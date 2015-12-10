#include <sstream>
#include <fstream>
#include <iostream>

#include "defs.h"
#include "leanmd.decl.h"
#include "Cell.h"
#include "ckmulticast.h"

Cell::Cell() : inbrs(NUM_NEIGHBORS), stepCount(1), updateCount(0), computesList(NUM_NEIGHBORS) {
  //load balancing to be called when AtSync is called
  usesAtSync = true;

  int numParticlesToAdd = 8;
  double halfH = H / 2;
  double HSquared = H * H;

  int myid = thisIndex.z+cellArrayDimZ*(thisIndex.y+thisIndex.x*cellArrayDimY); 
  myNumParts = 1;
  vec3 center((thisIndex.x * H), (thisIndex.y * H), (thisIndex.z * H));
  vec3 particlesToAdd[8];
  particlesToAdd[0] = vec3(center.x + halfH, center.y + halfH, center.z + halfH);
  particlesToAdd[1] = vec3(center.x + halfH, center.y + halfH, center.z - halfH);
  particlesToAdd[2] = vec3(center.x + halfH, center.y - halfH, center.z + halfH);
  particlesToAdd[3] = vec3(center.x + halfH, center.y - halfH, center.z - halfH);
  particlesToAdd[4] = vec3(center.x - halfH, center.y + halfH, center.z + halfH);
  particlesToAdd[5] = vec3(center.x - halfH, center.y + halfH, center.z - halfH);
  particlesToAdd[6] = vec3(center.x - halfH, center.y - halfH, center.z + halfH);
  particlesToAdd[7] = vec3(center.x - halfH, center.y - halfH, center.z - halfH);

  for(int i = 0;i < numParticlesToAdd;i++)
  {
    Particle p = Particle();
    p.pos = particlesToAdd[i];
    p.vel = vec3(0,0,0);
    p.acc = vec3(0,0,0);
    p.mass = PARTICLE_MASS;
    p.rho = RHO0;
    p.pressure = BOUNDARY_PRESSURE;
    /* Set as lower or top boundary particle (above and below z plane)*/
    if((p.pos.z > fluidMin.z && p.pos.z < fluidMax.z) && 
       (p.pos.x > fluidMin.x && p.pos.x < fluidMax.x) &&
       (p.pos.y > fluidMin.y && p.pos.y < fluidMax.y))
    {
      p.typeOfParticle = -1; // fluid Marker
      p.pressure = Eos(p.rho);

      particles.push_back(p);
    }
    /* */
    // else if(p.pos.x > fluidMin.x || p.pos.x < fluidMax.x)
    // {
    //   p.typeOfParticle = 0;
    //   particles.push_back(p);
    // }
    // else if(p.pos.y <= fluidMin.y || p.pos.y >= fluidMax.y)
    // {
    //   p.typeOfParticle = 0;
    //   particles.push_back(p);  
    // }
    else
    {
      p.typeOfParticle = 0; // Boundary
      particles.push_back(p);  
    }
  }


  energy[0] = energy[1] = 0;
  setMigratable(false);
}

//constructor for chare object migration
Cell::Cell(CkMigrateMessage *msg): CBase_Cell(msg) {
  usesAtSync = true;
  setMigratable(false);
  delete msg;
}  

Cell::~Cell() {}

//function to create my computes
void Cell::createComputes() {
  int x = thisIndex.x, y = thisIndex.y, z = thisIndex.z;
  int px1, py1, pz1, dx, dy, dz, px2, py2, pz2;

  /*  The computes X are inserted by a given cell:
   *
   *	^  X  X  X
   *	|  0  X  X
   *	y  0  0  0
   *	   x ---->
   */

  // for round robin insertion
  int currPe = CkMyPe();
  //CkPrintf("Creating computes...\n");

  for (int num = 0; num < inbrs; num++) {

    dx = num / (NBRS_Y * NBRS_Z)                - NBRS_X/2;
    dy = (num % (NBRS_Y * NBRS_Z)) / NBRS_Z     - NBRS_Y/2;
    dz = num % NBRS_Z                           - NBRS_Z/2;


    if (num >= inbrs / 2){
      px1 = x + KAWAY_X;
      py1 = y + KAWAY_Y;
      pz1 = z + KAWAY_Z;
      px2 = px1+dx;
      py2 = py1+dy;
      pz2 = pz1+dz;

      CkArrayIndex6D index(px1, py1, pz1, px2, py2, pz2);
      computeArray[index].insert((++currPe) % CkNumPes());
      computesList[num] = index;
    } else {
      // these computes will be created by pairing cells
      px1 = WRAP_X(x + dx) + KAWAY_X;
      py1 = WRAP_Y(y + dy) + KAWAY_Y;
      pz1 = WRAP_Z(z + dz) + KAWAY_Z;
      px2 = px1 - dx;
      py2 = py1 - dy;
      pz2 = pz1 - dz;
      CkArrayIndex6D index(px1, py1, pz1, px2, py2, pz2);
      computesList[num] = index;
    }
  } // end of for loop
  contribute(CkCallback(CkReductionTarget(Main,run),mainProxy));
}

//call multicast section creation
void Cell::createSection() {
  //knit the computes into a section
  mCastSecProxy = CProxySection_Compute::ckNew(computeArray.ckGetArrayID(), &computesList[0], computesList.size());

  //delegate the communication responsibility for this section to multicast library
  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  mCastSecProxy.ckSectionDelegate(mCastGrp);
  //mCastGrp->setReductionClient(mCastSecProxy, new CkCallback(CkReductionTarget(Cell,reduceForces), thisProxy(thisIndex.x, thisIndex.y, thisIndex.z)));
  // Reduction target to compute proxies are reduceForcesSPH
  mCastGrp->setReductionClient(mCastSecProxy, new CkCallback(CkReductionTarget(Cell,reduceForcesSPH), thisProxy(thisIndex.x, thisIndex.y, thisIndex.z)));

}

// Function to start interaction among particles in neighboring cells as well as its own particles
void Cell::sendPositions() {
  unsigned int len = particles.size();
  //create the particle and control message to be sent to computes
  ParticleDataMsg* msg = new (len) ParticleDataMsg(thisIndex.x, thisIndex.y, thisIndex.z, len);
  int id = thisIndex.x + thisIndex.y*cellArrayDimX + thisIndex.z*cellArrayDimX*cellArrayDimY;

  // Create a messahe with all the particles in the cell. Calculate the pressure of fluid particles
  // before sending the message.
  for(int i = 0; i < len; ++i)
  {
    if(particles[i].typeOfParticle < 0)
    {
      particles[i].pressure = Eos(particles[i].rho);
    }
    msg->part[i] = particles[i];
  }

  mCastSecProxy.calculateForces(msg);
}

void Cell::writeCell(int stepCount)
{
    int id = thisIndex.x + thisIndex.y*cellArrayDimX + thisIndex.z*cellArrayDimX*cellArrayDimY;

    std::stringstream ssParticles;
    ssParticles << "x,";
    ssParticles << "y,";
    ssParticles << "z,";
    ssParticles << "xVelocity,";
    ssParticles << "yVelocity,";
    ssParticles << "zVelocity,";
    ssParticles << "velMagnitude,";
    ssParticles << "mass,";
    ssParticles << "typeOfParticle";
    ssParticles << std::endl;

    for(int i = 0;i < particles.size();i++)
    {
      Particle p = particles[i];
      ssParticles << p.pos.x << ',';
      ssParticles << p.pos.y << ',';
      ssParticles << p.pos.z << ',';
      ssParticles << p.vel.x << ',';
      ssParticles << p.vel.y << ',';
      ssParticles << p.vel.z << ',';
      ssParticles << sqrt(dot(p.vel,p.vel)) << ',';
      ssParticles << p.mass << ',';
      ssParticles << p.typeOfParticle;
      ssParticles << std::endl;
    }

    std::ofstream fileNameParticles;
    std::stringstream ssFileNameParticles;
    ssFileNameParticles << "output/step." << stepCount << ".chare." << id;

    fileNameParticles.open(ssFileNameParticles.str().c_str());
    fileNameParticles<<ssParticles.str();
    fileNameParticles.close();
}

//send the atoms that have moved beyond my cell to neighbors
void Cell::migrateParticles()
{
  int id = thisIndex.x + thisIndex.y*cellArrayDimX + thisIndex.z*cellArrayDimX*cellArrayDimY;

  int x1, y1, z1;
  std::vector<std::vector<Particle> > outgoing;
  outgoing.resize(inbrs);
  //CkPrintf("Check 1 from chare %d",id);

  int size = particles.size();
  for(std::vector<Particle>::reverse_iterator iter = particles.rbegin(); iter != particles.rend(); iter++) 
  {
    migrateToCell(*iter, x1, y1, z1);
    if(x1!=0 || y1!=0 || z1!=0) 
    {
      outgoing[(x1+KAWAY_X)*NBRS_Y*NBRS_Z + (y1+KAWAY_Y)*NBRS_Z + (z1+KAWAY_Z)].push_back(wrapAround(*iter));
      std::swap(*iter, particles[size - 1]);
      size--;
    }
  }
  //CkPrintf("Check 2 from chare %d",id);

  particles.resize(size);
  for(int num = 0; num < inbrs; num++) 
  {
    x1 = num / (NBRS_Y * NBRS_Z)            - NBRS_X/2;
    y1 = (num % (NBRS_Y * NBRS_Z)) / NBRS_Z - NBRS_Y/2;
    z1 = num % NBRS_Z                       - NBRS_Z/2;
    cellArray(WRAP_X(thisIndex.x+x1), WRAP_Y(thisIndex.y+y1), WRAP_Z(thisIndex.z+z1)).receiveParticles(outgoing[num]);
  }
  //CkPrintf("Check 3 from chare %d",id);

}

//check if the particle is to be moved
void Cell::migrateToCell(Particle p, int &px, int &py, int &pz) {
  double x = thisIndex.x * CELL_SIZE_X + CELL_ORIGIN_X;
  double y = thisIndex.y * CELL_SIZE_Y + CELL_ORIGIN_Y;
  double z = thisIndex.z * CELL_SIZE_Z + CELL_ORIGIN_Z;
  px = py = pz = 0;

  if (p.pos.x < (x-CELL_SIZE_X)) px = -2;
  else if (p.pos.x < x) px = -1;
  else if (p.pos.x > (x+2*CELL_SIZE_X)) px = 2;
  else if (p.pos.x > (x+CELL_SIZE_X)) px = 1;

  if (p.pos.y < (y-CELL_SIZE_Y)) py = -2;
  else if (p.pos.y < y) py = -1;
  else if (p.pos.y > (y+2*CELL_SIZE_Y)) py = 2;
  else if (p.pos.y > (y+CELL_SIZE_Y)) py = 1;

  if (p.pos.z < (z-CELL_SIZE_Z)) pz = -2;
  else if (p.pos.z < z) pz = -1;
  else if (p.pos.z > (z+2*CELL_SIZE_Z)) pz = 2;
  else if (p.pos.z > (z+CELL_SIZE_Z)) pz = 1;
}

// // Function to update properties (i.e. acceleration, velocity and position) in particles
// void Cell::updateProperties(vec3 *forces) {
//   int i;
//   double powTen, powTwenty, realTimeDeltaVel, invMassParticle;
//   powTen = pow(10.0, 10);
//   powTwenty = pow(10.0, -20);
//   realTimeDeltaVel = DEFAULT_DELTA * powTwenty;
//   for(i = 0; i < particles.size(); i++) {
//     //calculate energy only in begining and end
//     if(stepCount == 1) {
//       energy[0] += (0.5 * particles[i].mass * dot(particles[i].vel, particles[i].vel) * powTen); // in milliJoules
//     } else if(stepCount == finalStepCount) { 
//       energy[1] += (0.5 * particles[i].mass * dot(particles[i].vel, particles[i].vel) * powTen);
//     }
//     // applying kinetic equations
//     invMassParticle = 1 / particles[i].mass;
//     particles[i].acc = forces[i] * invMassParticle; // in m/sec^2
//     particles[i].vel += particles[i].acc * realTimeDeltaVel; // in A/fs

//     limitVelocity(particles[i]);
//     particles[i].pos += particles[i].vel * DEFAULT_DELTA; // in A

//   }
// }

// Function to update properties (i.e. acceleration, velocity and position) in particles
void Cell::updatePropertiesSPH(vec4 *dVel_dRho) 
{
  int i;

  for(i = 0; i < particles.size(); i++) 
  {
    particles[i].acc = dVel_dRho[i].r;
    particles[i].dRho = dVel_dRho[i].l;

    particles[i].vel += dVel_dRho[i].r * DT;; 
    particles[i].pos += particles[i].vel * DT;
    particles[i].rho += dVel_dRho[i].l * DT;
  }
}

inline double velocityCheck(double inVelocity) {
  if(fabs(inVelocity) > MAX_VELOCITY) {
    if(inVelocity < 0.0 )
      return -MAX_VELOCITY;
    else
      return MAX_VELOCITY;
  } else {
    return inVelocity;
  }
}

void Cell::limitVelocity(Particle &p) {
  p.vel.x = velocityCheck(p.vel.x);
  p.vel.y = velocityCheck(p.vel.y);
  p.vel.z = velocityCheck(p.vel.z);
}

Particle& Cell::wrapAround(Particle &p) {
  if(p.pos.x < CELL_ORIGIN_X) p.pos.x += CELL_SIZE_X*cellArrayDimX;
  if(p.pos.y < CELL_ORIGIN_Y) p.pos.y += CELL_SIZE_Y*cellArrayDimY;
  if(p.pos.z < CELL_ORIGIN_Z) p.pos.z += CELL_SIZE_Z*cellArrayDimZ;
  if(p.pos.x > CELL_ORIGIN_X + CELL_SIZE_X*cellArrayDimX) p.pos.x -= CELL_SIZE_X*cellArrayDimX;
  if(p.pos.y > CELL_ORIGIN_Y + CELL_SIZE_Y*cellArrayDimY) p.pos.y -= CELL_SIZE_Y*cellArrayDimY;
  if(p.pos.z > CELL_ORIGIN_Z + CELL_SIZE_Z*cellArrayDimZ) p.pos.z -= CELL_SIZE_Z*cellArrayDimZ;

  return p;
}

//pack important data when I move/checkpoint
void Cell::pup(PUP::er &p) {
  CBase_Cell::pup(p);
  __sdag_pup(p);
  p | particles;
  p | stepCount;
  p | myNumParts;
  p | updateCount;
  p | stepTime;
  p | inbrs;
  p | numReadyCheckpoint;
  PUParray(p, energy, 2);

  p | computesList;

  p | mCastSecProxy;
  //adjust the multicast tree to give best performance after moving
  if (p.isUnpacking()){
    if(CkInRestarting()){
      createSection();
    }
    else{
      CkMulticastMgr *mg = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
      mg->resetSection(mCastSecProxy);
      //mg->setReductionClient(mCastSecProxy, new CkCallback(CkReductionTarget(Cell,reduceForces), thisProxy(thisIndex.x, thisIndex.y, thisIndex.z)));
      mg->setReductionClient(mCastSecProxy, new CkCallback(CkReductionTarget(Cell,reduceForcesSPH), thisProxy(thisIndex.x, thisIndex.y, thisIndex.z)));
    }
  }
}

