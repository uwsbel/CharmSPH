/* C++ includes */
#include <sstream>
#include <fstream>
#include <iostream>
/* Charm++ includes */
#include "ckmulticast.h"
#include "ckio.h"
/* CharmSPH includes */
#include "defs.h"
#include "charmsph.decl.h"
#include "Cell.h"

Cell::Cell() : inbrs(NUM_NEIGHBORS), stepCount(1), updateCount(0), computesList(NUM_NEIGHBORS) {
  //load balancing to be called when AtSync is called
  usesAtSync = true;

  double halfH = h / 2;
  double hSquared = h * h;
  double boundaryThickness = 3 * h;
  vec3 boundaryMin = domainMin + boundaryThickness;
  vec3 boundaryMax = domainMax - boundaryThickness;

  int myid = thisIndex.z + cellArrayDim.z * (thisIndex.y + cellArrayDim.y * thisIndex.x); 
  //std::cout << "Chare ID: " << myid << std::endl;
  myNumParts = 1;

  vec3 cellMin(thisIndex.x * cellSize.x, thisIndex.y * cellSize.y, thisIndex.z * cellSize.z);
  for (double px = 0.5 * mDist.x; px < cellSize.x; px += mDist.x) 
  {
    for (double py = 0.5 * mDist.y; py < cellSize.y; py += mDist.y) 
    {
      for (double pz = 0.5 * mDist.z; pz < cellSize.z; pz += mDist.z) 
      {
        Particle p = Particle();
        p.pos = vec3(px, py, pz) + cellMin;
        p.vel = vec3(0,0,0);
        p.acc = vec3(0,0,0);
        p.mass = particleMass;
        p.rho = RHO0;
        p.pressure = Eos(p.rho);
        /* Set as lower or top boundary particle (above and below z plane)*/
        if((p.pos.z < boundaryMin.z || p.pos.z > boundaryMax.z) || 
           (p.pos.x < boundaryMin.x || p.pos.x > boundaryMax.x) ||
           (p.pos.y < boundaryMin.y)) //|| p.pos.y > boundaryMax.y))
           //(p.pos.y < boundaryMin.y))
        {
          p.typeOfParticle = 0; // Boundary Marker
          numBoundaryMarkers++;
          //p.pressure = BOUNDARY_PRESSURE;
          particles.push_back(p);
        }
        else if((p.pos.z > (boundaryMin.z + 0.5 * h) && p.pos.z < (boundaryMax.z - 0.5 * h)) && 
                (p.pos.x > (boundaryMin.x + 0.5 * h) && p.pos.x < (boundaryMax.x / 2)) &&
                (p.pos.y > (boundaryMin.y + 0.5 * h) && p.pos.y < (boundaryMax.y / 2)))
        // else if(((p.pos.z > fluidMin.z) && (p.pos.z < (fluidMin.z + 10 * h))) && 
        //    ((p.pos.x > fluidMin.x) && (p.pos.x < (fluidMin.x + 10 * h))) &&
        //    ((p.pos.y > fluidMin.y) && (p.pos.y < (fluidMin.y + 10 * h))))
        {
          p.typeOfParticle = -1; // Fluid Marker
          numFluidMarkers++;
          particles.push_back(p);  
        }
      }
    }
  }
  // Copy particles into a second vector for Runge-Kutta Integration purposes.
  particles2 = particles;

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
  // CkPrintf("Creating computes...\n");
  // std::cout << "NumPe's:" << CkNumPes() << std::endl;
  // std::cout << "currPe:" << currPe << std::endl;
  // std::cout << "NBRS: " << NBRS_Y << std::endl;
  // std::cout << "inbrs: " << inbrs << std::endl;

  /* Single chare single compute tweak. Uncomment this and comment for loop for */
  // computesList.resize(1); // This is necessary because usually computesList is NUM_NEIGHBORS
  // CkArrayIndex6D index(1,1,1,1,1,1);
  // computeArray[index].insert((++currPe) % CkNumPes());
  // computesList[0] = index;
  
  /* Comment all of the following for loop to make the code serial. */
  for (int num = 0; num < inbrs; num++) 
  {
    /* The following computation gives us the following set of indeces (-1,-1,-1),(-1,-1,0)
     * ...,(0,0,0)...(1,1,1) for (dx,dy,dz).
     */
    dx = num / (NBRS_Y * NBRS_Z)                - NBRS_X/2;
    dy = (num % (NBRS_Y * NBRS_Z)) / NBRS_Z     - NBRS_Y/2;
    dz = num % NBRS_Z                           - NBRS_Z/2;

    if (num >= inbrs / 2)
    {
      px1 = x + KAWAY_X;
      py1 = y + KAWAY_Y;
      pz1 = z + KAWAY_Z;
      px2 = px1 + dx;
      py2 = py1 + dy;
      pz2 = pz1 + dz;
      if(px1 == 1 && py1 == 1 && pz1 == 1 && px2 == 1 && py2 == 1 && pz2 == 1)
      {
        std::cout << "1-All zeros" << std::endl;
      }
      CkArrayIndex6D index(px1, py1, pz1, px2, py2, pz2);
      computeArray[index].insert((++currPe) % CkNumPes());
      computesList[num] = index;
      numComputes++;
    } 
    else 
    {
      // these computes will be created by pairing cells
      px1 = WRAP_X(x + dx) + KAWAY_X;
      py1 = WRAP_Y(y + dy) + KAWAY_Y;
      pz1 = WRAP_Z(z + dz) + KAWAY_Z;
      px2 = px1 - dx;
      py2 = py1 - dy;
      pz2 = pz1 - dz;
      if(px1 == 1 && py1 == 1 && pz1 == 1 && px2 == 1 && py2 == 1 && pz2 == 1)
      {
        std::cout << "1All zeros" << std::endl;
      }
      CkArrayIndex6D index(px1, py1, pz1, px2, py2, pz2);
      computesList[num] = index;

    }
  } // end of for loop

  contribute(CkCallback(CkReductionTarget(Main,run),mainProxy));
}

//call multicast section creation
void Cell::createSection() 
{
  //knit the computes into a section
  mCastSecProxy = CProxySection_Compute::ckNew(computeArray.ckGetArrayID(), &computesList[0], computesList.size());

  //delegate the communication responsibility for this section to multicast library
  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  mCastSecProxy.ckSectionDelegate(mCastGrp);
  // Reduction target to compute proxies are reduceForcesSPH
  mCastGrp->setReductionClient(mCastSecProxy, new CkCallback(CkReductionTarget(Cell,reduceForcesSPH), thisProxy(thisIndex.x, thisIndex.y, thisIndex.z)));

}

// Function to start interaction among particles in neighboring cells as well as its own particles
void Cell::sendPositions(int iteration) 
{
  if(iteration == 1)
  {
    particles2 = particles;
  }

  unsigned int len = particles.size();
  //create the particle and control message to be sent to computes
  ParticleDataMsg* msg = new (len) ParticleDataMsg(thisIndex.x, thisIndex.y, thisIndex.z, len);

  // Create a message with all the particles in the cell. Calculate the pressure of fluid particles
  // before sending the message.
  
  int counter = 0;
  for(int i = 0; i < len; ++i)
  {
    msg->part[i] = particles2[i];
  }
  //CkPrintf("Finish sending positions at iteration %d.... Triggering calculateForces\n", iteration);
  mCastSecProxy.calculateForces(msg);
}

//send the atoms that have moved beyond my cell to neighbors
void Cell::migrateParticles(int step)
{
  int id = thisIndex.x + thisIndex.y*cellArrayDim.x + thisIndex.z*cellArrayDim.x*cellArrayDim.y;

  int x1, y1, z1;
  std::vector<std::vector<Particle> > outgoing;
  outgoing.resize(inbrs); // Resive to number of neighbor cells (27).
  ////CkPrintf("Check 1 from chare %d\n",id);

  int size = particles.size();
  for(std::vector<Particle>::reverse_iterator iter = particles.rbegin(); iter != particles.rend(); iter++) 
  {
    // x1, y1 and z1 have the neighbor indeces relative to current chare
    migrateToCell(*iter, x1, y1, z1); 
    if(x1!=0 || y1!=0 || z1!=0) 
    {
      outgoing[(x1+KAWAY_X)*NBRS_Y*NBRS_Z + (y1+KAWAY_Y)*NBRS_Z + (z1+KAWAY_Z)].push_back(wrapAround(*iter));
      //outgoing[(x1+KAWAY_X)*NBRS_Y*NBRS_Z + (y1+KAWAY_Y)*NBRS_Z + (z1+KAWAY_Z)].push_back((*iter));

      std::swap(*iter, particles[size - 1]);
      size--;
    }
  }
  ////CkPrintf("Check 2 from chare %d\n",id);

  particles.resize(size);
  for(int num = 0; num < inbrs; num++) 
  {
    x1 = num / (NBRS_Y * NBRS_Z)            - NBRS_X/2;
    y1 = (num % (NBRS_Y * NBRS_Z)) / NBRS_Z - NBRS_Y/2;
    z1 = num % NBRS_Z                       - NBRS_Z/2;
    cellArray(WRAP_X(thisIndex.x+x1), WRAP_Y(thisIndex.y+y1), WRAP_Z(thisIndex.z+z1)).receiveParticles(outgoing[num]);
  }
  ////CkPrintf("Check 3 from chare %d\n",id);

}

//check if the particle is to be moved
void Cell::migrateToCell(Particle p, int &px, int &py, int &pz) 
{
  // x, y, z give the coordinates of the bottom left part of the cell
  double x = thisIndex.x * cellSize.x + domainMin.x;
  double y = thisIndex.y * cellSize.y + domainMin.y;
  double z = thisIndex.z * cellSize.z + domainMin.z;
  px = py = pz = 0;

  if (p.pos.x < x) px = -1;
  else if (p.pos.x > (x+cellSize.x)) px = 1;

  if (p.pos.y < y) py = -1;
  else if (p.pos.y > (y+cellSize.y)) py = 1;

  if (p.pos.z < z) pz = -1;
  else if (p.pos.z > (z+cellSize.z)) pz = 1;
}


inline double velocityCheck(double inVelocity) 
{
  if(fabs(inVelocity) > maxVel) 
  {
    if(inVelocity < 0.0 ) return -maxVel;
    else return maxVel;
  } 
  else 
  {
    return inVelocity;
  }
}


// Function to update properties (i.e. acceleration, velocity and position) in particles
// add double dt arg, use particles 1 or 2
void Cell::updatePropertiesSPH(vec4 *dVel_dRho, int iteration) 
{
  int i;
  int typeOfParticle;
  if(iteration == 1)
  {
    for(i = 0; i < particles2.size(); i++) 
    {
      typeOfParticle = particles2[i].typeOfParticle;
      if(typeOfParticle == -1)
      {
        particles2[i].acc = dVel_dRho[i].r + gravity;
        particles2[i].pos += particles2[i].vel * 0.5 * dt;
        particles2[i].vel += particles2[i].acc * 0.5 * dt; 
      }
      particles2[i].rho += dVel_dRho[i].l * 0.5 * dt;
      particles2[i].pressure = Eos(particles2[i].rho);
    }   
  }
  else
  {
    for(i = 0; i < particles.size(); i++) 
    {
      typeOfParticle = particles[i].typeOfParticle;
      if(typeOfParticle == -1)
      {
        particles[i].acc = dVel_dRho[i].r + gravity;
        particles[i].pos += particles[i].vel * dt;
        particles[i].vel += particles[i].acc * dt; 

        particles[i].vel.x = velocityCheck(particles[i].vel.x);
        particles[i].vel.y = velocityCheck(particles[i].vel.y);
        particles[i].vel.z = velocityCheck(particles[i].vel.z);

      }
       particles[i].rho += dVel_dRho[i].l * dt; // With constant presure the density shouldnt beupdated
       particles[i].pressure = Eos(particles[i].rho);

    } 
  }
}

Particle& Cell::wrapAround(Particle &p) 
{
  if(p.pos.x < domainMin.x) p.pos.x += domainDim.x;
  if(p.pos.y < domainMin.y) p.pos.y += domainDim.y;
  if(p.pos.z < domainMin.z) p.pos.z += domainDim.z;

  if(p.pos.x > domainMax.x) p.pos.x -= domainDim.x;
  if(p.pos.y > domainMax.y) p.pos.y -= domainDim.y;
  if(p.pos.z > domainMax.z) p.pos.z -= domainDim.z;

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

void Cell::writeCell(int stepCount)
{
    int id = thisIndex.x + thisIndex.y*cellArrayDim.x + thisIndex.z*cellArrayDim.x*cellArrayDim.y;
    std::stringstream ssFluidParticles;
    std::stringstream ssBoundaryParticles;

    ssFluidParticles << "x,y,z,";
    ssFluidParticles << "xVelocity,yVelocity,zVelocity,";
    ssFluidParticles << "xAcc,yAcc,zAcc,";
    ssFluidParticles << "velMagnitude,density,pressure";
    ssFluidParticles << std::endl;
    ssBoundaryParticles << "x,y,z,";
    ssBoundaryParticles << "xVelocity,yVelocity,zVelocity,";
    ssBoundaryParticles << "xAcc,yAcc,zAcc,";
    ssBoundaryParticles << "velMagnitude,density,pressure";
    ssBoundaryParticles << std::endl;

    for(int i = 0;i < particles.size();i++) {
      Particle p = particles[i];
      if(p.typeOfParticle==-1) {
        ssFluidParticles << p.pos.x << ',';
        ssFluidParticles << p.pos.y << ',';
        ssFluidParticles << p.pos.z << ',';
        ssFluidParticles << p.vel.x << ',';
        ssFluidParticles << p.vel.y << ',';
        ssFluidParticles << p.vel.z << ',';
        ssFluidParticles << p.acc.x << ',';
        ssFluidParticles << p.acc.y << ',';
        ssFluidParticles << p.acc.z << ',';
        ssFluidParticles << sqrt(dot(p.vel,p.vel)) << ',';
        ssFluidParticles << p.rho << ',';
        ssFluidParticles << p.pressure << ',';
        ssFluidParticles << std::endl;
      }
      else if((writeBoundary==1) && p.typeOfParticle==0) {
        ssBoundaryParticles << p.pos.x << ',';
        ssBoundaryParticles << p.pos.y << ',';
        ssBoundaryParticles << p.pos.z << ',';
        ssBoundaryParticles << p.vel.x << ',';
        ssBoundaryParticles << p.vel.y << ',';
        ssBoundaryParticles << p.vel.z << ',';
        ssBoundaryParticles << p.acc.x << ',';
        ssBoundaryParticles << p.acc.y << ',';
        ssBoundaryParticles << p.acc.z << ',';
        ssBoundaryParticles << sqrt(dot(p.vel,p.vel)) << ',';
        ssBoundaryParticles << p.rho << ',';
        ssBoundaryParticles << p.pressure << ',';
        ssBoundaryParticles << std::endl;
      }
    }

    std::ofstream fileNameFluid, fileNameBoundary;
    std::stringstream ssFileNameFluid, ssFileNameBoundary;
    ssFileNameFluid << "output/" << simID << "/fluid/fluid." << stepCount << ".chare." << id;
    ssFileNameBoundary << "output/" << simID << "/boundary/boundary." << stepCount << ".chare." << id;

    fileNameFluid.open(ssFileNameFluid.str().c_str());
    fileNameFluid << ssFluidParticles.str();
    fileNameFluid.close();

    if(writeBoundary == 1){
      fileNameBoundary.open(ssFileNameBoundary.str().c_str());
      fileNameBoundary << ssBoundaryParticles.str();
      fileNameBoundary.close();
    }

}



