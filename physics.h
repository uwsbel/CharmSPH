#ifndef __PHYSICS_H__
#define __PHYSICS_H__

#include "ckmulticast.h"
#include "defs.h"

#define BLOCK_SIZE	512


/**
 * @brief Distance
 * @details 
 *          Distance between two particles, considering the periodic boundary condition
 * 
 * @param posRadA Position of Particle A
 * @param posRadB Position of Particle B
 * 
 * @return Distance vector (distance in x, distance in y, distance in z)
 */
inline vec3 Distance(vec3 a, vec3 b) {
  vec3 dist3 = a - b;
  // dist3.x -= ((dist3.x > 0.5f * domainDim.x) ? domainDim.x : 0);
  // dist3.x += ((dist3.x < -0.5f * domainDim.x) ? domainDim.x : 0);

  // dist3.y -= ((dist3.y > 0.5f * domainDim.y) ? domainDim.y : 0);
  // dist3.y += ((dist3.y < -0.5f * domainDim.y) ? domainDim.y : 0);

  // dist3.z -= ((dist3.z > 0.5f * domainDim.z) ? domainDim.z : 0);
  // dist3.z += ((dist3.z < -0.5f * domainDim.z) ? domainDim.z : 0);
  return dist3;
}

//3D SPH kernel function, W3_SplineA
inline double W3_Spline(double d) 
{ // d is positive. h is the sph particle radius (i.e. h in the document) d is the distance of 2 particles
  double q = fabs(d) / H;
  if (q < 1) {
    return (0.25f / (PI * H * H * H) * (pow(2 - q, 3) - 4 * pow(1 - q, 3)));
  }
  if (q < 2) {
    return (0.25f / (PI * H * H * H) * pow(2 - q, 3));
  }
  return 0;
}
//Gradient of the kernel function
// d: magnitude of the distance of the two particles
// dW * dist3 gives the gradiant of W3_Quadratic, where dist3 is the distance vector of the two particles, (dist3)a = pos_a - pos_b
inline vec3 GradW_Spline(vec3 d) { // d is positive. r is the sph particle radius (i.e. h in the document) d is the distance of 2 particles
  double q = sqrt(dot(d, d)) / H;
  if (q < 1) 
  {
    double constant1 = 0.75 * (INVPI) * pow(H, -5) * (3 * q - 4);
    return d * constant1;
  }
  if (q < 2) 
  {
    double constant1 = 0.75 * (INVPI) * pow(H, -5) * (-q + 4.0 - 4.0 / q);
    return d * constant1;
  }
  return vec3(0, 0, 0);
}

 //! Lennard-Jones boundary repulsion force
inline float calcLJForce(const float r)
{
  float force = 0.0f;
  float dcoeff =  abs(GRAVITY);// * (domainMax.y);
  int p1coeff = 12;
  int p2coeff = 6;
  if (r <= H)
  {
    force = dcoeff * (pow(H/r, p1coeff) - pow(H/r, p2coeff));   
  }

   return force;
}

inline double calcDRho(double rho_i, vec3 vel_i, 
                       double rho_j, vec3 vel_j, 
                       vec3 gradW)
{
  return (rho_i / rho_j) * PARTICLE_MASS * dot(vel_i - vel_j, gradW);
}

inline vec3 calcDVel(double rho_i, vec3 vel_i, double p_i, 
                     double rho_j, vec3 vel_j, double p_j, 
                     vec3 gradW, double r_ij_dot_gradW_overDist, double multViscosity)
{
  return gradW * -1 * PARTICLE_MASS * (p_i / (rho_i * rho_i) + p_j / (rho_j * rho_j)) +  (vel_i - vel_j) * PARTICLE_MASS * (8 * multViscosity) * MU * pow(rho_i + rho_j, -2) * r_ij_dot_gradW_overDist;
}

inline void calcPairForcesSPH(ParticleDataMsg* first, ParticleDataMsg* second, int stepCount, std::vector<vec4>& dVel_dRho1, std::vector<vec4>& dVel_dRho2)
{
  int i, j;
  int firstLen = first->lengthAll;
  int secondLen = second->lengthAll;

  dVel_dRho1.resize(firstLen);
  dVel_dRho2.resize(secondLen);

  vec3 pos_i, pos_j, vel_i, vel_j, r_ij;
  double p_i, p_j, rho_i, rho_j, absDist;
  int typeOfParticle_i, typeOfParticle_j;
  vec3 gradW;
  vec3 dVel_ij;
  double dRho_ij;

  int i1, j1;
  for(i1 = 0; i1 < firstLen; i1=i1+BLOCK_SIZE)
  {
    for(j1 = 0; j1 < secondLen; j1=j1+BLOCK_SIZE)
    {
      for(i = i1; i < i1+BLOCK_SIZE && i < firstLen; i++) 
      {
        pos_i            = first->part[i].pos;
        vel_i            = first->part[i].vel;
        p_i              = first->part[i].pressure;
        rho_i            = first->part[i].rho;
        typeOfParticle_i = first->part[i].typeOfParticle;
        for(j = j1; j < j1+BLOCK_SIZE && j < secondLen; j++) 
        {
  
          pos_j            = second->part[j].pos;
          vel_j            = second->part[j].vel;
          p_j              = second->part[j].pressure;
          rho_j            = second->part[j].rho;
          typeOfParticle_j = second->part[j].typeOfParticle;

          r_ij = Distance(pos_i, pos_j);
          absDist = magnitude(r_ij);
          
          if (absDist > PTP_CUT_OFF)
          {
            continue;
          }
          if (typeOfParticle_i >= 0  &&  typeOfParticle_j >= 0) {
            continue;
          }
          double multViscosity = 1.0;
          gradW = GradW_Spline(r_ij);
          double r_ij_dot_gradW = dot(r_ij, gradW);
          double r_ij_dot_gradW_overDist = r_ij_dot_gradW / (absDist * absDist + EPSILON * H * H);

          if ((typeOfParticle_i == -1)  &&  (typeOfParticle_j == -1)) 
          {
            dVel_ij = calcDVel(rho_i, vel_i, p_i, rho_j, vel_j, p_j, gradW, r_ij_dot_gradW_overDist,multViscosity);
          }
          //else
          //{
          // dVel_ij = (r_ij / absDist) * calcLJForce(absDist);
          //}

          dRho_ij = calcDRho(rho_i, vel_i, rho_j, vel_j, gradW);

          dVel_dRho1[i].r += dVel_ij; 
          dVel_dRho1[i].l += dRho_ij; 

          dVel_dRho2[j].r -= dVel_ij; 
          dVel_dRho2[j].l -= dRho_ij; 
        }
      }
    }
  }
}

inline void calcInternalForcesSPH(ParticleDataMsg* first, int stepCount, std::vector<vec4>& dVel_dRho)
{
  int i, j;
  int firstLen = first->lengthAll;

  //CkPrintf("Calculating Internal Forces...\n");
  //CkPrintf("NumParticles: %d\n", firstLen);


  dVel_dRho.resize(firstLen);

  vec3 pos_i, pos_j, vel_i, vel_j, r_ij;
  double p_i, p_j, rho_i, rho_j, absDist;
  int typeOfParticle_i, typeOfParticle_j;
  vec3 gradW;
  vec3 dVel_i = vec3(0,0,0);
  double dRho_i = 0.0;
  //CkPrintf("firstLen = %d\n",firstLen);

  for(i = 0; i < firstLen; i++)
  {
    pos_i            = first->part[i].pos;
    vel_i            = first->part[i].vel;
    p_i              = first->part[i].pressure;
    rho_i            = first->part[i].rho;
    typeOfParticle_i = first->part[i].typeOfParticle;
    for(j = 0; j < firstLen; j++) 
    {     
      if(i != j)
      {
        pos_j            = first->part[j].pos;
        vel_j            = first->part[j].vel;
        p_j              = first->part[j].pressure;
        rho_j            = first->part[j].rho;
        typeOfParticle_j = first->part[j].typeOfParticle;

        r_ij = Distance(pos_i, pos_j);
        absDist = magnitude(r_ij);


        if (absDist > PTP_CUT_OFF)
        {
          continue;
        }
        if (typeOfParticle_i >= 0  &&  typeOfParticle_j >= 0) 
        {
          continue;
        }
        double multViscosity = 1.0;
        gradW = GradW_Spline(r_ij);
        double r_ij_dot_gradW = dot(r_ij, gradW);
        double r_ij_dot_gradW_overDist = r_ij_dot_gradW / (absDist * absDist + EPSILON * H * H);

        if ((typeOfParticle_i == -1)  &&  (typeOfParticle_j == -1)) 
        {
          dVel_i = calcDVel(rho_i, vel_i, p_i, rho_j, vel_j, p_j, gradW, r_ij_dot_gradW_overDist,multViscosity);
        }
        //else
        //{
        //  dVel_i = (r_ij / absDist) * calcLJForce(absDist);
        //}

        dRho_i = calcDRho(rho_i, vel_i, rho_j, vel_j, gradW);

        dVel_dRho[i].r += dVel_i;
        dVel_dRho[i].l += dRho_i;
      }
    } 
  }
}
//function to calculate forces among atoms in a single list
inline double calcInternalForces(ParticleDataMsg* first, int stepCount, std::vector<vec3>& force1) {
  int i, j, ptpCutOffSqd;
  int firstLen = first->lengthAll;
  double powTwenty, powTen, firstx, firsty, firstz, rx, ry, rz, r, rsqd, fx, fy, fz, f, fr;
  vec3 firstpos, separation, force;
  double rSix, rTwelve;
  double energy = 0;
  int doEnergy = 0;
  if(stepCount == 1 || stepCount == finalStepCount)
    doEnergy = 1;
  force1.resize(firstLen);

  ptpCutOffSqd = PTP_CUT_OFF * PTP_CUT_OFF;
  powTen = pow(10.0, -10);
  powTwenty = pow(10.0, -20);
  for(i = 0; i < firstLen; i++){
    firstpos = first->part[i].pos;
    for(j = i+1; j < firstLen; j++) {
      // computing base values
      separation = firstpos - first->part[j].pos;
      rsqd = dot(separation, separation);
      if(rsqd > 1 && rsqd < ptpCutOffSqd){
        rsqd = rsqd * powTwenty;
        r = sqrt(rsqd);
        rSix = ((double)rsqd) * rsqd * rsqd;
        rTwelve = rSix * rSix;
        f = (double)( (12 * VDW_A) / rTwelve - (6 * VDW_B) / rSix);
        if(doEnergy)
          energy += (double)( VDW_A / rTwelve - VDW_B / rSix);

        fr = f / rsqd;
        force = separation * (fr * powTen);
        force1[i] += force;
        force1[j] -= force;
      }
    }
  }
  return energy;
}

//function to calculate forces among 2 lists of atoms
inline double calcPairForces(ParticleDataMsg* first, ParticleDataMsg* second, int stepCount,
                             std::vector<vec3>& force1, std::vector<vec3> &force2) {
  int i, j, ptpCutOffSqd, diff;
  int firstLen = first->lengthAll;
  int secondLen = second->lengthAll;
  double powTwenty, powTen, r, rsqd, f, fr;
  vec3 separation, force;
  double rSix, rTwelve;
  double energy = 0;
  int doEnergy = 0;
  if(stepCount == 1 || stepCount == finalStepCount)
    doEnergy = 1;

  force1.resize(firstLen);
  force2.resize(secondLen);
  //check for wrap around and adjust locations accordingly
  if (abs(first->x - second->x) > 1){
    diff = cellSize.x * cellArrayDim.x;
    if (second->x < first->x){
      diff = -1 * diff; 
    }
    for (i = 0; i < firstLen; i++)
      first->part[i].pos.x += diff;
  }
  if (abs(first->y - second->y) > 1){
    diff = cellSize.y * cellArrayDim.y;
    if (second->y < first->y)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->part[i].pos.y += diff;
  }
  if (abs(first->z - second->z) > 1){
    diff = cellSize.z * cellArrayDim.z;
    if (second->z < first->z)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->part[i].pos.z += diff;
  } 
  ptpCutOffSqd = PTP_CUT_OFF * PTP_CUT_OFF;
  powTen = pow(10.0, -10);
  powTwenty = pow(10.0, -20);

  int i1, j1;
  for(i1 = 0; i1 < firstLen; i1=i1+BLOCK_SIZE)
    for(j1 = 0; j1 < secondLen; j1=j1+BLOCK_SIZE)
      for(i = i1; i < i1+BLOCK_SIZE && i < firstLen; i++) {
        for(j = j1; j < j1+BLOCK_SIZE && j < secondLen; j++) {
          separation = first->part[i].pos - second->part[j].pos;
          rsqd = dot(separation, separation);
          if (rsqd > 1 && rsqd < ptpCutOffSqd) {
            rsqd = rsqd * powTwenty;
            r = sqrt(rsqd);
            rSix = ((double)rsqd) * rsqd * rsqd;
            rTwelve = rSix * rSix;
            f = (double)( (12 * VDW_A) / rTwelve - (6 * VDW_B) / rSix);
            if(doEnergy)
              energy += (double)( VDW_A / rTwelve - VDW_B / rSix); // in milliJoules
            fr = f / rsqd;
            force = separation * (fr * powTen);
            force1[i] += force;
            force2[j] -= force;
          }
        }
      }

  return energy;
}

#endif
