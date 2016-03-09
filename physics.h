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
  if (q < 1) 
  {
    return (0.25f / (PI * H * H * H) * (pow(2 - q, 3) - 4 * pow(1 - q, 3)));
  }
  if (q < 2) 
  {
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
    double constant1 = 0.75 * (INVPI) * pow(H, -5) * ((3 * q) - 4);
    return d * constant1;
  }
  if (q < 2) 
  {
    double constant1 = 0.75 * (INVPI) * pow(H, -5) * (-q + 4.0 - (4.0 / q));
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

        dVel_i = calcDVel(rho_i, vel_i, p_i, rho_j, vel_j, p_j, gradW, r_ij_dot_gradW_overDist,multViscosity);
        dRho_i = calcDRho(rho_i, vel_i, rho_j, vel_j, gradW);

        dVel_dRho[i].r += dVel_i;
        dVel_dRho[i].l += dRho_i;

      }
    } 
  }
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

          dVel_ij = calcDVel(rho_i, vel_i, p_i, rho_j, vel_j, p_j, gradW, r_ij_dot_gradW_overDist,multViscosity);

          dRho_ij = calcDRho(rho_i, vel_i, rho_j, vel_j, gradW);

          dVel_dRho1[i].r += dVel_ij; 
          dVel_dRho1[i].l += dRho_ij; 

          dVel_dRho2[j].r -= dVel_ij; 
          dVel_dRho2[j].l += dRho_ij; 
        }
      }
    }
  }
}


#endif
