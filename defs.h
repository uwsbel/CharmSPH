
#ifndef __DEFS__
#define __DEFS__

#include "pup.h"
#include "math.h"
#include <cstdio>

#define HYDROGEN_MASS           (1.67 * pow( 10.0,-24)) // in g
#define VDW_A                   (1.1328 * pow(10.0, -133)) // in (g m^2/s^2) m^12
#define VDW_B                   (2.23224 * pow(10.0, -76)) // in (g m^2/s^2) m^6

#define ENERGY_VAR              (1.0 * pow(10.0,-5))

//average of next two should be what you want as your atom density
//this should comply with the PERDIM parameter; for KAWAY 1 1 1, the maximum number
//of particles can be 10*10*10 = 1000 - 10 comes from PERDIM parameter, which is
//currently set to be 10, using a GAP of 3; as KWAYness increases, the maximum
//number of particles decreases - for 2 1 1, it is 500, for 2 2 1 it is 250; you
//can set them to have lower values but not higher; alternatively a host of
//paramters including PTP_CUT_OFF, PERDIM, GAP can be set to suitable values to
#define PARTICLES_PER_CELL_START        100 
#define PARTICLES_PER_CELL_END          250

#define DEFAULT_DELTA           1	// in femtoseconds

#define DEFAULT_FIRST_LDB       20
#define DEFAULT_LDB_PERIOD      20
#define DEFAULT_FT_PERIOD       100000

#define KAWAY_X                 1 //2 Original val
#define KAWAY_Y                 1 //2
#define KAWAY_Z                 1 //1
#define NBRS_X	                ((2 * KAWAY_X) + 1)
#define NBRS_Y                  ((2 * KAWAY_Y) + 1)
#define NBRS_Z                  ((2 * KAWAY_Z) + 1)
#define NUM_NEIGHBORS           (NBRS_X * NBRS_Y * NBRS_Z)
 
#define PI 3.1415926535897932384626433832795028841971693993751058
#define INVPI 0.3183098861837906715377675267450287240689192914809128


#define DT                      (5e-6)
#define H                       (0.05)
#define MarkDistMult            (1.0)
#define RHO0                    (1000)
#define PARTICLE_MASS           (H * H * H * RHO0) // 8 Particles per spatial cube/chare
#define MU                      (0.001)
#define GRAVITY                 (-9.81)
#define PRESSURE_CONSTANT       (0.5)
#define EPSILON                 (1e-2)
#define BASEPRES                (0)
#define MAXVEL                  (10)
#define BOUNDARY_PRESSURE       (1000) // Artificial Boundary Pressure
#define MULTVISCOSITY_FSI       (5.0)

#define MIN_X 0
#define MIN_Y 0
#define MIN_Z 0
#define MAX_X 1
#define MAX_Y 1
#define MAX_Z 1
#define FLUIDMIN_X 3 * H
#define FLUIDMIN_Y 3 * H
#define FLUIDMIN_Z 3 * H
#define FLUIDMAX_X 15 * H
#define FLUIDMAX_Y 15 * H
#define FLUIDMAX_Z 15 * H



#define CELLARRAY_DIM_X         3
#define CELLARRAY_DIM_Y         3
#define CELLARRAY_DIM_Z         3
#define PTP_CUT_OFF             H // cut off for atom to atom interactions
#define CELL_MARGIN             0  // constant diff between cutoff and cell size
#define CELL_SIZEX             (2 * PTP_CUT_OFF)/KAWAY_X // 
#define CELL_SIZEY             (2 * PTP_CUT_OFF)/KAWAY_Y
#define CELL_SIZEZ             (2 * PTP_CUT_OFF)/KAWAY_Z


//variables to control initial uniform placement of atoms;
//atoms should not be too close at startup for a stable system;  
//PERDIM * GAP should be less than (PTPCUTOFF+CELL_MARGIN);
//max particles per cell should not be greater thatn PERDIM^3 for 1 AWAY;
#define PERDIM                  10
#define GAP                     3 

#define CELL_ORIGIN_X           0
#define CELL_ORIGIN_Y	          0
#define CELL_ORIGIN_Z	          0

#define MIGRATE_STEPCOUNT	      5
#define DEFAULT_FINALSTEPCOUNT	1001
#define MAX_VELOCITY		        .1  //in A/fs

#define WRAP_X(a)		(((a) + cellArrayDim.x) % cellArrayDim.x)
#define WRAP_Y(a)		(((a) + cellArrayDim.y) % cellArrayDim.y)
#define WRAP_Z(a)		(((a) + cellArrayDim.z) % cellArrayDim.z)

struct vec3 {
  double x, y, z;

  vec3(double d = 0.0) : x(d), y(d), z(d) { }
  vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) { }
  void print()
  {
    printf("x = %f, y = %f, z = %f \n", x, y, z);
  }

  inline vec3& operator += (const vec3 &rhs) {
    x += rhs.x; y += rhs.y; z += rhs.z;
    return *this;
  }
  inline vec3 operator+ (const vec3& rhs) const {
    return vec3(x + rhs.x, y + rhs.y, z + rhs.z);
  }
  inline vec3 operator+ (const double rhs) const {
    return vec3(x + rhs, y + rhs, z + rhs);
  }
  inline vec3& operator -= (const vec3 &rhs) {
    return *this += (rhs * -1.0);
  }
  inline vec3 operator* (const double d) const {
    return vec3(d*x, d*y, d*z);
  }
  inline vec3 operator* (const vec3 a) const {
    return vec3(x * a.x, y * a.y, z * a.z);
  }
  inline vec3 operator/ (const double a) const {
    return vec3(x / a, y / a, z / a);
  }
  inline vec3 operator/ (const vec3 a) const {
    return vec3(x / a.x, y / a.y, z / a.z);
  }
  inline vec3 operator- (const vec3& rhs) const {
    return vec3(x - rhs.x, y - rhs.y, z - rhs.z);
  }
  inline vec3 operator- (const double rhs) const {
    return vec3(x - rhs, y - rhs, z - rhs);
  }
};
inline double dot(const vec3& a, const vec3& b) {
  return a.x*b.x + a.y*b.y + a.z*b.z;
}
inline double magnitude(const vec3& a){
  return sqrt(dot(a, a));
}

struct vec4 {
  vec3 r;
  double l;

  vec4(double d = 0.0) : r(d), l(d) { }
  vec4(double x_, double y_, double z_, double l_) : r(x_, y_, z_), l(l_) { }
  vec4(vec3 r_, double l_)
  { 
    r = r_;
    l = l_;
  }

  inline vec4& operator += (const vec4 &rhs) {
    r += rhs.r;
    l += rhs.l;
    return *this;
  }
  inline vec4 operator+ (const vec4& rhs) const {
    return vec4(r + rhs.r, l + rhs.l);
  }
  inline vec4& operator -= (const vec4 &rhs) {
    return *this += (rhs * -1.0);
  }
  inline vec4 operator* (const double d) const {
    return vec4(r * d, l * d);
  }
  inline vec4 operator* (const vec4 a) const {
    return vec4(r * a.r, l * a.l);
  }
  inline vec4 operator- (const vec4& rhs) const {
    return vec4(r - rhs.r, l - rhs.l);
  }
};



PUPbytes(vec3)
PUPbytes(vec4)

//class for keeping track of the properties for a particle
struct Particle {
  double mass;
  int typeOfParticle; // -1 if fluid, 0 if boundary, 1 if body

  double pressure;
  double rho, dRho;
  //   Position, acceleration, velocity
  vec3 pos, vel, acc;
};

PUPbytes(Particle);

#include "leanmd.decl.h"

/* Charm++ Globals */
extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CProxy_Cell cellArray;
extern /* readonly */ CProxy_Compute computeArray;
extern /* readonly */ CkGroupID mCastGrpID;

/* SPH Globals */
extern /* readonly */ vec3 cellArrayDim;
extern /* readonly */ vec3 domainMin;
extern /* readonly */ vec3 domainMax;
extern /* readonly */ vec3 domainDim;
extern /* readonly */ vec3 fluidMin;
extern /* readonly */ vec3 fluidMax;
extern /* readonly */ vec3 cellSize;
extern /* readonly */ vec3 mDist;

/* Charm++ Runtime System Globals */
extern /* readonly */ int finalStepCount;
extern /* readonly */ int checkptStrategy;
extern /* readonly */ std::string logs;

#endif
