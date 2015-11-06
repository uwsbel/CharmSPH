
#ifndef __DEFS__
#define __DEFS__

#include "pup.h"
#include "math.h"

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

#define KAWAY_X                 1
#define KAWAY_Y                 1
#define KAWAY_Z                 1
#define NBRS_X	                (2 * KAWAY_X+1)
#define NBRS_Y                  (2 * KAWAY_Y+1)
#define NBRS_Z                  (2 * KAWAY_Z+1)
#define NUM_NEIGHBORS           (NBRS_X * NBRS_Y * NBRS_Z)
 
#define PI 3.1415926535897932384626433832795028841971693993751058
#define INVPI 0.3183098861837906715377675267450287240689192914809128


#define DT                      (1e-5)
#define H                       (0.05)
#define RHO0                    (1000)
#define PARTICLE_MASS           (H * H * RHO0)
#define MU                      (0.001)
#define GRAVITY                 (-9.81)
#define PRESSURE_CONSTANT       (0.5)
#define EPSILON                 (1e-2)
#define BASEPRES                (0)
#define MAXVEL                  (10)
#define BOUNDARY_PRESSURE       (1000) // Artificial Boundary Pressure
#define MULTVISCOSITY_FSI       (5.0)

#define CELLARRAY_DIM_X         3
#define CELLARRAY_DIM_Y         3
#define CELLARRAY_DIM_Z         3
#define PTP_CUT_OFF             2 * H // cut off for atom to atom interactions
#define CELL_MARGIN             0.01 * H  // constant diff between cutoff and cell size
#define CELL_SIZE_X             (PTP_CUT_OFF + CELL_MARGIN)/KAWAY_X
#define CELL_SIZE_Y             (PTP_CUT_OFF + CELL_MARGIN)/KAWAY_Y
#define CELL_SIZE_Z             (PTP_CUT_OFF + CELL_MARGIN)/KAWAY_Z


//variables to control initial uniform placement of atoms;
//atoms should not be too close at startup for a stable system;  
//PERDIM * GAP should be less than (PTPCUTOFF+CELL_MARGIN);
//max particles per cell should not be greater thatn PERDIM^3 for 1 AWAY;
#define PERDIM                  10
#define GAP                     3 

#define CELL_ORIGIN_X           0
#define CELL_ORIGIN_Y	          0
#define CELL_ORIGIN_Z	          0

#define MIGRATE_STEPCOUNT	      20
#define DEFAULT_FINALSTEPCOUNT	1001
#define MAX_VELOCITY		        .1  //in A/fs

#define WRAP_X(a)		(((a) + cellArrayDimX) % cellArrayDimX)
#define WRAP_Y(a)		(((a) + cellArrayDimY) % cellArrayDimY)
#define WRAP_Z(a)		(((a) + cellArrayDimZ) % cellArrayDimZ)

struct vec3 {
  double x, y, z;

  vec3(double d = 0.0) : x(d), y(d), z(d) { }
  vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) { }

  inline vec3& operator += (const vec3 &rhs) {
    x += rhs.x; y += rhs.y; z += rhs.z;
    return *this;
  }
  inline vec3 operator+ (const vec3& rhs) const {
    return vec3(x + rhs.x, y + rhs.y, z + rhs.z);
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
  inline vec3 operator- (const vec3& rhs) const {
    return vec3(x - rhs.x, y - rhs.y, z - rhs.z);
  }
};
inline double dot(const vec3& a, const vec3& b) {
  return a.x*b.x + a.y*b.y + a.z*b.z;
}
inline double magnitude(const vec3& a){
  return sqrt(dot(a, a));
}


PUPbytes(vec3)

//class for keeping track of the properties for a particle
struct Particle {
  double mass;
  int typeOfParticle; // -1 if fluid, 0 if boundary, 1 if body

  double pressure;
  double rho, rhoHalf, dRho;
  //   Position, acceleration, velocity
  vec3 pos, vel, acc;
  vec3 halfVel;
};
PUPbytes(Particle);

#include "leanmd.decl.h"

extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CProxy_Cell cellArray;
extern /* readonly */ CProxy_Compute computeArray;
extern /* readonly */ CkGroupID mCastGrpID;

extern /* readonly */ int cellArrayDimX;
extern /* readonly */ int cellArrayDimY;
extern /* readonly */ int cellArrayDimZ;
extern /* readonly */ vec3 boundaryMin;
extern /* readonly */ vec3 boundaryMax;
extern /* readonly */ vec3 domainDim;
extern /* readonly */ vec3 fluidMin;
extern /* readonly */ vec3 fluidMax;
extern /* readonly */ int finalStepCount;
extern /* readonly */ int checkptStrategy;
extern /* readonly */ std::string logs;

#endif
