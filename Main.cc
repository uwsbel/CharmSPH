#include <string>
#include "time.h"
#include "ckmulticast.h"
#include "Main.h"
#include "Cell.h"
#include "Compute.h"
#include <iostream>

/* Charm++ Globals */
/* readonly */ CProxy_Main mainProxy;
/* readonly */ CProxy_Cell cellArray;
/* readonly */ CProxy_Compute computeArray;
/* readonly */ CkGroupID mCastGrpID;

/* SPH Globals */
/* readonly */ double h;
/* readonly */ double dt;
/* readonly */ double maxVel;
/* readonly */ double particleMass;
/* readonly */ int3 cellArrayDim;
/* readonly */ vec3 domainMin;
/* readonly */ vec3 domainMax;
/* readonly */ vec3 domainDim;
/* readonly */ vec3 fluidMin;
/* readonly */ vec3 fluidMax;
/* readonly */ vec3 cellSize;
/* readonly */ vec3 mDist;

/* Charm++ RTS Globals */
/* readonly */ int finalStepCount; 
/* readonly */ int firstLdbStep; 
/* readonly */ int ldbPeriod;
/* readonly */ int checkptFreq; 
/* readonly */ int checkptStrategy;
/* readonly */ std::string logs;   
/* readonly */ vec3 gravity;   

// Entry point of Charm++ application
Main::Main(CkArgMsg* m) 
{
  CkPrintf("\nLENNARD JONES MOLECULAR DYNAMICS START UP ...\n");

  initOutDirs(); // Initialize the fluid and boundary dirs where outputs goes
  setDefaultParams();

  for(int i = 0;i < m->argc;i++){
    std::cout << m->argv[i] << std::endl;
  }

  mainProxy = thisProxy;

  //branch factor for spanning tree of multicast
  int bFactor = 4;
  //creating the multicast spanning tree
  mCastGrpID = CProxy_CkMulticastMgr::ckNew(bFactor);

  int numPes = CkNumPes();
  int currPe = -1, pe;
  int cur_arg = 1;

  //read user parameters
  //number of cells in each dimension

  if (m->argc > cur_arg) 
  {
    // assume domainMin is (0, 0, 0)
    domainDim.x=atof(m->argv[cur_arg++]);
    domainDim.y=atof(m->argv[cur_arg++]);
    domainDim.z=atof(m->argv[cur_arg++]);
  }
  domainMin = vec3(0, 0, 0);
  domainMax = domainDim;

  cellSize.x = 4 * h;
  cellSize.y = 4 * h;
  cellSize.z = 4 * h;

  domainDim = domainMax - domainMin;

  //set variable values to a default set
  cellArrayDim.x = floor(domainDim.x / cellSize.x);
  cellArrayDim.y = floor(domainDim.y / cellSize.y);
  cellArrayDim.z = floor(domainDim.z / cellSize.z);

  cellSize.x = domainDim.x / cellArrayDim.x;
  cellSize.y = domainDim.y / cellArrayDim.y;
  cellSize.z = domainDim.z / cellArrayDim.z;

  // tune particle spacing based on cell size
  mDist = vec3(1, 1, 1) * MarkDistMult * h;

  int cNX = ceil(cellSize.x / mDist.x);
  int cNY = ceil(cellSize.y / mDist.y);
  int cNZ = ceil(cellSize.z / mDist.z);

  mDist.x = cellSize.x / cNX;
  mDist.y = cellSize.y / cNY;
  mDist.z = cellSize.z / cNZ;

  CkPrintf("mDist: \n");
  mDist.print();

  gravity = vec3(0, GRAVITY, 0);

  CkPrintf("\nInput Parameters...\n");
  
  //number of steps in simulation
  if (m->argc > cur_arg) 
  {
    finalStepCount=atoi(m->argv[cur_arg++]);
    CkPrintf("Final Step Count:%d\n",finalStepCount);
  }

  //step after which load balancing starts
  if (m->argc > cur_arg) 
  {
    firstLdbStep=atoi(m->argv[cur_arg++]);
    CkPrintf("First LB Step:%d\n",firstLdbStep);
  }

  //periodicity of load balancing
  if (m->argc > cur_arg) 
  {
    ldbPeriod=atoi(m->argv[cur_arg++]);
    CkPrintf("LB Period:%d\n",ldbPeriod);
  }

  //periodicity of checkpointing
  if (m->argc > cur_arg) 
  {
    checkptFreq=atoi(m->argv[cur_arg++]);
    CkPrintf("FT Period:%d\n",checkptFreq);
  }

  checkptStrategy = 1;
  //choose the checkpointing strategy use in disk checkpointing
  if (m->argc > cur_arg) 
  {
    checkptStrategy = 0;
    logs = m->argv[cur_arg];
  }

  printParams();

  cellArray = CProxy_Cell::ckNew();
  //initializing the 3D Patch array (with a uniform distribution) and 6D compute array
  int patchCount = 0;
  float ratio = ((float)CkNumPes() - 1)/(cellArrayDim.x*cellArrayDim.y*cellArrayDim.z);
  computeArray = CProxy_Compute::ckNew();

  for (int x=0; x<cellArrayDim.x; x++) {
    for (int y=0; y<cellArrayDim.y; y++) {
      for (int z=0; z<cellArrayDim.z; z++) {
        cellArray(x, y, z).insert((int)(patchCount++ * ratio));
        cellArray(x, y, z).createComputes();
      }
    }
  }
  cellArray.doneInserting();
  delete m;
}

//constructor for chare object migration
Main::Main(CkMigrateMessage* msg): CBase_Main(msg) 
{
}

void Main::setDefaultParams()
{
  /* Charm++ Default params */  
  finalStepCount = DEFAULT_FINALSTEPCOUNT;
  firstLdbStep = DEFAULT_FIRST_LDB;
  ldbPeriod = DEFAULT_LDB_PERIOD;
  checkptFreq = DEFAULT_FT_PERIOD;

  /* SPH Default params */  
  h = DEFAULT_H;
  dt = DEFAULT_DT;
  maxVel = DEFAULT_MAXVEL;
  particleMass = h * h * h * RHO0;
  domainMin = vec3(DEFAULT_MIN_X, DEFAULT_MIN_Y, DEFAULT_MIN_Z);
  domainMax = vec3(DEFAULT_MAX_X, DEFAULT_MAX_Y, DEFAULT_MAX_Z);
  domainDim = vec3(DEFAULT_MAX_X, DEFAULT_MAX_Y, DEFAULT_MAX_Z);
  fluidMin = vec3(DEFAULT_FLUIDMIN_X, DEFAULT_FLUIDMIN_Y, DEFAULT_FLUIDMIN_Z);
  fluidMax = vec3(DEFAULT_FLUIDMAX_X, DEFAULT_FLUIDMAX_Y, DEFAULT_FLUIDMAX_Z);
}

void Main::printParams()
{
  std::cout << "dt = " << dt << std::endl;
  std::cout << "h = " << h << std::endl;
  std::cout << "maxVel = " << maxVel << std::endl;
  std::cout << "particleMass = " << particleMass << std::endl;
  std::cout << "cutOffDist = " << PTP_CUT_OFF << std::endl;
  std::cout << "gravity = "; gravity.print();
  std::cout << "domainMin = "; domainMin.print();
  std::cout << "domainMax = "; domainMax.print();
  std::cout << "fluidMin = "; fluidMin.print();
  std::cout << "fluidMax = "; fluidMax.print();
}

/**
 * @brief initOutDirs
 * @details mkdir output directories if the don't exist and delete old output from them if any.
 */
void Main::initOutDirs()
{
  // Clear output directory
  const std::string outDir(" output");
  const std::string fluidDir(" output/fluid");
  const std::string boundaryDir(" output/boundary");
  const std::string mkOutputDir = std::string("mkdir ") + outDir + fluidDir + boundaryDir;
  system(mkOutputDir.c_str());
  const std::string rmOutCmd = std::string("rm ") + outDir + std::string("/*");
  const std::string rmFluidCmd = std::string("rm ") + fluidDir + std::string("/*");
  const std::string rmBoundaryCmd = std::string("rm ") + boundaryDir + std::string("/*");
  system(rmOutCmd.c_str());
  system(rmFluidCmd.c_str());
  system(rmBoundaryCmd.c_str());
}

//pup routine incase the main chare moves, pack important information
void Main::pup(PUP::er &p) 
{
  CBase_Main::pup(p);
  __sdag_pup(p);
}

#include "charmsph.def.h"
