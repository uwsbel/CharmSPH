#include <string>
#include "time.h"
#include "ckmulticast.h"
#include "Main.h"
#include "Cell.h"
#include "Compute.h"

/* readonly */ CProxy_Main mainProxy;
/* readonly */ CProxy_Cell cellArray;
/* readonly */ CProxy_Compute computeArray;
/* readonly */ CkGroupID mCastGrpID;

/* readonly */ int cellArrayDimX;
/* readonly */ int cellArrayDimY;
/* readonly */ int cellArrayDimZ;
/* readonly */ vec3 boundaryMin;
/* readonly */ vec3 boundaryMax;
/* readonly */ vec3 domainDim;
/* readonly */ vec3 fluidMin;
/* readonly */ vec3 fluidMax;
/* readonly */ int finalStepCount; 
/* readonly */ int firstLdbStep; 
/* readonly */ int ldbPeriod;
/* readonly */ int checkptFreq; 
/* readonly */ int checkptStrategy;
/* readonly */ std::string logs;

// Entry point of Charm++ application
Main::Main(CkArgMsg* m) {
  CkPrintf("\nLENNARD JONES MOLECULAR DYNAMICS START UP ...\n");

  const std::string out_dir("output");
  const std::string mkMyDir = std::string("mkdir ") + out_dir;
  system(mkMyDir.c_str());
  const std::string rmCmd = std::string("rm ") + out_dir + std::string("/*");
  system(rmCmd.c_str());


  //set variable values to a default set
  cellArrayDimX = CELLARRAY_DIM_X;
  cellArrayDimY = CELLARRAY_DIM_Y;
  cellArrayDimZ = CELLARRAY_DIM_Z;

  finalStepCount = DEFAULT_FINALSTEPCOUNT;
  firstLdbStep = DEFAULT_FIRST_LDB;
  ldbPeriod = DEFAULT_LDB_PERIOD;
  checkptFreq = DEFAULT_FT_PERIOD;

  mainProxy = thisProxy;

  //branch factor for spanning tree of multicast
  int bFactor = 4;
  //creating the multicast spanning tree
  mCastGrpID = CProxy_CkMulticastMgr::ckNew(bFactor);

  int numPes = CkNumPes();
  int currPe = -1, pe;
  int cur_arg = 1;

  CkPrintf("\nInput Parameters...\n");

  //read user parameters
  //number of cells in each dimension
  if (m->argc > cur_arg) {
    cellArrayDimX=atoi(m->argv[cur_arg++]);
    cellArrayDimY=atoi(m->argv[cur_arg++]);
    cellArrayDimZ=atoi(m->argv[cur_arg++]);
    CkPrintf("Cell Array Dimension X:%d Y:%d Z:%d of size %f %f %f\n",cellArrayDimX,cellArrayDimY,cellArrayDimZ,CELL_SIZE_X,CELL_SIZE_Y,CELL_SIZE_Z);
  }
  double halfH = H/2;
  boundaryMin = vec3(-halfH,-halfH,-halfH);
  boundaryMax = vec3((cellArrayDimX) * CELL_SIZE_X, (cellArrayDimY) * CELL_SIZE_Y, (cellArrayDimZ) * CELL_SIZE_Z) + halfH;
  boundaryMax.print();
  domainDim = boundaryMax - boundaryMin;
  double twoH = (2 * H);
  fluidMin = boundaryMax * 0.5;
  fluidMax = vec3(boundaryMax.x, boundaryMax.y, boundaryMax.z) - twoH;
  fluidMax.print();

  //number of steps in simulation
  if (m->argc > cur_arg) {
    finalStepCount=atoi(m->argv[cur_arg++]);
    CkPrintf("Final Step Count:%d\n",finalStepCount);
  }

  //step after which load balancing starts
  if (m->argc > cur_arg) {
    firstLdbStep=atoi(m->argv[cur_arg++]);
    CkPrintf("First LB Step:%d\n",firstLdbStep);
  }

  //periodicity of load balancing
  if (m->argc > cur_arg) {
    ldbPeriod=atoi(m->argv[cur_arg++]);
    CkPrintf("LB Period:%d\n",ldbPeriod);
  }

  //periodicity of checkpointing
  if (m->argc > cur_arg) {
    checkptFreq=atoi(m->argv[cur_arg++]);
    CkPrintf("FT Period:%d\n",checkptFreq);
  }

  checkptStrategy = 1;
  //choose the checkpointing strategy use in disk checkpointing
  if (m->argc > cur_arg) {
  	checkptStrategy = 0;
    logs = m->argv[cur_arg];
  }

  cellArray = CProxy_Cell::ckNew();
  //initializing the 3D Patch array (with a uniform distribution) and 6D compute array
  int patchCount = 0;
  float ratio = ((float)CkNumPes() - 1)/(cellArrayDimX*cellArrayDimY*cellArrayDimZ);
  computeArray = CProxy_Compute::ckNew();
  for (int x=0; x<cellArrayDimX; x++)
    for (int y=0; y<cellArrayDimY; y++)
      for (int z=0; z<cellArrayDimZ; z++) {
        cellArray(x, y, z).insert((int)(patchCount++ * ratio));
        cellArray(x, y, z).createComputes();
      }

  cellArray.doneInserting();
  CkPrintf("\nCells: %d X %d X %d .... created\n", cellArrayDimX, cellArrayDimY, cellArrayDimZ);

  delete m;
}

//constructor for chare object migration
Main::Main(CkMigrateMessage* msg): CBase_Main(msg) { 
}

//pup routine incase the main chare moves, pack important information
void Main::pup(PUP::er &p) {
  CBase_Main::pup(p);
  __sdag_pup(p);
}

#include "leanmd.def.h"
