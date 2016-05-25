/* C++ includes */
#include <iostream>
#include <sstream>
#include <fstream>
/* Charm++ includes */
#include "time.h"
#include "ckmulticast.h"
/* CharmSPH includes*/
#include "Main.h"
#include "Cell.h"
#include "Compute.h"


/* Charm++ Globals */
/* readonly */ CProxy_Main mainProxy;
/* readonly */ CProxy_Cell cellArray;
/* readonly */ CProxy_Compute computeArray;
/* readonly */ CkGroupID mCastGrpID;
/* readonly */ int numCells;
/* readonly */ int numComputes;

/* SPH Globals */
/* readonly */ double h;
/* readonly */ double dt;
/* readonly */ double maxVel;
/* readonly */ double particleMass;
/* readonly */ double cutOffDist;
/* readonly */ int cellSizeMult;
/* readonly */ int markDistMult;
/* readonly */ int numFluidMarkers;
/* readonly */ int numBoundaryMarkers;
/* readonly */ int writeAll;
/* readonly */ int writePeriod;
/* readonly */ int writeBoundary;
/* readonly */ std::string simID;
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
  mainProxy = thisProxy;
  setDefaultParams(); // Set default params. Keep these if args not set
  numCells = 0;
  numComputes = 0;
  numFluidMarkers = 0;
  numBoundaryMarkers = 0; 
  
  /**
   *  Parse input arguments
   *    * -x = x dimension
   *    * -y = y dimension
   *    * -z = z dimension
   *    * -t = t is the total number of time steps
   *    * -h = h is the particle interaction radius (cutoff radius)
   *    * -dt = dt delta t at every time step
   *    * -mv = mv is the estimate of the maximum velocity of the particles in the model.
   *    * -w = Write fluid.
   *    * -wp = Write period. After every wp steps we write output
   *    * -wb = Write boundary. 1 if you want to write the boundary, 0 if not.
   *    * -csm = Cell Size Multiplier. How much should we multiply
   *    * -mdm = Marker Distance Multiplier. How much should we multiply
   *    * -lbp = Load balance period
   */
  for(int i = 1;i < m->argc;i++){
    std::cout << "arg at " << i << " " << m->argv[i] << std::endl;

    if (i + 1 != m->argc){ // Check that we haven't finished parsing already
        if (strcmp(m->argv[i], "-x") == 0) { 
          domainDim.x = atof(m->argv[i + 1]);
        } 
        else if (strcmp(m->argv[i], "-y") == 0) {
          domainDim.y = atof(m->argv[i + 1]);
        } 
        else if (strcmp(m->argv[i], "-z") == 0) {
          domainDim.z = atof(m->argv[i + 1]);
        } 
        else if (strcmp(m->argv[i], "-t") == 0) {
          finalStepCount = atoi(m->argv[i + 1]);
        }
        else if (strcmp(m->argv[i], "-h") == 0) {
          h = atof(m->argv[i + 1]);
        }
        else if (strcmp(m->argv[i], "-dt") == 0) {
          dt = atof(m->argv[i + 1]);
        }
        else if (strcmp(m->argv[i], "-mv") == 0) {
          maxVel = atof(m->argv[i + 1]);
        }
        else if (strcmp(m->argv[i], "-w") == 0) {
          writeAll = atoi(m->argv[i + 1]);
        }
        else if (strcmp(m->argv[i], "-wp") == 0) {
          writePeriod = atoi(m->argv[i + 1]);
        }
        else if (strcmp(m->argv[i], "-wb") == 0) {
          writeBoundary = atoi(m->argv[i + 1]);
        }
        else if (strcmp(m->argv[i], "-csm") == 0){
          cellSizeMult = atof(m->argv[i + 1]);
        }
        else if (strcmp(m->argv[i], "-mdm") == 0){
          markDistMult = atof(m->argv[i + 1]);
        }
        else if (strcmp(m->argv[i], "-lbp") == 0) {
          ldbPeriod = atof(m->argv[i + 1]);
        }
    }
  }
  particleMass = h * h * h * RHO0;
  cutOffDist = 2 * h;
  setDimensions();
  gravity = vec3(0, GRAVITY, 0);
  simID = getSimulationID();
  initOutDirs(simID); // Initialize the fluid and boundary dirs where outputs goes

  int numPes = CkNumPes();
  int currPe = -1, pe;
  int cur_arg = 1;

  int bFactor = 4; //branch factor for spanning tree of multicast
  mCastGrpID = CProxy_CkMulticastMgr::ckNew(bFactor); //creating the multicast spanning tree
  
  cellArray = CProxy_Cell::ckNew();
  computeArray = CProxy_Compute::ckNew();
  //initializing the 3D Patch array (with a uniform distribution) and 6D compute array
  int patchCount = 0;
  float ratio = ((float)CkNumPes() - 1)/(cellArrayDim.x * cellArrayDim.y * cellArrayDim.z);
  for (int x=0; x<cellArrayDim.x; x++) {
    for (int y=0; y<cellArrayDim.y; y++) {
      for (int z=0; z<cellArrayDim.z; z++) {
        cellArray(x, y, z).insert((int)(patchCount++ * ratio));
        cellArray(x, y, z).createComputes();
        numCells++;
      }
    }
  }
  cellArray.doneInserting();
  delete m;
}

//constructor for chare object migration
Main::Main(CkMigrateMessage* msg): CBase_Main(msg) {}

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
  markDistMult = DEFAULT_MARKDISTMULT;
  cellSizeMult = DEFAULT_CELLSIZEMULT;

  writeAll = DEFAULT_WRITE;
  writePeriod = DEFAULT_WRITEPERIOD;
  writeBoundary = DEFAULT_WRITEBOUNDARY;
  domainMin = vec3(DEFAULT_MIN_X, DEFAULT_MIN_Y, DEFAULT_MIN_Z);
  domainMax = vec3(DEFAULT_MAX_X, DEFAULT_MAX_Y, DEFAULT_MAX_Z);
  domainDim = vec3(DEFAULT_MAX_X, DEFAULT_MAX_Y, DEFAULT_MAX_Z);
  fluidMin = vec3(DEFAULT_FLUIDMIN_X, DEFAULT_FLUIDMIN_Y, DEFAULT_FLUIDMIN_Z);
  fluidMax = vec3(DEFAULT_FLUIDMAX_X, DEFAULT_FLUIDMAX_Y, DEFAULT_FLUIDMAX_Z);
}

void Main::setDimensions()
{

  //number of cells in each dimension
  domainMin = vec3(0, 0, 0);
  domainMax = domainDim;

  cellSize.x = cellSizeMult * h;
  cellSize.y = cellSizeMult * h;
  cellSize.z = cellSizeMult * h;

  domainDim = domainMax - domainMin;

  //set variable values to a default set
  cellArrayDim.x = floor(domainDim.x / cellSize.x);
  cellArrayDim.y = floor(domainDim.y / cellSize.y);
  cellArrayDim.z = floor(domainDim.z / cellSize.z);

  cellSize.x = domainDim.x / cellArrayDim.x;
  cellSize.y = domainDim.y / cellArrayDim.y;
  cellSize.z = domainDim.z / cellArrayDim.z;

  // tune particle spacing based on cell size
  mDist = vec3(1, 1, 1) * markDistMult * h;

  int cNX = ceil(cellSize.x / mDist.x);
  int cNY = ceil(cellSize.y / mDist.y);
  int cNZ = ceil(cellSize.z / mDist.z);

  mDist.x = cellSize.x / cNX;
  mDist.y = cellSize.y / cNY;
  mDist.z = cellSize.z / cNZ;

  CkPrintf("mDist: \n");
  mDist.print();
}

/**
 * @brief getSimulationID
 * @details Put together the simulation ID using some of the main simulation parameters.
 *               charmsph_h_CellSizeMult_numCores_dt_t.json
 * @return simID
 */
std::string Main::getSimulationID()
{
  std::stringstream ssSimID;
  ssSimID << h << "_";
  ssSimID << cellSize.x / h << "_";
  ssSimID << CkNumPes() << "_";
  ssSimID << dt << "_";
  ssSimID << finalStepCount << "_";
  ssSimID << domainDim.x << "-" << domainDim.y << "-" << domainDim.z;

  return ssSimID.str();
}

void Main::printParams()
{
  std::cout << "************** SIMULATION PARAMETERS **************" << std::endl;
  std::cout << "Simulation ID = " << getSimulationID() << std::endl;
  std::cout << "dt = " << dt << std::endl;
  std::cout << "h = " << h << std::endl;
  std::cout << "maxVel = " << maxVel << std::endl;
  std::cout << "particleMass = " << particleMass << std::endl;
  std::cout << "cutOffDist = " << cutOffDist << std::endl;
  std::cout << "ldbPeriod = " << ldbPeriod << std::endl;
  std::cout << "write = " << writeAll << std::endl;
  std::cout << "write period = " << writePeriod << std::endl;
  std::cout << "write boundary = " << writeBoundary << std::endl;
  std::cout << "numCellChares = " << numCells << std::endl;
  std::cout << "numComputeChares = " << numComputes << std::endl;
  std::cout << "numFluidMarkers = " << numFluidMarkers << std::endl;
  std::cout << "numBoundaryMarkers = " << numBoundaryMarkers << std::endl;
  std::cout << "numPes = " << CkNumPes() << std::endl;
  std::cout << "cellArrayDim = "; cellArrayDim.print();
  std::cout << "gravity = "; gravity.print();
  std::cout << "cellSize = "; cellSize.print();
  std::cout << "domainMin = "; domainMin.print();
  std::cout << "domainMax = "; domainMax.print();
  std::cout << "fluidMin = "; fluidMin.print();
  std::cout << "fluidMax = "; fluidMax.print();
  std::cout << "***************************************************" << std::endl;
}

void Main::compileOutput()
{
  std::string cmd = "python compileOutput.py " + simID;
  system(cmd.c_str());
}

void Main::writeSimParams(int flag)
{
  std::ofstream simParamsFile;
  std::string filename;
  if(flag == 1){
    filename = "output/" + simID + "/SimParams1.json";
  }
  else{
    filename = "output/" + simID + "/SimParams2.json";
  }
  std::stringstream ssSimParams;

  ssSimParams << "{" << std::endl;
  ssSimParams << "\"simID\": " << "\"" << simID << "\""  << "," << std::endl;
  ssSimParams << "\"h\": " << h << "," << std::endl;
  ssSimParams << "\"dt\": " << dt << "," << std::endl;
  ssSimParams << "\"maxVel\": " << maxVel << "," << std::endl;
  ssSimParams << "\"particleMass\": " << particleMass << "," << std::endl;
  ssSimParams << "\"cutOffDist\": " << cutOffDist << "," << std::endl;
  ssSimParams << "\"numPes\": " << CkNumPes() << "," << std::endl;
  ssSimParams << "\"numFluidMarkers\": " << numFluidMarkers << "," << std::endl;
  ssSimParams << "\"numBoundaryMarkers\": " << numBoundaryMarkers << "," << std::endl;
  ssSimParams << "\"numCellChares\": " << numCells << "," << std::endl;
  ssSimParams << "\"numComputeChares\": " << numComputes << "," << std::endl;
  ssSimParams << "\"write\": " << writeAll << "," << std::endl;
  ssSimParams << "\"writePeriod\": " << writePeriod << "," << std::endl;
  ssSimParams << "\"writeBoundary\": " << writeBoundary << "," << std::endl;
  ssSimParams << "\"domainDim\": [" << domainDim.x << "," << domainDim.y << "," << domainDim.z << "]" << "," << std::endl;
  ssSimParams << "\"cellSize\": [" << cellSize.x << "," << cellSize.y << "," << cellSize.z << "]" << "," << std::endl;
  ssSimParams << "\"cellArrayDim\": [" << cellArrayDim.x << "," << cellArrayDim.y << "," << cellArrayDim.z << "]" << std::endl;


  ssSimParams << "}" << std::endl;

  simParamsFile.open(filename.c_str());
  simParamsFile << ssSimParams.str();
  simParamsFile.close();
}

/**
 * @brief initOutDirs
 * @details mkdir output directories if the don't exist and delete old output from them if any.
 */
void Main::initOutDirs(std::string simID)
{
  // Clear output directory
  const std::string outDir(" output");
  const std::string simDir(" output/" + simID);
  const std::string fluidDir(simDir + "/fluid");
  const std::string boundaryDir(simDir + "/boundary");
  const std::string mkOutputDir = std::string("mkdir ") + outDir + simDir + fluidDir + boundaryDir;
  system(mkOutputDir.c_str());
  const std::string rmOutCmd = std::string("rm ") + outDir + std::string("/*");
  const std::string rmSimCmd = std::string("rm ") + simDir + std::string("/*");
  const std::string rmFluidCmd = std::string("rm ") + fluidDir + std::string("/*");
  const std::string rmBoundaryCmd = std::string("rm ") + boundaryDir + std::string("/*");
  system(rmOutCmd.c_str());
  system(rmSimCmd.c_str());
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
