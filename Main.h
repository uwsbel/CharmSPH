#ifndef __MAIN_H__
#define __MAIN_H__

#include <string>
#include "defs.h"


//central controller chare
class Main : public CBase_Main {
  private:
    Main_SDAG_CODE
    double startBenchmarkTime;

  public:
    Main(CkArgMsg* msg);
    Main(CkMigrateMessage* msg);
    std::string getSimulationID();
    void compileOutput();
    void writeSimParams();
    void writeTimingResults(double totalTime);
    void setDimensions();
    void setDefaultParams();
    void printParams();
    void initOutDirs(std::string simID);
    void pup(PUP::er &p);
};
#endif
