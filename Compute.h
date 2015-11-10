#ifndef __COMPUTE_H__
#define __COMPUTE_H__

#include "defs.h"

//class representing the interaction agents between a couple of cells
class Compute : public CBase_Compute {
  private:
    Compute_SDAG_CODE
    int stepCount;  //current step number
    double energy[2]; //store potential energy
    //handles to differentiate the two multicast sections I am part of
    CkSectionInfo mcast1;     
    CkSectionInfo mcast2;

  public:
    Compute();
    Compute(CkMigrateMessage *msg);
    void pup(PUP::er &p);

    void selfInteractSPH(ParticleDataMsg *msg);
    void interactSPH(ParticleDataMsg *msg1, ParticleDataMsg *msg2);

    void selfInteract(ParticleDataMsg *msg);
    void interact(ParticleDataMsg *msg1, ParticleDataMsg *msg2);
};

#endif
