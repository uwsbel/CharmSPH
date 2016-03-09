#include "defs.h"
#include "Cell.h"
#include "Compute.h"
#include "physics.h"
#include "ckmulticast.h"
#include <algorithm>
using std::swap;

//compute - Default constructor
Compute::Compute() : stepCount(1) {
  energy[0] = energy[1] = 0;
  usesAtSync = true;
}

Compute::Compute(CkMigrateMessage *msg): CBase_Compute(msg)  { 
  usesAtSync = true;
  delete msg;
}

// interaction within cell. This function should be very similat to selfInteract.
void Compute::selfInteractSPH(ParticleDataMsg *msg){
  double energyP = 0;
  std::vector<vec4> dVel_dRho;

  calcInternalForcesSPH(msg, stepCount, dVel_dRho);

  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  CkGetSectionInfo(mcast1, msg);

  mCastGrp->contribute(sizeof(vec4) * (msg->lengthAll), &dVel_dRho[0], CkReduction::sum_double, mcast1);

  delete msg;
}

//interaction between two cells
void Compute::interactSPH(ParticleDataMsg *msg1, ParticleDataMsg *msg2){
  CkSectionInfo *handleA = &mcast1;
  CkSectionInfo *handleB = &mcast2;
  if (msg2->x * cellArrayDim.y * cellArrayDim.z + msg2->y * cellArrayDim.z + msg2->z <
      msg1->x * cellArrayDim.y * cellArrayDim.z + msg1->y * cellArrayDim.z + msg1->z)
  {
    swap(handleA, handleB);
  }

  std::vector<vec4> dVel_dRho1, dVel_dRho2;
  calcPairForcesSPH(msg1, msg2, stepCount, dVel_dRho1, dVel_dRho2);

  //contribute to force reduction
  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  CkGetSectionInfo(*handleA, msg1);
  mCastGrp->contribute(sizeof(vec4)*msg1->lengthAll, &dVel_dRho1[0], CkReduction::sum_double, *handleA);
  CkGetSectionInfo(*handleB, msg2);
  mCastGrp->contribute(sizeof(vec4)*msg2->lengthAll, &dVel_dRho2[0], CkReduction::sum_double, *handleB);

  delete msg1;
  delete msg2;
}

//pack important information if I am moving
void Compute::pup(PUP::er &p) {
  CBase_Compute::pup(p);
  __sdag_pup(p);
  p | stepCount;
  p | mcast1;
  p | mcast2;
  PUParray(p, energy, 2);
  if (p.isUnpacking() && CkInRestarting()) {
    mcast1.set_redNo(0); 
	if(mcast2.info.type != 0)
		mcast2.set_redNo(0);
  }
}
