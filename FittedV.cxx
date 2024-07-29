#ifndef _FittedV_cxx_
#define _FittedV_cxx_

#include "FittedV.h"

namespace hyperonreco {

TVector3 FittedV::GetArm1Dir() const {
  return MakeTVector3(Arm1Theta,Arm1Phi); 
} 

TVector3 FittedV::GetArm2Dir() const {
  return MakeTVector3(Arm2Theta,Arm2Phi); 
} 

LineWireTick FittedV::GetArm1_2D(int i_pl) const {
  return ProjectXYZWireTick(Vertex,GetArm1Dir(),i_pl);
}

LineWireTick FittedV::GetArm2_2D(int i_pl) const {
  return ProjectXYZWireTick(Vertex,GetArm2Dir(),i_pl);
}

LineWireTick2 FittedV::GetArm1_2D2(int i_pl) const {
  return ProjectXYZWireTick2(Vertex,GetArm1Dir(),i_pl);
}

LineWireTick2 FittedV::GetArm2_2D2(int i_pl) const {
  return ProjectXYZWireTick2(Vertex,GetArm2Dir(),i_pl);
}

}

#endif

