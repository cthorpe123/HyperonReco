#ifndef _FittedV_cxx_
#define _FittedV_cxx_

#include "FittedV.h"

namespace hyperonreco {

TVector3 FittedV::GetArm1Dir() const {
  return MakeTVector3(Arm1Theta,Arm1Phi); 
} 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

TVector3 FittedV::GetArm2Dir() const {
  return MakeTVector3(Arm2Theta,Arm2Phi); 
} 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

LineWireTick FittedV::GetArm1_2D(int i_pl) const {
  return ProjectXYZWireTick(Vertex,GetArm1Dir(),i_pl);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

LineWireTick FittedV::GetArm2_2D(int i_pl) const {
  return ProjectXYZWireTick(Vertex,GetArm2Dir(),i_pl);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

LineWireTick2 FittedV::GetArm1_2D2(int i_pl) const {
  return ProjectXYZWireTick2(Vertex,GetArm1Dir(),i_pl);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

LineWireTick2 FittedV::GetArm2_2D2(int i_pl) const {
  return ProjectXYZWireTick2(Vertex,GetArm2Dir(),i_pl);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::pair<TGraph*,TGraph*> FittedV::Make2DVGraphs(int i_pl) const {

   TGraph* g1 = new TGraph(2);
   TGraph* g2 = new TGraph(2);
   g1->SetName("fit_graph2D_1");
   g2->SetName("fit_graph2D_2");
   LineWireTick2 line_1 = GetArm1_2D2(i_pl);
   LineWireTick2 line_2 = GetArm2_2D2(i_pl);

   g1->SetPoint(0,line_1.Start.X(),line_1.Start.Y());
   g2->SetPoint(0,line_2.Start.X(),line_2.Start.Y());

   g1->SetPoint(1,line_1.Start.X()+line_1.Direction.X()*500,line_1.Start.Y()+line_1.Direction.Y()*500);
   g2->SetPoint(1,line_2.Start.X()+line_2.Direction.X()*500,line_2.Start.Y()+line_2.Direction.Y()*500);

   g1->SetLineColor(2); 
   g2->SetLineColor(3); 

   return std::make_pair(g1,g2);
  
}

}

#endif

