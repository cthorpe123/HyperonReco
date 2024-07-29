R__LOAD_LIBRARY(libVFitter.so);
R__LOAD_LIBRARY($HYP_TOP/lib/libHyperon.so);
R__LOAD_LIBRARY($HYP_TOP/lib/libParticleDict.so);

#include "VFitter.h"
#include "SpacePointVisualisation.h" 
#include "EventAssembler.h"


void DrawEvent(const std::string filename,const int run_to_view,const int subrun_to_view,const int event_to_view){

  EventAssembler E(false);
  E.SetFile(filename,"Background");

  int ievent=0;
  while(ievent<E.GetNEvents()){
   Event e = E.GetEvent(ievent);
    if(e.run == run_to_view && e.subrun == subrun_to_view && e.event == event_to_view)
      break;
      ievent++;
  }
   
  if(ievent == E.GetNEvents()){
    std::cout << "Event " << run_to_view << " " << subrun_to_view << " " << event_to_view << " not found, exiting" << std::endl;
    return;
  }

  Event e = E.GetEvent(ievent);
  //RecoParticle pfp = e.ShowerlikePrimaryDaughters.at(0);
  //hyperonreco::FittedV v = hyperonreco::MakeFittedVGuessShower(pfp);
  RecoParticle pfp = e.TracklikePrimaryDaughters.at(1);
  hyperonreco::FittedV v = hyperonreco::MakeFittedVGuessTrack(pfp);

  double roi_size_ch = 100; 
  double roi_size_tick = roi_size_ch*TickPerWire; 

  hyperonreco::VFitter fitter;
  fitter.SetROI(roi_size_ch,roi_size_tick,TVector3(pfp.X_NoSC,pfp.Y_NoSC,pfp.Z_NoSC));
  fitter.SetActivePlanes({0});
  fitter.SetOutlierCut(5); 
  fitter.SetFitTune(1.0); 
 
  fitter.AddData({pfp.HitChannels_Plane0,pfp.HitChannels_Plane1,pfp.HitChannels_Plane2},{pfp.HitTicks_Plane0,pfp.HitTicks_Plane1,pfp.HitTicks_Plane2},{pfp.HitWidths_Plane0,pfp.HitWidths_Plane1,pfp.HitWidths_Plane2});
  fitter.AddData({e.UnclaimedHitChannels_Plane0,e.UnclaimedHitChannels_Plane1,e.UnclaimedHitChannels_Plane2},{e.UnclaimedHitTicks_Plane0,e.UnclaimedHitTicks_Plane1,e.UnclaimedHitTicks_Plane2},{e.UnclaimedHitWidths_Plane0,e.UnclaimedHitWidths_Plane1,e.UnclaimedHitWidths_Plane2});

  fitter.SetGuess(v);
  fitter.DoFitGridSearch(v,1000000); 
  //fitter.DoFit(v); 

  std::cout << "Arm 1 Points: " << v.Arm1Points << "  Arm 2 Points: " << v.Arm2Points << std::endl;
/*
  for(int i=0;i<20;i++){
    roi_size_ch *= 1.02;
    roi_size_tick *= 1.02;
    fitter.RemoveOutliers(); 
    fitter.SetGuess(v);
    fitter.SetROI(roi_size_ch,roi_size_tick,TVector3(pfp.X_NoSC,pfp.Y_NoSC,pfp.Z_NoSC));
    fitter.DoFit(v); 
    std::cout << "Arm 1 Points: " << v.Arm1Points << "  Arm 2 Points: " << v.Arm2Points << std::endl;
  }
*/
  hyperonreco::KeepHitsInROI(TVector3(pfp.X,pfp.Y,pfp.Z),pfp.HitChannels_Plane0,pfp.HitTicks_Plane0,pfp.HitWidths_Plane0,roi_size_ch,roi_size_tick,0);
  hyperonreco::KeepHitsInROI(TVector3(pfp.X,pfp.Y,pfp.Z),pfp.HitChannels_Plane1,pfp.HitTicks_Plane1,pfp.HitWidths_Plane1,roi_size_ch,roi_size_tick,1);
  hyperonreco::KeepHitsInROI(TVector3(pfp.X,pfp.Y,pfp.Z),pfp.HitChannels_Plane2,pfp.HitTicks_Plane2,pfp.HitWidths_Plane2,roi_size_ch,roi_size_tick,2);

  std::vector<std::vector<TGraph*>> graph_2d_vec(3);
  graph_2d_vec.at(0).push_back(hyperonreco::Make2DGraph(pfp.HitChannels_Plane0,pfp.HitTicks_Plane0,pfp.HitWidths_Plane0,1,1));
  graph_2d_vec.at(1).push_back(hyperonreco::Make2DGraph(pfp.HitChannels_Plane1,pfp.HitTicks_Plane1,pfp.HitWidths_Plane1,1,1));
  graph_2d_vec.at(2).push_back(hyperonreco::Make2DGraph(pfp.HitChannels_Plane2,pfp.HitTicks_Plane2,pfp.HitWidths_Plane2,1,1));

  hyperonreco::KeepHitsInROI(TVector3(pfp.X,pfp.Y,pfp.Z),e.UnclaimedHitChannels_Plane0,e.UnclaimedHitTicks_Plane0,e.UnclaimedHitWidths_Plane0,roi_size_ch,roi_size_tick,0);
  hyperonreco::KeepHitsInROI(TVector3(pfp.X,pfp.Y,pfp.Z),e.UnclaimedHitChannels_Plane1,e.UnclaimedHitTicks_Plane1,e.UnclaimedHitWidths_Plane1,roi_size_ch,roi_size_tick,1);
  hyperonreco::KeepHitsInROI(TVector3(pfp.X,pfp.Y,pfp.Z),e.UnclaimedHitChannels_Plane2,e.UnclaimedHitTicks_Plane2,e.UnclaimedHitWidths_Plane2,roi_size_ch,roi_size_tick,2);

  graph_2d_vec.at(0).push_back(hyperonreco::Make2DGraph(e.UnclaimedHitChannels_Plane0,e.UnclaimedHitTicks_Plane0,e.UnclaimedHitWidths_Plane0,2,2));
  graph_2d_vec.at(1).push_back(hyperonreco::Make2DGraph(e.UnclaimedHitChannels_Plane1,e.UnclaimedHitTicks_Plane1,e.UnclaimedHitWidths_Plane1,2,2));
  graph_2d_vec.at(2).push_back(hyperonreco::Make2DGraph(e.UnclaimedHitChannels_Plane2,e.UnclaimedHitTicks_Plane2,e.UnclaimedHitWidths_Plane2,2,2));

  hyperonreco::Draw2DGraphs(graph_2d_vec,&v);

  
 



}
