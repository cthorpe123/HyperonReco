R__LOAD_LIBRARY(libVFitter.so);
R__LOAD_LIBRARY($HYP_TOP/lib/libHyperon.so);
R__LOAD_LIBRARY($HYP_TOP/lib/libParticleDict.so);

#include "VFitter.h"
#include "SpacePointVisualisation.h" 
#include "EventAssembler.h"


void DrawEvent(){

  const std::string filename = "HyperonTrees_Background.root";
  const int run_to_view = 6562;
  const int subrun_to_view = 102;
  const int event_to_view = 5111;

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
  RecoParticle pfp = e.TracklikePrimaryDaughters.at(0);
  hyperonreco::FittedV v = hyperonreco::MakeFittedVGuessTrack(pfp);

  double roi_size_ch = 50; 
  double roi_size_tick = roi_size_ch*TickPerWire; 

  std::vector<std::vector<hyperonreco::HitLite>> hits(3);
  hyperonreco::AddHits(hits,pfp.HitChannels,pfp.HitTicks,pfp.HitWidths,pfp.HitTrackIDs,pfp.HitPDGs);
  //hyperonreco::AddHits(hits,e.UnclaimedHitChannels,e.UnclaimedHitTicks,e.UnclaimedHitWidths,e.UnclaimedHitTrackIDs,e.UnclaimedHitPDGs);
  hyperonreco::KeepHitsInROI(TVector3(pfp.X_NoSC,pfp.Y_NoSC,pfp.Z_NoSC),hits,roi_size_ch,roi_size_tick);

  hyperonreco::VFitter fitter(true);
  fitter.SetEvent(run_to_view,subrun_to_view,event_to_view);   
  fitter.SetROI(roi_size_ch,roi_size_tick,TVector3(pfp.X_NoSC,pfp.Y_NoSC,pfp.Z_NoSC));
  fitter.AddData(hits);
  fitter.SetGuess(v);
  fitter.DoFitGridSearch3(v,50000); 

}
