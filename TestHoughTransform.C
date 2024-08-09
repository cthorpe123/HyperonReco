R__LOAD_LIBRARY(libVFitter.so);
R__LOAD_LIBRARY($HYP_TOP/lib/libHyperon.so);
R__LOAD_LIBRARY($HYP_TOP/lib/libParticleDict.so);

#include "EventAssembler.h"
#include "HoughTransformer.h"

void TestHoughTransform(){

  const std::string filename = "HyperonTrees.root";
  const int run_to_view = 7003;
  const int subrun_to_view = 1035;
  const int event_to_view = 51797;

  EventAssembler E(false);
  E.SetFile(filename,"Background");

  // Select the event we want to analyse
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
  std::cout << "Event has " << e.TracklikePrimaryDaughters.size() << " tracks and " << e.ShowerlikePrimaryDaughters.size() << " showers" << std::endl;

  // Select the PFP we want to use the hits from
  //RecoParticle pfp = e.ShowerlikePrimaryDaughters.at(0);
  RecoParticle pfp = e.TracklikePrimaryDaughters.at(1);

  double roi_size_ch = 100; 
  double roi_size_tick = roi_size_ch*TickPerWire; 

  std::vector<std::vector<hyperonreco::HitLite>> hits(3);
  hyperonreco::AddHits(hits,pfp.HitChannels,pfp.HitTicks,pfp.HitWidths,pfp.HitTrackIDs,pfp.HitPDGs);
  hyperonreco::AddHits(hits,e.UnclaimedHitChannels,e.UnclaimedHitTicks,e.UnclaimedHitWidths,e.UnclaimedHitTrackIDs,e.UnclaimedHitPDGs);
  hyperonreco::KeepHitsInROI(TVector3(pfp.X_NoSC,pfp.Y_NoSC,pfp.Z_NoSC),hits,roi_size_ch,roi_size_tick);

  std::vector<std::vector<hyperonreco::HoughTransformPoint>> clusters;
  for(int i_pl=0;i_pl<3;i_pl++){

    //std::cout << "i_pl=" << i_pl << std::endl;

    std::pair<double,double> vertex_ch_tick = WireTick(TVector3(pfp.X_NoSC,pfp.Y_NoSC,pfp.Z_NoSC),i_pl); 
    double or_ch = vertex_ch_tick.first;
    double or_ti = vertex_ch_tick.second;

    hyperonreco::HoughTransformer transformer(hits.at(i_pl),i_pl,vertex_ch_tick.first,vertex_ch_tick.second,true);
    transformer.SetEvent(run_to_view,subrun_to_view,event_to_view);   
    transformer.MakeTransform2();
    transformer.FindPeaks2();
    clusters.push_back(transformer.MakeClusters());
    std::pair<double,double> performance = transformer.GetPerformanceMetrics();
    std::cout << "Purity = " << performance.first << "  Completeness = " << performance.second << std::endl;

  }


}

