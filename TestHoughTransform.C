R__LOAD_LIBRARY(lib/libVFitter.so);
R__LOAD_LIBRARY($HYP_TOP/lib/libHyperon.so);
R__LOAD_LIBRARY($HYP_TOP/lib/libParticleDict.so);

#include "EventAssembler.h"
#include "HoughTransformer.h"

void TestHoughTransform(){

  const std::string filename = "GoodReco.root";

  EventAssembler E(false);
  E.SetFile(filename,"Background");

  // Select the event we want to analyse
  int ievent=0;
  while(ievent<E.GetNEvents()){

    Event e = E.GetEvent(ievent);
    ievent++;

    std::cout << "Event has " << e.TracklikePrimaryDaughters.size() << " tracks and " << e.ShowerlikePrimaryDaughters.size() << " showers" << std::endl;

    int pfp_ctr=0;
    for(RecoParticle pfp : e.ShowerlikePrimaryDaughters){

      // Select the PFP we want to use the hits from
      //RecoParticle pfp = e.TracklikePrimaryDaughters.at(1);

      double roi_size_ch = 100; 
      double roi_size_tick = roi_size_ch*TickPerWire; 

      std::vector<std::vector<hyperonreco::HitLite>> hits(3);
      hyperonreco::AddHits(hits,pfp.HitChannels,pfp.HitTicks,pfp.HitWidths,pfp.HitTrackIDs,pfp.HitPDGs);
      hyperonreco::AddHits(hits,e.UnclaimedHitChannels,e.UnclaimedHitTicks,e.UnclaimedHitWidths,e.UnclaimedHitTrackIDs,e.UnclaimedHitPDGs);
      hyperonreco::KeepHitsInROI(TVector3(pfp.X_NoSC,pfp.Y_NoSC,pfp.Z_NoSC),hits,roi_size_ch,roi_size_tick);

      std::vector<std::vector<hyperonreco::HoughTransformPoint>> clusters;
      for(int i_pl=0;i_pl<3;i_pl++){

        std::pair<double,double> vertex_ch_tick = WireTick(TVector3(pfp.X_NoSC,pfp.Y_NoSC,pfp.Z_NoSC),i_pl); 
        double or_ch = vertex_ch_tick.first;
        double or_ti = vertex_ch_tick.second;

        hyperonreco::HoughTransformer transformer(hits.at(i_pl),i_pl,vertex_ch_tick.first,vertex_ch_tick.second,true);
        transformer.SetEvent(e.run,e.subrun,e.event,pfp_ctr);   
        transformer.MakeTransform2();
        transformer.FindPeaks2();
        clusters.push_back(transformer.MakeClusters());

      } // i_pl

      pfp_ctr++;

    } // pfp

  } // ievent

}

