R__LOAD_LIBRARY(libVFitter.so);
R__LOAD_LIBRARY($HYP_TOP/lib/libHyperon.so);
R__LOAD_LIBRARY($HYP_TOP/lib/libParticleDict.so);

#include "EventAssembler.h"
//#include "VFitter.h"
//#include "SpacePointVisualisation.h" 
#include "HoughTransformer.h"
#include <ctime>    
#include <chrono>

// Apply Hough transform to make clusters on already reconstructed events, draw fitted Vs to resulting clusters

void DrawFitsGoodReco(){

  const std::string filename = "GoodReco.root";

  // number of pfps the transformer has been tested - use to get the average performance metrics
  int n_tests = 0;

  // Merge hit vectors together an remove hits outside roi 
  double roi_size_ch = 100; 
  double roi_size_tick = roi_size_ch*TickPerWire; 

  EventAssembler E(false);
  E.SetFile(filename,"Background");

  // Select the event we want to analyse
  int ievent=0;
  while(ievent<E.GetNEvents()){
    Event e = E.GetEvent(ievent);

    std::cout << ievent << "/" << E.GetNEvents() << std::endl;
    ievent++;

    if(e.run != 7001 || e.subrun != 964 || e.event != 48218) continue;


    std::vector<RecoParticle> pfps;
    bool has_proton=false,has_pion=false;

    for(RecoParticle pfp : e.TracklikePrimaryDaughters){
      if(pfp.TrackTruePDG == 2212 && pfp.TrackTrueOrigin == 2 && !has_proton){
        pfps.push_back(pfp);
        has_proton = true;
      }
      if(pfp.TrackTruePDG == -211 && pfp.TrackTrueOrigin == 2 && !has_pion){
        pfps.push_back(pfp);
        has_pion = true;
      }
      if(has_proton && has_pion) break;
    }

    if(!has_proton || !has_pion) continue;

    std::vector<std::vector<hyperonreco::HitLite>> hits(3);

    for(size_t i_pfp=0;i_pfp<pfps.size();i_pfp++){ 
      RecoParticle pfp = pfps.at(i_pfp);
      hyperonreco::AddHits(hits,pfp.HitChannels,pfp.HitTicks,pfp.HitWidths,pfp.HitTrackIDs,pfp.HitPDGs);
      hyperonreco::KeepHitsInROI(TVector3(pfp.X_NoSC,pfp.Y_NoSC,pfp.Z_NoSC),hits,roi_size_ch,roi_size_tick);
    }

    n_tests++;

    for(int i_pl=0;i_pl<3;i_pl++){

      if(hits.at(i_pl).size() < 7) continue;

      std::pair<double,double> vertex_ch_tick = WireTick(TVector3(pfps.at(0).X_NoSC,pfps.at(0).Y_NoSC,pfps.at(0).Z_NoSC),i_pl); 
      double or_ch = vertex_ch_tick.first;
      double or_ti = vertex_ch_tick.second;

      hyperonreco::HoughTransformer transformer(hits.at(i_pl),i_pl,vertex_ch_tick.first,vertex_ch_tick.second,true);
      transformer.SetEvent(e.run,e.subrun,e.event);   
      transformer.MakeTransform2();
      //transformer.FindPeaks2();
      transformer.MakeClusters();

    } // i_pl

  } // ievent

}

