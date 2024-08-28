R__LOAD_LIBRARY(lib/libVFitter.so);
R__LOAD_LIBRARY($HYP_TOP/lib/libHyperon.so);
R__LOAD_LIBRARY($HYP_TOP/lib/libParticleDict.so);

#include "EventAssembler.h"
//#include "VFitter.h"
#include "SpacePointVisualisation.h" 
#include "HoughTransformer.h"
#include "VFitter.h"
#include <ctime>    
#include <chrono>
#include "FitOrganiser.h"

using namespace hyperonreco;

// Apply Hough transform to make clusters on already reconstructed events, draw fitted Vs to resulting clusters

void DrawFitsGoodReco(){

  const std::string filename = "GoodReco.root";

  // Merge hit vectors together an remove hits outside roi 
  double roi_size_ch = 100; 
  double roi_size_tick = roi_size_ch*TickPerWire; 

  EventAssembler E(false);
  E.SetFile(filename,"Background");

  FitOrganiser organiser;
  TFile* f_Output = TFile::Open("rootfiles/FitResults.root","RECREATE");
  TTree* t_Output = new TTree("Output","Output");
  organiser.SetTreePtr(t_Output);

  // Select the event we want to analyse
  int ievent=0;
  while(ievent<E.GetNEvents() && ievent < 10){
    Event e = E.GetEvent(ievent);

    //std::cout << ievent << "/" << E.GetNEvents() << std::endl;
    ievent++;

    //if(e.run != 7003 || e.subrun != 682 || e.event != 34102) continue;

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

    std::vector<std::vector<HitLite>> hits(3);

    for(size_t i_pfp=0;i_pfp<pfps.size();i_pfp++){ 
      RecoParticle pfp = pfps.at(i_pfp);
      AddHits(hits,pfp.HitChannels,pfp.HitTicks,pfp.HitWidths,pfp.HitTrackIDs,pfp.HitPDGs);
      KeepHitsInROI(TVector3(pfp.X_NoSC,pfp.Y_NoSC,pfp.Z_NoSC),hits,roi_size_ch,roi_size_tick);
    }

    std::vector<std::vector<HoughTransformPoint>> clusters(3);

    for(size_t i_pl=0;i_pl<3;i_pl++){

      if(hits.at(i_pl).size() < 7) continue;

      std::pair<double,double> vertex_ch_tick = WireTick(TVector3(pfps.at(0).X_NoSC,pfps.at(0).Y_NoSC,pfps.at(0).Z_NoSC),i_pl); 
      double or_ch = vertex_ch_tick.first;
      double or_ti = vertex_ch_tick.second;

      HoughTransformer transformer(hits.at(i_pl),i_pl,vertex_ch_tick.first,vertex_ch_tick.second,true);
      transformer.SetVerbosity(1);
      transformer.SetEvent(e.run,e.subrun,e.event);   
      transformer.MakeTransform2();
      std::vector<HoughTransformPoint> cluster = transformer.MakeClusters();
      transformer.FindPeaks3();
      //clusters.at(i_pl) = cluster;

    } // i_pl

    //organiser.SetData(e.run,e.subrun,e.event,pfps.at(0),clusters,hits);
    //organiser.MakeFitList();

  } // ievent

  t_Output->Write();
  f_Output->Close();


}

