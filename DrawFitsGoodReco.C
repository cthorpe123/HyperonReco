R__LOAD_LIBRARY(libVFitter.so);
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

  // Select the event we want to analyse
  int ievent=0;
  while(ievent<E.GetNEvents() && ievent < 10){
    Event e = E.GetEvent(ievent);

    std::cout << ievent << "/" << E.GetNEvents() << std::endl;
    ievent++;

    //if(e.run != 6999 || e.subrun != 45 || e.event != 2285) continue;

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

    std::vector<std::vector<hyperonreco::HoughTransformPoint>> clusters(3);

    for(size_t i_pl=0;i_pl<3;i_pl++){

      if(hits.at(i_pl).size() < 7) continue;

      std::pair<double,double> vertex_ch_tick = WireTick(TVector3(pfps.at(0).X_NoSC,pfps.at(0).Y_NoSC,pfps.at(0).Z_NoSC),i_pl); 
      double or_ch = vertex_ch_tick.first;
      double or_ti = vertex_ch_tick.second;

      hyperonreco::HoughTransformer transformer(hits.at(i_pl),i_pl,vertex_ch_tick.first,vertex_ch_tick.second,true);
      transformer.SetEvent(e.run,e.subrun,e.event);   
      transformer.MakeTransform2();
      std::vector<hyperonreco::HoughTransformPoint> cluster = transformer.MakeClusters();
      clusters.at(i_pl) = cluster;

    } // i_pl

    //FindBestFit(e.run,e.subrun,e.event,clusters,pfps.at(0),hits);
    FitOrganiser organiser(e.run,e.subrun,e.event,pfps.at(0),clusters,hits);
    organiser.MakeFitList();


    /*
       int ctr=0;
       VFitter fitter(true);

       double bestfit = 1e10;
       int bestfit_ctr=-1;

       for(HoughTransformPoint cluster_Plane0 : clusters.at(0)){
       if(cluster_Plane0.Hits.size() < 5) continue;        
       for(HoughTransformPoint cluster_Plane1 : clusters.at(1)){
       if(cluster_Plane1.Hits.size() < 5) continue;        
       for(HoughTransformPoint cluster_Plane2 : clusters.at(2)){
       if(cluster_Plane2.Hits.size() < 5) continue;        

       fitter.SetEvent(e.run,e.subrun,e.event,ctr);   
       hyperonreco::FittedV v = hyperonreco::MakeFittedVGuessTrack(pfps.at(0));
       fitter.AddData(cluster_Plane0);
       fitter.AddData(cluster_Plane1);
       fitter.AddData(cluster_Plane2);
       fitter.SetGuess(v);
       fitter.DoFitGridSearch3(v,5000); 
       fitter.DrawFit(v,hits);
       fitter.Reset();

       if(v.Chi2/v.NDof/v.NDof < bestfit){
       bestfit = v.Chi2/v.NDof/v.NDof;
       bestfit_ctr = ctr;
       }

       ctr++;
       }
       }
       }

       std::cout << "Best fit: " << bestfit_ctr << std::endl;
       */

  } // ievent

}

