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

using namespace hyperonreco;

unsigned nChoosek( unsigned n, unsigned k )
{
  if (k > n) return 0;
  if (k * 2 > n) k = n-k;
  if (k == 0) return 1;

  int result = n;
  for( int i = 2; i <= k; ++i ) {
    result *= (n-i+1);
    result /= i;
  }
  return result;
}

void FindBestFit(int run,int subrun,int event,const std::vector<std::vector<HoughTransformPoint>>& clusters,const RecoParticle& pfp,const std::vector<std::vector<HitLite>>& hits){

  VFitter fitter(true);
  const std::vector<size_t> combos  = {2,3,4,5,6,7};
  const int min_hits = 5;

  double bestfit = 1e10;
  int bestfit_ctr=-1;

  // Flatten the input vector out, with pointers to the original data
  std::vector<const HoughTransformPoint*> clusters_pv;
  for(size_t i_pl=0;i_pl<3;i_pl++){
    for(const HoughTransformPoint& cluster : clusters.at(i_pl)){
      if(cluster.Hits.size() > min_hits){
        //std::cout << cluster.Plane << " " << cluster.Hits.size() << std::endl; 
        clusters_pv.push_back(&cluster); 
      }
    }
  }

  TRandom2* r = new TRandom2();

  int ctr=0;

  std::vector<std::vector<const HoughTransformPoint*>> already_tried;

  for(size_t combo : combos){

    if(combo > clusters_pv.size()) break;

    //std::cout << "combo=" <<  combo << " clusters_pv.size()=" <<  clusters_pv.size() << std::endl;

    // Make a copy of the vector we can erase elements from as they're used

    for(int c=0;c<nChoosek(clusters_pv.size(),combo)*2;c++){

      std::vector<const HoughTransformPoint*> clusters_pv_cp = clusters_pv;
      std::vector<const HoughTransformPoint*> clusters_pv_touse; 
      while(clusters_pv_touse.size() < combo){
        size_t pos = r->Uniform(-0.5,clusters_pv_cp.size()-0.5);
        clusters_pv_touse.push_back(clusters_pv_cp.at(pos));
        clusters_pv_cp.erase(clusters_pv_cp.begin()+pos);   
      }

      // Check clusters don't all belong to the same plane
      int pl=clusters_pv_touse.at(0)->Plane;
      bool two_planes = false;
      for(size_t i=1;i<clusters_pv_touse.size();i++)
        if(clusters_pv_touse.at(i)->Plane != pl) two_planes = true;

      if(!two_planes) continue;

      // Check we haven't already tried this combination
      bool this_combo_already_tried = false;
      for(std::vector<const HoughTransformPoint*> points : already_tried){  
        std::vector<bool> found(clusters_pv_touse.size(),false);
        for(size_t i_p=0;i_p<clusters_pv_touse.size();i_p++)
          if(std::find(points.begin(),points.end(),clusters_pv_touse.at(i_p)) != points.end())
            found.at(i_p) = true;

        // if all elements of found are true, we know we've already tried this combination
        if(std::find(found.begin(),found.end(),false) == found.end())
          this_combo_already_tried = true;

      }

      if(this_combo_already_tried) continue;

      //std::cout << std::endl << "Selected Clusters:" << std::endl;
      //std::vector<std::vector<HoughTransformPoint>> clusters_touse(3);
      for(size_t i=0;i<clusters_pv_touse.size();i++){
        //std::cout << clusters_pv_touse.at(i) << "  " <<  clusters_pv_touse.at(i)->Plane << " " << clusters_pv_touse.at(i)->Hits.size() << std::endl;
        //clusters_touse.at(clusters_pv_touse.at(i)->Plane).push_back(*(clusters_pv_touse.at(i)));
        fitter.AddData(*(clusters_pv_touse.at(i)));
      }


      fitter.SetEvent(run,subrun,event,ctr);
      hyperonreco::FittedV v = hyperonreco::MakeFittedVGuessTrack(pfp);
      fitter.SetGuess(v);
      fitter.DoFitGridSearch3(v,1000); 
      fitter.DrawFit(v,hits);
      fitter.Reset();

      std::cout << "Score = " << v.Chi2/v.NDof/v.NDof << std::endl;

      if(v.Chi2/v.NDof/v.NDof < bestfit && v.GetAsymmetry() < 0.9 && v.GetOpeningAngle() < 1.75){
        bestfit = v.Chi2/v.NDof/v.NDof;
        bestfit_ctr = ctr;
      }

      already_tried.push_back(clusters_pv_touse);     
      ctr++;

    }

  }

  for(size_t i=0;i<already_tried.at(bestfit_ctr).size();i++)
    fitter.AddData(*(already_tried.at(bestfit_ctr).at(i)));

  fitter.SetEvent(run,subrun,event,-1);
  hyperonreco::FittedV v = hyperonreco::MakeFittedVGuessTrack(pfp);
  fitter.SetGuess(v);
  fitter.DoFitGridSearch3(v,1000); 
  fitter.DrawFit(v);

}


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

    FindBestFit(e.run,e.subrun,e.event,clusters,pfps.at(0),hits);

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

