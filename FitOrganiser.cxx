#ifndef _FitOrganiser_cxx_
#define _FitOrganiser_cxx_

#include "FitOrganiser.h" 

namespace hyperonreco {

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

FitOrganiser::FitOrganiser(){

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FitOrganiser::SetData(int run,int subrun,int event,RecoParticle pfp,const std::vector<std::vector<HoughTransformPoint>>& clusters,std::vector<std::vector<HitLite>> allhits){

  FitResults.clear();

  Run = run;
  Subrun = subrun;
  Event = event;
  PFP = pfp;
  AllHits = allhits;

  // Flatten the input vector out, with pointers to the original data
  ClustersFlat.clear();
  for(size_t i_pl=0;i_pl<3;i_pl++){
    for(const HoughTransformPoint& cluster : clusters.at(i_pl)){
      if(cluster.Hits.size() > MinHits){
        //std::cout << cluster.Plane << " " << cluster.Hits.size() << std::endl; 
        ClustersFlat.push_back(&cluster); 
      }
    }
  }


}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FitOrganiser::MakeFitList(){

  //Fitter->Reset();

  for(size_t nclusters : NClusters){

    if(nclusters > ClustersFlat.size()) break;

    //std::cout << "combo=" <<  combo << " clusters_pv.size()=" <<  clusters_pv.size() << std::endl;

    for(unsigned int test=0;test<nChoosek(ClustersFlat.size(),nclusters)*2;test++){

      // Make a copy of the vector we can erase elements from as they're used
      std::vector<const HoughTransformPoint*> clusters_cp = ClustersFlat;
      std::vector<const HoughTransformPoint*> clusters_touse; 
      while(clusters_touse.size() < nclusters){
        size_t pos = RNG->Uniform(-0.5,clusters_cp.size()-0.5);
        clusters_touse.push_back(clusters_cp.at(pos));
        clusters_cp.erase(clusters_cp.begin()+pos);   
      }

      // Check clusters don't all belong to the same plane
      int pl=clusters_touse.at(0)->Plane;
      bool two_planes = false;
      for(size_t i=1;i<clusters_touse.size();i++)
        if(clusters_touse.at(i)->Plane != pl) two_planes = true;

      if(!two_planes || AlreadyTested(clusters_touse)) continue;

      VFitter Fitter(true); 

      for(size_t i=0;i<clusters_touse.size();i++) Fitter.AddData(*(clusters_touse.at(i)));

      //Fitter->SetEvent(Run,Subrun,Event,ctr);
      FittedV v;
      v.Vertex = TVector3(PFP.X_NoSC,PFP.Y_NoSC,PFP.Z_NoSC);
      Fitter.SetGuess(v);

      //bool converged = Fitter->DoFitGridSearch3(v,NThrows); 
      bool converged = Fitter.DoFit(v); 

      //Fitter->Reset();

      if(!converged || v.GetOpeningAngle() < OpeningAngleRange.first || v.GetOpeningAngle() > OpeningAngleRange.second) continue;

      FitResults[v.Chi2/pow(v.NDof,4)] = std::make_pair(clusters_touse,v); 

      if(t_Output != nullptr){
        t_run = Run;
        t_subrun = Subrun;
        t_event = Event;
        t_PFP = PFP.Index;
        
        for(int i_pl=0;i_pl<3;i_pl++){
          t_HitChannels.at(i_pl).clear();
          t_HitTicks.at(i_pl).clear();
          t_HitWidths.at(i_pl).clear();
        }      

        for(size_t i_c=0;i_c<clusters_touse.size();i_c++){
          HoughTransformPoint cluster = *(clusters_touse.at(i_c));
          for(HitLite hit : cluster.Hits){
            t_HitChannels.at(cluster.Plane).push_back(hit.Channel);
            t_HitTicks.at(cluster.Plane).push_back(hit.Tick);
            t_HitWidths.at(cluster.Plane).push_back(hit.Width);
          }
        }

        t_Score = v.Chi2;
        t_NHits = {t_HitChannels.at(0).size(),t_HitChannels.at(1).size(),t_HitChannels.at(2).size()};

        t_Output->Fill(); 

      }

    }

  }
  
  int ctr=0;
  std::map<double,std::pair<std::vector<const HoughTransformPoint*>,FittedV>>::iterator it;
  for(it = FitResults.begin();it != FitResults.end();it++){
    VFitter Fitter(true); 
    Fitter.SetEvent(Run,Subrun,Event,ctr);
    for(size_t i=0;i<it->second.first.size();i++) Fitter.AddData(*(it->second.first.at(i)));
    Fitter.DrawFit2(it->second.second,AllHits);
    Fitter.Reset();
    ctr++;
    if(ctr > 10) break;
  }  
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool FitOrganiser::AlreadyTested(const std::vector<const HoughTransformPoint*>& clusters){

  std::map<double,std::pair<std::vector<const HoughTransformPoint*>,FittedV>>::iterator it;

  for(it = FitResults.begin();it != FitResults.end();it++){

    const std::vector<const HoughTransformPoint*>& previous_fit_clusters = it->second.first;

    std::vector<bool> found(clusters.size(),false);

    if(previous_fit_clusters.size() != clusters.size()) continue;

    for(size_t i_p=0;i_p<clusters.size();i_p++)
      if(std::find(previous_fit_clusters.begin(),previous_fit_clusters.end(),clusters.at(i_p)) != previous_fit_clusters.end())
        found.at(i_p) = true;

    // if all elements of found are true, we know we've already tried this combination
    if(std::find(found.begin(),found.end(),false) == found.end())
      return true;

  }

  return false;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FitOrganiser::SetTreePtr(TTree* t){

  t_Output = t;
  InitialiseTree(); 

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FitOrganiser::InitialiseTree(){

  std::cout << "Output tree being set up" << std::endl;
  t_Output->Branch("run",&t_run);
  t_Output->Branch("subrun",&t_subrun);
  t_Output->Branch("event",&t_event);
  t_Output->Branch("PFP",&t_PFP);
  //t_Output->Branch("Fit",&t_Fit);
  t_Output->Branch("HitChannels",&t_HitChannels);
  t_Output->Branch("HitTicks",&t_HitTicks);
  t_Output->Branch("HitWidths",&t_HitWidths);
  t_Output->Branch("Score",&t_Score);
  t_Output->Branch("NHits",&t_NHits);

  t_HitChannels.resize(3);
  t_HitTicks.resize(3);
  t_HitWidths.resize(3);
  t_NHits.resize(3);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

}

#endif
