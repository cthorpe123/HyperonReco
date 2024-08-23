#ifndef _FitOrganiser_cxx_
#define _FitOrganiser_cxx_

#include "FitOrganiser.h" 

namespace hyperonreco {

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

FitOrganiser::FitOrganiser(int run,int subrun,int event,RecoParticle pfp,const std::vector<std::vector<HoughTransformPoint>>& clusters,std::vector<std::vector<HitLite>> allhits) :
    Run(run),Subrun(subrun),Event(event),PFP(pfp),AllHits(allhits)
{

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

      for(size_t i=0;i<clusters_touse.size();i++) Fitter->AddData(*(clusters_touse.at(i)));

      //Fitter->SetEvent(Run,Subrun,Event,ctr);
      FittedV v;
      v.Vertex = TVector3(PFP.X_NoSC,PFP.Y_NoSC,PFP.Z_NoSC);
      Fitter->SetGuess(v);
      bool converged = Fitter->DoFitGridSearch3(v,NThrows); 
      //bool converged = Fitter->DoFit(v); 
      Fitter->Reset();

      std::cout << "v.GetOpeningAngle()=" << v.GetOpeningAngle() << std::endl;

      if(!converged || v.GetOpeningAngle() < OpeningAngleRange.first || v.GetOpeningAngle() > OpeningAngleRange.second) continue;

      FitResults[v.Chi2/pow(v.NDof,4)] = std::make_pair(clusters_touse,v); 

    }

  }
  
  int ctr=0;
  std::map<double,std::pair<std::vector<const HoughTransformPoint*>,FittedV>>::iterator it;
  for(it = FitResults.begin();it != FitResults.end();it++){
    Fitter->SetEvent(Run,Subrun,Event,ctr);
    for(size_t i=0;i<it->second.first.size();i++) Fitter->AddData(*(it->second.first.at(i)));
    Fitter->DrawFit2(it->second.second,AllHits);
    Fitter->Reset();
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

}

#endif
