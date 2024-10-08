#ifndef _FitOrganiser_h_
#define _FitOrganiser_h_

// stdlib includes
#include <vector>

// root includes
#include "TTree.h"

// local includes
#include "VFitter.h"
#include "Objects.h"
#include "RecoParticle.h"

namespace hyperonreco {

  inline unsigned int nChoosek(unsigned int n,unsigned int k){
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    int result = n;
    for(unsigned int i = 2; i <= k; ++i ) {
      result *= (n-i+1);
      result /= i;
    }
    return result;
  }


  class FitOrganiser {

    public:

      FitOrganiser();
      void SetData(int run,int subrun,int event,RecoParticle pfp,const std::vector<std::vector<HoughTransformPoint>>& clusters,std::vector<std::vector<HitLite>> allhits = {});
      std::map<double,std::pair<std::vector<const HoughTransformPoint*>,FittedV>> FitResults;
      void MakeFitList();
      void SetTreePtr(TTree* t);

    private:

      int Run,Subrun,Event;
      RecoParticle PFP;
      int NThrows = 1000;
      size_t MinHits = 5;
      std::pair<double,double> OpeningAngleRange = {0.15,1.8};
      std::vector<int> NClusters = {3,4,5,6,7};

      std::vector<const HoughTransformPoint*> ClustersFlat; 
      std::vector<std::vector<HitLite>> AllHits;

      TRandom2* RNG = new TRandom2();

      bool AlreadyTested(const std::vector<const HoughTransformPoint*>& clusters);

      //VFitter* Fitter = new VFitter(true);

      // Recording the outputs      
      TTree* t_Output = nullptr;
      int t_run,t_subrun,t_event,t_PFP;
      int t_Fit;
      std::vector<std::vector<double>> t_HitChannels,t_HitWidths,t_HitTicks;
      double t_Score;
      std::vector<unsigned int> t_NHits;
      //const std::vector<std::string> Branches = {"run","subrun","event","PFP","Fit","HitChannels","HitWidths","HitTicks","Score","NHits"};

      void InitialiseTree();

   
  };


}

#endif
