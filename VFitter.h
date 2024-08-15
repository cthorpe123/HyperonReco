#ifndef _VFitter_h_
#define _VFitter_h_

// This object is compilable outside of LArSoft, comment out this line to do so
//#define _InsideLArSoft_

// C++ STL includes
#include <vector>
#include <stdexcept>
#include <iomanip>

#ifdef _InsideLArSoft_
// larsoft includes
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#endif

// root includes
#include "TVector3.h"
#include "TVector2.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "TRandom2.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"

// Local includes
//#include "RecoParticle.h"
#ifdef _InsideLArSoft_ 
#include "ubana/HyperonProduction/Alg/Position_To_Wire.h"
#else 
#include "Position_To_Wire.h" 
#endif
#include "FittedV.h"
#include "HoughTransformer.h"
#include "Objects.h"

namespace hyperonreco {

  class VFitter {

    public:

      VFitter(bool draw=false);

      void Reset();
      void AddData(const std::vector<std::vector<HitLite>>& hits);
      void AddData(const HoughTransformPoint& p);
      void SetGuess(const FittedV& fittedv);
      void SetEvent(int run,int subrun,int event,int combination=0);

      bool DoFitGridSearch3(FittedV& fittedv,int points);
      void SetActivePlanes(std::vector<int> i_pl);

      void SetROI(double roi_ch,double roi_tick,TVector3 center);
      void SetOutlierCut(double outliercut);
      void SetFitTune(double fittune);
      void RemoveOffset(std::vector<HitLite>& hits,TVector3 origin,int plane);
      void RestoreOffset(std::vector<HitLite>& hits,TVector3 origin,int plane);
      void DrawFit(const FittedV& v,const std::vector<std::vector<HitLite>>& allhits={}) const;
      void DrawFit2(const FittedV& v,const std::vector<std::vector<HitLite>>& allhits={}) const;



    private:

      int Run,Subrun,Event,Combination;
      std::string RSE;
      std::string RSEC;
      bool Draw = false;

      double ROIChannel = 100;
      double ROITick = 1500;
      TVector3 ROICenter;
      std::vector<int> ActivePlanes = {0,1,2};
      double OutlierCut = 50;
      double FitTune = 1;

      std::vector<std::vector<HitLite>> Hits;
 
      FittedV InitialGuess;
      bool InitialGuessSet = false;

      std::pair<double,int> FitScore2(const std::vector<std::vector<HitLite>>& hits,FittedV& fittedv,bool verbose=false);
      double HitLineSeparation2(const HitLite& hit,const LineWireTick2& line);

  };


}

#endif
