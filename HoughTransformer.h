#ifndef _HoughTransformer_h_
#define _HoughTransformer_h_

// This object is compilable outside of LArSoft, comment out this line to do so
//#define _InsideLArSoft_

// C++ STL includes
#include <vector>
#include <stdexcept>
#include <iomanip>
#include <algorithm>

// root includes
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "TLegend.h"

// local includes
#include "Position_To_Wire.h" 
#include "Objects.h"

namespace hyperonreco {

  const double PI = 3.1415;

  // 2D Gaussian of indep variables
  inline double Gaus2D(double mu_x,double mu_y,double sd_x,double sd_y,double x,double y){
    return (1.0/2/PI)/sqrt(sd_x*sd_y)*exp(-0.5*((mu_x-x)*(mu_x-x)/sd_x/sd_x + (mu_y-y)*(mu_y-y)/sd_y/sd_y));
  }

  inline double Gaus2D_GradX(double mu_x,double mu_y,double sd_x,double sd_y,double x,double y){
    return (1.0/2/PI)/sqrt(sd_x*sd_y)*(mu_x-x)/sd_x/sd_x*exp(-0.5*((mu_x-x)*(mu_x-x)/sd_x/sd_x + (mu_y-y)*(mu_y-y)/sd_y/sd_y));
  }
   
  inline double Gaus2D_GradY(double mu_x,double mu_y,double sd_x,double sd_y,double x,double y){
    return (1.0/2/PI)/sqrt(sd_x*sd_y)*(mu_y-y)/sd_y/sd_y*exp(-0.5*((mu_x-x)*(mu_x-x)/sd_x/sd_x + (mu_y-y)*(mu_y-y)/sd_y/sd_y));
  }

/*
  inline double Gaus2D_GradX2(double mu_x,double mu_y,double sd_x,double sd_y,double x,double y){
    return (1.0/2/PI)/sqrt(sd_x*sd_y)*(mu_x*mu_x - 2*mu_x*x - sd_x*sd_x + x*x)/sd_x/sd_x/sd_x/sd_x*exp(-0.5*((mu_x-x)*(mu_x-x)/sd_x/sd_x + (mu_y-y)*(mu_y-y)/sd_y/sd_y));
  }

  inline double Gaus2D_GradY2(double mu_x,double mu_y,double sd_x,double sd_y,double x,double y){
    return (1.0/2/PI)/sqrt(sd_x*sd_y)*(mu_y*mu_y - 2*mu_y*y - sd_y*sd_y + y*y)/sd_y/sd_y/sd_y/sd_y*exp(-0.5*((mu_x-x)*(mu_x-x)/sd_x/sd_x + (mu_y-y)*(mu_y-y)/sd_y/sd_y));
  }
*/

  inline bool SortHits(HitLite a,HitLite b){
    return a.Tick < b.Tick;
  }

  class HoughTransformer {

      public:

        HoughTransformer(std::vector<HitLite> hits,int plane,double origin_channel,double origin_tick,bool draw=false);
        ~HoughTransformer();
        void MakeTransform2();
        std::vector<HoughTransformPoint> FindPeaks2() const;
        std::vector<HoughTransformPoint> FindPeaks3() const;
        void DrawFits();
        void SubtractOffset(std::vector<HitLite>& hits) const;
        void RestoreOffset(std::vector<HitLite>& hits) const;
        std::vector<HoughTransformPoint> MakeClusters() const;
        
        void SetRBinSize(double size);
        void SetThetaBinSize(double size);
        void SetPeakSize(int size);
        void SetPointGrouping(int size);
        void SetMaxNeighbourDist(double dist);
        void SetChi2Cut(double chi2);
        void SetConvFloor(int fl);

        void SetEvent(int run,int subrun,int event,int pfp=0);
        void SetTuneID(int tuneid);
        void SetVerbosity(int verb);

        std::pair<double,double> GetPerformanceMetrics() const;

      private:

        int Run,Subrun,Event,Pfp;
        int TuneID = 0;
        std::string RSE;
        std::string RSEP;
        bool Draw;
        int Verbosity = 0;

        std::vector<HitLite> Hits;
        const int Plane;
        const double Origin_Channel;
        const double Origin_Tick;

        // Tuning parameters
        double RBinSize = 2.0;
        double ThetaBinSize = 0.2;
        int PeakSize = 1;
        int PointGrouping = 4;
        int ConvFloor = 2;
        double MaxNeighbourDist2 = 13*13;
        double Chi2Cut = 2.8;

        std::vector<HitLite> Hits_test;

        std::vector<HitLite> FindNearestNeighbours(int point,const std::vector<HitLite>& hits,size_t num) const;
        double Dist(const HitLite& hit,double r, double theta) const;
        std::tuple<double,double,double> GetLine(const std::vector<HitLite>& hits);

        //TH2D* h_Transform = nullptr;
        std::vector<HoughTransformPoint> Transform;
 
        const int Multiplier = 3;
        TH2D* MakeConvHistogram(std::vector<HoughTransformPoint> transform,bool guasskernel) const;
        TH2D* MakeConvHistogramDer(std::vector<HoughTransformPoint> transform,bool xy) const;


        ROOT::Math::Functor Func;
        std::unique_ptr<ROOT::Math::Minimizer> Minimizer = nullptr;
        
        void RemoveVerticalLines();
        std::vector<HoughTransformPoint> VerticalLines; 

  };

}

#endif
