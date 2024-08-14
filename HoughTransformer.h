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

// local includes
#include "Position_To_Wire.h" 
#include "Objects.h"

namespace hyperonreco {

  class HoughTransformer {

      public:

        HoughTransformer(std::vector<HitLite> hits,int plane,double origin_channel,double origin_tick,bool draw=false);
        ~HoughTransformer();
        void MakeTransform2();
        std::vector<HoughTransformPoint> FindPeaks() const;
        std::vector<HoughTransformPoint> FindPeaks2() const;
        void DrawFits();
        void SubtractOffset(std::vector<HitLite>& hits) const;
        void RestoreOffset(std::vector<HitLite>& hits) const;
        std::vector<HoughTransformPoint> MakeClusters() const;
        
        void SetRBinSize(double size);
        void SetThetaBinSize(double size);
        void SetPeakSize(int size);
        void SetPointGrouping(int size);
        void SetMaxNeighbourDist(double dist);
        void SetConvFloor(int fl);

        void SetEvent(int run,int subrun,int event,int pfp=0);
        void SetTuneID(int tuneid);

        std::pair<double,double> GetPerformanceMetrics() const;

      private:

        int Run,Subrun,Event,Pfp;
        int TuneID = 0;
        std::string RSE;
        std::string RSEP;
        bool Draw;

        std::vector<HitLite> Hits;
        const int Plane;
        const double Origin_Channel;
        const double Origin_Tick;

        // Tuning parameters
        double RBinSize = 3.7;
        double ThetaBinSize = 0.18;
        int PeakSize = 1;
        int PointGrouping = 4;
        int ConvFloor = 2;
        double MaxNeighbourDist2 = 11*11;

        std::vector<HitLite> Hits_test;

        std::vector<HitLite> FindNearestNeighbours(int point,const std::vector<HitLite>& hits,size_t num) const;
        double Dist(const HitLite& hit,double r, double theta) const;
        std::pair<double,double> GetLine(const std::vector<HitLite>& hits);

        TH2D* h_Transform = nullptr;
        std::vector<HoughTransformPoint> Transform;

        ROOT::Math::Functor Func;
        std::unique_ptr<ROOT::Math::Minimizer> Minimizer = nullptr;

  };

}

#endif
