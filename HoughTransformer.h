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
        void MakeTransform2();
        std::vector<HoughTransformPoint> FindPeaks() const;
        void DrawFits();
        void SubtractOffset(std::vector<HitLite>& hits);
        void RestoreOffset(std::vector<HitLite>& hits);
        std::vector<HoughTransformPoint> MakeClusters();
        
        void SetRBinSize(double size);
        void SetThetaBinSize(double size);
        void SetPeakSize(int size);
        void SetPointGrouping(int size);
        void SetMaxNeighbourDist(double dist);

        void SetEvent(int run,int subrun,int event);
        void SetTuneID(int tuneid);

      private:

        int Run,Subrun,Event;
        int TuneID = 0;
        std::string RSE;
        bool Draw;

        std::vector<HitLite> Hits;
        const int Plane;
        const double Origin_Channel;
        const double Origin_Tick;

        // Tuning parameters
        double RBinSize = 1;
        double ThetaBinSize = 0.1;
        int PeakSize = 2;
        int PointGrouping = 3;
        double MaxNeighbourDist2 = 5.0;

        std::vector<HitLite> Hits_test;

        std::vector<HitLite> FindNearestNeighbours(int point,const std::vector<HitLite>& hits,int num) const;
        double Dist(const HitLite& hit,double r, double theta) const;
        std::pair<double,double> GetLine(const std::vector<HitLite>& hits);

        TH2D* h_Transform = nullptr;
        std::vector<HoughTransformPoint> Transform;

        ROOT::Math::Functor Func;
        std::unique_ptr<ROOT::Math::Minimizer> Minimizer = nullptr;

  };

}

#endif
