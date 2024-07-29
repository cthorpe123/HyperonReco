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

namespace hyperonreco {

  struct HoughTransformPoint {

    HoughTransformPoint(int plane,double r,double theta) :
      Plane(plane),R(r),Theta(theta)
    {}

    int Plane;
    double R;
    double Theta;

    int Height;

    std::vector<double> Channels;
    std::vector<double> Ticks;
    std::vector<double> Widths;
    std::vector<int> PointNo; // int used to uniquely identify every point - used to avoid duplicate points later

  };

  class HoughTransformer {

      public:

        HoughTransformer(std::vector<double> channels,std::vector<double> ticks,std::vector<double> widths,int plane,double origin_channel,double origin_tick);
        void MakeTransform();
        void MakeTransform2();
        std::vector<HoughTransformPoint> FindPeaks() const;
        void DrawFits();
        void SubtractOffset(std::vector<double>& channels,std::vector<double>& ticks,std::vector<double>& widths);
        void RestoreOffset(std::vector<double>& channels,std::vector<double>& ticks,std::vector<double>& widths);
        std::vector<HoughTransformPoint> MakeClusters();
        
        void SetRBinSize(double size);
        void SetThetaBinSize(double size);
        void SetPeakSize(int size);
        void SetPointGrouping(int size);

        void SetEvent(int run,int subrun,int event);
        void SetTuneID(int tuneid);

      private:

        int Run,Subrun,Event;
        int TuneID = 0;

        std::vector<double> Channels_v;
        std::vector<double> Ticks_v;
        std::vector<double> Widths_v;
        const int Plane;
        const double Origin_Channel;
        const double Origin_Tick;

        double RBinSize = 1;
        double ThetaBinSize = 0.1;
        int PeakSize = 2;
        int PointGrouping = 3;

        std::vector<double> Channels_test;
        std::vector<double> Ticks_test;

        std::vector<int> FindNearestNeighbours(int point,const std::vector<double>& channels,const std::vector<double>& ticks,int num) const;

        double Dist(double channel,double tick,double r, double theta) const;
        std::pair<double,double> GetLine(std::vector<double> channels,std::vector<double> ticks);

        TH2D* h_Transform = nullptr;
        std::vector<HoughTransformPoint> Transform;

        ROOT::Math::Functor Func;
        std::unique_ptr<ROOT::Math::Minimizer> Minimizer = nullptr;

  };

}

#endif
