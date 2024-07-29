#ifndef _HoughTransformer_cxx_
#define _HoughTransformer_cxx_

#include "HoughTransformer.h"

namespace hyperonreco { 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

HoughTransformer::HoughTransformer(std::vector<double> channels,std::vector<double> ticks,std::vector<double> widths,int plane,double origin_channel,double origin_tick) : 
  Channels_v(channels),Ticks_v(ticks),Widths_v(widths),Plane(plane),Origin_Channel(origin_channel),Origin_Tick(origin_tick){

  Func = ROOT::Math::Functor( [&] (const double *coeff ){
      double val = 0.0;
      for(size_t i=0;i<Channels_test.size();i++)
          val += pow(Dist(Channels_test.at(i),Ticks_test.at(i),coeff[0],coeff[1]),2);
      return val;
  } , 2);
 
  Minimizer = std::unique_ptr<ROOT::Math::Minimizer>( ROOT::Math::Factory::CreateMinimizer( "Minuit2", "Migrad" ) );

  Minimizer->SetMaxFunctionCalls(1000);
  Minimizer->SetTolerance(0.01);
  Minimizer->SetFunction(Func);

  SubtractOffset(Channels_v,Ticks_v,Widths_v);

  system("mkdir Plots/");
 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HoughTransformer::SetEvent(int run,int subrun,int event){

  Run = run;
  Subrun = subrun;
  Event = event;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HoughTransformer::SetTuneID(int tuneid){

  TuneID = tuneid;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HoughTransformer::SetRBinSize(double size){

  RBinSize = size;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HoughTransformer::SetThetaBinSize(double size){

  ThetaBinSize = size;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HoughTransformer::SetPeakSize(int size){

  PeakSize = size;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HoughTransformer::SetPointGrouping(int size){

  PointGrouping = size;

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HoughTransformer::SubtractOffset(std::vector<double>& channels,std::vector<double>& ticks,std::vector<double>& widths){

  for(size_t i=0;i<channels.size();i++){
    channels.at(i) = channels.at(i) - Origin_Channel;
    ticks.at(i) = (ticks.at(i) - Origin_Tick)/TickPerWire;
    widths.at(i) /= TickPerWire;
  }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HoughTransformer::RestoreOffset(std::vector<double>& channels,std::vector<double>& ticks,std::vector<double>& widths){

  for(size_t i=0;i<channels.size();i++){
    channels.at(i) = channels.at(i) + Origin_Channel;
    ticks.at(i) = ticks.at(i)*TickPerWire + Origin_Tick;
    widths.at(i) *= TickPerWire;
  }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

double HoughTransformer::Dist(double channel,double tick,double r, double theta) const {

  return fabs(channel*cos(theta) + tick*sin(theta) - r);

}
  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::pair<double,double> HoughTransformer::GetLine(std::vector<double> channels,std::vector<double> ticks){

  Channels_test = channels;
  Ticks_test = ticks;

  Minimizer->SetVariable(0,"r",0,0.01);
  Minimizer->SetVariable(1,"theta",0,0.01);
  Minimizer->SetVariableLimits(1,-3.1415,3.1415);
  Minimizer->Minimize();

  return std::make_pair(Minimizer->X()[0],Minimizer->X()[1]);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HoughTransformer::MakeTransform(){

  std::cout << "Attempting Hough transform of " << Channels_v.size() << " points" << std::endl;

  std::vector<std::pair<double,double>> transform;
  double theta_min=3.1415,theta_max=-3.1415;
  double r_min=1e10,r_max=-1e10;

  for(size_t i=0;i<Channels_v.size();i++){
    for(size_t j=i+1;j<Channels_v.size();j++){
      for(size_t k=j+1;k<Channels_v.size();k++){

        if(i == j || i == k || j == k) continue;                           

        std::pair<double,double> line = GetLine({Channels_v.at(i),Channels_v.at(j),Channels_v.at(k)},{Ticks_v.at(i),Ticks_v.at(j),Ticks_v.at(k)});
        r_min = std::min(r_min,line.first);
        r_max = std::max(r_max,line.first);
        theta_min = std::min(theta_min,line.second);
        theta_max = std::max(theta_max,line.second);
        transform.push_back(line);

        Transform.push_back(HoughTransformPoint(Plane,line.first,line.second)); 
        Transform.back().Height++;
        Transform.back().Channels = {Channels_v.at(i),Channels_v.at(j),Channels_v.at(k)};  
        Transform.back().Ticks = {Ticks_v.at(i),Ticks_v.at(j),Ticks_v.at(k)};  
        Transform.back().Widths = {Widths_v.at(i),Widths_v.at(j),Widths_v.at(k)};  
        Transform.back().PointNo = {(int)i,(int)j,(int)k};
 
      }
    }
  }

  h_Transform = new TH2D("h_transform",";r;theta;",50,r_min,r_max,50,theta_min,theta_max);
  for(std::pair<double,double> result : transform)
    h_Transform->Fill(result.first,result.second);

  TCanvas* c = new TCanvas("c","c");
  h_Transform->Draw("colz");
  h_Transform->SetStats(0);
  c->Print(("HT_Plane" + std::to_string(Plane) + ".png").c_str());
  delete c;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Find nearest neighbours to a point

std::vector<int> HoughTransformer::FindNearestNeighbours(int point,const std::vector<double>& channels,const std::vector<double>& ticks,int num) const {

  std::vector<int> point_no_v;
  std::vector<double> dist_v;

  double ch = channels.at(point);
  double ti = ticks.at(point);

  point_no_v.push_back(0);
  dist_v.push_back((ch-channels.at(0))*(ch-channels.at(0)) + (ti-ticks.at(0))*(ti-ticks.at(0)));

  for(size_t i_ch=1;i_ch<channels.size();i_ch++){

    double dist = (ch-channels.at(i_ch))*(ch-channels.at(i_ch)) + (ti-ticks.at(i_ch))*(ti-ticks.at(i_ch));

    bool added = false;
    for(size_t i=0;i<point_no_v.size();i++){
      if(dist < dist_v.at(i)){
        dist_v.insert(dist_v.begin()+i,dist);
        point_no_v.insert(point_no_v.begin()+i,i_ch);
        break;
      }
    }

    if(!added){
      dist_v.push_back(dist);
      point_no_v.push_back(i_ch);
    } 

  }

  std::vector<int> nearest(point_no_v.begin()+1,point_no_v.begin()+num+1);
  std::vector<int> nearest_dist(dist_v.begin()+1,dist_v.begin()+num+1);

  // if most distant point is more than x cm away, return empty vector 
  if(nearest_dist.back() > 5*5) return {};
  else return nearest;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HoughTransformer::MakeTransform2(){

  std::cout << "Attempting Hough transform of " << Channels_v.size() << " points" << std::endl;

  std::vector<std::pair<double,double>> transform;
  double theta_min=3.1415,theta_max=-3.1415;
  double r_min=1e10,r_max=-1e10;

  for(size_t i=0;i<Channels_v.size();i++){
        
    std::vector<int> nearest_neighbors = FindNearestNeighbours(i,Channels_v,Ticks_v,PointGrouping); 

    if(!nearest_neighbors.size()) continue;

    std::vector<double> channels,ticks,widths;
    for(int n : nearest_neighbors){
      channels.push_back(Channels_v.at(n));
      ticks.push_back(Ticks_v.at(n));
      widths.push_back(Widths_v.at(n));
    }

    std::pair<double,double> line = GetLine(channels,ticks);

/*    
    std::cout << std::endl;
    std::cout << Channels_v.at(i) << "  " << Ticks_v.at(i) << std::endl; 
    std::cout << Channels_v.at(nearest_i) << "  " << Ticks_v.at(nearest_i) << std::endl; 
    std::cout << Channels_v.at(second_nearest_i) << "  " << Ticks_v.at(second_nearest_i) << std::endl; 
    std::cout << line.first << "  " << line.second << std::endl;
  */ 

    r_min = std::min(r_min,line.first);
    r_max = std::max(r_max,line.first);
    theta_min = std::min(theta_min,line.second);
    theta_max = std::max(theta_max,line.second);
    transform.push_back(line);

    Transform.push_back(HoughTransformPoint(Plane,line.first,line.second)); 
    Transform.back().Height++;
    Transform.back().Channels = channels;  
    Transform.back().Ticks = ticks;  
    Transform.back().Widths = widths;  
    Transform.back().PointNo = nearest_neighbors;

  }

  std::cout << "theta_min=" << theta_min <<  "  theta_max=" <<  theta_max << std::endl;
  std::cout << "theta bins=" << std::floor((theta_max-theta_min)/ThetaBinSize) << std::endl;

  //h_Transform = new TH2D("h_transform",";r;theta;",50,r_min,r_max,50,-3.1415,3.1415);
  h_Transform = new TH2D("h_transform",";r;theta;",std::floor((r_max-r_min)/RBinSize),r_min,r_max,std::floor((theta_max-theta_min)/ThetaBinSize),theta_min,theta_max);
  for(std::pair<double,double> result : transform)
    h_Transform->Fill(result.first,result.second);

  TCanvas* c = new TCanvas("c","c");
  h_Transform->Draw("colz");
  h_Transform->SetStats(0);
  c->Print(("HT_Plane" + std::to_string(Plane) + ".png").c_str());
  delete c;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<HoughTransformPoint> HoughTransformer::FindPeaks() const {
  
  std::vector<HoughTransformPoint> results;

  // Identify any bins taller than all of their neighbours
  for(int i_bx=1;i_bx<h_Transform->GetNbinsX()+1;i_bx++){
    for(int i_by=1;i_by<h_Transform->GetNbinsY()+1;i_by++){

      if(h_Transform->GetBinContent(i_bx,i_by) < 3) continue;

      double content = h_Transform->GetBinContent(i_bx,i_by);
      bool is_peak = true;
      for(int i=-PeakSize;i<=PeakSize;i++){
        for(int j=-PeakSize;j<=PeakSize;j++){
          if(i == 0 && j == 0) continue;
          if(i_bx+i < 1 || i_bx+i > h_Transform->GetNbinsX()) continue;
          if(i_by+j < 1 || i_by+j > h_Transform->GetNbinsY()) continue;
          if(content < h_Transform->GetBinContent(i_bx+i,i_by+j)) is_peak = false;
        }
      }

      if(is_peak){
        HoughTransformPoint result = HoughTransformPoint(Plane,h_Transform->GetXaxis()->GetBinCenter(i_bx),h_Transform->GetYaxis()->GetBinCenter(i_by));
        result.Height = content;
        bool added  = false;
        for(size_t i=0;i<results.size();i++){
          if(content > results.at(i).Height){
            results.insert(results.begin()+i,result);
            added = true;
            break;
          }
        }
        if(!added) results.push_back(result);
      }

    }
  }

  std::cout << "Found " << results.size() << " peaks" << std::endl;

  // Consolidate the results together, we want to know which data points 
  // contribute to each peak in Hough transform space

  for(HoughTransformPoint& peak : results){

    double r_min = peak.R - h_Transform->GetXaxis()->GetBinWidth(1)*(PeakSize-0.5);
    double r_max = peak.R + h_Transform->GetXaxis()->GetBinWidth(1)*(PeakSize+0.5);
    double theta_min = peak.Theta - h_Transform->GetYaxis()->GetBinWidth(1)*(PeakSize-0.5);
    double theta_max = peak.Theta + h_Transform->GetYaxis()->GetBinWidth(1)*(PeakSize+0.5);

    // Find all HT points within this ROI, make a list of all their points
    for(HoughTransformPoint point : Transform){
      if(point.R < r_min || point.R > r_max || point.Theta < theta_min || point.Theta > theta_max) continue;
      for(size_t i_p=0;i_p<point.PointNo.size();i_p++){
        if(std::find(peak.PointNo.begin(),peak.PointNo.end(),point.PointNo.at(i_p)) != peak.PointNo.end()) continue;
        peak.PointNo.push_back(point.PointNo.at(i_p));
        peak.Channels.push_back(point.Channels.at(i_p));
        peak.Ticks.push_back(point.Ticks.at(i_p));
        peak.Widths.push_back(point.Widths.at(i_p));
      }
    } 

  }

  return results;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HoughTransformer::DrawFits(){

  std::vector<HoughTransformPoint> peaks = FindPeaks();

  std::vector<double> channels,ticks;
  for(size_t i_p=0;i_p<Channels_v.size();i_p++){
    channels.push_back(Channels_v.at(i_p)); 
    ticks.push_back(Ticks_v.at(i_p)); 
  }

  TCanvas* c = new TCanvas("c","c");
 
  TGraphErrors* g = new TGraphErrors(channels.size(),&(channels[0]),&(ticks[0]),0,&(Widths_v[0]));
  g->SetName("graph2D");
  g->SetMarkerColor(1);
  g->SetMarkerSize(0.5);
  g->Draw("AP");

  c->Print(("data_Plane" + std::to_string(Plane) + ".png").c_str());

  std::vector<TF1*> f_v;
  int ctr = 0;
  for(HoughTransformPoint& peak : peaks){

    std::cout << "Drawing peak with r=" << peak.R << "  theta=" << peak.Theta << std::endl;

    f_v.push_back(new TF1(("f" + std::to_string(ctr)).c_str(),"-(cos([1])/sin([1]))*x + [0]/sin([1])",-100,100));
    f_v.back()->SetParameter(0,peak.R);
    f_v.back()->SetParameter(1,peak.Theta);
    f_v.back()->Draw("L same");
    f_v.back()->SetLineColor(ctr+1);

    ctr++;
    //if(ctr > 0) break;

  } 

  c->Print(("test_Plane" + std::to_string(Plane) + ".png").c_str());

  delete c;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<HoughTransformPoint> HoughTransformer::MakeClusters(){

  std::vector<HoughTransformPoint> peaks = FindPeaks();
  TMultiGraph *g = new TMultiGraph();
  std::vector<TGraph*> g_cluster_v;
  int ctr=1;
  for(HoughTransformPoint& peak : peaks){
    g_cluster_v.push_back(new TGraph(peak.Channels.size(),&(peak.Channels[0]),&(peak.Ticks[0])));
    g_cluster_v.back()->SetMarkerColor(ctr);
    g_cluster_v.back()->SetMarkerSize(0.8);
    g_cluster_v.back()->SetMarkerStyle(20);
    g->Add(g_cluster_v.back());
    ctr++;
  }

  TCanvas* c = new TCanvas("c","c");

  g->Draw("AP");
  c->Print(("Clusters_Plane" + std::to_string(Plane) + ".png").c_str());

  for(HoughTransformPoint& peak : peaks)
    RestoreOffset(peak.Channels,peak.Ticks,peak.Widths);

  delete c;
  delete g;

  return peaks;
    
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

}

#endif
