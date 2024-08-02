#ifndef _HoughTransformer_cxx_
#define _HoughTransformer_cxx_

#include "HoughTransformer.h"

namespace hyperonreco { 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

HoughTransformer::HoughTransformer(std::vector<HitLite> hits,int plane,double origin_channel,double origin_tick,bool draw) : 
  Hits(hits),Plane(plane),Origin_Channel(origin_channel),Origin_Tick(origin_tick){

  Func = ROOT::Math::Functor( [&] (const double *coeff ){
      double val = 0.0;
      for(size_t i=0;i<Hits_test.size();i++)
          val += pow(Dist(Hits_test.at(i),coeff[0],coeff[1]),2);
      return val;
  } , 2);
 
  Minimizer = std::unique_ptr<ROOT::Math::Minimizer>( ROOT::Math::Factory::CreateMinimizer( "Minuit2", "Migrad" ) );

  Minimizer->SetMaxFunctionCalls(1000);
  Minimizer->SetTolerance(0.01);
  Minimizer->SetFunction(Func);

  SubtractOffset(Hits);

  system("mkdir -p Plots/");

  // Ensure every hit is uniquely numbered
  for(size_t i_h=0;i_h<Hits.size();i_h++)
    Hits.at(i_h).Number = i_h;
 
  Draw = draw;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

HoughTransformer::~HoughTransformer(){

  if(h_Transform != nullptr) delete h_Transform;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HoughTransformer::SetEvent(int run,int subrun,int event,int pfp){

  Run = run;
  Subrun = subrun;
  Event = event;
  Pfp = pfp;
  RSE = std::to_string(run) + "_" + std::to_string(subrun) + "_" + std::to_string(event) + "_" + std::to_string(Pfp);

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

void HoughTransformer::SetMaxNeighbourDist(double dist){

  MaxNeighbourDist2 = dist*dist;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HoughTransformer::SubtractOffset(std::vector<HitLite>& hits) const {

  for(HitLite& hit : hits){

    hit.Channel = hit.Channel - Origin_Channel;
    hit.Tick = (hit.Tick - Origin_Tick)/TickPerWire;
    hit.Width /= TickPerWire;
  }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HoughTransformer::RestoreOffset(std::vector<HitLite>& hits) const {

  for(HitLite& hit : hits){
    hit.Channel = hit.Channel + Origin_Channel;
    hit.Tick = hit.Tick*TickPerWire + Origin_Tick;
    hit.Width *= TickPerWire;
  }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

double HoughTransformer::Dist(const HitLite& hit,double r, double theta) const {

  return fabs(hit.Channel*cos(theta) + hit.Tick*sin(theta) - r);

}
  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::pair<double,double> HoughTransformer::GetLine(const std::vector<HitLite>& hits){

  Hits_test = hits;

  Minimizer->SetVariable(0,"r",0,0.01);
  Minimizer->SetVariable(1,"theta",0,0.01);
  Minimizer->SetVariableLimits(1,-3.1415,3.1415);
  Minimizer->Minimize();

  return std::make_pair(Minimizer->X()[0],Minimizer->X()[1]);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Find nearest neighbours to a point

std::vector<HitLite> HoughTransformer::FindNearestNeighbours(int point,const std::vector<HitLite>& hits,int num) const {

  if(hits.size()-1 < num) return {}; 

  std::vector<HitLite> ordered_hits;
  std::vector<double> dist_v;

  double ch = hits.at(point).Channel;
  double ti = hits.at(point).Tick;
  int number = hits.at(point).Number;

  for(HitLite hit : hits){

    if(hit.Number == number) continue;

    double dist = (ch-hit.Channel)*(ch-hit.Channel) + (ti-hit.Tick)*(ti-hit.Tick);

    bool added = false;
    for(size_t i=0;i<ordered_hits.size();i++){
      if(dist < dist_v.at(i)){
        dist_v.insert(dist_v.begin()+i,dist);
        ordered_hits.insert(ordered_hits.begin()+i,hit);
        break;
      }
    }

    if(!added){
      dist_v.push_back(dist);
      ordered_hits.push_back(hit);
    } 

  }

  // if most distant point is more than x cm away, return empty vector 
  if(dist_v.at(num-1) > MaxNeighbourDist2) return {};

  std::vector<HitLite> nearest(ordered_hits.begin(),ordered_hits.begin()+num);

  return nearest; 

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HoughTransformer::MakeTransform2(){

  Transform.clear();
  if(h_Transform != nullptr){
    delete h_Transform;
    h_Transform = nullptr;
  }

  std::vector<std::pair<double,double>> transform;
  double theta_min=3.1415,theta_max=-3.1415;
  double r_min=1e10,r_max=-1e10;

  for(size_t i=0;i<Hits.size();i++){
      
    std::vector<HitLite> nearest_neighbors = FindNearestNeighbours(i,Hits,PointGrouping); 

    if(!nearest_neighbors.size()) continue;

    std::pair<double,double> line = GetLine(nearest_neighbors);
     
    /*
    std::cout << std::endl;
    for(HitLite hit : nearest_neighbors) std::cout << hit.Channel << "  " << hit.Tick << std::endl;
    std::cout << line.first << "  " << line.second << std::endl;
    */

    r_min = std::min(r_min,line.first);
    r_max = std::max(r_max,line.first);
    theta_min = std::min(theta_min,line.second);
    theta_max = std::max(theta_max,line.second);
    transform.push_back(line);

    Transform.push_back(HoughTransformPoint(Plane,line.first,line.second)); 
    Transform.back().Height++;
    Transform.back().Hits = nearest_neighbors;;

  }

  if(!Transform.size()) return;

  h_Transform = new TH2D("h_transform",";r;theta;",std::floor((r_max-r_min)/RBinSize),r_min,r_max,std::floor((theta_max-theta_min)/ThetaBinSize),theta_min,theta_max);
  for(std::pair<double,double> result : transform)
    h_Transform->Fill(result.first,result.second);
/*
  if(Draw){
    TCanvas* c = new TCanvas("c","c");
    h_Transform->Draw("colz");
    h_Transform->SetStats(0);
    c->Print(("Plots/Event_" + RSE + "_HT_Plane" + std::to_string(Plane) + "_Tune_" + std::to_string(TuneID) + ".png").c_str());
    delete c;
  }  
*/
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<HoughTransformPoint> HoughTransformer::FindPeaks() const {
 
  std::vector<HoughTransformPoint> results;

  if(h_Transform == nullptr) return results;

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

      for(HitLite hit : point.Hits){

        // Check hit hasn't already been saved 
        bool found = false;
        for(size_t i_h=0;i_h<peak.Hits.size();i_h++)
            if(hit.Number == peak.Hits.at(i_h).Number) found = true;
        if(found) continue;

        peak.Hits.push_back(hit);
      }

    } 

  }

  return results;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HoughTransformer::DrawFits(){

  std::vector<HoughTransformPoint> peaks = FindPeaks();

  std::vector<double> channels,ticks,widths;
  for(HitLite hit : Hits){
    channels.push_back(hit.Channel); 
    ticks.push_back(hit.Tick); 
    widths.push_back(hit.Width);
  }

  TCanvas* c = new TCanvas("c","c");

  TGraphErrors* g = new TGraphErrors(channels.size(),&(channels[0]),&(ticks[0]),0,&(widths[0]));
  g->SetName("graph2D");
  g->SetMarkerColor(1);
  g->SetMarkerSize(0.5);
  g->Draw("AP");

  c->Print(("Plots/Event_" + RSE + "_Data_Plane" + std::to_string(Plane) + "_Tune_" + std::to_string(TuneID) + ".png").c_str());

  std::vector<TF1*> f_v;
  int ctr = 0;
  for(HoughTransformPoint& peak : peaks){

    f_v.push_back(new TF1(("f" + std::to_string(ctr)).c_str(),"-(cos([1])/sin([1]))*x + [0]/sin([1])",-100,100));
    f_v.back()->SetParameter(0,peak.R);
    f_v.back()->SetParameter(1,peak.Theta);
    f_v.back()->Draw("L same");
    f_v.back()->SetLineColor(ctr+1);

    ctr++;
    //if(ctr > 0) break;

  } 

  c->Print(("Plots/Event_" + RSE + "_Fits_Plane" + std::to_string(Plane) + "_Tune_" + std::to_string(TuneID) + ".png").c_str());

  delete c;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<HoughTransformPoint> HoughTransformer::MakeClusters() const {

  std::vector<HoughTransformPoint> peaks = FindPeaks();

  for(HoughTransformPoint& peak : peaks)
    RestoreOffset(peak.Hits);

  if(Draw){

    TMultiGraph *g = new TMultiGraph();
    std::vector<TGraph*> g_cluster_v;
    int ctr=1;
    for(HoughTransformPoint& peak : peaks){

      std::vector<double> channels,ticks,widths;
      for(HitLite hit : peak.Hits){
        channels.push_back(hit.Channel); 
        ticks.push_back(hit.Tick); 
        widths.push_back(hit.Width);
      }

      g_cluster_v.push_back(new TGraph(channels.size(),&(channels[0]),&(ticks[0])));
      g_cluster_v.back()->SetMarkerColor(ctr);
      g_cluster_v.back()->SetMarkerSize(0.8);
      g_cluster_v.back()->SetMarkerStyle(20);
      g->Add(g_cluster_v.back());
      ctr++;
      //if(ctr > 4) break;
    }

    TCanvas* c = new TCanvas("c","c");

    g->Draw("AP");
    c->Print(("Plots/Event_" + RSE + "_Clusters_Plane" + std::to_string(Plane) + "_Tune_" + std::to_string(TuneID) + ".png").c_str());

    delete c;
    delete g;

  }

  return peaks;
    
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::pair<double,double> HoughTransformer::GetPerformanceMetrics() const {

  std::vector<HoughTransformPoint> clusters = MakeClusters();

  if(!clusters.size()) return std::make_pair(0,0);

  // Record how many hits belong to each trackid
  std::map<int,int> m_trackid_hits;
  std::map<int,double> m_trackid_completeness;
  for(HitLite hit : Hits){
    if(m_trackid_hits.find(hit.TrackID) == m_trackid_hits.end()){
      m_trackid_hits[hit.TrackID] = 0;
      m_trackid_completeness[hit.TrackID] = 0.0;
    }
    m_trackid_hits.at(hit.TrackID)++;
  }

  // Calculate the mean purity of the clusters
  //std::cout << "Cluster purities:" << std::endl;
  double mean_purity = 0;
  for(HoughTransformPoint cluster : clusters){
      std::pair<int,int> dominant_trackid = cluster.GetDominantTrackID();
      double purity = (double)dominant_trackid.second/cluster.Hits.size();
      //std::cout << "Purity=" << purity << std::endl;
      m_trackid_completeness.at(dominant_trackid.first) = std::max((double)dominant_trackid.second/m_trackid_hits.at(dominant_trackid.first),m_trackid_completeness.at(dominant_trackid.first));
      mean_purity += purity;
  }
  mean_purity /= clusters.size();

  //std::cout << "TrackID Completenesses:" << std::endl;
  std::map<int,double>::iterator it;
  double mean_completeness = 0;
  int ctr = 0;
  for(it = m_trackid_completeness.begin();it != m_trackid_completeness.end();it++){
    if(m_trackid_hits.at(it->first) < 5) continue;
    //std::cout << it->first << "  " << it->second << std::endl;
    ctr++;
    mean_completeness += it->second;
  }

  if(ctr == 0) return std::make_pair(0,0);

  mean_completeness /= ctr;

  //std::cout << "Mean purity = " << mean_purity << " mean completeness = " << mean_completeness << std::endl;

  return std::make_pair(mean_purity,mean_completeness);
 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

}

#endif
