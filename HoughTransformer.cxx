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
      //std::cout << coeff[0] << "  " << coeff[1] << "  " << val << std::endl;      
      return val;
  } , 2);
 
  Minimizer = std::unique_ptr<ROOT::Math::Minimizer>( ROOT::Math::Factory::CreateMinimizer( "Minuit2", "Migrad" ) );

  Minimizer->SetMaxFunctionCalls(5000);
  Minimizer->SetTolerance(0.001);
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
  RSE = std::to_string(run) + "_" + std::to_string(subrun) + "_" + std::to_string(event);
  RSEP = std::to_string(run) + "_" + std::to_string(subrun) + "_" + std::to_string(event) + "_" + std::to_string(Pfp) + "_Plane" + std::to_string(Plane);

  system(("mkdir -p Plots/Event_" + RSE).c_str());

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

void HoughTransformer::SetConvFloor(int fl){

  ConvFloor = fl;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HoughTransformer::SetMaxNeighbourDist(double dist){

  MaxNeighbourDist2 = dist*dist;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HoughTransformer::SetVerbosity(int verb){

  Verbosity = verb;

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

std::tuple<double,double,double> HoughTransformer::GetLine(const std::vector<HitLite>& hits){

  Hits_test = hits;

  Minimizer->SetVariable(0,"r",0,0.01);
  Minimizer->SetVariable(1,"theta",0,0.01);
  Minimizer->SetVariableLimits(1,-3.1415,3.1415);
  Minimizer->Minimize();

  //if(Minimizer->Status()) std::cout << "Fit failed" << std::endl; 

  return std::make_tuple(Minimizer->X()[0],Minimizer->X()[1],Minimizer->MinValue());

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Find nearest neighbours to a point

std::vector<HitLite> HoughTransformer::FindNearestNeighbours(int point,const std::vector<HitLite>& hits,size_t num) const {

  if(hits.size()-1 < num) return {}; 

  std::vector<HitLite> ordered_hits;
  std::vector<double> dist_v;

  double ch = hits.at(point).Channel;
  double ti = hits.at(point).Tick;
  size_t number = hits.at(point).Number;

  ordered_hits.push_back(hits.at(point));
  dist_v.push_back(0.0);

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
  if(dist_v.at(num) > MaxNeighbourDist2) return {};

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

  for(size_t i=0;i<Hits.size();i++){
      
    std::vector<HitLite> nearest_neighbors = FindNearestNeighbours(i,Hits,PointGrouping); 

    if(!nearest_neighbors.size()) continue;

    std::tuple<double,double,double> fit = GetLine(nearest_neighbors);
    std::pair<double,double> line = std::make_pair(std::get<0>(fit),std::get<1>(fit));

    //std::cout << std::get<2>(fit) << std::endl;


    //std::cout << "Before: " << line.first << "  " << line.second << "  a=" << -(cos(line.second)/sin(line.second)) << " b=" << line.first/sin(line.second) << std::endl;

    if(line.first < 0){
      line.first = abs(line.first); 
      if(line.second < 0) line.second += 3.1415;
      else if(line.second > 0) line.second -= 3.1415;
    }

    //std::cout << "After1: " << line.first << "  " << line.second << "  a=" << -(cos(line.second)/sin(line.second)) << " b=" << line.first/sin(line.second) << std::endl;

    //if(line.second < 3.1415) line.second += 2*3.1415;
    //if(line.second > 3.1415) line.second -= 2*3.1415;

    if(line.second < 0.0) line.second += 3.1415;

    //std::cout << "After2: " << line.first << "  " << line.second << "  a=" << -(cos(line.second)/sin(line.second)) << " b=" << line.first/sin(line.second) << std::endl;
 

    //std::cout << line.first << "  " << line.second << std::endl;
 
    std::vector<double> channel,tick;
    for(HitLite hit : nearest_neighbors){
      // std::cout << hit.Channel << "  " << hit.Tick << std::endl;
      channel.push_back(hit.Channel);
      tick.push_back(hit.Tick);
    }
   
    if(Draw && Verbosity == 2){

      TCanvas* c = new TCanvas("c","c");
      TGraph* g = new TGraph(channel.size(),&(channel[0]),&(tick[0]));
      TF1* f = new TF1("f","-(cos([1])/sin([1]))*x + [0]/sin([1])",-1000,1000);
      f->SetParameter(0,line.first);
      f->SetParameter(1,line.second);
      g->SetMarkerSize(0.4);
      g->SetMarkerStyle(20);
      g->Draw("AP");
      f->Draw("L same"); 

      TLegend* l = new TLegend(0.8,0.8,1.0,1.0);        
      l->AddEntry((TObject*)0,("R=" + std::to_string(line.first)).c_str(),"");
      l->AddEntry((TObject*)0,("Theta=" + std::to_string(line.second)).c_str(),"");
      l->AddEntry((TObject*)0,("Chi2=" + std::to_string(std::get<2>(fit))).c_str(),"");
      l->Draw();

      c->Print(("Plots/Event_" + RSE + "/Event_" + RSEP + "_Group_" + std::to_string(i) + ".png").c_str());

      delete g;
      delete f;
      delete c;
      delete l;

    }

    line.first = sqrt(line.first);

    transform.push_back(line);
    Transform.push_back(HoughTransformPoint(Plane,line.first,line.second)); 
    Transform.back().Height++;
    Transform.back().Hits = nearest_neighbors;;

  }

  if(Transform.size() < 4) return;

/*
  int nbins_r = std::max((r_max-r_min)/RBinSize,1.0);
  int nbins_theta = std::max((theta_max-theta_min)/ThetaBinSize,1.0);
  h_Transform = new TH2D("h_transform",";r;theta;",nbins_r,r_min,r_max,nbins_theta,theta_min,theta_max);
  for(std::pair<double,double> result : transform)
    h_Transform->Fill(result.first,result.second);
*/

  if(Draw && Verbosity > 0){
    TCanvas* c = new TCanvas("c","c");
 
    std::vector<double> r,theta;
    for(HoughTransformPoint point : Transform){
      r.push_back(point.R);
      theta.push_back(point.Theta);
    }

    TGraph* g = new TGraph(r.size(),&(r[0]),&(theta[0]));
    g->SetMarkerSize(0.25);
    g->SetMarkerStyle(20);
    g->Draw("AP");
    c->Print(("Plots/Event_" + RSE + "/Event_" + RSEP + "_HT_Scatter_Tune_" + std::to_string(TuneID) + ".png").c_str());

    delete c;
  }  


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

    double r_min = peak.R - h_Transform->GetXaxis()->GetBinWidth(1)*(PeakSize);
    double r_max = peak.R + h_Transform->GetXaxis()->GetBinWidth(1)*(PeakSize);
    double theta_min = peak.Theta - h_Transform->GetYaxis()->GetBinWidth(1)*(PeakSize);
    double theta_max = peak.Theta + h_Transform->GetYaxis()->GetBinWidth(1)*(PeakSize);

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

std::vector<HoughTransformPoint> HoughTransformer::FindPeaks2() const {

  std::vector<HoughTransformPoint> results;

  //double theta_min=3.1415,theta_max=-3.1415;
  double r_min=1e10,r_max=-1e10;

  for(HoughTransformPoint point : Transform){
    //std::cout << point.R << "  " << point.Theta << std::endl;
    r_min = std::min(r_min,point.R);
    r_max = std::max(r_max,point.R);
  }


  TH2D* h_Transform_Conv = new TH2D("h_Transform_Conv","",10*(r_max-r_min)/RBinSize,r_min,r_max,10*3.1415/ThetaBinSize,0.0,3.1415);

  for(int i=1;i<h_Transform_Conv->GetNbinsX()+1;i++){
    for(int j=1;j<h_Transform_Conv->GetNbinsY()+1;j++){

      double r = h_Transform_Conv->GetXaxis()->GetBinCenter(i);
      double theta = h_Transform_Conv->GetYaxis()->GetBinCenter(j);

      int points_inside = 0;
      for(HoughTransformPoint point : Transform){
        if(point.R > r - RBinSize/2 && point.R < r + RBinSize/2 && point.Theta > theta - ThetaBinSize/2 && point.Theta < theta + ThetaBinSize/2)
          points_inside++;        
      }
      if(points_inside > ConvFloor) h_Transform_Conv->SetBinContent(i,j,points_inside);
    }
  }

  if(Draw && Verbosity > 0){
    TCanvas* c = new TCanvas("c","c");
    h_Transform_Conv->Draw("colz");
    h_Transform_Conv->SetStats(0);
    c->Print(("Plots/Event_" + RSE + "/Event_" + RSEP + "_HT_Conv_Tune_" + std::to_string(TuneID) + ".png").c_str());
    c->Clear();
    delete c;
  }

  std::vector<std::vector<int>> bins_x,bins_y;

  while(true){   

    int max_bin_r = -1;
    int max_bin_theta = -1;
    int max_bin_content = 2;

    // Find the highest bin in the histogram
    for(int i=1;i<h_Transform_Conv->GetNbinsX()+1;i++){
      for(int j=1;j<h_Transform_Conv->GetNbinsY()+1;j++){
        //std::cout << h_Transform_Conv->GetBinContent(i,j) << std::endl;
        if(h_Transform_Conv->GetBinContent(i,j) > max_bin_content){
          max_bin_r = i;
          max_bin_theta = j;
          max_bin_content = h_Transform_Conv->GetBinContent(i,j);
        }
      }
    } 

    if(max_bin_content <= 2) break;

    bins_x.push_back(std::vector<int>());
    bins_y.push_back(std::vector<int>());

    results.push_back(HoughTransformPoint(Plane,max_bin_r,max_bin_theta)); 

    // Find contiguous regions of filled bins
    std::vector<std::pair<int,int>> BinsAddedLastPass;
    std::vector<std::pair<int,int>> BinsAddedThisPass;

    BinsAddedLastPass.push_back(std::make_pair(max_bin_r,max_bin_theta));

    int nfills_this_pass = 1;
    while(nfills_this_pass > 0){

      nfills_this_pass = 0;
      BinsAddedThisPass.clear();

      // iterate over the bins added in the last pass, check the bins that neighbour those
      for(size_t i_b=0;i_b<BinsAddedLastPass.size();i_b++){

        int current_bin_x = BinsAddedLastPass.at(i_b).first;
        int current_bin_y = BinsAddedLastPass.at(i_b).second;
        h_Transform_Conv->SetBinContent(current_bin_x,current_bin_y,1);

        // look at each of the eight bins surrounding the current one
        // if bin is occupied, and not already part of the cluster, add it
        for(int i=-1;i<=1;i++){
          for(int j=-1;j<=1;j++){
            if(i == 0  && j == 0) continue;
            if(h_Transform_Conv->GetBinContent(current_bin_x+i,current_bin_y+j) > 2){ 
              BinsAddedThisPass.push_back(std::make_pair(current_bin_x+i,current_bin_y+j));               
              bins_x.back().push_back(current_bin_x+i);   
              bins_y.back().push_back(current_bin_y+j);
              h_Transform_Conv->SetBinContent(current_bin_x+i,current_bin_y+j,1);
              nfills_this_pass++;   
            }
          }
        }

        // If the bin is the first/last in theta space, then the last/first bin in theta space also neighbours it as theta
        // has to wrap around
        for(int i=-1;i<=1;i++){
          if(current_bin_y == h_Transform_Conv->GetNbinsY() && h_Transform_Conv->GetBinContent(current_bin_x+i,1) > 2){
            BinsAddedThisPass.push_back(std::make_pair(current_bin_x+i,1));               
            bins_x.back().push_back(current_bin_x+i);   
            bins_y.back().push_back(1);
            h_Transform_Conv->SetBinContent(current_bin_x+i,1,1);
            nfills_this_pass++;   
          } 
          if(current_bin_y == 1 && h_Transform_Conv->GetBinContent(current_bin_x+i,h_Transform_Conv->GetNbinsY()) > 2){
            BinsAddedThisPass.push_back(std::make_pair(current_bin_x+i,h_Transform_Conv->GetNbinsY()));               
            bins_x.back().push_back(current_bin_x+i);   
            bins_y.back().push_back(h_Transform_Conv->GetNbinsY());
            h_Transform_Conv->SetBinContent(current_bin_x+i,h_Transform_Conv->GetNbinsY(),1);
            nfills_this_pass++;   
          } 
        }

      }

      BinsAddedLastPass = BinsAddedThisPass;

    }

  }

  if(Draw && Verbosity > 0){
    TCanvas* c = new TCanvas("c","c");
    h_Transform_Conv->Draw("colz");
    h_Transform_Conv->SetStats(0);
    c->Print(("Plots/Event_" + RSE +"/Event_" + RSEP + "_HT_Islands_Tune_" + std::to_string(TuneID) + ".png").c_str());
    c->Clear();
    delete c;
  }

  // For each peak, find all the points inside it, and get their corresponding hits
  
  std::vector<int> used_hits;
  for(HoughTransformPoint point : Transform){

    double r = point.R;
    double theta = point.Theta;
    int bin_x = h_Transform_Conv->GetXaxis()->FindBin(r);
    int bin_y = h_Transform_Conv->GetYaxis()->FindBin(theta);

    for(size_t p=0;p<bins_x.size();p++){
      if(std::find(bins_x.at(p).begin(),bins_x.at(p).end(),bin_x) != bins_x.at(p).end() &&
          std::find(bins_y.at(p).begin(),bins_y.at(p).end(),bin_y) != bins_y.at(p).end()){
        for(HitLite hit : point.Hits){
          if(std::find(used_hits.begin(),used_hits.end(),hit.Number) != used_hits.end()){
            //std::cout << "Hit already used in another cluster" << std::endl;
          }
          else {
            results.at(p).Hits.push_back(hit);
            used_hits.push_back(hit.Number);
          }
        }
        break; 
      }
    }

  } 

  // Remove any peaks with not hits
  // TODO: Why are there some peaks with no hits?
  std::vector<HoughTransformPoint> results_tmp;  
  for(HoughTransformPoint result : results)
      if(result.Hits.size()) results_tmp.push_back(result);
  results = results_tmp;

  for(HoughTransformPoint& result : results){
    result.RemoveDuplicateHits();
  } 

/*
  std::cout << "Lists of hits in each cluster" << std::endl;
  int ctr=0;
  used_hits.clear();
  for(HoughTransformPoint result : results){
    std::cout << "Cluster " << ctr << std::endl;
    for(HitLite hit : result.Hits){
      std::cout << hit.Number << std::endl;
      if(std::find(used_hits.begin(),used_hits.end(),hit.Number) != used_hits.end()){
        std::cout << "Hit already used in another cluster" << std::endl;
      }
      else used_hits.push_back(hit.Number);
    }
    ctr++;
  }
*/

  delete h_Transform_Conv; 

  return results;

} 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HoughTransformer::DrawFits(){

  std::vector<HoughTransformPoint> peaks = FindPeaks2();

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

  c->Print(("Plots/Event_" + RSE + "/Event_" + RSEP + "_Data_Tune_" + std::to_string(TuneID) + ".png").c_str());

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

  c->Print(("Plots/Event_" + RSE + "/Event_" + RSEP + "_Fits_Tune_" + std::to_string(TuneID) + ".png").c_str());

  delete c;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<HoughTransformPoint> HoughTransformer::MakeClusters() const {

  std::vector<HoughTransformPoint> peaks = FindPeaks2();

  if(Draw){

    TMultiGraph *g = new TMultiGraph();
    std::vector<TGraph*> g_cluster_v;

   std::vector<double> channels_allhits,ticks_allhits;
    for(HitLite hit : Hits){
        channels_allhits.push_back(hit.Channel);
        ticks_allhits.push_back(hit.Tick);
    }

    TGraph* g_allhits = new TGraph(channels_allhits.size(),&(channels_allhits[0]),&(ticks_allhits[0]));
    g_allhits->SetMarkerColor(1);
    g_allhits->SetMarkerSize(0.6);
    g_allhits->SetMarkerStyle(20);
    g->Add(g_allhits);

    int ctr=1;
    for(HoughTransformPoint& peak : peaks){

      std::vector<double> channels,ticks,widths;
      for(HitLite hit : peak.Hits){
        channels.push_back(hit.Channel); 
        ticks.push_back(hit.Tick); 
        widths.push_back(hit.Width);
      }

      g_cluster_v.push_back(new TGraph(channels.size(),&(channels[0]),&(ticks[0])));
      g_cluster_v.back()->SetMarkerColor(ctr+1);
      g_cluster_v.back()->SetMarkerSize(0.8);
      g_cluster_v.back()->SetMarkerStyle(20);
      g->Add(g_cluster_v.back());
      ctr++;
      //if(ctr > 4) break;
    }

    TCanvas* c = new TCanvas("c","c");

    g->Draw("AP");
    c->Print(("Plots/Event_" + RSE + "/Event_" + RSEP + "_Clusters_Tune_" + std::to_string(TuneID) + ".png").c_str());

    delete c;
    delete g;

  }

  for(HoughTransformPoint& peak : peaks)
    RestoreOffset(peak.Hits);

  return peaks;
    
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::pair<double,double> HoughTransformer::GetPerformanceMetrics() const {

  std::vector<HoughTransformPoint> clusters = MakeClusters();

  if(!clusters.size()) return std::make_pair(0,0);

  size_t hits = 0;
  for(HoughTransformPoint cluster : clusters){
    hits += cluster.Hits.size();
  }

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
  double mean_purity = 0;
  for(HoughTransformPoint cluster : clusters){
      std::pair<int,int> dominant_trackid = cluster.GetDominantTrackID();
      double purity = (double)dominant_trackid.second/cluster.Hits.size();
      //std::cout << "Purity=" << purity << std::endl;
      m_trackid_completeness.at(dominant_trackid.first) = std::max((double)dominant_trackid.second/m_trackid_hits.at(dominant_trackid.first),m_trackid_completeness.at(dominant_trackid.first));
      mean_purity += purity;
  }
  mean_purity /= clusters.size();

  std::map<int,double>::iterator it;
  double mean_completeness = 0;
  int ctr = 0;
  for(it = m_trackid_completeness.begin();it != m_trackid_completeness.end();it++){
    if(m_trackid_hits.at(it->first) < 5) continue;
    ctr++;
    mean_completeness += it->second;
  }

  if(ctr == 0) return std::make_pair(0,0);

  mean_completeness /= ctr;

  return std::make_pair(mean_purity,mean_completeness);
 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

}

#endif
