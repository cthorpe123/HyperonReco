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

void HoughTransformer::SetChi2Cut(double chi2){

   Chi2Cut = chi2;

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

  RemoveVerticalLines(); 

  Transform.clear();
  std::vector<std::pair<double,double>> transform;

  for(size_t i=0;i<Hits.size();i++){
      
    std::vector<HitLite> nearest_neighbors = FindNearestNeighbours(i,Hits,PointGrouping); 

    if(!nearest_neighbors.size()) continue;

    std::tuple<double,double,double> fit = GetLine(nearest_neighbors);

    //if(std::get<2>(fit) > Chi2Cut) continue; 

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
    Transform.back().Chi2 = std::get<2>(fit); 

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

std::vector<HoughTransformPoint> HoughTransformer::FindPeaks2() const {

  std::vector<HoughTransformPoint> results;

  //double theta_min=3.1415,theta_max=-3.1415;
  double r_min=1e10,r_max=-1e10;

  for(HoughTransformPoint point : Transform){
    //std::cout << point.R << "  " << point.Theta << std::endl;
    r_min = std::min(r_min,point.R);
    r_max = std::max(r_max,point.R);
  }

  TH2D* h_Transform_Conv = new TH2D("h_Transform_Conv","",20*(r_max-r_min)/RBinSize,r_min,r_max,20*3.1415/ThetaBinSize,0.0,3.1415);

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

  std::vector<HoughTransformPoint> peaks = FindPeaks3();

  // If there are vertical line clusters, add those back in now
  peaks.insert(peaks.end(),VerticalLines.begin(),VerticalLines.end());


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

      std::cout << "New cluster" << std::endl;

      std::vector<double> channels,ticks,widths;
      for(HitLite hit : peak.Hits){
        //std::cout << hit.Channel << std::endl;
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

std::vector<HoughTransformPoint> HoughTransformer::FindPeaks3() const {

  std::vector<HoughTransformPoint> results;

  std::vector<double> r,theta;
  for(HoughTransformPoint point : Transform){
    r.push_back(point.R);
    theta.push_back(point.Theta);
  }

  TGraph* g = new TGraph(r.size(),&(r[0]),&(theta[0]));
  g->SetMarkerSize(0.25);
  g->SetMarkerStyle(20);


  /*
  //double theta_min=3.1415,theta_max=-3.1415;
  double r_min=1e10,r_max=-1e10;

  for(HoughTransformPoint point : Transform){
    //std::cout << point.R << "  " << point.Theta << std::endl;
    r_min = std::min(r_min,point.R);
    r_max = std::max(r_max,point.R);
  }

  TH2D* h_Transform_Conv = new TH2D("h_Transform_Conv","",Multiplier*(r_max-r_min)/RBinSize,r_min-2*RBinSize,r_max+2*RBinSize,Multiplier*3.1415/ThetaBinSize,0.0,3.1415);

  for(int i=1;i<h_Transform_Conv->GetNbinsX()+1;i++){
    for(int j=1;j<h_Transform_Conv->GetNbinsY()+1;j++){

      double r = h_Transform_Conv->GetXaxis()->GetBinCenter(i);
      double theta = h_Transform_Conv->GetYaxis()->GetBinCenter(j);

      int points_inside = 0;
      for(HoughTransformPoint point : Transform){
        if(point.R > r - RBinSize/2 && point.R < r + RBinSize/2 && point.Theta > theta - ThetaBinSize/2 && point.Theta < theta + ThetaBinSize/2)
          points_inside++;        
      }

      h_Transform_Conv->SetBinContent(i,j,points_inside);

    }
  }
 */

  TH2D* h_Transform_Conv = MakeConvHistogram(Transform,true); 

  if(Draw && Verbosity > 0){
    TCanvas* c = new TCanvas("c","c");
    h_Transform_Conv->Draw("colz");
    h_Transform_Conv->SetStats(0);
    g->Draw("P same");
    c->Print(("Plots/Event_" + RSE + "/Event_" + RSEP + "_HT_Conv_Tune_" + std::to_string(TuneID) + ".png").c_str());
    c->Clear();
    delete c;
  }


/*
  TH2D* h_Transform_2Der_X = (TH2D*)h_Transform_Conv->Clone("h_Transform_2Der_X");
  TH2D* h_Transform_2Der_Y = (TH2D*)h_Transform_Conv->Clone("h_Transform_2Der_Y");
  h_Transform_2Der_X->Reset();
  h_Transform_2Der_Y->Reset();

  double x_bin_width = h_Transform_Conv->GetXaxis()->GetBinWidth(1);
  double y_bin_width = h_Transform_Conv->GetYaxis()->GetBinWidth(1);

  // Calculate the second dervivative at every bin
  for(int i=3;i<h_Transform_Conv->GetNbinsX()-1;i++){
    for(int j=3;j<h_Transform_Conv->GetNbinsY()-1;j++){
        double grad_low_x = (h_Transform_Conv->GetBinContent(i,j) - h_Transform_Conv->GetBinContent(i-1,j))/x_bin_width;
        double grad_high_x = (h_Transform_Conv->GetBinContent(i+1,j) - h_Transform_Conv->GetBinContent(i,j))/x_bin_width;
        double grad2_x = (grad_high_x - grad_low_x)/x_bin_width;
        h_Transform_2Der_X->SetBinContent(i,j,abs(grad2_x)); 

        double grad_low_y = (h_Transform_Conv->GetBinContent(i,j) - h_Transform_Conv->GetBinContent(i,j-1))/y_bin_width;
        double grad_high_y = (h_Transform_Conv->GetBinContent(i,j+1) - h_Transform_Conv->GetBinContent(i,j))/y_bin_width;
        double grad2_y = (grad_high_y - grad_low_y)/y_bin_width;
        h_Transform_2Der_Y->SetBinContent(i,j,abs(grad2_y)); 

    }
  }

  if(Draw && Verbosity > 0){
    TCanvas* c = new TCanvas("c","c");
    h_Transform_2Der_X->Draw("colz");
    h_Transform_2Der_X->SetStats(0);
    g->Draw("P same");
    c->Print(("Plots/Event_" + RSE + "/Event_" + RSEP + "_HT_Conv_2Der_X_Tune_" + std::to_string(TuneID) + ".png").c_str());
    c->Clear();

    h_Transform_2Der_Y->Draw("colz");
    h_Transform_2Der_Y->SetStats(0);
    g->Draw("P same");
    c->Print(("Plots/Event_" + RSE + "/Event_" + RSEP + "_HT_Conv_2Der_Y_Tune_" + std::to_string(TuneID) + ".png").c_str());
    c->Clear();

  }
*/

  /*
  TH2D* h_Transform_Grad_X = (TH2D*)h_Transform_Conv->Clone("h_TransformGrad_X");
  TH2D* h_Transform_Grad_Y = (TH2D*)h_Transform_Conv->Clone("h_TransformGrad_Y");
  h_Transform_Grad_X->Reset();
  h_Transform_Grad_Y->Reset();

  // Try calculating graident at every bin
  for(int i=1;i<h_Transform_Conv->GetNbinsX()+1;i++){

    for(int j=1;j<h_Transform_Conv->GetNbinsY()+1;j++){


        double grad_x = h_Transform_Conv->GetBinContent(i+1,j) - h_Transform_Conv->GetBinContent(i-1,j);
        double grad_y = h_Transform_Conv->GetBinContent(i,j+1) - h_Transform_Conv->GetBinContent(i,j-1);
        
        h_Transform_Grad_X->SetBinContent(i,j,grad_x);
        h_Transform_Grad_Y->SetBinContent(i,j,grad_y);

    }
  }

  if(Draw && Verbosity > 0){
    TCanvas* c = new TCanvas("c","c");
    h_Transform_Grad_X->Draw("colz");
    h_Transform_Grad_X->SetStats(0);
    g->Draw("P same");
    c->Print(("Plots/Event_" + RSE + "/Event_" + RSEP + "_HT_Conv_GradX_Tune_" + std::to_string(TuneID) + ".png").c_str());
    c->Clear();

    h_Transform_Grad_Y->Draw("colz");
    h_Transform_Grad_Y->SetStats(0);
    g->Draw("P same");
    c->Print(("Plots/Event_" + RSE + "/Event_" + RSEP + "_HT_Conv_GradY_Tune_" + std::to_string(TuneID) + ".png").c_str());
    c->Clear();

    delete c;
  }
*/

  TH2D* h_Transform_Grad_X = MakeConvHistogramDer(Transform,false);
  TH2D* h_Transform_Grad_Y = MakeConvHistogramDer(Transform,true);

  if(Draw && Verbosity > 0){
    TCanvas* c = new TCanvas("c","c");
    h_Transform_Grad_X->Draw("colz");
    h_Transform_Grad_X->SetStats(0);
    g->Draw("P same");
    c->Print(("Plots/Event_" + RSE + "/Event_" + RSEP + "_HT_Conv_GradX_Tune_" + std::to_string(TuneID) + ".png").c_str());
    c->Clear();

    h_Transform_Grad_Y->Draw("colz");
    h_Transform_Grad_Y->SetStats(0);
    g->Draw("P same");
    c->Print(("Plots/Event_" + RSE + "/Event_" + RSEP + "_HT_Conv_GradY_Tune_" + std::to_string(TuneID) + ".png").c_str());
    c->Clear();

    delete c;
  }

  TH2D* h_Transform_Grad_SignChange = (TH2D*)h_Transform_Grad_X->Clone("h_TransformGrad_SignChange");
  TH2D* h_Transform_Grad_X_SignChange = (TH2D*)h_Transform_Grad_X->Clone("h_TransformGrad_X_SignChange");
  TH2D* h_Transform_Grad_Y_SignChange = (TH2D*)h_Transform_Grad_X->Clone("h_TransformGrad_Y_SignChange");
  h_Transform_Grad_SignChange->Reset();
  h_Transform_Grad_X_SignChange->Reset();
  h_Transform_Grad_Y_SignChange->Reset();

  double _EPSILON_ = 0.001/Multiplier/RBinSize;

  // Find bins in which the sign of the gradient changes - X
  for(int j=1;j<h_Transform_Grad_X->GetNbinsY()+1;j++){
    bool up = false;
    //std::cout << "Starting new group" << std::endl;
    //std::cout << "_EPSILON_=" << _EPSILON_ << std::endl;
    for(int i=1;i<h_Transform_Grad_X->GetNbinsX()+1;i++){

      //std::cout << h_Transform_Grad_X->GetXaxis()->GetBinCenter(i) << "  " << h_Transform_Grad_X->GetYaxis()->GetBinCenter(j) << "   ";
      //std::cout << h_Transform_Grad_X->GetBinContent(i-1,j) << "  ";

      if(abs(h_Transform_Grad_X->GetBinContent(i-1,j)) < _EPSILON_ && h_Transform_Grad_X->GetBinContent(i,j) > _EPSILON_ && !up){
        h_Transform_Grad_X_SignChange->SetBinContent(i,j,1);
        h_Transform_Grad_SignChange->SetBinContent(i,j,1);
        up = true;
        //std::cout << "Found edge 1" ;
      }
      else if(h_Transform_Grad_X->GetBinContent(i-1,j) > _EPSILON_ && h_Transform_Grad_X->GetBinContent(i,j) < -_EPSILON_){
        up = false;
      }
      else if(h_Transform_Grad_X->GetBinContent(i-1,j) < -_EPSILON_ && h_Transform_Grad_X->GetBinContent(i,j) > _EPSILON_){
        h_Transform_Grad_X_SignChange->SetBinContent(i,j,1);
        h_Transform_Grad_SignChange->SetBinContent(i,j,1);
        up = true;
        //std::cout << "Found edge 2" ;
      }
      else if(h_Transform_Grad_X->GetBinContent(i-1,j) < -_EPSILON_ && abs(h_Transform_Grad_X->GetBinContent(i,j)) < _EPSILON_){
        h_Transform_Grad_X_SignChange->SetBinContent(i,j,1);
        h_Transform_Grad_SignChange->SetBinContent(i,j,1);
        up = false;
        //std::cout << "Found edge 3" ;
      }

      //std::cout << std::endl;        
    }
  }

  _EPSILON_ = 0.001/Multiplier/ThetaBinSize;

  // Find bins in which the sign of the gradient changes - Y
  for(int i=1;i<h_Transform_Grad_X->GetNbinsX()+1;i++){
    bool up = false;
    //std::cout << "_EPSILON_=" << _EPSILON_ << std::endl;
    //std::cout << "Starting new group" << std::endl;
    for(int j=1;j<h_Transform_Grad_X->GetNbinsY()+1;j++){

      //std::cout << h_Transform_Grad_Y->GetBinContent(i,j-1) << "  " << h_Transform_Grad_Y->GetBinContent(i,j);

      if(abs(h_Transform_Grad_Y->GetBinContent(i,j-1)) <_EPSILON_ && h_Transform_Grad_Y->GetBinContent(i,j) > _EPSILON_ && !up){
        h_Transform_Grad_Y_SignChange->SetBinContent(i,j,1);
        h_Transform_Grad_SignChange->SetBinContent(i,j,1);
        up = true;
        // std::cout << "Found edge" ;
      }
      else if(h_Transform_Grad_Y->GetBinContent(i,j-1) > _EPSILON_ && h_Transform_Grad_Y->GetBinContent(i,j) < -_EPSILON_){
        up = false;
      }
      else if(h_Transform_Grad_Y->GetBinContent(i,j-1) < -_EPSILON_ && h_Transform_Grad_Y->GetBinContent(i,j) > _EPSILON_){
        h_Transform_Grad_Y_SignChange->SetBinContent(i,j,1);
        h_Transform_Grad_SignChange->SetBinContent(i,j,1);
        up = true;
        // std::cout << "Found edge" ;
      }
      else if(h_Transform_Grad_Y->GetBinContent(i,j-1) < -_EPSILON_ && abs(h_Transform_Grad_Y->GetBinContent(i,j)) < _EPSILON_){
        h_Transform_Grad_Y_SignChange->SetBinContent(i,j,1);
        h_Transform_Grad_SignChange->SetBinContent(i,j,1);
        up = false;
        // std::cout << "Found edge" ;
      }

      //     std::cout << std::endl;        
    }
  }

  TH2D* h_Transform_Grad_SignChange_Flip = (TH2D*)h_Transform_Grad_SignChange->Clone("h_Transform_Grad_SignChange_Flip");
  for(int i=1;i<h_Transform_Grad_X->GetNbinsX()+1;i++){
    for(int j=1;j<h_Transform_Grad_X->GetNbinsY()+1;j++){
      if(h_Transform_Grad_SignChange_Flip->GetBinContent(i,j) > 0) h_Transform_Grad_SignChange_Flip->SetBinContent(i,j,0);
      else h_Transform_Grad_SignChange_Flip->SetBinContent(i,j,1);
    }
  }

  if(Draw && Verbosity > 0){
    TCanvas* c = new TCanvas("c","c",2000,2000);

    h_Transform_Grad_SignChange->Draw("colz");
    h_Transform_Grad_SignChange->SetStats(0);
    g->Draw("P same");
    c->Print(("Plots/Event_" + RSE + "/Event_" + RSEP + "_HT_Conv_Grad_SignChange_Tune_" + std::to_string(TuneID) + ".png").c_str());
    c->Clear();

    h_Transform_Grad_SignChange_Flip->Draw("colz");
    h_Transform_Grad_SignChange_Flip->SetStats(0);
    g->Draw("P same");
    c->Print(("Plots/Event_" + RSE + "/Event_" + RSEP + "_HT_Conv_Grad_SignChange_Flip_Tune_" + std::to_string(TuneID) + ".png").c_str());
    c->Clear();

    h_Transform_Grad_X_SignChange->Draw("colz");
    h_Transform_Grad_X_SignChange->SetStats(0);
    g->Draw("P same");
    c->Print(("Plots/Event_" + RSE + "/Event_" + RSEP + "_HT_Conv_GradX_SignChange_Tune_" + std::to_string(TuneID) + ".png").c_str());
    c->Clear();

    h_Transform_Grad_Y_SignChange->Draw("colz");
    h_Transform_Grad_Y_SignChange->SetStats(0);
    g->Draw("P same");
    c->Print(("Plots/Event_" + RSE + "/Event_" + RSEP + "_HT_Conv_GradY_SignChange_Tune_" + std::to_string(TuneID) + ".png").c_str());
    c->Clear();

    delete c;

  }

  std::vector<std::vector<int>> bins_x,bins_y;

  //std::cout << "Histogram bins = " << h_Transform_Conv->GetNbinsX()*h_Transform_Conv->GetNbinsY() << std::endl;

  // peak height of Guassian kernel from single point
  double min_height = 2*Gaus2D(0,0,RBinSize,ThetaBinSize,0,0)/Transform.size(); 
  //std::cout << "min_height=" << min_height << std::endl;

  // Now find the islands
  while(true){   

    int max_bin_r = -1;
    int max_bin_theta = -1;
    double max_bin_content = 0.0;

    // Find the highest bin in the histogram
    for(int i=1;i<h_Transform_Conv->GetNbinsX()+1;i++){
      for(int j=1;j<h_Transform_Conv->GetNbinsY()+1;j++){
        if(h_Transform_Conv->GetBinContent(i,j) > max_bin_content){
          max_bin_r = i;
          max_bin_theta = j;
          max_bin_content = h_Transform_Conv->GetBinContent(i,j);
        }
      }
    } 

    //std::cout << "max_bin_r=" << max_bin_r << "  max_bin_theta=" << max_bin_theta << "   max_bin_content=" << max_bin_content << std::endl;

    if(max_bin_content < min_height) break;

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
        h_Transform_Grad_SignChange_Flip->SetBinContent(current_bin_x,current_bin_y,0);
        h_Transform_Conv->SetBinContent(current_bin_x,current_bin_y,0);

        // look at each of the four bins surrounding the current one
        // if bin is occupied, and not already part of the cluster, add it

        if(h_Transform_Grad_SignChange_Flip->GetBinContent(current_bin_x+1,current_bin_y) > 0){ 
          BinsAddedThisPass.push_back(std::make_pair(current_bin_x+1,current_bin_y));               
          bins_x.back().push_back(current_bin_x+1);   
          bins_y.back().push_back(current_bin_y);
          h_Transform_Conv->SetBinContent(current_bin_x+1,current_bin_y,0);
          h_Transform_Grad_SignChange_Flip->SetBinContent(current_bin_x+1,current_bin_y,0);
          nfills_this_pass++;   
        }

        if(h_Transform_Grad_SignChange_Flip->GetBinContent(current_bin_x-1,current_bin_y) > 0){ 
          BinsAddedThisPass.push_back(std::make_pair(current_bin_x-1,current_bin_y));               
          bins_x.back().push_back(current_bin_x-1);   
          bins_y.back().push_back(current_bin_y);
          h_Transform_Conv->SetBinContent(current_bin_x-1,current_bin_y,0);
          h_Transform_Grad_SignChange_Flip->SetBinContent(current_bin_x-1,current_bin_y,0);
          nfills_this_pass++;   
        }

        if(h_Transform_Grad_SignChange_Flip->GetBinContent(current_bin_x,current_bin_y+1) > 0){ 
          BinsAddedThisPass.push_back(std::make_pair(current_bin_x,current_bin_y+1));               
          bins_x.back().push_back(current_bin_x);   
          bins_y.back().push_back(current_bin_y+1);
          h_Transform_Conv->SetBinContent(current_bin_x,current_bin_y+1,0);
          h_Transform_Grad_SignChange_Flip->SetBinContent(current_bin_x,current_bin_y+1,0);
          nfills_this_pass++;   
        }

        if(h_Transform_Grad_SignChange_Flip->GetBinContent(current_bin_x,current_bin_y-1) > 0){ 
          BinsAddedThisPass.push_back(std::make_pair(current_bin_x,current_bin_y-1));               
          bins_x.back().push_back(current_bin_x);   
          bins_y.back().push_back(current_bin_y-1);
          h_Transform_Conv->SetBinContent(current_bin_x,current_bin_y-1,0);
          h_Transform_Grad_SignChange_Flip->SetBinContent(current_bin_x,current_bin_y-1,0);
          nfills_this_pass++;   
        }

        // If the bin is the first/last in theta space, then the last/first bin in theta space also neighbours it as theta
        // has to wrap around
        /*
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
        */

      }

      //std::cout << BinsAddedThisPass.size() << " " << nfills_this_pass << std::endl;

      BinsAddedLastPass = BinsAddedThisPass;

    }

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

  delete h_Transform_Conv;
  delete h_Transform_Grad_X;
  delete h_Transform_Grad_Y;
  delete h_Transform_Grad_X_SignChange;
  delete h_Transform_Grad_Y_SignChange;

  return results;
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

TH2D* HoughTransformer::MakeConvHistogram(std::vector<HoughTransformPoint> transform,bool gausskernel) const {

  double r_min=1e10,r_max=-1e10;
  for(HoughTransformPoint point : transform){
    //std::cout << point.R << "  " << point.Theta << std::endl;
    r_min = std::min(r_min,point.R);
    r_max = std::max(r_max,point.R);
  }

  TH2D* h_Transform_Conv = new TH2D("h_Transform_Conv","",Multiplier*(r_max-r_min)/RBinSize,r_min-2*RBinSize,r_max+2*RBinSize,Multiplier*3.1415/ThetaBinSize,0.0,3.1415);

  if(!gausskernel){
    for(int i=1;i<h_Transform_Conv->GetNbinsX()+1;i++){
      for(int j=1;j<h_Transform_Conv->GetNbinsY()+1;j++){

        double r = h_Transform_Conv->GetXaxis()->GetBinCenter(i);
        double theta = h_Transform_Conv->GetYaxis()->GetBinCenter(j);

        int points_inside = 0;
        for(HoughTransformPoint point : transform){
          if(point.R > r - RBinSize/2 && point.R < r + RBinSize/2 && point.Theta > theta - ThetaBinSize/2 && point.Theta < theta + ThetaBinSize/2)
            points_inside++;        
        }

        h_Transform_Conv->SetBinContent(i,j,points_inside);

      }
    }
  }
  else {

    for(HoughTransformPoint point : transform){
      const double r = point.R;
      const double theta = point.Theta;
      //double width = 1+2*point.Chi2;
  //    std::cout << r << "  " << theta << std::endl;
      for(int i=1;i<h_Transform_Conv->GetNbinsX()+1;i++){
        for(int j=1;j<h_Transform_Conv->GetNbinsY()+1;j++){
//          std::cout << r << "  " << theta << "     " << h_Transform_Conv->GetXaxis()->GetBinCenter(i) << "  " << h_Transform_Conv->GetYaxis()->GetBinCenter(j) << "    " <<  Gaus2D(r,theta,RBinSize,ThetaBinSize,h_Transform_Conv->GetXaxis()->GetBinCenter(i), h_Transform_Conv->GetYaxis()->GetBinCenter(j)) << std::endl;
        double bin_r = h_Transform_Conv->GetXaxis()->GetBinCenter(i);
        double bin_theta = h_Transform_Conv->GetYaxis()->GetBinCenter(j);

        h_Transform_Conv->Fill(bin_r,bin_theta,Gaus2D(r,theta,RBinSize,ThetaBinSize,bin_r,bin_theta));
        }
      }
    }

  }

  h_Transform_Conv->Scale(1.0/transform.size());

  return h_Transform_Conv;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

TH2D* HoughTransformer::MakeConvHistogramDer(std::vector<HoughTransformPoint> transform,bool xy) const {

  double r_min=1e10,r_max=-1e10;
  for(HoughTransformPoint point : transform){
    //std::cout << point.R << "  " << point.Theta << std::endl;
    r_min = std::min(r_min,point.R);
    r_max = std::max(r_max,point.R);
  }

  std::string name = xy ? "h_Transform_Conv_Y" : "h_Transform_Conv_X";
  TH2D* h_Transform_Conv = new TH2D(name.c_str(),"",Multiplier*(r_max-r_min)/RBinSize,r_min-2*RBinSize,r_max+2*RBinSize,Multiplier*3.1415/ThetaBinSize,0.0,3.1415);
  for(HoughTransformPoint point : transform){
    const double r = point.R;
    const double theta = point.Theta;
    //double width = 1+2*point.Chi2;
    //    std::cout << r << "  " << theta << std::endl;
    for(int i=1;i<h_Transform_Conv->GetNbinsX()+1;i++){
      for(int j=1;j<h_Transform_Conv->GetNbinsY()+1;j++){
        //          std::cout << r << "  " << theta << "     " << h_Transform_Conv->GetXaxis()->GetBinCenter(i) << "  " << h_Transform_Conv->GetYaxis()->GetBinCenter(j) << "    " <<  Gaus2D(r,theta,RBinSize,ThetaBinSize,h_Transform_Conv->GetXaxis()->GetBinCenter(i), h_Transform_Conv->GetYaxis()->GetBinCenter(j)) << std::endl;
        double bin_r = h_Transform_Conv->GetXaxis()->GetBinCenter(i);
        double bin_theta = h_Transform_Conv->GetYaxis()->GetBinCenter(j);

        if(!xy) h_Transform_Conv->Fill(bin_r,bin_theta,Gaus2D_GradX(r,theta,RBinSize,ThetaBinSize,bin_r,bin_theta));
        else h_Transform_Conv->Fill(bin_r,bin_theta,Gaus2D_GradY(r,theta,RBinSize,ThetaBinSize,bin_r,bin_theta));
      }
    }
  }

  h_Transform_Conv->Scale(1.0/transform.size());

  return h_Transform_Conv;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HoughTransformer::RemoveVerticalLines(){

  std::map<double,std::vector<HitLite>> m_ch_ticks; 
  for(HitLite hit : Hits){
    //std::cout << hit.Channel << "   " << hit.Tick << std::endl; 
    if(m_ch_ticks.find(hit.Channel) == m_ch_ticks.end()) m_ch_ticks[hit.Channel] = {hit};
    else m_ch_ticks[hit.Channel].push_back(hit);
  } 

  std::map<double,std::vector<HitLite>>::iterator it;
  std::vector<size_t> bad_hits;
  for(it = m_ch_ticks.begin();it != m_ch_ticks.end();it++){
    std::sort(it->second.begin(),it->second.end(),SortHits);
    //std::cout << it->first << "  " << it->second.size() << std::endl;
    if(it->second.size() > 3){
      //std::cout << "Possible vertical line:" << std::endl;
      for(size_t i=0;i<it->second.size()-1;i++){
        //std::cout << it->second.at(i+1).Tick - it->second.at(i).Tick << std::endl;
        if(it->second.at(i+1).Tick - it->second.at(i).Tick <= 18/TickPerWire){
          bad_hits.push_back(it->second.at(i).Number);
          if(i == it->second.size()-2) bad_hits.push_back(it->second.at(i+1).Number); 
        }
      }        
    }
  }

  // TODO: atm can only clump all vertical hits into one cluster - should implement
  // something that can handle multiple
  VerticalLines.resize(1);

  // erase the vertical lines from the hit vector - write new alg for clustering them tomorrow
  for(size_t i=0;i<bad_hits.size();i++){
    for(size_t i_h=0;i_h<Hits.size();i_h++)
      if(bad_hits.at(i) == Hits.at(i_h).Number){
        VerticalLines.back().Hits.push_back(Hits.at(i_h));
        Hits.erase(Hits.begin()+i_h);
        break;
      }
  } 

//  for(HoughTransformPoint& point : VerticalLines)
//    SubtractOffset(point.Hits);
 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

}

#endif
