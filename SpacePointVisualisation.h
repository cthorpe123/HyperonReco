#include "TGraph2D.h"
#include "RecoParticle.h"
#include "VFitter.h"
#include "FittedV.h"
#include "Position_To_Wire.h"

namespace hyperonreco {

const int Markers[5] = {20,21,22,29,34};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

TGraph2D* MakeGraph(std::vector<double> X,std::vector<double> Y,std::vector<double> Z,int color,int marker){

  TGraph2D* g = new TGraph2D(X.size(),&(X[0]),&(Y[0]),&(Z[0]));
  g->SetName(("graph_"+std::to_string(color)).c_str());
  g->SetMarkerColor(color);
  g->SetMarkerStyle(Markers[marker]);
  g->SetMarkerSize(1.0);

  return g;

} 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

TGraph* Make2DGraph(std::vector<double> channel,std::vector<double> tick,std::vector<double> width,int color,int marker){

  TGraphErrors* g = new TGraphErrors(channel.size(),&(channel[0]),&(tick[0]),0,&(width[0]));
  g->SetName(("graph2D_"+std::to_string(color)).c_str());
  g->SetMarkerColor(color);
  g->SetMarkerStyle(Markers[marker]);
  g->SetMarkerSize(0.5);

  return g;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::pair<TGraph*,TGraph*> Make2DVGraphs(const FittedV& v,int i_pl,int color){

   TGraph* g1 = new TGraph(2);
   TGraph* g2 = new TGraph(2);
   g1->SetName("fit_graph2D_1");
   g2->SetName("fit_graph2D_2");
   LineWireTick2 line_1 = v.GetArm1_2D2(i_pl);
   LineWireTick2 line_2 = v.GetArm2_2D2(i_pl);

   g1->SetPoint(0,line_1.Start.X(),line_1.Start.Y());
   g2->SetPoint(0,line_2.Start.X(),line_2.Start.Y());

   g1->SetPoint(1,line_1.Start.X()+line_1.Direction.X()*200,line_1.Start.Y()+line_1.Direction.Y()*200);
   g2->SetPoint(1,line_2.Start.X()+line_2.Direction.X()*200,line_2.Start.Y()+line_2.Direction.Y()*200);

   g1->SetLineColor(color); 
   g2->SetLineColor(color); 

   return std::make_pair(g1,g2);
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DrawGraphs(const std::vector<TGraph2D*>& graph_v,std::pair<TGraph2D*,TGraph2D*> fit={nullptr,nullptr}){

  double min_x=1e6,max_x=-1e6;
  double min_y=1e6,max_y=-1e6;
  double min_z=1e6,max_z=-1e6;

  if(!graph_v.size()){
     std::cout << "Event contains no PFParticles, not drawing" << std::endl;
     return;
  }
      
  for(TGraph2D* g : graph_v){
    for(int i_p=0;i_p<g->GetN();i_p++){
      min_x = std::min(min_x,g->GetX()[i_p]);         
      max_x = std::max(max_x,g->GetX()[i_p]);         
      min_y = std::min(min_y,g->GetY()[i_p]);         
      max_y = std::max(max_y,g->GetY()[i_p]);         
      min_z = std::min(min_z,g->GetZ()[i_p]);         
      max_z = std::max(max_z,g->GetZ()[i_p]);         
    }       
  }
/*
  std::cout << "Box dimensions:" << std::endl;
  std::cout << "X: " << min_x << " " << max_x << std::endl;
  std::cout << "Y: " << min_y << " " << max_y << std::endl;
  std::cout << "Z: " << min_z << " " << max_z << std::endl;
*/
  TCanvas* c = new TCanvas("C","C",1200,800);

  TH3D* h = new TH3D("h_dummy","",1,min_x,max_x,1,min_y,max_y,1,min_z,max_z);
  h->Draw();
  h->SetStats(0);

  for(size_t i_g=0;i_g<graph_v.size();i_g++)
    graph_v.at(i_g)->Draw("P same");
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Draw2DGraphs(const std::vector<std::vector<TGraph*>> graph_v,FittedV* fit=nullptr){

  for(int i_pl=0;i_pl<3;i_pl++){

    double min_x=1e6,max_x=-1e6;
    double min_y=1e6,max_y=-1e6;

    if(!graph_v.at(i_pl).size()){
      std::cout << "Event contains no PFParticles, not drawing" << std::endl;
      return;
    }

    std::cout << "Event contains " << graph_v.at(i_pl).size() << " PFParticles" << std::endl;

    bool has_points=false;
    for(TGraph* g : graph_v.at(i_pl)){
      if(g->GetN() > 0) has_points = true;
      for(int i_p=0;i_p<g->GetN();i_p++){
        min_x = std::min(min_x,g->GetX()[i_p]);         
        max_x = std::max(max_x,g->GetX()[i_p]);         
        min_y = std::min(min_y,g->GetY()[i_p]);         
        max_y = std::max(max_y,g->GetY()[i_p]);         
      }       
    }

    if(!has_points) continue;

    std::pair<TGraph*,TGraph*> fit_arms = {nullptr,nullptr};
    if(fit != nullptr){
      fit_arms = Make2DVGraphs(*fit,i_pl,1);
      fit_arms.first->SetLineColor(2);
      fit_arms.second->SetLineColor(3);
      fit_arms.first->Draw("L same");
      fit_arms.second->Draw("L same");
      max_x = std::max(fit_arms.first->GetX()[0],max_x); 
      min_x = std::min(fit_arms.first->GetX()[0],min_x); 
      max_y = std::max(fit_arms.first->GetY()[0],max_y); 
      min_y = std::min(fit_arms.first->GetY()[0],min_y); 
    }

    double width_x = (max_x-min_x)/2;
    double width_y = (max_y-min_y)/2;

    TCanvas* c = new TCanvas("C_2d","C_2d");

    TH2D* h = new TH2D("h_2d_dummy",";Channel;Tick",1,min_x-width_x/10,max_x+width_x/10,1,min_y-width_y/10,max_y+width_y/10);
    h->Draw();
    h->SetStats(0);

    for(size_t i_g=0;i_g<graph_v.at(i_pl).size();i_g++)
      graph_v.at(i_pl).at(i_g)->Draw("P same");

    if(fit != nullptr){
      fit_arms.first->Draw("L same");
      fit_arms.second->Draw("L same");
    }

    c->Print(("Plane"+std::to_string(i_pl)+".png").c_str());

    delete h;
    c->Close();

  }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::pair<TGraph2D*,TGraph2D*> MakeVGraphs(const FittedV& v){

   TGraph2D* g1 = new TGraph2D(2);
   TGraph2D* g2 = new TGraph2D(2);
   
   g1->SetPoint(0,v.Vertex.X(),v.Vertex.Y(),v.Vertex.Z());
   g2->SetPoint(0,v.Vertex.X(),v.Vertex.Y(),v.Vertex.Z());

   TVector3 arm1end = v.Vertex + v.GetArm1Dir()*1000; 
   TVector3 arm2end = v.Vertex + v.GetArm2Dir()*1000; 

   g1->SetPoint(1,arm1end.X(),arm1end.Y(),arm1end.Z());
   g2->SetPoint(1,arm2end.X(),arm2end.Y(),arm2end.Z());

   return std::make_pair(g1,g2);
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

FittedV MakeFittedVGuessShower(const RecoParticle& shower){

  FittedV v;
  v.Vertex = TVector3(shower.X_NoSC,shower.Y_NoSC,shower.Z_NoSC);
  TVector3 showerdir(shower.ShowerDirectionX,shower.ShowerDirectionY,shower.ShowerDirectionZ);
  v.Arm1Theta = showerdir.Theta() + shower.ShowerOpeningAngle/5;
  v.Arm1Phi = showerdir.Phi() + shower.ShowerOpeningAngle/5;
  v.Arm2Theta = showerdir.Theta() - shower.ShowerOpeningAngle/5;
  v.Arm2Phi = showerdir.Phi() - shower.ShowerOpeningAngle/5;

  return v;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

FittedV MakeFittedVGuessTrack(const RecoParticle& track){

  FittedV v;

  v.Vertex = TVector3(track.X_NoSC,track.Y_NoSC,track.Z_NoSC);
  TVector3 trackdir(track.TrackDirectionX,track.TrackDirectionY,track.TrackDirectionZ);
  v.Arm1Theta = trackdir.Theta() + 0.3;
  v.Arm1Phi = trackdir.Phi() + 0.3;
  v.Arm2Theta = trackdir.Theta() - 0.3;
  v.Arm2Phi = trackdir.Phi() - 0.3;

  return v;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void KeepHitsInROI(TVector3 point,std::vector<double>& channels,std::vector<double>& ticks,std::vector<double>& widths,double roi_size_ch,double roi_size_tick,int i_pl){

  std::pair<double,double> vertex_ch_tick = WireTick(point,i_pl); 

  std::vector<double> channels_tmp;
  std::vector<double> ticks_tmp;
  std::vector<double> widths_tmp;

  for(size_t i_h=0;i_h<channels.size();i_h++){
    if(abs(channels.at(i_h) - vertex_ch_tick.first) < roi_size_ch && abs(ticks.at(i_h) - vertex_ch_tick.second) < roi_size_tick){
      channels_tmp.push_back(channels.at(i_h));
      ticks_tmp.push_back(ticks.at(i_h));
      widths_tmp.push_back(widths.at(i_h));
    }
  }

  channels = channels_tmp;
  ticks = ticks_tmp;
  widths = widths_tmp;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::pair<double,double>> FindPeaksTH2D(const TH2D* h){

  std::vector<double> peak_r;
  std::vector<double> peak_theta;
  std::vector<double> peak_height;
  std::vector<std::pair<double,double>> peaks;

  const int peak_size = 3;

  // Identify any bins taller than all of their neighbours
  for(int i_bx=peak_size-1;i_bx<h->GetNbinsX()-peak_size+1;i_bx++){
    for(int i_by=peak_size-1;i_by<h->GetNbinsY()-peak_size+1;i_by++){

      if(h->GetBinContent(i_bx,i_by) < 3) continue;

      double content = h->GetBinContent(i_bx,i_by);
      bool is_peak = true;
      for(int i=-peak_size;i<=peak_size;i++){
        for(int j=-peak_size;j<=peak_size;j++){
        if(i == 0 && j == 0) continue;
        if(content < h->GetBinContent(i_bx+i,i_by+j)) is_peak = false;
        }
      }

      if(is_peak){
        double height = h->GetBinContent(i_bx,i_by);
        if(!peak_height.size()){
          peak_height.push_back(height);
          peak_r.push_back(h->GetXaxis()->GetBinCenter(i_bx));
          peak_theta.push_back(h->GetYaxis()->GetBinCenter(i_by));
        }

        for(size_t i=0;i<peak_height.size();i++){
           if(height > peak_height.at(i)){
            peak_height.insert(peak_height.begin() + i, height);
            peak_r.insert(peak_r.begin() + i, h->GetXaxis()->GetBinCenter(i_bx));
            peak_theta.insert(peak_theta.begin() + i, h->GetYaxis()->GetBinCenter(i_by));
            break;
          }
        }
      }

    }
  }


  for(size_t i=0;i<peak_height.size();i++)
    peaks.push_back(std::make_pair(peak_r.at(i),peak_theta.at(i)));

  std::cout << "Found " << peaks.size() << " peaks in Hough transform" << std::endl;

  return peaks;

}

}

