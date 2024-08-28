#ifndef _VFitter_cxx_
#define _VFitter_cxx_

#include "VFitter.h"

namespace hyperonreco {

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

VFitter::VFitter(bool draw){

  Hits.resize(3);
  Draw = draw;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VFitter::SetEvent(int run,int subrun,int event,int combination){

  Run = run;
  Subrun = subrun;
  Event = event;
  Combination = combination;
  RSE = std::to_string(run) + "_" + std::to_string(subrun) + "_" + std::to_string(event);
  RSEC = std::to_string(run) + "_" + std::to_string(subrun) + "_" + std::to_string(event) + "_" + std::to_string(Combination);

  system(("mkdir -p Plots/Event_" + RSE).c_str());

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VFitter::SetROI(double roi_ch,double roi_tick,TVector3 center){

  ROIChannel = roi_ch;
  ROITick = roi_tick;
  ROICenter = center;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VFitter::RemoveOffset(std::vector<HitLite>& hits,TVector3 origin,int plane){

  std::pair<double,double> origin_ch_ti = WireTick(origin,plane);
  //std::cout << "Removing offset " << origin_ch_ti.first << "  " << origin_ch_ti.second << std::endl;

  for(HitLite& hit : hits){
    hit.Channel = (hit.Channel - origin_ch_ti.first)/A_w;
    hit.Tick = (hit.Tick - origin_ch_ti.second)/A_t;
    hit.Width /= A_t;
  }
/*
  if(Draw){

    std::vector<double> channels,ticks,widths;
    for(HitLite hit : hits){
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

    c->Print(("Plots/Event_" + RSE + "_Transposed_data_Plane" + std::to_string(plane) + ".png").c_str());
    c->Close();

  }
*/
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VFitter::RestoreOffset(std::vector<HitLite>& hits,TVector3 origin,int plane){

  std::pair<double,double> origin_ch_ti = WireTick(origin,plane);

  for(HitLite& hit : hits){
    hit.Channel = hit.Channel*A_w + origin_ch_ti.first;
    hit.Tick = hit.Tick*A_t + origin_ch_ti.second;
    hit.Width *= A_t;
  }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VFitter::SetFitTune(double fittune){

  FitTune = fittune;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VFitter::SetActivePlanes(std::vector<int> i_pl){

  ActivePlanes = i_pl;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VFitter::Reset(){

  for(size_t i_pl=0;i_pl<3;i_pl++){
      Hits.at(i_pl).clear();
  }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VFitter::AddData(const std::vector<std::vector<HitLite>>& hits){

  for(size_t i_pl=0;i_pl<3;i_pl++){
    if(std::find(ActivePlanes.begin(),ActivePlanes.end(),i_pl) == ActivePlanes.end()) continue;
    Hits.at(i_pl).insert(Hits.at(i_pl).end(),hits.at(i_pl).begin(),hits.at(i_pl).end());
  }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VFitter::AddData(const HoughTransformPoint& p){

  if(std::find(ActivePlanes.begin(),ActivePlanes.end(),p.Plane) == ActivePlanes.end()) return;

  Hits.at(p.Plane).insert(Hits.at(p.Plane).end(),p.Hits.begin(),p.Hits.end());

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VFitter::SetOutlierCut(double outliercut){

  OutlierCut = outliercut;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VFitter::SetGuess(const FittedV& fittedv){

  InitialGuess = fittedv;
  InitialGuessSet = true;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fit metric between 3D V and collection of hits 

std::pair<double,int> VFitter::FitScore2(const std::vector<std::vector<HitLite>>& hits,FittedV& fittedv,bool verbose){

  double score = 0.0;
  int ndof = 0;

  TVector3 arm1_direction_3d = fittedv.GetArm1Dir();
  TVector3 arm2_direction_3d = fittedv.GetArm2Dir();

  fittedv.Arm1Points = 0;
  fittedv.Arm2Points = 0;

  for(int i_pl=0;i_pl<kPlaneInvalid;i_pl++){

    //if(verbose) std::cout << "Plane " << i_pl << std::endl;

    LineWireTick2 line_1;
    LineWireTick2 line_2;
    if(i_pl == 0){
      line_1.Start = TVector2(-sin60*fittedv.Vertex.Y()+cos60*fittedv.Vertex.Z(),fittedv.Vertex.X());
      line_1.Direction = TVector2(-sin60*arm1_direction_3d.Y()+cos60*arm1_direction_3d.Z(),arm1_direction_3d.X());
      line_2.Start = TVector2(-sin60*fittedv.Vertex.Y()+cos60*fittedv.Vertex.Z(),fittedv.Vertex.X());
      line_2.Direction = TVector2(-sin60*arm2_direction_3d.Y()+cos60*arm2_direction_3d.Z(),arm2_direction_3d.X());
    }
    if(i_pl == 1){
      line_1.Start = TVector2(sin60*fittedv.Vertex.Y()+cos60*fittedv.Vertex.Z(),fittedv.Vertex.X());
      line_1.Direction = TVector2(sin60*arm1_direction_3d.Y()+cos60*arm1_direction_3d.Z(),arm1_direction_3d.X());
      line_2.Start = TVector2(sin60*fittedv.Vertex.Y()+cos60*fittedv.Vertex.Z(),fittedv.Vertex.X());
      line_2.Direction = TVector2(sin60*arm2_direction_3d.Y()+cos60*arm2_direction_3d.Z(),arm2_direction_3d.X());
    }
    if(i_pl == 2){
      line_1.Start = TVector2(fittedv.Vertex.Z(),fittedv.Vertex.X());
      line_1.Direction = TVector2(arm1_direction_3d.Z(),arm1_direction_3d.X());
      line_2.Start = TVector2(fittedv.Vertex.Z(),fittedv.Vertex.X());
      line_2.Direction = TVector2(arm2_direction_3d.Z(),arm2_direction_3d.X());
    }
    double mag = line_1.Direction.Mod();
    line_1.Direction /= mag;
    mag = line_2.Direction.Mod();
    line_2.Direction /= mag;
/*
    if(verbose && Draw){
      TCanvas* c = new TCanvas("c","c");
      std::vector<double> channels,ticks,widths;
      for(HitLite hit : Hits.at(i_pl)){
        channels.push_back(hit.Channel); 
        ticks.push_back(hit.Tick); 
        widths.push_back(hit.Width);
      }
      TGraphErrors* g = new TGraphErrors(channels.size(),&(channels[0]),&(ticks[0]),0,&(widths[0]));
      g->SetName("graph2D");
      g->SetMarkerColor(1);
      g->SetMarkerSize(0.5);
      g->Draw("AP");
      TGraph* g1 = new TGraph(2);
      TGraph* g2 = new TGraph(2);
      g1->SetName("fit_graph2D_1");
      g2->SetName("fit_graph2D_2");
      g1->SetPoint(0,line_1.Start.X(),line_1.Start.Y());
      g2->SetPoint(0,line_2.Start.X(),line_2.Start.Y());
      g1->SetPoint(1,line_1.Start.X()+line_1.Direction.X()*200,line_1.Start.Y()+line_1.Direction.Y()*200);
      g2->SetPoint(1,line_2.Start.X()+line_2.Direction.X()*200,line_2.Start.Y()+line_2.Direction.Y()*200);
      g1->SetLineColor(2); 
      g2->SetLineColor(3); 
      g1->Draw("L same");
      g2->Draw("L same");
      c->Print(("Plots/Event_" + RSE + "_Fit_Transposed_Data_Plane" + std::to_string(i_pl) + ".png").c_str());
    }
   */ 

    for(HitLite hit : hits.at(i_pl)){
      double arm1_sep = HitLineSeparation2(hit,line_1);
      double arm2_sep = HitLineSeparation2(hit,line_2);
      if(arm1_sep < arm2_sep){
        fittedv.Arm1Points++;
        score += arm1_sep; 
      }
      else {
        fittedv.Arm2Points++;
        score += arm2_sep; 
      }
      ndof++;
    }
  }

  fittedv.Chi2 = score;
  fittedv.NDof = ndof;

  return std::make_pair(score,ndof);

} 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Distance between a hit and a 2D line divided by width of hit 

double VFitter::HitLineSeparation2(const HitLite& hit,const LineWireTick2& line){

  double dot = (hit.Channel-line.Start.X())*line.Direction.X() + (hit.Tick-line.Start.Y())*line.Direction.Y();
  if(dot < 0) return (hit.Channel-line.Start.X())*(hit.Channel-line.Start.X()) + (hit.Tick-line.Start.Y())*(hit.Tick-line.Start.Y()); 
  else return (line.Start.X()-hit.Channel + dot*line.Direction.X())*(line.Start.X()-hit.Channel + dot*line.Direction.X()) + (line.Start.Y()-hit.Tick + dot*line.Direction.Y())*(line.Start.Y()-hit.Tick + dot*line.Direction.Y());

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Use minimizer to do fitting 

bool VFitter::DoFit(FittedV& fittedv){

  //std::cout << "Starting grid search 3" << std::endl;

  for(size_t i_pl=0;i_pl<3;i_pl++) RemoveOffset(Hits.at(i_pl),InitialGuess.Vertex,i_pl);

  FittedV bestfit_v;

  int points_suc = 0;

  ROOT::Math::Functor min = ROOT::Math::Functor( [&] (const double *coeff ){

      //std::cout << "Evaluating fit function" << std::endl;

      FittedV v;
      v.Vertex = fittedv.Vertex;
      v.Arm1Theta = coeff[0];
      v.Arm1Phi = coeff[1];
      v.Arm2Theta = coeff[2];
      v.Arm2Phi = coeff[3];
      v.Vertex = TVector3(coeff[4],coeff[5],coeff[6]);

      /*
         std::cout << "Vertex: " << v.Vertex.X() << "  " << v.Vertex.Y() << "  " << v.Vertex.Z() << std::endl;
         std::cout << "Arm 1: " << v.GetArm1Dir().X() << "  " << v.GetArm1Dir().Y() << "  " << v.GetArm1Dir().Z() << std::endl;
         std::cout << "Arm 2: " << v.GetArm2Dir().X() << "  " << v.GetArm2Dir().Y() << "  " << v.GetArm2Dir().Z() << std::endl;
         */

      std::pair<double,int> FitVal = FitScore2(Hits,v);

      return FitVal.first/1e8;

  } , 7);

  std::unique_ptr< ROOT::Math::Minimizer > fMinimizer = std::unique_ptr<ROOT::Math::Minimizer>
    ( ROOT::Math::Factory::CreateMinimizer( "Minuit2", "Migrad" ) );

  fMinimizer->SetMaxFunctionCalls(1000);
  fMinimizer->SetTolerance( 0.001 );
  fMinimizer->SetVariable(0,"Arm 1 Theta",InitialGuess.Arm1Theta,0.01);
  fMinimizer->SetVariable(1,"Arm 1 Phi",InitialGuess.Arm1Phi,0.01);
  fMinimizer->SetVariable(2,"Arm 2 Theta",InitialGuess.Arm2Theta,0.01);
  fMinimizer->SetVariable(3,"Arm 2 Phi",InitialGuess.Arm2Phi,0.01);
  fMinimizer->SetVariable(4,"Vertex X",0.0,0.05);
  fMinimizer->SetVariable(5,"Vertex Y",0.0,0.05);
  fMinimizer->SetVariable(6,"Vertex Z",0.0,0.05);

  fMinimizer->SetFunction(min);
  fMinimizer->Minimize();

  //if(fMinimizer->Status() != 0) std::cout << "Minimization did not converge" << std::endl;
  //else std::cout << "Minimization succeeded" << std::endl;

  points_suc++;

  bestfit_v.Arm1Theta = fMinimizer->X()[0];
  bestfit_v.Arm1Phi = fMinimizer->X()[1];
  bestfit_v.Arm2Theta = fMinimizer->X()[2];
  bestfit_v.Arm2Phi = fMinimizer->X()[3];
  bestfit_v.Vertex = TVector3(fMinimizer->X()[4],fMinimizer->X()[5],fMinimizer->X()[6]);

  std::pair<double,int> FitVal = FitScore2(Hits,fittedv);

  //std::cout << points_suc << "/" << points << std::endl;

  FitScore2(Hits,bestfit_v,true);

  fittedv = bestfit_v;
  fittedv.Vertex = bestfit_v.Vertex + InitialGuess.Vertex;
/*
  std::cout << "Best Fit:" << std::endl;
  std::cout << "Vertex: " << bestfit_v.Vertex.X() << "  " << bestfit_v.Vertex.Y() <<"  " << bestfit_v.Vertex.Z() << std::endl;
  std::cout << "Arm 1: theta=" << bestfit_v.Arm1Theta << "  phi=" << bestfit_v.Arm1Phi << std::endl;
  std::cout << "Arm 2: theta=" << bestfit_v.Arm2Theta << "  phi=" << bestfit_v.Arm2Phi << std::endl;
  std::cout << "Score=" << bestfit_v.Chi2 << "/"  << bestfit_v.NDof << "=" << bestfit_v.Chi2/bestfit_v.NDof << std::endl;
*/

  for(size_t i_pl=0;i_pl<3;i_pl++) RestoreOffset(Hits.at(i_pl),InitialGuess.Vertex,i_pl);

  //if(Draw) DrawFit(fittedv);

  return fMinimizer->Status() == 0;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hybrid - grid seach for vertex but use minimizer to find angles

bool VFitter::DoFitGridSearch3(FittedV& fittedv,int points){

  //std::cout << "Starting grid search 3" << std::endl;

  for(size_t i_pl=0;i_pl<3;i_pl++) RemoveOffset(Hits.at(i_pl),InitialGuess.Vertex,i_pl);
  TRandom2* r = new TRandom2();

  std::pair<double,double> vertex_range = {-3.0,3.0};

  double bestfit_score=1e10;
  FittedV bestfit_v;

  int points_suc = 0;

  for(int i=0;i<points;i++){
    if(i % 1000 == 0) std::cout << "Grid search " << (double)i/points*100 << "\% complete" << std::endl;
    fittedv.Vertex = TVector3(r->Uniform(vertex_range.first,vertex_range.second),r->Uniform(vertex_range.first,vertex_range.second),r->Uniform(vertex_range.first,vertex_range.second));

    ROOT::Math::Functor min = ROOT::Math::Functor( [&] (const double *coeff ){

        //std::cout << "Evaluating fit function" << std::endl;

        FittedV v;
        v.Vertex = fittedv.Vertex;
        v.Arm1Theta = coeff[0];
        v.Arm1Phi = coeff[1];
        v.Arm2Theta = coeff[2];
        v.Arm2Phi = coeff[3];

        /*
           std::cout << "Vertex: " << v.Vertex.X() << "  " << v.Vertex.Y() << "  " << v.Vertex.Z() << std::endl;
           std::cout << "Arm 1: " << v.GetArm1Dir().X() << "  " << v.GetArm1Dir().Y() << "  " << v.GetArm1Dir().Z() << std::endl;
           std::cout << "Arm 2: " << v.GetArm2Dir().X() << "  " << v.GetArm2Dir().Y() << "  " << v.GetArm2Dir().Z() << std::endl;
           */

        std::pair<double,int> FitVal = FitScore2(Hits,v);

        return FitVal.first/1e6;

    } , 4);

    std::unique_ptr< ROOT::Math::Minimizer > fMinimizer = std::unique_ptr<ROOT::Math::Minimizer>
      ( ROOT::Math::Factory::CreateMinimizer( "Minuit2", "Migrad" ) );

    fMinimizer->SetMaxFunctionCalls(1000);
    fMinimizer->SetTolerance( 0.001 );
    fMinimizer->SetVariable(0,"Arm 1 Theta",InitialGuess.Arm1Theta,0.01);
    fMinimizer->SetVariable(1,"Arm 1 Phi",InitialGuess.Arm1Phi,0.01);
    fMinimizer->SetVariable(2,"Arm 2 Theta",InitialGuess.Arm2Theta,0.01);
    fMinimizer->SetVariable(3,"Arm 2 Phi",InitialGuess.Arm2Phi,0.01);
    fMinimizer->SetFunction(min);
    fMinimizer->Minimize();

    if(fMinimizer->Status() != 0) continue;

    points_suc++;

    fittedv.Arm1Theta = fMinimizer->X()[0];
    fittedv.Arm1Phi = fMinimizer->X()[1];
    fittedv.Arm2Theta = fMinimizer->X()[2];
    fittedv.Arm2Phi = fMinimizer->X()[3];

    std::pair<double,int> FitVal = FitScore2(Hits,fittedv);

    if(FitVal.first < bestfit_score){
      bestfit_score = FitVal.first; 
      bestfit_v = fittedv;
      /*
      std::cout << "New bestfit" << std::endl;
      std::cout << "Vertex: " << bestfit_v.Vertex.X() << "  " << bestfit_v.Vertex.Y() <<"  " << bestfit_v.Vertex.Z() << std::endl;
      std::cout << "Arm 1: theta=" << bestfit_v.Arm1Theta << "  phi=" << bestfit_v.Arm1Phi << std::endl;
      std::cout << "Arm 2: theta=" << bestfit_v.Arm2Theta << "  phi=" << bestfit_v.Arm2Phi << std::endl;
      std::cout << "Score=" << fittedv.Chi2 << "/"  << fittedv.NDof << "=" << fittedv.Chi2/fittedv.NDof << std::endl;
      */
    }

  }

  //std::cout << points_suc << "/" << points << std::endl;

  FitScore2(Hits,bestfit_v,true);

  fittedv = bestfit_v;
  fittedv.Vertex = bestfit_v.Vertex + InitialGuess.Vertex;
/*
  std::cout << "Best Fit:" << std::endl;
  std::cout << "Vertex: " << bestfit_v.Vertex.X() << "  " << bestfit_v.Vertex.Y() <<"  " << bestfit_v.Vertex.Z() << std::endl;
  std::cout << "Arm 1: theta=" << bestfit_v.Arm1Theta << "  phi=" << bestfit_v.Arm1Phi << std::endl;
  std::cout << "Arm 2: theta=" << bestfit_v.Arm2Theta << "  phi=" << bestfit_v.Arm2Phi << std::endl;
  std::cout << "Score=" << bestfit_v.Chi2 << "/"  << bestfit_v.NDof << "=" << bestfit_v.Chi2/bestfit_v.NDof << std::endl;
*/
  delete r;

  for(size_t i_pl=0;i_pl<3;i_pl++) RestoreOffset(Hits.at(i_pl),InitialGuess.Vertex,i_pl);

  //if(Draw) DrawFit(fittedv);

  return true;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VFitter::DrawFit(const FittedV& v,const std::vector<std::vector<HitLite>>& allhits) const {

  TCanvas* c = new TCanvas("c","c",800,600);

  TPad *p_plot = new TPad("pad1","pad1",0,0,1,0.85);
  TPad *p_legend = new TPad("pad2","pad2",0,0.85,1,1);
  p_legend->SetBottomMargin(0);
  p_legend->SetTopMargin(0.1);
  p_plot->SetTopMargin(0.01);

  TLegend *l = new TLegend(0.1,0.0,0.9,1.0);
  l->SetBorderSize(0);
  l->SetNColumns(2);

  p_legend->Draw();
  p_legend->cd();
  l->Draw();
  c->cd();
  p_plot->Draw();
  p_plot->cd();

  bool has_allhits = allhits.size() == 3;

  for(int i_pl=0;i_pl<3;i_pl++){

    TMultiGraph* gm = new TMultiGraph();
    TGraphErrors* g_allhits = nullptr;

    std::vector<double> channel,tick,width;

    if(has_allhits){

      channel.clear(); 
      tick.clear();      
      width.clear();
        
      for(HitLite hit : allhits.at(i_pl)){
        channel.push_back(hit.Channel);
        tick.push_back(hit.Tick);
        width.push_back(hit.Width);
      }

      g_allhits = new TGraphErrors(channel.size(),&(channel[0]),&(tick[0]),0,&(width[0]));
      g_allhits->SetName("graph2D");
      g_allhits->SetMarkerColor(1);
      g_allhits->SetLineColor(1);
      g_allhits->SetMarkerStyle(20);
      g_allhits->SetMarkerSize(0.4);
      gm->Add(g_allhits);

    }

    channel.clear(); 
    tick.clear();      
    width.clear();

    for(HitLite hit : Hits.at(i_pl)){
      channel.push_back(hit.Channel);
      tick.push_back(hit.Tick);
      width.push_back(hit.Width);
    }

    TGraphErrors* g = new TGraphErrors(channel.size(),&(channel[0]),&(tick[0]),0,&(width[0]));
    g->SetName("graph2D");
    g->SetMarkerColor(2);
    g->SetLineColor(2);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.6);

    gm->Add(g);

    gm->Draw("AP");

    std::pair<TGraph*,TGraph*> fit_arms = v.Make2DVGraphs(i_pl); 

    fit_arms.first->Draw("L same");
    fit_arms.second->Draw("L same");
    l->Clear(); 
    //l->AddEntry((TObject*)0,("Chi2/ndof=" + std::to_string(v.Chi2) + " / " + std::to_string(v.NDof) + "^{2} = " + std::to_string(v.Chi2/v.NDof/v.NDof)).c_str(),"");
    l->AddEntry((TObject*)0,("Chi2=" + std::to_string(v.Chi2)).c_str(),"");
    l->AddEntry((TObject*)0,("NDof=" + std::to_string(v.NDof)).c_str(),"");
    l->AddEntry((TObject*)0,("Asymmetry=" + std::to_string(v.GetAsymmetry())).c_str(),"");
    l->AddEntry((TObject*)0,("Opening Angle=" + std::to_string(v.GetOpeningAngle())).c_str(),"");
    c->Print(("Plots/Event_" + RSE + "/Event_" + RSE + "_Plane" + std::to_string(i_pl) + "_" + std::to_string(Combination) + "_Fit.png").c_str());
    p_plot->Clear();

    delete g;
    //delete gm;

  }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VFitter::DrawFit2(const FittedV& v,const std::vector<std::vector<HitLite>>& allhits) const {

  bool has_allhits = allhits.size() == 3;

  TCanvas* c = new TCanvas("c","c",1400,600);
  TLegend *l = new TLegend(0.1,0.0,0.9,1.0);
  l->SetBorderSize(0);
  l->SetNColumns(2);

  TPad *p_plot_Plane0 = new TPad("pad_Plane0","pad_Plane0",0,0,0.333,0.85);
  TPad *p_plot_Plane1 = new TPad("pad_Plane1","pad_Plane1",0.333,0,0.666,0.85);
  TPad *p_plot_Plane2 = new TPad("pad_Plane2","pad_Plane2",0.666,0,1.0,0.85);

  TPad *p_legend = new TPad("pad2","pad2",0,0.85,1,1);
  p_legend->SetBottomMargin(0);
  p_legend->SetTopMargin(0.1);
  p_plot_Plane0->SetTopMargin(0.01);
  p_plot_Plane1->SetTopMargin(0.01);
  p_plot_Plane2->SetTopMargin(0.01);

  p_legend->Draw();
  p_legend->cd();
  l->Draw();
  c->cd();
  p_plot_Plane0->Draw();
  p_plot_Plane1->Draw();
  p_plot_Plane2->Draw();

  std::vector<TMultiGraph*> gm;  
  std::vector<std::pair<TGraph*,TGraph*>> fit_arms; 

  for(int i_pl=0;i_pl<3;i_pl++){

    if(i_pl == 0) p_plot_Plane0->cd();
    if(i_pl == 1) p_plot_Plane1->cd();
    if(i_pl == 2) p_plot_Plane2->cd();

    gm.push_back(new TMultiGraph(("Plane"+std::to_string(i_pl)).c_str(),";;"));
    
    TGraphErrors* g_allhits = nullptr;

    std::vector<double> channel,tick,width;

    if(has_allhits){

      channel.clear(); 
      tick.clear();      
      width.clear();
        
      for(HitLite hit : allhits.at(i_pl)){
        channel.push_back(hit.Channel);
        tick.push_back(hit.Tick);
        width.push_back(hit.Width);
      }

      g_allhits = new TGraphErrors(channel.size(),&(channel[0]),&(tick[0]),0,&(width[0]));
      g_allhits->SetName("graph2D");
      g_allhits->SetMarkerColor(1);
      g_allhits->SetLineColor(1);
      g_allhits->SetMarkerStyle(20);
      g_allhits->SetMarkerSize(0.4);
      gm.back()->Add(g_allhits);

    }

    channel.clear(); 
    tick.clear();      
    width.clear();

    for(HitLite hit : Hits.at(i_pl)){
      channel.push_back(hit.Channel);
      tick.push_back(hit.Tick);
      width.push_back(hit.Width);
    }

    TGraphErrors* g = new TGraphErrors(channel.size(),&(channel[0]),&(tick[0]),0,&(width[0]));
    g->SetName("graph2D");
    g->SetMarkerColor(2);
    g->SetLineColor(2);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.6);

    gm.back()->Add(g);
    gm.back()->Draw("AP");

    fit_arms.push_back(v.Make2DVGraphs(i_pl)); 

    fit_arms.back().first->Draw("L same");
    fit_arms.back().second->Draw("L same");

    c->cd();

  }
  
  //l->AddEntry((TObject*)0,("Chi2/ndof=" + std::to_string(v.Chi2) + " / " + std::to_string(v.NDof) + "^{2} = " + std::to_string(v.Chi2/v.NDof/v.NDof)).c_str(),"");
  l->AddEntry((TObject*)0,("Chi2=" + std::to_string(v.Chi2)).c_str(),"");
  l->AddEntry((TObject*)0,("NDof=" + std::to_string(v.NDof)).c_str(),"");
  l->AddEntry((TObject*)0,("Score=" + std::to_string(v.Chi2/pow(v.NDof,4))).c_str(),"");
  l->AddEntry((TObject*)0,("Asymmetry=" + std::to_string(v.GetAsymmetry())).c_str(),"");
  l->AddEntry((TObject*)0,("Opening Angle=" + std::to_string(v.GetOpeningAngle())).c_str(),"");

  c->Print(("Plots/Event_" + RSE + "/Event_" + RSE + "_" + std::to_string(Combination) + "_Fit.png").c_str());

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

}

#endif
