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

void VFitter::SetEvent(int run,int subrun,int event){

  Run = run;
  Subrun = subrun;
  Event = event;
  RSE = std::to_string(run) + "_" + std::to_string(subrun) + "_" + std::to_string(event);

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
  std::cout << "Removing offset " << origin_ch_ti.first << "  " << origin_ch_ti.second << std::endl;

  for(HitLite& hit : hits){
    hit.Channel = (hit.Channel - origin_ch_ti.first)/A_w;
    hit.Tick = (hit.Tick - origin_ch_ti.second)/A_t;
    hit.Width /= A_t;
  }

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
/*
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VFitter::RemoveOutliers(){

  std::vector<std::vector<double>> channel_v_tmp(3),tick_v_tmp(3),width_v_tmp(3);

  // Calculate the SD of the fit metrics
  double rms = 0.0;
  int points = 0;

  for(int i_pl=0;i_pl<kPlaneInvalid;i_pl++){

    LineWireTick2 line_1 = InitialGuess.GetArm1_2D2(i_pl); 
    LineWireTick2 line_2 = InitialGuess.GetArm2_2D2(i_pl); 
    std::pair<double,double> vertex_ch_tick = WireTick(ROICenter,i_pl); 

    for(size_t i_h=0;i_h<Channel_v.at(i_pl).size();i_h++){
      double channel = Channel_v.at(i_pl).at(i_h);
      double tick = Tick_v.at(i_pl).at(i_h);
      double width = Width_v.at(i_pl).at(i_h);
      if(abs(channel-vertex_ch_tick.first) > ROIChannel || abs(tick-vertex_ch_tick.second) > ROITick) continue;
      double sep = std::min(HitLineSeparation(channel,tick,width,line_1),HitLineSeparation(channel,tick,width,line_2));
      rms += sep*sep;
      points++;
    }

  }

  rms = sqrt(rms/points);  
  std::cout << "rms of displacements = " << rms << std::endl;

  for(int i_pl=0;i_pl<kPlaneInvalid;i_pl++){

    LineWireTick2 line_1 = InitialGuess.GetArm1_2D2(i_pl); 
    LineWireTick2 line_2 = InitialGuess.GetArm2_2D2(i_pl); 
    std::pair<double,double> vertex_ch_tick = WireTick(ROICenter,i_pl); 

    for(size_t i_h=0;i_h<Channel_v.at(i_pl).size();i_h++){
      double channel = Channel_v.at(i_pl).at(i_h);
      double tick = Tick_v.at(i_pl).at(i_h);
      double width = Width_v.at(i_pl).at(i_h);
      double sep = std::min(HitLineSeparation(channel,tick,width,line_1),HitLineSeparation(channel,tick,width,line_2));
      if(abs(channel-vertex_ch_tick.first) > ROIChannel || abs(tick-vertex_ch_tick.second) > ROITick || sep/rms < OutlierCut){
        channel_v_tmp.at(i_pl).push_back(channel);
        tick_v_tmp.at(i_pl).push_back(tick);
        width_v_tmp.at(i_pl).push_back(width);
      }
    }

  }

  Channel_v = channel_v_tmp;
  Tick_v = tick_v_tmp;
  Width_v = width_v_tmp;
}
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VFitter::SetOutlierCut(double outliercut){

  OutlierCut = outliercut;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _InsideLArSoft_
FittedV VFitter::DoFit(TVector3 start,std::map<art::Ptr<recob::SpacePoint>,art::Ptr<recob::Hit>> hitspacepointhap) const {

  FittedV result;

  std::cout << "Fitting V to shower with starting position " << start.X() << " " << start.Y() << " " << start.Z() << std::endl;

  return result;

}
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Distance between a hit and a 2D line divided by width of hit 

double VFitter::HitLineSeparation(const HitLite& hit,const LineWireTick2& line){

  TVector2 p(hit.Channel,hit.Tick);

  const TVector2& a = line.Start;
  const TVector2& n = line.Direction;

  double dot = (p-a).X()*n.X() + (p-a).Y()*n.Y();

  if(dot < 0){
    return (a-p).X()*(a-p).X() + (a-p).Y()*(a-p).Y()/hit.Width/hit.Width/TickPerWire/TickPerWire; 
  }

  TVector2 nearest_approach = (a-p)+dot*n;

  return nearest_approach.X()*nearest_approach.X() + nearest_approach.Y()*nearest_approach.Y()/hit.Width/hit.Width/TickPerWire/TickPerWire;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Distance between a hit and a 2D line divided by width of hit 

#ifdef _InsideLArSoft_
double VFitter::HitLineSeparation(art::Ptr<recob::Hit> hit,LineWireTick line){

  int hit_channel = hit->Channel();
  float hit_peak = hit->PeakTime();
  float hit_width = hit->EndTick() - hit->StartTick();

  return HitLineSeparation(hit_channel,hit_peak,hit_width,line);

}
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fit metric between 3D V and collection of hits 

std::pair<double,int> VFitter::FitScore(const std::vector<std::vector<HitLite>>& hits,FittedV& fittedv,bool verbose){

  double score = 0.0;
  int ndof = 0;
  fittedv.Arm1Points=0;
  fittedv.Arm2Points=0;

  for(int i_pl=0;i_pl<kPlaneInvalid;i_pl++){

    LineWireTick2 line_1 = fittedv.GetArm1_2D2(i_pl); 
    LineWireTick2 line_2 = fittedv.GetArm2_2D2(i_pl); 
    std::pair<double,double> vertex_ch_tick = WireTick(InitialGuess.Vertex,i_pl);
    
    //LineWireTick2 line_1;

    for(HitLite hit : hits.at(i_pl)){

      if(abs(hit.Channel-vertex_ch_tick.first) > ROIChannel || abs(hit.Tick-vertex_ch_tick.second) > ROITick) continue;
      double arm1_sep = HitLineSeparation(hit,line_1);
      double arm2_sep = HitLineSeparation(hit,line_2);

      if(arm1_sep < arm2_sep){
        score += arm1_sep;
        fittedv.Arm1Points++;
      }
      else {
        score += arm2_sep;
        fittedv.Arm2Points++;
      }

      ndof++;
    }
  }

  fittedv.Chi2 = score;
  fittedv.NDof = ndof;

  return std::make_pair(score,ndof);

} 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fit metric between 3D V and collection of hits 

#ifdef _InsideLArSoft_
std::pair<double,int> VFitter::FitScore(std::vector<art::Ptr<recob::Hit>> hit_v,FittedV fittedv){

  double score = 0.0;

  for(int i_pl=0;i_pl<kInvalid;i_pl++){
    fittedv.Arm1_2D.at(i_pl) = ProjectXYZWireTick(fittedv.Vertex,fittedv.Arm1Dir,i_pl);
    fittedv.Arm2_2D.at(i_pl) = ProjectXYZWireTick(fittedv.Vertex,fittedv.Arm2Dir,i_pl);

  }

  int ndof = 0;
  for(art::Ptr<recob::Hit> hit : hit_v){
    double arm1_sep = HitLineSeparation(hit,fittedv.Arm1_2D.at(hit->View()));
    double arm2_sep = HitLineSeparation(hit,fittedv.Arm2_2D.at(hit->View()));
    score += std::min(arm1_sep,arm2_sep);
    ndof++;
  }

  return std::make_pair(score,ndof);

} 
#endif


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Do Fit, returns bool indicating if fit was successful 

bool VFitter::DoFit(FittedV& fittedv){

  //std::cout << "Setting up fit" << std::endl;

  ROOT::Math::Functor min = ROOT::Math::Functor( [&] (const double *coeff ){

      //std::cout << "Evaluating fit function" << std::endl;

      FittedV v;
      v.Vertex = TVector3(coeff[0],coeff[1],coeff[2]);
      v.Arm1Theta = coeff[3];
      v.Arm1Phi = coeff[4];
      v.Arm2Theta = coeff[5];
      v.Arm2Phi = coeff[6];

      /*
         std::cout << "Vertex: " << v.Vertex.X() << "  " << v.Vertex.Y() << "  " << v.Vertex.Z() << std::endl;
         std::cout << "Arm 1: " << v.GetArm1Dir().X() << "  " << v.GetArm1Dir().Y() << "  " << v.GetArm1Dir().Z() << std::endl;
         std::cout << "Arm 2: " << v.GetArm2Dir().X() << "  " << v.GetArm2Dir().Y() << "  " << v.GetArm2Dir().Z() << std::endl;
         */

      std::pair<double,int> FitVal = FitScore(Hits,v);

      return FitVal.first/FitVal.second;

  } , 7);

  std::unique_ptr< ROOT::Math::Minimizer > fMinimizer = std::unique_ptr<ROOT::Math::Minimizer>
    ( ROOT::Math::Factory::CreateMinimizer( "Minuit2", "Migrad" ) );

  fMinimizer->SetMaxFunctionCalls(10000);
  fMinimizer->SetTolerance( 0.0001 );

  fMinimizer->SetVariable(0,"Vertex X",InitialGuess.Vertex.X(),0.3);
  fMinimizer->SetVariable(1,"Vertex Y",InitialGuess.Vertex.Y(),0.3);
  fMinimizer->SetVariable(2,"Vertex Z",InitialGuess.Vertex.Z(),0.3);
  fMinimizer->SetVariable(3,"Arm 1 Theta",InitialGuess.Arm1Theta,0.005);
  fMinimizer->SetVariable(4,"Arm 1 Phi",InitialGuess.Arm1Phi,0.005);
  fMinimizer->SetVariable(5,"Arm 2 Theta",InitialGuess.Arm2Theta,0.005);
  fMinimizer->SetVariable(6,"Arm 2 Phi",InitialGuess.Arm2Phi,0.005);

  fMinimizer->SetFunction(min);
  fMinimizer->Minimize();

  fittedv.Vertex = TVector3(fMinimizer->X()[0],fMinimizer->X()[1],fMinimizer->X()[2]);
  fittedv.Arm1Theta = fMinimizer->X()[3];
  fittedv.Arm1Phi = fMinimizer->X()[4];
  fittedv.Arm2Theta = fMinimizer->X()[5];
  fittedv.Arm2Phi = fMinimizer->X()[6];

  std::pair<double,int> FitVal = FitScore(Hits,fittedv);

  std::cout << "Score = " << FitVal.first << "/" << FitVal.second << "=" << FitVal.first/FitVal.second << std::endl;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool VFitter::DoFitGridSearch(FittedV& fittedv,int points){

  TRandom2* r = new TRandom2();

  std::pair<double,double> vertex_range = {-1.0,1.0};
  std::pair<double,double> arm1_theta_range = {-3.141,3.141};
  std::pair<double,double> arm2_theta_range = {-3.141,3.141};
  std::pair<double,double> arm1_phi_range = {-3.141,3.141};
  std::pair<double,double> arm2_phi_range = {-3.141,3.141};

  FittedV bestfit_v;
  double bestfit_score=1e10;

  for(int i=0;i<points;i++){
    if(i % 10000 == 0) std::cout << "Grid search iteration " << i << "/" << points << std::endl;
    fittedv.Vertex = InitialGuess.Vertex + TVector3(r->Uniform(vertex_range.first,vertex_range.second),r->Uniform(vertex_range.first,vertex_range.second),r->Uniform(vertex_range.first,vertex_range.second));
    fittedv.Arm1Theta = r->Uniform(arm1_theta_range.first,arm1_theta_range.second);
    fittedv.Arm2Theta = r->Uniform(arm2_theta_range.first,arm2_theta_range.second);
    fittedv.Arm1Phi = r->Uniform(arm1_phi_range.first,arm1_phi_range.second);
    fittedv.Arm2Phi = r->Uniform(arm2_phi_range.first,arm2_phi_range.second);
    std::pair<double,int> FitVal = FitScore(Hits,fittedv);
    if(FitVal.first < bestfit_score){
      bestfit_score = FitVal.first; 
      bestfit_v = fittedv;
    }
  }

  fittedv = bestfit_v;
  FitScore(Hits,fittedv,true); 

  delete r;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VFitter::SetGuess(const FittedV& fittedv){

  InitialGuess = fittedv;
  InitialGuessSet = true;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool VFitter::DoFitGridSearch2(FittedV& fittedv,int points){

  std::cout << "Starting grid search 2" << std::endl;
  for(size_t i_pl=0;i_pl<3;i_pl++) RemoveOffset(Hits.at(i_pl),InitialGuess.Vertex,i_pl);
  TRandom2* r = new TRandom2();

  std::pair<double,double> vertex_range = {-3.0,3.0};
  std::pair<double,double> arm1_theta_range = {-3.141,3.141};
  std::pair<double,double> arm2_theta_range = {-3.141,3.141};
  std::pair<double,double> arm1_phi_range = {-3.141,3.141};
  std::pair<double,double> arm2_phi_range = {-3.141,3.141};
 
  double bestfit_score=1e10;
  FittedV bestfit_v;

  for(int i=0;i<points;i++){
    if(i % 50000 == 0) std::cout << "Grid search " << (double)i/points*100 << "\% complete" << std::endl;
    fittedv.Vertex = TVector3(r->Uniform(vertex_range.first,vertex_range.second),r->Uniform(vertex_range.first,vertex_range.second),r->Uniform(vertex_range.first,vertex_range.second));
    fittedv.Arm1Theta = r->Uniform(arm1_theta_range.first,arm1_theta_range.second);
    fittedv.Arm2Theta = r->Uniform(arm2_theta_range.first,arm2_theta_range.second);
    fittedv.Arm1Phi = r->Uniform(arm1_phi_range.first,arm1_phi_range.second);
    fittedv.Arm2Phi = r->Uniform(arm2_phi_range.first,arm2_phi_range.second);

    std::pair<double,int> FitVal = FitScore2(Hits,fittedv);

    if(FitVal.first < bestfit_score){
      bestfit_score = FitVal.first; 
      bestfit_v = fittedv;
      std::cout << "New bestfit" << std::endl;
      std::cout << "Vertex: " << bestfit_v.Vertex.X() << "  " << bestfit_v.Vertex.Y() <<"  " << bestfit_v.Vertex.Z() << std::endl;
      std::cout << "Arm 1: theta=" << bestfit_v.Arm1Theta << "  phi=" << bestfit_v.Arm1Phi << std::endl;
      std::cout << "Arm 2: theta=" << bestfit_v.Arm2Theta << "  phi=" << bestfit_v.Arm2Phi << std::endl;
      std::cout << "Score=" << bestfit_score << std::endl;
    }
  }

  FitScore2(Hits,bestfit_v,true);

  fittedv = bestfit_v;
  fittedv.Vertex = bestfit_v.Vertex + InitialGuess.Vertex;
  //FitScore2(Channel_v,Tick_v,Width_v,fittedv,true); 

  delete r;

  for(size_t i_pl=0;i_pl<3;i_pl++) RestoreOffset(Hits.at(i_pl),InitialGuess.Vertex,i_pl);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fit metric between 3D V and collection of hits 

std::pair<double,int> VFitter::FitScore2(const std::vector<std::vector<HitLite>>& hits,FittedV& fittedv,bool verbose){

  double score = 0.0;
  int ndof = 0;

  TVector3 arm1_direction_3d = fittedv.GetArm1Dir();
  TVector3 arm2_direction_3d = fittedv.GetArm2Dir();

  for(int i_pl=0;i_pl<kPlaneInvalid;i_pl++){

    if(verbose) std::cout << "Plane " << i_pl << std::endl;

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
      line_1.Direction = TVector2(fittedv.Vertex.Z(),fittedv.Vertex.X());
      line_1.Direction = TVector2(arm1_direction_3d.Z(),arm1_direction_3d.X());
      line_2.Direction = TVector2(fittedv.Vertex.Z(),fittedv.Vertex.X());
      line_2.Direction = TVector2(arm2_direction_3d.Z(),arm2_direction_3d.X());
    }
    double mag = line_1.Direction.Mod();
    line_1.Direction /= mag;
    mag = line_2.Direction.Mod();
    line_2.Direction /= mag;
/*
    if(verbose){
      std::cout << "Line 1: " << line_1.Start.X() << " " << line_1.Start.Y() << " " << line_1.Direction.X() << "  " << line_1.Direction.Y() << std::endl;
      std::cout << "Line 2: " << line_2.Start.X() << " " << line_2.Start.Y() << " " << line_2.Direction.X() << "  " << line_2.Direction.Y() << std::endl;

      TCanvas* c = new TCanvas("c","c");
      TGraphErrors* g = new TGraphErrors(channel_v.at(i_pl).size(),&(channel_v.at(i_pl)[0]),&(tick_v.at(i_pl)[0]),0,&(width_v.at(i_pl)[0]));

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

      c->Print(("fit_transposed_data_Plane" + std::to_string(i_pl) + ".png").c_str());

    }
*/
    for(HitLite hit : hits.at(i_pl)){
      double arm1_sep = HitLineSeparation2(hit,line_1);
      double arm2_sep = HitLineSeparation2(hit,line_2);
      score += std::min(arm1_sep,arm2_sep);
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
// Hybrid - grid seach for vertex but use minimizer to find angles

bool VFitter::DoFitGridSearch3(FittedV& fittedv,int points){

  std::cout << "Starting grid search 3" << std::endl;

  for(size_t i_pl=0;i_pl<3;i_pl++) RemoveOffset(Hits.at(i_pl),InitialGuess.Vertex,i_pl);
  TRandom2* r = new TRandom2();

  std::pair<double,double> vertex_range = {-5.0,5.0};

  double bestfit_score=1e10;
  FittedV bestfit_v;

  for(int i=0;i<points;i++){
    if(i % 50000 == 0) std::cout << "Grid search " << (double)i/points*100 << "\% complete" << std::endl;
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

        return FitVal.first;

    } , 4);

    std::unique_ptr< ROOT::Math::Minimizer > fMinimizer = std::unique_ptr<ROOT::Math::Minimizer>
      ( ROOT::Math::Factory::CreateMinimizer( "Minuit2", "Migrad" ) );

    fMinimizer->SetMaxFunctionCalls(10000);
    fMinimizer->SetTolerance( 0.0001 );
    fMinimizer->SetVariable(0,"Arm 1 Theta",InitialGuess.Arm1Theta,0.005);
    fMinimizer->SetVariable(1,"Arm 1 Phi",InitialGuess.Arm1Phi,0.005);
    fMinimizer->SetVariable(2,"Arm 2 Theta",InitialGuess.Arm2Theta,0.005);
    fMinimizer->SetVariable(3,"Arm 2 Phi",InitialGuess.Arm2Phi,0.005);
    fMinimizer->SetFunction(min);
    fMinimizer->Minimize();

    if(fMinimizer->Status() != 0) continue;

    fittedv.Arm1Theta = fMinimizer->X()[0];
    fittedv.Arm1Phi = fMinimizer->X()[1];
    fittedv.Arm2Theta = fMinimizer->X()[2];
    fittedv.Arm2Phi = fMinimizer->X()[3];

    std::pair<double,int> FitVal = FitScore2(Hits,fittedv);

    if(FitVal.first < bestfit_score){
      bestfit_score = FitVal.first; 
      bestfit_v = fittedv;
      std::cout << "New bestfit" << std::endl;
      std::cout << "Vertex: " << bestfit_v.Vertex.X() << "  " << bestfit_v.Vertex.Y() <<"  " << bestfit_v.Vertex.Z() << std::endl;
      std::cout << "Arm 1: theta=" << bestfit_v.Arm1Theta << "  phi=" << bestfit_v.Arm1Phi << std::endl;
      std::cout << "Arm 2: theta=" << bestfit_v.Arm2Theta << "  phi=" << bestfit_v.Arm2Phi << std::endl;
      std::cout << "Score=" << bestfit_score << std::endl;
    }
  }

  FitScore2(Hits,bestfit_v,true);

  fittedv = bestfit_v;
  fittedv.Vertex = bestfit_v.Vertex + InitialGuess.Vertex;
  //FitScore2(Channel_v,Tick_v,Width_v,fittedv,true); 

  std::cout << "Best Fit:" << std::endl;
  std::cout << "Vertex: " << bestfit_v.Vertex.X() << "  " << bestfit_v.Vertex.Y() <<"  " << bestfit_v.Vertex.Z() << std::endl;
  std::cout << "Arm 1: theta=" << bestfit_v.Arm1Theta << "  phi=" << bestfit_v.Arm1Phi << std::endl;
  std::cout << "Arm 2: theta=" << bestfit_v.Arm2Theta << "  phi=" << bestfit_v.Arm2Phi << std::endl;
  std::cout << "Score=" << bestfit_score << std::endl;

  delete r;

  for(size_t i_pl=0;i_pl<3;i_pl++) RestoreOffset(Hits.at(i_pl),InitialGuess.Vertex,i_pl);

  if(Draw) DrawFit(fittedv);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VFitter::DrawFit(const FittedV& v) const {

  TCanvas* c = new TCanvas("c","c");

  for(int i_pl=0;i_pl<3;i_pl++){

    std::vector<double> channel,tick,width;
    for(HitLite hit : Hits.at(i_pl)){
      channel.push_back(hit.Channel);
      tick.push_back(hit.Tick);
      width.push_back(hit.Width);
    }

    TGraphErrors* g = new TGraphErrors(channel.size(),&(channel[0]),&(tick[0]),0,&(width[0]));
    g->SetName("graph2D");
    g->SetMarkerColor(1);
    g->SetLineColor(1);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.5);
    g->Draw("AP");

    std::pair<TGraph*,TGraph*> fit_arms = v.Make2DVGraphs(i_pl); 

    fit_arms.first->Draw("L same");
    fit_arms.second->Draw("L same");

    c->Print(("Plots/Event_" + RSE + "_Fit_Plane" + std::to_string(i_pl) + ".png").c_str());
    c->Clear();
    delete g;

  }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

}

#endif
