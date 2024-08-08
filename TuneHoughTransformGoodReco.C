R__LOAD_LIBRARY(libVFitter.so);
R__LOAD_LIBRARY($HYP_TOP/lib/libHyperon.so);
R__LOAD_LIBRARY($HYP_TOP/lib/libParticleDict.so);

#include "EventAssembler.h"
#include "VFitter.h"
#include "SpacePointVisualisation.h" 
#include "HoughTransformer.h"
#include <ctime>    
#include <chrono>

const int tunes = 1000;
std::vector<double> theta_bin_size;
std::vector<double> r_bin_size;
std::vector<int> peak_size;
std::vector<int> grouping;
std::vector<double> max_neighbour_dist;
std::vector<double> chi2_cut;
std::vector<double> score(tunes,0);

// Tune Hough transform algorithm on events with already reconstructed Lambdas

void TuneHoughTransformGoodReco(){

  const std::string filename = "GoodReco.root";

  // number of pfps the transformer has been tested - use to get the average performance metrics
  int n_tests = 0;

  // Prepare all the tunes
  auto start = std::chrono::system_clock::now();
  TRandom2* r = new TRandom2(std::chrono::system_clock::to_time_t(start));
  //TRandom2* r = new TRandom2();

  std::pair<double,double> theta_bin_size_range = {0.01,0.5};
  std::pair<double,double> r_bin_size_range = {0.1,20};
  std::pair<double,double> peak_size_range = {1,5};
  std::pair<double,double> grouping_range = {2.5,7.5};
  std::pair<double,double> max_neighbour_dist_range = {5.0,15.0};
  std::pair<double,double> chi2_cut_range = {0.1,2.0};

  for(int i=0;i<tunes;i++){
    theta_bin_size.push_back(r->Uniform(theta_bin_size_range.first,theta_bin_size_range.second)); 
    r_bin_size.push_back(r->Uniform(r_bin_size_range.first,r_bin_size_range.second)); 
    peak_size.push_back(r->Uniform(peak_size_range.first,peak_size_range.second)); 
    grouping.push_back(r->Uniform(grouping_range.first,grouping_range.second)); 
    max_neighbour_dist.push_back(r->Uniform(max_neighbour_dist_range.first,max_neighbour_dist_range.second)); 
    chi2_cut.push_back(r->Uniform(chi2_cut_range.first,chi2_cut_range.second));
  }
/*
  theta_bin_size.back() = 0.0155816;
  r_bin_size.back() = 11.8015;
  peak_size.back() = 1;
  grouping.back() = 6;
  max_neighbour_dist.back() = 12.9798;
  chi2_cut.back() = 1.12616;
*/
  // Merge hit vectors together an remove hits outside roi 
  double roi_size_ch = 100; 
  double roi_size_tick = roi_size_ch*TickPerWire; 

  EventAssembler E(false);
  E.SetFile(filename,"Background");

  // Select the event we want to analyse
  int ievent=0;
  while(ievent<E.GetNEvents()){
    Event e = E.GetEvent(ievent);

    std::cout << ievent << "/" << E.GetNEvents() << std::endl;

    std::tuple<int,int,int> event = std::make_tuple(e.run,e.subrun,e.event);

    ievent++;

    std::vector<RecoParticle> pfps;
    bool has_proton=false,has_pion=false;
    for(RecoParticle pfp : e.TracklikePrimaryDaughters){
      if(pfp.TrackTruePDG == 2212 && pfp.TrackTrueOrigin == 2 && !has_proton){
        pfps.push_back(pfp);
        has_proton = true;
      }
      if(pfp.TrackTruePDG == -211 && pfp.TrackTrueOrigin == 2 && !has_pion){
        pfps.push_back(pfp);
        has_pion = true;
      }
      if(has_proton && has_pion) break;
    }

    if(!has_proton || !has_pion) continue;

    //std::vector<RecoParticle> pfps = e.TracklikePrimaryDaughters;

    std::vector<std::vector<hyperonreco::HitLite>> hits(3);

    for(size_t i_pfp=0;i_pfp<pfps.size();i_pfp++){ 
      const RecoParticle& pfp = pfps.at(i_pfp);
      hyperonreco::AddHits(hits,pfp.HitChannels,pfp.HitTicks,pfp.HitWidths,pfp.HitTrackIDs,pfp.HitPDGs);
      hyperonreco::KeepHitsInROI(TVector3(pfp.X_NoSC,pfp.Y_NoSC,pfp.Z_NoSC),hits,roi_size_ch,roi_size_tick);
    }

    n_tests++;

    for(int i_pl=0;i_pl<3;i_pl++){

      if(hits.at(i_pl).size() < 7) continue;

      //std::cout << "i_pl=" << i_pl << std::endl;

      std::pair<double,double> vertex_ch_tick = WireTick(TVector3(pfps.at(0).X_NoSC,pfps.at(0).Y_NoSC,pfps.at(0).Z_NoSC),i_pl); 
      double or_ch = vertex_ch_tick.first;
      double or_ti = vertex_ch_tick.second;

      hyperonreco::HoughTransformer transformer(hits.at(i_pl),i_pl,vertex_ch_tick.first,vertex_ch_tick.second,false);

      for(int i=0;i<tunes;i++){
        /*
        std::cout << "Testing tune:" << std::endl;
        std::cout << "theta_bin_size=" << theta_bin_size.at(i) << std::endl;
        std::cout << "r_bin_size=" << r_bin_size.at(i) << std::endl;
        std::cout << "peak_size=" << peak_size.at(i) << std::endl;
        std::cout << "grouping=" << grouping.at(i) << std::endl;
        std::cout << "max_neighbour_dist=" << max_neighbour_dist.at(i) << std::endl;
        std::cout << "chi2_cut=" << chi2_cut.at(i) << std::endl;
        */
        transformer.SetTuneID(i);
        transformer.SetRBinSize(r_bin_size.at(i));
        transformer.SetThetaBinSize(theta_bin_size.at(i));
        transformer.SetPeakSize(peak_size.at(i));
        transformer.SetPointGrouping(grouping.at(i));
        transformer.SetMaxNeighbourDist(max_neighbour_dist.at(i));
        transformer.SetChi2Cut(chi2_cut.at(i));
        transformer.MakeTransform2();
        std::pair<double,double> performance = transformer.GetPerformanceMetrics();
        score.at(i) += performance.first*performance.second;      

      }

    } // i_pl

  } // ievent

  // Find the best tune
  double best_score = 0;
  int best_tune = -1;
  for(size_t i=0;i<tunes;i++){
    if(score.at(i) > best_score){
      best_score = score.at(i);
      best_tune = i;
    } 
  }

  std::cout << "Best score: " << best_score/n_tests << std::endl;
  std::cout << "Draw results of best tune" << std::endl;

  std::cout << "theta_bin_size=" << theta_bin_size.at(best_tune) << std::endl;
  std::cout << "r_bin_size=" << r_bin_size.at(best_tune) << std::endl;
  std::cout << "peak_size=" << peak_size.at(best_tune) << std::endl;
  std::cout << "grouping=" << grouping.at(best_tune) << std::endl;
  std::cout << "max_neighbour_dist=" << max_neighbour_dist.at(best_tune) << std::endl;
  std::cout << "chi2_cut=" << chi2_cut.at(best_tune) << std::endl;

  ievent=0;
  while(ievent<E.GetNEvents()){

    Event e = E.GetEvent(ievent);

    std::tuple<int,int,int> event = std::make_tuple(e.run,e.subrun,e.event);

    std::vector<RecoParticle> pfps;
    bool has_proton=false,has_pion=false;
    for(RecoParticle pfp : e.TracklikePrimaryDaughters){
      if(pfp.TrackTruePDG == 2212 && pfp.TrackTrueOrigin == 2 && !has_proton){
        pfps.push_back(pfp);
        has_proton = true;
      }
      if(pfp.TrackTruePDG == -211 && pfp.TrackTrueOrigin == 2 && !has_pion){
        pfps.push_back(pfp);
        has_pion = true;
      }
      if(has_proton && has_pion) break;
    }

    if(!has_proton || !has_pion) continue;

    ievent++;

    std::vector<std::vector<hyperonreco::HitLite>> hits(3);
    for(size_t i_pfp=0;i_pfp<pfps.size();i_pfp++){ 
      const RecoParticle& pfp = pfps.at(i_pfp);
      hyperonreco::AddHits(hits,pfp.HitChannels,pfp.HitTicks,pfp.HitWidths,pfp.HitTrackIDs,pfp.HitPDGs);
      hyperonreco::KeepHitsInROI(TVector3(pfp.X_NoSC,pfp.Y_NoSC,pfp.Z_NoSC),hits,roi_size_ch,roi_size_tick);
    }

    for(int i_pl=0;i_pl<3;i_pl++){

      if(hits.at(i_pl).size() < 7) continue;

      std::pair<double,double> vertex_ch_tick = WireTick(TVector3(pfps.at(0).X_NoSC,pfps.at(0).Y_NoSC,pfps.at(0).Z_NoSC),i_pl); 
      double or_ch = vertex_ch_tick.first;
      double or_ti = vertex_ch_tick.second;

      hyperonreco::HoughTransformer transformer(hits.at(i_pl),i_pl,vertex_ch_tick.first,vertex_ch_tick.second,true);
      transformer.SetTuneID(best_tune);
      transformer.SetRBinSize(r_bin_size.at(best_tune));
      transformer.SetThetaBinSize(theta_bin_size.at(best_tune));
      transformer.SetPeakSize(peak_size.at(best_tune));
      transformer.SetPointGrouping(grouping.at(best_tune));
      transformer.SetMaxNeighbourDist(max_neighbour_dist.at(best_tune));
      transformer.SetChi2Cut(chi2_cut.at(best_tune));
      transformer.SetEvent(e.run,e.subrun,e.event);   
      transformer.MakeTransform2();
      transformer.MakeClusters();

    } // i_pl

  } // ievent

  E.Close();


  gSystem->Exec("mkdir -p rootfiles");
  TFile* f_tunes = new TFile("rootfiles/Tunes.root","RECREATE");
  TTree* t_tunes = new TTree("tunes","tunes");

  double t_theta_bin_size;
  double t_r_bin_size;
  int t_peak_size;
  int t_grouping;
  double t_max_neighbour_dist;
  double t_chi2_cut;
  double t_score;
  t_tunes->Branch("theta_bin_size",&t_theta_bin_size);
  t_tunes->Branch("r_bin_size",&t_r_bin_size);
  t_tunes->Branch("peak_size",&t_peak_size);
  t_tunes->Branch("grouping",&t_grouping);
  t_tunes->Branch("max_neighbour_dist",&t_max_neighbour_dist);
  t_tunes->Branch("chi2_cut",&t_chi2_cut);
  t_tunes->Branch("score",&t_score);

  for(int i=0;i<tunes;i++){
    t_theta_bin_size = theta_bin_size.at(i);
    t_r_bin_size = r_bin_size.at(i);
    t_peak_size = peak_size.at(i);
    t_grouping = grouping.at(i);
    t_max_neighbour_dist = max_neighbour_dist.at(i);
    t_chi2_cut = chi2_cut.at(i);
    t_score = score.at(i)/n_tests;
    t_tunes->Fill();
  }

  t_tunes->Write("tunes");

  f_tunes->Close();

}

