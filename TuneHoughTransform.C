R__LOAD_LIBRARY(libVFitter.so);
R__LOAD_LIBRARY($HYP_TOP/lib/libHyperon.so);
R__LOAD_LIBRARY($HYP_TOP/lib/libParticleDict.so);

#include "EventAssembler.h"
#include "VFitter.h"
#include "SpacePointVisualisation.h" 
#include "HoughTransformer.h"

void TuneHoughTransform(){

  const std::string filename = "HyperonTrees.root";
  const int run_to_view = 7003;
  const int subrun_to_view = 217;
  const int event_to_view = 10897;

  EventAssembler E(false);
  E.SetFile(filename,"Background");

  // Select the event we want to analyse
  int ievent=0;
  while(ievent<E.GetNEvents()){
    Event e = E.GetEvent(ievent);
    if(e.run == run_to_view && e.subrun == subrun_to_view && e.event == event_to_view)
      break;
    ievent++;
  }

  if(ievent == E.GetNEvents()){
    std::cout << "Event " << run_to_view << " " << subrun_to_view << " " << event_to_view << " not found, exiting" << std::endl;
    return;
  }

  Event e = E.GetEvent(ievent);

  // Record which track IDs the decay products correspond to
  int proton_trackid = -1,pion_trackid = -1;
  for(SimParticle p : e.Decay){
    if(p.PDG == 2212) proton_trackid = p.TrackID;
    if(p.PDG == -211) pion_trackid = p.TrackID;
  }

  std::cout << "Proton TrackID: " << proton_trackid << "  Pion TrackID: " << pion_trackid << std::endl;

  // Select the PFP we want to use the hits from
  RecoParticle pfp = e.ShowerlikePrimaryDaughters.at(0);
  hyperonreco::FittedV v = hyperonreco::MakeFittedVGuessShower(pfp);
  //RecoParticle pfp = e.TracklikePrimaryDaughters.at(1);

  // Merge hit vectors together an remove hits outside roi 
  double roi_size_ch = 100; 
  double roi_size_tick = roi_size_ch*TickPerWire; 
  std::vector<std::vector<hyperonreco::HitLite>> hits(3);
  hyperonreco::AddHits(hits,pfp.HitChannels,pfp.HitTicks,pfp.HitWidths,pfp.HitTrackIDs);
  hyperonreco::AddHits(hits,e.UnclaimedHitChannels,e.UnclaimedHitTicks,e.UnclaimedHitWidths,e.UnclaimedHitTrackIDs);
  hyperonreco::KeepHitsInROI(TVector3(pfp.X_NoSC,pfp.Y_NoSC,pfp.Z_NoSC),hits,roi_size_ch,roi_size_tick);

  // Count how many hits belong to each true particle - needed to get cluster completeness
  std::map<int,std::vector<int>> m_trackid_hits;
  for(size_t i_pl=0;i_pl<3;i_pl++){
    for(hyperonreco::HitLite hit : hits.at(i_pl)){
      if(m_trackid_hits.find(hit.TrackID) == m_trackid_hits.end())
        m_trackid_hits[hit.TrackID] = std::vector<int>(3,0);
      m_trackid_hits.at(hit.TrackID).at(i_pl)++;
    }
  }

  double Best_Score = 0;
  int Best_Tune = -1;
  double Best_theta_bin_size;
  double Best_r_bin_size;
  int Best_peak_size;
  int Best_grouping;
  double Best_MaxNeighbourDist;

  // Run the hough transforms

  TRandom2* r = new TRandom2();

  std::pair<double,double> theta_bin_size_range = {0.01,0.5};
  std::pair<double,double> r_bin_size_range = {0.1,20};
  std::pair<double,double> peak_size_range = {1,5};
  std::pair<double,double> grouping_range = {3,8};
  std::pair<double,double> max_neighbour_dist_range = {0.5,10.0};

  for(int i=0;i<2000;i++){

    //std::vector<std::vector<hyperonreco::HoughTransformPoint>> clusters;

    double theta_bin_size = r->Uniform(theta_bin_size_range.first,theta_bin_size_range.second);
    double r_bin_size = r->Uniform(r_bin_size_range.first,r_bin_size_range.second);
    int peak_size = r->Uniform(peak_size_range.first,peak_size_range.second);
    int grouping = r->Uniform(grouping_range.first,grouping_range.second); 
    double max_neighbour_dist = r->Uniform(max_neighbour_dist_range.first,max_neighbour_dist_range.second);

    std::cout << std::endl;
    std::cout << "Tune " << i << std::endl;
    std::cout << "theta_bin_size=" << theta_bin_size << std::endl;
    std::cout << "r_bin_size=" << r_bin_size << std::endl;
    std::cout << "peak_size=" << peak_size << std::endl;
    std::cout << "grouping=" << grouping << std::endl;

    // Get the tune metrics
    double three_plane_score = 0.0; 

    std::cout << "max_neighbour_dist=" << max_neighbour_dist << std::endl;

    for(int i_pl=0;i_pl<3;i_pl++){

      //std::cout << "i_pl=" << i_pl << std::endl;

      std::pair<double,double> vertex_ch_tick = WireTick(TVector3(pfp.X_NoSC,pfp.Y_NoSC,pfp.Z_NoSC),i_pl); 
      double or_ch = vertex_ch_tick.first;
      double or_ti = vertex_ch_tick.second;

      hyperonreco::HoughTransformer transformer(hits.at(i_pl),i_pl,vertex_ch_tick.first,vertex_ch_tick.second);
      transformer.SetTuneID(i);
      transformer.SetRBinSize(r_bin_size);
      transformer.SetThetaBinSize(theta_bin_size);
      transformer.SetPeakSize(peak_size);
      transformer.SetPointGrouping(grouping);
      transformer.SetMaxNeighbourDist(max_neighbour_dist);

      transformer.MakeTransform2();
      //clusters.push_back(transformer.MakeClusters());
      std::pair<double,double> performance = transformer.GetPerformanceMetrics();
      three_plane_score += performance.first*performance.second;      

    }

    std::cout << "Score=" << three_plane_score << std::endl;

    if(three_plane_score > Best_Score){
      Best_Score = three_plane_score;
      Best_Tune = i; 
      Best_theta_bin_size = theta_bin_size; 
      Best_r_bin_size = r_bin_size; 
      Best_peak_size = peak_size; 
      Best_grouping = grouping; 
      Best_MaxNeighbourDist = max_neighbour_dist; 
    }

  }

  std::cout << "Best tune: " << Best_Tune << std::endl;

  std::vector<std::vector<hyperonreco::HoughTransformPoint>> clusters;
  for(int i_pl=0;i_pl<3;i_pl++){

    //std::cout << "i_pl=" << i_pl << std::endl;

    std::pair<double,double> vertex_ch_tick = WireTick(TVector3(pfp.X_NoSC,pfp.Y_NoSC,pfp.Z_NoSC),i_pl); 
    double or_ch = vertex_ch_tick.first;
    double or_ti = vertex_ch_tick.second;

    hyperonreco::HoughTransformer transformer(hits.at(i_pl),i_pl,vertex_ch_tick.first,vertex_ch_tick.second,true);
    transformer.SetTuneID(Best_Tune);
    transformer.SetRBinSize(Best_r_bin_size);
    transformer.SetThetaBinSize(Best_theta_bin_size);
    transformer.SetPeakSize(Best_peak_size);
    transformer.SetPointGrouping(Best_grouping);
    transformer.SetMaxNeighbourDist(Best_MaxNeighbourDist);

    transformer.SetEvent(run_to_view,subrun_to_view,event_to_view);   

    transformer.MakeTransform2();
    transformer.DrawFits();
    clusters.push_back(transformer.MakeClusters());

  }

}

