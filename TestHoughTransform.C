R__LOAD_LIBRARY(libVFitter.so);
R__LOAD_LIBRARY($HYP_TOP/lib/libHyperon.so);
R__LOAD_LIBRARY($HYP_TOP/lib/libParticleDict.so);

#include "EventAssembler.h"
#include "VFitter.h"
#include "SpacePointVisualisation.h" 
#include "HoughTransformer.h"

void TestHoughTransform(){

  const std::string filename = "HyperonTrees2.root";
  const int run_to_view = 7003;
  const int subrun_to_view = 1668;
  const int event_to_view = 83430;

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
  std::cout << "Event has " << e.TracklikePrimaryDaughters.size() << " tracks and " << e.ShowerlikePrimaryDaughters.size() << " showers" << std::endl;

  // Select the PFP we want to use the hits from
  RecoParticle pfp = e.ShowerlikePrimaryDaughters.at(0);
  hyperonreco::FittedV v = hyperonreco::MakeFittedVGuessShower(pfp);
  //RecoParticle pfp = e.TracklikePrimaryDaughters.at(1);
  //hyperonreco::FittedV v = hyperonreco::MakeFittedVGuessTrack(pfp);

  double roi_size_ch = 100; 
  double roi_size_tick = roi_size_ch*TickPerWire; 

  // Remove any hits outside the ROI
  hyperonreco::KeepHitsInROI(TVector3(pfp.X,pfp.Y,pfp.Z),pfp.HitChannels_Plane0,pfp.HitTicks_Plane0,pfp.HitWidths_Plane0,roi_size_ch,roi_size_tick,0);
  hyperonreco::KeepHitsInROI(TVector3(pfp.X,pfp.Y,pfp.Z),pfp.HitChannels_Plane1,pfp.HitTicks_Plane1,pfp.HitWidths_Plane1,roi_size_ch,roi_size_tick,1);
  hyperonreco::KeepHitsInROI(TVector3(pfp.X,pfp.Y,pfp.Z),pfp.HitChannels_Plane2,pfp.HitTicks_Plane2,pfp.HitWidths_Plane2,roi_size_ch,roi_size_tick,2);
  hyperonreco::KeepHitsInROI(TVector3(pfp.X,pfp.Y,pfp.Z),e.UnclaimedHitChannels_Plane0,e.UnclaimedHitTicks_Plane0,e.UnclaimedHitWidths_Plane0,roi_size_ch,roi_size_tick,0);
  hyperonreco::KeepHitsInROI(TVector3(pfp.X,pfp.Y,pfp.Z),e.UnclaimedHitChannels_Plane1,e.UnclaimedHitTicks_Plane1,e.UnclaimedHitWidths_Plane1,roi_size_ch,roi_size_tick,1);
  hyperonreco::KeepHitsInROI(TVector3(pfp.X,pfp.Y,pfp.Z),e.UnclaimedHitChannels_Plane2,e.UnclaimedHitTicks_Plane2,e.UnclaimedHitWidths_Plane2,roi_size_ch,roi_size_tick,2);

  // make the vectors of hits for the Hough transform
  std::vector<std::vector<double>> channels(3);
  std::vector<std::vector<double>> ticks(3);
  std::vector<std::vector<double>> widths(3);

  channels.at(0) = pfp.HitChannels_Plane0;
  channels.at(0).insert(channels.at(0).end(),e.UnclaimedHitChannels_Plane0.begin(),e.UnclaimedHitChannels_Plane0.end());
  ticks.at(0) = pfp.HitTicks_Plane0;
  ticks.at(0).insert(ticks.at(0).end(),e.UnclaimedHitTicks_Plane0.begin(),e.UnclaimedHitTicks_Plane0.end());
  widths.at(0) = pfp.HitWidths_Plane0;
  widths.at(0).insert(widths.at(0).end(),e.UnclaimedHitWidths_Plane0.begin(),e.UnclaimedHitWidths_Plane0.end());

  channels.at(1) = pfp.HitChannels_Plane1;
  channels.at(1).insert(channels.at(1).end(),e.UnclaimedHitChannels_Plane1.begin(),e.UnclaimedHitChannels_Plane1.end());
  ticks.at(1) = pfp.HitTicks_Plane1;
  ticks.at(1).insert(ticks.at(1).end(),e.UnclaimedHitTicks_Plane1.begin(),e.UnclaimedHitTicks_Plane1.end());
  widths.at(1) = pfp.HitWidths_Plane1;
  widths.at(1).insert(widths.at(1).end(),e.UnclaimedHitWidths_Plane1.begin(),e.UnclaimedHitWidths_Plane1.end());

  channels.at(2) = pfp.HitChannels_Plane2;
  channels.at(2).insert(channels.at(2).end(),e.UnclaimedHitChannels_Plane2.begin(),e.UnclaimedHitChannels_Plane2.end());
  ticks.at(2) = pfp.HitTicks_Plane2;
  ticks.at(2).insert(ticks.at(2).end(),e.UnclaimedHitTicks_Plane2.begin(),e.UnclaimedHitTicks_Plane2.end());
  widths.at(2) = pfp.HitWidths_Plane2;
  widths.at(2).insert(widths.at(2).end(),e.UnclaimedHitWidths_Plane2.begin(),e.UnclaimedHitWidths_Plane2.end());

  std::vector<std::vector<TGraph*>> graph_2d_vec(3);

  // Run the hough transforms
  std::vector<std::vector<hyperonreco::HoughTransformPoint>> clusters;
  for(int i_pl=0;i_pl<3;i_pl++){

    std::cout << "i_pl=" << i_pl << std::endl;

    std::pair<double,double> vertex_ch_tick = WireTick(v.Vertex,i_pl); 
    double or_ch = vertex_ch_tick.first;
    double or_ti = vertex_ch_tick.second;

    double theta_bin_size = 0.3;
    double r_bin_size = 2;
    int peak_size = 2;
    int grouping = 3; 

    hyperonreco::HoughTransformer transformer(channels.at(i_pl),ticks.at(i_pl),widths.at(i_pl),i_pl,vertex_ch_tick.first,vertex_ch_tick.second);
    transformer.SetRBinSize(r_bin_size);
    transformer.SetThetaBinSize(theta_bin_size);
    transformer.SetPeakSize(peak_size);
    transformer.SetPointGrouping(peak_size);
    
    transformer.MakeTransform2();
    transformer.DrawFits();
    clusters.push_back(transformer.MakeClusters());

  }
/*
  // Fit two tracks to the hit collections from the Hough transforms
  hyperonreco::VFitter fitter;

  fitter.SetGuess(v);
  fitter.SetROI(roi_size_ch,roi_size_tick,TVector3(pfp.X_NoSC,pfp.Y_NoSC,pfp.Z_NoSC));
  //fitter.SetActivePlanes();

  fitter.AddData({clusters.at(0).at(0).Channels,clusters.at(1).at(0).Channels,clusters.at(2).at(0).Channels},{clusters.at(0).at(0).Ticks,clusters.at(1).at(0).Ticks,clusters.at(2).at(0).Ticks},{clusters.at(0).at(0).Widths,clusters.at(1).at(0).Widths,clusters.at(2).at(0).Widths});
  fitter.AddData({clusters.at(0).at(1).Channels,clusters.at(1).at(1).Channels,clusters.at(2).at(1).Channels},{clusters.at(0).at(1).Ticks,clusters.at(1).at(1).Ticks,clusters.at(2).at(1).Ticks},{clusters.at(0).at(1).Widths,clusters.at(1).at(1).Widths,clusters.at(2).at(1).Widths});

  for(int i_pl=0;i_pl<3;i_pl++){
    graph_2d_vec.at(i_pl).push_back(hyperonreco::Make2DGraph(clusters.at(i_pl).at(0).Channels,clusters.at(i_pl).at(0).Ticks,clusters.at(i_pl).at(0).Widths,1,1));
    graph_2d_vec.at(i_pl).push_back(hyperonreco::Make2DGraph(clusters.at(i_pl).at(1).Channels,clusters.at(i_pl).at(1).Ticks,clusters.at(i_pl).at(1).Widths,1,1));
  }

  fitter.DoFitGridSearch2(v,100000); 
  fitter.SetGuess(v);
  fitter.DoFit(v); 
   
  //fitter.DoFitGridSearch(v,10000000); 
  hyperonreco::Draw2DGraphs(graph_2d_vec,&v);
*/
}

