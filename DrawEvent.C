R__LOAD_LIBRARY(lib/libVFitter.so);
R__LOAD_LIBRARY($HYP_TOP/lib/libHyperon.so);
R__LOAD_LIBRARY($HYP_TOP/lib/libParticleDict.so);

#include "VFitter.h"
#include "SpacePointVisualisation.h" 
#include "EventAssembler.h"

using namespace hyperonreco;

void DrawEvent(){

  const std::string filename = "HyperonTrees.root";
  const int run_to_view = 7001;
  const int subrun_to_view = 790;
  const int event_to_view = 39513;

  EventAssembler E(false);
  E.SetFile(filename,"Background");

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

  //FitOrganiser organiser;
  //TFile* f_Output = TFile::Open("rootfiles/FitResults.root","RECREATE");
  //TTree* t_Output = new TTree("Output","Output");
  //organiser.SetTreePtr(t_Output);

  Event e = E.GetEvent(ievent);
  //RecoParticle pfp = e.ShowerlikePrimaryDaughters.at(0);
  //FittedV v = MakeFittedVGuessShower(pfp);
  RecoParticle pfp = e.TracklikePrimaryDaughters.at(1);
  FittedV v = MakeFittedVGuessTrack(pfp);

  double roi_size_ch = 50; 
  double roi_size_tick = roi_size_ch*TickPerWire; 

  std::vector<std::vector<HitLite>> hits(3);
  AddHits(hits,pfp.HitChannels,pfp.HitTicks,pfp.HitWidths,pfp.HitTrackIDs,pfp.HitPDGs);
  AddHits(hits,e.UnclaimedHitChannels,e.UnclaimedHitTicks,e.UnclaimedHitWidths,e.UnclaimedHitTrackIDs,e.UnclaimedHitPDGs);
  KeepHitsInROI(TVector3(pfp.X_NoSC,pfp.Y_NoSC,pfp.Z_NoSC),hits,roi_size_ch,roi_size_tick);

  std::vector<std::vector<HoughTransformPoint>> clusters(3);

  for(size_t i_pl=0;i_pl<3;i_pl++){

    if(hits.at(i_pl).size() < 7) continue;

    std::pair<double,double> vertex_ch_tick = WireTick(TVector3(pfp.X_NoSC,pfp.Y_NoSC,pfp.Z_NoSC),i_pl); 
    double or_ch = vertex_ch_tick.first;
    double or_ti = vertex_ch_tick.second;

    HoughTransformer transformer(hits.at(i_pl),i_pl,vertex_ch_tick.first,vertex_ch_tick.second,true);
    transformer.SetVerbosity(1);
    transformer.SetEvent(e.run,e.subrun,e.event);   
    transformer.MakeTransform2();
    transformer.SetPointGrouping(5);
    transformer.SetRBinSize(0.3);
    transformer.SetThetaBinSize(0.03);
    std::vector<HoughTransformPoint> cluster = transformer.MakeClusters();
    //transformer.FindPeaks3();
    //clusters.at(i_pl) = cluster;

  } // i_pl

  //organiser.SetData(e.run,e.subrun,e.event,pfps,clusters,hits);
  //organiser.MakeFitList();

  /*


     VFitter fitter(true);
     fitter.SetEvent(run_to_view,subrun_to_view,event_to_view);   
     fitter.SetROI(roi_size_ch,roi_size_tick,TVector3(pfp.X_NoSC,pfp.Y_NoSC,pfp.Z_NoSC));
     fitter.AddData(hits);
     fitter.SetGuess(v);
     fitter.DoFitGridSearch3(v,50000); 
     */

}
