R__LOAD_LIBRARY(libVFitter.so);
R__LOAD_LIBRARY($HYP_TOP/lib/libHyperon.so);
R__LOAD_LIBRARY($HYP_TOP/lib/libParticleDict.so);

#include "EventAssembler.h"
#include "VFitter.h"
#include "SpacePointVisualisation.h" 
#include "HoughTransformer.h"

void DrawTrackIDs(){

  const std::string filename = "HyperonTrees.root";
  const int run_to_view = 7012;
  const int subrun_to_view = 1114;
  const int event_to_view = 55704;

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

  Event e = E.GetEvent(ievent);

  // Record which track IDs the decay products correspond to
  int proton_trackid = -1,pion_trackid = -1;
  for(SimParticle p : e.Decay){
    if(p.PDG == 2212) proton_trackid = p.TrackID;
    if(p.PDG == -211) pion_trackid = p.TrackID;
  }

  std::cout << "Proton TrackID: " << proton_trackid << "  Pion TrackID: " << pion_trackid << std::endl;

  // Select the PFP we want to use the hits from
  //RecoParticle pfp = e.TracklikePrimaryDaughters.at(1);
  //hyperonreco::FittedV v = hyperonreco::MakeFittedVGuessTrack(pfp);
  RecoParticle pfp = e.ShowerlikePrimaryDaughters.at(0);
  hyperonreco::FittedV v = hyperonreco::MakeFittedVGuessShower(pfp);

  double roi_size_ch = 100; 
  double roi_size_tick = roi_size_ch*TickPerWire; 

  // Merge hit vectors together an remove hits outside roi 
  std::vector<std::vector<hyperonreco::HitLite>> hits(3);
  hyperonreco::AddHits(hits,pfp.HitChannels,pfp.HitTicks,pfp.HitWidths,pfp.HitTrackIDs);
  hyperonreco::AddHits(hits,e.UnclaimedHitChannels,e.UnclaimedHitTicks,e.UnclaimedHitWidths,e.UnclaimedHitTrackIDs);
  hyperonreco::KeepHitsInROI(TVector3(pfp.X_NoSC,pfp.Y_NoSC,pfp.Z_NoSC),hits,roi_size_ch,roi_size_tick);

  // Make a map between each trackid and hit collections
  std::map<int,std::vector<std::vector<hyperonreco::HitLite>>> m_trackid_hits;
  for(size_t i_pl=0;i_pl<3;i_pl++){
    for(hyperonreco::HitLite hit : hits.at(i_pl)){
      if(m_trackid_hits.find(hit.TrackID) == m_trackid_hits.end()){
        m_trackid_hits[hit.TrackID] = std::vector<std::vector<hyperonreco::HitLite>>(3);
      }  
      m_trackid_hits.at(hit.TrackID).at(i_pl).push_back(hit);
    }
  } 

  TCanvas* c = new TCanvas("c","c");

  // Draw all of the hits for each trackid
  std::vector<std::vector<TGraph*>> graph_2d_vec(3);
  std::map<int,std::vector<std::vector<hyperonreco::HitLite>>>::iterator it;

  for(int i_pl=0;i_pl<3;i_pl++){

    for(it = m_trackid_hits.begin();it != m_trackid_hits.end();it++)
        graph_2d_vec.at(i_pl).push_back(hyperonreco::Make2DGraph(it->second.at(i_pl),it->first,20));

    TMultiGraph* mg = new TMultiGraph();
    for(TGraph* g : graph_2d_vec.at(i_pl))
      mg->Add(g);

    mg->Draw("AP"); 
    c->Print(("Plots/Event_" + std::to_string(e.run) + "_" + std::to_string(e.subrun) + "_" + std::to_string(e.event) +  "_Plane" + std::to_string(i_pl) + ".png").c_str());
    c->Clear();

  }

  hyperonreco::VFitter fitter(true);
  fitter.SetEvent(run_to_view,subrun_to_view,event_to_view);   
  fitter.SetROI(roi_size_ch,roi_size_tick,TVector3(pfp.X_NoSC,pfp.Y_NoSC,pfp.Z_NoSC));
  fitter.AddData(m_trackid_hits.at(proton_trackid));
  fitter.AddData(m_trackid_hits.at(pion_trackid));
  fitter.SetGuess(v);
  fitter.DoFitGridSearch3(v,10000); 

}

