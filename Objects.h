#ifndef _Objects_h_
#define _Objects_h_

#include "Position_To_Wire.h"

namespace hyperonreco {

  // Hit object that can live outside LArSoft
  struct HitLite {

    HitLite(int plane,double channel,double tick,double width,size_t number,int trackid=-1,int pdg=-1) :
      Plane(plane),Channel(channel),Tick(tick),Width(width),Number(number),TrackID(trackid),PDG(pdg)
    {}

    int Plane;
    double Channel; 
    double Tick;
    double Width;
    size_t Number;
    int TrackID;
    int PDG;

  }; 

  inline std::vector<HitLite> MakeHits(std::vector<double> channels,std::vector<double> ticks,std::vector<double> widths,std::vector<int> trackids,std::vector<int> pdgs,int plane){
    std::vector<HitLite> hits;
    for(size_t i_p=0;i_p<channels.size();i_p++){
      HitLite hit(plane,channels.at(i_p),ticks.at(i_p),widths.at(i_p),i_p,trackids.at(i_p),pdgs.at(i_p));
      hits.push_back(hit);
    }
    return hits;
  }

  inline void AddHits(std::vector<std::vector<HitLite>>& hit_container,std::vector<std::vector<double>> channels,std::vector<std::vector<double>> ticks,std::vector<std::vector<double>> widths,std::vector<std::vector<int>> trackids,std::vector<std::vector<int>> pdgs){
    for(size_t i_pl=0;i_pl<3;i_pl++){
      std::vector<HitLite> thisplane_hits = MakeHits(channels.at(i_pl),ticks.at(i_pl),widths.at(i_pl),trackids.at(i_pl),pdgs.at(i_pl),i_pl);
      for(size_t i_h=0;i_h<thisplane_hits.size();i_h++){
        hit_container.at(i_pl).push_back(thisplane_hits.at(i_h));
      }
    }
  }

  inline void KeepHitsInROI(TVector3 point,std::vector<std::vector<HitLite>>& hits,double roi_size_ch,double roi_size_tick){
    std::vector<std::vector<HitLite>> hits_tmp(3);
    for(int i_pl=0;i_pl<3;i_pl++){
      std::pair<double,double> vertex_ch_tick = WireTick(point,i_pl); 
      for(size_t i_h=0;i_h<hits.at(i_pl).size();i_h++){
        if(abs(hits.at(i_pl).at(i_h).Channel - vertex_ch_tick.first) < roi_size_ch && abs(hits.at(i_pl).at(i_h).Tick - vertex_ch_tick.second) < roi_size_tick)
          hits_tmp.at(i_pl).push_back(hits.at(i_pl).at(i_h));
      }  
    }
    hits = hits_tmp;
  }

  struct HoughTransformPoint {

   HoughTransformPoint(){} 

    HoughTransformPoint(int plane,double r,double theta) :
      Plane(plane),R(r),Theta(theta)
    {}

    int Plane;
    double R;
    double Theta;
    int Height;
    double Chi2;
    std::vector<HitLite> Hits;


    void RemoveDuplicateHits(){
      std::vector<HitLite> hits_tmp;
      for(HitLite hit : Hits){
        bool found = false;
        for(size_t i_h=0;i_h<hits_tmp.size();i_h++)
          if(hit.Number == hits_tmp.at(i_h).Number) found = true;
        if(!found) hits_tmp.push_back(hit);
      }       
      Hits = hits_tmp; 
    }

    // get trackid in the most hits in this cluster and how many hits it appears in
    std::pair<int,int> GetDominantTrackID(){

      std::map<int,int> m_trackid_hits;
      for(HitLite hit : Hits){
        if(m_trackid_hits.find(hit.TrackID) == m_trackid_hits.end())
          m_trackid_hits[hit.TrackID] = 0;
        m_trackid_hits.at(hit.TrackID)++;
      }

      std::map<int,int>::iterator it;
      int max_hits = 0;    
      int dominant_trackid = 0;
      for(it = m_trackid_hits.begin();it != m_trackid_hits.end();it++){
        if(it->second > max_hits){
          max_hits = it->second;
          dominant_trackid = it->first;
        } 
      }

      return std::make_pair(dominant_trackid,max_hits);

    }

  };


}

#endif

