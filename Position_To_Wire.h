#ifndef _Position_To_Wire_h_
#define _Position_To_Wire_h_

#include "TVector3.h"

// Convert 3D positions to channels/ticks
// See docdb 25505-v2 for explanation

const double A_w = 3.33328;
const double C_U = 338.140;
const double C_V = 2732.53;
const double C_Y = 4799.19;
const double A_t = 18.2148;
const double C_t = 818.351;

const double TickPerWire = A_t/A_w;

// cos(60) and sin(60)
const double cos60 = 0.5;
const double sin60 = sqrt(3)/2.0;

inline double U_wire(TVector3 pos) { return A_w*(-sin60*pos.Y()+cos60*pos.Z())+C_U; }
inline double V_wire(TVector3 pos) { return A_w*(sin60*pos.Y()+cos60*pos.Z())+C_V; }
inline double Y_wire(TVector3 pos) { return A_w*pos.Z() + C_Y; }
//inline double tick(TVector3 pos) { return A_t*pos.X() + C_t; }
inline double tick(TVector3 pos) { return A_t*pos.X() + C_t - 20; }

inline std::pair<double,double> WireTick(TVector3 pos,int plane){
  if(plane == 0) return std::make_pair(U_wire(pos),tick(pos));
  if(plane == 1) return std::make_pair(V_wire(pos),tick(pos));
  if(plane == 2) return std::make_pair(Y_wire(pos),tick(pos));
  else throw std::invalid_argument("Position_to_Wire: Invalid plane number");
}

inline double dUdt(TVector3 dir){ return A_w/A_t*(-sin60*dir.Y()/dir.X()+cos60*dir.Z()/dir.X()); }
inline double dVdt(TVector3 dir){ return A_w/A_t*(sin60*dir.Y()/dir.X()+cos60*dir.Z()/dir.X()); }
inline double dYdt(TVector3 dir){ return A_w/A_t*dir.Z()/dir.X(); }

inline double AngleU(TVector3 dir){
   bool invert = dir.X() < 0;
   double angle = (180/3.141)*atan(dUdt(dir));
   angle -= 2*(angle-45.0);
   if(invert && angle < 0) angle += 180;
   if(invert && angle > 0) angle -= 180;
   return angle;
}
inline double AngleV(TVector3 dir){
   bool invert = dir.X() < 0;
   double angle = (180/3.141)*atan(dVdt(dir));
   angle -= 2*(angle-45.0);
   if(invert && angle < 0) angle += 180;
   if(invert && angle > 0) angle -= 180;
   return angle;
}
inline double AngleY(TVector3 dir){
   bool invert = dir.X() < 0;
   double angle = (180/3.141)*atan(dYdt(dir));
   angle -= 2*(angle-45.0);
   if(invert && angle < 0) angle += 180;
   if(invert && angle > 0) angle -= 180;
   return angle;
}

inline double PointHitDistanceU(TVector3 point,double channel,double tick){

TVector3 a((tick-C_t)/A_t,0,2*(channel-C_U)/A_w);
TVector3 n(0.0,0.5,sqrt(3)/2);

return ((a-point)-((a-point).Dot(n))*n).Mag();

}

inline double PointHitDistanceV(TVector3 point,double channel,double tick){

TVector3 a((tick-C_t)/A_t,0,2*(channel-C_V)/A_w);
TVector3 n(0.0,0.5,-sqrt(3)/2);

return ((a-point)-((a-point).Dot(n))*n).Mag();

}

inline double PointHitDistanceY(TVector3 point,double channel,double tick){

TVector3 a((tick-C_t)/A_t,0,(channel-C_Y)/A_w);
TVector3 n(0,1.0,0);

return ((a-point)-((a-point).Dot(n))*n).Mag();

}

inline double PointHitDistance(TVector3 point,double channel,double tick,int plane){
  if(plane == 0) return PointHitDistanceU(point,channel,tick); 
  else if(plane == 1) return PointHitDistanceV(point,channel,tick); 
  else if(plane == 2) return PointHitDistanceY(point,channel,tick); 
  else return 1e10; 
}

#endif
