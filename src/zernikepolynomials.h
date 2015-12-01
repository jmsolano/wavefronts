/*
 * -------------------------------------------------------------------
 *
 *                        Wavefront-related code.
 *
 *        Copyright (c) 2015, Juan Manuel Solano-Altamirano
 *                                 <jmsolanoalt@gmail.com>
 *
 * -------------------------------------------------------------------
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ---------------------------------------------------------------------
 *
 */

#ifndef ZERNIKEPOLYNOMIALS_H_
#define ZERNIKEPOLYNOMIALS_H_

#include <cmath>

class ZernikePolynomials {
public:
   ZernikePolynomials();
   ~ZernikePolynomials();
   inline void ComputeNMFromJ(int j,int &n,int &m){
      int lm=j-1, ln=0;
      while ( lm>ln ) {++ln; lm-=ln;}
      n=ln; m=ln-2*lm;
   }
   int ComputePolynomialOrderFromN(int n);
   double RadialZernike(int n,int m,double rr);
   double Zernike(int n,int m,double rr,double phi);
   inline double Zernike(int j,double rr,double phi){
      int  n,m; ComputeNMFromJ(j,n,m); return Zernike(n,m,rr,phi);
   }
   inline double NormalizationFactor(int n,int m) {
      double ct=((m==0)?sqrtOneOverTwoPi:sqrtOneOverPi);
      return (ct*sqrt(2.0e0*double(n)+2.0e0));
   }
   inline double NormalZernike(int n,int m,double rr,double phi) {
      return (NormalizationFactor(n,m)*Zernike(n,m,rr,phi));
   }
   inline double NormalZernike(int j,double rr,double phi){
      int  n,m; ComputeNMFromJ(j,n,m); return NormalZernike(n,m,rr,phi);
   }
   inline double NormalZernikeFromCartesian(int n,int m,double x,double y){
      double rr=sqrt(x*x+y*y),pp=((rr>0.0e0)?atan2(y,x):0.0e0);
      return NormalZernike(n,m,rr,pp);
   }
   inline double NormalZernikeFromCartesian(int j,double x,double y){
      int  n,m; ComputeNMFromJ(j,n,m); return NormalZernikeFromCartesian(n,m,x,y);
   }
   inline double DNormalZernikeDr(int n,int m,double r,double p) {
      return (NormalizationFactor(n,m)*DZernikeDr(n,m,r,p));
   }
   inline double DNormalZernikeDp(int n,int m,double r,double p) {
      return (NormalizationFactor(n,m)*DZernikeDp(n,m,r,p));
   }
   inline void ComputeRPFromXY(double x,double y,double &r,double &p) {
      r=sqrt(x*x+y*y);
      p=((r>0.0e0?atan2(y,x):0.0e0));
   }
   double DZernikeDxFromCartesian(int n,int m,double x,double y);
   double DZernikeDyFromCartesian(int n,int m,double x,double y);
   inline double DNormalZernikeDxFromCartesian(int n,int m,double x,double y) {
      return (NormalizationFactor(n,m)*DZernikeDxFromCartesian(n,m,x,y));
   }
   inline double DNormalZernikeDyFromCartesian(int n,int m,double x,double y) {
      return (NormalizationFactor(n,m)*DZernikeDyFromCartesian(n,m,x,y));
   }
   inline double DNormalZernikeDxFromCartesian(int j,double x,double y) {
      int n,m; ComputeNMFromJ(j,n,m); return DNormalZernikeDxFromCartesian(n,m,x,y);
   }
   inline double DNormalZernikeDyFromCartesian(int j,double x,double y) {
      int n,m; ComputeNMFromJ(j,n,m); return DNormalZernikeDyFromCartesian(n,m,x,y);
   }
protected:
   static constexpr int RFNMAX=121;
   static double *factorial;
   static double *invfactorial;
   static bool initme;
   inline double Factorial(int n) {return factorial[n];}
   inline double OneOverFactorial(int n) { return invfactorial[n]; }
   inline double RadialCoefficient(int k,int n,int m){
      return (((k&1)?(-1.0e0):1.0e0)*(factorial[n-k]*\
            invfactorial[((n+m)>>1)-k])*invfactorial[k]*invfactorial[((n-m)>>1)-k]); //here i>>1 is i/2.
   }
   double DRadialZernikeDr(int n,int m,double r);
   double DZernikeDr(int n,int m,double r,double p);
   double DZernikeDp(int n,int m,double r,double p);
   static double pi;
   static double sqrtOneOverPi;
   static double sqrtOneOverTwoPi;
};

#endif  /* ZERNIKEPOLYNOMIALS_H_ */

