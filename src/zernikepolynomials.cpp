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

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <cstdlib>
#include "zernikepolynomials.h"
#include <string>
using std::string;

#ifndef DEBUG
#define DEBUG 0
#endif

#define CIRCULARHARMONICS_EPS_TO_ZERO 1.0e-4

double ZernikePolynomials::pi=acos(-1.0e0);
double ZernikePolynomials::sqrtOneOverPi=sqrt(1.0e0/acos(-1.0e0));
double ZernikePolynomials::sqrtOneOverTwoPi=sqrt(0.5e0/acos(-1.0e0));
bool ZernikePolynomials::initme=true;
double* ZernikePolynomials::factorial=NULL;
double* ZernikePolynomials::invfactorial=NULL;

ZernikePolynomials::ZernikePolynomials()
{
   if (initme) {
      initme=false;
      factorial=new double[RFNMAX];
      invfactorial=new double[RFNMAX];
      factorial[0]=1.0e0;
      invfactorial[0]=1.0e0;
      for (int i=1; i<RFNMAX; i++) {
         factorial[i]=double(i)*factorial[i-1];
         invfactorial[i]=1.0e0/factorial[i];
      }
   }
}

ZernikePolynomials::~ZernikePolynomials()
{
   delete[] factorial;
   delete[] invfactorial;
}

int ZernikePolynomials::ComputePolynomialOrderFromN(int n)
{
   return ((n+1)*(n+2)/2);
}

double ZernikePolynomials::RadialZernike(int n, int m,double rr)
{
   int am=((m<0)?(-m):m);
   if ((n-am)&1) {return 0.0;} //i&1 is a hack for testing oddness of i.
                               // 0&1 is not odd
                               // 1&1 is odd
                               // 2&1 is not odd
                               // 3&1 is odd, etc.
   int sumLimit=(((n-am)>>1)+1);
   double sum=0.0e0,nr=double(n),nm2k;
   for(int k=0; k<sumLimit; ++k) {
      nm2k=nr-2.0e0*double(k);
      sum+=RadialCoefficient(k,n,am)*pow(rr,nm2k);
   }
   return sum;
}

double ZernikePolynomials::DRadialZernikeDr(int n,int m,double r)
{
   int am=((m<0)?(-m):m);
   if ((n-am)&1) {return 0.0;} //i&1 is a hack for testing oddness of i.
                               // 0&1 is not odd
                               // 1&1 is odd
                               // 2&1 is not odd
                               // 3&1 is odd, etc.
   int sumLimit=(((n-am)>>1)+1);
   double sum=0.0,nr=double(n),nm2k;
   for(int k=0; k<sumLimit; ++k) {
      nm2k=nr-2.0e0*double(k);
      sum+=RadialCoefficient(k,n,am)*nm2k*pow(r,(nm2k-1.0e0));
   }
   return sum;
}

double ZernikePolynomials::Zernike(int n,int m,double rr,double phi)
{
   if ( m>0 ) { return RadialZernike(n,m,rr)*cos(double(m)*phi); }
   if ( m<0 ) { return RadialZernike(n,m,rr)*sin(-double(m)*phi); }
   return RadialZernike(n,0,rr);
}

double ZernikePolynomials::DZernikeDr(int n,int m,double r,double p)
{
   if ( m>0 ) { return (DRadialZernikeDr(n,m,r)*cos(double(m)*p)); }
   if ( m<0 ) { return (DRadialZernikeDr(n,m,r)*sin(double(-m)*p)); }
   return DRadialZernikeDr(n,m,r);
}

double ZernikePolynomials::DZernikeDp(int n,int m,double r,double p)
{
   double mr=double(m);
   if ( m>0 ) { return (-mr*RadialZernike(n,m,r)*sin(mr*p)); }
   if ( m<0 ) { return (-mr*RadialZernike(n,m,r)*cos(-mr*p)); }
   return 0.0e0;
}

double ZernikePolynomials::DZernikeDxFromCartesian(int n,int m,double x,double y)
{
   double r,p;
   ComputeRPFromXY(x,y,r,p);
   double oor=1.0e0/r;
   return ((x*oor*DZernikeDr(n,m,r,p))-(y*oor*oor*DZernikeDp(n,m,r,p)));
}

double ZernikePolynomials::DZernikeDyFromCartesian(int n,int m,double x,double y)
{
   double r,p;
   ComputeRPFromXY(x,y,r,p);
   double oor=1.0e0/r;
   return ((y*oor*DZernikeDr(n,m,r,p))+(x*oor*oor*DZernikeDp(n,m,r,p)));
}






