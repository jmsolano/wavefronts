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

#ifndef _MOCKWAVEFRONTGENERATOR_CPP_
#define _MOCKWAVEFRONTGENERATOR_CPP_

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <cstdlib>
#include <iomanip>
#include <string>
using std::string;
using std::to_string;
#include <fstream>
using std::ofstream;
#include <cmath>

#include "mockwavefrontgenerator.h"

#define MIN_NON_ZERO_EPS_FOR_R (1.0e-06)

MockWaveFunctionGenerator::MockWaveFunctionGenerator()
{
   x.clear();
   y.clear();
   z.clear();
   dx.clear();
   dy.clear();
   onSquare=false;
   rMin=0.0e0;
   rMax=1.0e0;
}

void MockWaveFunctionGenerator::SetupCoordinates(int n)
{
   x.clear();
   y.clear();
   z.clear();
   dx.clear();
   dy.clear();
   if ( onSquare ) {
      ComputeCoordinatesSquare(n);
   } else {
      ComputeCoordinatesCircle(n);
   }
   z.resize((x.size()),0.0e0);
   dx.resize((x.size()),0.0e0);
   dy.resize((x.size()),0.0e0);
}

void MockWaveFunctionGenerator::ComputeCoordinatesCircle(int n)
{
   if ( n<=0 ) {
      throw (string("You should use n>0!\n" )\
         +string(__FILE__)+string(": ")+to_string(__LINE__));
      return;
   }
   double dx=2.0e0/double(n-1);
   double cur2x=-1.0e0,cur2y,r2;
   double r2Max=rMax*rMax;
   double r2Min=rMin*rMin;
   for ( int i=0 ; i<n ; ++i ) {
      cur2y=-1.0e0;
      for ( int j=0 ; j<n ; ++j ) {
         r2=cur2x*cur2x+cur2y*cur2y;
         if ( (r2<=r2Max) && (r2>=r2Min) ) { 
            x.push_back(cur2x);
            y.push_back(cur2y);
         }
         cur2y+=dx;
      }
      cur2x+=dx;
   }
   cout << "x.size(): " << x.size() << endl;
}

void MockWaveFunctionGenerator::ComputeCoordinatesSquare(int n)
{
   if ( n<=0 ) {
      throw (string("You should use n>0!\n" )\
         +string(__FILE__)+string(": ")+to_string(__LINE__));
      return;
   }
   double dx=2.0e0*rMax/double(n-1);
   double cur2x=-rMax,cur2y;
   for ( int i=0 ; i<n ; ++i ) {
      cur2y=-rMax;
      for ( int j=0 ; j<n ; ++j ) {
         x.push_back(cur2x);
         y.push_back(cur2y);
         cur2y+=dx;
      }
      cur2x+=dx;
   }
}

int MockWaveFunctionGenerator::Size(void)
{
   return x.size();
}

void MockWaveFunctionGenerator::ComputeXYTiltWaveFront(double mx,double my)
{
   int nn=x.size();
   for ( int i=0 ; i<nn ; ++i ) {
      z[i]=mx*x[i]+my*y[i];
      dx[i]=mx;
      dy[i]=my;
   }
}

void MockWaveFunctionGenerator::GenerateXYTiltWaveFront(int n,double mx,double my)
{
   SetupCoordinates(n);
   ComputeXYTiltWaveFront(mx,my);
}


void MockWaveFunctionGenerator::GenerateGaussianWaveFront(int n,\
      double x0,double y0,double width,double height)
{
   SetupCoordinates(n);
   ComputeGaussianWaveFront(x0,y0,width,height);
}

void MockWaveFunctionGenerator::GenerateSuperGaussianWaveFront(int n,\
      double x0,double y0,int ord,double width,double height)
{
   SetupCoordinates(n);
   ComputeSuperGaussianWaveFront(x0,y0,ord,width,height);
}

void MockWaveFunctionGenerator::GenerateSquaredSuperGaussianWaveFront(int n,\
      double x0,double y0,int ord,double width,double height)
{
   SetupCoordinates(n);
   ComputeSquaredSuperGaussianWaveFront(x0,y0,ord,width,height);
}

void MockWaveFunctionGenerator::ComputeGaussianWaveFront(double x0,double y0,\
      double width,double height)
{
   int nn=x.size();
   double r2;
   double oow2=1.0e0/(width*width);
   double twoow2=2.0e0*oow2;
   for ( int i=0 ; i<nn ; ++i ) {
      r2=((x[i]-x0)*(x[i]-x0));
      r2+=((y[i]-y0)*(y[i]-y0));
      z[i]=height*exp(-(r2*oow2));
      dx[i]=-(twoow2*(x[i]-x0)*z[i]);
      dy[i]=-(twoow2*(y[i]-y0)*z[i]);
   }
}

void MockWaveFunctionGenerator::ComputeSuperGaussianWaveFront(double x0,\
      double y0,int order,double width,double height)
{
   int nn=x.size();
   double r,rn,ee,row;
   for ( int i=0 ; i<nn ; ++i ) {
      r=((x[i]-x0)*(x[i]-x0));
      r+=((y[i]-y0)*(y[i]-y0));
      r=sqrt(r);
      if ( r<MIN_NON_ZERO_EPS_FOR_R ) {
         z[i]=height;
         dx[i]=0.0e0;
         dy[i]=0.0e0;
      } else {
         row=r/width;
         rn=1.0e0;
         for ( int j=0 ; j<order ; ++j ) { rn*=row; }
         ee=exp(-(rn*0.5e0))*height;
         r=1.0e0/r; r*=r; //1/r^2
         z[i]=ee;
         ee*=(-0.5e0*double(order)*rn*r);
         dx[i]=x[i]*ee;
         dy[i]=y[i]*ee;
      }
   }
}

void MockWaveFunctionGenerator::ComputeSquaredSuperGaussianWaveFront(double x0,\
      double y0,int order,double width,double height)
{
   int nn=x.size();
   int efford=order;
   if ( order<1 ) {
      cerr << "Non valid order for square Gaussian beam!" << endl
           << "using N=2..." << endl;
#if DEBUG
      cerr << "file: " << __FILE__ << ", line: " << __LINE__ << endl;
#endif /* ( DEBUG ) */
      efford=2;
   }
   double xx,yy,ee,xow,yow,xn,yn;
   for ( int i=0 ; i<nn ; ++i ) {
      xx=((x[i]-x0)*(x[i]-x0));
      yy=((y[i]-y0)*(y[i]-y0));
      if ( sqrt(xx+yy)<MIN_NON_ZERO_EPS_FOR_R ) {
         z[i]=height;
         dx[i]=0.0e0;
         dy[i]=0.0e0;
      } else {
         xow=xx/(width*width);
         xn=1.0e0;
         for ( int j=0 ; j<efford ; ++j ) { xn*=xow; }
         yow=yy/(width*width);
         yn=1.0e0;
         for ( int j=0 ; j<efford ; ++j ) { yn*=yow; }
         ee=exp(-(xn+yn))*height;
         z[i]=ee;
         ee*=(-double(2*efford));
         dx[i]=xn*ee/(x[i]-x0);
         dy[i]=yn*ee/(y[i]-y0);
      }
   }
}

void MockWaveFunctionGenerator::GenerateTestF1WaveFront(int n)
{
   SetupCoordinates(n);
   ComputeTestF1WaveFront();
}

void MockWaveFunctionGenerator::GenerateTestF2WaveFront(int n)
{
   SetupCoordinates(n);
   ComputeTestF2WaveFront();
}

void MockWaveFunctionGenerator::ComputeTestF2WaveFront()
{
   size_t nn=x.size();
   double tmp1,tmp2,tmp3;
   double height=0.02e0;
   for ( size_t i=0 ; i<nn ; ++i ) {
      tmp1=x[i]*x[i]*x[i]; //x^3
      tmp2=y[i]*y[i]; tmp2*=tmp2; tmp2*=y[i]; //y^5
      z[i]=(-x[i]+20.0e0*tmp1+80.0e0*tmp2);
      z[i]*=(12.0e0*exp(1.0e0+4.0e0*x[i]));
      z[i]-=1.0e0;
      z[i]*=(exp(4.0e0*y[i]));
      tmp3=1.0e0-2.0e0*x[i]; tmp3*=tmp3; //(1-2x)^2
      z[i]+=(9.0e0*exp(4.0e0*x[i])*tmp3);
      tmp3=1.0e0+2.0e0*y[i]; tmp3*=tmp3; //(1+2y)^2
      z[i]*=(exp(-4.0e0*x[i]*(1.0e0+x[i])-tmp3)/15.0e0);
      z[i]*=height;
   }
   for ( size_t i=0 ; i<nn ; ++i ) {
      tmp1=x[i]*x[i]; //x^2
      tmp2=tmp1*tmp1; //x^4
      tmp3=y[i]*y[i]; tmp3*=tmp3; tmp3*=y[i]; //y^5
      dx[i]=(1.0e0-68.0e0*tmp1+160.0e0*tmp2+640.0e0*x[i]*tmp3);
      dx[i]*=(-3.0e0*exp(1.0e0+4.0e0*x[i]));
      dx[i]+=(1.0e0+2.0e0*x[i]);
      dx[i]*=(4.0e0*exp(4.0e0*y[i]));
      dx[i]-=(36.0e0*exp(4.0e0*x[i])*(1.0e0+8.0e0*(x[i]-1.0e0)*tmp1));
      tmp3=1.0e0+2.0e0*y[i]; tmp3*=tmp3; //(1+2y)^2
      dx[i]*=(exp(-4.0e0*x[i]*(1.0e0+x[i])-tmp3)/15.0e0);
      dx[i]*=height;
   }
   for ( size_t i=0 ; i<nn ; ++i ) {
      tmp1=x[i]*x[i]*x[i]; //x^3
      tmp2=y[i]*y[i];
      tmp3=tmp2*y[i]; //y^3
      tmp2*=tmp3; //y^5
      dy[i]=(-x[i]+20.0e0*tmp1-50.0e0*tmp3+80.0e0*tmp2);
      dy[i]*=(-12.0e0*exp(1.0e0+4.0e0*x[i]));
      dy[i]+=1.0e0;
      dy[i]*=(8.0e0*y[i]*exp(4.0e0*y[i]));
      tmp1=1.0e0-2.0e0*x[i]; tmp1*=tmp1; //(1-2x)^2
      tmp2=1.0e0+2.0e0*y[i];
      dy[i]-=(36.0e0*exp(4.0e0*x[i])*tmp1*tmp2);
      dy[i]*=(exp(-4.0e0*x[i]*(1.0e0+x[i])-tmp2*tmp2)/15.0e0);
      dy[i]*=height;
   }
}

void MockWaveFunctionGenerator::ComputeTestF1WaveFront()
{
   size_t nn=x.size();
   double tmp1,tmp2,tmp3;
   for ( size_t i=0 ; i<nn ; ++i ) {
      /*
      tmp1=1.0e0-2.0e0*x[i]; tmp1*=tmp1; //(1-2x)^2
      tmp2=1.0e0+2.0e0*y[i]; tmp2*=tmp2; //(1+2y)^2
      tmp3=4.0e0*x[i]*x[i]; //(2x)^2
      z[i]=3.0e0*tmp1*exp(-tmp3-tmp2);
      tmp1=tmp3; //(2x)^2
      tmp2=2.0e0*tmp1*x[i]; //(2x)^3
      tmp3=4.0e0*y[i]*y[i]; tmp3*=tmp3; tmp3*=(2.0e0*y[i]); //(2y)^5
      z[i]-=(10.0e0*(0.4e0*x[i]-tmp2-tmp3)*exp(-tmp1-4.0e0*y[i]*y[i]));
      tmp1=2.0e0*x[i]+1.0e0; tmp1*=tmp1; //(2x+1)^2
      tmp2=4.0e0*y[i]*y[i];
      z[i]-=(exp(-tmp1-tmp2)/3.0e0);
      z[i]*=0.2e0;
      // */
      //*
      tmp1=x[i]*x[i]*x[i]; //x^3
      tmp2=y[i]*y[i]; tmp2*=tmp2; tmp2*=y[i]; //y^5
      z[i]=(-x[i]+20.0e0*tmp1+80.0e0*tmp2);
      z[i]*=(12.0e0*exp(1.0e0+4.0e0*x[i]));
      z[i]-=1.0e0;
      z[i]*=(exp(4.0e0*y[i]));
      tmp3=1.0e0-2.0e0*x[i]; tmp3*=tmp3; //(1-2x)^2
      z[i]+=(9.0e0*exp(4.0e0*x[i])*tmp3);
      tmp3=1.0e0+2.0e0*y[i]; tmp3*=tmp3; //(1+2y)^2
      z[i]*=(exp(-4.0e0*x[i]*(1.0e0+x[i])-tmp3)/15.0e0);
      // */
   }
   for ( size_t i=0 ; i<nn ; ++i ) {
      tmp1=x[i]*x[i]; //x^2
      tmp2=tmp1*tmp1; //x^4
      tmp3=y[i]*y[i]; tmp3*=tmp3; tmp3*=y[i]; //y^5
      dx[i]=(1.0e0-68.0e0*tmp1+160.0e0*tmp2+640.0e0*x[i]*tmp3);
      dx[i]*=(-3.0e0*exp(1.0e0+4.0e0*x[i]));
      dx[i]+=(1.0e0+2.0e0*x[i]);
      dx[i]*=(4.0e0*exp(4.0e0*y[i]));
      dx[i]-=(36.0e0*exp(4.0e0*x[i])*(1.0e0+8.0e0*(x[i]-1.0e0)*tmp1));
      tmp3=1.0e0+2.0e0*y[i]; tmp3*=tmp3; //(1+2y)^2
      dx[i]*=(exp(-4.0e0*x[i]*(1.0e0+x[i])-tmp3)/15.0e0);
   }
   for ( size_t i=0 ; i<nn ; ++i ) {
      tmp1=x[i]*x[i]*x[i]; //x^3
      tmp2=y[i]*y[i];
      tmp3=tmp2*y[i]; //y^3
      tmp2*=tmp3; //y^5
      dy[i]=(-x[i]+20.0e0*tmp1-50.0e0*tmp3+80.0e0*tmp2);
      dy[i]*=(-12.0e0*exp(1.0e0+4.0e0*x[i]));
      dy[i]+=1.0e0;
      dy[i]*=(8.0e0*y[i]*exp(4.0e0*y[i]));
      tmp1=1.0e0-2.0e0*x[i]; tmp1*=tmp1; //(1-2x)^2
      tmp2=1.0e0+2.0e0*y[i];
      dy[i]-=(36.0e0*exp(4.0e0*x[i])*tmp1*tmp2);
      dy[i]*=(exp(-4.0e0*x[i]*(1.0e0+x[i])-tmp2*tmp2)/15.0e0);
   }
}

void MockWaveFunctionGenerator::Save(const char *fileName)
{
   ofstream ofil(fileName);
   int nn=x.size();
   ofil << std::scientific << std::setprecision(12);
   double tmp=x[0];
   for ( int i=0 ; i<nn ; ++i ) {
      if ( x[i]!=tmp ) {
         ofil << endl;
         tmp=x[i];
      }
      ofil << x[i] << " " << y[i] << " " << z[i];
      ofil << " " << dx[i] << " " << dy[i];
      ofil  << endl;
   }
   ofil.close();
}

void MockWaveFunctionGenerator::Print(void)
{
   int nn=x.size();
   cout << "     x             y                z                dx              dy" << endl;
   for ( int i=0 ; i<80 ; ++i ) { putchar('*'); } cout << endl;
   cout << std::scientific << std::setprecision(8);
   for ( int i =0 ; i<nn ; ++i ) {
      cout << std::setw(16)
           << x[i] << std::setw(16)
           << y[i] << std::setw(16) 
           << z[i] << std::setw(16)
           << dx[i] << std::setw(16)
           << dy[i] << endl;
   }
}


void MockWaveFunctionGenerator::CenterWavefrontAlongZ(void)
{
   int nn=z.size();
   if ( nn<=0 ) {
      cerr << "Error: non valid size of array!" << endl;
   }
   double sum=0.0e0;
   for ( int i=0 ; i<nn ; ++i ) { sum+=z[i]; }
   sum/=double(nn);
   for ( int i=0 ; i<nn ; ++i ) { z[i]-=sum; }
}




#endif  /* _MOCKWAVEFRONTGENERATOR_CPP_ */

