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

#ifndef _WAVEFRONTRECONSTRUCTOR_CPP_
#define _WAVEFRONTRECONSTRUCTOR_CPP_

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <cstdlib>
#include <iomanip>
#include <string>
using std::string;
#include <thread>
using std::thread;
#include <functional>
#include "wavefrontreconstructor.h"

WaveFrontReconstructor::WaveFrontReconstructor()
{
   N=J=0;
   slopes.clear();
   M.clear();
   R.clear();
   VSU.clear();
   coefficients.clear();
   phases.clear();
   haveMat=false;
   cpuTimeMatGen=0.0e0;
   cpuTimeRHSLSE=0.0e0;
}

void WaveFrontReconstructor::GenerateRHSOfLeastSquareEquation(void)
{
   arma::mat U,V;
   arma::vec Sv;
   if (!(arma::svd(U, Sv, V, M))){
      throw string("Error in svd decomposition!");
   }
   int nr=M.n_rows;
   nr-=(J);
   M.clear();
   arma::mat res(J,nr);
   res.zeros();
   arma::mat R=arma::join_rows(arma::pinv(arma::diagmat(Sv)),res);
   res.clear();
   arma::mat Ut=U.t();
   U.clear();
   arma::mat RUt=R*Ut;
   R.clear();
   VSU=(V*RUt);
   V.clear();
   RUt.clear();
   haveMat=true;
}

void WaveFrontReconstructor::ComputeCoefficients(void)
{
   coefficients=VSU*slopes;
}

void WaveFrontReconstructor::PrintCoefficients(void)
{
   int nn=coefficients.size();
   cout << std::scientific << std::setprecision(12) << endl;
   for ( int i=0 ; i<nn ; ++i ) {
      cout << coefficients[i] << endl;
   }
}

void WaveFrontReconstructor::CenterWavefrontAlongZ(vector<double> &zz)
{
   int nn=zz.size();
   if ( nn<=0 ) {
      cerr << "Error: non valid size of array!" << endl;
   }
   double sum=0.0e0;
   for ( int i=0 ; i<nn ; ++i ) { sum+=zz[i]; }
   sum/=double(nn);
   for ( int i=0 ; i<nn ; ++i ) { zz[i]-=sum; }
}

void WaveFrontReconstructor::SetSlopes(const vector<double> &dx,\
      const vector<double> &dy)
{
   size_t nn=dx.size();
   size_t p1=0,p2=nn;
   if ( slopes.size()!=(2*nn) ) { slopes.resize(2*nn); }
   thread t1(&WaveFrontReconstructor::CopyWholeArray,this,std::ref(dx),\
         std::ref(slopes),p1,nn);
   thread t2(&WaveFrontReconstructor::CopyWholeArray,this,std::ref(dy),\
         std::ref(slopes),p2,nn);
   t1.join();
   t2.join();
}

void WaveFrontReconstructor::CopyWholeArray(const vector<double> &vin,\
      arma::vec &vout,size_t posStart,size_t nElem)
{
   size_t nn=vout.size();
   if ( nn<(posStart+nElem) ) {
      cerr << "Error: Not enough space!" << endl;
      cerr << "File: " << __FILE__ << ", line: " << __LINE__ << endl;
   }
   nn=posStart+nElem;
   for ( size_t i=0 ; i<nElem ; ++i ) { vout[i+posStart]=vin[i]; }
}

void WaveFrontReconstructor::ComputeReconstructedWaveFront(const vector<double> &x,
      const vector<double> &y,const vector<double> &coeffs,
      vector<double> &z)
{
   size_t nn=x.size();
   if ( (nn)!=(y.size()) ) {
      cout << "Error: arrays x and y are of different sizes!!" << endl
         << "File: " << __FILE__ << ", line: " << __LINE__ << endl;
   }
   if ( z.size()!=nn ) {
      z.resize(nn);
   }
   if ( nn!=(phases.n_elem) ) {
      cout << "Error: the size of the output array is not\n"
        "equal to the internal size of phases!\n"
        "File: " << __FILE__ << ", line: " << __LINE__  << endl;
   }
   phases=R*coefficients;
   for ( size_t i=0 ; i<nn ; ++i ) { z[i]=phases[i]; }
}


#endif  /* _WAVEFRONTRECONSTRUCTOR_CPP_ */

