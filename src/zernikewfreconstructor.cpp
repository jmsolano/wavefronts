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
#include <ctime>
#include "zernikewfreconstructor.h"

ZernikeWaveFrontReconstructor::ZernikeWaveFrontReconstructor()
   : WaveFrontReconstructor()
{
   return;
}

ZernikeWaveFrontReconstructor::~ZernikeWaveFrontReconstructor()
{

}


void ZernikeWaveFrontReconstructor::FindWaveFrontCoefficientsFromCoordinatesAndSlopes(\
      const int principalOrder,const vector<double> &x,const vector<double> &y,\
      const vector<double> &dx,const vector<double> &dy,vector<double> &coeffs)
{
   N=principalOrder;
   J=zp.ComputePolynomialOrderFromN(N);
   int nn=x.size();
   slopes.resize(2*nn);
   for ( int i=0 ; i<nn ; ++i ) {
      slopes(i)=dx[i];
      slopes(i+nn)=dy[i];
   }
   clock_t c_init,c_end;
   if ( !haveMat ) {
      c_init=clock();
      GenerateMatrixM(x,y);
      GenerateRHSOfLeastSquareEquation();
      GenerateReconstructionMatrix(x,y);
      c_end=clock();
      cpuTimeMatGen=double(c_end-c_init)/double(CLOCKS_PER_SEC);
   }
   c_init=clock();
   ComputeCoefficients();
   nn=coefficients.n_elem;
   coeffs.resize(nn);
   for ( int i=0 ; i<nn ; ++i ) {
      coeffs[i]=coefficients.at(i);
   }
   c_end=clock();
   cpuTimeRHSLSE=double(c_end-c_init)/double(CLOCKS_PER_SEC);
}

/* // Only in very weird situation a re-implementation of this function
   // would be needed. First try to look for another solution ;)
   // Try to use GenerateReconstructionMatrix first.
void ZernikeWaveFrontReconstructor::ComputeReconstructedWaveFront(const vector<double> &x,\
         const vector<double> &y,const vector<double> &coeffs,vector<double> &z)
{
   double sum;
   size_t nn=x.size();
   if ( z.size()!=nn ) {
      z.resize(nn);
   }
   int jj=coeffs.size();
   for ( size_t i=0 ; i<nn ; ++i ) {
      sum=0.0e0;
      for ( int j=1 ; j<jj ; ++j ) {
         sum+=(coeffs[j]*\
               (zp.NormalZernikeFromCartesian((j+1),x[i],y[i])));
      }
      z[i]=sum;
   }
}
// */

void ZernikeWaveFrontReconstructor::GenerateMatrixM(const vector<double> &x,const vector<double> &y)
{
   int nn=x.size();
   M.resize((2*nn),J);
   for ( int j=1 ; j<=J ; ++j ) {
      for ( int i=0 ; i<nn ; ++i ) {
         M.at(i,j-1)=zp.DNormalZernikeDxFromCartesian(j,x[i],y[i]);
         M.at(i+nn,j-1)=zp.DNormalZernikeDyFromCartesian(j,x[i],y[i]);
      }
   }
}


void ZernikeWaveFrontReconstructor::GenerateReconstructionMatrix(\
      const vector<double> &x,const vector<double> &y)
{
   if ( !haveMat ) {
      cout << "Error: First generate the matrix M, and compute\n"
        "coefficients!" << endl;
      return;
   }
   int nn=x.size();
   int jj=J;
   phases.resize(nn);
   R.resize(nn,J);
   for ( int j=1 ; j<jj ; ++j ) {
      for ( int i=0 ; i<nn ; ++i ) {
         R.at(i,j)=(zp.NormalZernikeFromCartesian((j+1),x[i],y[i]));
      }
   }
}


