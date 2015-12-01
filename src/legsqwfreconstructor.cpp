#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <cstdlib>
#include <ctime>
#include "legendrepolynomial.h"
#include "legsqwfreconstructor.h"

LegendreSquareWaveFrontReconstructor::LegendreSquareWaveFrontReconstructor()
   : WaveFrontReconstructor()
{
   return;
}

LegendreSquareWaveFrontReconstructor::~LegendreSquareWaveFrontReconstructor()
{

}


void LegendreSquareWaveFrontReconstructor::FindWaveFrontCoefficientsFromCoordinatesAndSlopes(\
      const int principalOrder,const vector<double> &x,const vector<double> &y,\
      const vector<double> &dx,const vector<double> &dy,vector<double> &coeffs)
{
   N=principalOrder;
   J=N*N;
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
   //coefficients.print();
}

/* // Only in very weird situation a re-implementation of this function
   // would be needed. First try to look for another solution ;)
   // Try to use GenerateReconstructionMatrix first.
void LegendreSquareWaveFrontReconstructor::ComputeReconstructedWaveFront(const vector<double> &x,\
         const vector<double> &y,const vector<double> &coeffs,vector<double> &z)
{
   double sum;
   size_t nn=x.size();
   if ( z.size()!=nn ) {
      z.resize(nn);
   }
   int n1,n2;
   int jj=coeffs.size();
   double ooN=1.0e0/double(N);
   for ( size_t i=0 ; i<nn ; ++i ) {
      sum=0.0e0;
      for ( int j=0 ; j<jj ; ++j ) {
         n1=floor(j*ooN);
         n2=j%N;
         sum+=(coeffs[j]*\
               (legendrepolynomial::NormalPn(n1,x[i]))*\
               (legendrepolynomial::NormalPn(n2,y[i])));
      }
      z[i]=sum;
   }
}
// */

void LegendreSquareWaveFrontReconstructor::GenerateMatrixM(const vector<double> &x,const vector<double> &y)
{
   int nn=x.size();
   int n1,n2;
   double ooN=1.0e0/double(N);
   M.resize((2*nn),J);
   for ( int j=0 ; j<J ; ++j ) {
      n1=floor(double(j)*ooN);
      n2=j%N;
      for ( int i=0 ; i<nn ; ++i ) {
         M.at(i,j)=(legendrepolynomial::FirstDerivativeNormalPn(n1,x[i]))*\
                   (legendrepolynomial::NormalPn(n2,y[i]));
         M.at(i+nn,j)=(legendrepolynomial::NormalPn(n1,x[i]))*\
                        (legendrepolynomial::FirstDerivativeNormalPn(n2,y[i]));
      }
   }
}

void LegendreSquareWaveFrontReconstructor::GenerateReconstructionMatrix(\
      const vector<double> &x,const vector<double> &y)
{
   if ( !haveMat ) {
      cout << "Error: First generate the matrix M, and compute\n"
        "coefficients!" << endl;
      cout << "File: " << __FILE__ << ", line: " << __LINE__ << endl;
      return;
   }
   int nn=x.size();
   int n1,n2;
   int jj=J;
   phases.resize(nn);
   R.resize(nn,J);
   double ooN=1.0e0/double(N);
   for ( int j=0 ; j<jj ; ++j ) {
      for ( int i=0 ; i<nn ; ++i ) {
         n1=floor(j*ooN);
         n2=j%N;
         R.at(i,j)=((legendrepolynomial::NormalPn(n1,x[i]))*\
               (legendrepolynomial::NormalPn(n2,y[i])));
      }
   }
}



