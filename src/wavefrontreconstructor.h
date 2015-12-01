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

#ifndef WAVEFRONTRECONSTRUCTOR_H_
#define WAVEFRONTRECONSTRUCTOR_H_

#include <vector>
using std::vector;
#include <armadillo>

/* ************************************************************************** */
class WaveFrontReconstructor {
/* ************************************************************************** */
public:
   WaveFrontReconstructor(); //Initialization of base objects, such as
                             // arma::mat's, arma:vec's, etc.
                             // See private section.
   /** This is the main function for determining the coefficients of the expansion.
    * See reconstructor.cpp in test directory. This function must be implemented
    * within each particular concrete class since the details, including the number
    * of polynomial terms used in the expansion, depend on the specific polynomial
    * type.  */
   virtual void FindWaveFrontCoefficientsFromCoordinatesAndSlopes(const int principalOrder,\
         const vector<double> &x,const vector<double> &y,const vector<double> &dx,\
         const vector<double> &dy,vector<double> &coeffs) = 0;
   /** This function will compute the wavefront (and will store the results in
    * z). This function assumes that the estimation of the coefficients was
    * performed before calling it (the estimation is called by
    * FindWaveFrontCoefficientsFromCoordinatesAndSlopes).  */
   virtual void ComputeReconstructedWaveFront(const vector<double> &x,
         const vector<double> &y,const vector<double> &coeffs,
         vector<double> &z);
   /** Because the algorithm is in essence an integration, the wavefront is determined
    * up to an arbitrary constant. This function remedies somewhat this by centering
    * the resultant wavefront according to the average of zz.  */
   virtual void CenterWavefrontAlongZ(vector<double> &zz);
   void SetSlopes(const vector<double> &dx,const vector<double> &dy);
   /** As the name suggests, this will return the CPU time used for computing
    * all the matrices involved in the process (in seconds).  */
   double getCPUTimeMatrixGeneration(void) {return cpuTimeMatGen;}
   /** As the name suggests, this will return the CPU time employed 
    * to compute the coefficients  */
   double getCPUTimeCoefficientEstimation(void) {return cpuTimeRHSLSE;}
   /** Returns the actual number of terms used in the expansion.  */
   int PolynomialOrder(void) {return J;}
   void PrintCoefficients(void);
   virtual ~WaveFrontReconstructor() {};
/* ************************************************************************** */
protected:
   /** This function (pure virtual) will generate the matrix that will be
    * solved through the SVD algorithm. This matrix is the one with the 
    * derivatives of the polynomials at the pairs (x,y).  */
   virtual void GenerateMatrixM(const vector<double> &x,const vector<double> &y) = 0;
   /** This will generate the matrix used for the reconstruction, after the
    * Least-Squares method was applied. It contains the values of the
    * polynomials at the points (x,y).  */
   virtual void GenerateReconstructionMatrix(const vector<double> &x,\
         const vector<double> &y) = 0;
   /** Auxiliar, intermediate function for solving the reconstruction problem.
    * Usually there would be no need to re-implement this function.  */
   virtual void GenerateRHSOfLeastSquareEquation(void);
   /** Once the function GenerateReconstructionMatrix, and
    * FindWaveFrontCoefficientsFromCoordinatesAndSlopes had been called,
    * the coefficients would be simply computed with this function.
    * Also, usually there would be no need to re-implement this function.
    * */
   virtual void ComputeCoefficients(void);
   int J,N;
   arma::vec slopes;
   arma::mat M;
   arma::mat R;
   arma::mat VSU;
   arma::vec coefficients;
   arma::vec phases;
   bool haveMat;
   double cpuTimeMatGen,cpuTimeRHSLSE;
   void CopyWholeArray(const vector<double> &vin,arma::vec &vout,\
         size_t posStart,size_t nElem);
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* WAVEFRONTRECONSTRUCTOR_H_ */

