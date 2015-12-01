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

#ifndef ZERNIKEWFRECONSTRUCTOR_H_
#define ZERNIKEWFRECONSTRUCTOR_H_

#include "wavefrontreconstructor.h"
#include "zernikepolynomials.h"

/* ************************************************************************** */
class ZernikeWaveFrontReconstructor : public WaveFrontReconstructor {
/* ************************************************************************** */
public:
   ZernikeWaveFrontReconstructor();
   ~ZernikeWaveFrontReconstructor();
   /** See wavefrontreconstructor.h for more information  */
   void FindWaveFrontCoefficientsFromCoordinatesAndSlopes(const int principalOrder,\
         const vector<double> &x,const vector<double> &y,const vector<double> &dx,\
         const vector<double> &dy,vector<double> &coeffs);
/* ************************************************************************** */
protected:
   int N;
   /** See wavefrontreconstructor.h for more information  */
   void GenerateMatrixM(const vector<double> &x,const vector<double> &y);
   ZernikePolynomials zp;
   /** See wavefrontreconstructor.h for more information  */
   void GenerateReconstructionMatrix(const vector<double> &x,\
         const vector<double> &y);
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* ZERNIKEWFRECONSTRUCTOR_H_ */

