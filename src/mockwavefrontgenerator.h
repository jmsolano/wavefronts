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

#ifndef MOCKWAVEFRONTGENERATOR_H_
#define MOCKWAVEFRONTGENERATOR_H_

#include <vector>
using std::vector;
#include <string>

class MockWaveFunctionGenerator {
public:
   MockWaveFunctionGenerator();
   /** If the grid is nxn, this function returns n.  */
   int Size(void);
   vector<double> x,y,z,dx,dy;
   /** If false, it only saves the coordinates of points:
    * { (x,y) | rMin <= x^2+y^2 <= rMax  }  */
   void OnSquare(bool os) { onSquare=os; }
   /** Generates a simple tilted wave front with slopes mx and my.  */
   void GenerateXYTiltWaveFront(int n,double mx,double my);
   void GenerateXTiltWaveFront(int n,double mx) {return GenerateXYTiltWaveFront(n,mx,0.0e0);}
   void GenerateYTiltWaveFront(int n,double my) {return GenerateXYTiltWaveFront(n,0.0e0,my);}
   /** Returns a Gaussian-shaped wavefront, centered at (x0,y0) with width=width and
    * height=height. n is the number of points per direction of the grid (nxn).  */
   void GenerateGaussianWaveFront(int n, double x0,double y0,double width,double height=1.0e0);
   /** Generates a Super-Gaussian-shaped wavefront, centered at (x0,y0), etc. The order
    * of the superGaussian function is n (a Gaussian function is a super Gaussian function
    * of order 2).  */
   void GenerateSuperGaussianWaveFront(int n, double x0,double y0,int ord,double width,\
         double height=1.0e0);
   /** Generates z=height*exp(-(x^{2N}/width^{2N}+y^{2N}/width^{2N})), and  */
   void GenerateSquaredSuperGaussianWaveFront(int n, double x0,double y0,int ord,double width,\
         double height=1.0e0);
   void GenerateTestF1WaveFront(int n);
   void GenerateTestF2WaveFront(int n);
   void CenterWavefrontAlongZ(void);
   void Save(const char *fileName);
   void Save(std::string fileName) {Save(fileName.c_str());}
   void SetRMin(double rr) { rMin=rr; }
   void SetRMax(double rr) { rMax=rr; }
   void Print(void);
protected:
   void ComputeCoordinatesCircle(int n);
   void ComputeCoordinatesSquare(int n);
   void SetupCoordinates(int n);
   void ComputeXYTiltWaveFront(double mx,double my);
   void ComputeGaussianWaveFront(double x0, double y0, double width,double height);
   void ComputeTestF1WaveFront(void);
   void ComputeTestF2WaveFront(void);
   void ComputeSuperGaussianWaveFront(double x0,double y0,int order,\
         double width,double height);
   void ComputeSquaredSuperGaussianWaveFront(double x0,double y0,int order,\
         double width,double height);
   bool onSquare;
   double rMin,rMax;
};


#endif  /* MOCKWAVEFRONTGENERATOR_H_ */

