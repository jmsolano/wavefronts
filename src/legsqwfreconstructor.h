#ifndef LEGENDRESQUAREWFRECONSTRUCTOR_H_
#define LEGENDRESQUAREWFRECONSTRUCTOR_H_

#include "wavefrontreconstructor.h"

/* ************************************************************************** */
class LegendreSquareWaveFrontReconstructor : public WaveFrontReconstructor {
/* ************************************************************************** */
public:
   LegendreSquareWaveFrontReconstructor();
   ~LegendreSquareWaveFrontReconstructor();
   /** See wavefrontreconstructor.h for more information  */
   void FindWaveFrontCoefficientsFromCoordinatesAndSlopes(const int principalOrder,\
         const vector<double> &x,const vector<double> &y,const vector<double> &dx,\
         const vector<double> &dy,vector<double> &coeffs);
/* ************************************************************************** */
protected:
   /** See wavefrontreconstructor.h for more information  */
   void GenerateMatrixM(const vector<double> &x,const vector<double> &y);
   /** See wavefrontreconstructor.h for more information  */
   void GenerateReconstructionMatrix(const vector<double> &x,\
         const vector<double> &y);
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* LEGENDRESQUAREWFRECONSTRUCTOR_H_ */

