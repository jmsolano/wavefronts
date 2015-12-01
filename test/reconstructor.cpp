#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <cstdlib>
#include <fstream>
using std::ofstream;
#include <string>
using std::string;
#include <algorithm>    // std::min_element, std::max_element
using std::min_element;
using std::max_element;
#include "../src/mockwavefrontgenerator.h"
#include "../src/zernikewfreconstructor.h"
#include "../src/legsqwfreconstructor.h"

#define NUM_OF_POINTS_PER_DIRECTION 30

// Choose one of the following wavefront types (for the mock wavefront generator)
#ifndef TYPE_OF_MOCK_WAVEFRONT
//#define TYPE_OF_MOCK_WAVEFRONT "SingleCenteredGaussian"
//#define TYPE_OF_MOCK_WAVEFRONT "DoubleOffCenteredGaussians"
//#define TYPE_OF_MOCK_WAVEFRONT "TestFunction1"
//#define TYPE_OF_MOCK_WAVEFRONT "SquareTestFunction1"
//#define TYPE_OF_MOCK_WAVEFRONT "TestFunction2"
//#define TYPE_OF_MOCK_WAVEFRONT "SuperGaussian4"
#define TYPE_OF_MOCK_WAVEFRONT "SuperGaussian8"
//#define TYPE_OF_MOCK_WAVEFRONT "SuperGaussianSquared6"
//#define TYPE_OF_MOCK_WAVEFRONT "SuperGaussianSquared8"
#endif

//Set the polynomial order (based on the leading power of the
// polynomial set (e.g. in Z(n,m), n is the leading power).
#ifndef POLYNOMIAL_ORDER
#define POLYNOMIAL_ORDER 20
#endif

//Choose the type of polynomials that will be used for the reconstruction.
#ifndef BASE_POLYNOMIALS
#define BASE_POLYNOMIALS "Zernike"
//#define BASE_POLYNOMIALS "Legendre"
#endif

#ifndef USEGNUPLOTINTERP
#define USEGNUPLOTINTERP 0
#endif

double timeRecWF,timeMatOps,absNofPolTerms;
double coeffOfDet,sumTot,sumRes;
string gnuplotComments,gnuplotTitle,gnuplotZLabel;

void saveData3D(vector<double> &xx,vector<double> &yy,vector<double> &zz,\
      string fileName);
void mkPlot2DGnuplot(string fdat,string fpdf);
void mkPlot3DGnuplot(string fdat,string fpdf);
double getCoefficientOfDetermination(vector<double> &known,vector<double> &estim,\
      double &sumtot,double &sumres);

int main (int argc, char *argv[])
{
   MockWaveFunctionGenerator wf;
   //Generic wavefront reconstructor (based on wavefrontreconstructor.h)
   WaveFrontReconstructor *wfRec;

   int maxN=POLYNOMIAL_ORDER;
   int nPtsPerCoord=NUM_OF_POINTS_PER_DIRECTION;
   vector<double> coeffs;

   if ( string(BASE_POLYNOMIALS)==string("Zernike") ) {
      wfRec=new ZernikeWaveFrontReconstructor();
   } else if ( string(BASE_POLYNOMIALS)==string("Legendre") ) {
      wfRec=new LegendreSquareWaveFrontReconstructor();
      wf.SetRMax(0.99e0);
      wf.SetRMin(0.0e0);
      wf.OnSquare(true);
   } else {
      cerr << "Non valid reconstruction method!" << endl;
      return EXIT_FAILURE;
   }

   double gHeight=0.5e0;
   double gWidth=1.0e0/8.0e0;
   //Choose which type of wavefront will be used in the program
   //It is desinged in this way to use grep/sed for scripting
   //all types of wavefronts.
   if ( string(TYPE_OF_MOCK_WAVEFRONT)==string("DoubleOffCenteredGaussians") ) {
      double gHeight=0.5e0;
      wf.GenerateGaussianWaveFront(nPtsPerCoord,0.3e0,0.3e0,(1.0e0/8.0e0),gHeight);
      vector<double> ztmp=wf.z,dxtmp=wf.dx,dytmp=wf.dy;
      wf.GenerateGaussianWaveFront(nPtsPerCoord,0.3e0,-0.3e0,(1.0e0/8.0e0),gHeight);
      int kk=wf.z.size();
      for ( int i=0 ; i<kk ; ++i ) { 
         wf.z[i]+=ztmp[i];
         wf.dx[i]+=dxtmp[i];
         wf.dy[i]+=dytmp[i];
      }
   } else if ( string(TYPE_OF_MOCK_WAVEFRONT)==string("TestFunction1") ) {
      wf.GenerateTestF1WaveFront(nPtsPerCoord);
   } else if ( string(TYPE_OF_MOCK_WAVEFRONT)==string("TestFunction2") ) {
      wf.GenerateTestF2WaveFront(nPtsPerCoord);
   } else if ( string(TYPE_OF_MOCK_WAVEFRONT)==string("SingleCenteredGaussian") ) {
      wf.GenerateGaussianWaveFront(nPtsPerCoord,0.0e0,0.0e0,gWidth,gHeight);
   } else if ( string(TYPE_OF_MOCK_WAVEFRONT)==string("SuperGaussian4") ) {
      gWidth*=4.0e0; //gWidth=0.5;
      wf.GenerateSuperGaussianWaveFront(nPtsPerCoord,0.0e0,0.0e0,4,gWidth,gHeight);
   } else if ( string(TYPE_OF_MOCK_WAVEFRONT)==string("SuperGaussian6") ) {
      gHeight=0.02e0;
      gWidth*=4.0e0; //gWidth=0.5;
      wf.GenerateSuperGaussianWaveFront(nPtsPerCoord,0.0e0,0.0e0,6,gWidth,gHeight);
   } else if ( string(TYPE_OF_MOCK_WAVEFRONT)==string("SuperGaussian8") ) {
      gWidth*=4.0e0; //gWidth=0.5;
      wf.GenerateSuperGaussianWaveFront(nPtsPerCoord,0.0e0,0.0e0,8,gWidth,gHeight);
   } else if ( string(TYPE_OF_MOCK_WAVEFRONT)==string("SuperGaussianSquared6") ) {
      wf.OnSquare(true);
      gWidth*=8.0e0; //gWidth=1.0;
      wf.GenerateSquaredSuperGaussianWaveFront(nPtsPerCoord,0.0e0,0.0e0,3,gWidth,gHeight);
   } else if ( string(TYPE_OF_MOCK_WAVEFRONT)==string("SuperGaussianSquared8") ) {
      wf.OnSquare(true);
      gWidth*=6.0e0; //gWidth=1.0;
      wf.GenerateSquaredSuperGaussianWaveFront(nPtsPerCoord,0.0e0,0.0e0,4,gWidth,gHeight);
   } else if ( string(TYPE_OF_MOCK_WAVEFRONT)==string("SquareTestFunction1") ) {
      wf.OnSquare(true);
      wf.SetRMax(0.99);
      wf.GenerateTestF1WaveFront(nPtsPerCoord);
   }
   wf.CenterWavefrontAlongZ();

   //wf.Print();

   clock_t start=clock();   
   wfRec->FindWaveFrontCoefficientsFromCoordinatesAndSlopes(maxN,wf.x,wf.y,wf.dx,wf.dy,coeffs);
   clock_t temp=clock();
   cout << "Computation of coefficients finished!" << endl;

   //The actual number of polynomial terms used in the reconstruction.
   absNofPolTerms=wfRec->PolynomialOrder();

   cout << "Using " << absNofPolTerms << " " << BASE_POLYNOMIALS
        << " Polynomials" << endl;

   timeMatOps=(double(temp-start)/double(CLOCKS_PER_SEC));
   cout << "CPU time to process matrix operations: "
        << timeMatOps << endl;
   cout << "CPU time, after static matrix generation: "
        << wfRec->getCPUTimeCoefficientEstimation() << endl;

   vector<double> zrec;

   start=clock();
   wfRec->ComputeReconstructedWaveFront(wf.x,wf.y,coeffs,zrec);
   wfRec->CenterWavefrontAlongZ(zrec);
   temp=clock();
   timeRecWF=(double(temp-start)/double(CLOCKS_PER_SEC));
   cout << "CPU time to reconstruct wavefront: "
        << timeRecWF << endl;

   coeffOfDet=getCoefficientOfDetermination(wf.z,zrec,sumTot,sumRes);
   cout << "Coefficient of Determination: " << coeffOfDet;
   cout << " (1-(" << sumRes << "/" << sumTot << "))" << endl;

   string basename=TYPE_OF_MOCK_WAVEFRONT;
   basename+=BASE_POLYNOMIALS;
   basename+=string("POrd");
   basename+=std::to_string(POLYNOMIAL_ORDER);

   string recnam=basename+string("-reconst.dat");
   string recpdf=basename+string("-reconst.pdf");
   saveData3D(wf.x,wf.y,zrec,recnam);

   gnuplotTitle="Reconstructed WF (C. of D.: ";
   gnuplotTitle+=std::to_string(coeffOfDet);
   gnuplotTitle+=")";
   mkPlot3DGnuplot(recnam,recpdf);
   gnuplotTitle="";

   double sumor=0.0e0,sumrec=0.0e0;
   int k=0;
   for ( size_t i=0 ; i<zrec.size() ; ++i ) {
      if ( (sqrt(((wf.x[i])*(wf.x[i]))+((wf.y[i])*(wf.y[i])))<1.0e0) ) {
         sumor+=wf.z[i];
         sumrec+=zrec[i];
         ++k;
      }
   }

   double avrec=sumrec/double(k);
   double avor=sumor/double(k);
   double minOrig=*min_element(std::begin(wf.z),std::end(wf.z));
   double maxOrig=*max_element(std::begin(wf.z),std::end(wf.z));
   double rangeOrig=maxOrig-minOrig;
   double minRec=*min_element(std::begin(zrec),std::end(zrec));
   double maxRec=*max_element(std::begin(zrec),std::end(zrec));
   double rangeRec=maxRec-minRec;

   cout << "avor: " << avor << ", avrec: " << avrec << endl;
   cout << " min (or): " << minOrig;
   cout << ";  max (or): " << maxOrig;
   cout << ";  range(or): " << rangeOrig << endl;
   cout << "min (rec): " << minRec;
   cout << "; max (rec): " << maxRec;
   cout << "; range(rec): " << rangeRec << endl;

   for ( size_t i=0; i<zrec.size() ; ++i ) {
      zrec[i]+=(avor-avrec);
      zrec[i]=1.0e0-((zrec[i]+rangeOrig)/(wf.z[i]+rangeOrig));
      zrec[i]*=100.0e0;
      if ( zrec[i]<0.0e0 ) { zrec[i]=-zrec[i]; }
   }

   string dif3dnam=basename+string("-diff3d.dat");
   string dif3dpdf=basename+string("-diff3d.pdf");

   saveData3D(wf.x,wf.y,zrec,dif3dnam);
   gnuplotZLabel="Error (%)";
   gnuplotTitle="Relative error (C. of D.: ";
   gnuplotTitle+=std::to_string(coeffOfDet);
   gnuplotTitle+=")";
   mkPlot3DGnuplot(dif3dnam,dif3dpdf);

   string difnam=basename+string("-diffrec.dat");
   string difpdf=basename+string("-diffrec.pdf");
   ofstream ofil(difnam.c_str());
   for ( size_t i=0 ; i<zrec.size() ; ++i ) {
      ofil << zrec[i] << endl;
   }
   ofil.close();

   mkPlot2DGnuplot(difnam,difpdf);
   gnuplotTitle="";
   gnuplotZLabel="";
   
   
   //*
   string rawbasename=TYPE_OF_MOCK_WAVEFRONT;

   string fnam=rawbasename+string("-generated.dat");
   string pnam=rawbasename+string("-generated.pdf");

   wf.Save(fnam);

   gnuplotTitle="Target WF";
   mkPlot3DGnuplot(fnam,pnam);
   gnuplotTitle="";

   // */

   double maxval=-1.0e+50;
   for ( size_t i=0 ; i<zrec.size() ; ++i ) {
      if ( maxval<fabs(zrec[i]) ) { maxval=fabs(zrec[i]); }
   }

   cout << "Max(zrec): " << maxval << endl;
   delete wfRec;

   return 0;
}

void saveData3D(vector<double> &xx,vector<double> &yy,vector<double> &zz,string fileName)
{
   ofstream dfil(fileName.c_str());
   dfil << "#Number of maximum points per direction (used to build the grid): "
        << NUM_OF_POINTS_PER_DIRECTION << endl;
   dfil << "#Polynomial type: " << BASE_POLYNOMIALS << endl;
   dfil << "#Polynomial order: " << POLYNOMIAL_ORDER << endl;
   dfil << "#Absolute number of polynomials: " << absNofPolTerms << endl;
   dfil << "#CPU time to perform matrix decomposition: " << (1000.0e0*timeMatOps)
        << " milliseconds" << endl;
   dfil << "#CPU time to reconstruct wavefront: " << (1000.0e0*timeRecWF)
        << " milliseconds"<< endl;
   dfil << "#Mock wavefront type: " << TYPE_OF_MOCK_WAVEFRONT << endl;
   //dfil << << endl;
   int nn=xx.size();
   double tmp=xx[0];
   for ( int i=0 ; i<nn ; ++i ) {
      if ( xx[i]!=tmp ) {
         dfil << endl;
         tmp=xx[i];
      }
      dfil << xx[i] << " " << yy[i] << " " << zz[i] << endl;
   }
   dfil.close();
}

void mkPlot2DGnuplot(string fdat,string fpdf)
{
   //size_t pos=fpdf.find("pdf");
   string feps=fpdf.substr(0,(fpdf.length()-3));
   feps+="eps";
   string str="gnuplot -e \"set term postscript eps color enhanced fontscale 1.5";
   str+=string("; set title '")+gnuplotTitle+string("'");
   str+=string("; set ylabel '")+gnuplotZLabel+string("'");
   str+="; set xlabel 'Point index (in grid)'";
   //str+="; set output '|epstopdf --filter --outfile=";
   //str+=fpdf;
   str+="; set output '";
   str+=feps;
   str+="'";
   str+="; plot '";
   str+=fdat;
   str+="' w l notitle\"";
   system(str.c_str());
   cout << "Making plots..." << endl;
   str="dtkeps2pdf ";
   str+=feps;
   str+=(" > /dev/null 2>&1");
   system(str.c_str());
}

void mkPlot3DGnuplot(string fdat,string fpdf)
{
   string cmd="gnuplot -e \"set term postscript eps color enhanced fontscale 1.5";
   cmd+="; set xtics -1,0.4,1.0; set xlabel 'x'";
   cmd+="; set ytics -1,0.4,1.0; set ylabel 'y'";
   cmd+="; set xyplane relative 0.1";
   cmd+=string("; set zlabel '")+gnuplotZLabel+string("' rotate by 90");
   cmd+=string(";set title '")+gnuplotTitle+string("' offset character 0,-2");
   cmd+="; set palette rgbformulae 33,13,10";
#if USEGNUPLOTINTERP
   cmd+="; set dgrid3d 75,75 splines";
   cmd+="; set pm3d at s depthorder border lw 0.25 lt rgb '#000000'";
#else
   cmd+="; set pm3d at s depthorder border lw 0.25 lt rgb 'black'";
#endif
   cmd+="; set output '|epstopdf --filter --outfile=";
   cmd+=fpdf;
   cmd+="'; splot '";
   cmd+=fdat;
   cmd+="' u 1:2:3 w pm3d notitle\"";
   system(cmd.c_str());

}

double getCoefficientOfDetermination(vector<double> &known,vector<double> &estim,\
      double &st,double &sr)
{
   int nn=known.size();
   if ( nn!=int(estim.size()) ) {
      cerr << "Warning: vectors have not the same size" << endl;
      cerr << "file: " << __FILE__ << ", line: " << __LINE__ << endl;
      if ( nn<int(estim.size()) ) {
         nn=int(estim.size());
      }
   }
   double tot=0.0e0;
   double res=0.0e0;
   double av=0.0e0;
   for ( int i=0 ; i<nn ; ++i ) { av+=known[i]; }
   av/=double(nn);
   double tmp;
   for ( int i=0 ; i<nn ; ++i ) {
      tmp=known[i]-av;
      tot+=(tmp*tmp);
      tmp=known[i]-estim[i];
      res+=(tmp*tmp);
   }
   st=tot;
   sr=res;
   return (1.0e0-(res/tot));
}

