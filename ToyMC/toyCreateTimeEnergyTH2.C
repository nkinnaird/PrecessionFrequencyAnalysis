// Inspiration and various equations taken from Aaron's note: E989 Note 123: Characterization of an Energy Binned Precession Frequency Analysis, docDB 8789

// In order to run this first compile the macro as normal (root toyCreateTimeEnergyTH2.C+
// for i in {0..999}; do echo $i; done | xargs -i --max-procs=20 bash -c "root -b -q -l ../toyCreateTimeEnergyTH2.C\({}\)"
// where there is no + so the macro isn't recompiled, and where 999 is the number of jobs here, the max processes can be changed, and the path to the root macro might have to be changed depending on where this is ran.


#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TF2.h>
#include <TF1Convolution.h>
#include <TH1.h>
#include <TH2.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TImage.h>
#include <TRandom3.h>
#include <sstream>
#include <TVirtualFFT.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TVectorD.h>

#include "ratioAnalysisDefs.hh"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////

int numberOfJobs = 1000; // number of jobs in which to construct the full 2D histogram

// base global functions so that convoluted functions can access them
TF1* R_func; // resolution function
TF1* Acc_func; // acceptance
TF1* N_theoretical;
TF1* N_with_Ac; // N with acceptance
TF1* A_theoretical;
TF1* tdrift_energy_dependent; // drift time which has energy dependence

double resolutionModes[3990]; // x values of peaks in energy dependent resolution function, defined from E = 10 MeV to 4000 MeV

TF1* N_detected; // convolution functions - not directly used in full 2D function since everything is convoluted together at once
TF1* A_detected;

TF2* full2D_convolutedFunc;

/////////////////////////////////////////////////////////////////////////////////////

// variables 

double defaultPhase = 2.09;
int nBins = int(approxMaxTime/defaultBinWidth);
double histMaxTime = nBins*defaultBinWidth;

double function_Wa = blindingWa;


/////////////////////////////////////////////////////////////////////////////////////

double N_times_Ac(double* x, double* p)
{
  //x[0] is the true energy E
  return Acc_func->Eval(x[0]) * N_theoretical->Eval(x[0]);
}

double fullTheoreticalFunc(double* x, double* p)
{
  //x[0] is time t
  //x[1] is energy Ed
  return 1000.*exp(-x[0]/defaultLifetime) * N_with_Ac->Eval(x[1]) * (1 + A_theoretical->Eval(x[1]) * cos(function_Wa * x[0] - defaultPhase));
}

double resolutionGaus(double* x, double* p)
{
  //p[0] is the true energy E

  double sig = sqrt(.015*.015 + .045*.045/(p[0]/1000))*(p[0]);
  double R = (1/sqrt(2*pi*sig*sig)) * exp(-(x[0])*(x[0])/(2*sig*sig));

  return R;
}

double resolutionExpPlusGaus(double* x, double* p)
{
  //p[0] is the true energy E

  double sig = sqrt(.015*.015 + .045*.045/(p[0]/1000))*(p[0]);
  double R;

  // this needs to be improved if used - with (1-t)*func1 + t*func2 or something to splice them

  if(x[0] > -sig) R = (1/sqrt(2*pi*sig*sig)) * exp(-(x[0])*(x[0])/(2*sig*sig));
  else R = (1/sqrt(2*pi*sig*sig)) * exp(x[0]/(2*sig));

  // R = (1/sqrt(2*pi*sig*sig)) * exp(-(x[0])*(x[0])/(2*sig*sig));
  // if(x[0] < -sig) R += (1/sqrt(2*pi*sig*sig)) * exp(x[0]/(2*sig));

  return R;
}

double resolutionExpModGaus(double* x, double* p)
{
  // p[0] is the true energy E
  // p[1] is the mode value to shift the exp modified gaussian over - so the peak can be recentered at 0

  double sig = sqrt(.015*.015 + .045*.045/(p[0]/1000))*(p[0]);
  double lambda = .03; // value chosen by eye 
  double mode = p[1];
  double xval = x[0] + mode;

  double R;

  R = (lambda/2) * exp((lambda/2) * (lambda*sig*sig + 2*xval)) * erfc((lambda*sig*sig+xval)/(sqrt(2)*sig));
  if(isnan(R)) R = 0;

  return R;
}

double convoluteN(double* x, double* p)
{
  double Nd = 0;
  double Ed = x[0];

  for (int E = 10; E < 3100; ++E)
  {
     R_func->SetParameter(0, E);
     R_func->SetParameter(1, resolutionModes[E-10]);
     Nd += N_with_Ac->Eval(E) * R_func->Eval(Ed - E);
  }

  return Nd / p[0];
}

double convoluteA(double* x, double* p)
{
  double Ad = 0;
  double Ed = x[0];

  for (int E = 10; E < 3100; ++E)
  {
     R_func->SetParameter(0, E);
     R_func->SetParameter(1, resolutionModes[E-10]);
     Ad += (N_with_Ac->Eval(E) * A_theoretical->Eval(E)) * R_func->Eval(Ed - E);
  }

  return (Ad / N_detected->Eval(Ed)); // in order to get convolution of A with R by itself, need to divide out the N part (since it's all convoluted together)
}

double tdriftFunc(double* x, double* p)
{
  double energyFraction = x[0]/3100.;

  double tdrift_polynomial = (-0.255134 + 65.3034*energyFraction - 705.492*pow(energyFraction,2) + 5267.21*pow(energyFraction,3) - 
                              2.39865e4*pow(energyFraction,4) + 6.83481e4*pow(energyFraction,5) - 1.21761e5*pow(energyFraction,6) + 
                              1.31393e5*pow(energyFraction,7) - 7.8343e4*pow(energyFraction,8) + 1.97741e4*pow(energyFraction,9)) * 1e-3;
                              // polynomial values from Nandita, 1e-3 at end to convert from mrad to rad

  if(energyFraction > 1) tdrift_polynomial = 0; // no particles with energy > 3100 MeV

  return tdrift_polynomial;
} // I'm not sure how to get the convolution of R and the drift time / energy dependent phase by itself since it's within the cosine and everything is convoluted together

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

double fullFunctionConv(double* x, double* p)
{
  //x[0] is time t
  //x[1] is energy Ed

  double convFuncReturn = 0;
  double Ed = x[1];

  for (int E = 10; E < 3100; ++E)
  {
     R_func->SetParameter(0, E);
     R_func->SetParameter(1, resolutionModes[E-10]);
     convFuncReturn += R_func->Eval(Ed - E) * ( N_with_Ac->Eval(E) * (1 + A_theoretical->Eval(E) * cos(function_Wa * x[0] - function_Wa * tdrift_energy_dependent->Eval(E) - defaultPhase)) );
  }

  return 1000.*exp(-x[0]/defaultLifetime) * convFuncReturn;
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

void checkCanvases()
{
    cout << endl << "Creating canvases to check functions." << endl;

/////////////////////////////////////////////////////////////////////////////////////

    // plot R and various other resolution functions
  
      double plotEnergy = 1000; // set energy just for plotting
      R_func->SetParameter(0, plotEnergy);
      R_func->SetParameter(1, 0);

      double mode = R_func->GetMaximumX(-33,0);
      R_func->SetParameter(1,mode); // shift plot so peak R is centered at 0

        TF1* R_func_gaus = new TF1("R_func_gaus", resolutionGaus, -4000, 4000, 1); // Gaussian resolution function for comparison
        R_func_gaus->SetParameter(0,plotEnergy);
        R_func_gaus->SetNpx(8000);
        R_func_gaus->SetLineColor(2);

        // TF1* R_func_exp_plus_gaus = new TF1("R_func_exp_plus_gaus", resolutionExpPlusGaus, -4000, 4000, 1); // Gaussian spliced with an exponential tail - not fully worked out
        // R_func_exp_plus_gaus->SetParameter(0,plotEnergy);
        // R_func_exp_plus_gaus->SetNpx(8000);
        // R_func_exp_plus_gaus->SetLineColor(4);

      auto R_canvas = new TCanvas("R_canvas","R_canvas",200,10,1200,800);
      R_func->Draw();
      R_func_gaus->Draw("SAME");
      // R_func_exp_plus_gaus->Draw("SAME");

      R_func->GetYaxis()->SetRangeUser(0, 1.1*R_func_gaus->GetMaximum()); // fix plotting range so nothing is cut off
      R_canvas->Update();

/////////////////////////////////////////////////////////////////////////////////////

    // plot N and N_detected

    N_detected = new TF1("N_detected", convoluteN, 0, 3300, 1);
    N_detected->SetNpx(3300);
    N_detected->SetLineColor(2);
    N_detected->SetParameter(0, 1);
    // double norm_N_d = N_detected->Integral(10, 3300); // can normalize the plots if wanted
    // N_detected->SetParameter(0, norm_N_d); 

    auto N_canvas = new TCanvas("N_canvas","N_canvas",200,10,1200,800);
    N_theoretical->Draw();
    N_detected->Draw("SAME");

/////////////////////////////////////////////////////////////////////////////////////

    // plot A and A_detected

    A_detected = new TF1("A_detected", convoluteA, 0, 3300, 0);
    A_detected->SetNpx(3300);
    A_detected->SetLineColor(2);

    auto A_canvas = new TCanvas("A_canvas","A_canvas",200,10,1200,800);
    A_theoretical->Draw();
    A_detected->Draw("SAME");

/////////////////////////////////////////////////////////////////////////////////////

    // plot of the drift time with energy dependence

    auto t_canvas = new TCanvas("t_canvas","t_canvas",200,10,1200,800);
    tdrift_energy_dependent->Draw();

/////////////////////////////////////////////////////////////////////////////////////

    auto theoretical2D_func = new TF2("theoretical2D_func", fullTheoreticalFunc, 0, histMaxTime, 0, 3300, 0);
    theoretical2D_func->SetNpx(nBins);
    theoretical2D_func->SetNpy(3300);

    auto full_theoretical_canvas = new TCanvas("full_theoretical_canvas","full_theoretical_canvas",200,10,1200,800);
    theoretical2D_func->Draw("COLZ");

      // projections of 2D function

      auto timeProjection = new TH1F("timeProjection", "timeProjection", nBins, 0, histMaxTime);
      for (int xBin = 1; xBin <= nBins; ++xBin)
      {
        double timeSum = 0;
        for (int energyBin = 1; energyBin <= 3300; ++energyBin) // energy bin and value are the same
        {
          timeSum += theoretical2D_func->Eval(timeProjection->GetXaxis()->GetBinCenter(xBin), energyBin);
        }
        timeProjection->SetBinContent(xBin, timeSum);
      }
      auto timeProjection_canvas = new TCanvas("timeProjection_canvas","timeProjection_canvas",200,10,1200,800);
      timeProjection->Draw();


      auto energyProjection = new TH1F("energyProjection", "energyProjection", 3300, 0, 3300);
      for (int energyBin = 1; energyBin <= 3300; ++energyBin)
      {
        double energySum = 0;
        for (int xBin = 1; xBin <= nBins; ++xBin)
        {
          energySum += theoretical2D_func->Eval(timeProjection->GetXaxis()->GetBinCenter(xBin), energyBin);
        }
        energyProjection->SetBinContent(energyBin, energySum);
      }
      auto energyProjection_canvas = new TCanvas("energyProjection_canvas","energyProjection_canvas",200,10,1200,800);
      energyProjection->Draw(); // energy plot has a kink in it due to the acceptance, which is smeared out when convoluted with the resolution function

/////////////////////////////////////////////////////////////////////////////////////

    // could potentially plot 2D function of full convolution function, but it takes too long, even for the smaller number of bins here
/*
    auto full2Dfunction_small_hist = new TH2F("full2Dfunction_small_hist", "full2Dfunction_small_hist", nBins, 0, histMaxTime, 3290, 10, 3300);
    for (int i = 1; i <= nBins; ++i)
    {
      for (int j = 1; j <= 3290; ++j)
      {
        full2Dfunction_small_hist->SetBinContent(i, j, full2D_convolutedFunc->Eval(full2Dfunction_small_hist->GetXaxis()->GetBinCenter(i), full2Dfunction_small_hist->GetYaxis()->GetBinCenter(j)));
      }
    }

    auto full2D_func_hist_canv = new TCanvas("full2D_func_hist_canv","full2D_func_hist_canv",200,10,1200,800);
    full2Dfunction_small_hist->Draw("COLZ");
*/
/////////////////////////////////////////////////////////////////////////////////////

  return;
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

int toyCreateTimeEnergyTH2(int timeBinStart)
{

  TRandom3* randGenerator = new TRandom3(0); // random number generators not currently being used - careful with seeds - figure out whether they need to be set or not when used again
  gRandom->SetSeed(0);

  TFile* outputFile = new TFile(Form("timeEnergyHist%i.root",timeBinStart),"RECREATE");

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  // acceptance function

  Acc_func = new TF1("Acc_func", "(x<0)*0 + (x>=0 && x<2350)*(.75 / 2350.)*x + (x>=2350 && x<=3100)*.75 + (x>3100)*0", 0, 3300);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  // define the resolution function

  R_func = new TF1("R_func", resolutionExpModGaus, -4000, 4000, 2);
  R_func->SetNpx(8000);
  // p[0] is the true energy E
  // p[1] is the mode value to shift the exp modified gaussian over - so the peak can be recentered at 0

      R_func->SetParameter(1,0); // make sure resolution functions are unshifted to start in order to find the shifts
      for (int E = 10; E < 4000; ++E)
      {
        R_func->SetParameter(0,E);
        double peak_Xvalue = R_func->GetMaximumX(-33,0); // modes lie within this range, shifted slightly from 0

        resolutionModes[E-10] = peak_Xvalue;
      }

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  // define N and N multiplied by the acceptance

  N_theoretical = new TF1("N_theoretical", "(x<0)*0 + (x>=0 && x<=3100)*((x/3100 - 1)*(4*(x/3100)*(x/3100)-5*(x/3100)-5))/[0] + (x>3100)*0", 0, 3300);
  N_theoretical->SetNpx(3300);
  N_theoretical->SetParameter(0,1);
  double norm_N_th = N_theoretical->Integral(0, 3300);
  N_theoretical->SetParameter(0, norm_N_th);

  N_with_Ac = new TF1("N_with_Ac", N_times_Ac, 0, 3300, 0);
  N_with_Ac->SetNpx(3300);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  // define the Asymmetry

  A_theoretical = new TF1("A_theoretical", "(x<=0)*0 + (x>=0 && x<=3100)*(-8*(x/3100)*(x/3100)+(x/3100)+1) / (4*(x/3100)*(x/3100)-5*(x/3100)-5) + (x>3100)*0", 0, 3300);
  A_theoretical->SetNpx(3300);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  // define drift time with energy dependence

  tdrift_energy_dependent = new TF1("tdrift_energy_dependent", tdriftFunc, 0, 3300, 0); // defined as 0 above 3100 MeV
  tdrift_energy_dependent->SetNpx(3300);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  // define the full convolution function

  full2D_convolutedFunc = new TF2("full2D_convolutedFunc", fullFunctionConv, 0, histMaxTime, 0, 3300, 0); // can replace the function with fullTheoreticalFunc to create the 2D histogram with the perfect theoretical function (no resolution convolution)
  full2D_convolutedFunc->SetNpx(nBins);
  full2D_convolutedFunc->SetNpy(3300);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  // If you want to see plots of the individual functions comment this in.

  // checkCanvases();
  // return 1;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  // construct full 2D histogram by building it up in pieces (columns) and then hadding at the end

  int nTimeBinsTotal = nBins*8; // ~18.64 ns bins - probably good enough segmentation for random sampling later -> goes to 18.65 ns exactly
  int numberOfTimeBinsPerJob = nTimeBinsTotal/numberOfJobs;

  int nEnergyBinsTotal = (3290*2); // .5 MeV bins

  auto full2Dfunction_hist = new TH2F("full2Dfunction_hist", "full2Dfunction_hist", nTimeBinsTotal, 0, histMaxTime, nEnergyBinsTotal, 10, 3300);

  for (int i = 1+timeBinStart*numberOfTimeBinsPerJob; (i <= nTimeBinsTotal && i < (1 + timeBinStart*numberOfTimeBinsPerJob + numberOfTimeBinsPerJob)); ++i)
  {
    for (int j = 1; j <= nEnergyBinsTotal; ++j)
    {
      full2Dfunction_hist->SetBinContent(i, j, full2D_convolutedFunc->Eval(full2Dfunction_hist->GetXaxis()->GetBinCenter(i), full2Dfunction_hist->GetYaxis()->GetBinCenter(j)));
    }
  }

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

   outputFile->Write();
   delete outputFile;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


  return 1;

}
