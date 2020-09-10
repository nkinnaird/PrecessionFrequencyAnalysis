#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
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

#include "../ratioMacroHeaders/ratioAnalysisDefs.hh" // include paths for potential future grid jobs
#include "../ratioMacroHeaders/ratioAnalysisConfig.hh"

using namespace std;

class TruthFuncClass{

public:
  TruthFuncClass(int numTotalParameters)
  {
    numFunctionParameters = numTotalParameters;
  }

  double NCBOPart(double t, double* par){
    // par[5] = cbo frequency
    // par[6] = cbo lifetime
    // par[7] = N cbo amplitude
    // par[8] = N cbo phase

    double N_cbo;

    if(numFunctionParameters > 5) N_cbo = par[0] * (1 + exp(-t/par[6]) * par[7] * cos(par[5]*t + par[8]));
    else N_cbo = par[0]; // no cbo

    return N_cbo;
  };

  double ACBOPart(double t, double* par){
    // par[9] = A cbo amplitude
    // par[10] = A cbo phase

    double A_cbo;

    if(numFunctionParameters > 9) A_cbo = par[2] * (1 + exp(-t/par[6]) * par[9] * cos(par[5]*t + par[10]));
    else A_cbo = par[2]; // no cbo

    return A_cbo;
  };

  double PCBOPart(double t, double* par){
    // par[11] = Phase cbo amplitude
    // par[12] = Phase cbo phase

    double P_cbo;

    if(numFunctionParameters > 11) P_cbo = par[4] * (1 + exp(-t/par[6]) * par[11] * cos(par[5]*t + par[12]));
    else P_cbo = par[4]; // no cbo

    return P_cbo;
  };

  double Evaluate(double* x, double* par){
      double time = x[0];
      // double freq = defaultWa * (1 + par[3] * 1e-6);
      double freq = blindingWa * (1 + par[3] * 1e-6);  // want to use frequency corresponding to blinding software since ToyMC fits are unblinded
      double lifetime = par[1] * 1000.; // convert from us to ns for fit function

      // double value = NCBOPart(time, par) * exp(-time/lifetime) * (1 + ACBOPart(time, par) * cos(freq*time + PCBOPart(time, par)));
      double value = par[0] * exp(-time/lifetime) * (1 + par[2] * cos(freq*time + par[4]));
      return value;
  };

private:
  int numFunctionParameters;

}; // end TruthFuncClass


int makeRatioToyHists()
{
  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  gRandom->SetSeed(0); // used in GetRandom on TF1
  TRandom3* randGen_positrons = new TRandom3(0); // for number of entries
  TRandom3* randGen_FR = new TRandom3(0);
  TRandom3* randGen_UV = new TRandom3(0);

  double nPts = 1.5e5; // 1.5e9 for about the same precision as the 60 hr dataset - can increase the stats or run more jobs on the grid and hadd them

  // create output file that will hold plots
  TFile* outputFile = new TFile("ratioToyHists.root","RECREATE");

  // make top directory for output file
  auto topDir = gFile->mkdir("topDir");

/////////////////////////////////////////////////////////////////////////////////////

  int nBins = int(approxMaxTime/defaultBinWidth);
  double histMaxTime = nBins*defaultBinWidth;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  setDataset(1, 12); // set dataset case to 60h - for 'starting' parameters

/////////////////////////////////////////////////////////////////////////////////////

  int numParameters = 5; // 5, 9, 11, 13

  auto myTruthFuncClass = new TruthFuncClass(numParameters);
  auto truthFunc = new TF1("truthFunc", myTruthFuncClass, &TruthFuncClass::Evaluate, 0, histMaxTime, numParameters);
  truthFunc->SetNpx(10000);

  truthFunc->SetParameter(0, 1);
  truthFunc->SetParameter(1, defaultLifetime/1000.);
  truthFunc->SetParameter(2, startingAsymmetry);
  truthFunc->SetParameter(3, 0);
  truthFunc->SetParameter(4, startingPhase);

  // currently envelope uses exponential
  if(numParameters > 5){
    cout << "Have to update cbo frequency stuff in ToyMC before using it again." << endl;
    exit(1);
    truthFunc->SetParameter(5, 0);//startingCBOFreq); // cbo frequency
    truthFunc->SetParameter(6, 250000); // cbo lifetime

    truthFunc->SetParameter(7, 0.005); // N_cbo amplitude
    truthFunc->SetParameter(8, 0); // N_cbo phase
  }

  if(numParameters > 9){
    truthFunc->SetParameter(9, 0.005); // A_cbo amplitude
    truthFunc->SetParameter(10, pi/3); // A_cbo phase
  }

  if(numParameters > 11){
    truthFunc->SetParameter(11, 0.005); // Phi_cbo amplitude
    truthFunc->SetParameter(12, pi/2); // Phi_cbo phase
  }

  double Ntoy = nPts / ( truthFunc->Integral(0, histMaxTime) / defaultBinWidth );
  truthFunc->SetParameter(0, Ntoy);

  truthFunc->Write();

/////////////////////////////////////////////////////////////////////////////////////

  double nVal = 0.120;
  double fc = 1/defaultBinWidth;
  double const_fVW = (1-2*sqrt(nVal)) * fc;
  double VWperiod = 1./const_fVW;

  cout << "VW freq is: " << const_fVW << " period: " << VWperiod << endl;

  // exit(1);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  int totalIters = 1;

  TVectorD numPPMIters(1);
  numPPMIters[0] = totalIters;
  numPPMIters.Write("Iters");

  TDirectory* histIterDirs[totalIters];

  TH1F* toyFiveParamHists[totalIters];
  TH1F* toyUHists[totalIters];
  TH1F* toyVHists[totalIters];
  TH1F* toyNumHists[totalIters];
  TH1F* toyDenomHists[totalIters];

/////////////////////////////////////////////////////////////////////////////////////

  double gm2PeriodPPMOffsets[totalIters];
  double gm2PeriodGuesses[totalIters];

/////////////////////////////////////////////////////////////////////////////////////

for (int iter = 0; iter < totalIters; ++iter) // loop over totalIters
{
  gm2PeriodPPMOffsets[iter] = 0;
  // gm2PeriodGuesses[iter] = g2Period * (1 + 1e-6 * gm2PeriodPPMOffsets[iter]);
  gm2PeriodGuesses[iter] = 1/blindingFa * (1 + 1e-6 * gm2PeriodPPMOffsets[iter]); // want to use g2Period corresponding to blinding software since ToyMC fits are unblinded


  histIterDirs[iter] = topDir->mkdir(Form("Iter%d",iter));
  histIterDirs[iter]->cd();

  toyFiveParamHists[iter] = new TH1F("Toy_5_Param_Hist","Toy_5_Param_Hist",nBins,0,histMaxTime);
  toyUHists[iter] = new TH1F("Toy_U_Hist","Toy_U_Hist",nBins,0,histMaxTime);
  toyVHists[iter] = new TH1F("Toy_V_Hist","Toy_V_Hist",nBins,0,histMaxTime);
  // toyNumHists[iter] = new TH1F("Toy_Num_Hist","Toy_Num_Hist",nBins,0,histMaxTime);
  // toyDenomHists[iter] = new TH1F("Toy_Denom_Hist","Toy_Denom_Hist",nBins,0,histMaxTime);

  auto savDirs = histIterDirs[iter]->mkdir("SavedParameters");
  savDirs->cd();

  TVectorD parameterStore(2); // vector of doubles to store parameters used when creating histograms, for systematic studies
  // energy threshold, g-2 period used, possibly function parameters, etc.
  // should this just be a tntuple? or maybe a tlist with a variety of root objects? for tntuple see the ratioMacro code

  parameterStore[0] = gm2PeriodGuesses[iter];
  parameterStore[1] = 0; // supposed to be the energy threshold but that doesn't apply for this toy mc version
  parameterStore.Write("parameterStore");
}

    double randEntries = randGen_positrons->PoissonD(nPts);
    double tenIncrement = 0.1 * randEntries;
    double countdown = tenIncrement;
    int tenPercent = 0;

        for (double entry = 0; entry < randEntries; ++entry)
        {

          if(--countdown <= 0) // progress output
          {
            tenPercent++;
            cout << tenPercent << "0%" << " completed" << endl;
            countdown = tenIncrement;
          }

          double time_original = truthFunc->GetRandom();
         
              for (int iter = 0; iter < totalIters; ++iter)
              {
                double time = time_original + defaultBinWidth * (randGen_FR->Uniform() - 0.5) + VWperiod * (randGen_FR->Uniform() - 0.5);

                toyFiveParamHists[iter]->Fill(time);

                double randNum = randGen_UV->Uniform();

                double halfPeriodGuess = gm2PeriodGuesses[iter]/2;

                double totalChance = exp(halfPeriodGuess/defaultLifetime) + exp(-halfPeriodGuess/defaultLifetime) + 2.;
                double percentChanceUPlus = exp(halfPeriodGuess/defaultLifetime) / totalChance;
                double percentChanceUMinus = exp(-halfPeriodGuess/defaultLifetime) / totalChance;

                if     (randNum < percentChanceUPlus) toyUHists[iter]->Fill(time - halfPeriodGuess); // careful with the signs here - the U+ hist is filled with pulses shifted by t - T/2, and weighted by e^+T/2tau, and vice versa for U-
                else if(randNum < (percentChanceUPlus + percentChanceUMinus)) toyUHists[iter]->Fill(time + halfPeriodGuess);
                else if(randNum < 1) toyVHists[iter]->Fill(time);

              }
        }

        // for (int iter = 0; iter < totalIters; ++iter)
        // {
        //   toyNumHists[iter]->Add(toyUHists[iter], toyVHists[iter], -1, 1);
        //   toyDenomHists[iter]->Add(toyUHists[iter], toyVHists[iter]);
        // }

/////////////////////////////////////////////////////////////////////////////////////

  outputFile->Write();

  delete outputFile;

  return 1;

}
