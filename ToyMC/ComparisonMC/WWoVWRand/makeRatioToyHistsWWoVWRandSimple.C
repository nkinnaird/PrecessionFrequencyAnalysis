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

#include "../../../ratioMacroHeaders/ratioAnalysisDefs.hh" // include paths for potential future grid jobs
#include "../../../ratioMacroHeaders/ratioAnalysisConfig.hh"

using namespace std;

class TruthFuncClass{

public:
  TruthFuncClass(int numTotalParameters)
  {
    numFunctionParameters = numTotalParameters;
  }

  double Evaluate(double* x, double* par){
      double time = x[0];
      double freq = blindingWa * (1 + par[3] * 1e-6);  // want to use frequency corresponding to blinding software since ToyMC fits are unblinded
      double lifetime = par[1] * 1000.; // convert from us to ns for fit function

      double value = par[0] * exp(-time/lifetime) * (1 + par[2] * cos(freq*time + par[4]));
      return value;
  };

private:
  int numFunctionParameters;

}; // end TruthFuncClass


int makeRatioToyHistsWWoVWRandSimple()
{
  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  gRandom->SetSeed(0); // used in GetRandom on TF1
  TRandom3* randGen_positrons = new TRandom3(0); // for number of entries
  TRandom3* randGen_FR = new TRandom3(0);
  TRandom3* randGen_VW = new TRandom3(0);
  TRandom3* randGen_UV = new TRandom3(0);

  double nPts = 1.5e9; // 1.5e9 for about the same precision as the 60 hr dataset - can increase the stats or run more jobs on the grid and hadd them

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

  double Ntoy = nPts / ( truthFunc->Integral(0, histMaxTime) / defaultBinWidth );
  truthFunc->SetParameter(0, Ntoy);

  truthFunc->Write();

/////////////////////////////////////////////////////////////////////////////////////

  double nVal = 0.108;
  double fc = 1/defaultBinWidth;
  double const_fVW = (1-2*sqrt(nVal)) * fc;
  double VWperiod = 1./const_fVW;

  // cout << "VW freq is: " << const_fVW << " period: " << VWperiod << endl;

  // exit(1);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  int totalIters = 200; // number of random seeds basically

  TVectorD numPPMIters(1);
  numPPMIters[0] = totalIters;
  numPPMIters.Write("Iters");

  TDirectory* histIterDirs[totalIters];

  TH1F* toyFiveParamHists[totalIters];
  TH1F* toyUHists[totalIters];
  TH1F* toyVHists[totalIters];

  TH1F* toyFiveParamHists_withVWRand[totalIters];
  TH1F* toyUHists_withVWRand[totalIters];
  TH1F* toyVHists_withVWRand[totalIters];

/////////////////////////////////////////////////////////////////////////////////////

  double halfPeriodGuess = (1/blindingFa)/2;

  double totalChance = exp(halfPeriodGuess/defaultLifetime) + exp(-halfPeriodGuess/defaultLifetime) + 2.;
  double percentChanceUPlus = exp(halfPeriodGuess/defaultLifetime) / totalChance;
  double percentChanceUMinus = exp(-halfPeriodGuess/defaultLifetime) / totalChance;


/////////////////////////////////////////////////////////////////////////////////////

for (int iter = 0; iter < totalIters; ++iter) // loop over totalIters
{
  histIterDirs[iter] = topDir->mkdir(Form("Iter%d",iter));
  histIterDirs[iter]->cd();

  toyFiveParamHists[iter] = new TH1F("Toy_5_Param_Hist","Toy_5_Param_Hist",nBins,0,histMaxTime);
  toyUHists[iter] = new TH1F("Toy_U_Hist","Toy_U_Hist",nBins,0,histMaxTime);
  toyVHists[iter] = new TH1F("Toy_V_Hist","Toy_V_Hist",nBins,0,histMaxTime);

  toyFiveParamHists_withVWRand[iter] = new TH1F("Toy_5_Param_Hist_withVWRand","Toy_5_Param_Hist_withVWRand",nBins,0,histMaxTime);
  toyUHists_withVWRand[iter] = new TH1F("Toy_U_Hist_withVWRand","Toy_U_Hist_withVWRand",nBins,0,histMaxTime);
  toyVHists_withVWRand[iter] = new TH1F("Toy_V_Hist_withVWRand","Toy_V_Hist_withVWRand",nBins,0,histMaxTime);
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
                double time = time_original;
                double timeFR = defaultBinWidth * (randGen_FR->Uniform() - 0.5);
                double timeVW = VWperiod * (randGen_VW->Uniform() - 0.5);

                toyFiveParamHists[iter]->Fill(time + timeFR);
                toyFiveParamHists_withVWRand[iter]->Fill(time + timeFR + timeVW);

                double randNum = randGen_UV->Uniform();

                if     (randNum < percentChanceUPlus){
                  toyUHists[iter]->Fill(time + timeFR - halfPeriodGuess); // careful with the signs here - the U+ hist is filled with pulses shifted by t - T/2, and weighted by e^+T/2tau, and vice versa for U-
                  toyUHists_withVWRand[iter]->Fill(time + timeFR + timeVW - halfPeriodGuess);
                } 
                else if(randNum < (percentChanceUPlus + percentChanceUMinus)){
                  toyUHists[iter]->Fill(time + timeFR + halfPeriodGuess);
                  toyUHists_withVWRand[iter]->Fill(time + timeFR + timeVW + halfPeriodGuess);
                } 
                else if(randNum < 1){
                  toyVHists[iter]->Fill(time + timeFR);
                  toyVHists_withVWRand[iter]->Fill(time + timeFR + timeVW);
                } 

              }
        }

/////////////////////////////////////////////////////////////////////////////////////

  outputFile->Write();
  delete outputFile;

  return 1;

}
