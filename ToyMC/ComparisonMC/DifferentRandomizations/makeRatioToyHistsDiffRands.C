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


int makeRatioToyHistsDiffRands()
{
  clock_t t; // for time profiling
  t = clock();

  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  gRandom->SetSeed(0); // used in GetRandom on TF1
  TRandom3* randGen_positrons = new TRandom3(0); // for number of entries

  int randSeedStart = randGen_positrons->Integer(2147483647);

  TRandom3* randGen_1 = new TRandom3(randSeedStart);
  TRandom3* randGen_2 = new TRandom3(randSeedStart);
  TRandom3* randGen_3 = new TRandom3(randSeedStart);
  TRandom3* randGen_4 = new TRandom3(randSeedStart);
  TRandom3* randGen_5 = new TRandom3(randSeedStart);

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

  // double nVal = 0.120;
  // double nVal = 0.108;
  // double fc = 1/defaultBinWidth;
  // double const_fVW = (1-2*sqrt(nVal)) * fc;
  // double VWperiod = 1./const_fVW;

  // cout << "VW freq is: " << const_fVW << " period: " << VWperiod << endl;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  int totalIters = 100; // number of random seeds basically

  TVectorD numPPMIters(1);
  numPPMIters[0] = totalIters;
  numPPMIters.Write("Iters");

  TDirectory* histIterDirs[totalIters];

  TH1F* UVrands[totalIters];

  TH1F* toyFiveParamHists[totalIters];
  TH1F* toyUHists[totalIters];
  TH1F* toyVHists[totalIters];

  // TH1F* timeShifts_1[totalIters];

  TH1F* toyFiveParamHists_1[totalIters];
  TH1F* toyUHists_1[totalIters];
  TH1F* toyVHists_1[totalIters];

  // TH1F* timeShifts_2[totalIters];

  TH1F* toyFiveParamHists_2[totalIters];
  TH1F* toyUHists_2[totalIters];
  TH1F* toyVHists_2[totalIters];

  // TH1F* timeShifts_3[totalIters];

  TH1F* toyFiveParamHists_3[totalIters];
  TH1F* toyUHists_3[totalIters];
  TH1F* toyVHists_3[totalIters];

  // TH1F* timeShifts_4[totalIters];

  TH1F* toyFiveParamHists_4[totalIters];
  TH1F* toyUHists_4[totalIters];
  TH1F* toyVHists_4[totalIters];

  // TH1F* timeShifts_5[totalIters];

  TH1F* toyFiveParamHists_5[totalIters];
  TH1F* toyUHists_5[totalIters];
  TH1F* toyVHists_5[totalIters];

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

  UVrands[iter] = new TH1F("UVrands", "UV Random Numbers; ; Events", 100, 0, 1);

  toyFiveParamHists[iter] = new TH1F("Toy_5_Param_Hist","Toy_5_Param_Hist",nBins,0,histMaxTime);
  toyUHists[iter] = new TH1F("Toy_U_Hist","Toy_U_Hist",nBins,0,histMaxTime);
  toyVHists[iter] = new TH1F("Toy_V_Hist","Toy_V_Hist",nBins,0,histMaxTime);

  // timeShifts_1[iter] = new TH1F("TimeShifts_1", "Time Shifts 1; Time (ns); Events", 1000, -1000, 1000);

  toyFiveParamHists_1[iter] = new TH1F("Toy_5_Param_Hist_1","Toy_5_Param_Hist_1",nBins,0,histMaxTime);
  toyUHists_1[iter] = new TH1F("Toy_U_Hist_1","Toy_U_Hist_1",nBins,0,histMaxTime);
  toyVHists_1[iter] = new TH1F("Toy_V_Hist_1","Toy_V_Hist_1",nBins,0,histMaxTime);

  // timeShifts_2[iter] = new TH1F("TimeShifts_2", "Time Shifts 2; Time (ns); Events", 1000, -1000, 1000);

  toyFiveParamHists_2[iter] = new TH1F("Toy_5_Param_Hist_2","Toy_5_Param_Hist_2",nBins,0,histMaxTime);
  toyUHists_2[iter] = new TH1F("Toy_U_Hist_2","Toy_U_Hist_2",nBins,0,histMaxTime);
  toyVHists_2[iter] = new TH1F("Toy_V_Hist_2","Toy_V_Hist_2",nBins,0,histMaxTime);

  // timeShifts_3[iter] = new TH1F("TimeShifts_3", "Time Shifts 3; Time (ns); Events", 1000, -1000, 1000);

  toyFiveParamHists_3[iter] = new TH1F("Toy_5_Param_Hist_3","Toy_5_Param_Hist_3",nBins,0,histMaxTime);
  toyUHists_3[iter] = new TH1F("Toy_U_Hist_3","Toy_U_Hist_3",nBins,0,histMaxTime);
  toyVHists_3[iter] = new TH1F("Toy_V_Hist_3","Toy_V_Hist_3",nBins,0,histMaxTime);

  // timeShifts_4[iter] = new TH1F("TimeShifts_4", "Time Shifts 4; Time (ns); Events", 1000, -1000, 1000);

  toyFiveParamHists_4[iter] = new TH1F("Toy_5_Param_Hist_4","Toy_5_Param_Hist_4",nBins,0,histMaxTime);
  toyUHists_4[iter] = new TH1F("Toy_U_Hist_4","Toy_U_Hist_4",nBins,0,histMaxTime);
  toyVHists_4[iter] = new TH1F("Toy_V_Hist_4","Toy_V_Hist_4",nBins,0,histMaxTime);

  // timeShifts_5[iter] = new TH1F("TimeShifts_5", "Time Shifts 5; Time (ns); Events", 1000, -1000, 1000);

  toyFiveParamHists_5[iter] = new TH1F("Toy_5_Param_Hist_5","Toy_5_Param_Hist_5",nBins,0,histMaxTime);
  toyUHists_5[iter] = new TH1F("Toy_U_Hist_5","Toy_U_Hist_5",nBins,0,histMaxTime);
  toyVHists_5[iter] = new TH1F("Toy_V_Hist_5","Toy_V_Hist_5",nBins,0,histMaxTime);
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

                double timeRand_1 = 150. * (randGen_1->Uniform() - 0.5); // I've realized after the fact that this is dumb, calling Uniform 5 times when they all return the same number, fix it before rerunning this again
                double timeRand_2 = 300. * (randGen_2->Uniform() - 0.5);
                double timeRand_3 = 450. * (randGen_3->Uniform() - 0.5);
                double timeRand_4 = 600. * (randGen_4->Uniform() - 0.5);
                double timeRand_5 = 750. * (randGen_5->Uniform() - 0.5);

                // timeShifts_1[iter]->Fill(timeRand_1);
                // timeShifts_2[iter]->Fill(timeRand_2);
                // timeShifts_3[iter]->Fill(timeRand_3);
                // timeShifts_4[iter]->Fill(timeRand_4);
                // timeShifts_5[iter]->Fill(timeRand_5);

                toyFiveParamHists[iter]->Fill(time);
                toyFiveParamHists_1[iter]->Fill(time + timeRand_1);
                toyFiveParamHists_2[iter]->Fill(time + timeRand_2);
                toyFiveParamHists_3[iter]->Fill(time + timeRand_3);
                toyFiveParamHists_4[iter]->Fill(time + timeRand_4);
                toyFiveParamHists_5[iter]->Fill(time + timeRand_5);

                double randNum = randGen_UV->Uniform();
                UVrands[iter]->Fill(randNum);

                if     (randNum < percentChanceUPlus){
                  toyUHists[iter]->Fill(time - halfPeriodGuess); // careful with the signs here - the U+ hist is filled with pulses shifted by t - T/2, and weighted by e^+T/2tau, and vice versa for U-
                  toyUHists_1[iter]->Fill(time + timeRand_1 - halfPeriodGuess);
                  toyUHists_2[iter]->Fill(time + timeRand_2 - halfPeriodGuess);
                  toyUHists_3[iter]->Fill(time + timeRand_3 - halfPeriodGuess);
                  toyUHists_4[iter]->Fill(time + timeRand_4 - halfPeriodGuess);
                  toyUHists_5[iter]->Fill(time + timeRand_5 - halfPeriodGuess);
                } 
                else if(randNum < (percentChanceUPlus + percentChanceUMinus)){
                  toyUHists[iter]->Fill(time + halfPeriodGuess);
                  toyUHists_1[iter]->Fill(time + timeRand_1 + halfPeriodGuess);
                  toyUHists_2[iter]->Fill(time + timeRand_2 + halfPeriodGuess);
                  toyUHists_3[iter]->Fill(time + timeRand_3 + halfPeriodGuess);
                  toyUHists_4[iter]->Fill(time + timeRand_4 + halfPeriodGuess);
                  toyUHists_5[iter]->Fill(time + timeRand_5 + halfPeriodGuess);
                } 
                else if(randNum < 1){
                  toyVHists[iter]->Fill(time);
                  toyVHists_1[iter]->Fill(time + timeRand_1);
                  toyVHists_2[iter]->Fill(time + timeRand_2);
                  toyVHists_3[iter]->Fill(time + timeRand_3);
                  toyVHists_4[iter]->Fill(time + timeRand_4);
                  toyVHists_5[iter]->Fill(time + timeRand_5);
                } 

              }
        }

/////////////////////////////////////////////////////////////////////////////////////

  outputFile->Write();
  delete outputFile;

  t = clock() - t;
  printf ("It took me %f seconds.\n", ((float)t)/CLOCKS_PER_SEC);

  return 1;

}
