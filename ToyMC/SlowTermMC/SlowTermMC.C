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

#include "ratioAnalysisDefs.hh" // include paths if you want to run this on the grid
#include "ratioAnalysisConfig.hh"

using namespace std;

class TruthFuncClass{

public:
  TruthFuncClass()
  {}

  double slowPart(double time, double* par){
    double value = 0;
    double freq = blindingWa * (1 + par[3] * 1e-6);  // want to use frequency corresponding to blinding software since ToyMC fits are unblinded
    double lifetime = par[1] * 1000.; // convert from us to ns for fit function

    // value = -1e-3 * exp(-time/lifetime);
    // value = -1e-3 * exp(-time/lifetime) * (1 + par[2] * cos(freq*time + par[4]));
    value = -1e-3 * exp(-time/lifetime) * (1 + 0.2 * cos(freq*time + par[4]));

    return value;
  }

  double SlowTerm(double* x, double* par){
    return (1 + slowPart(x[0], par));
  }

  double Evaluate(double* x, double* par){
      double time = x[0];
      double freq = blindingWa * (1 + par[3] * 1e-6);  // want to use frequency corresponding to blinding software since ToyMC fits are unblinded
      double lifetime = par[1] * 1000.; // convert from us to ns for fit function

      double value = (1 + slowPart(time, par)) * par[0] * exp(-time/lifetime) * (1 + par[2] * cos(freq*time + par[4]));
      return value;
  };

}; // end TruthFuncClass


int SlowTermMC()
{
  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  gRandom->SetSeed(987); // used in GetRandom on TF1
  TRandom3* randGen_positrons = new TRandom3(444); // for number of entries
  TRandom3* randGen_UV = new TRandom3(123);

  double nPts = 1.5e9; // 1.5e9 for about the same precision as the 60 hr dataset - can increase the stats or run more jobs on the grid and hadd them

  // create output file that will hold plots
  TFile* outputFile = new TFile("slowMCHists.root","RECREATE");

/////////////////////////////////////////////////////////////////////////////////////

  int nBins = int(approxMaxTime/defaultBinWidth);
  double histMaxTime = nBins*defaultBinWidth;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  // setDataset(1, 12); // set dataset case to 60h - for 'starting' parameters - will print out some other things which won't apply to this Toy MC

/////////////////////////////////////////////////////////////////////////////////////




  auto myTruthFuncClass = new TruthFuncClass();
  auto truthFunc = new TF1("truthFunc", myTruthFuncClass, &TruthFuncClass::Evaluate, 0, histMaxTime, 5);
  truthFunc->SetNpx(10000);

  truthFunc->SetParameter(0, 1);
  truthFunc->SetParameter(1, defaultLifetime/1000.);
  truthFunc->SetParameter(2, 0.3688); // startingAsymmetry);
  truthFunc->SetParameter(3, 0);
  truthFunc->SetParameter(4, 2.076); // startingPhase);


  double Ntoy = nPts / ( truthFunc->Integral(0, histMaxTime) / defaultBinWidth );
  truthFunc->SetParameter(0, Ntoy);

  truthFunc->Write();



  auto justSlowTerm = new TF1("justSlowTerm", myTruthFuncClass, &TruthFuncClass::SlowTerm, 0, histMaxTime, 5);
  justSlowTerm->SetNpx(10000);

  justSlowTerm->SetParameter(0, 1);
  justSlowTerm->SetParameter(1, defaultLifetime/1000.);
  justSlowTerm->SetParameter(2, 0.3688); // startingAsymmetry);
  justSlowTerm->SetParameter(3, 0);
  justSlowTerm->SetParameter(4, 2.076); // startingPhase);

  justSlowTerm->Write();

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

          double halfPeriod = (1/blindingFa)/2;
          double totalChance = exp(halfPeriod/defaultLifetime) + exp(-halfPeriod/defaultLifetime) + 2.;
          double percentChanceUPlus = exp(halfPeriod/defaultLifetime) / totalChance;
          double percentChanceUMinus = exp(-halfPeriod/defaultLifetime) / totalChance;

/////////////////////////////////////////////////////////////////////////////////////


  TH1F* toyFiveParamHist = new TH1F("Toy_5_Param_Hist","Toy_5_Param_Hist",nBins,0,histMaxTime);
  TH1F* toyUHist = new TH1F("Toy_U_Hist","Toy_U_Hist",nBins,0,histMaxTime);
  TH1F* toyVHist = new TH1F("Toy_V_Hist","Toy_V_Hist",nBins,0,histMaxTime);


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
         
          double time = truthFunc->GetRandom();

          toyFiveParamHist->Fill(time);

          double randNum = randGen_UV->Uniform();

          if     (randNum < percentChanceUPlus) toyUHist->Fill(time - halfPeriod);
          else if(randNum < (percentChanceUPlus + percentChanceUMinus)) toyUHist->Fill(time + halfPeriod);
          else if(randNum < 1) toyVHist->Fill(time);

        }

/////////////////////////////////////////////////////////////////////////////////////

  outputFile->Write();

  delete outputFile;

  return 1;

}
