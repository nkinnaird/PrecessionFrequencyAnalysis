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

#include "../ratioMacroHeaders/ratioAnalysisDefs.hh"

using namespace std;

double gainSagFunc(double* x, double* p)
{
  double time = x[0];

  double alpha = 4.2e-7; // some constant
  double n0 = 20.; // number of pulses per fill - 500 CTAGs / 24 calos
  double energy = 20000 * (p[0]); // not really sure about this energy constant parameter

  double tauDrop = 1000.; // intitial gain drop time constant in ns
  double tauRelax = p[1]; // defaultLifetime; // this was tau mu for what Nandita had in her gain paper, but most recently with the mega boxes I think this is something different

  // double gainValue = 1 - alpha * energy * n0 * (tauDrop/(tauRelax-tauDrop)) * (exp(-time/tauRelax) - exp(-time/tauDrop));
  double gainValue = 1 - alpha * energy * n0 * 0.2 * (exp(-time/tauRelax) - exp(-time/tauDrop)); // just put in a constant (.2) instead of the tau combination term

  return gainValue;
}

class gainCorrClass {

public:
  gainCorrClass(TF1* gainSagFunc){
    sagFunc = gainSagFunc;
  }

  double gainCorrFunc(double* x, double* par){
    return 1./sagFunc->Eval(x[0]);
  }

private:
  TF1* sagFunc;

}; // end gainCorrClass class


/////////////////////////////////////////////////////////////////////////////////////

int makeToyHistFromTH2(string filePath)
{
  TRandom3* randGenerator = new TRandom3(444);
  gRandom->SetSeed(555);

  double nPts = 1e9;

  // pull in input file
  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  TFile* outputFile = new TFile("histFromPairs.root","RECREATE");

/////////////////////////////////////////////////////////////////////////////////////

  TH2F* initialHist = (TH2F*) inputFile->Get("full2Dfunction_hist");

/////////////////////////////////////////////////////////////////////////////////////

  auto topDir = gFile->mkdir("topDir");

/////////////////////////////////////////////////////////////////////////////////////

  int totalIters = 1;

  TVectorD numIters(1);
  numIters[0] = totalIters;
  numIters.Write("Iters");

  TDirectory* histIterDirs[totalIters];

  TH1F* toyFiveParamHists[totalIters];
  TH1F* toyUHists[totalIters];
  TH1F* toyVHists[totalIters];
  TH1F* toyNumHists[totalIters];
  TH1F* toyDenomHists[totalIters];
  TH1F* energies[totalIters];

/////////////////////////////////////////////////////////////////////////////////////

  TF1* gainSagFunctions[totalIters];
  TF1* gainCorrFunctions[totalIters];

/////////////////////////////////////////////////////////////////////////////////////

  double gainSagAmplitudeMultiplier[totalIters]; 
  double gainSagLifetime[totalIters];

  double gm2PeriodPPMOffsets[totalIters];
  double gm2PeriodGuesses[totalIters];
  double energyThresholds[totalIters];

  for (int iter = 0; iter < totalIters; ++iter)
  {
    // amplitude
    gainSagAmplitudeMultiplier[iter] = 1.;// - 1.*iter/totalIters; // realistic gain sags
    gainSagLifetime[iter] = 10000;
    // gainSagAmplitudeMultiplier[iter] = 1.*(iter-(totalIters/2.))/totalIters; // scan through 0 gain sag (negative to positive effect)
    // gainSagLifetime[iter] = 10000;

    // lifetime 
    // gainSagAmplitudeMultiplier[iter] = 1;
    // gainSagLifetime[iter] = 5000 + 2000*iter;


    gm2PeriodPPMOffsets[iter] = 0;
    gm2PeriodGuesses[iter] = g2Period * (1 + 1e-6 * gm2PeriodPPMOffsets[iter]);
    energyThresholds[iter] = defaultEThreshold;
  }

  gainSagAmplitudeMultiplier[0] = 0; // no gain sag for iteration 0
  gainSagLifetime[0] = 10000; // don't think it should matter what this number is exactly at long as the amplitude is 0

/////////////////////////////////////////////////////////////////////////////////////

  for (int iter = 0; iter < totalIters; ++iter)
  {
    histIterDirs[iter] = topDir->mkdir(Form("Iter%d", iter));
    histIterDirs[iter]->cd();

    toyFiveParamHists[iter] = new TH1F("Toy_5_Param_Hist","Toy_5_Param_Hist",nBins,0,histMaxTime);
    toyUHists[iter] = new TH1F("Toy_U_Hist","Toy_U_Hist",nBins,0,histMaxTime);
    toyVHists[iter] = new TH1F("Toy_V_Hist","Toy_V_Hist",nBins,0,histMaxTime);
    toyNumHists[iter] = new TH1F("Toy_Num_Hist","Toy_Num_Hist",nBins,0,histMaxTime);
    toyDenomHists[iter] = new TH1F("Toy_Denom_Hist","Toy_Denom_Hist",nBins,0,histMaxTime);

    energies[iter] = new TH1F("energies","energies",3100,0,3100);

/////////////////////////////////////////////////////////////////////////////////////

    auto savDirs = histIterDirs[iter]->mkdir("SavedParameters");
    savDirs->cd();

    gainSagFunctions[iter] = new TF1("gainSagFunction", gainSagFunc, 0, 1000000, 2);
    gainSagFunctions[iter]->SetNpx(10000);
    gainSagFunctions[iter]->SetParameter(0, gainSagAmplitudeMultiplier[iter]);
    gainSagFunctions[iter]->SetParameter(1, gainSagLifetime[iter]);
    gainSagFunctions[iter]->Write();

    auto corrClass = new gainCorrClass(gainSagFunctions[iter]); // if grabbing gainSagFunctions[0] (the no gain sag curve) then there is no correction performed
    gainCorrFunctions[iter] = new TF1("gainCorrFunction", corrClass, &gainCorrClass::gainCorrFunc, 0, 1000000, 0);
    gainCorrFunctions[iter]->SetNpx(10000);
    gainCorrFunctions[iter]->Write();


    TVectorD parameterStore(2); // vector of doubles to store parameters used when creating histograms, for systematic studies
    // energy threshold, g-2 period used, possibly function parameters, etc.
    // should this just be a tntuple? or maybe a tlist with a variety of root objects? for tntuple see the ratioMacro code

    parameterStore[0] = gm2PeriodGuesses[iter];
    parameterStore[1] = energyThresholds[iter];
    parameterStore.Write("parameterStore");
  }

/////////////////////////////////////////////////////////////////////////////////////

  double randEntries = randGenerator->Poisson(nPts);
  int tenIncrement = 0.1 * randEntries;
  int countdown = tenIncrement;
  int tenPercent = 0;

        for (int entry = 0; entry < randEntries; ++entry)
        {
          if(--countdown == 0) // progress output
          {
            tenPercent++;
            cout << tenPercent << "0%" << " completed" << endl;
            countdown = tenIncrement;
          }


          double time;
          double energy_raw;

          initialHist->GetRandom2(time, energy_raw);

          double randNum = randGenerator->Uniform(); // generate this random number above the iterations loop so histograms are filled consistently
          // could possibly just use multiple random number generators with the same seed that are only used when filling UV hists as it previously was

          for (int iter = 0; iter < totalIters; ++iter)
          {
              double gainSagMultiple = gainSagFunctions[iter]->Eval(time);
              double gainCorrMultiple = gainCorrFunctions[iter]->Eval(time);

              double gausSmearSag = randGenerator->Gaus(gainSagMultiple, 0.01*gainSagMultiple);

              double energy_measured;

              if(iter == 0) energy_measured = energy_raw; // don't modify the 0th iteration energy
              else energy_measured = energy_raw * gausSmearSag * gainCorrMultiple;

              energies[iter]->Fill(energy_measured);

              if(energy_measured > energyThresholds[iter])
              {
                    toyFiveParamHists[iter]->Fill(time);

                    double halfPeriodGuess = gm2PeriodGuesses[iter]/2;

                    double totalChance = exp(halfPeriodGuess/defaultLifetime) + exp(-halfPeriodGuess/defaultLifetime) + 2.;
                    double percentChanceUPlus = exp(halfPeriodGuess/defaultLifetime) / totalChance;
                    double percentChanceUMinus = exp(-halfPeriodGuess/defaultLifetime) / totalChance;

                    if     (randNum < percentChanceUPlus) toyUHists[iter]->Fill(time - halfPeriodGuess); // careful with the signs here - the U+ hist is filled with pulses shifted by t - T/2, and weighted by e^+T/2tau, and vice versa for U-
                    else if(randNum < (percentChanceUPlus + percentChanceUMinus)) toyUHists[iter]->Fill(time + halfPeriodGuess);
                    else if(randNum < 1) toyVHists[iter]->Fill(time);
              }

          }

        }

        for (int iter = 0; iter < totalIters; ++iter)
        {
          toyNumHists[iter]->Add(toyUHists[iter], toyVHists[iter], -1, 1);
          toyDenomHists[iter]->Add(toyUHists[iter], toyVHists[iter]);
        }

/////////////////////////////////////////////////////////////////////////////////////

  // this is just a place holder function so that I can use the same analysis code which gets this function for making truth pulls
  // the parameters don't have any direct meaning, but I can put in values that correspond to the original function for some
  // the R and lifetime parameters should still be able to be used in the pulls, the default phase should be slightly off, N I can just use the histogram entries I think (it will be close at least), and A I'm not sure about
  auto truthFunc = new TF1("truthFunc", "[0] + [1] + [2] + [3] + [4]", 0, histMaxTime);
  truthFunc->SetNpx(100);

  truthFunc->SetParameter(0, toyFiveParamHists[0]->GetEntries()*binWidth/defaultLifetime);
  truthFunc->SetParameter(1, defaultLifetime);
  truthFunc->SetParameter(2, 1); // random value
  truthFunc->SetParameter(3, 0);
  truthFunc->SetParameter(4, startingPhase);

  outputFile->cd();
  truthFunc->Write();

/////////////////////////////////////////////////////////////////////////////////////
  outputFile->Write();

  delete outputFile;


return 1;

}
