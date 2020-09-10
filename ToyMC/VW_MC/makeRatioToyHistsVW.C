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

#include "../../ratioMacroHeaders/ratioAnalysisDefs.hh" // include paths for potential future grid jobs
// #include "../../ratioMacroHeaders/ratioAnalysisConfig.hh"
#include "../../ratioMacroHeaders/ratioToyAnalysisConfig.hh"

using namespace std;

class TruthFuncClass{

public:
  TruthFuncClass()
  {}

  double CBOFreqPart(double t, double* par){
    // double frequency = par[5] * (1 + tM_Aconst*exp(-t/tM_Atau)/(tM_w0*t) + tM_Bconst*exp(-t/tM_Btau)/(tM_w0*t));
    double frequency = par[5];

    return (frequency * 1e-3); // convert from rad/us to rad/ns
  };

  double VWPart(double t, double* par){
    double VW_part;
    // double lifetime = par[10] * 1000.;  // convert from us to ns for the fit function
    double frequency;

      // double fcyc = 2*pi/(binWidth); // rad/ns
      // double fcbo = (1 + par[9]) * CBOFreqPart(t, par);
      // frequency = fcyc - 2*sqrt(2*fcyc*fcbo - fcbo*fcbo);
    
      frequency = par[9] * 1e-3; // convert to rad/ns for fit function from rad/us, which is the units the parameter is in
      // frequency = par[9];

    // VW_part = (1 + exp(-t/lifetime) * par[11] * cos(frequency*t + par[12]));
    VW_part = (1 + par[10] * cos(frequency*t + par[11]));
    
    return VW_part;
  };

  double NCBOPart(double t, double* par){
    double lifetime = par[6] * 1000.; // convert from us to ns for the fit function
    double N_cbo = (1 + exp(-t/lifetime) * par[7] * cos(CBOFreqPart(t, par)*t + par[8]));

    return N_cbo;
  };

  double VOPart(double t, double* par){
    double lifetime = 100 * 1000.; // convert from us to ns for the fit function
    double VO_part = (1 + exp(-t/lifetime) * 0.02 * cos(2*pi*2.2*1e-3*t + 1.5));

    return VO_part;
  };

  double Evaluate(double* x, double* par){
      double time = x[0];
      // double freq = defaultWa * (1 + par[3] * 1e-6);
      double freq = blindingWa * (1 + par[3] * 1e-6);  // want to use frequency corresponding to blinding software since ToyMC fits are unblinded
      double lifetime = par[1] * 1000.; // convert from us to ns for fit function

      double value = VWPart(time, par) * NCBOPart(time, par) * par[0] * exp(-time/lifetime) * (1 + par[2] * cos(freq*time + par[4]));
      // double value = VOPart(time, par) * NCBOPart(time, par) * par[0] * exp(-time/lifetime) * (1 + par[2] * cos(freq*time + par[4]));
      return value;
  };

private:


}; // end TruthFuncClass


int makeRatioToyHistsVW()
{
  clock_t t; // for time profiling - put after input
  t = clock();

/////////////////////////////////////////////////////////////////////////////////////

  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  gRandom->SetSeed(0);//444); // used in GetRandom on TF1
  TRandom3* randGen_positrons = new TRandom3(0);//43210);// for number of entries
  TRandom3* randGen_FR = new TRandom3(0);//987);
  TRandom3* randGen_UV = new TRandom3(0);//112233);

  int func_iterRandSeed = gRandom->Integer(1e8);
  int UV_iterRandSeed = randGen_UV->Integer(1e8);


  double nPts = 1.5e7; // 1.5e9 for about the same precision as the 60 hr dataset

  // create output file that will hold plots
  TFile* outputFile = new TFile("ratioToyHists.root","RECREATE");

  // make top directory for output file
  auto topDir = gFile->mkdir("topDir");

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  double functionStart = 10000; // can't be 0 because the frequency model blows up there

  int numParameters = 12; // 5, 9, 11, 13

  auto myTruthFuncClass = new TruthFuncClass();
  auto truthFunc = new TF1("truthFunc", myTruthFuncClass, &TruthFuncClass::Evaluate, functionStart, histMaxTime, numParameters);
  truthFunc->SetNpx(100000);

  truthFunc->SetParameter(0, 1);
  truthFunc->SetParameter(1, defaultLifetime/1000.);
  truthFunc->SetParameter(2, startingAsymmetry);
  truthFunc->SetParameter(3, 0);
  truthFunc->SetParameter(4, startingPhase);

  truthFunc->SetParameter(5, trackingCBOFreq);
  truthFunc->SetParameter(6, startingCBOTau);
  truthFunc->SetParameter(7, startingCBONAmp);
  truthFunc->SetParameter(8, startingCBONPhase);

  // truthFunc->SetParameter(9, startingVWFreq);
  // truthFunc->SetParameter(10, startingVWTau);
  // truthFunc->SetParameter(11, startingVWAmp);
  // truthFunc->SetParameter(12, startingVWPhase);

  truthFunc->SetParameter(9, startingVWFreq);
  truthFunc->SetParameter(10, startingVWAmp);
  truthFunc->SetParameter(11, startingVWPhase);

  TF1* temp_expFunc = new TF1("temp_expFunc", "[0]*exp(-x/[1])", 0, histMaxTime);

  temp_expFunc->SetParameter(0, 1); 
  temp_expFunc->SetParameter(1, defaultLifetime);

  double Ntoy = nPts / ( temp_expFunc->Integral(0, histMaxTime) / binWidth );

  // double Ntoy = nPts / ( truthFunc->Integral(functionStart, histMaxTime) / binWidth );
  truthFunc->SetParameter(0, Ntoy);


  for (int i = 0; i < truthFunc->GetNpar(); ++i) cout << "Truth func parameter " << i << ": " << truthFunc->GetParameter(i) << endl;

  truthFunc->Write();


  // cout << "Ntoy: " << Ntoy << endl;
  // outputFile->Write();
  // delete outputFile;
  // return 1;

/////////////////////////////////////////////////////////////////////////////////////

  int totalIters = 41;

  TVectorD numPPMIters(1);
  numPPMIters[0] = totalIters;
  numPPMIters.Write("Iters");

/////////////////////////////////////////////////////////////////////////////////////

  TF1* truthFuncCopies[totalIters];
  TRandom3* truthFuncGRandoms[totalIters];
  TRandom3* randGen_UV_iters[totalIters];

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

// Fast rotation stuff from James
/*
  // Momentum
  double pmumean = 3.094349631;
  double pmusd = 0.001 * pmumean; //0.1%
  double pmurang = pmumean + pmusd;
  double pmulo = pmumean-3*pmusd;
  double pmuhi = pmumean+3*pmusd;

  // Time & Freq
  double rMagic = 711.214761;
  double trevav = 2.*TMath::Pi()*rMagic/(0.99942*29.979246); //1.0e9*2.*TMath::Pi()*711.215/(0.99942*2.9979246e10);
  double freqav = 1.0/trevav;
  double freqrang = 0.020794769/pmurang;
  double freqsd = freqav - freqrang;
  double freqlo = freqav-freqsd*3.;
  double freqhi = freqav+freqsd*3.;

  // Muon pars
  double taumu = 2196.9803; //MuLan lifetime in ns
  double gamma = 29.4;

  // T0 Pars
  double t0mean = 0.0;
  double t0sd = 15;

  // Input information for realistic pulses
  Double_t inputT0Content[26] = { 0, 1033.210, 2878.228, 10258.30, 12250.92, 12841.32, 13284.13, 15350.55, 20590.40, 24575.64, 32693.72, 37859.77, 41697.41, 38892.98, 35129.15, 26937.26, 21107.01, 15276.75, 10996.30, 9667.896, 8634.686, 7011.070, 4428.044, 1623.616, 811.80, 0};
  Double_t inputRadBins[38] = { 708.0275, 708.2379, 708.4483, 708.6589, 708.8696, 709.0804, 709.2913, 709.5024, 709.7136, 709.9249, 710.1363, 710.3479, 710.5596, 710.7714, 710.9834, 711.1954, 711.4076, 711.6200, 711.8324, 712.0450, 712.2577, 712.4706, 712.6835, 712.8966, 713.1098, 713.3232, 713.5367, 713.7503, 713.9640, 714.1778, 714.3918, 714.6059, 714.8202, 715.0346, 715.2491, 715.4637, 715.6784, 715.8933};
  Double_t inputRadContent[38] = { 0.003010, 0.001786, -0.001229, -0.001294, 0.002887, 0.003798, 0.006922, 0.007412, 0.025320, 0.080977, 0.208575, 0.401118, 0.627794, 0.819469, 0.933363, 0.996938, 1.000000, 0.989440, 0.964585, 0.931940, 0.881744, 0.817763, 0.762939, 0.687833, 0.583787, 0.434679, 0.240783, 0.096240, 0.031282, 0.009576, 0.004438, 0.005673, 0.001774, -0.000256, -0.000667, -0.000568, 0.000937, 0.005438};

  auto frDir = topDir->mkdir("FastRotationHists");
  frDir->cd();

  // Input Spectra
  TH1D* hInputT0 = new TH1D("hInputT0","Input T_{0} Shape; Time [ns]",26,-65,65);
  for(int i = 1; i <= 26; i++) hInputT0->SetBinContent(i,inputT0Content[i-1]);
  TH1D* hInputRad = new TH1D("hInputRad","Input Radius; x_{e} [cm]", 37, inputRadBins);  // Values taken from Antoine's 60h FT analysis - could easily be off by half a bin
  for(int i = 1; i <= 37; i++){
    if(inputRadContent[i-1] > 0.008) hInputRad->SetBinContent(i,inputRadContent[i-1]);
  }

  // Initialise histograms
  TH1D* hFreqSpec = new TH1D("hFreqSpec","Frequency Spectrum (GHz)",101,freqlo,freqhi);
  TH1D* hMomSpec = new TH1D("hMomSpec","Momentum Spectrum (GeV/c)",101,pmulo,pmuhi);
  TH1D* hT0Pulse = new TH1D("hT0Pulse","T0Pulse (ns) ",201,-100.5,100.5);
  const int ndet = 1;
  TH1D* hSignalRaw[ndet];
  for(int idet=0;idet<ndet;idet++){
    hSignalRaw[idet] = new TH1D(Form("hSignalRaw%02d",idet),"Raw Signal;Time [#mus]",150000,-0.0005,149.9995); //only signal 00 is used
  }
  TH1D* hTimeDecay = new TH1D("hTimeDecay","Muon Decay Time ",150000,0.0,150.0);
  TH1D* hDetector = new TH1D("hDetector","Detector Number ",24,-0.5,23.5);
  TH1D* hnTurns = new TH1D("hnTurns","Number of turns",10000,0.0,6000.0);
  TH1D* hRadius = new TH1D("hRadius","Radius of muon orbit",650,711.2-6.5,711.2+6.5); 
*/

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  TDirectory* histIterDirs[totalIters];

  TH1F* toyFiveParamHists_raw[totalIters];
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

  toyFiveParamHists_raw[iter] = new TH1F("Toy_5_Param_Hist_Raw","Toy_5_Param_Hist_Raw",nBins*100,0,histMaxTime);
  toyFiveParamHists[iter] = new TH1F("Toy_5_Param_Hist","Toy_5_Param_Hist",nBins,0,histMaxTime);
  toyUHists[iter] = new TH1F("Toy_U_Hist","Toy_U_Hist",nBins,0,histMaxTime);
  toyVHists[iter] = new TH1F("Toy_V_Hist","Toy_V_Hist",nBins,0,histMaxTime);
  toyNumHists[iter] = new TH1F("Toy_Num_Hist","Toy_Num_Hist",nBins,0,histMaxTime);
  toyDenomHists[iter] = new TH1F("Toy_Denom_Hist","Toy_Denom_Hist",nBins,0,histMaxTime);

  auto savDirs = histIterDirs[iter]->mkdir("SavedParameters");
  savDirs->cd();

  TVectorD parameterStore(2); // vector of doubles to store parameters used when creating histograms, for systematic studies
  // energy threshold, g-2 period used, possibly function parameters, etc.
  // should this just be a tntuple? or maybe a tlist with a variety of root objects? for tntuple see the ratioMacro code

  parameterStore[0] = gm2PeriodGuesses[iter];
  parameterStore[1] = 0; // supposed to be the energy threshold but that doesn't apply for this toy mc version
  parameterStore.Write("parameterStore");

/////////////////////////////////////////////////////////////////////////////////////

  truthFuncCopies[iter] = new TF1();
  truthFunc->Copy(*truthFuncCopies[iter]);
  // truthFuncCopies[iter] = (TF1*) truthFunc->Clone(Form("truthFunc_copy_%i", iter)); // cloning doesn't work for some reason, evaluating the function later even after changing parameters returns the same value - probably something to do with transient integral data arrays or something

  truthFuncCopies[iter]->SetParameter(9, (9.98 + 0.001*iter)*blindingWa*1e3);
  truthFuncCopies[iter]->SetNpx(100000);
  truthFuncCopies[iter]->Write();

  truthFuncGRandoms[iter] = new TRandom3(func_iterRandSeed);
  randGen_UV_iters[iter] = new TRandom3(UV_iterRandSeed);

}

// cout << "pointer check: " << truthFuncCopies[0] << " " << truthFuncCopies[1] << endl;
// cout << "parameter check: " << truthFuncCopies[0]->GetParameter(10) << " " << truthFuncCopies[1]->GetParameter(10) << endl;
// cout << "function value check: " << truthFuncCopies[0]->Eval(50000) << " vs: " << truthFuncCopies[1]->Eval(50000) << " diff: " << truthFuncCopies[0]->Eval(50000)-truthFuncCopies[1]->Eval(50000) << endl;



    double randEntries = randGen_positrons->PoissonD(nPts);
    double tenIncrement = 0.1 * randEntries;
    double countdown = tenIncrement;
    int tenPercent = 0;

    double savedTime = 0;

        for (double entry = 0; entry < randEntries; ++entry)
        {

          if(--countdown <= 0) // progress output
          {
            tenPercent++;
            cout << tenPercent << "0%" << " completed" << endl;
            countdown = tenIncrement;
          }

          // double time_original = truthFunc->GetRandom();

/////////////////////////////////////////////////////////////////////////////////////
/*
    //Select a T0 for the muon from histogram of t0 pulse
    double t0 = hInputT0->GetRandom();
    
    //Now select "momentum" (actually a revolution time) from input radial distribution
    double rmu = hInputRad->GetRandom();
    double pmu = pmumean*rmu/rMagic;

    // Convert to time 
    double trev = trevav*pmu/pmumean;

    //Now choose decay time from wiggle plot
    double tdec = time_original;

    //With t0 and trev, grab the momentum (not relatistic) and frequency
    double frev = 1./trev;
    hMomSpec->Fill(pmu);
    hFreqSpec->Fill(frev);

    //TDec is the time of decay for this muon
    hTimeDecay->Fill(tdec / 1e3); // us

    //and determine the radius
    double murad = rMagic*trev/trevav;
    hRadius->Fill((float)murad);
      
    //Fill t0 pulse histogram
    double t0fill = t0;
    hT0Pulse->Fill(-t0fill);  //For use with fasrof

    //Calculate number of turns, detector
    double nturns = (tdec+t0)/trev; //number of turns; this assumes that positive t0 is already past detector 0
    hnTurns->Fill(nturns);
    double ttt = nturns-(int)nturns; // this is distance past detector 0 in turns. But the corresponding time is tdec
    double tfill = tdec-(ttt*trev); //this is the detector 0 signal
*/
/////////////////////////////////////////////////////////////////////////////////////
         
              for (int iter = 0; iter < totalIters; ++iter)
              {

                gRandom = truthFuncGRandoms[iter];
                
                double time_original = truthFuncCopies[iter]->GetRandom();


/////////////////////////////////////////////////////////////////////////////////////

                // double time = time_original + (binWidth * randGen_FR->Uniform() - binWidth/2.);
                double time = time_original;
                // double time = tfill;

                toyFiveParamHists[iter]->Fill(time);
                // toyFiveParamHists_raw[iter]->Fill(time);

                // double randNum = randGen_UV->Uniform();
                double randNum = randGen_UV_iters[iter]->Uniform();

                double halfPeriodGuess = gm2PeriodGuesses[iter]/2;

                double totalChance = exp(halfPeriodGuess/defaultLifetime) + exp(-halfPeriodGuess/defaultLifetime) + 2.;
                double percentChanceUPlus = exp(halfPeriodGuess/defaultLifetime) / totalChance;
                double percentChanceUMinus = exp(-halfPeriodGuess/defaultLifetime) / totalChance;

                if     (randNum < percentChanceUPlus) toyUHists[iter]->Fill(time - halfPeriodGuess); // careful with the signs here - the U+ hist is filled with pulses shifted by t - T/2, and weighted by e^+T/2tau, and vice versa for U-
                else if(randNum < (percentChanceUPlus + percentChanceUMinus)) toyUHists[iter]->Fill(time + halfPeriodGuess);
                else if(randNum < 1) toyVHists[iter]->Fill(time);

              }
        }

        for (int iter = 0; iter < totalIters; ++iter)
        {
          toyNumHists[iter]->Add(toyUHists[iter], toyVHists[iter], -1, 1);
          toyDenomHists[iter]->Add(toyUHists[iter], toyVHists[iter]);
        }

/////////////////////////////////////////////////////////////////////////////////////

    t = clock() - t;
    printf ("It took me %f seconds.\n", ((float)t)/CLOCKS_PER_SEC);


/////////////////////////////////////////////////////////////////////////////////////

  outputFile->Write();
  delete outputFile;

  return 1;

}
