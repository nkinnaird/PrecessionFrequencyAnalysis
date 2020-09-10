#include <algorithm>
#include <iostream>
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TMath.h"

void FastRotationToyMC(){

  //define a few constants

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

  // Initialise 5-parameter wiggle function
  TF1* fWiggle = new TF1("fWiggle","[0]*exp(-x/[1])*(1+[2]*cos([3]*x + [4]))",0,150000);
  fWiggle->SetParameters(1, 64372, 0.4, 0.00143944, 0.0);
  fWiggle->SetNpx(150000);
  
  // Main event loop
  TRandom3* re = new TRandom3(6280225);
  gRandom->SetSeed(12345); // Set seed for fWiggle->GetRandom()
  long long nEvents = 100000000;
  for(long long iev = 0; iev < nEvents; iev++){
  
    if(iev%1000000 ==0) cout << "Muon number " << iev << endl;

    //Select a T0 for the muon from histogram of t0 pulse
    double t0 = hInputT0->GetRandom();
    
    //Now select "momentum" (actually a revolution time) from input radial distribution
    double rmu = hInputRad->GetRandom();
    double pmu = pmumean*rmu/rMagic;

    // Convert to time 
    double trev = trevav*pmu/pmumean;

    //Now choose decay time from wiggle plot
    double tdec = fWiggle->GetRandom();

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
    int cdet = (int)(ndet*ttt); 

    hDetector->Fill((float)cdet);
    if(cdet > -1 && cdet < ndet) {
      hSignalRaw[cdet]->Fill(tfill/1e3);   //filling only detector 0  (us)
    }

  }

  // Write out histograms to file
  string outfile = "toyMCOutput.root";
  TFile *fout = new TFile(outfile.c_str(),"RECREATE");
  hDetector->Write();
  hMomSpec->Write();
  hnTurns->Write();
  hFreqSpec->Write();
  hT0Pulse->Write();
  hRadius->Write();
  hTimeDecay->Write();
  for(int idet=0;idet<ndet;idet++) {
    hSignalRaw[idet]->Write();
  }
  fout->Write();
  fout->Close();
}
