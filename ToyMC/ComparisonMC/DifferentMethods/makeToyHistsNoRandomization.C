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

#include <time.h>

using namespace std;

// Global variable definition and initialisation

// Total number of generated events
double nPts = 1e7;

// Number of sample iterations
int totalIters = 750;

// Total number of 100 MeV enegry bins (3.1 GeV max)
int eBins = 31;

// pi
static const double pi = 3.14159265358979323846;

// Approximum maximum histogram time
static const double approxMaxTime = 700000; // ns
// Histgram bin width (cyclotron period)
static const double binWidth = 149.15; // ns
static const double QbinWidth = 75; // ns
// Total number of bins
static const int nBins = int(approxMaxTime/binWidth);
static const int nQBins = int(approxMaxTime/QbinWidth);
// Maximum histogram time
static const double histMaxTime = nBins*binWidth;


/////////////////////////////////////////////////////////////////////////////////////

  // asymmetry function data from David

  vector<double> binCenters = {1025,1075,1125,1175,1225,1275,1325,1375,1425,1475,1525,1575,1625,1675,1725,1775,1825,1875,1925,1975,2025,2075,2125,2175,2225,2275,2325,2375,2425,2475,2525,2575,2625,2675,2725,2775,2825,2875,2925,2975};
  vector<double> A_60h = {0.00345991,0.013643,0.0245829,0.0362506,0.0490243,0.0625439,0.0764485,0.0910481,0.106377,0.122331,0.138858,0.156191,0.173996,0.192135,0.211069,0.230611,0.25058,0.271064,0.292124,0.313733,0.33603,0.358879,0.382346,0.406546,0.431054,0.455926,0.481466,0.507351,0.53356,0.559992,0.586049,0.611355,0.635609,0.658539,0.67937,0.69775,0.713558,0.726565,0.736108,0.742875};
  vector<double> A_HK = {0.00341471,0.0135044,0.0246224,0.0365415,0.0491475,0.0624716,0.0764926,0.0913398,0.106797,0.1226,0.139196,0.156396,0.174055,0.192553,0.211505,0.230849,0.250817,0.271251,0.292469,0.314358,0.336572,0.3593,0.382427,0.406379,0.43115,0.456207,0.481599,0.50755,0.534137,0.560397,0.586063,0.611518,0.63593,0.658607,0.679437,0.697708,0.71328,0.725815,0.735925,0.743152};
  vector<double> A_9d = {0.00328961,0.0133864,0.0245608,0.0365693,0.0489681,0.0621252,0.0761147,0.0908883,0.106353,0.122103,0.138547,0.155785,0.173598,0.191792,0.210475,0.23001,0.249958,0.270427,0.291451,0.313075,0.335401,0.358073,0.38131,0.405118,0.429543,0.454546,0.480038,0.506037,0.532116,0.558244,0.584374,0.609723,0.633986,0.656851,0.677732,0.696337,0.712076,0.725433,0.736151,0.742911};
  vector<double> A_EG = {0.00279735,0.0129582,0.024157,0.0361966,0.0488911,0.0622554,0.0764188,0.0912046,0.106639,0.12274,0.139354,0.156562,0.174329,0.192593,0.211408,0.230788,0.250743,0.271245,0.292199,0.313791,0.335934,0.358526,0.381776,0.405558,0.429813,0.454828,0.480393,0.506124,0.532004,0.558,0.584026,0.609238,0.63324,0.65588,0.676453,0.694782,0.711075,0.725051,0.736272,0.744989};

/////////////////////////////////////////////////////////////////////////////////////


// Main program
int makeToyHistsNoRandomization(string filePath)
{

  // Path name containing 2D time-energy MC histogram
  //  string filePath = "/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_03/srcs/gm2analyses/macros/RatioMacro/ToyMC/TimeEnergyPairsHist.root";
  // string filePath = "TimeEnergyPairsHist.root";

  TFile *inputFile = TFile::Open(filePath.c_str());
  if (inputFile == 0) {
    printf("Error: cannot open file\n");
    return 0;
  }
  // Get 2D histogram time-energy MC histogram
  TH2F* initialHist = (TH2F*) inputFile->Get("full2Dfunction_hist");

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////
  
  // Set batch mode to true for this macro so that nothing draws to the screen
  gROOT->SetBatch(kTRUE); 

  // Define random seeds for MC
  gRandom->SetSeed(0); // used in GetRandom on TF1
  TRandom3* randGen_positrons = new TRandom3(0); // for number of entries

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  // Create output file that will hold MC histograms
   TFile* outputFile = new TFile("toyHistsMethodsNoRand.root","RECREATE");

  // Make top directory for output file
  auto topDir = gFile->mkdir("topDir");

  /////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  // Define vector which store number of iterations
  TVectorD numIters(1);
  numIters[0] = totalIters;
  numIters.Write("Iters");
  TDirectory* histIterDirs[totalIters];

  // Initialise histograms below
  
  // Standard T-method histogram                  
  TH1F* toyFiveParamTwgtHists[totalIters];
  // A-weighted histogram                           
  TH1F* toyFiveParamAwgtHists[totalIters];
  TH1F* toyFiveParamAwgtHists_60h[totalIters];
  TH1F* toyFiveParamAwgtHists_HK[totalIters];
  TH1F* toyFiveParamAwgtHists_9d[totalIters];
  TH1F* toyFiveParamAwgtHists_EG[totalIters];
  // Q-method histogram                           
  TH1F* toyFiveParamQHists[totalIters];
  // Posistron energies and energy binned histogram
  TH1F* energies[totalIters];
  TH1F* toyEbins[totalIters][eBins];


  ///////////////////////////////////////////////////////////////////////////////////// 

  // Loop over number of MC sample iterations - now random seeds
  for (int iter = 0; iter < totalIters; ++iter){ // loop over totalIters
    
    // Make directory for each iteration and change directory to it
    histIterDirs[iter] = topDir->mkdir(Form("Iter%d",iter));
    histIterDirs[iter]->cd();

    // Standard T-method histogram                                                                     
    toyFiveParamTwgtHists[iter] = new TH1F("Toy_5_Param_Hist_T","Toy_5_Param_Hist_T; Time (ns); Events",nBins,0,histMaxTime);
    // A-weighted histogram                                                                            
    toyFiveParamAwgtHists[iter] = new TH1F("Toy_5_Param_Hist_A","Toy_5_Param_Hist_A; Time (ns); Events",nBins,0,histMaxTime);
    toyFiveParamAwgtHists_60h[iter] = new TH1F("Toy_5_Param_Hist_A_60h","Toy_5_Param_Hist_A_60h; Time (ns); Events",nBins,0,histMaxTime);
    toyFiveParamAwgtHists_HK[iter] = new TH1F("Toy_5_Param_Hist_A_HK","Toy_5_Param_Hist_A_HK; Time (ns); Events",nBins,0,histMaxTime);
    toyFiveParamAwgtHists_9d[iter] = new TH1F("Toy_5_Param_Hist_A_9d","Toy_5_Param_Hist_A_9d; Time (ns); Events",nBins,0,histMaxTime);
    toyFiveParamAwgtHists_EG[iter] = new TH1F("Toy_5_Param_Hist_A_EG","Toy_5_Param_Hist_A_EG; Time (ns); Events",nBins,0,histMaxTime);
    // Q-method histogram                                                                            
    toyFiveParamQHists[iter] = new TH1F("Toy_5_Param_Hist_Q","Toy_5_Param_Hist_Q; Time (ns); Energy sum (Mev)",nQBins,0,histMaxTime);
    // Histogram of all event energies
    energies[iter] = new TH1F("energies","energies; Energy (MeV); Events",700,0,7000);
    // Energy bin histograms
    for (int bin = 0; bin < eBins; ++bin){
      toyEbins[iter][bin] = new TH1F(Form("%d - %d MeV",bin*100,(bin+1)*100),Form("%d - %d MeV",bin*100,(bin+1)*100),nBins,0,histMaxTime);
    }
    
  }

  /////////////////////////////////////////////////////////////////////////////////////       
  /////////////////////////////////////////////////////////////////////////////////////       
  
  // Define energy limits
  double Tth = 1680; // T-method threshold [MeV]
  double Ath = 1100; // A-weighted threshold [MeV]
  double Qth = 250; // A-weighted threshold [MeV]
  double Emax = 3100; // Maximum energy [MeV]

  // Muon asymmetry function (lab frame)
  TF1 *muonA = new TF1("muonA", "(-8.0*(x/[0])*(x/[0]) + (x/[0]) + 1.1)/(4.0*(x/[0])*(x/[0]) - 5.0*(x/[0]) - 5.0)", 0, Emax);
  muonA->SetParameter(0,Emax);


  TH1F* muonAHist_60h = new TH1F("muonAHist_60h","muonAHist_60h",40,1000,3000);
  for (int bin = 1; bin <= muonAHist_60h->GetNbinsX(); ++bin) muonAHist_60h->SetBinContent(bin, A_60h.at(bin-1));
  
  TH1F* muonAHist_HK = new TH1F("muonAHist_HK","muonAHist_HK",40,1000,3000);
  for (int bin = 1; bin <= muonAHist_HK->GetNbinsX(); ++bin) muonAHist_HK->SetBinContent(bin, A_HK.at(bin-1));

  TH1F* muonAHist_9d = new TH1F("muonAHist_9d","muonAHist_9d",40,1000,3000);
  for (int bin = 1; bin <= muonAHist_9d->GetNbinsX(); ++bin) muonAHist_9d->SetBinContent(bin, A_9d.at(bin-1));

  TH1F* muonAHist_EG = new TH1F("muonAHist_EG","muonAHist_EG",40,1000,3000);
  for (int bin = 1; bin <= muonAHist_EG->GetNbinsX(); ++bin) muonAHist_EG->SetBinContent(bin, A_EG.at(bin-1));


  /////////////////////////////////////////////////////////////////////////////////////      
  /////////////////////////////////////////////////////////////////////////////////////
    


  clock_t startTime = clock();
  // int percentIncrement = 10;
  // int eventCountdown = 0.01*percentIncrement*randEntries;
  // double percentDone = 0;


// Loop over number of samples
for (int iter = 0; iter < totalIters; ++iter){

  if(iter > 0){
    double elapsedTime = double(clock() - startTime)/CLOCKS_PER_SEC;
    double timeRemainingSeconds = (double(totalIters) - iter)/iter * elapsedTime;
    cout << "Iters done: " << iter << "/" << totalIters << ". Time passed: " << int(elapsedTime/60/60) << " hours, " << int(elapsedTime/60)%60 << " minutes, and " << int(elapsedTime)%60 << " seconds." << " Estimated time remaining = " << int(timeRemainingSeconds/60/60) << " hours, " << int(timeRemainingSeconds/60)%60 << " minutes, and " << int(timeRemainingSeconds)%60 << " seconds." << endl;
  }

  // Initialise positron counter
  int64_t positronNum = 0;

  // Define random number of entries from Poisson distribution
  double randEntries = randGen_positrons->PoissonD(nPts);

  while(positronNum <= randEntries){

/////////////////////////////////////////////////////////////////////////////////////

    // if(--eventCountdown == 0)
    // {
    //   percentDone+=percentIncrement;
    //   double elapsedTime = double(clock() - startTime)/CLOCKS_PER_SEC;
    //   double timeRemainingSeconds = (100 - percentDone)/percentDone * elapsedTime;
    //   cout << percentDone << "% processed. Time passed: " << int(elapsedTime/60/60) << " hours, " << int(elapsedTime/60)%60 << " minutes, and " << int(elapsedTime)%60 << " seconds." << " Estimated time remaining = " << int(timeRemainingSeconds/60/60) << " hours, " << int(timeRemainingSeconds/60)%60 << " minutes, and " << int(timeRemainingSeconds)%60 << " seconds." << endl;
    //   eventCountdown = 0.01*percentIncrement*randEntries;
    // }

/////////////////////////////////////////////////////////////////////////////////////

    // Define time and energy variables
    double hit_time;
    double hit_energy_raw;
    
    // Get time and energy from 2D histogram
    initialHist->GetRandom2(hit_time, hit_energy_raw);


      double energy = hit_energy_raw; // energy
      double time = hit_time;
      
      /////////////////////////////////////////////////////////////////////////////////////      
      
      // Fill histogram of all event energies
      energies[iter]->Fill(energy);
      
      // Ensure energies are below the energy max (no pileup)
      if(energy < Emax){
	
      	// Fill E-bin histograms
      	int binNo = int(floor(energy/100));
      	toyEbins[iter][binNo]->Fill(time);
      	
      	// Fill A-weighted histogram above A-weighted threshold - with muon asymmetry value
      	if(energy > Ath){
      	  double Afunc = muonA->Eval(energy);
          toyFiveParamAwgtHists[iter]->Fill(time,Afunc);
          
          double Ahist_60h = muonAHist_60h->GetBinContent(muonAHist_60h->FindBin(energy));
          double Ahist_HK = muonAHist_HK->GetBinContent(muonAHist_HK->FindBin(energy));
          double Ahist_9d = muonAHist_9d->GetBinContent(muonAHist_9d->FindBin(energy));
          double Ahist_EG = muonAHist_EG->GetBinContent(muonAHist_EG->FindBin(energy));
          
          if(energy > 3000) Ahist_60h = muonAHist_60h->GetBinContent(muonAHist_60h->GetNbinsX());
          if(energy > 3000) Ahist_HK = muonAHist_60h->GetBinContent(muonAHist_HK->GetNbinsX());
          if(energy > 3000) Ahist_9d = muonAHist_60h->GetBinContent(muonAHist_9d->GetNbinsX());
          if(energy > 3000) Ahist_EG = muonAHist_60h->GetBinContent(muonAHist_EG->GetNbinsX());

          toyFiveParamAwgtHists_60h[iter]->Fill(time,Ahist_60h);
          toyFiveParamAwgtHists_HK[iter]->Fill(time,Ahist_HK);
          toyFiveParamAwgtHists_9d[iter]->Fill(time,Ahist_9d);
          toyFiveParamAwgtHists_EG[iter]->Fill(time,Ahist_EG);

          // if(energy > 3000) cout << "Ahist vs Afunc: " << Ahist << " " << Afunc << " " << Ahist-Afunc << endl;
      	} 
      	
      	 // Fill T-method and Ratio histogram above T-method threshold
      	if(energy > Tth){
      	   toyFiveParamTwgtHists[iter]->Fill(time);
      	} 

      	// Fill Q-method histograms above Q-method threshold
      	if(energy > Qth){
      	  toyFiveParamQHists[iter]->Fill(time,energy);
      	}

      } // End if E < Emax

      positronNum += 1;          
    } // End loop over positron number 

} // End sub loop over iterations - random seed


      
  outputFile->Write();
  delete outputFile;

  /////////////////////////////////////////////////////////////////////////////////////

  double runTime = double(clock() - startTime)/CLOCKS_PER_SEC;
  cout << "Total run time: " << int(runTime/60/60) << " hours, " << int(runTime/60)%60 << " minutes, and " << int(runTime)%60 << " seconds." << endl;

  /////////////////////////////////////////////////////////////////////////////////////    

  return 1;
  
}
 
