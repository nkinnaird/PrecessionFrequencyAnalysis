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
double nPts = 1.5e6;

// Was total number of sample iterations - now number of random seeds
int totalIters = 100;

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

// Muon lifetime [ns]
//static const double defaultLifetime = 64440; // ns
static const double defaultLifetime = 64400; // ns // number used in the production of the large 2D histogram                                                                                        

// Precession frequency values (blinding)
// 0.2291 MHz converted to units of 1/ns, 0.2291 is the value used in the blinding software                                             
static const double blindingFa = 0.2291 * 1e6 * 1/(1e9);
static const double blindingWa = 2*pi*blindingFa;


// Main program
// int makeToyHistsAlexModified()
int makeToyHistsAlexModified(string filePath)
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
  TRandom3* randGen_FR = new TRandom3(0);
  TRandom3* randGen_UV = new TRandom3(0);

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  // Create output file that will hold MC histograms
  // TFile* outputFile = new TFile("15e9toyEvents_500samplesGrid.root","RECREATE");
  // TFile* outputFile = new TFile("15e9toyEvents_50samplesGrid.root","RECREATE");
   TFile* outputFile = new TFile("toyHistsMethods.root","RECREATE");

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
  // Ratio histogram                                    
  TH1F* toyUHists[totalIters];
  TH1F* toyVHists[totalIters];
  // Q-method histogram                           
  TH1F* toyFiveParamQHists[totalIters];
  // Posistron energies and energy binned histogram
  TH1F* energies[totalIters];
  TH1F* toyEbins[totalIters][eBins];

  /////////////////////////////////////////////////////////////////////////////////////  

  double halfPeriodGuess = (1/blindingFa)/2;

  double totalChance = exp(halfPeriodGuess/defaultLifetime) + exp(-halfPeriodGuess/defaultLifetime) + 2.;
  double percentChanceUPlus = exp(halfPeriodGuess/defaultLifetime) / totalChance;
  double percentChanceUMinus = exp(-halfPeriodGuess/defaultLifetime) / totalChance;


  /////////////////////////////////////////////////////////////////////////////////////

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
    // Ratio histogram                                                                                 
    toyUHists[iter] = new TH1F("Toy_U_Hist","Toy_U_Hist; Time (ns); Events",nBins,0,histMaxTime);
    toyVHists[iter] = new TH1F("Toy_V_Hist","Toy_V_Hist; Time (ns); Events",nBins,0,histMaxTime);
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

  /////////////////////////////////////////////////////////////////////////////////////      
  /////////////////////////////////////////////////////////////////////////////////////
    
  // Define random number of entries from Poisson distribution
  double randEntries = randGen_positrons->PoissonD(nPts);


  clock_t startTime = clock();
  int percentIncrement = 10;
  int eventCountdown = 0.01*percentIncrement*randEntries;
  double percentDone = 0;


  // Initialise positron counter
  int64_t positronNum = 0;

  while(positronNum <= randEntries){

    /////////////////////////////////////////////////////////////////////////////////////

    if(--eventCountdown == 0)
    {
      percentDone+=percentIncrement;
      double elapsedTime = double(clock() - startTime)/CLOCKS_PER_SEC;
      double timeRemainingSeconds = (100 - percentDone)/percentDone * elapsedTime;
      cout << percentDone << "% processed. Time passed: " << int(elapsedTime/60/60) << " hours, " << int(elapsedTime/60)%60 << " minutes, and " << int(elapsedTime)%60 << " seconds." << " Estimated time remaining = " << int(timeRemainingSeconds/60/60) << " hours, " << int(timeRemainingSeconds/60)%60 << " minutes, and " << int(timeRemainingSeconds)%60 << " seconds." << endl;
      eventCountdown = 0.01*percentIncrement*randEntries;
    }

/////////////////////////////////////////////////////////////////////////////////////

    // Define time and energy variables
    double hit_time;
    double hit_energy_raw;
    
    // Get time and energy from 2D histogram
    initialHist->GetRandom2(hit_time, hit_energy_raw);


    // Loop over number of fills
    for (int iter = 0; iter < totalIters; ++iter){

      double energy = hit_energy_raw; // energy
      double time = hit_time + (binWidth * (randGen_FR->Uniform() - 0.5)); // smeared time
      
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
      	  double A = muonA->Eval(energy);
      	  toyFiveParamAwgtHists[iter]->Fill(time,A);
      	} 
      	
      	  // Fill T-method and Ratio histogram above T-method threshold
      	if(energy > Tth){
      	   toyFiveParamTwgtHists[iter]->Fill(time);
      	    
            double randNum = randGen_UV->Uniform();

      	    // Careful with the signs here - the U+ hist is filled with pulses shifted by t - T/2, 
      	    // and weighted by e^+T/2tau, and vice versa for U-
      	    if     (randNum < percentChanceUPlus) toyUHists[iter]->Fill(time - halfPeriodGuess);
      	    else if(randNum < (percentChanceUPlus + percentChanceUMinus)) toyUHists[iter]->Fill(time + halfPeriodGuess);
      	    else if(randNum < 1) toyVHists[iter]->Fill(time);
      	} // End if E > T-method threshold

      	// Fill Q-method histograms above Q-method threshold
      	if(energy > Qth){
      	  toyFiveParamQHists[iter]->Fill(time,energy);
      	}

      } // End if E < Emax
      
    } // End sub loop over iterations - random seed

    positronNum += 1;    
  } // End loop over positron number 


      
  outputFile->Write();
  delete outputFile;

  /////////////////////////////////////////////////////////////////////////////////////

  double runTime = double(clock() - startTime)/CLOCKS_PER_SEC;
  cout << "Total run time: " << int(runTime/60/60) << " hours, " << int(runTime/60)%60 << " minutes, and " << int(runTime)%60 << " seconds." << endl;

  /////////////////////////////////////////////////////////////////////////////////////    

  return 1;
  
}
 
