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
double nPts = 1.5e7; 
// 1.5e9 for about the same precision as the 60 hr dataset 
// Can increase the stats or run more jobs on the grid and hadd them

// Total number of sample iterations
int totalIters = 500;
//int totalIters = 50;

// Total number of 100 MeV enegry bins (3.1 GeV max)
int eBins = 31;

// pi
static const double pi = 3.14159265358979323846;

// Number of calorimeters
static const int nCalos = 24;

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

// Define number of fit parameters
static const int numTmethodParams = 5; // T/A-method
static const int numRatioParams = 3; // Ratio-method

// Fit start time
static const double fitStart = 30000; // ns
// Fit end time
static const double fitEnd = 350000; // ns
static const double ratioFitEnd = 350000; // ns
static const double QFitEnd = 218000; // ns

// Muon lifetime [ns]
//static const double defaultLifetime = 64440; // ns                                                                                         
static const double defaultLifetime = 64400; // ns                                                                                         

// Precession frequency values (blinding)
// 0.2291 MHz converted to units of 1/ns, 0.2291 is the value used in the blinding software                                             
static const double blindingFa = 0.2291 * 1e6 * 1/(1e9);
static const double blindingWa = 2*pi*blindingFa;

// Precession frequency values (no blinding)
// 0.2290735 is the average value in column 2 of Table 15 of the E821 Final Report, for the fa values for the different run periods
static const double defaultFa = 0.2290735 * 1e6 * 1/(1e9); 
static const double defaultWa = 2*pi*defaultFa;

// g-2 period
static const double g2Period = 1/defaultFa;

// Five parameter function
class wiggleFitClass{

public:
  wiggleFitClass(int numTotalParameters)
  {
    numFunctionParameters = numTotalParameters;
  }

  double fiveParFunc(double* x, double* par){
    double time = x[0];
    double freq = blindingWa * (1 + par[3] * 1e-6);  // want to use frequency corresponding to blinding software since ToyMC fits are unblinded
    double lifetime = par[1] * 1000.; // convert from us to ns for fit function

    double value = par[0] * exp(-time/lifetime) * (1 + par[2] * cos(freq*time + par[4]));

    return value;
  };

  double threeParFunc(double* x, double* par){
    double time = x[0];
    double freq = blindingWa * (1 + par[1] * 1e-6);  // want to use frequency corresponding to blinding software since ToyMC fits are unblinded

    double value = par[0] * cos(freq*time + par[2]);

    return value;
  };


private:
  int numFunctionParameters;

}; // end wiggleFitClass

// Main program
int makeToyHists()
{
  // Store program start time
  clock_t startTime = clock();

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  // Path name containing 2D time-energy MC histogram
  //  string filePath = "/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_03/srcs/gm2analyses/macros/RatioMacro/ToyMC/TimeEnergyPairsHist.root";
  string filePath = "TimeEnergyPairsHist.root";
  // Open file
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
  //  TFile* outputFile = new TFile("toyHists_test.root","RECREATE");
  //  TFile* outputFile = new TFile("/pnfs/GM2/scratch/users/alexkesh/toyEvents/22-11-19/15e8toyEvents_100samplesGrid.root","RECREATE");
  TFile* outputFile = new TFile("15e9toyEvents_500samplesGrid.root","RECREATE");
  //  TFile* outputFile = new TFile("15e9toyEvents_50samplesGrid.root","RECREATE");

  // Make top directory for output file
  auto topDir = gFile->mkdir("topDir");

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  // Build five-parameter fit function
  auto fiveParamFit = new wiggleFitClass(numTmethodParams);
  auto wiggleFit = new TF1("wiggleFit", fiveParamFit, &wiggleFitClass::fiveParFunc, 0, histMaxTime, numTmethodParams);
  // Set large number of points for fit function
  wiggleFit->SetNpx(10000);

  // Initiliase fit parameters
  wiggleFit->SetParameter(0, 1);
  wiggleFit->SetParameter(1, defaultLifetime/1000.);
  wiggleFit->SetParameter(2, 0.4);
  wiggleFit->SetParameter(3, 0);
  wiggleFit->SetParameter(4, pi);

  // Determine and set fit normalisation value
  double Ntoy = nPts / ( wiggleFit->Integral(0, histMaxTime) / binWidth );
  wiggleFit->SetParameter(0, Ntoy);

  // Write truth fit function to file
  wiggleFit->Write();

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
  TH1F* toyNumHists[totalIters];
  TH1F* toyDenomHists[totalIters];
  TH1F* toyRatioHists[totalIters];
  // Q-method histogram                           
  TH1F* toyFiveParamQHists[totalIters];
  // Posistron energies and energy binned histogram
  TH1F* energies[totalIters];
  TH1F* toyEbins[totalIters][eBins];

  /////////////////////////////////////////////////////////////////////////////////////  

  // g-2 period guesses per iteration
  double gm2PeriodPPMOffsets[totalIters];
  double gm2PeriodGuesses[totalIters];

  /////////////////////////////////////////////////////////////////////////////////////

  // Randomisation per iteration
  TRandom3* randIterUV[totalIters]; 

  ///////////////////////////////////////////////////////////////////////////////////// 

  // Loop over number of MC sample iterations
  for (int iter = 0; iter < totalIters; ++iter){ // loop over totalIters
    
    gm2PeriodPPMOffsets[iter] = 0;
    // gm2PeriodGuesses[iter] = g2Period * (1 + 1e-6 * gm2PeriodPPMOffsets[iter]);
    // g2Period corresponding to blinding software since ToyMC fits are unblinded
    gm2PeriodGuesses[iter] = 1/blindingFa * (1 + 1e-6 * gm2PeriodPPMOffsets[iter]); 
    
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
    toyNumHists[iter] = new TH1F("Toy_Num_Hist","Toy_Num_Hist; Time (ns); Events",nBins,0,histMaxTime);
    toyDenomHists[iter] = new TH1F("Toy_Denom_Hist","Toy_Denom_Hist; Time (ns); Events",nBins,0,histMaxTime);
    toyRatioHists[iter] = new TH1F("Toy_Ratio_Hist","Toy_Ratio_Hist; Time (ns); Events",nBins,0,histMaxTime);
    // Q-method histogram                                                                            
    toyFiveParamQHists[iter] = new TH1F("Toy_5_Param_Hist_Q","Toy_5_Param_Hist_Q; Time (ns); Energy sum (Mev)",nQBins,0,histMaxTime);
    // Histogram of all event energies
    energies[iter] = new TH1F("energies","energies; Energy (MeV); Events",700,0,7000);
    // Energy bin histograms
    for (int bin = 0; bin < eBins; ++bin){
      toyEbins[iter][bin] = new TH1F(Form("%d - %d MeV",bin*100,(bin+1)*100),Form("%d - %d MeV",bin*100,(bin+1)*100),nBins,0,histMaxTime);
    }
    
    /////////////////////////////////////////////////////////////////////////////////////     
    
    // Make directory of saved parameters and change directory to it
    auto savDirs = histIterDirs[iter]->mkdir("SavedParameters");
    savDirs->cd();
    
    TVectorD parameterStore(2); // vector of doubles to store parameters used when creating histograms, for systematic studies
    // energy threshold, g-2 period used, possibly function parameters, etc.
    // should this just be a tntuple? or maybe a tlist with a variety of root objects? for tntuple see the ratioMacro code
    
    parameterStore[0] = gm2PeriodGuesses[iter];
    parameterStore.Write("parameterStore");

    // Choose seed for iteration randomisation
    randIterUV[iter] = new TRandom3(31568);

  }

  /////////////////////////////////////////////////////////////////////////////////////       
  /////////////////////////////////////////////////////////////////////////////////////       
  
  // Initialise progress counter variables
  double targetPerc = 0; // Printout for percentage complete 
  int64_t eventCount = 0;
  int64_t totEvents = nPts*totalIters;
  clock_t refTime = clock();
  double tenPerTime = 0;
  
  /////////////////////////////////////////////////////////////////////////////////////       

  // Define vector to store hits per fill
  std::vector<std::pair<double, double> > trueHits;

  // Define energy limits
  double Tth = 1680; // T-method threshold [MeV]
  double Ath = 1100; // A-weighted threshold [MeV]
  double Qth = 250; // A-weighted threshold [MeV]
  double Emax = 3100; // Maximum energy [MeV]

  // Muon asymmetry function (lab frame)
  TF1 *muonA = new TF1("muonA", "(-8.0*(x/[0])*(x/[0]) + (x/[0]) + 1.1)/(4.0*(x/[0])*(x/[0]) - 5.0*(x/[0]) - 5.0)", 0, Emax);
  muonA->SetParameter(0,Emax);

  /////////////////////////////////////////////////////////////////////////////////////      

  // Begin filling histograms
  // Loop over number of sample iterations
  startTime = clock();
  for (int iter = 0; iter < totalIters; ++iter){
    
    // Define random number of entries from Poisson distribution
    double randEntries = randGen_positrons->PoissonD(nPts);

    // Initialise positron counter
    int64_t positronNum = 0;

    // Loop over number of fills
    while(positronNum <= randEntries){

      // Print the percentage of entries that have been processd and the time remaining                                              
      if(100*float(eventCount) / totEvents > targetPerc){
	// Calculate percentage of events processed
	double percent = int(100*float(eventCount) / totEvents);
	// Determine how long it took to process 10% of the events
	if (int(percent) == 10) {
	  refTime = clock();
	  tenPerTime = double(refTime - startTime)/CLOCKS_PER_SEC;
	}
	// Output the percentage completed for each 10% and time reminaing
	if (int(percent) % 10 == 0) {
	  if (int(percent) >= 10){
	    int seconds = int((100-percent)/10*tenPerTime);
	    int minutes = seconds / 60;
	    int hours = minutes / 60;
	    if (int(hours) > 0 && int(minutes%60) > 0) {
	      cout << "Processed " << int(percent)  << "%; Estimated time remaining = " << int(hours) << " hours, " << int(minutes%60) << " minutes & " << int(seconds%60) << " seconds." << endl;
	    }
	    else if (int(minutes%60) > 0) {
	      cout << "Processed " << int(percent)  << "%; Estimated time remaining = " << int(minutes%60) << " minutes & " << int(seconds%60) << " seconds." << endl;
	    }
	    else {
	      cout << "Processed " << int(percent)  << "%; Estimated time remaining = " << int(seconds%60) << " seconds." << endl;
	    }
	  }
	  else{
	    cout << "Processed " << int(percent)  << "%" << endl;
	  }
	}
	// Update percentage process by 1
	targetPerc += 1;
      }
      
      // Define time and energy variables
      double hit_time;
      double hit_energy_raw;
      
      // Get time and energy from 2D histogram
      initialHist->GetRandom2(hit_time, hit_energy_raw);

      /////////////////////////////////////////////////////////////////////////////////////            
            
      // Define time and energy variables for histogram filling
      double time_original = hit_time; // Event time
      double energy = hit_energy_raw; // Event enrgy
      double randNum = randGen_UV->Uniform(); // Random number generator for FR smearing 
      // Define smeared time
      double time = time_original + (binWidth * randGen_FR->Uniform() - binWidth/2.);
      
      /////////////////////////////////////////////////////////////////////////////////////      
      
      // Fill histogram of all event energies
      energies[iter]->Fill(energy);
      
      // Ensure energies are below the energy max (no pileup)
      if(energy < Emax){
	
	// Fill E-bin histograms                                                               
	int binNo = int(floor(energy/100));
	toyEbins[iter][binNo]->Fill(time);
	
	// Fill A-weighted histogram above A-weighted threshold
	if(energy > Ath){
	  
	  // Find value of muon asymmetry for given event energy
	  double A = muonA->Eval(energy);
	  // Fill histogram and weight by muon asymmetry
	  toyFiveParamAwgtHists[iter]->Fill(time,A);
	  
	} // End if E > A-weighted threshold
	
	  // Fill T-method and Ratio histogram above T-method threshold
	if(energy > Tth){
	  
	  // Fill T-method histograms
	    toyFiveParamTwgtHists[iter]->Fill(time);
	    
	    // Fill Ratio histograms
	    double halfPeriodGuess = gm2PeriodGuesses[iter]/2;
	    double totalChance = exp(halfPeriodGuess/defaultLifetime) + exp(-halfPeriodGuess/defaultLifetime) + 2.;
	    double percentChanceUPlus = exp(halfPeriodGuess/defaultLifetime) / totalChance;
	    double percentChanceUMinus = exp(-halfPeriodGuess/defaultLifetime) / totalChance;
	    // Careful with the signs here - the U+ hist is filled with pulses shifted by t - T/2, 
	    // and weighted by e^+T/2tau, and vice versa for U-
	    if     (randNum < percentChanceUPlus) toyUHists[iter]->Fill(time - halfPeriodGuess);
	    else if(randNum < (percentChanceUPlus + percentChanceUMinus)) toyUHists[iter]->Fill(time + halfPeriodGuess);
	    else if(randNum < 1) toyVHists[iter]->Fill(time);
	    
	    positronNum += 1;
	    eventCount += 1;
	    
	} // End if E > T-method threshold

	// Fill Q-method histograms above Q-method threshold
	if(energy > Qth){

	  // Fill Q-method histogram with positron energy
	  toyFiveParamQHists[iter]->Fill(time,energy);
	  
	}

      } // End if E < Emax
      
    } // End sub loop over number of fills
    
    // Calculate ratio numerator and denominator
    toyNumHists[iter]->Add(toyUHists[iter], toyVHists[iter], -1, 1);
    toyDenomHists[iter]->Add(toyUHists[iter], toyVHists[iter]);
    // Calculate ratio histogram and errors
    for (int i = 0; i < nBins; i++){
      double U = toyUHists[iter]->GetBinContent(i);
      double V = toyVHists[iter]->GetBinContent(i);
      double num = V-U;
      double den = V+U;
      double R = 0;
      double err = 0;
      if (den == 0){
	R = 0;
	err = 0;
      }
      else{
	R = (V-U)/(V+U);
	err = sqrt((1-pow(R,2))/(V+U));
      }
      toyRatioHists[iter]->SetBinContent(i,R);
      toyRatioHists[iter]->SetBinError(i,err);
    }
    
  } // End loop over sample iterations 

  // TGraph* asymPlot = new TGraph();
  // for (int bin = 0; bin < eBins; ++bin){
  
  //   auto eBinFit = new TF1("eBinFit", mywiggleFitClass, &wiggleFitClass::fiveParFunc, 0, histMaxTime, numTmethodParams);
     
  //   // Set fit parameters
  //   eBinFit->SetParameter(0, 1);
  //   eBinFit->SetParameter(1, defaultLifetime/1000);
  //   eBinFit->SetParameter(2, 0.4);
  //   eBinFit->SetParameter(3, 0);
  //   eBinFit->SetParameter(4, pi);
     
  //   // Choose good starting guess normalisation N by comparing integral of function and histogram     
  //   double normGuess = toyEbins[0][bin]->Integral(toyEbins[0][bin]->GetXaxis()->FindBin(0.0), toyEbins[0][bin]->GetXaxis()->FindBin(95000), "WIDTH") / eBinFit->Integral(0,95000);
  //   eBinFit->SetParameter(0, normGuess);  
     
  //   // Fit function
  //   toyEbins[0][bin]->Fit(eBinFit,"QNREML");
     
  //   cout << bin*100 << " - " << (bin+1)*100 << " MeV: A = " << eBinFit->GetParameter(2) << endl;
     
  //   double x = (bin*100 + (bin+1)*100)/2.0;
  //   double y = eBinFit->GetParameter(2);
     
  //   asymPlot->SetPoint(bin,x,y);
     
  // }
   
  // //asymPlot->Draw();
  // TF1* asymFit = new TF1("asymFit","pol5");
  // asymPlot->Fit(asymFit);
  // //asymFit->Draw("SAME");

  // for (int bin = 0; bin < eBins; ++bin){
  //   double x = (bin*100 + (bin+1)*100)/2.0;
  //   cout << bin*100 << " - " << (bin+1)*100 << " MeV: A = " << muonA->Eval(x) << ", " << asymFit->Eval(x) << endl;
  // }
      
  /////////////////////////////////////////////////////////////////////////////////////     
  
  //////////////////////////////////////////////////////////////////////////////////////////////
  // Fit histograms
  ////////////////////////////////////////////////////////////////////////////////////////////// 

  // for (int iter = 0; iter < totalIters; ++iter){

  //   cout << "\nIteration " << iter << endl;

  //   // Fit T-weighted histogram
  //   TH1F* hWiggleT = (TH1F*)toyFiveParamTwgtHists[iter]->Clone("hWiggleT");  
  //   auto fitTmethod = new TF1("fitTmethod", fiveParamFit, &wiggleFitClass::fiveParFunc, fitStart, fitEnd, numTmethodParams);
  //   fitTmethod->SetNpx(10000);
  
  //   // Set fit parameters
  //   fitTmethod->SetParameter(0, 1);
  //   fitTmethod->SetParameter(1, defaultLifetime/1000);
  //   fitTmethod->SetParameter(2, 0.4);
  //   fitTmethod->SetParameter(3, 0);
  //   fitTmethod->SetParameter(4, pi);
  
  //   // Choose good starting guess normalisation N by comparing integral of function and histogram     
  //   double normGuessT = hWiggleT->Integral(hWiggleT->GetXaxis()->FindBin(0.0), hWiggleT->GetXaxis()->FindBin(95000), "WIDTH") / fitTmethod->Integral(0,95000);
  //   fitTmethod->SetParameter(0, normGuessT);  

  //   // Fit function
  //   hWiggleT->Fit(fitTmethod,"QNREM");

  //   cout << "\nT-method" << endl;
  //   cout << "--------" << endl;
  //   cout << "N = " << fitTmethod->GetParameter(0) << endl;
  //   cout << "tau = " << fitTmethod->GetParameter(1) << endl;
  //   cout << "A = " << fitTmethod->GetParameter(2) << endl;
  //   cout << "R = " << fitTmethod->GetParameter(3) << endl;
  //   cout << "phi = " << fitTmethod->GetParameter(4) << endl;
  
  //   double TmethodPvalue = fitTmethod->GetProb();
  
  //   std::cout << "T-method fit p-value is: " << TmethodPvalue << " chi2/ndf: " << fitTmethod->GetChisquare()/fitTmethod->GetNDF() << std::endl;
  //   cout << "T-method R error: " << fitTmethod->GetParError(3) << " ppm "  << endl;

  //   ////////////////////////////////////////////////////////////////////////////////////////////// 

  //   // Fit A-weighted histogram
  //   TH1F* hWiggleA = (TH1F*)toyFiveParamAwgtHists[iter]->Clone("hWiggleA");
  //   auto fitAmethod = new TF1("fitAmethod", fiveParamFit, &wiggleFitClass::fiveParFunc, fitStart, fitEnd, numTmethodParams);
  //   fitAmethod->SetNpx(10000);
  //   // Set fit parameters
  //   fitAmethod->SetParameter(0, 1);
  //   fitAmethod->SetParameter(1, defaultLifetime/1000);
  //   fitAmethod->SetParameter(2, 0.4);
  //   fitAmethod->SetParameter(3, 0);
  //   fitAmethod->SetParameter(4, pi);
  
  //   // Choose good starting guess normalisation N by comparing integral of function and histogram     
  //   double normGuessA = hWiggleA->Integral(hWiggleA->GetXaxis()->FindBin(0.0), hWiggleA->GetXaxis()->FindBin(95000), "WIDTH") / fitAmethod->Integral(0,95000);
  //   fitAmethod->SetParameter(0, normGuessA);  
  
  //   // Fit function
  //   hWiggleA->Fit(fitAmethod,"QNREM");
  
  //   cout << "\nA-method" << endl;
  //   cout << "--------" << endl;
  //   cout << "N = " << fitAmethod->GetParameter(0) << endl;
  //   cout << "tau = " << fitAmethod->GetParameter(1) << endl;
  //   cout << "A = " << fitAmethod->GetParameter(2) << endl;
  //   cout << "R = " << fitAmethod->GetParameter(3) << endl;
  //   cout << "phi = " << fitAmethod->GetParameter(4) << endl;
  
  //   double AmethodPvalue = fitAmethod->GetProb();
  
  //   std::cout << "A-method fit p-value is: " << AmethodPvalue << " chi2/ndf: " << fitAmethod->GetChisquare()/fitAmethod->GetNDF() << std::endl;
  //   cout << "A-method R error: " << fitAmethod->GetParError(3) << " ppm "  << endl;

  //   ////////////////////////////////////////////////////////////////////////////////////////////// 

  //   // Fit ratio histogram
  //   TH1F* hWiggleR = (TH1F*)toyRatioHists[iter]->Clone("hWiggleR");  
  //   auto ratioFit = new wiggleFitClass(numRatioParams);
  //   auto fitRatio = new TF1("fitRatio", ratioFit, &wiggleFitClass::threeParFunc, fitStart, ratioFitEnd, numRatioParams);
  //   fitRatio->SetNpx(10000);

  //   // Set fit parameters
  //   fitRatio->SetParameter(0, 0.4);
  //   fitRatio->SetParameter(1, 0);
  //   fitRatio->SetParameter(2, pi);
  
  //   // Fit function
  //   hWiggleR->Fit(fitRatio,"QNREM");

  //   cout << "\nRatio-method" << endl;
  //   cout << "--------" << endl;
  //   cout << "A = " << fitRatio->GetParameter(0) << endl;
  //   cout << "R = " << fitRatio->GetParameter(1) << endl;
  //   cout << "phi = " << fitRatio->GetParameter(2) << endl;
  
  //   double ratioPvalue = fitRatio->GetProb();
  
  //   std::cout << "Ratio-method fit p-value is: " << ratioPvalue << " chi2/ndf: " << fitRatio->GetChisquare()/fitRatio->GetNDF() << std::endl;
  //   cout << "Ratio-method R error: " << fitRatio->GetParError(1) << " ppm "  << endl;

  //   // Fit Q-weighted histogram
  //   TH1F* hWiggleQ = (TH1F*)toyFiveParamQHists[iter]->Clone("hWiggleQ");  
  //   auto fitQmethod = new TF1("fitQmethod", fiveParamFit, &wiggleFitClass::fiveParFunc, fitStart, QFitEnd, numTmethodParams);
  //   fitQmethod->SetNpx(10000);
  
  //   // Set fit parameters
  //   fitQmethod->SetParameter(0, 1);
  //   fitQmethod->SetParameter(1, defaultLifetime/1000);
  //   fitQmethod->SetParameter(2, 0.4);
  //   fitQmethod->SetParameter(3, 0);
  //   fitQmethod->SetParameter(4, pi);
  
  //   // Choose good starting guess normalisation N by comparing integral of function and histogram     
  //   double normGuessQ = hWiggleQ->Integral(hWiggleQ->GetXaxis()->FindBin(0.0), hWiggleQ->GetXaxis()->FindBin(95000), "WIDTH") / fitQmethod->Integral(0,95000);
  //   fitQmethod->SetParameter(0, normGuessQ);  

  //   // Fit function
  //   hWiggleQ->Fit(fitQmethod,"QNREM");

  //   cout << "\nQ-method" << endl;
  //   cout << "--------" << endl;
  //   cout << "N = " << fitQmethod->GetParameter(0) << endl;
  //   cout << "tau = " << fitQmethod->GetParameter(1) << endl;
  //   cout << "A = " << fitQmethod->GetParameter(2) << endl;
  //   cout << "R = " << fitQmethod->GetParameter(3) << endl;
  //   cout << "phi = " << fitQmethod->GetParameter(4) << endl;
  
  //   double QmethodPvalue = fitQmethod->GetProb();
  
  //   std::cout << "Q-method fit p-value is: " << QmethodPvalue << " chi2/ndf: " << fitQmethod->GetChisquare()/fitQmethod->GetNDF() << std::endl;
  //   cout << "Q-method R error: " << fitQmethod->GetParError(3) << " ppm "  << endl;


  //  }

  /////////////////////////////////////////////////////////////////////////////////////

  outputFile->Write();

  delete outputFile;

  /////////////////////////////////////////////////////////////////////////////////////    

  return 1;
  
}
 
