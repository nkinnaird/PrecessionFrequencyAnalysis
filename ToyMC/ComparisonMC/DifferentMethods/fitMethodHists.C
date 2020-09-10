R__LOAD_LIBRARY(/cvmfs/gm2.opensciencegrid.org/prod/g-2/gm2util/v9_21_06/slf6.x86_64.e15.prof/lib/libgm2util_blinders.so)

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
#include <TFitResult.h>
#include <TSpectrum.h>
#include <TText.h>
#include <TLatex.h>
#include <TPaveText.h>

#include <time.h>

#include "gm2util/blinders/Blinders.hh"

#include "ratioAnalysisDefs.hh"
#include "ratioAnalysisConfig.hh"
#include "fiveParamFit.hh"
#include "threeParameterRatioFit.hh"

#include "residualPlots.hh"

using namespace std;

// ******************************************************************* //                                                                           
// ********** Global variable definition and initialisation ********** //
// ******************************************************************* //

// Total number of 100 MeV enegry bins (3.1 GeV max)
int eBins = 31;

// Fit start time
static const double toyFitStart = 5000; // 30000; // ns // need an earlier fit start time for more statistics
// Fit end time
// static const double toyFitEnd = 600000; // ns
// static const double QFitEnd = 206000; // ns // tuned so that the ratio of Q method R error to T method R error is about 1.547, which is what it is for the Run 1 data (needs to be tuned depending on the fit start and end times)
static const double toyFitEnd = 350000; // ns
static const double QFitEnd = 185000; // ns // tuned so that the ratio of Q method R error to T method R error is about 1.547, which is what it is for the Run 1 data (needs to be tuned depending on the fit start and end times)



bool saveFirstFits = true; // turn to false before fitting many files


//***********************************************************************//
//*************************** Main Program ******************************//
//***********************************************************************//

int fitMethodHists(std::string filePath)
{
   gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  // pull in input file
  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  // create output file that will hold plots
  TFile* outputFile = new TFile("methodFits.root","RECREATE");

  // make top directory for output file
  auto topDir = outputFile->mkdir("topDir");

  setDataset(1, 12); // set dataset case to 60h - for 'starting' parameters into fits


  /////////////////////////////////////////////////////////////////////////////////////

  bool doR = false, doQ = true, doE = false;
  
  int totalIters = (*(TVectorD*) inputFile->Get("Iters"))[0]; // total iterations in generated histograms )
  // int totalIters = 1;

  auto TMethod_dir = topDir->mkdir("T");
  TMethod_dir->cd();

  auto TmethodRPull = new TH1F("TmethodRPull", "TmethodRPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10);
  auto Tmethodchi2 = new TH1F("Tmethodchi2", "Tmethodchi2; #chi^{2}/d.o.f.; Events", 50, 0.5, 1.5);
  auto TmethodPvalues = new TH1F("TmethodPvalues", "TmethodPvalues; P-values; Events", 50, 0, 1);

  TNtuple* TMethodFitValues = new TNtuple("TMethodFitValues", "TMethodFitValues", "N:N_err:tau:tau_err:A:A_err:R:R_err:phi:phi_err");

  auto AMethod_dir = topDir->mkdir("A");
  AMethod_dir->cd();
 
  auto AmethodRPull = new TH1F("AmethodRPull", "AmethodRPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10);
  auto Amethodchi2 = new TH1F("Amethodchi2", "Amethodchi2; #chi^{2}/d.o.f.; Events", 50, 0.5, 1.5);
  auto AmethodPvalues = new TH1F("AmethodPvalues", "AmethodPvalues; P-values; Events", 50, 0, 1);

  TNtuple* AMethodFitValues = new TNtuple("AMethodFitValues", "AMethodFitValues", "N:N_err:tau:tau_err:A:A_err:R:R_err:phi:phi_err");
  
  auto RMethod_dir = topDir->mkdir("R");
  RMethod_dir->cd();

  auto RmethodRPull = new TH1F("RmethodRPull", "RmethodRPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10); 
  auto Rmethodchi2 = new TH1F("Rmethodchi2", "Rmethodchi2; #chi^{2}/d.o.f.; Events", 50, 0.5, 1.5);
  auto RmethodPvalues = new TH1F("RmethodPvalues", "RmethodPvalues; P-values; Events", 50, 0, 1);
  
  TNtuple* ratioFitValues = new TNtuple("ratioFitValues", "ratioFitValues_VW", "A:A_err:R:R_err:phi:phi_err");

  auto QMethod_dir = topDir->mkdir("Q");
  QMethod_dir->cd();

  auto QmethodRPull = new TH1F("QmethodRPull", "QmethodRPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10); 
  auto Qmethodchi2 = new TH1F("Qmethodchi2", "Qmethodchi2; #chi^{2}/d.o.f.; Events", 50, 0.5, 1.5);
  auto QmethodPvalues = new TH1F("QmethodPvalues", "QmethodPvalues; P-values; Events", 50, 0, 1);

  TNtuple* QMethodFitValues = new TNtuple("QMethodFitValues", "QMethodFitValues", "N:N_err:tau:tau_err:A:A_err:R:R_err:phi:phi_err");

  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////

  blinding::Blinders::fitType ftype = blinding::Blinders::kOmega_a;
  blinding::Blinders* myBlinder = new blinding::Blinders(ftype); // no blinding for ToyMC

  FiveParamFit fiveParamFitToyClass(myBlinder);
  ThreeParameterRatioFit ratioFitToyClass(myBlinder);

/////////////////////////////////////////////////////////////////////////////////////


  // Muon asymmetry function (lab frame)                                                                                                                                     
  double Emax = 3100; // Maximum energy [MeV]                                                                                                                                
  TF1 *muonA = new TF1("muonA", "(-8.0*(x/[0])*(x/[0]) + (x/[0]) + 1.1)/(4.0*(x/[0])*(x/[0]) - 5.0*(x/[0]) - 5.0)", 0, Emax);
  muonA->SetParameter(0,Emax);

 
  for (int iter = 0; iter < totalIters; ++iter){ // loop over iterations / different random seeds
   
    cout << "Iter: " << iter << endl;

    double RfromT = 0;
    double TmethodError = 0;
    double RfromA = 0;
    double AmethodError = 0;
    double RfromR = 0;
    double ratioError = 0;
    double RfromQ = 0;
    double QmethodError = 0;
    double RfromE = 0;
    double EbinnedError = 0;
    
    /////////////////////////////////////////////////////////////////////////////////////     
    
    //////////////////////////////////////////////////////////////////////////////////////////////
    // Fit histograms
    ////////////////////////////////////////////////////////////////////////////////////////////// 
      
    // Fit T-weighted histogram
    TH1F* hWiggleT = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_5_Param_Hist_T", iter)))->Clone();
    hWiggleT->SetName("hWiggleT");
    

    fiveParamFitToyClass.setFitRange(toyFitStart, toyFitEnd);
    fiveParamFitToyClass.fiveParameterFitMethod(hWiggleT);
    auto fitTmethod = hWiggleT->GetFunction("fiveParamFit");


    double Tchi2 = fitTmethod->GetChisquare()/fitTmethod->GetNDF();
    double TmethodPvalue = fitTmethod->GetProb();
      
    if (iter == 0){
      cout << "\nT-method" << endl;
      cout << "--------" << endl;
      cout << "N = " << fitTmethod->GetParameter(0) << endl;
      cout << "tau = " << fitTmethod->GetParameter(1) << endl;
      cout << "A = " << fitTmethod->GetParameter(2) << endl;
      cout << "R = " << fitTmethod->GetParameter(3) << endl;
      cout << "phi = " << fitTmethod->GetParameter(4) << endl;
      
      cout << "T-method fit p-value is: " << TmethodPvalue << " chi2/ndf: " << fitTmethod->GetChisquare()/fitTmethod->GetNDF() << std::endl;
      cout << "T-method R error: " << fitTmethod->GetParError(3) << " ppm "  << endl;

      if(saveFirstFits){
        TMethod_dir->cd();
        hWiggleT->Write();

        // TH1F* hWiggleT_FFT = doFFT(hWiggleT);

          ResidualPlots residualPlotsClass;
          residualPlotsClass.makeResidualPlots("TMethod", hWiggleT, "residual_fitRange", toyFitStart, toyFitEnd);
      }
    }

    RfromT = fitTmethod->GetParameter(3);
    TmethodError = fitTmethod->GetParError(3);
     

    // Calculate T-method fit pulls
    TmethodRPull->Fill( (RfromT - 0) / TmethodError );
    Tmethodchi2->Fill(Tchi2);
    TmethodPvalues->Fill(TmethodPvalue);

    TMethodFitValues->Fill(fitTmethod->GetParameter(0), fitTmethod->GetParError(0),
                           fitTmethod->GetParameter(1), fitTmethod->GetParError(1),
                           fitTmethod->GetParameter(2), fitTmethod->GetParError(2),
                           fitTmethod->GetParameter(3), fitTmethod->GetParError(3),
                           fitTmethod->GetParameter(4), fitTmethod->GetParError(4));

    delete hWiggleT;


    ////////////////////////////////////////////////////////////////////////////////////////////// 
      
    // Fit A-weighted histogram
    TH1F* hWiggleA = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_5_Param_Hist_A", iter)))->Clone();
    

    fiveParamFitToyClass.setFitRange(toyFitStart, toyFitEnd);
    fiveParamFitToyClass.fiveParameterFitMethod(hWiggleA);
    auto fitAmethod = hWiggleA->GetFunction("fiveParamFit");

    double Achi2 = fitAmethod->GetChisquare()/fitAmethod->GetNDF();
    double AmethodPvalue = fitAmethod->GetProb();
      
    if (iter == 0){
      cout << "\nA-method" << endl;
      cout << "--------" << endl;
      cout << "N = " << fitAmethod->GetParameter(0) << endl;
      cout << "tau = " << fitAmethod->GetParameter(1) << endl;
      cout << "A = " << fitAmethod->GetParameter(2) << endl;
      cout << "R = " << fitAmethod->GetParameter(3) << endl;
      cout << "phi = " << fitAmethod->GetParameter(4) << endl;
	
      cout << "A-method fit p-value is: " << AmethodPvalue << " chi2/ndf: " << fitAmethod->GetChisquare()/fitAmethod->GetNDF() << std::endl;
      cout << "A-method R error: " << fitAmethod->GetParError(3) << " ppm "  << endl;

      if(saveFirstFits){
        AMethod_dir->cd();
        hWiggleA->Write();
      }
    }
      
    RfromA = fitAmethod->GetParameter(3);
    AmethodError = fitAmethod->GetParError(3);


    // Calculate A-method fit pulls
    AmethodRPull->Fill( (RfromA - 0) / AmethodError );
    Amethodchi2->Fill(Achi2);
    AmethodPvalues->Fill(AmethodPvalue);

    AMethodFitValues->Fill(fitAmethod->GetParameter(0), fitAmethod->GetParError(0),
                           fitAmethod->GetParameter(1), fitAmethod->GetParError(1),
                           fitAmethod->GetParameter(2), fitAmethod->GetParError(2),
                           fitAmethod->GetParameter(3), fitAmethod->GetParError(3),
                           fitAmethod->GetParameter(4), fitAmethod->GetParError(4));

    delete hWiggleA;

    ////////////////////////////////////////////////////////////////////////////////////////////// 
      
    // Fit ratio histogram

    if(doR)
    {
      TH1F* toyUHist = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_U_Hist", iter)))->Clone();
      TH1F* toyVHist = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_V_Hist", iter)))->Clone();

      TH1F* toyNumHist = (TH1F*) toyVHist->Clone("numeratorHist");
      toyNumHist->Add(toyUHist, -1);

      TH1F* toyDenomHist = (TH1F*) toyVHist->Clone("denominatorHist");
      toyDenomHist->Add(toyUHist);

      TGraphErrors* toyRatioGraph = ratioFitToyClass.createRatioGraph(toyNumHist, toyDenomHist);
      
      ratioFitToyClass.setFitRange(toyFitStart, toyFitEnd);    
      ratioFitToyClass.fitMethod(toyRatioGraph);
      auto fitRatio = toyRatioGraph->GetFunction("threeParamRatioFit");


      double ratiochi2 = fitRatio->GetChisquare()/fitRatio->GetNDF();
      double ratioPvalue = fitRatio->GetProb();


      if (iter == 0){
        cout << "\nRatio-method" << endl;
        cout << "--------" << endl;
        cout << "A = " << fitRatio->GetParameter(0) << endl;
        cout << "R = " << fitRatio->GetParameter(1) << endl;
        cout << "phi = " << fitRatio->GetParameter(2) << endl;
    
        double ratioPvalue = fitRatio->GetProb();
    
        std::cout << "Ratio-method fit p-value is: " << ratioPvalue << " chi2/ndf: " << fitRatio->GetChisquare()/fitRatio->GetNDF() << std::endl;
        cout << "Ratio-method R error: " << fitRatio->GetParError(1) << " ppm "  << endl;

        if(saveFirstFits){
          RMethod_dir->cd();
          toyRatioGraph->Write("ToyRatioGraph");
        }
      }

      RfromR = fitRatio->GetParameter(1);
      ratioError = fitRatio->GetParError(1);

      // Calculate ratio method fit pulls
      RmethodRPull->Fill( (RfromR - 0) / ratioError );
      Rmethodchi2->Fill(ratiochi2);
      RmethodPvalues->Fill(ratioPvalue);

      ratioFitValues->Fill(fitRatio->GetParameter(0), fitRatio->GetParError(0),
                           fitRatio->GetParameter(1), fitRatio->GetParError(1),
                           fitRatio->GetParameter(2), fitRatio->GetParError(2));


      delete toyUHist;
      delete toyVHist;
      delete toyNumHist;
      delete toyDenomHist;
      delete toyRatioGraph; 
    }


/////////////////////////////////////////////////////////////////////////////////////

    // Fit Q-method histogram
    if (doQ){
      TH1F* hWiggleQ = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_5_Param_Hist_Q", iter)))->Clone();
      hWiggleQ->SetName("hWiggleQ");

      fiveParamFitToyClass.setFitRange(toyFitStart, QFitEnd);
      fiveParamFitToyClass.fiveParameterFitMethod(hWiggleQ);
      auto fitQmethod = hWiggleQ->GetFunction("fiveParamFit");


      double Qchi2 = fitQmethod->GetChisquare()/fitQmethod->GetNDF();
      double QmethodPvalue = fitQmethod->GetProb();
      
      if (iter == 0){
      	cout << "\nQ-method" << endl;
      	cout << "--------" << endl;
      	cout << "N = " << fitQmethod->GetParameter(0) << endl;
      	cout << "tau = " << fitQmethod->GetParameter(1) << endl;
      	cout << "A = " << fitQmethod->GetParameter(2) << endl;
      	cout << "R = " << fitQmethod->GetParameter(3) << endl;
      	cout << "phi = " << fitQmethod->GetParameter(4) << endl;
            
      	std::cout << "Q-method fit p-value is: " << QmethodPvalue << " chi2/ndf: " << fitQmethod->GetChisquare()/fitQmethod->GetNDF() << std::endl;
      	cout << "Q-method R error: " << fitQmethod->GetParError(3) << " ppm "  << endl;

        if(saveFirstFits){
          QMethod_dir->cd();
          hWiggleQ->Write();
        }
      }

      RfromQ = fitQmethod->GetParameter(3);
      QmethodError = fitQmethod->GetParError(3);
     
      // Calculate Q-method fit pulls
      QmethodRPull->Fill( (RfromQ - 0) / QmethodError );
      Qmethodchi2->Fill(Qchi2);
      QmethodPvalues->Fill(QmethodPvalue);


      QMethodFitValues->Fill(fitQmethod->GetParameter(0), fitQmethod->GetParError(0),
                             fitQmethod->GetParameter(1), fitQmethod->GetParError(1),
                             fitQmethod->GetParameter(2), fitQmethod->GetParError(2),
                             fitQmethod->GetParameter(3), fitQmethod->GetParError(3),
                             fitQmethod->GetParameter(4), fitQmethod->GetParError(4));

      delete hWiggleQ;
    }

    /////////////////////////////////////////////////////////////////////////////////////     

    // Fit E-binned histograms and average
    // not doing this at the moment so some of this may need to be updated
    if (doE){

      double avg = 0;
      double avgErr = 0;
      double avgNum = 0;
      double avgDen = 0;
      for (int bin = 0; bin < eBins; ++bin){

      	double eLow = double(bin*100);
      	double eHigh = double((bin+1)*100); 
      	double energy = (eLow+eHigh)/2.0;
      	double A = muonA->Eval(energy);
            
      	TH1F* toyEbin = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/%d - %d MeV",iter,bin*100,(bin+1)*100)))->Clone();      	

        fiveParamFitToyClass.setFitRange(toyFitStart, toyFitEnd); // need to be updated in some way
        fiveParamFitToyClass.fiveParameterFitMethod(toyEbin);
        auto eBinFit = toyEbin->GetFunction("fiveParamFit");


      	// Sum towards E-binned average
      	double par = eBinFit->GetParameter(3);
      	double parErr = eBinFit->GetParError(3);
      	// Set E-binned threshold
      	if (eLow > 500 && eLow < 2900) {
      	  avgNum += par/pow(parErr,2);
      	  avgDen += 1/pow(parErr,2);
      	}

        delete toyEbin;
           
      }
   
      // Calculate E-binned average
      avg = avgNum/avgDen;
      avgErr = 1/sqrt(avgDen);
      RfromE = avg;
      EbinnedError = avgErr;

      if (iter == 0){
      	cout << "\nE-binned" << endl;
      	cout << "--------" << endl;
      	cout << "R = " << RfromE << endl;
      	cout << "T-binned R error: " << EbinnedError << " ppm "  << endl;
      
        // if(saveFirstFits){
          // TMethod_dir->cd();
          // hWiggleT->Write();
        // }
      }
    }

/////////////////////////////////////////////////////////////////////////////////////

  } // end loop over iterations


/////////////////////////////////////////////////////////////////////////////


  outputFile->Write();
  delete outputFile;

  return 1;

}
