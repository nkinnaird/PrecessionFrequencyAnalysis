// 3-31-20: Old macro for saving some fit functions and lost muon spectra into a separate use file by someone else - probably Alex, don't really remember.

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2D.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TImage.h>
#include <sstream>
#include <TNtuple.h>
#include <TPaveStats.h>
#include <TMatrixD.h>

using namespace std;


int SaveFunction(std::string filePath)
{

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

   TFile* outputFile = new TFile("savedFunction.root","RECREATE");

   auto lostMuonDir = outputFile->mkdir("LostMuons");
   auto TmethodDir = outputFile->mkdir("TMethod");
   auto ratioMethodDir = outputFile->mkdir("Ratio");

/////////////////////////////////////////////////////////////////////////////////////

  lostMuonDir->cd();

  TH1F* lostMuonIntegral = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/LostMuons/Losses_triple_integral");
  lostMuonIntegral->Write();

/////////////////////////////////////////////////////////////////////////////////////

 	TmethodDir->cd();

	TF1* TMethod_fit = (TF1*) ((TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/TMethod/allTimesAdded_TMethod"))->GetFunction("TmethodFitFunc");
	TMethod_fit->Write();
	
	inputFile->Get("topDir/FitPasses/FitPass0/addedDir/TMethod/CorrelationMatrix")->Write("CorrelationMatrix");


/////////////////////////////////////////////////////////////////////////////////////

	ratioMethodDir->cd();

  TF1* fullRatio_fit = (TF1*) ((TGraphErrors*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/FullRatio/Added_Times_Full_Ratio_Graph"))->GetFunction("fullRatioFitFunc");
  fullRatio_fit->Write();

  inputFile->Get("topDir/FitPasses/FitPass0/addedDir/FullRatio/FullRatioFitCorrMatrix")->Write("CorrelationMatrix");


	return 1;

}
