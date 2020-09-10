// 3-31-20: Macro for plots for some integrated lost muons thing - don't remember exactly but it doesn't seem to do as much as the file title might suggest. I think it just gives the fraction of muons lost into the fill.

#include <iostream>
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

#include "ratioAnalysisDefs.hh"


int LostMuonSystematicEstimator(std::string filePath)
{

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

   TFile* outputFile = new TFile("lostMuonFraction.root","RECREATE");

   	TF1* fitFunction_TMethod = (TF1*) ((TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/TMethod/allTimesAdded_TMethod"))->GetFunction("TmethodFitFunc");   	

   	double fittedN0 = fitFunction_TMethod->GetParameter(0);
   	double fittedLifetime = fitFunction_TMethod->GetParameter(1);
   	double kappa_loss = fitFunction_TMethod->GetParameter(13);

   	cout << "fittedN0 : " << fittedN0 << endl;
   	cout << "fittedLifetime : " << fittedLifetime << endl;
   	cout << "kappa_loss : " << kappa_loss << endl;
   	
   	cout << "Don't forget about any arbitrary scaling factors in kappa loss." << endl;
	kappa_loss *= 1e-10;

   	TH1F* lostMuonsHist = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/LostMuons/Losses_triple");
   	TH1F* lostMuonsIntegralHist = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/LostMuons/Losses_triple_integral");

      lostMuonsHist->Write();
      lostMuonsIntegralHist->Write();

   	TH1F* totalHist = (TH1F*) lostMuonsHist->Clone("total"); // clone to get same bin structure 
   	totalHist->SetTitle("Total Muons");
      totalHist->SetName("totalMuons");
   	totalHist->Reset("ICES");

   	TH1F* fractionLost = (TH1F*) lostMuonsHist->Clone("fractionLost"); // clone to get same bin structure 
   	fractionLost->SetTitle("Fraction of Muons Lost");
      fractionLost->SetName("fractionLost");
   	fractionLost->Reset("ICES");

   	// TH1F* fractionLostIntegral = (TH1F*) lostMuonsHist->Clone("fractionLostIntegral"); // clone to get same bin structure 
   	// fractionLostIntegral->SetTitle("Integral of Fraction of Muons Lost");
   	// fractionLostIntegral->Reset("ICES");

   	// double lostIntegral = 0; for the fraction lost somehow need to combine with kappa_loss to get right number

   	for (int bin = 1; bin <= fractionLost->GetNbinsX(); ++bin)
   	{
   		double time = fractionLost->GetBinCenter(bin);
   		double totalInBin = fittedN0 * exp(-time/fittedLifetime) * (1 - kappa_loss * lostMuonsIntegralHist->GetBinContent(bin));
   		double lostInBin = lostMuonsHist->GetBinContent(bin);
   		double fractionInBin = lostInBin / totalInBin;

   		// if(time > 30000) lostIntegral += lostInBin;
   		
   		totalHist->SetBinContent(bin, totalInBin);
   		fractionLost->SetBinContent(bin, fractionInBin);
   		// fractionLostIntegral->SetBinContent(bin, lostIntegral / fittedN0);
   	}

		  // TF1* expFunc = new TF1("expFunc", "[0]*exp(-x/[1])", 0, histMaxTime);
		  // expFunc->SetParameter(0, totalHist->GetMaximum()); 
		  // expFunc->SetParameter(1, defaultLifetime);

		  // totalHist->Fit(expFunc, "QRLI");
		  // totalHist->GetFunction("expFunc")->SetLineColor(2);

      totalHist->Write();
      fractionLost->Write();

   	TCanvas* totalCanvas = new TCanvas("totalCanvas", "canvas",200,10,1000,800);
   		totalHist->Draw("hist");
   	
   	TCanvas* fractionCanvas = new TCanvas("fractionCanvas", "canvas",200,10,1000,800);
		   fractionLost->Draw("hist");

   	// TCanvas* fractionIntegralCanvas = new TCanvas("fractionIntegralCanvas", "canvas",200,10,1000,800);
		// fractionLostIntegral->Draw("hist");

      // delete outputFile;

  return 1;

}
