// 3-31-20: Macro for plots related to lost muons and accidental background subtraction. Don't remember exactly what it makes, should rerun to find out.

#include <iostream>
#include <fstream>
#include <string>
#include <functional>
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

#include "../ratioMacroHeaders/ratioAnalysisDefs.hh" // don't change these paths since this file is run on the grid (whereas other files can see my root include directory and the ../ can be removed)


int LMsubtractAccBG()
{

  int nBins = int(approxMaxTime/defaultBinWidth);
  double histMaxTime = nBins*defaultBinWidth;
  double histMinTime = 0;


  // TFile *inputFile = TFile::Open("/gm2/data/users/nkinnaird/Ratio/Endgame-FinalProduction/LostMuons/lostMuons-Endgame-spectra-check-2.root");
  // TH2F* deltaT12_vs_timeInFill = (TH2F*) inputFile->Get("topDir/Triples/triple_deltaT12_vs_timeInFill");
 
  TFile *inputFile = TFile::Open("/gm2/data/users/nkinnaird/Ratio/Endgame-FinalProduction/LostMuons/lostMuons-Endgame-spectra-check-3.root");
  TH2F* deltaT12_vs_timeInFill = (TH2F*) inputFile->Get("topDir/Info/TriplesInfo/triple_deltaT12_vs_timeInFill");
  
  auto inHist_canv = new TCanvas("inHist_canv","inHist_canv",0,0,800,600);
  inHist_canv->SetLogz();
  deltaT12_vs_timeInFill->Draw("COLZ");

  TH2F* accidentalBackground = (TH2F*) deltaT12_vs_timeInFill->Clone();
  accidentalBackground->Reset();

  // subtract off times below some cut off within which to grab the background
  double bgCutoffLow = 2; // ns
  double bgCutoffHigh = 4; // ns
  int yBinCutoffLow = accidentalBackground->GetYaxis()->FindBin(bgCutoffLow);
  int yBinCutoffHigh = accidentalBackground->GetYaxis()->FindBin(bgCutoffHigh) - 1;

  for (int yBin = yBinCutoffLow; yBin <= yBinCutoffHigh; ++yBin){
  	for (int xBin = 1; xBin <= accidentalBackground->GetXaxis()->GetNbins(); ++xBin){
  		accidentalBackground->SetBinContent(xBin, yBin, deltaT12_vs_timeInFill->GetBinContent(xBin, yBin));
	  }
  }

  auto bg_canv = new TCanvas("bg_canv","bg_canv",800,0,800,600);
  bg_canv->SetLogz();
  accidentalBackground->Draw("COLZ");

  // average accidental background 

  TH2F* accidentalBackground_average = (TH2F*) accidentalBackground->Clone();
  int yCounts = yBinCutoffHigh - yBinCutoffLow + 1;

  for (int xBin = 1; xBin <= accidentalBackground_average->GetXaxis()->GetNbins(); ++xBin){

  	double ySum = 0, yAverage = 0;

  		for (int yBin = yBinCutoffLow; yBin <= yBinCutoffHigh; ++yBin) ySum += accidentalBackground_average->GetBinContent(xBin, yBin);
  		yAverage = ySum/yCounts;
  		for (int yBin = 1; yBin <= accidentalBackground_average->GetYaxis()->GetNbins(); ++yBin) accidentalBackground_average->SetBinContent(xBin, yBin, yAverage);
  }

  auto bg_avg_canv = new TCanvas("bg_avg_canv","bg_avg_canv",0,600,800,600);
  bg_avg_canv->SetLogz();
  accidentalBackground_average->Draw("COLZ");

  // TH1F* accBg_avg_1D = new TH1F("Accidental_Background_Average_1D", "Accidental Background Average; time (ns); Events", nBins, histMinTime, histMaxTime);
  TH1F* accBg_avg_1D = new TH1F("Accidental_Background_Average_1D", "Accidental Background Average; time (ns); Events", 700, 0, 700000);

  int numYbinsInLossHist = accidentalBackground->GetYaxis()->FindBin(7.5) - accidentalBackground->GetYaxis()->FindBin(5);
  cout << "numYbinsInLossHist: " << numYbinsInLossHist << endl;

  for (int xBin = 1; xBin <= accidentalBackground_average->GetXaxis()->GetNbins(); ++xBin){
  	double stripeTime = accidentalBackground_average->GetXaxis()->GetBinCenter(xBin);
  	double stripeContent = accidentalBackground_average->GetBinContent(xBin, 1);

  	accBg_avg_1D->Fill(stripeTime, stripeContent*numYbinsInLossHist);
  }

  auto bg_1d_canv = new TCanvas("bg_1d_canv","bg_1d_canv",800,600,800,600);
  // ((TH1F*) inputFile->Get("topDir/Triples/Losses/Iter0/triple_losses_spectra"))->Draw("HIST");
  ((TH1F*) inputFile->Get("topDir/AdditionalCuts/Triples/Losses/triple_losses_spectra"))->Draw("HIST");
  accBg_avg_1D->SetLineColor(2);
  accBg_avg_1D->Draw("HISTSAME");


  // subtract off average from main histogram

  TH2F* deltaT12_vs_timeInFill_subtracted = (TH2F*) deltaT12_vs_timeInFill->Clone();

  deltaT12_vs_timeInFill_subtracted->Add(accidentalBackground_average, -1);

  for (int yBin = 1; yBin <= deltaT12_vs_timeInFill_subtracted->GetYaxis()->GetNbins(); ++yBin){
  	for (int xBin = 1; xBin <= deltaT12_vs_timeInFill_subtracted->GetXaxis()->GetNbins(); ++xBin){
  		if(deltaT12_vs_timeInFill_subtracted->GetBinContent(xBin, yBin) < 1) deltaT12_vs_timeInFill_subtracted->SetBinContent(xBin, yBin, 0);
      // if(deltaT12_vs_timeInFill_subtracted->GetBinContent(xBin, yBin) >= 1) deltaT12_vs_timeInFill_subtracted->SetBinContent(xBin, yBin, 0);
	  }
  }

  auto sub_canv = new TCanvas("sub_canv","sub_canv",800,600,800,600);
  sub_canv->SetLogz();
  deltaT12_vs_timeInFill_subtracted->Draw("COLZ");

return 1;
}
