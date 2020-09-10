// This macro is intended to compare pileup corrected times above 3.5 GeV. The incoming histogram needs to have the correct energy threshold set.
// ratioMacro.C should be ran in order to add the pileup histograms together - the fits do not need to be done though.
// Be sure to adjust the histogram paths "FitPass#" to the correct iteration (if more than one iterations is done with only one having the energy threshold of 3.5 GeV).

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
#include <TPaveStats.h>
#include <THStack.h>

#include "ratioAnalysisDefs.hh"
#include "pileupUtils.hh"
#include "plotUtils.hh"

using namespace std;


int pileupTimeCorrectedPlots(std::string filePath){

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  gStyle->SetOptStat(000000);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(101);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1);
  gStyle->SetLineColor(1);

  gStyle->SetPadRightMargin(.05);

  // check if dataset tag exists in file and if so append it to the file name and write it to the file

  string outputFileName = "pileupTimeCorrectedPlots";
  TNamed* tag = applyDatasetTag(inputFile, outputFileName);

  TFile* outputFile = new TFile((outputFileName + ".root").c_str(),"RECREATE");
  if(tag) tag->Write();

/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

  TH1F* times_Threshold = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/Input/Times_E_Threshold");
  TH1F* times_Threshold_Psubtracted = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/Input/Times_E_Threshold_Pileup_Subtracted");
  TH1F* times_Threshold_pileup = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/Pileup/added_pileupTimes_shadow_threshold");

  nsTOus(times_Threshold, "");
  nsTOus(times_Threshold_Psubtracted, "");
  nsTOus(times_Threshold_pileup, "");
  times_Threshold_Psubtracted->SetLineColor(2);

/////////////////////////////////////////////////////////////////////////////////////

  auto timeCompCanv = new TCanvas("timeCompCanv", "timeCompCanv", 200, 10, 1200, 800);

  auto legend = new TLegend(0.6,0.6,.9,0.8);
  legend->AddEntry(times_Threshold, "Cluster times", "l");
  legend->AddEntry(times_Threshold_Psubtracted, "Pileup corrected times", "l");

  THStack* histStack = new THStack("histStack","Temp Title");
  histStack->Add(times_Threshold);
  histStack->Add(times_Threshold_Psubtracted);

  histStack->Draw("nostack,hist");

  histStack->GetXaxis()->SetTitle("Time [#mus]");
  histStack->GetYaxis()->SetTitle("Counts above 3.5 GeV");

  legend->SetBorderSize(0);
  legend->Draw();

  histStack->GetXaxis()->SetRangeUser(30, 200);

  timeCompCanv->Write("PileupSubtractedTimesComparison");
  timeCompCanv->SaveAs(("Images/PileupSubtractedTimesComparison" + datasetTagForPlots + ".png").c_str());

/////////////////////////////////////////////////////////////////////////////////////

  auto timeRatioCanv = new TCanvas("timeRatioCanv", "timeRatioCanv", 200, 10, 1200, 800);

  TH1F* dividedHist = (TH1F*) (times_Threshold_Psubtracted->Clone("pileupCorrectedTimes"));
  dividedHist->Divide(times_Threshold);

  dividedHist->GetXaxis()->SetTitle("Time [#mus]");
  dividedHist->GetYaxis()->SetTitle("Pileup corrected times / cluster times (> 3.5 GeV)");

  dividedHist->SetLineColor(1);
  dividedHist->GetXaxis()->SetRangeUser(30, 200);


  TF1* scaleFactor = new TF1("scaleFactor", "[0]", 30, 200);
  scaleFactor->SetLineColor(2);
  dividedHist->Fit(scaleFactor, "RQ");

  dividedHist->Draw();

  timeRatioCanv->Write("PileupSubtractedTimesRatio");
  timeRatioCanv->SaveAs(("Images/PileupSubtractedTimesRatio" + datasetTagForPlots + ".png").c_str());

/////////////////////////////////////////////////////////////////////////////////////

  auto fftRatioCanv = new TCanvas("fftRatioCanv", "fftRatioCanv", 200, 10, 1200, 800);

  cout << "Can replace this bit with my doFFT method in plot utils now." << endl;

  TH1 *fftHistReal = 0;
  TVirtualFFT::SetTransform(0);
  fftHistReal = dividedHist->FFT(fftHistReal, "MAG");

  double histWidth = dividedHist->GetXaxis()->GetXmax()-dividedHist->GetXaxis()->GetXmin();

  TH1F* rescaledFFT = new TH1F((string(fftHistReal->GetName())+"_rescaled").c_str(), fftHistReal->GetTitle(), fftHistReal->GetNbinsX()/2, 0, fftHistReal->GetNbinsX()/(2*histWidth)); // convert to s and then to MHz
  rescaledFFT->GetXaxis()->SetTitle(fftHistReal->GetXaxis()->GetTitle());
  rescaledFFT->GetYaxis()->SetTitle(fftHistReal->GetYaxis()->GetTitle());

  rescaledFFT->GetXaxis()->SetTitle("Freq [MHz]"); // hasn't been converted to MHz yet but will be when rescaled
  rescaledFFT->GetYaxis()->SetTitle("FFT Mag [arb.]");

  for (int bin = 0; bin <= fftHistReal->GetNbinsX()/2; bin++) {
    rescaledFFT->SetBinContent(bin, fftHistReal->GetBinContent(bin));
    rescaledFFT->SetBinError(bin, fftHistReal->GetBinError(bin));
   }

  rescaledFFT->Draw("hist");

  fftRatioCanv->Write("PileupSubtractedTimesRatioFFT");
  fftRatioCanv->SaveAs(("Images/PileupSubtractedTimesRatioFFT" + datasetTagForPlots + ".png").c_str());

/////////////////////////////////////////////////////////////////////////////////////

/*
  auto pileupContamination_canv = new TCanvas("pileupContamination_canv", "pileupContamination_canv", 200, 10, 1200, 800);

  TH1F* dividedHist_contamination = (TH1F*) (times_Threshold_pileup->Clone("pileup"));
  dividedHist_contamination->Divide(times_Threshold);

  dividedHist_contamination->GetXaxis()->SetTitle("Time [#mus]");
  dividedHist_contamination->GetYaxis()->SetTitle("Pileup times / cluster times (> 1.7 GeV)");

  dividedHist_contamination->SetLineColor(1);
  dividedHist_contamination->GetXaxis()->SetRangeUser(30, 200);

  dividedHist_contamination->Draw("hist");

  pileupContamination_canv->Write("PileupTimesComparison");
  pileupContamination_canv->SaveAs(("Images/PileupTimesComparison" + datasetTagForPlots + ".png").c_str());
*/
/////////////////////////////////////////////////////////////////////////////////////



    return 1;

}
