// Macro to plot the time spectra of the pileup doublets over the measured times, in order to get an idea of the rate of doublets to singlets.


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


int pileupRateFractionPlot(std::string filePath){

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

  // string outputFileName = "pileupRatePlot";
  // TNamed* tag = applyDatasetTag(inputFile, outputFileName);

  // TFile* outputFile = new TFile((outputFileName + ".root").c_str(),"RECREATE");
  // if(tag) tag->Write();

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

  auto pileupContamination_canv = new TCanvas("pileupContamination_canv", "pileupContamination_canv", 200, 10, 1200, 800);

  TH1F* dividedHist_contamination = (TH1F*) (times_Threshold_pileup->Clone("pileup"));
  // dividedHist_contamination->Divide(times_Threshold);
  dividedHist_contamination->Divide(times_Threshold_Psubtracted);

  dividedHist_contamination->GetXaxis()->SetTitle("Time [#mus]");
  // dividedHist_contamination->GetYaxis()->SetTitle("Pileup times / cluster times (> 1.7 GeV)");
  dividedHist_contamination->GetYaxis()->SetTitle("Pileup times / corrected times (> 1.7 GeV)");

  dividedHist_contamination->SetLineColor(1);
  dividedHist_contamination->GetXaxis()->SetRangeUser(30, 200);

  dividedHist_contamination->Draw("hist");

  // pileupContamination_canv->Write("pileupRateFraction");
  pileupContamination_canv->SaveAs(("pileupRateFraction" + datasetTagForPlots + ".png").c_str());

/////////////////////////////////////////////////////////////////////////////////////

  // average the rate over the first X microseconds

  double averagingPeriod = g2Period/1000; // us
  int startingBin = dividedHist_contamination->FindBin(30.2876); // default fit starting value
  int endingBin = dividedHist_contamination->FindBin(30.2876 + averagingPeriod); // default fit starting value

  cout << "Starting bin: " << startingBin << " ending bin: " << endingBin << endl;

  double sumOfBinValues = 0;
  for (int i = startingBin; i < endingBin; ++i) sumOfBinValues += dividedHist_contamination->GetBinContent(i);
  double averageBinValue = sumOfBinValues/(endingBin-startingBin);

  cout << "Average bin value of: " << averageBinValue << " over the first " << averagingPeriod << " microseconds." << endl;
  cout << "Maximum value is: " << dividedHist_contamination->GetBinContent(dividedHist_contamination->GetMaximumBin()) << endl;

/////////////////////////////////////////////////////////////////////////////////////

    return 0;
}
