// 3-31-20: Macro to plot R value histograms for different random seeds for different datasets against each other.

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

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

bool saveImages = true;

int RandSeedFitsComparisonPlots()
{

  cout << "With the way I've written the projectGraphToHists method and constructed these histograms, the binning will be different between different datasets and between different sets of fits." << endl;

  gStyle->SetOptStat(000000);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1);
  gStyle->SetLineColor(1);
  gStyle->SetPadRightMargin(.05);
  gStyle->SetPadLeftMargin(.15);


  vector<string> dataset_files = {"/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/60h/RandSeeds/FitIterations/generalPlotsVsIter.root",
                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/HighKick/RandSeeds/FitIterations/generalPlotsVsIter.root",
                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/9d/RandSeeds/FitIterations/generalPlotsVsIter.root",
                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/Endgame/RandSeeds/FitIterations/generalPlotsVsIter.root"};


  auto combined_canv = new TCanvas("combined_canv", "combined_canv", 800, 10, 600, 400);
  // auto combined_canv = new TCanvas("combined_canv", "combined_canv", 10, 10, 1200, 800);

    auto legend = new TLegend(0.8,0.8,.975,0.95);
    legend->SetBorderSize(1);
    legend->SetFillStyle(1001);


  for (uint fileNum = 0; fileNum < dataset_files.size(); ++fileNum)
  {
    TFile* inputFile = TFile::Open(dataset_files.at(fileNum).c_str());
    if(inputFile == 0){
      printf("Error: cannot open file\n");
      return 0;
    }

    TH1F* R_hist = (TH1F*) ((TH1F*) inputFile->Get("topDir/Added/FullRatio/Hists/FullRatio_R_Vs_Iter_hist")->Clone());
    R_hist->SetDirectory(0); // decouple histogram from file/directory otherwise it disappears when closing the file

    R_hist->SetLineColor(dataset_colors.at(fileNum));


    if(fileNum == 0){
      R_hist->GetXaxis()->SetLimits(-21, -16);
      R_hist->GetYaxis()->SetRangeUser(0, 12);
      R_hist->Draw("HIST");
    } 
    else R_hist->Draw("SAME");

    legend->AddEntry(R_hist, dataset_names.at(fileNum).c_str(), "l");

    inputFile->Close();
  }

  legend->Draw("SAME");

  if(saveImages) combined_canv->SaveAs("randSeedFits_R_dataset_comparison.png");

	return 1;

}
