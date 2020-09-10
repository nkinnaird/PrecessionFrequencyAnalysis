
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
#include <TPaveStats.h>

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////


void makeDiffPlot(TH1F* hist, string histName, string title, string canvName)
{
    TH1F* histClone = (TH1F*) hist->Clone(histName.c_str());

    auto diff_canvas = new TCanvas(canvName.c_str(), "title", 200, 10, 1000, 800);

      histClone->SetTitle(title.c_str());
      histClone->Draw();

        diff_canvas->Update();

        TPaveStats* stats_diff = (TPaveStats*) histClone->GetListOfFunctions()->FindObject("stats");
        stats_diff->SetX1NDC(.70);
        stats_diff->SetX2NDC(.97);
        stats_diff->SetY1NDC(.63);
        stats_diff->SetY2NDC(.9);
        stats_diff->SetBorderSize(1);

        diff_canvas->Modified();

      diff_canvas->Write();
      diff_canvas->SaveAs(("Images/" + canvName + ".png").c_str());

      diff_canvas->DrawClonePad();

      delete histClone;
      delete diff_canvas;
}


void makeComparisonPlot(TH1F* firstHist, TH1F* secondHist, string firstHistName, string secondHistName, string title, string canvName)
{
    // clone histograms in order to change the names for plot stat boxes, without changing names of the histograms saved into the file
    TH1F* firstHistClone = (TH1F*) firstHist->Clone(firstHistName.c_str());
    TH1F* secondHistClone = (TH1F*) secondHist->Clone(secondHistName.c_str());

    auto compCanv = new TCanvas(canvName.c_str(), "title", 200, 10, 1000, 800);

      firstHistClone->SetTitle(title.c_str());
      firstHistClone->Draw();
      
      secondHistClone->SetLineColor(2);
      secondHistClone->Draw("sames");

        compCanv->Update();

        TPaveStats* stats_first = (TPaveStats*) firstHistClone->GetListOfFunctions()->FindObject("stats");
        stats_first->SetX1NDC(.70);
        stats_first->SetX2NDC(.97);
        stats_first->SetY1NDC(.63);
        stats_first->SetY2NDC(.9);
        stats_first->SetBorderSize(1);

        compCanv->Modified();

        TPaveStats* stats_second = (TPaveStats*) secondHistClone->GetListOfFunctions()->FindObject("stats");
        stats_second->SetX1NDC(.70);
        stats_second->SetX2NDC(.97);
        stats_second->SetY1NDC(.36);
        stats_second->SetY2NDC(.63);
        stats_second->SetTextColor(2);
        stats_second->SetLineColor(2);
        stats_second->SetBorderSize(1);

        compCanv->Modified();

      compCanv->Write();
      compCanv->SaveAs(("Images/" + canvName + ".png").c_str());

      compCanv->DrawClonePad(); // draw clone pad so that I can delete the associated histograms such that they do not get saved to the file, and so that the image on the screen is preserved

      delete firstHistClone;
      delete secondHistClone;
      delete compCanv;
}

int TRComparison(std::string fileList)
{
  // gStyle->SetOptStat(0);
  // gStyle->SetOptTitle(1);
  // gStyle->SetOptFit(2);
  // gStyle->SetMarkerStyle(20);
  // gStyle->SetMarkerColor(1);
  // gStyle->SetMarkerSize(1);
  // gStyle->SetLineColor(1);


  gStyle->SetOptStat(1111);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleStyle(0);


  gStyle->SetPadRightMargin(.05);


  // create output file that will hold plots
  TFile* outputFile = new TFile("TRComparison.root","RECREATE");

  // make top directory for output file
  auto topDir = outputFile->mkdir("topDir");
  auto individualRDist_Dir = topDir->mkdir("IndividualRDists");
  auto combinedDist_Dir = topDir->mkdir("CombinedRDists");
  auto canvasPlots_Dir = topDir->mkdir("Plots");

/////////////////////////////////////////////////////////////////////////////////////

  uint singleDistributionNumber = 16;
  TH1F* savedSingle_5Param_RDist = 0;
  TH1F* savedSingle_Ratio_RDist = 0;

/////////////////////////////////////////////////////////////////////////////////////

  combinedDist_Dir->cd();

  auto fiveParam_RMeans = new TH1F("fiveParam_RMeans", "Five Param Fit R Means; R Mean (ppm); Events", 100, -5, 5);
  auto fiveParam_RWidths = new TH1F("fiveParam_RWidths", "Five Param Fit R Widths; R Width (ppm); Events", 50, 0, .5);

  auto Ratio_RMeans = new TH1F("Ratio_RMeans", "Ratio Fit R Means; R Mean (ppm); Events", 100, -5, 5);
  auto Ratio_RWidths = new TH1F("Ratio_RWidths", "Ratio Fit R Widths; R Width (ppm); Events", 50, 0, .5);

  auto Rmean_Diffs = new TH1F("Rmean_Diffs", "Ratio Fit R Mean - Five Param Fit R Mean; R Mean Diff (ppm); Events", 100, -1, 1);

/////////////////////////////////////////////////////////////////////////////////////

  std::vector<string> fileVector;

  std::ifstream inList(fileList);
  string path;
  while(inList >> path) fileVector.push_back(path);

  if(int(fileVector.size() < singleDistributionNumber))
  {
    cout << "Trying to grab a single distribution but I didn't run enough files:" << endl;
    return -1;
  }

for (uint fileNum = 0; fileNum < fileVector.size(); ++fileNum)
{

  TFile* inputFile = TFile::Open(fileVector.at(fileNum).c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  cout << "Filenum: " << fileNum << " " << fileVector.at(fileNum) << endl;

  individualRDist_Dir->cd();

    // TH1F* toyFiveParamSingleJobRs = (TH1F*) ((TH1F*) inputFile->Get("topDir/ToyMC/Hist/Stats/toyHistFiveR"))->Clone(); // has extra Stats directory used in old comparison work
    // TH1F* toyRatioSingleJobRs = (TH1F*) ((TH1F*) inputFile->Get("topDir/ToyMC/Ratio/Stats/toyRatioR"))->Clone();
    TH1F* toyFiveParamSingleJobRs = (TH1F*) ((TH1F*) inputFile->Get("topDir/ToyMC/Hist/toyHistFiveR"))->Clone();
    TH1F* toyRatioSingleJobRs = (TH1F*) ((TH1F*) inputFile->Get("topDir/ToyMC/Ratio/toyRatioR"))->Clone();

    if(fileNum == singleDistributionNumber) // choose some number here such that the distributions aren't right on 0
    {
      savedSingle_5Param_RDist = (TH1F*) toyFiveParamSingleJobRs->Clone("tempClone_5param");
      savedSingle_Ratio_RDist = (TH1F*) toyRatioSingleJobRs->Clone("tempClone_ratio");
    }


    double fiveParMean = toyFiveParamSingleJobRs->GetMean();
    double ratioMean = toyRatioSingleJobRs->GetMean();

    fiveParam_RMeans->Fill(fiveParMean);
    Ratio_RMeans->Fill(ratioMean);

    fiveParam_RWidths->Fill(toyFiveParamSingleJobRs->GetRMS());
    Ratio_RWidths->Fill(toyRatioSingleJobRs->GetRMS());

    Rmean_Diffs->Fill(ratioMean - fiveParMean);

  inputFile->Close();
}

/////////////////////////////////////////////////////////////////////////////////////

  canvasPlots_Dir->cd();

    makeComparisonPlot(savedSingle_5Param_RDist, savedSingle_Ratio_RDist, "5 Param Fit", "Ratio Fit", "R For Many Random Seeds", "SingleDistributions");
    makeComparisonPlot(fiveParam_RMeans, Ratio_RMeans, "5 Param Fit", "Ratio Fit", "Distribution of Fitted R Means", "MeanDistributions");
    makeComparisonPlot(fiveParam_RWidths, Ratio_RWidths, "5 Param Fit", "Ratio Fit", "Distribution of Fitted R Widths", "WidthDistributions");

    makeDiffPlot(Rmean_Diffs, "R Mean Diffs", "Ratio Fit R Mean - Five Param Fit R Mean", "DiffDistribution");

/////////////////////////////////////////////////////////////////////////////////////

  outputFile->Write();
  // delete outputFile;

  return 1;

}
