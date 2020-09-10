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

int globalCanvX = 0;
int globalCanvY = 10;

/////////////////////////////////////////////////////////////////////////////////////

void makeDiffPlot(TH1F* hist, string histName, string title, string canvName)
{
    TH1F* histClone = (TH1F*) hist->Clone(histName.c_str());

    auto diff_canvas = new TCanvas(canvName.c_str(), "title", globalCanvX, globalCanvY, 600, 480);

      histClone->SetTitle(title.c_str());
      histClone->Draw();

        diff_canvas->Update();

        TPaveStats* stats_diff = (TPaveStats*) histClone->GetListOfFunctions()->FindObject("stats");
        stats_diff->SetX1NDC(.60);
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

      globalCanvX += 100;
}


void makeComparisonPlot(TH1F* firstHist, TH1F* secondHist, string firstHistName, string secondHistName, string title, string canvName)
{
    // clone histograms in order to change the names for plot stat boxes, without changing names of the histograms saved into the file
    TH1F* firstHistClone = (TH1F*) firstHist->Clone(firstHistName.c_str());
    TH1F* secondHistClone = (TH1F*) secondHist->Clone(secondHistName.c_str());

    auto compCanv = new TCanvas(canvName.c_str(), "title", globalCanvX, globalCanvY, 600, 480);

      firstHistClone->SetTitle(title.c_str());
      firstHistClone->Draw();
      
      secondHistClone->SetLineColor(2);
      secondHistClone->Draw("sames");

        compCanv->Update();

        TPaveStats* stats_first = (TPaveStats*) firstHistClone->GetListOfFunctions()->FindObject("stats");
        stats_first->SetX1NDC(.60);
        stats_first->SetX2NDC(.97);
        stats_first->SetY1NDC(.63);
        stats_first->SetY2NDC(.9);
        stats_first->SetBorderSize(1);
        stats_first->Draw("same");

        compCanv->Modified();

        TPaveStats* stats_second = (TPaveStats*) secondHistClone->GetListOfFunctions()->FindObject("stats");
        stats_second->SetX1NDC(.60);
        stats_second->SetX2NDC(.97);
        stats_second->SetY1NDC(.357);
        stats_second->SetY2NDC(.627);
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

      globalCanvX += 100;
}


int VWRandComparison(std::string fileList)
{
  // gStyle->SetOptStat(0);
  // gStyle->SetOptTitle(1);
  // gStyle->SetOptFit(2);
  // gStyle->SetMarkerStyle(20);
  // gStyle->SetMarkerColor(1);
  // gStyle->SetMarkerSize(1);
  // gStyle->SetLineColor(1);


  // gStyle->SetOptStat(1111);
  gStyle->SetOptStat(2211);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleStyle(0);


  gStyle->SetPadRightMargin(.05);


  // create output file that will hold plots
  TFile* outputFile = new TFile("VWRandComparison.root","RECREATE");

  // make top directory for output file
  auto NoVWRand_Dir = outputFile->mkdir("NoVWRand");
  auto individualRDist_Dir = NoVWRand_Dir->mkdir("IndividualRDists");
  auto combinedDist_Dir = NoVWRand_Dir->mkdir("CombinedRDists");

  auto VWRand_Dir = outputFile->mkdir("VWRand");
  auto individualRDist_VW_Dir = VWRand_Dir->mkdir("IndividualRDists");
  auto combinedDist_VW_Dir = VWRand_Dir->mkdir("CombinedRDists");

/////////////////////////////////////////////////////////////////////////////////////

  combinedDist_Dir->cd();

  auto fiveParam_RMeans = new TH1F("fiveParam_RMeans", "Five Param Fit R Means; R Mean (ppm); Events", 100, -5, 5);
  auto fiveParam_RWidths = new TH1F("fiveParam_RWidths", "Five Param Fit R Widths; R Width (ppm); Events", 50, 0, .5);

  auto Ratio_RMeans = new TH1F("Ratio_RMeans", "Ratio Fit R Means; R Mean (ppm); Events", 100, -5, 5);
  auto Ratio_RWidths = new TH1F("Ratio_RWidths", "Ratio Fit R Widths; R Width (ppm); Events", 50, 0, .5);

  auto RT_Rmean_Diffs = new TH1F("RT_Rmean_Diffs", "Ratio Fit R Mean - Five Param Fit R Mean; R Mean Diff (ppm); Events", 100, -1, 1);

/////////////////////////////////////////////////////////////////////////////////////

  combinedDist_VW_Dir->cd();

  auto fiveParam_RMeans_VW = new TH1F("fiveParam_RMeans_VW", "Five Param Fit R Means; R Mean (ppm); Events", 100, -5, 5);
  auto fiveParam_RWidths_VW = new TH1F("fiveParam_RWidths_VW", "Five Param Fit R Widths; R Width (ppm); Events", 50, 0, .5);

  auto Ratio_RMeans_VW = new TH1F("Ratio_RMeans_VW", "Ratio Fit R Means; R Mean (ppm); Events", 100, -5, 5);
  auto Ratio_RWidths_VW = new TH1F("Ratio_RWidths_VW", "Ratio Fit R Widths; R Width (ppm); Events", 50, 0, .5);

  auto RT_Rmean_Diffs_VW = new TH1F("RT_Rmean_Diffs_VW", "Ratio Fit R Mean - Five Param Fit R Mean; R Mean Diff (ppm); Events", 100, -1, 1);

/////////////////////////////////////////////////////////////////////////////////////

  auto VWRand_vw_NoVWRand_Dir = outputFile->mkdir("VWRand_Vs_NoVWRand");
  VWRand_vw_NoVWRand_Dir->cd();

  auto TT_Rmean_Diffs = new TH1F("TT_Rmean_Diffs", "Five Param R Mean (VW rand) - Five Param Fit R Mean; R Mean Diff (ppm); Events", 100, -1, 1);
  auto RR_Rmean_Diffs = new TH1F("RR_Rmean_Diffs", "Five Param R Mean (VW rand) - Five Param Fit R Mean; R Mean Diff (ppm); Events", 100, -1, 1);


/////////////////////////////////////////////////////////////////////////////////////

  uint singleDistributionNumber = 16;
  TH1F* savedSingle_5Param_RDist = 0;
  TH1F* savedSingle_Ratio_RDist = 0;
  TH1F* savedSingle_5Param_RDist_VW = 0;
  TH1F* savedSingle_Ratio_RDist_VW = 0;

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

/////////////////////////////////////////////////////////////////////////////////////

  individualRDist_Dir->cd();

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

    RT_Rmean_Diffs->Fill(ratioMean - fiveParMean);

/////////////////////////////////////////////////////////////////////////////////////

  individualRDist_VW_Dir->cd();

    TH1F* toyFiveParamSingleJobRs_VW = (TH1F*) ((TH1F*) inputFile->Get("topDir/ToyMC_VW/Hist/toyHistFiveR_VW"))->Clone();
    TH1F* toyRatioSingleJobRs_VW = (TH1F*) ((TH1F*) inputFile->Get("topDir/ToyMC_VW/Ratio/toyRatioR_VW"))->Clone();

    if(fileNum == singleDistributionNumber) // choose some number here such that the distributions aren't right on 0
    {
      savedSingle_5Param_RDist_VW = (TH1F*) toyFiveParamSingleJobRs_VW->Clone("tempClone_5param");
      savedSingle_Ratio_RDist_VW = (TH1F*) toyRatioSingleJobRs_VW->Clone("tempClone_ratio");
    }

    double fiveParMean_VW = toyFiveParamSingleJobRs_VW->GetMean();
    double ratioMean_VW = toyRatioSingleJobRs_VW->GetMean();

    fiveParam_RMeans_VW->Fill(fiveParMean_VW);
    Ratio_RMeans_VW->Fill(ratioMean_VW);

    fiveParam_RWidths_VW->Fill(toyFiveParamSingleJobRs_VW->GetRMS());
    Ratio_RWidths_VW->Fill(toyRatioSingleJobRs_VW->GetRMS());

    RT_Rmean_Diffs_VW->Fill(ratioMean_VW - fiveParMean_VW);

/////////////////////////////////////////////////////////////////////////////////////

    TT_Rmean_Diffs->Fill(fiveParMean_VW - fiveParMean);
    RR_Rmean_Diffs->Fill(ratioMean_VW - ratioMean);

/////////////////////////////////////////////////////////////////////////////////////

  inputFile->Close();
}

/////////////////////////////////////////////////////////////////////////////////////

  // R - T differences without VW randomization

  auto T_Vs_R_NoVWRand_Plots_Dir = outputFile->mkdir("T_Vs_R_NoVWRand_Plots");
  T_Vs_R_NoVWRand_Plots_Dir->cd();

    makeComparisonPlot(savedSingle_5Param_RDist, savedSingle_Ratio_RDist, "5 Param Fit", "Ratio Fit", "R For Many Random Seeds (no VW rand)", "SingleDistributions_RT_noVWRand");
    makeComparisonPlot(fiveParam_RMeans, Ratio_RMeans, "5 Param Fit", "Ratio Fit", "Distribution of Fitted R Means (no VW rand)", "MeanDistributions_RT_noVWRand");
    makeComparisonPlot(fiveParam_RWidths, Ratio_RWidths, "5 Param Fit", "Ratio Fit", "Distribution of Fitted R Widths (no VW rand)", "WidthDistributions_RT_noVWRand");

    makeDiffPlot(RT_Rmean_Diffs, "R Mean Diffs", "Ratio Fit R Mean - Five Param Fit R Mean (no VW rand)", "DiffDistribution_RT_noVWRand");

/////////////////////////////////////////////////////////////////////////////////////

    // R - T differences with VW randomization

    globalCanvX = 0;
    globalCanvY += 500;

  auto T_Vs_R_VWRand_Plots_Dir = outputFile->mkdir("T_Vs_R_VWRand_Plots");
  T_Vs_R_VWRand_Plots_Dir->cd();

    makeComparisonPlot(savedSingle_5Param_RDist_VW, savedSingle_Ratio_RDist_VW, "5 Param Fit", "Ratio Fit", "R For Many Random Seeds (VW rand)", "SingleDistributions_RT_VWRand");
    makeComparisonPlot(fiveParam_RMeans_VW, Ratio_RMeans_VW, "5 Param Fit", "Ratio Fit", "Distribution of Fitted R Means (VW rand)", "MeanDistributions_RT_VWRand");
    makeComparisonPlot(fiveParam_RWidths_VW, Ratio_RWidths_VW, "5 Param Fit", "Ratio Fit", "Distribution of Fitted R Widths (VW rand)", "WidthDistributions_RT_VWRand");

    makeDiffPlot(RT_Rmean_Diffs_VW, "R Mean Diffs", "Ratio Fit R Mean - 5 Param Fit R Mean (VW rand)", "DiffDistribution_RT_VWRand");

/////////////////////////////////////////////////////////////////////////////////////

    // T - T differences with and without VW randomization

    globalCanvX = 900;
    globalCanvY = 10;

  auto T_Vs_T_WWoVWRand_Plots_Dir = outputFile->mkdir("T_Vs_T_WWoVWRand_Plots");
  T_Vs_T_WWoVWRand_Plots_Dir->cd();

    makeComparisonPlot(savedSingle_5Param_RDist, savedSingle_5Param_RDist_VW, "5 Param Fit", "5 Param Fit (VW rand)", "R For Many Random Seeds", "SingleDistributions_TT_WWoVWRand");
    makeComparisonPlot(fiveParam_RMeans, fiveParam_RMeans_VW, "5 Param Fit", "5 Param Fit (VW rand)", "Distribution of Fitted R Means", "MeanDistributions_TT_WWoVWRandd");
    makeComparisonPlot(fiveParam_RWidths, fiveParam_RWidths_VW, "5 Param Fit", "5 Param Fit (VW rand)", "Distribution of Fitted R Widths", "WidthDistributions_TT_WWoVWRand");

    makeDiffPlot(TT_Rmean_Diffs, "R Mean Diffs", "5 Param Fit R Mean (VW rand) - 5 Param Fit R Mean", "DiffDistribution_TT_WWoVWRand");

/////////////////////////////////////////////////////////////////////////////////////

    // R - R differences with and without VW randomization

    globalCanvX = 900;
    globalCanvY += 500;

  auto R_Vs_R_WWoVWRand_Plots_Dir = outputFile->mkdir("R_Vs_R_WWoVWRand_Plots");
  R_Vs_R_WWoVWRand_Plots_Dir->cd();

    makeComparisonPlot(savedSingle_Ratio_RDist, savedSingle_Ratio_RDist_VW, "Ratio Fit", "Ratio Fit (VW rand)", "R For Many Random Seeds", "SingleDistributions_RR_WWoVWRand");
    makeComparisonPlot(Ratio_RMeans, Ratio_RMeans_VW, "Ratio Fit", "Ratio Fit (VW rand)", "Distribution of Fitted R Means", "MeanDistributions_RR_WWoVWRandd");
    makeComparisonPlot(Ratio_RWidths, Ratio_RWidths_VW, "Ratio Fit", "Ratio Fit (VW rand)", "Distribution of Fitted R Widths", "WidthDistributions_RR_WWoVWRand");

    makeDiffPlot(RR_Rmean_Diffs, "R Mean Diffs", "Ratio Fit R Mean (VW rand) - Ratio Fit R Mean", "DiffDistribution_RR_WWoVWRand");


/////////////////////////////////////////////////////////////////////////////////////

  outputFile->Write();

  return 1;

}
