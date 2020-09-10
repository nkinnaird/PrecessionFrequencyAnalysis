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

double calculateAverage(vector<double> inVec)
{
  double avg = 0;
  for (uint i = 0; i < inVec.size(); ++i) avg += inVec.at(i);
  avg /= inVec.size();
  return avg;
}


int DiffRandComparison(std::string fileList)
{

  gStyle->SetOptStat(2211);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleStyle(0);
  gStyle->SetPadRightMargin(.05);
  gStyle->SetPadLeftMargin(.15);

/////////////////////////////////////////////////////////////////////////////////////

  // first vector is the level of randomization, second vector are the R related parameters for that randomization

  vector<vector<double>> fiveParMeans(6);
  vector<vector<double>> fiveParWidths(6);
  vector<vector<double>> ratioMeans(6);
  vector<vector<double>> ratioWidths(6);

  vector<vector<double>> TTmeanDiffs(6); // Tmean - Tmean with no randomization
  vector<vector<double>> RRmeanDiffs(6); // Rmean - Rmean with no randomization
  vector<vector<double>> RTmeanDiffs(6); // Rmean - Tmean 

  vector<double> randomizations {0, 150, 300, 450, 600, 750};

/////////////////////////////////////////////////////////////////////////////////////

  std::vector<string> fileVector;

  std::ifstream inList(fileList);
  string path;
  while(inList >> path) fileVector.push_back(path);

  for (uint fileNum = 0; fileNum < fileVector.size(); ++fileNum)
  {
    TFile* inputFile = TFile::Open(fileVector.at(fileNum).c_str());
     if (inputFile == 0) {
        printf("Error: cannot open file\n");
        return 0;
     }

    cout << "Filenum: " << fileNum << " " << fileVector.at(fileNum) << endl;

/////////////////////////////////////////////////////////////////////////////////////

    for (uint i = 0; i < randomizations.size(); ++i)
    {
      TH1F* toyFiveParamSingleJobRs;
      TH1F* toyRatioSingleJobRs;

      if(i == 0){
        toyFiveParamSingleJobRs = (TH1F*) inputFile->Get("topDir/ToyMC/Hist/toyHistFiveR");
        toyRatioSingleJobRs = (TH1F*) inputFile->Get("topDir/ToyMC/Ratio/toyRatioR");
      }
      else{
        toyFiveParamSingleJobRs = (TH1F*) inputFile->Get(Form("topDir/ToyMC_%i/Hist/toyHistFiveR_%i", i, i));
        toyRatioSingleJobRs = (TH1F*) inputFile->Get(Form("topDir/ToyMC_%i/Ratio/toyRatioR_%i", i, i));
      }

      double fiveParMean = toyFiveParamSingleJobRs->GetMean();
      double fiveParWidth = toyFiveParamSingleJobRs->GetRMS();
      double ratioMean = toyRatioSingleJobRs->GetMean();
      double ratioWidth = toyRatioSingleJobRs->GetRMS();

      fiveParMeans.at(i).push_back(fiveParMean);
      fiveParWidths.at(i).push_back(fiveParWidth);
      ratioMeans.at(i).push_back(ratioMean);
      ratioWidths.at(i).push_back(ratioWidth);

      RTmeanDiffs.at(i).push_back(ratioMean - fiveParMean);
    }

/////////////////////////////////////////////////////////////////////////////////////

    inputFile->Close();
  }

/////////////////////////////////////////////////////////////////////////////////////

  // for (uint i = 0; i < fiveParWidths.size(); ++i)
  // {
  //   cout << "Rand seed level: " << i << " size: " << fiveParWidths.at(i).size() << endl;
  //   for (uint j = 0; j < fiveParWidths.at(i).size(); ++j)
  //   {
  //     cout << fiveParWidths.at(i).at(j) << endl;
  //   }
  //   cout << endl;
  // }

/////////////////////////////////////////////////////////////////////////////////////

  vector<double> fiveParMeanAvgs;
  vector<double> fiveParWidthAvgs;
  vector<double> ratioMeanAvgs;
  vector<double> ratioWidthAvgs;

  for (uint i = 0; i < fiveParMeans.size(); ++i) fiveParMeanAvgs.push_back(calculateAverage(fiveParMeans.at(i)));
  for (uint i = 0; i < fiveParWidths.size(); ++i) fiveParWidthAvgs.push_back(calculateAverage(fiveParWidths.at(i)));
  for (uint i = 0; i < ratioMeans.size(); ++i) ratioMeanAvgs.push_back(calculateAverage(ratioMeans.at(i)));
  for (uint i = 0; i < ratioWidths.size(); ++i) ratioWidthAvgs.push_back(calculateAverage(ratioWidths.at(i)));

/////////////////////////////////////////////////////////////////////////////////////


    TGraph* fiveParamMeanGraph = new TGraph();
    TGraph* fiveParamWidthGraph = new TGraph();
    TGraph* ratioMeanGraph = new TGraph();
    TGraph* ratioWidthGraph = new TGraph();

    for (uint i = 0; i < randomizations.size(); ++i)
    {
      fiveParamMeanGraph->SetPoint(i, randomizations.at(i), fiveParMeanAvgs.at(i));
      fiveParamWidthGraph->SetPoint(i, randomizations.at(i), fiveParWidthAvgs.at(i));
      ratioMeanGraph->SetPoint(i, randomizations.at(i), ratioMeanAvgs.at(i));
      ratioWidthGraph->SetPoint(i, randomizations.at(i), ratioWidthAvgs.at(i));
    }


    // fiveParamMeanGraph->GetXaxis()->SetTitle("Randomization in ns");
    // fiveParamMeanGraph->GetYaxis()->SetTitle("Five Parameter Average Mean (ppm)");
    // auto fiveParMeanCanv = new TCanvas("fiveParMean", "fiveParMean", 20, 20, 500, 400);
    // fiveParamMeanGraph->Draw("AP");

    fiveParamWidthGraph->GetXaxis()->SetTitle("Randomization in ns");
    fiveParamWidthGraph->GetYaxis()->SetTitle("Five Parameter Average Width (ppm)");
    fiveParamWidthGraph->GetYaxis()->SetTitleOffset(2);
    auto fiveParWidthCanv = new TCanvas("fiveParWidth", "fiveParWidth", 40, 20, 500, 400);
    fiveParamWidthGraph->Draw("AP");
    fiveParWidthCanv->SaveAs("Images/fiveParWidth.png");


    // ratioMeanGraph->GetXaxis()->SetTitle("Randomization in ns");
    // ratioMeanGraph->GetYaxis()->SetTitle("Ratio Average Mean (ppm)");
    // auto ratioMeanCanv = new TCanvas("ratioMean", "ratioMean", 60, 20, 500, 400);
    // ratioMeanGraph->Draw("AP");

    ratioWidthGraph->GetXaxis()->SetTitle("Randomization in ns");
    ratioWidthGraph->GetYaxis()->SetTitle("Ratio Average Width (ppm)");
    ratioWidthGraph->GetYaxis()->SetTitleOffset(2);
    auto ratioWidthCanv = new TCanvas("ratioWidth", "ratioWidth", 440, 20, 500, 400);
    ratioWidthGraph->Draw("AP");
    ratioWidthCanv->SaveAs("Images/ratioWidth.png");

/////////////////////////////////////////////////////////////////////////////////////

    vector<double> RTdiffWidths;

    for (uint i = 0; i < randomizations.size(); ++i)
    // for (uint i = 0; i < 1; ++i)
    {
      auto RT_Rmean_Diffs = new TH1F("RT_Rmean_Diffs", "Ratio Fit R Mean - Five Param Fit R Mean; R Mean Diff (ppm); Events", 100, -1, 1);
      for (uint j = 0; j < RTmeanDiffs.at(i).size(); ++j) RT_Rmean_Diffs->Fill(RTmeanDiffs.at(i).at(j));

        RTdiffWidths.push_back(RT_Rmean_Diffs->GetRMS());

      // auto RTmeanDiffCanv = new TCanvas("RTmeanDiffCanv", "RTmeanDiffCanv", 40, 20, 500, 400);
      // RT_Rmean_Diffs->Draw("HIST");

      //   RTmeanDiffCanv->Update();

      //   TPaveStats* stats_diff = (TPaveStats*) RT_Rmean_Diffs->GetListOfFunctions()->FindObject("stats");
      //   stats_diff->SetX1NDC(.60);
      //   stats_diff->SetX2NDC(.97);
      //   stats_diff->SetY1NDC(.63);
      //   stats_diff->SetY2NDC(.9);
      //   stats_diff->SetBorderSize(1);

      // RTmeanDiffCanv->DrawClonePad();

      // delete RTmeanDiffCanv;

    }

    TGraph* RTdiffWidthGraph = new TGraph();
    for (uint i = 0; i < randomizations.size(); ++i) RTdiffWidthGraph->SetPoint(i, randomizations.at(i), RTdiffWidths.at(i));

    RTdiffWidthGraph->GetXaxis()->SetTitle("Randomization in ns");
    RTdiffWidthGraph->GetYaxis()->SetTitle("Ratio - Five Parameter Mean Diff Width (ppm)");
    RTdiffWidthGraph->GetYaxis()->SetTitleOffset(2);
    auto RTdiffWidthCanv = new TCanvas("RTdiffWidth", "RTdiffWidth", 40, 200, 600, 480);
    RTdiffWidthGraph->Draw("AP");
    RTdiffWidthCanv->SaveAs("Images/RTdiffWidth.png");



/////////////////////////////////////////////////////////////////////////////////////

return 1;

}
