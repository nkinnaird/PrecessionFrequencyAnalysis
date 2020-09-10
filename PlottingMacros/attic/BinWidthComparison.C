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

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////

int BinWidthComparison()
{

  TFile *inputFile1 = TFile::Open("/gm2/data/users/nkinnaird/Ratio/60H-Gold-Report/BinWidth/output-BinWidth-148p9.root");
  TFile *inputFile2 = TFile::Open("/gm2/data/users/nkinnaird/Ratio/60H-Gold-Report/BinWidth/output-BinWidth-149.root");
  TFile *inputFile3 = TFile::Open("/gm2/data/users/nkinnaird/Ratio/60H-Gold-Report/BinWidth/output-BinWidth-149p1.root");
  TFile *inputFile4 = TFile::Open("/gm2/data/users/nkinnaird/Ratio/60H-Gold-Report/SingleIteration/SingleFit/output-Gold-SingleIteration.root");
  TFile *inputFile5 = TFile::Open("/gm2/data/users/nkinnaird/Ratio/60H-Gold-Report/BinWidth/output-BinWidth-149p3.root");
  TFile *inputFile6 = TFile::Open("/gm2/data/users/nkinnaird/Ratio/60H-Gold-Report/BinWidth/output-BinWidth-149p4.root");
  TFile *inputFile7 = TFile::Open("/gm2/data/users/nkinnaird/Ratio/60H-Gold-Report/BinWidth/output-BinWidth-149p5.root");
  TFile *inputFile8 = TFile::Open("/gm2/data/users/nkinnaird/Ratio/60H-Gold-Report/BinWidth/output-BinWidth-149p6.root");


   if (inputFile1 == 0 || inputFile2 == 0 || inputFile3 == 0 || inputFile4 == 0 || inputFile5 == 0 || inputFile6 == 0 || inputFile7 == 0 || inputFile8 == 0){
      printf("Error: cannot open file\n");
      return 0;
   }


/////////////////////////////////////////////////////////////////////////////////////
  // These only get set for the interactive root session (any generated canvases, etc.), but does not apply to the output root file - that comes from .rootlogon.C
  gStyle->SetOptStat(000000);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(1111);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1);
  gStyle->SetLineColor(1);
  gStyle->SetPadRightMargin(.05);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

std::vector<double> binWidthValues = {148.9, 149.0, 149.1, 149.2, 149.30, 149.4, 149.5, 149.6};
vector<pair<double, double> > data;

      TF1* fullRatioFitFunction1 = (TF1*) ((TGraphErrors*) inputFile1->Get("topDir/FitPasses/FitPass0/addedDir/FullRatio/Added_Times_Full_Ratio_Graph"))->GetFunction("fullRatioFitFunc")->Clone();
      TF1* fullRatioFitFunction2 = (TF1*) ((TGraphErrors*) inputFile2->Get("topDir/FitPasses/FitPass0/addedDir/FullRatio/Added_Times_Full_Ratio_Graph"))->GetFunction("fullRatioFitFunc")->Clone();
      TF1* fullRatioFitFunction3 = (TF1*) ((TGraphErrors*) inputFile3->Get("topDir/FitPasses/FitPass0/addedDir/FullRatio/Added_Times_Full_Ratio_Graph"))->GetFunction("fullRatioFitFunc")->Clone();
      TF1* fullRatioFitFunction4 = (TF1*) ((TGraphErrors*) inputFile4->Get("topDir/FitPasses/FitPass0/addedDir/FullRatio/Added_Times_Full_Ratio_Graph"))->GetFunction("fullRatioFitFunc")->Clone();
      TF1* fullRatioFitFunction5 = (TF1*) ((TGraphErrors*) inputFile5->Get("topDir/FitPasses/FitPass0/addedDir/FullRatio/Added_Times_Full_Ratio_Graph"))->GetFunction("fullRatioFitFunc")->Clone();
      TF1* fullRatioFitFunction6 = (TF1*) ((TGraphErrors*) inputFile6->Get("topDir/FitPasses/FitPass0/addedDir/FullRatio/Added_Times_Full_Ratio_Graph"))->GetFunction("fullRatioFitFunc")->Clone();
      TF1* fullRatioFitFunction7 = (TF1*) ((TGraphErrors*) inputFile7->Get("topDir/FitPasses/FitPass0/addedDir/FullRatio/Added_Times_Full_Ratio_Graph"))->GetFunction("fullRatioFitFunc")->Clone();
      TF1* fullRatioFitFunction8 = (TF1*) ((TGraphErrors*) inputFile8->Get("topDir/FitPasses/FitPass0/addedDir/FullRatio/Added_Times_Full_Ratio_Graph"))->GetFunction("fullRatioFitFunc")->Clone();


      data.emplace_back(fullRatioFitFunction1->GetParameter(1), fullRatioFitFunction1->GetParError(1));
      data.emplace_back(fullRatioFitFunction2->GetParameter(1), fullRatioFitFunction2->GetParError(1));
      data.emplace_back(fullRatioFitFunction3->GetParameter(1), fullRatioFitFunction3->GetParError(1));
      data.emplace_back(fullRatioFitFunction4->GetParameter(1), fullRatioFitFunction4->GetParError(1));
      data.emplace_back(fullRatioFitFunction5->GetParameter(1), fullRatioFitFunction5->GetParError(1));
      data.emplace_back(fullRatioFitFunction6->GetParameter(1), fullRatioFitFunction6->GetParError(1));
      data.emplace_back(fullRatioFitFunction7->GetParameter(1), fullRatioFitFunction7->GetParError(1));
      data.emplace_back(fullRatioFitFunction8->GetParameter(1), fullRatioFitFunction8->GetParError(1));

/////////////////////////////////////////////////////////////////////////////////////



    TGraphErrors* comparisonGraph = new TGraphErrors();

    for (uint pointNo = 0; pointNo < data.size(); ++pointNo)
    {
      comparisonGraph->SetPoint(pointNo, binWidthValues.at(pointNo), data.at(pointNo).first);
      comparisonGraph->SetPointError(pointNo, 0, data.at(pointNo).second);
    }


      comparisonGraph->SetName("BinWidthComparison");
      comparisonGraph->SetTitle("Full Ratio Fit R Vs Bin Width");
      comparisonGraph->GetXaxis()->SetTitle("Bin Width (ns)");
      comparisonGraph->GetYaxis()->SetTitle("R (ppm)");
      comparisonGraph->GetYaxis()->SetTitleOffset(2.0);

      // comparisonGraph->GetYaxis()->SetRangeUser(-.15, .15);

      auto binWidth_canv = new TCanvas("binWidth_canv","binWidth_canv",200,10,500,400);
      comparisonGraph->Draw("AP");

      double mean = comparisonGraph->GetMean(2);
      double rms = comparisonGraph->GetRMS(2);

      TPaveText* textBox = new TPaveText(0.75,0.7,0.975,0.9, "NDC");
      textBox->AddText(Form("Mean  %0.4g", mean));
      textBox->AddText(Form("RMS  %0.4g", rms));
      textBox->SetBorderSize(1);
      textBox->SetFillStyle(1001);
      textBox->Draw("SAME");

      binWidth_canv->SaveAs("Images/BinWidthComparison_R.png");


    TH1F* comparisonHist = projectGraphToHist(comparisonGraph);

      auto histCanvas = new TCanvas("binWidth_canv_hist", "binWidth_canv_hist",200,10,500,400);
        comparisonHist->Draw();
        textBox->Draw("SAME");

        histCanvas->SaveAs("Images/BinWidthComparison_R_hist.png");


/////////////////////////////////////////////////////////////////////////////////////

  return 1;
}
