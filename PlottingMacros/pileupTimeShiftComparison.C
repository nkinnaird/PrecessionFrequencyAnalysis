// 3-31-20: Old macro for plots for scans over fit start time with different pileup time shifts applied.

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
// #include <TRatioPlot.h> // wait till root v 6_08

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;


void comparePars(std::vector<float> fSs, std::vector<TF1*> fullRatioFits1, std::vector<TF1*> fullRatioFits2, std::vector<TF1*> fullRatioFits3, std::vector<TF1*> fullRatioFits4, std::vector<TF1*> fullRatioFits5)
{
  string xAxisTitle = "Fit Start Time (#mus)";

  // auto FullratioDir = inDir->mkdir("FullRatio");

    TGraph* FullRatio_R_diff_1 = new TGraph();
    TGraph* FullRatio_R_diff_2 = new TGraph();
    TGraph* FullRatio_R_diff_3 = new TGraph();
    TGraph* FullRatio_R_diff_4 = new TGraph();
    TGraph* FullRatio_R_diff_5 = new TGraph();

    FullRatio_R_diff_2->SetLineColor(2);
    FullRatio_R_diff_3->SetLineColor(3);
    FullRatio_R_diff_4->SetLineColor(4);
    FullRatio_R_diff_5->SetLineColor(6);

    FullRatio_R_diff_2->SetMarkerColor(2);
    FullRatio_R_diff_3->SetMarkerColor(3);
    FullRatio_R_diff_4->SetMarkerColor(4);
    FullRatio_R_diff_5->SetMarkerColor(6);

    FullRatio_R_diff_2->SetMarkerSize(1);
    FullRatio_R_diff_3->SetMarkerSize(1);
    FullRatio_R_diff_4->SetMarkerSize(1);
    FullRatio_R_diff_5->SetMarkerSize(1);

/////////////////////////////////////////////////////////////////////////////////////

int pointNo = 0;

    for (float fS : fSs)
    {
      double FullRatio_R_1 = fullRatioFits1.at(pointNo)->GetParameter(1);
      double FullRatio_R_2 = fullRatioFits2.at(pointNo)->GetParameter(1);
      double FullRatio_R_3 = fullRatioFits3.at(pointNo)->GetParameter(1);
      double FullRatio_R_4 = fullRatioFits4.at(pointNo)->GetParameter(1);
      double FullRatio_R_5 = fullRatioFits5.at(pointNo)->GetParameter(1);

        FullRatio_R_diff_1->SetPoint(pointNo, fS, FullRatio_R_1 - FullRatio_R_1);
        FullRatio_R_diff_2->SetPoint(pointNo, fS, FullRatio_R_2 - FullRatio_R_1);
        FullRatio_R_diff_3->SetPoint(pointNo, fS, FullRatio_R_3 - FullRatio_R_1);
        FullRatio_R_diff_4->SetPoint(pointNo, fS, FullRatio_R_4 - FullRatio_R_1);
        FullRatio_R_diff_5->SetPoint(pointNo, fS, FullRatio_R_5 - FullRatio_R_1);

      pointNo++;
    } // end loop over point numbers


  

  // FullratioDir->cd();

        nsTOus(FullRatio_R_diff_1, xAxisTitle);
        nsTOus(FullRatio_R_diff_2, xAxisTitle);
        nsTOus(FullRatio_R_diff_3, xAxisTitle);
        nsTOus(FullRatio_R_diff_4, xAxisTitle);
        nsTOus(FullRatio_R_diff_5, xAxisTitle);

      FullRatio_R_diff_1->SetName("FullRatio_R_diff");
      // FullRatio_R_diff_1->SetTitle("#Delta Ratio CBO R (ppm) Vs Fit Start Time");
      FullRatio_R_diff_1->SetTitle("R (PU phase shifted) - R (no shift) Vs Fit Start Time");
      FullRatio_R_diff_1->GetXaxis()->SetTitle(xAxisTitle.c_str());
      FullRatio_R_diff_1->GetYaxis()->SetTitle("#Delta R (ppm)");
      FullRatio_R_diff_1->GetYaxis()->SetRangeUser(-.15, .15);

      auto FullRatio_R_diff_canv = new TCanvas("FullRatio_R_diff_canv","FullRatio_R_diff_canv",200,10,1200,800);

      FullRatio_R_diff_1->Draw();
      FullRatio_R_diff_2->Draw("SAMEPL");
      FullRatio_R_diff_3->Draw("SAMEPL");
      FullRatio_R_diff_4->Draw("SAMEPL");
      FullRatio_R_diff_5->Draw("SAMEPL");

/////////////////////////////////////////////////////////////////////////////////////

        auto legend = new TLegend(0.7,0.63,0.9,0.83);

          legend->AddEntry(FullRatio_R_diff_1, "0 ns phase shift","l");
          legend->AddEntry(FullRatio_R_diff_2, "+10 ns phase shift","l");
          legend->AddEntry(FullRatio_R_diff_3, "+5 ns phase shift","l");
          legend->AddEntry(FullRatio_R_diff_4, "-5 ns phase shift","l");
          legend->AddEntry(FullRatio_R_diff_5, "-10 ns phase shift","l");
  
        legend->SetBorderSize(0);
        legend->Draw();

/////////////////////////////////////////////////////////////////////////////////////

        FullRatio_R_diff_canv->SaveAs("pileupTimeShiftComparison.png");
      // FullRatio_R_diff_canv->Write();
      // delete FullRatio_R_diff_canv;


}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

// have to make sure separate files are fitted in the same way
// filepath 1 is the non-time shifted file
int pileupTimeShiftComparison()
{
  // gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  TFile *inputFile1 = TFile::Open("output-fs-0nsShift.root");
  TFile *inputFile2 = TFile::Open("output-fs-p10nsShift.root");
  TFile *inputFile3 = TFile::Open("output-fs-p5nsShift.root");
  TFile *inputFile4 = TFile::Open("output-fs-n5nsShift.root");
  TFile *inputFile5 = TFile::Open("output-fs-n10nsShift.root");
   if (inputFile1 == 0 || inputFile2 == 0 || inputFile3 == 0 || inputFile4 == 0 || inputFile5 == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  // TFile* outputFile = new TFile("pileupTimeShiftComparison.root","RECREATE");
  // auto topDir = outputFile->mkdir("topDir");


/////////////////////////////////////////////////////////////////////////////////////
  // These only get set for the interactive root session (any generated canvases, etc.), but does not apply to the output root file - that comes from .rootlogon.C
  gStyle->SetOptStat(000000);
  gStyle->SetOptTitle(0);
  // gStyle->SetOptFit(0000);
  gStyle->SetOptFit(1111);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1);
  gStyle->SetLineColor(1);

  gStyle->SetPadRightMargin(.05);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
      // Plots over different fit conditions

  // auto addedDir = topDir->mkdir("Added");
  // auto calosDir = topDir->mkdir("Calos");


std::vector<float> fSs;
std::vector<TF1*> fullRatioFits1;
std::vector<TF1*> fullRatioFits2;
std::vector<TF1*> fullRatioFits3;
std::vector<TF1*> fullRatioFits4;
std::vector<TF1*> fullRatioFits5;


      TNtuple *firstTuple = (TNtuple*)inputFile1->Get(Form("topDir/FitPasses/FitPass0/FitConditions0"));
      float tP;
      firstTuple->SetBranchAddress("totalPasses", &tP);
      firstTuple->GetEntry(0);

// perform fit start scan for added data
    for (int fitPass = 0; fitPass < tP; ++fitPass)
    {
      float fS;

      TNtuple *ntuple = (TNtuple*)inputFile1->Get(Form("topDir/FitPasses/FitPass%d/FitConditions%d", fitPass, fitPass));
      ntuple->SetBranchAddress("fitStartTime", &fS);
      ntuple->GetEntry(0);


      TF1* fullRatioFitFunction1 = (TF1*) ((TGraphErrors*) inputFile1->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass)))->GetFunction("fullRatioFitFunc")->Clone();
      TF1* fullRatioFitFunction2 = (TF1*) ((TGraphErrors*) inputFile2->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass)))->GetFunction("fullRatioFitFunc")->Clone();
      TF1* fullRatioFitFunction3 = (TF1*) ((TGraphErrors*) inputFile3->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass)))->GetFunction("fullRatioFitFunc")->Clone();
      TF1* fullRatioFitFunction4 = (TF1*) ((TGraphErrors*) inputFile4->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass)))->GetFunction("fullRatioFitFunc")->Clone();
      TF1* fullRatioFitFunction5 = (TF1*) ((TGraphErrors*) inputFile5->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass)))->GetFunction("fullRatioFitFunc")->Clone();

      fSs.push_back(fS);
      fullRatioFits1.push_back(fullRatioFitFunction1);
      fullRatioFits2.push_back(fullRatioFitFunction2);
      fullRatioFits3.push_back(fullRatioFitFunction3);
      fullRatioFits4.push_back(fullRatioFitFunction4);
      fullRatioFits5.push_back(fullRatioFitFunction5);

    }

  comparePars(fSs, fullRatioFits1, fullRatioFits2, fullRatioFits3, fullRatioFits4, fullRatioFits5);

/////////////////////////////////////////////////////////////////////////////////////


      // delete outputFile;


  return 1;
}
