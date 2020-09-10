// 3-31-20: Macro for plots for comparison between T method and Ratio method fit start scans.

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


int FitStartComparisonTR(std::string filePath)
{
  // gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

/////////////////////////////////////////////////////////////////////////////////////
  // These only get set for the interactive root session (any generated canvases, etc.), but does not apply to the output root file - that comes from .rootlogon.C
  gStyle->SetOptStat(000000);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(2);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1);
  gStyle->SetLineColor(1);

  gStyle->SetPadRightMargin(.05);
  gStyle->SetPadLeftMargin(.2);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  vector<float> fSs;
  
  vector<TF1*> TmethodFits;
  vector<TF1*> ratioFits;


      TNtuple *firstTuple = (TNtuple*)inputFile->Get(Form("topDir/FitPasses/FitPass0/FitConditions0"));
      float tP;
      firstTuple->SetBranchAddress("totalPasses", &tP);
      firstTuple->GetEntry(0);

    for (int fitPass = 0; fitPass < tP; ++fitPass)
    {
      float fS;

      TNtuple *ntuple = (TNtuple*)inputFile->Get(Form("topDir/FitPasses/FitPass%d/FitConditions%d", fitPass, fitPass));
      ntuple->SetBranchAddress("fitStartTime", &fS);
      ntuple->GetEntry(0);

      TF1* TmethodFitFunction = (TF1*) ((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/TMethod/allTimesAdded_TMethod", fitPass)))->GetFunction("TmethodFitFunc")->Clone();
      TF1* fullRatioFitFunction = (TF1*) ((TGraphErrors*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass)))->GetFunction("fullRatioFitFunc")->Clone();

      fSs.push_back(fS);

      TmethodFits.push_back(TmethodFitFunction);
      ratioFits.push_back(fullRatioFitFunction);
    }

/////////////////////////////////////////////////////////////////////////////////////


  string xAxisTitle = "Fit Start Time [#mus]";


    TGraph* T_R_Diff = new TGraph();

    T_R_Diff->SetMarkerSize(0.75);
    T_R_Diff->SetLineColor(1);
    T_R_Diff->SetMarkerColor(1);

    TGraph* T_val = new TGraph();
    TGraph* R_val = new TGraph();

    T_val->SetMarkerColor(1);
    R_val->SetMarkerColor(2);

/////////////////////////////////////////////////////////////////////////////////////

int pointNo = 0;
double T_Rval_first = 0, R_Rval_first = 0;

    for (float fS : fSs)
    {
      double TMethod_R = TmethodFits.at(pointNo)->GetParameter(3);
      double FullRatio_R = ratioFits.at(pointNo)->GetParameter(1);

      if(pointNo == 0){
        T_Rval_first = TMethod_R;
        R_Rval_first = FullRatio_R;
      }

      T_R_Diff->SetPoint(pointNo, fS, 1000. * (TMethod_R - FullRatio_R));

      T_val->SetPoint(pointNo, fS, 1000. * (TMethod_R - T_Rval_first));
      R_val->SetPoint(pointNo, fS, 1000. * (FullRatio_R - T_Rval_first));

      pointNo++;
    } // end loop over point numbers


/////////////////////////////////////////////////////////////////////////////////////

        nsTOus(T_R_Diff, xAxisTitle);
        nsTOus(T_val, xAxisTitle);
        nsTOus(R_val, xAxisTitle);


      T_R_Diff->SetName("temp");
      T_R_Diff->SetTitle("#DeltaR (T-R) Vs Fit Start Time");
      T_R_Diff->GetXaxis()->SetTitle(xAxisTitle.c_str());
      T_R_Diff->GetYaxis()->SetTitle("#Delta R (T-R) [ppb]");
      T_R_Diff->GetYaxis()->SetRangeUser(-600, 600);

      auto TR_diff_compCanvas = new TCanvas("TR_diff_compCanvas","TR_diff_compCanvas",200,10,1200,800);

      T_R_Diff->Draw();

      TR_diff_compCanvas->Update();

      TLine *line = new TLine(TR_diff_compCanvas->GetUxmin(), 0, TR_diff_compCanvas->GetUxmax(), 0);
      line->SetLineColor(1);
      line->SetLineStyle(2);
      line->SetLineWidth(3);
      line->Draw();


        // auto legend_diff = new TLegend(0.7,0.675,0.975,0.975);

        //   legend_diff->AddEntry(T_R_Diff, "T-R (IFG 1 Multiplier) ","p");
  
        // legend_diff->SetBorderSize(1);
        // legend_diff->SetFillStyle(1001);
        // legend_diff->Draw();

      TR_diff_compCanvas->SaveAs("TR_Diff.png");

/////////////////////////////////////////////////////////////////////////////////////


      T_val->SetName("temp");
      T_val->SetTitle("#DeltaR (T & R) Vs Fit Start Time");
      T_val->GetXaxis()->SetTitle(xAxisTitle.c_str());
      T_val->GetYaxis()->SetTitle("#Delta R (T & R) [ppb]");
      // T_val->GetYaxis()->SetRangeUser(-600, 600);

      auto TR_val_canv = new TCanvas("TR_val_canv","TR_val_canv",200,10,1200,800);

      T_val->Draw();
      R_val->Draw("SAMEPL");

      TR_val_canv->Update();

      // TLine *line = new TLine(TR_val_canv->GetUxmin(), 0, TR_val_canv->GetUxmax(), 0);
      // line->SetLineColor(1);
      // line->SetLineStyle(2);
      // line->SetLineWidth(3);
      line->Draw();


        auto legend_val = new TLegend(0.7,0.675,0.975,0.875);

          legend_val->AddEntry(T_val, "T Method", "p");
          legend_val->AddEntry(R_val, "R Method", "p");
  
        legend_val->SetBorderSize(1);
        legend_val->SetFillStyle(1001);
        legend_val->Draw();

      TR_val_canv->SaveAs("TandR.png");

/////////////////////////////////////////////////////////////////////////////////////

  return 1;
}
