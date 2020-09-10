// 3-31-20: Macro for plots for scans fit start time, with different IFG parameters used. This was used to see how the R value oscillates as a function of fit start time when the wrong IFG is used, and to compare T method and Ratio method fits.

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


// have to make sure separate files are fitted in the same way
int IFG_FitStartComparison()
{
  // gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  // TFile *inputFile_0scale = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/60h/Gain/IFG/Amplitude/DifferentFitStartTimes/FitStartScans/output-IFG-FS-0scale.root");
  // TFile *inputFile_1scale = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/60h/Gain/IFG/Amplitude/DifferentFitStartTimes/FitStartScans/output-IFG-FS-1scale.root");
  // TFile *inputFile_2scale = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/60h/Gain/IFG/Amplitude/DifferentFitStartTimes/FitStartScans/output-IFG-FS-2scale.root");

  // TFile *inputFile_0scale = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/60h/Gain/IFG/Amplitude/DifferentFitStartTimes/FitStartScans/MoreFits/output-IFG-FS-0scale-moreFits.root");
  // TFile *inputFile_1scale = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/60h/Gain/IFG/Amplitude/DifferentFitStartTimes/FitStartScans/MoreFits/output-IFG-FS-1scale-moreFits.root");
  // TFile *inputFile_2scale = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/60h/Gain/IFG/Amplitude/DifferentFitStartTimes/FitStartScans/MoreFits/output-IFG-FS-2scale-moreFits.root");

  TFile *inputFile_0scale = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/60h/Gain/IFG/Amplitude-with-AdHoc/FitStartScans/FinerScans/output-60h-FS-0scale-withAdHoc.root");
  TFile *inputFile_1scale = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/60h/Gain/IFG/Amplitude-with-AdHoc/FitStartScans/FinerScans/output-60h-FS-1scale-withAdHoc.root");
  TFile *inputFile_2scale = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/60h/Gain/IFG/Amplitude-with-AdHoc/FitStartScans/FinerScans/output-60h-FS-2scale-withAdHoc.root");

   if (inputFile_0scale == 0 || inputFile_1scale == 0 || inputFile_2scale == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  // TFile* outputFile = new TFile("pileupTimeShiftComparison.root","RECREATE");
  // auto topDir = outputFile->mkdir("topDir");


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
  
  vector<TF1*> TmethodFits_0scale;
  vector<TF1*> TmethodFits_1scale;
  vector<TF1*> TmethodFits_2scale;

  vector<TF1*> ratioFits_0scale;
  vector<TF1*> ratioFits_1scale;
  vector<TF1*> ratioFits_2scale;


      TNtuple *firstTuple = (TNtuple*)inputFile_0scale->Get(Form("topDir/FitPasses/FitPass0/FitConditions0"));
      float tP;
      firstTuple->SetBranchAddress("totalPasses", &tP);
      firstTuple->GetEntry(0);

    for (int fitPass = 0; fitPass < tP; ++fitPass)
    {
      float fS;

      TNtuple *ntuple = (TNtuple*)inputFile_0scale->Get(Form("topDir/FitPasses/FitPass%d/FitConditions%d", fitPass, fitPass));
      ntuple->SetBranchAddress("fitStartTime", &fS);
      ntuple->GetEntry(0);

      TF1* TmethodFitFunction_0scale = (TF1*) ((TH1F*) inputFile_0scale->Get(Form("topDir/FitPasses/FitPass%d/addedDir/TMethod/allTimesAdded_TMethod", fitPass)))->GetFunction("TmethodFitFunc")->Clone();
      TF1* TmethodFitFunction_1scale = (TF1*) ((TH1F*) inputFile_1scale->Get(Form("topDir/FitPasses/FitPass%d/addedDir/TMethod/allTimesAdded_TMethod", fitPass)))->GetFunction("TmethodFitFunc")->Clone();
      TF1* TmethodFitFunction_2scale = (TF1*) ((TH1F*) inputFile_2scale->Get(Form("topDir/FitPasses/FitPass%d/addedDir/TMethod/allTimesAdded_TMethod", fitPass)))->GetFunction("TmethodFitFunc")->Clone();

      TF1* fullRatioFitFunction_0scale = (TF1*) ((TGraphErrors*) inputFile_0scale->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass)))->GetFunction("fullRatioFitFunc")->Clone();
      TF1* fullRatioFitFunction_1scale = (TF1*) ((TGraphErrors*) inputFile_1scale->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass)))->GetFunction("fullRatioFitFunc")->Clone();
      TF1* fullRatioFitFunction_2scale = (TF1*) ((TGraphErrors*) inputFile_2scale->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass)))->GetFunction("fullRatioFitFunc")->Clone();

      fSs.push_back(fS);

      TmethodFits_0scale.push_back(TmethodFitFunction_0scale);
      TmethodFits_1scale.push_back(TmethodFitFunction_1scale);
      TmethodFits_2scale.push_back(TmethodFitFunction_2scale);

      ratioFits_0scale.push_back(fullRatioFitFunction_0scale);
      ratioFits_1scale.push_back(fullRatioFitFunction_1scale);
      ratioFits_2scale.push_back(fullRatioFitFunction_2scale);

    }

/////////////////////////////////////////////////////////////////////////////////////


  string xAxisTitle = "Fit Start Time [#mus]";

    TGraph* TMethod_R_diff_0_scale = new TGraph();
    TGraph* TMethod_R_diff_1_scale = new TGraph();
    TGraph* TMethod_R_diff_2_scale = new TGraph();

    TMethod_R_diff_0_scale->SetLineColor(4);
    TMethod_R_diff_1_scale->SetLineColor(1);
    TMethod_R_diff_2_scale->SetLineColor(2);

    TMethod_R_diff_0_scale->SetMarkerColor(4);
    TMethod_R_diff_1_scale->SetMarkerColor(1);
    TMethod_R_diff_2_scale->SetMarkerColor(2);

    TMethod_R_diff_0_scale->SetMarkerSize(0.75);
    TMethod_R_diff_1_scale->SetMarkerSize(0.75);
    TMethod_R_diff_2_scale->SetMarkerSize(0.75);


    TGraph* FullRatio_R_diff_0_scale = new TGraph();
    TGraph* FullRatio_R_diff_1_scale = new TGraph();
    TGraph* FullRatio_R_diff_2_scale = new TGraph();

    FullRatio_R_diff_0_scale->SetLineStyle(2);
    FullRatio_R_diff_1_scale->SetLineStyle(2);
    FullRatio_R_diff_2_scale->SetLineStyle(2);

    FullRatio_R_diff_0_scale->SetLineColor(4);
    FullRatio_R_diff_1_scale->SetLineColor(1);
    FullRatio_R_diff_2_scale->SetLineColor(2);

    FullRatio_R_diff_0_scale->SetMarkerColor(4);
    FullRatio_R_diff_1_scale->SetMarkerColor(1);
    FullRatio_R_diff_2_scale->SetMarkerColor(2);

    FullRatio_R_diff_0_scale->SetMarkerSize(0.75);
    FullRatio_R_diff_1_scale->SetMarkerSize(0.75);
    FullRatio_R_diff_2_scale->SetMarkerSize(0.75);

/////////////////////////////////////////////////////////////////////////////////////

    TGraph* T_R_Diff_0_scale = new TGraph();
    TGraph* T_R_Diff_1_scale = new TGraph();
    TGraph* T_R_Diff_2_scale = new TGraph();

    T_R_Diff_0_scale->SetMarkerSize(0.75);
    T_R_Diff_1_scale->SetMarkerSize(0.75);
    T_R_Diff_2_scale->SetMarkerSize(0.75);

    T_R_Diff_0_scale->SetLineColor(4);
    T_R_Diff_1_scale->SetLineColor(1);
    T_R_Diff_2_scale->SetLineColor(2);

    T_R_Diff_0_scale->SetMarkerColor(4);
    T_R_Diff_1_scale->SetMarkerColor(1);
    T_R_Diff_2_scale->SetMarkerColor(2);

    T_R_Diff_0_scale->SetLineStyle(2);
    T_R_Diff_2_scale->SetLineStyle(2);


    TGraph* T_val_1_scale = new TGraph();
    TGraph* R_val_1_scale = new TGraph();

    T_val_1_scale->SetMarkerColor(1);
    R_val_1_scale->SetMarkerColor(2);

/////////////////////////////////////////////////////////////////////////////////////

int pointNo = 0;
double T_Rval_first = 0, R_Rval_first = 0;

    for (float fS : fSs)
    {
      double TMethod_R_0scale = TmethodFits_0scale.at(pointNo)->GetParameter(3);
      double TMethod_R_1scale = TmethodFits_1scale.at(pointNo)->GetParameter(3);
      double TMethod_R_2scale = TmethodFits_2scale.at(pointNo)->GetParameter(3);

        TMethod_R_diff_0_scale->SetPoint(pointNo, fS, 1000. * (TMethod_R_0scale - TMethod_R_1scale)); // 1000 multiple to put in units of ppb
        TMethod_R_diff_1_scale->SetPoint(pointNo, fS, 1000. * (TMethod_R_1scale - TMethod_R_1scale));
        TMethod_R_diff_2_scale->SetPoint(pointNo, fS, 1000. * (TMethod_R_2scale - TMethod_R_1scale));


      double FullRatio_R_0scale = ratioFits_0scale.at(pointNo)->GetParameter(1);
      double FullRatio_R_1scale = ratioFits_1scale.at(pointNo)->GetParameter(1);
      double FullRatio_R_2scale = ratioFits_2scale.at(pointNo)->GetParameter(1);

        FullRatio_R_diff_0_scale->SetPoint(pointNo, fS, 1000. * (FullRatio_R_0scale - FullRatio_R_1scale));
        FullRatio_R_diff_1_scale->SetPoint(pointNo, fS, 1000. * (FullRatio_R_1scale - FullRatio_R_1scale));
        FullRatio_R_diff_2_scale->SetPoint(pointNo, fS, 1000. * (FullRatio_R_2scale - FullRatio_R_1scale));

/////////////////////////////////////////////////////////////////////////////////////

        T_R_Diff_0_scale->SetPoint(pointNo, fS, 1000. * (TMethod_R_0scale - FullRatio_R_0scale));
        T_R_Diff_1_scale->SetPoint(pointNo, fS, 1000. * (TMethod_R_1scale - FullRatio_R_1scale));
        T_R_Diff_2_scale->SetPoint(pointNo, fS, 1000. * (TMethod_R_2scale - FullRatio_R_2scale));

/////////////////////////////////////////////////////////////////////////////////////

        if(pointNo == 0){
          T_Rval_first = TMethod_R_1scale;
          R_Rval_first = FullRatio_R_1scale;
        }

        T_val_1_scale->SetPoint(pointNo, fS, 1000. * (TMethod_R_1scale - T_Rval_first));
        R_val_1_scale->SetPoint(pointNo, fS, 1000. * (FullRatio_R_1scale - T_Rval_first));

/////////////////////////////////////////////////////////////////////////////////////

      pointNo++;
    } // end loop over point numbers


  
        nsTOus(TMethod_R_diff_0_scale, xAxisTitle);
        nsTOus(TMethod_R_diff_1_scale, xAxisTitle);
        nsTOus(TMethod_R_diff_2_scale, xAxisTitle);

        nsTOus(FullRatio_R_diff_0_scale, xAxisTitle);
        nsTOus(FullRatio_R_diff_1_scale, xAxisTitle);
        nsTOus(FullRatio_R_diff_2_scale, xAxisTitle);

      TMethod_R_diff_0_scale->SetName("temp");
      TMethod_R_diff_0_scale->SetTitle("#DeltaR Vs Fit Start Time");
      TMethod_R_diff_0_scale->GetXaxis()->SetTitle(xAxisTitle.c_str());
      TMethod_R_diff_0_scale->GetYaxis()->SetTitle("#Delta R [ppb]");
      TMethod_R_diff_0_scale->GetYaxis()->SetRangeUser(-400, 400);

      auto comparisonCanvas = new TCanvas("comparisonCanvas","comparisonCanvas",200,10,1200,800);

      TMethod_R_diff_0_scale->Draw();
      TMethod_R_diff_1_scale->Draw("SAMEPL");
      TMethod_R_diff_2_scale->Draw("SAMEPL");

      FullRatio_R_diff_0_scale->Draw("SAMEPL");
      FullRatio_R_diff_1_scale->Draw("SAMEPL");
      FullRatio_R_diff_2_scale->Draw("SAMEPL");

        auto legend = new TLegend(0.7,0.675,0.975,0.975);

          legend->AddEntry(TMethod_R_diff_0_scale, "T Method IFG 0 Amplitude Multiplier","l");
          legend->AddEntry(TMethod_R_diff_1_scale, "T Method IFG 1 Amplitude Multiplier","l");
          legend->AddEntry(TMethod_R_diff_2_scale, "T Method IFG 2 Amplitude Multiplier","l");

          legend->AddEntry(FullRatio_R_diff_0_scale, "Ratio IFG 0 Amplitude Multiplier","l");
          legend->AddEntry(FullRatio_R_diff_1_scale, "Ratio IFG 1 Amplitude Multiplier","l");
          legend->AddEntry(FullRatio_R_diff_2_scale, "Ratio IFG 2 Amplitude Multiplier","l");
  
        legend->SetBorderSize(1);
        legend->SetFillStyle(1001);
        legend->Draw();

/////////////////////////////////////////////////////////////////////////////////////

        comparisonCanvas->SaveAs("IFG_vs_FitStartTime.png");

/////////////////////////////////////////////////////////////////////////////////////

        nsTOus(T_R_Diff_0_scale, xAxisTitle);
        nsTOus(T_R_Diff_1_scale, xAxisTitle);
        nsTOus(T_R_Diff_2_scale, xAxisTitle);


        nsTOus(T_val_1_scale, xAxisTitle);
        nsTOus(R_val_1_scale, xAxisTitle);


      T_R_Diff_1_scale->SetName("temp");
      T_R_Diff_1_scale->SetTitle("#DeltaR (T-R) Vs Fit Start Time");
      T_R_Diff_1_scale->GetXaxis()->SetTitle(xAxisTitle.c_str());
      T_R_Diff_1_scale->GetYaxis()->SetTitle("#Delta R (T-R) [ppb]");
      T_R_Diff_1_scale->GetYaxis()->SetRangeUser(-600, 600);

      auto TR_diff_compCanvas = new TCanvas("TR_diff_compCanvas","TR_diff_compCanvas",200,10,1200,800);

      T_R_Diff_1_scale->Draw();
      // T_R_Diff_0_scale->Draw("SAMEPL");
      // T_R_Diff_2_scale->Draw("SAMEPL");

      // T_val_1_scale->Draw("SAMEPL");
      // R_val_1_scale->Draw("SAMEPL");

      TR_diff_compCanvas->Update();

      TLine *line = new TLine(TR_diff_compCanvas->GetUxmin(), 0, TR_diff_compCanvas->GetUxmax(), 0);
      line->SetLineColor(1);
      line->SetLineStyle(2);
      line->SetLineWidth(3);
      line->Draw();


        // auto legend_diff = new TLegend(0.7,0.675,0.975,0.975);

        //   legend_diff->AddEntry(T_R_Diff_1_scale, "T-R (IFG 1 Multiplier) ","p");
  
        // legend_diff->SetBorderSize(1);
        // legend_diff->SetFillStyle(1001);
        // legend_diff->Draw();

/////////////////////////////////////////////////////////////////////////////////////


      T_val_1_scale->SetName("temp");
      T_val_1_scale->SetTitle("#DeltaR (T & R) Vs Fit Start Time");
      T_val_1_scale->GetXaxis()->SetTitle(xAxisTitle.c_str());
      T_val_1_scale->GetYaxis()->SetTitle("#Delta R (T & R) [ppb]");
      // T_val_1_scale->GetYaxis()->SetRangeUser(-600, 600);

      auto TR_val_canv = new TCanvas("TR_val_canv","TR_val_canv",200,10,1200,800);

      T_val_1_scale->Draw();
      R_val_1_scale->Draw("SAMEPL");

      TR_val_canv->Update();

      // TLine *line = new TLine(TR_val_canv->GetUxmin(), 0, TR_val_canv->GetUxmax(), 0);
      // line->SetLineColor(1);
      // line->SetLineStyle(2);
      // line->SetLineWidth(3);
      line->Draw();


        auto legend_val = new TLegend(0.7,0.675,0.975,0.875);

          legend_val->AddEntry(T_val_1_scale, "T Method","p");
          legend_val->AddEntry(R_val_1_scale, "R Method","p");
  
        legend_val->SetBorderSize(1);
        legend_val->SetFillStyle(1001);
        legend_val->Draw();




/////////////////////////////////////////////////////////////////////////////////////

      // delete outputFile;


  return 1;
}
