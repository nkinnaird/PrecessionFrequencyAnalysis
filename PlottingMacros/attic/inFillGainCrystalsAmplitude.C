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
#include <sstream>
#include <TTree.h>
#include <TNtuple.h>
#include <TPaveStats.h>
#include <TVectorD.h>
#include <THStack.h>
#include <TLatex.h>

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////

int inFillGainCrystalsAmplitude()
{
  // older files

  // vector<string> amplitudeFiles = {"/gm2/data/users/nkinnaird/Ratio/9d-FinalProduction/SingleIteration/SingleFit/output-9d-SingleFit.root",
  //                                  "/gm2/data/users/nkinnaird/Ratio/9d-FinalProduction/Gain/Amp/Fits/output-9d-Gain-Amp-0.root", 
  //                                  "/gm2/data/users/nkinnaird/Ratio/9d-FinalProduction/Gain/Amp/Fits/output-9d-Gain-Amp-0p9.root", 
  //                                  "/gm2/data/users/nkinnaird/Ratio/9d-FinalProduction/Gain/Amp/Fits/output-9d-Gain-Amp-1p1.root", 
  //                                  "/gm2/data/users/nkinnaird/Ratio/9d-FinalProduction/Gain/Amp/Fits/output-9d-Gain-Amp-2.root"}; 

  // vector<string> tauFiles = {"/gm2/data/users/nkinnaird/Ratio/9d-FinalProduction/SingleIteration/SingleFit/output-9d-SingleFit.root",
  //                            "/gm2/data/users/nkinnaird/Ratio/9d-FinalProduction/Gain/Tau/Fits/output-9d-Gain-Tau-0p1.root", 
  //                            "/gm2/data/users/nkinnaird/Ratio/9d-FinalProduction/Gain/Tau/Fits/output-9d-Gain-Tau-0p9.root", 
  //                            "/gm2/data/users/nkinnaird/Ratio/9d-FinalProduction/Gain/Tau/Fits/output-9d-Gain-Tau-1p1.root", 
  //                            "/gm2/data/users/nkinnaird/Ratio/9d-FinalProduction/Gain/Tau/Fits/output-9d-Gain-Tau-2.root"}; 

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  // 60h Thesis
  // string ampDefaultFile = "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/60h/SingleIteration/SingleFits/output-60h-Thesis-SingleIter.root";
  // vector<string> amplitudeFiles = {"/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/60h/Gain/Amplitude/0/output-60h-Thesis-Gain-Amp-0.root",
  //                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/60h/Gain/Amplitude/0p25/output-60h-Thesis-Gain-Amp-0p25.root",
  //                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/60h/Gain/Amplitude/0p5/output-60h-Thesis-Gain-Amp-0p5.root",
  //                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/60h/Gain/Amplitude/0p75/output-60h-Thesis-Gain-Amp-0p75.root",
  //                                  ampDefaultFile,
  //                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/60h/Gain/Amplitude/1p25/output-60h-Thesis-Gain-Amp-1p25.root",
  //                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/60h/Gain/Amplitude/1p5/output-60h-Thesis-Gain-Amp-1p5.root",
  //                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/60h/Gain/Amplitude/1p75/output-60h-Thesis-Gain-Amp-1p75.root",
  //                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/60h/Gain/Amplitude/2/output-60h-Thesis-Gain-Amp-2.root"};

  // HighKick Thesis
  // string ampDefaultFile = "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/HighKick/SingleIteration/SingleFits/output-HighKick-Thesis-SingleIter-fixed-tau-cbo.root";
  // vector<string> amplitudeFiles = {"/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/HighKick/Gain/Amplitude/0/output-HighKick-Thesis-Gain-Amp-0.root",
  //                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/HighKick/Gain/Amplitude/0p25/output-HighKick-Thesis-Gain-Amp-0p25.root",
  //                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/HighKick/Gain/Amplitude/0p5/output-HighKick-Thesis-Gain-Amp-0p5.root",
  //                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/HighKick/Gain/Amplitude/0p75/output-HighKick-Thesis-Gain-Amp-0p75.root",
  //                                  ampDefaultFile,
  //                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/HighKick/Gain/Amplitude/1p25/output-HighKick-Thesis-Gain-Amp-1p25.root",
  //                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/HighKick/Gain/Amplitude/1p5/output-HighKick-Thesis-Gain-Amp-1p5.root",
  //                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/HighKick/Gain/Amplitude/1p75/output-HighKick-Thesis-Gain-Amp-1p75.root",
  //                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/HighKick/Gain/Amplitude/2/output-HighKick-Thesis-Gain-Amp-2.root"};


  // 9d Thesis
  // string ampDefaultFile = "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/9d/SingleIteration/SingleFits/output-9d-Thesis-SingleIter.root";
  // vector<string> amplitudeFiles = {"/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/9d/Gain/Amplitude/0/output-9d-Thesis-Gain-Amp-0.root",
  //                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/9d/Gain/Amplitude/0p25/output-9d-Thesis-Gain-Amp-0p25.root",
  //                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/9d/Gain/Amplitude/0p5/output-9d-Thesis-Gain-Amp-0p5.root",
  //                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/9d/Gain/Amplitude/0p75/output-9d-Thesis-Gain-Amp-0p75.root",
  //                                  ampDefaultFile,
  //                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/9d/Gain/Amplitude/1p25/output-9d-Thesis-Gain-Amp-1p25.root",
  //                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/9d/Gain/Amplitude/1p5/output-9d-Thesis-Gain-Amp-1p5.root",
  //                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/9d/Gain/Amplitude/1p75/output-9d-Thesis-Gain-Amp-1p75.root",
  //                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/9d/Gain/Amplitude/2/output-9d-Thesis-Gain-Amp-2.root"};

  // Endgame Thesis
  string ampDefaultFile = "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/Endgame/SingleIteration/SingleFits/output-Endgame-Thesis-SingleIter.root";
  vector<string> amplitudeFiles = {"/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/Endgame/Gain/Amplitude/0/output-Endgame-Thesis-Gain-Amp-0.root",
                                   "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/Endgame/Gain/Amplitude/0p25/output-Endgame-Thesis-Gain-Amp-0p25.root",
                                   "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/Endgame/Gain/Amplitude/0p5/output-Endgame-Thesis-Gain-Amp-0p5.root",
                                   "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/Endgame/Gain/Amplitude/0p75/output-Endgame-Thesis-Gain-Amp-0p75.root",
                                   ampDefaultFile,
                                   "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/Endgame/Gain/Amplitude/1p25/output-Endgame-Thesis-Gain-Amp-1p25.root",
                                   "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/Endgame/Gain/Amplitude/1p5/output-Endgame-Thesis-Gain-Amp-1p5.root",
                                   "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/Endgame/Gain/Amplitude/1p75/output-Endgame-Thesis-Gain-Amp-1p75.root",
                                   "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/Endgame/Gain/Amplitude/2/output-Endgame-Thesis-Gain-Amp-2.root"};
                                   // "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/Endgame/Gain/Amplitude/1p99/output-Endgame-Thesis-Gain-Amp-1p99.root"};


/////////////////////////////////////////////////////////////////////////////////////
  // These only get set for the interactive root session (any generated canvases, etc.), but does not apply to the output root file - that comes from .rootlogon.C
  gStyle->SetOptStat(000000);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(2);
  // gStyle->SetOptFit(11);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1);
  gStyle->SetLineColor(1);
  gStyle->SetPadRightMargin(.05);
  gStyle->SetPadLeftMargin(.15);

/////////////////////////////////////////////////////////////////////////////////////

  // vector<double> amplitudeFactors = {0, 0.5, 1, 1.5, 2};
  vector<double> amplitudeFactors = {0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2};

    TF1* lineFuncShort = new TF1("lineFuncShort", "[0] + (x * [1])", 0.5, 1.5);
    TF1* lineFuncLong = new TF1("lineFuncLong", "[0] + (x * [1])", -.2, 2.2);
    lineFuncShort->SetLineColor(2);
    lineFuncLong->SetLineColor(2);

/////////////////////////////////////////////////////////////////////////////////////

  TGraphErrors* R_vs_Amp_RMethod = new TGraphErrors();
  R_vs_Amp_RMethod->SetName("R_Vs_Gain_Amp_RMethod");
  R_vs_Amp_RMethod->SetTitle("Full Ratio Fit R Vs In-Fill Gain Amplitude Multiplier");
  R_vs_Amp_RMethod->GetXaxis()->SetTitle("In-Fill Gain Amplitude Multiplier");
  R_vs_Amp_RMethod->GetYaxis()->SetTitle("R (ppm)");
  R_vs_Amp_RMethod->GetYaxis()->SetTitleOffset(2.2);

  TGraphErrors* R_vs_Amp_TMethod = new TGraphErrors();
  R_vs_Amp_TMethod->SetName("R_Vs_Gain_Amp_TMethod");
  R_vs_Amp_TMethod->SetTitle("T Method Fit R Vs In-Fill Gain Amplitude Multiplier");
  R_vs_Amp_TMethod->GetXaxis()->SetTitle("In-Fill Gain Amplitude Multiplier");
  R_vs_Amp_TMethod->GetYaxis()->SetTitle("R (ppm)");
  R_vs_Amp_TMethod->GetYaxis()->SetTitleOffset(2.2);


  TGraphErrors* R_vs_Amp_compare_T = new TGraphErrors();
  R_vs_Amp_compare_T->SetName("R_Vs_Gain_Amp_compare_T");
  R_vs_Amp_compare_T->GetXaxis()->SetTitle("In-Fill Gain Amplitude Multiplier c_{A}");
  R_vs_Amp_compare_T->GetYaxis()->SetTitle("R - R_{1} (ppm)");
  R_vs_Amp_compare_T->GetYaxis()->SetTitleOffset(2.2);

  TGraphErrors* R_vs_Amp_compare_R = new TGraphErrors();
  R_vs_Amp_compare_R->SetName("R_Vs_Gain_Amp_compare_T");

/////////////////////////////////////////////////////////////////////////////////////

  TGraphErrors* Chi2_vs_Amp_RMethod = new TGraphErrors();
  Chi2_vs_Amp_RMethod->SetName("Chi2_Vs_Gain_Amp_RMethod");
  Chi2_vs_Amp_RMethod->SetTitle("Full Ratio Fit #chi^{2} Vs In-Fill Gain Amplitude Multiplier");
  Chi2_vs_Amp_RMethod->GetXaxis()->SetTitle("In-Fill Gain Amplitude Multiplier");
  Chi2_vs_Amp_RMethod->GetYaxis()->SetTitle("#chi^{2}");
  Chi2_vs_Amp_RMethod->GetYaxis()->SetTitleOffset(2.2);

  TGraphErrors* Chi2_vs_Amp_TMethod = new TGraphErrors();
  Chi2_vs_Amp_TMethod->SetName("Chi2_Vs_Gain_Amp_TMethod");
  Chi2_vs_Amp_TMethod->SetTitle("T Method Fit #chi^{2} Vs In-Fill Gain Amplitude Multiplier");
  Chi2_vs_Amp_TMethod->GetXaxis()->SetTitle("In-Fill Gain Amplitude Multiplier");
  Chi2_vs_Amp_TMethod->GetYaxis()->SetTitle("#chi^{2}");
  Chi2_vs_Amp_TMethod->GetYaxis()->SetTitleOffset(2.2);


  TGraphErrors* Chi2_vs_Amp_compare_T = new TGraphErrors();
  Chi2_vs_Amp_compare_T->SetName("Chi2_Vs_Gain_Amp_compare_T");
  Chi2_vs_Amp_compare_T->GetXaxis()->SetTitle("In-Fill Gain Amplitude Multiplier c_{A}");
  Chi2_vs_Amp_compare_T->GetYaxis()->SetTitle("#chi^{2} - #chi^{2}_{1}");
  Chi2_vs_Amp_compare_T->GetYaxis()->SetTitleOffset(2);

  TGraphErrors* Chi2_vs_Amp_compare_R = new TGraphErrors();
  Chi2_vs_Amp_compare_R->SetName("Chi2_Vs_Gain_Amp_compare_R");

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  TFile* tempFile = TFile::Open(ampDefaultFile.c_str()); 
  TF1* tempRatioFunc = (TF1*) ((TGraphErrors*) tempFile->Get("topDir/FitPasses/FitPass0/addedDir/FullRatio/Added_Times_Full_Ratio_Graph"))->GetFunction("fullRatioFitFunc")->Clone();
  TF1* tempTMethodFunc = (TF1*) ((TH1F*) tempFile->Get("topDir/FitPasses/FitPass0/addedDir/TMethod/allTimesAdded_TMethod"))->GetFunction("TmethodFitFunc")->Clone();
  double controlR_amp_RMethod = tempRatioFunc->GetParameter(1);
  double controlR_amp_TMethod = tempTMethodFunc->GetParameter(3);
  double controlR_chi2_RMethod = tempRatioFunc->GetChisquare();
  double controlR_chi2_TMethod = tempTMethodFunc->GetChisquare();
  tempFile->Close();

    for (uint pointNo = 0; pointNo < amplitudeFactors.size(); ++pointNo)
    {
      TFile* inputFile = TFile::Open(amplitudeFiles.at(pointNo).c_str()); 
        if(inputFile == 0){
          printf("Error: cannot open file\n");
          return 0;
        }


      TF1* fullRatioFitFunction = (TF1*) ((TGraphErrors*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/FullRatio/Added_Times_Full_Ratio_Graph"))->GetFunction("fullRatioFitFunc")->Clone();
      TF1* TMethodFitFunc = (TF1*) ((TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/TMethod/allTimesAdded_TMethod"))->GetFunction("TmethodFitFunc")->Clone();
      
      R_vs_Amp_RMethod->SetPoint(pointNo, amplitudeFactors.at(pointNo), fullRatioFitFunction->GetParameter(1));
        // R_vs_Amp_RMethod->SetPointError(pointNo, 0, fullRatioFitFunction->GetParError(1));
      Chi2_vs_Amp_RMethod->SetPoint(pointNo, amplitudeFactors.at(pointNo), fullRatioFitFunction->GetChisquare());

      R_vs_Amp_TMethod->SetPoint(pointNo, amplitudeFactors.at(pointNo), TMethodFitFunc->GetParameter(3));
        // R_vs_Amp_TMethod->SetPointError(pointNo, 0, TMethodFitFunc->GetParError(3));
      Chi2_vs_Amp_TMethod->SetPoint(pointNo, amplitudeFactors.at(pointNo), TMethodFitFunc->GetChisquare());

      cout << "Amp scan - multiplier: " << amplitudeFactors.at(pointNo) << " R Method R diff: " << fullRatioFitFunction->GetParameter(1) - controlR_amp_RMethod << " T Method R diff: " << TMethodFitFunc->GetParameter(3) - controlR_amp_TMethod << endl;


      R_vs_Amp_compare_T->SetPoint(pointNo, amplitudeFactors.at(pointNo), TMethodFitFunc->GetParameter(3) - controlR_amp_TMethod);
      R_vs_Amp_compare_R->SetPoint(pointNo, amplitudeFactors.at(pointNo), fullRatioFitFunction->GetParameter(1) - controlR_amp_RMethod);

      Chi2_vs_Amp_compare_T->SetPoint(pointNo, amplitudeFactors.at(pointNo), TMethodFitFunc->GetChisquare() - controlR_chi2_TMethod);
      Chi2_vs_Amp_compare_R->SetPoint(pointNo, amplitudeFactors.at(pointNo), fullRatioFitFunction->GetChisquare() - controlR_chi2_RMethod);

      inputFile->Close();
    }

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

      R_vs_Amp_RMethod->GetXaxis()->SetLimits(-0.2, 2.2);
      R_vs_Amp_RMethod->Fit(lineFuncLong, "QR"); // fit for the systematic effect on R

      auto ampCanvas_RMethod = new TCanvas("ampCanvas_RMethod","ampCanvas_RMethod",200,10,500,400);
      R_vs_Amp_RMethod->Draw("AP");
      ampCanvas_RMethod->Update();

      TPaveStats *statsBox_R_amp_RMethod = (TPaveStats*)ampCanvas_RMethod->GetPrimitive("stats");
      statsBox_R_amp_RMethod->SetBorderSize(1);
      statsBox_R_amp_RMethod->Draw("SAME");

      ampCanvas_RMethod->SaveAs("Images/R_Vs_IFG_Amp_crystals_RMethod.png");


      R_vs_Amp_TMethod->GetXaxis()->SetLimits(-0.2, 2.2);
      R_vs_Amp_TMethod->Fit(lineFuncLong, "QR"); // fit for the systematic effect on R

      auto ampCanvas_TMethod = new TCanvas("ampCanvas_TMethod","ampCanvas_TMethod",700,10,500,400);
      R_vs_Amp_TMethod->Draw("AP");
      ampCanvas_TMethod->Update();

      TPaveStats *statsBox_R_amp_TMethod = (TPaveStats*)ampCanvas_TMethod->GetPrimitive("stats");
      statsBox_R_amp_TMethod->SetBorderSize(1);
      statsBox_R_amp_TMethod->Draw("SAME");

      ampCanvas_TMethod->SaveAs("Images/R_Vs_IFG_Amp_crystals_TMethod.png");

/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

      Chi2_vs_Amp_RMethod->GetXaxis()->SetLimits(-0.2, 2.2);

        TF1* parabFunc = new TF1("parabFunc", "[2] * (x - [1]) * (x - [1]) + [0]", -0.2, 2.2);
        parabFunc->SetNpx(10000);
        parabFunc->SetParameter(0,4000);
        parabFunc->SetParameter(1,1);
        parabFunc->SetParameter(2,14);
        Chi2_vs_Amp_RMethod->Fit("parabFunc", "Q");
        Chi2_vs_Amp_RMethod->GetFunction("parabFunc")->SetLineColor(2);

      auto chi2_ampCanvas_RMethod = new TCanvas("chi2_ampCanvas_RMethod","chi2_ampCanvas_RMethod",200,410,500,400);
      Chi2_vs_Amp_RMethod->Draw("AP");
      chi2_ampCanvas_RMethod->Update();

      TPaveStats *statsBox_chi2_amp_RMethod = (TPaveStats*)chi2_ampCanvas_RMethod->GetPrimitive("stats");
      statsBox_chi2_amp_RMethod->SetBorderSize(1);
      statsBox_chi2_amp_RMethod->Draw("SAME");

      chi2_ampCanvas_RMethod->SaveAs("Images/Chi2_Vs_IFG_Amp_crystals_RMethod.png");


      Chi2_vs_Amp_TMethod->GetXaxis()->SetLimits(-0.2, 2.2);

        parabFunc->SetParameter(0,4000);
        parabFunc->SetParameter(1,1);
        parabFunc->SetParameter(2,14);
        Chi2_vs_Amp_TMethod->Fit("parabFunc", "Q");
        Chi2_vs_Amp_TMethod->GetFunction("parabFunc")->SetLineColor(2);

      auto chi2_ampCanvas_TMethod = new TCanvas("chi2_ampCanvas_TMethod","chi2_ampCanvas_TMethod",700,410,500,400);
      Chi2_vs_Amp_TMethod->Draw("AP");
      chi2_ampCanvas_TMethod->Update();

      TPaveStats *statsBox_chi2_amp_TMethod = (TPaveStats*)chi2_ampCanvas_TMethod->GetPrimitive("stats");
      statsBox_chi2_amp_TMethod->SetBorderSize(1);
      statsBox_chi2_amp_TMethod->Draw("SAME");

      chi2_ampCanvas_TMethod->SaveAs("Images/Chi2_Vs_IFG_Amp_crystals_TMethod.png");

/////////////////////////////////////////////////////////////////////////////////////

  gStyle->SetOptFit(0);

/////////////////////////////////////////////////////////////////////////////////////

      R_vs_Amp_compare_T->GetXaxis()->SetLimits(-0.2, 2.2);
        lineFuncLong->SetLineColor(1);
      R_vs_Amp_compare_T->Fit(lineFuncLong, "QR");
        lineFuncLong->SetLineColor(2);
      R_vs_Amp_compare_R->Fit(lineFuncLong, "QR");

      R_vs_Amp_compare_R->SetMarkerColor(2);

      auto ampCanvas_compare = new TCanvas("ampCanvas_compare","ampCanvas_compare",1300,10,500,400);
      R_vs_Amp_compare_T->Draw("AP");
      R_vs_Amp_compare_R->Draw("PSAME");

      TPaveText* textBox = new TPaveText(0.33,0.7,0.725,0.9, "NDC");
      textBox->AddText(Form("T Method #DeltaR/#Deltac_{A} = %0.4f", R_vs_Amp_compare_T->GetFunction("lineFuncLong")->GetParameter(1)));
      textBox->AddText(Form("R Method #DeltaR/#Deltac_{A} = %0.4f", R_vs_Amp_compare_R->GetFunction("lineFuncLong")->GetParameter(1)));
      textBox->SetBorderSize(1);
      textBox->SetFillStyle(1001);

       TList *listOfLines = textBox->GetListOfLines();
       TText *tconst = textBox->GetLineWith("R Method");
       tconst->SetTextColor(2);

      textBox->Draw("SAME");

      ampCanvas_compare->SaveAs("Images/R_Vs_IFG_Amp_crystals_compare.png");

/////////////////////////////////////////////////////////////////////////////////////

      Chi2_vs_Amp_compare_T->GetXaxis()->SetLimits(-0.2, 2.2);
        parabFunc->SetParameter(0,0);
        parabFunc->SetParameter(1,.5);
        parabFunc->SetParameter(2,14);
        parabFunc->SetLineColor(1);
      Chi2_vs_Amp_compare_T->Fit(parabFunc, "QR");

      Chi2_vs_Amp_compare_R->SetMarkerColor(2);

      auto chi2Canvas_compare = new TCanvas("chi2Canvas_compare","chi2Canvas_compare",1300,410,500,400);
      Chi2_vs_Amp_compare_T->Draw("AP");
      Chi2_vs_Amp_compare_R->Draw("PSAME");

    auto legend_chi2 = new TLegend(0.4,0.6,.7,0.8);
    legend_chi2->AddEntry(Chi2_vs_Amp_compare_T, "T Method", "p");
    legend_chi2->AddEntry(Chi2_vs_Amp_compare_R, "R Method", "p");
    legend_chi2->SetBorderSize(0);
    legend_chi2->Draw("SAME");

      chi2Canvas_compare->SaveAs("Images/Chi2_Vs_IFG_Amp_crystals_compare.png");


  return 1;
}
