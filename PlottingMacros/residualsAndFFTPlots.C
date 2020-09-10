// 3-31-20: Macro to make residual and associated FFT plots from fits. Is very similar to the residualPlots.hh file in the headers folder, but makes them a bit nicer and also makes some comparison plots between different types of fits.

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
#include "plotUtils.hh"


void drawLineOnCanv(double freq, TCanvas* canv, string label){

  canv->Update();

  TLine *line = new TLine(freq, canv->GetUymin(), freq, canv->GetUymax());
  line->SetLineColor(4);
  line->SetLineStyle(2);
  line->SetLineWidth(3);
  line->Draw();

  TText *t = new TText(freq, 1.01*canv->GetUymax(), label.c_str());
  t->SetTextAlign(12);
  t->SetTextColor(4);
  // t->SetTextFont(43);
  t->SetTextSize(.03);
  t->SetTextAngle(90);
  t->Draw();

  return;
}




int residualsAndFFTPlots(std::string filePath){

  // gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  string unneededString = "temp";
  applyDatasetTag(inputFile, unneededString);

/////////////////////////////////////////////////////////////////////////////////////

  TH1F* FiveParamFit_FFT = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/FiveParam/residual_fitRange/AddedTimes_fftHist_rescaled");
  TH1F* TMethodFit_FFT = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/TMethod/residual_fitRange/AddedTimes_fftHist_rescaled");
  TH1F* ThreeParamRatioFit_FFT = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/ThreeParamRatio/residual_fitRange/RatioAddedTimes_fftHist_rescaled");
  TH1F* ThreeParamRatioFit_FFT_earlyTimes = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/ThreeParamRatio/residual_otherRange/RatioAddedTimes_fftHist_rescaled");
  TH1F* FullRatioFit_FFT = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/FullRatio/residual_fitRange/FullRatioAddedTimes_fftHist_rescaled");

/////////////////////////////////////////////////////////////////////////////////////

  double FiveParam_max = FiveParamFit_FFT->GetMaximum();
  double TMethod_max = TMethodFit_FFT->GetMaximum();
  double ThreeParamRatio_max = ThreeParamRatioFit_FFT->GetMaximum();
  double ThreeParamRatio_max_earlyTimes = ThreeParamRatioFit_FFT_earlyTimes->GetMaximum();
  double FullRatio_max = FullRatioFit_FFT->GetMaximum();

  FiveParamFit_FFT->Scale(1./pow(10, int(log10(FiveParam_max))));
  // TMethodFit_FFT->Scale(1./pow(10, int(log10(FiveParam_max))));
  TMethodFit_FFT->Scale(1./pow(10, int(log10(TMethod_max))));
  ThreeParamRatioFit_FFT->Scale(1./pow(10, int(log10(ThreeParamRatio_max))));
  ThreeParamRatioFit_FFT_earlyTimes->Scale(1./pow(10, int(log10(ThreeParamRatio_max_earlyTimes))));
  FullRatioFit_FFT->Scale(1./pow(10, int(log10(FullRatio_max))));

  // FiveParam_max = FiveParamFit_FFT->GetMaximum();
  // TMethod_max = TMethodFit_FFT->GetMaximum();
  // FullRatio_max = FullRatioFit_FFT->GetMaximum();

/////////////////////////////////////////////////////////////////////////////////////

  // MHz
  double gm2Freq = 0.23;

  // 60h, Endgame
  double plot_cbo_freq = 0.37;
  double VWfreq = 2.29;

  // 9d, HighKick
  // double plot_cbo_freq = 0.415;
  // double VWfreq = 2.04;

/////////////////////////////////////////////////////////////////////////////////////

  // produce canvases for ratio cbo fit residuals and ffts

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(000000);
  gStyle->SetPadRightMargin(.05);

/////////////////////////////////////////////////////////////////////////////////////

  // make 5 parameter fit fft

  auto fftCanv_fivepar = new TCanvas("fftCanv_fivepar", "fftCanv_fivepar", 200, 10, 1000, 800);

  FiveParamFit_FFT->SetTitle("Five Parameter Fit Residual FFT");
  FiveParamFit_FFT->Draw("hist");

  drawLineOnCanv(plot_cbo_freq, fftCanv_fivepar, "cbo");
  drawLineOnCanv(gm2Freq, fftCanv_fivepar, "g-2");
  drawLineOnCanv(plot_cbo_freq-gm2Freq, fftCanv_fivepar, "cbo - g-2");
  drawLineOnCanv(plot_cbo_freq+gm2Freq, fftCanv_fivepar, "cbo + g-2");
  drawLineOnCanv(VWfreq, fftCanv_fivepar, "VW");

  fftCanv_fivepar->SaveAs(("Images/FFT_fiveParameter" + datasetTagForPlots + ".png").c_str());

  // make 3 param ratio fit fft

  auto fftCanv_threeParRatio = new TCanvas("fftCanv_threeParRatio", "fftCanv_threeParRatio", 200, 10, 1000, 800);

  ThreeParamRatioFit_FFT->SetTitle("Three Parameter Ratio Fit Residual FFT");
  ThreeParamRatioFit_FFT->Draw("hist");

  drawLineOnCanv(plot_cbo_freq, fftCanv_threeParRatio, "cbo");
  drawLineOnCanv(gm2Freq, fftCanv_threeParRatio, "g-2");
  drawLineOnCanv(plot_cbo_freq-gm2Freq, fftCanv_threeParRatio, "cbo - g-2");
  drawLineOnCanv(plot_cbo_freq+gm2Freq, fftCanv_threeParRatio, "cbo + g-2");
  drawLineOnCanv(VWfreq, fftCanv_threeParRatio, "VW");

  fftCanv_threeParRatio->SaveAs(("Images/FFT_threeParRatio" + datasetTagForPlots + ".png").c_str());

  // make 3 param ratio fit fft for early times

  auto fftCanv_threeParRatio_earlyTimes = new TCanvas("fftCanv_threeParRatio_earlyTimes", "fftCanv_threeParRatio_earlyTimes", 200, 10, 1000, 800);

  ThreeParamRatioFit_FFT_earlyTimes->SetTitle("Three Parameter Ratio Fit Residual FFT");
  ThreeParamRatioFit_FFT_earlyTimes->Draw("hist");

  drawLineOnCanv(plot_cbo_freq, fftCanv_threeParRatio_earlyTimes, "cbo");
  drawLineOnCanv(gm2Freq, fftCanv_threeParRatio_earlyTimes, "g-2");
  drawLineOnCanv(plot_cbo_freq-gm2Freq, fftCanv_threeParRatio_earlyTimes, "cbo - g-2");
  drawLineOnCanv(plot_cbo_freq+gm2Freq, fftCanv_threeParRatio_earlyTimes, "cbo + g-2");
  drawLineOnCanv(VWfreq, fftCanv_threeParRatio_earlyTimes, "VW");

  fftCanv_threeParRatio_earlyTimes->SaveAs(("Images/FFT_threeParRatio_earlyTimes" + datasetTagForPlots + ".png").c_str());

/////////////////////////////////////////////////////////////////////////////////////


  auto fullRatio_plots_canvas = new TCanvas("fullRatio_plots_canvas", "fullRatio_plots_canvas", 200, 10, 800, 800);

  TGraph* residualGraph = (TGraph*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/FullRatio/residual_fitRange/FullRatioAddedTimes_Residual");
  nsTOus(residualGraph, "Time (#mus)");
  residualGraph->SetTitle("Full Ratio Fit Residual");
  residualGraph->GetYaxis()->SetTitle("Fit Function - Ratio Value");
  residualGraph->Draw("AP");
  fullRatio_plots_canvas->SaveAs(("Images/fitResidual" + datasetTagForPlots + ".png").c_str());

  TGraph* pullGraph = (TGraph*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/FullRatio/residual_fitRange/FullRatioAddedTimes_Pull_Plot");
  nsTOus(pullGraph, "Time (#mus)");
  pullGraph->SetTitle("Full Ratio Fit Pull");
  pullGraph->GetYaxis()->SetTitle("Fit Residual / Fit Error");
  pullGraph->Draw("AP");
  fullRatio_plots_canvas->SaveAs(("Images/fitPull" + datasetTagForPlots + ".png").c_str());


  fullRatio_plots_canvas->SetCanvasSize(1000, 800);
  gStyle->SetOptStat("MR");
  gStyle->SetStatBorderSize(0);
  gStyle->SetStatW(0.37);
  gStyle->SetStatH(0.06);
  gStyle->SetStatX(0.94);
  gStyle->SetStatY(0.84);

  TH1F* pullHist = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/FullRatio/residual_fitRange/FullRatioAddedTimes_Projected_Pull");
  pullHist->SetTitle("Full Ratio Fit Projected Pull");
  pullHist->GetYaxis()->SetTitle("Entries");
  pullHist->GetXaxis()->SetTitle("Fit Residual / Fit Error");
  pullHist->Draw();

  fullRatio_plots_canvas->Modified();
  fullRatio_plots_canvas->SaveAs(("Images/fitPull_projected" + datasetTagForPlots + ".png").c_str());


  auto fftCanv = new TCanvas("fftCanv", "fftCanv", 200, 10, 1000, 800);
  gStyle->SetOptStat(0);

  FullRatioFit_FFT->SetTitle("Full Ratio Fit Residual FFT");
  FullRatioFit_FFT->Draw("hist");

  drawLineOnCanv(plot_cbo_freq, fftCanv, "cbo");
  drawLineOnCanv(gm2Freq, fftCanv, "g-2");
  drawLineOnCanv(plot_cbo_freq-gm2Freq, fftCanv, "cbo - g-2");
  drawLineOnCanv(plot_cbo_freq+gm2Freq, fftCanv, "cbo + g-2");
  drawLineOnCanv(VWfreq, fftCanv, "VW");

  fftCanv->SaveAs(("Images/FFT_fullRatioFit" + datasetTagForPlots + ".png").c_str());

/////////////////////////////////////////////////////////////////////////////////////

    // produce canvases for T method fit residuals and ffts

  auto Tmethod_plots_canvas = new TCanvas("Tmethod_plots_canvas", "Tmethod_plots_canvas", 200, 10, 800, 800);

  TGraph* residualGraph_TMethod = (TGraph*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/TMethod/residual_fitRange/AddedTimes_Residual");
  nsTOus(residualGraph_TMethod, "Time (#mus)");
  residualGraph_TMethod->SetTitle("T Method Fit Residual");
  residualGraph_TMethod->GetYaxis()->SetTitle("Fit Function - T Method Value");
  residualGraph_TMethod->Draw("AP");
  Tmethod_plots_canvas->SaveAs(("Images/fitResidual_TMethod" + datasetTagForPlots + ".png").c_str());

  TGraph* pullGraph_TMethod = (TGraph*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/TMethod/residual_fitRange/AddedTimes_Pull_Plot");
  nsTOus(pullGraph_TMethod, "Time (#mus)");
  pullGraph_TMethod->SetTitle("T Method Fit Pull");
  pullGraph_TMethod->GetYaxis()->SetTitle("Fit Residual / Fit Error");
  pullGraph_TMethod->Draw("AP");
  Tmethod_plots_canvas->SaveAs(("Images/fitPull_TMethod" + datasetTagForPlots + ".png").c_str());


  Tmethod_plots_canvas->SetCanvasSize(1000, 800);
  gStyle->SetOptStat("MR");
  gStyle->SetStatBorderSize(0);
  gStyle->SetStatW(0.37);
  gStyle->SetStatH(0.06);
  gStyle->SetStatX(0.94);
  gStyle->SetStatY(0.84);

  TH1F* pullHist_TMethod = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/TMethod/residual_fitRange/AddedTimes_Projected_Pull");
  pullHist_TMethod->SetTitle("T Method Fit Projected Pull");
  pullHist_TMethod->GetYaxis()->SetTitle("Entries");
  pullHist_TMethod->GetXaxis()->SetTitle("Fit Residual / Fit Error");
  pullHist_TMethod->Draw();

  Tmethod_plots_canvas->Modified();
  Tmethod_plots_canvas->SaveAs(("Images/fitPull_projected_TMethod" + datasetTagForPlots + ".png").c_str());


  auto fftCanv_TMethod = new TCanvas("fftCanv_TMethod", "fftCanv_TMethod", 200, 10, 1000, 800);
  gStyle->SetOptStat(0);

  TMethodFit_FFT->SetTitle("T Method Fit Residual FFT");
  TMethodFit_FFT->Draw("hist");

  drawLineOnCanv(plot_cbo_freq, fftCanv_TMethod, "cbo");
  drawLineOnCanv(gm2Freq, fftCanv_TMethod, "g-2");
  drawLineOnCanv(plot_cbo_freq-gm2Freq, fftCanv_TMethod, "cbo - g-2");
  drawLineOnCanv(plot_cbo_freq+gm2Freq, fftCanv_TMethod, "cbo + g-2");
  drawLineOnCanv(VWfreq, fftCanv_TMethod, "VW");

  fftCanv_TMethod->SaveAs(("Images/FFT_TMethodFit" + datasetTagForPlots + ".png").c_str());

/////////////////////////////////////////////////////////////////////////////////////
  
  // scale for comparison plots
  TMethodFit_FFT->Scale(pow(10, int(log10(TMethod_max))));
  TMethodFit_FFT->Scale(1./pow(10, int(log10(FiveParam_max))));

  TMethod_max = TMethodFit_FFT->GetMaximum();
  FullRatio_max = FullRatioFit_FFT->GetMaximum();

/////////////////////////////////////////////////////////////////////////////////////

  // make Residual FFT comparison plots after scaling appropriately

  auto Tmethod_Canv = new TCanvas("Tmethod_Canv", "Tmethod_Canv", 200, 10, 1000, 800);

  TMethodFit_FFT->SetLineColor(2);

  THStack* histStack_TMethod = new THStack("histStack_TMethod","Temp Title");
  histStack_TMethod->Add(FiveParamFit_FFT);
  histStack_TMethod->Add(TMethodFit_FFT);

  histStack_TMethod->SetTitle("Residual FFT Comparison");
  histStack_TMethod->Draw("nostack,hist");

  histStack_TMethod->GetXaxis()->SetTitle("Frequency [MHz]");
  histStack_TMethod->GetYaxis()->SetTitle("FFT Mag [arb.]");

  double userX = Tmethod_Canv->GetUxmax();
  double NDCx = GetNDCX(userX)-0.05;

  auto legend_TMethod = new TLegend(NDCx-.18,0.7,NDCx,0.8);
  legend_TMethod->AddEntry(FiveParamFit_FFT, "5 Parameter Fit", "l");
  legend_TMethod->AddEntry(TMethodFit_FFT, "T Method Fit", "l");

  legend_TMethod->SetTextSize(0.03);
  legend_TMethod->SetBorderSize(0);
  legend_TMethod->Draw();

  drawLineOnCanv(plot_cbo_freq, Tmethod_Canv, "cbo");
  drawLineOnCanv(gm2Freq, Tmethod_Canv, "g-2");
  drawLineOnCanv(plot_cbo_freq-gm2Freq, Tmethod_Canv, "cbo - g-2");
  drawLineOnCanv(plot_cbo_freq+gm2Freq, Tmethod_Canv, "cbo + g-2");
  drawLineOnCanv(VWfreq, Tmethod_Canv, "VW");

  Tmethod_Canv->Modified();
  Tmethod_Canv->SaveAs(("Images/FFTComparison_TMethod" + datasetTagForPlots + ".png").c_str());

/////////////////////////////////////////////////////////////////////////////////////

  auto FullRatio_Canv = new TCanvas("FullRatio_Canv", "FullRatio_Canv", 200, 10, 1000, 800);

  FullRatioFit_FFT->SetLineColor(2);
  FullRatioFit_FFT->Scale(TMethod_max/FullRatio_max);

  auto legend_FullRatio = new TLegend(NDCx-.18,0.7,NDCx,0.8);
  legend_FullRatio->AddEntry(FiveParamFit_FFT, "5 Parameter Fit", "l");
  legend_FullRatio->AddEntry(FullRatioFit_FFT, "Full Ratio Fit", "l");

  THStack* histStack_FullRatio = new THStack("histStack_FullRatio","Temp Title");
  histStack_FullRatio->Add(FiveParamFit_FFT);
  histStack_FullRatio->Add(FullRatioFit_FFT);

  histStack_FullRatio->SetTitle("Residual FFT Comparison");
  histStack_FullRatio->Draw("nostack,hist");

  histStack_FullRatio->GetXaxis()->SetTitle("Frequency [MHz]");
  histStack_FullRatio->GetYaxis()->SetTitle("FFT Mag [arb.]");

  legend_FullRatio->SetTextSize(0.03);
  legend_FullRatio->SetBorderSize(0);
  legend_FullRatio->Draw();

  drawLineOnCanv(plot_cbo_freq, FullRatio_Canv, "cbo");
  drawLineOnCanv(gm2Freq, FullRatio_Canv, "g-2");
  drawLineOnCanv(plot_cbo_freq-gm2Freq, FullRatio_Canv, "cbo - g-2");
  drawLineOnCanv(plot_cbo_freq+gm2Freq, FullRatio_Canv, "cbo + g-2");
  drawLineOnCanv(VWfreq, FullRatio_Canv, "VW");

  FullRatio_Canv->Modified();
  FullRatio_Canv->SaveAs(("Images/FFTComparison_FullRatio" + datasetTagForPlots + ".png").c_str());

/////////////////////////////////////////////////////////////////////////////////////

  return 1;

}
