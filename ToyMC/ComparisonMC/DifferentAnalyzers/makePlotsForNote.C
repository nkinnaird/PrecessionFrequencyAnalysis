#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TVectorD.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveStats.h>
#include <TGraphErrors.h>
#include <TROOT.h>

// #include <plotUtils.hh>

using namespace std;

void drawGM2Line(TCanvas* canv){

  canv->Update();

  TLine *line = new TLine(0.23, canv->GetUymin(), 0.23, canv->GetUymax());
  line->SetLineColor(2);
  line->SetLineStyle(2);
  line->Draw();

  double labelY = 1.01 * canv->GetUymax();

  TText *t = new TText(0.23, labelY, "g-2");
  t->SetTextColor(2);
  t->SetTextSize(.03);
  t->SetTextAngle(90);
  t->Draw();

  return;
}


int makePlotsForNote(string filePath){

  // gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  TFile* inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

   gStyle->SetOptTitle(0);
   // gStyle->SetStatX(0.9);
   gStyle->SetPadRightMargin(0.08);


   // fit plots

   // T

   // single fit

   TH1D* single_Tmethod_fit = (TH1D*) inputFile->Get("Aaron/aaron_T_wiggle_0");
   single_Tmethod_fit->SetName("T-Method Fit");
   single_Tmethod_fit->GetFunction("fiveParamFit")->SetLineColor(2);
   single_Tmethod_fit->GetXaxis()->SetTitle("Time [ns]");
   single_Tmethod_fit->GetYaxis()->SetTitle("Entries");

   auto canv_T = new TCanvas("canv_T", "canv_T", 20, 20, 800, 600);
   single_Tmethod_fit->Draw("HIST");
   single_Tmethod_fit->GetFunction("fiveParamFit")->Draw("SAME");
   canv_T->SetLogy();
   canv_T->Update();

   TPaveStats* T_stats = (TPaveStats*) single_Tmethod_fit->GetListOfFunctions()->FindObject("stats");
   T_stats->SetBorderSize(1);
   T_stats->Draw("SAME");
   T_stats->SetOptStat(10);

   canv_T->SaveAs("Example_TMethod_Fit.png");

   // p values

   TH1D* pValues_TMethod = (TH1D*) inputFile->Get("Aaron/aaron_T_Pvalues");
   pValues_TMethod->SetName("T-Method p Values");
   pValues_TMethod->GetXaxis()->SetTitle("p Value");
   pValues_TMethod->GetYaxis()->SetTitle("Entries");
   pValues_TMethod->GetYaxis()->SetRangeUser(0, 1.1*pValues_TMethod->GetMaximum());

   auto canv_T_pValues = new TCanvas("canv_T_pValues", "canv_T_pValues", 20, 200, 800, 600);
   pValues_TMethod->Draw("HIST");
   canv_T_pValues->Update();

   TPaveStats* T_pValue_stats = (TPaveStats*) pValues_TMethod->GetListOfFunctions()->FindObject("stats");
   T_pValue_stats->SetBorderSize(1);
   T_pValue_stats->SetOptStat(1110);
   T_pValue_stats->SetX1NDC(0.75);
   T_pValue_stats->SetX2NDC(0.975);
   T_pValue_stats->SetY1NDC(0.7);
   T_pValue_stats->SetY2NDC(0.9);
   T_pValue_stats->Draw("SAME");

   canv_T_pValues->SaveAs("PValues_TMethod.png");

   // pull plot

   TH1D* Rpull_TMethod = (TH1D*) inputFile->Get("Aaron/aaron_T_RPull");
   Rpull_TMethod->SetName("T-Method R Pull");
   Rpull_TMethod->GetYaxis()->SetTitle("Entries");

   auto canv_T_Rpull = new TCanvas("canv_T_Rpull", "canv_T_Rpull", 20, 300, 800, 600);
   Rpull_TMethod->Draw("HIST");
   canv_T_Rpull->Update();

   TPaveStats* T_Rpull_stats = (TPaveStats*) Rpull_TMethod->GetListOfFunctions()->FindObject("stats");
   T_Rpull_stats->SetBorderSize(1);
   T_Rpull_stats->SetOptStat(2210);
   T_Rpull_stats->SetX1NDC(0.65);
   T_Rpull_stats->SetX2NDC(0.975);
   T_Rpull_stats->SetY1NDC(0.7);
   T_Rpull_stats->SetY2NDC(0.9);
   T_Rpull_stats->Draw("SAME");

   canv_T_Rpull->SaveAs("Rpull_TMethod.png");   

   // FFT

   TH1D* FFT_TMethod = (TH1D*) inputFile->Get("Aaron/aaron_T_wiggle_residual_fitRange_0/aaron_T_wiggle_fftHist_rescaled");
   FFT_TMethod->SetName("T-Method Fit FFT");
   FFT_TMethod->Scale(1/FFT_TMethod->Integral("width"));

   auto canv_T_FFT = new TCanvas("canv_T_FFT", "canv_T_FFT", 20, 400, 800, 600);
   FFT_TMethod->Draw("HIST");
   canv_T_FFT->Update();
   FFT_TMethod->SetStats(0);
   drawGM2Line(canv_T_FFT);

   canv_T_FFT->SaveAs("FFT_TMethod.png");   



   // A

   // single fit

   TH1D* single_Amethod_fit = (TH1D*) inputFile->Get("David/david_A_wiggle_0");
   single_Amethod_fit->SetName("A-Method Fit");
   single_Amethod_fit->GetFunction("fiveParamFit")->SetLineColor(2);
   single_Amethod_fit->GetXaxis()->SetTitle("Time [ns]");
   single_Amethod_fit->GetYaxis()->SetTitle("Entries (Asymmetry Weighted)");

   auto canv_A = new TCanvas("canv_A", "canv_A", 200, 20, 800, 600);
   single_Amethod_fit->Draw("HIST");
   single_Amethod_fit->GetFunction("fiveParamFit")->Draw("SAME");
   canv_A->SetLogy();
   canv_A->Update();

   TPaveStats* A_stats = (TPaveStats*) single_Amethod_fit->GetListOfFunctions()->FindObject("stats");
   A_stats->SetBorderSize(1);
   A_stats->Draw("SAME");
   A_stats->SetOptStat(10);

   canv_A->SaveAs("Example_AMethod_Fit.png");

   // p values

   TH1D* pValues_AMethod = (TH1D*) inputFile->Get("David/david_A_Pvalues");
   pValues_AMethod->SetName("A-Method p Values");
   pValues_AMethod->GetXaxis()->SetTitle("p Value");
   pValues_AMethod->GetYaxis()->SetTitle("Entries");
   pValues_AMethod->GetYaxis()->SetRangeUser(0, 1.1*pValues_AMethod->GetMaximum());

   auto canv_A_pValues = new TCanvas("canv_A_pValues", "canv_A_pValues", 200, 200, 800, 600);
   pValues_AMethod->Draw("HIST");
   canv_A_pValues->Update();

   TPaveStats* A_pValue_stats = (TPaveStats*) pValues_AMethod->GetListOfFunctions()->FindObject("stats");
   A_pValue_stats->SetBorderSize(1);
   A_pValue_stats->SetOptStat(1110);
   A_pValue_stats->SetX1NDC(0.75);
   A_pValue_stats->SetX2NDC(0.975);
   A_pValue_stats->SetY1NDC(0.7);
   A_pValue_stats->SetY2NDC(0.9);
   A_pValue_stats->Draw("SAME");

   canv_A_pValues->SaveAs("PValues_AMethod.png");

   // pull plot

   TH1D* Rpull_AMethod = (TH1D*) inputFile->Get("David/david_A_RPull");
   Rpull_AMethod->SetName("A-Method R Pull");
   Rpull_AMethod->GetYaxis()->SetTitle("Entries");

   auto canv_A_Rpull = new TCanvas("canv_A_Rpull", "canv_A_Rpull", 200, 300, 800, 600);
   Rpull_AMethod->Draw("HIST");
   canv_A_Rpull->Update();

   TPaveStats* A_Rpull_stats = (TPaveStats*) Rpull_AMethod->GetListOfFunctions()->FindObject("stats");
   A_Rpull_stats->SetBorderSize(1);
   A_Rpull_stats->SetOptStat(2210);
   A_Rpull_stats->SetX1NDC(0.65);
   A_Rpull_stats->SetX2NDC(0.975);
   A_Rpull_stats->SetY1NDC(0.7);
   A_Rpull_stats->SetY2NDC(0.9);
   A_Rpull_stats->Draw("SAME");

   canv_A_Rpull->SaveAs("Rpull_AMethod.png");   

   // FFT

   TH1D* FFT_AMethod = (TH1D*) inputFile->Get("David/david_A_wiggle_residual_fitRange_0/david_A_wiggle_fftHist_rescaled");
   FFT_AMethod->SetName("A-Method Fit FFT");
   FFT_AMethod->Scale(1/FFT_AMethod->Integral("width"));

   auto canv_A_FFT = new TCanvas("canv_A_FFT", "canv_A_FFT", 200, 400, 800, 600);
   FFT_AMethod->Draw("HIST");
   canv_A_FFT->Update();
   FFT_AMethod->SetStats(0);
   drawGM2Line(canv_A_FFT);

   canv_A_FFT->SaveAs("FFT_AMethod.png");   


   // Q

   // single fit

   TH1D* single_Qmethod_fit = (TH1D*) inputFile->Get("Tim/tim_Q_wiggle_0");
   single_Qmethod_fit->SetName("Q-Method Fit");
   single_Qmethod_fit->GetFunction("fiveParamFit")->SetLineColor(2);
   single_Qmethod_fit->GetXaxis()->SetTitle("Time [ns]");
   single_Qmethod_fit->GetYaxis()->SetTitle("Entries (Energy Weighted)");

   auto canv_Q = new TCanvas("canv_Q", "canv_Q", 400, 20, 800, 600);
   single_Qmethod_fit->Draw("HIST");
   single_Qmethod_fit->GetFunction("fiveParamFit")->Draw("SAME");
   canv_Q->SetLogy();
   canv_Q->Update();

   TPaveStats* Q_stats = (TPaveStats*) single_Qmethod_fit->GetListOfFunctions()->FindObject("stats");
   Q_stats->SetBorderSize(1);
   Q_stats->Draw("SAME");
   Q_stats->SetOptStat(10);

   canv_Q->SaveAs("Example_QMethod_Fit.png");

   // p values

   TH1D* pValues_QMethod = (TH1D*) inputFile->Get("Tim/tim_Q_Pvalues");
   pValues_QMethod->SetName("Q-Method p Values");
   pValues_QMethod->GetXaxis()->SetTitle("p Value");
   pValues_QMethod->GetYaxis()->SetTitle("Entries");
   pValues_QMethod->GetYaxis()->SetRangeUser(0, 1.1*pValues_QMethod->GetMaximum());

   auto canv_Q_pValues = new TCanvas("canv_Q_pValues", "canv_Q_pValues", 400, 200, 800, 600);
   pValues_QMethod->Draw("HIST");
   canv_Q_pValues->Update();

   TPaveStats* Q_pValue_stats = (TPaveStats*) pValues_QMethod->GetListOfFunctions()->FindObject("stats");
   Q_pValue_stats->SetBorderSize(1);
   Q_pValue_stats->SetOptStat(1110);
   Q_pValue_stats->SetX1NDC(0.75);
   Q_pValue_stats->SetX2NDC(0.975);
   Q_pValue_stats->SetY1NDC(0.7);
   Q_pValue_stats->SetY2NDC(0.9);
   Q_pValue_stats->Draw("SAME");

   canv_Q_pValues->SaveAs("PValues_QMethod.png");

   // pull plot

   TH1D* Rpull_QMethod = (TH1D*) inputFile->Get("Tim/tim_Q_RPull");
   Rpull_QMethod->SetName("Q-Method R Pull");
   Rpull_QMethod->GetYaxis()->SetTitle("Entries");

   auto canv_Q_Rpull = new TCanvas("canv_Q_Rpull", "canv_Q_Rpull", 400, 300, 800, 600);
   Rpull_QMethod->Draw("HIST");
   canv_Q_Rpull->Update();

   TPaveStats* Q_Rpull_stats = (TPaveStats*) Rpull_QMethod->GetListOfFunctions()->FindObject("stats");
   Q_Rpull_stats->SetBorderSize(1);
   Q_Rpull_stats->SetOptStat(2210);
   Q_Rpull_stats->SetX1NDC(0.65);
   Q_Rpull_stats->SetX2NDC(0.975);
   Q_Rpull_stats->SetY1NDC(0.7);
   Q_Rpull_stats->SetY2NDC(0.9);
   Q_Rpull_stats->Draw("SAME");

   canv_Q_Rpull->SaveAs("Rpull_QMethod.png");   

   // FFT

   TH1D* FFT_QMethod = (TH1D*) inputFile->Get("Tim/tim_Q_wiggle_residual_fitRange_0/tim_Q_wiggle_fftHist_rescaled");
   FFT_QMethod->SetName("Q-Method Fit FFT");
   FFT_QMethod->Scale(1/FFT_QMethod->Integral("width"));

   auto canv_Q_FFT = new TCanvas("canv_Q_FFT", "canv_Q_FFT", 400, 400, 800, 600);
   FFT_QMethod->Draw("HIST");
   canv_Q_FFT->Update();
   FFT_QMethod->SetStats(0);
   drawGM2Line(canv_Q_FFT);

   canv_Q_FFT->SaveAs("FFT_QMethod.png");   



   // R

   // single fit

   gStyle->SetOptStat(0);

   TGraphErrors* single_Rmethod_fit = (TGraphErrors*) inputFile->Get("Nick/nick_R_graph_0");
   single_Rmethod_fit->SetName("R-Method Fit");
   single_Rmethod_fit->GetFunction("ratioFit")->SetLineColor(2);
   single_Rmethod_fit->GetXaxis()->SetTitle("Time [ns]");
   single_Rmethod_fit->GetYaxis()->SetTitle("Ratio");

   auto canv_R = new TCanvas("canv_R", "canv_R", 300, 20, 800, 600);
   single_Rmethod_fit->Draw("AP");
   single_Rmethod_fit->GetFunction("ratioFit")->Draw("SAME");
   canv_R->Update();

   TPaveStats* R_stats = (TPaveStats*) single_Rmethod_fit->GetListOfFunctions()->FindObject("stats");
   R_stats->SetBorderSize(1);
   R_stats->Draw("SAME");
   R_stats->SetOptStat(0);
   R_stats->SetOptFit(1111);
   R_stats->SetX1NDC(0.65);
   R_stats->SetX2NDC(0.975);
   R_stats->SetY1NDC(0.7);
   R_stats->SetY2NDC(0.9);

   canv_R->SaveAs("Example_RMethod_Fit.png");

   gStyle->SetOptStat(1);

   // p values

   TH1D* pValues_RMethod = (TH1D*) inputFile->Get("Nick/RatioSeeds/nick_R_Pvalues_0");
   pValues_RMethod->SetName("R-Method p Values");
   pValues_RMethod->GetXaxis()->SetTitle("p Value");
   pValues_RMethod->GetYaxis()->SetTitle("Entries");
   pValues_RMethod->GetYaxis()->SetRangeUser(0, 1.1*pValues_RMethod->GetMaximum());

   auto canv_R_pValues = new TCanvas("canv_R_pValues", "canv_R_pValues", 300, 200, 800, 600);
   pValues_RMethod->Draw("HIST");
   canv_R_pValues->Update();

   TPaveStats* R_pValue_stats = (TPaveStats*) pValues_RMethod->GetListOfFunctions()->FindObject("stats");
   R_pValue_stats->SetBorderSize(1);
   R_pValue_stats->SetOptStat(1110);
   R_pValue_stats->SetX1NDC(0.75);
   R_pValue_stats->SetX2NDC(0.975);
   R_pValue_stats->SetY1NDC(0.7);
   R_pValue_stats->SetY2NDC(0.9);
   R_pValue_stats->Draw("SAME");

   canv_R_pValues->SaveAs("PValues_RMethod.png");

   // pull plot

   TH1D* Rpull_RMethod = (TH1D*) inputFile->Get("Nick/RatioSeeds/nick_R_RPull_0");
   Rpull_RMethod->SetName("R-Method R Pull");
   Rpull_RMethod->GetYaxis()->SetTitle("Entries");

   auto canv_R_Rpull = new TCanvas("canv_R_Rpull", "canv_R_Rpull", 300, 300, 800, 600);
   Rpull_RMethod->Draw("HIST");
   canv_R_Rpull->Update();

   TPaveStats* R_Rpull_stats = (TPaveStats*) Rpull_RMethod->GetListOfFunctions()->FindObject("stats");
   R_Rpull_stats->SetBorderSize(1);
   R_Rpull_stats->SetOptStat(2210);
   R_Rpull_stats->SetX1NDC(0.65);
   R_Rpull_stats->SetX2NDC(0.975);
   R_Rpull_stats->SetY1NDC(0.7);
   R_Rpull_stats->SetY2NDC(0.9);
   R_Rpull_stats->Draw("SAME");

   canv_R_Rpull->SaveAs("Rpull_RMethod.png");   

   // FFT

   TH1D* FFT_RMethod = (TH1D*) inputFile->Get("Nick/nick_R_graph_residual_fitRange_0/nick_R_graph_fftHist_rescaled");
   FFT_RMethod->SetName("R-Method Fit FFT");
   FFT_RMethod->Scale(1/FFT_RMethod->Integral("width"));

   auto canv_R_FFT = new TCanvas("canv_R_FFT", "canv_R_FFT", 300, 400, 800, 600);
   FFT_RMethod->Draw("HIST");
   canv_R_FFT->Update();
   FFT_RMethod->SetStats(0);
   drawGM2Line(canv_R_FFT);

   canv_R_FFT->SaveAs("FFT_RMethod.png");


	return 1;
}
