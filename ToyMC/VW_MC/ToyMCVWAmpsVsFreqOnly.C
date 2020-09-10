#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TROOT.h>
#include <TImage.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

#include "ratioAnalysisDefs.hh" // include paths for potential future grid jobs

using namespace std;


int ToyMCVWAmpsVsFreqOnly(std::string filePath)
{
 
  // gStyle->SetOptStat(0);
  // gStyle->SetOptTitle(0);
  // gStyle->SetOptFit(2);
  // gStyle->SetMarkerStyle(20);
  // gStyle->SetMarkerColor(1);
  // gStyle->SetMarkerSize(1);
  // gStyle->SetLineColor(1);
  gStyle->SetPadRightMargin(.08);

  /////////////////////////////////////////////////////////////////////////////////////

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

cout << "Number of iterations is hardcoded here - make sure it's the right number." << endl;


  TGraph* Tmethod_pVal_graph = new TGraph();
  TGraph* Rmethod_pVal_graph = new TGraph();
  // TGraph* RstyleFunc_pVal_graph = new TGraph();

  TGraphErrors* Tmethod_VW_amps_graph = new TGraphErrors();
  TGraphErrors* Rmethod_VW_amps_graph = new TGraphErrors();
  // TGraphErrors* RstyleFunc_VW_amps_graph = new TGraphErrors();

  Tmethod_pVal_graph->SetMarkerColor(2);
  Tmethod_VW_amps_graph->SetLineColor(2);
  Tmethod_VW_amps_graph->SetMarkerColor(2);

  // RstyleFunc_pVal_graph->SetMarkerColor(4);
  // RstyleFunc_VW_amps_graph->SetLineColor(4);
  // RstyleFunc_VW_amps_graph->SetMarkerColor(4);

  // for (int iter = 0; iter < 21; ++iter)
  for (int iter = 0; iter < 37; ++iter)
  {
     TH1F* TmethodHist = (TH1F*) inputFile->Get(Form("topDir/ToyMC/Iter%i/5FitHist", iter));
     TF1* TmethodFit = TmethodHist->GetFunction("hist_fitFunc");

     TGraphErrors* ratioGraph = (TGraphErrors*) inputFile->Get(Form("topDir/ToyMC/Iter%i/Toy_Ratio_Graph", iter));
     TF1* ratioFit = ratioGraph->GetFunction("graph_fitFunc");

     // TGraphErrors* ratioStyleGraph = (TGraphErrors*) inputFile->Get(Form("topDir/ToyMC/Iter%i/Toy_Ratio_Graph_Full_Func", iter));
     // TF1* ratioStyleFit = ratioStyleGraph->GetFunction("ratioStyleFunction");

/////////////////////////////////////////////////////////////////////////////////////

     TF1* inputFunc = (TF1*) inputFile->Get(Form("topDir/ToyMC/Iter%i/truthFunc", iter));

     double funcVWFreq = inputFunc->GetParameter(2);

  	 Tmethod_VW_amps_graph->SetPoint(iter, funcVWFreq/(blindingWa*1e3), TmethodFit->GetParameter(3));
  	 Tmethod_VW_amps_graph->SetPointError(iter, 0, TmethodFit->GetParError(3));

  	 Rmethod_VW_amps_graph->SetPoint(iter, funcVWFreq/(blindingWa*1e3), ratioFit->GetParameter(1));
  	 Rmethod_VW_amps_graph->SetPointError(iter, 0, ratioFit->GetParError(1));

     // RstyleFunc_VW_amps_graph->SetPoint(iter, funcVWFreq/(blindingWa*1e3), ratioStyleFit->GetParameter(1));
     // RstyleFunc_VW_amps_graph->SetPointError(iter, 0, ratioStyleFit->GetParError(1));

/////////////////////////////////////////////////////////////////////////////////////

     Tmethod_pVal_graph->SetPoint(iter, funcVWFreq/(blindingWa*1e3), TmethodFit->GetProb());
     Rmethod_pVal_graph->SetPoint(iter, funcVWFreq/(blindingWa*1e3), ratioFit->GetProb());
     // RstyleFunc_pVal_graph->SetPoint(iter, funcVWFreq/(blindingWa*1e3), ratioStyleFit->GetProb());
  }

    auto myCanvas = new TCanvas("myCanvas","myCanvas",200,10,800,600);

    Rmethod_VW_amps_graph->GetXaxis()->SetTitle("#omega_{VW}/#omega_{a}");
    Rmethod_VW_amps_graph->GetYaxis()->SetTitle("Fitted  A_{VW}");
    // Rmethod_VW_amps_graph->GetYaxis()->SetRangeUser(0, 0.06);
    Rmethod_VW_amps_graph->GetYaxis()->SetRangeUser(0, 0.1);
    Rmethod_VW_amps_graph->Draw("AP");

    Tmethod_VW_amps_graph->Draw("PSAME");
    // RstyleFunc_VW_amps_graph->Draw("PSAME");

    auto legend = new TLegend(0.7,0.77,0.98,0.97);
    legend->AddEntry(Tmethod_VW_amps_graph,"5 Param Exp Func","p");
    legend->AddEntry(Rmethod_VW_amps_graph,"3 Param Ratio Func","p");
    // legend->AddEntry(RstyleFunc_VW_amps_graph,"Full Ratio Func","p");

    legend->Draw("SAME");

    myCanvas->SaveAs("Fitted_Avw_Vs_Wvw.png");


/////////////////////////////////////////////////////////////////////////////////////

    auto myCanvas2 = new TCanvas("myCanvas2","myCanvas2",200,10,800,600);

    Rmethod_pVal_graph->GetXaxis()->SetTitle("#omega_{VW}/#omega_{a}");
    Rmethod_pVal_graph->GetYaxis()->SetTitle("P value");
    Rmethod_pVal_graph->GetYaxis()->SetRangeUser(0, 1);
    Rmethod_pVal_graph->Draw("AP");

    Tmethod_pVal_graph->Draw("PSAME");
    // RstyleFunc_pVal_graph->Draw("PSAME");

    auto legend2 = new TLegend(0.7,0.77,0.98,0.97);
    legend2->AddEntry(Tmethod_pVal_graph,"5 Param Exp Func","p");
    legend2->AddEntry(Rmethod_pVal_graph,"3 Param Ratio Func","p");
    // legend2->AddEntry(RstyleFunc_pVal_graph,"Full Ratio Func","p");

    legend2->Draw("SAME");


return 1;
}
