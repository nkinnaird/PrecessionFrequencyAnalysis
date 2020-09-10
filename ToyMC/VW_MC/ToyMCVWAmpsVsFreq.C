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

#include "ratioAnalysisDefs.hh" // include paths for potential future grid jobs

using namespace std;


int ToyMCVWAmpsVsFreq(std::string filePath)
{

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

cout << "Number of iterations is hardcoded here - make sure it's the right number." << endl;


  TGraph* Tmethod_pVal_graph = new TGraph();
  TGraph* Rmethod_pVal_graph = new TGraph();

  TGraphErrors* Tmethod_VW_amps_graph = new TGraphErrors();
  TGraphErrors* Rmethod_VW_amps_graph = new TGraphErrors();

  Tmethod_pVal_graph->SetLineColor(2);
  Tmethod_VW_amps_graph->SetLineColor(2);
  Tmethod_VW_amps_graph->SetMarkerColor(2);

  for (int iter = 0; iter < 41; ++iter)
  {
  	 TH1F* TmethodHist = (TH1F*) inputFile->Get(Form("topDir/ToyMC/Iter%i/TmethodFitHist", iter));
  	 TF1* TmethodFit = TmethodHist->GetFunction("TmethodFitFunc");

  	 TGraphErrors* ratioCBOGraph = (TGraphErrors*) inputFile->Get(Form("topDir/ToyMC/Iter%i/Toy_Ratio_CBO_Graph", iter));
  	 TF1* ratioCBOFit = ratioCBOGraph->GetFunction("fullRatioFitFunc");

/////////////////////////////////////////////////////////////////////////////////////

     TF1* inputFunc = (TF1*) inputFile->Get(Form("topDir/ToyMC/Iter%i/truthFunc", iter));

     double funcVWFreq = inputFunc->GetParameter(9);

  	 Tmethod_VW_amps_graph->SetPoint(iter, funcVWFreq/(blindingWa*1e3), TmethodFit->GetParameter(10));
  	 Tmethod_VW_amps_graph->SetPointError(iter, 0, TmethodFit->GetParError(10));

  	 Rmethod_VW_amps_graph->SetPoint(iter, funcVWFreq/(blindingWa*1e3), ratioCBOFit->GetParameter(8));
  	 Rmethod_VW_amps_graph->SetPointError(iter, 0, ratioCBOFit->GetParError(8));

/////////////////////////////////////////////////////////////////////////////////////

     Tmethod_pVal_graph->SetPoint(iter, funcVWFreq/(blindingWa*1e3), TmethodFit->GetProb());
     Rmethod_pVal_graph->SetPoint(iter, funcVWFreq/(blindingWa*1e3), ratioCBOFit->GetProb());
  }

    auto myCanvas = new TCanvas("myCanvas","myCanvas",200,10,800,600);

    Rmethod_VW_amps_graph->GetXaxis()->SetTitle("#omega_{VW}/#omega_{a}");
    Rmethod_VW_amps_graph->GetYaxis()->SetTitle("A_{VW}");
    Rmethod_VW_amps_graph->Draw("AP");

    Tmethod_VW_amps_graph->Draw("PSAME");

    auto legend = new TLegend(0.2,0.6,0.4,0.8);
    legend->AddEntry(Rmethod_VW_amps_graph,"R Method Fit","p");
    legend->AddEntry(Tmethod_VW_amps_graph,"T Method Fit","p");

    legend->Draw("SAME");


/////////////////////////////////////////////////////////////////////////////////////

    auto myCanvas2 = new TCanvas("myCanvas2","myCanvas2",200,10,800,600);

    Rmethod_pVal_graph->GetXaxis()->SetTitle("#omega_{VW}/#omega_{a}");
    Rmethod_pVal_graph->GetYaxis()->SetTitle("P value");
    Rmethod_pVal_graph->Draw("AP");

    Tmethod_pVal_graph->Draw("PSAME");

    auto legend2 = new TLegend(0.2,0.6,0.4,0.8);
    legend2->AddEntry(Rmethod_pVal_graph,"R Method Fit","p");
    legend2->AddEntry(Tmethod_pVal_graph,"T Method Fit","p");

    legend2->Draw("SAME");


return 1;
}
