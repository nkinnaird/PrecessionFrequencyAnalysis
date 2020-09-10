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

using namespace std;


int ToyMCCompareVWAmplitudes()
{

  TFile *inputFile = TFile::Open("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_03/srcs/gm2analyses/macros/RatioMacro/ToyMC/VW_MC/toyOutputVW-50iters.root");



  TGraphErrors* Tmethod_VW_amps_graph = new TGraphErrors();
  TGraphErrors* Rmethod_VW_amps_graph = new TGraphErrors();

  Tmethod_VW_amps_graph->SetLineColor(2);
  Tmethod_VW_amps_graph->SetMarkerColor(2);

  for (int iter = 0; iter < 50; ++iter)
  {
  	 TH1F* TmethodHist = (TH1F*) inputFile->Get(Form("topDir/ToyMC/Iter%i/TmethoFitHist", iter)); // typo
  	 TF1* TmethodFit = TmethodHist->GetFunction("TmethodFitFunc");

  	 TGraphErrors* ratioCBOGraph = (TGraphErrors*) inputFile->Get(Form("topDir/ToyMC/Iter%i/Toy_Ratio_CBO_Graph", iter));
  	 TF1* ratioCBOFit = ratioCBOGraph->GetFunction("fullRatioFitFunc");


  	 cout << "amp paramters: " << TmethodFit->GetParameter(11) << " " << ratioCBOFit->GetParameter(9) << endl;

  	 Tmethod_VW_amps_graph->SetPoint(iter, iter+1, TmethodFit->GetParameter(11));
  	 Tmethod_VW_amps_graph->SetPointError(iter, 0, TmethodFit->GetParError(11));

  	 Rmethod_VW_amps_graph->SetPoint(iter, iter+1, ratioCBOFit->GetParameter(9));
  	 Rmethod_VW_amps_graph->SetPointError(iter, 0, ratioCBOFit->GetParError(9));
  }

    auto myCanvas = new TCanvas("myCanvas","myCanvas",200,10,800,600);

    Tmethod_VW_amps_graph->Draw("AP");
    Rmethod_VW_amps_graph->Draw("PSAME");






	// TGraph* Tmethod_VW_amps = (TGraph*) inputFile->Get("topDir/Added/TMethod/Graphs/TMethod_A_VW_Vs_Iter");
	// TGraph* Rmethod_VW_amps = (TGraph*) inputFile->Get("topDir/Added/RatioCBO/Graphs/RatioCBO_A_VW_Vs_Iter");

	// double seedNumber, T_pointAmp, R_pointAmp, averageAmpRatio = 0;

	// for (int pointNo = 0; pointNo < Rmethod_VW_amps->GetN(); ++pointNo)
	// {
	// 	Tmethod_VW_amps->GetPoint(pointNo, seedNumber, T_pointAmp);
	// 	Rmethod_VW_amps->GetPoint(pointNo, seedNumber, R_pointAmp);
	// 	cout << "Seed: " << seedNumber << " R amp: " << R_pointAmp << " T amp: " << T_pointAmp << " ratio: " << R_pointAmp/T_pointAmp << endl;
	// 	averageAmpRatio += R_pointAmp/T_pointAmp;
	// }

	// cout << endl << "Average amplitude ratio: " << averageAmpRatio/Rmethod_VW_amps->GetN() << endl;


return 1;
}
