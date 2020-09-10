// 3-31-20: Macro for plots comparing wiggle spectra between two different histogram files. Could be dusted off and updated to be more general.

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
#include <TVectorD.h>

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;


double fitFunc(double* x, double* par)
{
  double time = x[0];
  double funcValue = par[0] * exp(-time/par[1]) * ( 1 + par[2] * cos(par[3]*time + par[4]));
  return funcValue;
}



int CompareHists()
{
  string firstFile_string = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/Endgame/RandSeeds/addedHists-Endgame-RandSeeds.root";
  string secondFile_string = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/Endgame/RandSeeds/withAdHoc/addedHists-Endgame-RandSeeds-withAdHoc.root";


  TFile *firstFile = TFile::Open(firstFile_string.c_str());
  TFile *secondFile = TFile::Open(secondFile_string.c_str());
   if (firstFile == 0 || secondFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

/////////////////////////////////////////////////////////////////////////////////////
  // These only get set for the interactive root session (any generated canvases, etc.), but does not apply to the output root file - that comes from .rootlogon.C
  gStyle->SetOptStat(000000);
  gStyle->SetOptTitle(0);
  // gStyle->SetOptFit(2);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1);
  gStyle->SetLineColor(1);

  gStyle->SetPadRightMargin(.05);
  gStyle->SetPadLeftMargin(.2);

/////////////////////////////////////////////////////////////////////////////////////
  
  int totalHistogramIters = (*(TVectorD*) firstFile->Get("Iters"))[0];

  TH1F* averageHist_first = (TH1F*) ((TH1F*) firstFile->Get("topDir/Iter0/Added/Times_E_Threshold"))->Clone();
  TH1F* averageHist_second = (TH1F*) ((TH1F*) secondFile->Get("topDir/Iter0/Added/Times_E_Threshold"))->Clone();


  for (int iter = 1; iter < totalHistogramIters; ++iter)
  {
  	TH1F* timesThreshold_first = (TH1F*) ((TH1F*) firstFile->Get(Form("topDir/Iter%i/Added/Times_E_Threshold", iter)))->Clone();
  	averageHist_first->Add(timesThreshold_first);

  	TH1F* timesThreshold_second = (TH1F*) ((TH1F*) secondFile->Get(Form("topDir/Iter%i/Added/Times_E_Threshold", iter)))->Clone();
  	averageHist_second->Add(timesThreshold_second);
  }

  averageHist_first->Scale(1./totalHistogramIters);
  averageHist_second->Scale(1./totalHistogramIters);



  averageHist_second->Add(averageHist_first, -1);


  TF1* myFunc = new TF1("myFunc", fitFunc, 30000, 200000, 5);
  myFunc->SetLineColor(2);
  myFunc->SetNpx(10000);

  myFunc->SetParameter(0, 10000);
  myFunc->SetParameter(1, 32.2*1000.);
  myFunc->SetParameter(2, 0.2);
  myFunc->SetParameter(3, defaultWa);
  myFunc->SetParameter(4, 0);

  averageHist_second->Fit(myFunc, "R0");


  	auto canv_second = new TCanvas("canv_second", "canv_second", 200, 200, 1200, 800);
  	averageHist_second->Draw("HIST");
  	myFunc->Draw("SAME");





	return 1;
}
