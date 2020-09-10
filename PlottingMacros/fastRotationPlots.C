// 3-31-20: Macro for to make plots illustrating the fast rotation effect in data - specifically Figure 2.18 in my thesis. I don't remember which input file was used, but whatever it was it needed to have smaller bins than 149.2 ns.

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


int fastRotationPlots(std::string filePath)
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
  gStyle->SetOptTitle(1);
  gStyle->SetOptFit(1111);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1);
  gStyle->SetLineColor(1);
  gStyle->SetStatY(0.94);
  gStyle->SetPadRightMargin(.05);

/////////////////////////////////////////////////////////////////////////////////////


  auto timesHist_all = (TH1F*) inputFile->Get ("topDir/Iter0/Added/Times_E_Threshold");
  auto timesHist_calo = (TH1F*) inputFile->Get("topDir/Iter0/Calos/Calo1/Calo1_Times_E_Threshold");

  nsTOus(timesHist_all, "Time (#mus)");
  nsTOus(timesHist_calo, "Time (#mus)");


  auto myCanvas = new TCanvas("myCanvas","myCanvas",200,10,1000,800);
  myCanvas->Divide(2,2);

  myCanvas->cd(1);
  timesHist_calo->SetTitle("Calo 1 Times");
  timesHist_calo->GetXaxis()->SetRangeUser(5, 12);
  timesHist_calo->Draw();
  myCanvas->Update();


  myCanvas->cd(2);
  timesHist_calo->GetXaxis()->SetRangeUser(40, 47);
  timesHist_calo->Draw();
  myCanvas->Update();

  myCanvas->cd(3);
  timesHist_all->SetTitle("All Calo Times");
  timesHist_all->GetXaxis()->SetRangeUser(5, 12);
  timesHist_all->Draw();
  myCanvas->Update();

  myCanvas->cd(4);
  timesHist_all->GetXaxis()->SetRangeUser(40, 47);
  timesHist_all->Draw();
  myCanvas->Update();

  myCanvas->SaveAs("fastrotation.png");

/////////////////////////////////////////////////////////////////////////////////////

  auto myCanvasMore = new TCanvas("myCanvasMore","myCanvasMore",200,10,1000,800);

  TPad *pad1 = new TPad("pad1", "pad1", 0, .66, 1, 1);
  TPad *pad2 = new TPad("pad2", "pad2", 0, .33, .5, .66);
  TPad *pad3 = new TPad("pad3", "pad3", .5, .33, 1, .66);
  TPad *pad4 = new TPad("pad4", "pad4", 0, 0, .5, .33);
  TPad *pad5 = new TPad("pad5", "pad5", .5, 0, 1, .33);


  pad1->Draw();
  pad1->cd();
  timesHist_calo->SetTitle("Calo 1 Times");
  timesHist_calo->GetXaxis()->SetRangeUser(5, 50);
  timesHist_calo->Draw();
  myCanvasMore->Update();

  myCanvasMore->cd();
  pad2->Draw();
  pad2->cd();
  timesHist_calo->SetTitle("Calo 1 Times");
  timesHist_calo->GetXaxis()->SetRangeUser(5, 12);
  timesHist_calo->Draw();
  myCanvasMore->Update();

  myCanvasMore->cd();
  pad3->Draw();
  pad3->cd();
  timesHist_calo->GetXaxis()->SetRangeUser(40, 47);
  timesHist_calo->Draw();
  myCanvasMore->Update();

  myCanvasMore->cd();
  pad4->Draw();
  pad4->cd();
  timesHist_all->SetTitle("All Calo Times");
  timesHist_all->GetXaxis()->SetRangeUser(5, 12);
  timesHist_all->Draw();
  myCanvasMore->Update();

  myCanvasMore->cd();
  pad5->Draw();
  pad5->cd();
  timesHist_all->GetXaxis()->SetRangeUser(40, 47);
  timesHist_all->Draw();
  myCanvasMore->Update();

/////////////////////////////////////////////////////////////////////////////////////

  myCanvasMore->SaveAs("fastrotationMore.png");

/////////////////////////////////////////////////////////////////////////////////////

  return 1;

}
