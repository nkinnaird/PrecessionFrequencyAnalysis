#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TF2.h>
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
#include <TTree.h>
#include <THStack.h>
#include <TPaveStats.h>
#include <TVirtualHistPainter.h>



using namespace std;


int PileupScaleFactor()
{

  // TFile *inputFile = TFile::Open("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_04_00/srcs/gm2analyses/macros/RatioMacro/tempWork/pileupMove/output-8nsartificial-and-shadow.root");
  TFile *inputFile = TFile::Open("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_04_00/srcs/gm2analyses/macros/RatioMacro/tempWork/pileupMove/output.root");
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

   double lowE = 3100;
   double highE = 5000;


   TH1F* energiesTh = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/Input/Energy_Threshold")->Clone("energiesTh");
   TH1F* pileupEnergies = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/Pileup/added_pileupEnergies_shadow_threshold")->Clone("pileupEnergies");
  pileupEnergies->SetLineColor(2);

   for (int bin = 1; bin <= energiesTh->GetNbinsX(); ++bin)
   {
   	if(energiesTh->GetBinCenter(bin) < lowE || energiesTh->GetBinCenter(bin) > highE) energiesTh->SetBinContent(bin, 0);
   }
   for (int bin = 1; bin <= pileupEnergies->GetNbinsX(); ++bin)
   {
   	if(pileupEnergies->GetBinCenter(bin) < lowE || pileupEnergies->GetBinCenter(bin) > highE) pileupEnergies->SetBinContent(bin, 0);
   }

/////////////////////////////////////////////////////////////////////////////////////

  TH1F* dividedHist = (TH1F*) energiesTh->Clone("dividedHist");
  dividedHist->Divide(pileupEnergies);

  TF1* scaleFactor = new TF1("scaleFactor", "[0]", lowE, highE);
  dividedHist->Fit(scaleFactor, "R");

  auto canvas2 = new TCanvas("canvas2","canvas2",200,10,1200,1000);
  dividedHist->Draw("hist");
  scaleFactor->SetLineColor(2);
  scaleFactor->Draw("same");
  // canvas2->SaveAs("dividedHist.png");

/////////////////////////////////////////////////////////////////////////////////////

  gStyle->SetOptStat(0);

  double scaleNumber = scaleFactor->GetParameter(0);


  TH1F* scaledPileupEnergies = (TH1F*) pileupEnergies->Clone("scaledPileupEnergies");
  scaledPileupEnergies->Scale(scaleNumber);

  auto canvas3 = new TCanvas("canvas3","canvas3",200,10,1200,1000);
  energiesTh->SetTitle("Energy Comparison (scaled)");
  energiesTh->Draw();
  scaledPileupEnergies->Draw("histsame");

   auto legend3 = new TLegend(0.17,0.2,0.4,0.35);
   legend3->AddEntry(energiesTh,"Data (3250 - 5000 MeV)","l");
   legend3->AddEntry(scaledPileupEnergies, Form("Shadow pileup correction * %f", scaleNumber),"l");
   legend3->SetFillStyle(0);
   legend3->SetTextSize(0.018);
   legend3->Draw();

  // canvas3->SaveAs("scaledPileup.png");

/////////////////////////////////////////////////////////////////////////////////////

  auto canvas = new TCanvas("canvas","canvas",200,10,1200,1000);
  energiesTh->SetTitle("Energy Comparison");
  energiesTh->Draw();
  pileupEnergies->Draw("histsame");
  energiesTh->GetYaxis()->SetRangeUser(0, pileupEnergies->GetMaximum()*1.1);
  canvas->Update();

   auto legend1 = new TLegend(0.17,0.2,0.4,0.35);
   legend1->AddEntry(energiesTh,"Data (3250 - 5000 MeV)","l");
   legend1->AddEntry(pileupEnergies,"Shadow pileup correction","l");
   legend1->SetFillStyle(0);
   legend1->SetTextSize(0.018);
   legend1->Draw();


  // canvas->SaveAs("unscaledPileup.png");



  return 1;

 }
