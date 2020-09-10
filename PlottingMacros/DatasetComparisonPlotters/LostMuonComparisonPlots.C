// 3-31-20: Macro to take some lost muon spectra from different datasets and plot them against each other.

#include <iostream>
#include <fstream>
#include <string>
#include <functional>
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
#include <TVectorD.h>

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

bool saveImages = true;

int LostMuonComparisonPlots()
{

  gStyle->SetOptStat(000000);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1);
  gStyle->SetLineColor(1);
  gStyle->SetPadRightMargin(.05);
  gStyle->SetPadLeftMargin(.15);


  TFile *file_60h = TFile::Open("/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/60h/LostMuonsFiles/Plots/lostMuonPlots-60h-MainCuts.root");
  TFile *file_HighKick = TFile::Open("/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/HighKick/LostMuonsFiles/lostMuonPlots.root");
  TFile *file_9d = TFile::Open("/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/9d/LostMuonsFiles/Plots/lostMuonPlots-9d-MainCuts.root");
  TFile *file_Endgame = TFile::Open("/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/Endgame/LostMuonsFiles/Plots/lostMuonPlots-Endgame-MainCuts.root");

  TFile *file_Endgame_NegativeSide = TFile::Open("/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/Endgame/LostMuonsFiles/Plots/lostMuonPlots-Endgame-MainCuts-NegativeSide.root");

   if (file_60h == 0 || file_9d == 0 || file_Endgame == 0 || file_HighKick == 0 || file_Endgame_NegativeSide == 0){
      printf("Error: cannot open file\n");
      return 0;
   }


   TH1F* fractionalLosses_60h = (TH1F*) file_60h->Get("Fundamental/FractionalLosses");
   TH1F* fractionalLosses_HighKick = (TH1F*) file_HighKick->Get("Fundamental/FractionalLosses");
   TH1F* fractionalLosses_9d = (TH1F*) file_9d->Get("Fundamental/FractionalLosses");
   TH1F* fractionalLosses_Endgame = (TH1F*) file_Endgame->Get("Fundamental/FractionalLosses");

   fractionalLosses_60h->Scale(100);
   fractionalLosses_HighKick->Scale(100);
   fractionalLosses_9d->Scale(100);
   fractionalLosses_Endgame->Scale(100);

   fractionalLosses_HighKick->SetLineColor(4);
   fractionalLosses_9d->SetLineColor(2);
   fractionalLosses_Endgame->SetLineColor(8);


   fractionalLosses_60h->SetLineWidth(3);
   fractionalLosses_HighKick->SetLineWidth(3);
   fractionalLosses_9d->SetLineWidth(3);
   fractionalLosses_Endgame->SetLineWidth(3);

  auto canv_fractional_losses = new TCanvas("canv_fractional_losses","canv_fractional_losses",400,200,600,450);

  fractionalLosses_Endgame->GetYaxis()->SetTitle("Fractional Losses (%)");

  fractionalLosses_Endgame->Draw("HIST");
  fractionalLosses_60h->Draw("HISTSAME");
  fractionalLosses_HighKick->Draw("HISTSAME");
  fractionalLosses_9d->Draw("HISTSAME");

    auto legend_losses = new TLegend(0.2,0.6,.5,0.8);
    legend_losses->AddEntry(fractionalLosses_60h, "60h", "l");
    legend_losses->AddEntry(fractionalLosses_HighKick, "HighKick", "l");
    legend_losses->AddEntry(fractionalLosses_9d, "9d", "l");
    legend_losses->AddEntry(fractionalLosses_Endgame, "Endgame", "l");
    legend_losses->SetBorderSize(0);
    legend_losses->Draw("SAME");


   if(saveImages) canv_fractional_losses->SaveAs("fractionalLosses_dataset_comparison.png");

/////////////////////////////////////////////////////////////////////////////////////

  auto canv_fractional_losses_negative_side = new TCanvas("canv_fractional_losses_negative_side","canv_fractional_losses_negative_side",1000,200,600,450);

  fractionalLosses_Endgame->Draw("HIST");

  TH1F* fractionalLosses_Endgame_NegativeSide = (TH1F*) file_Endgame_NegativeSide->Get("Fundamental/FractionalLosses");
  fractionalLosses_Endgame_NegativeSide->Scale(100);

  // fractionalLosses_Endgame_NegativeSide->Scale((4.248 + 5*0.068) / 4.248);

  fractionalLosses_Endgame_NegativeSide->SetLineWidth(3);
  
  fractionalLosses_Endgame_NegativeSide->Draw("HISTSAME");

    auto legend_losses_neg = new TLegend(0.2,0.65,.7,0.8);
    legend_losses_neg->AddEntry(fractionalLosses_Endgame, "Endgame", "l");
    legend_losses_neg->AddEntry(fractionalLosses_Endgame_NegativeSide, "Endgame - negative side cut", "l");
    legend_losses_neg->SetBorderSize(0);
    legend_losses_neg->SetFillStyle(0);
    legend_losses_neg->Draw("SAME");

  if(saveImages) canv_fractional_losses_negative_side->SaveAs("fractionalLosses_Endgame_negativeSide_comparison.png");

/////////////////////////////////////////////////////////////////////////////////////


    // canv_fractional_losses->cd();

    // fractionalLosses_Endgame_NegativeSide->SetLineColor(8);
    // fractionalLosses_Endgame_NegativeSide->Draw("HISTSAME");



	return 1;

}
