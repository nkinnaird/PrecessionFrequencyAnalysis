// 3-31-20: Macro for various plots related to the lost muons.

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


int LostMuonPlots(std::string filePath)
{
  // gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(2);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1);
  gStyle->SetLineColor(1);
  gStyle->SetPadRightMargin(.08);


  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  // check if dataset tag exists in file and if so append it to the file name and write it to the file

  string outputFileName = "lostMuonPlots";
  TNamed* tag = applyDatasetTag(inputFile, outputFileName);

  TFile* outputFile = new TFile((outputFileName + ".root").c_str(),"RECREATE");
  if(tag) tag->Write();

/////////////////////////////////////////////////////////////////////////////////////

  // plots for quadruples

  auto quadruple_plots_dir = outputFile->mkdir("QuadruplePlots");
  quadruple_plots_dir->cd();

   TH1F* Losses_triple_minus_quads_accs = (TH1F*) inputFile->Get("topDir/AdditionalCuts/Triples/Losses_minus_quadruples_accidentals/triple_losses_spectra_minus_quadruples_accidentals");

   TH1F* Losses_quadruples = (TH1F*) inputFile->Get("topDir/AdditionalCuts/Quadruples/Losses/quadruple_losses_spectra_addedCuts");
   Losses_quadruples->SetLineColor(2);

   nsTOus(Losses_triple_minus_quads_accs, "time [#mus]");
   nsTOus(Losses_quadruples, "time [#mus]");

   auto canv_trips_vs_quads = new TCanvas("canv_trips_vs_quads","canv_trips_vs_quads",0,0,600,450);
   canv_trips_vs_quads->SetLogy();

   Losses_triple_minus_quads_accs->Draw("hist");
   Losses_quadruples->Draw("histsame");

   auto legend_trips_vs_quads = new TLegend(0.6,0.6,.9,0.8);
   legend_trips_vs_quads->AddEntry(Losses_triple_minus_quads_accs, "Triples", "l");
   legend_trips_vs_quads->AddEntry(Losses_quadruples, "Quadruples", "l");
   legend_trips_vs_quads->SetBorderSize(0);
   legend_trips_vs_quads->Draw("same");

   canv_trips_vs_quads->Write();
   canv_trips_vs_quads->SaveAs(("Images/Triples_vs_quadruples" + datasetTagForPlots + ".png").c_str());


   TH1F* triples_over_quadruples = (TH1F*) Losses_triple_minus_quads_accs->Clone("triples_over_quadruples");
   triples_over_quadruples->Divide(Losses_quadruples); // might want to do some type of rolling average here

   triples_over_quadruples->GetYaxis()->SetTitle("Triples / Quadruples");

   auto canv_trips_over_quads = new TCanvas("canv_trips_over_quads","canv_trips_over_quads",100,0,600,450);
   triples_over_quadruples->Draw("HIST");

/////////////////////////////////////////////////////////////////////////////////////

   // plots for accidentals

  auto accidental_plots_dir = outputFile->mkdir("AccidentalPlots");
  accidental_plots_dir->cd();

   TH1F* Losses_accidentals = (TH1F*) inputFile->Get("topDir/Info/TriplesInfo/AccidentalBackground/Accidental_Background_Average_1D");
   Losses_accidentals->SetLineColor(2);

   nsTOus(Losses_accidentals, "time [#mus]");

   auto canv_trips_vs_accs = new TCanvas("canv_trips_vs_accs","canv_trips_vs_accs",200,0,600,450);
   canv_trips_vs_accs->SetLogy();

   Losses_triple_minus_quads_accs->Draw("hist");
   Losses_accidentals->Draw("histsame");

   auto legend_trips_vs_accs = new TLegend(0.6,0.6,.9,0.8);
   legend_trips_vs_accs->AddEntry(Losses_triple_minus_quads_accs, "Triples", "l");
   legend_trips_vs_accs->AddEntry(Losses_accidentals, "Accidentals", "l");
   legend_trips_vs_accs->SetBorderSize(0);
   legend_trips_vs_accs->Draw("same");

   canv_trips_vs_accs->Write();
   canv_trips_vs_accs->SaveAs(("Images/Triples_vs_accidentals" + datasetTagForPlots + ".png").c_str());

   // 2d plot before accidental subtraction

   gStyle->SetPadRightMargin(.12);

   TH2F* deltaT12_beforeAccSubtraction = (TH2F*) inputFile->Get("topDir/Info/TriplesInfo/triple_deltaT12_vs_timeInFill");

   auto canv_deltaT12_beforeAccSubtraction = new TCanvas("canv_deltaT12_beforeAccSubtraction","canv_deltaT12_beforeAccSubtraction",400,0,600,450);
   canv_deltaT12_beforeAccSubtraction->SetLogz();

   nsTOus(deltaT12_beforeAccSubtraction, "time [#mus]");
   deltaT12_beforeAccSubtraction->Draw("colz");

   deltaT12_beforeAccSubtraction->GetYaxis()->SetRangeUser(2, 12);
   canv_deltaT12_beforeAccSubtraction->Update();

   canv_deltaT12_beforeAccSubtraction->Write();
   canv_deltaT12_beforeAccSubtraction->SaveAs(("Images/deltaT12_beforeAccSubtraction" + datasetTagForPlots + ".png").c_str());

   // 2d plot after accidental subtraction

   TH2F* deltaT12_afterAccSubtraction = (TH2F*) inputFile->Get("topDir/Info/TriplesInfo/AccidentalBackground/triple_deltaT12_vs_timeInFill_accBg_Subtracted");

   auto canv_deltaT12_afterAccSubtraction = new TCanvas("canv_deltaT12_afterAccSubtraction","canv_deltaT12_afterAccSubtraction",500,0,600,450);
   canv_deltaT12_afterAccSubtraction->SetLogz();

    for (int yBin = 1; yBin <= deltaT12_afterAccSubtraction->GetYaxis()->GetNbins(); ++yBin){
      for (int xBin = 1; xBin <= deltaT12_afterAccSubtraction->GetXaxis()->GetNbins(); ++xBin){
        if(deltaT12_afterAccSubtraction->GetBinContent(xBin, yBin) < 1) deltaT12_afterAccSubtraction->SetBinContent(xBin, yBin, 0); // set all bins less than 1 to 0 for visualization purposes
      }
    }

   nsTOus(deltaT12_afterAccSubtraction, "time [#mus]");
   deltaT12_afterAccSubtraction->Draw("colz");

   deltaT12_afterAccSubtraction->GetYaxis()->SetRangeUser(2, 12);
   canv_deltaT12_afterAccSubtraction->Update();

   canv_deltaT12_afterAccSubtraction->Write();
   canv_deltaT12_afterAccSubtraction->SaveAs(("Images/deltaT12_afterAccSubtraction" + datasetTagForPlots + ".png").c_str());

   gStyle->SetPadRightMargin(.08);

/////////////////////////////////////////////////////////////////////////////////////

   // plot with quadruples and accidentals

   outputFile->cd();

   auto canv_trips_all = new TCanvas("canv_trips_all","canv_trips_all",300,0,600,450);
   canv_trips_all->SetLogy();

   Losses_accidentals->SetLineColor(4);

   Losses_triple_minus_quads_accs->Draw("hist");
   Losses_quadruples->Draw("histsame");
   Losses_accidentals->Draw("histsame");

   auto legend_trips_vs_all = new TLegend(0.6,0.6,.9,0.8);
   legend_trips_vs_all->AddEntry(Losses_triple_minus_quads_accs, "Triples", "l");
   legend_trips_vs_all->AddEntry(Losses_quadruples, "Quadruples", "l");
   legend_trips_vs_all->AddEntry(Losses_accidentals, "Accidentals", "l");
   legend_trips_vs_all->SetBorderSize(0);
   legend_trips_vs_all->Draw("same");

   canv_trips_all->Write();
   canv_trips_all->SaveAs(("Images/Triples_vs_all" + datasetTagForPlots + ".png").c_str());

/////////////////////////////////////////////////////////////////////////////////////

   // deltaT13 plots

  auto deltaT13_plots_dir = outputFile->mkdir("DeltaT13");
  deltaT13_plots_dir->cd();

   // time in fill no cuts

   gStyle->SetPadRightMargin(.12);

   TH2F* deltaT13_timeInFill_noCuts = (TH2F*) inputFile->Get("topDir/Info/TriplesInfo/triple_deltaT13_vs_timeInFill");

   auto canv_deltaT13_timeInFill_noCuts = new TCanvas("canv_deltaT13_timeInFill_noCuts","canv_deltaT13_timeInFill_noCuts",600,0,600,450);
   canv_deltaT13_timeInFill_noCuts->SetLogz();

   nsTOus(deltaT13_timeInFill_noCuts, "time [#mus]");
   deltaT13_timeInFill_noCuts->Draw("colz");

   deltaT13_timeInFill_noCuts->GetYaxis()->SetRangeUser(8, 18);
   canv_deltaT13_timeInFill_noCuts->Update();

   canv_deltaT13_timeInFill_noCuts->Write();
   canv_deltaT13_timeInFill_noCuts->SaveAs(("Images/deltaT13_timeInFill_noCuts" + datasetTagForPlots + ".png").c_str());

   // late time deuteron population

   TH1F* deuteron_Population = (TH1F*) inputFile->Get("topDir/Info/TriplesInfo/triple_deltaT13_vs_energy_lateTimes");

   auto canv_deuteron_Population = new TCanvas("canv_deuteron_Population","canv_deuteron_Population",700,0,600,450);

   deuteron_Population->Draw("colz");

   deuteron_Population->GetYaxis()->SetRangeUser(8, 18);
   canv_deuteron_Population->Update();

   canv_deuteron_Population->Write();
   canv_deuteron_Population->SaveAs(("Images/deuteron_Population" + datasetTagForPlots + ".png").c_str());

   // final cut time in fill plot

   TH2F* deltaT13_timeInFill_finalCuts = (TH2F*) inputFile->Get("topDir/AdditionalCuts/Triples/triple_deltaT13_vs_timeInFill");

   auto canv_deltaT13_timeInFill_finalCuts = new TCanvas("canv_deltaT13_timeInFill_finalCuts","canv_deltaT13_timeInFill_finalCuts",0,200,600,450);
   canv_deltaT13_timeInFill_finalCuts->SetLogz();

   nsTOus(deltaT13_timeInFill_finalCuts, "time [#mus]");
   deltaT13_timeInFill_finalCuts->Draw("colz");

   deltaT13_timeInFill_finalCuts->GetYaxis()->SetRangeUser(8, 18);
   canv_deltaT13_timeInFill_finalCuts->Update();

   canv_deltaT13_timeInFill_finalCuts->Write();
   canv_deltaT13_timeInFill_finalCuts->SaveAs(("Images/deltaT13_timeInFill_finalCuts" + datasetTagForPlots + ".png").c_str());

/////////////////////////////////////////////////////////////////////////////////////

    // deltaT12 final plots

  auto deltaT12_plots_dir = outputFile->mkdir("DeltaT12");
  deltaT12_plots_dir->cd();

   TH2F* deltaT12_timeInFill_baseCuts = (TH2F*) inputFile->Get("topDir/BaseCuts/Triples/triple_deltaT12_vs_timeInFill");

   auto canv_deltaT12_timeInFill_baseCuts = new TCanvas("canv_deltaT12_timeInFill_baseCuts","canv_deltaT12_timeInFill_baseCuts",100,200,600,450);
   canv_deltaT12_timeInFill_baseCuts->SetLogz();

   nsTOus(deltaT12_timeInFill_baseCuts, "time [#mus]");
   deltaT12_timeInFill_baseCuts->Draw("colz");

   deltaT12_timeInFill_baseCuts->GetYaxis()->SetRangeUser(2, 12);
   canv_deltaT12_timeInFill_baseCuts->Update();

   canv_deltaT12_timeInFill_baseCuts->Write();
   canv_deltaT12_timeInFill_baseCuts->SaveAs(("Images/deltaT12_timeInFill_baseCuts" + datasetTagForPlots + ".png").c_str());


   TH2F* deltaT12_timeInFill_finalCuts = (TH2F*) inputFile->Get("topDir/AdditionalCuts/Triples/triple_deltaT12_vs_timeInFill");

   auto canv_deltaT12_timeInFill_finalCuts = new TCanvas("canv_deltaT12_timeInFill_finalCuts","canv_deltaT12_timeInFill_finalCuts",200,200,600,450);
   canv_deltaT12_timeInFill_finalCuts->SetLogz();

   nsTOus(deltaT12_timeInFill_finalCuts, "time [#mus]");
   deltaT12_timeInFill_finalCuts->Draw("colz");

   deltaT12_timeInFill_finalCuts->GetYaxis()->SetRangeUser(2, 12);
   canv_deltaT12_timeInFill_finalCuts->Update();

   canv_deltaT12_timeInFill_finalCuts->Write();
   canv_deltaT12_timeInFill_finalCuts->SaveAs(("Images/deltaT12_timeInFill_finalCuts" + datasetTagForPlots + ".png").c_str());

/////////////////////////////////////////////////////////////////////////////////////

   // 1D plots

  gStyle->SetPadRightMargin(.08);

  auto fundamental_plots_dir = outputFile->mkdir("Fundamental");
  fundamental_plots_dir->cd();  

   TH1F* losses_integral = (TH1F*) inputFile->Get("topDir/AdditionalCuts/Triples/Losses/triple_losses_spectra_integral");
   losses_integral->GetYaxis()->SetTitle("Integral of triples");

   auto canv_losses_integral = new TCanvas("canv_losses_integral","canv_losses_integral",300,200,600,450);
   
   nsTOus(losses_integral, "time [#mus]");
   losses_integral->Draw("colz");

   canv_losses_integral->Write();
   canv_losses_integral->SaveAs(("Images/losses_integral" + datasetTagForPlots + ".png").c_str());

   // fractional loss plot

   TH1F* fractionalLossHist = (TH1F*) ((TH1F*) inputFile->Get("topDir/AdditionalCuts/Triples/Losses/triple_losses_spectra_integral"))->Clone("FractionalLosses");
   fractionalLossHist->GetYaxis()->SetTitle("Fractional losses");

   // double input_k_loss = 8.974 * 1e-10; // 60h
   // double input_k_loss = 2.510 * 1e-10; // 9d
   // double input_k_loss = 2.345 * 1e-10; // Endgame
   // double input_k_loss = 4.248 * 1e-10; // Endgame - negative side
   // double input_k_loss = 5.651 * 1e-10; // HighKick
   for (int bin = 1; bin <= fractionalLossHist->GetNbinsX(); ++bin) fractionalLossHist->SetBinContent(bin, fractionalLossHist->GetBinContent(bin) * input_k_loss);

   auto canv_fractional_losses = new TCanvas("canv_fractional_losses","canv_fractional_losses",400,200,600,450);
   fractionalLossHist->Draw("hist");

   canv_fractional_losses->Write();
   canv_fractional_losses->SaveAs(("Images/losses_integral_fraction" + datasetTagForPlots + ".png").c_str());

   // delta T and energy distribution plots

   auto canv_deltaT = new TCanvas("canv_deltaT","canv_deltaT",500,200,600,450);
   ((TH1F*) inputFile->Get("topDir/Info/TriplesInfo/triple_deltaT12"))->Draw("hist");

   canv_deltaT->Write();
   canv_deltaT->SaveAs(("Images/triples_deltaT" + datasetTagForPlots + ".png").c_str());

   auto canv_Edep = new TCanvas("canv_Edep","canv_Edep",600,200,600,450);
   ((TH1F*) inputFile->Get("topDir/Info/TriplesInfo/triple_Edep1"))->Draw("hist");

   canv_Edep->Write();
   canv_Edep->SaveAs(("Images/triples_Edep" + datasetTagForPlots + ".png").c_str());

/////////////////////////////////////////////////////////////////////////////////////   

  outputFile->Write();
  // delete outputFile; // deleting the output file will delete some histograms drawn to the screen


return 1;

}
