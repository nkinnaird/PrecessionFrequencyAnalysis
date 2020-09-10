// 3-31-20: Macro for plots for pileup corrected vs uncorrected energy spectra.

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
#include "pileupUtils.hh"
#include "plotUtils.hh"

using namespace std;


void zoomInOnTail(TCanvas* canv)
{
  canv->SetLogy(0);

  THStack* stack = ((THStack*) canv->GetPrimitive("histStack"));
  stack->SetMinimum(-1000);
  stack->SetMaximum(2000);
  stack->SetTitle((string(stack->GetTitle()) + " Zoomed").c_str());

  canv->Update();
}


TCanvas* stackHists(TH1F* uncorrected, TH1F* corrected, string title)
{
  auto tempCanvas = new TCanvas("tempCanvas", "tempCanvas", 200, 10, 1200, 800);

  corrected->SetLineColor(2);

  auto legend = new TLegend(0.6,0.6,0.9,0.8);
  legend->AddEntry(uncorrected, "Pre-corrected", "l");
  legend->AddEntry(corrected, "Corrected", "l");


  THStack* histStack = new THStack("histStack","Temp Title"); // use a THStack because histogram attributes will be preserved when opening the canvas in a TBrowser
  histStack->Add(uncorrected);
  histStack->Add(corrected);

  histStack->SetTitle(("Energy Spectra - " + title).c_str());
  histStack->Draw("nostack,hist");

  histStack->GetXaxis()->SetTitle("Energy [MeV]");
  histStack->GetYaxis()->SetTitle("Events");


  legend->SetBorderSize(0);
  legend->Draw();

  tempCanvas->SetLogy();
  tempCanvas->Modified();

  return tempCanvas;
}



int energyCorrectedPlots(std::string filePath){

  // gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  // gStyle->SetOptFit(2);
  // gStyle->SetOptFit(1111);
  gStyle->SetOptFit(0);
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

  string outputFileName = "energyCorrectedPlots";
  TNamed* tag = applyDatasetTag(inputFile, outputFileName);

  TFile* outputFile = new TFile((outputFileName + ".root").c_str(),"RECREATE");
  if(tag) tag->Write();

/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

  auto addedDir = outputFile->mkdir("addedDir");
  addedDir->cd();

  TH1F* addedEnergySpectrum = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/Input/Energy");
  TH1F* correctedEnergySpectrum = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/Input/Energy_Pileup_Subtracted");

  TH1F* addedPileupEnergies = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/Pileup/added_pileupEnergies_shadow");

/////////////////////////////////////////////////////////////////////////////////////

  auto basicEnergyCanv = new TCanvas("basicEnergyCanv", "basicEnergyCanv", 200, 10, 1200, 800);

  TH1F* energySpecClone = (TH1F*) addedEnergySpectrum->Clone();

  energySpecClone->Draw("HIST");
  energySpecClone->GetXaxis()->SetRangeUser(0, 3200);

  basicEnergyCanv->SaveAs(("Images/basicEnergyHist" + datasetTagForPlots + ".png").c_str());

/////////////////////////////////////////////////////////////////////////////////////

  auto energyCompCanv = new TCanvas("energyCompCanv", "energyCompCanv", 200, 10, 1200, 800);

  addedPileupEnergies->SetLineColor(2);
  for (int bin = 1; bin <= addedPileupEnergies->GetNbinsX(); ++bin) addedPileupEnergies->SetBinContent(bin, abs(addedPileupEnergies->GetBinContent(bin)));

  auto legend = new TLegend(0.54,0.6,0.84,0.8);
  legend->AddEntry(addedEnergySpectrum, "Cluster Energies", "l");
  legend->AddEntry(addedPileupEnergies, "Pileup Energies", "l");

  THStack* histStack = new THStack("histStack","Temp Title");
  histStack->Add(addedEnergySpectrum);
  histStack->Add(addedPileupEnergies);

  histStack->SetTitle("Cluster Energies vs Pileup Energies");
  histStack->Draw("nostack,hist");

  histStack->GetXaxis()->SetTitle("Energy [MeV]");
  histStack->GetYaxis()->SetTitle("Events");

  legend->SetBorderSize(0);
  legend->Draw();

  energyCompCanv->SetLogy();

  // histStack->GetXaxis()->SetRangeUser(0, 5500);
  // histStack->SetMinimum(100);

  energyCompCanv->Modified();
  energyCompCanv->Write("ClusterEnergiesVsPileupEnergies");
  energyCompCanv->SaveAs(("Images/ClusterEnergiesVsPileupEnergies" + datasetTagForPlots + ".png").c_str());


/////////////////////////////////////////////////////////////////////////////////////

  TCanvas* myCanvas;

  myCanvas = stackHists(addedEnergySpectrum, correctedEnergySpectrum, "Added Calos");
  myCanvas->Write("AddedEnergies");
  myCanvas->SaveAs(("Images/AddedEnergies" + datasetTagForPlots + ".png").c_str());
  zoomInOnTail(myCanvas);
  myCanvas->SaveAs(("Images/AddedEnergiesZoomed" + datasetTagForPlots + ".png").c_str());

/////////////////////////////////////////////////////////////////////////////////////

  auto splitCanv = new TCanvas("splitCanv", "splitCanv", 200, 10, 1200, 800);
  splitCanv->Divide(2,2);

  auto splitCanvZoomed = new TCanvas("splitCanvZoomed", "splitCanvZoomed", 200, 10, 1200, 800);
  splitCanvZoomed->Divide(2,2);

  TH1F* calo1_EnergySpectrum = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/Calos/Calo1/Input/Calo1_Energy");
  TH1F* calo1_CorrectedEnergySpectrum = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/Calos/Calo1/Input/Calo1_Energy_Pileup_Subtracted");

  TH1F* calo7_EnergySpectrum = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/Calos/Calo7/Input/Calo7_Energy");
  TH1F* calo7_CorrectedEnergySpectrum = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/Calos/Calo7/Input/Calo7_Energy_Pileup_Subtracted");

  TH1F* calo13_EnergySpectrum = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/Calos/Calo13/Input/Calo13_Energy");
  TH1F* calo13_CorrectedEnergySpectrum = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/Calos/Calo13/Input/Calo13_Energy_Pileup_Subtracted");

  TH1F* calo23_EnergySpectrum = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/Calos/Calo23/Input/Calo23_Energy");
  TH1F* calo23_CorrectedEnergySpectrum = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/Calos/Calo23/Input/Calo23_Energy_Pileup_Subtracted");


  myCanvas = stackHists(calo1_EnergySpectrum, calo1_CorrectedEnergySpectrum, "Calo 1");
  splitCanv->cd(1);
  myCanvas->DrawClonePad();
  splitCanvZoomed->cd(1);
  zoomInOnTail(myCanvas);
  myCanvas->DrawClonePad();

  myCanvas = stackHists(calo7_EnergySpectrum, calo7_CorrectedEnergySpectrum, "Calo 7");
  splitCanv->cd(2);
  myCanvas->DrawClonePad();
  splitCanvZoomed->cd(2);
  zoomInOnTail(myCanvas);
  myCanvas->DrawClonePad();

  myCanvas = stackHists(calo13_EnergySpectrum, calo13_CorrectedEnergySpectrum, "Calo 13");
  splitCanv->cd(3);
  myCanvas->DrawClonePad();
  splitCanvZoomed->cd(3);
  zoomInOnTail(myCanvas);
  myCanvas->DrawClonePad();

  myCanvas = stackHists(calo23_EnergySpectrum, calo23_CorrectedEnergySpectrum, "Calo 23");
  splitCanv->cd(4);
  myCanvas->DrawClonePad();
  splitCanvZoomed->cd(4);
  zoomInOnTail(myCanvas);
  myCanvas->DrawClonePad();

  splitCanv->Update();
  splitCanv->Write("CaloEnergies");
  splitCanv->SaveAs(("Images/CaloEnergies" + datasetTagForPlots + ".png").c_str());

  splitCanvZoomed->Update();
  splitCanvZoomed->SaveAs(("Images/CaloEnergiesZoomed" + datasetTagForPlots + ".png").c_str());

/////////////////////////////////////////////////////////////////////////////////////

// plot pileup time spectra

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(2);
  gStyle->SetPadRightMargin(.05);


  // TH1F* addedTimeSpectrum = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/Input/Times_E_Threshold_Pileup_Subtracted");
  // nsTOus(addedTimeSpectrum, "Time (#mus)");


  auto timeCanvas = new TCanvas("timeCanvas", "timeCanvas", 200, 10, 1200, 600);
  timeCanvas->Divide(2,1);

  double xMaxBound = 300; // us

  TH1F* pileupTimeSpectrum = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/Pileup/added_pileupTimes_shadow_threshold");
  nsTOus(pileupTimeSpectrum, "Time (#mus)");
  pileupTimeSpectrum->GetXaxis()->SetRangeUser(0, xMaxBound);
  pileupTimeSpectrum->SetTitle("Pileup Time Spectrum (Above Threshold)");
  pileupTimeSpectrum->GetYaxis()->SetTitle("Entries");
  pileupTimeSpectrum->GetYaxis()->SetTitleOffset(1.8);

  timeCanvas->cd(1);
  pileupTimeSpectrum->Draw("hist");

  timeCanvas->cd(2);
  gPad->SetLogy();

  timeCanvas->Update();

  TF1* expFunc = new TF1("expFunc", "[0]*exp(-x/[1])", 25, xMaxBound);
  expFunc->SetParameter(0, 1e4); 
  expFunc->SetParameter(1, defaultLifetime/2.);
  expFunc->SetLineColor(2);

  pileupTimeSpectrum->Fit(expFunc, "QRL");

  pileupTimeSpectrum->Draw("hist");
  expFunc->Draw("SAME");

  timeCanvas->Write("PileupTimeSpectrum");
  timeCanvas->SaveAs(("Images/PileupTimeSpectrum" + datasetTagForPlots + ".png").c_str());

/////////////////////////////////////////////////////////////////////////////////////

    outputFile->Write();
    delete outputFile;

    return 1;

}
