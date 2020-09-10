// 3-31-20: Macro for plots for calorimeter time and energy spectra on top of each other. Hasn't been used in a long time and would need to be dusted off.

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


#include "../ratioMacroHeaders/ratioAnalysisDefs.hh"
#include "../ratioMacroHeaders/pileupUtils.hh"
#include "../ratioMacroHeaders/plotUtils.hh"

using namespace std;


void stackCaloHists(TH1F** hists, string id, string xLabel)
{

  THStack* histStack = new THStack("histStack","Temp Title"); // use a THStack because histogram attributes will be preserved when opening the canvas in a TBrowser
  
  string stackTitle = id;
  histStack->SetTitle(stackTitle.c_str());

  auto legend = new TLegend(0.8,0.2,0.9,0.9);

  for (uint i = 0; i < nCalos; ++i)
  {
    hists[i]->SetLineColor(53+2*i);
    histStack->Add(hists[i]);

    legend->AddEntry(hists[i],Form("Calo %d", i+1),"l");
  }


  string canvasName = "stacked_";
  canvasName.append(id);

  auto stackedCanvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 200, 10, 1200, 800);
  histStack->Draw("nostack,hist");

  histStack->GetXaxis()->SetTitle(xLabel.c_str());
  histStack->GetYaxis()->SetTitle("Events");

  legend->SetBorderSize(1);
  legend->Draw();


  stackedCanvas->Write();
  delete stackedCanvas;
}



int caloPileupCompare(std::string filePath){

  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  TFile* outputFile = new TFile("caloPileupCompare.root","RECREATE");
  auto topDir = outputFile->mkdir("topDir");

  auto calosDirectory = topDir->mkdir("Calos");
  calosDirectory->cd();

/////////////////////////////////////////////////////////////////////////////////////

  TH1F* scaleFactors = new TH1F("PileupScaleFactors", "PileupScaleFactors",100, .95, 1.05);

    TDirectory* caloDirs[nCalos];

    TH1F* caloEnergies[nCalos];
      TH1F* caloEnergiesPileupSubtracted[nCalos];

    TH1F* caloTimes[nCalos];
      TH1F* caloTimesPileupSubtracted[nCalos];
    TH1F* caloTimesThreshold[nCalos];
      TH1F* caloTimesThresholdPileupSubtracted[nCalos];


      TH1F* calo_pileupEnergies[nCalos];
      TH1F* calo_pileupTimes[nCalos];
      TH1F* calo_pileupTimes_threshold[nCalos];

      TH1F* calo_deltaTs[nCalos];

/////////////////////////////////////////////////////////////////////////////////////

      string dirString = "Iter0"; // hard set for now

/////////////////////////////////////////////////////////////////////////////////////

    // set bin info
    double binWidth = ((TH1F*) inputFile->Get(("topDir/" + dirString + "/Calos/Calo" + to_string(caloNum) + "/Calo" + to_string(caloNum) + "_Times").c_str()))->GetBinWidth(1);
    int nBins = ((TH1F*) inputFile->Get(("topDir/" + dirString + "/Calos/Calo" + to_string(caloNum) + "/Calo" + to_string(caloNum) + "_Times").c_str()))->GetNbinsX();
    double histMinTime = ((TH1F*) inputFile->Get(("topDir/" + dirString + "/Calos/Calo" + to_string(caloNum) + "/Calo" + to_string(caloNum) + "_Times").c_str()))->GetBinLowEdge(1);
    double histMaxTime = nBins*binWidth + histMinTime;

    for (int caloNum = 1; caloNum <= nCalos; ++caloNum)
    {
      caloDirs[caloNum-1] = calosDirectory->mkdir(Form("Calo%d",caloNum));
      caloDirs[caloNum-1]->cd();

      auto inputDir = caloDirs[caloNum-1]->mkdir("Input");
      inputDir->cd();

      caloEnergies[caloNum-1] = (TH1F*) ((TH1F*) inputFile->Get(("topDir/" + dirString + "/Calos/Calo" + to_string(caloNum) + "/Calo" + to_string(caloNum) + "_Energy").c_str()))->Clone();
        caloEnergiesPileupSubtracted[caloNum-1] = new TH1F( Form("Calo%d_Energy_Pileup_Subtracted", caloNum), "; Energy; Events", 1000, 0, 10000);

      caloTimes[caloNum-1] = (TH1F*) ((TH1F*) inputFile->Get(("topDir/" + dirString + "/Calos/Calo" + to_string(caloNum) + "/Calo" + to_string(caloNum) + "_Times").c_str()))->Clone();
        caloTimesPileupSubtracted[caloNum-1] = new TH1F(Form("Calo%d_Times_Pileup_Subtracted", caloNum), "; time (ns); Events", nBins, histMinTime, histMaxTime);
      caloTimesThreshold[caloNum-1] = (TH1F*) ((TH1F*) inputFile->Get(("topDir/" + dirString + "/Calos/Calo" + to_string(caloNum) + "/Calo" + to_string(caloNum) + "_Times_E_Threshold").c_str()))->Clone();
        caloTimesThresholdPileupSubtracted[caloNum-1] = new TH1F(Form("Calo%d_Times_E_Threshold_Pileup_Subtracted", caloNum), "; time (ns); Events", nBins, histMinTime, histMaxTime);

      nsTOus(caloTimes[caloNum-1], "Time (#mus)");
      nsTOus(caloTimesPileupSubtracted[caloNum-1], "Time (#mus)");
      nsTOus(caloTimesThreshold[caloNum-1], "Time (#mus)");
      nsTOus(caloTimesThresholdPileupSubtracted[caloNum-1], "Time (#mus)");

      // copy pileup histograms over

        auto pileupDir = caloDirs[caloNum-1]->mkdir("Pileup");
        pileupDir->cd();

          calo_deltaTs[caloNum-1] = (TH1F*) ((TH1F*) inputFile->Get(("topDir/" + dirString + "/PileupHists/Calos/Calo" + to_string(caloNum) + "/Time_Between_Hits_Calo" + to_string(caloNum)).c_str()))->Clone();

          calo_pileupEnergies[caloNum-1] = (TH1F*) ((TH1F*) inputFile->Get(("topDir/" + dirString + "/PileupHists/Calos/Calo" + to_string(caloNum) + "/pileupEnergies_shadow_Calo" + to_string(caloNum)).c_str()))->Clone();
          calo_pileupTimes[caloNum-1] = (TH1F*) ((TH1F*) inputFile->Get(("topDir/" + dirString + "/PileupHists/Calos/Calo" + to_string(caloNum) + "/pileupTimes_shadow_Calo" + to_string(caloNum)).c_str()))->Clone();
          calo_pileupTimes_threshold[caloNum-1] = (TH1F*) ((TH1F*) inputFile->Get(("topDir/" + dirString + "/PileupHists/Calos/Calo" + to_string(caloNum) + "/pileupTimes_shadow_threshold_Calo" + to_string(caloNum)).c_str()))->Clone();

          nsTOus(calo_pileupTimes[caloNum-1], "Time (#mus)");
          nsTOus(calo_pileupTimes_threshold[caloNum-1], "Time (#mus)");

/////////////////////////////////////////////////////////////////////////////////////

      auto autoScaleFactor_dir = pileupDir->mkdir("AutoScaleFactor");
      autoScaleFactor_dir->cd();

      double auto_pileup_scaleFactor = calcPileupScaleFactor(caloEnergies[caloNum-1], calo_pileupEnergies[caloNum-1]);
      scaleFactors->Fill(auto_pileup_scaleFactor);

      cout << "Calo: " << caloNum << " scale factor: " << auto_pileup_scaleFactor << endl;

    } // end calo loop

/////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////

      for (int caloNum = 1; caloNum <= nCalos; ++caloNum)
      {
          // create pileup subtracted histograms for each calo
              
              caloEnergiesPileupSubtracted[caloNum-1]->Add(caloEnergies[caloNum-1]);
              caloEnergiesPileupSubtracted[caloNum-1]->Add(calo_pileupEnergies[caloNum-1], -1);

              caloTimesPileupSubtracted[caloNum-1]->Add(caloTimes[caloNum-1]);
              caloTimesPileupSubtracted[caloNum-1]->Add(calo_pileupTimes[caloNum-1], -1);
              caloTimesThresholdPileupSubtracted[caloNum-1]->Add(caloTimesThreshold[caloNum-1]);
              caloTimesThresholdPileupSubtracted[caloNum-1]->Add(calo_pileupTimes_threshold[caloNum-1], -1);

      } // end calo loop

/////////////////////////////////////////////////////////////////////////////////////

  auto canvasDir = topDir->mkdir("Canvases");
  auto unnormalizedDir = canvasDir->mkdir("Unnormalized");
  auto normalizedDir = canvasDir->mkdir("Normalized");

  unnormalizedDir->cd();

  stackCaloHists(calo_deltaTs, "Delta_Ts", "#DeltaT (ns)");

  stackCaloHists(caloEnergies, "Energy", "Energy (MeV)");
  stackCaloHists(caloEnergiesPileupSubtracted, "Energy_Pileup_Subtracted", "Energy (MeV)");
  stackCaloHists(caloTimesPileupSubtracted, "Times_Pileup_Subtracted", "Time (#mus)");
  stackCaloHists(caloTimesThresholdPileupSubtracted, "Times_Threshold_Pileup_Subtracted", "Time (#mus)");

  stackCaloHists(calo_pileupEnergies, "Pileup_Energies", "Energy (MeV)");
  stackCaloHists(calo_pileupTimes, "Pileup_Times", "Time (#mus)");
  stackCaloHists(calo_pileupTimes_threshold, "Pileup_Times_Threshold", "Time (#mus)");

/////////////////////////////////////////////////////////////////////////////////////

  normalizedDir->cd();

  double normFactors[nCalos];
  for (int i = 0; i < 24; ++i)
  {
    double xmin = 500;
    double xmax = 7000;
    int binmin = caloEnergies[i]->GetXaxis()->FindBin(xmin);
    int binmax = caloEnergies[i]->GetXaxis()->FindBin(xmax);

    normFactors[i] = caloEnergies[i]->Integral(binmin, binmax); // should I be normalizing to this histogram or something else? original was normalizing to itself down below

    caloEnergies[i]->Scale(abs(1./normFactors[i]));
    caloEnergiesPileupSubtracted[i]->Scale(abs(1./normFactors[i]));
    caloTimesPileupSubtracted[i]->Scale(abs(1./normFactors[i]));
    caloTimesThresholdPileupSubtracted[i]->Scale(abs(1./normFactors[i]));
    calo_pileupEnergies[i]->Scale(abs(1./normFactors[i]));
    calo_pileupTimes[i]->Scale(abs(1./normFactors[i]));
    calo_pileupTimes_threshold[i]->Scale(abs(1./normFactors[i]));

    // caloEnergiesPileupSubtracted[i]->Scale(abs(1./caloEnergiesPileupSubtracted[i]->Integral()));
    // caloTimesPileupSubtracted[i]->Scale(abs(1./caloTimesPileupSubtracted[i]->Integral()));
    // caloTimesThresholdPileupSubtracted[i]->Scale(abs(1./caloTimesThresholdPileupSubtracted[i]->Integral()));
    // calo_pileupEnergies[i]->Scale(abs(1./calo_pileupEnergies[i]->Integral()));
    // calo_pileupTimes[i]->Scale(abs(1./calo_pileupTimes[i]->Integral()));
    // calo_pileupTimes_threshold[i]->Scale(abs(1./calo_pileupTimes_threshold[i]->Integral()));
  }

  stackCaloHists(caloEnergies, "Energy_Normalized", "Energy (MeV)");
  stackCaloHists(caloEnergiesPileupSubtracted, "Energy_Pileup_Subtracted_Normalized", "Energy (MeV)");
  stackCaloHists(caloTimesPileupSubtracted, "Times_Pileup_Subtracted_Normalized", "Time (#mus)");
  stackCaloHists(caloTimesThresholdPileupSubtracted, "Times_Threshold_Pileup_Subtracted_Normalized", "Time (#mus)");

  stackCaloHists(calo_pileupEnergies, "Pileup_Energies_Normalized", "Energy (MeV)");
  stackCaloHists(calo_pileupTimes, "Pileup_Times_Normalized", "Time (#mus)");
  stackCaloHists(calo_pileupTimes_threshold, "Pileup_Times_Threshold_Normalized", "Time (#mus)");


  // normalize to above 500 MeV - will need to store arrays of normalization factors for pileup energies and times (one set of 24), and calo energies and times (another set of 24)
  // will have to change method for this to work


/////////////////////////////////////////////////////////////////////////////////////


    outputFile->Write();
    delete outputFile;

    return 1;

}
