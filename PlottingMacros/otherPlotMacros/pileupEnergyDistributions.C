// Macro to put pileup energy distributions in a root file for figure making.

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

int pileupEnergyDistributions(){

  TFile *file_60h = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/60h/SingleIter/SingleFit-NewRange/output-60h-SingleIter-NewRange.root");
  TFile *file_HK = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/HighKick/SingleIter/SingleFit-NewRange/output-HighKick-SingleIter-NewRange.root");
  TFile *file_9d = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/9d/SingleIter/SingleFit-NewRange/output-9d-SingleIter-NewRange.root");
  TFile *file_EG = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/Endgame/SingleIter/SingleFit-NewRange/output-Endgame-SingleIter-NewRange.root");


   if (file_60h == 0 || file_HK == 0 || file_9d == 0 || file_EG == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  TFile* outputFile = new TFile("pileupDistributions_Kinnaird.root","RECREATE");

/////////////////////////////////////////////////////////////////////////////////////

  auto dir_60h = outputFile->mkdir("60h");
  dir_60h->cd();

  TH1F* clusterEnergySpectrum_60h = (TH1F*) file_60h->Get("topDir/FitPasses/FitPass0/addedDir/Input/Energy");
  TH1F* pileupEnergySpectrum_60h = (TH1F*) file_60h->Get("topDir/FitPasses/FitPass0/addedDir/Pileup/added_pileupEnergies_shadow");
  TH1F* correctedEnergySpectrum_60h = (TH1F*) file_60h->Get("topDir/FitPasses/FitPass0/addedDir/Input/Energy_Pileup_Subtracted");

  clusterEnergySpectrum_60h->Write("ClusterEnergies");
  pileupEnergySpectrum_60h->Write("PileupEnergies");
  correctedEnergySpectrum_60h->Write("PileupCorrectedClusterEnergies");


  auto dir_HK = outputFile->mkdir("HK");
  dir_HK->cd();

  TH1F* clusterEnergySpectrum_HK = (TH1F*) file_HK->Get("topDir/FitPasses/FitPass0/addedDir/Input/Energy");
  TH1F* pileupEnergySpectrum_HK = (TH1F*) file_HK->Get("topDir/FitPasses/FitPass0/addedDir/Pileup/added_pileupEnergies_shadow");
  TH1F* correctedEnergySpectrum_HK = (TH1F*) file_HK->Get("topDir/FitPasses/FitPass0/addedDir/Input/Energy_Pileup_Subtracted");

  clusterEnergySpectrum_HK->Write("ClusterEnergies");
  pileupEnergySpectrum_HK->Write("PileupEnergies");
  correctedEnergySpectrum_HK->Write("PileupCorrectedClusterEnergies");


  auto dir_9d = outputFile->mkdir("9d");
  dir_9d->cd();

  TH1F* clusterEnergySpectrum_9d = (TH1F*) file_9d->Get("topDir/FitPasses/FitPass0/addedDir/Input/Energy");
  TH1F* pileupEnergySpectrum_9d = (TH1F*) file_9d->Get("topDir/FitPasses/FitPass0/addedDir/Pileup/added_pileupEnergies_shadow");
  TH1F* correctedEnergySpectrum_9d = (TH1F*) file_9d->Get("topDir/FitPasses/FitPass0/addedDir/Input/Energy_Pileup_Subtracted");

  clusterEnergySpectrum_9d->Write("ClusterEnergies");
  pileupEnergySpectrum_9d->Write("PileupEnergies");
  correctedEnergySpectrum_9d->Write("PileupCorrectedClusterEnergies");


  auto dir_EG = outputFile->mkdir("EG");
  dir_EG->cd();

  TH1F* clusterEnergySpectrum_EG = (TH1F*) file_EG->Get("topDir/FitPasses/FitPass0/addedDir/Input/Energy");
  TH1F* pileupEnergySpectrum_EG = (TH1F*) file_EG->Get("topDir/FitPasses/FitPass0/addedDir/Pileup/added_pileupEnergies_shadow");
  TH1F* correctedEnergySpectrum_EG = (TH1F*) file_EG->Get("topDir/FitPasses/FitPass0/addedDir/Input/Energy_Pileup_Subtracted");

  clusterEnergySpectrum_EG->Write("ClusterEnergies");
  pileupEnergySpectrum_EG->Write("PileupEnergies");
  correctedEnergySpectrum_EG->Write("PileupCorrectedClusterEnergies");

/////////////////////////////////////////////////////////////////////////////////////

    outputFile->Write();
    delete outputFile;

    return 1;

}
