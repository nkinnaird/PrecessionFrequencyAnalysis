// 3-31-20: Old macro for plots related to the pileup construction. I don't remember exactly what it plots, but it would probably need to be dusted off to be reused.

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
#include <TVectorD.h>
#include <TFitResult.h>
#include <TSpectrum.h>
// #include <TRatioPlot.h>


#include "ratioAnalysisDefs.hh"
#include "pileupUtils.hh"


using namespace std;

// only energies looked at in this plotter so far
int pileupConstructionPlotter(std::string filePath)
{
  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  // pull in input file
  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  TFile* outputFile = new TFile("pileupConstruction.root","RECREATE");
  auto topDir = gFile->mkdir("topDir");

  int totalHistogramIters = (*(TVectorD*) inputFile->Get("Iters"))[0]; // total iterations in generated histograms (energy thresholds, etc.)
  int totalPasses = totalHistogramIters;

/////////////////////////////////////////////////////////////////////////////////////

for (int iter = 0; iter < totalPasses; ++iter)
{
  // make directory for specific fit conditions
  string dirString = Form("Iter%d", iter);
  auto iterDir = topDir->mkdir(dirString.c_str());
  iterDir->cd();

/////////////////////////////////////////////////////////////////////////////////////

  auto addedDir = iterDir->mkdir("addedDir");
  addedDir->cd();

    TH1F* allEnergies = (TH1F*) ((TH1F*) inputFile->Get(("topDir/" + dirString + "/Added/Energy").c_str()))->Clone();
    TH1F* allEnergiesPileupSubtracted = new TH1F( "Energy_Pileup_Subtracted", "Energy_Pileup_Subtracted; Energy; Events", 1000, 0, 10000);

    TH1F* allEnergies_early = (TH1F*) ((TH1F*) inputFile->Get(("topDir/" + dirString + "/Added/Energy_early").c_str()))->Clone();
    TH1F* allEnergies_late = (TH1F*) ((TH1F*) inputFile->Get(("topDir/" + dirString + "/Added/Energy_late").c_str()))->Clone();


    auto g2phaseDir = addedDir->mkdir("g2phase");
    g2phaseDir->cd();
    TH1F* allEnergies_g2phase[8];
      for (int g2phase = 0; g2phase < 8; ++g2phase) allEnergies_g2phase[g2phase] = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%i/Added/g2/Energies_g2phase%i", iter, g2phase)))->Clone();

    auto timeSliceDir = addedDir->mkdir("TimeSlice");
    timeSliceDir->cd();
    TH1F* allEnergies_timeSlice[10]; // this number needs to be automated somehow
      for (int timeSlice = 0; timeSlice < 10; ++timeSlice) allEnergies_timeSlice[timeSlice] = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%i/Added/TimeSlice/Energies_tS%i", iter, timeSlice)))->Clone();

/////////////////////////////////////////////////////////////////////////////////////

        auto added_pileupDir = addedDir->mkdir("Pileup");
        added_pileupDir->cd();

        TH1F* added_pileupEnergies_shadow = new TH1F("added_pileupEnergies_shadow", "added_pileupEnergies_shadow", 1000, 0, 10000);
        TH1F* added_pileupEnergies_early_shadow = new TH1F("added_pileupEnergies_early_shadow", "added_pileupEnergies_early_shadow", 1000, 0, 10000);
        TH1F* added_pileupEnergies_late_shadow = new TH1F("added_pileupEnergies_late_shadow", "added_pileupEnergies_late_shadow", 1000, 0, 10000);

          auto added_pileup_g2phaseDir = added_pileupDir->mkdir("g2phase");
          added_pileup_g2phaseDir->cd();
            TH1F* added_pileupEnergies_g2phase[8];
            for (int g2phase = 0; g2phase < 8; ++g2phase) added_pileupEnergies_g2phase[g2phase] = new TH1F(Form("added_pileupEnergies_g2phase%i",g2phase), Form("added_pileupEnergies_g2phase%i",g2phase), 1000, 0, 10000);

          auto added_pileup_timeSliceDir = added_pileupDir->mkdir("TimeSlice");
          added_pileup_timeSliceDir->cd();
            TH1F* added_pileupEnergies_timeSlice[10];
            for (int timeSlice = 0; timeSlice < 10; ++timeSlice) added_pileupEnergies_timeSlice[timeSlice] = new TH1F(Form("added_pileupEnergies_tS%i",timeSlice), Form("added_pileupEnergies_tS%i",timeSlice), 1000, 0, 10000);


/////////////////////////////////////////////////////////////////////////////////////

  // make calos directory where input data will live along with some analysis plots
  auto calosDirectory = iterDir->mkdir("Calos");
  calosDirectory->cd();

    // create individual calorimeter directories and histograms
    TDirectory* caloDirs[nCalos];

    TH1F* caloEnergies[nCalos];
    TH1F* caloEnergiesPileupSubtracted[nCalos];

    TH1F* caloEnergies_early[nCalos];
    TH1F* caloEnergies_late[nCalos];

    TH1F* caloEnergies_g2phase[nCalos][8];
    TH1F* caloEnergies_timeSlice[nCalos][10];

    TH1F* calo_pileupEnergies[nCalos];
    TH1F* calo_pileupEnergies_early[nCalos];
    TH1F* calo_pileupEnergies_late[nCalos];

    TH1F* calo_pileupEnergies_g2phase[nCalos][8];
    TH1F* calo_pileupEnergies_timeSlice[nCalos][10];

    for (int caloNum = 1; caloNum <= nCalos; ++caloNum)
    {
      caloDirs[caloNum-1] = calosDirectory->mkdir(Form("Calo%d",caloNum));
      caloDirs[caloNum-1]->cd();

      caloEnergies[caloNum-1] = (TH1F*) ((TH1F*) inputFile->Get(("topDir/" + dirString + "/Calos/Calo" + to_string(caloNum) + "/Calo" + to_string(caloNum) + "_Energy").c_str()))->Clone();
        caloEnergiesPileupSubtracted[caloNum-1] = new TH1F( Form("Calo%d_Energy_Pileup_Subtracted", caloNum), "; Energy; Events", 1000, 0, 10000);

      caloEnergies_early[caloNum-1] = (TH1F*) ((TH1F*) inputFile->Get(("topDir/" + dirString + "/Calos/Calo" + to_string(caloNum) + "/Calo" + to_string(caloNum) + "_Energy_early").c_str()))->Clone();
      caloEnergies_late[caloNum-1] = (TH1F*) ((TH1F*) inputFile->Get(("topDir/" + dirString + "/Calos/Calo" + to_string(caloNum) + "/Calo" + to_string(caloNum) + "_Energy_late").c_str()))->Clone();

          auto calo_g2phaseDir = caloDirs[caloNum-1]->mkdir("g2phase");
          calo_g2phaseDir->cd();
          for (int g2phase = 0; g2phase < 8; ++g2phase) caloEnergies_g2phase[caloNum-1][g2phase] = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%i/Calos/Calo%i/g2/Energies_g2phase%i", iter, caloNum, g2phase)))->Clone();

          auto calo_timeSliceDir = caloDirs[caloNum-1]->mkdir("TimeSlice");
          calo_timeSliceDir->cd();
          for (int timeSlice = 0; timeSlice < 10; ++timeSlice) caloEnergies_timeSlice[caloNum-1][timeSlice] = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%i/Calos/Calo%i/TimeSlice/Energies_tS%i", iter, caloNum, timeSlice)))->Clone();

      // copy pileup histograms over

        auto pileupDir = caloDirs[caloNum-1]->mkdir("Pileup");
        pileupDir->cd();

        calo_pileupEnergies[caloNum-1] = (TH1F*) ((TH1F*) inputFile->Get(("topDir/" + dirString + "/PileupHists/Calos/Calo" + to_string(caloNum) + "/pileupEnergies_shadow_Calo" + to_string(caloNum)).c_str()))->Clone();
        calo_pileupEnergies_early[caloNum-1] = (TH1F*) ((TH1F*) inputFile->Get(("topDir/" + dirString + "/PileupHists/Calos/Calo" + to_string(caloNum) + "/pileupEnergies_shadow_early_Calo" + to_string(caloNum)).c_str()))->Clone();
        calo_pileupEnergies_late[caloNum-1] = (TH1F*) ((TH1F*) inputFile->Get(("topDir/" + dirString + "/PileupHists/Calos/Calo" + to_string(caloNum) + "/pileupEnergies_shadow_late_Calo" + to_string(caloNum)).c_str()))->Clone();

            auto scaleFactDir = pileupDir->mkdir("ScaleFactor");
            scaleFactDir->cd();

            calcPileupScaleFactor(caloEnergies[caloNum-1], calo_pileupEnergies[caloNum-1]);

          auto calo_pileup_g2phaseDir = pileupDir->mkdir("g2phase");
          for (int g2phase = 0; g2phase < 8; ++g2phase)
          {
            calo_pileup_g2phaseDir->cd();
            calo_pileupEnergies_g2phase[caloNum-1][g2phase] = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%i/PileupHists/Calos/Calo%i/g2/pileupEnergies_shadow_g2phase%i", iter, caloNum, g2phase)))->Clone();
            auto g2phase_scaleFact_dir = calo_pileup_g2phaseDir->mkdir(Form("scaleFact%i",g2phase));
            g2phase_scaleFact_dir->cd();
            calcPileupScaleFactor(caloEnergies_g2phase[caloNum-1][g2phase], calo_pileupEnergies_g2phase[caloNum-1][g2phase]);
          }

          auto calo_pileup_timeSliceDir = pileupDir->mkdir("TimeSlice");
          for (int timeSlice = 0; timeSlice < 10; ++timeSlice)
          {
            calo_pileup_timeSliceDir->cd();
            calo_pileupEnergies_timeSlice[caloNum-1][timeSlice] = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%i/PileupHists/Calos/Calo%i/TimeSlice/pileupEnergies_shadow_tS%i", iter, caloNum, timeSlice)))->Clone();
            auto timeSlice_scaleFact_dir = calo_pileup_timeSliceDir->mkdir(Form("scaleFact%i",timeSlice));
            timeSlice_scaleFact_dir->cd();
            calcPileupScaleFactor(caloEnergies_timeSlice[caloNum-1][timeSlice], calo_pileupEnergies_timeSlice[caloNum-1][timeSlice]);
          }


    } // end calo loop

/////////////////////////////////////////////////////////////////////////////////////    

    // calculate best pileup scale factor by dividing energy histograms - do this here so that I can calculate a pileup scale factor before scaling and adding pileup hists

      auto autoScaleFactor_dir = added_pileupDir->mkdir("AutoScaleFactor");
      autoScaleFactor_dir->cd();

      TH1F* temp_added_calo_pileup_energies = new TH1F("temp_added_calo_pileup_energies", "temp_added_calo_pileup_energies", 1000, 0, 10000);
        for (int caloNum = 1; caloNum <= nCalos; ++caloNum) temp_added_calo_pileup_energies->Add(calo_pileupEnergies[caloNum-1]);
      double auto_pileup_scaleFactor = calcPileupScaleFactor(allEnergies, temp_added_calo_pileup_energies);
      delete temp_added_calo_pileup_energies; // so that this isn't saved into the file

/////////////////////////////////////////////////////////////////////////////////////

      // now apply pileup scale factor (calculated or explicit) and sum/subtract pileup spectra

      for (int caloNum = 1; caloNum <= nCalos; ++caloNum)
      {
          // scale pileup - comment this in and out for the time being - some stuff down below for g-2 phase and time slice

              // calo_pileupEnergies[caloNum-1]->Scale(auto_pileup_scaleFactor); // scaling the pileup results in a zeta that goes negative
              // calo_pileupEnergies_early[caloNum-1]->Scale(auto_pileup_scaleFactor);
              // calo_pileupEnergies_late[caloNum-1]->Scale(auto_pileup_scaleFactor);
 
          // create pileup subtracted histograms for each calo
              
              caloEnergiesPileupSubtracted[caloNum-1]->Add(caloEnergies[caloNum-1]);
              caloEnergiesPileupSubtracted[caloNum-1]->Add(calo_pileupEnergies[caloNum-1], -1);

          // construct added pileup spectra to later be subtracted off the added histograms

              added_pileupEnergies_shadow->Add(calo_pileupEnergies[caloNum-1]);
              added_pileupEnergies_early_shadow->Add(calo_pileupEnergies_early[caloNum-1]);
              added_pileupEnergies_late_shadow->Add(calo_pileupEnergies_late[caloNum-1]);

/////////////////////////////////////////////////////////////////////////////////////

              // add time slice and g2phase hists - careful with scaling stuff here

            for (int g2phase = 0; g2phase < 8; ++g2phase){
              // calo_pileupEnergies_g2phase[caloNum-1][g2phase]->Scale(auto_pileup_scaleFactor);
              added_pileupEnergies_g2phase[g2phase]->Add(calo_pileupEnergies_g2phase[caloNum-1][g2phase]);
            } 

            for (int timeSlice = 0; timeSlice < 10; ++timeSlice){
              // calo_pileupEnergies_timeSlice[caloNum-1][timeSlice]->Scale(auto_pileup_scaleFactor);
              added_pileupEnergies_timeSlice[timeSlice]->Add(calo_pileupEnergies_timeSlice[caloNum-1][timeSlice]);
            } 


      } // end calo loop


          // subtract pileup for added histograms

            allEnergiesPileupSubtracted->Add(allEnergies);
            allEnergiesPileupSubtracted->Add(added_pileupEnergies_shadow, -1);



/////////////////////////////////////////////////////////////////////////////////////
        
          for (int g2phase = 0; g2phase < 8; ++g2phase)
          {
            auto g2phase_scaleFact_dir = added_pileup_g2phaseDir->mkdir(Form("scaleFact%i",g2phase));
            g2phase_scaleFact_dir->cd();
            calcPileupScaleFactor(allEnergies_g2phase[g2phase], added_pileupEnergies_g2phase[g2phase]);
          }

          for (int timeSlice = 0; timeSlice < 10; ++timeSlice)
          {
            auto timeSlice_scaleFact_dir = added_pileup_timeSliceDir->mkdir(Form("scaleFact%i",timeSlice));
            timeSlice_scaleFact_dir->cd();
            calcPileupScaleFactor(allEnergies_timeSlice[timeSlice], added_pileupEnergies_timeSlice[timeSlice]);
          }


/////////////////////////////////////////////////////////////////////////////////////

          // calculate zeta efficiencies

          calcZetaEfficiency(addedDir, allEnergies_early, allEnergies_late, added_pileupEnergies_early_shadow, added_pileupEnergies_late_shadow);
          for (int caloNum = 1; caloNum <= nCalos; ++caloNum) calcZetaEfficiency(caloDirs[caloNum-1], caloEnergies_early[caloNum-1], caloEnergies_late[caloNum-1], calo_pileupEnergies_early[caloNum-1], calo_pileupEnergies_late[caloNum-1]);
            

} // end overall fit pass loops

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


      outputFile->Write();
      delete outputFile;

  return 1;
}
