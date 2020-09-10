// 3-31-20: Macro which produces lost muon histograms from passed in files, including various background subtractions. For an entire dataset this macro should be ran on every file and then hadded together.

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

#include "../ratioMacroHeaders/ratioAnalysisDefs.hh" // don't change these paths since this file is run on the grid (whereas other files can see my root include directory and the ../ can be removed)
#include "../ratioMacroHeaders/ratioLostMuonsConfig.hh" // don't change these paths since this file is run on the grid (whereas other files can see my root include directory and the ../ can be removed)

/////////////////////////////////////////////////////////////////////////////////////

// for storing config info into root file

const char* BoolToString(bool b)
{
  return b ? "true" : "false";
}

/////////////////////////////////////////////////////////////////////////////////////

int LostMuonHistograms(std::string filePath)
{

/////////////////////////////////////////////////////////////////////////////////////

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 1;
   }

/////////////////////////////////////////////////////////////////////////////////////

 // need to create unique output root file names for parallel processing - so grab the input file name which has a unique id

  size_t delPosition = filePath.find_last_of("/");
  string uniqueID = filePath.substr(delPosition+1);
  delPosition = uniqueID.find_first_of("_");
  uniqueID = uniqueID.substr(delPosition+1);
  delPosition = uniqueID.find_first_of("_");
  uniqueID = uniqueID.substr(delPosition+1);
  delPosition = uniqueID.find_last_of(".");
  uniqueID = uniqueID.substr(0, delPosition);

  string outputFileName = "lostMuons_hists_" + uniqueID + ".root";

  cout << "Output file name: " << outputFileName << endl;

  TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  auto topDir = gFile->mkdir("topDir");

/////////////////////////////////////////////////////////////////////////////////////

  int nBins = int(approxMaxTime/defaultBinWidth);
  double histMinTime = 0;
  double histMaxTime = nBins*defaultBinWidth;

  int timeBinsIn2D = nBins; // 700;
  double histMinTime2D = histMinTime; // 0;
  double histMaxTime2D = histMaxTime; // 700000;

/////////////////////////////////////////////////////////////////////////////////////

  // info dirs

  auto Info_dir = topDir->mkdir("Info");

  auto firstHitInfo_dir = Info_dir->mkdir("FirstHitInfo");
  firstHitInfo_dir->cd();

  TH1F* firstTimes = new TH1F("First_Times", "First Times; time (ns); Events", nBins, histMinTime, histMaxTime);
  TH1F* firstSizes = new TH1F("First_Cluster_Sizes", "First Cluster Sizes; Cluster size; Events", 10, 0, 10);
  TH1F* firstEnergyFracs = new TH1F("First_Energy_Fracs", "First Energy Fractions; Energy fraction; Events", 100, 0, 1);

  auto doublesInfo_dir = Info_dir->mkdir("DoublesInfo");
  doublesInfo_dir->cd();

  TH1F* double_deltaT12 = new TH1F("double_deltaT12", "#Deltat_{12}; #Deltat_{12} (ns); Events", 200, 0, 20);
  TH1F* double_Edep1 = new TH1F("double_Edep1", "E_{1}; Energy (MeV); Events", 500, 0, 500);
  TH1F* double_Edep2 = new TH1F("double_Edep2", "E_{2}; Energy (MeV); Events", 500, 0, 500);
  TH2F* double_calo_vs_deltaT12 = new TH2F("double_calo_vs_deltaT12", "Calo Num vs #Deltat_{12}; #Deltat_{12} (ns); Calo Num", 200, 0, 20, 24, 1, 25);
  
  TH2F* double_deltaT12_vs_energy = new TH2F("double_deltaT12_vs_energy", "Calo Num vs #Deltat_{12}; E_{1} (MeV); #Deltat_{12} (ns)", 100, 0, 500, 200, 0, 20);
  TH2F* double_deltaT12_vs_energy_lateTimes = new TH2F("double_deltaT12_vs_energy_lateTimes", "Calo Num vs #Deltat_{12} after 300 #mus; E_{1} (MeV); #Deltat_{12} (ns)", 100, 0, 500, 200, 0, 20);
  TH2F* double_deltaT12_vs_timeInFill = new TH2F("double_deltaT12_vs_timeInFill", "Calo Num vs #Deltat_{12}; time (ns); #Deltat_{12} (ns)", timeBinsIn2D, histMinTime2D, histMaxTime2D, 200, 0, 20);


  auto triplesInfo_dir = Info_dir->mkdir("TriplesInfo");
  triplesInfo_dir->cd();

  TH1F* triple_deltaT12 = new TH1F("triple_deltaT12", "#Deltat_{12}; #Deltat_{12} (ns); Events", 200, 0, 20);
  TH1F* triple_deltaT13 = new TH1F("triple_deltaT13", "#Deltat_{13}; #Deltat_{12} (ns); Events", 200, 0, 20);
  TH1F* triple_deltaT23 = new TH1F("triple_deltaT23", "#Deltat_{23}; #Deltat_{12} (ns); Events", 200, 0, 20);
  TH1F* triple_Edep1 = new TH1F("triple_Edep1", "E_{1}; Energy (MeV); Events", 500, 0, 500);
  TH1F* triple_Edep2 = new TH1F("triple_Edep2", "E_{2}; Energy (MeV); Events", 500, 0, 500);
  TH1F* triple_Edep3 = new TH1F("triple_Edep3", "E_{3}; Energy (MeV); Events", 500, 0, 500);

  TH2F* triple_deltaT12_vs_energy = new TH2F("triple_deltaT12_vs_energy", "Calo Num vs #Deltat_{12}; E_{1} (MeV); #Deltat_{12} (ns)", 100, 0, 500, 200, 0, 20);
  TH2F* triple_deltaT12_vs_energy_lateTimes = new TH2F("triple_deltaT12_vs_energy_lateTimes", "Calo Num vs #Deltat_{12} after 300 #mus; E_{1} (MeV); #Deltat_{12} (ns)", 100, 0, 500, 200, 0, 20);
  TH2F* triple_deltaT13_vs_energy = new TH2F("triple_deltaT13_vs_energy", "Calo Num vs #Deltat_{13}; E_{1} (MeV); #Deltat_{13} (ns)", 100, 0, 500, 200, 0, 20);
  TH2F* triple_deltaT13_vs_energy_lateTimes = new TH2F("triple_deltaT13_vs_energy_lateTimes", "Calo Num vs #Deltat_{13} after 300 #mus; E_{1} (MeV); #Deltat_{13} (ns)", 100, 0, 500, 200, 0, 20);
  TH2F* triple_deltaT23_vs_energy = new TH2F("triple_deltaT23_vs_energy", "Calo Num vs #Deltat_{23}; E_{1} (MeV); #Deltat_{23} (ns)", 100, 0, 500, 200, 0, 20);
  TH2F* triple_deltaT23_vs_energy_lateTimes = new TH2F("triple_deltaT23_vs_energy_lateTimes", "Calo Num vs #Deltat_{23} after 300 #mus; E_{1} (MeV); #Deltat_{23} (ns)", 100, 0, 500, 200, 0, 20);
  TH2F* triple_deltaT12_vs_timeInFill = new TH2F("triple_deltaT12_vs_timeInFill", "Calo Num vs #Deltat_{12}; time (ns); #Deltat_{12} (ns)", timeBinsIn2D, histMinTime2D, histMaxTime2D, 200, 0, 20);
  TH2F* triple_deltaT13_vs_timeInFill = new TH2F("triple_deltaT13_vs_timeInFill", "Calo Num vs #Deltat_{13}; time (ns); #Deltat_{13} (ns)", timeBinsIn2D, histMinTime2D, histMaxTime2D, 200, 0, 20);
  TH2F* triple_deltaT23_vs_timeInFill = new TH2F("triple_deltaT23_vs_timeInFill", "Calo Num vs #Deltat_{23}; time (ns); #Deltat_{23} (ns)", timeBinsIn2D, histMinTime2D, histMaxTime2D, 200, 0, 20);


  auto quadruplesInfo_dir = Info_dir->mkdir("QuadruplesInfo");
  quadruplesInfo_dir->cd();


/////////////////////////////////////////////////////////////////////////////////////

  // no cuts

  auto baseCuts_dir = topDir->mkdir("BaseCuts");

  auto listCuts_dir = baseCuts_dir->mkdir("Cuts"); // make this directory so that when hadding the many copies of the lostMuonCuts get put in their own directory
  listCuts_dir->cd();

    TVectorD lostMuonCuts(4); 
    lostMuonCuts[0] = deltaTCutLow;
    lostMuonCuts[1] = deltaTCutHigh;
    lostMuonCuts[2] = energyCutLow;
    lostMuonCuts[3] = energyCutHigh;
    lostMuonCuts.Write("LostMuonCuts");

  baseCuts_dir->cd();

  auto triples_baseCuts_dir = baseCuts_dir->mkdir("Triples");
  triples_baseCuts_dir->cd(); 

    TH2F* triple_deltaT12_vs_energy_baseCuts = new TH2F("triple_deltaT12_vs_energy", "Calo Num vs #Deltat_{12}; E_{1} (MeV); #Deltat_{12} (ns)", 100, 0, 500, 200, 0, 20);
    TH2F* triple_deltaT12_vs_energy_lateTimes_baseCuts = new TH2F("triple_deltaT12_vs_energy_lateTimes", "Calo Num vs #Deltat_{12} after 300 #mus; E_{1} (MeV); #Deltat_{12} (ns)", 100, 0, 500, 200, 0, 20);
    TH2F* triple_deltaT13_vs_energy_baseCuts = new TH2F("triple_deltaT13_vs_energy", "Calo Num vs #Deltat_{13}; E_{1} (MeV); #Deltat_{13} (ns)", 100, 0, 500, 200, 0, 20);
    TH2F* triple_deltaT13_vs_energy_lateTimes_baseCuts = new TH2F("triple_deltaT13_vs_energy_lateTimes", "Calo Num vs #Deltat_{13} after 300 #mus; E_{1} (MeV); #Deltat_{13} (ns)", 100, 0, 500, 200, 0, 20);
    TH2F* triple_deltaT23_vs_energy_baseCuts = new TH2F("triple_deltaT23_vs_energy", "Calo Num vs #Deltat_{23}; E_{1} (MeV); #Deltat_{23} (ns)", 100, 0, 500, 200, 0, 20);
    TH2F* triple_deltaT23_vs_energy_lateTimes_baseCuts = new TH2F("triple_deltaT23_vs_energy_lateTimes", "Calo Num vs #Deltat_{23} after 300 #mus; E_{1} (MeV); #Deltat_{23} (ns)", 100, 0, 500, 200, 0, 20);
    TH2F* triple_deltaT12_vs_timeInFill_baseCuts = new TH2F("triple_deltaT12_vs_timeInFill", "Calo Num vs #Deltat_{12}; time (ns); #Deltat_{12} (ns)", timeBinsIn2D, histMinTime2D, histMaxTime2D, 200, 0, 20);
    TH2F* triple_deltaT13_vs_timeInFill_baseCuts = new TH2F("triple_deltaT13_vs_timeInFill", "Calo Num vs #Deltat_{13}; time (ns); #Deltat_{13} (ns)", timeBinsIn2D, histMinTime2D, histMaxTime2D, 200, 0, 20);
    TH2F* triple_deltaT23_vs_timeInFill_baseCuts = new TH2F("triple_deltaT23_vs_timeInFill", "Calo Num vs #Deltat_{23}; time (ns); #Deltat_{23} (ns)", timeBinsIn2D, histMinTime2D, histMaxTime2D, 200, 0, 20);

  auto triple_losses_baseCuts_dir = triples_baseCuts_dir->mkdir("Losses");
  triple_losses_baseCuts_dir->cd();

    TH1F* triple_losses_spectra_baseCuts = new TH1F("triple_losses_spectra", "triple_losses_spectra; time (ns); Events", nBins, histMinTime, histMaxTime);
    TH1F* triple_losses_spectra_exp_baseCuts = new TH1F("triple_losses_spectra_exp", "triple_losses_spectra_exp; time (ns); Events * e^{t/#tau}", nBins, histMinTime, histMaxTime);
    TH1F* triple_losses_spectra_integral_baseCuts = new TH1F("triple_losses_spectra_integral", "triple_losses_spectra_integral; time (ns); Integral of lost muon function", nBins, histMinTime, histMaxTime);

  auto quadruples_baseCuts_dir = baseCuts_dir->mkdir("Quadruples");
  quadruples_baseCuts_dir->cd();

  auto quadruple_losses_baseCuts_dir = quadruples_baseCuts_dir->mkdir("Losses");
  quadruple_losses_baseCuts_dir->cd();

    TH1F* quadruple_losses_spectra_baseCuts = new TH1F("quadruple_losses_spectra_baseCuts", "quadruple_losses_spectra_baseCuts; time (ns); Events", nBins, histMinTime, histMaxTime);
    TH1F* quadruple_losses_spectra_triples_baseCuts = new TH1F("quadruple_losses_spectra_triples_baseCuts", "quadruple_losses_spectra_triples_baseCuts; time (ns); Events", nBins, histMinTime, histMaxTime);

/////////////////////////////////////////////////////////////////////////////////////

    // additional cuts

    auto addedCuts_dir = topDir->mkdir("AdditionalCuts");
    addedCuts_dir->cd();

    auto listAddedCuts_dir = addedCuts_dir->mkdir("Cuts"); // make this directory so that when hadding the many copies of the lostMuonCuts get put in their own directory
    listAddedCuts_dir->cd();

      TVectorD additionalCuts(5); 
      additionalCuts[0] = deltaT13Cut;
      additionalCuts[1] = accBgCutoffLow;
      additionalCuts[2] = accBgCutoffHigh;
      additionalCuts[3] = negativeSideCutLow;
      additionalCuts[4] = negativeSideCutHigh;
      additionalCuts.Write("AdditionalCuts");

    TNamed useNegSide("useNegativeSide", BoolToString(useNegativeSide)); useNegSide.Write();

    addedCuts_dir->cd();  

    auto triples_addedCuts_dir = addedCuts_dir->mkdir("Triples");
    triples_addedCuts_dir->cd(); 

    TH2F* triple_deltaT12_vs_energy_addedCuts = new TH2F("triple_deltaT12_vs_energy", "Calo Num vs #Deltat_{12}; E_{1} (MeV); #Deltat_{12} (ns)", 100, 0, 500, 200, 0, 20);
    TH2F* triple_deltaT12_vs_energy_lateTimes_addedCuts = new TH2F("triple_deltaT12_vs_energy_lateTimes", "Calo Num vs #Deltat_{12} after 300 #mus; E_{1} (MeV); #Deltat_{12} (ns)", 100, 0, 500, 200, 0, 20);
    TH2F* triple_deltaT13_vs_energy_addedCuts = new TH2F("triple_deltaT13_vs_energy", "Calo Num vs #Deltat_{13}; E_{1} (MeV); #Deltat_{13} (ns)", 100, 0, 500, 200, 0, 20);
    TH2F* triple_deltaT13_vs_energy_lateTimes_addedCuts = new TH2F("triple_deltaT13_vs_energy_lateTimes", "Calo Num vs #Deltat_{13} after 300 #mus; E_{1} (MeV); #Deltat_{13} (ns)", 100, 0, 500, 200, 0, 20);
    TH2F* triple_deltaT23_vs_energy_addedCuts = new TH2F("triple_deltaT23_vs_energy", "Calo Num vs #Deltat_{23}; E_{1} (MeV); #Deltat_{23} (ns)", 100, 0, 500, 200, 0, 20);
    TH2F* triple_deltaT23_vs_energy_lateTimes_addedCuts = new TH2F("triple_deltaT23_vs_energy_lateTimes", "Calo Num vs #Deltat_{23} after 300 #mus; E_{1} (MeV); #Deltat_{23} (ns)", 100, 0, 500, 200, 0, 20);
    TH2F* triple_deltaT12_vs_timeInFill_addedCuts = new TH2F("triple_deltaT12_vs_timeInFill", "Calo Num vs #Deltat_{12}; time (ns); #Deltat_{12} (ns)", timeBinsIn2D, histMinTime2D, histMaxTime2D, 200, 0, 20);
    TH2F* triple_deltaT13_vs_timeInFill_addedCuts = new TH2F("triple_deltaT13_vs_timeInFill", "Calo Num vs #Deltat_{13}; time (ns); #Deltat_{13} (ns)", timeBinsIn2D, histMinTime2D, histMaxTime2D, 200, 0, 20);
    TH2F* triple_deltaT23_vs_timeInFill_addedCuts = new TH2F("triple_deltaT23_vs_timeInFill", "Calo Num vs #Deltat_{23}; time (ns); #Deltat_{23} (ns)", timeBinsIn2D, histMinTime2D, histMaxTime2D, 200, 0, 20);

    auto triple_losses_addedCuts_dir = triples_addedCuts_dir->mkdir("Losses");
    triple_losses_addedCuts_dir->cd();

      TH1F* triple_losses_spectra_addedCuts = new TH1F("triple_losses_spectra", "triple_losses_spectra; time (ns); Events", nBins, histMinTime, histMaxTime);
      TH1F* triple_losses_spectra_exp_addedCuts = new TH1F("triple_losses_spectra_exp", "triple_losses_spectra_exp; time (ns); Events * e^{t/#tau}", nBins, histMinTime, histMaxTime);
      TH1F* triple_losses_spectra_integral_addedCuts = new TH1F("triple_losses_spectra_integral", "triple_losses_spectra_integral; time (ns); Integral of lost muon function", nBins, histMinTime, histMaxTime);

    auto quadruples_addedCuts_dir = addedCuts_dir->mkdir("Quadruples");
    quadruples_addedCuts_dir->cd();

    auto quadruple_losses_addedCuts_dir = quadruples_addedCuts_dir->mkdir("Losses");
    quadruple_losses_addedCuts_dir->cd();

    TH1F* quadruple_losses_spectra_addedCuts = new TH1F("quadruple_losses_spectra_addedCuts", "quadruple_losses_spectra_addedCuts; time (ns); Events", nBins, histMinTime, histMaxTime);
    TH1F* quadruple_losses_spectra_triples_addedCuts = new TH1F("quadruple_losses_spectra_triples_addedCuts", "quadruple_losses_spectra_triples_addedCuts; time (ns); Events", nBins, histMinTime, histMaxTime);

/////////////////////////////////////////////////////////////////////////////////////

  // fill lost muons histograms

  unsigned int coincidenceLevel_Get;
  std::vector<int>* caloNum_lostMuons_Get = 0;
  std::vector<double>* clusterEnergy_Get = 0;
  std::vector<double>* clusterTime_Get = 0;
  std::vector<double>* clusterX_Get = 0;
  std::vector<double>* clusterY_Get = 0;
  std::vector<double>* clusterSize_Get = 0;
  std::vector<double>* clusterEfrac_Get = 0;

  TTree* lostMuonTree = (TTree*)inputFile->Get("testCoincidenceFinder/t");
  Long64_t lostMuonEntries = lostMuonTree->GetEntries();

  lostMuonTree->SetBranchAddress("coincidenceLevel", &coincidenceLevel_Get);
  lostMuonTree->SetBranchAddress("caloNum", &caloNum_lostMuons_Get);
  lostMuonTree->SetBranchAddress("clusterEnergy", &clusterEnergy_Get);
  lostMuonTree->SetBranchAddress("clusterTime", &clusterTime_Get);
  lostMuonTree->SetBranchAddress("clusterX", &clusterX_Get);
  lostMuonTree->SetBranchAddress("clusterY", &clusterY_Get);
  lostMuonTree->SetBranchAddress("clusterSize", &clusterSize_Get);
  lostMuonTree->SetBranchAddress("clusterEfrac", &clusterEfrac_Get);

/////////////////////////////////////////////////////////////////////////////////////

    for(Long64_t entry=0; entry<lostMuonEntries; entry++) {
       lostMuonTree->GetEntry(entry);
       
       double time1 = clusterTime_Get->at(0) * 1.25;
       if(coincidenceLevel_Get == 1 || time1 < timeCutLowMuons) continue; // skip singlets and times below 15 micro seconds

/////////////////////////////////////////////////////////////////////////////////////

       // populate coincidence information

       vector<double> times, energies, sizes, eFracs, deltaTs;

       for (uint i = 0; i < coincidenceLevel_Get; ++i)
       {
        times.push_back(clusterTime_Get->at(i) * 1.25);
        energies.push_back(clusterEnergy_Get->at(i));
        sizes.push_back(clusterSize_Get->at(i));
        eFracs.push_back(clusterEfrac_Get->at(i));
       }
       for (uint i = 0; i < coincidenceLevel_Get-1; ++i) deltaTs.push_back(times.at(i+1) - times.at(i));

/////////////////////////////////////////////////////////////////////////////////////

        // fill first hit info for all coincidences barring singlets

       firstTimes->Fill(times.at(0));
       firstSizes->Fill(sizes.at(0));
       firstEnergyFracs->Fill(eFracs.at(0));

       for (uint i = 0; i < coincidenceLevel_Get; ++i) if(sizes.at(i) >= clusterSizeCut || eFracs.at(i) < clusterEFracCut) goto treeLoop; // skip coincidences where the hits aren't MIP-like

/////////////////////////////////////////////////////////////////////////////////////

       // doublets

       if(coincidenceLevel_Get == 2)
       {
        double_deltaT12->Fill(deltaTs.at(0));
        double_Edep1->Fill(energies.at(0));
        double_Edep2->Fill(energies.at(1));

        double_calo_vs_deltaT12->Fill(deltaTs.at(0), caloNum_lostMuons_Get->at(0));

        double_deltaT12_vs_energy->Fill(energies.at(0), deltaTs.at(0));
        if(times.at(0) > 300000) double_deltaT12_vs_energy_lateTimes->Fill(energies.at(0), deltaTs.at(0));

        double_deltaT12_vs_timeInFill->Fill(times.at(0), deltaTs.at(0));
       } // end doubles

/////////////////////////////////////////////////////////////////////////////////////

       // triplets

      if(coincidenceLevel_Get == 3)
       {
        double deltaT13 = times.at(2) - times.at(0);

        triple_deltaT12->Fill(deltaTs.at(0));
        triple_deltaT13->Fill(deltaT13);
        triple_deltaT23->Fill(deltaTs.at(1));
        triple_Edep1->Fill(energies.at(0));
        triple_Edep2->Fill(energies.at(1));
        triple_Edep3->Fill(energies.at(2));

        triple_deltaT12_vs_energy->Fill(energies.at(0), deltaTs.at(0));
        triple_deltaT13_vs_energy->Fill(energies.at(0), deltaT13);
        triple_deltaT23_vs_energy->Fill(energies.at(0), deltaTs.at(1));
        if(times.at(0) > 300000){
          triple_deltaT12_vs_energy_lateTimes->Fill(energies.at(0), deltaTs.at(0));
          triple_deltaT13_vs_energy_lateTimes->Fill(energies.at(0), deltaT13);
          triple_deltaT23_vs_energy_lateTimes->Fill(energies.at(0), deltaTs.at(1));
        }

        triple_deltaT12_vs_timeInFill->Fill(times.at(0), deltaTs.at(0));
        triple_deltaT13_vs_timeInFill->Fill(times.at(0), deltaT13);
        triple_deltaT23_vs_timeInFill->Fill(times.at(0), deltaTs.at(1));


        // apply base cuts

        bool passedBaseCuts = true;

        for (uint i = 0; i < coincidenceLevel_Get-1; ++i) if(deltaTs.at(i) < deltaTCutLow || deltaTs.at(i) > deltaTCutHigh) passedBaseCuts = false;
        for (uint i = 0; i < coincidenceLevel_Get; ++i) if(energies.at(i) < energyCutLow || energies.at(i) > energyCutHigh) passedBaseCuts = false;

        if(passedBaseCuts == true)
        {

          triple_deltaT12_vs_energy_baseCuts->Fill(energies.at(0), deltaTs.at(0));
          triple_deltaT13_vs_energy_baseCuts->Fill(energies.at(0), deltaT13);
          triple_deltaT23_vs_energy_baseCuts->Fill(energies.at(0), deltaTs.at(1));
          if(times.at(0) > 300000){
            triple_deltaT12_vs_energy_lateTimes_baseCuts->Fill(energies.at(0), deltaTs.at(0));
            triple_deltaT13_vs_energy_lateTimes_baseCuts->Fill(energies.at(0), deltaT13);
            triple_deltaT23_vs_energy_lateTimes_baseCuts->Fill(energies.at(0), deltaTs.at(1));
          }

          triple_deltaT12_vs_timeInFill_baseCuts->Fill(times.at(0), deltaTs.at(0));
          triple_deltaT13_vs_timeInFill_baseCuts->Fill(times.at(0), deltaT13);
          triple_deltaT23_vs_timeInFill_baseCuts->Fill(times.at(0), deltaTs.at(1));

          triple_losses_spectra_baseCuts->Fill(times.at(0));


            // apply additional cuts

            if(deltaT13 > deltaT13Cut) goto endCuts;
            if(useNegativeSide && (deltaT13 < negativeSideCutLow || deltaT13 > negativeSideCutHigh)) goto endCuts;

              triple_deltaT12_vs_energy_addedCuts->Fill(energies.at(0), deltaTs.at(0));
              triple_deltaT13_vs_energy_addedCuts->Fill(energies.at(0), deltaT13);
              triple_deltaT23_vs_energy_addedCuts->Fill(energies.at(0), deltaTs.at(1));
              if(times.at(0) > 300000){
                triple_deltaT12_vs_energy_lateTimes_addedCuts->Fill(energies.at(0), deltaTs.at(0));
                triple_deltaT13_vs_energy_lateTimes_addedCuts->Fill(energies.at(0), deltaT13);
                triple_deltaT23_vs_energy_lateTimes_addedCuts->Fill(energies.at(0), deltaTs.at(1));
              }

              triple_deltaT12_vs_timeInFill_addedCuts->Fill(times.at(0), deltaTs.at(0));
              triple_deltaT13_vs_timeInFill_addedCuts->Fill(times.at(0), deltaT13);
              triple_deltaT23_vs_timeInFill_addedCuts->Fill(times.at(0), deltaTs.at(1));

              triple_losses_spectra_addedCuts->Fill(times.at(0));

            endCuts: ;

        } // base cuts

       }  // end triples

/////////////////////////////////////////////////////////////////////////////////////

      // quadruplets - only make to subtract off from triplets 

      if(coincidenceLevel_Get == 4)
       {
        // fill loss spectra

        double deltaT13 = times.at(2) - times.at(0);
        double deltaT24 = times.at(3) - times.at(1);

        bool passedBaseCuts = true;

        for (uint i = 0; i < coincidenceLevel_Get-1; ++i) if(deltaTs.at(i) < deltaTCutLow || deltaTs.at(i) > deltaTCutHigh) passedBaseCuts = false;
        for (uint i = 0; i < coincidenceLevel_Get; ++i) if(energies.at(i) < energyCutLow || energies.at(i) > energyCutHigh) passedBaseCuts = false;

        if(passedBaseCuts == true){
            quadruple_losses_spectra_baseCuts->Fill(times.at(0));
            quadruple_losses_spectra_triples_baseCuts->Fill(times.at(0));
            quadruple_losses_spectra_triples_baseCuts->Fill(times.at(1));

            if(deltaT13 > deltaT13Cut || deltaT24 > deltaT13Cut) goto endCuts_quads;
            if(useNegativeSide && (deltaT13 < negativeSideCutLow || deltaT13 > negativeSideCutHigh || deltaT24 < negativeSideCutLow || deltaT24 > negativeSideCutHigh)) goto endCuts_quads;

            quadruple_losses_spectra_addedCuts->Fill(times.at(0));
            quadruple_losses_spectra_triples_addedCuts->Fill(times.at(0));
            quadruple_losses_spectra_triples_addedCuts->Fill(times.at(1));

            endCuts_quads: ;

        }
       }  // end quadruples


/////////////////////////////////////////////////////////////////////////////////////

      treeLoop: ; // goto point to skip non-MIP-like coincidences

     } // end for loop over lost muons tree

/////////////////////////////////////////////////////////////////////////////////////

    // make adjustments for additional cuts plots

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

    // subtract off quadruples (by subtracting off two triples)

    auto losses_minus_quads_dir = triples_addedCuts_dir->mkdir("Losses_minus_quadruples");
    losses_minus_quads_dir->cd();

    TH1F* triple_losses_spectra_minus_quadruples = (TH1F*) triple_losses_spectra_addedCuts->Clone("triple_losses_spectra_minus_quadruples");
    triple_losses_spectra_minus_quadruples->Add(quadruple_losses_spectra_triples_addedCuts, -1);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

    // subtract off accidentals (and make plots)

    auto accBgInfo_dir = triplesInfo_dir->mkdir("AccidentalBackground");
    accBgInfo_dir->cd();

    TH2F* accidentalBackground = (TH2F*) triple_deltaT12_vs_timeInFill->Clone("triple_deltaT12_vs_timeInFill_accBg");
    accidentalBackground->Reset();

    // grab times within some cut off for background
    int yBinCutoffLow = accidentalBackground->GetYaxis()->FindBin(accBgCutoffLow);
    int yBinCutoffHigh = accidentalBackground->GetYaxis()->FindBin(accBgCutoffHigh) - 1;

    for (int yBin = yBinCutoffLow; yBin <= yBinCutoffHigh; ++yBin){
      for (int xBin = 1; xBin <= accidentalBackground->GetXaxis()->GetNbins(); ++xBin){
        accidentalBackground->SetBinContent(xBin, yBin, triple_deltaT12_vs_timeInFill->GetBinContent(xBin, yBin));
      }
    }

    // average in time slices

    TH2F* accidentalBackground_average = (TH2F*) accidentalBackground->Clone("triple_deltaT12_vs_timeInFill_accBg_Avg");
    int yCounts = yBinCutoffHigh - yBinCutoffLow + 1;

    for (int xBin = 1; xBin <= accidentalBackground_average->GetXaxis()->GetNbins(); ++xBin){

      double ySum = 0, yAverage = 0;

        for (int yBin = yBinCutoffLow; yBin <= yBinCutoffHigh; ++yBin) ySum += accidentalBackground_average->GetBinContent(xBin, yBin);
        yAverage = ySum/yCounts;
        for (int yBin = 1; yBin <= accidentalBackground_average->GetYaxis()->GetNbins(); ++yBin) accidentalBackground_average->SetBinContent(xBin, yBin, yAverage);
    }

    // make 1D accidentals time spectrum

    TH1F* accBg_avg_1D = new TH1F("Accidental_Background_Average_1D", "Accidental Background Average; time (ns); Events", timeBinsIn2D, histMinTime2D, histMaxTime2D);

    int numYbinsInLossHist = accidentalBackground->GetYaxis()->FindBin(deltaTCutHigh) - accidentalBackground->GetYaxis()->FindBin(deltaTCutLow);

    for (int xBin = 1; xBin <= accidentalBackground_average->GetXaxis()->GetNbins(); ++xBin){
      double stripeTime = accidentalBackground_average->GetXaxis()->GetBinCenter(xBin);
      double stripeContent = accidentalBackground_average->GetBinContent(xBin, 1);

      accBg_avg_1D->Fill(stripeTime, stripeContent*numYbinsInLossHist);
    }

    // subtract off average from info 2D histogram for visualization purposes

    TH2F* triple_deltaT12_vs_timeInFill_subtracted = (TH2F*) triple_deltaT12_vs_timeInFill->Clone("triple_deltaT12_vs_timeInFill_accBg_Subtracted");
    triple_deltaT12_vs_timeInFill_subtracted->Add(accidentalBackground_average, -1);

    // subtract from 1D losses spectra (don't need to subtract from additional cuts 2D spectra since the changes are minimal and those are only for plotting purposes)

    auto losses_minus_accidentals_dir = triples_addedCuts_dir->mkdir("Losses_minus_accidentals");
    losses_minus_accidentals_dir->cd();

    TH1F* triple_losses_spectra_minus_accidentals = (TH1F*) triple_losses_spectra_addedCuts->Clone("triple_losses_spectra_minus_accidentals");
    triple_losses_spectra_minus_accidentals->Add(accBg_avg_1D, -1);

/////////////////////////////////////////////////////////////////////////////////////

    // make lost spectra with quadruple and accidental background subtraction

    auto losses_minus_quads_accidentals_dir = triples_addedCuts_dir->mkdir("Losses_minus_quadruples_accidentals");
    losses_minus_quads_accidentals_dir->cd();

    TH1F* triple_losses_spectra_minus_quadruples_accidentals = (TH1F*) triple_losses_spectra_minus_quadruples->Clone("triple_losses_spectra_minus_quadruples_accidentals");
    triple_losses_spectra_minus_quadruples_accidentals->Add(accBg_avg_1D, -1);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

    // fill out integral histograms - only for viewing purposes - the integral histograms are regenerated in ratioMacro.C so that the start time can be adjusted on the fly

    int integralStartBin = triple_losses_spectra_integral_baseCuts->FindBin(30200);

    // base cuts

    for (int bin = 1; bin <= nBins; ++bin) triple_losses_spectra_exp_baseCuts->SetBinContent(bin, triple_losses_spectra_baseCuts->GetBinContent(bin) * exp(triple_losses_spectra_baseCuts->GetBinCenter(bin)/defaultLifetime));
    for (int bin = integralStartBin; bin <= nBins; ++bin) triple_losses_spectra_integral_baseCuts->SetBinContent(bin, triple_losses_spectra_exp_baseCuts->Integral(integralStartBin, bin));

    // additional cuts

    for (int bin = 1; bin <= nBins; ++bin) triple_losses_spectra_exp_addedCuts->SetBinContent(bin, triple_losses_spectra_addedCuts->GetBinContent(bin) * exp(triple_losses_spectra_addedCuts->GetBinCenter(bin)/defaultLifetime));
    for (int bin = integralStartBin; bin <= nBins; ++bin) triple_losses_spectra_integral_addedCuts->SetBinContent(bin, triple_losses_spectra_exp_addedCuts->Integral(integralStartBin, bin));

    // quadruple subtraction integral

    losses_minus_quads_dir->cd();

    TH1F* triple_losses_spectra_minus_quadruples_exp = new TH1F("triple_losses_spectra_exp", "triple_losses_spectra_exp; time (ns); Events * e^{t/#tau}", nBins, histMinTime, histMaxTime);
    TH1F* triple_losses_spectra_minus_quadruples_integral = new TH1F("triple_losses_spectra_integral", "triple_losses_spectra_integral; time (ns); Integral of lost muon function", nBins, histMinTime, histMaxTime);

    for (int bin = 1; bin <= nBins; ++bin) triple_losses_spectra_minus_quadruples_exp->SetBinContent(bin, triple_losses_spectra_minus_quadruples->GetBinContent(bin) * exp(triple_losses_spectra_minus_quadruples->GetBinCenter(bin)/defaultLifetime));
    for (int bin = integralStartBin; bin <= nBins; ++bin) triple_losses_spectra_minus_quadruples_integral->SetBinContent(bin, triple_losses_spectra_minus_quadruples_exp->Integral(integralStartBin, bin));

    // accidental background subtraction integral

    losses_minus_accidentals_dir->cd();

    TH1F* triple_losses_spectra_minus_accidentals_exp = new TH1F("triple_losses_spectra_exp", "triple_losses_spectra_exp; time (ns); Events * e^{t/#tau}", nBins, histMinTime, histMaxTime);
    TH1F* triple_losses_spectra_minus_accidentals_integral = new TH1F("triple_losses_spectra_integral", "triple_losses_spectra_integral; time (ns); Integral of lost muon function", nBins, histMinTime, histMaxTime);

    for (int bin = 1; bin <= nBins; ++bin) triple_losses_spectra_minus_accidentals_exp->SetBinContent(bin, triple_losses_spectra_minus_accidentals->GetBinContent(bin) * exp(triple_losses_spectra_minus_accidentals->GetBinCenter(bin)/defaultLifetime));
    for (int bin = integralStartBin; bin <= nBins; ++bin) triple_losses_spectra_minus_accidentals_integral->SetBinContent(bin, triple_losses_spectra_minus_accidentals_exp->Integral(integralStartBin, bin));

    // quadruple and accidental background subtraction integral

    losses_minus_quads_accidentals_dir->cd();

    TH1F* triple_losses_spectra_minus_quadruples_accidentals_exp = new TH1F("triple_losses_spectra_exp", "triple_losses_spectra_exp; time (ns); Events * e^{t/#tau}", nBins, histMinTime, histMaxTime);
    TH1F* triple_losses_spectra_minus_quadruples_accidentals_integral = new TH1F("triple_losses_spectra_integral", "triple_losses_spectra_integral; time (ns); Integral of lost muon function", nBins, histMinTime, histMaxTime);

    for (int bin = 1; bin <= nBins; ++bin) triple_losses_spectra_minus_quadruples_accidentals_exp->SetBinContent(bin, triple_losses_spectra_minus_quadruples_accidentals->GetBinContent(bin) * exp(triple_losses_spectra_minus_quadruples_accidentals->GetBinCenter(bin)/defaultLifetime));
    for (int bin = integralStartBin; bin <= nBins; ++bin) triple_losses_spectra_minus_quadruples_accidentals_integral->SetBinContent(bin, triple_losses_spectra_minus_quadruples_accidentals_exp->Integral(integralStartBin, bin));


/////////////////////////////////////////////////////////////////////////////////////

   outputFile->Write();
   delete outputFile;

/////////////////////////////////////////////////////////////////////////////////////

   return 0;

}
