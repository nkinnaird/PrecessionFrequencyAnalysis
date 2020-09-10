// 3-31-20: Macro to make histograms from ROOT TTrees of clusters.
// To run with an input file use: root ClusterTreeToHistsOrdered.C+\(\"filepath\"\)

#include <iostream>
#include <iomanip>
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
#include "../ratioMacroHeaders/ratioHistogramsConfig.hh" // don't change these paths since this file is run on the grid (whereas other files can see my root include directory and the ../ can be removed)
#include "../ratioMacroHeaders/pileupUtils.hh"
// #include "../ratioMacroHeaders/ratioAnalysisConfig.hh" // include for tracker model parameters when applying time dependent randomization

/////////////////////////////////////////////////////////////////////////////////////

bool sortTuple(const tuple<double, uint, double, double, double, double, uint, uint, vector<uint>, vector<double>, vector<double>, vector<double>>& a,  
               const tuple<double, uint, double, double, double, double, uint, uint, vector<uint>, vector<double>, vector<double>, vector<double>>& b) 
{ 
    return (get<2>(a) < get<2>(b)); 
} 

/////////////////////////////////////////////////////////////////////////////////////

double getTVW(double t, double inputBinWidth)
{
  double period;
  double fc = 1/inputBinWidth;

/////////////////////////////////////////////////////////////////////////////////////
// time dependent 2fy and VW frequencies - need the analysis config file commented in

  // double fcbo = (1/(2*pi)) * (tM_w0 + tM_Aconst*exp(-t/tM_Atau)/t + tM_Bconst*exp(-t/tM_Btau)/t);
  // double two_fy = 2. * sqrt(2*fc*fcbo - fcbo*fcbo);

  // period = 1./two_fy;
  // return period;

  // double timeDep_fVW = fc - two_fy;
  // period = 1./timeDep_fVW;
  // return 1/timeDep_fVW;
/////////////////////////////////////////////////////////////////////////////////////

// constant VW frequency
  double const_fVW = (1-2*sqrt(nVal)) * fc;
  period = 1./const_fVW;

  return period;
}


double adHocGainCorrection(double* x, double* par)
{
  double time = x[0];
  double gainValue = par[0] * exp(-time/par[1]) * (1 + par[2] * cos(adHoc_w_a * time + adHoc_phi));

  return gainValue;
}


uint getCrystalRow(uint nHits, uint* xtalNumIn, float* energies){
  double largestEnergy = 0;
  int hitWithGreatestEnergy = -1;

  for (uint i = 0; i < nHits; ++i){
    if(energies[i] > largestEnergy){
      largestEnergy = energies[i];
      hitWithGreatestEnergy = i;
    }
  }
 
  return (1 + xtalNumIn[hitWithGreatestEnergy]/9);
}

uint getRunGroup(uint runNumber){
  uint runGroup = 1000; // some high number

  for (uint i = 0; i < runGroups.size()-1; ++i){
    if(runNumber >= runGroups.at(i) && runNumber < runGroups.at(i+1)){
      runGroup = i;
      break;
    }
  }

  return runGroup;
}

/////////////////////////////////////////////////////////////////////////////////////

int ClusterTreeToHistsOrdered(std::string filePath)
{
  cout << "Making histograms from cluster tree: " << filePath << endl;

  // check config before analyzing
  int configCheck = 0;

  if(varyRandSeeds) configCheck++;
  if(binWidthStep) configCheck++;
  if(binEdgeShiftStep) configCheck++;
  if(bunchNumScan) configCheck++;
  if(crystalRowScan) configCheck++;
  if(separateByRunGroup) configCheck++;
  if(gm2PeriodPPMStep) configCheck++;
  if(energyThresholdStep || makeEnergyBinnedHists) configCheck++;
  if(ADTStep) configCheck++;
  if(SDTStep) configCheck++;
  if(SGTStep) configCheck++;
  if(SDT2Step) configCheck++;
  if(pileupTimeShiftStep) configCheck++;
  if(pileupEnergyScalingStep) configCheck++;

  if(recorrectGainAtCrystals){
    if(crystalGainAmpFactorStep) configCheck++;
    if(crystalGainTauFactorStep) configCheck++;
  }

  if(applyAdHocGainCorrection){
    if(adHocAmplitudeStep) configCheck++;
    if(adHocLifetimeStep) configCheck++;
    if(adHocAsymmetryStep) configCheck++;
  }

  if(deltaTCutLowStep) configCheck++;
  if(deltaTCutHighStep) configCheck++;
  if(energyCutLowStep) configCheck++;
  if(energyCutHighStep) configCheck++;

  if(configCheck > 1){
    cout << "Histogram making configuration set incorrectly." << endl;
    return 1;
  }

  if(randomize_VW && !randomize_fc){
    cout << "Randomizing by VW but not randomizing by fc, which is odd." << endl;
    return 1;
  }

  if(crystalRowScan){
    cout << "Before doing a crystal row scan, need to reimplement and update the previous code." << endl;
    return 1;
  }

  if(setTimeShiftHalfGapTime && pileupTimeShiftStep){
    cout << "Setting the pileup time shift as half the gap time but trying to scan over the time shift." << endl;
    return 1;    
  }

/////////////////////////////////////////////////////////////////////////////////////

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 1;
   }

  clock_t t; // for time profiling - put after input
  t = clock();

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

  string outputFileName = "clusterHists_ordered_" + uniqueID + ".root";

  cout << "Output file name: " << outputFileName << endl;

  TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  auto topDir = gFile->mkdir("topDir");

/////////////////////////////////////////////////////////////////////////////////////

  uint totalIters = totalIterations;
  if(bunchNumScan) totalIters = 9;
  else if(crystalRowScan) totalIters = 7;
  else if(makeEnergyBinnedHists){
    binHistsByEnergy(makeEnergyBinnedHists); // call method to adjust energy threshold parameters, defined in ratioHistogramsConfig.hh
    totalIters = (energyBinMax - lowerEnergyThresholdStart) / energyThresholdStep;
  }
  else if(separateByRunGroup) totalIters = runGroups.size(); // size-1 individual run groups plus all run groups together for iteration 0

  cout << "Total iterations: " << totalIters << endl;

  TVectorD numIters(1);
  numIters[0] = totalIters;
  numIters.Write("Iters");

  TDirectory* histIterDirs[totalIters];

/////////////////////////////////////////////////////////////////////////////////////

  double binWidths[totalIters];
  double binEdgeShifts[totalIters];

  double lowerEnergyThresholds[totalIters];
  double upperEnergyThresholds[totalIters];

  double gm2PeriodPPMOffsets[totalIters];
  double gm2PeriodGuesses[totalIters];
  double weightingLifetimes[totalIters];

  double artificialDeadTimes[totalIters];
  double shadowDeadTimes[totalIters];
  double shadowGapTimes[totalIters];
  double shadowSecondDeadTimes[totalIters];

  double pileupTimeShifts[totalIters];
  double pileupEnergyScales[totalIters];

  double crystalGainAmpFactors[totalIters];
  double crystalGainTauFactors[totalIters];

/////////////////////////////////////////////////////////////////////////////////////

  TF1* adHocGainFunctions[totalIters];

  double adHocAmplitudes[totalIters];
  double adHocLifetimes[totalIters];
  double adHocAsymmetries[totalIters];

/////////////////////////////////////////////////////////////////////////////////////

  // lost muon cuts

  double deltaTCutLows[totalIters];
  double deltaTCutHighs[totalIters];
  double energyCutLows[totalIters];
  double energyCutHighs[totalIters];
  double clusterSizeCuts[totalIters];
  double clusterEFracCuts[totalIters];

/////////////////////////////////////////////////////////////////////////////////////

  for (uint iter = 0; iter < totalIters; ++iter)
  {
    crystalGainAmpFactors[iter] = crystalGainAmpFactorStart + iter * crystalGainAmpFactorStep;
    crystalGainTauFactors[iter] = crystalGainTauFactorStart + iter * crystalGainTauFactorStep;

    adHocAmplitudes[iter] = adHocAmplitudeStart + iter * adHocAmplitudeStep;
    adHocLifetimes[iter] = adHocLifetimeStart + iter * adHocLifetimeStep;
    adHocAsymmetries[iter] = adHocAsymmetryStart + iter * adHocAsymmetryStep;

/////////////////////////////////////////////////////////////////////////////////////

    binWidths[iter] = binWidthStart + iter * binWidthStep;
    binEdgeShifts[iter] = binEdgeShiftStart + iter * binEdgeShiftStep;

    gm2PeriodPPMOffsets[iter] = gm2PeriodPPMStart + iter * gm2PeriodPPMStep;
    gm2PeriodGuesses[iter] = g2Period * (1 + 1e-6 * gm2PeriodPPMOffsets[iter]);

    weightingLifetimes[iter] = weightingLifetimeStart + iter * weightingLifetimeStep;

    lowerEnergyThresholds[iter] = lowerEnergyThresholdStart + iter * energyThresholdStep;
    upperEnergyThresholds[iter] = upperEnergyThresholdStart + iter * energyThresholdStep;
    
    artificialDeadTimes[iter] = ADTStart + iter * ADTStep;
    shadowDeadTimes[iter] = SDTStart + iter * SDTStep;
    shadowGapTimes[iter] = SGTStart + iter * SGTStep;
    shadowSecondDeadTimes[iter] = SDT2Start + iter * SDT2Step;

    if(setTimeShiftHalfGapTime) pileupTimeShifts[iter] = 0.5 * shadowGapTimes[iter];
    else pileupTimeShifts[iter] = pileupTimeShiftStart + iter * pileupTimeShiftStep;

    pileupEnergyScales[iter] = pileupEnergyScalingStart + iter * pileupEnergyScalingStep;

/////////////////////////////////////////////////////////////////////////////////////

    deltaTCutLows[iter] = deltaTCutLow + iter * deltaTCutLowStep;
    deltaTCutHighs[iter] = deltaTCutHigh + iter * deltaTCutHighStep;
    energyCutLows[iter] = energyCutLow + iter * energyCutLowStep;
    energyCutHighs[iter] = energyCutHigh + iter * energyCutHighStep;
    clusterSizeCuts[iter] = clusterSizeCut;
    clusterEFracCuts[iter] = clusterEFracCut;

  }

/////////////////////////////////////////////////////////////////////////////////////

  // define histograms

  auto totalInfoDir = topDir->mkdir("TotalInfo");
  totalInfoDir->cd();

  // these hists will be filled by incoming hits before they are separated or adjusted due to the different iteration conditions

  TH1F* runHist = new TH1F("Runs", "Runs; Run Number; Events", 15000, 15000, 30000);
  TH1F* subRunHist = new TH1F("SubRuns", "SubRuns; SubRun Number; Events", 510, 0, 510);
  TH1F* numFills = new TH1F("Fills", "Fills; ; Events", 1, 0, 1);

/////////////////////////////////////////////////////////////////////////////////////

  TH1F* runGroupHist[totalIters];
  TH1F* runNumHist_perIter[totalIters];

  TH1F* bunchNumHist[totalIters];
  TH1F* crystalRowHist[totalIters];

  TH1F* timeShifts[totalIters];
  TH1F* UVrands[totalIters];

  TH1F* allEnergies[totalIters];
    TH1F* allEnergies_early[totalIters];
    TH1F* allEnergies_late[totalIters];
  TH1F* thresholdEnergies[totalIters];

  TH1F* allTimes[totalIters];
  TH1F* thresholdTimes[totalIters];

  TH1F* allTimesAdded_U[totalIters];
  TH1F* allTimesAdded_V[totalIters];

    TH1F* caloEnergies[totalIters][nCalos];
    TH1F* caloEnergiesThreshold[totalIters][nCalos];
    TH1F* caloTimes[totalIters][nCalos];
    TH1F* caloTimesThreshold[totalIters][nCalos];
    TH1F* caloTimesThresholdU[totalIters][nCalos];
    TH1F* caloTimesThresholdV[totalIters][nCalos];

  TH2F* timesAndEnergies[totalIters];
  TH2F* caloTimesAndEnergies[totalIters][nCalos];

        TH1F* caloEnergies_early[totalIters][nCalos];
        TH1F* caloEnergies_late[totalIters][nCalos];

/////////////////////////////////////////////////////////////////////////////////////

    // pileup histograms

     TH1F* pileupTimes_shadow[totalIters][nCalos];
     TH1F* pileupTimes_shadow_threshold[totalIters][nCalos];
     TH1F* pileupEnergies_shadow[totalIters][nCalos];
     TH1F* pileupEnergiesThreshold_shadow[totalIters][nCalos];
     TH1F* pileupTimes_U_shadow[totalIters][nCalos];
     TH1F* pileupTimes_V_shadow[totalIters][nCalos];

     TH1F* pileupErrorsNeither_shadow[totalIters][nCalos];
     TH1F* pileupErrorsBoth_shadow[totalIters][nCalos];

     TH1F* deltaT_all[totalIters];     
     TH1F* calo_deltaT[totalIters][nCalos];

        TH1F* pileupEnergies_shadow_early[totalIters][nCalos];
        TH1F* pileupEnergies_shadow_late[totalIters][nCalos];

/////////////////////////////////////////////////////////////////////////////////////

    // lost muons

    TH1F* DeltaT_triple[totalIters];
    TH1F* Edep_triple[totalIters];
    TH1F* Losses_triple[totalIters];

/////////////////////////////////////////////////////////////////////////////////////

for (uint iter = 0; iter < totalIters; ++iter)
{
  int nBins = int(approxMaxTime/binWidths[iter]);
  double histMaxTime = nBins*binWidths[iter] + binEdgeShifts[iter];
  double histMinTime = binEdgeShifts[iter];

/////////////////////////////////////////////////////////////////////////////////////

  histIterDirs[iter] = topDir->mkdir(Form("Iter%d", iter));
  histIterDirs[iter]->cd();

  auto savDirs = histIterDirs[iter]->mkdir("SavedParameters");
  savDirs->cd();

  TVectorD parameterStore(17); // vector of doubles to store parameters used when creating histograms, for systematic studies

  parameterStore[0] = gm2PeriodGuesses[iter];
  parameterStore[1] = lowerEnergyThresholds[iter];
  parameterStore[2] = shadowDeadTimes[iter];
  parameterStore[3] = shadowGapTimes[iter];
  parameterStore[4] = artificialDeadTimes[iter];
  parameterStore[5] = shadowSecondDeadTimes[iter];
  parameterStore[6] = pileupTimeShifts[iter];
  parameterStore[7] = pileupEnergyScales[iter];
  parameterStore[8] = weightingLifetimes[iter];
  parameterStore[9] = binWidths[iter];
  parameterStore[10] = binEdgeShifts[iter];
  parameterStore[11] = crystalGainAmpFactors[iter];
  parameterStore[12] = crystalGainTauFactors[iter];
  parameterStore[13] = adHocAmplitudes[iter];
  parameterStore[14] = adHocLifetimes[iter];
  parameterStore[15] = adHocAsymmetries[iter];
  parameterStore[16] = upperEnergyThresholds[iter];
  
  parameterStore.Write("parameterStore");

  // The only thing that is not technically correct about this saved function is that the w_a and phi in the function are updated based on the dataset, which is only set after the tree has been read for the first time, before applying to the clusters.
  adHocGainFunctions[iter] = new TF1("adHocGainCorrection", adHocGainCorrection, 0, histMaxTime, 3);
  adHocGainFunctions[iter]->SetNpx(10000);
  adHocGainFunctions[iter]->SetParameter(0, adHocAmplitudes[iter]);
  adHocGainFunctions[iter]->SetParameter(1, adHocLifetimes[iter]);
  adHocGainFunctions[iter]->SetParameter(2, adHocAsymmetries[iter]);
  adHocGainFunctions[iter]->Write();

/////////////////////////////////////////////////////////////////////////////////////  
/////////////////////////////////////////////////////////////////////////////////////

    auto iterInfoDir = histIterDirs[iter]->mkdir("IterInfo");
    iterInfoDir->cd();

      runGroupHist[iter] = new TH1F("RunGroup", "RunGroup; Run Group; Events", runGroups.size()-1, 0, runGroups.size()-1);
      runNumHist_perIter[iter] = new TH1F("RunsInIter", "Runs; Run Number; Events", 15000, 15000, 30000);

      bunchNumHist[iter] = new TH1F("BunchNum", "Bunch Number; Bunch Number; Events", 8, 0, 8);
      crystalRowHist[iter] = new TH1F("CrystalRow", "Crystal Row; Crystal Row; Events", 6, 1, 7);

      timeShifts[iter] = new TH1F("TimeShifts", "Time Shifts (all events); Time [ns]; Events", 1000, -500, 500);
      UVrands[iter] = new TH1F("UVrands", "UV Random Numbers; ; Events", 100, 0, 1);

    auto pileupDir = histIterDirs[iter]->mkdir("PileupHists");
    pileupDir->cd();

    auto otherDir = pileupDir->mkdir("Other");
    otherDir->cd();

    deltaT_all[iter] = new TH1F("Time_Between_Hits_All", "Time_Between_Hits_All; #DeltaT [ns]; Events", 2000, 0, 1000); // .5 ns bins
  
    auto calosDirPileup = pileupDir->mkdir("Calos");
    TDirectory* pileupCaloDirs[nCalos];

    for (int caloNum = 1; caloNum <= nCalos; ++caloNum)
    {
      pileupCaloDirs[caloNum-1] = calosDirPileup->mkdir(Form("Calo%d",caloNum));
      pileupCaloDirs[caloNum-1]->cd();

      pileupTimes_shadow[iter][caloNum-1] = new TH1F(Form("pileupTimes_shadow_Calo%d", caloNum), Form("pileupTimes_shadow_Calo%d", caloNum), nBins, histMinTime, histMaxTime);
      pileupTimes_shadow_threshold[iter][caloNum-1] = new TH1F(Form("pileupTimes_shadow_threshold_Calo%d", caloNum), Form("pileupTimes_shadow_threshold_Calo%d", caloNum), nBins, histMinTime, histMaxTime);

      pileupEnergies_shadow[iter][caloNum-1] = new TH1F(Form("pileupEnergies_shadow_Calo%d", caloNum), Form("pileupEnergies_shadow_Calo%d", caloNum), 1000, 0, 10000);
      pileupEnergies_shadow_early[iter][caloNum-1] = new TH1F(Form("pileupEnergies_shadow_early_Calo%d", caloNum), Form("pileupEnergies_shadow_early_Calo%d", caloNum), 1000, 0, 10000);
      pileupEnergies_shadow_late[iter][caloNum-1] = new TH1F(Form("pileupEnergies_shadow_late_Calo%d", caloNum), Form("pileupEnergies_shadow_late_Calo%d", caloNum), 1000, 0, 10000);
      pileupEnergiesThreshold_shadow[iter][caloNum-1] = new TH1F(Form("pileupEnergies_shadow_Calo%d_Threshold", caloNum), Form("pileupEnergies_shadow_Calo%d", caloNum), 1000, 0, 10000);

      pileupTimes_U_shadow[iter][caloNum-1] = new TH1F(Form("pileupTimes_U_shadow_Calo%d", caloNum), Form("pileupTimes_U_shadow_Calo%d", caloNum), nBins, histMinTime, histMaxTime);
      pileupTimes_V_shadow[iter][caloNum-1] = new TH1F(Form("pileupTimes_V_shadow_Calo%d", caloNum), Form("pileupTimes_V_shadow_Calo%d", caloNum), nBins, histMinTime, histMaxTime);

      pileupErrorsNeither_shadow[iter][caloNum-1] = new TH1F(Form("pileupErrorsNeither_shadow_Calo%d", caloNum), Form("pileupErrorsNeither_shadow_Calo%d", caloNum), nBins, histMinTime, histMaxTime);
      pileupErrorsBoth_shadow[iter][caloNum-1] = new TH1F(Form("pileupErrorsBoth_shadow_Calo%d", caloNum), Form("pileupErrorsBoth_shadow_Calo%d", caloNum), nBins, histMinTime, histMaxTime);

      calo_deltaT[iter][caloNum-1] = new TH1F(Form("Time_Between_Hits_Calo%d", caloNum), Form("Time_Between_Hits_Calo%d; #DeltaT [ns]; Events", caloNum), 2000, 0, 1000); // .5 ns bins
    }

/////////////////////////////////////////////////////////////////////////////////////

  auto addedDir = histIterDirs[iter]->mkdir("Added");
  addedDir->cd();

    allEnergies[iter] = new TH1F("Energy", "Energies; Energy [MeV]; Events", 1000, 0, 10000);
      allEnergies_early[iter] = new TH1F("Energy_early", "Energies Early; Energy [MeV]; Events", 1000, 0, 10000);
      allEnergies_late[iter] = new TH1F("Energy_late", "Energies Late; Energy [MeV]; Events", 1000, 0, 10000);
    thresholdEnergies[iter] = new TH1F("Energy_Threshold", "Energies; Energy [MeV]; Events", 1000, 0, 10000);

    allTimes[iter] = new TH1F("Times", "; time [ns]; Events", nBins, histMinTime, histMaxTime);
    thresholdTimes[iter] = new TH1F("Times_E_Threshold", "; time [ns]; Events", nBins, histMinTime, histMaxTime);

    if(!reduceMemory) timesAndEnergies[iter] = new TH2F("Times_And_Energies", "; time [ns]; Energy [MeV]; test", nBins, histMinTime, histMaxTime, 1000, 0, 10000 );

    allTimesAdded_U[iter] = new TH1F("allTimesAdded_U", "allTimesAdded_U", nBins, histMinTime, histMaxTime);
    allTimesAdded_V[iter] = new TH1F("allTimesAdded_V", "allTimesAdded_V", nBins, histMinTime, histMaxTime);

/////////////////////////////////////////////////////////////////////////////////////

  auto calosDirectory = histIterDirs[iter]->mkdir("Calos");
  calosDirectory->cd();

    TDirectory* caloDirs[nCalos];

    for (int caloNum = 1; caloNum <= nCalos; ++caloNum)
    {
      caloDirs[caloNum-1] = calosDirectory->mkdir(Form("Calo%d",caloNum));
      caloDirs[caloNum-1]->cd();

      caloEnergies[iter][caloNum-1] = new TH1F( Form("Calo%d_Energy", caloNum), "; Energy; Events", 1000, 0, 10000);
        caloEnergies_early[iter][caloNum-1] = new TH1F( Form("Calo%d_Energy_early", caloNum), "; Energy; Events", 1000, 0, 10000);
        caloEnergies_late[iter][caloNum-1] = new TH1F( Form("Calo%d_Energy_late", caloNum), "; Energy; Events", 1000, 0, 10000);
      caloEnergiesThreshold[iter][caloNum-1] = new TH1F( Form("Calo%d_Energy_Threshold", caloNum), "; Energy; Events", 1000, 0, 10000);

      caloTimes[iter][caloNum-1] = new TH1F( Form("Calo%d_Times", caloNum), "; time [ns]; Events", nBins, histMinTime, histMaxTime);
      caloTimesThreshold[iter][caloNum-1] = new TH1F( Form("Calo%d_Times_E_Threshold", caloNum), "; time [ns]; Events", nBins, histMinTime, histMaxTime);

      if(!reduceMemory) caloTimesAndEnergies[iter][caloNum-1] = new TH2F( Form("Calo%d_Times_And_Energies", caloNum), "; time [ns]; Energy [MeV]; test", nBins, histMinTime, histMaxTime, 1000, 0, 10000 );

      caloTimesThresholdU[iter][caloNum-1] = new TH1F( Form("Calo%d_Times_E_Threshold_U", caloNum), "; time [ns]; Events", nBins, histMinTime, histMaxTime);
      caloTimesThresholdV[iter][caloNum-1] = new TH1F( Form("Calo%d_Times_E_Threshold_V", caloNum), "; time [ns]; Events", nBins, histMinTime, histMaxTime);
    }

/////////////////////////////////////////////////////////////////////////////////////

    // lost muons

    auto lostMuonsDir = histIterDirs[iter]->mkdir("LostMuons");
    lostMuonsDir->cd();

    TVectorD lostMuonCuts(6); 
    lostMuonCuts[0] = deltaTCutLows[iter];
    lostMuonCuts[1] = deltaTCutHighs[iter];
    lostMuonCuts[2] = energyCutLows[iter];
    lostMuonCuts[3] = energyCutHighs[iter];
    lostMuonCuts[4] = clusterSizeCuts[iter];
    lostMuonCuts[5] = clusterEFracCuts[iter];
    lostMuonCuts.Write("LostMuonCuts");

    auto tripleDir = lostMuonsDir->mkdir("Triples");
    tripleDir->cd();

    DeltaT_triple[iter] = new TH1F("DeltaT_triple", "DeltaT_triple; time [ns]; Events", 1000, 0, 20);
    Edep_triple[iter] = new TH1F("Edep_triple", "Edep_triple; Energy [MeV]; Events", 1000, 0, 1000);
    Losses_triple[iter] = new TH1F("Losses_triple", "Losses_triple; time [ns]; Events", nBins, histMinTime, histMaxTime);

/////////////////////////////////////////////////////////////////////////////////////

} // end loop over hist iterations for defining histograms

/////////////////////////////////////////////////////////////////////////////////////

    double time_Get;
    double energy_Get;
    double x_Get;
    double y_Get;
    unsigned int nHit_Get;
    unsigned int caloNum_Get;
    unsigned int islandNum_Get;
    unsigned int event_Get;
    unsigned int bunchNum_Get;
    unsigned int midasSerialNum_Get;
    unsigned int subRun_Get;
    unsigned int run_Get;

    vector<unsigned int>* xtalNums_Get = 0;
    vector<double>* xtalEnergies_Get = 0;
    vector<double>* xtalTimes_Get = 0;
    vector<double>* xtalIFGAmplitudes_Get = 0;
    vector<double>* xtalIFGLifetimes_Get = 0;


  TTree* inputTree = (TTree*)inputFile->Get("clusterTree/clusters");
  Long64_t treeEntries = inputTree->GetEntries();

  inputTree->SetBranchAddress("time", &time_Get);
  inputTree->SetBranchAddress("energy", &energy_Get);
  inputTree->SetBranchAddress("x", &x_Get);
  inputTree->SetBranchAddress("y", &y_Get);
  inputTree->SetBranchAddress("nHit", &nHit_Get);
  inputTree->SetBranchAddress("caloNum", &caloNum_Get);
  inputTree->SetBranchAddress("islandNum", &islandNum_Get);
  inputTree->SetBranchAddress("eventNum", &event_Get);
  inputTree->SetBranchAddress("bunchNum", &bunchNum_Get);
  inputTree->SetBranchAddress("midasSerialNum", &midasSerialNum_Get);
  inputTree->SetBranchAddress("subRunNum", &subRun_Get);
  inputTree->SetBranchAddress("runNum", &run_Get);

  inputTree->SetBranchAddress("xtalNums", &xtalNums_Get);
  inputTree->SetBranchAddress("xtalEnergies", &xtalEnergies_Get);
  inputTree->SetBranchAddress("xtalTimes", &xtalTimes_Get);
  inputTree->SetBranchAddress("xtalIFGAmplitudes", &xtalIFGAmplitudes_Get);
  inputTree->SetBranchAddress("xtalIFGLifetimes", &xtalIFGLifetimes_Get);

/////////////////////////////////////////////////////////////////////////////////////

  // define randomizers above tree loop for later use
  TRandom3* randomizer_Time = new TRandom3(1); // the seeds don't matter, they are reset down below
  TRandom3* randomizer_UV = new TRandom3(1);

/////////////////////////////////////////////////////////////////////////////////////

  // loop through tree, fill vectors of hit info, and process them to fill pileup and main histograms

  vector<tuple<double, uint, double, double, double, double, uint, uint, vector<uint>, vector<double>, vector<double>, vector<double>>> hitsPerFillPerCalo;
  const int tuple_fill_index = 0, tuple_calo_index = 1, tuple_t_index = 2, tuple_E_index = 3, tuple_x_index = 4, tuple_y_index = 5, tuple_nHits_index = 6, tuple_bunch_index = 7, tuple_xtalNums_index = 8, tuple_xtalEnergies_index = 9, tuple_xtalGainAmps_index = 10, tuple_xtalGainTaus_index = 11;
  // 0 - fill ID, 1 - caloNum, 2 - time, 3 - energy, 4 - x, 5 - y, 6 - nHit, 7 - bunchNum, 8 - xtal nums, 9 - crystal hit energies, 10 - crystal gain amplitudes, 11 - crystal gain lifetimes

  set<double> fillSet; // set of fills - used for counting the number of fills
  int fillCaloNumber = 0;

  Long64_t treeEntry=0;
  while(treeEntry < treeEntries){

    if(fillCaloNumber % 1000 == 0) cout << "Processing fill/calo in file: " << fillCaloNumber << endl;
    fillCaloNumber++;

    hitsPerFillPerCalo.clear();

    // get first hit info of fill
    inputTree->GetEntry(treeEntry);
    double fill_ID = run_Get * 1e6 + subRun_Get*1e3 + event_Get;
    double thisFillID = fill_ID;
      fillSet.insert(fill_ID);
    uint caloNum = caloNum_Get, bunchNumber = bunchNum_Get, runNum = run_Get;

    // if(caloNum == 1) cout << "Proccessing FillID: " << std::setprecision(12) << fill_ID << endl; // << " " << caloNum << endl;

    // set the dataset (for time randomizing by dataset dependent VW period and setting correct ad hoc gain parameters)
    if(treeEntry == 0){
      int datasetCase = -1;

      totalInfoDir->cd();
      TNamed tag("datasetTag", ""); // write a tag to the file which says which dataset we're analyzing - the tag is contained within the Title of a TNamed object

      if(runNum > lowRunBound_60h && runNum < highRunBound_60h)                { datasetCase = 1; tag.SetTitle(datasetTag_60h.c_str()); } // 60h 
      else if(runNum > lowRunBound_9d && runNum < highRunBound_9d)             { datasetCase = 2; tag.SetTitle(datasetTag_9d.c_str()); } // 9d
      else if(runNum > lowRunBound_Endgame && runNum < highRunBound_Endgame)   { datasetCase = 3; tag.SetTitle(datasetTag_Endgame.c_str()); } // Endgame
      else if(runNum > lowRunBound_HighKick && runNum < highRunBound_HighKick) { datasetCase = 4; tag.SetTitle(datasetTag_HighKick.c_str()); } // HighKick
      else if(runNum > lowRunBound_Run2_C && runNum < highRunBound_Run2_C)     { datasetCase = 5; tag.SetTitle(datasetTag_Run2_C.c_str()); } // Run 2 C

      setDataset(datasetCase); // , trackerStationModel);
      tag.Write();
    }

/////////////////////////////////////////////////////////////////////////////////////

    // loop to fill vector of hits in fill/calo
    while(fill_ID == thisFillID && caloNum_Get == caloNum){
      double time = time_Get * 1.25; // clock ticks to ns

      if(time > timeCutLow && time < timeCutHigh && runNum > runCutLow && runNum < runCutHigh){
        hitsPerFillPerCalo.push_back(make_tuple(fill_ID, caloNum_Get, time, energy_Get, x_Get, y_Get, nHit_Get, bunchNum_Get, *xtalNums_Get, *xtalEnergies_Get, *xtalIFGAmplitudes_Get, *xtalIFGLifetimes_Get));
        runHist->Fill(runNum);
        subRunHist->Fill(subRun_Get);
      }

      // increment tree entry before jumping back to top of loop over hits in the fill/calo - if while condition fails then the tree entry will be properly set for the next fill/calo making sure to skip no hits
      treeEntry++;
      if(treeEntry >= treeEntries) break; // make sure I don't try to gather hits past the end of the tree for the last fill
      inputTree->GetEntry(treeEntry);
      fill_ID = runNum * 1e6 + subRun_Get*1e3 + event_Get;

    } // end of loop over fill/calo entries


/////////////////////////////////////////////////////////////////////////////////////
// Loop over iterations, make adjustements to hits, and fill histograms
/////////////////////////////////////////////////////////////////////////////////////

   // double crystalRow = getCrystalRow(get<tuple_nHits_index>(thisHit), get<tuple_xtalNums_index>(thisHit), get<tuple_xtalEnergies_index>(thisHit));

   for (uint iter = 0; iter < totalIters; ++iter) // iterations loop to scan over various parameters, energy threshold, pileup dead time, etc.
   {

    if(bunchNumScan && bunchNumber+1 != iter && iter != 0) continue;

    uint runGroup = getRunGroup(runNum);
    if(separateByRunGroup && runGroup+1 != iter && iter != 0) continue;

    // else if(crystalRowScan && thirdElement != iter && iter != 0) continue; // not sure how crystal row is supposed to work with my pileup method at the moment - need to think more about how iterations or clusters should be skipped

    // define iteration dependent quantities up here

    double halfPeriodGuess = gm2PeriodGuesses[iter]/2.;
    double totalChance = exp(halfPeriodGuess/weightingLifetimes[iter]) + exp(-halfPeriodGuess/weightingLifetimes[iter]) + 2.;
    double percentChanceUPlus = exp(halfPeriodGuess/weightingLifetimes[iter]) / totalChance;
    double percentChanceUMinus = exp(-halfPeriodGuess/weightingLifetimes[iter]) / totalChance;

/////////////////////////////////////////////////////////////////////////////////////

 // create adjusted hits vector to apply ADT or gain changes to on a per iteration basis

    auto adjustedHits = hitsPerFillPerCalo;
    if(adjustedHits.size() == 0) break;

    sort(adjustedHits.begin(), adjustedHits.end(), sortTuple); // in the Endgame dataset there is at least one time pair that is out of order, sort it here

/////////////////////////////////////////////////////////////////////////////////////

  // apply gain changes to crystals here

    if(recorrectGainAtCrystals){
      for(auto &hit : adjustedHits){

          double changedClusterEnergy = 0;
          
          for (uint crystalNum = 0; crystalNum < get<tuple_nHits_index>(hit); ++crystalNum)
          {
            double crystalTime = (get<tuple_t_index>(hit)/1.25) / 1000.; // units of clock ticks * 1000, convert to same units as the tau parameter (taken from the cluster time)
            double gainCorr = 1. / (1. - get<tuple_xtalGainAmps_index>(hit)[crystalNum] * exp(-crystalTime/get<tuple_xtalGainTaus_index>(hit)[crystalNum]));
            double newGainCorr = 1. / (1. - crystalGainAmpFactors[iter] * get<tuple_xtalGainAmps_index>(hit)[crystalNum] * exp(-crystalTime/(crystalGainTauFactors[iter] * get<tuple_xtalGainTaus_index>(hit)[crystalNum])));

            changedClusterEnergy += (get<tuple_xtalEnergies_index>(hit)[crystalNum]/gainCorr) * newGainCorr;
          }

          get<tuple_E_index>(hit) = changedClusterEnergy;
      } // loop over hits
    }

/////////////////////////////////////////////////////////////////////////////////////

  // apply global ad hoc gain correction here 

    if(applyAdHocGainCorrection){
      for(auto &hit : adjustedHits){
        get<tuple_E_index>(hit) = get<tuple_E_index>(hit) * (1 + adHocGainFunctions[iter]->Eval(get<tuple_t_index>(hit)));
      }
    }

/////////////////////////////////////////////////////////////////////////////////////

  // apply artificial dead time here - use same shadow pileup code but with a shadow gap time of 0 ns, and then changed adjustedHits vector

    if(artificialDeadTimes[iter] > 0){
        vector<pair<double, double> > preADT_hitPairs; // need times and energies of hits in this structure for the pileup code
        for(auto &hit : adjustedHits) preADT_hitPairs.emplace_back(get<tuple_t_index>(hit), get<tuple_E_index>(hit));

        uint indexA = 0;
        while (indexA < adjustedHits.size()-1)
        {
          vector<int> artificialShadowIndices, unneededTripleIndices;

          uint indSave = indexA;
          checkForShadowPulses(indexA, &preADT_hitPairs, artificialShadowIndices, unneededTripleIndices, artificialDeadTimes[iter], 0, 0, unseenPileupThreshold); // put in a 0 for the shadow window time when implementing ADT, put in a 0 for the ADT if you want to turn this off, put in a 0 for the second shadow deadtime

          if(artificialShadowIndices.size() > 0)
          {
            vector<pair<double, double> > pileupSinglets;
            pileupSinglets.emplace_back(get<tuple_t_index>(adjustedHits.at(indSave)), get<tuple_E_index>(adjustedHits.at(indSave)));

            for (int i = artificialShadowIndices.size()-1; i >= 0; --i) // combines all pulses in ADT window into a single pulse
            {
              int shadowIndex = artificialShadowIndices.at(i);
              pileupSinglets.emplace_back(get<tuple_t_index>(adjustedHits.at(shadowIndex)), get<tuple_E_index>(adjustedHits.at(shadowIndex)));

              adjustedHits.erase(adjustedHits.begin()+shadowIndex);
              preADT_hitPairs.erase(preADT_hitPairs.begin()+shadowIndex); // need to erase here too since preADT_hitPairs is defined only once above the loop
            } 

            double addedEnergy = pileupAddedEnergy(pileupSinglets, pileupEnergyScales[iter]);
            double weightedTime = energyWeightedTime(pileupSinglets, 0);

            auto tempTuple = adjustedHits.at(indSave);
            adjustedHits.erase(adjustedHits.begin()+indSave);
              get<tuple_t_index>(tempTuple) = weightedTime;
              get<tuple_E_index>(tempTuple) = addedEnergy;
            adjustedHits.insert(adjustedHits.begin()+indSave, tempTuple);

            preADT_hitPairs.erase(preADT_hitPairs.begin()+indSave); // need to erase here too since preADT_hitPairs is defined only once above the loop
            preADT_hitPairs.insert(preADT_hitPairs.begin()+indSave, make_pair(weightedTime, addedEnergy));

            indexA--; // turn this on or off depending on whether I want to re-check for pulses in the ADT of a pulse that has already been artificially combined
          } // end if artificial shadow indices size > 0

        } // end while
    } // end if ADT > 0

/////////////////////////////////////////////////////////////////////////////////////

  // fill histograms of delta T between hits (before I apply any time randomization and after the ADT)
           
      for (uint hit = 0; hit < adjustedHits.size(); ++hit)
      {    
        if(hit != adjustedHits.size()-1) // don't do this for the very last hit
        {
          double deltaT = get<tuple_t_index>(adjustedHits.at(hit+1)) - get<tuple_t_index>(adjustedHits.at(hit));
          deltaT_all[iter]->Fill(deltaT);
          calo_deltaT[iter][caloNum-1]->Fill(deltaT);
        }
      }

/////////////////////////////////////////////////////////////////////////////////////

   // construct pileup from shadow pulses in data in the same fill/calo

      vector<pair<double, double> > shadowNegativePileup, shadowPositivePileup; // energies and times of constructed pileup pulses which contribute either negatively or positively
      vector<uint> triggerIndicesNegative, triggerIndicesPositive; // trigger indices of corresponding additions to the constructed negative and positive pileup counts

      vector<pair<double, double> > hitPairs; // defined in the same way as was done in the ADT loop, but done again to include changes from said ADT loop
      for(auto &hit : adjustedHits) hitPairs.emplace_back(get<tuple_t_index>(hit), get<tuple_E_index>(hit));

      uint j = 0;
      while (j < adjustedHits.size()-1) // j is incremented in checkForShadowPulses()
      {
        std::vector<int> shadowIndices, tripleShadowIndices;

        uint indSave = j;
        checkForShadowPulses(j, &hitPairs, shadowIndices, tripleShadowIndices, shadowDeadTimes[iter], shadowGapTimes[iter], shadowSecondDeadTimes[iter], unseenPileupThreshold);

        if(shadowIndices.size() > 0)
        {
          vector<pair<double, double> > pileupSinglets;
          pileupSinglets.emplace_back(get<tuple_t_index>(adjustedHits.at(indSave)), get<tuple_E_index>(adjustedHits.at(indSave)));

          for (uint i = 0; i < shadowIndices.size(); ++i)
          {
            int shadowIndex = shadowIndices.at(i);
            pileupSinglets.emplace_back(get<tuple_t_index>(adjustedHits.at(shadowIndex)), get<tuple_E_index>(adjustedHits.at(shadowIndex)));
          }

          double addedEnergy = pileupAddedEnergy(pileupSinglets, pileupEnergyScales[iter]);
          double weightedTime = energyWeightedTime(pileupSinglets, shadowGapTimes[iter]);

          for (uint i = 0; i < pileupSinglets.size(); ++i){
            shadowNegativePileup.emplace_back(weightedTime, pileupSinglets.at(i).second);
            triggerIndicesNegative.push_back(indSave);
          } 
          shadowPositivePileup.emplace_back(weightedTime, addedEnergy);
          triggerIndicesPositive.push_back(indSave);

          // fill error histograms - currently only accounts for errors for doublets (0 and 1 entries in pileupSinglets)
          if(addedEnergy > lowerEnergyThresholds[iter] && addedEnergy < upperEnergyThresholds[iter])
          {
                 if((pileupSinglets.at(0).second > lowerEnergyThresholds[iter] && pileupSinglets.at(0).second < upperEnergyThresholds[iter]) && (pileupSinglets.at(1).second > lowerEnergyThresholds[iter] && pileupSinglets.at(1).second < upperEnergyThresholds[iter])) pileupErrorsBoth_shadow[iter][caloNum-1]->Fill(weightedTime); // both singlets are above threshold (or in the threshold range)
            else if((pileupSinglets.at(0).second < lowerEnergyThresholds[iter] || pileupSinglets.at(0).second > upperEnergyThresholds[iter]) && (pileupSinglets.at(1).second < lowerEnergyThresholds[iter] || pileupSinglets.at(1).second > upperEnergyThresholds[iter])) pileupErrorsNeither_shadow[iter][caloNum-1]->Fill(weightedTime); // nether singlet is above threshold (or in the threshold range)
          }

          // Triplets
          // Experimental section - doesn't quite work. Not sure how shadow dead times feed into what the combinatorics or amplitudes should be.

          if(performTripletCorrection && tripleShadowIndices.size() > 0){

            // cout << "Triplet pileup correction doesn't quite work right." << endl;
            // return 1;

            double pulse1_E = pileupSinglets.at(0).second;
            double pulse2_E = pileupSinglets.at(1).second; // energy of pulse(s) in first shadow window - only the first for now since there's only 1 pulse if SDT = ADT
            double pulse3_E = 0; // energy of pulse(s) in second shadow window

            for (uint i = 0; i < tripleShadowIndices.size(); ++i) { pulse3_E += get<tuple_E_index>(adjustedHits.at(tripleShadowIndices.at(i))); } // Add up all the energy in the second shadow window

            double pulse12_E = pulse1_E + pulse2_E;
            double pulse13_E = pulse1_E + pulse3_E;
            double pulse23_E = pulse2_E + pulse3_E;

            double pulse123_E = pulse1_E + pulse2_E + pulse3_E;

        // new tests

            shadowNegativePileup.emplace_back(weightedTime, pulse123_E); triggerIndicesNegative.push_back(indSave); // make sure to push corresponding trigger indices back everytime I add to the positive or negative pileup spectra
            shadowPositivePileup.emplace_back(weightedTime, pulse1_E); triggerIndicesPositive.push_back(indSave);
            shadowPositivePileup.emplace_back(weightedTime, pulse2_E); triggerIndicesPositive.push_back(indSave);
            shadowPositivePileup.emplace_back(weightedTime, pulse3_E); triggerIndicesPositive.push_back(indSave);

            shadowPositivePileup.emplace_back(weightedTime, pulse12_E); triggerIndicesPositive.push_back(indSave);
            shadowNegativePileup.emplace_back(weightedTime, pulse1_E); triggerIndicesNegative.push_back(indSave);
            shadowNegativePileup.emplace_back(weightedTime, pulse2_E); triggerIndicesNegative.push_back(indSave);
            shadowPositivePileup.emplace_back(weightedTime, pulse23_E); triggerIndicesPositive.push_back(indSave);
            shadowNegativePileup.emplace_back(weightedTime, pulse2_E); triggerIndicesNegative.push_back(indSave);
            shadowNegativePileup.emplace_back(weightedTime, pulse3_E); triggerIndicesNegative.push_back(indSave);

        // old code trying different combinations:

            // for(int n = 0; n < 3; n++) shadowNegativePileup.emplace_back(weightedTime, pulse1_E + pulse2_E + pulse3_E);
            // for(int n = 0; n < 3; n++) shadowPositivePileup.emplace_back(weightedTime, pulse1_E + pulse2_E);
            // for(int n = 0; n < 4; n++) shadowPositivePileup.emplace_back(weightedTime, pulse2_E + pulse3_E);
            // for(int n = 0; n < 1; n++) shadowNegativePileup.emplace_back(weightedTime, pulse1_E);
            // for(int n = 0; n < 2; n++) shadowNegativePileup.emplace_back(weightedTime, pulse2_E);
            // for(int n = 0; n < 2; n++) shadowNegativePileup.emplace_back(weightedTime, pulse3_E);

          } // end triplets stuff

        } // if shadow indices size > 0

      }  //end while

/////////////////////////////////////////////////////////////////////////////////////

  // apply randomization here, after gain corrections and pileup construction, but before filling pileup histograms - randomizers defined above the tree loop

    vector<double> randomizationTimeShifts (adjustedHits.size(), 0); // initialize all to 0 with the same size as the hits vector

    double T_fc = binWidths[iter]; // defaultBinWidth - if the bin width is changed to something significantly different than the FR period, then the randomization will not use the correct period unless forced to
    double T_VW = getTVW(0, T_fc); // not set up to work with time dependent randomization currently

    if(!varyRandSeeds && perFillRandomization) randomizer_Time->SetSeed(get<tuple_fill_index>(adjustedHits.at(0))); // set seed of randomizer using the fill ID
    else if(!varyRandSeeds && !perFillRandomization) randomizer_Time->SetSeed(get<tuple_fill_index>(adjustedHits.at(0)) * get<tuple_calo_index>(adjustedHits.at(0))); // set seed of randomizer using the fill ID and calo number
    else if(varyRandSeeds && perFillRandomization) randomizer_Time->SetSeed((iter+1) * get<tuple_fill_index>(adjustedHits.at(0))); // if varying the random seeds then multiply the fill id by the iteration number
    else if(varyRandSeeds && !perFillRandomization) randomizer_Time->SetSeed((iter+1) * get<tuple_fill_index>(adjustedHits.at(0)) * get<tuple_calo_index>(adjustedHits.at(0)));
    // else if(varyRandSeeds && !perFillRandomization) randomizer_Time->SetSeed((iter+1+100) * get<tuple_fill_index>(adjustedHits.at(0)) * get<tuple_calo_index>(adjustedHits.at(0))); // this for doing the next set of 100 seeds

    double timeRand;

    if(randomize_fc){
      if(perFillRandomization){ // per fill time randomization (same for all clusters in all calos in a fill)
        timeRand = T_fc * (randomizer_Time->Uniform() - 0.5);
        if(randomize_VW) timeRand += T_VW * (randomizer_Time->Uniform() - 0.5);
        else randomizer_Time->Uniform(); // still generate the random number even if I don't use it, so that the cluster f_c randomization is identical with and without the VW randomization (other option is to use a different randomizer for the VW, but then I'd have to regenerate the histograms I've already made)
        fill(randomizationTimeShifts.begin(), randomizationTimeShifts.end(), timeRand); // fill entire hits vector with same randomization
      }
      else{ // per cluster time randomization (different for every single cluster)
        for (uint i = 0; i < adjustedHits.size(); ++i){
          timeRand = T_fc * (randomizer_Time->Uniform() - 0.5);
          if(randomize_VW) timeRand += T_VW * (randomizer_Time->Uniform() - 0.5);
          else randomizer_Time->Uniform(); // still generate the random number even if I don't use it, so that the cluster f_c randomization is identical with and without the VW randomization (other option is to use a different randomizer for the VW, but then I'd have to regenerate the histograms I've already made)
          randomizationTimeShifts.at(i) = timeRand; // fill hits vector with different randomization for each cluster
        }
      }
    }

    // apply time randomization

    // cout << "Positive pileup size: " << shadowPositivePileup.size() << " positive trigger indices size: " << triggerIndicesPositive.size() << " negative pileup size: " << shadowNegativePileup.size() << " negative trigger indices size: " << triggerIndicesNegative.size() << endl;

    for(uint hit = 0; hit < adjustedHits.size(); hit++){
      get<tuple_t_index>(adjustedHits.at(hit)) = get<tuple_t_index>(adjustedHits.at(hit)) + randomizationTimeShifts.at(hit);
      timeShifts[iter]->Fill(randomizationTimeShifts.at(hit));
    }

    if(perFillRandomization){
      for (uint i = 0; i < shadowPositivePileup.size(); ++i)  shadowPositivePileup.at(i).first = shadowPositivePileup.at(i).first + randomizationTimeShifts.at(0);
      for (uint i = 0; i < shadowNegativePileup.size(); ++i)  shadowNegativePileup.at(i).first = shadowNegativePileup.at(i).first + randomizationTimeShifts.at(0);
    }
    else{
      for (uint i = 0; i < shadowPositivePileup.size(); ++i)  shadowPositivePileup.at(i).first = shadowPositivePileup.at(i).first + randomizationTimeShifts.at(triggerIndicesPositive.at(i));
      for (uint i = 0; i < shadowNegativePileup.size(); ++i)  shadowNegativePileup.at(i).first = shadowNegativePileup.at(i).first + randomizationTimeShifts.at(triggerIndicesNegative.at(i));
    }

    // do UV randomization - different random numbers for each cluster which are used later on when filling the histograms

    vector<double> randomizationUV;

    if(!varyRandSeeds) randomizer_UV->SetSeed(get<tuple_fill_index>(adjustedHits.at(0)) * get<tuple_calo_index>(adjustedHits.at(0)) * 123); // set seed of randomizer using the fill ID (but with a different number than the time randomizer, and a different seed for each calo (so UV is not the same from calo to calo))
    else randomizer_UV->SetSeed((iter+1) * get<tuple_fill_index>(adjustedHits.at(0)) * get<tuple_calo_index>(adjustedHits.at(0)) * 123);
    // else randomizer_UV->SetSeed((iter+1+100) * get<tuple_fill_index>(adjustedHits.at(0)) * get<tuple_calo_index>(adjustedHits.at(0)) * 123); // this for doing the next set of random seeds

    for (uint i = 0; i < adjustedHits.size(); ++i){
      double UV_rand = randomizer_UV->Uniform();
      randomizationUV.push_back(UV_rand);
      UVrands[iter]->Fill(UV_rand);
    }

/////////////////////////////////////////////////////////////////////////////////////

    // fill pileup histograms

          double pTime, pEnergy;

          for (uint i = 0; i < shadowNegativePileup.size(); ++i)
          {
            pTime = shadowNegativePileup.at(i).first + pileupTimeShifts[iter];
            pEnergy = shadowNegativePileup.at(i).second;

            pileupTimes_shadow[iter][caloNum-1]->Fill(pTime, -1);
            pileupEnergies_shadow[iter][caloNum-1]->Fill(pEnergy, -1);

                if(pTime < 250000) pileupEnergies_shadow_early[iter][caloNum-1]->Fill(pEnergy, -1);
                else if(pTime > 250000) pileupEnergies_shadow_late[iter][caloNum-1]->Fill(pEnergy, -1);

              if(pEnergy > lowerEnergyThresholds[iter] && pEnergy < upperEnergyThresholds[iter])
              {
                pileupTimes_shadow_threshold[iter][caloNum-1]->Fill(pTime, -1);
                pileupEnergiesThreshold_shadow[iter][caloNum-1]->Fill(pEnergy, -1);

                double randNum_pileup = randomizationUV.at(triggerIndicesNegative.at(i)); // grab the randomization with the pileup pulse trigger indices
                  if     (randNum_pileup < percentChanceUPlus) pileupTimes_U_shadow[iter][caloNum-1]->Fill(pTime - halfPeriodGuess, -1);
                  else if(randNum_pileup < (percentChanceUPlus + percentChanceUMinus)) pileupTimes_U_shadow[iter][caloNum-1]->Fill(pTime + halfPeriodGuess, -1);
                  else if(randNum_pileup < 1) pileupTimes_V_shadow[iter][caloNum-1]->Fill(pTime, -1);
              }

          }

          for (uint i = 0; i < shadowPositivePileup.size(); ++i)
          {            
            pTime = shadowPositivePileup.at(i).first + pileupTimeShifts[iter];
            pEnergy = shadowPositivePileup.at(i).second;

            pileupTimes_shadow[iter][caloNum-1]->Fill(pTime);
            pileupEnergies_shadow[iter][caloNum-1]->Fill(pEnergy);

                if(pTime < 250000) pileupEnergies_shadow_early[iter][caloNum-1]->Fill(pEnergy);
                if(pTime > 250000) pileupEnergies_shadow_late[iter][caloNum-1]->Fill(pEnergy);

              if(pEnergy > lowerEnergyThresholds[iter] && pEnergy < upperEnergyThresholds[iter])
              {
                pileupTimes_shadow_threshold[iter][caloNum-1]->Fill(pTime);
                pileupEnergiesThreshold_shadow[iter][caloNum-1]->Fill(pEnergy, -1);

                double randNum_pileup = randomizationUV.at(triggerIndicesPositive.at(i));
                  if     (randNum_pileup < percentChanceUPlus) pileupTimes_U_shadow[iter][caloNum-1]->Fill(pTime - halfPeriodGuess);
                  else if(randNum_pileup < (percentChanceUPlus + percentChanceUMinus)) pileupTimes_U_shadow[iter][caloNum-1]->Fill(pTime + halfPeriodGuess);
                  else if(randNum_pileup < 1) pileupTimes_V_shadow[iter][caloNum-1]->Fill(pTime);
              }
          }


/////////////////////////////////////////////////////////////////////////////////////
// Now fill main histograms (still under iterations loop)
/////////////////////////////////////////////////////////////////////////////////////

     for (uint thisHit = 0; thisHit < adjustedHits.size(); thisHit++) // loop over entries in that fill/calo
     {
       double time = get<tuple_t_index>(adjustedHits.at(thisHit));
       double energy = get<tuple_E_index>(adjustedHits.at(thisHit));

       // double crystalRow = getCrystalRow(get<tuple_nHits_index>(thisHit), get<tuple_xtalNums_index>(thisHit), get<tuple_xtalEnergies_index>(thisHit));
       // else if(crystalRowScan && crystalRow != iter && iter != 0) break;

          bunchNumHist[iter]->Fill(bunchNumber);
          if(separateByRunGroup) runGroupHist[iter]->Fill(runGroup);
          runNumHist_perIter[iter]->Fill(runNum);

          // crystalRowHist[iter]->Fill(crystalRow);

           allEnergies[iter]->Fill(energy);
           allTimes[iter]->Fill(time);

           caloEnergies[iter][caloNum-1]->Fill(energy);
           caloTimes[iter][caloNum-1]->Fill(time);

              if(time < 250000){
                allEnergies_early[iter]->Fill(energy);
                caloEnergies_early[iter][caloNum-1]->Fill(energy);
              } 
              else if(time > 250000){
                allEnergies_late[iter]->Fill(energy);
                caloEnergies_late[iter][caloNum-1]->Fill(energy);
              } 

           if(!reduceMemory){
             timesAndEnergies[iter]->Fill(time, energy);
             caloTimesAndEnergies[iter][caloNum-1]->Fill(time, energy);
           }

           if(energy > lowerEnergyThresholds[iter] && energy < upperEnergyThresholds[iter])
           {
            thresholdEnergies[iter]->Fill(energy);
            thresholdTimes[iter]->Fill(time);
            caloEnergiesThreshold[iter][caloNum-1]->Fill(energy);
            caloTimesThreshold[iter][caloNum-1]->Fill(time);

            // fill ratio histograms

            double randNum = randomizationUV.at(thisHit);

              if(randNum < percentChanceUPlus){
                 allTimesAdded_U[iter]->Fill(time - halfPeriodGuess); // careful with the signs here - the U+ hist is filled with pulses shifted by t - T/2, and weighted by e^+T/2tau, and vice versa for U-
                 caloTimesThresholdU[iter][caloNum-1]->Fill(time - halfPeriodGuess);
              }
              else if(randNum < (percentChanceUPlus + percentChanceUMinus)){
                 allTimesAdded_U[iter]->Fill(time + halfPeriodGuess);
                 caloTimesThresholdU[iter][caloNum-1]->Fill(time + halfPeriodGuess);
              }
              else if(randNum < 1){
                 allTimesAdded_V[iter]->Fill(time);
                 caloTimesThresholdV[iter][caloNum-1]->Fill(time);
              }

           } // end if above threshold

         } // end loop over vector of hits for main histogram filling

/////////////////////////////////////////////////////////////////////////////////////

     } // end loop over total iterations 

  } // end of loop over tree entries

  for (auto &i : fillSet) numFills->Fill(0); // count number of fills at end after loop through tree


/////////////////////////////////////////////////////////////////////////////////////

  // fill lost muons histograms

  unsigned int coincidenceLevel_Get;
  vector<int>* caloNum_lostMuons_Get = 0;
  vector<double>* clusterEnergy_Get = 0;
  vector<double>* clusterTime_Get = 0;
  vector<double>* clusterX_Get = 0;
  vector<double>* clusterY_Get = 0;
  vector<double>* clusterSize_Get = 0;
  vector<double>* clusterEfrac_Get = 0;

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

    for(Long64_t i=0; i<lostMuonEntries; i++) {
       lostMuonTree->GetEntry(i);

         if(coincidenceLevel_Get == 3)
         {
          double time1 = clusterTime_Get->at(0) * 1.25;
          double time2 = clusterTime_Get->at(1) * 1.25;
          double time3 = clusterTime_Get->at(2) * 1.25;
          double deltaT1 = time2 - time1;
          double deltaT2 = time3 - time2;

          double energy1 = clusterEnergy_Get->at(0);
          double energy2 = clusterEnergy_Get->at(1);
          double energy3 = clusterEnergy_Get->at(2);

          double size1 = clusterSize_Get->at(0);
          double size2 = clusterSize_Get->at(1);
          double size3 = clusterSize_Get->at(2);

          double eFrac1 = clusterEfrac_Get->at(0);
          double eFrac2 = clusterEfrac_Get->at(1);
          double eFrac3 = clusterEfrac_Get->at(2);

          for (uint iter = 0; iter < totalIters; ++iter)
          {
            if(time1 > timeCutLowMuons && time2 > timeCutLowMuons && time3 > timeCutLowMuons && size1 < clusterSizeCuts[iter] && size2 < clusterSizeCuts[iter] && size3 < clusterSizeCuts[iter] && eFrac1 > clusterEFracCuts[iter] && eFrac2 > clusterEFracCuts[iter] && eFrac3 > clusterEFracCuts[iter])
            {
              DeltaT_triple[iter]->Fill(deltaT1);
              DeltaT_triple[iter]->Fill(deltaT2);
              Edep_triple[iter]->Fill(energy1);
              Edep_triple[iter]->Fill(energy2);
              Edep_triple[iter]->Fill(energy3);

              if(deltaT1 > deltaTCutLows[iter] && deltaT1 < deltaTCutHighs[iter] && deltaT2 > deltaTCutLows[iter] && deltaT2 < deltaTCutHighs[iter]  && energy1 > energyCutLows[iter] && energy1 < energyCutHighs[iter] && energy2 > energyCutLows[iter] && energy2 < energyCutHighs[iter] && energy3 > energyCutLows[iter] && energy3 < energyCutHighs[iter])
              {
                Losses_triple[iter]->Fill(time1);
              }
             }
           } // end loop over iterations
         } // end triples
     } // end for loop over lost muons tree

/////////////////////////////////////////////////////////////////////////////////////

    t = clock() - t;
    printf ("It took me %f seconds.\n", ((float)t)/CLOCKS_PER_SEC);

/////////////////////////////////////////////////////////////////////////////////////

   outputFile->Write();
   delete outputFile;

/////////////////////////////////////////////////////////////////////////////////////

   return 0;

}
