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
#include <TVirtualFFT.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TVectorD.h>
#include <THStack.h>
#include <TPaveStats.h>
#include <TVirtualHistPainter.h>

#include <time.h>

#include "ratioAnalysisDefs.hh"
#include "ratioToyHistConfig.hh"
#include "pileupUtils.hh"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////

double gainSagFunc(double* x, double* p)
{
  double time = x[0];

  double alpha = 4.2e-7; // some constant
  double n0 = 20.; // number of pulses per fill - 500 CTAGs / 24 calos
  double energy = 20000 * (p[0]); // not really sure about this energy constant parameter

  double tauDrop = 1000.; // intitial gain drop time constant in ns
  double tauRelax = p[1]; // defaultLifetime; // this was tau mu for what Nandita had in her gain paper, but most recently with the mega boxes I think this is something different

  // double gainValue = 1 - alpha * energy * n0 * (tauDrop/(tauRelax-tauDrop)) * (exp(-time/tauRelax) - exp(-time/tauDrop));
  double gainValue = 1 - alpha * energy * n0 * 0.2 * (exp(-time/tauRelax) - exp(-time/tauDrop)); // just put in a constant (.2) instead of the tau combination term

  return gainValue;
}

class gainCorrClass {

public:
  gainCorrClass(TF1* gainSagFunc){
    sagFunc = gainSagFunc;
  }

  double gainCorrFunc(double* x, double* par){
    return 1./sagFunc->Eval(x[0]);
  }

private:
  TF1* sagFunc;

}; // end gainCorrClass class


/////////////////////////////////////////////////////////////////////////////////////

// int makeToyHistFromTH2PileupMod(int fN) // for running multiple instances of this module - the seed has to be set to 0 for this to work
int makeToyHistFromTH2PileupMod()
{

  // check config before analyzing - this isn't configured to work with energy changes (local gain) in any way yet
  int configCheck = 0;
  if(gm2PeriodPPMStep) configCheck++;
  if(energyThresholdStep) configCheck++;
  if(SDTStep) configCheck++;
  if(SGTStep) configCheck++;
  if(SDT2Step) configCheck++;
  if(pileupMultiplierStep) configCheck++;
  if(pileupTimeShiftStep) configCheck++;

  if(configCheck > 1){
    cout << "Histogram making configuration set incorrectly." << endl;
    return -1;
  }


/////////////////////////////////////////////////////////////////////////////////////

  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  clock_t startTime = clock();

  gRandom->SetSeed(toyRandSeed1); // for TH2 GetRandom
  TRandom3* randFill = new TRandom3(toyRandSeed2); // for randomized entries and entries per fill
  TRandom3* tripRand = new TRandom3(toyRandSeed3); // for triplet randomization

  TFile *inputFile = TFile::Open(inputFilePath.c_str());
     if (inputFile == 0) {
        printf("Error: cannot open file\n");
        return 0;
     }

  TFile* outputFile = new TFile("histFromPairs.root","RECREATE");
  // TFile* outputFile = new TFile(Form("histFromPairs-threeIters-%i.root",fN),"RECREATE");

  TH2F* initialHist = (TH2F*) inputFile->Get("full2Dfunction_hist");

/////////////////////////////////////////////////////////////////////////////////////

  auto topDir = gFile->mkdir("topDir");

/////////////////////////////////////////////////////////////////////////////////////

  TVectorD numIters(1);
  numIters[0] = totalIters;
  numIters.Write("Iters");

  TDirectory* histIterDirs[totalIters];

/////////////////////////////////////////////////////////////////////////////////////

  int subIters = 4; // 4 for truth, measured, measured - true pileup, measured - shadow pileup

  TH1F* toyFiveParamHists[totalIters][subIters];
  TH1F* toyUHists[totalIters][subIters];
  TH1F* toyVHists[totalIters][subIters];
  TH1F* toyNumHists[totalIters][subIters];
  TH1F* toyDenomHists[totalIters][subIters];
  TH1F* energies[totalIters][subIters];

  TH1F* deltaT_hits[totalIters][subIters];
  TH1F* deltaT_hits_threshold[totalIters][subIters];

/////////////////////////////////////////////////////////////////////////////////////

  TH1F* energies_early[totalIters][subIters];
  TH1F* energies_late[totalIters][subIters];

  int numTimeSlices = 10;
  std::vector<double> timeSliceLowEdges;
  for (int timeSlice = 0; timeSlice < numTimeSlices; ++timeSlice) timeSliceLowEdges.push_back(timeSlice*20000);

  TH1F* energies_timeSlice[totalIters][subIters][numTimeSlices];
  TH1F* energies_g2phase[totalIters][subIters][8];

/////////////////////////////////////////////////////////////////////////////////////

  TF1* gainSagFunctions[totalIters];
  TF1* gainCorrFunctions[totalIters];

/////////////////////////////////////////////////////////////////////////////////////

  double gainSagAmplitudeMultiplier[totalIters]; 
  double gainSagLifetime[totalIters];

  double gm2PeriodPPMOffsets[totalIters];
  double gm2PeriodGuesses[totalIters];
  double energyThresholds[totalIters];

  double shadowDeadTimes[totalIters];
  double shadowGapTimes[totalIters];
  double shadowSecondDeadTimes[totalIters];

  double pileupScaleFactors[totalIters];
  double pileupTimeShifts[totalIters];

/////////////////////////////////////////////////////////////////////////////////////

    TRandom3* randIterUV[totalIters]; // for individual iteration UV randomization
    TRandom3* randUVPileup[totalIters]; // for UV pileup filling

    int iterRandSeed1 = randFill->Integer(1e8);
    int iterRandSeed2 = randFill->Integer(1e8);

/////////////////////////////////////////////////////////////////////////////////////

  for (int iter = 0; iter < totalIters; ++iter)
  {
    // gain sag stuff - will need to revisit this at some point
    gainSagAmplitudeMultiplier[iter] = gainSagAmpFactor;
    gainSagLifetime[iter] = gainSagTau;
    // gainSagAmplitudeMultiplier[iter] = 1.*(iter-(totalIters/2.))/totalIters; // scan through 0 gain sag (negative to positive effect)

/////////////////////////////////////////////////////////////////////////////////////

    gm2PeriodPPMOffsets[iter] = gm2PeriodPPMStart + iter * gm2PeriodPPMStep;
    gm2PeriodGuesses[iter] = g2Period * (1 + 1e-6 * gm2PeriodPPMOffsets[iter]);
    
    energyThresholds[iter] = energyThresholdStart + iter * energyThresholdStep;

    shadowDeadTimes[iter] = SDTStart + iter * SDTStep;
    shadowGapTimes[iter] = SGTStart + iter * SGTStep;
    shadowSecondDeadTimes[iter] = SDT2Start + iter * SDT2Step;

    pileupScaleFactors[iter] = pileupMultiplierStart + iter * pileupMultiplierStep;
    pileupTimeShifts[iter] = pileupTimeShiftStart + iter * pileupTimeShiftStep;

    randIterUV[iter] = new TRandom3(iterRandSeed1);
    randUVPileup[iter] = new TRandom3(iterRandSeed2); // for UV pileup filling
  }

/////////////////////////////////////////////////////////////////////////////////////

    TH1::SetDefaultSumw2();

    double dT_maxtime = 100000;

  for (int iter = 0; iter < totalIters; ++iter)
  {
    histIterDirs[iter] = topDir->mkdir(Form("Iter%d", iter));
    histIterDirs[iter]->cd();

    for (int subIter = 0; subIter < subIters; ++subIter)
    {
      auto subIterDir = histIterDirs[iter]->mkdir(Form("SubIter%d", subIter));
      subIterDir->cd();

        toyFiveParamHists[iter][subIter] = new TH1F("Toy_5_Param_Hist","Toy_5_Param_Hist; time (ns); Events",nBins,0,histMaxTime);
        toyUHists[iter][subIter] = new TH1F("Toy_U_Hist","Toy_U_Hist; time (ns); Events",nBins,0,histMaxTime);
        toyVHists[iter][subIter] = new TH1F("Toy_V_Hist","Toy_V_Hist; time (ns); Events",nBins,0,histMaxTime);
        toyNumHists[iter][subIter] = new TH1F("Toy_Num_Hist","Toy_Num_Hist; time (ns); Events",nBins,0,histMaxTime);
        toyDenomHists[iter][subIter] = new TH1F("Toy_Denom_Hist","Toy_Denom_Hist; time (ns); Events",nBins,0,histMaxTime);

        energies[iter][subIter] = new TH1F("energies","energies; Energy (MeV); Events",1000,0,10000);

        energies_early[iter][subIter] = new TH1F("energies_early","energies_early; Energy (MeV); Events",1000,0,10000);
        energies_late[iter][subIter] = new TH1F("energies_late","energies_late; Energy (MeV); Events",1000,0,10000);

        deltaT_hits[iter][subIter] = new TH1F("Time_Between_Hits_All", "Time_Between_Hits_All; #DeltaT (ns); Events", 10000, 0, dT_maxtime);
        deltaT_hits_threshold[iter][subIter] = new TH1F("Time_Between_Hits_Above_Threshold", "Time_Between_Hits_Above_Threshold; #DeltaT (ns); Events", 10000, 0, dT_maxtime);

      auto TS_dir = subIterDir->mkdir("TimeSlice");
      TS_dir->cd();
        for (int timeSlice = 0; timeSlice < numTimeSlices; ++timeSlice) energies_timeSlice[iter][subIter][timeSlice] = new TH1F(Form("energies_tS%i",timeSlice), Form("energies_tS%i; Energy (MeV); Events",timeSlice), 1000, 0, 10000);
      auto g2_dir = subIterDir->mkdir("g2");
      g2_dir->cd();
        for (int g2phase = 0; g2phase < 8; ++g2phase) energies_g2phase[iter][subIter][g2phase] = new TH1F(Form("energies_g2phase%i",g2phase), Form("energies_g2phase%i; Energy (MeV); Events",g2phase), 1000, 0, 10000);
    }

/////////////////////////////////////////////////////////////////////////////////////

    auto savDirs = histIterDirs[iter]->mkdir("SavedParameters");
    savDirs->cd();

    gainSagFunctions[iter] = new TF1("gainSagFunction", gainSagFunc, 0, 1000000, 2);
    gainSagFunctions[iter]->SetNpx(10000);
    gainSagFunctions[iter]->SetParameter(0, gainSagAmplitudeMultiplier[iter]);
    gainSagFunctions[iter]->SetParameter(1, gainSagLifetime[iter]);
    gainSagFunctions[iter]->Write();

    auto corrClass = new gainCorrClass(gainSagFunctions[iter]); // if grabbing gainSagFunctions[0] (the no gain sag curve) then there is no correction performed
    gainCorrFunctions[iter] = new TF1("gainCorrFunction", corrClass, &gainCorrClass::gainCorrFunc, 0, 1000000, 0);
    gainCorrFunctions[iter]->SetNpx(10000);
    gainCorrFunctions[iter]->Write();


    TVectorD parameterStore(8); // vector of doubles to store parameters used when creating histograms, for systematic studies

    parameterStore[0] = gm2PeriodGuesses[iter];
    parameterStore[1] = energyThresholds[iter];
    parameterStore[2] = shadowDeadTimes[iter];
    parameterStore[3] = shadowGapTimes[iter];
    parameterStore[4] = artificialDeadtime;
    parameterStore[5] = pileupScaleFactors[iter];
    parameterStore[6] = shadowSecondDeadTimes[iter];
    parameterStore[7] = pileupTimeShifts[iter];
    parameterStore.Write("parameterStore");
  }

/////////////////////////////////////////////////////////////////////////////////////


TH1F* pileupTimes[totalIters];
TH1F* pileupTimes_threshold[totalIters];
TH1F* pileupEnergies[totalIters];
TH1F* pileupTimes_U[totalIters];
TH1F* pileupTimes_V[totalIters];

TH1F* pileupTimes_shadow[totalIters];
TH1F* pileupTimes_shadow_threshold[totalIters];
TH1F* pileupEnergies_shadow[totalIters];
TH1F* pileupTimes_U_shadow[totalIters];
TH1F* pileupTimes_V_shadow[totalIters];

  TH1F* errorsNeither[totalIters];
  TH1F* errorsBoth[totalIters];

  TH1F* errorsNeither_shadow[totalIters];
  TH1F* errorsBoth_shadow[totalIters];

TDirectory* pileupDirs[totalIters];

/////////////////////////////////////////////////////////////////////////////////////

TH1F* pileupEnergies_shadow_early[totalIters];
TH1F* pileupEnergies_shadow_late[totalIters];

TH1F* pileupEnergies_shadow_timeSlice[totalIters][numTimeSlices];
TH1F* pileupEnergies_shadow_g2phase[totalIters][8];

/////////////////////////////////////////////////////////////////////////////////////

  for (int iter = 0; iter < totalIters; ++iter)
  {
    pileupDirs[iter] = histIterDirs[iter]->mkdir("PileupPlots");
    pileupDirs[iter]->cd();

    auto truthDir = pileupDirs[iter]->mkdir("Truth");
    truthDir->cd();

      pileupTimes[iter] = new TH1F("pileupTimes", "pileupTimes; time (ns); Events", nBins, 0 , histMaxTime);
      pileupTimes_threshold[iter] = new TH1F("pileupTimes_threshold", "pileupTimes_threshold; time (ns); Events", nBins, 0 , histMaxTime);

      pileupEnergies[iter] = new TH1F("pileupEnergies", "pileupEnergies; Energy (MeV); Events", 1000, 0, 10000);
          
      pileupTimes_U[iter] = new TH1F("pileupTimes_U", "pileupTimes_U; time (ns); Events", nBins, 0 , histMaxTime);
      pileupTimes_V[iter] = new TH1F("pileupTimes_V", "pileupTimes_V; time (ns); Events", nBins, 0 , histMaxTime);

    auto errorsDir = truthDir->mkdir("Errors");
    errorsDir->cd();

      errorsNeither[iter] = new TH1F("errorsNeither", "errorsNeither; time (ns); Events", nBins, 0 , histMaxTime);
      errorsBoth[iter] = new TH1F("errorsBoth", "errorsBoth; time (ns); Events", nBins, 0 , histMaxTime);

    auto shadowDir = pileupDirs[iter]->mkdir("Shadow");
    shadowDir->cd();

      pileupTimes_shadow[iter] = new TH1F("pileupTimes_shadow", "pileupTimes_shadow; time (ns); Events", nBins, 0 , histMaxTime);
      pileupTimes_shadow_threshold[iter] = new TH1F("pileupTimes_shadow_threshold", "pileupTimes_shadow_threshold; time (ns); Events", nBins, 0 , histMaxTime);

      pileupEnergies_shadow[iter] = new TH1F("pileupEnergies_shadow", "pileupEnergies_shadow; Energy (MeV); Events", 1000, 0, 10000);
      pileupEnergies_shadow_early[iter] = new TH1F("pileupEnergies_shadow_early", "pileupEnergies_shadow_early; Energy (MeV); Events", 1000, 0, 10000);
      pileupEnergies_shadow_late[iter] = new TH1F("pileupEnergies_shadow_late", "pileupEnergies_shadow_late; Energy (MeV); Events", 1000, 0, 10000);
          
      pileupTimes_U_shadow[iter] = new TH1F("pileupTimes_U_shadow", "pileupTimes_U_shadow; time (ns); Events", nBins, 0 , histMaxTime);
      pileupTimes_V_shadow[iter] = new TH1F("pileupTimes_V_shadow", "pileupTimes_V_shadow; time (ns); Events", nBins, 0 , histMaxTime);

        pileupEnergies_shadow[iter]->SetLineColor(2);
        pileupTimes_shadow[iter]->SetLineColor(2);
        pileupTimes_shadow_threshold[iter]->SetLineColor(2);
        pileupTimes_U_shadow[iter]->SetLineColor(2);
        pileupTimes_V_shadow[iter]->SetLineColor(2);

    auto TS_dir = shadowDir->mkdir("TimeSlice");
    TS_dir->cd();
      for (int timeSlice = 0; timeSlice < numTimeSlices; ++timeSlice) pileupEnergies_shadow_timeSlice[iter][timeSlice] = new TH1F(Form("pileupEnergies_shadow_tS%i",timeSlice), Form("pileupEnergies_shadow_tS%i; Energy (MeV); Events",timeSlice), 1000, 0, 10000);
    auto g2_dir = shadowDir->mkdir("g2");
    g2_dir->cd();
      for (int g2phase = 0; g2phase < 8; ++g2phase) pileupEnergies_shadow_g2phase[iter][g2phase] = new TH1F(Form("pileupEnergies_shadow_tS%i",g2phase), Form("pileupEnergies_shadow_tS%i; Energy (MeV); Events",g2phase), 1000, 0, 10000);


    auto errorsDir_shadow = shadowDir->mkdir("Errors");
    errorsDir_shadow->cd();

      errorsNeither_shadow[iter] = new TH1F("errorsNeither_shadow", "errorsNeither_shadow; time (ns); Events", nBins, 0 , histMaxTime);
      errorsBoth_shadow[iter] = new TH1F("errorsBoth_shadow", "errorsBoth_shadow; time (ns); Events", nBins, 0 , histMaxTime);

  }

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  std::vector<std::pair<double, double> > trueHits;
  std::vector<std::pair<double, double> > measuredHits;
    std::vector<std::vector<std::pair<double, double> > > measuredHitsTruth;

    std::vector<std::pair<double, double> > trueNegativePileup;
    std::vector<std::pair<double, double> > truePositivePileup;

    std::vector<std::pair<double, double> > shadowNegativePileup;
    std::vector<std::pair<double, double> > shadowPositivePileup;

  std::map<int,long long> hitCombinations;

/////////////////////////////////////////////////////////////////////////////////////

  double randEntries = randFill->Poisson(nPts);
  double numFills = randFill->Poisson(nPts/nPerFillPerCalo);

    int countdown = 0.1 * numFills;
    int tenPercent = 0;

  for (double fillNum = 0; fillNum < numFills; ++fillNum)
  {
    if(--countdown == 0) // progress output
    {
      tenPercent++;
      cout << tenPercent << "0%" << " completed" << endl;
      countdown = 0.1 * numFills;
    }


    double hitsInThisFill = randFill->Poisson(nPerFillPerCalo);

    double hit_time;
    double hit_energy_raw;

      // fill a vector of hits corresponding to a "fill"

      for (int j = 0; j < hitsInThisFill; ++j)
      {
        if(monoEnergetic){
          hit_time = histMaxTime * randFill->Uniform();
          hit_energy_raw = 1000.;
          // hit_energy_raw = 3100. * randFill->Uniform();
        }
        else{
          initialHist->GetRandom2(hit_time, hit_energy_raw);          
        }

        trueHits.emplace_back(hit_time, hit_energy_raw);
      }
      std::sort(trueHits.begin(), trueHits.end());

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


    // for (int iter = 0; iter < totalIters; ++iter) // could put this up here if I want to change the ADT from iter to iter - will slow things down though and I need to be careful with random seeds and all that
    // {

      // loop through the fill hits and create the pileup spectrum

      uint k = 0;
      while (k < trueHits.size())
      {
        std::vector<std::pair<double, double> > trueHitsInMeasHit;

        if(k != trueHits.size()-1 && abs(trueHits.at(k).first - trueHits.at(k+1).first) < artificialDeadtime) // double pileup
        {
          std::vector<std::pair<double, double> > pileupSinglets;
          pileupSinglets.push_back(trueHits.at(k));
          pileupSinglets.push_back(trueHits.at(k+1));

            uint kk = k+2;
            while(kk != trueHits.size() && (abs(trueHits.at(kk-1).first - trueHits.at(kk).first) < artificialDeadtime/2.)) // triple and higher pileup - not clear what the factor at the end should be 
            // while(kk != trueHits.size() && (abs(trueHits.at(kk-1).first - trueHits.at(kk).first) < artificialDeadtime)) // triple and higher pileup
            {
              pileupSinglets.push_back(trueHits.at(kk));
              kk++;
            }

          double addedEnergy = pileupAddedEnergy(pileupSinglets);
          double weightedTime = energyWeightedTime(pileupSinglets, 0); // 0 for shadow gap time

          for (uint j = 0; j < pileupSinglets.size(); ++j){
            trueNegativePileup.emplace_back(weightedTime, pileupSinglets.at(j).second);
            trueHitsInMeasHit.push_back(pileupSinglets.at(j));
          } 
          
          truePositivePileup.emplace_back(weightedTime, addedEnergy);
          measuredHits.emplace_back(weightedTime, addedEnergy);

          k = kk;
        }
        else{
          measuredHits.push_back(trueHits.at(k));
          trueHitsInMeasHit.push_back(trueHits.at(k));
          k++;
        }

        measuredHitsTruth.push_back(trueHitsInMeasHit);  
      } // end while

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  for (int iter = 0; iter < totalIters; ++iter)
  {

    // construct pileup from shadow pulses in 'data'

    std::vector<int> shadowIndices;
    std::vector<int> tripleShadowIndices;

      uint j = 0;
      while (j < measuredHits.size()-1)
      {
        uint indSave = j;

        checkForShadowPulses(j, &measuredHits, shadowIndices, tripleShadowIndices, shadowDeadTimes[iter], shadowGapTimes[iter], shadowSecondDeadTimes[iter]);

        int nTriggerHits = measuredHitsTruth.at(indSave).size();
        int nShadowHits = 0;
        int nShadowHits2 = 0;


        if(shadowIndices.size() > 0) // size is 1 when ADT and SDT is the same
        {
          std::vector<std::pair<double, double> > pileupSinglets;
          pileupSinglets.push_back(measuredHits.at(indSave));

          for (uint i = 0; i < shadowIndices.size(); ++i)
          {
            int shadowIndex = shadowIndices.at(i);
            pileupSinglets.push_back(measuredHits.at(shadowIndex));

            nShadowHits += measuredHitsTruth.at(shadowIndex).size();
          }

          double addedEnergy = pileupAddedEnergy(pileupSinglets);
          double weightedTime = energyWeightedTime(pileupSinglets, shadowGapTimes[iter]);

          for (uint i = 0; i < pileupSinglets.size(); ++i) shadowNegativePileup.emplace_back(weightedTime, pileupSinglets.at(i).second);
          shadowPositivePileup.emplace_back(weightedTime, addedEnergy);

          // fill error histograms - currently only accounts for errors for doublets (0 and 1 entries in pileupSinglets)
          if(addedEnergy > energyThresholds[iter])
          {
                 if(pileupSinglets.at(0).second > energyThresholds[iter] && pileupSinglets.at(1).second > energyThresholds[iter]) errorsBoth_shadow[iter]->Fill(weightedTime);
            else if(pileupSinglets.at(0).second < energyThresholds[iter] && pileupSinglets.at(1).second < energyThresholds[iter]) errorsNeither_shadow[iter]->Fill(weightedTime);
          }

/////////////////////////////////////////////////////////////////////////////////////

          // Triplets - can comment this block out to only included doublets
          // Combinatorics were determined emperically

          if(performTripletCorrection && tripleShadowIndices.size() > 0) {

            double triggerEnergy = pileupSinglets.at(0).second;
            double addedEnergyDoublet = pileupSinglets.at(1).second; // energy of pulse(s) in first shadow window - only the first for now since there's only 1 pulse if SDT = ADT
            double addedEnergyTriplet = 0; // energy of pulse(s) in second shadow window

            for (uint i = 0; i < tripleShadowIndices.size(); ++i) {
               int tripletIndex = tripleShadowIndices.at(i);
               addedEnergyTriplet += measuredHits.at(tripletIndex).second;
               // nShadowHits2 += measuredHitsTruth.at(tripletIndex).size();
            }

            int randInt = tripRand->Integer(2);

/////////////////////////////////////////////////////////////////////////////////////
            cout << "Triplet stuff commented out for the time being." << endl;

        // seems to work but the numbers were just hacked in (when ADT = SDT2)

            // for(int n = 0; n < 1; n++) shadowNegativePileup.emplace_back(weightedTime, triggerEnergy + addedEnergyDoublet + addedEnergyTriplet);

            // for(int n = 0; n < 1; n++) shadowPositivePileup.emplace_back(weightedTime, triggerEnergy + addedEnergyDoublet);
            // for(int n = 0; n < 1; n++) shadowPositivePileup.emplace_back(weightedTime, triggerEnergy + addedEnergyTriplet);
            // for(int n = 0; n < 1; n++) shadowPositivePileup.emplace_back(weightedTime, addedEnergyDoublet + addedEnergyTriplet);
            // if(randInt > 0) shadowPositivePileup.emplace_back(weightedTime, addedEnergyDoublet + addedEnergyTriplet);

            // for(int n = 0; n < 1; n++) shadowNegativePileup.emplace_back(weightedTime, triggerEnergy);
            // for(int n = 0; n < 1; n++) shadowNegativePileup.emplace_back(weightedTime, addedEnergyDoublet);
            // for(int n = 0; n < 2; n++) shadowNegativePileup.emplace_back(weightedTime, addedEnergyTriplet);

        // what's seemed to work best for SDT2 = ADT/2 (for 5 ns ADT)

            // for(int n = 0; n < 3; n++) shadowNegativePileup.emplace_back(weightedTime, triggerEnergy + addedEnergyDoublet + addedEnergyTriplet);

            // for(int n = 0; n < 3; n++) shadowPositivePileup.emplace_back(weightedTime, triggerEnergy + addedEnergyDoublet);
            // for(int n = 0; n < 3; n++) shadowPositivePileup.emplace_back(weightedTime, triggerEnergy + addedEnergyTriplet);
            // for(int n = 0; n < 2; n++) shadowPositivePileup.emplace_back(weightedTime, addedEnergyDoublet + addedEnergyTriplet);

            // for(int n = 0; n < 2; n++) shadowNegativePileup.emplace_back(weightedTime, triggerEnergy);
            // for(int n = 0; n < 3; n++) shadowNegativePileup.emplace_back(weightedTime, addedEnergyDoublet);
            // for(int n = 0; n < 2; n++) shadowNegativePileup.emplace_back(weightedTime, addedEnergyTriplet);

        // play - doesn't work for 10 ns ADT but should it? the fixes that I would need to make seem possibly contradictory to 5 or 6 ns changes
/*
            for(int n = 0; n < 3; n++) shadowNegativePileup.emplace_back(weightedTime, triggerEnergy + addedEnergyDoublet + addedEnergyTriplet);

            for(int n = 0; n < 3; n++) shadowPositivePileup.emplace_back(weightedTime, triggerEnergy + addedEnergyDoublet);
            for(int n = 0; n < 3; n++) shadowPositivePileup.emplace_back(weightedTime, triggerEnergy + addedEnergyTriplet);
            for(int n = 0; n < 2; n++) shadowPositivePileup.emplace_back(weightedTime, addedEnergyDoublet + addedEnergyTriplet);
            // if(randInt > 0) shadowPositivePileup.emplace_back(weightedTime, addedEnergyDoublet + addedEnergyTriplet);

            for(int n = 0; n < 2; n++) shadowNegativePileup.emplace_back(weightedTime, triggerEnergy);
            for(int n = 0; n < 3; n++) shadowNegativePileup.emplace_back(weightedTime, addedEnergyDoublet);
            for(int n = 0; n < 2; n++) shadowNegativePileup.emplace_back(weightedTime, addedEnergyTriplet);
            // if(randInt > 0) shadowNegativePileup.emplace_back(weightedTime, addedEnergyTriplet);
*/
/////////////////////////////////////////////////////////////////////////////////////
/*
        // this works when only comparing deadtimes to first hit (hard window of 5 ns) when making pileup

            if(randInt > 0) for(int n = 0; n < 1; n++) shadowNegativePileup.emplace_back(weightedTime, triggerEnergy + addedEnergyDoublet + addedEnergyTriplet);
            for(int n = 0; n < 1; n++) shadowNegativePileup.emplace_back(weightedTime, triggerEnergy + addedEnergyDoublet + addedEnergyTriplet);

            if(randInt > 0) for(int n = 0; n < 1; n++) shadowPositivePileup.emplace_back(weightedTime, triggerEnergy + addedEnergyDoublet);
            for(int n = 0; n < 1; n++) shadowPositivePileup.emplace_back(weightedTime, triggerEnergy + addedEnergyDoublet);

            for(int n = 0; n < 2; n++) shadowPositivePileup.emplace_back(weightedTime, addedEnergyDoublet + addedEnergyTriplet);

            if(randInt > 0) for(int n = 0; n < 1; n++) shadowNegativePileup.emplace_back(weightedTime, triggerEnergy);   
            for(int n = 0; n < 1; n++) shadowNegativePileup.emplace_back(weightedTime, addedEnergyDoublet);
            for(int n = 0; n < 1; n++) shadowNegativePileup.emplace_back(weightedTime, addedEnergyTriplet);
*/

/////////////////////////////////////////////////////////////////////////////////////

        // original set
/*
            for(int n = 0; n < 3; n++) shadowNegativePileup.emplace_back(weightedTime, triggerEnergy + addedEnergyDoublet + addedEnergyTriplet); // 3
 
            for(int n = 0; n < 3; n++) shadowPositivePileup.emplace_back(weightedTime, triggerEnergy + addedEnergyDoublet);      // 3
            for(int n = 0; n < 4; n++) shadowPositivePileup.emplace_back(weightedTime, addedEnergyDoublet + addedEnergyTriplet); // 4
            // for(int n = 0; n < 1; n++) shadowPositivePileup.emplace_back(weightedTime, triggerEnergy + addedEnergyTriplet);   // 0

            for(int n = 0; n < 1; n++) shadowNegativePileup.emplace_back(weightedTime, triggerEnergy);      // 1
            for(int n = 0; n < 2; n++) shadowNegativePileup.emplace_back(weightedTime, addedEnergyDoublet); // 2
            for(int n = 0; n < 2; n++) shadowNegativePileup.emplace_back(weightedTime, addedEnergyTriplet); // 2
*/
          } // end triplet additions

/////////////////////////////////////////////////////////////////////////////////////

          pileupSinglets.clear();
        } // end shadow indices size if statement

/////////////////////////////////////////////////////////////////////////////////////
      // try running this even when shadow indices size are 0
            for (uint i = 0; i < tripleShadowIndices.size(); ++i) {
               int tripletIndex = tripleShadowIndices.at(i);
               nShadowHits2 += measuredHitsTruth.at(tripletIndex).size();
            }
/////////////////////////////////////////////////////////////////////////////////////

        shadowIndices.clear();
        tripleShadowIndices.clear();

        int hitCombo = nTriggerHits*100 + nShadowHits*10 + nShadowHits2;
        hitCombinations[hitCombo]++;

      }  //end while

/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

    // fill pileup and pileup shadow histograms

    // can make a separate loop over iterations for filling these pileup threshold histograms - which might be useful further down the line

    double halfPeriodGuess_temp = gm2PeriodGuesses[iter]/2.;
    double totalChance_temp = exp(halfPeriodGuess_temp/defaultLifetime) + exp(-halfPeriodGuess_temp/defaultLifetime) + 2.;
    double percentChanceUPlus_temp = exp(halfPeriodGuess_temp/defaultLifetime) / totalChance_temp;
    double percentChanceUMinus_temp = exp(-halfPeriodGuess_temp/defaultLifetime) / totalChance_temp;

    double randNum_temp;

    double pTime, pEnergy;

    warning here // fix pileup UV randomization as in data - and for triplets

/////////////////////////////////////////////////////////////////////////////////////

    for (uint i = 0; i < trueNegativePileup.size(); ++i)
    {
      pTime = trueNegativePileup.at(i).first;
      pEnergy = trueNegativePileup.at(i).second;

      pileupTimes[iter]->Fill(pTime, -1);
      pileupEnergies[iter]->Fill(pEnergy, -1);

        if(pEnergy > energyThresholds[iter])
        {
          pileupTimes_threshold[iter]->Fill(pTime, -1);

          randNum_temp = randUVPileup[iter]->Uniform();
            if     (randNum_temp < percentChanceUPlus_temp) pileupTimes_U[iter]->Fill(pTime - halfPeriodGuess_temp, -1);
            else if(randNum_temp < (percentChanceUPlus_temp + percentChanceUMinus_temp)) pileupTimes_U[iter]->Fill(pTime + halfPeriodGuess_temp, -1);
            else if(randNum_temp < 1) pileupTimes_V[iter]->Fill(pTime, -1);
        }
    }
    for (uint i = 0; i < truePositivePileup.size(); ++i)
    {
      pTime = truePositivePileup.at(i).first;
      pEnergy = truePositivePileup.at(i).second;

      pileupTimes[iter]->Fill(pTime);
      pileupEnergies[iter]->Fill(pEnergy);

        if(pEnergy > energyThresholds[iter])
        {
          pileupTimes_threshold[iter]->Fill(pTime);

          randNum_temp = randUVPileup[iter]->Uniform();
            if     (randNum_temp < percentChanceUPlus_temp) pileupTimes_U[iter]->Fill(pTime - halfPeriodGuess_temp);
            else if(randNum_temp < (percentChanceUPlus_temp + percentChanceUMinus_temp)) pileupTimes_U[iter]->Fill(pTime + halfPeriodGuess_temp);
            else if(randNum_temp < 1) pileupTimes_V[iter]->Fill(pTime);
        }
    }

/////////////////////////////////////////////////////////////////////////////////////

    for (uint i = 0; i < shadowNegativePileup.size(); ++i)
    {
      pTime = shadowNegativePileup.at(i).first + pileupTimeShifts[iter];
      pEnergy = shadowNegativePileup.at(i).second;

      pileupTimes_shadow[iter]->Fill(pTime, -1);
      pileupEnergies_shadow[iter]->Fill(pEnergy, -1);
      pileupEnergies_shadow_g2phase[iter][int(fmod(pTime,g2Period)/(g2Period/8))]->Fill(pEnergy, -1);

      if(pTime < 250000) pileupEnergies_shadow_early[iter]->Fill(pEnergy, -1);
      else if(pTime > 250000) pileupEnergies_shadow_late[iter]->Fill(pEnergy, -1);

      uint timeSlice = 0;
      while(timeSlice < timeSliceLowEdges.size() && pTime > timeSliceLowEdges.at(timeSlice)) timeSlice++;
      pileupEnergies_shadow_timeSlice[iter][timeSlice-1]->Fill(pEnergy, -1);

        if(pEnergy > energyThresholds[iter])
        {
          pileupTimes_shadow_threshold[iter]->Fill(pTime, -1);

          randNum_temp = randUVPileup[iter]->Uniform();
            if     (randNum_temp < percentChanceUPlus_temp) pileupTimes_U_shadow[iter]->Fill(pTime - halfPeriodGuess_temp, -1);
            else if(randNum_temp < (percentChanceUPlus_temp + percentChanceUMinus_temp)) pileupTimes_U_shadow[iter]->Fill(pTime + halfPeriodGuess_temp, -1);
            else if(randNum_temp < 1) pileupTimes_V_shadow[iter]->Fill(pTime, -1);
        }
    }
    for (uint i = 0; i < shadowPositivePileup.size(); ++i)
    {
      pTime = shadowPositivePileup.at(i).first + pileupTimeShifts[iter];
      pEnergy = shadowPositivePileup.at(i).second;

      pileupTimes_shadow[iter]->Fill(pTime);
      pileupEnergies_shadow[iter]->Fill(pEnergy);
      pileupEnergies_shadow_g2phase[iter][int(fmod(pTime,g2Period)/(g2Period/8))]->Fill(pEnergy);

      if(pTime < 250000) pileupEnergies_shadow_early[iter]->Fill(pEnergy);
      else if(pTime > 250000) pileupEnergies_shadow_late[iter]->Fill(pEnergy);

      uint timeSlice = 0;
      while(timeSlice < timeSliceLowEdges.size() && pTime > timeSliceLowEdges.at(timeSlice)) timeSlice++;
      pileupEnergies_shadow_timeSlice[iter][timeSlice-1]->Fill(pEnergy);

        if(pEnergy > energyThresholds[iter])
        {
          pileupTimes_shadow_threshold[iter]->Fill(pTime);

          randNum_temp = randUVPileup[iter]->Uniform();
            if     (randNum_temp < percentChanceUPlus_temp) pileupTimes_U_shadow[iter]->Fill(pTime - halfPeriodGuess_temp);
            else if(randNum_temp < (percentChanceUPlus_temp + percentChanceUMinus_temp)) pileupTimes_U_shadow[iter]->Fill(pTime + halfPeriodGuess_temp);
            else if(randNum_temp < 1) pileupTimes_V_shadow[iter]->Fill(pTime);
        }
    }


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
if(fillHitHists){
        for (int subIter = 0; subIter < 2; ++subIter) // only run for 2 iterations if I want to create one set of true hits histograms and one set of measured hits histograms (iters 2 and 3 are then filled down below)
        {
            std::vector<std::pair<double, double> > myHits;

            if(subIter == 0) myHits = trueHits;
            else myHits = measuredHits;

            for (uint hitJ = 0; hitJ < myHits.size(); ++hitJ) // loop through the new hits and fill histograms
            {
              double time = myHits.at(hitJ).first;
              double energy_raw = myHits.at(hitJ).second;

              double randNum = randIterUV[iter]->Uniform(); // this fills the UV histograms, and won't be the same for sub iterations 0 and 1 because true and measured hits have different sizes, but should be the same from iteration to iteration

/////////////////////////////////////////////////////////////////////////////////////
            
              double gainSagMultiple = gainSagFunctions[iter]->Eval(time);
              double gainCorrMultiple = gainCorrFunctions[iter]->Eval(time);


              double energy_measured;

              // double gausSmearSag = randGenerator->Gaus(gainSagMultiple, 0.01*gainSagMultiple);
              // if(iter == 0) energy_measured = energy_raw; // don't modify the 0th iteration energy
              // else energy_measured = energy_raw * gausSmearSag * gainCorrMultiple;

              // energy_measured = energy_raw * gainSagMultiple * gainCorrMultiple;
              energy_measured = energy_raw;

/////////////////////////////////////////////////////////////////////////////////////

              energies[iter][subIter]->Fill(energy_measured);

              if(time < 250000) energies_early[iter][subIter]->Fill(energy_measured);
              if(time > 250000) energies_late[iter][subIter]->Fill(energy_measured);

              energies_g2phase[iter][subIter][int(fmod(time,g2Period)/(g2Period/8))]->Fill(energy_measured);

              uint timeSlice = 0;
              while(timeSlice < timeSliceLowEdges.size() && time > timeSliceLowEdges.at(timeSlice)) timeSlice++;
              energies_timeSlice[iter][subIter][timeSlice-1]->Fill(energy_measured);

/////////////////////////////////////////////////////////////////////////////////////

              if(hitJ != myHits.size()-1) deltaT_hits[iter][subIter]->Fill(myHits.at(hitJ+1).first - time);

              if(energy_measured > energyThresholds[iter])
              {
                    toyFiveParamHists[iter][subIter]->Fill(time);

                    if(hitJ != myHits.size()-1) deltaT_hits_threshold[iter][subIter]->Fill(myHits.at(hitJ+1).first - time);

                    double halfPeriodGuess = gm2PeriodGuesses[iter]/2;

                    double totalChance = exp(halfPeriodGuess/defaultLifetime) + exp(-halfPeriodGuess/defaultLifetime) + 2.;
                    double percentChanceUPlus = exp(halfPeriodGuess/defaultLifetime) / totalChance;
                    double percentChanceUMinus = exp(-halfPeriodGuess/defaultLifetime) / totalChance;

                    if     (randNum < percentChanceUPlus) toyUHists[iter][subIter]->Fill(time - halfPeriodGuess); // careful with the signs here - the U+ hist is filled with pulses shifted by t - T/2, and weighted by e^+T/2tau, and vice versa for U-
                    else if(randNum < (percentChanceUPlus + percentChanceUMinus)) toyUHists[iter][subIter]->Fill(time + halfPeriodGuess);
                    else if(randNum < 1) toyVHists[iter][subIter]->Fill(time);
              }

            } // end loop over hits

        } // end loop over sub iterations
} // end if fillHitHists

        shadowNegativePileup.clear();
        shadowPositivePileup.clear();

      } // end loop over iterations


      trueHits.clear(); // these are the same from iter to iter
      measuredHits.clear();

      measuredHitsTruth.clear();

      trueNegativePileup.clear();
      truePositivePileup.clear();

    } // end outer for loop over entries/nFillPerCalo


/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

  // this code for filling sub iterations 2 and 3 with the same hits as sub iteration 1, but to which I will apply different pileup corrections to
  // subIter 2 has the true pileup correction, subIter 3 has the shadow pileup correction

  for (int iter = 0; iter < totalIters; ++iter)
  {

    for (int bin = 0; bin <= energies[iter][1]->GetNbinsX()+1; bin++)
    {
      energies[iter][2]->SetBinContent(bin, energies[iter][1]->GetBinContent(bin));
      energies[iter][3]->SetBinContent(bin, energies[iter][1]->GetBinContent(bin));

      // energies_early[iter][2]->SetBinContent(bin, energies_early[iter][1]->GetBinContent(bin));
      // energies_early[iter][3]->SetBinContent(bin, energies_early[iter][1]->GetBinContent(bin));
      
      // energies_late[iter][2]->SetBinContent(bin, energies_late[iter][1]->GetBinContent(bin));
      // energies_late[iter][3]->SetBinContent(bin, energies_late[iter][1]->GetBinContent(bin));   

      // for (int timeSlice = 0; timeSlice < numTimeSlices; ++timeSlice)
      // {
      //   energies_timeSlice[iter][2][timeSlice]->SetBinContent(bin, energies_timeSlice[iter][1][timeSlice]->GetBinContent(bin));
      //   energies_timeSlice[iter][3][timeSlice]->SetBinContent(bin, energies_timeSlice[iter][1][timeSlice]->GetBinContent(bin));
      // }
    } 

    energies[iter][2]->ResetStats();
    energies[iter][3]->ResetStats();

    // energies_early[iter][2]->ResetStats();
    // energies_early[iter][3]->ResetStats();

    // energies_late[iter][2]->ResetStats();
    // energies_late[iter][3]->ResetStats();

      // for (int timeSlice = 0; timeSlice < numTimeSlices; ++timeSlice)
      // {
      //   energies_timeSlice[iter][2][timeSlice]->ResetStats();
      //   energies_timeSlice[iter][3][timeSlice]->ResetStats();
      // }

/////////////////////////////////////////////////////////////////////////////////////

    for (int bin = 0; bin <= toyFiveParamHists[iter][1]->GetNbinsX()+1; bin++)
    {
      toyFiveParamHists[iter][2]->SetBinContent(bin, toyFiveParamHists[iter][1]->GetBinContent(bin));
      toyUHists[iter][2]->SetBinContent(bin, toyUHists[iter][1]->GetBinContent(bin));
      toyVHists[iter][2]->SetBinContent(bin, toyVHists[iter][1]->GetBinContent(bin));

      toyFiveParamHists[iter][3]->SetBinContent(bin, toyFiveParamHists[iter][1]->GetBinContent(bin));
      toyUHists[iter][3]->SetBinContent(bin, toyUHists[iter][1]->GetBinContent(bin));
      toyVHists[iter][3]->SetBinContent(bin, toyVHists[iter][1]->GetBinContent(bin));
    } 

    toyFiveParamHists[iter][2]->ResetStats();
    toyUHists[iter][2]->ResetStats();
    toyVHists[iter][2]->ResetStats();
    toyFiveParamHists[iter][3]->ResetStats();
    toyUHists[iter][3]->ResetStats();
    toyVHists[iter][3]->ResetStats();


    // apply pileup correction

    pileupTimes_threshold[iter]->Scale(pileupScaleFactors[iter]);
    pileupEnergies[iter]->Scale(pileupScaleFactors[iter]);
    pileupTimes_U[iter]->Scale(pileupScaleFactors[iter]);
    pileupTimes_V[iter]->Scale(pileupScaleFactors[iter]);

        toyFiveParamHists[iter][2]->Add(pileupTimes_threshold[iter], -1); // perfect correction
        energies[iter][2]->Add(pileupEnergies[iter], -1);
        toyUHists[iter][2]->Add(pileupTimes_U[iter], -1); // won't be the perfect correction for UV - because measured hits has a different size than the true hits leading to different random numbers when filling UV hists even though each iteration random number generator starts with the same seed
        toyVHists[iter][2]->Add(pileupTimes_V[iter], -1);

    // scale here instead of in runRatioToyMod since it's easier - might want to move it later
        pileupTimes_shadow_threshold[iter]->Scale(pileupScaleFactors[iter]);
        pileupEnergies_shadow[iter]->Scale(pileupScaleFactors[iter]);
        pileupTimes_U_shadow[iter]->Scale(pileupScaleFactors[iter]);
        pileupTimes_V_shadow[iter]->Scale(pileupScaleFactors[iter]);

        // errorsNeither_shadow // do I want to scale these errors?
        // errorsBoth_shadow

            toyFiveParamHists[iter][3]->Add(pileupTimes_shadow_threshold[iter], -1); // shadow correction            
            energies[iter][3]->Add(pileupEnergies_shadow[iter], -1);
            toyUHists[iter][3]->Add(pileupTimes_U_shadow[iter], -1);
            toyVHists[iter][3]->Add(pileupTimes_V_shadow[iter], -1);

        // energies_early[iter][3]->Add(pileupEnergies_shadow[iter], -1);


  // more time slice stuff

      // for (int timeSlice = 0; timeSlice < numTimeSlices; ++timeSlice)
      // {
      //   pileupEnergies_shadow_timeSlice[iter][timeSlice]->Scale(pileupScaleFactors[iter]); // should just be 1 for now
      //   energies_timeSlice[iter][3][timeSlice]->Add(pileupEnergies_shadow_timeSlice[iter][timeSlice], -1);
      // }

  
  } // end iter loop


  for (int iter = 0; iter < totalIters; ++iter)
  {
    for (int subIter = 0; subIter < subIters; ++subIter)
    {
      toyNumHists[iter][subIter]->Add(toyUHists[iter][subIter], toyVHists[iter][subIter], -1, 1);
      toyDenomHists[iter][subIter]->Add(toyUHists[iter][subIter], toyVHists[iter][subIter]);
    }  
  } // end iter loop

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  // replace errors for pileup subtracted histograms - really should move this into fitting stage
  // start with iter 3 5 param hist - ratio errors are modified in fitting stage

  for (int iter = 0; iter < totalIters; ++iter)
  {

  auto gammaDir = pileupDirs[iter]->mkdir("errorCompareDir");
  gammaDir->cd();

  TH1F* errorMultiplierByBin = pileupErrorMultiplier(errorsBoth_shadow[iter], errorsNeither_shadow[iter], toyFiveParamHists[iter][3]);
  TH1F* gammaHist = pileupErrorMultiplierApproximation(errorsBoth_shadow[iter], errorsNeither_shadow[iter], toyFiveParamHists[iter][3], 0);

    for (int bin = 1; bin <= toyFiveParamHists[iter][3]->GetNbinsX(); ++bin)
    {
      double newError = 0;

           if(iter == 0) newError = sqrt(toyFiveParamHists[iter][3]->GetBinContent(bin));
      else if(iter == 1) newError = sqrt(toyFiveParamHists[iter][3]->GetBinContent(bin) + 2.*errorsNeither_shadow[iter]->GetBinContent(bin) + 6.*errorsBoth_shadow[iter]->GetBinContent(bin));
      else if(iter == 2) newError = sqrt(toyFiveParamHists[iter][3]->GetBinContent(bin)) * gammaHist->GetBinContent(bin);

      toyFiveParamHists[iter][3]->SetBinError(bin, newError);
    }
  }


/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
for (int iter = 0; iter < totalIters; ++iter)
{
  auto diffDir = pileupDirs[iter]->mkdir("Diff");
  diffDir->cd();

/////////////////////////////////////////////////////////////////////////////////////

  THStack* timesHistsStack = new THStack("timesHistsStack","MC Pileup Times");
  timesHistsStack->Add(pileupTimes[iter]);
  timesHistsStack->Add(pileupTimes_shadow[iter]);

  auto canvas_Ptimes = new TCanvas("canvas_Ptimes","canvas_Ptimes",200,10,1200,1000);
  TPad* padT1 = new TPad("padT1", "padT1", 0.02, 0.5, 0.98, 0.98);
  padT1->Draw();
  padT1->cd();
    timesHistsStack->Draw("nostack,hist");

    timesHistsStack->GetXaxis()->SetTitle("time (ns)");
    timesHistsStack->GetYaxis()->SetTitle("Events");
    padT1->Modified();

       auto tHS_legend = new TLegend(0.7,0.7,0.9,0.9);
       tHS_legend->AddEntry(pileupTimes[iter],"True pileup","l");
       tHS_legend->AddEntry(pileupTimes_shadow[iter],"Shadow pileup correction","l");
       tHS_legend->SetBorderSize(1);
       tHS_legend->Draw();

  canvas_Ptimes->cd();
   TPad* padT2 = new TPad("padT2", "padT2", 0.02, 0.25, 0.98, 0.5);
  padT2->Draw();
  padT2->cd();
    TH1F* PtimesDiff = (TH1F*) pileupTimes_shadow[iter]->Clone("PtimesDiff");
    PtimesDiff->Add(pileupTimes[iter], -1);
    PtimesDiff->SetTitle("Residual (Shadow - Truth)");
    PtimesDiff->Draw("hist");

  canvas_Ptimes->cd();
   TPad* padT3 = new TPad("padT3", "padT3", 0.02, 0.02, 0.98, 0.25);
  padT3->Draw();
  padT3->cd();
    TH1F* PtimesSigDiff = (TH1F*) PtimesDiff->Clone("PtimesSigDiff");
      for (int bin = 1; bin <= PtimesSigDiff->GetNbinsX(); ++bin) if(PtimesSigDiff->GetBinError(bin) != 0) PtimesSigDiff->SetBinContent(bin, PtimesSigDiff->GetBinContent(bin)/PtimesSigDiff->GetBinError(bin));
    PtimesSigDiff->SetTitle("Residual (Shadow - Truth) / Bin Error");
    PtimesSigDiff->Draw("hist");

  pileupDirs[iter]->cd();
  canvas_Ptimes->Write();
  delete canvas_Ptimes;

/////////////////////////////////////////////////////////////////////////////////////

  diffDir->cd();

  THStack* timesHistsStack_threshold = new THStack("timesHistsStack_threshold","MC Pileup Times Above Threshold");
  timesHistsStack_threshold->Add(pileupTimes_threshold[iter]);
  timesHistsStack_threshold->Add(pileupTimes_shadow_threshold[iter]);

  auto canvas_Ptimes_threshold = new TCanvas("canvas_Ptimes_threshold","canvas_Ptimes_threshold",200,10,1200,1000);
  TPad* padT1_threshold = new TPad("padT1_threshold", "padT1_threshold", 0.02, 0.5, 0.98, 0.98);
  padT1_threshold->Draw();
  padT1_threshold->cd();
    timesHistsStack_threshold->Draw("nostack,hist");

    timesHistsStack_threshold->GetXaxis()->SetTitle("time (ns)");
    timesHistsStack_threshold->GetYaxis()->SetTitle("Events");
    padT1_threshold->Modified();

       auto tHS_threshold_legend = new TLegend(0.7,0.7,0.9,0.9);
       tHS_threshold_legend->AddEntry(pileupTimes_threshold[iter],"True pileup","l");
       tHS_threshold_legend->AddEntry(pileupTimes_shadow_threshold[iter],"Shadow pileup correction","l");
       tHS_threshold_legend->SetBorderSize(1);
       tHS_threshold_legend->Draw();

  canvas_Ptimes_threshold->cd();
   TPad* padT2_threshold = new TPad("padT2_threshold", "padT2_threshold", 0.02, 0.25, 0.98, 0.5);
  padT2_threshold->Draw();
  padT2_threshold->cd();
    TH1F* PtimesDiff_threshold = (TH1F*) pileupTimes_shadow_threshold[iter]->Clone("PtimesDiff_threshold");
    PtimesDiff_threshold->Add(pileupTimes_threshold[iter], -1);
    PtimesDiff_threshold->SetTitle("Residual (Shadow - Truth)");
    PtimesDiff_threshold->Draw("hist");

  canvas_Ptimes_threshold->cd();
   TPad* padT3_threshold = new TPad("padT3_threshold", "padT3_threshold", 0.02, 0.02, 0.98, 0.25);
  padT3_threshold->Draw();
  padT3_threshold->cd();
    TH1F* PtimesSigDiff_threshold = (TH1F*) PtimesDiff_threshold->Clone("PtimesSigDiff_threshold");
      for (int bin = 1; bin <= PtimesSigDiff_threshold->GetNbinsX(); ++bin) if(PtimesSigDiff_threshold->GetBinError(bin) != 0) PtimesSigDiff_threshold->SetBinContent(bin, PtimesSigDiff_threshold->GetBinContent(bin)/PtimesSigDiff_threshold->GetBinError(bin));
    PtimesSigDiff_threshold->SetTitle("Residual (Shadow - Truth) / Bin Error");
    PtimesSigDiff_threshold->Draw("hist");

  pileupDirs[iter]->cd();
  canvas_Ptimes_threshold->Write();
  delete canvas_Ptimes_threshold;

/////////////////////////////////////////////////////////////////////////////////////
  
/////////////////////////////////////////////////////////////////////////////////////

  diffDir->cd();

  THStack* energyHistsStack = new THStack("energyHistsStack","MC Pileup Energies");
  energyHistsStack->Add(pileupEnergies[iter]);
  energyHistsStack->Add(pileupEnergies_shadow[iter]);

  auto canvas_Penergies = new TCanvas("canvas_Penergies","canvas_Penergies",200,10,1200,1000);
  TPad* padE1 = new TPad("padE1", "padE1", 0.02, 0.5, 0.98, 0.98);
  padE1->Draw();
  padE1->cd();
    energyHistsStack->Draw("nostack,hist");
    
    energyHistsStack->GetXaxis()->SetTitle("Energy (MeV)");
    energyHistsStack->GetYaxis()->SetTitle("Events");
    padE1->Modified();

       auto eHS_legend = new TLegend(0.7,0.7,0.9,0.9);
       eHS_legend->AddEntry(pileupEnergies[iter],"True pileup","l");
       eHS_legend->AddEntry(pileupEnergies_shadow[iter],"Shadow pileup correction","l");
       eHS_legend->SetBorderSize(1);
       eHS_legend->Draw();

  canvas_Penergies->cd();
   TPad* padE2 = new TPad("padE2", "padE2", 0.02, 0.25, 0.98, 0.5);
  padE2->Draw();
  padE2->cd();
    TH1F* PenergiesDiff = (TH1F*) pileupEnergies_shadow[iter]->Clone("PenergiesDiff");
    PenergiesDiff->Add(pileupEnergies[iter], -1);
    PenergiesDiff->SetTitle("Residual (Shadow - Truth)");
    PenergiesDiff->Draw("hist");

  canvas_Penergies->cd();
   TPad* padE3 = new TPad("padE3", "padE3", 0.02, 0.02, 0.98, 0.25);
  padE3->Draw();
  padE3->cd();
    TH1F* PenergiesSigDiff = (TH1F*) PenergiesDiff->Clone("PenergiesSigDiff");
      for (int bin = 1; bin <= PenergiesSigDiff->GetNbinsX(); ++bin) if(PenergiesSigDiff->GetBinError(bin) != 0) PenergiesSigDiff->SetBinContent(bin, PenergiesSigDiff->GetBinContent(bin)/PenergiesSigDiff->GetBinError(bin));
    PenergiesSigDiff->SetTitle("Residual (Shadow - Truth) / Bin Error");
    PenergiesSigDiff->Draw("hist");

  pileupDirs[iter]->cd();
  canvas_Penergies->Write();
  delete canvas_Penergies;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

  // this should maybe be moved up if I want to apply this scale factor to the pileup instead of just calculate it - in that case would only want to apply it to iteration 3 or maybe even make another iteration 

  auto scaleDir = pileupDirs[iter]->mkdir("Scale");
  auto scaleTruthDir = scaleDir->mkdir("TruePileup");
  auto scaleShadowDir = scaleDir->mkdir("ShadowPileup");

  scaleTruthDir->cd();
  double scaleFactTruth = calcPileupScaleFactor(energies[iter][1], pileupEnergies[iter]);
  cout << "scale factor from truth should be 1 and is : " << scaleFactTruth << endl;

  scaleShadowDir->cd();
  double scaleFactShadow = calcPileupScaleFactor(energies[iter][1], pileupEnergies_shadow[iter]);
  cout << "scale factor from measured is (a scale factor might already have been applied before this keep in mind) : " << scaleFactShadow << endl;

    auto timeSliceDir = scaleShadowDir->mkdir("TimeSlice");
    timeSliceDir->cd();
      for (int timeSlice = 0; timeSlice < numTimeSlices; ++timeSlice)
      {
        auto tS_dir = timeSliceDir->mkdir(Form("TS%i",timeSlice));
        tS_dir->cd();
        calcPileupScaleFactor(energies_timeSlice[iter][1][timeSlice], pileupEnergies_shadow_timeSlice[iter][timeSlice]);
      }

    auto g2Dir = scaleShadowDir->mkdir("g2Phase");
    g2Dir->cd();
      for (int g2phase = 0; g2phase < 8; ++g2phase)
      {
        auto g2_dir = g2Dir->mkdir(Form("g2%i",g2phase));
        g2_dir->cd();
        calcPileupScaleFactor(energies_g2phase[iter][1][g2phase], pileupEnergies_shadow_g2phase[iter][g2phase]);
      }

/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  auto efficiencyDir = pileupDirs[iter]->mkdir("Efficiency");
  efficiencyDir->cd();

    TF1* efficiencyFactor = new TF1("efficiencyFactor", "[0]", 0, histMaxTime);
    efficiencyFactor->SetLineColor(2);

  auto noScaleFactDir = efficiencyDir->mkdir("NoScaleFactorApplied");
  noScaleFactDir->cd();

    TH1F* PtimesDivDiff = (TH1F*) PtimesDiff->Clone("PtimesDivDiff");
    PtimesDivDiff->Divide(pileupTimes[iter]);
    PtimesDivDiff->SetTitle("Residual (Shadow - Truth) / Truth");
    PtimesDivDiff->GetYaxis()->SetTitle("Efficiency (decimal)");
    PtimesDivDiff->Fit(efficiencyFactor, "RQ");

    TH1F* PtimesDivDiff_threshold = (TH1F*) PtimesDiff_threshold->Clone("PtimesDivDiff_threshold");
    PtimesDivDiff_threshold->Divide(pileupTimes_threshold[iter]);
    PtimesDivDiff_threshold->SetTitle("Residual (Shadow - Truth) / Truth");
    PtimesDivDiff_threshold->GetYaxis()->SetTitle("Efficiency (decimal)");
    PtimesDivDiff_threshold->Fit(efficiencyFactor, "RQ");

    efficiencyFactor->SetRange(0, 10000);

    TH1F* PenergiesDivDiff = (TH1F*) PenergiesDiff->Clone("PenergiesDivDiff");
    PenergiesDivDiff->Divide(pileupEnergies[iter]);
    PenergiesDivDiff->SetTitle("Residual (Shadow - Truth) / Truth");
    PenergiesDivDiff->GetYaxis()->SetTitle("Efficiency (decimal)");
    PenergiesDivDiff->Fit(efficiencyFactor, "RQ");

  auto scaleFactDir = efficiencyDir->mkdir("ScaleFactorApplied");
  scaleFactDir->cd();

    efficiencyFactor->SetRange(0, histMaxTime);

    TH1F* PtimesDivDiff_scaled = (TH1F*) pileupTimes_shadow[iter]->Clone("PtimesDivDiff_scaled");
    PtimesDivDiff_scaled->Scale(scaleFactShadow);
    PtimesDivDiff_scaled->Add(pileupTimes[iter], -1);
    PtimesDivDiff_scaled->Divide(pileupTimes[iter]);
    PtimesDivDiff_scaled->SetTitle("Residual (Shadow - Truth) / Truth");
    PtimesDivDiff_scaled->GetYaxis()->SetTitle("Efficiency (decimal)");
    PtimesDivDiff_scaled->Fit(efficiencyFactor, "RQ");

    TH1F* PtimesDivDiff_threshold_scaled = (TH1F*) pileupTimes_shadow_threshold[iter]->Clone("PtimesDivDiff_threshold_scaled");
    PtimesDivDiff_threshold_scaled->Scale(scaleFactShadow);
    PtimesDivDiff_threshold_scaled->Add(pileupTimes_threshold[iter], -1);
    PtimesDivDiff_threshold_scaled->Divide(pileupTimes_threshold[iter]);
    PtimesDivDiff_threshold_scaled->SetTitle("Residual (Shadow - Truth) / Truth");
    PtimesDivDiff_threshold_scaled->GetYaxis()->SetTitle("Efficiency (decimal)");
    PtimesDivDiff_threshold_scaled->Fit(efficiencyFactor, "RQ");

    efficiencyFactor->SetRange(0, 10000);

    TH1F* PenergiesDivDiff_scaled = (TH1F*) pileupEnergies_shadow[iter]->Clone("PenergiesDivDiff_scaled");
    PenergiesDivDiff_scaled->Scale(scaleFactShadow);
    PenergiesDivDiff_scaled->Add(pileupEnergies[iter], -1);
    PenergiesDivDiff_scaled->Divide(pileupEnergies[iter]);
    PenergiesDivDiff_scaled->SetTitle("Residual (Shadow - Truth) / Truth");
    PenergiesDivDiff_scaled->GetYaxis()->SetTitle("Efficiency (decimal)");
    PenergiesDivDiff_scaled->Fit(efficiencyFactor, "RQ");

/////////////////////////////////////////////////////////////////////////////////////

    calcZetaEfficiency(efficiencyDir, energies_early[iter][1], energies_late[iter][1], pileupEnergies_shadow_early[iter], pileupEnergies_shadow_late[iter]);


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  // print out hit combinations
/*
  cout << "hitCombinations:" << endl;
  for(auto & val : hitCombinations){
    cout << val.first << " : " << val.second << endl;
  }
*/
/////////////////////////////////////////////////////////////////////////////////////


} // end iter loop

/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

  // this is just a place holder function so that I can use the same analysis code which gets this function for making truth pulls
  // the parameters don't have any direct meaning, but I can put in values that correspond to the original function for some
  // the R and lifetime parameters should still be able to be used in the pulls, the default phase should be slightly off, N I can just use the histogram entries I think (it will be close at least), and A I'm not sure about
  auto truthFunc = new TF1("truthFunc", "[0] + [1] + [2] + [3] + [4]", 0, histMaxTime);
  truthFunc->SetNpx(100);

  truthFunc->SetParameter(0, toyFiveParamHists[0][0]->GetEntries()*binWidth/defaultLifetime);
  truthFunc->SetParameter(1, defaultLifetime);
  truthFunc->SetParameter(2, 1); // random value
  truthFunc->SetParameter(3, 0);
  truthFunc->SetParameter(4, startingPhase);

  outputFile->cd();
  truthFunc->Write();

/////////////////////////////////////////////////////////////////////////////////////

  outputFile->Write();
  delete outputFile;

/////////////////////////////////////////////////////////////////////////////////////

  clock_t endTime = clock();

  cout << "Elapsed time: " << double(endTime - startTime)/CLOCKS_PER_SEC << " seconds" << endl;

return 1;

}
