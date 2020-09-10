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

#include "../ratioMacroHeaders/ratioAnalysisDefs.hh"
#include "../ratioMacroHeaders/pileupUtils.hh"

using namespace std;

  int deadTime = 5; // set to 0 to remove pileup generation
  int shadowGapTime = 10;

  int nPerFillPerCalo = 150;
  double nPts = 1e7;

// Flags for how much we're going to cheat
bool useTruthTriggerPulse = false;
bool useTruthShadowPulses = false;

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

// int makeToyHistFromTH2Pileup(string filePath)
int makeToyHistFromTH2Pileup()
{
  string filePath = "/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_04_00/srcs/gm2analyses/macros/RatioMacro/ToyMC/TimeEnergyPairsHist.root";

  clock_t startTime = clock();

  gRandom->SetSeed(2342345); // for TH2 GetRandom
  // TRandom3* randGenerator = new TRandom3(47326); // older rng which has its use commented out currently
  TRandom3* randFill = new TRandom3(40586444); // for randomized entries and entries per fill
  TRandom3* randUVPileup = new TRandom3(3325679); // for UV pileup filling


  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  TFile* outputFile = new TFile("histFromPairs.root","RECREATE");

/////////////////////////////////////////////////////////////////////////////////////

  TH2F* initialHist = (TH2F*) inputFile->Get("full2Dfunction_hist");
  // TH2F* initialHist = (TH2F*) inputFile->Get("constHist");
  // TH2F* initialHist = (TH2F*) inputFile->Get("expHist");

/////////////////////////////////////////////////////////////////////////////////////

  auto topDir = gFile->mkdir("topDir");

/////////////////////////////////////////////////////////////////////////////////////

  int totalIters = 4;

  TVectorD numIters(1);
  numIters[0] = totalIters;
  numIters.Write("Iters");

  TDirectory* histIterDirs[totalIters];

  TH1F* toyFiveParamHists[totalIters];
  TH1F* toyUHists[totalIters];
  TH1F* toyVHists[totalIters];
  TH1F* toyNumHists[totalIters];
  TH1F* toyDenomHists[totalIters];
  TH1F* energies[totalIters];

  TH1F* deltaT_hits[totalIters];     
  TH1F* deltaT_hits_threshold[totalIters];     

/////////////////////////////////////////////////////////////////////////////////////

  TF1* gainSagFunctions[totalIters];
  TF1* gainCorrFunctions[totalIters];

/////////////////////////////////////////////////////////////////////////////////////

  double gainSagAmplitudeMultiplier[totalIters]; 
  double gainSagLifetime[totalIters];

  double gm2PeriodPPMOffsets[totalIters];
  double gm2PeriodGuesses[totalIters];
  double energyThresholds[totalIters];

/////////////////////////////////////////////////////////////////////////////////////

    TRandom3* randIterUV[totalIters]; // for individual iteration UV randomization

/////////////////////////////////////////////////////////////////////////////////////

  for (int iter = 0; iter < totalIters; ++iter)
  {
    // amplitude
    gainSagAmplitudeMultiplier[iter] = 0;//1.;// - 1.*iter/totalIters; // realistic gain sags
    gainSagLifetime[iter] = 10000;
    // gainSagAmplitudeMultiplier[iter] = 1.*(iter-(totalIters/2.))/totalIters; // scan through 0 gain sag (negative to positive effect)
    // gainSagLifetime[iter] = 10000;

    // lifetime 
    // gainSagAmplitudeMultiplier[iter] = 1;
    // gainSagLifetime[iter] = 5000 + 2000*iter;


    gm2PeriodPPMOffsets[iter] = 0;
    gm2PeriodGuesses[iter] = g2Period * (1 + 1e-6 * gm2PeriodPPMOffsets[iter]);
    energyThresholds[iter] = defaultEThreshold;


    randIterUV[iter] = new TRandom3(31568);
  }

  gainSagAmplitudeMultiplier[0] = 0; // no gain sag for iteration 0
  gainSagLifetime[0] = 10000; // don't think it should matter what this number is exactly at long as the amplitude is 0

/////////////////////////////////////////////////////////////////////////////////////

    TH1::SetDefaultSumw2();

    double dT_maxtime = 100000;

  for (int iter = 0; iter < totalIters; ++iter)
  {
    histIterDirs[iter] = topDir->mkdir(Form("Iter%d", iter));
    histIterDirs[iter]->cd();

    toyFiveParamHists[iter] = new TH1F("Toy_5_Param_Hist","Toy_5_Param_Hist; time (ns); Events",nBins,0,histMaxTime);
    toyUHists[iter] = new TH1F("Toy_U_Hist","Toy_U_Hist; time (ns); Events",nBins,0,histMaxTime);
    toyVHists[iter] = new TH1F("Toy_V_Hist","Toy_V_Hist; time (ns); Events",nBins,0,histMaxTime);
    toyNumHists[iter] = new TH1F("Toy_Num_Hist","Toy_Num_Hist; time (ns); Events",nBins,0,histMaxTime);
    toyDenomHists[iter] = new TH1F("Toy_Denom_Hist","Toy_Denom_Hist; time (ns); Events",nBins,0,histMaxTime);

    energies[iter] = new TH1F("energies","energies; Energy (MeV); Events",700,0,7000);

    deltaT_hits[iter] = new TH1F("Time_Between_Hits_All", "Time_Between_Hits_All; #DeltaT (ns); Events", 10000, 0, dT_maxtime);
    deltaT_hits_threshold[iter] = new TH1F("Time_Between_Hits_Above_Threshold", "Time_Between_Hits_Above_Threshold; #DeltaT (ns); Events", 10000, 0, dT_maxtime);

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


    TVectorD parameterStore(2); // vector of doubles to store parameters used when creating histograms, for systematic studies
    // energy threshold, g-2 period used, possibly function parameters, etc.
    // should this just be a tntuple? or maybe a tlist with a variety of root objects? for tntuple see the ratioMacro code

    parameterStore[0] = gm2PeriodGuesses[iter];
    parameterStore[1] = energyThresholds[iter];
    parameterStore.Write("parameterStore");
  }

/////////////////////////////////////////////////////////////////////////////////////

    auto pileupDir = topDir->mkdir("PileupPlots");
    pileupDir->cd();

    TH1F* pileupTimes = new TH1F("pileupTimes", "pileupTimes; time (ns); Events", nBins, 0 , histMaxTime);
    TH1F* pileupTimes_threshold = new TH1F("pileupTimes_threshold", "pileupTimes_threshold; time (ns); Events", nBins, 0 , histMaxTime);

    TH1F* pileupEnergies = new TH1F("pileupEnergies", "pileupEnergies; Energy (MeV); Events", 700, 0, 7000);
    TH1F* pileupEnergies_threshold = new TH1F("pileupEnergies_threshold", "pileupEnergies_threshold; Energy (MeV); Events", 700, 0, 7000);
        
    TH1F* pileupTimes_U = new TH1F("pileupTimes_U", "pileupTimes_U; time (ns); Events", nBins, 0 , histMaxTime);
    TH1F* pileupTimes_V = new TH1F("pileupTimes_V", "pileupTimes_V; time (ns); Events", nBins, 0 , histMaxTime);

      TH1F* pileupTimes_shadow = new TH1F("pileupTimes_shadow", "pileupTimes_shadow; time (ns); Events", nBins, 0 , histMaxTime);
      TH1F* pileupTimes_shadow_threshold = new TH1F("pileupTimes_shadow_threshold", "pileupTimes_shadow_threshold; time (ns); Events", nBins, 0 , histMaxTime);

      TH1F* pileupEnergies_shadow = new TH1F("pileupEnergies_shadow", "pileupEnergies_shadow; Energy (MeV); Events", 700, 0, 7000);
      TH1F* pileupEnergies_shadow_threshold = new TH1F("pileupEnergies_shadow_threshold", "pileupEnergies_shadow_threshold; Energy (MeV); Events", 700, 0, 7000);
          
      TH1F* pileupTimes_U_shadow = new TH1F("pileupTimes_U_shadow", "pileupTimes_U_shadow; time (ns); Events", nBins, 0 , histMaxTime);
      TH1F* pileupTimes_V_shadow = new TH1F("pileupTimes_V_shadow", "pileupTimes_V_shadow; time (ns); Events", nBins, 0 , histMaxTime);

        pileupEnergies_shadow->SetLineColor(2);
        pileupTimes_shadow->SetLineColor(2);
        pileupEnergies_shadow_threshold->SetLineColor(2);
        pileupTimes_shadow_threshold->SetLineColor(2);
        pileupTimes_U_shadow->SetLineColor(2);
        pileupTimes_V_shadow->SetLineColor(2);


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  double randEntries = randFill->Poisson(nPts);

  int tenIncrement = 0.1 * randEntries;
  int countdown = tenIncrement;
  int tenPercent = 0;

/////////////////////////////////////////////////////////////////////////////////////

  std::vector<std::pair<double, double> > trueHits;
  std::vector<std::pair<double, double> > measuredHits;
    std::vector<std::vector<std::pair<double, double> > > measuredHitsTruth; // vector of hits for each measured hit containing the true hit - if a measured hit is a pileup pulse then this contains the true pulses that make that pileup, otherwise it's just the single pulse


    std::vector<std::pair<double, double> > trueNegativePileup;
    std::vector<std::pair<double, double> > truePositivePileup;

    std::vector<std::pair<double, double> > shadowNegativePileup;
    std::vector<std::pair<double, double> > shadowPositivePileup;

/////////////////////////////////////////////////////////////////////////////////////

  int numFills = int(randFill->Poisson(nPts/nPerFillPerCalo));

  for (int entry = 0; entry < numFills; ++entry)
  {
    double hitsInThisFill = randFill->Poisson(nPerFillPerCalo);


    double hit_time;
    double hit_energy_raw;

      // fill a vector of hits corresponding to a "fill"

      for (int j = 0; j < hitsInThisFill; ++j)
      {
        initialHist->GetRandom2(hit_time, hit_energy_raw);
        trueHits.emplace_back(hit_time, hit_energy_raw);
      }
      std::sort(trueHits.begin(), trueHits.end());

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

      // loop through the fill hits and create the pileup spectrum

      uint k = 0;
      while (k < trueHits.size())
      {
        std::vector<std::pair<double, double> > trueHitsInMeasHit = {trueHits.at(k)};

        double t1 = trueHits.at(k).first;
        double addedEnergy = trueHits.at(k).second;

        if(k != trueHits.size()-1 && abs(trueHits.at(k).first - trueHits.at(k+1).first) < deadTime)
        {
          trueNegativePileup.push_back(trueHits.at(k));

            uint kk = k+1;
            while(kk != trueHits.size() && (abs(trueHits.at(k).first - trueHits.at(kk).first) < deadTime))
            {
              addedEnergy += trueHits.at(kk).second;
              trueNegativePileup.push_back(trueHits.at(kk));
              trueHitsInMeasHit.push_back(trueHits.at(kk));
              kk++;
            }

          truePositivePileup.emplace_back(t1, addedEnergy);
          measuredHits.emplace_back(t1, addedEnergy);
          measuredHitsTruth.push_back(trueHitsInMeasHit);

          k = kk;
        }
        else{
          measuredHits.push_back(trueHits.at(k));
          measuredHitsTruth.push_back(trueHitsInMeasHit);
          k++;
        }
      } // end while


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

// construct pileup from shadow pulses in 'data'

      std::vector<int> shadowIndices; // record shadow pileup indices for a single pulse, to loop over and make shadow pileup spectrum

      uint j = 0;
      while (j < measuredHits.size()-1)
      {
        uint indSave = j;

        checkForShadowPulses(j, &measuredHits, shadowIndices, deadTime, shadowGapTime);

        if(shadowIndices.size() > 0)
        {

/////////////////////////////////////////////////////////////////////////////////////

          double triggerTime = useTruthTriggerPulse ? measuredHitsTruth.at(indSave).at(0).first : measuredHits.at(indSave).first;
          double addedEnergy = useTruthTriggerPulse ? measuredHitsTruth.at(indSave).at(0).second : measuredHits.at(indSave).second;
          shadowNegativePileup.emplace_back(triggerTime, addedEnergy);

/////////////////////////////////////////////////////////////////////////////////////
    // with the above, a second pulse in the trigger window is ignored
    // below it is added in, but that is considering it improperly since the pileup of the 1st pulse should ignore that second pulse
    // there is then an open question though as to how to estimate the pileup of that second pulse, for which a measured hit doesn't exist, and for which the shadow pulse might lie within the shadow window of both trigger pulses
/*
          double triggerTime;
          double addedEnergy = 0.;
          
          if(useTruthTriggerPulse){
            triggerTime = measuredHitsTruth.at(indSave).at(0).first;

              for (uint k = 0; k < measuredHitsTruth.at(indSave).size(); ++k)
              {
                auto& truthHit = measuredHitsTruth.at(indSave).at(k);
                shadowNegativePileup.push_back(truthHit);
                addedEnergy += truthHit.second;
              }
          }
          else{
            triggerTime = measuredHits.at(indSave).first;
            addedEnergy += measuredHits.at(indSave).second;
            shadowNegativePileup.push_back(measuredHits.at(indSave));
          }
*/
/////////////////////////////////////////////////////////////////////////////////////

          for (uint i = 0; i < shadowIndices.size(); ++i)
          {
            int shadowIndex = shadowIndices.at(i);

            if(useTruthShadowPulses)
            {
              for(uint k = 0; k < measuredHitsTruth.at(shadowIndex).size(); k++)
              {
                auto& truthHit = measuredHitsTruth.at(shadowIndex).at(k);
                
                if( (truthHit.first - triggerTime > shadowGapTime) && (truthHit.first - triggerTime < shadowGapTime + deadTime)){ // extra dead time check for making sure all shadow measuredHitsTruth are in shadow window
                  addedEnergy += truthHit.second;
                  shadowNegativePileup.emplace_back(truthHit);
                }
              }
            } 
            else {
              addedEnergy += measuredHits.at(shadowIndex).second;
              shadowNegativePileup.push_back(measuredHits.at(shadowIndex));
            }
          }

/////////////////////////////////////////////////////////////////////////////////////

          shadowPositivePileup.emplace_back(measuredHits.at(indSave).first, addedEnergy);
          shadowIndices.clear();
        }

      }  //end while

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

    // fill pileup and pileup shadow histograms

    // can make a separate loop over iterations for filling these pileup threshold histograms - which might be useful further down the line

    double halfPeriodGuess_temp = g2Period/2.;
    double totalChance_temp = exp(halfPeriodGuess_temp/defaultLifetime) + exp(-halfPeriodGuess_temp/defaultLifetime) + 2.;
    double percentChanceUPlus_temp = exp(halfPeriodGuess_temp/defaultLifetime) / totalChance_temp;
    double percentChanceUMinus_temp = exp(-halfPeriodGuess_temp/defaultLifetime) / totalChance_temp;

    double randNum_temp;

/////////////////////////////////////////////////////////////////////////////////////

    for (uint i = 0; i < trueNegativePileup.size(); ++i)
    {
      pileupTimes->Fill(trueNegativePileup.at(i).first, -1);
      pileupEnergies->Fill(trueNegativePileup.at(i).second, -1);
        if(trueNegativePileup.at(i).second > defaultEThreshold)
        {
          pileupTimes_threshold->Fill(trueNegativePileup.at(i).first, -1);
          pileupEnergies_threshold->Fill(trueNegativePileup.at(i).second, -1);

          randNum_temp = randUVPileup->Uniform();
            if     (randNum_temp < percentChanceUPlus_temp) pileupTimes_U->Fill(trueNegativePileup.at(i).first - halfPeriodGuess_temp, -1);
            else if(randNum_temp < (percentChanceUPlus_temp + percentChanceUMinus_temp)) pileupTimes_U->Fill(trueNegativePileup.at(i).first + halfPeriodGuess_temp, -1);
            else if(randNum_temp < 1) pileupTimes_V->Fill(trueNegativePileup.at(i).first, -1);

        }
    }
    for (uint i = 0; i < truePositivePileup.size(); ++i)
    {
      pileupTimes->Fill(truePositivePileup.at(i).first);
      pileupEnergies->Fill(truePositivePileup.at(i).second);
        if(truePositivePileup.at(i).second > defaultEThreshold)
        {
          pileupTimes_threshold->Fill(truePositivePileup.at(i).first);
          pileupEnergies_threshold->Fill(truePositivePileup.at(i).second);

          randNum_temp = randUVPileup->Uniform();
            if     (randNum_temp < percentChanceUPlus_temp) pileupTimes_U->Fill(truePositivePileup.at(i).first - halfPeriodGuess_temp);
            else if(randNum_temp < (percentChanceUPlus_temp + percentChanceUMinus_temp)) pileupTimes_U->Fill(truePositivePileup.at(i).first + halfPeriodGuess_temp);
            else if(randNum_temp < 1) pileupTimes_V->Fill(truePositivePileup.at(i).first);

        }
    }

/////////////////////////////////////////////////////////////////////////////////////

    for (uint i = 0; i < shadowNegativePileup.size(); ++i)
    {
      pileupTimes_shadow->Fill(shadowNegativePileup.at(i).first, -1);
      pileupEnergies_shadow->Fill(shadowNegativePileup.at(i).second, -1);
        if(shadowNegativePileup.at(i).second > defaultEThreshold)
        {
          pileupTimes_shadow_threshold->Fill(shadowNegativePileup.at(i).first, -1);
          pileupEnergies_shadow_threshold->Fill(shadowNegativePileup.at(i).second, -1);

          randNum_temp = randUVPileup->Uniform();
            if     (randNum_temp < percentChanceUPlus_temp) pileupTimes_U_shadow->Fill(shadowNegativePileup.at(i).first - halfPeriodGuess_temp, -1);
            else if(randNum_temp < (percentChanceUPlus_temp + percentChanceUMinus_temp)) pileupTimes_U_shadow->Fill(shadowNegativePileup.at(i).first + halfPeriodGuess_temp, -1);
            else if(randNum_temp < 1) pileupTimes_V_shadow->Fill(shadowNegativePileup.at(i).first, -1);

        }
    }
    for (uint i = 0; i < shadowPositivePileup.size(); ++i)
    {
      pileupTimes_shadow->Fill(shadowPositivePileup.at(i).first);
      pileupEnergies_shadow->Fill(shadowPositivePileup.at(i).second);
        if(shadowPositivePileup.at(i).second > defaultEThreshold)
        {
          pileupTimes_shadow_threshold->Fill(shadowPositivePileup.at(i).first);
          pileupEnergies_shadow_threshold->Fill(shadowPositivePileup.at(i).second);

          randNum_temp = randUVPileup->Uniform();
            if     (randNum_temp < percentChanceUPlus_temp) pileupTimes_U_shadow->Fill(shadowPositivePileup.at(i).first - halfPeriodGuess_temp);
            else if(randNum_temp < (percentChanceUPlus_temp + percentChanceUMinus_temp)) pileupTimes_U_shadow->Fill(shadowPositivePileup.at(i).first + halfPeriodGuess_temp);
            else if(randNum_temp < 1) pileupTimes_V_shadow->Fill(shadowPositivePileup.at(i).first);

        }
    }


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/*
    // comment this in for faster running, where the hits are the same for each iteration
      for (uint hitJ = 0; hitJ < measuredHits.size(); ++hitJ) // loop through the new hits and fill histograms
      {
          if(--countdown == 0) // progress output
          {
            tenPercent++;
            cout << tenPercent << "0%" << " completed" << endl;
            countdown = tenIncrement;
          }


        double time = measuredHits.at(hitJ).first;
        double energy_raw = measuredHits.at(hitJ).second;


          double randNum = randGenerator->Uniform(); // generate this random number above the iterations loop so histograms are filled consistently
          // could possibly just use multiple random number generators with the same seed that are only used when filling UV hists as it previously was
*/
/////////////////////////////////////////////////////////////////////////////////////

          // for (int iter = 0; iter < totalIters; ++iter)
          for (int iter = 0; iter < 2; ++iter) // only run for 2 iterations if I want to create one set of true hits histograms and one set of measured hits histograms (iters 2 and 3 are then filled down below)
          {
        
/////////////////////////////////////////////////////////////////////////////////////

  // comment this in for when the hits are different from iteration to iteration (pileup vs no pileup) - the run time will be slower though

      std::vector<std::pair<double, double> > myHits;

      if(iter == 0) myHits = trueHits;
      else myHits = measuredHits;

      for (uint hitJ = 0; hitJ < myHits.size(); ++hitJ) // loop through the new hits and fill histograms
      {
          if(--countdown == 0) // progress output
          {
            tenPercent++;
            cout << tenPercent << "0%" << " completed" << endl;
            countdown = tenIncrement;
          }

        double time = myHits.at(hitJ).first;
        double energy_raw = myHits.at(hitJ).second;

        double randNum = randIterUV[iter]->Uniform(); // if hits loop is below iters loop, then I need individual random number generators with the same seed so that separate iterations are filled identically, but also so that from one outer loop to the next the random number keeps changing

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

              energies[iter]->Fill(energy_measured);

              if(hitJ != myHits.size()-1) deltaT_hits[iter]->Fill(myHits.at(hitJ+1).first - time);

              if(energy_measured > energyThresholds[iter])
              {
                    toyFiveParamHists[iter]->Fill(time);

                    if(hitJ != myHits.size()-1) deltaT_hits_threshold[iter]->Fill(myHits.at(hitJ+1).first - time);

                    double halfPeriodGuess = gm2PeriodGuesses[iter]/2;

                    double totalChance = exp(halfPeriodGuess/defaultLifetime) + exp(-halfPeriodGuess/defaultLifetime) + 2.;
                    double percentChanceUPlus = exp(halfPeriodGuess/defaultLifetime) / totalChance;
                    double percentChanceUMinus = exp(-halfPeriodGuess/defaultLifetime) / totalChance;

                    if     (randNum < percentChanceUPlus) toyUHists[iter]->Fill(time - halfPeriodGuess); // careful with the signs here - the U+ hist is filled with pulses shifted by t - T/2, and weighted by e^+T/2tau, and vice versa for U-
                    else if(randNum < (percentChanceUPlus + percentChanceUMinus)) toyUHists[iter]->Fill(time + halfPeriodGuess);
                    else if(randNum < 1) toyVHists[iter]->Fill(time);
              }

          } // end sub loop over iterations (or hits depending on what's commented in)


      } // end sub loop over pileup events (or iterations loop depending on what's commented in)

      trueHits.clear();
      measuredHits.clear();
      measuredHitsTruth.clear();

      trueNegativePileup.clear();
      truePositivePileup.clear();

      shadowNegativePileup.clear();
      shadowPositivePileup.clear();

    } // end outer for loop over entries/nFillPerCalo


/////////////////////////////////////////////////////////////////////////////////////

    // this code for filling iterations 2 and 3 with the same hits as iteration 1, but to which I will apply different pileup corrections to
    // using this will skip any other iteration specific modifiers for these iterations

    for (int bin = 0; bin <= energies[1]->GetNbinsX()+1; bin++)
    {
      energies[2]->SetBinContent(bin, energies[1]->GetBinContent(bin));
      energies[2]->ResetStats();

      energies[3]->SetBinContent(bin, energies[1]->GetBinContent(bin));
      energies[3]->ResetStats();
    } 

    for (int bin = 0; bin <= toyFiveParamHists[1]->GetNbinsX()+1; bin++)
    {
      toyFiveParamHists[2]->SetBinContent(bin, toyFiveParamHists[1]->GetBinContent(bin));
      toyUHists[2]->SetBinContent(bin, toyUHists[1]->GetBinContent(bin));
      toyVHists[2]->SetBinContent(bin, toyVHists[1]->GetBinContent(bin));

      toyFiveParamHists[2]->ResetStats();
      toyUHists[2]->ResetStats();
      toyVHists[2]->ResetStats();

      toyFiveParamHists[3]->SetBinContent(bin, toyFiveParamHists[1]->GetBinContent(bin));
      toyUHists[3]->SetBinContent(bin, toyUHists[1]->GetBinContent(bin));
      toyVHists[3]->SetBinContent(bin, toyVHists[1]->GetBinContent(bin));

      toyFiveParamHists[3]->ResetStats();
      toyUHists[3]->ResetStats();
      toyVHists[3]->ResetStats();
    } 



/////////////////////////////////////////////////////////////////////////////////////
    // apply pileup correction

        toyFiveParamHists[2]->Add(pileupTimes_threshold, -1); // perfect correction
        energies[2]->Add(pileupEnergies, -1);
        toyUHists[2]->Add(pileupTimes_U, -1); // won't be the perfect correction for UV - because measured hits has a different size than the true hits leading to different random numbers when filling UV hists even though each iteration random number generator starts with the same seed
        toyVHists[2]->Add(pileupTimes_V, -1);

            toyFiveParamHists[3]->Add(pileupTimes_shadow_threshold, -1); // shadow correction            
            energies[3]->Add(pileupEnergies_shadow, -1);
            toyUHists[3]->Add(pileupTimes_U_shadow, -1);
            toyVHists[3]->Add(pileupTimes_V_shadow, -1);

/////////////////////////////////////////////////////////////////////////////////////

        for (int iter = 0; iter < totalIters; ++iter)
        {
          toyNumHists[iter]->Add(toyUHists[iter], toyVHists[iter], -1, 1);
          toyDenomHists[iter]->Add(toyUHists[iter], toyVHists[iter]);
        }


/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

  pileupDir->cd();

/////////////////////////////////////////////////////////////////////////////////////

  THStack* timesHistsStack = new THStack("timesHistsStack","MC Pileup Times");
  timesHistsStack->Add(pileupTimes);
  timesHistsStack->Add(pileupTimes_shadow);

  auto canvas_Ptimes = new TCanvas("canvas_Ptimes","canvas_Ptimes",200,10,1200,1000);
  TPad* padT1 = new TPad("padT1", "padT1", 0.02, 0.5, 0.98, 0.98);
  padT1->Draw();
  padT1->cd();
    timesHistsStack->Draw("nostack,hist");

    timesHistsStack->GetXaxis()->SetTitle("time (ns)");
    timesHistsStack->GetYaxis()->SetTitle("Events");
    padT1->Modified();

       auto tHS_legend = new TLegend(0.7,0.7,0.9,0.9);
       tHS_legend->AddEntry(pileupTimes,"True pileup","l");
       tHS_legend->AddEntry(pileupTimes_shadow,"Shadow pileup correction","l");
       tHS_legend->SetBorderSize(1);
       tHS_legend->Draw();

  canvas_Ptimes->cd();
   TPad* padT2 = new TPad("padT2", "padT2", 0.02, 0.25, 0.98, 0.5);
  padT2->Draw();
  padT2->cd();
    TH1F* PtimesDiff = (TH1F*) pileupTimes_shadow->Clone("PtimesDiff");
    PtimesDiff->Add(pileupTimes, -1);
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

  canvas_Ptimes->Write();
  delete canvas_Ptimes;

/////////////////////////////////////////////////////////////////////////////////////

  THStack* timesHistsStack_threshold = new THStack("timesHistsStack_threshold","MC Pileup Times Above Threshold");
  timesHistsStack_threshold->Add(pileupTimes_threshold);
  timesHistsStack_threshold->Add(pileupTimes_shadow_threshold);

  auto canvas_Ptimes_threshold = new TCanvas("canvas_Ptimes_threshold","canvas_Ptimes_threshold",200,10,1200,1000);
  TPad* padT1_threshold = new TPad("padT1_threshold", "padT1_threshold", 0.02, 0.5, 0.98, 0.98);
  padT1_threshold->Draw();
  padT1_threshold->cd();
    timesHistsStack_threshold->Draw("nostack,hist");

    timesHistsStack_threshold->GetXaxis()->SetTitle("time (ns)");
    timesHistsStack_threshold->GetYaxis()->SetTitle("Events");
    padT1_threshold->Modified();

       auto tHS_threshold_legend = new TLegend(0.7,0.7,0.9,0.9);
       tHS_threshold_legend->AddEntry(pileupTimes_threshold,"True pileup","l");
       tHS_threshold_legend->AddEntry(pileupTimes_shadow_threshold,"Shadow pileup correction","l");
       tHS_threshold_legend->SetBorderSize(1);
       tHS_threshold_legend->Draw();

  canvas_Ptimes_threshold->cd();
   TPad* padT2_threshold = new TPad("padT2_threshold", "padT2_threshold", 0.02, 0.25, 0.98, 0.5);
  padT2_threshold->Draw();
  padT2_threshold->cd();
    TH1F* PtimesDiff_threshold = (TH1F*) pileupTimes_shadow_threshold->Clone("PtimesDiff_threshold");
    PtimesDiff_threshold->Add(pileupTimes_threshold, -1);
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

  canvas_Ptimes_threshold->Write();
  delete canvas_Ptimes_threshold;

/////////////////////////////////////////////////////////////////////////////////////
  
/////////////////////////////////////////////////////////////////////////////////////

  THStack* energyHistsStack = new THStack("energyHistsStack","MC Pileup Energies");
  energyHistsStack->Add(pileupEnergies);
  energyHistsStack->Add(pileupEnergies_shadow);

  auto canvas_Penergies = new TCanvas("canvas_Penergies","canvas_Penergies",200,10,1200,1000);
  TPad* padE1 = new TPad("padE1", "padE1", 0.02, 0.5, 0.98, 0.98);
  padE1->Draw();
  padE1->cd();
    energyHistsStack->Draw("nostack,hist");
    
    energyHistsStack->GetXaxis()->SetTitle("Energy (MeV)");
    energyHistsStack->GetYaxis()->SetTitle("Events");
    padE1->Modified();

       auto eHS_legend = new TLegend(0.7,0.7,0.9,0.9);
       eHS_legend->AddEntry(pileupEnergies,"True pileup","l");
       eHS_legend->AddEntry(pileupEnergies_shadow,"Shadow pileup correction","l");
       eHS_legend->SetBorderSize(1);
       eHS_legend->Draw();

/*
    canvas_Penergies->Update();
    TPaveStats* teststats = (TPaveStats*) pileupEnergies_shadow->GetListOfFunctions()->FindObject("stats");

    cout << "test stats pointer check: " << teststats << endl;

    // TPaveStats* teststats = (TPaveStats*) energyHistsStack->GetHistogram()->GetListOfFunctions()->FindObject("stats");
          teststats->SetX1NDC(.7);
          teststats->SetX2NDC(.9);
          teststats->SetY1NDC(.66);
          teststats->SetY2NDC(.9);
    canvas_Penergies->Modified();
    // teststats->Draw();
    // energyHistsStack->GetHistogram()->GetPainter()->PaintStat(1,0);
*/

  canvas_Penergies->cd();
   TPad* padE2 = new TPad("padE2", "padE2", 0.02, 0.25, 0.98, 0.5);
  padE2->Draw();
  padE2->cd();
    TH1F* PenergiesDiff = (TH1F*) pileupEnergies_shadow->Clone("PenergiesDiff");
    PenergiesDiff->Add(pileupEnergies, -1);
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

  canvas_Penergies->Write();
  delete canvas_Penergies;

/////////////////////////////////////////////////////////////////////////////////////

  THStack* energyHistsStack_threshold = new THStack("energyHistsStack_threshold","MC Pileup Energies Above Threshold");
  energyHistsStack_threshold->Add(pileupEnergies_threshold);
  energyHistsStack_threshold->Add(pileupEnergies_shadow_threshold);

  auto canvas_Penergies_threshold = new TCanvas("canvas_Penergies_threshold","canvas_Penergies_threshold",200,10,1200,1000);
  TPad* padE1_threshold = new TPad("padE1_threshold", "padE1_threshold", 0.02, 0.5, 0.98, 0.98);
  padE1_threshold->Draw();
  padE1_threshold->cd();
    energyHistsStack_threshold->Draw("nostack,hist");

    energyHistsStack_threshold->GetXaxis()->SetTitle("Energy (MeV)");
    energyHistsStack_threshold->GetYaxis()->SetTitle("Events");
    padE1_threshold->Modified();

       auto eHS_threshold_legend = new TLegend(0.7,0.7,0.9,0.9);
       eHS_threshold_legend->AddEntry(pileupEnergies_threshold,"True pileup","l");
       eHS_threshold_legend->AddEntry(pileupEnergies_shadow_threshold,"Shadow pileup correction","l");
       eHS_threshold_legend->SetBorderSize(1);
       eHS_threshold_legend->Draw();


  canvas_Penergies_threshold->cd();
   TPad* padE2_threshold = new TPad("padE2_threshold", "padE2_threshold", 0.02, 0.25, 0.98, 0.5);
  padE2_threshold->Draw();
  padE2_threshold->cd();
    TH1F* PenergiesDiff_threshold = (TH1F*) pileupEnergies_shadow_threshold->Clone("PenergiesDiff_threshold");
    PenergiesDiff_threshold->Add(pileupEnergies_threshold, -1);
    PenergiesDiff_threshold->SetTitle("Residual (Shadow - Truth)");
    PenergiesDiff_threshold->Draw("hist");

  canvas_Penergies_threshold->cd();
   TPad* padE3_threshold = new TPad("padE3_threshold", "padE3_threshold", 0.02, 0.02, 0.98, 0.25);
  padE3_threshold->Draw();
  padE3_threshold->cd();
    TH1F* PenergiesSigDiff_threshold = (TH1F*) PenergiesDiff_threshold->Clone("PenergiesSigDiff_threshold");
      for (int bin = 1; bin <= PenergiesSigDiff_threshold->GetNbinsX(); ++bin) if(PenergiesSigDiff_threshold->GetBinError(bin) != 0) PenergiesSigDiff_threshold->SetBinContent(bin, PenergiesSigDiff_threshold->GetBinContent(bin)/PenergiesSigDiff_threshold->GetBinError(bin));
    PenergiesSigDiff_threshold->SetTitle("Residual (Shadow - Truth) / Bin Error");
    PenergiesSigDiff_threshold->Draw("hist");

  canvas_Penergies_threshold->Write();
  delete canvas_Penergies_threshold;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
  
  TF1* expFunc = new TF1("expFunc", "[0]*exp(-x/[1])", 0, dT_maxtime);

  expFunc->SetParameter(1, histMaxTime/nPerFillPerCalo);
  cout << "Lifetime for deltaT hits is: " << histMaxTime/nPerFillPerCalo << endl;

  for (int iter = 0; iter < 2; ++iter)
  {
    expFunc->SetParameter(0, 1);
    double par0Start = deltaT_hits[iter]->GetEntries() / ( expFunc->Integral(0, dT_maxtime) / deltaT_hits[iter]->GetBinWidth(1) );
    expFunc->SetParameter(0, par0Start);

      deltaT_hits[iter]->Fit(expFunc, "QRLI");
      deltaT_hits[iter]->GetFunction("expFunc")->SetLineColor(2);

    expFunc->SetParameter(0, 1);
    par0Start = deltaT_hits_threshold[iter]->GetEntries() / ( expFunc->Integral(0, dT_maxtime) / deltaT_hits_threshold[iter]->GetBinWidth(1) );
    expFunc->SetParameter(0, par0Start);

      deltaT_hits_threshold[iter]->Fit(expFunc, "QRLI");
      deltaT_hits_threshold[iter]->GetFunction("expFunc")->SetLineColor(2);
  }

/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

  // this should maybe be moved up if I want to apply this scale factor to the pileup instead of just calculate it - in that case would only want to apply it to iteration 3 or maybe even make another iteration 

  auto scaleDir = pileupDir->mkdir("Scale");
  auto scaleTruthDir = scaleDir->mkdir("TruePileup");
  auto scaleShadowDir = scaleDir->mkdir("ShadowPileup");

  scaleTruthDir->cd();
  double scaleFactTruth = calcPileupScaleFactor(energies[1], pileupEnergies_threshold);
  cout << "scale factor from truth should be 1 and is : " << scaleFactTruth << endl;

  scaleShadowDir->cd();
  double scaleFactShadow = calcPileupScaleFactor(energies[1], pileupEnergies_shadow_threshold);
  cout << "scale factor from measured should hopefully be 1 and is : " << scaleFactShadow << endl;


/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

  // this is just a place holder function so that I can use the same analysis code which gets this function for making truth pulls
  // the parameters don't have any direct meaning, but I can put in values that correspond to the original function for some
  // the R and lifetime parameters should still be able to be used in the pulls, the default phase should be slightly off, N I can just use the histogram entries I think (it will be close at least), and A I'm not sure about
  auto truthFunc = new TF1("truthFunc", "[0] + [1] + [2] + [3] + [4]", 0, histMaxTime);
  truthFunc->SetNpx(100);

  truthFunc->SetParameter(0, toyFiveParamHists[0]->GetEntries()*binWidth/defaultLifetime);
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
