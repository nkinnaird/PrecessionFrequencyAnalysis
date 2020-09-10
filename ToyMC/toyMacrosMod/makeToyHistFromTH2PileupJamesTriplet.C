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

#include "ratioMacroHeaders/ratioAnalysisDefs.hh"
#include "ratioMacroHeaders/pileupUtils.hh"

using namespace std;

int deadTime = 5; // set to 0 to remove pileup generation
int shadowGapTime = 10;

/////////////////////////////////////////////////////////////////////////////////////

// Flags for how much we're going to cheat
bool useTruthTriggerPulse = false;
bool useTruthShadowPulses = false;

double nPts = 1e10;
int nPerFillPerCalo = 500;

// This draws randomly from the file if it's given or uniform in both if not
string filePath = "NONE";
//string filePath = "ConstantTimeEnergyPairsHist.root"; 
//string filePath = "/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_04_00/srcs/gm2analyses/macros/RatioMacro/ToyMC/TimeEnergyPairsHist.root";

// Choose single energy only (each particle has 550 MeV) - only valid with filepath "NONE"
bool monoEnergetic = false;

/////////////////////////////////////////////////////////////////////////////////////

int makeToyHistFromTH2PileupJamesTriplet(){

  if(shadowGapTime < deadTime) {
    cout << "SHADOW GAP IS SMALLER THAN DEADTIME - THIS WILL NOT WORK!" << endl;
    return -1;
  }
  // else if(shadowGapTime != deadTime){
  //cout << "SHADOW WINDOW DOESN'T START STRAIGHT AFTER DEADTIME - THIS CAN CAUSE PROBLEMS." << endl;
  //return -1;
  //}

  TRandom3* randGenerator = new TRandom3(3297);

  bool readInputFile = (filePath != "NONE");
  TFile* inputFile = 0;
  TH2F* initialHist = 0;
  if( readInputFile ) {
    if(monoEnergetic) {
      cout << "monoEnergetic flag doesn't work with filePath " << filePath << ". Quitting..." << endl;
      return -1;
    }
    inputFile = new TFile(filePath.c_str(),"READ");

    if(filePath == "ConstantTimeEnergyPairsHist.root") {
      initialHist = (TH2F*) inputFile->Get("constHist");
    } else {
      initialHist = (TH2F*) inputFile->Get("full2Dfunction_hist");
    }
  }

  TFile* outputFile = new TFile("histFromPairs.root","RECREATE");
  auto topDir = gFile->mkdir("topDir");

  /////////////////////////////////////////////////////////////////////////////////////

  TH1::SetDefaultSumw2();
  
  TDirectory* histIterDirs = topDir->mkdir(Form("Iter%d", 0));
  histIterDirs->cd();
  
  TH1D* energies = new TH1D("energies","energies; Energy (MeV); Events",1000,0,10000);
  TH1D* times = new TH1D("times","times; Time (ns); Events", nBins, 0 , histMaxTime);

  TH1D* trueEnergies = new TH1D("trueEnergies","trueEnergies; Energy (MeV); Events",1000,0,10000);
  TH1D* trueTimes = new TH1D("trueTimes","trueTimes; Time (ns); Events", nBins, 0 , histMaxTime);

  double dT_maxtime = 10000;
  TH1D* deltaT_true = new TH1D("Time_Between_Hits_True", "True_Time_Between_Hits; #DeltaT (ns); Events", 10000, 0, dT_maxtime);
  TH1D* deltaT_hits = new TH1D("Time_Between_Hits", "Time_Between_Hits; #DeltaT (ns); Events", 10000, 0, dT_maxtime);

  /////////////////////////////////////////////////////////////////////////////////////

  auto pileupDir = topDir->mkdir("PileupPlots");
  pileupDir->cd();

  TH1D* truePileupTimes = new TH1D("truePileupTimes", "truePileupTimes; time (ns); Events", nBins, 0 , histMaxTime);
  TH1D* truePileupEnergies = new TH1D("truePileupEnergies", "truePileupEnergies; Energy (MeV); Events", 1000,0,10000);
  TH1D* truePileupHits = new TH1D("truePileupHits", "truePileupHits;Hits; Events", 10, -0.5, 9.5);
  truePileupTimes->SetLineWidth(2);
  truePileupEnergies->SetLineWidth(2);
  truePileupHits->SetLineWidth(2);

  TH1D* pileupTimes_shadow = new TH1D("pileupTimes_shadow", "pileupTimes_shadow; time (ns); Events", nBins, 0 , histMaxTime);
  TH1D* pileupEnergies_shadow = new TH1D("pileupEnergies_shadow", "pileupEnergies_shadow; Energy (MeV); Events", 1000,0,10000);
  TH1D* pileupHits_shadow = new TH1D("pileupHits_shadow", "pileupHits_shadow;Hits; Events", 10, -0.5, 9.5);
  pileupTimes_shadow->SetLineColor(2);
  pileupEnergies_shadow->SetLineColor(2);
  pileupHits_shadow->SetLineColor(2);

  TH1D* pileupTripletTimes_shadow = new TH1D("pileupTripletTimes_shadow", "pileupTripletTimes_shadow; time (ns); Events", nBins, 0 , histMaxTime);
  TH1D* pileupTripletEnergies_shadow = new TH1D("pileupTripletEnergies_shadow", "pileupTripletEnergies_shadow; Energy (MeV); Events", 1000,0,10000);
  TH1D* pileupTripletHits_shadow = new TH1D("pileupTripletHits_shadow", "pileupTripletHits_shadow;Hits; Events", 10, -0.5, 9.5);
  pileupTripletTimes_shadow->SetLineColor(4);
  pileupTripletEnergies_shadow->SetLineColor(4);
  pileupTripletHits_shadow->SetLineColor(4);

  /////////////////////////////////////////////////////////////////////////////////////

  // Set number of positrons per fill and then calculate number of fills to get to nPts
  long long numFills = nPts/nPerFillPerCalo;

  std::vector<std::pair<double, double> > trueHits;
  std::vector<std::pair<double, double> > measuredHits;
  std::vector<std::vector<std::pair<double, double> > > measuredHitsTruth;

  std::vector<std::pair<double, double> > trueNegativePileup;
  std::vector<std::pair<double, double> > truePositivePileup;
  std::vector<int> trueNPileUp;

  std::vector<std::pair<double, double> > shadowNegativePileup;
  std::vector<std::pair<double, double> > shadowPositivePileup;
  std::vector<int> shadowNPileUp;

  std::vector<std::pair<double, double> > shadowNegativePileupTriplet;
  std::vector<std::pair<double, double> > shadowPositivePileupTriplet;

  std::map<int,long long> hitCombinations;
  std::map<int,long long> hitCombinationsTrueWindow;

  /////////////////////////////////////////////////////////////////////////////////////

  double lastTrueTimeLastFill = 0;
  double lastTimeLastFill = 0;
  
  double targetPercentage = 0;

  for (long long entry = 0; entry < numFills; ++entry) {

    if(100*double(entry) / numFills > targetPercentage){
      cout << "Processed " << Form("%.1f",100*double(entry)/numFills) << "%..." << endl;
      targetPercentage += 0.1;
    }

    double hitsInThisFill = randGenerator->Poisson(nPerFillPerCalo);

    double hit_time, hit_energy_raw;

    // fill a vector of hits corresponding to a "fill"
    for (int j = 0; j < hitsInThisFill; ++j){

      if( readInputFile ) {
        initialHist->GetRandom2(hit_time, hit_energy_raw);
      }  else {
	hit_time = (histMaxTime+300)*randGenerator->Uniform()-150; // Generate these before & after the timing window
	hit_energy_raw = monoEnergetic ? 550 : 3100*randGenerator->Uniform();
      }

      trueHits.emplace_back(hit_time, hit_energy_raw);
    }
    std::sort(trueHits.begin(), trueHits.end());

    for(uint hit = 0; hit < trueHits.size(); hit++){
      trueTimes->Fill(trueHits.at(hit).first);
      trueEnergies->Fill(trueHits.at(hit).second);
      if(hit < trueHits.size() - 1) deltaT_true->Fill(trueHits.at(hit+1).first - trueHits.at(hit).first);
      if(hit == 0) deltaT_true->Fill(trueHits.at(hit).first + (histMaxTime - lastTrueTimeLastFill));
      if(hit == trueHits.size() - 1) lastTrueTimeLastFill = trueHits.at(hit).first;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////
   
    // loop through the fill hits and create the pileup spectrum
    uint k = 0;
    while (k < trueHits.size()) {

      std::vector<std::pair<double, double> > trueHitsInMeasHit = {trueHits.at(k)};

      int nTrueHits = 0;

      //      int nPileUp = 1;
      double t1 = trueHits.at(k).first;
      //      double addedTime = trueHits.at(k).first;
      double addedEnergy = trueHits.at(k).second;

      if(k != trueHits.size()-1 && abs(trueHits.at(k).first - trueHits.at(k+1).first) < deadTime){

	trueNegativePileup.push_back(trueHits.at(k));
	nTrueHits++;

	uint kk = k+1;
	while(kk != trueHits.size() && (abs(trueHits.at(kk-1).first - trueHits.at(kk).first) < deadTime)){
	  //	  addedTime += trueHits.at(kk).first;
	  addedEnergy += trueHits.at(kk).second;
	  //  nPileUp++;
	  trueNegativePileup.push_back(trueHits.at(kk));
	  nTrueHits++;
	  trueHitsInMeasHit.push_back(trueHits.at(kk));
	  kk++;

	}
	
	//	double meanTime = addedTime / nPileUp;	
	//	truePositivePileup.emplace_back(meanTime, addedEnergy);
	//	measuredHits.emplace_back(meanTime, addedEnergy);

	trueNPileUp.push_back(nTrueHits);
	truePositivePileup.emplace_back(t1, addedEnergy);

	measuredHits.emplace_back(t1, addedEnergy);
	measuredHitsTruth.push_back(trueHitsInMeasHit);
	
       	k = kk;
      } else {
	measuredHits.push_back(trueHits.at(k));
	measuredHitsTruth.push_back(trueHitsInMeasHit);
	k++;
      }
    } // end while

    /////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////
    
    // construct pileup from shadow pulses in 'data'
    std::vector<int> shadowIndices;
    std::vector<int> doubleShadowIndices;

    double maxCheckTime = 1e8;//65000;
    uint j = 0;
    while (j < measuredHits.size()-1) {

      uint indSave = j;

      if(measuredHits.at(indSave).first > maxCheckTime){
	j++;
	continue;
      }

      checkForShadowPulsesTriplet(j, &measuredHits, shadowIndices, doubleShadowIndices, deadTime, shadowGapTime);

      int nTriggerHits = measuredHitsTruth.at(indSave).size();
      int nShadowHits = 0;
      int nShadowHits2 = 0;
      int nShadowHitsTrueWindow = 0;
      int nShadowHitsTrueWindow2 = 0;

      if(shadowIndices.size() > 0) {

	int nMeasHits = 0;

	double triggerTime = useTruthTriggerPulse ? measuredHitsTruth.at(indSave).at(0).first : measuredHits.at(indSave).first;
	double triggerEnergy = useTruthTriggerPulse ? measuredHitsTruth.at(indSave).at(0).second : measuredHits.at(indSave).second;
	shadowNegativePileup.emplace_back(triggerTime, triggerEnergy);
	nMeasHits++;
	//	cout << "-VE Meas:" << j << "=(" << triggerTime << ":" << triggerEnergy << ")" << endl;

	double addedEnergyDoublet = 0;
    	for (uint i = 0; i < shadowIndices.size(); ++i) {

    	  int shadowIndex = shadowIndices.at(i);

	  if(shadowIndices.size() != 1) {
	    std::cout << "shadowIndices.size() == " << shadowIndices.size() << endl;
	    return -1;
	  }

	  nShadowHits += measuredHitsTruth.at(shadowIndex).size();
	  for(uint hit = 0; hit < measuredHitsTruth.at(shadowIndex).size(); hit++){
	    auto& truthHit = measuredHitsTruth.at(shadowIndex).at(hit);
	    if(truthHit.first - triggerTime > shadowGapTime and truthHit.first - triggerTime < shadowGapTime + deadTime){
	      nShadowHitsTrueWindow++;
	    }
	  }

	  if(useTruthShadowPulses){
	    for(uint hit = 0; hit < measuredHitsTruth.at(shadowIndex).size(); hit++){
	      auto& truthHit = measuredHitsTruth.at(shadowIndex).at(hit);
	      if(truthHit.first - triggerTime > shadowGapTime and truthHit.first - triggerTime < shadowGapTime + deadTime){
		addedEnergyDoublet += truthHit.second;
		shadowNegativePileup.emplace_back(truthHit.first, truthHit.second);
		nMeasHits++;
		//	cout << "-VE Meas:" << j << "=(" << truthHit.first << ":" << truthHit.second << endl;
	      }
	    }
	  } else {
	    addedEnergyDoublet += measuredHits.at(shadowIndex).second;
	    shadowNegativePileup.push_back(measuredHits.at(shadowIndex));
	    nMeasHits++;
	    // cout << "-VE Meas:" << j << "=(" << measuredHits.at(shadowIndex).first << ":" << measuredHits.at(shadowIndex).second << ")" << endl;
	  }
    	}
	
	shadowNPileUp.push_back(nMeasHits);
	shadowPositivePileup.emplace_back(measuredHits.at(indSave).first, triggerEnergy + addedEnergyDoublet);
	//	cout << "+VE Meas:" << j << "=(" << measuredHits.at(indSave).first << ":" << triggerEnergy + addedEnergyDoublet << ")" << endl;
	//cout<< "========" << endl;
    
	// Try out triplets
	if(doubleShadowIndices.size() > 0) {

	  // Add up all the energy in the second shadow window
	  double addedEnergyTriplet = 0;
	  for (uint i = 0; i < doubleShadowIndices.size(); ++i) { addedEnergyTriplet += measuredHits.at(doubleShadowIndices.at(i)).second; }

	  // See combinatorix spreadsheet for these numbers - also found empirically
	  for(int n = 0; n < 1; n++) shadowNegativePileupTriplet.emplace_back(triggerTime, triggerEnergy + addedEnergyDoublet + addedEnergyTriplet);

	  for(int n = 0; n < 2; n++) shadowPositivePileupTriplet.emplace_back(triggerTime, triggerEnergy + addedEnergyDoublet);
	  for(int n = 0; n < 2; n++) shadowPositivePileupTriplet.emplace_back(triggerTime, addedEnergyDoublet + addedEnergyTriplet);

	  for(int n = 0; n < 1; n++) shadowNegativePileupTriplet.emplace_back(triggerTime, triggerEnergy);
	  for(int n = 0; n < 2; n++) shadowNegativePileupTriplet.emplace_back(triggerTime, addedEnergyDoublet);
	  for(int n = 0; n < 1; n++) shadowNegativePileupTriplet.emplace_back(triggerTime, addedEnergyTriplet);

  	  // Put same in non-test version
	  for(int n = 0; n < 1; n++) shadowNegativePileup.emplace_back(triggerTime, triggerEnergy + addedEnergyDoublet + addedEnergyTriplet);

	  for(int n = 0; n < 2; n++) shadowPositivePileup.emplace_back(triggerTime, triggerEnergy + addedEnergyDoublet);
	  for(int n = 0; n < 2; n++) shadowPositivePileup.emplace_back(triggerTime, addedEnergyDoublet + addedEnergyTriplet);

	  for(int n = 0; n < 1; n++) shadowNegativePileup.emplace_back(triggerTime, triggerEnergy);
	  for(int n = 0; n < 2; n++) shadowNegativePileup.emplace_back(triggerTime, addedEnergyDoublet);
	  for(int n = 0; n < 1; n++) shadowNegativePileup.emplace_back(triggerTime, addedEnergyTriplet);

	}

      }

      double triggerTime = useTruthTriggerPulse ? measuredHitsTruth.at(indSave).at(0).first : measuredHits.at(indSave).first;
      for (uint i = 0; i < doubleShadowIndices.size(); ++i) {
	int shadowIndex = doubleShadowIndices.at(i);
	nShadowHits2 += measuredHitsTruth.at(shadowIndex).size();
	for(uint hit = 0; hit < measuredHitsTruth.at(shadowIndex).size(); hit++){
	  auto& truthHit = measuredHitsTruth.at(shadowIndex).at(hit);
	  if(truthHit.first - triggerTime > 2*shadowGapTime and truthHit.first - triggerTime < 2*shadowGapTime + deadTime){
	    nShadowHitsTrueWindow2++;
	  }
	}
      }


      shadowIndices.clear();
      doubleShadowIndices.clear();

      // Increment maps
      int hitCombo = nTriggerHits*100 + nShadowHits*10 + nShadowHits2;
      int hitComboTrueWindow = nTriggerHits*100 + nShadowHitsTrueWindow*10 + nShadowHitsTrueWindow2;
      
      hitCombinations[hitCombo]++;
      hitCombinationsTrueWindow[hitComboTrueWindow]++;

    }  //end while

    //    cout << "shadowPulses = " << nShadowPulses << "\t2nd Shadow Pulses = "<< nDoubleShadowPulses << "\tBoth = " << nShadowDoubleShadowPulses << endl;

    if(false){
      cout << "Meas Times = {";
      for (uint hitJ = 0; hitJ < measuredHits.size(); ++hitJ){
	if(measuredHits.at(hitJ).first > maxCheckTime) break;
	cout << hitJ << ":" << measuredHits.at(hitJ).first << ", ";
      }
      cout << "}" << endl;
      cout << "True Times = {";
      for (uint hitJ = 0; hitJ < trueHits.size(); ++hitJ) {
	if(trueHits.at(hitJ).first > maxCheckTime) break;
	cout << hitJ << ":" << trueHits.at(hitJ).first << ", ";
      }
      cout << "}" << endl;
      return -1;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////
    
    // fill pileup and pileup shadow histograms
    for (uint i = 0; i < trueNegativePileup.size(); ++i) {
      if(trueNegativePileup.at(i).first < 0 or trueNegativePileup.at(i).first > histMaxTime) continue;
      truePileupTimes->Fill(trueNegativePileup.at(i).first, -1);
      truePileupEnergies->Fill(trueNegativePileup.at(i).second, -1);
      truePileupHits->Fill(1, -1);
    }

    for (uint i = 0; i < truePositivePileup.size(); ++i) {
      if(truePositivePileup.at(i).first < 0 or truePositivePileup.at(i).first > histMaxTime) continue;
      truePileupTimes->Fill(truePositivePileup.at(i).first);
      truePileupEnergies->Fill(truePositivePileup.at(i).second);
      //      truePileupHits->Fill(trueNPileUp.at(i));
    }
    /////////////////////////////////////////////////////////////////////////////////////

    for (uint i = 0; i < shadowNegativePileup.size(); ++i) {
      if(shadowNegativePileup.at(i).first < 0 or shadowNegativePileup.at(i).first > histMaxTime) continue;
      pileupTimes_shadow->Fill(shadowNegativePileup.at(i).first, -1);
      pileupEnergies_shadow->Fill(shadowNegativePileup.at(i).second, -1);
      pileupHits_shadow->Fill(1, -1);
    }
    for (uint i = 0; i < shadowPositivePileup.size(); ++i) {
      if(shadowPositivePileup.at(i).first < 0 or shadowPositivePileup.at(i).first > histMaxTime) continue;
      pileupTimes_shadow->Fill(shadowPositivePileup.at(i).first);
      pileupEnergies_shadow->Fill(shadowPositivePileup.at(i).second);
      //      pileupHits_shadow->Fill(shadowNPileUp.at(i));
    }

    for (uint i = 0; i < shadowNegativePileupTriplet.size(); ++i) {
      if(shadowNegativePileupTriplet.at(i).first < 0 or shadowNegativePileupTriplet.at(i).first > histMaxTime) continue;
      pileupTripletEnergies_shadow->Fill(shadowNegativePileupTriplet.at(i).second, -1);
    }
    for (uint i = 0; i < shadowPositivePileupTriplet.size(); ++i) {
      if(shadowPositivePileupTriplet.at(i).first < 0 or shadowPositivePileupTriplet.at(i).first > histMaxTime) continue;
      pileupTripletEnergies_shadow->Fill(shadowPositivePileupTriplet.at(i).second);
    }

    
    /////////////////////////////////////////////////////////////////////////////////////

    // comment this in for faster running, where the hits are the same for each iteration
    for (uint hitJ = 0; hitJ < measuredHits.size(); ++hitJ){

      double time = measuredHits.at(hitJ).first;
      double energy = measuredHits.at(hitJ).second;
    
      if(time < 0 or time > histMaxTime) continue;

      energies->Fill(energy);
      times->Fill(time);

      if(hitJ != measuredHits.size()-1) deltaT_hits->Fill(measuredHits.at(hitJ+1).first - time);
      if(hitJ == 0) deltaT_hits->Fill(time + (histMaxTime - lastTimeLastFill));
      if(hitJ == measuredHits.size() - 1) lastTimeLastFill = time;

    } // end sub loop over pileup events (or iterations loop depending on what's commented in)

    trueHits.clear();
    measuredHits.clear();
    measuredHitsTruth.clear();
    trueNegativePileup.clear();
    truePositivePileup.clear();
    trueNPileUp.clear();
    
    shadowNegativePileup.clear();
    shadowPositivePileup.clear();
    shadowNegativePileupTriplet.clear();
    shadowPositivePileupTriplet.clear();
    shadowNPileUp.clear();

  } // end outer for loop over entries/nFillPerCalo

  cout << "hitCombinations:" << endl;
  for(auto & val : hitCombinations){
    cout << val.first << " : " << val.second << endl;
  }

  cout << "hitCombinationsTrueWindow:" << endl;
  for(auto & val : hitCombinationsTrueWindow){
    cout << val.first << " : " << val.second << endl;
  }

  // Make sure truth is correct by adding back truth to measured and making sure we recover the right thing!
  histIterDirs->cd();
  auto corrEnergiesTruth = (TH1D*)energies->Clone("corrEnergiesTruth");
  corrEnergiesTruth->Add(truePileupEnergies, -1);

  auto corrTimesTruth = (TH1D*)times->Clone("corrTimesTruth");
  corrTimesTruth->Add(truePileupTimes, -1);

  auto corrEnergies = (TH1D*)energies->Clone("corrEnergies");
  corrEnergies->Add(pileupEnergies_shadow, -1);

  auto corrTimes = (TH1D*)times->Clone("corrTimes");
  corrTimes->Add(pileupTimes_shadow, -1);
  
  pileupDir->cd();


  // Figure out exactly what we need to add to pileup to get this right
  TH1D* pileupEnergyDiff = (TH1D*)truePileupEnergies->Clone("pileupEnergyDiff");
  pileupEnergyDiff->Add(pileupEnergies_shadow,-1);

  TCanvas* cEnergy = new TCanvas();
  energies->Draw();

  TCanvas* cEnergyDiff = new TCanvas();
  pileupEnergyDiff->Draw("HIST");
  pileupTripletEnergies_shadow->Draw("SAME HIST");
  cEnergyDiff->SaveAs("PileUpDiffTmp.C");

  TCanvas* cEnergyDiv = new TCanvas();
  TH1D* pileupEnergyDiffDiv = (TH1D*)pileupEnergyDiff->Clone("pileupEnergyDiffDiv");  
  pileupEnergyDiffDiv->Divide(pileupTripletEnergies_shadow);
  pileupEnergyDiffDiv->Draw();

  /////////////////////////////////////////////////////////////////////////////////////

  if(false){
  THStack* timesHistsStack = new THStack("timesHistsStack","MC Pileup Times");
  timesHistsStack->Add(truePileupTimes);
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
  tHS_legend->AddEntry(truePileupTimes,"True pileup","l");
  tHS_legend->AddEntry(pileupTimes_shadow,"Shadow pileup correction","l");
  tHS_legend->SetBorderSize(1);
  tHS_legend->Draw();

  canvas_Ptimes->cd();
  TPad* padT2 = new TPad("padT2", "padT2", 0.02, 0.25, 0.98, 0.5);
  padT2->Draw();
  padT2->cd();
  TH1D* PtimesDiff = (TH1D*) pileupTimes_shadow->Clone("PtimesDiff");
  PtimesDiff->Add(truePileupTimes, -1);
  PtimesDiff->SetTitle("Residual (Shadow - Truth)");
  PtimesDiff->Draw("hist");

  canvas_Ptimes->cd();
  TPad* padT3 = new TPad("padT3", "padT3", 0.02, 0.02, 0.98, 0.25);
  padT3->Draw();
  padT3->cd();
  TH1D* PtimesSigDiff = (TH1D*) PtimesDiff->Clone("PtimesSigDiff");
  for (int bin = 1; bin <= PtimesSigDiff->GetNbinsX(); ++bin) if(PtimesSigDiff->GetBinError(bin) != 0) PtimesSigDiff->SetBinContent(bin, PtimesSigDiff->GetBinContent(bin)/PtimesSigDiff->GetBinError(bin));
  PtimesSigDiff->SetTitle("Residual (Shadow - Truth) / Bin Error");
  PtimesSigDiff->Draw("hist");

  canvas_Ptimes->Write();
  //  delete canvas_Ptimes;

  /////////////////////////////////////////////////////////////////////////////////////

  THStack* hitsHistsStack = new THStack("hitsHistsStack","MC Pileup Hits");
  hitsHistsStack->Add(truePileupHits);
  hitsHistsStack->Add(pileupHits_shadow);

  auto canvas_Phits = new TCanvas("canvas_Phits","canvas_Phits",200,10,1200,1000);
  TPad* padH1 = new TPad("padH1", "padH1", 0.02, 0.5, 0.98, 0.98);
  padH1->Draw();
  padH1->cd();
  hitsHistsStack->Draw("nostack,hist");
    
  hitsHistsStack->GetXaxis()->SetTitle("Hits");
  hitsHistsStack->GetYaxis()->SetTitle("Events");
  padH1->Modified();

  auto hHS_legend = new TLegend(0.7,0.7,0.9,0.9);
  hHS_legend->AddEntry(truePileupHits,"True pileup","l");
  hHS_legend->AddEntry(pileupHits_shadow,"Shadow pileup correction","l");
  hHS_legend->SetBorderSize(1);
  hHS_legend->Draw();

  canvas_Phits->cd();
  TPad* padH2 = new TPad("padH2", "padH2", 0.02, 0.25, 0.98, 0.5);
  padH2->Draw();
  padH2->cd();
  TH1D* PhitsDiff = (TH1D*) pileupHits_shadow->Clone("PhitsDiff");
  PhitsDiff->Add(truePileupHits, -1);
  PhitsDiff->SetTitle("Residual (Shadow - Truth)");
  PhitsDiff->Draw("hist");

  canvas_Phits->cd();
  TPad* padH3 = new TPad("padH3", "padH3", 0.02, 0.02, 0.98, 0.25);
  padH3->Draw();
  padH3->cd();
  TH1D* PhitsSigDiff = (TH1D*) PhitsDiff->Clone("PhitsSigDiff");
  for (int bin = 1; bin <= PhitsSigDiff->GetNbinsX(); ++bin) if(PhitsSigDiff->GetBinError(bin) != 0) PhitsSigDiff->SetBinContent(bin, PhitsSigDiff->GetBinContent(bin)/PhitsSigDiff->GetBinError(bin));
  PhitsSigDiff->SetTitle("Residual (Shadow - Truth) / Bin Error");
  PhitsSigDiff->Draw("hist");

  canvas_Phits->Write();
  //  delete canvas_Phits;

  /////////////////////////////////////////////////////////////////////////////////////
  }

  THStack* energyHistsStack = new THStack("energyHistsStack","MC Pileup Energies");
  energyHistsStack->Add(truePileupEnergies);
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
  eHS_legend->AddEntry(truePileupEnergies,"True pileup","l");
  eHS_legend->AddEntry(pileupEnergies_shadow,"Shadow pileup correction","l");
  eHS_legend->SetBorderSize(1);
  eHS_legend->Draw();

  canvas_Penergies->cd();
  TPad* padE2 = new TPad("padE2", "padE2", 0.02, 0.25, 0.98, 0.5);
  padE2->Draw();
  padE2->cd();
  TH1D* PenergiesDiff = (TH1D*) pileupEnergies_shadow->Clone("PenergiesDiff");
  PenergiesDiff->Add(truePileupEnergies, -1);
  PenergiesDiff->SetTitle("Residual (Shadow - Truth)");
  PenergiesDiff->Draw("hist");

  canvas_Penergies->cd();
  TPad* padE3 = new TPad("padE3", "padE3", 0.02, 0.02, 0.98, 0.25);
  padE3->Draw();
  padE3->cd();
  TH1D* PenergiesSigDiff = (TH1D*) PenergiesDiff->Clone("PenergiesSigDiff");
  for (int bin = 1; bin <= PenergiesSigDiff->GetNbinsX(); ++bin) if(PenergiesSigDiff->GetBinError(bin) != 0) PenergiesSigDiff->SetBinContent(bin, PenergiesSigDiff->GetBinContent(bin)/PenergiesSigDiff->GetBinError(bin));
  PenergiesSigDiff->SetTitle("Residual (Shadow - Truth) / Bin Error");
  PenergiesSigDiff->Draw("hist");

  canvas_Penergies->Write();
  canvas_Penergies->SaveAs("EnergyCanvas.C");
  //  delete canvas_Penergies;
  return -1;
  /////////////////////////////////////////////////////////////////////////////////////

  TCanvas* c1 = new TCanvas("c1","c1");

  TF1* expFunc = new TF1("expFunc", "[0]*exp(-x/[1])", 100, dT_maxtime);
  //  TF1* expFunc = new TF1("expFunc", "[0]*exp(-x/[1])", deadTime*2, dT_maxtime);

  expFunc->SetParameter(1, histMaxTime/nPerFillPerCalo);
  cout << "Lifetime for deltaT hits is: " << histMaxTime/nPerFillPerCalo << endl;

  expFunc->SetParameter(0, 1);
  double par0Start = deltaT_hits->GetEntries() / ( expFunc->Integral(0, dT_maxtime) / deltaT_hits->GetBinWidth(1) );
  expFunc->SetParameter(0, par0Start);

  deltaT_hits->Fit(expFunc, "QRLI");
  deltaT_hits->GetFunction("expFunc")->SetLineColor(2);
  
  cout << "Fit values: All = " << expFunc->GetParameter(1) << " +/- " <<  expFunc->GetParError(1) << endl;

  expFunc->SetParameter(0, 1);
  par0Start = deltaT_true->GetEntries() / ( expFunc->Integral(0, dT_maxtime) / deltaT_true->GetBinWidth(1) );
  expFunc->SetParameter(0, par0Start);

  deltaT_true->Fit(expFunc, "QRLI");
  deltaT_true->GetFunction("expFunc")->SetLineColor(2);
  
  cout << "Fit values: True = " << expFunc->GetParameter(1) << " +/- " <<  expFunc->GetParError(1) << endl;

  outputFile->Write();
  //  delete outputFile;

  /////////////////////////////////////////////////////////////////////////////////////

  return 1;

}
