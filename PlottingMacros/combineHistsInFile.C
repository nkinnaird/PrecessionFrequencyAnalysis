#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TH1.h>
#include <TFile.h>


int combineHistsInFile()
{

  // TFile* inputFile_defaultConfig = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/60h/StatisticalComparisonFiles/DefaultConfig-NoTimeRand/addedHists-60h-defaultConfig-noTimeRand.root");
  // TFile* inputFile_compareConfig = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/60h/StatisticalComparisonFiles/ComparisonConfig-NoTimeRand/addedHists-60h-CompareConfig-noTimeRand.root");

  // with 5 ns ADT
  TFile* inputFile_defaultConfig = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/60h/StatisticalComparisonFiles/5ns-ADT/addedHists-60h-DefaultConfig-5ns-ADT.root");
  TFile* inputFile_compareConfig = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/60h/StatisticalComparisonFiles/5ns-ADT/addedHists-60h-ComparisonConfig-5ns-ADT.root");



  TFile* outputFile = new TFile("Kinnaird_histsForJosh.root","RECREATE"); // create output file that will hold plots

  TH1F* defaultHist = (TH1F*) inputFile_defaultConfig->Get("topDir/Iter0/Added/Times_E_Threshold");
  TH1F* compareHist = (TH1F*) inputFile_compareConfig->Get("topDir/Iter0/Added/Times_E_Threshold");

  defaultHist->Write("defaultHist");
  compareHist->Write("compareHist");

  outputFile->Write();
  delete outputFile;

  return 1;
}
