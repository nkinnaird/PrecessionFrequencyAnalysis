/////////////////////////////////////////////////////////////////////////////////////
// Macro to create some quick time and energy plots from root trees made using KinnairdTree in order to check runs.
// To run with an input file use: root -l runChecker.C+\(\"filepath\"\)
// The + compiles the program, and the \ escapes so bash knows to read the () and "" as actual characters (or something like that).
/////////////////////////////////////////////////////////////////////////////////////

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
// #include <TRatioPlot.h> // wait till root v 6_08

#include <time.h>

// #include "../../ratioMacroHeaders/ratioAnalysisDefs.hh"

using namespace std;

// struct for copying data from input tree - needs to be the same as what was used to generate the trees
  struct CaloData {
    unsigned int run;
    unsigned int event;
    int fillNum;   
    int caloNum;   
    int islandNum; 
    double time;   
    double energy; 
    double x;      
    double y;     
  };



int runChecker(std::string inString)
{
  // create output file that will hold plots
  TFile* outputFile = new TFile("runChecker.root","RECREATE");

/////////////////////////////////////////////////////////////////////////////////////

  // code to be able to run over single or multiple files

  std::vector<string> fileVector;

  std::size_t found = inString.find_last_of(".");
  std::string extension = inString.substr(found);

  if(extension == ".root") fileVector.push_back(inString);
  else if(extension == ".txt"){
    std::ifstream inList(inString);
    string path;
    while(inList >> path) fileVector.push_back(path);
  }
  else{
    printf("Not passing files in correctly. Should be either a single .root file or a .txt file with multiple .root files.\n");
    return 0;
  }

/////////////////////////////////////////////////////////////////////////////////////

      TH1F* allTimes = new TH1F( "Times", "; time (ns); Events", 7000, 0, 700000);
      TH1F* allEnergies = new TH1F( "Energies", "; Energy (~MeV); Events", 400, 0, 4000);

      TH1F* allTimes_threshold = new TH1F( "Times_Threshold", "; time (ns); Events", 7000, 0, 700000);
      TH1F* allEnergies_threshold = new TH1F( "Energies_Threshold", "; Energy (~MeV); Events", 400, 0, 4000);


      for (int fileNum = 0; fileNum < int(fileVector.size()); ++fileNum)
      {
        string line = fileVector.at(fileNum);
        cout << line << '\n';

              // pull in input file
              TFile *inputFile = TFile::Open(line.c_str());
               if (inputFile == 0) {
                  printf("Error: cannot open file\n");
                  continue;
               }
              // grab tree from input file, and link branch address to local variable
              TTree *inputTree = (TTree*)inputFile->Get("kinnairdTree/caloHits");
              if(inputTree == 0){
                  printf("No tree\n");
                  continue;                
              }

              CaloData caloData;
              inputTree->SetBranchAddress("caloHitsBranch", &caloData.run);
              Long64_t treeEntries = inputTree->GetEntries();

/////////////////////////////////////////////////////////////////////////////////////

              for(Long64_t i=0; i<treeEntries; i++) {
                 inputTree->GetEntry(i);
                 int run = caloData.run;
                 string runString = to_string(run);
                 double time = caloData.time * 1.25;
                 double energy = caloData.energy;

                 if(outputFile->FindKey(to_string(run).c_str())) // if the run directory exists already then fill hists
                 {
                   TH1F* runTimes = (TH1F*) outputFile->Get((runString+"/Times"+runString).c_str());
                   TH1F* runEnergies = (TH1F*) outputFile->Get((runString+"/Energies"+runString).c_str());
 
                   runTimes->Fill(time);
                   runEnergies->Fill(energy);
                 }
                 else // else make the hists and start filling
                 {
                   auto runDir = outputFile->mkdir(runString.c_str());
                   runDir->cd();
 
                   TH1F* runTimes = new TH1F( ("Times"+runString).c_str(), "; time (ns); Events", 7000, 0, 700000);
                   TH1F* runEnergies = new TH1F( ("Energies"+runString).c_str(), "; Energy (~MeV); Events", 400, 0, 4000);
 
                   runTimes->Fill(time);
                   runEnergies->Fill(energy);
                 }

                 allTimes->Fill(time);
                 allEnergies->Fill(energy);

                 if(energy > 1800)
                 {
                   allTimes_threshold->Fill(time);
                   allEnergies_threshold->Fill(energy);
                 }

              } // end loop over entries

/////////////////////////////////////////////////////////////////////////////////////

          } // end for loop over files

    outputFile->Write();

    delete outputFile;

  return 1;
}
