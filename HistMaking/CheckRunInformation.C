// 3-31-20: Older macro used to print out information relating to which runs, subruns, and fills are contained within the passed in file.

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
#include "../ratioMacroHeaders/ratioHistogramsConfig.hh" // don't change these paths since this file is run on the grid (whereas other files can see my root include directory and the ../ can be removed)
#include "../ratioMacroHeaders/pileupUtils.hh"


int CheckRunInformation(std::string filePath)
{

/////////////////////////////////////////////////////////////////////////////////////

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 1;
   }

/////////////////////////////////////////////////////////////////

 // need to create unique output root file names for parallel processing - so grab the input file name which has a unique id

  size_t found = filePath.find_first_of("_"); // using underscore here would not work if the directory path has any underscores in it or if the file names change
  string uniqueID = filePath.substr(found+1);
  found = uniqueID.find_first_of("_");
  uniqueID = uniqueID.substr(found+1);
  found = uniqueID.find_last_of(".");
  uniqueID = uniqueID.substr(0, found);

  // TFile* outputFile = new TFile(("run_info_" + uniqueID + ".root").c_str(),"RECREATE");

/////////////////////////////////////////////////////////////////////////////////////

    
    unsigned int event_Get;
    unsigned int subRun_Get;
    unsigned int run_Get;


  TTree* inputTree = (TTree*)inputFile->Get("clusterTree/clusters");
  Long64_t treeEntries = inputTree->GetEntries();


  inputTree->SetBranchAddress("eventNum", &event_Get);
  inputTree->SetBranchAddress("subRunNum", &subRun_Get);
  inputTree->SetBranchAddress("runNum", &run_Get);


/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

  set<int> subRunSet;
  set<double> fillSet;

  ofstream subRunInfo;
  subRunInfo.open("subRunInfo_" + uniqueID + ".txt");


  // loop through tree and fill maps containing the vectors of hits

    for(Long64_t i=0; i<treeEntries; i++) {
       inputTree->GetEntry(i);

       int subRun_ID_int = run_Get * 1e3 + subRun_Get;
       double fill_ID_int = run_Get * 1e6 + subRun_Get*1e3 + event_Get;

       subRunSet.insert(subRun_ID_int);
       fillSet.insert(fill_ID_int);

     }

     cout << "number of fills: " << fillSet.size() << endl << endl;
     cout << "number of run + subruns: " << subRunSet.size() << endl << endl;

     for (auto &i : subRunSet)
     {
       cout << i << endl;
       subRunInfo << i << endl;
     }

  subRunInfo.close();

/////////////////////////////////////////////////////////////////////////////////////

   return 0;

}
