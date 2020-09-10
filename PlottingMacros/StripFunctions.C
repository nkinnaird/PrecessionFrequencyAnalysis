// 3-31-20: Macro which strips functions out of a fit file from ratioMacro.C, and puts them in a similar directory structure as the original file. 
// This was written for my average fit start scans, where the large number of files with fits and various histograms were taking up too space.

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

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;

int StripFunctions(string fileList)
{
  vector<string> fileVector;
  std::ifstream inList(fileList);
  string path;
  while(inList >> path) fileVector.push_back(path);

  uint totalFiles = fileVector.size();
  uint filesToRunOver = totalFiles;

  // check if dataset tag exists in file and if so append it to the file name and write it to the file

  TFile* inputFile = TFile::Open(fileVector.at(0).c_str());

  string outputFileName = "outputFunctions";
  TNamed* tag = applyDatasetTag(inputFile, outputFileName);

  TFile* outputFile = new TFile((outputFileName + ".root").c_str(),"RECREATE");
  if(tag) tag->Write();

  inputFile->Close();

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


// loop through files 
  for (uint fileNum = 0; fileNum < filesToRunOver; fileNum++)
  {
    cout << "File num: " << fileNum << endl;

    TFile *thisFile = TFile::Open(fileVector.at(fileNum).c_str());
     if (thisFile == 0) {
        cout << "Error: cannot open file number " << fileNum << endl;
        return 0;
     }

    // in the future could do someone where the directories are also pulled from the input files, but for now just hard code them in
    auto fileDir = outputFile->mkdir(Form("File%i", fileNum));
    auto topDir = fileDir->mkdir("topDir");
    auto passesDir = topDir->mkdir("FitPasses");

      // get number of passes to loop over
      TNtuple *firstTuple = (TNtuple*)thisFile->Get(Form("topDir/FitPasses/FitPass0/FitConditions0"));
        if(firstTuple == 0){
          cout << "File num: " << fileNum << " failed somehow - no first fit conditions ntuple - exiting. Bad file: " << fileVector.at(fileNum) << endl;
          exit(1);
        }

      float totalFitPasses;
      firstTuple->SetBranchAddress("totalPasses", &totalFitPasses);
      firstTuple->GetEntry(0);

      for (int fitPass = 0; fitPass < totalFitPasses; ++fitPass)
      {
        auto passDir = passesDir->mkdir(Form("FitPass%i", fitPass));
        passDir->cd();

        TNtuple *ntuple = (TNtuple*)thisFile->Get(Form("topDir/FitPasses/FitPass%d/FitConditions%d", fitPass, fitPass));
          if(ntuple == 0){
            cout << "File num: " << fileNum << " fit pass: " << fitPass << " failed somehow - no fit conditions ntuple - exiting. Bad file: " << fileVector.at(fileNum) << endl;
            exit(1);
          }
        TNtuple *clonedNtuple = (TNtuple*) ntuple->CloneTree(); // just saving the tree doesn't save the data in the tree which is on disk and not in memory, need to clone it first
        clonedNtuple->Write();

        auto addedDir = passDir->mkdir("addedDir");
        auto TmethodDir = addedDir->mkdir("TMethod");
        auto fullRatioDir = addedDir->mkdir("FullRatio");

        TH1F* TmethodHist = (TH1F*) thisFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/TMethod/allTimesAdded_TMethod", fitPass));
        TF1* TMethod_func;
        if(TmethodHist) TMethod_func = (TF1*) TmethodHist->GetFunction("TmethodFitFunc");
        else {
          cout << "File num: " << fileNum << " fit pass: " << fitPass << " T method fit failed - exiting. Bad file: " << fileVector.at(fileNum) << endl;
          exit(1);
        }

        TmethodDir->cd();
        TMethod_func->Write();

        TGraphErrors* RmethodGraph = (TGraphErrors*) thisFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass));
        TF1* RMethod_func;
        if(RmethodGraph) RMethod_func = (TF1*) RmethodGraph->GetFunction("fullRatioFitFunc");
        else {
          cout << "File num: " << fileNum << " fit pass: " << fitPass << " R method fit failed - exiting. Bad file: " << fileVector.at(fileNum) << endl;
          exit(1);
        }

        fullRatioDir->cd();
        RMethod_func->Write();

        // clear objects from memory
        TmethodHist->Delete();
        RmethodGraph->Delete();
        ntuple->Delete();
        clonedNtuple->Delete();
      }

      delete thisFile;
  }


/////////////////////////////////////////////////////////////////////////////////////

      delete outputFile;


  return 1;
}
