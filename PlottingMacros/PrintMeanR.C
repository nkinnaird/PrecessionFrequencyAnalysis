// 3-31-20: Macro to print out R values from fits to many different random seeds. The input file should be that produced by generalIterPlotter.C.

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TF1.h>
#include <TH1.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TROOT.h>
#include <sstream>

using namespace std;

bool justPrintR = false;


int PrintMeanR(std::string filePath)
{

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  auto Tmethod_Rvalues = ((TH1F*) inputFile->Get("topDir/Added/TMethod/Hists/TMethod_R_Vs_Iter_hist"));
  auto Rmethod_Rvalues = ((TH1F*) inputFile->Get("topDir/Added/FullRatio/Hists/FullRatio_R_Vs_Iter_hist"));

  int numEntries = Tmethod_Rvalues->GetEntries();
  cout << "Number of seeds: " << numEntries << endl;

  cout << "Mean R value of T method random seeds (ppm): " << Tmethod_Rvalues->GetMean() << " width: " << Tmethod_Rvalues->GetRMS() << " error on mean: " << Tmethod_Rvalues->GetRMS()/sqrt(numEntries) << endl;
  cout << "Mean R value of R method random seeds (ppm): " << Rmethod_Rvalues->GetMean() << " width: " << Rmethod_Rvalues->GetRMS() << " error on mean: " << Rmethod_Rvalues->GetRMS()/sqrt(numEntries) << endl;
  cout << "Difference in means (R - T) (ppb): " << (Rmethod_Rvalues->GetMean() - Tmethod_Rvalues->GetMean())*1000. << endl;


	cout << endl;
	return 1;

}
