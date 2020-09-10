#include <iostream>
#include <iomanip>
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
#include <TVectorD.h>
#include <TFitResult.h>
#include <TSpectrum.h>
#include <TText.h>
#include <TPaveStats.h>

// #include "ratioAnalysisDefs.hh"
// #include "plotUtils.hh"

using namespace std;

bool compareReconCorrelations = true;


int CorrelationDifferences(std::string firstFileString, std::string secondFileString)
{
  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

/////////////////////////////////////////////////////////////////////////////////////

  // style setting for nice 2D matrix plots

  gStyle->SetPaintTextFormat("1.4f");
  gStyle->SetGridStyle(0);
  gStyle->SetOptStat(kFALSE);

/////////////////////////////////////////////////////////////////////////////////////

  TFile* firstFile = TFile::Open(firstFileString.c_str());
  TFile* secondFile = TFile::Open(secondFileString.c_str());
   if (firstFile == 0 || secondFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  string fileName;
  if(compareReconCorrelations) fileName = "correlationReconDifferences.root";
  else fileName = "correlationDifferences.root";

  TFile* outputFile = new TFile(fileName.c_str(),"RECREATE");

  auto histogramsDir = outputFile->mkdir("Histograms");
  auto canvasesDir = outputFile->mkdir("Canvases");

/////////////////////////////////////////////////////////////////////////////////////

  string pathToHists;
  if(compareReconCorrelations) pathToHists = "MethodAverage/Histograms";
  else pathToHists = "Histograms";


  TDirectory* firstHistogramsDir = firstFile->GetDirectory(pathToHists.c_str());
  TDirectory* secondHistogramsDir = secondFile->GetDirectory(pathToHists.c_str());

  TIter nextKeyFirst(firstHistogramsDir->GetListOfKeys());
  TIter nextKeySecond(secondHistogramsDir->GetListOfKeys());
  TKey* keyFirst;
  TKey* keySecond;


  while((keyFirst = (TKey*)nextKeyFirst()))
  {
    keySecond = (TKey*)nextKeySecond();

    string histName = keyFirst->GetName();
    if(histName.compare("CorrelationMatrix_R_R") != 0) continue; // only make plots and do print out for R v R correlations

    cout << "Histogram name: " << histName << endl;

    TH2D* firstCorr = (TH2D*) keyFirst->ReadObj();
    TH2D* secondCorr = (TH2D*) keySecond->ReadObj();

    TH2D* corrDiff = (TH2D*) firstCorr->Clone();
    corrDiff->Add(secondCorr,-1);


    histogramsDir->cd();

    corrDiff->SetName(("Diff_" + histName).c_str());
    corrDiff->GetZaxis()->SetTitle("Corr Diff");
    corrDiff->LabelsDeflate();
    corrDiff->GetXaxis()->SetLabelSize(0.04);
    corrDiff->GetXaxis()->SetLabelOffset(0.01);
    corrDiff->GetZaxis()->SetTitleOffset(1.7);
    corrDiff->GetYaxis()->SetLabelSize(0.04);
    corrDiff->SetTitle("");

    corrDiff->Write();


        TCanvas* c_corrMatDiff = new TCanvas("c_corrMatDiff","Correlation Matrix Diff",50,10,1100,900);

        gStyle->SetPaintTextFormat("1.4f");

        c_corrMatDiff->SetGrid();
        c_corrMatDiff->SetTopMargin(0.1);
        c_corrMatDiff->SetBottomMargin(0.12);
        c_corrMatDiff->SetRightMargin(0.18);
        c_corrMatDiff->SetLeftMargin(0.18);

        corrDiff->Draw("COLZTEXT");

        if(compareReconCorrelations) c_corrMatDiff->SetName(("Diff_Plot_Recon_" + histName).c_str());
        else c_corrMatDiff->SetName(("Diff_Plot_" + histName).c_str());

        c_corrMatDiff->SaveAs((string(c_corrMatDiff->GetName()) + ".png").c_str());

        canvasesDir->cd();
        c_corrMatDiff->Write();

/////////////////////////////////////////////////////////////////////////////////////

        // do print out of first input correlation

        cout.precision(4);

        cout << endl << endl << "First correlation matrix output: " << endl << endl;

        for (int i = 1; i <= firstCorr->GetNbinsX(); ++i)
        {
          for (int j = 1; j <= firstCorr->GetNbinsY(); ++j)
          {
            cout << setiosflags(ios::fixed) << "& " << firstCorr->GetBinContent(i, j) << " ";
          }
          cout << " \\\\ " <<endl;
        }


        cout << endl << endl << "Difference correlation matrix output: " << endl << endl;

        for (int i = 1; i <= corrDiff->GetNbinsX(); ++i)
        {
          for (int j = 1; j <= corrDiff->GetNbinsY(); ++j)
          {
            cout << setiosflags(ios::fixed) << "& " << corrDiff->GetBinContent(i, j) << " ";
          }
          cout << " \\\\ " <<endl;
        }
  }


/////////////////////////////////////////////////////////////////////////////////////

return 1;

}
