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

using namespace std;

int numSamples = 1000; // make sure this is set properly, 1000 or 851 for the samples I've done for the Run 1 datasets

double r_To_z(double r){
  return atanh(r);
}

double z_To_r(double z){
  return tanh(z);
}

double calcStatisticalError(double r){

  double var = 1/(double(numSamples) - 3);

  double plusError = z_To_r( r_To_z(r) + sqrt(var) ) - r;
  double minusError = r - z_To_r( r_To_z(r) - sqrt(var) );

  return (plusError + minusError)/2.;
}


void doAveraging(TFile* firstFile, TFile* secondFile, bool compareReconCorrelations)
{

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

    TH2D* corrAvg = (TH2D*) firstCorr->Clone();
    corrAvg->Add(secondCorr,1);

    corrAvg->Scale(0.5);

    if(compareReconCorrelations) corrAvg->SetName(("Avg_Recon_" + histName).c_str());
    else corrAvg->SetName(("Avg_" + histName).c_str());

    corrAvg->GetZaxis()->SetTitle("Corr Avg");
    corrAvg->LabelsDeflate();
    corrAvg->GetXaxis()->SetLabelSize(0.04);
    corrAvg->GetXaxis()->SetLabelOffset(0.01);
    corrAvg->GetZaxis()->SetTitleOffset(1.7);
    corrAvg->GetYaxis()->SetLabelSize(0.04);
    corrAvg->SetTitle("");

    corrAvg->Write();

      // make plot of correlations

        TCanvas* c_corrMatAvg = new TCanvas("c_corrMatAvg","Correlation Matrix Avg",50,10,1100,900);

        gStyle->SetPaintTextFormat("1.4f");

        c_corrMatAvg->SetGrid();
        c_corrMatAvg->SetTopMargin(0.1);
        c_corrMatAvg->SetBottomMargin(0.12);
        c_corrMatAvg->SetRightMargin(0.18);
        c_corrMatAvg->SetLeftMargin(0.18);

        corrAvg->Draw("COLZTEXT");

        if(compareReconCorrelations) c_corrMatAvg->SetName(("Avg_Recon_" + histName).c_str());
        else c_corrMatAvg->SetName(("Avg_" + histName).c_str());

        c_corrMatAvg->SaveAs((string(c_corrMatAvg->GetName()) + ".png").c_str());

/////////////////////////////////////////////////////////////////////////////////////

    // make histogram of the correlation error
    // calculated as the quadrature sum of the statistical and systematic errors, where the systematic error is calculated as the difference between the average and one of the inputs, and the statistical is calculated via the Fischer-Z transform

    TH2D* corrErrMat = (TH2D*) firstCorr->Clone();

    for (int i = 1; i <= corrErrMat->GetNbinsX(); ++i)
    {
      for (int j = 1; j <= corrErrMat->GetNbinsY(); ++j)
      {
        double corr = corrAvg->GetBinContent(i, j);
        double corrErr_syst = abs(corrAvg->GetBinContent(i, j) - firstCorr->GetBinContent(i, j));
        double corrErr_stat = calcStatisticalError(corr);
        double corrErr = sqrt(corrErr_stat*corrErr_stat + corrErr_syst*corrErr_syst);

        if(i==j) corrErr = 0;
        else if(abs(corr-1) < 1e-4) corrErr = 0;
        else if(i!=j && corrErr < 1e-4) corrErr = 1e-4;

        corrErrMat->SetBinContent(i, j, corrErr);

        // cout << "i: " << i << " j: " << j << " corr err: " << corrErr << " corr: " << corr << endl;
      }
    }

    if(compareReconCorrelations) corrErrMat->SetName(("Avg_Recon_" + histName + "_Err").c_str());
    else corrErrMat->SetName(("Avg_" + histName + "_Err").c_str());

    corrErrMat->GetZaxis()->SetTitle("Corr Err");
    corrErrMat->GetZaxis()->SetRangeUser(0, 1);
    corrErrMat->LabelsDeflate();
    corrErrMat->GetXaxis()->SetLabelSize(0.04);
    corrErrMat->GetXaxis()->SetLabelOffset(0.01);
    corrErrMat->GetZaxis()->SetTitleOffset(1.7);
    corrErrMat->GetYaxis()->SetLabelSize(0.04);
    corrErrMat->SetTitle("");

    corrErrMat->Write();

/////////////////////////////////////////////////////////////////////////////////////


        // do print out of correlations for latex tables
        cout.precision(4);

        cout << endl << endl << "Average correlation matrix output: " << endl << endl;

        for (int i = 1; i <= corrAvg->GetNbinsX(); ++i)
        {
          for (int j = 1; j <= corrAvg->GetNbinsY(); ++j)
          {
            double corr = corrAvg->GetBinContent(i, j);
            double corrErr = corrErrMat->GetBinContent(i, j);

            // cout << "stat err: " << corrErr_stat << " vs syst err: " << corrErr_syst << " diff: " << corrErr_stat - corrErr_syst << endl;

            cout << setiosflags(ios::fixed) << "& " << corr << " " << corrErr << " ";
          }
          cout << " \\\\ " <<endl;
        }
  }

}



int CorrelationAveragesAndErrors(std::string firstFileString, std::string secondFileString)
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

  TFile* outputFile = new TFile("correlationAverage.root","RECREATE");

  doAveraging(firstFile, secondFile, false);
  doAveraging(firstFile, secondFile, true);

  cout << endl << endl << "As a reminder, number of trials is set to: " << numSamples << ", make sure that is right." << endl;

/////////////////////////////////////////////////////////////////////////////////////

return 1;

}
