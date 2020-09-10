// 3-31-20: I don't remember what this macro does.

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
#include <TVectorD.h>
// #include <TRatioPlot.h> // wait till root v 6_08

#include "../../ratioMacroHeaders/ratioAnalysisDefs.hh"
#include "../../ratioMacroHeaders/copyFile.hh"

int numIters = 3;


void graphDiffs(double* xVals, double* yVals, string ID)
{
  TGraph* myGraph = new TGraph(numIters, xVals, yVals);

  myGraph->GetXaxis()->SetTitle("Iter");
  myGraph->GetYaxis()->SetTitle(("#Delta" + ID).c_str());

  myGraph->SetTitle(("#Delta" + ID).c_str());
  myGraph->SetName(ID.c_str());
  myGraph->Write();
}



int DifferentPileupErrors(std::string filePath)
{
    // create output file that will hold plots
  TFile* outputFile = new TFile("diffErrsCompare.root","RECREATE");

  // make top directory for output file
  auto topDir = outputFile->mkdir("topDir");
  topDir->cd();

  TFile *inputFile  = TFile::Open(filePath.c_str());
  if (inputFile == 0) {
     printf("Error: cannot open file\n");
     return 0;
  }

/////////////////////////////////////////////////////////////////////////////////////

  double pN[numIters];

  TF1* fiveParameterFits[numIters];

  double five_chi2s[numIters];
  double five_chi2ndfs[numIters];
  double five_Pvals[numIters];
  double five_Ns[numIters];
  double five_Taus[numIters];
  double five_As[numIters];
  double five_Rs[numIters];
  double five_sigmaRs[numIters];
  double five_Phis[numIters];

  double diff_five_chi2s[numIters];
  double diff_five_chi2ndfs[numIters];
  double diff_five_Pvals[numIters];
  double diff_five_Ns[numIters];
  double diff_five_Taus[numIters];
  double diff_five_As[numIters];
  double diff_five_Rs[numIters];
  double diff_five_sigmaRs[numIters];
  double diff_five_Phis[numIters];

/////////////////////////////////////////////////////////////////////////////////////

  TF1* ratioFits[numIters];

  double ratio_chi2s[numIters];
  double ratio_chi2ndfs[numIters];
  double ratio_Pvals[numIters];
  double ratio_As[numIters];
  double ratio_Rs[numIters];
  double ratio_sigmaRs[numIters];
  double ratio_Phis[numIters];

  double diff_ratio_chi2s[numIters];
  double diff_ratio_chi2ndfs[numIters];
  double diff_ratio_Pvals[numIters];
  double diff_ratio_As[numIters];
  double diff_ratio_Rs[numIters];
  double diff_ratio_sigmaRs[numIters];
  double diff_ratio_Phis[numIters];

/////////////////////////////////////////////////////////////////////////////////////

  for (int iter = 0; iter < numIters; ++iter)
  {
    pN[iter] = iter;
/*
    TH1F* fiveParamHist = (TH1F*) inputFile->Get(Form("topDir/ToyMCIter%i/SubIter3/Hist/Toy_5_Param_Hist", iter));
    fiveParameterFits[iter] = (TF1*) fiveParamHist->GetFunction("fiveParamFit");  

    five_chi2s[iter] = fiveParameterFits[iter]->GetChisquare();
    five_chi2ndfs[iter] = fiveParameterFits[iter]->GetChisquare()/fiveParameterFits[iter]->GetNDF();
    five_Pvals[iter] = fiveParameterFits[iter]->GetProb();
    five_Ns[iter] = fiveParameterFits[iter]->GetParameter(0);
    five_Taus[iter] = fiveParameterFits[iter]->GetParameter(1);
    five_As[iter] = fiveParameterFits[iter]->GetParameter(2);
    five_Rs[iter] = fiveParameterFits[iter]->GetParameter(3);
    five_sigmaRs[iter] = fiveParameterFits[iter]->GetParError(3);
    five_Phis[iter] = fiveParameterFits[iter]->GetParameter(4);

    diff_five_chi2s[iter] = five_chi2s[iter] - five_chi2s[0];
    diff_five_chi2ndfs[iter] = five_chi2ndfs[iter] - five_chi2ndfs[0];
    diff_five_Pvals[iter] = five_Pvals[iter] - five_Pvals[0];
    diff_five_Ns[iter] = five_Ns[iter] - five_Ns[0];
    diff_five_Taus[iter] = five_Taus[iter] - five_Taus[0];
    diff_five_As[iter] = five_As[iter] - five_As[0];
    diff_five_Rs[iter] = five_Rs[iter] - five_Rs[0];
    diff_five_sigmaRs[iter] = five_sigmaRs[iter] - five_sigmaRs[0];
    diff_five_Phis[iter] = five_Phis[iter] - five_Phis[0];
*/
/////////////////////////////////////////////////////////////////////////////////////

    // TGraph* ratioGraph = (TGraph*) inputFile->Get(Form("topDir/ToyMCIter%i/SubIter3/Ratio/Toy_Ratio_Graph", iter));
    // ratioFits[iter] = (TF1*) ratioGraph->GetFunction("ratioFitIntegral");

    TGraph* ratioGraph = (TGraph*) inputFile->Get(Form("topDir/FitPasses/FitPass%i/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", iter));
    ratioFits[iter] = (TF1*) ratioGraph->GetFunction("fullRatioFitFunc");

    ratio_chi2s[iter] = ratioFits[iter]->GetChisquare();
    ratio_chi2ndfs[iter] = ratioFits[iter]->GetChisquare()/ratioFits[iter]->GetNDF();
    ratio_Pvals[iter] = ratioFits[iter]->GetProb();
    ratio_As[iter] = ratioFits[iter]->GetParameter(0);
    ratio_Rs[iter] = ratioFits[iter]->GetParameter(1);
    ratio_sigmaRs[iter] = ratioFits[iter]->GetParError(1);
    ratio_Phis[iter] = ratioFits[iter]->GetParameter(2);

    diff_ratio_chi2s[iter] = ratio_chi2s[iter] - ratio_chi2s[0];
    diff_ratio_chi2ndfs[iter] = ratio_chi2ndfs[iter] - ratio_chi2ndfs[0];
    diff_ratio_Pvals[iter] = ratio_Pvals[iter] - ratio_Pvals[0];
    diff_ratio_As[iter] = ratio_As[iter] - ratio_As[0];
    diff_ratio_Rs[iter] = ratio_Rs[iter] - ratio_Rs[0];
    diff_ratio_sigmaRs[iter] = ratio_sigmaRs[iter] - ratio_sigmaRs[0];
    diff_ratio_Phis[iter] = ratio_Phis[iter] - ratio_Phis[0];

  }
/*
  auto fiveFitDir = topDir->mkdir("FiveParam");
  fiveFitDir->cd();

    graphDiffs(pN, diff_five_chi2s, "Chi2");
    graphDiffs(pN, diff_five_chi2ndfs, "Chi2NDF");
    graphDiffs(pN, diff_five_Pvals, "Pval");
    graphDiffs(pN, diff_five_Ns, "N");
    graphDiffs(pN, diff_five_Taus, "Tau");
    graphDiffs(pN, diff_five_As, "A");
    graphDiffs(pN, diff_five_Rs, "R");
    graphDiffs(pN, diff_five_sigmaRs, "sigmaR");
    graphDiffs(pN, diff_five_Phis, "Phi");
*/
  auto ratioFitDir = topDir->mkdir("RatioFit");
  ratioFitDir->cd();
  
    graphDiffs(pN, diff_ratio_chi2s, "Chi2");
    graphDiffs(pN, diff_ratio_chi2ndfs, "Chi2NDF");
    graphDiffs(pN, diff_ratio_Pvals, "Pval");
    graphDiffs(pN, diff_ratio_As, "A");
    graphDiffs(pN, diff_ratio_Rs, "R");
    graphDiffs(pN, diff_ratio_sigmaRs, "sigmaR");
    graphDiffs(pN, diff_ratio_Phis, "Phi");

/////////////////////////////////////////////////////////////////////////////////////

  outputFile->Write();
  delete outputFile;


  return 1;
 }
