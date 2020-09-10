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
#include <TFitResult.h>
#include <TSpectrum.h>
#include <TText.h>
#include <TPaveStats.h>

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;

double calculateAverage(vector<double> inVec)
{
  double avg = 0;
  for (uint i = 0; i < inVec.size(); ++i) avg += inVec.at(i);
  avg /= inVec.size();
  return avg;
}


int DiffRandCorrelation(std::string fileList)
{
  gStyle->SetOptStat(111111);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleStyle(0);
  gStyle->SetPadRightMargin(.05);
  gStyle->SetPadLeftMargin(.15);

  std::vector<string> fileVector;

  std::ifstream inList(fileList);
  string path;
  while(inList >> path) fileVector.push_back(path);

/////////////////////////////////////////////////////////////////////////////////////

  uint totalIters = 100;
  vector<double> randomizations {0, 150, 300, 450, 600, 750};

  // first array is the file number, second array is the randomization number, third array is the iteration number

  double fiveParRs[fileVector.size()][randomizations.size()][totalIters];
  double fiveParRerrs[fileVector.size()][randomizations.size()][totalIters];

  double ratioRs[fileVector.size()][randomizations.size()][totalIters];
  double ratioRerrs[fileVector.size()][randomizations.size()][totalIters];

/////////////////////////////////////////////////////////////////////////////////////


  for (uint fileNum = 0; fileNum < fileVector.size(); ++fileNum)
  {
    TFile* inputFile = TFile::Open(fileVector.at(fileNum).c_str());
     if (inputFile == 0) {
        printf("Error: cannot open file\n");
        return 0;
     }

    cout << "Filenum: " << fileNum << " " << fileVector.at(fileNum) << endl;

/////////////////////////////////////////////////////////////////////////////////////

    for (uint randNum = 0; randNum < randomizations.size(); ++randNum)
    // for (uint randNum = 0; randNum < 1; ++randNum)
    {
      string nameAddOn = "";
      if(randNum != 0) nameAddOn = "_" + to_string(randNum);

      TNtuple *fiveParamFitValues = (TNtuple*)inputFile->Get(("topDir/ToyMC" + nameAddOn + "/Hist/fiveParameterFitValues" + nameAddOn).c_str());
      
      float N_fivePar, N_err_fivePar, tau_fivePar, tau_err_fivePar, A_fivePar, A_err_fivePar, R_fivePar, R_err_fivePar, phi_fivePar, phi_err_fivePar;
      fiveParamFitValues->SetBranchAddress("N", &N_fivePar);
      fiveParamFitValues->SetBranchAddress("N_err", &N_err_fivePar);
      fiveParamFitValues->SetBranchAddress("tau", &tau_fivePar);
      fiveParamFitValues->SetBranchAddress("tau_err", &tau_err_fivePar);
      fiveParamFitValues->SetBranchAddress("A", &A_fivePar);
      fiveParamFitValues->SetBranchAddress("A_err", &A_err_fivePar);
      fiveParamFitValues->SetBranchAddress("R", &R_fivePar);
      fiveParamFitValues->SetBranchAddress("R_err", &R_err_fivePar);
      fiveParamFitValues->SetBranchAddress("phi", &phi_fivePar);
      fiveParamFitValues->SetBranchAddress("phi_err", &phi_err_fivePar);

      TNtuple *ratioFitValues = (TNtuple*)inputFile->Get(("topDir/ToyMC" + nameAddOn + "/Ratio/ratioFitValues" + nameAddOn).c_str());
      
      float A_ratio, A_err_ratio, R_ratio, R_err_ratio, phi_ratio, phi_err_ratio;
      ratioFitValues->SetBranchAddress("A", &A_ratio);
      ratioFitValues->SetBranchAddress("A_err", &A_err_ratio);
      ratioFitValues->SetBranchAddress("R", &R_ratio);
      ratioFitValues->SetBranchAddress("R_err", &R_err_ratio);
      ratioFitValues->SetBranchAddress("phi", &phi_ratio);
      ratioFitValues->SetBranchAddress("phi_err", &phi_err_ratio);

      for (uint iterNum = 0; iterNum < totalIters; ++iterNum)
      {
        fiveParamFitValues->GetEntry(iterNum);
        ratioFitValues->GetEntry(iterNum);

        fiveParRs[fileNum][randNum][iterNum] = R_fivePar;
        fiveParRerrs[fileNum][randNum][iterNum] = R_err_fivePar;

        ratioRs[fileNum][randNum][iterNum] = R_ratio;
        ratioRerrs[fileNum][randNum][iterNum] = R_err_ratio;
      } // end iter loop


      // for (int i = 0; i < totalIters; ++i)
      // {
      //   cout << "File num: " << fileNum << " randNum: " << randNum << " iterNum: " << i << " 5 R: " << fiveParRs[fileNum][randNum][i] << " Ratio R: " << ratioRs[fileNum][randNum][i] << endl;
      // }


    } // end rand num loop

    inputFile->Close();
  } // end file num loop


  cout << "Five parameter first error: " << fiveParRerrs[0][0][0] << endl;
  cout << "Ratio first error: " << ratioRerrs[0][0][0] << endl;


/////////////////////////////////////////////////////////////////////////////////////


  TH1F* corrHists[randomizations.size()];

  TGraphErrors* correlationGraphs[randomizations.size()][totalIters];
  TGraphErrors* correlationAvgGraphs[randomizations.size()];

  TGraphErrors* correlationGraphsOnlyRAvg[randomizations.size()][totalIters];


  int checknum = 0;

  for (uint randNum = 0; randNum < randomizations.size(); ++randNum)
  {
    corrHists[randNum] = new TH1F(("corrHist" + to_string(randNum)).c_str(), "corrHist; Correlation; Events", 100, 0.98, 1);

      correlationAvgGraphs[randNum] = new TGraphErrors();
      correlationAvgGraphs[randNum]->GetXaxis()->SetTitle("Mean Five Parameter R Value (ppm)");
      correlationAvgGraphs[randNum]->GetYaxis()->SetTitle("Mean Ratio R Value (ppm)");
      correlationAvgGraphs[randNum]->SetTitle("Average of Seeds");


    for (uint iterNum = 0; iterNum < totalIters; ++iterNum)
    {
      correlationGraphs[randNum][iterNum] = new TGraphErrors();
      correlationGraphs[randNum][iterNum]->GetXaxis()->SetTitle("Five Parameter R Value (ppm)");
      correlationGraphs[randNum][iterNum]->GetYaxis()->SetTitle("Ratio R Value (ppm)");
      correlationGraphs[randNum][iterNum]->SetTitle(Form("Single Seed %i",iterNum));

      correlationGraphsOnlyRAvg[randNum][iterNum] = new TGraphErrors();
      correlationGraphsOnlyRAvg[randNum][iterNum]->GetXaxis()->SetTitle("Five Parameter R Value (ppm)");
      correlationGraphsOnlyRAvg[randNum][iterNum]->GetYaxis()->SetTitle("Mean Ratio R Value (ppm)");
      correlationGraphsOnlyRAvg[randNum][iterNum]->SetTitle(Form("Seed %i",iterNum));


      for (uint fileNum = 0; fileNum < fileVector.size(); ++fileNum)
      {
        double fiveParR_perfile_perRand_perIter = fiveParRs[fileNum][randNum][iterNum];
        double ratioR_perfile_perRand_perIter = ratioRs[fileNum][randNum][iterNum];

        correlationGraphs[randNum][iterNum]->SetPoint(fileNum, fiveParR_perfile_perRand_perIter, ratioR_perfile_perRand_perIter);
        correlationGraphs[randNum][iterNum]->SetPointError(fileNum, fiveParRerrs[fileNum][randNum][iterNum], ratioRerrs[fileNum][randNum][iterNum]);

        correlationGraphsOnlyRAvg[randNum][iterNum]->SetPoint(fileNum, fiveParR_perfile_perRand_perIter, 0);
      } // end loop over file num

      double actualCorrelationFactor = correlationGraphs[randNum][iterNum]->GetCorrelationFactor();

      // cout << "checknum: " << checknum << " corr: " << actualCorrelationFactor << endl;

      corrHists[randNum]->Fill(actualCorrelationFactor);

      checknum++;
    } // end loop over iter num

/////////////////////////////////////////////////////////////////////////////////////

    for (uint fileNum = 0; fileNum < fileVector.size(); ++fileNum)
    {
      double averageFiveParamR_overSeeds = 0;
      double averageRatioParamR_overSeeds = 0;

      for (uint iterNum = 0; iterNum < totalIters; ++iterNum)
      {
        double fiveR, ratioR;
        correlationGraphs[randNum][iterNum]->GetPoint(fileNum, fiveR, ratioR);

        averageFiveParamR_overSeeds += fiveR;
        averageRatioParamR_overSeeds += ratioR;
      }
      
      averageFiveParamR_overSeeds /= totalIters;
      averageRatioParamR_overSeeds /= totalIters;

      correlationAvgGraphs[randNum]->SetPoint(fileNum, averageFiveParamR_overSeeds, averageRatioParamR_overSeeds);
   

      for (uint iterNum = 0; iterNum < totalIters; ++iterNum)
      {
        double fiveR, ratioR;
        correlationGraphsOnlyRAvg[randNum][iterNum]->GetPoint(fileNum, fiveR, ratioR);
        correlationGraphsOnlyRAvg[randNum][iterNum]->SetPoint(fileNum, fiveR, averageRatioParamR_overSeeds);
      }
    } // end loop over filenums/points

      // for (uint iterNum = 0; iterNum < totalIters; ++iterNum) cout << "corr: " << correlationGraphsOnlyRAvg[randNum][iterNum]->GetCorrelationFactor() << endl; // have to do this in a separate loop after I've run through all the points

    cout << endl;

  } // end loop over rand num


    // single seed correlation vs average correlation

    int plotIterNum = 22;
    int plotRandNum = 1; // 1 for 150 ns smearing

    correlationGraphs[plotRandNum][plotIterNum]->SetTitle("Single Seed Vs Average");
    correlationGraphs[plotRandNum][plotIterNum]->GetXaxis()->SetRangeUser(-4,4);
    correlationGraphs[plotRandNum][plotIterNum]->GetYaxis()->SetRangeUser(-4,4);

    auto corrCanv_both = new TCanvas("corrCanv_both", "corrCanv_both", 200, 100, 600, 600);
    correlationGraphs[plotRandNum][plotIterNum]->DrawClone("APX");
    correlationAvgGraphs[plotRandNum]->SetMarkerColor(2);
    correlationAvgGraphs[plotRandNum]->Draw("PSAME");

    TPaveText *textBoth = new TPaveText(0.22,0.68,0.62,0.83,"NDC");
    textBoth->SetTextSize(0.033);
    textBoth->AddText(Form("Single Seed Correlation = %1.3f",correlationGraphs[plotRandNum][plotIterNum]->GetCorrelationFactor()));
    textBoth->AddText(Form("Average Correlation = %1.3f",correlationAvgGraphs[plotRandNum]->GetCorrelationFactor()));
    ((TText*)textBoth->GetListOfLines()->Last())->SetTextColor(2);
    textBoth->Draw("SAME");

    corrCanv_both->SaveAs("Images/TR_Correlation_Single_Average.png");


    // correlation from average R method vs single seed of T method

    auto corrCanv_avgRmethod = new TCanvas("corrCanv_avgRmethod", "corrCanv_avgRmethod", 600, 100, 600, 600);
    correlationGraphsOnlyRAvg[plotRandNum][plotIterNum]->Draw("APX");
    correlationGraphsOnlyRAvg[plotRandNum][plotIterNum]->GetXaxis()->SetRangeUser(-4,4);
    correlationGraphsOnlyRAvg[plotRandNum][plotIterNum]->GetYaxis()->SetRangeUser(-4,4);

    TPaveText *text = new TPaveText(0.2,0.65,0.6,0.8,"NDC");
    text->SetTextSize(0.04);
    text->AddText(Form("Correlation = %1.3f",correlationGraphsOnlyRAvg[plotRandNum][plotIterNum]->GetCorrelationFactor()));
    text->Draw("SAME");
    
    corrCanv_avgRmethod->SaveAs("Images/TR_Correlation_AvgRMethod_SingleT.png");


    // average correlations vs mean of single seed correlations

    for (uint i = 0; i < randomizations.size(); ++i)
    {
      cout << "randomization level: " << randomizations.at(i) << " average correlation: " << correlationAvgGraphs[i]->GetCorrelationFactor() << " vs mean of random seeds correlation:  " << corrHists[i]->GetMean() << endl;

      if(i == 1){
        auto corrCanv_hist = new TCanvas(Form("corrCanv_hist_%i",i), Form("corrCanv_hist_%i",i), 200+100*i, 200, 600, 600);
        corrHists[i]->SetTitle("Correlation per Seed");
        corrHists[i]->SetName(Form("%3.0f ns rand",randomizations.at(i)));
        gStyle->SetStatBorderSize(1);
        gStyle->SetStatW(0.37);
        gStyle->SetStatX(0.94);
        corrHists[i]->GetXaxis()->SetNdivisions(405, false);
        corrHists[i]->Draw();

        corrCanv_hist->SaveAs("Images/CorrHist_PerRandSeed.png");
      }

    }

    cout << endl;

    for (uint i = 0; i < randomizations.size(); ++i)
    {
      cout << "randomization level: " << randomizations.at(i) << " correlation with only R method averaged: " << correlationGraphsOnlyRAvg[i][plotIterNum]->GetCorrelationFactor() << endl;
    }


/////////////////////////////////////////////////////////////////////////////////////

return 1;

}
