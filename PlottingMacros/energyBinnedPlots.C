// Macro for plots for different energy bins.

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TImage.h>
#include <sstream>
#include <TNtuple.h>
#include <TPaveStats.h>
#include <TVectorD.h>

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;

void makePlots(TDirectory* inDir, std::vector<TF1*> fits, string titleStr, vector<pair<float, float>> energyRanges)
{
  inDir->cd();

  string plotName;
  string xAxisTitle;

  xAxisTitle = "Cluster Energy [MeV]";
  plotName = "EBin";
  string fitTitle = inDir->GetName();

  int numParams = fits.at(0)->GetNpar();

/////////////////////////////////////////////////////////////////////////////////////

  TGraphErrors* Chi2NDF_Vs_Iter = new TGraphErrors();
  Chi2NDF_Vs_Iter->SetName((fitTitle + "_Chi2NDF_Vs_" + plotName).c_str());
  Chi2NDF_Vs_Iter->SetTitle((titleStr + " #chi^{2}/NDF Vs " + xAxisTitle).c_str());
  Chi2NDF_Vs_Iter->GetXaxis()->SetTitle(xAxisTitle.c_str());
  Chi2NDF_Vs_Iter->GetYaxis()->SetTitle("#chi^{2}/NDF");
  Chi2NDF_Vs_Iter->GetYaxis()->SetTitleOffset(2.0);

  TGraph* Pval_Vs_Iter = new TGraph();
  Pval_Vs_Iter->SetName((fitTitle + "_Pval_Vs_" + plotName).c_str());
  Pval_Vs_Iter->SetTitle((titleStr + " P Value Vs " + xAxisTitle).c_str());
  Pval_Vs_Iter->GetXaxis()->SetTitle(xAxisTitle.c_str());
  Pval_Vs_Iter->GetYaxis()->SetTitle("P value");
  Pval_Vs_Iter->GetYaxis()->SetTitleOffset(2.0);
  Pval_Vs_Iter->GetXaxis()->SetRangeUser(0, fits.size()+1);

  for (uint fitNum = 0; fitNum < fits.size(); ++fitNum)
  {
    double xPoint = energyRanges.at(fitNum).first;

    Chi2NDF_Vs_Iter->SetPoint(fitNum, xPoint, fits.at(fitNum)->GetChisquare()/fits.at(fitNum)->GetNDF());
    Chi2NDF_Vs_Iter->SetPointError(fitNum, 0, sqrt(2./fits.at(fitNum)->GetNDF()));

    Pval_Vs_Iter->SetPoint(fitNum, xPoint, fits.at(fitNum)->GetProb());
  }

    Chi2NDF_Vs_Iter->Write();
    Pval_Vs_Iter->Write();


    string canvasName = (fitTitle + "_Chi2NDF_Vs_" + plotName + "_Canv");
    auto chi2_canv = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,500,400);
    Chi2NDF_Vs_Iter->Draw("AP");

    drawHorizontalLine(chi2_canv, 1);

    chi2_canv->SaveAs(("Images/" + canvasName + datasetTagForPlots + ".png").c_str());

    TH1F* Chi2NDF_Vs_Iter_hist = projectGraphToHist(Chi2NDF_Vs_Iter);
    auto chi2_hist_canv = new TCanvas((canvasName + "_hist").c_str(),canvasName.c_str(),200,10,500,400);
    Chi2NDF_Vs_Iter_hist->Draw();

      double mean = Chi2NDF_Vs_Iter->GetMean(2);
      double rms = Chi2NDF_Vs_Iter->GetRMS(2);
      TPaveText* textBox = new TPaveText(0.75,0.7,0.975,0.9, "NDC");
      textBox->AddText(Form("Entries  %zu", fits.size()));
      textBox->AddText(Form("Mean  %0.4g", mean));
      textBox->AddText(Form("RMS  %0.4g", rms));
      textBox->SetBorderSize(1);
      textBox->SetFillStyle(1001);
      textBox->Draw("SAME");

    chi2_hist_canv->SaveAs(("Images/" + canvasName + "_hist" + datasetTagForPlots + ".png").c_str());

    canvasName = fitTitle + "_PVal_Vs_" + plotName + "_Canv";
    auto pval_canv = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,500,400);
    Pval_Vs_Iter->Draw("AP");
    pval_canv->SaveAs(("Images/" + canvasName + datasetTagForPlots + ".png").c_str());

/////////////////////////////////////////////////////////////////////////////////////

  for (int parNum = 0; parNum < numParams; ++parNum)
  {
    string paramString = fits.at(0)->GetParName(parNum);
    string nameString = removeDangerousCharacters(paramString);

    bool phaseParam = false;
    string strAddon = "";
    if(paramString.compare("R") == 0) strAddon = " [ppm]";
    else if(paramString.find("tau") != string::npos) strAddon = " [#mus]";
    else if(paramString.find("omega") != string::npos) strAddon = " [rad/#mus]";
    else if(paramString.find("phi") != string::npos) phaseParam = true;

    TGraphErrors* param_Vs_Iter = new TGraphErrors();
    param_Vs_Iter->SetName((fitTitle + "_" + nameString + "_Vs_" + plotName).c_str());
    param_Vs_Iter->SetTitle((titleStr + " " + paramString + " Vs " + xAxisTitle).c_str());
    param_Vs_Iter->GetYaxis()->SetTitle((paramString + strAddon).c_str());
    param_Vs_Iter->GetYaxis()->SetTitleOffset(2.0);
    param_Vs_Iter->GetXaxis()->SetTitle(xAxisTitle.c_str());
      
      for (uint fitNum = 0; fitNum < fits.size(); ++fitNum)
      {
        double xPoint = energyRanges.at(fitNum).first;

        double parameterValue = fits.at(fitNum)->GetParameter(parNum);
        double parameterError = fits.at(fitNum)->GetParError(parNum);
        
        if(phaseParam) normalizePhase(parameterValue);
        if(paramString.compare(A_string) == 0 && xPoint < 1000 && parameterValue > 0) flipParameter(parameterValue); // flip sign of asymmetry parameter so it goes negative at low energies - should maybe do something with the phase too

        param_Vs_Iter->SetPoint(fitNum, xPoint, parameterValue);
        param_Vs_Iter->SetPointError(fitNum, 0, parameterError);
      }
    
    if(paramString.compare(A_string) == 0){   
      param_Vs_Iter->Fit("pol5", "Q");
      param_Vs_Iter->GetFunction("pol5")->SetLineColor(2); 
    }
    else if(paramString.compare(kLoss_string) == 0){
      param_Vs_Iter->Fit("pol0", "Q", "", 1100, 3100);
      param_Vs_Iter->GetFunction("pol0")->SetLineColor(2);       
    }

    param_Vs_Iter->Write();

      string canvasName = (fitTitle + "_" + nameString + "_Vs_" + plotName + "_Canv");
      auto graphCanvas = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,500,400);
      param_Vs_Iter->Draw("AP");

      graphCanvas->Update();
      TPaveStats *statsBox = (TPaveStats*)graphCanvas->GetPrimitive("stats");
      if(statsBox != 0){
        statsBox->SetBorderSize(1);
        if(paramString.compare(A_string) == 0){
          statsBox->SetX1NDC(0.325);
          statsBox->SetX2NDC(0.625);
        }
        statsBox->Draw("SAME");
      }

    graphCanvas->SaveAs(("Images/" + canvasName + datasetTagForPlots + ".png").c_str());
      
  } // end loop over parameters

}

/////////////////////////////////////////////////////////////////////////////////////

int energyBinnedPlots(std::string filePath)
{
  // gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  // check if dataset tag exists in file and if so append it to the file name and write it to the file

  string outputFileName = "energyBinnedPlots";
  TNamed* tag = applyDatasetTag(inputFile, outputFileName);

  TFile* outputFile = new TFile((outputFileName + ".root").c_str(),"RECREATE");
  if(tag) tag->Write();
  auto topDir = outputFile->mkdir("topDir");

/////////////////////////////////////////////////////////////////////////////////////
  // These only get set for the interactive root session (any generated canvases, etc.), but does not apply to the output root file - that comes from .rootlogon.C
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(1111);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1);
  gStyle->SetLineColor(1);
  gStyle->SetPadRightMargin(.08);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
      // Plots over different fit conditions

  auto addedDir = topDir->mkdir("Added");

vector<pair<float, float>> energyRanges;

vector<TF1*> TmethodFits;
vector<TF1*> fullRatioFits;

      TNtuple *firstTuple = (TNtuple*)inputFile->Get(Form("topDir/FitPasses/FitPass0/FitConditions0"));
      float tP;
      firstTuple->SetBranchAddress("totalPasses", &tP);
      firstTuple->GetEntry(0);

    for (int fitPass = 0; fitPass < tP; ++fitPass)
    {
      TVectorD* histSavedParameters = (TVectorD*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/parameterStore", fitPass));
      float lowerThresh = (*histSavedParameters)[1];
      float upperThresh = (*histSavedParameters)[16];
      energyRanges.emplace_back(lowerThresh, upperThresh);

      TF1* TMethodFitFunction = (TF1*) ((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/TMethod/allTimesAdded_TMethod", fitPass)))->GetFunction("TmethodFitFunc")->Clone();
      TmethodFits.push_back(TMethodFitFunction);
      
      TGraphErrors* ratioGraph = (TGraphErrors*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass));
      TF1* fullRatioFitFunction = 0;
      if(ratioGraph != 0){ // check to make sure ratio fits exist before trying to make plots with them
        fullRatioFitFunction = (TF1*) ratioGraph->GetFunction("fullRatioFitFunc")->Clone();
        fullRatioFits.push_back(fullRatioFitFunction);
      }

    }

  auto TmethodDir = addedDir->mkdir("TMethod");
  makePlots(TmethodDir, TmethodFits, "T Method Fit", energyRanges);

  auto fullRatioDir = addedDir->mkdir("FullRatio");
  if(fullRatioFits.size() > 0) makePlots(fullRatioDir, fullRatioFits, "Full Ratio Fit", energyRanges);

/////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////

  delete outputFile;

  return 1;
}
