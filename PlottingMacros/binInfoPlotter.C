// 3-31-20: Macro for plots for scans over bin widths or bin edges.

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

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;

bool binEdgeShift = true; // true to plot against bin edges, false to plot against bin widths

void makePlots(TDirectory* inDir, vector<TF1*> fits, string titleStr, vector<pair<double, double> > binInfo, bool saveImages)
{
  inDir->cd();

  string plotName;
  string xAxisTitle;
  string addOn = " [ns]";

  if(binEdgeShift){
   plotName = "binEdgeShift";
   xAxisTitle = "Bin Edge Shift";
  }
  else{
   plotName = "binWidth";
   xAxisTitle = "Bin Width";
  }

  string fitTitle = inDir->GetName();

  int numParams = fits.at(0)->GetNpar();

/////////////////////////////////////////////////////////////////////////////////////

  TGraphErrors* Chi2NDF_Vs_Iter = new TGraphErrors();
  Chi2NDF_Vs_Iter->SetName((fitTitle + "_Chi2NDF_Vs_" + plotName).c_str());
  Chi2NDF_Vs_Iter->SetTitle((titleStr + " #chi^{2}/NDF Vs " + xAxisTitle).c_str());
  Chi2NDF_Vs_Iter->GetXaxis()->SetTitle((xAxisTitle + addOn).c_str());
  Chi2NDF_Vs_Iter->GetYaxis()->SetTitle("#chi^{2}/NDF");
  Chi2NDF_Vs_Iter->GetYaxis()->SetTitleOffset(2.0);

  TGraph* Pval_Vs_Iter = new TGraph();
  Pval_Vs_Iter->SetName((fitTitle + "_Pval_Vs_" + plotName).c_str());
  Pval_Vs_Iter->SetTitle((titleStr + " P Value Vs " + xAxisTitle).c_str());
  Pval_Vs_Iter->GetXaxis()->SetTitle((xAxisTitle + addOn).c_str());
  Pval_Vs_Iter->GetYaxis()->SetTitle("P value");
  Pval_Vs_Iter->GetYaxis()->SetTitleOffset(2.0);
  Pval_Vs_Iter->GetXaxis()->SetRangeUser(0, fits.size()+1);

  for (uint fitNum = 0; fitNum < fits.size(); ++fitNum)
  {
    double xPoint = (binEdgeShift) ? binInfo.at(fitNum).second : binInfo.at(fitNum).first;

    Chi2NDF_Vs_Iter->SetPoint(fitNum, xPoint, fits.at(fitNum)->GetChisquare()/fits.at(fitNum)->GetNDF());
    Chi2NDF_Vs_Iter->SetPointError(fitNum, 0, sqrt(2./fits.at(fitNum)->GetNDF()));

    Pval_Vs_Iter->SetPoint(fitNum, xPoint, fits.at(fitNum)->GetProb());
  }


  Chi2NDF_Vs_Iter->Write();
  Pval_Vs_Iter->Write();
/*
  if(saveImages){
    string canvasName = (fitTitle + "_Chi2NDF_Vs_" + plotName + "_Canv");
    auto tempCanv1 = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,500,400);
    Chi2NDF_Vs_Iter->GetXaxis()->SetLimits(xVector.front()-1, xVector.back()+1);
    Chi2NDF_Vs_Iter->Draw("AP");

      double mean = Chi2NDF_Vs_Iter->GetMean(2);
      double rms = Chi2NDF_Vs_Iter->GetRMS(2);
      TPaveText* textBox = new TPaveText(0.75,0.7,0.975,0.9, "NDC");
      textBox->AddText(Form("Mean  %0.4g", mean));
      textBox->AddText(Form("RMS  %0.4g", rms));
      textBox->SetBorderSize(1);
      textBox->SetFillStyle(1001);
      textBox->Draw("SAME");

    tempCanv1->SaveAs(("Images/" + canvasName + datasetTagForPlots + ".png").c_str());

    canvasName = fitTitle + "_PVal_Vs_" + plotName + "_Canv";
    auto tempCanv2 = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,500,400);
    Pval_Vs_Iter->GetXaxis()->SetLimits(xVector.front()-1, xVector.back()+1);
    Pval_Vs_Iter->Draw("AP");
    tempCanv2->SaveAs(("Images/" + canvasName + datasetTagForPlots + ".png").c_str());
  }
*/
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
    param_Vs_Iter->GetXaxis()->SetTitle((xAxisTitle + addOn).c_str());
      
      for (uint fitNum = 0; fitNum < fits.size(); ++fitNum)
      {
        double xPoint = (binEdgeShift) ? binInfo.at(fitNum).second : binInfo.at(fitNum).first;

        double parameterValue = fits.at(fitNum)->GetParameter(parNum);
        double parameterError = fits.at(fitNum)->GetParError(parNum);
        
        if(phaseParam) normalizePhase(parameterValue);

        param_Vs_Iter->SetPoint(fitNum, xPoint, parameterValue);
        param_Vs_Iter->SetPointError(fitNum, 0, parameterError);
      }

    
    if(paramString.compare("R") == 0){
        param_Vs_Iter->Fit("pol1", "Q"); // fit for the systematic effect on R
        param_Vs_Iter->GetFunction("pol1")->SetLineColor(2);
    }

    param_Vs_Iter->Write();

    if(saveImages && paramString.compare("R") == 0){
      string canvasName = (fitTitle + "_" + nameString + "_Vs_" + plotName + "_Canv");
      auto graphCanvas = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,600,400);
      param_Vs_Iter->Draw("AP");

      graphCanvas->Update();
      TPaveStats *statsBox = (TPaveStats*)graphCanvas->GetPrimitive("stats");
      statsBox->SetBorderSize(1);
      statsBox->Draw("SAME");

      graphCanvas->SaveAs(("Images/" + canvasName + datasetTagForPlots + ".png").c_str());
    }

  }

}

/////////////////////////////////////////////////////////////////////////////////////

int binInfoPlotter(std::string filePath)
{
  // gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  // check if dataset tag exists in file and if so append it to the file name and write it to the file

  string outputFileName = "plotsVsBinInfo";
  TNamed* tag = applyDatasetTag(inputFile, outputFileName);

  TFile* outputFile = new TFile((outputFileName + ".root").c_str(),"RECREATE");
  if(tag) tag->Write();
  auto topDir = outputFile->mkdir("topDir");


/////////////////////////////////////////////////////////////////////////////////////
  // These only get set for the interactive root session (any generated canvases, etc.), but does not apply to the output root file - that comes from .rootlogon.C
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(2);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1);
  gStyle->SetLineColor(1);
  gStyle->SetPadRightMargin(.08);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
      // Plots over different fit conditions

  auto addedDir = topDir->mkdir("Added");

std::vector<TF1*> TmethodFits;
std::vector<TF1*> fullRatioFits;

std::vector<pair<double, double> > binInfo;


      TNtuple *firstTuple = (TNtuple*)inputFile->Get(Form("topDir/FitPasses/FitPass0/FitConditions0"));
      float tP;
      firstTuple->SetBranchAddress("totalPasses", &tP);
      firstTuple->GetEntry(0);

    for (int fitPass = 0; fitPass < tP; ++fitPass)
    {
      TH1F* histForInfo = ((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/TMethod/allTimesAdded_TMethod", fitPass)));

      binInfo.emplace_back(histForInfo->GetBinWidth(1), histForInfo->GetBinLowEdge(1));

      TF1* TMethodFitFunction = (TF1*) ((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/TMethod/allTimesAdded_TMethod", fitPass)))->GetFunction("TmethodFitFunc")->Clone();
      TF1* fullRatioFitFunction = (TF1*) ((TGraphErrors*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass)))->GetFunction("fullRatioFitFunc")->Clone();

      TmethodFits.push_back(TMethodFitFunction);
      fullRatioFits.push_back(fullRatioFitFunction);
    }


  auto TmethodDir = addedDir->mkdir("TMethod");
  makePlots(TmethodDir, TmethodFits, "T Method Fit", binInfo, true);

  auto fullRatioDir = addedDir->mkdir("FullRatio");
  makePlots(fullRatioDir, fullRatioFits, "Full Ratio Fit", binInfo, true);



/////////////////////////////////////////////////////////////////////////////////////

  delete outputFile;

  return 1;
}
