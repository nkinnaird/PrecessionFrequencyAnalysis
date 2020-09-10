// 3-31-20: Macro for plots for scans over fit window. Was just done and used once, but may be needed again. Very similar to fitStartScan.C.

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

bool saveImagesDirectly = false;

void fitWindowScanPlots(TDirectory* inDir, std::vector<pair<float, float> > fitRanges, std::vector<TF1*> fits, string titleStr, bool saveImages)
{
  inDir->cd();

  int fitWindow = fitRanges.at(0).second - fitRanges.at(0).first;
  fitWindow /= 1000;
  cout << "Fit window is: " << fitWindow << " mus." << endl;

  string scanType = to_string(fitWindow) + " #mus Fit Window Start Time";
  string nameChar = "W_" + to_string(fitWindow) + "mus";

  string xAxisTitle = scanType + " [#mus]";
  string fitTitle = inDir->GetName();

  int numParams = fits.at(0)->GetNpar();

/////////////////////////////////////////////////////////////////////////////////////

  TGraphErrors* Chi2NDF_Vs_FS = new TGraphErrors();
  Chi2NDF_Vs_FS->SetName((fitTitle + "_Chi2NDF_Vs_F" + nameChar).c_str());
  // Chi2NDF_Vs_FS->SetTitle((titleStr + " #chi^{2}/NDF Vs " + scanType).c_str());
  Chi2NDF_Vs_FS->SetTitle(titleStr.c_str());
  Chi2NDF_Vs_FS->GetXaxis()->SetTitle(xAxisTitle.c_str());
  Chi2NDF_Vs_FS->GetYaxis()->SetTitle("#chi^{2}/NDF");

  TGraph* Pval_Vs_FS = new TGraph();
  Pval_Vs_FS->SetName((fitTitle + "_Pval_Vs_F" + nameChar).c_str());
  // Pval_Vs_FS->SetTitle((titleStr + " P Value Vs " + scanType).c_str());
  Pval_Vs_FS->SetTitle(titleStr.c_str());
  Pval_Vs_FS->GetXaxis()->SetTitle(xAxisTitle.c_str());
  Pval_Vs_FS->GetYaxis()->SetTitle("P value");

  for (uint fitNum = 0; fitNum < fits.size(); ++fitNum)
  {
    float fitTimePoint = fitRanges.at(fitNum).first;

    Chi2NDF_Vs_FS->SetPoint(fitNum, fitTimePoint, fits.at(fitNum)->GetChisquare()/fits.at(fitNum)->GetNDF());
    Chi2NDF_Vs_FS->SetPointError(fitNum, 0, sqrt(2./fits.at(fitNum)->GetNDF()));

    Pval_Vs_FS->SetPoint(fitNum, fitTimePoint, fits.at(fitNum)->GetProb());
  }

  nsTOus(Chi2NDF_Vs_FS, xAxisTitle);

  string canvName = (fitTitle + "_Chi2NDF_Vs_F" + nameChar + "_canv");
  auto chi2ndfcanv = new TCanvas(canvName.c_str(), canvName.c_str(), 200, 10, 600, 400);
  Chi2NDF_Vs_FS->Draw("AP");

  chi2ndfcanv->Write();
  if(saveImages) chi2ndfcanv->SaveAs(("Images/" + canvName + datasetTagForPlots + ".png").c_str());

  string Pval_canvName = (fitTitle + "_Pval_Vs_F" + nameChar + "_canv");
  auto pvalcanv = new TCanvas(Pval_canvName.c_str(), Pval_canvName.c_str(), 200, 10, 600, 400);

  nsTOus(Pval_Vs_FS, xAxisTitle);
  Pval_Vs_FS->Draw("AP");

  Pval_Vs_FS->Write();
  if(saveImages) pvalcanv->SaveAs(("Images/" + string(Pval_Vs_FS->GetName()) + datasetTagForPlots + ".png").c_str());

/////////////////////////////////////////////////////////////////////////////////////

  // figure out phase and asymmetry parameter numbers since those are needed for all plots when constructing the Kawall bands
  int asymmetryParamNum = -1;
  int phaseParamNum = -1;

  for (int parNum = 0; parNum < numParams; ++parNum)
  {
    string paramString = fits.at(0)->GetParName(parNum);
    if(paramString.compare("A") == 0) asymmetryParamNum = parNum;
    else if(paramString.compare("#phi") == 0) phaseParamNum = parNum;
  }

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

    TGraphErrors* param_Vs_FS = new TGraphErrors();
    param_Vs_FS->SetName((fitTitle + "_" + nameString + "_Vs_F" + nameChar).c_str());
    // param_Vs_FS->SetTitle((titleStr + " " + paramString + " Vs " + scanType).c_str());
    param_Vs_FS->SetTitle(titleStr.c_str());

    double param_at_start = fits.at(0)->GetParameter(parNum);
    double error_at_start = fits.at(0)->GetParError(parNum);

    double g2Phase_at_start = fits.at(0)->GetParameter(phaseParamNum);
    double asymmetry_at_start = fits.at(0)->GetParameter(asymmetryParamNum);

      if(phaseParam) normalizePhase(param_at_start);
      normalizePhase(g2Phase_at_start);

      int kawallBandPointNo = 0;
      for (uint fitNum = 0; fitNum < fits.size(); ++fitNum)
      {
        float fitStart = fitRanges.at(fitNum).first;
        float fitEnd = fitRanges.at(fitNum).second;
        float fitTimePoint = fitStart;

        double parameterValue = fits.at(fitNum)->GetParameter(parNum);
        double parameterError = fits.at(fitNum)->GetParError(parNum);

        double g2Phase = fits.at(fitNum)->GetParameter(phaseParamNum);
        double asymmetry = fits.at(fitNum)->GetParameter(asymmetryParamNum);

          if(phaseParam) {
            normalizePhase(parameterValue);
              if(fitNum > 0){
                double previousPhaseParameter = fits.at(fitNum-1)->GetParameter(parNum);
                normalizePhase(previousPhaseParameter);
                if(parameterValue < (previousPhaseParameter - pi)) parameterValue += 2*pi; // to clean up plots with phase parameter
              }
          }
          normalizePhase(g2Phase);


        param_Vs_FS->SetPoint(fitNum, fitTimePoint, parameterValue);
        param_Vs_FS->SetPointError(fitNum, 0, parameterError);

      }

    nsTOus(param_Vs_FS, xAxisTitle);

    if(paramString.compare("R") == 0){
      param_Vs_FS->Fit("pol1", "Q");
      param_Vs_FS->GetFunction("pol1")->SetLineColor(2);
    }

    string canvasName = (fitTitle + "_" + nameString + "_F" + nameChar + "_Canv");
    auto fitStart_canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 200, 10, 600, 400);
    param_Vs_FS->GetXaxis()->SetTitle(xAxisTitle.c_str());
    param_Vs_FS->GetYaxis()->SetTitle((paramString + strAddon).c_str());
    param_Vs_FS->GetYaxis()->SetTitleOffset(1.8);
    param_Vs_FS->Draw("AP");

    if(paramString.compare("R") == 0){
      fitStart_canvas->Update();
      TPaveStats *statsBox = (TPaveStats*)fitStart_canvas->GetPrimitive("stats");
      statsBox->SetBorderSize(1);
      statsBox->Draw("SAME");
    }

    fitStart_canvas->Write();
    if(saveImages) fitStart_canvas->SaveAs(("Images/" + canvasName + datasetTagForPlots + ".png").c_str());
    // delete fitStart_canvas;
    
    if(1 && paramString.compare("R") == 0) fitStart_canvas->SaveAs(("Images/" + canvasName + datasetTagForPlots + ".png").c_str());

  }

}

/////////////////////////////////////////////////////////////////////////////////////

int fitWindowScan(std::string filePath)
{
  // gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen
  // when setting this to true, for some reason there is an extra space in the axis title with the micro symbol in microseconds in the saved images - no idea why..
  // slower to run but might want to set to false when making final images

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  // check if dataset tag exists in file and if so append it to the file name and write it to the file

  string outputFileName = "fitWindowScan";
  TNamed* tag = applyDatasetTag(inputFile, outputFileName);

  TFile* outputFile = new TFile((outputFileName + ".root").c_str(),"RECREATE");
  if(tag) tag->Write();
  auto topDir = outputFile->mkdir("topDir");

/////////////////////////////////////////////////////////////////////////////////////
  // These only get set for the interactive root session (any generated canvases, etc.), but does not apply to the output root file - that comes from .rootlogon.C
  gStyle->SetOptStat(000000);
  // gStyle->SetOptTitle(0);
  // gStyle->SetOptFit(0000);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1);
  gStyle->SetLineColor(1);

  gStyle->SetPadRightMargin(.05);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
      // Plots over different fit conditions

  auto addedDir = topDir->mkdir("Added");


std::vector<pair<float, float> > fitRanges;
std::vector<TF1*> fiveParamFits;
std::vector<TF1*> TmethodFits;
std::vector<TF1*> threeParRatioFits;
std::vector<TF1*> fullRatioFits;


      TNtuple *firstTuple = (TNtuple*)inputFile->Get(Form("topDir/FitPasses/FitPass0/FitConditions0"));
      float tP;
      firstTuple->SetBranchAddress("totalPasses", &tP);
      firstTuple->GetEntry(0);

    for (int fitPass = 0; fitPass < tP; ++fitPass)
    {
      float fitStart;
      float fitEnd;

      TNtuple *ntuple = (TNtuple*)inputFile->Get(Form("topDir/FitPasses/FitPass%d/FitConditions%d", fitPass, fitPass));
      ntuple->SetBranchAddress("fitStartTime", &fitStart);
      ntuple->SetBranchAddress("fitEndTime", &fitEnd);
      ntuple->GetEntry(0);

      TF1* TMethodFitFunction = (TF1*) ((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/TMethod/allTimesAdded_TMethod", fitPass)))->GetFunction("TmethodFitFunc")->Clone();
      TF1* fullRatioFitFunction = (TF1*) ((TGraphErrors*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass)))->GetFunction("fullRatioFitFunc")->Clone();

      fitRanges.emplace_back(fitStart, fitEnd);
      TmethodFits.push_back(TMethodFitFunction);
      fullRatioFits.push_back(fullRatioFitFunction);
    }


  auto TmethodDir = addedDir->mkdir("TMethod");
  fitWindowScanPlots(TmethodDir, fitRanges, TmethodFits, "T Method Fit", saveImagesDirectly);

  auto fullRatioDir = addedDir->mkdir("FullRatio");
  fitWindowScanPlots(fullRatioDir, fitRanges, fullRatioFits, "Full Ratio Fit", saveImagesDirectly);

/////////////////////////////////////////////////////////////////////////////////////

  delete outputFile;

  return 1;
}
