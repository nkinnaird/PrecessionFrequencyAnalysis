// 3-31-20: Macro for plots for scans over general iterations. Can be used when anything is scanned over. 
// Two special cases include scans over random seed (which is basically just a general scan), or a scan vs any fixed parameter in the fits (just for ratio fits).
// The plots for the scan over the fixed parameter can be chosen by setting plotVsParam to the integer value corresponding the parameter of interest.

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

int plotVsParam = -1; // set to a positive int to plot against a specific parameter, otherwise set to a negative number to plot vs iteration/random seed, only applies to ratio fit right now - could make it apply to others but then the parameter number is different for the different fits
bool saveImagesDirectly = true;
bool makeHistogramPlotsFromGraphs = true;
bool vsRandSeed = true;

void iterationPlots(TDirectory* inDir, std::vector<TF1*> fits, string titleStr, int VsParamNum, bool saveImages)
{
  auto graphDir = inDir->mkdir("Graphs");
  auto histDir = inDir->mkdir("Hists");

  graphDir->cd();

  string plotName;
  string xAxisTitle;
  bool vsTimeParam = false;

  if(VsParamNum < 0){
    if(vsRandSeed) xAxisTitle = "Random Seed";
    else xAxisTitle = "Fit Iteration";
    plotName = "Iter";
  } 
  else{
    xAxisTitle = fits.at(0)->GetParName(VsParamNum);
    plotName = removeDangerousCharacters(xAxisTitle);
    if(xAxisTitle.find("tau") != string::npos) vsTimeParam = true;
  }

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
    double xPoint = (VsParamNum < 0) ? fitNum+1 : fits.at(fitNum)->GetParameter(VsParamNum);

    Chi2NDF_Vs_Iter->SetPoint(fitNum, xPoint, fits.at(fitNum)->GetChisquare()/fits.at(fitNum)->GetNDF());
    Chi2NDF_Vs_Iter->SetPointError(fitNum, 0, sqrt(2./fits.at(fitNum)->GetNDF()));

    Pval_Vs_Iter->SetPoint(fitNum, xPoint, fits.at(fitNum)->GetProb());
  }

  if(vsTimeParam){
    Chi2NDF_Vs_Iter->GetXaxis()->SetTitle((xAxisTitle + " [#mus]").c_str());
    Pval_Vs_Iter->GetXaxis()->SetTitle((xAxisTitle + " [#mus]").c_str());
  }

  graphDir->cd();
    Chi2NDF_Vs_Iter->Write();
    Pval_Vs_Iter->Write();

  histDir->cd();
    TH1F* Chi2NDF_Vs_Iter_hist = projectGraphToHist(Chi2NDF_Vs_Iter);
      Chi2NDF_Vs_Iter_hist->Write();
    TH1F* Pval_Vs_Iter_hist = projectGraphToHist(Pval_Vs_Iter);
      Pval_Vs_Iter_hist->Write();

  if(saveImages){
    string canvasName = (fitTitle + "_Chi2NDF_Vs_" + plotName + "_Canv");
    auto chi2_canv = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,500,400);
    if(VsParamNum < 0) Chi2NDF_Vs_Iter->GetXaxis()->SetLimits(0, fits.size()+1);
    Chi2NDF_Vs_Iter->Draw("AP");

      double mean = Chi2NDF_Vs_Iter->GetMean(2);
      double rms = Chi2NDF_Vs_Iter->GetRMS(2);
      TPaveText* textBox = new TPaveText(0.75,0.7,0.975,0.9, "NDC");
      textBox->AddText(Form("Mean  %0.4g", mean));
      textBox->AddText(Form("RMS  %0.4g", rms));
      textBox->SetBorderSize(1);
      textBox->SetFillStyle(1001);
      textBox->Draw("SAME");

    chi2_canv->SaveAs(("Images/" + canvasName + datasetTagForPlots + ".png").c_str());

    auto chi2_hist_canv = new TCanvas((canvasName + "_hist").c_str(),canvasName.c_str(),200,10,500,400);
    Chi2NDF_Vs_Iter_hist->Draw();
      textBox->Clear();
      textBox->AddText(Form("Entries  %zu", fits.size()));
      textBox->AddText(Form("Mean  %0.4g", mean));
      textBox->AddText(Form("RMS  %0.4g", rms));
      textBox->Draw("SAME");

    chi2_hist_canv->SaveAs(("Images/" + canvasName + "_hist" + datasetTagForPlots + ".png").c_str());

    canvasName = fitTitle + "_PVal_Vs_" + plotName + "_Canv";
    auto pval_canv = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,500,400);
    if(VsParamNum < 0) Pval_Vs_Iter->GetXaxis()->SetLimits(0, fits.size()+1);
    Pval_Vs_Iter->Draw("AP");
    pval_canv->SaveAs(("Images/" + canvasName + datasetTagForPlots + ".png").c_str());
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

    TGraphErrors* param_Vs_Iter = new TGraphErrors();
    param_Vs_Iter->SetName((fitTitle + "_" + nameString + "_Vs_" + plotName).c_str());
    param_Vs_Iter->SetTitle((titleStr + " " + paramString + " Vs " + xAxisTitle).c_str());
    param_Vs_Iter->GetYaxis()->SetTitle((paramString + strAddon).c_str());
    param_Vs_Iter->GetYaxis()->SetTitleOffset(2.0);
    param_Vs_Iter->GetXaxis()->SetTitle(xAxisTitle.c_str());
      
      for (uint fitNum = 0; fitNum < fits.size(); ++fitNum)
      {
        double xPoint = (VsParamNum < 0) ? fitNum+1 : fits.at(fitNum)->GetParameter(VsParamNum);

        double parameterValue = fits.at(fitNum)->GetParameter(parNum);
        double parameterError = fits.at(fitNum)->GetParError(parNum);
        
        if(phaseParam) normalizePhase(parameterValue);
        // if(phaseParam && parameterValue < 2) parameterValue += 2*pi;

        param_Vs_Iter->SetPoint(fitNum, xPoint, parameterValue);
        param_Vs_Iter->SetPointError(fitNum, 0, parameterError);
      }

    if(vsTimeParam) param_Vs_Iter->GetXaxis()->SetTitle((xAxisTitle + " [#mus]").c_str());
    
    if(paramString.compare("R") == 0 && VsParamNum >= 0){
      param_Vs_Iter->Fit("pol1", "Q"); // fit for the systematic effect on R
      param_Vs_Iter->GetFunction("pol1")->SetLineColor(2);

/////////////////////////////////////////////////////////////////////////////////////
/*
      // for fitting against the CBO lifetime which is non-linear
      TF1* expFunc = new TF1("expFunc", "[0] * exp(-x/[1]) + [2]", fits.front()->GetParameter(VsParamNum)/1000., fits.back()->GetParameter(VsParamNum)/1000.); // convert from ns to us because I've scaled the x axis in time
      expFunc->SetParameter(0, 1);
      expFunc->SetParameter(1, fits.front()->GetParameter(VsParamNum)/1000./3.);
      expFunc->SetParameter(2, fits.front()->GetParameter(1));
      param_Vs_Iter->Fit(expFunc, "QR"); // fit for the systematic effect on R
      param_Vs_Iter->GetFunction("expFunc")->SetLineColor(2);

      double fixedLifetime = 180.; // in us
      double lifetimeError = 16.;

      double R_fixed = expFunc->Eval(fixedLifetime);
      double R_above = expFunc->Eval(fixedLifetime + lifetimeError);
      double R_below = expFunc->Eval(fixedLifetime - lifetimeError);

      double errorR_above = R_above - R_fixed;
      double errorR_below = R_fixed - R_below;

      cout << "values of R: " << R_fixed << " " << R_above << " " << R_below << endl;
      cout << "systematic errors on R: " << errorR_above << " " << errorR_below << endl;
*/      
    }

  graphDir->cd();
    param_Vs_Iter->Write();

  histDir->cd();
    TH1F* param_Vs_Iter_hist = projectGraphToHist(param_Vs_Iter);
    if(paramString.compare(N_string) == 0 || paramString.compare(A_string) == 0 || paramString.compare(phi_string) == 0) param_Vs_Iter_hist->GetXaxis()->SetNdivisions(505);
      param_Vs_Iter_hist->Write();


    if(saveImages){
      string canvasName = (fitTitle + "_" + nameString + "_Vs_" + plotName + "_Canv");
      auto graphCanvas = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,500,400);
      if(VsParamNum < 0) param_Vs_Iter->GetXaxis()->SetLimits(0, fits.size()+1);
      param_Vs_Iter->Draw("AP");

      if(VsParamNum < 0){
        double mean = param_Vs_Iter->GetMean(2);
        double rms = param_Vs_Iter->GetRMS(2);

        TPaveText* textBox = new TPaveText(0.75,0.7,0.975,0.9, "NDC");
        textBox->AddText(Form("Mean  %0.4g", mean));
        textBox->AddText(Form("RMS  %0.4g", rms));
        textBox->SetBorderSize(1);
        textBox->SetFillStyle(1001);
        textBox->Draw("SAME");

        graphCanvas->SaveAs(("Images/" + canvasName + datasetTagForPlots + ".png").c_str());

        if(makeHistogramPlotsFromGraphs){
          auto histCanvas = new TCanvas((canvasName + "_hist").c_str(),canvasName.c_str(),200,10,500,400);
            param_Vs_Iter_hist->Draw();
            textBox->Clear();
            textBox->AddText(Form("Entries  %zu", fits.size()));
            textBox->AddText(Form("Mean  %0.4g", mean));
            textBox->AddText(Form("RMS  %0.4g", rms));
            textBox->Draw("SAME");

            histCanvas->SaveAs(("Images/" + canvasName + "_hist" + datasetTagForPlots + ".png").c_str());
        } // end if make histogram images from graphs

      }
      else if(VsParamNum >= 0 && paramString.compare("R") == 0){
        
        // draw options so shape of R vs parameter curve shows up
        double lowRange = (fits.front()->GetParameter(1) < fits.back()->GetParameter(1)) ? fits.front()->GetParameter(1) : fits.back()->GetParameter(1);
        double highRange = (fits.front()->GetParameter(1) > fits.back()->GetParameter(1)) ? fits.front()->GetParameter(1) : fits.back()->GetParameter(1);

        param_Vs_Iter->GetYaxis()->SetRangeUser(lowRange - .02, highRange + .02);
        param_Vs_Iter->Draw("XAP"); // X - draw without error bars

        graphCanvas->Update();
        TPaveStats *statsBox = (TPaveStats*)graphCanvas->GetPrimitive("stats");
        statsBox->SetBorderSize(1);
        statsBox->Draw("SAME");

        graphCanvas->SaveAs(("Images/" + canvasName + datasetTagForPlots + ".png").c_str());
      } // end else if
    } // end if save images

  } // end loop over parameters

}

/////////////////////////////////////////////////////////////////////////////////////

int generalIterPlotter(std::string filePath)
{
  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  vector<string> fileVector;

  size_t found = filePath.find_last_of(".");
  string extension = filePath.substr(found);

  if(extension == ".root") fileVector.push_back(filePath);
  else if(extension == ".txt"){
    ifstream inList(filePath);
    string path;
    while(inList >> path) fileVector.push_back(path);
  }
  else{
    printf("Not passing files in correctly. Should be either a single .root file or a .txt file with multiple .root files.\n");
    return 0;
  }

/////////////////////////////////////////////////////////////////////////////////////

  TFile *lastFile = TFile::Open(fileVector.back().c_str());
   if (lastFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  // check if dataset tag exists in file and if so append it to the file name and write it to the file

  string outputFileName = "generalPlotsVsIter";
  TNamed* tag = applyDatasetTag(lastFile, outputFileName);

  lastFile->Close();

/////////////////////////////////////////////////////////////////////////////////////

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
  auto calosDir = topDir->mkdir("Calos");


std::vector<TF1*> fiveParamFits;
std::vector<TF1*> TmethodFits;
std::vector<TF1*> threeParRatioFits;
std::vector<TF1*> fullRatioFits;


for (uint fileNum = 0; fileNum < fileVector.size(); ++fileNum)
{
  string fileName = fileVector.at(fileNum);

  // pull in input file
  TFile *inputFile = TFile::Open(fileName.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      continue;
   }

      TNtuple *firstTuple = (TNtuple*)inputFile->Get(Form("topDir/FitPasses/FitPass0/FitConditions0"));
      float tP;
      firstTuple->SetBranchAddress("totalPasses", &tP);
      firstTuple->GetEntry(0);

    for (int fitPass = 0; fitPass < tP; ++fitPass)
    {
      TF1* TMethodFitFunction = (TF1*) ((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/TMethod/allTimesAdded_TMethod", fitPass)))->GetFunction("TmethodFitFunc")->Clone();
      TmethodFits.push_back(TMethodFitFunction);

      TGraphErrors* ratioGraph = (TGraphErrors*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass));
      TF1* fullRatioFitFunction = 0;
      if(ratioGraph != 0){ // check to make sure ratio fits exist before trying to make plots with them
        fullRatioFitFunction = (TF1*) ratioGraph->GetFunction("fullRatioFitFunc")->Clone();
        fullRatioFits.push_back(fullRatioFitFunction);
      }      
    }

    inputFile->Close();

} // end loop over multiple files


if(TmethodFits.size() > 0){
  auto TmethodDir = addedDir->mkdir("TMethod");
  iterationPlots(TmethodDir, TmethodFits, "T Method Fit", plotVsParam, saveImagesDirectly);
}

if(fullRatioFits.size() > 0){
  auto fullRatioDir = addedDir->mkdir("FullRatio");
  iterationPlots(fullRatioDir, fullRatioFits, "Full Ratio Fit", plotVsParam, saveImagesDirectly);
}

/////////////////////////////////////////////////////////////////////////////////////

  delete outputFile;

  return 1;
}
