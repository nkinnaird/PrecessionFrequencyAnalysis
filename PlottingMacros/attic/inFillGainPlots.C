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
#include <sstream>
#include <TTree.h>
#include <TNtuple.h>
#include <TPaveStats.h>
#include <TVectorD.h>
#include <THStack.h>

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;

bool saveImagesDirectly = true; // only applies to one scan type at a time - hardcoded right now

void generalScan(TDirectory* inDir, std::tuple<string, string, string> tupleID, std::vector<float> xVector, std::vector<TF1*> fits, bool saveImages)
{
  inDir->cd();

  bool performFits = ((xVector.back() - xVector.front() != 0)) ? true : false;
  bool vsTimeParam = false;

    string fitName = std::get<0>(tupleID);
    string scanType = std::get<1>(tupleID);
    string titleStr = std::get<2>(tupleID);

    string scanName = removeSpaces(scanType);

    if(scanType.find("Lifetime") != string::npos) vsTimeParam = true;

    
/////////////////////////////////////////////////////////////////////////////////////

    TGraphErrors* Chi2_Vs_Val = new TGraphErrors();
      Chi2_Vs_Val->SetName((fitName + "_Chi2_Vs_" + scanName).c_str());
      Chi2_Vs_Val->SetTitle((titleStr + " #chi^{2} Vs " + scanType).c_str());
      Chi2_Vs_Val->GetXaxis()->SetTitle(scanType.c_str());
      Chi2_Vs_Val->GetYaxis()->SetTitle("#chi^{2}");
      Chi2_Vs_Val->GetYaxis()->SetTitleOffset(1.8);

    TGraphErrors* Chi2NDF_Vs_Val = new TGraphErrors();
      Chi2NDF_Vs_Val->SetName((fitName + "_Chi2NDF_Vs_" + scanName).c_str());
      Chi2NDF_Vs_Val->SetTitle((titleStr + " #chi^{2}/NDF Vs " + scanType).c_str());
      Chi2NDF_Vs_Val->GetXaxis()->SetTitle(scanType.c_str());
      Chi2NDF_Vs_Val->GetYaxis()->SetTitle("#chi^{2}/NDF");
      Chi2NDF_Vs_Val->GetYaxis()->SetTitleOffset(1.8);

    TGraphErrors* Pval_Vs_Val = new TGraphErrors();
      Pval_Vs_Val->SetName((fitName + "_Pval_Vs_" + scanName).c_str());
      Pval_Vs_Val->SetTitle((titleStr + " P Value Vs " + scanType).c_str());
      Pval_Vs_Val->GetXaxis()->SetTitle(scanType.c_str());
      Pval_Vs_Val->GetYaxis()->SetTitle("P value");
      Pval_Vs_Val->GetYaxis()->SetTitleOffset(1.8);


    int pointNo = 0;

    for (float xVal : xVector)
    {
      double chi2 = fits.at(pointNo)->GetChisquare();
      double chi2ndf = chi2/fits.at(pointNo)->GetNDF();
      double pVal = fits.at(pointNo)->GetProb();

      Chi2_Vs_Val->SetPoint(pointNo, xVal, chi2);
      Chi2NDF_Vs_Val->SetPoint(pointNo, xVal, chi2ndf);
      Pval_Vs_Val->SetPoint(pointNo, xVal, pVal);

      pointNo++;
    } // end loop over point numbers

    if(vsTimeParam){
      nsTOus(Chi2_Vs_Val, string(Chi2_Vs_Val->GetXaxis()->GetTitle()) + " (#mus)");
      nsTOus(Chi2NDF_Vs_Val, string(Chi2NDF_Vs_Val->GetXaxis()->GetTitle()) + " (#mus)");
      nsTOus(Pval_Vs_Val, string(Pval_Vs_Val->GetXaxis()->GetTitle()) + " (#mus)");
    }

      if(performFits){
        TF1* parabFunc = new TF1("parabFunc", "[2] * (x - [1]) * (x - [1]) + [0]", .8, 1.2);
        parabFunc->SetNpx(10000);
        parabFunc->SetParameter(0,fits.at(0)->GetChisquare());
        parabFunc->SetParameter(1,1);
        parabFunc->SetParameter(2,fits.at(0)->GetChisquare());
        Chi2_Vs_Val->Fit("parabFunc", "Q");
        Chi2_Vs_Val->GetFunction("parabFunc")->SetLineColor(2);
      }

      
      Chi2_Vs_Val->Write();
      Chi2NDF_Vs_Val->Write();
      Pval_Vs_Val->Write();

    if(saveImages){
      string canvasName = (string(Chi2_Vs_Val->GetName()) + "_Canv");
      auto graphCanvas = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,800,600);
      Chi2_Vs_Val->Draw("AP");

      graphCanvas->Update();
      TPaveStats *statsBox = (TPaveStats*)graphCanvas->GetPrimitive("stats");
      statsBox->SetBorderSize(1);
      statsBox->Draw("SAME");

      graphCanvas->SaveAs(("Images/" + canvasName + ".png").c_str());
    }

/////////////////////////////////////////////////////////////////////////////////////

  int numParams = fits.at(0)->GetNpar();

  for (int parNum = 0; parNum < numParams; ++parNum)
  {
    string paramString = fits.at(0)->GetParName(parNum);
    string nameString = removeDangerousCharacters(paramString);

    bool phaseParam = false;
    string strAddon = "";
    if(paramString.compare("R") == 0) strAddon = " (ppm)";
    else if(paramString.find("tau") != string::npos) strAddon = " (#mus)";
    else if(paramString.find("omega") != string::npos) strAddon = " (rad/#mus)";
    else if(paramString.find("phi") != string::npos) phaseParam = true;

    TGraphErrors* param_Vs_Val = new TGraphErrors();
    param_Vs_Val->SetName((fitName + "_" + nameString + "_Vs_" + scanName).c_str());
    param_Vs_Val->SetTitle((titleStr + " " + paramString + " Vs " + scanType).c_str());

      for (uint fitNum = 0; fitNum < fits.size(); ++fitNum)
      {
        double parameterValue = fits.at(fitNum)->GetParameter(parNum);
        double parameterError = fits.at(fitNum)->GetParError(parNum);
        
        if(phaseParam) normalizePhase(parameterValue);

        param_Vs_Val->SetPoint(fitNum, xVector.at(fitNum), parameterValue);
        param_Vs_Val->SetPointError(fitNum, 0, parameterError);
      }

    param_Vs_Val->GetXaxis()->SetTitle(scanType.c_str());
    param_Vs_Val->GetYaxis()->SetTitle((paramString + strAddon).c_str());
    param_Vs_Val->GetYaxis()->SetTitleOffset(1.8);

    if(vsTimeParam) nsTOus(param_Vs_Val, string(param_Vs_Val->GetXaxis()->GetTitle()) + " (#mus)");

    if(paramString.compare("R") == 0 && performFits){
      param_Vs_Val->Fit("pol1", "Q"); // fit for the systematic effect on R
      param_Vs_Val->GetFunction("pol1")->SetLineColor(2);
    }

    param_Vs_Val->Write();

    if(saveImages && paramString.compare("R") == 0){
      string canvasName = (fitName + "_" + nameString + "_Vs_" + scanName + "_Canv");
      auto graphCanvas = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,800,600);
      param_Vs_Val->Draw("AP");

      graphCanvas->Update();
      TPaveStats *statsBox = (TPaveStats*)graphCanvas->GetPrimitive("stats");
      statsBox->SetBorderSize(1);

      // statsBox->SetX2(graphCanvas->GetUxmax()); // can reposition and by extension resize the stats box with these commands - not the best though, there must be a better way
      // statsBox->SetY2(graphCanvas->GetUymax());

      statsBox->Draw("SAME");
      graphCanvas->SaveAs(("Images/" + canvasName + ".png").c_str());
    }

  }

}


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


int inFillGainPlots(std::string filePath)
{
  // gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  TFile* outputFile = new TFile("inFillGainPlots.root","RECREATE");
  auto topDir = outputFile->mkdir("topDir");


/////////////////////////////////////////////////////////////////////////////////////
  // These only get set for the interactive root session (any generated canvases, etc.), but does not apply to the output root file - that comes from .rootlogon.C
  gStyle->SetOptStat(000000);
  gStyle->SetOptTitle(1);
  gStyle->SetOptFit(2);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1);
  gStyle->SetLineColor(1);

  gStyle->SetPadRightMargin(.05);

  // gStyle->SetOptStat("m");


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  auto addedDir = topDir->mkdir("Added");

/////////////////////////////////////////////////////////////////////////////////////

// values scanned over

std::vector<float> amplitudes;
std::vector<float> lifetimes;

std::vector<TF1*> fullRatioFits;

/////////////////////////////////////////////////////////////////////////////////////

      TNtuple *firstTuple = (TNtuple*)inputFile->Get(Form("topDir/FitPasses/FitPass0/FitConditions0"));
      float tP;
      firstTuple->SetBranchAddress("totalPasses", &tP);
      firstTuple->GetEntry(0);

/////////////////////////////////////////////////////////////////////////////////////

    // loop through passes and store relevant objects

    for (int fitPass = 0; fitPass < tP; ++fitPass)
    {
      TF1* gainSagFunction = (TF1*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/gainSagFunction", fitPass));
      amplitudes.push_back(gainSagFunction->GetParameter(0));
      lifetimes.push_back(gainSagFunction->GetParameter(1));

      TF1* fullRatioFitFunction = (TF1*) ((TGraphErrors*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass)))->GetFunction("fullRatioFitFunc")->Clone();
      fullRatioFits.push_back(fullRatioFitFunction);
    }

/////////////////////////////////////////////////////////////////////////////////////

    auto added_FullRatio_dir = addedDir->mkdir("FullRatio");

    auto added_FullRatio_amp_dir = added_FullRatio_dir->mkdir("InFillGainAmplitude");
    std::tuple<string, string, string> ampInfo("FullRatio", "In Fill Gain Amplitude", "Full Ratio Fit");
    generalScan(added_FullRatio_amp_dir, ampInfo, amplitudes, fullRatioFits, saveImagesDirectly);

    auto added_FullRatio_tau_dir = added_FullRatio_dir->mkdir("InFillGainLifetime");
    std::tuple<string, string, string> tauInfo("FullRatio", "In Fill Gain Lifetime", "Full Ratio Fit");
    generalScan(added_FullRatio_tau_dir, tauInfo, lifetimes, fullRatioFits, false);

/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

      delete outputFile;


  return 1;
}
