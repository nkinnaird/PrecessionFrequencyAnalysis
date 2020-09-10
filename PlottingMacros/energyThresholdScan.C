// 3-31-20: Macro for plots for scans over energy threshold. Will find the optimum energy threshold for all calorimeters added together.

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
#include <TLegend.h>
#include <TLine.h>
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

void thresholdScanPlots(TDirectory* inDir, std::vector<float> eTs, std::vector<TF1*> fits, string titleStr)
{
  inDir->cd();

  string plotName;
  string xAxisTitle;
  bool vsTimeParam = false;

  xAxisTitle = "Energy Threshold [MeV]";
  plotName = "ETh";

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
    double xPoint = eTs.at(fitNum);

    Chi2NDF_Vs_Iter->SetPoint(fitNum, xPoint, fits.at(fitNum)->GetChisquare()/fits.at(fitNum)->GetNDF());
    Chi2NDF_Vs_Iter->SetPointError(fitNum, 0, sqrt(2./fits.at(fitNum)->GetNDF()));

    Pval_Vs_Iter->SetPoint(fitNum, xPoint, fits.at(fitNum)->GetProb());
  }

  Chi2NDF_Vs_Iter->Write();
  Pval_Vs_Iter->Write();

/////////////////////////////////////////////////////////////////////////////////////

  // figure out phase and asymmetry parameter numbers for the Kawall bands
  int asymmetryParamNum = -1;
  int phaseParamNum = -1;

  for (int parNum = 0; parNum < numParams; ++parNum)
  {
    string paramString = fits.at(0)->GetParName(parNum);
    if(paramString.compare("A") == 0) asymmetryParamNum = parNum;
    else if(paramString.compare("#phi") == 0) phaseParamNum = parNum;
  }

  double g2Phase_at_start = fits.at(0)->GetParameter(phaseParamNum);
  double asymmetry_at_start = fits.at(0)->GetParameter(asymmetryParamNum);

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
        double xPoint = eTs.at(fitNum);

        double parameterValue = fits.at(fitNum)->GetParameter(parNum);
        double parameterError = fits.at(fitNum)->GetParError(parNum);
        
        if(phaseParam) normalizePhase(parameterValue);

        param_Vs_Iter->SetPoint(fitNum, xPoint, parameterValue);
        param_Vs_Iter->SetPointError(fitNum, 0, parameterError);
      }

    param_Vs_Iter->Write();


    if(paramString.compare("R") == 0)
    {

      // make R plot (some code copied here to make Kawall bands but only for R)

    TGraph* KawallBandOne = new TGraph();
    TGraph* KawallBandTwo = new TGraph();
      int kawallBandPointNo = 0;

    double param_at_start = fits.at(0)->GetParameter(parNum);
    double error_at_start = fits.at(0)->GetParError(parNum);

      for (uint fitNum = 0; fitNum < fits.size(); ++fitNum)
      {
        double xPoint = eTs.at(fitNum);

        double g2Phase = fits.at(fitNum)->GetParameter(phaseParamNum);
        double asymmetry = fits.at(fitNum)->GetParameter(asymmetryParamNum);
        double analyzingPower = (2. * (asymmetry_at_start/asymmetry) * cos(g2Phase_at_start-g2Phase) - 1);

        double parameterError = fits.at(fitNum)->GetParError(parNum);

        double sigmaDiff = sqrt(pow(parameterError, 2) - pow(error_at_start, 2)*analyzingPower);
        // double sigmaDiff = sqrt(pow(parameterError, 2) - pow(error_at_start, 2));

          if(!isnan(sigmaDiff))
          {
            KawallBandOne->SetPoint(kawallBandPointNo, xPoint, param_at_start + sigmaDiff);
            KawallBandTwo->SetPoint(kawallBandPointNo, xPoint, param_at_start - sigmaDiff);
            kawallBandPointNo++;
          }
          else cout << "nan point in bands for R parameter at point " << fitNum << endl;
      }

      string canvasName = (fitTitle + "_R_Vs_" + plotName);

      auto R_canv = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,500,400);
      param_Vs_Iter->Draw("AP");
          KawallBandOne->Draw("SAME");
          KawallBandTwo->Draw("SAME");
      string graphName = param_Vs_Iter->GetName();
      R_canv->SaveAs(("Images/" + graphName + datasetTagForPlots + ".png").c_str());


/////////////////////////////////////////////////////////////////////////////////////

      // plot error on R as a function of threshold

      TGraphErrors* sigmaR = new TGraphErrors();

      for (uint fitNum = 0; fitNum < fits.size(); ++fitNum)
      {
        sigmaR->SetPoint(fitNum, eTs.at(fitNum), fits.at(fitNum)->GetParError(parNum));
        sigmaR->SetPointError(fitNum, 0, 0);
      }

        // TF1* parabFunc = new TF1("parabFunc", "[2] * (x - [1]) * (x - [1]) + [0]", 1600, 1800);
        // parabFunc->SetNpx(10000);
        // parabFunc->SetParameter(0, fits.at(0)->GetParError(0));
        // parabFunc->SetParameter(1, 1700);
        // parabFunc->SetParLimits(1, 1650, 1750);
        // // parabFunc->SetParameter(2, fits.at(0)->GetParError(0));
        // parabFunc->SetParameter(2, 4e-7);
        // parabFunc->SetParLimits(2, 3e-7, 5e-7);
        // sigmaR->Fit("parabFunc", "QR");
        // sigmaR->GetFunction("parabFunc")->SetLineColor(2);
        // cout << "Minimum: " << sigmaR->GetFunction("parabFunc")->GetMinimumX(1600, 1800) << endl;

        // TF1* myPolFunc = new TF1("myPolFunc", "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x + [7]*x*x*x*x*x*x*x", 1200, 2200);
        TF1* myPolFunc = new TF1("myPolFunc", "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x", 1200, 2200);
        
        // myPolFunc->SetParameter(0,      3.75756);
        // myPolFunc->SetParameter(1,   -0.0113081);
        // myPolFunc->SetParameter(2,   2.1821e-05); // these copied from TBrowser fit panel which seems to do just fine when fitting the function as opposed to just fitting from here
        // myPolFunc->SetParameter(3, -2.39895e-08);
        // myPolFunc->SetParameter(4,  1.55061e-11);
        // myPolFunc->SetParameter(5, -5.83892e-15);
        // myPolFunc->SetParameter(6,  1.18509e-18);
        // myPolFunc->SetParameter(7, -9.91002e-23);
        
        sigmaR->Fit("myPolFunc", "Q");
        sigmaR->Fit("myPolFunc", "Q");
        sigmaR->GetFunction("myPolFunc")->SetLineColor(2);

      double minimum = sigmaR->GetFunction("myPolFunc")->GetMinimumX(1600, 1800);
        cout << "Minimum: " << minimum << endl;

      canvasName = (fitTitle + "_sigmaR_Vs_" + plotName);

      sigmaR->GetYaxis()->SetTitle("#sigma_{R} [ppm]");
      sigmaR->SetName(canvasName.c_str());

      sigmaR->GetYaxis()->SetTitleOffset(2.0);
      sigmaR->GetXaxis()->SetTitle(xAxisTitle.c_str());

      sigmaR->Write();

      auto plotCanvas = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,500,400);
      sigmaR->Draw("AP");

      TPaveText* textBox = new TPaveText(0.75,0.7,0.975,0.8, "NDC");
      textBox->SetBorderSize(1);
      textBox->SetFillStyle(1001);
      textBox->AddText(Form("Minimum: %1.0f", minimum));
      textBox->Draw("SAME");

      plotCanvas->SaveAs(("Images/" + canvasName + datasetTagForPlots + ".png").c_str());

    }
    else if(paramString.compare(A_string) == 0)
    {
      // make A plot
      string canvasName = (fitTitle + "_A_Vs_" + plotName);

      auto A_canv = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,500,400);
      param_Vs_Iter->Draw("AP");
      string graphName = param_Vs_Iter->GetName();
      A_canv->SaveAs(("Images/" + graphName + datasetTagForPlots + ".png").c_str());
    }
    else if(paramString.compare(N_string) == 0)
    {
      // make N plot

      string canvasName = (fitTitle + "_N0_Vs_" + plotName);

      auto N_canv = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,500,400);
      param_Vs_Iter->Draw("AP");
      string graphName = param_Vs_Iter->GetName();
      N_canv->SaveAs(("Images/" + graphName + datasetTagForPlots + ".png").c_str());

      // plot NA^2

      TGraphErrors* NA2_graph = new TGraphErrors();

      for (uint fitNum = 0; fitNum < fits.size(); ++fitNum)
      {
        double NA2 = fits.at(fitNum)->GetParameter(0) * fits.at(fitNum)->GetParameter(2) * fits.at(fitNum)->GetParameter(2);

        NA2_graph->SetPoint(fitNum, eTs.at(fitNum), NA2);
        NA2_graph->SetPointError(fitNum, 0, 0);
      }
/*
      double firstNA2 = fits.at(0)->GetParameter(0) * fits.at(0)->GetParameter(2) * fits.at(0)->GetParameter(2);

        TF1* parabFunc = new TF1("parabFunc", "[2] * (x - [1]) * (x - [1]) + [0]", 1600, 1800);
        parabFunc->SetNpx(10000);
        parabFunc->SetParameter(0, firstNA2);
        parabFunc->SetParameter(1, 1700);
        parabFunc->SetParLimits(1, 1650, 1750);
        parabFunc->SetParameter(2, -firstNA2);
        // parabFunc->SetParameter(2, -1.026);
        NA2_graph->Fit("parabFunc", "QR");
        NA2_graph->GetFunction("parabFunc")->SetLineColor(2);
        cout << "Maxmimum: " << NA2_graph->GetFunction("parabFunc")->GetMaximumX(1600, 1800) << endl;
*/

        TF1* myPolFunc = new TF1("myPolFunc", "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x + [7]*x*x*x*x*x*x*x", 1200, 2200);

        NA2_graph->Fit("myPolFunc", "Q");
        NA2_graph->Fit("myPolFunc", "Q");
        NA2_graph->GetFunction("myPolFunc")->SetLineColor(2);

      double maximum = NA2_graph->GetFunction("myPolFunc")->GetMaximumX(1600, 1800);
        cout << "Maxmimum: " << maximum << endl;

      canvasName = (fitTitle + "_NAsq_Vs_" + plotName);

      NA2_graph->GetYaxis()->SetTitle("NA^{2}");
      NA2_graph->SetName(canvasName.c_str());

      NA2_graph->GetYaxis()->SetTitleOffset(2.0);
      NA2_graph->GetXaxis()->SetTitle(xAxisTitle.c_str());

      NA2_graph->Write();

      auto plotCanvas = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,500,400);
      NA2_graph->Draw("AP");

      TPaveText* textBox = new TPaveText(0.75,0.7,0.975,0.8, "NDC");
      textBox->SetBorderSize(1);
      textBox->SetFillStyle(1001);
      textBox->AddText(Form("Maximum: %1.0f", maximum));
      textBox->Draw("SAME");

      plotCanvas->SaveAs(("Images/" + canvasName + datasetTagForPlots + ".png").c_str());

    }

  }

}

/////////////////////////////////////////////////////////////////////////////////////

int energyThresholdScan(std::string filePath)
{
  // gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  // check if dataset tag exists in file and if so append it to the file name and write it to the file

  string outputFileName = "energyThresholdScan";
  TNamed* tag = applyDatasetTag(inputFile, outputFileName);

  TFile* outputFile = new TFile((outputFileName + ".root").c_str(),"RECREATE");
  if(tag) tag->Write();
  auto topDir = outputFile->mkdir("topDir");

/////////////////////////////////////////////////////////////////////////////////////
  // These only get set for the interactive root session (any generated canvases, etc.), but does not apply to the output root file - that comes from .rootlogon.C
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  // gStyle->SetOptFit(2);
  // gStyle->SetOptFit(1111);
  gStyle->SetOptFit(0);
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


std::vector<float> eTs;

std::vector<TF1*> fiveParamFits;
std::vector<TF1*> TmethodFits;
std::vector<TF1*> ratioFits;
std::vector<TF1*> ratioCBOFits;


      TNtuple *firstTuple = (TNtuple*)inputFile->Get(Form("topDir/FitPasses/FitPass0/FitConditions0"));
      float tP;
      firstTuple->SetBranchAddress("totalPasses", &tP);
      firstTuple->GetEntry(0);

    for (int fitPass = 0; fitPass < tP; ++fitPass)
    {
      TVectorD* histSavedParameters = (TVectorD*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/parameterStore", fitPass));
      float eT = (*histSavedParameters)[1];
      eTs.push_back(eT);

      TF1* fiveParamFitFunction = (TF1*) ((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FiveParam/allTimesAdded", fitPass)))->GetFunction("fiveParamFit")->Clone();
      TF1* TMethodFitFunction = (TF1*) ((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/TMethod/allTimesAdded_TMethod", fitPass)))->GetFunction("TmethodFitFunc")->Clone();
      TF1* threeParRatioFitFunction = (TF1*) ((TGraphErrors*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/ThreeParamRatio/Added_Times_Ratio_Graph", fitPass)))->GetFunction("threeParamRatioFit")->Clone();
      TF1* fullRatioFitFunction = (TF1*) ((TGraphErrors*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass)))->GetFunction("fullRatioFitFunc")->Clone();

      fiveParamFits.push_back(fiveParamFitFunction);
      TmethodFits.push_back(TMethodFitFunction);
      ratioFits.push_back(threeParRatioFitFunction);
      ratioCBOFits.push_back(fullRatioFitFunction);
    }

  auto fiveParamDir = addedDir->mkdir("FiveParameter");
  thresholdScanPlots(fiveParamDir, eTs, fiveParamFits, "5 Parameter Fit");

  auto TmethodDir = addedDir->mkdir("TMethod");
  thresholdScanPlots(TmethodDir, eTs, TmethodFits, "T Method Fit");

  auto threeParRatioDir = addedDir->mkdir("ThreeParamRatio");
  thresholdScanPlots(threeParRatioDir, eTs, ratioFits, "3 Parameter Ratio Fit");

  auto fullRatioDir = addedDir->mkdir("FullRatio");
  thresholdScanPlots(fullRatioDir, eTs, ratioCBOFits, "Full Ratio Fit");

/////////////////////////////////////////////////////////////////////////////////////


      delete outputFile;


  return 1;
}
