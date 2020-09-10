// 3-31-20: Macro for plots for scans over fit start or end time.

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

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;

bool fitEndScan = false;
bool saveImagesDirectly = true; // only applies to Ratio CBO fit since those are the plots I care most about

void fitStartScanPlots(TDirectory* inDir, std::vector<pair<float, float> > fitRanges, std::vector<TF1*> fits, string titleStr, bool saveImages, TH1F* fittedHistogram)
{
  inDir->cd();

  string scanType = "Fit Start Time";
  string nameChar = "S";
  if(fitEndScan){
    scanType = "Fit End Time";
    nameChar = "E";
  }

  string xAxisTitle = scanType + " [#mus]";
  string fitTitle = inDir->GetName();

  int numParams = fits.at(0)->GetNpar();

/////////////////////////////////////////////////////////////////////////////////////

  TGraphErrors* Chi2NDF_Vs_FS = new TGraphErrors();
  Chi2NDF_Vs_FS->SetName((fitTitle + "_Chi2NDF_Vs_F" + nameChar).c_str());
  Chi2NDF_Vs_FS->SetTitle((titleStr + " #chi^{2}/NDF Vs " + scanType).c_str());
  Chi2NDF_Vs_FS->GetXaxis()->SetTitle(xAxisTitle.c_str());
  Chi2NDF_Vs_FS->GetYaxis()->SetTitle("#chi^{2}/NDF");

  TGraph* Chi2BandOne = new TGraph();
    Chi2BandOne->SetName("KawallBandOne");
  TGraph* Chi2BandTwo = new TGraph();
    Chi2BandTwo->SetName("KawallBandTwo");

  TGraph* Pval_Vs_FS = new TGraph();
  Pval_Vs_FS->SetName((fitTitle + "_Pval_Vs_F" + nameChar).c_str());
  Pval_Vs_FS->SetTitle((titleStr + " P Value Vs " + scanType).c_str());
  Pval_Vs_FS->GetXaxis()->SetTitle(xAxisTitle.c_str());
  Pval_Vs_FS->GetYaxis()->SetTitle("P value");

  for (uint fitNum = 0; fitNum < fits.size(); ++fitNum)
  {
    float fitTimePoint = fitEndScan ? fitRanges.at(fitNum).second : fitRanges.at(fitNum).first;

    Chi2NDF_Vs_FS->SetPoint(fitNum, fitTimePoint, fits.at(fitNum)->GetChisquare()/fits.at(fitNum)->GetNDF());
    Chi2NDF_Vs_FS->SetPointError(fitNum, 0, sqrt(2./fits.at(fitNum)->GetNDF()));

    double chi2ndf_at_start = fits.at(0)->GetChisquare()/fits.at(0)->GetNDF();
    double sigmaDiff = sqrt(2. * (1./fits.at(fitNum)->GetNDF() - 1./fits.at(0)->GetNDF()));
      Chi2BandOne->SetPoint(fitNum, fitTimePoint, chi2ndf_at_start + sigmaDiff);
      Chi2BandTwo->SetPoint(fitNum, fitTimePoint, chi2ndf_at_start - sigmaDiff);

    Pval_Vs_FS->SetPoint(fitNum, fitTimePoint, fits.at(fitNum)->GetProb());
  }

  nsTOus(Chi2NDF_Vs_FS, xAxisTitle);
  nsTOus((TGraphErrors*)Chi2BandOne, xAxisTitle);
  nsTOus((TGraphErrors*)Chi2BandTwo, xAxisTitle);

  string canvName = (fitTitle + "_Chi2NDF_Vs_F" + nameChar + "_canv");
  auto chi2ndfcanv = new TCanvas(canvName.c_str(), canvName.c_str(), 200, 10, 600, 400);
  Chi2NDF_Vs_FS->Draw("AP");
  Chi2BandOne->Draw("SAME");
  Chi2BandTwo->Draw("SAME");

  chi2ndfcanv->Write();
  if(saveImages) chi2ndfcanv->SaveAs(("Images/" + canvName + datasetTagForPlots + ".png").c_str());

  string Pval_canvName = (fitTitle + "_Pval_Vs_F" + nameChar + "_canv");
  auto pvalcanv = new TCanvas(Pval_canvName.c_str(), Pval_canvName.c_str(), 200, 10, 600, 400);

  nsTOus(Pval_Vs_FS, xAxisTitle);
  Pval_Vs_FS->Draw("APL");

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

    TGraphErrors* param_Vs_FS = new TGraphErrors();
    param_Vs_FS->SetName((fitTitle + "_" + nameString + "_Vs_F" + nameChar).c_str());
    param_Vs_FS->SetTitle((titleStr + " " + paramString + " Vs " + scanType).c_str());

    TGraph* KawallBandOne = new TGraph();
      KawallBandOne->SetName("KawallBandOne");
    TGraph* KawallBandTwo = new TGraph();
      KawallBandTwo->SetName("KawallBandTwo");

    double param_at_start = fits.at(0)->GetParameter(parNum);
    double error_at_start = fits.at(0)->GetParError(parNum);

/////////////////////////////////////////////////////////////////////////////////////
    // tests with analytically calculating the Kawall band width - extra factor that I'm not sure how to deal with though

    double histBinWidth = fittedHistogram->GetBinWidth(1);
    double fitCounts_at_start = fittedHistogram->Integral(fitRanges.at(0).first/histBinWidth, fitRanges.at(0).second/histBinWidth);

    double tSA = pow((fitRanges.at(0).first / defaultLifetime + 1),2) + 1;

    TGraph* KawallBandOne_temp = new TGraph();
    TGraph* KawallBandTwo_temp = new TGraph();

    KawallBandOne_temp->SetLineColor(2);
    KawallBandOne_temp->SetMarkerColor(2);
    KawallBandTwo_temp->SetLineColor(2);
    KawallBandTwo_temp->SetMarkerColor(2);
/////////////////////////////////////////////////////////////////////////////////////



      if(phaseParam) normalizePhase(param_at_start);
      normalizePhase(g2Phase_at_start);

      int kawallBandPointNo = 0;
      for (uint fitNum = 0; fitNum < fits.size(); ++fitNum)
      {
        float fitStart = fitRanges.at(fitNum).first;
        float fitEnd = fitRanges.at(fitNum).second;
        float fitTimePoint = fitEndScan ? fitEnd : fitStart;

        double fitCounts = fittedHistogram->Integral(fitStart/histBinWidth, fitEnd/histBinWidth);


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

        // double analyzingPower = (2. * (asymmetry_at_start/asymmetry) * cos(g2Phase_at_start-g2Phase) - 1); // technically sigmaDiff is dependent on this value, but it's value is always very close to 1 so it doesn't really matter

        param_Vs_FS->SetPoint(fitNum, fitTimePoint, parameterValue);
        param_Vs_FS->SetPointError(fitNum, 0, parameterError);
          // double sigmaDiff = sqrt(pow(parameterError, 2) - pow(error_at_start, 2)*analyzingPower);
          double sigmaDiff = sqrt(pow(parameterError, 2) - pow(error_at_start, 2));

          if(!isnan(sigmaDiff))
          {
            KawallBandOne->SetPoint(kawallBandPointNo, fitTimePoint, param_at_start + sigmaDiff);
            KawallBandTwo->SetPoint(kawallBandPointNo, fitTimePoint, param_at_start - sigmaDiff);
            kawallBandPointNo++;
          }
          else cout << "nan point in bands for parameter: " << paramString << " at point " << fitNum << endl;


          // double sigmaDiffApprox = error_at_start * sqrt(fitCounts_at_start/fitCounts - 1);
          double tSB = pow((fitStart / defaultLifetime + 1),2) + 1;
          double sigmaDiffApprox = error_at_start * sqrt((fitCounts_at_start * tSB)/(fitCounts * tSA) - 1);


            KawallBandOne_temp->SetPoint(fitNum, fitTimePoint, param_at_start + sigmaDiffApprox);
            KawallBandTwo_temp->SetPoint(fitNum, fitTimePoint, param_at_start - sigmaDiffApprox);

      }

    nsTOus(param_Vs_FS, xAxisTitle);
    nsTOus(KawallBandOne, xAxisTitle);
    nsTOus(KawallBandTwo, xAxisTitle);

    nsTOus(KawallBandOne_temp, xAxisTitle);
    nsTOus(KawallBandTwo_temp, xAxisTitle);

    string canvasName = (fitTitle + "_" + nameString + "_F" + nameChar + "_Canv");
    auto fitStart_canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 200, 10, 600, 400);
    param_Vs_FS->GetXaxis()->SetTitle(xAxisTitle.c_str());
    param_Vs_FS->GetYaxis()->SetTitle((paramString + strAddon).c_str());
    // param_Vs_FS->GetYaxis()->SetTitleOffset(1.8);
    param_Vs_FS->Draw("AP");
    KawallBandOne->Draw("SAME");
    KawallBandTwo->Draw("SAME");

    // KawallBandOne_temp->Draw("SAME");
    // KawallBandTwo_temp->Draw("SAME");

    fitStart_canvas->Write();
    if(saveImages) fitStart_canvas->SaveAs(("Images/" + canvasName + datasetTagForPlots + ".png").c_str());
    // delete fitStart_canvas;
    
  }

}

/////////////////////////////////////////////////////////////////////////////////////

int fitStartScan(std::string filePath)
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

  string outputFileName = "fitStartScan";
  TNamed* tag = applyDatasetTag(inputFile, outputFileName);

  TFile* outputFile = new TFile((outputFileName + ".root").c_str(),"RECREATE");
  if(tag) tag->Write();
  auto topDir = outputFile->mkdir("topDir");

/////////////////////////////////////////////////////////////////////////////////////
  // These only get set for the interactive root session (any generated canvases, etc.), but does not apply to the output root file - that comes from .rootlogon.C
  // gStyle->SetOptStat(000000);
  // gStyle->SetOptTitle(0);
  // gStyle->SetOptFit(0000);
  // gStyle->SetMarkerStyle(20);
  // gStyle->SetMarkerColor(1);
  // gStyle->SetMarkerSize(1);
  // gStyle->SetLineColor(1);

  // gStyle->SetPadRightMargin(.05);

  SetGm2Style();

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

// perform fit start scan for added data
    for (int fitPass = 0; fitPass < tP; ++fitPass)
    {
      float fitStart;
      float fitEnd;

      TNtuple *ntuple = (TNtuple*)inputFile->Get(Form("topDir/FitPasses/FitPass%d/FitConditions%d", fitPass, fitPass));
      ntuple->SetBranchAddress("fitStartTime", &fitStart);
      ntuple->SetBranchAddress("fitEndTime", &fitEnd);
      ntuple->GetEntry(0);

      // TF1* fiveParamFitFunction = (TF1*) ((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/Hist/allTimesAdded", fitPass)))->GetFunction("fiveParamFit")->Clone();
      TF1* TMethodFitFunction = (TF1*) ((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/TMethod/allTimesAdded_TMethod", fitPass)))->GetFunction("TmethodFitFunc")->Clone();
      // TF1* theeParRatioFitFunction = (TF1*) ((TGraphErrors*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/Ratio/Added_Times_Ratio_Graph", fitPass)))->GetFunction("threeParamRatioFit")->Clone();
      TF1* fullRatioFitFunction = (TF1*) ((TGraphErrors*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass)))->GetFunction("fullRatioFitFunc")->Clone();
      // TF1* fullRatioFitFunction = (TF1*) ((TGraphErrors*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/RatioCBO/Added_Times_Ratio_CBO_Graph", fitPass)))->GetFunction("cboRatioFitFunc")->Clone();

      fitRanges.emplace_back(fitStart, fitEnd);
      // fiveParamFits.push_back(fiveParamFitFunction);
      TmethodFits.push_back(TMethodFitFunction);
      // threeParRatioFits.push_back(theeParRatioFitFunction);
      fullRatioFits.push_back(fullRatioFitFunction);
    }


    TH1F* fittedHistogram = (TH1F*) ((TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/Input/Times_E_Threshold_Pileup_Subtracted"))->Clone();


  // auto fiveParamDir = addedDir->mkdir("FiveParameter");
  // fitStartScanPlots(fiveParamDir, fitRanges, fiveParamFits, "5 Parameter Fit", false, fittedHistogram);

  auto TmethodDir = addedDir->mkdir("TMethod");
  fitStartScanPlots(TmethodDir, fitRanges, TmethodFits, "T Method Fit", saveImagesDirectly, fittedHistogram);

  // auto ratioDir = addedDir->mkdir("Ratio");
  // fitStartScanPlots(ratioDir, fitRanges, threeParRatioFits, "3 Parameter Ratio Fit", false);

  auto fullRatioDir = addedDir->mkdir("FullRatio");
  fitStartScanPlots(fullRatioDir, fitRanges, fullRatioFits, "Full Ratio Fit", saveImagesDirectly, fittedHistogram);

/////////////////////////////////////////////////////////////////////////////////////

  delete outputFile;

  return 1;
}
