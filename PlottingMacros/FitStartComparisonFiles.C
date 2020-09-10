// this macro compares fit start scans in two separate files on single seeds, typically with and without the ad hoc correction applied - see fitStartScan.C

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
#include <TLegend.h>

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;

bool fitEndScan = false;
bool saveImagesDirectly = true;

void fitStartScanPlots(TDirectory* inDir, vector<pair<float, float> > fitRanges, vector<TF1*> fits, vector<TF1*> fits_second, string titleStr, bool saveImages)
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
  if(fits_second.at(0)->GetNpar() != numParams){
    cout << "Number of parameters in fits are different - exiting." << endl;
    exit(0);
  }

/////////////////////////////////////////////////////////////////////////////////////

  TGraphErrors* Chi2NDF_Vs_FS = new TGraphErrors();
  Chi2NDF_Vs_FS->SetName((fitTitle + "_Chi2NDF_Vs_F" + nameChar).c_str());
  Chi2NDF_Vs_FS->SetTitle((titleStr + " #chi^{2}/NDF Vs " + scanType).c_str());
  Chi2NDF_Vs_FS->GetXaxis()->SetTitle(xAxisTitle.c_str());
  Chi2NDF_Vs_FS->GetYaxis()->SetTitle("#chi^{2}/NDF");

  TGraphErrors* Chi2NDF_Vs_FS_second = new TGraphErrors();
  // Chi2NDF_Vs_FS_second->SetLineStyle(2);
  Chi2NDF_Vs_FS_second->SetMarkerColor(2);
  Chi2NDF_Vs_FS_second->SetLineColor(2);

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

    Chi2NDF_Vs_FS_second->SetPoint(fitNum, fitTimePoint, fits_second.at(fitNum)->GetChisquare()/fits_second.at(fitNum)->GetNDF());
    Chi2NDF_Vs_FS_second->SetPointError(fitNum, 0, sqrt(2./fits_second.at(fitNum)->GetNDF()));

    double chi2ndf_at_start = fits.at(0)->GetChisquare()/fits.at(0)->GetNDF();
    double sigmaDiff = sqrt(2. * (1./fits.at(fitNum)->GetNDF() - 1./fits.at(0)->GetNDF()));
      Chi2BandOne->SetPoint(fitNum, fitTimePoint, chi2ndf_at_start + sigmaDiff);
      Chi2BandTwo->SetPoint(fitNum, fitTimePoint, chi2ndf_at_start - sigmaDiff);

    Pval_Vs_FS->SetPoint(fitNum, fitTimePoint, fits.at(fitNum)->GetProb());
  }

  nsTOus(Chi2NDF_Vs_FS, xAxisTitle);
  nsTOus(Chi2NDF_Vs_FS_second, xAxisTitle);
  nsTOus((TGraphErrors*)Chi2BandOne, xAxisTitle);
  nsTOus((TGraphErrors*)Chi2BandTwo, xAxisTitle);

  string canvName = (fitTitle + "_Chi2NDF_Vs_F" + nameChar + "_canv");
  auto chi2ndfcanv = new TCanvas(canvName.c_str(), canvName.c_str(), 200, 10, 600, 400);
  Chi2NDF_Vs_FS->Draw("ALX");
  Chi2NDF_Vs_FS_second->Draw("LXSAME");
  Chi2BandOne->Draw("SAME");
  Chi2BandTwo->Draw("SAME");

        auto legend_chi2 = new TLegend(0.7,0.675,0.975,0.875);

        legend_chi2->AddEntry(Chi2NDF_Vs_FS, "Uncorrected","l");
        legend_chi2->AddEntry(Chi2NDF_Vs_FS_second, "Corrected","l");
  
        legend_chi2->SetBorderSize(1);
        legend_chi2->SetFillStyle(1001);
        legend_chi2->Draw();

  chi2ndfcanv->Write();
  if(saveImages) chi2ndfcanv->SaveAs(("Images/" + canvName + "_compare.png").c_str());

  string Pval_canvName = (fitTitle + "_Pval_Vs_F" + nameChar + "_canv");
  auto pvalcanv = new TCanvas(Pval_canvName.c_str(), Pval_canvName.c_str(), 200, 10, 600, 400);

  nsTOus(Pval_Vs_FS, xAxisTitle);
  Pval_Vs_FS->Draw("APL");

  Pval_Vs_FS->Write();
  if(saveImages) pvalcanv->SaveAs(("Images/" + string(Pval_Vs_FS->GetName()) + "_compare.png").c_str());

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

    TGraphErrors* param_Vs_FS_second = new TGraphErrors();
    // param_Vs_FS_second->SetLineStyle(2);
    param_Vs_FS_second->SetMarkerColor(2);
    param_Vs_FS_second->SetLineColor(2);

    double param_at_start_second = fits_second.at(0)->GetParameter(parNum);

    TGraphErrors* param_Vs_FS_diff = new TGraphErrors();
    // param_Vs_FS_diff->SetLineStyle(4);
    param_Vs_FS_diff->SetMarkerColor(4);
    param_Vs_FS_diff->SetLineColor(4);


/////////////////////////////////////////////////////////////////////////////////////



      if(phaseParam) normalizePhase(param_at_start);
      normalizePhase(g2Phase_at_start);

      int kawallBandPointNo = 0;
      for (uint fitNum = 0; fitNum < fits.size(); ++fitNum)
      {
        float fitStart = fitRanges.at(fitNum).first;
        float fitEnd = fitRanges.at(fitNum).second;
        float fitTimePoint = fitEndScan ? fitEnd : fitStart;


        double parameterValue = fits.at(fitNum)->GetParameter(parNum);
        double parameterError = fits.at(fitNum)->GetParError(parNum);

        double g2Phase = fits.at(fitNum)->GetParameter(phaseParamNum);
        double asymmetry = fits.at(fitNum)->GetParameter(asymmetryParamNum);

                double parameterValue_second = fits_second.at(fitNum)->GetParameter(parNum);


          if(phaseParam) {
            normalizePhase(parameterValue);
              if(fitNum > 0){
                double previousPhaseParameter = fits.at(fitNum-1)->GetParameter(parNum);
                normalizePhase(previousPhaseParameter);
                if(parameterValue < (previousPhaseParameter - pi)) parameterValue += 2*pi; // to clean up plots with phase parameter
              }
          }
          normalizePhase(g2Phase);

        double analyzingPower = (2. * (asymmetry_at_start/asymmetry) * cos(g2Phase_at_start-g2Phase) - 1);

        param_Vs_FS->SetPoint(fitNum, fitTimePoint, parameterValue - param_at_start);
        param_Vs_FS->SetPointError(fitNum, 0, parameterError);
          double sigmaDiff = sqrt(pow(parameterError, 2) - pow(error_at_start, 2)*analyzingPower);

          if(!isnan(sigmaDiff))
          {
            KawallBandOne->SetPoint(kawallBandPointNo, fitTimePoint, sigmaDiff);
            KawallBandTwo->SetPoint(kawallBandPointNo, fitTimePoint, -sigmaDiff);
            kawallBandPointNo++;
          }
          else cout << "nan point in bands for parameter: " << paramString << " at point " << fitNum << endl;

        param_Vs_FS_second->SetPoint(fitNum, fitTimePoint, parameterValue_second - param_at_start_second);
        // param_Vs_FS_diff->SetPoint(fitNum, fitTimePoint, parameterValue_second - parameterValue);
        param_Vs_FS_diff->SetPoint(fitNum, fitTimePoint, (parameterValue_second - param_at_start_second) - (parameterValue - param_at_start)); // normalize to zero at start


      }

    nsTOus(param_Vs_FS, xAxisTitle);
    nsTOus(KawallBandOne, xAxisTitle);
    nsTOus(KawallBandTwo, xAxisTitle);

    nsTOus(param_Vs_FS_second, xAxisTitle);
    nsTOus(param_Vs_FS_diff, xAxisTitle);


    string canvasName = (fitTitle + "_" + nameString + "_F" + nameChar + "_Canv");
    auto fitStart_canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 200, 10, 600, 400);
    param_Vs_FS->GetXaxis()->SetTitle(xAxisTitle.c_str());
    param_Vs_FS->GetYaxis()->SetTitle(("#Delta" + paramString + strAddon).c_str());
    param_Vs_FS->GetYaxis()->SetTitleOffset(1.8);
    param_Vs_FS->Draw("ALX");
    param_Vs_FS_second->Draw("LSAME");
    param_Vs_FS_diff->Draw("LSAME");
    KawallBandOne->Draw("SAME");
    KawallBandTwo->Draw("SAME");

        auto legend_param = new TLegend(0.7,0.675,0.975,0.875);

        legend_param->AddEntry(param_Vs_FS, "Uncorrected","l");
        legend_param->AddEntry(param_Vs_FS_second, "Corrected","l");
        legend_param->AddEntry(param_Vs_FS_diff, "Difference","l");
  
        legend_param->SetBorderSize(1);
        legend_param->SetFillStyle(1001);
        legend_param->Draw();

    fitStart_canvas->Write();
    if(saveImages) fitStart_canvas->SaveAs(("Images/" + canvasName + "_compare.png").c_str());
    // delete fitStart_canvas;
    
  }

}

/////////////////////////////////////////////////////////////////////////////////////

int FitStartComparisonFiles()
{
  // gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen
  // when setting this to true, for some reason there is an extra space in the axis title with the micro symbol in microseconds in the saved images - no idea why..

  // 60h
  // string firstFile_string = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/60h/SingleIter/FitStartScan-NewRange/output-60h-FS-NewRange.root"; // typically uncorrected file
  // string secondFile_string = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/60h/Gain/IFG/Amplitude-with-AdHoc/FitStartScans/Ordinary/output-60h-FS-withAdHoc.root"; // typically corrected file
  // string secondFile_string = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/60h/Gain/AdHoc/FirstTests/Amplitude/FitStartScan-HighAmp/output-60h-FS-adHoc-highAmp.root"; // typically corrected file

  // Endgame
  string firstFile_string = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/Endgame/SingleIter/FitStartScan-NewRanges/output-Endgame-FS-NewRanges.root"; // typically uncorrected file
  string secondFile_string = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/Endgame/Gain/AdHoc/Amplitude/FitStartScan/output-Endgame-FS-withAdHoc-mostParamsFixed.root"; // typically corrected file

  TFile *firstFile = TFile::Open(firstFile_string.c_str());
  TFile *secondFile = TFile::Open(secondFile_string.c_str());
   if (firstFile == 0 || secondFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  TFile* outputFile = new TFile("fitStartScanComparison.root","RECREATE");
  auto topDir = outputFile->mkdir("topDir");

/////////////////////////////////////////////////////////////////////////////////////
  // These only get set for the interactive root session (any generated canvases, etc.), but does not apply to the output root file - that comes from .rootlogon.C
  gStyle->SetOptStat(000000);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0000);
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
std::vector<TF1*> TmethodFits_firstFile;
std::vector<TF1*> TmethodFits_secondFile;
std::vector<TF1*> fullRatioFits_firstFile;
std::vector<TF1*> fullRatioFits_secondFile;


      TNtuple *firstTuple = (TNtuple*)firstFile->Get(Form("topDir/FitPasses/FitPass0/FitConditions0"));
      float tP;
      firstTuple->SetBranchAddress("totalPasses", &tP);
      firstTuple->GetEntry(0);

      // make sure two files were fitted in the same way
      TNtuple *checkTuple = (TNtuple*)secondFile->Get(Form("topDir/FitPasses/FitPass0/FitConditions0"));
      float tP_check;
      checkTuple->SetBranchAddress("totalPasses", &tP_check);
      checkTuple->GetEntry(0);

      if(tP != tP_check){
        cout << "Total fit passes are different between files: " << tP << " vs " << tP_check << " - exiting." << endl;
        return 0;
      }



    for (int fitPass = 0; fitPass < tP; ++fitPass)
    {
      float fitStart, fitStart_check;
      float fitEnd, fitEnd_check;

      TNtuple *ntuple = (TNtuple*)firstFile->Get(Form("topDir/FitPasses/FitPass%d/FitConditions%d", fitPass, fitPass));
      ntuple->SetBranchAddress("fitStartTime", &fitStart);
      ntuple->SetBranchAddress("fitEndTime", &fitEnd);
      ntuple->GetEntry(0);

      TNtuple *ntuple_check = (TNtuple*)secondFile->Get(Form("topDir/FitPasses/FitPass%d/FitConditions%d", fitPass, fitPass));
      ntuple_check->SetBranchAddress("fitStartTime", &fitStart_check);
      ntuple_check->SetBranchAddress("fitEndTime", &fitEnd_check);
      ntuple_check->GetEntry(0);

      if(fitStart != fitStart_check || fitEnd != fitEnd_check){
        cout << "Fit ranges are different between files: " << fitStart << " vs " << fitStart_check << " and " << fitEnd << " vs " << fitEnd_check << " - exiting." << endl;
        return 0;
      }

      TF1* TMethodFitFunction_first = (TF1*) ((TH1F*) firstFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/TMethod/allTimesAdded_TMethod", fitPass)))->GetFunction("TmethodFitFunc")->Clone();
      TF1* TMethodFitFunction_second = (TF1*) ((TH1F*) secondFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/TMethod/allTimesAdded_TMethod", fitPass)))->GetFunction("TmethodFitFunc")->Clone();

      TF1* fullRatioFitFunction_first = (TF1*) ((TGraphErrors*) firstFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass)))->GetFunction("fullRatioFitFunc")->Clone();
      TF1* fullRatioFitFunction_second = (TF1*) ((TGraphErrors*) secondFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass)))->GetFunction("fullRatioFitFunc")->Clone();

      fitRanges.emplace_back(fitStart, fitEnd);
      TmethodFits_firstFile.push_back(TMethodFitFunction_first);
      TmethodFits_secondFile.push_back(TMethodFitFunction_second);
      fullRatioFits_firstFile.push_back(fullRatioFitFunction_first);
      fullRatioFits_secondFile.push_back(fullRatioFitFunction_second);
    }

  auto TmethodDir = addedDir->mkdir("TMethod");
  fitStartScanPlots(TmethodDir, fitRanges, TmethodFits_firstFile, TmethodFits_secondFile, "T Method Fit", saveImagesDirectly);

  auto fullRatioDir = addedDir->mkdir("FullRatio");
  fitStartScanPlots(fullRatioDir, fitRanges, fullRatioFits_firstFile, fullRatioFits_secondFile, "Full Ratio Fit", saveImagesDirectly);

/////////////////////////////////////////////////////////////////////////////////////

  delete outputFile;

  return 1;
}
