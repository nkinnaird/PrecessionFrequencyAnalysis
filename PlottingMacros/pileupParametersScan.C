// 3-31-20: Macro for plots for scans over various pileup paramters. The desired boolean value down below should be set to true whether it's the shadow dead time (SDT), shadow gap time (SGT), pileup multiplier (PM), pileup time shift (PTS), or pileup energy scale (PES).

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

bool includeErrorBars = false;
bool drawKawallBands = true;

bool saveImagesDirectly_SDT = false;
bool saveImagesDirectly_SGT = false;
bool saveImagesDirectly_PM  = false;
bool saveImagesDirectly_PTS = false;
bool saveImagesDirectly_PES = false;
bool saveImagesDirectly_ADT = true;

double histStackLegendY = 0;

void stackHists(std::vector<TH1F*> hists, string axisString, std::pair<double, double> xRange, std::pair<double, double> yRange, string plotName, vector<float> legendValues)
{

  THStack* histStack = new THStack("histStack","Temp Title"); // use a THStack because histogram attributes will be preserved when opening the canvas in a TBrowser
  histStack->SetTitle(hists.at(0)->GetName()); // or could use GetTitle

  // auto legend = new TLegend(0.65,histStackLegendY,.92,histStackLegendY+0.4);
  auto legend = new TLegend(0.65,histStackLegendY,.92,histStackLegendY+0.75);

  string legendLabel = plotName.substr(0,plotName.find("_"));

  for (uint i = 0; i < hists.size(); ++i)
  {
    hists.at(i)->SetLineColor(hists.size()-i);
    if(hists.size()-i == 5) hists.at(i)->SetLineColor(94); // replace yellow with orange
    nsTOus(hists.at(i), "doesntmatter");
    histStack->Add(hists.at(i));

    legend->AddEntry(hists.at(i), Form(string(legendLabel + " = %4.3g").c_str(), legendValues.at(i)), "l");
  }


  string canvasName = "stacked_";
  canvasName.append(hists.at(0)->GetName());

  auto stackedCanvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 200, 10, 500, 400);
  histStack->Draw("nostack,hist");

  histStack->GetXaxis()->SetTitle(axisString.c_str());
  histStack->GetYaxis()->SetTitle("Counts");

  histStack->GetXaxis()->SetRangeUser(xRange.first, xRange.second);
  histStack->SetMinimum(yRange.first); //SetRangeUser for Y doesn't want to work so do this instead
  histStack->SetMaximum(yRange.second);


  legend->SetBorderSize(1);
  legend->Draw("SAME");

  stackedCanvas->Write();

  stackedCanvas->SaveAs(("Images/" + plotName).c_str());

  // delete stackedCanvas;
}


void generalScan(TDirectory* inDir, std::tuple<string, string, string> tupleID, std::vector<float> xVector, std::vector<TF1*> fits, bool saveImages)
{
  inDir->cd();

  bool performFits = ((xVector.back() - xVector.front() != 0)) ? true : false;

    string fitName = std::get<0>(tupleID);
    string scanType = std::get<1>(tupleID);
    string titleStr = std::get<2>(tupleID);

    string scanName = removeSpaces(scanType);

    string axisPart = "";
    if(scanType.find("Time") != string::npos) axisPart = " [ns]";

    
/////////////////////////////////////////////////////////////////////////////////////

    TGraphErrors* Chi2_Vs_Val = new TGraphErrors();
      Chi2_Vs_Val->SetName((fitName + "_Chi2_Vs_" + scanName).c_str());
      Chi2_Vs_Val->SetTitle((titleStr + " #chi^{2} Vs " + scanType).c_str());
      Chi2_Vs_Val->GetXaxis()->SetTitle((scanType + axisPart).c_str());
      Chi2_Vs_Val->GetYaxis()->SetTitle("#chi^{2}");
      Chi2_Vs_Val->GetYaxis()->SetTitleOffset(2.4);

    TGraphErrors* Chi2NDF_Vs_Val = new TGraphErrors();
      Chi2NDF_Vs_Val->SetName((fitName + "_Chi2NDF_Vs_" + scanName).c_str());
      Chi2NDF_Vs_Val->SetTitle((titleStr + " #chi^{2}/NDF Vs " + scanType).c_str());
      Chi2NDF_Vs_Val->GetXaxis()->SetTitle((scanType + axisPart).c_str());
      Chi2NDF_Vs_Val->GetYaxis()->SetTitle("#chi^{2}/NDF");
      Chi2NDF_Vs_Val->GetYaxis()->SetTitleOffset(2.4);

    TGraphErrors* Pval_Vs_Val = new TGraphErrors();
      Pval_Vs_Val->SetName((fitName + "_Pval_Vs_" + scanName).c_str());
      Pval_Vs_Val->SetTitle((titleStr + " P Value Vs " + scanType).c_str());
      Pval_Vs_Val->GetXaxis()->SetTitle((scanType + axisPart).c_str());
      Pval_Vs_Val->GetYaxis()->SetTitle("P value");
      Pval_Vs_Val->GetYaxis()->SetTitleOffset(2.4);


    int pointNo = 0;

    for (float xVal : xVector)
    {
      double chi2 = fits.at(pointNo)->GetChisquare();
      double chi2ndf = chi2/fits.at(pointNo)->GetNDF();
      double pVal = fits.at(pointNo)->GetProb();

      Chi2_Vs_Val->SetPoint(pointNo, xVal, chi2);
      Chi2NDF_Vs_Val->SetPoint(pointNo, xVal, chi2ndf);
      Pval_Vs_Val->SetPoint(pointNo, xVal, pVal);

      Chi2_Vs_Val->SetPointError(pointNo, 0, sqrt(2*fits.at(pointNo)->GetNDF()));

      pointNo++;
    } // end loop over point numbers


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
      auto graphCanvas = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,500,400);
      if(includeErrorBars) Chi2_Vs_Val->Draw("AP");
      else{
        adjustGraphRanges(Chi2_Vs_Val);
        Chi2_Vs_Val->Draw("APX");
      } 

      graphCanvas->Update();
      TPaveStats *statsBox = (TPaveStats*)graphCanvas->GetPrimitive("stats");
      statsBox->SetBorderSize(1);

      statsBox->SetX1NDC(0.425);
      statsBox->SetX2NDC(0.725);

      // statsBox->SetX2(graphCanvas->GetUxmax()); // can reposition and by extension resize the stats box with these commands - not the best though, there must be a better way
      // statsBox->SetY2(graphCanvas->GetUymax());

      statsBox->Draw("SAME");
      graphCanvas->SaveAs(("Images/" + canvasName + datasetTagForPlots + ".png").c_str());
    }

/////////////////////////////////////////////////////////////////////////////////////

  int numParams = fits.at(0)->GetNpar();

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

    TGraphErrors* param_Vs_Val = new TGraphErrors();
    param_Vs_Val->SetName((fitName + "_" + nameString + "_Vs_" + scanName).c_str());
    param_Vs_Val->SetTitle((titleStr + " " + paramString + " Vs " + scanType).c_str());

    TGraph* KawallBandOne = new TGraph();
      KawallBandOne->SetName("KawallBandOne");
    TGraph* KawallBandTwo = new TGraph();
      KawallBandTwo->SetName("KawallBandTwo");

    double param_at_start = fits.at(0)->GetParameter(parNum);
    double error_at_start = fits.at(0)->GetParError(parNum);
      
    if(phaseParam) normalizePhase(param_at_start);

      int kawallBandPointNo = 0;

      for (uint fitNum = 0; fitNum < fits.size(); ++fitNum)
      {
        double parameterValue = fits.at(fitNum)->GetParameter(parNum);
        double parameterError = fits.at(fitNum)->GetParError(parNum);
        
        if(phaseParam) normalizePhase(parameterValue);

        param_Vs_Val->SetPoint(fitNum, xVector.at(fitNum), parameterValue);
        param_Vs_Val->SetPointError(fitNum, 0, parameterError);


          double sigmaDiff = sqrt(pow(parameterError, 2) - pow(error_at_start, 2));

          if(!isnan(sigmaDiff))
          {
            KawallBandOne->SetPoint(kawallBandPointNo, xVector.at(fitNum), param_at_start + sigmaDiff);
            KawallBandTwo->SetPoint(kawallBandPointNo, xVector.at(fitNum), param_at_start - sigmaDiff);
            kawallBandPointNo++;
          }
      }

    param_Vs_Val->GetXaxis()->SetTitle((scanType + axisPart).c_str());
    param_Vs_Val->GetYaxis()->SetTitle((paramString + strAddon).c_str());
    param_Vs_Val->GetYaxis()->SetTitleOffset(2.4);

    if(paramString.compare("R") == 0 && performFits){
      param_Vs_Val->Fit("pol1", "Q"); // fit for the systematic effect on R
      param_Vs_Val->GetFunction("pol1")->SetLineColor(2);
    }

    param_Vs_Val->Write();

    if(saveImages && paramString.compare("R") == 0){
      string canvasName = (fitName + "_" + nameString + "_Vs_" + scanName + "_Canv");
      auto graphCanvas = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,500,400);
      if(includeErrorBars) param_Vs_Val->Draw("AP");
      else{
        adjustGraphRanges(param_Vs_Val);
        param_Vs_Val->Draw("APX");
      } 

      if(drawKawallBands){
        KawallBandOne->Draw("SAME");
        KawallBandTwo->Draw("SAME");
      }
      
      graphCanvas->Update();
      TPaveStats *statsBox = (TPaveStats*)graphCanvas->GetPrimitive("stats");
      statsBox->SetBorderSize(1);

      // statsBox->SetX2(graphCanvas->GetUxmax()); // can reposition and by extension resize the stats box with these commands - not the best though, there must be a better way
      // statsBox->SetY2(graphCanvas->GetUymax());

      statsBox->Draw("SAME");
      graphCanvas->SaveAs(("Images/" + canvasName + datasetTagForPlots + ".png").c_str());
    }

  }

}


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


int pileupParametersScan(std::string filePath)
{
  // gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  if((saveImagesDirectly_SDT? 1:0) + (saveImagesDirectly_SGT? 1:0) + (saveImagesDirectly_PM? 1:0) + (saveImagesDirectly_PTS? 1:0) + (saveImagesDirectly_PES? 1:0) > 1)
  {
    printf("Scan booleans set incorrectly\n");
    return 0;
  }

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  // check if dataset tag exists in file and if so append it to the file name and write it to the file

  string outputFileName = "pileupParametersScan";
  TNamed* tag = applyDatasetTag(inputFile, outputFileName);

  TFile* outputFile = new TFile((outputFileName + ".root").c_str(),"RECREATE");
  if(tag) tag->Write();
  auto topDir = outputFile->mkdir("topDir");


/////////////////////////////////////////////////////////////////////////////////////
  // These only get set for the interactive root session (any generated canvases, etc.), but does not apply to the output root file - that comes from .rootlogon.C
  gStyle->SetOptStat(000000);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(2);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1);
  gStyle->SetLineColor(1);

  gStyle->SetPadRightMargin(.05);
  gStyle->SetPadLeftMargin(.2);

  // gStyle->SetOptStat("m");

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  auto addedDir = topDir->mkdir("Added");

/////////////////////////////////////////////////////////////////////////////////////

// values scanned over

// std::vector<float> ADTs; // can't iterate over this yet
std::vector<float> SDTs;
std::vector<float> SGTs;
std::vector<float> PSFs;
std::vector<float> pileupTimeShifts;
std::vector<float> pileupEnergyScales;
std::vector<float> ADTs;

std::vector<TF1*> fullRatioFits;
std::vector<TF1*> TMethodFits;

/////////////////////////////////////////////////////////////////////////////////////

  std::vector<TH1F*> pileupTimesAboveThreshold;
  std::vector<TH1F*> energiesPileupSubtracted;

/////////////////////////////////////////////////////////////////////////////////////

      TNtuple *firstTuple = (TNtuple*)inputFile->Get(Form("topDir/FitPasses/FitPass0/FitConditions0"));
      float tP;
      firstTuple->SetBranchAddress("totalPasses", &tP);
      firstTuple->GetEntry(0);

/////////////////////////////////////////////////////////////////////////////////////

    // loop through passes and store relevant objects

    for (int fitPass = 0; fitPass < tP; ++fitPass)
    {
      TVectorD* histSavedParameters = (TVectorD*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/parameterStore", fitPass));
      SDTs.push_back(float((*histSavedParameters)[2]));
      SGTs.push_back(float((*histSavedParameters)[3]));
      pileupTimeShifts.push_back(float((*histSavedParameters)[6]));
      pileupEnergyScales.push_back(float((*histSavedParameters)[7]));
      ADTs.push_back(float((*histSavedParameters)[4]));

      float pS;

      TNtuple *ntuple = (TNtuple*)inputFile->Get(Form("topDir/FitPasses/FitPass%d/FitConditions%d", fitPass, fitPass));
      ntuple->SetBranchAddress("pileupScaleFactor", &pS);
      ntuple->GetEntry(0);
      PSFs.push_back(pS);

      TF1* TmethodFitFunction = (TF1*) ((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/TMethod/allTimesAdded_TMethod", fitPass)))->GetFunction("TmethodFitFunc")->Clone();
      TMethodFits.push_back(TmethodFitFunction);


      TGraphErrors* ratioGraph = (TGraphErrors*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass));
      TF1* fullRatioFitFunction = 0;
      if(ratioGraph != 0){ // check to make sure ratio fits exist before trying to make plots with them
        fullRatioFitFunction = (TF1*) ratioGraph->GetFunction("fullRatioFitFunc")->Clone();
        fullRatioFits.push_back(fullRatioFitFunction);
      }
      
/////////////////////////////////////////////////////////////////////////////////////

      pileupTimesAboveThreshold.push_back((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/Pileup/added_pileupTimes_shadow_threshold", fitPass))->Clone());
      energiesPileupSubtracted.push_back((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/Input/Energy_Pileup_Subtracted", fitPass))->Clone());
    }

/////////////////////////////////////////////////////////////////////////////////////

    auto added_TMethod_dir = addedDir->mkdir("TMethod");

    auto added_TMethod_SDT_dir = added_TMethod_dir->mkdir("ShadowDeadTime");
    std::tuple<string, string, string> STDinfo_TMethod("TMethod", "Shadow Dead Time", "T Method Fit");
    generalScan(added_TMethod_SDT_dir, STDinfo_TMethod, SDTs, TMethodFits, saveImagesDirectly_SDT);

    auto added_TMethod_SGT_dir = added_TMethod_dir->mkdir("ShadowGapTime");
    std::tuple<string, string, string> SGTinfo_TMethod("TMethod", "Shadow Gap Time", "T Method Fit");
    generalScan(added_TMethod_SGT_dir, SGTinfo_TMethod, SGTs, TMethodFits, saveImagesDirectly_SGT);

    auto added_TMethod_PSF_dir = added_TMethod_dir->mkdir("PileupMultiplier");
    std::tuple<string, string, string> PSFinfo_TMethod("TMethod", "Pileup Multiplier", "T Method Fit");
    generalScan(added_TMethod_PSF_dir, PSFinfo_TMethod, PSFs, TMethodFits, saveImagesDirectly_PM);

    auto added_TMethod_PTS_dir = added_TMethod_dir->mkdir("PileupTimeShift");
    std::tuple<string, string, string> pTimeShiftInfo_TMethod("TMethod", "Pileup Time Shift", "T Method Fit");
    generalScan(added_TMethod_PTS_dir, pTimeShiftInfo_TMethod, pileupTimeShifts, TMethodFits, saveImagesDirectly_PTS);

    auto added_TMethod_PES_dir = added_TMethod_dir->mkdir("PileupEnergyScaling");
    std::tuple<string, string, string> pEScalingInfo_TMethod("TMethod", "Pileup Energy Scaling", "T Method Fit");
    generalScan(added_TMethod_PES_dir, pEScalingInfo_TMethod, pileupEnergyScales, TMethodFits, saveImagesDirectly_PES);

    auto added_TMethod_ADT_dir = added_TMethod_dir->mkdir("ArtificialDeadTime");
    std::tuple<string, string, string> ADTinfo_TMethod("TMethod", "Artificial Dead Time", "T Method Fit");
    generalScan(added_TMethod_ADT_dir, ADTinfo_TMethod, ADTs, TMethodFits, saveImagesDirectly_ADT);

/////////////////////////////////////////////////////////////////////////////////////

if(fullRatioFits.size() > 0) {

    auto added_FullRatio_dir = addedDir->mkdir("FullRatio");

    auto added_FullRatio_SDT_dir = added_FullRatio_dir->mkdir("ShadowDeadTime");
    std::tuple<string, string, string> STDinfo_ratio("FullRatio", "Shadow Dead Time", "Full Ratio Fit");
    generalScan(added_FullRatio_SDT_dir, STDinfo_ratio, SDTs, fullRatioFits, saveImagesDirectly_SDT);

    auto added_FullRatio_SGT_dir = added_FullRatio_dir->mkdir("ShadowGapTime");
    std::tuple<string, string, string> SGTinfo_ratio("FullRatio", "Shadow Gap Time", "Full Ratio Fit");
    generalScan(added_FullRatio_SGT_dir, SGTinfo_ratio, SGTs, fullRatioFits, saveImagesDirectly_SGT);

    auto added_FullRatio_PSF_dir = added_FullRatio_dir->mkdir("PileupMultiplier");
    std::tuple<string, string, string> PSFinfo_ratio("FullRatio", "Pileup Multiplier", "Full Ratio Fit");
    generalScan(added_FullRatio_PSF_dir, PSFinfo_ratio, PSFs, fullRatioFits, saveImagesDirectly_PM);

    auto added_FullRatio_PTS_dir = added_FullRatio_dir->mkdir("PileupTimeShift");
    std::tuple<string, string, string> pTimeShiftInfo_ratio("FullRatio", "Pileup Time Shift", "Full Ratio Fit");
    generalScan(added_FullRatio_PTS_dir, pTimeShiftInfo_ratio, pileupTimeShifts, fullRatioFits, saveImagesDirectly_PTS);

    auto added_FullRatio_PES_dir = added_FullRatio_dir->mkdir("PileupEnergyScaling");
    std::tuple<string, string, string> pEScalingInfo_ratio("FullRatio", "Pileup Energy Scaling", "Full Ratio Fit");
    generalScan(added_FullRatio_PES_dir, pEScalingInfo_ratio, pileupEnergyScales, fullRatioFits, saveImagesDirectly_PES);

    auto added_FullRatio_ADT_dir = added_FullRatio_dir->mkdir("ArtificialDeadTime");
    std::tuple<string, string, string> ADTinfo_ratio("FullRatio", "Artificial Dead Time", "Full Ratio Fit");
    generalScan(added_FullRatio_ADT_dir, ADTinfo_ratio, ADTs, fullRatioFits, saveImagesDirectly_ADT);

}
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

    auto PSFsVsDir = addedDir->mkdir("PSFsVs"); // should try and save psf errors for these plots
    PSFsVsDir->cd();

    TGraph* PSFs_vs_Iter = new TGraph();
      PSFs_vs_Iter->SetName("PSFs_vs_Iter");
      PSFs_vs_Iter->SetTitle("PSFs_vs_Iter");
        PSFs_vs_Iter->GetXaxis()->SetTitle("Iter");
        PSFs_vs_Iter->GetYaxis()->SetTitle("PSF");

    TGraph* PSFs_vs_SDT = new TGraph();
      PSFs_vs_SDT->SetName("PSFs_vs_SDT");
      PSFs_vs_SDT->SetTitle("PSFs_vs_SDT");
        PSFs_vs_SDT->GetXaxis()->SetTitle("SDT [ns]");
        PSFs_vs_SDT->GetYaxis()->SetTitle("PSF");

    TGraph* PSFs_vs_SGT = new TGraph();
      PSFs_vs_SGT->SetName("PSFs_vs_SGT");
      PSFs_vs_SGT->SetTitle("PSFs_vs_SGT");
        PSFs_vs_SGT->GetXaxis()->SetTitle("SGT [ns]");
        PSFs_vs_SGT->GetYaxis()->SetTitle("PSF");

    TGraph* PSFs_vs_ADT = new TGraph();
      PSFs_vs_ADT->SetName("PSFs_vs_ADT");
      PSFs_vs_ADT->SetTitle("PSFs_vs_ADT");
        PSFs_vs_ADT->GetXaxis()->SetTitle("ADT [ns]");
        PSFs_vs_ADT->GetYaxis()->SetTitle("PSF");

    for (uint i = 0; i < PSFs.size(); ++i)
    {
      PSFs_vs_Iter->SetPoint(i, i, PSFs.at(i));
      PSFs_vs_SDT->SetPoint(i, SDTs.at(i), PSFs.at(i));
      PSFs_vs_SGT->SetPoint(i, SGTs.at(i), PSFs.at(i));
      PSFs_vs_ADT->SetPoint(i, ADTs.at(i), PSFs.at(i));
    }

    PSFs_vs_Iter->Write("PSFs_vs_Iter");
    PSFs_vs_SDT->Write("PSFs_vs_SDT");
    PSFs_vs_SGT->Write("PSFs_vs_SGT");
    PSFs_vs_ADT->Write("PSFs_vs_ADT");

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  // make pileup energy and time comparison plots

  vector<float> valuesForLegend;

  string plotName = "noname";
  if(saveImagesDirectly_SDT){
    plotName = "SDT";
    valuesForLegend = SDTs;
  } 
  else if(saveImagesDirectly_SGT){
    plotName = "SGT";
    valuesForLegend = SGTs;
  } 
  else if(saveImagesDirectly_PM){
    plotName = "PM";
    valuesForLegend = PSFs;
  } 
  else if(saveImagesDirectly_PTS){
    plotName = "PTS";
    valuesForLegend = pileupTimeShifts;
  } 
  else if(saveImagesDirectly_PES){
    plotName = "PES";
    valuesForLegend = pileupEnergyScales;
  } 
  else if(saveImagesDirectly_ADT){
    plotName = "ADT";
    valuesForLegend = ADTs;
  }   

/////////////////////////////////////////////////////////////////////////////////////

    gStyle->SetPadLeftMargin(.15);

    auto PileupIterationComparison = addedDir->mkdir("PileupIterationComparison");
    PileupIterationComparison->cd();

    histStackLegendY = 0.2;
    stackHists(pileupTimesAboveThreshold, "Time (#mus)", make_pair(20, 60), make_pair(0, 45000), plotName+"_PileupTimeComparison" + datasetTagForPlots + ".png", valuesForLegend); // 20, 45, 80 for 60h, 9d, Endgame datasets
    histStackLegendY = 0.2;
    stackHists(energiesPileupSubtracted, "Energy (GeV)", make_pair(3, 7), make_pair(-2000, 2000), plotName+"_CorrEnergyComparison" + datasetTagForPlots + ".png", valuesForLegend);

/////////////////////////////////////////////////////////////////////////////////////

    // make energy ratio plot

    TH1F* energyRatioHist = (TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/Pileup/AutoScaleFactor/EnergyRatio");
    nsTOus(energyRatioHist, "Energy (GeV)");
    energyRatioHist->GetYaxis()->SetTitle("Cluster Energies / Pileup Energies");
    energyRatioHist->GetXaxis()->SetRangeUser(2, 9);

    double xMin, xMax;
    TF1* scaleFunc = energyRatioHist->GetFunction("scaleFactor");
    scaleFunc->GetRange(xMin, xMax);
    scaleFunc->SetRange(xMin/1000, xMax/1000);

    auto canv_energyRatio = new TCanvas("canv_energyRatio","canv_energyRatio",200,10,500,400);
    energyRatioHist->Draw("");

    canv_energyRatio->SaveAs(("Images/EnergyRatio" + datasetTagForPlots + ".png").c_str());

/////////////////////////////////////////////////////////////////////////////////////

      delete outputFile;


  return 1;
}
