// 3-31-20: Macro for plots for scans over bunch number.

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
#include <TLatex.h>

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;

bool saveImagesDirectly = true;

double lineFunc(double* x, double* p)
{
  if(x[0] < 1) TF1::RejectPoint(); // avoid fitting to first point which is added bunches
  return p[0];
}


void bunchPlots(TDirectory* inDir, std::vector<TF1*> fits, string titleStr, bool saveImages)
{
  auto graphDir = inDir->mkdir("Graphs");
  auto histDir = inDir->mkdir("Hists");

  graphDir->cd();

  string plotName;
  string xAxisTitle;
  bool vsTimeParam = false;

  xAxisTitle = "Bunch Number";
  plotName = "BunchNum";

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
    double xPoint = fitNum;

    Chi2NDF_Vs_Iter->SetPoint(fitNum, xPoint, fits.at(fitNum)->GetChisquare()/fits.at(fitNum)->GetNDF());
    Chi2NDF_Vs_Iter->SetPointError(fitNum, 0, sqrt(2./fits.at(fitNum)->GetNDF()));

    Pval_Vs_Iter->SetPoint(fitNum, xPoint, fits.at(fitNum)->GetProb());
  }

  if(vsTimeParam){
    nsTOus(Chi2NDF_Vs_Iter, xAxisTitle + " [#mus]");
    nsTOus(Pval_Vs_Iter, xAxisTitle + " [#mus]");
  }

  graphDir->cd();
    Chi2NDF_Vs_Iter->Write();
    Pval_Vs_Iter->Write();

  histDir->cd();
  
  auto tempGraph = (TGraph*) Chi2NDF_Vs_Iter->Clone(); // remove point for added bunch numbers before making histogram
  tempGraph->RemovePoint(0);
    TH1F* Chi2NDF_Vs_Iter_hist = projectGraphToHist(tempGraph);
      Chi2NDF_Vs_Iter_hist->Write();
    
    TH1F* Pval_Vs_Iter_hist = projectGraphToHist(Pval_Vs_Iter);
      Pval_Vs_Iter_hist->Write();

  if(saveImages){
    string canvasName = (fitTitle + "_Chi2NDF_Vs_" + plotName + "_Canv");
    auto chi2_canv = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,500,400);
    Chi2NDF_Vs_Iter->GetXaxis()->SetLimits(-1, fits.size());
    Chi2NDF_Vs_Iter->Draw("AP");

      double allBunches_result_chi2, tempX2;
      Chi2NDF_Vs_Iter->GetPoint(0, tempX2, allBunches_result_chi2);
      TLine *line = new TLine(-1, allBunches_result_chi2, 9, allBunches_result_chi2);
      line->SetLineColor(4);
      line->SetLineStyle(2);
      line->SetLineWidth(3);
      line->Draw();

    chi2_canv->SaveAs(("Images/" + canvasName + datasetTagForPlots + ".png").c_str());

    auto chi2_hist_canv = new TCanvas((canvasName + "_hist").c_str(),canvasName.c_str(),200,10,500,400);
    Chi2NDF_Vs_Iter_hist->Draw();

      TPaveText* textBox = new TPaveText(0.75,0.7,0.975,0.9, "NDC");
      textBox->SetBorderSize(1);
      textBox->SetFillStyle(1001);
      textBox->AddText(Form("Entries  %1.0f", Chi2NDF_Vs_Iter_hist->GetEntries()));
      textBox->AddText(Form("Mean  %0.4g", Chi2NDF_Vs_Iter_hist->GetMean()));
      textBox->AddText(Form("RMS  %0.4g", Chi2NDF_Vs_Iter_hist->GetRMS()));
      textBox->Draw("SAME");

    chi2_hist_canv->SaveAs(("Images/" + canvasName + "_hist" + datasetTagForPlots + ".png").c_str());

    canvasName = fitTitle + "_PVal_Vs_" + plotName + "_Canv";
    auto pval_canv = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,500,400);
    Pval_Vs_Iter->GetXaxis()->SetLimits(-1, fits.size());
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

    bool parFixed = false; // only for individual bunches
      
      for (uint fitNum = 0; fitNum < fits.size(); ++fitNum)
      {
        double xPoint = fitNum;

        double parameterValue = fits.at(fitNum)->GetParameter(parNum);
        double parameterError = fits.at(fitNum)->GetParError(parNum);
        
        if(phaseParam) normalizePhase(parameterValue);

        param_Vs_Iter->SetPoint(fitNum, xPoint, parameterValue);
        param_Vs_Iter->SetPointError(fitNum, 0, parameterError);

        if(fitNum > 0 && parameterError == 0) parFixed = true;
      }

    if(vsTimeParam) nsTOus(param_Vs_Iter, xAxisTitle + " [#mus]");

    TF1* fitLine = new TF1("fitLine", lineFunc, -1, 9, 1);
    fitLine->SetLineColor(2);
    fitLine->SetParameter(0, fits.at(0)->GetParameter(parNum));

    if(!parFixed) param_Vs_Iter->Fit(fitLine, "QR");


  graphDir->cd();
    param_Vs_Iter->Write();

  histDir->cd();

  auto tempGraph_param = (TGraph*) param_Vs_Iter->Clone(); // remove point for added bunch numbers before making histogram
  tempGraph_param->RemovePoint(0);
    TH1F* param_Vs_Iter_hist = projectGraphToHist(tempGraph_param);
      param_Vs_Iter_hist->Write();

    if(saveImages){
      string canvasName = (fitTitle + "_" + nameString + "_Vs_" + plotName + "_Canv");
      auto graphCanvas = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,500,400);
      param_Vs_Iter->GetXaxis()->SetLimits(-1, fits.size());
      param_Vs_Iter->Draw("AP");

        double allBunches_result, tempX;
        param_Vs_Iter->GetPoint(0, tempX, allBunches_result);
        TLine *line = new TLine(-1, allBunches_result, 9, allBunches_result);
        line->SetLineColor(4);
        line->SetLineStyle(2);
        line->SetLineWidth(3);
        line->Draw();

        // TLine *line_fitMinusSigma = new TLine(-1, fitLine->GetParameter(0) - fitLine->GetParError(0), 9, fitLine->GetParameter(0) - fitLine->GetParError(0));
        // line_fitMinusSigma->SetLineColor(2);
        // line_fitMinusSigma->SetLineStyle(2);
        // line_fitMinusSigma->SetLineWidth(3);
        // line_fitMinusSigma->Draw();

        // TLine *line_fitPlusSigma = new TLine(-1, fitLine->GetParameter(0) + fitLine->GetParError(0), 9, fitLine->GetParameter(0) + fitLine->GetParError(0));
        // line_fitPlusSigma->SetLineColor(2);
        // line_fitPlusSigma->SetLineStyle(2);
        // line_fitPlusSigma->SetLineWidth(3);
        // line_fitPlusSigma->Draw();

      if(!parFixed){
        graphCanvas->Update();
        TPaveStats *ps = (TPaveStats*)graphCanvas->GetPrimitive("stats");
        ps->SetBorderSize(1);
        graphCanvas->Modified();

        ps->Draw("SAME");
      }

        TPaveText* textBox = new TPaveText(0.25,0.8,0.55,1, "NDC");
        textBox->AddText(Form("All Bunches :   %0.4g", allBunches_result));
        textBox->SetBorderSize(0);
        textBox->SetTextColor(4);
        textBox->SetFillStyle(0);
        textBox->Draw("SAME");

        graphCanvas->SaveAs(("Images/" + canvasName + datasetTagForPlots + ".png").c_str());
      /*
      auto histCanvas = new TCanvas((canvasName + "_hist").c_str(),canvasName.c_str(),200,10,500,400);
        param_Vs_Iter_hist->Draw();

        TPaveText* textBox = new TPaveText(0.75,0.7,0.975,0.9, "NDC");
        textBox->SetBorderSize(1);
        textBox->SetFillStyle(1001);
        textBox->AddText(Form("Entries  %1.0f", param_Vs_Iter_hist->GetEntries()));
        textBox->AddText(Form("Mean  %0.4g", param_Vs_Iter_hist->GetMean()));
        textBox->AddText(Form("RMS  %0.4g", param_Vs_Iter_hist->GetRMS()));
        textBox->Draw("SAME");

        histCanvas->SaveAs(("Images/" + canvasName + "_hist" + datasetTagForPlots + ".png").c_str());
      */
    }

  }

}

/////////////////////////////////////////////////////////////////////////////////////

int bunchNumPlots(std::string filePath)
{
  
  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  // check if dataset tag exists in file and if so append it to the file name and write it to the file

  string outputFileName = "bunchNumPlots";
  TNamed* tag = applyDatasetTag(inputFile, outputFileName);

  TFile* outputFile = new TFile((outputFileName + ".root").c_str(),"RECREATE");
  if(tag) tag->Write();
  auto topDir = outputFile->mkdir("topDir");

/////////////////////////////////////////////////////////////////////////////////////
  // These only get set for the interactive root session (any generated canvases, etc.), but does not apply to the output root file - that comes from .rootlogon.C
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(1111);
  // gStyle->SetOptFit(11);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1);
  gStyle->SetLineColor(1);
  gStyle->SetPadRightMargin(.05);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
      // Plots over different fit conditions

  auto addedDir = topDir->mkdir("Added");

std::vector<TF1*> TmethodFits;
std::vector<TF1*> fullRatioFits;

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


  auto TmethodDir = addedDir->mkdir("TMethod");
  bunchPlots(TmethodDir, TmethodFits, "T Method Fit", saveImagesDirectly);

if(fullRatioFits.size() > 0) {
  auto ratioDir = addedDir->mkdir("FullRatio");
  bunchPlots(ratioDir, fullRatioFits, "Full Ratio Fit", saveImagesDirectly);
}

/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

  delete outputFile;

  return 1;
}
