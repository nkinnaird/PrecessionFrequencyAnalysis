// this macro compares average fit parameters between two different sets of fits to different random seeds - basically like the normal random seed parameter plotter but now the differences - see generalIterPlotter.C

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

int plotVsParam = -1; // set to a positive int to plot against a specific parameter, otherwise set to a negative number to plot vs iteration/random seed, only applies to ratio cbo fit right now - could make it apply to others but then the parameter number is different for the different fits
bool saveImagesDirectly = true; // this only applies to the ratio cbo fit right now

void iterationPlots(TDirectory* inDir, std::vector<TF1*> fits, std::vector<TF1*> fits_second, string titleStr, int VsParamNum, bool saveImages)
{
  auto graphDir = inDir->mkdir("Graphs");
  auto histDir = inDir->mkdir("Hists");

  graphDir->cd();

  string plotName;
  string xAxisTitle;
  bool vsTimeParam = false;

  if(VsParamNum < 0){
    xAxisTitle = "Random Seed";
    // xAxisTitle = "Fit Iteration";
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
  Chi2NDF_Vs_Iter->GetYaxis()->SetTitle("#Delta#chi^{2}/NDF");
  Chi2NDF_Vs_Iter->GetYaxis()->SetTitleOffset(2.0);

  TGraph* Pval_Vs_Iter = new TGraph();
  Pval_Vs_Iter->SetName((fitTitle + "_Pval_Vs_" + plotName).c_str());
  Pval_Vs_Iter->SetTitle((titleStr + " P Value Vs " + xAxisTitle).c_str());
  Pval_Vs_Iter->GetXaxis()->SetTitle(xAxisTitle.c_str());
  Pval_Vs_Iter->GetYaxis()->SetTitle("#DeltaP value");
  Pval_Vs_Iter->GetYaxis()->SetTitleOffset(2.0);
  Pval_Vs_Iter->GetXaxis()->SetRangeUser(0, fits.size()+1);

  for (uint fitNum = 0; fitNum < fits.size(); ++fitNum)
  {
    double xPoint = (VsParamNum < 0) ? fitNum+1 : fits.at(fitNum)->GetParameter(VsParamNum);

    Chi2NDF_Vs_Iter->SetPoint(fitNum, xPoint, fits_second.at(fitNum)->GetChisquare()/fits.at(fitNum)->GetNDF() - fits.at(fitNum)->GetChisquare()/fits.at(fitNum)->GetNDF());
    // Chi2NDF_Vs_Iter->SetPointError(fitNum, 0, sqrt(2./fits.at(fitNum)->GetNDF()));

    Pval_Vs_Iter->SetPoint(fitNum, xPoint, fits_second.at(fitNum)->GetProb() - fits.at(fitNum)->GetProb());
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

    chi2_canv->SaveAs(("Images/" + canvasName + "_compare.png").c_str());

    auto chi2_hist_canv = new TCanvas((canvasName + "_hist").c_str(),canvasName.c_str(),200,10,500,400);
    Chi2NDF_Vs_Iter_hist->Draw();
      textBox->Clear();
      textBox->AddText(Form("Entries  %zu", fits.size()));
      textBox->AddText(Form("Mean  %0.4g", mean));
      textBox->AddText(Form("RMS  %0.4g", rms));
      textBox->Draw("SAME");

    chi2_hist_canv->SaveAs(("Images/" + canvasName + "_hist_compare.png").c_str());

    canvasName = fitTitle + "_PVal_Vs_" + plotName + "_Canv";
    auto pval_canv = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,500,400);
    if(VsParamNum < 0) Pval_Vs_Iter->GetXaxis()->SetLimits(0, fits.size()+1);
    Pval_Vs_Iter->Draw("AP");
    pval_canv->SaveAs(("Images/" + canvasName + "_compare.png").c_str());
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
    param_Vs_Iter->GetYaxis()->SetTitle(("#Delta" + paramString + strAddon).c_str());
    param_Vs_Iter->GetYaxis()->SetTitleOffset(2.0);
    param_Vs_Iter->GetXaxis()->SetTitle(xAxisTitle.c_str());
      
      for (uint fitNum = 0; fitNum < fits.size(); ++fitNum)
      {
        double xPoint = (VsParamNum < 0) ? fitNum+1 : fits.at(fitNum)->GetParameter(VsParamNum);

        double parameterValue = fits.at(fitNum)->GetParameter(parNum);
        double parameterError = fits.at(fitNum)->GetParError(parNum);

        double parameterValue_second = fits_second.at(fitNum)->GetParameter(parNum);
        
        if(phaseParam) normalizePhase(parameterValue);
        // if(phaseParam && parameterValue < 2) parameterValue += 2*pi;

        param_Vs_Iter->SetPoint(fitNum, xPoint, parameterValue_second - parameterValue);
        // param_Vs_Iter->SetPointError(fitNum, 0, parameterError);
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

        graphCanvas->SaveAs(("Images/" + canvasName + "_compare.png").c_str());

      auto histCanvas = new TCanvas((canvasName + "_hist").c_str(),canvasName.c_str(),200,10,500,400);
        param_Vs_Iter_hist->Draw();
        textBox->Clear();
        textBox->AddText(Form("Entries  %zu", fits.size()));
        textBox->AddText(Form("Mean  %0.4g", mean));
        textBox->AddText(Form("RMS  %0.4g", rms));
        textBox->Draw("SAME");

        histCanvas->SaveAs(("Images/" + canvasName + "_hist_compare.png").c_str());
      }
      else if(VsParamNum >= 0 && paramString.compare("R") == 0){
        
        // draw options so shape of R vs parameter curve shows up
        param_Vs_Iter->GetYaxis()->SetRangeUser(fits.back()->GetParameter(1) - .02, fits.front()->GetParameter(1) + .02);
        param_Vs_Iter->Draw("XAP"); // X - draw without error bars

        graphCanvas->Update();
        TPaveStats *statsBox = (TPaveStats*)graphCanvas->GetPrimitive("stats");
        statsBox->SetBorderSize(1);
        statsBox->Draw("SAME");

        graphCanvas->SaveAs(("Images/" + canvasName + "_compare_compare.png").c_str());
      }
    }

  }

}

/////////////////////////////////////////////////////////////////////////////////////

int FitParameterComparison()
{
  // gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  // 60h
  // string firstFile_string = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/60h/RandSeeds/OneBinLater/output-60h-RandSeeds-OneBinLater.root"; // typically uncorrected file
  // string secondFile_string = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/60h/RandSeeds/With-Ad-Hoc/FitIterations/output-60h-RandSeeds-with-AdHoc.root"; // typically corrected file

  // HighKick
  string firstFile_string = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/HighKick/RandSeeds/OneBinLater/output-HighKick-RandSeeds-OneBinLater.root"; // typically uncorrected file
  string secondFile_string = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/HighKick/RandSeeds/withAdHoc/FitIterations/output-HighKick-RandSeeds-withAdHoc.root"; // typically corrected file


  // // 9d
  // string firstFile_string = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/9d/RandSeeds/OneBinLater/output-9d-RandSeeds-OneBinLater.root"; // typically uncorrected file
  // string secondFile_string = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/9d/RandSeeds/withAdHoc/FitIterations/output-9d-RandSeeds-withAdHoc.root"; // typically corrected file

  // // Endgame
  // string firstFile_string = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/Endgame/RandSeeds/OneBinLater/output-Endgame-RandSeeds-OneBinLater.root"; // typically uncorrected file
  // string secondFile_string = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/Endgame/RandSeeds/withAdHoc/FitIterations/output-Endgame-RandSeeds-withAdHoc.root"; // typically corrected file


  TFile *firstFile = TFile::Open(firstFile_string.c_str());
  TFile *secondFile = TFile::Open(secondFile_string.c_str());
   if (firstFile == 0 || secondFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }


  TFile* outputFile = new TFile("fitParameterComparison.root","RECREATE");
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


std::vector<TF1*> TmethodFits_firstFile;
std::vector<TF1*> fullRatioFits_firstFile;
std::vector<TF1*> TmethodFits_secondFile;
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
      TF1* TMethodFitFunction_firstFile = (TF1*) ((TH1F*) firstFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/TMethod/allTimesAdded_TMethod", fitPass)))->GetFunction("TmethodFitFunc")->Clone();
      TF1* TMethodFitFunction_secondFile = (TF1*) ((TH1F*) secondFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/TMethod/allTimesAdded_TMethod", fitPass)))->GetFunction("TmethodFitFunc")->Clone();
      
      TF1* fullRatioFitFunction_firstFile = (TF1*) ((TGraphErrors*) firstFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass)))->GetFunction("fullRatioFitFunc")->Clone();
      TF1* fullRatioFitFunction_secondFile = (TF1*) ((TGraphErrors*) secondFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass)))->GetFunction("fullRatioFitFunc")->Clone();

      TmethodFits_firstFile.push_back(TMethodFitFunction_firstFile);
      TmethodFits_secondFile.push_back(TMethodFitFunction_secondFile);
      fullRatioFits_firstFile.push_back(fullRatioFitFunction_firstFile);
      fullRatioFits_secondFile.push_back(fullRatioFitFunction_secondFile);
    }


  auto TmethodDir = addedDir->mkdir("TMethod");
  iterationPlots(TmethodDir, TmethodFits_firstFile, TmethodFits_secondFile, "T Method Fit", plotVsParam, saveImagesDirectly);

  auto fullRatioDir = addedDir->mkdir("FullRatio");
  iterationPlots(fullRatioDir, fullRatioFits_firstFile, fullRatioFits_secondFile, "Full Ratio Fit", plotVsParam, saveImagesDirectly);


/////////////////////////////////////////////////////////////////////////////////////

  delete outputFile;

  return 1;
}
