// 3-31-20: Macro to make plots for scans done over ad hoc gain parameters. The amplitude, lifetime, and asymmetry can be scanned over, and the related boolean value down below needs to be turned to true in order to save the plots directly to png format.
// The scans over the others will be included in the output root file, but won't be saved directly.

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
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
bool doParabFit = true; // initial value doesn't matter, gets set later
bool saveImagesDirectly_amplitude = true;
bool saveImagesDirectly_lifetime = false;
bool saveImagesDirectly_asymmetry = false;

pair<double, double> range_amplitude {-2e-3, 3e-3};
// pair<double, double> range_lifetime {,};
pair<double, double> range_asymmetry = {-0.1, 0.6};

TF1* lineFunc = new TF1("lineFunc", "[0] + (x * [1])", range_amplitude.first, range_amplitude.second);
TF1* parabFunc = new TF1("parabFunc", "[2] * (x - [1]) * (x - [1]) + [0]", range_amplitude.first, range_amplitude.second);

void generalScan(TDirectory* inDir, std::tuple<string, string, string> tupleID, std::vector<float> xVector, std::vector<TF1*> fits, bool saveImages, TMultiGraph* multiGraph_R, TMultiGraph* multiGraph_chisq)
{
  inDir->cd();

  bool performFits = ((xVector.back() - xVector.front() != 0)) ? true : false;

    string fitName = std::get<0>(tupleID);
    string scanType = std::get<1>(tupleID);
    string titleStr = std::get<2>(tupleID);

    string scanName = removeSpaces(scanType);

    string axisPart = "";
    if(scanType.find("Time") != string::npos) axisPart = " (ns)";

    
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

      if(performFits && doParabFit){
        parabFunc->SetParameter(0,fits.at(0)->GetChisquare());
        // parabFunc->SetParameter(1,0.6e-3);
        // parabFunc->SetParameter(2,5e6);
        // parabFunc->SetParameter(1,1e-3);
        // parabFunc->SetParameter(2,3e7);
        parabFunc->SetParameter(1,(xVector.at(0)+xVector.at(xVector.size()-1))/2);
        parabFunc->SetParameter(2, 5000./(xVector.at(xVector.size()-1)+xVector.at(0))); // 20 for asymmetry instead of 5000?
        Chi2_Vs_Val->Fit("parabFunc", "QR");
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

      multiGraph_chisq->Add((TGraphErrors*) Chi2_Vs_Val->Clone(("chisq_"+fitName).c_str()));

      graphCanvas->Update();

      if(graphCanvas->GetPrimitive("stats") != 0){
        TPaveStats *statsBox = (TPaveStats*)graphCanvas->GetPrimitive("stats");
        statsBox->SetBorderSize(1);

        statsBox->SetX1NDC(0.4);
        statsBox->SetX2NDC(0.7);

        // statsBox->SetX2(graphCanvas->GetUxmax()); // can reposition and by extension resize the stats box with these commands - not the best though, there must be a better way
        // statsBox->SetY2(graphCanvas->GetUymax());

        statsBox->Draw("SAME"); 
      }

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

      for (uint fitNum = 0; fitNum < fits.size(); ++fitNum)
      {
        double parameterValue = fits.at(fitNum)->GetParameter(parNum);
        double parameterError = fits.at(fitNum)->GetParError(parNum);
        
        if(phaseParam) normalizePhase(parameterValue);

        param_Vs_Val->SetPoint(fitNum, xVector.at(fitNum), parameterValue);
        param_Vs_Val->SetPointError(fitNum, 0, parameterError);
      }

    param_Vs_Val->GetXaxis()->SetTitle((scanType + axisPart).c_str());
    param_Vs_Val->GetYaxis()->SetTitle((paramString + strAddon).c_str());
    param_Vs_Val->GetYaxis()->SetTitleOffset(2.4);

    if(paramString.compare("R") == 0 && performFits){
      param_Vs_Val->Fit(lineFunc, "QR"); // fit for the systematic effect on R
      param_Vs_Val->GetFunction("lineFunc")->SetLineColor(2);
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

      multiGraph_R->Add((TGraphErrors*) param_Vs_Val->Clone(("R_"+fitName).c_str()));

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


int adHocGainScan(std::string filePath)
{
  // gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  if((saveImagesDirectly_amplitude? 1:0) + (saveImagesDirectly_lifetime? 1:0) + (saveImagesDirectly_asymmetry? 1:0) > 1)
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

  string outputFileName = "adHocGainScan";
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

  gStyle->SetPadRightMargin(.1);
  gStyle->SetPadLeftMargin(.2);

  // gStyle->SetOptStat("m");

  lineFunc->SetNpx(10000);
  parabFunc->SetNpx(10000);


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  auto addedDir = topDir->mkdir("Added");

/////////////////////////////////////////////////////////////////////////////////////

// values scanned over

// std::vector<float> ADTs; // can't iterate over this yet
std::vector<float> amplitudes;
std::vector<float> lifetimes;
std::vector<float> asymmetries;

std::vector<TF1*> fullRatioFits;
std::vector<TF1*> TMethodFits;

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
      amplitudes.push_back(float((*histSavedParameters)[13]));
      lifetimes.push_back(float((*histSavedParameters)[14]));
      asymmetries.push_back(float((*histSavedParameters)[15]));

      TF1* TmethodFitFunction = (TF1*) ((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/TMethod/allTimesAdded_TMethod", fitPass)))->GetFunction("TmethodFitFunc")->Clone();
      TMethodFits.push_back(TmethodFitFunction);

      TF1* fullRatioFitFunction = (TF1*) ((TGraphErrors*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass)))->GetFunction("fullRatioFitFunc")->Clone();
      fullRatioFits.push_back(fullRatioFitFunction);
    }

/////////////////////////////////////////////////////////////////////////////////////

    TMultiGraph* ampMultiGraph_R = new TMultiGraph();
    TMultiGraph* ampMultiGraph_chisq = new TMultiGraph();
    TMultiGraph* tauMultiGraph_R = new TMultiGraph();
    TMultiGraph* tauMultiGraph_chisq = new TMultiGraph();
    TMultiGraph* asymmetryMultiGraph_R = new TMultiGraph();
    TMultiGraph* asymmetryMultiGraph_chisq = new TMultiGraph();

/////////////////////////////////////////////////////////////////////////////////////

    auto added_TMethod_dir = addedDir->mkdir("TMethod");

      lineFunc->SetRange(range_amplitude.first, range_amplitude.second);
      parabFunc->SetRange(range_amplitude.first, range_amplitude.second);
      doParabFit = true;

    auto added_TMethod_amp_dir = added_TMethod_dir->mkdir("AdHoc_Amplitude");
    std::tuple<string, string, string> amp_info_TMethod("TMethod", "AdHoc Amplitude", "T Method Fit");
    generalScan(added_TMethod_amp_dir, amp_info_TMethod, amplitudes, TMethodFits, saveImagesDirectly_amplitude, ampMultiGraph_R, ampMultiGraph_chisq);

      // lineFunc->SetRange(lineFunc_shortRange.first, lineFunc_shortRange.second);
      // parabFunc->SetRange(parabFunc_shortRange.first, parabFunc_shortRange.second);
      // doParabFit = false;

    // auto added_TMethod_tau_dir = added_TMethod_dir->mkdir("AdHoc_Lifetime");
    // std::tuple<string, string, string> tau_info_TMethod("TMethod", "AdHoc Lifetime", "T Method Fit");
    // generalScan(added_TMethod_tau_dir, tau_info_TMethod, lifetimes, TMethodFits, saveImagesDirectly_lifetime, tauMultiGraph_R, tauMultiGraph_chisq);

      lineFunc->SetRange(range_asymmetry.first, range_asymmetry.second);
      parabFunc->SetRange(range_asymmetry.first, range_asymmetry.second);
      doParabFit = true;

    auto added_TMethod_asymmetry_dir = added_TMethod_dir->mkdir("AdHoc_Amsymmetry");
    std::tuple<string, string, string> asymmetry_info_TMethod("TMethod", "AdHoc Asymmetry", "T Method Fit");
    generalScan(added_TMethod_asymmetry_dir, asymmetry_info_TMethod, asymmetries, TMethodFits, saveImagesDirectly_asymmetry, asymmetryMultiGraph_R, asymmetryMultiGraph_chisq);

/////////////////////////////////////////////////////////////////////////////////////

    auto added_FullRatio_dir = addedDir->mkdir("FullRatio");

      lineFunc->SetRange(range_amplitude.first, range_amplitude.second);
      parabFunc->SetRange(range_amplitude.first, range_amplitude.second);
      doParabFit = true;

    auto added_FullRatio_amp_dir = added_FullRatio_dir->mkdir("AdHoc_Amplitude");
    std::tuple<string, string, string> amp_info_ratio("FullRatio", "AdHoc Amplitude", "Full Ratio Fit");
    generalScan(added_FullRatio_amp_dir, amp_info_ratio, amplitudes, fullRatioFits, saveImagesDirectly_amplitude, ampMultiGraph_R, ampMultiGraph_chisq);

      // lineFunc->SetRange(lineFunc_shortRange.first, lineFunc_shortRange.second);
      // parabFunc->SetRange(parabFunc_shortRange.first, parabFunc_shortRange.second);
      // doParabFit = false;

    // auto added_FullRatio_tau_dir = added_FullRatio_dir->mkdir("AdHoc_Lifetime");
    // std::tuple<string, string, string> tau_info_ratio("FullRatio", "AdHoc Lifetime", "Full Ratio Fit");
    // generalScan(added_FullRatio_tau_dir, tau_info_ratio, lifetimes, fullRatioFits, saveImagesDirectly_lifetime, tauMultiGraph_R, tauMultiGraph_chisq);

      lineFunc->SetRange(range_asymmetry.first, range_asymmetry.second);
      parabFunc->SetRange(range_asymmetry.first, range_asymmetry.second);
      doParabFit = true;

    auto added_FullRatio_asymmetry_dir = added_FullRatio_dir->mkdir("AdHoc_Asymmetry");
    std::tuple<string, string, string> asymmetry_info_ratio("FullRatio", "AdHoc Asymmetry", "Full Ratio Fit");
    generalScan(added_FullRatio_asymmetry_dir, asymmetry_info_ratio, asymmetries, fullRatioFits, saveImagesDirectly_asymmetry, asymmetryMultiGraph_R, asymmetryMultiGraph_chisq);

/////////////////////////////////////////////////////////////////////////////////////
/*
    gStyle->SetOptFit(0);

    auto compareDir = addedDir->mkdir("Compare");

    vector<string> legendStrings {"T Method", "R Method"};
    TList* grlist = 0;
    TObject *obj = 0;
    TGraph *gr = 0;
    int graphNum = 0;

/////////////////////////////////////////////////////////////////////////////////////

    lineFunc->SetRange(lineFunc_fullRange.first, lineFunc_fullRange.second);
    parabFunc->SetRange(parabFunc_fullRange.first, parabFunc_fullRange.second);

    auto compare_amp_dir = compareDir->mkdir("AdHoc_Amplitude");
    compare_amp_dir->cd();

    if(saveImagesDirectly_amplitude){

      auto legend_chi2_amp = new TLegend(0.25,0.6,.55,0.8);

          grlist = ampMultiGraph_chisq->GetListOfGraphs();
          TIter next_chisq_amp(grlist);
          
          while ((obj = next_chisq_amp())) {
            gr=(TGraph*)obj;
            gr->SetMarkerColor(graphNum+1);
            normalizeGraphToNominal(gr);

            if(graphNum == 0){ // only fit T method for now
              parabFunc->SetParameter(0,4000);
              parabFunc->SetParameter(1,1);
              parabFunc->SetParameter(2,14);
              gr->Fit("parabFunc", "QR");
              gr->GetFunction("parabFunc")->SetLineColor(graphNum+1);
            }

            legend_chi2_amp->AddEntry(gr, legendStrings.at(graphNum).c_str(), "p");

            graphNum++;
          } // end loop over graphs

        ampMultiGraph_chisq->Write("Amplitude_Multigraph_ChiSq");

        auto chi2Canvas_compare_amp = new TCanvas("chi2Canvas_compare_amp","chi2Canvas_compare_amp",1300,410,500,400);
        ampMultiGraph_chisq->Draw("APX");

        adjustGraphRanges(ampMultiGraph_chisq);
        ampMultiGraph_chisq->GetYaxis()->SetTitle("#chi^{2} - #chi^{2}_{nominal}");
        ampMultiGraph_chisq->GetXaxis()->SetTitle("AdHoc Amplitude");

        gPad->Modified();

        legend_chi2_amp->SetBorderSize(0);
        legend_chi2_amp->Draw("SAME");

        chi2Canvas_compare_amp->SaveAs("Images/AdHoc_Amplitude_Compare_Chisq.png");

/////////////////////////////////////////////////////////////////////////////////////

      auto legend_R_amp = new TLegend(0.25,0.6,.65,0.84);

        grlist = ampMultiGraph_R->GetListOfGraphs();
          TIter next_R(grlist);
          graphNum = 0;

          while ((obj = next_R())) {
            gr=(TGraph*)obj;
            gr->SetMarkerColor(graphNum+1);
            gr->SetLineColor(graphNum+1);
            normalizeGraphToNominal(gr);

              gr->Fit(lineFunc, "QR");
              gr->GetFunction("lineFunc")->SetLineColor(graphNum+1);

              double slope = lineFunc->GetParameter(1);

            legend_R_amp->AddEntry(gr, (legendStrings.at(graphNum) + " : " + to_string(slope)).c_str(), "pl");

            graphNum++;
          } // end loop over graphs

        ampMultiGraph_R->Write("Amplitude_Multigraph_R");

        auto Rcanvas_compare_amp = new TCanvas("Rcanvas_compare_amp","Rcanvas_compare_amp",1300,410,500,400);
        ampMultiGraph_R->Draw("APX");

        adjustGraphRanges(ampMultiGraph_R);
        ampMultiGraph_R->GetYaxis()->SetTitle("R - R_{nominal} (ppm)");
        ampMultiGraph_R->GetXaxis()->SetTitle("AdHoc Amplitude");

        gPad->Modified();

        legend_R_amp->SetBorderSize(0);
        legend_R_amp->Draw("SAME");

        Rcanvas_compare_amp->SaveAs("Images/AdHoc_Amplitude_Compare_R.png");
    }

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

    lineFunc->SetRange(lineFunc_shortRange.first, lineFunc_shortRange.second);
    parabFunc->SetRange(parabFunc_shortRange.first, parabFunc_shortRange.second);

    auto compare_tau_dir = compareDir->mkdir("AdHoc_Lifetime");
    compare_tau_dir->cd();

    if(saveImagesDirectly_lifetime){

      auto legend_chi2_tau = new TLegend(0.25,0.6,.55,0.8);

        grlist = tauMultiGraph_chisq->GetListOfGraphs();
          TIter next_chisq_tau(grlist);
          graphNum = 0;
          
          while ((obj = next_chisq_tau())) {
            gr=(TGraph*)obj;
            gr->SetMarkerColor(graphNum+1);
            normalizeGraphToNominal(gr);

            // if(graphNum == 0){ // only fit T method for now
            //   parabFunc->SetParameter(0,4000);
            //   parabFunc->SetParameter(1,1);
            //   parabFunc->SetParameter(2,14);
            //   gr->Fit("parabFunc", "QR");
            //   gr->GetFunction("parabFunc")->SetLineColor(graphNum+1);
            // }

            legend_chi2_tau->AddEntry(gr, legendStrings.at(graphNum).c_str(), "p");

            graphNum++;
          } // end loop over graphs

        tauMultiGraph_chisq->Write("Lifetime_Multigraph_ChiSq");

        auto chi2Canvas_compare_tau = new TCanvas("chi2Canvas_compare_tau","chi2Canvas_compare_tau",1300,410,500,400);
        tauMultiGraph_chisq->Draw("APX");

        adjustGraphRanges(tauMultiGraph_chisq);
        tauMultiGraph_chisq->GetYaxis()->SetTitle("#chi^{2} - #chi^{2}_{nominal}");
        tauMultiGraph_chisq->GetXaxis()->SetTitle("AdHoc Lifetime");

        gPad->Modified();

        legend_chi2_tau->SetBorderSize(0);
        legend_chi2_tau->Draw("SAME");

        chi2Canvas_compare_tau->SaveAs("Images/AdHoc_Lifetime_Compare_Chisq.png");

/////////////////////////////////////////////////////////////////////////////////////

      auto legend_R_tau = new TLegend(0.23,0.6,.64,0.84);

        grlist = tauMultiGraph_R->GetListOfGraphs();
          TIter next_R_tau(grlist);
          graphNum = 0;

          while ((obj = next_R_tau())) {
            gr=(TGraph*)obj;
            gr->SetMarkerColor(graphNum+1);
            gr->SetLineColor(graphNum+1);

            normalizeGraphToNominal(gr);

              gr->Fit(lineFunc, "QR");
              gr->GetFunction("lineFunc")->SetLineColor(graphNum+1);

              double slope = lineFunc->GetParameter(1);

            legend_R_tau->AddEntry(gr, (legendStrings.at(graphNum) + " : " + to_string(slope)).c_str(), "pl");

            graphNum++;
          } // end loop over graphs

        tauMultiGraph_R->Write("Lifetime_Multigraph_R");

        auto RCanvas_compare_tau = new TCanvas("RCanvas_compare_tau","RCanvas_compare_tau",1300,410,500,400);
        tauMultiGraph_R->Draw("APX");

        adjustGraphRanges(tauMultiGraph_R);
        tauMultiGraph_R->GetYaxis()->SetTitle("R - R_{nominal} (ppm)");
        tauMultiGraph_R->GetXaxis()->SetTitle("AdHoc Lifetime");

        gPad->Modified();

        legend_R_tau->SetBorderSize(0);
        legend_R_tau->Draw("SAME");

        RCanvas_compare_tau->SaveAs("Images/AdHoc_Lifetime_Compare_R.png");
    }
*/
/////////////////////////////////////////////////////////////////////////////////////

      delete outputFile;


  return 1;
}
