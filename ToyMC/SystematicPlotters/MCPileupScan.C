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

// #include <TRatioPlot.h> // wait till root v 6_08

#include "ratioAnalysisDefs.hh"

using namespace std;


void stackHists(std::vector<TH1F*> hists, string axisString, std::pair<double, double> xRange, std::pair<double, double> yRange)
{

  THStack* histStack = new THStack("histStack","Temp Title"); // use a THStack because histogram attributes will be preserved when opening the canvas in a TBrowser
  histStack->SetTitle(hists.at(0)->GetName()); // or could use GetTitle

  for (uint i = 0; i < hists.size(); ++i)
  {
    hists.at(i)->SetLineColor(hists.size()-i);
    histStack->Add(hists.at(i));
  }


  string canvasName = "stacked_";
  canvasName.append(hists.at(0)->GetName());

  auto stackedCanvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 200, 10, 1200, 1000);
  histStack->Draw("nostack,hist");

  histStack->GetXaxis()->SetTitle(axisString.c_str());
  histStack->GetYaxis()->SetTitle("Events");

  histStack->GetXaxis()->SetRangeUser(xRange.first, xRange.second);
  histStack->SetMinimum(yRange.first); //SetRangeUser for Y doesn't want to work so do this instead
  histStack->SetMaximum(yRange.second);

  stackedCanvas->Write();
  delete stackedCanvas;
}

void fiveParamScan(TDirectory* inDir, std::tuple<string, string, string> tupleID, std::vector<float> xVector, std::vector<TF1*> fiveFits)
{
  auto newDir = inDir->mkdir(std::get<0>(tupleID).c_str());
  newDir->cd();

/////////////////////////////////////////////////////////////////////////////////////


    TGraphErrors* Five_Chi2_Vs_Val = new TGraphErrors();
      Five_Chi2_Vs_Val->SetName(("Five_Chi2_Vs_" + std::get<0>(tupleID)).c_str());
      Five_Chi2_Vs_Val->SetTitle(("Five #chi^{2} Vs " + std::get<1>(tupleID)).c_str());

    TGraphErrors* Five_Chi2NDF_Vs_Val = new TGraphErrors();
      Five_Chi2NDF_Vs_Val->SetName(("Five_Chi2NDF_Vs_" + std::get<0>(tupleID)).c_str());
      Five_Chi2NDF_Vs_Val->SetTitle(("Five #chi^{2}/NDF Vs " + std::get<1>(tupleID)).c_str());

    TGraphErrors* Five_Pval_Vs_Val = new TGraphErrors();
      Five_Pval_Vs_Val->SetName(("Five_Pval_Vs_" + std::get<0>(tupleID)).c_str());
      Five_Pval_Vs_Val->SetTitle(("Five P Value Vs " + std::get<1>(tupleID)).c_str());

/////////////////////////////////////////////////////////////////////////////////////

    TGraphErrors* Five_N_Vs_Val = new TGraphErrors();
      Five_N_Vs_Val->SetName(("Five_N_Vs_" + std::get<0>(tupleID)).c_str());
      Five_N_Vs_Val->SetTitle(("Five N Vs " + std::get<1>(tupleID)).c_str());

    TGraphErrors* Five_Tau_Vs_Val = new TGraphErrors();
      Five_Tau_Vs_Val->SetName(("Five_Tau_Vs_" + std::get<0>(tupleID)).c_str());
      Five_Tau_Vs_Val->SetTitle(("Five #tau (ppm) Vs " + std::get<1>(tupleID)).c_str());

    TGraphErrors* Five_A_Vs_Val = new TGraphErrors();
      Five_A_Vs_Val->SetName(("Five_A_Vs_" + std::get<0>(tupleID)).c_str());
      Five_A_Vs_Val->SetTitle(("Five A Vs " + std::get<1>(tupleID)).c_str());

    TGraphErrors* Five_R_Vs_Val = new TGraphErrors();
      Five_R_Vs_Val->SetName(("Five_R_Vs_" + std::get<0>(tupleID)).c_str());
      Five_R_Vs_Val->SetTitle(("Five R (ppm) Vs " + std::get<1>(tupleID)).c_str());

    TGraphErrors* Five_Phi_Vs_Val = new TGraphErrors();
      Five_Phi_Vs_Val->SetName(("Five_Phi_Vs_" + std::get<0>(tupleID)).c_str());
      Five_Phi_Vs_Val->SetTitle(("Five #phi Vs " + std::get<1>(tupleID)).c_str());

/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////


int pointNo = 0;

    for (float xVal : xVector)
    {
      double chi2 = fiveFits.at(pointNo)->GetChisquare();
      double chi2ndf = chi2/fiveFits.at(pointNo)->GetNDF();
      double pVal = fiveFits.at(pointNo)->GetProb();

/////////////////////////////////////////////////////////////////////////////////////

      double Ratio_N = fiveFits.at(pointNo)->GetParameter(0);
      double Ratio_N_Error = fiveFits.at(pointNo)->GetParError(0);
      double Ratio_Tau = fiveFits.at(pointNo)->GetParameter(1);
      double Ratio_Tau_Error = fiveFits.at(pointNo)->GetParError(1);
      double Ratio_A = fiveFits.at(pointNo)->GetParameter(2);
      double Ratio_A_Error = fiveFits.at(pointNo)->GetParError(2);
      double Ratio_R = fiveFits.at(pointNo)->GetParameter(3);
      double Ratio_R_Error = fiveFits.at(pointNo)->GetParError(3);
      double Ratio_Phi = fiveFits.at(pointNo)->GetParameter(4);
      double Ratio_Phi_Error = fiveFits.at(pointNo)->GetParError(4);

        while(Ratio_Phi < 0) Ratio_Phi += 2*pi;

        Ratio_Phi = fmod(Ratio_Phi, 2*pi);

/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

        Five_Chi2_Vs_Val->SetPoint(pointNo, xVal, chi2);
        Five_Chi2NDF_Vs_Val->SetPoint(pointNo, xVal, chi2ndf);
        Five_Pval_Vs_Val->SetPoint(pointNo, xVal, pVal);

/////////////////////////////////////////////////////////////////////////////////////

        Five_N_Vs_Val->SetPoint(pointNo, xVal, Ratio_N);
        Five_N_Vs_Val->SetPointError(pointNo, 0, Ratio_N_Error);

        Five_Tau_Vs_Val->SetPoint(pointNo, xVal, Ratio_Tau);
        Five_Tau_Vs_Val->SetPointError(pointNo, 0, Ratio_Tau_Error);

        Five_A_Vs_Val->SetPoint(pointNo, xVal, Ratio_A);
        Five_A_Vs_Val->SetPointError(pointNo, 0, Ratio_A_Error);

        Five_R_Vs_Val->SetPoint(pointNo, xVal, Ratio_R);
        Five_R_Vs_Val->SetPointError(pointNo, 0, Ratio_R_Error);

        Five_Phi_Vs_Val->SetPoint(pointNo, xVal, Ratio_Phi);
        Five_Phi_Vs_Val->SetPointError(pointNo, 0, Ratio_Phi_Error);


/////////////////////////////////////////////////////////////////////////////////////

      pointNo++;

    } // end loop over point numbers

/////////////////////////////////////////////////////////////////////////////////////
        
        Five_Chi2_Vs_Val->GetXaxis()->SetTitle((std::get<1>(tupleID) + " " + std::get<2>(tupleID)).c_str());
        Five_Chi2_Vs_Val->GetYaxis()->SetTitle("#chi^{2}");
          TF1* parabFunc = new TF1("parabFunc", "[2] * (x - [1]) * (x - [1]) + [0]", .8, 1.2);
          parabFunc->SetParameter(0,fiveFits.at(0)->GetChisquare());
          parabFunc->SetParameter(1,1);
          parabFunc->SetParameter(2,fiveFits.at(0)->GetChisquare());
      Five_Chi2_Vs_Val->Fit("parabFunc", "Q");
      Five_Chi2_Vs_Val->GetFunction("parabFunc")->SetLineColor(2);
      Five_Chi2_Vs_Val->Write();

        Five_Chi2NDF_Vs_Val->GetXaxis()->SetTitle((std::get<1>(tupleID) + " " + std::get<2>(tupleID)).c_str());
        Five_Chi2NDF_Vs_Val->GetYaxis()->SetTitle("#chi^{2}/NDF");
      Five_Chi2NDF_Vs_Val->Write();

        Five_Pval_Vs_Val->GetXaxis()->SetTitle((std::get<1>(tupleID) + " " + std::get<2>(tupleID)).c_str());
        Five_Pval_Vs_Val->GetYaxis()->SetTitle("P value");
      Five_Pval_Vs_Val->Write();

/////////////////////////////////////////////////////////////////////////////////////

        Five_N_Vs_Val->GetXaxis()->SetTitle((std::get<1>(tupleID) + " " + std::get<2>(tupleID)).c_str());
        Five_N_Vs_Val->GetYaxis()->SetTitle("N");
      Five_N_Vs_Val->Write();

        Five_Tau_Vs_Val->GetXaxis()->SetTitle((std::get<1>(tupleID) + " " + std::get<2>(tupleID)).c_str());
        Five_Tau_Vs_Val->GetYaxis()->SetTitle("#tau");
      Five_Tau_Vs_Val->Write();

        Five_A_Vs_Val->GetXaxis()->SetTitle((std::get<1>(tupleID) + " " + std::get<2>(tupleID)).c_str());
        Five_A_Vs_Val->GetYaxis()->SetTitle("A");
      Five_A_Vs_Val->Write();

        Five_R_Vs_Val->GetXaxis()->SetTitle((std::get<1>(tupleID) + " " + std::get<2>(tupleID)).c_str());
        Five_R_Vs_Val->GetYaxis()->SetTitle("R (ppm)");
      Five_R_Vs_Val->Fit("pol1", "Q"); // fit for the systematic effect on R
      Five_R_Vs_Val->GetFunction("pol1")->SetLineColor(2);
      Five_R_Vs_Val->Write();

        Five_Phi_Vs_Val->GetXaxis()->SetTitle((std::get<1>(tupleID) + " " + std::get<2>(tupleID)).c_str());
        Five_Phi_Vs_Val->GetYaxis()->SetTitle("#phi");
      Five_Phi_Vs_Val->Write();


}



void generalScan(TDirectory* inDir, std::tuple<string, string, string> tupleID, std::vector<float> xVector, std::vector<TF1*> ratioFits)
{
  auto newDir = inDir->mkdir(std::get<0>(tupleID).c_str());
  newDir->cd();

/////////////////////////////////////////////////////////////////////////////////////

    TGraphErrors* Ratio_Chi2_Vs_Val = new TGraphErrors();
      Ratio_Chi2_Vs_Val->SetName(("Ratio_Chi2_Vs_" + std::get<0>(tupleID)).c_str());
      Ratio_Chi2_Vs_Val->SetTitle(("Ratio #chi^{2} Vs " + std::get<1>(tupleID)).c_str());

    TGraphErrors* Ratio_Chi2NDF_Vs_Val = new TGraphErrors();
      Ratio_Chi2NDF_Vs_Val->SetName(("Ratio_Chi2NDF_Vs_" + std::get<0>(tupleID)).c_str());
      Ratio_Chi2NDF_Vs_Val->SetTitle(("Ratio #chi^{2}/NDF Vs " + std::get<1>(tupleID)).c_str());

    TGraphErrors* Ratio_Pval_Vs_Val = new TGraphErrors();
      Ratio_Pval_Vs_Val->SetName(("Ratio_Pval_Vs_" + std::get<0>(tupleID)).c_str());
      Ratio_Pval_Vs_Val->SetTitle(("Ratio P Value Vs " + std::get<1>(tupleID)).c_str());

/////////////////////////////////////////////////////////////////////////////////////

    TGraphErrors* Ratio_A_Vs_Val = new TGraphErrors();
      Ratio_A_Vs_Val->SetName(("Ratio_A_Vs_" + std::get<0>(tupleID)).c_str());
      Ratio_A_Vs_Val->SetTitle(("Ratio A Vs " + std::get<1>(tupleID)).c_str());

    TGraphErrors* Ratio_R_Vs_Val = new TGraphErrors();
      Ratio_R_Vs_Val->SetName(("Ratio_R_Vs_" + std::get<0>(tupleID)).c_str());
      Ratio_R_Vs_Val->SetTitle(("Ratio R (ppm) Vs " + std::get<1>(tupleID)).c_str());

        // TGraph* Ratio_R_Err_Vs_Val = new TGraph();
        //   Ratio_R_Err_Vs_Val-(>SetName("Ratio_R_Err_Vs_" + std::get<0>(tupleID)).c_str());
        //   Ratio_R_Err_Vs_Val-(>SetTitle("Ratio R (ppm) Error Vs " + std::get<1>(tupleID)).c_str());

    TGraphErrors* Ratio_Phi_Vs_Val = new TGraphErrors();
      Ratio_Phi_Vs_Val->SetName(("Ratio_Phi_Vs_" + std::get<0>(tupleID)).c_str());
      Ratio_Phi_Vs_Val->SetTitle(("Ratio #phi Vs " + std::get<1>(tupleID)).c_str());

/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////


int pointNo = 0;

    for (float xVal : xVector)
    {
      double chi2 = ratioFits.at(pointNo)->GetChisquare();
      double chi2ndf = chi2/ratioFits.at(pointNo)->GetNDF();
      double pVal = ratioFits.at(pointNo)->GetProb();

/////////////////////////////////////////////////////////////////////////////////////

      double Ratio_A = ratioFits.at(pointNo)->GetParameter(0);
      double Ratio_A_Error = ratioFits.at(pointNo)->GetParError(0);
      double Ratio_R = ratioFits.at(pointNo)->GetParameter(1);
      double Ratio_R_Error = ratioFits.at(pointNo)->GetParError(1);
      double Ratio_Phi = ratioFits.at(pointNo)->GetParameter(2);
      double Ratio_Phi_Error = ratioFits.at(pointNo)->GetParError(2);

        while(Ratio_Phi < 0) Ratio_Phi += 2*pi;

        Ratio_Phi = fmod(Ratio_Phi, 2*pi);

/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

        Ratio_Chi2_Vs_Val->SetPoint(pointNo, xVal, chi2);
        Ratio_Chi2NDF_Vs_Val->SetPoint(pointNo, xVal, chi2ndf);
        Ratio_Pval_Vs_Val->SetPoint(pointNo, xVal, pVal);

/////////////////////////////////////////////////////////////////////////////////////

        Ratio_A_Vs_Val->SetPoint(pointNo, xVal, Ratio_A);
        Ratio_A_Vs_Val->SetPointError(pointNo, 0, Ratio_A_Error);

        Ratio_R_Vs_Val->SetPoint(pointNo, xVal, Ratio_R);
        Ratio_R_Vs_Val->SetPointError(pointNo, 0, Ratio_R_Error);

        Ratio_Phi_Vs_Val->SetPoint(pointNo, xVal, Ratio_Phi);
        Ratio_Phi_Vs_Val->SetPointError(pointNo, 0, Ratio_Phi_Error);


/////////////////////////////////////////////////////////////////////////////////////

      pointNo++;

    } // end loop over point numbers

/////////////////////////////////////////////////////////////////////////////////////

        Ratio_Chi2_Vs_Val->GetXaxis()->SetTitle((std::get<1>(tupleID) + " " + std::get<2>(tupleID)).c_str());
        Ratio_Chi2_Vs_Val->GetYaxis()->SetTitle("#chi^{2}");
          TF1* parabFunc = new TF1("parabFunc", "[2] * (x - [1]) * (x - [1]) + [0]", .8, 1.2);
          parabFunc->SetParameter(0,ratioFits.at(0)->GetChisquare());
          parabFunc->SetParameter(1,1);
          parabFunc->SetParameter(2,ratioFits.at(0)->GetChisquare());
      Ratio_Chi2_Vs_Val->Fit("parabFunc", "Q");
      Ratio_Chi2_Vs_Val->GetFunction("parabFunc")->SetLineColor(2);
      Ratio_Chi2_Vs_Val->Write();

        Ratio_Chi2NDF_Vs_Val->GetXaxis()->SetTitle((std::get<1>(tupleID) + " " + std::get<2>(tupleID)).c_str());
        Ratio_Chi2NDF_Vs_Val->GetYaxis()->SetTitle("#chi^{2}/NDF");
      Ratio_Chi2NDF_Vs_Val->Write();

        Ratio_Pval_Vs_Val->GetXaxis()->SetTitle((std::get<1>(tupleID) + " " + std::get<2>(tupleID)).c_str());
        Ratio_Pval_Vs_Val->GetYaxis()->SetTitle("P value");
      Ratio_Pval_Vs_Val->Write();

/////////////////////////////////////////////////////////////////////////////////////

        Ratio_A_Vs_Val->GetXaxis()->SetTitle((std::get<1>(tupleID) + " " + std::get<2>(tupleID)).c_str());
        Ratio_A_Vs_Val->GetYaxis()->SetTitle("A");
      Ratio_A_Vs_Val->Write();

        Ratio_R_Vs_Val->GetXaxis()->SetTitle((std::get<1>(tupleID) + " " + std::get<2>(tupleID)).c_str());
        Ratio_R_Vs_Val->GetYaxis()->SetTitle("R (ppm)");
      Ratio_R_Vs_Val->Fit("pol1", "Q"); // fit for the systematic effect on R
      Ratio_R_Vs_Val->GetFunction("pol1")->SetLineColor(2);
      Ratio_R_Vs_Val->Write();

        Ratio_Phi_Vs_Val->GetXaxis()->SetTitle((std::get<1>(tupleID) + " " + std::get<2>(tupleID)).c_str());
        Ratio_Phi_Vs_Val->GetYaxis()->SetTitle("#phi");
      Ratio_Phi_Vs_Val->Write();


}


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


int MCPileupScan(std::string filePath)
{
  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  TFile* outputFile = new TFile("MCPileupScan.root","RECREATE");
  auto topDir = outputFile->mkdir("topDir");

  auto ratioDir = topDir->mkdir("Ratio");
  auto fiveDir = topDir->mkdir("FiveParam");


/////////////////////////////////////////////////////////////////////////////////////
  // These only get set for the interactive root session (any generated canvases, etc.), but does not apply to the output root file - that comes from .rootlogon.C
  // gStyle->SetOptStat(000000);
  // gStyle->SetOptTitle(1);
  // // gStyle->SetOptFit(0000);
  // gStyle->SetOptFit(1111);
  // gStyle->SetMarkerStyle(20);
  // gStyle->SetMarkerColor(1);
  // gStyle->SetMarkerSize(.1);
  // gStyle->SetLineColor(1);


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
      // Plots over different fit conditions

  // auto addedDir = topDir->mkdir("Added");
  // auto calosDir = topDir->mkdir("Calos");

/////////////////////////////////////////////////////////////////////////////////////

// values scanned over

// std::vector<float> ADTs; // can't iterate over this yet
std::vector<float> SDTs;
std::vector<float> SGTs;
std::vector<float> PSFs;

std::vector<TF1*> ratioFits;
std::vector<TF1*> fiveFits;

/////////////////////////////////////////////////////////////////////////////////////

  std::vector<TH1F*> pileupTimesAboveThreshold;
  std::vector<TH1F*> energiesPileupSubtracted;



/////////////////////////////////////////////////////////////////////////////////////

int totalIters = 0;

// do this to get the total number of iterations in the file
TIter next(inputFile->GetListOfKeys());
TKey *key;

string rootFileName;

while ((key = (TKey*)next())) {
  string keyName = key->GetName();
  if(keyName.substr(keyName.size()-5) == ".root"){
    rootFileName = keyName;
    TVectorD* itersDouble = ((TVectorD*) inputFile->Get((keyName + "/Iters").c_str()));
    totalIters = (*itersDouble)[0];
  } 
}



/////////////////////////////////////////////////////////////////////////////////////

    // loop through passes and store relevant objects

    for (int fitPass = 0; fitPass < totalIters; ++fitPass)
    {
      TVectorD* histSavedParameters = (TVectorD*) inputFile->Get(Form("topDir/ToyMCIter%d/SubIter3/parameterStore", fitPass));
      SDTs.push_back(float((*histSavedParameters)[2]));
      SGTs.push_back(float((*histSavedParameters)[3]));
      PSFs.push_back(float((*histSavedParameters)[5]));


      TF1* ratioFitFunction = (TF1*) ((TGraphErrors*) inputFile->Get(Form("topDir/ToyMCIter%d/SubIter3/Ratio/Toy_Ratio_Graph", fitPass)))->GetFunction("ratioFitIntegral")->Clone();
      ratioFits.push_back(ratioFitFunction);

      // TF1* fiveParamFunc = (TF1*) ((TH1F*) inputFile->Get(Form("topDir/ToyMCIter%d/SubIter3/Hist/Toy_5_Param_Hist", fitPass)))->GetFunction("fiveParamFit")->Clone();
      // fiveFits.push_back(fiveParamFunc);

/////////////////////////////////////////////////////////////////////////////////////

      pileupTimesAboveThreshold.push_back((TH1F*) inputFile->Get((rootFileName + "/topDir/Iter" + to_string(fitPass) + "/PileupPlots/Shadow/pileupTimes_shadow_threshold").c_str())->Clone());
      // pileupTimesAboveThreshold.push_back((TH1F*) inputFile->Get((rootFileName + "/topDir/Iter" + to_string(fitPass) + "/PileupPlots/Truth/pileupTimes_threshold").c_str())->Clone());

      energiesPileupSubtracted.push_back((TH1F*) inputFile->Get((rootFileName + "/topDir/Iter" + to_string(fitPass) + "/SubIter3/energies").c_str())->Clone());
    }

/////////////////////////////////////////////////////////////////////////////////////

    std::tuple<string, string, string> STDinfo("SDT", "Shadow Dead Time", "(ns)");
    generalScan(ratioDir, STDinfo, SDTs, ratioFits);

    std::tuple<string, string, string> SGTinfo("SGT", "Shadow Gap Time", "(ns)");
    generalScan(ratioDir, SGTinfo, SGTs, ratioFits);

/////////////////////////////////////////////////////////////////////////////////////

    std::tuple<string, string, string> PSFinfo("PSF", "Pileup Scale Factor", "");
    generalScan(ratioDir, PSFinfo, PSFs, ratioFits);

    // fiveParamScan(fiveDir, PSFinfo, PSFs, fiveFits);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/*
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
        PSFs_vs_SDT->GetXaxis()->SetTitle("SDT (ns)");
        PSFs_vs_SDT->GetYaxis()->SetTitle("PSF");

    TGraph* PSFs_vs_SGT = new TGraph();
      PSFs_vs_SGT->SetName("PSFs_vs_SGT");
      PSFs_vs_SGT->SetTitle("PSFs_vs_SGT");
        PSFs_vs_SGT->GetXaxis()->SetTitle("SGT (ns)");
        PSFs_vs_SGT->GetYaxis()->SetTitle("PSF");

    for (uint i = 0; i < PSFs.size(); ++i)
    {
      PSFs_vs_Iter->SetPoint(i, i, PSFs.at(i));
      PSFs_vs_SDT->SetPoint(i, SDTs.at(i), PSFs.at(i));
      PSFs_vs_SGT->SetPoint(i, SGTs.at(i), PSFs.at(i));
    }

    PSFs_vs_Iter->Write("PSFs_vs_Iter");
    PSFs_vs_SDT->Write("PSFs_vs_SDT");
    PSFs_vs_SGT->Write("PSFs_vs_SGT");
*/
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


    auto PileupIterationComparison = topDir->mkdir("PileupIterationComparison");
    PileupIterationComparison->cd();

    stackHists(pileupTimesAboveThreshold, "Time (ns)", make_pair(20000, 60000), make_pair(0, 27000)); // right now these are hardcoded but I might want to grab max and min in Y over the given X range, problematic with the heavy drop off for energy though..
    stackHists(energiesPileupSubtracted, "Energy (MeV)", make_pair(3000, 6000), make_pair(-1000, 4000));


/////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////

      // delete outputFile;


  return 1;
}
