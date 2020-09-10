// Macro to create plots looking at the pileup rate errors.

#include <iostream>
#include <fstream>
#include <iomanip>
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
#include <TRandom3.h>
#include <sstream>
#include <TVirtualFFT.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TPaveStats.h>
#include <TVectorD.h>

#include "ratioAnalysisDefs.hh"
#include "pileupUtils.hh"
#include "plotUtils.hh"

using namespace std;


int pileupRateErrors(){

  TFile *file_60h = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/60h/SingleIter/SingleFit-NewRange/output-60h-SingleIter-NewRange.root");
  TFile *file_60h_manypoints = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/60h/SingleIter/SingleFit-100000pts/output.root");

  // TFile *file_HK = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/HighKick/SingleIter/SingleFit-NewRange/output-HighKick-SingleIter-NewRange.root");
  // TFile *file_9d = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/9d/SingleIter/SingleFit-NewRange/output-9d-SingleIter-NewRange.root");
  // TFile *file_EG = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/Endgame/SingleIter/SingleFit-NewRange/output-Endgame-SingleIter-NewRange.root");

  // TFile* outputFile = new TFile("pileupRateErrors.root","RECREATE");

/////////////////////////////////////////////////////////////////////////////////////

  // auto dir_60h = outputFile->mkdir("60h");
  // dir_60h->cd();

  // grab histogram for binning and whatnot

  auto TmethodHist = ((TH1F*) file_60h->Get("topDir/FitPasses/FitPass0/addedDir/TMethod/allTimesAdded_TMethod")->Clone("TmethodHist"));
  double histBinWidth = TmethodHist->GetBinWidth(1);
  
  TVectorD* histSavedParameters_60h = (TVectorD*) file_60h->Get("topDir/FitPasses/FitPass0/parameterStore");
  float SGT_60h = float((*histSavedParameters_60h)[3]);
  // float SGT_60h = 149.2;

  cout << "SGT_60h: " << SGT_60h << endl;

  double xmin, xmax;

  TF1* Tmethod_fitFunction_60h = (TF1*) TmethodHist->GetFunction("TmethodFitFunc")->Clone();
  Tmethod_fitFunction_60h->GetRange(xmin, xmax);

  cout << "Range: " << xmax - xmin << " integer of bin width: " << (xmax - xmin)/histBinWidth << endl;

  TF1* Tmethod_fitFunction_60h_manypoints = (TF1*) ((TH1F*) file_60h_manypoints->Get("topDir/FitPasses/FitPass0/addedDir/TMethod/allTimesAdded_TMethod"))->GetFunction("TmethodFitFunc")->Clone();
Tmethod_fitFunction_60h_manypoints->SetLineColor(2);

/////////////////////////////////////////////////////////////////////////////////////

  // auto funccanv = new TCanvas("funccanv","funccanv",100,10,1000,800);
  // Tmethod_fitFunction_60h->Draw("PL");
  // Tmethod_fitFunction_60h_manypoints->Draw("SAMEPL");
  // Tmethod_fitFunction_60h_manypoints->Draw("PL");

  /////////////////////////////////////////////////////////////////////////////////////


  TGraph* Tmethod_rateError_60h = new TGraph();
  TGraph* Tmethod_correctionFactor = new TGraph();

  TGraph* Tmethod_rateError_60h_manypoints = new TGraph();
  Tmethod_rateError_60h_manypoints->SetLineColor(2);
  Tmethod_rateError_60h_manypoints->SetMarkerColor(2);

/////////////////////////////////////////////////////////////////////////////////////


  double xPoint = (xmin + histBinWidth/2.);
  double xPointDivisions = 50;

  while(xPoint < xmax)
  // while(xPoint < 30560)
  // while(xPoint < 31560)
  {
    double Tmethod_val_at_t = Tmethod_fitFunction_60h->Eval(xPoint);;
    double Tmethod_val_at_tminus = Tmethod_fitFunction_60h->Eval(xPoint-SGT_60h/2);;
    double Tmethod_val_at_tplus = Tmethod_fitFunction_60h->Eval(xPoint+SGT_60h/2);;

    // double tminus_diff = Tmethod_val_at_tminus - Tmethod_val_at_t;
    // double tplus_diff = Tmethod_val_at_tplus - Tmethod_val_at_t;
    // double yPointNum =  (Tmethod_val_at_t-tminus_diff) * (Tmethod_val_at_t+tminus_diff) - pow(Tmethod_val_at_t,2);
    // double yPointNum =  Tmethod_val_at_tminus*Tmethod_val_at_tplus - pow(Tmethod_val_at_t,2);

    double yPointNum = Tmethod_val_at_t*Tmethod_val_at_t - Tmethod_val_at_tminus*Tmethod_val_at_tplus;
    double yPointDenom = pow(Tmethod_val_at_t,2);

    // cout << setprecision(12) <<  "A -eval: " << Tmethod_val_at_tminus << " +eval: " << Tmethod_val_at_tplus << " middleeval: " << Tmethod_val_at_t << " num: " << yPointNum << endl;

    double yPoint = yPointNum/yPointDenom;
    Tmethod_rateError_60h->SetPoint(Tmethod_rateError_60h->GetN(), xPoint, yPoint * 100);

/////////////////////////////////////////////////////////////////////////////////////

    double Tmethod_val_at_t_manypoints = Tmethod_fitFunction_60h_manypoints->Eval(xPoint);;
    double Tmethod_val_at_tminus_manypoints = Tmethod_fitFunction_60h_manypoints->Eval(xPoint-SGT_60h/2);;
    double Tmethod_val_at_tplus_manypoints = Tmethod_fitFunction_60h_manypoints->Eval(xPoint+SGT_60h/2);;

    double yPointNum_manypoints = Tmethod_val_at_t_manypoints*Tmethod_val_at_t_manypoints - Tmethod_val_at_tminus_manypoints*Tmethod_val_at_tplus_manypoints;
    double yPointDenom_manypoints = pow(Tmethod_val_at_t_manypoints,2);
    double yPoint_manypoints = yPointNum_manypoints/yPointDenom_manypoints;

    Tmethod_rateError_60h_manypoints->SetPoint(Tmethod_rateError_60h_manypoints->GetN(), xPoint, yPoint_manypoints * 100);

/////////////////////////////////////////////////////////////////////////////////////

    xPoint += histBinWidth/xPointDivisions;
  }

/////////////////////////////////////////////////////////////////////////////////////


  TGraph* Tmethod_rateError_60h_average = new TGraph();
  Tmethod_rateError_60h_average->SetLineColor(4);
  Tmethod_rateError_60h_average->SetMarkerColor(4);

  int numPoints = Tmethod_rateError_60h->GetN()/xPointDivisions;

  for (int i = 0; i < numPoints; ++i)
  {
    double yAverage = 0;
    double x,y;
    for (int j = 0; j < xPointDivisions; ++j){
      Tmethod_rateError_60h->GetPoint(int(i*xPointDivisions + j), x, y);
      yAverage += y;
    } 

    double xVal = (xmin + histBinWidth/2. + i*histBinWidth);

    Tmethod_rateError_60h_average->SetPoint(Tmethod_rateError_60h_average->GetN(), xVal, yAverage/xPointDivisions);
  }


/////////////////////////////////////////////////////////////////////////////////////


  nsTOus(Tmethod_rateError_60h, "Time [#mus]");
  Tmethod_rateError_60h->GetYaxis()->SetTitle("Pileup Rate Error [%]");

  nsTOus(Tmethod_rateError_60h_manypoints, "Time [#mus]");
  nsTOus(Tmethod_rateError_60h_average, "Time [#mus]");


  auto c_Tmethod_rateError_60h = new TCanvas("c_Tmethod_rateError_60h","c_Tmethod_rateError_60h",200,10,1000,800);
  Tmethod_rateError_60h->Draw("APL");
  Tmethod_rateError_60h_manypoints->Draw("SAMEPL");
  Tmethod_rateError_60h_average->Draw("SAMEPL");
  Tmethod_rateError_60h->GetXaxis()->SetRangeUser(30, 40);
  Tmethod_rateError_60h->GetYaxis()->SetRangeUser(-0.05, 0.03);


        auto myLegend = new TLegend(0.6,0.75,0.975,0.95);

        myLegend->AddEntry(Tmethod_rateError_60h, "Nominal fit - 10,000 points","pl");
        myLegend->AddEntry(Tmethod_rateError_60h_manypoints, "Fit with 100,000 points","pl");
        myLegend->AddEntry(Tmethod_rateError_60h_average, "Bin width average of nominal fit ","pl");
  
        myLegend->SetBorderSize(1);
        myLegend->SetFillStyle(1001);
        myLegend->Draw();


  c_Tmethod_rateError_60h->SaveAs("PileupRateError.png");

  Tmethod_rateError_60h->GetXaxis()->SetRangeUser(30, 70);
  c_Tmethod_rateError_60h->Modified();
  c_Tmethod_rateError_60h->SaveAs("PileupRateError_ZoomedOut.png");

  // Tmethod_rateError_60h->GetXaxis()->SetRangeUser(30, 70);
  // c_Tmethod_rateError_60h->Modified();
  // c_Tmethod_rateError_60h->SaveAs("PileupRateError_149p2SGT.png");


  // TH1F* Tmethod_rateError_60h_FFT = doFFT(Tmethod_rateError_60h_average);
  // auto tempcanv = new TCanvas("tempcanv","tempcanv",300,10,1000,800);
  // Tmethod_rateError_60h_FFT->Draw("HIST");
  // drawLineOnCanv(0.23, tempcanv);
  // drawLineOnCanv(2*0.23, tempcanv);


/////////////////////////////////////////////////////////////////////////////////////

    // outputFile->Write();
    // delete outputFile;

    return 1;

}
