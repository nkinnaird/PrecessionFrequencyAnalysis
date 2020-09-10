// this macro compares two different average fit start scans made with AverageFitStartScan.C
// was originally written to compare average fit start scans with and without the ad hoc correction applied

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
#include <TRandom3.h>
#include <sstream>
#include <TVirtualFFT.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TPaveStats.h>

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;


double fitFunc(double* x, double* par)
{
  double time = x[0];
  double funcValue = par[0] * exp(-time/par[1]) * cos(par[2]*time + par[3]);
  return funcValue;
}


void drawLines(TCanvas* inputCanv)
{
  inputCanv->Update();

  double canvUxmin = inputCanv->GetUxmin();
  double canvUxmax = inputCanv->GetUxmax();
  double canvUymin = inputCanv->GetUymin();
  double canvUymax = inputCanv->GetUymax();

  double startingUx = 30287.6/1000; // use units of us
  // int lineNum = 0;

  for (double Ux = startingUx; Ux < canvUxmax; Ux = Ux + g2Period/1000.)
  {
    TLine *thisLine = new TLine(Ux, canvUymin, Ux, canvUymax);
    thisLine->SetLineStyle(2);
    thisLine->SetLineWidth(2);
    // if(lineNum % 2 == 0) thisLine->SetLineColor(1);
    // else thisLine->SetLineColor(2);
    thisLine->SetLineColor(4);
    thisLine->Draw();
    // lineNum++;
  }
}

int AverageFitStartScanComparison()
{
  // gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  // Endgame short range
  // string firstFile_string = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/Endgame/RandSeeds/FitStartScans/averageFitStartScan.root"; // typically uncorrected file
  // string secondFile_string = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/Endgame/RandSeeds/withAdHoc/FitStartScans/averageFitStartScan.root"; // typically corrected file

  // Endgame longer range
  string firstFile_string = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/Endgame/RandSeeds/FitStartScans-LongerRange/averageFitStartScan.root"; // typically uncorrected file
  string secondFile_string = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/Endgame/RandSeeds/withAdHoc/FitStartScans-LongerRange/averageFitStartScan.root"; // typically corrected file

  TFile *firstFile = TFile::Open(firstFile_string.c_str());
  TFile *secondFile = TFile::Open(secondFile_string.c_str());
   if (firstFile == 0 || secondFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

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

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


  TGraph* TMethod_R_Average_firstFile = (TGraph*) firstFile->Get("TMethod_R_Average");
  TGraph* TMethod_R_Average_secondFile = (TGraph*) secondFile->Get("TMethod_R_Average");

  TGraph* TR_Diff_firstFile = (TGraph*) firstFile->Get("TRAverageDifference");

  TGraph* RMethod_R_Average_firstFile = (TGraph*) firstFile->Get("RMethod_R_Average");
  TGraph* RMethod_R_Average_secondFile = (TGraph*) secondFile->Get("RMethod_R_Average");

  TGraph* TR_Diff_secondFile = (TGraph*) secondFile->Get("TRAverageDifference");


  if(TMethod_R_Average_firstFile->GetN() != TMethod_R_Average_secondFile->GetN()){
    cout << "Files have different numbers of points - exiting." << endl;
    exit(0);
  }


  TGraph* TMethod_R_Differences = new TGraph();
  TGraph* RMethod_R_Differences = new TGraph();

  TGraph* TDiffs_minus_RDiffs = new TGraph();

  TGraph* Sum_of_Diffs = new TGraph(); // subtract first file TR diff from the T-T diff calculated here


  for (int pointNo = 0; pointNo < TMethod_R_Average_firstFile->GetN(); ++pointNo)
  {
    double T_x_first, T_y_first, T_x_second, T_y_second, R_x_first, R_y_first, R_x_second, R_y_second;

    TMethod_R_Average_firstFile->GetPoint(pointNo, T_x_first, T_y_first);
    TMethod_R_Average_secondFile->GetPoint(pointNo, T_x_second, T_y_second);

    RMethod_R_Average_firstFile->GetPoint(pointNo, R_x_first, R_y_first);
    RMethod_R_Average_secondFile->GetPoint(pointNo, R_x_second, R_y_second);

    // cout << "values check : " << T_y_first << " " << T_y_second << endl;

    if(T_x_first != T_x_second){
      cout << "Files have different fit start values - exiting." << endl;
      exit(0);    
    }

    double Tdiff = 1000. * (T_y_second - T_y_first);
    double Rdiff = 1000. * (R_y_second - R_y_first);

    TMethod_R_Differences->SetPoint(pointNo, T_x_first, Tdiff);
    RMethod_R_Differences->SetPoint(pointNo, R_x_first, Rdiff);

    TDiffs_minus_RDiffs->SetPoint(pointNo, T_x_first, Tdiff - Rdiff);


    double TRdiff_x_first, TRdiff_y_first;

    TR_Diff_firstFile->GetPoint(pointNo, TRdiff_x_first, TRdiff_y_first);

    Sum_of_Diffs->SetPoint(pointNo, TRdiff_x_first, TRdiff_y_first + Tdiff);

  }


  TMethod_R_Differences->SetName("TMethod_R_Differences");
  TMethod_R_Differences->GetXaxis()->SetTitle(TMethod_R_Average_firstFile->GetXaxis()->GetTitle());
  TMethod_R_Differences->GetYaxis()->SetTitle("#DeltaR (with-without) #LTT-Method#GT (ppb)");
  TMethod_R_Differences->GetYaxis()->SetTitleOffset(2);

  auto Tdiff_canv = new TCanvas("Tdiff_canv","Tdiff_canv",200,10,1200,800);
  TMethod_R_Differences->Draw("AL");

  Tdiff_canv->Update();

  TLine *zeroLine = new TLine(Tdiff_canv->GetUxmin(), 0, Tdiff_canv->GetUxmax(), 0);
  zeroLine->SetLineColor(1);
  zeroLine->SetLineStyle(2);
  zeroLine->SetLineWidth(3);
  zeroLine->Draw();

  drawLines(Tdiff_canv);

  // Tdiff_canv->SaveAs("TAverage_Diff_wwo_adHoc.png");

/////////////////////////////////////////////////////////////////////////////////////

  RMethod_R_Differences->SetName("RMethod_R_Differences");
  RMethod_R_Differences->GetXaxis()->SetTitle(RMethod_R_Average_firstFile->GetXaxis()->GetTitle());
  RMethod_R_Differences->GetYaxis()->SetTitle("#DeltaR (with-without) #LTR-Method#GT (ppb)");
  RMethod_R_Differences->GetYaxis()->SetTitleOffset(2);

  auto Rdiff_canv = new TCanvas("Rdiff_canv","Rdiff_canv",200,10,1200,800);
  RMethod_R_Differences->Draw("AL");

  Rdiff_canv->Update();
  // zeroLine->Draw();
  drawLines(Rdiff_canv);

  // Rdiff_canv->SaveAs("RAverage_Diff_wwo_adHoc.png");

/////////////////////////////////////////////////////////////////////////////////////

  TGraph* T_diff_clone = (TGraph*) TMethod_R_Differences->Clone();
  TGraph* R_diff_clone = (TGraph*) RMethod_R_Differences->Clone();


  T_diff_clone->GetYaxis()->SetTitle("#DeltaR (with-without) (ppb)");
  R_diff_clone->SetLineColor(2);

  auto TandR_diff_canv = new TCanvas("TandR_diff_canv","TandR_diff_canv",200,10,1200,800);
  T_diff_clone->Draw("AL");
  R_diff_clone->Draw("LSAME");
  zeroLine->Draw();
  drawLines(TandR_diff_canv);

        auto diffLegend = new TLegend(0.7,0.675,0.975,0.875);

          diffLegend->AddEntry(T_diff_clone, "T Method","l");
          diffLegend->AddEntry(R_diff_clone, "R Method","l");
  
        diffLegend->SetBorderSize(1);
        diffLegend->SetFillStyle(1001);
        diffLegend->Draw();

  // TandR_diff_canv->SaveAs("TandR_Average_Diff_wwo_adHoc.png");

/////////////////////////////////////////////////////////////////////////////////////

  gStyle->SetOptFit(1111);

  double fitRangeLow = 50;
  double fitRangeHigh = 95;

  TF1* diffFunction = new TF1("diffFunction", fitFunc, fitRangeLow, fitRangeHigh, 4);
  diffFunction->SetLineColor(2);
  diffFunction->SetNpx(10000);

  diffFunction->SetParameter(0, 100);
  diffFunction->SetParameter(1, 64.4);
  diffFunction->SetParameter(2, 1000.*defaultWa);
  diffFunction->SetParameter(3, 0);

  TDiffs_minus_RDiffs->Fit(diffFunction, "R");

  auto deltaT_minus_deltaR_canv = new TCanvas("deltaT_minus_deltaR_canv","deltaT_minus_deltaR_canv",200,10,1200,800);
  TDiffs_minus_RDiffs->Draw("AL");
  zeroLine->Draw();
  drawLines(deltaT_minus_deltaR_canv);

  deltaT_minus_deltaR_canv->Update();

  TPaveStats* statsBox = (TPaveStats*) TDiffs_minus_RDiffs->GetListOfFunctions()->FindObject("stats");
  statsBox->SetBorderSize(1);
  statsBox->Draw("SAME");


/////////////////////////////////////////////////////////////////////////////////////

  TGraph* playGraph = (TGraph*) TR_Diff_firstFile->Clone();
  // TGraph* playGraph = (TGraph*) TR_Diff_secondFile->Clone();

  TF1* playFunc = new TF1("playFunc", "[0]*cos([1]*x+[2])", 30, 40);
  // TF1* playFunc = new TF1("playFunc", "[0]*cos([1]*x+[2])", 60, 90);
  playFunc->SetLineColor(2);
  playFunc->SetNpx(10000);

  playFunc->SetParameter(0, 100);
  playFunc->SetParameter(1, 1000.*defaultWa);
  // playFunc->FixParameter(1, 1000.*defaultWa);
  playFunc->SetParameter(2, 0);

  playGraph->Fit(playFunc, "R");


          // TF1* playFunc = new TF1("playFunc", "[0]*cos([1]*x+[2])", 30, 40);
          // TF1* playFuncccc = new TF1("playFuncccc", "[0]*cos([1]*x+[2])", 60, 90);
          // playFuncccc->SetLineColor(2);
          // playFuncccc->SetNpx(10000);

          // playFuncccc->SetParameter(0, 100);
          // playFuncccc->SetParameter(1, 1000.*defaultWa);
          // playFuncccc->SetParameter(2, 0);

          // playGraph->Fit(playFuncccc, "R+");



  for (int pointNo = 0; pointNo < playGraph->GetN(); ++pointNo)
  {
    double x, y; 
    playGraph->GetPoint(pointNo, x, y);
    playGraph->SetPoint(pointNo, x, y - playFunc->Eval(x));
  }

  playGraph->GetFunction("playFunc")->SetBit(TF1::kNotDraw);

  auto playCanv = new TCanvas("playCanv","playCanv",200,10,1200,800);
  playGraph->Draw("AL");
  zeroLine->Draw();
  drawLines(playCanv);

  TPaveStats* statsBox_2 = (TPaveStats*) playGraph->GetListOfFunctions()->FindObject("stats");
  statsBox_2->SetBorderSize(1);
  statsBox_2->Draw("SAME");


  // do FFT of difference plot

  // TH1F* fft_of_TRdiff = doFFT(playGraph);
  // auto fftCanv = new TCanvas("fftCanv","fftCanv",200,10,1200,800);
  // fft_of_TRdiff->Draw("HIST");

  // drawLineOnCanv(0.23, fftCanv);
  // drawLineOnCanv(0.37-0.23, fftCanv);

/*

  TGraph* playGraph2 = (TGraph*) TR_Diff_secondFile->Clone();

  TF1* playFunc2 = new TF1("playFunc2", "[0]*cos([1]*x+[2])", 30, 100);
  playFunc2->SetLineColor(2);
  playFunc2->SetNpx(10000);

  playFunc2->SetParameter(0, 100);
  playFunc2->SetParameter(1, 1000.*defaultWa);
  playFunc2->SetParameter(2, 0);

  playGraph2->Fit(playFunc2, "R");

  for (int pointNo = 0; pointNo < playGraph2->GetN(); ++pointNo)
  {
    double x, y; 
    playGraph2->GetPoint(pointNo, x, y);
    playGraph2->SetPoint(pointNo, x, y);
  }

  playGraph2->GetFunction("playFunc2")->SetBit(TF1::kNotDraw);


  auto playCanv_2 = new TCanvas("playCanv_2","playCanv_2",200,10,1200,800);
  // TR_Diff_firstFile->Draw("AL");
  // playGraph2->Draw("LSAME");
  playGraph2->Draw("AL");
  zeroLine->Draw();
  drawLines(playCanv_2);

  TPaveStats* statsBox_3 = (TPaveStats*) playGraph2->GetListOfFunctions()->FindObject("stats");
  statsBox_3->SetBorderSize(1);
  statsBox_3->Draw("SAME");




  auto playCanv_3 = new TCanvas("playCanv_3","playCanv_3",200,10,1200,800);
  playGraph->Draw("AL");
  playGraph2->SetLineColor(2);
  playGraph2->Draw("LSAME");
  drawLines(playCanv_3);

  TGraph* playGraph3 = (TGraph*) playGraph->Clone();

  for (int pointNo = 0; pointNo < playGraph3->GetN(); ++pointNo)
  {
    double x, y; 
    playGraph3->GetPoint(pointNo, x, y);

    double x2, y2;
    playGraph2->GetPoint(pointNo, x2, y2);
    
    playGraph3->SetPoint(pointNo, x, y-y2);
  }

  auto playCanv_4 = new TCanvas("playCanv_4","playCanv_4",200,10,1200,800);
  playGraph3->Draw("AL");
  drawLines(playCanv_4);


*/
/////////////////////////////////////////////////////////////////////////////////////

/*
  Sum_of_Diffs->SetName("Sum_of_Diffs");
  Sum_of_Diffs->GetXaxis()->SetTitle(RMethod_R_Average_firstFile->GetXaxis()->GetTitle());
  Sum_of_Diffs->GetYaxis()->SetTitle("#DeltaR (T_{w}-T_{wo} plus T_{wo}-R_{wo}) (ppb)");
  Sum_of_Diffs->GetYaxis()->SetTitleOffset(2);

  auto sumDiff_canv = new TCanvas("sumDiff_canv","sumDiff_canv",200,10,1200,800);
  Sum_of_Diffs->Draw("AL");

  sumDiff_canv->Update();
  zeroLine->Draw();
  drawLines(sumDiff_canv);

  sumDiff_canv->SaveAs("ConsistencyCheck_TRAverageDiff.png");
*/

/////////////////////////////////////////////////////////////////////////////////////

  return 1;
}
