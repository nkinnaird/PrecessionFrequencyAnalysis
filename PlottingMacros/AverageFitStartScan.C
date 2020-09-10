// this macro plots the average fit start scan for many fit start scan fits to different random seeds
// do fits something like this in parallel (on gm2ucl):  #for i in `seq 0 99`; do echo $i; done | xargs -i --max-procs=25 bash -c 'mkdir Iter{} && cd Iter{}; root -l -b -q $RP/ratioMacroParallel.C+\(\"../../addedHists-Endgame-RandSeeds.root\",\"Iter{}\"\)'
// then create and pass in root file of fitted functions with StripFunctions.C

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

void lombScargle(TGraph* inGraph, double windowTime){
  ////////////////////////////////////////////////////////////////
  // Lomb Scargle (of sorts)
  ////////////////////////////////////////////////////////////////

  // Input params
  double omega = defaultWa * 1e3;
  double window = windowTime; // us (windows to scan in)

  // Output graph
  TGraph* powerVsTime = new TGraph();

  // Loop over all the time slices
  double t, f;
  int iStartPt = 0;

  while(iStartPt <= inGraph->GetN()-1){

    inGraph->GetPoint(iStartPt,t,f);
    double minTime = t;
    double maxTime = t + window;
    int iEndPt = inGraph->GetN() - 1;
    for(int iPt = iStartPt; iPt < inGraph->GetN(); iPt++){
      inGraph->GetPoint(iPt,t,f);
      if(t > maxTime) {
  iEndPt = iPt - 1;
  break;
      }
    }
    inGraph->GetPoint(iEndPt, maxTime, f);
    int N = iEndPt - iStartPt + 1;    

    // Calculate mean
    double sum = 0;
    for(int n = iStartPt; n < iEndPt; n++){
      inGraph->GetPoint(n, t, f);
      sum += f;
    }
    double mean = sum / N;
    
    // Calculate sum of sin and cosine terms to get tau
    double sumSin(0), sumCos(0);
    for(int n = iStartPt; n < iEndPt; n++){
      inGraph->GetPoint(n, t, f);
      sumSin += sin(TMath::TwoPi()*t);
      sumCos += cos(TMath::TwoPi()*t);
    }
    double tau = 1./(2*omega)*atan(sumSin/sumCos);
    
    // Calculate all the ingredients for the power
    double cosNum(0), cosDen(0), sinNum(0), sinDen(0);
    for(int n = iStartPt; n < iEndPt; n++){
      inGraph->GetPoint(n, t, f);
      cosNum += (f - mean)*cos(omega*(t - tau));
      cosDen += pow(cos(omega*(t - tau)),2);
      sinNum += (f - mean)*sin(omega*(t - tau));
      sinDen += pow(sin(omega*(t - tau)),2);
    }

    // Calculate sqrt of power (this is sort of amplitude)
    double power = sqrt(cosNum*cosNum/cosDen + sinNum*sinNum/sinDen);
    powerVsTime->SetPoint(powerVsTime->GetN(), 0.5*(minTime+maxTime), power);
    
    // Update start point
    iStartPt = iEndPt + 1;
  }

  TCanvas* c2 = new TCanvas();
  powerVsTime->SetMarkerStyle(20);
  // powerVsTime->SetTitle(";Time [us];Amplitude [arb. units]");
  powerVsTime->SetTitle(Form(";Time [#mus];Amplitude [arb. units], %0.4g #mus Window",window));
  powerVsTime->Draw("AP");

}

void drawLines(TCanvas* inputCanv)
{
  inputCanv->Update();

  double canvUxmin = inputCanv->GetUxmin();
  double canvUxmax = inputCanv->GetUxmax();
  double canvUymin = inputCanv->GetUymin();
  double canvUymax = inputCanv->GetUymax();

  double startingUx = 30287.6/1000; // use units of us

  for (double Ux = startingUx; Ux < canvUxmax; Ux = Ux + g2Period/1000.)
  {
    if(Ux > canvUxmin){
      TLine *thisLine = new TLine(Ux, canvUymin, Ux, canvUymax);
      thisLine->SetLineStyle(2);
      thisLine->SetLineWidth(2);
      thisLine->SetLineColor(4);
      thisLine->Draw();
    }
  }
}


int AverageFitStartScan(string functionsFile)
{
  // gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  TFile* outputFile = new TFile("averageFitStartScan.root","RECREATE");

  TFile *inputFile = TFile::Open(functionsFile.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  TIter next(inputFile->GetListOfKeys());
  TKey* key;

  bool containsTag = false;
  key = (TKey*)next(); // to skip datasetTag if it exists
  string firstKeyName = string(key->GetName());
  if(firstKeyName.compare("datasetTag") != 0){
    cout << "datasetTag missing, name of first key: " << firstKeyName << endl;
    next.Reset();
  }
  else{
    cout << "datasetTag contained within file, making sure to skip it: " << firstKeyName << endl;
    containsTag = true;
  }

  uint totalFiles = 0;
  while((key = (TKey*)next())) totalFiles++;

  cout << "Total files: " << totalFiles << endl;

  uint filesToRunOver = totalFiles;
  // uint filesToRunOver = 5;

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

  vector<float> fSs;

  vector<vector<double>> TMethod_R_Values; // first vector is over files, second vector is over fit start times
  vector<vector<double>> RMethod_R_Values; // first vector is over files, second vector is over fit start times


// get number of passes and the fit start times from the first file
      TNtuple *firstTuple = (TNtuple*)inputFile->Get(Form("File0/topDir/FitPasses/FitPass0/FitConditions0"));
      float totalFitPasses;
      firstTuple->SetBranchAddress("totalPasses", &totalFitPasses);
      firstTuple->GetEntry(0);
        for (int fitPass = 0; fitPass < totalFitPasses; ++fitPass)
        {
          float fS;
          TNtuple *ntuple = (TNtuple*)inputFile->Get(Form("File0/topDir/FitPasses/FitPass%d/FitConditions%d", fitPass, fitPass));
          ntuple->SetBranchAddress("fitStartTime", &fS);
          ntuple->GetEntry(0);
          fSs.push_back(fS);
        }



// loop through files and get R values
  uint fileNum = 0;
  next.Reset();
  if(containsTag) (key = (TKey*)next()); // to skip datasetTag 

  while((key = (TKey*)next()) && fileNum < filesToRunOver)
  {
    cout << "File num: " << fileNum << endl;
    string fileName = key->GetName();

    vector<double> Tmethod_file_R_Values;
    vector<double> Rmethod_file_R_Values;    

      for (int fitPass = 0; fitPass < totalFitPasses; ++fitPass)
      {
        string TmethodFunc_path = Form("/topDir/FitPasses/FitPass%d/addedDir/TMethod/TmethodFitFunc", fitPass);
        string RmethodFunc_path = Form("/topDir/FitPasses/FitPass%d/addedDir/FullRatio/fullRatioFitFunc", fitPass);

        double TMethod_R = ((TF1*) inputFile->Get((fileName + TmethodFunc_path).c_str()))->GetParameter(3);
        double RMethod_R = ((TF1*) inputFile->Get((fileName + RmethodFunc_path).c_str()))->GetParameter(1);

        Tmethod_file_R_Values.push_back(TMethod_R);
        Rmethod_file_R_Values.push_back(RMethod_R);
      }

      TMethod_R_Values.push_back(Tmethod_file_R_Values);
      RMethod_R_Values.push_back(Rmethod_file_R_Values);

      fileNum++;
  }


/////////////////////////////////////////////////////////////////////////////////////

  vector<double> TMethod_R_Averages(totalFitPasses, 0);
  vector<double> RMethod_R_Averages(totalFitPasses, 0);

    for (uint fitPass = 0; fitPass < totalFitPasses; ++fitPass)
    {
      for (uint fileNum = 0; fileNum < filesToRunOver; fileNum++)
      {
        TMethod_R_Averages.at(fitPass) += TMethod_R_Values.at(fileNum).at(fitPass);
        RMethod_R_Averages.at(fitPass) += RMethod_R_Values.at(fileNum).at(fitPass);
      }

      TMethod_R_Averages.at(fitPass) = TMethod_R_Averages.at(fitPass)/filesToRunOver;
      RMethod_R_Averages.at(fitPass) = RMethod_R_Averages.at(fitPass)/filesToRunOver;
    }


/////////////////////////////////////////////////////////////////////////////////////

    TGraph* Tmethod_Rval_average = new TGraph();
    TGraph* Rmethod_Rval_average = new TGraph();

    Tmethod_Rval_average->SetMarkerColor(1);
    Tmethod_Rval_average->SetLineColor(1);
    Rmethod_Rval_average->SetMarkerColor(2);
    Rmethod_Rval_average->SetLineColor(2);


    TGraph* T_R_Diff_average = new TGraph();

    T_R_Diff_average->SetMarkerSize(0.75);
    T_R_Diff_average->SetLineColor(1);
    T_R_Diff_average->SetMarkerColor(1);


/////////////////////////////////////////////////////////////////////////////////////

int pointNo = 0;
double T_Rval_first = 0, R_Rval_first = 0;

    for (float fS : fSs)
    {
      double TMethod_R = TMethod_R_Averages.at(pointNo);
      double FullRatio_R = RMethod_R_Averages.at(pointNo);

      if(pointNo == 0){
        T_Rval_first = TMethod_R;
        R_Rval_first = FullRatio_R;
      }

      // Tmethod_Rval_average->SetPoint(pointNo, fS, 1000 * (TMethod_R - T_Rval_first));
      // Rmethod_Rval_average->SetPoint(pointNo, fS, 1000 * (FullRatio_R - T_Rval_first));

      Tmethod_Rval_average->SetPoint(pointNo, fS, TMethod_R);
      Rmethod_Rval_average->SetPoint(pointNo, fS, FullRatio_R);

      T_R_Diff_average->SetPoint(pointNo, fS, 1000. * (TMethod_R - FullRatio_R));

      pointNo++;
    } // end loop over point numbers


/////////////////////////////////////////////////////////////////////////////////////

  string xAxisTitle = "Fit Start Time [#mus]";

        nsTOus(Tmethod_Rval_average, xAxisTitle);
        nsTOus(Rmethod_Rval_average, xAxisTitle);
        nsTOus(T_R_Diff_average, xAxisTitle);


      Tmethod_Rval_average->SetName("TMethod_R_Average");
      Tmethod_Rval_average->GetXaxis()->SetTitle(xAxisTitle.c_str());
      Tmethod_Rval_average->GetYaxis()->SetTitle("#LTR#GT [ppb]");
      Tmethod_Rval_average->GetYaxis()->SetTitleOffset(2);
      // Tmethod_Rval_average->GetYaxis()->SetRangeUser(-500, 100);
      // Tmethod_Rval_average->GetYaxis()->SetRangeUser(-800, 500);

      Rmethod_Rval_average->SetName("RMethod_R_Average");
      Rmethod_Rval_average->GetXaxis()->SetTitle(xAxisTitle.c_str());
      Rmethod_Rval_average->GetYaxis()->SetTitle("#LTR#GT [ppb]");
      Rmethod_Rval_average->GetYaxis()->SetTitleOffset(2);

      auto TR_val_average_canv = new TCanvas("TR_val_average_canv","TR_val_average_canv",200,10,1200,800);

      Tmethod_Rval_average->Draw("AL");
      Rmethod_Rval_average->Draw("SAMEL");

      TR_val_average_canv->Update();

      TLine *zeroLine = new TLine(TR_val_average_canv->GetUxmin(), 0, TR_val_average_canv->GetUxmax(), 0);
      zeroLine->SetLineColor(1);
      zeroLine->SetLineStyle(2);
      zeroLine->SetLineWidth(3);
      zeroLine->Draw();

      drawLines(TR_val_average_canv);

        auto legend_val_average = new TLegend(0.7,0.675,0.975,0.875);

          legend_val_average->AddEntry(Tmethod_Rval_average, "T Method Average","l");
          legend_val_average->AddEntry(Rmethod_Rval_average, "R Method Average","l");
  
        legend_val_average->SetBorderSize(1);
        legend_val_average->SetFillStyle(1001);
        legend_val_average->Draw();

      TR_val_average_canv->SaveAs("TandR_Average.png");

      outputFile->cd();
      Tmethod_Rval_average->Write();
      Rmethod_Rval_average->Write();
      TR_val_average_canv->Write();



/////////////////////////////////////////////////////////////////////////////////////

      T_R_Diff_average->SetName("TRAverageDifference");
      T_R_Diff_average->GetXaxis()->SetTitle(xAxisTitle.c_str());
      T_R_Diff_average->GetYaxis()->SetTitle("#Delta R (#LTT#GT-#LTR#GT) [ppb]");
      T_R_Diff_average->GetYaxis()->SetTitleOffset(2);
      T_R_Diff_average->GetYaxis()->SetRangeUser(-300,300);

      auto TR_diff_compCanvas = new TCanvas("TR_diff_compCanvas","TR_diff_compCanvas",200,10,1200,800);

      T_R_Diff_average->Draw("AL");
      TR_diff_compCanvas->Update();
      zeroLine->Draw();
      drawLines(TR_diff_compCanvas);

      TR_diff_compCanvas->SaveAs("TR_Diff_Average.png");
      TR_diff_compCanvas->SaveAs("TR_Diff_Average.C");

      T_R_Diff_average->Write();
      TR_diff_compCanvas->Write();

      TH1F* fft_of_diff = doFFT(T_R_Diff_average);
      fft_of_diff->Write("TR_Diff_Average_FFT");

      auto fft_diff_canv = new TCanvas("fft_diff_canv","fft_diff_canv",200,10,1200,800);
      fft_of_diff->Draw("HIST");

      drawLineOnCanv(0.23, fft_diff_canv);

      fft_diff_canv->SaveAs("TR_Diff_Average_FFT.png");

/////////////////////////////////////////////////////////////////////////////////////

      lombScargle(T_R_Diff_average, 4);
      lombScargle(T_R_Diff_average, 1*g2Period/1000);
      lombScargle(T_R_Diff_average, 6);
      lombScargle(T_R_Diff_average, 2*g2Period/1000);
      lombScargle(T_R_Diff_average, 10);

/////////////////////////////////////////////////////////////////////////////////////

      auto individualSeedsDir = outputFile->mkdir("IndividualSeeds");
      individualSeedsDir->cd();

      // make some single file fit start scans for comparison
      bool makeIndividualPlots = false;

      for (uint fileNum = 0; fileNum < filesToRunOver; ++fileNum)
      {
        TGraph* Tmethod_Rval_single = new TGraph();
        TGraph* Rmethod_Rval_single = new TGraph();

        Tmethod_Rval_single->SetMarkerColor(1);
        Tmethod_Rval_single->SetLineColor(1);
        Rmethod_Rval_single->SetMarkerColor(2);
        Rmethod_Rval_single->SetLineColor(2);

        for (uint fitPass = 0; fitPass < totalFitPasses; ++fitPass)
        {
          Tmethod_Rval_single->SetPoint(fitPass, fSs.at(fitPass), TMethod_R_Values.at(fileNum).at(fitPass));
          Rmethod_Rval_single->SetPoint(fitPass, fSs.at(fitPass), RMethod_R_Values.at(fileNum).at(fitPass));
        }

        nsTOus(Tmethod_Rval_single, xAxisTitle);
        nsTOus(Rmethod_Rval_single, xAxisTitle);

        Tmethod_Rval_single->SetName(Form("T_%i",fileNum));
        Tmethod_Rval_single->GetYaxis()->SetTitle("#LTR#GT [ppm]");
        Tmethod_Rval_single->GetYaxis()->SetTitleOffset(2);

        Rmethod_Rval_single->SetName(Form("R_%i",fileNum));
        Rmethod_Rval_single->GetYaxis()->SetTitle("#LTR#GT [ppm]");
        Rmethod_Rval_single->GetYaxis()->SetTitleOffset(2);

        Tmethod_Rval_single->Write();
        Rmethod_Rval_single->Write();


        if(fileNum < 4 && makeIndividualPlots)
        {
          auto TR_val_canv = new TCanvas(Form("TR_val_canv_%i",fileNum),"TR_val_canv",500,10,1200,800);

          Tmethod_Rval_single->Draw("AL");
          Rmethod_Rval_single->Draw("SAMEL");

          TR_val_canv->Update();

          drawLines(TR_val_canv);

          auto legend_val = new TLegend(0.7,0.675,0.975,0.875);

            legend_val->AddEntry(Tmethod_Rval_single, Form("T Method_%i",fileNum),"l");
            legend_val->AddEntry(Rmethod_Rval_single, Form("R Method_%i",fileNum),"l");
    
          legend_val->SetBorderSize(1);
          legend_val->SetFillStyle(1001);
          legend_val->Draw();

          TR_val_canv->SaveAs(Form("TandR_%i.png",fileNum));

          TR_val_canv->Write();
        }
      }


/////////////////////////////////////////////////////////////////////////////////////

  return 1;
}
