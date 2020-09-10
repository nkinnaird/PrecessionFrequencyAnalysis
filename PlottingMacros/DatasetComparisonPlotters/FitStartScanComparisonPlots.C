// 3-31-20: Macro to take fit start scan plots from different datasets, or at least the R parameter fit start scans, and plot them against each other.

#include <iostream>
#include <fstream>
#include <string>
#include <functional>
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
#include <TVectorD.h>
#include <TPaveText.h>

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

bool saveImages = true;

int FitStartScanComparisonPlots()
{

  gStyle->SetOptStat(000000);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1);
  gStyle->SetLineColor(1);
  gStyle->SetPadRightMargin(.05);
  gStyle->SetPadLeftMargin(.15);


  vector<string> dataset_files = {"/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/60h/SingleIteration/FitStartScans/fitStartScan-fixed-cbo-tau.root",
                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/HighKick/SingleIteration/FitStartScans/fitStartScan-fixed-cbo-phases-amps.root",
                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/9d/SingleIteration/FitStartScans/fitStartScan-fixed-cbo-phases-amps.root",
                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/Endgame/SingleIteration/FitStartScans/fitStartScan-fixed-cbo-A-A.root"};


  // auto combined_canv = new TCanvas("combined_canv", "combined_canv", 800, 10, 600, 400);
  auto combined_canv = new TCanvas("combined_canv", "combined_canv", 10, 10, 1200, 800);

    auto legend = new TLegend(0.8,0.8,.975,0.95);
    legend->SetBorderSize(1);
    legend->SetFillStyle(1001);


  for (uint fileNum = 0; fileNum < dataset_files.size(); ++fileNum)
  {
    TFile* inputFile = TFile::Open(dataset_files.at(fileNum).c_str());
    if(inputFile == 0){
      printf("Error: cannot open file\n");
      return 0;
    }

   TCanvas* individual_canv = (TCanvas*) inputFile->Get("topDir/Added/FullRatio/FullRatio_R_FS_Canv");

      TGraphErrors* R_graph = (TGraphErrors*) ((TGraphErrors*) individual_canv->GetListOfPrimitives()->FindObject("FullRatio_R_Vs_FS")->Clone());
      TGraph* K1 = (TGraph*) ((TGraph*) individual_canv->GetListOfPrimitives()->FindObject("KawallBandOne")->Clone());
      TGraph* K2 = (TGraph*) ((TGraph*) individual_canv->GetListOfPrimitives()->FindObject("KawallBandTwo")->Clone());

      double firstTime, firstVal;
      R_graph->GetPoint(0, firstTime, firstVal);

      for (int i = 0; i < R_graph->GetN(); ++i)
      {
        double time, val;
        R_graph->GetPoint(i, time, val);
        R_graph->SetPoint(i, time, val - firstVal);

        double K1_time, K1_val, K2_time, K2_val;
        K1->GetPoint(i, K1_time, K1_val);
        K2->GetPoint(i, K2_time, K2_val);

        K1->SetPoint(i, K1_time, K1_val - firstVal);
        K2->SetPoint(i, K2_time, K2_val - firstVal);
      }

    combined_canv->cd();

    R_graph->SetLineColor(dataset_colors.at(fileNum));
    R_graph->SetMarkerColor(dataset_colors.at(fileNum));

    R_graph->SetMarkerSize(0.8);
    R_graph->SetMarkerStyle(dataset_markers.at(fileNum));

    K1->SetLineColor(dataset_colors.at(fileNum));
    K1->SetMarkerColor(dataset_colors.at(fileNum));
    K2->SetLineColor(dataset_colors.at(fileNum));
    K2->SetMarkerColor(dataset_colors.at(fileNum));

    if(fileNum == 0){
      R_graph->GetYaxis()->SetTitle("#DeltaR (ppm)");
      R_graph->Draw("APZ");
    } 
    else R_graph->Draw("PZSAME");

      K1->Draw("SAME");
      K2->Draw("SAME");

    legend->AddEntry(R_graph, dataset_names.at(fileNum).c_str(), "pl");

    inputFile->Close();
  }

  legend->Draw("SAME");

  if(saveImages) combined_canv->SaveAs("fitStartScan_R_dataset_comparison.png");

/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
  
// code for when I was trying to put chi2 and R fit start scans vertically 

/*
  Float_t small = 1e-5;

  auto multiplePlots_canv = new TCanvas("multiplePlots_canv", "multiplePlots_canv", 10, 10, 800, 1000);
  // multiplePlots_canv->Divide(2,4,small,small);
  multiplePlots_canv->Divide(2,4,0,0);

  for (uint fileNum = 0; fileNum < dataset_files.size(); ++fileNum)
  {
    TFile* inputFile = TFile::Open(dataset_files.at(fileNum).c_str());

/////////////////////////////////////////////////////////////////////////////////////

    TCanvas* chi2_canv = (TCanvas*) inputFile->Get("topDir/Added/FullRatio/FullRatio_Chi2NDF_Vs_FS_canv");

    TGraphErrors* chi2_graph = (TGraphErrors*) ((TGraphErrors*) chi2_canv->GetListOfPrimitives()->FindObject("FullRatio_Chi2NDF_Vs_FS")->Clone());
    TGraph* K1_chi2 = (TGraph*) ((TGraph*) chi2_canv->GetListOfPrimitives()->FindObject("KawallBandOne")->Clone());
    TGraph* K2_chi2 = (TGraph*) ((TGraph*) chi2_canv->GetListOfPrimitives()->FindObject("KawallBandTwo")->Clone());
   
    int chi2_padNum = 2*fileNum+1;

    multiplePlots_canv->cd(chi2_padNum);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.01);

    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.01);
    
    // if(chi2_padNum == 1) gPad->SetTopMargin(0.1);
    // else if(chi2_padNum == 7) gPad->SetBottomMargin(0.1);

    chi2_graph->GetYaxis()->CenterTitle();
    chi2_graph->GetYaxis()->SetLabelFont(63);
    chi2_graph->GetYaxis()->SetLabelSize(14); // pixels
    chi2_graph->GetYaxis()->SetTitleFont(63);
    chi2_graph->GetYaxis()->SetTitleSize(14); // pixels
    chi2_graph->GetYaxis()->SetTitleOffset(4.5);

    chi2_graph->Draw("AP");
    K1_chi2->Draw("SAME");
    K2_chi2->Draw("SAME");

/////////////////////////////////////////////////////////////////////////////////////

    TCanvas* R_canv = (TCanvas*) inputFile->Get("topDir/Added/FullRatio/FullRatio_R_FS_Canv");

    TGraphErrors* R_graph = (TGraphErrors*) ((TGraphErrors*) R_canv->GetListOfPrimitives()->FindObject("FullRatio_R_Vs_FS")->Clone());
    TGraph* K1_R = (TGraph*) ((TGraph*) R_canv->GetListOfPrimitives()->FindObject("KawallBandOne")->Clone());
    TGraph* K2_R = (TGraph*) ((TGraph*) R_canv->GetListOfPrimitives()->FindObject("KawallBandTwo")->Clone());
   
    int R_padNum = 2*fileNum+2;

    multiplePlots_canv->cd(R_padNum);
    gPad->SetLeftMargin(0.2);
    gPad->SetRightMargin(0.01);

    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.01);

    // if(R_padNum == 2) gPad->SetTopMargin(0.1);
    // else if(R_padNum == 8) gPad->SetBottomMargin(0.1);

    R_graph->GetYaxis()->CenterTitle();
    R_graph->GetYaxis()->SetLabelFont(63);
    R_graph->GetYaxis()->SetLabelSize(14); // pixels
    R_graph->GetYaxis()->SetTitleFont(63);
    R_graph->GetYaxis()->SetTitleSize(14); // pixels
    R_graph->GetYaxis()->SetTitleOffset(6);

    R_graph->Draw("AP");
    K1_R->Draw("SAME");
    K2_R->Draw("SAME");


/////////////////////////////////////////////////////////////////////////////////////

    inputFile->Close();

  }
*/



  Float_t small = 1e-5;

  auto multiplePlots_canv = new TCanvas("multiplePlots_canv", "multiplePlots_canv", 10, 10, 1200, 600);
  // multiplePlots_canv->Divide(2,4,small,small);
  multiplePlots_canv->Divide(4,2,0,0);

  for (uint fileNum = 0; fileNum < dataset_files.size(); ++fileNum)
  {
    TFile* inputFile = TFile::Open(dataset_files.at(fileNum).c_str());

/////////////////////////////////////////////////////////////////////////////////////

    TCanvas* chi2_canv = (TCanvas*) inputFile->Get("topDir/Added/FullRatio/FullRatio_Chi2NDF_Vs_FS_canv");

    TGraphErrors* chi2_graph = (TGraphErrors*) ((TGraphErrors*) chi2_canv->GetListOfPrimitives()->FindObject("FullRatio_Chi2NDF_Vs_FS")->Clone());
    TGraph* K1_chi2 = (TGraph*) ((TGraph*) chi2_canv->GetListOfPrimitives()->FindObject("KawallBandOne")->Clone());
    TGraph* K2_chi2 = (TGraph*) ((TGraph*) chi2_canv->GetListOfPrimitives()->FindObject("KawallBandTwo")->Clone());
   
    int chi2_padNum = fileNum+1;

    multiplePlots_canv->cd(chi2_padNum);
    if(chi2_padNum == 1) gPad->SetLeftMargin(0.15);
    else gPad->SetLeftMargin(0.00);
      
    gPad->SetRightMargin(0.00);
    gPad->SetTopMargin(0.005);
    gPad->SetBottomMargin(0.01);
    
    // if(chi2_padNum == 1) gPad->SetTopMargin(0.1);
    // else if(chi2_padNum == 7) gPad->SetBottomMargin(0.1);

    chi2_graph->GetYaxis()->CenterTitle();
    chi2_graph->GetYaxis()->SetLabelFont(63);
    chi2_graph->GetYaxis()->SetLabelSize(14); // pixels
    chi2_graph->GetYaxis()->SetLabelOffset(0.01);
    chi2_graph->GetYaxis()->SetTitleFont(63);
    chi2_graph->GetYaxis()->SetTitleSize(14); // pixels
    chi2_graph->GetYaxis()->SetTitleOffset(2.8);
    chi2_graph->GetYaxis()->SetRangeUser(0.965, 1.055);


    chi2_graph->Draw("APZ");
    K1_chi2->Draw("SAME");
    K2_chi2->Draw("SAME");

    gPad->SetTickx(1);
    gPad->SetTicky(1);

    gPad->Update();

    TLine *line = new TLine(gPad->GetUxmin(), 1, gPad->GetUxmax(), 1);
    line->SetLineColor(2);
    line->SetLineWidth(4);
    line->SetLineStyle(2);
    line->Draw("SAME");

/////////////////////////////////////////////////////////////////////////////////////

    // TPaveText* datasetTitle = new TPaveText(GetNDCX(40),GetNDCY(1.06),GetNDCX(80),GetNDCY(1.2), "NDC");
    // datasetTitle->AddText(dataset_names.at(fileNum).c_str());
    // datasetTitle->SetBorderSize(0);
    // datasetTitle->SetTextColor(dataset_colors.at(fileNum));
    // datasetTitle->SetFillStyle(0);
    // datasetTitle->SetTextFont(63);
    // datasetTitle->SetTextSize(30);
    // datasetTitle->Draw("SAME");

/////////////////////////////////////////////////////////////////////////////////////

    TCanvas* R_canv = (TCanvas*) inputFile->Get("topDir/Added/FullRatio/FullRatio_R_FS_Canv");

    TGraphErrors* R_graph = (TGraphErrors*) ((TGraphErrors*) R_canv->GetListOfPrimitives()->FindObject("FullRatio_R_Vs_FS")->Clone());
    TGraph* K1_R = (TGraph*) ((TGraph*) R_canv->GetListOfPrimitives()->FindObject("KawallBandOne")->Clone());
    TGraph* K2_R = (TGraph*) ((TGraph*) R_canv->GetListOfPrimitives()->FindObject("KawallBandTwo")->Clone());
   
    int R_padNum = fileNum+5;

    multiplePlots_canv->cd(R_padNum);
    if(R_padNum == 5) gPad->SetLeftMargin(0.15);
    else gPad->SetLeftMargin(0.00);
    
    gPad->SetRightMargin(0.00);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.1);

    // if(R_padNum == 2) gPad->SetTopMargin(0.1);
    // else if(R_padNum == 8) gPad->SetBottomMargin(0.1);

    R_graph->GetYaxis()->CenterTitle();
    R_graph->GetYaxis()->SetLabelFont(63);
    R_graph->GetYaxis()->SetLabelSize(14); // pixels
    R_graph->GetYaxis()->SetLabelOffset(0.01);
    R_graph->GetYaxis()->SetTitleFont(63);
    R_graph->GetYaxis()->SetTitleSize(14); // pixels
    R_graph->GetYaxis()->SetTitleOffset(2.8);
    R_graph->GetYaxis()->SetRangeUser(-24.8, -15.2);

    R_graph->GetXaxis()->SetLabelFont(63);
    R_graph->GetXaxis()->SetLabelSize(14); // pixels
    R_graph->GetXaxis()->SetTitleFont(63);
    R_graph->GetXaxis()->SetTitleSize(14); // pixels
    R_graph->GetXaxis()->SetTitleOffset(2);

    gPad->SetTickx(1);
    gPad->SetTicky(1);

    R_graph->Draw("APZ");
    K1_R->Draw("SAME");
    K2_R->Draw("SAME");

/////////////////////////////////////////////////////////////////////////////////////

    inputFile->Close();

  }

multiplePlots_canv->cd();

  for (uint fileNum = 0; fileNum < dataset_files.size(); ++fileNum)
  {
    double change = 0.24 * fileNum;

    TPaveText* datasetTitle = new TPaveText(0.1+change,0.91,0.2+change,1, "NDC");
    datasetTitle->AddText(dataset_names.at(fileNum).c_str());
    datasetTitle->SetBorderSize(0);
    datasetTitle->SetTextColor(dataset_colors.at(fileNum));
    datasetTitle->SetFillStyle(0);
    datasetTitle->SetTextFont(63);
    datasetTitle->SetTextSize(30);
    datasetTitle->Draw("SAME");

  }

multiplePlots_canv->SaveAs("fitStartScan_R_fourDatasets.png");






	return 1;

}
