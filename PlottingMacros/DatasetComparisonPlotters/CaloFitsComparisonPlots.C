// 3-31-20: Macro to take per calorimeter plots from different datasets and plot them against each other.

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
#include <TPaveStats.h>

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

bool saveImages = true;

int CaloFitsComparisonPlots()
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


  vector<string> dataset_files = {"/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/60h/SingleIteration/CaloFits/perCaloPlots.root",
                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/HighKick/SingleIteration/CaloFits/perCaloPlots.root",
                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/9d/SingleIteration/CaloFits/perCaloPlots.root",
                                  "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/Endgame/SingleIteration/CaloFits/perCaloPlots.root"};

  if(dataset_files.size() > 4) return -1;

  Float_t small = 1e-5;

/////////////////////////////////////////////////////////////////////////////////////
/*
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

      TGraphErrors* R_graph = (TGraphErrors*) inputFile->Get("topDir/FitPass0/FullRatioFit/Individual_Graphs_and_Hists/FullRatioFit_R_Graph");

      R_graph->Fit("pol0", "Q0");
      double fitVal = R_graph->GetFunction("pol0")->GetParameter(0);

      for (int i = 0; i < R_graph->GetN(); ++i)
      {
        double xPoint, yPoint;
        R_graph->GetPoint(i, xPoint, yPoint);
        R_graph->SetPoint(i, xPoint, yPoint - fitVal);
      }

    combined_canv->cd();

    R_graph->SetLineColor(dataset_colors.at(fileNum));
    R_graph->SetMarkerColor(dataset_colors.at(fileNum));

    R_graph->SetMarkerSize(1);
    R_graph->SetMarkerStyle(dataset_markers.at(fileNum));


    if(fileNum == 0){
      R_graph->GetXaxis()->SetLimits(0, 25);
      R_graph->GetYaxis()->SetTitle("#DeltaR (ppm)");
      R_graph->Draw("APZ");
    } 
    else R_graph->Draw("PZSAME");


    legend->AddEntry(R_graph, dataset_names.at(fileNum).c_str(), "pl");

    inputFile->Close();
  }

  legend->Draw("SAME");

  combined_canv->Update();
  if(saveImages) combined_canv->SaveAs("caloFits_R_dataset_comparison.png");

/////////////////////////////////////////////////////////////////////////////////////
*/
  // four plots on one canvas

  auto multiplePlots_chi2_canv = new TCanvas("multiplePlots_chi2_canv", "multiplePlots_chi2_canv", 10, 10, 1000, 650);
  multiplePlots_chi2_canv->Divide(2,2,small,small);


  for (uint fileNum = 0; fileNum < dataset_files.size(); ++fileNum)
  {
    TFile* inputFile = TFile::Open(dataset_files.at(fileNum).c_str());
    if(inputFile == 0){
      printf("Error: cannot open file\n");
      return 0;
    }

    multiplePlots_chi2_canv->cd(fileNum+1);

    if(fileNum == 0){
      gPad->SetTopMargin(0.1);
      gPad->SetBottomMargin(0.005);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.02);
    }
    else if(fileNum == 1){
      gPad->SetTopMargin(0.1);
      gPad->SetBottomMargin(0.005);
      gPad->SetLeftMargin(0.02);
      gPad->SetRightMargin(0.15);
    }    
    else if(fileNum == 2){
      gPad->SetTopMargin(0.005);
      gPad->SetBottomMargin(0.1);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.02);
    }    
    else if(fileNum == 3){
      gPad->SetTopMargin(0.005);
      gPad->SetBottomMargin(0.1);
      gPad->SetLeftMargin(0.02);
      gPad->SetRightMargin(0.15);
    }


    TGraphErrors* chi2_graph = (TGraphErrors*) inputFile->Get("topDir/FitPass0/FullRatioFit/Individual_Graphs_and_Hists/FullRatioFit_chi2ndf_graph");
    chi2_graph->GetYaxis()->SetRangeUser(0.89, 1.11);
    chi2_graph->GetXaxis()->SetLimits(0, 25);

    chi2_graph->GetYaxis()->CenterTitle();
    chi2_graph->GetYaxis()->SetLabelFont(63);
    chi2_graph->GetYaxis()->SetLabelSize(15); // pixels
    chi2_graph->GetYaxis()->SetLabelOffset(0.01);
    chi2_graph->GetYaxis()->SetTitleFont(63);
    chi2_graph->GetYaxis()->SetTitleSize(15); // pixels
    chi2_graph->GetYaxis()->SetTitleOffset(2.8);

    chi2_graph->GetXaxis()->SetLabelFont(63);
    chi2_graph->GetXaxis()->SetLabelSize(15); // pixels
    chi2_graph->GetXaxis()->SetTitleFont(63);
    chi2_graph->GetXaxis()->SetTitleSize(15); // pixels
    chi2_graph->GetXaxis()->SetTitleOffset(2);

    if(fileNum == 1 || fileNum == 3) chi2_graph->GetYaxis()->SetLabelOffset(1);

    chi2_graph->Draw("AP");

    TLine *line = new TLine(0, 1, 25, 1);
    line->SetLineColor(2);
    line->SetLineStyle(2);
    line->Draw("SAME");

    TPaveText* textBox = new TPaveText(GetNDCX(1),GetNDCY(1.07),GetNDCX(10),GetNDCY(1.1), "NDC");
    textBox->AddText(dataset_names.at(fileNum).c_str());
    textBox->SetBorderSize(0);
    textBox->SetTextColor(dataset_colors.at(fileNum));
    textBox->SetFillStyle(0);
    textBox->Draw("SAME");

    multiplePlots_chi2_canv->Update();

    inputFile->Close();
  }
    
  // multiplePlots_chi2_canv->cd();
      // TBox* whiteBox = new TBox(1, 1, 5, 5);
      // whiteBox->SetFillStyle(1001);
      // whiteBox->SetFillColor(0);
      // whiteBox->Draw("SAME");

  if(saveImages) multiplePlots_chi2_canv->SaveAs("caloFits_chi2_fourDatsets.png");

  /////////////////////////////////////////////////////////////////////////////////////

  gStyle->SetOptFit(111);
              
  auto multiplePlots_R_canv = new TCanvas("multiplePlots_R_canv", "multiplePlots_R_canv", 710, 10, 1000, 650);
  // Float_t small = 1e-5;
  multiplePlots_R_canv->Divide(2,2,small,small);


  for (uint fileNum = 0; fileNum < dataset_files.size(); ++fileNum)
  {
    TFile* inputFile = TFile::Open(dataset_files.at(fileNum).c_str());
    if(inputFile == 0){
      printf("Error: cannot open file\n");
      return 0;
    }

    multiplePlots_R_canv->cd(fileNum+1);

    if(fileNum == 0){
      gPad->SetTopMargin(0.1);
      gPad->SetBottomMargin(0.005);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.02);
    }
    else if(fileNum == 1){
      gPad->SetTopMargin(0.1);
      gPad->SetBottomMargin(0.005);
      gPad->SetLeftMargin(0.02);
      gPad->SetRightMargin(0.15);
    }    
    else if(fileNum == 2){
      gPad->SetTopMargin(0.005);
      gPad->SetBottomMargin(0.1);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.02);
    }    
    else if(fileNum == 3){
      gPad->SetTopMargin(0.005);
      gPad->SetBottomMargin(0.1);
      gPad->SetLeftMargin(0.02);
      gPad->SetRightMargin(0.15);
    }


    TGraphErrors* R_graph = (TGraphErrors*) inputFile->Get("topDir/FitPass0/FullRatioFit/Individual_Graphs_and_Hists/FullRatioFit_R_Graph");
    R_graph->GetYaxis()->SetRangeUser(-44, 4);
    R_graph->GetXaxis()->SetLimits(0, 25);

    R_graph->GetYaxis()->CenterTitle();
    R_graph->GetYaxis()->SetLabelFont(63);
    R_graph->GetYaxis()->SetLabelSize(15); // pixels
    R_graph->GetYaxis()->SetLabelOffset(0.01);
    R_graph->GetYaxis()->SetTitleFont(63);
    R_graph->GetYaxis()->SetTitleSize(15); // pixels
    R_graph->GetYaxis()->SetTitleOffset(2.8);

    R_graph->GetXaxis()->SetLabelFont(63);
    R_graph->GetXaxis()->SetLabelSize(15); // pixels
    R_graph->GetXaxis()->SetTitleFont(63);
    R_graph->GetXaxis()->SetTitleSize(15); // pixels
    R_graph->GetXaxis()->SetTitleOffset(2);


    R_graph->Fit("pol0", "Q");
    R_graph->GetFunction("pol0")->SetLineColor(2);

    if(fileNum == 1 || fileNum == 3) R_graph->GetYaxis()->SetLabelOffset(1);

    R_graph->Draw("AP");

    multiplePlots_R_canv->Update();

    // TPaveText* textBox = new TPaveText(GetNDCX(2),GetNDCY(-2),GetNDCX(9),GetNDCY(4), "NDC");
    TPaveText* textBox = new TPaveText(GetNDCX(1),GetNDCY(-42),GetNDCX(9),GetNDCY(-36), "NDC");
    textBox->AddText(dataset_names.at(fileNum).c_str());
    textBox->SetBorderSize(0);
    textBox->SetTextColor(dataset_colors.at(fileNum));
    textBox->SetFillStyle(0);
    textBox->Draw("SAME");


    TPaveStats *ps = (TPaveStats*)gPad->GetPrimitive("stats");
    ps->SetBorderSize(1);

    ps->SetX1NDC(0.5);
    ps->SetX2NDC(0.8);
    ps->SetY1NDC(0.1);
    ps->SetY2NDC(0.3);

    gPad->Modified();

    // ps->SetX1NDC(GetNDCX(1));
    // ps->SetY1NDC(GetNDCY(-43));
    // ps->SetX2NDC(GetNDCX(10));
    // ps->SetY2NDC(GetNDCY(-35));

    ps->SetX1NDC(GetNDCX(11));
    ps->SetY1NDC(GetNDCY(-5));
    ps->SetX2NDC(GetNDCX(25));
    ps->SetY2NDC(GetNDCY(4));

    ps->Draw("SAME");

    inputFile->Close();
  }
    
  if(saveImages) multiplePlots_R_canv->SaveAs("caloFits_R_fourDatsets.png");



	return 1;

}
