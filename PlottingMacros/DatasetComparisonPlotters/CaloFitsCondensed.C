// 3-31-20: I believe this was a macro to take per calo plots and put them on the same canvas in a compact format for my thesis.

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
#include <TPaveStats.h>
#include <TGaxis.h>

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

bool saveImages = true;

int CaloFitsCondensed()
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

  // TGaxis::SetMaxDigits(2); // doesn't want to work

  Float_t small = 1e-5;

/*
  string fileString = "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/Endgame/SingleIteration/CaloFits/perCaloPlots.root";

  TFile* inputFile = TFile::Open(fileString.c_str());
  if(inputFile == 0){
    printf("Error: cannot open file\n");
    return 0;
  }


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  auto parametersCanv_first = new TCanvas("parametersCanv_first", "parametersCanv_first", 310, 10, 800, 600);
  parametersCanv_first->Divide(2,3,0.01,small);
  // parametersCanv_first->Divide(2,3,small,small);
  // parametersCanv_first->Divide(2,3);


  TGraphErrors* phi_graph = (TGraphErrors*) inputFile->Get("topDir/FitPass0/FullRatioFit/Individual_Graphs_and_Hists/FullRatioFit_phi_Graph");
  TGraphErrors* A_graph = (TGraphErrors*) inputFile->Get("topDir/FitPass0/FullRatioFit/Individual_Graphs_and_Hists/FullRatioFit_A_Graph");
  TGraphErrors* w_cbo_graph = (TGraphErrors*) inputFile->Get("topDir/FitPass0/FullRatioFit/Individual_Graphs_and_Hists/FullRatioFit_omega_cbo_Graph");
  TGraphErrors* tau_cbo_graph = (TGraphErrors*) inputFile->Get("topDir/FitPass0/FullRatioFit/Individual_Graphs_and_Hists/FullRatioFit_tau_cbo_Graph");
  TGraphErrors* A_cbo_N_graph = (TGraphErrors*) inputFile->Get("topDir/FitPass0/FullRatioFit/Individual_Graphs_and_Hists/FullRatioFit_A_cbo-N_Graph");
  TGraphErrors* phi_cbo_N_graph = (TGraphErrors*) inputFile->Get("topDir/FitPass0/FullRatioFit/Individual_Graphs_and_Hists/FullRatioFit_phi_cbo-N_Graph");

/////////////////////////////////////////////////////////////////////////////////////

  for (int pointNo = 0; pointNo < A_cbo_N_graph->GetN(); ++pointNo)
  {
  	double x,y,xErr,yErr;
  	A_cbo_N_graph->GetPoint(pointNo, x, y);
  	xErr = A_cbo_N_graph->GetErrorX(pointNo);
  	yErr = A_cbo_N_graph->GetErrorY(pointNo);
  	A_cbo_N_graph->SetPoint(pointNo, x, y*1000);
  	A_cbo_N_graph->SetPointError(pointNo, xErr, yErr*1000);
  }

/////////////////////////////////////////////////////////////////////////////////////

  parametersCanv_first->cd(1);
  	
  	  gPad->SetTopMargin(0.01);
      gPad->SetBottomMargin(0.005);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.02);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    phi_graph->GetYaxis()->CenterTitle();
    phi_graph->GetYaxis()->SetLabelFont(63);
    phi_graph->GetYaxis()->SetLabelSize(13); // pixels
    phi_graph->GetYaxis()->SetLabelOffset(0.01);
    phi_graph->GetYaxis()->SetTitleFont(63);
    phi_graph->GetYaxis()->SetTitleSize(14); // pixels
    phi_graph->GetYaxis()->SetTitleOffset(3.0);
   	phi_graph->Draw("AP");
  

  parametersCanv_first->cd(2);

  	  gPad->SetTopMargin(0.01);
      gPad->SetBottomMargin(0.005);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.02);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    A_graph->GetYaxis()->CenterTitle();
    A_graph->GetYaxis()->SetLabelFont(63);
    A_graph->GetYaxis()->SetLabelSize(13); // pixels
    A_graph->GetYaxis()->SetLabelOffset(0.01);
    A_graph->GetYaxis()->SetTitleFont(63);
    A_graph->GetYaxis()->SetTitleSize(14); // pixels
    A_graph->GetYaxis()->SetTitleOffset(3.0);
  	A_graph->Draw("AP");
  

  parametersCanv_first->cd(3);

  	  gPad->SetTopMargin(0.005);
      gPad->SetBottomMargin(0.005);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.02);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    w_cbo_graph->GetYaxis()->CenterTitle();
    w_cbo_graph->GetYaxis()->SetLabelFont(63);
    w_cbo_graph->GetYaxis()->SetLabelSize(13); // pixels
    w_cbo_graph->GetYaxis()->SetLabelOffset(0.01);
    w_cbo_graph->GetYaxis()->SetTitleFont(63);
    w_cbo_graph->GetYaxis()->SetTitleSize(14); // pixels
    w_cbo_graph->GetYaxis()->SetTitleOffset(3.0);
  	w_cbo_graph->Draw("AP");
  

  parametersCanv_first->cd(4);

  	  gPad->SetTopMargin(0.005);
      gPad->SetBottomMargin(0.005);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.02);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    tau_cbo_graph->GetYaxis()->CenterTitle();
    tau_cbo_graph->GetYaxis()->SetLabelFont(63);
    tau_cbo_graph->GetYaxis()->SetLabelSize(13); // pixels
    tau_cbo_graph->GetYaxis()->SetLabelOffset(0.01);
    tau_cbo_graph->GetYaxis()->SetTitleFont(63);
    tau_cbo_graph->GetYaxis()->SetTitleSize(14); // pixels
    tau_cbo_graph->GetYaxis()->SetTitleOffset(3.0);
  	tau_cbo_graph->Draw("AP");
  

  parametersCanv_first->cd(5);

  	  gPad->SetTopMargin(0.005);
      gPad->SetBottomMargin(0.15);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.02);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    A_cbo_N_graph->GetYaxis()->CenterTitle();
    A_cbo_N_graph->GetYaxis()->SetLabelFont(63);
    A_cbo_N_graph->GetYaxis()->SetLabelSize(13); // pixels
    A_cbo_N_graph->GetYaxis()->SetLabelOffset(0.01);
    A_cbo_N_graph->GetYaxis()->SetTitleFont(63);
    A_cbo_N_graph->GetYaxis()->SetTitleSize(14); // pixels
    A_cbo_N_graph->GetYaxis()->SetTitleOffset(3.0);
    A_cbo_N_graph->GetYaxis()->SetTitle((string(A_cbo_N_graph->GetYaxis()->GetTitle()) + " (#times 10^{-3})").c_str());
  	
    A_cbo_N_graph->GetXaxis()->SetLabelFont(63);
    A_cbo_N_graph->GetXaxis()->SetLabelSize(13); // pixels
    A_cbo_N_graph->GetXaxis()->SetTitleFont(63);
    A_cbo_N_graph->GetXaxis()->SetTitleSize(14); // pixels
    A_cbo_N_graph->GetXaxis()->SetTitleOffset(3);

  	A_cbo_N_graph->Draw("AP");
  

  parametersCanv_first->cd(6);

  	  gPad->SetTopMargin(0.005);
      gPad->SetBottomMargin(0.15);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.02);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    phi_cbo_N_graph->GetYaxis()->CenterTitle();
    phi_cbo_N_graph->GetYaxis()->SetLabelFont(63);
    phi_cbo_N_graph->GetYaxis()->SetLabelSize(13); // pixels
    phi_cbo_N_graph->GetYaxis()->SetLabelOffset(0.01);
    phi_cbo_N_graph->GetYaxis()->SetTitleFont(63);
    phi_cbo_N_graph->GetYaxis()->SetTitleSize(14); // pixels
    phi_cbo_N_graph->GetYaxis()->SetTitleOffset(3.0);

    phi_cbo_N_graph->GetXaxis()->SetLabelFont(63);
    phi_cbo_N_graph->GetXaxis()->SetLabelSize(13); // pixels
    phi_cbo_N_graph->GetXaxis()->SetTitleFont(63);
    phi_cbo_N_graph->GetXaxis()->SetTitleSize(14); // pixels
    phi_cbo_N_graph->GetXaxis()->SetTitleOffset(3);

  	phi_cbo_N_graph->Draw("AP");


/////////////////////////////////////////////////////////////////////////////////////

  	parametersCanv_first->SaveAs("Endgame_caloFitPars_first.png");

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  auto parametersCanv_second = new TCanvas("parametersCanv_second", "parametersCanv_second", 310, 10, 800, 600);
  parametersCanv_second->Divide(2,3,0.01,small);


  TGraphErrors* A_cbo_A_graph = (TGraphErrors*) inputFile->Get("topDir/FitPass0/FullRatioFit/Individual_Graphs_and_Hists/FullRatioFit_A_cbo-A_Graph");
  TGraphErrors* phi_cbo_A_graph = (TGraphErrors*) inputFile->Get("topDir/FitPass0/FullRatioFit/Individual_Graphs_and_Hists/FullRatioFit_phi_cbo-A_Graph");
  TGraphErrors* A_cbo_phi_graph = (TGraphErrors*) inputFile->Get("topDir/FitPass0/FullRatioFit/Individual_Graphs_and_Hists/FullRatioFit_A_cbo-phi_Graph");
  TGraphErrors* phi_cbo_phi_graph = (TGraphErrors*) inputFile->Get("topDir/FitPass0/FullRatioFit/Individual_Graphs_and_Hists/FullRatioFit_phi_cbo-phi_Graph");
  TGraphErrors* A_2cbo_N_graph = (TGraphErrors*) inputFile->Get("topDir/FitPass0/FullRatioFit/Individual_Graphs_and_Hists/FullRatioFit_A_2cbo-N_Graph");
  TGraphErrors* phi_2cbo_N_graph = (TGraphErrors*) inputFile->Get("topDir/FitPass0/FullRatioFit/Individual_Graphs_and_Hists/FullRatioFit_phi_2cbo-N_Graph");

/////////////////////////////////////////////////////////////////////////////////////

  for (int pointNo = 0; pointNo < A_cbo_A_graph->GetN(); ++pointNo)
  {
  	double x,y,xErr,yErr;
  	A_cbo_A_graph->GetPoint(pointNo, x, y);
  	xErr = A_cbo_A_graph->GetErrorX(pointNo);
  	yErr = A_cbo_A_graph->GetErrorY(pointNo);
  	A_cbo_A_graph->SetPoint(pointNo, x, y*1000);
  	A_cbo_A_graph->SetPointError(pointNo, xErr, yErr*1000);

    A_cbo_phi_graph->GetPoint(pointNo, x, y);
    xErr = A_cbo_phi_graph->GetErrorX(pointNo);
    yErr = A_cbo_phi_graph->GetErrorY(pointNo);
    A_cbo_phi_graph->SetPoint(pointNo, x, y*1000);
    A_cbo_phi_graph->SetPointError(pointNo, xErr, yErr*1000);

    A_2cbo_N_graph->GetPoint(pointNo, x, y);
    xErr = A_2cbo_N_graph->GetErrorX(pointNo);
    yErr = A_2cbo_N_graph->GetErrorY(pointNo);
    A_2cbo_N_graph->SetPoint(pointNo, x, y*1000);
    A_2cbo_N_graph->SetPointError(pointNo, xErr, yErr*1000);
  }

/////////////////////////////////////////////////////////////////////////////////////

  parametersCanv_second->cd(1);
  	
  	  gPad->SetTopMargin(0.01);
      gPad->SetBottomMargin(0.005);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.02);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    A_cbo_A_graph->GetYaxis()->CenterTitle();
    A_cbo_A_graph->GetYaxis()->SetLabelFont(63);
    A_cbo_A_graph->GetYaxis()->SetLabelSize(13); // pixels
    A_cbo_A_graph->GetYaxis()->SetLabelOffset(0.01);
    A_cbo_A_graph->GetYaxis()->SetTitleFont(63);
    A_cbo_A_graph->GetYaxis()->SetTitleSize(14); // pixels
    A_cbo_A_graph->GetYaxis()->SetTitleOffset(3.0);
    A_cbo_A_graph->GetYaxis()->SetTitle((string(A_cbo_A_graph->GetYaxis()->GetTitle()) + " (#times 10^{-3})").c_str());
   	A_cbo_A_graph->Draw("AP");
  

  parametersCanv_second->cd(2);

  	  gPad->SetTopMargin(0.01);
      gPad->SetBottomMargin(0.005);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.02);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    phi_cbo_A_graph->GetYaxis()->CenterTitle();
    phi_cbo_A_graph->GetYaxis()->SetLabelFont(63);
    phi_cbo_A_graph->GetYaxis()->SetLabelSize(13); // pixels
    phi_cbo_A_graph->GetYaxis()->SetLabelOffset(0.01);
    phi_cbo_A_graph->GetYaxis()->SetTitleFont(63);
    phi_cbo_A_graph->GetYaxis()->SetTitleSize(14); // pixels
    phi_cbo_A_graph->GetYaxis()->SetTitleOffset(3.0);
  	phi_cbo_A_graph->Draw("AP");
  

  parametersCanv_second->cd(3);

  	  gPad->SetTopMargin(0.005);
      gPad->SetBottomMargin(0.005);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.02);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    A_cbo_phi_graph->GetYaxis()->CenterTitle();
    A_cbo_phi_graph->GetYaxis()->SetLabelFont(63);
    A_cbo_phi_graph->GetYaxis()->SetLabelSize(13); // pixels
    A_cbo_phi_graph->GetYaxis()->SetLabelOffset(0.01);
    A_cbo_phi_graph->GetYaxis()->SetTitleFont(63);
    A_cbo_phi_graph->GetYaxis()->SetTitleSize(14); // pixels
    A_cbo_phi_graph->GetYaxis()->SetTitleOffset(3.0);
    A_cbo_phi_graph->GetYaxis()->SetTitle((string(A_cbo_phi_graph->GetYaxis()->GetTitle()) + " (#times 10^{-3})").c_str());
  	A_cbo_phi_graph->Draw("AP");
  

  parametersCanv_second->cd(4);

  	  gPad->SetTopMargin(0.005);
      gPad->SetBottomMargin(0.005);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.02);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    phi_cbo_phi_graph->GetYaxis()->CenterTitle();
    phi_cbo_phi_graph->GetYaxis()->SetLabelFont(63);
    phi_cbo_phi_graph->GetYaxis()->SetLabelSize(13); // pixels
    phi_cbo_phi_graph->GetYaxis()->SetLabelOffset(0.01);
    phi_cbo_phi_graph->GetYaxis()->SetTitleFont(63);
    phi_cbo_phi_graph->GetYaxis()->SetTitleSize(14); // pixels
    phi_cbo_phi_graph->GetYaxis()->SetTitleOffset(3.0);
  	phi_cbo_phi_graph->Draw("AP");
  

  parametersCanv_second->cd(5);

  	  gPad->SetTopMargin(0.005);
      gPad->SetBottomMargin(0.15);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.02);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    A_2cbo_N_graph->GetYaxis()->CenterTitle();
    A_2cbo_N_graph->GetYaxis()->SetLabelFont(63);
    A_2cbo_N_graph->GetYaxis()->SetLabelSize(13); // pixels
    A_2cbo_N_graph->GetYaxis()->SetLabelOffset(0.01);
    A_2cbo_N_graph->GetYaxis()->SetTitleFont(63);
    A_2cbo_N_graph->GetYaxis()->SetTitleSize(14); // pixels
    A_2cbo_N_graph->GetYaxis()->SetTitleOffset(3.0);
    A_2cbo_N_graph->GetYaxis()->SetTitle((string(A_2cbo_N_graph->GetYaxis()->GetTitle()) + " (#times 10^{-3})").c_str());
  	
    A_2cbo_N_graph->GetXaxis()->SetLabelFont(63);
    A_2cbo_N_graph->GetXaxis()->SetLabelSize(13); // pixels
    A_2cbo_N_graph->GetXaxis()->SetTitleFont(63);
    A_2cbo_N_graph->GetXaxis()->SetTitleSize(14); // pixels
    A_2cbo_N_graph->GetXaxis()->SetTitleOffset(3);

  	A_2cbo_N_graph->Draw("AP");
  

  parametersCanv_second->cd(6);

  	  gPad->SetTopMargin(0.005);
      gPad->SetBottomMargin(0.15);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.02);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    phi_2cbo_N_graph->GetYaxis()->CenterTitle();
    phi_2cbo_N_graph->GetYaxis()->SetLabelFont(63);
    phi_2cbo_N_graph->GetYaxis()->SetLabelSize(13); // pixels
    phi_2cbo_N_graph->GetYaxis()->SetLabelOffset(0.01);
    phi_2cbo_N_graph->GetYaxis()->SetTitleFont(63);
    phi_2cbo_N_graph->GetYaxis()->SetTitleSize(14); // pixels
    phi_2cbo_N_graph->GetYaxis()->SetTitleOffset(3.0);

    phi_2cbo_N_graph->GetXaxis()->SetLabelFont(63);
    phi_2cbo_N_graph->GetXaxis()->SetLabelSize(13); // pixels
    phi_2cbo_N_graph->GetXaxis()->SetTitleFont(63);
    phi_2cbo_N_graph->GetXaxis()->SetTitleSize(14); // pixels
    phi_2cbo_N_graph->GetXaxis()->SetTitleOffset(3);

  	phi_2cbo_N_graph->Draw("AP");


/////////////////////////////////////////////////////////////////////////////////////

  	parametersCanv_second->SaveAs("Endgame_caloFitPars_second.png");

*/
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


  string fileString_fitStartScan = "/gm2/data/users/nkinnaird/Ratio/ThesisAnalysis/Endgame/SingleIteration/FitStartScans/fitStartScan-fixed-cbo-A-A.root";

  TFile* inputFile_fitStartScan = TFile::Open(fileString_fitStartScan.c_str());
  if(inputFile_fitStartScan == 0){
    printf("Error: cannot open file\n");
    return 0;
  }


  auto fitStartCanv = new TCanvas("fitStartCanv", "fitStartCanv", 310, 10, 800, 600);
  fitStartCanv->Divide(2,3,0.01,small);

    // phi_canv->GetListOfPrimitives()->Print();

/////////////////////////////////////////////////////////////////////////////////////

    TCanvas* phi_canv = (TCanvas*) inputFile_fitStartScan->Get("topDir/Added/FullRatio/FullRatio_phi_FS_Canv");

    TGraphErrors* phi_graph_FS = (TGraphErrors*) ((TGraphErrors*) phi_canv->GetListOfPrimitives()->FindObject("FullRatio_phi_Vs_FS")->Clone());
    TGraph* K1_phi = (TGraph*) ((TGraph*) phi_canv->GetListOfPrimitives()->FindObject("KawallBandOne")->Clone());
    TGraph* K2_phi = (TGraph*) ((TGraph*) phi_canv->GetListOfPrimitives()->FindObject("KawallBandTwo")->Clone());

  fitStartCanv->cd(1);
    
      gPad->SetTopMargin(0.1);
      gPad->SetBottomMargin(0.005);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.02);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    phi_graph_FS->GetYaxis()->CenterTitle();
    phi_graph_FS->GetYaxis()->SetLabelFont(63);
    phi_graph_FS->GetYaxis()->SetLabelSize(13); // pixels
    phi_graph_FS->GetYaxis()->SetLabelOffset(0.01);
    phi_graph_FS->GetYaxis()->SetTitleFont(63);
    phi_graph_FS->GetYaxis()->SetTitleSize(14); // pixels
    phi_graph_FS->GetYaxis()->SetTitleOffset(3.0);
    phi_graph_FS->Draw("APZ");
    K1_phi->Draw("SAME");
    K2_phi->Draw("SAME");


    TCanvas* A_canv = (TCanvas*) inputFile_fitStartScan->Get("topDir/Added/FullRatio/FullRatio_A_FS_Canv");

    TGraphErrors* A_graph_FS = (TGraphErrors*) ((TGraphErrors*) A_canv->GetListOfPrimitives()->FindObject("FullRatio_A_Vs_FS")->Clone());
    TGraph* K1_A = (TGraph*) ((TGraph*) A_canv->GetListOfPrimitives()->FindObject("KawallBandOne")->Clone());
    TGraph* K2_A = (TGraph*) ((TGraph*) A_canv->GetListOfPrimitives()->FindObject("KawallBandTwo")->Clone());

  fitStartCanv->cd(2);

      gPad->SetTopMargin(0.1);
      gPad->SetBottomMargin(0.005);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.02);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    // A_graph_FS->GetYaxis()->SetNdivisions(710, false);
  // TGaxis::SetMaxDigits(2);

    A_graph_FS->GetYaxis()->CenterTitle();
    A_graph_FS->GetYaxis()->SetLabelFont(63);
    A_graph_FS->GetYaxis()->SetLabelSize(13); // pixels
    A_graph_FS->GetYaxis()->SetLabelOffset(0.01);
    A_graph_FS->GetYaxis()->SetTitleFont(63);
    A_graph_FS->GetYaxis()->SetTitleSize(14); // pixels
    A_graph_FS->GetYaxis()->SetTitleOffset(3.5);
    A_graph_FS->Draw("APZ");
    K1_A->Draw("SAME");
    K2_A->Draw("SAME");  

    // gPad->Modified();


    TCanvas* w_cbo_canv = (TCanvas*) inputFile_fitStartScan->Get("topDir/Added/FullRatio/FullRatio_omega_cbo_FS_Canv");

    TGraphErrors* w_cbo_graph_FS = (TGraphErrors*) ((TGraphErrors*) w_cbo_canv->GetListOfPrimitives()->FindObject("FullRatio_omega_cbo_Vs_FS")->Clone());
    TGraph* K1_w_cbo = (TGraph*) ((TGraph*) w_cbo_canv->GetListOfPrimitives()->FindObject("KawallBandOne")->Clone());
    TGraph* K2_w_cbo = (TGraph*) ((TGraph*) w_cbo_canv->GetListOfPrimitives()->FindObject("KawallBandTwo")->Clone());

  fitStartCanv->cd(3);

      gPad->SetTopMargin(0.005);
      gPad->SetBottomMargin(0.005);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.02);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    w_cbo_graph_FS->GetYaxis()->CenterTitle();
    w_cbo_graph_FS->GetYaxis()->SetLabelFont(63);
    w_cbo_graph_FS->GetYaxis()->SetLabelSize(13); // pixels
    w_cbo_graph_FS->GetYaxis()->SetLabelOffset(0.01);
    w_cbo_graph_FS->GetYaxis()->SetTitleFont(63);
    w_cbo_graph_FS->GetYaxis()->SetTitleSize(14); // pixels
    w_cbo_graph_FS->GetYaxis()->SetTitleOffset(3.0);
    w_cbo_graph_FS->Draw("APZ");
    K1_w_cbo->Draw("SAME");
    K2_w_cbo->Draw("SAME");  



    TCanvas* phi_cbo_N_canv = (TCanvas*) inputFile_fitStartScan->Get("topDir/Added/FullRatio/FullRatio_phi_cbo-N_FS_Canv");

    TGraphErrors* phi_cbo_N_graph_FS = (TGraphErrors*) ((TGraphErrors*) phi_cbo_N_canv->GetListOfPrimitives()->FindObject("FullRatio_phi_cbo-N_Vs_FS")->Clone());
    TGraph* K1_phi_cbo_N = (TGraph*) ((TGraph*) phi_cbo_N_canv->GetListOfPrimitives()->FindObject("KawallBandOne")->Clone());
    TGraph* K2_phi_cbo_N = (TGraph*) ((TGraph*) phi_cbo_N_canv->GetListOfPrimitives()->FindObject("KawallBandTwo")->Clone());

  fitStartCanv->cd(4);

      gPad->SetTopMargin(0.005);
      gPad->SetBottomMargin(0.005);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.02);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    phi_cbo_N_graph_FS->GetYaxis()->CenterTitle();
    phi_cbo_N_graph_FS->GetYaxis()->SetLabelFont(63);
    phi_cbo_N_graph_FS->GetYaxis()->SetLabelSize(13); // pixels
    phi_cbo_N_graph_FS->GetYaxis()->SetLabelOffset(0.01);
    phi_cbo_N_graph_FS->GetYaxis()->SetTitleFont(63);
    phi_cbo_N_graph_FS->GetYaxis()->SetTitleSize(14); // pixels
    phi_cbo_N_graph_FS->GetYaxis()->SetTitleOffset(3.0);
    phi_cbo_N_graph_FS->Draw("APZ");
    K1_phi_cbo_N->Draw("SAME");
    K2_phi_cbo_N->Draw("SAME");  



    TCanvas* A_cbo_N_canv = (TCanvas*) inputFile_fitStartScan->Get("topDir/Added/FullRatio/FullRatio_A_cbo-N_FS_Canv");

    TGraphErrors* A_cbo_N_graph_FS = (TGraphErrors*) ((TGraphErrors*) A_cbo_N_canv->GetListOfPrimitives()->FindObject("FullRatio_A_cbo-N_Vs_FS")->Clone());
    TGraph* K1_A_cbo_N = (TGraph*) ((TGraph*) A_cbo_N_canv->GetListOfPrimitives()->FindObject("KawallBandOne")->Clone());
    TGraph* K2_A_cbo_N = (TGraph*) ((TGraph*) A_cbo_N_canv->GetListOfPrimitives()->FindObject("KawallBandTwo")->Clone());

  fitStartCanv->cd(5);

      gPad->SetTopMargin(0.005);
      gPad->SetBottomMargin(0.15);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.02);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    A_cbo_N_graph_FS->GetYaxis()->CenterTitle();
    A_cbo_N_graph_FS->GetYaxis()->SetLabelFont(63);
    A_cbo_N_graph_FS->GetYaxis()->SetLabelSize(13); // pixels
    A_cbo_N_graph_FS->GetYaxis()->SetLabelOffset(0.01);
    A_cbo_N_graph_FS->GetYaxis()->SetTitleFont(63);
    A_cbo_N_graph_FS->GetYaxis()->SetTitleSize(14); // pixels
    A_cbo_N_graph_FS->GetYaxis()->SetTitleOffset(3.0);
    A_cbo_N_graph_FS->GetYaxis()->SetTitle((string(A_cbo_N_graph_FS->GetYaxis()->GetTitle()) + " (#times 10^{-3})").c_str());
    
    A_cbo_N_graph_FS->GetXaxis()->SetLabelFont(63);
    A_cbo_N_graph_FS->GetXaxis()->SetLabelSize(13); // pixels
    A_cbo_N_graph_FS->GetXaxis()->SetTitleFont(63);
    A_cbo_N_graph_FS->GetXaxis()->SetTitleSize(14); // pixels
    A_cbo_N_graph_FS->GetXaxis()->SetTitleOffset(3);

    A_cbo_N_graph_FS->Draw("APZ");
    K1_A_cbo_N->Draw("SAME");
    K2_A_cbo_N->Draw("SAME");  



    TCanvas* A_cbo_phi_canv = (TCanvas*) inputFile_fitStartScan->Get("topDir/Added/FullRatio/FullRatio_A_cbo-phi_FS_Canv");

    TGraphErrors* A_cbo_phi_graph_FS = (TGraphErrors*) ((TGraphErrors*) A_cbo_phi_canv->GetListOfPrimitives()->FindObject("FullRatio_A_cbo-phi_Vs_FS")->Clone());
    TGraph* K1_A_cbo_phi = (TGraph*) ((TGraph*) A_cbo_phi_canv->GetListOfPrimitives()->FindObject("KawallBandOne")->Clone());
    TGraph* K2_A_cbo_phi = (TGraph*) ((TGraph*) A_cbo_phi_canv->GetListOfPrimitives()->FindObject("KawallBandTwo")->Clone());

  fitStartCanv->cd(6);

      gPad->SetTopMargin(0.005);
      gPad->SetBottomMargin(0.15);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.02);
    gPad->SetTickx(1);
    gPad->SetTicky(1);

    A_cbo_phi_graph_FS->GetYaxis()->CenterTitle();
    A_cbo_phi_graph_FS->GetYaxis()->SetLabelFont(63);
    A_cbo_phi_graph_FS->GetYaxis()->SetLabelSize(13); // pixels
    A_cbo_phi_graph_FS->GetYaxis()->SetLabelOffset(0.01);
    A_cbo_phi_graph_FS->GetYaxis()->SetTitleFont(63);
    A_cbo_phi_graph_FS->GetYaxis()->SetTitleSize(14); // pixels
    A_cbo_phi_graph_FS->GetYaxis()->SetTitleOffset(3.0);
    A_cbo_phi_graph_FS->GetYaxis()->SetTitle((string(A_cbo_phi_graph_FS->GetYaxis()->GetTitle()) + " (#times 10^{-3})").c_str());
    
    A_cbo_phi_graph_FS->GetXaxis()->SetLabelFont(63);
    A_cbo_phi_graph_FS->GetXaxis()->SetLabelSize(13); // pixels
    A_cbo_phi_graph_FS->GetXaxis()->SetTitleFont(63);
    A_cbo_phi_graph_FS->GetXaxis()->SetTitleSize(14); // pixels
    A_cbo_phi_graph_FS->GetXaxis()->SetTitleOffset(3);

    A_cbo_phi_graph_FS->Draw("APZ");
    K1_A_cbo_phi->Draw("SAME");
    K2_A_cbo_phi->Draw("SAME");  

/////////////////////////////////////////////////////////////////////////////////////

        fitStartCanv->SaveAs("Endgame_fitStartScan_pars.png");


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

return 1;

}
