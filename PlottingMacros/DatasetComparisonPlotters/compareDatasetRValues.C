// 3-31-20: Macro to plot single point R values from different datasets against each other.

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TLegend.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TImage.h>
#include <TLine.h>

#include <Eigen/Dense>

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;

bool saveImages = true;

int compareDatasetRValues()
{
	double correlationCoefficient = 1;

	// means are in ppm

	double Tmethod_mean_60h = -20.3703; // different software blinding string
	double Tmethod_mean_HighKick = -16.6131;
	double Tmethod_mean_9d = -17.4557;
	double Tmethod_mean_Endgame = -16.8398;

  double Rmethod_mean_60h = -20.5474; // different software blinding string
  double Rmethod_mean_HighKick = -16.8131;
  double Rmethod_mean_9d = -17.5192;
  double Rmethod_mean_Endgame = -16.9319;

	vector<double> Tmethod_datasetMeans = {Tmethod_mean_60h, Tmethod_mean_HighKick, Tmethod_mean_9d, Tmethod_mean_Endgame};
  vector<double> Rmethod_datasetMeans = {Rmethod_mean_60h, Rmethod_mean_HighKick, Rmethod_mean_9d, Rmethod_mean_Endgame};

	// errors are statistical and then systematic
	// units in ppb

	pair<double, double> Tmethod_errors_60h = {1360.0, 0};
	pair<double, double> Tmethod_errors_HighKick = {1150.0, 0};
	pair<double, double> Tmethod_errors_9d = {930.0, 0};
	pair<double, double> Tmethod_errors_Endgame = {650.0, 0};

  pair<double, double> Rmethod_errors_60h = {1360.0, 0};
  pair<double, double> Rmethod_errors_HighKick = {1150.0, 0};
  pair<double, double> Rmethod_errors_9d = {930.0, 0};
  pair<double, double> Rmethod_errors_Endgame = {650.0, 0};

	vector<pair<double, double> > Tmethod_datasetErrors = {Tmethod_errors_60h, Tmethod_errors_HighKick, Tmethod_errors_9d, Tmethod_errors_Endgame};
  vector<pair<double, double> > Rmethod_datasetErrors = {Rmethod_errors_60h, Rmethod_errors_HighKick, Rmethod_errors_9d, Rmethod_errors_Endgame};


/////////////////////////////////////////////////////////////////////////////////////

  	// make plot of dataset means against each other

  gStyle->SetOptStat(000000);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1);
  gStyle->SetLineColor(1);
  gStyle->SetPadRightMargin(.05);
  gStyle->SetPadLeftMargin(.15);


  	TGraphErrors* Tmethod_datasetGraph = new TGraphErrors();
    TGraphErrors* Rmethod_datasetGraph = new TGraphErrors();
      Rmethod_datasetGraph->SetMarkerColor(2);

 	  int pointNo = 0;
  	for (uint i = 0; i < Tmethod_datasetMeans.size(); ++i)
  	{
  		Tmethod_datasetGraph->SetPoint(pointNo, i+1, Tmethod_datasetMeans.at(i));
  		Tmethod_datasetGraph->SetPointError(pointNo, 0, sqrt(pow(Tmethod_datasetErrors.at(i).first,2) + pow(Tmethod_datasetErrors.at(i).second,2))/1000);

      Rmethod_datasetGraph->SetPoint(pointNo, i+1+0.1, Rmethod_datasetMeans.at(i));
      Rmethod_datasetGraph->SetPointError(pointNo, 0, sqrt(pow(Rmethod_datasetErrors.at(i).first,2) + pow(Rmethod_datasetErrors.at(i).second,2))/1000);

      pointNo++;
  	}

   Tmethod_datasetGraph->GetHistogram()->GetXaxis()->Set(Tmethod_datasetMeans.size(),0.5,Tmethod_datasetMeans.size()+0.5);

   for(uint i=0; i < Tmethod_datasetMeans.size(); i++){
      Tmethod_datasetGraph->GetHistogram()->GetXaxis()->SetBinLabel(i+1, dataset_names.at(i).c_str());
   }

   Tmethod_datasetGraph->GetYaxis()->SetTitle("R (ppm, blinded)");
   Tmethod_datasetGraph->GetXaxis()->SetLabelSize(0.08);

    auto canv = new TCanvas("canv","canv",200,10,800,600);
    Tmethod_datasetGraph->Draw("AP");
    Rmethod_datasetGraph->Draw("PSAME");

    canv->Update();

  TLine *line = new TLine(1.5, canv->GetUymin(), 1.5, canv->GetUymax());
  line->SetLineColor(4);
  line->SetLineStyle(2);
  line->SetLineWidth(3);
  line->Draw("SAME");


  auto myLegend = new TLegend(0.6,0.2,.9,0.4);
  myLegend->AddEntry(Tmethod_datasetGraph, "T Method", "p");
  myLegend->AddEntry(Rmethod_datasetGraph, "R Method", "p");
  myLegend->SetBorderSize(0);
  myLegend->Draw("SAME");

  if(saveImages) canv->SaveAs("R_dataset_comparison.png");


	return 1;
}
