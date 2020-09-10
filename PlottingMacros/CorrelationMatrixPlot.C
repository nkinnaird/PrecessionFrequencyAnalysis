// 3-31-20: Macro for correlation matrix plots. Also prints out correlations in a latex table style format.

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2D.h>
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
#include <TMatrixD.h>

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;

int fitPass = 0; // choose which fit pass to make the correlation plot of using this integer


int CorrelationMatrixPlot(std::string filePath)
{

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

    string unneededString = "temp";
	applyDatasetTag(inputFile, unneededString);

	cout.precision(3);
	gStyle->SetMarkerSize(.7);
	gStyle->SetPadLeftMargin(.1);
	gStyle->SetPaintTextFormat("1.3f");
	gStyle->SetGridStyle(0);
	gStyle->SetOptStat(kFALSE);

	// Get correlation matrix
	TMatrixD* corrMatrix_Rmethod = (TMatrixD*) inputFile->Get(Form("topDir/FitPasses/FitPass%i/addedDir/FullRatio/FullRatioFitCorrMatrix", fitPass));

if(corrMatrix_Rmethod != 0)
{
    TF1* fitFunction_Rmethod = (TF1*) ((TGraphErrors*) inputFile->Get(Form("topDir/FitPasses/FitPass%i/addedDir/FullRatio/Added_Times_Full_Ratio_Graph", fitPass)))->GetFunction("fullRatioFitFunc");

	// Labels for covariance matrix                                                                                                                                                                                       
	const int nPars = fitFunction_Rmethod->GetNpar();
	char* labels[nPars];
	for (int parNum = 0; parNum < nPars; ++parNum) labels[parNum] = Form("%s", fitFunction_Rmethod->GetParName(parNum));

	// Draw correlation matrix                                                                                                                                                                                            
	TCanvas* c_corrMat= new TCanvas("c_corrMat","Correlation Matrix",50,10,700,700);
	c_corrMat->SetRightMargin(0.14);

	TH2D* hCorrMatrix_Rmethod = new TH2D("hCorrMatrix_Rmethod","Correlation Matrix", nPars, 0, nPars, nPars, 0, nPars);

	cout << endl << endl << "Full ratio fit correlation matrix: " << endl << endl;

	for(int i = 0; i < nPars; i++){
	  for(int j = 0; j < nPars; j++){
	    hCorrMatrix_Rmethod->Fill(labels[i], labels[j], (*corrMatrix_Rmethod)[i][j]);
	    // if(i == j) hCorrMatrix_Rmethod->SetBinContent(i, j, 1); // to put 1 in for fixed parameters

	    cout << setiosflags(ios::fixed) << "& " << (*corrMatrix_Rmethod)[i][j] << " "; // also print correlation matrix elements in a nice way to the screen
	  }
	  cout << " \\\\ " <<endl;
	}

	// hCorrMatrix_Rmethod->LabelsOption("h");

	hCorrMatrix_Rmethod->LabelsDeflate();
	hCorrMatrix_Rmethod->GetZaxis()->SetRangeUser(-1,1);
	hCorrMatrix_Rmethod->GetXaxis()->SetLabelSize(0.04);
	hCorrMatrix_Rmethod->GetXaxis()->SetLabelOffset(0.02);
	hCorrMatrix_Rmethod->GetYaxis()->SetLabelSize(0.04);
	hCorrMatrix_Rmethod->SetTitle("");
	hCorrMatrix_Rmethod->Draw("COLZTEXT");

	c_corrMat->SetGrid();
	c_corrMat->SetTopMargin(0.1);
	c_corrMat->SetBottomMargin(0.1);

	c_corrMat->SaveAs(("Images/CorrelationMatrixFullRatioFit" + datasetTagForPlots + ".png").c_str()); // save ratio fit correlation matrix plot
}
else
{
	cout << "No R method fit correlation matrix for fit pass: " << fitPass << endl;
}

/////////////////////////////////////////////////////////////////////////////////////

	TMatrixD* corrMatrix_TMethod = (TMatrixD*) inputFile->Get(Form("topDir/FitPasses/FitPass%i/addedDir/TMethod/TMethodFitCorrMatrix",fitPass));

if(corrMatrix_TMethod != 0)
{	
	TF1* fitFunction_TMethod = (TF1*) ((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass%i/addedDir/TMethod/allTimesAdded_TMethod",fitPass)))->GetFunction("TmethodFitFunc");
	
	const int nPars_TMethod = fitFunction_TMethod->GetNpar();
	char* labels_TMethod[nPars_TMethod];
	for (int parNum = 0; parNum < nPars_TMethod; ++parNum) labels_TMethod[parNum] = Form("%s", fitFunction_TMethod->GetParName(parNum));

	// Draw correlation matrix for TMethod                                                                                                                                                                                         
	TCanvas* c_corrMat_Tmethod = new TCanvas("c_corrMat_Tmethod","Correlation Matrix",800,10,700,700);
	c_corrMat_Tmethod->SetRightMargin(0.14);

	TH2D* hCorrMatrix_TMethod = new TH2D("hCorrMatrix_TMethod","Correlation Matrix", nPars_TMethod, 0, nPars_TMethod, nPars_TMethod, 0, nPars_TMethod);

	cout << endl << endl << "T Method fit correlation matrix: " << endl << endl;

	for(int i = 0; i < nPars_TMethod; i++){
	  for(int j = 0; j < nPars_TMethod; j++){
	    hCorrMatrix_TMethod->Fill(labels_TMethod[i], labels_TMethod[j], (*corrMatrix_TMethod)[i][j]);
	    cout << setiosflags(ios::fixed) << "& " << (*corrMatrix_TMethod)[i][j] << " ";
	  }
	  cout << " \\\\ " <<endl;
	}

	hCorrMatrix_TMethod->LabelsDeflate();
	hCorrMatrix_TMethod->GetZaxis()->SetRangeUser(-1,1);
	hCorrMatrix_TMethod->GetXaxis()->SetLabelSize(0.04);
	hCorrMatrix_TMethod->GetXaxis()->SetLabelOffset(0.02);
	hCorrMatrix_TMethod->GetYaxis()->SetLabelSize(0.04);
	hCorrMatrix_TMethod->SetTitle("");
	hCorrMatrix_TMethod->Draw("COLZTEXT");

	c_corrMat_Tmethod->SetGrid();
	c_corrMat_Tmethod->SetTopMargin(0.1);
	c_corrMat_Tmethod->SetBottomMargin(0.1);

	c_corrMat_Tmethod->SaveAs(("Images/CorrelationMatrixTMethod" + datasetTagForPlots + ".png").c_str()); // save T method fit correlation matrix plot
}
else
{
	cout << "No T method fit correlation matrix for fit pass: " << fitPass << endl;
}
/////////////////////////////////////////////////////////////////////////////////////


	return 1;

}
