// 3-31-20: Macro for plots for comparing calorimeter VW parameters between T method and ratio method fits. Hasn't been used in a while and would need to be dusted off.

#include <iostream>
#include <fstream>
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
#include <sstream>
#include <TNtuple.h>

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;

// macro only set up for VW params at the moment - I could generalize it to more parameters like my other macros, but would need to be careful of any parameter mismatches between the two types of fits

int compareCaloParamsVW(std::string filePath)
{
  // gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  gStyle->SetPadRightMargin(.05);
  gStyle->SetOptFit(111);
  gStyle->SetOptTitle(0);

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

    int Tmethod_VW_Tau_paramNum = -1;


	TF1* TMethodFits[nCalos];
	TF1* fullRatioFits[nCalos];

	  for (int caloNum = 1; caloNum <= nCalos; ++caloNum)
	  {
	  	TF1* TMethodFitFunction = (TF1*) ((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass0/Calos/Calo%d/TMethod/Calo%d_Threshold_Clone", caloNum, caloNum)))->GetFunction("TmethodFitFunc")->Clone();
	    TF1* fullRatioFitFunction = (TF1*) ((TGraphErrors*) inputFile->Get(Form("topDir/FitPasses/FitPass0/Calos/Calo%d/FullRatio/Full_Ratio_TGraph_Calo_%d", caloNum, caloNum)))->GetFunction("fullRatioFitFunc")->Clone();

	  	TMethodFits[caloNum-1] = TMethodFitFunction;
	    fullRatioFits[caloNum-1] = fullRatioFitFunction;

	    if(caloNum ==1){
	    	for (int parNum = 0; parNum < TMethodFitFunction->GetNpar(); ++parNum) if(vw_tau_string.compare(TMethodFitFunction->GetParName(parNum)) == 0) Tmethod_VW_Tau_paramNum = parNum;
	    }
	  }

/////////////////////////////////////////////////////////////////////////////////////

	TGraphErrors* VWtau_Tmethod = new TGraphErrors();
	TGraphErrors* VWamp_Tmethod = new TGraphErrors();
	TGraphErrors* VWphase_Tmethod = new TGraphErrors();

	TGraphErrors* VWtau_Rmethod = new TGraphErrors();
	TGraphErrors* VWamp_Rmethod = new TGraphErrors();
	TGraphErrors* VWphase_Rmethod = new TGraphErrors();

	TGraphErrors* VWtau_Differences = new TGraphErrors();
	TGraphErrors* VWamp_Differences = new TGraphErrors();
	TGraphErrors* VWphase_Differences = new TGraphErrors();

/////////////////////////////////////////////////////////////////////////////////////

	  for (int caloNum = 1; caloNum <= nCalos; ++caloNum)
	  {

	  	double Tmethod_VWtau, Tmethod_VWtau_Err;
	  	double Tmethod_VWamp, Tmethod_VWamp_Err;
	  	double Tmethod_VWphi, Tmethod_VWphi_Err;

	  	double Rmethod_VWtau, Rmethod_VWtau_Err;
	  	double Rmethod_VWamp, Rmethod_VWamp_Err;
	  	double Rmethod_VWphi, Rmethod_VWphi_Err;

	  	Tmethod_VWtau = TMethodFits[caloNum-1]->GetParameter(Tmethod_VW_Tau_paramNum);
	  	Tmethod_VWtau_Err = TMethodFits[caloNum-1]->GetParError(Tmethod_VW_Tau_paramNum);

	  	Tmethod_VWamp = TMethodFits[caloNum-1]->GetParameter(Tmethod_VW_Tau_paramNum+1);
	  	Tmethod_VWamp_Err = TMethodFits[caloNum-1]->GetParError(Tmethod_VW_Tau_paramNum+1);
	  	
	  	Tmethod_VWphi = TMethodFits[caloNum-1]->GetParameter(Tmethod_VW_Tau_paramNum+2);
	  	Tmethod_VWphi_Err = TMethodFits[caloNum-1]->GetParError(Tmethod_VW_Tau_paramNum+2);

	  	Rmethod_VWtau = fullRatioFits[caloNum-1]->GetParameter(Tmethod_VW_Tau_paramNum-2);
	  	Rmethod_VWtau_Err = fullRatioFits[caloNum-1]->GetParError(Tmethod_VW_Tau_paramNum-2);

	  	Rmethod_VWamp = fullRatioFits[caloNum-1]->GetParameter(Tmethod_VW_Tau_paramNum-1);
	  	Rmethod_VWamp_Err = fullRatioFits[caloNum-1]->GetParError(Tmethod_VW_Tau_paramNum-1);

	  	Rmethod_VWphi = fullRatioFits[caloNum-1]->GetParameter(Tmethod_VW_Tau_paramNum);
	  	Rmethod_VWphi_Err = fullRatioFits[caloNum-1]->GetParError(Tmethod_VW_Tau_paramNum);

		normalizePhase(Tmethod_VWphi);
		normalizePhase(Rmethod_VWphi);

/////////////////////////////////////////////////////////////////////////////////////

	  	VWtau_Tmethod->SetPoint(caloNum-1, caloNum, Tmethod_VWtau);
	  	VWtau_Tmethod->SetPointError(caloNum-1, 0, Tmethod_VWtau_Err);

	  	VWamp_Tmethod->SetPoint(caloNum-1, caloNum, Tmethod_VWamp);
	  	VWamp_Tmethod->SetPointError(caloNum-1, 0, Tmethod_VWamp_Err);

	  	VWphase_Tmethod->SetPoint(caloNum-1, caloNum, Tmethod_VWphi);
	  	VWphase_Tmethod->SetPointError(caloNum-1, 0, Tmethod_VWphi_Err);

	  	VWtau_Rmethod->SetPoint(caloNum-1, caloNum, Rmethod_VWtau);
	  	VWtau_Rmethod->SetPointError(caloNum-1, 0, Rmethod_VWtau_Err);

	  	VWamp_Rmethod->SetPoint(caloNum-1, caloNum, Rmethod_VWamp);
	  	VWamp_Rmethod->SetPointError(caloNum-1, 0, Rmethod_VWamp_Err);

	  	VWphase_Rmethod->SetPoint(caloNum-1, caloNum, Rmethod_VWphi);
	  	VWphase_Rmethod->SetPointError(caloNum-1, 0, Rmethod_VWphi_Err);

/////////////////////////////////////////////////////////////////////////////////////

		double tauDiffError, ampDiffError, phiDiffError;

		tauDiffError = sqrt(pow(Tmethod_VWtau_Err,2) + pow(Rmethod_VWtau_Err,2));
		ampDiffError = sqrt(pow(Tmethod_VWamp_Err,2) + pow(Rmethod_VWamp_Err,2));
		phiDiffError = sqrt(pow(Tmethod_VWphi_Err,2) + pow(Rmethod_VWphi_Err,2));

	  	VWtau_Differences->SetPoint(caloNum-1, caloNum, Rmethod_VWtau-Tmethod_VWtau);
	  	VWtau_Differences->SetPointError(caloNum-1, 0, tauDiffError);

	  	VWamp_Differences->SetPoint(caloNum-1, caloNum, Rmethod_VWamp-Tmethod_VWamp);
	  	VWamp_Differences->SetPointError(caloNum-1, 0, ampDiffError);

	  	double phiDiff = Rmethod_VWphi-Tmethod_VWphi;
	  	// normalizePhase(phiDiff); // maybe not the way to make the plot look better

	  	VWphase_Differences->SetPoint(caloNum-1, caloNum, phiDiff);
	  	VWphase_Differences->SetPointError(caloNum-1, 0, phiDiffError);
	  }

/////////////////////////////////////////////////////////////////////////////////////

	VWtau_Tmethod->SetLineColor(2);
	VWamp_Tmethod->SetLineColor(2);
	VWphase_Tmethod->SetLineColor(2);
	VWtau_Tmethod->SetMarkerColor(2);
	VWamp_Tmethod->SetMarkerColor(2);
	VWphase_Tmethod->SetMarkerColor(2);


	VWtau_Tmethod->SetTitle("Tmethod #tau_{VW}");
    VWtau_Tmethod->GetXaxis()->SetTitle("CaloNum");
    VWtau_Tmethod->GetYaxis()->SetTitle("#tau_{VW} (#mus)");

	VWamp_Tmethod->SetTitle("Tmethod A_{VW}");
    VWamp_Tmethod->GetXaxis()->SetTitle("CaloNum");
    VWamp_Tmethod->GetYaxis()->SetTitle("A_{VW}");

	VWphase_Tmethod->SetTitle("Tmethod #phi_{VW}");
    VWphase_Tmethod->GetXaxis()->SetTitle("CaloNum");
    VWphase_Tmethod->GetYaxis()->SetTitle("#phi_{VW}");


	VWtau_Rmethod->SetTitle("Rmethod #tau_{VW}");
    VWtau_Rmethod->GetXaxis()->SetTitle("CaloNum");
    VWtau_Rmethod->GetYaxis()->SetTitle("#tau_{VW} (#mus)");

	VWamp_Rmethod->SetTitle("Rmethod A_{VW}");
    VWamp_Rmethod->GetXaxis()->SetTitle("CaloNum");
    VWamp_Rmethod->GetYaxis()->SetTitle("A_{VW}");

	VWphase_Rmethod->SetTitle("Rmethod #phi_{VW}");
    VWphase_Rmethod->GetXaxis()->SetTitle("CaloNum");
    VWphase_Rmethod->GetYaxis()->SetTitle("#phi_{VW}");


/////////////////////////////////////////////////////////////////////////////////////

	VWtau_Differences->SetTitle("#tau_{VW} Diff");
    VWtau_Differences->GetXaxis()->SetTitle("CaloNum");
    VWtau_Differences->GetYaxis()->SetTitle("#tau_{VW} (R - T) (#mus)");

	VWamp_Differences->SetTitle("A_{VW} Diff");
    VWamp_Differences->GetXaxis()->SetTitle("CaloNum");
    VWamp_Differences->GetYaxis()->SetTitle("A_{VW} (R - T)");

	VWphase_Differences->SetTitle("#phi_{VW} Diff");
    VWphase_Differences->GetXaxis()->SetTitle("CaloNum");
    VWphase_Differences->GetYaxis()->SetTitle("#phi_{VW} (R - T)");


/////////////////////////////////////////////////////////////////////////////////////

  	auto legend_tau = new TLegend(0.6,0.6,.9,0.8);
  	legend_tau->AddEntry(VWtau_Tmethod, "T method", "lp");
  	legend_tau->AddEntry(VWtau_Rmethod, "R method", "lp");
  	legend_tau->SetBorderSize(0);

    auto tau_canv = new TCanvas("tau_canv","tau_canv",200,10,500,400);
    VWtau_Tmethod->Draw("AP");
    VWtau_Rmethod->Draw("PSAME");
 	legend_tau->Draw("SAME");


  	auto legend_amp = new TLegend(0.6,0.6,.9,0.8);
  	legend_amp->AddEntry(VWamp_Tmethod, "T method", "lp");
  	legend_amp->AddEntry(VWamp_Rmethod, "R method", "lp");
  	legend_amp->SetBorderSize(0);

    auto amp_canv = new TCanvas("amp_canv","amp_canv",700,10,500,400);
    VWamp_Tmethod->Draw("AP");
    VWamp_Rmethod->Draw("PSAME");
    legend_amp->Draw("SAME");
    amp_canv->SaveAs("Images/Avw_noTimeRand.png");


  	auto legend_phi = new TLegend(0.6,0.6,.9,0.8);
  	legend_phi->AddEntry(VWphase_Tmethod, "T method", "lp");
  	legend_phi->AddEntry(VWphase_Rmethod, "R method", "lp");
  	legend_phi->SetBorderSize(0);

    auto phi_canv = new TCanvas("phi_canv","phi_canv",1200,10,500,400);
    VWphase_Tmethod->Draw("AP");
    VWphase_Rmethod->Draw("PSAME");
    legend_phi->Draw("SAME");
    phi_canv->SaveAs("Images/Phivw_noTimeRand.png");

/////////////////////////////////////////////////////////////////////////////////////

    auto tauDiff_canv = new TCanvas("tauDiff_canv","tauDiff_canv",200,410,500,400);
    VWtau_Differences->Draw("AP");

    auto ampDiff_canv = new TCanvas("ampDiff_canv","ampDiff_canv",700,410,500,400);
    VWamp_Differences->Draw("AP");
    ampDiff_canv->SaveAs("Images/Avw_Diff_noTimeRand.png");

    auto phiDiff_canv = new TCanvas("phiDiff_canv","phiDiff_canv",1200,410,500,400);
    VWphase_Differences->Draw("AP");
    phiDiff_canv->SaveAs("Images/Phivw_Diff_noTimeRand.png");



  return 1;

}
