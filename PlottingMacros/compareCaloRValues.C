// 3-31-20: Macro for plots for comparing calo R values between two different sets of fits. Would need to be dusted off.

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

int compareCaloRValues()
{
  // gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  gStyle->SetPadRightMargin(.05);
  gStyle->SetOptFit(111);
  gStyle->SetOptTitle(0);



  TFile *inputFile1 = TFile::Open("/gm2/data/users/nkinnaird/Ratio/9d-Fresh/SingleIteration/Fits/output-allParams-cbo-w-fixed.root");
  TFile *inputFile2 = TFile::Open("/gm2/data/users/nkinnaird/Ratio/9d-Fresh/TimeRandomize-VW/SingleIter/CaloFits/output-caloFits-withoutVW.root");
   if (inputFile1 == 0 || inputFile2 == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }


	TGraphErrors* Rdiff_graph = new TGraphErrors();



	  for (int caloNum = 1; caloNum <= nCalos; ++caloNum)
	  {
	    TF1* fullRatioFitFunction_1 = (TF1*) ((TGraphErrors*) inputFile1->Get(Form("topDir/FitPasses/FitPass0/Calos/Calo%d/FullRatio/Full_Ratio_TGraph_Calo_%d", caloNum, caloNum)))->GetFunction("fullRatioFitFunc")->Clone();
	    TF1* fullRatioFitFunction_2 = (TF1*) ((TGraphErrors*) inputFile2->Get(Form("topDir/FitPasses/FitPass0/Calos/Calo%d/FullRatio/Full_Ratio_TGraph_Calo_%d", caloNum, caloNum)))->GetFunction("fullRatioFitFunc")->Clone();

	    double R1 = fullRatioFitFunction_1->GetParameter(1);
	    double R1_Err = fullRatioFitFunction_1->GetParError(1);

	    double R2 = fullRatioFitFunction_2->GetParameter(1);
	    double R2_Err = fullRatioFitFunction_2->GetParError(1);

	    double R_diff = R1 - R2;
	    double errDiff = sqrt(R1_Err*R1_Err + R2_Err*R2_Err);

	    cout << "Calo " << caloNum << " with VW: " << R1 << " +- " << R1_Err << " without VW with fVW Rand: " << R2 << " +- " << R2_Err << " Diff: " << R_diff << " err: " << errDiff << endl;

	    Rdiff_graph->SetPoint(caloNum-1, caloNum, R_diff);
	    Rdiff_graph->SetPointError(caloNum-1, 0, errDiff);
	  }


	Rdiff_graph->SetTitle("R Diff");
    Rdiff_graph->GetXaxis()->SetTitle("CaloNum");
    Rdiff_graph->GetYaxis()->SetTitle("R (with VW - without VW with fVW Rand)");

    auto RDiff_canv = new TCanvas("RDiff_canv","RDiff_canv",700,410,500,400);
    Rdiff_graph->Draw("AP");

return 1;

}
