#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
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
#include <sstream>
#include <TTree.h>
#include <TNtuple.h>
#include <TPaveStats.h>
#include <TVectorD.h>
#include <THStack.h>

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;

int IFGAmpCompareRA()
{
  // gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  SetGm2Style(); // call later to overwrite some things


  TFile *inputFile1 = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/9d/Gain/IFG/Amplitude/OneBinLater/inFillGainCrystalsScan.root");
  // TFile *inputFile2 = TFile::Open("/gm2/app/users/msorbara/gm2DevEuropa_v9_36_00/srcs/gm2analyses/europa/SystematicsTable/9d_A/fit_9d_Awgh_ifgScan.root");
  TFile *inputFile2 = TFile::Open("/gm2/app/users/sweigart/run1-run/run1/analysis/multiplier-gain/9dy-weighted/results_9dy.root");
   if (inputFile1 == 0 || inputFile2 == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  TMultiGraph* ampMultiGraph_R = (TMultiGraph*) inputFile1->Get("topDir/Added/Compare/IFG_Amplitude/Amplitude_Multigraph_R");
  // TGraph* ampGraph_A = (TGraph*) inputFile2->Get("ifgScan/R");
  TGraph* ampGraph_A = (TGraph*) inputFile2->Get("9dy_a_r_amp_i_g");


	pair<double, double> lineFunc_fullRange {-0.2, 2.2};
	TF1* lineFunc2 = new TF1("lineFunc2", "[0] + (x * [1])", lineFunc_fullRange.first, lineFunc_fullRange.second);

/////////////////////////////////////////////////////////////////////////////////////

    TList* grlist = 0;
    vector<string> legendStrings {"T Method", "R Method"};
    TObject *obj = 0;
    TGraph *gr = 0;


      auto legend_R_amp = new TLegend(0.45,0.2,.75,0.35);

        grlist = ampMultiGraph_R->GetListOfGraphs();
          TIter next_R(grlist);

          obj = next_R(); // skip T Method

          obj = next_R(); // R method
          gr=(TGraph*)obj;


              int pointNo = 0;
              while(pointNo < gr->GetN())
              {
                 double x, y;
                 gr->GetPoint(pointNo, x, y);

                 // if(fmod(x,0.25) != 0){
                 //  gr->RemovePoint(pointNo);
                 //  continue;
                 // }

                if(abs(remainder(x,0.1)) > 1e-3){ // use this instead of fmod since there is some precision issues
                  gr->RemovePoint(pointNo);
                  continue;
                 }

                 gr->SetPoint(pointNo, x, y*1000);
                 pointNo++;
              }

              gr->SetMarkerColor(2);
              gr->SetMarkerStyle(20);
              gr->SetMarkerSize(1.2);

              gr->Fit(lineFunc2, "QR");
              gr->GetFunction("lineFunc2")->SetLineColor(2);

              double slope = lineFunc2->GetParameter(1);

            // legend_R_amp->AddEntry(gr, Form("%s : %2.1f ppb/unit amp", legendStrings.at(1).c_str(), slope), "pl");

/////////////////////////////////////////////////////////////////////////////////////

            // A method from Matteo

        // normalizeGraphToNominal(ampGraph_A);

              int pointNoA = 0;
              while(pointNoA < ampGraph_A->GetN())
              {
                 double x, y;
                 ampGraph_A->GetPoint(pointNoA, x, y);

                 // if(fmod(x,0.25) != 0){
                 //  ampGraph_A->RemovePoint(pointNoA);
                 //  continue;
                 // }

                 if(abs(remainder(x,0.1)) > 1e-3){ // use this instead of fmod since there is some precision issues
                  ampGraph_A->RemovePoint(pointNoA);
                  continue;
                 }

                 pointNoA++;
              }

              // for (int i = 0; i < ampGraph_A->GetN(); ++i)
              // {
              //    double x, y;
              //    ampGraph_A->GetPoint(i, x, y);
              //    ampGraph_A->SetPoint(i, x, y*1000);
              // }

  ampGraph_A->SetMarkerStyle(20);
  ampGraph_A->SetMarkerSize(1.2);


              ampGraph_A->Fit(lineFunc2, "QR");
              ampGraph_A->GetFunction("lineFunc2")->SetLineColor(1);

              double slopeA = lineFunc2->GetParameter(1);

/////////////////////////////////////////////////////////////////////////////////////

  TF1* davidsFit = (TF1*) inputFile2->Get("9dy_a_r_amp_i_f");

  double davidFitSlope = davidsFit->GetParameter(1);
  cout << " David fit slope = " << davidFitSlope << endl;


/////////////////////////////////////////////////////////////////////////////////////

            // legend_R_amp->AddEntry(ampGraph_A, Form("%s : %2.1f ppb/unit amp", "A Method", slopeA), "pl");
            legend_R_amp->AddEntry(ampGraph_A, Form("%s : %2.1f ppb/unit amp", "A Method", davidFitSlope), "pl");
            legend_R_amp->AddEntry(gr, Form("%s : %2.1f ppb/unit amp", legendStrings.at(1).c_str(), slope), "pl");

/////////////////////////////////////////////////////////////////////////////////////

  // SetGm2Style();

        auto Rcanvas_compare_amp = new TCanvas("Rcanvas_compare_amp","Rcanvas_compare_amp",1300,410,500,400);

        ampGraph_A->Draw("APX");
        gr->Draw("PXSAME");

        adjustGraphRanges(ampGraph_A);
        ampGraph_A->GetYaxis()->SetTitle("#DeltaR [ppb]");
        ampGraph_A->GetXaxis()->SetTitle("IFG Amplitude Multiplier");

        gPad->Modified();

        legend_R_amp->SetBorderSize(0);
        legend_R_amp->Draw("SAME");


		  Rcanvas_compare_amp->SaveAs("RvsA_IFG_amp.png");
		  Rcanvas_compare_amp->SaveAs("RvsA_IFG_amp.eps");
		  Rcanvas_compare_amp->SaveAs("RvsA_IFG_amp.pdf");
		  Rcanvas_compare_amp->SaveAs("RvsA_IFG_amp.C");



return 1;

 }
