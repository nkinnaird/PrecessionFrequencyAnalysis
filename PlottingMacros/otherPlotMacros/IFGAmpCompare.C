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

int IFGAmpCompare()
{
  SetGm2Style(); // call later to overwrite some things

  TFile *inputFile = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/9d/Gain/IFG/Amplitude/OneBinLater/inFillGainCrystalsScan.root");
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  TMultiGraph* ampMultiGraph_R = (TMultiGraph*) inputFile->Get("topDir/Added/Compare/IFG_Amplitude/Amplitude_Multigraph_R");


	pair<double, double> lineFunc_fullRange {-0.2, 2.2};
	TF1* lineFunc2 = new TF1("lineFunc2", "[0] + (x * [1])", lineFunc_fullRange.first, lineFunc_fullRange.second);


    TList* grlist = 0;
    int graphNum = 0;
    vector<string> legendStrings {"T Method", "R Method"};
    TObject *obj = 0;
    TGraph *gr = 0;


      auto legend_R_amp = new TLegend(0.45,0.2,.75,0.35);

        grlist = ampMultiGraph_R->GetListOfGraphs();
          TIter next_R(grlist);
          graphNum = 0;

          while ((obj = next_R())) {
            gr=(TGraph*)obj;

              gr->SetMarkerStyle(20);
              gr->SetMarkerSize(1.2);

              for (int i = 0; i < gr->GetN(); ++i)
              {
              	double x, y;
      			gr->GetPoint(i, x, y);
      			gr->SetPoint(i, x, y*1000);
              }

              gr->Fit(lineFunc2, "QR");
              gr->GetFunction("lineFunc2")->SetLineColor(graphNum+1);

              double slope = lineFunc2->GetParameter(1);

            legend_R_amp->AddEntry(gr, Form("%s : %2.1f ppb/unit amp", legendStrings.at(graphNum).c_str(), slope), "pl");

            graphNum++;
          } // end loop over graphs

  // SetGm2Style();

        auto Rcanvas_compare_amp = new TCanvas("Rcanvas_compare_amp","Rcanvas_compare_amp",1300,410,500,400);
        ampMultiGraph_R->Draw("APX");

        adjustGraphRanges(ampMultiGraph_R);
        ampMultiGraph_R->GetYaxis()->SetTitle("#DeltaR [ppb]");
        ampMultiGraph_R->GetXaxis()->SetTitle("IFG Amplitude Multiplier");

        gPad->Modified();

        legend_R_amp->SetBorderSize(0);
        legend_R_amp->Draw("SAME");


		  Rcanvas_compare_amp->SaveAs("RvsT_IFG_amp.png");
		  Rcanvas_compare_amp->SaveAs("RvsT_IFG_amp.eps");
		  Rcanvas_compare_amp->SaveAs("RvsT_IFG_amp.pdf");
		  Rcanvas_compare_amp->SaveAs("RvsT_IFG_amp.C");



return 1;

 }
