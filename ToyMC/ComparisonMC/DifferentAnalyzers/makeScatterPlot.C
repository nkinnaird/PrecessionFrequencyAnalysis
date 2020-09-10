#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TVectorD.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveStats.h>
#include <TGraphErrors.h>
#include <TROOT.h>

// #include <plotUtils.hh>

using namespace std;


int makeScatterPlot(string filePath){

  TFile* inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

   gStyle->SetPadRightMargin(0.05);
   gStyle->SetPadTopMargin(0.05);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   TGraph* scatterPlot = (TGraph*) inputFile->Get("Graphs/Graph_Aaron_A_R_vs_David_A_R");
   scatterPlot->GetXaxis()->SetTitle("Aaron A-Method R [ppm]");
   scatterPlot->GetYaxis()->SetTitle("David A-Method R [ppm]");

   auto canv_scatter = new TCanvas("canv_scatter", "canv_scatter", 20, 20, 700, 600);

   scatterPlot->Draw("AP");

   canv_scatter->Update();

   double corr = scatterPlot->GetCorrelationFactor();

    TPaveText *textCorr = new TPaveText(0.22,0.73,0.62,0.88,"NDC");
    textCorr->SetTextSize(0.05);
    textCorr->AddText(Form("Correlation = %1.4f",corr));
    textCorr->Draw("SAME");

    TF1* line = new TF1("straightline", "[0]*x+[1]", -3, 3);
    scatterPlot->Fit(line, "QR");
    line->SetLineColor(2);

    line->Draw("SAME");

    canv_scatter->SaveAs("ScatterPlot.png");


   return 1;

}
