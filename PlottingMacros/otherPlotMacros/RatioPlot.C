
#include <iostream>
#include <fstream>
#include <string>
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
#include <sstream>
#include <TPaveStats.h>

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;


int RatioPlot()
{

  TFile *inputFile = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/Endgame/SingleIter/SingleFit-NewRange/output-Endgame-SingleIter-NewRange.root");

   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  SetGm2Style();

/////////////////////////////////////////////////////////////////////////////////////


  auto ratioGraph = (TGraph*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/FullRatio/Added_Times_Full_Ratio_Graph");
  nsTOus(ratioGraph, "Time [#mus]");


  auto myCanvas = new TCanvas("myCanvas","myCanvas",200,10,1100,600);

  ratioGraph->GetXaxis()->SetRangeUser(30, 60);
  ratioGraph->GetYaxis()->SetRangeUser(-0.6, 0.6);
  ratioGraph->GetYaxis()->SetTitle("Ratio");
  ratioGraph->Draw();


  myCanvas->SaveAs("ratioPlot.png");
  myCanvas->SaveAs("ratioPlot.eps");
  myCanvas->SaveAs("ratioPlot.pdf");
  myCanvas->SaveAs("ratioPlot.C");


/////////////////////////////////////////////////////////////////////////////////////

  return 1;

}
