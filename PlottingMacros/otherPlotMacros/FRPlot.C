
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


int FRPlot()
{

  TFile *firstFile = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/Misc/OtherImages/FRPlot/addedHists-9d-5nsBins-noFRrand.root");
  TFile *secondFile = TFile::Open("/gm2/data/users/nkinnaird/Ratio/FinalProductions/Misc/OtherImages/FRPlot/addedHists-9d-5nsBins-withFRrand.root");

   if (firstFile == 0 || secondFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  SetGm2Style();

/////////////////////////////////////////////////////////////////////////////////////


  auto timesHist_calo_first = (TH1F*) firstFile->Get("topDir/Iter0/Calos/Calo1/Calo1_Times_E_Threshold");
  auto timesHist_calo_second = (TH1F*) secondFile->Get("topDir/Iter0/Calos/Calo1/Calo1_Times_E_Threshold");

  nsTOus(timesHist_calo_first, "Time [#mus]");
  nsTOus(timesHist_calo_second, "Time [#mus]");


  auto myCanvas = new TCanvas("myCanvas","myCanvas",200,10,1100,600);

  timesHist_calo_first->SetTitle("Calo 1 Times");
  timesHist_calo_first->GetXaxis()->SetRangeUser(28, 44);
  timesHist_calo_first->GetYaxis()->SetTitle("Clusters / 5 ns");
  timesHist_calo_first->Draw();

  timesHist_calo_second->SetLineColor(2);
  timesHist_calo_second->Draw("SAME");


  TLegend* legend = new TLegend(0.6,0.73,0.8,0.88);
  legend->SetBorderSize(0);  // no border
  legend->AddEntry(timesHist_calo_first, "Without time randomization", "F");
  legend->AddEntry(timesHist_calo_second, "With time randomization", "F");
  legend->Draw("SAME");

  // myCanvas->Update();

  myCanvas->SaveAs("fastrotationSimple.png");
  myCanvas->SaveAs("fastrotationSimple.eps");
  myCanvas->SaveAs("fastrotationSimple.pdf");
  myCanvas->SaveAs("fastrotationSimple.C");


/////////////////////////////////////////////////////////////////////////////////////

  return 1;

}
