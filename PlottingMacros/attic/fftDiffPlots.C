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
#include <TRandom3.h>
#include <sstream>
#include <TVirtualFFT.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TPaveStats.h>
// #include <TRatioPlot.h> // wait till root v 6_08

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;

void drawLineOnCanv(double freq, TCanvas* canv, string label){

  TLine *line = new TLine(freq, canv->GetUymin(), freq, canv->GetUymax());
  line->SetLineColor(2);
  line->SetLineStyle(2);
  line->Draw();

  TText *t = new TText(freq, canv->GetUymax()+0.02, label.c_str());
  // t->SetTextAlign(22);
  t->SetTextColor(2);
  // t->SetTextFont(43);
  t->SetTextSize(.02);
  t->SetTextAngle(90);
  t->Draw();

  return;
}

void makeFFTCanvas(TH1F* fft){

  TH1F* fft_clone = (TH1F*) fft->Clone("fft_clone");
  // fft_clone->SetStats(false); // doesn't work in TBrowser

  double yMax = fft_clone->GetMaximum(); // scale y axis to get rid of 10^# text in axis label
  fft_clone->Scale(1./pow(10, int(log10(yMax))));

  auto canv = new TCanvas("FFT_Canvas", "FFT Canvas", 200, 10, 1200, 800);
  fft_clone->Draw("hist");
  canv->Update();

  // MHz
  double plot_cbo_freq = 0.37;
  double gm2Freq = 0.23;
  double VWfreq = 2.3;

  drawLineOnCanv(plot_cbo_freq, canv, "cbo");
  drawLineOnCanv(gm2Freq, canv, "g-2");
  drawLineOnCanv(plot_cbo_freq-gm2Freq, canv, "cbo - g-2");
  drawLineOnCanv(plot_cbo_freq+gm2Freq, canv, "cbo + g-2");
  drawLineOnCanv(VWfreq, canv, "VW");

  canv->Write();
  delete canv;

/////////////////////////////////////////////////////////////////////////////////////

  auto zoomed_canv = new TCanvas("FFT_Canvas_Zoomed", "FFT Canvas Zoomed", 200, 10, 1200, 800);
  fft_clone->GetXaxis()->SetRangeUser(0.0, 1.0);
  fft_clone->Draw("hist");
  zoomed_canv->Update();

  drawLineOnCanv(plot_cbo_freq, zoomed_canv, "cbo");
  drawLineOnCanv(gm2Freq, zoomed_canv, "g-2");
  drawLineOnCanv(plot_cbo_freq-gm2Freq, zoomed_canv, "cbo - g-2");
  drawLineOnCanv(plot_cbo_freq+gm2Freq, zoomed_canv, "cbo + g-2");

  zoomed_canv->Write();
  delete zoomed_canv;

/////////////////////////////////////////////////////////////////////////////////////

  delete fft_clone;
  return;
}

int fftDiffPlots(std::string filePath)
{
  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  TFile* outputFile = new TFile("fftDiffPlots.root","RECREATE");
  auto topDir = outputFile->mkdir("topDir");

  auto fitRange_Dir = topDir->mkdir("FitRange");
  auto otherRange_Dir = topDir->mkdir("OtherRange");

  auto fitRange_FFT_dir = fitRange_Dir->mkdir("FFTs");
  auto otherRange_FFT_dir = otherRange_Dir->mkdir("FFTs");

  auto fitRange_diff_dir = fitRange_Dir->mkdir("Diffs");
  auto otherRange_diff_dir = otherRange_Dir->mkdir("Diffs");

/////////////////////////////////////////////////////////////////////////////////////
  // These only get set for the interactive root session (any generated canvases, etc.), but does not apply to the output root file - that comes from .rootlogon.C
  gStyle->SetOptStat(000000);
  gStyle->SetOptTitle(1);
  // gStyle->SetOptFit(0000);
  gStyle->SetOptFit(1111);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(.1);
  gStyle->SetLineColor(1);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
      // Plots over different fit conditions

std::vector<TH1F*> ffts_fit_range;
std::vector<TH1F*> ffts_other_range;



      TNtuple *firstTuple = (TNtuple*)inputFile->Get(Form("topDir/FitPasses/FitPass0/FitConditions0"));
      float tP;
      firstTuple->SetBranchAddress("totalPasses", &tP);
      firstTuple->GetEntry(0);

    for (int fitPass = 0; fitPass < tP; ++fitPass)
    {
      fitRange_FFT_dir->cd();
      ffts_fit_range.push_back((TH1F*) ((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/RatioCBO/residual_fitRange/RatioCBOAddedTimes_fftHist_rescaled", fitPass)))->Clone(Form("FFT_fitRange_FitPass%d", fitPass)));

      otherRange_FFT_dir->cd();
      ffts_other_range.push_back((TH1F*) ((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/addedDir/RatioCBO/residual_otherRange/RatioCBOAddedTimes_fftHist_rescaled", fitPass)))->Clone(Form("FFT_otherRange_FitPass%d", fitPass)));
    }

/////////////////////////////////////////////////////////////////////////////////////

    for (int fitPass = 0; fitPass < tP; ++fitPass)
    {
      fitRange_diff_dir->cd();

      TH1F* diffHist_fitRange = (TH1F*) ffts_fit_range.at(0)->Clone(Form("Diff%d", fitPass));
      diffHist_fitRange->Reset("ICES");

      for(int bin = 1; bin <= diffHist_fitRange->GetNbinsX(); bin++){
        // diffHist_fitRange->SetBinContent(bin, fabs(ffts_fit_range.at(fitPass)->GetBinContent(bin) - ffts_fit_range.back()->GetBinContent(bin)));
        diffHist_fitRange->SetBinContent(bin, ffts_fit_range.at(fitPass)->GetBinContent(bin) - ffts_fit_range.back()->GetBinContent(bin));
      }

      auto canvasDir_fitRange = fitRange_diff_dir->mkdir(Form("Canvases%d", fitPass));
      canvasDir_fitRange->cd();
      makeFFTCanvas(diffHist_fitRange);

/////////////////////////////////////////////////////////////////////////////////////

      otherRange_diff_dir->cd();

      TH1F* diffHist_otherRange = (TH1F*) ffts_other_range.at(0)->Clone(Form("Diff%d", fitPass));
      diffHist_otherRange->Reset("ICES");

      for(int bin = 1; bin <= diffHist_otherRange->GetNbinsX(); bin++){
        diffHist_otherRange->SetBinContent(bin, fabs(ffts_other_range.at(fitPass)->GetBinContent(bin) - ffts_other_range.back()->GetBinContent(bin)));
      }

      auto canvasDir_otherRange = otherRange_diff_dir->mkdir(Form("Canvases%d", fitPass));
      canvasDir_otherRange->cd();
      makeFFTCanvas(diffHist_otherRange);
    }


/////////////////////////////////////////////////////////////////////////////////////

      outputFile->Write();
      delete outputFile;


  return 1;
}
