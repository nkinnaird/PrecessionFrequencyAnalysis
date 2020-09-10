// note this module assumes an inpute file from fits to a toy mc histFromPairs.root file with 4 iterations contained within it - truth, 'measured', 'measuree' - true pileup, and 'measured' - shadow pileup

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
#include <TVectorD.h>
// #include <TRatioPlot.h> // wait till root v 6_08

#include "../../ratioMacroHeaders/ratioAnalysisDefs.hh"
#include "../../ratioMacroHeaders/copyFile.hh"


int PileupPlots(std::string filePath)
{
    // create output file that will hold plots
  TFile* outputFile = new TFile("pileupMCPlots.root","RECREATE");

  // make top directory for output file
  auto topDir = outputFile->mkdir("topDir");
  topDir->cd();

  TFile *inputFile  = TFile::Open(filePath.c_str());
  if (inputFile == 0) {
     printf("Error: cannot open file\n");
     return 0;
  }

/////////////////////////////////////////////////////////////////////////////////////

  // CopyFile(filePath.c_str());


TIter next(inputFile->GetListOfKeys());
TKey *key;

string genFileName;
while ((key = (TKey*)next())) {
  string name = key->GetName();
  if(name.size() > 5 && name.substr(name.size()-5) == ".root"){
    genFileName = name;
    cout << "gen file name is: " << genFileName << endl;
    break;
  }
}


/////////////////////////////////////////////////////////////////////////////////////

  auto pileupDir = topDir->mkdir("BasePileupPlots");


  auto compareDir = topDir->mkdir("ComparePlots");
  auto fiveFitDir = compareDir->mkdir("FiveParam");


  auto ratioFitDir = compareDir->mkdir("Ratio");
  ratioFitDir->cd();

  TGraph* ratio_A_diff_graph = new TGraph();
  ratio_A_diff_graph->SetName("ratio_A_diff_graph");
  ratio_A_diff_graph->SetTitle("ratio_A_diff_graph; PointNo; A(Iter# - Iter0)");

    TGraph* ratio_A_err_diff_graph = new TGraph();
    ratio_A_err_diff_graph->SetName("ratio_A_err_diff_graph");
    ratio_A_err_diff_graph->SetTitle("ratio_A_err_diff_graph; PointNo; A_err(Iter# - Iter0)");

  TGraph* ratio_R_diff_graph = new TGraph();
  ratio_R_diff_graph->SetName("ratio_R_diff_graph");
  ratio_R_diff_graph->SetTitle("ratio_R_diff_graph; PointNo; R(Iter# - Iter0)");

    TGraph* ratio_R_err_diff_graph = new TGraph();
    ratio_R_err_diff_graph->SetName("ratio_R_err_diff_graph");
    ratio_R_err_diff_graph->SetTitle("ratio_R_err_diff_graph; PointNo; R_err(Iter# - Iter0)");

  TGraph* ratio_Phi_diff_graph = new TGraph();
  ratio_Phi_diff_graph->SetName("ratio_Phi_diff_graph");
  ratio_Phi_diff_graph->SetTitle("ratio_Phi_diff_graph; PointNo; #phi(Iter# - Iter0)");

    TGraph* ratio_Phi_err_diff_graph = new TGraph();
    ratio_Phi_err_diff_graph->SetName("ratio_Phi_err_diff_graph");
    ratio_Phi_err_diff_graph->SetTitle("ratio_Phi_err_diff_graph; PointNo; #phi_err(Iter# - Iter0)");
/////////////////////////////////////////////////////////////////////////////////////

    // five param fit comparisons just because I can

    TH1F* true_fiveParamHist = (TH1F*) inputFile->Get("topDir/ToyMCIter0/Hist/Iter0/Toy_5_Param_Hist");
      TF1* true_fittedFiveParamFunction = (TF1*) true_fiveParamHist->GetFunction("fiveParamFit");   
    TH1F* measured_fiveParamHist = (TH1F*) inputFile->Get("topDir/ToyMCIter1/Hist/Iter1/Toy_5_Param_Hist");
      TF1* measured_fittedFiveParamFunction = (TF1*) measured_fiveParamHist->GetFunction("fiveParamFit");   
    TH1F* meas_truePileup_fiveParamHist = (TH1F*) inputFile->Get("topDir/ToyMCIter2/Hist/Iter2/Toy_5_Param_Hist");
      TF1* meas_truePileup_fittedFiveParamFunction = (TF1*) meas_truePileup_fiveParamHist->GetFunction("fiveParamFit");   
    TH1F* meas_shadowPileup_fiveParamHist = (TH1F*) inputFile->Get("topDir/ToyMCIter3/Hist/Iter3/Toy_5_Param_Hist");
      TF1* meas_shadowPileup_fittedFiveParamFunction = (TF1*) meas_shadowPileup_fiveParamHist->GetFunction("fiveParamFit");   

    double five_true_N = true_fittedFiveParamFunction->GetParameter(0);
      double five_true_N_err = true_fittedFiveParamFunction->GetParError(0);
    double five_true_tau = true_fittedFiveParamFunction->GetParameter(1);
      double five_true_tau_err = true_fittedFiveParamFunction->GetParError(1);
    double five_true_A = true_fittedFiveParamFunction->GetParameter(2);
      double five_true_A_err = true_fittedFiveParamFunction->GetParError(2);
    double five_true_R = true_fittedFiveParamFunction->GetParameter(3);
      double five_true_R_err = true_fittedFiveParamFunction->GetParError(3);
    double five_true_phi = true_fittedFiveParamFunction->GetParameter(4);
      double five_true_phi_err = true_fittedFiveParamFunction->GetParError(4);

    double five_measured_N = measured_fittedFiveParamFunction->GetParameter(0);
      double five_measured_N_err = measured_fittedFiveParamFunction->GetParError(0);
    double five_measured_tau = measured_fittedFiveParamFunction->GetParameter(1);
      double five_measured_tau_err = measured_fittedFiveParamFunction->GetParError(1);
    double five_measured_A = measured_fittedFiveParamFunction->GetParameter(2);
      double five_measured_A_err = measured_fittedFiveParamFunction->GetParError(2);
    double five_measured_R = measured_fittedFiveParamFunction->GetParameter(3);
      double five_measured_R_err = measured_fittedFiveParamFunction->GetParError(3);
    double five_measured_phi = measured_fittedFiveParamFunction->GetParameter(4);
      double five_measured_phi_err = measured_fittedFiveParamFunction->GetParError(4);

    double five_meas_truePileup_N = meas_truePileup_fittedFiveParamFunction->GetParameter(0);
      double five_meas_truePileup_N_err = meas_truePileup_fittedFiveParamFunction->GetParError(0);
    double five_meas_truePileup_tau = meas_truePileup_fittedFiveParamFunction->GetParameter(1);
      double five_meas_truePileup_tau_err = meas_truePileup_fittedFiveParamFunction->GetParError(1);
    double five_meas_truePileup_A = meas_truePileup_fittedFiveParamFunction->GetParameter(2);
      double five_meas_truePileup_A_err = meas_truePileup_fittedFiveParamFunction->GetParError(2);
    double five_meas_truePileup_R = meas_truePileup_fittedFiveParamFunction->GetParameter(3);
      double five_meas_truePileup_R_err = meas_truePileup_fittedFiveParamFunction->GetParError(3);
    double five_meas_truePileup_phi = meas_truePileup_fittedFiveParamFunction->GetParameter(4);
      double five_meas_truePileup_phi_err = meas_truePileup_fittedFiveParamFunction->GetParError(4);

    double five_meas_shadowPileup_N = meas_shadowPileup_fittedFiveParamFunction->GetParameter(0);
      double five_meas_shadowPileup_N_err = meas_shadowPileup_fittedFiveParamFunction->GetParError(0);
    double five_meas_shadowPileup_tau = meas_shadowPileup_fittedFiveParamFunction->GetParameter(1);
      double five_meas_shadowPileup_tau_err = meas_shadowPileup_fittedFiveParamFunction->GetParError(1);
    double five_meas_shadowPileup_A = meas_shadowPileup_fittedFiveParamFunction->GetParameter(2);
      double five_meas_shadowPileup_A_err = meas_shadowPileup_fittedFiveParamFunction->GetParError(2);
    double five_meas_shadowPileup_R = meas_shadowPileup_fittedFiveParamFunction->GetParameter(3);
      double five_meas_shadowPileup_R_err = meas_shadowPileup_fittedFiveParamFunction->GetParError(3);
    double five_meas_shadowPileup_phi = meas_shadowPileup_fittedFiveParamFunction->GetParameter(4);
      double five_meas_shadowPileup_phi_err = meas_shadowPileup_fittedFiveParamFunction->GetParError(4);


/////////////////////////////////////////////////////////////////////////////////////

    // 3 param ratio fit comparisons - don't need to go higher for these studies

    TF1* true_fittedRatioFunction = (TF1*) ((TGraph*) inputFile->Get("topDir/ToyMCIter0/Ratio/Iter0/Toy_Ratio_Graph"))->GetFunction("ratioFitIntegral");
    TF1* measured_fittedRatioFunction = (TF1*) ((TGraph*) inputFile->Get("topDir/ToyMCIter1/Ratio/Iter1/Toy_Ratio_Graph"))->GetFunction("ratioFitIntegral");
    TF1* meas_truePileup_fittedRatioFunction = (TF1*) ((TGraph*) inputFile->Get("topDir/ToyMCIter2/Ratio/Iter2/Toy_Ratio_Graph"))->GetFunction("ratioFitIntegral");
    TF1* meas_shadowPileup_fittedRatioFunction = (TF1*) ((TGraph*) inputFile->Get("topDir/ToyMCIter3/Ratio/Iter3/Toy_Ratio_Graph"))->GetFunction("ratioFitIntegral");

    double ratio_true_A = true_fittedRatioFunction->GetParameter(0);
      double ratio_true_A_err = true_fittedRatioFunction->GetParError(0);
    double ratio_true_R = true_fittedRatioFunction->GetParameter(1);
      double ratio_true_R_err = true_fittedRatioFunction->GetParError(1);
    double ratio_true_phi = true_fittedRatioFunction->GetParameter(2);
      double ratio_true_phi_err = true_fittedRatioFunction->GetParError(2);

    double ratio_measured_A = measured_fittedRatioFunction->GetParameter(0);
      double ratio_measured_A_err = measured_fittedRatioFunction->GetParError(0);
    double ratio_measured_R = measured_fittedRatioFunction->GetParameter(1);
      double ratio_measured_R_err = measured_fittedRatioFunction->GetParError(1);
    double ratio_measured_phi = measured_fittedRatioFunction->GetParameter(2);
      double ratio_measured_phi_err = measured_fittedRatioFunction->GetParError(2);

    double ratio_meas_truePileup_A = meas_truePileup_fittedRatioFunction->GetParameter(0);
      double ratio_meas_truePileup_A_err = meas_truePileup_fittedRatioFunction->GetParError(0);
    double ratio_meas_truePileup_R = meas_truePileup_fittedRatioFunction->GetParameter(1);
      double ratio_meas_truePileup_R_err = meas_truePileup_fittedRatioFunction->GetParError(1);
    double ratio_meas_truePileup_phi = meas_truePileup_fittedRatioFunction->GetParameter(2);
      double ratio_meas_truePileup_phi_err = meas_truePileup_fittedRatioFunction->GetParError(2);

    double ratio_meas_shadowPileup_A = meas_shadowPileup_fittedRatioFunction->GetParameter(0);
      double ratio_meas_shadowPileup_A_err = meas_shadowPileup_fittedRatioFunction->GetParError(0);
    double ratio_meas_shadowPileup_R = meas_shadowPileup_fittedRatioFunction->GetParameter(1);
      double ratio_meas_shadowPileup_R_err = meas_shadowPileup_fittedRatioFunction->GetParError(1);
    double ratio_meas_shadowPileup_phi = meas_shadowPileup_fittedRatioFunction->GetParameter(2);
      double ratio_meas_shadowPileup_phi_err = meas_shadowPileup_fittedRatioFunction->GetParError(2);

/////////////////////////////////////////////////////////////////////////////////////

      ratio_A_diff_graph->SetPoint(0, 0, ratio_true_A - ratio_true_A);
      ratio_A_diff_graph->SetPoint(1, 1, ratio_measured_A - ratio_true_A);
      ratio_A_diff_graph->SetPoint(2, 2, ratio_meas_truePileup_A - ratio_true_A);
      ratio_A_diff_graph->SetPoint(3, 3, ratio_meas_shadowPileup_A - ratio_true_A);
        ratio_A_diff_graph->Write();

        ratio_A_err_diff_graph->SetPoint(0, 0, ratio_true_A_err - ratio_true_A_err);
        ratio_A_err_diff_graph->SetPoint(1, 1, ratio_measured_A_err - ratio_true_A_err);
        ratio_A_err_diff_graph->SetPoint(2, 2, ratio_meas_truePileup_A_err - ratio_true_A_err);
        ratio_A_err_diff_graph->SetPoint(3, 3, ratio_meas_shadowPileup_A_err - ratio_true_A_err);
          ratio_A_err_diff_graph->Write();

      ratio_R_diff_graph->SetPoint(0, 0, ratio_true_R - ratio_true_R);
      ratio_R_diff_graph->SetPoint(1, 1, ratio_measured_R - ratio_true_R);
      ratio_R_diff_graph->SetPoint(2, 2, ratio_meas_truePileup_R - ratio_true_R);
      ratio_R_diff_graph->SetPoint(3, 3, ratio_meas_shadowPileup_R - ratio_true_R);
        ratio_R_diff_graph->Write();

        ratio_R_err_diff_graph->SetPoint(0, 0, ratio_true_R_err - ratio_true_R_err);
        ratio_R_err_diff_graph->SetPoint(1, 1, ratio_measured_R_err - ratio_true_R_err);
        ratio_R_err_diff_graph->SetPoint(2, 2, ratio_meas_truePileup_R_err - ratio_true_R_err);
        ratio_R_err_diff_graph->SetPoint(3, 3, ratio_meas_shadowPileup_R_err - ratio_true_R_err);
          ratio_R_err_diff_graph->Write();

      ratio_Phi_diff_graph->SetPoint(0, 0, ratio_true_phi - ratio_true_phi);
      ratio_Phi_diff_graph->SetPoint(1, 1, ratio_measured_phi - ratio_true_phi);
      ratio_Phi_diff_graph->SetPoint(2, 2, ratio_meas_truePileup_phi - ratio_true_phi);
      ratio_Phi_diff_graph->SetPoint(3, 3, ratio_meas_shadowPileup_phi - ratio_true_phi);
        ratio_Phi_diff_graph->Write();

        ratio_Phi_err_diff_graph->SetPoint(0, 0, ratio_true_phi_err - ratio_true_phi_err);
        ratio_Phi_err_diff_graph->SetPoint(1, 1, ratio_measured_phi_err - ratio_true_phi_err);
        ratio_Phi_err_diff_graph->SetPoint(2, 2, ratio_meas_truePileup_phi_err - ratio_true_phi_err);
        ratio_Phi_err_diff_graph->SetPoint(3, 3, ratio_meas_shadowPileup_phi_err - ratio_true_phi_err);
          ratio_Phi_err_diff_graph->Write();


/////////////////////////////////////////////////////////////////////////////////////


    auto dT_dir = topDir->mkdir("DeltaT");
    dT_dir->cd();

    TH1F* deltaT_true_all = (TH1F*) inputFile->Get((genFileName + "/topDir/Iter0/Time_Between_Hits_All").c_str())->Clone();
      TF1* exp_function = new TF1("exp_func", "[0]*exp([1]*x)+[2]");
      exp_function->SetLineColor(2);
      deltaT_true_all->Fit("exp_func");

/////////////////////////////////////////////////////////////////////////////////////

  outputFile->Write();
  delete outputFile;


  return 1;
 }
