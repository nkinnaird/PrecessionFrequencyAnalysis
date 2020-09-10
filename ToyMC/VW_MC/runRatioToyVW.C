R__LOAD_LIBRARY(/cvmfs/gm2.opensciencegrid.org/prod/g-2/gm2util/v9_16_00/slf6.x86_64.e15.prof/lib/libgm2util_blinders.so)

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
#include <TFitResult.h>
#include <TSpectrum.h>
#include <TText.h>
#include <Math/Minimizer.h>

#include "gm2util/blinders/Blinders.hh"

#include "ratioAnalysisDefs.hh"
// #include "ratioAnalysisConfig.hh"
#include "ratioToyAnalysisConfig.hh"
#include "fiveParamFit.hh"
#include "TmethodFit.hh"
#include "ratioFit.hh"
#include "ratioCBOFit.hh"
#include "residualPlots.hh"
#include "copyFile.hh"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////

int runRatioToyVW(std::string filePath)
{
  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);

  // cout << "Remember to comment out VW lifetime stuff in fit classes again." << endl;
  // return 0;

  // pull in input file
  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  // create output file that will hold plots
  TFile* outputFile = new TFile("toyOutputVW.root","RECREATE");

  // make top directory for output file
  auto topDir = outputFile->mkdir("topDir");

/////////////////////////////////////////////////////////////////////////////////////

  double toyFitStart = 30000;
  double toyFitEnd = 650000;

  int totalIters = (*(TVectorD*) inputFile->Get("Iters"))[0]; // total iterations in generated histograms (energy thresholds, etc.)
  TVectorD* histSavedParameters = (TVectorD*) inputFile->Get("topDir/Iter0/SavedParameters/parameterStore");

/////////////////////////////////////////////////////////////////////////////////////
 
  auto toyMCDir = topDir->mkdir("ToyMC");

/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

  blinding::Blinders::fitType ftype = blinding::Blinders::kOmega_a;
  blinding::Blinders* myBlinder = new blinding::Blinders(ftype); // no blinding for ToyMC

  ResidualPlots residualPlotsClass;

  FiveParamFit fiveParamFitToyClass(myBlinder);
  fiveParamFitToyClass.setFitRange(toyFitStart, toyFitEnd);

  TmethodFit TmethodFitClass(myBlinder);
  TmethodFitClass.setFitParamBools(includeChangingCBOFreq, includeChangingVWFreq, TmethodFit_Include2CBO, TmethodFit_IncludeACBO, TmethodFit_IncludePhiCBO, TmethodFit_IncludeVW, TmethodFit_IncludeLM);
  TmethodFitClass.prepareFitFunction();

  TmethodFitClass.setFitRange(toyFitStart, toyFitEnd);
  TH1F* emptyLossHist1 = 0; 
  TH1F* emptyLossHist2 = 0;
  TmethodFitClass.setLostMuonHists(emptyLossHist1, emptyLossHist2);

  RatioFit ratioFitToyClass(myBlinder);
  ratioFitToyClass.setFitRange(toyFitStart, toyFitEnd);

  RatioCBOFit ratioCBOFitClass(myBlinder);
  ratioCBOFitClass.setFitParamBools(includeChangingCBOFreq, includeChangingVWFreq, ratioFit_Include2CBO, ratioFit_IncludeACBO, ratioFit_IncludePhiCBO, ratioFit_IncludeVW, ratioFit_IncludeLM);
  ratioCBOFitClass.prepareFitFunction();
  ratioCBOFitClass.setFitRange(toyFitStart, toyFitEnd);
  ratioCBOFitClass.classForFunc->setHalfPeriod((*histSavedParameters)[0]/2.);
  ratioCBOFitClass.setLostMuonIntegral(emptyLossHist1);

/////////////////////////////////////////////////////////////////////////////////////

  TF1* truthFunc = (TF1*) inputFile->Get("truthFunc"); // could maybe put this outside the file loop, but could cause problems if truth func changes even though it shouldn't
  for (int i = 0; i < truthFunc->GetNpar(); ++i) cout << "Truth func parameter " << i << ": " << truthFunc->GetParameter(i) << endl;

/////////////////////////////////////////////////////////////////////////////////////

  // for (int iter = 0; iter < totalIters; ++iter) // loop over iterations / different random seeds
  // for (int iter = 0; iter < 2; ++iter)
  // for (int iter = 17; iter < 18; ++iter)
  for (int iter = 19; iter < 20; ++iter)
  {
    cout << "Iter: " << iter << endl;
    auto iterDir = toyMCDir->mkdir(Form("Iter%i",iter));
    iterDir->cd();

/////////////////////////////////////////////////////////////////////////////////////

  TF1* iter_truthFunc = (TF1*) inputFile->Get(Form("topDir/Iter%d/SavedParameters/truthFunc", iter));
  for (int i = 0; i < iter_truthFunc->GetNpar(); ++i) cout << "Truth func parameter " << i << ": " << iter_truthFunc->GetParameter(i) << endl;
  iter_truthFunc->Write();

/////////////////////////////////////////////////////////////////////////////////////

      // TH1F* toyFiveParamHist = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_5_Param_Hist_Raw", iter)))->Clone("5FitHist");
      TH1F* toyFiveParamHist = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_5_Param_Hist", iter)))->Clone("5FitHist");

            fiveParamFitToyClass.fiveParameterFitMethod(toyFiveParamHist);
            auto toyHistFiveFit = toyFiveParamHist->GetFunction("fiveParamFit");

              std::cout << "Toy hist 5 param fit p-value is: " << toyHistFiveFit->GetProb() << " chi2/ndf: " << toyHistFiveFit->GetChisquare()/toyHistFiveFit->GetNDF() << " R error: " << toyHistFiveFit->GetParError(3) << " ppm "  << endl;

              toyFiveParamHist->Write();

              residualPlotsClass.makeResidualPlots("FiveParamFit", toyFiveParamHist, "FiveParamFit_residual", toyFitStart, toyFitEnd, true);

/////////////////////////////////////////////////////////////////////////////////////
      
      // TH1F* TmethodFitHist = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_5_Param_Hist_Raw", iter)))->Clone("TmethodFitHist");
      TH1F* TmethodFitHist = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_5_Param_Hist", iter)))->Clone("TmethodFitHist");

            TmethodFitClass.clearFunctionParams();
            TmethodFitClass.setFiveParamFunction(toyHistFiveFit);
            TmethodFitClass.setStartingFitParameters();

            TmethodFitClass.doTMethodFit(TmethodFitHist);
            auto Tmethod_Fit = TmethodFitHist->GetFunction("TmethodFitFunc");

              std::cout << "T method fit p-value is: " << Tmethod_Fit->GetProb() << " chi2/ndf: " << Tmethod_Fit->GetChisquare()/Tmethod_Fit->GetNDF() << " R error: " << Tmethod_Fit->GetParError(3) << " ppm "  << endl;

              TmethodFitHist->Write();

              residualPlotsClass.makeResidualPlots("TmethodFit", TmethodFitHist, "TmethodFit_residual", toyFitStart, toyFitEnd, true);

/////////////////////////////////////////////////////////////////////////////////////

      // TH1F* toyUHist = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_U_Hist", iter)))->Clone();
      // TH1F* toyVHist = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_V_Hist", iter)))->Clone();
      TH1F* toyNumHist = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_Num_Hist", iter)))->Clone();
      TH1F* toyDenomHist = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_Denom_Hist", iter)))->Clone();

            TGraphErrors* toyRatioGraph = ratioFitToyClass.createRatioGraph(toyNumHist, toyDenomHist);
            ratioFitToyClass.ratioFitMethod(toyRatioGraph);

            auto toyRatioFit = toyRatioGraph->GetFunction("threeParamRatioFit");

              std::cout << "Toy Ratio fit p-value is: " << toyRatioFit->GetProb() << " chi2/ndf: " << toyRatioFit->GetChisquare()/toyRatioFit->GetNDF() << " R error: " << toyRatioFit->GetParError(1) << " ppm " << endl;

              toyRatioGraph->SetTitle("Toy Ratio Graph; Time (ns); R (unitless)");
              toyRatioGraph->Write("Toy_Ratio_Graph");

              residualPlotsClass.makeResidualPlots("ThreeParamRatio", toyRatioGraph, "ThreeParamRatio_residual", toyFitStart, toyFitEnd, true);

/////////////////////////////////////////////////////////////////////////////////////

            auto toyRatioCBOGraph = ratioFitToyClass.createRatioGraph(toyNumHist, toyDenomHist);

            ratioCBOFitClass.clearFunctionParams();
            ratioCBOFitClass.setThreeParamRatioFunction(toyRatioFit);

            // fix VW frequency and lifetime from T method fit for check
            if(setStartingParamsFromTMethod) ratioCBOFitClass.setTMethodFunction(Tmethod_Fit, TmethodFitClass.getParamIndicesMap());
            ratioCBOFitClass.setFixedParameters(ratioFit_added_fixedParameters);

            ratioCBOFitClass.setStartingFitParameters();

            ratioCBOFitClass.doFullRatioFit(toyRatioCBOGraph);
            auto toyRatioCBOFit = toyRatioCBOGraph->GetFunction("fullRatioFitFunc");

              std::cout << "R method fit p-value is: " << toyRatioCBOFit->GetProb() << " chi2/ndf: " << toyRatioCBOFit->GetChisquare()/toyRatioCBOFit->GetNDF() << " R error: " << toyRatioCBOFit->GetParError(1) << " ppm "  << endl;

              toyRatioCBOGraph->Write("Toy_Ratio_CBO_Graph");

              residualPlotsClass.makeResidualPlots("FullRatio", toyRatioCBOGraph, "FullRatio_residual", toyFitStart, toyFitEnd, true);

/////////////////////////////////////////////////////////////////////////////////////

            delete toyFiveParamHist;
            delete TmethodFitHist;

            // delete toyUHist;
            // delete toyVHist;
            delete toyNumHist;
            delete toyDenomHist;

} // end loop over iterations

/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

  outputFile->Write();
  delete outputFile;

  return 1;

}
