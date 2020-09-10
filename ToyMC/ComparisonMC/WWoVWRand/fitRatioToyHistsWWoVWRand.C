R__LOAD_LIBRARY(/cvmfs/gm2.opensciencegrid.org/prod/g-2/gm2util/v9_21_06/slf6.x86_64.e15.prof/lib/libgm2util_blinders.so)

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

#include "gm2util/blinders/Blinders.hh"

#include "ratioAnalysisDefs.hh"
#include "ratioAnalysisConfig.hh"
#include "fiveParamFit.hh"
#include "threeParameterRatioFit.hh"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////

int fitRatioToyHistsWWoVWRand(std::string filePath)
{
  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen


  // pull in input file
  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  // create output file that will hold plots
  TFile* outputFile = new TFile("toyOutputCompareWWoVWRand.root","RECREATE");

  // make top directory for output file
  auto topDir = outputFile->mkdir("topDir");

/////////////////////////////////////////////////////////////////////////////////////

  TF1* truthFunc = (TF1*) inputFile->Get("truthFunc"); // could maybe put this outside the file loop, but could cause problems if truth func changes even though it shouldn't
  truthFunc->Write();

  for (int i = 0; i < truthFunc->GetNpar(); ++i) cout << "Truth func parameter " << i << ": " << truthFunc->GetParameter(i) << endl;

/////////////////////////////////////////////////////////////////////////////////////

  setDataset(1, 12); // set dataset case to 60h - for 'starting' parameters

  double toyFitStart = 30200;
  double toyFitEnd = 650000;

  int totalIters = (*(TVectorD*) inputFile->Get("Iters"))[0]; // total iterations in generated histograms (energy thresholds, etc.)

/////////////////////////////////////////////////////////////////////////////////////
 
  auto toyMC_Dir = topDir->mkdir("ToyMC");
  auto toyHist_Dir = toyMC_Dir->mkdir("Hist");
  auto toyRatio_Dir = toyMC_Dir->mkdir("Ratio");

  toyHist_Dir->cd();

  auto toy5ParamHistPVal = new TH1F("toy5ParamHistPVal", "toy5ParamHistPVal; p value; Events", 100, 0, 1);
  auto toy5ParamHistChi2NDF = new TH1F("toy5ParamHistChi2NDF", "toy5ParamHistChi2NDF; #chi^{2}/ndf; Events", 50, .9, 1.1);
  auto toyHistFiveR = new TH1F("toyHistFiveR", "Toy MC Five Param Fit R Vs Random Seed; R (ppm); Events", 100, -5, 5);

  TNtuple* fiveParameterFitValues = new TNtuple("fiveParameterFitValues", "fiveParameterFitValues", "N:N_err:tau:tau_err:A:A_err:R:R_err:phi:phi_err");

  toyRatio_Dir->cd();

  auto toyRatioPVal = new TH1F("toyRatioPVal", "toyRatioPVal; p value; Events", 100, 0, 1);
  auto toyRatioChi2NDF = new TH1F("toyRatioChi2NDF", "toyRatioChi2NDF; #chi^{2}/ndf; Events", 50, .9, 1.1);
  auto toyRatioR = new TH1F("toyRatioR", "Toy Ratio Fit R Vs Random Seed; R (ppm); Events", 100, -5, 5);

  TNtuple* ratioFitValues = new TNtuple("ratioFitValues", "ratioFitValues_VW", "A:A_err:R:R_err:phi:phi_err");

/////////////////////////////////////////////////////////////////////////////////////

  // construct VW randomization histograms directories

  auto toyMC_VW_Dir = topDir->mkdir("ToyMC_VW");
  auto toyHist_VW_Dir = toyMC_VW_Dir->mkdir("Hist");
  auto toyRatio_VW_Dir = toyMC_VW_Dir->mkdir("Ratio");

  toyHist_VW_Dir->cd();

  auto toy5ParamHistPVal_VW = new TH1F("toy5ParamHistPVal_VW", "toy5ParamHistPVal_VW; p value; Events", 100, 0, 1);
  auto toy5ParamHistChi2NDF_VW = new TH1F("toy5ParamHistChi2NDF_VW", "toy5ParamHistChi2NDF_VW; #chi^{2}/ndf; Events", 50, .9, 1.1);
  auto toyHistFiveR_VW = new TH1F("toyHistFiveR_VW", "Toy MC Five Param Fit R Vs Random Seed; R (ppm); Events", 100, -5, 5);

  TNtuple* fiveParameterFitValues_VW = new TNtuple("fiveParameterFitValues_VW", "fiveParameterFitValues", "N:N_err:tau:tau_err:A:A_err:R:R_err:phi:phi_err");

  toyRatio_VW_Dir->cd();

  auto toyRatioPVal_VW = new TH1F("toyRatioPVal_VW", "toyRatioPVal_VW; p value; Events", 100, 0, 1);
  auto toyRatioChi2NDF_VW = new TH1F("toyRatioChi2NDF_VW", "toyRatioChi2NDF_VW; #chi^{2}/ndf; Events", 50, .9, 1.1);
  auto toyRatioR_VW = new TH1F("toyRatioR_VW", "Toy Ratio Fit R Vs Random Seed; R (ppm); Events", 100, -5, 5);

  TNtuple* ratioFitValues_VW = new TNtuple("ratioFitValues_VW", "ratioFitValues_VW", "A:A_err:R:R_err:phi:phi_err");

/////////////////////////////////////////////////////////////////////////////////////

  blinding::Blinders::fitType ftype = blinding::Blinders::kOmega_a;
  blinding::Blinders* myBlinder = new blinding::Blinders(ftype); // no blinding for ToyMC

  FiveParamFit fiveParamFitToyClass(myBlinder);
  fiveParamFitToyClass.setFitRange(toyFitStart, toyFitEnd);

  ThreeParameterRatioFit ratioFitToyClass(myBlinder);
  ratioFitToyClass.setFitRange(toyFitStart, toyFitEnd);

/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

  for (int iter = 0; iter < totalIters; ++iter) // loop over iterations / different random seeds
  {
    cout << "Iter: " << iter << endl;

      TH1F* toyFiveParamHist = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_5_Param_Hist", iter)))->Clone();

            fiveParamFitToyClass.fiveParameterFitMethod(toyFiveParamHist);

            auto toyHistFiveFit = toyFiveParamHist->GetFunction("fiveParamFit");

            double histPValue = toyHistFiveFit->GetProb();
            toy5ParamHistPVal->Fill(histPValue);
            toy5ParamHistChi2NDF->Fill(toyHistFiveFit->GetChisquare()/toyHistFiveFit->GetNDF());

            toyHistFiveR->Fill(toyHistFiveFit->GetParameter(3));

            if(iter == 0){
              std::cout << "Toy hist 5 param fit p-value is: " << histPValue << " chi2/ndf: " << toyHistFiveFit->GetChisquare()/toyHistFiveFit->GetNDF() << std::endl;
              cout << "Toy hist 5 param fit R error: " << toyHistFiveFit->GetParError(3) << " ppm "  << endl;
            }

            fiveParameterFitValues->Fill(toyHistFiveFit->GetParameter(0), toyHistFiveFit->GetParError(0),
                                         toyHistFiveFit->GetParameter(1), toyHistFiveFit->GetParError(1),
                                         toyHistFiveFit->GetParameter(2), toyHistFiveFit->GetParError(2),
                                         toyHistFiveFit->GetParameter(3), toyHistFiveFit->GetParError(3),
                                         toyHistFiveFit->GetParameter(4), toyHistFiveFit->GetParError(4));

            delete toyFiveParamHist;

/////////////////////////////////////////////////////////////////////////////////////

            // now do 5 parameter stuff for histograms with VW randomization

      TH1F* toyFiveParamHist_VW = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_5_Param_Hist_withVWRand", iter)))->Clone();

            fiveParamFitToyClass.fiveParameterFitMethod(toyFiveParamHist_VW);

            auto toyHistFiveFit_VW = toyFiveParamHist_VW->GetFunction("fiveParamFit");

            double histPValue_VW = toyHistFiveFit_VW->GetProb();
            toy5ParamHistPVal_VW->Fill(histPValue_VW);
            toy5ParamHistChi2NDF_VW->Fill(toyHistFiveFit_VW->GetChisquare()/toyHistFiveFit_VW->GetNDF());

            toyHistFiveR_VW->Fill(toyHistFiveFit_VW->GetParameter(3));

            if(iter == 0){
              std::cout << "Toy hist 5 param fit p-value is: " << histPValue_VW << " chi2/ndf: " << toyHistFiveFit_VW->GetChisquare()/toyHistFiveFit_VW->GetNDF() << std::endl;
              cout << "Toy hist 5 param fit R error: " << toyHistFiveFit_VW->GetParError(3) << " ppm "  << endl;
            }

            fiveParameterFitValues_VW->Fill(toyHistFiveFit_VW->GetParameter(0), toyHistFiveFit_VW->GetParError(0),
                                            toyHistFiveFit_VW->GetParameter(1), toyHistFiveFit_VW->GetParError(1),
                                            toyHistFiveFit_VW->GetParameter(2), toyHistFiveFit_VW->GetParError(2),
                                            toyHistFiveFit_VW->GetParameter(3), toyHistFiveFit_VW->GetParError(3),
                                            toyHistFiveFit_VW->GetParameter(4), toyHistFiveFit_VW->GetParError(4));

            delete toyFiveParamHist_VW;

/////////////////////////////////////////////////////////////////////////////////////

      TH1F* toyUHist = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_U_Hist", iter)))->Clone();
      TH1F* toyVHist = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_V_Hist", iter)))->Clone();

      TH1F* toyNumHist = (TH1F*) toyVHist->Clone("numeratorHist");
      toyNumHist->Add(toyUHist, -1);

      TH1F* toyDenomHist = (TH1F*) toyVHist->Clone("denominatorHist");
      toyDenomHist->Add(toyUHist);

            toyRatio_Dir->cd();

            TGraphErrors* toyRatioGraph = ratioFitToyClass.createRatioGraph(toyNumHist, toyDenomHist);
            ratioFitToyClass.fitMethod(toyRatioGraph);

            auto toyRatioFit = toyRatioGraph->GetFunction("threeParamRatioFit");

            if(iter == 0)
            {
              std::cout << "Toy Ratio fit p-value is: " << toyRatioFit->GetProb() << " chi2/ndf: " << toyRatioFit->GetChisquare()/toyRatioFit->GetNDF() << std::endl;
              cout << "Toy Ratio fit R error: " << toyRatioFit->GetParError(1) << " ppm " << endl;
             }

            toyRatioPVal->Fill(toyRatioFit->GetProb());
            toyRatioChi2NDF->Fill(toyRatioFit->GetChisquare()/toyRatioFit->GetNDF());

            toyRatioR->Fill(toyRatioFit->GetParameter(1));

            ratioFitValues->Fill(toyRatioFit->GetParameter(0), toyRatioFit->GetParError(0),
                                 toyRatioFit->GetParameter(1), toyRatioFit->GetParError(1),
                                 toyRatioFit->GetParameter(2), toyRatioFit->GetParError(2));


            delete toyUHist;
            delete toyVHist;
            delete toyNumHist;
            delete toyDenomHist;
            delete toyRatioGraph;

/////////////////////////////////////////////////////////////////////////////////////

      // now do ratio stuff for histograms with VW randomization
        
      TH1F* toyUHist_VW = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_U_Hist_withVWRand", iter)))->Clone();
      TH1F* toyVHist_VW = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_V_Hist_withVWRand", iter)))->Clone();

      TH1F* toyNumHist_VW = (TH1F*) toyVHist_VW->Clone("numeratorHist");
      toyNumHist_VW->Add(toyUHist_VW, -1);

      TH1F* toyDenomHist_VW = (TH1F*) toyVHist_VW->Clone("denominatorHist");
      toyDenomHist_VW->Add(toyUHist_VW);

            toyRatio_VW_Dir->cd();

            TGraphErrors* toyRatioGraph_VW = ratioFitToyClass.createRatioGraph(toyNumHist_VW, toyDenomHist_VW);
            ratioFitToyClass.fitMethod(toyRatioGraph_VW);

            auto toyRatioFit_VW = toyRatioGraph_VW->GetFunction("threeParamRatioFit");

            if(iter == 0){
              std::cout << "Toy Ratio fit p-value is: " << toyRatioFit_VW->GetProb() << " chi2/ndf: " << toyRatioFit_VW->GetChisquare()/toyRatioFit_VW->GetNDF() << std::endl;
              cout << "Toy Ratio fit R error: " << toyRatioFit_VW->GetParError(1) << " ppm " << endl;
             }

            toyRatioPVal_VW->Fill(toyRatioFit_VW->GetProb());
            toyRatioChi2NDF_VW->Fill(toyRatioFit_VW->GetChisquare()/toyRatioFit_VW->GetNDF());

            toyRatioR_VW->Fill(toyRatioFit_VW->GetParameter(1));

            ratioFitValues_VW->Fill(toyRatioFit_VW->GetParameter(0), toyRatioFit_VW->GetParError(0),
                                    toyRatioFit_VW->GetParameter(1), toyRatioFit_VW->GetParError(1),
                                    toyRatioFit_VW->GetParameter(2), toyRatioFit_VW->GetParError(2));

            delete toyUHist_VW;
            delete toyVHist_VW;
            delete toyNumHist_VW;
            delete toyDenomHist_VW;
            delete toyRatioGraph_VW;

} // end loop over iterations


/////////////////////////////////////////////////////////////////////////////////////

  outputFile->Write();
  delete outputFile;

  return 1;

}
