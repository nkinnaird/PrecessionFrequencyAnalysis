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

#include "gm2util/blinders/Blinders.hh"

#include "ratioAnalysisDefs.hh"
#include "ratioAnalysisConfig.hh"
#include "fiveParamFit.hh"
#include "threeParameterRatioFit.hh"
#include "fullRatioFit.hh"
#include "residualPlots.hh"
#include "copyFile.hh"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////

int runRatioToyCompare(std::string filePath)
{
  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen


  // pull in input file
  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  // create output file that will hold plots
  TFile* outputFile = new TFile("toyOutputCompare.root","RECREATE");

  // make top directory for output file
  auto topDir = outputFile->mkdir("topDir");

/////////////////////////////////////////////////////////////////////////////////////

  setDataset(1, 12); // set dataset case to 60h - for 'starting' parameters

  double toyFitStart = 30000;
  double toyFitEnd = 650000;

  int totalIters = (*(TVectorD*) inputFile->Get("Iters"))[0]; // total iterations in generated histograms (energy thresholds, etc.)

/////////////////////////////////////////////////////////////////////////////////////
 
  auto toyMCDir = topDir->mkdir("ToyMC");
  auto toyHistDir = toyMCDir->mkdir("Hist");
  auto toyRatioDir = toyMCDir->mkdir("Ratio");

  auto comparisonDir = toyMCDir->mkdir("Comparison");

/////////////////////////////////////////////////////////////////////////////////////

  toyMCDir->cd();

  auto toy5ParamHistPVal = new TH1F("toy5ParamHistPVal", "toy5ParamHistPVal; p value; Events", 100, 0, 1);
  auto toyRatioPVal = new TH1F("toyRatioPVal", "toyRatioPVal; p value; Events", 100, 0, 1);

  auto toy5ParamHistChi2NDF = new TH1F("toy5ParamHistChi2NDF", "toy5ParamHistChi2NDF; #chi^{2}/ndf; Events", 50, .9, 1.1);
  auto toyRatioChi2NDF = new TH1F("toyRatioChi2NDF", "toyRatioChi2NDF; #chi^{2}/ndf; Events", 50, .9, 1.1);

/////////////////////////////////////////////////////////////////////////////////////

  auto toyHistStatsDir = toyHistDir->mkdir("Stats");
  toyHistStatsDir->cd();

  auto toyHistFiveNPull = new TH1F("toyHistFiveNPull", "toyHistFiveNPull; (N_{fit}-N_{true})/#sigma_{N_{fit}}; Events", 50, -10, 10);
  auto toyHistFiveTauPull = new TH1F("toyHistFiveTauPull", "toyHistFiveTauPull; (#tau_{fit}-#tau_{true})/#sigma_{#tau_{fit}}; Events", 50, -10, 10);
  auto toyHistFiveAPull = new TH1F("toyHistFiveAPull", "toyHistFiveAPull; (A_{fit}-A_{true})/#sigma_{A_{fit}}; Events", 50, -10, 10);
  auto toyHistFiveRPull = new TH1F("toyHistFiveRPull", "toyHistFiveRPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10);
  auto toyHistFivePhasePull = new TH1F("toyHistFivePhasePull", "toyHistFivePhasePull; (#phi_{fit}-#phi_{true})/#sigma_{#phi_{fit}}; Events", 50, -10, 10);

  auto toyHistFiveR = new TH1F("toyHistFiveR", "Toy MC Five Param Fit R Vs Random Seed; R (ppm); Events", 100, -5, 5);

/////////////////////////////////////////////////////////////////////////////////////

  auto toyRatioStatsDir = toyRatioDir->mkdir("Stats");
  toyRatioStatsDir->cd();

  auto toyRatioAPull = new TH1F("toyRatioAPull", "toyRatioAPull; (A_{fit}-A_{true})/#sigma_{A_{fit}}; Events", 50, -10, 10);
  auto toyRatioRPull = new TH1F("toyRatioRPull", "toyRatioRPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10);
  auto toyRatioPhasePull = new TH1F("toyRatioPhasePull", "toyRatioPhasePull; (#phi_{fit}-#phi_{true})/#sigma_{#phi_{fit}}; Events", 50, -20, 10); 

  auto toyRatioRDiff = new TH1F("toyRatioRDiff", "toyRatioRDiff; (R_{fit}-R_{true}) (ppm); Events", 50, -25, 25);

  auto toyRatioR = new TH1F("toyRatioR", "Toy Ratio Fit R Vs Random Seed; R (ppm); Events", 100, -5, 5);

/////////////////////////////////////////////////////////////////////////////////////

  comparisonDir->cd();

  // auto RDiff = new TH1F("RDiff", "RDiff; (R_{ratio}-R_{5 param}) (ppm); Events", 50, -25, 25);



/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

  blinding::Blinders::fitType ftype = blinding::Blinders::kOmega_a;
  blinding::Blinders* myBlinder = new blinding::Blinders(ftype); // no blinding for ToyMC

  FiveParamFit fiveParamFitToyClass(myBlinder);
  fiveParamFitToyClass.setFitRange(toyFitStart, toyFitEnd);

  ThreeParameterRatioFit ratioFitToyClass(myBlinder);
  ratioFitToyClass.setFitRange(toyFitStart, toyFitEnd);

  ResidualPlots residualPlotsToyClass;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  TF1* truthFunc = (TF1*) inputFile->Get("truthFunc"); // could maybe put this outside the file loop, but could cause problems if truth func changes even though it shouldn't

  for (int i = 0; i < truthFunc->GetNpar(); ++i) cout << "Truth func parameter " << i << ": " << truthFunc->GetParameter(i) << endl;

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

            toyHistFiveNPull    ->Fill( (toyHistFiveFit->GetParameter(0) - truthFunc->GetParameter(0)) / toyHistFiveFit->GetParError(0) );
            toyHistFiveTauPull  ->Fill( (toyHistFiveFit->GetParameter(1) - truthFunc->GetParameter(1)) / toyHistFiveFit->GetParError(1) );
            toyHistFiveAPull    ->Fill( (toyHistFiveFit->GetParameter(2) - truthFunc->GetParameter(2)) / toyHistFiveFit->GetParError(2) );
            toyHistFiveRPull    ->Fill( (toyHistFiveFit->GetParameter(3) - truthFunc->GetParameter(3)) / toyHistFiveFit->GetParError(3) );
            toyHistFivePhasePull->Fill( (toyHistFiveFit->GetParameter(4) - truthFunc->GetParameter(4)) / toyHistFiveFit->GetParError(4) );

            toyHistFiveR->Fill(toyHistFiveFit->GetParameter(3));

            if(iter == 0){
              std::cout << "Toy hist 5 param fit p-value is: " << histPValue << " chi2/ndf: " << toyHistFiveFit->GetChisquare()/toyHistFiveFit->GetNDF() << std::endl;
              cout << "Toy hist 5 param fit R error: " << toyHistFiveFit->GetParError(3) << " ppm "  << endl;

              auto histIterNumDir = toyHistDir->mkdir("Iter0Results");
              histIterNumDir->cd();

              truthFunc->Write();
              toyFiveParamHist->Write();

              // residualPlotsToyClass.makeResidualPlots("ToyHist", toyFiveParamHist, "residual_fitRange", toyFitStart, toyFitEnd, true);
            }

            delete toyFiveParamHist;

    /////////////////////////////////////////////////////////////////////////////////////

      TH1F* toyUHist = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_U_Hist", iter)))->Clone();
      TH1F* toyVHist = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_V_Hist", iter)))->Clone();
      // TH1F* toyNumHist = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_Num_Hist", iter)))->Clone();
      // TH1F* toyDenomHist = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_Denom_Hist", iter)))->Clone();

      TH1F* toyNumHist = (TH1F*) toyVHist->Clone("numeratorHist");
      toyNumHist->Add(toyUHist, -1);

      TH1F* toyDenomHist = (TH1F*) toyVHist->Clone("denominatorHist");
      toyDenomHist->Add(toyUHist);

            toyRatioDir->cd();

            TGraphErrors* toyRatioGraph = ratioFitToyClass.createRatioGraph(toyNumHist, toyDenomHist);
            ratioFitToyClass.fitMethod(toyRatioGraph);

            auto toyRatioFit = toyRatioGraph->GetFunction("threeParamRatioFit");

            if(iter == 0)
            {
              std::cout << "Toy Ratio fit p-value is: " << toyRatioFit->GetProb() << " chi2/ndf: " << toyRatioFit->GetChisquare()/toyRatioFit->GetNDF() << std::endl;
              cout << "Toy Ratio fit R error: " << toyRatioFit->GetParError(1) << " ppm " << endl;

              auto ratioIterNumDir = toyRatioDir->mkdir("Iter0Results");
              ratioIterNumDir->cd();

              toyRatioGraph->SetTitle("Toy Ratio Graph; Time (ns); R (unitless)");
              toyRatioGraph->Write("Toy_Ratio_Graph");

              toyUHist->Write();
              toyVHist->Write();
              toyNumHist->Write();
              toyDenomHist->Write();

              // residualPlotsToyClass.makeResidualPlots("Toy_Ratio", toyRatioGraph, "residual_fitRange", toyFitStart, toyFitEnd, true);
             }

            toyRatioPVal->Fill(toyRatioFit->GetProb());
            toyRatioChi2NDF->Fill(toyRatioFit->GetChisquare()/toyRatioFit->GetNDF());

            toyRatioAPull    ->Fill( (toyRatioFit->GetParameter(0) - truthFunc->GetParameter(2)) / toyRatioFit->GetParError(0) );
            toyRatioRPull    ->Fill( (toyRatioFit->GetParameter(1) - truthFunc->GetParameter(3)) / toyRatioFit->GetParError(1) );
            toyRatioPhasePull->Fill( (toyRatioFit->GetParameter(2) - truthFunc->GetParameter(4)) / toyRatioFit->GetParError(2) );   

            toyRatioRDiff->Fill(toyRatioFit->GetParameter(1) - truthFunc->GetParameter(3));

            toyRatioR->Fill(toyRatioFit->GetParameter(1));

            delete toyUHist;
            delete toyVHist;
            delete toyNumHist;
            delete toyDenomHist;

} // end loop over iterations

/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

  outputFile->Write();
  delete outputFile;

  return 1;

}
