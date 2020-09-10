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

using namespace std;

int fitRatioToyHistsDiffRands(std::string filePath)
{
  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen


  // pull in input file
  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  // create output file that will hold plots
  TFile* outputFile = new TFile("toyOutputCompareDiffRands.root","RECREATE");

  // make top directory for output file
  auto topDir = outputFile->mkdir("topDir");

/////////////////////////////////////////////////////////////////////////////////////

  setDataset(1, 12); // set dataset case to 60h - for 'starting' parameters

  double toyFitStart = 30200;
  // double toyFitStart = 5000;
  double toyFitEnd = 650000;

  int totalIters = (*(TVectorD*) inputFile->Get("Iters"))[0]; // total iterations in generated histograms (energy thresholds, etc.)

/////////////////////////////////////////////////////////////////////////////////////

  blinding::Blinders::fitType ftype = blinding::Blinders::kOmega_a;
  blinding::Blinders* myBlinder = new blinding::Blinders(ftype); // no blinding for ToyMC

  FiveParamFit fiveParamFitToyClass(myBlinder);
  fiveParamFitToyClass.setFitRange(toyFitStart, toyFitEnd);

  ThreeParameterRatioFit ratioFitToyClass(myBlinder);
  ratioFitToyClass.setFitRange(toyFitStart, toyFitEnd);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  TF1* truthFunc = (TF1*) inputFile->Get("truthFunc");
  truthFunc->Write();
  // for (int i = 0; i < truthFunc->GetNpar(); ++i) cout << "Truth func parameter " << i << ": " << truthFunc->GetParameter(i) << endl;

/////////////////////////////////////////////////////////////////////////////////////

uint numRands = 6;

for (uint i = 0; i < numRands; ++i)
// for (uint i = 0; i < 1; ++i)
{

  string nameAddOn = "";
  if(i != 0) nameAddOn = "_" + to_string(i);

  auto toyMC_Dir = topDir->mkdir(("ToyMC" + nameAddOn).c_str());
  auto toyHist_Dir = toyMC_Dir->mkdir("Hist");
  auto toyRatio_Dir = toyMC_Dir->mkdir("Ratio");

  toyHist_Dir->cd();

  auto toy5ParamHistPVal = new TH1F(("toy5ParamHistPVal" + nameAddOn).c_str(), "toy5ParamHistPVal; p value; Events", 100, 0, 1);
  auto toy5ParamHistChi2NDF = new TH1F(("toy5ParamHistChi2NDF" + nameAddOn).c_str(), "toy5ParamHistChi2NDF; #chi^{2}/ndf; Events", 50, .9, 1.1);
  auto toyHistFiveR = new TH1F(("toyHistFiveR" + nameAddOn).c_str(), "Toy MC Five Param Fit R Vs Random Seed; R (ppm); Events", 100, -5, 5);

    TNtuple* fiveParameterFitValues = new TNtuple(("fiveParameterFitValues" + nameAddOn).c_str(), "fiveParameterFitValues", "N:N_err:tau:tau_err:A:A_err:R:R_err:phi:phi_err");

  toyRatio_Dir->cd();

  auto toyRatioPVal = new TH1F(("toyRatioPVal" + nameAddOn).c_str(), "toyRatioPVal; p value; Events", 100, 0, 1);
  auto toyRatioChi2NDF = new TH1F(("toyRatioChi2NDF" + nameAddOn).c_str(), "toyRatioChi2NDF; #chi^{2}/ndf; Events", 50, .9, 1.1);
  auto toyRatioR = new TH1F(("toyRatioR" + nameAddOn).c_str(), "Toy Ratio Fit R Vs Random Seed; R (ppm); Events", 100, -5, 5);

    TNtuple* ratioFitValues = new TNtuple(("ratioFitValues" + nameAddOn).c_str(), "ratioFitValues", "A:A_err:R:R_err:phi:phi_err");


  for (int iter = 0; iter < totalIters; ++iter) // loop over iterations / different random seeds
  {
    // cout << "Rand: " << i << " Iter: " << iter << endl;

      TH1F* toyFiveParamHist = (TH1F*) ((TH1F*) inputFile->Get(("topDir/Iter" + to_string(iter) + "/Toy_5_Param_Hist" + nameAddOn).c_str()))->Clone();

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

      TH1F* toyUHist = (TH1F*) ((TH1F*) inputFile->Get(("topDir/Iter" + to_string(iter) + "/Toy_U_Hist" + nameAddOn).c_str()))->Clone();
      TH1F* toyVHist = (TH1F*) ((TH1F*) inputFile->Get(("topDir/Iter" + to_string(iter) + "/Toy_V_Hist" + nameAddOn).c_str()))->Clone();

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

            ratioFitValues->Fill(toyRatioFit->GetParameter(0), toyRatioFit->GetParError(0),
                                 toyRatioFit->GetParameter(1), toyRatioFit->GetParError(1),
                                 toyRatioFit->GetParameter(2), toyRatioFit->GetParError(2));


            toyRatioPVal->Fill(toyRatioFit->GetProb());
            toyRatioChi2NDF->Fill(toyRatioFit->GetChisquare()/toyRatioFit->GetNDF());

            toyRatioR->Fill(toyRatioFit->GetParameter(1));

            delete toyUHist;
            delete toyVHist;
            delete toyNumHist;
            delete toyDenomHist;
            delete toyRatioGraph;

  } // end loop over iterations

    cout << "End rand: " << i << endl;

} // end loop over randomizations


/////////////////////////////////////////////////////////////////////////////////////

  outputFile->Write();
  delete outputFile;

  return 1;

}
