// R__LOAD_LIBRARY(/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_02_01/build_slf6.x86_64/gm2util/lib/libgm2util_blinders.so) // these paths are wrong/need to be updated
// R__ADD_INCLUDE_PATH(/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_02_01/srcs/gm2analyses/macros/RatioMacro/ratioMacroHeaders) // this didn't work for some reason (besides path issue after latest pull)

R__LOAD_LIBRARY(/cvmfs/gm2.opensciencegrid.org/prod/g-2/gm2util/v9_04_00/slf6.x86_64.e15.prof/lib/libgm2util_blinders.so)

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
// #include <TRatioPlot.h> // wait till root v 6_08

#include "gm2util/blinders/Blinders.hh"

#include "../ratioMacroHeaders/ratioAnalysisDefs.hh"
#include "../ratioMacroHeaders/fiveParamFit.hh"
#include "../ratioMacroHeaders/ratioFit.hh"
#include "../ratioMacroHeaders/ratioCBOFit.hh"
#include "../ratioMacroHeaders/residualPlots.hh"

#include "../ratioMacroHeaders/copyFile.hh"

using namespace std;

int toyMC(string filePath, TDirectory* inDir, FitCondStruct toyFitConditions){

  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
// generate histograms for many iterations

  auto toyMCDir = inDir->mkdir(("ToyMC" + toyFitConditions.directoryString).c_str());
  auto toyHistDir = toyMCDir->mkdir("Hist");
  auto toyRatioDir = toyMCDir->mkdir("Ratio");
  auto toyRatioCBODir = toyMCDir->mkdir("RatioCBO");

  toyMCDir->cd();

  auto toy5ParamHistPVal = new TH1F("toy5ParamHistPVal", "toy5ParamHistPVal; p value; Events", 100, 0, 1);
  auto toyRatioPVal = new TH1F("toyRatioPVal", "toyRatioPVal; p value; Events", 100, 0, 1);
  auto toyRatioCBOPVal = new TH1F("toyRatioCBOPVal", "toyRatioCBOPVal; p value; Events", 100, 0, 1);

/////////////////////////////////////////////////////////////////////////////////////

  auto toyHistStatsDir = toyHistDir->mkdir("Stats");
  toyHistStatsDir->cd();

  auto toyHistFiveNPull = new TH1F("toyHistFiveNPull", "toyHistFiveNPull; (N_{fit}-N_{true})/#sigma_{N_{fit}}; Events", 50, -10, 10);
  auto toyHistFiveTauPull = new TH1F("toyHistFiveTauPull", "toyHistFiveTauPull; (#tau_{fit}-#tau_{true})/#sigma_{#tau_{fit}}; Events", 50, -10, 10);
  auto toyHistFiveAPull = new TH1F("toyHistFiveAPull", "toyHistFiveAPull; (A_{fit}-A_{true})/#sigma_{A_{fit}}; Events", 50, -10, 10);
  auto toyHistFiveRPull = new TH1F("toyHistFiveRPull", "toyHistFiveRPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10);
  auto toyHistFivePhasePull = new TH1F("toyHistFivePhasePull", "toyHistFivePhasePull; (#phi_{fit}-#phi_{true})/#sigma_{#phi_{fit}}; Events", 50, -10, 10);

  toyHistDir->cd();
  auto toyHistFivePhase = new TH1F("toyHistFivePhase", "toyHistFivePhase; #phi_{fit}; Events", 100, TMath::Pi()-.1, TMath::Pi()+.1);

/////////////////////////////////////////////////////////////////////////////////////

  auto toyRatioStatsDir = toyRatioDir->mkdir("Stats");
  toyRatioStatsDir->cd();

  auto toyRatioAPull = new TH1F("toyRatioAPull", "toyRatioAPull; (A_{fit}-A_{true})/#sigma_{A_{fit}}; Events", 50, -10, 10);
  auto toyRatioRPull = new TH1F("toyRatioRPull", "toyRatioRPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10);
  auto toyRatioPhasePull = new TH1F("toyRatioPhasePull", "toyRatioPhasePull; (#phi_{fit}-#phi_{true})/#sigma_{#phi_{fit}}; Events", 50, -20, 10); 

    auto toyRatioADiv = new TH1F("toyRatioADiv", "toyRatioADiv; (A_{fit}-A_{true})/A_{true} (ppm); Events", 50, -100, 100);
    auto toyRatioPhaseDiv = new TH1F("toyRatioPhaseDiv", "toyRatioPhaseDiv; (#phi_{fit}-#phi_{true})/#phi_{true} (ppm); Events", 50, -100, 100);

    auto toyRatioRDiff = new TH1F("toyRatioRDiff", "toyRatioRDiff; (R_{fit}-R_{true}) (ppm); Events", 50, -25, 25);

    auto toyRatioChi2NDF = new TH1F("toyRatioChi2NDF", "toyRatioChi2NDF; #chi^{2}/ndf; Events", 50, .9, 1.1);

  toyRatioDir->cd();
  auto toyRatioPhase = new TH1F("toyRatioPhase", "toyRatioPhase; #phi_{fit}; Events", 100, TMath::Pi()-.1, TMath::Pi()+.1);

/////////////////////////////////////////////////////////////////////////////////////

  auto toyRatioCBOStatsDir = toyRatioCBODir->mkdir("Stats");
  toyRatioCBOStatsDir->cd();

  auto tRCBOAPull = new TH1F("tRCBOAPull", "tRCBOAPull; (A_{fit}-A_{true})/#sigma_{A_{fit}}; Events", 50, -20, 10);
  auto tRCBORPull = new TH1F("tRCBORPull", "tRCBORPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10);
  auto tRCBOPhasePull = new TH1F("tRCBOPhasePull", "tRCBOPhasePull; (#phi_{fit}-#phi_{true})/#sigma_{#phi_{fit}}; Events", 50, -20, 10); 

  auto toyRCBO_cboFreq_Pull = new TH1F("toyRCBO_cboFreq_Pull", "toyRCBO_cboFreq_Pull; (#omega_{cbo-fit}-#omega_{cbo-true})/#sigma_{#omega_{cbo-fit}}; Events", 50, -10, 10);
  auto toyRCBO_cboTau_Pull = new TH1F("toyRCBO_cboTau_Pull", "toyRCBO_cboTau_Pull; (#tau_{cbo-fit}-#tau_{cbo-true})/#sigma_{#tau_{cbo-fit}}; Events", 50, -10, 10); 

  auto tRCBONampPull = new TH1F("tRCBONampPull", "tRCBONampPull; (A_{N_{cbo-fit}}-A_{N_{cbo-true}})/#sigma_{A_{N_{cbo-fit}}}; Events", 50, -10, 10);
  auto tRCBONphasePull = new TH1F("tRCBONphasePull", "tRCBONphasePull; (#phi_{N_{cbo-fit}}-#phi_{N_{cbo-true}})/#sigma_{#phi_{N_{cbo-fit}}}; Events", 50, -10, 10); 

  auto tRCBOAampPull = new TH1F("tRCBOAampPull", "tRCBOAampPull; (A_{A_{cbo-fit}}-A_{A_{cbo-true}})/#sigma_{A_{A_{cbo-fit}}}; Events", 50, -10, 10);
  auto tRCBOAphasePull = new TH1F("tRCBOAphasePull", "tRCBOAphasePull; (#phi_{A_{cbo-fit}}-#phi_{A_{cbo-true}})/#sigma_{#phi_{A_{cbo-fit}}}; Events", 50, -10, 10); 

  auto tRCBOPampPull = new TH1F("tRCBOPampPull", "tRCBOPampPull; (A_{#phi_{cbo-fit}}-A_{#phi_{cbo-true}})/#sigma_{A_{#phi_{cbo-fit}}}; Events", 50, -10, 10);
  auto tRCBOPphasePull = new TH1F("tRCBOPphasePull", "tRCBOPphasePull; (#phi_{#phi_{cbo-fit}}-#phi_{#phi_{cbo-true}})/#sigma_{#phi_{#phi_{cbo-fit}}}; Events", 50, -10, 10); 

    auto tRCBOADiv = new TH1F("tRCBOADiv", "tRCBOADiv; (A_{fit}-A_{true})/A_{true} (ppm); Events", 50, -1000, 100);
    auto tRCBOPhaseDiv = new TH1F("tRCBOPhaseDiv", "tRCBOPhaseDiv; (#phi_{fit}-#phi_{true})/#phi_{true} (ppm); Events", 50, -1000, 100);

    auto tRCBORDiff = new TH1F("tRCBORDiff", "tRCBORDiff; (R_{fit}-R_{true}) (ppm); Events", 50, -25, 25);

    auto tRCBOChi2NDF = new TH1F("tRCBOChi2NDF", "tRCBOChi2NDF; #chi^{2}/ndf; Events", 50, .9, 1.1);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  std::vector<string> fileVector;

	std::size_t found = filePath.find_last_of(".");
	std::string extension = filePath.substr(found);

if(extension == ".root") fileVector.push_back(filePath);
else if(extension == ".txt"){
  std::ifstream inList(filePath);
  string path;
  while(inList >> path) fileVector.push_back(path);
}
else{
	printf("Not passing files in correctly.\n");
    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
// do some stuff here to initialize single objects used in all iterations, with some knowledge of mc conditions

  double pts;
  double toyFitEnd;

  // pull in first input file
   TFile *firstFile  = TFile::Open(fileVector.front().c_str());
   if (firstFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  TH1F* initialHist = (TH1F*) firstFile->Get(("topDir/" + toyFitConditions.directoryString + "/Toy_5_Param_Hist").c_str());

  pts = initialHist->GetEntries();
  toyFitEnd = 200000 + log10(int(pts/1e7))*100000; // function to change fit end time based on number of events, 200 us + some timefor events greater than 1e7 events - should probably clean this up
  // toyFitEnd = 300000;

  delete initialHist;

  firstFile->Close();

/////////////////////////////////////////////////////////////////////////////////////

  blinding::Blinders::fitType ftype = blinding::Blinders::kOmega_a;
  blinding::Blinders* myBlinder = new blinding::Blinders(ftype); // no blinding for ToyMC
  int numberOfCBOFitPars = 7; // 7, 9, 11
  double cbo_freq = 0;

  FiveParamFit fiveParamFitToyClass(toyFitConditions.fitStartTime, toyFitEnd, myBlinder);
  RatioFit ratioFitToyClass(toyFitConditions.fitStartTime, toyFitEnd, myBlinder);
  RatioCBOFit ratioCBOFitClass(toyFitConditions.fitStartTime, toyFitEnd, myBlinder, numberOfCBOFitPars);

  ResidualPlots residualPlotsToyClass;
    double wholeRange_min = 5000; // use these numbers to cut out the splash and spike at the beginning and end of the data (copied from ratioMacro, can be changed for this toy)
    double wholeRange_max = 650000;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


for (int fileNum = 0; fileNum < int(fileVector.size()); ++fileNum)
{

  TFile* inputFile = TFile::Open(fileVector.at(fileNum).c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

   cout << "Filenum: " << fileNum << endl;

   TVectorD* histSavedParameters = (TVectorD*) inputFile->Get(("topDir/" + toyFitConditions.directoryString + "/SavedParameters/parameterStore").c_str());

      TF1* truthFunc = (TF1*) inputFile->Get("truthFunc"); // could maybe put this outside the file loop, but could cause problems if truth func changes even though it shouldn't

      TH1F* toyFiveParamHist = (TH1F*) ((TH1F*) inputFile->Get(("topDir/" + toyFitConditions.directoryString + "/Toy_5_Param_Hist").c_str()))->Clone();

            fiveParamFitToyClass.fiveParameterFitMethod(toyFiveParamHist);

            auto toyHistFiveFit = toyFiveParamHist->GetFunction("fiveParamFit");

            double histPValue = toyHistFiveFit->GetProb();

            toy5ParamHistPVal->Fill(histPValue);

            toyHistFiveNPull    ->Fill( (toyHistFiveFit->GetParameter(0) - truthFunc->GetParameter(0)) / toyHistFiveFit->GetParError(0) );
            toyHistFiveTauPull  ->Fill( (toyHistFiveFit->GetParameter(1) - truthFunc->GetParameter(1)) / toyHistFiveFit->GetParError(1) );
            toyHistFiveAPull    ->Fill( (toyHistFiveFit->GetParameter(2) - truthFunc->GetParameter(2)) / toyHistFiveFit->GetParError(2) );
            toyHistFiveRPull    ->Fill( (toyHistFiveFit->GetParameter(3) - truthFunc->GetParameter(3)) / toyHistFiveFit->GetParError(3) );
            toyHistFivePhasePull->Fill( (toyHistFiveFit->GetParameter(4) - truthFunc->GetParameter(4)) / toyHistFiveFit->GetParError(4) );

            toyHistFivePhase->Fill(toyHistFiveFit->GetParameter(4));

            if(fileNum == 0){
              std::cout << "Toy hist 5 param fit p-value is: " << histPValue << " chi2/ndf: " << toyHistFiveFit->GetChisquare()/toyHistFiveFit->GetNDF() << std::endl;
              cout << "Toy hist 5 param fit R error: " << toyHistFiveFit->GetParError(3) << " ppm "  << endl;

              auto histIterNumDir = toyHistDir->mkdir(toyFitConditions.directoryString.c_str());
              histIterNumDir->cd();

              truthFunc->Write();
              toyFiveParamHist->Write();

              TH1F* fiveFitFFT = residualPlotsToyClass.makeResidualPlots("ToyHist", toyFiveParamHist, "residual_fitRange", toyFitConditions.fitStartTime, toyFitEnd);
              cbo_freq = residualPlotsToyClass.findPeaks(fiveFitFFT);
              residualPlotsToyClass.makeResidualPlots("ToyHist", toyFiveParamHist, "residual_wholeRange", wholeRange_min, wholeRange_max);
            }

            delete toyFiveParamHist;

    /////////////////////////////////////////////////////////////////////////////////////

      TH1F* toyUHist = (TH1F*) ((TH1F*) inputFile->Get(("topDir/" + toyFitConditions.directoryString + "/Toy_U_Hist").c_str()))->Clone();
      TH1F* toyVHist = (TH1F*) ((TH1F*) inputFile->Get(("topDir/" + toyFitConditions.directoryString + "/Toy_V_Hist").c_str()))->Clone();
      TH1F* toyNumHist = (TH1F*) ((TH1F*) inputFile->Get(("topDir/" + toyFitConditions.directoryString + "/Toy_Num_Hist").c_str()))->Clone();
      TH1F* toyDenomHist = (TH1F*) ((TH1F*) inputFile->Get(("topDir/" + toyFitConditions.directoryString + "/Toy_Denom_Hist").c_str()))->Clone();

            toyRatioDir->cd();

            TGraphErrors* toyRatioGraph = ratioFitToyClass.createRatioGraph(toyNumHist, toyDenomHist);
            ratioFitToyClass.ratioFitMethod(toyRatioGraph);

            auto toyRatioFit = toyRatioGraph->GetFunction("threeParamRatioFit");

            if(fileNum == 0)
            {
              std::cout << "Toy Ratio fit p-value is: " << toyRatioFit->GetProb() << " chi2/ndf: " << toyRatioFit->GetChisquare()/toyRatioFit->GetNDF() << std::endl;
              cout << "Toy Ratio fit R error: " << toyRatioFit->GetParError(1) << " ppm " << endl;

              auto ratioIterNumDir = toyRatioDir->mkdir(toyFitConditions.directoryString.c_str());
              ratioIterNumDir->cd();

              toyRatioGraph->SetTitle("Toy Ratio Graph; Time (ns); R (unitless)");
              toyRatioGraph->Write("Toy_Ratio_Graph");

              toyUHist->Write();
              toyVHist->Write();
              toyNumHist->Write();
              toyDenomHist->Write();

              residualPlotsToyClass.makeResidualPlots("Toy_Ratio", toyRatioGraph, "residual_fitRange", toyFitConditions.fitStartTime, toyFitEnd);
              residualPlotsToyClass.makeResidualPlots("Toy_Ratio", toyRatioGraph, "residual_wholeRange", wholeRange_min+g2Period, wholeRange_max);
             }

            toyRatioPVal->Fill(toyRatioFit->GetProb());

            toyRatioAPull    ->Fill( (toyRatioFit->GetParameter(0) - truthFunc->GetParameter(2)) / toyRatioFit->GetParError(0) );
            toyRatioRPull    ->Fill( (toyRatioFit->GetParameter(1) - truthFunc->GetParameter(3)) / toyRatioFit->GetParError(1) );
            toyRatioPhasePull->Fill( (toyRatioFit->GetParameter(2) - truthFunc->GetParameter(4)) / toyRatioFit->GetParError(2) );   

            toyRatioADiv->Fill( 1e6 * (toyRatioFit->GetParameter(0) - truthFunc->GetParameter(2)) / truthFunc->GetParameter(2) );
            toyRatioPhaseDiv->Fill( 1e6 * (toyRatioFit->GetParameter(2) - truthFunc->GetParameter(4)) / truthFunc->GetParameter(4) );

            toyRatioRDiff->Fill(toyRatioFit->GetParameter(1) - truthFunc->GetParameter(3));

            toyRatioChi2NDF->Fill(toyRatioFit->GetChisquare()/toyRatioFit->GetNDF());

            toyRatioPhase->Fill(toyRatioFit->GetParameter(2));


    /////////////////////////////////////////////////////////////////////////////////////
            //cbo testing area

            toyRatioCBODir->cd();

            ratioCBOFitClass.setHalfPeriod((*histSavedParameters)[0]/2.); // set the half period for the cbo fitting

            TGraphErrors* toyRatioCBOGraph = ratioFitToyClass.createRatioGraph(toyNumHist, toyDenomHist);
            ratioCBOFitClass.setCBOFreq(cbo_freq);
            ratioCBOFitClass.setThreeParamRatioFunction(toyRatioFit);
            ratioCBOFitClass.ratioCBOFitMethod(toyRatioCBOGraph);

            auto ratioCBOToyFit = toyRatioCBOGraph->GetFunction("cboRatioFitFunc");

            if(fileNum == 0)
            {
              cout << "ratioCBOToyFit p value: " << ratioCBOToyFit->GetProb() << " chi2/ndf: " << ratioCBOToyFit->GetChisquare()/ratioCBOToyFit->GetNDF() << endl;

              auto ratioCBOIterNumDir = toyRatioCBODir->mkdir(toyFitConditions.directoryString.c_str());
              ratioCBOIterNumDir->cd();

              toyRatioCBOGraph->SetTitle("Toy Ratio CBO Graph; Time (ns); R (unitless)");
              toyRatioCBOGraph->Write("Toy_Ratio_CBO_Graph");

              toyUHist->Write();
              toyVHist->Write();
              toyNumHist->Write();
              toyDenomHist->Write();

              residualPlotsToyClass.makeResidualPlots("Toy_Ratio_CBO", toyRatioCBOGraph, "residual_fitRange", toyFitConditions.fitStartTime, toyFitEnd);
              residualPlotsToyClass.makeResidualPlots("Toy_Ratio_CBO", toyRatioCBOGraph, "residual_wholeRange", wholeRange_min+g2Period, wholeRange_max);
            }


            toyRatioCBOPVal->Fill(ratioCBOToyFit->GetProb());

            tRCBOAPull    ->Fill( (ratioCBOToyFit->GetParameter(0) - truthFunc->GetParameter(2)) / ratioCBOToyFit->GetParError(0) );
            tRCBORPull    ->Fill( (ratioCBOToyFit->GetParameter(1) - truthFunc->GetParameter(3)) / ratioCBOToyFit->GetParError(1) );
            tRCBOPhasePull->Fill( (ratioCBOToyFit->GetParameter(2) - truthFunc->GetParameter(4)) / ratioCBOToyFit->GetParError(2) ); 

            toyRCBO_cboFreq_Pull->Fill( (ratioCBOToyFit->GetParameter(3) - truthFunc->GetParameter(5)) / ratioCBOToyFit->GetParError(3));
            toyRCBO_cboTau_Pull ->Fill( (ratioCBOToyFit->GetParameter(4) - truthFunc->GetParameter(6)) / ratioCBOToyFit->GetParError(4)); 

            tRCBONampPull  ->Fill( (ratioCBOToyFit->GetParameter(5) - truthFunc->GetParameter(7)) / ratioCBOToyFit->GetParError(5) );
            tRCBONphasePull->Fill( (ratioCBOToyFit->GetParameter(6) - truthFunc->GetParameter(8)) / ratioCBOToyFit->GetParError(6) );
            tRCBOAampPull  ->Fill( (ratioCBOToyFit->GetParameter(7) - truthFunc->GetParameter(9)) / ratioCBOToyFit->GetParError(7) );
            tRCBOAphasePull->Fill( (ratioCBOToyFit->GetParameter(8) - truthFunc->GetParameter(10)) / ratioCBOToyFit->GetParError(8) );
            tRCBOPampPull  ->Fill( (ratioCBOToyFit->GetParameter(9) - truthFunc->GetParameter(11)) / ratioCBOToyFit->GetParError(9) );
            tRCBOPphasePull->Fill( (ratioCBOToyFit->GetParameter(10) - truthFunc->GetParameter(12)) / ratioCBOToyFit->GetParError(10) );

            tRCBOADiv->Fill( 1e6 * (ratioCBOToyFit->GetParameter(0) - truthFunc->GetParameter(2)) / truthFunc->GetParameter(2) );
            tRCBOPhaseDiv->Fill( 1e6 * (ratioCBOToyFit->GetParameter(2) - truthFunc->GetParameter(4)) / truthFunc->GetParameter(4) );

            tRCBORDiff->Fill(ratioCBOToyFit->GetParameter(1) - truthFunc->GetParameter(3));

            tRCBOChi2NDF->Fill(ratioCBOToyFit->GetChisquare()/ratioCBOToyFit->GetNDF());


    /////////////////////////////////////////////////////////////////////////////////////


            delete toyUHist;
            delete toyVHist;
            delete toyNumHist;
            delete toyDenomHist;

            delete truthFunc;


    inputFile->Close();


} //end file path loop

  return 1;

}


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

int runRatioToy(std::string filePath)
{
    // create output file that will hold plots
  TFile* outputFile = new TFile("toyOutputNew.root","RECREATE");

  // make top directory for output file
  auto topDir = outputFile->mkdir("topDir");

/////////////////////////////////////////////////////////////////////////////////////

  double toyFitStart = g2Period;
  // double toyFitStart = 30000;

  std::size_t found = filePath.find_last_of(".");
  std::string extension = filePath.substr(found);

  string dirString;

  if(extension == ".root") // if I'm running over a single root file, which possibly has multiple iterations
  {
    CopyFile(filePath.c_str()); // copy the input file into the output file (only if it's a single input file)

    found = filePath.find_last_of("/"); // do this to get the filename which corresponds to the directory name of the copied file
    std::string fileName;
    if(found < filePath.size()) fileName = filePath.substr(found);
    else fileName = filePath;

    TVectorD* itersDouble = ((TVectorD*) outputFile->Get((filePath + "/Iters").c_str()));
    int totalIters = (*itersDouble)[0];

/////////////////////////////////////////////////////////////////////////////////////

    for (int iter = 0; iter < totalIters; ++iter)
    {
      dirString = "Iter" + to_string(iter);
      FitCondStruct toyFitConditions(dirString, toyFitStart, 0, 0); // fit end range isn't set here because it's automatically set later based on the number of events in the toy, pileup scale factor is dealt with in makeToyHistFromTH2Pileup for now as well
      toyMC(filePath, topDir, toyFitConditions);
    }

  }
  else // if I'm running over a list of root files, where each file has a single iteration
  {
    dirString = "Iter0";
    FitCondStruct toyFitConditions(dirString, toyFitStart, 0, 0);
    toyMC(filePath, topDir, toyFitConditions);
  }

/////////////////////////////////////////////////////////////////////////////////////

  outputFile->Write();
  delete outputFile;

  return 1;

}
