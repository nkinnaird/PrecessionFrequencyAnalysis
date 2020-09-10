
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

#include "gm2util/blinders/Blinders.hh"

#include "ratioAnalysisDefs.hh"
#include "ratioToyAnalysisConfig.hh"
#include "fiveParamFit.hh"
#include "threeParameterRatioFit.hh"
// #include "residualPlots.hh"
// #include "ratioAnalysisConfig.hh"
#include "plotUtils.hh"


using namespace std;

/////////////////////////////////////////////////////////////////////////////////////

int FitSlowMC(std::string filePath)
{

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  // create output file that will hold plots
  TFile* outputFile = new TFile("slowMCFits.root","RECREATE");
  // outputFile->cd();

  blinding::Blinders::fitType ftype = blinding::Blinders::kOmega_a;
  blinding::Blinders* myBlinder = new blinding::Blinders(ftype); // no blinding for ToyMC

  FiveParamFit fiveParamFitClass(myBlinder);

  ThreeParameterRatioFit threeParRatioFitClass(myBlinder);
  threeParRatioFitClass.setBinWidth(defaultBinWidth);

/////////////////////////////////////////////////////////////////////////////////////

    // setDataset(3, 12); // set dataset case to Endgame (or whatever was used to generate Toy MC hists) - for 'starting' parameters - will print out some other things which won't apply to this Toy MC

/////////////////////////////////////////////////////////////////////////////////////

  int totalFits = 500;

  vector<double> T_Rs, R_Rs;

  vector<double> fitStarts;
  double approxStart = 10000; 
  double realStart = defaultBinWidth * int(approxStart/defaultBinWidth);
  double approxEnd = 300000;
  double realEnd = defaultBinWidth * int(approxEnd/defaultBinWidth);


  for (int fitNum = 0; fitNum < totalFits; ++fitNum)
  {
    fitStarts.push_back(realStart + fitNum * defaultBinWidth);
  }



/////////////////////////////////////////////////////////////////////////////////////


  TH1F* fiveParamHistogram = (TH1F*) inputFile->Get("Toy_5_Param_Hist")->Clone("fiveParamHistogram");

  int nBins = int(approxMaxTime/defaultBinWidth);
  double histMaxTime = nBins*defaultBinWidth;

  TH1F* Uhist = (TH1F*) inputFile->Get("Toy_U_Hist")->Clone("Uhist");
  TH1F* Vhist = (TH1F*) inputFile->Get("Toy_V_Hist")->Clone("Vhist");

  TH1F* toyNumHist = new TH1F("Toy_Num_Hist","Toy_Num_Hist",nBins,0,histMaxTime);
  TH1F* toyDenomHist = new TH1F("Toy_Denom_Hist","Toy_Denom_Hist",nBins,0,histMaxTime);

  toyNumHist->Add(Uhist, Vhist, -1, 1);
  toyDenomHist->Add(Uhist, Vhist);

  auto ratioGraph = threeParRatioFitClass.createRatioGraph(toyNumHist, toyDenomHist);
  ratioGraph->SetName("ratioGraph");



  for (int fitNum = 0; fitNum < totalFits; ++fitNum)
  {

    // if(fitNum == 6 || fitNum == 46 || fitNum == 48 || fitNum == 55 || fitNum == 55 || fitNum == 55 || fitNum == 55 || fitNum == 55){
    //   T_Rs.push_back(0);
    //   R_Rs.push_back(0);
    //   continue;
    // } 

    fiveParamFitClass.setFitRange(fitStarts.at(fitNum), realEnd);
    threeParRatioFitClass.setFitRange(fitStarts.at(fitNum), realEnd);

    fiveParamFitClass.fiveParameterFitMethod(fiveParamHistogram);
    auto fiveParamFittedFunction = fiveParamHistogram->GetFunction("fiveParamFit");
    // cout << endl << "Added 5 Param Fit p-value is: " << fiveParamFittedFunction->GetProb() << " chi2/ndf: " << fiveParamFittedFunction->GetChisquare()/fiveParamFittedFunction->GetNDF() << " R error: " << fiveParamFittedFunction->GetParError(3) << endl;

    T_Rs.push_back(fiveParamFittedFunction->GetParameter(3));

/////////////////////////////////////////////////////////////////////////////////////

    threeParRatioFitClass.fitMethod(ratioGraph);
    auto theRatioFit = ratioGraph->GetFunction("threeParamRatioFit");
    // cout << "Added 3 Param Ratio Fit p-value is: " << theRatioFit->GetProb() << " chi2/ndf: " << theRatioFit->GetChisquare()/theRatioFit->GetNDF() << " R error: " << theRatioFit->GetParError(1) << endl;

    R_Rs.push_back(theRatioFit->GetParameter(1));

    cout << "Fit: " << fitNum << " T p-value: " << fiveParamFittedFunction->GetProb() << " R p-value: " << theRatioFit->GetProb() << " Error on R: " << fiveParamFittedFunction->GetParError(3) << endl;

  } // end of fit loop



  outputFile->cd();


    TGraph* graph_T_Rs = new TGraph();
    TGraph* graph_R_Rs = new TGraph();
    TGraph* graph_diffs = new TGraph();

    graph_R_Rs->SetLineColor(2);


    for (uint pointNo = 0; pointNo < T_Rs.size(); ++pointNo)
    {
      graph_T_Rs->SetPoint(pointNo, fitStarts.at(pointNo), T_Rs.at(pointNo));
      graph_R_Rs->SetPoint(pointNo, fitStarts.at(pointNo), R_Rs.at(pointNo));

      graph_diffs->SetPoint(pointNo, fitStarts.at(pointNo), T_Rs.at(pointNo) - R_Rs.at(pointNo));
    }

    auto myCanv = new TCanvas("name", "title", 200, 200, 1000, 800);
    graph_T_Rs->Draw("AL");
    graph_R_Rs->Draw("LSAME");
    drawg2Lines(myCanv, false);
    drawHorizontalLine(myCanv, 0);

    auto myCanv2 = new TCanvas("name2", "title", 200, 200, 1000, 800);
    graph_diffs->Draw("AL");
    drawg2Lines(myCanv2, false);
    drawHorizontalLine(myCanv2, 0);


    graph_T_Rs->Write("T_Rs");
    graph_R_Rs->Write("R_Rs");
    graph_diffs->Write("diffs");


/////////////////////////////////////////////////////////////////////////////////////

  // outputFile->Write();
  delete outputFile;

  return 1;

}
