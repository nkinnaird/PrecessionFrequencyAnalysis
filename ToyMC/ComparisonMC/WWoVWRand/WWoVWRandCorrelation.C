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
#include <TPaveStats.h>

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;

double calculateAverage(vector<double> inVec)
{
  double avg = 0;
  for (uint i = 0; i < inVec.size(); ++i) avg += inVec.at(i);
  avg /= inVec.size();
  return avg;
}


int WWoVWRandCorrelation(std::string fileList)
{
  // gStyle->SetOptStat(111111);
  // gStyle->SetOptTitle(1);
  // gStyle->SetTitleStyle(0);
  // gStyle->SetPadRightMargin(.05);
  // gStyle->SetPadLeftMargin(.15);

  std::vector<string> fileVector;

  std::ifstream inList(fileList);
  string path;
  while(inList >> path) fileVector.push_back(path);

/////////////////////////////////////////////////////////////////////////////////////

  uint totalIters = 100; // this I should pull from the file


  double fiveParRs[totalIters][fileVector.size()]; // put iteration loop first so that I could access the arrays more easily via pointers
  double fiveParRerrs[totalIters][fileVector.size()];

    double fiveParRs_Avg[fileVector.size()];

  double fiveParRs_VW[totalIters][fileVector.size()];
  double fiveParRerrs_VW[totalIters][fileVector.size()];

    double fiveParRs_VW_Avg[fileVector.size()];

  double ratioRs[totalIters][fileVector.size()];
  double ratioRerrs[totalIters][fileVector.size()];

    double ratioRs_Avg[fileVector.size()];

  double ratioRs_VW[totalIters][fileVector.size()];
  double ratioRerrs_VW[totalIters][fileVector.size()];

    double ratioRs_VW_Avg[fileVector.size()];

/////////////////////////////////////////////////////////////////////////////////////


  for (uint fileNum = 0; fileNum < fileVector.size(); ++fileNum)
  {
    TFile* inputFile = TFile::Open(fileVector.at(fileNum).c_str());
     if (inputFile == 0) {
        printf("Error: cannot open file\n");
        return 0;
     }

    cout << "Filenum: " << fileNum << " " << fileVector.at(fileNum) << endl;

/////////////////////////////////////////////////////////////////////////////////////


      TNtuple *fiveParamFitValues = (TNtuple*)inputFile->Get("topDir/ToyMC/Hist/fiveParameterFitValues");
      
      float N_fivePar, N_err_fivePar, tau_fivePar, tau_err_fivePar, A_fivePar, A_err_fivePar, R_fivePar, R_err_fivePar, phi_fivePar, phi_err_fivePar;
      fiveParamFitValues->SetBranchAddress("N", &N_fivePar);
      fiveParamFitValues->SetBranchAddress("N_err", &N_err_fivePar);
      fiveParamFitValues->SetBranchAddress("tau", &tau_fivePar);
      fiveParamFitValues->SetBranchAddress("tau_err", &tau_err_fivePar);
      fiveParamFitValues->SetBranchAddress("A", &A_fivePar);
      fiveParamFitValues->SetBranchAddress("A_err", &A_err_fivePar);
      fiveParamFitValues->SetBranchAddress("R", &R_fivePar);
      fiveParamFitValues->SetBranchAddress("R_err", &R_err_fivePar);
      fiveParamFitValues->SetBranchAddress("phi", &phi_fivePar);
      fiveParamFitValues->SetBranchAddress("phi_err", &phi_err_fivePar);

      TNtuple *fiveParamFitValues_VW = (TNtuple*)inputFile->Get("topDir/ToyMC_VW/Hist/fiveParameterFitValues_VW");
      
      float N_fivePar_VW, N_err_fivePar_VW, tau_fivePar_VW, tau_err_fivePar_VW, A_fivePar_VW, A_err_fivePar_VW, R_fivePar_VW, R_err_fivePar_VW, phi_fivePar_VW, phi_err_fivePar_VW;
      fiveParamFitValues_VW->SetBranchAddress("N", &N_fivePar_VW);
      fiveParamFitValues_VW->SetBranchAddress("N_err", &N_err_fivePar_VW);
      fiveParamFitValues_VW->SetBranchAddress("tau", &tau_fivePar_VW);
      fiveParamFitValues_VW->SetBranchAddress("tau_err", &tau_err_fivePar_VW);
      fiveParamFitValues_VW->SetBranchAddress("A", &A_fivePar_VW);
      fiveParamFitValues_VW->SetBranchAddress("A_err", &A_err_fivePar_VW);
      fiveParamFitValues_VW->SetBranchAddress("R", &R_fivePar_VW);
      fiveParamFitValues_VW->SetBranchAddress("R_err", &R_err_fivePar_VW);
      fiveParamFitValues_VW->SetBranchAddress("phi", &phi_fivePar_VW);
      fiveParamFitValues_VW->SetBranchAddress("phi_err", &phi_err_fivePar_VW);


      TNtuple *ratioFitValues = (TNtuple*)inputFile->Get("topDir/ToyMC/Ratio/ratioFitValues");
      
      float A_ratio, A_err_ratio, R_ratio, R_err_ratio, phi_ratio, phi_err_ratio;
      ratioFitValues->SetBranchAddress("A", &A_ratio);
      ratioFitValues->SetBranchAddress("A_err", &A_err_ratio);
      ratioFitValues->SetBranchAddress("R", &R_ratio);
      ratioFitValues->SetBranchAddress("R_err", &R_err_ratio);
      ratioFitValues->SetBranchAddress("phi", &phi_ratio);
      ratioFitValues->SetBranchAddress("phi_err", &phi_err_ratio);


      TNtuple *ratioFitValues_VW = (TNtuple*)inputFile->Get("topDir/ToyMC_VW/Ratio/ratioFitValues_VW");
      
      float A_ratio_VW, A_err_ratio_VW, R_ratio_VW, R_err_ratio_VW, phi_ratio_VW, phi_err_ratio_VW;
      ratioFitValues_VW->SetBranchAddress("A", &A_ratio_VW);
      ratioFitValues_VW->SetBranchAddress("A_err", &A_err_ratio_VW);
      ratioFitValues_VW->SetBranchAddress("R", &R_ratio_VW);
      ratioFitValues_VW->SetBranchAddress("R_err", &R_err_ratio_VW);
      ratioFitValues_VW->SetBranchAddress("phi", &phi_ratio_VW);
      ratioFitValues_VW->SetBranchAddress("phi_err", &phi_err_ratio_VW);



      for (uint iterNum = 0; iterNum < totalIters; ++iterNum)
      {
        fiveParamFitValues->GetEntry(iterNum);
        fiveParamFitValues_VW->GetEntry(iterNum);
        ratioFitValues->GetEntry(iterNum);
        ratioFitValues_VW->GetEntry(iterNum);

        fiveParRs[iterNum][fileNum] = R_fivePar;
        fiveParRerrs[iterNum][fileNum] = R_err_fivePar;

        fiveParRs_VW[iterNum][fileNum] = R_fivePar_VW;
        fiveParRerrs_VW[iterNum][fileNum] = R_err_fivePar_VW;

        ratioRs[iterNum][fileNum] = R_ratio;
        ratioRerrs[iterNum][fileNum] = R_err_ratio;

        ratioRs_VW[iterNum][fileNum] = R_ratio_VW;
        ratioRerrs_VW[iterNum][fileNum] = R_err_ratio_VW;
      } // end iter loop


    inputFile->Close();
  } // end file num loop


  cout << "Five parameter first error: " << fiveParRerrs[0][0] << endl;
  cout << "Ratio first error: " << ratioRerrs[0][0] << endl;

/////////////////////////////////////////////////////////////////////////////////////

  for (uint fileNum = 0; fileNum < fileVector.size(); ++fileNum)
  {
    fiveParRs_Avg[fileNum] = 0;
    fiveParRs_VW_Avg[fileNum] = 0;
    ratioRs_Avg[fileNum] = 0;
    ratioRs_VW_Avg[fileNum] = 0;

    for (uint seedNum = 0; seedNum < totalIters; ++seedNum)
    {
      fiveParRs_Avg[fileNum] += fiveParRs[seedNum][fileNum];
      fiveParRs_VW_Avg[fileNum] += fiveParRs_VW[seedNum][fileNum];
      ratioRs_Avg[fileNum] += ratioRs[seedNum][fileNum];
      ratioRs_VW_Avg[fileNum] += ratioRs_VW[seedNum][fileNum];

    }

      fiveParRs_Avg[fileNum] /= totalIters;
      fiveParRs_VW_Avg[fileNum] /= totalIters;
      ratioRs_Avg[fileNum] /= totalIters;
      ratioRs_VW_Avg[fileNum] /= totalIters;
  }




/////////////////////////////////////////////////////////////////////////////////////


  // make 1D arrays of either single seed R values or average R values, which I can then pass into the constructors of the T graph

  vector<string> types = {"Single T Avg", "Single R Avg", "Single T VW Avg", "Single R VW Avg", "Avg T", "Avg R", "Avg T VW", "Avg R VW"};
  char* labels[types.size()];
  for (uint type = 0; type < types.size(); ++type) labels[type] = Form("%s", types.at(type).c_str());

  double correlations[types.size()][types.size()];
  vector<double*> typePointers = {fiveParRs[0], ratioRs[0], fiveParRs_VW[0], ratioRs_VW[0], fiveParRs_Avg, ratioRs_Avg, fiveParRs_VW_Avg, ratioRs_VW_Avg}; // default to 0 seed and replace in the loop


  for (uint typeA = 0; typeA < types.size(); ++typeA)
  {
    for (uint typeB = 0; typeB < types.size(); ++typeB)
    {

      correlations[typeA][typeB] = 0;

      for (uint seedNum = 0; seedNum < totalIters; ++seedNum)
      {

        typePointers.at(0) = fiveParRs[seedNum];
        typePointers.at(1) = ratioRs[seedNum];
        typePointers.at(2) = fiveParRs_VW[seedNum];
        typePointers.at(3) = ratioRs_VW[seedNum];


        TGraph* corrGraph = new TGraph(fileVector.size(), typePointers.at(typeA), typePointers.at(typeB));
        correlations[typeA][typeB] += corrGraph->GetCorrelationFactor();

        // cout << "Seed: " << seedNum << " type A: " << types.at(typeA) << " typeB: " << types.at(typeB) << " corr: " << corrGraph->GetCorrelationFactor() << endl;

        delete corrGraph;
      }

      correlations[typeA][typeB] /= totalIters;

      cout << " Final | Type A: " << types.at(typeA) << " | Type B: " << types.at(typeB) << " | Corr: " << correlations[typeA][typeB] << endl;

    }
  }



/////////////////////////////////////////////////////////////////////////////////////


  // Draw correlation matrix                                                                                                                                                                                            
  TCanvas* c_corrMat= new TCanvas("c_corrMat","Correlation Matrix",50,10,800,800);

  // gStyle->SetMarkerSize(.7);
  gStyle->SetPaintTextFormat("1.4f");
  gStyle->SetGridStyle(0);
  gStyle->SetOptStat(kFALSE);

  c_corrMat->SetGrid();
  c_corrMat->SetTopMargin(0.1);
  c_corrMat->SetBottomMargin(0.15);
  c_corrMat->SetRightMargin(0.2);
  c_corrMat->SetLeftMargin(0.2);


  TH2D* corrMatrixHist = new TH2D("corrMatrixHist","Correlation Matrix", types.size(), 0, types.size(), types.size(), 0, types.size());

  for(uint i = 0; i < types.size(); i++){
    for(uint j = 0; j < types.size(); j++){
      corrMatrixHist->Fill(labels[i], labels[j], correlations[i][j]);
    }
  }

  corrMatrixHist->LabelsDeflate();
  corrMatrixHist->GetZaxis()->SetRangeUser(0.95,1);
  corrMatrixHist->GetXaxis()->SetLabelSize(0.04);
  corrMatrixHist->GetXaxis()->SetLabelOffset(0.02);
  corrMatrixHist->GetYaxis()->SetLabelSize(0.04);
  corrMatrixHist->SetTitle("");
  corrMatrixHist->Draw("COLZTEXT");

  c_corrMat->SaveAs("TypeCorrelationMatrixPlot.png");

/////////////////////////////////////////////////////////////////////////////////////

return 1;

}
