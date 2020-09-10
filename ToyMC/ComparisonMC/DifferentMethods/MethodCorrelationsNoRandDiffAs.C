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

/////////////////////////////////////////////////////////////////////////////////////

// some dataset error numbers for allowed deviation plots

static vector<double> dataset_singleseed_errors_TMethod = {1.33, 1.13, 0.91, 0.75}; // approximate T method errors for the run 1 datasets
static vector<double> dataset_singleseed_errors_AMethod = {0.9 * 1.33, 0.9 * 1.13, 0.9 * 0.91, 0.9 * 0.75}; // A method errors approximately 10% better than the T method
// static vector<double> dataset_singleseed_errors_RMethod = {1.33, 1.13, 0.91, 0.75}; // R method errors are basically the same as the T method
static vector<double> dataset_singleseed_errors_QMethod = {1.55 * 1.33, 1.55 * 1.13, 1.55 * 0.91, 1.55 * 0.75}; // Q method errors approximately 1.55 times worse than the T method

// vector<vector<double>> datasetErrors = {dataset_singleseed_errors_TMethod, dataset_singleseed_errors_AMethod, dataset_singleseed_errors_RMethod, dataset_singleseed_errors_QMethod};
vector<vector<double>> datasetErrors = {dataset_singleseed_errors_TMethod, dataset_singleseed_errors_AMethod, dataset_singleseed_errors_QMethod};

/////////////////////////////////////////////////////////////////////////////////////

uint totalIters = 0;

double calculateAverage(vector<double> inVec)
{
  double avg = 0;
  for (uint i = 0; i < inVec.size(); ++i) avg += inVec.at(i);
  avg /= inVec.size();
  return avg;
}


int MethodCorrelationsNoRandDiffAs(std::string filePath)
{
  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

/////////////////////////////////////////////////////////////////////////////////////

  // style setting for nice 2D matrix plots

  gStyle->SetGridStyle(0);
  gStyle->SetOptStat(kFALSE);

/////////////////////////////////////////////////////////////////////////////////////

  TFile* inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  TNtuple *tempTree = (TNtuple*)inputFile->Get("topDir/T/TMethodFitValues");
  totalIters = tempTree->GetEntries();
  cout << "Num tree entries: " << totalIters << endl;

  // inputFile->Close();

/////////////////////////////////////////////////////////////////////////////////////

  TFile* outputFile = new TFile("methodFitCorrelationsNoRandDiffAs.root","RECREATE");

  auto histogramsDir = outputFile->mkdir("Histograms");
  auto canvasesDir = outputFile->mkdir("Canvases");
  auto deviationsDir = outputFile->mkdir("Deviations");

/////////////////////////////////////////////////////////////////////////////////////

  // dynamically allocate memory so that the stack isn't overflown - https://www.techiedelight.com/dynamic-memory-allocation-in-c-for-2d-3d-array/

  double** Tmethod_Pars = new double*[5];
  double** Amethod_Pars = new double*[5];
  double** Amethod_Pars_60h = new double*[5];
  double** Amethod_Pars_HK = new double*[5];
  double** Amethod_Pars_9d = new double*[5];
  double** Amethod_Pars_EG = new double*[5];
  double** Qmethod_Pars = new double*[5];

  for (int i = 0; i < 5; ++i)
  {
    Tmethod_Pars[i] = new double[totalIters];
    Amethod_Pars[i] = new double[totalIters];
    Amethod_Pars_60h[i] = new double[totalIters];
    Amethod_Pars_HK[i] = new double[totalIters];
    Amethod_Pars_9d[i] = new double[totalIters];
    Amethod_Pars_EG[i] = new double[totalIters];
    Qmethod_Pars[i] = new double[totalIters];
  }

/////////////////////////////////////////////////////////////////////////////////////


      TNtuple *TMethodFitValues = (TNtuple*)inputFile->Get("topDir/T/TMethodFitValues");
      
      float N_Tm, N_err_Tm, tau_Tm, tau_err_Tm, A_Tm, A_err_Tm, R_Tm, R_err_Tm, phi_Tm, phi_err_Tm;
      TMethodFitValues->SetBranchAddress("N", &N_Tm);
      TMethodFitValues->SetBranchAddress("N_err", &N_err_Tm);
      TMethodFitValues->SetBranchAddress("tau", &tau_Tm);
      TMethodFitValues->SetBranchAddress("tau_err", &tau_err_Tm);
      TMethodFitValues->SetBranchAddress("A", &A_Tm);
      TMethodFitValues->SetBranchAddress("A_err", &A_err_Tm);
      TMethodFitValues->SetBranchAddress("R", &R_Tm);
      TMethodFitValues->SetBranchAddress("R_err", &R_err_Tm);
      TMethodFitValues->SetBranchAddress("phi", &phi_Tm);
      TMethodFitValues->SetBranchAddress("phi_err", &phi_err_Tm);

      TNtuple *AMethodFitValues = (TNtuple*)inputFile->Get("topDir/A/AMethodFitValues");
      
      float N_Am, N_err_Am, tau_Am, tau_err_Am, A_Am, A_err_Am, R_Am, R_err_Am, phi_Am, phi_err_Am;
      AMethodFitValues->SetBranchAddress("N", &N_Am);
      AMethodFitValues->SetBranchAddress("N_err", &N_err_Am);
      AMethodFitValues->SetBranchAddress("tau", &tau_Am);
      AMethodFitValues->SetBranchAddress("tau_err", &tau_err_Am);
      AMethodFitValues->SetBranchAddress("A", &A_Am);
      AMethodFitValues->SetBranchAddress("A_err", &A_err_Am);
      AMethodFitValues->SetBranchAddress("R", &R_Am);
      AMethodFitValues->SetBranchAddress("R_err", &R_err_Am);
      AMethodFitValues->SetBranchAddress("phi", &phi_Am);
      AMethodFitValues->SetBranchAddress("phi_err", &phi_err_Am);


      TNtuple *AMethodFitValues_60h = (TNtuple*)inputFile->Get("topDir/A_60h/AMethodFitValues_60h");
      
      float N_Am_60h, N_err_Am_60h, tau_Am_60h, tau_err_Am_60h, A_Am_60h, A_err_Am_60h, R_Am_60h, R_err_Am_60h, phi_Am_60h, phi_err_Am_60h;
      AMethodFitValues_60h->SetBranchAddress("N", &N_Am_60h);
      AMethodFitValues_60h->SetBranchAddress("N_err", &N_err_Am_60h);
      AMethodFitValues_60h->SetBranchAddress("tau", &tau_Am_60h);
      AMethodFitValues_60h->SetBranchAddress("tau_err", &tau_err_Am_60h);
      AMethodFitValues_60h->SetBranchAddress("A", &A_Am_60h);
      AMethodFitValues_60h->SetBranchAddress("A_err", &A_err_Am_60h);
      AMethodFitValues_60h->SetBranchAddress("R", &R_Am_60h);
      AMethodFitValues_60h->SetBranchAddress("R_err", &R_err_Am_60h);
      AMethodFitValues_60h->SetBranchAddress("phi", &phi_Am_60h);
      AMethodFitValues_60h->SetBranchAddress("phi_err", &phi_err_Am_60h);


      TNtuple *AMethodFitValues_HK = (TNtuple*)inputFile->Get("topDir/A_HK/AMethodFitValues_HK");
      
      float N_Am_HK, N_err_Am_HK, tau_Am_HK, tau_err_Am_HK, A_Am_HK, A_err_Am_HK, R_Am_HK, R_err_Am_HK, phi_Am_HK, phi_err_Am_HK;
      AMethodFitValues_HK->SetBranchAddress("N", &N_Am_HK);
      AMethodFitValues_HK->SetBranchAddress("N_err", &N_err_Am_HK);
      AMethodFitValues_HK->SetBranchAddress("tau", &tau_Am_HK);
      AMethodFitValues_HK->SetBranchAddress("tau_err", &tau_err_Am_HK);
      AMethodFitValues_HK->SetBranchAddress("A", &A_Am_HK);
      AMethodFitValues_HK->SetBranchAddress("A_err", &A_err_Am_HK);
      AMethodFitValues_HK->SetBranchAddress("R", &R_Am_HK);
      AMethodFitValues_HK->SetBranchAddress("R_err", &R_err_Am_HK);
      AMethodFitValues_HK->SetBranchAddress("phi", &phi_Am_HK);
      AMethodFitValues_HK->SetBranchAddress("phi_err", &phi_err_Am_HK);


      TNtuple *AMethodFitValues_9d = (TNtuple*)inputFile->Get("topDir/A_9d/AMethodFitValues_9d");
      
      float N_Am_9d, N_err_Am_9d, tau_Am_9d, tau_err_Am_9d, A_Am_9d, A_err_Am_9d, R_Am_9d, R_err_Am_9d, phi_Am_9d, phi_err_Am_9d;
      AMethodFitValues_9d->SetBranchAddress("N", &N_Am_9d);
      AMethodFitValues_9d->SetBranchAddress("N_err", &N_err_Am_9d);
      AMethodFitValues_9d->SetBranchAddress("tau", &tau_Am_9d);
      AMethodFitValues_9d->SetBranchAddress("tau_err", &tau_err_Am_9d);
      AMethodFitValues_9d->SetBranchAddress("A", &A_Am_9d);
      AMethodFitValues_9d->SetBranchAddress("A_err", &A_err_Am_9d);
      AMethodFitValues_9d->SetBranchAddress("R", &R_Am_9d);
      AMethodFitValues_9d->SetBranchAddress("R_err", &R_err_Am_9d);
      AMethodFitValues_9d->SetBranchAddress("phi", &phi_Am_9d);
      AMethodFitValues_9d->SetBranchAddress("phi_err", &phi_err_Am_9d);


      TNtuple *AMethodFitValues_EG = (TNtuple*)inputFile->Get("topDir/A_EG/AMethodFitValues_EG");
      
      float N_Am_EG, N_err_Am_EG, tau_Am_EG, tau_err_Am_EG, A_Am_EG, A_err_Am_EG, R_Am_EG, R_err_Am_EG, phi_Am_EG, phi_err_Am_EG;
      AMethodFitValues_EG->SetBranchAddress("N", &N_Am_EG);
      AMethodFitValues_EG->SetBranchAddress("N_err", &N_err_Am_EG);
      AMethodFitValues_EG->SetBranchAddress("tau", &tau_Am_EG);
      AMethodFitValues_EG->SetBranchAddress("tau_err", &tau_err_Am_EG);
      AMethodFitValues_EG->SetBranchAddress("A", &A_Am_EG);
      AMethodFitValues_EG->SetBranchAddress("A_err", &A_err_Am_EG);
      AMethodFitValues_EG->SetBranchAddress("R", &R_Am_EG);
      AMethodFitValues_EG->SetBranchAddress("R_err", &R_err_Am_EG);
      AMethodFitValues_EG->SetBranchAddress("phi", &phi_Am_EG);
      AMethodFitValues_EG->SetBranchAddress("phi_err", &phi_err_Am_EG);


      TNtuple *QMethodFitValues = (TNtuple*)inputFile->Get("topDir/Q/QMethodFitValues");
      
      float N_Qm, N_err_Qm, tau_Qm, tau_err_Qm, A_Qm, A_err_Qm, R_Qm, R_err_Qm, phi_Qm, phi_err_Qm;
      QMethodFitValues->SetBranchAddress("N", &N_Qm);
      QMethodFitValues->SetBranchAddress("N_err", &N_err_Qm);
      QMethodFitValues->SetBranchAddress("tau", &tau_Qm);
      QMethodFitValues->SetBranchAddress("tau_err", &tau_err_Qm);
      QMethodFitValues->SetBranchAddress("A", &A_Qm);
      QMethodFitValues->SetBranchAddress("A_err", &A_err_Qm);
      QMethodFitValues->SetBranchAddress("R", &R_Qm);
      QMethodFitValues->SetBranchAddress("R_err", &R_err_Qm);
      QMethodFitValues->SetBranchAddress("phi", &phi_Qm);
      QMethodFitValues->SetBranchAddress("phi_err", &phi_err_Qm);



      for (uint iterNum = 0; iterNum < totalIters; ++iterNum)
      {
        TMethodFitValues->GetEntry(iterNum);
        AMethodFitValues->GetEntry(iterNum);
        AMethodFitValues_60h->GetEntry(iterNum);
        AMethodFitValues_HK->GetEntry(iterNum);
        AMethodFitValues_9d->GetEntry(iterNum);
        AMethodFitValues_EG->GetEntry(iterNum);
        QMethodFitValues->GetEntry(iterNum);

        Tmethod_Pars[0][iterNum] = N_Tm;
        Tmethod_Pars[1][iterNum] = tau_Tm;
        Tmethod_Pars[2][iterNum] = A_Tm;
        Tmethod_Pars[3][iterNum] = R_Tm;
        Tmethod_Pars[4][iterNum] = phi_Tm;

        Amethod_Pars[0][iterNum] = N_Am;
        Amethod_Pars[1][iterNum] = tau_Am;
        Amethod_Pars[2][iterNum] = A_Am;
        Amethod_Pars[3][iterNum] = R_Am;
        Amethod_Pars[4][iterNum] = phi_Am;

        Amethod_Pars_60h[0][iterNum] = N_Am_60h;
        Amethod_Pars_60h[1][iterNum] = tau_Am_60h;
        Amethod_Pars_60h[2][iterNum] = A_Am_60h;
        Amethod_Pars_60h[3][iterNum] = R_Am_60h;
        Amethod_Pars_60h[4][iterNum] = phi_Am_60h;

        Amethod_Pars_HK[0][iterNum] = N_Am_HK;
        Amethod_Pars_HK[1][iterNum] = tau_Am_HK;
        Amethod_Pars_HK[2][iterNum] = A_Am_HK;
        Amethod_Pars_HK[3][iterNum] = R_Am_HK;
        Amethod_Pars_HK[4][iterNum] = phi_Am_HK;

        Amethod_Pars_9d[0][iterNum] = N_Am_9d;
        Amethod_Pars_9d[1][iterNum] = tau_Am_9d;
        Amethod_Pars_9d[2][iterNum] = A_Am_9d;
        Amethod_Pars_9d[3][iterNum] = R_Am_9d;
        Amethod_Pars_9d[4][iterNum] = phi_Am_9d;

        Amethod_Pars_EG[0][iterNum] = N_Am_EG;
        Amethod_Pars_EG[1][iterNum] = tau_Am_EG;
        Amethod_Pars_EG[2][iterNum] = A_Am_EG;
        Amethod_Pars_EG[3][iterNum] = R_Am_EG;
        Amethod_Pars_EG[4][iterNum] = phi_Am_EG;

        Qmethod_Pars[0][iterNum] = N_Qm;
        Qmethod_Pars[1][iterNum] = tau_Qm;
        Qmethod_Pars[2][iterNum] = A_Qm;
        Qmethod_Pars[3][iterNum] = R_Qm;
        Qmethod_Pars[4][iterNum] = phi_Qm;

        if(iterNum == 0){
          cout << "T Method first error: " << R_err_Tm << endl;
          cout << "A Method first error: " << R_err_Am << endl;
          cout << "Q Method first error: " << R_err_Qm << endl;
        }
      } // end iter loop

/////////////////////////////////////////////////////////////////////////////////////

 // construct correlation matrices for each individual parameter among the different methods

  vector<string> parameter_strings = {"N", "#tau", "A", "R", "#phi"};
  vector<string> methodTypes = {"T", "A", "A_60h", "A_HK", "A_9d", "A_EG", "Q"}; // make sure the pointers down below correspond to the correct labels


  for (uint firstPar = 0; firstPar < parameter_strings.size(); ++firstPar)
  {
    for (uint secondPar = firstPar; secondPar < parameter_strings.size(); ++secondPar)
    {

      // tried to put everything below here in a separate method to make things easier to code, but had trouble dealing with 2D and 3D arrays when attempting to pass them to another method

        char* labelsX[methodTypes.size()];
        char* labelsY[methodTypes.size()];
        for (uint type = 0; type < methodTypes.size(); ++type){
          labelsX[type] = Form("%s : %s", methodTypes.at(type).c_str(), parameter_strings.at(firstPar).c_str());
          labelsY[type] = Form("%s : %s", methodTypes.at(type).c_str(), parameter_strings.at(secondPar).c_str());
        } 

        vector<double*> typePointersA = {Tmethod_Pars[firstPar], Amethod_Pars[firstPar], Amethod_Pars_60h[firstPar], Amethod_Pars_HK[firstPar], Amethod_Pars_9d[firstPar], Amethod_Pars_EG[firstPar], Qmethod_Pars[firstPar]};
        vector<double*> typePointersB = {Tmethod_Pars[secondPar], Amethod_Pars[secondPar], Amethod_Pars_60h[secondPar], Amethod_Pars_HK[secondPar], Amethod_Pars_9d[secondPar], Amethod_Pars_EG[secondPar], Qmethod_Pars[secondPar]};

        double correlations[methodTypes.size()][methodTypes.size()];

        double corrMax = -1, corrMin = 1;

        for (uint typeA = 0; typeA < methodTypes.size(); ++typeA)
        {
          for (uint typeB = 0; typeB < methodTypes.size(); ++typeB)
          {
            correlations[typeA][typeB] = 0;

              TGraph* corrGraph = new TGraph(totalIters, typePointersA.at(typeA), typePointersB.at(typeB));
              correlations[typeA][typeB] = corrGraph->GetCorrelationFactor();

              // TCanvas* c_corrMattest = new TCanvas("c_corrMattest","Correlation Matrix",50,10,925,900);
              // corrGraph->Draw();
              // corrGraph->Write("testGraph");

              delete corrGraph;

            if(correlations[typeA][typeB] > corrMax) corrMax = correlations[typeA][typeB];
            if(correlations[typeA][typeB] < corrMin) corrMin = correlations[typeA][typeB];
          }
        }


      /////////////////////////////////////////////////////////////////////////////////////

        // Draw correlation matrix

        TH2D* corrMatrixHist = new TH2D("corrMatrixHist","Correlation Matrix", methodTypes.size(), 0, methodTypes.size(), methodTypes.size(), 0, methodTypes.size());

        for(uint i = 0; i < methodTypes.size(); i++){
          for(uint j = 0; j < methodTypes.size(); j++){
            corrMatrixHist->Fill(labelsX[i], labelsY[j], correlations[i][j]);
          }
        }

        corrMatrixHist->LabelsDeflate();
        corrMatrixHist->GetZaxis()->SetRangeUser(corrMin,corrMax);
        corrMatrixHist->GetXaxis()->SetLabelSize(0.04);
        corrMatrixHist->GetXaxis()->SetLabelOffset(0.02);
        corrMatrixHist->GetZaxis()->SetTitleOffset(1.9);
        corrMatrixHist->GetYaxis()->SetLabelSize(0.04);
        corrMatrixHist->SetTitle("");

        corrMatrixHist->SetName(Form("CorrelationMatrix_%s_%s", parameter_strings.at(firstPar).c_str(), parameter_strings.at(secondPar).c_str()));
        // corrMatrixHist->SetTitle(Form("CorrelationMatrix_%s_%s", parameter_strings.at(firstPar).c_str(), parameter_strings.at(secondPar).c_str()));


        TCanvas* c_corrMat = new TCanvas("c_corrMat","Correlation Matrix",50,10,925,900);

        gStyle->SetPaintTextFormat("1.4f");

        c_corrMat->SetGrid();
        c_corrMat->SetTopMargin(0.1);
        c_corrMat->SetBottomMargin(0.15);
        c_corrMat->SetRightMargin(0.22);
        c_corrMat->SetLeftMargin(0.23);

        corrMatrixHist->Draw("COLZTEXT");

        c_corrMat->SetName(Form("CorrelationMatrix_Plot_%s_%s", parameter_strings.at(firstPar).c_str(), parameter_strings.at(secondPar).c_str()));
        c_corrMat->SetTitle(Form("CorrelationMatrix_Plot_%s_%s", parameter_strings.at(firstPar).c_str(), parameter_strings.at(secondPar).c_str()));

        corrMatrixHist->GetZaxis()->SetTitle("Correlation");
        gPad->Modified();
        gPad->Update();

        c_corrMat->SaveAs(Form("MethodType_CorrelationMatrixPlot_%s_%s.png", parameter_strings.at(firstPar).c_str(), parameter_strings.at(secondPar).c_str()));

        canvasesDir->cd();
        c_corrMat->Write();

        histogramsDir->cd();
        corrMatrixHist->Write();

/////////////////////////////////////////////////////////////////////////////////////

        // produce ppb 1 sigma allowed difference plots for the Run 1 datasets
/*
        if(parameter_strings.at(firstPar).compare("R") == 0 && parameter_strings.at(secondPar).compare("R") == 0)
        {
          for (uint datasetNum = 0; datasetNum < dataset_singleseed_errors_TMethod.size(); ++datasetNum)
          {
            TH2D* ppbDevMatrixHist = new TH2D("name","title", methodTypes.size(), 0, methodTypes.size(), methodTypes.size(), 0, methodTypes.size());

            ppbDevMatrixHist->LabelsDeflate();
            // ppbDevMatrixHist->GetZaxis()->SetRangeUser(corrMin,corrMax);
            ppbDevMatrixHist->GetXaxis()->SetLabelSize(0.04);
            ppbDevMatrixHist->GetXaxis()->SetLabelOffset(0.02);
            ppbDevMatrixHist->GetYaxis()->SetLabelSize(0.04);

            ppbDevMatrixHist->SetName(Form("AllowedDeviation_%s_%s_%s", parameter_strings.at(firstPar).c_str(), parameter_strings.at(secondPar).c_str(), dataset_names.at(datasetNum).c_str()));
            // ppbDevMatrixHist->SetTitle(Form("AllowedDeviation_%s_%s_%s", parameter_strings.at(firstPar).c_str(), parameter_strings.at(secondPar).c_str(), dataset_names.at(datasetNum).c_str()));
            ppbDevMatrixHist->SetTitle("");

            for(uint i = 0; i < methodTypes.size(); i++){
              for(uint j = 0; j < methodTypes.size(); j++){

                double errorA = datasetErrors.at(i%4).at(datasetNum);
                double errorB = datasetErrors.at(j%4).at(datasetNum);

                double squareSigma = pow(errorA,2) + pow(errorB,2) - 2 * correlations[i][j] * errorA * errorB;
                if(abs(squareSigma) < 0.001) squareSigma = 0;

                double oneSigma = sqrt(squareSigma);

                // cout << "Dataset: " << dataset_names.at(datasetNum) <<  " Method type i: " << methodTypes.at(i) << " j: " << methodTypes.at(j) << " with errors: " << errorA << " : " << errorB << " with correlation coefficient " << correlations[i][j] << " allowed deviation: " << oneSigma << " square of it: " << (pow(errorA,2) + pow(errorB,2) - 2 * correlations[i][j] * errorA * errorB) << endl;

                ppbDevMatrixHist->Fill(labelsX[i], labelsY[j], oneSigma);
              }
            }

            deviationsDir->cd();

            ppbDevMatrixHist->Write();

            TCanvas* c_devMat = new TCanvas("c_devMat","1#sigma Deviation",50,10,925,900);

            c_devMat->SetGrid();
            c_devMat->SetTopMargin(0.1);
            c_devMat->SetBottomMargin(0.15);
            c_devMat->SetRightMargin(0.22);
            c_devMat->SetLeftMargin(0.23);

            c_devMat->SetName(Form("AllowedDeviation_Plot_%s_%s_%s", parameter_strings.at(firstPar).c_str(), parameter_strings.at(secondPar).c_str(),dataset_names.at(datasetNum).c_str()));
            c_devMat->SetTitle(Form("AllowedDeviation_Plot_%s_%s", parameter_strings.at(firstPar).c_str(), parameter_strings.at(secondPar).c_str()));

            ppbDevMatrixHist->Draw("COLZTEXT");

            ppbDevMatrixHist->GetZaxis()->SetTitleOffset(1.9);
            ppbDevMatrixHist->GetZaxis()->SetTitle(Form("%s 1#sigma (ppm)", dataset_names.at(datasetNum).c_str()));
            gPad->Modified();
            gPad->Update();

            c_devMat->SaveAs(Form("MethodType_AllowedDeviation_Plot_%s_%s_%s.png", parameter_strings.at(firstPar).c_str(), parameter_strings.at(secondPar).c_str(),dataset_names.at(datasetNum).c_str()));
            c_devMat->Write();

            delete c_devMat;
            delete ppbDevMatrixHist;
          } // end dataset loop
        } // end R R std deviation plot construction
*/

/////////////////////////////////////////////////////////////////////////////////////

        delete c_corrMat;
        delete corrMatrixHist;

    } // end loop over parameter 2
  } // end loop over parameter 1


/////////////////////////////////////////////////////////////////////////////////////

return 1;

}
