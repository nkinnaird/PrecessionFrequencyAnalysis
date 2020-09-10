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
static vector<double> dataset_singleseed_errors_RMethod = {1.33, 1.13, 0.91, 0.75}; // R method errors are basically the same as the T method
static vector<double> dataset_singleseed_errors_QMethod = {1.55 * 1.33, 1.55 * 1.13, 1.55 * 0.91, 1.55 * 0.75}; // Q method errors approximately 1.55 times worse than the T method

vector<vector<double>> datasetErrors = {dataset_singleseed_errors_TMethod, dataset_singleseed_errors_AMethod, dataset_singleseed_errors_RMethod, dataset_singleseed_errors_QMethod};

/////////////////////////////////////////////////////////////////////////////////////

uint totalIters = 0;
uint totalFiles = 0;

double calculateAverage(vector<double> inVec)
{
  double avg = 0;
  for (uint i = 0; i < inVec.size(); ++i) avg += inVec.at(i);
  avg /= inVec.size();
  return avg;
}


int MethodCorrelations(std::string fileList)
{
  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

/////////////////////////////////////////////////////////////////////////////////////

  // style setting for nice 2D matrix plots

  gStyle->SetGridStyle(0);
  gStyle->SetOptStat(kFALSE);

/////////////////////////////////////////////////////////////////////////////////////

  std::vector<string> fileVector;

  std::ifstream inList(fileList);
  string path;
  while(inList >> path) fileVector.push_back(path);

  totalFiles = fileVector.size();
  cout << "Num files: " << totalFiles << endl;

/////////////////////////////////////////////////////////////////////////////////////

  TFile* firstFile = TFile::Open(fileVector.at(0).c_str());
   if (firstFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  TNtuple *tempTree = (TNtuple*)firstFile->Get("topDir/T/TMethodFitValues");
  totalIters = tempTree->GetEntries();
  cout << "Num tree entries: " << totalIters << endl;

  firstFile->Close();

/////////////////////////////////////////////////////////////////////////////////////

  TFile* outputFile = new TFile("methodFitCorrelations.root","RECREATE");

  auto histogramsDir = outputFile->mkdir("Histograms");
  auto canvasesDir = outputFile->mkdir("Canvases");
  auto deviationsDir = outputFile->mkdir("Deviations");

/////////////////////////////////////////////////////////////////////////////////////

  // dynamically allocate memory so that the stack isn't overflown - https://www.techiedelight.com/dynamic-memory-allocation-in-c-for-2d-3d-array/

  double*** Tmethod_Pars = new double**[5];
  double*** Amethod_Pars = new double**[5];
  double*** Rmethod_Pars = new double**[5];
  double*** Qmethod_Pars = new double**[5];

  double** Tmethod_Pars_Avg = new double*[5];
  double** Amethod_Pars_Avg = new double*[5];
  double** Rmethod_Pars_Avg = new double*[5];
  double** Qmethod_Pars_Avg = new double*[5];

  for (int i = 0; i < 5; ++i)
  {
    Tmethod_Pars[i] = new double*[totalIters];
    Amethod_Pars[i] = new double*[totalIters];
    Rmethod_Pars[i] = new double*[totalIters];
    Qmethod_Pars[i] = new double*[totalIters];

    Tmethod_Pars_Avg[i] = new double[totalFiles];
    Amethod_Pars_Avg[i] = new double[totalFiles];
    Rmethod_Pars_Avg[i] = new double[totalFiles];
    Qmethod_Pars_Avg[i] = new double[totalFiles];

    for (uint j = 0; j < totalIters; ++j)
    {
      Tmethod_Pars[i][j] = new double[totalFiles];
      Amethod_Pars[i][j] = new double[totalFiles];
      Rmethod_Pars[i][j] = new double[totalFiles];
      Qmethod_Pars[i][j] = new double[totalFiles];
    }
  }

/////////////////////////////////////////////////////////////////////////////////////


  for (uint fileNum = 0; fileNum < totalFiles; ++fileNum)
  {
    TFile* inputFile = TFile::Open(fileVector.at(fileNum).c_str());
     if (inputFile == 0) {
        printf("Error: cannot open file\n");
        return 0;
     }

    cout << "Filenum: " << fileNum << " " << fileVector.at(fileNum) << endl;

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


      TNtuple *RMethodFitValues = (TNtuple*)inputFile->Get("topDir/R/ratioFitValues");
      
      float A_Rm, A_err_Rm, R_Rm, R_err_Rm, phi_Rm, phi_err_Rm;
      RMethodFitValues->SetBranchAddress("A", &A_Rm);
      RMethodFitValues->SetBranchAddress("A_err", &A_err_Rm);
      RMethodFitValues->SetBranchAddress("R", &R_Rm);
      RMethodFitValues->SetBranchAddress("R_err", &R_err_Rm);
      RMethodFitValues->SetBranchAddress("phi", &phi_Rm);
      RMethodFitValues->SetBranchAddress("phi_err", &phi_err_Rm);


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
        RMethodFitValues->GetEntry(iterNum);
        QMethodFitValues->GetEntry(iterNum);

        Tmethod_Pars[0][iterNum][fileNum] = N_Tm;
        Tmethod_Pars[1][iterNum][fileNum] = tau_Tm;
        Tmethod_Pars[2][iterNum][fileNum] = A_Tm;
        Tmethod_Pars[3][iterNum][fileNum] = R_Tm;
        Tmethod_Pars[4][iterNum][fileNum] = phi_Tm;

        Amethod_Pars[0][iterNum][fileNum] = N_Am;
        Amethod_Pars[1][iterNum][fileNum] = tau_Am;
        Amethod_Pars[2][iterNum][fileNum] = A_Am;
        Amethod_Pars[3][iterNum][fileNum] = R_Am;
        Amethod_Pars[4][iterNum][fileNum] = phi_Am;

        Rmethod_Pars[0][iterNum][fileNum] = 0;
        Rmethod_Pars[1][iterNum][fileNum] = 0;
        Rmethod_Pars[2][iterNum][fileNum] = A_Rm;
        Rmethod_Pars[3][iterNum][fileNum] = R_Rm;
        Rmethod_Pars[4][iterNum][fileNum] = phi_Rm;

        Qmethod_Pars[0][iterNum][fileNum] = N_Qm;
        Qmethod_Pars[1][iterNum][fileNum] = tau_Qm;
        Qmethod_Pars[2][iterNum][fileNum] = A_Qm;
        Qmethod_Pars[3][iterNum][fileNum] = R_Qm;
        Qmethod_Pars[4][iterNum][fileNum] = phi_Qm;

        if(fileNum == 0 && iterNum == 0){
          cout << "T Method first error: " << R_err_Tm << endl;
          cout << "A Method first error: " << R_err_Am << endl;
          cout << "R Method first error: " << R_err_Rm << endl;
          cout << "Q Method first error: " << R_err_Qm << endl;
        }
      } // end iter loop

    inputFile->Close();
  } // end file num loop


/////////////////////////////////////////////////////////////////////////////////////

  // calculate average parameter values

  for (uint fileNum = 0; fileNum < totalFiles; ++fileNum)
  {
    for (int parNum = 0; parNum < 5; ++parNum)
    {
      Tmethod_Pars_Avg[parNum][fileNum] = 0;
      Amethod_Pars_Avg[parNum][fileNum] = 0;
      Rmethod_Pars_Avg[parNum][fileNum] = 0;
      Qmethod_Pars_Avg[parNum][fileNum] = 0;

      for (uint seedNum = 0; seedNum < totalIters; ++seedNum)
      {
        Tmethod_Pars_Avg[parNum][fileNum] += Tmethod_Pars[parNum][seedNum][fileNum];
        Amethod_Pars_Avg[parNum][fileNum] += Amethod_Pars[parNum][seedNum][fileNum];
        Rmethod_Pars_Avg[parNum][fileNum] += Rmethod_Pars[parNum][seedNum][fileNum];
        Qmethod_Pars_Avg[parNum][fileNum] += Qmethod_Pars[parNum][seedNum][fileNum];          
      }   

      Tmethod_Pars_Avg[parNum][fileNum] /= totalIters;
      Amethod_Pars_Avg[parNum][fileNum] /= totalIters;
      Rmethod_Pars_Avg[parNum][fileNum] /= totalIters;
      Qmethod_Pars_Avg[parNum][fileNum] /= totalIters;

    } // end parameter number loop

  } // end file vector loop

/////////////////////////////////////////////////////////////////////////////////////

 // construct correlation matrices for each individual parameter among the different methods

  vector<string> parameter_strings = {"N", "#tau", "A", "R", "#phi"};
  vector<string> methodTypes = {"T_{single seed avg}", "A_{single seed avg}", "R_{single seed avg}", "Q_{single seed avg}", "T_{Avg}", "A_{Avg}", "R_{Avg}", "Q_{Avg}"}; // make sure the pointers down below correspond to the correct labels


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

        vector<double*> typePointersA = {Tmethod_Pars[firstPar][0], Amethod_Pars[firstPar][0], Rmethod_Pars[firstPar][0], Qmethod_Pars[firstPar][0], Tmethod_Pars_Avg[firstPar], Amethod_Pars_Avg[firstPar], Rmethod_Pars_Avg[firstPar], Qmethod_Pars_Avg[firstPar]}; // default to 0 seed and replace in the loop
        vector<double*> typePointersB = {Tmethod_Pars[secondPar][0], Amethod_Pars[secondPar][0], Rmethod_Pars[secondPar][0], Qmethod_Pars[secondPar][0], Tmethod_Pars_Avg[secondPar], Amethod_Pars_Avg[secondPar], Rmethod_Pars_Avg[secondPar], Qmethod_Pars_Avg[secondPar]}; // default to 0 seed and replace in the loop

        double correlations[methodTypes.size()][methodTypes.size()];

        double corrMax = -1, corrMin = 1;

        for (uint typeA = 0; typeA < methodTypes.size(); ++typeA)
        {
          for (uint typeB = 0; typeB < methodTypes.size(); ++typeB)
          {
            correlations[typeA][typeB] = 0;

            for (uint seedNum = 0; seedNum < totalIters; ++seedNum)
            {
              typePointersA.at(0) = Tmethod_Pars[firstPar][seedNum];
              typePointersA.at(1) = Amethod_Pars[firstPar][seedNum];
              typePointersA.at(2) = Rmethod_Pars[firstPar][seedNum];
              typePointersA.at(3) = Qmethod_Pars[firstPar][seedNum];

              typePointersB.at(0) = Tmethod_Pars[secondPar][seedNum];
              typePointersB.at(1) = Amethod_Pars[secondPar][seedNum];
              typePointersB.at(2) = Rmethod_Pars[secondPar][seedNum];
              typePointersB.at(3) = Qmethod_Pars[secondPar][seedNum];


              TGraph* corrGraph = new TGraph(totalFiles, typePointersA.at(typeA), typePointersB.at(typeB));
              correlations[typeA][typeB] += corrGraph->GetCorrelationFactor();

              delete corrGraph;
            }

            correlations[typeA][typeB] /= totalIters;

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


/////////////////////////////////////////////////////////////////////////////////////

        delete c_corrMat;
        delete corrMatrixHist;

    } // end loop over parameter 2
  } // end loop over parameter 1


/////////////////////////////////////////////////////////////////////////////////////

return 1;

}
