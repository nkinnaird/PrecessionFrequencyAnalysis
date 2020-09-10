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

static vector<double> dataset_singleseed_errors_TMethod = {1.335, 1.132, 0.911, 0.746};
static vector<double> dataset_singleseed_errors_RMethod = {1.337, 1.133, 0.914, 0.745};
static vector<double> dataset_singleseed_errors_TMethod_VW = {1.358, 1.156, 0.930, 0.758};
static vector<double> dataset_singleseed_errors_RMethod_VW = {1.360, 1.157, 0.933, 0.758};

vector<vector<double>> datasetErrors = {dataset_singleseed_errors_TMethod, dataset_singleseed_errors_RMethod, dataset_singleseed_errors_TMethod_VW, dataset_singleseed_errors_RMethod_VW};

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


int WWoVWRandCorrelationAll(std::string fileList)
{
  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

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

  TNtuple *tempTree = (TNtuple*)firstFile->Get("topDir/ToyMC/Hist/fiveParameterFitValues");
  totalIters = tempTree->GetEntries();
  cout << "Num tree entries: " << totalIters << endl;

  firstFile->Close();
      
/////////////////////////////////////////////////////////////////////////////////////

  TFile* outputFile = new TFile("WWoVWRandCorrelations.root","RECREATE");

  auto histogramsDir = outputFile->mkdir("Histograms");
  auto canvasesDir = outputFile->mkdir("Canvases");
  auto deviationsDir = outputFile->mkdir("Deviations");

/////////////////////////////////////////////////////////////////////////////////////

  // dynamically allocate memory so that the stack isn't overflown - https://www.techiedelight.com/dynamic-memory-allocation-in-c-for-2d-3d-array/

  double*** fivePar_Pars = new double**[5];
  double*** fivePar_Pars_VW = new double**[5];
  double*** ratio_Pars = new double**[5];
  double*** ratio_Pars_VW = new double**[5];

  double** fivePar_Pars_Avg = new double*[5];
  double** fivePar_Pars_VW_Avg = new double*[5];
  double** ratio_Pars_Avg = new double*[5];
  double** ratio_Pars_VW_Avg = new double*[5];

  for (int i = 0; i < 5; ++i)
  {
    fivePar_Pars[i] = new double*[totalIters];
    fivePar_Pars_VW[i] = new double*[totalIters];
    ratio_Pars[i] = new double*[totalIters];
    ratio_Pars_VW[i] = new double*[totalIters];

    fivePar_Pars_Avg[i] = new double[totalFiles];
    fivePar_Pars_VW_Avg[i] = new double[totalFiles];
    ratio_Pars_Avg[i] = new double[totalFiles];
    ratio_Pars_VW_Avg[i] = new double[totalFiles];

    for (uint j = 0; j < totalIters; ++j)
    {
      fivePar_Pars[i][j] = new double[totalFiles];
      fivePar_Pars_VW[i][j] = new double[totalFiles];
      ratio_Pars[i][j] = new double[totalFiles];
      ratio_Pars_VW[i][j] = new double[totalFiles];
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

        fivePar_Pars[0][iterNum][fileNum] = N_fivePar;
        fivePar_Pars[1][iterNum][fileNum] = tau_fivePar;
        fivePar_Pars[2][iterNum][fileNum] = A_fivePar;
        fivePar_Pars[3][iterNum][fileNum] = R_fivePar;
        fivePar_Pars[4][iterNum][fileNum] = phi_fivePar;

        fivePar_Pars_VW[0][iterNum][fileNum] = N_fivePar_VW;
        fivePar_Pars_VW[1][iterNum][fileNum] = tau_fivePar_VW;
        fivePar_Pars_VW[2][iterNum][fileNum] = A_fivePar_VW;
        fivePar_Pars_VW[3][iterNum][fileNum] = R_fivePar_VW;
        fivePar_Pars_VW[4][iterNum][fileNum] = phi_fivePar_VW;

        ratio_Pars[0][iterNum][fileNum] = 0;
        ratio_Pars[1][iterNum][fileNum] = 0;
        ratio_Pars[2][iterNum][fileNum] = A_ratio;
        ratio_Pars[3][iterNum][fileNum] = R_ratio;
        ratio_Pars[4][iterNum][fileNum] = phi_ratio;

        ratio_Pars_VW[0][iterNum][fileNum] = 0;
        ratio_Pars_VW[1][iterNum][fileNum] = 0;
        ratio_Pars_VW[2][iterNum][fileNum] = A_ratio_VW;
        ratio_Pars_VW[3][iterNum][fileNum] = R_ratio_VW;
        ratio_Pars_VW[4][iterNum][fileNum] = phi_ratio_VW;

        if(fileNum == 0 && iterNum == 0){
          cout << "Five parameter first error: " << R_err_fivePar << endl;
          cout << "Ratio first error: " << R_err_ratio << endl;
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
      fivePar_Pars_Avg[parNum][fileNum] = 0;
      fivePar_Pars_VW_Avg[parNum][fileNum] = 0;
      ratio_Pars_Avg[parNum][fileNum] = 0;
      ratio_Pars_VW_Avg[parNum][fileNum] = 0;

      for (uint seedNum = 0; seedNum < totalIters; ++seedNum)
      {
        fivePar_Pars_Avg[parNum][fileNum] += fivePar_Pars[parNum][seedNum][fileNum];
        fivePar_Pars_VW_Avg[parNum][fileNum] += fivePar_Pars_VW[parNum][seedNum][fileNum];
        ratio_Pars_Avg[parNum][fileNum] += ratio_Pars[parNum][seedNum][fileNum];
        ratio_Pars_VW_Avg[parNum][fileNum] += ratio_Pars_VW[parNum][seedNum][fileNum];          
      }   

      fivePar_Pars_Avg[parNum][fileNum] /= totalIters;
      fivePar_Pars_VW_Avg[parNum][fileNum] /= totalIters;
      ratio_Pars_Avg[parNum][fileNum] /= totalIters;
      ratio_Pars_VW_Avg[parNum][fileNum] /= totalIters;

    } // end parameter number loop

  } // end file vector loop

/////////////////////////////////////////////////////////////////////////////////////

 // construct correlation matrices for each individual parameter among the different methods

  vector<string> parameter_strings = {"N", "#tau", "A", "R", "#phi"};
  vector<string> methodTypes = {"T_{single seed avg}", "R_{single seed avg}", "T_{VW single seed avg}", "R_{VW single seed avg}", "T_{Avg}", "R_{Avg}", "T_{VW Avg}", "R_{VW Avg}"}; // make sure the pointers down below correspond to the correct labels


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

        vector<double*> typePointersA = {fivePar_Pars[firstPar][0], ratio_Pars[firstPar][0], fivePar_Pars_VW[firstPar][0], ratio_Pars_VW[firstPar][0], fivePar_Pars_Avg[firstPar], ratio_Pars_Avg[firstPar], fivePar_Pars_VW_Avg[firstPar], ratio_Pars_VW_Avg[firstPar]}; // default to 0 seed and replace in the loop
        vector<double*> typePointersB = {fivePar_Pars[secondPar][0], ratio_Pars[secondPar][0], fivePar_Pars_VW[secondPar][0], ratio_Pars_VW[secondPar][0], fivePar_Pars_Avg[secondPar], ratio_Pars_Avg[secondPar], fivePar_Pars_VW_Avg[secondPar], ratio_Pars_VW_Avg[secondPar]}; // default to 0 seed and replace in the loop

        double correlations[methodTypes.size()][methodTypes.size()];

        double corrMax = -1, corrMin = 1;

        for (uint typeA = 0; typeA < methodTypes.size(); ++typeA)
        {
          for (uint typeB = 0; typeB < methodTypes.size(); ++typeB)
          {
            correlations[typeA][typeB] = 0;

            for (uint seedNum = 0; seedNum < totalIters; ++seedNum)
            {
              typePointersA.at(0) = fivePar_Pars[firstPar][seedNum];
              typePointersA.at(1) = ratio_Pars[firstPar][seedNum];
              typePointersA.at(2) = fivePar_Pars_VW[firstPar][seedNum];
              typePointersA.at(3) = ratio_Pars_VW[firstPar][seedNum];

              typePointersB.at(0) = fivePar_Pars[secondPar][seedNum];
              typePointersB.at(1) = ratio_Pars[secondPar][seedNum];
              typePointersB.at(2) = fivePar_Pars_VW[secondPar][seedNum];
              typePointersB.at(3) = ratio_Pars_VW[secondPar][seedNum];


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
        gStyle->SetGridStyle(0);
        gStyle->SetOptStat(kFALSE);

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
