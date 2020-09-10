#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TImage.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TVectorD.h>
#include <TFitResult.h>
#include <TText.h>
#include <TPaveStats.h>

// #include "plotUtils.hh"

using namespace std;

uint totalSamples = 0;

/////////////////////////////////////////////////////////////////////////////////////

// some dataset error numbers for allowed deviation plots

// static vector<double> dataset_singleseed_errors_TMethod = {1.33, 1.13, 0.91, 0.75}; // approximate T method errors for the run 1 datasets
// static vector<double> dataset_singleseed_errors_AMethod = {0.9 * 1.33, 0.9 * 1.13, 0.9 * 0.91, 0.9 * 0.75}; // A method errors approximately 10% better than the T method
// static vector<double> dataset_singleseed_errors_RMethod = {1.33, 1.13, 0.91, 0.75}; // R method errors are basically the same as the T method
// static vector<double> dataset_singleseed_errors_QMethod = {1.55 * 1.33, 1.55 * 1.13, 1.55 * 0.91, 1.55 * 0.75}; // Q method errors approximately 1.55 times worse than the T method

// vector<vector<double>> datasetErrors = {dataset_singleseed_errors_TMethod, dataset_singleseed_errors_AMethod, dataset_singleseed_errors_RMethod, dataset_singleseed_errors_QMethod};
// vector<vector<double>> datasetErrors = {dataset_singleseed_errors_TMethod, dataset_singleseed_errors_AMethod, dataset_singleseed_errors_QMethod};

/////////////////////////////////////////////////////////////////////////////////////

void readInTree(TNtuple* inTree, double** parameters){

  float N, tau, A, R, phi;
  inTree->SetBranchAddress("N", &N);
  inTree->SetBranchAddress("tau", &tau);
  inTree->SetBranchAddress("A", &A);
  inTree->SetBranchAddress("R", &R);
  inTree->SetBranchAddress("phi", &phi);


  for (uint sampleNum = 0; sampleNum < totalSamples; ++sampleNum)
  {
    inTree->GetEntry(sampleNum);

    parameters[0][sampleNum] = N;
    parameters[1][sampleNum] = tau;
    parameters[2][sampleNum] = A;
    parameters[3][sampleNum] = R;
    parameters[4][sampleNum] = phi;
  } // end loop over samples
}


int calculateAnalyzerCorrelations(std::string filePath)
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

  TNtuple *tempTree = (TNtuple*)inputFile->Get("Nick/nick_T_fitParams");
  totalSamples = tempTree->GetEntries();
  cout << "Num samples: " << totalSamples << endl;

/////////////////////////////////////////////////////////////////////////////////////

  TFile* outputFile = new TFile("analyzerCorrelations.root","RECREATE");

  auto graphsDir = outputFile->mkdir("Graphs");
  auto histogramsDir = outputFile->mkdir("Histograms");
  auto canvasesDir = outputFile->mkdir("Canvases");
  auto deviationsDir = outputFile->mkdir("Deviations");

/////////////////////////////////////////////////////////////////////////////////////

  // dynamically allocate memory so that the stack isn't overflown - https://www.techiedelight.com/dynamic-memory-allocation-in-c-for-2d-3d-array/

  double** nick_T_pars = new double*[5]; // 5 different parameters
  double** nick_R_pars = new double*[5]; 
  double** david_T_pars = new double*[5]; 
  double** david_A_pars = new double*[5]; 
  double** aaron_T_pars = new double*[5]; 
  double** aaron_A_pars = new double*[5]; 
  double** matteo_T_pars = new double*[5]; 
  double** matteo_A_pars = new double*[5]; 
  double** bingzhi_T_pars = new double*[5]; 
  double** bingzhi_A_pars = new double*[5]; 
  double** tim_Q_pars = new double*[5]; 

  // averaged and combined parameters

  double** RW_T_pars = new double*[5]; 
  double** RW_A_pars = new double*[5];

  // double ** RE_TA_pars = new double*[5];

  // double ** RW_TA_pars = new double*[5];
  // double ** RW_TR_pars = new double*[5];
  // double ** RW_AR_pars = new double*[5];
  // double ** RW_TAR_pars = new double*[5];


  for (int i = 0; i < 5; ++i)
  {
    nick_T_pars[i] = new double[totalSamples];
    nick_R_pars[i] = new double[totalSamples];
    david_T_pars[i] = new double[totalSamples];
    david_A_pars[i] = new double[totalSamples];
    aaron_T_pars[i] = new double[totalSamples];
    aaron_A_pars[i] = new double[totalSamples];
    matteo_T_pars[i] = new double[totalSamples];
    matteo_A_pars[i] = new double[totalSamples];
    bingzhi_T_pars[i] = new double[totalSamples];
    bingzhi_A_pars[i] = new double[totalSamples];
    tim_Q_pars[i] = new double[totalSamples];

    RW_T_pars[i] = new double[totalSamples];
    RW_A_pars[i] = new double[totalSamples];

    // RE_TA_pars[i] = new double[totalSamples];

    // RW_TA_pars[i] = new double[totalSamples];
    // RW_TR_pars[i] = new double[totalSamples];
    // RW_AR_pars[i] = new double[totalSamples];
    // RW_TAR_pars[i] = new double[totalSamples];
  }

/////////////////////////////////////////////////////////////////////////////////////

      readInTree((TNtuple*)inputFile->Get("Nick/nick_T_fitParams"), nick_T_pars);
      readInTree((TNtuple*)inputFile->Get("David/david_T_fitParams"), david_T_pars);
      readInTree((TNtuple*)inputFile->Get("David/david_A_fitParams"), david_A_pars);
      readInTree((TNtuple*)inputFile->Get("Aaron/aaron_T_fitParams"), aaron_T_pars);
      readInTree((TNtuple*)inputFile->Get("Aaron/aaron_A_fitParams"), aaron_A_pars);
      readInTree((TNtuple*)inputFile->Get("Matteo/matteo_T_fitParams"), matteo_T_pars);
      readInTree((TNtuple*)inputFile->Get("Matteo/matteo_A_fitParams"), matteo_A_pars);
      readInTree((TNtuple*)inputFile->Get("Bingzhi/bingzhi_T_fitParams"), bingzhi_T_pars);
      readInTree((TNtuple*)inputFile->Get("Bingzhi/bingzhi_A_fitParams"), bingzhi_A_pars);
      readInTree((TNtuple*)inputFile->Get("Tim/tim_Q_fitParams"), tim_Q_pars);

        // do R method separately 

        TNtuple *nick_R_fitParams = (TNtuple*)inputFile->Get("Nick/nick_R_fitParams_average");
        
        float A_nick_R, R_nick_R, phi_nick_R;
        nick_R_fitParams->SetBranchAddress("A", &A_nick_R);
        nick_R_fitParams->SetBranchAddress("R", &R_nick_R);
        nick_R_fitParams->SetBranchAddress("phi", &phi_nick_R);

        for (uint sampleNum = 0; sampleNum < totalSamples; ++sampleNum)
        {
          nick_R_fitParams->GetEntry(sampleNum);

          nick_R_pars[0][sampleNum] = 0;
          nick_R_pars[1][sampleNum] = 0;
          nick_R_pars[2][sampleNum] = A_nick_R;
          nick_R_pars[3][sampleNum] = R_nick_R;
          nick_R_pars[4][sampleNum] = phi_nick_R;
        } // end iter loop

/////////////////////////////////////////////////////////////////////////////////////

 // make average/combined parameters

  for (uint sampleNum = 0; sampleNum < totalSamples; ++sampleNum)
  {
    for (int parNum = 0; parNum < 5; ++parNum)
    {
      RW_T_pars[parNum][sampleNum] = (nick_T_pars[parNum][sampleNum] + aaron_T_pars[parNum][sampleNum] + matteo_T_pars[parNum][sampleNum] + bingzhi_T_pars[parNum][sampleNum])/4.;
      RW_A_pars[parNum][sampleNum] = (aaron_A_pars[parNum][sampleNum] + matteo_A_pars[parNum][sampleNum] + bingzhi_A_pars[parNum][sampleNum])/3.;

      // RE_TA_pars[parNum][sampleNum] = (david_T_pars[parNum][sampleNum] + david_A_pars[parNum][sampleNum])/2.;

      // RW_TA_pars[parNum][sampleNum] = (nick_T_pars[parNum][sampleNum] + aaron_T_pars[parNum][sampleNum] + matteo_T_pars[parNum][sampleNum] + bingzhi_T_pars[parNum][sampleNum] + aaron_A_pars[parNum][sampleNum] + matteo_A_pars[parNum][sampleNum] + bingzhi_A_pars[parNum][sampleNum])/7.;
      // RW_TR_pars[parNum][sampleNum] = (nick_T_pars[parNum][sampleNum] + aaron_T_pars[parNum][sampleNum] + matteo_T_pars[parNum][sampleNum] + bingzhi_T_pars[parNum][sampleNum] + nick_R_pars[parNum][sampleNum])/5.;
      // RW_AR_pars[parNum][sampleNum] = (aaron_A_pars[parNum][sampleNum] + matteo_A_pars[parNum][sampleNum] + bingzhi_A_pars[parNum][sampleNum] + nick_R_pars[parNum][sampleNum])/4.;
      // RW_TAR_pars[parNum][sampleNum] = (nick_T_pars[parNum][sampleNum] + aaron_T_pars[parNum][sampleNum] + matteo_T_pars[parNum][sampleNum] + bingzhi_T_pars[parNum][sampleNum] + aaron_A_pars[parNum][sampleNum] + matteo_A_pars[parNum][sampleNum] + bingzhi_A_pars[parNum][sampleNum] + nick_R_pars[parNum][sampleNum])/8.;
    }
  }


/////////////////////////////////////////////////////////////////////////////////////

 // construct correlation matrices for each individual parameter among the different methods

  vector<string> parameter_strings = {"N", "#tau", "A", "R", "#phi"};
  vector<string> methodTypes = {"Nick_T", "Nick_R", "David_T", "David_A", "Aaron_T", "Aaron_A", "Matteo_T", "Matteo_A", "Bingzhi_T", "Bingzhi_A", "Tim_Q"}; // make sure the pointers down below correspond to the correct labels

  for (uint firstPar = 0; firstPar < parameter_strings.size(); ++firstPar)
  {
    for (uint secondPar = firstPar; secondPar < parameter_strings.size(); ++secondPar)
    {

      // tried to put everything below here in a separate method to make things easier to code, but had trouble dealing with 2D and 3D arrays when attempting to pass them to another method

        vector<double*> typePointersA = {nick_T_pars[firstPar], nick_R_pars[firstPar], david_T_pars[firstPar], david_A_pars[firstPar], aaron_T_pars[firstPar], aaron_A_pars[firstPar], matteo_T_pars[firstPar], matteo_A_pars[firstPar], bingzhi_T_pars[firstPar], bingzhi_A_pars[firstPar], tim_Q_pars[firstPar]};
        vector<double*> typePointersB = {nick_T_pars[secondPar], nick_R_pars[secondPar], david_T_pars[secondPar], david_A_pars[secondPar], aaron_T_pars[secondPar], aaron_A_pars[secondPar], matteo_T_pars[secondPar], matteo_A_pars[secondPar], bingzhi_T_pars[secondPar], bingzhi_A_pars[secondPar], tim_Q_pars[secondPar]};

        double correlations[methodTypes.size()][methodTypes.size()];

        double corrMax = -1, corrMin = 1;

        for (uint typeA = 0; typeA < methodTypes.size(); ++typeA)
        {
          for (uint typeB = 0; typeB < methodTypes.size(); ++typeB)
          {
            correlations[typeA][typeB] = 0;

              TGraph* corrGraph = new TGraph(totalSamples, typePointersA.at(typeA), typePointersB.at(typeB));
              correlations[typeA][typeB] = corrGraph->GetCorrelationFactor();

              graphsDir->cd();

              corrGraph->SetTitle(Form(";%s_%s;%s_%s", methodTypes.at(typeA).c_str(), parameter_strings.at(firstPar).c_str(), methodTypes.at(typeB).c_str(), parameter_strings.at(secondPar).c_str()));
              corrGraph->Write(Form("Graph_%s_%s_vs_%s_%s", methodTypes.at(typeA).c_str(), parameter_strings.at(firstPar).c_str(), methodTypes.at(typeB).c_str(), parameter_strings.at(secondPar).c_str()));

              // TCanvas* c_corrMattest = new TCanvas("c_corrMattest","Correlation Matrix",50,10,925,900);
              // corrGraph->Draw();
              // c_corrMattest->SaveAs();

              delete corrGraph;

            if(correlations[typeA][typeB] > corrMax) corrMax = correlations[typeA][typeB];
            if(correlations[typeA][typeB] < corrMin) corrMin = correlations[typeA][typeB];
          }
        }

        // generate labels down here, because otherwise I run into some weird character pointer issue

        char* labelsX[methodTypes.size()];
        char* labelsY[methodTypes.size()];
        for (uint type = 0; type < methodTypes.size(); ++type){
          labelsX[type] = Form("%s : %s", methodTypes.at(type).c_str(), parameter_strings.at(firstPar).c_str());
          labelsY[type] = Form("%s : %s", methodTypes.at(type).c_str(), parameter_strings.at(secondPar).c_str());
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
        corrMatrixHist->GetXaxis()->SetLabelOffset(0.01);
        corrMatrixHist->GetZaxis()->SetTitleOffset(1.7);
        corrMatrixHist->GetYaxis()->SetLabelSize(0.04);
        corrMatrixHist->SetTitle("");

        corrMatrixHist->SetName(Form("CorrelationMatrix_%s_%s", parameter_strings.at(firstPar).c_str(), parameter_strings.at(secondPar).c_str()));
        // corrMatrixHist->SetTitle(Form("CorrelationMatrix_%s_%s", parameter_strings.at(firstPar).c_str(), parameter_strings.at(secondPar).c_str()));


        TCanvas* c_corrMat = new TCanvas("c_corrMat","Correlation Matrix",50,10,1100,900);

        gStyle->SetPaintTextFormat("1.4f");

        c_corrMat->SetGrid();
        c_corrMat->SetTopMargin(0.1);
        c_corrMat->SetBottomMargin(0.12);
        c_corrMat->SetRightMargin(0.18);
        c_corrMat->SetLeftMargin(0.18);

        corrMatrixHist->Draw("COLZTEXT");

        c_corrMat->SetName(Form("CorrelationMatrix_Plot_%s_%s", parameter_strings.at(firstPar).c_str(), parameter_strings.at(secondPar).c_str()));
        c_corrMat->SetTitle(Form("CorrelationMatrix_Plot_%s_%s", parameter_strings.at(firstPar).c_str(), parameter_strings.at(secondPar).c_str()));

        corrMatrixHist->GetZaxis()->SetTitle("Correlation");
        gPad->Modified();
        gPad->Update();

        if(parameter_strings.at(firstPar).compare("R") == 0 && parameter_strings.at(secondPar).compare("R") == 0) 
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


/////////////////////////////////////////////////////////////////////////////////////
*/
        delete c_corrMat;
        delete corrMatrixHist;

    } // end loop over parameter 2
  } // end loop over parameter 1

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  auto methodAverage_dir = outputFile->mkdir("MethodAverage");
  auto methodAverage_graphs_dir = methodAverage_dir->mkdir("Graphs");
  auto methodAverage_histograms_dir = methodAverage_dir->mkdir("Histograms");
  auto methodAverage_canvases_dir = methodAverage_dir->mkdir("Canvases");


 // construct correlation matrices for the different methods, averaged over the different analyzers
  // only do it for the R parameter

  // vector<string> parameter_strings = {"N", "#tau", "A", "R", "#phi"};
  vector<string> methodTypes_avg = {"RE_T", "RE_A", "RW_T", "RW_A", "RW_R", "Q"}; // make sure the pointers down below correspond to the correct labels
  // vector<string> methodTypes_avg = {"RE_T", "RE_A", "RE_TA", "RW_T", "RW_A", "RW_R", "RW_TA", "RW_TR", "RW_AR", "RW_TAR", "Q"}; // make sure the pointers down below correspond to the correct labels





  for (uint firstPar = 0; firstPar < parameter_strings.size(); ++firstPar)
  {
    for (uint secondPar = firstPar; secondPar < parameter_strings.size(); ++secondPar)
    {
        vector<double*> typePointersA = {david_T_pars[firstPar], david_A_pars[firstPar], RW_T_pars[firstPar], RW_A_pars[firstPar], nick_R_pars[firstPar], tim_Q_pars[firstPar]};
        vector<double*> typePointersB = {david_T_pars[secondPar], david_A_pars[secondPar], RW_T_pars[secondPar], RW_A_pars[secondPar], nick_R_pars[secondPar], tim_Q_pars[secondPar]};

        // vector<double*> typePointersA = {david_T_pars[firstPar], david_A_pars[firstPar], RE_TA_pars[firstPar], RW_T_pars[firstPar], RW_A_pars[firstPar], nick_R_pars[firstPar], RW_TA_pars[firstPar], RW_TR_pars[firstPar], RW_AR_pars[firstPar], RW_TAR_pars[firstPar], tim_Q_pars[firstPar]};
        // vector<double*> typePointersB = {david_T_pars[secondPar], david_A_pars[secondPar], RE_TA_pars[secondPar], RW_T_pars[secondPar], RW_A_pars[secondPar], nick_R_pars[secondPar], RW_TA_pars[secondPar], RW_TR_pars[secondPar], RW_AR_pars[secondPar], RW_TAR_pars[secondPar], tim_Q_pars[secondPar]};


        double correlations[methodTypes_avg.size()][methodTypes_avg.size()];

        double corrMax = -1, corrMin = 1;

        for (uint typeA = 0; typeA < methodTypes_avg.size(); ++typeA)
        {
          for (uint typeB = 0; typeB < methodTypes_avg.size(); ++typeB)
          {
            correlations[typeA][typeB] = 0;

              TGraph* corrGraph = new TGraph(totalSamples, typePointersA.at(typeA), typePointersB.at(typeB));
              correlations[typeA][typeB] = corrGraph->GetCorrelationFactor();

              methodAverage_graphs_dir->cd();

              corrGraph->SetTitle(Form(";%s_%s;%s_%s", methodTypes_avg.at(typeA).c_str(), parameter_strings.at(firstPar).c_str(), methodTypes_avg.at(typeB).c_str(), parameter_strings.at(secondPar).c_str()));
              corrGraph->Write(Form("Graph_%s_%s_vs_%s_%s", methodTypes_avg.at(typeA).c_str(), parameter_strings.at(firstPar).c_str(), methodTypes_avg.at(typeB).c_str(), parameter_strings.at(secondPar).c_str()));

              // TCanvas* c_corrMattest = new TCanvas("c_corrMattest","Correlation Matrix",50,10,925,900);
              // corrGraph->Draw();
              // c_corrMattest->SaveAs();

              delete corrGraph;

            if(correlations[typeA][typeB] > corrMax) corrMax = correlations[typeA][typeB];
            if(correlations[typeA][typeB] < corrMin) corrMin = correlations[typeA][typeB];
          }
        }

        // generate labels down here, because otherwise I run into some weird character pointer issue

        char* labelsX[methodTypes_avg.size()];
        char* labelsY[methodTypes_avg.size()];
        for (uint type = 0; type < methodTypes_avg.size(); ++type){
          labelsX[type] = Form("%s : %s", methodTypes_avg.at(type).c_str(), parameter_strings.at(firstPar).c_str());
          labelsY[type] = Form("%s : %s", methodTypes_avg.at(type).c_str(), parameter_strings.at(secondPar).c_str());
        } 

      /////////////////////////////////////////////////////////////////////////////////////

        // Draw correlation matrix

        TH2D* corrMatrixHist = new TH2D("corrMatrixHist","Correlation Matrix", methodTypes_avg.size(), 0, methodTypes_avg.size(), methodTypes_avg.size(), 0, methodTypes_avg.size());

        for(uint i = 0; i < methodTypes_avg.size(); i++){
          for(uint j = 0; j < methodTypes_avg.size(); j++){
            corrMatrixHist->Fill(labelsX[i], labelsY[j], correlations[i][j]);
          }
        }

        corrMatrixHist->LabelsDeflate();
        corrMatrixHist->GetZaxis()->SetRangeUser(corrMin,corrMax);
        corrMatrixHist->GetXaxis()->SetLabelSize(0.04);
        corrMatrixHist->GetXaxis()->SetLabelOffset(0.01);
        corrMatrixHist->GetZaxis()->SetTitleOffset(1.7);
        corrMatrixHist->GetYaxis()->SetLabelSize(0.04);
        corrMatrixHist->SetTitle("");

        corrMatrixHist->SetName(Form("CorrelationMatrix_%s_%s", parameter_strings.at(firstPar).c_str(), parameter_strings.at(secondPar).c_str()));
        // corrMatrixHist->SetTitle(Form("CorrelationMatrix_%s_%s", parameter_strings.at(firstPar).c_str(), parameter_strings.at(secondPar).c_str()));


        TCanvas* c_corrMat = new TCanvas("c_corrMat","Correlation Matrix",50,10,1100,900);

        gStyle->SetPaintTextFormat("1.4f");

        c_corrMat->SetGrid();
        c_corrMat->SetTopMargin(0.1);
        c_corrMat->SetBottomMargin(0.12);
        c_corrMat->SetRightMargin(0.18);
        c_corrMat->SetLeftMargin(0.18);

        corrMatrixHist->Draw("COLZTEXT");

        c_corrMat->SetName(Form("CorrelationMatrix_Plot_%s_%s", parameter_strings.at(firstPar).c_str(), parameter_strings.at(secondPar).c_str()));
        c_corrMat->SetTitle(Form("CorrelationMatrix_Plot_%s_%s", parameter_strings.at(firstPar).c_str(), parameter_strings.at(secondPar).c_str()));

        corrMatrixHist->GetZaxis()->SetTitle("Correlation");
        gPad->Modified();
        gPad->Update();

        if(parameter_strings.at(firstPar).compare("R") == 0 && parameter_strings.at(secondPar).compare("R") == 0) 
          c_corrMat->SaveAs(Form("MethodType_Average_CorrelationMatrixPlot_%s_%s.png", parameter_strings.at(firstPar).c_str(), parameter_strings.at(secondPar).c_str()));

        methodAverage_canvases_dir->cd();
        c_corrMat->Write();

        methodAverage_histograms_dir->cd();
        corrMatrixHist->Write();

        delete c_corrMat;
        delete corrMatrixHist;

    } // end loop over parameter 2
  } // end loop over parameter 1


/////////////////////////////////////////////////////////////////////////////////////

return 1;

}
