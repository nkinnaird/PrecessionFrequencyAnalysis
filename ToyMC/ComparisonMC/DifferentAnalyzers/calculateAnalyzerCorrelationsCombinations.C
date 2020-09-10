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

#include <algorithm>
#include <Eigen/Dense>


// #include "plotUtils.hh"

using namespace std;

uint totalSamples = 0;

bool clipCoefficients = false;
bool matrixPrintout = false;

/////////////////////////////////////////////////////////////////////////////////////

// make sure input parameters and correlations have the right order
double getWeightedAverage(string combo, int sampleNum, vector<double**> pars, TH2D* corrMat, vector<double> corrIndices){

  if(sampleNum != 0) matrixPrintout = false;

  if(matrixPrintout) cout << endl << "Combination: " << combo << endl;

  if(pars.size() != corrIndices.size()){
    cout << "Not passing in the right number of methods or correlations." << endl;
    exit(0);
  }  

  vector<double> Rs (pars.size(), 0);
  vector<double> Rerrs (pars.size(), 0);

  for (uint methodNum = 0; methodNum < pars.size(); ++methodNum)
  {
    Rs.at(methodNum) = pars.at(methodNum)[0][sampleNum];
    Rerrs.at(methodNum) = pars.at(methodNum)[1][sampleNum];

    if(matrixPrintout) cout << "Method: " << methodNum << " R: " << Rs.at(methodNum) << " err: " << Rerrs.at(methodNum) << endl;
  }

  Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(pars.size(),pars.size());

  for (uint i = 0; i < pars.size(); ++i)
  {
    for (uint j = 0; j < pars.size(); ++j)
    {
      double corr = corrMat->GetBinContent(corrIndices.at(i), corrIndices.at(j));

      if(clipCoefficients){
        double clipCorr = min(Rerrs.at(i), Rerrs.at(j))/max(Rerrs.at(i), Rerrs.at(j));
        if(matrixPrintout) cout << "Checking clipping. corr: " << corr << " clipCorr: " << clipCorr << " diff: " << corr-clipCorr << endl;
        if(corr > clipCorr) corr = clipCorr;
      }

      cov(i,j) = corr * Rerrs.at(i) * Rerrs.at(j);
      if(matrixPrintout) cout << "Matrix element i,j = (" << i << "," << j << ")= " << corr << " * " << Rerrs.at(i) << " * " << Rerrs.at(j) << endl;
    }
  }

  if(matrixPrintout) cout << endl << "Cov: " << endl << cov << endl;

  Eigen::MatrixXd covInverse = cov.inverse();

  if(matrixPrintout) cout << endl << "Cov inverse: " << endl << covInverse << endl << endl;

  double weightedAverage = 0;

  for (uint i = 0; i < pars.size(); ++i)
  {
    double weight = covInverse.row(i).sum() / covInverse.sum();
    if(matrixPrintout) cout << "weight i : " << i << " row sum: " << covInverse.row(i).sum() << " mat sum: " << covInverse.sum() << " weight: " << weight << endl;

    weightedAverage += weight * Rs.at(i);
  }

  if(matrixPrintout) cout << endl << "weightedAverage: " << weightedAverage << endl;

  return weightedAverage;
}


/////////////////////////////////////////////////////////////////////////////////////

void readInTree(TNtuple* inTree, double** parameters){

  float R, R_err;
  inTree->SetBranchAddress("R", &R);
  inTree->SetBranchAddress("R_err", &R_err);


  for (uint sampleNum = 0; sampleNum < totalSamples; ++sampleNum)
  {
    inTree->GetEntry(sampleNum);

    parameters[0][sampleNum] = R;
    parameters[1][sampleNum] = R_err;
  } // end loop over samples
}


int calculateAnalyzerCorrelationsCombinations(std::string filePath)
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

  if(clipCoefficients) cout << "Clipping coefficients." << endl;
  else cout << "Not clipping coefficients." << endl;

/////////////////////////////////////////////////////////////////////////////////////

  TFile* outputFile = new TFile("analyzerCorrelations.root","RECREATE");

  auto graphsDir = outputFile->mkdir("Graphs");
  auto histogramsDir = outputFile->mkdir("Histograms");
  auto canvasesDir = outputFile->mkdir("Canvases");
  auto deviationsDir = outputFile->mkdir("Deviations");

/////////////////////////////////////////////////////////////////////////////////////

  // dynamically allocate memory so that the stack isn't overflown - https://www.techiedelight.com/dynamic-memory-allocation-in-c-for-2d-3d-array/

  double** nick_T_pars = new double*[2]; // 2 different parameters, 1 for R and 1 for the error on R
  double** nick_R_pars = new double*[2]; 
  double** david_T_pars = new double*[2]; 
  double** david_A_pars = new double*[2]; 
  double** aaron_T_pars = new double*[2]; 
  double** aaron_A_pars = new double*[2]; 
  double** matteo_T_pars = new double*[2]; 
  double** matteo_A_pars = new double*[2]; 
  double** bingzhi_T_pars = new double*[2]; 
  double** bingzhi_A_pars = new double*[2]; 
  double** tim_Q_pars = new double*[2]; 

  // averaged and combined parameters

  double** RW_T_pars = new double*[2]; 
  double** RW_A_pars = new double*[2];

  double ** RE_TA_pars = new double*[2];

  double ** RW_TA_pars = new double*[2];
  double ** RW_TR_pars = new double*[2];
  double ** RW_AR_pars = new double*[2];
  double ** RW_TAR_pars = new double*[2];


  for (int i = 0; i < 2; ++i)
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

    RE_TA_pars[i] = new double[totalSamples];

    RW_TA_pars[i] = new double[totalSamples];
    RW_TR_pars[i] = new double[totalSamples];
    RW_AR_pars[i] = new double[totalSamples];
    RW_TAR_pars[i] = new double[totalSamples];
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
        
        float R_nick_R, R_nick_R_err;
        nick_R_fitParams->SetBranchAddress("R", &R_nick_R);
        nick_R_fitParams->SetBranchAddress("R_err", &R_nick_R);

        for (uint sampleNum = 0; sampleNum < totalSamples; ++sampleNum)
        {
          nick_R_fitParams->GetEntry(sampleNum);

          nick_R_pars[0][sampleNum] = R_nick_R;
          // nick_R_pars[1][sampleNum] = R_nick_R_err; // bug in my code where I was storing the R value also into the R_error tree entry, instead just use my T method R error
          nick_R_pars[1][sampleNum] = nick_T_pars[1][sampleNum];

        } // end iter loop

/////////////////////////////////////////////////////////////////////////////////////

 // read in previously made file of correlation coefficients so I can use them to get the weights

  TFile* prevFile = TFile::Open("../../analyzerCorrelations.root");
   if (prevFile == 0) {
      cout << "cannot open previously made correlations file" << endl;
      return 0;
   }

  TH2D* previousReconCorrelations = (TH2D*) prevFile->Get("MethodAverage/Histograms/CorrelationMatrix_R_R");

 // make average/combined parameters

  for (uint sampleNum = 0; sampleNum < totalSamples; ++sampleNum)
  {
    for (int parNum = 0; parNum < 2; ++parNum)
    {
      RW_T_pars[parNum][sampleNum] = (nick_T_pars[parNum][sampleNum] + aaron_T_pars[parNum][sampleNum] + matteo_T_pars[parNum][sampleNum] + bingzhi_T_pars[parNum][sampleNum])/4.;
      RW_A_pars[parNum][sampleNum] = (aaron_A_pars[parNum][sampleNum] + matteo_A_pars[parNum][sampleNum] + bingzhi_A_pars[parNum][sampleNum])/3.;

/////////////////////////////////////////////////////////////////////////////////////
// simple average

      // RE_TA_pars[parNum][sampleNum] = (david_T_pars[parNum][sampleNum] + david_A_pars[parNum][sampleNum])/2.;

      // RW_TA_pars[parNum][sampleNum] = (nick_T_pars[parNum][sampleNum] + aaron_T_pars[parNum][sampleNum] + matteo_T_pars[parNum][sampleNum] + bingzhi_T_pars[parNum][sampleNum] + aaron_A_pars[parNum][sampleNum] + matteo_A_pars[parNum][sampleNum] + bingzhi_A_pars[parNum][sampleNum])/7.;
      // RW_TR_pars[parNum][sampleNum] = (nick_T_pars[parNum][sampleNum] + aaron_T_pars[parNum][sampleNum] + matteo_T_pars[parNum][sampleNum] + bingzhi_T_pars[parNum][sampleNum] + nick_R_pars[parNum][sampleNum])/5.;
      // RW_AR_pars[parNum][sampleNum] = (aaron_A_pars[parNum][sampleNum] + matteo_A_pars[parNum][sampleNum] + bingzhi_A_pars[parNum][sampleNum] + nick_R_pars[parNum][sampleNum])/4.;
      // RW_TAR_pars[parNum][sampleNum] = (nick_T_pars[parNum][sampleNum] + aaron_T_pars[parNum][sampleNum] + matteo_T_pars[parNum][sampleNum] + bingzhi_T_pars[parNum][sampleNum] + aaron_A_pars[parNum][sampleNum] + matteo_A_pars[parNum][sampleNum] + bingzhi_A_pars[parNum][sampleNum] + nick_R_pars[parNum][sampleNum])/8.;

/////////////////////////////////////////////////////////////////////////////////////

    } // end parNum loop

/////////////////////////////////////////////////////////////////////////////////////    
// weighted average - don't worry about combining the errors in any way - so do the averaging outside the parNum loop

      RE_TA_pars[0][sampleNum] = getWeightedAverage("RE_TA", sampleNum, {david_T_pars, david_A_pars}, previousReconCorrelations, {1,2});

      RW_TA_pars[0][sampleNum] = getWeightedAverage("RW_TA", sampleNum, {RW_T_pars, RW_A_pars}, previousReconCorrelations, {3,4});
      RW_TR_pars[0][sampleNum] = getWeightedAverage("RW_TR", sampleNum, {RW_T_pars, nick_R_pars}, previousReconCorrelations, {3,5});
      RW_AR_pars[0][sampleNum] = getWeightedAverage("RW_AR", sampleNum, {RW_A_pars, nick_R_pars}, previousReconCorrelations, {4,5});

      RW_TAR_pars[0][sampleNum] = getWeightedAverage("RW_TAR", sampleNum, {RW_T_pars, RW_A_pars, nick_R_pars}, previousReconCorrelations, {3,4,5});
  }


// exit(0);

/////////////////////////////////////////////////////////////////////////////////////

 // construct correlation matrices for each individual parameter among the different methods

  vector<string> parameter_strings = {"R"};
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

  // vector<string> methodTypes_avg = {"RE_T", "RE_A", "RW_T", "RW_A", "RW_R", "Q"}; // make sure the pointers down below correspond to the correct labels
  vector<string> methodTypes_avg = {"RE_T", "RE_A", "RE_TA", "RW_T", "RW_A", "RW_R", "RW_TA", "RW_TR", "RW_AR", "RW_TAR", "Q"}; // make sure the pointers down below correspond to the correct labels



  for (uint firstPar = 0; firstPar < parameter_strings.size(); ++firstPar)
  {
    for (uint secondPar = firstPar; secondPar < parameter_strings.size(); ++secondPar)
    {
        // vector<double*> typePointersA = {david_T_pars[firstPar], david_A_pars[firstPar], RW_T_pars[firstPar], RW_A_pars[firstPar], nick_R_pars[firstPar], tim_Q_pars[firstPar]};
        // vector<double*> typePointersB = {david_T_pars[secondPar], david_A_pars[secondPar], RW_T_pars[secondPar], RW_A_pars[secondPar], nick_R_pars[secondPar], tim_Q_pars[secondPar]};

        vector<double*> typePointersA = {david_T_pars[firstPar], david_A_pars[firstPar], RE_TA_pars[firstPar], RW_T_pars[firstPar], RW_A_pars[firstPar], nick_R_pars[firstPar], RW_TA_pars[firstPar], RW_TR_pars[firstPar], RW_AR_pars[firstPar], RW_TAR_pars[firstPar], tim_Q_pars[firstPar]};
        vector<double*> typePointersB = {david_T_pars[secondPar], david_A_pars[secondPar], RE_TA_pars[secondPar], RW_T_pars[secondPar], RW_A_pars[secondPar], nick_R_pars[secondPar], RW_TA_pars[secondPar], RW_TR_pars[secondPar], RW_AR_pars[secondPar], RW_TAR_pars[secondPar], tim_Q_pars[secondPar]};


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
