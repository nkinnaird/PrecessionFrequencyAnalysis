// 3-31-20: Macro to calculate combined correlated error from several datasets. 
// The Run 1 dataset numbers are hardcoded in and can be changed to whatever.
// See this paper for the technique I used to calculate the correlated error: https://arxiv.org/pdf/1507.08210.pdf

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TLegend.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TImage.h>
#include <TLine.h>

#include <Eigen/Dense>

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;

bool saveImages = true;


void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}


int calculateCorrelatedErrors()
{
	double correlationCoefficient = 1;

	// means are in ppm

	double mean_60h = -20.5562; // different software blinding string
	double mean_HighKick = -17.4755;
	double mean_9d = -17.7182;
	double mean_Endgame = -17.3406;

	vector<double> datasetMeans = {mean_60h, mean_HighKick, mean_9d, mean_Endgame};
	// vector<double> datasetMeans = {mean_HighKick, mean_9d, mean_Endgame};

	// errors are statistical and then systematic
	// units in ppb

	pair<double, double> errors_60h = {1358.1, 89.8};
	pair<double, double> errors_HighKick = {1411.2, 80.8};
	pair<double, double> errors_9d = {903.3, 72.9};
	pair<double, double> errors_Endgame = {639.3, 102.8};

	vector<pair<double, double> > datasetErrors = {errors_60h, errors_HighKick, errors_9d, errors_Endgame};
	// vector<pair<double, double> > datasetErrors = {errors_HighKick, errors_9d, errors_Endgame};

/////////////////////////////////////////////////////////////////////////////////////

	// basic checks for comparison with matrix algebra

	double totalStatisticalError = 0;
	double totalErrorUncorrelated = 0;

	for (uint i = 0; i < datasetErrors.size(); ++i)
	{
		totalStatisticalError += 1/(pow(datasetErrors.at(i).first,2));
		totalErrorUncorrelated += 1/(pow(datasetErrors.at(i).first,2) + pow(datasetErrors.at(i).second,2));
	}

	totalStatisticalError = sqrt(1/totalStatisticalError);
	totalErrorUncorrelated = sqrt(1/totalErrorUncorrelated);

	cout << "total statistical error: " << totalStatisticalError << endl;
	cout << "total uncorrelated error: " << totalErrorUncorrelated << endl;

/////////////////////////////////////////////////////////////////////////////////////

  	Eigen::MatrixXd fullErrorMatrix = Eigen::MatrixXd::Zero(datasetMeans.size(),datasetMeans.size());


  	for (uint i = 0; i < datasetErrors.size(); ++i)
  	{
  		for (uint j = 0; j < datasetErrors.size(); ++j)
  		{
  			if(i == j) fullErrorMatrix(i,j) = pow(datasetErrors.at(i).first,2) + pow(datasetErrors.at(i).second,2);
  			else fullErrorMatrix(i,j) = correlationCoefficient * datasetErrors.at(i).second * datasetErrors.at(j).second;
  		}
  	}

  	cout << endl << "fullErrorMatrix: " << endl << fullErrorMatrix << endl;


  	vector<double> datasetWeightFactors;

  	for (uint i = 0; i < datasetErrors.size(); ++i)
  	{
  		Eigen::MatrixXd weightMatrix = fullErrorMatrix;

  		Eigen::VectorXd columnVector = weightMatrix.col(i);

  		weightMatrix = weightMatrix.colwise() - columnVector;

  		removeRow(weightMatrix, i);
  		removeColumn(weightMatrix, i);

  		// cout << endl << "weightMatrix " << i << ": " << endl << weightMatrix << endl;

  		double weightFactor = weightMatrix.determinant();
  		datasetWeightFactors.push_back(weightFactor);

  	}

  	double sumOfWeights = 0;
  	for (uint i = 0; i < datasetWeightFactors.size(); ++i) sumOfWeights += datasetWeightFactors.at(i);


  	double totalError = sqrt(fullErrorMatrix.determinant() / sumOfWeights);
  	cout << "total correlated error: " << totalError << endl;


  	double totalStatisticalPart = 0;
  	for (uint i = 0; i < datasetErrors.size(); ++i) totalStatisticalPart += pow(datasetWeightFactors.at(i),2) * pow(datasetErrors.at(i).first,2);
  	totalStatisticalPart = sqrt(totalStatisticalPart);
    totalStatisticalPart /= sumOfWeights;

    cout << "total statistical part of error: " << totalStatisticalPart << endl;

/////////////////////////////////////////////////////////////////////////////////////

  	double weightedMean = 0;

  	for (uint i = 0; i < datasetMeans.size(); ++i) weightedMean += datasetMeans.at(i) * datasetWeightFactors.at(i);
  	weightedMean /= sumOfWeights;

  	cout << "weighted mean: " << weightedMean << endl;

/////////////////////////////////////////////////////////////////////////////////////

  	// make plot of dataset means against each other

  gStyle->SetOptStat(000000);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1);
  gStyle->SetLineColor(1);
  gStyle->SetPadRightMargin(.05);
  gStyle->SetPadLeftMargin(.15);


  	TGraphErrors* datasetGraph = new TGraphErrors();
 	int pointNo = 0;

  	for (uint i = 0; i < datasetMeans.size(); ++i)
  	{
  		datasetGraph->SetPoint(pointNo, i+1, datasetMeans.at(i));
  		datasetGraph->SetPointError(pointNo, 0, sqrt(pow(datasetErrors.at(i).first,2) + pow(datasetErrors.at(i).second,2))/1000);

  		pointNo++;
  	}

   datasetGraph->GetHistogram()->GetXaxis()->Set(datasetMeans.size(),0.5,datasetMeans.size()+0.5);

   for(uint i=0; i < datasetMeans.size(); i++){
      datasetGraph->GetHistogram()->GetXaxis()->SetBinLabel(i+1, dataset_names.at(i).c_str());
   }

   datasetGraph->GetYaxis()->SetTitle("R (ppm, blinded)");
   datasetGraph->GetXaxis()->SetLabelSize(0.08);

    auto canv = new TCanvas("canv","canv",200,10,800,600);
    datasetGraph->Draw("AP");

    canv->Update();

  TLine *line = new TLine(1.5, canv->GetUymin(), 1.5, canv->GetUymax());
  line->SetLineColor(4);
  line->SetLineStyle(2);
  line->SetLineWidth(3);
  line->Draw("SAME");

  if(saveImages) canv->SaveAs("R_dataset_comparison.png");


	return 1;
}
