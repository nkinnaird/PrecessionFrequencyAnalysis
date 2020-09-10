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

vector<vector<double>> Import(const string &p_file) {
  // the input data file
  ifstream infile;
  infile.open(p_file, ios::in);

  // the input data
  vector<vector<double>> values;

  // the line of input data
  string line;

  // discard the header
  getline(infile, line);

  // read in the input data
  while (getline(infile, line)) {
    // the variables to parse the line of input data
    stringstream stream(line);
    string entry;
    vector<double> entries;

    // discard the header
    getline(stream, entry, ',');

    // parse the line of input data
    while (getline(stream, entry, ',')) {
      // the value
      double value = (entry == "" || entry == "\r") ? 0.0 : stod(entry);

      // include the value
      entries.push_back(value);
    }

    // include the data
    values.push_back(move(entries));
  }

  // close the input data
  infile.close();

  // return the input data
  return values;
} // end function Import


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


int combineRun1Errors()
{

  const auto errors_60h = Import("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/PlottingMacros/DatasetComparisonPlotters/Run1Combination/Data/1a.csv");

uint numAnalysisTypes = errors_60h.at(0).size();
uint numIndividualSystematics = errors_60h.size();

vector<double> errors_60h_systematic; // vector of analysis type systematics

  // for (uint j = 0; j < numAnalysisTypes; ++j) // columns
  // {
  //   double columnError = 0;

  //   for (uint i = 0; i < numIndividualSystematics; ++i) // rows
  //   {
  //     columnError +=  pow(errors_60h.at(i).at(j),2); // this is assuming all the systematics in one analysis are uncorrelated with each other
  //   }
  //   errors_60h_systematic.push_back(sqrt(columnError));
  // }

  // for (uint i = 0; i < numAnalysisTypes; ++i)
  // {
  //   cout << "Analysis " << i << " total systematic error: " << errors_60h_systematic.at(i) << endl;
  // }


errors_60h_systematic.push_back(100);
errors_60h_systematic.push_back(200);

  vector<double> errors_60h_statistical = {200, 200};
  // vector<double> errors_60h_statistical = {1358.17, 1359.81}; // my T and R methods
  // vector<double> errors_60h_statistical = {0, 0}; // set to 0 to just look at the systematics part
  // vector<double> errors_60h_statistical;
  // for (uint i = 0; i < numAnalysisTypes; ++i) errors_60h_statistical.push_back(0);

  uint numberOfAnalysesToCombine = errors_60h_statistical.size();

	double corr_anaType_statistical = 0; // will want to pull these in with another csv file
  double corr_anaType_systematic = 1;


/////////////////////////////////////////////////////////////////////////////////////

  	Eigen::MatrixXd fullErrorMatrix = Eigen::MatrixXd::Zero(numberOfAnalysesToCombine,numberOfAnalysesToCombine);


  	for (uint i = 0; i < numberOfAnalysesToCombine; ++i)
  	{
  		for (uint j = 0; j < numberOfAnalysesToCombine; ++j)
  		{
  			if(i == j) fullErrorMatrix(i,j) = pow(errors_60h_statistical.at(i),2) + pow(errors_60h_systematic.at(i),2);
  			else fullErrorMatrix(i,j) = corr_anaType_statistical * errors_60h_statistical.at(i) * errors_60h_statistical.at(j) + corr_anaType_systematic * errors_60h_systematic.at(i) * errors_60h_systematic.at(j);
  		}
  	}

  	cout << endl << "fullErrorMatrix: " << endl << endl << fullErrorMatrix << endl << endl;


  	vector<double> analysisWeightFactors;

  	for (uint i = 0; i < numberOfAnalysesToCombine; ++i)
  	{
  		Eigen::MatrixXd weightMatrix = fullErrorMatrix;

  		Eigen::VectorXd columnVector = weightMatrix.col(i);

  		weightMatrix = weightMatrix.colwise() - columnVector;

  		removeRow(weightMatrix, i);
  		removeColumn(weightMatrix, i);

  		// cout << endl << "weightMatrix " << i << ": " << endl << weightMatrix << endl;

  		double weightFactor = weightMatrix.determinant();

      cout << "weight factor " << i << ": " << weightFactor << endl;

  		analysisWeightFactors.push_back(weightFactor);

  	}

    cout << endl;


  	double sumOfWeights = 0;
  	for (uint i = 0; i < analysisWeightFactors.size(); ++i) sumOfWeights += analysisWeightFactors.at(i);
    cout << "Sum of weights: " << sumOfWeights << endl;

    for (uint i = 0; i < analysisWeightFactors.size(); ++i){
      cout << "Normalized weights: " << i << " " << analysisWeightFactors.at(i)/sumOfWeights << endl;
    } 




    double errorMatrixDeterminant = fullErrorMatrix.determinant();
    cout << "Error matrix determinant: " << errorMatrixDeterminant << endl;

    if(errorMatrixDeterminant / sumOfWeights < 0) cout << "Variance was negative." << endl;

  	double totalError = sqrt(abs(errorMatrixDeterminant / sumOfWeights));
  	cout << "Total error: " << totalError << endl;


  	// double totalStatisticalPart = 0;
  	// for (uint i = 0; i < numberOfAnalysesToCombine; ++i) totalStatisticalPart += pow(analysisWeightFactors.at(i),2) * pow(datasetErrors.at(i).first,2);
  	// totalStatisticalPart = sqrt(totalStatisticalPart);
   //  totalStatisticalPart /= sumOfWeights;

   //  cout << "total statistical part of error: " << totalStatisticalPart << endl;


/////////////////////////////////////////////////////////////////////////////////////

	return 1;
}
