// 3-31-20: Macro which printed out fit results in a format which was at one time the input to a combination analysis I believe by Alex.

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2D.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TImage.h>
#include <sstream>
#include <TNtuple.h>
#include <TPaveStats.h>
#include <TMatrixD.h>
#include <TDatime.h>
#include <TVectorD.h>

using namespace std;


int OutputForCombination(std::string filePath)
{

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }


string BU = "BU";
string dataset = "60hour";
string analysisMethod;

TDatime* newDate = new TDatime();
string date = to_string(newDate->GetMonth()) + "-" + to_string(newDate->GetDay()) + "-" + to_string(newDate->GetYear());

double usedBinWidth = ((TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/Input/Times"))->GetBinWidth(1);

/////////////////////////////////////////////////////////////////////////////////////

if(inputFile->Get("topDir/FitPasses/FitPass0/addedDir/FullRatio/Added_Times_Full_Ratio_Graph") != 0)
{
    TF1* fullRatio_fit = (TF1*) ((TGraphErrors*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/FullRatio/Added_Times_Full_Ratio_Graph"))->GetFunction("fullRatioFitFunc");
	const int nPars = fullRatio_fit->GetNpar();

      analysisMethod = "R";
      string fileName = dataset + "_" + BU + "_" + analysisMethod + "_" + date + ".txt";

	  ofstream ratioFile;
	  ratioFile.open (fileName.c_str());

	  ratioFile << dataset << "\n";
	  ratioFile << analysisMethod << "\n";
	  ratioFile << BU << "\n";
	  ratioFile << nPars << "\n";

/*
	TVectorD* histSavedParameters = (TVectorD*) inputFile->Get("topDir/FitPasses/FitPass0/parameterStore");
    float energyThreshold = (*histSavedParameters)[1];
    
      ratioFile << energyThreshold << "\n";
*/

	  ratioFile << usedBinWidth << "\n";

      ratioFile << fullRatio_fit->GetXmin() << "\n";
      ratioFile << fullRatio_fit->GetXmax() << "\n";

      ratioFile << "pluscosplus" << "\n";

	for (int parNum = 0; parNum < nPars; ++parNum){
		string paramString = fullRatio_fit->GetParName(parNum);
		ratioFile << paramString << " " << fullRatio_fit->GetParameter(parNum) << " " << fullRatio_fit->GetParError(parNum);

		    if(paramString.find("tau") != string::npos) ratioFile << " \u03BCs";
		    else if(paramString.find("omega") != string::npos) ratioFile << " rad/\u03BCs";

		ratioFile << "\n";
	} 

	TMatrixD* corrMatrix = (TMatrixD*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/FullRatio/FullRatioFitCorrMatrix");

	for(int i = 0; i < nPars; i++){
	  for(int j = 0; j < nPars; j++){
	    ratioFile << (*corrMatrix)[i][j] << "\n";
	  }
	}


	ratioFile.close();
}

/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

if(inputFile->Get("topDir/FitPasses/FitPass0/addedDir/TMethod/allTimesAdded_TMethod") != 0)
{
	TF1* TMethod_fit = (TF1*) ((TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/TMethod/allTimesAdded_TMethod"))->GetFunction("TmethodFitFunc");
	const int nPars = TMethod_fit->GetNpar();

      analysisMethod = "T";
      string fileName = dataset + "_" + BU + "_" + analysisMethod + "_" + date + ".txt";

	  ofstream TmethodFile;
	  TmethodFile.open (fileName.c_str());

	  TmethodFile << dataset << "\n";
	  TmethodFile << analysisMethod << "\n";
	  TmethodFile << BU << "\n";
	  TmethodFile << nPars << "\n";

/*
	TVectorD* histSavedParameters = (TVectorD*) inputFile->Get("topDir/FitPasses/FitPass0/parameterStore");
    float energyThreshold = (*histSavedParameters)[1];
    
      TmethodFile << energyThreshold << "\n";
*/

	  TmethodFile << usedBinWidth << "\n";

      TmethodFile << TMethod_fit->GetXmin() << "\n";
      TmethodFile << TMethod_fit->GetXmax() << "\n";

      TmethodFile << "pluscosplus" << "\n";

	for (int parNum = 0; parNum < nPars; ++parNum){
		string paramString = TMethod_fit->GetParName(parNum);
		TmethodFile << paramString << " " << TMethod_fit->GetParameter(parNum) << " " << TMethod_fit->GetParError(parNum);

		    if(paramString.find("tau") != string::npos) TmethodFile << " \u03BCs";
		    else if(paramString.find("omega") != string::npos) TmethodFile << " rad/\u03BCs";

		TmethodFile << "\n";
	} 

	TMatrixD* corrMatrix = (TMatrixD*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/TMethod/CorrelationMatrix");

	for(int i = 0; i < nPars; i++){
	  for(int j = 0; j < nPars; j++){
	    TmethodFile << (*corrMatrix)[i][j] << "\n";
	  }
	}


	TmethodFile.close();
}

/////////////////////////////////////////////////////////////////////////////////////

return 1;

}
