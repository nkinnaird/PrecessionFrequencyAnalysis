// 3-31-20: Macro to calculate some average and max errors from the IFG parameters. Input files are those given down below, having been taken from various elogs from https://muon.npl.washington.edu/elog/g2/Production/

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>

using namespace std;

void calculateLTDPError(string fileName)
{
	ifstream inputFile(fileName);

	if(!inputFile.is_open()){
		cout << "Trouble opening file: " << fileName << endl;
		exit(1);
	}

  	string line;
  	getline(inputFile, line); // get header line

  	vector<double> taus, tauErrs;

  	 while(!inputFile.eof()){
  		double calo, xtal, tau, tauErr;
		inputFile >> calo >> xtal >> tau >> tauErr;
		if(inputFile.eof()) break;
		taus.push_back(tau);
		tauErrs.push_back(tauErr);
  	 } 


  	 double maxTauSigma = 0;
  	 double averageTauSigma = 0; // errors converted to uncertainties in percent

  	 for (uint i = 0; i < taus.size(); ++i)
  	 {
  	   	double tauUncertainty = tauErrs.at(i) / taus.at(i);
  	 	averageTauSigma += tauUncertainty;
  	 	if(tauUncertainty > maxTauSigma) maxTauSigma = tauUncertainty;	
  	 }

  	 averageTauSigma /= taus.size();

  	 cout << "Average lifetime relative uncertainty is: " << averageTauSigma << " maximum uncertainty: " << maxTauSigma << endl;

  	 inputFile.close();
}


void calculateAverageErrors(string fileName)
{
	ifstream inputFile(fileName);

	if(!inputFile.is_open()){
		cout << "Trouble opening file: " << fileName << endl;
		exit(1);
	}

  	string line;
  	getline(inputFile, line); // get header line

  	vector<double> amps, taus, ampErrs, tauErrs, covATs;

  	 while(!inputFile.eof()){
  		double calo, xtal, amp, tau, ampErr, tauErr, covAA, covTT, covAT;
		inputFile >> calo >> xtal >> amp >> tau >> ampErr >> tauErr >> covAA >> covTT >> covAT;
		// cout << " " << calo << " " << xtal << " " << amp << " " << tau << " " << ampErr << " " << tauErr << " " << covAA << " " << covTT << " " << covAT << endl;
		if(inputFile.eof()) break;
		amps.push_back(amp); 
		taus.push_back(tau);
		ampErrs.push_back(ampErr);
		tauErrs.push_back(tauErr);
		covATs.push_back(covAT);
  	 } 


  	 double maxAmpSigma = 0, maxTauSigma = 0, maxCorrAT = 0;
  	 double averageAmpSigma = 0, averageTauSigma = 0, averageCorrAT = 0; // errors converted to uncertainties in percent

  	 for (uint i = 0; i < amps.size(); ++i)
  	 {
  	 	double ampUncertainty = ampErrs.at(i) / amps.at(i);
  	 	averageAmpSigma += ampUncertainty;
  	 	if(ampUncertainty > maxAmpSigma) maxAmpSigma = ampUncertainty;

  	   	double tauUncertainty = tauErrs.at(i) / taus.at(i);
  	 	averageTauSigma += tauUncertainty;
  	 	if(tauUncertainty > maxTauSigma) maxTauSigma = tauUncertainty;	

	  	double corrAT = covATs.at(i) / (ampErrs.at(i) * tauErrs.at(i));
	  	averageCorrAT += corrAT;
	  	if(abs(corrAT) > abs(maxCorrAT)) maxCorrAT = corrAT;	
  	 }


  	 averageAmpSigma /= amps.size();
  	 averageTauSigma /= taus.size();
  	 averageCorrAT /= covATs.size();

  	 cout << "Average amplitude relative uncertainty is: " << averageAmpSigma << " maximum uncertainty: " << maxAmpSigma << endl;
  	 cout << "Average lifetime relative uncertainty is: " << averageTauSigma << " maximum uncertainty: " << maxTauSigma << endl;
  	 cout << "Average correlation is: " << averageCorrAT << " maximum correlation: " << maxCorrAT << endl;	

  	 inputFile.close();

}



int IFGErrors()
{
	int totalCrystals = 1296;

	string parameters60h = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/Misc/IFGErrors/parameters60h.txt";
	string parametersHighKick = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/Misc/IFGErrors/infill_highkick_withstdp_taultdp_cov.txt";
	string parameters9d = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/Misc/IFGErrors/parameters9d.txt";
	string parametersEndgame = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/Misc/IFGErrors/infill_endgame_withstdp_taultdp_cov.txt";
	
  	cout << endl << "60h IFG Parameter Errors: " << endl;
	calculateAverageErrors(parameters60h);
  	cout << endl;

  	cout << endl << "HighKick IFG Parameter Errors: " << endl;
	calculateAverageErrors(parametersHighKick);
  	cout << endl;

  	cout << endl << "9d IFG Parameter Errors: " << endl;
	calculateAverageErrors(parameters9d);
  	cout << endl;

  	cout << endl << "Endgame IFG Parameter Errors: " << endl;
	calculateAverageErrors(parametersEndgame);
  	cout << endl;


	string tauLTDP_parameters = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/Misc/IFGErrors/tau_ltdp_run1.txt";

  	cout << endl << "Tau LTDP Error: " << endl;
	calculateLTDPError(tauLTDP_parameters);
  	cout << endl;


	return 1;
}
