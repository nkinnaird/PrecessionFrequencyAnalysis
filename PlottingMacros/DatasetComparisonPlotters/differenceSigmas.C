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

using namespace std;

double allowedDeviation(double errA, double errB, double corr)
{
	return sqrt(errA*errA + errB*errB - 2 * corr * errA*errB);
}


int differenceSigmas()
{
	double corr = 0.9983;

	std::vector<double> Ts;
	std::vector<double> Rs;
	std::vector<double> Errs;


	// Ts.push_back(5);
	// Ts.push_back(6);

	// Rs.push_back(5.05);
	// Rs.push_back(6.05);

	// Errs.push_back(1);
	// Errs.push_back(0.55);


	Ts.push_back(-20.3016);
	Ts.push_back(-16.6644);
	Ts.push_back(-17.5373);
	Ts.push_back(-17.3222);

	Rs.push_back(-20.4661);
	Rs.push_back(-16.8295);
	Rs.push_back(-17.542);
	Rs.push_back(-17.3856);

	Errs.push_back(1.3582);
	Errs.push_back(1.1561);
	Errs.push_back(0.9301);
	Errs.push_back(0.7584);

/////////////////////////////////////////////////////////////////////////////////////

	std::vector<double> RTdiffs;
	std::vector<double> allowedDiffs;

	for (uint i = 0; i < Ts.size(); ++i)
	{
		RTdiffs.push_back(Rs.at(i) - Ts.at(i));
		allowedDiffs.push_back(allowedDeviation(Errs.at(i), Errs.at(i), corr));
	}


	for (uint i = 0; i < Ts.size(); ++i)
	{
		cout << "i: " << i << " R-T diff: " << RTdiffs.at(i) << " 1sigma allowed: " << allowedDiffs.at(i) << " number deviations: " << RTdiffs.at(i)/allowedDiffs.at(i) << endl;
	}


	double sumInvW2 = 0;

	for (uint i = 0; i < Ts.size(); ++i)
	{
		sumInvW2 += 1/pow(Errs.at(i),2);
	}

	double variance = 1/sumInvW2;
	double totalError = sqrt(variance);

	cout << "Total error: " << totalError << endl;

	double totalAllowedDiff = allowedDeviation(totalError, totalError, corr);

	cout << "Total allowed 1sigma: " << totalAllowedDiff << endl;


	double totalT = 0, totalR = 0;

	for (uint i = 0; i < Ts.size(); ++i)
	{
		totalT += Ts.at(i)/pow(Errs.at(i),2);
		totalR += Rs.at(i)/pow(Errs.at(i),2);
	}

	totalT *= variance;
	totalR *= variance;


	cout << "Combined T: " << totalT << " combined R: " << totalR << endl;

	double totalRTDiff = totalR - totalT;

	cout << "Total RT diff: " << totalRTDiff << " number deviations: " << totalRTDiff/totalAllowedDiff << endl;


	cout << "Careful about the blinding, can't calculate the total T and total R necessarily and then get the R-T difference, but I can get the combined R-T difference from the individual ones." << endl;


	double combinedRTDiffFromSum = 0;

	for (uint i = 0; i < Ts.size(); ++i)
	{
		combinedRTDiffFromSum += RTdiffs.at(i)/pow(Errs.at(i),2);
	}

	combinedRTDiffFromSum *= variance;

	cout << "Avoiding the blinding issue: total RT diff: " << combinedRTDiffFromSum << " number deviations: " << combinedRTDiffFromSum/totalAllowedDiff << endl;


	cout << "The numbers end up being exactly the same - so maybe for determining the number of deviations it doesn't care about any constant offset from the blinding." << endl;



	return 1;
}
