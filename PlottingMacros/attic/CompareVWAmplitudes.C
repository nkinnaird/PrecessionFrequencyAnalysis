#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TROOT.h>
#include <TImage.h>
#include <sstream>

using namespace std;


int CompareVWAmplitudes()
{

  TFile *inputFile = TFile::Open("/gm2/data/users/nkinnaird/Ratio/9d-FinalProduction/RandomSeeds/tests_with_VW/fixed-vw-w-tau/generalPlotsVsIter.root");
  // TFile *inputFile = TFile::Open("/gm2/data/users/nkinnaird/Ratio/60h-FinalProduction/RandSeeds/FitIterations/generalPlotsVsIter.root");


	TGraph* Tmethod_VW_amps = (TGraph*) inputFile->Get("topDir/Added/TMethod/Graphs/TMethod_A_VW_Vs_Iter");
	TGraph* Rmethod_VW_amps = (TGraph*) inputFile->Get("topDir/Added/RatioCBO/Graphs/RatioCBO_A_VW_Vs_Iter");

	double seedNumber, T_pointAmp, R_pointAmp, averageAmpRatio = 0;

	for (int pointNo = 0; pointNo < Rmethod_VW_amps->GetN(); ++pointNo)
	{
		Tmethod_VW_amps->GetPoint(pointNo, seedNumber, T_pointAmp);
		Rmethod_VW_amps->GetPoint(pointNo, seedNumber, R_pointAmp);
		cout << "Seed: " << seedNumber << " R amp: " << R_pointAmp << " T amp: " << T_pointAmp << " ratio: " << R_pointAmp/T_pointAmp << endl;
		averageAmpRatio += R_pointAmp/T_pointAmp;
	}

	cout << endl << "Average amplitude ratio: " << averageAmpRatio/Rmethod_VW_amps->GetN() << endl;


return 1;
}
