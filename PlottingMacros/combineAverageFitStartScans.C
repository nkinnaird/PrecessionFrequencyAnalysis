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
#include <TPaveStats.h>

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;


void lombScargle(TGraph* inGraph, double windowTime){
  ////////////////////////////////////////////////////////////////
  // Lomb Scargle (of sorts)
  ////////////////////////////////////////////////////////////////

  // Input params
  double omega = defaultWa * 1e3;
  double window = windowTime; // us (windows to scan in)

  // Output graph
  TGraph* powerVsTime = new TGraph();

  // Loop over all the time slices
  double t, f;
  int iStartPt = 0;

  while(iStartPt <= inGraph->GetN()-1){

    inGraph->GetPoint(iStartPt,t,f);
    double minTime = t;
    double maxTime = t + window;
    int iEndPt = inGraph->GetN() - 1;
    for(int iPt = iStartPt; iPt < inGraph->GetN(); iPt++){
      inGraph->GetPoint(iPt,t,f);
      if(t > maxTime) {
  iEndPt = iPt - 1;
  break;
      }
    }
    inGraph->GetPoint(iEndPt, maxTime, f);
    int N = iEndPt - iStartPt + 1;    

    // Calculate mean
    double sum = 0;
    for(int n = iStartPt; n < iEndPt; n++){
      inGraph->GetPoint(n, t, f);
      sum += f;
    }
    double mean = sum / N;
    
    // Calculate sum of sin and cosine terms to get tau
    double sumSin(0), sumCos(0);
    for(int n = iStartPt; n < iEndPt; n++){
      inGraph->GetPoint(n, t, f);
      sumSin += sin(TMath::TwoPi()*t);
      sumCos += cos(TMath::TwoPi()*t);
    }
    double tau = 1./(2*omega)*atan(sumSin/sumCos);
    
    // Calculate all the ingredients for the power
    double cosNum(0), cosDen(0), sinNum(0), sinDen(0);
    for(int n = iStartPt; n < iEndPt; n++){
      inGraph->GetPoint(n, t, f);
      cosNum += (f - mean)*cos(omega*(t - tau));
      cosDen += pow(cos(omega*(t - tau)),2);
      sinNum += (f - mean)*sin(omega*(t - tau));
      sinDen += pow(sin(omega*(t - tau)),2);
    }

    // Calculate sqrt of power (this is sort of amplitude)
    double power = sqrt(cosNum*cosNum/cosDen + sinNum*sinNum/sinDen);
    powerVsTime->SetPoint(powerVsTime->GetN(), 0.5*(minTime+maxTime), power);
    
    // Update start point
    iStartPt = iEndPt + 1;
  }

  TCanvas* c2 = new TCanvas();
  powerVsTime->SetMarkerStyle(20);
  // powerVsTime->SetTitle(";Time [us];Amplitude [arb. units]");
  powerVsTime->SetTitle(Form(";Time [#mus];Amplitude [arb. units], %0.4g #mus Window",window));
  powerVsTime->Draw("AP");

}



int combineAverageFitStartScans(){

	string file1 = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/Endgame/RandSeeds/FitStartScans-LongerRange/LombScargle/averageFitStartScan.root";
	string file2 = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/Endgame/RandSeeds/FitStartScans-LaterRange/FFT/averageFitStartScan.root";

    TFile *firstFile = TFile::Open(file1.c_str());
	TFile *secondFile = TFile::Open(file2.c_str());
	if (firstFile == 0 || secondFile == 0) {
	  printf("Error: cannot open file\n");
	  return 0;
	}


	TGraph* firstDiffGraph = (TGraph*) firstFile->Get("TRAverageDifference");

	TGraph* secondDiffGraph = (TGraph*) secondFile->Get("TRAverageDifference");
	secondDiffGraph->SetLineColor(2);
	secondDiffGraph->SetMarkerColor(2);

	auto canv1 = new TCanvas();
	firstDiffGraph->Draw();
	secondDiffGraph->Draw("PLSAME");

	double maxPointX, maxPointY;
	firstDiffGraph->GetPoint(firstDiffGraph->GetN()-1, maxPointX, maxPointY);

	TGraph* firstGraphCopy = (TGraph*) firstDiffGraph->Clone();

	for (int i = 0; i < secondDiffGraph->GetN(); ++i)
	{
		double secondGraphX, secondGraphY;
		secondDiffGraph->GetPoint(i, secondGraphX, secondGraphY);

		if(secondGraphX <= maxPointX) continue;
		else firstGraphCopy->SetPoint(firstGraphCopy->GetN(), secondGraphX, secondGraphY);
	}

	auto canv2 = new TCanvas();
	firstGraphCopy->Draw();

/////////////////////////////////////////////////////////////////////////////////////

      lombScargle(firstGraphCopy, 4);
      lombScargle(firstGraphCopy, 6);
      lombScargle(firstGraphCopy, 2*g2Period/1000);
      lombScargle(firstGraphCopy, 10);

/////////////////////////////////////////////////////////////////////////////////////

	return 0;
}
