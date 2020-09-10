#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <sstream>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveStats.h>

using namespace std;

vector<string> datasets {"60h", "HK", "9d", "EG"};
vector<string> analyses {"Nick T", "Nick R", "David T", "David A", "Aaron T", "Aaron A", "Matteo T", "Matteo A", "Bingzhi T", "Bingzhi A", "Tim Q"};

vector<double> Rvalues_60h {-28.802286, -28.966786, -28.211122, -28.228822, -28.619947, -28.637307, -28.884793, -28.481292, -28.739823, -28.422823, -29.206209};
vector<double> Rvalues_HK  {-27.044193, -27.209293, -27.209341, -26.946590, -27.004863, -26.965693, -27.080591, -27.021291, -27.001927, -27.091047, -24.946396};
vector<double> Rvalues_9d  {-27.917094, -27.921794, -28.020191, -27.553151, -27.898924, -27.572324, -27.888992, -27.599792, -27.893477, -27.744007, -26.279397};
vector<double> Rvalues_EG  {-27.701994, -27.765394, -27.715241, -27.590211, -27.714364, -27.669444, -27.877192, -27.727592, -27.665847, -27.694497, -27.990498};

vector<double> Rerrors_60h {1.358170,	1.359810, 1.337677, 1.207920, 1.330790, 1.218350, 1.332700, 1.193750, 1.331350, 1.206100, 2.058500};
vector<double> Rerrors_HK  {1.156130,	1.157420, 1.133584, 1.023270, 1.127740, 1.030180, 1.120300, 1.012020, 1.128100, 1.022280, 1.747800};
vector<double> Rerrors_9d  {0.930105,	0.932710, 0.912601, 0.823995, 0.907900, 0.830450, 0.906720, 0.814609, 0.908410, 0.822380, 1.403200};
vector<double> Rerrors_EG  {0.758400,	0.757600, 0.747400, 0.675823, 0.743690, 0.679850, 0.743493, 0.667947, 0.744140, 0.672860, 1.269000};


double calcAllowedDiff(double err1, double err2, double corr){
	return sqrt(err1*err1 + err2*err2 - 2*corr*err1*err2);
}


int calcRDiffs(){

  TFile* outputFile = new TFile("datasetDiffs.root","RECREATE");

  TH1F* allSigmas = new TH1F("allSigmas", "allSigmas", 100, -5, 5);

  vector<string> fileList; // make sure file order is 60h, HK, 9d, EG

   fileList.push_back("/gm2/data/users/nkinnaird/Ratio/ToyMC/DifferentAnalyzers/60h/Average/correlationAverage_60h.root");
   fileList.push_back("/gm2/data/users/nkinnaird/Ratio/ToyMC/DifferentAnalyzers/HK/Average/correlationAverage_HK.root");
   fileList.push_back("/gm2/data/users/nkinnaird/Ratio/ToyMC/DifferentAnalyzers/9d/Average/correlationAverage_9d.root");
   fileList.push_back("/gm2/data/users/nkinnaird/Ratio/ToyMC/DifferentAnalyzers/EG/Average/correlationAverage_EG.root");

	vector<vector<double>> Rvalues, Rerrors;

	Rvalues.push_back(Rvalues_60h);
	Rerrors.push_back(Rerrors_60h);

	Rvalues.push_back(Rvalues_HK);
	Rerrors.push_back(Rerrors_HK);

	Rvalues.push_back(Rvalues_9d);
	Rerrors.push_back(Rerrors_9d);

	Rvalues.push_back(Rvalues_EG);
	Rerrors.push_back(Rerrors_EG);


for (uint fileNum = 0; fileNum < fileList.size(); ++fileNum)
{
	TFile* inputFile = TFile::Open(fileList.at(fileNum).c_str());
	if(inputFile == 0){
		cout << "Can't open file" << endl;
		return 0;
	}

	TH2D* analyzerCorrelations = (TH2D*) inputFile->Get("Avg_CorrelationMatrix_R_R");
	// TH2D* reconCorrelations = (TH2D*) inputFile->Get("Avg_Recon_CorrelationMatrix_R_R");

	TH2D* allowedDiffHist = (TH2D*) analyzerCorrelations->Clone();
	TH2D* sigmaDiffHist = (TH2D*) analyzerCorrelations->Clone();

	allowedDiffHist->SetName((datasets.at(fileNum) + "_allowed_diff").c_str());
	sigmaDiffHist->SetName((datasets.at(fileNum) + "_sigma_diff").c_str());

	allowedDiffHist->GetZaxis()->SetTitle("Allowed Diff [ppm]");
	sigmaDiffHist->GetZaxis()->SetTitle("# Sigmas Diff");

	cout << endl << datasets.at(fileNum) << endl;

	for (int binY = 1; binY <= analyzerCorrelations->GetNbinsY(); ++binY){
		for (int binX = 1; binX <= analyzerCorrelations->GetNbinsX(); ++binX){
			
			double corr = analyzerCorrelations->GetBinContent(binX, binY);
			double Rdiff = Rvalues.at(fileNum).at(binX-1) - Rvalues.at(fileNum).at(binY-1);
			double allowedDiff = 0; 
			double sigmaDiff = 0;

			if(binX != binY){
				allowedDiff = calcAllowedDiff(Rerrors.at(fileNum).at(binX-1), Rerrors.at(fileNum).at(binY-1), corr);
				sigmaDiff = Rdiff/allowedDiff;
			}

			allowedDiffHist->SetBinContent(binX, binY, allowedDiff);
			sigmaDiffHist->SetBinContent(binX, binY, sigmaDiff);

			if(binX > binY) allSigmas->Fill(sigmaDiff);

            cout << std::setprecision(4) << setiosflags(ios::fixed) << std::noshowpos << "& " << allowedDiff;
            cout << std::setprecision(2) << setiosflags(ios::fixed) << std::showpos << " $" << sigmaDiff << "$ ";

			// cout << analyses.at(binX-1) << " " << analyses.at(binY-1) << " : corr: " << corr << " R diff: " << Rdiff << " allowedDiff: " << allowedDiff << " sigmaDiff: " << sigmaDiff << endl;
		}
		cout << " \\\\ " << endl;
	}

	outputFile->cd();

	allowedDiffHist->Write();
	sigmaDiffHist->Write();

	inputFile->Close();

} // end fileNum loop

	gStyle->SetOptStat(1110);
	gStyle->SetOptTitle(0);

	allSigmas->GetXaxis()->SetTitle("#DeltaR/#sigma_{allowed}");
	allSigmas->GetYaxis()->SetTitle("Entries");

    TCanvas* canv = new TCanvas("canv","canv",50,10,600,500);
    canv->SetRightMargin(0.1);
    allSigmas->Draw("HIST");

    canv->Update();

    TPaveStats* statsBox = (TPaveStats*) allSigmas->GetListOfFunctions()->FindObject("stats");
   	statsBox->SetBorderSize(1);

   	canv->SaveAs("AllSigmas.png");

    allSigmas->Write();

	// outputFile->Close();

	return 1;
}
