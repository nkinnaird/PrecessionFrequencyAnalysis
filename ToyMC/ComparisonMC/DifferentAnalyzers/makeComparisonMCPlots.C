#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TVectorD.h>
#include <TLegend.h>
#include <TLine.h>

using namespace std;

// choose which dataset plots to make
int datasetNum = 2; // 60h - 0, HK - 1, 9d - 2, EG - 3
vector<string> datasetStrings {"60h", "HK", "9d", "EG"};

double defaultLifetime = 64440;
double pi = 3.14159265358979323846;
double blindingFa = 0.2291 * 1e6 * 1/(1e9); // 0.2291 MHz converted to units of 1/ns, 0.2291 is the value used in the blinding software
double blindingWa = 2*pi*blindingFa;

double approximateEndTime = 700000;

/////////////////////////////////////////////////////////////////////////////////////

void checkFile(TFile* file, string path){
  if(file == 0){
    cout << endl << "Error: cannot open file: " << path << endl;
    exit(0);
  }
}

// set up according to David's file naming structure
vector<vector<TH1D*>> readEastFiles(vector<string> fileList){

  vector<vector<TH1D*>> inputEastHistograms; // outer vector is file, inner vector is N, A, Phi

  for (uint file = 0; file < fileList.size(); ++file)
  {
    TFile* inputFile;

    string path = fileList.at(file);

    inputFile = TFile::Open(path.c_str());
    checkFile(inputFile, path);

    string d_datasetString = path.substr(path.find_last_of(".")-3, 3);

    vector<TH1D*> fileHists;
    fileHists.push_back((TH1D*) (inputFile->Get((d_datasetString + "_t_h_no_mu").c_str())->Clone()));
    fileHists.push_back((TH1D*) (inputFile->Get((d_datasetString + "_t_h_as_mu").c_str())->Clone()));
    fileHists.push_back((TH1D*) (inputFile->Get((d_datasetString + "_t_h_ph_mu").c_str())->Clone()));

    for (uint i = 0; i < fileHists.size(); ++i) fileHists.at(i)->SetDirectory(0);

    inputEastHistograms.push_back(fileHists);
    inputFile->Close();
  }

  return inputEastHistograms;
}

// set up according to Matteo's file naming structure
vector<vector<TH1D*>> readWestFiles(vector<string> fileList){

  vector<vector<TH1D*>> inputWestHistograms; // outer vector is file, inner vector is N, A, Phi

  for (uint file = 0; file < fileList.size(); ++file)
  {
    TFile* inputFile;

    string path = fileList.at(file);

    inputFile = TFile::Open(path.c_str());
    checkFile(inputFile, path);

    vector<TH1D*> fileHists;
    fileHists.push_back((TH1D*) (inputFile->Get("hN")->Clone( (datasetStrings.at(file) + "_hN").c_str() )));
    fileHists.push_back((TH1D*) (inputFile->Get("hAsy")->Clone( (datasetStrings.at(file) + "_hAsy").c_str() )));
    fileHists.push_back((TH1D*) (inputFile->Get("hPhi")->Clone( (datasetStrings.at(file) + "_hPhi").c_str() )));

    for (uint i = 0; i < fileHists.size(); ++i) fileHists.at(i)->SetDirectory(0);

    inputWestHistograms.push_back(fileHists);
    inputFile->Close();
  }

  return inputWestHistograms;
}


// set up according to Josh's file naming structure
vector<TH2I*> readComparisonFiles(vector<string> fileList){

  vector<TH2I*> comparisonHists;

  for (uint file = 0; file < fileList.size(); ++file)
  {
    TFile* inputFile;

    string path = fileList.at(file);

    inputFile = TFile::Open(path.c_str());
    checkFile(inputFile, path);

    comparisonHists.push_back((TH2I*) (inputFile->Get("farline/evwEnergyEvW")->Clone((datasetStrings.at(file) + "_evwEnergyEvW").c_str())));
    comparisonHists.at(file)->SetDirectory(0);
    inputFile->Close();
  }

  return comparisonHists;
}

int getAnalyzerNumBins(double anaBinWidth, double anaBinEdge){
  return int((approximateEndTime - anaBinEdge)/anaBinWidth);
}

/////////////////////////////////////////////////////////////////////////////////////


class MCGeneratorFunc{

public:

  MCGeneratorFunc(){};

  void setHistograms(TH1D* N_hist_in, TH1D* A_hist_in, TH1D* Phi_hist_in){
    N_hist = N_hist_in;
    A_hist = A_hist_in;
    Phi_hist = Phi_hist_in;
  };

  double getAsymmetry(double energy){
    return A_hist->Interpolate(energy);
  };

  TH1D* getNHist(){ return N_hist; };
  TH1D* getAHist(){ return A_hist; };
  TH1D* getPhiHist(){ return Phi_hist; };

  double full2DFunc(double* x, double* p){
    double time = x[0];
    double energy = x[1];

    double N0 = N_hist->Interpolate(energy);
    double A = A_hist->Interpolate(energy);
    double phi = Phi_hist->Interpolate(energy);
    
    return N0 * exp(-time/defaultLifetime) * (1 + A * cos(blindingWa * time + phi));
  };

private:

  TH1D* N_hist;
  TH1D* A_hist;
  TH1D* Phi_hist;
};

/////////////////////////////////////////////////////////////////////////////////////

int makeComparisonMCPlots(){
  
  TH1::SetDefaultSumw2();

/////////////////////////////////////////////////////////////////////////////////////

  // Read in files containing energy bin input functions

  // Recon East - From David

  vector<string> fileList_east_eBin; // make sure file order is 60h, HK, 9d, EG

   fileList_east_eBin.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/SweigartEnergyBinFunctions/300-MeV-LowRange/fornick_60h.root");
   fileList_east_eBin.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/SweigartEnergyBinFunctions/300-MeV-LowRange/fornick_hik.root");
   fileList_east_eBin.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/SweigartEnergyBinFunctions/300-MeV-LowRange/fornick_9dy.root");
   fileList_east_eBin.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/SweigartEnergyBinFunctions/300-MeV-LowRange/fornick_end.root");

  vector<vector<TH1D*>> eastHistograms = readEastFiles(fileList_east_eBin);

  // Recon West - From Matteo

  vector<string> fileList_west_eBin; // make sure file order is 60h, HK, 9d, EG

   fileList_west_eBin.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/MatteoEnergyBinFunctions/fit_ebinned_60h_Like_60MeV_binAna.root");
   fileList_west_eBin.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/MatteoEnergyBinFunctions/fit_ebinned_HK_Like_60MeV_binAna.root");
   fileList_west_eBin.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/MatteoEnergyBinFunctions/fit_ebinned_9d_Like_60MeV_binAna.root");
   fileList_west_eBin.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/MatteoEnergyBinFunctions/fit_ebinned_EG_Like_60MeV_binAna.root");

  vector<vector<TH1D*>> westHistograms = readWestFiles(fileList_west_eBin);


/////////////////////////////////////////////////////////////////////////////////////

  // make 2D function for hit generation

  vector<MCGeneratorFunc*> funcClass_east;
  vector<MCGeneratorFunc*> funcClass_west;

  vector<TF2*> dataMCFunc_east;
  vector<TF2*> dataMCFunc_west;

  for (uint dataset = 0; dataset < datasetStrings.size(); ++dataset)
  {
    funcClass_east.push_back(new MCGeneratorFunc());
    funcClass_west.push_back(new MCGeneratorFunc());

    funcClass_east.at(dataset)->setHistograms(eastHistograms.at(dataset).at(0), eastHistograms.at(dataset).at(1), eastHistograms.at(dataset).at(2));
    funcClass_west.at(dataset)->setHistograms(westHistograms.at(dataset).at(0), westHistograms.at(dataset).at(1), westHistograms.at(dataset).at(2));


    double minEnergy_east = funcClass_east.at(dataset)->getNHist()->GetXaxis()->GetXmin();
    double maxEnergy_east = funcClass_east.at(dataset)->getNHist()->GetXaxis()->GetXmax();
    double eBinWidth_east = funcClass_east.at(dataset)->getNHist()->GetBinWidth(1);
    int numEnergyBins_east = funcClass_east.at(dataset)->getNHist()->GetNbinsX();
    int nPy_east = numEnergyBins_east * 1;

    double minEnergy_west = funcClass_west.at(dataset)->getNHist()->GetXaxis()->GetXmin();
    double maxEnergy_west = funcClass_west.at(dataset)->getNHist()->GetXaxis()->GetXmax();
    double eBinWidth_west = funcClass_west.at(dataset)->getNHist()->GetBinWidth(1);
    int numEnergyBins_west = funcClass_west.at(dataset)->getNHist()->GetNbinsX();
    int nPy_west = numEnergyBins_west * 1;

    cout << endl << "Generating 2D functions for dataset: " << datasetStrings.at(dataset) << endl;

    cout << "East minEnergy_east: " << minEnergy_east << " maxEnergy_east: " << maxEnergy_east <<  " width: " << eBinWidth_east << " num Energy bins: " << numEnergyBins_east << endl;
    cout << "West minEnergy_east: " << minEnergy_west << " maxEnergy_west: " << maxEnergy_west <<  " width: " << eBinWidth_west << " num Energy bins: " << numEnergyBins_west << endl;

    double timePointWidth = 149.2/2;
    int nPx = getAnalyzerNumBins(timePointWidth, 0);
    double functionEndRange = nPx * timePointWidth;

    cout << "East: Number time function points: " << nPx << " function time end range: " << functionEndRange << " energy points: " << nPy_east << endl;
    cout << "West: Number time function points: " << nPx << " function time end range: " << functionEndRange << " energy points: " << nPy_east << endl;

    dataMCFunc_east.push_back(new TF2(Form("dataMCFunc_east_%s", datasetStrings.at(dataset).c_str()), funcClass_east.at(dataset), &MCGeneratorFunc::full2DFunc, 0, functionEndRange, minEnergy_east, maxEnergy_east, 0));    
    dataMCFunc_east.at(dataset)->SetNpx(nPx);
    dataMCFunc_east.at(dataset)->SetNpy(nPy_east);

    dataMCFunc_west.push_back(new TF2(Form("dataMCFunc_west_%s", datasetStrings.at(dataset).c_str()), funcClass_west.at(dataset), &MCGeneratorFunc::full2DFunc, 0, functionEndRange, minEnergy_west, maxEnergy_west, 0));    
    dataMCFunc_west.at(dataset)->SetNpx(nPx);
    dataMCFunc_west.at(dataset)->SetNpy(nPy_west);
  }




  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  // make 2D function plot

  auto canv_2dfunc = new TCanvas("canv_2dfunc", "canv_2dfunc", 10,10,800,600);
  dataMCFunc_east.at(datasetNum)->GetXaxis()->SetTitle("Time [ns]");
  dataMCFunc_east.at(datasetNum)->GetYaxis()->SetTitle("Energy [MeV]");
  dataMCFunc_east.at(datasetNum)->Draw("COLZ");
  canv_2dfunc->SetLogz();
  canv_2dfunc->Update();

  canv_2dfunc->SaveAs("2D_function.png");

  // make energy bin function plots

  // N 

  TH1D* plotN_east = (TH1D*) eastHistograms.at(datasetNum).at(0)->Clone();
  TH1D* plotN_west = (TH1D*) westHistograms.at(datasetNum).at(0)->Clone();

  plotN_east->Scale(1/plotN_east->Integral("width"));
  plotN_west->Scale(1/plotN_west->Integral("width"));

  plotN_east->GetXaxis()->SetTitle("Energy [MeV]");
  plotN_east->GetXaxis()->SetMaxDigits(4);
  plotN_east->GetYaxis()->SetTitle("N(E) (Normalized)");

  plotN_west->SetLineColor(2);

  auto canv_N = new TCanvas("canv_N", "canv_N", 100,10,800,600);
  canv_N->SetRightMargin(0.05);
  canv_N->SetTopMargin(0.1);

  plotN_east->Draw("HIST");
  plotN_west->Draw("HISTSAME");

  auto legend_N = new TLegend(0.65,0.675,0.925,0.825);

  legend_N->AddEntry(plotN_east, "ReconEast","l");
  legend_N->AddEntry(plotN_west, "ReconWest","l");

  legend_N->SetBorderSize(0);
  legend_N->SetFillStyle(1001);
  legend_N->Draw();

  canv_N->SaveAs("ReconEastvWest_N.png");

  // A 

  TH1D* plotA_east = (TH1D*) eastHistograms.at(datasetNum).at(1)->Clone();
  TH1D* plotA_west = (TH1D*) westHistograms.at(datasetNum).at(1)->Clone();

  plotA_east->GetXaxis()->SetTitle("Energy [MeV]");
  plotA_east->GetXaxis()->SetMaxDigits(4);
  plotA_east->GetYaxis()->SetTitle("A(E)");

  plotA_west->SetLineColor(2);

  auto canv_A = new TCanvas("canv_A", "canv_A", 200,10,800,600);
  canv_A->SetRightMargin(0.05);
  canv_A->SetTopMargin(0.1);

  plotA_east->Draw("HIST");
  plotA_west->Draw("HISTSAME");

  auto legend_A = new TLegend(0.275,0.675,0.55,0.825);

  legend_A->AddEntry(plotA_east, "ReconEast","l");
  legend_A->AddEntry(plotA_west, "ReconWest","l");

  legend_A->SetBorderSize(0);
  legend_A->SetFillStyle(1001);
  legend_A->Draw();

  canv_A->SaveAs("ReconEastvWest_A.png");

  // Phi 

  TH1D* plotPhi_east = (TH1D*) eastHistograms.at(datasetNum).at(2)->Clone();
  TH1D* plotPhi_west = (TH1D*) westHistograms.at(datasetNum).at(2)->Clone();

  plotPhi_east->GetXaxis()->SetTitle("Energy [MeV]");
  plotPhi_east->GetXaxis()->SetMaxDigits(4);
  plotPhi_east->GetYaxis()->SetTitle("#phi(E)");

  plotPhi_west->SetLineColor(2);

  auto canv_Phi = new TCanvas("canv_Phi", "canv_Phi", 300,10,800,600);
  canv_Phi->SetRightMargin(0.05);
  canv_Phi->SetTopMargin(0.1);

  plotPhi_east->Draw("HIST");
  plotPhi_west->Draw("HISTSAME");

  auto legend_Phi = new TLegend(0.65,0.675,0.925,0.825);

  legend_Phi->AddEntry(plotPhi_east, "ReconEast","l");
  legend_Phi->AddEntry(plotPhi_west, "ReconWest","l");

  legend_Phi->SetBorderSize(0);
  legend_Phi->SetFillStyle(1001);
  legend_Phi->Draw();

  canv_Phi->SaveAs("ReconEastvWest_Phi.png");



/////////////////////////////////////////////////////////////////////////////////////

  // input Josh's functions for energy comparison

  vector<string> fileList_comparison; // make sure file order is 60h, HK, 9d, EG

   fileList_comparison.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/JoshEnergyComparison/results_EvW_60h_Feb5_CorrectThreshold.root");
   fileList_comparison.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/JoshEnergyComparison/results_EvW_HighKick_Feb5.root");
   fileList_comparison.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/JoshEnergyComparison/results_EvW_9day_Feb3_CombinedRecoveryAndOriginal.root");
   fileList_comparison.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/JoshEnergyComparison/results_EvW_EndGame_Feb3_BeforeRecovery.root");

  vector<TH2I*> comparisonHistograms = readComparisonFiles(fileList_comparison);

  int reconComparisonEnergyBins = comparisonHistograms.at(0)->GetNbinsX(); // these binning parameters should all be the same, in X and Y and for the different datasets
  double reconComparisonEnergyBinWidth = comparisonHistograms.at(0)->GetXaxis()->GetBinWidth(1);

  cout << "Num recon comparison energy bins: " << reconComparisonEnergyBins << " with width: " << reconComparisonEnergyBinWidth << endl;

  TH1D* EtoWProjections[4][reconComparisonEnergyBins];
  TH1D* WtoEProjections[4][reconComparisonEnergyBins];

  for (int dataset = 0; dataset < 4; ++dataset)
  {
    for (int bin = 1; bin <= reconComparisonEnergyBins; ++bin)
    {
      int eHigh = int(bin * reconComparisonEnergyBinWidth);
      int eLow = int(eHigh - reconComparisonEnergyBinWidth);
      EtoWProjections[dataset][bin-1] = comparisonHistograms.at(dataset)->ProjectionY(Form("%s_Proj_EtoW_%i_%i",datasetStrings.at(dataset).c_str(),eLow,eHigh), bin, bin);
      WtoEProjections[dataset][bin-1] = comparisonHistograms.at(dataset)->ProjectionX(Form("%s_Proj_WtoE_%i_%i",datasetStrings.at(dataset).c_str(),eLow,eHigh), bin, bin);
    }
  }



  // make energy comparison plots

  auto canv_comparison = new TCanvas("canv_comparison", "canv_comparison", 400,10,800,600);

  comparisonHistograms.at(datasetNum)->Draw("COLZ");
  comparisonHistograms.at(datasetNum)->GetXaxis()->SetMaxDigits(4);
  comparisonHistograms.at(datasetNum)->GetYaxis()->SetMaxDigits(4);
  canv_comparison->SetLogz();
  canv_comparison->Update();

  canv_comparison->SaveAs("ReconEastvWest_Energies.png");

  // projection

  auto canv_projection = new TCanvas("canv_projection", "canv_projection", 500,10,800,600);
  canv_projection->SetRightMargin(0.1);

  double energyForPlot = 2000;
  int plotBin = energyForPlot/reconComparisonEnergyBinWidth;

  EtoWProjections[datasetNum][plotBin]->Draw("HIST");
  EtoWProjections[datasetNum][plotBin]->GetXaxis()->SetMaxDigits(4);
  EtoWProjections[datasetNum][plotBin]->GetYaxis()->SetTitle("Entries");
  canv_projection->SetLogy();
  canv_projection->Update();

  cout << "max and min: " << canv_projection->GetUymin() << " " << canv_projection->GetUymax() << endl;

  // for some reason the canvas uy min and max values aren't right

  // TLine *line = new TLine(energyForPlot, canv_projection->GetUymin(), energyForPlot, canv_projection->GetUymax());
  TLine *line = new TLine(energyForPlot, 0, energyForPlot, 2e6);
  line->SetLineColor(2);
  line->SetLineStyle(2);
  line->Draw();

  auto legend_REenergy = new TLegend(0.525,0.675,0.825,0.825);

  legend_REenergy->AddEntry(line, "ReconEast Energy","l");

  legend_REenergy->SetBorderSize(0);
  legend_REenergy->SetFillStyle(1001);
  legend_REenergy->Draw();

  canv_projection->SaveAs("ReconEastvWest_Projection.png");

  


 	return 1;

}
