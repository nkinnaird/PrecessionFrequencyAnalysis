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
#include <TRandomGen.h>

#include <chrono> 

using namespace std;

bool runOverAllDatasets = false; // set to true to run over all 4 datasets, the output will be put into samples 0 through 3, datasetStatistics will overwrite nPts, and numberIndependentSamples will be set to 4, numberRandomSeeds still needs to be specified
vector<double> datasetStatistics {5.04e9, 7.01e9, 1.08e10, 2.21e10};
double statsReduction = 1; // change to reduce the amount of stats used when running over all the datasets - for testing

double nPts = 1.08e10;

// decide the dataset to simulate - unless runOverAllDatasets is set to true
int datasetNum = 2; // 60h - 0, HK - 1, 9d - 2, EG - 3
vector<string> datasetStrings {"60h", "HK", "9d", "EG"};

bool eastToWest = false; // set to false to go from west to east

int numberIndependentSamples = 1;
int numberRandomSeeds = 100; // just for the Ratio method

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
vector<vector<TH1D*>> readEastFiles(vector<string> fileList, int grid){

  vector<vector<TH1D*>> inputEastHistograms; // outer vector is file, inner vector is N, A, Phi

  for (uint file = 0; file < fileList.size(); ++file)
  {
    TFile* inputFile;

    string path = fileList.at(file);
    if(grid == 1) path = path.substr(path.find_last_of("/")+1);

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
vector<vector<TH1D*>> readWestFiles(vector<string> fileList, int grid){

  vector<vector<TH1D*>> inputWestHistograms; // outer vector is file, inner vector is N, A, Phi

  for (uint file = 0; file < fileList.size(); ++file)
  {
    TFile* inputFile;

    string path = fileList.at(file);
    if(grid == 1) path = path.substr(path.find_last_of("/")+1);

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
vector<TH2I*> readComparisonFiles(vector<string> fileList, int grid){

  vector<TH2I*> comparisonHists;

  for (uint file = 0; file < fileList.size(); ++file)
  {
    TFile* inputFile;

    string path = fileList.at(file);
    if(grid == 1) path = path.substr(path.find_last_of("/")+1);

    inputFile = TFile::Open(path.c_str());
    checkFile(inputFile, path);

    comparisonHists.push_back((TH2I*) (inputFile->Get("farline/evwEnergyEvW")->Clone((datasetStrings.at(file) + "_evwEnergyEvW").c_str())));
    comparisonHists.at(file)->SetDirectory(0);
    inputFile->Close();
  }

  return comparisonHists;
}

/////////////////////////////////////////////////////////////////////////////////////

// TUUID implementation in order to generate a unique random seed - code copied from TRandom.SetSeed(0) - this just generates a unique random seed for the random generator used later, which is NOT TRandom
uint getUniqueSeed(){
  TUUID u;
  UChar_t uuid[16];
  u.GetUUID(uuid);
  uint seed = uint(uuid[3])*16777216 + uint(uuid[2])*65536 + uint(uuid[1])*256 + uint(uuid[0]);
  return seed;
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

// set grid to 1 so that input file paths are updated properly - grid submission scripts must tar up the macro and input files accordingly
int makeAnalyzerHists(int grid = 0){

  cout << endl << endl;

  if(runOverAllDatasets){
    cout << "Running over all datasets, ouput will be placed into samples 0 through 3 for the 60h, HighKick, 9d, and EndGame datasets respectively." << endl;
    numberIndependentSamples = 4;
  }
  else cout << "Running approximately " << nPts << " events for " << numberIndependentSamples << " samples for the dataset: " << datasetStrings.at(datasetNum) << endl;

  cout << "Running " << numberRandomSeeds << " ratio random seeds." << endl;

  if(eastToWest) cout << "Going from east to west." << endl;
  else cout << "Going from west to east." << endl;

  if(grid == 0) cout << "Running locally." << endl;
  else if(grid == 1) cout << "Grid implementation set to true, input file paths will be updated accordingly." << endl;
  else{
    cout << "Wrong macro input" << endl;
    return 0; 
  }

	bool timeCheck = true;
	clock_t t;
  t = clock();

/////////////////////////////////////////////////////////////////////////////////////

  // use TRandomMixMax instead of TRandom3 which has some deficiencies for statistical studies, use TRandomMixMax instead of the Ranlux generators since it is a lot faster
  // TRandomMixMax will produce a different set of random numbers for every input seed, it does something weird with an input seed of 0
  // something is off with the way that SetSeed works for TRandomGen and it's derived classes like TRandomMixMax, either a mistake or on purpose so that the different generators can be implemented differently
  // it doesn't properly change the internal variable fSeed so GetSeed calls on the generator return the default value of 65539 implemented in TRandom
  // setting the seed to 0 also doesn't implement the unique seed based on the TUUID, therefore I've added my own implementation to get a unique seed with which to start the generators

  // gRandom is used in GetRandom function calls in TF1, TH1, TH2, etc.
  // use it for the other randomization calls in this code as well
  if(gRandom) delete gRandom;
  gRandom = new TRandomMixMax(getUniqueSeed());
  // gRandom = new TRandomMixMax(123456); // static seed for tests

  cout << "gRandom generator: ";
  gRandom->Print();
  
/////////////////////////////////////////////////////////////////////////////////////

  TH1::SetDefaultSumw2();

/////////////////////////////////////////////////////////////////////////////////////

  // Read in files containing energy bin input functions

  // Recon East - From David

  vector<string> fileList_east_eBin; // make sure file order is 60h, HK, 9d, EG

   fileList_east_eBin.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/SweigartEnergyBinFunctions/300-MeV-LowRange/fornick_60h.root");
   fileList_east_eBin.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/SweigartEnergyBinFunctions/300-MeV-LowRange/fornick_hik.root");
   fileList_east_eBin.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/SweigartEnergyBinFunctions/300-MeV-LowRange/fornick_9dy.root");
   fileList_east_eBin.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/SweigartEnergyBinFunctions/300-MeV-LowRange/fornick_end.root");

  vector<vector<TH1D*>> eastHistograms = readEastFiles(fileList_east_eBin, grid);

  // Recon West - From Matteo

  vector<string> fileList_west_eBin; // make sure file order is 60h, HK, 9d, EG

   fileList_west_eBin.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/MatteoEnergyBinFunctions/fit_ebinned_60h_Like_60MeV_binAna.root");
   fileList_west_eBin.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/MatteoEnergyBinFunctions/fit_ebinned_HK_Like_60MeV_binAna.root");
   fileList_west_eBin.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/MatteoEnergyBinFunctions/fit_ebinned_9d_Like_60MeV_binAna.root");
   fileList_west_eBin.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/MatteoEnergyBinFunctions/fit_ebinned_EG_Like_60MeV_binAna.root");

  vector<vector<TH1D*>> westHistograms = readWestFiles(fileList_west_eBin, grid);


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


/////////////////////////////////////////////////////////////////////////////////////

  // input Josh's functions for energy comparison

  vector<string> fileList_comparison; // make sure file order is 60h, HK, 9d, EG

   fileList_comparison.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/JoshEnergyComparison/results_EvW_60h_Feb5_CorrectThreshold.root");
   fileList_comparison.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/JoshEnergyComparison/results_EvW_HighKick_Feb5.root");
   fileList_comparison.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/JoshEnergyComparison/results_EvW_9day_Feb3_CombinedRecoveryAndOriginal.root");
   fileList_comparison.push_back("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_21_06/srcs/gm2analyses/macros/RatioMacro/ToyMC/ComparisonMC/OtherAnalyzerFunctions/JoshEnergyComparison/results_EvW_EndGame_Feb3_BeforeRecovery.root");

  vector<TH2I*> comparisonHistograms = readComparisonFiles(fileList_comparison, grid);

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

/////////////////////////////////////////////////////////////////////////////////////

  TFile* outputFile = new TFile("analyzerHists.root","RECREATE");

  TVectorD numSamples(1);
  numSamples[0] = numberIndependentSamples;
  numSamples.Write("Samples");

  TVectorD numSeeds(1);
  numSeeds[0] = numberRandomSeeds;
  numSeeds.Write("Seeds");

  // write input information

  auto inputDir = gFile->mkdir("Input");
  inputDir->cd();

  for (uint dataset = 0; dataset < datasetStrings.size(); ++dataset)
  {
    auto dir_dataset = inputDir->mkdir(datasetStrings.at(dataset).c_str());
      dir_dataset->cd();
      dataMCFunc_east.at(dataset)->Write();
      dataMCFunc_west.at(dataset)->Write();

    auto east_dir = dir_dataset->mkdir("East");
      east_dir->cd();
      for (uint i = 0; i < eastHistograms.at(dataset).size(); ++i) eastHistograms.at(dataset).at(i)->Write();

    auto west_dir = dir_dataset->mkdir("West");
      west_dir->cd();
      for (uint i = 0; i < westHistograms.at(dataset).size(); ++i) westHistograms.at(dataset).at(i)->Write();

    auto vs_dir = dir_dataset->mkdir("EvW");
      vs_dir->cd();
      comparisonHistograms.at(dataset)->Write();
      // for (int bin = 1; bin <= reconComparisonEnergyBins; ++bin) EtoWProjections[dataset][bin-1]->Write();
      // for (int bin = 1; bin <= reconComparisonEnergyBins; ++bin) WtoEProjections[dataset][bin-1]->Write();
  }


  auto topDir = gFile->mkdir("topDir"); // generated histograms go under here


/////////////////////////////////////////////////////////////////////////////////////

  // define analyzer parameters and histograms

  // Nick 
  double nick_BinWidth = 149.2;
  double nick_BinEdge = 0;
  int nick_NumBins = getAnalyzerNumBins(nick_BinWidth, nick_BinEdge);
  double nick_T_LowThresh = 1700;
  double nick_T_HighThresh = 1e6;

  TH1I* nick_T_Hist[numberIndependentSamples];
  TH1I* nick_R_U_Hists[numberIndependentSamples][numberRandomSeeds];
  TH1I* nick_R_V_Hists[numberIndependentSamples][numberRandomSeeds];    

  // David 
  double david_BinWidth = 149.2;
  double david_BinEdge = 0;
  int david_NumBins = getAnalyzerNumBins(david_BinWidth, david_BinEdge);
  double david_T_LowThresh = 1700;
  double david_T_HighThresh = 1e6;
  double david_A_LowThresh = 1000;
  double david_A_HighThresh = 3000;

  TH1I* david_T_Hist[numberIndependentSamples];
  TH1D* david_A_Hist[numberIndependentSamples];

  // Aaron 
  double aaron_BinWidth = 149.19;
  double aaron_BinEdge = 53.62;
  int aaron_NumBins = getAnalyzerNumBins(aaron_BinWidth, aaron_BinEdge);
  double aaron_T_LowThresh = 1700;
  double aaron_T_HighThresh = 6000;
  double aaron_A_LowThresh = 1000;
  double aaron_A_HighThresh = 3020;

  TH1I* aaron_T_Hist[numberIndependentSamples];
  TH1D* aaron_A_Hist[numberIndependentSamples];

  // Matteo 
  double matteo_BinWidth = 149.19;
  double matteo_BinEdge = 0;
  int matteo_NumBins = getAnalyzerNumBins(matteo_BinWidth, matteo_BinEdge);
  double matteo_T_LowThresh = 1680;
  double matteo_T_HighThresh = 7020;
  double matteo_A_LowThresh = 1080;
  double matteo_A_HighThresh = 3020;

  TH1I* matteo_T_Hist[numberIndependentSamples];
  TH1D* matteo_A_Hist[numberIndependentSamples];

  // Bingzhi 
  double bingzhi_BinWidth = 149.2;
  double bingzhi_BinEdge = 0;
  int bingzhi_NumBins = getAnalyzerNumBins(bingzhi_BinWidth, bingzhi_BinEdge);
  double bingzhi_T_LowThresh = 1700;
  double bingzhi_T_HighThresh = 9300;
  double bingzhi_A_LowThresh = 1000;
  double bingzhi_A_HighThresh = 3100;

  TH1I* bingzhi_T_Hist[numberIndependentSamples];
  TH1D* bingzhi_A_Hist[numberIndependentSamples];

  // Tim 
  double tim_BinWidth = 150;
  double tim_BinEdge = 0;
  int tim_NumBins = getAnalyzerNumBins(tim_BinWidth, tim_BinEdge);
  double tim_Q_LowThresh = 300;
  double tim_Q_HighThresh = 1e6;

  TH1D* tim_Q_Hist[numberIndependentSamples];

/////////////////////////////////////////////////////////////////////////////////////

  TDirectory* sampleNumDir[numberIndependentSamples];

  for (int sampleNum = 0; sampleNum < numberIndependentSamples; ++sampleNum)
  {
    // Make directory for each sample number and change directory to it
    sampleNumDir[sampleNum] = topDir->mkdir(Form("SampleNum%d",sampleNum));
    sampleNumDir[sampleNum]->cd();

    nick_T_Hist[sampleNum] = new TH1I("nick_T_Hist", "nick_T_Hist", nick_NumBins, nick_BinEdge, nick_NumBins*nick_BinWidth + nick_BinEdge);

      auto ratioDir = sampleNumDir[sampleNum]->mkdir("Ratio");
      ratioDir->cd();

      for (int seedNum = 0; seedNum < numberRandomSeeds; ++seedNum)
      {
        nick_R_U_Hists[sampleNum][seedNum] = new TH1I(Form("nick_R_U_Hist_%i",seedNum), Form("nick_R_U_Hist_%i",seedNum), nick_NumBins, nick_BinEdge, nick_NumBins*nick_BinWidth + nick_BinEdge);
        nick_R_V_Hists[sampleNum][seedNum] = new TH1I(Form("nick_R_V_Hist_%i",seedNum), Form("nick_R_V_Hist_%i",seedNum), nick_NumBins, nick_BinEdge, nick_NumBins*nick_BinWidth + nick_BinEdge);
      }

    sampleNumDir[sampleNum]->cd();

    david_T_Hist[sampleNum] = new TH1I("david_T_Hist", "david_T_Hist", david_NumBins, david_BinEdge, david_NumBins*david_BinWidth + david_BinEdge);
    david_A_Hist[sampleNum] = new TH1D("david_A_Hist", "david_A_Hist", david_NumBins, david_BinEdge, david_NumBins*david_BinWidth + david_BinEdge);

    aaron_T_Hist[sampleNum] = new TH1I("aaron_T_Hist", "aaron_T_Hist", aaron_NumBins, aaron_BinEdge, aaron_NumBins*aaron_BinWidth + aaron_BinEdge);
    aaron_A_Hist[sampleNum] = new TH1D("aaron_A_Hist", "aaron_A_Hist", aaron_NumBins, aaron_BinEdge, aaron_NumBins*aaron_BinWidth + aaron_BinEdge);

    matteo_T_Hist[sampleNum] = new TH1I("matteo_T_Hist", "matteo_T_Hist", matteo_NumBins, matteo_BinEdge, matteo_NumBins*matteo_BinWidth + matteo_BinEdge);
    matteo_A_Hist[sampleNum] = new TH1D("matteo_A_Hist", "matteo_A_Hist", matteo_NumBins, matteo_BinEdge, matteo_NumBins*matteo_BinWidth + matteo_BinEdge);

    bingzhi_T_Hist[sampleNum] = new TH1I("bingzhi_T_Hist", "bingzhi_T_Hist", bingzhi_NumBins, bingzhi_BinEdge, bingzhi_NumBins*bingzhi_BinWidth + bingzhi_BinEdge);
    bingzhi_A_Hist[sampleNum] = new TH1D("bingzhi_A_Hist", "bingzhi_A_Hist", bingzhi_NumBins, bingzhi_BinEdge, bingzhi_NumBins*bingzhi_BinWidth + bingzhi_BinEdge);

    tim_Q_Hist[sampleNum] = new TH1D("tim_Q_Hist", "tim_Q_Hist", tim_NumBins, tim_BinEdge, tim_NumBins*tim_BinWidth + tim_BinEdge);

    // auto unique_Dir = sampleNumDir[sampleNum]->mkdir("UniqueHits");
    // unique_Dir->cd();
  }

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  // calculate probabilities for ratio histogram filling

  double gm2HalfPeriod = (1/blindingFa)/2;

  double totalChance = exp(gm2HalfPeriod/defaultLifetime) + exp(-gm2HalfPeriod/defaultLifetime) + 2.;
  double percentChanceUPlus = exp(gm2HalfPeriod/defaultLifetime) / totalChance;
  double percentChanceUMinus = exp(-gm2HalfPeriod/defaultLifetime) / totalChance;

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

cout << endl << endl;

clock_t startTime = clock();

for (int sampleNum = 0; sampleNum < numberIndependentSamples; ++sampleNum){

  cout << endl << "Beginning sample " << sampleNum << endl;

  if(sampleNum > 0){
    double elapsedTime = double(clock() - startTime)/CLOCKS_PER_SEC;
    double timeRemainingSeconds = (double(numberIndependentSamples) - sampleNum)/sampleNum * elapsedTime;
    cout << "Samples done: " << sampleNum << "/" << numberIndependentSamples << ". Time passed: " << int(elapsedTime/60/60) << " hours, " << int(elapsedTime/60)%60 << " minutes, and " << int(elapsedTime)%60 << " seconds." << " Estimated time remaining = " << int(timeRemainingSeconds/60/60) << " hours, " << int(timeRemainingSeconds/60)%60 << " minutes, and " << int(timeRemainingSeconds)%60 << " seconds." << endl;
    if(runOverAllDatasets) cout << "Note that runOverAllDatasets is set to true, so the running will take longer than estimated here." << endl;
  }

  double randEntries = gRandom->PoissonD(nPts);

  if(runOverAllDatasets){
    datasetNum = sampleNum;
    randEntries = gRandom->PoissonD(datasetStatistics.at(datasetNum)/statsReduction);
  } 

  cout << "Simulating dataset: " << datasetStrings.at(datasetNum) << endl;
  cout << "Simulating " << randEntries << " entries." << endl;

  for (double entry = 0; entry < randEntries; ++entry)
  {
    if(fmod(entry, 1000000) == 0) cout << "At entry num: " << entry << endl;

    double time, energy;
    double energyEast, energyWest, asymmetryEast, asymmetryWest;

    if(eastToWest){
      dataMCFunc_east.at(datasetNum)->GetRandom2(time, energy);
      energyEast = energy;
      asymmetryEast = funcClass_east.at(datasetNum)->getAsymmetry(energyEast);

      int projectionNum = energyEast/reconComparisonEnergyBinWidth;
      energyWest = EtoWProjections[datasetNum][projectionNum]->GetRandom();
      asymmetryWest = funcClass_west.at(datasetNum)->getAsymmetry(energyWest);
    }
    else{
      dataMCFunc_west.at(datasetNum)->GetRandom2(time, energy);
      energyWest = energy;
      asymmetryWest = funcClass_west.at(datasetNum)->getAsymmetry(energyWest);

      int projectionNum = energyWest/reconComparisonEnergyBinWidth;
      energyEast = WtoEProjections[datasetNum][projectionNum]->GetRandom();
      asymmetryEast = funcClass_east.at(datasetNum)->getAsymmetry(energyEast);
    }

    // cout << "Time: " << time << " energy: " << energy << " east energy: " << energyEast << " west energy: " << energyWest << " east asymmetry: " << asymmetryEast << " west asymmetry: " << asymmetryWest << endl;

/////////////////////////////////////////////////////////////////////////////////////
    // Recon West

      // fill Nick hists
      if(energyWest > nick_T_LowThresh && energyWest < nick_T_HighThresh){
        nick_T_Hist[sampleNum]->Fill(time, 1);

        // do Ratio method randomization
        for (int seedNum = 0; seedNum < numberRandomSeeds; ++seedNum)
        {
          double ratioProb = gRandom->Uniform();

          if     (ratioProb < percentChanceUPlus) nick_R_U_Hists[sampleNum][seedNum]->Fill(time - gm2HalfPeriod, 1);
          else if(ratioProb < (percentChanceUPlus + percentChanceUMinus)) nick_R_U_Hists[sampleNum][seedNum]->Fill(time + gm2HalfPeriod, 1);
          else if(ratioProb < 1) nick_R_V_Hists[sampleNum][seedNum]->Fill(time, 1);
        }
      }

      // fill Aaron hists
      if(energyWest > aaron_T_LowThresh && energyWest < aaron_T_HighThresh) aaron_T_Hist[sampleNum]->Fill(time, 1);
      if(energyWest > aaron_A_LowThresh && energyWest < aaron_A_HighThresh) aaron_A_Hist[sampleNum]->Fill(time, asymmetryWest);

      // fill Matteo hists
      if(energyWest > matteo_T_LowThresh && energyWest < matteo_T_HighThresh) matteo_T_Hist[sampleNum]->Fill(time, 1);
      if(energyWest > matteo_A_LowThresh && energyWest < matteo_A_HighThresh) matteo_A_Hist[sampleNum]->Fill(time, asymmetryWest);

      // fill Bingzhi hists
      if(energyWest > bingzhi_T_LowThresh && energyWest < bingzhi_T_HighThresh) bingzhi_T_Hist[sampleNum]->Fill(time, 1);
      if(energyWest > bingzhi_A_LowThresh && energyWest < bingzhi_A_HighThresh) bingzhi_A_Hist[sampleNum]->Fill(time, asymmetryWest);

/////////////////////////////////////////////////////////////////////////////////////
    // Recon East

      // fill David hists
      if(energyEast > david_T_LowThresh && energyEast < david_T_HighThresh) david_T_Hist[sampleNum]->Fill(time, 1);
      if(energyEast > david_A_LowThresh && energyEast < david_A_HighThresh) david_A_Hist[sampleNum]->Fill(time, asymmetryEast);

/////////////////////////////////////////////////////////////////////////////////////
    // Q

    // fill Tim hists
    if(energy > tim_Q_LowThresh && energy < tim_Q_HighThresh) tim_Q_Hist[sampleNum]->Fill(time, energy);

  } // end loop over hits

} // end loop over samples

/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////


  outputFile->Write();
  delete outputFile;

	if(timeCheck){
    	t = clock() - t;
    	printf ("It took me %f seconds.\n", ((float)t)/CLOCKS_PER_SEC);
 	}

 	return 1;

}
