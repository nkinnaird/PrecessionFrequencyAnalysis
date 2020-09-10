#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TH1.h>
#include <TFile.h>
#include <TLine.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TVectorD.h>
#include <TFitResult.h>
#include <TF1.h>
#include <TROOT.h>
#include <TGraphErrors.h>

using namespace std;

bool saveAllFits = true; // turn to false before fitting many files

double defaultLifetime = 64440;
double pi = 3.14159265358979323846;
double blindingFa = 0.2291 * 1e6 * 1/(1e9); // 0.2291 MHz converted to units of 1/ns, 0.2291 is the value used in the blinding software
double blindingWa = 2*pi*blindingFa;


double plot_cbo_freq = -1;
double plot_VW_freq = -1;
#include <TSpectrum.h>
#include <TVirtualFFT.h>
#include <TCanvas.h>
#include <TText.h>
#include "residualPlots.hh"


TGraphErrors* createRatioGraph(TH1D* num, TH1D* denom){

        int pointN = 0;

        TGraphErrors* ratioGraph = new TGraphErrors();

          for (int bin = 1; bin <= num->GetNbinsX(); ++bin)
          {
            if(denom->GetBinContent(bin) > 0){

              double ratio = num->GetBinContent(bin) / denom->GetBinContent(bin);
              double ratioErr = sqrt( (1 - ratio*ratio) / denom->GetBinContent(bin) ); // error: dR = sqrt( (1-R^2) / (U+V) )

              if(abs(ratio) <= 1 && ratioErr != 0){ // both U and V should be non-zero and postiive, which this takes care of (they can go negative sometimes when subtracting off pileup leading to a ratio greater than 1)

                ratioGraph->SetPoint(pointN, num->GetBinCenter(bin), ratio);
                ratioGraph->SetPointError(pointN, 0, ratioErr);

                pointN++;
              }
            }
          }

        return ratioGraph;
}

double ratioFunc(double *x, double *p){
  double w_a = blindingWa * (1 + p[1] * 1e-6);
  return p[0] * cos(w_a * x[0] + p[2]);
}

void setParameters(TH1D* histToFit, TF1* fiveParFit){
  double startingN0 = histToFit->Integral("WIDTH") / defaultLifetime;
  fiveParFit->SetParameters(startingN0, defaultLifetime, 0.37, 0, 2.1);
}

double fiveParFunc(double *x, double *p){
  double w_a = blindingWa * (1 + p[3] * 1e-6);
  return p[0] * exp(-x[0]/p[1]) * (1 + p[2] * cos(w_a * x[0] + p[4]));
}


TF1* fitAndFillHists(int sampleNum, TH1D* wiggle, TF1* fiveParamFit, double fitStart, double fitEnd, TH1F* RPull_hist, TH1F* chi2_ndf_hist, TH1F* pValues_hist, TNtuple* fitParamsTree)
{
    fiveParamFit->SetRange(fitStart, fitEnd);
    setParameters(wiggle, fiveParamFit);

    wiggle->Fit(fiveParamFit, "QR");

    if(sampleNum == 0 || saveAllFits){
      cout << endl << string(wiggle->GetName()) << endl << "--------" << endl;
      cout << "N = " << fiveParamFit->GetParameter(0) << " +- " << fiveParamFit->GetParError(0) << endl;
      cout << "tau = " << fiveParamFit->GetParameter(1) << " +- " << fiveParamFit->GetParError(1) << endl;
      cout << "A = " << fiveParamFit->GetParameter(2) <<  " +- " << fiveParamFit->GetParError(2) << endl;
      cout << "R = " << fiveParamFit->GetParameter(3) << " +- " << fiveParamFit->GetParError(3) << endl;
      cout << "phi = " << fiveParamFit->GetParameter(4) << " +- " << fiveParamFit->GetParError(4) << endl;
      cout << "p-value is: " << fiveParamFit->GetProb() << " chi2/ndf: " << fiveParamFit->GetChisquare()/fiveParamFit->GetNDF() << endl;

        wiggle->Write((string(wiggle->GetName()) + "_" + to_string(sampleNum)).c_str());

        ResidualPlots residualPlotsClass;
        residualPlotsClass.makeResidualPlots(wiggle->GetName(), wiggle, string(wiggle->GetName()) + "_residual_fitRange_" + to_string(sampleNum), fitStart, fitEnd); 
  	}

     
    RPull_hist->Fill((fiveParamFit->GetParameter(3) - 0) / fiveParamFit->GetParError(3));
    chi2_ndf_hist->Fill(fiveParamFit->GetChisquare()/fiveParamFit->GetNDF());
    pValues_hist->Fill(fiveParamFit->GetProb());

    fitParamsTree->Fill(fiveParamFit->GetParameter(0), fiveParamFit->GetParError(0),
                           fiveParamFit->GetParameter(1), fiveParamFit->GetParError(1),
                           fiveParamFit->GetParameter(2), fiveParamFit->GetParError(2),
                           fiveParamFit->GetParameter(3), fiveParamFit->GetParError(3),
                           fiveParamFit->GetParameter(4), fiveParamFit->GetParError(4));

  return fiveParamFit;
}

macro not finished and can be deleted since I went back to submitting the datasets separately

int fitAnalyzerHistsDatasets(string filePath)
{
  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  vector<string> fileVector;

  size_t found = filePath.find_last_of(".");
  string extension = filePath.substr(found);

  if(extension == ".root") fileVector.push_back(filePath);
  else if(extension == ".txt"){
    ifstream inList(filePath);
    string path;
    while(inList >> path) fileVector.push_back(path);
  }
  else{
    printf("Not passing files in correctly. Should be either a single .root file or a .txt file with multiple .root files.\n");
    return 0;
  }



  // create output file that will hold plots
  TFile* outputFile = new TFile("analyzerFitsDatasets.root","RECREATE");

  int totalFiles = fileVector.size();
  int totalSamples = 4;


  int totalSeeds = (*(TVectorD*) inputFile->Get("Seeds"))[0];

  cout << "Number of files: " << totalFiles << endl;
  cout << "Number of samples: " << totalSamples << endl;
  cout << "Number of seeds: " << totalSeeds << endl;

/////////////////////////////////////////////////////////////////////////////////////

  // define analyzer parameters

  double nick_fit_start = 30287.6;
  double nick_fit_start_EG50 = 49982.0;
  double nick_fit_end =  650064.4;

  double david_fit_start = 30287.6;
  double david_fit_start_EG50 = 49982.0;
  double david_fit_end =  649915.2;

  double aaron_fit_start = 30190.0;
  double aaron_fit_start = 49883.08;
  double aaron_fit_end =  649925.26;

  double matteo_fit_start = 30136.4;
  double matteo_fit_start_EG50 = 49829.5;
  double matteo_fit_end =  650021.0;

  double bingzhi_fit_start = 30287.6;
  double bingzhi_fit_start_EG50 = 49982.0
  double bingzhi_fit_end =  671400.0;

  double tim_fit_start = 30000.0;
  double tim_fit_start_EG50 = 49976.2;
  double tim_fit_end =  215500.0;

/////////////////////////////////////////////////////////////////////////////////////

  auto nick_dir = outputFile->mkdir("Nick");
  nick_dir->cd();

    auto nick_T_RPull = new TH1F("nick_T_RPull", "nick_T_RPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10);
    auto nick_T_chi2_ndf = new TH1F("nick_T_chi2_ndf", "nick_T_chi2_ndf; #chi^{2}/d.o.f.; Events", 50, 0.5, 1.5);
    auto nick_T_Pvalues = new TH1F("nick_T_Pvalues", "nick_T_Pvalues; P-values; Events", 50, 0, 1);

    TNtuple* nick_T_fitParams = new TNtuple("nick_T_fitParams", "nick_T_fitParams", "N:N_err:tau:tau_err:A:A_err:R:R_err:phi:phi_err");

    auto nick_R_RPull_average = new TH1F("nick_R_RPull_average", "nick_R_RPull_average; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10);
    auto nick_R_chi2_ndf_average = new TH1F("nick_R_chi2_ndf_average", "nick_R_chi2_ndf_average; #chi^{2}/d.o.f.; Events", 50, 0.5, 1.5);
    auto nick_R_Pvalues_average = new TH1F("nick_R_Pvalues_average", "nick_R_Pvalues_average; P-values; Events", 50, 0, 1);

    TNtuple* nick_R_fitParams_average = new TNtuple("nick_R_fitParams_average", "nick_R_fitParams_average", "A:A_err:R:R_err:phi:phi_err");

    auto ratio_dir = nick_dir->mkdir("RatioSeeds");
    ratio_dir->cd();

      TH1F* nick_R_RPulls[totalSamples];
      TH1F* nick_R_chi2_ndfs[totalSamples];
      TH1F* nick_R_Pvalues[totalSamples];

      TNtuple* nick_R_fitParams[totalSamples];

      for (int sampleNum = 0; sampleNum < totalSamples; ++sampleNum)
      {
        nick_R_RPulls[sampleNum] = new TH1F(Form("nick_R_RPull_%i",sampleNum), "nick_R_RPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10);
        nick_R_chi2_ndfs[sampleNum] = new TH1F(Form("nick_R_chi2_ndf_%i",sampleNum), "nick_R_chi2_ndf; #chi^{2}/d.o.f.; Events", 50, 0.5, 1.5);
        nick_R_Pvalues[sampleNum] = new TH1F(Form("nick_R_Pvalues_%i",sampleNum), "nick_R_Pvalues; P-values; Events", 50, 0, 1);

        nick_R_fitParams[sampleNum] = new TNtuple(Form("nick_R_fitParams_%i",sampleNum), "nick_R_fitParams", "A:A_err:R:R_err:phi:phi_err");
      }


  auto david_dir = outputFile->mkdir("David");
  david_dir->cd();

    auto david_T_RPull = new TH1F("david_T_RPull", "david_T_RPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10);
    auto david_T_chi2_ndf = new TH1F("david_T_chi2_ndf", "david_T_chi2_ndf; #chi^{2}/d.o.f.; Events", 50, 0.5, 1.5);
    auto david_T_Pvalues = new TH1F("david_T_Pvalues", "david_T_Pvalues; P-values; Events", 50, 0, 1);

    TNtuple* david_T_fitParams = new TNtuple("david_T_fitParams", "david_T_fitParams", "N:N_err:tau:tau_err:A:A_err:R:R_err:phi:phi_err");

    auto david_A_RPull = new TH1F("david_A_RPull", "david_A_RPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10);
    auto david_A_chi2_ndf = new TH1F("david_A_chi2_ndf", "david_A_chi2_ndf; #chi^{2}/d.o.f.; Events", 50, 0.5, 1.5);
    auto david_A_Pvalues = new TH1F("david_A_Pvalues", "david_A_Pvalues; P-values; Events", 50, 0, 1);

    TNtuple* david_A_fitParams = new TNtuple("david_A_fitParams", "david_A_fitParams", "N:N_err:tau:tau_err:A:A_err:R:R_err:phi:phi_err");


  auto aaron_dir = outputFile->mkdir("Aaron");
  aaron_dir->cd();

    auto aaron_T_RPull = new TH1F("aaron_T_RPull", "aaron_T_RPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10);
    auto aaron_T_chi2_ndf = new TH1F("aaron_T_chi2_ndf", "aaron_T_chi2_ndf; #chi^{2}/d.o.f.; Events", 50, 0.5, 1.5);
    auto aaron_T_Pvalues = new TH1F("aaron_T_Pvalues", "aaron_T_Pvalues; P-values; Events", 50, 0, 1);

    TNtuple* aaron_T_fitParams = new TNtuple("aaron_T_fitParams", "aaron_T_fitParams", "N:N_err:tau:tau_err:A:A_err:R:R_err:phi:phi_err");

    auto aaron_A_RPull = new TH1F("aaron_A_RPull", "aaron_A_RPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10);
    auto aaron_A_chi2_ndf = new TH1F("aaron_A_chi2_ndf", "aaron_A_chi2_ndf; #chi^{2}/d.o.f.; Events", 50, 0.5, 1.5);
    auto aaron_A_Pvalues = new TH1F("aaron_A_Pvalues", "aaron_A_Pvalues; P-values; Events", 50, 0, 1);

    TNtuple* aaron_A_fitParams = new TNtuple("aaron_A_fitParams", "aaron_A_fitParams", "N:N_err:tau:tau_err:A:A_err:R:R_err:phi:phi_err");


  auto matteo_dir = outputFile->mkdir("Matteo");
  matteo_dir->cd();

    auto matteo_T_RPull = new TH1F("matteo_T_RPull", "matteo_T_RPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10);
    auto matteo_T_chi2_ndf = new TH1F("matteo_T_chi2_ndf", "matteo_T_chi2_ndf; #chi^{2}/d.o.f.; Events", 50, 0.5, 1.5);
    auto matteo_T_Pvalues = new TH1F("matteo_T_Pvalues", "matteo_T_Pvalues; P-values; Events", 50, 0, 1);

    TNtuple* matteo_T_fitParams = new TNtuple("matteo_T_fitParams", "matteo_T_fitParams", "N:N_err:tau:tau_err:A:A_err:R:R_err:phi:phi_err");

    auto matteo_A_RPull = new TH1F("matteo_A_RPull", "matteo_A_RPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10);
    auto matteo_A_chi2_ndf = new TH1F("matteo_A_chi2_ndf", "matteo_A_chi2_ndf; #chi^{2}/d.o.f.; Events", 50, 0.5, 1.5);
    auto matteo_A_Pvalues = new TH1F("matteo_A_Pvalues", "matteo_A_Pvalues; P-values; Events", 50, 0, 1);

    TNtuple* matteo_A_fitParams = new TNtuple("matteo_A_fitParams", "matteo_A_fitParams", "N:N_err:tau:tau_err:A:A_err:R:R_err:phi:phi_err");


  auto bingzhi_dir = outputFile->mkdir("Bingzhi");
  bingzhi_dir->cd();

    auto bingzhi_T_RPull = new TH1F("bingzhi_T_RPull", "bingzhi_T_RPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10);
    auto bingzhi_T_chi2_ndf = new TH1F("bingzhi_T_chi2_ndf", "bingzhi_T_chi2_ndf; #chi^{2}/d.o.f.; Events", 50, 0.5, 1.5);
    auto bingzhi_T_Pvalues = new TH1F("bingzhi_T_Pvalues", "bingzhi_T_Pvalues; P-values; Events", 50, 0, 1);

    TNtuple* bingzhi_T_fitParams = new TNtuple("bingzhi_T_fitParams", "bingzhi_T_fitParams", "N:N_err:tau:tau_err:A:A_err:R:R_err:phi:phi_err");

    auto bingzhi_A_RPull = new TH1F("bingzhi_A_RPull", "bingzhi_A_RPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10);
    auto bingzhi_A_chi2_ndf = new TH1F("bingzhi_A_chi2_ndf", "bingzhi_A_chi2_ndf; #chi^{2}/d.o.f.; Events", 50, 0.5, 1.5);
    auto bingzhi_A_Pvalues = new TH1F("bingzhi_A_Pvalues", "bingzhi_A_Pvalues; P-values; Events", 50, 0, 1);

    TNtuple* bingzhi_A_fitParams = new TNtuple("bingzhi_A_fitParams", "bingzhi_A_fitParams", "N:N_err:tau:tau_err:A:A_err:R:R_err:phi:phi_err");


  auto tim_dir = outputFile->mkdir("Tim");
  tim_dir->cd();

    auto tim_Q_RPull = new TH1F("tim_Q_RPull", "tim_Q_RPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10);
    auto tim_Q_chi2_ndf = new TH1F("tim_Q_chi2_ndf", "tim_Q_chi2_ndf; #chi^{2}/d.o.f.; Events", 50, 0.5, 1.5);
    auto tim_Q_Pvalues = new TH1F("tim_Q_Pvalues", "tim_Q_Pvalues; P-values; Events", 50, 0, 1);

    TNtuple* tim_Q_fitParams = new TNtuple("tim_Q_fitParams", "tim_Q_fitParams", "N:N_err:tau:tau_err:A:A_err:R:R_err:phi:phi_err");


/////////////////////////////////////////////////////////////////////////////////////

  // setup five parameter fit function

  TF1* fiveParamFit = new TF1("fiveParamFit", fiveParFunc, 0, 1, 5);
  fiveParamFit->SetNpx(10000);
  fiveParamFit->SetLineColor(2);

  TF1* ratioFit = new TF1("ratioFit", ratioFunc, 0, 1, 3);
  ratioFit->SetNpx(10000);
  ratioFit->SetLineColor(2);

  for (int sampleNum = 0; sampleNum < totalSamples; ++sampleNum)
  {
    cout << endl << "Sample number: " << sampleNum << "/" << totalSamples << endl;

  	// nick
  	nick_dir->cd();

  	TH1D* nick_T_wiggle = (TH1D*) (inputFile->Get(Form("topDir/SampleNum%d/nick_T_Hist", sampleNum)))->Clone("nick_T_wiggle");
    fitAndFillHists(sampleNum, nick_T_wiggle, fiveParamFit, nick_fit_start, nick_fit_end, nick_T_RPull, nick_T_chi2_ndf, nick_T_Pvalues, nick_T_fitParams);
    delete nick_T_wiggle;

      // fit different random seed ratio hists

      // ratio_dir->cd();

      double nick_R_A_total = 0, nick_R_A_err_total = 0, nick_R_R_total = 0, nick_R_R_err_total = 0, nick_R_phi_total = 0, nick_R_phi_err_total = 0;

      for (int seedNum = 0; seedNum < totalSeeds; ++seedNum)
      {
        TH1D* nick_R_U_wiggle = (TH1D*) (inputFile->Get(Form("topDir/SampleNum%d/Ratio/nick_R_U_Hist_%i", sampleNum, seedNum)))->Clone();
        TH1D* nick_R_V_wiggle = (TH1D*) (inputFile->Get(Form("topDir/SampleNum%d/Ratio/nick_R_V_Hist_%i", sampleNum, seedNum)))->Clone();

        TH1D* numeratorHist = (TH1D*) nick_R_V_wiggle->Clone("numeratorHist");
        numeratorHist->Add(nick_R_U_wiggle, -1);

        TH1D* denominatorHist = (TH1D*) nick_R_V_wiggle->Clone("denominatorHist");
        denominatorHist->Add(nick_R_U_wiggle);

        TGraphErrors* nick_R_graph = createRatioGraph(numeratorHist, denominatorHist);
        nick_R_graph->SetName("nick_R_graph");

        ratioFit->SetRange(nick_fit_start, nick_fit_end);
        ratioFit->SetParameters(0.37, 0, 2.1);

        nick_R_graph->Fit(ratioFit, "QR");

        if((sampleNum == 0 || saveAllFits) && seedNum == 0){
          cout << endl << string(nick_R_graph->GetName()) << endl << "--------" << endl << "A = " << ratioFit->GetParameter(0) << endl << "R = " << ratioFit->GetParameter(1) << endl << "phi = " << ratioFit->GetParameter(2) << endl;
          cout << "p-value is: " << ratioFit->GetProb() << " chi2/ndf: " << ratioFit->GetChisquare()/ratioFit->GetNDF() << endl;
          cout << "R error: " << ratioFit->GetParError(1) << " ppm "  << endl;

          nick_R_graph->Write((string(nick_R_graph->GetName()) + "_" + to_string(sampleNum)).c_str());

          ResidualPlots residualPlotsClass_ratio;
          residualPlotsClass_ratio.makeResidualPlots(nick_R_graph->GetName(), nick_R_graph, string(nick_R_graph->GetName()) + "_residual_fitRange_" + to_string(sampleNum), nick_fit_start, nick_fit_end); 
        }

        nick_R_RPulls[sampleNum]->Fill( (ratioFit->GetParameter(1) - 0) / ratioFit->GetParError(1) );
        nick_R_chi2_ndfs[sampleNum]->Fill(ratioFit->GetChisquare()/ratioFit->GetNDF());
        nick_R_Pvalues[sampleNum]->Fill(ratioFit->GetProb());

        nick_R_fitParams[sampleNum]->Fill(ratioFit->GetParameter(0), ratioFit->GetParError(0),
                                          ratioFit->GetParameter(1), ratioFit->GetParError(1),
                                          ratioFit->GetParameter(2), ratioFit->GetParError(2));


        nick_R_A_total += ratioFit->GetParameter(0);
        nick_R_A_err_total += ratioFit->GetParError(0);
        nick_R_R_total += ratioFit->GetParameter(1);
        nick_R_R_err_total += ratioFit->GetParameter(1);
        nick_R_phi_total += ratioFit->GetParameter(2);
        nick_R_phi_err_total += ratioFit->GetParameter(2);

        delete nick_R_U_wiggle;
        delete nick_R_V_wiggle;
        delete numeratorHist;
        delete denominatorHist;
        delete nick_R_graph;

      } // end loop over ratio random seeds

      // fill seed-averaged ratio histograms and parameters

      nick_R_RPull_average->Fill(nick_R_RPulls[sampleNum]->GetMean());
      nick_R_chi2_ndf_average->Fill(nick_R_chi2_ndfs[sampleNum]->GetMean());

      double averageChi2 = nick_R_chi2_ndfs[sampleNum]->GetMean() * ratioFit->GetNDF();
      double averagePValue = TMath::Prob(averageChi2, ratioFit->GetNDF());
      nick_R_Pvalues_average->Fill(averagePValue);

      nick_R_fitParams_average->Fill(nick_R_A_total/totalSeeds, nick_R_A_err_total/totalSeeds,
                                     nick_R_R_total/totalSeeds, nick_R_R_err_total/totalSeeds,
                                     nick_R_phi_total/totalSeeds, nick_R_phi_err_total/totalSeeds);


/////////////////////////////////////////////////////////////////////////////////////

  	// david
  	david_dir->cd();

    TH1D* david_T_wiggle = (TH1D*) (inputFile->Get(Form("topDir/SampleNum%d/david_T_Hist", sampleNum)))->Clone("david_T_wiggle");
    fitAndFillHists(sampleNum, david_T_wiggle, fiveParamFit, david_fit_start, david_fit_end, david_T_RPull, david_T_chi2_ndf, david_T_Pvalues, david_T_fitParams);
    delete david_T_wiggle;

  	TH1D* david_A_wiggle = (TH1D*) (inputFile->Get(Form("topDir/SampleNum%d/david_A_Hist", sampleNum)))->Clone("david_A_wiggle");
    fitAndFillHists(sampleNum, david_A_wiggle, fiveParamFit, david_fit_start, david_fit_end, david_A_RPull, david_A_chi2_ndf, david_A_Pvalues, david_A_fitParams);
    
    delete david_A_wiggle;

/////////////////////////////////////////////////////////////////////////////////////

  	// aaron
  	aaron_dir->cd();

  	TH1D* aaron_T_wiggle = (TH1D*) (inputFile->Get(Form("topDir/SampleNum%d/aaron_T_Hist", sampleNum)))->Clone("aaron_T_wiggle");
    fitAndFillHists(sampleNum, aaron_T_wiggle, fiveParamFit, aaron_fit_start, aaron_fit_end, aaron_T_RPull, aaron_T_chi2_ndf, aaron_T_Pvalues, aaron_T_fitParams);
    delete aaron_T_wiggle;

  	TH1D* aaron_A_wiggle = (TH1D*) (inputFile->Get(Form("topDir/SampleNum%d/aaron_A_Hist", sampleNum)))->Clone("aaron_A_wiggle");
    fitAndFillHists(sampleNum, aaron_A_wiggle, fiveParamFit, aaron_fit_start, aaron_fit_end, aaron_A_RPull, aaron_A_chi2_ndf, aaron_A_Pvalues, aaron_A_fitParams);
    delete aaron_A_wiggle;

/////////////////////////////////////////////////////////////////////////////////////

  	// matteo
  	matteo_dir->cd();

  	TH1D* matteo_T_wiggle = (TH1D*) (inputFile->Get(Form("topDir/SampleNum%d/matteo_T_Hist", sampleNum)))->Clone("matteo_T_wiggle");
    fitAndFillHists(sampleNum, matteo_T_wiggle, fiveParamFit, matteo_fit_start, matteo_fit_end, matteo_T_RPull, matteo_T_chi2_ndf, matteo_T_Pvalues, matteo_T_fitParams);
    delete matteo_T_wiggle;

  	TH1D* matteo_A_wiggle = (TH1D*) (inputFile->Get(Form("topDir/SampleNum%d/matteo_A_Hist", sampleNum)))->Clone("matteo_A_wiggle");
    fitAndFillHists(sampleNum, matteo_A_wiggle, fiveParamFit, matteo_fit_start, matteo_fit_end, matteo_A_RPull, matteo_A_chi2_ndf, matteo_A_Pvalues, matteo_A_fitParams);
    delete matteo_A_wiggle;

/////////////////////////////////////////////////////////////////////////////////////

  	// bingzhi
  	bingzhi_dir->cd();

  	TH1D* bingzhi_T_wiggle = (TH1D*) (inputFile->Get(Form("topDir/SampleNum%d/bingzhi_T_Hist", sampleNum)))->Clone("bingzhi_T_wiggle");
    fitAndFillHists(sampleNum, bingzhi_T_wiggle, fiveParamFit, bingzhi_fit_start, bingzhi_fit_end, bingzhi_T_RPull, bingzhi_T_chi2_ndf, bingzhi_T_Pvalues, bingzhi_T_fitParams);
    delete bingzhi_T_wiggle;

  	TH1D* bingzhi_A_wiggle = (TH1D*) (inputFile->Get(Form("topDir/SampleNum%d/bingzhi_A_Hist", sampleNum)))->Clone("bingzhi_A_wiggle");
    fitAndFillHists(sampleNum, bingzhi_A_wiggle, fiveParamFit, bingzhi_fit_start, bingzhi_fit_end, bingzhi_A_RPull, bingzhi_A_chi2_ndf, bingzhi_A_Pvalues, bingzhi_A_fitParams);
    delete bingzhi_A_wiggle;

/////////////////////////////////////////////////////////////////////////////////////

  	// tim
  	tim_dir->cd();

  	TH1D* tim_Q_wiggle = (TH1D*) (inputFile->Get(Form("topDir/SampleNum%d/tim_Q_Hist", sampleNum)))->Clone("tim_Q_wiggle");
    fitAndFillHists(sampleNum, tim_Q_wiggle, fiveParamFit, tim_fit_start, tim_fit_end, tim_Q_RPull, tim_Q_chi2_ndf, tim_Q_Pvalues, tim_Q_fitParams);
    delete tim_Q_wiggle;


  } // end loop over samples

/////////////////////////////////////////////////////////////////////////////////////


  outputFile->Write();
  delete outputFile;

  return 1;
}

