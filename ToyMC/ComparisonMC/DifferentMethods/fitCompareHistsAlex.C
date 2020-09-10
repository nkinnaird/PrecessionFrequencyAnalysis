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
#include <TVectorD.h>
#include <TFitResult.h>
#include <TSpectrum.h>
#include <TText.h>
#include <TLatex.h>
#include <TPaveText.h>

#include <time.h>

using namespace std;

// ******************************************************************* //                                                                           
// ********** Global variable definition and initialisation ********** //
// ******************************************************************* //

// Total number of 100 MeV enegry bins (3.1 GeV max)
int eBins = 31;

// pi
static const double pi = 3.14159265358979323846;

// Approximum maximum histogram time
static const double approxMaxTime = 700000; // ns
// Histgram bin width (cyclotron period)
static const double binWidth = 149.15;
// Total number of bins
static const int nBins = int(approxMaxTime/binWidth);
// Maximum histogram time
static const double histMaxTime = nBins*binWidth;

// Define number of fit parameters
static const int numTmethodParams = 5; // T/A-method
static const int numRatioParams = 3; // Ratio-method

// Fit start time
static const double fitStart = 30000; // ns
// Fit end time
static const double fitEnd = 350000; // ns
static const double ratioFitEnd = 350000; // ns
static const double QFitEnd = 218000; // ns

// Muon lifetime [ns]
static const double defaultLifetime = 64400; // ns                                                                                        

// Precession frequency values (blinding)
// 0.2291 MHz converted to units of 1/ns, 0.2291 is the value used in the blinding software                                             
static const double blindingFa = 0.2291 * 1e6 * 1/(1e9);
static const double blindingWa = 2*pi*blindingFa;

// Precession frequency values (no blinding)
// 0.2290735 is the average value in column 2 of Table 15 of the E821 Final Report, for the fa values for the different run periods
static const double defaultFa = 0.2290735 * 1e6 * 1/(1e9); 
static const double defaultWa = 2*pi*defaultFa;

// g-2 period
static const double g2Period = 1/defaultFa;

TH1F* RescaleAxis(TH1* input, Double_t Scale) {
  int bins = input->GetNbinsX();
  TAxis* xaxis = input->GetXaxis();
  double* ba = new double[bins+1];
  xaxis->GetLowEdge(ba);
  ba[bins] = ba[bins-1] + xaxis->GetBinWidth(bins);
  for (int i = 0; i < bins+1; i++) {
    ba[i] *= Scale;
  }
  TH1F* out = new TH1F(Form("%s_Rescaled",input->GetName()), input->GetTitle(), bins, ba);
  for (int i = 0; i <= bins; i++) {
    out->SetBinContent(i, input->GetBinContent(i));
    out->SetBinError(i, input->GetBinError(i));
  }
  return out;
}


//***********************************************************************//
//************** Gauss-Jordan Elimination Function **********************//
//************** Taken from Numerical Recipes for C++ ******************//
//**********************************************************************//

template <class T> void SWAP ( T& a, T& b )
{
  T c(a); a=b; b=c;
}

void gaussj(double** a, int n, double* b) {
  
  int i,icol=0,irow=0,j,k,l,ll;
  double big,dum,pivinv;
  int indxr[n], indxc[n], ipiv[n];
  
  for (j=0;j<n;j++) ipiv[j]=0;
  for (i=0;i<n;i++) {
    big=0.0;
    for (j=0;j<n;j++)
      if (ipiv[j] != 1)
	for (k=0;k<n;k++) {
	  if (ipiv[k] == 0) {
	    if (fabs(a[j][k]) >= big) {
	      big=fabs(a[j][k]);
	      irow=j;
	      icol=k;
	    }
	  }
	}
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l]);
      SWAP(b[irow],b[icol]);
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if (a[icol][icol] == 0.0) cout << "gaussj: Singular Matrix" << endl;
    pivinv=1.0/a[icol][icol];
    a[icol][icol]=1.0;
    for (l=0;l<n;l++) a[icol][l] *= pivinv;
    b[icol] *= pivinv;
    for (ll=0;ll<n;ll++)
      if (ll != icol) {
	dum=a[ll][icol];
	a[ll][icol]=0.0;
	for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
	b[ll] -= b[icol]*dum;
      }
  }
  for (l=n-1;l>=0;l--) {
    if (indxr[l] != indxc[l])
      for (k=0;k<n;k++)
	SWAP(a[k][indxr[l]],a[k][indxc[l]]);
  }
}


//***********************************************************************//
//************************ Fit Function Class ***************************//
//***********************************************************************//

// Five parameter function
class wiggleFitClass{

public:
  wiggleFitClass(int numTotalParameters)
  {
    numFunctionParameters = numTotalParameters;
  }

  double fiveParFunc(double* x, double* par){
    double time = x[0];
    double freq = blindingWa * (1 + par[3] * 1e-6);  // want to use frequency corresponding to blinding software since ToyMC fits are unblinded
    double lifetime = par[1] * 1000.; // convert from us to ns for fit function

    double value = par[0] * exp(-time/lifetime) * (1 + par[2] * cos(freq*time + par[4]));

    return value;
  };

  double threeParFunc(double* x, double* par){
    double time = x[0];
    double freq = blindingWa * (1 + par[1] * 1e-6);  // want to use frequency corresponding to blinding software since ToyMC fits are unblinded

    double value = par[0] * cos(freq*time + par[2]);

    return value;
  };


private:
  int numFunctionParameters;

}; // end wiggleFitClass

TH1F* SetupFFT(TH1* h, double xmin, double xmax){
  double timeTot = xmax - xmin;
  double binW = h->GetBinWidth(1);
  int nBins = timeTot / binW;
  TH1F* hout = new TH1F("","",nBins, xmin, xmin + (nBins*binW));
  int binCount = 0;
  for (int i(0); i < h->GetXaxis()->GetNbins(); i++){
    if (h->GetBinCenter(i) < xmin ) continue;
    if (binCount > nBins) break;
    binCount++;
    double cont = h->GetBinContent(i);
    double err = h->GetBinError(i);
    hout->SetBinContent(binCount, cont);
    hout->SetBinError(binCount, err);
  }
  return hout;
}

//***********************************************************************//
//*************************** Main Program ******************************//
//***********************************************************************//

int fitCompareHists(std::string filePath)
{
  //  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  // pull in input file
  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  // create output file that will hold plots
  TFile* outputFile = new TFile("histCompare.root","RECREATE");

  // make top directory for output file
  auto topDir = outputFile->mkdir("topDir");

  /////////////////////////////////////////////////////////////////////////////////////

  int analysisNum = 5; 
  
  int totalIters = (*(TVectorD*) inputFile->Get("Iters"))[0]; // total iterations in generated histograms (energy thresholds, etc.)

  //  totalIters = 1;

  TGraph* asymPlot[totalIters];
  TF1* asymFit[totalIters];
  TGraphErrors* corrGraph[analysisNum][analysisNum];
  for (int i = 0; i <analysisNum; i++){
    for (int j = 0; j < analysisNum; j++){
      corrGraph[i][j]= new TGraphErrors();
    }
  }
  // TF1s for comparison fits
  TF1* corrFit[analysisNum][analysisNum];

  auto TmethodRPull = new TH1F("TmethodRPull", "TmethodRPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10);
  auto Tmethodchi2 = new TH1F("Tmethodchi2", "Tmethodchi2; #chi^{2}/d.o.f.; Events", 50, 0.5, 1.5);
  auto TmethodPvalues = new TH1F("TmethodPvalues", "TmethodPvalues; P-values; Events", 50, 0, 1);
  auto AmethodRPull = new TH1F("AmethodRPull", "AmethodRPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10);
  auto Amethodchi2 = new TH1F("Amethodchi2", "Amethodchi2; #chi^{2}/d.o.f.; Events", 50, 0.5, 1.5);
  auto AmethodPvalues = new TH1F("AmethodPvalues", "AmethodPvalues; P-values; Events", 50, 0, 1);
  auto RmethodRPull = new TH1F("RmethodRPull", "RmethodRPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10); 
  auto Rmethodchi2 = new TH1F("Rmethodchi2", "Rmethodchi2; #chi^{2}/d.o.f.; Events", 50, 0.5, 1.5);
  auto RmethodPvalues = new TH1F("RmethodPvalues", "RmethodPvalues; P-values; Events", 50, 0, 1);
  auto QmethodRPull = new TH1F("QmethodRPull", "QmethodRPull; (R_{fit}-R_{true})/#sigma_{R_{fit}}; Events", 50, -10, 10); 
  auto Qmethodchi2 = new TH1F("Qmethodchi2", "Qmethodchi2; #chi^{2}/d.o.f.; Events", 50, 0.5, 1.5);
  auto QmethodPvalues = new TH1F("QmethodPvalues", "QmethodPvalues; P-values; Events", 50, 0, 1);

  // TH1F* TmethodResiduals[totalIters];
  // TH1F* TmethodResidualsFFT = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_5_Param_Hist_T", 0)))->Clone();;
  // TCanvas* residualPlot = new TCanvas("residuals","residuals",800,600);
  // TCanvas* residualFFT = new TCanvas("residualFFT","residualFFT",800,600);

  string analysisName[analysisNum];
  string analysisID[analysisNum];
  analysisName[0] = "T-method";
  analysisName[1] = "A-method";
  analysisName[2] = "Ratio";
  if (analysisNum > 3) analysisName[3] = "Q-method";
  if (analysisNum > 4) analysisName[4] = "E-binned";
  analysisID[0] = "T";
  analysisID[1] = "A";
  analysisID[2] = "R";
  if (analysisNum > 3) analysisID[3] = "Q";
  if (analysisNum > 4) analysisID[4] = "E";

  double eThreshold[analysisNum];
  eThreshold[0] = 1680; // T
  eThreshold[1] = 1100; // A
  eThreshold[2] = 1680; // R
  eThreshold[3] = 250;  // Q
  eThreshold[4] = 500;  // E

  double R[analysisNum][totalIters];
  double Rerr[analysisNum][totalIters];
  
  /////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////
  
  auto fiveParamFit = new wiggleFitClass(numTmethodParams);
  auto ratioFit = new wiggleFitClass(numRatioParams);

  // Initialise progress counter variables                                                                                                               
  double targetPerc = 0; // Printout for percentage complete                                                                                               
  clock_t refTime = clock();
  double tenPerTime = 0;

  double entries[analysisNum];
  // double Tentries = 0;
  // double Aentries = 0;
  // double Rentries = 0;
  // double Qentries = 0;
  // double Eentries = 0;

  // Muon asymmetry function (lab frame)                                                                                                                                     
  double Emax = 3100; // Maximum energy [MeV]                                                                                                                                
  TF1 *muonA = new TF1("muonA", "(-8.0*(x/[0])*(x/[0]) + (x/[0]) + 1.1)/(4.0*(x/[0])*(x/[0]) - 5.0*(x/[0]) - 5.0)", 0, Emax);
  muonA->SetParameter(0,Emax);

  clock_t startTime = clock();
 
  for (int iter = 0; iter < totalIters; ++iter){ // loop over iterations / different random seeds
 
    // //    cout << "\nIteration: " << iter << endl;    
  
    double RfromT = 0;
    double TmethodError = 0;
    double RfromA = 0;
    double AmethodError = 0;
    double RfromR = 0;
    double ratioError = 0;
    double RfromQ = 0;
    double QmethodError = 0;
    double RfromE = 0;
    double EbinnedError = 0;
    
    /////////////////////////////////////////////////////////////////////////////////////     
    
    //////////////////////////////////////////////////////////////////////////////////////////////
    // Fit histograms
    ////////////////////////////////////////////////////////////////////////////////////////////// 
      
    // Fit T-weighted histogram
    TH1F* hWiggleT = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_5_Param_Hist_T", iter)))->Clone();
    hWiggleT->SetName("hWiggleT");
    // TmethodResiduals[iter] = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_5_Param_Hist_T", iter)))->Clone();
    // TmethodResiduals[iter]->SetName("TmethodResiduals");
    // TH1F* TmethodResidualsFFT = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_5_Param_Hist_T", iter)))->Clone();
    // TmethodResiduals[iter]->Reset();
    // TmethodResidualsFFT->Reset();
    // TmethodResidualsFFT->SetName("TmethodResidualsFFT");
    //    Tentries = hWiggleT->GetEntries();
    entries[0] = hWiggleT->GetEntries();
    auto fitTmethod = new TF1("fitTmethod", fiveParamFit, &wiggleFitClass::fiveParFunc, fitStart, fitEnd, numTmethodParams);
    fitTmethod->SetNpx(10000);
  
    // Set fit parameters
    fitTmethod->SetParameter(0, 1);
    fitTmethod->SetParameter(1, defaultLifetime/1000);
    fitTmethod->SetParameter(2, 0.4);
    fitTmethod->SetParameter(3, 0);
    fitTmethod->SetParameter(4, pi);
      
    // Choose good starting guess normalisation N by comparing integral of function and histogram     
    double normGuessT = hWiggleT->Integral(hWiggleT->GetXaxis()->FindBin(0.0), hWiggleT->GetXaxis()->FindBin(95000), "WIDTH") / fitTmethod->Integral(0,95000);
    fitTmethod->SetParameter(0, normGuessT);  
      
    // Fit function
    hWiggleT->Fit(fitTmethod,"QNR");

    double Tchi2 = fitTmethod->GetChisquare()/fitTmethod->GetNDF();
    double TmethodPvalue = fitTmethod->GetProb();
      
    if (totalIters == 1){

      cout << "\nT-method" << endl;
      cout << "--------" << endl;
      cout << "N = " << fitTmethod->GetParameter(0) << endl;
      cout << "tau = " << fitTmethod->GetParameter(1) << endl;
      cout << "A = " << fitTmethod->GetParameter(2) << endl;
      cout << "R = " << fitTmethod->GetParameter(3) << endl;
      cout << "phi = " << fitTmethod->GetParameter(4) << endl;
      
      std::cout << "T-method fit p-value is: " << TmethodPvalue << " chi2/ndf: " << fitTmethod->GetChisquare()/fitTmethod->GetNDF() << std::endl;
      cout << "T-method R error: " << fitTmethod->GetParError(3) << " ppm "  << endl;

    }

    RfromT = fitTmethod->GetParameter(3);
    TmethodError = fitTmethod->GetParError(3);
     
    R[0][iter] = RfromT;
    Rerr[0][iter] = TmethodError;

    // Calculate T-method fit pulls
    TmethodRPull->Fill( (RfromT - 0) / TmethodError );
    Tmethodchi2->Fill(Tchi2);
    TmethodPvalues->Fill(TmethodPvalue);

    // Calculate residuals
    // residualPlot->cd();
    // int nBins = hWiggleT->GetNbinsX();
    // for (int i = 1; i <= nBins; i++){
    // 	double time = hWiggleT->GetBinCenter(i);
    // 	double data = hWiggleT->GetBinContent(i);
    // 	double fit = fitTmethod->Eval(time);
    // 	double residual = data - fit;
    // 	TmethodResiduals[iter]->SetBinContent(i,residual);
    // }
    // TmethodResiduals[iter]->SetAxisRange(fitStart, fitEnd);
    // TmethodResiduals[iter]->Draw();
    // residualFFT->cd();
    // // Do FFT                                                                                                                                                              
    // TH1 *hm2 = 0;
    // TVirtualFFT::SetTransform(0); 

    // TH1F* fftHist = SetupFFT(TmethodResiduals[iter],fitStart, fitEnd);
    // hm2 = fftHist->FFT(hm2,"MAG");
    // //      TH1F* fft2 = RescaleAxis(hm2,1./(TmethodResiduals[iter]->GetXaxis()->GetXmax() - TmethodResiduals[iter]->GetXaxis()->GetXmin()));
    // TH1F* fft2 = RescaleAxis(hm2,1./(fitEnd-fitStart));
    // fft2->GetXaxis()->SetRangeUser(0.1,35);
    // fft2->SetTitle("FFT of Residuals;Frequency (MHz);Magnitude [Arb Units]");
    // fft2->SetStats(0);
    // fft2->Draw("HIST");
    // residualFFT->Update();
    // return 1;

    ////////////////////////////////////////////////////////////////////////////////////////////// 
      
    // Fit A-weighted histogram
    TH1F* hWiggleA = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_5_Param_Hist_A", iter)))->Clone();
    //    Aentries = hWiggleA->GetEntries();
    entries[1] = hWiggleT->GetEntries();
    auto fitAmethod = new TF1("fitAmethod", fiveParamFit, &wiggleFitClass::fiveParFunc, fitStart, fitEnd, numTmethodParams);
    fitAmethod->SetNpx(10000);
    // Set fit parameters
    fitAmethod->SetParameter(0, 1);
    fitAmethod->SetParameter(1, defaultLifetime/1000);
    fitAmethod->SetParameter(2, 0.4);
    fitAmethod->SetParameter(3, 0);
    fitAmethod->SetParameter(4, pi);
      
    // Choose good starting guess normalisation N by comparing integral of function and histogram     
    double normGuessA = hWiggleA->Integral(hWiggleA->GetXaxis()->FindBin(0.0), hWiggleA->GetXaxis()->FindBin(95000), "WIDTH") / fitAmethod->Integral(0,95000);
    fitAmethod->SetParameter(0, normGuessA);  
      
    // Fit function
    hWiggleA->Fit(fitAmethod,"QNR");

    double Achi2 = fitAmethod->GetChisquare()/fitAmethod->GetNDF();
    double AmethodPvalue = fitAmethod->GetProb();
      
    if (totalIters == 1){
      
      cout << "\nA-method" << endl;
      cout << "--------" << endl;
      cout << "N = " << fitAmethod->GetParameter(0) << endl;
      cout << "tau = " << fitAmethod->GetParameter(1) << endl;
      cout << "A = " << fitAmethod->GetParameter(2) << endl;
      cout << "R = " << fitAmethod->GetParameter(3) << endl;
      cout << "phi = " << fitAmethod->GetParameter(4) << endl;
	
      std::cout << "A-method fit p-value is: " << AmethodPvalue << " chi2/ndf: " << fitAmethod->GetChisquare()/fitAmethod->GetNDF() << std::endl;
      cout << "A-method R error: " << fitAmethod->GetParError(3) << " ppm "  << endl;
	
    }
      
    RfromA = fitAmethod->GetParameter(3);
    AmethodError = fitAmethod->GetParError(3);

    R[1][iter] = RfromA;
    Rerr[1][iter] = AmethodError;

    // Calculate A-method fit pulls
    AmethodRPull->Fill( (RfromA - 0) / AmethodError );
    Amethodchi2->Fill(Achi2);
    AmethodPvalues->Fill(AmethodPvalue);

    ////////////////////////////////////////////////////////////////////////////////////////////// 
      
    // Fit ratio histogram
    TH1F* ratioNum = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_Num_Hist", iter)))->Clone();
    int nBins = ratioNum->GetNbinsX();
    TH1F* ratioDen = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_Denom_Hist", iter)))->Clone();
    TH1F* hWiggleR = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_Ratio_Hist", iter)))->Clone();
    //    Rentries = hWiggleR->GetEntries();
    entries[2] = hWiggleT->GetEntries();
    //      hWiggleR->Reset();
    // Calculate ratio histogram and errors 
    for (int i = 0; i < nBins; i++){
      double U = ratioNum->GetBinContent(i);
      double V = ratioDen->GetBinContent(i);
      double num = V-U;
      double den = V+U;
      double R = 0;
      double err = 0;
      if (den == 0){
	R = 0;
	err = 0;
      }
      else{
	R = (V-U)/(V+U);
	err = sqrt((1-pow(R,2))/(V+U));
      }
      //	hWiggleR->SetBinContent(i,R);
      //	hWiggleR->SetBinError(i,err);
    }
    auto ratioFit = new wiggleFitClass(numRatioParams);
    auto fitRatio = new TF1("fitRatio", ratioFit, &wiggleFitClass::threeParFunc, fitStart, ratioFitEnd, numRatioParams);
    fitRatio->SetNpx(10000);
      
    // Set fit parameters
    fitRatio->SetParameter(0, 0.4);
    fitRatio->SetParameter(1, 0);
    fitRatio->SetParameter(2, pi);

    // Fit function
    hWiggleR->Fit(fitRatio,"QNR");

    double ratiochi2 = fitRatio->GetChisquare()/fitRatio->GetNDF();
    double ratioPvalue = fitRatio->GetProb();


    if (totalIters == 1){
	
      cout << "\nRatio-method" << endl;
      cout << "--------" << endl;
      cout << "A = " << fitRatio->GetParameter(0) << endl;
      cout << "R = " << fitRatio->GetParameter(1) << endl;
      cout << "phi = " << fitRatio->GetParameter(2) << endl;
	
      double ratioPvalue = fitRatio->GetProb();
	
      std::cout << "Ratio-method fit p-value is: " << ratioPvalue << " chi2/ndf: " << fitRatio->GetChisquare()/fitRatio->GetNDF() << std::endl;
      cout << "Ratio-method R error: " << fitRatio->GetParError(1) << " ppm "  << endl;

    }

    RfromR = fitRatio->GetParameter(1);
    ratioError = fitRatio->GetParError(1);

    R[2][iter] = RfromR;
    Rerr[2][iter] = ratioError;

    // Calculate ratio method fit pulls
    RmethodRPull->Fill( (RfromR - 0) / ratioError );
    Rmethodchi2->Fill(ratiochi2);
    RmethodPvalues->Fill(ratioPvalue);

    if (analysisNum > 3){

      ////////////////////////////////////////////////////////////////////////////////////////////// 

      // Fit Q-method histogram
      TH1F* hWiggleQ = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/Toy_5_Param_Hist_Q", iter)))->Clone();
      hWiggleQ->SetName("hWiggleQ");
      entries[3] = hWiggleT->GetEntries();
      //      Qentries = hWiggleQ->GetEntries();
      auto fitQmethod = new TF1("fitQmethod", fiveParamFit, &wiggleFitClass::fiveParFunc, fitStart, QFitEnd, numTmethodParams);
      fitQmethod->SetNpx(10000);
  
      // Set fit parameters
      fitQmethod->SetParameter(0, 1);
      fitQmethod->SetParameter(1, defaultLifetime/1000);
      fitQmethod->SetParameter(2, 0.4);
      fitQmethod->SetParameter(3, 0);
      fitQmethod->SetParameter(4, pi);
      
      // Choose good starting guess normalisation N by comparing integral of function and histogram     
      double normGuessQ = hWiggleQ->Integral(hWiggleQ->GetXaxis()->FindBin(0.0), hWiggleQ->GetXaxis()->FindBin(95000), "WIDTH") / fitQmethod->Integral(0,95000);
      fitQmethod->SetParameter(0, normGuessQ);  
      
      // Fit function
      hWiggleQ->Fit(fitQmethod,"QNR");

      double Qchi2 = fitQmethod->GetChisquare()/fitQmethod->GetNDF();
      double QmethodPvalue = fitQmethod->GetProb();
      
      if (totalIters == 1){

	cout << "\nQ-method" << endl;
	cout << "--------" << endl;
	cout << "N = " << fitQmethod->GetParameter(0) << endl;
	cout << "tau = " << fitQmethod->GetParameter(1) << endl;
	cout << "A = " << fitQmethod->GetParameter(2) << endl;
	cout << "R = " << fitQmethod->GetParameter(3) << endl;
	cout << "phi = " << fitQmethod->GetParameter(4) << endl;
      
	std::cout << "Q-method fit p-value is: " << QmethodPvalue << " chi2/ndf: " << fitQmethod->GetChisquare()/fitQmethod->GetNDF() << std::endl;
	cout << "Q-method R error: " << fitQmethod->GetParError(3) << " ppm "  << endl;

      }

      RfromQ = fitQmethod->GetParameter(3);
      QmethodError = fitQmethod->GetParError(3);
     
      R[3][iter] = RfromQ;
      Rerr[3][iter] = QmethodError;

      // Calculate Q-method fit pulls
      QmethodRPull->Fill( (RfromQ - 0) / QmethodError );
      Qmethodchi2->Fill(Qchi2);
      QmethodPvalues->Fill(QmethodPvalue);

    }

    /////////////////////////////////////////////////////////////////////////////////////     

    if (analysisNum > 4){

      // Fit E-binned histograms and average
      double avg = 0;
      double avgErr = 0;
      double avgNum = 0;
      double avgDen = 0;
      // asymPlot[iter] = new TGraph(); 
      entries[4] = 0;
      for (int bin = 0; bin < eBins; ++bin){

	double eLow = double(bin*100);
	double eHigh = double((bin+1)*100); 
	double energy = (eLow+eHigh)/2.0;
	double A = muonA->Eval(energy);
      
	TH1F* toyEbin = (TH1F*) ((TH1F*) inputFile->Get(Form("topDir/Iter%d/%d - %d MeV",iter,bin*100,(bin+1)*100)))->Clone();
	entries[4] += toyEbin->GetEntries();
	auto eBinFit = new TF1("eBin", fiveParamFit, &wiggleFitClass::fiveParFunc, fitStart, fitEnd, numTmethodParams);
	eBinFit->SetNpx(10000);
      
	// Set fit parameters
	eBinFit->SetParameter(0, 1);
	eBinFit->SetParameter(1, defaultLifetime/1000);
	eBinFit->SetParameter(2, A);
	eBinFit->SetParameter(3, 0);
	eBinFit->SetParameter(4, pi);
      
	// Choose good starting guess normalisation N by comparing integral of function and histogram     
	double normGuess = toyEbin->Integral(toyEbin->GetXaxis()->FindBin(0.0), toyEbin->GetXaxis()->FindBin(95000), "WIDTH") / eBinFit->Integral(0,95000);
	eBinFit->SetParameter(0, normGuess);  
      
	// Fit function
	toyEbin->Fit(eBinFit,"QNREM");

	// Sum towards E-binned average
	double par = eBinFit->GetParameter(3);
	double parErr = eBinFit->GetParError(3);
	// Set E-binned threshold
	if (eLow > eThreshold[4] && eLow < 2900) {
	  avgNum += par/pow(parErr,2);
	  avgDen += 1/pow(parErr,2);
	}
     
	// double x = (bin*100 + (bin+1)*100)/2.0;
	// double y = eBinFit->GetParameter(2);
	//      asymPlot[iter]->SetPoint(bin,x,y);
      
      }
   
      // Calculate E-binned average
      avg = avgNum/avgDen;
      avgErr = 1/sqrt(avgDen);
      RfromE = avg;
      EbinnedError = avgErr;
      R[4][iter] = RfromE;
      Rerr[4][iter] = EbinnedError;

      if (totalIters == 1){

	cout << "\nE-binned" << endl;
	cout << "--------" << endl;
	cout << "R = " << RfromE << endl;
	cout << "T-binned R error: " << EbinnedError << " ppm "  << endl;

      }

      // asymFit[iter] = new TF1(Form("asymFit_Iter%d",iter),"pol5");
      // asymPlot[iter]->Fit(asymFit[iter],"Q");
   
    }

    /////////////////////////////////////////////////////////////////////////////////////

    for (int i = 0; i < analysisNum; i++){
      for (int j = 0; j < analysisNum; j++){
	corrGraph[i][j]->SetPoint(iter,R[i][iter],R[j][iter]);
	corrGraph[i][j]->SetPointError(iter,Rerr[i][iter],Rerr[j][iter]);
      }
    }

    // Print the percentage of entries that have been processd and the time remaining                                                                                     
    if(100*float(iter) / totalIters > targetPerc){
      // Calculate percentage of events processed                                                                                                                          
      double percent = int(100*float(iter) / totalIters);
      // Determine how long it took to process 10% of the events                                                                                                          
      if (int(percent) == 10) {
	refTime = clock();
	tenPerTime = double(refTime - startTime)/CLOCKS_PER_SEC;
      }
      // Output the percentage completed for each 10% and time reminaing                                                                                                  
      if (int(percent) % 10 == 0) {
	if (int(percent) >= 10){
	  int seconds = int((100-percent)/10*tenPerTime);
	  int minutes = seconds / 60;
	  int hours = minutes / 60;
	  if (int(hours) > 0 && int(minutes%60) > 0) {
	    cout << "Processed " << int(percent)  << "%; Estimated time remaining = " << int(hours) << " hours, " << int(minutes%60) << " minutes & " << int(seconds%60) << " seconds." << endl;
	  }
	  else if (int(minutes%60) > 0) {
	    cout << "Processed " << int(percent)  << "%; Estimated time remaining = " << int(minutes%60) << " minutes & " << int(seconds%60) << " seconds." << endl;
	  }
	  else {
	    cout << "Processed " << int(percent)  << "%; Estimated time remaining = " << int(seconds%60) << " seconds." << endl;
	  }
	}
	else{
	  cout << "Processed " << int(percent)  << "%" << endl;
	}
      }
      // Update percentage process by 1                                                                                                                                   
      targetPerc += 1;
    }
      
  } // end loop over iterations

  if (totalIters > 1){

    TCanvas* corrPlot[analysisNum][analysisNum];    
    for (int i = 0; i < analysisNum; i++){
      for (int j = 0; j < analysisNum; j++){
	// Print correlation coefficients
	cout << Form("%s vs. %s: correlation = ",analysisName[i].c_str(),analysisName[j].c_str()) << corrGraph[i][j]->GetCorrelationFactor() << endl;
	if (i != j && j > i){
	  corrPlot[i][j] = new TCanvas(Form("%svs%splot",analysisID[i].c_str(),analysisID[j].c_str()),
				       Form("%svs%splot",analysisID[i].c_str(),analysisID[j].c_str()),800,800);
	  corrPlot[i][j]->cd();
	  corrGraph[i][j]->Draw("Ap*E");
	  corrGraph[i][j]->SetTitle(Form("%s vs. %s",analysisName[i].c_str(),analysisName[j].c_str()));
	  corrGraph[i][j]->GetXaxis()->SetTitle(Form("%s: R [ppm]",analysisName[i].c_str()));
	  corrGraph[i][j]->GetYaxis()->SetTitle(Form("%s: R [ppm]",analysisName[j].c_str()));
	  TPaveText *text = new TPaveText(0.15,0.7,0.48,0.85,"NDC");
	  text->SetTextSize(0.022);
	  text->AddText(Form("No. of samples = %d" ,totalIters));
	  text->AddText(Form("%s threshold = %.2f GeV",analysisName[i].c_str(),eThreshold[i]/1000));
	  if (i == 0 or i == 1 or i == 4){ // If T, A or E
	    text->AddText(Form("%s events ~ %e",analysisName[i].c_str(),entries[i]));
	  }
	  else if (i == 3){
	    text->AddText(Form("%s entries ~ %e",analysisName[i].c_str(),entries[i]));
	  }
	  text->AddText(Form("%s threshold = %.2f GeV",analysisName[j].c_str(),eThreshold[j]/1000));
	  if (j == 0 or j == 1 or j == 4){ // If T, A or E
	    text->AddText(Form("%s events ~ %e",analysisName[j].c_str(),entries[j]));
	  }
	  else if (j == 3){
	    text->AddText(Form("%s entries ~ %e",analysisName[j].c_str(),entries[j]));
	  }
	  text->AddText(Form("Correlation = %2.2f %%",corrGraph[i][j]->GetCorrelationFactor()*100));
	  text->SetFillColor(0);
	  text->Draw("SAME");
	  corrFit[i][j] = new TF1(Form("%svs%s",analysisName[i].c_str(),analysisName[j].c_str()),"pol1");
	  double xmin = TMath::MinElement(corrGraph[i][j]->GetN(),corrGraph[i][j]->GetX());
	  double xmax = TMath::MaxElement(corrGraph[i][j]->GetN(),corrGraph[i][j]->GetX());
	  corrFit[i][j]->SetRange(xmin,xmax);
	  corrGraph[i][j]->Fit(corrFit[i][j],"QR");
	  corrFit[i][j]->SetLineColor(kBlue);
	  corrFit[i][j]->DrawF1(xmin,xmax,"SAME");				       
	  TImage *img = TImage::Create();
	  img->FromPad(corrPlot[i][j]);
	  img->WriteImage(Form("Plots/%svs%s.png",analysisID[i].c_str(),analysisID[j].c_str()));
	  delete img;  
	}
      }
    }
    
    // Plot master correlation matrix                                                                                               
    double matrixBins[analysisNum];
    for (int j = 0; j < analysisNum; j++){
      matrixBins[j] = double(j+1);
    }
    TCanvas *corrMatPlot = new TCanvas("corrMatPlot","Statistical correlation matrix of R",800,800);
    corrMatPlot->SetLeftMargin(0.15);
    corrMatPlot->SetBottomMargin(0.15);
    corrMatPlot->SetGrid();
    TH2D *corrMat = new TH2D("corrMat","Statistical correlation matrix of only R parameters",analysisNum,0,analysisNum,analysisNum,0,analysisNum);
    corrMat->SetStats(0);
    double corr[analysisNum][analysisNum];
    for (int j = 1; j <= analysisNum; j++) {
      corrMat->GetXaxis()->SetBinLabel(j,Form("%s",analysisName[j-1].c_str()));
      corrMat->GetYaxis()->SetBinLabel(j,Form("%s",analysisName[j-1].c_str()));
      for (int k = 1; k <= analysisNum; k++) {
	corr[j-1][k-1] = corrGraph[j-1][k-1]->GetCorrelationFactor();
	corrMat->SetBinContent(j,k,corr[j-1][k-1]);
      }
    }
    gStyle->SetPalette(55);
    gStyle->SetPaintTextFormat("1.2f");
    corrMat->GetZaxis()->SetRangeUser(-1.0,1.0);
    corrMatPlot->cd();
    corrMat->Draw("TEXT COLZ");
    
    TImage *img = TImage::Create();
    img->FromPad(corrMatPlot);
    img->WriteImage("Plots/corrMatPlot.png");
    delete img;  

    // Calculate analysis weighted average values and correlated errors
    double analysisMean[analysisNum];
    double meanWeight[analysisNum];
    double analysisMeanErr[analysisNum];
    double analysisAvg[analysisNum];
    double errWeight[analysisNum];
    double analysisAvgErr[analysisNum];
    cout << "\nAveraging for the same analysis type..." << endl;
    for (int i = 0 ; i < analysisNum; i++){
      cout << "\n" << analysisName[i] << endl;
      analysisMean[i] = 0;
      meanWeight[i] = 0;
      analysisMeanErr[i] = 0;
      analysisAvgErr[i] = 0;
      analysisAvg[i] = 0;
      errWeight[i] = 0;
      analysisAvgErr[i] = 0;
      double analysisAvgNum = 0;
      double analysisAvgDen = 0;
      double meanNum = 0;
      double meanDen = 0;
      for (int j = 0; j < totalIters; j++){
	meanNum += R[i][j];
	meanDen += 1;
	analysisAvgNum += R[i][j]/pow(Rerr[i][j],2);
	analysisAvgDen += 1/pow(Rerr[i][j],2);
	errWeight[i] += 1/pow(Rerr[i][j],2);
      }
      analysisMean[i] = meanNum/meanDen;
      meanWeight[i] = 1/meanDen;
      analysisAvg[i] = analysisAvgNum/analysisAvgDen;
      errWeight[i] = 1.0/errWeight[i];
      // Calculate correlated errors
      for (int j = 0; j < totalIters; j++){
	for (int k = 0; k < totalIters; k++){
	  double varj = Rerr[i][j]*Rerr[i][j];
	  double vark = Rerr[i][k]*Rerr[i][k];
	  double covjk = Rerr[i][j]*Rerr[i][k];
	  analysisMeanErr[i] += meanWeight[i]*covjk*meanWeight[i];
	  analysisAvgErr[i] += errWeight[i]*1/(varj)*covjk*1/(vark)*errWeight[i];
	}
      }
      analysisMeanErr[i] = sqrt(analysisMeanErr[i]);
      analysisAvgErr[i] = sqrt(analysisAvgErr[i]);
      cout << "--" << endl;
      cout << "    Mean average = " << analysisMean[i] << " ± " << analysisMeanErr[i] << endl;
      cout << "Weighted average = " << analysisAvg[i] << " ± " << analysisAvgErr[i] << endl;
      cout << "--" << endl;
    }

    cout << "\nAveraging the different analyses..." << endl;
    double Rmean = 0;
    double RmeanWeight = 0;
    double RmeanErr = 0;
    double Ravg = 0;
    double RerrWeight = 0;
    double RavgErr = 0;
    double analysisAvgNum = 0;
    double analysisAvgDen = 0;
    double meanNum = 0;
    double meanDen = 0;
    for (int j = 0; j < analysisNum; j++){
      meanNum += analysisAvg[j];
      meanDen += 1;
      analysisAvgNum += analysisAvg[j]/pow(analysisAvgErr[j],2);
      analysisAvgDen += 1/pow(analysisAvgErr[j],2);
      RerrWeight += 1/pow(analysisAvgErr[j],2);
    }
    Rmean = meanNum/meanDen;
    RmeanWeight = 1/meanDen;
    Ravg = analysisAvgNum/analysisAvgDen;
    RerrWeight = 1.0/RerrWeight;
    // Calculate correlated errors
    for (int j = 0; j < analysisNum; j++){
      for (int k = 0; k < analysisNum; k++){
	double varj = analysisAvgErr[j]*analysisAvgErr[j];
	double vark = analysisAvgErr[k]*analysisAvgErr[k];
	double covjk = analysisAvgErr[k]*corr[j][k]*analysisAvgErr[k];
	RmeanErr += RmeanWeight*covjk*RmeanWeight;
	RavgErr += RerrWeight*1/(varj)*covjk*1/(vark)*RerrWeight;
      }
    }
    RmeanErr = sqrt(RmeanErr);
    RavgErr = sqrt(RavgErr);

    // Turn to chi^2 minimisation
    // Define variables for minimisation
    double chisq = 0.0;
    double chisq_dof = 0.0;
    TMatrixD covMat(analysisNum,analysisNum);
    TMatrixD invCovMat(analysisNum,analysisNum);
    double Rfit = 0;
    double RfitErr = 0;
    
     // Determine covariance matrix
    for (int j = 0; j < analysisNum; j++){
      double errj = analysisAvgErr[j];
      for (int k = 0; k < analysisNum; k++){
     	double errk = analysisAvgErr[k];
     	double corrjk = corr[j][k];
     	covMat(j,k) = errj*corrjk*errk;
      }
    }
    
    // Define inverse covariance matrix as TMatrixD for easy inversion            
    for (int j = 0; j < analysisNum; j++){
      for (int k = 0; k < analysisNum; k++){
     	invCovMat(j,k) = covMat(j,k);
      }
    }
    // Invert master covariance matrix  
    invCovMat.Invert();  
    
    int fitNum = 1;
    double* gaussVec = new double [fitNum];
    double** gaussMat = new double* [fitNum];
    for (int i = 0; i < fitNum; i++) gaussMat[i] = new double [fitNum];    
    // Re-initialise variables for minimisation
    for (int j = 0; j < fitNum; j++){
      gaussVec[j] = 0.0 ;
      for (int k = 0; k < fitNum; k++){
	gaussMat[j][k] = 0.0;
      }
    }
    // Populate vector and matrix for minimisation
    for (int j = 0; j < analysisNum; j++){
      for (int k = 0; k < analysisNum; k++){
	//
	gaussVec[0] += analysisAvg[k]*invCovMat(j,k);
	gaussMat[0][0] += invCovMat(j,k);
	//
      }
    }
    
    // Perform minimisation via Gauss-Jordan Elimination
    gaussj(gaussMat, fitNum, gaussVec);
    
    // Calculate chi^2 minimum
    chisq = 0.0;
    for (int j = 0; j < analysisNum; j++){
      for (int k = 0; k < analysisNum; k++){
	chisq += (analysisAvg[j]-gaussVec[0])*invCovMat(j,k)*(analysisAvg[k]-gaussVec[0]);
      }
    }
    // Find no. of degrees of freedom                          
    double nDoF = analysisNum - 1;
    // Find global chi^2 per degree of freedom       
    chisq_dof = chisq/nDoF;
    
    Rfit = gaussVec[0];
    RfitErr = sqrt(gaussMat[0][0]);
    
    // // if (itNum >= 50){
    // //   cout << "---------------------!!ERROR!!---------------------" << endl;
    // //   cout <<  "chi^2 minimisation does not converge!" << endl;
    // //   cout <<  "Exiting..." << endl;
    // //   cout << "---------------------------------------------------" << endl;
    // //   std::exit(0);
    // // }
    
     // Inflate error by chi^2/d.o.f.               
     if (chisq_dof > 1.0) RfitErr *= chisq_dof;
    
     cout.precision(6); 
    
     // Print fit results
     cout << "--" << endl;
     cout << "    Mean average = " << Rmean << " ± " << RmeanErr << endl;
     cout << "Weighted average = " << Ravg << " ± " << RavgErr << endl;
     cout << "             Fit = " << Rfit << " ± " << RfitErr << endl;
     cout << "--" << endl;


     // Print global chi^2 of the fit                                                                                                                          
     cout << "\n_________________ Total chi^2/d.o.f.: ____________________" << endl;
     cout << "chi^2 = " << chisq << endl;
     cout << "chi^2/d.o.f. = " << chisq_dof << endl;
     cout << "------------------------------------------------------------" << endl;
    
  }

/////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////

  outputFile->Write();
  delete outputFile;

  return 1;

}
