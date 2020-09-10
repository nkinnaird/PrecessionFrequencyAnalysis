// 3-31-20: Five parameter fit header file. This fit needs to be done before the T Method fit using TmethodFit.hh.

#ifndef FIVEPARAMFIT_HH
#define FIVEPARAMFIT_HH

class FiveParamFit{

public:

  FiveParamFit(blinding::Blinders* inBlinder);
  void fiveParameterFitMethod(TObject* obj);

  void setFitRange(double fitStartTime, double fitEndTime){
   xmin = fitStartTime;
   xmax = fitEndTime;
  };

private:
  
  double xmin;
  double xmax;
  blinding::Blinders* blinder;

  double fiveParamFunc(double* x, double* par){
    double time = x[0];
    double freq = blinder->paramToFreq(par[3]) * 1e-3; // 1e-3 at end to convert to inv ns
    double lifetime = par[1] * 1000.; // convert from us to ns for fit function

    double t0 = 0; // defaultLifetime + startTimeOfFit;
    double fitVal = par[0] * exp(-time/lifetime) * (1 + par[2] * cos(freq*(time-t0) + par[4]));
    
    return fitVal;
  };

};


FiveParamFit::FiveParamFit(blinding::Blinders* inBlinder){
  blinder = inBlinder;
}

void FiveParamFit::fiveParameterFitMethod(TObject* obj){

/////////////////////////////////////////////////////////////////////////////////////
  bool isGraph = false;
  bool isHist = false;

    if(obj->InheritsFrom("TGraph")) isGraph = true;
    else if (obj->InheritsFrom("TH1")) isHist = true;

  double startingN0 = 0;

  if(isHist) startingN0 = ((TH1F*)obj)->Integral("WIDTH") / defaultLifetime;
  else if(isGraph){ 

    double sum = 0;
    double x, y;

    for (int pointI = 0; pointI < ((TGraph*)obj)->GetN(); ++pointI)
    {
      ((TGraph*)obj)->GetPoint(pointI, x, y);
      sum += y;
    }

    startingN0 = sum*(((TGraph*)obj)->GetX()[1] - ((TGraph*)obj)->GetX()[0]) / defaultLifetime;
 }
/////////////////////////////////////////////////////////////////////////////////////

  auto fiveParamFitFunc = new TF1("fiveParamFit", this, &FiveParamFit::fiveParamFunc, xmin, xmax, 5);
  fiveParamFitFunc->SetLineColor(2);
  fiveParamFitFunc->SetNpx(10000);
  fiveParamFitFunc->SetParNames(N_string.c_str(), tau_string.c_str(), A_string.c_str(), R_string.c_str(), phi_string.c_str());

  fiveParamFitFunc->SetParameters(startingN0, startingTau, startingAsymmetry, startingR, startingPhase);

/////////////////////////////////////////////////////////////////////////////////////

  // first fit to the main phase and N (seems to work best for the time being)

  fiveParamFitFunc->FixParameter(1, startingTau);
  fiveParamFitFunc->FixParameter(2, startingAsymmetry);
  fiveParamFitFunc->FixParameter(3, startingR);
  fiveParamFitFunc->SetParLimits(4, -2*pi, 2*pi);

  TFitResultPtr fitResult;

  if(isHist) fitResult = ((TH1F*)obj)->Fit(fiveParamFitFunc, "QR0");
  else if(isGraph) fitResult = ((TGraph*)obj)->Fit(fiveParamFitFunc, "QR");

  // then fit to everything else

  fiveParamFitFunc->ReleaseParameter(0);

  fiveParamFitFunc->SetParLimits(1, 64, 65);
  // fiveParamFitFunc->ReleaseParameter(1);

  // fiveParamFitFunc->SetParLimits(2, 0.3, 0.5);
  fiveParamFitFunc->ReleaseParameter(2);
  
  fiveParamFitFunc->SetParLimits(3, -100, 100);
  // fiveParamFitFunc->ReleaseParameter(3);

  if(isHist) fitResult = ((TH1F*)obj)->Fit(fiveParamFitFunc, "QRS0");
  else if(isGraph) fitResult = ((TGraph*)obj)->Fit(fiveParamFitFunc, "QRS");

/////////////////////////////////////////////////////////////////////////////////////

  int fitStatus = -1;
  fitStatus = fitResult;

  // if(fitStatus != 0)
  // {
  //   cout << "Five param fit returned non-zero fit status. Exiting with fit status: " << fitStatus << endl;
  //   exit(1);
  // }

  // fitResult->GetCorrelationMatrix().Write("CorrelationMatrix"); // write the correlation matrix to the file

}

#endif //! FIVEPARAMFIT_HH
