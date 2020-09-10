// 3-31-20: Three parameter ratio fit header file. This fit needs to be done before the full ratio fit using fullRatioFit.hh.

#ifndef THREEPARAMETERRATIOFIT_HH
#define THREEPARAMETERRATIOFIT_HH

/////////////////////////////////////////////////////////////////////////////////////

class ThreeParamRatioFuncClass {

public:
  ThreeParamRatioFuncClass(double xmin, double xmax, blinding::Blinders* inBlinder, double inputBinWidth)
  {
    threeParBackEndFunc = new TF1("threeParBackEndFunc", this, &ThreeParamRatioFuncClass::threeParFunc, xmin, xmax, 3); // a second function is used to that the main fit function can be an integral or not
    threeParBackEndFunc->SetNpx(10000);

    blinder = inBlinder;
    binWidth = inputBinWidth;
  }

  double threeParFunc(double* x, double* par){
    double time = x[0];
    double freq = blinder->paramToFreq(par[1]) * 1e-3; // 1e-3 at end to convert to inv ns

    double t0 = 0; // defaultLifetime + startTimeOfFit;

    double fitVal = par[0] * cos(freq*(time-t0) + par[2]);
    return fitVal;
  }

  double Evaluate(double *x, double *p) {
    // threeParBackEndFunc->SetParameters(p);
    // return threeParBackEndFunc->Integral( x[0] - 0.5*binWidth, x[0] + 0.5*binWidth) / binWidth; // comment these in to turn on the integral functionality
    return threeParBackEndFunc->EvalPar(x, p);

  }

private:
  TF1* threeParBackEndFunc;
  blinding::Blinders* blinder;
  double binWidth;

}; // end ThreeParamRatioFuncClass class

/////////////////////////////////////////////////////////////////////////////////////


class ThreeParameterRatioFit{

public:

  ThreeParameterRatioFit(blinding::Blinders* inBlinder);
  TGraphErrors* createRatioGraph(TH1F* num, TH1F* denom);
  void fitMethod(TGraphErrors* ratioGraph);
  void denominatorFitMethod(TH1F* denom, double binWidth);

  void setBinWidth(double inputBinWidth){
    binWidth = inputBinWidth;
  };

  void setFitRange(double fitStartTime, double fitEndTime){
   xmin = fitStartTime;
   xmax = fitEndTime;   
  };

private :

  double xmin;
  double xmax;
  blinding::Blinders* blinder;
  double binWidth;

}; // end ThreeParameterRatioFit class


ThreeParameterRatioFit::ThreeParameterRatioFit(blinding::Blinders* inBlinder){
  blinder = inBlinder;
}

TGraphErrors* ThreeParameterRatioFit::createRatioGraph(TH1F* num, TH1F* denom){

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

void ThreeParameterRatioFit::fitMethod(TGraphErrors* ratioGraph){
        
        auto classForFunc = new ThreeParamRatioFuncClass(xmin, xmax, blinder, binWidth);
        auto threeParamRatioFit = new TF1("threeParamRatioFit", classForFunc, &ThreeParamRatioFuncClass::Evaluate, xmin, xmax, 3);
        threeParamRatioFit->SetNpx(10000);
        threeParamRatioFit->SetLineColor(2);
        threeParamRatioFit->SetParNames(A_string.c_str(), R_string.c_str(), phi_string.c_str());

        threeParamRatioFit->SetParameters(startingAsymmetry, startingR, startingPhase);

/////////////////////////////////////////////////////////////////////////////////////
        
        // try fitting phase first and then A and R

        threeParamRatioFit->FixParameter(0, startingAsymmetry);
        threeParamRatioFit->FixParameter(1, startingR);
        threeParamRatioFit->SetParLimits(2, -2*pi, 2*pi);
        // threeParamRatioFit->ReleaseParameter(2);

        TFitResultPtr fitResult;
        fitResult = ratioGraph->Fit(threeParamRatioFit, "QR");

        // threeParamRatioFit->SetParLimits(0, 0.3, 0.5);
        threeParamRatioFit->SetParLimits(1, -100, +100);
        threeParamRatioFit->ReleaseParameter(0);
        // threeParamRatioFit->ReleaseParameter(1);

        fitResult = ratioGraph->Fit(threeParamRatioFit, "QRS");


/////////////////////////////////////////////////////////////////////////////////////        

        int fitStatus = -1;
        fitStatus = fitResult;

        if(fitStatus != 0)
        {
          cout << "3 param Ratio fit returned non-zero fit status. Exiting with fit status: " << fitStatus << endl;
          exit(1);
        }

        // fitResult->GetCorrelationMatrix().Write("CorrelationMatrix"); // write the correlation matrix to the file
     
}

void ThreeParameterRatioFit::denominatorFitMethod(TH1F* denom, double binWidth){

  TF1* expFunc = new TF1("expFunc", "[0]*exp(-x/[1])", xmin, xmax);

  int fitStartBin = denom->FindBin(xmin);
  int fitEndBin = denom->FindBin(xmax);

  expFunc->SetParameter(0, 1);
  expFunc->SetParameter(1, defaultLifetime);

  double par0 = denom->Integral(fitStartBin, fitEndBin) / ( expFunc->Integral(xmin, xmax) / binWidth );

  expFunc->SetParameter(0, par0); 

  denom->Fit(expFunc, "QRLI");
  denom->GetFunction("expFunc")->SetLineColor(2);

}

#endif //! THREEPARAMETERRATIOFIT_HH
