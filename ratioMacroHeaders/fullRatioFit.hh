// 3-31-20: Full ratio fit header file. Do fits using threeParameterRatioFit.hh as input to this.

#ifndef FULLRATIOFIT_HH
#define FULLRATIOFIT_HH

class FullRatioFit{

public:

  FullRatioFit(blinding::Blinders* inBlinder);
  
  TF1* prepareFitFunction();
  void setStartingFitParameters(TF1* fullRatioFitFunc, double amplitudeMultiplier);

  void clearUpdatedParLimits(){
    updatedParLimits.clear();
  };

  double* doFullRatioFit(TF1* fullRatioFitFunc, TGraphErrors* ratioGraph);
  double* fitAllParameters(TF1* fullRatioFitFunc, TGraphErrors* ratioGraph, double* pars);

  void setBinWidth(double inputBinWidth){
    binWidth = inputBinWidth;
  };

  void setHalfPeriod(double halfPerGuess){
    halfPeriodGuess = halfPerGuess;
  };

  void setFitRange(double fitStartTime, double fitEndTime){
    fitStart = fitStartTime;
    fitEnd = fitEndTime;
  };

  void setThreeParamRatioFunction(TF1* ratioFit){
    threeParamRatioFit = ratioFit;
  };

  void setTMethodFunction(TF1* TmethodFit, map<string, int> TmethodParamMap){
    TmethodFitFunction = TmethodFit;
    TmethodParameterIndices = TmethodParamMap;
  };

  void setFullRatioFunction(TF1* fullRatioFunc){
    input_FullRatioFitFunction = fullRatioFunc;
  };

  void setupFitParams(bool changingCBOFrequency, bool changingVWFrequency, bool N2_cbo, bool A_cbo, bool Phi_cbo, bool incVW, bool incLM){
    changingCBOFreq = changingCBOFrequency;
    changingVWFreq = changingVWFrequency;
    include2CBO = N2_cbo;
    includeACBO = A_cbo;
    includePhiCBO = Phi_cbo;
    includeVW = incVW;
    includeLM = incLM;

    parameterNameList.insert(parameterNameList.begin(), {A_string, R_string, phi_string, cbo_w_string, cbo_tau_string, cbo_N_amp_string, cbo_N_phi_string});
    if(include2CBO) parameterNameList.insert(parameterNameList.end(), {cbo_N2_amp_string, cbo_N2_phi_string});
    if(includeACBO) parameterNameList.insert(parameterNameList.end(), {cbo_A_amp_string, cbo_A_phi_string});
    if(includePhiCBO) parameterNameList.insert(parameterNameList.end(), {cbo_Phi_amp_string, cbo_Phi_phi_string});
    if(includeVW) parameterNameList.insert(parameterNameList.end(), {vw_w_string, vw_tau_string, vw_amp_string, vw_phi_string});
    // if(includeVW) parameterNameList.insert(parameterNameList.end(), {vw_w_string, vw_amp_string, vw_phi_string});
    if(includeLM) parameterNameList.push_back(kLoss_string);

    // parameterNameList.insert(parameterNameList.end(), "C");
    // parameterNameList.insert(parameterNameList.end(), "C_{T}");
    // parameterNameList.insert(parameterNameList.end(), "C_{#phi}");

    for (uint parNum = 0; parNum < parameterNameList.size(); ++parNum) parameterIndices[parameterNameList.at(parNum)] = parNum; // set indexes for each parameter

    numFitParams = int(parameterIndices.size());

  };

  void setFixedParameters(vector<string> inputFixedParameters){
    fitFixedParameters = inputFixedParameters;
  };

  bool isParamFree(string paramString){
    if(find(fitFixedParameters.begin(), fitFixedParameters.end(), paramString) != fitFixedParameters.end()) return false;
    else return true;
  };

  void setLostMuonIntegral(TH1F* integralHist){
    lostMuonIntegral = integralHist;
  };

  int getNumFitParams(){
    return numFitParams;
  };

  void setScanParameterChange(double inputChange){
    scannedParameterValueChange = inputChange;
  };

  void printParamsAndLimits(TF1* fullRatioFitFunc){
      for (int parNum = 0; parNum < fullRatioFitFunc->GetNpar(); ++parNum){
        double low, high = 0;
        double value = fullRatioFitFunc->GetParameter(parNum);
        fullRatioFitFunc->GetParLimits(parNum, low, high);
        cout << "Parameter: " << parNum << " name: " << fullRatioFitFunc->GetParName(parNum) << " value: " << value << " limits: " << low << " - " << high;
          if((low != 0 && high != 0) && (value > high || value < low)) cout << " OUTSIDE LIMIT" << endl;
          else cout << endl;
      }
  };

  void fixAllParameters(TF1* fullRatioFitFunc){
    for (int parNum = 0; parNum < fullRatioFitFunc->GetNpar(); ++parNum) fullRatioFitFunc->FixParameter(parNum, fullRatioFitFunc->GetParameter(parNum));
  };

  void fillDiffParLimits(string parName, double low, double high){
    updatedParLimits.emplace_back(parName, make_pair(low, high));
  };
  
  // update limits dynamically for individual fits using passed in vector which checks if those parameters are included in the updatedParLimits vector
  void setDiffParLimits(TF1* fullRatioFitFunc, vector<string> parsToUpdate){

    for (uint index = 0; index < parsToUpdate.size(); ++index)
    {
      auto it = std::find_if( updatedParLimits.begin(), updatedParLimits.end(),
      [&parsToUpdate, index](const pair<string, pair<double, double> >& element){ return parsToUpdate.at(index) == element.first; } );

      if(it != updatedParLimits.end() && isParamFree(parsToUpdate.at(index))) fullRatioFitFunc->SetParLimits(parameterIndices[parsToUpdate.at(index)], (*it).second.first, (*it).second.second);
    }
  };

private:
  
  blinding::Blinders* blinder;
  double binWidth = 0;
  double fitStart, fitEnd;

  vector<string> parameterNameList; // initialized in setupFitParams
  std::map<string, int> parameterIndices; // map of strings to parameter indices for more flexible fitting with respect to additional parameters
  int numFitParams = 0;

  TF1* threeParamRatioFit;
  TF1* TmethodFitFunction = 0;
  map<string, int> TmethodParameterIndices;
  TF1* input_FullRatioFitFunction = 0;


  vector<string> fitFixedParameters = ratioFit_added_fixedParameters; // default initial fixed parameters to R method added fixed parameters

  double halfPeriodGuess;

  bool changingCBOFreq = false;
  bool changingVWFreq = false;

  bool include2CBO = false;
  bool includeACBO = false;
  bool includePhiCBO = false;

  bool includeVW = false;

  bool includeLM = false;
  TH1F* lostMuonIntegral;

  double scannedParameterValueChange = 0;

  vector< pair<string, pair<double, double> > > updatedParLimits;

/////////////////////////////////////////////////////////////////////////////////////

 // adjust errors to Pearsons - doesn't currently work, needs a fit function to the denominator which does not currently exist
 // regarding fitting the denominator, a similar function to the full ratio fit could be used, but in my experience including parameters which don't exist in the data in a fit (R value for the denominator), causes issues
 // in that case I think it would probably be better to define a new fit function entirely, drawing from the other bits where possible
/*
  void setPearsonErrors(TGraphErrors* ratioGraph, TF1* ratioFitFunction, TF* denominatorFitFunction)
  {
    for (int pointNo = 0; pointNo < ratioGraph->GetN(); ++pointNo)
    {
      double x, y;
      ratioGraph->GetPoint(pointNo, x, y);
      double errY = ratioGraph->GetErrorY(pointNo);
      double RfuncValue = ratioFitFunction->Eval(x);
      double denomFuncValue = denominatorFitFunction->Eval(x);

      // error formula goes as dR = sqrt( (1-R^2) / (U+V) )
      // R can be replaced by RfuncValue
      // for U+V, need to fit the denominator and then supply that function here

      double newError = sqrt((1 - pow(RfuncValue,2)) / (denomFuncValue));

      cout << "x point: " << x << " original error: " << errY << " new error: " << newError << endl; // print out the before and after errors to make sure they are very close to one another

      // ratioGraph->SetPointError(pointNo, 0, newError); // comment this in once everything is working as it should, and then condense the code and remove extraneous bits
    }
  };
*/

/////////////////////////////////////////////////////////////////////////////////////

   // methods in fit function

  double LMPart(double t, double* par){
    double LM_part = 1;

    if(includeLM)
    {
      int timeBin = lostMuonIntegral->FindBin(t);
      double integralPart = 0;
      integralPart = lostMuonIntegral->GetBinContent(timeBin);
      // for (int bin = fitStartBin; bin <= timeBin; ++bin) integralPart += lostMuonHistogram->GetBinContent(bin) * exp(lostMuonHistogram->GetBinCenter(bin)/defaultLifetime); // this is too slow at the moment, could also switch out defaultLifetime for par[1] potentially, also might need to update with latest bin stuff

      LM_part = 1 - par[parameterIndices[kLoss_string]] * integralPart * 1e-10; // the 1e-10 is just an arbitrary scaling factor
    }

    return LM_part;
  };

  double CBOFreqPart(double t, double* par){
    double frequency = par[parameterIndices[cbo_w_string]];
    if(changingCBOFreq) frequency = par[parameterIndices[cbo_w_string]] * (1 + tM_Aconst*exp(-t/tM_Atau)/(tM_w0*t) + tM_Bconst*exp(-t/tM_Btau)/(tM_w0*t));

    return (frequency * 1e-3); // convert from rad/us to rad/ns
  };

  double VWPart(double t, double* par){
    double VW_part = 1;

    if(includeVW){
      double lifetime = par[parameterIndices[vw_tau_string]] * 1000.;  // convert from us to ns for the fit function
      double frequency;

      if(!changingVWFreq) frequency = par[parameterIndices[vw_w_string]] * 1e-3; // convert to rad/ns for fit function from rad/us, which is the units the parameter is in
      else{
        double fcyc = 2*pi/binWidth; // rad/ns
        double fcbo = (1 + par[parameterIndices[vw_w_string]]) * CBOFreqPart(t, par);
        frequency = fcyc - 2*sqrt(2*fcyc*fcbo - fcbo*fcbo);
      }

       VW_part = (1 + exp(-t/lifetime) * par[parameterIndices[vw_amp_string]] * cos(frequency*t + par[parameterIndices[vw_phi_string]]));
    }

    return VW_part;
  };

  double NCBO2Part(double t, double* par){
    double N2_cbo = 1;

    if(include2CBO){
      double lifetime = par[parameterIndices[cbo_tau_string]] * 1000. * 0.5; // convert from us to ns for the fit function, half the cbo lifetime for the 2 cbo term
      N2_cbo = (1 + exp(-t/lifetime) * par[parameterIndices[cbo_N2_amp_string]] * cos(2*CBOFreqPart(t, par)*t + par[parameterIndices[cbo_N2_phi_string]]));
    }

    return N2_cbo;
  };

  double NCBOPart(double t, double* par){
    double lifetime = par[parameterIndices[cbo_tau_string]] * 1000.; // convert from us to ns for the fit function
    double N_cbo = (1 + exp(-t/lifetime) * par[parameterIndices[cbo_N_amp_string]] * cos(CBOFreqPart(t, par)*t + par[parameterIndices[cbo_N_phi_string]])); // there's no constant out front since it will divide out in the ratio
    // double N_cbo = (1 + (exp(-t/lifetime) * par[parameterIndices[cbo_N_amp_string]] + par[parameterIndices["C"]]) * cos(CBOFreqPart(t, par)*t + par[parameterIndices[cbo_N_phi_string]])); // envelope with constant term + exponential on the amplitude
    // double N_cbo = (1 + (exp(-t/lifetime) * par[parameterIndices[cbo_N_amp_string]] * (1 + par[parameterIndices["C"]] * cos(t/par[parameterIndices["C_{T}"]] - par[parameterIndices["C_{#phi}"]])) ) * cos(CBOFreqPart(t, par)*t + par[parameterIndices[cbo_N_phi_string]])); // envelope with exponential term and extra oscillatory term

    return N_cbo;
  };

  double ACBOPart(double t, double* par){
    double A_cbo;

    if(includeACBO){
      double lifetime = 1.0 * par[parameterIndices[cbo_tau_string]] * 1000.; // convert from us to ns for the fit function
      A_cbo = par[parameterIndices[A_string]] * (1 + exp(-t/lifetime) * par[parameterIndices[cbo_A_amp_string]] * cos(CBOFreqPart(t, par)*t + par[parameterIndices[cbo_A_phi_string]]));
    }
    else A_cbo = par[parameterIndices[A_string]]; // no A cbo

    return A_cbo;
  };

  double PCBOPart(double t, double* par){
    double P_cbo;

    if(includePhiCBO){
      double lifetime = 1.0 * par[parameterIndices[cbo_tau_string]] * 1000.; // convert from us to ns for the fit function
      P_cbo = par[parameterIndices[phi_string]] + exp(-t/lifetime) * par[parameterIndices[cbo_Phi_amp_string]] * cos(CBOFreqPart(t, par)*t + par[parameterIndices[cbo_Phi_phi_string]]);
    }
    else P_cbo = par[parameterIndices[phi_string]]; // no Phi cbo

    return P_cbo;
  };

  double fullRatioFunc(double* x, double* par){
    double time = x[0];
    double freq = blinder->paramToFreq(par[parameterIndices[R_string]]) * 1e-3; // 1e-3 at end to convert to inv ns

    double t0 = 0; // defaultLifetime + startTimeOfFit;

    double f0 = LMPart(time, par) * VWPart(time, par) * NCBO2Part(time, par) * NCBOPart(time, par) * (1 + ACBOPart(time, par) * cos(freq*(time-t0) + PCBOPart(time, par)));
    double fplus  = LMPart(time+halfPeriodGuess, par) * VWPart(time+halfPeriodGuess, par) * NCBO2Part(time+halfPeriodGuess, par) * NCBOPart(time+halfPeriodGuess, par) * (1 + ACBOPart(time+halfPeriodGuess, par) * cos(freq*((time-t0)+halfPeriodGuess) + PCBOPart(time+halfPeriodGuess, par)));
    double fminus = LMPart(time-halfPeriodGuess, par) * VWPart(time-halfPeriodGuess, par) * NCBO2Part(time-halfPeriodGuess, par) * NCBOPart(time-halfPeriodGuess, par) * (1 + ACBOPart(time-halfPeriodGuess, par) * cos(freq*((time-t0)-halfPeriodGuess) + PCBOPart(time-halfPeriodGuess, par)));

    double fitVal = (2*f0 - fplus - fminus) / (2*f0 + fplus + fminus);
    return fitVal;
  };

}; // end FullRatioFit class


FullRatioFit::FullRatioFit(blinding::Blinders* inBlinder){
  blinder = inBlinder;
}

TF1* FullRatioFit::prepareFitFunction(){

  // set up fit function

        auto fullRatioFitFunc = new TF1("fullRatioFitFunc", this, &FullRatioFit::fullRatioFunc, fitStart, fitEnd, int(parameterIndices.size()));
        
        fullRatioFitFunc->SetNpx(10000);
        fullRatioFitFunc->SetLineColor(2);

        for (int parNum = 0; parNum < fullRatioFitFunc->GetNpar(); ++parNum){
          fullRatioFitFunc->SetParName(parNum, parameterNameList[parNum].c_str());
          fullRatioFitFunc->FixParameter(parNum, 0); // fix everything to 0 at the start
        } 
        
        return fullRatioFitFunc;
}

void FullRatioFit::setStartingFitParameters(TF1* fullRatioFitFunc, double amplitudeMultiplier = 1.){ // amplitude multiplier is used to help per calo fits converge which have higher amplitudes than the calorimeter sum fits

        // pull in the fit results from the 3 param ratio fit
        fullRatioFitFunc->FixParameter(parameterIndices[A_string], threeParamRatioFit->GetParameter(parameterIndices[A_string]));
        fullRatioFitFunc->FixParameter(parameterIndices[R_string], threeParamRatioFit->GetParameter(parameterIndices[R_string]));
        fullRatioFitFunc->FixParameter(parameterIndices[phi_string], threeParamRatioFit->GetParameter(parameterIndices[phi_string]));

    if(TmethodFitFunction != 0 && input_FullRatioFitFunction == 0) // set parameters from input T Method function
    {
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_w_string], TmethodFitFunction->GetParameter(TmethodParameterIndices[cbo_w_string]));
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_tau_string], TmethodFitFunction->GetParameter(TmethodParameterIndices[cbo_tau_string]));
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_N_amp_string], amplitudeMultiplier * TmethodFitFunction->GetParameter(TmethodParameterIndices[cbo_N_amp_string]));
          // fullRatioFitFunc->FixParameter(parameterIndices[cbo_N_amp_string], TmethodFitFunction->GetParameter(TmethodParameterIndices[cbo_N_amp_string]));
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_N_phi_string], TmethodFitFunction->GetParameter(TmethodParameterIndices[cbo_N_phi_string]));

        if(include2CBO){
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_N2_amp_string], amplitudeMultiplier * TmethodFitFunction->GetParameter(TmethodParameterIndices[cbo_N2_amp_string]));
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_N2_phi_string], TmethodFitFunction->GetParameter(TmethodParameterIndices[cbo_N2_phi_string]));
        }
        if(includeACBO){
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_A_amp_string], amplitudeMultiplier * TmethodFitFunction->GetParameter(TmethodParameterIndices[cbo_A_amp_string]));
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_A_phi_string], TmethodFitFunction->GetParameter(TmethodParameterIndices[cbo_A_phi_string]));
        }
        if(includePhiCBO){
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_Phi_amp_string], amplitudeMultiplier * TmethodFitFunction->GetParameter(TmethodParameterIndices[cbo_Phi_amp_string]));
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_Phi_phi_string], TmethodFitFunction->GetParameter(TmethodParameterIndices[cbo_Phi_phi_string]));
        }

        if(includeVW){
          fullRatioFitFunc->FixParameter(parameterIndices[vw_w_string],   TmethodFitFunction->GetParameter(TmethodParameterIndices[vw_w_string]));
          fullRatioFitFunc->FixParameter(parameterIndices[vw_tau_string], TmethodFitFunction->GetParameter(TmethodParameterIndices[vw_tau_string]));
          // if(isParamFree(vw_amp_string)) 
          // fullRatioFitFunc->FixParameter(parameterIndices[vw_amp_string], 0.545938*TmethodFitFunction->GetParameter(TmethodParameterIndices[vw_amp_string])); // the factor out front was determined to be the reduction on the VW amplitude in the ratio from the 9d from the average of 50 random seeds
          // else 
          fullRatioFitFunc->FixParameter(parameterIndices[vw_amp_string], TmethodFitFunction->GetParameter(TmethodParameterIndices[vw_amp_string]));
          fullRatioFitFunc->FixParameter(parameterIndices[vw_phi_string], TmethodFitFunction->GetParameter(TmethodParameterIndices[vw_phi_string]));
        }

        if(includeLM){
          fullRatioFitFunc->FixParameter(parameterIndices[kLoss_string], TmethodFitFunction->GetParameter(TmethodParameterIndices[kLoss_string]));
        }
    } 
    else if(input_FullRatioFitFunction != 0) // sets parameters to first ratio fit function - only used when passParamsForward is false and setStartingParamsFromFirstRMethodFit is true
    {
      cout << "setting starting parameters from first ratio method fit" << endl;

          fullRatioFitFunc->FixParameter(parameterIndices[cbo_w_string], input_FullRatioFitFunction->GetParameter(parameterIndices[cbo_w_string]));
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_tau_string], input_FullRatioFitFunction->GetParameter(parameterIndices[cbo_tau_string]));
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_N_amp_string], input_FullRatioFitFunction->GetParameter(parameterIndices[cbo_N_amp_string]));
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_N_phi_string], input_FullRatioFitFunction->GetParameter(parameterIndices[cbo_N_phi_string]));

        if(include2CBO){
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_N2_amp_string], input_FullRatioFitFunction->GetParameter(parameterIndices[cbo_N2_amp_string]));
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_N2_phi_string], input_FullRatioFitFunction->GetParameter(parameterIndices[cbo_N2_phi_string]));
        }
        if(includeACBO){
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_A_amp_string], input_FullRatioFitFunction->GetParameter(parameterIndices[cbo_A_amp_string]));
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_A_phi_string], input_FullRatioFitFunction->GetParameter(parameterIndices[cbo_A_phi_string]));
        }
        if(includePhiCBO){
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_Phi_amp_string], input_FullRatioFitFunction->GetParameter(parameterIndices[cbo_Phi_amp_string]));
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_Phi_phi_string], input_FullRatioFitFunction->GetParameter(parameterIndices[cbo_Phi_phi_string]));
        }

        if(includeVW){
          fullRatioFitFunc->FixParameter(parameterIndices[vw_w_string],   input_FullRatioFitFunction->GetParameter(parameterIndices[vw_w_string]));
          fullRatioFitFunc->FixParameter(parameterIndices[vw_tau_string], input_FullRatioFitFunction->GetParameter(parameterIndices[vw_tau_string]));
          fullRatioFitFunc->FixParameter(parameterIndices[vw_amp_string], input_FullRatioFitFunction->GetParameter(parameterIndices[vw_amp_string]));
          fullRatioFitFunc->FixParameter(parameterIndices[vw_phi_string], input_FullRatioFitFunction->GetParameter(parameterIndices[vw_phi_string]));
        }

        if(includeLM){
          // fullRatioFitFunc->FixParameter(parameterIndices[kLoss_string], input_FullRatioFitFunction->GetParameter(parameterIndices[kLoss_string]));
          fullRatioFitFunc->FixParameter(parameterIndices[kLoss_string], TmethodFitFunction->GetParameter(TmethodParameterIndices[kLoss_string])); // hack-y thing so that I can fix some parameters to the first ratio fit function but not k_loss
        }
    } 
    else // set parameters from defaults in ratioAnalysisConfig.hh
    {
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_w_string], startingCBOFreq);
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_tau_string], startingCBOTau);
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_N_amp_string], amplitudeMultiplier * startingCBONAmp);
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_N_phi_string], startingCBONPhase);

        if(include2CBO){
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_N2_amp_string], amplitudeMultiplier * startingN2CBO_amp);
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_N2_phi_string], startingN2CBO_phi);
        }
        if(includeACBO){
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_A_amp_string], amplitudeMultiplier * startingACBO_amp);
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_A_phi_string], startingACBO_phi);
        }
        if(includePhiCBO){
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_Phi_amp_string], amplitudeMultiplier * startingPhiCBO_amp);
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_Phi_phi_string], startingPhiCBO_phi);
        }

        if(includeVW){
          fullRatioFitFunc->FixParameter(parameterIndices[vw_w_string],   startingVWFreq);
          fullRatioFitFunc->FixParameter(parameterIndices[vw_tau_string], startingVWTau);
          fullRatioFitFunc->FixParameter(parameterIndices[vw_amp_string], startingVWAmp);
          fullRatioFitFunc->FixParameter(parameterIndices[vw_phi_string], startingVWPhase);
        }

        if(includeLM){
          fullRatioFitFunc->FixParameter(parameterIndices[kLoss_string], startingLMAmp);
        }
    } 
// 
    // fullRatioFitFunc->FixParameter(parameterIndices["C"], 0);
    // fullRatioFitFunc->FixParameter(parameterIndices["C"], TmethodFitFunction->GetParameter(TmethodParameterIndices["C"]));

        // fix parameter value for that which is being scanned over
        // if(!isParamFree(parameterToScan)) fullRatioFitFunc->FixParameter(parameterIndices[parameterToScan],  fullRatioFitFunc->GetParameter(parameterIndices[parameterToScan]) + scannedParameterValueChange);
        if(!isParamFree(parameterToScan)) fullRatioFitFunc->FixParameter(parameterIndices[parameterToScan],  parameterScanStart + scannedParameterValueChange);

}

double* FullRatioFit::doFullRatioFit(TF1* fullRatioFitFunc, TGraphErrors* ratioGraph){

        // fit cbo

          fixAllParameters(fullRatioFitFunc);

          if(isParamFree(cbo_N_phi_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_N_phi_string], -4*pi, 4*pi);
          setDiffParLimits(fullRatioFitFunc, {cbo_N_phi_string});

        ratioGraph->Fit(fullRatioFitFunc, "QR"); // fit to cbo phase first - helps with calo fits

          if(isParamFree(cbo_w_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_w_string], .9*fullRatioFitFunc->GetParameter(parameterIndices[cbo_w_string]), 1.1*fullRatioFitFunc->GetParameter(parameterIndices[cbo_w_string]));
          if(isParamFree(cbo_N_phi_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_N_phi_string], -4*pi, 4*pi);
          setDiffParLimits(fullRatioFitFunc, {cbo_w_string, cbo_N_phi_string});

        ratioGraph->Fit(fullRatioFitFunc, "QR"); // fit to cbo phase and frequency

          if(isParamFree(cbo_tau_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_tau_string], 10, 600); // typically fix for calo fits or fit start scans
          if(isParamFree(cbo_N_amp_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_N_amp_string], 0, .05);
          // if(isParamFree(cbo_tau_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_tau_string], 10, 100); // typically fix for calo fits or fit start scans
          // if(isParamFree(cbo_N_amp_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_N_amp_string], 0, .003);
          setDiffParLimits(fullRatioFitFunc, {cbo_w_string, cbo_N_phi_string, cbo_tau_string, cbo_N_amp_string});

// extra cbo envelope parameters
//hasn't been fully updated with new set starting params method
          // fullRatioFitFunc->SetParLimits(parameterIndices["C"], -0.1, 0.002);

          // fullRatioFitFunc->FixParameter(parameterIndices["C_{T}"], 114500);
          // fullRatioFitFunc->SetParLimits(parameterIndices["C_{T}"], 10000, 300000);

          // fullRatioFitFunc->SetParameter(parameterIndices["C_{#phi}"], 0.482);
          // fullRatioFitFunc->SetParLimits(parameterIndices["C_{#phi}"], -2*pi, 2*pi);

        ratioGraph->Fit(fullRatioFitFunc, "QR"); // fit to cbo parameters

/////////////////////////////////////////////////////////////////////////////////////

        // fit 2 CBO N term
        if(include2CBO){
          fixAllParameters(fullRatioFitFunc);
          if(isParamFree(cbo_N2_amp_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_N2_amp_string], 0, .01);
          if(isParamFree(cbo_N2_phi_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_N2_phi_string], -4*pi, 4*pi);
          setDiffParLimits(fullRatioFitFunc, {cbo_N2_amp_string, cbo_N2_phi_string});
          ratioGraph->Fit(fullRatioFitFunc, "QR");
        }

        // fit A CBO term
        if(includeACBO){
          fixAllParameters(fullRatioFitFunc);
          if(isParamFree(cbo_A_amp_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_A_amp_string], 0, .02);
          if(isParamFree(cbo_A_phi_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_A_phi_string], -3*pi, 3*pi);
          setDiffParLimits(fullRatioFitFunc, {cbo_A_amp_string, cbo_A_phi_string});
          ratioGraph->Fit(fullRatioFitFunc, "QR");
        }

        // fit Phi CBO term
        if(includePhiCBO){
          fixAllParameters(fullRatioFitFunc);
          if(isParamFree(cbo_Phi_amp_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_Phi_amp_string], 0, .02);
          if(isParamFree(cbo_Phi_phi_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_Phi_phi_string], -4*pi, 4*pi);
          setDiffParLimits(fullRatioFitFunc, {cbo_Phi_amp_string, cbo_Phi_phi_string});
          ratioGraph->Fit(fullRatioFitFunc, "QR");
        }

/////////////////////////////////////////////////////////////////////////////////////

        // fit VW
        if(includeVW){
          fixAllParameters(fullRatioFitFunc);

          if(includeChangingVWFreq){ if(isParamFree(vw_w_string)) fullRatioFitFunc->SetParLimits(parameterIndices[vw_w_string], 0, 0.03); }
          else{ if(isParamFree(vw_w_string)) fullRatioFitFunc->SetParLimits(parameterIndices[vw_w_string], .8*fullRatioFitFunc->GetParameter(parameterIndices[vw_w_string]), 1.2*fullRatioFitFunc->GetParameter(parameterIndices[vw_w_string])); }

          if(isParamFree(vw_phi_string)) fullRatioFitFunc->SetParLimits(parameterIndices[vw_phi_string], -2*pi, 2*pi);
          setDiffParLimits(fullRatioFitFunc, {vw_w_string, vw_phi_string});

          ratioGraph->Fit(fullRatioFitFunc, "QR"); // fit to freq and phase

          fixAllParameters(fullRatioFitFunc);

          if(isParamFree(vw_tau_string)) fullRatioFitFunc->SetParLimits(parameterIndices[vw_tau_string], 1, 75); // 50
          if(isParamFree(vw_amp_string)) fullRatioFitFunc->SetParLimits(parameterIndices[vw_amp_string], 0, 1); // 0.5
          setDiffParLimits(fullRatioFitFunc, {vw_tau_string, vw_amp_string});

          ratioGraph->Fit(fullRatioFitFunc, "QR"); // fit to VW amp and lifetime
          
          if(includeChangingVWFreq){ if(isParamFree(vw_w_string)) fullRatioFitFunc->SetParLimits(parameterIndices[vw_w_string], 0, 0.03); }
          else{ if(isParamFree(vw_w_string)) fullRatioFitFunc->SetParLimits(parameterIndices[vw_w_string], .8*fullRatioFitFunc->GetParameter(parameterIndices[vw_w_string]), 1.2*fullRatioFitFunc->GetParameter(parameterIndices[vw_w_string])); }

          if(isParamFree(vw_phi_string)) fullRatioFitFunc->SetParLimits(parameterIndices[vw_phi_string], -2*pi, 2*pi);
          setDiffParLimits(fullRatioFitFunc, {vw_w_string, vw_phi_string, vw_tau_string, vw_amp_string});

          ratioGraph->Fit(fullRatioFitFunc, "QR"); // fit to VW parameters
        }

/////////////////////////////////////////////////////////////////////////////////////

  // fit all parameters at once

  fitAllParameters(fullRatioFitFunc, ratioGraph, fullRatioFitFunc->GetParameters());
  return fullRatioFitFunc->GetParameters();

}

double* FullRatioFit::fitAllParameters(TF1* fullRatioFitFunc, TGraphErrors* ratioGraph, double* pars){

        TFitResultPtr fitResult;
        int fitStatus = -1;

        fullRatioFitFunc->SetParameters(pars); // any parameters that were fixed before will still be fixed

        // explicitly list out all parameter limits even though some of the lines of code are redundant
          
          // if(isParamFree(A_string)) fullRatioFitFunc->SetParLimits(parameterIndices[A_string], 0.3, 0.5);
          // if(isParamFree(R_string)) fullRatioFitFunc->SetParLimits(parameterIndices[R_string], -100, 100);
          // if(isParamFree(phi_string)) fullRatioFitFunc->SetParLimits(parameterIndices[phi_string], 2, 2.5);
          if(isParamFree(A_string)) fullRatioFitFunc->ReleaseParameter(parameterIndices[A_string]);
          if(isParamFree(R_string)) fullRatioFitFunc->ReleaseParameter(parameterIndices[R_string]);
          if(isParamFree(phi_string)) fullRatioFitFunc->ReleaseParameter(parameterIndices[phi_string]);

          // main cbo parameters
          if(isParamFree(cbo_w_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_w_string], .9*fullRatioFitFunc->GetParameter(parameterIndices[cbo_w_string]), 1.1*fullRatioFitFunc->GetParameter(parameterIndices[cbo_w_string]));
          if(isParamFree(cbo_tau_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_tau_string], 10, 600); // typically fix for calo fits or fit start scans
          if(isParamFree(cbo_N_amp_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_N_amp_string], 0, .05);
          // if(isParamFree(cbo_tau_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_tau_string], 10, 100); // typically fix for calo fits or fit start scans
          // if(isParamFree(cbo_N_amp_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_N_amp_string], 0, .003);
          fullRatioFitFunc->FixParameter(parameterIndices[cbo_N_phi_string], fmod(fullRatioFitFunc->GetParameter(parameterIndices[cbo_N_phi_string]), 2*pi));
          if(isParamFree(cbo_N_phi_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_N_phi_string], -3*pi, 3*pi);

          // 2 cbo parameters
          if(include2CBO){
            if(isParamFree(cbo_N2_amp_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_N2_amp_string], 0, .01);
            fullRatioFitFunc->SetParameter(parameterIndices[cbo_N2_phi_string], fmod(fullRatioFitFunc->GetParameter(parameterIndices[cbo_N2_phi_string]), 2*pi));
            if(isParamFree(cbo_N2_phi_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_N2_phi_string], -3*pi, 3*pi);
          }

          // A cbo parameters
          if(includeACBO){
            if(isParamFree(cbo_A_amp_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_A_amp_string], 0, .02);
            fullRatioFitFunc->SetParameter(parameterIndices[cbo_A_phi_string], fmod(fullRatioFitFunc->GetParameter(parameterIndices[cbo_A_phi_string]), 2*pi));
            if(isParamFree(cbo_A_phi_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_A_phi_string], -3*pi, 3*pi);
          }

          // phi cbo parameters
          if(includePhiCBO){
            if(isParamFree(cbo_Phi_amp_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_Phi_amp_string], 0, .02);
            fullRatioFitFunc->SetParameter(parameterIndices[cbo_Phi_phi_string], fmod(fullRatioFitFunc->GetParameter(parameterIndices[cbo_Phi_phi_string]), 2*pi));
            if(isParamFree(cbo_Phi_phi_string)) fullRatioFitFunc->SetParLimits(parameterIndices[cbo_Phi_phi_string], -3*pi, 3*pi);
          }

          // vw parameters
          if(includeVW){
            if(includeChangingVWFreq){ if(isParamFree(vw_w_string)) fullRatioFitFunc->SetParLimits(parameterIndices[vw_w_string], 0, 0.03); }
            else{ if(isParamFree(vw_w_string)) fullRatioFitFunc->SetParLimits(parameterIndices[vw_w_string], .8*fullRatioFitFunc->GetParameter(parameterIndices[vw_w_string]), 1.2*fullRatioFitFunc->GetParameter(parameterIndices[vw_w_string])); }
            if(isParamFree(vw_tau_string)) fullRatioFitFunc->SetParLimits(parameterIndices[vw_tau_string], 1, 75); // 50
            if(isParamFree(vw_amp_string)) fullRatioFitFunc->SetParLimits(parameterIndices[vw_amp_string], 0, 1); // 0.5
            fullRatioFitFunc->FixParameter(parameterIndices[vw_phi_string], fmod(fullRatioFitFunc->GetParameter(parameterIndices[vw_phi_string]), 2*pi));
            if(isParamFree(vw_phi_string)) fullRatioFitFunc->SetParLimits(parameterIndices[vw_phi_string], -2*pi, 2*pi);
          }

          // lost muons
          if(includeLM && isParamFree(kLoss_string)){ fullRatioFitFunc->SetParLimits(parameterIndices[kLoss_string], -5, 15); }


          // other parameters - not setup with newest set stating params method

          // fullRatioFitFunc->SetParLimits(parameterIndices["C"], -0.1, 0.002);
          // // fullRatioFitFunc->SetParLimits(parameterIndices["C_{T}"], 10000, 300000);
          // fullRatioFitFunc->SetParLimits(parameterIndices["C_{#phi}"], -2*pi, 2*pi);

/////////////////////////////////////////////////////////////////////////////////////

          // set updated parameter limits for all parameters
          for (uint i = 0; i < updatedParLimits.size(); ++i) if(isParamFree(updatedParLimits.at(i).first)) fullRatioFitFunc->SetParLimits(parameterIndices[updatedParLimits.at(i).first], updatedParLimits.at(i).second.first, updatedParLimits.at(i).second.second);

/////////////////////////////////////////////////////////////////////////////////////

        // printParamsAndLimits(fullRatioFitFunc);

        fitResult = ratioGraph->Fit(fullRatioFitFunc, RMethod_lastFitOptions.c_str());
        // setPearsonErrors(ratioGraph, fullRatioFitFunc); // for doing Pearson chi2 fit, no P option for T Graphs - could do things manually though - doesn't currently work
        // fitResult = ratioGraph->Fit(fullRatioFitFunc, RMethod_lastFitOptions.c_str());


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

        fitStatus = fitResult;

        // if(fitStatus != 0)
        // {
        //   cout << "Ratio CBO fit returned non-zero fit status. Exiting with fit status: " << fitStatus << endl;
        //   exit(1);
        // }

        fitResult->GetCorrelationMatrix().Write("FullRatioFitCorrMatrix"); // write the correlation matrix to the file

    return fullRatioFitFunc->GetParameters();

}

#endif //! FULLRATIOFIT_HH
