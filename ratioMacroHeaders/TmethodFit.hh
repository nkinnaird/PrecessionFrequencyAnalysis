// 3-31-20: T Method fit header file. Do fits using fiveParamFit.hh as input to this.

#ifndef TMETHODFIT_HH
#define TMETHODFIT_HH

class TmethodFit{

public:

  TmethodFit(blinding::Blinders* inBlinder);

  TF1* prepareFitFunction();
  void setStartingFitParameters(TF1* TmethodFitFunc, double amplitudeMultiplier);

  void clearUpdatedParLimits(){
    updatedParLimits.clear();  
  };

  double* doTMethodFit(TF1* TmethodFitFunc, TH1F* histToFit);
  double* fitAllTMethodParameters(TF1* TmethodFitFunc, TH1F* histToFit, double* pars);

  void setBinWidth(double inputBinWidth){
    binWidth = inputBinWidth;
  };

  void setFitRange(double fitStartTime, double fitEndTime){
    fitStart = fitStartTime;
    fitEnd = fitEndTime;
  };

  void setFiveParamFunction(TF1* input_fiveFit){
    inputFiveParamFit = input_fiveFit;
  };

  void setTMethodFunction(TF1* input_TmethodFit){
    input_TmethodFitFunction = input_TmethodFit;
  };

  void setupFitParams(bool changingCBOFrequency, bool changingVWFrequency, bool N2_cbo, bool A_cbo, bool Phi_cbo, bool incVW, bool incLM, bool incBR){
    changingCBOFreq = changingCBOFrequency;
    changingVWFreq = changingVWFrequency;
    include2CBO = N2_cbo;
    includeACBO = A_cbo;
    includePhiCBO = Phi_cbo;
    includeVW = incVW;
    includeLM = incLM;
    includeBR = incBR;

    parameterNameList.insert(parameterNameList.begin(), {N_string, tau_string, A_string, R_string, phi_string, cbo_w_string, cbo_tau_string, cbo_N_amp_string, cbo_N_phi_string});
    if(include2CBO) parameterNameList.insert(parameterNameList.end(), {cbo_N2_amp_string, cbo_N2_phi_string});
    if(includeACBO) parameterNameList.insert(parameterNameList.end(), {cbo_A_amp_string, cbo_A_phi_string});
    if(includePhiCBO) parameterNameList.insert(parameterNameList.end(), {cbo_Phi_amp_string, cbo_Phi_phi_string});
    if(includeVW) parameterNameList.insert(parameterNameList.end(), {vw_w_string, vw_tau_string, vw_amp_string, vw_phi_string});
    if(includeLM) parameterNameList.push_back(kLoss_string);
    if(includeBR) parameterNameList.insert(parameterNameList.end(), {br_tau_string, br_amp_string});

    // parameterNameList.insert(parameterNameList.end(), "C");

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

  void setLostMuonHists(TH1F* lostMuonHist, TH1F* integralHist){
    lostMuonHistogram = lostMuonHist;
    lostMuonIntegral = integralHist;
  };

  int getNumFitParams(){
    return numFitParams;
  };

  map<string, int> getParamIndicesMap(){
    return parameterIndices;
  };

  void printParamsAndLimits(TF1* TmethodFitFunc){
      for (int parNum = 0; parNum < TmethodFitFunc->GetNpar(); ++parNum){
        double low, high = 0;
        double value = TmethodFitFunc->GetParameter(parNum);
        TmethodFitFunc->GetParLimits(parNum, low, high);
        cout << "Parameter: " << parNum << " name: " << TmethodFitFunc->GetParName(parNum) << " value: " << value << " limits: " << low << " - " << high;
          if((value > high || value < low) && !(low == 0 && high == 0)) cout << " OUTSIDE LIMIT" << endl;
          else cout << endl;
      }
  };

  void fixAllParameters(TF1* TmethodFitFunc){
    for (int parNum = 0; parNum < TmethodFitFunc->GetNpar(); ++parNum) TmethodFitFunc->FixParameter(parNum, TmethodFitFunc->GetParameter(parNum));
  };

  void fillDiffParLimits(string parName, double low, double high){
    updatedParLimits.emplace_back(parName, make_pair(low, high));
  };
  
  // update limits dynamically for individual fits using passed in vector which checks if those parameters are included in the updatedParLimits vector
  void setDiffParLimits(TF1* TmethodFitFunc, vector<string> parsToUpdate){
    for (uint index = 0; index < parsToUpdate.size(); ++index)
    {
      auto it = std::find_if( updatedParLimits.begin(), updatedParLimits.end(),
      [&parsToUpdate, index](const pair<string, pair<double, double> >& element){ return parsToUpdate.at(index) == element.first; } );

      if(it != updatedParLimits.end() && isParamFree(parsToUpdate.at(index))) TmethodFitFunc->SetParLimits(parameterIndices[parsToUpdate.at(index)], (*it).second.first, (*it).second.second);
    }
  };

private:

  blinding::Blinders* blinder;
  double binWidth = 0;
  double fitStart, fitEnd;

  vector<string> parameterNameList; // initialized in setupFitParams
  map<string, int> parameterIndices; // map of strings to parameter indices for more flexible fitting with respect to additional parameters
  int numFitParams = 0;

  TF1* inputFiveParamFit = 0;
  TF1* input_TmethodFitFunction = 0;

  vector<string> fitFixedParameters = TmethodFit_added_fixedParameters; // default initial fixed parameters to T method added fixed parameters

  bool changingCBOFreq = false;
  bool changingVWFreq = false;
  bool include2CBO = false;
  bool includeACBO = false;
  bool includePhiCBO = false;
  bool includeVW = false;
  bool includeBR = false;
  
  bool includeLM = false;
  TH1F* lostMuonHistogram;
  TH1F* lostMuonIntegral;
  bool useDefaultIntegral = true;

  vector< pair<string, pair<double, double> > > updatedParLimits;

/////////////////////////////////////////////////////////////////////////////////////

  // methods in fit function

  double BRPart(double t, double* par){
    double BR_part = 1;
    
    if(includeBR){
      double lifetime = par[parameterIndices[br_tau_string]] * 1000.; // convert from us to ns for the fit function, do this in if statement so that parameterIndices map doesn't add br_tau_string
      BR_part = (1 + exp(-t/lifetime) * par[parameterIndices[br_amp_string]]);
    } 

    return BR_part;
  };

  double LMPart(double t, double* par){
    double muonLossPart = 1;

    if(includeLM)
    {
      int timeBin = lostMuonIntegral->FindBin(t);
      double integralPart = 0;

      if(useDefaultIntegral) integralPart = lostMuonIntegral->GetBinContent(timeBin);
      else{ // this is extremely slow and doesn't appear to change the final fitted parameters at all - tightening up parameter limits may speed this up considerably
        int fitStartBin = int(startTimeOfFit/binWidth) + 1; // this bin should be the same for each iteration, that being the start bin of the first fit
        double lifetime = par[parameterIndices[tau_string]] * 1000.;
        for (int bin = fitStartBin; bin <= timeBin; ++bin) integralPart += lostMuonHistogram->GetBinContent(bin) * exp(lostMuonHistogram->GetBinCenter(bin)/lifetime); // calculate the lost muon integral directly with the variable muon lifetime parameter
      }

      muonLossPart = 1 - par[parameterIndices[kLoss_string]] * integralPart * 1e-10; // the 1e-10 is just an arbitrary scaling factor
    }

    return muonLossPart;
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
    double N_cbo = (1 + exp(-t/lifetime) * par[parameterIndices[cbo_N_amp_string]] * cos(CBOFreqPart(t, par)*t + par[parameterIndices[cbo_N_phi_string]]));
    // double N_cbo = (1 + (exp(-t/lifetime) * par[parameterIndices[cbo_N_amp_string]] + par[parameterIndices["C"]]) * cos(CBOFreqPart(t, par)*t + par[parameterIndices[cbo_N_phi_string]])); // envelope with constant term + exponential on the amplitude

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

  double totalFunction(double* x, double* par){
    
    double time = x[0];
    double freq = blinder->paramToFreq(par[parameterIndices[R_string]]) * 1e-3; // 1e-3 at end to convert to inv ns
    double lifetime = par[parameterIndices[tau_string]] * 1000.; // convert from us to ns for fit function

    double t0 = 0; // defaultLifetime + startTimeOfFit;

    double fitVal = BRPart(time, par) * LMPart(time, par) * VWPart(time, par) * NCBO2Part(time, par) * NCBOPart(time, par) * par[parameterIndices[N_string]] * exp(-time/lifetime) * (1 + ACBOPart(time, par) * cos(freq*(time-t0) + PCBOPart(time, par)));
    return fitVal;
  };

/////////////////////////////////////////////////////////////////////////////////////

}; // end TmethodFit class


TmethodFit::TmethodFit(blinding::Blinders* inBlinder){
  blinder = inBlinder;
}

TF1* TmethodFit::prepareFitFunction(){

        auto TmethodFitFunc = new TF1("TmethodFitFunc", this, &TmethodFit::totalFunction, fitStart, fitEnd, int(parameterIndices.size()));
        TmethodFitFunc->SetNpx(10000);
        TmethodFitFunc->SetLineColor(2);

        for (int parNum = 0; parNum < TmethodFitFunc->GetNpar(); ++parNum){
          TmethodFitFunc->SetParName(parNum, parameterNameList[parNum].c_str());
          TmethodFitFunc->FixParameter(parNum, 0); // fix everything to 0 at the start
        }

        return TmethodFitFunc;
}

void TmethodFit::setStartingFitParameters(TF1* TmethodFitFunc, double amplitudeMultiplier = 1.){ // amplitude multiplier is used to help per calo fits converge which have higher amplitudes than the calorimeter sum fits

      useDefaultIntegral = true;

        // pull in the fit results from the 5 param fit
        TmethodFitFunc->FixParameter(parameterIndices[N_string], inputFiveParamFit->GetParameter(parameterIndices[N_string]));
        TmethodFitFunc->FixParameter(parameterIndices[tau_string], inputFiveParamFit->GetParameter(parameterIndices[tau_string]));
        TmethodFitFunc->FixParameter(parameterIndices[A_string], inputFiveParamFit->GetParameter(parameterIndices[A_string]));
        TmethodFitFunc->FixParameter(parameterIndices[R_string], inputFiveParamFit->GetParameter(parameterIndices[R_string]));
        TmethodFitFunc->FixParameter(parameterIndices[phi_string], inputFiveParamFit->GetParameter(parameterIndices[phi_string]));

    if(input_TmethodFitFunction != 0)
    {
          TmethodFitFunc->FixParameter(parameterIndices[cbo_w_string], input_TmethodFitFunction->GetParameter(parameterIndices[cbo_w_string]));
          TmethodFitFunc->FixParameter(parameterIndices[cbo_tau_string], input_TmethodFitFunction->GetParameter(parameterIndices[cbo_tau_string]));
          // TmethodFitFunc->FixParameter(parameterIndices[cbo_N_amp_string], amplitudeMultiplier * input_TmethodFitFunction->GetParameter(parameterIndices[cbo_N_amp_string]));
          TmethodFitFunc->FixParameter(parameterIndices[cbo_N_amp_string], input_TmethodFitFunction->GetParameter(parameterIndices[cbo_N_amp_string]));
          TmethodFitFunc->FixParameter(parameterIndices[cbo_N_phi_string], input_TmethodFitFunction->GetParameter(parameterIndices[cbo_N_phi_string]));

        if(include2CBO){
          TmethodFitFunc->FixParameter(parameterIndices[cbo_N2_amp_string], amplitudeMultiplier * input_TmethodFitFunction->GetParameter(parameterIndices[cbo_N2_amp_string]));
          TmethodFitFunc->FixParameter(parameterIndices[cbo_N2_phi_string], input_TmethodFitFunction->GetParameter(parameterIndices[cbo_N2_phi_string]));
        }
        if(includeACBO){
          TmethodFitFunc->FixParameter(parameterIndices[cbo_A_amp_string], amplitudeMultiplier * input_TmethodFitFunction->GetParameter(parameterIndices[cbo_A_amp_string]));
          TmethodFitFunc->FixParameter(parameterIndices[cbo_A_phi_string], input_TmethodFitFunction->GetParameter(parameterIndices[cbo_A_phi_string]));
        }
        if(includePhiCBO){
          TmethodFitFunc->FixParameter(parameterIndices[cbo_Phi_amp_string], amplitudeMultiplier * input_TmethodFitFunction->GetParameter(parameterIndices[cbo_Phi_amp_string]));
          TmethodFitFunc->FixParameter(parameterIndices[cbo_Phi_phi_string], input_TmethodFitFunction->GetParameter(parameterIndices[cbo_Phi_phi_string]));
        }

        if(includeVW){
          TmethodFitFunc->FixParameter(parameterIndices[vw_w_string],   input_TmethodFitFunction->GetParameter(parameterIndices[vw_w_string]));
          TmethodFitFunc->FixParameter(parameterIndices[vw_tau_string], input_TmethodFitFunction->GetParameter(parameterIndices[vw_tau_string]));
          TmethodFitFunc->FixParameter(parameterIndices[vw_amp_string], input_TmethodFitFunction->GetParameter(parameterIndices[vw_amp_string]));
          TmethodFitFunc->FixParameter(parameterIndices[vw_phi_string], input_TmethodFitFunction->GetParameter(parameterIndices[vw_phi_string]));
        }

        if(includeLM){
          TmethodFitFunc->FixParameter(parameterIndices[kLoss_string], input_TmethodFitFunction->GetParameter(parameterIndices[kLoss_string]));
        }

        if(includeBR){
          TmethodFitFunc->FixParameter(parameterIndices[br_tau_string], input_TmethodFitFunction->GetParameter(parameterIndices[br_amp_string]));
          TmethodFitFunc->FixParameter(parameterIndices[br_amp_string], input_TmethodFitFunction->GetParameter(parameterIndices[br_amp_string]));        
        }
    } 
    else
    {
          TmethodFitFunc->FixParameter(parameterIndices[cbo_w_string], startingCBOFreq);
          TmethodFitFunc->FixParameter(parameterIndices[cbo_tau_string], startingCBOTau);
          TmethodFitFunc->FixParameter(parameterIndices[cbo_N_amp_string], amplitudeMultiplier * startingCBONAmp);
          TmethodFitFunc->FixParameter(parameterIndices[cbo_N_phi_string], startingCBONPhase);

        if(include2CBO){
          TmethodFitFunc->FixParameter(parameterIndices[cbo_N2_amp_string], amplitudeMultiplier * startingN2CBO_amp);
          TmethodFitFunc->FixParameter(parameterIndices[cbo_N2_phi_string], startingN2CBO_phi);
        }
        if(includeACBO){
          TmethodFitFunc->FixParameter(parameterIndices[cbo_A_amp_string], amplitudeMultiplier * startingACBO_amp);
          TmethodFitFunc->FixParameter(parameterIndices[cbo_A_phi_string], startingACBO_phi);
        }
        if(includePhiCBO){
          TmethodFitFunc->FixParameter(parameterIndices[cbo_Phi_amp_string], amplitudeMultiplier * startingPhiCBO_amp);
          TmethodFitFunc->FixParameter(parameterIndices[cbo_Phi_phi_string], startingPhiCBO_phi);
        }

        if(includeVW){
          TmethodFitFunc->FixParameter(parameterIndices[vw_w_string],   startingVWFreq);
          TmethodFitFunc->FixParameter(parameterIndices[vw_tau_string], startingVWTau);
          TmethodFitFunc->FixParameter(parameterIndices[vw_amp_string], startingVWAmp);
          TmethodFitFunc->FixParameter(parameterIndices[vw_phi_string], startingVWPhase);
        }

        if(includeLM){
          TmethodFitFunc->FixParameter(parameterIndices[kLoss_string], startingLMAmp);
        }

        if(includeBR){
          TmethodFitFunc->FixParameter(parameterIndices[br_tau_string], startingBRTau);
          // TmethodFitFunc->FixParameter(parameterIndices[br_amp_string], 0);//startingBRAmp);
          TmethodFitFunc->FixParameter(parameterIndices[br_amp_string], startingBRAmp);
        }
    } 

    // TmethodFitFunc->FixParameter(parameterIndices["C"], 0);

        // fix parameter value for that which is being scanned over - not fully implemented for T method fit
        // if(!isParamFree(parameterToScan)) TmethodFitFunc->FixParameter(parameterIndices[parameterToScan],  TmethodFitFunc->GetParameter(parameterIndices[parameterToScan]) + scannedParameterValueChange);

}

double* TmethodFit::doTMethodFit(TF1* TmethodFitFunc, TH1F* histToFit){

        // fit lost muons, N, and tau

        fixAllParameters(TmethodFitFunc);
          
          if(isParamFree(N_string)) TmethodFitFunc->ReleaseParameter(parameterIndices[N_string]);
          // if(isParamFree(tau_string)) TmethodFitFunc->SetParLimits(parameterIndices[tau_string], 64.2, 64.6);
          if(isParamFree(tau_string)) TmethodFitFunc->ReleaseParameter(parameterIndices[tau_string]);
          if(includeLM && isParamFree(kLoss_string)){ TmethodFitFunc->SetParLimits(parameterIndices[kLoss_string], -35, 35); }
          setDiffParLimits(TmethodFitFunc, {N_string, tau_string, kLoss_string});

        histToFit->Fit(TmethodFitFunc, "QR"); // fit to muon loss term, lifetime, and N

/////////////////////////////////////////////////////////////////////////////////////

        // fit cbo

        fixAllParameters(TmethodFitFunc);

          if(isParamFree(cbo_N_phi_string)) TmethodFitFunc->SetParLimits(parameterIndices[cbo_N_phi_string], -4*pi, 4*pi);
          setDiffParLimits(TmethodFitFunc, {cbo_N_phi_string});

        histToFit->Fit(TmethodFitFunc, "QR"); // fit to cbo phase first only - helps with calo fits

          if(isParamFree(cbo_w_string)) TmethodFitFunc->SetParLimits(parameterIndices[cbo_w_string], .9*TmethodFitFunc->GetParameter(parameterIndices[cbo_w_string]), 1.1*TmethodFitFunc->GetParameter(parameterIndices[cbo_w_string]));
          if(isParamFree(cbo_N_phi_string)) TmethodFitFunc->SetParLimits(parameterIndices[cbo_N_phi_string], -4*pi, 4*pi);
          setDiffParLimits(TmethodFitFunc, {cbo_w_string, cbo_N_phi_string});

        histToFit->Fit(TmethodFitFunc, "QR"); // fit to cbo phase and frequency

          if(isParamFree(cbo_tau_string)) TmethodFitFunc->SetParLimits(parameterIndices[cbo_tau_string], 10, 500);
          if(isParamFree(cbo_N_amp_string)) TmethodFitFunc->SetParLimits(parameterIndices[cbo_N_amp_string], 0, .05);
          setDiffParLimits(TmethodFitFunc, {cbo_w_string, cbo_N_phi_string, cbo_tau_string, cbo_N_amp_string});

          // extra cbo envelope parameters
          // TmethodFitFunc->SetParLimits(parameterIndices["C"], -0.4, 0.4);

        histToFit->Fit(TmethodFitFunc, "QR"); // fit to all cbo parameters

/////////////////////////////////////////////////////////////////////////////////////

        // for beam relaxation term

        if(includeBR){
          fixAllParameters(TmethodFitFunc);
          if(isParamFree(br_tau_string)) TmethodFitFunc->SetParLimits(parameterIndices[br_tau_string], 0, 100);
          if(isParamFree(br_amp_string)) TmethodFitFunc->SetParLimits(parameterIndices[br_amp_string], -.05, .05);
          setDiffParLimits(TmethodFitFunc, {br_tau_string, br_amp_string});
          histToFit->Fit(TmethodFitFunc, "QR");
        }

/////////////////////////////////////////////////////////////////////////////////////

        // fit 2 CBO term
        if(include2CBO){
          fixAllParameters(TmethodFitFunc);
          if(isParamFree(cbo_N2_amp_string)) TmethodFitFunc->SetParLimits(parameterIndices[cbo_N2_amp_string], 0, .01);
          if(isParamFree(cbo_N2_phi_string)) TmethodFitFunc->SetParLimits(parameterIndices[cbo_N2_phi_string], -3*pi, 3*pi);
          setDiffParLimits(TmethodFitFunc, {cbo_N2_amp_string, cbo_N2_phi_string});
          histToFit->Fit(TmethodFitFunc, "QR");
        }

        // fit A CBO term
        if(includeACBO){
          fixAllParameters(TmethodFitFunc);
          if(isParamFree(cbo_A_amp_string)) TmethodFitFunc->SetParLimits(parameterIndices[cbo_A_amp_string], 0, .02);
          if(isParamFree(cbo_A_phi_string)) TmethodFitFunc->SetParLimits(parameterIndices[cbo_A_phi_string], -3*pi, 3*pi);
          setDiffParLimits(TmethodFitFunc, {cbo_A_amp_string, cbo_A_phi_string});
          histToFit->Fit(TmethodFitFunc, "QR");
        }

        // fit Phi CBO term
        if(includePhiCBO){
          fixAllParameters(TmethodFitFunc);
          if(isParamFree(cbo_Phi_amp_string)) TmethodFitFunc->SetParLimits(parameterIndices[cbo_Phi_amp_string], 0, .01);
          if(isParamFree(cbo_Phi_phi_string)) TmethodFitFunc->SetParLimits(parameterIndices[cbo_Phi_phi_string], -3*pi, 3*pi);
          setDiffParLimits(TmethodFitFunc, {cbo_Phi_amp_string, cbo_Phi_phi_string});
          histToFit->Fit(TmethodFitFunc, "QR");
        }

/////////////////////////////////////////////////////////////////////////////////////

        // fit VW
        if(includeVW){
          fixAllParameters(TmethodFitFunc);

          if(includeChangingVWFreq){ if(isParamFree(vw_w_string)) TmethodFitFunc->SetParLimits(parameterIndices[vw_w_string], 0, 0.03); }
          else{ if(isParamFree(vw_w_string)) TmethodFitFunc->SetParLimits(parameterIndices[vw_w_string], .9*TmethodFitFunc->GetParameter(parameterIndices[vw_w_string]), 1.1*TmethodFitFunc->GetParameter(parameterIndices[vw_w_string])); }          

          if(isParamFree(vw_phi_string)) TmethodFitFunc->SetParLimits(parameterIndices[vw_phi_string], -2*pi, 2*pi);
          setDiffParLimits(TmethodFitFunc, {vw_w_string, vw_phi_string});

          histToFit->Fit(TmethodFitFunc, "QR"); // fit to VW phase and frequency

          fixAllParameters(TmethodFitFunc);

          if(isParamFree(vw_tau_string)) TmethodFitFunc->SetParLimits(parameterIndices[vw_tau_string], 0, 75); // 50
          if(isParamFree(vw_amp_string)) TmethodFitFunc->SetParLimits(parameterIndices[vw_amp_string], 0, 0.1);
          setDiffParLimits(TmethodFitFunc, {vw_tau_string, vw_amp_string});

          histToFit->Fit(TmethodFitFunc, "QR"); // fit to VW amp and lifetime

          if(includeChangingVWFreq){ if(isParamFree(vw_w_string)) TmethodFitFunc->SetParLimits(parameterIndices[vw_w_string], 0, 0.03); }
          else{ if(isParamFree(vw_w_string)) TmethodFitFunc->SetParLimits(parameterIndices[vw_w_string], .9*TmethodFitFunc->GetParameter(parameterIndices[vw_w_string]), 1.1*TmethodFitFunc->GetParameter(parameterIndices[vw_w_string])); }          
          
          if(isParamFree(vw_phi_string)) TmethodFitFunc->SetParLimits(parameterIndices[vw_phi_string], -2*pi, 2*pi);
          setDiffParLimits(TmethodFitFunc, {vw_w_string, vw_phi_string, vw_tau_string, vw_amp_string});

          histToFit->Fit(TmethodFitFunc, "QR"); // fit to VW parameters
        }

/////////////////////////////////////////////////////////////////////////////////////

  fitAllTMethodParameters(TmethodFitFunc, histToFit, TmethodFitFunc->GetParameters());
  return TmethodFitFunc->GetParameters();

}


double* TmethodFit::fitAllTMethodParameters(TF1* TmethodFitFunc, TH1F* histToFit, double* pars){

        TFitResultPtr fitResult;
        int fitStatus = -1;

        TmethodFitFunc->SetParameters(pars);

        // explicitly list out all parameter limits even though some of the lines of code are redundant

          if(isParamFree(N_string)) TmethodFitFunc->ReleaseParameter(parameterIndices[N_string]);
          // if(isParamFree(tau_string)) TmethodFitFunc->SetParLimits(parameterIndices[tau_string], 64.2, 64.6);
          // if(isParamFree(A_string)) TmethodFitFunc->SetParLimits(parameterIndices[A_string], 0.3, 0.5);
          // if(isParamFree(R_string)) TmethodFitFunc->SetParLimits(parameterIndices[R_string], -100, 100);
          // if(isParamFree(phi_string)) TmethodFitFunc->SetParLimits(parameterIndices[phi_string], 2, 2.5);
          if(isParamFree(tau_string)) TmethodFitFunc->ReleaseParameter(parameterIndices[tau_string]);
          if(isParamFree(A_string)) TmethodFitFunc->ReleaseParameter(parameterIndices[A_string]);
          if(isParamFree(R_string)) TmethodFitFunc->ReleaseParameter(parameterIndices[R_string]);
          if(isParamFree(phi_string)) TmethodFitFunc->ReleaseParameter(parameterIndices[phi_string]);

          if(isParamFree(cbo_w_string)) TmethodFitFunc->SetParLimits(parameterIndices[cbo_w_string], .9*TmethodFitFunc->GetParameter(parameterIndices[cbo_w_string]), 1.1*TmethodFitFunc->GetParameter(parameterIndices[cbo_w_string]));
          if(isParamFree(cbo_tau_string)) TmethodFitFunc->SetParLimits(parameterIndices[cbo_tau_string], 10, 500);
          if(isParamFree(cbo_N_amp_string)) TmethodFitFunc->SetParLimits(parameterIndices[cbo_N_amp_string], 0, .05);
          TmethodFitFunc->FixParameter(parameterIndices[cbo_N_phi_string], fmod(TmethodFitFunc->GetParameter(parameterIndices[cbo_N_phi_string]), 2*pi));          
          if(isParamFree(cbo_N_phi_string)) TmethodFitFunc->SetParLimits(parameterIndices[cbo_N_phi_string], -3*pi, 3*pi);

          // 2 cbo parameters
          if(include2CBO){
            if(isParamFree(cbo_N2_amp_string)) TmethodFitFunc->SetParLimits(parameterIndices[cbo_N2_amp_string], 0, .01);
            TmethodFitFunc->FixParameter(parameterIndices[cbo_N2_phi_string], fmod(TmethodFitFunc->GetParameter(parameterIndices[cbo_N2_phi_string]), 2*pi));
            if(isParamFree(cbo_N2_phi_string)) TmethodFitFunc->SetParLimits(parameterIndices[cbo_N2_phi_string], -3*pi, 3*pi);
          }

          // A cbo parameters
          if(includeACBO){
            if(isParamFree(cbo_A_amp_string)) TmethodFitFunc->SetParLimits(parameterIndices[cbo_A_amp_string], 0, .02);
            TmethodFitFunc->FixParameter(parameterIndices[cbo_A_phi_string], fmod(TmethodFitFunc->GetParameter(parameterIndices[cbo_A_phi_string]), 2*pi));
            if(isParamFree(cbo_A_phi_string)) TmethodFitFunc->SetParLimits(parameterIndices[cbo_A_phi_string], -3*pi, 3*pi);
          }

          // phi cbo parameters
          if(includePhiCBO){
            if(isParamFree(cbo_Phi_amp_string)) TmethodFitFunc->SetParLimits(parameterIndices[cbo_Phi_amp_string], 0, .01);
            TmethodFitFunc->FixParameter(parameterIndices[cbo_Phi_phi_string], fmod(TmethodFitFunc->GetParameter(parameterIndices[cbo_Phi_phi_string]), 2*pi));
            if(isParamFree(cbo_Phi_phi_string)) TmethodFitFunc->SetParLimits(parameterIndices[cbo_Phi_phi_string], -3*pi, 3*pi);
          }

          if(includeVW){
            if(includeChangingVWFreq){ if(isParamFree(vw_w_string)) TmethodFitFunc->SetParLimits(parameterIndices[vw_w_string], 0, 0.03); }
            else{ if(isParamFree(vw_w_string)) TmethodFitFunc->SetParLimits(parameterIndices[vw_w_string], .9*TmethodFitFunc->GetParameter(parameterIndices[vw_w_string]), 1.1*TmethodFitFunc->GetParameter(parameterIndices[vw_w_string])); }          
            if(isParamFree(vw_tau_string)) TmethodFitFunc->SetParLimits(parameterIndices[vw_tau_string], 0, 75); // 50
            if(isParamFree(vw_amp_string)) TmethodFitFunc->SetParLimits(parameterIndices[vw_amp_string], 0, 0.1);
            TmethodFitFunc->FixParameter(parameterIndices[vw_phi_string], fmod(TmethodFitFunc->GetParameter(parameterIndices[vw_phi_string]), 2*pi));
            if(isParamFree(vw_phi_string)) TmethodFitFunc->SetParLimits(parameterIndices[vw_phi_string], -2*pi, 2*pi);
          }

          if(includeLM && isParamFree(kLoss_string)){ TmethodFitFunc->SetParLimits(parameterIndices[kLoss_string], -35, 35); }

          if(includeBR){
            if(isParamFree(br_tau_string)) TmethodFitFunc->SetParLimits(parameterIndices[br_tau_string], 0, 100);
            if(isParamFree(br_amp_string)) TmethodFitFunc->SetParLimits(parameterIndices[br_amp_string], -.05, .05);
          }


          // TmethodFitFunc->SetParLimits(parameterIndices["C"], -0.4, 0.4);

/////////////////////////////////////////////////////////////////////////////////////

          // set updated parameter limits for all parameters
          for (uint i = 0; i < updatedParLimits.size(); ++i) if(isParamFree(updatedParLimits.at(i).first)) TmethodFitFunc->SetParLimits(parameterIndices[updatedParLimits.at(i).first], updatedParLimits.at(i).second.first, updatedParLimits.at(i).second.second);

/////////////////////////////////////////////////////////////////////////////////////

        // printParamsAndLimits(TmethodFitFunc);

        fitResult = histToFit->Fit(TmethodFitFunc, TMethod_lastFitOptions.c_str());
        // fitResult = histToFit->Fit(TmethodFitFunc, (TMethod_lastFitOptions + "P").c_str()); // do this to include Pearson errors in the chi2 fit - could write something so that the option P could be included in the config, and it checks if that option exists here

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

        fitStatus = fitResult;

        // if(fitStatus != 0)
        // {
        //   cout << "T method fit returned non-zero fit status. Exiting with fit status: " << fitStatus << endl;
        //   exit(1);
        // }

        fitResult->GetCorrelationMatrix().Write("TMethodFitCorrMatrix"); // write the correlation matrix to the file

    return TmethodFitFunc->GetParameters();

}

#endif //! TMETHODFIT_HH
