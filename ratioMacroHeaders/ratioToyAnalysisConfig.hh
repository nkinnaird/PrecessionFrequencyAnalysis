// 3-31-20: Analysis config header file used with some of the Toy MC. Similar in many respects to ratioAnalysisConfig.hh but simpler.

#ifndef RATIOTOYANALYSISCONFIG_HH
#define RATIOTOYANALYSISCONFIG_HH

/////////////////////////////////////////////////////////////////////////////////////

    static const string TMethod_lastFitOptions = "RSM";
    static const string RMethod_lastFitOptions = "RSM";

    static const bool setStartingParamsFromTMethod = true;

/////////////////////////////////////////////////////////////////////////////////////

    static const bool includeChangingCBOFreq = false;
    static const bool includeChangingVWFreq = false;
    static const string vw_w_string = (includeChangingVWFreq) ? "#kappa_{VW}" : "#omega_{VW}";

/////////////////////////////////////////////////////////////////////////////////////

    static const vector<string> ratioFit_added_fixedParameters {vw_w_string, vw_tau_string, vw_phi_string};
    // static const vector<string> ratioFit_added_fixedParameters {cbo_w_string, cbo_tau_string, cbo_N_amp_string, cbo_N_phi_string, vw_w_string, vw_tau_string, vw_phi_string};

    static const bool ratioFit_Include2CBO = false;
    static const bool ratioFit_IncludeACBO = false;
    static const bool ratioFit_IncludePhiCBO = false;
    static const bool ratioFit_IncludeVW = true;
    static const bool ratioFit_IncludeLM = false;

/////////////////////////////////////////////////////////////////////////////////////

    static const vector<string> TmethodFit_added_fixedParameters {};

    static const bool TmethodFit_Include2CBO = false;
    static const bool TmethodFit_IncludeACBO = false;
    static const bool TmethodFit_IncludePhiCBO = false;
    static const bool TmethodFit_IncludeVW = true;
    static const bool TmethodFit_IncludeLM = false;

/////////////////////////////////////////////////////////////////////////////////////

    static const string parameterToScan = "nothing"; // set to some random string when not scanning the parameter

/////////////////////////////////////////////////////////////////////////////////////

    // no starting N0, that's grabbed from an integral of the histogram
    static const double startingTau = defaultLifetime/1000.; // us for parameter, converted from ns
    static const double startingR = 0; 
    static const double startingAsymmetry = 0.37;
    static const double startingPhase = 2.1;

    static const double startingCBOTau = 200; // us

    static const double startingCBONAmp = 0.004;
    static const double startingCBONPhase = -2;

    static const double startingN2CBO_amp = 1e-4;
    static const double startingN2CBO_phi = 4.363;

    static const double startingACBO_amp = 5e-4;
    static const double startingACBO_phi = -0.87;

    static const double startingPhiCBO_amp = 5e-4;
    static const double startingPhiCBO_phi = -1.675;

    static const double startingVWFreq = (includeChangingVWFreq) ? 0.01295 : 8.8*blindingWa*1e3; // fudge fact : rad/us
    static const double startingVWTau = (includeChangingVWFreq) ? 26.05 : 50; // us
    static const double startingVWAmp = (includeChangingVWFreq) ? 0.00226 : 0.05;
    static const double startingVWPhase = (includeChangingVWFreq) ? 2.747 : 1;

    static const double startingLMAmp = 6.467;

    // New tracker model
    // Station 12
    static const double trackingCBOFreq = 2.339; // rad/us
    static const double tM_w0 = 2.3389 * 1e-3; // convert from rad/us to rad/ns
    static const double tM_Aconst = 2.90; // rad
    static const double tM_Atau = 81.8e3; // us to ns
    static const double tM_Bconst = 5.12;
    static const double tM_Btau = 7.7e3;


    static const double plot_cbo_freq = 0.37; // for residual plots - doesn't affect fitting
    static const double plot_VW_freq = 2.3;

/////////////////////////////////////////////////////////////////////////////////////

#endif //! RATIOTOYANALYSISCONFIG_HH
