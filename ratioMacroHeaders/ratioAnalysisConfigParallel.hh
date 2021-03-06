// config file for what fit conditions and the like are being used in the analysis
// 3-31-20: This config file has Parallel in the name since it is the config file used with ratioMacroParallel.C
// See notes at the top of that macro.

#ifndef RATIOANALYSISCONFIGPARALLEL_HH
#define RATIOANALYSISCONFIGPARALLEL_HH

/////////////////////////////////////////////////////////////////////////////////////

    static const bool timeCheck = true; // print out time taken to do fits

    static const bool doTMethodFits = true;
    static const bool doRMethodFits = true;

    static const string TMethod_lastFitOptions = "QRSM";
    static const string RMethod_lastFitOptions = "QRSM";

    // if setStartingParamsFromTMethod and passParamsForward are set to true, the ratio fit fixed parameters in all fits will be fixed to the first T method fit results, otherwise the ratio fit fixed parameters are set from the most recent T method fit results each iteration
    static const bool setStartingParamsFromTMethod = true;
    static const bool passParamsForward = false; // only applies to added fits - passes the fit parameters from one fit on to the next
    static const bool setStartingParamsFromFirstRMethodFit = false; // if this is set to true, then the ratio fit parameters will be set to the first ratio fit results - only applies when passParamsForward is false

    static const bool fitAddedCalos = true;
    static const vector<int> fitIndividualCalos{};
    // static const vector<int> fitIndividualCalos{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
    // static const vector<int> fitIndividualCalos{11, 16};

/////////////////////////////////////////////////////////////////////////////////////

    static const int trackerStationModel = 12;

    static bool includeChangingCBOFreq = true;
    static bool includeChangingVWFreq = true;
    static string vw_w_string = (includeChangingVWFreq) ? "#kappa_{VW}" : "#omega_{VW}";

/////////////////////////////////////////////////////////////////////////////////////

    // adjust which parameters are fixed in the ratio fit here - not implemented for all parameters yet though minimal changes can be included to do so
    // static vector<string> ratioFit_added_fixedParameters {kLoss_string};
    static vector<string> ratioFit_added_fixedParameters {vw_w_string, vw_tau_string, kLoss_string};
    // static vector<string> ratioFit_added_fixedParameters {cbo_tau_string, vw_w_string, vw_tau_string, kLoss_string}; // highkick and 9d
    // static vector<string> ratioFit_added_fixedParameters {cbo_tau_string, vw_w_string, vw_tau_string, kLoss_string}; // endgame - fit start scan, or HighKick single fit since the CBO lifetime is small
    // static vector<string> ratioFit_added_fixedParameters {cbo_tau_string, cbo_N2_amp_string, cbo_N2_phi_string, cbo_A_amp_string, cbo_A_phi_string, cbo_Phi_amp_string, cbo_Phi_phi_string, vw_w_string, vw_tau_string, kLoss_string}; // 9d - fit start scan (and highkick)

    static const vector<string> ratioFit_calos_fixedParameters {vw_w_string, vw_tau_string, kLoss_string};
    // static const vector<string> ratioFit_calos_fixedParameters {cbo_tau_string, vw_w_string, vw_tau_string, kLoss_string};
    
    static vector<string> ratioFit_scan_fixedParameters {vw_w_string, vw_tau_string, kLoss_string};    
    // static vector<string> ratioFit_scan_fixedParameters {cbo_tau_string, vw_w_string, vw_tau_string, kLoss_string};
    // static vector<string> ratioFit_scan_fixedParameters {cbo_tau_string, kLoss_string}; // for bunch num fits
    // static vector<string> ratioFit_scan_fixedParameters {cbo_tau_string, cbo_A_phi_string, cbo_Phi_phi_string, cbo_N2_amp_string, cbo_N2_phi_string, vw_w_string, vw_tau_string, kLoss_string}; // 60h
    // static vector<string> ratioFit_scan_fixedParameters {cbo_tau_string, cbo_A_amp_string, cbo_A_phi_string, cbo_Phi_phi_string, cbo_N2_amp_string, cbo_N2_phi_string, vw_w_string, vw_tau_string, kLoss_string}; // endgame
    // static vector<string> ratioFit_scan_fixedParameters {cbo_tau_string, cbo_A_amp_string, cbo_A_phi_string, cbo_Phi_amp_string, cbo_Phi_phi_string, cbo_N2_amp_string, cbo_N2_phi_string, vw_w_string, vw_tau_string, kLoss_string}; // 9d, and highkick, or all higher order cbo params

    static const bool ratioFit_Include2CBO = false;
    static const bool ratioFit_IncludeACBO = false;
    static const bool ratioFit_IncludePhiCBO = false;
    static const bool ratioFit_IncludeVW = false;
    static const bool ratioFit_IncludeLM = true;
    static const bool ratioFit_IncludeBR = false;

/////////////////////////////////////////////////////////////////////////////////////

    static const vector<string> TmethodFit_added_fixedParameters {};
    // static const vector<string> TmethodFit_added_fixedParameters {cbo_tau_string, cbo_N2_amp_string, cbo_Phi_amp_string}; // highkick tests - fit start scans

    // static const vector<string> TmethodFit_calos_fixedParameters {cbo_w_string, vw_w_string, vw_tau_string, kLoss_string};    
    static const vector<string> TmethodFit_calos_fixedParameters {vw_w_string, vw_tau_string};
    
    static const vector<string> TmethodFit_scan_fixedParameters {vw_w_string, vw_tau_string, br_tau_string, br_amp_string}; // extra vector so that the first fit is allowed to float - only applies to fit start scan at the moment
    // static const vector<string> TmethodFit_scan_fixedParameters {vw_w_string, vw_tau_string, cbo_A_phi_string, cbo_Phi_phi_string, cbo_N2_amp_string, cbo_N2_phi_string, br_tau_string, br_amp_string}; // 60h
    // static const vector<string> TmethodFit_scan_fixedParameters {vw_w_string, vw_tau_string, cbo_A_amp_string, cbo_A_phi_string, cbo_Phi_phi_string, cbo_N2_amp_string, cbo_N2_phi_string, br_tau_string, br_amp_string}; // endgame
    // static const vector<string> TmethodFit_scan_fixedParameters {vw_w_string, vw_tau_string, cbo_A_amp_string, cbo_A_phi_string, cbo_Phi_phi_string, cbo_N2_amp_string, cbo_N2_phi_string, br_tau_string, br_amp_string}; // 9d
    // static const vector<string> TmethodFit_scan_fixedParameters {cbo_tau_string, vw_w_string, vw_tau_string, cbo_A_amp_string, cbo_A_phi_string, cbo_Phi_amp_string, cbo_Phi_phi_string, cbo_N2_amp_string, cbo_N2_phi_string, br_tau_string, br_amp_string}; // highkick
    // static const vector<string> TmethodFit_scan_fixedParameters {cbo_A_amp_string, cbo_A_phi_string, cbo_Phi_amp_string, cbo_Phi_phi_string, cbo_N2_amp_string, cbo_N2_phi_string, vw_w_string, vw_tau_string, }; // all higher order cbo params


    static const bool TmethodFit_Include2CBO = false;
    static const bool TmethodFit_IncludeACBO = false;
    static const bool TmethodFit_IncludePhiCBO = false;
    static const bool TmethodFit_IncludeVW = false;
    static const bool TmethodFit_IncludeLM = true;
    static const bool TmethodFit_IncludeBR = false;

/////////////////////////////////////////////////////////////////////////////////////

    // choose a parameter to fix and scan over - only implemented for ratio method - pass params forward should be false
    static const string parameterToScan = "nothing"; // set to some random string when not scanning the parameter
    static const double parameterScanStart = 0;
    static const double parameterScanStep = 0;
    // static const string parameterToScan = kLoss_string;
    // static const double parameterScanStart = 5.651 - 5 * 0.697; // 60h : 8.974 +- 0.338, 9d : 2.510 +- 0.170, Endgame : 2.345 +- 0.038, HighKick : 5.651 +- 0.697
    // static const double parameterScanStep = 0.697;
    // static const string parameterToScan = cbo_tau_string;
    // static const double parameterScanStart = 98;
    // static const double parameterScanStep = 0.1;

/////////////////////////////////////////////////////////////////////////////////////

    static const double defaultFitStart = 29094; // 29094; // 30287.6; // 30200; // 32525.6; half a g-2 period later
    static const double defaultFitEnd = 650064.4; // 4357 bins // 650000;

/////////////////////////////////////////////////////////////////////////////////////

	static const bool analyzeAllHistogramIters = false; // analyze multiple iterations in generated histograms, totalAnalysesPasses will then be overridden in ratioMacro.C
    static const string analyzeSpecificHistIter = "Iter0";
    static const int totalAnalysesPasses = 150.;

    static const double startTimeOfFit = defaultFitStart;
    static const double fitStartStep = 3*149.2;

    static const double endTimeOfFit = defaultFitEnd; // 215500;
    static const double fitEndStep = 0;

    static const double caloFitEndTime = 400000; // fit with a different end point for individual calos because of the lack of stats - can't scan over this yet - not sure if I really need to

    static bool useAutoPileupFactor = false; // apply automatic scaling factor to pileup from divided energy histograms, will override pileupMultiplier below
    static const bool useFixedPileupMultiplier = false; // applies the pileup scaling factor down below from each dataset
    static const double pileupMultiplierStart = 1;
    static const double pileupMultiplerStep = 0;

/////////////////////////////////////////////////////////////////////////////////////

    static const int ignoreCalos[0] = {}; // calorimeters to ignore when analyzing the added histograms
    // static const int ignoreCalos[20] = {5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24}; // keep 1-4
    // static const int ignoreCalos[20] = {2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 20, 21, 22, 23, 24}; // keep 1, 7, 13, 19

/////////////////////////////////////////////////////////////////////////////////////

    static bool useSeparateLostMuonsFile = true;
    string lostMuonsFilePath;

/////////////////////////////////////////////////////////////////////////////////////

    static const bool blindResults = true;
    // blinding strings down below for the individual datasets

/////////////////////////////////////////////////////////////////////////////////////

    // starting guesses for fit parameters based on dataset

    static string myBlindingString;

    static double fixedPileupMultiplier;

    // no starting N0, that's grabbed from an integral of the histogram
    static double startingTau; // us
    static double startingR;
    static double startingAsymmetry;
    static double startingPhase;

    static double startingCBOFreq; // rad/us
    static double startingCBOTau; // us

    static double startingCBONAmp;
    static double startingCBONPhase;

    static double startingN2CBO_amp;
    static double startingN2CBO_phi;

    static double startingACBO_amp;
    static double startingACBO_phi;

    static double startingPhiCBO_amp;
    static double startingPhiCBO_phi;

    static double startingVWFreq; // rad/us or fudge fact multiplying rad/us
    static double startingVWTau; // us
    static double startingVWAmp;
    static double startingVWPhase;

    static double startingLMAmp;

    static double startingBRTau;
    static double startingBRAmp;

    static double tM_w0;
    static double tM_Aconst;
    static double tM_Atau;
    static double tM_Bconst;
    static double tM_Btau;

    static double plot_cbo_freq; // for residual plots - doesn't affect fitting
    static double plot_VW_freq;


/////////////////////////////////////////////////////////////////////////////////////


void setDataset(int datasetCase, int trackerCase){

  switch(datasetCase){

    // 60H dataset
    case 1:

    cout << "Analyzing 60h dataset - setting config accordingly." << endl;

    lostMuonsFilePath = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/60h/LostMuons/LostMuons-60h-FP.root";
    if(useSeparateLostMuonsFile) cout << "Using separate lost muons file: " << lostMuonsFilePath << endl;

    myBlindingString = "See no evil."; // 60 hr blinding string

    fixedPileupMultiplier = 1.; // 1.03211; // Final Production 60h with ADT = SDT = 5 ns

    // no starting N0, that's grabbed from an integral of the histogram
    startingTau = defaultLifetime/1000.; // us for parameter, converted from ns
    startingR = 0;
    startingAsymmetry = 0.3637;
    startingPhase = 2.091;

    startingCBOTau = 184.2; // us

    startingCBONAmp = 0.004;
    startingCBONPhase = -2.035;

    startingN2CBO_amp = 1e-4;
    startingN2CBO_phi = 3.833;

    startingACBO_amp = 6e-4;
    startingACBO_phi = -1.077;

    startingPhiCBO_amp = 5e-4;
    startingPhiCBO_phi = -1.299;

    startingVWFreq = (includeChangingVWFreq) ? 0.01282 : 0; // fudge fact : rad/us // have not copied here the non-changing frequency starting parameters since I included the changing model
    startingVWTau = (includeChangingVWFreq) ? 22.80 : 0; // us
    startingVWAmp = (includeChangingVWFreq) ? 0.00272 : 0;
    startingVWPhase = (includeChangingVWFreq) ? 0.269 : 0;

    startingLMAmp = 8.088;

    startingBRTau = 0;
    startingBRAmp = 0;

    // https://muon.npl.washington.edu/elog/g2/Tracking+Analysis/199 - location of tracker parameters
    if(trackerCase == 12){ // Station 12
      tM_w0 = 2.3376 * 1e-3; // convert from rad/us to rad/ns
      tM_Aconst = 2.79; // rad
      tM_Atau = 61.1e3; // us to ns
      tM_Bconst = 5.63;
      tM_Btau = 6.07e3;
    }
    else if(trackerCase == 18){ // Station 18
      tM_w0 = 2.3375 * 1e-3; // convert from rad/us to rad/ns
      tM_Aconst = 2.79; // rad
      tM_Atau = 59.6e3; // us to ns
      tM_Bconst = 5.25;
      tM_Btau = 6.57e3;
    }
    else{
        cout << "bad tracker config" << endl;
        exit(-1);
    }

    startingCBOFreq = tM_w0 * 1e3; // rad/us

    plot_cbo_freq = 0.37; // for residual plots - doesn't affect fitting
    plot_VW_freq = 2.3;

    break;

/////////////////////////////////////////////////////////////////////////////////////

    // 9d dataset
    case 2:

    cout << "Analyzing 9d dataset - setting config accordingly." << endl;

    ratioFit_added_fixedParameters.push_back(cbo_tau_string); // fix cbo lifetime in 9d fits
    ratioFit_scan_fixedParameters.push_back(cbo_tau_string);

    lostMuonsFilePath = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/9d/LostMuons/LostMuons-9d-FP.root";
    if(useSeparateLostMuonsFile) cout << "Using separate lost muons file: " << lostMuonsFilePath << endl;

    myBlindingString = "Never wait for inspiration!"; // non-60h blinding string

    fixedPileupMultiplier = 1.; // 1.03387; // Thesis with with ADT = SDT = 5 ns

    // no starting N0, that's grabbed from an integral of the histogram
    startingTau = defaultLifetime/1000.; // us for parameter, converted from ns
    startingR = 0;
    startingAsymmetry = 0.3639;
    startingPhase = 2.0803;

    startingCBOTau = 225.6; // us

    startingCBONAmp = 0.0035;
    startingCBONPhase = 3.97;

    startingN2CBO_amp = 1.2e-4;
    startingN2CBO_phi = 5.511;

    startingACBO_amp = 2.2e-4;
    startingACBO_phi = 1.87;

    startingPhiCBO_amp = 7.8e-4;
    startingPhiCBO_phi = 4.66;

    startingVWFreq = (includeChangingVWFreq) ? 0.0147 : 0; // fudge fact : rad/us
    startingVWTau = (includeChangingVWFreq) ? 21.18 : 0; // us
    startingVWAmp = (includeChangingVWFreq) ? 0.00377 : 0;
    startingVWPhase = (includeChangingVWFreq) ? 1.359 : 0;

    startingLMAmp = 2.831;

    startingBRTau = 0;
    startingBRAmp = 0;

    // https://muon.npl.washington.edu/elog/g2/Tracking+Analysis/199 - location of tracker parameters
    if(trackerCase == 12){ // Station 12
      tM_w0 = 2.6094 * 1e-3; // convert from rad/us to rad/ns
      tM_Aconst = 2.80; // rad
      tM_Atau = 56.6e3; // us to ns
      tM_Bconst = 6.18;
      tM_Btau = 6.32e3;
    }
    else if(trackerCase == 18){ // Station 18
      tM_w0 = 2.6095 * 1e-3; // convert from rad/us to rad/ns
      tM_Aconst = 2.80; // rad
      tM_Atau = 57.6e3; // us to ns
      tM_Bconst = 5.72;
      tM_Btau = 6.99e3;
    }
    else{
        cout << "bad tracker config" << endl;
        exit(-1);
    }

    startingCBOFreq = tM_w0 * 1e3; // rad/us

    plot_cbo_freq = 0.415; // for residual plots - doesn't affect fitting
    plot_VW_freq = 2.04;

    break;

/////////////////////////////////////////////////////////////////////////////////////

    // Endgame dataset
    case 3:

    cout << "Analyzing Endgame dataset - setting config accordingly." << endl;

    lostMuonsFilePath = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/Endgame/LostMuons/LostMuons-Endgame-FP.root";
    if(useSeparateLostMuonsFile) cout << "Using separate lost muons file: " << lostMuonsFilePath << endl;

    myBlindingString = "Never wait for inspiration!"; // non-60h blinding string

    fixedPileupMultiplier = 1.; // 1.03819; // Thesis with with ADT = SDT = 5 ns

    // no starting N0, that's grabbed from an integral of the histogram
    startingTau = defaultLifetime/1000.; // us for parameter, converted from ns
    startingR = 0;
    startingAsymmetry = 0.3688;
    startingPhase = 2.076;

    startingCBOTau = 187.3; // us

    startingCBONAmp = 0.0034;
    startingCBONPhase = 0.239;

    startingN2CBO_amp = 1.2e-4;
    startingN2CBO_phi = 2.894;

    startingACBO_amp = 2.3e-4;
    startingACBO_phi = -2.240;

    startingPhiCBO_amp = 3.4e-4;
    startingPhiCBO_phi = 1.795;

    startingVWFreq = (includeChangingVWFreq) ? 0.0116 : 0; // fudge fact : rad/us // have not copied here the non-changing frequency starting parameters since I included the changing model
    startingVWTau = (includeChangingVWFreq) ? 21.71 : 0; // us
    startingVWAmp = (includeChangingVWFreq) ? 0.00283 : 0;
    startingVWPhase = (includeChangingVWFreq) ? 0.220 : 0;

    startingLMAmp = 2.638; // 1.709 - base cuts

    startingBRTau = 71.52;//20;
    startingBRAmp = -0.002544;//0.025;

    // https://muon.npl.washington.edu/elog/g2/Tracking+Analysis/199 - location of tracker parameters
    if(trackerCase == 12){ // Station 12
      tM_w0 = 2.3354 * 1e-3; // convert from rad/us to rad/ns
      tM_Aconst = 6.82; // rad
      tM_Atau = 78.3e3; // us to ns
      tM_Bconst = 5.42;
      tM_Btau = 6.54e3;
    }
    else if(trackerCase == 18){ // Station 18
      tM_w0 = 2.3356 * 1e-3; // convert from rad/us to rad/ns
      tM_Aconst = 6.85; // rad
      tM_Atau = 79.8e3; // us to ns
      tM_Bconst = 5.04;
      tM_Btau = 7.34e3;
    }
    else{
        cout << "bad tracker config" << endl;
        exit(-1);
    }

    startingCBOFreq = tM_w0 * 1e3; // rad/us

    plot_cbo_freq = 0.37; // for residual plots - doesn't affect fitting
    plot_VW_freq = 2.32;

    break;

/////////////////////////////////////////////////////////////////////////////////////

    // HighKick dataset
    case 4:

    cout << "Analyzing HighKick dataset - setting config accordingly." << endl;

    ratioFit_added_fixedParameters.push_back(cbo_tau_string); // fix cbo lifetime in HighKick fits
    ratioFit_scan_fixedParameters.push_back(cbo_tau_string);

    lostMuonsFilePath = "/gm2/data/users/nkinnaird/Ratio/FinalProductions/HighKick/LostMuons/LostMuons-HighKick-FP.root";
    if(useSeparateLostMuonsFile) cout << "Using separate lost muons file: " << lostMuonsFilePath << endl;

    myBlindingString = "Never wait for inspiration!"; // non-60h blinding string

    fixedPileupMultiplier = 1.; // 1.03413; // Thesis with with ADT = SDT = 5 ns

    // no starting N0, that's grabbed from an integral of the histogram
    startingTau = defaultLifetime/1000.; // us for parameter, converted from ns
    startingR = 0;
    startingAsymmetry = 0.3624;
    startingPhase = 2.083;

    startingCBOTau = 116.4; // us

    startingCBONAmp = 0.0034;
    startingCBONPhase = 3.696;

    startingN2CBO_amp = 1.5e-4;
    startingN2CBO_phi = 5.5;

    startingACBO_amp = 6e-4;
    startingACBO_phi = -2.06;

    startingPhiCBO_amp = 2.9e-4;
    startingPhiCBO_phi = 0.36;

    startingVWFreq = (includeChangingVWFreq) ? 0.0151 : 0; // fudge fact : rad/us
    startingVWTau = (includeChangingVWFreq) ? 15.46 : 0; // us
    startingVWAmp = (includeChangingVWFreq) ? 0.00790 : 0;
    startingVWPhase = (includeChangingVWFreq) ? 1.085 : 0;

    startingLMAmp = 4.012; // ? - base cuts

    startingBRTau = 0;
    startingBRAmp = 0;

    // https://muon.npl.washington.edu/elog/g2/Tracking+Analysis/199 - location of tracker parameters
    if(trackerCase == 12){ // Station 12
      tM_w0 = 2.6137 * 1e-3; // convert from rad/us to rad/ns
      tM_Aconst = 2.88; // rad
      tM_Atau = 49.2e3; // us to ns
      tM_Bconst = 6.27;
      tM_Btau = 6.18e3;
    }
    else if(trackerCase == 18){ // Station 18
      tM_w0 = 2.6134 * 1e-3; // convert from rad/us to rad/ns
      tM_Aconst = 2.93; // rad
      tM_Atau = 44.6e3; // us to ns
      tM_Bconst = 5.79;
      tM_Btau = 6.43e3;
    }
    else{
        cout << "bad tracker config" << endl;
        exit(-1);
    }

    startingCBOFreq = tM_w0 * 1e3; // rad/us

    plot_cbo_freq = 0.415; // for residual plots - doesn't affect fitting
    plot_VW_freq = 2.04;

    break;

/////////////////////////////////////////////////////////////////////////////////////

    // Run 2 stuff
    case 5:

    cout << "Analyzing Run 2 data - setting config accordingly." << endl;

    // some things here are kind of hacky for first looks

    useAutoPileupFactor = true;

    includeChangingCBOFreq = false;
    includeChangingVWFreq = false;
    vw_w_string = "#omega_{VW}";

    useSeparateLostMuonsFile = false;
    lostMuonsFilePath = "";
    if(useSeparateLostMuonsFile) cout << "Using separate lost muons file: " << lostMuonsFilePath << endl;

    myBlindingString = "Keep moving forward."; // Run 2 blinding string

    fixedPileupMultiplier = 1.;

    // no starting N0, that's grabbed from an integral of the histogram
    startingTau = defaultLifetime/1000.; // us for parameter, converted from ns
    startingR = 0;
    startingAsymmetry = 0.3688;
    startingPhase = 2.16;

    startingCBOTau = 255; // us

    startingCBONAmp = 0.003;
    startingCBONPhase = 0.3;

    startingN2CBO_amp = 1.2e-4;
    startingN2CBO_phi = 2.894;

    startingACBO_amp = 2.3e-4;
    startingACBO_phi = -2.240;

    startingPhiCBO_amp = 3.4e-4;
    startingPhiCBO_phi = 1.795;

    startingVWFreq = (includeChangingVWFreq) ? 0.0116 : 14; // fudge fact : rad/us // have not copied here the non-changing frequency starting parameters since I included the changing model
    startingVWTau = (includeChangingVWFreq) ? 21.71 : 20; // us
    startingVWAmp = (includeChangingVWFreq) ? 0.00283 : 0.001;
    startingVWPhase = (includeChangingVWFreq) ? 0.220 : 4.4;

    startingLMAmp = 2.638; // 1.709 - base cuts

    // https://muon.npl.washington.edu/elog/g2/Tracking+Analysis/199 - location of tracker parameters
    if(trackerCase == 12){ // Station 12
      tM_w0 = 2.3354 * 1e-3; // convert from rad/us to rad/ns
      tM_Aconst = 6.82; // rad
      tM_Atau = 78.3e3; // us to ns
      tM_Bconst = 5.42;
      tM_Btau = 6.54e3;
    }
    else if(trackerCase == 18){ // Station 18
      tM_w0 = 2.3356 * 1e-3; // convert from rad/us to rad/ns
      tM_Aconst = 6.85; // rad
      tM_Atau = 79.8e3; // us to ns
      tM_Bconst = 5.04;
      tM_Btau = 7.34e3;
    }
    else{
        cout << "bad tracker config" << endl;
        exit(-1);
    }

    startingCBOFreq = 2.32; // tM_w0 * 1e3; // rad/us

    plot_cbo_freq = 0.37; // for residual plots - doesn't affect fitting
    plot_VW_freq = 2.27;

    break;


/////////////////////////////////////////////////////////////////////////////////////

    default:

    cout << "bad dataset definition" << endl;
    exit(-1);

    }
}

/////////////////////////////////////////////////////////////////////////////////////

#endif //! RATIOANALYSISCONFIGPARALLEL_HH
