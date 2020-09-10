// 3-31-20: Header file with basic analysis definitions common to many parts of the analysis.

#ifndef RATIOANALYSISDEFS_HH
#define RATIOANALYSISDEFS_HH

/////////////////////////////////////////////////////////////////////////////////////

    static const double pi = 3.14159265358979323846;

    static const int nCalos = 24;

    static const double approxMaxTime = 700000; 
    static const double defaultBinWidth = 149.2;

    static const double defaultLifetime = 64440; // ns

    static const double blindingFa = 0.2291 * 1e6 * 1/(1e9); // 0.2291 MHz converted to units of 1/ns, 0.2291 is the value used in the blinding software
    static const double blindingWa = 2*pi*blindingFa;

    static const double defaultFa = 0.2290735 * 1e6 * 1/(1e9); // 0.2290735 is the average value in column 2 of Table 15 of the E821 Final Report, for the fa values for the different run periods
    static const double defaultWa = 2*pi*defaultFa;

    static const double g2Period = 1/defaultFa;

/////////////////////////////////////////////////////////////////////////////////////

    // dataset run and fill info

    // Run 1

    static const string datasetTag_60h = "60h";
    static const int lowRunBound_60h = 15900;
    static const int highRunBound_60h = 16000;
    static const double numFills_60h = 1505362;

    static const string datasetTag_HighKick = "HK";
    static const int lowRunBound_HighKick = 16100;
    static const int highRunBound_HighKick = 16260;
    static const double numFills_HighKick = 1958449;
    
    static const string datasetTag_9d = "9d";
    static const int lowRunBound_9d = 16350;
    static const int highRunBound_9d = 16550;
    static const double numFills_9d = 3329228;

    static const string datasetTag_Endgame = "EG";
    static const int lowRunBound_Endgame = 16904;
    static const int highRunBound_Endgame = 17550;
    static const double numFills_Endgame = 7327570;

/////////////////////////////////////////////////////////////////////////////////////

    // Run 2

    static const string datasetTag_Run2_C = "2C";
    static const int lowRunBound_Run2_C = 24680;
    static const int highRunBound_Run2_C = 25870;
    static const double numFills_Run2_C = -1; // needs to be updated down the line once the final dataset is available    

/////////////////////////////////////////////////////////////////////////////////////

    // parameter name strings used in the T method and Ratio CBO fits

    static const string N_string           = "N_{0}";
    static const string tau_string         = "#tau";
    static const string A_string           = "A";
    static const string R_string           = "R";
    static const string phi_string         = "#phi";
    
    static const string cbo_w_string       = "#omega_{cbo}";
    static const string cbo_tau_string     = "#tau_{cbo}";
    static const string cbo_N_amp_string   = "A_{cbo-N}";
    static const string cbo_N_phi_string   = "#phi_{cbo-N}";
    
    static const string cbo_N2_amp_string  = "A_{2cbo-N}";
    static const string cbo_N2_phi_string  = "#phi_{2cbo-N}";
    static const string cbo_A_amp_string   = "A_{cbo-A}";
    static const string cbo_A_phi_string   = "#phi_{cbo-A}";
    static const string cbo_Phi_amp_string = "A_{cbo-#phi}";
    static const string cbo_Phi_phi_string = "#phi_{cbo-#phi}";

    // static const string vw_w_string        = "#omega_{VW}"; // defined in ratioAnalysisConfig so it can be adjusted depending on whether the frequency is changing or not
    static const string vw_tau_string      = "#tau_{VW}";
    static const string vw_amp_string      = "A_{VW}";
    static const string vw_phi_string      = "#phi_{VW}";

    static const string kLoss_string       = "#kappa_{loss}";

    static const string br_tau_string      = "#tau_{br}";
    static const string br_amp_string      = "A_{br}";

/////////////////////////////////////////////////////////////////////////////////////

  // struct for fit conditions - located here so that data and ToyMC can access it equally - might want to move it somewhere else at some point
  struct FitCondStruct{
    string directoryString;
    double fitStartTime;
    double fitEndTime;
    double pileupScaleFactor;
    pair<string, double> parScan;

    FitCondStruct() : directoryString(), fitStartTime(), fitEndTime(), pileupScaleFactor(), parScan() {}

    FitCondStruct(string dirString, double fitLowerBound, double fitUpperBound, double pScaleFactor, pair<string, double> parToScan)
    : directoryString(dirString)
    , fitStartTime(fitLowerBound)
    , fitEndTime(fitUpperBound)
    , pileupScaleFactor(pScaleFactor)
    , parScan(parToScan)
    {}
  };

/////////////////////////////////////////////////////////////////////////////////////

#endif //! RATIOANALYSISDEFS_HH
