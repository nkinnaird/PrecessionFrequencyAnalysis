// 3-31-20: Histograms config file for what conditions are being used when constructing the histograms in makeToyHistFromTH2.C.
// This probably needs to be dusted off and cleaned up the next time said macro is used.

#ifndef RATIOTOYHISTCONFIG_HH
#define RATIOTOYHISTCONFIG_HH

/////////////////////////////////////////////////////////////////////////////////////

    static const string inputFilePath = "/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_10_00/srcs/gm2analyses/macros/RatioMacro/ToyMC/TimeEnergyPairsHist.root";

/////////////////////////////////////////////////////////////////////////////////////
	
    // random seeds in the histogram making
	static const int toyRandSeed1 = 45111;
    static const int toyRandSeed2 = 34134;
    static const int toyRandSeed3 = 55233;

/////////////////////////////////////////////////////////////////////////////////////

    static const bool fillHitHists = true; // bool for whether to make hit histograms besides just pileup histograms

    static const int totalIters = 10; // number of histogram iterations

    static const int nPerFillPerCalo = 150;
    static const double nPts = 1e9;

/////////////////////////////////////////////////////////////////////////////////////

    static const bool monoEnergetic = false; // Choose single energy only - input file will be opened but ignored - this does some other stuff too which needs to be cleaned up

/////////////////////////////////////////////////////////////////////////////////////

    static const double gm2PeriodPPMStart = 0;
    static const double gm2PeriodPPMStep = 0;

    static const double energyThresholdStart = 2000;
    static const double energyThresholdStep = 0;

/////////////////////////////////////////////////////////////////////////////////////

    static const double artificialDeadtime = 6.;

    static const double SDTStart = artificialDeadtime;
    static const double SDTStep = 0;

    static const double SGTStart = 2. * artificialDeadtime;
    static const double SGTStep = 2;

    static const bool performTripletCorrection = false;
    static const double SDT2Start = artificialDeadtime;
    static const double SDT2Step = 0;

    static const double pileupMultiplierStart = 1;
    static const double pileupMultiplierStep = 0;

    static const double pileupTimeShiftStart = 0;
    static const double pileupTimeShiftStep = 0;

/////////////////////////////////////////////////////////////////////////////////////

    static const double gainSagAmpFactor = 0;
    static const double gainSagTau = 10000; // doesn't matter what this is if the constant out front is 0

/////////////////////////////////////////////////////////////////////////////////////

#endif //! RATIOTOYHISTCONFIG_HH
