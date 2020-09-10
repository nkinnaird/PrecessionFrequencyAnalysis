// 3-31-20: Older config file used with the ClusterTreeToHistsPileup.C macro.

#ifndef RATIOHISTOGRAMCONFIG_HH
#define RATIOHISTOGRAMCONFIG_HH

// config file for what conditions and the like are being used when constructing the histograms in ClusterTreeToHistsPileup.C or ClusterTreeToHistsOrderec.C

/////////////////////////////////////////////////////////////////////////////////////
	
    // random seeds in the histogram making - set to 0 before running on the grid
    static const bool useHashAsSeed = true;
	static const int histRandSeed1 = 0;
    static const int gRandomSeed = 0;
    static const bool varyRandSeeds = false; // boolean for filling histograms with varying random seeds - each iteration will have a different random seed

    static const bool perFillRandomization = false; // do randomization at the fill level instead of the cluster level - doesn't work yet with varying the random seeds

    static const bool randomize_fc = true; // randomize out fc component
    static const bool randomizeByDefaultBinWidth = true; // set to true to randomize by the default bin width of 149.2 ns, set to false to randomize by the choice of bin width

    static const bool randomize_VW = true; // randomize out VW component - can also be used to randomize out 2fy or anything else

/////////////////////////////////////////////////////////////////////////////////////

    static double nVal;
    void setnVal(int datasetCase){
        switch(datasetCase){
            case 1: nVal = 0.108; // 60h
            break;
            case 2: nVal = 0.120; // 9d
            break;
            case 3: nVal = 0.108; // Endgame
            break;
            case 4: nVal = 0.120; // HighKick
            break;
            default: cout << "bad dataset definition" << endl;
            exit(-1);
        }
    }

/////////////////////////////////////////////////////////////////////////////////////

    static const uint totalIterations = 1; // number of histogram iterations
    static const bool reduceMemory = true; // set to false to produce many 2D time energy histograms - default to true because they take up a lot of memory

    static const bool bunchNumScan = false; // if true will overwrite totalIterations to 9
    static const bool crystalRowScan = false; // if true will overwrite totalIterations to 7

/////////////////////////////////////////////////////////////////////////////////////

    static const double binWidthStart = defaultBinWidth;
    static const double binWidthStep = 0;

    static const double binEdgeShiftStart = 0;
    static const double binEdgeShiftStep = 0;

/////////////////////////////////////////////////////////////////////////////////////

    static const double gm2PeriodPPMStart = 0;
    static const double gm2PeriodPPMStep = 0;

    static const double weightingLifetimeStart = defaultLifetime;
    static const double weightingLifetimeStep = 0;

    static const double energyThresholdStart = 1700;
    static const double energyThresholdStep = 0;

/////////////////////////////////////////////////////////////////////////////////////

    static const double ADTStart = 5.;
    static const double ADTStep = 0;

    static const double SDTStart = ADTStart;
    static const double SDTStep = 0;

    static const double SGTStart = 2. * ADTStart;
    static const double SGTStep = 0;

    static const bool performTripletCorrection = false; // currently doesn't work quite right - would need updates and improvement
    static const double SDT2Start = 0;
    static const double SDT2Step = 0;

    static const double pileupTimeShiftStart = 0;
    static const double pileupTimeShiftStep = 0;

    static const double pileupEnergyScalingStart = 1;
    static const double pileupEnergyScalingStep = 0;

/////////////////////////////////////////////////////////////////////////////////////

    // my global gain correction - hasn't been used in a while and not really the right way to do things, but it could be nice in the future for faster scans
    // static const double gainSagAmpStart = 0; // 0.03 (had offsets) // eyeballed average // -0.1 - 21 iterations
    // static const double gainSagAmpStep = 0; // 0.01
    // static const double gainSagTauStart = 6700.; // eyeballed average // 1000 - 15 iterations
    // static const double gainSagTauStep = 0; // 1000

    // crystal hit gain correction factors

    static const bool recorrectGainAtCrystals = false;

    static const double crystalGainAmpFactorStart = 1;
    static const double crystalGainAmpFactorStep = 0;

    static const double crystalGainTauFactorStart = 1;
    static const double crystalGainTauFactorStep = 0;

/////////////////////////////////////////////////////////////////////////////////////

    static const double timeCutLow = 25000;
    static const double timeCutHigh = 655000; // 660000;

    static const uint runCutLow = 1; // avoid setting to 0 so it doesn't throw warnings
    static const uint runCutHigh = 1e8; // some high number for including all runs

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

    // lost muons cuts

    static const double timeCutLowMuons = 15000;

    static const double deltaTCutLow = 5;
    static const double deltaTCutLowStep = 0;
    static const double deltaTCutHigh = 7.5;
    static const double deltaTCutHighStep = 0;

    static const double energyCutLow = 100;
    static const double energyCutLowStep = 0;
    static const double energyCutHigh = 250;
    static const double energyCutHighStep = 0;

    static const double clusterSizeCut = 3;
    static const double clusterEFracCut = 0.8;

/////////////////////////////////////////////////////////////////////////////////////

#endif //! RATIOHISTOGRAMCONFIG_HH
