// config file for what conditions are being used when constructing histograms with ClusterTreeToHistsOrdered.C

#ifndef RATIOHISTOGRAMCONFIG_HH
#define RATIOHISTOGRAMCONFIG_HH

/////////////////////////////////////////////////////////////////////////////////////
	
    static const bool varyRandSeeds = false; // boolean for filling histograms with varying random seeds - each iteration will have a different random seed
    static const bool perFillRandomization = false; // do randomization at the fill level instead of the cluster level

    static const bool randomize_fc = true; // randomize out fc component - the width of the randomization comes directly from the bin width
    static const bool randomize_VW = true; // randomize out VW component - can also be used to randomize out 2fy or anything else

/////////////////////////////////////////////////////////////////////////////////////

    static double nVal;
    static double adHoc_w_a = defaultWa; // defaults which are overwritten based on nominal fit results without the correction
    static double adHoc_phi = 2.084;

    void setDataset(int datasetCase){
        switch(datasetCase){

            // 60h
            case 1: 
                nVal = 0.108;
                adHoc_w_a = defaultWa; // just leave these frequencies the default for now
                adHoc_phi = 2.09064;
            break;

            // 9d
            case 2: 
                nVal = 0.120; 
                adHoc_w_a = defaultWa;
                adHoc_phi = 2.08031;
            break;

            // Endgame
            case 3: 
                nVal = 0.108;
                adHoc_w_a = defaultWa;
                adHoc_phi = 2.07567;
            break;

            // HighKick
            case 4: 
                nVal = 0.120;
                adHoc_w_a = defaultWa;
                adHoc_phi = 2.08277;
            break;

            // Run 2 C
            case 5:
                nVal = 0.108;
                adHoc_w_a = defaultWa;
                adHoc_phi = 2.166;
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

    static const double binWidthStart = defaultBinWidth; // 149.19;
    static const double binWidthStep = 0;

    static const double binEdgeShiftStart = 0; // 53.62;
    static const double binEdgeShiftStep = 0;

/////////////////////////////////////////////////////////////////////////////////////

    // ratio histogram settings

    static const double gm2PeriodPPMStart = 0;
    static const double gm2PeriodPPMStep = 0;

    static const double weightingLifetimeStart = defaultLifetime;
    static const double weightingLifetimeStep = 0;

/////////////////////////////////////////////////////////////////////////////////////

    // energy threshold settings

    static double lowerEnergyThresholdStart = 1700; // 1700;
    static double upperEnergyThresholdStart = 1e9; // 1e9; // set default to some very high number

    static double energyThresholdStep = 0; // will step both lower and upper energy thresholds - will have to change this if I ever want to just step over the upper threshold

    static const bool makeEnergyBinnedHists = false; // set to true to make energy binned hists, will overwrite energy threshold parameters and totalIterations 
    static double energyBinMax = 3100;
    // totalIterations gets set to int( (energyBinMax - lowerEnergyThresoldStart) / energyThresholdStep)

    void binHistsByEnergy(bool inBool){ // has to be a method in order to use an if statement, is called in ClusterTreeToHistsOrdered.C with makeEnergyBinnedHists as the passed in parameter
        if(inBool){
            energyThresholdStep = 100;
            lowerEnergyThresholdStart = 500;
            upperEnergyThresholdStart = lowerEnergyThresholdStart + energyThresholdStep;
        }
    }

/////////////////////////////////////////////////////////////////////////////////////

    static const double unseenPileupThreshold = 0; // ignore pulses in pileup construction below this threshold, implemented in pileupUtils.hh, can make this scannable if desired

    static const double ADTStart = 0; // 5.;
    static const double ADTStep = 0;

    static const double SDTStart = 4.06; // ADTStart;
    static const double SDTStep = 0;

    static const double SGTStart = 10.; // 2. * ADTStart;
    static const double SGTStep = 0;

    static const bool performTripletCorrection = false; // currently doesn't work quite right - would need updates and improvement
    static const double SDT2Start = 0; // 4.06; // 5.;
    static const double SDT2Step = 0;

    static const bool setTimeShiftHalfGapTime = true;
    static const double pileupTimeShiftStart = 0;
    static const double pileupTimeShiftStep = 0;

    static const double pileupEnergyScalingStart = 1;
    static const double pileupEnergyScalingStep = 0;

/////////////////////////////////////////////////////////////////////////////////////

    // crystal hit gain correction factors

    static const bool recorrectGainAtCrystals = false;

    static const double crystalGainAmpFactorStart = 1;
    static const double crystalGainAmpFactorStep = 0;

    static const double crystalGainTauFactorStart = 1;
    static const double crystalGainTauFactorStep = 0;

    // ad hoc gain correction

    static const bool applyAdHocGainCorrection = false;

    static const double adHocAmplitudeStart = 6e-4; // -1e-3 - default scan start // 8.5e-4 - 60h // 1.12e-3 - Endgame // 5.1e-4 - HighKick // 8.4e-4 - 9d // 5.1e-4 - 2C
    static const double adHocAmplitudeStep = 0; // 1e-4 - default step  // usually 31 iters for tests

    static const double adHocLifetimeStart = defaultLifetime;
    static const double adHocLifetimeStep = 0;

    static const double adHocAsymmetryStart = 0.2; // 0.2 - default for all // 0.28 - 60h // 0.15 - Endgame // 0.24 - 9d or HighKick
    static const double adHocAsymmetryStep = 0; // 0.02 - default step  // usually 21 iters for tests


/////////////////////////////////////////////////////////////////////////////////////

    static const double timeCutLow = 25000;
    static const double timeCutHigh = 655000;

    static const uint runCutLow = 1; // avoid setting to 0 so it doesn't throw warnings
    static const uint runCutHigh = 1e8; // some high number for including all runs

    static const bool separateByRunGroup = false;
    static const vector<uint> runGroups = {17065, 17156, 17240, 17291, 17320, 17384, 17410, 17453, 17527}; // last entry is the upper bound on the last group

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

    // lost muons cuts - in order to use a more sophisticated lost muon spectrum, use LostMuonHistograms.C

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
