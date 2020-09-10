// 3-31-20: Config file used with LostMuonHistograms.C file.

#ifndef RATIOLOSTMUONSCONFIG_HH
#define RATIOLOSTMUONSCONFIG_HH

// config file for what conditions and the like are being used when constructing the histograms in LostMuonHistograms.C

/////////////////////////////////////////////////////////////////////////////////////

    // base cuts

    static const double timeCutLowMuons = 15000;

    static const double clusterSizeCut = 3;
    static const double clusterEFracCut = 0.8;

    static const double deltaTCutLow = 5; // ns
    static const double deltaTCutHigh = 7.5; // ns

    static const double energyCutLow = 100; // MeV
    static const double energyCutHigh = 250; // MeV

/////////////////////////////////////////////////////////////////////////////////////

    // additional cuts

    static const double deltaT13Cut = 14.4; // ns - for deuteron subtraction - set to high number to remove the cut

    static const double accBgCutoffLow = 2; // ns
    static const double accBgCutoffHigh = 4; // ns

    static const bool useNegativeSide = false;
    static const double negativeSideCutLow = deltaTCutLow * 2; // ns - 10 by default
    static const double negativeSideCutHigh = deltaTCutHigh * 2 - (deltaTCutHigh - deltaTCutLow); // ns - 12.5 by default

/////////////////////////////////////////////////////////////////////////////////////

#endif //! RATIOLOSTMUONSCONFIG_HH
