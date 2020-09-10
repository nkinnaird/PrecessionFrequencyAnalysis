// 3-31-20: Header file for various pileup utils. Main methods include code for pileup construction using the shadow method and modifying histogram or graph errors due to correlated pileup subtraction.

#include <THStack.h>

pair<double, double> tupleToPair(tuple<double, double, uint> tupleIn){
  pair<double, double> pairOut(get<0>(tupleIn), get<1>(tupleIn));
  return pairOut;
}

vector<pair<double, double> >* vectTupleToPair(vector<tuple<double, double, uint> >* tupleVector){
  vector<pair<double, double> >* vectorOfPairs = new vector<pair<double, double> >; // since this is newed the object/pointer will need to be deleted later
  for (uint i = 0; i < tupleVector->size(); ++i) vectorOfPairs->push_back(tupleToPair(tupleVector->at(i)));
  return vectorOfPairs;
}

double pileupAddedEnergy(vector<pair<double, double> > singlets, double energyScaling){
  
  double addedEnergy = 0;
  for (uint i = 0; i < singlets.size(); ++i) addedEnergy += singlets.at(i).second;
  
  return (energyScaling * addedEnergy);
}


double energyWeightedTime(vector<pair<double, double> > singlets, double shadowGaptime){
  
  double weightedTime = 0;
  double addedEnergy = pileupAddedEnergy(singlets, 1); // for the energy weighted time use the exact added energy with no scaling - could update this though since I've implemented the per iteration scaling ability
  
  for (uint i = 0; i < singlets.size(); ++i){
    if(i == 0) weightedTime += singlets.at(i).first * singlets.at(i).second;
    else weightedTime += (singlets.at(i).first - shadowGaptime) * singlets.at(i).second; // subtract off the gap time for shadow pulses (this isn't set up for triplets yet where the gap time would need to change)
  }

  weightedTime = weightedTime/addedEnergy;

/////////////////////////////////////////////////////////////////////////////////////

 // most energetic cluster time model (more akin to what the reconstruction does) - again not built for triplets with a different gap time
/*
  double mostEnergeticCluster_time = singlets.at(0).first;
  double mostEnergeticCluster_energy = singlets.at(0).second;
  
  for (uint i = 1; i < singlets.size(); ++i){
    double singlet_time = singlets.at(i).first - shadowGaptime;
    double singlet_energy = singlets.at(i).second;
    
    if(singlet_energy > mostEnergeticCluster_energy){
      mostEnergeticCluster_time = singlet_time;
      mostEnergeticCluster_energy = singlet_energy;
    }
  }

  return mostEnergeticCluster_time;
*/
/////////////////////////////////////////////////////////////////////////////////////
  
  // return singlets.front().first; // return the time of the first singlet
  // return singlets.back().first - shadowGaptime; // return the time of the last singlet

/////////////////////////////////////////////////////////////////////////////////////

  return weightedTime;
}

void checkForShadowPulses(uint &checkIndex, vector<pair<double, double> >* measuredHits, vector<int> &shadowIndices, vector<int> &tripletShadowIndices, double deadtime, double shadowGaptime, double secondDeadtime, double unseenPileupThreshold) 
{
  // if(shadowGaptime < 2. * deadtime || shadowGaptime < 2. * secondDeadtime )
  // {
  //   cout << "Shadow gap time is too small relative to deadtime. gap time: " << shadowGaptime << " deadtime: " << deadtime << " second deadtime: " << secondDeadtime << endl;
  //   exit(1);
  // }
  
  // these booleans are only used when the shadow dead time windows are larger than the ADT (for scans) and more than one pulse is found within the windows - which the code can't really handle
  bool alreadyFoundDoublet = false;
  bool alreadyFoundTriplet = false;

  double t1 = measuredHits->at(checkIndex).first;
  double E1 = measuredHits->at(checkIndex).second;

  uint shadowCheckIndex = checkIndex+1;
  while(shadowCheckIndex < measuredHits->size()) 
  {
    double t2 = measuredHits->at(shadowCheckIndex).first;
    double E2 = measuredHits->at(shadowCheckIndex).second;
    
    if(t2 < t1){
      cout << "Time hits should be ordered and are not - exiting. t1: " << t1 << " t2: " << t2  << endl;
      exit(1);
    } 
    else if( (t2 - t1) > (2.*shadowGaptime + deadtime + secondDeadtime)){ // past the windows
      checkIndex++;
      return;
    } 
    else if( (t2 - t1) > (2.*shadowGaptime + deadtime)){ // in the second dead-time window
      if(!alreadyFoundTriplet) tripletShadowIndices.push_back(shadowCheckIndex);
      shadowCheckIndex++;
      if(shadowCheckIndex == measuredHits->size()) checkIndex++;
      alreadyFoundTriplet = true;
      continue;
    }
    else if( (t2 - t1) > (1.*shadowGaptime + deadtime)){ // in the second gap time window
      shadowCheckIndex++;
      if(shadowCheckIndex == measuredHits->size()) checkIndex++;
      continue;
    }
    else if( (t2 - t1) >= shadowGaptime ){ // in the window
      if(!alreadyFoundDoublet && E1 > unseenPileupThreshold && E2 > unseenPileupThreshold) shadowIndices.push_back(shadowCheckIndex); // ignore pileup pulses when one is below some threshold
      shadowCheckIndex++;
      if(shadowCheckIndex == measuredHits->size()) checkIndex++;
      alreadyFoundDoublet = true;
      continue;
    } 
    else if( (t2 - t1) < shadowGaptime ){ // before the window
      shadowCheckIndex++;
      if(shadowCheckIndex == measuredHits->size()) checkIndex++;
      continue;
    } 
    else{
      cout << "Shouldn't get here, pileup if statement logic is off - exiting. Times and energies are t1: " << t1 << " t2: " << t2 << " E1: " << E1 << " E2: " << E2 << endl;
      exit(1);
    } 
  }
  
  return;
}

double calcPileupScaleFactor(TH1F* inEnergies, TH1F* inPileup)
{
   double plotRangeLow = 3000;
   double plotRangeHigh = 9000; // 5500

   double fitRangeLow = 3500; // energy range to divide over
   double fitRangeHigh = 4500; // 5500

   TH1F* energiesTh = (TH1F*) inEnergies->Clone("ThresholdEnergiesClone");
   TH1F* pileupEnergies = (TH1F*) inPileup->Clone("UnscaledPileupEnergies");


   for (int bin = 0; bin <= energiesTh->GetNbinsX()+1; ++bin) // set bin content outside of the range to be 0 so that plots are better visually
   {
   	if(energiesTh->GetBinCenter(bin) < plotRangeLow || energiesTh->GetBinCenter(bin) > plotRangeHigh)
   	{
   	  energiesTh->SetBinContent(bin, 0);
   	  pileupEnergies->SetBinContent(bin, 0);
   	  energiesTh->SetBinError(bin, 0);
   	  pileupEnergies->SetBinError(bin, 0);
   	}
   }

/////////////////////////////////////////////////////////////////////////////////////

  TH1F* dividedHist = (TH1F*) energiesTh->Clone("EnergyRatio");
  dividedHist->Divide(pileupEnergies);

  TF1* scaleFactor = new TF1("scaleFactor", "[0]", fitRangeLow, fitRangeHigh);
  scaleFactor->SetLineColor(2);
  dividedHist->Fit(scaleFactor, "RQ");

  dividedHist->SetTitle("Measured Energies / Pileup Energies");
  dividedHist->GetYaxis()->SetTitle("Ratio");


  double scaleNumber = scaleFactor->GetParameter(0);

  TH1F* scaledPileupEnergies = (TH1F*) pileupEnergies->Clone("scaledPileupEnergies");
  scaledPileupEnergies->Scale(scaleNumber);

  pileupEnergies->SetLineColor(2);
  scaledPileupEnergies->SetLineColor(2);

/////////////////////////////////////////////////////////////////////////////////////

  THStack* stack_unscaled = new THStack("stack_unscaled","Energy Comparison"); // use a THStack because histogram attributes will be preserved when opening the canvas in a TBrowser
  stack_unscaled->Add(energiesTh);
  stack_unscaled->Add(pileupEnergies);   

  auto canvas_unscaled = new TCanvas("canvas_unscaled","canvas_unscaled",200,10,1200,1000);
  stack_unscaled->Draw("nostack,hist");
  stack_unscaled->GetXaxis()->SetTitle("Energy (MeV)");
  stack_unscaled->GetYaxis()->SetTitle("Events");

   auto legend_unscaled = new TLegend(0.17,0.2,0.4,0.35);
   legend_unscaled->AddEntry(energiesTh,"Data","l");
   legend_unscaled->AddEntry(pileupEnergies,"Shadow pileup correction","l");
   legend_unscaled->SetFillStyle(0);
   legend_unscaled->SetTextSize(0.018);
   legend_unscaled->Draw();

   canvas_unscaled->Write();
   delete canvas_unscaled;

/////////////////////////////////////////////////////////////////////////////////////

  THStack* stack_scaled = new THStack("stack_scaled","Energy Comparison (scaled)"); // use a THStack because histogram attributes will be preserved when opening the canvas in a TBrowser
  stack_scaled->Add(energiesTh);
  stack_scaled->Add(scaledPileupEnergies);

  auto canvas_scaled = new TCanvas("canvas_scaled","canvas_scaled",200,10,1200,1000);
  stack_scaled->Draw("nostack,hist");
  stack_scaled->GetXaxis()->SetTitle("Energy (MeV)");
  stack_scaled->GetYaxis()->SetTitle("Events");

   auto legend_scaled = new TLegend(0.17,0.2,0.4,0.35);
   legend_scaled->AddEntry(energiesTh,"Data","l");
   legend_scaled->AddEntry(scaledPileupEnergies, Form("Shadow pileup correction * %f", scaleNumber),"l");
   legend_scaled->SetFillStyle(0);
   legend_scaled->SetTextSize(0.018);
   legend_scaled->Draw();

   canvas_scaled->Write();
   delete canvas_scaled;

/////////////////////////////////////////////////////////////////////////////////////

   delete energiesTh; // never really look at these, delete them so they aren't saved into the file for now
   delete pileupEnergies;
   delete scaledPileupEnergies;

  return scaleNumber;

 }

double calcZetaEfficiency(TDirectory* inDir, TH1F* earlyEnergies, TH1F* lateEnergies, TH1F* earlyPileupEnergies, TH1F* latePileupEnergies)
{
  // take in early and late pileup energies, make sure they cross zero at about the same bin, and then divide early pileup energies by late pileup energies, and scale the late pileup energies by either the value in the bin where the pileup crosses zero or the average value in the surrounding bins
  // crosses zero in the 2390 - 2400 bin
  // instead of using the lambda value below use this scale factor
  // also make a plot of the early energies / late energies with and without the pieup correction applied to both early and late energies



  auto efficiencyDir = inDir->mkdir("EarlyVsLateEfficiency");
  efficiencyDir->cd();

/////////////////////////////////////////////////////////////////////////////////////

  // either divide or subtract early pileup energies by late pileup energies

  // 0 - 10000
  // 1000 bins, 10 MeV bins
  // bin 239 is where it crosses 0 - in data at least - not in mc

  int binZeroCrossing = 1;
  int binLow = int(2000/earlyPileupEnergies->GetBinWidth(1));

  for (int bin = binLow; bin <= earlyPileupEnergies->GetNbinsX(); ++bin)
  {
    if(earlyPileupEnergies->GetBinContent(bin) < 0 && earlyPileupEnergies->GetBinContent(bin+1) > 0)
    {
      binZeroCrossing = bin;
      cout << "The bin where it crosses 0 is : " << binZeroCrossing << " with content: " << earlyPileupEnergies->GetBinContent(bin) << " and the next bin with content: " << earlyPileupEnergies->GetBinContent(bin+1) << endl;
      break;
    }
  }

  double testScaleFact = 0;
  for (int i = -2; i <= 2; ++i) // average over 5 bins, 50 MeV
  {
    testScaleFact += earlyEnergies->GetBinContent(binZeroCrossing+i) / lateEnergies->GetBinContent(binZeroCrossing+i);
  }
  testScaleFact /= 5.;


  cout << "average scaling fact: " << testScaleFact << endl;
  cout << "new scaling factor: " << earlyEnergies->GetBinContent(binZeroCrossing) / lateEnergies->GetBinContent(binZeroCrossing) << endl;
  cout << "new scaling factor+1: " << earlyEnergies->GetBinContent(binZeroCrossing+1) / lateEnergies->GetBinContent(binZeroCrossing+1) << endl;
  cout << "new scaling factor+2: " << earlyEnergies->GetBinContent(binZeroCrossing+2) / lateEnergies->GetBinContent(binZeroCrossing+2) << endl;
  cout << "new scaling factor-1: " << earlyEnergies->GetBinContent(binZeroCrossing-1) / lateEnergies->GetBinContent(binZeroCrossing-1) << endl;
  cout << "new scaling factor-2: " << earlyEnergies->GetBinContent(binZeroCrossing-2) / lateEnergies->GetBinContent(binZeroCrossing-2) << endl;

/////////////////////////////////////////////////////////////////////////////////////

    TH1F* uncorrectedEnergyRatio = (TH1F*) earlyEnergies->Clone("uncorrectedEnergyRatio");
    TH1F* uncorrectedEnergyRatioDenom = (TH1F*) lateEnergies->Clone("uncorrectedEnergyRatioDenom");
    uncorrectedEnergyRatioDenom->Scale(testScaleFact);
    uncorrectedEnergyRatio->Divide(uncorrectedEnergyRatioDenom);
    delete uncorrectedEnergyRatioDenom;

    TH1F* correctedEnergyRatio = (TH1F*) earlyEnergies->Clone("correctedEnergyRatio");
    TH1F* correctedEnergyRatioDenom = (TH1F*) lateEnergies->Clone("correctedEnergyRatioDenom");
    correctedEnergyRatio->Add(earlyPileupEnergies, -1);
    correctedEnergyRatioDenom->Scale(testScaleFact);
    TH1F* scaledLatePileup = (TH1F*) latePileupEnergies->Clone("scaledLatePileup");
    scaledLatePileup->Scale(testScaleFact);
    correctedEnergyRatioDenom->Add(scaledLatePileup, -1);
    correctedEnergyRatio->Divide(correctedEnergyRatioDenom);
    delete correctedEnergyRatioDenom;





/////////////////////////////////////////////////////////////////////////////////////

    TF1* efficiencyFactor = new TF1("efficiencyFactor", "[0]", 0, 10000);
    efficiencyFactor->SetLineColor(2);

    efficiencyFactor->SetRange(1000, 2000); // use the same function to fit and get lambda

    TH1F* lambdaHist = (TH1F*) earlyEnergies->Clone("lambdaHist");
    // lambdaHist->Add(earlyPileupEnergies, -1);
    lambdaHist->Divide(lateEnergies);
    lambdaHist->Fit(efficiencyFactor, "RQ");
    lambdaHist->SetTitle("Lambda");
    double lambda = efficiencyFactor->GetParameter(0);
    cout << "lambda is: " << lambda << endl;    

      // int scalePasses = 1;
      // double scaleSteps = 0.05;
      //   for (int scaleNum = 0; scaleNum < scalePasses; ++scaleNum)
      //   {
            // double lambdaScale = lambda;// + scaleSteps * (scaleNum - scalePasses/2.);
            double lambdaScale = testScaleFact;

            TH1F* energies_late_scaled = (TH1F*) lateEnergies->Clone("energies_late_scaled");
            energies_late_scaled->Scale(lambdaScale);

            // TH1F* zetaHist = (TH1F*) earlyEnergies->Clone(Form("zetaHist_%f",lambdaScale));
            TH1F* zetaHist = (TH1F*) earlyEnergies->Clone("zetaHist");
            zetaHist->Add(earlyPileupEnergies, -1);
            zetaHist->Add(energies_late_scaled, -1);

            TH1F* denom_Early_to_Late = (TH1F*) earlyEnergies->Clone("denom_Early_to_Late");
            denom_Early_to_Late->Add(energies_late_scaled, -1);

            zetaHist->Divide(denom_Early_to_Late);

            efficiencyFactor->SetRange(3500, 4500);
            zetaHist->Fit(efficiencyFactor, "RQ");
            zetaHist->GetYaxis()->SetTitle("Efficiency (decimal)");
            zetaHist->SetTitle("Zeta Inefficiency");
            // cout << "zeta is: " << efficiencyFactor->GetParameter(0) << endl;

            zetaHist->Write();
            delete zetaHist;

            delete energies_late_scaled; // never really look a these plots anyways, delete them so they aren't saved into the file
            delete denom_Early_to_Late; 
        // }

    inDir->cd();

    return efficiencyFactor->GetParameter(0);
}

TH1F* pileupErrorMultiplier(TH1F* bothErrors, TH1F* neitherErrors, TH1F* baseHistogram) // create histogram of pileup multiplier to scale errors by
{
  int nBins = baseHistogram->GetNbinsX();
  double histMinTime = baseHistogram->GetBinLowEdge(1);
  double histMaxTime = baseHistogram->GetBinLowEdge(nBins) + baseHistogram->GetBinWidth(1);

  TH1F* errorMultiplierByBin = new TH1F("errorMultiplierByBin", "errorMultiplierByBin; time (ns); Error Multiplier", nBins, histMinTime, histMaxTime);

    for (int bin = 1; bin <= errorMultiplierByBin->GetNbinsX(); ++bin)
    {
      if(baseHistogram->GetBinContent(bin) == 0) continue;

      double multip = 0;
      multip = sqrt(1 + (2.*neitherErrors->GetBinContent(bin) + 6.*bothErrors->GetBinContent(bin)) / baseHistogram->GetBinContent(bin) );
      errorMultiplierByBin->SetBinContent(bin, multip);
    }

    return errorMultiplierByBin;
}

TH1F* pileupErrorMultiplierApproximation(TH1F* bothErrors, TH1F* neitherErrors, TH1F* baseHistogram, double fitStart) // create histogram of pileup multiplier to scale errors by, but averaged over the first g-2 cycle and then extrapolated outwards
{
  int nBins = baseHistogram->GetNbinsX();
  double histBinWidth = baseHistogram->GetBinWidth(1);
  double histMinTime = baseHistogram->GetBinLowEdge(1);
  double histMaxTime = baseHistogram->GetBinLowEdge(nBins) + histBinWidth;

  TH1F* errorModHist = new TH1F("errorModHist", "errorModHist; time (ns); Error Multiplier (Approx)", nBins, histMinTime, histMaxTime); // using a histogram here for ease of use with the times

  // calculate a gamma for each bin in the first g-2 cycle and then average them

  int binsToAverageOver = int(g2Period/histBinWidth);

  double gamma = 0;

  double binStart = (fitStart/histBinWidth) + 1;
  double binStartCenter = errorModHist->GetBinCenter(binStart);

  for (int bin = binStart; bin < binStart + binsToAverageOver; ++bin) gamma += (2 * neitherErrors->GetBinContent(bin) + 6 * bothErrors->GetBinContent(bin)) / baseHistogram->GetBinContent(bin); 
  
  gamma = gamma/binsToAverageOver;

  for (int bin = binStart; bin <= nBins; ++bin)
  {
    double binCenter = errorModHist->GetBinCenter(bin);
    double sigmaMultiplier = sqrt(1 + gamma * exp(-(binCenter - binStartCenter)/defaultLifetime));
    errorModHist->SetBinContent(bin, sigmaMultiplier);
  }

  return errorModHist;
}

void modifyGraphErrors(TGraphErrors* inGraph, TH1F* errorMultiplierHist)
{  
  double histBinWidth = errorMultiplierHist->GetBinWidth(1);

  for (int i = 0; i < inGraph->GetN(); ++i)
  {
    double time, y;
    inGraph->GetPoint(i, time, y);

    int binNum = (time/histBinWidth + 1); // time and bin centers are the same, except for points later in the fill where the ratio is sometimes not formed

    inGraph->SetPointError(i, 0, inGraph->GetErrorY(i) * errorMultiplierHist->GetBinContent(binNum));
  }

  return;
}

void modifyHistErrors(TH1F* inHist, TH1F* errorMultiplierHist)
{  
  for (int binNum = 1; binNum <= inHist->GetNbinsX(); ++binNum) inHist->SetBinError(binNum, sqrt(inHist->GetBinContent(binNum)) * errorMultiplierHist->GetBinContent(binNum));
  return;
}
