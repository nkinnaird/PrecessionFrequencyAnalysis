// 3-31-20: Macro for plots for scans over calorimeter number.

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TImage.h>
#include <sstream>
#include <TNtuple.h>
#include <TLine.h>

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;

bool saveImagesDirectly = false; // this is also set elsewhere in the file

void groupPlots(TDirectory* inDir, TF1** fits, string titleStr)
{
  auto individualDir = inDir->mkdir("Individual_Graphs_and_Hists");
  individualDir->cd();

  string xAxisTitle = "Calorimeter Number";
  string fitTitle = inDir->GetName();

  int numParams = fits[0]->GetNpar();

/////////////////////////////////////////////////////////////////////////////////////

  auto pValue_hist = new TH1F((fitTitle + "_pValue_hist").c_str(), (fitTitle + "_pValue_hist; P value; Events").c_str(), 100, 0, 1);
  auto chi2ndf_hist = new TH1F((fitTitle + "_chi2ndf_hist").c_str(), (fitTitle + "_chi2ndf_hist; #chi^{2}/ndf; Events").c_str(), 100, .8, 1.2);

  auto timeShifts = new TH1F((fitTitle + "_time_shifts").c_str(), (fitTitle + "_time_shifts; Time Shift (ns); Events").c_str(), 100, 0, 1000);

  auto chi2ndf_graph = new TGraphErrors();
    chi2ndf_graph->SetName((fitTitle + "_chi2ndf_graph").c_str());
    chi2ndf_graph->SetTitle((titleStr + " #chi^{2}/NDF Vs Calo Num").c_str());
      chi2ndf_graph->GetXaxis()->SetTitle(xAxisTitle.c_str());
      chi2ndf_graph->GetYaxis()->SetTitle("#chi^{2}/NDF");
      // chi2ndf_graph->GetYaxis()->SetTitleOffset(1.8);
      chi2ndf_graph->GetXaxis()->SetRangeUser(0, 25);

  for (int caloNum = 1; caloNum <= nCalos; ++caloNum)
  {
    pValue_hist->Fill(fits[caloNum-1]->GetProb());
    chi2ndf_hist->Fill(fits[caloNum-1]->GetChisquare()/fits[caloNum-1]->GetNDF());
    chi2ndf_graph->SetPoint(caloNum-1, caloNum, fits[caloNum-1]->GetChisquare()/fits[caloNum-1]->GetNDF());
    chi2ndf_graph->SetPointError(caloNum-1, 0, sqrt(2./fits[caloNum-1]->GetNDF()));
  }

  if(saveImagesDirectly){
    string canvasName = (fitTitle + "_Chi2NDF_Vs_Calo_Canv");
    auto chi2ndfcanv = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,500,400);
    chi2ndf_graph->GetXaxis()->SetLimits(0, 25); // needs to be done after the points are added to the graph
    chi2ndf_graph->GetYaxis()->SetRangeUser(0.9, 1.1);

      // TF1* constFunc = new TF1("constFunc", "[0]", 0, 25);
      // constFunc->FixParameter(0, 1);
      // constFunc->SetLineColor(2);
      // chi2ndf_graph->Fit("constFunc", "Q");

    chi2ndf_graph->Draw("AP");

      TLine *line = new TLine(0, 1, 25, 1);
      line->SetLineColor(2);
      line->SetLineStyle(2);
      line->Draw("SAME");

    chi2ndfcanv->SaveAs(("Images/" + canvasName + datasetTagForPlots + ".png").c_str());
  }

   chi2ndf_graph->Write();
   pValue_hist->Write();
   chi2ndf_hist->Write();

/////////////////////////////////////////////////////////////////////////////////////

  std::vector<TGraphErrors*> parameter_graphs;
  parameter_graphs.push_back(chi2ndf_graph);

  for (int parNum = 0; parNum < numParams; ++parNum)
  {
    string paramString = fits[0]->GetParName(parNum);
    string nameString = removeDangerousCharacters(paramString);

    bool phaseParam = false;
    string strAddon = "";
    if(paramString.compare("R") == 0) strAddon = " [ppm]";
    else if(paramString.find("tau") != string::npos) strAddon = " [#mus]";
    else if(paramString.find("omega") != string::npos) strAddon = " [rad/#mus]";
    else if(paramString.find("phi") != string::npos) phaseParam = true;

    auto paramGraph = new TGraphErrors();
    paramGraph->SetName((fitTitle + "_" + nameString + "_Graph").c_str());
    paramGraph->SetTitle((titleStr + " " + paramString + " Vs Calo Num").c_str());
    paramGraph->GetXaxis()->SetTitle(xAxisTitle.c_str());
    paramGraph->GetYaxis()->SetTitle((paramString + strAddon).c_str());
    // paramGraph->GetYaxis()->SetTitleOffset(2);

  	  for (int caloNum = 1; caloNum <= nCalos; ++caloNum)
      {
        double parameterValue = fits[caloNum-1]->GetParameter(parNum);
        double parameterError = fits[caloNum-1]->GetParError(parNum);

        if(phaseParam) normalizePhase(parameterValue);
        if(paramString.compare("#phi") == 0) timeShifts->Fill((pi-parameterValue)/defaultWa);

        paramGraph->SetPoint(caloNum-1, caloNum, parameterValue);
        paramGraph->SetPointError(caloNum-1, 0, parameterError);
      }

    paramGraph->Write();
    parameter_graphs.push_back(paramGraph);

    if(saveImagesDirectly){
      string canvasName = (fitTitle + "_" + nameString + "_Vs_Calo_Canv");
      auto graphCanvas = new TCanvas(canvasName.c_str(),canvasName.c_str(),200,10,500,400);
      paramGraph->GetXaxis()->SetLimits(0, 25);

      if(paramString.compare("R") == 0){
        paramGraph->Fit("pol0", "Q");
        paramGraph->GetFunction("pol0")->SetLineColor(2);
      }

      paramGraph->Draw("AP");
      graphCanvas->SaveAs(("Images/" + canvasName + datasetTagForPlots + ".png").c_str());
    }
  }

  timeShifts->Write();

/////////////////////////////////////////////////////////////////////////////////////

  inDir->cd();

  int split;
  if(parameter_graphs.size() % 4 == 0) split = 4;
  else split = 6;

  int numCanvases = 1 + int((numParams)/split);

  for (int canvasNum = 0; canvasNum < numCanvases; ++canvasNum)
  {
    string canvasName = fitTitle + Form("_Pars_Vs_CaloNum_%d",canvasNum+1);
    auto canvas_of_calo_plots = new TCanvas(canvasName.c_str(),"Canvas",200,10,(500 * split/2),800);
    canvas_of_calo_plots->Divide(split/2, 2);

    for (int padNum = 0; padNum < split; ++padNum)
    {
      canvas_of_calo_plots->cd(padNum+1);

      uint graphIndex = canvasNum*split + padNum;
      if(graphIndex >= parameter_graphs.size()) break;

      parameter_graphs.at(graphIndex)->Draw("AP");
    }

    canvas_of_calo_plots->Write();
    delete canvas_of_calo_plots;
  }

/////////////////////////////////////////////////////////////////////////////////////

  return;
}


int perCaloPlots(std::string filePath)
{
  gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  // gStyle->SetPadRightMargin(.05);
  // gStyle->SetOptFit(111);
  // gStyle->SetOptTitle(0);

  SetGm2Style();

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  // check if dataset tag exists in file and if so append it to the file name and write it to the file

  string outputFileName = "perCaloPlots";
  TNamed* tag = applyDatasetTag(inputFile, outputFileName);

  TFile* outputFile = new TFile((outputFileName + ".root").c_str(),"RECREATE");
  if(tag) tag->Write();
  auto topDir = outputFile->mkdir("topDir");

/////////////////////////////////////////////////////////////////////////////////////

  TNtuple *firstTuple = (TNtuple*)inputFile->Get(Form("topDir/FitPasses/FitPass0/FitConditions0"));
  float tP;
  firstTuple->SetBranchAddress("totalPasses", &tP);
  firstTuple->GetEntry(0);


  for (int fitPass = 0; fitPass < tP; ++fitPass)
  {
    auto fitPassDir = topDir->mkdir(Form("FitPass%d", fitPass));

  	TF1* fiveParamFits[nCalos];
    TF1* TMethodFits[nCalos];
  	TF1* threeParRatioFits[nCalos];
  	TF1* fullRatioFits[nCalos];

	  for (int caloNum = 1; caloNum <= nCalos; ++caloNum)
	  {
      // TF1* fiveParamFitFunction = (TF1*) ((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/Calos/Calo%d/Hist/Calo%d_Threshold_Clone", fitPass, caloNum, caloNum)))->GetFunction("fiveParamFit")->Clone();
      TF1* TMethodFitFunction = (TF1*) ((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/Calos/Calo%d/TMethod/Calo%d_Threshold_Clone", fitPass, caloNum, caloNum)))->GetFunction("TmethodFitFunc")->Clone();
	    // TF1* threeParRatioFitFunction = (TF1*) ((TGraphErrors*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/Calos/Calo%d/ThreeParamRatio/Ratio_TGraph_Calo_%d", fitPass, caloNum, caloNum)))->GetFunction("threeParamRatioFit")->Clone();
	    TF1* fullRatioFitFunction = (TF1*) ((TGraphErrors*) inputFile->Get(Form("topDir/FitPasses/FitPass%d/Calos/Calo%d/FullRatio/Full_Ratio_TGraph_Calo_%d", fitPass, caloNum, caloNum)))->GetFunction("fullRatioFitFunc")->Clone();

	    // fiveParamFits[caloNum-1] = fiveParamFitFunction;
      TMethodFits[caloNum-1] = TMethodFitFunction;
	    // threeParRatioFits[caloNum-1] = threeParRatioFitFunction;
	    fullRatioFits[caloNum-1] = fullRatioFitFunction;
	  }

    saveImagesDirectly = true;

  	// auto fiveParamDir = fitPassDir->mkdir("FiveParamFit", "5 Parameter Fit");
  	// groupPlots(fiveParamDir, fiveParamFits, "5 Parameter Fit");

    auto TmethodDir = fitPassDir->mkdir("TMethodFit", "T Method Fit");
    groupPlots(TmethodDir, TMethodFits, "T Method Fit");

   //  auto threeParRatioDir = fitPassDir->mkdir("ThreeParRatioFit");
   //  groupPlots(threeParRatioDir, threeParRatioFits, "3 Parameter Ratio Fit");

    auto fullRatioDir = fitPassDir->mkdir("FullRatioFit");
    groupPlots(fullRatioDir, fullRatioFits, "Full Ratio Fit");

  }


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  delete outputFile;

  return 1;

}
