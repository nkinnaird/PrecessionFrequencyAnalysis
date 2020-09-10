// 3-31-20: Macro for modulo plots for T method and Ratio method fits. I remember having some issues with canvases the last time I used it - in some cases canvases might be updated after they are saved, such that the image which shows up on screen isn't the same as that which is saved.

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
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TImage.h>
#include <TRandom3.h>
#include <sstream>
#include <TVirtualFFT.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TPaveStats.h>
// #include <TRatioPlot.h> // wait till root v 6_08

#include "ratioAnalysisDefs.hh"
#include "plotUtils.hh"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////

class modulusFunctionClass {

  public:
    modulusFunctionClass(){}

    double Evaluate(double* x, double* p) {
      return (modFunc->Eval(x[0]+100000*p[0]) - 2*p[1]);
    }

    void setFunction(TF1* inputFunc){
      modFunc = inputFunc;
    }

  private:
    TF1* modFunc;

};

/////////////////////////////////////////////////////////////////////////////////////

TF1* global_TMethod_Func;

double histFuncModulus(double* x, double* p)
{
  // return (global_TMethod_Func->Eval(x[0]+100000*p[0]));
  return (global_TMethod_Func->Eval(x[0]+100*p[0])); // 100 instead of 100000 since I'm now passing in a microsecond function
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


int moduloPlots(std::string filePath)
{
  // gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

  // check if dataset tag exists in file and if so append it to the file name and write it to the file

  string outputFileName = "moduloPlots";
  TNamed* tag = applyDatasetTag(inputFile, outputFileName);

  TFile* outputFile = new TFile((outputFileName + ".root").c_str(),"RECREATE");
  if(tag) tag->Write();
  auto topDir = outputFile->mkdir("topDir");

/////////////////////////////////////////////////////////////////////////////////////
  // These only get set for the interactive root session (any generated canvases, etc.), but does not apply to the output root file - that comes from .rootlogon.C
  gStyle->SetOptStat(000000);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(1111);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerColor(1);
  gStyle->SetMarkerSize(1);
  gStyle->SetLineColor(1);
  gStyle->SetStatY(0.94);

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

// create plot of 24 calo fits on a single canvas
/*
    auto Calos_24_Dir = topDir->mkdir("Calos_24");
    Calos_24_Dir->cd();

        auto caloFiveParamHists = new TCanvas("caloFiveParamHists","Canvas",200,10,1200,800);
        caloFiveParamHists->Divide(6,4);

        auto caloRatioGraphs = new TCanvas("caloRatioGraphs","Canvas",200,10,1200,800);
        caloRatioGraphs->Divide(6,4);

        auto caloRatioCBOGraphs = new TCanvas("caloRatioCBOGraphs","Canvas",200,10,1200,800);
        caloRatioCBOGraphs->Divide(6,4);
  for (int caloNum = 1; caloNum <= nCalos; ++caloNum)
  { 
      TH1F* caloTimes = ((TH1F*) inputFile->Get(("topDir/FitPasses/FitPass0/Calos/Calo" + to_string(caloNum) + "/Hist/Calo" + to_string(caloNum) + "_Threshold_Clone").c_str())->Clone(("Calo" + to_string(caloNum)).c_str()));
      caloTimes->SetTitle(Form("5 Param Hist Calo %d", caloNum));

            caloFiveParamHists->cd();
            caloFiveParamHists->cd(caloNum);
            caloTimes->GetXaxis()->SetRangeUser(50000,70000);
            caloTimes->Draw();

    TGraphErrors* ratioGraph = ((TGraphErrors*) inputFile->Get(("topDir/FitPasses/FitPass0/Calos/Calo" + to_string(caloNum) + "/Ratio/Ratio_TGraph_Calo_" + to_string(caloNum)).c_str())->Clone(("Calo" + to_string(caloNum)).c_str()));

        caloRatioGraphs->cd();
        caloRatioGraphs->cd(caloNum);
        ratioGraph->GetXaxis()->SetRangeUser(50000,70000);
        ratioGraph->Draw("AP");

    ratioGraph->GetFunction("threeParamRatioFit")->Draw("SAME");


    TGraphErrors* ratioCBOGraph = ((TGraphErrors*) inputFile->Get(("topDir/FitPasses/FitPass0/Calos/Calo" + to_string(caloNum) + "/RatioCBO/Ratio_CBO_TGraph_Calo_" + to_string(caloNum)).c_str())->Clone(("Calo" + to_string(caloNum)).c_str()));

        caloRatioCBOGraphs->cd();
        caloRatioCBOGraphs->cd(caloNum);
        ratioCBOGraph->GetXaxis()->SetRangeUser(50000,70000);
        ratioCBOGraph->Draw("AP");

    ratioCBOGraph->GetFunction("fullRatioFitFunc")->Draw("SAME");

  } // end calo for loop

  caloFiveParamHists->Write();
  caloRatioGraphs->Write();
  caloRatioCBOGraphs->Write();

  delete caloFiveParamHists;
  delete caloRatioGraphs;
  delete caloRatioCBOGraphs;
*/
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  // modulo plots

  auto microSecondClass = new microSecondFunctionClass();
  auto modulusClass = new modulusFunctionClass();

/////////////////////////////////////////////////////////////////////////////////////

  // don't really need the 3 parameter ratio modulo plot, but I define some objects in this block of code

  double xmin;
  double xmax;
  int numMod;

    auto modulusDir = topDir->mkdir("Modulus");
    modulusDir->cd();

/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

      auto fullRatio_moduloPlot = new TCanvas("fullRatio_moduloPlot","Canvas",200,10,1200,800);

      TGraphErrors* fullRatio_moduloGraph = ((TGraphErrors*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/FullRatio/Added_Times_Full_Ratio_Graph")->Clone("fullRatio_moduloGraph"));
      // TGraphErrors* fullRatio_moduloGraph = ((TGraphErrors*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/RatioCBO/Added_Times_Ratio_CBO_Graph")->Clone("fullRatio_moduloGraph"));

      TF1* fullRatio_fitFunction = (TF1*) fullRatio_moduloGraph->GetFunction("fullRatioFitFunc");
      // TF1* fullRatio_fitFunction = (TF1*) fullRatio_moduloGraph->GetFunction("cboRatioFitFunc");
      fullRatio_fitFunction->GetRange(xmin, xmax);
      numMod = int(xmax/100000) + 1;

      // remove times before and after fit for display purposes
      for (int pointNum = 0; pointNum < fullRatio_moduloGraph->GetN(); ++pointNum)
      {
        double pointTime, pointRatio;
        fullRatio_moduloGraph->GetPoint(pointNum, pointTime, pointRatio);

        if(pointTime < xmin || pointTime > xmax){
          fullRatio_moduloGraph->RemovePoint(pointNum);
          pointNum--;
          continue;
        } 

        fullRatio_moduloGraph->SetPoint(pointNum, fmod(pointTime,100000), pointRatio-2*int(pointTime/100000));
      }

      nsTOus(fullRatio_moduloGraph, "");

      fullRatio_moduloGraph->GetXaxis()->SetRangeUser(0, 100);
      fullRatio_moduloGraph->GetYaxis()->SetRangeUser(-13, 1);
      fullRatio_moduloGraph->GetXaxis()->SetTitle("Time [#mus] % 100 #mus");
      fullRatio_moduloGraph->GetYaxis()->SetTitle("Ratio Value (Counts above 1.7 GeV)");
      fullRatio_moduloGraph->Draw("AP");



      cout << "Full Ratio Fit R: " << fullRatio_fitFunction->GetParameter(1) << endl;

      microSecondClass->setFunction(fullRatio_fitFunction);
      auto ratioCBOFunc_microsecond = new TF1("ratioCBOFunc_microsecond", microSecondClass, &microSecondFunctionClass::Evaluate, fullRatio_fitFunction->GetXmin()/1000., 100, fullRatio_fitFunction->GetNpar());
      ratioCBOFunc_microsecond->SetLineColor(2);
      ratioCBOFunc_microsecond->SetNpx(10000);
      ratioCBOFunc_microsecond->DrawClone("SAME");

      modulusClass->setFunction(fullRatio_fitFunction);
      auto ratioCBOFunc_modulus = new TF1("ratioCBOFunc_modulus", modulusClass, &modulusFunctionClass::Evaluate, 0, 100000, 2);
      ratioCBOFunc_microsecond->SetRange(0, 100);

      for (int i = 1; i < numMod; ++i)
      {
        if(i == numMod-1) ratioCBOFunc_microsecond->SetRange(0, 50); // for last plotted function set the end to 50 us to correspond to fit end time of 650 us

        ratioCBOFunc_modulus->SetParameter(0, i);
        ratioCBOFunc_modulus->SetParameter(1, i);
        microSecondClass->setFunction(ratioCBOFunc_modulus);
        ratioCBOFunc_microsecond->DrawClone("SAME");
      }

      // I don't understand why but I need to paint the stats box before it exists in order to grab and redraw it - even though the stats box appears when the canvas is created
      // actually there is a different and probably better way to do this by updating the canvas and then grabbing the stats primitive directly from the canvas - SetOptFit in that case needs to be non-zero
      fullRatio_moduloGraph->PaintStats(fullRatio_fitFunction);
      TPaveStats* ratioCBO_stats = (TPaveStats*) fullRatio_moduloGraph->GetListOfFunctions()->FindObject("stats");
      ratioCBO_stats->SetBorderSize(1);
      ratioCBO_stats->Draw("SAME");

      fullRatio_moduloGraph->SetTitle("Full Ratio Fit");

      // having trouble with replacing axes labels and the like - just draw a box over it and careful with the hardcoded positions
      TBox* whiteBox = new TBox(-7.15, 1, -.2, -13);
      whiteBox->SetFillStyle(1001);
      whiteBox->SetFillColor(0);
      whiteBox->Draw("SAME");

      TText* zeroText = new TText();
      zeroText->SetTextAlign(32); // number reflect options relating to y axis somehow
      zeroText->SetTextFont(42);
      for (int i=0;i<numMod;i++) {
        zeroText->DrawText(-1, -i*2, "0"); // coordinates are xy of graph
      }

      fullRatio_moduloPlot->Modified();

      fullRatio_moduloPlot->SaveAs(("Images/fullRatio_moduloPlot" + datasetTagForPlots + ".png").c_str());
      fullRatio_moduloPlot->Write();




/*
  gStyle->SetOptFit(0);

      auto simplifiedCanv_ratio = (TCanvas*) fullRatio_moduloPlot->DrawClonePad();

      TGraphErrors *graph_copy = (TGraphErrors*)simplifiedCanv_ratio->FindObject("fullRatio_moduloGraph");
      graph_copy->GetYaxis()->SetTitleOffset(1.2);

      // simplifiedCanv_ratio->ls();

      TPaveStats *st = (TPaveStats*)simplifiedCanv_ratio->FindObject("stats");
      delete st;

      simplifiedCanv_ratio->SetRightMargin(0.05);
      simplifiedCanv_ratio->SetLeftMargin(0.1);
      simplifiedCanv_ratio->SetTopMargin(0.1);
      simplifiedCanv_ratio->SetBottomMargin(0.1);
      simplifiedCanv_ratio->Update();

      int datasetNum = 2;
      // double xTextRange = 10; // 60h
      double xTextRange = 6; // 9d
      // double xTextRange = 18; // HighKick, Endgame

      gPad->Update();
      double yMax = gPad->GetUymax();

      TPaveText* datasetName = new TPaveText(0,yMax,xTextRange,yMax+2);
      datasetName->AddText(dataset_names.at(datasetNum).c_str());
      datasetName->SetBorderSize(0);
      datasetName->SetTextColor(dataset_colors.at(datasetNum));
      datasetName->SetFillStyle(0);
      datasetName->Draw("SAME");

      double chi2 = fullRatio_fitFunction->GetChisquare();
      int ndf = fullRatio_fitFunction->GetNDF();
      double fracError = fullRatio_fitFunction->GetParError(1);

      TPaveText* tpt1 = new TPaveText(20,yMax,60,yMax+2);
      tpt1->AddText(Form("#chi^{2} / ndf = %.1f / %d",chi2,ndf));
      tpt1->SetBorderSize(0);
      tpt1->SetFillStyle(-1);
      tpt1->Draw("SAME");

      TPaveText* tpt2 = new TPaveText(65,yMax,100,yMax+2);
      tpt2->AddText(Form("#delta#omega_{a}/#omega_{a} = %.3f ppm",fracError));
      tpt2->SetBorderSize(0);
      tpt2->SetFillStyle(-1);
      tpt2->Draw("SAME");

      simplifiedCanv_ratio->Update();
      simplifiedCanv_ratio->Modified();
      simplifiedCanv_ratio->SaveAs(("Images/fullRatio_moduloPlot_noStats" + datasetTagForPlots + ".png").c_str());
*/
/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////

  // gStyle->SetOptFit(1111);
  gStyle->SetStatH(0.13);


      auto Tmethod_moduloPlot = new TCanvas("Tmethod_moduloPlot","Canvas",200,10,1200,800);

      auto moduloHist = ((TH1F*) inputFile->Get("topDir/FitPasses/FitPass0/addedDir/TMethod/allTimesAdded_TMethod")->Clone("moduloHist"));
      double histBinWidth = moduloHist->GetBinWidth(1);
      cout << "Counts in fit: " << moduloHist->Integral(xmin/histBinWidth, xmax/histBinWidth) << endl;

      // remove times before and after fit for display purposes
      for (int bin = 1; bin < moduloHist->GetNbinsX(); ++bin){
        if(moduloHist->GetBinCenter(bin) < xmin || moduloHist->GetBinCenter(bin) > xmax) moduloHist->SetBinContent(bin, 0);
      }

      int numHists = 7;
      TH1F* modHists[numHists];
      for (int i = 0; i < numHists; ++i)
      {
            modHists[i] = new TH1F(Form("modHist%d",i), "T Method Fit", 100000/histBinWidth, 0, 100000); // bins in this modulo plot won't be exactly like the non-modulated histogram, but that's fine
      }

      for (int bin = 1; bin <= moduloHist->GetNbinsX(); ++bin)
      {
        double time = moduloHist->GetBinCenter(bin);
        double moduloTime = fmod(time,100000);

        int histNum = int(time/100000);
        if(histNum >= numHists) break;

        int binCont = moduloHist->GetBinContent(bin);
        int modBin = int(moduloTime/histBinWidth)+1;
        modHists[histNum]->SetBinContent(modBin, binCont);
      }

      auto TmethodFitFunc = (TF1*) moduloHist->GetFunction("TmethodFitFunc")->Clone();
      TmethodFitFunc->GetRange(xmin, xmax);
      numMod = int(xmax/100000) + 1;

      modHists[0]->GetListOfFunctions()->Add(TmethodFitFunc); // so I can see the fit parameters in the stats box

      for (int i = 0; i < numHists; ++i)
      {
        nsTOus(modHists[i], "");
        if(i == 0) modHists[i]->Draw("E");
        else modHists[i]->Draw("E SAME");
      }

/////////////////////////////////////////////////////////////////////////////////////

      // in the ratio modulo plot I plot the microsecond function which is filled from the modulo plot, whereas here it's the reverse - might want to fix that someday

      microSecondClass->setFunction(TmethodFitFunc);
      auto TmethodFunc_microsecond = new TF1("TmethodFunc_microsecond", microSecondClass, &microSecondFunctionClass::Evaluate, TmethodFitFunc->GetXmin()/1000., TmethodFitFunc->GetXmax()/1000., TmethodFitFunc->GetNpar());
      TmethodFunc_microsecond->SetLineColor(2);
      TmethodFunc_microsecond->SetNpx(10000);
      TmethodFunc_microsecond->DrawClone("SAME");

      global_TMethod_Func = TmethodFunc_microsecond;

      for (int i = 1; i < numMod; ++i)
      {
        double fitEnd = 100; 
        if(i == numMod-1) fitEnd = 50; // for last plotted function set the end to 50 us to correspond to fit end time of 650 us

        TF1* TMethod_mod_func = new TF1("TMethod_mod_func", histFuncModulus, 0, fitEnd, 1);
        TMethod_mod_func->SetNpx(10000);
        TMethod_mod_func->SetLineColor(2);

        TMethod_mod_func->SetParameter(0, i);
        TMethod_mod_func->DrawClone("SAME");
      }

/////////////////////////////////////////////////////////////////////////////////////

      modHists[0]->GetXaxis()->SetTitle("Time [#mus] % 100 #mus");
      modHists[0]->GetYaxis()->SetTitle("Counts above 1.7 GeV");

      modHists[0]->SetMinimum(30);

      Tmethod_moduloPlot->SetLogy(); // this will be overridden when opening root files because of .rootlogon.C gRoot->ForceStyle()

      Tmethod_moduloPlot->Update();

      TPaveStats* TMethod_stats = (TPaveStats*) modHists[0]->GetListOfFunctions()->FindObject("stats");
      TMethod_stats->SetBorderSize(1);
      TMethod_stats->Draw("SAME");

      Tmethod_moduloPlot->SaveAs(("Images/TMethod_moduloPlot" + datasetTagForPlots + ".png").c_str());
      Tmethod_moduloPlot->Write();
/////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////

  return 1;
}
