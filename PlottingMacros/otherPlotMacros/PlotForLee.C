#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TF2.h>
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
#include <TTree.h>
#include <THStack.h>
#include <TPaveStats.h>
#include <TVirtualHistPainter.h>

using namespace std;

const int nFunc = 5;
TF1* global_func[nFunc];
double microSecondFunction(double* x, double* p){
  return (global_func[int(p[0])]->Eval(x[0]*1000));
}

int PlotForLee(){

  TFile *inputFile = TFile::Open("/gm2/app/users/nkinnaird/RatioAnalysis/gm2Dev_v9_04_00/srcs/gm2analyses/macros/RatioMacro/Data/Run1/60hr-work/60hr-pileup/scale-0p8643/plotsExtra.root");
  if (inputFile == 0) {
    printf("Error: cannot open file\n");
    return 0;
  }

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // Get canvas and graph within it 
  TCanvas* ratioCBOFit_moduloCanvas = (TCanvas*) inputFile->Get("topDir/Modulus/ratioCBO_moduloPlot")->Clone("modulo_canvas"); // get canvas
  TGraphErrors* moduloGraph = ((TGraphErrors*) ratioCBOFit_moduloCanvas->GetPrimitive("ratioCBO_moduloGraph")); // get graph
  moduloGraph->SetTitle("Ratio Method - 60h Dataset;Time [#mus] % 100 #mus;V-U/V+U Ratio");

  /////////////////////////////////////////////////////////////////////////////////////

  // change from ns to us
  for(int pointNum = 0; pointNum < moduloGraph->GetN(); ++pointNum)  {
    double pointTime, pointRatio;
    moduloGraph->GetPoint(pointNum, pointTime, pointRatio);
    // if(pointTime > 100e3) moduloGraph->RemovePoint(pointNum);
    // else                  
    moduloGraph->SetPoint(pointNum, pointTime/1000., pointRatio);
  }
  ratioCBOFit_moduloCanvas->Draw();
  moduloGraph->GetYaxis()->SetLabelColor(0);
  moduloGraph->GetXaxis()->SetLimits(0,100);

  /////////////////////////////////////////////////////////////////////////////////////

  TList* listOfFuncs = new TList();

  TIter myList(ratioCBOFit_moduloCanvas->GetListOfPrimitives());
  while(TObject* obj = myList())  {
    // cout << "primitive name: " << obj->GetName() << endl;
    if(obj->InheritsFrom("TF1"))    {
      // cout << "primitive: " << obj->GetName() << " inherits from TF1 so removing" << endl;
      ratioCBOFit_moduloCanvas->GetListOfPrimitives()->Remove(obj);
      listOfFuncs->Add(obj);
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////
  
  int funcNum = 0;
  global_func[funcNum] = moduloGraph->GetFunction("cboRatioFitFunc");
  TF1* cboFunc_us_main = new TF1("microSecondFunction", microSecondFunction, 30, 100, 1);
  cboFunc_us_main->SetParameter(0,funcNum++);
  cboFunc_us_main->SetNpx(10000);
  cboFunc_us_main->SetLineColor(2);
  cboFunc_us_main->Draw("SAME");

  TIter funcListIter(listOfFuncs);
  while(TObject* obj = funcListIter())  {
    // cout << "function name: " << obj->GetName() << endl;
    global_func[funcNum] = (TF1*) obj;
    TF1* cboFunc_us = new TF1("microSecondFunctiontest", microSecondFunction, 0, 100, 1);
    cboFunc_us->SetParameter(0,funcNum++);
    cboFunc_us->SetNpx(10000);
    cboFunc_us->SetLineColor(2);
    cboFunc_us->Draw("SAME");
  }

  /////////////////////////////////////////////////////////////////////////////////////

  TText* t = new TText();
  t->SetTextAlign(32); // number reflect options relating to y axis somehow
  t->SetTextFont(42);
  for (int i=0;i<6;i++) {
    t->DrawText(-1, -i*2, "0"); // coordinates are xy of graph
  }

  double chi2 = moduloGraph->GetFunction("cboRatioFitFunc")->GetChisquare();
  int ndf = moduloGraph->GetFunction("cboRatioFitFunc")->GetNDF();
  double fracError = moduloGraph->GetFunction("cboRatioFitFunc")->GetParError(1);

  gPad->Update();
  double yMax = gPad->GetUymax();
  TPaveText* tpt1 = new TPaveText(10,yMax,50,yMax+1);
  tpt1->AddText(Form("#chi^{2} / ndf = %.1f / %d",chi2,ndf));
  tpt1->SetBorderSize(0);
  tpt1->SetFillStyle(-1);
  tpt1->Draw();

  TPaveText* tpt2 = new TPaveText(55,yMax,95,yMax+1);
  tpt2->AddText(Form("#delta#omega_{a}/#omega_{a} = %.3f ppm",fracError));
  tpt2->SetBorderSize(0);
  tpt2->SetFillStyle(-1);
  tpt2->Draw();

/////////////////////////////////////////////////////////////////////////////////////


  ratioCBOFit_moduloCanvas->SaveAs("ratioFit_grantPlot.png");



  return 0;

}
