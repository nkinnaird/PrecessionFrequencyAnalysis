// 3-31-20: Header file for various useful plotting methods.

#include <TMultiGraph.h>
#include <TVirtualFFT.h>
#include <TLine.h>

#include <iostream>
#include "gm2Style.h"
#include "TROOT.h"

static vector<string> dataset_names = {"60h", "HighKick", "9d", "Endgame"};
static vector<int> dataset_colors = {1, 4, 2, 8};
static vector<int> dataset_markers = {20, 21, 22, 23};

static string datasetTagForPlots = "";

// check if output file from fitting has the dataset tag, and if so adjust the output file name and associated plot names
TNamed* applyDatasetTag(TFile* inputFile, string &outputFileName)
{
  if(inputFile->GetListOfKeys()->Contains("datasetTag")){
    TNamed* tag = (TNamed*) inputFile->Get("datasetTag");
    datasetTagForPlots = string("_") + tag->GetTitle();
    outputFileName.append(datasetTagForPlots);
    return tag;
  }
  else return 0;
}

// copied from gm2papers ROOT style file, couldn't get it to work directly
void SetGm2Style ()
{
  static TStyle* Gm2Style = 0;
  std::cout << "\nApplying g-2 style settings...\n" << std::endl ;
  if ( Gm2Style==0 ) Gm2Style = gm2Style();
  gROOT->SetStyle("gm2");
  gROOT->ForceStyle();
}

TStyle* gm2Style() 
{
  TStyle *gm2Style = new TStyle("gm2","gm2 style");

  // use plain black on white colors
  Int_t icol=0; // WHITE
  gm2Style->SetFrameBorderMode(icol);
  gm2Style->SetFrameFillColor(icol);
  gm2Style->SetCanvasBorderMode(icol);
  gm2Style->SetCanvasColor(icol);
  gm2Style->SetPadBorderMode(icol);
  gm2Style->SetPadColor(icol);
  gm2Style->SetStatColor(icol);

  // set the paper & margin sizes
  gm2Style->SetPaperSize(20,26);

  // set margin sizes
  gm2Style->SetPadTopMargin(0.05);
  gm2Style->SetPadRightMargin(0.05);
  gm2Style->SetPadBottomMargin(0.16);
  gm2Style->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  gm2Style->SetTitleXOffset(1.4);
  gm2Style->SetTitleYOffset(1.4);

  // use large fonts
  //Int_t font=72; // Helvetica italics
  Int_t font=42; // Helvetica
  Double_t tsize=0.05;
  gm2Style->SetTextFont(font);

  gm2Style->SetLegendFont(font);
  gm2Style->SetLegendTextSize(0.04);

  gm2Style->SetTextSize(tsize);
  gm2Style->SetLabelFont(font,"x");
  gm2Style->SetTitleFont(font,"x");
  gm2Style->SetLabelFont(font,"y");
  gm2Style->SetTitleFont(font,"y");
  gm2Style->SetLabelFont(font,"z");
  gm2Style->SetTitleFont(font,"z");
  
  gm2Style->SetLabelSize(tsize,"x");
  gm2Style->SetTitleSize(tsize,"x");
  gm2Style->SetLabelSize(tsize,"y");
  gm2Style->SetTitleSize(tsize,"y");
  gm2Style->SetLabelSize(tsize,"z");
  gm2Style->SetTitleSize(tsize,"z");

  // use bold lines and markers
  gm2Style->SetMarkerStyle(20);
  gm2Style->SetMarkerSize(1.2);
  gm2Style->SetHistLineWidth(2.);
  gm2Style->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars (as recommended in ATLAS figure guidelines)
  gm2Style->SetErrorX(0.0001);
  // get rid of error bar caps
  gm2Style->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  gm2Style->SetOptTitle(0);
  gm2Style->SetOptStat(0);
  gm2Style->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gm2Style->SetPadTickX(1);
  gm2Style->SetPadTickY(1);

  return gm2Style;

}


// pass in a function in units of ns and this will return it in us (not the fit parameters though, just the function to be drawn)
class microSecondFunctionClass {

  public:
    microSecondFunctionClass(){}

    double Evaluate(double* x, double* p) {
      return microSecondFunction->Eval(x[0]*1000);
    }

    void setFunction(TF1* inputFunc){
      microSecondFunction = inputFunc;
    }

  private:
    TF1* microSecondFunction;

};

// void drawLineOnCanv(double freq, TCanvas* canv, string label){
void drawLineOnCanv(double freq, TCanvas* canv){ // draws a vertical line

  canv->Update();

  TLine *line = new TLine(freq, canv->GetUymin(), freq, canv->GetUymax());
  line->SetLineColor(2);
  line->SetLineStyle(2);
  line->Draw();

  // double labelY = 1.01 * canv->GetUymax();

  // TText *t = new TText(freq, labelY, label.c_str());
  // // t->SetTextAlign(22);
  // t->SetTextColor(2);
  // // t->SetTextFont(43);
  // t->SetTextSize(.02);
  // t->SetTextAngle(90);
  // t->Draw();
}

void drawHorizontalLine(TCanvas* inputCanv, double lineY)
{
  inputCanv->Update();
  
  TLine *thisLine = new TLine(inputCanv->GetUxmin(), lineY, inputCanv->GetUxmax(), lineY);
  thisLine->SetLineColor(2);
  thisLine->SetLineStyle(2);
  thisLine->SetLineWidth(3);
  thisLine->Draw();
}


// draw vertical lines on canvas every g-2 period
void drawg2Lines(TCanvas* inputCanv, bool microseconds)
{
  // set bool microseconds to true to scale in microseconds, otherwise it draws in nanoseconds

  inputCanv->Update();

  double canvUxmin = inputCanv->GetUxmin();
  double canvUxmax = inputCanv->GetUxmax();
  double canvUymin = inputCanv->GetUymin();
  double canvUymax = inputCanv->GetUymax();

  double referenceZeroCrossing = 30287.6;

  double startingUx = referenceZeroCrossing - int(referenceZeroCrossing/g2Period) * g2Period; // closest point to 0 that is an integer multiple of g-2 periods down from the reference zero crossing
  double step = g2Period;

  if(microseconds){
    startingUx = startingUx/1000.;
    step = step/1000.;
  }

  for (double Ux = startingUx; Ux < canvUxmax; Ux = Ux + step)
  {
    if(Ux > canvUxmin){
      TLine *thisLine = new TLine(Ux, canvUymin, Ux, canvUymax);
      thisLine->SetLineStyle(2);
      thisLine->SetLineWidth(2);
      thisLine->SetLineColor(4);
      thisLine->Draw();
    }
  }
}


 // can also be used to convert from MeV to GeV
void nsTOus(TObject* obj, string xTitle)
{
	if(obj->InheritsFrom("TGraph"))
	{
      for (int point = 0; point < ((TGraph*)obj)->GetN(); ++point)
      {
      	double x, y;
      	((TGraph*)obj)->GetPoint(point, x, y);
      	((TGraph*)obj)->SetPoint(point, x/1000., y);
      }

	  ((TGraph*)obj)->GetXaxis()->SetTitle(xTitle.c_str());
	}
  else if (obj->InheritsFrom("TH1"))
  {
    ((TH1F*)obj)->GetXaxis()->Set(((TH1F*)obj)->GetXaxis()->GetNbins(), ((TH1F*)obj)->GetXaxis()->GetXmin()/1000., ((TH1F*)obj)->GetXaxis()->GetXmax()/1000.);
    ((TH1F*)obj)->ResetStats();
    ((TH1F*)obj)->GetXaxis()->SetTitle(xTitle.c_str());
  }
}

void nsTOusY(TObject* inGraph)
{
  for (int point = 0; point < ((TGraph*)inGraph)->GetN(); ++point)
  {
    double x, y;
    ((TGraph*)inGraph)->GetPoint(point, x, y);
    ((TGraph*)inGraph)->SetPoint(point, x, y/1000.);
    if(inGraph->InheritsFrom("TGraphErrors")) ((TGraphErrors*)inGraph)->SetPointError(point, ((TGraphErrors*)inGraph)->GetErrorX(point), ((TGraphErrors*)inGraph)->GetErrorY(point)/1000.);
  }
}

void normalizePhase(double &phase)
{
  while(phase < 0) phase += 2*pi;
  phase = fmod(phase, 2*pi);
}

void flipParameter(double &parameter)
{
  parameter = -parameter;
}

string removeDangerousCharacters(string inputStr)
{
  string str = inputStr;
  
  str.erase(std::remove(str.begin(), str.end(), '#'), str.end());
  str.erase(std::remove(str.begin(), str.end(), '{'), str.end());
  str.erase(std::remove(str.begin(), str.end(), '}'), str.end());

  return str;
}

string removeSpaces(string inputStr)
{
  string str = inputStr;
  str.erase(std::remove(str.begin(), str.end(), ' '), str.end());

  return str;
}


//pass in user coordinates to get out the normalized pad coordinates - useful for placing legends and such consistently
double GetNDCX(double x) {
  gPad->Update();//this is necessary!
  return (x - gPad->GetX1())/(gPad->GetX2()-gPad->GetX1());
}

double GetNDCY(double y) {
  gPad->Update();//this is necessary!
  return (y - gPad->GetY1())/(gPad->GetY2()-gPad->GetY1());
}


// create a histogram from a graph by projecting onto the Y-axis
TH1F* projectGraphToHist(TGraph* inGraph)
{
  TH1F* hist = new TH1F();

  int nPts = inGraph->GetN();

  double minVal = inGraph->GetYaxis()->GetXmin();
  double maxVal = inGraph->GetYaxis()->GetXmax();
  double range = maxVal-minVal;

  hist->SetName((string(inGraph->GetName()) + "_hist").c_str());
  hist->SetTitle(inGraph->GetTitle());
  hist->SetBins(nPts, minVal-.1*range, maxVal+.1*range);

  hist->GetXaxis()->SetTitle(inGraph->GetYaxis()->GetTitle());
  hist->GetYaxis()->SetTitle("Entries");

  double x,y;

  for (int pointNo = 0; pointNo < nPts; ++pointNo)
  {
    inGraph->GetPoint(pointNo, x, y);
    hist->Fill(y);
  }

  return hist;
}

void normalizeGraphToNominal(TGraph* inGraph)
{
  double xPoint, yPoint;
  double xDiff = 1000, xOne = 0, yOne = 0; // closest x point to x = 1, and corresponding y point

  for (int pointNo = 0; pointNo < inGraph->GetN(); ++pointNo) // find the x and y point for the nominal multiplier of 1
  {
    inGraph->GetPoint(pointNo, xPoint, yPoint);

    if(abs(xPoint - 1) < xDiff){
      xDiff = abs(xPoint - 1);
      xOne = xPoint;
      yOne = yPoint;
    }
  }

  // normalize all points to x = 1 value
  for (int pointNo = 0; pointNo < inGraph->GetN(); ++pointNo) 
  {
    inGraph->GetPoint(pointNo, xPoint, yPoint);
    inGraph->SetPoint(pointNo, xPoint, yPoint - yOne);
  }
}


void adjustGraphRanges(TObject* inGraph) // object so TGraph or TGraph errors can be passed in
{
  double xPoint, yPoint;
  double xMin = 0, xMax = 0, yMin = 0, yMax = 0;

  TGraph* graphToCheck = 0;

  if(inGraph->InheritsFrom("TMultiGraph")){
    TList* grlist = ((TMultiGraph*) inGraph)->GetListOfGraphs();
    TIter next(grlist);
    TObject *obj;
    
    obj = next();
    graphToCheck=(TGraph*)obj; // set limits from first graph in TMultiGraph
  }
  else graphToCheck = (TGraph*) inGraph;

  for (int pointNo = 0; pointNo < graphToCheck->GetN(); ++pointNo)
  {
    graphToCheck->GetPoint(pointNo, xPoint, yPoint);
    
    if(pointNo == 0)
    {
      xMin = xPoint;
      xMax = xPoint;
      yMin = yPoint;
      yMax = yPoint;
    }
    else{
      if(xPoint < xMin) xMin = xPoint;
      if(xPoint > xMax) xMax = xPoint;
      if(yPoint < yMin) yMin = yPoint;
      if(yPoint > yMax) yMax = yPoint;
    }
  }

  double xRange = xMax-xMin;
  double yRange = yMax-yMin;
  double percentAdditionalRange = 0.05;

  if(inGraph->InheritsFrom("TMultiGraph")){
    ((TMultiGraph*) inGraph)->GetXaxis()->SetLimits(xMin-percentAdditionalRange*xRange, xMax+percentAdditionalRange*xRange);
    ((TMultiGraph*) inGraph)->GetYaxis()->SetRangeUser(yMin-percentAdditionalRange*yRange, yMax+percentAdditionalRange*yRange); 
  }
  else{
    graphToCheck->GetXaxis()->SetLimits(xMin-percentAdditionalRange*xRange, xMax+percentAdditionalRange*xRange);
    graphToCheck->GetYaxis()->SetRangeUser(yMin-percentAdditionalRange*yRange, yMax+percentAdditionalRange*yRange); 
  }

}

// rescale the x axis in the FFT to remove the mirroring
TH1F* rescaleFFTAxis(TH1* inputFFT, TH1* originalHist) {

  double histWidth = originalHist->GetXaxis()->GetXmax()-originalHist->GetXaxis()->GetXmin();

  // TH1F* rescaledFFT = new TH1F(inputFFT->GetName(), inputFFT->GetTitle(), inputFFT->GetNbinsX()/2, 0, inputFFT->GetNbinsX()/(2*histWidth*1e-3)); // convert to s and then to MHz - for input histogram in nanoseconds
  TH1F* rescaledFFT = new TH1F(inputFFT->GetName(), inputFFT->GetTitle(), inputFFT->GetNbinsX()/2, 0, inputFFT->GetNbinsX()/(2*histWidth)); // convert to s and then to MHz - for input histogram in microseconds
  rescaledFFT->GetXaxis()->SetTitle(inputFFT->GetXaxis()->GetTitle());
  rescaledFFT->GetYaxis()->SetTitle(inputFFT->GetYaxis()->GetTitle());

  for (int bin = 0; bin <= inputFFT->GetNbinsX()/2; bin++) {
    rescaledFFT->SetBinContent(bin, inputFFT->GetBinContent(bin));
    rescaledFFT->SetBinError(bin, inputFFT->GetBinError(bin));
   }

  return rescaledFFT;
}


// careful - default right now is that the input object to this method is in units of microseconds
TH1F* doFFT(TObject* obj)
{
  TH1F* originalHist = 0;

  // if the input object is a graph turn it into a histogram
  if(obj->InheritsFrom("TGraph")){
        double inBinWidth = ((TGraph*)obj)->GetX()[2] - ((TGraph*)obj)->GetX()[1];
        double binEdgeShift =  fmod(((TGraph*)obj)->GetX()[1], inBinWidth) - inBinWidth/2;


        double xMin, yMin, xMax, yMax;
        ((TGraph*) obj)->GetPoint(0, xMin, yMin);
        ((TGraph*) obj)->GetPoint(((TGraph*) obj)->GetN()-1, xMax, yMax);


        double minRangeBin = int((xMin-binEdgeShift)/inBinWidth)+1;
        double maxRangeBin = int((xMax-binEdgeShift)/inBinWidth)+1;

        int binsInRange = 1 + maxRangeBin - minRangeBin;

        double histMinTime = (minRangeBin - 1) * inBinWidth + binEdgeShift;
        double histMaxTime = maxRangeBin * inBinWidth + binEdgeShift;

        originalHist = new TH1F("histname", "histtitle", binsInRange, histMinTime, histMaxTime);


        for (int pointNo = 0; pointNo < ((TGraph*) obj)->GetN(); ++pointNo)
        {
          double x,y;
          ((TGraph*) obj)->GetPoint(pointNo, x, y);
          originalHist->Fill(x, y);
        }
  }
  else if (obj->InheritsFrom("TH1")) originalHist = (TH1F*) obj;

  // return originalHist;
  

  TH1 *fftHistReal = 0;
  TVirtualFFT::SetTransform(0);
  fftHistReal = originalHist->FFT(fftHistReal, "MAG");
  fftHistReal->SetTitle("FFT");
  fftHistReal->GetXaxis()->SetTitle("Freq [MHz]"); // hasn't been converted to MHz yet but will be when rescaled
  fftHistReal->GetYaxis()->SetTitle("FFT Mag [arb.]");

  TH1F* scaledFFT = rescaleFFTAxis(fftHistReal, originalHist);
  delete fftHistReal;

  return scaledFFT;
}


