// 3-31-20: Header file for producing residual plots from fits done with ratioMacro.C. 
// This macro will make associated FFT plots.

#ifndef RESIDUALPLOTS_HH
#define RESIDUALPLOTS_HH


class ResidualPlots{

public:

  ResidualPlots();

  TH1F* makeResidualPlots(string plotName, TObject* obj, string dirName, double minRange, double maxRange);
  double findPeaks(TH1F* input);

private:

  TH1F* rescaleFFTAxis(TH1* input, TH1* residualHistogram);
  void makeFFTCanvas(TH1F* fft);
  void drawLineOnCanv(double freq, TCanvas* canv, string label);

}; // end ResidualPlots class

ResidualPlots::ResidualPlots() {}

double ResidualPlots::findPeaks(TH1F* input){

  int numberOfPeaks = 1; // 1 for just the cbo frequency - careful with this failing

  TSpectrum* spec = new TSpectrum(numberOfPeaks);
  spec->Search(input);

  // for (int i = 0; i < spec->GetNPeaks(); ++i)
  // {
  //   cout << "peak: " << i << " with xposition: " << spec->GetPositionX()[i] << endl;
  // }

  double cboPeak = spec->GetPositionX()[0]; // f_cbo
  cboPeak = cboPeak * 1e-3 * 2 * pi; // convert to 1/ns and to w_cbo

  return cboPeak;
}

void ResidualPlots::drawLineOnCanv(double freq, TCanvas* canv, string label){

  canv->Update();

  TLine *line = new TLine(freq, canv->GetUymin(), freq, canv->GetUymax());
  line->SetLineColor(2);
  line->SetLineStyle(2);
  line->Draw();

  double labelY = 1.01 * canv->GetUymax();

  TText *t = new TText(freq, labelY, label.c_str());
  // t->SetTextAlign(22);
  t->SetTextColor(2);
  // t->SetTextFont(43);
  t->SetTextSize(.02);
  t->SetTextAngle(90);
  t->Draw();

  return;
}

void ResidualPlots::makeFFTCanvas(TH1F* fft){

  TH1F* fft_clone = (TH1F*) fft->Clone("fft_clone");
  // fft_clone->SetStats(false); // doesn't work in TBrowser

  double yMax = fft_clone->GetMaximum(); // scale y axis to get rid of 10^# text in axis label
  fft_clone->Scale(1./pow(10, int(log10(yMax))));

  auto canv = new TCanvas("FFT_Canvas", "FFT Canvas", 200, 10, 1200, 800);
  fft_clone->Draw("hist");
  canv->Update();

  // MHz
  double gm2Freq = 0.23;

  drawLineOnCanv(plot_cbo_freq, canv, "cbo");
  drawLineOnCanv(gm2Freq, canv, "g-2");
  drawLineOnCanv(plot_cbo_freq-gm2Freq, canv, "cbo - g-2");
  drawLineOnCanv(plot_cbo_freq+gm2Freq, canv, "cbo + g-2");
  drawLineOnCanv(plot_VW_freq, canv, "VW");

  canv->Write();
  delete canv;

/////////////////////////////////////////////////////////////////////////////////////

  auto zoomed_canv = new TCanvas("FFT_Canvas_Zoomed", "FFT Canvas Zoomed", 200, 10, 1200, 800);
  fft_clone->GetXaxis()->SetRangeUser(0.0, 1.0);
  fft_clone->Draw("hist");
  zoomed_canv->Update();

  drawLineOnCanv(plot_cbo_freq, zoomed_canv, "cbo");
  drawLineOnCanv(gm2Freq, zoomed_canv, "g-2");
  drawLineOnCanv(plot_cbo_freq-gm2Freq, zoomed_canv, "cbo - g-2");
  drawLineOnCanv(plot_cbo_freq+gm2Freq, zoomed_canv, "cbo + g-2");

  zoomed_canv->Write();
  delete zoomed_canv;

/////////////////////////////////////////////////////////////////////////////////////

  delete fft_clone;
  return;
}



TH1F* ResidualPlots::rescaleFFTAxis(TH1* input, TH1* residualHistogram) {

  double histWidth = residualHistogram->GetXaxis()->GetXmax()-residualHistogram->GetXaxis()->GetXmin();

  TH1F* rescaledFFT = new TH1F((string(input->GetName())+"_rescaled").c_str(), input->GetTitle(), input->GetNbinsX()/2, 0, input->GetNbinsX()/(2*histWidth*1e-3)); // convert to s and then to MHz
  rescaledFFT->GetXaxis()->SetTitle(input->GetXaxis()->GetTitle());
  rescaledFFT->GetYaxis()->SetTitle(input->GetYaxis()->GetTitle());

  for (int bin = 0; bin <= input->GetNbinsX()/2; bin++) {
    rescaledFFT->SetBinContent(bin, input->GetBinContent(bin));
    rescaledFFT->SetBinError(bin, input->GetBinError(bin));
   }

  return rescaledFFT;

}

TH1F* ResidualPlots::makeResidualPlots(string plotName, TObject* obj, string dirName, double minRange, double maxRange){

      auto aboveDir = gDirectory;
      auto residualPlots_dir = aboveDir->mkdir(dirName.c_str()); // take current directory and make a new directory to put these residual plots in
      residualPlots_dir->cd();

      bool isGraph = false;
      bool isHist = false;

      if(obj->InheritsFrom("TGraph")) isGraph = true;
      else if (obj->InheritsFrom("TH1")) isHist = true;

      TGraph* pullGraph = new TGraph();
      TGraph* residualGraph = new TGraph();
      TH1F* pullHist = new TH1F( (plotName + "_Projected_Pull").c_str(), (plotName + "_Projected_Pull").c_str(), 100, -5, 5);

/////////////////////////////////////////////////////////////////////////////////////

      double inBinWidth = 0;
      int minRangeBin = 0, maxRangeBin = 0, binsInRange = 0;
      double histMinTime = 0, histMaxTime = 0;


      if(obj->InheritsFrom("TGraph")){
        inBinWidth = ((TGraphErrors*)obj)->GetX()[2] - ((TGraphErrors*)obj)->GetX()[1];

        double binEdgeShift =  fmod(((TGraphErrors*)obj)->GetX()[1], inBinWidth) - inBinWidth/2;

        minRangeBin = int((minRange-binEdgeShift)/inBinWidth)+1;
        maxRangeBin = int((maxRange-binEdgeShift)/inBinWidth)+1;

        binsInRange = 1 + maxRangeBin - minRangeBin;

        histMinTime = (minRangeBin - 1) * inBinWidth + binEdgeShift;
        histMaxTime = maxRangeBin * inBinWidth + binEdgeShift;

      } 
      else if (obj->InheritsFrom("TH1")){
        inBinWidth = ((TH1F*)obj)->GetBinWidth(1);

        minRangeBin = ((TH1F*)obj)->FindBin(minRange);
        maxRangeBin = ((TH1F*)obj)->FindBin(maxRange);

        binsInRange = 1 + maxRangeBin - minRangeBin;

        histMinTime = ((TH1F*)obj)->GetBinLowEdge(minRangeBin);
        histMaxTime = ((TH1F*)obj)->GetBinLowEdge(maxRangeBin) + inBinWidth;
      } 

/////////////////////////////////////////////////////////////////////////////////////

      TH1F* residualHist = new TH1F( (plotName + "_residualHist").c_str(), (plotName + "_residualHist").c_str(), binsInRange, histMinTime, histMaxTime);

      double time, height, error = 0;
      double residual, pull = 0;

      if(isGraph){
        for (int i = 0; i < ((TGraphErrors*)obj)->GetN(); ++i){
          ((TGraphErrors*)obj)->GetPoint(i, time, height);
          error = ((TGraphErrors*)obj)->GetErrorY(i);

          if (time >= minRange && time <= maxRange){
            residual = ((TGraphErrors*)obj)->GetFunction( ((TGraphErrors*)obj)->GetListOfFunctions()->First()->GetName() )->Eval(time) - height;
            residualGraph->SetPoint(residualGraph->GetN(), time, residual);
            residualHist->SetBinContent(residualHist->FindBin(time), residual);

            pull = residual/error;
            pullGraph->SetPoint(pullGraph->GetN(), time, pull);
            pullHist->Fill(pull);
          }
        }
      }

      else if(isHist){
        for (int bin = 1; bin <= ((TH1F*)obj)->GetNbinsX(); ++bin){
          time = ((TH1F*)obj)->GetBinCenter(bin);
          height = ((TH1F*)obj)->GetBinContent(bin);
          error = ((TH1F*)obj)->GetBinError(bin);

          if (time >= minRange && time <= maxRange){
            residual = ((TH1F*)obj)->GetFunction( ((TH1F*)obj)->GetListOfFunctions()->First()->GetName() )->Eval(time) - height;
            residualGraph->SetPoint(residualGraph->GetN(), time, residual);
            residualHist->SetBinContent(residualHist->FindBin(time), residual);

            if(height >= 1){
              pull = residual/error;
              pullGraph->SetPoint(pullGraph->GetN(), time, pull);
              pullHist->Fill(pull);
            }
          }
        }
      }


      pullGraph->SetName((plotName + "_Pull_Plot").c_str());
      pullGraph->SetTitle((plotName + "_Pull_Plot").c_str());
      pullGraph->GetXaxis()->SetTitle("Time [ns]");
      pullGraph->GetYaxis()->SetTitle("Pull");
      pullGraph->Write();

      pullHist->GetXaxis()->SetTitle("Pull");
      pullHist->GetYaxis()->SetTitle("Events");

      residualGraph->SetName((plotName + "_Residual").c_str());
      residualGraph->SetTitle((plotName + "_Residual").c_str());
      residualGraph->GetXaxis()->SetTitle("Time [ns]");
      residualGraph->GetYaxis()->SetTitle("Residual");
      residualGraph->Write();

      residualHist->GetXaxis()->SetTitle("Time [ns]");
      residualHist->GetYaxis()->SetTitle("Residual");

/////////////////////////////////////////////////////////////////////////////////////

      TH1 *fftHistReal = 0;
      TVirtualFFT::SetTransform(0);
      fftHistReal = residualHist->FFT(fftHistReal, "MAG");
      fftHistReal->SetTitle((plotName + "_Residual_FFT").c_str());
      fftHistReal->SetName((plotName + "_fftHist").c_str());
      fftHistReal->GetXaxis()->SetTitle("Freq [MHz]"); // hasn't been converted to MHz yet but will be when rescaled
      fftHistReal->GetYaxis()->SetTitle("FFT Mag [arb.]");

      TH1F* scaledFFT = rescaleFFTAxis(fftHistReal, residualHist);
      delete fftHistReal;

/////////////////////////////////////////////////////////////////////////////////////

      makeFFTCanvas(scaledFFT);

/////////////////////////////////////////////////////////////////////////////////////

      aboveDir->cd(); // change back to above directory
      return scaledFFT;
}

#endif //! RESIDUALPLOTS_HH
