// 3-31-20: Macro to print out parameters from T method and Ratio method fits. Also prints them out in a latex table style format.

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2D.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TImage.h>
#include <sstream>
#include <TNtuple.h>
#include <TPaveStats.h>
#include <TMatrixD.h>

using namespace std;

bool justPrintR = false;
bool printInLatexFormat = false;

bool printAllFits = false;
uint fitNum = 0;


int PrintParameters(std::string filePath)
{

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

/////////////////////////////////////////////////////////////////////////////////////


if(!printAllFits){

       // print number of counts in T Method fit range

      auto TmethodHist = ((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass%i/addedDir/TMethod/allTimesAdded_TMethod",fitNum))->Clone("TmethodHist"));
      double histBinWidth = TmethodHist->GetBinWidth(1);
      
      double xmin, xmax;

      TF1* Tmethod_fitFunction = (TF1*) TmethodHist->GetFunction("TmethodFitFunc");
      Tmethod_fitFunction->GetRange(xmin, xmax);

      cout << "Counts in fit: " << TmethodHist->Integral(xmin/histBinWidth, xmax/histBinWidth) << endl;


    if(inputFile->Get(Form("topDir/FitPasses/FitPass%i/addedDir/TMethod/allTimesAdded_TMethod",fitNum)) != 0)
    {
      TF1* TMethod_fit = (TF1*) ((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass%i/addedDir/TMethod/allTimesAdded_TMethod",fitNum)))->GetFunction("TmethodFitFunc");
      cout << endl << "T Method fit #chi^{2}: " << TMethod_fit->GetChisquare() << " NDF: " << TMethod_fit->GetNDF() << " #chi^{2}/NDF: " << TMethod_fit->GetChisquare()/TMethod_fit->GetNDF() << " p value: " << TMethod_fit->GetProb() << endl;

      if(justPrintR) cout << endl << "T Method fit R: " << TMethod_fit->GetParameter(3) << " error: " << TMethod_fit->GetParError(3) << endl;
      else{
        cout << endl << "T Method parameter results: " << endl << endl;
        for (int parNum = 0; parNum < TMethod_fit->GetNpar(); ++parNum) cout << TMethod_fit->GetParName(parNum) << " value: " << TMethod_fit->GetParameter(parNum) << " error: " << TMethod_fit->GetParError(parNum) << endl; 

          if(printInLatexFormat){
            cout << endl << "And in LaTex format: " << endl << endl;
            for (int parNum = 0; parNum < TMethod_fit->GetNpar(); ++parNum) cout << "$\\SI{" << TMethod_fit->GetParameter(parNum) << "}{}$ & $\\SI{" << TMethod_fit->GetParError(parNum) << "}{}$ \\\\" << endl;
          }
      }
    }

    if(inputFile->Get(Form("topDir/FitPasses/FitPass%i/addedDir/FullRatio/Added_Times_Full_Ratio_Graph",fitNum)) != 0)
    {
      TF1* fullRatio_fit = (TF1*) ((TGraphErrors*) inputFile->Get(Form("topDir/FitPasses/FitPass%i/addedDir/FullRatio/Added_Times_Full_Ratio_Graph",fitNum)))->GetFunction("fullRatioFitFunc");
      cout << endl << "Full Ratio fit #chi^{2}: " << fullRatio_fit->GetChisquare() << " NDF: " << fullRatio_fit->GetNDF() << " #chi^{2}/NDF: " << fullRatio_fit->GetChisquare()/fullRatio_fit->GetNDF() << " p value: " << fullRatio_fit->GetProb() << endl;

      if(justPrintR) cout << endl << "Full Ratio fit R: " << fullRatio_fit->GetParameter(1) << " error: " << fullRatio_fit->GetParError(1) << endl;
      else{
        cout << endl << "Full Ratio fit parameter results: " << endl << endl;
        for (int parNum = 0; parNum < fullRatio_fit->GetNpar(); ++parNum) cout << fullRatio_fit->GetParName(parNum) << " value: " << fullRatio_fit->GetParameter(parNum) << " error: " << fullRatio_fit->GetParError(parNum) << endl;

          if(printInLatexFormat){
            cout << endl << "And in LaTex format: " << endl << endl;
            for (int parNum = 0; parNum < fullRatio_fit->GetNpar(); ++parNum) cout << "$\\SI{" << fullRatio_fit->GetParameter(parNum) << "}{}$ & $\\SI{" << fullRatio_fit->GetParError(parNum) << "}{}$ \\\\" << endl;        
          }
      }
    }
}

else{

      TNtuple *firstTuple = (TNtuple*)inputFile->Get(Form("topDir/FitPasses/FitPass0/FitConditions0"));
      float tP;
      firstTuple->SetBranchAddress("totalPasses", &tP);
      firstTuple->GetEntry(0);

    for (int fitPass = 0; fitPass < tP; ++fitPass)
    {
      cout << endl << "Fit num: " << fitPass << endl;

      if(inputFile->Get(Form("topDir/FitPasses/FitPass%i/addedDir/TMethod/allTimesAdded_TMethod",fitPass)) != 0)
      {
        TF1* TMethod_fit = (TF1*) ((TH1F*) inputFile->Get(Form("topDir/FitPasses/FitPass%i/addedDir/TMethod/allTimesAdded_TMethod",fitPass)))->GetFunction("TmethodFitFunc");
        cout << endl << "T Method fit #chi^{2}: " << TMethod_fit->GetChisquare() << " NDF: " << TMethod_fit->GetNDF() << " #chi^{2}/NDF: " << TMethod_fit->GetChisquare()/TMethod_fit->GetNDF() << " p value: " << TMethod_fit->GetProb() << endl;

        if(justPrintR) cout << endl << "T Method fit R: " << TMethod_fit->GetParameter(3) << " error: " << TMethod_fit->GetParError(3) << endl;
        else{
          cout << endl << "T Method parameter results: " << endl << endl;
          for (int parNum = 0; parNum < TMethod_fit->GetNpar(); ++parNum) cout << TMethod_fit->GetParName(parNum) << " value: " << TMethod_fit->GetParameter(parNum) << " error: " << TMethod_fit->GetParError(parNum) << endl; 
        }
      }

      if(inputFile->Get(Form("topDir/FitPasses/FitPass%i/addedDir/FullRatio/Added_Times_Full_Ratio_Graph",fitPass)) != 0)
      {
        TF1* fullRatio_fit = (TF1*) ((TGraphErrors*) inputFile->Get(Form("topDir/FitPasses/FitPass%i/addedDir/FullRatio/Added_Times_Full_Ratio_Graph",fitPass)))->GetFunction("fullRatioFitFunc");
        cout << endl << "Full Ratio fit #chi^{2}: " << fullRatio_fit->GetChisquare() << " NDF: " << fullRatio_fit->GetNDF() << " #chi^{2}/NDF: " << fullRatio_fit->GetChisquare()/fullRatio_fit->GetNDF() << " p value: " << fullRatio_fit->GetProb() << endl;

        if(justPrintR) cout << endl << "Full Ratio fit R: " << fullRatio_fit->GetParameter(1) << " error: " << fullRatio_fit->GetParError(1) << endl;
        else{
          cout << endl << "Full Ratio fit parameter results: " << endl << endl;
          for (int parNum = 0; parNum < fullRatio_fit->GetNpar(); ++parNum) cout << fullRatio_fit->GetParName(parNum) << " value: " << fullRatio_fit->GetParameter(parNum) << " error: " << fullRatio_fit->GetParError(parNum) << endl;
        }
      }

    }
}

/////////////////////////////////////////////////////////////////////////////////////

	cout << endl;
	return 1;

}
