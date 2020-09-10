// 3-31-20: A macro which is supposed to loop through a ROOT file and print all TGraphs and histograms as pngs. Hasn't been used in a while so I'm not sure if it would need updating or not.

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
#include <TCanvas.h>
#include <TClass.h>
#include <TFile.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TImage.h>
#include <sstream>
#include <TKey.h>

// this class recursively prints root objects as pngs into a sub-directory rootPrint (this dir has to be created beforehand)
// objects with the same name (but different directories in the root file) will be overwritten with the last one

void readdir(TDirectory *dir) {
   TDirectory *dirsav = gDirectory;
   TIter next(dir->GetListOfKeys());
   TKey *key;
   while ((key = (TKey*)next())) {
      if (key->IsFolder()) {

          if(gROOT->GetClass(key->GetClassName())->InheritsFrom("TTree")){
            cout << key->GetName() << " is an object inheriting from TTree, so skip." << endl;
            continue;
          }

         dir->cd(key->GetName());
         TDirectory *subdir = gDirectory;
         readdir(subdir);
         dirsav->cd();
         continue;
      }

      TCanvas c1;

      TClass *cl = gROOT->GetClass(key->GetClassName());

      string name;

      if (cl->InheritsFrom("TCanvas")){
         TCanvas *h = (TCanvas*)(key->ReadObj());
         h->Draw();
         name = h->GetName();
         h->Print(("rootPrint/"+name+".png").c_str());
      }

      else if (cl->InheritsFrom("TH1")){
         TH1 *h = (TH1*)key->ReadObj();
         h->Draw();
         name = h->GetName();
         c1.Print(("rootPrint/"+name+".png").c_str());
      }

      else if (cl->InheritsFrom("TGraph")){
         TGraph *h = (TGraph*)key->ReadObj();
         h->Draw("AP");
         name = h->GetName();
         c1.Print(("rootPrint/"+name+".png").c_str());
      }

   }


}

int printPlots(std::string filePath) {

  TFile *inputFile = TFile::Open(filePath.c_str());
   if (inputFile == 0) {
      printf("Error: cannot open file\n");
      return 0;
   }

   gROOT->SetBatch(kTRUE); // set batch mode to true for this macro so that nothing draws to the screen
   
   gStyle->SetOptStat(0);

   readdir(inputFile);

   return 1;

}

