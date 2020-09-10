// Macro to subtract histograms in two root files in order to compare the differences between them. The base code was taken from the hadd.C script before modification.

#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

TList *FileList;
TFile *Target;
bool nonZeroFile = false;

void RecursiveSubtract(TDirectory *target, TList *sourcelist);

int SubtractRootFiles(string filePath1, string filePath2) {

   // Prepare the files to me merged
  TFile *inputFile1 = TFile::Open(filePath1.c_str());
  TFile *inputFile2 = TFile::Open(filePath2.c_str());

   if (inputFile1 == 0 || inputFile2 == 0) {
      printf("Error: cannot open file\n");
      return 1;
   }

   Target = TFile::Open( "differenceHists.root", "RECREATE" );
   FileList = new TList();
   FileList->Add( inputFile1 );
   FileList->Add( inputFile2  );
   RecursiveSubtract( Target, FileList );

   if(nonZeroFile) std::cout << "Difference histograms found to be non-zero." << std::endl;
   else std::cout << "All difference histograms were found to be zero." << std::endl;

   return 0;
}
void RecursiveSubtract(TDirectory *target, TList *sourcelist) {
   //  cout << "Target path: " << target->GetPath() << endl;
   TString path( (char*)strstr( target->GetPath(), ":" ) );
   path.Remove( 0, 2 );
   TFile *first_source = (TFile*)sourcelist->First();
   first_source->cd( path );
   TDirectory *current_sourcedir = gDirectory;
   //gain time, do not add the objects in the list in memory
   Bool_t status = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);
   // loop over all keys in this directory
   TChain *globChain = 0;
   TIter nextkey( current_sourcedir->GetListOfKeys() );
   TKey *key, *oldkey=0;
   while ( (key = (TKey*)nextkey())) {
      //keep only the highest cycle number for each key
      if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;
      // read object from first source file
      first_source->cd( path );
      TObject *obj = key->ReadObj();
      if (obj->IsA()->InheritsFrom(TH1::Class())) {
         // descendant of TH1 -> merge it
         //      cout << "Merging histogram " << obj->GetName() << endl;
         TH1 *h1 = (TH1*)obj;
         // loop over all source files and add the content of the
         // correspondant histogram to the one pointed to by "h1"
         TFile *nextsource = (TFile*)sourcelist->After( first_source );
         while ( nextsource ) {
            // make sure we are at the correct directory level by cd'ing to path
            nextsource->cd( path );
            TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(h1->GetName());
            if (key2) {
               TH1 *h2 = (TH1*)key2->ReadObj();
               h1->Add(h2, -1);

               // set bin errors to zero, since adding them doesn't do that
               std:vector<double> zeroErrors(h1->GetNbinsX(), 0);
               h1->SetError(&zeroErrors[0]);

               // cout to the terminal if a subtracted histogram is non-zero, so that the user doesn't necessarily need to open the root file and check the histograms manually
               for (int bin = 0; bin < h1->GetNbinsX(); ++bin){
                  if(h1->GetBinContent(bin) != 0){
                    std::cout << "Histogram " << h1->GetName() << " is non-zero." << std::endl;
                    nonZeroFile = true;
                    break; 
                  }
               }

               delete h2;
            }
            nextsource = (TFile*)sourcelist->After( nextsource );
         }
      }
      else if ( obj->IsA()->InheritsFrom( TDirectory::Class() ) ) {
         // it's a subdirectory
         cout << "Entering subdirectory: " << obj->GetName() << endl;
         // create a new subdir of same name and title in the target file
         target->cd();
         TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );
         // newdir is now the starting point of another round of merging
         // newdir still knows its depth within the target file via
         // GetPath(), so we can still figure out where we are in the recursion
         RecursiveSubtract( newdir, sourcelist );
      } 
      // now write the merged histogram (which is "in" obj) to the target file
      // note that this will just store obj in the current directory level,
      // which is not persistent until the complete directory itself is stored
      // by "target->Write()" below
      if ( obj ) {
         target->cd();
         if(obj->IsA()->InheritsFrom(TH1::Class())) obj->Write( key->GetName() ); // only write histograms to the file, drop everything else
      }
   } // while ( ( TKey *key = (TKey*)nextkey() ) )
   // save modifications to target file
   target->SaveSelf(kTRUE);
   TH1::AddDirectory(status);
}
