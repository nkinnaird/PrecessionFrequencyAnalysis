// This file a modified variant of : $ROOTSYS/tutorials/copyFiles.C
// 3-31-20: It copies an input ROOT file into a subdirectory of an output ROOT file.

#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"


void CopyDir(TDirectory *source) {
   //copy all objects and subdirs of directory source as a subdir of the current directory
   // source->ls();
   TDirectory *savdir = gDirectory;
   string sourcename = string(source->GetName());

   // use this block of code to avoid creating ../../ and other extraneous directories when running on a file located in a different folder
   TDirectory *adir = new TDirectory();
   size_t slashChar = sourcename.find_last_of("/");
   if(slashChar >= sourcename.size()) adir = savdir->mkdir(sourcename.c_str());
   else adir = savdir->mkdir(sourcename.substr(slashChar+1).c_str());

   // TDirectory *adir = savdir->mkdir(sourcename.c_str());
   adir->cd();

   //loop on all entries of this directory
   TKey *key;
   TIter nextkey(source->GetListOfKeys());
   while ((key = (TKey*)nextkey())) {
      const char *classname = key->GetClassName();
      const char *keyname = key->GetName();
      TClass *cl = gROOT->GetClass(classname);
      if (!cl) continue;
      if (cl->InheritsFrom(TDirectory::Class())) {
         source->cd(keyname);
         TDirectory *subdir = gDirectory;
         adir->cd();
         CopyDir(subdir);
         adir->cd();
      } else if (cl->InheritsFrom(TTree::Class())) {
         TTree *T = (TTree*)source->Get(keyname);
         adir->cd();
         TTree *newT = T->CloneTree(-1,"fast");
         newT->Write();
      } else {
         source->cd();
         TObject *obj = key->ReadObj();
         adir->cd();
         if(!adir->GetListOfKeys()->Contains(keyname)) obj->Write(keyname); // need the GetName here otherwise the TVectorDs don't save with the proper names, also don't copy duplicate objects
         delete obj;
     }
  }
  adir->SaveSelf(kTRUE);
  savdir->cd();
}

void CopyFile(const char *fname) {
   //Copy all objects and subdirs of file fname as a subdir of the current directory
   TDirectory *target = gDirectory;
   TFile* f = TFile::Open(fname);
   target->cd(); // opening a file changes gDirectory to that file, so we need to switch back to the output file directory

   // TNamed n("inputFilePath", fname); // store relative path of file being ran on
   // n.Write();

   CopyDir(f);
   f->Close();
   target->cd();
}
