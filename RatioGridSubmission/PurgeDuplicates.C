// 4-1-20: This macro purges duplicate cycles in a ROOT file, typically introduced by hadd on objects which can't add directly, like functions, vectors, etc.

#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"

void RecursivePurgeDir(TDirectory *source) {
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
      subdir->Purge();
      RecursivePurgeDir(subdir);
    } 
  }
}


void PurgeDuplicates(const char *fname) {
   TFile* f = TFile::Open(fname, "UPDATE");
   f->Purge();
   RecursivePurgeDir(f);
   f->Write();
   f->Close();
}
