//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jul  3 09:30:40 2016 by ROOT version 5.34/36
// from TTree tree/tree
// found on file: map00335.root
//////////////////////////////////////////////////////////

#ifndef spulse_h
#define spulse_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class spulse {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           det;
   Int_t           fr;
   Int_t           str;
   ULong64_t       ts;
   Double_t        er;
   Double_t        e;
   Int_t           np;
   Double_t        ep[2];   //[np]
   ULong64_t       tsp[2];   //[np]
   Double_t        eh[2];   //[np]
   Double_t        dwave[1500];
   Short_t         cr[1500];
   UShort_t        sample[1500];
   Int_t           bd;
   UChar_t         ch;
   Int_t           nevt;

   // List of branches
   TBranch        *b_det;   //!
   TBranch        *b_fr;   //!
   TBranch        *b_str;   //!
   TBranch        *b_ts;   //!
   TBranch        *b_er;   //!
   TBranch        *b_e;   //!
   TBranch        *b_np;   //!
   TBranch        *b_ep;   //!
   TBranch        *b_tsp;   //!
   TBranch        *b_eh;   //!
   TBranch        *b_dwave;   //!
   TBranch        *b_cr;   //!
   TBranch        *b_sample;   //!
   TBranch        *b_bd;   //!
   TBranch        *b_ch;   //!
   TBranch        *b_nevt;   //!

   spulse(TTree *tree=0);
   virtual ~spulse();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef spulse_cxx
spulse::spulse(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("map00335.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("map00335.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

spulse::~spulse()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t spulse::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t spulse::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void spulse::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("det", &det, &b_det);
   fChain->SetBranchAddress("fr", &fr, &b_fr);
   fChain->SetBranchAddress("str", &str, &b_str);
   fChain->SetBranchAddress("ts", &ts, &b_ts);
   fChain->SetBranchAddress("er", &er, &b_er);
   fChain->SetBranchAddress("e", &e, &b_e);
   fChain->SetBranchAddress("np", &np, &b_np);
   fChain->SetBranchAddress("ep", ep, &b_ep);
   fChain->SetBranchAddress("tsp", tsp, &b_tsp);
   fChain->SetBranchAddress("eh", eh, &b_eh);
   fChain->SetBranchAddress("dwave", dwave, &b_dwave);
   fChain->SetBranchAddress("cr", cr, &b_cr);
   fChain->SetBranchAddress("sample", sample, &b_sample);
   fChain->SetBranchAddress("bd", &bd, &b_bd);
   fChain->SetBranchAddress("ch", &ch, &b_ch);
   fChain->SetBranchAddress("nevt", &nevt, &b_nevt);
   Notify();
}

Bool_t spulse::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void spulse::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t spulse::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef spulse_cxx
