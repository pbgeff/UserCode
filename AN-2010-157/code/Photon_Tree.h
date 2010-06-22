//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr 23 09:16:06 2010 by ROOT version 5.17/02
// from TTree Photon_Tree/Photon Tree
// found on file: CMS_phojet_data.root
//////////////////////////////////////////////////////////

#ifndef Photon_Tree_h
#define Photon_Tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class Photon_Tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Int_t           run;
   Int_t           event;
   Int_t           num_pho;
   Float_t         pho_et[64];   //[num_pho]
   Float_t         pho_px[64];   //[num_pho]
   Float_t         pho_py[64];   //[num_pho]
   Float_t         pho_pz[64];   //[num_pho]
   Float_t         pho_eta[64];   //[num_pho]
   Float_t         pho_E[64];   //[num_pho]
   Float_t         pho_theta[64];   //[num_pho]
   Float_t         pho_phi[64];   //[num_pho]
   Int_t           pho_hasPixelSeed[64];   //[num_pho]
   Float_t         pho_hadOverEM[64];   //[num_pho]
   Float_t         pho_isoEcalRecHitDR04[64];   //[num_pho]
   Float_t         pho_isoHcalRecHitDR04[64];   //[num_pho]
   Float_t         pho_isoHollowTrkConeDR04[64];   //[num_pho]
   Float_t         pho_maxEnergyXtal[64];   //[num_pho]
   Float_t         pho_e3x3[64];   //[num_pho]
   Float_t         pho_sigmaIetaIeta[64];   //[num_pho]
   Float_t         pho_scPhiWidth[64];   //[num_pho]
   Float_t         pho_scEtaWidth[64];   //[num_pho]
   Float_t         pho_r9[64];   //[num_pho]
   Int_t           pho_istight[64];   //[num_pho]
   Int_t           pho_isloose[64];   //[num_pho]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_num_pho;   //!
   TBranch        *b_pho_et;   //!
   TBranch        *b_pho_px;   //!
   TBranch        *b_pho_py;   //!
   TBranch        *b_pho_pz;   //!
   TBranch        *b_pho_eta;   //!
   TBranch        *b_pho_E;   //!
   TBranch        *b_pho_theta;   //!
   TBranch        *b_pho_phi;   //!
   TBranch        *b_pho_hasPixelSeed;   //!
   TBranch        *b_pho_hadOverEM;   //!
   TBranch        *b_pho_isoEcalRecHitDR04;   //!
   TBranch        *b_pho_isoHcalRecHitDR04;   //!
   TBranch        *b_pho_isoHollowTrkConeDR04;   //!
   TBranch        *b_pho_maxEnergyXtal;   //!
   TBranch        *b_pho_e3x3;   //!
   TBranch        *b_pho_sigmaIetaIeta;   //!
   TBranch        *b_pho_scPhiWidth;   //!
   TBranch        *b_pho_scEtaWidth;   //!
   TBranch        *b_pho_r9;   //!
   TBranch        *b_pho_istight;   //!
   TBranch        *b_pho_isloose;   //!

   Photon_Tree(TTree *tree=0);
   virtual ~Photon_Tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Photon_Tree_cxx
Photon_Tree::Photon_Tree(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("CMS_phojet_data.root");
      if (!f) {
         f = new TFile("CMS_phojet_data.root");
      }
      tree = (TTree*)gDirectory->Get("Photon_Tree");

   }
   Init(tree);
}

Photon_Tree::~Photon_Tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Photon_Tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Photon_Tree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Photon_Tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("num_pho", &num_pho, &b_num_pho);
   fChain->SetBranchAddress("pho_et", pho_et, &b_pho_et);
   fChain->SetBranchAddress("pho_px", pho_px, &b_pho_px);
   fChain->SetBranchAddress("pho_py", pho_py, &b_pho_py);
   fChain->SetBranchAddress("pho_pz", pho_pz, &b_pho_pz);
   fChain->SetBranchAddress("pho_eta", pho_eta, &b_pho_eta);
   fChain->SetBranchAddress("pho_E", pho_E, &b_pho_E);
   fChain->SetBranchAddress("pho_theta", pho_theta, &b_pho_theta);
   fChain->SetBranchAddress("pho_phi", pho_phi, &b_pho_phi);
   fChain->SetBranchAddress("pho_hasPixelSeed", pho_hasPixelSeed, &b_pho_hasPixelSeed);
   fChain->SetBranchAddress("pho_hadOverEM", pho_hadOverEM, &b_pho_hadOverEM);
   fChain->SetBranchAddress("pho_isoEcalRecHitDR04", pho_isoEcalRecHitDR04, &b_pho_isoEcalRecHitDR04);
   fChain->SetBranchAddress("pho_isoHcalRecHitDR04", pho_isoHcalRecHitDR04, &b_pho_isoHcalRecHitDR04);
   fChain->SetBranchAddress("pho_isoHollowTrkConeDR04", pho_isoHollowTrkConeDR04, &b_pho_isoHollowTrkConeDR04);
   fChain->SetBranchAddress("pho_maxEnergyXtal", pho_maxEnergyXtal, &b_pho_maxEnergyXtal);
   fChain->SetBranchAddress("pho_e3x3", pho_e3x3, &b_pho_e3x3);
   fChain->SetBranchAddress("pho_sigmaIetaIeta", pho_sigmaIetaIeta, &b_pho_sigmaIetaIeta);
   fChain->SetBranchAddress("pho_scPhiWidth", pho_scPhiWidth, &b_pho_scPhiWidth);
   fChain->SetBranchAddress("pho_scEtaWidth", pho_scEtaWidth, &b_pho_scEtaWidth);
   fChain->SetBranchAddress("pho_r9", pho_r9, &b_pho_r9);
   fChain->SetBranchAddress("pho_istight", pho_istight, &b_pho_istight);
   fChain->SetBranchAddress("pho_isloose", pho_isloose, &b_pho_isloose);
   Notify();
}

Bool_t Photon_Tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Photon_Tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Photon_Tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Photon_Tree_cxx
