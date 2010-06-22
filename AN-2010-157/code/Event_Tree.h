//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr 29 11:43:24 2010 by ROOT version 5.17/02
// from TTree Event_Tree/Event Tree
// found on file: CMS_bjets_data.root
//////////////////////////////////////////////////////////

#ifndef Event_Tree_h
#define Event_Tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class Event_Tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Int_t           run;
   Int_t           event;
   Int_t           lumiblock;
   Float_t         met;
   Float_t         metx;
   Float_t         mety;
   Float_t         met_phi;
   Float_t         rawmet;
   Float_t         rawmet_phi;
   Float_t         vtxmet;
   Float_t         vtxmet_phi;
   Int_t           num_vert;
   Int_t           num_jet;
   Float_t         jet_et[64];   //[num_jet]
   Float_t         jet_px[64];   //[num_jet]
   Float_t         jet_py[64];   //[num_jet]
   Float_t         jet_pz[64];   //[num_jet]
   Float_t         jet_eta[64];   //[num_jet]
   Float_t         jet_E[64];   //[num_jet]
   Float_t         jet_theta[64];   //[num_jet]
   Float_t         jet_phi[64];   //[num_jet]
   Float_t         jet_emf[64];   //[num_jet]
   Float_t         jet_fHPD[64];   //[num_jet]
   Float_t         jet_fRBX[64];   //[num_jet]
   Int_t           jet_n90Hits[64];   //[num_jet]
   Float_t         jet_btag_jetProb[64];   //[num_jet]
   Float_t         jet_scalefac[64];   //[num_jet]
   Int_t           jet_remove[64];   //[num_jet]
   Float_t         ele1_et;
   Float_t         ele1_px;
   Float_t         ele1_py;
   Float_t         ele1_pz;
   Float_t         ele1_E;
   Float_t         ele1_eta;
   Float_t         ele1_phi;
   Int_t           ele1_pass;
   Float_t         ele2_et;
   Float_t         ele2_px;
   Float_t         ele2_py;
   Float_t         ele2_pz;
   Float_t         ele2_E;
   Float_t         ele2_eta;
   Float_t         ele2_phi;
   Int_t           ele2_pass;
   Int_t           num_ele;
   Float_t         mu1_et;
   Float_t         mu1_px;
   Float_t         mu1_py;
   Float_t         mu1_pz;
   Float_t         mu1_pt;
   Float_t         mu1_E;
   Float_t         mu1_eta;
   Float_t         mu1_phi;
   Float_t         mu1_emE;
   Float_t         mu1_hadE;
   Int_t           mu1_pass;
   Float_t         mu2_et;
   Float_t         mu2_px;
   Float_t         mu2_py;
   Float_t         mu2_pz;
   Float_t         mu2_pt;
   Float_t         mu2_E;
   Float_t         mu2_eta;
   Float_t         mu2_phi;
   Float_t         mu2_emE;
   Float_t         mu2_hadE;
   Int_t           mu2_pass;
   Int_t           num_mu;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumiblock;   //!
   TBranch        *b_met;   //!
   TBranch        *b_metx;   //!
   TBranch        *b_mety;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_rawmet;   //!
   TBranch        *b_rawmet_phi;   //!
   TBranch        *b_vtxmet;   //!
   TBranch        *b_vtxmet_phi;   //!
   TBranch        *b_num_vert;   //!
   TBranch        *b_num_jet;   //!
   TBranch        *b_jet_et;   //!
   TBranch        *b_jet_px;   //!
   TBranch        *b_jet_py;   //!
   TBranch        *b_jet_pz;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_E;   //!
   TBranch        *b_jet_theta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_emf;   //!
   TBranch        *b_jet_fHPD;   //!
   TBranch        *b_jet_fRBX;   //!
   TBranch        *b_jet_n90Hits;   //!
   TBranch        *b_jet_btag_jetProb;   //!
   TBranch        *b_jet_scalefac;   //!
   TBranch        *b_jet_remove;   //!
   TBranch        *b_ele1_et;   //!
   TBranch        *b_ele1_px;   //!
   TBranch        *b_ele1_py;   //!
   TBranch        *b_ele1_pz;   //!
   TBranch        *b_ele1_E;   //!
   TBranch        *b_ele1_eta;   //!
   TBranch        *b_ele1_phi;   //!
   TBranch        *b_ele1_pass;   //!
   TBranch        *b_ele2_et;   //!
   TBranch        *b_ele2_px;   //!
   TBranch        *b_ele2_py;   //!
   TBranch        *b_ele2_pz;   //!
   TBranch        *b_ele2_E;   //!
   TBranch        *b_ele2_eta;   //!
   TBranch        *b_ele2_phi;   //!
   TBranch        *b_ele2_pass;   //!
   TBranch        *b_num_ele;   //!
   TBranch        *b_mu1_et;   //!
   TBranch        *b_mu1_px;   //!
   TBranch        *b_mu1_py;   //!
   TBranch        *b_mu1_pz;   //!
   TBranch        *b_mu1_pt;   //!
   TBranch        *b_mu1_E;   //!
   TBranch        *b_mu1_eta;   //!
   TBranch        *b_mu1_phi;   //!
   TBranch        *b_mu1_emE;   //!
   TBranch        *b_mu1_hadE;   //!
   TBranch        *b_mu1_pass;   //!
   TBranch        *b_mu2_et;   //!
   TBranch        *b_mu2_px;   //!
   TBranch        *b_mu2_py;   //!
   TBranch        *b_mu2_pz;   //!
   TBranch        *b_mu2_pt;   //!
   TBranch        *b_mu2_E;   //!
   TBranch        *b_mu2_eta;   //!
   TBranch        *b_mu2_phi;   //!
   TBranch        *b_mu2_emE;   //!
   TBranch        *b_mu2_hadE;   //!
   TBranch        *b_mu2_pass;   //!
   TBranch        *b_num_mu;   //!

   Event_Tree(TTree *tree=0);
   virtual ~Event_Tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Event_Tree_cxx
Event_Tree::Event_Tree(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("CMS_bjets_data.root");
      if (!f) {
         f = new TFile("CMS_bjets_data.root");
      }
      tree = (TTree*)gDirectory->Get("Event_Tree");

   }
   Init(tree);
}

Event_Tree::~Event_Tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Event_Tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Event_Tree::LoadTree(Long64_t entry)
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

void Event_Tree::Init(TTree *tree)
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
   fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("metx", &metx, &b_metx);
   fChain->SetBranchAddress("mety", &mety, &b_mety);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("rawmet", &rawmet, &b_rawmet);
   fChain->SetBranchAddress("rawmet_phi", &rawmet_phi, &b_rawmet_phi);
   fChain->SetBranchAddress("vtxmet", &vtxmet, &b_vtxmet);
   fChain->SetBranchAddress("vtxmet_phi", &vtxmet_phi, &b_vtxmet_phi);
   fChain->SetBranchAddress("num_vert", &num_vert, &b_num_vert);
   fChain->SetBranchAddress("num_jet", &num_jet, &b_num_jet);
   fChain->SetBranchAddress("jet_et", jet_et, &b_jet_et);
   fChain->SetBranchAddress("jet_px", jet_px, &b_jet_px);
   fChain->SetBranchAddress("jet_py", jet_py, &b_jet_py);
   fChain->SetBranchAddress("jet_pz", jet_pz, &b_jet_pz);
   fChain->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_E", jet_E, &b_jet_E);
   fChain->SetBranchAddress("jet_theta", jet_theta, &b_jet_theta);
   fChain->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_emf", jet_emf, &b_jet_emf);
   fChain->SetBranchAddress("jet_fHPD", jet_fHPD, &b_jet_fHPD);
   fChain->SetBranchAddress("jet_fRBX", jet_fRBX, &b_jet_fRBX);
   fChain->SetBranchAddress("jet_n90Hits", jet_n90Hits, &b_jet_n90Hits);
   fChain->SetBranchAddress("jet_btag_jetProb", jet_btag_jetProb, &b_jet_btag_jetProb);
   fChain->SetBranchAddress("jet_scalefac", jet_scalefac, &b_jet_scalefac);
   fChain->SetBranchAddress("jet_remove", jet_remove, &b_jet_remove);
   fChain->SetBranchAddress("ele1_et", &ele1_et, &b_ele1_et);
   fChain->SetBranchAddress("ele1_px", &ele1_px, &b_ele1_px);
   fChain->SetBranchAddress("ele1_py", &ele1_py, &b_ele1_py);
   fChain->SetBranchAddress("ele1_pz", &ele1_pz, &b_ele1_pz);
   fChain->SetBranchAddress("ele1_E", &ele1_E, &b_ele1_E);
   fChain->SetBranchAddress("ele1_eta", &ele1_eta, &b_ele1_eta);
   fChain->SetBranchAddress("ele1_phi", &ele1_phi, &b_ele1_phi);
   fChain->SetBranchAddress("ele1_pass", &ele1_pass, &b_ele1_pass);
   fChain->SetBranchAddress("ele2_et", &ele2_et, &b_ele2_et);
   fChain->SetBranchAddress("ele2_px", &ele2_px, &b_ele2_px);
   fChain->SetBranchAddress("ele2_py", &ele2_py, &b_ele2_py);
   fChain->SetBranchAddress("ele2_pz", &ele2_pz, &b_ele2_pz);
   fChain->SetBranchAddress("ele2_E", &ele2_E, &b_ele2_E);
   fChain->SetBranchAddress("ele2_eta", &ele2_eta, &b_ele2_eta);
   fChain->SetBranchAddress("ele2_phi", &ele2_phi, &b_ele2_phi);
   fChain->SetBranchAddress("ele2_pass", &ele2_pass, &b_ele2_pass);
   fChain->SetBranchAddress("num_ele", &num_ele, &b_num_ele);
   fChain->SetBranchAddress("mu1_et", &mu1_et, &b_mu1_et);
   fChain->SetBranchAddress("mu1_px", &mu1_px, &b_mu1_px);
   fChain->SetBranchAddress("mu1_py", &mu1_py, &b_mu1_py);
   fChain->SetBranchAddress("mu1_pz", &mu1_pz, &b_mu1_pz);
   fChain->SetBranchAddress("mu1_pt", &mu1_pt, &b_mu1_pt);
   fChain->SetBranchAddress("mu1_E", &mu1_E, &b_mu1_E);
   fChain->SetBranchAddress("mu1_eta", &mu1_eta, &b_mu1_eta);
   fChain->SetBranchAddress("mu1_phi", &mu1_phi, &b_mu1_phi);
   fChain->SetBranchAddress("mu1_emE", &mu1_emE, &b_mu1_emE);
   fChain->SetBranchAddress("mu1_hadE", &mu1_hadE, &b_mu1_hadE);
   fChain->SetBranchAddress("mu1_pass", &mu1_pass, &b_mu1_pass);
   fChain->SetBranchAddress("mu2_et", &mu2_et, &b_mu2_et);
   fChain->SetBranchAddress("mu2_px", &mu2_px, &b_mu2_px);
   fChain->SetBranchAddress("mu2_py", &mu2_py, &b_mu2_py);
   fChain->SetBranchAddress("mu2_pz", &mu2_pz, &b_mu2_pz);
   fChain->SetBranchAddress("mu2_pt", &mu2_pt, &b_mu2_pt);
   fChain->SetBranchAddress("mu2_E", &mu2_E, &b_mu2_E);
   fChain->SetBranchAddress("mu2_eta", &mu2_eta, &b_mu2_eta);
   fChain->SetBranchAddress("mu2_phi", &mu2_phi, &b_mu2_phi);
   fChain->SetBranchAddress("mu2_emE", &mu2_emE, &b_mu2_emE);
   fChain->SetBranchAddress("mu2_hadE", &mu2_hadE, &b_mu2_hadE);
   fChain->SetBranchAddress("mu2_pass", &mu2_pass, &b_mu2_pass);
   fChain->SetBranchAddress("num_mu", &num_mu, &b_num_mu);
   Notify();
}

Bool_t Event_Tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Event_Tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Event_Tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Event_Tree_cxx
