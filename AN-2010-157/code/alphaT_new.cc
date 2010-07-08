#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TH3F.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include "TString.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "Event_Tree.C"
#include "Photon_Tree.C"
using namespace std;


float smear_jet_et(float et, int type);
float smear_phi(float phi, int type);

double deltaphi(double phi1, double phi2)
{
  double result = phi1-phi2;
  while (result>TMath::Pi()) result -= 2*TMath::Pi();
  while (result<=-TMath::Pi()) result += 2*TMath::Pi();
  return result;
}

double findphi(double x, double y)
{
  //to compute missing Ht phi
  double pi = 3.141592653589;
  double phi, phi1, phi2;
  double c;
  c = sqrt(x*x+y*y);
  phi2= asin(y/c);
  phi1= acos(x/c);
  if(phi2>0 && phi1<pi/2.) phi=phi1;
  if(phi2>0 && phi1>pi/2.) phi=phi1;
  if(phi2<0 && phi1<pi/2.) phi=phi2;
  if(phi2<0 && phi1>pi/2.) phi=-phi1;
  return phi;
}

void event_dump(Event_Tree *evt2, int num_removed, float* ht_vec, float ht, TLorentzVector DiMuon);

const float pi = 3.141592653589;

void plot_maker(int sample, int n_jet=-1, int mode=-1, int test=0, int soft=0, int which_file=1, int jet_smear_flag=0, int z_smear=0)
{

  gRandom->SetSeed(0);

  int debug=0; // Debugging level. 0=none, 100=trace, 300=detailed trace

  //mode==-1 don't cut on Z mass window (for photon study)
  //mode==0 Z mass cut [81,101]; no sideband subtraction
  //mode==1 [71,81]&[101,111] sideband subtraction
  //mode==2  [61,71]&[111,121] sideband subtraction

  //test==0 usual sample
  //test==1 N-jet sample where any jet can be removed
  //test==2 N-jet sample where any jet can be removed, alternate values

  //soft==0 usual sample
  //soft==1 add in soft jets
  //soft==2 add in soft jets, alternate values  
  //soft==? add in soft jets; account for Z mass; Not implemented right now
  //soft==3 add in soft jets and then tighten Z and jet cuts

  //z_smear==0 usual
  //z_smear==1 smear Z pt 
  //z_smear==2 smear Z phi
  //z_smear==3 smear Z pt and phi

  //jet_smear_flag==0 usual
  //jet_smear_flag==1 smear pt of jets 
  //jet_smear_flag==2 smear pt of jets, alternate values

  //n_jet==-1 >=2-jet with no Z
  //n_jet==-2 2-jet with no Z
  //n_jet==-3 3-jet with no Z
  //n_jet==-4 4-jet with no Z
  //n_jet==0 Z+Njet where N>=1
  //n_jet==1 Z+1-jet 
  //n_jet==2 Z+2-jet 
  //n_jet==3 Z+3-jet 
  //n_jet==4 Z+4-jet 

  int pho_sample = 0; //Select photons instead of Z's if set to 1
                      //Photons in samples 6,

  int jet_quality = 0; //Impose jet quality cuts if set to 1
                       //Important for early CMS data

  int bjet_sample = 0; //Require >=1 bjet if set to 1


  //input flat ntuples and output histograms
  TFile* input,*hist_file; 
  //if(sample==5) {input = new TFile("CMS_MC.root"); hist_file = new TFile("hists_cms_mc.root","RECREATE");} //big Z sample
  if(sample==5) {input = new TFile("CMS_MC_oldZmumu_800pb.root"); hist_file = new TFile("hists_cms_mc.root","RECREATE");} //800 pb-1 sample
  if(sample==6)  {input = new TFile("CMS_phojet_MC.root"); hist_file = new TFile("hists_cms_mc_pho.root","RECREATE"); pho_sample=1; jet_quality=1;}  //gamma+jets MinBias MC sample
  if(sample==7)  {input = new TFile(/*"CMS_phojet_data.root"*//*"CMS_phojet_data_Commissioning10-SD_EG-v9.root"*//*"CMS_phojet_data_EG_Jun9thReReco_v1.root"*//*"CMS_phojet_data_EG_12nb.root"*//*"CMS_phojet_data_MinBias_TuneD6T_7TeV-pythia6_Spring10-START3X_V26B-v2_GEN-SIM-RECO.root"*//*"CMS_phojet_15_data_EG_12nb.root"*//*"CMS_phojet_15_data_EG_Run2010A-PromptReco-v4.root"*/"CMS_phojet_15_data_EG_26nb.root"); hist_file = new TFile("hists_cms_data_pho.root","RECREATE"); pho_sample=1; jet_quality=1;}  //gamma+jets MinBias data sample
  if(sample==-1)  {input = new TFile("CMS_jets_MC.root"/*"CMS_QCD_MC_Pthat30toInf.root"*/); hist_file = new TFile("hists_cms_mc_jets.root","RECREATE"); jet_quality=1;}  //dijet MinBias MC sample
  if(sample==-2)  {input = new TFile("CMS_jets_data.root"); hist_file = new TFile("hists_cms_data_jets.root","RECREATE"); jet_quality=1;}  //dijet MinBias data sample
  if(sample==-3)  {input = new TFile(/*"CMS_jets_data_new.root"*//*"CMS_jets_data_Commissioning10-SD_JetMETTau-v9.root"*//*"CMS_jets_data_JMT_Jun9thReReco_v1.root"*/"CMS_jets_data_JetMETTau_12nb.root"/*"CMS_JPTjets_data_JetMETTau_12nb.root"*//*"CMS_PFjets_data_JetMETTau_12nb.root"*//*"CMS_PFjets_pt10_data_JMT_Jun9thReReco_v1.root"*/); hist_file = new TFile("hists_cms_data_jets.root","RECREATE"); jet_quality=1;}  //dijet data sample
  if(sample==-4)  {input = new TFile("CMS_bjets_data.root"); hist_file = new TFile("hists_cms_data_jets.root","RECREATE"); jet_quality=1; bjet_sample=1;}  //dijet MinBias data sample; >=1 bjet 
  if(sample==-5)  {
    //               input = new TFile("CMS_QCD_MC_v9.root"); 
    //               input = new TFile("CMS_QCD_MC_Pt100to250.root"); 
    //               input = new TFile("CMS_QCD_MC_Pt250to500.root"); 
    //               input = new TFile("CMS_QCD_MC_Pt500to1000.root");
    //               input = new TFile("CMS_QCD_MC_Pt1000toInf.root");
                   input = new TFile("CMS_QCD_MC_Pthat30toInf.root");
    //               input = new TFile("CMS_QCD_MC_Pthat15toInf.root");
	  	     hist_file = new TFile("hists_cms_MC_jets.root","RECREATE"); jet_quality=1;}  //QCD MC sample

  TFile* graph1, *graph2;
  if(which_file==1) graph1 = new TFile("graph_1.root","RECREATE");  
  if(which_file==2) graph2 = new TFile("graph_2.root","RECREATE");  

  TH1F *h_alphaT = new TH1F("h_alphaT","#alpha_{T}",60,0.,3.);
  TH1F *h_alphaT_ht_20_40 = new TH1F("h_alphaT_ht_20_40","#alpha_{T}, 20<H_{T}<40",60,0.,3.);
  TH1F *h_alphaT_ht_40_up = new TH1F("h_alphaT_ht_40_up","#alpha_{T}, 40<H_{T}",60,0.,3.);
  TH1F *h_alphaT_ht_40_80 = new TH1F("h_alphaT_ht_40_80","#alpha_{T}, 40<H_{T}<80",60,0.,3.);
  TH1F *h_alphaT_ht_100_up = new TH1F("h_alphaT_ht_100_up","#alpha_{T}, 100<H_{T}",40,0.,2.);
  TH1F *h_alphaT_ht_80_up = new TH1F("h_alphaT_ht_80_up","#alpha_{T}, 80<H_{T}",60,0.,3.);
  TH1F *h_alphaT_ht_60_up = new TH1F("h_alphaT_ht_60_up","#alpha_{T}, 60<H_{T}",40,0.,2.);
  TH1F *h_alphaT_ht_30_50 = new TH1F("h_alphaT_ht_30_50","#alpha_{T}, 30<H_{T}<50",60,0.,3.);
  TH1F *h_alphaT_ht_50_80 = new TH1F("h_alphaT_ht_50_80","#alpha_{T}, 50<H_{T}<80",60,0.,3.);
  TH1F *h_alphaT_ht_80_120 = new TH1F("h_alphaT_ht_80_120","#alpha_{T}, 80<H_{T}<120",60,0.,3.);
  TH1F *h_alphaT_ht_120_up = new TH1F("h_alphaT_ht_120_up","#alpha_{T}, 120<H_{T}",60,0.,3.);
  TH1F *h_alphaT_ht_120_160 = new TH1F("h_alphaT_ht_120_160","#alpha_{T}, 120<H_{T}<160",60,0.,3.);
  TH1F *h_alphaT_ht_160_up = new TH1F("h_alphaT_ht_160_up","#alpha_{T}, 160<H_{T}",60,0.,3.);
  TH1F *h_alphaT_ht_60_100 = new TH1F("h_alphaT_ht_60_100","#alpha_{T}, 60<H_{T}<100",40,0.,2.);
  TH1F *h_alphaT_ht_100_150 = new TH1F("h_alphaT_ht_100_150","#alpha_{T}, 100<H_{T}<150",40,0.,2.);
  TH1F *h_alphaT_ht_150_200 = new TH1F("h_alphaT_ht_150_200","#alpha_{T} = #frac{E_{T,ZorJ1}}{M_{T,j1,Z}}, ==1 jet, 150<H_{T}<200",20,0.,1.);
  TH1F *h_alphaT_ht_200_250 = new TH1F("h_alphaT_ht_200_250","#alpha_{T} = #frac{E_{T,ZorJ1}}{M_{T,j1,Z}}, ==1 jet, 200<H_{T}<250",20,0.,1.);
  TH1F *h_alphaT_ht_250_300 = new TH1F("h_alphaT_ht_250_300","#alpha_{T} = #frac{E_{T,ZorJ1}}{M_{T,j1,Z}}, ==1 jet, 250<H_{T}<300",20,0.,1.);
  TH1F *h_alphaT_ht_300_500 = new TH1F("h_alphaT_ht_300_500","#alpha_{T} = #frac{E_{T,ZorJ1}}{M_{T,j1,Z}}, ==1 jet, 300<H_{T}<500",20,0.,1.);
  TH1F *h_alphaT_ht_500up = new TH1F("h_alphaT_ht_500up","#alpha_{T} = #frac{E_{T,ZorJ1}}{M_{T,j1,Z}}, ==1 jet, H_{T}>500",20,0.,1.);

  TH1F *h_jet_pt = new TH1F("h_jet_pt","Jet pt",40,0,400.);
  TH1F *h_jet_pt_wide = new TH1F("h_jet_pt_wide","Jet pt",500,0,4000.);
  TH1F *h_jet_phi = new TH1F("h_jet_phi","Jet phi",20,0,2*pi);
  TH1F *h_jet_eta = new TH1F("h_jet_eta","Jet Eta",16,-4,4.);
  TH1F *h_jet_num = new TH1F("h_jet_num","Jet Multiplicity",25,0,25);
  TH1F *h_jet_fRBX = new TH1F("h_jet_fRBX","Jet fRBX",40,0,1.0);
  TH1F *h_jet_fHPD = new TH1F("h_jet_fHPD","Jet fHPD",40,0,1.0);
  TH1F *h_jet_n90Hits = new TH1F("h_jet_n90Hits","Jet n90Hits",100,0,100);
  TH1F *h_jet_emf = new TH1F("h_jet_emf","Jet EMF",100,-1.0,1.0);
  TH1F *h_jet_btag_jetProb = new TH1F("h_jet_btag_jetProb","Jet btag_jetProb",100,0,1.0);
  TH1F *h_jet_pt_btag = new TH1F("h_jet_pt_btag","Jet pt, b-Tagged jets",50,0,400.);

  TH1F *h_pho_pt = new TH1F("h_pho_pt","Photon Pt",100,0,25./*,40,0,200*/);
  TH1F *h_pho_phi = new TH1F("h_pho_phi","Photon phi",20,-pi,pi);
  TH1F *h_pho_eta = new TH1F("h_pho_eta","Photon Eta",16,-4.,4.);
  TH1F *h_pho_E1overE3x3 = new TH1F("h_pho_E1overE3x3","Photon maxEnergyXtal/e3x3",50,0,1.);
  TH1F *h_pho_sigmaIetaIeta = new TH1F("h_pho_sigmaIetaIeta","Photon sigmaIetaIeta",50,0,0.1);
  TH1F *h_pho_hasPixelSeed = new TH1F("h_pho_hasPixelSeed","Photon hasPixelSeed",2,0,2.); 
  TH1F *h_pho_hadOverEM = new TH1F("h_pho_hadOverEM","Photon hadOverEM",50,0,0.2);
  TH1F *h_pho_isoEcalRecHitDR04 = new TH1F("h_pho_isoEcalRecHitDR04","Photon isoEcalRecHitDR04",50,0,10.);
  TH1F *h_pho_isoHcalRecHitDR04 = new TH1F("h_pho_isoHcalRecHitDR04","Photon isoHcalRecHitDR04",50,0,10.);
  TH1F *h_pho_isoHollowTrkConeDR04 = new TH1F("h_pho_isoHollowTrkConeDR04","Photon isoHollowTrkConeDR04",50,0,10.);

  TH2F *h_phopt_vs_jetpt_overlap = new TH2F("h_phopt_vs_jetpt_overlap","Jet Pt vs Photon Pt, overlap",1000,0,100.,1000,0,100);
  TH1F *h_pho_isoEcalRecHitDR04_overlap = new TH1F("h_pho_isoEcalRecHitDR04_overlap","Photon isoEcalRecHitDR04, overlapping jet",50,0,10.);
  TH1F *h_jet_emf_overlap = new TH1F("h_jet_emf_overlap","Jet EMF, overlapping photon",60,-1.0,1.1);

  TH1F *h_ht = new TH1F("h_ht","HT",100,0,2000.);
  //TH1F *h_mht = new TH1F("h_mht","Missing HT",100,0,1000.);
  TH1F *h_mht = new TH1F("h_mht","Missing HT",40,0,200.);  
  TH1F *h_met = new TH1F("h_met","Missing ET",40,0,200.);
  TH1F *h_genht = new TH1F("h_genht","genHT",100,0,2000.);
  TH2F *h_ht_vs_htplusmet = new TH2F("h_ht_vs_htplusmet","HT vs HT+MET",1000,0,1000,1000,0,1000);
  TProfile* p_ht_vs_htplusmet = new TProfile("p_ht_vs_htplusmet","HT vs HT+MET",100,0,1000);
  TH2F *h_ht_vs_genht = new TH2F("h_ht_vs_genht","HT vs genHT",1000,0,1000,1000,0,1000);
  TProfile* p_ht_vs_genht = new TProfile("p_ht_vs_genht","HT vs genHT",100,0,1000);
  TH1F *h_Npv = new TH1F("h_Npv","Number of Primary Verticies",20,0,20.);
  TH1F *h_Npv_good = new TH1F("h_Npv_good","Number of good Primary Verticies",20,0,20.);

  TH1F *h_phojet_balance = new TH1F("h_phojet_balance","Photon Pt/Jet Pt",30,0,3.);
  TH1F *h_dijet_balance = new TH1F("h_dijet_balance","Jet Pt/Jet Pt",30,0,3.);
  TH1F *h_dphi_phomet =  new TH1F("h_dphi_phomet","dPhi(Photon,MET)",30,0,pi);
  TH1F *h_dphi_phomet_met0to20 =  new TH1F("h_dphi_phomet_met0to20","dPhi(Photon,MET), MET<20",30,0,pi);
  TH1F *h_dphi_phomet_met20toinf =  new TH1F("h_dphi_phomet_met20toinf","dPhi(Photon,MET), MET>20",30,0,pi);


  TH1F *h_measuredht = new TH1F("h_measuredht","measured HT",60,0,600.);
  TH1F *h_measuredht_failalphat = new TH1F("h_measuredht_failalphat","measured HT, #alpha_{T}>0.55",60,0,600.);


  Event_Tree evt = Event_Tree((TTree*)input->Get("Event_Tree")); 
  Photon_Tree pho_tr = Photon_Tree((TTree*)input->Get("Photon_Tree")); 
  Event_Tree *evt2;

  //Number of events to loop over
  Int_t nentries = (Int_t)evt.fChain->GetEntriesFast();
  cout<<"\nNumber of entries: "<<nentries<<"\n";
 
  float mumass = 0.105658;
  float elemass = 0.000511;
  int ijet = 0; 
  int num_ele_evts=0;
  int num_mu_evts=0;
  int dump_counter=0;
  int num_evts=0;
  float weight1, weight2; //weight events based on Z mass.  used for sideband subtraction
  //weight1 is for 1D hists and weight2 is for 2D hists

  
  float leadjet_pt_min=30; //minimum pt of lead jet
  float jet_pt_min=30; //minimum pt of any jet (to ensure only 1 jet) 
  float z_pt_min=30; //minimum Z pt

  //z_pt_min=10; jet_pt_min=10; leadjet_pt_min=10; //lower cuts to 10 GeV as a check
  //z_pt_min=15; jet_pt_min=15; leadjet_pt_min=15; //lower cuts to 15 GeV as a check
  //z_pt_min=20; jet_pt_min=20; leadjet_pt_min=20; //lower cuts to 20 GeV as a check
  //z_pt_min=40; jet_pt_min=40; leadjet_pt_min=40; //cuts to 40 GeV
  //z_pt_min=50; jet_pt_min=50; leadjet_pt_min=50; //cuts to 50 GeV 
  //z_pt_min=75; jet_pt_min=75; leadjet_pt_min=75; //cuts to 75 GeV 
  //z_pt_min=100; jet_pt_min=100; leadjet_pt_min=100; //cuts to 100 GeV 
  if(sample==-3 || sample==-5) {z_pt_min=40; jet_pt_min=40; leadjet_pt_min=40;} //cuts to 40 GeV cout<<"\n"<<z_pt_min<<" "<<jet_pt_min;
  if(sample==-2) {z_pt_min=20; jet_pt_min=20; leadjet_pt_min=20;} //cuts to 20 GeV
  if(sample==7 && mode==-2) {z_pt_min=10; jet_pt_min=10; leadjet_pt_min=10; }
  if(sample==7 && mode==-1) {z_pt_min=18; jet_pt_min=18; leadjet_pt_min=18; }
  //z_pt_min=18; jet_pt_min=18; leadjet_pt_min=18;
  //z_pt_min=50; jet_pt_min=50; leadjet_pt_min=50;

  float btag_cut = 0.5;
  int min_num_jets;
  int max_num_jets = 100;
  if(n_jet==-1) {min_num_jets=1;  max_num_jets = 100;} //N-Jet case (no Z or photon)
  if(n_jet==-2) {min_num_jets=1;  max_num_jets = 1;} // 2-jet case (no Z or photon)
  if(n_jet==-3) {min_num_jets=2;  max_num_jets = 2;} // 3-jet case (no Z or photon)
  if(n_jet==-4) {min_num_jets=3;  max_num_jets = 3;} // 4-jet case (no Z or photon)
  if(n_jet==-5) {min_num_jets=2;  max_num_jets = 100;} // >=3 Jet case (no Z or photon)
  if(n_jet==0)  { min_num_jets=1;  max_num_jets = 100; } // >=1 jet
  if(n_jet==1)  { min_num_jets=1;  max_num_jets=1; } //exactly 1 jet
  if(n_jet==2)  { min_num_jets=2;  max_num_jets = 2; } //exactly 2 jets 
  if(n_jet==3)  { min_num_jets=3;  max_num_jets = 3; } //exactly 3 jets 
  if(n_jet==4)  { min_num_jets=4;  max_num_jets = 4; } //exactly 4 jets 
  if(n_jet==5)  { min_num_jets=2;  max_num_jets = 100; } // >=2 jet
  if(n_jet>5) {cout<<"\n\nInvalid Option \n\n"; exit(1);}


  //arrays to hold the number of events above alphaT=0.55 and total events 
  //and arrays to hold the ratio and error
  int max_htbins = 60;
  double num_above_smear[8][max_htbins], num_total_smear[8][max_htbins], num_total_smear_noscale[8][max_htbins];
  double ratio_smear[8][max_htbins]  , error_smear[8][max_htbins];
  double x_sum_smear[8][max_htbins], x2_sum_smear[8][max_htbins];
  for(int j=0; j<8; j++) {
    for(int i=0; i<max_htbins; i++) {
      num_above_smear[j][i]=num_total_smear[j][i]=num_total_smear_noscale[j][i]=0;
      ratio_smear[j][i]=error_smear[j][i]=0;
      x_sum_smear[j][i]=x2_sum_smear[j][i]=0;
    }
  }

  //arrays to hold the number of events above alphaT=0.55 and total events for eta bins
  //and arrays to hold the ratio and error
  float num_above_eta[8][4][3], num_total_eta[8][4][3]; // indicies are [smear][ht][eta]
  float ratio_eta[8][4][3]  , error_eta[8][4][3];
  float x_sum_eta[8][4][3], x2_sum_eta[8][4][3];
  float xerr_eta[8][4][3], x_eta[8][4][3];
  for(int j=0; j<8; j++) {
    for(int i=0; i<4; i++) {
      for(int k=0; k<3; k++){
	num_above_eta[j][i][k]=num_total_eta[j][i][k]=0;
	ratio_eta[j][i][k]=error_eta[j][i][k]=0;
	x_sum_eta[j][i][k]=x2_sum_eta[j][i][k]=0;
	x_eta[j][i][k]=xerr_eta[j][i][k]=0;
      }
    }
  }

  //Variables for failure fraction vs HT bins
  float offset;
  float bin_width; 
  float max_ht; 


  //Main Event loop
  for(int ia = 0; ia<nentries;ia++){
    if (debug>100) cout << "Entry #"<<ia<<"\n";

    //Only use every 9th event to get 100 pb-1
    //if(ia%8) continue;

    evt.GetEntry(ia);
    pho_tr.GetEntry(ia);

    //if(evt.num_vert>1) continue; //only take events with certain number of verticies
    //if(evt.run>135802) continue; 

    //Only look at events with genHT in specific range
    //if(evt.rawmet<120 || evt.rawmet>150) continue; 
    //if(evt.rawmet<170 || evt.rawmet>200) continue; 
    //if(evt.rawmet<270 || evt.rawmet>300) continue; 
    //if(evt.rawmet<520 || evt.rawmet>550) continue;

    TLorentzVector DiMuon;
    DiMuon.SetXYZM(0,0,0,0);
    int pho_index;

    if(pho_sample==0 && mode>=0){ //select a Z
      if(evt.num_mu!=2 && evt.num_ele!=2) continue; //require 2 mu or ele
      if((evt.num_mu==2 && evt.num_ele>0) || (evt.num_mu>0 && evt.num_ele==2)) continue;  //veto leptons besides the 2 that make the Z
      TLorentzVector Muon1LV;
      if(evt.num_ele==2) Muon1LV.SetXYZM(evt.ele1_px,evt.ele1_py,evt.ele1_pz,elemass);
      if(evt.num_mu==2) Muon1LV.SetXYZM(evt.mu1_px,evt.mu1_py,evt.mu1_pz,mumass);
      TLorentzVector Muon2LV;
      if(evt.num_ele==2) Muon2LV.SetXYZM(evt.ele2_px,evt.ele2_py,evt.ele2_pz,elemass);
      if(evt.num_mu==2) Muon2LV.SetXYZM(evt.mu2_px,evt.mu2_py,evt.mu2_pz,mumass);
      DiMuon = Muon1LV+Muon2LV;
      if(DiMuon.Pt()< z_pt_min) continue;
    }
    if(pho_sample==1){ //select a photon
     if(pho_tr.num_pho<1) continue;
     //How should photons be selected?  Exactly 1 above some threshold? 
     //loop through photons and make sure there is only 1 with pT above threshold
     int pho_pass=0; int many_pho=0;
     int removed_pho=0;
     for(int i=0; i<pho_tr.num_pho; i++) {
       //Apply quality cuts?
       if(pho_tr.pho_et[i]>z_pt_min && !pho_tr.pho_istight[i]) { removed_pho=1; continue;}
       if(pho_tr.pho_et[i]>z_pt_min && pho_tr.pho_istight[i]) {
	 if(fabs(pho_tr.pho_eta[i]) > 1.4/*3.0*/) { removed_pho=1; continue;}
	 if(pho_tr.pho_maxEnergyXtal[i]/pho_tr.pho_e3x3[i] > 0.9) { removed_pho=1; continue;}
	 if(pho_tr.pho_sigmaIetaIeta[i]<0.002) { removed_pho=1; continue;} 
	 
	 //if(pho_tr.pho_sigmaIetaIeta[i]>0.015) { removed_pho=1; continue;} //Only cut her for central photons
	 if(pho_pass==0) {
	   //pho_tr.pho_et[i] += (pho_tr.pho_isoEcalRecHitDR04[i]+pho_tr.pho_isoHcalRecHitDR04[i]);  //Add in isolation energy to the photon
	   pho_tr.pho_px[i] = pho_tr.pho_et[i]*cos(pho_tr.pho_phi[i]);
	   pho_tr.pho_py[i] = pho_tr.pho_et[i]*sin(pho_tr.pho_phi[i]);
	                  //pho_tr.pho_et[i] *= 2; pho_tr.pho_px[i]*= 2;  pho_tr.pho_py[i]*= 2;  pho_tr.pho_pz[i]*= 2;
           	          DiMuon.SetXYZM(pho_tr.pho_px[i],pho_tr.pho_py[i],pho_tr.pho_pz[i],0); pho_pass=1; pho_index=i; //save photon
	 }
	 else many_pho=1;
       }
     }
     if(many_pho || pho_pass==0 || removed_pho) continue; //too many or not enough good photons; Or photons not passing quality cuts
   }


    TLorentzVector DiMuon_smear[8];
    for(int i=0; i<8; i++) DiMuon_smear[i] = DiMuon;

    int smear_flag = 0;
    float jet_smear_et[evt.num_jet][8]; 
    float jet_smear_phi[evt.num_jet][8]; 
    for(ijet=0; ijet<evt.num_jet; ijet++) { //loop over jets to smear them
      if(evt.jet_remove[ijet]==1) continue; //skip removed jets 
      jet_smear_et[ijet][0]= smear_jet_et(evt.jet_et[ijet],0); //saved true jet et
      //jet_smear_et[ijet][0] *= 1.2; //change jet energy scale
      //evt.jet_et[ijet] *= 0.8; //change jet energy scale
      if(jet_smear_flag){  //Jet smearing
	float temp_pt=evt.jet_et[ijet];
	if(temp_pt<jet_pt_min) continue;
	if(jet_smear_flag==1) {
	  if(gRandom->Uniform(0,1)<0.9) evt.jet_et[ijet] = smear_jet_et(evt.jet_et[ijet],1);     //smear jet with 10% gaussian 90% of the time
	  else { evt.jet_et[ijet] = smear_jet_et(evt.jet_et[ijet],8);  smear_flag++; }                        //smear jet downward with 50% gaussian 10% of the time
	}
	if(jet_smear_flag==2) {
	  if(gRandom->Uniform(0,1)<0.5) evt.jet_et[ijet] = smear_jet_et(evt.jet_et[ijet],1);     //smear jet with 10% gaussian 50% of the time
	  else { evt.jet_et[ijet] = smear_jet_et(evt.jet_et[ijet],8);  smear_flag++; }                        //smear jet downward with 50% gaussian 50% of the time
	}
	//cout<<"\n"<<temp_pt<<" "<<evt.jet_et[ijet];
	if(evt.jet_et[ijet]<jet_pt_min) evt.jet_et[ijet]=temp_pt;   //reset jet to original pt if smeared below threshold   
	//if(evt.jet_et[ijet]<jet_pt_min) continue;                          
      }  //Jet smearing
      /*
      jet_smear_et[ijet][1]= smear_jet_et(evt.jet_et[ijet],1); //saved smeared jet et
      jet_smear_et[ijet][2]= smear_jet_et(evt.jet_et[ijet],2); //saved smeared jet et
      jet_smear_et[ijet][3]= smear_jet_et(evt.jet_et[ijet],3); //saved smeared jet et
      jet_smear_et[ijet][4]= smear_jet_et(evt.jet_et[ijet],4); //saved smeared jet et
      jet_smear_et[ijet][5]= smear_jet_et(evt.jet_et[ijet],5); //saved smeared jet et
      jet_smear_et[ijet][6]= smear_jet_et(evt.jet_et[ijet],6); //saved smeared jet et
      jet_smear_et[ijet][7]= smear_jet_et(evt.jet_et[ijet],0); //saved true jet et
      for(int i=0; i<8; i++) jet_smear_phi[ijet][i] = evt.jet_phi[ijet]; //save jet phi
      */
    } //loop over jets to smear them
    //if(!smear_flag) continue; 


    /*
    //For event dump
    if(evt.event==43711304 && evt.run==136066){
      //////// //event dump
      float ht_vect[2] = {0,0};
      evt2=&evt;
      event_dump(evt2,0,ht_vect,0,DiMuon);
      //////////
      dump_counter++;
      //if(dump_counter>5) exit(0);
    }
    */


    int temp_numjets = 0;
    int num_bjets = 0;
    float temp_emf=-999;
    float temp_phoet, temp_jetet, temp_pho_ecaliso, temp_overlap_jetpt[2];
    if(jet_quality==1){ //impose jet quality cuts
      int removed_jet = 0;
      for(ijet=0; ijet<evt.num_jet; ijet++)  {//loop over jets  
	if(pho_sample==1){  //make sure photon isn't same as jet
	  float dphi = deltaphi(pho_tr.pho_phi[pho_index],evt.jet_phi[ijet]);
	  float deta = fabs(evt.jet_eta[ijet] - pho_tr.pho_eta[pho_index]);
	  float dR = sqrt(pow(dphi,2)+pow(deta,2));
	  //cout<<"\n"<<pho_tr.pho_phi[pho_index]<<" "<<evt.jet_phi[ijet]<<" "<<evt.jet_eta[ijet]<<" "<<pho_tr.pho_eta[pho_index];
	  //cout<<"\n"<<dR;
	  if(dR < 0.2) {
	    evt.jet_remove[ijet]=1;  /*removed_jet=1;*/ 
	    //h_phopt_vs_jetpt_overlap->Fill(pho_tr.pho_et[pho_index],evt.jet_et[ijet],weight2);
	    temp_phoet = pho_tr.pho_et[pho_index]; temp_jetet = evt.jet_et[ijet]; 
	    //h_pho_isoEcalRecHitDR04_overlap->Fill(pho_tr.pho_isoEcalRecHitDR04[pho_index],weight1);
	    temp_pho_ecaliso = pho_tr.pho_isoEcalRecHitDR04[pho_index];
	    temp_emf = evt.jet_emf[ijet];
	    temp_overlap_jetpt[0] = evt.jet_px[ijet]; temp_overlap_jetpt[1] = evt.jet_py[ijet]; 
            //h_jet_emf_overlap->Fill(evt.jet_emf[ijet],weight1);
	    //pho_tr.pho_phi[pho_index] = evt.jet_phi[ijet]; pho_tr.pho_et[pho_index] = evt.jet_et[ijet]; pho_tr.pho_eta[pho_index] = evt.jet_eta[ijet]; //set photon equal to jet
	    //DiMuon.SetXYZM(evt.jet_px[ijet],evt.jet_py[ijet],evt.jet_pz[ijet],0);
	    continue;
	  } //photon overlaps jet
	}
	if(evt.jet_et[ijet]<jet_pt_min)   { evt.jet_remove[ijet]=1; continue;} //remove jets below pt cut
	if(evt.jet_et[ijet]>jet_pt_min) { 
	  if(fabs(evt.jet_eta[ijet]) > 3.0 || evt.jet_fHPD[ijet]>0.98 || evt.jet_n90Hits[ijet]<=1 || (evt.jet_emf[ijet]<0.01 && fabs(evt.jet_eta[ijet])<2.6 )) { evt.jet_remove[ijet]=1; removed_jet=1; continue;}  //calo jets
	  //if(fabs(evt.jet_eta[ijet]) > 2.0 || evt.jet_fHPD[ijet]>0.98 || evt.jet_n90Hits[ijet]<=1 /*|| evt.jet_emf[ijet]<0.01*/) { evt.jet_remove[ijet]=1; removed_jet=1; continue;}  //Jpt jets	
	  //if(fabs(evt.jet_eta[ijet]) > 3.0 || evt.jet_fHPD[ijet]/evt.jet_E[ijet]>=0.99 || evt.jet_fHPD[ijet]/evt.jet_E[ijet]>=0.99 || evt.jet_emf[ijet]/evt.jet_E[ijet]>=0.99 || ((evt.jet_n90Hits[ijet]<=0.01 || evt.jet_btag_jetProb[ijet]<=0.01) && fabs(evt.jet_eta[ijet])<2.4 )) { evt.jet_remove[ijet]=1; removed_jet=1; continue;}  //PF jets
	}
	if(evt.jet_remove[ijet]==0) temp_numjets++;	
	if(evt.jet_btag_jetProb[ijet]>btag_cut) num_bjets++;
      }
      if(removed_jet) continue; //throw out event if any jets fail quality cuts
      if(bjet_sample && num_bjets==0) continue;  //throw out event if no b-jets
    }


    
    if(mode==-1) { //Set weightings for histograms; No Z mass cut
      weight1=weight2=1; 
      if(sample==-5) { //QCD stitching
	if(ia<=9447272) weight1=weight2 = 1.0; 
	if(ia>=9447273 && ia<=14315150) weight1=weight2 = 0.0526694; 
	if(ia>=14315151 && ia<=18550848) weight1=weight2 = 0.00183962; 
	if(ia>=18550849) weight1=weight2 = 0.000075; 
      }
    }
   
    if(mode==-2) weight1=weight2=1; 
 
    /*
    if(mode==-1) { //Set weightings for histograms; No Z mass cut
      weight1=weight2=1; 
      if(sample==-5) { //QCD stitching 
	if(ia<=9447272) weight1=weight2 = 13333.3; 
	if(ia>=9447273 && ia<=14315150) weight1=weight2 = 702.2587; 
	if(ia>=14315151 && ia<=18550848) weight1=weight2 = 24.5283;
 	if(ia>=18550849) weight1=weight2 = 1.0; 
      }
    }
    */


    //Fill jet histograms
    //if(temp_numjets!=min_num_jets) continue;
    //if(temp_numjets<=4) continue;
    for(ijet=0; ijet<evt.num_jet; ijet++)  {//loop over jets  
      if((temp_numjets<2 && pho_sample==0) || (temp_numjets<1 && pho_sample==1)) continue;
      if(evt.jet_remove[ijet]) continue;
      h_jet_pt->Fill(evt.jet_et[ijet],weight1);
      h_jet_pt_wide->Fill(evt.jet_et[ijet],weight1);
      h_jet_phi->Fill(evt.jet_phi[ijet],weight1);
      h_jet_eta->Fill(evt.jet_eta[ijet],weight1);
      h_jet_fRBX->Fill(evt.jet_fRBX[ijet],weight1);
      h_jet_fHPD->Fill(evt.jet_fHPD[ijet],weight1);
      //h_jet_fRBX->Fill(evt.jet_fRBX[ijet]/evt.jet_E[ijet],weight1);
      //h_jet_fHPD->Fill(evt.jet_fHPD[ijet]/evt.jet_E[ijet],weight1);
      h_jet_n90Hits->Fill(evt.jet_n90Hits[ijet],weight1);
      h_jet_emf->Fill(evt.jet_emf[ijet],weight1);
      h_jet_btag_jetProb->Fill(evt.jet_btag_jetProb[ijet],weight1);
      //h_jet_emf->Fill(evt.jet_emf[ijet]/evt.jet_E[ijet],weight1);
      //h_jet_btag_jetProb->Fill(evt.jet_btag_jetProb[ijet]/evt.jet_E[ijet],weight1);
      if(evt.jet_btag_jetProb[ijet]>btag_cut) h_jet_pt_btag->Fill(evt.jet_et[ijet],weight1);
    }
    //Fill photon histograms
    if(pho_sample==1) {
      if(temp_numjets<1) continue;
      h_pho_pt->Fill(DiMuon.Pt(),weight1);
      h_pho_phi->Fill(DiMuon.Phi(),weight1);
      h_pho_eta->Fill(pho_tr.pho_eta[pho_index],weight1);
      h_pho_E1overE3x3->Fill(pho_tr.pho_maxEnergyXtal[pho_index]/pho_tr.pho_e3x3[pho_index],weight1);
      h_pho_sigmaIetaIeta->Fill(pho_tr.pho_sigmaIetaIeta[pho_index],weight1);
      h_pho_hasPixelSeed->Fill(pho_tr.pho_hasPixelSeed[pho_index],weight1);
      h_pho_hadOverEM->Fill(pho_tr.pho_hadOverEM[pho_index],weight1);
      h_pho_isoEcalRecHitDR04->Fill(pho_tr.pho_isoEcalRecHitDR04[pho_index],weight1);
      h_pho_isoHcalRecHitDR04->Fill(pho_tr.pho_isoHcalRecHitDR04[pho_index],weight1);
      h_pho_isoHollowTrkConeDR04->Fill(pho_tr.pho_isoHollowTrkConeDR04[pho_index],weight1);
      if(temp_emf>-900) {
	//h_jet_emf_overlap->Fill(temp_emf,weight1);
	//h_phopt_vs_jetpt_overlap->Fill(temp_phoet,temp_jetet,weight2);
	h_pho_isoEcalRecHitDR04_overlap->Fill(temp_pho_ecaliso,weight1);
      }
    }
    h_jet_num->Fill(temp_numjets,weight1);


    //cout<<"\n"<<DiMuon.Pt()<<" "<<DiMuon.Phi()<<" ";

     //smear Z
    if(z_smear==1 || z_smear==3){ //smear Z pt
      float temp_pt = 0;
      float pt_old = DiMuon.Pt();
      //temp_pt = smear_jet_et(pt_old,1); //saved smeared et
      temp_pt = smear_jet_et(pt_old,2); //saved smeared et
      //temp_pt = smear_jet_et(pt_old,3); //saved smeared et
      //temp_pt = smear_jet_et(pt_old,4); //saved smeared et
      //temp_pt = smear_jet_et(pt_old,5); //saved smeared et
      //temp_pt = smear_jet_et(pt_old,6); //saved smeared et
      float px,py,pz,mass;
      mass = DiMuon.M();
      pz = DiMuon.Pz(); px = DiMuon.Px(); py = DiMuon.Py();
      px *= temp_pt/pt_old; py *= temp_pt/pt_old; 
      DiMuon.SetXYZM(px,py,pz,mass);
      DiMuon_smear[0] = DiMuon;
    }
    if(z_smear==2 || z_smear==3) { //smear Z phi
      float temp_phi, phi_old;
      phi_old = DiMuon.Phi();
      temp_phi = smear_phi(phi_old,1); //saved smeared phi 
      //temp_phi = smear_phi(phi_old,2); //saved smeared phi 
      //temp_phi = smear_phi(phi_old,3); //saved smeared phi 
      //temp_phi = smear_phi(phi_old,4); //saved smeared phi 
      float temp_pt,temp_eta,temp_E;
      temp_pt = DiMuon.Pt();
      temp_eta = DiMuon.Eta();
      temp_E = DiMuon.E();
      DiMuon.SetPtEtaPhiE(temp_pt,temp_eta,temp_phi,temp_E);
      DiMuon_smear[0] = DiMuon;
    }

    //cout<<DiMuon.Pt()<<" "<<DiMuon.Phi();


    //Z mass cuts    
    if(mode==0) {
      if(DiMuon.M()>81 && DiMuon.M()<101) weight1=weight2=1; //main mass window
      else continue;
    }
    if(mode==1) {
      if(DiMuon.M()>81 && DiMuon.M()<101) weight1=weight2=1; //main mass window
      else if((DiMuon.M()>71 && DiMuon.M()<81) || (DiMuon.M()>101 && DiMuon.M()<111)) { //sidebands
	weight1=-1; weight2=0;
      }
      else continue;
    }
    if(mode==2) {
      if(DiMuon.M()>81 && DiMuon.M()<101) weight1=weight2=1; //main mass window
      else if((DiMuon.M()>61 && DiMuon.M()<71) || (DiMuon.M()>111 && DiMuon.M()<121)) { //sidebands
	weight1=-1; weight2=0;
      }
      else continue;
    }



 
    //jet removal for n-jet events
    if(test) { 
      int num_removed = 0; 
      for(ijet=0; ijet<evt.num_jet; ijet++) { //loop through all jets
	if(evt.jet_remove[ijet]==1 || evt.jet_et[ijet]<jet_pt_min) continue;
	//if(gRandom->Uniform()>0.80) { evt.jet_remove[ijet]=1; num_removed++;} //remove 20% of all jets 
	//if(gRandom->Uniform()>0.75) { evt.jet_remove[ijet]=1; num_removed++;} //remove 25% of all jets 
	if(test==1) if(gRandom->Uniform()>0.67) { evt.jet_remove[ijet]=1; num_removed++;} //remove 33% of all jets
	//if(gRandom->Uniform()>0.55) { evt.jet_remove[ijet]=1; num_removed++;} //remove 45% of all jets
	//if((evt.jet_et[ijet]-jet_pt_min)/(200.-2*jet_pt_min)+0.5<gRandom->Uniform())  { evt.jet_remove[ijet]=1; num_removed++;} //remove 50% of jets at threshold to 0% of jets at 100 GeV
	if(test==2) if((evt.jet_et[ijet]-jet_pt_min)/(400.-2*jet_pt_min)+0.5<gRandom->Uniform())  { evt.jet_remove[ijet]=1; num_removed++;} //remove 50% of jets at threshold to 0% of jets at 200 GeV
	//if(test==2) if((evt.jet_et[ijet]-jet_pt_min)/(200.-jet_pt_min)*(1.0-0.35)+0.35<gRandom->Uniform())  { evt.jet_remove[ijet]=1; num_removed++;} //remove 25% of jets at threshold to 0% of jets at 200 GeV
	//if((evt.jet_et[ijet]-jet_pt_min)/(100.-2*jet_pt_min)+0.5<gRandom->Uniform())  { evt.jet_remove[ijet]=1; num_removed++;} //remove 50% of jets at threshold to 0% of jets at 50 GeV
	//if((evt.jet_et[ijet]-jet_pt_min)/(2*500.-2*jet_pt_min)+0.5<gRandom->Uniform())  { evt.jet_remove[ijet]=1; num_removed++;} //remove 50% of jets at threshold to 0% of jets at 500 GeV
	//if(fabs(evt.jet_phi[ijet]) < 1.5 && fabs(evt.jet_eta[ijet])<1.0) { evt.jet_remove[ijet]=1; num_removed++;}  //remove jets with eta in [-1,1] and phi in [0,1.5] 
	//if(evt.jet_et[ijet]<100 && evt.jet_et[ijet]>40 &&  0.5<gRandom->Uniform()) { evt.jet_remove[ijet]=1; num_removed++;} //remove 50% of 40<pt<100 jets
      }
      if(num_removed!=1) continue; //Only save events with exactly 1 above threshold jet removed
    }




    //N-Jet code
    int num_passing_jets=0; 
    float jet_pt[evt.num_jet+1][2];  //second index is for px and py
    float jet_phi[evt.num_jet+1];
    float eta=-999;
    float njet_ht = 0;
    float njet_mht = 0;
    float njet_htvec[2] = {0,0}; 
    float njet_deltaht = 10000;
    float njet_2ndbestdeltaht = 10000;
    float njet_alphat=-999.;
    njet_ht += DiMuon.Pt(); //include Z as one of the jets
    njet_htvec[0] += DiMuon.Px();  njet_htvec[1] += DiMuon.Py(); 
    jet_pt[num_passing_jets][0] = DiMuon.Px();
    jet_pt[num_passing_jets][1] = DiMuon.Py();
    num_passing_jets++;
    if(n_jet<0) num_passing_jets--; //For exclusive jets;  This command overwrites the Z
    for(ijet=0; ijet<evt.num_jet; ijet++) { //loop over jets to save good jets
      if(evt.jet_et[ijet]>jet_pt_min && evt.jet_remove[ijet]==0) { 
	jet_pt[num_passing_jets][0] = evt.jet_px[ijet]; //cout<<"\nBeginning: "<<jet_pt[num_passing_jets][0];
	jet_pt[num_passing_jets][1] = evt.jet_py[ijet]; //cout<<" "<<jet_pt[num_passing_jets][1]<<" "<<evt.jet_et[ijet]<<" "<<num_passing_jets;
	jet_phi[num_passing_jets] = evt.jet_phi[ijet];
	if(eta<-900) eta = evt.jet_eta[ijet];
	num_passing_jets++;
	njet_ht += evt.jet_et[ijet]; 
	njet_htvec[0] += evt.jet_px[ijet];  njet_htvec[1] += evt.jet_py[ijet]; 
      }
    }
    h_ht->Fill(njet_ht,weight1);
    h_mht->Fill(sqrt(pow(njet_htvec[0],2)+pow(njet_htvec[1],2)),weight1);
    h_met->Fill(evt.met);
    h_genht->Fill(evt.rawmet,weight1);
    h_ht_vs_htplusmet->Fill(njet_ht+evt.met,njet_ht,weight2);
    p_ht_vs_htplusmet->Fill(njet_ht+evt.met,njet_ht,weight1);
    h_ht_vs_genht->Fill(evt.rawmet,njet_ht,weight2);
    p_ht_vs_genht->Fill(evt.rawmet,njet_ht,weight1); 
    /*
    //To study the strange properties of photons in the early data
    //They have lower energy than the jets in pho+1jet events
    if(num_passing_jets==2) { h_phojet_balance->Fill(sqrt(pow(jet_pt[0][0],2)+pow(jet_pt[0][1],2))/sqrt(pow(jet_pt[1][0],2)+pow(jet_pt[1][1],2)),weight1);
                              if(gRandom->Uniform(0,1)>0.5) h_dijet_balance->Fill(sqrt(pow(jet_pt[0][0],2)+pow(jet_pt[0][1],2))/sqrt(pow(jet_pt[1][0],2)+pow(jet_pt[1][1],2)),weight1);
                              else h_dijet_balance->Fill(sqrt(pow(jet_pt[1][0],2)+pow(jet_pt[1][1],2))/sqrt(pow(jet_pt[0][0],2)+pow(jet_pt[0][1],2)),weight1);
			      float met_vec[2] = {evt.met*cos(evt.met_phi),evt.met*sin(evt.met_phi)};
			      met_vec[0] = met_vec[0] + temp_overlap_jetpt[0] - DiMuon.Px();
			      met_vec[1] = met_vec[1] + temp_overlap_jetpt[1] - DiMuon.Py();
			      evt.met = sqrt(pow(met_vec[0],2)+pow(met_vec[1],2));
			      evt.met_phi = findphi(met_vec[0],met_vec[1]);
			      h_met->Fill(evt.met,weight1);
			      //if(evt.met>35){
			      //h_jet_emf_overlap->Fill(temp_emf,weight1);
			      //h_phopt_vs_jetpt_overlap->Fill(temp_phoet,temp_jetet,weight2);
			      //h_dphi_phomet->Fill(fabs(deltaphi(evt.met_phi,pho_tr.pho_phi[pho_index])),weight1);
			      //}
			      h_dphi_phomet->Fill(fabs(deltaphi(evt.met_phi,pho_tr.pho_phi[pho_index])),weight1);
			      if(evt.met<20) h_dphi_phomet_met0to20->Fill(fabs(deltaphi(evt.met_phi,DiMuon.Phi())),weight1);
			      if(evt.met>20) h_dphi_phomet_met20toinf->Fill(fabs(deltaphi(evt.met_phi,DiMuon.Phi())),weight1);
    }
    */
    //cout<<"\n"<<njet_mht;
    //cout<<"\n"<<num_passing_jets;
    njet_mht = sqrt(pow(njet_htvec[0],2) + pow(njet_htvec[1],2)); 
    if(num_passing_jets<=min_num_jets) continue; //only consider events with a Z and minimum number of jets
    if(num_passing_jets>max_num_jets+1) continue; //only consider events with a Z and less than or equal to maximum number of jets
    num_evts++;
    //Smear a jet ****
    /*
    if(jet_smear_flag) {
      float temp_pt[2], smeared_pt;
      int smear_index=0; //Smear leading jet for now
      temp_pt[0] = jet_pt[smear_index][0]; temp_pt[1] = jet_pt[smear_index][1]; 
      smeared_pt = smear_jet_et(sqrt(pow(temp_pt[0],2)+pow(temp_pt[1],2)),6); //compute new pt
      smeared_pt = sqrt(pow(temp_pt[0],2)+pow(temp_pt[1],2))*0.75;  //scale the pT of lead jet
      if(smeared_pt<jet_pt_min) continue;  //Make sure smeared_pt is above jet threshold
      //smeared_pt = 10;
      jet_pt[smear_index][0] = smeared_pt*cos(jet_phi[smear_index]);
      jet_pt[smear_index][1] = smeared_pt*sin(jet_phi[smear_index]);     
      njet_ht = njet_ht - sqrt(pow(temp_pt[0],2)+pow(temp_pt[1],2)) + smeared_pt; //recompute HT
      njet_htvec[0] = njet_htvec[0] - temp_pt[0] + jet_pt[smear_index][0];        
      njet_htvec[1] = njet_htvec[1] - temp_pt[1] + jet_pt[smear_index][1];
      njet_mht = sqrt(pow(njet_htvec[0],2) + pow(njet_htvec[1],2));               //recompute MHT
    }
    */
    //End of smear a jet ****
    //soft jet stuff ****
    if(soft) { //Add soft jets right after N-jet system has been found and before computing DeltaHT
      int nsoft; //number of soft jets to add
      int softadded = 0; //number of soft jets added 
      //nsoft = (int) gRandom->Exp(1.5);
      //nsoft = (int) gRandom->Exp(2);
      if(soft==1) nsoft = (int) gRandom->Exp(3); //1,3,4
      if(soft==2) nsoft = (int) gRandom->Exp(4);
      //nsoft = (int) gRandom->Exp(10); //6
      //nsoft = (int) gRandom->Exp(pow((ht_smear[i]-60)/30.,3)); //5
      while(softadded<nsoft) { //add multiple soft jets
	//determine soft jet pT and phi
	float soft_phi = gRandom->Uniform(0,2*pi);
	float soft_pt;
	//soft_pt = leadjet_pt_min; //1
	//soft_pt = gRandom->Exp(4);
	if(soft==1 && sample==7) soft_pt = gRandom->Exp(10);
	if(soft==1 && sample==-3) soft_pt = gRandom->Exp(30);
	if(soft==2 && sample==7) soft_pt = gRandom->Exp(15); //3,6
	if(soft==2 && sample==-3) soft_pt = gRandom->Exp(20);
	//soft_pt = gRandom->Exp(30); 
	//soft_pt = gRandom->Gaus(30,10); //4
	//soft_pt = gRandom->Gaus(pow((ht_smear[i]-60)/26.,4)/*ht_smear[i]/4.*/,1);  if(ht_smear[i]>120) soft_pt=30; //5
	if(soft_pt>jet_pt_min) continue; //make sure jet is soft
	else softadded++;
	float temp_htvec[2] = {0,0};
	float temp_ht=0;
	float temp_mass;
	temp_mass = njet_ht + njet_mht;
	TLorentzVector object[num_passing_jets+1];
	object[0].SetXYZM(-njet_htvec[0],-njet_htvec[1],0,0); //MHT
	for(int j=0; j< num_passing_jets; j++) { //Save passing jets as Lorentz vectors
	  object[j+1].SetXYZM(jet_pt[j][0],jet_pt[j][1],0,0);
	}
	TLorentzVector s,t;
	s.SetXYZM(soft_pt*cos(soft_phi),soft_pt*sin(soft_phi),0,0);
	//t finds the boost of the system after the soft jet
	t.SetXYZM(s.Px(),s.Py(),s.Pz(),temp_mass);
	for(int j=1; j< num_passing_jets+1; j++) { //boost jets and save their momenta
	  object[j].Boost(-t.BoostVector());
	  jet_pt[j-1][0] = object[j].Px();
	  jet_pt[j-1][1] = object[j].Py(); 
	  temp_htvec[0] += jet_pt[j-1][0];
	  temp_htvec[1] += jet_pt[j-1][1];
	  temp_ht += object[j].Pt();
	}
	njet_ht = temp_ht;
	njet_htvec[0] = temp_htvec[0];
	njet_htvec[1] = temp_htvec[1];
	njet_mht = sqrt(pow(njet_htvec[0],2) + pow(njet_htvec[1],2)); 
      } //add multiple soft jets
      if(soft/*==3*/) { //tighten Z and jet cuts
	for(int j=0; j< num_passing_jets; j++) { 
	  if(sqrt(pow(jet_pt[j][0],2)+pow(jet_pt[j][1],2)) < jet_pt_min/*+10*/) continue; //jet fails tighter pT
	}
      }
    }
    //End of soft jet stuff ****
    //Now Z and all good jets are stored in jet_pt which has num_passing_jets entries 
    //Try all combinations of 2 pseudojets to minimize deltaHT 
    float temp_deltaht;
    for(int n1=1; n1<= num_passing_jets/2; n1++) { //number of jets in first pseudojet 
      int end_flag = 0;
      int jet_index[n1];
      for(int a=0; a<n1; a++) jet_index[a] = n1-1-a; 
      jet_index[0] -= 1;
      while(end_flag==0){ //search all combinations for lowest deltaHT value
	jet_index[0]++;
	if(jet_index[0] == num_passing_jets){
	  int carry=0;
	  for(int a=1; a<n1; a++) {
	    if(jet_index[a] < num_passing_jets-1-a) { //find closest digit than can be incremented
	      carry = a;
	      jet_index[a]++;
	    }
	  }
	  if(carry==0) {end_flag=1; continue;} //no more combinations
	  for(int a=0; a<carry; a++) { //reset all previous digits; avoid double counting
	    jet_index[a] = jet_index[carry]+(carry-a); 
	  }
	}
	//Compute deltaHT for given combination
	float jet1_ht[2] = {0,0};
	//cout<<"\n"<<num_passing_jets<<" ";
	for(int a=0; a<n1; a++) { 
	  //cout<<jet_index[a]<<" ";
	  jet1_ht[0] += jet_pt[jet_index[a]][0];
	  jet1_ht[1] += jet_pt[jet_index[a]][1];
	}
	float jet2_ht[2];
	jet2_ht[0] = njet_htvec[0] - jet1_ht[0];
	jet2_ht[1] = njet_htvec[1] - jet1_ht[1];
	temp_deltaht = fabs(sqrt(pow(jet1_ht[0],2)+pow(jet1_ht[1],2)) - sqrt(pow(jet2_ht[0],2)+pow(jet2_ht[1],2)));
	if(temp_deltaht < njet_deltaht) { //save new deltaHT if it is smaller than previous smallest value 
	  njet_2ndbestdeltaht = njet_deltaht; //save 2nd best value
	  njet_deltaht = temp_deltaht;   	    
	}
      }
    }
    //njet_deltaht = njet_2ndbestdeltaht;  //use second lowest deltaHT instead
    njet_alphat = 0.5*(njet_ht-njet_deltaht)/sqrt(pow(njet_ht,2)-pow(njet_mht,2));
    offset = jet_pt_min*(min_num_jets) + z_pt_min; 
    bin_width = 10.0;
    if(sample==-5) bin_width = 10;
    max_ht = max_htbins*bin_width+offset;
    //if(sample==0) { bin_width = 10.0; max_ht = max_htbins*bin_width+offset; }
    //cout<<"\n"<<njet_ht<<" "<<njet_mht<<" "<<njet_deltaht<<" "<<njet_alphat;
    //njet_ht+= evt.met; //include MET into HT **********
    //njet_ht = evt.rawmet; //Use genHT for HT   **********
    h_alphaT->Fill(njet_alphat,weight1);
    if(njet_ht>80) h_alphaT_ht_80_up->Fill(njet_alphat,weight1);
    if(njet_ht>40) h_alphaT_ht_40_up->Fill(njet_alphat,weight1);
    if(njet_ht>120) h_alphaT_ht_120_up->Fill(njet_alphat,weight1);
    if(njet_ht>160) h_alphaT_ht_160_up->Fill(njet_alphat,weight1);
    if(njet_ht>40 && njet_ht<80) h_alphaT_ht_40_80->Fill(njet_alphat,weight1);
    if(njet_ht>20 && njet_ht<40) h_alphaT_ht_20_40->Fill(njet_alphat,weight1); 
    if(njet_ht>30 && njet_ht<50) h_alphaT_ht_30_50->Fill(njet_alphat,weight1);
    if(njet_ht>50 && njet_ht<80) h_alphaT_ht_50_80->Fill(njet_alphat,weight1);
    if(njet_ht>80 && njet_ht<120) h_alphaT_ht_80_120->Fill(njet_alphat,weight1);
    if(njet_ht>120 && njet_ht<160) h_alphaT_ht_120_160->Fill(njet_alphat,weight1);
    if(njet_ht<max_ht && njet_ht>offset) num_total_smear[0][(int)((njet_ht-offset)/bin_width)]+=weight1;
    if(njet_ht<max_ht && njet_ht>offset) if(weight1>0) num_total_smear_noscale[0][(int)((njet_ht-offset)/bin_width)]+= 1; 
    if(njet_ht<max_ht && njet_ht>offset) x_sum_smear[0][(int)((njet_ht-offset)/bin_width)]+= weight1*njet_ht;
    if(njet_ht<max_ht && njet_ht>offset) x2_sum_smear[0][(int)((njet_ht-offset)/bin_width)]+= weight1*pow(njet_ht,2);
    if(njet_alphat>0.55 && njet_ht<max_ht && njet_ht>offset) num_above_smear[0][(int)((njet_ht-offset)/bin_width)]+=weight1; 
    /* //for fraction vs Z pt
       if(ht_smear[i]<max_ht && ht_smear[i]>offset) num_total_smear[i][(int)((DiMuon.Pt()-offset)/bin_width)]+=weight1;
       if(ht_smear[i]<max_ht && ht_smear[i]>offset) x_sum_smear[i][(int)((DiMuon.Pt()-offset)/bin_width)]+= weight1*DiMuon.Pt();
       if(ht_smear[i]<max_ht && ht_smear[i]>offset) x2_sum_smear[i][(int)((DiMuon.Pt()-offset)/bin_width)]+= weight1*DiMuon.Pt()*DiMuon.Pt();
       if(alphaT>0.55 && ht_smear[i]<max_ht && ht_smear[i]>offset) num_above_smear[i][(int)((DiMuon.Pt()-offset)/bin_width)]+=weight1; 
    */
     //for flatness in eta
    eta= fabs(eta);
    if(fabs(eta)<3. && njet_ht>80 && njet_ht<120 ) {   //cout<<"\n"<<x2_sum_eta[0][0][2]<<" "<<x_sum_eta[0][0][2]<<" "<<num_total_eta[0][0][2]<<" "<<weight1*fabs(eta)<<" "<<(int)fabs(eta); 
    //    if(fabs(eta)<3. && njet_ht>120 && njet_ht<180 ) {
      if(njet_alphat>0.55) num_above_eta[0][0][(int)fabs(eta)]+=weight1; num_total_eta[0][0][(int)fabs(eta)]+=weight1; x_sum_eta[0][0][(int)fabs(eta)]+=weight1*fabs(eta);  x2_sum_eta[0][0][(int)fabs(eta)]+=weight1*eta*eta;
    }
    
    if(eta<3. && njet_ht>120 && njet_ht<160) { //cout<<"\n"<<x2_sum_eta[0][1][2]<<" "<<x_sum_eta[0][1][2]<<" "<<num_total_eta[0][1][2]<<" "<<num_above_eta[0][1][2]<<" "<<weight1*eta<<" "<<(int)eta;
    //  if(eta<3. && njet_ht>180 ) {
      if(njet_alphat>0.55) num_above_eta[0][1][(int)eta]+=weight1; num_total_eta[0][1][(int)eta]+=weight1; x_sum_eta[0][1][(int)eta]+=weight1*eta;  x2_sum_eta[0][1][(int)eta]+=weight1*eta*eta;
    }
    if(eta<3. && njet_ht>160) {
//  if(eta<3. && njet_ht>200) {
      if(njet_alphat>0.55) num_above_eta[0][2][(int)eta]+=weight1; num_total_eta[0][2][(int)eta]+=weight1; x_sum_eta[0][2][(int)eta]+=weight1*eta;  x2_sum_eta[0][2][(int)eta]+=weight1*eta*eta;
    }
    if(eta<3. && njet_ht>160 && njet_ht<300) {
      if(njet_alphat>0.55) num_above_eta[0][3][(int)eta]+=weight1; num_total_eta[0][3][(int)eta]+=weight1; x_sum_eta[0][3][(int)eta]+=weight1*eta;  x2_sum_eta[0][3][(int)eta]+=weight1*eta*eta;
    }
        

    h_measuredht->Fill(njet_ht,weight1);
    if(njet_alphat>0.55) h_measuredht_failalphat->Fill(njet_ht,weight1);
    h_Npv->Fill(evt.num_vert,weight1);
    h_Npv_good->Fill(evt.num_mu,weight1); //using num_mu temporarily to hold this information

    /*
    //For event dump
    if(njet_alphat>0.55 && njet_ht>180){
      //////// //event dump
      float ht_vect[2] = {njet_htvec[0],njet_htvec[1]};
      evt2=&evt;
      event_dump(evt2,0,ht_vect,njet_ht,DiMuon);
      //////////
      dump_counter++;
      //if(dump_counter>5) exit(0);
    }
    */    



  }//end of event loop


  //cout<<"\n\n"<<num_evts<<"\n\n";

  //Check the failure fraction entries
  //cout<<"\n\n";
  //for(int a=0; a<8; a++) cout<<num_total_smear[0][a]<<" "<<num_above_smear[0][a]<<" "<<x_sum_smear[0][a]<<" "<<x2_sum_smear[0][a]<<"\n";
  //exit(0);



  TCanvas* c1= new TCanvas();  

  gStyle->SetOptStat(1111111);
  gStyle->SetOptTitle(1);


  /*
   //don't need this for early data
  //compute eta relate ratios and errors
  for(int i=0; i<1; i++) { //loop through different smearings 
    for(int j=0; j<3; j++) { //loop through HT bins
      for(int k=0; k<3; k++) { //loop through eta bins
	//if(num_above_eta[i][j][k]<1) num_above_eta[i][j][k]=1;
	ratio_eta[i][j][k] = num_above_eta[i][j][k]/num_total_eta[i][j][k]; 
	//error_eta[i][j][k] = ratio_eta[i][j][k] *  sqrt(1./num_above_eta[i][j][k] + 1./num_total_eta[i][j][k]);
	error_eta[i][j][k] = sqrt(ratio_eta[i][j][k]*(1-ratio_eta[i][j][k])/num_total_eta[i][j][k]);  //binomial errors
	xerr_eta[i][j][k] = sqrt(x2_sum_eta[i][j][k]/num_total_eta[i][j][k]-pow(x_sum_eta[i][j][k]/num_total_eta[i][j][k],2));
	x_eta[i][j][k] = x_sum_eta[i][j][k]/num_total_eta[i][j][k];
	//error_eta[i][j][k] *= sqrt(8.); //Increase error for 1/8th sample size (100pb-1)
      } //loop through eta bins
  
      
      //if(i==0){
      if(0){ //What is this for?  Bins of HT?
	TGraphAsymmErrors *g_etavsfracabove55 = new TGraphAsymmErrors(3,x_eta[i][j],ratio_eta[i][j],xerr_eta[i][j],xerr_eta[i][j],error_eta[i][j],error_eta[i][j]);
	g_etavsfracabove55->SetTitle("Fraction with #alpha_{T}>0.55 vs Eta"); 
	gStyle->SetOptFit();
	c1->SetLogy(0); 
	g_etavsfracabove55->Draw("ALP"); if(j==0) c1->Print("g_etavsfracabove55_a_1.pdf"); if(j==1) c1->Print("g_etavsfracabove55_b_1.pdf"); 
	if(j==2) c1->Print("g_etavsfracabove55_c_1.pdf"); if(j==3) c1->Print("g_etavsfracabove55_d_1.pdf"); 
	c1->SetLogy(); 
	g_etavsfracabove55->Draw("ALP"); if(j==0) c1->Print("g_etavsfracabove55_a_2.pdf"); if(j==1) c1->Print("g_etavsfracabove55_b_2.pdf"); 
	if(j==2) c1->Print("g_etavsfracabove55_c_2.pdf"); if(j==3) c1->Print("g_etavsfracabove55_d_2.pdf"); 
	gStyle->SetOptFit(0);
      }
      
         
    } //loop through HT bins

  } //loop through different smearings 

  TGraphAsymmErrors *g0_etavsfracabove55 = new TGraphAsymmErrors(3,x_eta[0][0],ratio_eta[0][0],xerr_eta[0][0],xerr_eta[0][0],error_eta[0][0],error_eta[0][0]);
  g0_etavsfracabove55->SetTitle("Fraction(#alpha_{T}>0.55)  vs |#eta|");
  g0_etavsfracabove55->SetMarkerStyle(22); g0_etavsfracabove55->SetMarkerColor(1); g0_etavsfracabove55->SetMarkerSize(2);
  TGraphAsymmErrors *g1_etavsfracabove55 = new TGraphAsymmErrors(3,x_eta[0][1],ratio_eta[0][1],xerr_eta[0][1],xerr_eta[0][1],error_eta[0][1],error_eta[0][1]);
  g1_etavsfracabove55->SetMarkerStyle(22); g1_etavsfracabove55->SetMarkerColor(2);  g1_etavsfracabove55->SetMarkerSize(2);
  TGraphAsymmErrors *g2_etavsfracabove55 = new TGraphAsymmErrors(3,x_eta[0][2],ratio_eta[0][2],xerr_eta[0][2],xerr_eta[0][2],error_eta[0][2],error_eta[0][2]);      
  g2_etavsfracabove55->SetMarkerStyle(22); g2_etavsfracabove55->SetMarkerColor(3);  g2_etavsfracabove55->SetMarkerSize(2);
  TGraphAsymmErrors *g3_etavsfracabove55 = new TGraphAsymmErrors(3,x_eta[0][3],ratio_eta[0][3],xerr_eta[0][3],xerr_eta[0][3],error_eta[0][3],error_eta[0][3]);      
  
  TLegend *leg = new TLegend(.15,0.72,.45,0.87);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  leg->AddEntry(g0_etavsfracabove55,"H_{T} in [80,100]","p");
  leg->AddEntry(g1_etavsfracabove55,"H_{T} in [100,120]","p");
  leg->AddEntry(g2_etavsfracabove55,"H_{T} in [120,160]","p");

  c1->SetLogy(0); 	
  g0_etavsfracabove55->Draw("AP"); g1_etavsfracabove55->Draw("P"); g2_etavsfracabove55->Draw("P"); //g3_etavsfracabove55->Draw("P"); 
  g0_etavsfracabove55->GetYaxis()->SetRangeUser(0,0.15);  
  leg->Draw("same");
  c1->Print("g_etavsfracabove55_1.pdf");  c1->SaveSource("g_etavsfracabove55_1.C");
  c1->SetLogy(); 	
  g0_etavsfracabove55->Draw("AP"); //g1_etavsfracabove55->Draw("P"); g2_etavsfracabove55->Draw("P"); //g3_etavsfracabove55->Draw("P");
  //g0_etavsfracabove55->GetYaxis()->SetRangeUser(0,0.01);  //leg->Draw("same");
  c1->Print("g_etavsfracabove55_2.pdf");  c1->SaveSource("g_etavsfracabove55_2.C");
  
  */



  
  double num_above_smear_small[2][7]={{0,0,0,0,0,0,0},{0,0,0,0,0,0,0}};
  double num_total_smear_small[2][7]={{0,0,0,0,0,0,0},{0,0,0,0,0,0,0}};
  double ratio_smear_small[7][7], yup_error_smear_small[7][7], ydown_error_smear_small[7][7];
  double x_smear_small[7][7], xerr_smear_small[7][7]; //[smearing][ht bin]
  double x_sum_smear_small[2][7]={{0,0,0,0,0,0,0},{0,0,0,0,0,0,0}};
  double x2_sum_smear_small[2][7]={{0,0,0,0,0,0,0},{0,0,0,0,0,0,0}};
  /*
  double num_above_smear_small[2][7]={{0,0,0,0,0,0,0},{0,0,0,0,0,0,0}};
  double num_total_smear_small[2][7]={{0,0,0,0,0,0,0},{0,0,0,0,0,0,0}};
  double ratio_smear_small[7][11], yup_error_smear_small[7][11], ydown_error_smear_small[7][11];
  double x_smear_small[7][11], xerr_smear_small[7][11]; //[smearing][ht bin]
  double x_sum_smear_small[2][7]={{0,0,0,0,0,0,0},{0,0,0,0,0,0,0}};
  double x2_sum_smear_small[2][7]={{0,0,0,0,0,0,0},{0,0,0,0,0,0,0}};
  */
  int ndiv;
  ndiv= 8; //15;
  if(sample==1 || sample==3) ndiv=4;
  if(sample==0) ndiv=4;
  if(sample==2 || sample==4) ndiv=2;
  if(sample==6 || sample==7) ndiv=12;
  if(sample==-1) ndiv= 10; //5;
  if(sample==-2) ndiv= 10; //6;
  if(sample==-3) ndiv= 20;
  if(sample==-4) ndiv=4;
  if(sample==-5) ndiv=max_htbins-1;
  double ratio_smear_small2[7][ndiv], yup_error_smear_small2[7][ndiv], ydown_error_smear_small2[7][ndiv];
  double x_smear_small2[7][ndiv], xerr_smear_small2[7][ndiv];
  int first_nonzero_bin, last_nonzero_bin, last_nonzero_bin_num; 

  //compute ratios and errors
  for(int j=0; j<1 /*7*/; j++) { //loop through different smearings 
    first_nonzero_bin = -10; 
    for(int i=0; i<max_htbins; i++) {
      if(num_above_smear[j][i]>0 && first_nonzero_bin < 0 && i<ndiv) first_nonzero_bin = x_sum_smear[j][i]/num_total_smear[j][i];
      if(num_above_smear[j][i]>0 && i<ndiv) { last_nonzero_bin = x_sum_smear[j][i]/num_total_smear[j][i]; last_nonzero_bin_num = i+1; }
      if(i<ndiv) { ratio_smear_small2[j][i]=0; yup_error_smear_small2[j][i]=1; ydown_error_smear_small2[j][i]=0;}
      if(i<6) { ratio_smear_small[j][i]=0; yup_error_smear_small[j][i]=1; ydown_error_smear_small[j][i]=0;}
      
      if(num_total_smear[j][i]<=0){ 
	ratio_smear_small2[j][i]=0;
	xerr_smear_small2[j][i] = 0;
	//yup_error_smear_small2[j][i]=0;
	ydown_error_smear_small2[j][i]=0.0;
	x_smear_small2[j][i] = offset + ((double) i + 0.5)*bin_width; //cout<<"\n"<<offset<<" "<<z_pt_min<<" "<<jet_pt_min;
      }
      
      if(num_total_smear[j][i]>0){ 
	ratio_smear[j][i]= num_above_smear[j][i]/num_total_smear[j][i];
	//if(num_above_smear[j][i]<0.5)	num_above_smear[j][i]=1;
	//error_smear[j][i]= num_above_smear[j][i]/num_total_smear[j][i] * sqrt(1./num_above_smear[j][i] + 1./num_total_smear[j][i]);
	error_smear[j][i]= sqrt(ratio_smear[j][i]*(1-ratio_smear[j][i])/num_total_smear_noscale[j][i]); //binomial errors
	//if(num_above_smear[j][i]<=0) {  //If no events in the bin fail, compute error as if 1 event failed
	//  float temp_ratio = 1.0/num_total_smear[j][i];
	//  error_smear[j][i]= sqrt(temp_ratio*(1-temp_ratio)/num_total_smear_noscale[j][i]);
	//}
	if(i<5){
	  xerr_smear_small[j][i] = sqrt(x2_sum_smear[j][i]/num_total_smear[j][i]-pow(x_sum_smear[j][i]/num_total_smear[j][i],2));
	  x_smear_small[j][i] = x_sum_smear[j][i]/num_total_smear[j][i];
	  ratio_smear_small[j][i] = ratio_smear[j][i];
	  yup_error_smear_small[j][i] = error_smear[j][i];
	  if(ratio_smear[j][i]-error_smear[j][i]<0) ydown_error_smear_small[j][i] = ratio_smear[j][i];
	  else ydown_error_smear_small[j][i]= error_smear[j][i];
	}
	if(i>=5 && i<=8) {
	  if(num_above_smear[j][i]>0.51) num_above_smear_small[0][j] += num_above_smear[j][i];
	  num_total_smear_small[0][j] += num_total_smear[j][i];
	  x_sum_smear_small[0][j] += x_sum_smear[j][i];
	  x2_sum_smear_small[0][j] += x2_sum_smear[j][i];
	}
	if(i>=9 && i<=13) {
	  if(num_above_smear[j][i]>0.51) num_above_smear_small[1][j] += num_above_smear[j][i];
	  num_total_smear_small[1][j] += num_total_smear[j][i];
	  x_sum_smear_small[1][j] += x_sum_smear[j][i];
	  x2_sum_smear_small[1][j] += x2_sum_smear[j][i];
	}
	/*
	if(i<9){
	  xerr_smear_small[j][i] = sqrt(x2_sum_smear[j][i]/num_total_smear[j][i]-pow(x_sum_smear[j][i]/num_total_smear[j][i],2));
	  x_smear_small[j][i] = x_sum_smear[j][i]/num_total_smear[j][i];
	  ratio_smear_small[j][i] = ratio_smear[j][i];
	  yup_error_smear_small[j][i] = error_smear[j][i];
	  if(ratio_smear[j][i]-error_smear[j][i]<0) ydown_error_smear_small[j][i] = ratio_smear[j][i];
	  else ydown_error_smear_small[j][i]= error_smear[j][i];
	}
	if(i>=9 && i<=12) {
	  if(num_above_smear[j][i]>0.51) num_above_smear_small[0][j] += num_above_smear[j][i];
	  num_total_smear_small[0][j] += num_total_smear[j][i];
	  x_sum_smear_small[0][j] += x_sum_smear[j][i];
	  x2_sum_smear_small[0][j] += x2_sum_smear[j][i];
	}
	if(i>=13 && i<=17) {
	  if(num_above_smear[j][i]>0.51) num_above_smear_small[1][j] += num_above_smear[j][i];
	  num_total_smear_small[1][j] += num_total_smear[j][i];
	  x_sum_smear_small[1][j] += x_sum_smear[j][i];
	  x2_sum_smear_small[1][j] += x2_sum_smear[j][i];
	}	
	*/
	if(i<ndiv){
	  //cout<<"\n"<<x2_sum_smear[j][i]/num_total_smear[j][i]-pow(x_sum_smear[j][i]/num_total_smear[j][i],2)<<" "<<x2_sum_smear[j][i]<<" "<<x_sum_smear[j][i];
	  //cout<<" "<<num_above_smear[j][i]<<" "<<num_total_smear[j][i]<<" "<<ratio_smear[j][i];
	  xerr_smear_small2[j][i] = sqrt(x2_sum_smear[j][i]/num_total_smear[j][i]-pow(x_sum_smear[j][i]/num_total_smear[j][i],2));
	  x_smear_small2[j][i] = x_sum_smear[j][i]/num_total_smear[j][i];
	  ratio_smear_small2[j][i] = ratio_smear[j][i];
	  yup_error_smear_small2[j][i] = error_smear[j][i];
	  if(num_above_smear[j][i]<=0) yup_error_smear_small2[j][i] = (1.0/num_total_smear[j][i])/(1+1.0/num_total_smear[j][i]);
	  if(ratio_smear[j][i]-error_smear[j][i]<0) ydown_error_smear_small2[j][i] = ratio_smear[j][i];
	  else ydown_error_smear_small2[j][i]= error_smear[j][i];
	}
      }
    }


    /*
    x_smear_small[j][9] = x_sum_smear_small[0][j]/num_total_smear_small[0][j];
    xerr_smear_small[j][9] = sqrt(x2_sum_smear_small[0][j]/num_total_smear_small[0][j]-pow(x_sum_smear_small[0][j]/num_total_smear_small[0][j],2));
    if(num_above_smear_small[0][j]<0.5) num_above_smear_small[0][j]=0.001; //0.5;
    if(num_total_smear_small[0][j]>0){
      ratio_smear_small[j][9]= num_above_smear_small[0][j]/num_total_smear_small[0][j];
      //yup_error_smear_small[j][9]= num_above_smear_small[0][j]/num_total_smear_small[0][j] * sqrt(1./num_above_smear_small[0][j] + 1./num_total_smear_small[0][j]);
      yup_error_smear_small[j][9]= sqrt(ratio_smear_small[j][9]*(1-ratio_smear_small[j][9])/num_total_smear_small[0][j]);
      if(ratio_smear_small[j][9]-yup_error_smear_small[j][9]<0) ydown_error_smear_small[j][9] = ratio_smear_small[j][9];
      else ydown_error_smear_small[j][9] = yup_error_smear_small[j][9];
    }

    //if(j==0) cout<<"\n\n"<<ratio_smear_small[j][9]<<" "<<yup_error_smear_small[j][9]<<" "<<ydown_error_smear_small[j][9]<<"\n\n";
    //exit(1);

    x_smear_small[j][10] = x_sum_smear_small[1][j]/num_total_smear_small[1][j];
    xerr_smear_small[j][10] = sqrt(x2_sum_smear_small[1][j]/num_total_smear_small[1][j]-pow(x_sum_smear_small[1][j]/num_total_smear_small[1][j],2));
    if(num_above_smear_small[1][j]<0.5) num_above_smear_small[1][j]=0.001; //0.5;
    if(num_total_smear_small[1][j]>0){
      ratio_smear_small[j][10]= num_above_smear_small[1][j]/num_total_smear_small[1][j];
      //yup_error_smear_small[j][10]= num_above_smear_small[1][j]/num_total_smear_small[1][j] * sqrt(1./num_above_smear_small[1][j] + 1./num_total_smear_small[1][j]);
      yup_error_smear_small[j][10]= sqrt(ratio_smear_small[j][10]*(1-ratio_smear_small[j][10])/num_total_smear_small[1][j]);
      if(ratio_smear_small[j][10]-yup_error_smear_small[j][10]<0) ydown_error_smear_small[j][10] = ratio_smear_small[j][10];
      else ydown_error_smear_small[j][10] = yup_error_smear_small[j][10];
    }
    */
    
    x_smear_small[j][5] = x_sum_smear_small[0][j]/num_total_smear_small[0][j];
    xerr_smear_small[j][5] = sqrt(x2_sum_smear_small[0][j]/num_total_smear_small[0][j]-pow(x_sum_smear_small[0][j]/num_total_smear_small[0][j],2));
    if(num_above_smear_small[0][j]<0.5) num_above_smear_small[0][j]= 0.001; //0.5;
    if(num_total_smear_small[0][j]>0){
      ratio_smear_small[j][5]= num_above_smear_small[0][j]/num_total_smear_small[0][j];
      //yup_error_smear_small[j][5]= ratio_smear_small[j][5] * sqrt(1./num_above_smear_small[0][j] + 1./num_total_smear_small[0][j]);
      yup_error_smear_small[j][5]= sqrt(ratio_smear_small[j][5]*(1-ratio_smear_small[j][5])/num_total_smear_small[0][j]);
      if(ratio_smear_small[j][5]-yup_error_smear_small[j][5]<0) ydown_error_smear_small[j][5] = ratio_smear_small[j][5];
      else ydown_error_smear_small[j][5] = yup_error_smear_small[j][5];
    }

    x_smear_small[j][6] = x_sum_smear_small[1][j]/num_total_smear_small[1][j];
    xerr_smear_small[j][6] = sqrt(x2_sum_smear_small[1][j]/num_total_smear_small[1][j]-pow(x_sum_smear_small[1][j]/num_total_smear_small[1][j],2));
    if(num_above_smear_small[1][j]<0.5) num_above_smear_small[1][j]= 0.001; //0.5;
    if(num_total_smear_small[1][j]>0){
      ratio_smear_small[j][6]= num_above_smear_small[1][j]/num_total_smear_small[1][j];
      //yup_error_smear_small[j][6]= ratio_smear_small[j][6] * sqrt(1./num_above_smear_small[1][j] + 1./num_total_smear_small[1][j]);
      yup_error_smear_small[j][6]= sqrt(ratio_smear_small[j][6]*(1-ratio_smear_small[j][6])/num_total_smear_small[1][j]);
      if(ratio_smear_small[j][6]-yup_error_smear_small[j][6]<0) ydown_error_smear_small[j][6] = ratio_smear_small[j][6];
      else ydown_error_smear_small[j][6] = yup_error_smear_small[j][6];
    }
    
    
    TGraphAsymmErrors *g_htvsfracabove55 = new TGraphAsymmErrors(6,x_smear_small[j],ratio_smear_small[j],xerr_smear_small[j],xerr_smear_small[j],ydown_error_smear_small[j],yup_error_smear_small[j]);
    //TGraphAsymmErrors *g_htvsfracabove55 = new TGraphAsymmErrors(8,x_smear_small[j],ratio_smear_small[j],xerr_smear_small[j],xerr_smear_small[j],ydown_error_smear_small[j],yup_error_smear_small[j]);
    char name1[64], name2[64], name3[64], name4[64], name5[64], name6[64];
    if(j==0) {
      g_htvsfracabove55->SetTitle("Fraction with #alpha_{T}>0.55 vs HT"); 
      strcpy(name1,"g_htvsfracabove55_aa1.pdf"); strcpy(name2,"g_htvsfracabove55_aa2.pdf"); strcpy(name3,"g_htvsfracabove55_ab1.pdf"); strcpy(name4,"g_htvsfracabove55_ab2.pdf"); 
      strcpy(name5,"g_htvsfracabove55_ac1.pdf"); strcpy(name6,"g_htvsfracabove55_ac2.pdf");
    }
    g_htvsfracabove55->SetMarkerColor(4);
    g_htvsfracabove55->SetMarkerStyle(21);
    g_htvsfracabove55->Write();
    gStyle->SetOptFit();
    TF1* func1 = new TF1("func1","x*[1]+[0]",0,150); //func1->SetParameters(0.3,-0.0005);
    TF1* func2 = new TF1("func2","expo",0,250); // func2->SetParameters(0.3,-0.0005);
    TF1 *func3 = new TF1("func3", "(x<60)*(x*[1]+[0]) + (x>60)*[2]", 0, 150); //func3->SetParameters(0.3,-0.0005,0.01);
    c1->SetLogy(0); 
    //g_htvsfracabove55->Draw("ALP"); g_htvsfracabove55->Fit(func1); c1->Print(name1);
    //g_htvsfracabove55->Draw("AP"); if(test==1) g_htvsfracabove55->Fit(func2,"","",60,240); else g_htvsfracabove55->Fit(func2); if(func2->GetParameter(0)>50){c1->Clear(); c1->Update(); c1->Print(name3);}  else c1->Print(name3);
    //g_htvsfracabove55->Draw("ALP"); g_htvsfracabove55->Fit(func3); c1->Print(name5);
    c1->SetLogy(); 
    //g_htvsfracabove55->Draw("ALP"); g_htvsfracabove55->Fit(func1); c1->Print(name2);
    //g_htvsfracabove55->Draw("AP"); /*g_htvsfracabove55->GetYaxis()->SetRangeUser(0.00001,0.12);*/  if(test==1) g_htvsfracabove55->Fit(func2,"","",40,200); else g_htvsfracabove55->Fit(func2); if(func2->GetParameter(0)>50){c1->Clear(); c1->Update(); c1->Print(name4);} else c1->Print(name4);
    //g_htvsfracabove55->Draw("ALP"); g_htvsfracabove55->Fit(func3); c1->Print(name6);
    gStyle->SetOptFit(0);

    if(1) {
      TGraphAsymmErrors *g2_htvsfracabove55 = new TGraphAsymmErrors(last_nonzero_bin_num,x_smear_small2[j],ratio_smear_small2[j],xerr_smear_small2[j],xerr_smear_small2[j],ydown_error_smear_small2[j],yup_error_smear_small2[j]);
      char name1[64], name2[64], name3[64], name4[64], name5[64], name6[64];
      if(j==0) {
	g2_htvsfracabove55->SetTitle("Fraction with #alpha_{T}>0.55 vs HT"); 
	strcpy(name1,"g_htvsfracabove55_ba1.pdf"); strcpy(name2,"g_htvsfracabove55_ba2.pdf"); strcpy(name3,"g_htvsfracabove55_bb1.pdf"); strcpy(name4,"g_htvsfracabove55_bb2.pdf");
	strcpy(name5,"g_htvsfracabove55_bc1.pdf"); strcpy(name6,"g_htvsfracabove55_bc2.pdf"); 
      }
   

      g2_htvsfracabove55->SetMarkerColor(4);
      g2_htvsfracabove55->SetMarkerStyle(21);
      g2_htvsfracabove55->Write();
      gStyle->SetOptFit();
      TF1* func1 = new TF1("func1","x*[1]+[0]",0,150); //func1->SetParameters(0.3,-0.0005);
      TF1* func2 = new TF1("func2","expo",first_nonzero_bin,last_nonzero_bin); //func2->SetParameters(0.3,-0.0005);
      TF1 *func3 = new TF1("func3", "(x<60)*(x*[1]+[0]) + (x>60)*[2]", 0, 150); //func3->SetParameters(0.3,-0.0005,0.01);
      c1->SetLogy(0);
      //g2_htvsfracabove55->Draw("ALP"); g2_htvsfracabove55->Fit(func1); c1->Print(name1);
      g2_htvsfracabove55->Draw("AP");   /*else g2_htvsfracabove55->Fit(func2);*/  /*if(func2->GetParameter(0)>80){c1->Clear(); c1->Update(); c1->Print(name3);}  else*/ c1->Print(name3);
      //g2_htvsfracabove55->Draw("ALP"); g2_htvsfracabove55->Fit(func3); c1->Print(name5);
      c1->SetLogy(); 
      //g2_htvsfracabove55->Draw("ALP"); g2_htvsfracabove55->Fit(func1); c1->Print(name2);
      gStyle->SetOptFit(0);
      g2_htvsfracabove55->Draw("AP");   if(1) g2_htvsfracabove55->Fit(func2,"","",first_nonzero_bin,last_nonzero_bin);  g2_htvsfracabove55->GetYaxis()->SetRangeUser(0.0005,1);  c1->Print(name4);
      if(j==0) { c1->SaveSource("g_htvsfracabove55_bb2.C"); if(which_file==1) graph1->cd(); if(which_file==2) graph2->cd();  g2_htvsfracabove55->Write(); }
      //g2_htvsfracabove55->Draw("ALP"); g2_htvsfracabove55->Fit(func3); c1->Print(name6);
      gStyle->SetOptFit(0);
    }
  } //loop through different smearings

  //cout<<"\n\n"<<last_nonzero_bin<<" "<<last_nonzero_bin_num<<"\n\n";
  //if(which_file==1) {graph1->Write(); graph1->Close();} 
  //if(which_file==2) {graph2->Write(); graph2->Close();}
  hist_file->cd();

  //  if(n_jet==-2 || n_jet==1) {graph1->Close(); graph1->Close();}
  // if(n_jet==-5 || n_jet==5) graph2->Close();
  //hist_file->Write();

  //  TCanvas* c1= new TCanvas();  

  //  gStyle->SetOptStat(1111111);
  //  gStyle->SetOptTitle(1);

  //exit(1);


c1->SaveSource("g_htvsfracabove55_bb2.C");


  c1->SetLogy(0); h_alphaT->Draw(); c1->Print("alphat1.pdf"); 
  c1->SetLogy(); h_alphaT->Draw(); c1->Print("alphat2.pdf"); 
  //c1->SaveSource("alphat3.C");
  /*
  c1->SetLogy(0); h_alphaT_ht_60_100->Draw(); c1->Print("alphat_ht_60_1001.pdf"); 
  c1->SetLogy(); h_alphaT_ht_60_100->Draw(); c1->Print("alphat_ht_60_1002.pdf");
  //c1->SaveSource("alphat1.C");
  c1->SetLogy(0); h_alphaT_ht_100_150->Draw(); c1->Print("alphat_ht_100_1501.pdf"); 
  c1->SetLogy(); h_alphaT_ht_100_150->Draw(); c1->Print("alphat_ht_100_1502.pdf");
  //c1->SaveSource("alphat2.C");
  c1->SetLogy(0); h_alphaT_ht_100_up->Draw(); c1->Print("alphat_ht_100_up1.pdf"); 
  c1->SetLogy(); h_alphaT_ht_100_up->Draw(); c1->Print("alphat_ht_100_up2.pdf");
  //c1->SaveSource("alphat_JES1_ht100up.C");
  //  c1->SetLogy(); h_alphaT_ht_60_up->Draw(); c1->SaveSource("alphat1.C");
  //  c1->SetLogy(); h_alphaT_ht_40_80->Draw(); c1->SaveSource("alphat1.C");
  //  c1->SetLogy(); h_alphaT_ht_80_120->Draw(); c1->SaveSource("alphat2.C");
  */
  
  c1->SetLogy(); h_alphaT_ht_20_40->Draw(); c1->Print("alphat_ht_20_402.pdf");  c1->SaveSource("alphat_ht_20_40.C");
  c1->SetLogy(); h_alphaT_ht_40_80->Draw(); c1->Print("alphat_ht_40_802.pdf");  c1->SaveSource("alphat_ht_40_80.C");
  c1->SetLogy(); h_alphaT_ht_40_up->Draw(); c1->Print("alphat_ht_40_up2.pdf");  c1->SaveSource("alphat_ht_40_up.C");
  c1->SetLogy(); h_alphaT_ht_80_up->Draw(); c1->Print("alphat_ht_80_up2.pdf");  c1->SaveSource("alphat_ht_80_up.C");
  c1->SetLogy(); h_alphaT_ht_30_50->Draw(); c1->Print("alphat_ht_30_502.pdf");  c1->SaveSource("alphat_ht_30_50.C");
  c1->SetLogy(); h_alphaT_ht_50_80->Draw(); c1->Print("alphat_ht_50_802.pdf");  c1->SaveSource("alphat_ht_50_80.C");
  c1->SetLogy(); h_alphaT_ht_80_120->Draw(); c1->Print("alphat_ht_80_1202.pdf");  c1->SaveSource("alphat_ht_80_120.C");
  c1->SetLogy(); h_alphaT_ht_120_up->Draw(); c1->Print("alphat_ht_120_up2.pdf");  c1->SaveSource("alphat_ht_120_up.C");
  c1->SetLogy(); h_alphaT_ht_120_160->Draw(); c1->Print("alphat_ht_120_1602.pdf");  c1->SaveSource("alphat_ht_120_160.C");
  c1->SetLogy(); h_alphaT_ht_160_up->Draw(); c1->Print("alphat_ht_160_up2.pdf");  c1->SaveSource("alphat_ht_160_up.C");

  
  /*
  c1->SetLogy(0); h_alphaT_ht_150_200->Draw(); c1->Print("alphat_ht_150_2001.pdf"); 
  c1->SetLogy(); h_alphaT_ht_150_200->Draw(); c1->Print("alphat_ht_150_2002.pdf");
  c1->SetLogy(0); h_alphaT_ht_200_250->Draw(); c1->Print("alphat_ht_200_2501.pdf"); 
  c1->SetLogy(); h_alphaT_ht_200_250->Draw(); c1->Print("alphat_ht_200_2502.pdf");
  c1->SetLogy(0); h_alphaT_ht_250_300->Draw(); c1->Print("alphat_ht_250_3001.pdf"); 
  c1->SetLogy(); h_alphaT_ht_250_300->Draw(); c1->Print("alphat_ht_250_3002.pdf");
  c1->SetLogy(0); h_alphaT_ht_300_500->Draw(); c1->Print("alphat_ht_300_5001.pdf"); 
  c1->SetLogy(); h_alphaT_ht_300_500->Draw(); c1->Print("alphat_ht_300_5002.pdf");
  c1->SetLogy(0); h_alphaT_ht_500up->Draw(); c1->Print("alphat_ht_500up1.pdf"); 
  c1->SetLogy(); h_alphaT_ht_500up->Draw(); c1->Print("alphat_ht_500up2.pdf");
  */  

  c1->SetLogy();  h_jet_pt->Draw(); c1->Print("jet_pt.pdf");
  c1->SetLogy();  h_jet_pt_wide->Draw(); c1->Print("jet_pt_wide.pdf");
  c1->SetLogy(0); h_jet_phi->Draw(); c1->Print("jet_phi.pdf");
  c1->SetLogy();  h_jet_num->Draw(); c1->Print("jet_num.pdf"); 
  c1->SetLogy(0); h_jet_eta->Draw(); c1->Print("jet_eta.pdf"); 
  c1->SetLogy(); h_jet_fRBX->Draw(); c1->Print("jet_fRBX.pdf"); 
  c1->SetLogy(); h_jet_fHPD->Draw(); c1->Print("jet_fHPD.pdf"); 
  c1->SetLogy();  h_jet_n90Hits->Draw(); c1->Print("jet_n90Hits.pdf"); 
  c1->SetLogy(); h_jet_emf->Draw(); c1->Print("jet_emf.pdf"); 
  c1->SetLogy();  h_jet_btag_jetProb->Draw(); c1->Print("jet_btag_jetProb.pdf"); 
  c1->SetLogy();  h_jet_pt_btag->Draw(); c1->Print("jet_pt_btag.pdf");

  c1->SetLogy();  h_pho_pt->Draw(); c1->Print("pho_pt.pdf"); 
  c1->SetLogy();  h_pho_pt->Draw(); c1->Print("pho_pt.pdf");
  c1->SetLogy(0); h_pho_phi->Draw(); c1->Print("pho_phi.pdf");
  c1->SetLogy(0); h_pho_eta->Draw(); c1->Print("pho_eta.pdf");
  c1->SetLogy(0); h_pho_E1overE3x3->Draw(); c1->Print("pho_E1overE3x3.pdf");
  c1->SetLogy(0); h_pho_sigmaIetaIeta->Draw(); c1->Print("pho_sigmaIetaIeta.pdf");
  c1->SetLogy(0); h_pho_hasPixelSeed->Draw(); c1->Print("pho_hasPixelSeed.pdf");
  c1->SetLogy(0); h_pho_hadOverEM->Draw(); c1->Print("pho_hadOverEM.pdf");
  c1->SetLogy(0); h_pho_isoEcalRecHitDR04->Draw(); c1->Print("pho_isoEcalRecHitDR04.pdf");
  c1->SetLogy(0); h_pho_isoHcalRecHitDR04->Draw(); c1->Print("pho_isoHcalRecHitDR04.pdf");
  c1->SetLogy(0); h_pho_isoHollowTrkConeDR04->Draw(); c1->Print("pho_isoHollowTrkConeDR04.pdf");

  c1->SetLogy();  h_ht->Draw(); c1->Print("ht.pdf"); 
  c1->SetLogy();  h_mht->Draw(); c1->Print("mht.pdf"); 
  c1->SetLogy();  h_met->Draw(); c1->Print("met.pdf"); 
  c1->SetLogy();  h_genht->Draw(); c1->Print("ht_generator.pdf"); 
  c1->SetLogy(0);  h_ht_vs_htplusmet->Draw(); c1->Print("ht_vs_htplusmet.pdf"); 
  c1->SetLogy(0);  p_ht_vs_htplusmet->Draw(); c1->Print("ht_vs_htplusmet_profile.pdf"); 
  c1->SetLogy(0);  h_ht_vs_genht->Draw(); c1->Print("ht_vs_genht.pdf"); 
  c1->SetLogy(0);  p_ht_vs_genht->Draw(); c1->Print("ht_vs_genht_profile.pdf"); 
  c1->SetLogy();  h_Npv->Draw(); c1->Print("Npv.pdf"); 
  c1->SetLogy();  h_Npv_good->Draw(); c1->Print("Npv_good.pdf"); 

  c1->SetLogy();  h_measuredht->Draw(); c1->Print("ht_measured.pdf"); 
  c1->SetLogy();  h_measuredht_failalphat->Draw(); c1->Print("ht_measured_failalphat.pdf"); 
  /*
  c1->SetLogy(0);  h_phopt_vs_jetpt_overlap->Draw(); c1->Print("phopt_vs_jetpt_overlap.pdf"); 
  c1->SetLogy(0); h_pho_isoEcalRecHitDR04_overlap->Draw(); c1->Print("pho_isoEcalRecHitDR04_overlap.pdf"); 
  c1->SetLogy(0);  h_jet_emf_overlap->Draw(); c1->Print("jet_emf_overlap.pdf"); 
  c1->SetLogy(0);  h_phojet_balance->Draw(); c1->Print("phojet_balance.pdf"); c1->SaveSource("phojet_balance.C");
  c1->SetLogy(0);  h_dijet_balance->Draw(); c1->Print("dijet_balance.pdf"); 
  c1->SetLogy(0);  h_dphi_phomet->Draw(); c1->Print("dphi_phomet.pdf");
  c1->SetLogy(0);  h_dphi_phomet_met0to20->Draw(); c1->Print("dphi_phomet_met0to20.pdf"); 
  c1->SetLogy(0);  h_dphi_phomet_met20toinf->Draw(); c1->Print("dphi_phomet_met20toinf.pdf");
  */

  hist_file->Close();
  input->Close();
  
 
  //system("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=mu_plots/data_tighter_mu_alphaT.pdf *.pdf"); //run from command line
  //system("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=pho1jet.pdf pho_*.pdf g_htvsfracabove55_bb2.pdf"); 

}




/////////////////////////////////////////////////////////////////////////////
float smear_jet_et(float et, int type) {
  
  float smear_et;
  smear_et=0;
  
  float rms, dummy;

  if(type==0){ //no change
    smear_et = et; 
  }  

  if(type==1){ //10% gaussian
    rms=et*0.1; 
    smear_et = et + gRandom->Gaus(0,rms);
  }
    
  if(type==2){ //30% gaussian
    rms=et*0.3; 
    smear_et = et + gRandom->Gaus(0,rms);
  }
    
  if(type==3){ //10% gaussian Left, 30% gaussian Right
    rms=et*0.1; 
    dummy = gRandom->Gaus(0,rms);
    if(dummy<=0) smear_et = et + dummy;
    if(dummy>0) smear_et = et + fabs(gRandom->Gaus(0,3.0*rms));
  }
  
  if(type==4){ //30% gaussian Left, 10% gaussian Right
    rms=et*0.3; 
    dummy = gRandom->Gaus(0,rms);
    if(dummy<=0) smear_et = et + dummy;
    if(dummy>0) smear_et = et + fabs(gRandom->Gaus(0,rms/3.0));
  }  
  
  if(type==5){ //Remove 10% of jets
    dummy = gRandom->Uniform();
    if(dummy<=0.1) smear_et = 0;
    if(dummy>0.1) smear_et = et;
  } 

  if(type==6){ //50% gaussian Left, 10% gaussian Right
    rms=et*0.5; 
    dummy = gRandom->Gaus(0,rms);
    if(dummy<=0) smear_et = et + dummy;
    if(dummy>0) smear_et = et + fabs(gRandom->Gaus(0,rms/5.0));
  } 

  if(type==7){ //50% gaussian
    rms=et*0.5; 
    smear_et = et + gRandom->Gaus(0,rms);
  }

  if(type==8){ //50% gaussian Left
    rms=et*0.5; 
    dummy = fabs(gRandom->Gaus(0,rms));
    smear_et = et - dummy;
  }


  if(type==10){ //5% gaussian 80% of the time, 30% gaussian Left 20% of the time
    if(gRandom->Uniform(0,1)<0.1) {
      rms=et*0.05; 
      dummy = gRandom->Gaus(0,rms);
      smear_et = et + dummy;
    }
    else {
      rms=et*0.80; 
      dummy = fabs(gRandom->Gaus(0,rms));
      smear_et = et + dummy;
    }
  } 


  /*
  //change jet pt spectrum
  double x = smear_et;
  if(smear_et<90 && gRandom->Uniform() > exp(-pow((x-90),2)/500)) smear_et=0;
  */ 

  return smear_et;
}



/////////////////////////////////////////////////////////////////////////////
float smear_phi(float phi, int type) {
  
  float smear_phi;
  
  float rms, dummy;

  if(type==0){ //no change
    smear_phi = phi; 
  }  

  if(type==1){ //gaussian of width 0.4 rad
    rms = 0.4; 
    smear_phi = phi + gRandom->Gaus(0,rms);
  }
   
  if(type==2){ //gaussian of width 1.0 rad
    rms = 1.0; 
    smear_phi = phi + gRandom->Gaus(0,rms);
  }
   
  if(type==3){ //uniformly distributed in phi
    smear_phi = gRandom->Uniform(0,2*pi);
  }
    
  if(type==4){ //opposite direction in phi
    smear_phi = phi + pi;
  }

  return smear_phi;
}



//////////////////////////////////////////////////////////////////////////////
//dump objects from an event
void event_dump(Event_Tree* evt2, int num_removed, float* ht_vec, float ht, TLorentzVector DiMuon)
{
  int n;
  float mht,mht_phi;
  mht=sqrt(ht_vec[0]*ht_vec[0]+ht_vec[1]*ht_vec[1]);
  mht_phi= findphi(-ht_vec[0],-ht_vec[1]);

  cout<<"\nRun, Event = "<<evt2->run<<", "<<evt2->event<<"\n";
  cout<<"MET "<<evt2->met<<", MET Phi "<<evt2->met_phi<<", HT "<<ht<<", MHT "<<mht<<", MHT Phi "<<mht_phi<<"\n";

  cout<<"Z pt "<<DiMuon.Pt()<<", Z Phi "<<DiMuon.Phi()<<", Z Eta "<<DiMuon.Eta()<<", Z Mass "<<DiMuon.M()<<"\n";

  n=0;
  cout<<"There are "<<evt2->num_jet-num_removed<<" jets: \n";
  for(int ijet=0; ijet<evt2->num_jet; ijet++) { //loop over jets
    n++;
    //if(evt2->jet_remove[ijet]==1) continue; //skip removed jets
    cout<<"Jet "<<n<<": Px "<<evt2->jet_px[ijet]<<", Py "<<evt2->jet_py[ijet]<<", Pt "<<evt2->jet_et[ijet]<<", eta "<<evt2->jet_eta[ijet]<<", phi "<<evt2->jet_phi[ijet]<<", emf "<<evt2->jet_emf[ijet]<<  ",  removed "<<evt2->jet_remove[ijet]<< " \n";
  }


  cout<<"There are "<<evt2->num_mu<<" muons: \n";
  if(evt2->num_mu>0) 
    cout<<"Muon "<<1<<": Px "<<evt2->mu1_px<<", Py "<<evt2->mu1_py<<", Pt "<<evt2->mu1_pt<<", eta "<<evt2->mu1_eta<<", phi "<<evt2->mu1_phi<<", type "<<evt2->mu1_pass<<" \n";
  if(evt2->num_mu>1) 
    cout<<"Muon "<<2<<": Px "<<evt2->mu2_px<<", Py "<<evt2->mu2_py<<", Pt "<<evt2->mu2_pt<<", eta "<<evt2->mu2_eta<<", phi "<<evt2->mu2_phi<<", type "<<evt2->mu2_pass<<" \n";


  cout<<"There are "<<evt2->num_ele<<" electrons: \n";
  if(evt2->num_ele>0) 
    cout<<"Electron "<<1<<": Px "<<evt2->ele1_px<<", Py "<<evt2->ele1_py<<", Pt "<<evt2->ele1_et<<", eta "<<evt2->ele1_eta<<", phi "<<evt2->ele1_phi<<", type "<<evt2->ele1_pass<<" \n";
  if(evt2->num_ele>1) 
    cout<<"Electron "<<2<<": Px "<<evt2->ele2_px<<", Py "<<evt2->ele2_py<<", Pt "<<evt2->ele2_et<<", eta "<<evt2->ele2_eta<<", phi "<<evt2->ele2_phi<<", type "<<evt2->ele2_pass<<" \n";

  //exit(1);
}
