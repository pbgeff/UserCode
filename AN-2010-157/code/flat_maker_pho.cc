#include "flat_maker_pho.hh" //Version 24 of cfA
#include "inJSON.h" //For good run list
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include "TString.h"
#include "TVector3.h"
#include "TLorentzVector.h"
using namespace std;


double pi = 3.141592653589;

int flat_maker(){

  Int_t debug = 0;

  TString output_filename;

  TChain *chainB = new TChain("event");
  TChain *chainV = new TChain("event");

  //Data
  chainB->Add("/data01/cfA/EG_Run2010A-PromptReco-v1_RECO_UCSB0206_v24/cfA*.root/configurableAnalysis/eventB");
  chainV->Add("/data01/cfA/EG_Run2010A-PromptReco-v1_RECO_UCSB0206_v24/cfA*.root/configurableAnalysis/eventV");
  int MC_flag = 0;
  output_filename="CMS_phojet_data_EG_Run2010A-PromptReco-v1_RECO.root";
/*
  chainB->Add("/data01/cfA/MinimumBias_Commissioning10-SD_EG-v9_RECO_UCSB0214_v24/cfA*.root/configurableAnalysis/eventB");
  chainV->Add("/data01/cfA/MinimumBias_Commissioning10-SD_EG-v9_RECO_UCSB0214_v24/cfA*.root/configurableAnalysis/eventV");
  int MC_flag = 0;
  output_filename="CMS_phojet_data_MinimumBias_Commissioning10-SD_EG-v9_RECO.root";
*/
/*
  chainB->Add("/data01/cfA/EG_Run2010A-PromptReco-v2_RECO_UCSB0215_v24/cfA*.root/configurableAnalysis/eventB");
  chainV->Add("/data01/cfA/EG_Run2010A-PromptReco-v2_RECO_UCSB0215_v24/cfA*.root/configurableAnalysis/eventV");
  int MC_flag = 0;
  output_filename="CMS_phojet_data_EG_Run2010A-PromptReco-v2_RECO.root";
*/


  //MC
  //chainB->Add("/data01/cfA/MinBias_Spring10-START3X_V26A_356ReReco-v1_GEN-SIM-RECO_UCSB0162_v20/cfA*.root/configurableAnalysis/eventB");
  //chainV->Add("/data01/cfA/MinBias_Spring10-START3X_V26A_356ReReco-v1_GEN-SIM-RECO_UCSB0162_v20/cfA*.root/configurableAnalysis/eventV");
  //int MC_flag = 1;
  //output_filename="CMS_phojet_MC.root";


  // Initializes both Trees, i.e. set all the branch addresses, etc.
  InitializeB(chainB);
  InitializeV(chainV);


  TFile *outFile = new TFile(output_filename.Data(),"RECREATE");

  outFile->cd();



  //declare ntuple
  TTree* tree = new TTree("Event_Tree","Event Tree");    
  tree->Branch("run",&_run,"run/I");
  tree->Branch("event",&_event,"event/I");
  tree->Branch("lumiblock",&_lumiblock,"lumiblock/I");	  
  tree->Branch("met",&_met,"met/F");
  tree->Branch("metx",&_metx,"metx/F");
  tree->Branch("mety",&_mety,"mety/F");
  tree->Branch("met_phi",&_metphi,"met_phi/F");
  tree->Branch("rawmet",&_rawmet,"rawmet/F");
  tree->Branch("rawmet_phi",&_rawmetphi,"rawmet_phi/F");
  tree->Branch("vtxmet",&_vtxmet,"vtxmet/F");
  tree->Branch("vtxmet_phi",&_vtxmetphi,"vtxmet_phi/F");
  tree->Branch("num_vert",&_num_vert,"num_vert/I");
  tree->Branch("num_jet",&_num_jet,"num_jet/I");
  tree->Branch("jet_et",_jet_et,"jet_et[num_jet]/F");
  tree->Branch("jet_px",_jet_px,"jet_px[num_jet]/F");
  tree->Branch("jet_py",_jet_py,"jet_py[num_jet]/F");
  tree->Branch("jet_pz",_jet_pz,"jet_pz[num_jet]/F");
  tree->Branch("jet_eta",_jet_eta,"jet_eta[num_jet]/F");
  tree->Branch("jet_E",_jet_E,"jet_E[num_jet]/F");
  tree->Branch("jet_theta",_jet_theta,"jet_theta[num_jet]/F");
  tree->Branch("jet_phi",_jet_phi,"jet_phi[num_jet]/F");
  tree->Branch("jet_emf",_jet_emf,"jet_emf[num_jet]/F");
  tree->Branch("jet_fHPD",_jet_fHPD,"jet_fHPD[num_jet]/F");
  tree->Branch("jet_fRBX",_jet_fRBX,"jet_fRBX[num_jet]/F");
  tree->Branch("jet_n90Hits",_jet_n90Hits,"jet_n90Hits[num_jet]/I");
  tree->Branch("jet_btag_jetProb",_jet_btag_jetProb,"jet_btag_jetProb[num_jet]/F");
  tree->Branch("jet_scalefac",_jet_scalefac,"jet_scalefac[num_jet]/F");
  tree->Branch("jet_remove",_jet_remove,"jet_remove[num_jet]/I");
  tree->Branch("ele1_et",&_ele1_et,"ele1_et/F");
  tree->Branch("ele1_px",&_ele1_px,"ele1_px/F");
  tree->Branch("ele1_py",&_ele1_py,"ele1_py/F");
  tree->Branch("ele1_pz",&_ele1_pz,"ele1_pz/F");
  tree->Branch("ele1_E",&_ele1_E,"ele1_E/F");
  tree->Branch("ele1_eta",&_ele1_eta,"ele1_eta/F");
  tree->Branch("ele1_phi",&_ele1_phi,"ele1_phi/F");
  tree->Branch("ele1_pass",&_ele1_pass,"ele1_pass/I");
  tree->Branch("ele2_et",&_ele2_et,"ele2_et/F");
  tree->Branch("ele2_px",&_ele2_px,"ele2_px/F");
  tree->Branch("ele2_py",&_ele2_py,"ele2_py/F");
  tree->Branch("ele2_pz",&_ele2_pz,"ele2_pz/F");
  tree->Branch("ele2_E",&_ele2_E,"ele2_E/F");
  tree->Branch("ele2_eta",&_ele2_eta,"ele2_eta/F");
  tree->Branch("ele2_phi",&_ele2_phi,"ele2_phi/F");
  tree->Branch("ele2_pass",&_ele2_pass,"ele2_pass/I");
  tree->Branch("num_ele",&_num_ele,"num_ele/I");
  tree->Branch("mu1_et",&_mu1_et,"mu1_et/F");
  tree->Branch("mu1_px",&_mu1_px,"mu1_px/F");
  tree->Branch("mu1_py",&_mu1_py,"mu1_py/F");
  tree->Branch("mu1_pz",&_mu1_pz,"mu1_pz/F");
  tree->Branch("mu1_pt",&_mu1_pt,"mu1_pt/F");
  tree->Branch("mu1_E",&_mu1_E,"mu1_E/F");
  tree->Branch("mu1_eta",&_mu1_eta,"mu1_eta/F");
  tree->Branch("mu1_phi",&_mu1_phi,"mu1_phi/F");
  tree->Branch("mu1_emE",&_mu1_emE,"mu1_emE/F");
  tree->Branch("mu1_hadE",&_mu1_hadE,"mu1_hadE/F");
  /*
  tree->Branch("mu1_trkZatCES",&_mu1_trkZatCES,"mu1_trkZatCES/F");
  tree->Branch("mu1_D0C",&_mu1_D0C,"mu1_D0C/F");
  tree->Branch("mu1_isol",&_mu1_isol,"mu1_isol/F");
  tree->Branch("mu1_CMUStub",&_mu1_CMUStub,"mu1_CMUStub/F");
  tree->Branch("mu1_CMPStub",&_mu1_CMPStub,"mu1_CMPStub/F");
  tree->Branch("mu1_CMXStub",&_mu1_CMXStub,"mu1_CMXStub/F");
  tree->Branch("mu1_BMUStub",&_mu1_BMUStub,"mu1_BMUStub/F");
  tree->Branch("mu1_CmuDx",&_mu1_CmuDx,"mu1_CmuDx/F");
  tree->Branch("mu1_CmpDx",&_mu1_CmpDx,"mu1_CmpDx/F");
  tree->Branch("mu1_CmxDx",&_mu1_CmxDx,"mu1_CmxDx/F");
  tree->Branch("mu1_NumAxSeg",&_mu1_NumAxSeg,"mu1_NumAxSeg/F");
  tree->Branch("mu1_NumStSeg",&_mu1_NumStSeg,"mu1_NumStSeg/F");
  tree->Branch("mu1_trkchisqproblog10",&_mu1_trkchisqproblog10,"mu1_trkchisqproblog10/F");
  tree->Branch("mu1_trkZ0",&_mu1_trkZ0,"mu1_trkZ0/F");
  tree->Branch("mu1_SIHits",&_mu1_SIHits,"mu1_SIHits/F");
  tree->Branch("mu1_trkisol",&_mu1_trkisol,"mu1_trkisol/F");
  */
  tree->Branch("mu1_pass",&_mu1_pass,"mu1_pass/I");
  tree->Branch("mu2_et",&_mu2_et,"mu2_et/F");
  tree->Branch("mu2_px",&_mu2_px,"mu2_px/F");
  tree->Branch("mu2_py",&_mu2_py,"mu2_py/F");
  tree->Branch("mu2_pz",&_mu2_pz,"mu2_pz/F");
  tree->Branch("mu2_pt",&_mu2_pt,"mu2_pt/F");
  tree->Branch("mu2_E",&_mu2_E,"mu2_E/F");
  tree->Branch("mu2_eta",&_mu2_eta,"mu2_eta/F");
  tree->Branch("mu2_phi",&_mu2_phi,"mu2_phi/F");
  tree->Branch("mu2_emE",&_mu2_emE,"mu2_emE/F");
  tree->Branch("mu2_hadE",&_mu2_hadE,"mu2_hadE/F");
  /*
  tree->Branch("mu2_trkZatCES",&_mu2_trkZatCES,"mu2_trkZatCES/F");
  tree->Branch("mu2_D0C",&_mu2_D0C,"mu2_D0C/F");
  tree->Branch("mu2_isol",&_mu2_isol,"mu2_isol/F");
  tree->Branch("mu2_CMUStub",&_mu2_CMUStub,"mu2_CMUStub/F");
  tree->Branch("mu2_CMPStub",&_mu2_CMPStub,"mu2_CMPStub/F");
  tree->Branch("mu2_CMXStub",&_mu2_CMXStub,"mu2_CMXStub/F");
  tree->Branch("mu2_BMUStub",&_mu2_BMUStub,"mu2_BMUStub/F");
  tree->Branch("mu2_CmuDx",&_mu2_CmuDx,"mu2_CmuDx/F");
  tree->Branch("mu2_CmpDx",&_mu2_CmpDx,"mu2_CmpDx/F");
  tree->Branch("mu2_CmxDx",&_mu2_CmxDx,"mu2_CmxDx/F");
  tree->Branch("mu2_NumAxSeg",&_mu2_NumAxSeg,"mu2_NumAxSeg/F");
  tree->Branch("mu2_NumStSeg",&_mu2_NumStSeg,"mu2_NumStSeg/F");
  tree->Branch("mu2_trkchisqproblog10",&_mu2_trkchisqproblog10,"mu2_trkchisqproblog10/F");
  tree->Branch("mu2_trkZ0",&_mu2_trkZ0,"mu2_trkZ0/F");
  tree->Branch("mu2_SIHits",&_mu2_SIHits,"mu2_SIHits/F");
  tree->Branch("mu2_trkisol",&_mu2_trkisol,"mu2_trkisol/F");
  */
  tree->Branch("mu2_pass",&_mu2_pass,"mu2_pass/I");
  tree->Branch("num_mu",&_num_mu,"num_mu/I");


  //declare photon ntuple
  TTree* pho_tree = new TTree("Photon_Tree","Photon Tree");    
  pho_tree->Branch("run",&_run,"run/I");
  pho_tree->Branch("event",&_event,"event/I");	  
  pho_tree->Branch("num_pho",&_num_pho,"num_pho/I");
  pho_tree->Branch("pho_et",_pho_et,"pho_et[num_pho]/F");
  pho_tree->Branch("pho_px",_pho_px,"pho_px[num_pho]/F");
  pho_tree->Branch("pho_py",_pho_py,"pho_py[num_pho]/F");
  pho_tree->Branch("pho_pz",_pho_pz,"pho_pz[num_pho]/F");
  pho_tree->Branch("pho_eta",_pho_eta,"pho_eta[num_pho]/F");
  pho_tree->Branch("pho_E",_pho_E,"pho_E[num_pho]/F");
  pho_tree->Branch("pho_theta",_pho_theta,"pho_theta[num_pho]/F");
  pho_tree->Branch("pho_phi",_pho_phi,"pho_phi[num_pho]/F");
  pho_tree->Branch("pho_hasPixelSeed",_pho_hasPixelSeed,"pho_hasPixelSeed[num_pho]/I");
  pho_tree->Branch("pho_hadOverEM",_pho_hadOverEM,"pho_hadOverEM[num_pho]/F");
  pho_tree->Branch("pho_isoEcalRecHitDR04",_pho_isoEcalRecHitDR04,"pho_isoEcalRecHitDR04[num_pho]/F");
  pho_tree->Branch("pho_isoHcalRecHitDR04",_pho_isoHcalRecHitDR04,"pho_isoHcalRecHitDR04[num_pho]/F");
  pho_tree->Branch("pho_isoHollowTrkConeDR04",_pho_isoHollowTrkConeDR04,"pho_isoHollowTrkConeDR04[num_pho]/F");
  pho_tree->Branch("pho_maxEnergyXtal",_pho_maxEnergyXtal,"pho_maxEnergyXtal[num_pho]/F");
  pho_tree->Branch("pho_e3x3",_pho_e3x3,"pho_e3x3[num_pho]/F");

  pho_tree->Branch("pho_sigmaIetaIeta",_pho_sigmaIetaIeta,"pho_sigmaIetaIeta[num_pho]/F");
  pho_tree->Branch("pho_scPhiWidth",_pho_scPhiWidth,"pho_scPhiWidth[num_pho]/F");
  pho_tree->Branch("pho_scEtaWidth",_pho_scEtaWidth,"pho_scEtaWidth[num_pho]/F");
  pho_tree->Branch("pho_r9",_pho_r9,"pho_r9[num_pho]/F");

  pho_tree->Branch("pho_istight",_pho_istight,"pho_istight[num_pho]/I");
  pho_tree->Branch("pho_isloose",_pho_isloose,"pho_isloose[num_pho]/I");





  //Number of events to loop over
  Int_t nentries = (Int_t)chainB->GetEntries();
  cout<<"The number of entries is: "<<nentries<<endl;

  float beamx = 0.0322;
  float beamy = 0.;

  //Main event loop
  for(int ia = 0; ia<nentries; ia++){ 

    if(ia && ia%100000==0) cout<<"\nProcessed entries: "<<ia;


    bool LooseMuon[1000];
    bool TightMuon[1000];

    bool LooseE[1000];
    bool TightE[1000];

    bool LooseJet[10000];
    bool GoodJet[10000];

    bool LoosePho[1000];
    bool TightPho[1000];

    chainB->GetEntry(ia);
    chainV->GetEntry(ia);

    //Physics Declared
    //if(!experimentType && !MC_flag) continue;

    //Trigger Bits
    if(!L1Bit_0 || (!L1Bit_40 && !L1Bit_41)) continue;
    if(!HLT_Photon10_L1R) continue;
    if(L1Bit_36 || L1Bit_37 || L1Bit_38 || L1Bit_39) continue; 

    //Require vertex
    if(Npv==0) continue;
    else if(pv_ndof->at(0)<5 || fabs(pv_z->at(0))>15 ) continue;

    //Reject Monster Events
    int num_highPurity = 0;
    if(Ntracks>10) {
      for (unsigned int k=0; k<Ntracks; k++) {
	if(tracks_highPurity->at(k)>0.5) num_highPurity++;
      }
      if((float)num_highPurity/Ntracks < 0.2) continue;
    }


    //Check good run list
    if(!MC_flag && !inJSON(run,lumiblock)) continue;


    _run = run;
    _event = event;
    _lumiblock = lumiblock;
    _met = mets_AK5_et->at(0);
    _metx = mets_AK5_et->at(0)*cos(mets_AK5_phi->at(0));
    _mety = mets_AK5_et->at(0)*sin(mets_AK5_phi->at(0));
    _metphi = mets_AK5_phi->at(0);
    _rawmet = pfmets_et->at(0);
    _rawmetphi = pfmets_phi->at(0);
    _vtxmet =  tcmets_et->at(0); 
    _vtxmetphi = tcmets_phi->at(0);
    if(_metphi<0) _metphi+=2*pi;
    _num_vert = Npv;


    // Clear old info from previous event
    _jet_et[0] = _jet_px[0] = _jet_py[0] = _jet_pz[0] = _jet_eta[0] = -999;
    _jet_E[0] = _jet_theta[0] = _jet_phi[0] = _jet_emf[0] = -999;
    _num_jet = 0;
    _ele1_et = _ele1_px = _ele1_py = _ele1_pz = _ele1_E = -999;
    _ele2_et = _ele2_px = _ele2_py = _ele2_pz = _ele2_E = -999;
    _ele1_pass = _ele2_pass = 0;
    _num_ele = 0;
    _mu1_et = _mu1_px = _mu1_py = _mu1_pz = _mu1_pt = _mu1_E = -999;
    _mu2_et = _mu2_px = _mu2_py = _mu2_pz = _mu2_pt = _mu2_E = -999;
    _mu1_pass = _mu2_pass = 0;
    _num_mu = 0;
    _pho_et[0] = _pho_px[0] = _pho_py[0] = _pho_pz[0] = _pho_E[0] = -999;
    _pho_eta[0] = _pho_theta[0] = _pho_phi[0] = -999;
    _pho_istight[0] = 0; _pho_isloose[0] = 0;
    _num_pho=0;



    //loop through jets
    for (unsigned int k=0; k<jets_AK5_pt->size(); k++) {
      if(jets_AK5_et->at(k)<10) continue;
      //if(fabs(jets_AK5_eta->at(k))>3.0) continue;
      //write out jet data
      float et = jets_AK5_et->at(k);
      float phi = jets_AK5_phi->at(k);
      float eta = jets_AK5_eta->at(k);
      _jet_theta[_num_jet] = 2*atan(exp(fabs(eta)));
      _jet_et[_num_jet]= et;
      _jet_px[_num_jet]= et*cos(phi);
      _jet_py[_num_jet]= et*sin(phi);
      _jet_pz[_num_jet] = et*(1./tan(_jet_theta[_num_jet]));
      if (eta<0) _jet_pz[_num_jet] *= -1;
      _jet_eta[_num_jet] = eta;
      if(phi<0) phi+=2*pi;
      _jet_phi[_num_jet] = phi;
      _jet_E[_num_jet]=sqrt(_jet_et[_num_jet]*_jet_et[_num_jet]
				+_jet_pz[_num_jet]*_jet_pz[_num_jet]);
      _jet_emf[_num_jet]= jets_AK5_emf->at(k);
      _jet_fHPD[_num_jet]= jets_AK5_fHPD->at(k);
      _jet_fRBX[_num_jet]= jets_AK5_fRBX->at(k);
      _jet_n90Hits[_num_jet]= jets_AK5_n90Hits->at(k);
      _jet_btag_jetProb[_num_jet]= jets_AK5_btag_jetProb->at(k); 
      _jet_scalefac[_num_jet]= 1;
      _jet_remove[_num_jet]= 0;
      _num_jet++;
    }//loop through jets


    
    /* //Don't look at electrons or muons

    //loop through electrons
    for (unsigned int k=0; k<els_pt->size(); k++) {
      // Identify each electron as loose or tight.
      // First calculate variables that I need.
      //Correct d0 for PV not beamspot 
      float d0; 
      if(Npv) d0 = els_d0dum->at(k)-pv_x->at(0)*sin(els_tk_phi->at(k))+pv_y->at(0)*cos(els_tk_phi->at(k));
      else d0 = els_d0dum->at(k);
      // First check loose
      LooseE[k]=true;
      if (els_robustLooseId->at(k) == 0) LooseE[k]=false;
      if (els_pt->at(k) < 10.) LooseE[k]=false;
      if (fabs(els_eta->at(k)) > 2.5) LooseE[k]=false;
      if (fabs(d0)>=0.2) LooseE[k]=false;
      // Then check tight
      TightE[k]=true;
      if (els_robustTightId->at(k) == 0) TightE[k]=false;
      if (els_pt->at(k) < 20.) TightE[k]=false;
      if (fabs(els_eta->at(k)) > 2.5) TightE[k]=false;
      if ((els_ecalIso->at(k)+els_hcalIso->at(k)+els_tIso->at(k))/els_et->at(k)>=0.1) TightE[k] = false;
      if (fabs(d0)>=0.2) TightE[k]=false;
      
      int tightE=0; int looseE=0; int triggerE=0;
      if(TightE[k]) tightE=1; if(LooseE[k]) looseE=1;
      //if(TightE[k] || LooseE[k]) triggerE=1;
      
      if (tightE || looseE || triggerE) { // Electron is good enough to save
	_num_ele++;
	if(_num_ele==1) {
	  _ele1_px=els_px->at(k);
	  _ele1_py=els_py->at(k);
	  _ele1_pz=els_pz->at(k);
	  _ele1_E=els_energy->at(k);
	  _ele1_eta=els_eta->at(k);
	  _ele1_phi=els_phi->at(k);
	  if(_ele1_phi<0) _ele1_phi+=2*pi;
	  _ele1_et=els_et->at(k);
	  if(looseE)_ele1_pass+=2; if(tightE)_ele1_pass+=4; if(triggerE)_ele1_pass+=8; //if(passes_minimal==1)_ele1_pass+=1; 
	}
	if(_num_ele==2) {
	  _ele2_px=els_px->at(k);
	  _ele2_py=els_py->at(k);
	  _ele2_pz=els_pz->at(k);
	  _ele2_E=els_energy->at(k);
	  _ele2_eta=els_eta->at(k);
	  _ele2_phi=els_phi->at(k);
	  if(_ele2_phi<0) _ele2_phi+=2*pi;
	  _ele2_et=els_et->at(k);
	  if(looseE)_ele2_pass+=2; if(tightE)_ele2_pass+=4; if(triggerE)_ele2_pass+=8; //if(passes_minimal==1)_ele2_pass+=1; 
	}
      } // Electron is good enough to save
    } //loop through electrons


    



    //loop through muons
    for (unsigned int k=0; k<mus_pt->size(); k++) {
      //check if muons are loose or tight 
      // First calculate variables that I need.
      //Correct d0 for PV not beamspot
      float d0; 
      if(Npv) d0 = mus_tk_d0dum->at(k)-pv_x->at(0)*sin(mus_tk_phi->at(k))+pv_y->at(0)*cos(mus_tk_phi->at(k));
      else d0 = mus_tk_d0dum->at(k);
      // First check loose
      LooseMuon[k]=true;
      if (mus_id_GlobalMuonPromptTight->at(k) == 0) LooseMuon[k]=false;
      if (mus_pt->at(k) < 10.) LooseMuon[k]=false;
      if (fabs(mus_eta->at(k)) > 2.4) LooseMuon[k]=false;
      if (fabs(d0)>=0.2) LooseMuon[k]=false;
      if (mus_cm_ndof->at(k)==0 || (mus_cm_chi2->at(k)/mus_cm_ndof->at(k)>=10)) LooseMuon[k]=false;
      if (mus_tk_numvalhits->at(k) < 11) LooseMuon[k]=false;
      // Then check tight
      TightMuon[k]=true;
      if (mus_id_GlobalMuonPromptTight->at(k) == 0) TightMuon[k]=false;
      if (mus_pt->at(k) < 20.) TightMuon[k]=false;
      if (fabs(mus_eta->at(k)) > 2.1) TightMuon[k]=false;
      if ((mus_ecalIso->at(k)+mus_hcalIso->at(k)+mus_tIso->at(k))/mus_pt->at(k)>=0.1) TightMuon[k] = false;
      if (fabs(d0)>=0.2) TightMuon[k]=false;
      if (mus_cm_ndof->at(k)==0 || mus_cm_chi2->at(k)/mus_cm_ndof->at(k) >= 10) TightMuon[k]=false;
      if (mus_tk_numvalhits->at(k) < 11) TightMuon[k]=false;
      if (mus_hcalvetoDep->at(k) >= 6.) TightMuon[k]=false;
      if (mus_ecalvetoDep->at(k) >= 4.) TightMuon[k]=false;

      int tightmu=0; int loosemu=0; int triggermu=0;
      if(TightMuon[k]) tightmu=1; if(LooseMuon[k]) loosemu=1;
      //if(TightMuon[k] || LooseMuon[k]) triggermu=1;

      if (tightmu || loosemu || triggermu) { // Muon is good enough to save
	_num_mu++;
	if(_num_mu==1) {
	  _mu1_px=mus_px->at(k);
	  _mu1_py=mus_py->at(k);
	  _mu1_pz=mus_pz->at(k);
	  _mu1_eta=mus_eta->at(k);
	  _mu1_phi=mus_phi->at(k);
	  if(_mu1_phi<0) _mu1_phi+=2*pi;
	  _mu1_pt=mus_pt->at(k);
	  _mu1_E=mus_energy->at(k);
	  _mu1_trkZatCES = -999;
	  _mu1_D0C = 0;
	  _mu1_emE = mus_iso03_emEt->at(k);
	  _mu1_hadE = mus_iso03_hadEt->at(k);
	  _mu1_isol = mus_cIso->at(k);
	  _mu1_CMUStub = 0;
	  _mu1_CMPStub = 0;
	  _mu1_CMXStub = 0;
	  _mu1_BMUStub = 0;
	  _mu1_CmuDx = -999;
	  _mu1_CmpDx = -999;
	  _mu1_CmxDx = -999;
	  _mu1_NumAxSeg = 0;
	  _mu1_NumStSeg = 0;
	  _mu1_trkchisqproblog10 = 0;
	  _mu1_trkZ0 = -999;
	  _mu1_SIHits = 0;
	  _mu1_trkisol = mus_tIso->at(k);
	  if(loosemu)_mu1_pass+=2; if(tightmu)_mu1_pass+=4; if(triggermu)_mu1_pass+=8; //if(passes_minimal==1)_mu1_pass+=1; 
	}
	if(_num_mu==2) {
	  _mu2_px=mus_px->at(k);
	  _mu2_py=mus_py->at(k);
	  _mu2_pz=mus_pz->at(k);
	  _mu2_eta=mus_eta->at(k);
	  _mu2_phi=mus_phi->at(k);
	  if(_mu2_phi<0) _mu2_phi+=2*pi;
	  _mu2_pt=mus_pt->at(k);
	  _mu2_E=mus_energy->at(k);
	  _mu2_trkZatCES = -999;
	  _mu2_D0C = 0;
	  _mu2_emE = mus_iso03_emEt->at(k);
	  _mu2_hadE = mus_iso03_hadEt->at(k);
	  _mu2_isol = mus_cIso->at(k);
	  _mu2_CMUStub = 0;
	  _mu2_CMPStub = 0;
	  _mu2_CMXStub = 0;
	  _mu2_BMUStub = 0;
	  _mu2_CmuDx = -999;
	  _mu2_CmpDx = -999;
	  _mu2_CmxDx = -999;
	  _mu2_NumAxSeg = 0;
	  _mu2_NumStSeg = 0;
	  _mu2_trkchisqproblog10 = 0;
	  _mu2_trkZ0 = -999;
	  _mu2_SIHits = 0;
	  _mu2_trkisol = mus_tIso->at(k);
	  if(loosemu)_mu2_pass+=2; if(tightmu)_mu2_pass+=4; if(triggermu)_mu2_pass+=8; //if(passes_minimal==1)_mu2_pass+=1;
	}
      } // Muon is good enough to save
    } //loop through muons

    */ //Don't look at electrons or muons



    //loop through photons
    for (unsigned int k=0; k<photons_et->size(); k++) {
      if(photons_et->at(k) <10) continue;
      //if(fabs(photons_eta->at(k))>3.0 || photons_isEBGap->at(k) || photons_isEEGap->at(k) || photons_isEBEEGap->at(k)) continue; 
      //if(photons_isEBPho->at(k) && photons_maxEnergyXtal->at(k)/photons_e3x3->at(k)>0.9) continue;  //spike rejection

      // First check loose
      LoosePho[k]=true;
      if(photons_isEBGap->at(k) || photons_isEEGap->at(k) || photons_isEBEEGap->at(k)) LoosePho[k]=false;
      if(photons_hadOverEM->at(k)>0.1) LoosePho[k]=false;
      if(photons_isoEcalRecHitDR04->at(k)>8.0+0.004*photons_et->at(k)) LoosePho[k]=false;
      if(photons_isoHcalRecHitDR04->at(k)>8.0+0.002*photons_et->at(k)) LoosePho[k]=false;
      if(photons_isoHollowTrkConeDR04->at(k)>8.0+0.002*photons_et->at(k)) LoosePho[k]=false;

      // Then check tight
      TightPho[k]=true;
      if(photons_isEBGap->at(k) || photons_isEEGap->at(k) || photons_isEBEEGap->at(k)) TightPho[k]=false;
      if(photons_hadOverEM->at(k)>0.05) TightPho[k]=false;
      if(photons_isoEcalRecHitDR04->at(k)>4.2+0.002*photons_et->at(k)) TightPho[k]=false;
      if(photons_isoHcalRecHitDR04->at(k)>4.0+0.001*photons_et->at(k)) TightPho[k]=false;
      if(photons_isoHollowTrkConeDR04->at(k)>4.0+0.001*photons_et->at(k)) TightPho[k]=false;
      if(photons_hasPixelSeed->at(k)) TightPho[k]=false;

      int tightPho=0; int loosePho=0; int triggerPho=0;
      if(TightPho[k]) tightPho=1; if(LoosePho[k]) loosePho=1;
      triggerPho=1;
      //if(TightPho[k] || LoosePho[k]) triggerPho=1;

     if (tightPho || loosePho || triggerPho) { // Photon is good enough to save
       float et = photons_et->at(k);
       float phi = photons_phi->at(k);
       float eta = photons_eta->at(k);
       _pho_eta[_num_pho] = eta;
       _pho_theta[_num_pho] = 2*atan(exp(fabs(eta)));
       _pho_phi[_num_pho] = phi;
       _pho_et[_num_pho] = et;
       _pho_px[_num_pho] = et*cos(phi);
       _pho_py[_num_pho] = et*sin(phi);
       _pho_pz[_num_pho] = et*(1./tan(_pho_theta[_num_pho]));
       _pho_E[_num_pho] = photons_energy->at(k);
       _pho_hadOverEM[_num_pho] = photons_hadOverEM->at(k);
       _pho_isoEcalRecHitDR04[_num_pho] = photons_isoEcalRecHitDR04->at(k);
       _pho_isoHcalRecHitDR04[_num_pho] = photons_isoHcalRecHitDR04->at(k);
       _pho_isoHollowTrkConeDR04[_num_pho] = photons_isoHollowTrkConeDR04->at(k);
       _pho_hasPixelSeed[_num_pho] = photons_hasPixelSeed->at(k);
       _pho_maxEnergyXtal[_num_pho] = photons_maxEnergyXtal->at(k);
       _pho_e3x3[_num_pho] = photons_e3x3->at(k);
       _pho_sigmaIetaIeta[_num_pho] = photons_sigmaIetaIeta->at(k);
       _pho_scPhiWidth[_num_pho] = photons_scPhiWidth->at(k);
       _pho_scEtaWidth[_num_pho] = photons_scEtaWidth->at(k);
       _pho_r9[_num_pho] = photons_r9->at(k);
       _pho_istight[_num_pho] = tightPho;
       _pho_isloose[_num_pho] = loosePho;
       _num_pho++;
     } // Photon is good enough to save
    } //loop through photons



    // Save this event
    if (/*_num_mu==0 && _num_ele==0 &&*/ _num_jet>0 && _num_pho>0)  { //save event if at least one photon and at leaste one jet
      tree->Fill();
      pho_tree->Fill();
    }

  } //Main event loop

  outFile->Write();
  outFile->Close();

  return 0;
}



//Strange spikes in sigmaIetaIeta
//t->Draw("pho_et","(pho_maxEnergyXtal/pho_e3x3)<0.9 && pho_istight && (pho_sigmaIetaIeta>0.01948 || pho_sigmaIetaIeta<0.01943) && (pho_sigmaIetaIeta>0.01941 || pho_sigmaIetaIeta<0.01937) && (pho_sigmaIetaIeta>0.0247 || pho_sigmaIetaIeta<0.0246) ")
