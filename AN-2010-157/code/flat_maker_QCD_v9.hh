#include "TTree.h"
#if !defined (__CINT__) || defined (__MAKECINT__)
#endif
#include "TMath.h"
#include <vector>
#include "TChain.h"

Int_t _run;
Int_t _event;
Int_t _lumiblock;
Float_t _met;
Float_t _metx;
Float_t _mety;
Float_t _metphi;
Float_t _rawmet;
Float_t _rawmetphi;
Float_t _vtxmet;
Float_t _vtxmetphi;
Int_t _num_vert;

Float_t _jet_et[64], temp_jet_et[64];
Float_t _jet_px[64], temp_jet_px[64];
Float_t _jet_py[64], temp_jet_py[64];
Float_t _jet_pz[64], temp_jet_pz[64];
Float_t _jet_eta[64], temp_jet_eta[64];
Float_t _jet_E[64], temp_jet_E[64];
Float_t _jet_theta[64], temp_jet_theta[64];
Float_t _jet_phi[64], temp_jet_phi[64];
Float_t _jet_emf[64], temp_jet_emf[64];
Float_t _jet_fHPD[64], temp_jet_fHPD[64];
Float_t _jet_fRBX[64], temp_jet_fRBX[64];
Int_t _jet_n90Hits[64], temp_jet_n90Hits[64];
Float_t _jet_btag_jetProb[64], temp_jet_btag_jetProb[64];
Float_t _jet_scalefac[64], temp_jet_scalefac[64];
Int_t _jet_remove[64];
float jet_sort[64][2];
Int_t _num_jet;

Float_t _ele1_et;
Float_t _ele1_px;
Float_t _ele1_py;
Float_t _ele1_pz;
Float_t _ele1_E;
Float_t _ele1_eta;
Float_t _ele1_phi;
Int_t _ele1_pass;
Float_t _ele2_et;
Float_t _ele2_px;
Float_t _ele2_py;
Float_t _ele2_pz;
Float_t _ele2_E;
Float_t _ele2_eta;
Float_t _ele2_phi;
Int_t _ele2_pass;
Int_t _num_ele;

Float_t _mu1_et;
Float_t _mu1_px;
Float_t _mu1_py;
Float_t _mu1_pz;
Float_t _mu1_pt;
Float_t _mu1_E;
Float_t _mu1_eta;
Float_t _mu1_phi;
Float_t _mu1_trkZatCES;
Float_t _mu1_D0C;
Float_t _mu1_emE;
Float_t _mu1_hadE;
Float_t _mu1_isol;
Float_t _mu1_CMUStub;
Float_t _mu1_CMPStub;
Float_t _mu1_CMXStub;
Float_t _mu1_BMUStub;
Float_t _mu1_CmuDx;
Float_t _mu1_CmpDx;
Float_t _mu1_CmxDx;
Float_t _mu1_NumAxSeg;
Float_t _mu1_NumStSeg;
Float_t _mu1_trkchisqproblog10;
Float_t _mu1_trkZ0;
Float_t _mu1_SIHits;
Float_t _mu1_trkisol;
Int_t _mu1_pass;
Float_t _mu2_et;
Float_t _mu2_px;
Float_t _mu2_py;
Float_t _mu2_pz;
Float_t _mu2_pt;
Float_t _mu2_E;
Float_t _mu2_eta;
Float_t _mu2_phi;
Float_t _mu2_trkZatCES;
Float_t _mu2_D0C;
Float_t _mu2_emE;
Float_t _mu2_hadE;
Float_t _mu2_isol;
Float_t _mu2_CMUStub;
Float_t _mu2_CMPStub;
Float_t _mu2_CMXStub;
Float_t _mu2_BMUStub;
Float_t _mu2_CmuDx;
Float_t _mu2_CmpDx;
Float_t _mu2_CmxDx;
Float_t _mu2_NumAxSeg;
Float_t _mu2_NumStSeg;
Float_t _mu2_trkchisqproblog10;
Float_t _mu2_trkZ0;
Float_t _mu2_SIHits;
Float_t _mu2_trkisol;
Int_t _mu2_pass;
Int_t _num_mu;

Float_t _pho_et[64], temp_pho_et[64];
Float_t _pho_px[64], temp_pho_px[64];
Float_t _pho_py[64], temp_pho_py[64];
Float_t _pho_pz[64], temp_pho_pz[64];
Float_t _pho_eta[64], temp_pho_eta[64];
Float_t _pho_E[64], temp_pho_E[64];
Float_t _pho_theta[64], temp_pho_theta[64];
Float_t _pho_phi[64], temp_pho_phi[64];
Float_t _pho_hadOverEM[64];
Float_t _pho_isoEcalRecHitDR04[64];
Float_t _pho_isoHcalRecHitDR04[64];
Float_t _pho_isoHollowTrkConeDR04[64];
Int_t _pho_hasPixelSeed[64];
Float_t _pho_maxEnergyXtal[64];
Float_t _pho_e3x3[64];
Int_t   _pho_istight[64]; 
Int_t   _pho_isloose[64]; 
float pho_sort[64][2];
Int_t _num_pho;


double deltaphi(double phi1, double phi2);
double findphi(double x, double y);

void et_sort(int num_jet){
  int imax;
  float etmax;
  //loop through jets to order by et
  for(int i=0; i<num_jet; i++) {
    float etmax=0;
    for(int j=0; j<num_jet; j++) {
	if(jet_sort[j][0]>etmax) {
	  etmax=jet_sort[j][0];
	  imax=j;
	}
    }
    jet_sort[imax][1] = i+0.1;
    jet_sort[imax][0] = -1.;
  }

  for(int i= 0; i<num_jet; i++){ //save jets in correct order
    int n;
    n = (int)jet_sort[i][1];
    _jet_et[n] = temp_jet_et[i];
    _jet_px[n] = temp_jet_px[i];
    _jet_py[n] = temp_jet_py[i];
    _jet_pz[n] = temp_jet_pz[i];
    _jet_eta[n] = temp_jet_eta[i];
    _jet_E[n] = temp_jet_E[i];
    _jet_theta[n] = temp_jet_theta[i];
    _jet_phi[n] = temp_jet_phi[i];
    _jet_emf[n] = temp_jet_emf[i];
    _jet_scalefac[n] = temp_jet_scalefac[i];  
    _jet_remove[n] = 0;
  } //save jets in correct order
}


void overlap(){
  double dR;
 
  if(_ele1_et>0) { // 1st ele
    for(int i=0; i<_num_jet; i++) { //check jets
      dR= sqrt(pow(_ele1_eta-_jet_eta[i],2) + pow(deltaphi(_ele1_phi,_jet_phi[i]),2));
      if(dR < 0.25) _jet_remove[i] = 1;
    } //check jets
  } //1st ele

  if(_ele2_et>0) { //2nd ele
    for(int i=0; i<_num_jet; i++) { //check jets
      dR= sqrt(pow(_ele2_eta-_jet_eta[i],2) + pow(deltaphi(_ele2_phi,_jet_phi[i]),2));
      if(dR < 0.25) _jet_remove[i] = 1;
    } //check jets
  } //2nd ele  
}


void met_correct(){ //correct MET for muons
  double pi = 3.141592653589;
  
  _metx -= (_mu1_pt-_mu1_emE-_mu1_hadE)*cos(_mu1_phi);
  _mety -= (_mu1_pt-_mu1_emE-_mu1_hadE)*sin(_mu1_phi);

  if(_num_mu==2){
    _metx -= (_mu2_pt-_mu2_emE-_mu2_hadE)*cos(_mu2_phi);
    _mety -= (_mu2_pt-_mu2_emE-_mu2_hadE)*sin(_mu2_phi);
  }

  _met = sqrt(pow(_metx,2)+pow(_mety,2));
  _metphi = findphi(_metx,_mety);
  if(_metphi<0) _metphi+=2*pi;

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

double deltaphi(double phi1, double phi2)
{
  double result = phi1-phi2;
  while (result>TMath::Pi()) result -= 2*TMath::Pi();
  while (result<=-TMath::Pi()) result += 2*TMath::Pi();
  return result;
}


//Version 9 of cfA

   UInt_t          NbeamSpot;
   vector<float>   *beamSpot_x;
   vector<float>   *beamSpot_y;
   vector<float>   *beamSpot_z;
   vector<float>   *beamSpot_x0Error;
   vector<float>   *beamSpot_y0Error;
   vector<float>   *beamSpot_z0Error;
   vector<float>   *beamSpot_sigmaZ;
   vector<float>   *beamSpot_sigmaZ0Error;
   vector<float>   *beamSpot_dxdz;
   vector<float>   *beamSpot_dxdzError;
   vector<float>   *beamSpot_dydz;
   vector<float>   *beamSpot_dydzError;
   vector<float>   *beamSpot_beamWidthX;
   vector<float>   *beamSpot_beamWidthY;
   vector<float>   *beamSpot_beamWidthXError;
   vector<float>   *beamSpot_beamWidthYError;
   UInt_t          Nels;
   vector<float>   *els_energy;
   vector<float>   *els_et;
   vector<float>   *els_eta;
   vector<float>   *els_phi;
   vector<float>   *els_pt;
   vector<float>   *els_px;
   vector<float>   *els_py;
   vector<float>   *els_pz;
   vector<float>   *els_status;
   vector<float>   *els_theta;
   vector<float>   *els_gen_id;
   vector<float>   *els_gen_phi;
   vector<float>   *els_gen_pt;
   vector<float>   *els_gen_pz;
   vector<float>   *els_gen_px;
   vector<float>   *els_gen_py;
   vector<float>   *els_gen_eta;
   vector<float>   *els_gen_theta;
   vector<float>   *els_gen_et;
   vector<float>   *els_gen_mother_id;
   vector<float>   *els_gen_mother_phi;
   vector<float>   *els_gen_mother_pt;
   vector<float>   *els_gen_mother_pz;
   vector<float>   *els_gen_mother_px;
   vector<float>   *els_gen_mother_py;
   vector<float>   *els_gen_mother_eta;
   vector<float>   *els_gen_mother_theta;
   vector<float>   *els_gen_mother_et;
   vector<float>   *els_tightId;
   vector<float>   *els_looseId;
   vector<float>   *els_robustTightId;
   vector<float>   *els_robustLooseId;
   vector<float>   *els_robustHighEnergyId;
   vector<float>   *els_cIso;
   vector<float>   *els_tIso;
   vector<float>   *els_ecalIso;
   vector<float>   *els_hcalIso;
   vector<float>   *els_chi2;
   vector<float>   *els_charge;
   vector<float>   *els_caloEnergy;
   vector<float>   *els_hadOverEm;
   vector<float>   *els_eOverPIn;
   vector<float>   *els_eSeedOverPOut;
   vector<float>   *els_eSCraw;
   vector<float>   *els_eSeed;
   vector<float>   *els_sigmaEtaEta;
   vector<float>   *els_sigmaIEtaIEta;
   vector<float>   *els_scE1x5;
   vector<float>   *els_scE2x5Max;
   vector<float>   *els_scE5x5;
   vector<float>   *els_dEtaIn;
   vector<float>   *els_dPhiIn;
   vector<float>   *els_dEtaOut;
   vector<float>   *els_dPhiOut;
   vector<float>   *els_numvalhits;
   vector<float>   *els_numlosthits;
   vector<float>   *els_basicClustersSize;
   vector<float>   *els_tk_pt;
   vector<float>   *els_tk_phi;
   vector<float>   *els_tk_eta;
   vector<float>   *els_d0dum;
   vector<float>   *els_dz;
   vector<float>   *els_vx;
   vector<float>   *els_vy;
   vector<float>   *els_vz;
   vector<float>   *els_ndof;
   vector<float>   *els_ptError;
   vector<float>   *els_d0dumError;
   vector<float>   *els_dzError;
   vector<float>   *els_etaError;
   vector<float>   *els_phiError;
   vector<float>   *els_tk_charge;
   vector<float>   *els_ctf_tk_id;
   vector<float>   *els_ctf_tk_charge;
   vector<float>   *els_ctf_tk_eta;
   vector<float>   *els_ctf_tk_phi;
   vector<float>   *els_fbrem;
   vector<float>   *els_shFracInnerHits;
   vector<float>   *els_dr03EcalRecHitSumEt;
   vector<float>   *els_dr03HcalTowerSumEt;
   vector<float>   *els_dr03HcalDepth1TowerSumEt;
   vector<float>   *els_dr03HcalDepth2TowerSumEt;
   vector<float>   *els_dr03TkSumPt;
   vector<float>   *els_dr04EcalRecHitSumEt;
   vector<float>   *els_dr04HcalTowerSumEt;
   vector<float>   *els_dr04HcalDepth1TowerSumEt;
   vector<float>   *els_dr04HcalDepth2TowerSumEt;
   vector<float>   *els_dr04TkSumPt;
   vector<float>   *els_cpx;
   vector<float>   *els_cpy;
   vector<float>   *els_cpz;
   vector<float>   *els_vpx;
   vector<float>   *els_vpy;
   vector<float>   *els_vpz;
   vector<float>   *els_cx;
   vector<float>   *els_cy;
   vector<float>   *els_cz;
   UInt_t          Njets_AK5;
   vector<float>   *jets_AK5_energy;
   vector<float>   *jets_AK5_et;
   vector<float>   *jets_AK5_eta;
   vector<float>   *jets_AK5_phi;
   vector<float>   *jets_AK5_pt;
   vector<float>   *jets_AK5_px;
   vector<float>   *jets_AK5_py;
   vector<float>   *jets_AK5_pz;
   vector<float>   *jets_AK5_status;
   vector<float>   *jets_AK5_theta;
   vector<float>   *jets_AK5_parton_Id;
   vector<float>   *jets_AK5_parton_motherId;
   vector<float>   *jets_AK5_parton_pt;
   vector<float>   *jets_AK5_parton_phi;
   vector<float>   *jets_AK5_parton_eta;
   vector<float>   *jets_AK5_parton_Energy;
   vector<float>   *jets_AK5_parton_mass;
   vector<float>   *jets_AK5_parton_motherID;
   vector<float>   *jets_AK5_gen_et;
   vector<float>   *jets_AK5_gen_pt;
   vector<float>   *jets_AK5_gen_eta;
   vector<float>   *jets_AK5_gen_phi;
   vector<float>   *jets_AK5_gen_mass;
   vector<float>   *jets_AK5_gen_Energy;
   vector<float>   *jets_AK5_gen_Id;
   vector<float>   *jets_AK5_gen_motherID;
   vector<float>   *jets_AK5_gen_threeCharge;
   vector<float>   *jets_AK5_partonFlavour;
   vector<float>   *jets_AK5_btag_TC_highPur;
   vector<float>   *jets_AK5_btag_TC_highEff;
   vector<float>   *jets_AK5_btag_jetProb;
   vector<float>   *jets_AK5_btag_jetBProb;
   vector<float>   *jets_AK5_btag_softEle;
   vector<float>   *jets_AK5_btag_softMuon;
   vector<float>   *jets_AK5_btag_softMuonNoIP;
   vector<float>   *jets_AK5_btag_secVertex;
   vector<float>   *jets_AK5_chgEmE;
   vector<float>   *jets_AK5_chgHadE;
   vector<float>   *jets_AK5_chgMuE;
   vector<float>   *jets_AK5_chg_Mult;
   vector<float>   *jets_AK5_neutralEmE;
   vector<float>   *jets_AK5_neutralHadE;
   vector<float>   *jets_AK5_neutral_Mult;
   vector<float>   *jets_AK5_mu_Mult;
   vector<float>   *jets_AK5_emf;
   vector<float>   *jets_AK5_ehf;
   vector<float>   *jets_AK5_n60;
   vector<float>   *jets_AK5_n90;
   vector<float>   *jets_AK5_area;
   vector<float>   *jets_AK5_mass;
   UInt_t          Njets_SC5;
   vector<float>   *jets_SC5_energy;
   vector<float>   *jets_SC5_et;
   vector<float>   *jets_SC5_eta;
   vector<float>   *jets_SC5_phi;
   vector<float>   *jets_SC5_pt;
   vector<float>   *jets_SC5_px;
   vector<float>   *jets_SC5_py;
   vector<float>   *jets_SC5_pz;
   vector<float>   *jets_SC5_status;
   vector<float>   *jets_SC5_theta;
   vector<float>   *jets_SC5_parton_Id;
   vector<float>   *jets_SC5_parton_motherId;
   vector<float>   *jets_SC5_parton_pt;
   vector<float>   *jets_SC5_parton_phi;
   vector<float>   *jets_SC5_parton_eta;
   vector<float>   *jets_SC5_parton_Energy;
   vector<float>   *jets_SC5_parton_mass;
   vector<float>   *jets_SC5_parton_motherID;
   vector<float>   *jets_SC5_gen_et;
   vector<float>   *jets_SC5_gen_pt;
   vector<float>   *jets_SC5_gen_eta;
   vector<float>   *jets_SC5_gen_phi;
   vector<float>   *jets_SC5_gen_mass;
   vector<float>   *jets_SC5_gen_Energy;
   vector<float>   *jets_SC5_gen_Id;
   vector<float>   *jets_SC5_gen_motherID;
   vector<float>   *jets_SC5_gen_threeCharge;
   vector<float>   *jets_SC5_partonFlavour;
   vector<float>   *jets_SC5_btag_TC_highPur;
   vector<float>   *jets_SC5_btag_TC_highEff;
   vector<float>   *jets_SC5_btag_jetProb;
   vector<float>   *jets_SC5_btag_jetBProb;
   vector<float>   *jets_SC5_btag_softEle;
   vector<float>   *jets_SC5_btag_softMuon;
   vector<float>   *jets_SC5_btag_softMuonNoIP;
   vector<float>   *jets_SC5_btag_secVertex;
   vector<float>   *jets_SC5_chgEmE;
   vector<float>   *jets_SC5_chgHadE;
   vector<float>   *jets_SC5_chgMuE;
   vector<float>   *jets_SC5_chg_Mult;
   vector<float>   *jets_SC5_neutralEmE;
   vector<float>   *jets_SC5_neutralHadE;
   vector<float>   *jets_SC5_neutral_Mult;
   vector<float>   *jets_SC5_mu_Mult;
   vector<float>   *jets_SC5_emf;
   vector<float>   *jets_SC5_ehf;
   vector<float>   *jets_SC5_n60;
   vector<float>   *jets_SC5_n90;
   vector<float>   *jets_SC5_area;
   vector<float>   *jets_SC5_mass;
   UInt_t          Nmc_doc;
   vector<float>   *mc_doc_id;
   vector<float>   *mc_doc_pt;
   vector<float>   *mc_doc_px;
   vector<float>   *mc_doc_py;
   vector<float>   *mc_doc_pz;
   vector<float>   *mc_doc_eta;
   vector<float>   *mc_doc_phi;
   vector<float>   *mc_doc_theta;
   vector<float>   *mc_doc_energy;
   vector<float>   *mc_doc_status;
   vector<float>   *mc_doc_charge;
   vector<float>   *mc_doc_mother_id;
   vector<float>   *mc_doc_grandmother_id;
   vector<float>   *mc_doc_ggrandmother_id;
   vector<float>   *mc_doc_mother_pt;
   vector<float>   *mc_doc_vertex_x;
   vector<float>   *mc_doc_vertex_y;
   vector<float>   *mc_doc_vertex_z;
   vector<float>   *mc_doc_mass;
   vector<float>   *mc_doc_numOfDaughters;
   vector<float>   *mc_doc_numOfMothers;
   UInt_t          Nmc_electrons;
   vector<float>   *mc_electrons_id;
   vector<float>   *mc_electrons_pt;
   vector<float>   *mc_electrons_px;
   vector<float>   *mc_electrons_py;
   vector<float>   *mc_electrons_pz;
   vector<float>   *mc_electrons_eta;
   vector<float>   *mc_electrons_phi;
   vector<float>   *mc_electrons_theta;
   vector<float>   *mc_electrons_status;
   vector<float>   *mc_electrons_energy;
   vector<float>   *mc_electrons_charge;
   vector<float>   *mc_electrons_mother_id;
   vector<float>   *mc_electrons_mother_pt;
   vector<float>   *mc_electrons_grandmother_id;
   vector<float>   *mc_electrons_ggrandmother_id;
   vector<float>   *mc_electrons_vertex_x;
   vector<float>   *mc_electrons_vertex_y;
   vector<float>   *mc_electrons_vertex_z;
   vector<float>   *mc_electrons_mass;
   vector<float>   *mc_electrons_numOfDaughters;
   UInt_t          Nmc_mus;
   vector<float>   *mc_mus_id;
   vector<float>   *mc_mus_pt;
   vector<float>   *mc_mus_px;
   vector<float>   *mc_mus_py;
   vector<float>   *mc_mus_pz;
   vector<float>   *mc_mus_eta;
   vector<float>   *mc_mus_phi;
   vector<float>   *mc_mus_theta;
   vector<float>   *mc_mus_status;
   vector<float>   *mc_mus_energy;
   vector<float>   *mc_mus_charge;
   vector<float>   *mc_mus_mother_id;
   vector<float>   *mc_mus_mother_pt;
   vector<float>   *mc_mus_grandmother_id;
   vector<float>   *mc_mus_ggrandmother_id;
   vector<float>   *mc_mus_vertex_x;
   vector<float>   *mc_mus_vertex_y;
   vector<float>   *mc_mus_vertex_z;
   vector<float>   *mc_mus_mass;
   vector<float>   *mc_mus_numOfDaughters;
   UInt_t          Nmets_AK5;
   vector<float>   *mets_AK5_et;
   vector<float>   *mets_AK5_phi;
   vector<float>   *mets_AK5_ex;
   vector<float>   *mets_AK5_ey;
   vector<float>   *mets_AK5_gen_et;
   vector<float>   *mets_AK5_gen_phi;
   vector<float>   *mets_AK5_sign;
   vector<float>   *mets_AK5_sumEt;
   vector<float>   *mets_AK5_unCPhi;
   vector<float>   *mets_AK5_unCPt;
   UInt_t          Nmets_SC5;
   vector<float>   *mets_SC5_et;
   vector<float>   *mets_SC5_phi;
   vector<float>   *mets_SC5_ex;
   vector<float>   *mets_SC5_ey;
   vector<float>   *mets_SC5_gen_et;
   vector<float>   *mets_SC5_gen_phi;
   vector<float>   *mets_SC5_sign;
   vector<float>   *mets_SC5_sumEt;
   vector<float>   *mets_SC5_unCPhi;
   vector<float>   *mets_SC5_unCPt;
   UInt_t          Nmus;
   vector<float>   *mus_energy;
   vector<float>   *mus_et;
   vector<float>   *mus_eta;
   vector<float>   *mus_phi;
   vector<float>   *mus_pt;
   vector<float>   *mus_px;
   vector<float>   *mus_py;
   vector<float>   *mus_pz;
   vector<float>   *mus_status;
   vector<float>   *mus_theta;
   vector<float>   *mus_gen_id;
   vector<float>   *mus_gen_phi;
   vector<float>   *mus_gen_pt;
   vector<float>   *mus_gen_pz;
   vector<float>   *mus_gen_px;
   vector<float>   *mus_gen_py;
   vector<float>   *mus_gen_eta;
   vector<float>   *mus_gen_theta;
   vector<float>   *mus_gen_et;
   vector<float>   *mus_gen_mother_id;
   vector<float>   *mus_gen_mother_phi;
   vector<float>   *mus_gen_mother_pt;
   vector<float>   *mus_gen_mother_pz;
   vector<float>   *mus_gen_mother_px;
   vector<float>   *mus_gen_mother_py;
   vector<float>   *mus_gen_mother_eta;
   vector<float>   *mus_gen_mother_theta;
   vector<float>   *mus_gen_mother_et;
   vector<float>   *mus_tkHits;
   vector<float>   *mus_cIso;
   vector<float>   *mus_tIso;
   vector<float>   *mus_ecalIso;
   vector<float>   *mus_hcalIso;
   vector<float>   *mus_ecalvetoDep;
   vector<float>   *mus_hcalvetoDep;
   vector<float>   *mus_calEnergyEm;
   vector<float>   *mus_calEnergyHad;
   vector<float>   *mus_calEnergyHo;
   vector<float>   *mus_calEnergyEmS9;
   vector<float>   *mus_calEnergyHadS9;
   vector<float>   *mus_calEnergyHoS9;
   vector<float>   *mus_iso03_sumPt;
   vector<float>   *mus_iso03_emEt;
   vector<float>   *mus_iso03_hadEt;
   vector<float>   *mus_iso03_hoEt;
   vector<float>   *mus_iso03_nTracks;
   vector<float>   *mus_iso05_sumPt;
   vector<float>   *mus_iso05_emEt;
   vector<float>   *mus_iso05_hadEt;
   vector<float>   *mus_iso05_hoEt;
   vector<float>   *mus_iso05_nTracks;
   vector<float>   *mus_charge;
   vector<float>   *mus_cm_chi2;
   vector<float>   *mus_cm_ndof;
   vector<float>   *mus_cm_chg;
   vector<float>   *mus_cm_pt;
   vector<float>   *mus_cm_px;
   vector<float>   *mus_cm_py;
   vector<float>   *mus_cm_pz;
   vector<float>   *mus_cm_eta;
   vector<float>   *mus_cm_phi;
   vector<float>   *mus_cm_theta;
   vector<float>   *mus_cm_d0dum;
   vector<float>   *mus_cm_dz;
   vector<float>   *mus_cm_vx;
   vector<float>   *mus_cm_vy;
   vector<float>   *mus_cm_vz;
   vector<float>   *mus_cm_numvalhits;
   vector<float>   *mus_cm_numlosthits;
   vector<float>   *mus_cm_d0dumErr;
   vector<float>   *mus_cm_dzErr;
   vector<float>   *mus_cm_ptErr;
   vector<float>   *mus_cm_etaErr;
   vector<float>   *mus_cm_phiErr;
   vector<float>   *mus_tk_id;
   vector<float>   *mus_tk_chi2;
   vector<float>   *mus_tk_ndof;
   vector<float>   *mus_tk_chg;
   vector<float>   *mus_tk_pt;
   vector<float>   *mus_tk_px;
   vector<float>   *mus_tk_py;
   vector<float>   *mus_tk_pz;
   vector<float>   *mus_tk_eta;
   vector<float>   *mus_tk_phi;
   vector<float>   *mus_tk_theta;
   vector<float>   *mus_tk_d0dum;
   vector<float>   *mus_tk_dz;
   vector<float>   *mus_tk_vx;
   vector<float>   *mus_tk_vy;
   vector<float>   *mus_tk_vz;
   vector<float>   *mus_tk_numvalhits;
   vector<float>   *mus_tk_numlosthits;
   vector<float>   *mus_tk_d0dumErr;
   vector<float>   *mus_tk_dzErr;
   vector<float>   *mus_tk_ptErr;
   vector<float>   *mus_tk_etaErr;
   vector<float>   *mus_tk_phiErr;
   vector<float>   *mus_stamu_chi2;
   vector<float>   *mus_stamu_ndof;
   vector<float>   *mus_stamu_chg;
   vector<float>   *mus_stamu_pt;
   vector<float>   *mus_stamu_px;
   vector<float>   *mus_stamu_py;
   vector<float>   *mus_stamu_pz;
   vector<float>   *mus_stamu_eta;
   vector<float>   *mus_stamu_phi;
   vector<float>   *mus_stamu_theta;
   vector<float>   *mus_stamu_d0dum;
   vector<float>   *mus_stamu_dz;
   vector<float>   *mus_stamu_vx;
   vector<float>   *mus_stamu_vy;
   vector<float>   *mus_stamu_vz;
   vector<float>   *mus_stamu_numvalhits;
   vector<float>   *mus_stamu_numlosthits;
   vector<float>   *mus_stamu_d0dumErr;
   vector<float>   *mus_stamu_dzErr;
   vector<float>   *mus_stamu_ptErr;
   vector<float>   *mus_stamu_etaErr;
   vector<float>   *mus_stamu_phiErr;
   vector<float>   *mus_num_matches;
   vector<float>   *mus_id_All;
   vector<float>   *mus_id_AllGlobalMuons;
   vector<float>   *mus_id_AllStandAloneMuons;
   vector<float>   *mus_id_AllTrackerMuons;
   vector<float>   *mus_id_TrackerMuonArbitrated;
   vector<float>   *mus_id_AllArbitrated;
   vector<float>   *mus_id_GlobalMuonPromptTight;
   vector<float>   *mus_id_TMLastStationLoose;
   vector<float>   *mus_id_TMLastStationTight;
   vector<float>   *mus_id_TM2DCompatibilityLoose;
   vector<float>   *mus_id_TM2DCompatibilityTight;
   vector<float>   *mus_id_TMOneStationLoose;
   vector<float>   *mus_id_TMOneStationTight;
   vector<float>   *mus_id_TMLastStationOptimizedLowPtLoose;
   vector<float>   *mus_id_TMLastStationOptimizedLowPtTight;
   UInt_t          Nphotons;
   vector<float>   *photons_energy;
   vector<float>   *photons_et;
   vector<float>   *photons_eta;
   vector<float>   *photons_phi;
   vector<float>   *photons_pt;
   vector<float>   *photons_px;
   vector<float>   *photons_py;
   vector<float>   *photons_pz;
   vector<float>   *photons_status;
   vector<float>   *photons_theta;
   vector<float>   *photons_hadOverEM;
   vector<float>   *photons_scEnergy;
   vector<float>   *photons_scRawEnergy;
   vector<float>   *photons_scEta;
   vector<float>   *photons_scPhi;
   vector<float>   *photons_scEtaWidth;
   vector<float>   *photons_scPhiWidth;
   vector<float>   *photons_tIso;
   vector<float>   *photons_ecalIso;
   vector<float>   *photons_hcalIso;
   vector<float>   *photons_isoEcalRecHitDR04;
   vector<float>   *photons_isoHcalRecHitDR04;
   vector<float>   *photons_isoSolidTrkConeDR04;
   vector<float>   *photons_isoHollowTrkConeDR04;
   vector<float>   *photons_nTrkSolidConeDR04;
   vector<float>   *photons_nTrkHollowConeDR04;
   vector<float>   *photons_isoEcalRecHitDR03;
   vector<float>   *photons_isoHcalRecHitDR03;
   vector<float>   *photons_isoSolidTrkConeDR03;
   vector<float>   *photons_isoHollowTrkConeDR03;
   vector<float>   *photons_nTrkSolidConeDR03;
   vector<float>   *photons_nTrkHollowConeDR03;
   vector<float>   *photons_isAlsoElectron;
   vector<float>   *photons_hasPixelSeed;
   vector<float>   *photons_isConverted;
   vector<float>   *photons_isEBGap;
   vector<float>   *photons_isEEGap;
   vector<float>   *photons_isEBEEGap;
   vector<float>   *photons_isEBPho;
   vector<float>   *photons_isEEPho;
   vector<float>   *photons_isLoosePhoton;
   vector<float>   *photons_isTightPhoton;
   vector<float>   *photons_r9;
   vector<float>   *photons_gen_et;
   vector<float>   *photons_gen_eta;
   vector<float>   *photons_gen_phi;
   vector<float>   *photons_gen_id;
   UInt_t          Npv;
   vector<float>   *pv_x;
   vector<float>   *pv_y;
   vector<float>   *pv_z;
   vector<float>   *pv_xErr;
   vector<float>   *pv_yErr;
   vector<float>   *pv_zErr;
   UInt_t          Ntaus;
   vector<float>   *taus_energy;
   vector<float>   *taus_et;
   vector<float>   *taus_eta;
   vector<float>   *taus_phi;
   vector<float>   *taus_pt;
   vector<float>   *taus_px;
   vector<float>   *taus_py;
   vector<float>   *taus_pz;
   vector<float>   *taus_status;
   vector<float>   *taus_theta;
   vector<float>   *taus_charge;
   vector<float>   *taus_emf;
   vector<float>   *taus_hcalTotOverPLead;
   vector<float>   *taus_hcalMaxOverPLead;
   vector<float>   *taus_hcal3x3OverPLead;
   vector<float>   *taus_ecalStripSumEOverPLead;
   vector<float>   *taus_elecPreIdOutput;
   vector<float>   *taus_elecPreIdDecision;
   vector<float>   *taus_leadPFChargedHadrCand_pt;
   vector<float>   *taus_leadPFChargedHadrCand_charge;
   vector<float>   *taus_leadPFChargedHadrCand_eta;
   vector<float>   *taus_leadPFChargedHadrCand_ECAL_eta;
   vector<float>   *taus_leadPFChargedHadrCand_phi;
   vector<float>   *taus_isoPFGammaCandsEtSum;
   vector<float>   *taus_isoPFChargedHadrCandsPtSum;
   vector<float>   *taus_leadingTrackFinding;
   vector<float>   *taus_leadingTrackPtCut;
   vector<float>   *taus_trackIsolation;
   vector<float>   *taus_ecalIsolation;
   vector<float>   *taus_byIsolation;
   vector<float>   *taus_againstElectron;
   vector<float>   *taus_againstMuon;
   vector<float>   *taus_taNC_quarter;
   vector<float>   *taus_taNC_one;
   vector<float>   *taus_taNC_half;
   vector<float>   *taus_taNC_tenth;
   vector<float>   *taus_muDecision;
   vector<float>   *taus_Nprongs;
   UInt_t          Ntcmets;
   vector<float>   *tcmets_et;
   vector<float>   *tcmets_phi;
   vector<float>   *tcmets_ex;
   vector<float>   *tcmets_ey;
   vector<float>   *tcmets_sumEt;
   UInt_t          Ntracks;
   vector<float>   *tracks_chi2;
   vector<float>   *tracks_ndof;
   vector<float>   *tracks_chg;
   vector<float>   *tracks_pt;
   vector<float>   *tracks_px;
   vector<float>   *tracks_py;
   vector<float>   *tracks_pz;
   vector<float>   *tracks_eta;
   vector<float>   *tracks_phi;
   vector<float>   *tracks_theta;
   vector<float>   *tracks_d0dum;
   vector<float>   *tracks_dz;
   vector<float>   *tracks_vx;
   vector<float>   *tracks_vy;
   vector<float>   *tracks_vz;
   vector<float>   *tracks_numvalhits;
   vector<float>   *tracks_numlosthits;
   vector<float>   *tracks_d0dumErr;
   vector<float>   *tracks_dzErr;
   vector<float>   *tracks_ptErr;
   vector<float>   *tracks_etaErr;
   vector<float>   *tracks_phiErr;
   vector<float>   *tracks_Nrechits;
   vector<float>   *tracks_innerHitX;
   vector<float>   *tracks_innerHitY;
   vector<float>   *tracks_innerHitZ;
   vector<float>   *tracks_outerHitX;
   vector<float>   *tracks_outerHitY;
   vector<float>   *tracks_outerHitZ;
   vector<float>   *tracks_highPurity;
   UInt_t          run;
   UInt_t          event;

   // List of branches
   TBranch        *b_NbeamSpot;   //!
   TBranch        *b_beamSpot_x;   //!
   TBranch        *b_beamSpot_y;   //!
   TBranch        *b_beamSpot_z;   //!
   TBranch        *b_beamSpot_x0Error;   //!
   TBranch        *b_beamSpot_y0Error;   //!
   TBranch        *b_beamSpot_z0Error;   //!
   TBranch        *b_beamSpot_sigmaZ;   //!
   TBranch        *b_beamSpot_sigmaZ0Error;   //!
   TBranch        *b_beamSpot_dxdz;   //!
   TBranch        *b_beamSpot_dxdzError;   //!
   TBranch        *b_beamSpot_dydz;   //!
   TBranch        *b_beamSpot_dydzError;   //!
   TBranch        *b_beamSpot_beamWidthX;   //!
   TBranch        *b_beamSpot_beamWidthY;   //!
   TBranch        *b_beamSpot_beamWidthXError;   //!
   TBranch        *b_beamSpot_beamWidthYError;   //!
   TBranch        *b_Nels;   //!
   TBranch        *b_els_energy;   //!
   TBranch        *b_els_et;   //!
   TBranch        *b_els_eta;   //!
   TBranch        *b_els_phi;   //!
   TBranch        *b_els_pt;   //!
   TBranch        *b_els_px;   //!
   TBranch        *b_els_py;   //!
   TBranch        *b_els_pz;   //!
   TBranch        *b_els_status;   //!
   TBranch        *b_els_theta;   //!
   TBranch        *b_els_gen_id;   //!
   TBranch        *b_els_gen_phi;   //!
   TBranch        *b_els_gen_pt;   //!
   TBranch        *b_els_gen_pz;   //!
   TBranch        *b_els_gen_px;   //!
   TBranch        *b_els_gen_py;   //!
   TBranch        *b_els_gen_eta;   //!
   TBranch        *b_els_gen_theta;   //!
   TBranch        *b_els_gen_et;   //!
   TBranch        *b_els_gen_mother_id;   //!
   TBranch        *b_els_gen_mother_phi;   //!
   TBranch        *b_els_gen_mother_pt;   //!
   TBranch        *b_els_gen_mother_pz;   //!
   TBranch        *b_els_gen_mother_px;   //!
   TBranch        *b_els_gen_mother_py;   //!
   TBranch        *b_els_gen_mother_eta;   //!
   TBranch        *b_els_gen_mother_theta;   //!
   TBranch        *b_els_gen_mother_et;   //!
   TBranch        *b_els_tightId;   //!
   TBranch        *b_els_looseId;   //!
   TBranch        *b_els_robustTightId;   //!
   TBranch        *b_els_robustLooseId;   //!
   TBranch        *b_els_robustHighEnergyId;   //!
   TBranch        *b_els_cIso;   //!
   TBranch        *b_els_tIso;   //!
   TBranch        *b_els_ecalIso;   //!
   TBranch        *b_els_hcalIso;   //!
   TBranch        *b_els_chi2;   //!
   TBranch        *b_els_charge;   //!
   TBranch        *b_els_caloEnergy;   //!
   TBranch        *b_els_hadOverEm;   //!
   TBranch        *b_els_eOverPIn;   //!
   TBranch        *b_els_eSeedOverPOut;   //!
   TBranch        *b_els_eSCraw;   //!
   TBranch        *b_els_eSeed;   //!
   TBranch        *b_els_sigmaEtaEta;   //!
   TBranch        *b_els_sigmaIEtaIEta;   //!
   TBranch        *b_els_scE1x5;   //!
   TBranch        *b_els_scE2x5Max;   //!
   TBranch        *b_els_scE5x5;   //!
   TBranch        *b_els_dEtaIn;   //!
   TBranch        *b_els_dPhiIn;   //!
   TBranch        *b_els_dEtaOut;   //!
   TBranch        *b_els_dPhiOut;   //!
   TBranch        *b_els_numvalhits;   //!
   TBranch        *b_els_numlosthits;   //!
   TBranch        *b_els_basicClustersSize;   //!
   TBranch        *b_els_tk_pt;   //!
   TBranch        *b_els_tk_phi;   //!
   TBranch        *b_els_tk_eta;   //!
   TBranch        *b_els_d0dum;   //!
   TBranch        *b_els_dz;   //!
   TBranch        *b_els_vx;   //!
   TBranch        *b_els_vy;   //!
   TBranch        *b_els_vz;   //!
   TBranch        *b_els_ndof;   //!
   TBranch        *b_els_ptError;   //!
   TBranch        *b_els_d0dumError;   //!
   TBranch        *b_els_dzError;   //!
   TBranch        *b_els_etaError;   //!
   TBranch        *b_els_phiError;   //!
   TBranch        *b_els_tk_charge;   //!
   TBranch        *b_els_ctf_tk_id;   //!
   TBranch        *b_els_ctf_tk_charge;   //!
   TBranch        *b_els_ctf_tk_eta;   //!
   TBranch        *b_els_ctf_tk_phi;   //!
   TBranch        *b_els_fbrem;   //!
   TBranch        *b_els_shFracInnerHits;   //!
   TBranch        *b_els_dr03EcalRecHitSumEt;   //!
   TBranch        *b_els_dr03HcalTowerSumEt;   //!
   TBranch        *b_els_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_els_dr03HcalDepth2TowerSumEt;   //!
   TBranch        *b_els_dr03TkSumPt;   //!
   TBranch        *b_els_dr04EcalRecHitSumEt;   //!
   TBranch        *b_els_dr04HcalTowerSumEt;   //!
   TBranch        *b_els_dr04HcalDepth1TowerSumEt;   //!
   TBranch        *b_els_dr04HcalDepth2TowerSumEt;   //!
   TBranch        *b_els_dr04TkSumPt;   //!
   TBranch        *b_els_cpx;   //!
   TBranch        *b_els_cpy;   //!
   TBranch        *b_els_cpz;   //!
   TBranch        *b_els_vpx;   //!
   TBranch        *b_els_vpy;   //!
   TBranch        *b_els_vpz;   //!
   TBranch        *b_els_cx;   //!
   TBranch        *b_els_cy;   //!
   TBranch        *b_els_cz;   //!
   TBranch        *b_Njets_AK5;   //!
   TBranch        *b_jets_AK5_energy;   //!
   TBranch        *b_jets_AK5_et;   //!
   TBranch        *b_jets_AK5_eta;   //!
   TBranch        *b_jets_AK5_phi;   //!
   TBranch        *b_jets_AK5_pt;   //!
   TBranch        *b_jets_AK5_px;   //!
   TBranch        *b_jets_AK5_py;   //!
   TBranch        *b_jets_AK5_pz;   //!
   TBranch        *b_jets_AK5_status;   //!
   TBranch        *b_jets_AK5_theta;   //!
   TBranch        *b_jets_AK5_parton_Id;   //!
   TBranch        *b_jets_AK5_parton_motherId;   //!
   TBranch        *b_jets_AK5_parton_pt;   //!
   TBranch        *b_jets_AK5_parton_phi;   //!
   TBranch        *b_jets_AK5_parton_eta;   //!
   TBranch        *b_jets_AK5_parton_Energy;   //!
   TBranch        *b_jets_AK5_parton_mass;   //!
   TBranch        *b_jets_AK5_parton_motherID;   //!
   TBranch        *b_jets_AK5_gen_et;   //!
   TBranch        *b_jets_AK5_gen_pt;   //!
   TBranch        *b_jets_AK5_gen_eta;   //!
   TBranch        *b_jets_AK5_gen_phi;   //!
   TBranch        *b_jets_AK5_gen_mass;   //!
   TBranch        *b_jets_AK5_gen_Energy;   //!
   TBranch        *b_jets_AK5_gen_Id;   //!
   TBranch        *b_jets_AK5_gen_motherID;   //!
   TBranch        *b_jets_AK5_gen_threeCharge;   //!
   TBranch        *b_jets_AK5_partonFlavour;   //!
   TBranch        *b_jets_AK5_btag_TC_highPur;   //!
   TBranch        *b_jets_AK5_btag_TC_highEff;   //!
   TBranch        *b_jets_AK5_btag_jetProb;   //!
   TBranch        *b_jets_AK5_btag_jetBProb;   //!
   TBranch        *b_jets_AK5_btag_softEle;   //!
   TBranch        *b_jets_AK5_btag_softMuon;   //!
   TBranch        *b_jets_AK5_btag_softMuonNoIP;   //!
   TBranch        *b_jets_AK5_btag_secVertex;   //!
   TBranch        *b_jets_AK5_chgEmE;   //!
   TBranch        *b_jets_AK5_chgHadE;   //!
   TBranch        *b_jets_AK5_chgMuE;   //!
   TBranch        *b_jets_AK5_chg_Mult;   //!
   TBranch        *b_jets_AK5_neutralEmE;   //!
   TBranch        *b_jets_AK5_neutralHadE;   //!
   TBranch        *b_jets_AK5_neutral_Mult;   //!
   TBranch        *b_jets_AK5_mu_Mult;   //!
   TBranch        *b_jets_AK5_emf;   //!
   TBranch        *b_jets_AK5_ehf;   //!
   TBranch        *b_jets_AK5_n60;   //!
   TBranch        *b_jets_AK5_n90;   //!
   TBranch        *b_jets_AK5_area;   //!
   TBranch        *b_jets_AK5_mass;   //!
   TBranch        *b_Njets_SC5;   //!
   TBranch        *b_jets_SC5_energy;   //!
   TBranch        *b_jets_SC5_et;   //!
   TBranch        *b_jets_SC5_eta;   //!
   TBranch        *b_jets_SC5_phi;   //!
   TBranch        *b_jets_SC5_pt;   //!
   TBranch        *b_jets_SC5_px;   //!
   TBranch        *b_jets_SC5_py;   //!
   TBranch        *b_jets_SC5_pz;   //!
   TBranch        *b_jets_SC5_status;   //!
   TBranch        *b_jets_SC5_theta;   //!
   TBranch        *b_jets_SC5_parton_Id;   //!
   TBranch        *b_jets_SC5_parton_motherId;   //!
   TBranch        *b_jets_SC5_parton_pt;   //!
   TBranch        *b_jets_SC5_parton_phi;   //!
   TBranch        *b_jets_SC5_parton_eta;   //!
   TBranch        *b_jets_SC5_parton_Energy;   //!
   TBranch        *b_jets_SC5_parton_mass;   //!
   TBranch        *b_jets_SC5_parton_motherID;   //!
   TBranch        *b_jets_SC5_gen_et;   //!
   TBranch        *b_jets_SC5_gen_pt;   //!
   TBranch        *b_jets_SC5_gen_eta;   //!
   TBranch        *b_jets_SC5_gen_phi;   //!
   TBranch        *b_jets_SC5_gen_mass;   //!
   TBranch        *b_jets_SC5_gen_Energy;   //!
   TBranch        *b_jets_SC5_gen_Id;   //!
   TBranch        *b_jets_SC5_gen_motherID;   //!
   TBranch        *b_jets_SC5_gen_threeCharge;   //!
   TBranch        *b_jets_SC5_partonFlavour;   //!
   TBranch        *b_jets_SC5_btag_TC_highPur;   //!
   TBranch        *b_jets_SC5_btag_TC_highEff;   //!
   TBranch        *b_jets_SC5_btag_jetProb;   //!
   TBranch        *b_jets_SC5_btag_jetBProb;   //!
   TBranch        *b_jets_SC5_btag_softEle;   //!
   TBranch        *b_jets_SC5_btag_softMuon;   //!
   TBranch        *b_jets_SC5_btag_softMuonNoIP;   //!
   TBranch        *b_jets_SC5_btag_secVertex;   //!
   TBranch        *b_jets_SC5_chgEmE;   //!
   TBranch        *b_jets_SC5_chgHadE;   //!
   TBranch        *b_jets_SC5_chgMuE;   //!
   TBranch        *b_jets_SC5_chg_Mult;   //!
   TBranch        *b_jets_SC5_neutralEmE;   //!
   TBranch        *b_jets_SC5_neutralHadE;   //!
   TBranch        *b_jets_SC5_neutral_Mult;   //!
   TBranch        *b_jets_SC5_mu_Mult;   //!
   TBranch        *b_jets_SC5_emf;   //!
   TBranch        *b_jets_SC5_ehf;   //!
   TBranch        *b_jets_SC5_n60;   //!
   TBranch        *b_jets_SC5_n90;   //!
   TBranch        *b_jets_SC5_area;   //!
   TBranch        *b_jets_SC5_mass;   //!
   TBranch        *b_Nmc_doc;   //!
   TBranch        *b_mc_doc_id;   //!
   TBranch        *b_mc_doc_pt;   //!
   TBranch        *b_mc_doc_px;   //!
   TBranch        *b_mc_doc_py;   //!
   TBranch        *b_mc_doc_pz;   //!
   TBranch        *b_mc_doc_eta;   //!
   TBranch        *b_mc_doc_phi;   //!
   TBranch        *b_mc_doc_theta;   //!
   TBranch        *b_mc_doc_energy;   //!
   TBranch        *b_mc_doc_status;   //!
   TBranch        *b_mc_doc_charge;   //!
   TBranch        *b_mc_doc_mother_id;   //!
   TBranch        *b_mc_doc_grandmother_id;   //!
   TBranch        *b_mc_doc_ggrandmother_id;   //!
   TBranch        *b_mc_doc_mother_pt;   //!
   TBranch        *b_mc_doc_vertex_x;   //!
   TBranch        *b_mc_doc_vertex_y;   //!
   TBranch        *b_mc_doc_vertex_z;   //!
   TBranch        *b_mc_doc_mass;   //!
   TBranch        *b_mc_doc_numOfDaughters;   //!
   TBranch        *b_mc_doc_numOfMothers;   //!
   TBranch        *b_Nmc_electrons;   //!
   TBranch        *b_mc_electrons_id;   //!
   TBranch        *b_mc_electrons_pt;   //!
   TBranch        *b_mc_electrons_px;   //!
   TBranch        *b_mc_electrons_py;   //!
   TBranch        *b_mc_electrons_pz;   //!
   TBranch        *b_mc_electrons_eta;   //!
   TBranch        *b_mc_electrons_phi;   //!
   TBranch        *b_mc_electrons_theta;   //!
   TBranch        *b_mc_electrons_status;   //!
   TBranch        *b_mc_electrons_energy;   //!
   TBranch        *b_mc_electrons_charge;   //!
   TBranch        *b_mc_electrons_mother_id;   //!
   TBranch        *b_mc_electrons_mother_pt;   //!
   TBranch        *b_mc_electrons_grandmother_id;   //!
   TBranch        *b_mc_electrons_ggrandmother_id;   //!
   TBranch        *b_mc_electrons_vertex_x;   //!
   TBranch        *b_mc_electrons_vertex_y;   //!
   TBranch        *b_mc_electrons_vertex_z;   //!
   TBranch        *b_mc_electrons_mass;   //!
   TBranch        *b_mc_electrons_numOfDaughters;   //!
   TBranch        *b_Nmc_mus;   //!
   TBranch        *b_mc_mus_id;   //!
   TBranch        *b_mc_mus_pt;   //!
   TBranch        *b_mc_mus_px;   //!
   TBranch        *b_mc_mus_py;   //!
   TBranch        *b_mc_mus_pz;   //!
   TBranch        *b_mc_mus_eta;   //!
   TBranch        *b_mc_mus_phi;   //!
   TBranch        *b_mc_mus_theta;   //!
   TBranch        *b_mc_mus_status;   //!
   TBranch        *b_mc_mus_energy;   //!
   TBranch        *b_mc_mus_charge;   //!
   TBranch        *b_mc_mus_mother_id;   //!
   TBranch        *b_mc_mus_mother_pt;   //!
   TBranch        *b_mc_mus_grandmother_id;   //!
   TBranch        *b_mc_mus_ggrandmother_id;   //!
   TBranch        *b_mc_mus_vertex_x;   //!
   TBranch        *b_mc_mus_vertex_y;   //!
   TBranch        *b_mc_mus_vertex_z;   //!
   TBranch        *b_mc_mus_mass;   //!
   TBranch        *b_mc_mus_numOfDaughters;   //!
   TBranch        *b_Nmets_AK5;   //!
   TBranch        *b_mets_AK5_et;   //!
   TBranch        *b_mets_AK5_phi;   //!
   TBranch        *b_mets_AK5_ex;   //!
   TBranch        *b_mets_AK5_ey;   //!
   TBranch        *b_mets_AK5_gen_et;   //!
   TBranch        *b_mets_AK5_gen_phi;   //!
   TBranch        *b_mets_AK5_sign;   //!
   TBranch        *b_mets_AK5_sumEt;   //!
   TBranch        *b_mets_AK5_unCPhi;   //!
   TBranch        *b_mets_AK5_unCPt;   //!
   TBranch        *b_Nmets_SC5;   //!
   TBranch        *b_mets_SC5_et;   //!
   TBranch        *b_mets_SC5_phi;   //!
   TBranch        *b_mets_SC5_ex;   //!
   TBranch        *b_mets_SC5_ey;   //!
   TBranch        *b_mets_SC5_gen_et;   //!
   TBranch        *b_mets_SC5_gen_phi;   //!
   TBranch        *b_mets_SC5_sign;   //!
   TBranch        *b_mets_SC5_sumEt;   //!
   TBranch        *b_mets_SC5_unCPhi;   //!
   TBranch        *b_mets_SC5_unCPt;   //!
   TBranch        *b_Nmus;   //!
   TBranch        *b_mus_energy;   //!
   TBranch        *b_mus_et;   //!
   TBranch        *b_mus_eta;   //!
   TBranch        *b_mus_phi;   //!
   TBranch        *b_mus_pt;   //!
   TBranch        *b_mus_px;   //!
   TBranch        *b_mus_py;   //!
   TBranch        *b_mus_pz;   //!
   TBranch        *b_mus_status;   //!
   TBranch        *b_mus_theta;   //!
   TBranch        *b_mus_gen_id;   //!
   TBranch        *b_mus_gen_phi;   //!
   TBranch        *b_mus_gen_pt;   //!
   TBranch        *b_mus_gen_pz;   //!
   TBranch        *b_mus_gen_px;   //!
   TBranch        *b_mus_gen_py;   //!
   TBranch        *b_mus_gen_eta;   //!
   TBranch        *b_mus_gen_theta;   //!
   TBranch        *b_mus_gen_et;   //!
   TBranch        *b_mus_gen_mother_id;   //!
   TBranch        *b_mus_gen_mother_phi;   //!
   TBranch        *b_mus_gen_mother_pt;   //!
   TBranch        *b_mus_gen_mother_pz;   //!
   TBranch        *b_mus_gen_mother_px;   //!
   TBranch        *b_mus_gen_mother_py;   //!
   TBranch        *b_mus_gen_mother_eta;   //!
   TBranch        *b_mus_gen_mother_theta;   //!
   TBranch        *b_mus_gen_mother_et;   //!
   TBranch        *b_mus_tkHits;   //!
   TBranch        *b_mus_cIso;   //!
   TBranch        *b_mus_tIso;   //!
   TBranch        *b_mus_ecalIso;   //!
   TBranch        *b_mus_hcalIso;   //!
   TBranch        *b_mus_ecalvetoDep;   //!
   TBranch        *b_mus_hcalvetoDep;   //!
   TBranch        *b_mus_calEnergyEm;   //!
   TBranch        *b_mus_calEnergyHad;   //!
   TBranch        *b_mus_calEnergyHo;   //!
   TBranch        *b_mus_calEnergyEmS9;   //!
   TBranch        *b_mus_calEnergyHadS9;   //!
   TBranch        *b_mus_calEnergyHoS9;   //!
   TBranch        *b_mus_iso03_sumPt;   //!
   TBranch        *b_mus_iso03_emEt;   //!
   TBranch        *b_mus_iso03_hadEt;   //!
   TBranch        *b_mus_iso03_hoEt;   //!
   TBranch        *b_mus_iso03_nTracks;   //!
   TBranch        *b_mus_iso05_sumPt;   //!
   TBranch        *b_mus_iso05_emEt;   //!
   TBranch        *b_mus_iso05_hadEt;   //!
   TBranch        *b_mus_iso05_hoEt;   //!
   TBranch        *b_mus_iso05_nTracks;   //!
   TBranch        *b_mus_charge;   //!
   TBranch        *b_mus_cm_chi2;   //!
   TBranch        *b_mus_cm_ndof;   //!
   TBranch        *b_mus_cm_chg;   //!
   TBranch        *b_mus_cm_pt;   //!
   TBranch        *b_mus_cm_px;   //!
   TBranch        *b_mus_cm_py;   //!
   TBranch        *b_mus_cm_pz;   //!
   TBranch        *b_mus_cm_eta;   //!
   TBranch        *b_mus_cm_phi;   //!
   TBranch        *b_mus_cm_theta;   //!
   TBranch        *b_mus_cm_d0dum;   //!
   TBranch        *b_mus_cm_dz;   //!
   TBranch        *b_mus_cm_vx;   //!
   TBranch        *b_mus_cm_vy;   //!
   TBranch        *b_mus_cm_vz;   //!
   TBranch        *b_mus_cm_numvalhits;   //!
   TBranch        *b_mus_cm_numlosthits;   //!
   TBranch        *b_mus_cm_d0dumErr;   //!
   TBranch        *b_mus_cm_dzErr;   //!
   TBranch        *b_mus_cm_ptErr;   //!
   TBranch        *b_mus_cm_etaErr;   //!
   TBranch        *b_mus_cm_phiErr;   //!
   TBranch        *b_mus_tk_id;   //!
   TBranch        *b_mus_tk_chi2;   //!
   TBranch        *b_mus_tk_ndof;   //!
   TBranch        *b_mus_tk_chg;   //!
   TBranch        *b_mus_tk_pt;   //!
   TBranch        *b_mus_tk_px;   //!
   TBranch        *b_mus_tk_py;   //!
   TBranch        *b_mus_tk_pz;   //!
   TBranch        *b_mus_tk_eta;   //!
   TBranch        *b_mus_tk_phi;   //!
   TBranch        *b_mus_tk_theta;   //!
   TBranch        *b_mus_tk_d0dum;   //!
   TBranch        *b_mus_tk_dz;   //!
   TBranch        *b_mus_tk_vx;   //!
   TBranch        *b_mus_tk_vy;   //!
   TBranch        *b_mus_tk_vz;   //!
   TBranch        *b_mus_tk_numvalhits;   //!
   TBranch        *b_mus_tk_numlosthits;   //!
   TBranch        *b_mus_tk_d0dumErr;   //!
   TBranch        *b_mus_tk_dzErr;   //!
   TBranch        *b_mus_tk_ptErr;   //!
   TBranch        *b_mus_tk_etaErr;   //!
   TBranch        *b_mus_tk_phiErr;   //!
   TBranch        *b_mus_stamu_chi2;   //!
   TBranch        *b_mus_stamu_ndof;   //!
   TBranch        *b_mus_stamu_chg;   //!
   TBranch        *b_mus_stamu_pt;   //!
   TBranch        *b_mus_stamu_px;   //!
   TBranch        *b_mus_stamu_py;   //!
   TBranch        *b_mus_stamu_pz;   //!
   TBranch        *b_mus_stamu_eta;   //!
   TBranch        *b_mus_stamu_phi;   //!
   TBranch        *b_mus_stamu_theta;   //!
   TBranch        *b_mus_stamu_d0dum;   //!
   TBranch        *b_mus_stamu_dz;   //!
   TBranch        *b_mus_stamu_vx;   //!
   TBranch        *b_mus_stamu_vy;   //!
   TBranch        *b_mus_stamu_vz;   //!
   TBranch        *b_mus_stamu_numvalhits;   //!
   TBranch        *b_mus_stamu_numlosthits;   //!
   TBranch        *b_mus_stamu_d0dumErr;   //!
   TBranch        *b_mus_stamu_dzErr;   //!
   TBranch        *b_mus_stamu_ptErr;   //!
   TBranch        *b_mus_stamu_etaErr;   //!
   TBranch        *b_mus_stamu_phiErr;   //!
   TBranch        *b_mus_num_matches;   //!
   TBranch        *b_mus_id_All;   //!
   TBranch        *b_mus_id_AllGlobalMuons;   //!
   TBranch        *b_mus_id_AllStandAloneMuons;   //!
   TBranch        *b_mus_id_AllTrackerMuons;   //!
   TBranch        *b_mus_id_TrackerMuonArbitrated;   //!
   TBranch        *b_mus_id_AllArbitrated;   //!
   TBranch        *b_mus_id_GlobalMuonPromptTight;   //!
   TBranch        *b_mus_id_TMLastStationLoose;   //!
   TBranch        *b_mus_id_TMLastStationTight;   //!
   TBranch        *b_mus_id_TM2DCompatibilityLoose;   //!
   TBranch        *b_mus_id_TM2DCompatibilityTight;   //!
   TBranch        *b_mus_id_TMOneStationLoose;   //!
   TBranch        *b_mus_id_TMOneStationTight;   //!
   TBranch        *b_mus_id_TMLastStationOptimizedLowPtLoose;   //!
   TBranch        *b_mus_id_TMLastStationOptimizedLowPtTight;   //!
   TBranch        *b_Nphotons;   //!
   TBranch        *b_photons_energy;   //!
   TBranch        *b_photons_et;   //!
   TBranch        *b_photons_eta;   //!
   TBranch        *b_photons_phi;   //!
   TBranch        *b_photons_pt;   //!
   TBranch        *b_photons_px;   //!
   TBranch        *b_photons_py;   //!
   TBranch        *b_photons_pz;   //!
   TBranch        *b_photons_status;   //!
   TBranch        *b_photons_theta;   //!
   TBranch        *b_photons_hadOverEM;   //!
   TBranch        *b_photons_scEnergy;   //!
   TBranch        *b_photons_scRawEnergy;   //!
   TBranch        *b_photons_scEta;   //!
   TBranch        *b_photons_scPhi;   //!
   TBranch        *b_photons_scEtaWidth;   //!
   TBranch        *b_photons_scPhiWidth;   //!
   TBranch        *b_photons_tIso;   //!
   TBranch        *b_photons_ecalIso;   //!
   TBranch        *b_photons_hcalIso;   //!
   TBranch        *b_photons_isoEcalRecHitDR04;   //!
   TBranch        *b_photons_isoHcalRecHitDR04;   //!
   TBranch        *b_photons_isoSolidTrkConeDR04;   //!
   TBranch        *b_photons_isoHollowTrkConeDR04;   //!
   TBranch        *b_photons_nTrkSolidConeDR04;   //!
   TBranch        *b_photons_nTrkHollowConeDR04;   //!
   TBranch        *b_photons_isoEcalRecHitDR03;   //!
   TBranch        *b_photons_isoHcalRecHitDR03;   //!
   TBranch        *b_photons_isoSolidTrkConeDR03;   //!
   TBranch        *b_photons_isoHollowTrkConeDR03;   //!
   TBranch        *b_photons_nTrkSolidConeDR03;   //!
   TBranch        *b_photons_nTrkHollowConeDR03;   //!
   TBranch        *b_photons_isAlsoElectron;   //!
   TBranch        *b_photons_hasPixelSeed;   //!
   TBranch        *b_photons_isConverted;   //!
   TBranch        *b_photons_isEBGap;   //!
   TBranch        *b_photons_isEEGap;   //!
   TBranch        *b_photons_isEBEEGap;   //!
   TBranch        *b_photons_isEBPho;   //!
   TBranch        *b_photons_isEEPho;   //!
   TBranch        *b_photons_isLoosePhoton;   //!
   TBranch        *b_photons_isTightPhoton;   //!
   TBranch        *b_photons_r9;   //!
   TBranch        *b_photons_gen_et;   //!
   TBranch        *b_photons_gen_eta;   //!
   TBranch        *b_photons_gen_phi;   //!
   TBranch        *b_photons_gen_id;   //!
   TBranch        *b_Npv;   //!
   TBranch        *b_pv_x;   //!
   TBranch        *b_pv_y;   //!
   TBranch        *b_pv_z;   //!
   TBranch        *b_pv_xErr;   //!
   TBranch        *b_pv_yErr;   //!
   TBranch        *b_pv_zErr;   //!
   TBranch        *b_Ntaus;   //!
   TBranch        *b_taus_energy;   //!
   TBranch        *b_taus_et;   //!
   TBranch        *b_taus_eta;   //!
   TBranch        *b_taus_phi;   //!
   TBranch        *b_taus_pt;   //!
   TBranch        *b_taus_px;   //!
   TBranch        *b_taus_py;   //!
   TBranch        *b_taus_pz;   //!
   TBranch        *b_taus_status;   //!
   TBranch        *b_taus_theta;   //!
   TBranch        *b_taus_charge;   //!
   TBranch        *b_taus_emf;   //!
   TBranch        *b_taus_hcalTotOverPLead;   //!
   TBranch        *b_taus_hcalMaxOverPLead;   //!
   TBranch        *b_taus_hcal3x3OverPLead;   //!
   TBranch        *b_taus_ecalStripSumEOverPLead;   //!
   TBranch        *b_taus_elecPreIdOutput;   //!
   TBranch        *b_taus_elecPreIdDecision;   //!
   TBranch        *b_taus_leadPFChargedHadrCand_pt;   //!
   TBranch        *b_taus_leadPFChargedHadrCand_charge;   //!
   TBranch        *b_taus_leadPFChargedHadrCand_eta;   //!
   TBranch        *b_taus_leadPFChargedHadrCand_ECAL_eta;   //!
   TBranch        *b_taus_leadPFChargedHadrCand_phi;   //!
   TBranch        *b_taus_isoPFGammaCandsEtSum;   //!
   TBranch        *b_taus_isoPFChargedHadrCandsPtSum;   //!
   TBranch        *b_taus_leadingTrackFinding;   //!
   TBranch        *b_taus_leadingTrackPtCut;   //!
   TBranch        *b_taus_trackIsolation;   //!
   TBranch        *b_taus_ecalIsolation;   //!
   TBranch        *b_taus_byIsolation;   //!
   TBranch        *b_taus_againstElectron;   //!
   TBranch        *b_taus_againstMuon;   //!
   TBranch        *b_taus_taNC_quarter;   //!
   TBranch        *b_taus_taNC_one;   //!
   TBranch        *b_taus_taNC_half;   //!
   TBranch        *b_taus_taNC_tenth;   //!
   TBranch        *b_taus_muDecision;   //!
   TBranch        *b_taus_Nprongs;   //!
   TBranch        *b_Ntcmets;   //!
   TBranch        *b_tcmets_et;   //!
   TBranch        *b_tcmets_phi;   //!
   TBranch        *b_tcmets_ex;   //!
   TBranch        *b_tcmets_ey;   //!
   TBranch        *b_tcmets_sumEt;   //!
   TBranch        *b_Ntracks;   //!
   TBranch        *b_tracks_chi2;   //!
   TBranch        *b_tracks_ndof;   //!
   TBranch        *b_tracks_chg;   //!
   TBranch        *b_tracks_pt;   //!
   TBranch        *b_tracks_px;   //!
   TBranch        *b_tracks_py;   //!
   TBranch        *b_tracks_pz;   //!
   TBranch        *b_tracks_eta;   //!
   TBranch        *b_tracks_phi;   //!
   TBranch        *b_tracks_theta;   //!
   TBranch        *b_tracks_d0dum;   //!
   TBranch        *b_tracks_dz;   //!
   TBranch        *b_tracks_vx;   //!
   TBranch        *b_tracks_vy;   //!
   TBranch        *b_tracks_vz;   //!
   TBranch        *b_tracks_numvalhits;   //!
   TBranch        *b_tracks_numlosthits;   //!
   TBranch        *b_tracks_d0dumErr;   //!
   TBranch        *b_tracks_dzErr;   //!
   TBranch        *b_tracks_ptErr;   //!
   TBranch        *b_tracks_etaErr;   //!
   TBranch        *b_tracks_phiErr;   //!
   TBranch        *b_tracks_Nrechits;   //!
   TBranch        *b_tracks_innerHitX;   //!
   TBranch        *b_tracks_innerHitY;   //!
   TBranch        *b_tracks_innerHitZ;   //!
   TBranch        *b_tracks_outerHitX;   //!
   TBranch        *b_tracks_outerHitY;   //!
   TBranch        *b_tracks_outerHitZ;   //!
   TBranch        *b_tracks_highPurity;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!

   // Declaration of leaf types
   Double_t        AlCa_EcalEta_8E29;
   Double_t        AlCa_EcalPhiSym;
   Double_t        AlCa_EcalPi0_8E29;
   Double_t        AlCa_HcalPhiSym;
   Double_t        AlCa_RPCMuonNoHits;
   Double_t        AlCa_RPCMuonNormalisation;
   Double_t        HLTAnalyzerEndpath;
   Double_t        HLT_BTagIP_Jet50U;
   Double_t        HLT_BTagMu_Jet10U;
   Double_t        HLT_BackwardBSC;
   Double_t        HLT_CSCBeamHalo;
   Double_t        HLT_CSCBeamHaloOverlapRing1;
   Double_t        HLT_CSCBeamHaloOverlapRing2;
   Double_t        HLT_CSCBeamHaloRing2or3;
   Double_t        HLT_DiJetAve15U_8E29;
   Double_t        HLT_DiJetAve30U_8E29;
   Double_t        HLT_DoubleEle5_SW_L1R;
   Double_t        HLT_DoubleLooseIsoTau15;
   Double_t        HLT_DoubleMu0;
   Double_t        HLT_DoubleMu3;
   Double_t        HLT_DoublePhoton10_L1R;
   Double_t        HLT_DoublePhoton5_Jpsi_L1R;
   Double_t        HLT_DoublePhoton5_Upsilon_L1R;
   Double_t        HLT_DoublePhoton5_eeRes_L1R;
   Double_t        HLT_Ele10_LW_EleId_L1R;
   Double_t        HLT_Ele10_LW_L1R;
   Double_t        HLT_Ele15_LW_L1R;
   Double_t        HLT_Ele15_SC10_LW_L1R;
   Double_t        HLT_Ele15_SiStrip_L1R;
   Double_t        HLT_Ele20_LW_L1R;
   Double_t        HLT_ForwardBSC;
   Double_t        HLT_FwdJet20U;
   Double_t        HLT_HT100U;
   Double_t        HLT_IsoMu3;
   Double_t        HLT_IsoTrack_8E29;
   Double_t        HLT_Jet15U;
   Double_t        HLT_Jet30U;
   Double_t        HLT_Jet50U;
   Double_t        HLT_L1DoubleEG5;
   Double_t        HLT_L1DoubleMuOpen;
   Double_t        HLT_L1Jet6U;
   Double_t        HLT_L1MET20;
   Double_t        HLT_L1Mu;
   Double_t        HLT_L1Mu14_L1ETM30;
   Double_t        HLT_L1Mu14_L1SingleEG10;
   Double_t        HLT_L1Mu14_L1SingleJet6U;
   Double_t        HLT_L1Mu20;
   Double_t        HLT_L1MuOpen;
   Double_t        HLT_L1SingleEG5;
   Double_t        HLT_L1SingleEG8;
   Double_t        HLT_L2Mu11;
   Double_t        HLT_L2Mu9;
   Double_t        HLT_MET100;
   Double_t        HLT_MET45;
   Double_t        HLT_MinBiasEcal;
   Double_t        HLT_MinBiasHcal;
   Double_t        HLT_MinBiasPixel;
   Double_t        HLT_MinBiasPixel_Trk5;
   Double_t        HLT_Mu3;
   Double_t        HLT_Mu5;
   Double_t        HLT_Mu9;
   Double_t        HLT_Photon10_L1R;
   Double_t        HLT_Photon15_L1R;
   Double_t        HLT_Photon15_LooseEcalIso_L1R;
   Double_t        HLT_Photon15_TrackIso_L1R;
   Double_t        HLT_Photon20_L1R;
   Double_t        HLT_Photon30_L1R_8E29;
   Double_t        HLT_QuadJet15U;
   Double_t        HLT_SingleLooseIsoTau20;
   Double_t        HLT_StoppedHSCP_8E29;
   Double_t        HLT_TrackerCosmics;
   Double_t        HLT_ZeroBias;
   Double_t        HLTriggerFinalPath;
   Double_t        HLTriggerFirstPath;

   // List of branches
   TBranch        *b_AlCa_EcalEta_8E29;   //!
   TBranch        *b_AlCa_EcalPhiSym;   //!
   TBranch        *b_AlCa_EcalPi0_8E29;   //!
   TBranch        *b_AlCa_HcalPhiSym;   //!
   TBranch        *b_AlCa_RPCMuonNoHits;   //!
   TBranch        *b_AlCa_RPCMuonNormalisation;   //!
   TBranch        *b_HLTAnalyzerEndpath;   //!
   TBranch        *b_HLT_BTagIP_Jet50U;   //!
   TBranch        *b_HLT_BTagMu_Jet10U;   //!
   TBranch        *b_HLT_BackwardBSC;   //!
   TBranch        *b_HLT_CSCBeamHalo;   //!
   TBranch        *b_HLT_CSCBeamHaloOverlapRing1;   //!
   TBranch        *b_HLT_CSCBeamHaloOverlapRing2;   //!
   TBranch        *b_HLT_CSCBeamHaloRing2or3;   //!
   TBranch        *b_HLT_DiJetAve15U_8E29;   //!
   TBranch        *b_HLT_DiJetAve30U_8E29;   //!
   TBranch        *b_HLT_DoubleEle5_SW_L1R;   //!
   TBranch        *b_HLT_DoubleLooseIsoTau15;   //!
   TBranch        *b_HLT_DoubleMu0;   //!
   TBranch        *b_HLT_DoubleMu3;   //!
   TBranch        *b_HLT_DoublePhoton10_L1R;   //!
   TBranch        *b_HLT_DoublePhoton5_Jpsi_L1R;   //!
   TBranch        *b_HLT_DoublePhoton5_Upsilon_L1R;   //!
   TBranch        *b_HLT_DoublePhoton5_eeRes_L1R;   //!
   TBranch        *b_HLT_Ele10_LW_EleId_L1R;   //!
   TBranch        *b_HLT_Ele10_LW_L1R;   //!
   TBranch        *b_HLT_Ele15_LW_L1R;   //!
   TBranch        *b_HLT_Ele15_SC10_LW_L1R;   //!
   TBranch        *b_HLT_Ele15_SiStrip_L1R;   //!
   TBranch        *b_HLT_Ele20_LW_L1R;   //!
   TBranch        *b_HLT_ForwardBSC;   //!
   TBranch        *b_HLT_FwdJet20U;   //!
   TBranch        *b_HLT_HT100U;   //!
   TBranch        *b_HLT_IsoMu3;   //!
   TBranch        *b_HLT_IsoTrack_8E29;   //!
   TBranch        *b_HLT_Jet15U;   //!
   TBranch        *b_HLT_Jet30U;   //!
   TBranch        *b_HLT_Jet50U;   //!
   TBranch        *b_HLT_L1DoubleEG5;   //!
   TBranch        *b_HLT_L1DoubleMuOpen;   //!
   TBranch        *b_HLT_L1Jet6U;   //!
   TBranch        *b_HLT_L1MET20;   //!
   TBranch        *b_HLT_L1Mu;   //!
   TBranch        *b_HLT_L1Mu14_L1ETM30;   //!
   TBranch        *b_HLT_L1Mu14_L1SingleEG10;   //!
   TBranch        *b_HLT_L1Mu14_L1SingleJet6U;   //!
   TBranch        *b_HLT_L1Mu20;   //!
   TBranch        *b_HLT_L1MuOpen;   //!
   TBranch        *b_HLT_L1SingleEG5;   //!
   TBranch        *b_HLT_L1SingleEG8;   //!
   TBranch        *b_HLT_L2Mu11;   //!
   TBranch        *b_HLT_L2Mu9;   //!
   TBranch        *b_HLT_MET100;   //!
   TBranch        *b_HLT_MET45;   //!
   TBranch        *b_HLT_MinBiasEcal;   //!
   TBranch        *b_HLT_MinBiasHcal;   //!
   TBranch        *b_HLT_MinBiasPixel;   //!
   TBranch        *b_HLT_MinBiasPixel_Trk5;   //!
   TBranch        *b_HLT_Mu3;   //!
   TBranch        *b_HLT_Mu5;   //!
   TBranch        *b_HLT_Mu9;   //!
   TBranch        *b_HLT_Photon10_L1R;   //!
   TBranch        *b_HLT_Photon15_L1R;   //!
   TBranch        *b_HLT_Photon15_LooseEcalIso_L1R;   //!
   TBranch        *b_HLT_Photon15_TrackIso_L1R;   //!
   TBranch        *b_HLT_Photon20_L1R;   //!
   TBranch        *b_HLT_Photon30_L1R_8E29;   //!
   TBranch        *b_HLT_QuadJet15U;   //!
   TBranch        *b_HLT_SingleLooseIsoTau20;   //!
   TBranch        *b_HLT_StoppedHSCP_8E29;   //!
   TBranch        *b_HLT_TrackerCosmics;   //!
   TBranch        *b_HLT_ZeroBias;   //!
   TBranch        *b_HLTriggerFinalPath;   //!
   TBranch        *b_HLTriggerFirstPath;   //!

//void Initialize(TTree *t1,TTree *t2)
void InitializeB(TChain *fChain)
{

   // Set object pointer
   beamSpot_x = 0;
   beamSpot_y = 0;
   beamSpot_z = 0;
   beamSpot_x0Error = 0;
   beamSpot_y0Error = 0;
   beamSpot_z0Error = 0;
   beamSpot_sigmaZ = 0;
   beamSpot_sigmaZ0Error = 0;
   beamSpot_dxdz = 0;
   beamSpot_dxdzError = 0;
   beamSpot_dydz = 0;
   beamSpot_dydzError = 0;
   beamSpot_beamWidthX = 0;
   beamSpot_beamWidthY = 0;
   beamSpot_beamWidthXError = 0;
   beamSpot_beamWidthYError = 0;
   els_energy = 0;
   els_et = 0;
   els_eta = 0;
   els_phi = 0;
   els_pt = 0;
   els_px = 0;
   els_py = 0;
   els_pz = 0;
   els_status = 0;
   els_theta = 0;
   els_gen_id = 0;
   els_gen_phi = 0;
   els_gen_pt = 0;
   els_gen_pz = 0;
   els_gen_px = 0;
   els_gen_py = 0;
   els_gen_eta = 0;
   els_gen_theta = 0;
   els_gen_et = 0;
   els_gen_mother_id = 0;
   els_gen_mother_phi = 0;
   els_gen_mother_pt = 0;
   els_gen_mother_pz = 0;
   els_gen_mother_px = 0;
   els_gen_mother_py = 0;
   els_gen_mother_eta = 0;
   els_gen_mother_theta = 0;
   els_gen_mother_et = 0;
   els_tightId = 0;
   els_looseId = 0;
   els_robustTightId = 0;
   els_robustLooseId = 0;
   els_robustHighEnergyId = 0;
   els_cIso = 0;
   els_tIso = 0;
   els_ecalIso = 0;
   els_hcalIso = 0;
   els_chi2 = 0;
   els_charge = 0;
   els_caloEnergy = 0;
   els_hadOverEm = 0;
   els_eOverPIn = 0;
   els_eSeedOverPOut = 0;
   els_eSCraw = 0;
   els_eSeed = 0;
   els_sigmaEtaEta = 0;
   els_sigmaIEtaIEta = 0;
   els_scE1x5 = 0;
   els_scE2x5Max = 0;
   els_scE5x5 = 0;
   els_dEtaIn = 0;
   els_dPhiIn = 0;
   els_dEtaOut = 0;
   els_dPhiOut = 0;
   els_numvalhits = 0;
   els_numlosthits = 0;
   els_basicClustersSize = 0;
   els_tk_pt = 0;
   els_tk_phi = 0;
   els_tk_eta = 0;
   els_d0dum = 0;
   els_dz = 0;
   els_vx = 0;
   els_vy = 0;
   els_vz = 0;
   els_ndof = 0;
   els_ptError = 0;
   els_d0dumError = 0;
   els_dzError = 0;
   els_etaError = 0;
   els_phiError = 0;
   els_tk_charge = 0;
   els_ctf_tk_id = 0;
   els_ctf_tk_charge = 0;
   els_ctf_tk_eta = 0;
   els_ctf_tk_phi = 0;
   els_fbrem = 0;
   els_shFracInnerHits = 0;
   els_dr03EcalRecHitSumEt = 0;
   els_dr03HcalTowerSumEt = 0;
   els_dr03HcalDepth1TowerSumEt = 0;
   els_dr03HcalDepth2TowerSumEt = 0;
   els_dr03TkSumPt = 0;
   els_dr04EcalRecHitSumEt = 0;
   els_dr04HcalTowerSumEt = 0;
   els_dr04HcalDepth1TowerSumEt = 0;
   els_dr04HcalDepth2TowerSumEt = 0;
   els_dr04TkSumPt = 0;
   els_cpx = 0;
   els_cpy = 0;
   els_cpz = 0;
   els_vpx = 0;
   els_vpy = 0;
   els_vpz = 0;
   els_cx = 0;
   els_cy = 0;
   els_cz = 0;
   jets_AK5_energy = 0;
   jets_AK5_et = 0;
   jets_AK5_eta = 0;
   jets_AK5_phi = 0;
   jets_AK5_pt = 0;
   jets_AK5_px = 0;
   jets_AK5_py = 0;
   jets_AK5_pz = 0;
   jets_AK5_status = 0;
   jets_AK5_theta = 0;
   jets_AK5_parton_Id = 0;
   jets_AK5_parton_motherId = 0;
   jets_AK5_parton_pt = 0;
   jets_AK5_parton_phi = 0;
   jets_AK5_parton_eta = 0;
   jets_AK5_parton_Energy = 0;
   jets_AK5_parton_mass = 0;
   jets_AK5_parton_motherID = 0;
   jets_AK5_gen_et = 0;
   jets_AK5_gen_pt = 0;
   jets_AK5_gen_eta = 0;
   jets_AK5_gen_phi = 0;
   jets_AK5_gen_mass = 0;
   jets_AK5_gen_Energy = 0;
   jets_AK5_gen_Id = 0;
   jets_AK5_gen_motherID = 0;
   jets_AK5_gen_threeCharge = 0;
   jets_AK5_partonFlavour = 0;
   jets_AK5_btag_TC_highPur = 0;
   jets_AK5_btag_TC_highEff = 0;
   jets_AK5_btag_jetProb = 0;
   jets_AK5_btag_jetBProb = 0;
   jets_AK5_btag_softEle = 0;
   jets_AK5_btag_softMuon = 0;
   jets_AK5_btag_softMuonNoIP = 0;
   jets_AK5_btag_secVertex = 0;
   jets_AK5_chgEmE = 0;
   jets_AK5_chgHadE = 0;
   jets_AK5_chgMuE = 0;
   jets_AK5_chg_Mult = 0;
   jets_AK5_neutralEmE = 0;
   jets_AK5_neutralHadE = 0;
   jets_AK5_neutral_Mult = 0;
   jets_AK5_mu_Mult = 0;
   jets_AK5_emf = 0;
   jets_AK5_ehf = 0;
   jets_AK5_n60 = 0;
   jets_AK5_n90 = 0;
   jets_AK5_area = 0;
   jets_AK5_mass = 0;
   jets_SC5_energy = 0;
   jets_SC5_et = 0;
   jets_SC5_eta = 0;
   jets_SC5_phi = 0;
   jets_SC5_pt = 0;
   jets_SC5_px = 0;
   jets_SC5_py = 0;
   jets_SC5_pz = 0;
   jets_SC5_status = 0;
   jets_SC5_theta = 0;
   jets_SC5_parton_Id = 0;
   jets_SC5_parton_motherId = 0;
   jets_SC5_parton_pt = 0;
   jets_SC5_parton_phi = 0;
   jets_SC5_parton_eta = 0;
   jets_SC5_parton_Energy = 0;
   jets_SC5_parton_mass = 0;
   jets_SC5_parton_motherID = 0;
   jets_SC5_gen_et = 0;
   jets_SC5_gen_pt = 0;
   jets_SC5_gen_eta = 0;
   jets_SC5_gen_phi = 0;
   jets_SC5_gen_mass = 0;
   jets_SC5_gen_Energy = 0;
   jets_SC5_gen_Id = 0;
   jets_SC5_gen_motherID = 0;
   jets_SC5_gen_threeCharge = 0;
   jets_SC5_partonFlavour = 0;
   jets_SC5_btag_TC_highPur = 0;
   jets_SC5_btag_TC_highEff = 0;
   jets_SC5_btag_jetProb = 0;
   jets_SC5_btag_jetBProb = 0;
   jets_SC5_btag_softEle = 0;
   jets_SC5_btag_softMuon = 0;
   jets_SC5_btag_softMuonNoIP = 0;
   jets_SC5_btag_secVertex = 0;
   jets_SC5_chgEmE = 0;
   jets_SC5_chgHadE = 0;
   jets_SC5_chgMuE = 0;
   jets_SC5_chg_Mult = 0;
   jets_SC5_neutralEmE = 0;
   jets_SC5_neutralHadE = 0;
   jets_SC5_neutral_Mult = 0;
   jets_SC5_mu_Mult = 0;
   jets_SC5_emf = 0;
   jets_SC5_ehf = 0;
   jets_SC5_n60 = 0;
   jets_SC5_n90 = 0;
   jets_SC5_area = 0;
   jets_SC5_mass = 0;
   mc_doc_id = 0;
   mc_doc_pt = 0;
   mc_doc_px = 0;
   mc_doc_py = 0;
   mc_doc_pz = 0;
   mc_doc_eta = 0;
   mc_doc_phi = 0;
   mc_doc_theta = 0;
   mc_doc_energy = 0;
   mc_doc_status = 0;
   mc_doc_charge = 0;
   mc_doc_mother_id = 0;
   mc_doc_grandmother_id = 0;
   mc_doc_ggrandmother_id = 0;
   mc_doc_mother_pt = 0;
   mc_doc_vertex_x = 0;
   mc_doc_vertex_y = 0;
   mc_doc_vertex_z = 0;
   mc_doc_mass = 0;
   mc_doc_numOfDaughters = 0;
   mc_doc_numOfMothers = 0;
   mc_electrons_id = 0;
   mc_electrons_pt = 0;
   mc_electrons_px = 0;
   mc_electrons_py = 0;
   mc_electrons_pz = 0;
   mc_electrons_eta = 0;
   mc_electrons_phi = 0;
   mc_electrons_theta = 0;
   mc_electrons_status = 0;
   mc_electrons_energy = 0;
   mc_electrons_charge = 0;
   mc_electrons_mother_id = 0;
   mc_electrons_mother_pt = 0;
   mc_electrons_grandmother_id = 0;
   mc_electrons_ggrandmother_id = 0;
   mc_electrons_vertex_x = 0;
   mc_electrons_vertex_y = 0;
   mc_electrons_vertex_z = 0;
   mc_electrons_mass = 0;
   mc_electrons_numOfDaughters = 0;
   mc_mus_id = 0;
   mc_mus_pt = 0;
   mc_mus_px = 0;
   mc_mus_py = 0;
   mc_mus_pz = 0;
   mc_mus_eta = 0;
   mc_mus_phi = 0;
   mc_mus_theta = 0;
   mc_mus_status = 0;
   mc_mus_energy = 0;
   mc_mus_charge = 0;
   mc_mus_mother_id = 0;
   mc_mus_mother_pt = 0;
   mc_mus_grandmother_id = 0;
   mc_mus_ggrandmother_id = 0;
   mc_mus_vertex_x = 0;
   mc_mus_vertex_y = 0;
   mc_mus_vertex_z = 0;
   mc_mus_mass = 0;
   mc_mus_numOfDaughters = 0;
   mets_AK5_et = 0;
   mets_AK5_phi = 0;
   mets_AK5_ex = 0;
   mets_AK5_ey = 0;
   mets_AK5_gen_et = 0;
   mets_AK5_gen_phi = 0;
   mets_AK5_sign = 0;
   mets_AK5_sumEt = 0;
   mets_AK5_unCPhi = 0;
   mets_AK5_unCPt = 0;
   mets_SC5_et = 0;
   mets_SC5_phi = 0;
   mets_SC5_ex = 0;
   mets_SC5_ey = 0;
   mets_SC5_gen_et = 0;
   mets_SC5_gen_phi = 0;
   mets_SC5_sign = 0;
   mets_SC5_sumEt = 0;
   mets_SC5_unCPhi = 0;
   mets_SC5_unCPt = 0;
   mus_energy = 0;
   mus_et = 0;
   mus_eta = 0;
   mus_phi = 0;
   mus_pt = 0;
   mus_px = 0;
   mus_py = 0;
   mus_pz = 0;
   mus_status = 0;
   mus_theta = 0;
   mus_gen_id = 0;
   mus_gen_phi = 0;
   mus_gen_pt = 0;
   mus_gen_pz = 0;
   mus_gen_px = 0;
   mus_gen_py = 0;
   mus_gen_eta = 0;
   mus_gen_theta = 0;
   mus_gen_et = 0;
   mus_gen_mother_id = 0;
   mus_gen_mother_phi = 0;
   mus_gen_mother_pt = 0;
   mus_gen_mother_pz = 0;
   mus_gen_mother_px = 0;
   mus_gen_mother_py = 0;
   mus_gen_mother_eta = 0;
   mus_gen_mother_theta = 0;
   mus_gen_mother_et = 0;
   mus_tkHits = 0;
   mus_cIso = 0;
   mus_tIso = 0;
   mus_ecalIso = 0;
   mus_hcalIso = 0;
   mus_ecalvetoDep = 0;
   mus_hcalvetoDep = 0;
   mus_calEnergyEm = 0;
   mus_calEnergyHad = 0;
   mus_calEnergyHo = 0;
   mus_calEnergyEmS9 = 0;
   mus_calEnergyHadS9 = 0;
   mus_calEnergyHoS9 = 0;
   mus_iso03_sumPt = 0;
   mus_iso03_emEt = 0;
   mus_iso03_hadEt = 0;
   mus_iso03_hoEt = 0;
   mus_iso03_nTracks = 0;
   mus_iso05_sumPt = 0;
   mus_iso05_emEt = 0;
   mus_iso05_hadEt = 0;
   mus_iso05_hoEt = 0;
   mus_iso05_nTracks = 0;
   mus_charge = 0;
   mus_cm_chi2 = 0;
   mus_cm_ndof = 0;
   mus_cm_chg = 0;
   mus_cm_pt = 0;
   mus_cm_px = 0;
   mus_cm_py = 0;
   mus_cm_pz = 0;
   mus_cm_eta = 0;
   mus_cm_phi = 0;
   mus_cm_theta = 0;
   mus_cm_d0dum = 0;
   mus_cm_dz = 0;
   mus_cm_vx = 0;
   mus_cm_vy = 0;
   mus_cm_vz = 0;
   mus_cm_numvalhits = 0;
   mus_cm_numlosthits = 0;
   mus_cm_d0dumErr = 0;
   mus_cm_dzErr = 0;
   mus_cm_ptErr = 0;
   mus_cm_etaErr = 0;
   mus_cm_phiErr = 0;
   mus_tk_id = 0;
   mus_tk_chi2 = 0;
   mus_tk_ndof = 0;
   mus_tk_chg = 0;
   mus_tk_pt = 0;
   mus_tk_px = 0;
   mus_tk_py = 0;
   mus_tk_pz = 0;
   mus_tk_eta = 0;
   mus_tk_phi = 0;
   mus_tk_theta = 0;
   mus_tk_d0dum = 0;
   mus_tk_dz = 0;
   mus_tk_vx = 0;
   mus_tk_vy = 0;
   mus_tk_vz = 0;
   mus_tk_numvalhits = 0;
   mus_tk_numlosthits = 0;
   mus_tk_d0dumErr = 0;
   mus_tk_dzErr = 0;
   mus_tk_ptErr = 0;
   mus_tk_etaErr = 0;
   mus_tk_phiErr = 0;
   mus_stamu_chi2 = 0;
   mus_stamu_ndof = 0;
   mus_stamu_chg = 0;
   mus_stamu_pt = 0;
   mus_stamu_px = 0;
   mus_stamu_py = 0;
   mus_stamu_pz = 0;
   mus_stamu_eta = 0;
   mus_stamu_phi = 0;
   mus_stamu_theta = 0;
   mus_stamu_d0dum = 0;
   mus_stamu_dz = 0;
   mus_stamu_vx = 0;
   mus_stamu_vy = 0;
   mus_stamu_vz = 0;
   mus_stamu_numvalhits = 0;
   mus_stamu_numlosthits = 0;
   mus_stamu_d0dumErr = 0;
   mus_stamu_dzErr = 0;
   mus_stamu_ptErr = 0;
   mus_stamu_etaErr = 0;
   mus_stamu_phiErr = 0;
   mus_num_matches = 0;
   mus_id_All = 0;
   mus_id_AllGlobalMuons = 0;
   mus_id_AllStandAloneMuons = 0;
   mus_id_AllTrackerMuons = 0;
   mus_id_TrackerMuonArbitrated = 0;
   mus_id_AllArbitrated = 0;
   mus_id_GlobalMuonPromptTight = 0;
   mus_id_TMLastStationLoose = 0;
   mus_id_TMLastStationTight = 0;
   mus_id_TM2DCompatibilityLoose = 0;
   mus_id_TM2DCompatibilityTight = 0;
   mus_id_TMOneStationLoose = 0;
   mus_id_TMOneStationTight = 0;
   mus_id_TMLastStationOptimizedLowPtLoose = 0;
   mus_id_TMLastStationOptimizedLowPtTight = 0;
   photons_energy = 0;
   photons_et = 0;
   photons_eta = 0;
   photons_phi = 0;
   photons_pt = 0;
   photons_px = 0;
   photons_py = 0;
   photons_pz = 0;
   photons_status = 0;
   photons_theta = 0;
   photons_hadOverEM = 0;
   photons_scEnergy = 0;
   photons_scRawEnergy = 0;
   photons_scEta = 0;
   photons_scPhi = 0;
   photons_scEtaWidth = 0;
   photons_scPhiWidth = 0;
   photons_tIso = 0;
   photons_ecalIso = 0;
   photons_hcalIso = 0;
   photons_isoEcalRecHitDR04 = 0;
   photons_isoHcalRecHitDR04 = 0;
   photons_isoSolidTrkConeDR04 = 0;
   photons_isoHollowTrkConeDR04 = 0;
   photons_nTrkSolidConeDR04 = 0;
   photons_nTrkHollowConeDR04 = 0;
   photons_isoEcalRecHitDR03 = 0;
   photons_isoHcalRecHitDR03 = 0;
   photons_isoSolidTrkConeDR03 = 0;
   photons_isoHollowTrkConeDR03 = 0;
   photons_nTrkSolidConeDR03 = 0;
   photons_nTrkHollowConeDR03 = 0;
   photons_isAlsoElectron = 0;
   photons_hasPixelSeed = 0;
   photons_isConverted = 0;
   photons_isEBGap = 0;
   photons_isEEGap = 0;
   photons_isEBEEGap = 0;
   photons_isEBPho = 0;
   photons_isEEPho = 0;
   photons_isLoosePhoton = 0;
   photons_isTightPhoton = 0;
   photons_r9 = 0;
   photons_gen_et = 0;
   photons_gen_eta = 0;
   photons_gen_phi = 0;
   photons_gen_id = 0;
   pv_x = 0;
   pv_y = 0;
   pv_z = 0;
   pv_xErr = 0;
   pv_yErr = 0;
   pv_zErr = 0;
   taus_energy = 0;
   taus_et = 0;
   taus_eta = 0;
   taus_phi = 0;
   taus_pt = 0;
   taus_px = 0;
   taus_py = 0;
   taus_pz = 0;
   taus_status = 0;
   taus_theta = 0;
   taus_charge = 0;
   taus_emf = 0;
   taus_hcalTotOverPLead = 0;
   taus_hcalMaxOverPLead = 0;
   taus_hcal3x3OverPLead = 0;
   taus_ecalStripSumEOverPLead = 0;
   taus_elecPreIdOutput = 0;
   taus_elecPreIdDecision = 0;
   taus_leadPFChargedHadrCand_pt = 0;
   taus_leadPFChargedHadrCand_charge = 0;
   taus_leadPFChargedHadrCand_eta = 0;
   taus_leadPFChargedHadrCand_ECAL_eta = 0;
   taus_leadPFChargedHadrCand_phi = 0;
   taus_isoPFGammaCandsEtSum = 0;
   taus_isoPFChargedHadrCandsPtSum = 0;
   taus_leadingTrackFinding = 0;
   taus_leadingTrackPtCut = 0;
   taus_trackIsolation = 0;
   taus_ecalIsolation = 0;
   taus_byIsolation = 0;
   taus_againstElectron = 0;
   taus_againstMuon = 0;
   taus_taNC_quarter = 0;
   taus_taNC_one = 0;
   taus_taNC_half = 0;
   taus_taNC_tenth = 0;
   taus_muDecision = 0;
   taus_Nprongs = 0;
   tcmets_et = 0;
   tcmets_phi = 0;
   tcmets_ex = 0;
   tcmets_ey = 0;
   tcmets_sumEt = 0;
   tracks_chi2 = 0;
   tracks_ndof = 0;
   tracks_chg = 0;
   tracks_pt = 0;
   tracks_px = 0;
   tracks_py = 0;
   tracks_pz = 0;
   tracks_eta = 0;
   tracks_phi = 0;
   tracks_theta = 0;
   tracks_d0dum = 0;
   tracks_dz = 0;
   tracks_vx = 0;
   tracks_vy = 0;
   tracks_vz = 0;
   tracks_numvalhits = 0;
   tracks_numlosthits = 0;
   tracks_d0dumErr = 0;
   tracks_dzErr = 0;
   tracks_ptErr = 0;
   tracks_etaErr = 0;
   tracks_phiErr = 0;
   tracks_Nrechits = 0;
   tracks_innerHitX = 0;
   tracks_innerHitY = 0;
   tracks_innerHitZ = 0;
   tracks_outerHitX = 0;
   tracks_outerHitY = 0;
   tracks_outerHitZ = 0;
   tracks_highPurity = 0;

   fChain->SetBranchAddress("NbeamSpot", &NbeamSpot, &b_NbeamSpot);
   fChain->SetBranchAddress("beamSpot_x", &beamSpot_x, &b_beamSpot_x);
   fChain->SetBranchAddress("beamSpot_y", &beamSpot_y, &b_beamSpot_y);
   fChain->SetBranchAddress("beamSpot_z", &beamSpot_z, &b_beamSpot_z);
   fChain->SetBranchAddress("beamSpot_x0Error", &beamSpot_x0Error, &b_beamSpot_x0Error);
   fChain->SetBranchAddress("beamSpot_y0Error", &beamSpot_y0Error, &b_beamSpot_y0Error);
   fChain->SetBranchAddress("beamSpot_z0Error", &beamSpot_z0Error, &b_beamSpot_z0Error);
   fChain->SetBranchAddress("beamSpot_sigmaZ", &beamSpot_sigmaZ, &b_beamSpot_sigmaZ);
   fChain->SetBranchAddress("beamSpot_sigmaZ0Error", &beamSpot_sigmaZ0Error, &b_beamSpot_sigmaZ0Error);
   fChain->SetBranchAddress("beamSpot_dxdz", &beamSpot_dxdz, &b_beamSpot_dxdz);
   fChain->SetBranchAddress("beamSpot_dxdzError", &beamSpot_dxdzError, &b_beamSpot_dxdzError);
   fChain->SetBranchAddress("beamSpot_dydz", &beamSpot_dydz, &b_beamSpot_dydz);
   fChain->SetBranchAddress("beamSpot_dydzError", &beamSpot_dydzError, &b_beamSpot_dydzError);
   fChain->SetBranchAddress("beamSpot_beamWidthX", &beamSpot_beamWidthX, &b_beamSpot_beamWidthX);
   fChain->SetBranchAddress("beamSpot_beamWidthY", &beamSpot_beamWidthY, &b_beamSpot_beamWidthY);
   fChain->SetBranchAddress("beamSpot_beamWidthXError", &beamSpot_beamWidthXError, &b_beamSpot_beamWidthXError);
   fChain->SetBranchAddress("beamSpot_beamWidthYError", &beamSpot_beamWidthYError, &b_beamSpot_beamWidthYError);
   fChain->SetBranchAddress("Nels", &Nels, &b_Nels);
   fChain->SetBranchAddress("els_energy", &els_energy, &b_els_energy);
   fChain->SetBranchAddress("els_et", &els_et, &b_els_et);
   fChain->SetBranchAddress("els_eta", &els_eta, &b_els_eta);
   fChain->SetBranchAddress("els_phi", &els_phi, &b_els_phi);
   fChain->SetBranchAddress("els_pt", &els_pt, &b_els_pt);
   fChain->SetBranchAddress("els_px", &els_px, &b_els_px);
   fChain->SetBranchAddress("els_py", &els_py, &b_els_py);
   fChain->SetBranchAddress("els_pz", &els_pz, &b_els_pz);
   fChain->SetBranchAddress("els_status", &els_status, &b_els_status);
   fChain->SetBranchAddress("els_theta", &els_theta, &b_els_theta);
   fChain->SetBranchAddress("els_gen_id", &els_gen_id, &b_els_gen_id);
   fChain->SetBranchAddress("els_gen_phi", &els_gen_phi, &b_els_gen_phi);
   fChain->SetBranchAddress("els_gen_pt", &els_gen_pt, &b_els_gen_pt);
   fChain->SetBranchAddress("els_gen_pz", &els_gen_pz, &b_els_gen_pz);
   fChain->SetBranchAddress("els_gen_px", &els_gen_px, &b_els_gen_px);
   fChain->SetBranchAddress("els_gen_py", &els_gen_py, &b_els_gen_py);
   fChain->SetBranchAddress("els_gen_eta", &els_gen_eta, &b_els_gen_eta);
   fChain->SetBranchAddress("els_gen_theta", &els_gen_theta, &b_els_gen_theta);
   fChain->SetBranchAddress("els_gen_et", &els_gen_et, &b_els_gen_et);
   fChain->SetBranchAddress("els_gen_mother_id", &els_gen_mother_id, &b_els_gen_mother_id);
   fChain->SetBranchAddress("els_gen_mother_phi", &els_gen_mother_phi, &b_els_gen_mother_phi);
   fChain->SetBranchAddress("els_gen_mother_pt", &els_gen_mother_pt, &b_els_gen_mother_pt);
   fChain->SetBranchAddress("els_gen_mother_pz", &els_gen_mother_pz, &b_els_gen_mother_pz);
   fChain->SetBranchAddress("els_gen_mother_px", &els_gen_mother_px, &b_els_gen_mother_px);
   fChain->SetBranchAddress("els_gen_mother_py", &els_gen_mother_py, &b_els_gen_mother_py);
   fChain->SetBranchAddress("els_gen_mother_eta", &els_gen_mother_eta, &b_els_gen_mother_eta);
   fChain->SetBranchAddress("els_gen_mother_theta", &els_gen_mother_theta, &b_els_gen_mother_theta);
   fChain->SetBranchAddress("els_gen_mother_et", &els_gen_mother_et, &b_els_gen_mother_et);
   fChain->SetBranchAddress("els_tightId", &els_tightId, &b_els_tightId);
   fChain->SetBranchAddress("els_looseId", &els_looseId, &b_els_looseId);
   fChain->SetBranchAddress("els_robustTightId", &els_robustTightId, &b_els_robustTightId);
   fChain->SetBranchAddress("els_robustLooseId", &els_robustLooseId, &b_els_robustLooseId);
   fChain->SetBranchAddress("els_robustHighEnergyId", &els_robustHighEnergyId, &b_els_robustHighEnergyId);
   fChain->SetBranchAddress("els_cIso", &els_cIso, &b_els_cIso);
   fChain->SetBranchAddress("els_tIso", &els_tIso, &b_els_tIso);
   fChain->SetBranchAddress("els_ecalIso", &els_ecalIso, &b_els_ecalIso);
   fChain->SetBranchAddress("els_hcalIso", &els_hcalIso, &b_els_hcalIso);
   fChain->SetBranchAddress("els_chi2", &els_chi2, &b_els_chi2);
   fChain->SetBranchAddress("els_charge", &els_charge, &b_els_charge);
   fChain->SetBranchAddress("els_caloEnergy", &els_caloEnergy, &b_els_caloEnergy);
   fChain->SetBranchAddress("els_hadOverEm", &els_hadOverEm, &b_els_hadOverEm);
   fChain->SetBranchAddress("els_eOverPIn", &els_eOverPIn, &b_els_eOverPIn);
   fChain->SetBranchAddress("els_eSeedOverPOut", &els_eSeedOverPOut, &b_els_eSeedOverPOut);
   fChain->SetBranchAddress("els_eSCraw", &els_eSCraw, &b_els_eSCraw);
   fChain->SetBranchAddress("els_eSeed", &els_eSeed, &b_els_eSeed);
   fChain->SetBranchAddress("els_sigmaEtaEta", &els_sigmaEtaEta, &b_els_sigmaEtaEta);
   fChain->SetBranchAddress("els_sigmaIEtaIEta", &els_sigmaIEtaIEta, &b_els_sigmaIEtaIEta);
   fChain->SetBranchAddress("els_scE1x5", &els_scE1x5, &b_els_scE1x5);
   fChain->SetBranchAddress("els_scE2x5Max", &els_scE2x5Max, &b_els_scE2x5Max);
   fChain->SetBranchAddress("els_scE5x5", &els_scE5x5, &b_els_scE5x5);
   fChain->SetBranchAddress("els_dEtaIn", &els_dEtaIn, &b_els_dEtaIn);
   fChain->SetBranchAddress("els_dPhiIn", &els_dPhiIn, &b_els_dPhiIn);
   fChain->SetBranchAddress("els_dEtaOut", &els_dEtaOut, &b_els_dEtaOut);
   fChain->SetBranchAddress("els_dPhiOut", &els_dPhiOut, &b_els_dPhiOut);
   fChain->SetBranchAddress("els_numvalhits", &els_numvalhits, &b_els_numvalhits);
   fChain->SetBranchAddress("els_numlosthits", &els_numlosthits, &b_els_numlosthits);
   fChain->SetBranchAddress("els_basicClustersSize", &els_basicClustersSize, &b_els_basicClustersSize);
   fChain->SetBranchAddress("els_tk_pt", &els_tk_pt, &b_els_tk_pt);
   fChain->SetBranchAddress("els_tk_phi", &els_tk_phi, &b_els_tk_phi);
   fChain->SetBranchAddress("els_tk_eta", &els_tk_eta, &b_els_tk_eta);
   fChain->SetBranchAddress("els_d0dum", &els_d0dum, &b_els_d0dum);
   fChain->SetBranchAddress("els_dz", &els_dz, &b_els_dz);
   fChain->SetBranchAddress("els_vx", &els_vx, &b_els_vx);
   fChain->SetBranchAddress("els_vy", &els_vy, &b_els_vy);
   fChain->SetBranchAddress("els_vz", &els_vz, &b_els_vz);
   fChain->SetBranchAddress("els_ndof", &els_ndof, &b_els_ndof);
   fChain->SetBranchAddress("els_ptError", &els_ptError, &b_els_ptError);
   fChain->SetBranchAddress("els_d0dumError", &els_d0dumError, &b_els_d0dumError);
   fChain->SetBranchAddress("els_dzError", &els_dzError, &b_els_dzError);
   fChain->SetBranchAddress("els_etaError", &els_etaError, &b_els_etaError);
   fChain->SetBranchAddress("els_phiError", &els_phiError, &b_els_phiError);
   fChain->SetBranchAddress("els_tk_charge", &els_tk_charge, &b_els_tk_charge);
   fChain->SetBranchAddress("els_ctf_tk_id", &els_ctf_tk_id, &b_els_ctf_tk_id);
   fChain->SetBranchAddress("els_ctf_tk_charge", &els_ctf_tk_charge, &b_els_ctf_tk_charge);
   fChain->SetBranchAddress("els_ctf_tk_eta", &els_ctf_tk_eta, &b_els_ctf_tk_eta);
   fChain->SetBranchAddress("els_ctf_tk_phi", &els_ctf_tk_phi, &b_els_ctf_tk_phi);
   fChain->SetBranchAddress("els_fbrem", &els_fbrem, &b_els_fbrem);
   fChain->SetBranchAddress("els_shFracInnerHits", &els_shFracInnerHits, &b_els_shFracInnerHits);
   fChain->SetBranchAddress("els_dr03EcalRecHitSumEt", &els_dr03EcalRecHitSumEt, &b_els_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("els_dr03HcalTowerSumEt", &els_dr03HcalTowerSumEt, &b_els_dr03HcalTowerSumEt);
   fChain->SetBranchAddress("els_dr03HcalDepth1TowerSumEt", &els_dr03HcalDepth1TowerSumEt, &b_els_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("els_dr03HcalDepth2TowerSumEt", &els_dr03HcalDepth2TowerSumEt, &b_els_dr03HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("els_dr03TkSumPt", &els_dr03TkSumPt, &b_els_dr03TkSumPt);
   fChain->SetBranchAddress("els_dr04EcalRecHitSumEt", &els_dr04EcalRecHitSumEt, &b_els_dr04EcalRecHitSumEt);
   fChain->SetBranchAddress("els_dr04HcalTowerSumEt", &els_dr04HcalTowerSumEt, &b_els_dr04HcalTowerSumEt);
   fChain->SetBranchAddress("els_dr04HcalDepth1TowerSumEt", &els_dr04HcalDepth1TowerSumEt, &b_els_dr04HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("els_dr04HcalDepth2TowerSumEt", &els_dr04HcalDepth2TowerSumEt, &b_els_dr04HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("els_dr04TkSumPt", &els_dr04TkSumPt, &b_els_dr04TkSumPt);
   fChain->SetBranchAddress("els_cpx", &els_cpx, &b_els_cpx);
   fChain->SetBranchAddress("els_cpy", &els_cpy, &b_els_cpy);
   fChain->SetBranchAddress("els_cpz", &els_cpz, &b_els_cpz);
   fChain->SetBranchAddress("els_vpx", &els_vpx, &b_els_vpx);
   fChain->SetBranchAddress("els_vpy", &els_vpy, &b_els_vpy);
   fChain->SetBranchAddress("els_vpz", &els_vpz, &b_els_vpz);
   fChain->SetBranchAddress("els_cx", &els_cx, &b_els_cx);
   fChain->SetBranchAddress("els_cy", &els_cy, &b_els_cy);
   fChain->SetBranchAddress("els_cz", &els_cz, &b_els_cz);
   fChain->SetBranchAddress("Njets_AK5", &Njets_AK5, &b_Njets_AK5);
   fChain->SetBranchAddress("jets_AK5_energy", &jets_AK5_energy, &b_jets_AK5_energy);
   fChain->SetBranchAddress("jets_AK5_et", &jets_AK5_et, &b_jets_AK5_et);
   fChain->SetBranchAddress("jets_AK5_eta", &jets_AK5_eta, &b_jets_AK5_eta);
   fChain->SetBranchAddress("jets_AK5_phi", &jets_AK5_phi, &b_jets_AK5_phi);
   fChain->SetBranchAddress("jets_AK5_pt", &jets_AK5_pt, &b_jets_AK5_pt);
   fChain->SetBranchAddress("jets_AK5_px", &jets_AK5_px, &b_jets_AK5_px);
   fChain->SetBranchAddress("jets_AK5_py", &jets_AK5_py, &b_jets_AK5_py);
   fChain->SetBranchAddress("jets_AK5_pz", &jets_AK5_pz, &b_jets_AK5_pz);
   fChain->SetBranchAddress("jets_AK5_status", &jets_AK5_status, &b_jets_AK5_status);
   fChain->SetBranchAddress("jets_AK5_theta", &jets_AK5_theta, &b_jets_AK5_theta);
   fChain->SetBranchAddress("jets_AK5_parton_Id", &jets_AK5_parton_Id, &b_jets_AK5_parton_Id);
   fChain->SetBranchAddress("jets_AK5_parton_motherId", &jets_AK5_parton_motherId, &b_jets_AK5_parton_motherId);
   fChain->SetBranchAddress("jets_AK5_parton_pt", &jets_AK5_parton_pt, &b_jets_AK5_parton_pt);
   fChain->SetBranchAddress("jets_AK5_parton_phi", &jets_AK5_parton_phi, &b_jets_AK5_parton_phi);
   fChain->SetBranchAddress("jets_AK5_parton_eta", &jets_AK5_parton_eta, &b_jets_AK5_parton_eta);
   fChain->SetBranchAddress("jets_AK5_parton_Energy", &jets_AK5_parton_Energy, &b_jets_AK5_parton_Energy);
   fChain->SetBranchAddress("jets_AK5_parton_mass", &jets_AK5_parton_mass, &b_jets_AK5_parton_mass);
   fChain->SetBranchAddress("jets_AK5_parton_motherID", &jets_AK5_parton_motherID, &b_jets_AK5_parton_motherID);
   fChain->SetBranchAddress("jets_AK5_gen_et", &jets_AK5_gen_et, &b_jets_AK5_gen_et);
   fChain->SetBranchAddress("jets_AK5_gen_pt", &jets_AK5_gen_pt, &b_jets_AK5_gen_pt);
   fChain->SetBranchAddress("jets_AK5_gen_eta", &jets_AK5_gen_eta, &b_jets_AK5_gen_eta);
   fChain->SetBranchAddress("jets_AK5_gen_phi", &jets_AK5_gen_phi, &b_jets_AK5_gen_phi);
   fChain->SetBranchAddress("jets_AK5_gen_mass", &jets_AK5_gen_mass, &b_jets_AK5_gen_mass);
   fChain->SetBranchAddress("jets_AK5_gen_Energy", &jets_AK5_gen_Energy, &b_jets_AK5_gen_Energy);
   fChain->SetBranchAddress("jets_AK5_gen_Id", &jets_AK5_gen_Id, &b_jets_AK5_gen_Id);
   fChain->SetBranchAddress("jets_AK5_gen_motherID", &jets_AK5_gen_motherID, &b_jets_AK5_gen_motherID);
   fChain->SetBranchAddress("jets_AK5_gen_threeCharge", &jets_AK5_gen_threeCharge, &b_jets_AK5_gen_threeCharge);
   fChain->SetBranchAddress("jets_AK5_partonFlavour", &jets_AK5_partonFlavour, &b_jets_AK5_partonFlavour);
   fChain->SetBranchAddress("jets_AK5_btag_TC_highPur", &jets_AK5_btag_TC_highPur, &b_jets_AK5_btag_TC_highPur);
   fChain->SetBranchAddress("jets_AK5_btag_TC_highEff", &jets_AK5_btag_TC_highEff, &b_jets_AK5_btag_TC_highEff);
   fChain->SetBranchAddress("jets_AK5_btag_jetProb", &jets_AK5_btag_jetProb, &b_jets_AK5_btag_jetProb);
   fChain->SetBranchAddress("jets_AK5_btag_jetBProb", &jets_AK5_btag_jetBProb, &b_jets_AK5_btag_jetBProb);
   fChain->SetBranchAddress("jets_AK5_btag_softEle", &jets_AK5_btag_softEle, &b_jets_AK5_btag_softEle);
   fChain->SetBranchAddress("jets_AK5_btag_softMuon", &jets_AK5_btag_softMuon, &b_jets_AK5_btag_softMuon);
   fChain->SetBranchAddress("jets_AK5_btag_softMuonNoIP", &jets_AK5_btag_softMuonNoIP, &b_jets_AK5_btag_softMuonNoIP);
   fChain->SetBranchAddress("jets_AK5_btag_secVertex", &jets_AK5_btag_secVertex, &b_jets_AK5_btag_secVertex);
   fChain->SetBranchAddress("jets_AK5_chgEmE", &jets_AK5_chgEmE, &b_jets_AK5_chgEmE);
   fChain->SetBranchAddress("jets_AK5_chgHadE", &jets_AK5_chgHadE, &b_jets_AK5_chgHadE);
   fChain->SetBranchAddress("jets_AK5_chgMuE", &jets_AK5_chgMuE, &b_jets_AK5_chgMuE);
   fChain->SetBranchAddress("jets_AK5_chg_Mult", &jets_AK5_chg_Mult, &b_jets_AK5_chg_Mult);
   fChain->SetBranchAddress("jets_AK5_neutralEmE", &jets_AK5_neutralEmE, &b_jets_AK5_neutralEmE);
   fChain->SetBranchAddress("jets_AK5_neutralHadE", &jets_AK5_neutralHadE, &b_jets_AK5_neutralHadE);
   fChain->SetBranchAddress("jets_AK5_neutral_Mult", &jets_AK5_neutral_Mult, &b_jets_AK5_neutral_Mult);
   fChain->SetBranchAddress("jets_AK5_mu_Mult", &jets_AK5_mu_Mult, &b_jets_AK5_mu_Mult);
   fChain->SetBranchAddress("jets_AK5_emf", &jets_AK5_emf, &b_jets_AK5_emf);
   fChain->SetBranchAddress("jets_AK5_ehf", &jets_AK5_ehf, &b_jets_AK5_ehf);
   fChain->SetBranchAddress("jets_AK5_n60", &jets_AK5_n60, &b_jets_AK5_n60);
   fChain->SetBranchAddress("jets_AK5_n90", &jets_AK5_n90, &b_jets_AK5_n90);
   fChain->SetBranchAddress("jets_AK5_area", &jets_AK5_area, &b_jets_AK5_area);
   fChain->SetBranchAddress("jets_AK5_mass", &jets_AK5_mass, &b_jets_AK5_mass);
   fChain->SetBranchAddress("Njets_SC5", &Njets_SC5, &b_Njets_SC5);
   fChain->SetBranchAddress("jets_SC5_energy", &jets_SC5_energy, &b_jets_SC5_energy);
   fChain->SetBranchAddress("jets_SC5_et", &jets_SC5_et, &b_jets_SC5_et);
   fChain->SetBranchAddress("jets_SC5_eta", &jets_SC5_eta, &b_jets_SC5_eta);
   fChain->SetBranchAddress("jets_SC5_phi", &jets_SC5_phi, &b_jets_SC5_phi);
   fChain->SetBranchAddress("jets_SC5_pt", &jets_SC5_pt, &b_jets_SC5_pt);
   fChain->SetBranchAddress("jets_SC5_px", &jets_SC5_px, &b_jets_SC5_px);
   fChain->SetBranchAddress("jets_SC5_py", &jets_SC5_py, &b_jets_SC5_py);
   fChain->SetBranchAddress("jets_SC5_pz", &jets_SC5_pz, &b_jets_SC5_pz);
   fChain->SetBranchAddress("jets_SC5_status", &jets_SC5_status, &b_jets_SC5_status);
   fChain->SetBranchAddress("jets_SC5_theta", &jets_SC5_theta, &b_jets_SC5_theta);
   fChain->SetBranchAddress("jets_SC5_parton_Id", &jets_SC5_parton_Id, &b_jets_SC5_parton_Id);
   fChain->SetBranchAddress("jets_SC5_parton_motherId", &jets_SC5_parton_motherId, &b_jets_SC5_parton_motherId);
   fChain->SetBranchAddress("jets_SC5_parton_pt", &jets_SC5_parton_pt, &b_jets_SC5_parton_pt);
   fChain->SetBranchAddress("jets_SC5_parton_phi", &jets_SC5_parton_phi, &b_jets_SC5_parton_phi);
   fChain->SetBranchAddress("jets_SC5_parton_eta", &jets_SC5_parton_eta, &b_jets_SC5_parton_eta);
   fChain->SetBranchAddress("jets_SC5_parton_Energy", &jets_SC5_parton_Energy, &b_jets_SC5_parton_Energy);
   fChain->SetBranchAddress("jets_SC5_parton_mass", &jets_SC5_parton_mass, &b_jets_SC5_parton_mass);
   fChain->SetBranchAddress("jets_SC5_parton_motherID", &jets_SC5_parton_motherID, &b_jets_SC5_parton_motherID);
   fChain->SetBranchAddress("jets_SC5_gen_et", &jets_SC5_gen_et, &b_jets_SC5_gen_et);
   fChain->SetBranchAddress("jets_SC5_gen_pt", &jets_SC5_gen_pt, &b_jets_SC5_gen_pt);
   fChain->SetBranchAddress("jets_SC5_gen_eta", &jets_SC5_gen_eta, &b_jets_SC5_gen_eta);
   fChain->SetBranchAddress("jets_SC5_gen_phi", &jets_SC5_gen_phi, &b_jets_SC5_gen_phi);
   fChain->SetBranchAddress("jets_SC5_gen_mass", &jets_SC5_gen_mass, &b_jets_SC5_gen_mass);
   fChain->SetBranchAddress("jets_SC5_gen_Energy", &jets_SC5_gen_Energy, &b_jets_SC5_gen_Energy);
   fChain->SetBranchAddress("jets_SC5_gen_Id", &jets_SC5_gen_Id, &b_jets_SC5_gen_Id);
   fChain->SetBranchAddress("jets_SC5_gen_motherID", &jets_SC5_gen_motherID, &b_jets_SC5_gen_motherID);
   fChain->SetBranchAddress("jets_SC5_gen_threeCharge", &jets_SC5_gen_threeCharge, &b_jets_SC5_gen_threeCharge);
   fChain->SetBranchAddress("jets_SC5_partonFlavour", &jets_SC5_partonFlavour, &b_jets_SC5_partonFlavour);
   fChain->SetBranchAddress("jets_SC5_btag_TC_highPur", &jets_SC5_btag_TC_highPur, &b_jets_SC5_btag_TC_highPur);
   fChain->SetBranchAddress("jets_SC5_btag_TC_highEff", &jets_SC5_btag_TC_highEff, &b_jets_SC5_btag_TC_highEff);
   fChain->SetBranchAddress("jets_SC5_btag_jetProb", &jets_SC5_btag_jetProb, &b_jets_SC5_btag_jetProb);
   fChain->SetBranchAddress("jets_SC5_btag_jetBProb", &jets_SC5_btag_jetBProb, &b_jets_SC5_btag_jetBProb);
   fChain->SetBranchAddress("jets_SC5_btag_softEle", &jets_SC5_btag_softEle, &b_jets_SC5_btag_softEle);
   fChain->SetBranchAddress("jets_SC5_btag_softMuon", &jets_SC5_btag_softMuon, &b_jets_SC5_btag_softMuon);
   fChain->SetBranchAddress("jets_SC5_btag_softMuonNoIP", &jets_SC5_btag_softMuonNoIP, &b_jets_SC5_btag_softMuonNoIP);
   fChain->SetBranchAddress("jets_SC5_btag_secVertex", &jets_SC5_btag_secVertex, &b_jets_SC5_btag_secVertex);
   fChain->SetBranchAddress("jets_SC5_chgEmE", &jets_SC5_chgEmE, &b_jets_SC5_chgEmE);
   fChain->SetBranchAddress("jets_SC5_chgHadE", &jets_SC5_chgHadE, &b_jets_SC5_chgHadE);
   fChain->SetBranchAddress("jets_SC5_chgMuE", &jets_SC5_chgMuE, &b_jets_SC5_chgMuE);
   fChain->SetBranchAddress("jets_SC5_chg_Mult", &jets_SC5_chg_Mult, &b_jets_SC5_chg_Mult);
   fChain->SetBranchAddress("jets_SC5_neutralEmE", &jets_SC5_neutralEmE, &b_jets_SC5_neutralEmE);
   fChain->SetBranchAddress("jets_SC5_neutralHadE", &jets_SC5_neutralHadE, &b_jets_SC5_neutralHadE);
   fChain->SetBranchAddress("jets_SC5_neutral_Mult", &jets_SC5_neutral_Mult, &b_jets_SC5_neutral_Mult);
   fChain->SetBranchAddress("jets_SC5_mu_Mult", &jets_SC5_mu_Mult, &b_jets_SC5_mu_Mult);
   fChain->SetBranchAddress("jets_SC5_emf", &jets_SC5_emf, &b_jets_SC5_emf);
   fChain->SetBranchAddress("jets_SC5_ehf", &jets_SC5_ehf, &b_jets_SC5_ehf);
   fChain->SetBranchAddress("jets_SC5_n60", &jets_SC5_n60, &b_jets_SC5_n60);
   fChain->SetBranchAddress("jets_SC5_n90", &jets_SC5_n90, &b_jets_SC5_n90);
   fChain->SetBranchAddress("jets_SC5_area", &jets_SC5_area, &b_jets_SC5_area);
   fChain->SetBranchAddress("jets_SC5_mass", &jets_SC5_mass, &b_jets_SC5_mass);
   fChain->SetBranchAddress("Nmc_doc", &Nmc_doc, &b_Nmc_doc);
   fChain->SetBranchAddress("mc_doc_id", &mc_doc_id, &b_mc_doc_id);
   fChain->SetBranchAddress("mc_doc_pt", &mc_doc_pt, &b_mc_doc_pt);
   fChain->SetBranchAddress("mc_doc_px", &mc_doc_px, &b_mc_doc_px);
   fChain->SetBranchAddress("mc_doc_py", &mc_doc_py, &b_mc_doc_py);
   fChain->SetBranchAddress("mc_doc_pz", &mc_doc_pz, &b_mc_doc_pz);
   fChain->SetBranchAddress("mc_doc_eta", &mc_doc_eta, &b_mc_doc_eta);
   fChain->SetBranchAddress("mc_doc_phi", &mc_doc_phi, &b_mc_doc_phi);
   fChain->SetBranchAddress("mc_doc_theta", &mc_doc_theta, &b_mc_doc_theta);
   fChain->SetBranchAddress("mc_doc_energy", &mc_doc_energy, &b_mc_doc_energy);
   fChain->SetBranchAddress("mc_doc_status", &mc_doc_status, &b_mc_doc_status);
   fChain->SetBranchAddress("mc_doc_charge", &mc_doc_charge, &b_mc_doc_charge);
   fChain->SetBranchAddress("mc_doc_mother_id", &mc_doc_mother_id, &b_mc_doc_mother_id);
   fChain->SetBranchAddress("mc_doc_grandmother_id", &mc_doc_grandmother_id, &b_mc_doc_grandmother_id);
   fChain->SetBranchAddress("mc_doc_ggrandmother_id", &mc_doc_ggrandmother_id, &b_mc_doc_ggrandmother_id);
   fChain->SetBranchAddress("mc_doc_mother_pt", &mc_doc_mother_pt, &b_mc_doc_mother_pt);
   fChain->SetBranchAddress("mc_doc_vertex_x", &mc_doc_vertex_x, &b_mc_doc_vertex_x);
   fChain->SetBranchAddress("mc_doc_vertex_y", &mc_doc_vertex_y, &b_mc_doc_vertex_y);
   fChain->SetBranchAddress("mc_doc_vertex_z", &mc_doc_vertex_z, &b_mc_doc_vertex_z);
   fChain->SetBranchAddress("mc_doc_mass", &mc_doc_mass, &b_mc_doc_mass);
   fChain->SetBranchAddress("mc_doc_numOfDaughters", &mc_doc_numOfDaughters, &b_mc_doc_numOfDaughters);
   fChain->SetBranchAddress("mc_doc_numOfMothers", &mc_doc_numOfMothers, &b_mc_doc_numOfMothers);
   fChain->SetBranchAddress("Nmc_electrons", &Nmc_electrons, &b_Nmc_electrons);
   fChain->SetBranchAddress("mc_electrons_id", &mc_electrons_id, &b_mc_electrons_id);
   fChain->SetBranchAddress("mc_electrons_pt", &mc_electrons_pt, &b_mc_electrons_pt);
   fChain->SetBranchAddress("mc_electrons_px", &mc_electrons_px, &b_mc_electrons_px);
   fChain->SetBranchAddress("mc_electrons_py", &mc_electrons_py, &b_mc_electrons_py);
   fChain->SetBranchAddress("mc_electrons_pz", &mc_electrons_pz, &b_mc_electrons_pz);
   fChain->SetBranchAddress("mc_electrons_eta", &mc_electrons_eta, &b_mc_electrons_eta);
   fChain->SetBranchAddress("mc_electrons_phi", &mc_electrons_phi, &b_mc_electrons_phi);
   fChain->SetBranchAddress("mc_electrons_theta", &mc_electrons_theta, &b_mc_electrons_theta);
   fChain->SetBranchAddress("mc_electrons_status", &mc_electrons_status, &b_mc_electrons_status);
   fChain->SetBranchAddress("mc_electrons_energy", &mc_electrons_energy, &b_mc_electrons_energy);
   fChain->SetBranchAddress("mc_electrons_charge", &mc_electrons_charge, &b_mc_electrons_charge);
   fChain->SetBranchAddress("mc_electrons_mother_id", &mc_electrons_mother_id, &b_mc_electrons_mother_id);
   fChain->SetBranchAddress("mc_electrons_mother_pt", &mc_electrons_mother_pt, &b_mc_electrons_mother_pt);
   fChain->SetBranchAddress("mc_electrons_grandmother_id", &mc_electrons_grandmother_id, &b_mc_electrons_grandmother_id);
   fChain->SetBranchAddress("mc_electrons_ggrandmother_id", &mc_electrons_ggrandmother_id, &b_mc_electrons_ggrandmother_id);
   fChain->SetBranchAddress("mc_electrons_vertex_x", &mc_electrons_vertex_x, &b_mc_electrons_vertex_x);
   fChain->SetBranchAddress("mc_electrons_vertex_y", &mc_electrons_vertex_y, &b_mc_electrons_vertex_y);
   fChain->SetBranchAddress("mc_electrons_vertex_z", &mc_electrons_vertex_z, &b_mc_electrons_vertex_z);
   fChain->SetBranchAddress("mc_electrons_mass", &mc_electrons_mass, &b_mc_electrons_mass);
   fChain->SetBranchAddress("mc_electrons_numOfDaughters", &mc_electrons_numOfDaughters, &b_mc_electrons_numOfDaughters);
   fChain->SetBranchAddress("Nmc_mus", &Nmc_mus, &b_Nmc_mus);
   fChain->SetBranchAddress("mc_mus_id", &mc_mus_id, &b_mc_mus_id);
   fChain->SetBranchAddress("mc_mus_pt", &mc_mus_pt, &b_mc_mus_pt);
   fChain->SetBranchAddress("mc_mus_px", &mc_mus_px, &b_mc_mus_px);
   fChain->SetBranchAddress("mc_mus_py", &mc_mus_py, &b_mc_mus_py);
   fChain->SetBranchAddress("mc_mus_pz", &mc_mus_pz, &b_mc_mus_pz);
   fChain->SetBranchAddress("mc_mus_eta", &mc_mus_eta, &b_mc_mus_eta);
   fChain->SetBranchAddress("mc_mus_phi", &mc_mus_phi, &b_mc_mus_phi);
   fChain->SetBranchAddress("mc_mus_theta", &mc_mus_theta, &b_mc_mus_theta);
   fChain->SetBranchAddress("mc_mus_status", &mc_mus_status, &b_mc_mus_status);
   fChain->SetBranchAddress("mc_mus_energy", &mc_mus_energy, &b_mc_mus_energy);
   fChain->SetBranchAddress("mc_mus_charge", &mc_mus_charge, &b_mc_mus_charge);
   fChain->SetBranchAddress("mc_mus_mother_id", &mc_mus_mother_id, &b_mc_mus_mother_id);
   fChain->SetBranchAddress("mc_mus_mother_pt", &mc_mus_mother_pt, &b_mc_mus_mother_pt);
   fChain->SetBranchAddress("mc_mus_grandmother_id", &mc_mus_grandmother_id, &b_mc_mus_grandmother_id);
   fChain->SetBranchAddress("mc_mus_ggrandmother_id", &mc_mus_ggrandmother_id, &b_mc_mus_ggrandmother_id);
   fChain->SetBranchAddress("mc_mus_vertex_x", &mc_mus_vertex_x, &b_mc_mus_vertex_x);
   fChain->SetBranchAddress("mc_mus_vertex_y", &mc_mus_vertex_y, &b_mc_mus_vertex_y);
   fChain->SetBranchAddress("mc_mus_vertex_z", &mc_mus_vertex_z, &b_mc_mus_vertex_z);
   fChain->SetBranchAddress("mc_mus_mass", &mc_mus_mass, &b_mc_mus_mass);
   fChain->SetBranchAddress("mc_mus_numOfDaughters", &mc_mus_numOfDaughters, &b_mc_mus_numOfDaughters);
   fChain->SetBranchAddress("Nmets_AK5", &Nmets_AK5, &b_Nmets_AK5);
   fChain->SetBranchAddress("mets_AK5_et", &mets_AK5_et, &b_mets_AK5_et);
   fChain->SetBranchAddress("mets_AK5_phi", &mets_AK5_phi, &b_mets_AK5_phi);
   fChain->SetBranchAddress("mets_AK5_ex", &mets_AK5_ex, &b_mets_AK5_ex);
   fChain->SetBranchAddress("mets_AK5_ey", &mets_AK5_ey, &b_mets_AK5_ey);
   fChain->SetBranchAddress("mets_AK5_gen_et", &mets_AK5_gen_et, &b_mets_AK5_gen_et);
   fChain->SetBranchAddress("mets_AK5_gen_phi", &mets_AK5_gen_phi, &b_mets_AK5_gen_phi);
   fChain->SetBranchAddress("mets_AK5_sign", &mets_AK5_sign, &b_mets_AK5_sign);
   fChain->SetBranchAddress("mets_AK5_sumEt", &mets_AK5_sumEt, &b_mets_AK5_sumEt);
   fChain->SetBranchAddress("mets_AK5_unCPhi", &mets_AK5_unCPhi, &b_mets_AK5_unCPhi);
   fChain->SetBranchAddress("mets_AK5_unCPt", &mets_AK5_unCPt, &b_mets_AK5_unCPt);
   fChain->SetBranchAddress("Nmets_SC5", &Nmets_SC5, &b_Nmets_SC5);
   fChain->SetBranchAddress("mets_SC5_et", &mets_SC5_et, &b_mets_SC5_et);
   fChain->SetBranchAddress("mets_SC5_phi", &mets_SC5_phi, &b_mets_SC5_phi);
   fChain->SetBranchAddress("mets_SC5_ex", &mets_SC5_ex, &b_mets_SC5_ex);
   fChain->SetBranchAddress("mets_SC5_ey", &mets_SC5_ey, &b_mets_SC5_ey);
   fChain->SetBranchAddress("mets_SC5_gen_et", &mets_SC5_gen_et, &b_mets_SC5_gen_et);
   fChain->SetBranchAddress("mets_SC5_gen_phi", &mets_SC5_gen_phi, &b_mets_SC5_gen_phi);
   fChain->SetBranchAddress("mets_SC5_sign", &mets_SC5_sign, &b_mets_SC5_sign);
   fChain->SetBranchAddress("mets_SC5_sumEt", &mets_SC5_sumEt, &b_mets_SC5_sumEt);
   fChain->SetBranchAddress("mets_SC5_unCPhi", &mets_SC5_unCPhi, &b_mets_SC5_unCPhi);
   fChain->SetBranchAddress("mets_SC5_unCPt", &mets_SC5_unCPt, &b_mets_SC5_unCPt);
   fChain->SetBranchAddress("Nmus", &Nmus, &b_Nmus);
   fChain->SetBranchAddress("mus_energy", &mus_energy, &b_mus_energy);
   fChain->SetBranchAddress("mus_et", &mus_et, &b_mus_et);
   fChain->SetBranchAddress("mus_eta", &mus_eta, &b_mus_eta);
   fChain->SetBranchAddress("mus_phi", &mus_phi, &b_mus_phi);
   fChain->SetBranchAddress("mus_pt", &mus_pt, &b_mus_pt);
   fChain->SetBranchAddress("mus_px", &mus_px, &b_mus_px);
   fChain->SetBranchAddress("mus_py", &mus_py, &b_mus_py);
   fChain->SetBranchAddress("mus_pz", &mus_pz, &b_mus_pz);
   fChain->SetBranchAddress("mus_status", &mus_status, &b_mus_status);
   fChain->SetBranchAddress("mus_theta", &mus_theta, &b_mus_theta);
   fChain->SetBranchAddress("mus_gen_id", &mus_gen_id, &b_mus_gen_id);
   fChain->SetBranchAddress("mus_gen_phi", &mus_gen_phi, &b_mus_gen_phi);
   fChain->SetBranchAddress("mus_gen_pt", &mus_gen_pt, &b_mus_gen_pt);
   fChain->SetBranchAddress("mus_gen_pz", &mus_gen_pz, &b_mus_gen_pz);
   fChain->SetBranchAddress("mus_gen_px", &mus_gen_px, &b_mus_gen_px);
   fChain->SetBranchAddress("mus_gen_py", &mus_gen_py, &b_mus_gen_py);
   fChain->SetBranchAddress("mus_gen_eta", &mus_gen_eta, &b_mus_gen_eta);
   fChain->SetBranchAddress("mus_gen_theta", &mus_gen_theta, &b_mus_gen_theta);
   fChain->SetBranchAddress("mus_gen_et", &mus_gen_et, &b_mus_gen_et);
   fChain->SetBranchAddress("mus_gen_mother_id", &mus_gen_mother_id, &b_mus_gen_mother_id);
   fChain->SetBranchAddress("mus_gen_mother_phi", &mus_gen_mother_phi, &b_mus_gen_mother_phi);
   fChain->SetBranchAddress("mus_gen_mother_pt", &mus_gen_mother_pt, &b_mus_gen_mother_pt);
   fChain->SetBranchAddress("mus_gen_mother_pz", &mus_gen_mother_pz, &b_mus_gen_mother_pz);
   fChain->SetBranchAddress("mus_gen_mother_px", &mus_gen_mother_px, &b_mus_gen_mother_px);
   fChain->SetBranchAddress("mus_gen_mother_py", &mus_gen_mother_py, &b_mus_gen_mother_py);
   fChain->SetBranchAddress("mus_gen_mother_eta", &mus_gen_mother_eta, &b_mus_gen_mother_eta);
   fChain->SetBranchAddress("mus_gen_mother_theta", &mus_gen_mother_theta, &b_mus_gen_mother_theta);
   fChain->SetBranchAddress("mus_gen_mother_et", &mus_gen_mother_et, &b_mus_gen_mother_et);
   fChain->SetBranchAddress("mus_tkHits", &mus_tkHits, &b_mus_tkHits);
   fChain->SetBranchAddress("mus_cIso", &mus_cIso, &b_mus_cIso);
   fChain->SetBranchAddress("mus_tIso", &mus_tIso, &b_mus_tIso);
   fChain->SetBranchAddress("mus_ecalIso", &mus_ecalIso, &b_mus_ecalIso);
   fChain->SetBranchAddress("mus_hcalIso", &mus_hcalIso, &b_mus_hcalIso);
   fChain->SetBranchAddress("mus_ecalvetoDep", &mus_ecalvetoDep, &b_mus_ecalvetoDep);
   fChain->SetBranchAddress("mus_hcalvetoDep", &mus_hcalvetoDep, &b_mus_hcalvetoDep);
   fChain->SetBranchAddress("mus_calEnergyEm", &mus_calEnergyEm, &b_mus_calEnergyEm);
   fChain->SetBranchAddress("mus_calEnergyHad", &mus_calEnergyHad, &b_mus_calEnergyHad);
   fChain->SetBranchAddress("mus_calEnergyHo", &mus_calEnergyHo, &b_mus_calEnergyHo);
   fChain->SetBranchAddress("mus_calEnergyEmS9", &mus_calEnergyEmS9, &b_mus_calEnergyEmS9);
   fChain->SetBranchAddress("mus_calEnergyHadS9", &mus_calEnergyHadS9, &b_mus_calEnergyHadS9);
   fChain->SetBranchAddress("mus_calEnergyHoS9", &mus_calEnergyHoS9, &b_mus_calEnergyHoS9);
   fChain->SetBranchAddress("mus_iso03_sumPt", &mus_iso03_sumPt, &b_mus_iso03_sumPt);
   fChain->SetBranchAddress("mus_iso03_emEt", &mus_iso03_emEt, &b_mus_iso03_emEt);
   fChain->SetBranchAddress("mus_iso03_hadEt", &mus_iso03_hadEt, &b_mus_iso03_hadEt);
   fChain->SetBranchAddress("mus_iso03_hoEt", &mus_iso03_hoEt, &b_mus_iso03_hoEt);
   fChain->SetBranchAddress("mus_iso03_nTracks", &mus_iso03_nTracks, &b_mus_iso03_nTracks);
   fChain->SetBranchAddress("mus_iso05_sumPt", &mus_iso05_sumPt, &b_mus_iso05_sumPt);
   fChain->SetBranchAddress("mus_iso05_emEt", &mus_iso05_emEt, &b_mus_iso05_emEt);
   fChain->SetBranchAddress("mus_iso05_hadEt", &mus_iso05_hadEt, &b_mus_iso05_hadEt);
   fChain->SetBranchAddress("mus_iso05_hoEt", &mus_iso05_hoEt, &b_mus_iso05_hoEt);
   fChain->SetBranchAddress("mus_iso05_nTracks", &mus_iso05_nTracks, &b_mus_iso05_nTracks);
   fChain->SetBranchAddress("mus_charge", &mus_charge, &b_mus_charge);
   fChain->SetBranchAddress("mus_cm_chi2", &mus_cm_chi2, &b_mus_cm_chi2);
   fChain->SetBranchAddress("mus_cm_ndof", &mus_cm_ndof, &b_mus_cm_ndof);
   fChain->SetBranchAddress("mus_cm_chg", &mus_cm_chg, &b_mus_cm_chg);
   fChain->SetBranchAddress("mus_cm_pt", &mus_cm_pt, &b_mus_cm_pt);
   fChain->SetBranchAddress("mus_cm_px", &mus_cm_px, &b_mus_cm_px);
   fChain->SetBranchAddress("mus_cm_py", &mus_cm_py, &b_mus_cm_py);
   fChain->SetBranchAddress("mus_cm_pz", &mus_cm_pz, &b_mus_cm_pz);
   fChain->SetBranchAddress("mus_cm_eta", &mus_cm_eta, &b_mus_cm_eta);
   fChain->SetBranchAddress("mus_cm_phi", &mus_cm_phi, &b_mus_cm_phi);
   fChain->SetBranchAddress("mus_cm_theta", &mus_cm_theta, &b_mus_cm_theta);
   fChain->SetBranchAddress("mus_cm_d0dum", &mus_cm_d0dum, &b_mus_cm_d0dum);
   fChain->SetBranchAddress("mus_cm_dz", &mus_cm_dz, &b_mus_cm_dz);
   fChain->SetBranchAddress("mus_cm_vx", &mus_cm_vx, &b_mus_cm_vx);
   fChain->SetBranchAddress("mus_cm_vy", &mus_cm_vy, &b_mus_cm_vy);
   fChain->SetBranchAddress("mus_cm_vz", &mus_cm_vz, &b_mus_cm_vz);
   fChain->SetBranchAddress("mus_cm_numvalhits", &mus_cm_numvalhits, &b_mus_cm_numvalhits);
   fChain->SetBranchAddress("mus_cm_numlosthits", &mus_cm_numlosthits, &b_mus_cm_numlosthits);
   fChain->SetBranchAddress("mus_cm_d0dumErr", &mus_cm_d0dumErr, &b_mus_cm_d0dumErr);
   fChain->SetBranchAddress("mus_cm_dzErr", &mus_cm_dzErr, &b_mus_cm_dzErr);
   fChain->SetBranchAddress("mus_cm_ptErr", &mus_cm_ptErr, &b_mus_cm_ptErr);
   fChain->SetBranchAddress("mus_cm_etaErr", &mus_cm_etaErr, &b_mus_cm_etaErr);
   fChain->SetBranchAddress("mus_cm_phiErr", &mus_cm_phiErr, &b_mus_cm_phiErr);
   fChain->SetBranchAddress("mus_tk_id", &mus_tk_id, &b_mus_tk_id);
   fChain->SetBranchAddress("mus_tk_chi2", &mus_tk_chi2, &b_mus_tk_chi2);
   fChain->SetBranchAddress("mus_tk_ndof", &mus_tk_ndof, &b_mus_tk_ndof);
   fChain->SetBranchAddress("mus_tk_chg", &mus_tk_chg, &b_mus_tk_chg);
   fChain->SetBranchAddress("mus_tk_pt", &mus_tk_pt, &b_mus_tk_pt);
   fChain->SetBranchAddress("mus_tk_px", &mus_tk_px, &b_mus_tk_px);
   fChain->SetBranchAddress("mus_tk_py", &mus_tk_py, &b_mus_tk_py);
   fChain->SetBranchAddress("mus_tk_pz", &mus_tk_pz, &b_mus_tk_pz);
   fChain->SetBranchAddress("mus_tk_eta", &mus_tk_eta, &b_mus_tk_eta);
   fChain->SetBranchAddress("mus_tk_phi", &mus_tk_phi, &b_mus_tk_phi);
   fChain->SetBranchAddress("mus_tk_theta", &mus_tk_theta, &b_mus_tk_theta);
   fChain->SetBranchAddress("mus_tk_d0dum", &mus_tk_d0dum, &b_mus_tk_d0dum);
   fChain->SetBranchAddress("mus_tk_dz", &mus_tk_dz, &b_mus_tk_dz);
   fChain->SetBranchAddress("mus_tk_vx", &mus_tk_vx, &b_mus_tk_vx);
   fChain->SetBranchAddress("mus_tk_vy", &mus_tk_vy, &b_mus_tk_vy);
   fChain->SetBranchAddress("mus_tk_vz", &mus_tk_vz, &b_mus_tk_vz);
   fChain->SetBranchAddress("mus_tk_numvalhits", &mus_tk_numvalhits, &b_mus_tk_numvalhits);
   fChain->SetBranchAddress("mus_tk_numlosthits", &mus_tk_numlosthits, &b_mus_tk_numlosthits);
   fChain->SetBranchAddress("mus_tk_d0dumErr", &mus_tk_d0dumErr, &b_mus_tk_d0dumErr);
   fChain->SetBranchAddress("mus_tk_dzErr", &mus_tk_dzErr, &b_mus_tk_dzErr);
   fChain->SetBranchAddress("mus_tk_ptErr", &mus_tk_ptErr, &b_mus_tk_ptErr);
   fChain->SetBranchAddress("mus_tk_etaErr", &mus_tk_etaErr, &b_mus_tk_etaErr);
   fChain->SetBranchAddress("mus_tk_phiErr", &mus_tk_phiErr, &b_mus_tk_phiErr);
   fChain->SetBranchAddress("mus_stamu_chi2", &mus_stamu_chi2, &b_mus_stamu_chi2);
   fChain->SetBranchAddress("mus_stamu_ndof", &mus_stamu_ndof, &b_mus_stamu_ndof);
   fChain->SetBranchAddress("mus_stamu_chg", &mus_stamu_chg, &b_mus_stamu_chg);
   fChain->SetBranchAddress("mus_stamu_pt", &mus_stamu_pt, &b_mus_stamu_pt);
   fChain->SetBranchAddress("mus_stamu_px", &mus_stamu_px, &b_mus_stamu_px);
   fChain->SetBranchAddress("mus_stamu_py", &mus_stamu_py, &b_mus_stamu_py);
   fChain->SetBranchAddress("mus_stamu_pz", &mus_stamu_pz, &b_mus_stamu_pz);
   fChain->SetBranchAddress("mus_stamu_eta", &mus_stamu_eta, &b_mus_stamu_eta);
   fChain->SetBranchAddress("mus_stamu_phi", &mus_stamu_phi, &b_mus_stamu_phi);
   fChain->SetBranchAddress("mus_stamu_theta", &mus_stamu_theta, &b_mus_stamu_theta);
   fChain->SetBranchAddress("mus_stamu_d0dum", &mus_stamu_d0dum, &b_mus_stamu_d0dum);
   fChain->SetBranchAddress("mus_stamu_dz", &mus_stamu_dz, &b_mus_stamu_dz);
   fChain->SetBranchAddress("mus_stamu_vx", &mus_stamu_vx, &b_mus_stamu_vx);
   fChain->SetBranchAddress("mus_stamu_vy", &mus_stamu_vy, &b_mus_stamu_vy);
   fChain->SetBranchAddress("mus_stamu_vz", &mus_stamu_vz, &b_mus_stamu_vz);
   fChain->SetBranchAddress("mus_stamu_numvalhits", &mus_stamu_numvalhits, &b_mus_stamu_numvalhits);
   fChain->SetBranchAddress("mus_stamu_numlosthits", &mus_stamu_numlosthits, &b_mus_stamu_numlosthits);
   fChain->SetBranchAddress("mus_stamu_d0dumErr", &mus_stamu_d0dumErr, &b_mus_stamu_d0dumErr);
   fChain->SetBranchAddress("mus_stamu_dzErr", &mus_stamu_dzErr, &b_mus_stamu_dzErr);
   fChain->SetBranchAddress("mus_stamu_ptErr", &mus_stamu_ptErr, &b_mus_stamu_ptErr);
   fChain->SetBranchAddress("mus_stamu_etaErr", &mus_stamu_etaErr, &b_mus_stamu_etaErr);
   fChain->SetBranchAddress("mus_stamu_phiErr", &mus_stamu_phiErr, &b_mus_stamu_phiErr);
   fChain->SetBranchAddress("mus_num_matches", &mus_num_matches, &b_mus_num_matches);
   fChain->SetBranchAddress("mus_id_All", &mus_id_All, &b_mus_id_All);
   fChain->SetBranchAddress("mus_id_AllGlobalMuons", &mus_id_AllGlobalMuons, &b_mus_id_AllGlobalMuons);
   fChain->SetBranchAddress("mus_id_AllStandAloneMuons", &mus_id_AllStandAloneMuons, &b_mus_id_AllStandAloneMuons);
   fChain->SetBranchAddress("mus_id_AllTrackerMuons", &mus_id_AllTrackerMuons, &b_mus_id_AllTrackerMuons);
   fChain->SetBranchAddress("mus_id_TrackerMuonArbitrated", &mus_id_TrackerMuonArbitrated, &b_mus_id_TrackerMuonArbitrated);
   fChain->SetBranchAddress("mus_id_AllArbitrated", &mus_id_AllArbitrated, &b_mus_id_AllArbitrated);
   fChain->SetBranchAddress("mus_id_GlobalMuonPromptTight", &mus_id_GlobalMuonPromptTight, &b_mus_id_GlobalMuonPromptTight);
   fChain->SetBranchAddress("mus_id_TMLastStationLoose", &mus_id_TMLastStationLoose, &b_mus_id_TMLastStationLoose);
   fChain->SetBranchAddress("mus_id_TMLastStationTight", &mus_id_TMLastStationTight, &b_mus_id_TMLastStationTight);
   fChain->SetBranchAddress("mus_id_TM2DCompatibilityLoose", &mus_id_TM2DCompatibilityLoose, &b_mus_id_TM2DCompatibilityLoose);
   fChain->SetBranchAddress("mus_id_TM2DCompatibilityTight", &mus_id_TM2DCompatibilityTight, &b_mus_id_TM2DCompatibilityTight);
   fChain->SetBranchAddress("mus_id_TMOneStationLoose", &mus_id_TMOneStationLoose, &b_mus_id_TMOneStationLoose);
   fChain->SetBranchAddress("mus_id_TMOneStationTight", &mus_id_TMOneStationTight, &b_mus_id_TMOneStationTight);
   fChain->SetBranchAddress("mus_id_TMLastStationOptimizedLowPtLoose", &mus_id_TMLastStationOptimizedLowPtLoose, &b_mus_id_TMLastStationOptimizedLowPtLoose);
   fChain->SetBranchAddress("mus_id_TMLastStationOptimizedLowPtTight", &mus_id_TMLastStationOptimizedLowPtTight, &b_mus_id_TMLastStationOptimizedLowPtTight);
   fChain->SetBranchAddress("Nphotons", &Nphotons, &b_Nphotons);
   fChain->SetBranchAddress("photons_energy", &photons_energy, &b_photons_energy);
   fChain->SetBranchAddress("photons_et", &photons_et, &b_photons_et);
   fChain->SetBranchAddress("photons_eta", &photons_eta, &b_photons_eta);
   fChain->SetBranchAddress("photons_phi", &photons_phi, &b_photons_phi);
   fChain->SetBranchAddress("photons_pt", &photons_pt, &b_photons_pt);
   fChain->SetBranchAddress("photons_px", &photons_px, &b_photons_px);
   fChain->SetBranchAddress("photons_py", &photons_py, &b_photons_py);
   fChain->SetBranchAddress("photons_pz", &photons_pz, &b_photons_pz);
   fChain->SetBranchAddress("photons_status", &photons_status, &b_photons_status);
   fChain->SetBranchAddress("photons_theta", &photons_theta, &b_photons_theta);
   fChain->SetBranchAddress("photons_hadOverEM", &photons_hadOverEM, &b_photons_hadOverEM);
   fChain->SetBranchAddress("photons_scEnergy", &photons_scEnergy, &b_photons_scEnergy);
   fChain->SetBranchAddress("photons_scRawEnergy", &photons_scRawEnergy, &b_photons_scRawEnergy);
   fChain->SetBranchAddress("photons_scEta", &photons_scEta, &b_photons_scEta);
   fChain->SetBranchAddress("photons_scPhi", &photons_scPhi, &b_photons_scPhi);
   fChain->SetBranchAddress("photons_scEtaWidth", &photons_scEtaWidth, &b_photons_scEtaWidth);
   fChain->SetBranchAddress("photons_scPhiWidth", &photons_scPhiWidth, &b_photons_scPhiWidth);
   fChain->SetBranchAddress("photons_tIso", &photons_tIso, &b_photons_tIso);
   fChain->SetBranchAddress("photons_ecalIso", &photons_ecalIso, &b_photons_ecalIso);
   fChain->SetBranchAddress("photons_hcalIso", &photons_hcalIso, &b_photons_hcalIso);
   fChain->SetBranchAddress("photons_isoEcalRecHitDR04", &photons_isoEcalRecHitDR04, &b_photons_isoEcalRecHitDR04);
   fChain->SetBranchAddress("photons_isoHcalRecHitDR04", &photons_isoHcalRecHitDR04, &b_photons_isoHcalRecHitDR04);
   fChain->SetBranchAddress("photons_isoSolidTrkConeDR04", &photons_isoSolidTrkConeDR04, &b_photons_isoSolidTrkConeDR04);
   fChain->SetBranchAddress("photons_isoHollowTrkConeDR04", &photons_isoHollowTrkConeDR04, &b_photons_isoHollowTrkConeDR04);
   fChain->SetBranchAddress("photons_nTrkSolidConeDR04", &photons_nTrkSolidConeDR04, &b_photons_nTrkSolidConeDR04);
   fChain->SetBranchAddress("photons_nTrkHollowConeDR04", &photons_nTrkHollowConeDR04, &b_photons_nTrkHollowConeDR04);
   fChain->SetBranchAddress("photons_isoEcalRecHitDR03", &photons_isoEcalRecHitDR03, &b_photons_isoEcalRecHitDR03);
   fChain->SetBranchAddress("photons_isoHcalRecHitDR03", &photons_isoHcalRecHitDR03, &b_photons_isoHcalRecHitDR03);
   fChain->SetBranchAddress("photons_isoSolidTrkConeDR03", &photons_isoSolidTrkConeDR03, &b_photons_isoSolidTrkConeDR03);
   fChain->SetBranchAddress("photons_isoHollowTrkConeDR03", &photons_isoHollowTrkConeDR03, &b_photons_isoHollowTrkConeDR03);
   fChain->SetBranchAddress("photons_nTrkSolidConeDR03", &photons_nTrkSolidConeDR03, &b_photons_nTrkSolidConeDR03);
   fChain->SetBranchAddress("photons_nTrkHollowConeDR03", &photons_nTrkHollowConeDR03, &b_photons_nTrkHollowConeDR03);
   fChain->SetBranchAddress("photons_isAlsoElectron", &photons_isAlsoElectron, &b_photons_isAlsoElectron);
   fChain->SetBranchAddress("photons_hasPixelSeed", &photons_hasPixelSeed, &b_photons_hasPixelSeed);
   fChain->SetBranchAddress("photons_isConverted", &photons_isConverted, &b_photons_isConverted);
   fChain->SetBranchAddress("photons_isEBGap", &photons_isEBGap, &b_photons_isEBGap);
   fChain->SetBranchAddress("photons_isEEGap", &photons_isEEGap, &b_photons_isEEGap);
   fChain->SetBranchAddress("photons_isEBEEGap", &photons_isEBEEGap, &b_photons_isEBEEGap);
   fChain->SetBranchAddress("photons_isEBPho", &photons_isEBPho, &b_photons_isEBPho);
   fChain->SetBranchAddress("photons_isEEPho", &photons_isEEPho, &b_photons_isEEPho);
   fChain->SetBranchAddress("photons_isLoosePhoton", &photons_isLoosePhoton, &b_photons_isLoosePhoton);
   fChain->SetBranchAddress("photons_isTightPhoton", &photons_isTightPhoton, &b_photons_isTightPhoton);
   fChain->SetBranchAddress("photons_r9", &photons_r9, &b_photons_r9);
   fChain->SetBranchAddress("photons_gen_et", &photons_gen_et, &b_photons_gen_et);
   fChain->SetBranchAddress("photons_gen_eta", &photons_gen_eta, &b_photons_gen_eta);
   fChain->SetBranchAddress("photons_gen_phi", &photons_gen_phi, &b_photons_gen_phi);
   fChain->SetBranchAddress("photons_gen_id", &photons_gen_id, &b_photons_gen_id);
   fChain->SetBranchAddress("Npv", &Npv, &b_Npv);
   fChain->SetBranchAddress("pv_x", &pv_x, &b_pv_x);
   fChain->SetBranchAddress("pv_y", &pv_y, &b_pv_y);
   fChain->SetBranchAddress("pv_z", &pv_z, &b_pv_z);
   fChain->SetBranchAddress("pv_xErr", &pv_xErr, &b_pv_xErr);
   fChain->SetBranchAddress("pv_yErr", &pv_yErr, &b_pv_yErr);
   fChain->SetBranchAddress("pv_zErr", &pv_zErr, &b_pv_zErr);
   fChain->SetBranchAddress("Ntaus", &Ntaus, &b_Ntaus);
   fChain->SetBranchAddress("taus_energy", &taus_energy, &b_taus_energy);
   fChain->SetBranchAddress("taus_et", &taus_et, &b_taus_et);
   fChain->SetBranchAddress("taus_eta", &taus_eta, &b_taus_eta);
   fChain->SetBranchAddress("taus_phi", &taus_phi, &b_taus_phi);
   fChain->SetBranchAddress("taus_pt", &taus_pt, &b_taus_pt);
   fChain->SetBranchAddress("taus_px", &taus_px, &b_taus_px);
   fChain->SetBranchAddress("taus_py", &taus_py, &b_taus_py);
   fChain->SetBranchAddress("taus_pz", &taus_pz, &b_taus_pz);
   fChain->SetBranchAddress("taus_status", &taus_status, &b_taus_status);
   fChain->SetBranchAddress("taus_theta", &taus_theta, &b_taus_theta);
   fChain->SetBranchAddress("taus_charge", &taus_charge, &b_taus_charge);
   fChain->SetBranchAddress("taus_emf", &taus_emf, &b_taus_emf);
   fChain->SetBranchAddress("taus_hcalTotOverPLead", &taus_hcalTotOverPLead, &b_taus_hcalTotOverPLead);
   fChain->SetBranchAddress("taus_hcalMaxOverPLead", &taus_hcalMaxOverPLead, &b_taus_hcalMaxOverPLead);
   fChain->SetBranchAddress("taus_hcal3x3OverPLead", &taus_hcal3x3OverPLead, &b_taus_hcal3x3OverPLead);
   fChain->SetBranchAddress("taus_ecalStripSumEOverPLead", &taus_ecalStripSumEOverPLead, &b_taus_ecalStripSumEOverPLead);
   fChain->SetBranchAddress("taus_elecPreIdOutput", &taus_elecPreIdOutput, &b_taus_elecPreIdOutput);
   fChain->SetBranchAddress("taus_elecPreIdDecision", &taus_elecPreIdDecision, &b_taus_elecPreIdDecision);
   fChain->SetBranchAddress("taus_leadPFChargedHadrCand_pt", &taus_leadPFChargedHadrCand_pt, &b_taus_leadPFChargedHadrCand_pt);
   fChain->SetBranchAddress("taus_leadPFChargedHadrCand_charge", &taus_leadPFChargedHadrCand_charge, &b_taus_leadPFChargedHadrCand_charge);
   fChain->SetBranchAddress("taus_leadPFChargedHadrCand_eta", &taus_leadPFChargedHadrCand_eta, &b_taus_leadPFChargedHadrCand_eta);
   fChain->SetBranchAddress("taus_leadPFChargedHadrCand_ECAL_eta", &taus_leadPFChargedHadrCand_ECAL_eta, &b_taus_leadPFChargedHadrCand_ECAL_eta);
   fChain->SetBranchAddress("taus_leadPFChargedHadrCand_phi", &taus_leadPFChargedHadrCand_phi, &b_taus_leadPFChargedHadrCand_phi);
   fChain->SetBranchAddress("taus_isoPFGammaCandsEtSum", &taus_isoPFGammaCandsEtSum, &b_taus_isoPFGammaCandsEtSum);
   fChain->SetBranchAddress("taus_isoPFChargedHadrCandsPtSum", &taus_isoPFChargedHadrCandsPtSum, &b_taus_isoPFChargedHadrCandsPtSum);
   fChain->SetBranchAddress("taus_leadingTrackFinding", &taus_leadingTrackFinding, &b_taus_leadingTrackFinding);
   fChain->SetBranchAddress("taus_leadingTrackPtCut", &taus_leadingTrackPtCut, &b_taus_leadingTrackPtCut);
   fChain->SetBranchAddress("taus_trackIsolation", &taus_trackIsolation, &b_taus_trackIsolation);
   fChain->SetBranchAddress("taus_ecalIsolation", &taus_ecalIsolation, &b_taus_ecalIsolation);
   fChain->SetBranchAddress("taus_byIsolation", &taus_byIsolation, &b_taus_byIsolation);
   fChain->SetBranchAddress("taus_againstElectron", &taus_againstElectron, &b_taus_againstElectron);
   fChain->SetBranchAddress("taus_againstMuon", &taus_againstMuon, &b_taus_againstMuon);
   fChain->SetBranchAddress("taus_taNC_quarter", &taus_taNC_quarter, &b_taus_taNC_quarter);
   fChain->SetBranchAddress("taus_taNC_one", &taus_taNC_one, &b_taus_taNC_one);
   fChain->SetBranchAddress("taus_taNC_half", &taus_taNC_half, &b_taus_taNC_half);
   fChain->SetBranchAddress("taus_taNC_tenth", &taus_taNC_tenth, &b_taus_taNC_tenth);
   fChain->SetBranchAddress("taus_muDecision", &taus_muDecision, &b_taus_muDecision);
   fChain->SetBranchAddress("taus_Nprongs", &taus_Nprongs, &b_taus_Nprongs);
   fChain->SetBranchAddress("Ntcmets", &Ntcmets, &b_Ntcmets);
   fChain->SetBranchAddress("tcmets_et", &tcmets_et, &b_tcmets_et);
   fChain->SetBranchAddress("tcmets_phi", &tcmets_phi, &b_tcmets_phi);
   fChain->SetBranchAddress("tcmets_ex", &tcmets_ex, &b_tcmets_ex);
   fChain->SetBranchAddress("tcmets_ey", &tcmets_ey, &b_tcmets_ey);
   fChain->SetBranchAddress("tcmets_sumEt", &tcmets_sumEt, &b_tcmets_sumEt);
   fChain->SetBranchAddress("Ntracks", &Ntracks, &b_Ntracks);
   fChain->SetBranchAddress("tracks_chi2", &tracks_chi2, &b_tracks_chi2);
   fChain->SetBranchAddress("tracks_ndof", &tracks_ndof, &b_tracks_ndof);
   fChain->SetBranchAddress("tracks_chg", &tracks_chg, &b_tracks_chg);
   fChain->SetBranchAddress("tracks_pt", &tracks_pt, &b_tracks_pt);
   fChain->SetBranchAddress("tracks_px", &tracks_px, &b_tracks_px);
   fChain->SetBranchAddress("tracks_py", &tracks_py, &b_tracks_py);
   fChain->SetBranchAddress("tracks_pz", &tracks_pz, &b_tracks_pz);
   fChain->SetBranchAddress("tracks_eta", &tracks_eta, &b_tracks_eta);
   fChain->SetBranchAddress("tracks_phi", &tracks_phi, &b_tracks_phi);
   fChain->SetBranchAddress("tracks_theta", &tracks_theta, &b_tracks_theta);
   fChain->SetBranchAddress("tracks_d0dum", &tracks_d0dum, &b_tracks_d0dum);
   fChain->SetBranchAddress("tracks_dz", &tracks_dz, &b_tracks_dz);
   fChain->SetBranchAddress("tracks_vx", &tracks_vx, &b_tracks_vx);
   fChain->SetBranchAddress("tracks_vy", &tracks_vy, &b_tracks_vy);
   fChain->SetBranchAddress("tracks_vz", &tracks_vz, &b_tracks_vz);
   fChain->SetBranchAddress("tracks_numvalhits", &tracks_numvalhits, &b_tracks_numvalhits);
   fChain->SetBranchAddress("tracks_numlosthits", &tracks_numlosthits, &b_tracks_numlosthits);
   fChain->SetBranchAddress("tracks_d0dumErr", &tracks_d0dumErr, &b_tracks_d0dumErr);
   fChain->SetBranchAddress("tracks_dzErr", &tracks_dzErr, &b_tracks_dzErr);
   fChain->SetBranchAddress("tracks_ptErr", &tracks_ptErr, &b_tracks_ptErr);
   fChain->SetBranchAddress("tracks_etaErr", &tracks_etaErr, &b_tracks_etaErr);
   fChain->SetBranchAddress("tracks_phiErr", &tracks_phiErr, &b_tracks_phiErr);
   fChain->SetBranchAddress("tracks_Nrechits", &tracks_Nrechits, &b_tracks_Nrechits);
   fChain->SetBranchAddress("tracks_innerHitX", &tracks_innerHitX, &b_tracks_innerHitX);
   fChain->SetBranchAddress("tracks_innerHitY", &tracks_innerHitY, &b_tracks_innerHitY);
   fChain->SetBranchAddress("tracks_innerHitZ", &tracks_innerHitZ, &b_tracks_innerHitZ);
   fChain->SetBranchAddress("tracks_outerHitX", &tracks_outerHitX, &b_tracks_outerHitX);
   fChain->SetBranchAddress("tracks_outerHitY", &tracks_outerHitY, &b_tracks_outerHitY);
   fChain->SetBranchAddress("tracks_outerHitZ", &tracks_outerHitZ, &b_tracks_outerHitZ);
   fChain->SetBranchAddress("tracks_highPurity", &tracks_highPurity, &b_tracks_highPurity);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);

}

   //t2
void InitializeV(TChain *fChain)
{
   fChain->SetBranchAddress("AlCa_EcalEta_8E29", &AlCa_EcalEta_8E29, &b_AlCa_EcalEta_8E29);
   fChain->SetBranchAddress("AlCa_EcalPhiSym", &AlCa_EcalPhiSym, &b_AlCa_EcalPhiSym);
   fChain->SetBranchAddress("AlCa_EcalPi0_8E29", &AlCa_EcalPi0_8E29, &b_AlCa_EcalPi0_8E29);
   fChain->SetBranchAddress("AlCa_HcalPhiSym", &AlCa_HcalPhiSym, &b_AlCa_HcalPhiSym);
   fChain->SetBranchAddress("AlCa_RPCMuonNoHits", &AlCa_RPCMuonNoHits, &b_AlCa_RPCMuonNoHits);
   fChain->SetBranchAddress("AlCa_RPCMuonNormalisation", &AlCa_RPCMuonNormalisation, &b_AlCa_RPCMuonNormalisation);
   fChain->SetBranchAddress("HLTAnalyzerEndpath", &HLTAnalyzerEndpath, &b_HLTAnalyzerEndpath);
   fChain->SetBranchAddress("HLT_BTagIP_Jet50U", &HLT_BTagIP_Jet50U, &b_HLT_BTagIP_Jet50U);
   fChain->SetBranchAddress("HLT_BTagMu_Jet10U", &HLT_BTagMu_Jet10U, &b_HLT_BTagMu_Jet10U);
   fChain->SetBranchAddress("HLT_BackwardBSC", &HLT_BackwardBSC, &b_HLT_BackwardBSC);
   fChain->SetBranchAddress("HLT_CSCBeamHalo", &HLT_CSCBeamHalo, &b_HLT_CSCBeamHalo);
   fChain->SetBranchAddress("HLT_CSCBeamHaloOverlapRing1", &HLT_CSCBeamHaloOverlapRing1, &b_HLT_CSCBeamHaloOverlapRing1);
   fChain->SetBranchAddress("HLT_CSCBeamHaloOverlapRing2", &HLT_CSCBeamHaloOverlapRing2, &b_HLT_CSCBeamHaloOverlapRing2);
   fChain->SetBranchAddress("HLT_CSCBeamHaloRing2or3", &HLT_CSCBeamHaloRing2or3, &b_HLT_CSCBeamHaloRing2or3);
   fChain->SetBranchAddress("HLT_DiJetAve15U_8E29", &HLT_DiJetAve15U_8E29, &b_HLT_DiJetAve15U_8E29);
   fChain->SetBranchAddress("HLT_DiJetAve30U_8E29", &HLT_DiJetAve30U_8E29, &b_HLT_DiJetAve30U_8E29);
   fChain->SetBranchAddress("HLT_DoubleEle5_SW_L1R", &HLT_DoubleEle5_SW_L1R, &b_HLT_DoubleEle5_SW_L1R);
   fChain->SetBranchAddress("HLT_DoubleLooseIsoTau15", &HLT_DoubleLooseIsoTau15, &b_HLT_DoubleLooseIsoTau15);
   fChain->SetBranchAddress("HLT_DoubleMu0", &HLT_DoubleMu0, &b_HLT_DoubleMu0);
   fChain->SetBranchAddress("HLT_DoubleMu3", &HLT_DoubleMu3, &b_HLT_DoubleMu3);
   fChain->SetBranchAddress("HLT_DoublePhoton10_L1R", &HLT_DoublePhoton10_L1R, &b_HLT_DoublePhoton10_L1R);
   fChain->SetBranchAddress("HLT_DoublePhoton5_Jpsi_L1R", &HLT_DoublePhoton5_Jpsi_L1R, &b_HLT_DoublePhoton5_Jpsi_L1R);
   fChain->SetBranchAddress("HLT_DoublePhoton5_Upsilon_L1R", &HLT_DoublePhoton5_Upsilon_L1R, &b_HLT_DoublePhoton5_Upsilon_L1R);
   fChain->SetBranchAddress("HLT_DoublePhoton5_eeRes_L1R", &HLT_DoublePhoton5_eeRes_L1R, &b_HLT_DoublePhoton5_eeRes_L1R);
   fChain->SetBranchAddress("HLT_Ele10_LW_EleId_L1R", &HLT_Ele10_LW_EleId_L1R, &b_HLT_Ele10_LW_EleId_L1R);
   fChain->SetBranchAddress("HLT_Ele10_LW_L1R", &HLT_Ele10_LW_L1R, &b_HLT_Ele10_LW_L1R);
   fChain->SetBranchAddress("HLT_Ele15_LW_L1R", &HLT_Ele15_LW_L1R, &b_HLT_Ele15_LW_L1R);
   fChain->SetBranchAddress("HLT_Ele15_SC10_LW_L1R", &HLT_Ele15_SC10_LW_L1R, &b_HLT_Ele15_SC10_LW_L1R);
   fChain->SetBranchAddress("HLT_Ele15_SiStrip_L1R", &HLT_Ele15_SiStrip_L1R, &b_HLT_Ele15_SiStrip_L1R);
   fChain->SetBranchAddress("HLT_Ele20_LW_L1R", &HLT_Ele20_LW_L1R, &b_HLT_Ele20_LW_L1R);
   fChain->SetBranchAddress("HLT_ForwardBSC", &HLT_ForwardBSC, &b_HLT_ForwardBSC);
   fChain->SetBranchAddress("HLT_FwdJet20U", &HLT_FwdJet20U, &b_HLT_FwdJet20U);
   fChain->SetBranchAddress("HLT_HT100U", &HLT_HT100U, &b_HLT_HT100U);
   fChain->SetBranchAddress("HLT_IsoMu3", &HLT_IsoMu3, &b_HLT_IsoMu3);
   fChain->SetBranchAddress("HLT_IsoTrack_8E29", &HLT_IsoTrack_8E29, &b_HLT_IsoTrack_8E29);
   fChain->SetBranchAddress("HLT_Jet15U", &HLT_Jet15U, &b_HLT_Jet15U);
   fChain->SetBranchAddress("HLT_Jet30U", &HLT_Jet30U, &b_HLT_Jet30U);
   fChain->SetBranchAddress("HLT_Jet50U", &HLT_Jet50U, &b_HLT_Jet50U);
   fChain->SetBranchAddress("HLT_L1DoubleEG5", &HLT_L1DoubleEG5, &b_HLT_L1DoubleEG5);
   fChain->SetBranchAddress("HLT_L1DoubleMuOpen", &HLT_L1DoubleMuOpen, &b_HLT_L1DoubleMuOpen);
   fChain->SetBranchAddress("HLT_L1Jet6U", &HLT_L1Jet6U, &b_HLT_L1Jet6U);
   fChain->SetBranchAddress("HLT_L1MET20", &HLT_L1MET20, &b_HLT_L1MET20);
   fChain->SetBranchAddress("HLT_L1Mu", &HLT_L1Mu, &b_HLT_L1Mu);
   fChain->SetBranchAddress("HLT_L1Mu14_L1ETM30", &HLT_L1Mu14_L1ETM30, &b_HLT_L1Mu14_L1ETM30);
   fChain->SetBranchAddress("HLT_L1Mu14_L1SingleEG10", &HLT_L1Mu14_L1SingleEG10, &b_HLT_L1Mu14_L1SingleEG10);
   fChain->SetBranchAddress("HLT_L1Mu14_L1SingleJet6U", &HLT_L1Mu14_L1SingleJet6U, &b_HLT_L1Mu14_L1SingleJet6U);
   fChain->SetBranchAddress("HLT_L1Mu20", &HLT_L1Mu20, &b_HLT_L1Mu20);
   fChain->SetBranchAddress("HLT_L1MuOpen", &HLT_L1MuOpen, &b_HLT_L1MuOpen);
   fChain->SetBranchAddress("HLT_L1SingleEG5", &HLT_L1SingleEG5, &b_HLT_L1SingleEG5);
   fChain->SetBranchAddress("HLT_L1SingleEG8", &HLT_L1SingleEG8, &b_HLT_L1SingleEG8);
   fChain->SetBranchAddress("HLT_L2Mu11", &HLT_L2Mu11, &b_HLT_L2Mu11);
   fChain->SetBranchAddress("HLT_L2Mu9", &HLT_L2Mu9, &b_HLT_L2Mu9);
   fChain->SetBranchAddress("HLT_MET100", &HLT_MET100, &b_HLT_MET100);
   fChain->SetBranchAddress("HLT_MET45", &HLT_MET45, &b_HLT_MET45);
   fChain->SetBranchAddress("HLT_MinBiasEcal", &HLT_MinBiasEcal, &b_HLT_MinBiasEcal);
   fChain->SetBranchAddress("HLT_MinBiasHcal", &HLT_MinBiasHcal, &b_HLT_MinBiasHcal);
   fChain->SetBranchAddress("HLT_MinBiasPixel", &HLT_MinBiasPixel, &b_HLT_MinBiasPixel);
   fChain->SetBranchAddress("HLT_MinBiasPixel_Trk5", &HLT_MinBiasPixel_Trk5, &b_HLT_MinBiasPixel_Trk5);
   fChain->SetBranchAddress("HLT_Mu3", &HLT_Mu3, &b_HLT_Mu3);
   fChain->SetBranchAddress("HLT_Mu5", &HLT_Mu5, &b_HLT_Mu5);
   fChain->SetBranchAddress("HLT_Mu9", &HLT_Mu9, &b_HLT_Mu9);
   fChain->SetBranchAddress("HLT_Photon10_L1R", &HLT_Photon10_L1R, &b_HLT_Photon10_L1R);
   fChain->SetBranchAddress("HLT_Photon15_L1R", &HLT_Photon15_L1R, &b_HLT_Photon15_L1R);
   fChain->SetBranchAddress("HLT_Photon15_LooseEcalIso_L1R", &HLT_Photon15_LooseEcalIso_L1R, &b_HLT_Photon15_LooseEcalIso_L1R);
   fChain->SetBranchAddress("HLT_Photon15_TrackIso_L1R", &HLT_Photon15_TrackIso_L1R, &b_HLT_Photon15_TrackIso_L1R);
   fChain->SetBranchAddress("HLT_Photon20_L1R", &HLT_Photon20_L1R, &b_HLT_Photon20_L1R);
   fChain->SetBranchAddress("HLT_Photon30_L1R_8E29", &HLT_Photon30_L1R_8E29, &b_HLT_Photon30_L1R_8E29);
   fChain->SetBranchAddress("HLT_QuadJet15U", &HLT_QuadJet15U, &b_HLT_QuadJet15U);
   fChain->SetBranchAddress("HLT_SingleLooseIsoTau20", &HLT_SingleLooseIsoTau20, &b_HLT_SingleLooseIsoTau20);
   fChain->SetBranchAddress("HLT_StoppedHSCP_8E29", &HLT_StoppedHSCP_8E29, &b_HLT_StoppedHSCP_8E29);
   fChain->SetBranchAddress("HLT_TrackerCosmics", &HLT_TrackerCosmics, &b_HLT_TrackerCosmics);
   fChain->SetBranchAddress("HLT_ZeroBias", &HLT_ZeroBias, &b_HLT_ZeroBias);
   fChain->SetBranchAddress("HLTriggerFinalPath", &HLTriggerFinalPath, &b_HLTriggerFinalPath);
   fChain->SetBranchAddress("HLTriggerFirstPath", &HLTriggerFirstPath, &b_HLTriggerFirstPath);

}
