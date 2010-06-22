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


//Version 24 of cfA

   // Declaration of leaf types
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
   UInt_t          NconvTk;
   vector<float>   *convTk_tracksPin_x1;
   vector<float>   *convTk_tracksPin_x2;
   vector<float>   *convTk_tracksPin_y1;
   vector<float>   *convTk_tracksPin_y2;
   vector<float>   *convTk_tracksPin_z1;
   vector<float>   *convTk_tracksPin_z2;
   vector<float>   *convTk_tracksPout_x1;
   vector<float>   *convTk_tracksPout_x2;
   vector<float>   *convTk_tracksPout_y1;
   vector<float>   *convTk_tracksPout_y2;
   vector<float>   *convTk_tracksPout_z1;
   vector<float>   *convTk_tracksPout_z2;
   vector<float>   *convTk_ecalImpactPosition_x1;
   vector<float>   *convTk_ecalImpactPosition_x2;
   vector<float>   *convTk_ecalImpactPosition_y1;
   vector<float>   *convTk_ecalImpactPosition_y2;
   vector<float>   *convTk_ecalImpactPosition_z1;
   vector<float>   *convTk_ecalImpactPosition_z2;
   vector<float>   *convTk_tracksInnerPosition_x1;
   vector<float>   *convTk_tracksInnerPosition_x2;
   vector<float>   *convTk_tracksInnerPosition_y1;
   vector<float>   *convTk_tracksInnerPosition_y2;
   vector<float>   *convTk_tracksInnerPosition_z1;
   vector<float>   *convTk_tracksInnerPosition_z2;
   vector<float>   *convTk_tracks_px1;
   vector<float>   *convTk_tracks_py1;
   vector<float>   *convTk_tracks_pz1;
   vector<float>   *convTk_tracks_phi1;
   vector<float>   *convTk_tracks_d01;
   vector<float>   *convTk_tracks_dz1;
   vector<float>   *convTk_tracks_charge1;
   vector<float>   *convTk_tracks_px2;
   vector<float>   *convTk_tracks_py2;
   vector<float>   *convTk_tracks_pz2;
   vector<float>   *convTk_tracks_phi2;
   vector<float>   *convTk_tracks_d02;
   vector<float>   *convTk_tracks_dz2;
   vector<float>   *convTk_tracks_charge2;
   vector<float>   *convTk_bcMatchingWithTracks_x1;
   vector<float>   *convTk_bcMatchingWithTracks_y1;
   vector<float>   *convTk_bcMatchingWithTracks_z1;
   vector<float>   *convTk_bcMatchingWithTracks_x2;
   vector<float>   *convTk_bcMatchingWithTracks_y2;
   vector<float>   *convTk_bcMatchingWithTracks_z2;
   vector<float>   *convTk_bcMatchingWithTracks_energy1;
   vector<float>   *convTk_bcMatchingWithTracks_energy2;
   vector<float>   *convTk_EoverP;
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
   vector<float>   *els_tk_pz;
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
   vector<float>   *els_n_inner_layer;
   vector<float>   *els_n_outer_layer;
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
   UInt_t          NhcalNoiseRBX;
   vector<float>   *hcalNoiseRBX_idnumber;
   vector<float>   *hcalNoiseRBX_allChargeTotal;
   vector<float>   *hcalNoiseRBX_allChargeHighest2TS;
   vector<float>   *hcalNoiseRBX_allChargeHighest3TS;
   vector<float>   *hcalNoiseRBX_totalZeros;
   vector<float>   *hcalNoiseRBX_maxZeros;
   vector<float>   *hcalNoiseRBX_recHitEnergy;
   vector<float>   *hcalNoiseRBX_minRecHitTime;
   vector<float>   *hcalNoiseRBX_maxRecHitTime;
   vector<float>   *hcalNoiseRBX_numRecHits;
   vector<float>   *hcalNoiseRBX_caloTowerHadE;
   vector<float>   *hcalNoiseRBX_caloTowerEmE;
   vector<float>   *hcalNoiseRBX_caloTowerTotalE;
   vector<float>   *hcalNoiseRBX_caloTowerEmFraction;
   UInt_t          NhcalNoiseSummary;
   vector<float>   *hcalNoiseSummary_passLooseNoiseFilter;
   vector<float>   *hcalNoiseSummary_passTightNoiseFilter;
   vector<float>   *hcalNoiseSummary_passHighLevelNoiseFilter;
   vector<float>   *hcalNoiseSummary_noiseFilterStatus;
   vector<float>   *hcalNoiseSummary_noiseType;
   vector<float>   *hcalNoiseSummary_eventEMEnergy;
   vector<float>   *hcalNoiseSummary_eventHadEnergy;
   vector<float>   *hcalNoiseSummary_eventTrackEnergy;
   vector<float>   *hcalNoiseSummary_eventEMFraction;
   vector<float>   *hcalNoiseSummary_eventChargeFraction;
   vector<float>   *hcalNoiseSummary_min10GeVHitTime;
   vector<float>   *hcalNoiseSummary_max10GeVHitTime;
   vector<float>   *hcalNoiseSummary_rms10GeVHitTime;
   vector<float>   *hcalNoiseSummary_min25GeVHitTime;
   vector<float>   *hcalNoiseSummary_max25GeVHitTime;
   vector<float>   *hcalNoiseSummary_rms25GeVHitTime;
   vector<float>   *hcalNoiseSummary_num10GeVHits;
   vector<float>   *hcalNoiseSummary_num25GeVHits;
   vector<float>   *hcalNoiseSummary_minE2TS;
   vector<float>   *hcalNoiseSummary_minE10TS;
   vector<float>   *hcalNoiseSummary_minE2Over10TS;
   vector<float>   *hcalNoiseSummary_maxZeros;
   vector<float>   *hcalNoiseSummary_maxHPDHits;
   vector<float>   *hcalNoiseSummary_maxRBXHits;
   vector<float>   *hcalNoiseSummary_minHPDEMF;
   vector<float>   *hcalNoiseSummary_minRBXEMF;
   vector<float>   *hcalNoiseSummary_numProblematicRBXs;
   UInt_t          NhybridBBC;
   vector<float>   *hybridBBC_energy;
   vector<float>   *hybridBBC_x;
   vector<float>   *hybridBBC_y;
   vector<float>   *hybridBBC_z;
   vector<float>   *hybridBBC_rho;
   vector<float>   *hybridBBC_phi;
   vector<float>   *hybridBBC_eta;
   vector<float>   *hybridBBC_theta;
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
   vector<float>   *jets_AK5_btag_secVertexHighPur;
   vector<float>   *jets_AK5_btag_secVertexHighEff;
   vector<float>   *jets_AK5_btag_secVertexCombined;
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
   vector<float>   *jets_AK5_etaetaMoment;
   vector<float>   *jets_AK5_etaphiMoment;
   vector<float>   *jets_AK5_phiphiMoment;
   vector<float>   *jets_AK5_n90Hits;
   vector<float>   *jets_AK5_fHPD;
   vector<float>   *jets_AK5_fRBX;
   vector<float>   *jets_AK5_hitsInN90;
   vector<float>   *jets_AK5_nECALTowers;
   vector<float>   *jets_AK5_nHCALTowers;
   vector<float>   *jets_AK5_fSubDetector1;
   vector<float>   *jets_AK5_fSubDetector2;
   vector<float>   *jets_AK5_fSubDetector3;
   vector<float>   *jets_AK5_fSubDetector4;
   vector<float>   *jets_AK5_area;
   vector<float>   *jets_AK5_mass;
   UInt_t          Njets_AK5JPT;
   vector<float>   *jets_AK5JPT_energy;
   vector<float>   *jets_AK5JPT_et;
   vector<float>   *jets_AK5JPT_eta;
   vector<float>   *jets_AK5JPT_phi;
   vector<float>   *jets_AK5JPT_pt;
   vector<float>   *jets_AK5JPT_px;
   vector<float>   *jets_AK5JPT_py;
   vector<float>   *jets_AK5JPT_pz;
   vector<float>   *jets_AK5JPT_status;
   vector<float>   *jets_AK5JPT_theta;
   vector<float>   *jets_AK5JPT_parton_Id;
   vector<float>   *jets_AK5JPT_parton_motherId;
   vector<float>   *jets_AK5JPT_parton_pt;
   vector<float>   *jets_AK5JPT_parton_phi;
   vector<float>   *jets_AK5JPT_parton_eta;
   vector<float>   *jets_AK5JPT_parton_Energy;
   vector<float>   *jets_AK5JPT_parton_mass;
   vector<float>   *jets_AK5JPT_gen_et;
   vector<float>   *jets_AK5JPT_gen_pt;
   vector<float>   *jets_AK5JPT_gen_eta;
   vector<float>   *jets_AK5JPT_gen_phi;
   vector<float>   *jets_AK5JPT_gen_mass;
   vector<float>   *jets_AK5JPT_gen_Energy;
   vector<float>   *jets_AK5JPT_gen_Id;
   vector<float>   *jets_AK5JPT_gen_motherID;
   vector<float>   *jets_AK5JPT_gen_threeCharge;
   vector<float>   *jets_AK5JPT_partonFlavour;
   vector<float>   *jets_AK5JPT_btag_TC_highPur;
   vector<float>   *jets_AK5JPT_btag_TC_highEff;
   vector<float>   *jets_AK5JPT_btag_jetProb;
   vector<float>   *jets_AK5JPT_btag_jetBProb;
   vector<float>   *jets_AK5JPT_btag_softEle;
   vector<float>   *jets_AK5JPT_btag_softMuon;
   vector<float>   *jets_AK5JPT_btag_secVertexHighPur;
   vector<float>   *jets_AK5JPT_btag_secVertexHighEff;
   vector<float>   *jets_AK5JPT_btag_secVertexCombined;
   vector<float>   *jets_AK5JPT_chgEmE;
   vector<float>   *jets_AK5JPT_chgHadE;
   vector<float>   *jets_AK5JPT_chgMuE;
   vector<float>   *jets_AK5JPT_chg_Mult;
   vector<float>   *jets_AK5JPT_neutralEmE;
   vector<float>   *jets_AK5JPT_neutralHadE;
   vector<float>   *jets_AK5JPT_neutral_Mult;
   vector<float>   *jets_AK5JPT_mu_Mult;
   vector<float>   *jets_AK5JPT_emf;
   vector<float>   *jets_AK5JPT_ehf;
   vector<float>   *jets_AK5JPT_n60;
   vector<float>   *jets_AK5JPT_n90;
   vector<float>   *jets_AK5JPT_etaetaMoment;
   vector<float>   *jets_AK5JPT_etaphiMoment;
   vector<float>   *jets_AK5JPT_phiphiMoment;
   vector<float>   *jets_AK5JPT_n90Hits;
   vector<float>   *jets_AK5JPT_fHPD;
   vector<float>   *jets_AK5JPT_fRBX;
   vector<float>   *jets_AK5JPT_hitsInN90;
   vector<float>   *jets_AK5JPT_nECALTowers;
   vector<float>   *jets_AK5JPT_nHCALTowers;
   vector<float>   *jets_AK5JPT_fSubDetector1;
   vector<float>   *jets_AK5JPT_fSubDetector2;
   vector<float>   *jets_AK5JPT_fSubDetector3;
   vector<float>   *jets_AK5JPT_fSubDetector4;
   vector<float>   *jets_AK5JPT_area;
   vector<float>   *jets_AK5JPT_mass;
   UInt_t          Njets_AK5PF;
   vector<float>   *jets_AK5PF_energy;
   vector<float>   *jets_AK5PF_et;
   vector<float>   *jets_AK5PF_eta;
   vector<float>   *jets_AK5PF_phi;
   vector<float>   *jets_AK5PF_pt;
   vector<float>   *jets_AK5PF_px;
   vector<float>   *jets_AK5PF_py;
   vector<float>   *jets_AK5PF_pz;
   vector<float>   *jets_AK5PF_status;
   vector<float>   *jets_AK5PF_theta;
   vector<float>   *jets_AK5PF_parton_Id;
   vector<float>   *jets_AK5PF_parton_motherId;
   vector<float>   *jets_AK5PF_parton_pt;
   vector<float>   *jets_AK5PF_parton_phi;
   vector<float>   *jets_AK5PF_parton_eta;
   vector<float>   *jets_AK5PF_parton_Energy;
   vector<float>   *jets_AK5PF_parton_mass;
   vector<float>   *jets_AK5PF_gen_et;
   vector<float>   *jets_AK5PF_gen_pt;
   vector<float>   *jets_AK5PF_gen_eta;
   vector<float>   *jets_AK5PF_gen_phi;
   vector<float>   *jets_AK5PF_gen_mass;
   vector<float>   *jets_AK5PF_gen_Energy;
   vector<float>   *jets_AK5PF_gen_Id;
   vector<float>   *jets_AK5PF_gen_motherID;
   vector<float>   *jets_AK5PF_gen_threeCharge;
   vector<float>   *jets_AK5PF_partonFlavour;
   vector<float>   *jets_AK5PF_btag_TC_highPur;
   vector<float>   *jets_AK5PF_btag_TC_highEff;
   vector<float>   *jets_AK5PF_btag_jetProb;
   vector<float>   *jets_AK5PF_btag_jetBProb;
   vector<float>   *jets_AK5PF_btag_softEle;
   vector<float>   *jets_AK5PF_btag_softMuon;
   vector<float>   *jets_AK5PF_btag_secVertexHighPur;
   vector<float>   *jets_AK5PF_btag_secVertexHighEff;
   vector<float>   *jets_AK5PF_btag_secVertexCombined;
   vector<float>   *jets_AK5PF_chgEmE;
   vector<float>   *jets_AK5PF_chgHadE;
   vector<float>   *jets_AK5PF_chgMuE;
   vector<float>   *jets_AK5PF_chg_Mult;
   vector<float>   *jets_AK5PF_neutralEmE;
   vector<float>   *jets_AK5PF_neutralHadE;
   vector<float>   *jets_AK5PF_neutral_Mult;
   vector<float>   *jets_AK5PF_mu_Mult;
   vector<float>   *jets_AK5PF_emf;
   vector<float>   *jets_AK5PF_ehf;
   vector<float>   *jets_AK5PF_n60;
   vector<float>   *jets_AK5PF_n90;
   vector<float>   *jets_AK5PF_etaetaMoment;
   vector<float>   *jets_AK5PF_etaphiMoment;
   vector<float>   *jets_AK5PF_phiphiMoment;
   vector<float>   *jets_AK5PF_n90Hits;
   vector<float>   *jets_AK5PF_fHPD;
   vector<float>   *jets_AK5PF_fRBX;
   vector<float>   *jets_AK5PF_hitsInN90;
   vector<float>   *jets_AK5PF_nECALTowers;
   vector<float>   *jets_AK5PF_nHCALTowers;
   vector<float>   *jets_AK5PF_fSubDetector1;
   vector<float>   *jets_AK5PF_fSubDetector2;
   vector<float>   *jets_AK5PF_fSubDetector3;
   vector<float>   *jets_AK5PF_fSubDetector4;
   vector<float>   *jets_AK5PF_area;
   vector<float>   *jets_AK5PF_mass;
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
   vector<float>   *mus_isTrackerMuon;
   vector<float>   *mus_isStandAloneMuon;
   vector<float>   *mus_isCaloMuon;
   vector<float>   *mus_isGlobalMuon;
   vector<float>   *mus_isElectron;
   vector<float>   *mus_isConvertedPhoton;
   vector<float>   *mus_isPhoton;
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
   UInt_t          Npfmets;
   vector<float>   *pfmets_et;
   vector<float>   *pfmets_phi;
   vector<float>   *pfmets_ex;
   vector<float>   *pfmets_ey;
   vector<float>   *pfmets_gen_et;
   vector<float>   *pfmets_gen_phi;
   vector<float>   *pfmets_sign;
   vector<float>   *pfmets_sumEt;
   vector<float>   *pfmets_unCPhi;
   vector<float>   *pfmets_unCPt;
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
   vector<float>   *photons_maxEnergyXtal;
   vector<float>   *photons_e1x5;
   vector<float>   *photons_e2x5;
   vector<float>   *photons_e3x3;
   vector<float>   *photons_e5x5;
   vector<float>   *photons_sigmaEtaEta;
   vector<float>   *photons_sigmaIetaIeta;
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
   vector<float>   *pv_chi2;
   vector<float>   *pv_ndof;
   vector<float>   *pv_tracksSize;
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
   vector<float>   *taus_taNC;
   vector<float>   *taus_byIsoUsingLeadingPi;
   vector<float>   *taus_tkIsoUsingLeadingPi;
   vector<float>   *taus_ecalIsoUsingLeadingPi;
   vector<float>   *taus_signalPFChargedHadrCandsSize;
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
   vector<float>   *tracks_trkExptHitsInner;
   vector<float>   *tracks_trkExptHitsOuter;
   vector<float>   *tracks_trks_nlayerslost;
   vector<float>   *tracks_trks_nlayers;
   vector<float>   *tracks_trksvalidpixelhits;
   vector<float>   *tracks_trkslostpixelhits;
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
   vector<float>   *tracks_outerPx;
   vector<float>   *tracks_outerPy;
   vector<float>   *tracks_outerPz;
   vector<float>   *tracks_algo;
   vector<float>   *tracks_highPurity;
   UInt_t          run;
   UInt_t          event;
   UInt_t          lumiblock;
   UInt_t          experimentType;
   UInt_t          bunchCrossing;
   UInt_t          orbitNumber;

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
   TBranch        *b_NconvTk;   //!
   TBranch        *b_convTk_tracksPin_x1;   //!
   TBranch        *b_convTk_tracksPin_x2;   //!
   TBranch        *b_convTk_tracksPin_y1;   //!
   TBranch        *b_convTk_tracksPin_y2;   //!
   TBranch        *b_convTk_tracksPin_z1;   //!
   TBranch        *b_convTk_tracksPin_z2;   //!
   TBranch        *b_convTk_tracksPout_x1;   //!
   TBranch        *b_convTk_tracksPout_x2;   //!
   TBranch        *b_convTk_tracksPout_y1;   //!
   TBranch        *b_convTk_tracksPout_y2;   //!
   TBranch        *b_convTk_tracksPout_z1;   //!
   TBranch        *b_convTk_tracksPout_z2;   //!
   TBranch        *b_convTk_ecalImpactPosition_x1;   //!
   TBranch        *b_convTk_ecalImpactPosition_x2;   //!
   TBranch        *b_convTk_ecalImpactPosition_y1;   //!
   TBranch        *b_convTk_ecalImpactPosition_y2;   //!
   TBranch        *b_convTk_ecalImpactPosition_z1;   //!
   TBranch        *b_convTk_ecalImpactPosition_z2;   //!
   TBranch        *b_convTk_tracksInnerPosition_x1;   //!
   TBranch        *b_convTk_tracksInnerPosition_x2;   //!
   TBranch        *b_convTk_tracksInnerPosition_y1;   //!
   TBranch        *b_convTk_tracksInnerPosition_y2;   //!
   TBranch        *b_convTk_tracksInnerPosition_z1;   //!
   TBranch        *b_convTk_tracksInnerPosition_z2;   //!
   TBranch        *b_convTk_tracks_px1;   //!
   TBranch        *b_convTk_tracks_py1;   //!
   TBranch        *b_convTk_tracks_pz1;   //!
   TBranch        *b_convTk_tracks_phi1;   //!
   TBranch        *b_convTk_tracks_d01;   //!
   TBranch        *b_convTk_tracks_dz1;   //!
   TBranch        *b_convTk_tracks_charge1;   //!
   TBranch        *b_convTk_tracks_px2;   //!
   TBranch        *b_convTk_tracks_py2;   //!
   TBranch        *b_convTk_tracks_pz2;   //!
   TBranch        *b_convTk_tracks_phi2;   //!
   TBranch        *b_convTk_tracks_d02;   //!
   TBranch        *b_convTk_tracks_dz2;   //!
   TBranch        *b_convTk_tracks_charge2;   //!
   TBranch        *b_convTk_bcMatchingWithTracks_x1;   //!
   TBranch        *b_convTk_bcMatchingWithTracks_y1;   //!
   TBranch        *b_convTk_bcMatchingWithTracks_z1;   //!
   TBranch        *b_convTk_bcMatchingWithTracks_x2;   //!
   TBranch        *b_convTk_bcMatchingWithTracks_y2;   //!
   TBranch        *b_convTk_bcMatchingWithTracks_z2;   //!
   TBranch        *b_convTk_bcMatchingWithTracks_energy1;   //!
   TBranch        *b_convTk_bcMatchingWithTracks_energy2;   //!
   TBranch        *b_convTk_EoverP;   //!
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
   TBranch        *b_els_tk_pz;   //!
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
   TBranch        *b_els_n_inner_layer;   //!
   TBranch        *b_els_n_outer_layer;   //!
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
   TBranch        *b_NhcalNoiseRBX;   //!
   TBranch        *b_hcalNoiseRBX_idnumber;   //!
   TBranch        *b_hcalNoiseRBX_allChargeTotal;   //!
   TBranch        *b_hcalNoiseRBX_allChargeHighest2TS;   //!
   TBranch        *b_hcalNoiseRBX_allChargeHighest3TS;   //!
   TBranch        *b_hcalNoiseRBX_totalZeros;   //!
   TBranch        *b_hcalNoiseRBX_maxZeros;   //!
   TBranch        *b_hcalNoiseRBX_recHitEnergy;   //!
   TBranch        *b_hcalNoiseRBX_minRecHitTime;   //!
   TBranch        *b_hcalNoiseRBX_maxRecHitTime;   //!
   TBranch        *b_hcalNoiseRBX_numRecHits;   //!
   TBranch        *b_hcalNoiseRBX_caloTowerHadE;   //!
   TBranch        *b_hcalNoiseRBX_caloTowerEmE;   //!
   TBranch        *b_hcalNoiseRBX_caloTowerTotalE;   //!
   TBranch        *b_hcalNoiseRBX_caloTowerEmFraction;   //!
   TBranch        *b_NhcalNoiseSummary;   //!
   TBranch        *b_hcalNoiseSummary_passLooseNoiseFilter;   //!
   TBranch        *b_hcalNoiseSummary_passTightNoiseFilter;   //!
   TBranch        *b_hcalNoiseSummary_passHighLevelNoiseFilter;   //!
   TBranch        *b_hcalNoiseSummary_noiseFilterStatus;   //!
   TBranch        *b_hcalNoiseSummary_noiseType;   //!
   TBranch        *b_hcalNoiseSummary_eventEMEnergy;   //!
   TBranch        *b_hcalNoiseSummary_eventHadEnergy;   //!
   TBranch        *b_hcalNoiseSummary_eventTrackEnergy;   //!
   TBranch        *b_hcalNoiseSummary_eventEMFraction;   //!
   TBranch        *b_hcalNoiseSummary_eventChargeFraction;   //!
   TBranch        *b_hcalNoiseSummary_min10GeVHitTime;   //!
   TBranch        *b_hcalNoiseSummary_max10GeVHitTime;   //!
   TBranch        *b_hcalNoiseSummary_rms10GeVHitTime;   //!
   TBranch        *b_hcalNoiseSummary_min25GeVHitTime;   //!
   TBranch        *b_hcalNoiseSummary_max25GeVHitTime;   //!
   TBranch        *b_hcalNoiseSummary_rms25GeVHitTime;   //!
   TBranch        *b_hcalNoiseSummary_num10GeVHits;   //!
   TBranch        *b_hcalNoiseSummary_num25GeVHits;   //!
   TBranch        *b_hcalNoiseSummary_minE2TS;   //!
   TBranch        *b_hcalNoiseSummary_minE10TS;   //!
   TBranch        *b_hcalNoiseSummary_minE2Over10TS;   //!
   TBranch        *b_hcalNoiseSummary_maxZeros;   //!
   TBranch        *b_hcalNoiseSummary_maxHPDHits;   //!
   TBranch        *b_hcalNoiseSummary_maxRBXHits;   //!
   TBranch        *b_hcalNoiseSummary_minHPDEMF;   //!
   TBranch        *b_hcalNoiseSummary_minRBXEMF;   //!
   TBranch        *b_hcalNoiseSummary_numProblematicRBXs;   //!
   TBranch        *b_NhybridBBC;   //!
   TBranch        *b_hybridBBC_energy;   //!
   TBranch        *b_hybridBBC_x;   //!
   TBranch        *b_hybridBBC_y;   //!
   TBranch        *b_hybridBBC_z;   //!
   TBranch        *b_hybridBBC_rho;   //!
   TBranch        *b_hybridBBC_phi;   //!
   TBranch        *b_hybridBBC_eta;   //!
   TBranch        *b_hybridBBC_theta;   //!
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
   TBranch        *b_jets_AK5_btag_secVertexHighPur;   //!
   TBranch        *b_jets_AK5_btag_secVertexHighEff;   //!
   TBranch        *b_jets_AK5_btag_secVertexCombined;   //!
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
   TBranch        *b_jets_AK5_etaetaMoment;   //!
   TBranch        *b_jets_AK5_etaphiMoment;   //!
   TBranch        *b_jets_AK5_phiphiMoment;   //!
   TBranch        *b_jets_AK5_n90Hits;   //!
   TBranch        *b_jets_AK5_fHPD;   //!
   TBranch        *b_jets_AK5_fRBX;   //!
   TBranch        *b_jets_AK5_hitsInN90;   //!
   TBranch        *b_jets_AK5_nECALTowers;   //!
   TBranch        *b_jets_AK5_nHCALTowers;   //!
   TBranch        *b_jets_AK5_fSubDetector1;   //!
   TBranch        *b_jets_AK5_fSubDetector2;   //!
   TBranch        *b_jets_AK5_fSubDetector3;   //!
   TBranch        *b_jets_AK5_fSubDetector4;   //!
   TBranch        *b_jets_AK5_area;   //!
   TBranch        *b_jets_AK5_mass;   //!
   TBranch        *b_Njets_AK5JPT;   //!
   TBranch        *b_jets_AK5JPT_energy;   //!
   TBranch        *b_jets_AK5JPT_et;   //!
   TBranch        *b_jets_AK5JPT_eta;   //!
   TBranch        *b_jets_AK5JPT_phi;   //!
   TBranch        *b_jets_AK5JPT_pt;   //!
   TBranch        *b_jets_AK5JPT_px;   //!
   TBranch        *b_jets_AK5JPT_py;   //!
   TBranch        *b_jets_AK5JPT_pz;   //!
   TBranch        *b_jets_AK5JPT_status;   //!
   TBranch        *b_jets_AK5JPT_theta;   //!
   TBranch        *b_jets_AK5JPT_parton_Id;   //!
   TBranch        *b_jets_AK5JPT_parton_motherId;   //!
   TBranch        *b_jets_AK5JPT_parton_pt;   //!
   TBranch        *b_jets_AK5JPT_parton_phi;   //!
   TBranch        *b_jets_AK5JPT_parton_eta;   //!
   TBranch        *b_jets_AK5JPT_parton_Energy;   //!
   TBranch        *b_jets_AK5JPT_parton_mass;   //!
   TBranch        *b_jets_AK5JPT_gen_et;   //!
   TBranch        *b_jets_AK5JPT_gen_pt;   //!
   TBranch        *b_jets_AK5JPT_gen_eta;   //!
   TBranch        *b_jets_AK5JPT_gen_phi;   //!
   TBranch        *b_jets_AK5JPT_gen_mass;   //!
   TBranch        *b_jets_AK5JPT_gen_Energy;   //!
   TBranch        *b_jets_AK5JPT_gen_Id;   //!
   TBranch        *b_jets_AK5JPT_gen_motherID;   //!
   TBranch        *b_jets_AK5JPT_gen_threeCharge;   //!
   TBranch        *b_jets_AK5JPT_partonFlavour;   //!
   TBranch        *b_jets_AK5JPT_btag_TC_highPur;   //!
   TBranch        *b_jets_AK5JPT_btag_TC_highEff;   //!
   TBranch        *b_jets_AK5JPT_btag_jetProb;   //!
   TBranch        *b_jets_AK5JPT_btag_jetBProb;   //!
   TBranch        *b_jets_AK5JPT_btag_softEle;   //!
   TBranch        *b_jets_AK5JPT_btag_softMuon;   //!
   TBranch        *b_jets_AK5JPT_btag_secVertexHighPur;   //!
   TBranch        *b_jets_AK5JPT_btag_secVertexHighEff;   //!
   TBranch        *b_jets_AK5JPT_btag_secVertexCombined;   //!
   TBranch        *b_jets_AK5JPT_chgEmE;   //!
   TBranch        *b_jets_AK5JPT_chgHadE;   //!
   TBranch        *b_jets_AK5JPT_chgMuE;   //!
   TBranch        *b_jets_AK5JPT_chg_Mult;   //!
   TBranch        *b_jets_AK5JPT_neutralEmE;   //!
   TBranch        *b_jets_AK5JPT_neutralHadE;   //!
   TBranch        *b_jets_AK5JPT_neutral_Mult;   //!
   TBranch        *b_jets_AK5JPT_mu_Mult;   //!
   TBranch        *b_jets_AK5JPT_emf;   //!
   TBranch        *b_jets_AK5JPT_ehf;   //!
   TBranch        *b_jets_AK5JPT_n60;   //!
   TBranch        *b_jets_AK5JPT_n90;   //!
   TBranch        *b_jets_AK5JPT_etaetaMoment;   //!
   TBranch        *b_jets_AK5JPT_etaphiMoment;   //!
   TBranch        *b_jets_AK5JPT_phiphiMoment;   //!
   TBranch        *b_jets_AK5JPT_n90Hits;   //!
   TBranch        *b_jets_AK5JPT_fHPD;   //!
   TBranch        *b_jets_AK5JPT_fRBX;   //!
   TBranch        *b_jets_AK5JPT_hitsInN90;   //!
   TBranch        *b_jets_AK5JPT_nECALTowers;   //!
   TBranch        *b_jets_AK5JPT_nHCALTowers;   //!
   TBranch        *b_jets_AK5JPT_fSubDetector1;   //!
   TBranch        *b_jets_AK5JPT_fSubDetector2;   //!
   TBranch        *b_jets_AK5JPT_fSubDetector3;   //!
   TBranch        *b_jets_AK5JPT_fSubDetector4;   //!
   TBranch        *b_jets_AK5JPT_area;   //!
   TBranch        *b_jets_AK5JPT_mass;   //!
   TBranch        *b_Njets_AK5PF;   //!
   TBranch        *b_jets_AK5PF_energy;   //!
   TBranch        *b_jets_AK5PF_et;   //!
   TBranch        *b_jets_AK5PF_eta;   //!
   TBranch        *b_jets_AK5PF_phi;   //!
   TBranch        *b_jets_AK5PF_pt;   //!
   TBranch        *b_jets_AK5PF_px;   //!
   TBranch        *b_jets_AK5PF_py;   //!
   TBranch        *b_jets_AK5PF_pz;   //!
   TBranch        *b_jets_AK5PF_status;   //!
   TBranch        *b_jets_AK5PF_theta;   //!
   TBranch        *b_jets_AK5PF_parton_Id;   //!
   TBranch        *b_jets_AK5PF_parton_motherId;   //!
   TBranch        *b_jets_AK5PF_parton_pt;   //!
   TBranch        *b_jets_AK5PF_parton_phi;   //!
   TBranch        *b_jets_AK5PF_parton_eta;   //!
   TBranch        *b_jets_AK5PF_parton_Energy;   //!
   TBranch        *b_jets_AK5PF_parton_mass;   //!
   TBranch        *b_jets_AK5PF_gen_et;   //!
   TBranch        *b_jets_AK5PF_gen_pt;   //!
   TBranch        *b_jets_AK5PF_gen_eta;   //!
   TBranch        *b_jets_AK5PF_gen_phi;   //!
   TBranch        *b_jets_AK5PF_gen_mass;   //!
   TBranch        *b_jets_AK5PF_gen_Energy;   //!
   TBranch        *b_jets_AK5PF_gen_Id;   //!
   TBranch        *b_jets_AK5PF_gen_motherID;   //!
   TBranch        *b_jets_AK5PF_gen_threeCharge;   //!
   TBranch        *b_jets_AK5PF_partonFlavour;   //!
   TBranch        *b_jets_AK5PF_btag_TC_highPur;   //!
   TBranch        *b_jets_AK5PF_btag_TC_highEff;   //!
   TBranch        *b_jets_AK5PF_btag_jetProb;   //!
   TBranch        *b_jets_AK5PF_btag_jetBProb;   //!
   TBranch        *b_jets_AK5PF_btag_softEle;   //!
   TBranch        *b_jets_AK5PF_btag_softMuon;   //!
   TBranch        *b_jets_AK5PF_btag_secVertexHighPur;   //!
   TBranch        *b_jets_AK5PF_btag_secVertexHighEff;   //!
   TBranch        *b_jets_AK5PF_btag_secVertexCombined;   //!
   TBranch        *b_jets_AK5PF_chgEmE;   //!
   TBranch        *b_jets_AK5PF_chgHadE;   //!
   TBranch        *b_jets_AK5PF_chgMuE;   //!
   TBranch        *b_jets_AK5PF_chg_Mult;   //!
   TBranch        *b_jets_AK5PF_neutralEmE;   //!
   TBranch        *b_jets_AK5PF_neutralHadE;   //!
   TBranch        *b_jets_AK5PF_neutral_Mult;   //!
   TBranch        *b_jets_AK5PF_mu_Mult;   //!
   TBranch        *b_jets_AK5PF_emf;   //!
   TBranch        *b_jets_AK5PF_ehf;   //!
   TBranch        *b_jets_AK5PF_n60;   //!
   TBranch        *b_jets_AK5PF_n90;   //!
   TBranch        *b_jets_AK5PF_etaetaMoment;   //!
   TBranch        *b_jets_AK5PF_etaphiMoment;   //!
   TBranch        *b_jets_AK5PF_phiphiMoment;   //!
   TBranch        *b_jets_AK5PF_n90Hits;   //!
   TBranch        *b_jets_AK5PF_fHPD;   //!
   TBranch        *b_jets_AK5PF_fRBX;   //!
   TBranch        *b_jets_AK5PF_hitsInN90;   //!
   TBranch        *b_jets_AK5PF_nECALTowers;   //!
   TBranch        *b_jets_AK5PF_nHCALTowers;   //!
   TBranch        *b_jets_AK5PF_fSubDetector1;   //!
   TBranch        *b_jets_AK5PF_fSubDetector2;   //!
   TBranch        *b_jets_AK5PF_fSubDetector3;   //!
   TBranch        *b_jets_AK5PF_fSubDetector4;   //!
   TBranch        *b_jets_AK5PF_area;   //!
   TBranch        *b_jets_AK5PF_mass;   //!
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
   TBranch        *b_mus_isTrackerMuon;   //!
   TBranch        *b_mus_isStandAloneMuon;   //!
   TBranch        *b_mus_isCaloMuon;   //!
   TBranch        *b_mus_isGlobalMuon;   //!
   TBranch        *b_mus_isElectron;   //!
   TBranch        *b_mus_isConvertedPhoton;   //!
   TBranch        *b_mus_isPhoton;   //!
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
   TBranch        *b_Npfmets;   //!
   TBranch        *b_pfmets_et;   //!
   TBranch        *b_pfmets_phi;   //!
   TBranch        *b_pfmets_ex;   //!
   TBranch        *b_pfmets_ey;   //!
   TBranch        *b_pfmets_gen_et;   //!
   TBranch        *b_pfmets_gen_phi;   //!
   TBranch        *b_pfmets_sign;   //!
   TBranch        *b_pfmets_sumEt;   //!
   TBranch        *b_pfmets_unCPhi;   //!
   TBranch        *b_pfmets_unCPt;   //!
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
   TBranch        *b_photons_maxEnergyXtal;   //!
   TBranch        *b_photons_e1x5;   //!
   TBranch        *b_photons_e2x5;   //!
   TBranch        *b_photons_e3x3;   //!
   TBranch        *b_photons_e5x5;   //!
   TBranch        *b_photons_sigmaEtaEta;   //!
   TBranch        *b_photons_sigmaIetaIeta;   //!
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
   TBranch        *b_pv_chi2;   //!
   TBranch        *b_pv_ndof;   //!
   TBranch        *b_pv_tracksSize;   //!
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
   TBranch        *b_taus_taNC;   //!
   TBranch        *b_taus_byIsoUsingLeadingPi;   //!
   TBranch        *b_taus_tkIsoUsingLeadingPi;   //!
   TBranch        *b_taus_ecalIsoUsingLeadingPi;   //!
   TBranch        *b_taus_signalPFChargedHadrCandsSize;   //!
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
   TBranch        *b_tracks_trkExptHitsInner;   //!
   TBranch        *b_tracks_trkExptHitsOuter;   //!
   TBranch        *b_tracks_trks_nlayerslost;   //!
   TBranch        *b_tracks_trks_nlayers;   //!
   TBranch        *b_tracks_trksvalidpixelhits;   //!
   TBranch        *b_tracks_trkslostpixelhits;   //!
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
   TBranch        *b_tracks_outerPx;   //!
   TBranch        *b_tracks_outerPy;   //!
   TBranch        *b_tracks_outerPz;   //!
   TBranch        *b_tracks_algo;   //!
   TBranch        *b_tracks_highPurity;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumiblock;   //!
   TBranch        *b_experimentType;   //!
   TBranch        *b_bunchCrossing;   //!
   TBranch        *b_orbitNumber;   //!

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
   Double_t        L1Bit_0;
   Double_t        L1Bit_1;
   Double_t        L1Bit_10;
   Double_t        L1Bit_100;
   Double_t        L1Bit_101;
   Double_t        L1Bit_102;
   Double_t        L1Bit_103;
   Double_t        L1Bit_104;
   Double_t        L1Bit_105;
   Double_t        L1Bit_106;
   Double_t        L1Bit_107;
   Double_t        L1Bit_108;
   Double_t        L1Bit_109;
   Double_t        L1Bit_11;
   Double_t        L1Bit_110;
   Double_t        L1Bit_111;
   Double_t        L1Bit_112;
   Double_t        L1Bit_113;
   Double_t        L1Bit_114;
   Double_t        L1Bit_115;
   Double_t        L1Bit_116;
   Double_t        L1Bit_117;
   Double_t        L1Bit_118;
   Double_t        L1Bit_119;
   Double_t        L1Bit_12;
   Double_t        L1Bit_120;
   Double_t        L1Bit_121;
   Double_t        L1Bit_122;
   Double_t        L1Bit_123;
   Double_t        L1Bit_124;
   Double_t        L1Bit_125;
   Double_t        L1Bit_126;
   Double_t        L1Bit_127;
   Double_t        L1Bit_13;
   Double_t        L1Bit_14;
   Double_t        L1Bit_15;
   Double_t        L1Bit_16;
   Double_t        L1Bit_17;
   Double_t        L1Bit_18;
   Double_t        L1Bit_19;
   Double_t        L1Bit_2;
   Double_t        L1Bit_20;
   Double_t        L1Bit_21;
   Double_t        L1Bit_22;
   Double_t        L1Bit_23;
   Double_t        L1Bit_24;
   Double_t        L1Bit_25;
   Double_t        L1Bit_26;
   Double_t        L1Bit_27;
   Double_t        L1Bit_28;
   Double_t        L1Bit_29;
   Double_t        L1Bit_3;
   Double_t        L1Bit_30;
   Double_t        L1Bit_31;
   Double_t        L1Bit_32;
   Double_t        L1Bit_33;
   Double_t        L1Bit_34;
   Double_t        L1Bit_35;
   Double_t        L1Bit_36;
   Double_t        L1Bit_37;
   Double_t        L1Bit_38;
   Double_t        L1Bit_39;
   Double_t        L1Bit_4;
   Double_t        L1Bit_40;
   Double_t        L1Bit_41;
   Double_t        L1Bit_42;
   Double_t        L1Bit_43;
   Double_t        L1Bit_44;
   Double_t        L1Bit_45;
   Double_t        L1Bit_46;
   Double_t        L1Bit_47;
   Double_t        L1Bit_48;
   Double_t        L1Bit_49;
   Double_t        L1Bit_5;
   Double_t        L1Bit_50;
   Double_t        L1Bit_51;
   Double_t        L1Bit_52;
   Double_t        L1Bit_53;
   Double_t        L1Bit_54;
   Double_t        L1Bit_55;
   Double_t        L1Bit_56;
   Double_t        L1Bit_57;
   Double_t        L1Bit_58;
   Double_t        L1Bit_59;
   Double_t        L1Bit_6;
   Double_t        L1Bit_60;
   Double_t        L1Bit_61;
   Double_t        L1Bit_62;
   Double_t        L1Bit_63;
   Double_t        L1Bit_64;
   Double_t        L1Bit_65;
   Double_t        L1Bit_66;
   Double_t        L1Bit_67;
   Double_t        L1Bit_68;
   Double_t        L1Bit_69;
   Double_t        L1Bit_7;
   Double_t        L1Bit_70;
   Double_t        L1Bit_71;
   Double_t        L1Bit_72;
   Double_t        L1Bit_73;
   Double_t        L1Bit_74;
   Double_t        L1Bit_75;
   Double_t        L1Bit_76;
   Double_t        L1Bit_77;
   Double_t        L1Bit_78;
   Double_t        L1Bit_79;
   Double_t        L1Bit_8;
   Double_t        L1Bit_80;
   Double_t        L1Bit_81;
   Double_t        L1Bit_82;
   Double_t        L1Bit_83;
   Double_t        L1Bit_84;
   Double_t        L1Bit_85;
   Double_t        L1Bit_86;
   Double_t        L1Bit_87;
   Double_t        L1Bit_88;
   Double_t        L1Bit_89;
   Double_t        L1Bit_9;
   Double_t        L1Bit_90;
   Double_t        L1Bit_91;
   Double_t        L1Bit_92;
   Double_t        L1Bit_93;
   Double_t        L1Bit_94;
   Double_t        L1Bit_95;
   Double_t        L1Bit_96;
   Double_t        L1Bit_97;
   Double_t        L1Bit_98;
   Double_t        L1Bit_99;

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
   TBranch        *b_L1Bit_0;   //!
   TBranch        *b_L1Bit_1;   //!
   TBranch        *b_L1Bit_10;   //!
   TBranch        *b_L1Bit_100;   //!
   TBranch        *b_L1Bit_101;   //!
   TBranch        *b_L1Bit_102;   //!
   TBranch        *b_L1Bit_103;   //!
   TBranch        *b_L1Bit_104;   //!
   TBranch        *b_L1Bit_105;   //!
   TBranch        *b_L1Bit_106;   //!
   TBranch        *b_L1Bit_107;   //!
   TBranch        *b_L1Bit_108;   //!
   TBranch        *b_L1Bit_109;   //!
   TBranch        *b_L1Bit_11;   //!
   TBranch        *b_L1Bit_110;   //!
   TBranch        *b_L1Bit_111;   //!
   TBranch        *b_L1Bit_112;   //!
   TBranch        *b_L1Bit_113;   //!
   TBranch        *b_L1Bit_114;   //!
   TBranch        *b_L1Bit_115;   //!
   TBranch        *b_L1Bit_116;   //!
   TBranch        *b_L1Bit_117;   //!
   TBranch        *b_L1Bit_118;   //!
   TBranch        *b_L1Bit_119;   //!
   TBranch        *b_L1Bit_12;   //!
   TBranch        *b_L1Bit_120;   //!
   TBranch        *b_L1Bit_121;   //!
   TBranch        *b_L1Bit_122;   //!
   TBranch        *b_L1Bit_123;   //!
   TBranch        *b_L1Bit_124;   //!
   TBranch        *b_L1Bit_125;   //!
   TBranch        *b_L1Bit_126;   //!
   TBranch        *b_L1Bit_127;   //!
   TBranch        *b_L1Bit_13;   //!
   TBranch        *b_L1Bit_14;   //!
   TBranch        *b_L1Bit_15;   //!
   TBranch        *b_L1Bit_16;   //!
   TBranch        *b_L1Bit_17;   //!
   TBranch        *b_L1Bit_18;   //!
   TBranch        *b_L1Bit_19;   //!
   TBranch        *b_L1Bit_2;   //!
   TBranch        *b_L1Bit_20;   //!
   TBranch        *b_L1Bit_21;   //!
   TBranch        *b_L1Bit_22;   //!
   TBranch        *b_L1Bit_23;   //!
   TBranch        *b_L1Bit_24;   //!
   TBranch        *b_L1Bit_25;   //!
   TBranch        *b_L1Bit_26;   //!
   TBranch        *b_L1Bit_27;   //!
   TBranch        *b_L1Bit_28;   //!
   TBranch        *b_L1Bit_29;   //!
   TBranch        *b_L1Bit_3;   //!
   TBranch        *b_L1Bit_30;   //!
   TBranch        *b_L1Bit_31;   //!
   TBranch        *b_L1Bit_32;   //!
   TBranch        *b_L1Bit_33;   //!
   TBranch        *b_L1Bit_34;   //!
   TBranch        *b_L1Bit_35;   //!
   TBranch        *b_L1Bit_36;   //!
   TBranch        *b_L1Bit_37;   //!
   TBranch        *b_L1Bit_38;   //!
   TBranch        *b_L1Bit_39;   //!
   TBranch        *b_L1Bit_4;   //!
   TBranch        *b_L1Bit_40;   //!
   TBranch        *b_L1Bit_41;   //!
   TBranch        *b_L1Bit_42;   //!
   TBranch        *b_L1Bit_43;   //!
   TBranch        *b_L1Bit_44;   //!
   TBranch        *b_L1Bit_45;   //!
   TBranch        *b_L1Bit_46;   //!
   TBranch        *b_L1Bit_47;   //!
   TBranch        *b_L1Bit_48;   //!
   TBranch        *b_L1Bit_49;   //!
   TBranch        *b_L1Bit_5;   //!
   TBranch        *b_L1Bit_50;   //!
   TBranch        *b_L1Bit_51;   //!
   TBranch        *b_L1Bit_52;   //!
   TBranch        *b_L1Bit_53;   //!
   TBranch        *b_L1Bit_54;   //!
   TBranch        *b_L1Bit_55;   //!
   TBranch        *b_L1Bit_56;   //!
   TBranch        *b_L1Bit_57;   //!
   TBranch        *b_L1Bit_58;   //!
   TBranch        *b_L1Bit_59;   //!
   TBranch        *b_L1Bit_6;   //!
   TBranch        *b_L1Bit_60;   //!
   TBranch        *b_L1Bit_61;   //!
   TBranch        *b_L1Bit_62;   //!
   TBranch        *b_L1Bit_63;   //!
   TBranch        *b_L1Bit_64;   //!
   TBranch        *b_L1Bit_65;   //!
   TBranch        *b_L1Bit_66;   //!
   TBranch        *b_L1Bit_67;   //!
   TBranch        *b_L1Bit_68;   //!
   TBranch        *b_L1Bit_69;   //!
   TBranch        *b_L1Bit_7;   //!
   TBranch        *b_L1Bit_70;   //!
   TBranch        *b_L1Bit_71;   //!
   TBranch        *b_L1Bit_72;   //!
   TBranch        *b_L1Bit_73;   //!
   TBranch        *b_L1Bit_74;   //!
   TBranch        *b_L1Bit_75;   //!
   TBranch        *b_L1Bit_76;   //!
   TBranch        *b_L1Bit_77;   //!
   TBranch        *b_L1Bit_78;   //!
   TBranch        *b_L1Bit_79;   //!
   TBranch        *b_L1Bit_8;   //!
   TBranch        *b_L1Bit_80;   //!
   TBranch        *b_L1Bit_81;   //!
   TBranch        *b_L1Bit_82;   //!
   TBranch        *b_L1Bit_83;   //!
   TBranch        *b_L1Bit_84;   //!
   TBranch        *b_L1Bit_85;   //!
   TBranch        *b_L1Bit_86;   //!
   TBranch        *b_L1Bit_87;   //!
   TBranch        *b_L1Bit_88;   //!
   TBranch        *b_L1Bit_89;   //!
   TBranch        *b_L1Bit_9;   //!
   TBranch        *b_L1Bit_90;   //!
   TBranch        *b_L1Bit_91;   //!
   TBranch        *b_L1Bit_92;   //!
   TBranch        *b_L1Bit_93;   //!
   TBranch        *b_L1Bit_94;   //!
   TBranch        *b_L1Bit_95;   //!
   TBranch        *b_L1Bit_96;   //!
   TBranch        *b_L1Bit_97;   //!
   TBranch        *b_L1Bit_98;   //!
   TBranch        *b_L1Bit_99;   //!


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
   convTk_tracksPin_x1 = 0;
   convTk_tracksPin_x2 = 0;
   convTk_tracksPin_y1 = 0;
   convTk_tracksPin_y2 = 0;
   convTk_tracksPin_z1 = 0;
   convTk_tracksPin_z2 = 0;
   convTk_tracksPout_x1 = 0;
   convTk_tracksPout_x2 = 0;
   convTk_tracksPout_y1 = 0;
   convTk_tracksPout_y2 = 0;
   convTk_tracksPout_z1 = 0;
   convTk_tracksPout_z2 = 0;
   convTk_ecalImpactPosition_x1 = 0;
   convTk_ecalImpactPosition_x2 = 0;
   convTk_ecalImpactPosition_y1 = 0;
   convTk_ecalImpactPosition_y2 = 0;
   convTk_ecalImpactPosition_z1 = 0;
   convTk_ecalImpactPosition_z2 = 0;
   convTk_tracksInnerPosition_x1 = 0;
   convTk_tracksInnerPosition_x2 = 0;
   convTk_tracksInnerPosition_y1 = 0;
   convTk_tracksInnerPosition_y2 = 0;
   convTk_tracksInnerPosition_z1 = 0;
   convTk_tracksInnerPosition_z2 = 0;
   convTk_tracks_px1 = 0;
   convTk_tracks_py1 = 0;
   convTk_tracks_pz1 = 0;
   convTk_tracks_phi1 = 0;
   convTk_tracks_d01 = 0;
   convTk_tracks_dz1 = 0;
   convTk_tracks_charge1 = 0;
   convTk_tracks_px2 = 0;
   convTk_tracks_py2 = 0;
   convTk_tracks_pz2 = 0;
   convTk_tracks_phi2 = 0;
   convTk_tracks_d02 = 0;
   convTk_tracks_dz2 = 0;
   convTk_tracks_charge2 = 0;
   convTk_bcMatchingWithTracks_x1 = 0;
   convTk_bcMatchingWithTracks_y1 = 0;
   convTk_bcMatchingWithTracks_z1 = 0;
   convTk_bcMatchingWithTracks_x2 = 0;
   convTk_bcMatchingWithTracks_y2 = 0;
   convTk_bcMatchingWithTracks_z2 = 0;
   convTk_bcMatchingWithTracks_energy1 = 0;
   convTk_bcMatchingWithTracks_energy2 = 0;
   convTk_EoverP = 0;
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
   els_tk_pz = 0;
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
   els_n_inner_layer = 0;
   els_n_outer_layer = 0;
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
   hcalNoiseRBX_idnumber = 0;
   hcalNoiseRBX_allChargeTotal = 0;
   hcalNoiseRBX_allChargeHighest2TS = 0;
   hcalNoiseRBX_allChargeHighest3TS = 0;
   hcalNoiseRBX_totalZeros = 0;
   hcalNoiseRBX_maxZeros = 0;
   hcalNoiseRBX_recHitEnergy = 0;
   hcalNoiseRBX_minRecHitTime = 0;
   hcalNoiseRBX_maxRecHitTime = 0;
   hcalNoiseRBX_numRecHits = 0;
   hcalNoiseRBX_caloTowerHadE = 0;
   hcalNoiseRBX_caloTowerEmE = 0;
   hcalNoiseRBX_caloTowerTotalE = 0;
   hcalNoiseRBX_caloTowerEmFraction = 0;
   hcalNoiseSummary_passLooseNoiseFilter = 0;
   hcalNoiseSummary_passTightNoiseFilter = 0;
   hcalNoiseSummary_passHighLevelNoiseFilter = 0;
   hcalNoiseSummary_noiseFilterStatus = 0;
   hcalNoiseSummary_noiseType = 0;
   hcalNoiseSummary_eventEMEnergy = 0;
   hcalNoiseSummary_eventHadEnergy = 0;
   hcalNoiseSummary_eventTrackEnergy = 0;
   hcalNoiseSummary_eventEMFraction = 0;
   hcalNoiseSummary_eventChargeFraction = 0;
   hcalNoiseSummary_min10GeVHitTime = 0;
   hcalNoiseSummary_max10GeVHitTime = 0;
   hcalNoiseSummary_rms10GeVHitTime = 0;
   hcalNoiseSummary_min25GeVHitTime = 0;
   hcalNoiseSummary_max25GeVHitTime = 0;
   hcalNoiseSummary_rms25GeVHitTime = 0;
   hcalNoiseSummary_num10GeVHits = 0;
   hcalNoiseSummary_num25GeVHits = 0;
   hcalNoiseSummary_minE2TS = 0;
   hcalNoiseSummary_minE10TS = 0;
   hcalNoiseSummary_minE2Over10TS = 0;
   hcalNoiseSummary_maxZeros = 0;
   hcalNoiseSummary_maxHPDHits = 0;
   hcalNoiseSummary_maxRBXHits = 0;
   hcalNoiseSummary_minHPDEMF = 0;
   hcalNoiseSummary_minRBXEMF = 0;
   hcalNoiseSummary_numProblematicRBXs = 0;
   hybridBBC_energy = 0;
   hybridBBC_x = 0;
   hybridBBC_y = 0;
   hybridBBC_z = 0;
   hybridBBC_rho = 0;
   hybridBBC_phi = 0;
   hybridBBC_eta = 0;
   hybridBBC_theta = 0;
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
   jets_AK5_btag_secVertexHighPur = 0;
   jets_AK5_btag_secVertexHighEff = 0;
   jets_AK5_btag_secVertexCombined = 0;
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
   jets_AK5_etaetaMoment = 0;
   jets_AK5_etaphiMoment = 0;
   jets_AK5_phiphiMoment = 0;
   jets_AK5_n90Hits = 0;
   jets_AK5_fHPD = 0;
   jets_AK5_fRBX = 0;
   jets_AK5_hitsInN90 = 0;
   jets_AK5_nECALTowers = 0;
   jets_AK5_nHCALTowers = 0;
   jets_AK5_fSubDetector1 = 0;
   jets_AK5_fSubDetector2 = 0;
   jets_AK5_fSubDetector3 = 0;
   jets_AK5_fSubDetector4 = 0;
   jets_AK5_area = 0;
   jets_AK5_mass = 0;
   jets_AK5JPT_energy = 0;
   jets_AK5JPT_et = 0;
   jets_AK5JPT_eta = 0;
   jets_AK5JPT_phi = 0;
   jets_AK5JPT_pt = 0;
   jets_AK5JPT_px = 0;
   jets_AK5JPT_py = 0;
   jets_AK5JPT_pz = 0;
   jets_AK5JPT_status = 0;
   jets_AK5JPT_theta = 0;
   jets_AK5JPT_parton_Id = 0;
   jets_AK5JPT_parton_motherId = 0;
   jets_AK5JPT_parton_pt = 0;
   jets_AK5JPT_parton_phi = 0;
   jets_AK5JPT_parton_eta = 0;
   jets_AK5JPT_parton_Energy = 0;
   jets_AK5JPT_parton_mass = 0;
   jets_AK5JPT_gen_et = 0;
   jets_AK5JPT_gen_pt = 0;
   jets_AK5JPT_gen_eta = 0;
   jets_AK5JPT_gen_phi = 0;
   jets_AK5JPT_gen_mass = 0;
   jets_AK5JPT_gen_Energy = 0;
   jets_AK5JPT_gen_Id = 0;
   jets_AK5JPT_gen_motherID = 0;
   jets_AK5JPT_gen_threeCharge = 0;
   jets_AK5JPT_partonFlavour = 0;
   jets_AK5JPT_btag_TC_highPur = 0;
   jets_AK5JPT_btag_TC_highEff = 0;
   jets_AK5JPT_btag_jetProb = 0;
   jets_AK5JPT_btag_jetBProb = 0;
   jets_AK5JPT_btag_softEle = 0;
   jets_AK5JPT_btag_softMuon = 0;
   jets_AK5JPT_btag_secVertexHighPur = 0;
   jets_AK5JPT_btag_secVertexHighEff = 0;
   jets_AK5JPT_btag_secVertexCombined = 0;
   jets_AK5JPT_chgEmE = 0;
   jets_AK5JPT_chgHadE = 0;
   jets_AK5JPT_chgMuE = 0;
   jets_AK5JPT_chg_Mult = 0;
   jets_AK5JPT_neutralEmE = 0;
   jets_AK5JPT_neutralHadE = 0;
   jets_AK5JPT_neutral_Mult = 0;
   jets_AK5JPT_mu_Mult = 0;
   jets_AK5JPT_emf = 0;
   jets_AK5JPT_ehf = 0;
   jets_AK5JPT_n60 = 0;
   jets_AK5JPT_n90 = 0;
   jets_AK5JPT_etaetaMoment = 0;
   jets_AK5JPT_etaphiMoment = 0;
   jets_AK5JPT_phiphiMoment = 0;
   jets_AK5JPT_n90Hits = 0;
   jets_AK5JPT_fHPD = 0;
   jets_AK5JPT_fRBX = 0;
   jets_AK5JPT_hitsInN90 = 0;
   jets_AK5JPT_nECALTowers = 0;
   jets_AK5JPT_nHCALTowers = 0;
   jets_AK5JPT_fSubDetector1 = 0;
   jets_AK5JPT_fSubDetector2 = 0;
   jets_AK5JPT_fSubDetector3 = 0;
   jets_AK5JPT_fSubDetector4 = 0;
   jets_AK5JPT_area = 0;
   jets_AK5JPT_mass = 0;
   jets_AK5PF_energy = 0;
   jets_AK5PF_et = 0;
   jets_AK5PF_eta = 0;
   jets_AK5PF_phi = 0;
   jets_AK5PF_pt = 0;
   jets_AK5PF_px = 0;
   jets_AK5PF_py = 0;
   jets_AK5PF_pz = 0;
   jets_AK5PF_status = 0;
   jets_AK5PF_theta = 0;
   jets_AK5PF_parton_Id = 0;
   jets_AK5PF_parton_motherId = 0;
   jets_AK5PF_parton_pt = 0;
   jets_AK5PF_parton_phi = 0;
   jets_AK5PF_parton_eta = 0;
   jets_AK5PF_parton_Energy = 0;
   jets_AK5PF_parton_mass = 0;
   jets_AK5PF_gen_et = 0;
   jets_AK5PF_gen_pt = 0;
   jets_AK5PF_gen_eta = 0;
   jets_AK5PF_gen_phi = 0;
   jets_AK5PF_gen_mass = 0;
   jets_AK5PF_gen_Energy = 0;
   jets_AK5PF_gen_Id = 0;
   jets_AK5PF_gen_motherID = 0;
   jets_AK5PF_gen_threeCharge = 0;
   jets_AK5PF_partonFlavour = 0;
   jets_AK5PF_btag_TC_highPur = 0;
   jets_AK5PF_btag_TC_highEff = 0;
   jets_AK5PF_btag_jetProb = 0;
   jets_AK5PF_btag_jetBProb = 0;
   jets_AK5PF_btag_softEle = 0;
   jets_AK5PF_btag_softMuon = 0;
   jets_AK5PF_btag_secVertexHighPur = 0;
   jets_AK5PF_btag_secVertexHighEff = 0;
   jets_AK5PF_btag_secVertexCombined = 0;
   jets_AK5PF_chgEmE = 0;
   jets_AK5PF_chgHadE = 0;
   jets_AK5PF_chgMuE = 0;
   jets_AK5PF_chg_Mult = 0;
   jets_AK5PF_neutralEmE = 0;
   jets_AK5PF_neutralHadE = 0;
   jets_AK5PF_neutral_Mult = 0;
   jets_AK5PF_mu_Mult = 0;
   jets_AK5PF_emf = 0;
   jets_AK5PF_ehf = 0;
   jets_AK5PF_n60 = 0;
   jets_AK5PF_n90 = 0;
   jets_AK5PF_etaetaMoment = 0;
   jets_AK5PF_etaphiMoment = 0;
   jets_AK5PF_phiphiMoment = 0;
   jets_AK5PF_n90Hits = 0;
   jets_AK5PF_fHPD = 0;
   jets_AK5PF_fRBX = 0;
   jets_AK5PF_hitsInN90 = 0;
   jets_AK5PF_nECALTowers = 0;
   jets_AK5PF_nHCALTowers = 0;
   jets_AK5PF_fSubDetector1 = 0;
   jets_AK5PF_fSubDetector2 = 0;
   jets_AK5PF_fSubDetector3 = 0;
   jets_AK5PF_fSubDetector4 = 0;
   jets_AK5PF_area = 0;
   jets_AK5PF_mass = 0;
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
   mus_isTrackerMuon = 0;
   mus_isStandAloneMuon = 0;
   mus_isCaloMuon = 0;
   mus_isGlobalMuon = 0;
   mus_isElectron = 0;
   mus_isConvertedPhoton = 0;
   mus_isPhoton = 0;
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
   pfmets_et = 0;
   pfmets_phi = 0;
   pfmets_ex = 0;
   pfmets_ey = 0;
   pfmets_gen_et = 0;
   pfmets_gen_phi = 0;
   pfmets_sign = 0;
   pfmets_sumEt = 0;
   pfmets_unCPhi = 0;
   pfmets_unCPt = 0;
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
   photons_maxEnergyXtal = 0;
   photons_e1x5 = 0;
   photons_e2x5 = 0;
   photons_e3x3 = 0;
   photons_e5x5 = 0;
   photons_sigmaEtaEta = 0;
   photons_sigmaIetaIeta = 0;
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
   pv_chi2 = 0;
   pv_ndof = 0;
   pv_tracksSize = 0;
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
   taus_taNC = 0;
   taus_byIsoUsingLeadingPi = 0;
   taus_tkIsoUsingLeadingPi = 0;
   taus_ecalIsoUsingLeadingPi = 0;
   taus_signalPFChargedHadrCandsSize = 0;
   taus_muDecision = 0;
   taus_Nprongs = 0;
   tcmets_et = 0;
   tcmets_phi = 0;
   tcmets_ex = 0;
   tcmets_ey = 0;
   tcmets_sumEt = 0;
   tracks_chi2 = 0;
   tracks_trkExptHitsInner = 0;
   tracks_trkExptHitsOuter = 0;
   tracks_trks_nlayerslost = 0;
   tracks_trks_nlayers = 0;
   tracks_trksvalidpixelhits = 0;
   tracks_trkslostpixelhits = 0;
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
   tracks_outerPx = 0;
   tracks_outerPy = 0;
   tracks_outerPz = 0;
   tracks_algo = 0;
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
   fChain->SetBranchAddress("NconvTk", &NconvTk, &b_NconvTk);
   fChain->SetBranchAddress("convTk_tracksPin_x1", &convTk_tracksPin_x1, &b_convTk_tracksPin_x1);
   fChain->SetBranchAddress("convTk_tracksPin_x2", &convTk_tracksPin_x2, &b_convTk_tracksPin_x2);
   fChain->SetBranchAddress("convTk_tracksPin_y1", &convTk_tracksPin_y1, &b_convTk_tracksPin_y1);
   fChain->SetBranchAddress("convTk_tracksPin_y2", &convTk_tracksPin_y2, &b_convTk_tracksPin_y2);
   fChain->SetBranchAddress("convTk_tracksPin_z1", &convTk_tracksPin_z1, &b_convTk_tracksPin_z1);
   fChain->SetBranchAddress("convTk_tracksPin_z2", &convTk_tracksPin_z2, &b_convTk_tracksPin_z2);
   fChain->SetBranchAddress("convTk_tracksPout_x1", &convTk_tracksPout_x1, &b_convTk_tracksPout_x1);
   fChain->SetBranchAddress("convTk_tracksPout_x2", &convTk_tracksPout_x2, &b_convTk_tracksPout_x2);
   fChain->SetBranchAddress("convTk_tracksPout_y1", &convTk_tracksPout_y1, &b_convTk_tracksPout_y1);
   fChain->SetBranchAddress("convTk_tracksPout_y2", &convTk_tracksPout_y2, &b_convTk_tracksPout_y2);
   fChain->SetBranchAddress("convTk_tracksPout_z1", &convTk_tracksPout_z1, &b_convTk_tracksPout_z1);
   fChain->SetBranchAddress("convTk_tracksPout_z2", &convTk_tracksPout_z2, &b_convTk_tracksPout_z2);
   fChain->SetBranchAddress("convTk_ecalImpactPosition_x1", &convTk_ecalImpactPosition_x1, &b_convTk_ecalImpactPosition_x1);
   fChain->SetBranchAddress("convTk_ecalImpactPosition_x2", &convTk_ecalImpactPosition_x2, &b_convTk_ecalImpactPosition_x2);
   fChain->SetBranchAddress("convTk_ecalImpactPosition_y1", &convTk_ecalImpactPosition_y1, &b_convTk_ecalImpactPosition_y1);
   fChain->SetBranchAddress("convTk_ecalImpactPosition_y2", &convTk_ecalImpactPosition_y2, &b_convTk_ecalImpactPosition_y2);
   fChain->SetBranchAddress("convTk_ecalImpactPosition_z1", &convTk_ecalImpactPosition_z1, &b_convTk_ecalImpactPosition_z1);
   fChain->SetBranchAddress("convTk_ecalImpactPosition_z2", &convTk_ecalImpactPosition_z2, &b_convTk_ecalImpactPosition_z2);
   fChain->SetBranchAddress("convTk_tracksInnerPosition_x1", &convTk_tracksInnerPosition_x1, &b_convTk_tracksInnerPosition_x1);
   fChain->SetBranchAddress("convTk_tracksInnerPosition_x2", &convTk_tracksInnerPosition_x2, &b_convTk_tracksInnerPosition_x2);
   fChain->SetBranchAddress("convTk_tracksInnerPosition_y1", &convTk_tracksInnerPosition_y1, &b_convTk_tracksInnerPosition_y1);
   fChain->SetBranchAddress("convTk_tracksInnerPosition_y2", &convTk_tracksInnerPosition_y2, &b_convTk_tracksInnerPosition_y2);
   fChain->SetBranchAddress("convTk_tracksInnerPosition_z1", &convTk_tracksInnerPosition_z1, &b_convTk_tracksInnerPosition_z1);
   fChain->SetBranchAddress("convTk_tracksInnerPosition_z2", &convTk_tracksInnerPosition_z2, &b_convTk_tracksInnerPosition_z2);
   fChain->SetBranchAddress("convTk_tracks_px1", &convTk_tracks_px1, &b_convTk_tracks_px1);
   fChain->SetBranchAddress("convTk_tracks_py1", &convTk_tracks_py1, &b_convTk_tracks_py1);
   fChain->SetBranchAddress("convTk_tracks_pz1", &convTk_tracks_pz1, &b_convTk_tracks_pz1);
   fChain->SetBranchAddress("convTk_tracks_phi1", &convTk_tracks_phi1, &b_convTk_tracks_phi1);
   fChain->SetBranchAddress("convTk_tracks_d01", &convTk_tracks_d01, &b_convTk_tracks_d01);
   fChain->SetBranchAddress("convTk_tracks_dz1", &convTk_tracks_dz1, &b_convTk_tracks_dz1);
   fChain->SetBranchAddress("convTk_tracks_charge1", &convTk_tracks_charge1, &b_convTk_tracks_charge1);
   fChain->SetBranchAddress("convTk_tracks_px2", &convTk_tracks_px2, &b_convTk_tracks_px2);
   fChain->SetBranchAddress("convTk_tracks_py2", &convTk_tracks_py2, &b_convTk_tracks_py2);
   fChain->SetBranchAddress("convTk_tracks_pz2", &convTk_tracks_pz2, &b_convTk_tracks_pz2);
   fChain->SetBranchAddress("convTk_tracks_phi2", &convTk_tracks_phi2, &b_convTk_tracks_phi2);
   fChain->SetBranchAddress("convTk_tracks_d02", &convTk_tracks_d02, &b_convTk_tracks_d02);
   fChain->SetBranchAddress("convTk_tracks_dz2", &convTk_tracks_dz2, &b_convTk_tracks_dz2);
   fChain->SetBranchAddress("convTk_tracks_charge2", &convTk_tracks_charge2, &b_convTk_tracks_charge2);
   fChain->SetBranchAddress("convTk_bcMatchingWithTracks_x1", &convTk_bcMatchingWithTracks_x1, &b_convTk_bcMatchingWithTracks_x1);
   fChain->SetBranchAddress("convTk_bcMatchingWithTracks_y1", &convTk_bcMatchingWithTracks_y1, &b_convTk_bcMatchingWithTracks_y1);
   fChain->SetBranchAddress("convTk_bcMatchingWithTracks_z1", &convTk_bcMatchingWithTracks_z1, &b_convTk_bcMatchingWithTracks_z1);
   fChain->SetBranchAddress("convTk_bcMatchingWithTracks_x2", &convTk_bcMatchingWithTracks_x2, &b_convTk_bcMatchingWithTracks_x2);
   fChain->SetBranchAddress("convTk_bcMatchingWithTracks_y2", &convTk_bcMatchingWithTracks_y2, &b_convTk_bcMatchingWithTracks_y2);
   fChain->SetBranchAddress("convTk_bcMatchingWithTracks_z2", &convTk_bcMatchingWithTracks_z2, &b_convTk_bcMatchingWithTracks_z2);
   fChain->SetBranchAddress("convTk_bcMatchingWithTracks_energy1", &convTk_bcMatchingWithTracks_energy1, &b_convTk_bcMatchingWithTracks_energy1);
   fChain->SetBranchAddress("convTk_bcMatchingWithTracks_energy2", &convTk_bcMatchingWithTracks_energy2, &b_convTk_bcMatchingWithTracks_energy2);
   fChain->SetBranchAddress("convTk_EoverP", &convTk_EoverP, &b_convTk_EoverP);
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
   fChain->SetBranchAddress("els_tk_pz", &els_tk_pz, &b_els_tk_pz);
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
   fChain->SetBranchAddress("els_n_inner_layer", &els_n_inner_layer, &b_els_n_inner_layer);
   fChain->SetBranchAddress("els_n_outer_layer", &els_n_outer_layer, &b_els_n_outer_layer);
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
   fChain->SetBranchAddress("NhcalNoiseRBX", &NhcalNoiseRBX, &b_NhcalNoiseRBX);
   fChain->SetBranchAddress("hcalNoiseRBX_idnumber", &hcalNoiseRBX_idnumber, &b_hcalNoiseRBX_idnumber);
   fChain->SetBranchAddress("hcalNoiseRBX_allChargeTotal", &hcalNoiseRBX_allChargeTotal, &b_hcalNoiseRBX_allChargeTotal);
   fChain->SetBranchAddress("hcalNoiseRBX_allChargeHighest2TS", &hcalNoiseRBX_allChargeHighest2TS, &b_hcalNoiseRBX_allChargeHighest2TS);
   fChain->SetBranchAddress("hcalNoiseRBX_allChargeHighest3TS", &hcalNoiseRBX_allChargeHighest3TS, &b_hcalNoiseRBX_allChargeHighest3TS);
   fChain->SetBranchAddress("hcalNoiseRBX_totalZeros", &hcalNoiseRBX_totalZeros, &b_hcalNoiseRBX_totalZeros);
   fChain->SetBranchAddress("hcalNoiseRBX_maxZeros", &hcalNoiseRBX_maxZeros, &b_hcalNoiseRBX_maxZeros);
   fChain->SetBranchAddress("hcalNoiseRBX_recHitEnergy", &hcalNoiseRBX_recHitEnergy, &b_hcalNoiseRBX_recHitEnergy);
   fChain->SetBranchAddress("hcalNoiseRBX_minRecHitTime", &hcalNoiseRBX_minRecHitTime, &b_hcalNoiseRBX_minRecHitTime);
   fChain->SetBranchAddress("hcalNoiseRBX_maxRecHitTime", &hcalNoiseRBX_maxRecHitTime, &b_hcalNoiseRBX_maxRecHitTime);
   fChain->SetBranchAddress("hcalNoiseRBX_numRecHits", &hcalNoiseRBX_numRecHits, &b_hcalNoiseRBX_numRecHits);
   fChain->SetBranchAddress("hcalNoiseRBX_caloTowerHadE", &hcalNoiseRBX_caloTowerHadE, &b_hcalNoiseRBX_caloTowerHadE);
   fChain->SetBranchAddress("hcalNoiseRBX_caloTowerEmE", &hcalNoiseRBX_caloTowerEmE, &b_hcalNoiseRBX_caloTowerEmE);
   fChain->SetBranchAddress("hcalNoiseRBX_caloTowerTotalE", &hcalNoiseRBX_caloTowerTotalE, &b_hcalNoiseRBX_caloTowerTotalE);
   fChain->SetBranchAddress("hcalNoiseRBX_caloTowerEmFraction", &hcalNoiseRBX_caloTowerEmFraction, &b_hcalNoiseRBX_caloTowerEmFraction);
   fChain->SetBranchAddress("NhcalNoiseSummary", &NhcalNoiseSummary, &b_NhcalNoiseSummary);
   fChain->SetBranchAddress("hcalNoiseSummary_passLooseNoiseFilter", &hcalNoiseSummary_passLooseNoiseFilter, &b_hcalNoiseSummary_passLooseNoiseFilter);
   fChain->SetBranchAddress("hcalNoiseSummary_passTightNoiseFilter", &hcalNoiseSummary_passTightNoiseFilter, &b_hcalNoiseSummary_passTightNoiseFilter);
   fChain->SetBranchAddress("hcalNoiseSummary_passHighLevelNoiseFilter", &hcalNoiseSummary_passHighLevelNoiseFilter, &b_hcalNoiseSummary_passHighLevelNoiseFilter);
   fChain->SetBranchAddress("hcalNoiseSummary_noiseFilterStatus", &hcalNoiseSummary_noiseFilterStatus, &b_hcalNoiseSummary_noiseFilterStatus);
   fChain->SetBranchAddress("hcalNoiseSummary_noiseType", &hcalNoiseSummary_noiseType, &b_hcalNoiseSummary_noiseType);
   fChain->SetBranchAddress("hcalNoiseSummary_eventEMEnergy", &hcalNoiseSummary_eventEMEnergy, &b_hcalNoiseSummary_eventEMEnergy);
   fChain->SetBranchAddress("hcalNoiseSummary_eventHadEnergy", &hcalNoiseSummary_eventHadEnergy, &b_hcalNoiseSummary_eventHadEnergy);
   fChain->SetBranchAddress("hcalNoiseSummary_eventTrackEnergy", &hcalNoiseSummary_eventTrackEnergy, &b_hcalNoiseSummary_eventTrackEnergy);
   fChain->SetBranchAddress("hcalNoiseSummary_eventEMFraction", &hcalNoiseSummary_eventEMFraction, &b_hcalNoiseSummary_eventEMFraction);
   fChain->SetBranchAddress("hcalNoiseSummary_eventChargeFraction", &hcalNoiseSummary_eventChargeFraction, &b_hcalNoiseSummary_eventChargeFraction);
   fChain->SetBranchAddress("hcalNoiseSummary_min10GeVHitTime", &hcalNoiseSummary_min10GeVHitTime, &b_hcalNoiseSummary_min10GeVHitTime);
   fChain->SetBranchAddress("hcalNoiseSummary_max10GeVHitTime", &hcalNoiseSummary_max10GeVHitTime, &b_hcalNoiseSummary_max10GeVHitTime);
   fChain->SetBranchAddress("hcalNoiseSummary_rms10GeVHitTime", &hcalNoiseSummary_rms10GeVHitTime, &b_hcalNoiseSummary_rms10GeVHitTime);
   fChain->SetBranchAddress("hcalNoiseSummary_min25GeVHitTime", &hcalNoiseSummary_min25GeVHitTime, &b_hcalNoiseSummary_min25GeVHitTime);
   fChain->SetBranchAddress("hcalNoiseSummary_max25GeVHitTime", &hcalNoiseSummary_max25GeVHitTime, &b_hcalNoiseSummary_max25GeVHitTime);
   fChain->SetBranchAddress("hcalNoiseSummary_rms25GeVHitTime", &hcalNoiseSummary_rms25GeVHitTime, &b_hcalNoiseSummary_rms25GeVHitTime);
   fChain->SetBranchAddress("hcalNoiseSummary_num10GeVHits", &hcalNoiseSummary_num10GeVHits, &b_hcalNoiseSummary_num10GeVHits);
   fChain->SetBranchAddress("hcalNoiseSummary_num25GeVHits", &hcalNoiseSummary_num25GeVHits, &b_hcalNoiseSummary_num25GeVHits);
   fChain->SetBranchAddress("hcalNoiseSummary_minE2TS", &hcalNoiseSummary_minE2TS, &b_hcalNoiseSummary_minE2TS);
   fChain->SetBranchAddress("hcalNoiseSummary_minE10TS", &hcalNoiseSummary_minE10TS, &b_hcalNoiseSummary_minE10TS);
   fChain->SetBranchAddress("hcalNoiseSummary_minE2Over10TS", &hcalNoiseSummary_minE2Over10TS, &b_hcalNoiseSummary_minE2Over10TS);
   fChain->SetBranchAddress("hcalNoiseSummary_maxZeros", &hcalNoiseSummary_maxZeros, &b_hcalNoiseSummary_maxZeros);
   fChain->SetBranchAddress("hcalNoiseSummary_maxHPDHits", &hcalNoiseSummary_maxHPDHits, &b_hcalNoiseSummary_maxHPDHits);
   fChain->SetBranchAddress("hcalNoiseSummary_maxRBXHits", &hcalNoiseSummary_maxRBXHits, &b_hcalNoiseSummary_maxRBXHits);
   fChain->SetBranchAddress("hcalNoiseSummary_minHPDEMF", &hcalNoiseSummary_minHPDEMF, &b_hcalNoiseSummary_minHPDEMF);
   fChain->SetBranchAddress("hcalNoiseSummary_minRBXEMF", &hcalNoiseSummary_minRBXEMF, &b_hcalNoiseSummary_minRBXEMF);
   fChain->SetBranchAddress("hcalNoiseSummary_numProblematicRBXs", &hcalNoiseSummary_numProblematicRBXs, &b_hcalNoiseSummary_numProblematicRBXs);
   fChain->SetBranchAddress("NhybridBBC", &NhybridBBC, &b_NhybridBBC);
   fChain->SetBranchAddress("hybridBBC_energy", &hybridBBC_energy, &b_hybridBBC_energy);
   fChain->SetBranchAddress("hybridBBC_x", &hybridBBC_x, &b_hybridBBC_x);
   fChain->SetBranchAddress("hybridBBC_y", &hybridBBC_y, &b_hybridBBC_y);
   fChain->SetBranchAddress("hybridBBC_z", &hybridBBC_z, &b_hybridBBC_z);
   fChain->SetBranchAddress("hybridBBC_rho", &hybridBBC_rho, &b_hybridBBC_rho);
   fChain->SetBranchAddress("hybridBBC_phi", &hybridBBC_phi, &b_hybridBBC_phi);
   fChain->SetBranchAddress("hybridBBC_eta", &hybridBBC_eta, &b_hybridBBC_eta);
   fChain->SetBranchAddress("hybridBBC_theta", &hybridBBC_theta, &b_hybridBBC_theta);
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
   fChain->SetBranchAddress("jets_AK5_btag_secVertexHighPur", &jets_AK5_btag_secVertexHighPur, &b_jets_AK5_btag_secVertexHighPur);
   fChain->SetBranchAddress("jets_AK5_btag_secVertexHighEff", &jets_AK5_btag_secVertexHighEff, &b_jets_AK5_btag_secVertexHighEff);
   fChain->SetBranchAddress("jets_AK5_btag_secVertexCombined", &jets_AK5_btag_secVertexCombined, &b_jets_AK5_btag_secVertexCombined);
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
   fChain->SetBranchAddress("jets_AK5_etaetaMoment", &jets_AK5_etaetaMoment, &b_jets_AK5_etaetaMoment);
   fChain->SetBranchAddress("jets_AK5_etaphiMoment", &jets_AK5_etaphiMoment, &b_jets_AK5_etaphiMoment);
   fChain->SetBranchAddress("jets_AK5_phiphiMoment", &jets_AK5_phiphiMoment, &b_jets_AK5_phiphiMoment);
   fChain->SetBranchAddress("jets_AK5_n90Hits", &jets_AK5_n90Hits, &b_jets_AK5_n90Hits);
   fChain->SetBranchAddress("jets_AK5_fHPD", &jets_AK5_fHPD, &b_jets_AK5_fHPD);
   fChain->SetBranchAddress("jets_AK5_fRBX", &jets_AK5_fRBX, &b_jets_AK5_fRBX);
   fChain->SetBranchAddress("jets_AK5_hitsInN90", &jets_AK5_hitsInN90, &b_jets_AK5_hitsInN90);
   fChain->SetBranchAddress("jets_AK5_nECALTowers", &jets_AK5_nECALTowers, &b_jets_AK5_nECALTowers);
   fChain->SetBranchAddress("jets_AK5_nHCALTowers", &jets_AK5_nHCALTowers, &b_jets_AK5_nHCALTowers);
   fChain->SetBranchAddress("jets_AK5_fSubDetector1", &jets_AK5_fSubDetector1, &b_jets_AK5_fSubDetector1);
   fChain->SetBranchAddress("jets_AK5_fSubDetector2", &jets_AK5_fSubDetector2, &b_jets_AK5_fSubDetector2);
   fChain->SetBranchAddress("jets_AK5_fSubDetector3", &jets_AK5_fSubDetector3, &b_jets_AK5_fSubDetector3);
   fChain->SetBranchAddress("jets_AK5_fSubDetector4", &jets_AK5_fSubDetector4, &b_jets_AK5_fSubDetector4);
   fChain->SetBranchAddress("jets_AK5_area", &jets_AK5_area, &b_jets_AK5_area);
   fChain->SetBranchAddress("jets_AK5_mass", &jets_AK5_mass, &b_jets_AK5_mass);
   fChain->SetBranchAddress("Njets_AK5JPT", &Njets_AK5JPT, &b_Njets_AK5JPT);
   fChain->SetBranchAddress("jets_AK5JPT_energy", &jets_AK5JPT_energy, &b_jets_AK5JPT_energy);
   fChain->SetBranchAddress("jets_AK5JPT_et", &jets_AK5JPT_et, &b_jets_AK5JPT_et);
   fChain->SetBranchAddress("jets_AK5JPT_eta", &jets_AK5JPT_eta, &b_jets_AK5JPT_eta);
   fChain->SetBranchAddress("jets_AK5JPT_phi", &jets_AK5JPT_phi, &b_jets_AK5JPT_phi);
   fChain->SetBranchAddress("jets_AK5JPT_pt", &jets_AK5JPT_pt, &b_jets_AK5JPT_pt);
   fChain->SetBranchAddress("jets_AK5JPT_px", &jets_AK5JPT_px, &b_jets_AK5JPT_px);
   fChain->SetBranchAddress("jets_AK5JPT_py", &jets_AK5JPT_py, &b_jets_AK5JPT_py);
   fChain->SetBranchAddress("jets_AK5JPT_pz", &jets_AK5JPT_pz, &b_jets_AK5JPT_pz);
   fChain->SetBranchAddress("jets_AK5JPT_status", &jets_AK5JPT_status, &b_jets_AK5JPT_status);
   fChain->SetBranchAddress("jets_AK5JPT_theta", &jets_AK5JPT_theta, &b_jets_AK5JPT_theta);
   fChain->SetBranchAddress("jets_AK5JPT_parton_Id", &jets_AK5JPT_parton_Id, &b_jets_AK5JPT_parton_Id);
   fChain->SetBranchAddress("jets_AK5JPT_parton_motherId", &jets_AK5JPT_parton_motherId, &b_jets_AK5JPT_parton_motherId);
   fChain->SetBranchAddress("jets_AK5JPT_parton_pt", &jets_AK5JPT_parton_pt, &b_jets_AK5JPT_parton_pt);
   fChain->SetBranchAddress("jets_AK5JPT_parton_phi", &jets_AK5JPT_parton_phi, &b_jets_AK5JPT_parton_phi);
   fChain->SetBranchAddress("jets_AK5JPT_parton_eta", &jets_AK5JPT_parton_eta, &b_jets_AK5JPT_parton_eta);
   fChain->SetBranchAddress("jets_AK5JPT_parton_Energy", &jets_AK5JPT_parton_Energy, &b_jets_AK5JPT_parton_Energy);
   fChain->SetBranchAddress("jets_AK5JPT_parton_mass", &jets_AK5JPT_parton_mass, &b_jets_AK5JPT_parton_mass);
   fChain->SetBranchAddress("jets_AK5JPT_gen_et", &jets_AK5JPT_gen_et, &b_jets_AK5JPT_gen_et);
   fChain->SetBranchAddress("jets_AK5JPT_gen_pt", &jets_AK5JPT_gen_pt, &b_jets_AK5JPT_gen_pt);
   fChain->SetBranchAddress("jets_AK5JPT_gen_eta", &jets_AK5JPT_gen_eta, &b_jets_AK5JPT_gen_eta);
   fChain->SetBranchAddress("jets_AK5JPT_gen_phi", &jets_AK5JPT_gen_phi, &b_jets_AK5JPT_gen_phi);
   fChain->SetBranchAddress("jets_AK5JPT_gen_mass", &jets_AK5JPT_gen_mass, &b_jets_AK5JPT_gen_mass);
   fChain->SetBranchAddress("jets_AK5JPT_gen_Energy", &jets_AK5JPT_gen_Energy, &b_jets_AK5JPT_gen_Energy);
   fChain->SetBranchAddress("jets_AK5JPT_gen_Id", &jets_AK5JPT_gen_Id, &b_jets_AK5JPT_gen_Id);
   fChain->SetBranchAddress("jets_AK5JPT_gen_motherID", &jets_AK5JPT_gen_motherID, &b_jets_AK5JPT_gen_motherID);
   fChain->SetBranchAddress("jets_AK5JPT_gen_threeCharge", &jets_AK5JPT_gen_threeCharge, &b_jets_AK5JPT_gen_threeCharge);
   fChain->SetBranchAddress("jets_AK5JPT_partonFlavour", &jets_AK5JPT_partonFlavour, &b_jets_AK5JPT_partonFlavour);
   fChain->SetBranchAddress("jets_AK5JPT_btag_TC_highPur", &jets_AK5JPT_btag_TC_highPur, &b_jets_AK5JPT_btag_TC_highPur);
   fChain->SetBranchAddress("jets_AK5JPT_btag_TC_highEff", &jets_AK5JPT_btag_TC_highEff, &b_jets_AK5JPT_btag_TC_highEff);
   fChain->SetBranchAddress("jets_AK5JPT_btag_jetProb", &jets_AK5JPT_btag_jetProb, &b_jets_AK5JPT_btag_jetProb);
   fChain->SetBranchAddress("jets_AK5JPT_btag_jetBProb", &jets_AK5JPT_btag_jetBProb, &b_jets_AK5JPT_btag_jetBProb);
   fChain->SetBranchAddress("jets_AK5JPT_btag_softEle", &jets_AK5JPT_btag_softEle, &b_jets_AK5JPT_btag_softEle);
   fChain->SetBranchAddress("jets_AK5JPT_btag_softMuon", &jets_AK5JPT_btag_softMuon, &b_jets_AK5JPT_btag_softMuon);
   fChain->SetBranchAddress("jets_AK5JPT_btag_secVertexHighPur", &jets_AK5JPT_btag_secVertexHighPur, &b_jets_AK5JPT_btag_secVertexHighPur);
   fChain->SetBranchAddress("jets_AK5JPT_btag_secVertexHighEff", &jets_AK5JPT_btag_secVertexHighEff, &b_jets_AK5JPT_btag_secVertexHighEff);
   fChain->SetBranchAddress("jets_AK5JPT_btag_secVertexCombined", &jets_AK5JPT_btag_secVertexCombined, &b_jets_AK5JPT_btag_secVertexCombined);
   fChain->SetBranchAddress("jets_AK5JPT_chgEmE", &jets_AK5JPT_chgEmE, &b_jets_AK5JPT_chgEmE);
   fChain->SetBranchAddress("jets_AK5JPT_chgHadE", &jets_AK5JPT_chgHadE, &b_jets_AK5JPT_chgHadE);
   fChain->SetBranchAddress("jets_AK5JPT_chgMuE", &jets_AK5JPT_chgMuE, &b_jets_AK5JPT_chgMuE);
   fChain->SetBranchAddress("jets_AK5JPT_chg_Mult", &jets_AK5JPT_chg_Mult, &b_jets_AK5JPT_chg_Mult);
   fChain->SetBranchAddress("jets_AK5JPT_neutralEmE", &jets_AK5JPT_neutralEmE, &b_jets_AK5JPT_neutralEmE);
   fChain->SetBranchAddress("jets_AK5JPT_neutralHadE", &jets_AK5JPT_neutralHadE, &b_jets_AK5JPT_neutralHadE);
   fChain->SetBranchAddress("jets_AK5JPT_neutral_Mult", &jets_AK5JPT_neutral_Mult, &b_jets_AK5JPT_neutral_Mult);
   fChain->SetBranchAddress("jets_AK5JPT_mu_Mult", &jets_AK5JPT_mu_Mult, &b_jets_AK5JPT_mu_Mult);
   fChain->SetBranchAddress("jets_AK5JPT_emf", &jets_AK5JPT_emf, &b_jets_AK5JPT_emf);
   fChain->SetBranchAddress("jets_AK5JPT_ehf", &jets_AK5JPT_ehf, &b_jets_AK5JPT_ehf);
   fChain->SetBranchAddress("jets_AK5JPT_n60", &jets_AK5JPT_n60, &b_jets_AK5JPT_n60);
   fChain->SetBranchAddress("jets_AK5JPT_n90", &jets_AK5JPT_n90, &b_jets_AK5JPT_n90);
   fChain->SetBranchAddress("jets_AK5JPT_etaetaMoment", &jets_AK5JPT_etaetaMoment, &b_jets_AK5JPT_etaetaMoment);
   fChain->SetBranchAddress("jets_AK5JPT_etaphiMoment", &jets_AK5JPT_etaphiMoment, &b_jets_AK5JPT_etaphiMoment);
   fChain->SetBranchAddress("jets_AK5JPT_phiphiMoment", &jets_AK5JPT_phiphiMoment, &b_jets_AK5JPT_phiphiMoment);
   fChain->SetBranchAddress("jets_AK5JPT_n90Hits", &jets_AK5JPT_n90Hits, &b_jets_AK5JPT_n90Hits);
   fChain->SetBranchAddress("jets_AK5JPT_fHPD", &jets_AK5JPT_fHPD, &b_jets_AK5JPT_fHPD);
   fChain->SetBranchAddress("jets_AK5JPT_fRBX", &jets_AK5JPT_fRBX, &b_jets_AK5JPT_fRBX);
   fChain->SetBranchAddress("jets_AK5JPT_hitsInN90", &jets_AK5JPT_hitsInN90, &b_jets_AK5JPT_hitsInN90);
   fChain->SetBranchAddress("jets_AK5JPT_nECALTowers", &jets_AK5JPT_nECALTowers, &b_jets_AK5JPT_nECALTowers);
   fChain->SetBranchAddress("jets_AK5JPT_nHCALTowers", &jets_AK5JPT_nHCALTowers, &b_jets_AK5JPT_nHCALTowers);
   fChain->SetBranchAddress("jets_AK5JPT_fSubDetector1", &jets_AK5JPT_fSubDetector1, &b_jets_AK5JPT_fSubDetector1);
   fChain->SetBranchAddress("jets_AK5JPT_fSubDetector2", &jets_AK5JPT_fSubDetector2, &b_jets_AK5JPT_fSubDetector2);
   fChain->SetBranchAddress("jets_AK5JPT_fSubDetector3", &jets_AK5JPT_fSubDetector3, &b_jets_AK5JPT_fSubDetector3);
   fChain->SetBranchAddress("jets_AK5JPT_fSubDetector4", &jets_AK5JPT_fSubDetector4, &b_jets_AK5JPT_fSubDetector4);
   fChain->SetBranchAddress("jets_AK5JPT_area", &jets_AK5JPT_area, &b_jets_AK5JPT_area);
   fChain->SetBranchAddress("jets_AK5JPT_mass", &jets_AK5JPT_mass, &b_jets_AK5JPT_mass);
   fChain->SetBranchAddress("Njets_AK5PF", &Njets_AK5PF, &b_Njets_AK5PF);
   fChain->SetBranchAddress("jets_AK5PF_energy", &jets_AK5PF_energy, &b_jets_AK5PF_energy);
   fChain->SetBranchAddress("jets_AK5PF_et", &jets_AK5PF_et, &b_jets_AK5PF_et);
   fChain->SetBranchAddress("jets_AK5PF_eta", &jets_AK5PF_eta, &b_jets_AK5PF_eta);
   fChain->SetBranchAddress("jets_AK5PF_phi", &jets_AK5PF_phi, &b_jets_AK5PF_phi);
   fChain->SetBranchAddress("jets_AK5PF_pt", &jets_AK5PF_pt, &b_jets_AK5PF_pt);
   fChain->SetBranchAddress("jets_AK5PF_px", &jets_AK5PF_px, &b_jets_AK5PF_px);
   fChain->SetBranchAddress("jets_AK5PF_py", &jets_AK5PF_py, &b_jets_AK5PF_py);
   fChain->SetBranchAddress("jets_AK5PF_pz", &jets_AK5PF_pz, &b_jets_AK5PF_pz);
   fChain->SetBranchAddress("jets_AK5PF_status", &jets_AK5PF_status, &b_jets_AK5PF_status);
   fChain->SetBranchAddress("jets_AK5PF_theta", &jets_AK5PF_theta, &b_jets_AK5PF_theta);
   fChain->SetBranchAddress("jets_AK5PF_parton_Id", &jets_AK5PF_parton_Id, &b_jets_AK5PF_parton_Id);
   fChain->SetBranchAddress("jets_AK5PF_parton_motherId", &jets_AK5PF_parton_motherId, &b_jets_AK5PF_parton_motherId);
   fChain->SetBranchAddress("jets_AK5PF_parton_pt", &jets_AK5PF_parton_pt, &b_jets_AK5PF_parton_pt);
   fChain->SetBranchAddress("jets_AK5PF_parton_phi", &jets_AK5PF_parton_phi, &b_jets_AK5PF_parton_phi);
   fChain->SetBranchAddress("jets_AK5PF_parton_eta", &jets_AK5PF_parton_eta, &b_jets_AK5PF_parton_eta);
   fChain->SetBranchAddress("jets_AK5PF_parton_Energy", &jets_AK5PF_parton_Energy, &b_jets_AK5PF_parton_Energy);
   fChain->SetBranchAddress("jets_AK5PF_parton_mass", &jets_AK5PF_parton_mass, &b_jets_AK5PF_parton_mass);
   fChain->SetBranchAddress("jets_AK5PF_gen_et", &jets_AK5PF_gen_et, &b_jets_AK5PF_gen_et);
   fChain->SetBranchAddress("jets_AK5PF_gen_pt", &jets_AK5PF_gen_pt, &b_jets_AK5PF_gen_pt);
   fChain->SetBranchAddress("jets_AK5PF_gen_eta", &jets_AK5PF_gen_eta, &b_jets_AK5PF_gen_eta);
   fChain->SetBranchAddress("jets_AK5PF_gen_phi", &jets_AK5PF_gen_phi, &b_jets_AK5PF_gen_phi);
   fChain->SetBranchAddress("jets_AK5PF_gen_mass", &jets_AK5PF_gen_mass, &b_jets_AK5PF_gen_mass);
   fChain->SetBranchAddress("jets_AK5PF_gen_Energy", &jets_AK5PF_gen_Energy, &b_jets_AK5PF_gen_Energy);
   fChain->SetBranchAddress("jets_AK5PF_gen_Id", &jets_AK5PF_gen_Id, &b_jets_AK5PF_gen_Id);
   fChain->SetBranchAddress("jets_AK5PF_gen_motherID", &jets_AK5PF_gen_motherID, &b_jets_AK5PF_gen_motherID);
   fChain->SetBranchAddress("jets_AK5PF_gen_threeCharge", &jets_AK5PF_gen_threeCharge, &b_jets_AK5PF_gen_threeCharge);
   fChain->SetBranchAddress("jets_AK5PF_partonFlavour", &jets_AK5PF_partonFlavour, &b_jets_AK5PF_partonFlavour);
   fChain->SetBranchAddress("jets_AK5PF_btag_TC_highPur", &jets_AK5PF_btag_TC_highPur, &b_jets_AK5PF_btag_TC_highPur);
   fChain->SetBranchAddress("jets_AK5PF_btag_TC_highEff", &jets_AK5PF_btag_TC_highEff, &b_jets_AK5PF_btag_TC_highEff);
   fChain->SetBranchAddress("jets_AK5PF_btag_jetProb", &jets_AK5PF_btag_jetProb, &b_jets_AK5PF_btag_jetProb);
   fChain->SetBranchAddress("jets_AK5PF_btag_jetBProb", &jets_AK5PF_btag_jetBProb, &b_jets_AK5PF_btag_jetBProb);
   fChain->SetBranchAddress("jets_AK5PF_btag_softEle", &jets_AK5PF_btag_softEle, &b_jets_AK5PF_btag_softEle);
   fChain->SetBranchAddress("jets_AK5PF_btag_softMuon", &jets_AK5PF_btag_softMuon, &b_jets_AK5PF_btag_softMuon);
   fChain->SetBranchAddress("jets_AK5PF_btag_secVertexHighPur", &jets_AK5PF_btag_secVertexHighPur, &b_jets_AK5PF_btag_secVertexHighPur);
   fChain->SetBranchAddress("jets_AK5PF_btag_secVertexHighEff", &jets_AK5PF_btag_secVertexHighEff, &b_jets_AK5PF_btag_secVertexHighEff);
   fChain->SetBranchAddress("jets_AK5PF_btag_secVertexCombined", &jets_AK5PF_btag_secVertexCombined, &b_jets_AK5PF_btag_secVertexCombined);
   fChain->SetBranchAddress("jets_AK5PF_chgEmE", &jets_AK5PF_chgEmE, &b_jets_AK5PF_chgEmE);
   fChain->SetBranchAddress("jets_AK5PF_chgHadE", &jets_AK5PF_chgHadE, &b_jets_AK5PF_chgHadE);
   fChain->SetBranchAddress("jets_AK5PF_chgMuE", &jets_AK5PF_chgMuE, &b_jets_AK5PF_chgMuE);
   fChain->SetBranchAddress("jets_AK5PF_chg_Mult", &jets_AK5PF_chg_Mult, &b_jets_AK5PF_chg_Mult);
   fChain->SetBranchAddress("jets_AK5PF_neutralEmE", &jets_AK5PF_neutralEmE, &b_jets_AK5PF_neutralEmE);
   fChain->SetBranchAddress("jets_AK5PF_neutralHadE", &jets_AK5PF_neutralHadE, &b_jets_AK5PF_neutralHadE);
   fChain->SetBranchAddress("jets_AK5PF_neutral_Mult", &jets_AK5PF_neutral_Mult, &b_jets_AK5PF_neutral_Mult);
   fChain->SetBranchAddress("jets_AK5PF_mu_Mult", &jets_AK5PF_mu_Mult, &b_jets_AK5PF_mu_Mult);
   fChain->SetBranchAddress("jets_AK5PF_emf", &jets_AK5PF_emf, &b_jets_AK5PF_emf);
   fChain->SetBranchAddress("jets_AK5PF_ehf", &jets_AK5PF_ehf, &b_jets_AK5PF_ehf);
   fChain->SetBranchAddress("jets_AK5PF_n60", &jets_AK5PF_n60, &b_jets_AK5PF_n60);
   fChain->SetBranchAddress("jets_AK5PF_n90", &jets_AK5PF_n90, &b_jets_AK5PF_n90);
   fChain->SetBranchAddress("jets_AK5PF_etaetaMoment", &jets_AK5PF_etaetaMoment, &b_jets_AK5PF_etaetaMoment);
   fChain->SetBranchAddress("jets_AK5PF_etaphiMoment", &jets_AK5PF_etaphiMoment, &b_jets_AK5PF_etaphiMoment);
   fChain->SetBranchAddress("jets_AK5PF_phiphiMoment", &jets_AK5PF_phiphiMoment, &b_jets_AK5PF_phiphiMoment);
   fChain->SetBranchAddress("jets_AK5PF_n90Hits", &jets_AK5PF_n90Hits, &b_jets_AK5PF_n90Hits);
   fChain->SetBranchAddress("jets_AK5PF_fHPD", &jets_AK5PF_fHPD, &b_jets_AK5PF_fHPD);
   fChain->SetBranchAddress("jets_AK5PF_fRBX", &jets_AK5PF_fRBX, &b_jets_AK5PF_fRBX);
   fChain->SetBranchAddress("jets_AK5PF_hitsInN90", &jets_AK5PF_hitsInN90, &b_jets_AK5PF_hitsInN90);
   fChain->SetBranchAddress("jets_AK5PF_nECALTowers", &jets_AK5PF_nECALTowers, &b_jets_AK5PF_nECALTowers);
   fChain->SetBranchAddress("jets_AK5PF_nHCALTowers", &jets_AK5PF_nHCALTowers, &b_jets_AK5PF_nHCALTowers);
   fChain->SetBranchAddress("jets_AK5PF_fSubDetector1", &jets_AK5PF_fSubDetector1, &b_jets_AK5PF_fSubDetector1);
   fChain->SetBranchAddress("jets_AK5PF_fSubDetector2", &jets_AK5PF_fSubDetector2, &b_jets_AK5PF_fSubDetector2);
   fChain->SetBranchAddress("jets_AK5PF_fSubDetector3", &jets_AK5PF_fSubDetector3, &b_jets_AK5PF_fSubDetector3);
   fChain->SetBranchAddress("jets_AK5PF_fSubDetector4", &jets_AK5PF_fSubDetector4, &b_jets_AK5PF_fSubDetector4);
   fChain->SetBranchAddress("jets_AK5PF_area", &jets_AK5PF_area, &b_jets_AK5PF_area);
   fChain->SetBranchAddress("jets_AK5PF_mass", &jets_AK5PF_mass, &b_jets_AK5PF_mass);
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
   fChain->SetBranchAddress("mus_isTrackerMuon", &mus_isTrackerMuon, &b_mus_isTrackerMuon);
   fChain->SetBranchAddress("mus_isStandAloneMuon", &mus_isStandAloneMuon, &b_mus_isStandAloneMuon);
   fChain->SetBranchAddress("mus_isCaloMuon", &mus_isCaloMuon, &b_mus_isCaloMuon);
   fChain->SetBranchAddress("mus_isGlobalMuon", &mus_isGlobalMuon, &b_mus_isGlobalMuon);
   fChain->SetBranchAddress("mus_isElectron", &mus_isElectron, &b_mus_isElectron);
   fChain->SetBranchAddress("mus_isConvertedPhoton", &mus_isConvertedPhoton, &b_mus_isConvertedPhoton);
   fChain->SetBranchAddress("mus_isPhoton", &mus_isPhoton, &b_mus_isPhoton);
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
   fChain->SetBranchAddress("Npfmets", &Npfmets, &b_Npfmets);
   fChain->SetBranchAddress("pfmets_et", &pfmets_et, &b_pfmets_et);
   fChain->SetBranchAddress("pfmets_phi", &pfmets_phi, &b_pfmets_phi);
   fChain->SetBranchAddress("pfmets_ex", &pfmets_ex, &b_pfmets_ex);
   fChain->SetBranchAddress("pfmets_ey", &pfmets_ey, &b_pfmets_ey);
   fChain->SetBranchAddress("pfmets_gen_et", &pfmets_gen_et, &b_pfmets_gen_et);
   fChain->SetBranchAddress("pfmets_gen_phi", &pfmets_gen_phi, &b_pfmets_gen_phi);
   fChain->SetBranchAddress("pfmets_sign", &pfmets_sign, &b_pfmets_sign);
   fChain->SetBranchAddress("pfmets_sumEt", &pfmets_sumEt, &b_pfmets_sumEt);
   fChain->SetBranchAddress("pfmets_unCPhi", &pfmets_unCPhi, &b_pfmets_unCPhi);
   fChain->SetBranchAddress("pfmets_unCPt", &pfmets_unCPt, &b_pfmets_unCPt);
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
   fChain->SetBranchAddress("photons_maxEnergyXtal", &photons_maxEnergyXtal, &b_photons_maxEnergyXtal);
   fChain->SetBranchAddress("photons_e1x5", &photons_e1x5, &b_photons_e1x5);
   fChain->SetBranchAddress("photons_e2x5", &photons_e2x5, &b_photons_e2x5);
   fChain->SetBranchAddress("photons_e3x3", &photons_e3x3, &b_photons_e3x3);
   fChain->SetBranchAddress("photons_e5x5", &photons_e5x5, &b_photons_e5x5);
   fChain->SetBranchAddress("photons_sigmaEtaEta", &photons_sigmaEtaEta, &b_photons_sigmaEtaEta);
   fChain->SetBranchAddress("photons_sigmaIetaIeta", &photons_sigmaIetaIeta, &b_photons_sigmaIetaIeta);
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
   fChain->SetBranchAddress("pv_chi2", &pv_chi2, &b_pv_chi2);
   fChain->SetBranchAddress("pv_ndof", &pv_ndof, &b_pv_ndof);
   fChain->SetBranchAddress("pv_tracksSize", &pv_tracksSize, &b_pv_tracksSize);
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
   fChain->SetBranchAddress("taus_taNC", &taus_taNC, &b_taus_taNC);
   fChain->SetBranchAddress("taus_byIsoUsingLeadingPi", &taus_byIsoUsingLeadingPi, &b_taus_byIsoUsingLeadingPi);
   fChain->SetBranchAddress("taus_tkIsoUsingLeadingPi", &taus_tkIsoUsingLeadingPi, &b_taus_tkIsoUsingLeadingPi);
   fChain->SetBranchAddress("taus_ecalIsoUsingLeadingPi", &taus_ecalIsoUsingLeadingPi, &b_taus_ecalIsoUsingLeadingPi);
   fChain->SetBranchAddress("taus_signalPFChargedHadrCandsSize", &taus_signalPFChargedHadrCandsSize, &b_taus_signalPFChargedHadrCandsSize);
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
   fChain->SetBranchAddress("tracks_trkExptHitsInner", &tracks_trkExptHitsInner, &b_tracks_trkExptHitsInner);
   fChain->SetBranchAddress("tracks_trkExptHitsOuter", &tracks_trkExptHitsOuter, &b_tracks_trkExptHitsOuter);
   fChain->SetBranchAddress("tracks_trks_nlayerslost", &tracks_trks_nlayerslost, &b_tracks_trks_nlayerslost);
   fChain->SetBranchAddress("tracks_trks_nlayers", &tracks_trks_nlayers, &b_tracks_trks_nlayers);
   fChain->SetBranchAddress("tracks_trksvalidpixelhits", &tracks_trksvalidpixelhits, &b_tracks_trksvalidpixelhits);
   fChain->SetBranchAddress("tracks_trkslostpixelhits", &tracks_trkslostpixelhits, &b_tracks_trkslostpixelhits);
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
   fChain->SetBranchAddress("tracks_outerPx", &tracks_outerPx, &b_tracks_outerPx);
   fChain->SetBranchAddress("tracks_outerPy", &tracks_outerPy, &b_tracks_outerPy);
   fChain->SetBranchAddress("tracks_outerPz", &tracks_outerPz, &b_tracks_outerPz);
   fChain->SetBranchAddress("tracks_algo", &tracks_algo, &b_tracks_algo);
   fChain->SetBranchAddress("tracks_highPurity", &tracks_highPurity, &b_tracks_highPurity);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
   fChain->SetBranchAddress("experimentType", &experimentType, &b_experimentType);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing);
   fChain->SetBranchAddress("orbitNumber", &orbitNumber, &b_orbitNumber);


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
   fChain->SetBranchAddress("L1Bit_0", &L1Bit_0, &b_L1Bit_0);
   fChain->SetBranchAddress("L1Bit_1", &L1Bit_1, &b_L1Bit_1);
   fChain->SetBranchAddress("L1Bit_10", &L1Bit_10, &b_L1Bit_10);
   fChain->SetBranchAddress("L1Bit_100", &L1Bit_100, &b_L1Bit_100);
   fChain->SetBranchAddress("L1Bit_101", &L1Bit_101, &b_L1Bit_101);
   fChain->SetBranchAddress("L1Bit_102", &L1Bit_102, &b_L1Bit_102);
   fChain->SetBranchAddress("L1Bit_103", &L1Bit_103, &b_L1Bit_103);
   fChain->SetBranchAddress("L1Bit_104", &L1Bit_104, &b_L1Bit_104);
   fChain->SetBranchAddress("L1Bit_105", &L1Bit_105, &b_L1Bit_105);
   fChain->SetBranchAddress("L1Bit_106", &L1Bit_106, &b_L1Bit_106);
   fChain->SetBranchAddress("L1Bit_107", &L1Bit_107, &b_L1Bit_107);
   fChain->SetBranchAddress("L1Bit_108", &L1Bit_108, &b_L1Bit_108);
   fChain->SetBranchAddress("L1Bit_109", &L1Bit_109, &b_L1Bit_109);
   fChain->SetBranchAddress("L1Bit_11", &L1Bit_11, &b_L1Bit_11);
   fChain->SetBranchAddress("L1Bit_110", &L1Bit_110, &b_L1Bit_110);
   fChain->SetBranchAddress("L1Bit_111", &L1Bit_111, &b_L1Bit_111);
   fChain->SetBranchAddress("L1Bit_112", &L1Bit_112, &b_L1Bit_112);
   fChain->SetBranchAddress("L1Bit_113", &L1Bit_113, &b_L1Bit_113);
   fChain->SetBranchAddress("L1Bit_114", &L1Bit_114, &b_L1Bit_114);
   fChain->SetBranchAddress("L1Bit_115", &L1Bit_115, &b_L1Bit_115);
   fChain->SetBranchAddress("L1Bit_116", &L1Bit_116, &b_L1Bit_116);
   fChain->SetBranchAddress("L1Bit_117", &L1Bit_117, &b_L1Bit_117);
   fChain->SetBranchAddress("L1Bit_118", &L1Bit_118, &b_L1Bit_118);
   fChain->SetBranchAddress("L1Bit_119", &L1Bit_119, &b_L1Bit_119);
   fChain->SetBranchAddress("L1Bit_12", &L1Bit_12, &b_L1Bit_12);
   fChain->SetBranchAddress("L1Bit_120", &L1Bit_120, &b_L1Bit_120);
   fChain->SetBranchAddress("L1Bit_121", &L1Bit_121, &b_L1Bit_121);
   fChain->SetBranchAddress("L1Bit_122", &L1Bit_122, &b_L1Bit_122);
   fChain->SetBranchAddress("L1Bit_123", &L1Bit_123, &b_L1Bit_123);
   fChain->SetBranchAddress("L1Bit_124", &L1Bit_124, &b_L1Bit_124);
   fChain->SetBranchAddress("L1Bit_125", &L1Bit_125, &b_L1Bit_125);
   fChain->SetBranchAddress("L1Bit_126", &L1Bit_126, &b_L1Bit_126);
   fChain->SetBranchAddress("L1Bit_127", &L1Bit_127, &b_L1Bit_127);
   fChain->SetBranchAddress("L1Bit_13", &L1Bit_13, &b_L1Bit_13);
   fChain->SetBranchAddress("L1Bit_14", &L1Bit_14, &b_L1Bit_14);
   fChain->SetBranchAddress("L1Bit_15", &L1Bit_15, &b_L1Bit_15);
   fChain->SetBranchAddress("L1Bit_16", &L1Bit_16, &b_L1Bit_16);
   fChain->SetBranchAddress("L1Bit_17", &L1Bit_17, &b_L1Bit_17);
   fChain->SetBranchAddress("L1Bit_18", &L1Bit_18, &b_L1Bit_18);
   fChain->SetBranchAddress("L1Bit_19", &L1Bit_19, &b_L1Bit_19);
   fChain->SetBranchAddress("L1Bit_2", &L1Bit_2, &b_L1Bit_2);
   fChain->SetBranchAddress("L1Bit_20", &L1Bit_20, &b_L1Bit_20);
   fChain->SetBranchAddress("L1Bit_21", &L1Bit_21, &b_L1Bit_21);
   fChain->SetBranchAddress("L1Bit_22", &L1Bit_22, &b_L1Bit_22);
   fChain->SetBranchAddress("L1Bit_23", &L1Bit_23, &b_L1Bit_23);
   fChain->SetBranchAddress("L1Bit_24", &L1Bit_24, &b_L1Bit_24);
   fChain->SetBranchAddress("L1Bit_25", &L1Bit_25, &b_L1Bit_25);
   fChain->SetBranchAddress("L1Bit_26", &L1Bit_26, &b_L1Bit_26);
   fChain->SetBranchAddress("L1Bit_27", &L1Bit_27, &b_L1Bit_27);
   fChain->SetBranchAddress("L1Bit_28", &L1Bit_28, &b_L1Bit_28);
   fChain->SetBranchAddress("L1Bit_29", &L1Bit_29, &b_L1Bit_29);
   fChain->SetBranchAddress("L1Bit_3", &L1Bit_3, &b_L1Bit_3);
   fChain->SetBranchAddress("L1Bit_30", &L1Bit_30, &b_L1Bit_30);
   fChain->SetBranchAddress("L1Bit_31", &L1Bit_31, &b_L1Bit_31);
   fChain->SetBranchAddress("L1Bit_32", &L1Bit_32, &b_L1Bit_32);
   fChain->SetBranchAddress("L1Bit_33", &L1Bit_33, &b_L1Bit_33);
   fChain->SetBranchAddress("L1Bit_34", &L1Bit_34, &b_L1Bit_34);
   fChain->SetBranchAddress("L1Bit_35", &L1Bit_35, &b_L1Bit_35);
   fChain->SetBranchAddress("L1Bit_36", &L1Bit_36, &b_L1Bit_36);
   fChain->SetBranchAddress("L1Bit_37", &L1Bit_37, &b_L1Bit_37);
   fChain->SetBranchAddress("L1Bit_38", &L1Bit_38, &b_L1Bit_38);
   fChain->SetBranchAddress("L1Bit_39", &L1Bit_39, &b_L1Bit_39);
   fChain->SetBranchAddress("L1Bit_4", &L1Bit_4, &b_L1Bit_4);
   fChain->SetBranchAddress("L1Bit_40", &L1Bit_40, &b_L1Bit_40);
   fChain->SetBranchAddress("L1Bit_41", &L1Bit_41, &b_L1Bit_41);
   fChain->SetBranchAddress("L1Bit_42", &L1Bit_42, &b_L1Bit_42);
   fChain->SetBranchAddress("L1Bit_43", &L1Bit_43, &b_L1Bit_43);
   fChain->SetBranchAddress("L1Bit_44", &L1Bit_44, &b_L1Bit_44);
   fChain->SetBranchAddress("L1Bit_45", &L1Bit_45, &b_L1Bit_45);
   fChain->SetBranchAddress("L1Bit_46", &L1Bit_46, &b_L1Bit_46);
   fChain->SetBranchAddress("L1Bit_47", &L1Bit_47, &b_L1Bit_47);
   fChain->SetBranchAddress("L1Bit_48", &L1Bit_48, &b_L1Bit_48);
   fChain->SetBranchAddress("L1Bit_49", &L1Bit_49, &b_L1Bit_49);
   fChain->SetBranchAddress("L1Bit_5", &L1Bit_5, &b_L1Bit_5);
   fChain->SetBranchAddress("L1Bit_50", &L1Bit_50, &b_L1Bit_50);
   fChain->SetBranchAddress("L1Bit_51", &L1Bit_51, &b_L1Bit_51);
   fChain->SetBranchAddress("L1Bit_52", &L1Bit_52, &b_L1Bit_52);
   fChain->SetBranchAddress("L1Bit_53", &L1Bit_53, &b_L1Bit_53);
   fChain->SetBranchAddress("L1Bit_54", &L1Bit_54, &b_L1Bit_54);
   fChain->SetBranchAddress("L1Bit_55", &L1Bit_55, &b_L1Bit_55);
   fChain->SetBranchAddress("L1Bit_56", &L1Bit_56, &b_L1Bit_56);
   fChain->SetBranchAddress("L1Bit_57", &L1Bit_57, &b_L1Bit_57);
   fChain->SetBranchAddress("L1Bit_58", &L1Bit_58, &b_L1Bit_58);
   fChain->SetBranchAddress("L1Bit_59", &L1Bit_59, &b_L1Bit_59);
   fChain->SetBranchAddress("L1Bit_6", &L1Bit_6, &b_L1Bit_6);
   fChain->SetBranchAddress("L1Bit_60", &L1Bit_60, &b_L1Bit_60);
   fChain->SetBranchAddress("L1Bit_61", &L1Bit_61, &b_L1Bit_61);
   fChain->SetBranchAddress("L1Bit_62", &L1Bit_62, &b_L1Bit_62);
   fChain->SetBranchAddress("L1Bit_63", &L1Bit_63, &b_L1Bit_63);
   fChain->SetBranchAddress("L1Bit_64", &L1Bit_64, &b_L1Bit_64);
   fChain->SetBranchAddress("L1Bit_65", &L1Bit_65, &b_L1Bit_65);
   fChain->SetBranchAddress("L1Bit_66", &L1Bit_66, &b_L1Bit_66);
   fChain->SetBranchAddress("L1Bit_67", &L1Bit_67, &b_L1Bit_67);
   fChain->SetBranchAddress("L1Bit_68", &L1Bit_68, &b_L1Bit_68);
   fChain->SetBranchAddress("L1Bit_69", &L1Bit_69, &b_L1Bit_69);
   fChain->SetBranchAddress("L1Bit_7", &L1Bit_7, &b_L1Bit_7);
   fChain->SetBranchAddress("L1Bit_70", &L1Bit_70, &b_L1Bit_70);
   fChain->SetBranchAddress("L1Bit_71", &L1Bit_71, &b_L1Bit_71);
   fChain->SetBranchAddress("L1Bit_72", &L1Bit_72, &b_L1Bit_72);
   fChain->SetBranchAddress("L1Bit_73", &L1Bit_73, &b_L1Bit_73);
   fChain->SetBranchAddress("L1Bit_74", &L1Bit_74, &b_L1Bit_74);
   fChain->SetBranchAddress("L1Bit_75", &L1Bit_75, &b_L1Bit_75);
   fChain->SetBranchAddress("L1Bit_76", &L1Bit_76, &b_L1Bit_76);
   fChain->SetBranchAddress("L1Bit_77", &L1Bit_77, &b_L1Bit_77);
   fChain->SetBranchAddress("L1Bit_78", &L1Bit_78, &b_L1Bit_78);
   fChain->SetBranchAddress("L1Bit_79", &L1Bit_79, &b_L1Bit_79);
   fChain->SetBranchAddress("L1Bit_8", &L1Bit_8, &b_L1Bit_8);
   fChain->SetBranchAddress("L1Bit_80", &L1Bit_80, &b_L1Bit_80);
   fChain->SetBranchAddress("L1Bit_81", &L1Bit_81, &b_L1Bit_81);
   fChain->SetBranchAddress("L1Bit_82", &L1Bit_82, &b_L1Bit_82);
   fChain->SetBranchAddress("L1Bit_83", &L1Bit_83, &b_L1Bit_83);
   fChain->SetBranchAddress("L1Bit_84", &L1Bit_84, &b_L1Bit_84);
   fChain->SetBranchAddress("L1Bit_85", &L1Bit_85, &b_L1Bit_85);
   fChain->SetBranchAddress("L1Bit_86", &L1Bit_86, &b_L1Bit_86);
   fChain->SetBranchAddress("L1Bit_87", &L1Bit_87, &b_L1Bit_87);
   fChain->SetBranchAddress("L1Bit_88", &L1Bit_88, &b_L1Bit_88);
   fChain->SetBranchAddress("L1Bit_89", &L1Bit_89, &b_L1Bit_89);
   fChain->SetBranchAddress("L1Bit_9", &L1Bit_9, &b_L1Bit_9);
   fChain->SetBranchAddress("L1Bit_90", &L1Bit_90, &b_L1Bit_90);
   fChain->SetBranchAddress("L1Bit_91", &L1Bit_91, &b_L1Bit_91);
   fChain->SetBranchAddress("L1Bit_92", &L1Bit_92, &b_L1Bit_92);
   fChain->SetBranchAddress("L1Bit_93", &L1Bit_93, &b_L1Bit_93);
   fChain->SetBranchAddress("L1Bit_94", &L1Bit_94, &b_L1Bit_94);
   fChain->SetBranchAddress("L1Bit_95", &L1Bit_95, &b_L1Bit_95);
   fChain->SetBranchAddress("L1Bit_96", &L1Bit_96, &b_L1Bit_96);
   fChain->SetBranchAddress("L1Bit_97", &L1Bit_97, &b_L1Bit_97);
   fChain->SetBranchAddress("L1Bit_98", &L1Bit_98, &b_L1Bit_98);
   fChain->SetBranchAddress("L1Bit_99", &L1Bit_99, &b_L1Bit_99);

}

