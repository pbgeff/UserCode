
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooGlobalFunc.h"
#include "RooFitResult.h"
#include "RooRandom.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/MCMCInterval.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/ProposalHelper.h"
#include "RooStats/SimpleInterval.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/PointSetInterval.h"
#include "RooRealVar.h"
#include <sstream>
#include <algorithm>

using namespace RooFit;
using namespace RooStats;


// declarations
void SetConstants(RooWorkspace * w, RooStats::ModelConfig * mc);
void SetConstant(const RooArgSet * vars, Bool_t value );


void new_RA4(){
  
  // let's time this challenging example
  TStopwatch t;
  t.Start();

  // set RooFit random seed for reproducible results
  RooRandom::randomGenerator()->SetSeed(4357);

  // make model
  RooWorkspace* wspace = new RooWorkspace("wspace");

  wspace->factory("Gaussian::sigCons(prime_SigEff[0,-5,5], nom_SigEff[0,-5,5], 1)");
  wspace->factory("expr::SigEff('1.0*pow(1.20,@0)',prime_SigEff)"); // // 1+-20%, 1.20=exp(20%)

  wspace->factory("Poisson::on(non[0,50], sum::splusb(prod::SigUnc(s[0,0,50],SigEff),mainb[8.8,0,50],dilep[0.9,0,20],tau[2.3,0,20],QCD[0.,0,10],MC[0.1,0,4]))");

  wspace->factory("Gaussian::mcCons(prime_rho[0,-5,5], nom_rho[0,-5,5], 1)");
  wspace->factory("expr::rho('1.0*pow(1.39,@0)',prime_rho)"); // // 1+-39%
  wspace->factory("Poisson::off(noff[0,200], prod::rhob(mainb,rho,mu_plus_e[0.74,0.01,10],1.08))");
  wspace->factory("Gaussian::mcCons2(mu_plus_enom[0.74,0.01,4], mu_plus_e, sigmatwo[.05])");

  wspace->factory("Gaussian::dilep_pred(dilep_nom[0.9,0,20], dilep, sigma3[2.2])");
  wspace->factory("Gaussian::tau_pred(tau_nom[2.3,0,20], tau, sigma4[0.5])");
  wspace->factory("Gaussian::QCD_pred(QCD_nom[0.0,0,10], QCD, sigma5[1.0])");
  wspace->factory("Gaussian::MC_pred(MC_nom[0.1,0.01,4], MC, sigma7[0.14])");

  wspace->factory("PROD::model(on,off,mcCons,mcCons2,sigCons,dilep_pred,tau_pred,QCD_pred,MC_pred)");

  RooArgSet obs(*wspace->var("non"), *wspace->var("noff"), *wspace->var("mu_plus_enom"), *wspace->var("dilep_nom"), *wspace->var("tau_nom"), "obs");
  obs.add(*wspace->var("QCD_nom"));  obs.add(*wspace->var("MC_nom"));
  RooArgSet globalObs(*wspace->var("nom_SigEff"), *wspace->var("nom_rho"), "global_obs");
  // fix global observables to their nominal values
  wspace->var("nom_SigEff")->setConstant();
  wspace->var("nom_rho")->setConstant();

  RooArgSet poi(*wspace->var("s"), "poi");
  RooArgSet nuis(*wspace->var("mainb"), *wspace->var("prime_rho"), *wspace->var("prime_SigEff"), *wspace->var("mu_plus_e"), *wspace->var("dilep"), *wspace->var("tau"), "nuis");
  nuis.add(*wspace->var("QCD"));  nuis.add(*wspace->var("MC"));


  wspace->factory("Uniform::prior_poi({s})");
  wspace->factory("Uniform::prior_nuis({mainb,mu_plus_e,dilep,tau,QCD,MC})");
  wspace->factory("PROD::prior(prior_poi,prior_nuis)");

  wspace->var("non")->setVal(8); //observed
  //wspace->var("non")->setVal(12); //expected observation
  wspace->var("noff")->setVal(7); //observed events in control region
  wspace->var("mu_plus_enom")->setVal(0.74);
  wspace->var("dilep_nom")->setVal(0.9);
  wspace->var("tau_nom")->setVal(2.3);
  wspace->var("QCD")->setVal(0.0);
  wspace->var("MC")->setVal(0.1);


  RooDataSet * data = new RooDataSet("data","",obs);
  data->add(obs);
  wspace->import(*data);


  /////////////////////////////////////////////////////
  // Now the statistical tests
  // model config
  ModelConfig* pSbModel = new ModelConfig("SbModel");
  pSbModel->SetWorkspace(*wspace);
  pSbModel->SetPdf(*wspace->pdf("model"));
  pSbModel->SetPriorPdf(*wspace->pdf("prior"));
  pSbModel->SetParametersOfInterest(poi);
  pSbModel->SetNuisanceParameters(nuis);
  pSbModel->SetObservables(obs);
  pSbModel->SetGlobalObservables(globalObs);
  wspace->import(*pSbModel);

  // set all but obs, poi and nuisance to const
  SetConstants(wspace, pSbModel);
  wspace->import(*pSbModel);


  Double_t poiValueForBModel = 0.0;
  ModelConfig* pBModel = new ModelConfig(*(RooStats::ModelConfig *)wspace->obj("SbModel"));
  pBModel->SetName("BModel");
  pBModel->SetWorkspace(*wspace);
  wspace->import(*pBModel);


  RooAbsReal * pNll = pSbModel->GetPdf()->createNLL(*data);
  RooAbsReal * pProfile = pNll->createProfile(RooArgSet());
  pProfile->getVal(); // this will do fit and set POI and nuisance parameters to fitted values
  RooArgSet * pPoiAndNuisance = new RooArgSet();
  //if(pSbModel->GetNuisanceParameters())
  //  pPoiAndNuisance->add(*pSbModel->GetNuisanceParameters());
  pPoiAndNuisance->add(*pSbModel->GetParametersOfInterest());
  cout << "\nWill save these parameter points that correspond to the fit to data" << endl;
  pPoiAndNuisance->Print("v");
  pSbModel->SetSnapshot(*pPoiAndNuisance);
  delete pProfile;
  delete pNll;
  delete pPoiAndNuisance;


  pNll = pBModel->GetPdf()->createNLL(*data);
  pProfile = pNll->createProfile(poi);
  ((RooRealVar *)poi.first())->setVal(poiValueForBModel);
  pProfile->getVal(); // this will do fit and set nuisance parameters to profiled values
  pPoiAndNuisance = new RooArgSet();
  //if(pBModel->GetNuisanceParameters())
  //  pPoiAndNuisance->add(*pBModel->GetNuisanceParameters());
  pPoiAndNuisance->add(*pBModel->GetParametersOfInterest());
  cout << "\nShould use these parameter points to generate pseudo data for bkg only" << endl;
  pPoiAndNuisance->Print("v");
  pBModel->SetSnapshot(*pPoiAndNuisance);
  delete pProfile;
  delete pNll;
  delete pPoiAndNuisance;


  // inspect workspace
  wspace->Print();

  // save workspace to file
  wspace->writeToFile("tight.root");
  //wspace->writeToFile("tight_median.root");


  // clean up
  delete wspace;
  delete data;
  delete pSbModel;
  delete pBModel;

}


////////////////////////////////////////////////////////////

// helper functions

void SetConstants(RooWorkspace * pWs, RooStats::ModelConfig * pMc){
  //
  // Fix all variables in the PDF except observables, POI and
  // nuisance parameters. Note that global observables are fixed.
  // If you need global observables floated, you have to set them
  // to float separately.
  //
  
  pMc->SetWorkspace(*pWs);
  
  RooAbsPdf * pPdf = pMc->GetPdf(); // we do not own this
  
  RooArgSet * pVars = pPdf->getVariables(); // we do own this
  
  RooArgSet * pFloated = new RooArgSet(*pMc->GetObservables());
  pFloated->add(*pMc->GetParametersOfInterest());
  pFloated->add(*pMc->GetNuisanceParameters());
  
  TIterator * pIter = pVars->createIterator(); // we do own this
  
  for(TObject * pObj = pIter->Next(); pObj; pObj = pIter->Next() ){
    std::string _name = pObj->GetName();
    RooRealVar * pFloatedObj = (RooRealVar *)pFloated->find(_name.c_str());
    if (pFloatedObj){
      ((RooRealVar *)pObj)->setConstant(kFALSE);
    }
    else{
      ((RooRealVar *)pObj)->setConstant(kTRUE);
    }
    //pObj->Print();
  }
  
  delete pIter;
  delete pVars;
  delete pFloated;
  
  return;
}


void SetConstant(const RooArgSet * vars, Bool_t value ){
  //
  // Set the constant attribute for all vars in the set
  //

  TIterator * pIter = vars->createIterator(); // we do own this

  for(TObject * pObj = pIter->Next(); pObj; pObj = pIter->Next() ){
    ((RooRealVar *)pObj)->setConstant(value);
    //pObj->Print();
  }

  delete pIter;

  return;
}


