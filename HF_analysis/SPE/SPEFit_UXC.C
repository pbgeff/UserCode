#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"

#include <vector>
#include <utility>
#include <algorithm>
#include <map>
#include <fstream>
#include <sstream>

Int_t Nev;


Double_t FitFun(Double_t *x, Double_t *par);

void SPEFit(char * fLEDname, char * fPEDname, int run, int LED_amp, double cutmax = 250.0)
{

  //set plotting styles
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetStatColor(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

    //set file names
    stringstream out_fname;
    stringstream out_fname1;
    out_fname<<"SPEconstants_Run_"<<run<<".txt";
    out_fname1<<"SPEspec_Run_"<<run<<".txt";

    ofstream  constants_file(out_fname.str().c_str(),ios_base::trunc); 
    //ofstream  constants_file1(out_fname1.str().c_str(),ios_base::trunc); 
    constants_file<<"Run "<<run<<endl;
    constants_file<<"type SPE"<<endl;
    constants_file<<"LED_amplitude "<<LED_amp<<endl<<endl;

    constants_file<<endl<<"LED_amplitude Depth Phi Eta Ped_mean Ped_mean_err Ped_RMS  Ped_RMS_err SPEPeak_RMS SPEPeak_RMS_err Gain Gain_err Normalized_Chi2 MeanPE_fit MeanPE_fit_err MeanPE_estimate PE5flag"<<endl;

    out_fname.str("");
    out_fname<<"SPEdistributions_Run_"<<run<<".txt";

    out_fname.str("");
    out_fname<<"SPEextra_Run_"<<run<<".txt";
    //ofstream  extra_file(out_fname.str().c_str(),ios_base::trunc); 


    double scale = 1.0;
    scale = 2.6; //Need to scale up HF charge
    double fC2electrons = 6240.; //convert fC to #electrons

    char spename[128], pedname[128], spehistname[128];
 
    TFile *tfLED = new TFile(fLEDname);
    TFile *tfPED = new TFile(fPEDname);
    


    //const int NnewBins = 106;
    //double binsX[NnewBins] = {0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122,124,126,128,130,132,134,136,138,140,142,144,146,148,150,152,154,156,158,160,162,164,166,168,170,180,190,200,210,220,230,240,250,266,282,298,316,336,356,378,404,430,456,482,500};

    const int NnewBins = 80;
    double binsX[NnewBins] = {0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,96,99,102,105,108,111,114,117,120,123,126,129,132,135,138,141,144,147,150,153,156,159,162,165,168,171,174,177,180,190,200,210,220,230,240,250,266,282,298,316,336,356,378,404,430,456,482,500};	  
    TH1F* hspe = new TH1F("hspe","hspe",NnewBins-1,binsX);


    int NDepth = 2; //number of depths
    int MinDepth = 1;
    int MaxDepth = 2;
    int MinEta = 29; 
    int MaxEta = 41;
    int MinPhi = 41;
    int MaxPhi = 53;
   
 
    TCanvas *Carray[NDepth+1][MaxPhi+1];
    bool drawflag[NDepth+1][MaxPhi+1];  
    TH1F *LED[NDepth+1][MaxEta+1][MaxPhi+1];
    TH1F *PED[NDepth+1][MaxEta+1][MaxPhi+1];

    for(int iDepth = MinDepth; iDepth <= MaxDepth; iDepth++){
      for(int iPhi = MinPhi; iPhi <= MaxPhi; iPhi++){

	bool nonNull = false;

	for(int iEta = MinEta; iEta <= MaxEta; iEta++){

	  sprintf(spename,"Analyzer/CommonDir/ResPlotDir/Histo_for_Depth_%d_Eta_%d_Phi_%d",iDepth,iEta,iPhi);
	  LED[iDepth][iEta][iPhi]=(TH1F *)tfLED->Get(spename);
	  if(LED[iDepth][iEta][iPhi]) nonNull = true;
      
	  sprintf(spename,"Analyzer/CommonDir/ResPlotDir/Histo_for_Depth_%d_Eta_%d_Phi_%d",iDepth,iEta,iPhi);
	  PED[iDepth][iEta][iPhi]=(TH1F *)tfPED->Get(spename);
	}

	drawflag[iDepth][iPhi] = false;
	char canvname[16];
	sprintf(canvname, "c_%d_%d", iDepth,iPhi);
	if(nonNull){ //only create canvas if distributions exist
	  Carray[iDepth][iPhi] = new TCanvas(canvname,canvname,1200,700);
	  Carray[iDepth][iPhi]->Divide(5,3);
	}

      }
    }



    int HV=0;

    for(int iDepth = MinDepth; iDepth <= MaxDepth; iDepth++){
      for(int iPhi = MinPhi; iPhi <= MaxPhi; iPhi++){
	for(int iEta = MinEta; iEta <= MaxEta; iEta++){

	  //cout<<iDepth<<" "<<iPhi<<" "<<iEta<<endl;

	  if(!LED[iDepth][iEta][iPhi]) continue;

	  sprintf(spehistname,"led %d %d %d",iDepth,iEta,iPhi);
	  TH1F *hspe_temp = (TH1F *)LED[iDepth][iEta][iPhi]->Clone(spehistname);
	  sprintf(spehistname,"ped %d %d %d",iDepth,iEta,iPhi);
	  TH1F *hped = (TH1F *)PED[iDepth][iEta][iPhi]->Clone(spehistname);
	  hspe->Reset();
	  sprintf (spehistname, "SumLED_Depth_%d_Eta_%d_Phi_%d",iDepth,iEta,iPhi);
	  hspe->SetTitle(spehistname);

	  //combine bins of original SPE histogram
	  for(int ib=1; ib<=hspe_temp->GetNbinsX(); ib++) {
	    double bin_center = hspe_temp->GetBinCenter(ib);
	    if(bin_center>hspe->GetXaxis()->GetXmax()) continue;
	    int newbin = hspe->FindBin(bin_center);
	    double new_content = hspe->GetBinContent(newbin) + hspe_temp->GetBinContent(ib);
	    double new_error = sqrt(pow(hspe->GetBinError(newbin),2)+pow(hspe_temp->GetBinError(ib),2));
	    hspe->SetBinContent(newbin,new_content);
	    hspe->SetBinError(newbin,new_error);
	  }
	  TH1F* hspe_unscaled = (TH1F*)hspe->Clone("hspe_unscaled");
	  //renormalize bins of new SPE histogram
	  for(int ib=1; ib<=hspe->GetNbinsX(); ib++) {
	    double new_content = hspe->GetBinContent(ib)/hspe->GetXaxis()->GetBinWidth(ib)*hspe_temp->GetXaxis()->GetBinWidth(1);
	    double new_error = hspe->GetBinError(ib)/hspe->GetXaxis()->GetBinWidth(ib)*hspe_temp->GetXaxis()->GetBinWidth(1);
	    hspe->SetBinContent(ib,new_content);
	    hspe->SetBinError(ib,new_error);
	  }
	  
	  if(hspe_temp->Integral()==0) continue;
	  else drawflag[iDepth][iPhi] = true;	  

	  Nev = hspe_temp->Integral()*hspe_temp->GetXaxis()->GetBinWidth(1); 
	  TF1 *fped = new TF1("fped","gaus",0, 80);
	  hped->Fit(fped,"NQR");
	  double pploc = fped->GetParameter(1), ppwidth = fped->GetParameter(2);
	  hspe->Fit(fped, "NQ", "", pploc - 3*ppwidth, pploc + ppwidth);  
	  
	  //estimate SPE peak location
	  int max_SPE_bin, maxbin, Nbins;
	  double max_SPE_height=0, minheight, max_SPE_location;
	  bool minflag = false;
	  maxbin=hspe->FindBin(fped->GetParameter(1)); //location of pedestal peak
	  minheight=hspe->GetBinContent(maxbin); //initialize minheight
	  Nbins = hspe->GetNbinsX();
	  for(int j=maxbin+1; j<Nbins-1; j++) { //start from pedestal peak and loop through bins
	    if(hspe->GetBinContent(j) > minheight && !minflag) minflag=true; //only look for SPE peak when minflag=true
	    if(hspe->GetBinContent(j) < minheight )  minheight = hspe->GetBinContent(j);
	    if(minflag && hspe->GetBinContent(j) > max_SPE_height){
	      max_SPE_bin = j;
	      max_SPE_location = hspe->GetBinCenter(max_SPE_bin);
	      max_SPE_height = hspe->GetBinContent(j);
	    }
	  } //start from pedestal peak and loop through bins
	  //find minimum bin between pedestal and SPE peaks
	  hspe->GetXaxis()->SetRange(maxbin,max_SPE_bin);
	  int minbin = hspe->GetMinimumBin(); 
	  double minbin_location = hspe->GetBinCenter(minbin);
	  hspe->GetXaxis()->SetRange(1,Nbins);	    
	  
	  TF1 *fit = new TF1("fit", FitFun, 0, 500, 5);
	    
	  double mu = - log(fped->Integral(0,100)/Nev);
	  if(mu<0) mu=0.01;
	  double gain_est = max_SPE_location-1.0*fped->GetParameter(1);
	  if(max_SPE_bin > (minbin+1)) fit->SetParameters(mu, 20, 1, gain_est, gain_est*0.5);
	  else fit->SetParameters(mu, 20, 1, 2.1*fped->GetParameter(2), 10); //case of no clear minimum; start looking for SPE peak at 2sigma away from pedestal peak
	  fit->SetParLimits(0, 0, 10);
	  fit->FixParameter(1, fped->GetParameter(1));
	  fit->FixParameter(2, fped->GetParameter(2));
	  fit->SetParLimits(3, fped->GetParameter(2)*2, 350);
	  fit->SetParLimits(4, fped->GetParameter(2)*1.01, 250);
	  

	  double maxfitrange = 500.;    
	  double minfitrange = 0.;
	  hspe->Fit(fit, "MNQL", "", minfitrange, maxfitrange);
	  maxfitrange = fped->GetParameter(1)+4*fit->GetParameter(3)+fit->GetParameter(4);
	  if(500<maxfitrange) maxfitrange = 500;
	  hspe->Fit(fit, "MNQL", "", minfitrange, maxfitrange);

	  //calculate NDOF of fit excluding bins with 0 entries
	  int myNDOF=-3; //three free parameters
	  for(int j=hspe->FindBin(minfitrange); j<=hspe->FindBin(maxfitrange); j++) { //loop through fitted spe bins
	    if(hspe->GetBinContent(j)) myNDOF++;
	  } //loop through fitted spe bins


	  //calculate means and integrals of the fit and data
	  double fint, fint_error, hint, favg, havg;
	  int temp_lowbin, temp_highbin;
	  temp_lowbin = hspe->FindBin(minfitrange);
	  temp_highbin = hspe->FindBin(maxfitrange);
	  hspe_unscaled->GetXaxis()->SetRangeUser(minfitrange, maxfitrange);
	  havg = hspe_unscaled->GetMean();
	  hint = hspe->Integral(temp_lowbin,temp_highbin,"width");
	  double min_frange = hspe->GetBinLowEdge(temp_lowbin);
	  favg = fit->Mean(min_frange, maxfitrange);
	  fint = fit->Integral(min_frange, maxfitrange);
	  //fint_error = fit->IntegralError(min_frange, maxfitrange);
	  
	  double PE5int = 0; //integral of events with >=5 PE
	  double PE5loc =  fped->GetParameter(1)+ 5*fit->GetParameter(3);
	  if(PE5loc>500) PE5int = 0;
	  else {
	    int PE5bin =  hspe_temp->FindBin(PE5loc);
	    temp_highbin = hspe_temp->FindBin(maxfitrange)-1;
	    PE5int =  hspe_temp->Integral(PE5bin,temp_highbin,"width");
	  }
	  int PE5flag = 0;
	  if(PE5int/hint>0.05) PE5flag = 1; //set flag if more than 5% of events in the fit correspond to >=5PE
	//=========================================    
	  //for(int i1=1;i1<hspe->GetNbinsX();i1++){
	    //constants_file1<<HV<<"\t"<<iDepth<<"\t"<<iEta<<"\t"<<iPhi<<"\t"<<2.6*hspe->GetBinCenter(i1)<<"\t"<<hspe->GetBinContent(i1)<<"\t"<<fit->Eval(hspe->GetBinCenter(i1))<<"\n";
          //}
        //=========================================    

	  //printf("%d\n",myNDOF);
	  //output calibrations constants
	  //constants_file<<endl<<"LED_amplitude HV Spigot Channel Ped_mean Ped_mean_err Ped_RMS  Ped_RMS_err SPEPeak_RMS SPEPeak_RMS_err Gain Gain_err Normalized_Chi2 MeanPE_fit MeanPE_fit_err MeanPE_estimate PE5flag"<<endl;
	  constants_file<<LED_amp<<" "<<iDepth<<" "<<iPhi<<" "<<iEta<<" "<<scale*fped->GetParameter(1)<<" "<<scale*fped->GetParError(1)<<" "<<scale*fped->GetParameter(2)<<" "<<scale*fped->GetParError(2)<<" "<<scale*fit->GetParameter(4)<<" "<<scale*fit->GetParError(4)<<" "<<scale*fit->GetParameter(3)*fC2electrons<<" "<<scale*fit->GetParError(3)*fC2electrons<<" "<<fit->GetChisquare()/myNDOF/*fit->GetNDF()*/<<" "<<fit->GetParameter(0)<<" "<<fit->GetParError(0)<<" "<<mu<<" "<<PE5flag<<endl;
	    

	  /*
	  if(iDepth==2 && iPhi==53 && iEta==36){
	    cout<<iDepth<<" "<<iPhi<<" "<<iEta<<" "<<gain_est<<" "<<fit->GetParameter(3)<<endl;
	    cout<<LED_amp<<" "<<iDepth<<" "<<iPhi<<" "<<iEta<<" "<<scale*fped->GetParameter(1)<<" "<<scale*fped->GetParError(1)<<" "<<scale*fped->GetParameter(2)<<" "<<scale*fped->GetParError(2)<<" "<<scale*fit->GetParameter(4)<<" "<<scale*fit->GetParError(4)<<" "<<scale*fit->GetParameter(3)*fC2electrons<<" "<<scale*fit->GetParError(3)*fC2electrons<<" "<<fit->GetChisquare()/myNDOF<<" "<<fit->GetParameter(0)<<" "<<fit->GetParError(0)<<" "<<mu<<" "<<PE5flag<<endl;
	  }
	  */

	  Carray[iDepth][iPhi]->cd(iEta-MinEta+1);
	  gPad->SetBorderMode(0);
	  gPad->SetBorderSize(0);
	  gPad->SetRightMargin(0.01);
	  gPad->SetBottomMargin(0.1);
	  gPad->SetLogy(true);
	  hspe->GetXaxis()->SetRangeUser(0, 200 /*300*//*508*/);
	  hspe->SetLineColor(kBlue);
	  hspe->DrawClone("hist");
	  fit->SetLineWidth(2);
	  fit->Draw("same");

	}
    
	if(drawflag[iDepth][iPhi]) { //draw plots of fit if data for the HV is present
	  stringstream plot_name;
	  plot_name<<"Plots/SPEFits_Run_"<<run<<"_Depth"<<iDepth<<"_Phi"<<iPhi<<".pdf";
	  Carray[iDepth][iPhi]->SaveAs(plot_name.str().c_str());
	  plot_name.str( std::string() );
	}

      }
    }

    constants_file.close();
    //constants_file1.close();
}

Double_t FitFun(Double_t *x, Double_t *par){
        
  Double_t sum,xx,A0,C0,r0,sigma0,mean1,sigma1,A1,C1,r1,mean2,sigma2,A2,C2,r2,mean3,sigma3,A3,C3,r3,mean4,sigma4,A4,C4,r4;
        
  const Double_t k0=2.0, k1=1.6, k2=2.0, k3=2.0, k4=2.0;
        
  xx=x[0];  
  sigma0 = par[2];
  A0 = Nev/TMath::Exp(par[0]);
  r0 = ((xx-par[1])/sigma0);
  C0 = 1/(sigma0* TMath::Exp(-k0*k0/2)/k0 +
          sigma0*sqrt(2*3.14159)*0.5*(1+TMath::Erf(k0/1.41421)));
  //sum = 1/(sqrt(2*3.14159)*par[2])*A0*TMath::Exp(-0.5*r0*r0);
  if(r0 < k0) sum = C0*A0*TMath::Exp(-0.5*r0*r0);
  else sum = C0*A0*TMath::Exp(0.5*k0*k0-k0*r0);
        
  mean1 = par[1]+par[3];
  sigma1 = par[4];
  A1 = A0*par[0];
  C1 = 1/(sigma1* TMath::Exp(-k1*k1/2)/k1 +
          sigma1*sqrt(2*3.14159)*0.5*(1+TMath::Erf(k1/1.41421)));
  r1 = ((xx-mean1)/sigma1);
  if(r1 < k1 && xx>par[1]-3.0*par[2]) sum += C1*A1*TMath::Exp(-0.5*r1*r1);
  //else if(r1 < k1 && xx<par[1]) sum += C1*A1*TMath::Exp(-0.5*r1*r1);
  else if(r1 < k1 && xx<par[1]) sum += 0;
  else sum += C1*A1*TMath::Exp(0.5*k1*k1-k1*r1);

                
  mean2 = 2*par[3]+par[1];
  sigma2 = sqrt(2*sigma1*sigma1 - pow(par[2],2));
  A2 = A0*par[0]*par[0]/2;
  C2 = 1/(sigma2* TMath::Exp(-k2*k2/2)/k2 +
          sigma2*sqrt(2*3.14159)*0.5*(1+TMath::Erf(k2/1.41421)));
  r2 = ((xx-mean2)/sigma2);
  if(r2 < k2 && xx>par[1]-3.0*par[2]) sum += C2*A2*TMath::Exp(-0.5*r2*r2);
  //else if(r2 < k2 && xx<par[1]) sum += C1*A1*TMath::Exp(-0.5*r1*r1);
  else if(r2 < k2 && xx<par[1]) sum += 0;
  else sum += C2*A2*TMath::Exp(0.5*k2*k2-k2*r2);

 
  mean3 = 3*par[3]+par[1];
  sigma3 = sqrt(3*sigma1*sigma1 - 2*pow(par[2],2));
  A3 = A0*par[0]*par[0]*par[0]/6;
  //C3 = 1/(sigma3*sqrt(2*3.14159));
  C3 = 1/(sigma3* TMath::Exp(-k3*k3/2)/k3 + 
	  sigma3*sqrt(2*3.14159)*0.5*(1+TMath::Erf(k3/1.41421)));
    r3 = ((xx-mean3)/sigma3);
  if(r3 < k3 && xx>par[1]-3.0*par[2]) sum += C3*A3*TMath::Exp(-0.5*r3*r3);
  else if(r3 < k3 && xx<par[1]) sum += 0;
  else sum += C3*A3*TMath::Exp(0.5*k3*k3-k3*r3); 
  

  mean4 = 4*par[3]+par[1];
  sigma4 = sqrt(4*sigma1*sigma1 - 3*pow(par[2],2));
  A4 = A0*par[0]*par[0]*par[0]*par[0]/24;
  //C4 = 1/(sigma4*sqrt(2*3.14159));
  C4 = 1/(sigma4* TMath::Exp(-k4*k4/2)/k4 + 
	  sigma4*sqrt(2*3.14159)*0.5*(1+TMath::Erf(k4/1.41421)));
  r4 = ((xx-mean4)/sigma4);
  if(r4 < k4 && xx>par[1]-3.0*par[2]) sum += C4*A4*TMath::Exp(-0.5*r4*r4);
  else if(r4 < k4 && xx<par[1]) sum += 0;
  else sum += C4*A4*TMath::Exp(0.5*k4*k4-k4*r4);
  

  return sum;
}

