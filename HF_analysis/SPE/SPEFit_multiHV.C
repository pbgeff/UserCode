#include "TFile.h"
#include "TH1.h"
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

void SPEFit(char * fname, int run, int LED_amp, double cutmax = 250.0)
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
    out_fname<<"SPEconstants_Run_"<<run<<".txt";

    ofstream  constants_file(out_fname.str().c_str(),ios_base::trunc); 
    constants_file<<"Run "<<run<<endl;
    constants_file<<"type SPE"<<endl;
    constants_file<<"LED_amplitude "<<LED_amp<<endl<<endl;

    constants_file<<endl<<"LED_amplitude HV Spigot Channel Ped_mean Ped_mean_err Ped_RMS  Ped_RMS_err SPEPeak_RMS SPEPeak_RMS_err Gain Gain_err Normalized_Chi2 MeanPE_fit MeanPE_fit_err MeanPE_estimate PE5flag"<<endl;

    out_fname.str("");
    out_fname<<"SPEdistributions_Run_"<<run<<".txt";

    out_fname.str("");
    out_fname<<"SPEextra_Run_"<<run<<".txt";
    ofstream  extra_file(out_fname.str().c_str(),ios_base::trunc); 

    //extra_file<<endl<<"LED_amplitude HV Spigot Channel PedSubtracted_mean Gain Gain_err MeanPE_fit MeanPE_fit_err MeanPE_estimate"<<endl;


    double scale = 1.0;
    scale = 2.6; //Need to scale up HF charge
    double fC2electrons = 6240.; //convert fC to #electrons

    int QIECh[]={18,24,13,19,20,14,22,16,17,23,15,21,6,12,1,7,8,2,10,4,5,11,3,9};
    int cde[]={1,2,2,1,2,1,2,1,1,2,1,2,1,2,1,2,2,1,2,1,1,2,1,2};
    int ceta[]={39,39,40,40,37,37,38,38,36,36,35,35,34,34,33,33,31,31,32,32,30,30,29,29};

    char spename[128], pedname[128], spehistname[128];
    bool drawflag;   
 
    TFile *tf = new TFile(fname);
    
    TCanvas *c1 = new TCanvas("c1","c1",1200,700);
    c1->Divide(6,4);
    c1->SetBorderMode(0);
    c1->SetBorderSize(0);
    TCanvas *c2 = new TCanvas("c2","c2",1200,700);
    c2->Divide(6,4);
    c2->SetBorderMode(0);
    c2->SetBorderSize(0);  

    //TCanvas *c3 = new TCanvas("c3","c3",500,500);

    const int NnewBins = 106;
    double binsX[NnewBins] = {0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122,124,126,128,130,132,134,136,138,140,142,144,146,148,150,152,154,156,158,160,162,164,166,168,170,180,190,200,210,220,230,240,250,266,282,298,316,336,356,378,404,430,456,482,500};	  
    TH1F* hspe = new TH1F("hspe","hspe",NnewBins-1,binsX);


    for(int HV=600; HV<=750; HV+=50) {
      drawflag=false;

      for (int iSpig = 0; iSpig < 2; iSpig++) {
	for(int i = 0; i < 24; i++) {

	    //if(iSpig!=0 || QIECh[i]!=24 || HV!=700) continue; //for debugging
	    //cout<<"Spig "<<iSpig<<" Ch "<<i<<endl;

	    sprintf (spename, "Analyzer/QIEsSumLED%d_QiECh_%d_DBOX_%d_eta_%d_D_%d", HV, QIECh[i], iSpig, ceta[i], cde[i]);
	    TH1F *hspe_temp = (TH1F *)tf->Get(spename);
	    sprintf (pedname, "Analyzer/QIEsSumPED%d_QiECh_%d_DBOX_%d_eta_%d_D_%d", HV, QIECh[i], iSpig, ceta[i], cde[i]);
	    TH1F *hped = (TH1F *)tf->Get(pedname);


	    hspe->Reset();
	    sprintf (spehistname, "SumLED%d_QiECh_%d_DBOX_%d_eta_%d_D_%d", HV, QIECh[i], iSpig, ceta[i], cde[i]);
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
	    else drawflag=true;	  
    
            Nev = hspe->Integral(1,hspe->GetNbinsX(),"width");
	    //cout<<"Nev "<<Nev<<endl;
	    //cout<<"hspe_temp Integral "<<hspe_temp->Integral(1,hspe_temp->FindBin(499),"width")<<endl;

	    
	    TF1 *fped = new TF1("fped","gaus",0, 80);
	    hped->Fit(fped,"NQR");
	    double pploc = fped->GetParameter(1), ppwidth = fped->GetParameter(2);
	    //cout<<"Ped only: ped mean "<<fped->GetParameter(1)<<", ped width "<<fped->GetParameter(2)<<" normalization "<<fped->GetParameter(0)<<endl;
	    hspe->Fit(fped, "NQ", "", pploc - 3*ppwidth, pploc + ppwidth);
	    //cout<<"SPE distribution: ped mean "<<fped->GetParameter(1)<<", ped width "<<fped->GetParameter(2)<<" normalization "<<fped->GetParameter(0)<<endl;
	    
	    
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
	    hspe->Fit(fit, "MNQL", "", 0, maxfitrange);
	    //double maxfitrange = fped->GetParameter(1)+2*fit->GetParameter(3); //testing
	    if(500<maxfitrange) maxfitrange = 500;
	    double minfitrange = 0.;

	    //cout<<fit->GetParameter(3)<<" "<<fped->GetParameter(2)*2<<endl;
	    if(fit->GetParameter(3)>fped->GetParameter(2)*2.001){ //If fitted gain is greater than minimum allowed value
	      minfitrange = fped->GetParameter(1)+ 0.75*fit->GetParameter(3); 
	      //cout<<"minfitrange "<<minfitrange<<endl;
	      fit->SetParLimits(0, 0, fit->GetParameter(0)*1.1);
	      fit->SetParLimits(4, fped->GetParameter(2)*1.01, fit->GetParameter(4)*1.01);
	      hspe->Fit(fit, "MNQL", "", minfitrange, maxfitrange);
	      minfitrange = fped->GetParameter(1)+ 0.75*fit->GetParameter(3); 
	      //cout<<"minfitrange "<<minfitrange<<endl;
	      hspe->Fit(fit, "MNQL", "", minfitrange, maxfitrange);
	    }

	    //calculate NDOF of fit excluding bins with 0 entries
	    int myNDOF=-3; //three free parameters
	    for(int j=hspe->FindBin(minfitrange); j<=hspe->FindBin(maxfitrange); j++) { //loop through fitted spe bins
	      if(hspe->GetBinContent(j)) myNDOF++;
            } //loop through fitted spe bins


	    //cout<<"estimate of gain "<<gain_est<<", fit "<<fit->GetParameter(3)<<endl;
	    //cout<<"SPE width "<<fit->GetParameter(4)<<endl;
	    //cout<<"Fit normalization constant: estimate "<<mu<<" fit "<<fit->GetParameter(0)<<endl;
	    //cout<<spename<<endl;


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


	    //output calibrations constants
	    //constants_file<<endl<<"LED_amplitude HV Spigot Channel Ped_mean Ped_mean_err Ped_RMS  Ped_RMS_err SPEPeak_RMS SPEPeak_RMS_err Gain Gain_err Normalized_Chi2 MeanPE_fit MeanPE_fit_err MeanPE_estimate PE5flag"<<endl;
	    constants_file<<LED_amp<<" "<<HV<<" "<<iSpig<<" "<<QIECh[i]<<" "<<scale*fped->GetParameter(1)<<" "<<scale*fped->GetParError(1)<<" "<<scale*fped->GetParameter(2)<<" "<<scale*fped->GetParError(2)<<" "<<scale*fit->GetParameter(4)<<" "<<scale*fit->GetParError(4)<<" "<<scale*fit->GetParameter(3)*fC2electrons<<" "<<scale*fit->GetParError(3)*fC2electrons<<" "<<fit->GetChisquare()/myNDOF/*fit->GetNDF()*/<<" "<<fit->GetParameter(0)<<" "<<fit->GetParError(0)<<" "<<mu<<" "<<PE5flag<<endl;

//	    cout<<LED_amp<<" "<<HV<<" "<<iSpig<<" "<<QIECh[i]<<" "<<scale*fped->GetParameter(1)<<" "<<scale*fped->GetParError(1)<<" "<<scale*fped->GetParameter(2)<<" "<<scale*fped->GetParError(2)<<" "<<scale*fit->GetParameter(4)<<" "<<scale*fit->GetParError(4)<<" "<<scale*fit->GetParameter(3)*fC2electrons<<" "<<scale*fit->GetParError(3)*fC2electrons<<" "<<fit->GetChisquare()/fit->GetNDF()<<" "<<fit->GetChisquare()<<" "<<fit->GetNDF()<<" "<<myNDOF<<" "<<fit->GetParameter(0)<<" "<<fit->GetParError(0)<<" "<<mu<<endl;
	   
	    //extra_file<<endl<<"LED_amplitude HV Spigot Channel SignalAvg_inFitRange FitAvg_inFitRange SignalInt_inFitRange FitInt_inFitRange PEge5Int Gain(fC) Gain_err MeanPE_fit MeanPE_fit_err MeanPE_estimate"<<endl;
	    //extra_file<<LED_amp<<" "<<HV<<" "<<iSpig<<" "<<QIECh[i]<<" "<<scale*havg<<" "<<scale*favg<<" "<<hint<<" "<<fint<<" "<<PE5int<<" "<<scale*fit->GetParameter(3)<<" "<<scale*fit->GetParError(3)<<" "<<fit->GetParameter(0)<<" "<<fit->GetParError(0)<<" "<<mu<<" "<<PE5flag<<endl;


	    //cout<<"# Spigot Channel Ped_mean Ped_RMS SPE_PeakRMS Gain Normalized_Chi2 Avg_PE"<<endl;
	    //cout<<iSpig<<" "<<i<<" "<<scale*fped->GetParameter(1)<<" "<<scale*fped->GetParameter(2)<<" "<<scale*fit->GetParameter(4)<<" "<<scale*fit->GetParameter(3)*fC2electrons<<" "<<fit->GetChisquare()/fit->GetNDF()<<" "<<fit->GetParameter(0)<<endl<<endl;
	    
	    
	  
	    if(iSpig==0) c1->cd(i+1);
	    else if(iSpig==1) c2->cd(i+1);
	    gPad->SetBorderMode(0);
	    gPad->SetBorderSize(0);
	    gPad->SetRightMargin(0.01);
	    gPad->SetBottomMargin(0.1);
	    gPad->SetLogy(true);
	    hspe->GetXaxis()->SetRangeUser(0, /*300*/508);
	    hspe->SetLineColor(kBlue);
	    hspe->DrawClone("hist");
	    fit->SetLineWidth(2);
	    fit->Draw("same");
	   
	    /*
	    c3->cd();
	    gPad->SetBorderMode(0);
	    gPad->SetBorderSize(0);
	    gPad->SetRightMargin(0.01);
	    gPad->SetBottomMargin(0.1);
	    gPad->SetLogy(true);
	    hspe->GetXaxis()->SetRangeUser(0, 300);
	    hspe->Draw();
	    fit->SetLineWidth(1);
	    fit->SetLineColor(2);
	    fit->Draw("same");
            stringstream plot_name;
	    plot_name<<"Plots/Channels/Spig"<<iSpig<<"_Chan"<<QIECh[i]<<"_HV"<<HV<<"_LED"<<LED_amp<<"_Run_"<<run<<"_SPEFit.pdf";
	    //c3->SaveAs(plot_name.str().c_str());
	    plot_name.str( std::string() );
	    */
	  }
      }

      if(drawflag) { //draw plots of fit if data for the HV is present
	stringstream plot_name;
	plot_name<<"Plots/SPEFits_Spigot0_Run_"<<run<<"_HV"<<HV<<".pdf";
	c1->SaveAs(plot_name.str().c_str());
	plot_name.str( std::string() );
	plot_name<<"Plots/SPEFits_Spigot1_Run_"<<run<<"_HV"<<HV<<".pdf";
	c2->SaveAs(plot_name.str().c_str());
      }

    } //HV loop

    delete c1;
    delete c2;
    //delete c3;

    constants_file.close();
}

Double_t FitFun(Double_t *x, Double_t *par)
{
// Spectra fit function: Pedestal Gaussian + asymmetric 1PE + 2PE +3PE peaks
        
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

