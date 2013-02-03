#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"

#include <vector>
#include <utility>
#include <algorithm>
#include <map>
#include <fstream>
#include <sstream>

Int_t Nev;

Double_t FitFun(Double_t *x, Double_t *par);

void SPEFit(bool isHF=true, char * fname = "output_700V.root", int run=-1, int HV=-1, int LED_amp=-1, double cutmin = 30.0, double cutmax = 250.0)
{

    //set file names
    stringstream out_fname;
    out_fname<<"SPEconstants_Run_"<<run<<".txt";

    ofstream  constants_file(out_fname.str().c_str(),ios_base::trunc); 
    constants_file<<"Run "<<run<<endl;
    constants_file<<"type SPE"<<endl;
    constants_file<<"HV "<<HV<<endl;
    constants_file<<"LED_amplitude "<<LED_amp<<endl<<endl;

    constants_file<<endl<<"Spigot Channel Ped_mean Ped_RMS SPEPeak_RMS Gain Normalized_Chi2 MeanPE_fit MeanPE_estimate"<<endl;

    out_fname.str("");
    out_fname<<"SPEdistributions_Run_"<<run<<".txt";

    //ofstream  plots_file(out_fname.str().c_str(),ios_base::trunc);
    //plots_file<<"Run "<<run<<endl;
    //plots_file<<"type SPE"<<endl<<endl;
    //plots_file<<"HV "<<HV<<endl<<endl;

    double scale = 1.0;
    if(isHF) scale = 2.6; //Need to scale up HF charge
    double fC2electrons = 6240.; //convert fC to #electrons

    char spename[128], pedname[128];
    
    TFile *tf = new TFile(fname);
    
    TCanvas *c1 = new TCanvas("c1","c1",1200,700);
    c1->Divide(6,4);
    c1->SetBorderMode(0);
    c1->SetBorderSize(0);
    TCanvas *c2 = new TCanvas("c2","c2",1200,700);
    c2->Divide(6,4);
    c2->SetBorderMode(0);
    c2->SetBorderSize(0);  

    for (int iSpig = 0; iSpig < 2; iSpig++) {
      for(int i = 1; i <= 24; i++)
	{
	  if(i < 10) sprintf(spename, "kplotter/FED_700/Spigot 0%d/LED Plots/Sum/QIEsSumLED_S0%d_Ch0%d", iSpig, iSpig, i);
	  else sprintf(spename, "kplotter/FED_700/Spigot 0%d/LED Plots/Sum/QIEsSumLED_S0%d_Ch%d", iSpig, iSpig, i);
	  TH1 *hspe = (TH1*)tf->Get(spename);
	  if(i < 10) sprintf(pedname, "kplotter/FED_700/Spigot 0%d/Pedestal Plots/Sum/QIEsSumPedestal_S0%d_Ch0%d", iSpig, iSpig, i);
	  else sprintf(pedname, "kplotter/FED_700/Spigot 0%d/Pedestal Plots/Sum/QIEsSumPedestal_S0%d_Ch%d", iSpig, iSpig, i);
	  TH1 *hped = (TH1*)tf->Get(pedname);
	  
	  Nev = hspe->Integral()*hspe->GetXaxis()->GetBinWidth(1);
	  //cout<<"#Events in LED distribution "<<hspe->Integral()<<endl;

	  TF1 *fped = new TF1("fped","gaus",0, 50);
	  hped->Fit(fped,"NQR");
	  double pploc = fped->GetParameter(1), ppwidth = fped->GetParameter(2);
	  cout<<"Ped only: ped mean "<<fped->GetParameter(1)<<", ped width "<<fped->GetParameter(2)<<" normalization "<<fped->GetParameter(0)<<endl;
	  hspe->Fit(fped, "NQ", "", pploc - 3*ppwidth, pploc + ppwidth);
	  cout<<"SPE distribution: ped mean "<<fped->GetParameter(1)<<", ped width "<<fped->GetParameter(2)<<" normalization "<<fped->GetParameter(0)<<endl;
	  
	  
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
	  //fit->SetParameters(1, 20, 1, 40, 2);
	  //fit->SetParameters(1, 20, 1, 22, 2);
	  //fit->SetParameters(0.5, 20, 1, max_SPE_location-0.0*fped->GetParameter(1), 2);

double mu = - log(fped->Integral(0,100)/Nev);

	  if(max_SPE_bin > (minbin+1)) fit->SetParameters(/*0.5*/mu, 20, 1, max_SPE_location-1.0*fped->GetParameter(1), 10);
	  else fit->SetParameters(/*0.5*/mu, 20, 1, 2.0*fped->GetParameter(2), 10); //case of no clear minimum; start looking for SPE peak at 2sigma away from pedestal peak
	  fit->FixParameter(1, fped->GetParameter(1));
	  fit->FixParameter(2, fped->GetParameter(2));
          //fit->SetParameter(1, fped->GetParameter(1));
          //fit->SetParameter(2, fped->GetParameter(2));

	  fit->SetParLimits(3, ppwidth*3, cutmax-pploc);
	  //hspe->Fit(fit, "MNQB", "", cutmin, cutmax);
	 
          fit->SetParLimits(4, fped->GetParameter(2), 100);

 
	  hspe->Fit(fit, "MNQB", "", fped->GetParameter(1)+fped->GetParameter(2), cutmax);
	  double maxfitrange = fped->GetParameter(1)+3*fped->GetParameter(3)+fped->GetParameter(4);
	  //if(cutmax<maxfitrange) maxfitrange = cutmax;
	  if(500<maxfitrange) maxfitrange = 500;
	  hspe->Fit(fit, "MNQB", "", fped->GetParameter(1)+fped->GetParameter(2), maxfitrange);
	  /*
	  hspe->Fit(fit, "MNQB", "", fped->GetParameter(1)-fped->GetParameter(2), cutmax);
	  double maxfitrange = fped->GetParameter(1)+3*fped->GetParameter(3)+fped->GetParameter(4);
	  if(cutmax<maxfitrange) maxfitrange = cutmax;
	  hspe->Fit(fit, "MNQB", "", fped->GetParameter(1)-fped->GetParameter(2), maxfitrange);
	  */

	  cout<<"estimate of gain "<<max_SPE_location-fped->GetParameter(1)<<", fit "<<fit->GetParameter(3)<<endl;
	  cout<<"Fit normalization constant: estimate "<<mu<<" fit "<<fit->GetParameter(0)<<endl;
	  cout<<spename<<endl;

	  //output calibrations constants
	  //constants_file<<"# Spigot Channel Ped_mean Ped_RMS SPE_PeakRMS Gain Normalized_Chi2 Avg_PE_fit Avg_PE_est"<<endl;
	  constants_file<<iSpig<<" "<<i<<" "<<scale*fped->GetParameter(1)<<" "<<scale*fped->GetParameter(2)<<" "<<scale*fit->GetParameter(4)<<" "<<scale*fit->GetParameter(3)*fC2electrons<<" "<<fit->GetChisquare()/fit->GetNDF()<<" "<<fit->GetParameter(0)<<" "<<mu<<endl;
	  
	  cout<<"# Spigot Channel Ped_mean Ped_RMS SPE_PeakRMS Gain Normalized_Chi2 Avg_PE"<<endl;
	  cout<<iSpig<<" "<<i<<" "<<scale*fped->GetParameter(1)<<" "<<scale*fped->GetParameter(2)<<" "<<scale*fit->GetParameter(4)<<" "<<scale*fit->GetParameter(3)*fC2electrons<<" "<<fit->GetChisquare()/fit->GetNDF()<<" "<<fit->GetParameter(0)<<endl<<endl;

	  //dump plot and fit to text file
	  //plots_file<<"Spigot "<<iSpig<<" Channel "<<i<<endl; 
	  //plots_file<<"Energy Sum"<<endl;
	  Nbins = hspe->GetBin(cutmax);
	  //for(int ia=1; ia < Nbins; ia++) plots_file<<scale*hspe->GetBinCenter(ia)<<" "<<hspe->GetBinContent(ia)<<endl;
	  //plots_file<<endl<<"SPE fit"<<endl;
	  //for(int ia=1; ia < Nbins; ia++) plots_file<<scale*hspe->GetBinCenter(ia)<<" "<<fit->Eval(hspe->GetBinCenter(ia))<<endl;
	  //plots_file<<endl;
	  
	  
	  if(iSpig==0) c1->cd(i);
	  else if(iSpig==1) c2->cd(i);
	  gPad->SetBorderMode(0);
	  gPad->SetBorderSize(0);
	  gPad->SetRightMargin(0.01);
	  gPad->SetBottomMargin(0.1);
	  gPad->SetLogy(true);
	  hspe->GetXaxis()->SetRangeUser(0, /*300*/200);
	  hspe->Draw();
	  fit->Draw("same");
	  
	}
    }

    stringstream plot_name;
    plot_name<<"Plots/SPEFits_Spigot0_Run_"<<run<<".pdf";
    c1->SaveAs(plot_name.str().c_str());
    plot_name.str( std::string() );
    plot_name<<"Plots/SPEFits_Spigot1_Run_"<<run<<".pdf";
    c2->SaveAs(plot_name.str().c_str());

    //plots_file.close();    
    constants_file.close();

}

Double_t FitFun(Double_t *x, Double_t *par)
{
// Spectra fit function: Pedestal Gaussian + asymmetric 1PE + 2PE +3PE peaks
        
  Double_t sum,xx,A0,C0,r0,sigma0,mean1,sigma1,A1,C1,r1,mean2,sigma2,A2,C2,r2,mean3,sigma3,A3,C3,r3;
        
  const Double_t k0=2.0, k1=1.6, k2=2.0;
        
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
  C3 = 1/(sigma3*sqrt(2*3.14159));
  r3 = ((xx-mean3)/sigma3);
  if(xx>par[1]-3.0*par[2]) sum += C3*A3*TMath::Exp(-0.5*r3*r3);

  return sum;
}

