TGraphErrors* LP_get5fbObs_NLO_HT500(){

Int_t nl = 200;
Double_t xl[200]; Double_t yl[200];
for (int i=0; i<nl; i++) { xl[i]=0.; yl[i]=0.; }
Int_t NN = 0;
  xl[NN] = 0;     yl[NN] = 540; NN++; 
  xl[NN] = 200;   yl[NN] = 480; NN++; 
  xl[NN] = 500;   yl[NN] = 480; NN++; 
  xl[NN] = 700;   yl[NN] = 470; NN++; 
  xl[NN] = 1000;  yl[NN] = 360; NN++; 
  xl[NN] = 1200;  yl[NN] = 280; NN++;   
  xl[NN] = 1500;  yl[NN] = 260; NN++; 
  xl[NN] = 1750;  yl[NN] = 240; NN++;   
  xl[NN] = 2000;  yl[NN] = 230; NN++; 
  Double_t exl[21];
  Double_t eyl[21];

  TGraphErrors* gr1 = new TGraphErrors(NN,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);

  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);

  return gr1;
}


TGraphErrors* LP_get5fbUp_NLO_HT500(){

Int_t nl = 200;
Double_t xl[200]; Double_t yl[200];
for (int i=0; i<nl; i++) { xl[i]=0.; yl[i]=0.; }
Int_t NN = 0;
  xl[NN] = 0;     yl[NN] = 560; NN++; 
  xl[NN] = 200;   yl[NN] = 520;     NN++; 
  xl[NN] = 500;   yl[NN] = 520;     NN++; 
  xl[NN] = 700;   yl[NN] = 510;  NN++; 
  xl[NN] = 1000;  yl[NN] = 400;     NN++; 
  xl[NN] = 1200;  yl[NN] = 320;  NN++;   
  xl[NN] = 1500;  yl[NN] = 280;     NN++; 
  xl[NN] = 1750;  yl[NN] = 280;     NN++;   
  xl[NN] = 2000;  yl[NN] = 270;   NN++; 
  Double_t exl[21];
  Double_t eyl[21];

  TGraphErrors* gr1 = new TGraphErrors(NN,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);

  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);

  return gr1;

}


TGraphErrors* LP_get5fbLow_NLO_HT500(){

Int_t nl = 200;
Double_t xl[200]; Double_t yl[200];
for (int i=0; i<nl; i++) { xl[i]=0.; yl[i]=0.; }
Int_t NN = 0;
  xl[NN] = 0;     yl[NN] = 500; NN++; 
  xl[NN] = 200;   yl[NN] = 440;     NN++; 
  xl[NN] = 500;   yl[NN] = 420;  NN++; 
  xl[NN] = 700;   yl[NN] = 380;     NN++; 
  xl[NN] = 1000;  yl[NN] = 320; NN++; 
  xl[NN] = 1200;  yl[NN] = 240;  NN++; 
  xl[NN] = 1500;  yl[NN] = 200;     NN++; 
  xl[NN] = 1750;  yl[NN] = 200;     NN++;   
  xl[NN] = 2000;  yl[NN] = 190;   NN++; 
  Double_t exl[21];
  Double_t eyl[21];

  TGraphErrors* gr1 = new TGraphErrors(NN,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);

  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);

  return gr1;

}


TGraphErrors* LP_get5fbExp_NLO_HT500(){

Int_t nl = 200;
Double_t xl[200]; Double_t yl[200];
for (int i=0; i<nl; i++) { xl[i]=0.; yl[i]=0.; }
Int_t NN = 0;
  xl[NN] = 0;     yl[NN] = 540; NN++; 
  xl[NN] = 200;   yl[NN] = 480;     NN++; 
  xl[NN] = 500;   yl[NN] = 480;     NN++; 
  xl[NN] = 700;   yl[NN] = 460;  NN++; 
  xl[NN] = 1000;  yl[NN] = 380;  NN++; 
  xl[NN] = 1200;  yl[NN] = 280;     NN++;   
  xl[NN] = 1500;  yl[NN] = 260;  NN++; 
  xl[NN] = 1750;  yl[NN] = 240;     NN++; 
  xl[NN] = 2000;  yl[NN] = 230;   NN++; 
  Double_t exl[21];
  Double_t eyl[21];

  TGraphErrors* gr1 = new TGraphErrors(NN,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);

  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);

  return gr1;

}


//This graph encloses the +/-1 sigma expected band.
//May need to adjust slightly so the filled region matches +/-1 sigma curves.
TGraphErrors* LP_get5fbEnclosed_NLO_HT500(){

Int_t nl = 200;
Double_t xl[200]; Double_t yl[200];
for (int i=0; i<nl; i++) { xl[i]=0.; yl[i]=0.; }
Int_t NN = 0;
  xl[NN] = 0;     yl[NN] = 500; NN++; 
  xl[NN] = 200;   yl[NN] = 440;     NN++; 
  xl[NN] = 500;   yl[NN] = 420;  NN++; 
  xl[NN] = 700;   yl[NN] = 380;     NN++; 
  xl[NN] = 1000;  yl[NN] = 320; NN++; 
  xl[NN] = 1200;  yl[NN] = 240;  NN++; 
  xl[NN] = 1300;  yl[NN] = 215;  NN++; 
  xl[NN] = 1400;  yl[NN] = 200;  NN++; 
  xl[NN] = 1500;  yl[NN] = 200;     NN++; 
  xl[NN] = 1750;  yl[NN] = 200;     NN++;   
  xl[NN] = 2000;  yl[NN] = 190;   NN++; 
  xl[NN] = 2000;  yl[NN] = 270;   NN++; 
  xl[NN] = 1750;  yl[NN] = 280;     NN++;   
  xl[NN] = 1400;  yl[NN] = 285;  NN++; 
  xl[NN] = 1300;  yl[NN] = 300;  NN++; 
  xl[NN] = 1200;  yl[NN] = 320;  NN++;   
  xl[NN] = 1000;  yl[NN] = 400;     NN++; 
  xl[NN] = 800;  yl[NN] = 480;     NN++; 
  xl[NN] = 700;   yl[NN] = 510;  NN++; 
  xl[NN] = 500;   yl[NN] = 520;     NN++; 
  xl[NN] = 200;   yl[NN] = 520;     NN++; 
  xl[NN] = 0;     yl[NN] = 560; NN++; 
  xl[NN] = 0;     yl[NN] = 500; NN++; 
  Double_t exl[21];
  Double_t eyl[21];

  TGraphErrors* gr1 = new TGraphErrors(NN,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);

  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);

  return gr1;
}


TGraphErrors* LS_get5fbObs_NLO_HT500(){

Int_t nl = 200;
Double_t xl[200]; Double_t yl[200];
for (int i=0; i<nl; i++) { xl[i]=0.; yl[i]=0.; }
Int_t NN = 0;
  xl[NN] = 0;  yl[NN] = 560; NN++; 
  xl[NN] = 200;  yl[NN] = 490; NN++; 
  xl[NN] = 500;  yl[NN] = 460; NN++; 
  xl[NN] = 700;  yl[NN] = 445; NN++; 
  xl[NN] = 1000;  yl[NN] = 290; NN++; 
  xl[NN] = 1200;  yl[NN] = 245; NN++;   
  xl[NN] = 1500;  yl[NN] = 215; NN++; 
  xl[NN] = 1750;  yl[NN] = 208; NN++;   
  xl[NN] = 2000;  yl[NN] = 200; NN++; 
  Double_t exl[21];
  Double_t eyl[21];

  TGraphErrors* gr1 = new TGraphErrors(NN,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);

  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);

  return gr1;
}


TGraphErrors* LS_get5fbUp_NLO_HT500(){

Int_t nl = 200;
Double_t xl[200]; Double_t yl[200];
for (int i=0; i<nl; i++) { xl[i]=0.; yl[i]=0.; }
Int_t NN = 0;
  xl[NN] = 0;  yl[NN] = 600; NN++; 
  xl[NN] = 200;  yl[NN] = 500+5; NN++; 
  xl[NN] = 500;  yl[NN] = 470+5; NN++; 
  xl[NN] = 700;  yl[NN] = 460+10; NN++; 
  xl[NN] = 1000;  yl[NN] = 330+10; NN++; 
  xl[NN] = 1200;  yl[NN] = 280+5; NN++;   
  xl[NN] = 1500;  yl[NN] = 245+5; NN++; 
  xl[NN] = 1750;  yl[NN] = 242+5; NN++;   
  xl[NN] = 2000;  yl[NN] = 235+5; NN++; 
  Double_t exl[21];
  Double_t eyl[21];

  TGraphErrors* gr1 = new TGraphErrors(NN,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);

  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);

  return gr1;

}


TGraphErrors* LS_get5fbLow_NLO_HT500(){

Int_t nl = 200;
Double_t xl[200]; Double_t yl[200];
for (int i=0; i<nl; i++) { xl[i]=0.; yl[i]=0.; }
Int_t NN = 0;
  xl[NN] = 0;  yl[NN] = 500; NN++; 
  xl[NN] = 200;  yl[NN] = 425; NN++; 
  xl[NN] = 500;  yl[NN] = 400+5; NN++; 
  xl[NN] = 700;  yl[NN] = 345+10; NN++; 
  xl[NN] = 1000;  yl[NN] = 235; NN++; 
  xl[NN] = 1200;  yl[NN] = 190; NN++; 
  xl[NN] = 1500;  yl[NN] = 175+5; NN++; 
  xl[NN] = 1750;  yl[NN] = 172+5; NN++;   
  xl[NN] = 2000;  yl[NN] = 170+5; NN++; 
  Double_t exl[21];
  Double_t eyl[21];

  TGraphErrors* gr1 = new TGraphErrors(NN,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);

  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);

  return gr1;

}


TGraphErrors* LS_get5fbExp_NLO_HT500(){

Int_t nl = 200;
Double_t xl[200]; Double_t yl[200];
for (int i=0; i<nl; i++) { xl[i]=0.; yl[i]=0.; }
Int_t NN = 0;
  xl[NN] = 0;  yl[NN] = 540; NN++; 
  xl[NN] = 200;  yl[NN] = 460+5; NN++; 
  xl[NN] = 500;  yl[NN] = 440; NN++; 
  xl[NN] = 700;  yl[NN] = 425; NN++; 
  xl[NN] = 1000;  yl[NN] = 285; NN++; 
  xl[NN] = 1200;  yl[NN] = 245; NN++;   
  xl[NN] = 1500;  yl[NN] = 220; NN++; 
  xl[NN] = 1750;  yl[NN] = 212; NN++; 
  xl[NN] = 2000;  yl[NN] = 205; NN++; 
  Double_t exl[21];
  Double_t eyl[21];

  TGraphErrors* gr1 = new TGraphErrors(NN,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);

  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);

  return gr1;

}





