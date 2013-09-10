#include <iostream>
#include <TROOT.h>
#include <TRandom.h>
#include <TH1.h>
#include <TH1F.h>
#include <vector>

using namespace std;

void sig(float back, float eback, float data);

void dirty()

{
  sig(0.7,0.05,2.);
}

void sig(float back, float eback, float data)

{
 int nevt = 10000000;

 float xmax = data*5.;
 int nbins = (int)xmax+1;

 TCanvas* c = new TCanvas();
 c->cd();

 TH1D* h = new TH1D("h","h",nbins,0.,xmax);

 TRandom ran;

 TRandom reso;

 for (int i=0; i < nevt; i++) {

  float u_back = reso.Gaus(back,eback);
  int x = ran.Poisson(u_back);
  float y=x;
   h->Fill(y); 
 }
 
 h->Draw();
 //c->Update();

 float s=0;
 for (int j=(int)data; j <xmax; j++) {
  int bin = h -> FindBin(j);
  float uu = h -> GetBinContent(bin);
  s += uu;
 }

 cout << "Background : " << back << " +/- " << eback << endl;
 cout << "Data : " << data << endl; 

 cout << "Fluctuation probability : ";
 cout << s/nevt << endl;

 float proba = s/nevt;

 // --- Deviation in terms of sigma :

 //delete h;

 for (int k=0; k < 1000; k++) {
  float xx = (float)(k+1)/100.;
  float vv = TMath::Erf(xx);
     // -- searching for xx where 1 - Erf(xx) = Proba *2
  if (vv > 1.-proba*2.) {
    cout << "corresponds to " << xx*sqrt(2) << " sigma " << endl;
    break;
  }
 }

}
