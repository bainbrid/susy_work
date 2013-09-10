#include <TF1.h>
#include <TH1F.h>
#include <TMath.h>
#include <TROOT.h>

Double_t myfunction(Double_t *x, Double_t *par)
{
  Float_t xx =x[0];
  //Double_t f = TMath::Abs(par[0]*sin(par[1]*xx)/xx);
  Double_t f = TMath::Abs(par[0]+TMath::Exp(par[1]+par[2]*xx));
  return f;
}

void myfunc()
{
  TF1 *f1 = new TF1("myfunc",myfunction,0,10,3);
  f1->SetParameters(2,1,-1);
  f1->SetParNames("constant","const","decay");
  f1->Draw();
}

void fit()
{
  TH1F* h1 = new TH1F("h1","test",100,0,10);
  TF1* f1 = (TF1*)gROOT->GetFunction("myfunction");
  h1->FillRandom("myfunc",20000);
  f1->SetParameters(2,1,-1);
  h1->Fit("myfunction");
  h1->Draw();
}

