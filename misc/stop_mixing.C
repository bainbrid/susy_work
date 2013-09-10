#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TLegend.h"

int stop_mixing() {

  TCanvas* c1 = new TCanvas("c1","");
  c1->SetFillColor(0);

  TStyle* tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
  tdrStyle->SetTitleSize(0.045, "XYZ");
  tdrStyle->cd();
  
  const Int_t n = 26;
  
  Double_t x[n] = {0.,2.,4.,6.,8.,10.,12.,14.,16.,18.,20.,22.,24.,26.,
		   28.,30.,32.,34.,36.,38.,40.,42.,44.,46.,48.,50.};

  Double_t y1[n] = {1.00,1.00,1.00,1.00,1.00,1.00,0.99,0.99,0.97,0.95,0.90,0.84,0.75,
		    0.65,0.55,0.45,0.35,0.27,0.21,0.16,0.13,0.10,0.07,0.06,0.05,0.04};
  
  Double_t y2[n] = {1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,
 		    0.99,0.99,0.98,0.98,0.97,0.96,0.94,0.92,0.90,0.87,0.83,0.79,0.75};
  
  for ( int i = 0; i < n; ++i ) {
    y1[i] = 1. - y1[i];
    y2[i] = 1. - y2[i];
  }

  TGraph* gr1 = new TGraph(n,x,y1);
  gr1->SetTitle("");
  gr1->GetYaxis()->SetTitle("BF(stop #rightarrow b,W*,LSP)");
  //gr1->GetYaxis()->SetTitle("BF(stop #rightarrow charm,LSP)");
  gr1->GetXaxis()->SetTitle("#DeltaM(stop,LSP)");
  gr1->GetYaxis()->SetRangeUser(0.,1.);
  gr1->GetXaxis()->SetRangeUser(0.,50.);
  gr1->SetLineColor(2);
  gr1->SetLineStyle(1);
  gr1->SetLineWidth(3);
  gr1->Draw("AL");

  TGraph* gr2 = new TGraph(n,x,y2);
  gr2->SetTitle("");
  gr2->GetYaxis()->SetTitle("BF(stop #rightarrow b,W*,LSP)");
  //gr2->GetYaxis()->SetTitle("BF(stop #rightarrow charm,LSP)");
  gr2->GetXaxis()->SetTitle("#DeltaM(stop,LSP)");
  gr2->GetYaxis()->SetRangeUser(0.,1.);
  gr2->GetXaxis()->SetRangeUser(0.,50.);
  gr2->SetLineColor(4);
  gr2->SetLineStyle(1);
  gr2->SetLineWidth(3);
  gr2->Draw("Lsame");

  TLegend* leg = new TLegend(0.20,0.20,0.40,0.40,NULL,"brNDC");
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->SetTextFont(62);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetFillColor(19);
  leg->SetFillStyle(0);
  leg->AddEntry(gr2,"tan #beta = 30","l");
  leg->AddEntry(gr1,"tan #beta = 10","l");
  leg->Draw("same");

  return 0;

}
