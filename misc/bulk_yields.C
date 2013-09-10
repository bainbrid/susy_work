#include "TMultiGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveText.h"
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;
void setTDRStyle();

void qcd() {
   
  setTDRStyle();

  std::vector<double> bins;
  bins.push_back(275.);
  bins.push_back(325.);
  for ( int i = 0; i < 6; ++i ) { bins.push_back(375.+100.*i); }

  std::vector<double> widths;
  widths.push_back(25.);
  widths.push_back(25.);
  for ( int i = 0; i < 6; ++i ) { widths.push_back(50.); }

  std::vector<double> effs;
  effs.push_back(1.13);
  effs.push_back(0.83);
  effs.push_back(0.72);
  effs.resize(8,0.84);
//   effs.push_back(0.84);
//   effs.push_back(1.13);
//   effs.push_back(1.32);
//   effs.push_back(1.41);
//   effs.push_back(2.09);

  std::vector<double> errs;
  errs.push_back(0.34);
  errs.push_back(0.24);
  errs.push_back(0.19);
  errs.resize(8,0.14);
//   errs.push_back(0.14);
//   errs.push_back(0.24);
//   errs.push_back(0.35);
//   errs.push_back(0.50);
//   errs.push_back(0.80);

  std::vector<double> entries;
  entries.push_back(46272.);
  entries.push_back(43592.);
  entries.push_back(58746.);
  entries.push_back(78644.);
  entries.push_back(29516.);
  entries.push_back(11668.);
  entries.push_back(5416.);
  entries.push_back(5368.);

  std::vector<double> prescales;
  prescales.push_back(2000.);
  prescales.push_back(1000.);
  prescales.push_back(500.);
  prescales.push_back(125.);
  prescales.push_back(125.);
  prescales.push_back(125.);
  prescales.push_back(125.);
  prescales.push_back(125.);
  
  TCanvas* c = new TCanvas("tmp","tmp",900,600);
  
  TGraphAsymmErrors* gr1 = new TGraphAsymmErrors(8);
  gr1->SetName("Graph1");
  gr1->SetTitle("");
  gr1->SetFillColor(1);
  gr1->SetLineWidth(2);
  gr1->SetMarkerStyle(24);
  gr1->SetMarkerColor(1);
  gr1->SetMarkerSize(2.);

  TGraphAsymmErrors* gr2 = new TGraphAsymmErrors(8);
  gr2->SetName("Graph2");
  gr2->SetTitle("");
  gr2->SetFillColor(1);
  gr2->SetLineWidth(2);
  gr2->SetMarkerStyle(20);
  gr1->SetMarkerColor(1);
  gr2->SetMarkerSize(1.5);

  for ( int i = 0; i < bins.size(); ++i ) {
    double x = bins[i]+widths[i];
    double y = entries[i]*prescales[i];
    double ex = widths[i];
    double ey = sqrt(entries[i])*prescales[i];
    
    gr1->SetPoint(i+1,x,y);
    gr1->SetPointError(i+1,0.,0.,ey,ey);

    y = y / effs[i];
    ey = y * sqrt( (ey/y)*(ey/y) + (errs[i]/effs[i])*(errs[i]/effs[i]) );
    
    gr2->SetPoint(i+1,x,y);
    gr2->SetPointError(i+1,ex,ex,ey,ey);

    std::cout << " HT= " << bins[i] << " yield= " << y << " +/- " << ey << std::endl; 
 
  }
  
  TMultiGraph* mg = new TMultiGraph();
  mg->Add(gr1,"p");
  mg->Add(gr2,"p");
  mg->Draw("ap");
  mg->GetXaxis()->SetTitle("H_{T} (GeV)");
  mg->GetYaxis()->SetTitle("Bulk yields");
  mg->GetYaxis()->SetRangeUser(1.e5,1.e9);
  mg->GetXaxis()->SetRangeUser(275.,975.);

  TF1* fit = new TF1("fit","expo",275.,875.);
  gr2->Fit(fit,"R");
  fit->Draw("same");
  
  TLegend* leg = new TLegend( 0.45, 0.70, 0.85, 0.80 );
  leg->SetFillColor(0);
  leg->SetLineColor(0); 
  leg->SetShadowColor(0); 
  leg->SetTextSize(0.035);
  leg->AddEntry(gr1,"Raw counts weighted by prescales","p");
  leg->AddEntry(gr2,"Trigger efficiency corrected","p");
  leg->Draw("same");

  std::stringstream sss;
  sss << "CMS, 1.5 fb^{-1}, #sqrt{s} = 8 TeV";
  TLatex* tex = new TLatex(0.17,0.88,sss.str().c_str());
  tex->SetNDC();
  tex->SetTextSize(0.035);
  tex->Draw();

  c->SetLogy();

  return;
  
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

void setTDRStyle() {

  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  // For the histo:
  tdrStyle->SetHistFillColor(63);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  //  tdrStyle->SetEndErrorSize(0);
  tdrStyle->SetErrorX(0.);
  //  tdrStyle->SetErrorMarker(20);
  
  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

  // Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.13);
  tdrStyle->SetPadRightMargin(0.05);

  // For the Global title:

  //  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);


  // For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.05);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  //@@ was: tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetNdivisions(510, "YZ"); 
  tdrStyle->SetNdivisions(1005, "X");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  // Postscript options:
  // tdrStyle->SetPaperSize(15.,15.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);


  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();

}

