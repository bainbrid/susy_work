#include "TMultiGraph.h"
#include "TLegend.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TStyle.h"
#include "TBox.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

// options
bool inclusive = false;
bool weighted = true;

typedef unsigned int uint;
typedef vector<double> vdouble;
typedef vector<string> vstring;
typedef vector<vdouble> vvdouble;

const uint nbins = 8;

void setTDRStyle();
void init( uint, vdouble& );
void init( uint, uint, vvdouble& );
uint extractData( vdouble& b, vvdouble& v, vvdouble& e, vstring&, bool incl );
void extractNumbers( vvdouble& v, vvdouble& e, vdouble& b, vdouble& m, vdouble& s );

void closure() {
   
  setTDRStyle();

  vdouble markers;
  markers.push_back(24);
  markers.push_back(25);
  markers.push_back(26);
  markers.push_back(28);
  markers.push_back(30);

  // Bins
  vdouble bins;
  bins.push_back(275.);
  bins.push_back(325.);
  for ( uint i = 0; i < 6; ++i ) { bins.push_back(375.+100.*i); }
  bins.push_back(975.);
  
  // Bin errors
  vdouble widths;
  widths.push_back(25.);
  widths.push_back(25.);
  for ( uint i = 0; i < 6; ++i ) { widths.push_back(50.); }
  
  // Extract data from histograms
  vvdouble values; 
  vvdouble errors;
  vstring titles;
  const uint nhistos = extractData( bins, values, errors, titles, inclusive );

  // Extract central values and errors for each HT bin
  vdouble b1; for ( uint i = 0; i < nbins; ++i ) { b1.push_back(i); }
  vdouble m1; 
  vdouble s1;
  extractNumbers( values, errors, b1, m1, s1 ); 

  // Extract central values and errors for each HT region
  vdouble b2; b2.push_back(0); b2.push_back(4); b2.push_back(6);
  //vdouble b2; b2.push_back(0); b2.push_back(5); 
  vdouble m2; 
  vdouble s2;
  extractNumbers( values, errors, b2, m2, s2 ); 
  
  for ( uint i = 0; i < m2.size(); ++i ) { 
    double min = fabs( m2[i] - s2[i] );
    double max = fabs( m2[i] + s2[i] );
    if (weighted) { s2[i] *= 3.; }
    m2[i] = 0.;
  }
  
  // Plot numbers
  TLegend* leg = new TLegend( 0.16, 0.91-0.04*(nhistos+2), 0.35, 0.91 );
  leg->SetFillColor(0);
  leg->SetLineColor(0); 
  leg->SetShadowColor(0); 
  leg->SetTextSize(0.027);

  TMultiGraph* mg = new TMultiGraph();
  int offset = 0;
  for ( uint i = 0; i < nhistos+2; ++i ) { 
    offset++;

    if ( i == 0 ) {
      
      vdouble x; for ( uint j = 0; j < b2.size(); ++j ) { x.push_back( bins[b2[j]] ); }
      vdouble xel(b2.size(),0.);
      vdouble xeh; for ( uint j = 0; j < b2.size(); ++j ) { xeh.push_back( ((j<b2.size()-1)?x[j+1]-x[j]:bins[nbins]-x[j]) ); }

      for ( uint j = 0; j < b2.size(); ++j ) { std::cout << " HT: " << x[j] << " error: " << 100.*s2[j] << endl; }
      
      TGraphAsymmErrors* tmp = new TGraphAsymmErrors( b2.size(),
						      &x.front(),
						      &m2.front(),
						      &xel.front(),
						      &xeh.front(),
						      &s2.front(),
						      &s2.front() );
      tmp->SetMarkerStyle(20);
      tmp->SetMarkerSize(1.5);
      tmp->SetFillColor(kGray);
      tmp->SetTitle("");
      mg->Add(tmp,"2");
      stringstream ss; ss << "Systematic uncertainty (" << (weighted?"#pm3#sigma)":"r.m.s.)");
      leg->AddEntry(tmp,ss.str().c_str(),"f");

  } else if ( i > 0 && i <= nhistos ) {
      
      //continue;
      
      vdouble e(widths.size(),0.);
      vdouble x(bins); for ( uint j = 0; j < x.size(); ++j ) { x[j] += widths[j] + 5.*(offset%2?1*offset/2:-1*offset/2); }
      
      TGraphErrors* tmp = new TGraphErrors( nbins,
					    &x.front(),
					    &values[i-1].front(),
					    &e.front(),
					    &errors[i-1].front() );
      tmp->SetMarkerStyle(markers[i-1]);
      tmp->SetMarkerSize(2.0);
      tmp->SetLineWidth(2);
      tmp->SetTitle("");
      mg->Add(tmp,"pZ");
      leg->AddEntry(tmp,titles[i-1].c_str(),"p");

    } else if ( i > nhistos ) {

      vdouble e(widths.size(),0.);
      vdouble x(bins); for ( uint j = 0; j < x.size(); ++j ) { x[j] += widths[j]; }
      
      TGraphErrors* tmp = new TGraphErrors( nbins,
					    &x.front(),
					    &m1.front(),
					    &e.front(),
					    &s1.front() );
      tmp->SetMarkerStyle(20);
      tmp->SetMarkerSize(1.5);
      tmp->SetLineWidth(2);
      tmp->SetTitle("");
      mg->Add(tmp,"pZ");
      stringstream ss; ss << "Weighted average of all " << nhistos << " closure tests";
			 leg->AddEntry(tmp,ss.str().c_str(),"p");

    }

  }

  TCanvas* c = new TCanvas("tmp","tmp",900,600);
  mg->Draw("ap");
  mg->GetXaxis()->SetTitle("H_{T} (GeV)");
  mg->GetYaxis()->SetTitle("( N_{obs} - N_{pred} ) / N_{pred}");
  mg->GetYaxis()->SetRangeUser(-1.5,3.5);
  mg->GetXaxis()->SetRangeUser(270.,980.);
  mg->GetXaxis()->SetNdivisions(510);
  leg->Draw("same");
  c->Update();
  
  return;
  
}

// -----------------------------------------------------------------------------
//
void init( uint size, vdouble& v ){ 
  v.clear(); 
  v.resize(size,0.); 
}

void init( uint size1, uint size2, vvdouble& v ){ 
  v.clear(); 
  v.resize(size1,vdouble(size2,0.)); 
}

// -----------------------------------------------------------------------------
//
uint extractData( vdouble& bins, 
		  vvdouble& values, 
		  vvdouble& errors,
		  vstring& titles,
		  bool incl ) {
  
  if ( incl ) {

    uint nhistos = 5;
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);
    
    // Darren's input files
    TFile* file[2];
    file[0] = new TFile("RA1_List2_4_6_5_Spread.root");
    file[1] = new TFile("RA1_List2_3_8_8_Spread.root");
    
    // Extract histogram contents
    for ( uint i = 0; i < 8; ++i ) {
      stringstream ss; ss << "ht_" << bins[i];
      values[0][i] = ((TH1D*)file[1]->Get(ss.str().c_str()))->GetBinContent(1);
      values[1][i] = ((TH1D*)file[1]->Get(ss.str().c_str()))->GetBinContent(2);
      values[2][i] = ((TH1D*)file[1]->Get(ss.str().c_str()))->GetBinContent(3);
      values[3][i] = ((TH1D*)file[0]->Get(ss.str().c_str()))->GetBinContent(1);
      values[4][i] = ((TH1D*)file[0]->Get(ss.str().c_str()))->GetBinContent(2);
      errors[0][i] = ((TH1D*)file[1]->Get(ss.str().c_str()))->GetBinError(1);
      errors[1][i] = ((TH1D*)file[1]->Get(ss.str().c_str()))->GetBinError(2);
      errors[2][i] = ((TH1D*)file[1]->Get(ss.str().c_str()))->GetBinError(3);
      errors[3][i] = ((TH1D*)file[0]->Get(ss.str().c_str()))->GetBinError(1);
      errors[4][i] = ((TH1D*)file[0]->Get(ss.str().c_str()))->GetBinError(2);
    }

    // Missing entries for photon tests
    values[3][0] = 100.;
    values[3][1] = 100.;
    values[4][0] = 100.;
    values[4][1] = 100.;

    titles.push_back("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (using #mu + jets sample)");
    titles.push_back("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (using #mu#mu + jets sample)");
    titles.push_back("0 b-tags #rightarrow 1 b-tag (using #mu + jets sample with no #alpha_{T} cut)");
    titles.push_back("#mu + jets #rightarrow #mu#mu + jets (both samples with no #alpha_{T} cut)");
    titles.push_back("#mu + jets (sample with no #alpha_{T} cut) #rightarrow #gamma + jets");

//     titles.push_back("#mu + jets, #alpha_{T} < 0.55 #rightarrow #mu + jets, #alpha_{T} > 0.55");
//     titles.push_back("#mu#mu + jets, #alpha_{T} < 0.55 #rightarrow #mu#mu + jets, #alpha_{T} > 0.55");
//     titles.push_back("#mu + jets, 0 b-tags, no #alpha_{T} cut #rightarrow #mu + jets, 1 b-tag, no #alpha_{T} cut");
//     titles.push_back("#mu + jets, no #alpha_{T} cut #rightarrow #mu#mu + jets, no #alpha_{T} cut");
//     titles.push_back("#mu + jets, no #alpha_{T} cut #rightarrow #gamma + jets, no #alpha_{T} cut");

    return nhistos;

  } else {

    uint nhistos = 4;
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);
    
    // Darren's input files
    TFile* file = new TFile("RA1_List_3_Spread.root");
    
    // Extract histogram contents
    for ( uint i = 0; i < 8; ++i ) {
      stringstream ss; ss << "ht_" << bins[i];
      values[0][i] = ((TH1D*)file->Get(ss.str().c_str()))->GetBinContent(3);
      values[1][i] = ((TH1D*)file->Get(ss.str().c_str()))->GetBinContent(4);
      values[2][i] = ((TH1D*)file->Get(ss.str().c_str()))->GetBinContent(1);
      values[3][i] = ((TH1D*)file->Get(ss.str().c_str()))->GetBinContent(2);
      errors[0][i] = ((TH1D*)file->Get(ss.str().c_str()))->GetBinError(3);
      errors[1][i] = ((TH1D*)file->Get(ss.str().c_str()))->GetBinError(4);
      errors[2][i] = ((TH1D*)file->Get(ss.str().c_str()))->GetBinError(1);
      errors[3][i] = ((TH1D*)file->Get(ss.str().c_str()))->GetBinError(2);
    }

    titles.push_back("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (using #mu + jets sample with 1 b-tag)");
    titles.push_back("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (using #mu + jets sample with #geq 2 b-tags)");
    titles.push_back("0 b-tags #rightarrow 1 b-tag (using #mu + jets sample with no #alpha_{T} cut)");
    titles.push_back("#mu + jets #rightarrow #mu#mu + jets (samples with 1 b-tag and no #alpha_{T} cut)");

//     titles.push_back("#mu + jets, 1 b-tag, #alpha_{T} < 0.55 #rightarrow #mu + jets, 1 b-tag, #alpha_{T} > 0.55");
//     titles.push_back("#mu#mu + jets, 1 b-tag, #alpha_{T} < 0.55 #rightarrow #mu#mu + jets, 1 b-tag, #alpha_{T} > 0.55");
//     titles.push_back("#mu + jets, 0 b-tags #rightarrow #mu + jets, 1 b-tag");
//     titles.push_back("#mu + jets, 1 b-tag, no #alpha_{T} cut #rightarrow #mu#mu + jets, 1 b-tag, no #alpha_{T} cut");

    return nhistos;

  }
  
}

// -----------------------------------------------------------------------------
//
void extractNumbers( vvdouble& v, vvdouble& e, vdouble& b, vdouble& m, vdouble& s ) {

  // Init
  m.clear(); m.resize(!v.empty()?v[0].size():0,0.);
  s.clear(); s.resize(!v.empty()?v[0].size():0,0.);
  vdouble n; init(b.size(),n);
  vdouble n2; init(b.size(),n2);
  vdouble d; init(b.size(),d);
  
  // Check sizes
  sort(b.begin(),b.end());
  if ( b.empty() || b.back() >= v[0].size() ) { cout << "Problem!" << endl; return; } 
  
  for ( uint i = 0; i < v.size(); ++i ) {
    for ( uint j = 0; j < b.size(); ++j ) {
      for ( uint k = b[j]; k < (j<b.size()-1?b[j+1]:v[i].size()); ++k ) {
	//cout << " test " << i << " " << j << " " << k << endl;
 	if ( e[i][k] > 0. && v[i][k] < 100. ) {
	  double w = weighted ? 1./(e[i][k]*e[i][k]) : 1.;
 	  n[j] += v[i][k] * w;
 	  d[j] += w;
 	  n2[j] += v[i][k]*v[i][k];
	}
      }
    }
  }
  for ( uint j = 0; j < b.size(); ++j ) {
    if ( d[j] > 0 ) {
      m[j] = n[j] / d[j];
      if (weighted) { s[j] = 1. / TMath::Sqrt(d[j]); } 
      else { s[j] = TMath::Sqrt(n2[j]/d[j]); }
    }
  }
  
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


