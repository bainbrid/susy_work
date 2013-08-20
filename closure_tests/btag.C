#include "TMultiGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TStyle.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveText.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

typedef unsigned int uint;
typedef vector<double> vdouble;
typedef vector<string> vstring;
typedef vector<vdouble> vvdouble;

// options
/*
  0=incl, 
  1=excl, 
  2=btag, 

  3=paper, 
  4=7tevmc, 
  5=7tevmclatest, 

  6=mumu-gamma-only, 

  7=8tev-pre-app, 
  8=8tev-updated, 

  9=8tev,kine-only, 
  10=8tev,alphat
*/

int set_option = 3; 
int syst_option = 0; // 0=weighted,3sigma, 1=unweighted,|mean+std.dev.|, 2=68%CL
bool weighted = true;
bool draw_weighted = false;
bool detailed_syst = false;
bool detailed_fit = false;
bool fixed_syst_bands = true;
std::string lumi = "4.98";
std::string energy = "7";
double confidence_level = 0.90;

vdouble bins;
vdouble widths;
vdouble syst; 
vdouble markers;

void setTDRStyle();
void init( uint, vdouble& );
void init( uint, uint, vvdouble& );
int extractData( vvdouble&, vvdouble&, vstring& );
void extractNumbers( vvdouble&, vvdouble&, vdouble&, vdouble&, vdouble&, vdouble& );
void fitResults( TGraphErrors*, double&, double&, double&, int& );
double fabs( double x ) { return x<0.?x*-1.:x; }

typedef std::pair<int,double> Pair;
typedef std::vector<Pair> Pairs;
bool ascending( const Pair& x, const Pair& y ) { return  x.second < y.second; }

void btag() {
   
  setTDRStyle();

  // Bins and widths
  bins.push_back(275.);
  bins.push_back(325.);
  for ( uint i = 0; i < 6; ++i ) { bins.push_back(375.+100.*i); }
  bins.push_back(975.);
  
  widths.push_back(25.);
  widths.push_back(25.);
  for ( uint i = 0; i < 6; ++i ) { widths.push_back(50.); }

  markers.push_back(24);
  markers.push_back(25);
  markers.push_back(26);
  markers.push_back(28);
  markers.push_back(30);
  markers.push_back(32);
  
  // Extract data from histograms
  vvdouble values; 
  vvdouble errors;
  vstring titles;
  const int nhistos = extractData( values, errors, titles );

  // Extract central values, errors, and bands for each HT bin
  vdouble n1; for ( uint i = 0; i < bins.size()-1; ++i ) { n1.push_back(i); }
  vdouble m1; 
  vdouble s1;
  vdouble u1; 
  extractNumbers( values, errors, n1, m1, s1, u1 ); 

  // Extract central values, errors, and bands for each HT region
  vdouble n2; n2.push_back(0); n2.push_back(4); n2.push_back(6);
  vdouble m2; 
  vdouble s2;
  vdouble u2; 
  extractNumbers( values, errors, n2, m2, s2, u2 ); 
  
  vdouble m; 
  m.resize(m2.size(),0.);

  TLegend* leg = new TLegend( 0.16, 0.85-0.045*(nhistos+(draw_weighted?2:1)), 0.35, 0.85 );
  leg->SetFillColor(0);
  leg->SetLineColor(0); 
  leg->SetShadowColor(0); 
  leg->SetTextSize(0.035);

  TMultiGraph* mg = new TMultiGraph();
  int offset = 0;
  for ( int i = 0; i < nhistos+2; ++i ) { 
    offset++;

    if ( i == 0 ) { // Systematic bands for HT regions
      
      vdouble x; for ( uint j = 0; j < n2.size(); ++j ) { x.push_back( bins[n2[j]] ); } 
      vdouble xel(n2.size(),0.);
      vdouble xeh; for ( uint j = 0; j < n2.size(); ++j ) { xeh.push_back( (j<n2.size()-1)?x[j+1]-x[j]:bins[bins.size()-1]-x[j] ); }
      
      TGraphAsymmErrors* tmp = new TGraphAsymmErrors( n2.size(),
						      &x.front(),
						      &m.front(),
						      &xel.front(),
						      &xeh.front(),
						      &u2.front(),
						      &u2.front() );
      tmp->SetMarkerStyle(20);
      tmp->SetMarkerSize(1.5);
      tmp->SetFillColor(kGray);
      tmp->SetTitle("");
      mg->Add(tmp,"2");

      stringstream ss; 
      if ( detailed_syst ) { 
 	ss << "Systematic uncertainty ("
 	   << int(100.*s2[0]) << "%, "
 	   << int(100.*s2[1]) << "%, and "
 	   << int(100.*s2[2]) << "%)";
      } else {
	ss << "Systematic uncertainty";
      }
      leg->AddEntry(tmp,ss.str().c_str(),"f");

    } else if ( i > 0 && i <= nhistos ) { // All individual closure tests
      
      //continue;
      
      vdouble e(widths.size(),0.);
      vdouble x(bins); for ( uint j = 0; j < x.size(); ++j ) { x[j] += widths[j] + 5.*(offset%2?1*offset/2:-1*offset/2); }
      
      TGraphErrors* rob1 = new TGraphErrors( bins.size()-1,
					     &x.front(),
					     &values[i-1].front(),
					     &e.front(),
					     &errors[i-1].front() );
      rob1->SetMarkerStyle(markers[i-1]);
      rob1->SetMarkerSize(2.0);
      rob1->SetLineWidth(2);
      rob1->SetTitle("");
      mg->Add(rob1,"pZ");

      std::cout << "Plot #: " << i-1 << " Nbins: " << bins.size() << " Title: " << titles[i-1] << std::endl;
      for ( uint ii = 0; ii < bins.size()-1; ++ii ) { 
	std::cout << std::fixed
		  << std::setprecision(0)
		  << " HT: " << x[ii]
	  //<< " width: " << e[ii]
		  << std::setprecision(2)
		  << " val: " << values[i-1][ii]
		  << " err: " << errors[i-1][ii]
		  << std::endl;
      }

      double fmean,ferr,fchi2; int fdof;
      fitResults(rob1,fmean,ferr,fchi2,fdof);
      
      stringstream ss; 
      ss << titles[i-1];
      if ( detailed_fit ) { 
	ss  << " (A = " << std::setprecision(2) << fmean
	    << "#pm"  << std::setprecision(2) << ferr
	    << ", #chi^{2}/dof = "  << std::setprecision(1) << fchi2
	    << "/" << std::setprecision(0) << fdof
	    << ")";
      }
      leg->AddEntry(rob1,ss.str().c_str(),"p");

    } else if ( draw_weighted && i > nhistos ) { // Bin-averaged result

      vdouble e(widths.size(),0.);
      vdouble x; x.reserve(bins.size()-1); 
      std::copy( bins.begin(), bins.end()-1, back_inserter(x) ); 
      for ( uint j = 0; j < x.size(); ++j ) { x[j] += widths[j]; }
      
      TGraphErrors* tmp2 = new TGraphErrors( bins.size(),
					     &x.front(),
					     &m1.front(),
					     &e.front(),
					     &s1.front() );
      tmp2->SetMarkerStyle(20);
      tmp2->SetMarkerSize(1.5);
      tmp2->SetLineWidth(2);
      tmp2->SetTitle("");
      mg->Add(tmp2,"pZ");
      stringstream ss; ss << "Weighted average of all " << nhistos << " closure tests";
      leg->AddEntry(tmp2,ss.str().c_str(),"p");

      for ( uint ii = 0; ii < m1.size(); ++ii ) { 
	std::cout << std::fixed
		  << std::setprecision(0)
		  << " HT: " << bins[ii]
		  << std::setprecision(2)
		  << " mean: " << m1[ii]
		  << " std. dev.: " << s1[ii]/3.
		  << std::endl;
      }

      double fmean,ferr,fchi2; int fdof;
      fitResults(tmp2,fmean,ferr,fchi2,fdof);
      
    }

  }

  TCanvas* c = new TCanvas("tmp","tmp",900,600);
  mg->Draw("ap");
  mg->GetXaxis()->SetTitle("H_{T} (GeV)");
  mg->GetYaxis()->SetTitle("( N_{obs} - N_{pred} ) / N_{pred}");
  //mg->GetYaxis()->SetRangeUser(-1.5,2.5);
  mg->GetYaxis()->SetRangeUser(-2.,4.);
  mg->GetXaxis()->SetRangeUser(270.,980.);
  mg->GetXaxis()->SetNdivisions(510);
  leg->Draw("same");
  std::stringstream sss;
  sss << "CMS, L_{int} = " << lumi << " fb^{-1}, #sqrt{s} = " << energy << " TeV";
  TLatex* tex = new TLatex(0.17,0.88,sss.str().c_str());
  tex->SetNDC();
  tex->SetTextSize(0.035);
  tex->Draw();
  c->Update();
  
  return;
  
}

// -----------------------------------------------------------------------------
//
void init( uint size, vdouble& v ) { 
  v.clear(); 
  v.resize(size,0.); 
}

void init( uint size1, uint size2, vvdouble& v ) { 
  v.clear(); 
  v.resize(size1,vdouble(size2,0.)); 
}

// -----------------------------------------------------------------------------
//
void extractNumbers( vvdouble& v, 
		     vvdouble& e, 
		     vdouble& n, 
		     vdouble& m, 
		     vdouble& s,
		     vdouble& u ) {

  // Init
  m.clear(); m.resize(!v.empty()?v[0].size():0,0.);
  s.clear(); s.resize(!v.empty()?v[0].size():0,0.);
  u.clear(); u.resize(!v.empty()?v[0].size():0,0.);

  vdouble numer; init(n.size(),numer);
  vdouble numer2; init(n.size(),numer2);
  vdouble denom; init(n.size(),denom);
  
  // Check sizes
  sort(n.begin(),n.end());
  if ( n.empty() || n.back() >= v[0].size() ) { cout << "Problem!" << endl; return; } 

  // Increment counters
  //   vvdouble means;
  //   vvdouble ;
  for ( uint i = 0; i < v.size(); ++i ) {
    for ( uint j = 0; j < n.size(); ++j ) {
      for ( uint k = n[j]; k < (j<n.size()-1?n[j+1]:v[i].size()); ++k ) {
	//cout << " test " << i << " " << j << " " << k << endl;
 	if ( e[i][k] > 0. && v[i][k] < 100. ) {
	  double w = weighted ? 1./(e[i][k]*e[i][k]) : 1.;
 	  numer[j] += v[i][k] * w;
 	  denom[j] += w;
 	  numer2[j] += v[i][k]*v[i][k];
	}
      }
    }
  }
  
  // Weighted mean/variance or not?
  if ( syst_option == 1 ) { weighted = false; }
  else { weighted = true; }
  
  // Calc mean and error
  for ( uint j = 0; j < n.size(); ++j ) {
    if ( denom[j] > 0 ) {
      m[j] = numer[j] / denom[j];
      if (weighted) { s[j] = TMath::Sqrt( 1. / denom[j] ); } 
      //else { s[j] = TMath::Sqrt( numer2[j]/denom[j] ); } // r.m.s.
      else { s[j] = TMath::Sqrt( numer2[j]/denom[j] - m[j]*m[j] ); } // std. dev.
    } else { m[j] = -1000.; s[j] = 0.; }
  }  

  // Calc bands
  if ( syst_option == 0 ) {

    // weighted, 3.*std.dev.

    u.resize(s.size(),0.);
    for ( uint i = 0; i < s.size(); ++i ) { s[i] *= 3.; u[i] = s[i]; }

  } else if ( syst_option == 1 ) {

    // unweighted, | mean + std.dev.|

    u.resize(s.size(),0.);
    for ( uint i = 0; i < s.size(); ++i ) { u[i] = fabs( m[i]>0.?m[i]+s[i]:m[i]-s[i] ); }
    
  } else if ( syst_option == 2 ) {

    // calc 68% CL

    s.clear();
    s.resize(n.size(),0.);

    u.clear();
    u.resize(n.size(),0.);

    // Iterate through HT regions
    for ( uint j = 0; j < n.size(); ++j ) {
      vdouble vals;
      // Iterate through sets of closure tests
      for ( uint i = 0; i < v.size(); ++i ) {
	// Iterate through HT bins in a given region
	for ( uint k = n[j]; k < (j<n.size()-1?n[j+1]:v[i].size()); ++k ) { 
	  if ( fabs(v[i][k]) > 999. ) { continue; }
	  vals.push_back(fabs(v[i][k])); 
	  std::cout << " j: " << j 
		    << " i: " << i
		    << " k: " << k
		    << " v: " << v[i][k]
		    << " vals: " << vals.back()
		    << std::endl;
	}
      }
      std::sort( vals.begin(), vals.end() );
      for ( uint i = 0; i < vals.size(); ++i ) {
	std::cout << " i: " << i
		  << " vals: " << vals[i]
		  << std::endl;
      }
      //int index = int(confidence_level*vals.size()) < vals.size() ? int(confidence_level*vals.size()) : vals.size()-1;
      int index = vals.size()-1;
      uint temp = int( confidence_level * vals.size() );
      if ( temp < vals.size() ) {
	index = temp;
	double tmp3 = confidence_level * vals.size() - double(index);
	if ( tmp3 > 0. )  index++;
	index--;
      }
      s[j] = vals[index];
      u[j] = vals[index];
      std::cout << " j: " << j
		<< " size: " << vals.size()
		<< " idx: " << confidence_level*vals.size()
		<< " index: " << index
		<< " val: " << vals[index]
		<< std::endl;
    }
    
    //       Pairs sig;
    //       for ( int ii = 0; ii < means.size(); ++ii ) { sig.push_back( std::make_pair(ii,s[ii]>0.?fabs(m[ii]/s[ii]):0.) ); }
    //       std::sort( sig.begin(), sig.end(), ascending );
      
    //       int kk = 0;
    //       Pairs::const_iterator ii = sig.begin();
    //       Pairs::const_iterator jj = sig.end();
    //       for ( ; ii != jj; ++ii  ) { 
    // 	std::cout << " index " << kk
    // 		  << " was " << ii->first
    // 		  << " sig " << ii->second
    // 		  << " mean " << m[ii->first]
    // 		  << " err " << s[ii->first]
    // 		  << " " << std::endl;
    // 	kk++;
    //       }
      
    //       int index = int ( confidence_level * m.size() ) + 1;
    //       if ( index >= m.size() ) { index = m.size()-1; }

    
  } else {
    std::cout << "UNKNOWN OPTION!!!" << std::endl;
  }
  
  for ( uint j = 0; j < n.size(); ++j ) {
    std::cout << " bin: " << j
	      << " HT: " << bins[n[j]]
	      << " sum: " << numer[j]
	      << " sum2: " << numer2[j]
	      << " n: " << denom[j]
	      << " mean: " << m[j]
	      << " err: " << s[j]
	      << " uncert: " << u[j]
	      << std::endl;
  }
  

  // Use fixed systematic bands
  if ( fixed_syst_bands ) {
    syst.resize(u.size(),0.);
    u.clear();
    for ( uint i = 0; i < syst.size(); ++i ) { u.push_back(syst[i]); }
  }

}  
// -----------------------------------------------------------------------------
//
void fitResults( TGraphErrors* tmp, double& tmp1, double& tmp2, double& tmp3, int& tmp4 ) {

  TGraphErrors* temp = new TGraphErrors( *tmp );
  TF1* fit = 0;


  double x = 0., y = 0.; (((TGraphErrors*)(temp))->GetPoint(0,x,y));
  if ( x < 275. ) { fit = new TF1("fit","pol0",370.,980.); }
  else            { fit = new TF1("fit","pol0",270.,980.); }
  temp->Fit(fit,"QR");

  //   double a275 = fit->GetParameter(0) + 275.*fit->GetParameter(1);
  //   double e275 = fit->GetParameter(0) + 275.*fit->GetParameter(1);
  
  std::cout << " name & " << std::setprecision(2) << fit->GetParameter(0)
	    << " & " << std::setprecision(2) << fit->GetParError(0)
	    << " & " << std::setprecision(1) << fabs(fit->GetParameter(0)/fit->GetParError(0))
    // 	    << " & " << std::setprecision(2) << fit->GetParameter(1)
    // 	    << " & " << std::setprecision(2) << fit->GetParError(1)
	    << " & " << std::setprecision(2) << fit->GetChisquare() 
	    << " & " << std::setprecision(0) << fit->GetNDF() 
	    << " & " << std::setprecision(2) << fit->GetProb() 
	    << " \\\\ "
	    << std::endl;
  
  //   std::cout << std::setprecision(2)
  // 	    << " mean: " << fit->GetParameter(0)
  // 	    << " error: " << fit->GetParError(0)
  // 	    << " chi2: " << fit->GetChisquare() 
  // 	    << std::setprecision(0)
  // 	    << " ndof: " << fit->GetNDF() 
  // 	    << std::setprecision(2)
  // 	    << " p-value: " << fit->GetProb() 
  // 	    << std::endl;
  
  tmp1 = fit->GetParameter(0);
  tmp2 = fit->GetParError(0);
  tmp3 = fit->GetChisquare();
  tmp4 = fit->GetNDF();
  
}

// -----------------------------------------------------------------------------

void setTDRStyle() {

  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

  //   tdrStyle->SetOptFit(kFALSE);
  //   tdrStyle->SetOptStat(kFALSE);

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

// -----------------------------------------------------------------------------
//
int extractData( vvdouble& values, 
		 vvdouble& errors,
		 vstring& titles ) {
  
  if ( set_option == 0 ) {

    int nhistos = 5;
    
    syst.clear(); 
    syst.push_back(0.10);
    syst.push_back(0.20);
    syst.push_back(0.40);

    init(nhistos,bins.size()-1,values);
    init(nhistos,bins.size()-1,errors);
    
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
    titles.push_back("0 b tags #rightarrow 1 b tag (using #mu + jets sample with no #alpha_{T} cut)");
    titles.push_back("#mu + jets #rightarrow #mu#mu + jets (both samples with no #alpha_{T} cut)");
    titles.push_back("#mu + jets (sample with no #alpha_{T} cut) #rightarrow #gamma + jets");

    //     titles.push_back("#mu + jets, #alpha_{T} < 0.55 #rightarrow #mu + jets, #alpha_{T} > 0.55");
    //     titles.push_back("#mu#mu + jets, #alpha_{T} < 0.55 #rightarrow #mu#mu + jets, #alpha_{T} > 0.55");
    //     titles.push_back("#mu + jets, 0 b tags, no #alpha_{T} cut #rightarrow #mu + jets, 1 b tag, no #alpha_{T} cut");
    //     titles.push_back("#mu + jets, no #alpha_{T} cut #rightarrow #mu#mu + jets, no #alpha_{T} cut");
    //     titles.push_back("#mu + jets, no #alpha_{T} cut #rightarrow #gamma + jets, no #alpha_{T} cut");

    return nhistos;

    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
  } else if ( set_option == 1 ) {

    int nhistos = 4;
    
    syst.clear();
    syst.push_back(0.10);
    syst.push_back(0.20);
    syst.push_back(0.40);
    
    init(nhistos,bins.size()-1,values);
    init(nhistos,bins.size()-1,errors);
    
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

    titles.push_back("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (using #mu + jets sample with 1 b tag)");
    titles.push_back("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (using #mu + jets sample with #geq 2 b tags)");
    titles.push_back("0 b tags #rightarrow 1 b tag (using #mu + jets sample with no #alpha_{T} cut)");
    titles.push_back("#mu + jets #rightarrow #mu#mu + jets (samples with 1 b tag and no #alpha_{T} cut)");

    //     titles.push_back("#mu + jets, 1 b tag, #alpha_{T} < 0.55 #rightarrow #mu + jets, 1 b tag, #alpha_{T} > 0.55");
    //     titles.push_back("#mu#mu + jets, 1 b tag, #alpha_{T} < 0.55 #rightarrow #mu#mu + jets, 1 b tag, #alpha_{T} > 0.55");
    //     titles.push_back("#mu + jets, 0 b tags #rightarrow #mu + jets, 1 b tag");
    //     titles.push_back("#mu + jets, 1 b tag, no #alpha_{T} cut #rightarrow #mu#mu + jets, 1 b tag, no #alpha_{T} cut");

    return nhistos;

    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
  } else if ( set_option == 2 ) {

    int nhistos = 4;
    
    syst.clear();
    syst.push_back(0.10);
    syst.push_back(0.20);
    syst.push_back(0.40);

    init(nhistos,bins.size()-1,values);
    init(nhistos,bins.size()-1,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu, 0->1, no aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (using #mu + jets sample with no #alpha_{T} cut)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.07545189);
    grae->SetPointError(1,0,0,0.2161792,0.2161792);
    grae->SetPoint(2,350,0.02341433);
    grae->SetPointError(2,0,0,0.3263656,0.3263656);
    grae->SetPoint(3,425,0.05799095);
    grae->SetPointError(3,0,0,0.05314723,0.05314723);
    grae->SetPoint(4,525,0.1261303);
    grae->SetPointError(4,0,0,0.07825716,0.07825716);
    grae->SetPoint(5,625,0.02131099);
    grae->SetPointError(5,0,0,0.1180257,0.1180257);
    grae->SetPoint(6,725,-0.05882136);
    grae->SetPointError(6,0,0,0.17784,0.17784);
    grae->SetPoint(7,825,-0.1930679);
    grae->SetPointError(7,0,0,0.2431383,0.2431383);
    grae->SetPoint(8,925,0.04710806);
    grae->SetPointError(8,0,0,0.2413123,0.2413123);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, 1->2, no aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (using #mu + jets sample with no #alpha_{T} cut)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,-0.09088821);
    grae->SetPointError(1,0,0,0.2529462,0.2529462);
    grae->SetPoint(2,350,-0.09028696);
    grae->SetPointError(2,0,0,0.3757164,0.3757164);
    grae->SetPoint(3,425,0.03644863);
    grae->SetPointError(3,0,0,0.08315404,0.08315404);
    grae->SetPoint(4,525,0.04415105);
    grae->SetPointError(4,0,0,0.112311,0.112311);
    grae->SetPoint(5,625,0.05825824);
    grae->SetPointError(5,0,0,0.1792305,0.1792305);
    grae->SetPoint(6,725,0.131543);
    grae->SetPointError(6,0,0,0.278603,0.278603);
    grae->SetPoint(7,825,-0.02169071);
    grae->SetPointError(7,0,0,0.4068863,0.4068863);
    grae->SetPoint(8,925,0.116807);
    grae->SetPointError(8,0,0,0.3787215,0.3787215);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, 0 no aT -> 1 w/ aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags, no #alpha_{T} cut #rightarrow 1 b tag, with #alpha_{T} cut (using #mu + jets sample)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.07545189);
    grae->SetPointError(1,0,0,0.166154,0.166154);
    grae->SetPoint(2,350,0.02341433);
    grae->SetPointError(2,0,0,0.256357,0.256357);
    grae->SetPoint(3,425,0.06873148);
    grae->SetPointError(3,0,0,0.1187567,0.1187567);
    grae->SetPoint(4,525,0.212922);
    grae->SetPointError(4,0,0,0.201002,0.201002);
    grae->SetPoint(5,625,-0.6010028);
    grae->SetPointError(5,0,0,0.3623177,0.3623177);
    grae->SetPoint(6,725,0.1037815);
    grae->SetPointError(6,0,0,0.5709566,0.5709566);
    grae->SetPoint(7,825,0.4608314);
    grae->SetPointError(7,0,0,1.039041,1.039041);
    grae->SetPoint(8,925,0.8461076);
    grae->SetPointError(8,0,0,1.866444,1.866444);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, 1 no aT -> 2 w/ aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag, no #alpha_{T} cut #rightarrow 2 b tags, with #alpha_{T} cut (using #mu + jets sample)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,-0.09088821);
    grae->SetPointError(1,0,0,0.196108,0.196108);
    grae->SetPoint(2,350,-0.09028696);
    grae->SetPointError(2,0,0,0.2923519,0.2923519);
    grae->SetPoint(3,425,0.144123);
    grae->SetPointError(3,0,0,0.2028936,0.2028936);
    grae->SetPoint(4,525,0.135388);
    grae->SetPointError(4,0,0,0.3097176,0.3097176);
    grae->SetPoint(5,625,0.2220949);
    grae->SetPointError(5,0,0,0.4833377,0.4833377);
    grae->SetPoint(6,725,1.026469);
    grae->SetPointError(6,0,0,1.291714,1.291714);
    grae->SetPoint(7,825,-1);
    grae->SetPointError(7,0,0,2.896569,2.896569);
    grae->SetPoint(8,925,-1);
    grae->SetPointError(8,0,0,2.220931,2.220931);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    return nhistos;

    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
  } else if ( set_option == 3 ) {

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.10);
    syst.push_back(0.20);
    syst.push_back(0.40);

    init(nhistos,bins.size()-1,values);
    init(nhistos,bins.size()-1,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,-0.02592075);
    grae->SetPointError(1,0,0,0.06408704,0.06408704);
    grae->SetPoint(2,350,-0.01530341);
    grae->SetPointError(2,0,0,0.09168221,0.09168221);
    grae->SetPoint(3,425,0.07170955);
    grae->SetPointError(3,0,0,0.06258445,0.06258445);
    grae->SetPoint(4,525,-0.0682347);
    grae->SetPointError(4,0,0,0.1029296,0.1029296);
    grae->SetPoint(5,625,-0.2259311);
    grae->SetPointError(5,0,0,0.1757231,0.1757231);
    grae->SetPoint(6,725,0.3348285);
    grae->SetPointError(6,0,0,0.3169539,0.3169539);
    grae->SetPoint(7,825,0.2588858);
    grae->SetPointError(7,0,0,0.5014861,0.5014861);
    grae->SetPoint(8,925,-0.701126);
    grae->SetPointError(8,0,0,0.7159149,0.7159149);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.07519879);
    grae->SetPointError(1,0,0,0.2170091,0.2170091);
    grae->SetPoint(2,350,0.02346486);
    grae->SetPointError(2,0,0,0.3281894,0.3281894);
    grae->SetPoint(3,425,0.05780617);
    grae->SetPointError(3,0,0,0.05482093,0.05482093);
    grae->SetPoint(4,525,0.1261081);
    grae->SetPointError(4,0,0,0.07941886,0.07941886);
    grae->SetPoint(5,625,0.02135401);
    grae->SetPointError(5,0,0,0.1187902,0.1187902);
    grae->SetPoint(6,725,-0.05868353);
    grae->SetPointError(6,0,0,0.1783514,0.1783514);
    grae->SetPoint(7,825,-0.1931447);
    grae->SetPointError(7,0,0,0.243516,0.243516);
    grae->SetPoint(8,925,0.04742296);
    grae->SetPointError(8,0,0,0.2417228,0.2417228);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,-0.09107025);
    grae->SetPointError(1,0,0,0.2536254,0.2536254);
    grae->SetPoint(2,350,-0.09040525);
    grae->SetPointError(2,0,0,0.3773004,0.3773004);
    grae->SetPoint(3,425,0.03650087);
    grae->SetPointError(3,0,0,0.08423514,0.08423514);
    grae->SetPoint(4,525,0.04394681);
    grae->SetPointError(4,0,0,0.1131024,0.1131024);
    grae->SetPoint(5,625,0.05812638);
    grae->SetPointError(5,0,0,0.1797252,0.1797252);
    grae->SetPoint(6,725,0.1316546);
    grae->SetPointError(6,0,0,0.2789421,0.2789421);
    grae->SetPoint(7,825,-0.02185926);
    grae->SetPointError(7,0,0,0.4070769,0.4070769);
    grae->SetPoint(8,925,0.1164298);
    grae->SetPointError(8,0,0,0.3789139,0.3789139);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.05145438);
    grae->SetPointError(1,0,0,0.1312231,0.1312231);
    grae->SetPoint(2,350,-0.1429071);
    grae->SetPointError(2,0,0,0.1868665,0.1868665);
    grae->SetPoint(3,425,-0.0417508);
    grae->SetPointError(3,0,0,0.08235534,0.08235534);
    grae->SetPoint(4,525,-0.1234049);
    grae->SetPointError(4,0,0,0.1219376,0.1219376);
    grae->SetPoint(5,625,0.223824);
    grae->SetPointError(5,0,0,0.1895285,0.1895285);
    grae->SetPoint(6,725,0.4025885);
    grae->SetPointError(6,0,0,0.2932544,0.2932544);
    grae->SetPoint(7,825,-0.524625);
    grae->SetPointError(7,0,0,0.4604737,0.4604737);
    grae->SetPoint(8,925,0.6670811);
    grae->SetPointError(8,0,0,0.4608829,0.4608829);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#mu#mu + jets #rightarrow #gamma + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,0,0);
    grae->SetPointError(1,0,0,0,0);
    grae->SetPoint(2,0,0);
    grae->SetPointError(2,0,0,0,0);
    grae->SetPoint(3,425,-0.1411698);
    grae->SetPointError(3,0,0,0.08648166,0.08648166);
    grae->SetPoint(4,525,-0.02322981);
    grae->SetPointError(4,0,0,0.1315553,0.1315553);
    grae->SetPoint(5,625,-0.3074706);
    grae->SetPointError(5,0,0,0.1873511,0.1873511);
    grae->SetPoint(6,725,-0.1557285);
    grae->SetPointError(6,0,0,0.2854318,0.2854318);
    grae->SetPoint(7,825,1.073652);
    grae->SetPointError(7,0,0,0.9323202,0.9323202);
    grae->SetPoint(8,925,-0.276403);
    grae->SetPointError(8,0,0,0.4515533,0.4515533);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    return nhistos;

    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
  } else if ( set_option == 4 ) {

    energy = "8";
    lumi = "0.56";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.20);
    syst.push_back(0.40);
    syst.push_back(0.60);

    init(nhistos,bins.size()-1,values);
    init(nhistos,bins.size()-1,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.06051738);
    grae->SetPointError(1,0,0,0.09525032,0.09525032);
    grae->SetPoint(2,350,0.04926848);
    grae->SetPointError(2,0,0,0.1353258,0.1353258);
    grae->SetPoint(3,425,-0.03563764);
    grae->SetPointError(3,0,0,0.13925,0.13925);
    grae->SetPoint(4,525,0.1102946);
    grae->SetPointError(4,0,0,0.247521,0.247521);
    grae->SetPoint(5,625,-0.2204474);
    grae->SetPointError(5,0,0,0.4310368,0.4310368);
    grae->SetPoint(6,725,0.1303485);
    grae->SetPointError(6,0,0,0.647188,0.647188);
    grae->SetPoint(7,825,0.2483287);
    grae->SetPointError(7,0,0,1.17989,1.17989);
    grae->SetPoint(8,925,-1);
    grae->SetPointError(8,0,0,1.863572,1.863572);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.08590583);
    grae->SetPointError(1,0,0,0.08300547,0.08300547);
    grae->SetPoint(2,350,0.1482044);
    grae->SetPointError(2,0,0,0.1121917,0.1121917);
    grae->SetPoint(3,425,0.01777188);
    grae->SetPointError(3,0,0,0.107452,0.107452);
    grae->SetPoint(4,525,-0.06281969);
    grae->SetPointError(4,0,0,0.157514,0.157514);
    grae->SetPoint(5,625,-0.3042731);
    grae->SetPointError(5,0,0,0.2648389,0.2648389);
    grae->SetPoint(6,725,0.1305306);
    grae->SetPointError(6,0,0,0.3787774,0.3787774);
    grae->SetPoint(7,825,-0.03162577);
    grae->SetPointError(7,0,0,0.5261531,0.5261531);
    grae->SetPoint(8,925,-0.3932323);
    grae->SetPointError(8,0,0,0.4961877,0.4961877);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.07775372);
    grae->SetPointError(1,0,0,0.1268628,0.1268628);
    grae->SetPoint(2,350,0.2059022);
    grae->SetPointError(2,0,0,0.1641371,0.1641371);
    grae->SetPoint(3,425,-0.1083611);
    grae->SetPointError(3,0,0,0.1719001,0.1719001);
    grae->SetPoint(4,525,0.020751);
    grae->SetPointError(4,0,0,0.2465672,0.2465672);
    grae->SetPoint(5,625,0.2107079);
    grae->SetPointError(5,0,0,0.4596019,0.4596019);
    grae->SetPoint(6,725,0.7616245);
    grae->SetPointError(6,0,0,0.6821133,0.6821133);
    grae->SetPoint(7,825,-0.09948001);
    grae->SetPointError(7,0,0,0.822696,0.822696);
    grae->SetPoint(8,925,1.715365);
    grae->SetPointError(8,0,0,1.935599,1.935599);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.05858088);
    grae->SetPointError(1,0,0,0.1218501,0.1218501);
    grae->SetPoint(2,350,0.0577255);
    grae->SetPointError(2,0,0,0.1633545,0.1633545);
    grae->SetPoint(3,425,-0.10601);
    grae->SetPointError(3,0,0,0.1667387,0.1667387);
    grae->SetPoint(4,525,-0.227977);
    grae->SetPointError(4,0,0,0.2652359,0.2652359);
    grae->SetPoint(5,625,1.120917);
    grae->SetPointError(5,0,0,0.6269889,0.6269889);
    grae->SetPoint(6,725,-0.476071);
    grae->SetPointError(6,0,0,0.6069773,0.6069773);
    grae->SetPoint(7,825,1.080867);
    grae->SetPointError(7,0,0,1.153193,1.153193);
    grae->SetPoint(8,925,-1);
    grae->SetPointError(8,0,0,1.085947,1.085947);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#mu#mu + jets #rightarrow #gamma + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,0,0);
    grae->SetPointError(1,0,0,0,0);
    grae->SetPoint(2,0,0);
    grae->SetPointError(2,0,0,0,0);
    grae->SetPoint(3,425,-0.1455919);
    grae->SetPointError(3,0,0,0.1806481,0.1806481);
    grae->SetPoint(4,525,-0.3547269);
    grae->SetPointError(4,0,0,0.2846392,0.2846392);
    grae->SetPoint(5,625,0.7052895);
    grae->SetPointError(5,0,0,0.5310966,0.5310966);
    grae->SetPoint(6,725,-0.45608);
    grae->SetPointError(6,0,0,0.7556384,0.7556384);
    grae->SetPoint(7,825,3.042455);
    grae->SetPointError(7,0,0,4.659075,4.659075);
    grae->SetPoint(8,925,-1);
    grae->SetPointError(8,0,0,0,0);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    return nhistos;

    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
  } else if ( set_option == 5 ) {

    // http://www.hep.ph.ic.ac.uk/~db1110/22May_801pb_Documents/

    energy = "8";
    lumi = "0.80";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.20);
    syst.push_back(0.40);
    syst.push_back(0.60);

    init(nhistos,bins.size()-1,values);
    init(nhistos,bins.size()-1,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,-0.01615104);
    grae->SetPointError(1,0,0,0.07714497,0.07714497);
    grae->SetPoint(2,350,0.03785635);
    grae->SetPointError(2,0,0,0.1040441,0.1040441);
    grae->SetPoint(3,425,-0.0259905);
    grae->SetPointError(3,0,0,0.1028576,0.1028576);
    grae->SetPoint(4,525,0.1719211);
    grae->SetPointError(4,0,0,0.1738626,0.1738626);
    grae->SetPoint(5,625,-0.1925034);
    grae->SetPointError(5,0,0,0.2868012,0.2868012);
    grae->SetPoint(6,725,0.3034191);
    grae->SetPointError(6,0,0,0.4517342,0.4517342);
    grae->SetPoint(7,825,-0.470608);
    grae->SetPointError(7,0,0,0.8294227,0.8294227);
    grae->SetPoint(8,925,-1);
    grae->SetPointError(8,0,0,1.276578,1.276578);
       
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,-0.02706759);
    grae->SetPointError(1,0,0,0.06895844,0.06895844);
    grae->SetPoint(2,350,-0.05126671);
    grae->SetPointError(2,0,0,0.08654936,0.08654936);
    grae->SetPoint(3,425,-0.1244027);
    grae->SetPointError(3,0,0,0.08234624,0.08234624);
    grae->SetPoint(4,525,-0.1470136);
    grae->SetPointError(4,0,0,0.1159392,0.1159392);
    grae->SetPoint(5,625,-0.4010268);
    grae->SetPointError(5,0,0,0.1817466,0.1817466);
    grae->SetPoint(6,725,0.02873572);
    grae->SetPointError(6,0,0,0.2516495,0.2516495);
    grae->SetPoint(7,825,-0.004536864);
    grae->SetPointError(7,0,0,0.3523288,0.3523288);
    grae->SetPoint(8,925,-0.4306654);
    grae->SetPointError(8,0,0,0.3482062,0.3482062);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.09245505);
    grae->SetPointError(1,0,0,0.0927446,0.0927446);
    grae->SetPoint(2,350,0.1241869);
    grae->SetPointError(2,0,0,0.1201659,0.1201659);
    grae->SetPoint(3,425,-0.04148397);
    grae->SetPointError(3,0,0,0.1240865,0.1240865);
    grae->SetPoint(4,525,-0.2374894);
    grae->SetPointError(4,0,0,0.1785591,0.1785591);
    grae->SetPoint(5,625,0.1602992);
    grae->SetPointError(5,0,0,0.3026645,0.3026645);
    grae->SetPoint(6,725,0.3835603);
    grae->SetPointError(6,0,0,0.3778018,0.3778018);
    grae->SetPoint(7,825,0.3158481);
    grae->SetPointError(7,0,0,0.5401025,0.5401025);
    grae->SetPoint(8,925,0.9641363);
    grae->SetPointError(8,0,0,0.897725,0.897725);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.1298237);
    grae->SetPointError(1,0,0,0.09511788,0.09511788);
    grae->SetPoint(2,350,0.2364315);
    grae->SetPointError(2,0,0,0.1268917,0.1268917);
    grae->SetPoint(3,425,0.1075587);
    grae->SetPointError(3,0,0,0.1259463,0.1259463);
    grae->SetPoint(4,525,-0.07831348);
    grae->SetPointError(4,0,0,0.1882171,0.1882171);
    grae->SetPoint(5,625,1.089043);
    grae->SetPointError(5,0,0,0.4240721,0.4240721);
    grae->SetPoint(6,725,-0.4112483);
    grae->SetPointError(6,0,0,0.406656,0.406656);
    grae->SetPoint(7,825,1.060315);
    grae->SetPointError(7,0,0,0.7873214,0.7873214);
    grae->SetPoint(8,925,-0.7456088);
    grae->SetPointError(8,0,0,0.6942266,0.6942266);

    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#mu#mu + jets #rightarrow #gamma + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,0,0);
    grae->SetPointError(1,0,0,0,0);
    grae->SetPoint(2,0,0);
    grae->SetPointError(2,0,0,0,0);
    grae->SetPoint(3,425,-0.1455919);
    grae->SetPointError(3,0,0,0.1806481,0.1806481);
    grae->SetPoint(4,525,-0.3547269);
    grae->SetPointError(4,0,0,0.2846392,0.2846392);
    grae->SetPoint(5,625,0.7052895);
    grae->SetPointError(5,0,0,0.5310966,0.5310966);
    grae->SetPoint(6,725,-0.45608);
    grae->SetPointError(6,0,0,0.7556384,0.7556384);
    grae->SetPoint(7,825,3.042455);
    grae->SetPointError(7,0,0,4.659075,4.659075);
    grae->SetPoint(8,925,-1);
    grae->SetPointError(8,0,0,0,0);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    return nhistos;

    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
  } else if ( set_option == 6 ) {

    int nhistos = 1;
    
    syst.clear();
    syst.push_back(0.20);
    syst.push_back(0.20);
    syst.push_back(0.20);

    init(nhistos,bins.size()-1,values);
    init(nhistos,bins.size()-1,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#mu#mu + jets #rightarrow #gamma + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,0,0);
    grae->SetPointError(1,0,0,0,0);
    grae->SetPoint(2,0,0);
    grae->SetPointError(2,0,0,0,0);
    grae->SetPoint(3,425,-0.1411698);
    grae->SetPointError(3,0,0,0.08648166,0.08648166);
    grae->SetPoint(4,525,-0.02322981);
    grae->SetPointError(4,0,0,0.1315553,0.1315553);
    grae->SetPoint(5,625,-0.3074706);
    grae->SetPointError(5,0,0,0.1873511,0.1873511);
    grae->SetPoint(6,725,-0.1557285);
    grae->SetPointError(6,0,0,0.2854318,0.2854318);
    grae->SetPoint(7,825,1.073652);
    grae->SetPointError(7,0,0,0.9323202,0.9323202);
    grae->SetPoint(8,925,-0.276403);
    grae->SetPointError(8,0,0,0.4515533,0.4515533);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    return nhistos;

    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
  } else if ( set_option == 7 ) {

    energy = "8";
    lumi = "1.6";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.20);
    syst.push_back(0.40);
    syst.push_back(0.60);
    
    init(nhistos,bins.size()-1,values);
    init(nhistos,bins.size()-1,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,-0.2768671);
    grae->SetPointError(1,0,0,0.2307345,0.2307345);
    grae->SetPoint(2,350,0.04872682);
    grae->SetPointError(2,0,0,0.2663714,0.2663714);
    grae->SetPoint(3,425,-0.2189002);
    grae->SetPointError(3,0,0,0.3141957,0.3141957);
    grae->SetPoint(4,525,1.22434);
    grae->SetPointError(4,0,0,0.6002669,0.6002669);
    grae->SetPoint(5,625,-0.05519551);
    grae->SetPointError(5,0,0,0.3348582,0.3348582);
    grae->SetPoint(6,725,0.05164244);
    grae->SetPointError(6,0,0,0.463439,0.463439);
    grae->SetPoint(7,825,-0.06994708);
    grae->SetPointError(7,0,0,0.74134,0.74134);
    grae->SetPoint(8,925,-1);
    grae->SetPointError(8,0,0,1.181829,1.181829);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,-0.07026767);
    grae->SetPointError(1,0,0,0.1520701,0.1520701);
    grae->SetPoint(2,350,-0.1515203);
    grae->SetPointError(2,0,0,0.1622125,0.1622125);
    grae->SetPoint(3,425,0.0502432);
    grae->SetPointError(3,0,0,0.1802333,0.1802333);
    grae->SetPoint(4,525,-0.2411856);
    grae->SetPointError(4,0,0,0.2254115,0.2254115);
    grae->SetPoint(5,625,0.2888916);
    grae->SetPointError(5,0,0,0.2507603,0.2507603);
    grae->SetPoint(6,725,0.7491793);
    grae->SetPointError(6,0,0,0.4094419,0.4094419);
    grae->SetPoint(7,825,0.2132012);
    grae->SetPointError(7,0,0,0.4146076,0.4146076);
    grae->SetPoint(8,925,-0.4169336);
    grae->SetPointError(8,0,0,0.445368,0.445368);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.1351485);
    grae->SetPointError(1,0,0,0.1719981,0.1719981);
    grae->SetPoint(2,350,0.2532158);
    grae->SetPointError(2,0,0,0.2127321,0.2127321);
    grae->SetPoint(3,425,-0.07493383);
    grae->SetPointError(3,0,0,0.2190089,0.2190089);
    grae->SetPoint(4,525,-0.3612194);
    grae->SetPointError(4,0,0,0.2670054,0.2670054);
    grae->SetPoint(5,625,0.05574149);
    grae->SetPointError(5,0,0,0.3459646,0.3459646);
    grae->SetPoint(6,725,-0.03813885);
    grae->SetPointError(6,0,0,0.4168869,0.4168869);
    grae->SetPoint(7,825,-0.6015736);
    grae->SetPointError(7,0,0,0.666715,0.666715);
    grae->SetPoint(8,925,0.02737636);
    grae->SetPointError(8,0,0,0.8577577,0.8577577);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.2083239);
    grae->SetPointError(1,0,0,0.1779821,0.1779821);
    grae->SetPoint(2,350,0.09422682);
    grae->SetPointError(2,0,0,0.2246674,0.2246674);
    grae->SetPoint(3,425,0.4772324);
    grae->SetPointError(3,0,0,0.2618881,0.2618881);
    grae->SetPoint(4,525,-0.274191);
    grae->SetPointError(4,0,0,0.289073,0.289073);
    grae->SetPoint(5,625,0.05885632);
    grae->SetPointError(5,0,0,0.3666106,0.3666106);
    grae->SetPoint(6,725,0.06785131);
    grae->SetPointError(6,0,0,0.6577643,0.6577643);
    grae->SetPoint(7,825,-0.2988956);
    grae->SetPointError(7,0,0,0.7079096,0.7079096);
    grae->SetPoint(8,925,-0.06152218);
    grae->SetPointError(8,0,0,0.6335361,0.6335361);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#mu#mu + jets #rightarrow #gamma + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,0,0);
    grae->SetPointError(1,0,0,0,0);
    grae->SetPoint(2,0,0);
    grae->SetPointError(2,0,0,0,0);
    grae->SetPoint(3,0,0);
    grae->SetPointError(3,0,0,0,0);
    grae->SetPoint(4,0,0);
    grae->SetPointError(4,0,0,0,0);
    grae->SetPoint(5,625,-0.4433618);
    grae->SetPointError(5,0,0,0.409085,0.409085);
    grae->SetPoint(6,725,0.1833945);
    grae->SetPointError(6,0,0,0.7837883,0.7837883);
    grae->SetPoint(7,825,-0.484331);
    grae->SetPointError(7,0,0,0.9713683,0.9713683);
    grae->SetPoint(8,925,-0.2750741);
    grae->SetPointError(8,0,0,0.9380708,0.9380708);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    return nhistos;

    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
  } else if ( set_option == 8 ) {

    energy = "8";
    lumi = "3.9";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.20);
    syst.push_back(0.30);
    syst.push_back(0.70);
    
    init(nhistos,bins.size()-1,values);
    init(nhistos,bins.size()-1,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,-0.01912947);
    grae->SetPointError(1,0,0,0.08938719,0.08938719);
    grae->SetPoint(2,350,-0.1222202);
    grae->SetPointError(2,0,0,0.07609002,0.07609002);
    grae->SetPoint(3,425,-0.03527486);
    grae->SetPointError(3,0,0,0.06820642,0.06820642);
    grae->SetPoint(4,525,-0.1020959);
    grae->SetPointError(4,0,0,0.111265,0.111265);
    grae->SetPoint(5,625,-0.2678577);
    grae->SetPointError(5,0,0,0.1854963,0.1854963);
    grae->SetPoint(6,725,-0.4320945);
    grae->SetPointError(6,0,0,0.3378501,0.3378501);
    grae->SetPoint(7,825,0.1432385);
    grae->SetPointError(7,0,0,0.4316852,0.4316852);
    grae->SetPoint(8,925,-0.7577979);
    grae->SetPointError(8,0,0,0.6817167,0.6817167);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.1675455);
    grae->SetPointError(1,0,0,0.08771491,0.08771491);
    grae->SetPoint(2,350,-0.006801409);
    grae->SetPointError(2,0,0,0.06597502,0.06597502);
    grae->SetPoint(3,425,0.03028101);
    grae->SetPointError(3,0,0,0.06517107,0.06517107);
    grae->SetPoint(4,525,-0.01456659);
    grae->SetPointError(4,0,0,0.09853544,0.09853544);
    grae->SetPoint(5,625,0.08209643);
    grae->SetPointError(5,0,0,0.1506132,0.1506132);
    grae->SetPoint(6,725,0.2929802);
    grae->SetPointError(6,0,0,0.2183663,0.2183663);
    grae->SetPoint(7,825,-0.03940762);
    grae->SetPointError(7,0,0,0.3181101,0.3181101);
    grae->SetPoint(8,925,-0.5127236);
    grae->SetPointError(8,0,0,0.3083436,0.3083436);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.0734074);
    grae->SetPointError(1,0,0,0.09296023,0.09296023);
    grae->SetPoint(2,350,-0.03892464);
    grae->SetPointError(2,0,0,0.09608791,0.09608791);
    grae->SetPoint(3,425,0.1032287);
    grae->SetPointError(3,0,0,0.09616556,0.09616556);
    grae->SetPoint(4,525,0.2036508);
    grae->SetPointError(4,0,0,0.1467392,0.1467392);
    grae->SetPoint(5,625,-0.03596342);
    grae->SetPointError(5,0,0,0.207465,0.207465);
    grae->SetPoint(6,725,0.02051861);
    grae->SetPointError(6,0,0,0.2957373,0.2957373);
    grae->SetPoint(7,825,-0.06235633);
    grae->SetPointError(7,0,0,0.4414847,0.4414847);
    grae->SetPoint(8,925,0.4326466);
    grae->SetPointError(8,0,0,0.5975056,0.5975056);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.1412435);
    grae->SetPointError(1,0,0,0.09142167,0.09142167);
    grae->SetPoint(2,350,-0.06178053);
    grae->SetPointError(2,0,0,0.1022141,0.1022141);
    grae->SetPoint(3,425,-0.06250405);
    grae->SetPointError(3,0,0,0.1150601,0.1150601);
    grae->SetPoint(4,525,0.2754829);
    grae->SetPointError(4,0,0,0.1728933,0.1728933);
    grae->SetPoint(5,625,0.1990504);
    grae->SetPointError(5,0,0,0.2511346,0.2511346);
    grae->SetPoint(6,725,0.3355871);
    grae->SetPointError(6,0,0,0.3278686,0.3278686);
    grae->SetPoint(7,825,0.6662954);
    grae->SetPointError(7,0,0,0.5921539,0.5921539);
    grae->SetPoint(8,925,-0.2322451);
    grae->SetPointError(8,0,0,0.5188978,0.5188978);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#mu#mu + jets #rightarrow #gamma + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,0,0);
    grae->SetPointError(1,0,0,0,0);
    grae->SetPoint(2,0,0);
    grae->SetPointError(2,0,0,0,0);
    grae->SetPoint(3,425,-0.0745309);
    grae->SetPointError(3,0,0,0.1208496,0.1208496);
    grae->SetPoint(4,525,0.2166894);
    grae->SetPointError(4,0,0,0.1821702,0.1821702);
    grae->SetPoint(5,625,0.05906951);
    grae->SetPointError(5,0,0,0.2648079,0.2648079);
    grae->SetPoint(6,725,-0.01750412);
    grae->SetPointError(6,0,0,0.3328937,0.3328937);
    grae->SetPoint(7,825,0.4635275);
    grae->SetPointError(7,0,0,0.5991121,0.5991121);
    grae->SetPoint(8,925,-0.07715072);
    grae->SetPointError(8,0,0,0.6037971,0.6037971);

    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    return nhistos;


    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
  } else if ( set_option == 9 ) {

    energy = "8";
    lumi = "3.9";

    int nhistos = 3;
    
    syst.clear();
    syst.push_back(0.20);
    syst.push_back(0.40);
    syst.push_back(0.60);
    
    init(nhistos,bins.size()-1,values);
    init(nhistos,bins.size()-1,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    //     // mu: no aT -> with aT
    //     grae = new TGraphAsymmErrors(8);
    //     grae->SetName("Graph");
    //     grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    //     grae->SetFillColor(1);
    //     grae->SetLineWidth(3);
    //     grae->SetMarkerStyle(20);
    //     grae->SetMarkerSize(1.5);
    //    grae->SetPoint(1,300,-0.2640229);
    //    grae->SetPointError(1,0,0,0.2835055,0.2835055);
    //    grae->SetPoint(2,350,-0.1107953);
    //    grae->SetPointError(2,0,0,0.326635,0.326635);
    //    grae->SetPoint(3,425,0.2781013);
    //    grae->SetPointError(3,0,0,0.2572024,0.2572024);
    //    grae->SetPoint(4,525,0.8828158);
    //    grae->SetPointError(4,0,0,0.5013902,0.5013902);
    //    grae->SetPoint(5,625,0.04499625);
    //    grae->SetPointError(5,0,0,0.5721765,0.5721765);
    //    grae->SetPoint(6,725,-0.1100501);
    //    grae->SetPointError(6,0,0,0.5635164,0.5635164);
    //    grae->SetPoint(7,825,8.209098);
    //    grae->SetPointError(7,0,0,14.97828,14.97828);
    //    grae->SetPoint(8,925,-1);
    //    grae->SetPointError(8,0,0,1.45329,1.45329);
    
    //     for ( uint i = 0; i < 8; ++i ) {
    //       double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
    //       errors[index][i] = grae->GetErrorYhigh(i+1);
    //     }
    //     titles.push_back( grae->GetTitle() );
    //     index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.109504);
    grae->SetPointError(1,0,0,0.08403678,0.08403678);
    grae->SetPoint(2,350,-0.05021044);
    grae->SetPointError(2,0,0,0.06424859,0.06424859);
    grae->SetPoint(3,425,0.0004423723);
    grae->SetPointError(3,0,0,0.06327076,0.06327076);
    grae->SetPoint(4,525,-0.01307169);
    grae->SetPointError(4,0,0,0.09691845,0.09691845);
    grae->SetPoint(5,625,0.09159275);
    grae->SetPointError(5,0,0,0.1482814,0.1482814);
    grae->SetPoint(6,725,0.3601117);
    grae->SetPointError(6,0,0,0.2206155,0.2206155);
    grae->SetPoint(7,825,-0.0105901);
    grae->SetPointError(7,0,0,0.3182493,0.3182493);
    grae->SetPoint(8,925,-0.5004471);
    grae->SetPointError(8,0,0,0.3076232,0.3076232);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.08576673);
    grae->SetPointError(1,0,0,0.09047675,0.09047675);
    grae->SetPoint(2,350,-0.01624375);
    grae->SetPointError(2,0,0,0.0950765,0.0950765);
    grae->SetPoint(3,425,0.1214643);
    grae->SetPointError(3,0,0,0.09529702,0.09529702);
    grae->SetPoint(4,525,0.2140894);
    grae->SetPointError(4,0,0,0.145538,0.145538);
    grae->SetPoint(5,625,-0.02372244);
    grae->SetPointError(5,0,0,0.2048419,0.2048419);
    grae->SetPoint(6,725,-0.02159368);
    grae->SetPointError(6,0,0,0.2891405,0.2891405);
    grae->SetPoint(7,825,-0.04972092);
    grae->SetPointError(7,0,0,0.4416311,0.4416311);
    grae->SetPoint(8,925,0.457147);
    grae->SetPointError(8,0,0,0.6061101,0.6061101);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.1507975);
    grae->SetPointError(1,0,0,0.09006486,0.09006486);
    grae->SetPoint(2,350,-0.05207282);
    grae->SetPointError(2,0,0,0.1005548,0.1005548);
    grae->SetPoint(3,425,-0.04101989);
    grae->SetPointError(3,0,0,0.1127603,0.1127603);
    grae->SetPoint(4,525,0.2891903);
    grae->SetPointError(4,0,0,0.1703224,0.1703224);
    grae->SetPoint(5,625,0.2190631);
    grae->SetPointError(5,0,0,0.2511198,0.2511198);
    grae->SetPoint(6,725,0.4755818);
    grae->SetPointError(6,0,0,0.3451707,0.3451707);
    grae->SetPoint(7,825,0.7018587);
    grae->SetPointError(7,0,0,0.6042459,0.6042459);
    grae->SetPoint(8,925,-0.1291913);
    grae->SetPointError(8,0,0,0.5131332,0.5131332);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    //     // mumu->gamma, no aT
    //     grae = new TGraphAsymmErrors(8);
    //     grae->SetName("Graph");
    //     grae->SetTitle("#mu#mu + jets #rightarrow #gamma + jets");
    //     grae->SetFillColor(1);
    //     grae->SetLineWidth(3);
    //     grae->SetMarkerStyle(20);
    //     grae->SetMarkerSize(1.5);
    //    grae->SetPoint(1,0,0);
    //    grae->SetPointError(1,0,0,0,0);
    //    grae->SetPoint(2,0,0);
    //    grae->SetPointError(2,0,0,0,0);
    //    grae->SetPoint(3,0,0);
    //    grae->SetPointError(3,0,0,0,0);
    //    grae->SetPoint(4,0,0);
    //    grae->SetPointError(4,0,0,0,0);
    //    grae->SetPoint(5,625,-0.4433618);
    //    grae->SetPointError(5,0,0,0.409085,0.409085);
    //    grae->SetPoint(6,725,0.1833945);
    //    grae->SetPointError(6,0,0,0.7837883,0.7837883);
    //    grae->SetPoint(7,825,-0.484331);
    //    grae->SetPointError(7,0,0,0.9713683,0.9713683);
    //    grae->SetPoint(8,925,-0.2750741);
    //    grae->SetPointError(8,0,0,0.9380708,0.9380708);
    
    //     for ( uint i = 0; i < 8; ++i ) {
    //       double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
    //       if ( x < 275. ) { y = -1000.; }
    //       values[index][i] = y;
    //       errors[index][i] = grae->GetErrorYhigh(i+1);
    //     }
    //     titles.push_back( grae->GetTitle() );
    //     index++;

    return nhistos;

    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
  } else if ( set_option == 10 ) {

    energy = "8";
    lumi = "3.9";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.20);
    syst.push_back(0.40);
    syst.push_back(0.60);
    
    init(nhistos,bins.size()-1,values);
    init(nhistos,bins.size()-1,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    //     // mu, with aT -> mu, no aT 
    //     grae = new TGraphAsymmErrors(8);
    //     grae->SetName("Graph");
    //     grae->SetTitle("#alpha_{T} > 0.55 #rightarrow #alpha_{T} < 0.55 (#mu + jets)");
    //     grae->SetFillColor(1);
    //     grae->SetLineWidth(3);
    //     grae->SetMarkerStyle(20);
    //     grae->SetMarkerSize(1.5);
    //    grae->SetPoint(1,300,-0.01912947);
    //    grae->SetPointError(1,0,0,0.08938719,0.08938719);
    //    grae->SetPoint(2,350,-0.1222202);
    //    grae->SetPointError(2,0,0,0.07609002,0.07609002);
    //    grae->SetPoint(3,425,-0.03527486);
    //    grae->SetPointError(3,0,0,0.06820642,0.06820642);
    //    grae->SetPoint(4,525,-0.1020959);
    //    grae->SetPointError(4,0,0,0.111265,0.111265);
    //    grae->SetPoint(5,625,-0.2678577);
    //    grae->SetPointError(5,0,0,0.1854963,0.1854963);
    //    grae->SetPoint(6,725,-0.4320945);
    //    grae->SetPointError(6,0,0,0.3378501,0.3378501);
    //    grae->SetPoint(7,825,0.1432385);
    //    grae->SetPointError(7,0,0,0.4316852,0.4316852);
    //    grae->SetPoint(8,925,-0.7577979);
    //    grae->SetPointError(8,0,0,0.6817167,0.6817167);
    
    //     for ( uint i = 0; i < 8; ++i ) {
    //       double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
    //       errors[index][i] = grae->GetErrorYhigh(i+1);
    //     }
    //     titles.push_back( grae->GetTitle() );
    //     index++;

    // mu, with aT -> mu, no aT (=0b) 
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets, 0 b tags)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.01897203);
    grae->SetPointError(1,0,0,0.07547812,0.07547812);
    grae->SetPoint(2,350,-0.1487618);
    grae->SetPointError(2,0,0,0.0964829,0.0964829);
    grae->SetPoint(3,425,0.04572239);
    grae->SetPointError(3,0,0,0.079521,0.079521);
    grae->SetPoint(4,525,-0.2035728);
    grae->SetPointError(4,0,0,0.1298842,0.1298842);
    grae->SetPoint(5,625,-0.1882691);
    grae->SetPointError(5,0,0,0.2114145,0.2114145);
    grae->SetPoint(6,725,-0.4067594);
    grae->SetPointError(6,0,0,0.4121097,0.4121097);
    grae->SetPoint(7,825,-0.1073573);
    grae->SetPointError(7,0,0,0.4811773,0.4811773);
    grae->SetPoint(8,925,-0.7421818);
    grae->SetPointError(8,0,0,0.715033,0.715033);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, with aT -> mu, no aT (>0b) 
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets, #geq 1 b tags)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,-0.02540268);
    grae->SetPointError(1,0,0,0.1351748,0.1351748);
    grae->SetPoint(2,350,-0.04726063);
    grae->SetPointError(2,0,0,0.2239728,0.2239728);
    grae->SetPoint(3,425,0.06647788);
    grae->SetPointError(3,0,0,0.211705,0.211705);
    grae->SetPoint(4,525,0.1550522);
    grae->SetPointError(4,0,0,0.3979539,0.3979539);
    grae->SetPoint(5,625,-0.1405223);
    grae->SetPointError(5,0,0,0.6016138,0.6016138);
    grae->SetPoint(6,725,-1);
    grae->SetPointError(6,0,0,1.043344,1.043344);
    grae->SetPoint(7,825,2.517561);
    grae->SetPointError(7,0,0,5.300113,5.300113);
    grae->SetPoint(8,925,-1);
    grae->SetPointError(8,0,0,4.44563,4.44563);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu: no aT -> with aT XXX
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.1895047);
    grae->SetPointError(1,0,0,0.1807637,0.1807637);
    grae->SetPoint(2,350,0.3722965);
    grae->SetPointError(2,0,0,0.29114,0.29114);
    grae->SetPoint(3,425,-0.3631829);
    grae->SetPointError(3,0,0,0.2931261,0.2931261);
    grae->SetPoint(4,525,-0.3998448);
    grae->SetPointError(4,0,0,0.4199226,0.4199226);
    grae->SetPoint(5,625,-0.8023914);
    grae->SetPointError(5,0,0,0.7934083,0.7934083);
    grae->SetPoint(6,725,1.097085);
    grae->SetPointError(6,0,0,1.857348,1.857348);
    grae->SetPoint(7,825,-1);
    grae->SetPointError(7,0,0,0,0);
    grae->SetPoint(8,925,-1);
    grae->SetPointError(8,0,0,11.33996,11.33996);

    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, with aT -> mumu, no aT 
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} > 0.55 (#mu + jets) #rightarrow #alpha_{T} < 0.55 (#mu#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.1635007);
    grae->SetPointError(1,0,0,0.1189439,0.1189439);
    grae->SetPoint(2,350,0.06885518);
    grae->SetPointError(2,0,0,0.1242743,0.1242743);
    grae->SetPoint(3,425,-0.02822481);
    grae->SetPointError(3,0,0,0.1297572,0.1297572);
    grae->SetPoint(4,525,0.4205114);
    grae->SetPointError(4,0,0,0.2154006,0.2154006);
    grae->SetPoint(5,625,0.6377286);
    grae->SetPointError(5,0,0,0.3777515,0.3777515);
    grae->SetPoint(6,725,1.351777);
    grae->SetPointError(6,0,0,0.8767194,0.8767194);
    grae->SetPoint(7,825,0.4575222);
    grae->SetPointError(7,0,0,0.6673845,0.6673845);
    grae->SetPoint(8,925,2.169893);
    grae->SetPointError(8,0,0,2.977805,2.977805);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu, no aT -> gamma, with aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 (#mu#mu + jets) #rightarrow #alpha_{T} > 0.55 (#gamma + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,0,0);
    grae->SetPointError(1,0,0,0,0);
    grae->SetPoint(2,0,0);
    grae->SetPointError(2,0,0,0,0);
    grae->SetPoint(3,425,-0.0745309);
    grae->SetPointError(3,0,0,0.1208496,0.1208496);
    grae->SetPoint(4,525,0.2166894);
    grae->SetPointError(4,0,0,0.1821702,0.1821702);
    grae->SetPoint(5,625,0.05906951);
    grae->SetPointError(5,0,0,0.2648079,0.2648079);
    grae->SetPoint(6,725,-0.01750412);
    grae->SetPointError(6,0,0,0.3328937,0.3328937);
    grae->SetPoint(7,825,0.4635275);
    grae->SetPointError(7,0,0,0.5991121,0.5991121);
    grae->SetPoint(8,925,-0.07715072);
    grae->SetPointError(8,0,0,0.6037971,0.6037971);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    return nhistos;

    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
  } else if ( set_option == 11 ) { // 2011, aT acceptance

    energy = "7";
    lumi = "5.0";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.20);
    syst.push_back(0.40);
    syst.push_back(0.60);
    
    init(nhistos,bins.size()-1,values);
    init(nhistos,bins.size()-1,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu, with aT -> mu, no aT (=0b) 
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets, 0 b tags)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,0,0);
    grae->SetPointError(1,0,0,0,0);
    grae->SetPoint(2,0,0);
    grae->SetPointError(2,0,0,0,0);
    grae->SetPoint(3,425,0.0644265);
    grae->SetPointError(3,0,0,0.06956349,0.06956349);
    grae->SetPoint(4,525,-0.1121516);
    grae->SetPointError(4,0,0,0.1180477,0.1180477);
    grae->SetPoint(5,625,-0.28351);
    grae->SetPointError(5,0,0,0.1991521,0.1991521);
    grae->SetPoint(6,725,0.005680904);
    grae->SetPointError(6,0,0,0.3208864,0.3208864);
    grae->SetPoint(7,825,0.06888876);
    grae->SetPointError(7,0,0,0.4895633,0.4895633);
    grae->SetPoint(8,925,-1);
    grae->SetPointError(8,0,0,0.9221425,0.9221425);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, with aT -> mu, no aT (>0b) 
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets, #geq 1 b tags)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,0,0);
    grae->SetPointError(1,0,0,0,0);
    grae->SetPoint(2,0,0);
    grae->SetPointError(2,0,0,0,0);
    grae->SetPoint(3,425,0.1279024);
    grae->SetPointError(3,0,0,0.1907289,0.1907289);
    grae->SetPoint(4,525,-0.05648134);
    grae->SetPointError(4,0,0,0.2678968,0.2678968);
    grae->SetPoint(5,625,0.3720487);
    grae->SetPointError(5,0,0,0.4448949,0.4448949);
    grae->SetPoint(6,725,0.8280071);
    grae->SetPointError(6,0,0,1.107324,1.107324);
    grae->SetPoint(7,825,-1);
    grae->SetPointError(7,0,0,2.478996,2.478996);
    grae->SetPoint(8,925,-1);
    grae->SetPointError(8,0,0,2.361392,2.361392);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu: no aT -> with aT XXX
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,0,0);
    grae->SetPointError(1,0,0,0,0);
    grae->SetPoint(2,0,0);
    grae->SetPointError(2,0,0,0,0);
    grae->SetPoint(3,425,-0.0699633);
    grae->SetPointError(3,0,0,0.2094057,0.2094057);
    grae->SetPoint(4,525,-0.002593247);
    grae->SetPointError(4,0,0,0.3709601,0.3709601);
    grae->SetPoint(5,625,-0.4440163);
    grae->SetPointError(5,0,0,0.5176479,0.5176479);
    grae->SetPoint(6,725,-0.5304814);
    grae->SetPointError(6,0,0,0.9670614,0.9670614);
    grae->SetPoint(7,825,-1);
    grae->SetPointError(7,0,0,6.879462,6.879462);
    grae->SetPoint(8,925,-0.2125803);
    grae->SetPointError(8,0,0,1.127312,1.127312);

    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, with aT -> mumu, no aT 
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} > 0.55 (#mu + jets) #rightarrow #alpha_{T} < 0.55 (#mu#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,0,0);
    grae->SetPointError(1,0,0,0,0);
    grae->SetPoint(2,0,0);
    grae->SetPointError(2,0,0,0,0);
    grae->SetPoint(3,425,-0.09714969);
    grae->SetPointError(3,0,0,0.1001485,0.1001485);
    grae->SetPoint(4,525,-0.05980608);
    grae->SetPointError(4,0,0,0.1696416,0.1696416);
    grae->SetPoint(5,625,0.5411776);
    grae->SetPointError(5,0,0,0.3333413,0.3333413);
    grae->SetPoint(6,725,0.1809098);
    grae->SetPointError(6,0,0,0.4118215,0.4118215);
    grae->SetPoint(7,825,-0.6185918);
    grae->SetPointError(7,0,0,0.6944006,0.6944006);
    grae->SetPoint(8,925,-1);
    grae->SetPointError(8,0,0,0,0);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu, no aT -> gamma, with aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 (#mu#mu + jets) #rightarrow #alpha_{T} > 0.55 (#gamma + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,0,0);
    grae->SetPointError(1,0,0,0,0);
    grae->SetPoint(2,0,0);
    grae->SetPointError(2,0,0,0,0);
    grae->SetPoint(3,425,0.01719715);
    grae->SetPointError(3,0,0,0.1829103,0.1829103);
    grae->SetPoint(4,525,-0.09022204);
    grae->SetPointError(4,0,0,0.298347,0.298347);
    grae->SetPoint(5,625,0.005966493);
    grae->SetPointError(5,0,0,0.4429791,0.4429791);
    grae->SetPoint(6,725,-0.612747);
    grae->SetPointError(6,0,0,0.8840258,0.8840258);
    grae->SetPoint(7,825,-1);
    grae->SetPointError(7,0,0,4.07869,4.07869);
    grae->SetPoint(8,925,0.1471287);
    grae->SetPointError(8,0,0,0.9906887,0.9906887);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    return nhistos;

  } else if ( set_option == 20 ) { // pdf

    energy = "7";
    lumi = "5.0";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.20);
    syst.push_back(0.40);
    syst.push_back(0.60);
    
    init(nhistos,bins.size()-1,values);
    init(nhistos,bins.size()-1,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu, with aT -> mu, no aT (=0b) 
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets, 0 b tags)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,0,0);
    grae->SetPointError(1,0,0,0,0);
    grae->SetPoint(2,0,0);
    grae->SetPointError(2,0,0,0,0);
    grae->SetPoint(3,425,0.0644265);
    grae->SetPointError(3,0,0,0.06956349,0.06956349);
    grae->SetPoint(4,525,-0.1121516);
    grae->SetPointError(4,0,0,0.1180477,0.1180477);
    grae->SetPoint(5,625,-0.28351);
    grae->SetPointError(5,0,0,0.1991521,0.1991521);
    grae->SetPoint(6,725,0.005680904);
    grae->SetPointError(6,0,0,0.3208864,0.3208864);
    grae->SetPoint(7,825,0.06888876);
    grae->SetPointError(7,0,0,0.4895633,0.4895633);
    grae->SetPoint(8,925,-1);
    grae->SetPointError(8,0,0,0.9221425,0.9221425);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, with aT -> mu, no aT (>0b) 
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets, #geq 1 b tags)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,0,0);
    grae->SetPointError(1,0,0,0,0);
    grae->SetPoint(2,0,0);
    grae->SetPointError(2,0,0,0,0);
    grae->SetPoint(3,425,0.1279024);
    grae->SetPointError(3,0,0,0.1907289,0.1907289);
    grae->SetPoint(4,525,-0.05648134);
    grae->SetPointError(4,0,0,0.2678968,0.2678968);
    grae->SetPoint(5,625,0.3720487);
    grae->SetPointError(5,0,0,0.4448949,0.4448949);
    grae->SetPoint(6,725,0.8280071);
    grae->SetPointError(6,0,0,1.107324,1.107324);
    grae->SetPoint(7,825,-1);
    grae->SetPointError(7,0,0,2.478996,2.478996);
    grae->SetPoint(8,925,-1);
    grae->SetPointError(8,0,0,2.361392,2.361392);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu: no aT -> with aT XXX
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,0,0);
    grae->SetPointError(1,0,0,0,0);
    grae->SetPoint(2,0,0);
    grae->SetPointError(2,0,0,0,0);
    grae->SetPoint(3,425,-0.0699633);
    grae->SetPointError(3,0,0,0.2094057,0.2094057);
    grae->SetPoint(4,525,-0.002593247);
    grae->SetPointError(4,0,0,0.3709601,0.3709601);
    grae->SetPoint(5,625,-0.4440163);
    grae->SetPointError(5,0,0,0.5176479,0.5176479);
    grae->SetPoint(6,725,-0.5304814);
    grae->SetPointError(6,0,0,0.9670614,0.9670614);
    grae->SetPoint(7,825,-1);
    grae->SetPointError(7,0,0,6.879462,6.879462);
    grae->SetPoint(8,925,-0.2125803);
    grae->SetPointError(8,0,0,1.127312,1.127312);

    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, with aT -> mumu, no aT 
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} > 0.55 (#mu + jets) #rightarrow #alpha_{T} < 0.55 (#mu#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,0,0);
    grae->SetPointError(1,0,0,0,0);
    grae->SetPoint(2,0,0);
    grae->SetPointError(2,0,0,0,0);
    grae->SetPoint(3,425,-0.09714969);
    grae->SetPointError(3,0,0,0.1001485,0.1001485);
    grae->SetPoint(4,525,-0.05980608);
    grae->SetPointError(4,0,0,0.1696416,0.1696416);
    grae->SetPoint(5,625,0.5411776);
    grae->SetPointError(5,0,0,0.3333413,0.3333413);
    grae->SetPoint(6,725,0.1809098);
    grae->SetPointError(6,0,0,0.4118215,0.4118215);
    grae->SetPoint(7,825,-0.6185918);
    grae->SetPointError(7,0,0,0.6944006,0.6944006);
    grae->SetPoint(8,925,-1);
    grae->SetPointError(8,0,0,0,0);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu, no aT -> gamma, with aT
    grae = new TGraphAsymmErrors(8);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 (#mu#mu + jets) #rightarrow #alpha_{T} > 0.55 (#gamma + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,0,0);
    grae->SetPointError(1,0,0,0,0);
    grae->SetPoint(2,0,0);
    grae->SetPointError(2,0,0,0,0);
    grae->SetPoint(3,425,0.01719715);
    grae->SetPointError(3,0,0,0.1829103,0.1829103);
    grae->SetPoint(4,525,-0.09022204);
    grae->SetPointError(4,0,0,0.298347,0.298347);
    grae->SetPoint(5,625,0.005966493);
    grae->SetPointError(5,0,0,0.4429791,0.4429791);
    grae->SetPoint(6,725,-0.612747);
    grae->SetPointError(6,0,0,0.8840258,0.8840258);
    grae->SetPoint(7,825,-1);
    grae->SetPointError(7,0,0,4.07869,4.07869);
    grae->SetPoint(8,925,0.1471287);
    grae->SetPointError(8,0,0,0.9906887,0.9906887);
    
    for ( uint i = 0; i < 8; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    return nhistos;

  }

  return 0;

}

