#include "TMultiGraph.h"
#include "TLatex.h"
#include "TBox.h"
#include "TLegend.h"
#include "TLine.h"
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
#include "Math/DistFunc.h"

using namespace std;

typedef unsigned int uint;
typedef vector<double> vdouble;
typedef vector<string> vstring;
typedef vector<vdouble> vvdouble;

float round( double val, int places ) {
  return (int( val * pow(10.,places) )*1.)/(1.*pow(10.,places));
}

// options
/*
  2011 = 2011 paper
  0 = ICHEP PAS
  1 = 5/fb "top up", inclusive, 532p4
  2 = 5/fb "top up", 2-3 jets, 532p4
  3 = 5/fb "top up", >3 jets, 532p4
  4 = 8/fb "top up", inclusive, 532p4
  5 = 8/fb "top up", 2-3 jets, 532p4
  6 = 8/fb "top up", >3 jets, 532p4
  7 = 5/fb, inclusive, 553p2, nominal
  8 = 5/fb, inclusive, 553p2, LO WJets
  9 = 5/fb, inclusive, 553p2, LO WJets/DY/Zinv
  10 = 11/fb, inclusive, 553p2, nominal
  11 = 11/fb, inclusive, 553p2, LO WJets/DY/Zinv

  110 = 11/fb, inclusive, 553p2, Darren's corrections
  111 = 11/fb, inclusive, 553p2, NNLO XS

  12 = 5/fb, inclusive, 553p2, nominal, 8 bins
  13 = 5/fb, inclusive, 553p2, LO WJets/DY/Zinv, 8 bins
  14 = 11/fb, inclusive, 553p2, LO WJets/DY/Zinv, 8 bins

  15 = 11/fb, 2-3, 553p2, LO WJets/DY/Zinv, 8 bins
  16 = 11/fb, >=4, 553p2, LO WJets/DY/Zinv, 8 bins
  17 = 5/fb, 2-3, 553p2, LO WJets/DY/Zinv, 8 bins
  18 = 5/fb, >=4, 553p2, LO WJets/DY/Zinv, 8 bins

  19 = 11/fb, >=4, 553p2, LO WJets/DY/Zinv, 8 bins, no mumum/gamma

  20 = 11/fb, 2-3, 553p2, LO WJets/DY/Zinv, 8 bins, 8 tests
  21 = 11/fb, >=4, 553p2, LO WJets/DY/Zinv, 8 bins, 8 tests

  22 = 11/fb, 2-3, 553p2, Darren's corrections, 8 bins, 8 tests
  23 = 11/fb, >=4, 553p2, Darren's corrections, 8 bins, 8 tests
  230 = 11/fb, >=4, 553p2, Darren's corrections, 8 bins, JUST 2,3 -> >=4, mu+jets

  24 = 11/fb, 2-3, 553p2, effect of varying W+20% / TTbar-20%, 8 bins, 8 tests
  25 = 11/fb, >=4, 553p2, effect of varying W+20% / TTbar-20%, 8 bins, 8 tests


  31 = 11/fb, 553p2, alphaT acceptance

  42 = 18/fb, 2-3, 553p2, Darren's corrections, 8 bins, 8 tests
  43 = 18/fb, >=4, 553p2, Darren's corrections, 8 bins, 8 tests

  119 = 11/fb, diff jet multiplicity bins, diff b-tab bins, LO V+jets, 10 bins
  120 = 11/fb, diff jet multiplicity bins, incl. b-tag bin, LO V+jets, 10 bins

  130 = Different TTbar samples

*/

int set_option = 22; uint nbins = 8;

int syst_option = 5; // 0 = std. dev., 1 = sqrt( (mean-1)^ 2 + var ), 2 = X% coverage, 3 = paris, 4 = louis, 5 = sample variance
bool weighted = true;
bool draw_weighted = false;
bool draw_syst = true;
bool detailed_syst = false;
bool detailed_fit = false;
bool fixed_syst_bands = true;
std::string lumi = "4.98";
std::string energy = "7";
double confidence_level = 0.682689; //0.9973 //0.9545 
std::string njet = "";

vdouble bins;
vdouble widths;
vdouble syst; 
vdouble markers;

void setTDRStyle();
void init( uint, vdouble& );
void init( uint, uint, vvdouble& );
int extractData( vvdouble&, vvdouble&, vstring& );
void extractNumbers( vstring&, vvdouble&, vvdouble&, vdouble&, vdouble&, vdouble&, vdouble&, bool );
void fitResults( TGraphErrors*, double&, double&, double&, int&, std::string );
double fabs( double x ) { return x<0.?x*-1.:x; }

typedef std::pair<int,double> Pair;
typedef std::vector<Pair> Pairs;
bool ascending( const Pair& x, const Pair& y ) { return  x.second < y.second; }

void btag2012() {


  std::cout << " 1sigma: " << 1.-(1.-ROOT::Math::normal_cdf(1.))*2.
	    << " 2sigma: " << 1.-(1.-ROOT::Math::normal_cdf(2.))*2.
	    << " 3sigma: " << 1.-(1.-ROOT::Math::normal_cdf(3.))*2.
	    << std::endl;

  setTDRStyle();

  // Bins and widths
  bins.push_back(275.);
  bins.push_back(325.);
  for ( uint i = 0; i < nbins-2; ++i ) { bins.push_back(375.+100.*i); }
  bins.push_back(375.+100.*(nbins-2));
  
  widths.push_back(25.);
  widths.push_back(25.);
  for ( uint i = 0; i < nbins-2; ++i ) { widths.push_back(50.); }

  markers.push_back(24);
  markers.push_back(25);
  markers.push_back(26);
  markers.push_back(5);
  markers.push_back(20);
  markers.push_back(32);
  markers.push_back(27);
  markers.push_back(31);
  
  // Extract data from histograms
  vvdouble values; 
  vvdouble errors;
  vstring titles;
  const int nhistos = extractData( values, errors, titles );

  // Extract central values, errors, and bands for each HT bin
  vdouble n1; for ( uint i = 0; i < nbins; ++i ) { n1.push_back(i); }
  vdouble m1; 
  vdouble s1;
  vdouble u1; 
  extractNumbers( titles, values, errors, n1, m1, s1, u1, false ); 

  // Extract central values, errors, and bands for each HT region
  vdouble n2; n2.push_back(0); n2.push_back(1); n2.push_back(2); n2.push_back(4); n2.push_back(6); 
  vdouble m2; 
  vdouble s2;
  vdouble u2; 
  extractNumbers( titles, values, errors, n2, m2, s2, u2, true ); 
  
  vdouble m; 
  m.resize(m2.size(),0.);

  TLegend* leg = new TLegend( 0.15, 0.91-0.045*(nhistos+(draw_weighted?2:1)), 0.35, 0.91 );
  leg->SetFillColor(0);
  leg->SetLineColor(0); 
  leg->SetShadowColor(0); 
  leg->SetTextSize(0.035);

  TMultiGraph* mg = new TMultiGraph();
  int offset = 0;
  for ( int i = 0; i < nhistos+2; ++i ) { 
    offset++;

    if ( i > 0 && i <= nhistos ) { // All individual closure tests
      
      //continue;
      
      vdouble e(widths.size(),0.);
      vdouble x(bins); for ( uint j = 0; j < x.size(); ++j ) { x[j] += widths[j] + 5.*(offset%2?1*offset/2:-1*offset/2); }
      
      TGraphErrors* rob1 = new TGraphErrors( nbins,
					     &x.front(),
					     &values[i-1].front(),
					     &e.front(),
					     &errors[i-1].front() );
      rob1->SetMarkerStyle(markers[i-1]);
      if ( rob1->GetMarkerStyle() != 20 ) { 
	rob1->SetMarkerSize(2.0);
      } else {
	rob1->SetMarkerSize(1.0);
      }
      rob1->SetLineWidth(2);
      rob1->SetTitle("");
      mg->Add(rob1,"pZ");

//       std::cout << "Plot #: " << i-1
// 		<< " Nbins: " << bins.size()
// 		<< " Title: " << titles[i-1] << std::endl;
//       for ( uint ii = 0; ii < nbins; ++ii ) { 
// 	std::cout << std::fixed
// 		  << std::setprecision(0)
// 		  << " HT: " << bins[ii]
// 	  //<< " width: " << e[ii]
// 		  << std::setprecision(2)
// 		  << " val: " << values[i-1][ii]
// 		  << " err: " << errors[i-1][ii]
// 		  << " pull: " << (errors[i-1][ii]>0.?(values[i-1][ii]/errors[i-1][ii]):0.)
// 		  << std::endl;
//       }

      double fmean,ferr,fchi2; int fdof;
      fitResults(rob1,fmean,ferr,fchi2,fdof,titles[i-1]);
      
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

    } else if ( draw_syst && i == 0 ) { // Systematic bands for HT regions

      vdouble x; for ( uint j = 0; j < n2.size(); ++j ) { x.push_back( bins[n2[j]] ); } 
      vdouble xel(n2.size(),0.);
      vdouble xeh; for ( uint j = 0; j < n2.size(); ++j ) { xeh.push_back( (j<n2.size()-1)?x[j+1]-x[j]:bins[nbins]-x[j] ); }
      
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
	ss << "Systematic (";
	//<< int(100.*confidence_level) << "% coverage yields: ";
	for ( uint iii = 0; iii < u2.size(); ++iii ) {
	  if ( iii < u2.size()-1 ) { ss << int(100.*u2[iii]) << "%, "; }
	  else                     { ss << int(100.*u2[iii]) << "%)"; }
	}
      } else {
	ss << "Systematic uncertainty";
      }
      leg->AddEntry(tmp,ss.str().c_str(),"f");

    } else if ( draw_weighted && i > nhistos ) { // Bin-averaged result
      
      vdouble e(widths.size(),0.);
      vdouble x; x.reserve(nbins); 
      std::copy( bins.begin(), bins.end()-1, back_inserter(x) ); 
      for ( uint j = 0; j < x.size(); ++j ) { x[j] += widths[j]; }
	
      TGraphErrors* tmp2 = new TGraphErrors( nbins,
					     &x.front(),
					     &m1.front(),
					     &e.front(),
					     &s1.front() );
      tmp2->SetMarkerStyle(20);
      tmp2->SetMarkerSize(1.5);
      tmp2->SetLineWidth(2);
      tmp2->SetTitle("");
      mg->Add(tmp2,"p");
      stringstream ss; ss << "Weighted average of all " << nhistos << " closure tests";
      leg->AddEntry(tmp2,ss.str().c_str(),"p");

      for ( uint ii = 0; ii < m1.size(); ++ii ) { 
	std::cout << std::fixed
		  << std::setprecision(0)
		  << " HT: " << bins[ii]
		  << std::setprecision(2)
		  << " mean: " << m1[ii]
		  << " std. dev.: " << s1[ii]
		  << " std. dev.: " << u1[ii]
		  << std::endl;
      }

      double fmean,ferr,fchi2; int fdof;
      fitResults(tmp2,fmean,ferr,fchi2,fdof,"");

    }

  }

  TCanvas* c = new TCanvas("tmp","tmp",900,600);
  mg->Draw("ap");
  mg->GetXaxis()->SetTitle("H_{T} (GeV)");
  mg->GetYaxis()->SetTitle("( N_{obs} - N_{pred} ) / N_{pred}");
  //mg->GetYaxis()->SetRangeUser(-1.25,2.75);
  mg->GetYaxis()->SetRangeUser(-2.,4.);
  mg->GetXaxis()->SetRangeUser(bins[0],bins[nbins]);
  mg->GetXaxis()->SetNdivisions(510);
  leg->Draw("same");
  std::stringstream sss;
  sss << "#splitline{CMS Preliminary}{L_{int} = " << lumi << " fb^{-1}, #sqrt{s} = " << energy << " TeV}";
  float tsize = 0.04;
  TLatex* tex = new TLatex(0.62,0.91-0.06,sss.str().c_str());
  tex->SetNDC();
  tex->SetTextSize(tsize);
  tex->Draw();

  std::stringstream sss1;
  if       ( njet == "le3" ) { sss1 << "2 #leq n_{jet} #leq 3"; }
  else if  ( njet == "ge4" ) { sss1 << "n_{jet} #geq 4"; }
  else if  ( njet == "ge2" ) { sss1 << "n_{jet} #geq 2"; }
  else                       { sss1 << ""; }
  TLatex* tex1 = new TLatex(0.62,0.76,sss1.str().c_str());
  tex1->SetNDC();
  tex1->SetTextSize(tsize);
  tex1->Draw();

  if ( !draw_syst ) { 
    TLine* line = new TLine(bins[0],0.,bins[nbins],0.); 
    line->SetLineColor(15);
    line->Draw(); 
  }
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
void extractNumbers( vstring& t, 
		     vvdouble& v, 
		     vvdouble& e, 
		     vdouble& n, 
		     vdouble& m, 
		     vdouble& s,
		     vdouble& u,
		     bool do_it ) {

  // Init
  m.clear(); m.resize(!v.empty()?v[0].size():0,0.);
  s.clear(); s.resize(!v.empty()?v[0].size():0,0.);
  u.clear(); u.resize(!v.empty()?v[0].size():0,0.);

  // Check sizes
  sort(n.begin(),n.end());
  if ( n.empty() || n.back() >= v[0].size() ) { cout << "Problem!" << endl; return; } 
  
  // Calc bands
  if ( syst_option == 0 ) {

    // weighted mean and variance

    vdouble numer; init(n.size(),numer);
    vdouble denom; init(n.size(),denom);
    for ( uint j = 0; j < n.size(); ++j ) {
      for ( uint k = n[j]; k < (j<n.size()-1?n[j+1]:v[0].size()); ++k ) {
	for ( uint i = 0; i < v.size(); ++i ) {
	  //cout << " test: " << i << " region: " << j << " bin: " << k << " val: " << v[i][k] << " err: " << e[i][k];
	  if ( e[i][k] > 0. && v[i][k] < 100. ) {
	    double val = v[i][k];
	    double err = e[i][k];
	    double w = weighted ? 1./(err*err) : 1.;
	    numer[j] += val*w;
	    denom[j] += w;
	    //cout << " w: " << w << " numer: " << numer[j] << " denom: " << denom[j];
	  }
	  //cout << endl;
	}
      }
    }
    for ( uint j = 0; j < n.size(); ++j ) {
      if ( denom[j] > 0 ) {
	m[j] = numer[j] / denom[j];
	s[j] = TMath::Sqrt( 1. / denom[j] ); 
      } else { m[j] = -1000.; s[j] = 0.; }
    }  
    
    float coverage = 1. - (1.-confidence_level)/2.;
    float sigma = ROOT::Math::normal_quantile(coverage,1.);
    std::cout << " CL: " << confidence_level 
	      << " coverage: " << coverage
	      << " sigma: " << sigma
	      << std::endl;
    
    s.resize(n.size(),0.);
    u.resize(n.size(),0.);
    for ( uint i = 0; i < n.size(); ++i ) { u[i] = sigma*s[i]; }

  } else if ( syst_option == 1 ) {

    // syst = sqrt( mean^2 + var )

    vdouble numer; init(n.size(),numer);
    vdouble denom; init(n.size(),denom);
    for ( uint j = 0; j < n.size(); ++j ) {
      for ( uint k = n[j]; k < (j<n.size()-1?n[j+1]:v[0].size()); ++k ) {
	for ( uint i = 0; i < v.size(); ++i ) {
	  //cout << " test: " << i << " region: " << j << " bin: " << k << " val: " << v[i][k] << " err: " << e[i][k];
	  if ( e[i][k] > 0. && v[i][k] < 100. ) {
	    double err = e[i][k];
	    double w = weighted ? 1./(err*err) : 1.;
	    numer[j] += v[i][k] * w;
	    denom[j] += w;
	    //cout << " w: " << w << " numer: " << numer[j] << " denom: " << denom[j];
	  }
	  cout << endl;
	}
      }
    }
    for ( uint j = 0; j < n.size(); ++j ) {
      if ( denom[j] > 0 ) {
	m[j] = numer[j] / denom[j];
	s[j] = TMath::Sqrt( 1. / denom[j] ); 
      } else { m[j] = -1000.; s[j] = 0.; }
    }  
    
    s.resize(n.size(),0.);
    u.resize(n.size(),0.);
    for ( uint i = 0; i < n.size(); ++i ) { u[i] = sqrt( m[i]*m[i]  + s[i]*s[i] ); }
    
  } else if ( syst_option == 3 ) {

    // paris

    vdouble numer; init(n.size(),numer);
    vdouble denom; init(n.size(),denom);
    for ( uint j = 0; j < n.size(); ++j ) {
      for ( uint k = n[j]; k < (j<n.size()-1?n[j+1]:v[0].size()); ++k ) {
	for ( uint i = 0; i < v.size(); ++i ) {
	  //cout << " test: " << i << " region: " << j << " bin: " << k << " val: " << v[i][k] << " err: " << e[i][k];
	  double err = sqrt( e[i][k]*e[i][k] + v[i][k]*v[i][k] );
	  if ( err > 0. && v[i][k] < 100. ) {
	    double w = weighted ? 1./(err*err) : 1.;
	    numer[j] += v[i][k] * w;
	    denom[j] += w;
	    //cout << " w: " << w << " numer: " << numer[j] << " denom: " << denom[j];
	  }
	  cout << endl;
	}
      }
    }
    for ( uint j = 0; j < n.size(); ++j ) {
      if ( denom[j] > 0 ) {
	m[j] = numer[j] / denom[j];
	s[j] = TMath::Sqrt( 1. / denom[j] ); 
      } else { m[j] = -1000.; s[j] = 0.; }
    }  

    float coverage = 1. - (1.-confidence_level)/2.;
    float sigma = ROOT::Math::normal_quantile(coverage,1.);
    std::cout << " CL: " << confidence_level 
	      << " coverage: " << coverage
	      << " sigma: " << sigma
	      << std::endl;
    
    s.resize(n.size(),0.);
    u.resize(n.size(),0.);
    for ( uint i = 0; i < n.size(); ++i ) { u[i] = sigma*s[i]; }
    
  } else if ( syst_option == 4 ) {

    // louis

    vdouble var; init(n.size(),var);
    vdouble num; init(n.size(),num);
    vdouble numer; init(n.size(),numer);
    vdouble denom; init(n.size(),denom);
    for ( uint j = 0; j < n.size(); ++j ) {
      for ( uint k = n[j]; k < (j<n.size()-1?n[j+1]:v[0].size()); ++k ) {
	for ( uint i = 0; i < v.size(); ++i ) {
	  //cout << " test: " << i << " region: " << j << " bin: " << k << " val: " << v[i][k] << " err: " << e[i][k];
	  double err = e[i][k];
	  if ( err > 0. && v[i][k] < 100. ) {
	    double w = weighted ? 1./(err*err) : 1.;
	    numer[j] += v[i][k] * w;
	    denom[j] += w;
	    //cout << " w: " << w << " numer: " << numer[j] << " denom: " << denom[j];
	  }
	  //cout << endl;
	  if ( v[i][k] > -999. ) {
	    double tmp = ( v[i][k]*v[i][k] - e[i][k]*e[i][k] );
	    var[j] += tmp;
	    num[j]++;
	    cout << " region: " << j
		 << " bin: " << k
		 << " test: " << i
		 << " v: " << v[i][k]
		 << " vv: " << v[i][k]*v[i][k]
		 << " e: " << e[i][k]
		 << " ee: " << e[i][k]*e[i][k]
		 << " diff: " << tmp
		 << " pos?: " << (tmp>0.?1.:0.)
		 << " var: " << var[j]
		 << endl;
	  }
	}
      }
    }
    for ( uint j = 0; j < n.size(); ++j ) {
      //if ( denom[j] > 0 ) {
      if ( num[j] > 0 ) {
	//m[j] = numer[j] / denom[j];
	s[j] = TMath::Sqrt( (var[j]>0.?var[j]:0.)/num[j] ); 
      } else { m[j] = -1000.; s[j] = 0.; }
    }  

    float coverage = 1. - (1.-confidence_level)/2.;
    float sigma = ROOT::Math::normal_quantile(coverage,1.);
    std::cout << " CL: " << confidence_level 
	      << " coverage: " << coverage
	      << " sigma: " << sigma
	      << std::endl;
    
    s.resize(n.size(),0.);
    u.resize(n.size(),0.);
    for ( uint i = 0; i < n.size(); ++i ) { u[i] = sigma*s[i]; }
    
  } else if ( syst_option == 5 ) {

    // sample variance

    // Calculate (normalised) weights
    vvdouble weights; init( v.size(), v[0].size(), weights );
    for ( uint j = 0; j < n.size(); ++j ) {
      double total = 0.;
      for ( uint k = n[j]; k < (j<n.size()-1?n[j+1]:v[0].size()); ++k ) {
	for ( uint i = 0; i < v.size(); ++i ) {
	  if ( v[i][k] < -999. ) { continue; }
	  double err = e[i][k]; //err = sqrt( v[i][k]*v[i][k] + e[i][k]*e[i][k] ); // paris
	  weights[i][k] = 1./(err*err);
	  total += weights[i][k];
	}
      }
      for ( uint k = n[j]; k < (j<n.size()-1?n[j+1]:v[0].size()); ++k ) {
	for ( uint i = 0; i < v.size(); ++i ) {
	  weights[i][k] /= total;
	}
      }
    }
    
    // Calculate weighted mean
    vdouble num; init(n.size(),num);
    vdouble means; init( n.size(), means );
    vdouble vars; init( n.size(), vars );
    for ( uint j = 0; j < n.size(); ++j ) {
      for ( uint k = n[j]; k < (j<n.size()-1?n[j+1]:v[0].size()); ++k ) {
	for ( uint i = 0; i < v.size(); ++i ) {
	  if ( v[i][k] < -999. ) { continue; }
	  means[j] += ( v[i][k] * weights[i][k] );
	  //vars[j] += ( weights[i][k] * ( v[i][k]*v[i][k] - e[i][k]*e[i][k] ) ); // louis
	  vars[j] += ( weights[i][k] * ( v[i][k] - means[j] ) * ( v[i][k] - means[j] ) );
	  num[j] += weights[i][k];
	}
      }
    }

    // Set values
    for ( uint j = 0; j < n.size(); ++j ) {
      m[j] = means[j]/num[j];
      //s[j] = TMath::Sqrt( (vars[j]>0?vars[j]:0.)/num[j] ); // louis
      s[j] = TMath::Sqrt( vars[j]/num[j] ); 
    }  

    float coverage = 1. - (1.-confidence_level)/2.;
    float sigma = ROOT::Math::normal_quantile(coverage,1.);
    std::cout << " CL: " << confidence_level 
	      << " coverage: " << coverage
	      << " sigma: " << sigma
	      << std::endl;
    
    s.resize(n.size(),0.);
    u.resize(n.size(),0.);
    for ( uint i = 0; i < n.size(); ++i ) { u[i] = sigma*s[i]; }
    //for ( uint i = 0; i < n.size(); ++i ) { u[i] = sqrt( m[i]*m[i] + s[i]*s[i] ); } // paris
    
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
	if ( tmp3 > 0.5 )  index++;
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

    
  } //else if ( syst_option == 6 && do_it ) 
  
  if ( do_it ) {

    // Use Gaussian fits
    
    std::vector<TF1*> fit(n.size(),0);
    std::vector<TCanvas*> c(n.size(),0);
    std::vector<TLegend*> leg(n.size(),0);
    for ( uint j = 0; j < n.size(); ++j ) { // regions
      
      leg[j] = new TLegend( 0.5, 0.82-0.04*v.size(), 0.6, 0.82 );
      leg[j]->SetFillColor(0);
      leg[j]->SetLineColor(0); 
      leg[j]->SetShadowColor(0); 
      leg[j]->SetTextSize(0.03);
      
      std::vector<TH1F*> h(v.size(),0);
      for ( uint i = 0; i < v.size(); ++i ) { // Test
	
	h[i] = new TH1F("","",500,-1.,4.); 
	h[i]->Sumw2();
	
	for ( uint k = n[j]; k < (j<n.size()-1?n[j+1]:v[0].size()); ++k ) { // HT bins
	  if ( v[i][k] < -999. ) { continue; }
	  h[i]->Fill( v[i][k], 1.);///(e[i][k]*e[i][k]) );
	}
	
	h[i]->SetLineColor(i+2);
	h[i]->SetFillColor(i+2);
	if ( i ) { h[i]->Add( h[i-1] ); }
	leg[j]->AddEntry(h[i],std::string("  "+t[i]).c_str(),"f");
	
      }
      
      std::stringstream ss; ss << "Region" << j;
      c[j] = new TCanvas(ss.str().c_str(),ss.str().c_str(),600,600);
      for ( uint i = v.size(); i > 0; --i ) { 
	if ( i == v.size() ) { 
	  h[i-1]->Draw("hist"); 

	  float band = u[j]; //h[i-1]->GetRMS();
	  TBox* box = new TBox( -1.*band, 0., 1.*band, 1.2*h[i-1]->GetMaximum() );
	  box->SetFillColor(kGray);
	  box->SetLineColor(18);
	  box->Draw("same");
	  h[i-1]->Draw("hist same"); 
	  h[i-1]->SetMaximum( 1.2*h[i-1]->GetMaximum() );
	  h[i-1]->GetXaxis()->SetNdivisions(1005); 
	  h[i-1]->GetYaxis()->SetNdivisions(505); 
	  h[i-1]->GetXaxis()->SetTitle("( N_{obs} - N_{pred} ) / N_{pred}"); 
	  //h[i-1]->GetYaxis()->SetTitle("Number of tests");// ("a.u."); 
	  //h[i-1]->GetYaxis()->SetTitle("#Sigma_{i}^{N} 1/#sigma_{i}^{2}"); 
	  h[i-1]->GetYaxis()->SetTitle("#Sigma 1/#sigma^{2}"); 
	}
	else { h[i-1]->Draw("hist same"); }

      }

      c[j]->Draw();
      leg[j]->Draw("same");
      
      TH1F* tmp = new TH1F(*h.back());
      fit[j] = new TF1("fit","gaus",tmp->GetXaxis()->GetXmin(),tmp->GetXaxis()->GetXmax()); 
      fit[j]->FixParameter( 1, 0. );
      fit[j]->FixParameter( 2, u[j] );
      //fit[j]->SetParLimits(0,0.,5.*tmp->GetMaximum());
      //fit[j]->SetParLimits(1,-5.*tmp->GetMean(),5.*tmp->GetMean());
      //fit[j]->SetParLimits(2,0.,5.*tmp->GetRMS());
      tmp->Fit(fit[j],"NB");
      fit[j]->SetLineColor(1);
      fit[j]->SetLineWidth(1.5);
      //tmp->SetFillColor(1);
      //tmp->Draw("hist same");
      fit[j]->Draw("c same");

      {
	std::stringstream sss;
	sss << "CMS Preliminary, L_{int} = " << lumi << " fb^{-1}, #sqrt{s} = " << energy << " TeV";
	TLatex* tex = new TLatex(0.4,0.90,sss.str().c_str());
	tex->SetNDC();
	tex->SetTextSize(0.03);
	tex->Draw();
      }

      {
	std::stringstream sss;
	int low = n[j];
	int high = j < n.size() - 1 ? n[j+1] : v[0].size();
	sss << bins[low] << "< H_{T} < " << bins[high] << " GeV";
	TLatex* tex = new TLatex(0.5,0.84,sss.str().c_str());
	tex->SetNDC();
	tex->SetTextSize(0.03);
	tex->Draw();
      }

      {
	std::stringstream sss;
	sss << "His mean: " 
	    << std::fixed
	    << std::setprecision(2) 
	    << tmp->GetMean()
	    << " #pm " 
	    << tmp->GetMeanError()
	    << " rms: " 
	    << tmp->GetRMS()
	    << " #pm "
	    << tmp->GetRMSError();

	TLatex* tex = new TLatex(0.4,0.40,sss.str().c_str());
	tex->SetNDC();
	tex->SetTextSize(0.03);
	//tex->Draw();
      }

      {
	std::stringstream sss;
	sss << "Fit mean: " 
	    << std::fixed
	    << std::setprecision(2) 
	    << fit[j]->GetParameter(1)
	    << " #pm " 
	    << fit[j]->GetParError(1)
	    << " sigma: " 
	    << fit[j]->GetParameter(2)
	    << " #pm " 
	    << fit[j]->GetParError(2);

	TLatex* tex = new TLatex(0.4,0.36,sss.str().c_str());
	tex->SetNDC();
	tex->SetTextSize(0.03);
	//tex->Draw();
      }

      gPad->RedrawAxis();

      //m[j] = fabs(fit[j]->GetParameter(1));
      //s[j] = fabs(fit[j]->GetParameter(2));
//       m[j] = fabs(tmp->GetMean());
//       s[j] = fabs(tmp->GetRMS());
      
      
      std::stringstream ssss;
      ssss << c[j]->GetTitle() << ".pdf";
      cout << ssss.str() << endl;
      c[j]->SaveAs(ssss.str().c_str());

    }

    float coverage = 1. - (1.-confidence_level)/2.;
    float sigma = ROOT::Math::normal_quantile(coverage,1.);
    std::cout << " CL: " << confidence_level 
	      << " coverage: " << coverage
	      << " sigma: " << sigma
	      << std::endl;
    
    s.resize(n.size(),0.);
    u.resize(n.size(),0.);
    for ( uint i = 0; i < n.size(); ++i ) { u[i] = sigma*s[i]; }

//   } else {
//     std::cout << "UNKNOWN OPTION!!!" << std::endl;
   }
  
//   for ( uint j = 0; j < n.size(); ++j ) {
//     std::cout << " bin: " << j
// 	      << " HT: " << bins[n[j]]
// 	      << " sum: " << numer[j]
// 	      << " sum2: " << numer2[j]
// 	      << " n: " << denom[j]
// 	      << " mean: " << m[j]
// 	      << " err: " << s[j]
// 	      << " uncert: " << u[j]
// 	      << std::endl;
//   }
  
    for ( uint ii = 0; ii < n.size(); ++ii ) { 
      std::cout	<< " n: " << n[ii]
		<< " m: " << m[ii]
		<< " s: " << s[ii]
		<< " u: " << u[ii]
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
void fitResults( TGraphErrors* tmp, double& tmp1, double& tmp2, double& tmp3, int& tmp4, std::string title ) {

  TGraphErrors* temp0 = new TGraphErrors( *tmp );
  TF1* fit0 = 0;

  TGraphErrors* temp1 = new TGraphErrors( *tmp );
  TF1* fit1 = 0;

  double x0 = 0., y0 = 0.; (((TGraphErrors*)(temp0))->GetPoint(0,x0,y0));
  if ( x0 < bins[0]-5. ) { fit0 = new TF1("fit0","pol0",bins[2]-5.,bins[nbins]+5.); }
  else                   { fit0 = new TF1("fit0","pol0",bins[0]-5.,bins[nbins]+5.); }
  temp0->Fit(fit0,"QR");

  double x1 = 0., y1 = 0.; (((TGraphErrors*)(temp1))->GetPoint(0,x1,y1));
  if ( x1 < bins[0]-5. ) { fit1 = new TF1("fit1","pol1",bins[2]-5.,bins[nbins]+5.); }
  else                   { fit1 = new TF1("fit1","pol1",bins[0]-5.,bins[nbins]+5.); }
  temp1->Fit(fit1,"QR");

  double a275 = fit1->GetParameter(0) + 275.*fit1->GetParameter(1);
  double e275 = fit1->GetParameter(0) + 275.*fit1->GetParameter(1);
  

//   \begin{tabular}{ ll|rccc|rc }
//     \hline
//                                    &                   & \multicolumn{4}{c|}{Constant fit} & \multicolumn{2}{c}{Linear fit}                                \\
//     Sample                         & Symbol            & Best fit value                    & $\chi^2$ & d.o.f. & $p$-value & Slope ($10^{-4}$) & $p$-value \\
//     \hline                                                                                                            
//     \mj                            & Inverted triangle & $-0.03 \pm 0.02$                  & 17.3     & 7      & 0.02      & $0.0 \pm 1.0$     & 0.01      \\ 
//     \mj (outlier removed)$^{\dag}$ & Inverted triangle & $-0.04 \pm 0.01$                  & 6.1      & 6      & 0.42      & $-1.4 \pm 1.1$    & 0.49      \\ 
//     \gj                            & Diamond           & $ 0.12 \pm 0.05$                  & 2.42     & 5      & 0.79      & $6.0 \pm 4.7$     & 0.94      \\ 
//     \mmj                           & Asterisk          & $-0.04 \pm 0.07$                  & 9.76     & 7      & 0.20      & $4.9 \pm 4.4$     & 0.20      \\ 
//     \hline


  std::cout << "& $" 
	    << std::setprecision(2) << fit0->GetParameter(0)
	    << " \\pm " 
	    << std::setprecision(2) << fit0->GetParError(0)
	    << "$ & " 
	    << std::setprecision(3) << fit0->GetChisquare() 
	    << " & " 
	    << std::setprecision(0) << fit0->GetNDF() 
	    << " & " 
	    << std::setprecision(2) << fit0->GetProb() 
	    << " & $"
	    << std::setprecision(2) << fit1->GetParameter(1)
	    << " \\pm " 
	    << std::setprecision(2) << fit1->GetParError(1)
	    << "$ & " 
	    << std::setprecision(3) << fit1->GetChisquare() 
	    << " & " 
	    << std::setprecision(0) << fit1->GetNDF() 
	    << " & " 
	    << std::setprecision(2) << fit1->GetProb() 
	    << " \\\\"
	    << std::endl;

  tmp1 = fit0->GetParameter(0);
  tmp2 = fit0->GetParError(0);
  tmp3 = fit0->GetChisquare();
  tmp4 = fit0->GetNDF();
  
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
  
  if ( set_option == 100 ) {
    
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
    
  } else if ( set_option == 0 ) {
    
    // ICHEP PAS

    energy = "8";
    lumi = "3.9";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.20);
    syst.push_back(0.30);
    syst.push_back(0.70);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
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
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
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
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
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
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
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
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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

    for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 1 ) {

    // 5/fb "top up", inclusive

    energy = "8";
    lumi = "5.0";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.20);
    syst.push_back(0.30);
    syst.push_back(0.40);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.0292667);
   grae->SetPointError(1,0,0,0.06189248,0.06189248);
   grae->SetPoint(2,350,-0.07421998);
   grae->SetPointError(2,0,0,0.0539984,0.0539984);
   grae->SetPoint(3,425,-0.06773209);
   grae->SetPointError(3,0,0,0.05558163,0.05558163);
   grae->SetPoint(4,525,-0.1157101);
   grae->SetPointError(4,0,0,0.09051336,0.09051336);
   grae->SetPoint(5,625,-0.2348931);
   grae->SetPointError(5,0,0,0.1530066,0.1530066);
   grae->SetPoint(6,725,-0.208552);
   grae->SetPointError(6,0,0,0.2746149,0.2746149);
   grae->SetPoint(7,825,0.04282867);
   grae->SetPointError(7,0,0,0.362482,0.362482);
   grae->SetPoint(8,925,-0.8292535);
   grae->SetPointError(8,0,0,0.5857774,0.5857774);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1928643);
   grae->SetPointError(1,0,0,0.06223478,0.06223478);
   grae->SetPoint(2,350,0.1118001);
   grae->SetPointError(2,0,0,0.05573281,0.05573281);
   grae->SetPoint(3,425,0.1838906);
   grae->SetPointError(3,0,0,0.05707422,0.05707422);
   grae->SetPoint(4,525,0.1420095);
   grae->SetPointError(4,0,0,0.07725208,0.07725208);
   grae->SetPoint(5,625,0.1325877);
   grae->SetPointError(5,0,0,0.1188084,0.1188084);
   grae->SetPoint(6,725,0.4611449);
   grae->SetPointError(6,0,0,0.1993745,0.1993745);
   grae->SetPoint(7,825,0.1787153);
   grae->SetPointError(7,0,0,0.2478256,0.2478256);
   grae->SetPoint(8,925,-0.2787331);
   grae->SetPointError(8,0,0,0.2551869,0.2551869);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.01685955);
   grae->SetPointError(1,0,0,0.07202448,0.07202448);
   grae->SetPoint(2,350,-0.03066361);
   grae->SetPointError(2,0,0,0.07974449,0.07974449);
   grae->SetPoint(3,425,-0.06451732);
   grae->SetPointError(3,0,0,0.07624193,0.07624193);
   grae->SetPoint(4,525,0.02759051);
   grae->SetPointError(4,0,0,0.1076957,0.1076957);
   grae->SetPoint(5,625,-0.0523672);
   grae->SetPointError(5,0,0,0.1631205,0.1631205);
   grae->SetPoint(6,725,-0.1871218);
   grae->SetPointError(6,0,0,0.2370259,0.2370259);
   grae->SetPoint(7,825,-0.07837915);
   grae->SetPointError(7,0,0,0.3342762,0.3342762);
   grae->SetPoint(8,925,0.03482808);
   grae->SetPointError(8,0,0,0.4127586,0.4127586);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.09870746);
   grae->SetPointError(1,0,0,0.06552383,0.06552383);
   grae->SetPoint(2,350,-0.03546442);
   grae->SetPointError(2,0,0,0.08323372,0.08323372);
   grae->SetPoint(3,425,0.06143153);
   grae->SetPointError(3,0,0,0.08481386,0.08481386);
   grae->SetPoint(4,525,0.01069821);
   grae->SetPointError(4,0,0,0.1337018,0.1337018);
   grae->SetPoint(5,625,0.249106);
   grae->SetPointError(5,0,0,0.1901936,0.1901936);
   grae->SetPoint(6,725,0.0490647);
   grae->SetPointError(6,0,0,0.2675704,0.2675704);
   grae->SetPoint(7,825,0.1958008);
   grae->SetPointError(7,0,0,0.3804563,0.3804563);
   grae->SetPoint(8,925,0.1082358);
   grae->SetPointError(8,0,0,0.3744375,0.3744375);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,-0.04030034);
   grae->SetPointError(3,0,0,0.08787241,0.08787241);
   grae->SetPoint(4,525,-0.1275339);
   grae->SetPointError(4,0,0,0.1401435,0.1401435);
   grae->SetPoint(5,625,0.03829255);
   grae->SetPointError(5,0,0,0.1960235,0.1960235);
   grae->SetPoint(6,725,-0.2289197);
   grae->SetPointError(6,0,0,0.2915995,0.2915995);
   grae->SetPoint(7,825,-0.1800639);
   grae->SetPointError(7,0,0,0.4032558,0.4032558);
   grae->SetPoint(8,925,0.3705466);
   grae->SetPointError(8,0,0,0.5012573,0.5012573);

    for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 2 ) {

    // 5/fb "top up", 2-3 jets

    energy = "8";
    lumi = "5.0";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.25);
    syst.push_back(0.45);
    syst.push_back(0.55);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.05593175);
   grae->SetPointError(1,0,0,0.07519763,0.07519763);
   grae->SetPoint(2,350,-0.04975588);
   grae->SetPointError(2,0,0,0.05732012,0.05732012);
   grae->SetPoint(3,425,-0.05126801);
   grae->SetPointError(3,0,0,0.06239574,0.06239574);
   grae->SetPoint(4,525,-0.1605633);
   grae->SetPointError(4,0,0,0.1104441,0.1104441);
   grae->SetPoint(5,625,-0.1694169);
   grae->SetPointError(5,0,0,0.1924076,0.1924076);
   grae->SetPoint(6,725,-0.232026);
   grae->SetPointError(6,0,0,0.3729678,0.3729678);
   grae->SetPoint(7,825,-0.07911233);
   grae->SetPointError(7,0,0,0.4718303,0.4718303);
   grae->SetPoint(8,925,-1);
   grae->SetPointError(8,0,0,0.9148121,0.9148121);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.2283006);
   grae->SetPointError(1,0,0,0.08386912,0.08386912);
   grae->SetPoint(2,350,0.1106251);
   grae->SetPointError(2,0,0,0.06645799,0.06645799);
   grae->SetPoint(3,425,0.2564887);
   grae->SetPointError(3,0,0,0.07095924,0.07095924);
   grae->SetPoint(4,525,0.2480567);
   grae->SetPointError(4,0,0,0.09946865,0.09946865);
   grae->SetPoint(5,625,0.1889247);
   grae->SetPointError(5,0,0,0.1598149,0.1598149);
   grae->SetPoint(6,725,0.2325977);
   grae->SetPointError(6,0,0,0.263721,0.263721);
   grae->SetPoint(7,825,0.1812302);
   grae->SetPointError(7,0,0,0.3298413,0.3298413);
   grae->SetPoint(8,925,-0.5321255);
   grae->SetPointError(8,0,0,0.4564076,0.4564076);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.03168183);
   grae->SetPointError(1,0,0,0.1012294,0.1012294);
   grae->SetPoint(2,350,-0.06112091);
   grae->SetPointError(2,0,0,0.1039506,0.1039506);
   grae->SetPoint(3,425,-0.112151);
   grae->SetPointError(3,0,0,0.1000182,0.1000182);
   grae->SetPoint(4,525,-0.03154164);
   grae->SetPointError(4,0,0,0.1621475,0.1621475);
   grae->SetPoint(5,625,-0.2880484);
   grae->SetPointError(5,0,0,0.2991294,0.2991294);
   grae->SetPoint(6,725,-0.4879241);
   grae->SetPointError(6,0,0,0.5440645,0.5440645);
   grae->SetPoint(7,825,0.04123313);
   grae->SetPointError(7,0,0,0.6156035,0.6156035);
   grae->SetPoint(8,925,0.3493286);
   grae->SetPointError(8,0,0,1.022466,1.022466);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1373553);
   grae->SetPointError(1,0,0,0.07067155,0.07067155);
   grae->SetPoint(2,350,0.04466097);
   grae->SetPointError(2,0,0,0.08831794,0.08831794);
   grae->SetPoint(3,425,0.06311979);
   grae->SetPointError(3,0,0,0.08989551,0.08989551);
   grae->SetPoint(4,525,-0.08023821);
   grae->SetPointError(4,0,0,0.1503083,0.1503083);
   grae->SetPoint(5,625,0.180593);
   grae->SetPointError(5,0,0,0.2119664,0.2119664);
   grae->SetPoint(6,725,0.4398447);
   grae->SetPointError(6,0,0,0.3763019,0.3763019);
   grae->SetPoint(7,825,0.3246639);
   grae->SetPointError(7,0,0,0.4672136,0.4672136);
   grae->SetPoint(8,925,0.1115713);
   grae->SetPointError(8,0,0,0.4565074,0.4565074);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,-0.02059047);
   grae->SetPointError(3,0,0,0.09320259,0.09320259);
   grae->SetPoint(4,525,-0.1297196);
   grae->SetPointError(4,0,0,0.1597241,0.1597241);
   grae->SetPoint(5,625,-0.03832861);
   grae->SetPointError(5,0,0,0.222412,0.222412);
   grae->SetPoint(6,725,-0.1812387);
   grae->SetPointError(6,0,0,0.3610356,0.3610356);
   grae->SetPoint(7,825,0.05063377);
   grae->SetPointError(7,0,0,0.4993848,0.4993848);
   grae->SetPoint(8,925,0.1890704);
   grae->SetPointError(8,0,0,0.5764542,0.5764542);

    for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 3 ) {

    // 5/fb "top up", >3 jets

    energy = "8";
    lumi = "5.0";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.40);
    syst.push_back(0.50);
    syst.push_back(0.70);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.05958225);
   grae->SetPointError(1,0,0,0.08010016,0.08010016);
   grae->SetPoint(2,350,-0.1597706);
   grae->SetPointError(2,0,0,0.1334475,0.1334475);
   grae->SetPoint(3,425,-0.1141818);
   grae->SetPointError(3,0,0,0.1164997,0.1164997);
   grae->SetPoint(4,525,-0.03340484);
   grae->SetPointError(4,0,0,0.1542016,0.1542016);
   grae->SetPoint(5,625,-0.3166521);
   grae->SetPointError(5,0,0,0.2477168,0.2477168);
   grae->SetPoint(6,725,-0.2364497);
   grae->SetPointError(6,0,0,0.3987308,0.3987308);
   grae->SetPoint(7,825,0.200159);
   grae->SetPointError(7,0,0,0.576704,0.576704);
   grae->SetPoint(8,925,-0.6932969);
   grae->SetPointError(8,0,0,0.754413,0.754413);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.06638716);
   grae->SetPointError(1,0,0,0.1020102,0.1020102);
   grae->SetPoint(2,350,0.146388);
   grae->SetPointError(2,0,0,0.1147124,0.1147124);
   grae->SetPoint(3,425,0.08264517);
   grae->SetPointError(3,0,0,0.1055916,0.1055916);
   grae->SetPoint(4,525,0.1404841);
   grae->SetPointError(4,0,0,0.1291087,0.1291087);
   grae->SetPoint(5,625,0.1018502);
   grae->SetPointError(5,0,0,0.1803489,0.1803489);
   grae->SetPoint(6,725,0.4767506);
   grae->SetPointError(6,0,0,0.2841773,0.2841773);
   grae->SetPoint(7,825,0.2193412);
   grae->SetPointError(7,0,0,0.3762187,0.3762187);
   grae->SetPoint(8,925,-0.1842508);
   grae->SetPointError(8,0,0,0.308465,0.308465);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.01630152);
   grae->SetPointError(1,0,0,0.09007696,0.09007696);
   grae->SetPoint(2,350,0.008089063);
   grae->SetPointError(2,0,0,0.1238612,0.1238612);
   grae->SetPoint(3,425,0.04414519);
   grae->SetPointError(3,0,0,0.1188969,0.1188969);
   grae->SetPoint(4,525,0.1160613);
   grae->SetPointError(4,0,0,0.1476922,0.1476922);
   grae->SetPoint(5,625,0.03794876);
   grae->SetPointError(5,0,0,0.2048531,0.2048531);
   grae->SetPoint(6,725,-0.1795167);
   grae->SetPointError(6,0,0,0.278857,0.278857);
   grae->SetPoint(7,825,-0.1043569);
   grae->SetPointError(7,0,0,0.4158672,0.4158672);
   grae->SetPoint(8,925,-0.1726461);
   grae->SetPointError(8,0,0,0.4391614,0.4391614);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.1115897);
   grae->SetPointError(1,0,0,0.1766693,0.1766693);
   grae->SetPoint(2,350,-0.6058388);
   grae->SetPointError(2,0,0,0.2873025,0.2873025);
   grae->SetPoint(3,425,0.01192187);
   grae->SetPointError(3,0,0,0.2458318,0.2458318);
   grae->SetPoint(4,525,0.4056748);
   grae->SetPointError(4,0,0,0.3098359,0.3098359);
   grae->SetPoint(5,625,0.6015155);
   grae->SetPointError(5,0,0,0.461605,0.461605);
   grae->SetPoint(6,725,-0.4463881);
   grae->SetPointError(6,0,0,0.4707705,0.4707705);
   grae->SetPoint(7,825,-0.177848);
   grae->SetPointError(7,0,0,0.6908435,0.6908435);
   grae->SetPoint(8,925,0.233931);
   grae->SetPointError(8,0,0,0.6445523,0.6445523);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,-0.2017726);
   grae->SetPointError(3,0,0,0.256438,0.256438);
   grae->SetPoint(4,525,-0.08462917);
   grae->SetPointError(4,0,0,0.2826813,0.2826813);
   grae->SetPoint(5,625,0.3841786);
   grae->SetPointError(5,0,0,0.4430207,0.4430207);
   grae->SetPoint(6,725,-0.4240053);
   grae->SetPointError(6,0,0,0.5214629,0.5214629);
   grae->SetPoint(7,825,-0.5540355);
   grae->SetPointError(7,0,0,0.746766,0.746766);
   grae->SetPoint(8,925,0.9160421);
   grae->SetPointError(8,0,0,1.140142,1.140142);

    for ( uint i = 0; i < nbins; ++i ) {
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

    // 8/fb "top up", inclusive

    energy = "8";
    lumi = "8.0";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.15);
    syst.push_back(0.30);
    syst.push_back(0.90);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.0427723);
   grae->SetPointError(1,0,0,0.05972027,0.05972027);
   grae->SetPoint(2,350,-0.08086641);
   grae->SetPointError(2,0,0,0.04796461,0.04796461);
   grae->SetPoint(3,425,-0.06679328);
   grae->SetPointError(3,0,0,0.04832904,0.04832904);
   grae->SetPoint(4,525,-0.09650421);
   grae->SetPointError(4,0,0,0.07851997,0.07851997);
   grae->SetPoint(5,625,-0.2620246);
   grae->SetPointError(5,0,0,0.131442,0.131442);
   grae->SetPoint(6,725,-0.13468);
   grae->SetPointError(6,0,0,0.2410398,0.2410398);
   grae->SetPoint(7,825,-0.08330379);
   grae->SetPointError(7,0,0,0.3198348,0.3198348);
   grae->SetPoint(8,925,-0.8791075);
   grae->SetPointError(8,0,0,0.518717,0.518717);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1551209);
   grae->SetPointError(1,0,0,0.05930979,0.05930979);
   grae->SetPoint(2,350,0.1198618);
   grae->SetPointError(2,0,0,0.050379,0.050379);
   grae->SetPoint(3,425,0.1534372);
   grae->SetPointError(3,0,0,0.05108691,0.05108691);
   grae->SetPoint(4,525,0.1203351);
   grae->SetPointError(4,0,0,0.06795124,0.06795124);
   grae->SetPoint(5,625,0.1015318);
   grae->SetPointError(5,0,0,0.1033242,0.1033242);
   grae->SetPoint(6,725,0.3257863);
   grae->SetPointError(6,0,0,0.1660006,0.1660006);
   grae->SetPoint(7,825,0.2755312);
   grae->SetPointError(7,0,0,0.2312355,0.2312355);
   grae->SetPoint(8,925,-0.2226956);
   grae->SetPointError(8,0,0,0.2313517,0.2313517);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.01273126);
   grae->SetPointError(1,0,0,0.06786036,0.06786036);
   grae->SetPoint(2,350,-0.01046181);
   grae->SetPointError(2,0,0,0.0712126,0.0712126);
   grae->SetPoint(3,425,-0.008207413);
   grae->SetPointError(3,0,0,0.06790508,0.06790508);
   grae->SetPoint(4,525,0.04967546);
   grae->SetPointError(4,0,0,0.0952294,0.0952294);
   grae->SetPoint(5,625,-0.03662366);
   grae->SetPointError(5,0,0,0.1418481,0.1418481);
   grae->SetPoint(6,725,-0.2132366);
   grae->SetPointError(6,0,0,0.2142033,0.2142033);
   grae->SetPoint(7,825,-0.1316644);
   grae->SetPointError(7,0,0,0.2971728,0.2971728);
   grae->SetPoint(8,925,-0.1320009);
   grae->SetPointError(8,0,0,0.3585166,0.3585166);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.07356922);
   grae->SetPointError(1,0,0,0.06111266,0.06111266);
   grae->SetPoint(2,350,-0.02891513);
   grae->SetPointError(2,0,0,0.07725978,0.07725978);
   grae->SetPoint(3,425,0.05603726);
   grae->SetPointError(3,0,0,0.07850519,0.07850519);
   grae->SetPoint(4,525,0.0212776);
   grae->SetPointError(4,0,0,0.1240002,0.1240002);
   grae->SetPoint(5,625,0.2172916);
   grae->SetPointError(5,0,0,0.1730176,0.1730176);
   grae->SetPoint(6,725,-0.07215984);
   grae->SetPointError(6,0,0,0.2458828,0.2458828);
   grae->SetPoint(7,825,0.1860796);
   grae->SetPointError(7,0,0,0.3508879,0.3508879);
   grae->SetPoint(8,925,0.2918501);
   grae->SetPointError(8,0,0,0.3575249,0.3575249);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,0.2026566);
   grae->SetPointError(3,0,0,0.08733804,0.08733804);
   grae->SetPoint(4,525,0.09605987);
   grae->SetPointError(4,0,0,0.1361478,0.1361478);
   grae->SetPoint(5,625,0.3271275);
   grae->SetPointError(5,0,0,0.2015542,0.2015542);
   grae->SetPoint(6,725,-0.1636186);
   grae->SetPointError(6,0,0,0.280972,0.280972);
   grae->SetPoint(7,825,-0.02011504);
   grae->SetPointError(7,0,0,0.392038,0.392038);
   grae->SetPoint(8,925,0.8790799);
   grae->SetPointError(8,0,0,0.6063338,0.6063338);

    for ( uint i = 0; i < nbins; ++i ) {
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

    // 8/fb "top up", 2-3 jets

    energy = "8";
    lumi = "8.0";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.15);
    syst.push_back(0.30);
    syst.push_back(0.90);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.06774757);
   grae->SetPointError(1,0,0,0.07286714,0.07286714);
   grae->SetPoint(2,350,-0.0521656);
   grae->SetPointError(2,0,0,0.05001001,0.05001001);
   grae->SetPoint(3,425,-0.04557748);
   grae->SetPointError(3,0,0,0.05391019,0.05391019);
   grae->SetPoint(4,525,-0.1093819);
   grae->SetPointError(4,0,0,0.09419403,0.09419403);
   grae->SetPoint(5,625,-0.2451693);
   grae->SetPointError(5,0,0,0.1624391,0.1624391);
   grae->SetPoint(6,725,-0.4767087);
   grae->SetPointError(6,0,0,0.3413534,0.3413534);
   grae->SetPoint(7,825,-0.2983007);
   grae->SetPointError(7,0,0,0.4364511,0.4364511);
   grae->SetPoint(8,925,-1);
   grae->SetPointError(8,0,0,0.7821514,0.7821514);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1792736);
   grae->SetPointError(1,0,0,0.08029543,0.08029543);
   grae->SetPoint(2,350,0.1219996);
   grae->SetPointError(2,0,0,0.06037335,0.06037335);
   grae->SetPoint(3,425,0.2294154);
   grae->SetPointError(3,0,0,0.06366705,0.06366705);
   grae->SetPoint(4,525,0.1823164);
   grae->SetPointError(4,0,0,0.08408354,0.08408354);
   grae->SetPoint(5,625,0.2586192);
   grae->SetPointError(5,0,0,0.1360979,0.1360979);
   grae->SetPoint(6,725,0.3275219);
   grae->SetPointError(6,0,0,0.2307485,0.2307485);
   grae->SetPoint(7,825,0.2725715);
   grae->SetPointError(7,0,0,0.3007453,0.3007453);
   grae->SetPoint(8,925,-0.4187433);
   grae->SetPointError(8,0,0,0.4086704,0.4086704);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.0004741051);
   grae->SetPointError(1,0,0,0.09589927,0.09589927);
   grae->SetPoint(2,350,-0.06700659);
   grae->SetPointError(2,0,0,0.09252062,0.09252062);
   grae->SetPoint(3,425,-0.06471827);
   grae->SetPointError(3,0,0,0.08841442,0.08841442);
   grae->SetPoint(4,525,-0.02859128);
   grae->SetPointError(4,0,0,0.1417515,0.1417515);
   grae->SetPoint(5,625,-0.1554906);
   grae->SetPointError(5,0,0,0.239135,0.239135);
   grae->SetPoint(6,725,-0.2250501);
   grae->SetPointError(6,0,0,0.4156814,0.4156814);
   grae->SetPoint(7,825,-0.2560571);
   grae->SetPointError(7,0,0,0.5510231,0.5510231);
   grae->SetPoint(8,925,-0.2290178);
   grae->SetPointError(8,0,0,0.7873435,0.7873435);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1040311);
   grae->SetPointError(1,0,0,0.06567307,0.06567307);
   grae->SetPoint(2,350,0.0422922);
   grae->SetPointError(2,0,0,0.08190388,0.08190388);
   grae->SetPoint(3,425,0.06518751);
   grae->SetPointError(3,0,0,0.08341801,0.08341801);
   grae->SetPoint(4,525,-0.05500957);
   grae->SetPointError(4,0,0,0.1397402,0.1397402);
   grae->SetPoint(5,625,0.1378851);
   grae->SetPointError(5,0,0,0.1926447,0.1926447);
   grae->SetPoint(6,725,0.2453775);
   grae->SetPointError(6,0,0,0.3238319,0.3238319);
   grae->SetPoint(7,825,0.2981313);
   grae->SetPointError(7,0,0,0.4334659,0.4334659);
   grae->SetPoint(8,925,0.296861);
   grae->SetPointError(8,0,0,0.4363215,0.4363215);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,0.2304226);
   grae->SetPointError(3,0,0,0.09326034,0.09326034);
   grae->SetPoint(4,525,0.1002445);
   grae->SetPointError(4,0,0,0.1551307,0.1551307);
   grae->SetPoint(5,625,0.2164782);
   grae->SetPointError(5,0,0,0.2224116,0.2224116);
   grae->SetPoint(6,725,-0.1347363);
   grae->SetPointError(6,0,0,0.3508925,0.3508925);
   grae->SetPoint(7,825,0.1844014);
   grae->SetPointError(7,0,0,0.5004343,0.5004343);
   grae->SetPoint(8,925,0.6497535);
   grae->SetPointError(8,0,0,0.6577595,0.6577595);

    for ( uint i = 0; i < nbins; ++i ) {
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

    // 8/fb "top up", >= 3 jets

    energy = "8";
    lumi = "8.0";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.15);
    syst.push_back(0.30);
    syst.push_back(0.90);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.04686251);
   grae->SetPointError(1,0,0,0.07291975,0.07291975);
   grae->SetPoint(2,350,-0.1801224);
   grae->SetPointError(2,0,0,0.122429,0.122429);
   grae->SetPoint(3,425,-0.1295565);
   grae->SetPointError(3,0,0,0.1018122,0.1018122);
   grae->SetPoint(4,525,-0.07061857);
   grae->SetPointError(4,0,0,0.1354875,0.1354875);
   grae->SetPoint(5,625,-0.2843153);
   grae->SetPointError(5,0,0,0.2148417,0.2148417);
   grae->SetPoint(6,725,0.03176696);
   grae->SetPointError(6,0,0,0.3550344,0.3550344);
   grae->SetPoint(7,825,0.1671163);
   grae->SetPointError(7,0,0,0.4985956,0.4985956);
   grae->SetPoint(8,925,-0.7796283);
   grae->SetPointError(8,0,0,0.6878265,0.6878265);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.09299544);
   grae->SetPointError(1,0,0,0.09711739,0.09711739);
   grae->SetPoint(2,350,0.1301034);
   grae->SetPointError(2,0,0,0.1013426,0.1013426);
   grae->SetPoint(3,425,0.01832605);
   grae->SetPointError(3,0,0,0.0915348,0.0915348);
   grae->SetPoint(4,525,0.1483839);
   grae->SetPointError(4,0,0,0.1134731,0.1134731);
   grae->SetPoint(5,625,-0.02970365);
   grae->SetPointError(5,0,0,0.155483,0.155483);
   grae->SetPoint(6,725,0.1204834);
   grae->SetPointError(6,0,0,0.2201471,0.2201471);
   grae->SetPoint(7,825,0.2635532);
   grae->SetPointError(7,0,0,0.3423838,0.3423838);
   grae->SetPoint(8,925,-0.1704005);
   grae->SetPointError(8,0,0,0.2759024,0.2759024);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.01187622);
   grae->SetPointError(1,0,0,0.08194207,0.08194207);
   grae->SetPoint(2,350,0.05930813);
   grae->SetPointError(2,0,0,0.1105245,0.1105245);
   grae->SetPoint(3,425,0.1142039);
   grae->SetPointError(3,0,0,0.1064149,0.1064149);
   grae->SetPoint(4,525,0.1257219);
   grae->SetPointError(4,0,0,0.130468,0.130468);
   grae->SetPoint(5,625,0.04911519);
   grae->SetPointError(5,0,0,0.1832748,0.1832748);
   grae->SetPoint(6,725,-0.2329434);
   grae->SetPointError(6,0,0,0.2603626,0.2603626);
   grae->SetPoint(7,825,-0.1024017);
   grae->SetPointError(7,0,0,0.3711564,0.3711564);
   grae->SetPoint(8,925,-0.2277692);
   grae->SetPointError(8,0,0,0.3911081,0.3911081);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.1293121);
   grae->SetPointError(1,0,0,0.1651975,0.1651975);
   grae->SetPoint(2,350,-0.5071591);
   grae->SetPointError(2,0,0,0.256843,0.256843);
   grae->SetPoint(3,425,-0.03178712);
   grae->SetPointError(3,0,0,0.2211298,0.2211298);
   grae->SetPoint(4,525,0.3579974);
   grae->SetPointError(4,0,0,0.2759461,0.2759461);
   grae->SetPoint(5,625,0.5976744);
   grae->SetPointError(5,0,0,0.4182476,0.4182476);
   grae->SetPoint(6,725,-0.4657772);
   grae->SetPointError(6,0,0,0.4370661,0.4370661);
   grae->SetPoint(7,825,-0.0481418);
   grae->SetPointError(7,0,0,0.6125295,0.6125295);
   grae->SetPoint(8,925,0.4185124);
   grae->SetPointError(8,0,0,0.6138237,0.6138237);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,-0.0186351);
   grae->SetPointError(3,0,0,0.2436365,0.2436365);
   grae->SetPoint(4,525,0.1247645);
   grae->SetPointError(4,0,0,0.2762266,0.2762266);
   grae->SetPoint(5,625,0.8255158);
   grae->SetPointError(5,0,0,0.5234025,0.5234025);
   grae->SetPoint(6,725,-0.314747);
   grae->SetPointError(6,0,0,0.488046,0.488046);
   grae->SetPoint(7,825,-0.3366879);
   grae->SetPointError(7,0,0,0.680088,0.680088);
   grae->SetPoint(8,925,1.553071);
   grae->SetPointError(8,0,0,1.480891,1.480891);

    for ( uint i = 0; i < nbins; ++i ) {
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

    // 5/fb, 553p2 inclusive

    energy = "8";
    lumi = "5.0";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.15);
    syst.push_back(0.30);
    syst.push_back(0.90);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.04551591);
   grae->SetPointError(1,0,0,0.04388825,0.04388825);
   grae->SetPoint(2,350,-0.06173592);
   grae->SetPointError(2,0,0,0.05665237,0.05665237);
   grae->SetPoint(3,425,-0.001766697);
   grae->SetPointError(3,0,0,0.06312409,0.06312409);
   grae->SetPoint(4,525,-0.0952039);
   grae->SetPointError(4,0,0,0.104175,0.104175);
   grae->SetPoint(5,625,-0.07428241);
   grae->SetPointError(5,0,0,0.170718,0.170718);
   grae->SetPoint(6,725,0.05191541);
   grae->SetPointError(6,0,0,0.291749,0.291749);
   grae->SetPoint(7,825,-0.3560347);
   grae->SetPointError(7,0,0,0.4496893,0.4496893);
   grae->SetPoint(8,925,-1);
   grae->SetPointError(8,0,0,1.079336,1.079336);
   grae->SetPoint(9,1025,-1);
   grae->SetPointError(9,0,0,1.253105,1.253105);
   grae->SetPoint(10,1125,-0.1264823);
   grae->SetPointError(10,0,0,0.977621,0.977621);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.2978023);
   grae->SetPointError(1,0,0,0.03595918,0.03595918);
   grae->SetPoint(2,350,0.1949233);
   grae->SetPointError(2,0,0,0.04596902,0.04596902);
   grae->SetPoint(3,425,0.1968553);
   grae->SetPointError(3,0,0,0.04719326,0.04719326);
   grae->SetPoint(4,525,0.2197711);
   grae->SetPointError(4,0,0,0.06930008,0.06930008);
   grae->SetPoint(5,625,0.1505478);
   grae->SetPointError(5,0,0,0.1036013,0.1036013);
   grae->SetPoint(6,725,0.341073);
   grae->SetPointError(6,0,0,0.170557,0.170557);
   grae->SetPoint(7,825,0.4941215);
   grae->SetPointError(7,0,0,0.2405517,0.2405517);
   grae->SetPoint(8,925,-0.2477387);
   grae->SetPointError(8,0,0,0.328646,0.328646);
   grae->SetPoint(9,1025,-0.1989215);
   grae->SetPointError(9,0,0,0.3816686,0.3816686);
   grae->SetPoint(10,1125,-0.2423733);
   grae->SetPointError(10,0,0,0.3459498,0.3459498);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.03530544);
   grae->SetPointError(1,0,0,0.04579474,0.04579474);
   grae->SetPoint(2,350,0.02770325);
   grae->SetPointError(2,0,0,0.06467179,0.06467179);
   grae->SetPoint(3,425,-0.06656379);
   grae->SetPointError(3,0,0,0.06658049,0.06658049);
   grae->SetPoint(4,525,0.05679914);
   grae->SetPointError(4,0,0,0.09474903,0.09474903);
   grae->SetPoint(5,625,0.07350577);
   grae->SetPointError(5,0,0,0.1471161,0.1471161);
   grae->SetPoint(6,725,0.0840413);
   grae->SetPointError(6,0,0,0.2207806,0.2207806);
   grae->SetPoint(7,825,-0.3653033);
   grae->SetPointError(7,0,0,0.3021225,0.3021225);
   grae->SetPoint(8,925,0.6984769);
   grae->SetPointError(8,0,0,0.724342,0.724342);
   grae->SetPoint(9,1025,-0.6610242);
   grae->SetPointError(9,0,0,0.8315556,0.8315556);
   grae->SetPoint(10,1125,-0.3788335);
   grae->SetPointError(10,0,0,0.7068793,0.7068793);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.05894689);
   grae->SetPointError(1,0,0,0.04878277,0.04878277);
   grae->SetPoint(2,350,0.04404579);
   grae->SetPointError(2,0,0,0.06473565,0.06473565);
   grae->SetPoint(3,425,0.04736159);
   grae->SetPointError(3,0,0,0.06727939,0.06727939);
   grae->SetPoint(4,525,-0.09854121);
   grae->SetPointError(4,0,0,0.1006002,0.1006002);
   grae->SetPoint(5,625,0.2741295);
   grae->SetPointError(5,0,0,0.1508801,0.1508801);
   grae->SetPoint(6,725,0.02511093);
   grae->SetPointError(6,0,0,0.2217853,0.2217853);
   grae->SetPoint(7,825,-0.1470127);
   grae->SetPointError(7,0,0,0.3034018,0.3034018);
   grae->SetPoint(8,925,1.709319);
   grae->SetPointError(8,0,0,1.129369,1.129369);
   grae->SetPoint(9,1025,-0.2064802);
   grae->SetPointError(9,0,0,0.6072,0.6072);
   grae->SetPoint(10,1125,-0.4214689);
   grae->SetPointError(10,0,0,0.4972042,0.4972042);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,-0.08415214);
   grae->SetPointError(3,0,0,0.06862587,0.06862587);
   grae->SetPoint(4,525,-0.2158789);
   grae->SetPointError(4,0,0,0.106198,0.106198);
   grae->SetPoint(5,625,0.1025569);
   grae->SetPointError(5,0,0,0.1564233,0.1564233);
   grae->SetPoint(6,725,-0.1441609);
   grae->SetPointError(6,0,0,0.239534,0.239534);
   grae->SetPoint(7,825,-0.2279739);
   grae->SetPointError(7,0,0,0.350117,0.350117);
   grae->SetPoint(8,925,1.191221);
   grae->SetPointError(8,0,0,0.9295743,0.9295743);
   grae->SetPoint(9,1025,-0.3088652);
   grae->SetPointError(9,0,0,0.8115539,0.8115539);
   grae->SetPoint(10,1125,1.398296);
   grae->SetPointError(10,0,0,2.247238,2.247238);

    for ( uint i = 0; i < nbins; ++i ) {
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

    // 5/fb, 553p2 inclusive, LO WJets

    energy = "8";
    lumi = "5.0";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.15);
    syst.push_back(0.30);
    syst.push_back(0.90);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.0539077);
   grae->SetPointError(1,0,0,0.04346965,0.04346965);
   grae->SetPoint(2,350,-0.06861896);
   grae->SetPointError(2,0,0,0.05644557,0.05644557);
   grae->SetPoint(3,425,-0.005978784);
   grae->SetPointError(3,0,0,0.06296756,0.06296756);
   grae->SetPoint(4,525,-0.09726743);
   grae->SetPointError(4,0,0,0.1040286,0.1040286);
   grae->SetPoint(5,625,-0.07380554);
   grae->SetPointError(5,0,0,0.1706697,0.1706697);
   grae->SetPoint(6,725,0.0571735);
   grae->SetPointError(6,0,0,0.2923091,0.2923091);
   grae->SetPoint(7,825,-0.3496427);
   grae->SetPointError(7,0,0,0.4505994,0.4505994);
   grae->SetPoint(8,925,-1);
   grae->SetPointError(8,0,0,1.075849,1.075849);
   grae->SetPoint(9,1025,-1);
   grae->SetPointError(9,0,0,1.255776,1.255776);
   grae->SetPoint(10,1125,-0.1284646);
   grae->SetPointError(10,0,0,0.9762018,0.9762018);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1640901);
   grae->SetPointError(1,0,0,0.03382819,0.03382819);
   grae->SetPoint(2,350,0.07375513);
   grae->SetPointError(2,0,0,0.04373548,0.04373548);
   grae->SetPoint(3,425,0.07619107);
   grae->SetPointError(3,0,0,0.04490463,0.04490463);
   grae->SetPoint(4,525,0.09553914);
   grae->SetPointError(4,0,0,0.065695,0.065695);
   grae->SetPoint(5,625,0.03708835);
   grae->SetPointError(5,0,0,0.09916304,0.09916304);
   grae->SetPoint(6,725,0.209668);
   grae->SetPointError(6,0,0,0.1594543,0.1594543);
   grae->SetPoint(7,825,0.348466);
   grae->SetPointError(7,0,0,0.2209454,0.2209454);
   grae->SetPoint(8,925,-0.3221379);
   grae->SetPointError(8,0,0,0.323393,0.323393);
   grae->SetPoint(9,1025,-0.2642547);
   grae->SetPointError(9,0,0,0.3774117,0.3774117);
   grae->SetPoint(10,1125,-0.3119179);
   grae->SetPointError(10,0,0,0.3412287,0.3412287);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.06747191);
   grae->SetPointError(1,0,0,0.04541952,0.04541952);
   grae->SetPoint(2,350,-0.008940374);
   grae->SetPointError(2,0,0,0.06393545,0.06393545);
   grae->SetPoint(3,425,-0.1004893);
   grae->SetPointError(3,0,0,0.06602354,0.06602354);
   grae->SetPoint(4,525,0.01971263);
   grae->SetPointError(4,0,0,0.09363279,0.09363279);
   grae->SetPoint(5,625,0.0333879);
   grae->SetPointError(5,0,0,0.1451688,0.1451688);
   grae->SetPoint(6,725,0.04158869);
   grae->SetPointError(6,0,0,0.2176237,0.2176237);
   grae->SetPoint(7,825,-0.3897394);
   grae->SetPointError(7,0,0,0.301169,0.301169);
   grae->SetPoint(8,925,0.6317201);
   grae->SetPointError(8,0,0,0.6932742,0.6932742);
   grae->SetPoint(9,1025,-0.6798079);
   grae->SetPointError(9,0,0,0.823416,0.823416);
   grae->SetPoint(10,1125,-0.4112713);
   grae->SetPointError(10,0,0,0.7020888,0.7020888);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.0525007);
   grae->SetPointError(1,0,0,0.04721084,0.04721084);
   grae->SetPoint(2,350,-0.06674759);
   grae->SetPointError(2,0,0,0.06264574,0.06264574);
   grae->SetPoint(3,425,-0.06435969);
   grae->SetPointError(3,0,0,0.06503731,0.06503731);
   grae->SetPoint(4,525,-0.1923244);
   grae->SetPointError(4,0,0,0.09846368,0.09846368);
   grae->SetPoint(5,625,0.1431219);
   grae->SetPointError(5,0,0,0.1418041,0.1418041);
   grae->SetPoint(6,725,-0.08293527);
   grae->SetPointError(6,0,0,0.2144613,0.2144613);
   grae->SetPoint(7,825,-0.2386887);
   grae->SetPointError(7,0,0,0.2972145,0.2972145);
   grae->SetPoint(8,925,1.408328);
   grae->SetPointError(8,0,0,0.9427036,0.9427036);
   grae->SetPoint(9,1025,-0.2962134);
   grae->SetPointError(9,0,0,0.5934361,0.5934361);
   grae->SetPoint(10,1125,-0.4910262);
   grae->SetPointError(10,0,0,0.4888562,0.4888562);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,-0.08415214);
   grae->SetPointError(3,0,0,0.06862587,0.06862587);
   grae->SetPoint(4,525,-0.2158432);
   grae->SetPointError(4,0,0,0.1061995,0.1061995);
   grae->SetPoint(5,625,0.1025569);
   grae->SetPointError(5,0,0,0.1564233,0.1564233);
   grae->SetPoint(6,725,-0.1441609);
   grae->SetPointError(6,0,0,0.239534,0.239534);
   grae->SetPoint(7,825,-0.2279739);
   grae->SetPointError(7,0,0,0.350117,0.350117);
   grae->SetPoint(8,925,1.191221);
   grae->SetPointError(8,0,0,0.9295743,0.9295743);
   grae->SetPoint(9,1025,-0.3088652);
   grae->SetPointError(9,0,0,0.8115539,0.8115539);
   grae->SetPoint(10,1125,1.398296);
   grae->SetPointError(10,0,0,2.247238,2.247238);

    for ( uint i = 0; i < nbins; ++i ) {
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

    // 5/fb, 553p2 inclusive, LO WJets/DY/Zinv

    energy = "8";
    lumi = "5.0";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.15);
    syst.push_back(0.30);
    syst.push_back(0.90);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.04369646);
   grae->SetPointError(1,0,0,0.03454606,0.03454606);
   grae->SetPoint(2,350,-0.08614341);
   grae->SetPointError(2,0,0,0.0462547,0.0462547);
   grae->SetPoint(3,425,-0.04620612);
   grae->SetPointError(3,0,0,0.05196707,0.05196707);
   grae->SetPoint(4,525,-0.00544409);
   grae->SetPointError(4,0,0,0.08371275,0.08371275);
   grae->SetPoint(5,625,-0.2060631);
   grae->SetPointError(5,0,0,0.1419464,0.1419464);
   grae->SetPoint(6,725,0.04983341);
   grae->SetPointError(6,0,0,0.2387854,0.2387854);
   grae->SetPoint(7,825,0.278086);
   grae->SetPointError(7,0,0,0.3740809,0.3740809);
   grae->SetPoint(8,925,-1);
   grae->SetPointError(8,0,0,0.940454,0.940454);
   grae->SetPoint(9,1025,-1);
   grae->SetPointError(9,0,0,1.156233,1.156233);
   grae->SetPoint(10,1125,-0.3168086);
   grae->SetPointError(10,0,0,0.9090339,0.9090339);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1588869);
   grae->SetPointError(1,0,0,0.03374388,0.03374388);
   grae->SetPoint(2,350,0.06941408);
   grae->SetPointError(2,0,0,0.04365203,0.04365203);
   grae->SetPoint(3,425,0.07169922);
   grae->SetPointError(3,0,0,0.04481655,0.04481655);
   grae->SetPoint(4,525,0.09121781);
   grae->SetPointError(4,0,0,0.06557311,0.06557311);
   grae->SetPoint(5,625,0.03271968);
   grae->SetPointError(5,0,0,0.09898286,0.09898286);
   grae->SetPoint(6,725,0.2053112);
   grae->SetPointError(6,0,0,0.1590989,0.1590989);
   grae->SetPoint(7,825,0.3454099);
   grae->SetPointError(7,0,0,0.2203969,0.2203969);
   grae->SetPoint(8,925,-0.3235124);
   grae->SetPointError(8,0,0,0.3232874,0.3232874);
   grae->SetPoint(9,1025,-0.2662276);
   grae->SetPointError(9,0,0,0.3772493,0.3772493);
   grae->SetPoint(10,1125,-0.3131689);
   grae->SetPointError(10,0,0,0.3411448,0.3411448);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.06904826);
   grae->SetPointError(1,0,0,0.04539424,0.04539424);
   grae->SetPoint(2,350,-0.01029602);
   grae->SetPointError(2,0,0,0.06390437,0.06390437);
   grae->SetPoint(3,425,-0.1018287);
   grae->SetPointError(3,0,0,0.06600035,0.06600035);
   grae->SetPoint(4,525,0.0184776);
   grae->SetPointError(4,0,0,0.09359805,0.09359805);
   grae->SetPoint(5,625,0.03199971);
   grae->SetPointError(5,0,0,0.1450944,0.1450944);
   grae->SetPoint(6,725,0.04025967);
   grae->SetPointError(6,0,0,0.2175314,0.2175314);
   grae->SetPoint(7,825,-0.3919949);
   grae->SetPointError(7,0,0,0.3009637,0.3009637);
   grae->SetPoint(8,925,0.6304989);
   grae->SetPointError(8,0,0,0.6927251,0.6927251);
   grae->SetPoint(9,1025,-0.6804437);
   grae->SetPointError(9,0,0,0.8231403,0.8231403);
   grae->SetPoint(10,1125,-0.4118773);
   grae->SetPointError(10,0,0,0.7020024,0.7020024);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1047921);
   grae->SetPointError(1,0,0,0.04939415,0.04939415);
   grae->SetPoint(2,350,0.09096862);
   grae->SetPointError(2,0,0,0.06564988,0.06564988);
   grae->SetPoint(3,425,0.09359867);
   grae->SetPointError(3,0,0,0.06824156,0.06824156);
   grae->SetPoint(4,525,-0.0556755);
   grae->SetPointError(4,0,0,0.1015223,0.1015223);
   grae->SetPoint(5,625,0.3364184);
   grae->SetPointError(5,0,0,0.1556905,0.1556905);
   grae->SetPoint(6,725,0.07384903);
   grae->SetPointError(6,0,0,0.2254088,0.2254088);
   grae->SetPoint(7,825,-0.1090939);
   grae->SetPointError(7,0,0,0.305895,0.305895);
   grae->SetPoint(8,925,1.802086);
   grae->SetPointError(8,0,0,1.187526,1.187526);
   grae->SetPoint(9,1025,-0.1839468);
   grae->SetPointError(9,0,0,0.609126,0.609126);
   grae->SetPoint(10,1125,-0.4025936);
   grae->SetPointError(10,0,0,0.4985755,0.4985755);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,0.07513038);
   grae->SetPointError(3,0,0,0.07150761,0.07150761);
   grae->SetPoint(4,525,-0.07993755);
   grae->SetPointError(4,0,0,0.1085003,0.1085003);
   grae->SetPoint(5,625,0.2949428);
   grae->SetPointError(5,0,0,0.1690054,0.1690054);
   grae->SetPoint(6,725,0.005900243);
   grae->SetPointError(6,0,0,0.2466829,0.2466829);
   grae->SetPoint(7,825,-0.09193787);
   grae->SetPointError(7,0,0,0.3563153,0.3563153);
   grae->SetPoint(8,925,1.555202);
   grae->SetPointError(8,0,0,1.152673,1.152673);
   grae->SetPoint(9,1025,-0.1954904);
   grae->SetPointError(9,0,0,0.8130912,0.8130912);
   grae->SetPoint(10,1125,1.821677);
   grae->SetPointError(10,0,0,2.798859,2.798859);

    for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 10 ) {

    // 11/fb, 553p2 inclusive, nominal NNLO XS

    energy = "8";
    lumi = "11.6";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.15);
    syst.push_back(0.30);
    syst.push_back(0.90);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.04472765);
   grae->SetPointError(1,0,0,0.02797026,0.02797026);
   grae->SetPoint(2,350,-0.0693895);
   grae->SetPointError(2,0,0,0.03351583,0.03351583);
   grae->SetPoint(3,425,-0.03873836);
   grae->SetPointError(3,0,0,0.03744866,0.03744866);
   grae->SetPoint(4,525,-0.008860317);
   grae->SetPointError(4,0,0,0.05923141,0.05923141);
   grae->SetPoint(5,625,-0.1551768);
   grae->SetPointError(5,0,0,0.09540817,0.09540817);
   grae->SetPoint(6,725,0.1525483);
   grae->SetPointError(6,0,0,0.162338,0.162338);
   grae->SetPoint(7,825,0.0320218);
   grae->SetPointError(7,0,0,0.2604509,0.2604509);
   grae->SetPoint(8,925,-0.4707327);
   grae->SetPointError(8,0,0,0.5000436,0.5000436);
   grae->SetPoint(9,1025,-0.7285871);
   grae->SetPointError(9,0,0,0.7098185,0.7098185);
   grae->SetPoint(10,1125,-0.6702097);
   grae->SetPointError(10,0,0,0.7494358,0.7494358);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.2716923);
   grae->SetPointError(1,0,0,0.02679599,0.02679599);
   grae->SetPoint(2,350,0.233476);
   grae->SetPointError(2,0,0,0.03334648,0.03334648);
   grae->SetPoint(3,425,0.1850367);
   grae->SetPointError(3,0,0,0.0341338,0.0341338);
   grae->SetPoint(4,525,0.2357377);
   grae->SetPointError(4,0,0,0.04964064,0.04964064);
   grae->SetPoint(5,625,0.0908104);
   grae->SetPointError(5,0,0,0.07011198,0.07011198);
   grae->SetPoint(6,725,0.1887614);
   grae->SetPointError(6,0,0,0.109259,0.109259);
   grae->SetPoint(7,825,0.3545164);
   grae->SetPointError(7,0,0,0.167739,0.167739);
   grae->SetPoint(8,925,-0.03675048);
   grae->SetPointError(8,0,0,0.220042,0.220042);
   grae->SetPoint(9,1025,-0.2067498);
   grae->SetPointError(9,0,0,0.271977,0.271977);
   grae->SetPoint(10,1125,0.239281);
   grae->SetPointError(10,0,0,0.2597612,0.2597612);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,-0.02697954);
   grae->SetPointError(1,0,0,0.03290404,0.03290404);
   grae->SetPoint(2,350,0.001190671);
   grae->SetPointError(2,0,0,0.04456932,0.04456932);
   grae->SetPoint(3,425,0.03265944);
   grae->SetPointError(3,0,0,0.04704095,0.04704095);
   grae->SetPoint(4,525,0.01139806);
   grae->SetPointError(4,0,0,0.06609412,0.06609412);
   grae->SetPoint(5,625,0.08528368);
   grae->SetPointError(5,0,0,0.1011062,0.1011062);
   grae->SetPoint(6,725,0.05384737);
   grae->SetPointError(6,0,0,0.1530462,0.1530462);
   grae->SetPoint(7,825,-0.1873718);
   grae->SetPointError(7,0,0,0.2197093,0.2197093);
   grae->SetPoint(8,925,0.08770927);
   grae->SetPointError(8,0,0,0.3593508,0.3593508);
   grae->SetPoint(9,1025,-0.361775);
   grae->SetPointError(9,0,0,0.5162023,0.5162023);
   grae->SetPoint(10,1125,-0.6005455);
   grae->SetPointError(10,0,0,0.4530595,0.4530595);
   
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.06229846);
   grae->SetPointError(1,0,0,0.03906769,0.03906769);
   grae->SetPoint(2,350,0.02375957);
   grae->SetPointError(2,0,0,0.05083061,0.05083061);
   grae->SetPoint(3,425,0.05077076);
   grae->SetPointError(3,0,0,0.05294899,0.05294899);
   grae->SetPoint(4,525,-0.02639442);
   grae->SetPointError(4,0,0,0.0785245,0.0785245);
   grae->SetPoint(5,625,0.113949);
   grae->SetPointError(5,0,0,0.1105075,0.1105075);
   grae->SetPoint(6,725,-0.0887129);
   grae->SetPointError(6,0,0,0.1683529,0.1683529);
   grae->SetPoint(7,825,-0.01935264);
   grae->SetPointError(7,0,0,0.2383442,0.2383442);
   grae->SetPoint(8,925,1.785808);
   grae->SetPointError(8,0,0,0.8626186,0.8626186);
   grae->SetPoint(9,1025,-0.02005198);
   grae->SetPointError(9,0,0,0.446597,0.446597);
   grae->SetPoint(10,1125,-0.2352875);
   grae->SetPointError(10,0,0,0.3653169,0.3653169);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,-0.1029378);
   grae->SetPointError(3,0,0,0.05757584,0.05757584);
   grae->SetPoint(4,525,-0.1849123);
   grae->SetPointError(4,0,0,0.0887351,0.0887351);
   grae->SetPoint(5,625,0.0006743843);
   grae->SetPointError(5,0,0,0.1313864,0.1313864);
   grae->SetPoint(6,725,-0.2121653);
   grae->SetPointError(6,0,0,0.2062823,0.2062823);
   grae->SetPoint(7,825,-0.2253278);
   grae->SetPointError(7,0,0,0.2998814,0.2998814);
   grae->SetPoint(8,925,1.273781);
   grae->SetPointError(8,0,0,0.8023161,0.8023161);
   grae->SetPoint(9,1025,-0.1910069);
   grae->SetPointError(9,0,0,0.6994237,0.6994237);
   grae->SetPoint(10,1125,1.896447);
   grae->SetPointError(10,0,0,2.511828,2.511828);

    for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 110 ) {

    // 11/fb, 553p2 inclusive, darren's k-factors for WJets/DY/Zinv

    energy = "8";
    lumi = "11.6";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.15);
    syst.push_back(0.30);
    syst.push_back(0.90);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
       grae->SetPoint(1,300,-0.05326944);
   grae->SetPointError(1,0,0,0.02721819,0.02721819);
   grae->SetPoint(2,350,-0.07678704);
   grae->SetPointError(2,0,0,0.03326405,0.03326405);
   grae->SetPoint(3,425,-0.04376405);
   grae->SetPointError(3,0,0,0.03726134,0.03726134);
   grae->SetPoint(4,525,-0.009469896);
   grae->SetPointError(4,0,0,0.05908797,0.05908797);
   grae->SetPoint(5,625,-0.1523521);
   grae->SetPointError(5,0,0,0.09535056,0.09535056);
   grae->SetPoint(6,725,0.1696597);
   grae->SetPointError(6,0,0,0.163568,0.163568);
   grae->SetPoint(7,825,0.05731179);
   grae->SetPointError(7,0,0,0.2635855,0.2635855);
   grae->SetPoint(8,925,-0.5975212);
   grae->SetPointError(8,0,0,0.358875,0.358875);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.06407818);
   grae->SetPointError(1,0,0,0.02476339,0.02476339);
   grae->SetPoint(2,350,0.03396877);
   grae->SetPointError(2,0,0,0.03086003,0.03086003);
   grae->SetPoint(3,425,-0.002208437);
   grae->SetPointError(3,0,0,0.03190943,0.03190943);
   grae->SetPoint(4,525,0.04049819);
   grae->SetPointError(4,0,0,0.04577476,0.04577476);
   grae->SetPoint(5,625,-0.0660106);
   grae->SetPointError(5,0,0,0.06658213,0.06658213);
   grae->SetPoint(6,725,0.02110587);
   grae->SetPointError(6,0,0,0.1021119,0.1021119);
   grae->SetPoint(7,825,0.1324085);
   grae->SetPointError(7,0,0,0.1501536,0.1501536);
   grae->SetPoint(8,925,-0.1482618);
   grae->SetPointError(8,0,0,0.1346768,0.1346768);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.08320381);
   grae->SetPointError(1,0,0,0.03232637,0.03232637);
   grae->SetPoint(2,350,-0.05148706);
   grae->SetPointError(2,0,0,0.04377343,0.04377343);
   grae->SetPoint(3,425,-0.02413611);
   grae->SetPointError(3,0,0,0.04611211,0.04611211);
   grae->SetPoint(4,525,-0.05069055);
   grae->SetPointError(4,0,0,0.06467226,0.06467226);
   grae->SetPoint(5,625,0.007152619);
   grae->SetPointError(5,0,0,0.09814455,0.09814455);
   grae->SetPoint(6,725,-0.02827604);
   grae->SetPointError(6,0,0,0.1485783,0.1485783);
   grae->SetPoint(7,825,-0.2238742);
   grae->SetPointError(7,0,0,0.2172615,0.2172615);
   grae->SetPoint(8,925,-0.322685);
   grae->SetPointError(8,0,0,0.2346623,0.2346623);
   
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1218363);
   grae->SetPointError(1,0,0,0.03924644,0.03924644);
   grae->SetPoint(2,350,0.08337729);
   grae->SetPointError(2,0,0,0.05112376,0.05112376);
   grae->SetPoint(3,425,0.1113075);
   grae->SetPointError(3,0,0,0.05340018,0.05340018);
   grae->SetPoint(4,525,0.03340948);
   grae->SetPointError(4,0,0,0.0789309,0.0789309);
   grae->SetPoint(5,625,0.1844923);
   grae->SetPointError(5,0,0,0.1124332,0.1124332);
   grae->SetPoint(6,725,-0.03216693);
   grae->SetPointError(6,0,0,0.1689611,0.1689611);
   grae->SetPoint(7,825,0.03774525);
   grae->SetPointError(7,0,0,0.2406252,0.2406252);
   grae->SetPoint(8,925,0.4975177);
   grae->SetPointError(8,0,0,0.2616271,0.2616271);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,-0.0503512);
   grae->SetPointError(3,0,0,0.05432958,0.05432958);
   grae->SetPoint(4,525,-0.1424387);
   grae->SetPointError(4,0,0,0.08251005,0.08251005);
   grae->SetPoint(5,625,0.1313328);
   grae->SetPointError(5,0,0,0.1241772,0.1241772);
   grae->SetPoint(6,725,-0.1408899);
   grae->SetPointError(6,0,0,0.1885196,0.1885196);
   grae->SetPoint(7,825,-0.2416665);
   grae->SetPointError(7,0,0,0.2663151,0.2663151);
   grae->SetPoint(8,925,0.2821308);
   grae->SetPointError(8,0,0,0.2932455,0.2932455);

    for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 111 ) {

    // 11/fb, 553p2 inclusive, NNLO

    energy = "8";
    lumi = "11.6";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.15);
    syst.push_back(0.30);
    syst.push_back(0.90);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.03102826);
   grae->SetPointError(1,0,0,0.0308242,0.0308242);
   grae->SetPoint(2,350,-0.05226311);
   grae->SetPointError(2,0,0,0.03376251,0.03376251);
   grae->SetPoint(3,425,-0.04372183);
   grae->SetPointError(3,0,0,0.03675145,0.03675145);
   grae->SetPoint(4,525,-0.01132264);
   grae->SetPointError(4,0,0,0.05825812,0.05825812);
   grae->SetPoint(5,625,-0.1558879);
   grae->SetPointError(5,0,0,0.09312599,0.09312599);
   grae->SetPoint(6,725,0.1465856);
   grae->SetPointError(6,0,0,0.1576263,0.1576263);
   grae->SetPoint(7,825,0.005220447);
   grae->SetPointError(7,0,0,0.2538092,0.2538092);
   grae->SetPoint(8,925,-0.6027043);
   grae->SetPointError(8,0,0,0.3488301,0.3488301);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.288736);
   grae->SetPointError(1,0,0,0.02811059,0.02811059);
   grae->SetPoint(2,350,0.196813);
   grae->SetPointError(2,0,0,0.03213419,0.03213419);
   grae->SetPoint(3,425,0.1434895);
   grae->SetPointError(3,0,0,0.03249779,0.03249779);
   grae->SetPoint(4,525,0.1655189);
   grae->SetPointError(4,0,0,0.0464262,0.0464262);
   grae->SetPoint(5,625,0.06945767);
   grae->SetPointError(5,0,0,0.06579184,0.06579184);
   grae->SetPoint(6,725,0.1417698);
   grae->SetPointError(6,0,0,0.1009796,0.1009796);
   grae->SetPoint(7,825,0.3075983);
   grae->SetPointError(7,0,0,0.1543874,0.1543874);
   grae->SetPoint(8,925,0.008681426);
   grae->SetPointError(8,0,0,0.1339139,0.1339139);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.003636446);
   grae->SetPointError(1,0,0,0.03225708,0.03225708);
   grae->SetPoint(2,350,-0.007688387);
   grae->SetPointError(2,0,0,0.0449688,0.0449688);
   grae->SetPoint(3,425,0.04602087);
   grae->SetPointError(3,0,0,0.04590317,0.04590317);
   grae->SetPoint(4,525,0.04548716);
   grae->SetPointError(4,0,0,0.06425633,0.06425633);
   grae->SetPoint(5,625,0.09814623);
   grae->SetPointError(5,0,0,0.09697728,0.09697728);
   grae->SetPoint(6,725,0.05175404);
   grae->SetPointError(6,0,0,0.1465807,0.1465807);
   grae->SetPoint(7,825,-0.1511205);
   grae->SetPointError(7,0,0,0.2092724,0.2092724);
   grae->SetPoint(8,925,-0.2720603);
   grae->SetPointError(8,0,0,0.2297656,0.2297656);
   
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1120709);
   grae->SetPointError(1,0,0,0.03154878,0.03154878);
   grae->SetPoint(2,350,0.06042295);
   grae->SetPointError(2,0,0,0.03983448,0.03983448);
   grae->SetPoint(3,425,0.04495716);
   grae->SetPointError(3,0,0,0.0402217,0.0402217);
   grae->SetPoint(4,525,-0.06420767);
   grae->SetPointError(4,0,0,0.05805213,0.05805213);
   grae->SetPoint(5,625,0.1276312);
   grae->SetPointError(5,0,0,0.08490047,0.08490047);
   grae->SetPoint(6,725,-0.1289624);
   grae->SetPointError(6,0,0,0.124976,0.124976);
   grae->SetPoint(7,825,-0.0963107);
   grae->SetPointError(7,0,0,0.1776318,0.1776318);
   grae->SetPoint(8,925,0.006458404);
   grae->SetPointError(8,0,0,0.1603616,0.1603616);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,-0.04373164);
   grae->SetPointError(3,0,0,0.0434185,0.0434185);
   grae->SetPoint(4,525,-0.1692159);
   grae->SetPointError(4,0,0,0.06530448,0.06530448);
   grae->SetPoint(5,625,0.1639127);
   grae->SetPointError(5,0,0,0.1056159,0.1056159);
   grae->SetPoint(6,725,-0.1381286);
   grae->SetPointError(6,0,0,0.1563199,0.1563199);
   grae->SetPoint(7,825,-0.1375594);
   grae->SetPointError(7,0,0,0.2254803,0.2254803);
   grae->SetPoint(8,925,0.0174795);
   grae->SetPointError(8,0,0,0.2332812,0.2332812);

    for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 11 ) {

    // 11/fb, 553p2 inclusive, LO WJets/DY/Zinv

    energy = "8";
    lumi = "11.6";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.15);
    syst.push_back(0.30);
    syst.push_back(0.90);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.05293503);
   grae->SetPointError(1,0,0,0.02743293,0.02743293);
   grae->SetPoint(2,350,-0.07700453);
   grae->SetPointError(2,0,0,0.03335601,0.03335601);
   grae->SetPoint(3,425,-0.04357718);
   grae->SetPointError(3,0,0,0.03735359,0.03735359);
   grae->SetPoint(4,525,-0.01013511);
   grae->SetPointError(4,0,0,0.05923718,0.05923718);
   grae->SetPoint(5,625,-0.1534659);
   grae->SetPointError(5,0,0,0.09559843,0.09559843);
   grae->SetPoint(6,725,0.1656907);
   grae->SetPointError(6,0,0,0.1637356,0.1637356);
   grae->SetPoint(7,825,0.05251158);
   grae->SetPointError(7,0,0,0.2633199,0.2633199);
   grae->SetPoint(8,925,-0.4716888);
   grae->SetPointError(8,0,0,0.5005874,0.5005874);
   grae->SetPoint(9,1025,-0.7253268);
   grae->SetPointError(9,0,0,0.7180253,0.7180253);
   grae->SetPoint(10,1125,-0.6716884);
   grae->SetPointError(10,0,0,0.7500792,0.7500792);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1351431);
   grae->SetPointError(1,0,0,0.02542201,0.02542201);
   grae->SetPoint(2,350,0.1038273);
   grae->SetPointError(2,0,0,0.03171034,0.03171034);
   grae->SetPoint(3,425,0.06134457);
   grae->SetPointError(3,0,0,0.03266981,0.03266981);
   grae->SetPoint(4,525,0.1059089);
   grae->SetPointError(4,0,0,0.04715596,0.04715596);
   grae->SetPoint(5,625,-0.02124414);
   grae->SetPointError(5,0,0,0.06786822,0.06786822);
   grae->SetPoint(6,725,0.06835945);
   grae->SetPointError(6,0,0,0.1045784,0.1045784);
   grae->SetPoint(7,825,0.2201414);
   grae->SetPointError(7,0,0,0.1571294,0.1571294);
   grae->SetPoint(8,925,-0.1338325);
   grae->SetPointError(8,0,0,0.2150535,0.2150535);
   grae->SetPoint(9,1025,-0.2739195);
   grae->SetPointError(9,0,0,0.2707102,0.2707102);
   grae->SetPoint(10,1125,0.1234567);
   grae->SetPointError(10,0,0,0.2478308,0.2478308);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.06106963);
   grae->SetPointError(1,0,0,0.03266276,0.03266276);
   grae->SetPoint(2,350,-0.03584204);
   grae->SetPointError(2,0,0,0.04415882,0.04415882);
   grae->SetPoint(3,425,-0.006394232);
   grae->SetPointError(3,0,0,0.04655098,0.04655098);
   grae->SetPoint(4,525,-0.02518445);
   grae->SetPointError(4,0,0,0.06553118,0.06553118);
   grae->SetPoint(5,625,0.04329015);
   grae->SetPointError(5,0,0,0.09986558,0.09986558);
   grae->SetPoint(6,725,0.01120909);
   grae->SetPointError(6,0,0,0.1513383,0.1513383);
   grae->SetPoint(7,825,-0.2215862);
   grae->SetPointError(7,0,0,0.2187708,0.2187708);
   grae->SetPoint(8,925,0.04411602);
   grae->SetPointError(8,0,0,0.3542408,0.3542408);
   grae->SetPoint(9,1025,-0.3981753);
   grae->SetPointError(9,0,0,0.513788,0.513788);
   grae->SetPoint(10,1125,-0.6218312);
   grae->SetPointError(10,0,0,0.4506003,0.4506003);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1082205);
   grae->SetPointError(1,0,0,0.03932382,0.03932382);
   grae->SetPoint(2,350,0.06971762);
   grae->SetPointError(2,0,0,0.05119974,0.05119974);
   grae->SetPoint(3,425,0.09734531);
   grae->SetPointError(3,0,0,0.05342853,0.05342853);
   grae->SetPoint(4,525,0.01954224);
   grae->SetPointError(4,0,0,0.0790247,0.0790247);
   grae->SetPoint(5,625,0.1688629);
   grae->SetPointError(5,0,0,0.1122823,0.1122823);
   grae->SetPoint(6,725,-0.04527908);
   grae->SetPointError(6,0,0,0.1692812,0.1692812);
   grae->SetPoint(7,825,0.02455099);
   grae->SetPointError(7,0,0,0.2405246,0.2405246);
   grae->SetPoint(8,925,1.880537);
   grae->SetPointError(8,0,0,0.9026659,0.9026659);
   grae->SetPoint(9,1025,0.008009069);
   grae->SetPointError(9,0,0,0.4479901,0.4479901);
   grae->SetPoint(10,1125,-0.2105265);
   grae->SetPointError(10,0,0,0.3655223,0.3655223);    

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,0.05319363);
   grae->SetPointError(3,0,0,0.05878618,0.05878618);
   grae->SetPoint(4,525,-0.04402408);
   grae->SetPointError(4,0,0,0.08943532,0.08943532);
   grae->SetPoint(5,625,0.1755666);
   grae->SetPointError(5,0,0,0.1366582,0.1366582);
   grae->SetPoint(6,725,-0.07386437);
   grae->SetPointError(6,0,0,0.207088,0.207088);
   grae->SetPoint(7,825,-0.08880057);
   grae->SetPointError(7,0,0,0.3003513,0.3003513);
   grae->SetPoint(8,925,1.651309);
   grae->SetPointError(8,0,0,0.9811638,0.9811638);
   grae->SetPoint(9,1025,-0.05799386);
   grae->SetPointError(9,0,0,0.6985753,0.6985753);
   grae->SetPoint(10,1125,2.406475);
   grae->SetPointError(10,0,0,3.103424,3.103424);

    for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 12 ) {

    // 5/fb, 553p2 inclusive, standard, 8 bins 

    energy = "8";
    lumi = "5.1";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.15);
    syst.push_back(0.30);
    syst.push_back(0.90);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.03521584);
   grae->SetPointError(1,0,0,0.03501903,0.03501903);
   grae->SetPoint(2,350,-0.07867912);
   grae->SetPointError(2,0,0,0.04643089,0.04643089);
   grae->SetPoint(3,425,-0.04130978);
   grae->SetPointError(3,0,0,0.0520864,0.0520864);
   grae->SetPoint(4,525,-0.003675133);
   grae->SetPointError(4,0,0,0.0837503,0.0837503);
   grae->SetPoint(5,625,-0.2074532);
   grae->SetPointError(5,0,0,0.1417806,0.1417806);
   grae->SetPoint(6,725,0.03796963);
   grae->SetPointError(6,0,0,0.2372576,0.2372576);
   grae->SetPoint(7,825,0.2535281);
   grae->SetPointError(7,0,0,0.3682382,0.3682382);
   grae->SetPoint(8,925,-0.8245295);
   grae->SetPointError(8,0,0,0.569061,0.569061);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.2978023);
   grae->SetPointError(1,0,0,0.03595918,0.03595918);
   grae->SetPoint(2,350,0.1949233);
   grae->SetPointError(2,0,0,0.04596902,0.04596902);
   grae->SetPoint(3,425,0.1968553);
   grae->SetPointError(3,0,0,0.04719326,0.04719326);
   grae->SetPoint(4,525,0.2197711);
   grae->SetPointError(4,0,0,0.06930008,0.06930008);
   grae->SetPoint(5,625,0.1505478);
   grae->SetPointError(5,0,0,0.1036013,0.1036013);
   grae->SetPoint(6,725,0.341073);
   grae->SetPointError(6,0,0,0.170557,0.170557);
   grae->SetPoint(7,825,0.4941215);
   grae->SetPointError(7,0,0,0.2405517,0.2405517);
   grae->SetPoint(8,925,-0.2351277);
   grae->SetPointError(8,0,0,0.201622,0.201622);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.03530544);
   grae->SetPointError(1,0,0,0.04579474,0.04579474);
   grae->SetPoint(2,350,0.02770325);
   grae->SetPointError(2,0,0,0.06467179,0.06467179);
   grae->SetPoint(3,425,-0.06656379);
   grae->SetPointError(3,0,0,0.06658049,0.06658049);
   grae->SetPoint(4,525,0.05679914);
   grae->SetPointError(4,0,0,0.09474903,0.09474903);
   grae->SetPoint(5,625,0.07350577);
   grae->SetPointError(5,0,0,0.1471161,0.1471161);
   grae->SetPoint(6,725,0.0840413);
   grae->SetPointError(6,0,0,0.2207806,0.2207806);
   grae->SetPoint(7,825,-0.3653033);
   grae->SetPointError(7,0,0,0.3021225,0.3021225);
   grae->SetPoint(8,925,-0.03785047);
   grae->SetPointError(8,0,0,0.3713945,0.3713945);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.05894573);
   grae->SetPointError(1,0,0,0.0487812,0.0487812);
   grae->SetPoint(2,350,0.0442105);
   grae->SetPointError(2,0,0,0.06473845,0.06473845);
   grae->SetPoint(3,425,0.04757992);
   grae->SetPointError(3,0,0,0.06728443,0.06728443);
   grae->SetPoint(4,525,-0.09829882);
   grae->SetPointError(4,0,0,0.1006065,0.1006065);
   grae->SetPoint(5,625,0.2746089);
   grae->SetPointError(5,0,0,0.1509137,0.1509137);
   grae->SetPoint(6,725,0.02531396);
   grae->SetPointError(6,0,0,0.2217891,0.2217891);
   grae->SetPoint(7,825,-0.1470402);
   grae->SetPointError(7,0,0,0.3033662,0.3033662);
   grae->SetPoint(8,925,0.2532282);
   grae->SetPointError(8,0,0,0.3183834,0.3183834);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,-0.08415214);
   grae->SetPointError(3,0,0,0.06862587,0.06862587);
   grae->SetPoint(4,525,-0.2158789);
   grae->SetPointError(4,0,0,0.106198,0.106198);
   grae->SetPoint(5,625,0.1029119);
   grae->SetPointError(5,0,0,0.1564342,0.1564342);
   grae->SetPoint(6,725,-0.1443433);
   grae->SetPointError(6,0,0,0.2395165,0.2395165);
   grae->SetPoint(7,825,-0.2281211);
   grae->SetPointError(7,0,0,0.3500768,0.3500768);
   grae->SetPoint(8,925,0.3185155);
   grae->SetPointError(8,0,0,0.4191179,0.4191179);

    for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 13 ) {

    // 5/fb, 553p2 inclusive, mixture, 8 bins 

    energy = "8";
    lumi = "5.1";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.15);
    syst.push_back(0.30);
    syst.push_back(0.90);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.04370573);
   grae->SetPointError(1,0,0,0.03454504,0.03454504);
   grae->SetPoint(2,350,-0.08593754);
   grae->SetPointError(2,0,0,0.04625758,0.04625758);
   grae->SetPoint(3,425,-0.04615293);
   grae->SetPointError(3,0,0,0.05196745,0.05196745);
   grae->SetPoint(4,525,-0.005285702);
   grae->SetPointError(4,0,0,0.08371743,0.08371743);
   grae->SetPoint(5,625,-0.2058542);
   grae->SetPointError(5,0,0,0.1419541,0.1419541);
   grae->SetPoint(6,725,0.05009956);
   grae->SetPointError(6,0,0,0.2388121,0.2388121);
   grae->SetPoint(7,825,0.2781935);
   grae->SetPointError(7,0,0,0.3740974,0.3740974);
   grae->SetPoint(8,925,-0.8244216);
   grae->SetPointError(8,0,0,0.5701106,0.5701106);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1588869);
   grae->SetPointError(1,0,0,0.03374388,0.03374388);
   grae->SetPoint(2,350,0.06941408);
   grae->SetPointError(2,0,0,0.04365203,0.04365203);
   grae->SetPoint(3,425,0.07169922);
   grae->SetPointError(3,0,0,0.04481655,0.04481655);
   grae->SetPoint(4,525,0.09121781);
   grae->SetPointError(4,0,0,0.06557311,0.06557311);
   grae->SetPoint(5,625,0.03271968);
   grae->SetPointError(5,0,0,0.09898286,0.09898286);
   grae->SetPoint(6,725,0.2053112);
   grae->SetPointError(6,0,0,0.1590989,0.1590989);
   grae->SetPoint(7,825,0.3454099);
   grae->SetPointError(7,0,0,0.2203969,0.2203969);
   grae->SetPoint(8,925,-0.3074322);
   grae->SetPointError(8,0,0,0.1987152,0.1987152);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.06904826);
   grae->SetPointError(1,0,0,0.04539424,0.04539424);
   grae->SetPoint(2,350,-0.01029602);
   grae->SetPointError(2,0,0,0.06390437,0.06390437);
   grae->SetPoint(3,425,-0.1018287);
   grae->SetPointError(3,0,0,0.06600035,0.06600035);
   grae->SetPoint(4,525,0.0184776);
   grae->SetPointError(4,0,0,0.09359805,0.09359805);
   grae->SetPoint(5,625,0.03199971);
   grae->SetPointError(5,0,0,0.1450944,0.1450944);
   grae->SetPoint(6,725,0.04025967);
   grae->SetPointError(6,0,0,0.2175314,0.2175314);
   grae->SetPoint(7,825,-0.3919949);
   grae->SetPointError(7,0,0,0.3009637,0.3009637);
   grae->SetPoint(8,925,-0.08433129);
   grae->SetPointError(8,0,0,0.3664553,0.3664553);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1047848);
   grae->SetPointError(1,0,0,0.04939226,0.04939226);
   grae->SetPoint(2,350,0.09115966);
   grae->SetPointError(2,0,0,0.06565371,0.06565371);
   grae->SetPoint(3,425,0.09385497);
   grae->SetPointError(3,0,0,0.06824826,0.06824826);
   grae->SetPoint(4,525,-0.05539108);
   grae->SetPointError(4,0,0,0.1015307,0.1015307);
   grae->SetPoint(5,625,0.336974);
   grae->SetPointError(5,0,0,0.1557343,0.1557343);
   grae->SetPoint(6,725,0.07434181);
   grae->SetPointError(6,0,0,0.2254519,0.2254519);
   grae->SetPoint(7,825,-0.1086907);
   grae->SetPointError(7,0,0,0.3059287,0.3059287);
   grae->SetPoint(8,925,0.2951331);
   grae->SetPointError(8,0,0,0.3249491,0.3249491);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,0.07513038);
   grae->SetPointError(3,0,0,0.07150761,0.07150761);
   grae->SetPoint(4,525,-0.07993755);
   grae->SetPointError(4,0,0,0.1085003,0.1085003);
   grae->SetPoint(5,625,0.2953368);
   grae->SetPointError(5,0,0,0.1690295,0.1690295);
   grae->SetPoint(6,725,0.005897748);
   grae->SetPointError(6,0,0,0.2466826,0.2466826);
   grae->SetPoint(7,825,-0.09171677);
   grae->SetPointError(7,0,0,0.3563247,0.3563247);
   grae->SetPoint(8,925,0.5429077);
   grae->SetPointError(8,0,0,0.4712034,0.4712034);

    for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 14 ) {

    // 11/fb, 553p2 inclusive, mixture, 8 bins 

    energy = "8";
    lumi = "11.7";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.15);
    syst.push_back(0.30);
    syst.push_back(0.90);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.05315763);
   grae->SetPointError(1,0,0,0.02742961,0.02742961);
   grae->SetPoint(2,350,-0.07687933);
   grae->SetPointError(2,0,0,0.03335667,0.03335667);
   grae->SetPoint(3,425,-0.04356434);
   grae->SetPointError(3,0,0,0.03735331,0.03735331);
   grae->SetPoint(4,525,-0.0101238);
   grae->SetPointError(4,0,0,0.05923594,0.05923594);
   grae->SetPoint(5,625,-0.1533897);
   grae->SetPointError(5,0,0,0.09559697,0.09559697);
   grae->SetPoint(6,725,0.1656802);
   grae->SetPointError(6,0,0,0.1637282,0.1637282);
   grae->SetPoint(7,825,0.05223907);
   grae->SetPointError(7,0,0,0.2632692,0.2632692);
   grae->SetPoint(8,925,-0.5981793);
   grae->SetPointError(8,0,0,0.3594511,0.3594511);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1351431);
   grae->SetPointError(1,0,0,0.02542201,0.02542201);
   grae->SetPoint(2,350,0.1038273);
   grae->SetPointError(2,0,0,0.03171034,0.03171034);
   grae->SetPoint(3,425,0.06134457);
   grae->SetPointError(3,0,0,0.03266981,0.03266981);
   grae->SetPoint(4,525,0.1059089);
   grae->SetPointError(4,0,0,0.04715596,0.04715596);
   grae->SetPoint(5,625,-0.02124414);
   grae->SetPointError(5,0,0,0.06786822,0.06786822);
   grae->SetPoint(6,725,0.06835945);
   grae->SetPointError(6,0,0,0.1045784,0.1045784);
   grae->SetPoint(7,825,0.2201414);
   grae->SetPointError(7,0,0,0.1571294,0.1571294);
   grae->SetPoint(8,925,-0.09482368);
   grae->SetPointError(8,0,0,0.1375655,0.1375655);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.06106963);
   grae->SetPointError(1,0,0,0.03266276,0.03266276);
   grae->SetPoint(2,350,-0.03584204);
   grae->SetPointError(2,0,0,0.04415882,0.04415882);
   grae->SetPoint(3,425,-0.006394232);
   grae->SetPointError(3,0,0,0.04655098,0.04655098);
   grae->SetPoint(4,525,-0.02518445);
   grae->SetPointError(4,0,0,0.06553118,0.06553118);
   grae->SetPoint(5,625,0.04329015);
   grae->SetPointError(5,0,0,0.09986558,0.09986558);
   grae->SetPoint(6,725,0.01120909);
   grae->SetPointError(6,0,0,0.1513383,0.1513383);
   grae->SetPoint(7,825,-0.2215862);
   grae->SetPointError(7,0,0,0.2187708,0.2187708);
   grae->SetPoint(8,925,-0.2982789);
   grae->SetPointError(8,0,0,0.2372008,0.2372008);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1082156);
   grae->SetPointError(1,0,0,0.03932267,0.03932267);
   grae->SetPoint(2,350,0.06990228);
   grae->SetPointError(2,0,0,0.0512016,0.0512016);
   grae->SetPoint(3,425,0.09745933);
   grae->SetPointError(3,0,0,0.05343038,0.05343038);
   grae->SetPoint(4,525,0.01976576);
   grae->SetPointError(4,0,0,0.07902915,0.07902915);
   grae->SetPoint(5,625,0.1693077);
   grae->SetPointError(5,0,0,0.1123009,0.1123009);
   grae->SetPoint(6,725,-0.04509585);
   grae->SetPointError(6,0,0,0.1692792,0.1692792);
   grae->SetPoint(7,825,0.02493317);
   grae->SetPointError(7,0,0,0.240549,0.240549);
   grae->SetPoint(8,925,0.4818484);
   grae->SetPointError(8,0,0,0.2602567,0.2602567);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,0.05319363);
   grae->SetPointError(3,0,0,0.05878618,0.05878618);
   grae->SetPoint(4,525,-0.04402655);
   grae->SetPointError(4,0,0,0.08943522,0.08943522);
   grae->SetPoint(5,625,0.1760582);
   grae->SetPointError(5,0,0,0.1366745,0.1366745);
   grae->SetPoint(6,725,-0.0739753);
   grae->SetPointError(6,0,0,0.2070785,0.2070785);
   grae->SetPoint(7,825,-0.08859914);
   grae->SetPointError(7,0,0,0.300348,0.300348);
   grae->SetPoint(8,925,0.693386);
   grae->SetPointError(8,0,0,0.427462,0.427462);

    for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 15 ) {

    // 11/fb, 553p2, 2-3, mixture, 8 bins 

    energy = "8";
    lumi = "11.7";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.15);
    syst.push_back(0.30);
    syst.push_back(0.90);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.0546802);
   grae->SetPointError(1,0,0,0.03025781,0.03025781);
   grae->SetPoint(2,350,-0.06052873);
   grae->SetPointError(2,0,0,0.03694661,0.03694661);
   grae->SetPoint(3,425,-0.03560808);
   grae->SetPointError(3,0,0,0.04184187,0.04184187);
   grae->SetPoint(4,525,-0.07427342);
   grae->SetPointError(4,0,0,0.07445687,0.07445687);
   grae->SetPoint(5,625,-0.1416568);
   grae->SetPointError(5,0,0,0.1269872,0.1269872);
   grae->SetPoint(6,725,0.02794661);
   grae->SetPointError(6,0,0,0.2256453,0.2256453);
   grae->SetPoint(7,825,-0.3371197);
   grae->SetPointError(7,0,0,0.3771194,0.3771194);
   grae->SetPoint(8,925,-0.5331764);
   grae->SetPointError(8,0,0,0.4782493,0.4782493);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.160398);
   grae->SetPointError(1,0,0,0.02885296,0.02885296);
   grae->SetPoint(2,350,0.1576586);
   grae->SetPointError(2,0,0,0.03643444,0.03643444);
   grae->SetPoint(3,425,0.1234905);
   grae->SetPointError(3,0,0,0.03802028,0.03802028);
   grae->SetPoint(4,525,0.1485403);
   grae->SetPointError(4,0,0,0.06239629,0.06239629);
   grae->SetPoint(5,625,0.07797849);
   grae->SetPointError(5,0,0,0.09611793,0.09611793);
   grae->SetPoint(6,725,0.04869653);
   grae->SetPointError(6,0,0,0.1604718,0.1604718);
   grae->SetPoint(7,825,0.1536973);
   grae->SetPointError(7,0,0,0.2358036,0.2358036);
   grae->SetPoint(8,925,-0.04202397);
   grae->SetPointError(8,0,0,0.2119533,0.2119533);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.04383);
   grae->SetPointError(1,0,0,0.04206461,0.04206461);
   grae->SetPoint(2,350,-0.07827005);
   grae->SetPointError(2,0,0,0.0562285,0.0562285);
   grae->SetPoint(3,425,-0.009007773);
   grae->SetPointError(3,0,0,0.0608713,0.0608713);
   grae->SetPoint(4,525,-0.04370875);
   grae->SetPointError(4,0,0,0.1065396,0.1065396);
   grae->SetPoint(5,625,-0.272785);
   grae->SetPointError(5,0,0,0.1858908,0.1858908);
   grae->SetPoint(6,725,-0.4405969);
   grae->SetPointError(6,0,0,0.3530987,0.3530987);
   grae->SetPoint(7,825,-0.4784584);
   grae->SetPointError(7,0,0,0.5250312,0.5250312);
   grae->SetPoint(8,925,-0.04425773);
   grae->SetPointError(8,0,0,0.5377488,0.5377488);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1247631);
   grae->SetPointError(1,0,0,0.0417249,0.0417249);
   grae->SetPoint(2,350,0.09700737);
   grae->SetPointError(2,0,0,0.05406695,0.05406695);
   grae->SetPoint(3,425,0.07092422);
   grae->SetPointError(3,0,0,0.05589118,0.05589118);
   grae->SetPoint(4,525,-0.07717243);
   grae->SetPointError(4,0,0,0.08811558,0.08811558);
   grae->SetPoint(5,625,0.1369423);
   grae->SetPointError(5,0,0,0.1289094,0.1289094);
   grae->SetPoint(6,725,0.15845);
   grae->SetPointError(6,0,0,0.215432,0.215432);
   grae->SetPoint(7,825,0.11716);
   grae->SetPointError(7,0,0,0.2999268,0.2999268);
   grae->SetPoint(8,925,0.3686635);
   grae->SetPointError(8,0,0,0.2973786,0.2973786);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,0.05530473);
   grae->SetPointError(3,0,0,0.06207301,0.06207301);
   grae->SetPoint(4,525,-0.06996815);
   grae->SetPointError(4,0,0,0.1014,0.1014);
   grae->SetPoint(5,625,0.1249143);
   grae->SetPointError(5,0,0,0.1564522,0.1564522);
   grae->SetPoint(6,725,0.07826812);
   grae->SetPointError(6,0,0,0.2575003,0.2575003);
   grae->SetPoint(7,825,0.007436524);
   grae->SetPointError(7,0,0,0.3805916,0.3805916);
   grae->SetPoint(8,925,0.3296416);
   grae->SetPointError(8,0,0,0.424438,0.424438);

    for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 16 ) {

    // 11/fb, 553p2, >=4, mixture, 8 bins 

    energy = "8";
    lumi = "11.7";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.15);
    syst.push_back(0.30);
    syst.push_back(0.90);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.04791268);
   grae->SetPointError(1,0,0,0.05400354,0.05400354);
   grae->SetPoint(2,350,-0.1359035);
   grae->SetPointError(2,0,0,0.06898996,0.06898996);
   grae->SetPoint(3,425,-0.06227768);
   grae->SetPointError(3,0,0,0.07573852,0.07573852);
   grae->SetPoint(4,525,0.09921159);
   grae->SetPointError(4,0,0,0.09619947,0.09619947);
   grae->SetPoint(5,625,-0.1633746);
   grae->SetPointError(5,0,0,0.1440042,0.1440042);
   grae->SetPoint(6,725,0.2864719);
   grae->SetPointError(6,0,0,0.2403701,0.2403701);
   grae->SetPoint(7,825,0.4843173);
   grae->SetPointError(7,0,0,0.4305148,0.4305148);
   grae->SetPoint(8,925,-0.6627245);
   grae->SetPointError(8,0,0,0.5436823,0.5436823);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.06741467);
   grae->SetPointError(1,0,0,0.05010283,0.05010283);
   grae->SetPoint(2,350,-0.01486992);
   grae->SetPointError(2,0,0,0.06325661,0.06325661);
   grae->SetPoint(3,425,-0.01327023);
   grae->SetPointError(3,0,0,0.06641844,0.06641844);
   grae->SetPoint(4,525,0.137712);
   grae->SetPointError(4,0,0,0.07732728,0.07732728);
   grae->SetPoint(5,625,-0.06690647);
   grae->SetPointError(5,0,0,0.1025216,0.1025216);
   grae->SetPoint(6,725,-0.06720695);
   grae->SetPointError(6,0,0,0.1400445,0.1400445);
   grae->SetPoint(7,825,0.2784852);
   grae->SetPointError(7,0,0,0.2284207,0.2284207);
   grae->SetPoint(8,925,-0.1361178);
   grae->SetPointError(8,0,0,0.1907261,0.1907261);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.06340631);
   grae->SetPointError(1,0,0,0.04887925,0.04887925);
   grae->SetPoint(2,350,0.06423014);
   grae->SetPointError(2,0,0,0.07115013,0.07115013);
   grae->SetPoint(3,425,0.05479067);
   grae->SetPointError(3,0,0,0.07248438,0.07248438);
   grae->SetPoint(4,525,0.009008043);
   grae->SetPointError(4,0,0,0.08464941,0.08464941);
   grae->SetPoint(5,625,0.1961726);
   grae->SetPointError(5,0,0,0.1263657,0.1263657);
   grae->SetPoint(6,725,0.07517212);
   grae->SetPointError(6,0,0,0.1748874,0.1748874);
   grae->SetPoint(7,825,-0.1847061);
   grae->SetPointError(7,0,0,0.2494827,0.2494827);
   grae->SetPoint(8,925,-0.3286244);
   grae->SetPointError(8,0,0,0.2748807,0.2748807);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.01123228);
   grae->SetPointError(1,0,0,0.1093682,0.1093682);
   grae->SetPoint(2,350,-0.1922941);
   grae->SetPointError(2,0,0,0.1576935,0.1576935);
   grae->SetPoint(3,425,0.2561119);
   grae->SetPointError(3,0,0,0.1717025,0.1717025);
   grae->SetPoint(4,525,0.4465904);
   grae->SetPointError(4,0,0,0.1928888,0.1928888);
   grae->SetPoint(5,625,0.2438495);
   grae->SetPointError(5,0,0,0.2269505,0.2269505);
   grae->SetPoint(6,725,-0.3347452);
   grae->SetPointError(6,0,0,0.2961806,0.2961806);
   grae->SetPoint(7,825,-0.1383484);
   grae->SetPointError(7,0,0,0.4131128,0.4131128);
   grae->SetPoint(8,925,0.7123682);
   grae->SetPointError(8,0,0,0.5258358,0.5258358);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,0.04807188);
   grae->SetPointError(3,0,0,0.1790766,0.1790766);
   grae->SetPoint(4,525,0.08529442);
   grae->SetPointError(4,0,0,0.1891667,0.1891667);
   grae->SetPoint(5,625,0.329845);
   grae->SetPointError(5,0,0,0.2840337,0.2840337);
   grae->SetPoint(6,725,-0.3556287);
   grae->SetPointError(6,0,0,0.3657536,0.3657536);
   grae->SetPoint(7,825,-0.2513208);
   grae->SetPointError(7,0,0,0.5037379,0.5037379);
   grae->SetPoint(8,925,2.045471);
   grae->SetPointError(8,0,0,1.68273,1.68273);

    for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 17 ) {

    // 5/fb, 553p2, >=4, mixture, 8 bins 

    energy = "8";
    lumi = "5.1";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.15);
    syst.push_back(0.30);
    syst.push_back(0.90);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.05794724);
   grae->SetPointError(1,0,0,0.03869255,0.03869255);
   grae->SetPoint(2,350,-0.08765549);
   grae->SetPointError(2,0,0,0.05158851,0.05158851);
   grae->SetPoint(3,425,-0.04894145);
   grae->SetPointError(3,0,0,0.05853597,0.05853597);
   grae->SetPoint(4,525,-0.08933298);
   grae->SetPointError(4,0,0,0.1056509,0.1056509);
   grae->SetPoint(5,625,-0.1907294);
   grae->SetPointError(5,0,0,0.1920857,0.1920857);
   grae->SetPoint(6,725,-0.002123071);
   grae->SetPointError(6,0,0,0.3513309,0.3513309);
   grae->SetPoint(7,825,0.05266484);
   grae->SetPointError(7,0,0,0.483124,0.483124);
   grae->SetPoint(8,925,-1);
   grae->SetPointError(8,0,0,0.8500192,0.8500192);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1936001);
   grae->SetPointError(1,0,0,0.0395771,0.0395771);
   grae->SetPoint(2,350,0.1222037);
   grae->SetPointError(2,0,0,0.05096365,0.05096365);
   grae->SetPoint(3,425,0.141719);
   grae->SetPointError(3,0,0,0.05324682,0.05324682);
   grae->SetPoint(4,525,0.1534589);
   grae->SetPointError(4,0,0,0.08864727,0.08864727);
   grae->SetPoint(5,625,0.05650573);
   grae->SetPointError(5,0,0,0.1453681,0.1453681);
   grae->SetPoint(6,725,0.1871141);
   grae->SetPointError(6,0,0,0.2586053,0.2586053);
   grae->SetPoint(7,825,0.2489897);
   grae->SetPointError(7,0,0,0.3277943,0.3277943);
   grae->SetPoint(8,925,-0.4760935);
   grae->SetPointError(8,0,0,0.3289043,0.3289043);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.05254709);
   grae->SetPointError(1,0,0,0.05985427,0.05985427);
   grae->SetPoint(2,350,-0.02945259);
   grae->SetPointError(2,0,0,0.08238973,0.08238973);
   grae->SetPoint(3,425,-0.1391575);
   grae->SetPointError(3,0,0,0.08748967,0.08748967);
   grae->SetPoint(4,525,-0.1088257);
   grae->SetPointError(4,0,0,0.1534943,0.1534943);
   grae->SetPoint(5,625,-0.1844267);
   grae->SetPointError(5,0,0,0.2814228,0.2814228);
   grae->SetPoint(6,725,-0.6294387);
   grae->SetPointError(6,0,0,0.573563,0.573563);
   grae->SetPoint(7,825,-0.06138404);
   grae->SetPointError(7,0,0,0.6280453,0.6280453);
   grae->SetPoint(8,925,0.7750279);
   grae->SetPointError(8,0,0,1.291212,1.291212);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1141267);
   grae->SetPointError(1,0,0,0.05257708,0.05257708);
   grae->SetPoint(2,350,0.1195105);
   grae->SetPointError(2,0,0,0.06938741,0.06938741);
   grae->SetPoint(3,425,0.07689184);
   grae->SetPointError(3,0,0,0.07126281,0.07126281);
   grae->SetPoint(4,525,-0.1350528);
   grae->SetPointError(4,0,0,0.1129661,0.1129661);
   grae->SetPoint(5,625,0.3209068);
   grae->SetPointError(5,0,0,0.1805772,0.1805772);
   grae->SetPoint(6,725,0.4394803);
   grae->SetPointError(6,0,0,0.3214426,0.3214426);
   grae->SetPoint(7,825,0.1106235);
   grae->SetPointError(7,0,0,0.3733761,0.3733761);
   grae->SetPoint(8,925,0.2308933);
   grae->SetPointError(8,0,0,0.3750249,0.3750249);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,0.09470716);
   grae->SetPointError(3,0,0,0.07559551,0.07559551);
   grae->SetPoint(4,525,-0.0850152);
   grae->SetPointError(4,0,0,0.1228471,0.1228471);
   grae->SetPoint(5,625,0.2252092);
   grae->SetPointError(5,0,0,0.1904266,0.1904266);
   grae->SetPoint(6,725,0.1789309);
   grae->SetPointError(6,0,0,0.3105728,0.3105728);
   grae->SetPoint(7,825,0.1972394);
   grae->SetPointError(7,0,0,0.4567579,0.4567579);
   grae->SetPoint(8,925,0.2701397);
   grae->SetPointError(8,0,0,0.4828748,0.4828748);

    for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 18 ) {

    // 5/fb, 553p2, >=4, mixture, 8 bins 

    energy = "8";
    lumi = "5.1";

    int nhistos = 5;
    
    syst.clear();
    syst.push_back(0.15);
    syst.push_back(0.30);
    syst.push_back(0.90);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.001290323);
   grae->SetPointError(1,0,0,0.06889915,0.06889915);
   grae->SetPoint(2,350,-0.0781273);
   grae->SetPointError(2,0,0,0.09848256,0.09848256);
   grae->SetPoint(3,425,-0.02603525);
   grae->SetPointError(3,0,0,0.1082178,0.1082178);
   grae->SetPoint(4,525,0.1450656);
   grae->SetPointError(4,0,0,0.1384182,0.1384182);
   grae->SetPoint(5,625,-0.2261092);
   grae->SetPointError(5,0,0,0.2098671,0.2098671);
   grae->SetPoint(6,725,0.07374993);
   grae->SetPointError(6,0,0,0.3246947,0.3246947);
   grae->SetPoint(7,825,0.5598106);
   grae->SetPointError(7,0,0,0.6213942,0.6213942);
   grae->SetPoint(8,925,-0.6211451);
   grae->SetPointError(8,0,0,0.7662211,0.7662211);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.03986036);
   grae->SetPointError(1,0,0,0.06541399,0.06541399);
   grae->SetPoint(2,350,-0.04029156);
   grae->SetPointError(2,0,0,0.0898567,0.0898567);
   grae->SetPoint(3,425,0.01404287);
   grae->SetPointError(3,0,0,0.09293667,0.09293667);
   grae->SetPoint(4,525,0.1556405);
   grae->SetPointError(4,0,0,0.1119939,0.1119939);
   grae->SetPoint(5,625,0.004789268);
   grae->SetPointError(5,0,0,0.146456,0.146456);
   grae->SetPoint(6,725,-0.0230194);
   grae->SetPointError(6,0,0,0.197564,0.197564);
   grae->SetPoint(7,825,0.5241694);
   grae->SetPointError(7,0,0,0.3499499,0.3499499);
   grae->SetPoint(8,925,-0.1878883);
   grae->SetPointError(8,0,0,0.2741571,0.2741571);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.06991643);
   grae->SetPointError(1,0,0,0.06867174,0.06867174);
   grae->SetPoint(2,350,0.06399421);
   grae->SetPointError(2,0,0,0.1035835,0.1035835);
   grae->SetPoint(3,425,0.005787923);
   grae->SetPointError(3,0,0,0.1035099,0.1035099);
   grae->SetPoint(4,525,0.12509);
   grae->SetPointError(4,0,0,0.122883,0.122883);
   grae->SetPoint(5,625,0.1026806);
   grae->SetPointError(5,0,0,0.1749495,0.1749495);
   grae->SetPoint(6,725,0.1255689);
   grae->SetPointError(6,0,0,0.2469546,0.2469546);
   grae->SetPoint(7,825,-0.4618831);
   grae->SetPointError(7,0,0,0.3526918,0.3526918);
   grae->SetPoint(8,925,-0.241596);
   grae->SetPointError(8,0,0,0.40197,0.40197);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.07162308);
   grae->SetPointError(1,0,0,0.1402287,0.1402287);
   grae->SetPoint(2,350,-0.2166778);
   grae->SetPointError(2,0,0,0.2063398,0.2063398);
   grae->SetPoint(3,425,0.09872895);
   grae->SetPointError(3,0,0,0.2210173,0.2210173);
   grae->SetPoint(4,525,0.2704021);
   grae->SetPointError(4,0,0,0.242335,0.242335);
   grae->SetPoint(5,625,0.4121727);
   grae->SetPointError(5,0,0,0.3153333,0.3153333);
   grae->SetPoint(6,725,-0.3306592);
   grae->SetPointError(6,0,0,0.3850166,0.3850166);
   grae->SetPoint(7,825,-0.5974379);
   grae->SetPointError(7,0,0,0.6170791,0.6170791);
   grae->SetPoint(8,925,0.3737041);
   grae->SetPointError(8,0,0,0.6260212,0.6260212);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,-0.08633713);
   grae->SetPointError(3,0,0,0.2211158,0.2211158);
   grae->SetPoint(4,525,-0.03019842);
   grae->SetPointError(4,0,0,0.2304689,0.2304689);
   grae->SetPoint(5,625,0.5144609);
   grae->SetPointError(5,0,0,0.3724645,0.3724645);
   grae->SetPoint(6,725,-0.3151908);
   grae->SetPointError(6,0,0,0.436886,0.436886);
   grae->SetPoint(7,825,-0.6209835);
   grae->SetPointError(7,0,0,0.6914305,0.6914305);
   grae->SetPoint(8,925,1.476955);
   grae->SetPointError(8,0,0,1.619158,1.619158);

    for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 119 ) {

    energy = "8";
    lumi = "11.7";

    int nhistos = 6;
    
    syst.clear();
    syst.push_back(0.15);
    syst.push_back(0.30);
    syst.push_back(0.90);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets, 0 b-tags");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.0853453);
   grae->SetPointError(1,0,0,0.05307479,0.05307479);
   grae->SetPoint(2,350,0.004181127);
   grae->SetPointError(2,0,0,0.06688684,0.06688684);
   grae->SetPoint(3,425,-0.09644717);
   grae->SetPointError(3,0,0,0.06762834,0.06762834);
   grae->SetPoint(4,525,-0.1593712);
   grae->SetPointError(4,0,0,0.07880235,0.07880235);
   grae->SetPoint(5,625,0.02315273);
   grae->SetPointError(5,0,0,0.1139954,0.1139954);
   grae->SetPoint(6,725,0.5167808);
   grae->SetPointError(6,0,0,0.2063214,0.2063214);
   grae->SetPoint(7,825,-0.1365756);
   grae->SetPointError(7,0,0,0.22034,0.22034);
   grae->SetPoint(8,925,0.005151018);
   grae->SetPointError(8,0,0,0.3171721,0.3171721);
   grae->SetPoint(9,1025,-0.3373959);
   grae->SetPointError(9,0,0,0.3798285,0.3798285);
   grae->SetPoint(10,1125,0.06581228);
   grae->SetPointError(10,0,0,0.3177888,0.3177888);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets, #geq1 b-tags");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.079952);
   grae->SetPointError(1,0,0,0.04267817,0.04267817);
   grae->SetPoint(2,350,-0.1173302);
   grae->SetPointError(2,0,0,0.05969023,0.05969023);
   grae->SetPoint(3,425,-0.1591889);
   grae->SetPointError(3,0,0,0.06093977,0.06093977);
   grae->SetPoint(4,525,-0.09722469);
   grae->SetPointError(4,0,0,0.08594708,0.08594708);
   grae->SetPoint(5,625,0.016354);
   grae->SetPointError(5,0,0,0.1419724,0.1419724);
   grae->SetPoint(6,725,0.5406711);
   grae->SetPointError(6,0,0,0.2815224,0.2815224);
   grae->SetPoint(7,825,-0.09389022);
   grae->SetPointError(7,0,0,0.2891103,0.2891103);
   grae->SetPoint(8,925,-0.1315193);
   grae->SetPointError(8,0,0,0.5301807,0.5301807);
   grae->SetPoint(9,1025,0.6352715);
   grae->SetPointError(9,0,0,1.004588,1.004588);
   grae->SetPoint(10,1125,0.6698436);
   grae->SetPointError(10,0,0,1.006486,1.006486);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#gamma + jets, 0 b-tags");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,0,0);
   grae->SetPointError(1,0,0,0,0);
   grae->SetPoint(2,0,0);
   grae->SetPointError(2,0,0,0,0);
   grae->SetPoint(3,425,0.1079193);
   grae->SetPointError(3,0,0,0.1075442,0.1075442);
   grae->SetPoint(4,525,0.1707858);
   grae->SetPointError(4,0,0,0.1364307,0.1364307);
   grae->SetPoint(5,625,-0.05785309);
   grae->SetPointError(5,0,0,0.1962089,0.1962089);
   grae->SetPoint(6,725,0.2184542);
   grae->SetPointError(6,0,0,0.3218333,0.3218333);
   grae->SetPoint(7,825,-0.1745394);
   grae->SetPointError(7,0,0,0.4719121,0.4719121);
   grae->SetPoint(8,925,-0.4530208);
   grae->SetPointError(8,0,0,0.7738448,0.7738448);
   grae->SetPoint(9,1025,-0.4863946);
   grae->SetPointError(9,0,0,1.26975,1.26975);
   grae->SetPoint(10,1125,-1);
   grae->SetPointError(10,0,0,2.467818,2.467818);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#gamma + jets, #geq1 b-tags");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,0,0);
   grae->SetPointError(1,0,0,0,0);
   grae->SetPoint(2,0,0);
   grae->SetPointError(2,0,0,0,0);
   grae->SetPoint(3,425,-0.05868979);
   grae->SetPointError(3,0,0,0.2115641,0.2115641);
   grae->SetPoint(4,525,0.7703572);
   grae->SetPointError(4,0,0,0.3772374,0.3772374);
   grae->SetPoint(5,625,-0.3743926);
   grae->SetPointError(5,0,0,0.4263535,0.4263535);
   grae->SetPoint(6,725,0.1005332);
   grae->SetPointError(6,0,0,0.5689451,0.5689451);
   grae->SetPoint(7,825,0.8183261);
   grae->SetPointError(7,0,0,1.010632,1.010632);
   grae->SetPoint(8,925,-0.2054984);
   grae->SetPointError(8,0,0,1.218292,1.218292);
   grae->SetPoint(9,1025,-1);
   grae->SetPointError(9,0,0,1.782282,1.782282);
   grae->SetPoint(10,1125,-1);
   grae->SetPointError(10,0,0,0,0);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu#mu + jets, 0 b-tags");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.0531029);
   grae->SetPointError(1,0,0,0.1653923,0.1653923);
   grae->SetPoint(2,350,-0.3250542);
   grae->SetPointError(2,0,0,0.2455483,0.2455483);
   grae->SetPoint(3,425,-0.1302368);
   grae->SetPointError(3,0,0,0.2629382,0.2629382);
   grae->SetPoint(4,525,0.5422405);
   grae->SetPointError(4,0,0,0.3597825,0.3597825);
   grae->SetPoint(5,625,0.2499893);
   grae->SetPointError(5,0,0,0.3543421,0.3543421);
   grae->SetPoint(6,725,-0.1856812);
   grae->SetPointError(6,0,0,0.5119455,0.5119455);
   grae->SetPoint(7,825,-0.5724379);
   grae->SetPointError(7,0,0,0.7610719,0.7610719);
   grae->SetPoint(8,925,-1);
   grae->SetPointError(8,0,0,1.588015,1.588015);
   grae->SetPoint(9,1025,0.8925247);
   grae->SetPointError(9,0,0,2.273729,2.273729);
   grae->SetPoint(10,1125,-1);
   grae->SetPointError(10,0,0,1.894768,1.894768);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu#mu + jets, #geq1 b-tags");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.143902);
   grae->SetPointError(1,0,0,0.282208,0.282208);
   grae->SetPoint(2,350,-0.3835099);
   grae->SetPointError(2,0,0,0.3695968,0.3695968);
   grae->SetPoint(3,425,-0.1252243);
   grae->SetPointError(3,0,0,0.3512171,0.3512171);
   grae->SetPoint(4,525,-0.2119804);
   grae->SetPointError(4,0,0,0.4420018,0.4420018);
   grae->SetPoint(5,625,-0.3498493);
   grae->SetPointError(5,0,0,0.6704157,0.6704157);
   grae->SetPoint(6,725,-0.6303498);
   grae->SetPointError(6,0,0,1.054553,1.054553);
   grae->SetPoint(7,825,-1);
   grae->SetPointError(7,0,0,2.006929,2.006929);
   grae->SetPoint(8,925,0.1349797);
   grae->SetPointError(8,0,0,1.234297,1.234297);
   grae->SetPoint(9,1025,-1);
   grae->SetPointError(9,0,0,0,0);
   grae->SetPoint(10,1125,-1);
   grae->SetPointError(10,0,0,0,0);

    for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 120 ) {

    energy = "8";
    lumi = "11.7";

    int nhistos = 3;
    
    syst.clear();
    syst.push_back(0.15);
    syst.push_back(0.30);
    syst.push_back(0.90);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.03906917);
   grae->SetPointError(1,0,0,0.03353366,0.03353366);
   grae->SetPoint(2,350,-0.03533564);
   grae->SetPointError(2,0,0,0.04383138,0.04383138);
   grae->SetPoint(3,425,-0.1056004);
   grae->SetPointError(3,0,0,0.04400892,0.04400892);
   grae->SetPoint(4,525,-0.07988218);
   grae->SetPointError(4,0,0,0.05497124,0.05497124);
   grae->SetPoint(5,625,0.03079012);
   grae->SetPointError(5,0,0,0.08253733,0.08253733);
   grae->SetPoint(6,725,0.5489098);
   grae->SetPointError(6,0,0,0.1549391,0.1549391);
   grae->SetPoint(7,825,-0.03631155);
   grae->SetPointError(7,0,0,0.162666,0.162666);
   grae->SetPoint(8,925,-0.098666);
   grae->SetPointError(8,0,0,0.2480876,0.2480876);
   grae->SetPoint(9,1025,-0.2785468);
   grae->SetPointError(9,0,0,0.3058,0.3058);
   grae->SetPoint(10,1125,-0.05397235);
   grae->SetPointError(10,0,0,0.2634776,0.2634776);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#gamma + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,0,0);
   grae->SetPointError(1,0,0,0,0);
   grae->SetPoint(2,0,0);
   grae->SetPointError(2,0,0,0,0);
   grae->SetPoint(3,425,0.09278029);
   grae->SetPointError(3,0,0,0.09664966,0.09664966);
   grae->SetPoint(4,525,0.2739804);
   grae->SetPointError(4,0,0,0.1275042,0.1275042);
   grae->SetPoint(5,625,-0.1090606);
   grae->SetPointError(5,0,0,0.1778259,0.1778259);
   grae->SetPoint(6,725,0.2388484);
   grae->SetPointError(6,0,0,0.2877924,0.2877924);
   grae->SetPoint(7,825,0.1027996);
   grae->SetPointError(7,0,0,0.4100296,0.4100296);
   grae->SetPoint(8,925,-0.3708192);
   grae->SetPointError(8,0,0,0.6611959,0.6611959);
   grae->SetPoint(9,1025,-0.6806142);
   grae->SetPointError(9,0,0,1.133479,1.133479);
   grae->SetPoint(10,1125,-1);
   grae->SetPointError(10,0,0,2.386582,2.386582);
    

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,2.793741e-06);
   grae->SetPointError(1,0,0,0.1419687,0.1419687);
   grae->SetPoint(2,350,-0.3251202);
   grae->SetPointError(2,0,0,0.2083702,0.2083702);
   grae->SetPoint(3,425,-0.08794268);
   grae->SetPointError(3,0,0,0.213224,0.213224);
   grae->SetPoint(4,525,0.3507572);
   grae->SetPointError(4,0,0,0.2762535,0.2762535);
   grae->SetPoint(5,625,0.1011718);
   grae->SetPointError(5,0,0,0.3006639,0.3006639);
   grae->SetPoint(6,725,-0.28048);
   grae->SetPointError(6,0,0,0.4577014,0.4577014);
   grae->SetPoint(7,825,-0.6510373);
   grae->SetPointError(7,0,0,0.7113626,0.7113626);
   grae->SetPoint(8,925,-0.2464512);
   grae->SetPointError(8,0,0,0.9219985,0.9219985);
   grae->SetPoint(9,1025,0.503977);
   grae->SetPointError(9,0,0,1.695403,1.695403);
   grae->SetPoint(10,1125,0.8699347);
   grae->SetPointError(10,0,0,1.68778,1.68778);

    for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 19 ) {

    // 11/fb, 553p2, >=4, mixture, 8 bins, no gamma/mumu

    energy = "8";
    lumi = "11.7";

    int nhistos = 4;
    
    syst.clear();
    syst.push_back(0.15);
    syst.push_back(0.30);
    syst.push_back(0.90);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.04791268);
   grae->SetPointError(1,0,0,0.05400354,0.05400354);
   grae->SetPoint(2,350,-0.1359035);
   grae->SetPointError(2,0,0,0.06898996,0.06898996);
   grae->SetPoint(3,425,-0.06227768);
   grae->SetPointError(3,0,0,0.07573852,0.07573852);
   grae->SetPoint(4,525,0.09921159);
   grae->SetPointError(4,0,0,0.09619947,0.09619947);
   grae->SetPoint(5,625,-0.1633746);
   grae->SetPointError(5,0,0,0.1440042,0.1440042);
   grae->SetPoint(6,725,0.2864719);
   grae->SetPointError(6,0,0,0.2403701,0.2403701);
   grae->SetPoint(7,825,0.4843173);
   grae->SetPointError(7,0,0,0.4305148,0.4305148);
   grae->SetPoint(8,925,-0.6627245);
   grae->SetPointError(8,0,0,0.5436823,0.5436823);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.06741467);
   grae->SetPointError(1,0,0,0.05010283,0.05010283);
   grae->SetPoint(2,350,-0.01486992);
   grae->SetPointError(2,0,0,0.06325661,0.06325661);
   grae->SetPoint(3,425,-0.01327023);
   grae->SetPointError(3,0,0,0.06641844,0.06641844);
   grae->SetPoint(4,525,0.137712);
   grae->SetPointError(4,0,0,0.07732728,0.07732728);
   grae->SetPoint(5,625,-0.06690647);
   grae->SetPointError(5,0,0,0.1025216,0.1025216);
   grae->SetPoint(6,725,-0.06720695);
   grae->SetPointError(6,0,0,0.1400445,0.1400445);
   grae->SetPoint(7,825,0.2784852);
   grae->SetPointError(7,0,0,0.2284207,0.2284207);
   grae->SetPoint(8,925,-0.1361178);
   grae->SetPointError(8,0,0,0.1907261,0.1907261);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.06340631);
   grae->SetPointError(1,0,0,0.04887925,0.04887925);
   grae->SetPoint(2,350,0.06423014);
   grae->SetPointError(2,0,0,0.07115013,0.07115013);
   grae->SetPoint(3,425,0.05479067);
   grae->SetPointError(3,0,0,0.07248438,0.07248438);
   grae->SetPoint(4,525,0.009008043);
   grae->SetPointError(4,0,0,0.08464941,0.08464941);
   grae->SetPoint(5,625,0.1961726);
   grae->SetPointError(5,0,0,0.1263657,0.1263657);
   grae->SetPoint(6,725,0.07517212);
   grae->SetPointError(6,0,0,0.1748874,0.1748874);
   grae->SetPoint(7,825,-0.1847061);
   grae->SetPointError(7,0,0,0.2494827,0.2494827);
   grae->SetPoint(8,925,-0.3286244);
   grae->SetPointError(8,0,0,0.2748807,0.2748807);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.01123228);
   grae->SetPointError(1,0,0,0.1093682,0.1093682);
   grae->SetPoint(2,350,-0.1922941);
   grae->SetPointError(2,0,0,0.1576935,0.1576935);
   grae->SetPoint(3,425,0.2561119);
   grae->SetPointError(3,0,0,0.1717025,0.1717025);
   grae->SetPoint(4,525,0.4465904);
   grae->SetPointError(4,0,0,0.1928888,0.1928888);
   grae->SetPoint(5,625,0.2438495);
   grae->SetPointError(5,0,0,0.2269505,0.2269505);
   grae->SetPoint(6,725,-0.3347452);
   grae->SetPointError(6,0,0,0.2961806,0.2961806);
   grae->SetPoint(7,825,-0.1383484);
   grae->SetPointError(7,0,0,0.4131128,0.4131128);
   grae->SetPoint(8,925,0.7123682);
   grae->SetPointError(8,0,0,0.5258358,0.5258358);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

//     // mumu->gamma, no aT
//     grae = new TGraphAsymmErrors(nbins);
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
//    grae->SetPoint(3,425,0.04807188);
//    grae->SetPointError(3,0,0,0.1790766,0.1790766);
//    grae->SetPoint(4,525,0.08529442);
//    grae->SetPointError(4,0,0,0.1891667,0.1891667);
//    grae->SetPoint(5,625,0.329845);
//    grae->SetPointError(5,0,0,0.2840337,0.2840337);
//    grae->SetPoint(6,725,-0.3556287);
//    grae->SetPointError(6,0,0,0.3657536,0.3657536);
//    grae->SetPoint(7,825,-0.2513208);
//    grae->SetPointError(7,0,0,0.5037379,0.5037379);
//    grae->SetPoint(8,925,2.045471);
//    grae->SetPointError(8,0,0,1.68273,1.68273);

//     for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 20 ) {

    // 11/fb, 553p2, 2-3, mixture, 8 bins 

    energy = "8";
    lumi = "11.7";

    int nhistos = 8;
    
    syst.clear();
    syst.push_back(0.10);
    syst.push_back(0.20);
    syst.push_back(0.20);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.0546802);
   grae->SetPointError(1,0,0,0.03025781,0.03025781);
   grae->SetPoint(2,350,-0.06052873);
   grae->SetPointError(2,0,0,0.03694661,0.03694661);
   grae->SetPoint(3,425,-0.03560808);
   grae->SetPointError(3,0,0,0.04184187,0.04184187);
   grae->SetPoint(4,525,-0.07427342);
   grae->SetPointError(4,0,0,0.07445687,0.07445687);
   grae->SetPoint(5,625,-0.1416568);
   grae->SetPointError(5,0,0,0.1269872,0.1269872);
   grae->SetPoint(6,725,0.02794661);
   grae->SetPointError(6,0,0,0.2256453,0.2256453);
   grae->SetPoint(7,825,-0.3371197);
   grae->SetPointError(7,0,0,0.3771194,0.3771194);
   grae->SetPoint(8,925,-0.5331764);
   grae->SetPointError(8,0,0,0.4782493,0.4782493);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.160398);
   grae->SetPointError(1,0,0,0.02885296,0.02885296);
   grae->SetPoint(2,350,0.1576586);
   grae->SetPointError(2,0,0,0.03643444,0.03643444);
   grae->SetPoint(3,425,0.1234905);
   grae->SetPointError(3,0,0,0.03802028,0.03802028);
   grae->SetPoint(4,525,0.1485403);
   grae->SetPointError(4,0,0,0.06239629,0.06239629);
   grae->SetPoint(5,625,0.07797849);
   grae->SetPointError(5,0,0,0.09611793,0.09611793);
   grae->SetPoint(6,725,0.04869653);
   grae->SetPointError(6,0,0,0.1604718,0.1604718);
   grae->SetPoint(7,825,0.1536973);
   grae->SetPointError(7,0,0,0.2358036,0.2358036);
   grae->SetPoint(8,925,-0.04202397);
   grae->SetPointError(8,0,0,0.2119533,0.2119533);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.04383);
   grae->SetPointError(1,0,0,0.04206461,0.04206461);
   grae->SetPoint(2,350,-0.07827005);
   grae->SetPointError(2,0,0,0.0562285,0.0562285);
   grae->SetPoint(3,425,-0.009007773);
   grae->SetPointError(3,0,0,0.0608713,0.0608713);
   grae->SetPoint(4,525,-0.04370875);
   grae->SetPointError(4,0,0,0.1065396,0.1065396);
   grae->SetPoint(5,625,-0.272785);
   grae->SetPointError(5,0,0,0.1858908,0.1858908);
   grae->SetPoint(6,725,-0.4405969);
   grae->SetPointError(6,0,0,0.3530987,0.3530987);
   grae->SetPoint(7,825,-0.4784584);
   grae->SetPointError(7,0,0,0.5250312,0.5250312);
   grae->SetPoint(8,925,-0.04425773);
   grae->SetPointError(8,0,0,0.5377488,0.5377488);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1247631);
   grae->SetPointError(1,0,0,0.0417249,0.0417249);
   grae->SetPoint(2,350,0.09700737);
   grae->SetPointError(2,0,0,0.05406695,0.05406695);
   grae->SetPoint(3,425,0.07092422);
   grae->SetPointError(3,0,0,0.05589118,0.05589118);
   grae->SetPoint(4,525,-0.07717243);
   grae->SetPointError(4,0,0,0.08811558,0.08811558);
   grae->SetPoint(5,625,0.1369423);
   grae->SetPointError(5,0,0,0.1289094,0.1289094);
   grae->SetPoint(6,725,0.15845);
   grae->SetPointError(6,0,0,0.215432,0.215432);
   grae->SetPoint(7,825,0.11716);
   grae->SetPointError(7,0,0,0.2999268,0.2999268);
   grae->SetPoint(8,925,0.3686635);
   grae->SetPointError(8,0,0,0.2973786,0.2973786);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,0.05530473);
   grae->SetPointError(3,0,0,0.06207301,0.06207301);
   grae->SetPoint(4,525,-0.06996815);
   grae->SetPointError(4,0,0,0.1014,0.1014);
   grae->SetPoint(5,625,0.1249143);
   grae->SetPointError(5,0,0,0.1564522,0.1564522);
   grae->SetPoint(6,725,0.07826812);
   grae->SetPointError(6,0,0,0.2575003,0.2575003);
   grae->SetPoint(7,825,0.007436524);
   grae->SetPointError(7,0,0,0.3805916,0.3805916);
   grae->SetPoint(8,925,0.3296416);
   grae->SetPointError(8,0,0,0.424438,0.424438);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.03906917);
   grae->SetPointError(1,0,0,0.03353366,0.03353366);
   grae->SetPoint(2,350,-0.03533564);
   grae->SetPointError(2,0,0,0.04383138,0.04383138);
   grae->SetPoint(3,425,-0.1056004);
   grae->SetPointError(3,0,0,0.04400892,0.04400892);
   grae->SetPoint(4,525,-0.07988218);
   grae->SetPointError(4,0,0,0.05497124,0.05497124);
   grae->SetPoint(5,625,0.03079012);
   grae->SetPointError(5,0,0,0.08253733,0.08253733);
   grae->SetPoint(6,725,0.5489098);
   grae->SetPointError(6,0,0,0.1549391,0.1549391);
   grae->SetPoint(7,825,-0.03631155);
   grae->SetPointError(7,0,0,0.162666,0.162666);
   grae->SetPoint(8,925,-0.098666);
   grae->SetPointError(8,0,0,0.2480876,0.2480876);
   grae->SetPoint(9,1025,-0.2785468);
   grae->SetPointError(9,0,0,0.3058,0.3058);
   grae->SetPoint(10,1125,-0.05397235);
   grae->SetPointError(10,0,0,0.2634776,0.2634776);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#gamma + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,0,0);
   grae->SetPointError(1,0,0,0,0);
   grae->SetPoint(2,0,0);
   grae->SetPointError(2,0,0,0,0);
   grae->SetPoint(3,425,0.09278029);
   grae->SetPointError(3,0,0,0.09664966,0.09664966);
   grae->SetPoint(4,525,0.2739804);
   grae->SetPointError(4,0,0,0.1275042,0.1275042);
   grae->SetPoint(5,625,-0.1090606);
   grae->SetPointError(5,0,0,0.1778259,0.1778259);
   grae->SetPoint(6,725,0.2388484);
   grae->SetPointError(6,0,0,0.2877924,0.2877924);
   grae->SetPoint(7,825,0.1027996);
   grae->SetPointError(7,0,0,0.4100296,0.4100296);
   grae->SetPoint(8,925,-0.3708192);
   grae->SetPointError(8,0,0,0.6611959,0.6611959);
   grae->SetPoint(9,1025,-0.6806142);
   grae->SetPointError(9,0,0,1.133479,1.133479);
   grae->SetPoint(10,1125,-1);
   grae->SetPointError(10,0,0,2.386582,2.386582);
    

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,2.793741e-06);
   grae->SetPointError(1,0,0,0.1419687,0.1419687);
   grae->SetPoint(2,350,-0.3251202);
   grae->SetPointError(2,0,0,0.2083702,0.2083702);
   grae->SetPoint(3,425,-0.08794268);
   grae->SetPointError(3,0,0,0.213224,0.213224);
   grae->SetPoint(4,525,0.3507572);
   grae->SetPointError(4,0,0,0.2762535,0.2762535);
   grae->SetPoint(5,625,0.1011718);
   grae->SetPointError(5,0,0,0.3006639,0.3006639);
   grae->SetPoint(6,725,-0.28048);
   grae->SetPointError(6,0,0,0.4577014,0.4577014);
   grae->SetPoint(7,825,-0.6510373);
   grae->SetPointError(7,0,0,0.7113626,0.7113626);
   grae->SetPoint(8,925,-0.2464512);
   grae->SetPointError(8,0,0,0.9219985,0.9219985);
   grae->SetPoint(9,1025,0.503977);
   grae->SetPointError(9,0,0,1.695403,1.695403);
   grae->SetPoint(10,1125,0.8699347);
   grae->SetPointError(10,0,0,1.68778,1.68778);

    for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 21 ) {

    // 11/fb, 553p2, >=4, mixture, 8 bins

    energy = "8";
    lumi = "11.7";

    int nhistos = 8;
    
    syst.clear();
    syst.push_back(0.10);
    syst.push_back(0.20);
    syst.push_back(0.30);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.04791268);
   grae->SetPointError(1,0,0,0.05400354,0.05400354);
   grae->SetPoint(2,350,-0.1359035);
   grae->SetPointError(2,0,0,0.06898996,0.06898996);
   grae->SetPoint(3,425,-0.06227768);
   grae->SetPointError(3,0,0,0.07573852,0.07573852);
   grae->SetPoint(4,525,0.09921159);
   grae->SetPointError(4,0,0,0.09619947,0.09619947);
   grae->SetPoint(5,625,-0.1633746);
   grae->SetPointError(5,0,0,0.1440042,0.1440042);
   grae->SetPoint(6,725,0.2864719);
   grae->SetPointError(6,0,0,0.2403701,0.2403701);
   grae->SetPoint(7,825,0.4843173);
   grae->SetPointError(7,0,0,0.4305148,0.4305148);
   grae->SetPoint(8,925,-0.6627245);
   grae->SetPointError(8,0,0,0.5436823,0.5436823);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.06741467);
   grae->SetPointError(1,0,0,0.05010283,0.05010283);
   grae->SetPoint(2,350,-0.01486992);
   grae->SetPointError(2,0,0,0.06325661,0.06325661);
   grae->SetPoint(3,425,-0.01327023);
   grae->SetPointError(3,0,0,0.06641844,0.06641844);
   grae->SetPoint(4,525,0.137712);
   grae->SetPointError(4,0,0,0.07732728,0.07732728);
   grae->SetPoint(5,625,-0.06690647);
   grae->SetPointError(5,0,0,0.1025216,0.1025216);
   grae->SetPoint(6,725,-0.06720695);
   grae->SetPointError(6,0,0,0.1400445,0.1400445);
   grae->SetPoint(7,825,0.2784852);
   grae->SetPointError(7,0,0,0.2284207,0.2284207);
   grae->SetPoint(8,925,-0.1361178);
   grae->SetPointError(8,0,0,0.1907261,0.1907261);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;
 
      // mu, no aT: 1->2
      grae = new TGraphAsymmErrors(nbins);
      grae->SetName("Graph");
      grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
      grae->SetFillColor(1);
      grae->SetLineWidth(3);
      grae->SetMarkerStyle(20);
      grae->SetMarkerSize(1.5);
     grae->SetPoint(1,300,-0.06340631);
     grae->SetPointError(1,0,0,0.04887925,0.04887925);
     grae->SetPoint(2,350,0.06423014);
     grae->SetPointError(2,0,0,0.07115013,0.07115013);
     grae->SetPoint(3,425,0.05479067);
     grae->SetPointError(3,0,0,0.07248438,0.07248438);
     grae->SetPoint(4,525,0.009008043);
     grae->SetPointError(4,0,0,0.08464941,0.08464941);
     grae->SetPoint(5,625,0.1961726);
     grae->SetPointError(5,0,0,0.1263657,0.1263657);
     grae->SetPoint(6,725,0.07517212);
     grae->SetPointError(6,0,0,0.1748874,0.1748874);
     grae->SetPoint(7,825,-0.1847061);
     grae->SetPointError(7,0,0,0.2494827,0.2494827);
     grae->SetPoint(8,925,-0.3286244);
     grae->SetPointError(8,0,0,0.2748807,0.2748807);
   
      for ( uint i = 0; i < nbins; ++i ) {
        double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
        errors[index][i] = grae->GetErrorYhigh(i+1);
      }
      titles.push_back( grae->GetTitle() );
      index++;
 
      // mu->mumu, no aT
      grae = new TGraphAsymmErrors(nbins);
      grae->SetName("Graph");
      grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
      grae->SetFillColor(1);
      grae->SetLineWidth(3);
      grae->SetMarkerStyle(20);
      grae->SetMarkerSize(1.5);
     grae->SetPoint(1,300,-0.01123228);
     grae->SetPointError(1,0,0,0.1093682,0.1093682);
     grae->SetPoint(2,350,-0.1922941);
     grae->SetPointError(2,0,0,0.1576935,0.1576935);
     grae->SetPoint(3,425,0.2561119);
     grae->SetPointError(3,0,0,0.1717025,0.1717025);
     grae->SetPoint(4,525,0.4465904);
     grae->SetPointError(4,0,0,0.1928888,0.1928888);
     grae->SetPoint(5,625,0.2438495);
     grae->SetPointError(5,0,0,0.2269505,0.2269505);
     grae->SetPoint(6,725,-0.3347452);
     grae->SetPointError(6,0,0,0.2961806,0.2961806);
     grae->SetPoint(7,825,-0.1383484);
     grae->SetPointError(7,0,0,0.4131128,0.4131128);
     grae->SetPoint(8,925,0.7123682);
     grae->SetPointError(8,0,0,0.5258358,0.5258358);
 
      for ( uint i = 0; i < nbins; ++i ) {
        double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
        errors[index][i] = grae->GetErrorYhigh(i+1);
      }
      titles.push_back( grae->GetTitle() );
      index++;
 
      // mumu->gamma, no aT
      grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,-0.03868952);
   grae->SetPointError(3,0,0,0.1660727,0.1660727);
   grae->SetPoint(4,525,0.0774636);
   grae->SetPointError(4,0,0,0.1781558,0.1781558);
   grae->SetPoint(5,625,0.1110814);
   grae->SetPointError(5,0,0,0.2383428,0.2383428);
   grae->SetPoint(6,725,-0.4092417);
   grae->SetPointError(6,0,0,0.3349802,0.3349802);
   grae->SetPoint(7,825,-0.4450173);
   grae->SetPointError(7,0,0,0.4535583,0.4535583);
   grae->SetPoint(8,925,0.06917934);
   grae->SetPointError(8,0,0,0.445336,0.445336);

//      grae->SetPoint(1,0,0);
//      grae->SetPointError(1,0,0,0,0);
//      grae->SetPoint(2,0,0);
//      grae->SetPointError(2,0,0,0,0);
//      grae->SetPoint(3,425,0.04807188);
//      grae->SetPointError(3,0,0,0.1790766,0.1790766);
//      grae->SetPoint(4,525,0.08529442);
//      grae->SetPointError(4,0,0,0.1891667,0.1891667);
//      grae->SetPoint(5,625,0.329845);
//      grae->SetPointError(5,0,0,0.2840337,0.2840337);
//      grae->SetPoint(6,725,-0.3556287);
//      grae->SetPointError(6,0,0,0.3657536,0.3657536);
//      grae->SetPoint(7,825,-0.2513208);
//      grae->SetPointError(7,0,0,0.5037379,0.5037379);
//      grae->SetPoint(8,925,2.045471);
//      grae->SetPointError(8,0,0,1.68273,1.68273);
 
      for ( uint i = 0; i < nbins; ++i ) {
        double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
        if ( x < 275. ) { y = -1000.; }
        values[index][i] = y;
        errors[index][i] = grae->GetErrorYhigh(i+1);
      }
      titles.push_back( grae->GetTitle() );
      index++;
 
      grae = new TGraphAsymmErrors(nbins);
      grae->SetName("Graph");
      grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu + jets)");
      grae->SetFillColor(1);
      grae->SetLineWidth(3);
      grae->SetMarkerStyle(20);
      grae->SetMarkerSize(1.5);
     grae->SetPoint(1,300,0.03906917);
     grae->SetPointError(1,0,0,0.03353366,0.03353366);
     grae->SetPoint(2,350,-0.03533564);
     grae->SetPointError(2,0,0,0.04383138,0.04383138);
     grae->SetPoint(3,425,-0.1056004);
     grae->SetPointError(3,0,0,0.04400892,0.04400892);
     grae->SetPoint(4,525,-0.07988218);
     grae->SetPointError(4,0,0,0.05497124,0.05497124);
     grae->SetPoint(5,625,0.03079012);
     grae->SetPointError(5,0,0,0.08253733,0.08253733);
     grae->SetPoint(6,725,0.5489098);
     grae->SetPointError(6,0,0,0.1549391,0.1549391);
     grae->SetPoint(7,825,-0.03631155);
     grae->SetPointError(7,0,0,0.162666,0.162666);
     grae->SetPoint(8,925,-0.098666);
     grae->SetPointError(8,0,0,0.2480876,0.2480876);
     grae->SetPoint(9,1025,-0.2785468);
     grae->SetPointError(9,0,0,0.3058,0.3058);
     grae->SetPoint(10,1125,-0.05397235);
     grae->SetPointError(10,0,0,0.2634776,0.2634776);
   
      for ( uint i = 0; i < nbins; ++i ) {
        double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
        errors[index][i] = grae->GetErrorYhigh(i+1);
      }
      titles.push_back( grae->GetTitle() );
      index++;
 
      grae = new TGraphAsymmErrors(nbins);
      grae->SetName("Graph");
      grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#gamma + jets)");
      grae->SetFillColor(1);
      grae->SetLineWidth(3);
      grae->SetMarkerStyle(20);
      grae->SetMarkerSize(1.5);
     grae->SetPoint(1,0,0);
     grae->SetPointError(1,0,0,0,0);
     grae->SetPoint(2,0,0);
     grae->SetPointError(2,0,0,0,0);
     grae->SetPoint(3,425,0.09278029);
     grae->SetPointError(3,0,0,0.09664966,0.09664966);
     grae->SetPoint(4,525,0.2739804);
     grae->SetPointError(4,0,0,0.1275042,0.1275042);
     grae->SetPoint(5,625,-0.1090606);
     grae->SetPointError(5,0,0,0.1778259,0.1778259);
     grae->SetPoint(6,725,0.2388484);
     grae->SetPointError(6,0,0,0.2877924,0.2877924);
     grae->SetPoint(7,825,0.1027996);
     grae->SetPointError(7,0,0,0.4100296,0.4100296);
     grae->SetPoint(8,925,-0.3708192);
     grae->SetPointError(8,0,0,0.6611959,0.6611959);
     grae->SetPoint(9,1025,-0.6806142);
     grae->SetPointError(9,0,0,1.133479,1.133479);
     grae->SetPoint(10,1125,-1);
     grae->SetPointError(10,0,0,2.386582,2.386582);
   
 
      for ( uint i = 0; i < nbins; ++i ) {
        double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
        if ( x < 275. ) { y = -1000.; }
        values[index][i] = y;
        errors[index][i] = grae->GetErrorYhigh(i+1);
      }
      titles.push_back( grae->GetTitle() );
      index++;
 
      grae = new TGraphAsymmErrors(nbins);
      grae->SetName("Graph");
      grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu#mu + jets)");
      grae->SetFillColor(1);
      grae->SetLineWidth(3);
      grae->SetMarkerStyle(20);
      grae->SetMarkerSize(1.5);
     grae->SetPoint(1,300,2.793741e-06);
     grae->SetPointError(1,0,0,0.1419687,0.1419687);
     grae->SetPoint(2,350,-0.3251202);
     grae->SetPointError(2,0,0,0.2083702,0.2083702);
     grae->SetPoint(3,425,-0.08794268);
     grae->SetPointError(3,0,0,0.213224,0.213224);
     grae->SetPoint(4,525,0.3507572);
     grae->SetPointError(4,0,0,0.2762535,0.2762535);
     grae->SetPoint(5,625,0.1011718);
     grae->SetPointError(5,0,0,0.3006639,0.3006639);
     grae->SetPoint(6,725,-0.28048);
     grae->SetPointError(6,0,0,0.4577014,0.4577014);
     grae->SetPoint(7,825,-0.6510373);
     grae->SetPointError(7,0,0,0.7113626,0.7113626);
     grae->SetPoint(8,925,-0.2464512);
     grae->SetPointError(8,0,0,0.9219985,0.9219985);
     grae->SetPoint(9,1025,0.503977);
     grae->SetPointError(9,0,0,1.695403,1.695403);
     grae->SetPoint(10,1125,0.8699347);
     grae->SetPointError(10,0,0,1.68778,1.68778);
 
      for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 22 ) {

    // 11/fb, 553p2, 2-3, mixture, 8 bins 

    energy = "8";
    lumi = "11.7";

    int nhistos = 8;
    
    syst.clear();
    syst.push_back(0.10);
    syst.push_back(0.10);
    syst.push_back(0.10);
    syst.push_back(0.20);
    syst.push_back(0.20);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

     // mu: no aT -> with aT
     grae = new TGraphAsymmErrors(nbins);
     grae->SetName("Graph");
     grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
     grae->SetFillColor(1);
     grae->SetLineWidth(3);
     grae->SetMarkerStyle(20);
     grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,-0.0548582);
    grae->SetPointError(1,0,0,0.03007567,0.03007567);
    grae->SetPoint(2,350,-0.0603661);
    grae->SetPointError(2,0,0,0.03686175,0.03686175);
    grae->SetPoint(3,425,-0.03576732);
    grae->SetPointError(3,0,0,0.0417654,0.0417654);
    grae->SetPoint(4,525,-0.07323264);
    grae->SetPointError(4,0,0,0.07437727,0.07437727);
    grae->SetPoint(5,625,-0.1400356);
    grae->SetPointError(5,0,0,0.1269372,0.1269372);
    grae->SetPoint(6,725,0.03079586);
    grae->SetPointError(6,0,0,0.225785,0.225785);
    grae->SetPoint(7,825,-0.3353326);
    grae->SetPointError(7,0,0,0.3770813,0.3770813);
    grae->SetPoint(8,925,-0.5314985);
    grae->SetPointError(8,0,0,0.4784117,0.4784117);
  
     for ( uint i = 0; i < nbins; ++i ) {
       double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
       errors[index][i] = grae->GetErrorYhigh(i+1);
     }
     titles.push_back( grae->GetTitle() );
     index++;
 
     // mu, no aT: 0->1
     grae = new TGraphAsymmErrors(nbins);
     grae->SetName("Graph");
     grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
     grae->SetFillColor(1);
     grae->SetLineWidth(3);
     grae->SetMarkerStyle(20);
     grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.07932807);
    grae->SetPointError(1,0,0,0.02792049,0.02792049);
    grae->SetPoint(2,350,0.07738148);
    grae->SetPointError(2,0,0,0.03516747,0.03516747);
    grae->SetPoint(3,425,0.04996524);
    grae->SetPointError(3,0,0,0.03683893,0.03683893);
    grae->SetPoint(4,525,0.07075566);
    grae->SetPointError(4,0,0,0.06009775,0.06009775);
    grae->SetPoint(5,625,0.02415148);
    grae->SetPointError(5,0,0,0.0937526,0.0937526);
    grae->SetPoint(6,725,-0.002330288);
    grae->SetPointError(6,0,0,0.1568307,0.1568307);
    grae->SetPoint(7,825,0.07128817);
    grae->SetPointError(7,0,0,0.2259437,0.2259437);
    grae->SetPoint(8,925,-0.09203304);
    grae->SetPointError(8,0,0,0.207548,0.207548);
  
     for ( uint i = 0; i < nbins; ++i ) {
       double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
       errors[index][i] = grae->GetErrorYhigh(i+1);
     }
     titles.push_back( grae->GetTitle() );
     index++;
 
     // mu, no aT: 1->2
     grae = new TGraphAsymmErrors(nbins);
     grae->SetName("Graph");
     grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
     grae->SetFillColor(1);
     grae->SetLineWidth(3);
     grae->SetMarkerStyle(20);
     grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,-0.06626417);
    grae->SetPointError(1,0,0,0.04162607,0.04162607);
    grae->SetPoint(2,350,-0.09220085);
    grae->SetPointError(2,0,0,0.05580856,0.05580856);
    grae->SetPoint(3,425,-0.02420883);
    grae->SetPointError(3,0,0,0.06036907,0.06036907);
    grae->SetPoint(4,525,-0.06380463);
    grae->SetPointError(4,0,0,0.1054617,0.1054617);
    grae->SetPoint(5,625,-0.2993131);
    grae->SetPointError(5,0,0,0.1840868,0.1840868);
    grae->SetPoint(6,725,-0.4615446);
    grae->SetPointError(6,0,0,0.3499595,0.3499595);
    grae->SetPoint(7,825,-0.4882799);
    grae->SetPointError(7,0,0,0.5213155,0.5213155);
    grae->SetPoint(8,925,-0.1116731);
    grae->SetPointError(8,0,0,0.5229513,0.5229513);
  
     for ( uint i = 0; i < nbins; ++i ) {
       double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
       errors[index][i] = grae->GetErrorYhigh(i+1);
     }
     titles.push_back( grae->GetTitle() );
     index++;
 
     // mu->mumu, no aT
     grae = new TGraphAsymmErrors(nbins);
     grae->SetName("Graph");
     grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
     grae->SetFillColor(1);
     grae->SetLineWidth(3);
     grae->SetMarkerStyle(20);
     grae->SetMarkerSize(1.5);
    grae->SetPoint(1,300,0.1352768);
    grae->SetPointError(1,0,0,0.041633,0.041633);
    grae->SetPoint(2,350,0.1079514);
    grae->SetPointError(2,0,0,0.05397643,0.05397643);
    grae->SetPoint(3,425,0.08140191);
    grae->SetPointError(3,0,0,0.05580835,0.05580835);
    grae->SetPoint(4,525,-0.06897097);
    grae->SetPointError(4,0,0,0.08786318,0.08786318);
    grae->SetPoint(5,625,0.1460586);
    grae->SetPointError(5,0,0,0.1289715,0.1289715);
    grae->SetPoint(6,725,0.1676168);
    grae->SetPointError(6,0,0,0.2154547,0.2154547);
    grae->SetPoint(7,825,0.1248903);
    grae->SetPointError(7,0,0,0.3000854,0.3000854);
    grae->SetPoint(8,925,0.377553);
    grae->SetPointError(8,0,0,0.2980829,0.2980829);
 
     for ( uint i = 0; i < nbins; ++i ) {
       double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
       errors[index][i] = grae->GetErrorYhigh(i+1);
     }
     titles.push_back( grae->GetTitle() );
     index++;
 
     // mumu->gamma, no aT
     grae = new TGraphAsymmErrors(nbins);
     grae->SetName("Graph");
     grae->SetTitle("#mu#mu + jets #rightarrow #gamma + jets");
     grae->SetFillColor(1);
     grae->SetLineWidth(3);
     grae->SetMarkerStyle(20);
     grae->SetMarkerSize(0.1);
    grae->SetPoint(1,0,0);
    grae->SetPointError(1,0,0,0,0);
    grae->SetPoint(2,0,0);
    grae->SetPointError(2,0,0,0,0);
    grae->SetPoint(3,425,-0.04982457);
    grae->SetPointError(3,0,0,0.05732224,0.05732224);
    grae->SetPoint(4,525,-0.190063);
    grae->SetPointError(4,0,0,0.093251,0.093251);
    grae->SetPoint(5,625,0.1463898);
    grae->SetPointError(5,0,0,0.1457974,0.1457974);
    grae->SetPoint(6,725,0.009716137);
    grae->SetPointError(6,0,0,0.2324358,0.2324358);
    grae->SetPoint(7,825,-0.07715776);
    grae->SetPointError(7,0,0,0.336785,0.336785);
    grae->SetPoint(8,925,0.5076362);
    grae->SetPointError(8,0,0,0.4158233,0.4158233);
 
     for ( uint i = 0; i < nbins; ++i ) {
       double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
       if ( x < 275. ) { y = -1000.; }
       values[index][i] = y;
       errors[index][i] = grae->GetErrorYhigh(i+1);
     }
     titles.push_back( grae->GetTitle() );
     index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.0005695891);
   grae->SetPointError(1,0,0,0.02565151,0.02565151);
   grae->SetPoint(2,350,-0.02616328);
   grae->SetPointError(2,0,0,0.03152215,0.03152215);
   grae->SetPoint(3,425,-0.08750849);
   grae->SetPointError(3,0,0,0.03225158,0.03225158);
   grae->SetPoint(4,525,-0.06502171);
   grae->SetPointError(4,0,0,0.03943689,0.03943689);
   grae->SetPoint(5,625,-0.05162088);
   grae->SetPointError(5,0,0,0.0556754,0.0556754);
   grae->SetPoint(6,725,0.2697853);
   grae->SetPointError(6,0,0,0.09060845,0.09060845);
   grae->SetPoint(7,825,0.04631919);
   grae->SetPointError(7,0,0,0.1210327,0.1210327);
   grae->SetPoint(8,925,-0.1085483);
   grae->SetPointError(8,0,0,0.1081994,0.1081994);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#gamma + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,0,0);
   grae->SetPointError(1,0,0,0,0);
   grae->SetPoint(2,0,0);
   grae->SetPointError(2,0,0,0,0);
   grae->SetPoint(3,425,0.0765702);
   grae->SetPointError(3,0,0,0.07466225,0.07466225);
   grae->SetPoint(4,525,0.1210882);
   grae->SetPointError(4,0,0,0.09441587,0.09441587);
   grae->SetPoint(5,625,0.08951763);
   grae->SetPointError(5,0,0,0.1408225,0.1408225);
   grae->SetPoint(6,725,0.2678406);
   grae->SetPointError(6,0,0,0.2231914,0.2231914);
   grae->SetPoint(7,825,0.3643913);
   grae->SetPointError(7,0,0,0.3319153,0.3319153);
   grae->SetPoint(8,925,0.5954869);
   grae->SetPointError(8,0,0,0.4613771,0.4613771);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.11258);
   grae->SetPointError(1,0,0,0.1101806,0.1101806);
   grae->SetPoint(2,350,-0.2749154);
   grae->SetPointError(2,0,0,0.1606634,0.1606634);
   grae->SetPoint(3,425,0.08022565);
   grae->SetPointError(3,0,0,0.1643356,0.1643356);
   grae->SetPoint(4,525,0.4796398);
   grae->SetPointError(4,0,0,0.2141343,0.2141343);
   grae->SetPoint(5,625,0.04870089);
   grae->SetPointError(5,0,0,0.2354536,0.2354536);
   grae->SetPoint(6,725,-0.2629874);
   grae->SetPointError(6,0,0,0.3537635,0.3537635);
   grae->SetPoint(7,825,-0.183167);
   grae->SetPointError(7,0,0,0.4867338,0.4867338);
   grae->SetPoint(8,925,0.1247106);
   grae->SetPointError(8,0,0,0.4412594,0.4412594);

    for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 23 ) {

    // 11/fb, 553p2, >=4, mixture, 8 bins

    energy = "8";
    lumi = "11.7";

    int nhistos = 8;
    
    syst.clear();
    syst.push_back(0.10);
    syst.push_back(0.10);
    syst.push_back(0.10);
    syst.push_back(0.20);
    syst.push_back(0.30);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.04661271);
   grae->SetPointError(1,0,0,0.05311414,0.05311414);
   grae->SetPoint(2,350,-0.1354332);
   grae->SetPointError(2,0,0,0.06862723,0.06862723);
   grae->SetPoint(3,425,-0.06136346);
   grae->SetPointError(3,0,0,0.07538184,0.07538184);
   grae->SetPoint(4,525,0.1012961);
   grae->SetPointError(4,0,0,0.0958097,0.0958097);
   grae->SetPoint(5,625,-0.1610576);
   grae->SetPointError(5,0,0,0.1433355,0.1433355);
   grae->SetPoint(6,725,0.2924913);
   grae->SetPointError(6,0,0,0.2395603,0.2395603);
   grae->SetPoint(7,825,0.4952695);
   grae->SetPointError(7,0,0,0.4328884,0.4328884);
   grae->SetPoint(8,925,-0.6620565);
   grae->SetPointError(8,0,0,0.541967,0.541967);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.008869397);
   grae->SetPointError(1,0,0,0.04930047,0.04930047);
   grae->SetPoint(2,350,-0.06934603);
   grae->SetPointError(2,0,0,0.06235936,0.06235936);
   grae->SetPoint(3,425,-0.06747649);
   grae->SetPointError(3,0,0,0.06556353,0.06556353);
   grae->SetPoint(4,525,0.07208692);
   grae->SetPointError(4,0,0,0.07531432,0.07531432);
   grae->SetPoint(5,625,-0.1117422);
   grae->SetPointError(5,0,0,0.101114,0.101114);
   grae->SetPoint(6,725,-0.1101752);
   grae->SetPointError(6,0,0,0.1378534,0.1378534);
   grae->SetPoint(7,825,0.1806107);
   grae->SetPointError(7,0,0,0.2180237,0.2180237);
   grae->SetPoint(8,925,-0.1954219);
   grae->SetPointError(8,0,0,0.1875318,0.1875318);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;
 
      // mu, no aT: 1->2
      grae = new TGraphAsymmErrors(nbins);
      grae->SetName("Graph");
      grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
      grae->SetFillColor(1);
      grae->SetLineWidth(3);
      grae->SetMarkerStyle(20);
      grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.09252062);
   grae->SetPointError(1,0,0,0.04821615,0.04821615);
   grae->SetPoint(2,350,0.03879051);
   grae->SetPointError(2,0,0,0.07018167,0.07018167);
   grae->SetPoint(3,425,0.02822489);
   grae->SetPointError(3,0,0,0.07150649,0.07150649);
   grae->SetPoint(4,525,-0.02417513);
   grae->SetPointError(4,0,0,0.08326054,0.08326054);
   grae->SetPoint(5,625,0.1531042);
   grae->SetPointError(5,0,0,0.1233281,0.1233281);
   grae->SetPoint(6,725,0.03005917);
   grae->SetPointError(6,0,0,0.1709384,0.1709384);
   grae->SetPoint(7,825,-0.1851216);
   grae->SetPointError(7,0,0,0.2475068,0.2475068);
   grae->SetPoint(8,925,-0.3471344);
   grae->SetPointError(8,0,0,0.2721103,0.2721103);
   
      for ( uint i = 0; i < nbins; ++i ) {
        double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
        errors[index][i] = grae->GetErrorYhigh(i+1);
      }
      titles.push_back( grae->GetTitle() );
      index++;
 
      // mu->mumu, no aT
      grae = new TGraphAsymmErrors(nbins);
      grae->SetName("Graph");
      grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
      grae->SetFillColor(1);
      grae->SetLineWidth(3);
      grae->SetMarkerStyle(20);
      grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.009092265);
   grae->SetPointError(1,0,0,0.1085735,0.1085735);
   grae->SetPoint(2,350,-0.1746513);
   grae->SetPointError(2,0,0,0.156561,0.156561);
   grae->SetPoint(3,425,0.2813321);
   grae->SetPointError(3,0,0,0.1721546,0.1721546);
   grae->SetPoint(4,525,0.4749682);
   grae->SetPointError(4,0,0,0.1946408,0.1946408);
   grae->SetPoint(5,625,0.2684642);
   grae->SetPointError(5,0,0,0.2271014,0.2271014);
   grae->SetPoint(6,725,-0.3219431);
   grae->SetPointError(6,0,0,0.2939868,0.2939868);
   grae->SetPoint(7,825,-0.1213883);
   grae->SetPointError(7,0,0,0.4120549,0.4120549);
   grae->SetPoint(8,925,0.7373891);
   grae->SetPointError(8,0,0,0.5304941,0.5304941);
 
      for ( uint i = 0; i < nbins; ++i ) {
        double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
        errors[index][i] = grae->GetErrorYhigh(i+1);
      }
      titles.push_back( grae->GetTitle() );
      index++;
 
      // mumu->gamma, no aT
      grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,-0.0458997);
   grae->SetPointError(3,0,0,0.1654557,0.1654557);
   grae->SetPoint(4,525,0.06978332);
   grae->SetPointError(4,0,0,0.1774263,0.1774263);
   grae->SetPoint(5,625,0.1042074);
   grae->SetPointError(5,0,0,0.2374312,0.2374312);
   grae->SetPoint(6,725,-0.4127053);
   grae->SetPointError(6,0,0,0.3345253,0.3345253);
   grae->SetPoint(7,825,-0.4473453);
   grae->SetPointError(7,0,0,0.4531042,0.4531042);
   grae->SetPoint(8,925,0.06191244);
   grae->SetPointError(8,0,0,0.4439505,0.4439505);
 
      for ( uint i = 0; i < nbins; ++i ) {
        double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
        if ( x < 275. ) { y = -1000.; }
        values[index][i] = y;
        errors[index][i] = grae->GetErrorYhigh(i+1);
      }
      titles.push_back( grae->GetTitle() );
      index++;
 

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.0005695891);
   grae->SetPointError(1,0,0,0.02565151,0.02565151);
   grae->SetPoint(2,350,-0.02616328);
   grae->SetPointError(2,0,0,0.03152215,0.03152215);
   grae->SetPoint(3,425,-0.08750849);
   grae->SetPointError(3,0,0,0.03225158,0.03225158);
   grae->SetPoint(4,525,-0.06502171);
   grae->SetPointError(4,0,0,0.03943689,0.03943689);
   grae->SetPoint(5,625,-0.05162088);
   grae->SetPointError(5,0,0,0.0556754,0.0556754);
   grae->SetPoint(6,725,0.2697853);
   grae->SetPointError(6,0,0,0.09060845,0.09060845);
   grae->SetPoint(7,825,0.04631919);
   grae->SetPointError(7,0,0,0.1210327,0.1210327);
   grae->SetPoint(8,925,-0.1085483);
   grae->SetPointError(8,0,0,0.1081994,0.1081994);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#gamma + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,0,0);
   grae->SetPointError(1,0,0,0,0);
   grae->SetPoint(2,0,0);
   grae->SetPointError(2,0,0,0,0);
   grae->SetPoint(3,425,0.0765702);
   grae->SetPointError(3,0,0,0.07466225,0.07466225);
   grae->SetPoint(4,525,0.1210882);
   grae->SetPointError(4,0,0,0.09441587,0.09441587);
   grae->SetPoint(5,625,0.08951763);
   grae->SetPointError(5,0,0,0.1408225,0.1408225);
   grae->SetPoint(6,725,0.2678406);
   grae->SetPointError(6,0,0,0.2231914,0.2231914);
   grae->SetPoint(7,825,0.3643913);
   grae->SetPointError(7,0,0,0.3319153,0.3319153);
   grae->SetPoint(8,925,0.5954869);
   grae->SetPointError(8,0,0,0.4613771,0.4613771);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.11258);
   grae->SetPointError(1,0,0,0.1101806,0.1101806);
   grae->SetPoint(2,350,-0.2749154);
   grae->SetPointError(2,0,0,0.1606634,0.1606634);
   grae->SetPoint(3,425,0.08022565);
   grae->SetPointError(3,0,0,0.1643356,0.1643356);
   grae->SetPoint(4,525,0.4796398);
   grae->SetPointError(4,0,0,0.2141343,0.2141343);
   grae->SetPoint(5,625,0.04870089);
   grae->SetPointError(5,0,0,0.2354536,0.2354536);
   grae->SetPoint(6,725,-0.2629874);
   grae->SetPointError(6,0,0,0.3537635,0.3537635);
   grae->SetPoint(7,825,-0.183167);
   grae->SetPointError(7,0,0,0.4867338,0.4867338);
   grae->SetPoint(8,925,0.1247106);
   grae->SetPointError(8,0,0,0.4412594,0.4412594);
 
      for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 42 ) {

    // 18/fb, 553p2, 2-3, mixture, 8 bins 

    energy = "8";
    lumi = "18.7";
    njet = "le3";

    int nhistos = 8;
    
    syst.clear();
    syst.push_back(0.10);
    syst.push_back(0.10);
    syst.push_back(0.10);
    syst.push_back(0.10);
    syst.push_back(0.20);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

     // mu: no aT -> with aT
     grae = new TGraphAsymmErrors(nbins);
     grae->SetName("Graph");
     grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
     grae->SetFillColor(1);
     grae->SetLineWidth(3);
     grae->SetMarkerStyle(20);
     grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.02830745);
   grae->SetPointError(1,0,0,0.03040306,0.03040306);
   grae->SetPoint(2,350,-0.04842127);
   grae->SetPointError(2,0,0,0.03088054,0.03088054);
   grae->SetPoint(3,425,-0.03626645);
   grae->SetPointError(3,0,0,0.03443865,0.03443865);
   grae->SetPoint(4,525,-0.04149093);
   grae->SetPointError(4,0,0,0.06098234,0.06098234);
   grae->SetPoint(5,625,-0.1590971);
   grae->SetPointError(5,0,0,0.1041572,0.1041572);
   grae->SetPoint(6,725,-0.002870784);
   grae->SetPointError(6,0,0,0.1773236,0.1773236);
   grae->SetPoint(7,825,-0.2292653);
   grae->SetPointError(7,0,0,0.2982348,0.2982348);
   grae->SetPoint(8,925,-0.4067336);
   grae->SetPointError(8,0,0,0.3713284,0.3713284);
  
     for ( uint i = 0; i < nbins; ++i ) {
       double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
       errors[index][i] = grae->GetErrorYhigh(i+1);
     }
     titles.push_back( grae->GetTitle() );
     index++;
 
     // mu, no aT: 0->1
     grae = new TGraphAsymmErrors(nbins);
     grae->SetName("Graph");
     grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
     grae->SetFillColor(1);
     grae->SetLineWidth(3);
     grae->SetMarkerStyle(20);
     grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1441884);
   grae->SetPointError(1,0,0,0.02600097,0.02600097);
   grae->SetPoint(2,350,0.111872);
   grae->SetPointError(2,0,0,0.02939634,0.02939634);
   grae->SetPoint(3,425,0.07578086);
   grae->SetPointError(3,0,0,0.03036036,0.03036036);
   grae->SetPoint(4,525,0.06823984);
   grae->SetPointError(4,0,0,0.04907832,0.04907832);
   grae->SetPoint(5,625,0.05555008);
   grae->SetPointError(5,0,0,0.07595524,0.07595524);
   grae->SetPoint(6,725,0.01552418);
   grae->SetPointError(6,0,0,0.1232152,0.1232152);
   grae->SetPoint(7,825,0.007762816);
   grae->SetPointError(7,0,0,0.1857697,0.1857697);
   grae->SetPoint(8,925,-0.06301305);
   grae->SetPointError(8,0,0,0.1713697,0.1713697);
  
     for ( uint i = 0; i < nbins; ++i ) {
       double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
       errors[index][i] = grae->GetErrorYhigh(i+1);
     }
     titles.push_back( grae->GetTitle() );
     index++;
 
     // mu, no aT: 1->2
     grae = new TGraphAsymmErrors(nbins);
     grae->SetName("Graph");
     grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
     grae->SetFillColor(1);
     grae->SetLineWidth(3);
     grae->SetMarkerStyle(20);
     grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.05582098);
   grae->SetPointError(1,0,0,0.03425992,0.03425992);
   grae->SetPoint(2,350,-0.05473331);
   grae->SetPointError(2,0,0,0.04497547,0.04497547);
   grae->SetPoint(3,425,-0.007118279);
   grae->SetPointError(3,0,0,0.04847374,0.04847374);
   grae->SetPoint(4,525,-0.04600748);
   grae->SetPointError(4,0,0,0.08566305,0.08566305);
   grae->SetPoint(5,625,-0.2428315);
   grae->SetPointError(5,0,0,0.1469793,0.1469793);
   grae->SetPoint(6,725,-0.3994991);
   grae->SetPointError(6,0,0,0.2666527,0.2666527);
   grae->SetPoint(7,825,-0.354267);
   grae->SetPointError(7,0,0,0.4304296,0.4304296);
   grae->SetPoint(8,925,-0.05473629);
   grae->SetPointError(8,0,0,0.4099286,0.4099286);
  
     for ( uint i = 0; i < nbins; ++i ) {
       double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
       errors[index][i] = grae->GetErrorYhigh(i+1);
     }
     titles.push_back( grae->GetTitle() );
     index++;
 
     // mu->mumu, no aT
     grae = new TGraphAsymmErrors(nbins);
     grae->SetName("Graph");
     grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
     grae->SetFillColor(1);
     grae->SetLineWidth(3);
     grae->SetMarkerStyle(20);
     grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1021524);
   grae->SetPointError(1,0,0,0.0283057,0.0283057);
   grae->SetPoint(2,350,0.1173886);
   grae->SetPointError(2,0,0,0.03462356,0.03462356);
   grae->SetPoint(3,425,0.09098752);
   grae->SetPointError(3,0,0,0.03522444,0.03522444);
   grae->SetPoint(4,525,-0.04424214);
   grae->SetPointError(4,0,0,0.05316926,0.05316926);
   grae->SetPoint(5,625,0.1580795);
   grae->SetPointError(5,0,0,0.08061528,0.08061528);
   grae->SetPoint(6,725,0.09982629);
   grae->SetPointError(6,0,0,0.1231553,0.1231553);
   grae->SetPoint(7,825,0.1141393);
   grae->SetPointError(7,0,0,0.1763231,0.1763231);
   grae->SetPoint(8,925,-0.0312421);
   grae->SetPointError(8,0,0,0.1537213,0.1537213);
 
     for ( uint i = 0; i < nbins; ++i ) {
       double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
       errors[index][i] = grae->GetErrorYhigh(i+1);
     }
     titles.push_back( grae->GetTitle() );
     index++;
 
     // mumu->gamma, no aT
     grae = new TGraphAsymmErrors(nbins);
     grae->SetName("Graph");
     grae->SetTitle("#mu#mu + jets #rightarrow #gamma + jets");
     grae->SetFillColor(1);
     grae->SetLineWidth(3);
     grae->SetMarkerStyle(20);
     grae->SetMarkerSize(0.1);
   grae->SetPoint(1,0,0);
   grae->SetPointError(1,0,0,0,0);
   grae->SetPoint(2,0,0);
   grae->SetPointError(2,0,0,0,0);
   grae->SetPoint(3,425,-0.06521067);
   grae->SetPointError(3,0,0,0.03795643,0.03795643);
   grae->SetPoint(4,525,-0.2052821);
   grae->SetPointError(4,0,0,0.0604652,0.0604652);
   grae->SetPoint(5,625,0.01829117);
   grae->SetPointError(5,0,0,0.09725351,0.09725351);
   grae->SetPoint(6,725,-0.07583478);
   grae->SetPointError(6,0,0,0.1530836,0.1530836);
   grae->SetPoint(7,825,0.1457793);
   grae->SetPointError(7,0,0,0.2484594,0.2484594);
   grae->SetPoint(8,925,-0.06185581);
   grae->SetPointError(8,0,0,0.2428195,0.2428195);
 
     for ( uint i = 0; i < nbins; ++i ) {
       double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
       if ( x < 275. ) { y = -1000.; }
       values[index][i] = y;
       errors[index][i] = grae->GetErrorYhigh(i+1);
     }
     titles.push_back( grae->GetTitle() );
     index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.00389501);
   grae->SetPointError(1,0,0,0.02694342,0.02694342);
   grae->SetPoint(2,350,-0.01971203);
   grae->SetPointError(2,0,0,0.02587963,0.02587963);
   grae->SetPoint(3,425,-0.07007172);
   grae->SetPointError(3,0,0,0.02682278,0.02682278);
   grae->SetPoint(4,525,-0.02113186);
   grae->SetPointError(4,0,0,0.03236048,0.03236048);
   grae->SetPoint(5,625,-0.03637465);
   grae->SetPointError(5,0,0,0.04439642,0.04439642);
   grae->SetPoint(6,725,0.1353426);
   grae->SetPointError(6,0,0,0.06781349,0.06781349);
   grae->SetPoint(7,825,-0.01657415);
   grae->SetPointError(7,0,0,0.09644335,0.09644335);
   grae->SetPoint(8,925,-0.1097893);
   grae->SetPointError(8,0,0,0.08621139,0.08621139);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#gamma + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,0,0);
   grae->SetPointError(1,0,0,0,0);
   grae->SetPoint(2,0,0);
   grae->SetPointError(2,0,0,0,0);
   grae->SetPoint(3,425,0.04467433);
   grae->SetPointError(3,0,0,0.0676834,0.0676834);
   grae->SetPoint(4,525,0.09998912);
   grae->SetPointError(4,0,0,0.08713809,0.08713809);
   grae->SetPoint(5,625,0.03466702);
   grae->SetPointError(5,0,0,0.123327,0.123327);
   grae->SetPoint(6,725,0.2994604);
   grae->SetPointError(6,0,0,0.1965639,0.1965639);
   grae->SetPoint(7,825,0.3119324);
   grae->SetPointError(7,0,0,0.3058743,0.3058743);
   grae->SetPoint(8,925,0.4550672);
   grae->SetPointError(8,0,0,0.3492483,0.3492483);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.1048325);
   grae->SetPointError(1,0,0,0.06775633,0.06775633);
   grae->SetPoint(2,350,-0.1713302);
   grae->SetPointError(2,0,0,0.09955166,0.09955166);
   grae->SetPoint(3,425,0.0326672);
   grae->SetPointError(3,0,0,0.09449645,0.09449645);
   grae->SetPoint(4,525,0.3420602);
   grae->SetPointError(4,0,0,0.1259856,0.1259856);
   grae->SetPoint(5,625,-0.008494672);
   grae->SetPointError(5,0,0,0.1445688,0.1445688);
   grae->SetPoint(6,725,-0.190231);
   grae->SetPointError(6,0,0,0.2207145,0.2207145);
   grae->SetPoint(7,825,0.1533652);
   grae->SetPointError(7,0,0,0.309377,0.309377);
   grae->SetPoint(8,925,0.006665431);
   grae->SetPointError(8,0,0,0.269327,0.269327);
 
      for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 43 ) {

    // 18/fb, 553p2, >=4, mixture, 8 bins

    energy = "8";
    lumi = "18.7";
    njet = "ge4";

    int nhistos = 8;
    
    syst.clear();
    syst.push_back(0.10);
    syst.push_back(0.10);
    syst.push_back(0.10);
    syst.push_back(0.20);
    syst.push_back(0.20);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.0166241);
   grae->SetPointError(1,0,0,0.05123598,0.05123598);
   grae->SetPoint(2,350,-0.1122772);
   grae->SetPointError(2,0,0,0.05313585,0.05313585);
   grae->SetPoint(3,425,-0.1077604);
   grae->SetPointError(3,0,0,0.05827368,0.05827368);
   grae->SetPoint(4,525,-0.003629828);
   grae->SetPointError(4,0,0,0.07319543,0.07319543);
   grae->SetPoint(5,625,-0.1385773);
   grae->SetPointError(5,0,0,0.1098369,0.1098369);
   grae->SetPoint(6,725,0.08996313);
   grae->SetPointError(6,0,0,0.1811404,0.1811404);
   grae->SetPoint(7,825,0.3519611);
   grae->SetPointError(7,0,0,0.332177,0.332177);
   grae->SetPoint(8,925,-0.5810519);
   grae->SetPointError(8,0,0,0.398623,0.398623);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.09726359);
   grae->SetPointError(1,0,0,0.05344529,0.05344529);
   grae->SetPoint(2,350,-0.04665348);
   grae->SetPointError(2,0,0,0.0488755,0.0488755);
   grae->SetPoint(3,425,-0.04243294);
   grae->SetPointError(3,0,0,0.05284108,0.05284108);
   grae->SetPoint(4,525,0.03985102);
   grae->SetPointError(4,0,0,0.05719635,0.05719635);
   grae->SetPoint(5,625,-0.1080215);
   grae->SetPointError(5,0,0,0.07571437,0.07571437);
   grae->SetPoint(6,725,-0.1505843);
   grae->SetPointError(6,0,0,0.1059535,0.1059535);
   grae->SetPoint(7,825,-0.08698489);
   grae->SetPointError(7,0,0,0.161541,0.161541);
   grae->SetPoint(8,925,-0.2434888);
   grae->SetPointError(8,0,0,0.1443955,0.1443955);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;
 
      // mu, no aT: 1->2
      grae = new TGraphAsymmErrors(nbins);
      grae->SetName("Graph");
      grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
      grae->SetFillColor(1);
      grae->SetLineWidth(3);
      grae->SetMarkerStyle(20);
      grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.01225817);
   grae->SetPointError(1,0,0,0.03787374,0.03787374);
   grae->SetPoint(2,350,0.08875529);
   grae->SetPointError(2,0,0,0.05406094,0.05406094);
   grae->SetPoint(3,425,0.0645663);
   grae->SetPointError(3,0,0,0.05484707,0.05484707);
   grae->SetPoint(4,525,0.07742724);
   grae->SetPointError(4,0,0,0.0645151,0.0645151);
   grae->SetPoint(5,625,0.05368706);
   grae->SetPointError(5,0,0,0.0918347,0.0918347);
   grae->SetPoint(6,725,0.199296);
   grae->SetPointError(6,0,0,0.1392207,0.1392207);
   grae->SetPoint(7,825,0.05237678);
   grae->SetPointError(7,0,0,0.2088297,0.2088297);
   grae->SetPoint(8,925,-0.09543222);
   grae->SetPointError(8,0,0,0.2088065,0.2088065);
   
      for ( uint i = 0; i < nbins; ++i ) {
        double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
        errors[index][i] = grae->GetErrorYhigh(i+1);
      }
      titles.push_back( grae->GetTitle() );
      index++;
 
      // mu->mumu, no aT
      grae = new TGraphAsymmErrors(nbins);
      grae->SetName("Graph");
      grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
      grae->SetFillColor(1);
      grae->SetLineWidth(3);
      grae->SetMarkerStyle(20);
      grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.009531155);
   grae->SetPointError(1,0,0,0.06982161,0.06982161);
   grae->SetPoint(2,350,-0.05543449);
   grae->SetPointError(2,0,0,0.1013211,0.1013211);
   grae->SetPoint(3,425,0.2115203);
   grae->SetPointError(3,0,0,0.101195,0.101195);
   grae->SetPoint(4,525,0.3103753);
   grae->SetPointError(4,0,0,0.1150605,0.1150605);
   grae->SetPoint(5,625,0.1915855);
   grae->SetPointError(5,0,0,0.1454524,0.1454524);
   grae->SetPoint(6,725,-0.2155625);
   grae->SetPointError(6,0,0,0.1950038,0.1950038);
   grae->SetPoint(7,825,0.3066664);
   grae->SetPointError(7,0,0,0.3029229,0.3029229);
   grae->SetPoint(8,925,0.0954879);
   grae->SetPointError(8,0,0,0.2471279,0.2471279);
 
      for ( uint i = 0; i < nbins; ++i ) {
        double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
        errors[index][i] = grae->GetErrorYhigh(i+1);
      }
      titles.push_back( grae->GetTitle() );
      index++;
 
      // mumu->gamma, no aT
      grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,-0.07595482);
   grae->SetPointError(3,0,0,0.1058347,0.1058347);
   grae->SetPoint(4,525,-0.03039107);
   grae->SetPointError(4,0,0,0.1198184,0.1198184);
   grae->SetPoint(5,625,-0.0241874);
   grae->SetPointError(5,0,0,0.1621843,0.1621843);
   grae->SetPoint(6,725,-0.4240992);
   grae->SetPointError(6,0,0,0.2319989,0.2319989);
   grae->SetPoint(7,825,0.007294223);
   grae->SetPointError(7,0,0,0.3245537,0.3245537);
   grae->SetPoint(8,925,-0.3509596);
   grae->SetPointError(8,0,0,0.3066588,0.3066588);
 
      for ( uint i = 0; i < nbins; ++i ) {
        double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
        if ( x < 275. ) { y = -1000.; }
        values[index][i] = y;
        errors[index][i] = grae->GetErrorYhigh(i+1);
      }
      titles.push_back( grae->GetTitle() );
      index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.00389501);
   grae->SetPointError(1,0,0,0.02694342,0.02694342);
   grae->SetPoint(2,350,-0.01971203);
   grae->SetPointError(2,0,0,0.02587963,0.02587963);
   grae->SetPoint(3,425,-0.07007172);
   grae->SetPointError(3,0,0,0.02682278,0.02682278);
   grae->SetPoint(4,525,-0.02113186);
   grae->SetPointError(4,0,0,0.03236048,0.03236048);
   grae->SetPoint(5,625,-0.03637465);
   grae->SetPointError(5,0,0,0.04439642,0.04439642);
   grae->SetPoint(6,725,0.1353426);
   grae->SetPointError(6,0,0,0.06781349,0.06781349);
   grae->SetPoint(7,825,-0.01657415);
   grae->SetPointError(7,0,0,0.09644335,0.09644335);
   grae->SetPoint(8,925,-0.1097893);
   grae->SetPointError(8,0,0,0.08621139,0.08621139);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#gamma + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,0,0);
   grae->SetPointError(1,0,0,0,0);
   grae->SetPoint(2,0,0);
   grae->SetPointError(2,0,0,0,0);
   grae->SetPoint(3,425,0.04467433);
   grae->SetPointError(3,0,0,0.0676834,0.0676834);
   grae->SetPoint(4,525,0.09998912);
   grae->SetPointError(4,0,0,0.08713809,0.08713809);
   grae->SetPoint(5,625,0.03466702);
   grae->SetPointError(5,0,0,0.123327,0.123327);
   grae->SetPoint(6,725,0.2994604);
   grae->SetPointError(6,0,0,0.1965639,0.1965639);
   grae->SetPoint(7,825,0.3119324);
   grae->SetPointError(7,0,0,0.3058743,0.3058743);
   grae->SetPoint(8,925,0.4550672);
   grae->SetPointError(8,0,0,0.3492483,0.3492483);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.1048325);
   grae->SetPointError(1,0,0,0.06775633,0.06775633);
   grae->SetPoint(2,350,-0.1713302);
   grae->SetPointError(2,0,0,0.09955166,0.09955166);
   grae->SetPoint(3,425,0.0326672);
   grae->SetPointError(3,0,0,0.09449645,0.09449645);
   grae->SetPoint(4,525,0.3420602);
   grae->SetPointError(4,0,0,0.1259856,0.1259856);
   grae->SetPoint(5,625,-0.008494672);
   grae->SetPointError(5,0,0,0.1445688,0.1445688);
   grae->SetPoint(6,725,-0.190231);
   grae->SetPointError(6,0,0,0.2207145,0.2207145);
   grae->SetPoint(7,825,0.1533652);
   grae->SetPointError(7,0,0,0.309377,0.309377);
   grae->SetPoint(8,925,0.006665431);
   grae->SetPointError(8,0,0,0.269327,0.269327);
 
      for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 230 ) {

    // 11/fb, 553p2, >=4, mixture, 8 bins

    energy = "8";
    lumi = "11.7";

    int nhistos = 1;
    
    syst.clear();
    syst.push_back(0.10);
    syst.push_back(0.10);
    syst.push_back(0.10);
    syst.push_back(0.20);
    syst.push_back(0.30);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.0005695891);
   grae->SetPointError(1,0,0,0.02565151,0.02565151);
   grae->SetPoint(2,350,-0.02616328);
   grae->SetPointError(2,0,0,0.03152215,0.03152215);
   grae->SetPoint(3,425,-0.08750849);
   grae->SetPointError(3,0,0,0.03225158,0.03225158);
   grae->SetPoint(4,525,-0.06502171);
   grae->SetPointError(4,0,0,0.03943689,0.03943689);
   grae->SetPoint(5,625,-0.05162088);
   grae->SetPointError(5,0,0,0.0556754,0.0556754);
   grae->SetPoint(6,725,0.2697853);
   grae->SetPointError(6,0,0,0.09060845,0.09060845);
   grae->SetPoint(7,825,0.04631919);
   grae->SetPointError(7,0,0,0.1210327,0.1210327);
   grae->SetPoint(8,925,-0.1085483);
   grae->SetPointError(8,0,0,0.1081994,0.1081994);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    return nhistos;

    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
  } else if ( set_option == 24 ) {

    // 11/fb, 553p2, 2-3, mixture, 8 bins 

    energy = "8";
    lumi = "11.7";

    int nhistos = 8;
    
    syst.clear();
    syst.push_back(0.10);
    syst.push_back(0.10);
    syst.push_back(0.10);
    syst.push_back(0.20);
    syst.push_back(0.20);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.03693739);
   grae->SetPointError(1,0,0,0.03112361,0.03112361);
   grae->SetPoint(2,350,-0.04563713);
   grae->SetPointError(2,0,0,0.03727786,0.03727786);
   grae->SetPoint(3,425,-0.02657793);
   grae->SetPointError(3,0,0,0.04207249,0.04207249);
   grae->SetPoint(4,525,-0.07445496);
   grae->SetPointError(4,0,0,0.07449139,0.07449139);
   grae->SetPoint(5,625,-0.1474447);
   grae->SetPointError(5,0,0,0.1268713,0.1268713);
   grae->SetPoint(6,725,0.003398276);
   grae->SetPointError(6,0,0,0.2233636,0.2233636);
   grae->SetPoint(7,825,-0.3510622);
   grae->SetPointError(7,0,0,0.375176,0.375176);
   grae->SetPoint(8,925,-0.5363633);
   grae->SetPointError(8,0,0,0.4769889,0.4769889);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.3522908);
   grae->SetPointError(1,0,0,0.03136335,0.03136335);
   grae->SetPoint(2,350,0.345388);
   grae->SetPointError(2,0,0,0.03965973,0.03965973);
   grae->SetPoint(3,425,0.3042249);
   grae->SetPointError(3,0,0,0.0411169,0.0411169);
   grae->SetPoint(4,525,0.3062347);
   grae->SetPointError(4,0,0,0.06708324,0.06708324);
   grae->SetPoint(5,625,0.225768);
   grae->SetPointError(5,0,0,0.1023557,0.1023557);
   grae->SetPoint(6,725,0.1767982);
   grae->SetPointError(6,0,0,0.1692398,0.1692398);
   grae->SetPoint(7,825,0.2457454);
   grae->SetPointError(7,0,0,0.2460841,0.2460841);
   grae->SetPoint(8,925,0.04219919);
   grae->SetPointError(8,0,0,0.2175458,0.2175458);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 1->2
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.05613523);
   grae->SetPointError(1,0,0,0.04329188,0.04329188);
   grae->SetPoint(2,350,0.02845411);
   grae->SetPointError(2,0,0,0.05790894,0.05790894);
   grae->SetPoint(3,425,0.1107588);
   grae->SetPointError(3,0,0,0.06337875,0.06337875);
   grae->SetPoint(4,525,0.08321115);
   grae->SetPointError(4,0,0,0.1111175,0.1111175);
   grae->SetPoint(5,625,-0.181872);
   grae->SetPointError(5,0,0,0.1897713,0.1897713);
   grae->SetPoint(6,725,-0.3641274);
   grae->SetPointError(6,0,0,0.3603639,0.3603639);
   grae->SetPoint(7,825,-0.3952155);
   grae->SetPointError(7,0,0,0.5386013,0.5386013);
   grae->SetPoint(8,925,0.06495327);
   grae->SetPointError(8,0,0,0.5614005,0.5614005);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.06049229);
   grae->SetPointError(1,0,0,0.04117406,0.04117406);
   grae->SetPoint(2,350,0.03063481);
   grae->SetPointError(2,0,0,0.05326427,0.05326427);
   grae->SetPoint(3,425,0.009478347);
   grae->SetPointError(3,0,0,0.05514037,0.05514037);
   grae->SetPoint(4,525,-0.1176822);
   grae->SetPointError(4,0,0,0.08773292,0.08773292);
   grae->SetPoint(5,625,0.09230731);
   grae->SetPointError(5,0,0,0.1272211,0.1272211);
   grae->SetPoint(6,725,0.117822);
   grae->SetPointError(6,0,0,0.2121976,0.2121976);
   grae->SetPoint(7,825,0.08806207);
   grae->SetPointError(7,0,0,0.2973869,0.2973869);
   grae->SetPoint(8,925,0.3460223);
   grae->SetPointError(8,0,0,0.2945599,0.2945599);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mumu->gamma, no aT
    grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,-0.0331756);
   grae->SetPointError(3,0,0,0.05783626,0.05783626);
   grae->SetPoint(4,525,-0.1769477);
   grae->SetPointError(4,0,0,0.09387518,0.09387518);
   grae->SetPoint(5,625,0.1596395);
   grae->SetPointError(5,0,0,0.1469274,0.1469274);
   grae->SetPoint(6,725,0.01674581);
   grae->SetPointError(6,0,0,0.2332761,0.2332761);
   grae->SetPoint(7,825,-0.0674598);
   grae->SetPointError(7,0,0,0.3381702,0.3381702);
   grae->SetPoint(8,925,0.5286271);
   grae->SetPointError(8,0,0,0.4215346,0.4215346);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1985665);
   grae->SetPointError(1,0,0,0.02836791,0.02836791);
   grae->SetPoint(2,350,0.1701261);
   grae->SetPointError(2,0,0,0.03395709,0.03395709);
   grae->SetPoint(3,425,0.102229);
   grae->SetPointError(3,0,0,0.03471205,0.03471205);
   grae->SetPoint(4,525,0.1366804);
   grae->SetPointError(4,0,0,0.04186897,0.04186897);
   grae->SetPoint(5,625,0.1545897);
   grae->SetPointError(5,0,0,0.05923602,0.05923602);
   grae->SetPoint(6,725,0.549097);
   grae->SetPointError(6,0,0,0.1050738,0.1050738);
   grae->SetPoint(7,825,0.2699094);
   grae->SetPointError(7,0,0,0.1316974,0.1316974);
   grae->SetPoint(8,925,0.06293647);
   grae->SetPointError(8,0,0,0.1127411,0.1127411);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#gamma + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,0,0);
   grae->SetPointError(1,0,0,0,0);
   grae->SetPoint(2,0,0);
   grae->SetPointError(2,0,0,0,0);
   grae->SetPoint(3,425,0.0765702);
   grae->SetPointError(3,0,0,0.07466225,0.07466225);
   grae->SetPoint(4,525,0.1210882);
   grae->SetPointError(4,0,0,0.09441587,0.09441587);
   grae->SetPoint(5,625,0.08951763);
   grae->SetPointError(5,0,0,0.1408225,0.1408225);
   grae->SetPoint(6,725,0.2678406);
   grae->SetPointError(6,0,0,0.2231914,0.2231914);
   grae->SetPoint(7,825,0.3643913);
   grae->SetPointError(7,0,0,0.3319153,0.3319153);
   grae->SetPoint(8,925,0.5954869);
   grae->SetPointError(8,0,0,0.4613771,0.4613771);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.08146677);
   grae->SetPointError(1,0,0,0.1133193,0.1133193);
   grae->SetPoint(2,350,-0.2576047);
   grae->SetPointError(2,0,0,0.1638928,0.1638928);
   grae->SetPoint(3,425,0.1362042);
   grae->SetPointError(3,0,0,0.1715621,0.1715621);
   grae->SetPoint(4,525,0.5536013);
   grae->SetPointError(4,0,0,0.2274302,0.2274302);
   grae->SetPoint(5,625,0.09573878);
   grae->SetPointError(5,0,0,0.2428634,0.2428634);
   grae->SetPoint(6,725,-0.2297499);
   grae->SetPointError(6,0,0,0.3585898,0.3585898);
   grae->SetPoint(7,825,-0.1605828);
   grae->SetPointError(7,0,0,0.4923044,0.4923044);
   grae->SetPoint(8,925,0.1883984);
   grae->SetPointError(8,0,0,0.4587862,0.4587862);

    for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 25 ) {

    // 11/fb, 553p2, >=4, mixture, 8 bins

    energy = "8";
    lumi = "11.7";

    int nhistos = 8;
    
    syst.clear();
    syst.push_back(0.10);
    syst.push_back(0.10);
    syst.push_back(0.10);
    syst.push_back(0.20);
    syst.push_back(0.30);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.05773176);
   grae->SetPointError(1,0,0,0.05854098,0.05854098);
   grae->SetPoint(2,350,-0.1403123);
   grae->SetPointError(2,0,0,0.0689494,0.0689494);
   grae->SetPoint(3,425,-0.07014288);
   grae->SetPointError(3,0,0,0.07546122,0.07546122);
   grae->SetPoint(4,525,0.07831582);
   grae->SetPointError(4,0,0,0.0946931,0.0946931);
   grae->SetPoint(5,625,-0.1871356);
   grae->SetPointError(5,0,0,0.1415266,0.1415266);
   grae->SetPoint(6,725,0.2393693);
   grae->SetPointError(6,0,0,0.2307421,0.2307421);
   grae->SetPoint(7,825,0.3697157);
   grae->SetPointError(7,0,0,0.3963437,0.3963437);
   grae->SetPoint(8,925,-0.6710476);
   grae->SetPointError(8,0,0,0.5340043,0.5340043);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu, no aT: 0->1
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.2766688);
   grae->SetPointError(1,0,0,0.05474697,0.05474697);
   grae->SetPoint(2,350,0.1722806);
   grae->SetPointError(2,0,0,0.06656115,0.06656115);
   grae->SetPoint(3,425,0.1734265);
   grae->SetPointError(3,0,0,0.07043922,0.07043922);
   grae->SetPoint(4,525,0.3508758);
   grae->SetPointError(4,0,0,0.08387267,0.08387267);
   grae->SetPoint(5,625,0.1223192);
   grae->SetPointError(5,0,0,0.106697,0.106697);
   grae->SetPoint(6,725,0.1223089);
   grae->SetPointError(6,0,0,0.1451113,0.1451113);
   grae->SetPoint(7,825,0.4943273);
   grae->SetPointError(7,0,0,0.2529994,0.2529994);
   grae->SetPoint(8,925,0.0155342);
   grae->SetPointError(8,0,0,0.194547,0.194547);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;
 
      // mu, no aT: 1->2
      grae = new TGraphAsymmErrors(nbins);
      grae->SetName("Graph");
      grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
      grae->SetFillColor(1);
      grae->SetLineWidth(3);
      grae->SetMarkerStyle(20);
      grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.05517262);
   grae->SetPointError(1,0,0,0.04851318,0.04851318);
   grae->SetPoint(2,350,0.08173446);
   grae->SetPointError(2,0,0,0.07096752,0.07096752);
   grae->SetPoint(3,425,0.07113314);
   grae->SetPointError(3,0,0,0.07223998,0.07223998);
   grae->SetPoint(4,525,0.01944286);
   grae->SetPointError(4,0,0,0.08392015,0.08392015);
   grae->SetPoint(5,625,0.2113582);
   grae->SetPointError(5,0,0,0.1257937,0.1257937);
   grae->SetPoint(6,725,0.08759045);
   grae->SetPointError(6,0,0,0.1731041,0.1731041);
   grae->SetPoint(7,825,-0.1400369);
   grae->SetPointError(7,0,0,0.2481987,0.2481987);
   grae->SetPoint(8,925,-0.301244);
   grae->SetPointError(8,0,0,0.2726859,0.2726859);
   
      for ( uint i = 0; i < nbins; ++i ) {
        double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
        errors[index][i] = grae->GetErrorYhigh(i+1);
      }
      titles.push_back( grae->GetTitle() );
      index++;
 
      // mu->mumu, no aT
      grae = new TGraphAsymmErrors(nbins);
      grae->SetName("Graph");
      grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
      grae->SetFillColor(1);
      grae->SetLineWidth(3);
      grae->SetMarkerStyle(20);
      grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.1883982);
   grae->SetPointError(1,0,0,0.1074085,0.1074085);
   grae->SetPoint(2,350,-0.3464716);
   grae->SetPointError(2,0,0,0.1570951,0.1570951);
   grae->SetPoint(3,425,0.03935144);
   grae->SetPointError(3,0,0,0.160514,0.160514);
   grae->SetPoint(4,525,0.2041592);
   grae->SetPointError(4,0,0,0.172688,0.172688);
   grae->SetPoint(5,625,0.03535927);
   grae->SetPointError(5,0,0,0.2118097,0.2118097);
   grae->SetPoint(6,725,-0.4444772);
   grae->SetPointError(6,0,0,0.2943726,0.2943726);
   grae->SetPoint(7,825,-0.2812203);
   grae->SetPointError(7,0,0,0.4062696,0.4062696);
   grae->SetPoint(8,925,0.5058076);
   grae->SetPointError(8,0,0,0.4643104,0.4643104);
 
      for ( uint i = 0; i < nbins; ++i ) {
        double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
        errors[index][i] = grae->GetErrorYhigh(i+1);
      }
      titles.push_back( grae->GetTitle() );
      index++;
 
      // mumu->gamma, no aT
      grae = new TGraphAsymmErrors(nbins);
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
   grae->SetPoint(3,425,0.01912165);
   grae->SetPointError(3,0,0,0.1712809,0.1712809);
   grae->SetPoint(4,525,0.139158);
   grae->SetPointError(4,0,0,0.1843831,0.1843831);
   grae->SetPoint(5,625,0.1649464);
   grae->SetPointError(5,0,0,0.2458821,0.2458821);
   grae->SetPoint(6,725,-0.3828802);
   grae->SetPointError(6,0,0,0.3383693,0.3383693);
   grae->SetPoint(7,825,-0.4265511);
   grae->SetPointError(7,0,0,0.4571255,0.4571255);
   grae->SetPoint(8,925,0.1399153);
   grae->SetPointError(8,0,0,0.4603481,0.4603481);
 
      for ( uint i = 0; i < nbins; ++i ) {
        double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
        if ( x < 275. ) { y = -1000.; }
        values[index][i] = y;
        errors[index][i] = grae->GetErrorYhigh(i+1);
      }
      titles.push_back( grae->GetTitle() );
      index++;
 

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,0.1985665);
   grae->SetPointError(1,0,0,0.02836791,0.02836791);
   grae->SetPoint(2,350,0.1701261);
   grae->SetPointError(2,0,0,0.03395709,0.03395709);
   grae->SetPoint(3,425,0.102229);
   grae->SetPointError(3,0,0,0.03471205,0.03471205);
   grae->SetPoint(4,525,0.1366804);
   grae->SetPointError(4,0,0,0.04186897,0.04186897);
   grae->SetPoint(5,625,0.1545897);
   grae->SetPointError(5,0,0,0.05923602,0.05923602);
   grae->SetPoint(6,725,0.549097);
   grae->SetPointError(6,0,0,0.1050738,0.1050738);
   grae->SetPoint(7,825,0.2699094);
   grae->SetPointError(7,0,0,0.1316974,0.1316974);
   grae->SetPoint(8,925,0.06293647);
   grae->SetPointError(8,0,0,0.1127411,0.1127411);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#gamma + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,0,0);
   grae->SetPointError(1,0,0,0,0);
   grae->SetPoint(2,0,0);
   grae->SetPointError(2,0,0,0,0);
   grae->SetPoint(3,425,0.0765702);
   grae->SetPointError(3,0,0,0.07466225,0.07466225);
   grae->SetPoint(4,525,0.1210882);
   grae->SetPointError(4,0,0,0.09441587,0.09441587);
   grae->SetPoint(5,625,0.08951763);
   grae->SetPointError(5,0,0,0.1408225,0.1408225);
   grae->SetPoint(6,725,0.2678406);
   grae->SetPointError(6,0,0,0.2231914,0.2231914);
   grae->SetPoint(7,825,0.3643913);
   grae->SetPointError(7,0,0,0.3319153,0.3319153);
   grae->SetPoint(8,925,0.5954869);
   grae->SetPointError(8,0,0,0.4613771,0.4613771);

    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y));
      if ( x < 275. ) { y = -1000.; }
      values[index][i] = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.08146677);
   grae->SetPointError(1,0,0,0.1133193,0.1133193);
   grae->SetPoint(2,350,-0.2576047);
   grae->SetPointError(2,0,0,0.1638928,0.1638928);
   grae->SetPoint(3,425,0.1362042);
   grae->SetPointError(3,0,0,0.1715621,0.1715621);
   grae->SetPoint(4,525,0.5536013);
   grae->SetPointError(4,0,0,0.2274302,0.2274302);
   grae->SetPoint(5,625,0.09573878);
   grae->SetPointError(5,0,0,0.2428634,0.2428634);
   grae->SetPoint(6,725,-0.2297499);
   grae->SetPointError(6,0,0,0.3585898,0.3585898);
   grae->SetPoint(7,825,-0.1605828);
   grae->SetPointError(7,0,0,0.4923044,0.4923044);
   grae->SetPoint(8,925,0.1883984);
   grae->SetPointError(8,0,0,0.4587862,0.4587862);

      for ( uint i = 0; i < nbins; ++i ) {
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
  } else if ( set_option == 130 ) {

    // TTbar samples 

    energy = "8";
    lumi = "11.7";

    int nhistos = 4;
    
    syst.clear();
    syst.push_back(0.10);
    syst.push_back(0.10);
    syst.push_back(0.10);
    syst.push_back(0.20);
    syst.push_back(0.30);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

  Int_t nbins1 = 8;
  Double_t bins1[8] = {300.,350.,425.,525.,625.,725.,825.,925.};

  TH1F *Muon_btag_eq0_category_eq2_and_3 = new TH1F("Muon_btag_eq0_category_eq2_and_3","",nbins1,bins1);
  Muon_btag_eq0_category_eq2_and_3->SetBinContent(1,0.9931067);
  Muon_btag_eq0_category_eq2_and_3->SetBinContent(2,0.9931463);
  Muon_btag_eq0_category_eq2_and_3->SetBinContent(3,0.9939956);
  Muon_btag_eq0_category_eq2_and_3->SetBinContent(4,1.008251);
  Muon_btag_eq0_category_eq2_and_3->SetBinContent(5,0.9920727);
  Muon_btag_eq0_category_eq2_and_3->SetBinContent(6,1.025841);
  Muon_btag_eq0_category_eq2_and_3->SetBinContent(7,1.000836);
  Muon_btag_eq0_category_eq2_and_3->SetBinContent(8,0.9630669);
  Muon_btag_eq0_category_eq2_and_3->SetBinError(1,0.05116221);
  Muon_btag_eq0_category_eq2_and_3->SetBinError(2,0.0367912);
  Muon_btag_eq0_category_eq2_and_3->SetBinError(3,0.02956191);
  Muon_btag_eq0_category_eq2_and_3->SetBinError(4,0.06134912);
  Muon_btag_eq0_category_eq2_and_3->SetBinError(5,0.1072359);
  Muon_btag_eq0_category_eq2_and_3->SetBinError(6,0.1476388);
  Muon_btag_eq0_category_eq2_and_3->SetBinError(7,0.1906373);
  Muon_btag_eq0_category_eq2_and_3->SetBinError(8,0.2126528);

  TH1F *Muon_btag_gr0_category_eq2_and_3 = new TH1F("Muon_btag_gr0_category_eq2_and_3","",nbins1,bins1);
  Muon_btag_gr0_category_eq2_and_3->SetBinContent(1,1.005801);
  Muon_btag_gr0_category_eq2_and_3->SetBinContent(2,1.039622);
  Muon_btag_gr0_category_eq2_and_3->SetBinContent(3,1.070993);
  Muon_btag_gr0_category_eq2_and_3->SetBinContent(4,0.923012);
  Muon_btag_gr0_category_eq2_and_3->SetBinContent(5,0.9483998);
  Muon_btag_gr0_category_eq2_and_3->SetBinContent(6,1.031156);
  Muon_btag_gr0_category_eq2_and_3->SetBinContent(7,0.871369);
  Muon_btag_gr0_category_eq2_and_3->SetBinContent(8,1.156432);
  Muon_btag_gr0_category_eq2_and_3->SetBinError(1,0.02861841);
  Muon_btag_gr0_category_eq2_and_3->SetBinError(2,0.03269337);
  Muon_btag_gr0_category_eq2_and_3->SetBinError(3,0.03756606);
  Muon_btag_gr0_category_eq2_and_3->SetBinError(4,0.06562057);
  Muon_btag_gr0_category_eq2_and_3->SetBinError(5,0.1267605);
  Muon_btag_gr0_category_eq2_and_3->SetBinError(6,0.1864612);
  Muon_btag_gr0_category_eq2_and_3->SetBinError(7,0.1863392);
  Muon_btag_gr0_category_eq2_and_3->SetBinError(8,0.2607779);

  TH1F *Muon_btag_eq0_category_greq4 = new TH1F("Muon_btag_eq0_category_greq4","",nbins1,bins1);
  Muon_btag_eq0_category_greq4->SetBinContent(1,0.9921936);
  Muon_btag_eq0_category_greq4->SetBinContent(2,0.9920093);
  Muon_btag_eq0_category_greq4->SetBinContent(3,0.9909783);
  Muon_btag_eq0_category_greq4->SetBinContent(4,1.0409);
  Muon_btag_eq0_category_greq4->SetBinContent(5,1.012983);
  Muon_btag_eq0_category_greq4->SetBinContent(6,0.9361136);
  Muon_btag_eq0_category_greq4->SetBinContent(7,1.029574);
  Muon_btag_eq0_category_greq4->SetBinContent(8,0.9005623);
  Muon_btag_eq0_category_greq4->SetBinError(1,0.09571537);
  Muon_btag_eq0_category_greq4->SetBinError(2,0.09049988);
  Muon_btag_eq0_category_greq4->SetBinError(3,0.05452576);
  Muon_btag_eq0_category_greq4->SetBinError(4,0.07158675);
  Muon_btag_eq0_category_greq4->SetBinError(5,0.1135813);
  Muon_btag_eq0_category_greq4->SetBinError(6,0.1239103);
  Muon_btag_eq0_category_greq4->SetBinError(7,0.186428);
  Muon_btag_eq0_category_greq4->SetBinError(8,0.169362);

  TH1F *Muon_btag_gr0_category_greq4 = new TH1F("Muon_btag_gr0_category_greq4","",nbins1,bins1);
  Muon_btag_gr0_category_greq4->SetBinContent(1,1.009486);
  Muon_btag_gr0_category_greq4->SetBinContent(2,0.9679167);
  Muon_btag_gr0_category_greq4->SetBinContent(3,1.026647);
  Muon_btag_gr0_category_greq4->SetBinContent(4,0.9557768);
  Muon_btag_gr0_category_greq4->SetBinContent(5,0.9665126);
  Muon_btag_gr0_category_greq4->SetBinContent(6,0.8542372);
  Muon_btag_gr0_category_greq4->SetBinContent(7,0.7857504);
  Muon_btag_gr0_category_greq4->SetBinContent(8,1.389679);
  Muon_btag_gr0_category_greq4->SetBinError(1,0.02788042);
  Muon_btag_gr0_category_greq4->SetBinError(2,0.03294511);
  Muon_btag_gr0_category_greq4->SetBinError(3,0.03741779);
  Muon_btag_gr0_category_greq4->SetBinError(4,0.05727073);
  Muon_btag_gr0_category_greq4->SetBinError(5,0.09880627);
  Muon_btag_gr0_category_greq4->SetBinError(6,0.102831);
  Muon_btag_gr0_category_greq4->SetBinError(7,0.1199039);
  Muon_btag_gr0_category_greq4->SetBinError(8,0.2248223);

  // HT BINNING

//   Muon_btag_eq0_category_eq2_and_3->SetBins(nbins1,bins1);
//   Muon_btag_gr0_category_eq2_and_3->SetBins(nbins1,bins1);
//   Muon_btag_eq0_category_greq4->SetBins(nbins1,bins1);
//   Muon_btag_gr0_category_greq4->SetBins(nbins1,bins1);

  uint index = 0;
    TGraphAsymmErrors* grae = 0;

    grae = new TGraphAsymmErrors(Muon_btag_eq0_category_eq2_and_3);
    cout << " test " << Muon_btag_eq0_category_eq2_and_3->GetXaxis()->GetNbins()
	 << " " << grae->GetXaxis()->GetNbins()
	 << endl;

    grae->SetName("Graph");
    grae->SetTitle("#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    grae = new TGraphAsymmErrors(Muon_btag_gr0_category_eq2_and_3);
    grae->SetName("Graph");
    grae->SetTitle("0 b tags #rightarrow 1 b tag (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;
 
    grae = new TGraphAsymmErrors(Muon_btag_eq0_category_greq4);
    grae->SetName("Graph");
    grae->SetTitle("1 b tag #rightarrow 2 b tags (#mu + jets)");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;
    
    // mu->mumu, no aT
    grae = new TGraphAsymmErrors(Muon_btag_gr0_category_greq4);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets #rightarrow #mu#mu + jets");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;
    
    return nhistos;


    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
  } else if ( set_option == 31 ) {

    // 11/fb, alphaT acceptance

    energy = "8";
    lumi = "11.7";

    int nhistos = 3;
    
    syst.clear();
    syst.push_back(0.10);
    syst.push_back(0.20);
    syst.push_back(0.30);
    
    init(nhistos,nbins,values);
    init(nhistos,nbins,errors);

    uint index = 0;
    TGraphAsymmErrors* grae = 0;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets, 2 #leq N_{jet} #leq 3");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.0548582);
   grae->SetPointError(1,0,0,0.03007567,0.03007567);
   grae->SetPoint(2,350,-0.0603661);
   grae->SetPointError(2,0,0,0.03686175,0.03686175);
   grae->SetPoint(3,425,-0.03576732);
   grae->SetPointError(3,0,0,0.0417654,0.0417654);
   grae->SetPoint(4,525,-0.07323264);
   grae->SetPointError(4,0,0,0.07437727,0.07437727);
   grae->SetPoint(5,625,-0.1400356);
   grae->SetPointError(5,0,0,0.1269372,0.1269372);
   grae->SetPoint(6,725,0.03079586);
   grae->SetPointError(6,0,0,0.225785,0.225785);
   grae->SetPoint(7,825,-0.3353326);
   grae->SetPointError(7,0,0,0.3770813,0.3770813);
   grae->SetPoint(8,925,-0.5314985);
   grae->SetPointError(8,0,0,0.4784117,0.4784117);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu + jets, N_{jet} #geq 4");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.04661271);
   grae->SetPointError(1,0,0,0.05311414,0.05311414);
   grae->SetPoint(2,350,-0.1354332);
   grae->SetPointError(2,0,0,0.06862723,0.06862723);
   grae->SetPoint(3,425,-0.06136346);
   grae->SetPointError(3,0,0,0.07538184,0.07538184);
   grae->SetPoint(4,525,0.1012961);
   grae->SetPointError(4,0,0,0.0958097,0.0958097);
   grae->SetPoint(5,625,-0.1610576);
   grae->SetPointError(5,0,0,0.1433355,0.1433355);
   grae->SetPoint(6,725,0.2924913);
   grae->SetPointError(6,0,0,0.2395603,0.2395603);
   grae->SetPoint(7,825,0.4952695);
   grae->SetPointError(7,0,0,0.4328884,0.4328884);
   grae->SetPoint(8,925,-0.6620565);
   grae->SetPointError(8,0,0,0.541967,0.541967);
    
    for ( uint i = 0; i < nbins; ++i ) {
      double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
      errors[index][i] = grae->GetErrorYhigh(i+1);
    }
    titles.push_back( grae->GetTitle() );
    index++;
 

    // mu: no aT -> with aT
    grae = new TGraphAsymmErrors(nbins);
    grae->SetName("Graph");
    grae->SetTitle("#mu#mu + jets, 2 #leq N_{jet} #leq 3");
    grae->SetFillColor(1);
    grae->SetLineWidth(3);
    grae->SetMarkerStyle(20);
    grae->SetMarkerSize(1.5);
   grae->SetPoint(1,300,-0.005709641);
   grae->SetPointError(1,0,0,0.08605346,0.08605346);
   grae->SetPoint(2,350,0.1633132);
   grae->SetPointError(2,0,0,0.1321284,0.1321284);
   grae->SetPoint(3,425,-0.1149851);
   grae->SetPointError(3,0,0,0.1332386,0.1332386);
   grae->SetPoint(4,525,-0.01579951);
   grae->SetPointError(4,0,0,0.2413322,0.2413322);
   grae->SetPoint(5,625,0.1833385);
   grae->SetPointError(5,0,0,0.4501614,0.4501614);
   grae->SetPoint(6,725,-0.2603421);
   grae->SetPointError(6,0,0,0.7373136,0.7373136);
   grae->SetPoint(7,825,-1);
   grae->SetPointError(7,0,0,0,0);
   grae->SetPoint(8,925,4.169746);
   grae->SetPointError(8,0,0,10.17011,10.17011);
   
      for ( uint i = 0; i < nbins; ++i ) {
        double x = 0., y = 0.; (((TGraphAsymmErrors*)(grae))->GetPoint(i+1,x,y)); values[index][i] = y; x = y;
        errors[index][i] = grae->GetErrorYhigh(i+1);
      }
      titles.push_back( grae->GetTitle() );
      index++;

    return nhistos;

    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
  } else if ( set_option == 2011 ) {

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

  }

  return 0;

}

