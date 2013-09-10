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
#include <math.h>

using namespace std;
typedef unsigned int uint;
void setTDRStyle();

void qcd2012() {
   
  bool normalise = false;

  setTDRStyle();

  typedef std::vector<double> doubles;
  
  uint ncat = 11;

  std::vector<std::string> cats;
  cats.push_back("2012, 11.7/fb, ge2, QCD MC");
  cats.push_back("2012, 11.7/fb, le3, QCD MC");
  cats.push_back("2012, 11.7/fb, ge4, QCD MC");
  cats.push_back("2012, 11.7/fb, ge2, Darren");
  cats.push_back("2012, 11.7/fb, le3, Darren");
  cats.push_back("2012, 11.7/fb, ge4, Darren");
  cats.push_back("2012, 11.7/fb, ge2, Zhaoxia");
  cats.push_back("2012, 11.7/fb, le3, Zhaoxia");
  cats.push_back("2012, 11.7/fb, ge4, Zhaoxia");
  cats.push_back("2012, 3.9/fb, ge2");
  cats.push_back("2011, 5.0/fb, ge2");
  
  doubles bins;
  //bins.push_back(275.);
  //bins.push_back(325.);
  //for ( int i = 0; i < 6; ++i ) { bins.push_back(375.+100.*i); }
  bins.push_back(298);
  bins.push_back(348);
  bins.push_back(416);
  bins.push_back(517);
  bins.push_back(617);
  bins.push_back(719);
  bins.push_back(819);
  bins.push_back(1044);

  doubles widths;
  widths.push_back(25.);
  widths.push_back(25.);
  for ( int i = 0; i < 6; ++i ) { widths.push_back(50.); }

  std::vector<doubles> values;
  values.resize( ncat, doubles(0,0.) );

  std::vector<doubles> errors;
  errors.resize( ncat, doubles(0,0.) );

  // 2012, HCP, ge2, QCD
  uint icat = 0;
  values[icat].push_back(644595861.);
  values[icat].push_back(274899534.);
  values[icat].push_back(195734721.);
  values[icat].push_back(64548912.);
  values[icat].push_back(23380976.);
  values[icat].push_back(9311077.);
  values[icat].push_back(4064315.);
  values[icat].push_back(3944842.);
  errors[icat].push_back(787304.);
  errors[icat].push_back(410284.);
  errors[icat].push_back(239708.);
  errors[icat].push_back(97840.);
  errors[icat].push_back(42730.);
  errors[icat].push_back(20372.);
  errors[icat].push_back(10643.);
  errors[icat].push_back(7449.);

  // 2012, HCP, le3, QCD
  icat++;
  values[icat].push_back(521996475.);
  values[icat].push_back(225770454.);
  values[icat].push_back(157659684.);
  values[icat].push_back(42608311.);
  values[icat].push_back(13789640.);
  values[icat].push_back(5145155.);
  values[icat].push_back(2134882.);
  values[icat].push_back(1980795.);
  errors[icat].push_back(664993.);
  errors[icat].push_back(334154.);
  errors[icat].push_back(181395.);
  errors[icat].push_back(60079.);
  errors[icat].push_back(25575.);
  errors[icat].push_back(10802.);
  errors[icat].push_back(5084.);
  errors[icat].push_back(2748.);

  // 2012, HCP, ge4, QCD
  icat++;
  values[icat].push_back(122599386.);
  values[icat].push_back(49162577.);
  values[icat].push_back(38108534.);
  values[icat].push_back(21940600.);
  values[icat].push_back(9582404.);
  values[icat].push_back(4167039.);
  values[icat].push_back(1929433.);
  values[icat].push_back(1962930.);
  errors[icat].push_back(421463.);
  errors[icat].push_back(238063.);
  errors[icat].push_back(156702.);
  errors[icat].push_back(77222.);
  errors[icat].push_back(34232.);
  errors[icat].push_back(17272.);
  errors[icat].push_back(9350.);
  errors[icat].push_back(6923.);

  // 2012, HCP, ge2, Darren
  icat++;
  values[icat].push_back(653500000.); 
  values[icat].push_back(294800000.); 
  values[icat].push_back(214600000.); 
  values[icat].push_back(72190000.); 
  values[icat].push_back(26470000.); 
  values[icat].push_back(10860000.); 
  values[icat].push_back(4741000.);  
  values[icat].push_back(4718000.);  
  errors[icat].push_back(1202585.);
  errors[icat].push_back(571151.);
  errors[icat].push_back(344609.);
  errors[icat].push_back(141325.);
  errors[icat].push_back(85562.);
  errors[icat].push_back(54837.);
  errors[icat].push_back(36223.);
  errors[icat].push_back(36122.);

  // 2012, HCP, le3, Darren
  icat++;
  values[icat].push_back(559500000.);
  values[icat].push_back(252400000.);
  values[icat].push_back(180600000.);
  values[icat].push_back(51650000.);
  values[icat].push_back(17060000.);
  values[icat].push_back(6499000.);
  values[icat].push_back(2674000.);
  values[icat].push_back(2501000.);
  errors[icat].push_back(1112701.);
  errors[icat].push_back(528494.);
  errors[icat].push_back(316123.);
  errors[icat].push_back(119511.);
  errors[icat].push_back(68677.);
  errors[icat].push_back(42396.);
  errors[icat].push_back(27183.);
  errors[icat].push_back(26298.);

  // 2012, HCP, ge4, Darren
  icat++;
  values[icat].push_back(93940000.); 
  values[icat].push_back(42330000.); 
  values[icat].push_back(33950000.);
  values[icat].push_back(20540000.); 
  values[icat].push_back(9410000.);  
  values[icat].push_back(4363000.);  
  values[icat].push_back(2067000.);  
  values[icat].push_back(2217000.);  
  errors[icat].push_back(456188.);
  errors[icat].push_back(216580.);
  errors[icat].push_back(137192.);
  errors[icat].push_back(75431.);
  errors[icat].push_back(51032.);
  errors[icat].push_back(34780.);
  errors[icat].push_back(23942.);
  errors[icat].push_back(24763.);

  // 2012, HCP, ge2, Zhaoxia
  icat++;
  values[icat].push_back(630453600.); 
  values[icat].push_back(286166200.); 
  values[icat].push_back(209611400.); 
  values[icat].push_back(69777150.); 
  values[icat].push_back(26101500.); 
  values[icat].push_back(20182300.); 
  values[icat].push_back(4745175.);  
  values[icat].push_back(4776350.);  
  errors[icat].push_back(1181263.323734);
  errors[icat].push_back(562792.359579);  
  errors[icat].push_back(340602.025831);  
  errors[icat].push_back(98254.312882);	   
  errors[icat].push_back(60093.832878);	   
  errors[icat].push_back(52839.982494);	   
  errors[icat].push_back(25625.390241);	   
  errors[icat].push_back(25701.459103);   

  // 2012, HCP, le3, Zhaoxia
  icat++;
  values[icat].push_back(487992800.);
  values[icat].push_back(202369400.);
  values[icat].push_back(134976100.);
  values[icat].push_back(36965375.);
  values[icat].push_back(12292400.);
  values[icat].push_back(8301900.);
  values[icat].push_back(1925125.);
  values[icat].push_back(1768325.);
  errors[icat].push_back(1038846.437160);   
  errors[icat].push_back(472978.096744);    
  errors[icat].push_back(273134.234398);    
  errors[icat].push_back(71457.570453);         
  errors[icat].push_back(41210.101917);         
  errors[icat].push_back(33859.322941); 	        
  errors[icat].push_back(16303.469723); 	        
  errors[icat].push_back(15630.159148);       

  // 2012, HCP, ge4, Zhaoxia
  icat++;
  values[icat].push_back(142460800.);
  values[icat].push_back(83796800.);
  values[icat].push_back(74635300.);
  values[icat].push_back(32811775.);
  values[icat].push_back(13809100.);
  values[icat].push_back(11880400.);
  values[icat].push_back(2820050.);
  values[icat].push_back(3008025.);
  errors[icat].push_back(562299.848835);
  errors[icat].push_back(305003.213098);
  errors[icat].push_back(203488.156904);
  errors[icat].push_back(67436.826920);
  errors[icat].push_back(43737.812588);
  errors[icat].push_back(40566.118868);
  errors[icat].push_back(19770.116338);
  errors[icat].push_back(20402.527417);

  // 2012, ICHEP, ge2
  icat++;
  values[icat].push_back(231496000.);	
  values[icat].push_back(103615000.);
  values[icat].push_back(76347400.);
  values[icat].push_back(25456300.);
  values[icat].push_back(9467480.);
  values[icat].push_back(3855680.);
  values[icat].push_back(1729150.);
  values[icat].push_back(1750550.);
  for ( uint j = 0; j < values[icat].size(); ++j ) {
    errors[icat].push_back( sqrt(values[icat][j]) );
  }
  
  // 2011, 5/fb, ge2
  icat++;
  values[icat].push_back(2.792e+08);
  values[icat].push_back(1.214e+08);
  values[icat].push_back(8.544e+07);
  values[icat].push_back(2.842e+07);
  values[icat].push_back(9.953e+06);
  values[icat].push_back(3.954e+06);
  values[icat].push_back(1.679e+06);
  values[icat].push_back(1.563e+06); 
  for ( uint j = 0; j < values[icat].size(); ++j ) {
    errors[icat].push_back( sqrt(values[icat][j]) );
  }

  // Calc norm, max, min
  double norm = 0.;
  double max = 1.e0;
  double min = 1.e10;
  for ( uint i = 0; i < ncat; ++i ) {
    double temp = 0.;
    for ( uint j = 0; j < values[i].size(); ++j ) {
      values[i][j] /= (widths[j]/100.);
      double tmp = values[i][j];
      temp += tmp; 
      max = ( tmp > max ? tmp : max ); 
      min = ( tmp < min ? tmp : min ); 
    }
    norm = ( temp > norm ? temp : norm ); 
  }

  // Then normalise
  if (normalise) {
    for ( uint i = 0; i < ncat; ++i ) {
      for ( uint j = 0; j < values[i].size(); ++j ) {
	values[i][j] /= norm;
	errors[i][j] /= norm;
      }
    }
  }
  
//   for ( uint j = 0; j < values[0].size(); ++j ) {
//     cout << (values[0][j]/1.) / (values[3][j]/1.) << endl;
//   }
  
  std::vector<doubles> effs;
  std::vector<doubles> errs;

  if (1) {
    effs.resize( ncat, doubles(bins.size(),1.) );
    errs.resize( ncat, doubles(bins.size(),0.) );
  } else {
    
    effs.resize( ncat, doubles(0,0.) );
    errs.resize( ncat, doubles(0,0.) );
    
    effs[0].push_back(1.03);
    effs[0].push_back(1.02);
    effs[0].push_back(1.03);
    effs[0].push_back(0.97);
    effs[0].push_back(1.17);
    effs[0].push_back(1.27);
    effs[0].push_back(1.55);
    effs[0].push_back(1.36);
    errs[0].push_back(0.10);
    errs[0].push_back(0.09);
    errs[0].push_back(0.09);
    errs[0].push_back(0.09);
    errs[0].push_back(0.13);
    errs[0].push_back(0.20);
    errs[0].push_back(0.30);
    errs[0].push_back(0.37);

    effs[1].push_back(1.00);
    effs[1].push_back(0.95);
    effs[1].push_back(1.02);
    effs[1].push_back(0.95);
    effs[1].push_back(1.09);
    effs[1].push_back(1.19);
    effs[1].push_back(1.38);
    effs[1].push_back(1.14);
    errs[1].push_back(0.13);
    errs[1].push_back(0.12);
    errs[1].push_back(0.10);
    errs[1].push_back(0.11);
    errs[1].push_back(0.17);
    errs[1].push_back(0.26);
    errs[1].push_back(0.39);
    errs[1].push_back(0.48);

    effs[2].push_back(1.07);
    effs[2].push_back(1.13);
    effs[2].push_back(1.06);
    effs[2].push_back(1.01);
    effs[2].push_back(1.28);
    effs[2].push_back(1.36);
    effs[2].push_back(1.78);
    effs[2].push_back(1.67);
    errs[2].push_back(0.17);
    errs[2].push_back(0.14);
    errs[2].push_back(0.16);
    errs[2].push_back(0.14);
    errs[2].push_back(0.21);
    errs[2].push_back(0.29);
    errs[2].push_back(0.46);
    errs[2].push_back(0.58);

  }
  
  TCanvas* c = new TCanvas("tmp","tmp",900,600);
  
  std::vector<int> style;
  style.push_back(20);
  style.push_back(20);
  style.push_back(20);
  style.push_back(21);
  style.push_back(21);
  style.push_back(21);
  style.push_back(24);
  style.push_back(24);
  style.push_back(24);
  style.push_back(25);
  style.push_back(26);

  std::vector<int> col;
  col.push_back(1);
  col.push_back(2);
  col.push_back(4);
  col.push_back(1);
  col.push_back(2);
  col.push_back(4);
  col.push_back(1);
  col.push_back(2);
  col.push_back(4);
  col.push_back(1);
  col.push_back(1);

  std::vector<int> size;
  size.push_back(1.0);
  size.push_back(1.0);
  size.push_back(1.0);
  size.push_back(1.0);
  size.push_back(1.0);
  size.push_back(1.0);
  size.push_back(2.0);
  size.push_back(2.0);
  size.push_back(2.0);
  size.push_back(2.0);
  size.push_back(2.0);

  std::vector<TGraphAsymmErrors*> gr;
  for ( uint i = 0; i < ncat; ++i ) {
    gr.push_back( new TGraphAsymmErrors(8) );
    gr.back()->SetName( cats[i].c_str() );
    gr.back()->SetTitle( cats[i].c_str() );
    gr.back()->SetMarkerStyle( style[i] );
    gr.back()->SetMarkerColor( col[i] );
    gr.back()->SetMarkerSize( size[i] );
    gr.back()->SetLineStyle(1);
    gr.back()->SetLineColor( col[i] );
    gr.back()->SetLineWidth(2);
  }

  for ( uint i = 0; i < ncat; ++i ) {
    TGraphAsymmErrors* temp = gr[i];

    for ( uint j = 0; j < bins.size(); ++j ) {
      double x = bins[j];
      double y = effs[i][j] > 0. ? values[i][j] / effs[i][j] : 0.;
      double ey = errors[i][j];
      ey = y*sqrt( (ey/y)*(ey/y) + (errs[i][j]/effs[i][j])*(errs[i][j]/effs[i][j]) );
      
      temp->SetPoint(j+1,x,y);
      temp->SetPointError(j+1,0.,0.,ey,ey);
      
    }
  }
  
  TMultiGraph* mg = new TMultiGraph();
  for ( uint i = 0; i < ncat; ++i ) { mg->Add(gr[i],"pl"); }
  mg->Draw("ap");
  mg->GetXaxis()->SetTitle("Mean H_{T} (GeV)");
  mg->GetYaxis()->SetTitle("Bulk yields / 100 GeV");
  cout << max << " " << min << endl;
  if (normalise) { max /= norm; min /= norm; }
  cout << max << " " << min << endl;
  max = pow(10.,int(log10(max)+1));
  min = pow(10.,int(log10(min)-0));
  cout << max << " " << min << endl;
  mg->GetYaxis()->SetRangeUser(min,max); 
  mg->GetXaxis()->SetRangeUser(275.,1075.);
  
  TLegend* leg = new TLegend( 0.6, 0.9-0.05*ncat, 0.8, 0.9 );
  leg->SetFillColor(0);
  leg->SetLineColor(0); 
  leg->SetShadowColor(0); 
  leg->SetTextSize(0.04);
  for ( uint i = 0; i < ncat; ++i ) { leg->AddEntry(gr[i],cats[i].c_str(),"p"); }
  leg->Draw("same");

  std::stringstream sss;
  sss << "CMS Preliminary, 11.7 fb^{-1}, #sqrt{s} = 8 TeV";
  TLatex* tex = new TLatex(0.5,0.87,sss.str().c_str());
  tex->SetNDC();
  tex->SetTextSize(0.04);
  //tex->Draw();

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
  tdrStyle->SetOptFit(0); // To display the mean and RMS:   SetOptStat("mr");
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

