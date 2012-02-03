
/***************************************************************************** 
  * Project: RooFit                                                           * 
  *                                                                           * 
  * This code was autogenerated by RooClassFactory                            * 
  *****************************************************************************/ 

 // Your description goes here... 

 #include "Riostream.h" 

 #include "eRooWPlus.h" 
 #include "RooAbsReal.h" 
 #include "RooAbsCategory.h" 
 #include <math.h> 
 #include "TMath.h" 

 ClassImp(eRooWPlus) 

   eRooWPlus::eRooWPlus(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _fL,
                 	RooAbsReal& _fR
                        ) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   fL("fL","fL",this,_fL),
   fR("fR","fR",this,_fR)

 { 
 } 


 eRooWPlus::eRooWPlus(const eRooWPlus& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   fL("fL",this,other.fL),
   fR("fR",this,other.fR)

 { 

   rbin=10;

   //   dir = "lp_postmht_mt50_chargeID_sig/";
   dir = "RECO_ElPolPlots_PostMHT50MT50/";

   pHist1 = dir+"RECO_ICVarPFPlus_VertexReweighting_LH";
   pHist2 = dir+"RECO_ICVarPFPlus_VertexReweighting_RH";
   pHist3 = dir+"RECO_ICVarPFPlus_VertexReweighting_LO";

   fW=new TFile("/vols/cms02/georgia/WPol/eWPol_WJets-madgraph_fall10_pile_up.root");

   mc1 = (TH1D*)fW->Get(pHist1); mc1->Rebin(rbin);
   mc2 = (TH1D*)fW->Get(pHist2); mc2->Rebin(rbin);
   mc3 = (TH1D*)fW->Get(pHist3); mc3->Rebin(rbin);

 } 

 Double_t eRooWPlus::evaluate() const 
 { 

   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
   Double_t f;
   Double_t xx = x;

   // Helicity templates
   Double_t m1val = mc1->GetBinContent(mc1->FindBin(xx));
   Double_t m2val = mc2->GetBinContent(mc2->FindBin(xx));
   Double_t m3val = mc3->GetBinContent(mc3->FindBin(xx));

   // Signal PDF
   //f = (fL * m1val) + (fR * m2val) + ((1-fL-fR) * m3val); 
   f = (0.5*fL*(m1val-m2val)) + (0.5*(1-fR)*(m1val+m2val)) + (fR*m3val);

   return f; 

 } 



