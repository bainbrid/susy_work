/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOWPLUSEXT
#define ROOWPLUSEXT

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 // ROOT headers
#include "TFile.h"
#include "TH1D.h"
#include "TString.h" 
class RooWPlusExt : public RooAbsPdf {
public:
  RooWPlusExt() {} ; 
  RooWPlusExt(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _f1,
	      RooAbsReal& _f2,
	      RooAbsReal& _fs,
	      RooAbsReal& _fm);
  RooWPlusExt(const RooWPlusExt& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooWPlusExt(*this,newname); }
  inline virtual ~RooWPlusExt() { }

protected:

  RooRealProxy x ;
  RooRealProxy f1 ;
  RooRealProxy f2 ;
  RooRealProxy fs ;
  RooRealProxy fm ;
  
  Double_t evaluate() const ;

private:

 Int_t rbin;
  TString dir;

  TFile *fW; // use W file to extract templates
  TH1D *mc1, *mc2, *mc3;
  
  TString pHist1, pHist2, pHist3;

  ClassDef(RooWPlusExt,1) // Your description goes here...
};
 
#endif
