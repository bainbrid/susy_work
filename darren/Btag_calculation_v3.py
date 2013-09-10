#!/usr/bin/env python
from ROOT import *
import ROOT as r
import logging,itertools
import os,fnmatch,sys
import glob, errno
from time import strftime
from optparse import OptionParser
import array, ast
from math import *
from plottingUtils import *
import time

class Btag_Calc(object):
        
   def __init__(self,settings,samples,btag_measure,number,alphaT_check,lumi_dict,analysis,analysis_category):
        print "In btag Calc taking samples"

        print "\n\n Making Predictions for %s \n\n" %number
        self.Keep_AlphaT = alphaT_check
        self.settings = settings
        self.Analysis = analysis
        self.analysis_category = analysis_category
        self.btag_multiplicity = btag_measure
        self.Btag_Rate(btag_measure)
        
        #Use to produce Btag - Mistage plots
        #self.DiMuon_Fit(self.Btag_Efficiencies,"Muon_Z0","Mistag")
        #self.DiMuon_Fit(self.Btag_Efficiencies,"Muon_Z2","Mistag")
        #self.DiMuon_Fit(self.Btag_Efficiencies,"Muon_Z0","Btag")
        #self.DiMuon_Fit(self.Btag_Efficiencies,"Muon_Z2","Btag")
        self.lumi_dict = lumi_dict

   def DiMuon_Fit(self,dictionary,sample,tag):
        
        c1= r.TCanvas("Yields", "Yields",0,0,900,600)
        c1.cd()
        fit = r.TF1("fit","pol0",float(self.settings["dirs"][0]) - 50.0,float(self.settings["dirs"][:-1])+50.0)
        data = r.TGraphAsymmErrors(8)
        i = 0
        for num,entry in sorted(dictionary[sample].iteritems()):
                i+=1
                data.SetPoint(i,float(num),dictionary[sample][num]["%s_Efficiency"%tag ])
                data.SetPointEYhigh(i,dictionary[sample][num]["%s_Error"%tag])
                data.SetPointEYlow(i,dictionary[sample][num]["%s_Error"%tag])

        data.GetXaxis().SetRangeUser(float(self.settings["dirs"][0])- 50.,float(self.settings["dirs"][:-1])+50.0)
        if tag == "Mistag" :data.GetYaxis().SetRangeUser(0,0.1)
        if tag == "Btag" : data.GetYaxis().SetRangeUser(0.5,1.0)
        data.SetTitle("")
        data.GetXaxis().SetTitle("H_{T} (GeV)")
        data.GetYaxis().SetTitle("%s Efficiency"%tag)
        data.GetYaxis().SetTitleOffset(1.4)
        data.GetXaxis().SetTitleOffset(0.8)
        data.SetLineWidth(3)
        data.SetMarkerStyle(20)
        data.SetMarkerSize(1.5)
        gStyle.SetOptStat(0)
        data.Draw("AP")
        #data.Fit(fit)
        c1.SaveAs("Uncorrected_%s_%s_Fit.png"%(sample,tag) )
        #print fit.GetParameter(0)
        #for num,entry in sorted(dictionary[sample].iteritems()):
        #        self.Btag_Efficiencies[sample][num]["Mistag_Efficiency"] = fit.GetParameter(0)

   def Make_Dict(self,settings,samples,number):
        
        htbins = settings["dirs"]

        table_entries = "{"
        for key,fi in sorted(samples.iteritems()):
              file = r.TFile("%s.root"%fi[0])
              
              DirKeys = file.GetListOfKeys()  
              i = 0
              for dir in settings['dirs']:
                fixed_dir = dir
                for alphat in settings['AlphaTSlices']:
                   dir = fixed_dir
                   lower = alphat.split('_')[0]
                   higher = alphat.split('_')[1]
                   table_entries += "\t\"%s_%d\"  : "%(key,i)
                   i+=1
                   table_entries += "{\"HT\":\"%s\","%(dir.split('_')[0])
                   checkht = dir
                   for entry in DirKeys:
                        subdirect = file.FindObjectAny(entry.GetName())
                        sample_dir = fi[1]+dir 
                        subdirect.GetName()
                         
                        if sample_dir == subdirect.GetName() and fi[2] != "Data":  
                            for subkey in subdirect.GetListOfKeys():
                                if fi[3] =="Had" or fi[3] == "Photon" or (str(lower) == "0.55" and self.Keep_AlphaT == "True") :
                                    if subkey.GetName() == "Matched_vs_Matched_noB_vs_c_alphaT_%s" %self.analysis_category: 
                                      plot = file.Get(sample_dir+"/"+subkey.GetName())
                                else: 
                                    if subkey.GetName() == "Matched_vs_Matched_noB_vs_c_%s" %self.analysis_category: 
                                      plot = file.Get(sample_dir+"/"+subkey.GetName())
                                
                   if fi[2] != "Data": 
                        table_entries += self.Make_Prediction(plot,fi[3],fi[2],number,dir.split('_')[0],lower)
                        del plot
                   else: table_entries += self.Data_Yield(fi[0],fi[1],dir,lower,higher,fi[2],fi[3])
              file.Close()
        table_entries += "}"
        return_dict = ast.literal_eval(table_entries)
        return return_dict
 
   
   def Data_Yield(self,file,dir_path,dir,lower,higher,sample_type,category):
      for histName in self.settings['plots']:
          histName = str(histName+self.analysis_category)
          checkht = dir
          dir = dir_path+dir
          normal =  GetSumHist(File = ["%s.root"%file], Directories = [dir], Hist = histName, Col = r.kBlack, Norm = None, LegendText = "nBtag")  
          normal.HideOverFlow()
          if self.Keep_AlphaT == "True" and str(lower)=="0.55": 
            err = r.Double(0.0)
            if category != "Had": normal.hObj.IntegralAndError(int(float(lower)/0.01)+1,int(float(higher)/0.01),err)
            else: normal.hObj.IntegralAndError(int(float(0.55)/0.01)+1,int(float(higher)/0.01),err)
            table_string =" \"Yield\": %.3e ,\"Error\":\"%s\",\"SampleType\":\"%s\",\"Category\":\"%s\",\"AlphaT\":%s},\n"%((normal.hObj.Integral(int(float(0.55)/0.01)+1,int(float(10)/0.01)) if category =="Had" else (normal.hObj.Integral(int(float(lower)/0.01)+1,int(float(higher)/0.01)))),err,sample_type,category,lower)
          else:
            err = r.Double(0.0)
            if category == "Had" :normal.hObj.IntegralAndError(int(float(lower)/0.01)+1,int(float(higher)/0.01),err)
            else: normal.hObj.IntegralAndError(1,2000,err)
            table_string =" \"Yield\": %.3e ,\"Error\":\"%s\",\"SampleType\":\"%s\",\"Category\":\"%s\",\"AlphaT\":%s},\n"%((normal.hObj.Integral(int(float(0.55)/0.01)+1,int(float(higher)/0.01)) if category =="Had" else (normal.hObj.Integral())),err,sample_type,category,lower)
          normal.a.Close()
      return table_string    

   def Make_Prediction(self,plot,sample,category,btag_number,htbin,alphaT):  
      
        def bcombo(b, s,charm, e, m,c, hist):

              Nb = b;
              Ns = s;
              Nc = charm;

              #here you set the upper limits for the loop...
              Nbmax = 4 #hist.GetNbinsX()
              Nsmax = hist.GetNbinsY()
              Ncmax = 4 #hist.GetNbinsZ()

              #this is the result to return...
              final_yield = 0.0
              final_error = 0.0
              for x in range(Nb,Nbmax):
                for y in range(Ns,Nsmax):
                  for z in range(Nc,Ncmax):
                    if hist.GetBinContent(x+1,y+1,z+1) == 0 : continue
                    final_yield += hist.GetBinContent(x+1, y+1,z+1) * TMath.Binomial(x,b) * pow(e,b) * pow(1.0 - e, x-b) * TMath.Binomial(y,s) * pow(m,s) * pow(1.0 - m, y-s) * TMath.Binomial(z,charm) * pow(c,charm) * pow(1.0 - c, z-charm)

                    final_error += pow(hist.GetBinError(x+1, y+1,z+1) * TMath.Binomial(x,b) * pow(e,b) * pow(1.0 - e, x-b) * TMath.Binomial(y,s) * pow(m,s) * pow(1.0 - m, y-s) * TMath.Binomial(z,charm) * pow(c,charm) * pow(1.0 - c, z-charm) ,2)
 
              return(final_yield,final_error)

        prediction_dictionary = {"Zero_btags":[0],"One_btag":[1],"Two_btags":[2],"Three_btags":[3],"More_Than_Zero_btag":[1,2,3,4],"More_Than_One_btag":[2,3,4],"Inclusive":[0,1,2,3,4],"More_Than_Three_btag":[4]}
        temp_yield = 0
        temp_yield_error = 0

        Formula_List = prediction_dictionary[btag_number]


        for entry in  Formula_List:

           sample_type = sample
           Znum = ""
           if category in ["Zinv","WJets","DY"]: Znum = "_Z0"
           if category in ["TTbar","DiBoson","SingleTop"] : Znum = "_Z2"
           if sample == "DiMuon": sample_type = "Muon"
           if sample == "Photon": Znum = ""
           if sample in ["OSSF","OSOF","SSSF","SSOF"] : sample_type = "DiLepton"
           
           btag_eff = self.Btag_Efficiencies[sample_type+Znum][htbin]['Btag_Efficiency'] # +  self.Btag_Efficiencies[sample_type+Znum][htbin]['Btag_Error']
           mistag_eff = self.Btag_Efficiencies[sample_type+Znum][htbin]['Mistag_Efficiency']  #+  self.Btag_Efficiencies[sample_type+Znum][htbin]['Mistag_Error']
           ctag_eff = self.Btag_Efficiencies[sample_type+Znum][htbin]['Ctag_Efficiency']  #+  self.Btag_Efficiencies[sample_type+Znum][htbin]['Ctag_Error']

           for j in range(0,entry+1):
             for k in range(0,entry+1):
               for l in range (0,entry+1):
                 if j + k + l == entry:
                   temp_yield += bcombo(j,k,l,btag_eff,mistag_eff,ctag_eff,plot)[0]
                   temp_yield_error += bcombo(j,k,l,btag_eff,mistag_eff,ctag_eff,plot)[1]
           
        yield_pred = temp_yield
        error_pred = sqrt(temp_yield_error)
        Luminosity = self.lumi_dict[sample]
        table_string =" \"Yield\": %.3e ,\"Error\":\"%s\",\"SampleType\":\"%s\",\"Category\":\"%s\",\"AlphaT\":%s},\n"%(yield_pred*(10*Luminosity),error_pred*(Luminosity*10),category,sample,alphaT)
        return table_string
  

   def Btag_Rate(self,btag_measurement):
        self.Btag_Efficiencies = {'Had_Z0':{},'Muon_Z0':{},'DiMuon':{},'Had_Z2':{},'Muon_Z2':{},'Photon':{},'DiLepton_Z0':{},'DiLepton_Z2':{}}
        self.bins = tuple(self.settings["bins"])
        dict_entries = ('Btag_Efficiency','Mistag_Efficiency','Ctag_Efficiency','Btag_Error','Mistag_Error','Ctag_Error')
        for key in self.Btag_Efficiencies:
             self.Btag_Efficiencies[key] = dict.fromkeys(self.bins)
             for a in self.bins: 
                        self.Btag_Efficiencies[key][a] = dict.fromkeys(dict_entries,0)

        htbins = self.settings["dirs"]
        for key,fi in sorted(btag_measurement.iteritems()):
           file = r.TFile.Open(fi[0])
           DirKeys = file.GetListOfKeys()
          
           for num,bin in enumerate(htbins):
             if fi[2] == "Photon_" and bin in ["150_200","200_275","275_325","325_375"] : continue
             for entry in DirKeys:
                subdirect = file.FindObjectAny(entry.GetName())
                dir = fi[2]+bin
                subdirect.GetName()
                if dir == subdirect.GetName():
                   for subkey in [ "GenJetPt_nBgen_all", "GenJetPt_noB_nBgen_all", "GenJetPt_c_nBgen_all", "Btagged_GenJetPt_nBgen_SFb_Medium_all", "Btagged_GenJetPt_noB_nBgen_SFlight_Medium_all", "Btagged_GenJetPt_c_nBgen_SFlight_Medium_all"  ]:
                       #=========================================#
                       if subkey == "GenJetPt_nBgen_all":
                           err = r.Double(0.0)     
                           plot = file.Get(dir+"/"+subkey)
                           self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]["Btag_Efficiency"] = plot.Integral()
                           plot.IntegralAndError(1,10000,err)
                           try:self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]["Btag_Error"] = err/plot.Integral()
                           except ZeroDivisionError : self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]["Btag_Error"] = 0

                       if subkey == "Btagged_GenJetPt_nBgen_SFb_Medium_all":
                       #if subkey == "Btagged_GenJetPt_nBgen_Medium_all":

                           err = r.Double(0.0)
                           bplot = file.Get(dir+"/"+subkey)
                           try: self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]['Btag_Efficiency'] = bplot.Integral()/(self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]['Btag_Efficiency'])
                           except ZeroDivisionError : self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]['Btag_Efficiency'] = 0
                           bplot.IntegralAndError(1,10000,err)
                           try: self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]["Btag_Error"] =  self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]['Btag_Efficiency']*sqrt(pow(self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]["Btag_Error"],2)+pow(err/bplot.Integral(),2))
                           except ZeroDivisionError : self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]["Btag_Error"] = 0


                       #========================================#
                       if subkey == "GenJetPt_noB_nBgen_all":
                           mistag_plot = file.Get(dir+"/"+subkey)
                           err = r.Double(0.0)
                           self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]['Mistag_Efficiency'] =mistag_plot.Integral()
                           plot.IntegralAndError(1,10000,err)
                           self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]["Mistag_Error"] = err/mistag_plot.Integral()
                       
                       if subkey == "Btagged_GenJetPt_noB_nBgen_SFlight_Medium_all":
                       #if subkey == "Btagged_GenJetPt_noB_nBgen_Medium_all":
                                
                           aplot = file.Get(dir+"/"+subkey)
                           err = r.Double(0.0)
                           self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]['Mistag_Efficiency'] = aplot.Integral()/(self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]['Mistag_Efficiency'])
                           aplot.IntegralAndError(1,10000,err)
                           try : self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]["Mistag_Error"] =  self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]['Mistag_Efficiency']*sqrt(pow(self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]["Mistag_Error"],2)+pow(err/aplot.Integral(),2))
                           except ZeroDivisionError : self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]["Mistag_Error"] = 0


                        #=======================================#

                       if subkey == "GenJetPt_c_nBgen_all":
                           ctag_plot = file.Get(dir+"/"+subkey)
                           err = r.Double(0.0)
                           self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]['Ctag_Efficiency'] =ctag_plot.Integral()
                           plot.IntegralAndError(1,10000,err)
                           self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]["Ctag_Error"] = err/ctag_plot.Integral()

                       if subkey == "Btagged_GenJetPt_c_nBgen_SFlight_Medium_all":
                       #if subkey == "Btagged_GenJetPt_c_nBgen_Medium_all":

                           cplot = file.Get(dir+"/"+subkey)
                           err = r.Double(0.0)
                           self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]['Ctag_Efficiency'] = cplot.Integral()/(self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]['Ctag_Efficiency'])
                           aplot.IntegralAndError(1,10000,err)
                           self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]["Ctag_Error"] =  self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]['Ctag_Efficiency']*sqrt(pow(self.Btag_Efficiencies[fi[1]][bin.split('_')[0]]["Ctag_Error"],2)+pow(err/cplot.Integral(),2))


                        #=======================================#


                       
 
