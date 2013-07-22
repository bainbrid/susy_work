#!/usr/bin/env python
from ROOT import *
import ROOT as r
import array

hists = ["HT_vs_AlphaT_","LeadJet_vs_SecondJet_"]
#hists = ["HT_vs_AlphaT_"]

jetmult = ["2","3","all"]
#slices = ["Slice_one_","Slice_two_","Slice_three_"]
slices = [""]
dirs = ["200_275","275_325","325_375","375_475","475_575","575_675","675_775","775_875","875"]
#dirs = ["200_275","275_325","325_375","375_475"]
EWK = ["Had_EWK.root"][0]
SMS = ["T1_700_100","T2_700_100", "T1_700_100","T2_500_300","T2cc_175_165","T2cc_175_95"][2:] 
Lumi = 200.0

syst_dict = {"175":0.1,"200":0.1,"275":0.1,"325":0.1,"375":0.1,"475":0.1,"575":0.2,"675":0.2,"775":0.2,"875":0.3}
syst_dict_small = {"175":0.05,"200":0.05,"275":0.05,"325":0.05,"375":0.05,"475":0.05,"575":0.2,"675":0.2,"775":0.2,"875":0.3}

#Fill the root files

for sig in SMS:
  
  file = r.TFile("./%s_ratio.root"%sig ,"RECREATE")
  for model in slices:
    for dir in dirs:
      file.cd("/")
      file.mkdir(str(model)+str(dir))

  sig_model = r.TFile.Open("./"+sig+".root")
  back_model = r.TFile.Open("./"+EWK)
  for model in slices:
    for dir in dirs: 
      for plot in hists: 
       for jet in jetmult:
         
         normal_syst = syst_dict[dir.split("_")[0]]
         small_syst = syst_dict_small[dir.split("_")[0]]

         sig_plot = sig_model.Get(model+dir+"/"+plot+jet)
         back_plot = back_model.Get(model+dir+"/"+plot+jet)
         sig_plot.Scale(Lumi)
         back_plot.Scale(Lumi)


         slice_dict = {}
         ht_list = [175.,225.,275.,325.,375.,475.,575.,675.,775.,875.]

         if plot == "LeadJet_vs_SecondJet_":
         #slice_dict = {"1":{"Low":0.55,"High":1.0},"2":{"Low":1.0,"High":3.0},"3":{"Low":3.0,"High":10.0}  }
            slice_dict = {"1":{"Low":50.,"High":60.0},"2":{"Low":60.0,"High":70.0},"3":{"Low":70.0,"High":80.0},"4":{"Low":80.0,"High":90.0},"5":{"Low":90.0,"High":100.0},"6":{"Low":100.0,"High":1000.0}}


         if plot == "HT_vs_AlphaT_":
            #slice_dict = {"1":{"Low":0.55,"High":1.0},"2":{"Low":1.0,"High":3.0},"3":{"Low":3.0,"High":10.0}  }

            slice_dict = {"1":{"Low":0.55,"High":0.60},"2":{"Low":0.60,"High":0.7},"3":{"Low":0.7,"High":0.8},"4":{"Low":0.8,"High":0.9},"5":{"Low":0.9,"High":1.0},"6":{"Low":1.0,"High":2.0},"7":{"Low":2.0,"High":3.0},"8":{"Low":3.0,"High":10.0}  }

            slice_list = [0.55,0.60,0.7,0.8,0.9,1.0,2.0,3.0,10.0]
            #slice_dict = {"1":{"Low":0.55,"High":1.0},"2":{"Low":1.0,"High":1.5},"3":{"Low":1.5,"High":2.0} ,"4":{"Low":2.0,"High":2.5} ,"5":{"Low":2.5,"High":3.0},"6":{"Low":3.0,"High":3.5},"7":{"Low":3.5,"High":4.0},"8":{"Low":4.0,"High":4.5},"9":{"Low":4.5,"High":5.0},"10":{"Low":5.0,"High":10.0}}


         if plot == "HT_vs_SecondJetPt_":

            slice_dict = {"1":{"Low":50.,"High":60.0},"2":{"Low":60.0,"High":70.0},"3":{"Low":70.0,"High":80.0},"4":{"Low":80.0,"High":90.0},"5":{"Low":90.0,"High":100.0},"6":{"Low":100.0,"High":1000.0}}

         if plot =="HT_vs_AlphaT_":
            
            s_b = r.TH2D("HT_vs_AlphaT_","",len(ht_list)-1,array.array('d',ht_list),len(slice_list)-1,array.array('d',slice_list))
            s_sqrtb = r.TH2D("HT_vs_AlphaT_","",len(ht_list)-1,array.array('d',ht_list),len(slice_list)-1,array.array('d',slice_list))
            s_sqrtb_syst = r.TH2D("HT_vs_AlphaT_","",len(ht_list)-1,array.array('d',ht_list),len(slice_list)-1,array.array('d',slice_list))
            s_sqrtb_syst_small = r.TH2D("HT_vs_AlphaT_","",len(ht_list)-1,array.array('d',ht_list),len(slice_list)-1,array.array('d',slice_list))
 
         else :

           s_b = sig_plot.Clone()
           s_sqrtb = sig_plot.Clone()
           s_sqrtb_syst = sig_plot.Clone()
           s_sqrtb_syst_small = sig_plot.Clone()

         s_b.SetName("%ss_b_%s"%(plot,jet))
         s_sqrtb.SetName("%ss_sqrtb_%s"%(plot,jet))
         s_sqrtb_syst.SetName("%ss_sqrtb_syst_%s"%(plot,jet))
         s_sqrtb_syst_small.SetName("%ss_sqrtb_syst_small_%s"%(plot,jet))

         print sig,dir,plot
         file.cd(model+dir)

         
         for xbin in range(1,sig_plot.GetNbinsX()+1):
           int_dict = {}
           int_x = xbin
           if plot == "LeadJet_vs_SecondJet_": int_x = sig_plot.GetNbinsX()+1

           for entry in sorted(slice_dict.iterkeys()):

              if plot == "LeadJet_vs_SecondJet_" or plot == "HT_vs_SecondJet_":
                int_dict[entry+"_sig"] = sig_plot.Integral(xbin,int_x,sig_plot.GetYaxis().FindBin(slice_dict[entry]["Low"]),sig_plot.GetNbinsY()+1)
                int_dict[entry+"_back"] = back_plot.Integral(xbin,int_x,back_plot.GetYaxis().FindBin(slice_dict[entry]["Low"]),back_plot.GetNbinsY()+1)

              if plot == "HT_vs_AlphaT_":
                int_dict[entry+"_sig"] = sig_plot.Integral(xbin,int_x,sig_plot.GetYaxis().FindBin(slice_dict[entry]["Low"]), sig_plot.GetYaxis().FindBin(slice_dict[entry]["High"])-1)
                int_dict[entry+"_back"] = back_plot.Integral(xbin,int_x,back_plot.GetYaxis().FindBin(slice_dict[entry]["Low"]), sig_plot.GetYaxis().FindBin(slice_dict[entry]["High"])-1)

           print sig_plot.GetYaxis().FindBin(slice_dict[entry]["Low"]), sig_plot.GetYaxis().FindBin(slice_dict[entry]["High"])-1
           print int_dict 


           for ybin in range(1,s_b.GetNbinsY()+1):
              binedge = s_b.GetYaxis().GetBinLowEdge(ybin)
              slice_to_use = 0
              for entry in sorted(slice_dict.iterkeys()):
                if binedge >= slice_dict[entry]["Low"] and binedge < slice_dict[entry]["High"]:
                  slice_to_use = entry
             
              if slice_to_use == 0 : continue


              try:
                  if  int_dict[str(slice_to_use)+"_back"] != 0 and int_dict[str(slice_to_use)+"_sig"] != 0: s_b.SetBinContent(xbin,ybin,int_dict[str(slice_to_use)+"_sig"]/int_dict[str(slice_to_use)+"_back"])
                  else : s_b.SetBinContent(xbin,ybin,0) 
              except ZeroDivisionError: s_b.SetBinContent(xbin,ybin,0) 
              
              try:
                  if  int_dict[str(slice_to_use)+"_back"] != 0 and int_dict[str(slice_to_use)+"_sig"] != 0: s_sqrtb.SetBinContent(xbin,ybin,int_dict[str(slice_to_use)+"_sig"]/sqrt(int_dict[str(slice_to_use)+"_back"]))
                  else : s_sqrtb.SetBinContent(xbin,ybin,0) 
              except ZeroDivisionError: s_sqrtb.SetBinContent(xbin,ybin,0) 

              try:
                  if  int_dict[str(slice_to_use)+"_back"] != 0 and int_dict[str(slice_to_use)+"_sig"] != 0 : s_sqrtb_syst.SetBinContent(xbin,ybin,int_dict[str(slice_to_use)+"_sig"]/sqrt((int_dict[str(slice_to_use)+"_back"]+pow(normal_syst*int_dict[str(slice_to_use)+"_back"] ,2))))
                  else : s_sqrtb_syst.SetBinContent(xbin,ybin,0) 
              except ZeroDivisionError: s_sqrtb_syst.SetBinContent(xbin,ybin,0) 

              try:
                  if  int_dict[str(slice_to_use)+"_back"] != 0 and int_dict[str(slice_to_use)+"_sig"] != 0  : s_sqrtb_syst_small.SetBinContent(xbin,ybin,int_dict[str(slice_to_use)+"_sig"]/sqrt((int_dict[str(slice_to_use)+"_back"]+pow(small_syst*int_dict[str(slice_to_use)+"_back"],2))))
                  else : s_sqrtb_syst_small.SetBinContent(xbin,ybin,0) 
              except ZeroDivisionError: s_sqrtb_syst_small.SetBinContent(xbin,ybin,0) 

         s_b.Write("",r.TObject.kOverwrite)
         s_sqrtb.Write("",r.TObject.kOverwrite)
         s_sqrtb_syst.Write("",r.TObject.kOverwrite)
         s_sqrtb_syst_small.Write("",r.TObject.kOverwrite)

         file.cd("../")
  file.Close()
