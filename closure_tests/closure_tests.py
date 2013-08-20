from subprocess import Popen, PIPE, STDOUT
from setTDRStyle import setTDRStyle
from datetime import datetime
import pickle
import pprint
import os
from ROOT import *
import numpy as np
np.set_printoptions(precision=3)
from collections import OrderedDict

################################################################################
# LISTS AND DICTS ##############################################################
################################################################################

master_files_list = [
    "Btag_dimuon_to_dimuon_with_without_alphaT_2.C",
    "Btag_dimuon_to_dimuon_with_without_alphaT_3.C",
    "Btag_dimuon_to_dimuon_with_without_alphaT_all.C",
    "Btag_gamma_to_dimuon_2.C",
    "Btag_gamma_to_dimuon_3.C",
    "Btag_gamma_to_dimuon_all.C",
    "Btag_mu_greater_one_mu_greater_one_with_without_alphaT_2.C",
    "Btag_mu_greater_one_mu_greater_one_with_without_alphaT_3.C",
    "Btag_mu_greater_one_mu_greater_one_with_without_alphaT_all.C",
    "Btag_mu_one_mu_greater_one_no_alphaT_Cut_2.C",
    "Btag_mu_one_mu_greater_one_no_alphaT_Cut_3.C",              
    "Btag_mu_one_mu_greater_one_no_alphaT_Cut_all.C",            
    "Btag_mu_one_mu_one_with_without_alphaT_Cut_2.C",            
    "Btag_mu_one_mu_one_with_without_alphaT_Cut_3.C",            
    "Btag_mu_one_mu_one_with_without_alphaT_Cut_all.C",          
    "Btag_mu_one_mu_two_no_alphaT_Cut_2.C",                      
    "Btag_mu_one_mu_two_no_alphaT_Cut_3.C",                      
    "Btag_mu_one_mu_two_no_alphaT_Cut_all.C",                    
    "Btag_mu_to_dimuon_alphaT_Cut_2.C",                          
    "Btag_mu_to_dimuon_alphaT_Cut_3.C",                          
    "Btag_mu_to_dimuon_alphaT_Cut_all.C",                        
    "Btag_mu_to_dimuon_no_alphaT_Cut_2.C",                       
    "Btag_mu_to_dimuon_no_alphaT_Cut_3.C",                       
    "Btag_mu_to_dimuon_no_alphaT_Cut_all.C",                     
    "Btag_mu_to_mu_with_without_alphaT_2.C",                     
    "Btag_mu_to_mu_with_without_alphaT_3.C",                     
    "Btag_mu_to_mu_with_without_alphaT_all.C",                   
    "Btag_mu_zero_mu_greater_one_no_alphaT_Cut_2.C",             
    "Btag_mu_zero_mu_greater_one_no_alphaT_Cut_3.C",             
    "Btag_mu_zero_mu_greater_one_no_alphaT_Cut_all.C",           
    "Btag_mu_zero_mu_one_no_alphaT_Cut_2.C",                     
    "Btag_mu_zero_mu_one_no_alphaT_Cut_3.C",                     
    "Btag_mu_zero_mu_one_no_alphaT_Cut_all.C",                   
    "Btag_mu_zero_mu_zero_with_without_alphaT_Cut_2.C",          
    "Btag_mu_zero_mu_zero_with_without_alphaT_Cut_3.C",          
    "Btag_mu_zero_mu_zero_with_without_alphaT_Cut_all.C",        
    "Btag_mu_zero_to_dimuon_zero_alphaT_Cut_2.C",                
    "Btag_mu_zero_to_dimuon_zero_alphaT_Cut_3.C",                
    "Btag_mu_zero_to_dimuon_zero_alphaT_Cut_all.C",              
    "Btag_mumu_zero_mmuu_one_no_alphaT_Cut_2.C",                 
    "Btag_mumu_zero_mmuu_one_no_alphaT_Cut_3.C",                 
    "Btag_mumu_zero_mmuu_one_no_alphaT_Cut_all.C",               
    "Btag_one_gamma_to_dimuon_2.C",                              
    "Btag_one_gamma_to_dimuon_3.C",                              
    "Btag_one_gamma_to_dimuon_all.C",                            
    "Btag_zero_gamma_to_inclusive_dimuon_2.C",                   
    "Btag_zero_gamma_to_inclusive_dimuon_3.C",                   
    "Btag_zero_gamma_to_inclusive_dimuon_all.C",                 
    "JetCat_dimuon_to_dimuon.C",                                 
    "JetCat_gamma_to_gamma.C",                                   
    "JetCat_muon_to_muon.C",
    ]

master_titles_dict = {
    "Btag_mu_to_mu_with_without_alphaT_all.C":"#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + #geq 2 jets)",
    "Btag_mu_to_mu_with_without_alphaT_2.C":"#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + 2-3 jets)",
    "Btag_mu_to_mu_with_without_alphaT_3.C":"#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + #geq 4 jets)",
    "Btag_dimuon_to_dimuon_with_without_alphaT_all.C":"#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu#mu + #geq 2 jets)",
    "Btag_dimuon_to_dimuon_with_without_alphaT_2.C":"#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu#mu + 2-3 jets)",
    "Btag_dimuon_to_dimuon_with_without_alphaT_3.C":"#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu#mu + #geq 4 jets)",
    "Btag_mu_zero_mu_one_no_alphaT_Cut_all.C":"0 b tags #rightarrow 1 b tag (#mu + jets)",
    "Btag_mu_zero_mu_one_no_alphaT_Cut_2.C":"0 b tags #rightarrow 1 b tag (#mu + jets)",
    "Btag_mu_zero_mu_one_no_alphaT_Cut_3.C":"0 b tags #rightarrow 1 b tag (#mu + jets)",
    "Btag_mu_one_mu_greater_one_no_alphaT_Cut_all.C":"1 b tag #rightarrow #geq 2 b tags (#mu + jets)",
    "Btag_mu_one_mu_two_no_alphaT_Cut_all.C":"1 b tag #rightarrow 2 b tags (#mu + jets)",
    "Btag_mu_one_mu_two_no_alphaT_Cut_2.C":"1 b tag #rightarrow 2 b tags (#mu + jets)",
    "Btag_mu_one_mu_two_no_alphaT_Cut_3.C":"1 b tag #rightarrow 2 b tags (#mu + jets)",
    "Btag_mu_to_dimuon_no_alphaT_Cut_all.C":"#mu + jets #rightarrow #mu#mu + jets",
    "Btag_mu_to_dimuon_no_alphaT_Cut_2.C":"#mu + jets #rightarrow #mu#mu + jets",
    "Btag_mu_to_dimuon_no_alphaT_Cut_3.C":"#mu + jets #rightarrow #mu#mu + jets",
    "Btag_gamma_to_dimuon_no_alphaT_Cut_all.C":"#mu#mu + jets #rightarrow #gamma + jets", 
    "Btag_gamma_to_dimuon_no_alphaT_Cut_2.C":"#mu#mu + jets #rightarrow #gamma + jets", 
    "Btag_gamma_to_dimuon_no_alphaT_Cut_3.C":"#mu#mu + jets #rightarrow #gamma + jets", 
    "Btag_gamma_to_dimuon_all.C":"#mu#mu + jets #rightarrow #gamma + jets",
    }

################################################################################
# CLASS ########################################################################
################################################################################

class data :

    def __init__(self,
                 paths=[],
                 filenames=[],
                 label="",
                 file=False,
                 energy="8",
                 lumi="19.3",
                 prelim=True,
                 ) :
        # Data
        if len(paths) < len(filenames) :
            paths.extend( paths[-1:]*(len(filenames)- len(paths)) )
        self.paths_ = [ x.append("/")  if x[-1] != "/" else x for x in paths ]
        self.names_ = filenames
        self.files_ = [ x[0]+x[1] for x in zip(paths,filenames) ]
        self.titles_ = master_titles_dict
        self.timestamp_ = datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S')
        if len(label) > 0 : self.label_ = label
        else : self.label_ = "closure_tests_" + self.timestamp_
        self.file_ = file
        self.stream_ = open(self.label_+".log","w")
        # Histogram
        self.nbins_ = 11
        self.bins_ = [200.,275.,325.,375.,475.,575.,675.,775.,875.,975.,1075.,1175.]
        self.regions_ = [0,1,2,3,5,7,11]
        self.markers_ = [24,25,26,5,20,32,27,31]
        self.size_ = [ 1.6 if x != 20 else 1.4 for x in self.markers_ ]
        print self.markers_
        print self.size_
        self.offset_ = [ 5.*x-(len(self.files_)/2) for x in range(len(self.files_)) ]
        self.energy_ = energy
        self.lumi_ = lumi
        self.prelim_ = prelim

    def init(self) :
        if self.file_ is True :
            self.data_ = self.read_from_file(label+".pkl")
        else :
            self.wget(self.files_)
            self.data_ = self.parse_files(self.names_)
            self.data_ = self.patch(self.data_)
            self.write_to_file(self.label_+".pkl")

    def print_pretty(self,dict) :
        f = self.stream_
        for key,value in dict.iteritems() :
            print >>f, key
            for idx,tuple in enumerate(value) :
                print >>f, " ",tuple 
                
    def verbose(self) :
        f = self.stream_
        print >>f, " Read from file?",self.file_
        print >>f, " Timestamp: '%s'"%(self.timestamp_)
        print >>f, " Label: '%s'"%(self.label_)
        print >>f, " Paths:"
        for i in self.paths_ : print >>f, "  ",i 
        print >>f," Names:"
        for i in self.names_ : print >>f, "  ",i 
        print >>f, " Files:"
        for i in self.files_ : print >>f, "  ",i 
        print >>f, " Nbins:", self.nbins_
        print >>f, " Binning:", self.bins_
        print >>f, " Regions", self.regions_
        print >>f, " Contents:"
        self.print_pretty(self.data_)

    def write_to_file(self,file) :
        f = self.stream_
        print >>f, "Writing to file %s"%(file)
        output = open(file,'w')
        pickle.dump(self.data_,output)
            
    def read_from_file(self,file) :
        f = self.stream_
        print >>f, "Reading from file '%s'"%(file)
        output = open(file,'r')
        return pickle.load(output)

    def wget(self,files) :
        f = self.stream_
        if not os.path.exists(self.label_) :
            os.makedirs(self.label_)
            for i in files :
                process = Popen(['wget','-P',self.label_,i], stdout=PIPE,stderr=PIPE)
                stdout,stderr = process.communicate()
                print >>f, stdout, stderr
        else :
            print >>f, "Directory '%s' already exists!"%(self.label_)
            quit()

    def parse_files(self,filenames) :
        dict = OrderedDict()
        #gROOT.SetBatch(kTRUE)
        for file in filenames :
            gROOT.ProcessLine('.x '+self.label_+"/"+file)
            list = []
            for point in range(grae.GetN()) : 
                if point == 0 : continue # GetPoint(0) is always empty!
                x = Double(0.)
                y = Double(0.)
                grae.GetPoint(point,x,y)
                xel = grae.GetErrorXlow(point)
                xeh = grae.GetErrorXhigh(point)
                yel = grae.GetErrorYlow(point)
                yeh = grae.GetErrorYhigh(point)
                list.append( (point-1,float(x),float(y),xel,xeh,yel,yeh) )
                dict[file] = list
        gROOT.SetBatch(kFALSE)
        return dict

    def patch(self,data) :
        for key,value in data.iteritems() :
            if data[key][0][1] == 225. :
                data[key][0] = (data[key][0][0],
                                237.5,
                                data[key][0][2],
                                data[key][0][3],
                                data[key][0][4],
                                data[key][0][5],
                                data[key][0][6])
            if data[key][2][1] == 375. :
                data[key][2] = (data[key][0][0],
                                350.,
                                data[key][0][2],
                                data[key][0][3],
                                data[key][0][4],
                                data[key][0][5],
                                data[key][0][6])
        return data

    def fill(self,input,index=None,systematics=False) :
        data = map(list,zip(*input))
        x = np.array( data[1] )
        if index is not None : x = x + self.offset_[index] # marker offsets
        y = np.array( data[2] )
        xeh = np.array( data[3] )
        xel = np.array( data[4] )
        yeh = np.array( data[5] )
        yel = np.array( data[6] )
        plot = TGraphAsymmErrors(len(x),x,y,xel,xeh,yel,yeh)
        return plot

    def plot(self) :
        if len(self.data_) == 0 :
            print "No data!"
            return

        setTDRStyle()
        
        leg = TLegend( 0.15, 0.91-(0.045*(len(self.data_)+1)), 0.35, 0.91 )
        leg.SetFillColor(0)
        leg.SetLineColor(0) 
        leg.SetShadowColor(0) 
        leg.SetTextSize(0.035)
        
        mg = TMultiGraph()

        plot = self.fill(self.systematics(),systematics=True)
        plot.SetLineColor(kGray)
        plot.SetFillColor(kGray)
        mg.Add(plot,"2")
        leg.AddEntry(plot,"Systematic uncertainty","f")

        for idx,key in enumerate(self.data_.keys()) :
            print self.titles_[key]
            plot = self.fill(self.data_[key],index=idx)
            plot.SetTitle("")
            plot.SetMarkerStyle(self.markers_[idx])
            plot.SetMarkerSize(self.size_[idx])
            plot.SetLineColor(1)
            plot.SetLineWidth(2)
            mg.Add(plot,"pZ")
            leg.AddEntry(plot,self.titles_[key],"p")

        canvas = TCanvas("tmp","tmp",900,600)
        mg.Draw("ap")
        mg.GetXaxis().SetTitle("H_{T} (GeV)")
        mg.GetYaxis().SetTitle("( N_{obs} - N_{pred} ) / N_{pred}")
        mg.GetYaxis().SetRangeUser(-2.,4.)
        mg.GetXaxis().SetRangeUser(self.bins_[0],self.bins_[self.nbins_])
        mg.GetXaxis().SetNdivisions(510)
        leg.Draw("same")
        prelim = "CMS Preliminary" if self.prelim_ else "CMS"
        str1 = "#splitline{%s}{L_{int} = %s fb^{-1}, #sqrt{s} = %s TeV}" % ( prelim, self.lumi_, self.energy_ )
        txt1 = TLatex(0.62,(0.91-0.06),str1)
        txt1.SetNDC()
        txt1.SetTextSize(0.04)
        txt1.Draw()

        str2 = ""
        temp = self.data_.keys()[0]
        if temp.count("_all") > 0 :
            str2 = "n_{jet} #geq 2"
        elif temp.count("_2") > 0 :
            str2 = "2 #leq n_{jet} #leq 3"
        elif temp.count("_3") > 0 :
            str2 = "n_{jet} #geq 4"
        else :
            str2 = "n_{jet} ERROR!"
        txt2 = TLatex(0.62,0.76,str2)
        txt2.SetNDC()
        txt2.SetTextSize(0.04)
        txt2.Draw()

        canvas.Update()
        raw_input("Press any key to continue...")
        canvas.SaveAs("temp.pdf")
        
    def systematics(self) :
#        for data in self.data_ :
#            total = 
#            for point in data :
#                if point[1] >= self.bins_[0] :
#                    total = total + 1./(total)
#            data = map(list,zip(*plot))
#            for point in data[]
        syst = [0.1,0.1,0.1,0.1,0.2,0.3]
        list = []
        for idx in range(len(self.regions_)) :
            if idx == len(self.regions_) - 1 : continue
            curr = self.regions_[idx]
            next = self.regions_[idx+1]
            xe = (self.bins_[next]-self.bins_[curr]) / 2.
            x = self.bins_[curr] + xe
            list.append( (idx,x,0.,xe,xe,syst[idx],syst[idx]) )
        return list
        
################################################################################
# CONFIGS ######################################################################
################################################################################

conf_8tev_20fb_broken = data(
    paths = [
    "https://clucas.web.cern.ch/clucas/RA1/Parked/Analysis/ControlLook_NewBinning_24Jun13/RA1_Documents_24_Jun/ClosureTests/"
    ],
    filenames = [
    "Btag_mu_to_mu_with_without_alphaT_all.C",
    "Btag_mu_zero_mu_one_no_alphaT_Cut_all.C",
    "Btag_mu_one_mu_greater_one_no_alphaT_Cut_all.C",
    "Btag_mu_to_dimuon_no_alphaT_Cut_all.C",
    "Btag_gamma_to_dimuon_all.C",
    ],
    )

conf_8tev_20fb_latest = data(
    paths=[
    "http://www.hep.ph.ic.ac.uk/~db1110/RA1_Documents_23_Jul/ClosureTests/"
    ],
    filenames=[
    "Btag_gamma_to_dimuon_no_alphaT_Cut_all.C",
    "Btag_mu_one_mu_two_no_alphaT_Cut_all.C",
    "Btag_mu_to_dimuon_no_alphaT_Cut_all.C",
    "Btag_mu_to_mu_with_without_alphaT_all.C",
    "Btag_mu_zero_mu_one_no_alphaT_Cut_all.C",
    ],
    )

conf_8tev_20fb_alphat = data(
    paths=[
    "http://www.hep.ph.ic.ac.uk/~db1110/RA1_Documents_23_Jul/ClosureTests/"
    ],
    filenames=[
    "Btag_mu_to_mu_with_without_alphaT_all.C",
    "Btag_mu_to_mu_with_without_alphaT_2.C",
    "Btag_mu_to_mu_with_without_alphaT_3.C",
    "Btag_dimuon_to_dimuon_with_without_alphaT_all.C",
    "Btag_dimuon_to_dimuon_with_without_alphaT_2.C",
    "Btag_dimuon_to_dimuon_with_without_alphaT_3.C",
    ],
    )

conf_8tev_20fb_mhtmet_2 = data(
    paths=[
    "http://www.hep.ph.ic.ac.uk/~db1110/RA1_Documents_ISR_PU_MHTMETSideband/ClosureTests_FullRange/"
    ],
    filenames=[
    "Btag_mu_to_mu_with_without_alphaT_2.C",
    "Btag_mu_zero_mu_one_no_alphaT_Cut_2.C",
    "Btag_mu_one_mu_two_no_alphaT_Cut_2.C",
    "Btag_mu_to_dimuon_no_alphaT_Cut_2.C",
    "Btag_gamma_to_dimuon_no_alphaT_Cut_2.C",
    ],
    )

conf_8tev_20fb_mhtmet_3 = data(
    paths=[
    "http://www.hep.ph.ic.ac.uk/~db1110/RA1_Documents_ISR_PU_MHTMETSideband/ClosureTests_FullRange/"
    ],
    filenames=[
    "Btag_mu_to_mu_with_without_alphaT_3.C",
    "Btag_mu_zero_mu_one_no_alphaT_Cut_3.C",
    "Btag_mu_one_mu_two_no_alphaT_Cut_3.C",
    "Btag_mu_to_dimuon_no_alphaT_Cut_3.C",
    "Btag_gamma_to_dimuon_no_alphaT_Cut_3.C",
    ],
    )

################################################################################
# EXECUTE ######################################################################
################################################################################

choice = 4

temp = None
if   choice == 0 : temp = conf_7tev_5fb_paper
elif choice == 1 : temp = conf_8tev_20fb_broken
elif choice == 2 : temp = conf_8tev_20fb_latest
elif choice == 3 : temp = conf_8tev_20fb_alphat
elif choice == 4 : temp = conf_8tev_20fb_mhtmet_2

temp.init()
temp.verbose()
temp.plot()

