import sys
#import ROOT as r
import numpy as np
import pickle

regions = {"signal":("Had","QCD"),"muon":("Muon","QCD_OneMuon"),}
lumis = {"signal":18.583,"bulk":19.294,"muon":19.255,}
processes = ["Data","SM","EWK","QCD",]
ht_bins = [200,275,325,375,475,575,675,775,875,None]
alphat_bins = [0.51,0.52,0.53,0.54,0.55,0.60,None]
mht_met_bins = [0.,1.25,None]

def get_file( dir=".",
              region="signal",
              process="SM" ) :
    file = dir+"/"+regions[region][0]+"_"+process+".root"
    f = r.TFile(file)
    if f.IsOpen() == False :
        print "File problem!",dir,region,process,file
        return None
    else :
        return f

def get_histo( file,
               region="signal",
               process="SM",
               ht_bin=(200,275),
               draw=False ) :
    name = "MHTovMET_vs_AlphaT_all"
    if file.IsOpen() == False :
        print "File problem!",file,region,process,ht_bin,draw,file
        return None
    bin = str(ht_bin[0])+"_"+str(ht_bin[1])
    if ht_bin[1] is not None : bin = str(ht_bin[0])
    histo = regions[region][1]+"_"+bin+"/"+name
    h = file.Get(histo)
    if type(h) != type(r.TH2D()) :
        #print "Histogram problem",file.GetName(),region,process,ht_bin,draw,histo
        return None
    c = None
    if draw :
        c = r.TCanvas("canvas","canvas")
        c.cd()
        c.SetLogz(1)
        h.Draw("COLZ")
    return h

def get_binning( histo, axis="X" ) :
    if histo is None : return None
    a = None
    if   axis == "X" : a = histo.GetXaxis()
    elif axis == "Y" : a = histo.GetYaxis()
    elif axis == "Z" : a = histo.GetZaxis()
    if a == None : return None
    if a.GetNbins() == 0 : return None
    bins = np.zeros((a.GetNbins()))
    for i in range(0,a.GetNbins()) :
        bins[i] = a.GetBinLowEdge(i+1)
    bins[i] = a.GetBinUpEdge(a.GetNbins())
    return bins

def get_contents( histo ) :
    if histo is None : return None
    xbins = get_binning(histo,"X")
    ybins = get_binning(histo,"Y")
    val = np.zeros( (histo.GetNbinsX()+2,histo.GetNbinsY()+2) )
    err = np.zeros( (histo.GetNbinsX()+2,histo.GetNbinsY()+2) )
    for i in range(0,histo.GetNbinsX()+1) :
        for j in range(0,histo.GetNbinsY()+1) :
            val[i][j] = histo.GetBinContent(i,j)
            err[i][j] = histo.GetBinError(i,j)
    return (xbins,ybins,val,err)

def dump_pickle() :
    data = {}
    for i in regions.keys() :
        for j in processes :
            for k in range(0,len(ht_bins)-1) :
                f = get_file( dir=".",
                              region=i,
                              process=j )
                h = get_histo( f,
                               region=i,
                               process=j,
                               ht_bin=(ht_bins[k],ht_bins[k+1]),
                               draw=False )
                if h is not None :
                    x = get_contents(h)
                    data[(i,j,k)] = x
                    print "Indices",i,j,k
                    print "X axis",len(x[0]),x[0]
                    print "Y axis",len(x[1]),x[1]
                    print "Values:",x[2].shape,x[2]
                    print "Errors",x[3].shape,x[3]
    pickle.dump(data,open("qcd_bkgd_est.pkl","w"))

def load_pickle() :
    return pickle.load(open("qcd_bkgd_est.pkl","r"))

def rebin(tuple) :
    def rebin_(old,new) :
        index = np.zeros((len(tuple[0])))
        for idx,val in enumerate(old) : 
            bin = None
            if val >= new[i] and new[i+1] is None : x = k
            elif ii >= xbins[i] and kk < xbins[i+1] : x = k
            
        return index
    xbins = rebin_(tuple[0],mht_met_bins)
    ybins = rebin_(tuple[1],alphat_bins)
    print xbins
    print ybins


#    val = np.zeros(len(mht_met_bins)-1,len(alphat_bins)-1)
#    err = np.zeros(len(mht_met_bins)-1,len(alphat_bins)-1)

def test(data) :
    for i in regions.keys() :
        for j in processes :
            if (i,j,0) in data.keys() : 
                x = data[(i,j,0)]
                rebin(x)
                print "Indices",i,j,0
                print "X axis",len(x[0]),x[0]
                print "Y axis",len(x[1]),x[1]
                print "Values:",x[2].shape,x[2]
                print "Errors",x[3].shape,x[3]
            else :
                print i,j,0

#dump_pickle()
data = load_pickle()
test(data)
#input("Please any key to continue...")

    

