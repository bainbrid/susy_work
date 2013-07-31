import sys
import ROOT as r
import numpy as np
import pickle
from copy import copy
np.set_printoptions(precision=3)

################################################################################
# OPTIONS ######################################################################
################################################################################

pois_h = [1.15,1.36,1.53,1.73,1.98,2.21,2.42,2.61,2.80,3.00]
pois_l = [0.00,1.00,2.00,2.14,2.30,2.49,2.68,2.86,3.03,3.19]
regions = {"signal":("Had","QCD"),"muon":("Muon","QCD_OneMuon"),}
lumis = {"signal":18.583,"bulk":19.294,"muon":19.255,}
processes = ["Data","SM","EWK","QCD",]
ht_bins = [200,275,325,375,475,575,675,775,875,None]
mht_met_bins = [0.,1.25,None]
pre_alphat_bins = [0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.60,None]
alphat_bins = [0.51,0.52,0.53,0.54,0.55,0.60,None]
signal_eff = [0.012,0.013,0.032,0.053,0.123,0.310,0.512,0.679,0.795,0.795]
muon_eff = 0.88

################################################################################
# METHODS ######################################################################
################################################################################

def trigger_effs( path ) :

    mht_met_options = {"with":(0,0.,"_mht"),
                       "without":(1,1.25,""),}
    ht_bin_options = {200:"200_275_73_73_36",
                      275:"275_325_73_73_36",
                      325:"325_375_86_86_43",
                      375:"375_475_100_100_50",
                      475:"475_575_100_100_50",
                      575:"575_675_100_100_50",
                      675:"675_775_100_100_50",
                      775:"775_875_100_100_50",
                      875:"875_100_100_50",}
    alpha_t_options = {"0.515000":(0,0.51),
                       "0.525000":(1,0.52),
                       "0.535000":(2,0.53),
                       "0.545000":(3,0.54),
                       "0.555000":(4,0.55),
                       "0.565000":(5,0.56),
                       "0.575000":(6,0.57),
                       "0.585000":(7,0.58),
                       "0.595000":(8,0.59),
                       "0.645000":(9,0.60),}

    data = {}
    for k in range(0,len(ht_bins)-1) :
        effs = np.zeros((len(mht_met_options),len(alpha_t_options)))
        errh = np.zeros((len(mht_met_options),len(alpha_t_options)))
        errl = np.zeros((len(mht_met_options),len(alpha_t_options)))
        for key,val in mht_met_options.iteritems():
            option = None
            if ht_bins[k] in ht_bin_options.keys() :
                name = "text_HT"+ht_bin_options[ht_bins[k]]+val[2]+"_AlphaT_ge2j.txt"
                print name
            else : continue
            file = open(path+name)
            for line in file.readlines() :
                entries = line.split()
                if entries[2] in alpha_t_options.keys() :
                    alphat_bin = alpha_t_options[entries[2]][0]
                    mht_met_bin = val[0]
                    effs[mht_met_bin][alphat_bin] = entries[3]
                    errh[mht_met_bin][alphat_bin] = entries[5]
                    errl[mht_met_bin][alphat_bin] = entries[7]
            file.close()

        mht_met_bins = [ val[1] for key,val in mht_met_options.iteritems() ]
        alpha_t_bins = [ val[1] for key,val in alpha_t_options.iteritems() ]
        data[k] = (mht_met_bins,alpha_t_bins,effs,errh,errl)
        if True :
            print "HT:",ht_bins[k]
            print "effs:"
            print view(effs)
            print "errh:"
            print view(errh)
            print "errl:"
            print view(errl)
        else :
            print "Problem:",k

    return data

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
    if ht_bin[1] is None : bin = str(ht_bin[0])
    histo = regions[region][1]+"_"+bin+"/"+name
    h = file.Get(histo)
    if type(h) != type(r.TH2D()) :
        #print "Histogram problem",file.GetName(),region,process,ht_bin,draw,histo
        return None
    c = None
    if draw is True :
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
    if a is None : return None
    if a.GetNbins() == 0 : return None
    bins = np.zeros((a.GetNbins()+1))
    for i in range(a.GetNbins()) :
        bins[i] = a.GetBinLowEdge(i+1)
    bins[-1] = a.GetBinUpEdge(a.GetNbins())
    return bins

def get_contents( histo ) :
    if histo is None : return None
    xbins = get_binning(histo,"X")
    ybins = get_binning(histo,"Y")
    val = np.zeros( (histo.GetNbinsX(),histo.GetNbinsY()) )
    errh = np.zeros( (histo.GetNbinsX(),histo.GetNbinsY()) )
    errl = np.zeros( (histo.GetNbinsX(),histo.GetNbinsY()) )
    for i in range(histo.GetNbinsX()) :
        for j in range(histo.GetNbinsY()) :
            val[i][j] = histo.GetBinContent(i+1,j+1)
            errh[i][j] = histo.GetBinError(i+1,j+1)
            np.copyto(errh,errl)
    return (xbins,ybins,val,errh,errl)

def preprocess(histo,ii,jj,kk) :
    if type(histo) != type(r.TH2D()) :
        return None
    for i in range(0,histo.GetNbinsX()+2) :
        for j in range(0,histo.GetNbinsY()+2) :
            val = histo.GetBinContent(i,j)
            err = histo.GetBinError(i,j)
            if err > 0. :
                entries = (val*val)/(err*err)
                weight = val/entries
                error = err
                index = int(entries)
                if index < 10 :
                    if pois_h[index] > pois_l[index] : # symmetric errors!
                        error = pois_h[index]*weight
                    else :
                        error = pois_l[index]*weight
                    histo.SetBinError(i,j,error)
            err = histo.GetBinError(i,j)
#            @@
            if False : # debug only
                histo.SetBinContent(i,j,1.) 
                histo.SetBinError(i,j,1.) 
    return histo

def get_bins(histo,xbins,ybins) :
    x = copy(xbins)
    y = copy(ybins)
    if x[0] is None : x[0] = histo.GetXaxis().GetXmin()
    if y[0] is None : y[0] = histo.GetYaxis().GetXmin()
    if x[-1] is None : x[-1] = histo.GetXaxis().GetXmax()
    if y[-1] is None : y[-1] = histo.GetYaxis().GetXmax()
    return x,y

def rebin(histo,xbins,ybins,xoverflow=True,yoverflow=True) :
    if type(histo) != type(r.TH2D()) : return None
    x,y = get_bins(histo,xbins,ybins)
    h = r.TH2D(histo.GetName(),histo.GetTitle(),len(x)-1,np.array(x),len(y)-1,np.array(y))
    xaxis = histo.GetXaxis()
    yaxis = histo.GetYaxis()
    for xbin in range(len(x)-1) :
        for ybin in range(len(y)-1) :
            xindex = 0 if xoverflow==True and xbin==(len(x)-2) else 1
            yindex = 0 if yoverflow==True and ybin==(len(y)-2) else 1
            error = r.Double(0.)
            integral = histo.IntegralAndError(xaxis.FindBin(x[xbin]),
                                              xaxis.FindBin(x[xbin+1])-xindex,
                                              yaxis.FindBin(y[ybin]),
                                              yaxis.FindBin(y[ybin+1])-yindex,
                                              error)
            h.SetBinContent(xbin+1,ybin+1,integral)
            h.SetBinError(xbin+1,ybin+1,error)
    return h

def view(arr) :
    return np.flipud(copy(arr).transpose())

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
                    h = rebin(h,mht_met_bins,pre_alphat_bins)
                    h = preprocess(h,i,j,k)
                    h = rebin(h,mht_met_bins,alphat_bins)
                    x = get_contents(h)
                    data[(i,j,k)] = x
                    if False :
                        print "Region:",i,"Sample:",j,"HT:",ht_bins[k]
                        print "X axis binning:",x[0]
                        print "Y axis binning:",x[1]
                        print "Yields:"
                        print view(x[2])
                        print "Errors:"
                        print view(x[3])
    pickle.dump(data,open("qcd_bkgd_est.pkl","w"))

def load_pickle() :
    return pickle.load(open("qcd_bkgd_est.pkl","r"))

def yields(data) :
    for i in regions.keys() :
        for j in processes :
            for k in range(0,len(ht_bins)-1) :
                if (i,j,k) in data.keys() : 
                    x = data[(i,j,k)]
                    if True :
                        print "Region:",i,"Sample:",j,"HT:",ht_bins[k]
                        print "X axis binning:",x[0]
                        print "Y axis binning:",x[1]
                        print "Yields:"
                        print view(x[2])
                        print "Errh:"
                        print view(x[3])
                        print "Errl:"
                        print view(x[4])
                        print "Fractions:"
                        print view(x[3]/x[2])
                else :
                    print "Problem:",i,j,k

def ewk_tf(data) :
    tf = {}
    for k in range(0,len(ht_bins)-1) :
        signal = ("signal","EWK",k)
        muon = ("muon","EWK",k)
        if signal in data.keys() and muon in data.keys() :
            xbins = data[signal][0]
            ybins = data[signal][1]
            val = data[signal][2]/data[muon][2]
            errh_sig = data[signal][3]/data[signal][2]
            errl_sig = data[signal][4]/data[signal][2]
            errh_mu = data[muon][3]/data[muon][2]
            errl_mu = data[muon][4]/data[muon][2]
            errh = np.sqrt(errh_sig*errh_sig+errh_mu*errh_mu)
            errl = np.sqrt(errl_sig*errl_sig+errl_mu*errl_mu)
            tf[k] = (xbins,ybins,val,errh,errl)
            if True :
                print "HT:",ht_bins[k]
                print "X axis binning:",xbins
                print "Y axis binning:",ybins
                print "Yields:"
                print view(val)
                print "Errh:"
                print view(errh)
                print "Errl:"
                print view(errl)
                print "Fractions:"
                print view(errh/val)
            else :
                print "Problem:","EWK",k
    return tf

def ewk_pred(data,ewk_tf) :
    pred = {}
    for k in range(0,len(ht_bins)-1) :
        obs = ("muon","Data",k)
        if obs in data.keys() and k in ewk_tf.keys() :
            xbins = data[obs][0]
            ybins = data[obs][1]
            val = data[obs][2] * ewk_tf[k][2]
            errh_obs = data[obs][3]/data[obs][2]
            errl_obs = data[obs][4]/data[obs][2]
            errh_tf = ewk_tf[k][3]/ewk_tf[k][2]
            errl_tf = ewk_tf[k][4]/ewk_tf[k][2]
            errh = np.sqrt(errh_obs*errh_obs+errh_tf*errh_tf)
            errl = np.sqrt(errl_obs*errl_obs+errl_tf*errl_tf)
            pred[k] = (xbins,ybins,val,errh,errl)
            if True :
                print "HT:",ht_bins[k]
                print "X axis binning:",xbins
                print "Y axis binning:",ybins
                print "Yields:"
                print view(val)
                print "Errh:"
                print view(errh)
                print "Errl:"
                print view(errl)
                print "Fractions:"
                print view(errh/val)
            else :
                print "Problem:","EWK",k
    return pred

def qcd_only(data,ewk_pred) :
    qcd = {}
    for k in range(0,len(ht_bins)-1) :
        obs = ("signal","Data",k)
        if obs in data.keys() and k in ewk_pred.keys() :
            xbins = data[obs][0]
            ybins = data[obs][1]
            val = data[obs][2] - ewk_pred[k][2]
            errh_obs = data[obs][3]
            errl_obs = data[obs][4]
            errh_pred = ewk_pred[k][3]
            errl_pred = ewk_pred[k][4]
            errh = np.sqrt(errh_obs*errh_obs+errh_pred*errh_pred)
            errl = np.sqrt(errl_obs*errl_obs+errl_pred*errl_pred)
            qcd[k] = (xbins,ybins,val,errh,errl)
            if True :
                print "HT:",ht_bins[k]
                print "X axis binning:",xbins
                print "Y axis binning:",ybins
                print "Yields:"
                print view(val)
                print "Errh:"
                print view(errh)
                print "Errl:"
                print view(errl)
                print "Fractions:"
                print view(errh/val)
            else :
                print "Problem:","EWK",k
    return qcd

################################################################################
# EXECUTE ######################################################################
################################################################################

effs = trigger_effs("./files/")

#dump_pickle()
#data = load_pickle()
#yields(data)
#tf = ewk_tf(data)
#ewk = ewk_pred(data,tf)
#qcd = qcd_only(data,ewk)

#input("Please any key to continue...")
