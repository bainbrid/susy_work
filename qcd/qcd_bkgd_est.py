import sys
import ROOT as r
import numpy as np
import pickle
import math
from copy import copy
np.set_printoptions(precision=3)

################################################################################
# OPTIONS ######################################################################
################################################################################

pois_h = [1.15,1.36,1.53,1.73,1.98,2.21,2.42,2.61,2.80,3.00]
pois_l = [0.00,1.00,2.00,2.14,2.30,2.49,2.68,2.86,3.03,3.19]

regions = {"signal":("Had","QCD"),"muon":("Muon","QCD_OneMuon"),}
lumis = {"signal":18.583,"muon":19.255,"bulk":19.294,}
processes = ["Data","EWK","QCD",]

ht_bins = [200,275]#,325,375,475,575,675,775,875]
mht_met_bins = [0.,1.25]
alphat_bins = [0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.60]

muon_eff = (0.88,0.03)

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

################################################################################
# METHODS ######################################################################
################################################################################

def trigger_effs( path, verbose=False ) :
    data = {}
    for k in range(0,len(ht_bins)-1) :
        effs = np.zeros((len(mht_met_options),len(alpha_t_options)))
        errh = np.zeros((len(mht_met_options),len(alpha_t_options)))
        errl = np.zeros((len(mht_met_options),len(alpha_t_options)))
        for key,val in mht_met_options.iteritems():
            option = None
            if ht_bins[k] in ht_bin_options.keys() :
                name = "text_HT"+ht_bin_options[ht_bins[k]]+val[2]+"_AlphaT_ge2j.txt"
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
        data[k] = (mht_met_bins,alphat_bins,effs,errh,errl)
        if verbose == True :
            print "HT:",ht_bins[k]
            print "effs:"
            print view(effs)
            print "errh:"
            print view(errh)
            print "errl:"
            print view(errl)
            print "Fractions:"
            print view(errl/effs)
    return data

def trigger_hack( data, verbose=False ) :
    bin = 4
    if len(data) < bin+1 : return data
    tmp = data[bin]
    for k in range(0,len(ht_bins)-1) :
        if k > bin :
            data[k] = (data[k][0],data[k][1],tmp[2],tmp[3],tmp[4])
        if verbose == True :
            print "HT:",ht_bins[k]
            print "effs:"
            print view(data[k][2])
            print "errh:"
            print view(data[k][3])
            print "errl:"
            print view(data[k][4])
            print "Fractions:"
            print view(data[k][4]/data[k][2])
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
    h = file.Get(histo).Clone(histo)
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

def get_contents( stath, statl, systh, systl ) :
    if stath is None or statl is None or systh is None or systl is None : return None
    xbins = get_binning(stath,"X")
    ybins = get_binning(stath,"Y")
    val  = np.zeros( (stath.GetNbinsX(),stath.GetNbinsY()) )
    errh = np.zeros( (stath.GetNbinsX(),stath.GetNbinsY()) )
    errl = np.zeros( (stath.GetNbinsX(),stath.GetNbinsY()) )
    sysh = np.zeros( (stath.GetNbinsX(),stath.GetNbinsY()) )
    sysl = np.zeros( (stath.GetNbinsX(),stath.GetNbinsY()) )
    toth = np.zeros( (stath.GetNbinsX(),stath.GetNbinsY()) )
    totl = np.zeros( (stath.GetNbinsX(),stath.GetNbinsY()) )
    for i in range(stath.GetNbinsX()) :
        for j in range(stath.GetNbinsY()) :
            val[i][j] = stath.GetBinContent(i+1,j+1)
            errh[i][j] = stath.GetBinError(i+1,j+1)
            errl[i][j] = statl.GetBinError(i+1,j+1)
            sysh[i][j] = systh.GetBinError(i+1,j+1)
            sysl[i][j] = systl.GetBinError(i+1,j+1)
    toth = np.sqrt( errh*errh + sysh*sysh )
    totl = np.sqrt( errl*errl + sysl*sysl )
    return (xbins,ybins,val,errh,errl,sysh,sysl,toth,totl)

def debug(histo,ii,jj,kk) :
    if type(histo) != type(r.TH2D()) :
        return None
    for i in range(0,histo.GetNbinsX()+2) :
        for j in range(0,histo.GetNbinsY()+2) :
            histo.SetBinContent(i,j,1.) 
            histo.SetBinError(i,j,1.) 
    return histo

def poisson_errors(histo,errh) :
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
                    err = pois_h[index] if errh == True else pois_l[index]
                    histo.SetBinError(i,j,err*weight)
    return histo

def rebin( histo, xbins, ybins, xoverflow=False, yoverflow=False ) :
    if type(histo) != type(r.TH2D()) : return None
    x = copy(xbins)
    y = copy(ybins)
    x.append( histo.GetXaxis().GetXmax() )
    y.append( histo.GetYaxis().GetXmax() )
    h = r.TH2D(histo.GetName(),histo.GetTitle(),len(x)-1,np.array(x),len(y)-1,np.array(y))
    xaxis = histo.GetXaxis()
    yaxis = histo.GetYaxis()
    for xbin in range(len(x)-1) :
        for ybin in range(len(y)-1) :
            xindex = 0 if xoverflow is True and xbin == (len(x)-2) else 1
            yindex = 0 if yoverflow is True and ybin == (len(y)-2) else 1
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

def dump_pickle( effs=None, verbose=False ) :
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
                    #h = debug(h,i,j,k) 
                    stath = h
                    statl = h.Clone()
                    stath = rebin(stath,mht_met_bins,alphat_bins,True,True)
                    statl = rebin(statl,mht_met_bins,alphat_bins,True,True)
                    systh = stath.Clone()
                    systl = statl.Clone()
                    #stath = poisson_errors(stath,True,i,j,k)
                    #statl = poisson_errors(stath,False,i,j,k)
                    for ii in range(stath.GetNbinsX()) :
                        for jj in range(stath.GetNbinsY()) :
                            if i == "signal" and j == "Data" :
                                eff = 1.
                                sysh = 0.
                                sysl = 0.
                                if effs is not None :
                                    eff = effs[k][2][ii][jj]
                                    sysh = effs[k][3][ii][jj]
                                    sysl = effs[k][4][ii][jj]
                                val = stath.GetBinContent(ii+1,jj+1)/eff
                                errh = stath.GetBinError(ii+1,jj+1)/eff
                                errl = statl.GetBinError(ii+1,jj+1)/eff
                                stath.SetBinContent(ii+1,jj+1,val)
                                statl.SetBinContent(ii+1,jj+1,val)
                                stath.SetBinError(ii+1,jj+1,errh)
                                statl.SetBinError(ii+1,jj+1,errl)
                                systh.SetBinError(ii+1,jj+1,val*sysh)
                                systl.SetBinError(ii+1,jj+1,val*sysl)
                            elif i == "muon" and j == "Data" :
                                eff = muon_eff[0]
                                sysh = muon_eff[1]
                                sysl = muon_eff[1]
                                val = stath.GetBinContent(ii+1,jj+1)/eff
                                errh = stath.GetBinError(ii+1,jj+1)/eff
                                errl = statl.GetBinError(ii+1,jj+1)/eff
                                stath.SetBinContent(ii+1,jj+1,val)
                                statl.SetBinContent(ii+1,jj+1,val)
                                stath.SetBinError(ii+1,jj+1,errh)
                                statl.SetBinError(ii+1,jj+1,errl)
                                systh.SetBinError(ii+1,jj+1,val*sysh)
                                systl.SetBinError(ii+1,jj+1,val*sysl)
                            else :
                                systh.SetBinError(ii+1,jj+1,0.)
                                systl.SetBinError(ii+1,jj+1,0.)
                    x = get_contents(stath,stath,systh,systl)
                    data[(i,j,k)] = x
                    if verbose == True :
                        print "Region:",i,"Sample:",j,"HT:",ht_bins[k]
                        print "X axis binning:",x[0]
                        print "Y axis binning:",x[1]
                        print "Yields:"
                        print view(x[2])
                        print "Errh"
                        print view(x[3])
                        print "Errl"
                        print view(x[4])
                        print "Sysh"
                        print view(x[5])
                        print "Sysl"
                        print view(x[6])
                        print "Toth"
                        print view(x[7])
                        print "Totl"
                        print view(x[8])
    pickle.dump(data,open("qcd_bkgd_est.pkl","w"))

def load_pickle() :
    return pickle.load(open("qcd_bkgd_est.pkl","r"))

def yields( data, verbose=False ) :
    for i in regions.keys() :
        for j in processes :
            for k in range(0,len(ht_bins)-1) :
                if (i,j,k) in data.keys() : 
                    x = data[(i,j,k)]
                    if verbose == True :
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

def ewk_tf( data, verbose=False ) :
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
            if verbose == True :
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

def ewk_pred( data, ewk_tf, verbose=False ) :
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
            if verbose == True :
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

def qcd_only( data, ewk_pred, verbose=False ) :
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
            if verbose == True :
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
effs = trigger_hack(effs)

dump_pickle(effs,True)
#data = load_pickle()
#yields(data,True)
#tf = ewk_tf(data)
#ewk = ewk_pred(data,tf)
#qcd = qcd_only(data,ewk)

#input("Please any key to continue...")
