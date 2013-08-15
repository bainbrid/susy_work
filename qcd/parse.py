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

# (filename prefix, dir prefix)
regions = {"signal":("Had","QCD"),
           "muon":("Muon","QCD_OneMuon"),}

processes = ["Data","EWK","QCD","SM","Signal","SM_Signal"]

directory = "./root_files_2"

ht_bins = [200,275]#,325,375,475,575,675,775,875]

mht_met_bins = [0.,1.25,2.50,3.75,5.00]

alphat_bins = [0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.60,0.65,0.70]

overflow = True

################################################################################
# METHODS ######################################################################
################################################################################

def view(arr) :
    return np.flipud(copy(arr).transpose())

def get_file( dir=".",
              region="signal",
              process="SM",
              filename=None) :
    f = None
    if filename is not None :
        f = r.TFile(filename)
    else : 
        file = dir+"/"+regions[region][0]+"_"+process+".root"
        f = r.TFile(file)
    if f.IsOpen() == False :
        print "File problem!",dir,region,process,file
        return None
    else :
        return f

def get_histo( file,
               region="signal",
               process="Data",
               ht_bin=(200,275) ) :
    name = "MHTovMET_vs_AlphaT_all"
    if file.IsOpen() == False :
        print "File problem!",file,region,process,ht_bin,file
        return None
    bin = str(ht_bin[0])+"_"+str(ht_bin[1])
    if ht_bin[1] is None : bin = str(ht_bin[0])
    histo = regions[region][1]+"_"+bin+"/"+name
    h = file.Get(histo).Clone(file.GetName()+"/"+histo)
    if type(h) != type(r.TH2D()) :
        #print "Histogram problem",file.GetName(),region,process,ht_bin,histo
        return None
    return h

#def plot_histo( histo,
#                dir=None,
#                region=None,
#                process=None,
#                ht_bin=None ) :
#
#    file = dir+"/"+regions[region][0]+"_"+process+".root"
#    name = "MHTovMET_vs_AlphaT_all"
#    bin = str(ht_bin[0])+"_"+str(ht_bin[1])
#    if ht_bin[1] is None : bin = str(ht_bin[0])
#    his = regions[region][1]+"_"+bin+"/"+name
#
#    c = r.TCanvas(his,his)
#
#    if a is None : return None
#    if a.GetNbins() == 0 : return None
#    bins = np.zeros((a.GetNbins()+1))
#    for i in range(a.GetNbins()) :
#        bins[i] = a.GetBinLowEdge(i+1)
#    bins[-1] = a.GetBinUpEdge(a.GetNbins())
#
#    if draw is True :
#        c = 
#        c.cd()
#        c.SetLogz(1)
#        h.Draw("COLZ")
#    return c

def get_binning( histo, axis="X" ) :
    if histo is None : return None
    a = None
    if   axis == "X" : a = histo.GetXaxis()
    elif axis == "Y" : a = histo.GetYaxis()
    elif axis == "Z" : a = histo.GetZaxis()
    if a is None : return None
    if a.GetNbins() == 0 : return None
    bins = np.zeros((a.GetNbins()))#+1))
    for i in range(a.GetNbins()) :
        bins[i] = a.GetBinLowEdge(i+1)
    #bins[-1] = a.GetBinUpEdge(a.GetNbins())
    return bins

def get_contents( stath, statl, systh, systl ) :
    if stath is None or statl is None or systh is None or systl is None : return None
    xbins = get_binning(stath,"X")
    ybins = get_binning(stath,"Y")
    val  = np.zeros( (stath.GetNbinsX(),stath.GetNbinsY()) )
    errh = np.zeros_like( val )
    errl = np.zeros_like( val )
    sysh = np.zeros_like( val )
    sysl = np.zeros_like( val )
    toth = np.zeros_like( val )
    totl = np.zeros_like( val )
    for i in range(stath.GetNbinsX()) :
        for j in range(stath.GetNbinsY()) :
            val[i][j] = stath.GetBinContent(i+1,j+1)
            errh[i][j] = stath.GetBinError(i+1,j+1)
            errl[i][j] = statl.GetBinError(i+1,j+1)
#            sysh[i][j] = systh.GetBinError(i+1,j+1)
#            sysl[i][j] = systl.GetBinError(i+1,j+1)
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

def dump_pickle( verbose=False ) :
    data = {}
    for i in regions.keys() :
        for j in processes :
            for k in range(0,len(ht_bins)-1) :
                f = get_file( dir=directory,
                              region=i,
                              process=j )
                if f is None : continue
                h = get_histo( f,
                               region=i,
                               process=j,
                               ht_bin=(ht_bins[k],ht_bins[k+1]) )
                if h is not None :
                    #h = debug(h,i,j,k) 
                    stath = h
                    statl = h.Clone()
                    stath = rebin(stath,mht_met_bins,alphat_bins,overflow,overflow)
                    statl = rebin(statl,mht_met_bins,alphat_bins,overflow,overflow)
                    systh = stath.Clone()
                    systl = statl.Clone()
#                    c = plot_histo( stath,
#                                    dir=directory,
#                                    region=i,
#                                    process=j,
#                                    ht_bin=(ht_bins[k],ht_bins[k+1]) )
#                    for ii in range(stath.GetNbinsX()) :
#                        for jj in range(stath.GetNbinsY()) :
#                            if i == "signal" and j == "Data" :
#
#                                eff = 1.
#                                eff_stath = 0.
#                                eff_statl = 0.
#                                eff_systh = 0.
#                                eff_systl = 0.
#                                if effs is not None :
#                                    eff = effs[k][2][ii][jj]
#                                    eff_stath = effs[k][3][ii][jj]
#                                    eff_statl = effs[k][4][ii][jj]
#                                    eff_systh = effs[k][5][ii][jj]
#                                    eff_systl = effs[k][6][ii][jj]
#
#                                val = stath.GetBinContent(ii+1,jj+1)/eff
#                                val_stath = pois_h[0]
#                                val_statl = pois_l[0]
#                                if stath.GetBinContent(ii+1,jj+1) > 0. :
#                                    val_stath = stath.GetBinError(ii+1,jj+1)/stath.GetBinContent(ii+1,jj+1)
#                                    val_statl = statl.GetBinError(ii+1,jj+1)/stath.GetBinContent(ii+1,jj+1)
#                                stath.SetBinContent(ii+1,jj+1,val)
#                                statl.SetBinContent(ii+1,jj+1,val)
#                                stath.SetBinError(ii+1,jj+1, math.sqrt(val_stath*val_stath+eff_stath*eff_stath) )
#                                statl.SetBinError(ii+1,jj+1, math.sqrt(val_statl*val_statl+eff_statl*eff_statl) )
#                                systh.SetBinError(ii+1,jj+1,eff_systh)
#                                systl.SetBinError(ii+1,jj+1,eff_systl)
#
#                            elif i == "muon" and j == "Data" :
#                                eff = muon_eff[0]
#                                eff_stath = muon_eff[1]
#                                eff_statl = muon_eff[1]
#                                eff_systh = muon_eff[2]
#                                eff_systl = muon_eff[2]
#                                val = stath.GetBinContent(ii+1,jj+1)/eff
#                                val_stath = pois_h[0]
#                                val_statl = pois_l[0]
#                                if stath.GetBinContent(ii+1,jj+1) > 0. :
#                                    val_stath = stath.GetBinError(ii+1,jj+1)/stath.GetBinContent(ii+1,jj+1)
#                                    val_statl = statl.GetBinError(ii+1,jj+1)/stath.GetBinContent(ii+1,jj+1)
#                                stath.SetBinContent(ii+1,jj+1,val)
#                                statl.SetBinContent(ii+1,jj+1,val)
#                                stath.SetBinError(ii+1,jj+1, math.sqrt(val_stath*val_stath+eff_stath*eff_stath) )
#                                statl.SetBinError(ii+1,jj+1, math.sqrt(val_statl*val_statl+eff_statl*eff_statl) )
#                                systh.SetBinError(ii+1,jj+1,eff_systh)
#                                systl.SetBinError(ii+1,jj+1,eff_systl)
#
#                            elif j != "Data" :
#                                lumi = regions[i][2]
#                                lumi_stath = 0.
#                                lumi_statl = 0.
#                                lumi_systh = regions[i][3]
#                                lumi_systl = regions[i][3]
#                                val = stath.GetBinContent(ii+1,jj+1)*lumi
#                                val_stath = pois_h[0]
#                                val_statl = pois_l[0]
#                                if stath.GetBinContent(ii+1,jj+1) > 0. :
#                                    val_stath = stath.GetBinError(ii+1,jj+1)/stath.GetBinContent(ii+1,jj+1)
#                                    val_statl = statl.GetBinError(ii+1,jj+1)/stath.GetBinContent(ii+1,jj+1)
#                                stath.SetBinContent(ii+1,jj+1,val)
#                                statl.SetBinContent(ii+1,jj+1,val)
#                                stath.SetBinError(ii+1,jj+1, math.sqrt(val_stath*val_stath+lumi_stath*lumi_stath) )
#                                statl.SetBinError(ii+1,jj+1, math.sqrt(val_statl*val_statl+lumi_statl*lumi_statl) )
#                                systh.SetBinError(ii+1,jj+1,lumi_systh)
#                                systl.SetBinError(ii+1,jj+1,lumi_systl)
#
#                            else :
#                                errh = pois_h[0]
#                                errl = pois_l[0]
#                                if stath.GetBinContent(ii+1,jj+1) > 0. :
#                                    errh = stath.GetBinError(ii+1,jj+1)/stath.GetBinContent(ii+1,jj+1)
#                                    errl = statl.GetBinError(ii+1,jj+1)/stath.GetBinContent(ii+1,jj+1)
#                                stath.SetBinError(ii+1,jj+1,errh)
#                                statl.SetBinError(ii+1,jj+1,errl)
#                                systh.SetBinError(ii+1,jj+1,0.)
#                                systl.SetBinError(ii+1,jj+1,0.)
                    x = get_contents(stath,stath,systh,systl)
                    tuple = (i,j,k)
                    data[tuple] = x
                    if verbose == True :
                        print tuple
                        #summary(data[tuple])
    pickle.dump(data,open("data.pkl","w"))
    return data

################################################################################
# EXECUTE ######################################################################
################################################################################

data = dump_pickle(True)

#input("Please any key to continue...")

#print data[("signal","Data",0)][0]
#print data[("signal","Data",0)][1]
#print view( data[("signal","Data",0)][2] )
