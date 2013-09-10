import sys
import ROOT as r
import numpy as np
import pickle
import math
from copy import copy
np.set_printoptions(precision=3)
from uncertainties import unumpy

################################################################################
# OPTIONS ######################################################################
################################################################################

pois_h = [1.15,1.36,1.53,1.73,1.98,2.21,2.42,2.61,2.80,3.00]
pois_l = [0.00,1.00,2.00,2.14,2.30,2.49,2.68,2.86,3.03,3.19]

# (filename prefix, dir prefix)
regions = {"signal":("Had","QCD"),
           "muon":("Muon","QCD_OneMuon"),}

processes = ["Data","EWK","QCD","SM","Signal"]

directory = "../root_files_5"

#ht_bins = [275,325,None]
ht_bins = [150,200,275,325,375,475,575,675,775,875,975,1075,None]

mht_met_bins = [0.,1.25,2.50,3.75,5.00]

#alphat_bins = [0.51,0.55,0.60,0.65,0.70]
alphat_bins = [0.50,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,
               0.60,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,
               0.70]

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
    #print file.GetName()+"/"+histo
    h = file.Get(histo)#.Clone(file.GetName()+"/"+histo)
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

def weight( h ) :
    if h is None : return (None,None,None)
    err = r.Double(0.)
    val = h.IntegralAndError(1,h.GetNbinsX(),1,h.GetNbinsY(),err)
    weight = err*err/val if val > 0. else 1.
    #temp = h.Integral()/h.GetEntries() if h.GetEntries() > 0 else 1.
    #print val,err,weight,val/weight,h.Integral(),h.GetEntries(),temp
    return (val,err,weight,val/weight)

def get_contents( h ) :
    if h is None : return None
    values = np.zeros((h.GetNbinsX(),h.GetNbinsY()))
    errors = np.zeros((h.GetNbinsX(),h.GetNbinsY()))
    tuple = weight(h) # val,err,wei
    for i in range(h.GetNbinsX()) :
        for j in range(h.GetNbinsY()) :
            val = h.GetBinContent(i+1,j+1)
            err = h.GetBinError(i+1,j+1) 
            values[i][j] = val
            errors[i][j] = err
            wei = err*err/val if val > 0. else tuple[2]
            entries = int(abs(val/wei if wei > 0. else val))
            if entries < 10 : 
                max_err = pois_h[entries] 
                if pois_l[entries] > pois_h[entries] : 
                    max_err = pois_l[entries]
                errors[i][j] = max_err*wei
#                print "i=%i, j=%i, val=%g, err=%g, wei=%g, (val,err,wei,num)=(%g,%g,%g,%g), entries=%i, pois=%g, (val,err)=(%g,%g)"\
#                    %(i,j,val,err,wei,tuple[0],tuple[1],tuple[2],tuple[3],entries,max_err,val,max_err*wei)
    temp = unumpy.uarray(values,errors)
    return temp

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
    weights = {}
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
                    # h = debug(h,i,j,k) 
#                    h = rebin(h,mht_met_bins,alphat_bins,overflow,overflow)
                    tuple = (i,j,ht_bins[k])
                    weights[tuple] = weight(h)
                    data[tuple] = get_contents(h)
                    if verbose == True : 
                        print tuple#, weights[tuple][0], weights[tuple][1], weights[tuple][2]
                        # summary(data[tuple])
    pickle.dump(data,open("data.pkl","w"))
    #print [ x[1][2] for x in weights.iteritems() ]
    return data

################################################################################
# EXECUTE ######################################################################
################################################################################

data = dump_pickle(True)

#input("Please any key to continue...")

#print data[("signal","Data",0)][0]
#print data[("signal","Data",0)][1]
#print view( data[("signal","Data",0)][2] )
