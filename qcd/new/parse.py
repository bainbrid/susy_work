import sys
import ROOT as r
import numpy as np
import pickle
import math
from copy import copy
np.set_printoptions(precision=3)
from uncertainties import unumpy
from utils import view

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
alphat_bins = [0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.60,]
               #0.61,0.62,0.63,0.64,0.65]#,0.66,0.67,0.68,0.69,0.70]

#mht_met_bins = [0.,5.]
#alphat_bins = [0.,10.]

################################################################################
# METHODS ######################################################################
################################################################################

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
#    temp = h.Integral()/h.GetEntries() if h.GetEntries() > 0 else 1.
#    print "val=%f, err=%f, weight=%f, val/weight=%f, h.Integral()=%f, h.GetEntries()=%f, temp=%f"\
#        %(val,err,weight,val/weight,h.Integral(),h.GetEntries(),temp)
    return (val,err,weight,val/weight)

def correlation_factor( h ) :
    if h is None : return (None,None,None)
    return h.GetCorrelationFactor()

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
    mhtmet = []
    alphat = []
    for i in range(h.GetNbinsX()) :
        mhtmet.append(h.GetXaxis().GetBinLowEdge(i+1))
    #mhtmet.append(h.GetXaxis().GetBinUpEdge(h.GetNbinsX()))
    for j in range(h.GetNbinsY()) :
        alphat.append(h.GetYaxis().GetBinLowEdge(j+1))
    #alphat.append(h.GetYaxis().GetBinUpEdge(h.GetNbinsY()))
    temp = unumpy.uarray(values,errors)
    return temp#(mhtmet,alphat,temp)

def debug(histo,ii,jj,kk) :
    if type(histo) != type(r.TH2D()) :
        return None
    for i in range(0,histo.GetNbinsX()+2) :
        for j in range(0,histo.GetNbinsY()+2) :
            value = 1.
            #value = j + i*(histo.GetNbinsY()+2) # 1.
            #value = i + j*(histo.GetNbinsX()+2) # 1.
            histo.SetBinContent(i,j,value) 
            histo.SetBinError(i,j,value) 
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

def rebin( histo, xbins, ybins, xoverflow=True, yoverflow=True ) :
    if type(histo) != type(r.TH2D()) : return None
    x = copy(xbins)
    y = copy(ybins)
    h = r.TH2D(histo.GetName()+"_copy1",histo.GetTitle()+"_copy1",len(x)-1,np.array(x),len(y)-1,np.array(y))
    for xbin in range(len(x)+1) :
        for ybin in range(len(y)+1) :

            xlow = xbin
            xhigh = histo.GetNbinsX()+1
            if xbin > 0 :
                xlow = histo.GetXaxis().FindBin( x[xbin-1] + 1.e-9 )
            if xbin < len(x) :
                xhigh = histo.GetXaxis().FindBin( x[xbin] + 1.e-9 ) - 1

            ylow = ybin
            yhigh = histo.GetNbinsY()+1
            if ybin > 0 :
                ylow = histo.GetYaxis().FindBin( y[ybin-1] + 1.e-9 )
            if ybin < len(y) :
                yhigh = histo.GetYaxis().FindBin( y[ybin] + 1.e-9 ) - 1

            error = r.Double(0.)
            integral = histo.IntegralAndError(xlow,xhigh,ylow,yhigh,error)

            h.SetBinContent(xbin,ybin,integral)
            h.SetBinError(xbin,ybin,error)

    if xoverflow : x.append( x[-1] + histo.GetBinWidth(histo.GetNbinsX()) )
    if yoverflow : y.append( y[-1] + histo.GetBinWidth(histo.GetNbinsY()) )
    hnew = r.TH2D(histo.GetName()+"_copy2",histo.GetTitle()+"_copy2",len(x)-1,np.array(x),len(y)-1,np.array(y))
    for xbin in range(len(x)+1) :
        for ybin in range(len(y)+1) :
            if xbin < h.GetNbinsX()+2 and ybin < h.GetNbinsY()+2 :
                hnew.SetBinContent( xbin, ybin, h.GetBinContent(xbin,ybin) )
                hnew.SetBinError( xbin, ybin, h.GetBinError(xbin,ybin) )
               
#    if xoverflow :
#        for ybin in range( 0, h.GetNbinsY()+2 ) :
#            val1 = h.GetBinContent( h.GetNbinsX(), ybin )
#            err1 = h.GetBinError( h.GetNbinsX(), ybin )
#            val2 = h.GetBinContent( h.GetNbinsX()+1, ybin )
#            err2 = h.GetBinError( h.GetNbinsX()+1, ybin )
#            h.SetBinContent( h.GetNbinsX(), ybin, val1 + val2 )
#            h.SetBinError( h.GetNbinsX(), ybin, math.sqrt( err1*err1 + err2*err2 ) )
#            h.SetBinContent( h.GetNbinsX()+1, ybin, 0. )
#            h.SetBinError( h.GetNbinsX()+1, ybin, 0. )

#    if yoverflow :
#        for xbin in range( 0, h.GetNbinsY()+2 ) :
#            val1 = h.GetBinContent( xbin, h.GetNbinsY() )
#            err1 = h.GetBinError( xbin, h.GetNbinsY() )
#            val2 = h.GetBinContent( xbin, h.GetNbinsY()+1 )
#            err2 = h.GetBinError( xbin, h.GetNbinsY()+1 )
#            h.SetBinContent( xbin, h.GetNbinsY(), val1 + val2 )
#            h.SetBinError( xbin, h.GetNbinsY(), math.sqrt( err1*err1 + err2*err2 ) )
#            h.SetBinContent( xbin, h.GetNbinsY()+1, 0. )
#            h.SetBinError( xbin, h.GetNbinsY()+1, 0. )

    #hnew.Draw()
    #input("")
    return hnew

def dump_pickle( verbose=False ) :
    correlation = {}
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
                    tuple = (i,j,ht_bins[k])
                    correlation[tuple] = correlation_factor(h)
                    #h.GetXaxis().SetRangeUser(0.,2.50)
                    h.GetYaxis().SetRangeUser(0.51,0.70)
                    corr = correlation_factor(h)
#                    if tuple == ('signal','QCD',200) :
#                        h.Draw()
#                        input()
                    #h = debug(h,i,j,k) 
                    h = rebin(h,mht_met_bins,alphat_bins,True,True)
                    weights[tuple] = weight(h)
                    data[tuple] = get_contents(h)
                    if verbose == True : 
                        #print tuple, weights[tuple][2]
                        print tuple, correlation[tuple],corr
                        # summary(data[tuple])
    temp = (weights,data)
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
#print view( data[("signal","QCD",275)] )
#print data[("signal","EWK",200)] 
#print view( data[("signal","EWK",200)] )
