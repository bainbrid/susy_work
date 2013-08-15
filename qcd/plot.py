import numpy as np
from parse import *

dir = "./root_files_2"
region = "signal"
process = "QCD"
ht_bin = (200,275)

#f = get_file( filename="T2cc_175_165.root" )
f = get_file( dir, region, process )

if f is None : 
    print "Problem with file"
    quit()

#h = debug(h,i,j,k) 
histo = get_histo( f, region, process, ht_bin )

if histo is None :
    print "Problem with histo"
    quit()

mht_met_bins = [0.,1.25,]#list(get_binning(h,"X"))
#alphat_bins = [0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.60]
alphat_bins = [0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.60]
#his = h.Clone()
#his = rebin(his,mht_met_bins,alphat_bins)

#c1 = r.TCanvas("temp1","temp1")
#c1.cd()
#c1.SetLogz()
#his.GetYaxis().SetRangeUser(0.51,0.6)
#his.Draw("text colZ")
#
#c2 = r.TCanvas("temp2","temp2")
#c2.cd()
#prof2 = his.ProfileY()
#prof2.GetXaxis().SetTitle("alphaT")
#prof2.GetYaxis().SetTitle("Mean MHT/MET")
#prof2.Draw()

#c = r.TCanvas("projection","projection")
#c.cd()
#tmp1 = his.ProjectionX("tmp1",0,1)
##tmp2 = tmp1.Clone()
##tmp2.Scale(2.)
#tmp2 = his.ProjectionX("tmp2",1,2)
#print tmp1.GetBinContent(1)
#print tmp2.GetBinContent(1)
#gr1 =  r.TGraph(tmp1)
#gr2 =  r.TGraph(tmp2)
#mg = r.TMultiGraph()
#mg.Add(gr1,"l")
#mg.Add(gr2,"l")
#mg.Draw("al")
#c.Update()
#tmp1.Draw()
#tmp2.Draw("same")

def projection( h ) :
    c = r.TCanvas()
    c.cd()
    mg = r.TMultiGraph()
    l = r.TLegend(0.2,0.6,0.4,0.8) # 0.15, 0.91-0.045*(nhistos+(draw_weighted?2:1)), 0.35, 0.91 );
    gr = []
    for i in range( h.GetYaxis().GetNbins() ) :
        name = str(h.GetName()) + "_" + str(i)
        proj = []
        if i == 0 :
            tmp = h.ProjectionX(name)
            #tmp.Scale(1./tmp.Integral())
            #gr = r.TGraph(tmp)
            #gr.GetXaxis().SetTitle("MHT/MET")
            #gr.GetYaxis().SetTitle("a.u.")
            #gr.SetLineWidth(1*i)
            #mg.Add(gr,"l")
            tmp.SetLineWidth(1*(i+1))
            tmp.DrawNormalized("L")
            tmp.GetYaxis().SetRangeUser(0.,1.)
        else :
            tmp = h.ProjectionX(name,i,i+1)
            #tmp.Scale(1./float(tmp.Integral()))
            #gr = r.TGraph(tmp)
            #gr.SetLineWidth(1*i)
            #mg.Add(gr,"l")
            tmp.SetLineWidth(1*(i+1))
            tmp.DrawNormalized("L same")
        #l.AddEntry(gr)
    #mg.Draw("a")
    c.Update()
    input("Please any key to continue...")
#    l.SetFillColor(0)
#    l.SetLineColor(0) 
#    l.SetShadowColor(0)
#    l.SetTextSize(0.035)
#    l.Draw("same")

#projection(his)

#r.gStyle.SetOptStat(1111)
r.gStyle.SetOptFit(1111)

def plot_eff_vs_mhtmet( h ) :
    mhtmet = [0.,1.25,2.5,3.75]
    alphat = [0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.60]
    #mhtmet = [0.,0.25,0.5,0.75,1.,1.25,1.5,1.75,2.,2.25,2.5,2.75,3.,3.25,3.5,3.75,4.]
    #alphat = [0.51,0.55]
    his = rebin(h.Clone(),mhtmet,alphat)
    numer = his.ProjectionY("numer",1,2)
    denom = his.ProjectionY("denom",2,3)
    total = denom.Clone()
    #total.Add(numer)
    eff = r.TGraphAsymmErrors(numer,total,"pois")
    c = r.TCanvas()
    c.cd()
    eff.GetXaxis().SetRangeUser(0.51,0.6)
    eff.GetYaxis().SetRangeUser(0.,1.)
    eff.Draw("ap")
    c.Update()
    input("Please any key to continue...")
#    l.SetFillColor(0)
#    l.SetLineColor(0) 
#    l.SetShadowColor(0)
#    l.SetTextSize(0.035)
#    l.Draw("same")



#c3 = r.TCanvas("temp3","temp3")
#c3.cd()
#prof3 = his.ProfileX()
#prof3.GetXaxis().SetTitle("MHT/MET")
#prof3.GetYaxis().SetTitle("Mean alphaT")
#prof3.Draw()

plot_eff_vs_mhtmet(histo)

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
