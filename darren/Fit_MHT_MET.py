#!/usr/bin/env python
from ROOT import *
import ROOT as r
import logging,itertools
import os,fnmatch,sys

c1= r.TCanvas("Yields", "Yields",0,0,900,600)
c1.SetHighLightColor(2)
c1.SetFillColor(0)
c1.SetBorderMode(0)
c1.SetBorderSize(2)
c1.SetTickx(1)
c1.SetTicky(1)
c1.SetFrameBorderMode(0)
c1.SetFrameBorderMode(0)
c1.cd(1)
r.gStyle.SetOptFit(1)


#htbin = [150,200,275,325,375,475,575,675,775,875,975,1075,1175]

htbin = [178.2,235.2,297.5,347.5,416.4,517.3,618.4,716.9,819.9,919.,1019.,1289.,1300.]
# bands for PU ISR Only

"""
signal_list_zero = [0,0.8839,0.8559,0.8065,0.8105,0.7750,0.7842,0.6828,0.6565,0.6511,0.7215,0.6340]
signal_list_zero_err = [0,0.0042,0.0067,0.0088,0.0091,0.014,0.0219,0.0298,0.0415,0.0565,0.0805,0.062]

sideband_list_zero = [0,0.9031,0.8503,0.8861,0.8398,0.8303,0.7209,0.6861,0.7851,0.4644,0.603,0.4853]
sideband_list_zero_err = [0,0.0071,0.0114,0.0165,0.0171,0.0287,0.0419,0.0635,0.0986,0.1092,0.102,0.1122]

signal_list_two = [0,1.065,1.051,1.028,0.8865,0.9075,0.8814,0.851,0.6054,0.5871,0.4149,0.4102]
signal_list_two_err = [0,0.054,0.028,0.038,0.0344,0.0416,0.0584,0.085,0.1072,0.1509,0.2343,0.1968]

sideband_list_two = [0,1.085,1.011,0.8585,0.9708,0.7952,0.7945,0.7548,0.9791,1.221,0.8663,0.4283]
sideband_list_two_err = [0,0.025,0.027,0.0347,0.04,0.0527,0.0845,0.1261,0.2341,0.433,0.5875,0.3034]

dimuon_list_zero = [0,0.9383,0.980,0.9453,0.8899,0.7469,0.8890,0.5851,0.7341,0.7344,0.7033,0.5331]
dimuon_list_zero_err = [0,0.0130,0.02,0.0276,0.0272,0.0380,0.063,0.0732,0.1173,0.1739,0.2399,0.1683]

dimuon_sideband_list_zero = [0,0.8725,0.9001,0.8825,0.7921,0.7176,0.556,0.5608,0.7564,0.3604,0.5994,0.8981]
dimuon_sideband_list_zero_err = [0,0.0235,0.0389,0.0554,0.0503,0.0856,0.1182,0.1880,0.3332,0.2697,0.5994,0.5072]
"""


# bands for PU ISR Using Table yields with error propagation

#signal_list_two = [0.,1.0538,1.0539,1.0304,0.9071,0.8946,0.8935]#,0.9993,0.7555,0.6025,0.3731,0.4551]
#signal_list_two_err = [0.,0.05669,0.03078,0.04094,0.03739,0.0431,0.0614]#,0.09779,0.1211,0.1412,0.1545,0.1541]

#sideband_list_two = [0.,1.0112,1.0435,0.8779,1.0501,0.9052,0.8040]#,0.80211,1.2800,1.2487,0.9299,0.5116]
#sideband_list_two_err = [0.,0.0762,0.0405,0.0512,0.0575,0.0665,0.0914]#,0.1379,0.2779,0.492,0.4763,0.3708]

#signal_list_two =   [1.19,1.12,1.05,1.01,0.95,0.95,0.86]#,0.89,0.72,0.73,0.41,0.58  ]
#signal_list_two_err = [0.05,0.02,0.02,0.03,0.03,0.04,0.05]#,0.08,0.11,0.15,0.16,0.17]

#sideband_list_two = [1.08,1.11,1.04,0.89,0.99,0.88,0.82]#,0.86,1.08,1.75,1.04,0.92 ]
#sideband_list_two_err = [0.07,0.03,0.03,0.04,0.04,0.06,0.09]#,0.13,0.23,0.46,0.47,0.47]

#Gr4 cat only
fitnum = 12
fitlength = 12
signal_list_two = [0.0,1.05,1.05,1.03,0.91,0.95,0.90,1.00,0.77,0.66,0.43,0.53][:fitnum]
signal_list_two_err = [0.0,0.06,0.03,0.04,0.04,0.05,0.06,0.10,0.12,0.15,0.18,0.18][:fitnum]


sideband_list_two = [0.0,1.01,1.04,0.88,1.05,0.91,0.81,0.81,0.97,0.7,0.70,0.54,][:fitnum]
sideband_list_two_err = [0.0,0.08,0.04,0.05,0.06,0.07,0.09,0.14,0.27,0.48,0.47,0.42][:fitnum]

dimuon_list_zero = [0.96,0.95,0.98,0.96,0.90,0.75,0.89,0.69,0.77,0.65,0.82,0.48]
dimuon_list_zero_err = [0.03,0.02,0.02,0.03,0.03,0.04,0.06,0.08,0.12,0.15,0.23,0.14]

dimuon_sideband_list_zero = [0.92,0.89,0.93,0.93,0.78,0.77,0.72,0.93,0.70,0.46,0.28,1.22]
dimuon_sideband_list_zero_err = [0.05,0.03,0.04,0.06,0.05,0.09,0.13,0.22,0.25,0.27,0.38,0.43]

signal_list_zero = [0.92,0.89,0.85,0.81,0.81,0.78,0.79,0.69,0.68,0.65,0.75,0.66]
signal_list_zero_err = [0.01,0.01,0.01,0.01,0.01,0.02,0.02,0.03,0.04,0.06,0.09,0.07]

sideband_list_zero = [0.95,0.91,0.85,0.89,0.85,0.84,0.74,0.70,0.82,0.61,0.60,0.54]
sideband_list_zero_err = [0.02,0.01,0.02,0.02,0.02,0.03,0.04,0.06,0.10,0.11,0.14,0.11]

signal = r.TGraphErrors(12)
side = r.TGraphErrors(12)

signal_ttbar = r.TGraphErrors(12)
side_ttbar = r.TGraphErrors(12)

signal_dimuon = r.TGraphErrors(12)
side_dimuon = r.TGraphErrors(12)

signal.SetLineColor(4)
signal_ttbar.SetLineColor(4)
signal_dimuon.SetLineColor(4)

for n,entry in enumerate(signal_list_two):

     signal.SetPoint(n+1,htbin[n],signal_list_zero[n])
     signal.SetPointError(n+1,0,signal_list_zero_err[n])

     side.SetPoint(n+1,htbin[n],sideband_list_zero[n])
     side.SetPointError(n+1,0,sideband_list_zero_err[n])

     signal_ttbar.SetPoint(n+1,htbin[n],signal_list_two[n])
     signal_ttbar.SetPointError(n+1,0,signal_list_two_err[n])

     side_ttbar.SetPoint(n+1,htbin[n],sideband_list_two[n])
     side_ttbar.SetPointError(n+1,0,sideband_list_two_err[n])

     signal_dimuon.SetPoint(n+1,htbin[n],dimuon_list_zero[n])
     signal_dimuon.SetPointError(n+1,0,dimuon_list_zero_err[n])

     side_dimuon.SetPoint(n+1,htbin[n],dimuon_sideband_list_zero[n])
     side_dimuon.SetPointError(n+1,0,dimuon_sideband_list_zero_err[n])


signal.GetYaxis().SetRangeUser(0.0,1.9)
signal.GetXaxis().SetRangeUser(150,1500)

signal.Draw("AP")

fit_signal = r.TF1("fit","pol1", htbin[0], htbin[fitlength])
signal.Fit(fit_signal,"R")
fit_signal.SetLineColor(4)
fit_signal.Draw("SAME")

c1.SaveAs("./signal_zero_wjets.png")

side.GetYaxis().SetRangeUser(0.0,1.9)
side.GetXaxis().SetRangeUser(150,1500)

side.Draw("AP")

fit_side = r.TF1("fit","pol1", htbin[0], htbin[fitlength])
side.Fit(fit_side,"R")


fit_side.Draw("SAME")
signal.Draw("SAMEP")
fit_signal.Draw("SAME")

c1.SaveAs("./sideband_zero_wjets_both.png")

signal_dimuon.GetYaxis().SetRangeUser(0.0,1.9)
signal_dimuon.GetXaxis().SetRangeUser(150,1500)

signal_dimuon.Draw("AP")

fit_signal_dimuon = r.TF1("fit","pol1", htbin[0], htbin[fitlength])
signal_dimuon.Fit(fit_signal_dimuon,"R")
fit_signal_dimuon.SetLineColor(4)
fit_signal_dimuon.Draw("SAME")

c1.SaveAs("./signal_zero_dimuon.png")

side_dimuon.GetYaxis().SetRangeUser(0.0,1.9)
side_dimuon.GetXaxis().SetRangeUser(150,1500)

side_dimuon.Draw("AP")

fit_side_dimuon = r.TF1("fit","pol1", htbin[0], htbin[fitlength])
side_dimuon.Fit(fit_side_dimuon,"R")

signal_dimuon.Draw("SAMEP")
fit_signal_dimuon.Draw("SAME")
fit_side_dimuon.Draw("SAME")

c1.SaveAs("./sideband_zero_dimuon_both.png")

signal_ttbar.GetXaxis().SetRangeUser(150,1500)
signal_ttbar.GetYaxis().SetRangeUser(0.0,1.9)

signal_ttbar.Draw("AP")

fit_signal_ttbar = r.TF1("fit","pol1", htbin[0], htbin[fitlength])
signal_ttbar.Fit(fit_signal_ttbar,"R")
fit_signal_ttbar.SetLineColor(4)
fit_signal_ttbar.Draw("SAME")

c1.SaveAs("./signal_two_ttbar.png")


side_ttbar.GetXaxis().SetRangeUser(150,1500)
side_ttbar.GetYaxis().SetRangeUser(0.0,1.9)
side_ttbar.Draw("AP")

fit_side_ttbar = r.TF1("fit","pol1", htbin[0], htbin[fitlength])
side_ttbar.Fit(fit_side_ttbar,"R")


fit_side_ttbar.Draw("SAME")
signal_ttbar.Draw("SAMEP")
fit_signal_ttbar.Draw("SAME")

c1.SaveAs("./sideband_two_ttbar_both.png")

print htbin[fitnum]
