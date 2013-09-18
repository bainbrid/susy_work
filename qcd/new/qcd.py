import sys
import numpy as np
import pickle
import math
import utils
import trigger
from copy import copy
from uncertainties import ufloat
from uncertainties import UFloat
from uncertainties import unumpy

################################################################################
# CONFIG #######################################################################
################################################################################

# (filename prefix, dir prefix, lumi, % uncert)
regions = {"signal":("Had","QCD",18.583*10.,0.04),
           "muon":("Muon","QCD_OneMuon",19.255*10.,0.04),}

processes = ["Data","EWK","QCD","SM","Signal"]

# (eff, stat uncert, syst uncert)
muon_eff = (0.88,0.00,0.03)

# (syst uncert, debug mode, debug eff, debug syst uncert)
trigger_syst = (0.03,False,1.,0.)

pickle_file = "data.pkl"

ht_bins = [200,275,375,475,575,675,775,875,975,1075]
mhtmet_bins = [0.,1.25,2.50,3.75,5.00]
#alphat_bins = [0.,10.]
alphat_bins = [0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.60,]
               #0.61,0.62,0.63,0.64,0.65]#,0.66,0.67,0.68,0.69,0.70]

index = 0
ht_bin = ht_bins[index]

add_signal = False

use_data = True

qcd_debug = False

################################################################################
# UTILITY METHODS ##############################################################
################################################################################

def load_pickle( verbose=False ) :
    data = pickle.load(open(pickle_file,"r"))
    if verbose == True :
        for i in regions.keys() :
            for j in processes :
                for k in range(0,len(ht_bins)) :
                    if (i,j,ht_bins[k]) in data.keys() : 
                        x = data[(i,j,ht_bins[k])]
                        print "Region= '%s', sample= '%s', HT= %i "%(i,j,ht_bins[k])
                    else : print "Problem:",i,j,ht_bins[k]
    return data

def ewk_tf( signal, control ) :
    return utils.safe_div(signal,control)

def ewk_pred( control, ewk_tf ) :
    return control * ewk_tf

def qcd_only( had_data, ewk_pred ) :
    return had_data - ewk_pred

def qcd_ratio( qcd_only, method="" ) :
    ratio = None
    if method == "abcd" :
        q = qcd_only
        r = utils.safe_div(q[:-1,:],q[1:,:])
        r = np.delete(r,r.shape[1]-1,1) # delete last row 
        r = np.insert(r,0,ufloat(np.nan,np.nan),1) # insert first row
        r = np.insert(r,0,ufloat(np.nan,np.nan),0) # insert first column
        ratio = r
    elif method == "double" :
        q = qcd_only
        r = utils.safe_div(q[:-1,:],q[1:,:])
        rr = utils.safe_div(r[:,1:-1]*r[:,1:-1],r[:,:-2])
        rr = np.insert(rr,(0,0),ufloat(np.nan,np.nan),1) # insert first two rows
        rr = np.insert(rr,0,ufloat(np.nan,np.nan),0) # insert first column
        ratio = rr
    else : # "raw"
        q = qcd_only
        r = utils.safe_div(q[:-1,:],q[1:,:])
        r = np.insert(r,0,ufloat(np.nan,np.nan),0) # insert first column
        ratio = r
    return ratio

def qcd_abcd_ratio( qcd_only ) :
    return qcd_ratio( qcd_only, "abcd" ) 

def qcd_double_ratio( qcd_only ) :
    return qcd_ratio( qcd_only, "double" ) 

def qcd_pred( qcd_only, qcd_ratio ) :
    r = qcd_ratio
    q = qcd_only
    pr = r*q
    pr = np.delete(pr,0,0) # delete first column 
    pr = np.insert(pr,pr.shape[0],ufloat(np.nan,np.nan),0) # insert last column
    return pr

################################################################################
# EXECUTE ######################################################################
################################################################################

if __name__=="__main__":

    data = load_pickle(True)

    had_ewk = data[("signal","EWK",ht_bin)]*regions["signal"][2]
    had_qcd = data[("signal","QCD",ht_bin)]*regions["signal"][2]

    mu_ewk = data[("muon","EWK",ht_bin)]*regions["muon"][2]
    mu_qcd = data[("muon","QCD",ht_bin)]*regions["muon"][2]

    if qcd_debug : 
        mu_qcd.fill( ufloat(0.,1.15) ) # Set QCD to zero for muon sample
        for i in range(had_qcd.shape[0]) :
            for j in range(had_qcd.shape[1]) :
                ii = had_qcd.shape[0]-i
                jj = had_qcd.shape[1]-j
                #exp = 8. # MHT/MET and AlphaT are independent
                #exp = jj/2. # MHT/MET and AlphaT are independent
                #exp = ii # MHT/MET and AlphaT are independent
                #exp = jj/2. + i/2. # MHT/MET and AlphaT are independent
                #exp = jj /  2. - abs(i-1.) # MHT/MET and AlphaT are independent
                exp = jj / ( 2. + abs(i-1.) / 4. ) - 1. # MHT/MET and AlphaT are dependent
                val = pow(10.,exp) 
                #val = 0.75*math.exp(-0.5*pow((mhtmet_bins[i]-2.4)/1.06,2.))*math.exp(70.0-120.0*alphat_bins[j])*1850.
                #val = 0.75*math.exp(-0.5*pow((mhtmet_bins[i]-(1.5+(alphat_bins[j]-0.5)*10.))/1.06,2.))*math.exp(70.0-120.0*alphat_bins[j])*580.
                had_qcd[i,j] = ufloat(val,math.sqrt(val)) 
                #had_qcd[i,j] = ufloat(exp,0.) 
                if had_qcd[i,j] < 1. : had_qcd[i,j] = ufloat(0.,1.15)
                #print j,i,had_qcd.shape[1]-j,(had_qcd.shape[1]-j)/2.,val

    had_sm = had_ewk + had_qcd
    mu_sm = mu_ewk + mu_qcd

    had_sig = np.empty_like(had_ewk)
    had_sig.fill( ufloat(0.,0.) )
    mu_sig = np.empty_like(had_ewk)
    mu_sig.fill( ufloat(0.,0.) )
    if add_signal is True :
        had_sig = data[("signal","Signal",ht_bin)]*regions["signal"][2]
        #mu_sig = data[("muon","Signal",ht_bin)]*regions["muon"][2]
    had_sm = had_sm + had_sig
    mu_sm = mu_sm + mu_sig

    mu_ewk_low_ht = data[("muon","EWK",150)]*regions["muon"][2]
    mu_qcd_low_ht = data[("muon","QCD",150)]*regions["muon"][2]
    mu_sm_low_ht = mu_ewk_low_ht + mu_qcd_low_ht

    had_data = data[("signal","Data",ht_bin)]
    mu_data = data[("muon","Data",ht_bin)]
    if use_data is False :
        had_data = had_sm
        mu_data = mu_sm

    from table import *
    tab = Table()

    effs = trigger.parse("../trigger_files/")
    diff,cumu = trigger.efficiencies(effs,ht_bin,alphat_bins,mhtmet_bins)
    diff[:,-1] = cumu[:,-1] # use cumu effs for final alphat bin
    #print diff[:,-1]
    #print cumu[:,-1]

    mu_eff_val = np.zeros_like(mu_data) 
    mu_eff_val.fill(muon_eff[0])
    mu_eff_err = np.zeros_like(mu_eff_val) 
    mu_eff_err.fill(muon_eff[1])
    mu_sys_val = np.ones_like(mu_eff_val)
    mu_sys_err = np.zeros_like(mu_eff_val) 
    mu_sys_err.fill(muon_eff[2])
    mu_effs = unumpy.uarray(mu_eff_val,mu_eff_err)
    mu_syst = unumpy.uarray(mu_sys_val,mu_sys_err)

    tab.newpage("Trigger efficiencies")
    tab.alphat_bins(effs[index][1])
    tab.mhtmet_bins(effs[index][2])
    tab.add_table(effs[index][3],"Had trigger effs")
    tab.alphat_bins(alphat_bins)
    tab.mhtmet_bins(mhtmet_bins)
    tab.add_table(diff,"Rebinned had trigger effs")
    #tab.add_table(effs[index][5],"Had trigger syst. uncert.")
    #tab.add_table(mu_effs*mu_syst,"Muon trigger effs (total uncert.)")
    tab.alphat_bins(effs[index][1])
    tab.mhtmet_bins(effs[index][2])

    tab.alphat_bins(alphat_bins)
    tab.mhtmet_bins(mhtmet_bins)

    had_data_title = "Had data (taken from SM MC)"
    had_data_raw = copy(had_data)
    if use_data is True :
        had_data = had_data/diff
        had_data_title = "Had data (trigger eff corrected)"

    tab.newpage("Hadronic trigger and data")
    tab.add_table(diff,"Rebinned had trigger effs")
    tab.add_table(had_data_raw,"Had data raw")
    tab.add_table(had_data,had_data_title)

    tab.newpage("Hadronic data and MC yields")
    tab.add_table(had_data,had_data_title)
    tab.add_table(had_sm,"Had SM MC")
    tab.add_table(utils.safe_div(had_data,had_sm),"Had data / SM MC")
    tab.add_table(utils.safe_div(had_data,had_ewk),"Had data / EWK MC")

    tab.newpage("Hadronic data and MC yields")
    had_data_raw = copy(had_data)
    if use_data is True :
        had_data = had_data/diff
        tab.add_table(had_data,"Had data (trigger eff corrected)")
    else :
        tab.add_table(had_data,"Had data (taken from SM MC)")
    tab.add_table(had_sm,"Had SM MC")
    tab.add_table(utils.safe_div(had_data,had_sm),"Had data / SM MC")

    tab.newpage("Hadronic EWK and QCD yields from MC")
    tab.add_table(had_ewk,"Had EWK MC")
    tab.add_table(had_qcd,"Had QCD MC")
    tab.add_table(utils.safe_div(had_qcd,had_sm),"Fraction QCD/SM MC")

    tab.newpage("Muon data and MC yields")
    mu_data_raw = copy(mu_data)
    if use_data is True :
        mu_data = mu_data/(mu_effs*mu_syst)
        tab.add_table(mu_data,"Mu data (trigger eff corrected)")
    else :
        tab.add_table(mu_data,"Mu data (taken from SM MC)")
    tab.add_table(mu_sm,"Mu SM MC")
    tab.add_table(utils.safe_div(mu_data,mu_sm),"Mu data / SM MC")

    tab.newpage("Muon EWK and QCD yields from MC")
    tab.add_table(mu_ewk,"Mu EWK MC")
    tab.add_table(mu_qcd,"Mu QCD MC")
    tab.add_table(utils.safe_div(mu_qcd,mu_sm),"Fraction QCD/SM MC")

    tab.newpage("Muon EWK and QCD yields from MC at low HT")
    tab.add_table(mu_ewk_low_ht,"Mu EWK MC")
    tab.add_table(mu_qcd_low_ht,"Mu QCD MC")
    tab.add_table(utils.safe_div(mu_qcd_low_ht,mu_sm_low_ht),"Fraction QCD/SM MC")

    tab.newpage("Had and muon SIG yields from MC")
    tab.add_table(had_sig,"Had SIG MC")
    tab.add_table(mu_sig,"Mu SIG MC")
    tab.add_table(utils.safe_div(mu_sig,had_sig),"Fraction Muon/Had SIG MC")
    
    tf = ewk_tf(had_ewk,mu_ewk)
    tab.newpage("Transfer factors")
    tab.add_table(had_ewk,"Had EWK MC")
    tab.add_table(mu_ewk,"Mu EWK MC")
    tab.add_table(tf,"Transfer factors")

    ewk = ewk_pred(mu_data,tf)
    tab.newpage("EWK predictions")
    tab.add_table(tf,"Transfer factors")
    tab.add_table(mu_data,"Mu data (trigger eff corrected)")
    tab.add_table(ewk,"EWK prediction")
    
    qcd = qcd_only(had_data,ewk)
    tab.newpage("'QCD only' yields (Data - EWK pred)")
    tab.add_table(had_data,"Had data (trigger eff corrected)")
    tab.add_table(ewk,"EWK prediction")
    tab.add_table(qcd,"QCD only")

    ratio_raw = qcd_ratio(qcd)
    ratio_abcd = qcd_abcd_ratio(qcd)
    ratio_double = qcd_double_ratio(qcd)

    tab.newpage("QCD ratios (ABCD and 'double ratio' methods)")
    #tab.add_table(qcd,"QCD only")
    tab.add_table(ratio_raw,"QCD ratios (raw)")
    tab.add_table(ratio_abcd,"QCD ratios ABCD")
#    tab.add_table(ratio_double,"QCD ratios 'double ratio'")
    #tab.add_table(ratio_double/ratio_abcd,"QCD ratios 'double ratio' / ABCD")

    pred_abcd = qcd_pred(qcd,ratio_abcd)
    #if use_data is True : pred_abcd *= diff 

#    pred_double = qcd_pred(qcd,ratio_double)

    tab.newpage("QCD predictions (ABCD method)")
    tab.add_table(qcd,"QCD only")
    tab.add_table(ratio_abcd,"QCD ratios")
    tab.add_table(pred_abcd,"QCD pred ABCD (trigger corr)")

#    tab.newpage("QCD predictions ('double ratio' method)")
#    tab.add_table(qcd,"QCD only")
#    tab.add_table(ratio_double,"QCD ratios")
#    tab.add_table(pred_double,"QCD pred 'double ratio' (trigger corr, stat. uncert.)")
    
#    tab.newpage("Comparison of QCD predictions (ABCD vs 'double ratio')")
#    tab.add_table(pred_abcd,"QCD pred ABCD (trigger corr, stat. uncert.)")
#    tab.add_table(pred_double,"QCD pred 'double ratio' (trigger corr, stat. uncert.)")
#    tab.add_table(utils.safe_div(pred_double,pred_abcd),"'double ratio' / ABCD (trigger corr, stat. uncert.)")

    tab.newpage("QCD prediction (ABCD) vs hadronic yields")
    tab.add_table(pred_abcd,"QCD pred ABCD (trigger corr)")
    tab.add_table(had_data_raw,"Had data (no trigger corr)")
    tab.add_table(utils.safe_div(pred_abcd,had_data_raw),"QCD pred ABCD / had data (no trigger corr)")

    tab.newpage("QCD prediction (ABCD) vs QCD only")
    tab.add_table(pred_abcd,"QCD pred ABCD (trigger corr)")
    tab.add_table(qcd,"QCD only (no trigger corr)")
    tab.add_table(utils.safe_div(qcd,pred_abcd),"QCD only / QCD pred")

    tab.newpage("QCD prediction (ABCD) vs EWK pred")
    tab.add_table(pred_abcd,"QCD pred ABCD (trigger corr)")
    tab.add_table(ewk,"EWK pred (no trigger corr)")
    tab.add_table(utils.safe_div(pred_abcd,ewk),"QCD pred ABCD / EWK pred")

    tab.newpage("'Signal': Data - ( EWK pred + QCD pred)")
    tab.add_table(had_data,"Had data (trigger eff corrected)")
    tab.add_table(pred_abcd,"QCD pred ABCD (trigger corr)")
    tab.add_table(had_data-ewk-pred_abcd,"Data - ( EWK pred + QCD pred)")

    if ( add_signal ) :
        tab.newpage("Fraction of signal recovered")
        excess = had_data - ewk - pred_abcd
        tab.add_table(excess,"Excess: Data - (EWK pred + QCD pred)")
        tab.add_table(had_sig,"Signal")
        tab.add_table(utils.safe_div(had_sig-excess,had_sig),"(Signal-Excess)/Signal")
        #tab.add_table((had_sig-excess)/had_sig,"(Signal-Excess)/Signal")
        #tab.add_table(excess-excess,"(Signal-Excess)/Signal")

    tab.newpage("Closure")
    tab.add_table(had_qcd,"Had QCD MC")
    tab.add_table(ratio_abcd,"QCD ratios")
    tab.add_table(utils.safe_div(qcd,pred_abcd),"QCD only / QCD pred")
    #tab.add_table(utils.safe_div(had_qcd,pred_abcd),"QCD only / QCD pred")

#    tab.newpage("QCD prediction ('double ratio') vs hadronic yields")
#    tab.add_table(pred_double,"QCD pred double (trigger corr)")
#    tab.add_table(had_data_raw,"Had data (no trigger corr)")
#    tab.add_table(utils.safe_div(pred_double,had_data_raw),"QCD pred double / had data (no trigger corr)")
    
    from array import *
    import ROOT as r

    if False :
        his_title = "(obs-pred)/pred for ABCD method"
        his_contents = utils.safe_div(had_data_raw-pred_abcd,pred_abcd).flatten().tolist()
        if use_double_ratio == True :
            his_title = "(obs-pred)/pred for 'double ratio' method" 
            his_contents = utils.safe_div(had_data_raw-pred_double,pred_double).flatten().tolist()
            
        his = r.TH1D(his_title,his_title,100,-2.,2.)
        # his = r.TH1D(his_title,his_title,100,-5.,5.)
        for i in his_contents : 
            if i.s != 0. : 
                # his.Fill(i.n,abs((i.n-1)/i.s))
                his.Fill( i.n )
                his.Fill( (i.n-1.)/i.s )
        his.Draw()
        input("Please any key to continue...")

    if False :
        his_title = "Pred (ABCD method) "
        his_contents = utils.safe_div(had_data_raw-pred_abcd,pred_abcd).flatten().tolist()
        if use_double_ratio == True :
            his_title = "(obs-pred)/pred for 'double ratio' method" 
            his_contents = utils.safe_div(had_data_raw-pred_double,pred_double).flatten().tolist()
            
        his = r.TH1D(his_title,his_title,100,-2.,2.)
        # his = r.TH1D(his_title,his_title,100,-5.,5.)
        for i in his_contents : 
            if i.s != 0. : 
                # his.Fill(i.n,abs((i.n-1)/i.s))
                his.Fill( i.n )
                his.Fill( (i.n-1.)/i.s )
        his.Draw()
        input("Please any key to continue...")



