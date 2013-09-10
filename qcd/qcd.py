import sys
import numpy as np
import pickle
import math
from copy import copy
from table import *
tab = Table()
np.set_printoptions(precision=3)

################################################################################
# OPTIONS ######################################################################
################################################################################

pois_h = [1.15,1.36,1.53,1.73,1.98,2.21,2.42,2.61,2.80,3.00]
pois_l = [0.00,1.00,2.00,2.14,2.30,2.49,2.68,2.86,3.03,3.19]

# (filename prefix, dir prefix, lumi, % uncert)
regions = {"signal":("Had","QCD",18.583*10.,0.04),
           "muon":("Muon","QCD_OneMuon",19.255*10.,0.04),}

processes = ["Data","EWK","QCD","SM","Signal","SM_Signal"]

pickle_file = "data.pkl"

first_bin = 6
last_bin = None

use_data = False

double_ratio = False

# getter
def entry( list, str ) :
    if len(list) == 0 :
        return None
    strs = [ x[-1] for x in list ]
    if str in strs :
        return strs.index(str),alphat_bins[strs.index(str)]
    else :
        return None,(None,)*len(list[0])

# (bin lower bound, alphat bin, filename)
ht_bins = [(200,6,"200_275_73_73_36"),
           (275,4,"275_325_73_73_36"),
#           (325,2,"325_375_86_86_43"),
#           (375,1,"375_475_100_100_50"),
#           (475,0,"475_575_100_100_50"),
#           (575,0,"575_675_100_100_50"),
#           (675,0,"675_775_100_100_50"),
#           (775,0,"775_875_100_100_50"),
#           (875,0,"875_100_100_50"),
           ]

# (bin lower bound, bin centre)
alphat_bins = [(0.51,"0.515000"), 
               (0.52,"0.525000"), 
               (0.53,"0.535000"), 
               (0.54,"0.545000"), 
               (0.55,"0.555000"), 
               (0.56,"0.565000"), 
               (0.57,"0.575000"), 
               (0.58,"0.585000"), 
               (0.59,"0.595000"), 
#               (0.60,"0.645000"), # cumu
               (0.60,"0.650000"), # diff
               (0.65,"0.750000"), # diff
               (0.70,"0.745000"), # cumu
               ]

# (bin lower bound, filename)
mhtmet_bins = [(0.00,"_mht"),  
               (1.25,""),
               (2.50,""),
               (3.75,""),
               (5.00,""),
               ]

# (eff, stat uncert, syst uncert)
muon_eff = (0.88,0.00,0.03)

# (syst uncert, debug mode, debug eff, debug syst uncert)
trigger_syst = (0.03,False,1.,0.)

################################################################################
# METHODS ######################################################################
################################################################################

def view(arr) :
    return np.flipud(copy(arr).transpose())

def summary(x) :
    print "X axis binning:",x[0]
    print "Y axis binning:",x[1]
    print "Values:"
    print view(x[2])
    print "Stath:"
    print view(x[3])
    print "Statl:"
    print view(x[4])
    print "Systh:"
    print view(x[5])
    print "Systl:"
    print view(x[6])
    print "Toth:"
    print view(x[7])
    print "Totl:"
    print view(x[8])

def trigger_effs( path, verbose=False, final_bin=4 ) :
    data = {}
    for iht in range(len(ht_bins)) :
        effs = np.zeros((len(mhtmet_bins),len(alphat_bins))) 
        errh = np.zeros_like(effs)
        errl = np.zeros_like(effs)
        sysh = np.zeros_like(effs)
        sysl = np.zeros_like(effs)
        toth = np.zeros_like(effs)
        totl = np.zeros_like(effs)
        for imht in range(len(mhtmet_bins)) :
            name = "text_HT"+ht_bins[iht][2]+mhtmet_bins[imht][1]+"_AlphaT_ge2j.txt"
            file = open(path+name)
            for line in file.readlines() :
                entries = line.split()
                tuple = entry(alphat_bins,entries[2])
                if tuple[0] is not None :
                    if trigger_syst[1] == True :
                        effs[imht][tuple[0]] = trigger_syst[2]
                        errh[imht][tuple[0]] = trigger_syst[3]
                        errl[imht][tuple[0]] = trigger_syst[3]
                    else :
                        effs[imht][tuple[0]] = entries[3]
                        errh[imht][tuple[0]] = entries[5]
                        errl[imht][tuple[0]] = entries[7]
            file.close()
        sysh.fill(trigger_syst[0])
        sysl.fill(trigger_syst[0])
        toth = np.sqrt(errh*errh+sysh*sysh)
        totl = np.sqrt(errl*errl+sysl*sysl)
        data[iht] = (mhtmet_bins,alphat_bins,effs,errh,errl,sysh,sysl,toth,totl)
        if final_bin > -1 and iht > final_bin : 
            data[iht] = data[final_bin]
        if verbose == True :
            print "TRIGGER EFFS, HT:",ht_bins[iht][0]
            summary(data[iht])
    return data

def load_pickle( verbose=False ) :
    data = pickle.load(open(pickle_file,"r"))
    for i in regions.keys() :
        for j in processes :
            for k in range(0,len(ht_bins)-1) :
                if (i,j,k) in data.keys() : 
                    x = data[(i,j,k)]
                    if verbose == True :
                        print "YIELDS, region:",i," sample:",j," HT:",ht_bins[k][0]
                        summary(x)
                else :
                    print "Problem:",i,j,k
    return data

def ewk_tf( data, verbose=False ) :
    tf = {}
    for k in range(0,len(ht_bins)-1) :
        signal = ("signal","EWK",k)
        muon = ("muon","EWK",k)
        if signal in data.keys() and muon in data.keys() :
            xbins = data[signal][0]
            ybins = data[signal][1]
            sig = data[signal][2]
            mu = data[muon][2]
            val = sig/mu
            sig_stath = data[signal][3]
            sig_statl = data[signal][4]
            mu_stath = data[muon][3]
            mu_statl = data[muon][4]
            stath = np.sqrt(sig_stath*sig_stath + mu_stath*mu_stath)
            statl = np.sqrt(sig_statl*sig_statl + mu_statl*mu_statl)
            sig_systh = data[signal][5]
            sig_systl = data[signal][6]
            mu_systh = data[muon][5]
            mu_systl = data[muon][6]
            systh = np.sqrt(sig_systh*sig_systh + mu_systh*mu_systh)
            systl = np.sqrt(sig_systl*sig_systl + mu_systl*mu_systl)
            toth = np.sqrt( stath*stath + systh*systh )
            totl = np.sqrt( statl*statl + systl*systl )
            tf[k] = (xbins,ybins,val,stath,statl,systh,systl,toth,totl)
            if verbose == True :
                print "TRANSFER FACTORS, HT:",ht_bins[k][0]
                summary(tf[k])
        else :
            print "Problem:","EWK",k
    return tf

def calc( val1, val2, multiply=True ) :
    val = val1[0]*val2[0]
    if multiply == False : 
        if val2[0] > 0. : val = val1[0]/val2[0]
        else : return None
    val_stath = pois_h[0]
    val_statl = pois_l[0]
    if val1 > 0. and val2 > 0. :
        val1_stath = val1[1] / val1[0]
        val1_statl = val1[2] / val1[0]
        val2_stath = val2[1] / val2[0]
        val2_statl = val2[2] / val2[0]
        stath = math.sqrt( val1_stath*val1_stath + val2_stath*val2_stath )
        statl = math.sqrt( val1_statl*val1_statl + val2_statl*val2_statl )
    return None

def ewk_pred( data, ewk_tf, verbose=False ) :
    pred = {}
    for k in range(0,len(ht_bins)-1) :
        obs = None
        if use_data is True : 
            obs = ("muon","Data",k)
        else :
            obs = ("muon","SM",k)
        if obs in data.keys() and k in ewk_tf.keys() :
            xbins = data[obs][0]
            ybins = data[obs][1]
            val = data[obs][2] * ewk_tf[k][2]
            obs_stath = data[obs][3]/data[obs][2]
            obs_statl = data[obs][4]/data[obs][2]
            obs_systh = 0.
            obs_systl = 0.
            tf_stath = ewk_tf[k][3]/ewk_tf[k][2]
            tf_statl = ewk_tf[k][4]/ewk_tf[k][2]
            tf_systh = ewk_tf[k][5]/ewk_tf[k][2]
            tf_systl = ewk_tf[k][6]/ewk_tf[k][2]
            stath = np.sqrt(obs_stath*obs_stath+tf_stath*tf_stath)
            statl = np.sqrt(obs_statl*obs_statl+tf_statl*tf_statl)
            systh = np.zeros_like(stath) #np.sqrt(obs_systh*obs_systh+tf_systh*tf_systh)
            systl = np.zeros_like(stath) #np.sqrt(obs_systl*obs_systl+tf_systl*tf_systl)
            toth = np.sqrt( stath*stath + systh*systh )
            totl = np.sqrt( statl*statl + systl*systl )
            pred[k] = (xbins,ybins,val,stath,statl,systh,systl,toth,totl)
            if verbose == True :
                print "EWK PRED, HT:",ht_bins[k][0]
                summary(pred[k])
        else :
            print "Problem:","EWK",k
    return pred

def qcd_only( data, ewk_pred, verbose=False ) :
    qcd = {}
    for k in range(0,len(ht_bins)-1) :
        obs = None
        if use_data is True : 
            obs = ("signal","Data",k)
        else :
            obs = ("signal","SM",k)
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
            qcd[k] = (xbins,ybins,val,errh,errl,errh,errl,errh,errl)
            if verbose == True :
                print "QCD ONLY, HT:",ht_bins[k][0]
                summary(qcd[k])
        else :
            print "Problem:","EWK",k
    return qcd

def qcd_ratio( qcd_only, verbose=False ) :
    ratio = {}
    for k in range(0,len(ht_bins)-1) :
        ratio[k] = (qcd_only[k][0],
                    qcd_only[k][1],
                    np.insert(qcd[0][2][:-1,:]/qcd[0][2][1:,:],0,0.,0),
                    qcd_only[k][3][0,:],
                    qcd_only[k][4][0,:],
                    qcd_only[k][5][0,:],
                    qcd_only[k][6][0,:],
                    qcd_only[k][7][0,:],
                    qcd_only[k][8][0,:])
        if verbose == True :
            print "QCD RATIO, HT:",ht_bins[k][0]
            summary(ratio[k])
    return ratio

def qcd_pred( qcd_only, qcd_ratio, verbose=False ) :
    pred = {}
    for k in range(0,len(ht_bins)-1) :
        pr = None
        r = qcd_ratio[k][2]
        q = qcd_only[k][2]
        z = np.zeros_like(r)
        if double_ratio is True :
            rr = r[:,1:-1]*r[:,1:-1]/r[:,:-2]
            rr = np.hstack((z[:,0:2],rr))
            pr = np.delete(rr*q,0,0)
            pr = np.vstack((pr,z[0:1,:]))
        else :
            pr = r[:,:-1]*q[:,1:]
            pr = np.hstack((z[:,0:1],pr))
            pr = np.delete(pr,0,0)
            pr = np.vstack((pr,z[0:1,:]))
        pred[k] = (qcd_only[k][0],
                   qcd_only[k][1],
                   pr,
                   np.insert(qcd_only[k][3][1,1:],0,0.),
                   np.insert(qcd_only[k][4][1,1:],0,0.),
                   np.insert(qcd_only[k][5][1,1:],0,0.),
                   np.insert(qcd_only[k][6][1,1:],0,0.),
                   np.insert(qcd_only[k][7][1,1:],0,0.),
                   np.insert(qcd_only[k][8][1,1:],0,0.))
        if use_data is True :
            pred[k] = broadcast( pred[k], effs[k] )
        if verbose == True :
            print "QCD PRED, HT:",ht_bins[k][0]
            summary(pred[k])
    return pred

################################################################################
# EXECUTE ######################################################################
################################################################################

if __name__=="__main__":

    had_ewk = ("signal","EWK",0)
    had_qcd = ("signal","QCD",0)
    had_sm = ("signal","SM",0)
    had_sig = ("signal","Signal",0)
    had_sm_sig = ("signal","SM_Signal",0)

    mu_ewk = ("muon","EWK",0)
    mu_qcd = ("muon","QCD",0)
    mu_sm = ("muon","SM",0)
    
    had_data = ("signal","Data",0)
    mu_data = ("muon","Data",0)
    if use_data is False :
        had_data = ("signal","SM",0)
        mu_data = ("muon","SM",0)

    effs = trigger_effs("./trigger_files/")
    data = load_pickle()

    tab.add_table(data[had_data],"Had data (raw)")
    tab.add_table(data[had_sm],"Had SM MC (raw)")
    # tab.add_table(data[had_qcd]/data[had_sm],"Fraction QCD/SM MC (raw)")

    tab.newpage()
    tab.add_table(data[had_sm],"Had SM (raw)")
    tab.add_table(data[had_ewk],"Had EWK (raw)")
    tab.add_table(data[had_qcd],"Had QCD (raw)")
    
    tab.newpage()
    tab.add_table(data[mu_data],"Muon data (raw)")
    tab.add_table(data[mu_sm],"Muon SM MC (raw)")
    # tab.add_table(data[mu_qcd]/data[mu_sm],"Fraction QCD/SM MC (raw)")
    
    tab.newpage()
    tab.add_table(data[mu_sm],"Muon SM (raw)")
    tab.add_table(data[mu_ewk],"Muon EWK (raw)")
    tab.add_table(data[mu_qcd],"Muon QCD (raw)")
    
    tab.newpage()
    tab.add_table(data[mu_sm],"Muon SM (raw)")
    
    print "Trigger effs"
    print view(effs[0][2][:,first_bin:last_bin])
    
    print "Muon data (raw)"
    print view(data[mu_data][2][:,first_bin:last_bin])
    print "Muon EWK (raw)"
    print view(data[mu_ewk][2][:,first_bin:last_bin])
    print "Muon EWK ERR (raw)"
    print view(data[mu_ewk][3][:,first_bin:last_bin])
    print "Muon QCD (raw)"
    print view(data[mu_qcd][2][:,first_bin:last_bin])
    print "Muon QCD ERR (raw)"
    print view(data[mu_qcd][3][:,first_bin:last_bin])
    
    print "Had data (raw)"
    print view(data[had_data][2][:,first_bin:last_bin])
    print "Had SM (raw)"
    print view(data[had_sm][2][:,first_bin:last_bin])
    print "Had EWK (raw)"
    print view(data[had_ewk][2][:,first_bin:last_bin])
    print "Had EWK ERR (raw)"
    print view(data[had_ewk][3][:,first_bin:last_bin])
    print "Had QCD (raw)"
    print view(data[had_qcd][2][:,first_bin:last_bin])
    print "Had QCD ERR (raw)"
    print view(data[had_qcd][3][:,first_bin:last_bin])
    print "Had Signal (raw)"
    print view(data[had_sig][2][:,first_bin:last_bin])
    print "Had Signal ERR (raw)"
    print view(data[had_sig][3][:,first_bin:last_bin])
    
    # Hack to deal with lumis and trigger effs
    def scale( tmp, factor ) : return ( tmp[0], tmp[1], factor*tmp[2], factor*tmp[3], factor*tmp[4], tmp[5], tmp[6], factor*tmp[7], factor*tmp[8] )
    def broadcast( tmp, arr ) : return ( tmp[0], tmp[1], tmp[2]*arr[2], tmp[3]*arr[2], tmp[4]*arr[2], tmp[5], tmp[6], tmp[7]*arr[2], tmp[8]*arr[2] )
    def broadcast_divide( tmp, arr ) : return ( tmp[0], tmp[1], tmp[2]/arr[2], tmp[3]/arr[2], tmp[4]/arr[2], tmp[5], tmp[6], tmp[7]/arr[2], tmp[8]/arr[2] )
    
    data[had_ewk] = scale( data[had_ewk], 18.583*10. )
    data[had_qcd] = scale( data[had_qcd], 18.583*10. )
    # data[had_sm] = scale( data[had_sm], 18.583*10. )
    data[had_sig] = scale( data[had_sig], 18.583*10. )
    data[had_sm_sig] = scale( data[had_sm_sig], 18.583*10. )
    
    data[mu_ewk] = scale( data[mu_ewk], 19.255*10. )
    data[mu_qcd] = scale( data[mu_qcd], 19.255*10. )
    data[mu_sm] = scale( data[mu_sm], 19.255*10. )
    
    if use_data is True :
        data[had_data] = broadcast_divide( data[had_data], effs[0] )
        data[mu_data] = scale( data[mu_data], 1./muon_eff[0] )
    
    tf = ewk_tf(data)
    ewk = ewk_pred(data,tf)
    qcd = qcd_only(data,ewk)
    ratio = qcd_ratio(qcd)
    pred = qcd_pred(qcd,ratio)
    
    print "Muon data"
    print view(data[mu_data][2][:,first_bin:last_bin])
    print "Muon EWK"
    print view(data[mu_ewk][2][:,first_bin:last_bin])
    print "Muon QCD"
    print view(data[mu_qcd][2][:,first_bin:last_bin])
    
    print "Muon SM/Data"
    print view(data[mu_sm][2][:,first_bin:last_bin]/data[mu_data][2][:,first_bin:last_bin])
    
    print "Had data"
    print view(data[had_data][2][:,first_bin:last_bin])
    print "Had SM"
    print view(data[had_sm][2][:,first_bin:last_bin])
    print "Had EWK"
    print view(data[had_ewk][2][:,first_bin:last_bin])
    print "Had QCD"
    print view(data[had_qcd][2][:,first_bin:last_bin])
    print "Had Signal"
    print view(data[had_sig][2][:,first_bin:last_bin])
    
    print "Transfer factors"
    print view(tf[0][2][:,first_bin:last_bin])
    
    print "EWK pred"
    print view(ewk[0][2][:,first_bin:last_bin])
    
    tmp = qcd[0][2][:,first_bin:last_bin]
    if use_data is True :
        tmp = broadcast(qcd[0],effs[0])[2][:,first_bin:last_bin]
    print "QCD only"
    print view(tmp)
    
    print "QCD ratio"
    print view(ratio[0][2][:,first_bin:last_bin])
    
    print "QCD pred"
    print view(pred[0][2][:,first_bin:last_bin])
    
    print "QCD pred/only"
    pred_over_qcd = pred[0][2][:,first_bin:last_bin]/tmp
    print view(pred_over_qcd)

    # tab.add_table(data[had_sm])
    # tab.add_table(data[had_sm])
    # tab.add_table(data[had_sm])
    del tab

    # from array import *
    # import ROOT as r
    # his = r.TH1D("","",20,0.,2.)
    # str = "Double Ratio" if double_ratio is True else "ABCD"
    # his = r.TH1D(str,str,20,0.,2.)
    # for i in pred_over_qcd.flatten().tolist() : his.Fill(i)
    # his.Draw()
    # input("Please any key to continue...")

