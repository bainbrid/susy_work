import sys
import numpy as np
import pickle
import math
import utils
from copy import copy
from uncertainties import ufloat
from uncertainties import UFloat
from uncertainties import unumpy

################################################################################
# CONFIG #######################################################################
################################################################################

# bin lower bound, alphat bin, filename string
ht_bins = [
    (200,6,"200_275_73_73_36"),
    (275,4,"275_325_73_73_36"),
    (325,2,"325_375_86_86_43"),
    (375,1,"375_475_100_100_50"),
    (475,0,"475_575_100_100_50"),
    (575,0,"575_675_100_100_50"),
    (675,0,"675_775_100_100_50"),
    (775,0,"775_875_100_100_50"),
    (875,0,"875_100_100_50"),
#    (875,0,"875_975_100_100_50"),
#    (975,0,"975_1075_100_100_50"),
#    (1075,0,"1075_100_100_50"),
    ]

# bin lower bound and "bin centre" string for diff/cumu effs
alphat_bins = [
    (0.50,"0.505000","0.500000"), 
    (0.51,"0.515000","0.510000"), 
    (0.52,"0.525000","0.520000"), 
    (0.53,"0.535000","0.530000"), 
    (0.54,"0.545000","0.540000"), 
    (0.55,"0.555000","0.550000"), 
    (0.56,"0.565000","0.560000"), 
    (0.57,"0.575000","0.570000"), 
    (0.58,"0.585000","0.580000"), 
    (0.59,"0.595000","0.590000"), 
    (0.60,"0.650000","0.645000"),
    (0.70,"0.750000","0.745000"),
    (0.80,"0.850000","0.845000"),
    (0.90,"0.950000","0.945000"),
    (1.00,"1.050000","1.045000"),
    ]

# bin lower bound, filename string
mhtmet_bins = [
    (0.00,"_mht"),  
    (1.25,""),
    (2.50,""),
    (3.75,""),
    (5.00,""),
    ]

# alphat systematic uncertainty
alphat_syst = 0.03

# debug mode on/off, debug eff, debug syst uncert
debug_mode = (False,1.,0.)

# muon trigger eff, stat uncert, syst uncert
muon_eff = (0.88,0.00,0.03)

# utility method to pull correct the tuple from the alphat_bins list, keyed by strings in the text file
def alphat_bin( list, str, diff=True ) :
    if len(list) == 0 :
        return None
    strs = [ x[1] for x in list ] if diff is True else [ x[2] for x in list ]
    check = False 
    if str[-1] == "Differential" and diff is True : check = True
    elif str[-1] == "Cumulative" and diff is False : check = True
    if check is True and str[2] in strs :
        return strs.index(str),alphat_bins[strs.index(str)]
    else :
        return None,(None,)*len(list[0])

def parse( path, final_bin=-1 ) :
    data = {}
    for iht in range(len(ht_bins)) :
        diff_vals = np.zeros((len(mhtmet_bins),len(alphat_bins))) 
        diff_errs = np.zeros_like(diff_vals) 
        cumu_vals = np.zeros_like(diff_vals) 
        cumu_errs = np.zeros_like(diff_vals) 
        for imht in range(len(mhtmet_bins)) :
            name = "text_HT"+ht_bins[iht][2]+mhtmet_bins[imht][1]+"_AlphaT_ge2j.txt"
            file = open(path+name)
            for line in file.readlines() :
                entries = line.split()
                strs = []
                if entries[-1] == "Differential" : strs = [ x[1] for x in alphat_bins ] 
                if entries[-1] == "Cumu" : strs = [ x[2] for x in alphat_bins ] 
                if len(strs) > 0 and entries[2] in strs :
                    alphat_bin = strs.index(entries[2])
                    err_max = entries[5] if entries[5] > entries[7] else entries[7]
                    if entries[-1] == "Differential" :
                        diff_vals[imht][alphat_bin] = entries[3] if debug_mode[0] is False else debug_mode[1]
                        diff_errs[imht][alphat_bin] = err_max if debug_mode[0] is False else debug_mode[2] 
                    elif entries[-1] == "Cumu" : 
                        cumu_vals[imht][alphat_bin] = entries[3] if debug_mode[0] is False else debug_mode[1]
                        cumu_errs[imht][alphat_bin] = err_max if debug_mode[0] is False else debug_mode[2]
#                        if alphat_bin == len(alphat_bins)-1 : # use cumu effs for last open bin
#                            diff_vals[imht][alphat_bin] = entries[3] if debug_mode[0] is False else debug_mode[1]
#                            diff_errs[imht][alphat_bin] = err_max if debug_mode[0] is False else debug_mode[2] 
            file.close()
        ht_binning = [ x[0] for x in ht_bins ]
        alphat_binning = [ x[0] for x in alphat_bins ]
        mhtmet_binning = [ x[0] for x in mhtmet_bins ]
        diff = unumpy.uarray(diff_vals,diff_errs)
        cumu = unumpy.uarray(cumu_vals,cumu_errs)
        errs = np.zeros_like(diff_vals) 
        errs.fill(alphat_syst)
        syst = unumpy.uarray( np.ones_like(diff_vals), errs )
        data[iht] = (ht_binning,alphat_binning,mhtmet_binning,diff,cumu,syst)
        if final_bin > -1 and iht > final_bin : 
            data[iht] = data[final_bin]
    return data

def efficiency( effs, ht_val, alphat_val, mhtmet_val, diff=True ) :

    ht_binning = effs[0][0]
    ht_bin = None
    for bin in range(len(ht_binning)) :
        if bin == len(ht_binning)-1 and ht_val >= ht_binning[bin] : ht_bin = bin
        elif ht_val >= ht_binning[bin] and ht_val < ht_binning[bin+1] : ht_bin = bin
    if ht_bin is None : return None

    alphat_binning = effs[ht_bin][1]
    alphat_bin = None
    for bin in range(len(alphat_binning)) :
        if bin == len(alphat_binning)-1 and alphat_val >= alphat_binning[bin] : alphat_bin = bin
        elif alphat_val >= alphat_binning[bin] and alphat_val < alphat_binning[bin+1] : alphat_bin = bin
    if alphat_bin is None : return None

    mhtmet_binning = effs[ht_bin][2]
    mhtmet_bin = None
    for bin in range(len(mhtmet_binning)) :
        if bin == len(mhtmet_binning)-1 and mhtmet_val >= mhtmet_binning[bin] : mhtmet_bin = bin
        elif mhtmet_val >= mhtmet_binning[bin] and mhtmet_val < mhtmet_binning[bin+1] : mhtmet_bin = bin
    if mhtmet_bin is None : return None

    index = 3 if diff is True else 4
    lower = effs[ht_bin][index][mhtmet_bin][alphat_bin]
    upper = effs[ht_bin][index][mhtmet_bin][alphat_bin+1]
    diff = upper - lower

    val = None
    if alphat_bin == len(alphat_binning)-1 or lower == upper : val = lower
    else : 
        bin_lower = alphat_binning[alphat_bin]
        bin_upper = alphat_binning[alphat_bin+1]
        val = lower + diff * ( (alphat_val-bin_lower) / (bin_upper-bin_lower) )
    return val

def efficiencies( input, ht_val, alphat_binning, mhtmet_binning ) :
    vals = np.zeros((len(mhtmet_binning),len(alphat_binning))) 
    errs = np.zeros_like(vals) 
    diff = unumpy.uarray(vals,errs)
    cumu = unumpy.uarray(vals,errs)
    for imht in range(len(mhtmet_binning)) :
        mhtmet_val = mhtmet_binning[imht]
        for ialphat in range(len(alphat_binning)) :
            alphat_val = alphat_binning[ialphat]
            diff[imht][ialphat] = efficiency( input, ht_val, alphat_val, mhtmet_val, True ) 
            cumu[imht][ialphat] = efficiency( input, ht_val, alphat_val, mhtmet_val, False ) 
#            if ialphat == len(alphat_binning)-1 :
#                diff[imht][ialphat] = cumu[imht][ialphat]
    return diff,cumu

################################################################################
# EXECUTE ######################################################################
################################################################################

if __name__=="__main__":

    ht = [200,275,325,375,475,575,675,775,875]
    alphat = [0.50,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,
              0.60,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,
              0.70]
    
    effs = parse("../trigger_files/")

    from table import *
    tab = Table()
    tab.mhtmet_bins(effs[0][2])
    
    for index in range(len(ht)) :

        diff,cumu = efficiencies(effs,ht[index],alphat,effs[index][2])

        str = "trigger efficiencies for " 
        if index == len(ht)-1 : 
            str += ("$H_{\\textrm{T}} > %i \, \\textrm{GeV}$"%ht[index])
        else :
            str += ("$%i < H_{\\textrm{T}} < %i \, \\textrm{GeV}$"%(ht[index],ht[index+1]))

        tab.newpage("Differential "+str)
        tab.alphat_bins(effs[index][1])
        tab.add_table(effs[index][3],"Hadronic trigger differential efficiencies (stat. uncert. only)")
        tab.alphat_bins(alphat)
        tab.add_table(diff,"Rebinned hadronic trigger differential efficiencies (stat. uncert. only)")

        tab.newpage("Cumulative "+str)
        tab.alphat_bins(effs[index][1])
        tab.add_table(effs[index][4],"Hadronic trigger cumulative efficiencies (stat. uncert. only)",cumu=True)
        tab.alphat_bins(alphat)
        tab.add_table(cumu,"Rebinned hadronic trigger cumulative efficiencies (stat. uncert. only)",cumu=True)
