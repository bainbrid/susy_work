import numpy as np
import math
import qcd
import uncertainties

data = qcd.load_pickle()

had_data = data[("signal","Data",0)]
had_ewk = data[("signal","EWK",0)]
had_qcd = data[("signal","QCD",0)]
had_sig = data[("signal","Signal",0)]

mu_data = data[("muon","Data",0)]
mu_ewk = data[("muon","EWK",0)]
mu_qcd = data[("muon","QCD",0)]

