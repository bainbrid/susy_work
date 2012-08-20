#!/usr/bin/env python

import setupSUSY
from libFrameworkSUSY import *
#from libbryn import *
from libHadronic import *
from libOneLepton import *
from icf.core import PSet,Analysis
from time import strftime
from batchGolden import *
from ra1objectid.vbtfElectronId_cff import *
from ra1objectid.vbtfMuonId_cff import *
from ra1objectid.ra3PhotonId_cff import *
from ra1objectid.ra3PhotonId2012_cff import *

vbtfMuonId_cff = Muon_IDFilter( vbtfmuonidps.ps()  )
vbtfElectronIdFilter = Electron_IDFilter( vbtfelectronidWP95ps.ps() )
ra3PhotonIdFilter    = Photon_IDFilter2012( ra3photonid2012ps.ps() )
CustomEleID = Electron_Egamma_Veto()
CustomMuID = OL_TightMuID(mu_2012_had.ps())
#  Change the settings from golden to use the lowest scaled bin.
default_common.Jets.PtCut=50.*(275./375.)
cutTree,blah,blah2,l = MakeDataTree(100.*(275./375.), Muon = None)

def addCutFlowData(a) :
  a.AddMuonFilter("PreCC",CustomMuID)
  a.AddPhotonFilter("PreCC",ra3PhotonIdFilter)
  a.AddElectronFilter("PreCC",CustomEleID)
  a+=cutTree

# AK5 Calo

conf_ak5_caloData = deepcopy(defaultConfig)
conf_ak5_caloData.Ntuple = deepcopy(ak5_calo)
conf_ak5_caloData.XCleaning = deepcopy(default_cc)
conf_ak5_caloData.Common = deepcopy(default_common)
# conf_ak5_calo.Common.print_out()
anal_ak5_caloData=Analysis("AK5Calo")
addCutFlowData(anal_ak5_caloData)


from data.Run2012.FNAL.SingleMu_Run2012C_PromptReco_v2_V17_0_taus_0_doTypeIMetReco_1_zmengJob287 import *
SingleMu_Run2012C_PromptReco_v2_V17_0_taus_0_doTypeIMetReco_1_zmengJob287.File=SingleMu_Run2012C_PromptReco_v2_V17_0_taus_0_doTypeIMetReco_1_zmengJob287.File[0:1]

from data.Run2012.FNAL.DoubleMu_Run2012C_PromptReco_v1_V17_0_taus_0_doTypeIMetReco_1_yeshaqJob312 import *
from data.Run2012.FNAL.DoubleMu_Run2012C_PromptReco_v2_V17_0_taus_0_doTypeIMetReco_1_yeshaqJob312 import *

DoubleMu_2012C = [DoubleMu_Run2012C_PromptReco_v1_V17_0_taus_0_doTypeIMetReco_1_yeshaqJob312, DoubleMu_Run2012C_PromptReco_v2_V17_0_taus_0_doTypeIMetReco_1_yeshaqJob312]

testData = DoubleMu_2012C

#outDir = "../Split_Jsons_"+strftime("%d_%b")+"/Data37/"
outDir = "../../results_"+strftime("%d_%b")+"//Data37"
ensure_dir(outDir)
anal_ak5_caloData.Run(outDir,conf_ak5_caloData,testData)

