#!/usr/bin/env python
import setupSUSY
from libFrameworkSUSY import *

from icf.core import PSet,Analysis

# This is the default configuration
from icf.config import defaultConfig

# This imports our demo sample
import samples.icf_demo as icf

import samples.QCD_EMBCE_PSet as qcd
# import default configuration
conf=defaultConfig

filter_skinny=PSet(
    OutputDir="/tmp/as1604/",
    Branches = [
        " drop * ",
        " keep beamSpotPosition ",
        " keep hltHandleValid ",
        " keep triggered ",
        " keep metnohfP4Calo ",
        " keep metP4TC ",
        " keep metP4PF ",
        " keep metP4IC5 ",
        " keep metP4Calo ",
        " keep metP4AK5 ",
        " keep ak5JetGenJetMatchIndexPat ",
        " keep ak5JetGenJetP4Pat ",
        " keep genPdgId ",
        " keep genMother ",
        " keep genMotherStored ",
        " keep genHasMother ",
        " keep genStatus ",
        " keep genP4 ",
        " keep genHandleValid ",
        " keep photonIDTightPat ",
        " keep photonIDLoosePat ",
        " keep photonCaloIsoPat ",
        " keep photonHcalIsoPat ",
        " keep photonEcalIsoPat ",
        " keep photonTrkIsoPat ",
        " keep photonSigmaIetaIetaPat ",
        " keep photonR9Pat ",
        " keep photonHasPixelSeedPat ",
        " keep photonHadronicOverEmPat ",
        " keep photonTrkSumPtHollowConeDR04Pat ",
        " keep photonHcalDepth2TowSumEtConeDR04Pat ",
        " keep photonHcalDepth1TowSumEtConeDR04Pat ",
        " keep photonEcalRecHitEtConeDR04Pat ",
        " keep photonP4Pat ",
        " keep ak5JetSimpleSecondaryVertexBJetTagsPat ",
        " keep ak5JetTrkCountingHighEffBJetTagsPat ",
        " keep ak5JetJetIDtightPat ",
        " keep ak5JetJetIDminimalPat ",
        " keep ak5JetJetIDloosePat ",
        " keep ak5JetEmEnergyFractionPat ",
        " keep ak5JetCorrFactorPat ",
        " keep ak5JetCorrectedP4Pat ",
        " keep tauTauIdbyTaNCfrTenthPercentPat ",
        " keep tauTauIdbyTaNCfrQuarterPercentPat ",
        " keep tauTauIdbyTaNCfrHalfPercentPat ",
        " keep tauTauIdbyTaNCfrOnePercentPat ",
        " keep tauTauIdbyTaNCPat ",
        " keep tauHcalIsoPat ",
        " keep tauEcalIsoPat ",
        " keep tauTrkIsoPat ",
        " keep tauChargePat ",
        " keep tauP4Pat ",
        " keep muonOuterTrackNumberOfValidHitsPat ",
        " keep muonInnerTrackNumberOfValidHitsPat ",
        " keep muonGlobalTracknumberOfValidHitsPat ",
        " keep muonInnerTrackNormalizedChi2Pat ",
        " keep muonGlobalTracknormalizedChi2Pat ",
        " keep muonGlobalTrackDxyErrorPat ",
        " keep muonGlobalTrackDxyPat ",
        " keep muonInnerTrackDzPat ",
        " keep muonInnerTrackDxyBSPat ",
        " keep muonInnerTrackDxyErrorPat ",
        " keep muonInnerTrackDxyPat ",
        " keep muonVertexNdofPat ",
        " keep muonVertexChi2Pat ",
        " keep muonIDGlobalMuonPromptTightPat ",
        " keep muonIsGlobalMuonPat ",
        " keep muonIsTrackerMuonPat ",
        " keep muonTMLastStationAngTightPat ",
        " keep muonPhotonIsoPat ",
        " keep muonNeutralHadronIsoPat ",
        " keep muonChargedHadronIsoPat ",
        " keep muonHcalIsoPat ",
        " keep muonEcalIsoPat ",
        " keep muonTrackIsoPat ",
        " keep muonChargePat ",
        " keep muonP4Pat ",
        " keep electronGsfTrackTrackerExpectedHitsInnerPat ",
        " keep electronConversionPartnerTrackTrackerExpectedHitsInnerPat ",
        " keep electronConversionDistPat ",
        " keep electronConversionDCotPat ",
        " keep electronGsfTrackDxyPat ",
        " keep electronFbremPat ",
        " keep electronESuperClusterOverPPat ",
        " keep electronHasValidHitInFirstPixelBarrelPat ",
        " keep electronSigmaIetaIetaPat ",
        " keep electronDeltaEtaSuperClusterTrackAtVtxPat ",
        " keep electronDeltaPhiSuperClusterTrackAtVtxPat ",
        " keep electronHcalOverEcalPat ",
        " keep electronMvaPat ",
        " keep electronScPixChargePat ",
        " keep electronClosestCtfTrackChargePat ",
        " keep electronEcalEnergyPat ",
        " keep electronGsfTrackPtPat ",
        " keep electronEIDTightPat ",
        " keep electronEIDLoosePat ",
        " keep electronPhotonIsoPat ",
        " keep electronNeutralHadronIsoPat ",
        " keep electronChargedHadronIsoPat ",
        " keep electronDr03HcalTowerSumEtPat ",
        " keep electronDr03EcalRecHitSumEtPat ",
        " keep electronDr03TkSumPtPat ",
        " keep electronHcalIsoPat ",
        " keep electronEcalIsoPat ",
        " keep electronTrackIsoPat ",
        " keep electronChargePat ",
        " keep electronP4Pat ",
        " keep physicsDeclared ",
        " keep vertexPosition ",
        " keep vertexIsFake ",
        " keep vertexNdof ",
        " keep tracksNEtaGT1p5HighPurityTracks ",
        " keep tracksNEta0p9to1p5HighPurityTracks ",
        " keep tracksNEtaLT0p9HighPurityTracks ",
        " keep tracksNEtaGT1p5AllTracks ",
        " keep tracksNEta0p9to1p5AllTracks ",
        " keep tracksNEtaLT0p9AllTracks ",
        " keep lumiSection ",
        " keep genpthat ",
        " keep event ",
        " keep run ",
        ]
)

# Create the analysis
a=Analysis("eWPolSlim",filter_skinny)
s=[
    qcd.bce_20_30,
]

a.Run("../results",conf,s)
