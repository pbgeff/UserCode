#
#  SUSY-PAT configuration file
#
#  PAT configuration for the SUSY group - 52X series
#  More information here:
#  https://twiki.cern.ch/twiki/bin/view/CMS/SusyPatLayer1DefV12
#

#switch between MC and data
#switch between fastsim and fullsim/data

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("PAT")
options = VarParsing ('standard')

globalTags = {"PromptReco": "GR_P_V40_AN1::All", "Aug24ReReco": "FT_53_V10_AN2::All" , "Aug06ReReco": "FT_53_V6C_AN2::All", "July13ReReco": "FT_53_V6_AN2::All",
	      "pre526FastSim": "START53_V7F::All", "526andLaterFastSim": "START53_V7F::All", "MC": "START53_V7F::All"}
## ************************************************************
## The following line must be configured properly
## datasetType must be one of the options above in "globalTags"
## ************************************************************
datasetType = "MC"

## options for testing
options.output='configurableAnalysis.root'
options.files='file:/cmsdata/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_PU_S10_START53_V7A_AODSIM/ECDEFDB7-AAE1-E111-B576-003048C68A88.root'
maxEvents=10

## determine if we are running on an MC dataset
isMC = False
if datasetType=="MC" or datasetType=="pre526FastSim" or datasetType=="526andLaterFastSim":
	isMC = True

## Print out the cfA configuration information
print "Using global tag " + globalTags[datasetType] + " selected from datasetType=" + datasetType
print "isMC: " + str(isMC)

#-- Message Logger ------------------------------------------------------------
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
    limit = cms.untracked.int32(-1),
    reportEvery = cms.untracked.int32(100)
)
process.MessageLogger.suppressError = cms.untracked.vstring('patElectronsPF', 'patMuonsPF')

#-- Source information ------------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.files)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(maxEvents) )
#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(
#  '190645:10-190645:110',
#)
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.Geometry.GeometryIdeal_cff")
## The geometry sequence now generates a deprecation warning
## and so has been replaced by the one above
#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.BFieldColl = cms.EDProducer('BFieldProducer')
process.JetCorrColl = cms.EDProducer('JetCorrProducer')

#Need this for L1 triggers with CMSSW >= 381
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
process.patTrigger.addL1Algos = cms.bool( True )

## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
     #verbose = cms.untracked.bool(True),
     fileName = cms.untracked.string('patTuple.root'),
     SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
     outputCommands = cms.untracked.vstring('drop *', "keep *_BFieldColl_*_*_","keep *_JetCorrColl_*_*_", *patEventContent 
		 )
)

#-- SUSYPAT and GlobalTag Settings -----------------------------------------------------------
from PhysicsTools.Configuration.SUSY_pattuple_cff import addDefaultSUSYPAT, getSUSY_pattuple_outputCommands

process.GlobalTag.globaltag = globalTags[datasetType]

if isMC:
	addDefaultSUSYPAT(process,True,'HLT',['L1FastJet','L2Relative','L3Absolute'],'',['AK5PF'])
else:
	addDefaultSUSYPAT(process,False,'HLT',['L1FastJet','L2Relative','L3Absolute','L2L3Residual'],'',['AK5PF'])

process.pfNoTauPF.enable = cms.bool(False)
SUSY_pattuple_outputCommands = getSUSY_pattuple_outputCommands( process )

# make loose clones of the original electron collection
process.pfRelaxedElectronsPF = process.pfIsolatedElectronsPF.clone()
process.pfRelaxedElectronsPF.isolationCut = 9999.0
process.patElectronsPF.pfElectronSource  = "pfRelaxedElectronsPF"
process.pfElectronSequencePF.replace(process.pfIsolatedElectronsPF,
                                     process.pfIsolatedElectronsPF +
                                     process.pfRelaxedElectronsPF)

# make loose clones of the original muon collection
process.pfRelaxedMuonsPF = process.pfIsolatedMuonsPF.clone()
process.pfRelaxedMuonsPF.isolationCut = 9999.0
process.patMuonsPF.pfMuonSource  = "pfRelaxedMuonsPF"
process.pfMuonSequencePF.replace(process.pfIsolatedMuonsPF,
                                 process.pfIsolatedMuonsPF +
                                 process.pfRelaxedMuonsPF)

############################## END SUSYPAT specifics ####################################


################### Add Type-I PFMET (for default RECO-PF jets) ########################
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")

# NOTE: use "ak5PFL1FastL2L3" for MC / "ak5PFL1FastL2L3Residual" for Data
if isMC:
	process.pfJetMETcorr.jetCorrLabel = "ak5PFL1FastL2L3"
	process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_mc
else:
	process.pfJetMETcorr.jetCorrLabel = "ak5PFL1FastL2L3Residual"
	process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_data

##fix for Fastsim
if datasetType=="pre526FastSim":
	process.pfCandsNotInJet.bottomCollection = cms.InputTag("FSparticleFlow")
	process.pfJetMETcorr.skipMuons = cms.bool(False)


## generate typeI corrected pfMET (no Type0 correction)
process.patPFMETs = process.patMETs.clone(
    metSource = cms.InputTag('pfMet'),
    addMuonCorrections = cms.bool(False),
    #genMETSource = cms.InputTag('genMetTrue'),
    #addGenMET = cms.bool(True)
    )

process.pfType1CorrectedMet.applyType0Corrections = cms.bool(False)

process.patPFMETsTypeIcorrected = process.patPFMETs.clone(
    metSource = cms.InputTag('pfType1CorrectedMet'),
    )

## generate typeI corrected pfMET (with Type0 correction)
process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
process.pfType1Type0PFCandidateCorrectedMet = process.pfType1CorrectedMet.clone(
	applyType0Corrections = cms.bool(False),
	srcType1Corrections = cms.VInputTag(
	cms.InputTag('pfMETcorrType0'),
	cms.InputTag('pfJetMETcorr', 'type1')        
	)
    )

process.patPFMETsTypeIType0PFCandcorrected = process.patPFMETs.clone(
    metSource = cms.InputTag('pfType1Type0PFCandidateCorrectedMet'),
   )

## generate typeI corrected MET with x/y shift
process.pfType1XYCorrectedMet = process.pfType1CorrectedMet.clone(
	applyType0Corrections = cms.bool(False),
	srcType1Corrections = cms.VInputTag(
	cms.InputTag('pfJetMETcorr', 'type1'),        
	cms.InputTag('pfMEtSysShiftCorr')
	)
   )

process.patPFMETsTypeIXYcorrected = process.patPFMETs.clone(
    metSource = cms.InputTag('pfType1XYCorrectedMet'),
   )

#Turn on trigger info
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
switchOnTrigger(process, triggerProducer='patTrigger', triggerEventProducer='patTriggerEvent', sequence='patDefaultSequence', hltProcess="HLT")

process.load("Workspace.ConfigurableAnalysis.configurableAnalysis_ForPattuple_cff")
process.TFileService.fileName = cms.string(options.output)
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')
process.load('RecoMET.METAnalyzers.CSCHaloFilter_cfi')
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
## For AOD and RECO recommendation to use recovered rechits
process.EcalDeadCellTriggerPrimitiveFilter.tpDigiCollection = cms.InputTag("ecalTPSkimNA")
process.load('RecoMET.METFilters.inconsistentMuonPFCandidateFilter_cfi')
process.load('RecoMET.METFilters.greedyMuonPFCandidateFilter_cfi')
process.load('RecoMET.METFilters.eeNoiseFilter_cfi')
process.load('MyAnalyzers.TriggerFilter.triggerFilter_cfi')

process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                     applyfilter = cms.untracked.bool(True),
                                     debugOn = cms.untracked.bool(False),
                                     numtrack = cms.untracked.uint32(10),
                                     thresh = cms.untracked.double(0.25)
)

# The section below is for the filter on Boundary Energy. Available in AOD in CMSSW>44x
process.load('RecoMET.METFilters.EcalDeadCellBoundaryEnergyFilter_cfi')
process.EcalDeadCellBoundaryEnergyFilter.taggingMode = cms.bool(False)
process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyDeadCellsEB=cms.untracked.double(10)
process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyDeadCellsEE=cms.untracked.double(10)
process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyGapEB=cms.untracked.double(100)
process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyGapEE=cms.untracked.double(100)
process.EcalDeadCellBoundaryEnergyFilter.enableGap=cms.untracked.bool(False)
process.EcalDeadCellBoundaryEnergyFilter.limitDeadCellToChannelStatusEB = cms.vint32(12,14)
process.EcalDeadCellBoundaryEnergyFilter.limitDeadCellToChannelStatusEE = cms.vint32(12,14)
# End of Boundary Energy filter configuration 

## The HCAL laser filter _____________________________________________________||
process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
process.hcalLaserEventFilter.vetoByRunEventNumber=cms.untracked.bool(False)
process.hcalLaserEventFilter.vetoByHBHEOccupancy=cms.untracked.bool(True)

## The EE bad SuperCrystal filter ____________________________________________||
process.load('RecoMET.METFilters.eeBadScFilter_cfi')

#compute rho for 2011 effective area Egamma isolation corrections
from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation2011 = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation2011.Rho_EtaMax = cms.double(2.5)
#compute rho for 2012 effective area Egamma isolation corrections
process.kt6PFJetsForIsolation2012 = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation2012.Rho_EtaMax = cms.double(4.4)
process.kt6PFJetsForIsolation2012.voronoiRfact = cms.double(0.9)

#-- Output module configuration -----------------------------------------------
process.out.fileName = "SUSYPAT.root" 
process.out.splitLevel = cms.untracked.int32(99)  # Turn on split level (smaller files???)
process.out.overrideInputFileSplitLevels = cms.untracked.bool(True)
process.out.dropMetaData = cms.untracked.string('DROPPED')   # Get rid of metadata related to dropped collections
process.out.outputCommands = cms.untracked.vstring('keep *',"keep *_HBHENoiseFilterResultProducer_*_*","keep *_BFieldColl_*_*","keep *_JetCorrectionColl_*_*", *SUSY_pattuple_outputCommands )
process.out.outputCommands.append('keep *_patPFMETsTypeIcorrected_*_PAT')

#-- Execution path ------------------------------------------------------------
# Full path
#This is to run on full sim or data
process.hcallaserfilter = cms.Path(process.hcalLaserEventFilter)
process.eebadscfilter = cms.Path(process.eeBadScFilter)
process.ecaltpfilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
#Run both TP and BE filters; doesn't work right now
#process.ecaltpfilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter*process.EcalDeadCellBoundaryEnergyFilter)
process.scrapingveto = cms.Path(process.scrapingVeto)
process.greedymuonfilter = cms.Path(process.greedyMuonPFCandidateFilter)
process.inconsistentPFmuonfilter = cms.Path(process.inconsistentMuonPFCandidateFilter)
process.eenoisefilter = cms.Path(process.eeNoiseFilter)
process.passprescalePFHT350filter = cms.Path( process.pfht350PassPrescaleFilter )
## Adding more HT "active trigger" variables
process.triggerFilterHT250 = process.pfht350PassPrescaleFilter.clone(
    HLTPaths = cms.vstring('HLT_HT250_v[0-9]')
    )
process.passprescaleHT250filter = cms.Path( process.triggerFilterHT250 )

process.triggerFilterHT300 = process.pfht350PassPrescaleFilter.clone(
    HLTPaths = cms.vstring('HLT_HT300_v[0-9]')
    )
process.passprescaleHT300filter = cms.Path( process.triggerFilterHT300 )

process.triggerFilterHT350 = process.pfht350PassPrescaleFilter.clone(
    HLTPaths = cms.vstring('HLT_HT350_v[0-9]')
    )
process.passprescaleHT350filter = cms.Path( process.triggerFilterHT350 )

process.triggerFilterHT400 = process.pfht350PassPrescaleFilter.clone(
    HLTPaths = cms.vstring('HLT_HT400_v[0-9]')
    )
process.passprescaleHT400filter = cms.Path( process.triggerFilterHT400 )

process.triggerFilterHT450 = process.pfht350PassPrescaleFilter.clone(
    HLTPaths = cms.vstring('HLT_HT450_v[0-9]')
    )
process.passprescaleHT450filter = cms.Path( process.triggerFilterHT450 )

if datasetType=="526andLaterFastSim" or datasetType=="pre526FastSim":
        process.p = cms.Path( process.goodOfflinePrimaryVertices + process.BFieldColl + process.susyPatDefaultSequence + process.JetCorrColl)
	process.p += process.kt6PFJetsForIsolation2012
else:
	process.csctighthalofilter = cms.Path(process.CSCTightHaloFilter)
        process.p = cms.Path(process.goodOfflinePrimaryVertices + process.HBHENoiseFilterResultProducer + process.BFieldColl + process.susyPatDefaultSequence + process.JetCorrColl)

process.p += process.kt6PFJetsForIsolation2011

process.p += process.pfMEtSysShiftCorrSequence
process.p += process.type0PFMEtCorrection
process.p += process.producePFMETCorrections

process.p += process.patPFMETsTypeIcorrected
process.p += process.pfType1Type0PFCandidateCorrectedMet
process.p += process.patPFMETsTypeIType0PFCandcorrected
process.p += process.pfType1XYCorrectedMet
process.p += process.patPFMETsTypeIXYcorrected


process.trackingfailturefilter = cms.Path(process.trackingFailureFilter)
process.outpath = cms.EndPath(cms.ignore(process.configurableAnalysis))
#process.outpath = cms.EndPath(process.out) #Output the SUSYPAT.root file

#-- Dump config ------------------------------------------------------------
file = open('SusyPAT_cfg.py','w')
file.write(str(process.dumpPython()))
file.close()
