import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

from BprimeTobHAnalysisv1.BprimeTobHAnalysis.BpBpToBHBHinc.BprimeBprimeTobHbHinc_M_800_cfi import * 
#from BprimeTobHAnalysisv1.BprimeTobHAnalysis.Data.JetHT_Run2012BCD_cfi import * 

options = VarParsing('python')

options.register('outFilename', 'bprimeTobH.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
    )
options.register('reportEvery', 1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=1000)"
    )
options.register('jetPtMin', 30.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum jet Pt"
    )
options.register('jetPtMax', 1.E6,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum jet Pt"
    )
options.register('bJetPtMin', 100.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum b jet Pt"
    )
options.register('fatJetPtMin', 300.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum fat jet Pt"
    )
options.register('fatJetPtMax', 1.E6,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum fat jet Pt"
    )
options.register('fatJetMassMin', 100.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum fat jet mass"
    )
options.register('fatJetMassMax', 150.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum fat jet mass"
    )
options.register('fatJetPrunedMassMin', 75.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum fat jet pruned mass"
    )
options.register('fatJetPrunedMassMax', 1.E6,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum fat jet pruned mass"
    )
options.register('fatJetTau2ByTau1Min', 0.5,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum fat jet tau2/tau1"
    )
options.register('subjet1CSVDiscMin', 0.679,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum subjet1 b discriminator"
    )
options.register('subjet1CSVDiscMax', 1.000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum subjet1 b discriminator"
    )
options.register('subjet2CSVDiscMin', 0.679,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum subjet2 b discriminator"
    )
options.register('subjet2CSVDiscMax', 1.000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum subjet2 b discriminator"
    )
options.register('hTMin', 1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum HT"
    )
options.register('hTMax', 1.E6,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum HT"
    )

options.setDefault('maxEvents', -50000) 

options.parseArguments()

process = cms.Process("BprimebH")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cout = cms.untracked.PSet(
    threshold = cms.untracked.string('INFO'), 
    ) 
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) ) # Leave it this way. 

process.source = cms.Source("EmptySource")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outFilename) 
    )

process.BprimebH = cms.EDAnalyzer('BprimeTobHAnalysis',
    MaxEvents           = cms.int32(options.maxEvents),
    ReportEvery         = cms.int32(options.reportEvery),  
    InputTTree          = cms.string('ntuple/tree'),
    InputFiles          = cms.vstring(FileNames), 
    HLTPaths            = cms.vint32(3225,4136,4137,5089,5537,5538), 
    JetPtMin            = cms.double(options.jetPtMin),
    JetPtMax            = cms.double(options.jetPtMax),
    JetAbsEtaMax        = cms.double(2.4),
    BJetPtMin           = cms.double(options.bJetPtMin),
    FatJetPtMin         = cms.double(options.fatJetPtMin),
    FatJetPtMax         = cms.double(options.fatJetPtMax),
    FatJetAbsEtaMax     = cms.double(2.4),
    FatJetMassMin       = cms.double(options.fatJetMassMin),
    FatJetMassMax       = cms.double(options.fatJetMassMax),
    FatJetPrunedMassMin = cms.double(options.fatJetPrunedMassMin),
    FatJetPrunedMassMax = cms.double(options.fatJetPrunedMassMax),
    FatJetTau2ByTau1Min = cms.double(options.fatJetTau2ByTau1Min),
    Subjet1CSVDiscMin   = cms.double(options.subjet1CSVDiscMin),
    Subjet1CSVDiscMax   = cms.double(options.subjet1CSVDiscMax),
    Subjet2CSVDiscMin   = cms.double(options.subjet2CSVDiscMin),
    Subjet2CSVDiscMax   = cms.double(options.subjet2CSVDiscMax),
    HTMin               = cms.double(options.hTMin),
    HTMax               = cms.double(options.hTMax), 
    ) 

process.p = cms.Path(process.BprimebH)

