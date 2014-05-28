import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

from inputFiles_cfi import * 

options = VarParsing('python')

options.register('outFilename', 'bprimeTobH_BkgABCDMethod.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
    )
options.register('reportEvery', 1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=1000)"
    )
options.register('ttreedir', 'ntuple',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Name of ROOT TTree dir: Either 'ntuple' or 'skim' or 'bVeto'"
    )
options.register('doPUReweighting', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Do pileup reweighting"
)
options.register('pileupmchist', 'pileup_mc',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Name of Histogram for MC pileup weights"
    )
options.register('pileupdatahist', 'pileup_data',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Name of Histogram for data pileup weights"
    )
options.register('JetPtMin', 50.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum b jet Pt"
    )
options.register('bJetPtMin', 80.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum b jet Pt"
    )
options.register('bjetCSVDiscMin', 0.679,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum b jet CSV discriminator"
    )
options.register('bjetCSVDiscMax', 1.000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum b jet CSV discriminator"
    )
options.register('fatJetPtMin', 300.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum fat jet pt"
    )
options.register('fatJetPtMax', 1.E6,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum fat jet pt"
    )
options.register('fatJetPrunedMassMin', 90.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum fat jet pruned mass"
    )
options.register('fatJetPrunedMassMax', 140.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum fat jet pruned mass"
    )
options.register('dRSubjetsMin', 0.3,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum dR(subjet1, subjet2)"
    )
options.register('dRSubjetsMax', 999,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum dR(subjet1, subjet2)"
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
options.register('ApplyJEC', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Apply JEC" 
    )
options.register('ApplyBTagSF', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Apply b-tagging scale factors" 
    )
options.register('JESShift', 0.0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "JES shift in unit of sigmas" 
    )
options.register('JERShift', 0.0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "JER shift in unit of sigmas" 
    )
options.register('SFbShiftHtag', 0.0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "SFb shift (for Higgs-tagging) in unit of sigmas" 
    )
options.register('SFlShiftHtag', 0.0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "SFl shift (for Higgs-tagging) in unit of sigmas" 
    )
options.register('SFbShift', 0.0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "SFb shift in unit of sigmas" 
    )
options.register('SFlShift', 0.0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "SFl shift in unit of sigmas" 
    )

if options.SFbShiftHtag != 0.0 and options.SFlShiftHtag != 0.0: 
  print "SFbshiftHtag = ",  options.SFbShiftHtag, " and SFlshiftHtag = ", options.SFlShiftHtag
  print "Warning: must be varied independently."

if options.SFbShift != 0.0 and options.SFlShift != 0.0: 
  print "SFbshift = ",  options.SFbShift, " and SFlshift = ", options.SFlShift
  print "Warning: must be varied independently."

options.setDefault('maxEvents', -1000) 

options.parseArguments()

process = cms.Process("ABCD")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring(
      'ABCD',
      ),
    ABCD = cms.untracked.PSet(
      threshold = cms.untracked.string('INFO'),  
      ), 
    suppressInfo = cms.untracked.vstring('ABCD'),
    ) 
#process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) ) # Leave it this way. 

process.source = cms.Source("EmptySource")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outFilename) 
    )

from BpbH.BprimeTobH.TriggerSelector_cfi import * 
from BpbH.BprimeTobH.HiggsJetSelector_cfi import * 
from BpbH.BprimeTobH.HTSelector_cfi import * 
from BpbH.BprimeTobHAnalysis.EventSelector_cfi import * 
from BpbH.BprimeTobHAnalysis.JMEUncertUntilParameters_cfi import * 


process.ABCD = cms.EDAnalyzer('BackgroundEstimationABCD',
    MaxEvents           = cms.int32(options.maxEvents),
    ReportEvery         = cms.int32(options.reportEvery),  
    InputTTree          = cms.string(options.ttreedir+'/tree'),
    InputFiles          = cms.vstring(FileNames), 
    #InputFiles          = cms.vstring(FileNames_TTbar), 
    #InputFiles          = cms.vstring(FileNames_BpBp800), 
    #InputFiles          = cms.vstring(SkimmedFileNames_BpBp500), 
    #InputFiles          = cms.vstring(SkimmedFileNames_BpBp1000), 
    #InputFiles          = cms.vstring(SkimmedFileNames_QCD300to470), 
    #InputFiles          = cms.vstring(SkimmedFileNames_QCD1800), 
    #InputFiles          = cms.vstring(FileNames_QCD_HT_100To250), 
    #InputFiles          = cms.vstring(FileNames_QCD_HT_250To500), 
    #InputFiles          = cms.vstring(FileNames_QCD_HT_500To1000), 
    #InputFiles          = cms.vstring(FileNames_QCD_HT_1000ToInf), 
    HLTPaths            = defaultTriggerSelectionParameters.clone(), 
    DoPUReweighting     = cms.bool(options.doPUReweighting),
    File_PUDistMC       = cms.string('pileup_Data_Summer12_53X_S10.root'),
    File_PUDistData     = cms.string('pileup_Data_Summer12_53X_S10.root'),
    Hist_PUDistMC       = cms.string(options.pileupmchist),
    Hist_PUDistData     = cms.string(options.pileupdatahist),
    HTAK5Min		= cms.double(900),
    HTAK5Max		= cms.double(100000),
    HJetPtMin		= cms.double(300),
    HJetPtMax		= cms.double(100000),
    HJetAbsEtaMax	= cms.double(2.4),
    HJetMassMin 	= cms.double(100),	
    HJetMassMax 	= cms.double(150),	
    HJetPrunedMassMin = cms.double(90),	
    HJetPrunedMassMax = cms.double(140),	
    HJetTau2ByTau1Min = cms.double(-1),
    HJetTau2ByTau1Max = cms.double(0.5),
    DRSubjetsMin        = cms.double(options.dRSubjetsMin),
    DRSubjetsMax        = cms.double(options.dRSubjetsMax),
    Subjet1CSVDiscMin   = cms.double(options.subjet1CSVDiscMin),
    Subjet1CSVDiscMax   = cms.double(options.subjet1CSVDiscMax),
    Subjet2CSVDiscMin   = cms.double(options.subjet2CSVDiscMin),
    Subjet2CSVDiscMax   = cms.double(options.subjet2CSVDiscMax),
    HJetSBMassMin       = cms.double(-1),	
    HJetSBMassMax       = cms.double(80),
    JetPtMin           = cms.double(options.JetPtMin),
    BJetPtMin           = cms.double(options.bJetPtMin),
    BJetCSVDiscMin  	  = cms.double(options.bjetCSVDiscMin),
    BJetCSVDiscMax   	  = cms.double(options.bjetCSVDiscMax),
    FatJetPtMin         = cms.double(options.fatJetPtMin),
    FatJetPtMax         = cms.double(options.fatJetPtMax),
    FatJetPrunedMassMin = cms.double(options.fatJetPrunedMassMin),
    FatJetPrunedMassMax = cms.double(options.fatJetPrunedMassMax),
    bVetoJetPtMin	      = cms.double(30),
    bVetoJetPtMax	      = cms.double(100000),
    bVetoJetAbsEtaMin	  = cms.double(0.0),
    bVetoJetAbsEtaMax	  = cms.double(2.4),
    bVetoJetCSVDiscMin 	= cms.double(0.244),
    bVetoJetCSVDiscMax 	= cms.double(2),
    numbJetMin		      = cms.int32(1),
    numHiggsJetMin    	= cms.int32(1),
    JetSelParams        = defaultJetSelectionParameters.clone(
      jetPtMin            = cms.double(30),
      ), 
    FatJetSelParams     = defaultFatJetSelectionParameters.clone(), 
    HiggsJetSelParams   = defaultHiggsJetSelectionParameters.clone(
      subjet1CSVDiscMin = cms.double(0.679),
      subjet2CSVDiscMin = cms.double(0.679),
      ), 
    HTSelParams         = defaultHTSelectionParameters.clone(
      HTMin = cms.double(950), 
      ),
    EvtSelParams        = defaultEventSelectionParameters.clone(),
    JMEParams           = defaultJMEUncertUntilParameters.clone(), 
    ApplyJEC            = cms.bool(options.ApplyJEC), 
    ApplyBTagSF         = cms.bool(options.ApplyBTagSF), 
    JESShift            = cms.double(options.JESShift), 
    JERShift            = cms.double(options.JERShift), 
    SFbShiftHtag        = cms.double(options.SFbShiftHtag), 
    SFlShiftHtag        = cms.double(options.SFlShiftHtag), 
    SFbShift            = cms.double(options.SFbShift), 
    SFlShift            = cms.double(options.SFlShift), 
    BuildMinTree        = cms.bool(True),
    ) 

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )
process.p = cms.Path(process.ABCD)

