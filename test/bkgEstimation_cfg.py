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
options.register('fatJetMassMin', 0.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum fat jet mass"
    )
options.register('fatJetMassMax', 1.E6,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum fat jet mass"
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
options.register('ApplyJEC', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Apply JEC" 
    )
options.register('ApplyBTagSF', False,
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
if options.SFbShift != 0.0 and options.SFlShift != 0.0: 
  print "SFbshift = ",  options.SFbShift, " and SFlshift = ", options.SFlShift
  print "Warning: must be varied independently."

options.setDefault('maxEvents', -1000) 

options.parseArguments()

process = cms.Process("ABCD")

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

from BpbH.BprimeTobH.TriggerSelector_cfi import * 
from BpbH.BprimeTobH.HiggsJetSelector_cfi import * 
from BpbH.BprimeTobH.HTSelector_cfi import * 
from BpbH.BprimeTobHAnalysis.EventSelector_cfi import * 
from BpbH.BprimeTobHAnalysis.JMEUncertUntilParameters_cfi import * 


process.ABCD = cms.EDAnalyzer('BackgroundEstimationABCD',
    MaxEvents           = cms.int32(options.maxEvents),
    ReportEvery         = cms.int32(options.reportEvery),  
    InputTTree          = cms.string(options.ttreedir+'/tree'),
    #InputFiles          = cms.vstring(FileNames), 
    InputFiles          = cms.vstring(FileNames_BpBp800), 
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
    Hist_PUDistMC       = cms.string('pileup_mc'),
    Hist_PUDistData     = cms.string('pileup_data'),

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
    Subjet1CSVDiscMin   = cms.double(0.679),
    Subjet1CSVDiscMax   = cms.double(2),
    Subjet2CSVDiscMin   = cms.double(0.679),
    Subjet2CSVDiscMax   = cms.double(2),
    dRSubjetsMin 	= cms.double(0.3),
    dRSubjetsMax 	= cms.double(1000),

    HJetSBMassMin = cms.double(-1),	
    HJetSBMassMax = cms.double(80),
	
    JetPtMin		= cms.double(50),
    JetPtMax		= cms.double(100000),
    JetAbsEtaMin	= cms.double(-1),
    JetAbsEtaMax	= cms.double(2.4),
    BJetPtMin		= cms.double(80),
    BJetCSVDiscMin  	  = cms.double(options.bjetCSVDiscMin),
    BJetCSVDiscMax   	  = cms.double(options.bjetCSVDiscMax),
    FatJetPtMin         = cms.double(options.fatJetPtMin),
    FatJetPtMax         = cms.double(options.fatJetPtMax),
    FatJetMassMin       = cms.double(options.fatJetMassMin),
    FatJetMassMax       = cms.double(options.fatJetMassMax),
    FatJetPrunedMassMin = cms.double(options.fatJetPrunedMassMin),
    FatJetPrunedMassMax = cms.double(options.fatJetPrunedMassMax),

    bVetoJetPtMin	= cms.double(30),
    bVetoJetPtMax	= cms.double(100000),
    bVetoJetAbsEtaMin	= cms.double(-1),
    bVetoJetAbsEtaMax	= cms.double(2.4),
    bVetoJetCSVDiscMin 	= cms.double(0.244),
    bVetoJetCSVDiscMax 	= cms.double(2),

    numbJetMin		= cms.int32(1),
    numHiggsJetMin	= cms.int32(1),

    JetSelParams        = defaultJetSelectionParameters.clone(), 
    FatJetSelParams     = defaultFatJetSelectionParameters.clone(), 
    HiggsJetSelParams   = defaultHiggsJetSelectionParameters.clone(
      subjet1CSVDiscMin = cms.double(0.679),
      subjet2CSVDiscMin = cms.double(0.679),
      ), 
    HTSelParams         = defaultHTSelectionParameters.clone(
      HTMin = cms.double(900), 
      ),
    EvtSelParams        = defaultEventSelectionParameters.clone(),
    JMEParams           = defaultJMEUncertUntilParameters.clone(), 
    ApplyJEC            = cms.bool(options.ApplyJEC), 
    ApplyBTagSF         = cms.bool(options.ApplyBTagSF), 
    JESShift            = cms.double(options.JESShift), 
    JERShift            = cms.double(options.JERShift), 
    SFbShift            = cms.double(options.SFbShift), 
    SFlShift            = cms.double(options.SFlShift), 
    BuildMinTree        = cms.bool(True),
    ) 

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )
process.p = cms.Path(process.ABCD)

