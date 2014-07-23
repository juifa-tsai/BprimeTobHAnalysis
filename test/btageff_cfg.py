import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

from inputFiles_cfi import * 

options = VarParsing('python')

options.register('outFilename', 'BtagEff.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
    )
options.register('reportEvery', 1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=1000)"
    )
options.register('jetPtMin', 50.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum jet Pt"
    )
options.register('jetPtMax', 1.E6,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum jet Pt"
    )
options.register('bJetPtMin', 80.,
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
options.register('fatJetMassMin', 0.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum fat jet mass"
    )
options.register('fatJetMassMax', 100000.,
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
options.register('fatJetTau2ByTau1Max', 0.5,
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
options.register('hTMin', 900,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum HT"
    )
options.register('hTMax', 1.E6,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Maximum HT"
    )

options.register('doPUReweighting', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Do pileup reweighting"
    )

options.setDefault('maxEvents', -1000) 

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

from BpbH.BprimeTobHAnalysis.EventSelector_cfi import * 

process.BTagEff = cms.EDAnalyzer('BTagEff',
    MaxEvents            = cms.int32(options.maxEvents),
    ReportEvery          = cms.int32(options.reportEvery),  
    InputTTree           = cms.string('ntuple/tree'),
    #InputFiles           = cms.vstring(FileNames), 
    InputFiles           = cms.vstring(FileNames_BpBp500), 
    #InputFiles           = cms.vstring(FileNames_TTbar), 
    HLTPaths             = defaultTriggerSelectionParameters.clone(), 
    DoPUReweighting      = cms.bool(options.doPUReweighting),
    File_PUDistMC        = cms.string('pileup_Data_Summer12_53X_S10.root'),
    File_PUDistData      = cms.string('pileup_Data_Summer12_53X_S10.root'),
    Hist_PUDistMC        = cms.string('pileup_mc'),
    Hist_PUDistData      = cms.string('pileup_data'),
    EvtSelParams         = defaultEventSelectionParameters.clone(
      CutLevels          = cms.vstring('Trigger', 'Vertex', 'MinHiggsjets', 'MaxHiggsjets', 'MinNJets', 'MaxNJets', 'HT'),
      FatJetSelParams          = defaultFatJetSelectionParameters.clone(
        fatJetPtMin        = cms.double(300), 
        ), 
      HiggsJetSelParams        = defaultHiggsJetSelectionParameters.clone(
        fatJetPtMin        = cms.double(300), 
        ), 
      ), 
    FatJetSelParams     = defaultFatJetSelectionParameters.clone(
      jettype = cms.string('CA8JET'), 
      fatJetPtMin = cms.double(30.),
      fatJetTau2ByTau1Max = cms.double(1.),
      ), 
    FatBJetSelParams     = defaultBTaggedFatJetSelectionParameters.clone(
      fatJetPtMin = cms.double(30.),
      fatJetTau2ByTau1Max = cms.double(1.),
      fatJetCSVDiscMin    = cms.double(-2000), 
      fatJetCSVDiscMax    = cms.double(1.000), 
      subjet1CSVDiscMin   = cms.double(-2000),
      subjet1CSVDiscMax   = cms.double(1.000),
      subjet2CSVDiscMin   = cms.double(-2000),
      subjet2CSVDiscMax   = cms.double(1.000),
      ), 
    JetSelParams        = defaultJetSelectionParameters.clone(), 
    BJetSelParams       = defaultBJetSelectionParameters.clone(), 
    ) 

process.p = cms.Path(process.BTagEff)

