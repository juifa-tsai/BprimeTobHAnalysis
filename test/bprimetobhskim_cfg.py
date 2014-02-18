import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

from inputFiles_cfi import * 

options = VarParsing('python')

options.register('outFilename', 'evtSkim.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
    )
options.register('reportEvery', 1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=1000)"
    )
options.setDefault('maxEvents', -1000) 

options.parseArguments()

process = cms.Process("BprimebH")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger = cms.Service("MessageLogger",
#    destinations  =  cms.untracked.vstring('debugmessages', 'cout'
#      ), 
#    categories    = cms.untracked.vstring('BprimebHSkim'),
#    debugModules  = cms.untracked.vstring('*'),
#    debugmessages = cms.untracked.PSet(
#      threshold =  cms.untracked.string('DEBUG'),
#      ) 
#    ) 
process.MessageLogger.cout = cms.untracked.PSet(
    threshold =  cms.untracked.string('INFO'),
    ) 
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) ) # Leave it this way. 

process.source = cms.Source("EmptySource")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outFilename) 
    )

process.skim = cms.EDAnalyzer('BprimebHSkim',
    MaxEvents            = cms.int32(options.maxEvents),
    ReportEvery          = cms.int32(options.reportEvery),  
    InputTTree           = cms.string('ntuple/tree'),
    #InputFiles           = cms.vstring(FileNames), 
    InputFiles          = cms.vstring(FileNamesTTbar), 
    ) 

process.p = cms.Path(process.skim)
