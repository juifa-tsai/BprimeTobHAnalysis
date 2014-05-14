import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

from inputFiles_cfi import * 

options = VarParsing('python')

options.register('outFilename', 'bprimeTobH_PDFSyst.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
    )
options.register('reportEvery', 1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=1000)"
    )
options.register('ttreedir', 'BprimebH',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Name of ROOT TTree dir: Either 'ntuple' or 'skim' or 'bVeto'"
    )

options.setDefault('maxEvents', -1000) 

options.parseArguments()

process = cms.Process("PDFSyst")

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

FileNames = ['bprimeTobHAnalysis_PDFtree.root']

process.PDFSyst = cms.EDAnalyzer('PdfSystematicsAnalyzer',
    MaxEvents           = cms.int32(options.maxEvents),
    ReportEvery         = cms.int32(options.reportEvery),  
    InputTTree          = cms.string(options.ttreedir+'/pdftree'),
    InputFiles          = cms.vstring(FileNames), 
    ) 

process.p = cms.Path(process.PDFSyst)

