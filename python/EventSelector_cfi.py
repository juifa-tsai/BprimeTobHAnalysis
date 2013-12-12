import FWCore.ParameterSet.Config as cms

from BpbH.BprimeTobH.TriggerSelector_cfi import *
from BpbH.BprimeTobH.JetSelector_cfi import *
from BpbH.BprimeTobH.BJetSelector_cfi import *
from BpbH.BprimeTobH.FatJetSelector_cfi import * 
from BpbH.BprimeTobH.HiggsJetSelector_cfi import * 
from BpbH.BprimeTobH.HTSelector_cfi import * 

defaultEventSelectionParameters = cms.PSet (
    CutLevels                = cms.vstring('Trigger', 'Vertex', 
      'MinBjets', 'MaxBjets', 'MinHiggsjets', 'MaxHiggsjets', 'HT'),  
    RemoveJetHiggsJetOverlap = cms.bool(True), 
    MinNJets                 = cms.int32(0),  
    MaxNJets                 = cms.int32(10000),  
    MinNBJets                = cms.int32(2),  
    MaxNBJets                = cms.int32(10000),  
    MinNFatJets              = cms.int32(0),  
    MaxNFatJets              = cms.int32(10000),  
    MinNHiggsJets            = cms.int32(1),  
    MaxNHiggsJets            = cms.int32(10000),   
    TrigSelParams            = defaultTriggerSelectionParameters.clone(), 
    JetSelParams             = defaultJetSelectionParameters.clone(), 
    BJetSelParams            = defaultBJetSelectionParameters.clone(), 
    FatJetSelParams          = defaultFatJetSelectionParameters.clone(), 
    HiggsJetSelParams        = defaultHiggsJetSelectionParameters.clone(), 
    HTSelParams              = defaultHTSelectionParameters.clone(),
    )
