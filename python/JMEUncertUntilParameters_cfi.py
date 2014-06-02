import FWCore.ParameterSet.Config as cms

defaultJMEUncertUntilParameters = cms.PSet(
    FilenameJECAK5MC = cms.untracked.string('/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_13_patch3_Bpbh/src/BpbH/BprimeTobHAnalysis/test/Summer13_V4_Uncertainties/Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt'), 
    FilenameJECCA8MC = cms.untracked.string('/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_13_patch3_Bpbh/src/BpbH/BprimeTobHAnalysis/test/Summer13_V4_Uncertainties/Summer13_V4_DATA_UncertaintySources_AK7PFchs.txt'), 
    #https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    JEREta = cms.untracked.vdouble(0.5,1.1,1.7,2.3,5.0),
    JERNominal = cms.untracked.vdouble(1.052,1.057,1.096,1.134,1.288),
    JERSigmaSym = cms.untracked.vdouble(0.012,0.012,0.017,0.035,0.127),
    JERSigmaNeg = cms.untracked.vdouble(0.061,0.055,0.062,0.085,0.153),
    JERSigmaPos = cms.untracked.vdouble(0.062,0.056,0.063,0.087,0.155),
    )
