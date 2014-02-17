#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TString.h>
#include <iostream> 

TString dir = "/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_13_patch3_Bpbh/src/BpbH/BprimeTobHAnalysis/test/OnLxplus/LXBATCH_Jobs_18Dec_JetCorr/" ; 

void evtCounts () {

  TString fnames [27] = 
  {"BprimeBprimeToBHBHinc_M-1000_TuneZ2star_8TeV-madgraph__Summer12_DR53X-PU_S10_START53_V7C-v1__AODSIM", 
    "BprimeBprimeToBHBHinc_M-1200_TuneZ2star_8TeV-madgraph__Summer12_DR53X-PU_S10_START53_V7C-v1__AODSIM", 
    "BprimeBprimeToBHBHinc_M-1500_TuneZ2star_8TeV-madgraph__Summer12_DR53X-PU_S10_START53_V7C-v1__AODSIM", 
    "BprimeBprimeToBHBHinc_M-450_TuneZ2star_8TeV-madgraph__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM", 
    "BprimeBprimeToBHBHinc_M-500_TuneZ2star_8TeV-madgraph__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM", 
    "BprimeBprimeToBHBHinc_M-550_TuneZ2star_8TeV-madgraph__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM", 
    "BprimeBprimeToBHBHinc_M-600_TuneZ2star_8TeV-madgraph__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM", 
    "BprimeBprimeToBHBHinc_M-650_TuneZ2star_8TeV-madgraph__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM", 
    "BprimeBprimeToBHBHinc_M-700_TuneZ2star_8TeV-madgraph__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM", 
    "BprimeBprimeToBHBHinc_M-750_TuneZ2star_8TeV-madgraph__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM", 
    "BprimeBprimeToBHBHinc_M-800_TuneZ2star_8TeV-madgraph__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM", 
    "JetHT__Run2012B-22Jan2013-v1__AOD", 
    "JetHT__Run2012B_missingPart-22Jan2013-v1__AOD", 
    "JetHT__Run2012C-22Jan2013-v1__AOD", 
    "JetHT__Run2012D-22Jan2013-v1__AOD", 
    "JetHT__Run2012D_missingPart-22Jan2013-v1__AOD", 
    "Jet__Run2012A-22Jan2013-v1__AOD", 
    "Jet__Run2012A_missingPart-22Jan2013-v1__AOD", 
    "QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM", 
    "QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM", 
    "QCD_Pt-170to300_TuneZ2star_8TeV_pythia6__Summer12_DR53X-PU_S10_START53_V7A-v2__AODSIM", 
    "QCD_Pt-1800_TuneZ2star_8TeV_pythia6__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM", 
    "QCD_Pt-300to470_TuneZ2star_8TeV_pythia6__Summer12_DR53X-PU_S10_START53_V7A-v2__AODSIM", 
    "QCD_Pt-470to600_TuneZ2star_8TeV_pythia6__Summer12_DR53X-PU_S10_START53_V7A-v2__AODSIM", 
    "QCD_Pt-600to800_TuneZ2star_8TeV_pythia6__Summer12_DR53X-PU_S10_START53_V7A-v2__AODSIM", 
    "QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6__Summer12_DR53X-PU_S10_START53_V7A-v2__AODSIM", 
    "TTJets_HadronicMGDecays_8TeV-madgraph__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM" } ; 

  std::cout << " Dataset " << " Events after full selection\n" ; 
  for (int ii = 0; ii < 27; ++ii) {
    TFile* fin = TFile::Open(dir+fnames[ii]+".root", "READ") ; 
    TH1D* hist = (TH1D*)fin->Get("BprimebH/h_cutflow") ;
    std::cout << fnames[ii] << " : " <<  hist->GetBinContent(7) << std::endl ;
  }

}
