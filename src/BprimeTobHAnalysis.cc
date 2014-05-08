// -*- C++ -*-
//
// Package:    BprimeTobHAnalysis
// Class:      BprimeTobHAnalysis
// 
/**\class BprimeTobHAnalysis BprimeTobHAnalysis.cc Bprime_kit/BprimeTobHAnalysis/src/BprimeTobHAnalysis.cc

Description: 
Analyzer class for Bprime -> b Higgs studies 
- National Taiwan University - 

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Eleni Petrakou,27 2-020,+41227674870,
//         Created:  Tue Jul 16 19:48:47 CEST 2013
// Second Author:    Devdatta Majumder 
// $Id$
//
//

// system include files
#include <memory>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <assert.h>
#include <vector>
#include <map>
#include <string>

// Root headers 
#include <TLorentzVector.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "BpbH/BprimeTobH/interface/format.h"
#include "BpbH/BprimeTobH/interface/TriggerBooking.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h" 

#include "BpbH/BprimeTobHAnalysis/interface/EventSelector.h"
#include "BpbH/BprimeTobHAnalysis/interface/JMEUncertUtil.h"
#include "BpbH/BprimeTobHAnalysis/interface/ApplyBTagSF.h"
#include "BpbH/BprimeTobHAnalysis/interface/ApplyHiggsTagSF.h"
#include "BpbH/BprimeTobHAnalysis/interface/HiggsBRscaleFactors.h" 

//
// class declaration
//

class BprimeTobHAnalysis : public edm::EDAnalyzer {
  public:
    explicit BprimeTobHAnalysis(const edm::ParameterSet&);
    ~BprimeTobHAnalysis();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    void CreateHistos(const TString&) ; 
    void AddHisto(const TString&, const TString&, const TString&, const int&, const double&, const double&) ; 
    void AddHisto(const TString&, const TString&, const TString&, const int&, const double&, const double&, const int&, const double&, const double&) ; 
    template <class Type>
      void FillHisto(const TString& name, const Type value, const double weight);
    template <class Type1, class Type2>
      void FillHisto(const TString& name, const Type1 value1, const Type2 value2, const double weight);

    edm::LumiReWeighting LumiWeights_; 

    // ----------member data ---------------------------

    //// Configurables 

    int                             maxEvents_; 
    const int                       reportEvery_; 
    const std::string               inputTTree_;
    const std::vector<std::string>  inputFiles_;
    const edm::ParameterSet         hltPaths_; 
    const bool                      doPUReweighting_ ;
    const std::string               file_PUDistMC_ ;
    const std::string               file_PUDistData_ ;
    const std::string               hist_PUDistMC_ ;
    const std::string               hist_PUDistData_ ;

    const double jetPtMin_ ; 
    const double jetPtMax_ ; 
    const double jetAbsEtaMax_ ;
    const double bjetPtMin_ ; 
		const double bjetCSVDiscMin_;
		const double bjetCSVDiscMax_;
    const double fatJetPtMin_ ;
    const double fatJetPtMax_ ; 
    const double fatJetMassMin_ ;
    const double fatJetMassMax_ ; 
    const double fatJetPrunedMassMin_ ;
    const double fatJetPrunedMassMax_ ; 
    const double dRSubjetsMin_ ; 
    const double dRSubjetsMax_ ; 
    const double subj1CSVDiscMin_ ; 
    const double subj1CSVDiscMax_ ; 
    const double subj2CSVDiscMin_ ; 
    const double subj2CSVDiscMax_ ; 
    const double HTAK5Min_ ; 
    const double HTAK5Max_ ; 

    const edm::ParameterSet jetSelParams_ ; 
    const edm::ParameterSet fatJetSelParams_ ; 
    const edm::ParameterSet higgsJetSelParams_ ; 
    const edm::ParameterSet HTSelParams_ ; 
    const edm::ParameterSet evtSelParams_ ; 
    const edm::ParameterSet jmeParams_; 
    const bool   applyJEC_ ; 
    const bool   applyBTagSF_ ; 
    const double jesShift_;
    const double jerShift_; 
    const double SFbShift_;
    const double SFlShift_;
    const bool   doTrigEff_;
    const bool   doABCDPlots_; 
    const bool   fillBDTTrees_; 

    TChain*            chain_;

    EvtInfoBranches    EvtInfo;
    VertexInfoBranches VtxInfo;
    GenInfoBranches    GenInfo;
    JetInfoBranches    GenJetInfo;
    JetInfoBranches    JetInfo;
    JetInfoBranches    FatJetInfo;
    JetInfoBranches    SubJetInfo;
    LepInfoBranches    LepInfo;

    bool isData_ ; 
    double evtwt_ ; 
    double puweight_ ; 

    edm::Service<TFileService> fs; 
    TH1D *h_bprimeDecayModes, *h_cutflow ; 
    std::map<TString, TH1D*> hmap_1d ;  
    std::map<TString, TH2D*> hmap_2d ;  

    EventSelector *eSelector_ ; 
    std::vector<std::string> cutLevels_; 

};

//
// constructors and destructor
//
BprimeTobHAnalysis::BprimeTobHAnalysis(const edm::ParameterSet& iConfig) : 
  maxEvents_(iConfig.getParameter<int>("MaxEvents")), 
  reportEvery_(iConfig.getParameter<int>("ReportEvery")),
  inputTTree_(iConfig.getParameter<std::string>("InputTTree")),
  inputFiles_(iConfig.getParameter<std::vector<std::string> >("InputFiles")),
  hltPaths_(iConfig.getParameter<edm::ParameterSet>("HLTPaths")),
  doPUReweighting_(iConfig.getParameter<bool>("DoPUReweighting")),
  file_PUDistMC_(iConfig.getParameter<std::string>("File_PUDistMC")),
  file_PUDistData_(iConfig.getParameter<std::string>("File_PUDistData")),
  hist_PUDistMC_(iConfig.getParameter<std::string>("Hist_PUDistMC")),
  hist_PUDistData_(iConfig.getParameter<std::string>("Hist_PUDistData")),
  jetPtMin_(iConfig.getParameter<double>("JetPtMin")),
  jetPtMax_(iConfig.getParameter<double>("JetPtMax")),
  jetAbsEtaMax_(iConfig.getParameter<double>("JetAbsEtaMax")),
  bjetPtMin_(iConfig.getParameter<double>("BJetPtMin")),
	bjetCSVDiscMin_(iConfig.getParameter<double>("BJetCSVDiscMin")),
	bjetCSVDiscMax_(iConfig.getParameter<double>("BJetCSVDiscMax")),
  fatJetPtMin_(iConfig.getParameter<double>("FatJetPtMin")),
  fatJetPtMax_(iConfig.getParameter<double>("FatJetPtMax")), 
  fatJetMassMin_(iConfig.getParameter<double>("FatJetMassMin")),
  fatJetMassMax_(iConfig.getParameter<double>("FatJetMassMax")), 
  fatJetPrunedMassMin_(iConfig.getParameter<double>("FatJetPrunedMassMin")),
  fatJetPrunedMassMax_(iConfig.getParameter<double>("FatJetPrunedMassMax")),
  dRSubjetsMin_(iConfig.getParameter<double>("DRSubjetsMin")),
  dRSubjetsMax_(iConfig.getParameter<double>("DRSubjetsMax")),
  subj1CSVDiscMin_(iConfig.getParameter<double>("Subjet1CSVDiscMin")),
  subj1CSVDiscMax_(iConfig.getParameter<double>("Subjet1CSVDiscMax")),
  subj2CSVDiscMin_(iConfig.getParameter<double>("Subjet2CSVDiscMin")),
  subj2CSVDiscMax_(iConfig.getParameter<double>("Subjet2CSVDiscMax")),
  HTAK5Min_(iConfig.getParameter<double>("HTAK5Min")), 
  HTAK5Max_(iConfig.getParameter<double>("HTAK5Max")),
  jetSelParams_(iConfig.getParameter<edm::ParameterSet>("JetSelParams")), 
  fatJetSelParams_(iConfig.getParameter<edm::ParameterSet>("FatJetSelParams")), 
  higgsJetSelParams_(iConfig.getParameter<edm::ParameterSet>("HiggsJetSelParams")), 
  HTSelParams_(iConfig.getParameter<edm::ParameterSet>("HTSelParams")), 
  evtSelParams_(iConfig.getParameter<edm::ParameterSet>("EvtSelParams")), 
  jmeParams_(iConfig.getParameter<edm::ParameterSet>("JMEParams")),
  applyJEC_(iConfig.getParameter<bool>("ApplyJEC")),
  applyBTagSF_(iConfig.getParameter<bool>("ApplyBTagSF")),
  jesShift_(iConfig.getParameter<double>("JESShift")),
  jerShift_(iConfig.getParameter<double>("JERShift")),
  SFbShift_(iConfig.getParameter<double>("SFbShift")),
  SFlShift_(iConfig.getParameter<double>("SFlShift")),
  doTrigEff_(iConfig.getParameter<double>("DoTrigEff")),
  doABCDPlots_(iConfig.getParameter<double>("DoABCDPlots")),
  fillBDTTrees_(iConfig.getParameter<double>("FillBDTTrees")), 
  isData_(0),
  evtwt_(1), 
  puweight_(1)  
{ 

  if (doPUReweighting_) LumiWeights_ = edm::LumiReWeighting(file_PUDistMC_, file_PUDistData_, hist_PUDistMC_, hist_PUDistData_) ;

}


BprimeTobHAnalysis::~BprimeTobHAnalysis() { 
  delete chain_;
}

// ------------ method called once each job just before starting event loop  ------------
void BprimeTobHAnalysis::beginJob() { 

  h_bprimeDecayModes = fs->make<TH1D>("h_bprimeDecayModes", "b' decay modes", 8, 0, 8) ; 
  h_bprimeDecayModes -> GetXaxis() -> SetBinLabel(1,"N(b')") ; 
  h_bprimeDecayModes -> GetXaxis() -> SetBinLabel(2,"b'b' -> bHbH") ; 
  h_bprimeDecayModes -> GetXaxis() -> SetBinLabel(3,"b'b' -> bHbZ") ; 
  h_bprimeDecayModes -> GetXaxis() -> SetBinLabel(4,"b'b' -> bHtW") ; 
  h_bprimeDecayModes -> GetXaxis() -> SetBinLabel(5,"b'b' -> bZbZ") ; 
  h_bprimeDecayModes -> GetXaxis() -> SetBinLabel(6,"b'b' -> bZtW") ; 
  h_bprimeDecayModes -> GetXaxis() -> SetBinLabel(7,"b'b' -> tWtW") ; 
  h_bprimeDecayModes -> GetXaxis() -> SetBinLabel(8,"Others") ; 

  int ncutflowbins(12) ; 
  h_cutflow = fs->make<TH1D>("h_cutflow" ,"Cut flow" ,ncutflowbins, 0 ,ncutflowbins) ;  
  h_cutflow -> Sumw2() ; 

  h_cutflow -> GetXaxis() -> SetBinLabel(1 ,"AllEvents") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(2 ,"SkimmedEvents") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(3 ,"TriggerSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(4 ,"FatJetSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(5 ,"HiggsJetSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(6 ,"JetSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(7 ,"BJetsSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(8 ,"HTSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(9 ,"HTSel_1H_1b") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(10,"HTSel_1H_2b") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(11,"HTSel_2H_1b") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(12,"HTSel_2H_2b") ; 

  for (int ii = 1; ii <= h_cutflow->GetNbinsX(); ++ii) 
    CreateHistos(h_cutflow->GetXaxis()->GetBinLabel(ii)) ; 

  chain_ = new TChain(inputTTree_.c_str());

  for(unsigned i=0; i<inputFiles_.size(); ++i) {
    chain_->Add(inputFiles_.at(i).c_str());

    TFile *f = TFile::Open(inputFiles_.at(i).c_str(),"READ");
    if ( TString(inputTTree_).Contains("skim") ) {
      TH1F* h_events = (TH1F*)f->Get("skim/h_cutflow") ;  
      if ( maxEvents_ < 0 ) {
        h_cutflow->Fill("AllEvents", h_events->GetBinContent(1)) ; 
        h_cutflow->Fill("SkimmedEvents", h_events->GetBinContent(4)) ; 
        h_cutflow->SetBinError(1, h_events->GetBinError(1)) ; 
        h_cutflow->SetBinError(2, h_events->GetBinError(4)) ; 
      }
    }
    else if ( TString(inputTTree_).Contains("bVeto") ) {
      TH1F* h_events = (TH1F*)f->Get("bVeto/h_cutflow") ;  
      if ( maxEvents_ < 0 ) {
        h_cutflow->Fill("AllEvents", h_events->GetBinContent(1)) ; 
        h_cutflow->Fill("SkimmedEvents", h_events->GetBinContent(6)) ; 
        h_cutflow->SetBinError(1, h_events->GetBinError(1)) ; 
        h_cutflow->SetBinError(2, h_events->GetBinError(6)) ; 
      }
    }
    f->Close();
  }

  EvtInfo.Register(chain_);
  VtxInfo.Register(chain_);
  GenInfo.Register(chain_);
  GenJetInfo.Register(chain_,"GenJetInfo");
  JetInfo.Register(chain_,"JetInfo");
  FatJetInfo.Register(chain_,"FatJetInfo");
  SubJetInfo.Register(chain_,"SubJetInfo");
  LepInfo.Register(chain_);
  eSelector_ = new EventSelector(evtSelParams_,EvtInfo,VtxInfo,JetInfo,FatJetInfo,SubJetInfo);  
  cutLevels_ = eSelector_->getCutLevels();

  if(maxEvents_<0 || maxEvents_>chain_->GetEntries()) maxEvents_ = chain_->GetEntries();

  return ;  

}

void BprimeTobHAnalysis::CreateHistos(const TString& cutname) {

  AddHisto(cutname ,"_nPVtx_NoPUWt"                ,"N(PV), No PU weight"                          ,50     ,-0.5     ,49.5    ) ; 
  AddHisto(cutname ,"_nPVtx_PUWt"                  ,"N(PV)"                                        ,50     ,-0.5     ,49.5    ) ; 

  AddHisto(cutname ,"_nJets"                       ,"N(AK5 jets)"                                  ,20     ,-0.5     ,19.5    ) ; 

  AddHisto(cutname ,"_nAllAK5"                     ,"N(All AK5 jets)"                              ,20     ,-0.5     ,19.5    ) ; 
  AddHisto(cutname ,"_AK5Jet_1_Pt"                 ,";p_{T} (leading AK5 jet)[GeV]; Events"        ,100    ,0.       ,2000.   ) ;
  AddHisto(cutname ,"_AK5Jet_2_Pt"                 ,";p_{T} (2nd AK5 jet)[GeV]; Events"            ,100    ,0.       ,2000.   ) ;
  AddHisto(cutname ,"_AK5Jet_3_Pt"                 ,";p_{T} (3rd AK5 jet)[GeV]; Events"            ,100    ,0.       ,2000.   ) ;
  AddHisto(cutname ,"_AK5Jet_4_Pt"                 ,";p_{T} (4th AK5 jet)[GeV]; Events"            ,100    ,0.       ,2000.   ) ;
  AddHisto(cutname ,"_AK5Jet_1_Eta"                ,";#eta (leading AK5 jet); Events"              ,50     ,-4.      ,4.      ) ;
  AddHisto(cutname ,"_AK5Jet_2_Eta"                ,";#eta (2nd AK5 jet); Events"                  ,50     ,-4.      ,4.      ) ;
  AddHisto(cutname ,"_AK5Jet_3_Eta"                ,";#eta (3rd AK5 jet); Events"                  ,50     ,-4.      ,4.      ) ;
  AddHisto(cutname ,"_AK5Jet_4_Eta"                ,";#eta (4th AK5 jet); Events"                  ,50     ,-4.      ,4.      ) ;

  AddHisto(cutname ,"_nFatJets"                    ,"N(fat jets)"                                  ,20     ,-0.5     ,19.5    ) ; 
  AddHisto(cutname ,"_FatJetsNoMassCut_Pt"         ,";p_{T}(fat jets)| No mass cut;"               ,1000   ,0.       ,1000.   ) ; 
  AddHisto(cutname ,"_FatJets_Pt"                  ,"p_{T}(fat jets)"                              ,1000   ,0.       ,1000.   ) ; 
  AddHisto(cutname ,"_FatJets_Eta"                 ,"#eta(fat jets)"                               ,50     ,-4.      ,4.      ) ; 
  AddHisto(cutname ,"_FatJets_tau2ByTau1"          ,"Fat jet #tau2/#tau1"                          ,20     ,0.       ,1.      ) ; 
  AddHisto(cutname ,"_FatJets_tau3ByTau2"          ,"Fat jet #tau2/#tau1"                          ,20     ,0.       ,1.      ) ; 
  AddHisto(cutname ,"_FatJets_tau3ByTau1"          ,"Fat jet #tau2/#tau1"                          ,20     ,0.       ,1.      ) ; 
  AddHisto(cutname ,"_FatJets_CombinedSVBJetTags"  ,"Fat jet CSV discriminator"                    ,20     ,0.       ,1.      ) ; 
  AddHisto(cutname ,"_FatJets_Mass"                ,"Fat jet mass [GeV]"                           ,200    ,0.       ,2000.   ) ; 
  AddHisto(cutname ,"_FatJets_MassPruned"          ,"Fat jet pruned mass [GeV]"                    ,200    ,0.       ,2000.   ) ; 
  AddHisto(cutname ,"_FatJets_DRSubjets"           ,";#DeltaR_{#eta,#phi}(subj1, subj2);"      ,100    ,0.       ,2.      ) ; 
  AddHisto(cutname ,"_FatJets_DYPhiSubjets"        ,";#DeltaR_{#y,#phi}(subj1, subj2);"        ,100    ,0.       ,2.      ) ; 

  AddHisto(cutname ,"_SubJet1_Pt"                  ,"SubJet1 p_{T} [GeV]"                          ,200    ,0.       ,2000.   ) ;
  AddHisto(cutname ,"_SubJet1_Eta"                 ,"SubJet1 #eta"                                 ,50     ,-4.      ,4.      ) ;
  AddHisto(cutname ,"_SubJet1_Mass"                ,"SubJet1 mass [GeV]"                           ,200    ,0.       ,2000.   ) ; 
  AddHisto(cutname ,"_SubJet1_CombinedSVBJetTags"  ,"SubJet1 CSV discriminator"                    ,20     ,0.       ,1.      ) ; 

  AddHisto(cutname ,"_SubJet2_Pt"                  ,"SubJet2 p_{T} [GeV]"                          ,200    ,0.       ,2000.   ) ;
  AddHisto(cutname ,"_SubJet2_Eta"                 ,"SubJet2 #eta"                                 ,50     ,-4.      ,4.      ) ;
  AddHisto(cutname ,"_SubJet2_Mass"                ,"SubJet2 mass [GeV]"                           ,200    ,0.       ,2000.   ) ; 
  AddHisto(cutname ,"_SubJet2_CombinedSVBJetTags"  ,"SubJet2 CSV discriminator"                    ,20     ,0.       ,1.      ) ; 

  AddHisto(cutname ,"_nHJets"                      ,"N(Higgs jets)"                                ,20     ,-0.5     ,19.5    ) ; 
  AddHisto(cutname ,"_HiggsJet_Pt"                 ,"p_{T} (Higgs jet)[GeV]"                       ,200    ,0.       ,2000.   ) ;
  AddHisto(cutname ,"_HiggsJet_Eta"                ,"#eta (Higgs jet)"                             ,50     ,-4.      ,4.      ) ;
  AddHisto(cutname ,"_HiggsJet_Mass"               ,"Mass (Higgs jet) [GeV]"                       ,40     ,0.       ,400.    ) ; 

  AddHisto(cutname ,"_nHJetsNoMDyCut"              ,"N(Higgs jets)"                                ,20     ,-0.5     ,19.5    ) ; 
  AddHisto(cutname ,"_HiggsJetNoMDyCut_Pt"         ,"p_{T} (Higgs jet)[GeV]"                       ,200    ,0.       ,2000.   ) ;
  AddHisto(cutname ,"_HiggsJetNoMDyCut_Eta"        ,"#eta (Higgs jet)"                             ,50     ,-4.      ,4.      ) ;
  AddHisto(cutname ,"_HiggsJetNoMDyCut_Mass"       ,"Mass (Higgs jet) [GeV]"                       ,40     ,0.       ,400.    ) ; 
  AddHisto(cutname ,"_HiggsJetNoMDyCut_DRsubjets"  ,"#DeltaR_{y,#phi}(subjets)"                    ,100    ,0.       ,2.      ) ; 

  AddHisto(cutname ,"_nBJets"                      ,"N(b jets)"                                    ,20     ,-0.5     ,19.5    ) ; 
  AddHisto(cutname ,"_BJet_Pt"                     ,"p_{T} (b jet)[GeV]"                           ,200    ,0.       ,2000.   ) ;
  AddHisto(cutname ,"_BJet_Eta"                    ,"#eta (b jet)"                                 ,50     ,-4.      ,4.      ) ;

  AddHisto(cutname ,"_HT"                          ,"H_{T}[GeV]"                                   ,400   ,0.       ,4000.   ) ;
  AddHisto(cutname ,"_HTAK5"                       ,"H_{T} (AK5 jets, no Higgs jet overlap) [GeV]" ,400   ,0.       ,4000.   ) ;
  AddHisto(cutname ,"_HTAK5_leading4"              ,"H_{T} (leading 4 AK5 jets) [GeV]"             ,400   ,0.       ,4000.   ) ;
  AddHisto(cutname ,"_HTCA8_leading2_AK5_leading2" ,"H_{T} (leading 2 CA8 and AK5 jets) [GeV]"     ,400   ,0.       ,4000.   ) ;
  AddHisto(cutname ,"_HTAllAK5"                    ,"H_{T} (AK5 jets) [GeV]"                       ,400   ,0.       ,4000.   ) ;
  AddHisto(cutname ,"_Nbprimes"                    ,";N (b' candidates); Events"                   ,5     ,-0.5     ,4.5     ) ;
  AddHisto(cutname ,"_bprimePt"                    ,"b' p_{T} [GeV]"                               ,200   ,0.       ,2000.   ) ;
  AddHisto(cutname ,"_bprimeMass"                  ,"b' mass [GeV]"                                ,200   ,0.       ,2000.   ) ;

  if ( doABCDPlots_ ) {
    AddHisto(cutname ,"_nAntiHjets"                ,";N(Anti Higgs-tagged jets);;"                 ,20     ,-0.5     ,19.5    ) ; 
    AddHisto(cutname ,"_AntiHjet_Pt"               ,";p_{T} (Anti Higgs-tagged jet)[GeV];;"        ,200    ,0.       ,2000.   ) ;
    AddHisto(cutname ,"_AntiHjet_Eta"              ,";#eta (Anti Higgs-tagged jet)''"              ,50     ,-4.      ,4.      ) ;
    AddHisto(cutname ,"_AntiHjet_Mass"             ,";Mass (Anti Higgs-tagged jet) [GeV];;"        ,40     ,0.       ,400.    ) ; 

    AddHisto(cutname ,"_Htag_HT"                   ,";HT [GeV];;;"      ,200 ,0. ,2000. ,2 ,-.5, 1.5) ; 
    AddHisto(cutname ,"_Htag_HTAK5"                ,";HT(AK5) [GeV];;;" ,200 ,0. ,2000. ,2 ,-.5, 1.5) ; 
    AddHisto(cutname ,"_2b_Htag_HT"                ,";HT [GeV];;;"      ,200 ,0. ,2000. ,2 ,-.5, 1.5) ; 
    AddHisto(cutname ,"_2b_Htag_HTAK5"             ,";HT(AK5) [GeV];;;" ,200 ,0. ,2000. ,2 ,-.5, 1.5) ; 
  }

  return ; 

}

void BprimeTobHAnalysis::AddHisto(const TString& cutname, const TString& histname, const TString& histtitle, const int& nbins, const double& min, const double& max) { 
  TH1D* h1d ; 
  h1d = fs->make<TH1D>(cutname+histname, cutname+histtitle, nbins, min, max);  
  h1d -> Sumw2() ; 
  hmap_1d[cutname+histname] = h1d ; 
  return ; 
}

void BprimeTobHAnalysis::AddHisto(const TString& cutname, const TString& histname, const TString& histtitle, const int& nbinsx, const double& minx, const double& maxx, const int& nbinsy, const double& miny, const double& maxy) { 
  TH2D* h2d ; 
  h2d = fs->make<TH2D>(cutname+histname, cutname+histtitle, nbinsx, minx, maxx, nbinsy, miny, maxy);  
  h2d -> Sumw2() ; 
  hmap_2d[cutname+histname] = h2d ; 
  return ; 
}

template <class Type>
void BprimeTobHAnalysis::FillHisto(const TString& name, const Type value, const double weight){
  hmap_1d[name]->Fill(double(value),weight);
  return ; 
}

template <class Type1, class Type2>
void BprimeTobHAnalysis::FillHisto(const TString& name, const Type1 value1, const Type2 value2, const double weight){
  hmap_2d[name]->Fill(double(value1),double(value2),weight);
  return ; 
}

// ------------ method called for each event  ------------
void BprimeTobHAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) { 
  using namespace edm;
  using namespace std;

  if(chain_ == 0) return;

  JetID jetIDTight(JetID::FIRSTDATA,JetID::TIGHT, JetInfo) ; 
  JetID fatjetIDLoose(JetID::FIRSTDATA,JetID::LOOSE, FatJetInfo) ; 
  pat::strbitset retak5 = jetIDTight.getBitTemplate() ;
  pat::strbitset retca8 = fatjetIDLoose.getBitTemplate() ;

  FatJetSelector jetSelCA8(fatJetSelParams_) ; 
  FatJetSelector jetSelHiggs(higgsJetSelParams_) ; 
  JetSelector jetSelAK5(jetSelParams_) ; 
  pat::strbitset retjetidak5 = jetSelAK5.getBitTemplate() ; 
  pat::strbitset retjetidca8 = jetSelCA8.getBitTemplate() ; 

  edm::LogInfo("StartingAnalysisLoop") << "Starting analysis loop\n";

  for(int entry=0; entry<maxEvents_; entry++) {

    if((entry%reportEvery_) == 0) edm::LogInfo("Event") << entry << " of " << maxEvents_ ; 

    //// Event variables 
    bool passHLT(false), passBSel(false) ; 
    int nGoodVtxs(0) ;
    JetCollection fatjets_corr, HjetsNoMassDyCuts_corr, Hjets_corr, AntiHjets_corr, ak5jets_corr, allAK5jets_corr, ak5jets_leading4, bjets, jets_CA8_leading2_AK5_leading2 ; 
    HT HTAK5, HTAllAK5, HTAK5_leading4, HTCA8_leading2_AK5_leading2, MyHT ; 
    std::vector<TLorentzVector>p4bprimes ; 

    chain_->GetEntry(entry);

    isData_   = EvtInfo.McFlag ? 0 : 1; 
    if ( !isData_ ) evtwt_    = EvtInfo.Weight ; 
    if ( doPUReweighting_ && !isData_ ) puweight_ = LumiWeights_.weight(EvtInfo.TrueIT[0]) ; 

    if ( !isData_ ) { //// Gen info 
      std::vector<int> gen_bp_indices, gen_higgs_indices ; 
      for (int igen = 0; igen < GenInfo.Size; ++igen) { 
        if ( GenInfo.Status[igen] == 3 && abs(GenInfo.PdgID[igen]) == 7  ) { //// bprime found 
          h_bprimeDecayModes -> AddBinContent(1) ; 
          gen_bp_indices.push_back(igen) ; 
        }
        if ( GenInfo.Status[igen] == 3 && abs(GenInfo.PdgID[igen]) == 25 && GenInfo.nDa[igen] >= 2  ) gen_higgs_indices.push_back(igen) ; 
      }
      int nbprimeBH(0), nbprimeBZ(0), nbprimeTW(0) ; 
      for (std::vector<int>::const_iterator ibp = gen_bp_indices.begin(); ibp != gen_bp_indices.end(); ++ibp) {
        int bprimeDau0(abs(GenInfo.Da0PdgID[*ibp])), bprimeDau1(abs(GenInfo.Da1PdgID[*ibp])) ;
        if ( (bprimeDau0 == 5 && bprimeDau1 == 25) || (bprimeDau0 == 25 && bprimeDau1 == 5) ) ++nbprimeBH ; 
        else if ( (bprimeDau0 == 5 && bprimeDau1 == 23) || (bprimeDau0 == 23 && bprimeDau1 == 5) ) ++nbprimeBZ ; 
        else if ( (bprimeDau0 == 6 && bprimeDau1 == 24) || (bprimeDau0 == 24 && bprimeDau1 == 6) ) ++nbprimeTW ; 
      }
      if      ( nbprimeBH == 2 )                   h_bprimeDecayModes -> AddBinContent(2) ; 
      else if ( nbprimeBH == 1 && nbprimeBZ == 1 ) h_bprimeDecayModes -> AddBinContent(3) ; 
      else if ( nbprimeBH == 1 && nbprimeTW == 1 ) h_bprimeDecayModes -> AddBinContent(4) ; 
      else if ( nbprimeBZ == 2 )                   h_bprimeDecayModes -> AddBinContent(5) ; 
      else if ( nbprimeBZ == 1 && nbprimeTW == 1 ) h_bprimeDecayModes -> AddBinContent(6) ; 
      else if ( nbprimeTW == 2 )                   h_bprimeDecayModes -> AddBinContent(7) ; 
      else                                         h_bprimeDecayModes -> AddBinContent(8) ; 
      //// Higgs BR reweighting
      for (std::vector<int>::const_iterator ihig = gen_higgs_indices.begin(); ihig != gen_higgs_indices.end(); ++ihig) {
        int higgsDau0(abs(GenInfo.Da0PdgID[*ihig])), higgsDau1(abs(GenInfo.Da1PdgID[*ihig])) ;
        int higgsDau = higgsDau0+higgsDau1 ; 
        if(higgsDau==10) evtwt_ *= HiggsBRscaleFactors::higgsBBSf;
        else if(higgsDau==30) evtwt_ *= HiggsBRscaleFactors::higgsTauTauSf;
        else if(higgsDau==26) evtwt_ *= HiggsBRscaleFactors::higgsMuMuSf;
        else if(higgsDau==8)  evtwt_ *= HiggsBRscaleFactors::higgsCCSf;
        else if(higgsDau==6)  evtwt_ *= HiggsBRscaleFactors::higgsSSSf;
        else if(higgsDau==16) evtwt_ *= HiggsBRscaleFactors::higgsTTSf;
        else if(higgsDau==42) evtwt_ *= HiggsBRscaleFactors::higgsGGSf;
        else if(higgsDau==44) evtwt_ *= HiggsBRscaleFactors::higgsGammaGammaSf;
        else if(higgsDau==45) evtwt_ *= HiggsBRscaleFactors::higgsZGammaSf;
        else if(higgsDau==48) evtwt_ *= HiggsBRscaleFactors::higgsWWSf;
        else if(higgsDau==46) evtwt_ *= HiggsBRscaleFactors::higgsZZSf; 
        else evtwt_ *= 1 ; 
      } //// Higgs BR reweighting  
    } //// Gen info

    nGoodVtxs = 0 ;
    VertexSelector vtxSel(VtxInfo) ; 
    nGoodVtxs = vtxSel.NGoodVtxs(); 
    if (nGoodVtxs < 1)  { edm::LogInfo("NoGoodPrimaryVertex") << " No good primary vertex " ; continue ; }
    if ( TString(inputTTree_).Contains("ntuple") ) h_cutflow -> Fill("AllEvents", 1) ; 
    FillHisto(TString("AllEvents")+TString("_nPVtx_NoPUWt"), nGoodVtxs, evtwt_) ; 
    FillHisto(TString("AllEvents")+TString("_nPVtx_PUWt"), nGoodVtxs, evtwt_*puweight_) ; 
    FillHisto(TString("AllEvents")+TString("_nJets"), JetInfo.Size, evtwt_*puweight_) ; 

    TriggerSelector trigSel(hltPaths_) ; 
    passHLT = trigSel.getTrigDecision(EvtInfo) ; 
    if ( !doTrigEff_ && !passHLT ) continue ; 
    if ( !doTrigEff_ ) {
      h_cutflow -> Fill("TriggerSel", 1) ; 
      FillHisto(TString("TriggerSel")+TString("_nPVtx_NoPUWt"), nGoodVtxs, evtwt_) ; 
      FillHisto(TString("TriggerSel")+TString("_nPVtx_PUWt"), nGoodVtxs, evtwt_*puweight_) ; 
    }

    evtwt_ *= puweight_ ; 

    int nfatjets(0) ; 
    JetCollection fatjets ; 
    for (int ifatjet = 0; ifatjet < FatJetInfo.Size; ++ifatjet) { 
      if (jetSelCA8(FatJetInfo, ifatjet, SubJetInfo, retjetidca8) == 0) continue ; //// pt > 150 GeV, |eta| < 2.4 tau2/tau1 < 0.5 
      if (FatJetInfo.Pt[ifatjet] > 300) {
        FillHisto(TString("TriggerSel")+TString("_FatJetsNoMassCut_Pt"), FatJetInfo.Pt[ifatjet], evtwt_) ; 
        ++nfatjets ;
      }
      Jet thisjet(FatJetInfo, ifatjet) ; 
      fatjets.push_back(thisjet) ; 
    } //// Loop over fat jets 
    FillHisto(TString("AllEvents")+TString("_nFatJets"), nfatjets, evtwt_) ; 
    nfatjets = 0 ; 

    //// Apply JEC and b-tagging SFs for CA8 jets in MC  
    if ( !isData_ ) {
      JMEUncertUtil* jmeUtil_jer = new JMEUncertUtil(jmeParams_, fatjets, "JERCA8MC", jerShift_) ; 
      JetCollection ca8jets_jer = jmeUtil_jer->GetModifiedJetColl() ; 
      delete jmeUtil_jer ; 
      fatjets.clear() ; 

      if ( abs(jesShift_) > 1E-6 ) {
        JMEUncertUtil* jmeUtil_jes = new JMEUncertUtil(jmeParams_, ca8jets_jer, "JESCA8MC", jesShift_) ; 
        JetCollection ca8jets_jes = jmeUtil_jes->GetModifiedJetColl() ; 
        delete jmeUtil_jes ; 

        for (JetCollection::const_iterator ijet = ca8jets_jes.begin(); ijet != ca8jets_jes.end(); ++ijet) {
          Jet thisjet(*ijet) ; 
          fatjets.push_back(thisjet) ; 
        }
      }
      else {
        for (JetCollection::const_iterator ijet = ca8jets_jer.begin(); ijet != ca8jets_jer.end(); ++ijet) {
          Jet thisjet(*ijet) ; 
          fatjets.push_back(thisjet) ; 
        }
      }
    } //// Apply JEC and b-tagging SFs for CA8 jets in MC 

    for (JetCollection::const_iterator ifat = fatjets.begin(); ifat != fatjets.end(); ++ifat) {
      TLorentzVector fatjet_p4;
      fatjet_p4.SetPtEtaPhiM(ifat->Pt(), ifat->Eta(), ifat->Phi (), ifat->Mass ()); 
      int iSubJet1 = ifat->Jet_SubJet1Idx();
      int iSubJet2 = ifat->Jet_SubJet2Idx();
      TLorentzVector subj1_p4, subj2_p4;
      subj1_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet1], SubJetInfo.Eta[iSubJet1], SubJetInfo.Phi[iSubJet1], SubJetInfo.Mass[iSubJet1]);
      subj2_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet2], SubJetInfo.Eta[iSubJet2], SubJetInfo.Phi[iSubJet2], SubJetInfo.Mass[iSubJet2]); 
      double subjet_dR = subj1_p4.DeltaR(subj2_p4);
      double subjet_dy = subj1_p4.Rapidity() - subj2_p4.Rapidity() ;
      double subjet_dphi = subj1_p4.DeltaPhi(subj2_p4); ;
      double subjet_dyphi = sqrt( subjet_dy*subjet_dy + subjet_dphi*subjet_dphi ) ;

      FillHisto(TString("TriggerSel")+TString("_FatJets_Pt")                 ,fatjet_p4.Pt() ,evtwt_)  ;  
      FillHisto(TString("TriggerSel")+TString("_FatJets_Eta")                ,fatjet_p4.Eta() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_FatJets_CombinedSVBJetTags") ,ifat->CombinedSVBJetTags() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_FatJets_tau2ByTau1")         ,ifat->tau2()/ifat->tau1() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_FatJets_tau3ByTau2")         ,ifat->tau3()/ifat->tau2() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_FatJets_tau3ByTau1")         ,ifat->tau3()/ifat->tau1() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_FatJets_Mass")               ,fatjet_p4.Mag() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_FatJets_MassPruned")         ,ifat->MassPruned() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_FatJets_DRSubjets")          ,subjet_dR ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_FatJets_DYPhiSubjets")       ,subjet_dyphi ,evtwt_)  ; 

      FillHisto(TString("TriggerSel")+TString("_SubJet1_Pt") ,subj1_p4.Pt() ,evtwt_)  ;  
      FillHisto(TString("TriggerSel")+TString("_SubJet1_Eta") ,subj1_p4.Eta() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_SubJet1_Mass") ,subj1_p4.Mag() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_SubJet1_CombinedSVBJetTags") ,SubJetInfo.CombinedSVBJetTags[iSubJet1] ,evtwt_)  ; 

      FillHisto(TString("TriggerSel")+TString("_SubJet2_Pt") ,subj2_p4.Pt() ,evtwt_)  ;  
      FillHisto(TString("TriggerSel")+TString("_SubJet2_Eta") ,subj2_p4.Eta() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_SubJet2_Mass") ,subj2_p4.Mag() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_SubJet2_CombinedSVBJetTags") ,SubJetInfo.CombinedSVBJetTags[iSubJet2] ,evtwt_)  ; 

      if ( fatjet_p4.Pt() > fatJetPtMin_ && fatjet_p4.Pt() < fatJetPtMax_ ) {  

        Jet thisjet(*ifat) ; 

        //// Selecting fat jets with mass and dy cuts 
        if ( fatjet_p4.Mag() > fatJetMassMin_ && fatjet_p4.Mag() < fatJetMassMax_ 
            && ifat->MassPruned() > fatJetPrunedMassMin_ && ifat->MassPruned() < fatJetPrunedMassMax_ 
            && subjet_dyphi > dRSubjetsMin_ && subjet_dyphi < dRSubjetsMax_ ) {  
          fatjets_corr.push_back(thisjet) ; 
          ++nfatjets ; 
          if (nfatjets < 2) jets_CA8_leading2_AK5_leading2.push_back(thisjet) ;  
        } //// Selecting fat jets with mass and dy cuts 

        //// Higgs tagging
        if ( SubJetInfo.CombinedSVBJetTags[iSubJet1] > subj1CSVDiscMin_ 
            && SubJetInfo.CombinedSVBJetTags[iSubJet1] < subj1CSVDiscMax_ 
            && SubJetInfo.CombinedSVBJetTags[iSubJet2] > subj2CSVDiscMin_ 
            && SubJetInfo.CombinedSVBJetTags[iSubJet2] < subj2CSVDiscMax_) { 

          HjetsNoMassDyCuts_corr.push_back(thisjet) ; 

          if ( fatjet_p4.Mag() > fatJetMassMin_ && fatjet_p4.Mag() < fatJetMassMax_ 
              && ifat->MassPruned() > fatJetPrunedMassMin_ && ifat->MassPruned() < fatJetPrunedMassMax_ 
              && subjet_dyphi > dRSubjetsMin_ && subjet_dyphi < dRSubjetsMax_ ) {  
            FillHisto(TString("TriggerSel")+TString("_HiggsJet_Pt")   ,fatjet_p4.Pt() ,evtwt_)  ; 
            FillHisto(TString("TriggerSel")+TString("_HiggsJet_Eta")  ,fatjet_p4.Eta() ,evtwt_)  ;
            FillHisto(TString("TriggerSel")+TString("_HiggsJet_Mass") ,fatjet_p4.Mag() ,evtwt_)  ;

            Hjets_corr.push_back(thisjet) ; 

            if ( !isData_ ) { //// Apply Higgs-tagging scale factor 
              ApplyHiggsTagSF* higgsTagSF = new ApplyHiggsTagSF(double(SubJetInfo.Pt[iSubJet1]), double(SubJetInfo.Pt[iSubJet2]), 
                  double(SubJetInfo.Eta[iSubJet1]), double(SubJetInfo.Eta[iSubJet2]),
                  SubJetInfo.GenFlavor[iSubJet1], SubJetInfo.GenFlavor[iSubJet2], 
                  SubJetInfo.CombinedSVBJetTags[iSubJet1], SubJetInfo.CombinedSVBJetTags[iSubJet2]) ; 
              evtwt_ *= higgsTagSF->GetHiggsTagSF() ;
              delete higgsTagSF ; 
            }

          } //// Selecting Higgs jets with mass and dy cuts 
        } //// Higgs tagging 
        //// Anti-Higgs tagging : Attn: Hardcoding anti-Higgs jets to have subjet CSV disc. > CSVL 
        else if ( doABCDPlots_ 
            && (SubJetInfo.CombinedSVBJetTags[iSubJet1] < subj1CSVDiscMin_ && SubJetInfo.CombinedSVBJetTags[iSubJet1] > 0.244 ) 
            && (SubJetInfo.CombinedSVBJetTags[iSubJet2] < subj2CSVDiscMin_ && SubJetInfo.CombinedSVBJetTags[iSubJet2] > 0.244 ) 
            && fatjet_p4.Mag() > fatJetMassMin_ && fatjet_p4.Mag() < fatJetMassMax_ 
            && ifat->MassPruned() > fatJetPrunedMassMin_ && ifat->MassPruned() < fatJetPrunedMassMax_ 
            && subjet_dyphi > dRSubjetsMin_ && subjet_dyphi < dRSubjetsMax_ ) {  
          FillHisto(TString("TriggerSel")+TString("_AntiHjet_Pt")   ,fatjet_p4.Pt() ,evtwt_)  ; 
          FillHisto(TString("TriggerSel")+TString("_AntiHjet_Eta")  ,fatjet_p4.Eta() ,evtwt_)  ;
          FillHisto(TString("TriggerSel")+TString("_AntiHjet_Mass") ,fatjet_p4.Mag() ,evtwt_)  ;
          AntiHjets_corr.push_back(thisjet) ; 
        } //// Anit-Higgs tagging 

      } //// Applying fat jet pt cut 
    } //// Looping over all fat jets with pt > 150 GeV, |eta| < 2.4, loose jet ID, tau2/tau1 < 0.5 

    FillHisto(TString("TriggerSel")+TString("_nHJetsNoMDyCut")     ,HjetsNoMassDyCuts_corr.size() ,evtwt_)  ; 
    for (JetCollection::const_iterator ijet = HjetsNoMassDyCuts_corr.begin(); ijet != HjetsNoMassDyCuts_corr.end(); ++ijet ) {
      FillHisto(TString("TriggerSel")+TString("_HiggsJetNoMDyCut_Pt")   ,ijet->Pt() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_HiggsJetNoMDyCut_Eta")  ,ijet->Eta() ,evtwt_)  ;
      FillHisto(TString("TriggerSel")+TString("_HiggsJetNoMDyCut_Mass") ,ijet->Mass() ,evtwt_)  ;
      int iSubJet1(ijet->Jet_SubJet1Idx()), iSubJet2(ijet->Jet_SubJet2Idx()); 
      TLorentzVector subj1_p4, subj2_p4;
      subj1_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet1], SubJetInfo.Eta[iSubJet1], SubJetInfo.Phi[iSubJet1], SubJetInfo.Mass[iSubJet1]);
      subj2_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet2], SubJetInfo.Eta[iSubJet2], SubJetInfo.Phi[iSubJet2], SubJetInfo.Mass[iSubJet2]); 
      double subjet_dy = subj1_p4.Rapidity() - subj2_p4.Rapidity() ;
      double subjet_dphi = subj1_p4.DeltaPhi(subj2_p4); ;
      double subjet_dyphi = sqrt( subjet_dy*subjet_dy + subjet_dphi*subjet_dphi ) ;
      FillHisto(TString("TriggerSel")+TString("_HiggsJetNoMDyCut_DRsubjets"), subjet_dyphi, evtwt_) ; 
    }

    JetCollection ak5jetsNotHJets ;
    JetCollection allAK5jets ;
    for (int ijet = 0; ijet < JetInfo.Size; ++ijet) {
      retjetidak5.set(false) ;
      if (jetSelAK5(JetInfo, ijet,retjetidak5) == 0) continue ; 
      Jet thisjet(JetInfo, ijet) ;
      allAK5jets.push_back(thisjet) ; 
      bool isJetNotHiggs(false) ; 
      if (Hjets_corr.size() > 0) {
        for (JetCollection::const_iterator ihig = Hjets_corr.begin(); ihig != Hjets_corr.end(); ++ihig) {
          if (thisjet.DeltaR(*ihig) < 1.2) { //// Attn. hard-coded 
            isJetNotHiggs = false ; 
            break ; 
          }
          else isJetNotHiggs = true ; 
        } 
      }
      else isJetNotHiggs = true ; 
      if (isJetNotHiggs) ak5jetsNotHJets.push_back(thisjet) ; 
    }

    //// Apply JEC and b-tagging SFs to MC 
    //DMif ( !isData_ ) {
    //DM  if ( applyJEC_ ) { //// Apply JEC for MC   

    //DM    //// Only AK5 jets not overlapping with Higgs jets 
    //DM    JMEUncertUtil* jmeUtil_jer = new JMEUncertUtil(jmeParams_, ak5jetsNotHJets, "JERAK5MC", jerShift_) ; 
    //DM    JetCollection ak5jets_jer = jmeUtil_jer->GetModifiedJetColl() ; 
    //DM    delete jmeUtil_jer ; 
    //DM    ak5jetsNotHJets.clear() ; 

    //DM    if ( abs(jesShift_) > 1E-6 ) {
    //DM      JMEUncertUtil* jmeUtil_jes = new JMEUncertUtil(jmeParams_, ak5jets_jer, "JESAK5MC", jesShift_) ; 
    //DM      JetCollection ak5jets_jes = jmeUtil_jes->GetModifiedJetColl() ; 
    //DM      delete jmeUtil_jes ; 

    //DM      for (JetCollection::const_iterator ijet = ak5jets_jes.begin(); ijet != ak5jets_jes.end(); ++ijet) {
    //DM        Jet thisjet(*ijet) ; 
    //DM        ak5jetsNotHJets.push_back(thisjet) ; 
    //DM      }
    //DM    }
    //DM    else {
    //DM      for (JetCollection::const_iterator ijet = ak5jets_jer.begin(); ijet != ak5jets_jer.end(); ++ijet) {
    //DM        Jet thisjet(*ijet) ; 
    //DM        ak5jetsNotHJets.push_back(thisjet) ; 
    //DM      }
    //DM    }

    //DM    //// All AK5 jets 
    //DM    jmeUtil_jer = new JMEUncertUtil(jmeParams_, allAK5jets, "JERAK5MC", jerShift_) ; 
    //DM    JetCollection allAK5jets_jer = jmeUtil_jer->GetModifiedJetColl() ; 
    //DM    delete jmeUtil_jer ; 
    //DM    allAK5jets.clear() ; 

    //DM    if ( abs(jesShift_) > 1E-6 ) {
    //DM      JMEUncertUtil* jmeUtil_jes = new JMEUncertUtil(jmeParams_, allAK5jets_jer, "JESAK5MC", jesShift_) ; 
    //DM      JetCollection allAK5jets_jec = jmeUtil_jes->GetModifiedJetColl() ; 
    //DM      delete jmeUtil_jes ; 

    //DM      for (JetCollection::const_iterator ijet = allAK5jets_jec.begin(); ijet != allAK5jets_jec.end(); ++ijet) {
    //DM        Jet thisjet(*ijet) ; 
    //DM        allAK5jets.push_back(thisjet) ; 
    //DM      }
    //DM    }
    //DM    else {
    //DM      for (JetCollection::const_iterator ijet = allAK5jets_jer.begin(); ijet != allAK5jets_jer.end(); ++ijet) {
    //DM        Jet thisjet(*ijet) ; 
    //DM        allAK5jets.push_back(thisjet) ; 
    //DM      }
    //DM    }

    //DM  } //// Apply JEC for MC 

    //DM  if ( applyBTagSF_  ) { //// Apply b-tagging SF for MC  
    //DM    ApplyBTagSF * btagsf =  new ApplyBTagSF(ak5jetsNotHJets, 0.679, "CSVM", SFbShift_, SFlShift_) ;  
    //DM    ak5jetsNotHJets.clear() ; 
    //DM    ak5jetsNotHJets =  btagsf->getBtaggedJetsWithSF () ; 
    //DM    delete btagsf ; 
    //DM  } //// Apply b-tagging SF for MC 
    //DM} //// Apply JEC and b-tagging SFs to MC 

    for (JetCollection::const_iterator ijet = ak5jetsNotHJets.begin(); ijet != ak5jetsNotHJets.end(); ++ijet) {
      if (ijet->Pt() > jetPtMin_ && ijet->Pt() < jetPtMax_) { 
        ak5jets_corr.push_back(*ijet) ; 
        if ( ijet - ak5jetsNotHJets.begin() < 4 ) ak5jets_leading4.push_back(*ijet) ; 
        if ( ijet - ak5jetsNotHJets.begin() < 2 ) jets_CA8_leading2_AK5_leading2.push_back(*ijet) ;  
      }
      if (ijet->Pt() > bjetPtMin_ && ijet->CombinedSVBJetTags() > bjetCSVDiscMin_) bjets.push_back(*ijet) ; //// Select b-tagged AK5 jets 
    }

    for (JetCollection::const_iterator ijet = allAK5jets.begin(); ijet != allAK5jets.end(); ++ijet) {
      if (ijet->Pt() > jetPtMin_ && ijet->Pt() < jetPtMax_) allAK5jets_corr.push_back(*ijet) ; 
    }

    HTAK5.setJetCollection(ak5jets_corr) ; 
    HTAK5.buildHT() ; 

    HTAllAK5.setJetCollection(allAK5jets_corr) ; 
    HTAllAK5.buildHT() ; 

    HTAK5_leading4.setJetCollection(ak5jets_leading4) ; 
    HTAK5_leading4.buildHT() ; 

    HTCA8_leading2_AK5_leading2.setJetCollection(jets_CA8_leading2_AK5_leading2) ; 
    HTCA8_leading2_AK5_leading2.buildHT() ; 

    MyHT.setJetCollection(Hjets_corr) ; 
    MyHT.setJetCollection(bjets) ; 
    MyHT.buildHT() ; 

    if ( !doTrigEff_ ) {
      FillHisto(TString("TriggerSel")+TString("_nJets"), ak5jets_corr.size(), evtwt_) ; 
      FillHisto(TString("TriggerSel")+TString("_nFatJets"), fatjets_corr.size(), evtwt_) ; 
      FillHisto(TString("TriggerSel")+TString("_nAllAK5"), allAK5jets_corr.size() , evtwt_) ; 
      FillHisto(TString("TriggerSel")+TString("_nBJets"), bjets.size(), evtwt_) ; 
      FillHisto(TString("TriggerSel")+TString("_nHJets"), Hjets_corr.size(), evtwt_) ; 
      if ( doABCDPlots_ ) FillHisto(TString("TriggerSel")+TString("_nAntiHjets"), AntiHjets_corr.size(), evtwt_) ; 
      FillHisto(TString("TriggerSel")+TString("_HTCA8_leading2_AK5_leading2"), HTCA8_leading2_AK5_leading2.getHT(), evtwt_) ; 
      FillHisto(TString("TriggerSel")+TString("_HTAK5_leading4"), HTAK5_leading4.getHT(), evtwt_) ; 
      FillHisto(TString("TriggerSel")+TString("_HTAK5"), HTAK5.getHT(), evtwt_) ; 
      FillHisto(TString("TriggerSel")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt_) ; 
      FillHisto(TString("TriggerSel")+TString("_HT"), MyHT.getHT(), evtwt_) ; 
      for (JetCollection::const_iterator ib = bjets.begin(); ib != bjets.end(); ++ib) { 
        FillHisto(TString("TriggerSel")+TString("_BJet_Pt"), ib->Pt(), evtwt_) ; 
        FillHisto(TString("TriggerSel")+TString("_BJet_Eta"), ib->Eta(), evtwt_) ; 
      }

      if (allAK5jets_corr.size() > 0) FillHisto(TString("TriggerSel")+TString("_AK5Jet_1_Pt"), (allAK5jets_corr.at(0)).Pt(), evtwt_) ; 
      else FillHisto(TString("TriggerSel")+TString("_AK5Jet_1_Pt"), 0., evtwt_) ; 
      if (allAK5jets_corr.size() > 1) FillHisto(TString("TriggerSel")+TString("_AK5Jet_2_Pt"), (allAK5jets_corr.at(1)).Pt(), evtwt_) ; 
      else FillHisto(TString("TriggerSel")+TString("_AK5Jet_2_Pt"), 0., evtwt_) ; 
      if (allAK5jets_corr.size() > 2) FillHisto(TString("TriggerSel")+TString("_AK5Jet_3_Pt"), (allAK5jets_corr.at(2)).Pt(), evtwt_) ; 
      else FillHisto(TString("TriggerSel")+TString("_AK5Jet_3_Pt"), 0., evtwt_) ; 
      if (allAK5jets_corr.size() > 3) FillHisto(TString("TriggerSel")+TString("_AK5Jet_4_Pt"), (allAK5jets_corr.at(3)).Pt(), evtwt_) ; 
      else FillHisto(TString("TriggerSel")+TString("_AK5Jet_4_Pt"), 0., evtwt_) ; 
    }

    //// Event selection  
    passBSel = false ; 
    if (fatjets_corr.size() >= 1) {
      h_cutflow -> Fill("FatJetSel", 1) ; 
      FillHisto(TString("FatJetSel")+TString("_nJets"), ak5jets_corr.size(), evtwt_) ; 
      FillHisto(TString("FatJetSel")+TString("_nAllAK5"), allAK5jets_corr.size() , evtwt_*puweight_) ; 
      FillHisto(TString("FatJetSel")+TString("_nBJets"), bjets.size(), evtwt_) ; 
      FillHisto(TString("FatJetSel")+TString("_nHJets"), Hjets_corr.size(), evtwt_) ; 
      FillHisto(TString("FatJetSel")+TString("_HTCA8_leading2_AK5_leading2"), HTCA8_leading2_AK5_leading2.getHT(), evtwt_) ; 
      FillHisto(TString("FatJetSel")+TString("_HTAK5_leading4"), HTAK5_leading4.getHT(), evtwt_) ; 
      FillHisto(TString("FatJetSel")+TString("_HTAK5"), HTAK5.getHT(), evtwt_) ; 
      FillHisto(TString("FatJetSel")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt_) ; 
      FillHisto(TString("FatJetSel")+TString("_HT"), MyHT.getHT(), evtwt_) ; 
      for (JetCollection::const_iterator ib = bjets.begin(); ib != bjets.end(); ++ib) { 
        FillHisto(TString("FatJetSel")+TString("_BJet_Pt"), ib->Pt(), evtwt_) ; 
        FillHisto(TString("FatJetSel")+TString("_BJet_Eta"), ib->Eta(), evtwt_) ; 
      }
      FillHisto(TString("FatJetSel")+TString("_nHJetsNoMDyCut")     ,HjetsNoMassDyCuts_corr.size() ,evtwt_)  ; 
      for (JetCollection::const_iterator ijet = HjetsNoMassDyCuts_corr.begin(); ijet != HjetsNoMassDyCuts_corr.end(); ++ijet ) {
        FillHisto(TString("FatJetSel")+TString("_HiggsJetNoMDyCut_Pt")   ,ijet->Pt() ,evtwt_)  ; 
        FillHisto(TString("FatJetSel")+TString("_HiggsJetNoMDyCut_Eta")  ,ijet->Eta() ,evtwt_)  ;
        FillHisto(TString("FatJetSel")+TString("_HiggsJetNoMDyCut_Mass") ,ijet->Mass() ,evtwt_)  ;
      }
      for (JetCollection::const_iterator ijet = fatjets_corr.begin(); ijet != fatjets_corr.end(); ++ijet ) {
        int iSubJet1(ijet->Jet_SubJet1Idx()), iSubJet2(ijet->Jet_SubJet2Idx()); 
        TLorentzVector fatjet_p4, subj1_p4, subj2_p4;
        fatjet_p4.SetPtEtaPhiM(ijet->Pt(), ijet->Eta(), ijet->Phi (), ijet->Mass ()); 
        subj1_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet1], SubJetInfo.Eta[iSubJet1], SubJetInfo.Phi[iSubJet1], SubJetInfo.Mass[iSubJet1]);
        subj2_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet2], SubJetInfo.Eta[iSubJet2], SubJetInfo.Phi[iSubJet2], SubJetInfo.Mass[iSubJet2]); 
        double subjet_dR = subj1_p4.DeltaR(subj2_p4);
        double subjet_dy = subj1_p4.Rapidity() - subj2_p4.Rapidity() ;
        double subjet_dphi = subj1_p4.DeltaPhi(subj2_p4); ;
        double subjet_dyphi = sqrt( subjet_dy*subjet_dy + subjet_dphi*subjet_dphi ) ;

        FillHisto(TString("FatJetSel")+TString("_FatJets_Pt")                 ,fatjet_p4.Pt() ,evtwt_)  ;  
        FillHisto(TString("FatJetSel")+TString("_FatJets_Eta")                ,fatjet_p4.Eta() ,evtwt_)  ; 
        FillHisto(TString("FatJetSel")+TString("_FatJets_CombinedSVBJetTags") ,ijet->CombinedSVBJetTags() ,evtwt_)  ; 
        FillHisto(TString("FatJetSel")+TString("_FatJets_tau2ByTau1")         ,ijet->tau2()/ijet->tau1() ,evtwt_)  ; 
        FillHisto(TString("FatJetSel")+TString("_FatJets_tau3ByTau2")         ,ijet->tau3()/ijet->tau2() ,evtwt_)  ; 
        FillHisto(TString("FatJetSel")+TString("_FatJets_tau3ByTau1")         ,ijet->tau3()/ijet->tau1() ,evtwt_)  ; 
        FillHisto(TString("FatJetSel")+TString("_FatJets_Mass")               ,fatjet_p4.Mag() ,evtwt_)  ; 
        FillHisto(TString("FatJetSel")+TString("_FatJets_MassPruned")         ,ijet->MassPruned() ,evtwt_)  ; 
        FillHisto(TString("FatJetSel")+TString("_FatJets_DRSubjets")          ,subjet_dR ,evtwt_)  ; 
        FillHisto(TString("FatJetSel")+TString("_FatJets_DYPhiSubjets")       ,subjet_dyphi ,evtwt_)  ; 

        FillHisto(TString("FatJetSel")+TString("_SubJet1_Pt") ,subj1_p4.Pt() ,evtwt_)  ;  
        FillHisto(TString("FatJetSel")+TString("_SubJet1_Eta") ,subj1_p4.Eta() ,evtwt_)  ; 
        FillHisto(TString("FatJetSel")+TString("_SubJet1_Mass") ,subj1_p4.Mag() ,evtwt_)  ; 
        FillHisto(TString("FatJetSel")+TString("_SubJet1_CombinedSVBJetTags") ,SubJetInfo.CombinedSVBJetTags[iSubJet1] ,evtwt_)  ; 

        FillHisto(TString("FatJetSel")+TString("_SubJet2_Pt") ,subj2_p4.Pt() ,evtwt_)  ;  
        FillHisto(TString("FatJetSel")+TString("_SubJet2_Eta") ,subj2_p4.Eta() ,evtwt_)  ; 
        FillHisto(TString("FatJetSel")+TString("_SubJet2_Mass") ,subj2_p4.Mag() ,evtwt_)  ; 
        FillHisto(TString("FatJetSel")+TString("_SubJet2_CombinedSVBJetTags") ,SubJetInfo.CombinedSVBJetTags[iSubJet2] ,evtwt_)  ; 

        FillHisto(TString("FatJetSel")+TString("_HiggsJetNoMDyCut_DRsubjets"), subjet_dyphi, evtwt_) ; 
      }
      for (JetCollection::const_iterator ih = Hjets_corr.begin(); ih != Hjets_corr.end(); ++ih) { 
        FillHisto(TString("FatJetSel")+TString("_HiggsJet_Pt"), ih->Pt(), evtwt_) ; 
        FillHisto(TString("FatJetSel")+TString("_HiggsJet_Eta"), ih->Eta(), evtwt_) ; 
        FillHisto(TString("FatJetSel")+TString("_HiggsJet_Mass") ,ih->Mass() ,evtwt_)  ;
      }
      if (allAK5jets_corr.size() > 0) FillHisto(TString("FatJetSel")+TString("_AK5Jet_1_Pt"), (allAK5jets_corr.at(0)).Pt(), evtwt_) ; 
      else FillHisto(TString("FatJetSel")+TString("_AK5Jet_1_Pt"), 0., evtwt_) ; 
      if (allAK5jets_corr.size() > 1) FillHisto(TString("FatJetSel")+TString("_AK5Jet_2_Pt"), (allAK5jets_corr.at(1)).Pt(), evtwt_) ; 
      else FillHisto(TString("FatJetSel")+TString("_AK5Jet_2_Pt"), 0., evtwt_) ; 
      if (allAK5jets_corr.size() > 2) FillHisto(TString("FatJetSel")+TString("_AK5Jet_3_Pt"), (allAK5jets_corr.at(2)).Pt(), evtwt_) ; 
      else FillHisto(TString("FatJetSel")+TString("_AK5Jet_3_Pt"), 0., evtwt_) ; 
      if (allAK5jets_corr.size() > 3) FillHisto(TString("FatJetSel")+TString("_AK5Jet_4_Pt"), (allAK5jets_corr.at(3)).Pt(), evtwt_) ; 
      else FillHisto(TString("FatJetSel")+TString("_AK5Jet_4_Pt"), 0., evtwt_) ; 

      HTSelector htsel(HTSelParams_) ; 
      pat::strbitset retht = htsel.getBitTemplate() ; 
      retht.set(false) ; 

      if ( allAK5jets_corr.size() > 0 && bjets.size() >= 1 ) { //// Special "N-1" selection for Higgs mass and dy, after all other cuts.
        FillHisto(TString("BJetsSel")+TString("_nHJets")     ,Hjets_corr.size() ,evtwt_)  ; 
        FillHisto(TString("BJetsSel")+TString("_nHJetsNoMDyCut")     ,HjetsNoMassDyCuts_corr.size() ,evtwt_)  ; 
        for (JetCollection::const_iterator ijet = HjetsNoMassDyCuts_corr.begin(); ijet != HjetsNoMassDyCuts_corr.end(); ++ijet ) {
          FillHisto(TString("BJetsSel")+TString("_HiggsJetNoMDyCut_Pt")   ,ijet->Pt() ,evtwt_)  ; 
          FillHisto(TString("BJetsSel")+TString("_HiggsJetNoMDyCut_Eta")  ,ijet->Eta() ,evtwt_)  ;
          FillHisto(TString("BJetsSel")+TString("_HiggsJetNoMDyCut_Mass") ,ijet->Mass() ,evtwt_)  ;
          int iSubJet1(ijet->Jet_SubJet1Idx()), iSubJet2(ijet->Jet_SubJet2Idx()); 
          TLorentzVector subj1_p4, subj2_p4;
          subj1_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet1], SubJetInfo.Eta[iSubJet1], SubJetInfo.Phi[iSubJet1], SubJetInfo.Mass[iSubJet1]);
          subj2_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet2], SubJetInfo.Eta[iSubJet2], SubJetInfo.Phi[iSubJet2], SubJetInfo.Mass[iSubJet2]); 
          double subjet_dy = subj1_p4.Rapidity() - subj2_p4.Rapidity() ;
          double subjet_dphi = subj1_p4.DeltaPhi(subj2_p4); ;
          double subjet_dyphi = sqrt( subjet_dy*subjet_dy + subjet_dphi*subjet_dphi ) ;

          FillHisto(TString("BJetsSel")+TString("_SubJet1_Pt") ,subj1_p4.Pt() ,evtwt_)  ;  
          FillHisto(TString("BJetsSel")+TString("_SubJet1_Eta") ,subj1_p4.Eta() ,evtwt_)  ; 
          FillHisto(TString("BJetsSel")+TString("_SubJet1_Mass") ,subj1_p4.Mag() ,evtwt_)  ; 
          FillHisto(TString("BJetsSel")+TString("_SubJet1_CombinedSVBJetTags") ,SubJetInfo.CombinedSVBJetTags[iSubJet1] ,evtwt_)  ; 

          FillHisto(TString("BJetsSel")+TString("_SubJet2_Pt") ,subj2_p4.Pt() ,evtwt_)  ;  
          FillHisto(TString("BJetsSel")+TString("_SubJet2_Eta") ,subj2_p4.Eta() ,evtwt_)  ; 
          FillHisto(TString("BJetsSel")+TString("_SubJet2_Mass") ,subj2_p4.Mag() ,evtwt_)  ; 
          FillHisto(TString("BJetsSel")+TString("_SubJet2_CombinedSVBJetTags") ,SubJetInfo.CombinedSVBJetTags[iSubJet2] ,evtwt_)  ; 

          FillHisto(TString("BJetsSel")+TString("_HiggsJetNoMDyCut_DRsubjets"), subjet_dyphi, evtwt_) ; 
        }
        if ( htsel(HTAllAK5, retht) != 0 ) { //// Special "N-1" selection for Higgs mass and dy, after all other cuts.
          FillHisto(TString("HTSel")+TString("_nHJetsNoMDyCut")     ,HjetsNoMassDyCuts_corr.size() ,evtwt_)  ; 
          for (JetCollection::const_iterator ijet = HjetsNoMassDyCuts_corr.begin(); ijet != HjetsNoMassDyCuts_corr.end(); ++ijet ) {
            FillHisto(TString("HTSel")+TString("_HiggsJetNoMDyCut_Pt")   ,ijet->Pt() ,evtwt_)  ; 
            FillHisto(TString("HTSel")+TString("_HiggsJetNoMDyCut_Eta")  ,ijet->Eta() ,evtwt_)  ;
            FillHisto(TString("HTSel")+TString("_HiggsJetNoMDyCut_Mass") ,ijet->Mass() ,evtwt_)  ;
            int iSubJet1(ijet->Jet_SubJet1Idx()), iSubJet2(ijet->Jet_SubJet2Idx()); 
            TLorentzVector subj1_p4, subj2_p4;
            subj1_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet1], SubJetInfo.Eta[iSubJet1], SubJetInfo.Phi[iSubJet1], SubJetInfo.Mass[iSubJet1]);
            subj2_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet2], SubJetInfo.Eta[iSubJet2], SubJetInfo.Phi[iSubJet2], SubJetInfo.Mass[iSubJet2]); 
            double subjet_dy = subj1_p4.Rapidity() - subj2_p4.Rapidity() ;
            double subjet_dphi = subj1_p4.DeltaPhi(subj2_p4); ;
            double subjet_dyphi = sqrt( subjet_dy*subjet_dy + subjet_dphi*subjet_dphi ) ;
            FillHisto(TString("HTSel")+TString("_HiggsJetNoMDyCut_DRsubjets"), subjet_dyphi, evtwt_) ; 
          }
        }
      } //// Special "N-1" selection for Higgs mass and dy, after all other cuts.

      if ( doABCDPlots_ ) {
        if (bjets.size() >= 1 ) { 
          FillHisto(TString("BJetsSel")+TString("_nAntiHjets"), AntiHjets_corr.size(), evtwt_) ; 
          if ( Hjets_corr.size() >= 1 ) {
            FillHisto(TString("BJetsSel")+TString("_Htag_HTAK5"), HTAK5.getHT(), 1., evtwt_) ; 
            FillHisto(TString("BJetsSel")+TString("_Htag_HT"), MyHT.getHT(), 1., evtwt_) ; 
          }
          else if ( AntiHjets_corr.size() > 0 ) {
            FillHisto(TString("BJetsSel")+TString("_Htag_HTAK5"), MyHT.getHT(), 0., evtwt_) ; 
            FillHisto(TString("BJetsSel")+TString("_Htag_HT"), MyHT.getHT(), 0., evtwt_) ; 
          }
        }
        if (bjets.size() >= 2 ) { 
          FillHisto(TString("BJetsSel")+TString("_nAntiHjets"), AntiHjets_corr.size(), evtwt_) ; 
          if ( Hjets_corr.size() >= 1 ) {
            FillHisto(TString("BJetsSel")+TString("_2b_Htag_HTAK5"), HTAK5.getHT(), 1., evtwt_) ; 
            FillHisto(TString("BJetsSel")+TString("_2b_Htag_HT"), MyHT.getHT(), 1., evtwt_) ; 
          }
          else if ( AntiHjets_corr.size() > 0 ) {
            FillHisto(TString("BJetsSel")+TString("_2b_Htag_HTAK5"), MyHT.getHT(), 0., evtwt_) ; 
            FillHisto(TString("BJetsSel")+TString("_2b_Htag_HT"), MyHT.getHT(), 0., evtwt_) ; 
          }
        }
      }

      if (Hjets_corr.size() >= 1) {
        h_cutflow -> Fill("HiggsJetSel", 1) ; 
        FillHisto(TString("HiggsJetSel")+TString("_nJets"), ak5jets_corr.size(), evtwt_) ; 
        FillHisto(TString("HiggsJetSel")+TString("_nAllAK5"), allAK5jets_corr.size() , evtwt_) ; 
        FillHisto(TString("HiggsJetSel")+TString("_nBJets"), bjets.size(), evtwt_) ; 
        FillHisto(TString("HiggsJetSel")+TString("_nHJets"), Hjets_corr.size(), evtwt_) ; 
        FillHisto(TString("HiggsJetSel")+TString("_HTCA8_leading2_AK5_leading2"), HTCA8_leading2_AK5_leading2.getHT(), evtwt_) ; 
        FillHisto(TString("HiggsJetSel")+TString("_HTAK5_leading4"), HTAK5_leading4.getHT(), evtwt_) ; 
        FillHisto(TString("HiggsJetSel")+TString("_HTAK5"), HTAK5.getHT(), evtwt_) ; 
        FillHisto(TString("HiggsJetSel")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt_) ; 
        FillHisto(TString("HiggsJetSel")+TString("_HT"), MyHT.getHT(), evtwt_) ; 
        for (JetCollection::const_iterator ib = bjets.begin(); ib != bjets.end(); ++ib) { 
          FillHisto(TString("HiggsJetSel")+TString("_BJet_Pt"), ib->Pt(), evtwt_) ; 
          FillHisto(TString("HiggsJetSel")+TString("_BJet_Eta"), ib->Eta(), evtwt_) ; 
        }
        for (JetCollection::const_iterator ih = Hjets_corr.begin(); ih != Hjets_corr.end(); ++ih) { 
          FillHisto(TString("HiggsJetSel")+TString("_HiggsJet_Pt"), ih->Pt(), evtwt_) ; 
          FillHisto(TString("HiggsJetSel")+TString("_HiggsJet_Eta"), ih->Eta(), evtwt_) ; 
          FillHisto(TString("HiggsJetSel")+TString("_HiggsJet_Mass") ,ih->Mass() ,evtwt_)  ;
        }
        if (allAK5jets_corr.size() > 0) FillHisto(TString("HiggsJetSel")+TString("_AK5Jet_1_Pt"), (allAK5jets_corr.at(0)).Pt(), evtwt_) ; 
        else FillHisto(TString("HiggsJetSel")+TString("_AK5Jet_1_Pt"), 0., evtwt_) ; 
        if (allAK5jets_corr.size() > 1) FillHisto(TString("HiggsJetSel")+TString("_AK5Jet_2_Pt"), (allAK5jets_corr.at(1)).Pt(), evtwt_) ; 
        else FillHisto(TString("HiggsJetSel")+TString("_AK5Jet_2_Pt"), 0., evtwt_) ; 
        if (allAK5jets_corr.size() > 2) FillHisto(TString("HiggsJetSel")+TString("_AK5Jet_3_Pt"), (allAK5jets_corr.at(2)).Pt(), evtwt_) ; 
        else FillHisto(TString("HiggsJetSel")+TString("_AK5Jet_3_Pt"), 0., evtwt_) ; 
        if (allAK5jets_corr.size() > 3) FillHisto(TString("HiggsJetSel")+TString("_AK5Jet_4_Pt"), (allAK5jets_corr.at(3)).Pt(), evtwt_) ; 
        else FillHisto(TString("HiggsJetSel")+TString("_AK5Jet_4_Pt"), 0., evtwt_) ; 

        if (allAK5jets_corr.size() > 0) { 
          h_cutflow -> Fill("JetSel", 1) ; 
          FillHisto(TString("JetSel")+TString("_nJets"), ak5jets_corr.size(), evtwt_) ; 
          FillHisto(TString("JetSel")+TString("_nAllAK5"), allAK5jets_corr.size() , evtwt_) ; 
          FillHisto(TString("JetSel")+TString("_nBJets"), bjets.size(), evtwt_) ; 
          FillHisto(TString("JetSel")+TString("_nHJets"), Hjets_corr.size(), evtwt_) ; 
          FillHisto(TString("JetSel")+TString("_HTCA8_leading2_AK5_leading2"), HTCA8_leading2_AK5_leading2.getHT(), evtwt_) ; 
          FillHisto(TString("JetSel")+TString("_HTAK5_leading4"), HTAK5_leading4.getHT(), evtwt_) ; 
          FillHisto(TString("JetSel")+TString("_HTAK5"), HTAK5.getHT(), evtwt_) ; 
          FillHisto(TString("JetSel")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt_) ; 
          FillHisto(TString("JetSel")+TString("_HT"), MyHT.getHT(), evtwt_) ; 
          for (JetCollection::const_iterator ib = bjets.begin(); ib != bjets.end(); ++ib) { 
            FillHisto(TString("JetSel")+TString("_BJet_Pt"), ib->Pt(), evtwt_) ; 
            FillHisto(TString("JetSel")+TString("_BJet_Eta"), ib->Eta(), evtwt_) ; 
          }
          FillHisto(TString("JetSel")+TString("_nHJetsNoMDyCut")     ,HjetsNoMassDyCuts_corr.size() ,evtwt_)  ; 
          for (JetCollection::const_iterator ijet = HjetsNoMassDyCuts_corr.begin(); ijet != HjetsNoMassDyCuts_corr.end(); ++ijet ) {
            FillHisto(TString("JetSel")+TString("_HiggsJetNoMDyCut_Pt")   ,ijet->Pt() ,evtwt_)  ; 
            FillHisto(TString("JetSel")+TString("_HiggsJetNoMDyCut_Eta")  ,ijet->Eta() ,evtwt_)  ;
            FillHisto(TString("JetSel")+TString("_HiggsJetNoMDyCut_Mass") ,ijet->Mass() ,evtwt_)  ;
            int iSubJet1(ijet->Jet_SubJet1Idx()), iSubJet2(ijet->Jet_SubJet2Idx()); 
            TLorentzVector subj1_p4, subj2_p4;
            subj1_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet1], SubJetInfo.Eta[iSubJet1], SubJetInfo.Phi[iSubJet1], SubJetInfo.Mass[iSubJet1]);
            subj2_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet2], SubJetInfo.Eta[iSubJet2], SubJetInfo.Phi[iSubJet2], SubJetInfo.Mass[iSubJet2]); 
            double subjet_dy = subj1_p4.Rapidity() - subj2_p4.Rapidity() ;
            double subjet_dphi = subj1_p4.DeltaPhi(subj2_p4); ;
            double subjet_dyphi = sqrt( subjet_dy*subjet_dy + subjet_dphi*subjet_dphi ) ;
            FillHisto(TString("JetSel")+TString("_HiggsJetNoMDyCut_DRsubjets"), subjet_dyphi, evtwt_) ; 
          }
          for (JetCollection::const_iterator ih = Hjets_corr.begin(); ih != Hjets_corr.end(); ++ih) { 
            FillHisto(TString("JetSel")+TString("_HiggsJet_Pt"), ih->Pt(), evtwt_) ; 
            FillHisto(TString("JetSel")+TString("_HiggsJet_Eta"), ih->Eta(), evtwt_) ; 
            FillHisto(TString("JetSel")+TString("_HiggsJet_Mass") ,ih->Mass() ,evtwt_)  ;
          }
          if (allAK5jets_corr.size() > 0) FillHisto(TString("JetSel")+TString("_AK5Jet_1_Pt"), (allAK5jets_corr.at(0)).Pt(), evtwt_) ; 
          else FillHisto(TString("JetSel")+TString("_AK5Jet_1_Pt"), 0., evtwt_) ; 
          if (allAK5jets_corr.size() > 1) FillHisto(TString("JetSel")+TString("_AK5Jet_2_Pt"), (allAK5jets_corr.at(1)).Pt(), evtwt_) ; 
          else FillHisto(TString("JetSel")+TString("_AK5Jet_2_Pt"), 0., evtwt_) ; 
          if (allAK5jets_corr.size() > 2) FillHisto(TString("JetSel")+TString("_AK5Jet_3_Pt"), (allAK5jets_corr.at(2)).Pt(), evtwt_) ; 
          else FillHisto(TString("JetSel")+TString("_AK5Jet_3_Pt"), 0., evtwt_) ; 
          if (allAK5jets_corr.size() > 3) FillHisto(TString("JetSel")+TString("_AK5Jet_4_Pt"), (allAK5jets_corr.at(3)).Pt(), evtwt_) ; 
          else FillHisto(TString("JetSel")+TString("_AK5Jet_4_Pt"), 0., evtwt_) ; 

          if (bjets.size() >= 1 ) { 
            h_cutflow -> Fill("BJetsSel", 1) ;
            FillHisto(TString("BJetsSel")+TString("_nJets"), ak5jets_corr.size(), evtwt_) ; 
            FillHisto(TString("BJetsSel")+TString("_nAllAK5"), allAK5jets_corr.size() , evtwt_*puweight_) ; 
            FillHisto(TString("BJetsSel")+TString("_HTCA8_leading2_AK5_leading2"), HTCA8_leading2_AK5_leading2.getHT(), evtwt_) ; 
            FillHisto(TString("BJetsSel")+TString("_HTAK5_leading4"), HTAK5_leading4.getHT(), evtwt_) ; 
            FillHisto(TString("BJetsSel")+TString("_HTAK5"), HTAK5.getHT(), evtwt_) ; 
            FillHisto(TString("BJetsSel")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt_) ; 
            FillHisto(TString("BJetsSel")+TString("_HT"), MyHT.getHT(), evtwt_) ; 
            for (JetCollection::const_iterator ib = bjets.begin(); ib != bjets.end(); ++ib) { 
              FillHisto(TString("BJetsSel")+TString("_BJet_Pt"), ib->Pt(), evtwt_) ; 
              FillHisto(TString("BJetsSel")+TString("_BJet_Eta"), ib->Eta(), evtwt_) ; 
            }
            for (JetCollection::const_iterator ih = Hjets_corr.begin(); ih != Hjets_corr.end(); ++ih) { 
              FillHisto(TString("BJetsSel")+TString("_HiggsJet_Pt"), ih->Pt(), evtwt_) ; 
              FillHisto(TString("BJetsSel")+TString("_HiggsJet_Eta"), ih->Eta(), evtwt_) ; 
              FillHisto(TString("BJetsSel")+TString("_HiggsJet_Mass") ,ih->Mass() ,evtwt_)  ;
            }
            if (allAK5jets_corr.size() > 0) FillHisto(TString("BJetsSel")+TString("_AK5Jet_1_Pt"), (allAK5jets_corr.at(0)).Pt(), evtwt_) ; 
            else FillHisto(TString("BJetsSel")+TString("_AK5Jet_1_Pt"), 0., evtwt_) ; 
            if (allAK5jets_corr.size() > 1) FillHisto(TString("BJetsSel")+TString("_AK5Jet_2_Pt"), (allAK5jets_corr.at(1)).Pt(), evtwt_) ; 
            else FillHisto(TString("BJetsSel")+TString("_AK5Jet_2_Pt"), 0., evtwt_) ; 
            if (allAK5jets_corr.size() > 2) FillHisto(TString("BJetsSel")+TString("_AK5Jet_3_Pt"), (allAK5jets_corr.at(2)).Pt(), evtwt_) ; 
            else FillHisto(TString("BJetsSel")+TString("_AK5Jet_3_Pt"), 0., evtwt_) ; 
            if (allAK5jets_corr.size() > 3) FillHisto(TString("BJetsSel")+TString("_AK5Jet_4_Pt"), (allAK5jets_corr.at(3)).Pt(), evtwt_) ; 
            else FillHisto(TString("BJetsSel")+TString("_AK5Jet_4_Pt"), 0., evtwt_) ; 

            if ( htsel(HTAllAK5, retht) != 0 ) { //// If event passes HT 
              h_cutflow -> Fill("HTSel", 1) ; 
              FillHisto(TString("HTSel")+TString("_nJets"), ak5jets_corr.size(), evtwt_) ; 
              FillHisto(TString("HTSel")+TString("_nAllAK5"), allAK5jets_corr.size() , evtwt_*puweight_) ; 
              FillHisto(TString("HTSel")+TString("_nBJets"), bjets.size(), evtwt_) ; 
              FillHisto(TString("HTSel")+TString("_nHJets"), Hjets_corr.size(), evtwt_) ; 
              FillHisto(TString("HTSel")+TString("_HTCA8_leading2_AK5_leading2"), HTCA8_leading2_AK5_leading2.getHT(), evtwt_) ; 
              FillHisto(TString("HTSel")+TString("_HTAK5_leading4"), HTAK5_leading4.getHT(), evtwt_) ; 
              FillHisto(TString("HTSel")+TString("_HTAK5"), HTAK5.getHT(), evtwt_) ; 
              FillHisto(TString("HTSel")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt_) ; 
              FillHisto(TString("HTSel")+TString("_HT"), MyHT.getHT(), evtwt_) ; 
              for (JetCollection::const_iterator ib = bjets.begin(); ib != bjets.end(); ++ib) { 
                FillHisto(TString("HTSel")+TString("_BJet_Pt"), ib->Pt(), evtwt_) ; 
                FillHisto(TString("HTSel")+TString("_BJet_Eta"), ib->Eta(), evtwt_) ; 
              }
              for (JetCollection::const_iterator ih = Hjets_corr.begin(); ih != Hjets_corr.end(); ++ih) { 
                FillHisto(TString("HTSel")+TString("_HiggsJet_Pt"), ih->Pt(), evtwt_) ; 
                FillHisto(TString("HTSel")+TString("_HiggsJet_Eta"), ih->Eta(), evtwt_) ; 
                FillHisto(TString("HTSel")+TString("_HiggsJet_Mass") ,ih->Mass() ,evtwt_)  ;
              }
              if (allAK5jets_corr.size() > 0) FillHisto(TString("HTSel")+TString("_AK5Jet_1_Pt"), (allAK5jets_corr.at(0)).Pt(), evtwt_) ; 
              else FillHisto(TString("HTSel")+TString("_AK5Jet_1_Pt"), 0., evtwt_) ; 
              if (allAK5jets_corr.size() > 1) FillHisto(TString("HTSel")+TString("_AK5Jet_2_Pt"), (allAK5jets_corr.at(1)).Pt(), evtwt_) ; 
              else FillHisto(TString("HTSel")+TString("_AK5Jet_2_Pt"), 0., evtwt_) ; 
              if (allAK5jets_corr.size() > 2) FillHisto(TString("HTSel")+TString("_AK5Jet_3_Pt"), (allAK5jets_corr.at(2)).Pt(), evtwt_) ; 
              else FillHisto(TString("HTSel")+TString("_AK5Jet_3_Pt"), 0., evtwt_) ; 
              if (allAK5jets_corr.size() > 3) FillHisto(TString("HTSel")+TString("_AK5Jet_4_Pt"), (allAK5jets_corr.at(3)).Pt(), evtwt_) ; 
              else FillHisto(TString("HTSel")+TString("_AK5Jet_4_Pt"), 0., evtwt_) ; 

              //// Reconstruct b' candidates
              for (JetCollection::const_iterator ihig = Hjets_corr.begin(); ihig != Hjets_corr.end(); ++ihig) { 
                unsigned int closestBIndex(bjets.size()) ;
                double deltaR(TMath::Pi()) ; 
                for (JetCollection::const_iterator ib = bjets.begin(); ib != bjets.end(); ++ib) { 
                  if ( ihig->DeltaR(*ib) < deltaR) {
                    deltaR = ihig->DeltaR(*ib) ; 
                    closestBIndex = ib - bjets.begin() ; 
                  }
                }
                if (/*deltaR < TMath::Pi() &&*/ closestBIndex < bjets.size()) {
                  p4bprimes.push_back(ihig->p4() + (bjets.at(closestBIndex)).p4()) ; 
                }
              } //// Reconstruct b' candidates
              FillHisto(TString("HTSel")+TString("_Nbprimes"), p4bprimes.size(), evtwt_) ; 
              for (std::vector<TLorentzVector>::const_iterator ibprime = p4bprimes.begin(); ibprime != p4bprimes.end(); ++ibprime) {
                FillHisto(TString("HTSel")+TString("_bprimePt") ,ibprime->Pt() ,evtwt_)  ; 
                FillHisto(TString("HTSel")+TString("_bprimeMass") ,ibprime->Mag() ,evtwt_)  ;
              } //// Fill all bprime candidates 

              if (Hjets_corr.size() == 1 && bjets.size() == 1) { //// 1 HJet and 1 bjet 
                h_cutflow -> Fill("HTSel_1H_1b", 1) ; 
                FillHisto(TString("HTSel_1H_1b")+TString("_nJets"), ak5jets_corr.size(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_1b")+TString("_nAllAK5"), allAK5jets_corr.size() , evtwt_*puweight_) ; 
                FillHisto(TString("HTSel_1H_1b")+TString("_nBJets"), bjets.size(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_1b")+TString("_nHJets"), Hjets_corr.size(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_1b")+TString("_HTCA8_leading2_AK5_leading2"), HTCA8_leading2_AK5_leading2.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_1b")+TString("_HTAK5_leading4"), HTAK5_leading4.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_1b")+TString("_HTAK5"), HTAK5.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_1b")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_1b")+TString("_HT"), MyHT.getHT(), evtwt_) ; 
                for (JetCollection::const_iterator ib = bjets.begin(); ib != bjets.end(); ++ib) { 
                  FillHisto(TString("HTSel_1H_1b")+TString("_BJet_Pt"), ib->Pt(), evtwt_) ; 
                  FillHisto(TString("HTSel_1H_1b")+TString("_BJet_Eta"), ib->Eta(), evtwt_) ; 
                }
                for (JetCollection::const_iterator ih = Hjets_corr.begin(); ih != Hjets_corr.end(); ++ih) { 
                  FillHisto(TString("HTSel_1H_1b")+TString("_HiggsJet_Pt"), ih->Pt(), evtwt_) ; 
                  FillHisto(TString("HTSel_1H_1b")+TString("_HiggsJet_Eta"), ih->Eta(), evtwt_) ; 
                  FillHisto(TString("HTSel_1H_1b")+TString("_HiggsJet_Mass") ,ih->Mass() ,evtwt_)  ;
                }
                if (allAK5jets_corr.size() > 0) FillHisto(TString("HTSel_1H_1b")+TString("_AK5Jet_1_Pt"), (allAK5jets_corr.at(0)).Pt(), evtwt_) ; 
                else FillHisto(TString("HTSel_1H_1b")+TString("_AK5Jet_1_Pt"), 0., evtwt_) ; 
                if (allAK5jets_corr.size() > 1) FillHisto(TString("HTSel_1H_1b")+TString("_AK5Jet_2_Pt"), (allAK5jets_corr.at(1)).Pt(), evtwt_) ; 
                else FillHisto(TString("HTSel_1H_1b")+TString("_AK5Jet_2_Pt"), 0., evtwt_) ; 
                if (allAK5jets_corr.size() > 2) FillHisto(TString("HTSel_1H_1b")+TString("_AK5Jet_3_Pt"), (allAK5jets_corr.at(2)).Pt(), evtwt_) ; 
                else FillHisto(TString("HTSel_1H_1b")+TString("_AK5Jet_3_Pt"), 0., evtwt_) ; 
                if (allAK5jets_corr.size() > 3) FillHisto(TString("HTSel_1H_1b")+TString("_AK5Jet_4_Pt"), (allAK5jets_corr.at(3)).Pt(), evtwt_) ; 
                else FillHisto(TString("HTSel_1H_1b")+TString("_AK5Jet_4_Pt"), 0., evtwt_) ; 
                FillHisto(TString("HTSel_1H_1b")+TString("_HTCA8_leading2_AK5_leading2"), HTCA8_leading2_AK5_leading2.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_1b")+TString("_HTAK5_leading4"), HTAK5_leading4.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_1b")+TString("_HTAK5"), HTAK5.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_1b")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_1b")+TString("_HT"), MyHT.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_1b")+TString("_Nbprimes"), p4bprimes.size(), evtwt_) ; 
                for (std::vector<TLorentzVector>::const_iterator ibprime = p4bprimes.begin(); ibprime != p4bprimes.end(); ++ibprime) {
                  FillHisto(TString("HTSel_1H_1b")+TString("_bprimePt") ,ibprime->Pt() ,evtwt_)  ; 
                  FillHisto(TString("HTSel_1H_1b")+TString("_bprimeMass") ,ibprime->Mag() ,evtwt_)  ;
                } //// Fill all bprime candidates 
              } //// 1 HJet and 1 bjet 

              if (Hjets_corr.size() == 1 && bjets.size() >= 2) { //// 1 HJet and >=2 bjet 
                h_cutflow -> Fill("HTSel_1H_2b", 1) ; 
                FillHisto(TString("HTSel_1H_2b")+TString("_nJets"), ak5jets_corr.size(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_2b")+TString("_nAllAK5"), allAK5jets_corr.size() , evtwt_*puweight_) ; 
                FillHisto(TString("HTSel_1H_2b")+TString("_nBJets"), bjets.size(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_2b")+TString("_nHJets"), Hjets_corr.size(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_2b")+TString("_HTCA8_leading2_AK5_leading2"), HTCA8_leading2_AK5_leading2.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_2b")+TString("_HTAK5_leading4"), HTAK5_leading4.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_2b")+TString("_HTAK5"), HTAK5.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_2b")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_2b")+TString("_HT"), MyHT.getHT(), evtwt_) ; 
                for (JetCollection::const_iterator ib = bjets.begin(); ib != bjets.end(); ++ib) { 
                  FillHisto(TString("HTSel_1H_2b")+TString("_BJet_Pt"), ib->Pt(), evtwt_) ; 
                  FillHisto(TString("HTSel_1H_2b")+TString("_BJet_Eta"), ib->Eta(), evtwt_) ; 
                }
                for (JetCollection::const_iterator ih = Hjets_corr.begin(); ih != Hjets_corr.end(); ++ih) { 
                  FillHisto(TString("HTSel_1H_2b")+TString("_HiggsJet_Pt"), ih->Pt(), evtwt_) ; 
                  FillHisto(TString("HTSel_1H_2b")+TString("_HiggsJet_Eta"), ih->Eta(), evtwt_) ; 
                  FillHisto(TString("HTSel_1H_2b")+TString("_HiggsJet_Mass") ,ih->Mass() ,evtwt_)  ;
                }
                if (allAK5jets_corr.size() > 0) FillHisto(TString("HTSel_1H_2b")+TString("_AK5Jet_1_Pt"), (allAK5jets_corr.at(0)).Pt(), evtwt_) ; 
                else FillHisto(TString("HTSel_1H_2b")+TString("_AK5Jet_1_Pt"), 0., evtwt_) ; 
                if (allAK5jets_corr.size() > 1) FillHisto(TString("HTSel_1H_2b")+TString("_AK5Jet_2_Pt"), (allAK5jets_corr.at(1)).Pt(), evtwt_) ; 
                else FillHisto(TString("HTSel_1H_2b")+TString("_AK5Jet_2_Pt"), 0., evtwt_) ; 
                if (allAK5jets_corr.size() > 2) FillHisto(TString("HTSel_1H_2b")+TString("_AK5Jet_3_Pt"), (allAK5jets_corr.at(2)).Pt(), evtwt_) ; 
                else FillHisto(TString("HTSel_1H_2b")+TString("_AK5Jet_3_Pt"), 0., evtwt_) ; 
                if (allAK5jets_corr.size() > 3) FillHisto(TString("HTSel_1H_2b")+TString("_AK5Jet_4_Pt"), (allAK5jets_corr.at(3)).Pt(), evtwt_) ; 
                else FillHisto(TString("HTSel_1H_2b")+TString("_AK5Jet_4_Pt"), 0., evtwt_) ; 
                FillHisto(TString("HTSel_1H_2b")+TString("_HTCA8_leading2_AK5_leading2"), HTCA8_leading2_AK5_leading2.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_2b")+TString("_HTAK5_leading4"), HTAK5_leading4.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_2b")+TString("_HTAK5"), HTAK5.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_2b")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_2b")+TString("_HT"), MyHT.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_1H_2b")+TString("_Nbprimes"), p4bprimes.size(), evtwt_) ; 
                for (std::vector<TLorentzVector>::const_iterator ibprime = p4bprimes.begin(); ibprime != p4bprimes.end(); ++ibprime) {
                  FillHisto(TString("HTSel_1H_2b")+TString("_bprimePt") ,ibprime->Pt() ,evtwt_)  ; 
                  FillHisto(TString("HTSel_1H_2b")+TString("_bprimeMass") ,ibprime->Mag() ,evtwt_)  ;
                } //// Fill all bprime candidates 
              } //// 1 HJet and >=2 bjet 

              if (Hjets_corr.size() >= 2 && bjets.size() == 1) { //// >= 2 HJet and 1 bjet 
                h_cutflow -> Fill("HTSel_2H_1b", 1) ; 
                FillHisto(TString("HTSel_2H_1b")+TString("_nJets"), ak5jets_corr.size(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_1b")+TString("_nAllAK5"), allAK5jets_corr.size() , evtwt_*puweight_) ; 
                FillHisto(TString("HTSel_2H_1b")+TString("_nBJets"), bjets.size(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_1b")+TString("_nHJets"), Hjets_corr.size(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_1b")+TString("_HTCA8_leading2_AK5_leading2"), HTCA8_leading2_AK5_leading2.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_1b")+TString("_HTAK5_leading4"), HTAK5_leading4.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_1b")+TString("_HTAK5"), HTAK5.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_1b")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_1b")+TString("_HT"), MyHT.getHT(), evtwt_) ; 
                for (JetCollection::const_iterator ib = bjets.begin(); ib != bjets.end(); ++ib) { 
                  FillHisto(TString("HTSel_2H_1b")+TString("_BJet_Pt"), ib->Pt(), evtwt_) ; 
                  FillHisto(TString("HTSel_2H_1b")+TString("_BJet_Eta"), ib->Eta(), evtwt_) ; 
                }
                for (JetCollection::const_iterator ih = Hjets_corr.begin(); ih != Hjets_corr.end(); ++ih) { 
                  FillHisto(TString("HTSel_2H_1b")+TString("_HiggsJet_Pt"), ih->Pt(), evtwt_) ; 
                  FillHisto(TString("HTSel_2H_1b")+TString("_HiggsJet_Eta"), ih->Eta(), evtwt_) ; 
                  FillHisto(TString("HTSel_2H_1b")+TString("_HiggsJet_Mass") ,ih->Mass() ,evtwt_)  ;
                }
                if (allAK5jets_corr.size() > 0) FillHisto(TString("HTSel_2H_1b")+TString("_AK5Jet_1_Pt"), (allAK5jets_corr.at(0)).Pt(), evtwt_) ; 
                else FillHisto(TString("HTSel_2H_1b")+TString("_AK5Jet_1_Pt"), 0., evtwt_) ; 
                if (allAK5jets_corr.size() > 1) FillHisto(TString("HTSel_2H_1b")+TString("_AK5Jet_2_Pt"), (allAK5jets_corr.at(1)).Pt(), evtwt_) ; 
                else FillHisto(TString("HTSel_2H_1b")+TString("_AK5Jet_2_Pt"), 0., evtwt_) ; 
                if (allAK5jets_corr.size() > 2) FillHisto(TString("HTSel_2H_1b")+TString("_AK5Jet_3_Pt"), (allAK5jets_corr.at(2)).Pt(), evtwt_) ; 
                else FillHisto(TString("HTSel_2H_1b")+TString("_AK5Jet_3_Pt"), 0., evtwt_) ; 
                if (allAK5jets_corr.size() > 3) FillHisto(TString("HTSel_2H_1b")+TString("_AK5Jet_4_Pt"), (allAK5jets_corr.at(3)).Pt(), evtwt_) ; 
                else FillHisto(TString("HTSel_2H_1b")+TString("_AK5Jet_4_Pt"), 0., evtwt_) ; 
                FillHisto(TString("HTSel_2H_1b")+TString("_HTCA8_leading2_AK5_leading2"), HTCA8_leading2_AK5_leading2.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_1b")+TString("_HTAK5_leading4"), HTAK5_leading4.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_1b")+TString("_HTAK5"), HTAK5.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_1b")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_1b")+TString("_HT"), MyHT.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_1b")+TString("_Nbprimes"), p4bprimes.size(), evtwt_) ; 
                for (std::vector<TLorentzVector>::const_iterator ibprime = p4bprimes.begin(); ibprime != p4bprimes.end(); ++ibprime) {
                  FillHisto(TString("HTSel_2H_1b")+TString("_bprimePt") ,ibprime->Pt() ,evtwt_)  ; 
                  FillHisto(TString("HTSel_2H_1b")+TString("_bprimeMass") ,ibprime->Mag() ,evtwt_)  ;
                } //// Fill all bprime candidates 
              } //// >= 2 HJet and 1 bjet 

              if (Hjets_corr.size() >= 2 && bjets.size() >= 2) { //// >= 2 HJet and >=2 bjet  
                h_cutflow -> Fill("HTSel_2H_2b", 1) ; 
                FillHisto(TString("HTSel_2H_2b")+TString("_nJets"), ak5jets_corr.size(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_2b")+TString("_nAllAK5"), allAK5jets_corr.size() , evtwt_*puweight_) ; 
                FillHisto(TString("HTSel_2H_2b")+TString("_nBJets"), bjets.size(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_2b")+TString("_nHJets"), Hjets_corr.size(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_2b")+TString("_HTCA8_leading2_AK5_leading2"), HTCA8_leading2_AK5_leading2.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_2b")+TString("_HTAK5_leading4"), HTAK5_leading4.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_2b")+TString("_HTAK5"), HTAK5.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_2b")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_2b")+TString("_HT"), MyHT.getHT(), evtwt_) ; 
                for (JetCollection::const_iterator ib = bjets.begin(); ib != bjets.end(); ++ib) { 
                  FillHisto(TString("HTSel_2H_2b")+TString("_BJet_Pt"), ib->Pt(), evtwt_) ; 
                  FillHisto(TString("HTSel_2H_2b")+TString("_BJet_Eta"), ib->Eta(), evtwt_) ; 
                }
                for (JetCollection::const_iterator ih = Hjets_corr.begin(); ih != Hjets_corr.end(); ++ih) { 
                  FillHisto(TString("HTSel_2H_2b")+TString("_HiggsJet_Pt"), ih->Pt(), evtwt_) ; 
                  FillHisto(TString("HTSel_2H_2b")+TString("_HiggsJet_Eta"), ih->Eta(), evtwt_) ; 
                  FillHisto(TString("HTSel_2H_2b")+TString("_HiggsJet_Mass") ,ih->Mass() ,evtwt_)  ;
                }
                if (allAK5jets_corr.size() > 0) FillHisto(TString("HTSel_2H_2b")+TString("_AK5Jet_1_Pt"), (allAK5jets_corr.at(0)).Pt(), evtwt_) ; 
                else FillHisto(TString("HTSel_2H_2b")+TString("_AK5Jet_1_Pt"), 0., evtwt_) ; 
                if (allAK5jets_corr.size() > 1) FillHisto(TString("HTSel_2H_2b")+TString("_AK5Jet_2_Pt"), (allAK5jets_corr.at(1)).Pt(), evtwt_) ; 
                else FillHisto(TString("HTSel_2H_2b")+TString("_AK5Jet_2_Pt"), 0., evtwt_) ; 
                if (allAK5jets_corr.size() > 2) FillHisto(TString("HTSel_2H_2b")+TString("_AK5Jet_3_Pt"), (allAK5jets_corr.at(2)).Pt(), evtwt_) ; 
                else FillHisto(TString("HTSel_2H_2b")+TString("_AK5Jet_3_Pt"), 0., evtwt_) ; 
                if (allAK5jets_corr.size() > 3) FillHisto(TString("HTSel_2H_2b")+TString("_AK5Jet_4_Pt"), (allAK5jets_corr.at(3)).Pt(), evtwt_) ; 
                else FillHisto(TString("HTSel_2H_2b")+TString("_AK5Jet_4_Pt"), 0., evtwt_) ; 
                FillHisto(TString("HTSel_2H_2b")+TString("_HTCA8_leading2_AK5_leading2"), HTCA8_leading2_AK5_leading2.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_2b")+TString("_HTAK5_leading4"), HTAK5_leading4.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_2b")+TString("_HTAK5"), HTAK5.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_2b")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_2b")+TString("_HT"), MyHT.getHT(), evtwt_) ; 
                FillHisto(TString("HTSel_2H_2b")+TString("_Nbprimes"), p4bprimes.size(), evtwt_) ; 
                for (std::vector<TLorentzVector>::const_iterator ibprime = p4bprimes.begin(); ibprime != p4bprimes.end(); ++ibprime) {
                  FillHisto(TString("HTSel_2H_2b")+TString("_bprimePt") ,ibprime->Pt() ,evtwt_)  ; 
                  FillHisto(TString("HTSel_2H_2b")+TString("_bprimeMass") ,ibprime->Mag() ,evtwt_)  ;
                } //// Fill all bprime candidates 
              } //// >= 2 HJet and >=2 bjet  
            } //// If event passes HT 
            passBSel = true ; 
          } //// If at least one b-jets 
        } //// If at least one AK5 jet 
      } //// If at least one Higgs jet 
    } //// If at least one fat jet 

    if ( doTrigEff_ && passBSel && passHLT) { 
      h_cutflow -> Fill("TriggerSel", 1) ; 
      FillHisto(TString("TriggerSel")+TString("_nPVtx_NoPUWt"), nGoodVtxs, evtwt_) ; 
      FillHisto(TString("TriggerSel")+TString("_nPVtx_PUWt"), nGoodVtxs, evtwt_*puweight_) ; 
      FillHisto(TString("TriggerSel")+TString("_nJets"), ak5jets_corr.size(), evtwt_) ; 
      FillHisto(TString("TriggerSel")+TString("_nFatJets"), fatjets_corr.size(), evtwt_) ; 
      FillHisto(TString("TriggerSel")+TString("_nAllAK5"), allAK5jets_corr.size() , evtwt_) ; 
      FillHisto(TString("TriggerSel")+TString("_nBJets"), bjets.size(), evtwt_) ; 
      FillHisto(TString("TriggerSel")+TString("_nHJets"), Hjets_corr.size(), evtwt_) ; 
      FillHisto(TString("TriggerSel")+TString("_HTCA8_leading2_AK5_leading2"), HTCA8_leading2_AK5_leading2.getHT(), evtwt_) ; 
      FillHisto(TString("TriggerSel")+TString("_HTAK5_leading4"), HTAK5_leading4.getHT(), evtwt_) ; 
      FillHisto(TString("TriggerSel")+TString("_HTAK5"), HTAK5.getHT(), evtwt_) ; 
      FillHisto(TString("TriggerSel")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt_) ; 
      FillHisto(TString("TriggerSel")+TString("_HT"), MyHT.getHT(), evtwt_) ; 
      for (JetCollection::const_iterator ib = bjets.begin(); ib != bjets.end(); ++ib) { 
        FillHisto(TString("TriggerSel")+TString("_BJet_Pt"), ib->Pt(), evtwt_) ; 
        FillHisto(TString("TriggerSel")+TString("_BJet_Eta"), ib->Eta(), evtwt_) ; 
      }

      if (allAK5jets_corr.size() > 0) FillHisto(TString("TriggerSel")+TString("_AK5Jet_1_Pt"), (allAK5jets_corr.at(0)).Pt(), evtwt_) ; 
      else FillHisto(TString("TriggerSel")+TString("_AK5Jet_1_Pt"), 0., evtwt_) ; 
      if (allAK5jets_corr.size() > 1) FillHisto(TString("TriggerSel")+TString("_AK5Jet_2_Pt"), (allAK5jets_corr.at(1)).Pt(), evtwt_) ; 
      else FillHisto(TString("TriggerSel")+TString("_AK5Jet_2_Pt"), 0., evtwt_) ; 
      if (allAK5jets_corr.size() > 2) FillHisto(TString("TriggerSel")+TString("_AK5Jet_3_Pt"), (allAK5jets_corr.at(2)).Pt(), evtwt_) ; 
      else FillHisto(TString("TriggerSel")+TString("_AK5Jet_3_Pt"), 0., evtwt_) ; 
      if (allAK5jets_corr.size() > 3) FillHisto(TString("TriggerSel")+TString("_AK5Jet_4_Pt"), (allAK5jets_corr.at(3)).Pt(), evtwt_) ; 
      else FillHisto(TString("TriggerSel")+TString("_AK5Jet_4_Pt"), 0., evtwt_) ; 
    }

  } //// entry loop 

}

// ------------ method called once each job just after ending the event loop  ------------
void BprimeTobHAnalysis::endJob() { 
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void BprimeTobHAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BprimeTobHAnalysis);
