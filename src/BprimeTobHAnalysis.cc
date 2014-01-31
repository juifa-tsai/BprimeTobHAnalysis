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

// Root headers 
#include <TLorentzVector.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH1I.h>
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
    template <class Type>
      void FillHisto(const TString& name, const Type value, const double weight);

    edm::LumiReWeighting LumiWeights_; 

    // ----------member data ---------------------------

    //// Configurables 

    int                             maxEvents_; 
    const int                       reportEvery_; 
    const std::string               inputTTree_;
    const std::vector<std::string>  inputFiles_;
    //const std::vector<string>       hltPaths_; 
    const edm::ParameterSet         hltPaths_; 
    const int                       doPUReweighting_ ;
    const std::string               file_PUDistMC_ ;
    const std::string               file_PUDistData_ ;
    const std::string               hist_PUDistMC_ ;
    const std::string               hist_PUDistData_ ;

    const double jetPtMin_ ; 
    const double jetPtMax_ ; 
    const double jetAbsEtaMax_ ;
    const double bjetPtMin_ ; 
    const double fatJetPtMin_ ; 
    const double fatJetPtMax_ ; 
    const double fatJetAbsEtaMax_ ;
    const double fatJetMassMin_ ;
    const double fatJetMassMax_ ; 
    const double fatJetPrunedMassMin_ ;
    const double fatJetPrunedMassMax_ ; 
    const double fatJetTau2ByTau1Max_ ; 
    const double subjet1CSVDiscMin_ ; 
    const double subjet1CSVDiscMax_ ; 
    const double subjet2CSVDiscMin_ ; 
    const double subjet2CSVDiscMax_ ; 
    const double HTMin_ ; 
    const double HTMax_ ; 
    const edm::ParameterSet jetSelParams_ ; 
    const edm::ParameterSet fatJetSelParams_ ; 
    const edm::ParameterSet higgsJetSelParams_ ; 
    const edm::ParameterSet HTSelParams_ ; 
    const edm::ParameterSet evtSelParams_ ; 
    const edm::ParameterSet jmeParams_; 
    const double jesShift_;
    const double jerShift_; 
    const double SFbShift_;
    const double SFlShift_;

    TChain*            chain_;

    EvtInfoBranches    EvtInfo;
    VertexInfoBranches VtxInfo;
    GenInfoBranches    GenInfo;
    JetInfoBranches    GenJetInfo;
    JetInfoBranches    JetInfo;
    JetInfoBranches    FatJetInfo;
    JetInfoBranches    SubJetInfo;
    LepInfoBranches    LepInfo;

    edm::Service<TFileService> fs; 

    bool isData_ ; 
    double evtwt_ ; 
    double puweight_ ; 

    TH1D* h_cutflow ; 

    std::map<TString, TH1D*> hmap_1d ;  

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
  //hltPaths_(iConfig.getParameter<std::vector<string> >("HLTPaths")),
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
  fatJetPtMin_(iConfig.getParameter<double>("FatJetPtMin")),
  fatJetPtMax_(iConfig.getParameter<double>("FatJetPtMax")),
  fatJetAbsEtaMax_(iConfig.getParameter<double>("FatJetAbsEtaMax")),
  fatJetMassMin_(iConfig.getParameter<double>("FatJetMassMin")),
  fatJetMassMax_(iConfig.getParameter<double>("FatJetMassMax")), 
  fatJetPrunedMassMin_(iConfig.getParameter<double>("FatJetPrunedMassMin")),
  fatJetPrunedMassMax_(iConfig.getParameter<double>("FatJetPrunedMassMax")),
  fatJetTau2ByTau1Max_(iConfig.getParameter<double>("FatJetTau2ByTau1Max")),
  subjet1CSVDiscMin_(iConfig.getParameter<double>("Subjet1CSVDiscMin")),
  subjet1CSVDiscMax_(iConfig.getParameter<double>("Subjet1CSVDiscMax")),
  subjet2CSVDiscMin_(iConfig.getParameter<double>("Subjet2CSVDiscMin")),
  subjet2CSVDiscMax_(iConfig.getParameter<double>("Subjet2CSVDiscMax")),
  HTMin_(iConfig.getParameter<double>("HTMin")), 
  HTMax_(iConfig.getParameter<double>("HTMax")),
  jetSelParams_(iConfig.getParameter<edm::ParameterSet>("JetSelParams")), 
  fatJetSelParams_(iConfig.getParameter<edm::ParameterSet>("FatJetSelParams")), 
  higgsJetSelParams_(iConfig.getParameter<edm::ParameterSet>("HiggsJetSelParams")), 
  HTSelParams_(iConfig.getParameter<edm::ParameterSet>("HTSelParams")), 
  evtSelParams_(iConfig.getParameter<edm::ParameterSet>("EvtSelParams")), 
  jmeParams_(iConfig.getParameter<edm::ParameterSet>("JMEParams")),
  jesShift_(iConfig.getParameter<double>("JESShift")),
  jerShift_(iConfig.getParameter<double>("JERShift")),
  SFbShift_(iConfig.getParameter<double>("SFbShift")),
  SFlShift_(iConfig.getParameter<double>("SFlShift")),
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

  chain_ = new TChain(inputTTree_.c_str());

  for(unsigned i=0; i<inputFiles_.size(); ++i) {
    chain_->Add(inputFiles_.at(i).c_str());

    TFile *f = TFile::Open(inputFiles_.at(i).c_str(),"READ");
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

  h_cutflow = fs->make<TH1D>("h_cutflow" ,"Cut flow" ,20 ,0. ,20.) ;  
  h_cutflow -> Sumw2() ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(1,"AllEvents") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(2,"TriggerSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(3,"FatJetSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(4,"HiggsJetSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(5,"JetSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(6,"BJetsSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(7,"HTSel") ; 

  for (int ii = 1; ii <= h_cutflow->GetNbinsX(); ++ii) 
    CreateHistos(h_cutflow->GetXaxis()->GetBinLabel(ii)) ; 

  return ;  

}

void BprimeTobHAnalysis::CreateHistos(const TString& cutname) {

  AddHisto(cutname ,"_nPVtx_NoPUWt"                ,"N(PV), No PU weight"                          ,50     ,-0.5     ,49.5    ) ; 
  AddHisto(cutname ,"_nPVtx_PUWt"                  ,"N(PV)"                                        ,50     ,-0.5     ,49.5    ) ; 
  AddHisto(cutname ,"_nJets"                       ,"N(AK5 jets)"                                  ,20     ,-0.5     ,19.5    ) ; 
  AddHisto(cutname ,"_nAllAK5"                     ,"N(All AK5 jets)"                              ,20     ,-0.5     ,19.5    ) ; 
  AddHisto(cutname ,"_nBJets"                      ,"N(b jets)"                                    ,20     ,-0.5     ,19.5    ) ; 
  AddHisto(cutname ,"_nFatJets"                    ,"N(fat jets)"                                  ,20     ,-0.5     ,19.5    ) ; 
  AddHisto(cutname ,"_nHJets"                      ,"N(Higgs jets)"                                ,20     ,-0.5     ,19.5    ) ; 
  AddHisto(cutname ,"_FatJets_Pt"                  ,"p_{T}(fat jets)"                              ,1000   ,0.       ,1000.   ) ; 
  AddHisto(cutname ,"_FatJets_Eta"                 ,"#eta(fat jets)"                               ,50     ,-4.      ,4.      ) ; 
  AddHisto(cutname ,"_FatJets_Mass"                ,"Fat jet mass [GeV]"                           ,100    ,0.       ,2000.   ) ; 
  AddHisto(cutname ,"_FatJets_MassPruned"          ,"Fat jet pruned mass [GeV]"                    ,100    ,0.       ,2000.   ) ; 
  AddHisto(cutname ,"_FatJets_tau2ByTau1"          ,"Fat jet #tau2/#tau1"                          ,20     ,0.       ,1.      ) ; 
  AddHisto(cutname ,"_FatJets_tau3ByTau2"          ,"Fat jet #tau2/#tau1"                          ,20     ,0.       ,1.      ) ; 
  AddHisto(cutname ,"_FatJets_tau3ByTau1"          ,"Fat jet #tau2/#tau1"                          ,20     ,0.       ,1.      ) ; 
  AddHisto(cutname ,"_FatJets_CombinedSVBJetTags"  ,"Fat jet CSV discriminator"                    ,20     ,0.       ,1.      ) ; 

  AddHisto(cutname ,"_SubJet1_Pt"                  ,"SubJet1 p_{T} [GeV]"                          ,100    ,0.       ,2000.   ) ;
  AddHisto(cutname ,"_SubJet1_Eta"                 ,"SubJet1 #eta"                                 ,50     ,-4.      ,4.      ) ;
  AddHisto(cutname ,"_SubJet1_Mass"                ,"SubJet1 mass [GeV]"                           ,100    ,0.       ,2000.   ) ; 
  AddHisto(cutname ,"_SubJet1_CombinedSVBJetTags"  ,"SubJet1 CSV discriminator"                    ,20     ,0.       ,1.      ) ; 

  AddHisto(cutname ,"_SubJet2_Pt"                  ,"SubJet2 p_{T} [GeV]"                          ,100    ,0.       ,2000.   ) ;
  AddHisto(cutname ,"_SubJet2_Eta"                 ,"SubJet2 #eta"                                 ,50     ,-4.      ,4.      ) ;
  AddHisto(cutname ,"_SubJet2_Mass"                ,"SubJet2 mass [GeV]"                           ,100    ,0.       ,2000.   ) ; 
  AddHisto(cutname ,"_SubJet2_CombinedSVBJetTags"  ,"SubJet2 CSV discriminator"                    ,20     ,0.       ,1.      ) ; 

  AddHisto(cutname ,"_HiggsJet_Pt"                 ,"p_{T} (Higgs jet)[GeV]"                       ,100    ,0.       ,2000.   ) ;
  AddHisto(cutname ,"_HiggsJet_Eta"                ,"#eta (Higgs jet)"                             ,50     ,-4.      ,4.      ) ;
  AddHisto(cutname ,"_HiggsJet_Mass"               ,"Mass (Higgs jet) [GeV]"                       ,100    ,0.       ,200.    ) ; 

  AddHisto(cutname ,"_BJet_Pt"                     ,"p_{T} (b jet)[GeV]"                           ,100    ,0.       ,2000.   ) ;
  AddHisto(cutname ,"_BJet_Eta"                    ,"#eta (b jet)"                                 ,50     ,-4.      ,4.      ) ;

  AddHisto(cutname ,"_HTCA8_leading2_AK5_leading2" ,"H_{T} (leading 2 CA8 and AK5 jets) [GeV]"     ,200   ,0.       ,4000.   ) ;
  AddHisto(cutname ,"_HTAK5_leading4"              ,"H_{T} (leading 4 AK5 jets) [GeV]"             ,200   ,0.       ,4000.   ) ;
  AddHisto(cutname ,"_HTAK5"                       ,"H_{T} (AK5 jets, no Higgs jet overlap) [GeV]" ,200   ,0.       ,4000.   ) ;
  AddHisto(cutname ,"_HTAllAK5"                    ,"H_{T} (AK5 jets) [GeV]"                       ,200   ,0.       ,4000.   ) ;
  AddHisto(cutname ,"_HT"                          ,"H_{T}[GeV]"                                   ,200   ,0.       ,4000.   ) ;
  AddHisto(cutname ,"_bprimePt"                    ,"b' p_{T} [GeV]"                               ,100   ,0.       ,2000.   ) ;
  AddHisto(cutname ,"_bprimeMass"                  ,"b' mass [GeV]"                                ,40    ,0.       ,2000.   ) ;

  return ; 

}

void BprimeTobHAnalysis::AddHisto(const TString& cutname, const TString& histname, const TString& histtitle, const int& nbins, const double& min, const double& max) { 

  TH1D* h1d ; 
  h1d = fs->make<TH1D>(cutname+histname, cutname+histtitle, nbins, min, max);  
  h1d -> Sumw2() ; 
  hmap_1d[cutname+histname] = h1d ; 

  return ; 

}

template <class Type>
void BprimeTobHAnalysis::FillHisto(const TString& name, const Type value, const double weight){
  hmap_1d[name]->Fill(double(value),weight);

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

  JetSelector jetSelAK5(jetSelParams_) ; 
  FatJetSelector jetSelCA8(fatJetSelParams_) ; 
  FatJetSelector jetSelHiggs(higgsJetSelParams_) ; 
  pat::strbitset retjetidak5 = jetSelAK5.getBitTemplate() ; 
  pat::strbitset retjetidca8 = jetSelCA8.getBitTemplate() ; 

  ofstream fout("Evt_NoJets.txt") ; 
  if ( isData_ ) {
    fout << "EvtInfo.RunNo " << " EvtInfo.LumiNo " << " EvtInfo.EvtNo " << std::endl ;
  }

  edm::LogInfo("StartingAnalysisLoop") << "Starting analysis loop\n";

  for(int entry=0; entry<maxEvents_; entry++) {

    if((entry%reportEvery_) == 0) edm::LogInfo("Event") << entry << " of " << maxEvents_ ; 

    //// Event variables 
    std::vector<TLorentzVector>p4_fatJets ; 
    std::vector<TLorentzVector>p4_higgsJets ; 
    std::vector<TLorentzVector>p4_bJets ; 
    std::vector<TLorentzVector>bprimes ; 

    bool passHLT(false) ; 

    int nGoodVtxs(0) ;
    int njets(0) ; 
    int nfatjets(0) ; 
    double H_TCA8_leading2_AK5_leading2(0) ; 
    double H_TAK5_leading4(0) ; 
    double H_TAK5(0) ; 
    double H_TAllAK5(0) ; 
    double H_T(0) ; 

    chain_->GetEntry(entry);

    isData_   = EvtInfo.McFlag ? 0 : 1; 
    if ( !isData_ ) evtwt_    = EvtInfo.Weight ; 
    if ( doPUReweighting_ && !isData_ ) puweight_ = LumiWeights_.weight(EvtInfo.TrueIT[0]) ; 

    //// Higgs BR reweighting
    if ( !isData_ ) {

      int nbprimeBH(0), nbprimeBZ(0), nbprimeTW(0) ; 

      for (int igen=0; igen < GenInfo.Size; ++igen) {
        if ( GenInfo.Status[igen] == 3 &&
            TMath::Abs(GenInfo.PdgID[igen]) == 25 && 
            GenInfo.nDa[igen] >= 2 ) { //// H found
          int higgsDau0 = abs(GenInfo.Da0PdgID[igen]) ;
          int higgsDau1 = abs(GenInfo.Da1PdgID[igen]) ;
          int higgsDau = higgsDau0+higgsDau1 ; 
          if(higgsDau==10) evtwt_ *= HiggsBRscaleFactors::higgsBBSf;
          if(higgsDau==30) evtwt_ *= HiggsBRscaleFactors::higgsTauTauSf;
          if(higgsDau==26) evtwt_ *= HiggsBRscaleFactors::higgsMuMuSf;
          if(higgsDau==8)  evtwt_ *= HiggsBRscaleFactors::higgsCCSf;
          if(higgsDau==6)  evtwt_ *= HiggsBRscaleFactors::higgsSSSf;
          if(higgsDau==16) evtwt_ *= HiggsBRscaleFactors::higgsTTSf;
          if(higgsDau==42) evtwt_ *= HiggsBRscaleFactors::higgsGGSf;
          if(higgsDau==44) evtwt_ *= HiggsBRscaleFactors::higgsGammaGammaSf;
          if(higgsDau==45) evtwt_ *= HiggsBRscaleFactors::higgsZGammaSf;
          if(higgsDau==48) evtwt_ *= HiggsBRscaleFactors::higgsWWSf;
          if(higgsDau==46) evtwt_ *= HiggsBRscaleFactors::higgsZZSf; 
        }
      }

    }

    nGoodVtxs = 0 ;
    VertexSelector vtxSel(VtxInfo) ; 
    nGoodVtxs = vtxSel.NGoodVtxs(); 
    if (nGoodVtxs < 1)  { edm::LogInfo("NoGoodPrimaryVertex") << " No good primary vertex " ; continue ; }

    FillHisto(TString("AllEvents")+TString("_nPVtx_NoPUWt"), nGoodVtxs, evtwt_) ; 
    FillHisto(TString("AllEvents")+TString("_nPVtx_PUWt"), nGoodVtxs, evtwt_*puweight_) ; 
    FillHisto(TString("AllEvents")+TString("_nJets"), JetInfo.Size, evtwt_*puweight_) ; 
    FillHisto(TString("AllEvents")+TString("_nFatJets"), FatJetInfo.Size, evtwt_*puweight_) ; 
    h_cutflow -> Fill("AllEvents", 1) ; 

    TriggerSelector trigSel(hltPaths_) ; 
    passHLT = trigSel.getTrigDecision(EvtInfo) ; 

    if ( !passHLT ) continue ; 

    if ( isData_ ) {
      if ( JetInfo.Size == 0 ) fout << EvtInfo.RunNo << " " << EvtInfo.LumiNo << " " << EvtInfo.EvtNo << std::endl ; 
    }

    FillHisto(TString("TriggerSel")+TString("_nPVtx_NoPUWt"), nGoodVtxs, evtwt_) ; 
    FillHisto(TString("TriggerSel")+TString("_nPVtx_PUWt"), nGoodVtxs, evtwt_*puweight_) ; 
    h_cutflow -> Fill("TriggerSel", 1) ; 

    evtwt_ *= puweight_ ; 

    JetCollection allAK5jets_corr ; //// All AK5 jets, including those overlapping with Higgs jets 
    JetCollection ak5jets_corr ; 
    JetCollection ak5jets_leading4 ; 
    JetCollection bjets ; 
    JetCollection fatjets ; 
    JetCollection higgsjets ; 
    JetCollection jets_CA8_leading2_AK5_leading2 ; 

    for (int ifatjet=0; ifatjet < FatJetInfo.Size; ++ifatjet) {


      //DM 31Jan if ( FatJetInfo.Pt[ifatjet] < fatJetPtMin_ 
      //DM 31Jan     || FatJetInfo.Pt[ifatjet] > fatJetPtMax_ ) continue; //// apply jet pT cut
      //DM 31Jan if ( fabs(FatJetInfo.Eta[ifatjet]) > fatJetAbsEtaMax_ ) continue; //// apply jet eta cut
      //DM 31Jan retca8.set(false);
      //DM 31Jan if ( fatjetIDLoose(FatJetInfo, ifatjet,retca8) == 0 ) continue; //// apply loose jet ID
      if (jetSelCA8(FatJetInfo, ifatjet, SubJetInfo, retjetidca8) == 0) continue ; 

      TLorentzVector fatjet_p4;
      fatjet_p4.SetPtEtaPhiM(FatJetInfo.Pt[ifatjet], FatJetInfo.Eta[ifatjet], 
          FatJetInfo.Phi[ifatjet], FatJetInfo.Mass[ifatjet]);

      FillHisto(TString("TriggerSel")+TString("_FatJets_Pt")                 ,fatjet_p4.Pt() ,evtwt_)  ;  
      FillHisto(TString("TriggerSel")+TString("_FatJets_Eta")                ,fatjet_p4.Eta() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_FatJets_Mass")               ,fatjet_p4.Mag() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_FatJets_MassPruned")         ,FatJetInfo.MassPruned[ifatjet] ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_FatJets_CombinedSVBJetTags") ,FatJetInfo.CombinedSVBJetTags[ifatjet] ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_FatJets_tau2ByTau1")         ,FatJetInfo.tau2[ifatjet]/FatJetInfo.tau1[ifatjet] ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_FatJets_tau3ByTau2")         ,FatJetInfo.tau3[ifatjet]/FatJetInfo.tau2[ifatjet] ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_FatJets_tau3ByTau1")         ,FatJetInfo.tau3[ifatjet]/FatJetInfo.tau1[ifatjet] ,evtwt_)  ; 

      //// Get subjets of fat jets 

      int iSubJet1 = FatJetInfo.Jet_SubJet1Idx[ifatjet];
      int iSubJet2 = FatJetInfo.Jet_SubJet2Idx[ifatjet];

      if( SubJetInfo.Pt[iSubJet1]==0. || SubJetInfo.Pt[iSubJet2]==0. ) 
        continue; //// skip fat jets for which one of the subjets has pT=0

      TLorentzVector subjet1_p4, subjet2_p4;
      subjet1_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet1], SubJetInfo.Eta[iSubJet1], 
          SubJetInfo.Phi[iSubJet1], SubJetInfo.Mass[iSubJet1]);
      subjet2_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet2], SubJetInfo.Eta[iSubJet2], 
          SubJetInfo.Phi[iSubJet2], SubJetInfo.Mass[iSubJet2]);

      double subjet_dR = subjet1_p4.DeltaR(subjet2_p4);
      double subjet_dy = subjet1_p4.Rapidity() - subjet2_p4.Rapidity() ;
      double subjet_dphi = subjet1_p4.DeltaPhi(subjet2_p4); ;
      double subjet_dyphi = sqrt( subjet_dy*subjet_dy + subjet_dphi*subjet_dphi ) ;

      //if( subjet_dyphi < max(double(FatJetInfo.Mass[ifatjet]/FatJetInfo.Pt[ifatjet]), 0.4) ) 
      //  continue; //// skip infrared unsafe configurations and overlap between the subjets 

      FillHisto(TString("TriggerSel")+TString("_SubJet1_Pt") ,subjet1_p4.Pt() ,evtwt_)  ;  
      FillHisto(TString("TriggerSel")+TString("_SubJet1_Eta") ,subjet1_p4.Eta() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_SubJet1_Mass") ,subjet1_p4.Mag() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_SubJet1_CombinedSVBJetTags") ,SubJetInfo.CombinedSVBJetTags[iSubJet1] ,evtwt_)  ; 

      FillHisto(TString("TriggerSel")+TString("_SubJet2_Pt") ,subjet2_p4.Pt() ,evtwt_)  ;  
      FillHisto(TString("TriggerSel")+TString("_SubJet2_Eta") ,subjet2_p4.Eta() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_SubJet2_Mass") ,subjet2_p4.Mag() ,evtwt_)  ; 
      FillHisto(TString("TriggerSel")+TString("_SubJet2_CombinedSVBJetTags") ,SubJetInfo.CombinedSVBJetTags[iSubJet2] ,evtwt_)  ; 

      //// Selecting fat jets 
      if (fatjet_p4.Mag() > fatJetMassMin_ 
          && fatjet_p4.Mag() < fatJetMassMax_ 
          && FatJetInfo.tau2[ifatjet]/FatJetInfo.tau1[ifatjet] < fatJetTau2ByTau1Max_) {

        Jet thisjet(FatJetInfo, ifatjet) ; 
        fatjets.push_back(thisjet) ; 
        if (nfatjets < 2) jets_CA8_leading2_AK5_leading2.push_back(thisjet) ;  

        p4_fatJets.push_back(fatjet_p4) ; 

        ++nfatjets ; 

        //// Higgs tagging
        if ( SubJetInfo.CombinedSVBJetTags[iSubJet1] > subjet1CSVDiscMin_ 
            && SubJetInfo.CombinedSVBJetTags[iSubJet1] < subjet1CSVDiscMax_ 
            && SubJetInfo.CombinedSVBJetTags[iSubJet2] > subjet2CSVDiscMin_ 
            && SubJetInfo.CombinedSVBJetTags[iSubJet2] < subjet2CSVDiscMax_) { 
          FillHisto(TString("TriggerSel")+TString("_HiggsJet_Pt")   ,fatjet_p4.Pt() ,evtwt_)  ; 
          FillHisto(TString("TriggerSel")+TString("_HiggsJet_Eta")  ,fatjet_p4.Eta() ,evtwt_)  ;
          FillHisto(TString("TriggerSel")+TString("_HiggsJet_Mass") ,fatjet_p4.Mag() ,evtwt_)  ;

          higgsjets.push_back(thisjet) ; 

          p4_higgsJets.push_back(fatjet_p4) ; 

          if ( !isData_ ) { //// Apply Higgs-tagging scale factor 
            ApplyHiggsTagSF* higgsTagSF = new ApplyHiggsTagSF(double(SubJetInfo.Pt[iSubJet1]), double(SubJetInfo.Pt[iSubJet2]), 
                double(SubJetInfo.Eta[iSubJet1]), double(SubJetInfo.Eta[iSubJet2]),
                SubJetInfo.GenFlavor[iSubJet1], SubJetInfo.GenFlavor[iSubJet2], 
                SubJetInfo.CombinedSVBJetTags[iSubJet1], SubJetInfo.CombinedSVBJetTags[iSubJet2]) ; 
            evtwt_ *= higgsTagSF->GetHiggsTagSF() ;
            delete higgsTagSF ; 
          }

        } //// Higgs tagging 

      } //// Selecting fat jets 

    } //// Loop over fat jets 

    JetCollection myjets ;
    JetCollection allMyJets ;
    for (int ijet = 0; ijet < JetInfo.Size; ++ijet) {
      retjetidak5.set(false) ;
      if (jetSelAK5(JetInfo, ijet,retjetidak5) == 0) continue ; 
      Jet thisjet(JetInfo, ijet) ;
      allMyJets.push_back(thisjet) ; 
      bool isJetNotHiggs(false) ; 
      TLorentzVector jet_p4;
      jet_p4.SetPtEtaPhiM(JetInfo.Pt[ijet], JetInfo.Eta[ijet], JetInfo.Phi[ijet], JetInfo.Mass[ijet]); 
      for (std::vector<TLorentzVector>::const_iterator ihig = p4_higgsJets.begin(); ihig != p4_higgsJets.end(); ++ihig) {
        if (jet_p4.DeltaR(*ihig) < 1.2) { //// Attn. hard-coded 
          isJetNotHiggs = false ; 
          break ; 
        }
        else {
          isJetNotHiggs = true ; 
        } 
      } 
      if (isJetNotHiggs) myjets.push_back(thisjet) ; 
    }

    JetCollection ak5jets_jec, allAK5jets_jec ; 
    if ( !isData_) {
      //// Only AK5 jets not overlapping with Higgs jets 
      JMEUncertUtil* jmeUtil_jer = new JMEUncertUtil(jmeParams_, EvtInfo, myjets, "JER", jerShift_) ; 
      JetCollection ak5jets_jer = jmeUtil_jer->GetModifiedJetColl() ; 
      delete jmeUtil_jer ; 

      JMEUncertUtil* jmeUtil_jes = new JMEUncertUtil(jmeParams_, EvtInfo, ak5jets_jer, "JESAK5MC", jesShift_) ; 
      JetCollection ak5jets_jec = jmeUtil_jes->GetModifiedJetColl() ; 
      delete jmeUtil_jes ; 

      //// All AK5 jets 
      jmeUtil_jer = new JMEUncertUtil(jmeParams_, EvtInfo, allMyJets, "JER", jerShift_) ; 
      JetCollection allAK5jets_jer = jmeUtil_jer->GetModifiedJetColl() ; 
      delete jmeUtil_jer ; 

      jmeUtil_jes = new JMEUncertUtil(jmeParams_, EvtInfo, allAK5jets_jer, "JESAK5MC", jesShift_) ; 
      JetCollection allAK5jets_jec = jmeUtil_jes->GetModifiedJetColl() ; 
      delete jmeUtil_jes ; 

    }
    else {
      //// Only AK5 jets not overlapping with Higgs jets 
      JMEUncertUtil* jmeUtil_jes = new JMEUncertUtil(jmeParams_, EvtInfo, myjets, "JESAK5DATA", jesShift_) ; 
      ak5jets_jec = jmeUtil_jes->GetModifiedJetColl() ; 
      delete jmeUtil_jes ; 

      //// All AK5 jets 
      jmeUtil_jes = new JMEUncertUtil(jmeParams_, EvtInfo, allMyJets, "JESAK5DATA", jesShift_) ; 
      allAK5jets_jec = jmeUtil_jes->GetModifiedJetColl() ; 
      delete jmeUtil_jes ; 
    }

    for (JetCollection::const_iterator ijet = ak5jets_jec.begin(); ijet != ak5jets_jec.end(); ++ijet) {
      if (ijet->Pt() < jetPtMin_ || ijet->Pt() > jetPtMax_) continue ; 
      ak5jets_corr.push_back(*ijet) ; 

      if (njets < 4) ak5jets_leading4.push_back(*ijet) ; 
      if (njets < 2) jets_CA8_leading2_AK5_leading2.push_back(*ijet) ;  
      ++njets ; 

    }

    for (JetCollection::const_iterator ijet = allAK5jets_jec.begin(); ijet != allAK5jets_jec.end(); ++ijet) {
      if (ijet->Pt() < jetPtMin_ || ijet->Pt() > jetPtMax_) continue ; 
      allAK5jets_corr.push_back(*ijet) ; 
    }

    JetCollection ak5jets_btag ; 
    if ( !isData_ ) {
      ApplyBTagSF * btagsf =  new ApplyBTagSF(ak5jets_corr, 0.679, "CSVM", SFbShift_, SFlShift_) ;  
      ak5jets_btag =  btagsf->getBtaggedJetsWithSF () ; 
      delete btagsf ; 
    }
    else {
      for (JetCollection::const_iterator ijet = allAK5jets_jec.begin(); ijet != allAK5jets_jec.end(); ++ijet) {
        ak5jets_btag.push_back(*ijet) ; 
      }
    }

    for (JetCollection::const_iterator ijet = ak5jets_btag.begin(); ijet != ak5jets_btag.end(); ++ijet) { 

      if (ijet->Pt() > bjetPtMin_ && ijet->CombinedSVBJetTags() > 0.679) { //// Select b-tagged AK5 jets  
        bjets.push_back(*ijet) ; 
        TLorentzVector bjet_p4;
        bjet_p4.SetPtEtaPhiM(ijet->Pt(), ijet->Eta(), ijet->Phi(), ijet->Mass()); 
        p4_bJets.push_back(bjet_p4) ; 
      } //// Select b-tagged AK5 jets 
    }

    HT HTAK5 ; 
    HTAK5.setJetCollection(ak5jets_corr) ; 
    HTAK5.buildHT() ; 
    H_TAK5 = HTAK5.getHT() ; 

    HT HTAllAK5 ; 
    HTAllAK5.setJetCollection(allAK5jets_corr) ; 
    HTAllAK5.buildHT() ; 
    H_TAllAK5 = HTAllAK5.getHT() ; 

    HT HTAK5_leading4 ; 
    HTAK5_leading4.setJetCollection(ak5jets_leading4) ; 
    HTAK5_leading4.buildHT() ; 
    H_TAK5_leading4 = HTAK5_leading4.getHT() ; 

    HT HTCA8_leading2_AK5_leading2 ; 
    HTCA8_leading2_AK5_leading2.setJetCollection(jets_CA8_leading2_AK5_leading2) ; 
    HTCA8_leading2_AK5_leading2.buildHT() ; 
    H_TCA8_leading2_AK5_leading2 = HTCA8_leading2_AK5_leading2.getHT() ; 

    HT MyHT ; 
    MyHT.setJetCollection(higgsjets) ; 
    MyHT.setJetCollection(bjets) ; 
    MyHT.buildHT() ; 
    H_T = MyHT.getHT() ; 

    FillHisto(TString("TriggerSel")+TString("_nJets"), njets, evtwt_) ; 
    FillHisto(TString("TriggerSel")+TString("_nFatJets"), fatjets.size(), evtwt_*puweight_) ; 
    FillHisto(TString("TriggerSel")+TString("_nAllAK5"), allAK5jets_corr.size() , evtwt_*puweight_) ; 
    FillHisto(TString("TriggerSel")+TString("_nBJets"), p4_bJets.size(), evtwt_) ; 
    FillHisto(TString("TriggerSel")+TString("_HTCA8_leading2_AK5_leading2"), H_TCA8_leading2_AK5_leading2, evtwt_) ; 
    FillHisto(TString("TriggerSel")+TString("_HTAK5_leading4"), H_TAK5_leading4, evtwt_) ; 
    FillHisto(TString("TriggerSel")+TString("_HTAK5"), H_TAK5, evtwt_) ; 
    FillHisto(TString("TriggerSel")+TString("_HTAllAK5"), H_TAllAK5, evtwt_) ; 
    FillHisto(TString("TriggerSel")+TString("_HT"), H_T, evtwt_) ; 
    for (std::vector<TLorentzVector>::const_iterator ib = p4_bJets.begin(); ib != p4_bJets.end(); ++ib) { 
      FillHisto(TString("TriggerSel")+TString("_BJet_Pt"), ib->Pt(), evtwt_) ; 
      FillHisto(TString("TriggerSel")+TString("_BJet_Eta"), ib->Eta(), evtwt_) ; 
    }

    if (p4_fatJets.size() >= 1) {
      FillHisto(TString("FatJetSel")+TString("_nJets"), njets, evtwt_) ; 
      FillHisto(TString("FatJetSel")+TString("_nAllAK5"), allAK5jets_corr.size() , evtwt_*puweight_) ; 
      FillHisto(TString("FatJetSel")+TString("_nBJets"), p4_bJets.size(), evtwt_) ; 
      FillHisto(TString("FatJetSel")+TString("_HTCA8_leading2_AK5_leading2"), H_TCA8_leading2_AK5_leading2, evtwt_) ; 
      FillHisto(TString("FatJetSel")+TString("_HTAK5_leading4"), H_TAK5_leading4, evtwt_) ; 
      FillHisto(TString("FatJetSel")+TString("_HTAK5"), H_TAK5, evtwt_) ; 
      FillHisto(TString("FatJetSel")+TString("_HTAllAK5"), H_TAllAK5, evtwt_) ; 
      FillHisto(TString("FatJetSel")+TString("_HT"), H_T, evtwt_) ; 
      for (std::vector<TLorentzVector>::const_iterator ib = p4_bJets.begin(); ib != p4_bJets.end(); ++ib) { 
        FillHisto(TString("FatJetSel")+TString("_BJet_Pt"), ib->Pt(), evtwt_) ; 
        FillHisto(TString("FatJetSel")+TString("_BJet_Eta"), ib->Eta(), evtwt_) ; 
      }
      for (std::vector<TLorentzVector>::const_iterator ih = p4_higgsJets.begin(); ih != p4_higgsJets.end(); ++ih) { 
        FillHisto(TString("FatJetSel")+TString("_HiggsJet_Pt"), ih->Pt(), evtwt_) ; 
        FillHisto(TString("FatJetSel")+TString("_HiggsJet_Eta"), ih->Eta(), evtwt_) ; 
        FillHisto(TString("FatJetSel")+TString("_HiggsJet_Mass") ,ih->Mag() ,evtwt_)  ;
      }
      h_cutflow -> Fill("FatJetSel", 1) ; 
    }

    if (higgsjets.size() >= 1) {
      FillHisto(TString("HiggsJetSel")+TString("_nJets"), njets, evtwt_) ; 
      FillHisto(TString("HiggsJetSel")+TString("_nAllAK5"), allAK5jets_corr.size() , evtwt_*puweight_) ; 
      FillHisto(TString("HiggsJetSel")+TString("_nBJets"), p4_bJets.size(), evtwt_) ; 
      FillHisto(TString("HiggsJetSel")+TString("_HTCA8_leading2_AK5_leading2"), H_TCA8_leading2_AK5_leading2, evtwt_) ; 
      FillHisto(TString("HiggsJetSel")+TString("_HTAK5_leading4"), H_TAK5_leading4, evtwt_) ; 
      FillHisto(TString("HiggsJetSel")+TString("_HTAK5"), H_TAK5, evtwt_) ; 
      FillHisto(TString("HiggsJetSel")+TString("_HTAllAK5"), H_TAllAK5, evtwt_) ; 
      FillHisto(TString("HiggsJetSel")+TString("_HT"), H_T, evtwt_) ; 
      for (std::vector<TLorentzVector>::const_iterator ib = p4_bJets.begin(); ib != p4_bJets.end(); ++ib) { 
        FillHisto(TString("HiggsJetSel")+TString("_BJet_Pt"), ib->Pt(), evtwt_) ; 
        FillHisto(TString("HiggsJetSel")+TString("_BJet_Eta"), ib->Eta(), evtwt_) ; 
      }
      for (std::vector<TLorentzVector>::const_iterator ih = p4_higgsJets.begin(); ih != p4_higgsJets.end(); ++ih) { 
        FillHisto(TString("HiggsJetSel")+TString("_HiggsJet_Pt"), ih->Pt(), evtwt_) ; 
        FillHisto(TString("HiggsJetSel")+TString("_HiggsJet_Eta"), ih->Eta(), evtwt_) ; 
        FillHisto(TString("HiggsJetSel")+TString("_HiggsJet_Mass") ,ih->Mag() ,evtwt_)  ;
      }
      h_cutflow -> Fill("HiggsJetSel", 1) ; 

      if (p4_bJets.size() >= 2 ) { 
        FillHisto(TString("BJetsSel")+TString("_nJets"), njets, evtwt_) ; 
        FillHisto(TString("BJetsSel")+TString("_nAllAK5"), allAK5jets_corr.size() , evtwt_*puweight_) ; 
        FillHisto(TString("BJetsSel")+TString("_HTCA8_leading2_AK5_leading2"), H_TCA8_leading2_AK5_leading2, evtwt_) ; 
        FillHisto(TString("BJetsSel")+TString("_HTAK5_leading4"), H_TAK5_leading4, evtwt_) ; 
        FillHisto(TString("BJetsSel")+TString("_HTAK5"), H_TAK5, evtwt_) ; 
        FillHisto(TString("BJetsSel")+TString("_HTAllAK5"), H_TAllAK5, evtwt_) ; 
        FillHisto(TString("BJetsSel")+TString("_HT"), H_T, evtwt_) ; 
        for (std::vector<TLorentzVector>::const_iterator ib = p4_bJets.begin(); ib != p4_bJets.end(); ++ib) { 
          FillHisto(TString("BJetsSel")+TString("_BJet_Pt"), ib->Pt(), evtwt_) ; 
          FillHisto(TString("BJetsSel")+TString("_BJet_Eta"), ib->Eta(), evtwt_) ; 
        }
        for (std::vector<TLorentzVector>::const_iterator ih = p4_higgsJets.begin(); ih != p4_higgsJets.end(); ++ih) { 
          FillHisto(TString("BJetsSel")+TString("_HiggsJet_Pt"), ih->Pt(), evtwt_) ; 
          FillHisto(TString("BJetsSel")+TString("_HiggsJet_Eta"), ih->Eta(), evtwt_) ; 
          FillHisto(TString("BJetsSel")+TString("_HiggsJet_Mass") ,ih->Mag() ,evtwt_)  ;
        }
        h_cutflow -> Fill("BJetsSel", 1) ;

        HTSelector htsel(HTSelParams_) ; 
        pat::strbitset retht = htsel.getBitTemplate() ; 
        retht.set(false) ; 
        if ( htsel(HTAllAK5, retht) == 0 ) continue ; 

        FillHisto(TString("HTSel")+TString("_nJets"), njets, evtwt_) ; 
        FillHisto(TString("HTSel")+TString("_nAllAK5"), allAK5jets_corr.size() , evtwt_*puweight_) ; 
        FillHisto(TString("HTSel")+TString("_nBJets"), p4_bJets.size(), evtwt_) ; 
        FillHisto(TString("HTSel")+TString("_nHJets"), p4_higgsJets.size(), evtwt_) ; 
        FillHisto(TString("HTSel")+TString("_HTCA8_leading2_AK5_leading2"), H_TCA8_leading2_AK5_leading2, evtwt_) ; 
        FillHisto(TString("HTSel")+TString("_HTAK5_leading4"), H_TAK5_leading4, evtwt_) ; 
        FillHisto(TString("HTSel")+TString("_HTAK5"), H_TAK5, evtwt_) ; 
        FillHisto(TString("HTSel")+TString("_HTAllAK5"), H_TAllAK5, evtwt_) ; 
        FillHisto(TString("HTSel")+TString("_HT"), H_T, evtwt_) ; 
        for (std::vector<TLorentzVector>::const_iterator ib = p4_bJets.begin(); ib != p4_bJets.end(); ++ib) { 
          FillHisto(TString("HTSel")+TString("_BJet_Pt"), ib->Pt(), evtwt_) ; 
          FillHisto(TString("HTSel")+TString("_BJet_Eta"), ib->Eta(), evtwt_) ; 
        }
        for (std::vector<TLorentzVector>::const_iterator ih = p4_higgsJets.begin(); ih != p4_higgsJets.end(); ++ih) { 
          FillHisto(TString("HTSel")+TString("_HiggsJet_Pt"), ih->Pt(), evtwt_) ; 
          FillHisto(TString("HTSel")+TString("_HiggsJet_Eta"), ih->Eta(), evtwt_) ; 
          FillHisto(TString("HTSel")+TString("_HiggsJet_Mass") ,ih->Mag() ,evtwt_)  ;
        }
        h_cutflow -> Fill("HTSel", 1) ; 

        //// Reconstruct b' candidates
        for (std::vector<TLorentzVector>::const_iterator ihig = p4_higgsJets.begin(); ihig != p4_higgsJets.end(); ++ihig) { 
          const TLorentzVector* closestB_p4 ;
          double deltaR(TMath::Pi()) ; 
          for (std::vector<TLorentzVector>::const_iterator ib = p4_bJets.begin(); ib != p4_bJets.end(); ++ib) { 
            if ( ihig->DeltaR(*ib) < deltaR) {
              deltaR = ihig->DeltaR(*ib) ; 
              closestB_p4 = &(*ib) ; 
            }
          }
          if (deltaR < TMath::Pi()) {
            bprimes.push_back(*ihig + *closestB_p4) ; 
            FillHisto(TString("HTSel")+TString("_bprimePt") ,(*ihig + *closestB_p4).Pt() ,evtwt_)  ; 
            FillHisto(TString("HTSel")+TString("_bprimeMass") ,(*ihig + *closestB_p4).Mag() ,evtwt_)  ;
          }
        } //// Reconstruct b' candidates

      } //// If at least two b-jets 
    } //// If at least one Higgs jet 

  } //// entry loop 

  fout.close() ; 

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
