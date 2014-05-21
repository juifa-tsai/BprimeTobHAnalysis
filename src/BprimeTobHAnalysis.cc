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

    void isolateCollection( JetCollection control, JetCollection input, JetCollection& output, double dR=1.2 );
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

    const double bJetPtMin_ ; 
    const double bJetCSVDiscMin_;
    const double bJetCSVDiscMax_;
    const double fatJetPtMin_ ;
    const double fatJetPtMax_ ; 
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
    const bool   doPDFTree_;
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

    edm::Service<TFileService> fs; 
    TH1D *h_bprimeDecayModes, *h_cutflow ; 
    std::map<TString, TH1D*> hmap_1d ;  
    std::map<TString, TH2D*> hmap_2d ;  

    TTree* pdfTree_ ; 
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
  bJetPtMin_(iConfig.getParameter<double>("BJetPtMin")),
  bJetCSVDiscMin_(iConfig.getParameter<double>("BJetCSVDiscMin")),
  bJetCSVDiscMax_(iConfig.getParameter<double>("BJetCSVDiscMax")),
  fatJetPtMin_(iConfig.getParameter<double>("FatJetPtMin")),
  fatJetPtMax_(iConfig.getParameter<double>("FatJetPtMax")), 
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
  doPDFTree_(iConfig.getParameter<double>("DoPDFTree")),
  fillBDTTrees_(iConfig.getParameter<double>("FillBDTTrees"))  
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

  int ncutflowbins(13) ; 
  h_cutflow = fs->make<TH1D>("h_cutflow" ,"Cut flow" ,ncutflowbins, 0 ,ncutflowbins) ;  
  h_cutflow -> Sumw2() ; 

  h_cutflow -> GetXaxis() -> SetBinLabel(1 ,"AllEvents") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(2 ,"SkimmedEvents") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(3 ,"TriggerSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(4 ,"VertexSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(5 ,"FatJetSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(6 ,"HiggsJetSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(7 ,"JetSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(8 ,"BJetsSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(9 ,"HTSel") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(10,"HTSel_1H_1b") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(11,"HTSel_1H_2b") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(12,"HTSel_2H_1b") ; 
  h_cutflow -> GetXaxis() -> SetBinLabel(13,"HTSel_2H_2b") ; 

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

void BprimeTobHAnalysis::isolateCollection( JetCollection control, JetCollection input, JetCollection& output, double dR ){
  for( JetCollection::const_iterator ii = input.begin(); ii != input.end(); ++ii ){
    bool isControl(false); 
    for( JetCollection::const_iterator ic = control.begin(); ic != control.end(); ++ic ){
      if( ic->DeltaR(*ii) < dR ){
        isControl = true; 
        break; 
      }else {
        isControl = false; 
      } 
    }
    if( !isControl ){
      output.push_back(*ii);
    } 
  }
}

void BprimeTobHAnalysis::CreateHistos(const TString& cutname) {

  AddHisto(cutname ,"_nPVtx_NoPUWt"                ,"N(PV), No PU weight"                          ,50     ,-0.5     ,49.5    ) ; 
  AddHisto(cutname ,"_nPVtx_PUWt"                  ,"N(PV)"                                        ,50     ,-0.5     ,49.5    ) ; 

  AddHisto(cutname ,"_nFatJets"                    ,"N(fat jets)"                                  ,20     ,-0.5     ,19.5    ) ; 
  AddHisto(cutname ,"_FatJets_Pt"                  ,"p_{T}(fat jets)"                              ,1000   ,0.       ,1000.   ) ; 
  AddHisto(cutname ,"_FatJets_Eta"                 ,"#eta(fat jets)"                               ,50     ,-4.      ,4.      ) ; 
  AddHisto(cutname ,"_FatJets_tau2ByTau1"          ,"Fat jet #tau2/#tau1"                          ,20     ,0.       ,1.      ) ; 
  AddHisto(cutname ,"_FatJets_tau3ByTau2"          ,"Fat jet #tau2/#tau1"                          ,20     ,0.       ,1.      ) ; 
  AddHisto(cutname ,"_FatJets_tau3ByTau1"          ,"Fat jet #tau2/#tau1"                          ,20     ,0.       ,1.      ) ; 
  AddHisto(cutname ,"_FatJets_CombinedSVBJetTags"  ,"Fat jet CSV discriminator"                    ,20     ,0.       ,1.      ) ; 
  AddHisto(cutname ,"_FatJets_Mass"                ,"Fat jet mass [GeV]"                           ,200    ,0.       ,2000.   ) ; 
  AddHisto(cutname ,"_FatJets_MassPruned"          ,"Fat jet pruned mass [GeV]"                    ,200    ,0.       ,2000.   ) ; 
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

  AddHisto(cutname ,"_nBJets"                      ,"N(b jets)"                                    ,20     ,-0.5     ,19.5    ) ; 
  AddHisto(cutname ,"_BJet_Pt"                     ,"p_{T} (b jet)[GeV]"                           ,200    ,0.       ,2000.   ) ;
  AddHisto(cutname ,"_BJet_Eta"                    ,"#eta (b jet)"                                 ,50     ,-4.      ,4.      ) ;

  AddHisto(cutname ,"_HT"                          ,"H_{T}[GeV]"                                   ,400   ,0.       ,4000.   ) ;
  AddHisto(cutname ,"_HTAllAK5"                    ,"H_{T} (AK5 jets) [GeV]"                       ,400   ,0.       ,4000.   ) ;
  AddHisto(cutname ,"_Nbprimes"                    ,";N (b' candidates); Events"                   ,5     ,-0.5     ,4.5     ) ;
  AddHisto(cutname ,"_bprimePt"                    ,"b' p_{T} [GeV]"                               ,200   ,0.       ,2000.   ) ;
  AddHisto(cutname ,"_bprimeMass"                  ,"b' mass [GeV]"                                ,200   ,0.       ,2000.   ) ;

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

  FatJetSelector jetSelCA8(fatJetSelParams_) ; 
  FatJetSelector jetSelHiggs(higgsJetSelParams_) ; 
  JetSelector jetSelAK5(jetSelParams_) ; 
  pat::strbitset retjetidak5 = jetSelAK5.getBitTemplate() ; 
  pat::strbitset retjetidca8 = jetSelCA8.getBitTemplate() ; 

  //// Event variables 
  bool passHLT(false) ; 
  int nGoodVtxs(0) ;
  bool isdata(0);
  double evtwt(1); 
  double puwt(1); 
  if (doPDFTree_) { 
    pdfTree_ = fs->make<TTree>("pdftree", "") ; 
    pdfTree_ -> Branch ("McFlag"   , &EvtInfo.McFlag  , "EvtInfo.McFlag/I"   ) ; 
    pdfTree_ -> Branch ("PDFid1"   , &EvtInfo.PDFid1  , "EvtInfo.PDFid1/I"   ) ; 
    pdfTree_ -> Branch ("PDFid2"   , &EvtInfo.PDFid1  , "EvtInfo.PDFid2/I"   ) ; 
    pdfTree_ -> Branch ("PDFx1"    , &EvtInfo.PDFx1   , "EvtInfo.PDFx1/F"    ) ; 
    pdfTree_ -> Branch ("PDFx2"    , &EvtInfo.PDFx1   , "EvtInfo.PDFx2/F"    ) ; 
    pdfTree_ -> Branch ("PDFv1"    , &EvtInfo.PDFv1   , "EvtInfo.PDFv1/F"    ) ; 
    pdfTree_ -> Branch ("PDFv2"    , &EvtInfo.PDFv2   , "EvtInfo.PDFv2/F"    ) ; 
    pdfTree_ -> Branch ("qScale"   , &EvtInfo.qScale  , "EvtInfo.qScale/F"   ) ; 
    pdfTree_ -> Branch ("alphaQCD" , &EvtInfo.alphaQCD, "EvtInfo.alphaQCD/F" ) ; 
    pdfTree_ -> Branch ("alphaQED" , &EvtInfo.alphaQED, "EvtInfo.alphaQED/F" ) ; 
    pdfTree_ -> Branch ("evtwt"    , &evtwt           , "evtwt/D"            ) ; 
    pdfTree_ -> Branch ("puwt"     , &puwt            , "puwt/D"             ) ; 
  }

  edm::LogInfo("StartingAnalysisLoop") << "Starting analysis loop\n";

  for(int entry=0; entry<maxEvents_; entry++) {
    if((entry%reportEvery_) == 0) edm::LogInfo("Event") << entry << " of " << maxEvents_ ; 

    HT HTAllAK5, MyHT ; 
    std::vector<TLorentzVector>p4bprimes ; 

    chain_->GetEntry(entry);

    isdata   = EvtInfo.McFlag ? 0 : 1; 
    if ( !isdata ) evtwt = EvtInfo.Weight ; 
    if ( doPUReweighting_ && !isdata ) puwt = LumiWeights_.weight(EvtInfo.TrueIT[0]) ; 

    if ( !isdata ) { //// Gen info 
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
        if(higgsDau==10) evtwt *= HiggsBRscaleFactors::higgsBBSf;
        else if(higgsDau==30) evtwt *= HiggsBRscaleFactors::higgsTauTauSf;
        else if(higgsDau==26) evtwt *= HiggsBRscaleFactors::higgsMuMuSf;
        else if(higgsDau==8)  evtwt *= HiggsBRscaleFactors::higgsCCSf;
        else if(higgsDau==6)  evtwt *= HiggsBRscaleFactors::higgsSSSf;
        else if(higgsDau==16) evtwt *= HiggsBRscaleFactors::higgsTTSf;
        else if(higgsDau==42) evtwt *= HiggsBRscaleFactors::higgsGGSf;
        else if(higgsDau==44) evtwt *= HiggsBRscaleFactors::higgsGammaGammaSf;
        else if(higgsDau==45) evtwt *= HiggsBRscaleFactors::higgsZGammaSf;
        else if(higgsDau==48) evtwt *= HiggsBRscaleFactors::higgsWWSf;
        else if(higgsDau==46) evtwt *= HiggsBRscaleFactors::higgsZZSf; 
        else evtwt *= 1 ; 
      } //// Higgs BR reweighting  
    } //// Gen info

    if ( TString(inputTTree_).Contains("ntuple") ) h_cutflow -> Fill("AllEvents", 1) ; 

    TriggerSelector trigSel(hltPaths_) ; 
    passHLT = trigSel.getTrigDecision(EvtInfo) ; 
    if ( !doTrigEff_ && !passHLT ) evtwt = 0 ; 
    else if ( doTrigEff_ || (!doTrigEff_  && passHLT) ) {
      h_cutflow -> Fill("TriggerSel", 1) ; 
      FillHisto(TString("TriggerSel")+TString("_nPVtx_NoPUWt"), nGoodVtxs, evtwt) ; 
      FillHisto(TString("TriggerSel")+TString("_nPVtx_PUWt"), nGoodVtxs, evtwt*puwt) ; 

      VertexSelector vtxSel(VtxInfo) ; 
      nGoodVtxs = vtxSel.NGoodVtxs(); 
      if (nGoodVtxs < 1)  { edm::LogInfo("NoGoodPrimaryVertex") << " No good primary vertex " ; continue ; }
      FillHisto(TString("VertexSel")+TString("_nPVtx_NoPUWt"), nGoodVtxs, evtwt) ; 
      FillHisto(TString("VertexSel")+TString("_nPVtx_PUWt"), nGoodVtxs, evtwt*puwt) ; 

      evtwt *= puwt ; 

      JetCollection fatjets ; 
      for (int ifatjet = 0; ifatjet < FatJetInfo.Size; ++ifatjet) { 
        if (jetSelCA8(FatJetInfo, ifatjet, SubJetInfo, retjetidca8) == 0) continue ; //// pt > 150 GeV, |eta| < 2.4 tau2/tau1 < 0.5 
        Jet thisjet(FatJetInfo, ifatjet) ; 
        fatjets.push_back(thisjet) ; 
      } //// Loop over fat jets 

      //// Apply JEC and b-tagging SFs for CA8 jets in MC  
      if ( !isdata && applyJEC_ ) {
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

      FillHisto(TString("VertexSel")+TString("_nFatJets"), fatjets.size(), evtwt) ; 
      JetCollection selectedFatJets, HiggsJets, AllHiggsJets, AllHiggsAntiHiggsJets; 
      JetCollection allAK5Jets, cleanedAK5Jets, allBJets, selectedBJets ; 
      for (JetCollection::const_iterator ifat = fatjets.begin(); ifat != fatjets.end(); ++ifat) {
        Jet thisjet(*ifat) ; 
        if ( thisjet.Pt() < fatJetPtMin_ || thisjet.Pt() > fatJetPtMax_ ) continue ; //// Apply fat jet pT cut  
        Jet subjet1(SubJetInfo, ifat->Jet_SubJet1Idx()) ;
        Jet subjet2(SubJetInfo, ifat->Jet_SubJet2Idx()) ;
        double subjet_dyphi = subjet1.DeltaR(subjet2) ; 
        if ( ifat->MassPruned() > fatJetPrunedMassMin_ && ifat->MassPruned() < fatJetPrunedMassMax_ 
            && subjet_dyphi > dRSubjetsMin_ && subjet_dyphi < dRSubjetsMax_ ) { //// Selecting fat jets with mass and dy cuts   
          selectedFatJets.push_back(thisjet) ; 
        } //// Selecting fat jets with mass and dy cuts 
        if (subjet1.CombinedSVBJetTags() < 0.244 || subjet2.CombinedSVBJetTags() < 0.244) continue ;  
        if (subjet1.CombinedSVBJetTags() <= subj1CSVDiscMin_ && subjet2.CombinedSVBJetTags() <= subj2CSVDiscMin_) { //// subjet disc
          AllHiggsAntiHiggsJets.push_back(thisjet) ; 
        }
        else if (subjet1.CombinedSVBJetTags() > subj1CSVDiscMin_ && subjet2.CombinedSVBJetTags() > subj2CSVDiscMin_) {
          AllHiggsAntiHiggsJets.push_back(thisjet) ; 
          AllHiggsJets.push_back(thisjet);
          if (thisjet.MassPruned() > fatJetPrunedMassMin_ && thisjet.MassPruned() < fatJetPrunedMassMax_ 
              && subjet_dyphi >= dRSubjetsMin_ && subjet_dyphi <= dRSubjetsMax_  ) { //// fat jet pruned mass 
            if ( !isdata && applyBTagSF_ ) { //// Apply Higgs-tagging scale factor 
              ApplyHiggsTagSF* higgsTagSF = new ApplyHiggsTagSF(double(subjet1.Pt()), double(subjet2.Pt()), 
                  double(subjet1.Eta()), double(subjet2.Eta()),
                  subjet1.GenFlavor(), subjet2.GenFlavor(), 
                  subjet1.CombinedSVBJetTags(), subjet2.CombinedSVBJetTags()) ; 
              evtwt *= higgsTagSF->GetHiggsTagSF() ;
              delete higgsTagSF ; 
            } //// Apply Higgs-tagging scale factor  
            HiggsJets.push_back(thisjet);
          } //// fat jet pruned mass   
        } //// Subjet CSV 

        FillHisto(TString("VertexSel")+TString("_FatJets_Pt")                 ,ifat->Pt() ,evtwt)  ;  
        FillHisto(TString("VertexSel")+TString("_FatJets_Eta")                ,ifat->Eta() ,evtwt)  ; 
        FillHisto(TString("VertexSel")+TString("_FatJets_CombinedSVBJetTags") ,ifat->CombinedSVBJetTags() ,evtwt)  ; 
        FillHisto(TString("VertexSel")+TString("_FatJets_tau2ByTau1")         ,ifat->tau2()/ifat->tau1() ,evtwt)  ; 
        FillHisto(TString("VertexSel")+TString("_FatJets_tau3ByTau2")         ,ifat->tau3()/ifat->tau2() ,evtwt)  ; 
        FillHisto(TString("VertexSel")+TString("_FatJets_tau3ByTau1")         ,ifat->tau3()/ifat->tau1() ,evtwt)  ; 
        FillHisto(TString("VertexSel")+TString("_FatJets_Mass")               ,ifat->Mass() ,evtwt)  ; 
        FillHisto(TString("VertexSel")+TString("_FatJets_MassPruned")         ,ifat->MassPruned() ,evtwt)  ; 
        FillHisto(TString("VertexSel")+TString("_FatJets_DYPhiSubjets")       ,subjet_dyphi ,evtwt)  ; 

        FillHisto(TString("VertexSel")+TString("_SubJet1_Pt") ,subjet1.Pt() ,evtwt)  ;  
        FillHisto(TString("VertexSel")+TString("_SubJet1_Eta") ,subjet1.Eta() ,evtwt)  ; 
        FillHisto(TString("VertexSel")+TString("_SubJet1_Mass") ,subjet1.Mass() ,evtwt)  ; 
        FillHisto(TString("VertexSel")+TString("_SubJet1_CombinedSVBJetTags") ,subjet1.CombinedSVBJetTags() ,evtwt)  ; 

        FillHisto(TString("VertexSel")+TString("_SubJet2_Pt") ,subjet2.Pt() ,evtwt)  ;  
        FillHisto(TString("VertexSel")+TString("_SubJet2_Eta") ,subjet2.Eta() ,evtwt)  ; 
        FillHisto(TString("VertexSel")+TString("_SubJet2_Mass") ,subjet2.Mass() ,evtwt)  ; 
        FillHisto(TString("VertexSel")+TString("_SubJet2_CombinedSVBJetTags") ,subjet2.CombinedSVBJetTags() ,evtwt)  ; 

      } //// Looping over all fat jets with pt > 150 GeV, |eta| < 2.4, loose jet ID, tau2/tau1 < 0.5 


      fatjets.clear() ; 

      for (int ijet = 0; ijet < JetInfo.Size; ++ijet) {
        retjetidak5.set(false) ;
        if (jetSelAK5(JetInfo, ijet,retjetidak5) == 0) continue ; 
        Jet thisjet(JetInfo, ijet) ;
        allAK5Jets.push_back(thisjet) ; 
      }

      if ( !isdata && applyJEC_ ) { //// Apply JEC for MC   

        //// All AK5 jets 
        JMEUncertUtil* jmeUtil_jer = new JMEUncertUtil(jmeParams_, allAK5Jets, "JERAK5MC", jerShift_) ; 
        JetCollection allAK5JetsJER = jmeUtil_jer->GetModifiedJetColl() ; 
        delete jmeUtil_jer ; 
        allAK5Jets.clear() ; 

        if ( abs(jesShift_) > 1E-6 ) {
          JMEUncertUtil* jmeUtil_jes = new JMEUncertUtil(jmeParams_, allAK5JetsJER, "JESAK5MC", jesShift_) ; 
          JetCollection allAK5JetsJES = jmeUtil_jes->GetModifiedJetColl() ; 
          delete jmeUtil_jes ; 

          for (JetCollection::const_iterator ijet = allAK5JetsJES.begin(); ijet != allAK5JetsJES.end(); ++ijet) {
            Jet thisjet(*ijet) ; 
            allAK5Jets.push_back(thisjet) ; 
          }
        }
        else {
          for (JetCollection::const_iterator ijet = allAK5JetsJER.begin(); ijet != allAK5JetsJER.end(); ++ijet) {
            Jet thisjet(*ijet) ; 
            allAK5Jets.push_back(thisjet) ; 
          }
        }

      } //// Apply JEC for MC 

      isolateCollection (AllHiggsAntiHiggsJets, allAK5Jets, cleanedAK5Jets) ; 

      if ( !isdata && applyBTagSF_  ) { //// Apply b-tagging SF for MC  
        ApplyBTagSF * btagsf =  new ApplyBTagSF(cleanedAK5Jets, 0.679, SFbShift_, SFlShift_) ;  
        allBJets.clear() ; 
        allBJets =  btagsf->getBtaggedJetsWithSF () ; 
        delete btagsf ; 
      } //// Apply b-tagging SF for MC 
      else {
        for (JetCollection::const_iterator ijet = cleanedAK5Jets.begin(); ijet != cleanedAK5Jets.end(); ++ijet) {
          if (ijet->CombinedSVBJetTags() > bJetCSVDiscMin_) allBJets.push_back(*ijet) ; 
        }
      }

      for (JetCollection::const_iterator ijet = allBJets.begin(); ijet != allBJets.end(); ++ijet) {
        if (ijet->Pt() > bJetPtMin_ /*&& ijet->CombinedSVBJetTags() > bJetCSVDiscMin_*/ ) selectedBJets.push_back(*ijet) ; 
      }

      cleanedAK5Jets.clear() ; 
      allBJets.clear() ; 

      HTAllAK5.setJetCollection(allAK5Jets) ; 
      HTAllAK5.buildHT() ; 

      MyHT.setJetCollection(HiggsJets) ; 
      MyHT.setJetCollection(selectedBJets) ; 
      MyHT.buildHT() ; 

      HTSelector htsel(HTSelParams_) ; 
      pat::strbitset retht = htsel.getBitTemplate() ; 
      retht.set(false) ; 

      if ( HiggsJets.size() >= 1
          && selectedBJets.size() >= 1
          && HTAllAK5.getHT() >= HTAK5Min_ ) {
        h_cutflow -> Fill("HTSel", 1) ; 
        FillHisto(TString("HTSel")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt) ; 
      }
      else evtwt = 0. ; 
    } /// Event passes HLT 
    if (doPDFTree_) pdfTree_->Fill() ; 

    /*
       if ( !doTrigEff_ ) {
       FillHisto(TString("TriggerSel")+TString("_nFatJets"), selectedFatJets.size(), evtwt) ; 
       FillHisto(TString("TriggerSel")+TString("_nAllAK5"), allAK5Jets.size() , evtwt) ; 
       FillHisto(TString("TriggerSel")+TString("_nBJets"), selectedBJets.size(), evtwt) ; 
       FillHisto(TString("TriggerSel")+TString("_nHJets"), HiggsJets.size(), evtwt) ; 
       FillHisto(TString("TriggerSel")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt) ; 
       FillHisto(TString("TriggerSel")+TString("_HT"), MyHT.getHT(), evtwt) ; 
       for (JetCollection::const_iterator ib = selectedBJets.begin(); ib != selectedBJets.end(); ++ib) { 
       FillHisto(TString("TriggerSel")+TString("_BJet_Pt"), ib->Pt(), evtwt) ; 
       FillHisto(TString("TriggerSel")+TString("_BJet_Eta"), ib->Eta(), evtwt) ; 
       }

       }

    //// Event selection  
    passBSel = false ; 
    if (selectedFatJets.size() >= 1) {
    h_cutflow -> Fill("FatJetSel", 1) ; 
    FillHisto(TString("FatJetSel")+TString("_nAllAK5"), allAK5Jets.size() , evtwt*puwt) ; 
    FillHisto(TString("FatJetSel")+TString("_nBJets"), selectedBJets.size(), evtwt) ; 
    FillHisto(TString("FatJetSel")+TString("_nHJets"), HiggsJets.size(), evtwt) ; 
    FillHisto(TString("FatJetSel")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt) ; 
    FillHisto(TString("FatJetSel")+TString("_HT"), MyHT.getHT(), evtwt) ; 
    for (JetCollection::const_iterator ib = selectedBJets.begin(); ib != selectedBJets.end(); ++ib) { 
    FillHisto(TString("FatJetSel")+TString("_BJet_Pt"), ib->Pt(), evtwt) ; 
    FillHisto(TString("FatJetSel")+TString("_BJet_Eta"), ib->Eta(), evtwt) ; 
    }
    for (JetCollection::const_iterator ijet = selectedFatJets.begin(); ijet != selectedFatJets.end(); ++ijet ) {
    int iSubJet1(ijet->Jet_SubJet1Idx()), iSubJet2(ijet->Jet_SubJet2Idx()); 
    TLorentzVector fatjet_p4, subj1_p4, subj2_p4;
    fatjet_p4.SetPtEtaPhiM(ijet->Pt(), ijet->Eta(), ijet->Phi (), ijet->Mass ()); 
    subj1_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet1], SubJetInfo.Eta[iSubJet1], SubJetInfo.Phi[iSubJet1], SubJetInfo.Mass[iSubJet1]);
    subj2_p4.SetPtEtaPhiM(SubJetInfo.Pt[iSubJet2], SubJetInfo.Eta[iSubJet2], SubJetInfo.Phi[iSubJet2], SubJetInfo.Mass[iSubJet2]); 
    double subjet_dy = subj1_p4.Rapidity() - subj2_p4.Rapidity() ;
    double subjet_dphi = subj1_p4.DeltaPhi(subj2_p4); ;
    double subjet_dyphi = sqrt( subjet_dy*subjet_dy + subjet_dphi*subjet_dphi ) ;

    FillHisto(TString("FatJetSel")+TString("_FatJets_Pt")                 ,fatjet_p4.Pt() ,evtwt)  ;  
    FillHisto(TString("FatJetSel")+TString("_FatJets_Eta")                ,fatjet_p4.Eta() ,evtwt)  ; 
    FillHisto(TString("FatJetSel")+TString("_FatJets_CombinedSVBJetTags") ,ijet->CombinedSVBJetTags() ,evtwt)  ; 
    FillHisto(TString("FatJetSel")+TString("_FatJets_tau2ByTau1")         ,ijet->tau2()/ijet->tau1() ,evtwt)  ; 
    FillHisto(TString("FatJetSel")+TString("_FatJets_tau3ByTau2")         ,ijet->tau3()/ijet->tau2() ,evtwt)  ; 
    FillHisto(TString("FatJetSel")+TString("_FatJets_tau3ByTau1")         ,ijet->tau3()/ijet->tau1() ,evtwt)  ; 
    FillHisto(TString("FatJetSel")+TString("_FatJets_Mass")               ,fatjet_p4.Mag() ,evtwt)  ; 
    FillHisto(TString("FatJetSel")+TString("_FatJets_MassPruned")         ,ijet->MassPruned() ,evtwt)  ; 
    FillHisto(TString("FatJetSel")+TString("_FatJets_DYPhiSubjets")       ,subjet_dyphi ,evtwt)  ; 

    FillHisto(TString("FatJetSel")+TString("_SubJet1_Pt") ,subj1_p4.Pt() ,evtwt)  ;  
    FillHisto(TString("FatJetSel")+TString("_SubJet1_Eta") ,subj1_p4.Eta() ,evtwt)  ; 
    FillHisto(TString("FatJetSel")+TString("_SubJet1_Mass") ,subj1_p4.Mag() ,evtwt)  ; 
    FillHisto(TString("FatJetSel")+TString("_SubJet1_CombinedSVBJetTags") ,SubJetInfo.CombinedSVBJetTags[iSubJet1] ,evtwt)  ; 

    FillHisto(TString("FatJetSel")+TString("_SubJet2_Pt") ,subj2_p4.Pt() ,evtwt)  ;  
    FillHisto(TString("FatJetSel")+TString("_SubJet2_Eta") ,subj2_p4.Eta() ,evtwt)  ; 
    FillHisto(TString("FatJetSel")+TString("_SubJet2_Mass") ,subj2_p4.Mag() ,evtwt)  ; 
    FillHisto(TString("FatJetSel")+TString("_SubJet2_CombinedSVBJetTags") ,SubJetInfo.CombinedSVBJetTags[iSubJet2] ,evtwt)  ; 

    }
    for (JetCollection::const_iterator ih = HiggsJets.begin(); ih != HiggsJets.end(); ++ih) { 
      FillHisto(TString("FatJetSel")+TString("_HiggsJet_Pt"), ih->Pt(), evtwt) ; 
      FillHisto(TString("FatJetSel")+TString("_HiggsJet_Eta"), ih->Eta(), evtwt) ; 
      FillHisto(TString("FatJetSel")+TString("_HiggsJet_Mass") ,ih->Mass() ,evtwt)  ;
    }

    if (HiggsJets.size() >= 1) {
      h_cutflow -> Fill("HiggsJetSel", 1) ; 
      FillHisto(TString("HiggsJetSel")+TString("_nAllAK5"), allAK5Jets.size() , evtwt) ; 
      FillHisto(TString("HiggsJetSel")+TString("_nBJets"), selectedBJets.size(), evtwt) ; 
      FillHisto(TString("HiggsJetSel")+TString("_nHJets"), HiggsJets.size(), evtwt) ; 
      FillHisto(TString("HiggsJetSel")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt) ; 
      FillHisto(TString("HiggsJetSel")+TString("_HT"), MyHT.getHT(), evtwt) ; 
      for (JetCollection::const_iterator ib = selectedBJets.begin(); ib != selectedBJets.end(); ++ib) { 
        FillHisto(TString("HiggsJetSel")+TString("_BJet_Pt"), ib->Pt(), evtwt) ; 
        FillHisto(TString("HiggsJetSel")+TString("_BJet_Eta"), ib->Eta(), evtwt) ; 
      }
      for (JetCollection::const_iterator ih = HiggsJets.begin(); ih != HiggsJets.end(); ++ih) { 
        FillHisto(TString("HiggsJetSel")+TString("_HiggsJet_Pt"), ih->Pt(), evtwt) ; 
        FillHisto(TString("HiggsJetSel")+TString("_HiggsJet_Eta"), ih->Eta(), evtwt) ; 
        FillHisto(TString("HiggsJetSel")+TString("_HiggsJet_Mass") ,ih->Mass() ,evtwt)  ;
      }

      if (allAK5Jets.size() > 0) { 
        h_cutflow -> Fill("JetSel", 1) ; 
        FillHisto(TString("JetSel")+TString("_nAllAK5"), allAK5Jets.size() , evtwt) ; 
        FillHisto(TString("JetSel")+TString("_nBJets"), selectedBJets.size(), evtwt) ; 
        FillHisto(TString("JetSel")+TString("_nHJets"), HiggsJets.size(), evtwt) ; 
        FillHisto(TString("JetSel")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt) ; 
        FillHisto(TString("JetSel")+TString("_HT"), MyHT.getHT(), evtwt) ; 
        for (JetCollection::const_iterator ib = selectedBJets.begin(); ib != selectedBJets.end(); ++ib) { 
          FillHisto(TString("JetSel")+TString("_BJet_Pt"), ib->Pt(), evtwt) ; 
          FillHisto(TString("JetSel")+TString("_BJet_Eta"), ib->Eta(), evtwt) ; 
        }
        for (JetCollection::const_iterator ih = HiggsJets.begin(); ih != HiggsJets.end(); ++ih) { 
          FillHisto(TString("JetSel")+TString("_HiggsJet_Pt"), ih->Pt(), evtwt) ; 
          FillHisto(TString("JetSel")+TString("_HiggsJet_Eta"), ih->Eta(), evtwt) ; 
          FillHisto(TString("JetSel")+TString("_HiggsJet_Mass") ,ih->Mass() ,evtwt)  ;
        }

        if (selectedBJets.size() >= 1 ) { 
          h_cutflow -> Fill("BJetsSel", 1) ;
          FillHisto(TString("BJetsSel")+TString("_nAllAK5"), allAK5Jets.size() , evtwt*puwt) ; 
          FillHisto(TString("BJetsSel")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt) ; 
          FillHisto(TString("BJetsSel")+TString("_HT"), MyHT.getHT(), evtwt) ; 
          for (JetCollection::const_iterator ib = selectedBJets.begin(); ib != selectedBJets.end(); ++ib) { 
            FillHisto(TString("BJetsSel")+TString("_BJet_Pt"), ib->Pt(), evtwt) ; 
            FillHisto(TString("BJetsSel")+TString("_BJet_Eta"), ib->Eta(), evtwt) ; 
          }
          for (JetCollection::const_iterator ih = HiggsJets.begin(); ih != HiggsJets.end(); ++ih) { 
            FillHisto(TString("BJetsSel")+TString("_HiggsJet_Pt"), ih->Pt(), evtwt) ; 
            FillHisto(TString("BJetsSel")+TString("_HiggsJet_Eta"), ih->Eta(), evtwt) ; 
            FillHisto(TString("BJetsSel")+TString("_HiggsJet_Mass") ,ih->Mass() ,evtwt)  ;
          }

          if ( htsel(HTAllAK5, retht) != 0 ) { //// If event passes HT 
            h_cutflow -> Fill("HTSel", 1) ; 
            FillHisto(TString("HTSel")+TString("_nAllAK5"), allAK5Jets.size() , evtwt*puwt) ; 
            FillHisto(TString("HTSel")+TString("_nBJets"), selectedBJets.size(), evtwt) ; 
            FillHisto(TString("HTSel")+TString("_nHJets"), HiggsJets.size(), evtwt) ; 
            FillHisto(TString("HTSel")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt) ; 
            FillHisto(TString("HTSel")+TString("_HT"), MyHT.getHT(), evtwt) ; 
            for (JetCollection::const_iterator ib = selectedBJets.begin(); ib != selectedBJets.end(); ++ib) { 
              FillHisto(TString("HTSel")+TString("_BJet_Pt"), ib->Pt(), evtwt) ; 
              FillHisto(TString("HTSel")+TString("_BJet_Eta"), ib->Eta(), evtwt) ; 
            }
            for (JetCollection::const_iterator ih = HiggsJets.begin(); ih != HiggsJets.end(); ++ih) { 
              FillHisto(TString("HTSel")+TString("_HiggsJet_Pt"), ih->Pt(), evtwt) ; 
              FillHisto(TString("HTSel")+TString("_HiggsJet_Eta"), ih->Eta(), evtwt) ; 
              FillHisto(TString("HTSel")+TString("_HiggsJet_Mass") ,ih->Mass() ,evtwt)  ;
            }

            //// Reconstruct b' candidates
            for (JetCollection::const_iterator ihig = HiggsJets.begin(); ihig != HiggsJets.end(); ++ihig) { 
              unsigned int closestBIndex(selectedBJets.size()) ;
              double deltaR(TMath::Pi()) ; 
              for (JetCollection::const_iterator ib = selectedBJets.begin(); ib != selectedBJets.end(); ++ib) { 
                if ( ihig->DeltaR(*ib) < deltaR) {
                  deltaR = ihig->DeltaR(*ib) ; 
                  closestBIndex = ib - selectedBJets.begin() ; 
                }
              }
              if (closestBIndex < selectedBJets.size()) {
                p4bprimes.push_back(ihig->p4() + (selectedBJets.at(closestBIndex)).p4()) ; 
              }
            } //// Reconstruct b' candidates
            FillHisto(TString("HTSel")+TString("_Nbprimes"), p4bprimes.size(), evtwt) ; 
            for (std::vector<TLorentzVector>::const_iterator ibprime = p4bprimes.begin(); ibprime != p4bprimes.end(); ++ibprime) {
              FillHisto(TString("HTSel")+TString("_bprimePt") ,ibprime->Pt() ,evtwt)  ; 
              FillHisto(TString("HTSel")+TString("_bprimeMass") ,ibprime->Mag() ,evtwt)  ;
            } //// Fill all bprime candidates 

            if (HiggsJets.size() == 1 && selectedBJets.size() == 1) { //// 1 HJet and 1 bjet 
              h_cutflow -> Fill("HTSel_1H_1b", 1) ; 
              FillHisto(TString("HTSel_1H_1b")+TString("_nAllAK5"), allAK5Jets.size() , evtwt*puwt) ; 
              FillHisto(TString("HTSel_1H_1b")+TString("_nBJets"), selectedBJets.size(), evtwt) ; 
              FillHisto(TString("HTSel_1H_1b")+TString("_nHJets"), HiggsJets.size(), evtwt) ; 
              FillHisto(TString("HTSel_1H_1b")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt) ; 
              FillHisto(TString("HTSel_1H_1b")+TString("_HT"), MyHT.getHT(), evtwt) ; 
              for (JetCollection::const_iterator ib = selectedBJets.begin(); ib != selectedBJets.end(); ++ib) { 
                FillHisto(TString("HTSel_1H_1b")+TString("_BJet_Pt"), ib->Pt(), evtwt) ; 
                FillHisto(TString("HTSel_1H_1b")+TString("_BJet_Eta"), ib->Eta(), evtwt) ; 
              }
              for (JetCollection::const_iterator ih = HiggsJets.begin(); ih != HiggsJets.end(); ++ih) { 
                FillHisto(TString("HTSel_1H_1b")+TString("_HiggsJet_Pt"), ih->Pt(), evtwt) ; 
                FillHisto(TString("HTSel_1H_1b")+TString("_HiggsJet_Eta"), ih->Eta(), evtwt) ; 
                FillHisto(TString("HTSel_1H_1b")+TString("_HiggsJet_Mass") ,ih->Mass() ,evtwt)  ;
              }
              FillHisto(TString("HTSel_1H_1b")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt) ; 
              FillHisto(TString("HTSel_1H_1b")+TString("_HT"), MyHT.getHT(), evtwt) ; 
              FillHisto(TString("HTSel_1H_1b")+TString("_Nbprimes"), p4bprimes.size(), evtwt) ; 
              for (std::vector<TLorentzVector>::const_iterator ibprime = p4bprimes.begin(); ibprime != p4bprimes.end(); ++ibprime) {
                FillHisto(TString("HTSel_1H_1b")+TString("_bprimePt") ,ibprime->Pt() ,evtwt)  ; 
                FillHisto(TString("HTSel_1H_1b")+TString("_bprimeMass") ,ibprime->Mag() ,evtwt)  ;
              } //// Fill all bprime candidates 
            } //// 1 HJet and 1 bjet 

            if (HiggsJets.size() == 1 && selectedBJets.size() >= 2) { //// 1 HJet and >=2 bjet 
              h_cutflow -> Fill("HTSel_1H_2b", 1) ; 
              FillHisto(TString("HTSel_1H_2b")+TString("_nAllAK5"), allAK5Jets.size() , evtwt*puwt) ; 
              FillHisto(TString("HTSel_1H_2b")+TString("_nBJets"), selectedBJets.size(), evtwt) ; 
              FillHisto(TString("HTSel_1H_2b")+TString("_nHJets"), HiggsJets.size(), evtwt) ; 
              FillHisto(TString("HTSel_1H_2b")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt) ; 
              FillHisto(TString("HTSel_1H_2b")+TString("_HT"), MyHT.getHT(), evtwt) ; 
              for (JetCollection::const_iterator ib = selectedBJets.begin(); ib != selectedBJets.end(); ++ib) { 
                FillHisto(TString("HTSel_1H_2b")+TString("_BJet_Pt"), ib->Pt(), evtwt) ; 
                FillHisto(TString("HTSel_1H_2b")+TString("_BJet_Eta"), ib->Eta(), evtwt) ; 
              }
              for (JetCollection::const_iterator ih = HiggsJets.begin(); ih != HiggsJets.end(); ++ih) { 
                FillHisto(TString("HTSel_1H_2b")+TString("_HiggsJet_Pt"), ih->Pt(), evtwt) ; 
                FillHisto(TString("HTSel_1H_2b")+TString("_HiggsJet_Eta"), ih->Eta(), evtwt) ; 
                FillHisto(TString("HTSel_1H_2b")+TString("_HiggsJet_Mass") ,ih->Mass() ,evtwt)  ;
              }
              FillHisto(TString("HTSel_1H_2b")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt) ; 
              FillHisto(TString("HTSel_1H_2b")+TString("_HT"), MyHT.getHT(), evtwt) ; 
              FillHisto(TString("HTSel_1H_2b")+TString("_Nbprimes"), p4bprimes.size(), evtwt) ; 
              for (std::vector<TLorentzVector>::const_iterator ibprime = p4bprimes.begin(); ibprime != p4bprimes.end(); ++ibprime) {
                FillHisto(TString("HTSel_1H_2b")+TString("_bprimePt") ,ibprime->Pt() ,evtwt)  ; 
                FillHisto(TString("HTSel_1H_2b")+TString("_bprimeMass") ,ibprime->Mag() ,evtwt)  ;
              } //// Fill all bprime candidates 
            } //// 1 HJet and >=2 bjet 

            if (HiggsJets.size() >= 2 && selectedBJets.size() == 1) { //// >= 2 HJet and 1 bjet 
              h_cutflow -> Fill("HTSel_2H_1b", 1) ; 
              FillHisto(TString("HTSel_2H_1b")+TString("_nAllAK5"), allAK5Jets.size() , evtwt*puwt) ; 
              FillHisto(TString("HTSel_2H_1b")+TString("_nBJets"), selectedBJets.size(), evtwt) ; 
              FillHisto(TString("HTSel_2H_1b")+TString("_nHJets"), HiggsJets.size(), evtwt) ; 
              FillHisto(TString("HTSel_2H_1b")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt) ; 
              FillHisto(TString("HTSel_2H_1b")+TString("_HT"), MyHT.getHT(), evtwt) ; 
              for (JetCollection::const_iterator ib = selectedBJets.begin(); ib != selectedBJets.end(); ++ib) { 
                FillHisto(TString("HTSel_2H_1b")+TString("_BJet_Pt"), ib->Pt(), evtwt) ; 
                FillHisto(TString("HTSel_2H_1b")+TString("_BJet_Eta"), ib->Eta(), evtwt) ; 
              }
              for (JetCollection::const_iterator ih = HiggsJets.begin(); ih != HiggsJets.end(); ++ih) { 
                FillHisto(TString("HTSel_2H_1b")+TString("_HiggsJet_Pt"), ih->Pt(), evtwt) ; 
                FillHisto(TString("HTSel_2H_1b")+TString("_HiggsJet_Eta"), ih->Eta(), evtwt) ; 
                FillHisto(TString("HTSel_2H_1b")+TString("_HiggsJet_Mass") ,ih->Mass() ,evtwt)  ;
              }
              FillHisto(TString("HTSel_2H_1b")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt) ; 
              FillHisto(TString("HTSel_2H_1b")+TString("_HT"), MyHT.getHT(), evtwt) ; 
              FillHisto(TString("HTSel_2H_1b")+TString("_Nbprimes"), p4bprimes.size(), evtwt) ; 
              for (std::vector<TLorentzVector>::const_iterator ibprime = p4bprimes.begin(); ibprime != p4bprimes.end(); ++ibprime) {
                FillHisto(TString("HTSel_2H_1b")+TString("_bprimePt") ,ibprime->Pt() ,evtwt)  ; 
                FillHisto(TString("HTSel_2H_1b")+TString("_bprimeMass") ,ibprime->Mag() ,evtwt)  ;
              } //// Fill all bprime candidates 
            } //// >= 2 HJet and 1 bjet 

            if (HiggsJets.size() >= 2 && selectedBJets.size() >= 2) { //// >= 2 HJet and >=2 bjet  
              h_cutflow -> Fill("HTSel_2H_2b", 1) ; 
              FillHisto(TString("HTSel_2H_2b")+TString("_nAllAK5"), allAK5Jets.size() , evtwt*puwt) ; 
              FillHisto(TString("HTSel_2H_2b")+TString("_nBJets"), selectedBJets.size(), evtwt) ; 
              FillHisto(TString("HTSel_2H_2b")+TString("_nHJets"), HiggsJets.size(), evtwt) ; 
              FillHisto(TString("HTSel_2H_2b")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt) ; 
              FillHisto(TString("HTSel_2H_2b")+TString("_HT"), MyHT.getHT(), evtwt) ; 
              for (JetCollection::const_iterator ib = selectedBJets.begin(); ib != selectedBJets.end(); ++ib) { 
                FillHisto(TString("HTSel_2H_2b")+TString("_BJet_Pt"), ib->Pt(), evtwt) ; 
                FillHisto(TString("HTSel_2H_2b")+TString("_BJet_Eta"), ib->Eta(), evtwt) ; 
              }
              for (JetCollection::const_iterator ih = HiggsJets.begin(); ih != HiggsJets.end(); ++ih) { 
                FillHisto(TString("HTSel_2H_2b")+TString("_HiggsJet_Pt"), ih->Pt(), evtwt) ; 
                FillHisto(TString("HTSel_2H_2b")+TString("_HiggsJet_Eta"), ih->Eta(), evtwt) ; 
                FillHisto(TString("HTSel_2H_2b")+TString("_HiggsJet_Mass") ,ih->Mass() ,evtwt)  ;
              }
              FillHisto(TString("HTSel_2H_2b")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt) ; 
              FillHisto(TString("HTSel_2H_2b")+TString("_HT"), MyHT.getHT(), evtwt) ; 
              FillHisto(TString("HTSel_2H_2b")+TString("_Nbprimes"), p4bprimes.size(), evtwt) ; 
              for (std::vector<TLorentzVector>::const_iterator ibprime = p4bprimes.begin(); ibprime != p4bprimes.end(); ++ibprime) {
                FillHisto(TString("HTSel_2H_2b")+TString("_bprimePt") ,ibprime->Pt() ,evtwt)  ; 
                FillHisto(TString("HTSel_2H_2b")+TString("_bprimeMass") ,ibprime->Mag() ,evtwt)  ;
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
    FillHisto(TString("TriggerSel")+TString("_nPVtx_NoPUWt"), nGoodVtxs, evtwt) ; 
    FillHisto(TString("TriggerSel")+TString("_nPVtx_PUWt"), nGoodVtxs, evtwt*puwt) ; 
    FillHisto(TString("TriggerSel")+TString("_nFatJets"), selectedFatJets.size(), evtwt) ; 
    FillHisto(TString("TriggerSel")+TString("_nAllAK5"), allAK5Jets.size() , evtwt) ; 
    FillHisto(TString("TriggerSel")+TString("_nBJets"), selectedBJets.size(), evtwt) ; 
    FillHisto(TString("TriggerSel")+TString("_nHJets"), HiggsJets.size(), evtwt) ; 
    FillHisto(TString("TriggerSel")+TString("_HTAllAK5"), HTAllAK5.getHT(), evtwt) ; 
    FillHisto(TString("TriggerSel")+TString("_HT"), MyHT.getHT(), evtwt) ; 
    for (JetCollection::const_iterator ib = selectedBJets.begin(); ib != selectedBJets.end(); ++ib) { 
      FillHisto(TString("TriggerSel")+TString("_BJet_Pt"), ib->Pt(), evtwt) ; 
      FillHisto(TString("TriggerSel")+TString("_BJet_Eta"), ib->Eta(), evtwt) ; 
    }

  }
  */ 

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
