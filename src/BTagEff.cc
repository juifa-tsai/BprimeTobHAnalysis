// -*- C++ -*-
//
// Package:    BTagEff
// Class:      BTagEff
// 
/**\class BTagEff BTagEff.cc Bprime_kit/BTagEff/src/BTagEff.cc

Description: 
Analyzer class for Bprime -> b Higgs studies 
- National Taiwan University - 

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Devdatta Majumder 
//         Created:  
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
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH1I.h>
#include <TEfficiency.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "BpbH/BprimeTobH/interface/TriggerBooking.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h" 

#include "BpbH/BprimeTobHAnalysis/interface/EventSelector.h"

//
// class declaration
//

class BTagEff : public edm::EDAnalyzer{
  public:
    explicit BTagEff(const edm::ParameterSet&);
    ~BTagEff();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();

		double NearestJet( JetCollection* jets, const Jet& myJet, const Jet* nearestJet );
		void isolateCollection( JetCollection control, JetCollection input, JetCollection* output, double dR=1.2 );

    edm::LumiReWeighting LumiWeights_; 

    // ----------member data ---------------------------

    //// Configurables 

    int                             maxEvents_; 
    const int                       reportEvery_; 
    const std::string               inputTTree_;
    const std::vector<std::string>  inputFiles_;
    const edm::ParameterSet         hltPaths_; 
    const int                       doPUReweighting_;
    const std::string               file_PUDistMC_;
    const std::string               file_PUDistData_;
    const std::string               hist_PUDistMC_;
    const std::string               hist_PUDistData_;

    const edm::ParameterSet         fatJetSelParams_ ; 
    const edm::ParameterSet         fatbjetSelParams_ ; 
    const edm::ParameterSet         jetSelParams_ ; 
    const edm::ParameterSet         bjetSelParams_ ; 
    const edm::ParameterSet         evtSelParams_; 

    TChain*            chain_;

    GenInfoBranches    GenInfo;
    EvtInfoBranches    EvtInfo;
    VertexInfoBranches VtxInfo;
    JetInfoBranches    JetInfo;
    JetInfoBranches    FatJetInfo;
    JetInfoBranches    SubJetInfo;

    GenInfoBranches    newGenInfo;
    JetInfoBranches    newJetInfo;
    JetInfoBranches    newFatJetInfo;
    JetInfoBranches    newSubJetInfo;

    edm::Service<TFileService> fs; 

    EventSelector *eSelector_ ; 
    std::vector<std::string> cutLevels_; 

    bool isData_; 
    bool McFlag_; 
    double evtwt_; 
    double puweight_;
    double evtwtPu_; 

    TH1D*   h_cutflow ; 
    TH1D*   h_nAK5_bFlav ; 
    TH1D*   h_nAK5_cFlav ; 
    TH1D*   h_nAK5_lFlav ;  
    TH1D*   h_nAK5_0Flav ;  
    TH2D*   h2_AK5_pt_eta_bFlav ; 
    TH2D*   h2_AK5_pt_eta_cFlav ; 
    TH2D*   h2_AK5_pt_eta_lFlav ;  
    TH2D*   h2_AK5_pt_eta_0Flav ;  
    TH2D*   h2_AK5_pt_eta_bFlav_btagged ; 
    TH2D*   h2_AK5_pt_eta_cFlav_btagged ; 
    TH2D*   h2_AK5_pt_eta_lFlav_btagged ;  
    TH2D*   h2_AK5_pt_eta_0Flav_btagged ;  
    TH2D*   h2_AK5_btagEff_pt_eta_bFlav ; 
    TH2D*   h2_AK5_btagEff_pt_eta_cFlav ; 
    TH2D*   h2_AK5_btagEff_pt_eta_lFlav ;  
    TH2D*   h2_AK5_btagEff_pt_eta_0Flav ;  
    TH1D*   h_AK5_btagEff_pt_bFlav ; 
    TH1D*   h_AK5_btagEff_pt_cFlav ; 
    TH1D*   h_AK5_btagEff_pt_lFlav ;  
    TH1D*   h_AK5_btagEff_pt_0Flav ;  
    TH1D*   h_AK5_nbtagged_bFlav ; 
    TH1D*   h_AK5_nbtagged_cFlav ; 
    TH1D*   h_AK5_nbtagged_lFlav ;  
    TH1D*   h_AK5_nbtagged_0Flav ;  

    TH1D*   h_nCA8 ; 
    TH1D*   h_nCA8_bFlav ; 
    TH1D*   h_nCA8_cFlav ; 
    TH1D*   h_nCA8_lFlav ;  
    TH1D*   h_nCA8_0Flav ;  
    TH2D*   h2_CA8_pt_eta_bFlav ; 
    TH2D*   h2_CA8_pt_eta_cFlav ; 
    TH2D*   h2_CA8_pt_eta_lFlav ;  
    TH2D*   h2_CA8_pt_eta_0Flav ;  
    TH2D*   h2_CA8_pt_eta_bFlav_btagged ; 
    TH2D*   h2_CA8_pt_eta_cFlav_btagged ; 
    TH2D*   h2_CA8_pt_eta_lFlav_btagged ;  
    TH2D*   h2_CA8_pt_eta_0Flav_btagged ;  
    TH2D*   h2_CA8_btagEff_pt_eta_bFlav ; 
    TH2D*   h2_CA8_btagEff_pt_eta_cFlav ; 
    TH2D*   h2_CA8_btagEff_pt_eta_lFlav ;  
    TH2D*   h2_CA8_btagEff_pt_eta_0Flav ;  
    TH1D*   h_CA8_btagEff_pt_bFlav ; 
    TH1D*   h_CA8_btagEff_pt_cFlav ; 
    TH1D*   h_CA8_btagEff_pt_lFlav ;  
    TH1D*   h_CA8_btagEff_pt_0Flav ;  
    TH1D*   h_CA8_nbtagged_bFlav ; 
    TH1D*   h_CA8_nbtagged_cFlav ; 
    TH1D*   h_CA8_nbtagged_lFlav ;  
    TH1D*   h_CA8_nbtagged_0Flav ;  

    TH2D*   h2_pt_AK5_CA8_bFlav ; 
    TH2D*   h2_pt_AK5_CA8_cFlav ; 
    TH2D*   h2_pt_AK5_CA8_lFlav ; 
    TH2D*   h2_pt_AK5_CA8_0Flav ; 

    TH1D*   h_DR_AK5_CA8_bFlav ; 
    TH1D*   h_DR_AK5_CA8_cFlav ; 
    TH1D*   h_DR_AK5_CA8_lFlav ; 
    TH1D*   h_DR_AK5_CA8_0Flav ; 

};

//
// constructors and destructor
//
BTagEff::BTagEff(const edm::ParameterSet& iConfig) : 
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
  fatJetSelParams_(iConfig.getParameter<edm::ParameterSet>("FatJetSelParams")), 
  fatbjetSelParams_(iConfig.getParameter<edm::ParameterSet>("FatBJetSelParams")), 
  jetSelParams_(iConfig.getParameter<edm::ParameterSet>("JetSelParams")), 
  bjetSelParams_(iConfig.getParameter<edm::ParameterSet>("BJetSelParams")), 
  evtSelParams_(iConfig.getParameter<edm::ParameterSet>("EvtSelParams")),
  isData_(0),
  evtwt_(1), 
  puweight_(1)  
{ 

  if( doPUReweighting_) LumiWeights_ = edm::LumiReWeighting(file_PUDistMC_, file_PUDistData_, hist_PUDistMC_, hist_PUDistData_);

}


BTagEff::~BTagEff(){ 
  delete chain_;
}

// ------------ method called once each job just before starting event loop  ------------
void BTagEff::beginJob(){ 
  chain_  = new TChain(inputTTree_.c_str());

  for(unsigned i=0; i<inputFiles_.size(); ++i){
    chain_->Add(inputFiles_.at(i).c_str());
    TFile *f = TFile::Open(inputFiles_.at(i).c_str(),"READ");
    f->Close();
  }

  EvtInfo.Register(chain_);
  VtxInfo.Register(chain_);
  GenInfo.Register(chain_);
  JetInfo.Register(chain_,"JetInfo");
  FatJetInfo.Register(chain_,"FatJetInfo");
  SubJetInfo.Register(chain_,"SubJetInfo");

  if(  maxEvents_<0 || maxEvents_>chain_->GetEntries()) maxEvents_ = chain_->GetEntries();

  eSelector_ = new EventSelector(evtSelParams_,EvtInfo,VtxInfo,JetInfo,FatJetInfo,SubJetInfo);  
  cutLevels_ = eSelector_->getCutLevels();

  h_cutflow = fs->make<TH1D>("h_cutflow" ,"Cut flow" ,(int)cutLevels_.size() ,0. ,(int)cutLevels_.size() ); 
  for (int ii = 0; ii < (int)cutLevels_.size(); ++ii) {
    h_cutflow -> GetXaxis() -> SetBinLabel(ii+1, (cutLevels_.at(ii)).c_str()) ;  
  }

  h2_AK5_pt_eta_bFlav = fs->make<TH2D>("h2_AK5_pt_eta_bFlav" ,";p_{T} (b flavour jet) [GeV]; #eta (jet); Events"     ,1000 ,0. ,1000. ,80 ,-4. ,4.) ; 
  h2_AK5_pt_eta_cFlav = fs->make<TH2D>("h2_AK5_pt_eta_cFlav" ,";p_{T} (c flavour jet) [GeV]; #eta (jet); Events"     ,1000 ,0. ,1000. ,80 ,-4. ,4.) ; 
  h2_AK5_pt_eta_lFlav = fs->make<TH2D>("h2_AK5_pt_eta_lFlav" ,";p_{T} (light flavour jet) [GeV]; #eta (jet); Events" ,1000 ,0. ,1000. ,80 ,-4. ,4.) ;  
  h2_AK5_pt_eta_0Flav = fs->make<TH2D>("h2_AK5_pt_eta_0Flav" ,";p_{T} (flavourless jet) [GeV]; #eta (jet); Events"   ,1000 ,0. ,1000. ,80 ,-4. ,4.) ;  
  h_nAK5_bFlav = fs->make<TH1D>("h_nAK5_bFlav" , ";N(b-tagged jets) ;Events" ,20 ,-0.5 ,19.5) ; 
  h_nAK5_cFlav = fs->make<TH1D>("h_nAK5_cFlav" , ";N(b-tagged jets) ;Events" ,20 ,-0.5 ,19.5) ; 
  h_nAK5_lFlav = fs->make<TH1D>("h_nAK5_lFlav" , ";N(b-tagged jets) ;Events" ,20 ,-0.5 ,19.5) ;  
  h_nAK5_0Flav = fs->make<TH1D>("h_nAK5_0Flav" , ";N(b-tagged jets) ;Events" ,20 ,-0.5 ,19.5) ;  
  h2_AK5_pt_eta_bFlav_btagged = fs->make<TH2D>("h2_AK5_pt_eta_bFlav_btagged" , ";p_{T} (b flavour jet) [GeV]; #eta (jet); Events"     ,1000 ,0. ,1000. ,80 ,-4. ,4.) ; 
  h2_AK5_pt_eta_cFlav_btagged = fs->make<TH2D>("h2_AK5_pt_eta_cFlav_btagged" , ";p_{T} (c flavour jet) [GeV]; #eta (jet); Events"     ,1000 ,0. ,1000. ,80 ,-4. ,4.) ; 
  h2_AK5_pt_eta_lFlav_btagged = fs->make<TH2D>("h2_AK5_pt_eta_lFlav_btagged" , ";p_{T} (light flavour jet) [GeV]; #eta (jet); Events" ,1000 ,0. ,1000. ,80 ,-4. ,4.) ;  
  h2_AK5_pt_eta_0Flav_btagged = fs->make<TH2D>("h2_AK5_pt_eta_0Flav_btagged" , ";p_{T} (flavourless jet) [GeV]; #eta (jet); Events"   ,1000 ,0. ,1000. ,80 ,-4. ,4.) ;  
  h2_AK5_btagEff_pt_eta_bFlav = fs->make<TH2D>("h2_AK5_btagEff_pt_eta_bFlav" , ";p_{T} (b flavour jet) [GeV]; #eta (jet); CSVM b-tagging efficiency"     ,1000 ,0. ,1000. ,80 ,-4. ,4.) ; 
  h2_AK5_btagEff_pt_eta_cFlav = fs->make<TH2D>("h2_AK5_btagEff_pt_eta_cFlav" , ";p_{T} (c flavour jet) [GeV]; #eta (jet); CSVM b-tagging efficiency"     ,1000 ,0. ,1000. ,80 ,-4. ,4.) ; 
  h2_AK5_btagEff_pt_eta_lFlav = fs->make<TH2D>("h2_AK5_btagEff_pt_eta_lFlav" , ";p_{T} (light flavour jet) [GeV]; #eta (jet); CSVM b-tagging efficiency" ,1000 ,0. ,1000. ,80 ,-4. ,4.) ;  
  h2_AK5_btagEff_pt_eta_0Flav = fs->make<TH2D>("h2_AK5_btagEff_pt_eta_0Flav" , ";p_{T} (flvourless jet) [GeV]; #eta (jet); CSVM b-tagging efficiency"    ,1000 ,0. ,1000. ,80 ,-4. ,4.) ;  
  h_AK5_btagEff_pt_bFlav = fs->make<TH1D>("h_AK5_btagEff_pt_bFlav" , ";p_{T} (b flavour jet) [GeV]; CSVM b-tagging efficiency"     ,20 ,0. ,1000.) ; 
  h_AK5_btagEff_pt_cFlav = fs->make<TH1D>("h_AK5_btagEff_pt_cFlav" , ";p_{T} (c flavour jet) [GeV]; CSVM b-tagging efficiency"     ,20 ,0. ,1000.) ; 
  h_AK5_btagEff_pt_lFlav = fs->make<TH1D>("h_AK5_btagEff_pt_lFlav" , ";p_{T} (light flavour jet) [GeV]; CSVM b-tagging efficiency" ,20 ,0. ,1000.) ;  
  h_AK5_btagEff_pt_0Flav = fs->make<TH1D>("h_AK5_btagEff_pt_0Flav" , ";p_{T} (flvourless jet) [GeV]; CSVM b-tagging efficiency"    ,20 ,0. ,1000.) ;  
  h_AK5_nbtagged_bFlav = fs->make<TH1D>("h_AK5_nbtagged_bFlav" , ";N(b-tagged jets) ;Events" ,20 ,-0.5 ,19.5) ; 
  h_AK5_nbtagged_cFlav = fs->make<TH1D>("h_AK5_nbtagged_cFlav" , ";N(b-tagged jets) ;Events" ,20 ,-0.5 ,19.5) ; 
  h_AK5_nbtagged_lFlav = fs->make<TH1D>("h_AK5_nbtagged_lFlav" , ";N(b-tagged jets) ;Events" ,20 ,-0.5 ,19.5) ;  
  h_AK5_nbtagged_0Flav = fs->make<TH1D>("h_AK5_nbtagged_0Flav" , ";N(b-tagged jets) ;Events" ,20 ,-0.5 ,19.5) ;  

  h_nCA8       = fs->make<TH1D>("h_nCA8"       , ";N(CA8 jets)          ;Events" ,20 ,-0.5 ,19.5) ; 
  h_nCA8_bFlav = fs->make<TH1D>("h_nCA8_bFlav" , ";N(b-tagged CA8 jets) ;Events" ,20 ,-0.5 ,19.5) ; 
  h_nCA8_cFlav = fs->make<TH1D>("h_nCA8_cFlav" , ";N(b-tagged CA8 jets) ;Events" ,20 ,-0.5 ,19.5) ; 
  h_nCA8_lFlav = fs->make<TH1D>("h_nCA8_lFlav" , ";N(b-tagged CA8 jets) ;Events" ,20 ,-0.5 ,19.5) ;  
  h_nCA8_0Flav = fs->make<TH1D>("h_nCA8_0Flav" , ";N(b-tagged CA8 jets) ;Events" ,20 ,-0.5 ,19.5) ;  
  h2_CA8_pt_eta_bFlav = fs->make<TH2D>("h2_CA8_pt_eta_bFlav" ,";p_{T} (b flavour CA8 jet) [GeV]; #eta (CA8 jet); Events"     ,1000 ,0. ,1000. ,80 ,-4. ,4.) ; 
  h2_CA8_pt_eta_cFlav = fs->make<TH2D>("h2_CA8_pt_eta_cFlav" ,";p_{T} (c flavour CA8 jet) [GeV]; #eta (CA8 jet); Events"     ,1000 ,0. ,1000. ,80 ,-4. ,4.) ; 
  h2_CA8_pt_eta_lFlav = fs->make<TH2D>("h2_CA8_pt_eta_lFlav" ,";p_{T} (light flavour CA8 jet) [GeV]; #eta (CA8 jet); Events" ,1000 ,0. ,1000. ,80 ,-4. ,4.) ;  
  h2_CA8_pt_eta_0Flav = fs->make<TH2D>("h2_CA8_pt_eta_0Flav" ,";p_{T} (flavourless CA8 jet) [GeV]; #eta (CA8 jet); Events"   ,1000 ,0. ,1000. ,80 ,-4. ,4.) ;  
  h2_CA8_pt_eta_bFlav_btagged = fs->make<TH2D>("h2_CA8_pt_eta_bFlav_btagged" , ";p_{T} (b flavour CA8 jet) [GeV]; #eta (CA8 jet); Events"     ,1000 ,0. ,1000. ,80 ,-4. ,4.) ; 
  h2_CA8_pt_eta_cFlav_btagged = fs->make<TH2D>("h2_CA8_pt_eta_cFlav_btagged" , ";p_{T} (c flavour CA8 jet) [GeV]; #eta (CA8 jet); Events"     ,1000 ,0. ,1000. ,80 ,-4. ,4.) ; 
  h2_CA8_pt_eta_lFlav_btagged = fs->make<TH2D>("h2_CA8_pt_eta_lFlav_btagged" , ";p_{T} (light flavour CA8 jet) [GeV]; #eta (CA8 jet); Events" ,1000 ,0. ,1000. ,80 ,-4. ,4.) ;  
  h2_CA8_pt_eta_0Flav_btagged = fs->make<TH2D>("h2_CA8_pt_eta_0Flav_btagged" , ";p_{T} (flavourless CA8 jet) [GeV]; #eta (CA8 jet); Events"   ,1000 ,0. ,1000. ,80 ,-4. ,4.) ;  
  h2_CA8_btagEff_pt_eta_bFlav = fs->make<TH2D>("h2_CA8_btagEff_pt_eta_bFlav" , ";p_{T} (b flavour CA8 jet) [GeV]; #eta (CA8 jet); CSVM b-tagging efficiency"     ,1000 ,0. ,1000. ,80 ,-4. ,4.) ; 
  h2_CA8_btagEff_pt_eta_cFlav = fs->make<TH2D>("h2_CA8_btagEff_pt_eta_cFlav" , ";p_{T} (c flavour CA8 jet) [GeV]; #eta (CA8 jet); CSVM b-tagging efficiency"     ,1000 ,0. ,1000. ,80 ,-4. ,4.) ; 
  h2_CA8_btagEff_pt_eta_lFlav = fs->make<TH2D>("h2_CA8_btagEff_pt_eta_lFlav" , ";p_{T} (light flavour CA8 jet) [GeV]; #eta (CA8 jet); CSVM b-tagging efficiency" ,1000 ,0. ,1000. ,80 ,-4. ,4.) ;  
  h2_CA8_btagEff_pt_eta_0Flav = fs->make<TH2D>("h2_CA8_btagEff_pt_eta_0Flav" , ";p_{T} (flvourless CA8 jet) [GeV]; #eta (CA8 jet); CSVM b-tagging efficiency"    ,1000 ,0. ,1000. ,80 ,-4. ,4.) ;  
  h_CA8_btagEff_pt_bFlav = fs->make<TH1D>("h_CA8_btagEff_pt_bFlav" , ";p_{T} (b flavour CA8 jet) [GeV]; CSVM b-tagging efficiency"     ,20 ,0. ,1000.) ; 
  h_CA8_btagEff_pt_cFlav = fs->make<TH1D>("h_CA8_btagEff_pt_cFlav" , ";p_{T} (c flavour CA8 jet) [GeV]; CSVM b-tagging efficiency"     ,20 ,0. ,1000.) ; 
  h_CA8_btagEff_pt_lFlav = fs->make<TH1D>("h_CA8_btagEff_pt_lFlav" , ";p_{T} (light flavour CA8 jet) [GeV]; CSVM b-tagging efficiency" ,20 ,0. ,1000.) ;  
  h_CA8_btagEff_pt_0Flav = fs->make<TH1D>("h_CA8_btagEff_pt_0Flav" , ";p_{T} (flvourless CA8 jet) [GeV]; CSVM b-tagging efficiency"    ,20 ,0. ,1000.) ;  
  h_CA8_nbtagged_bFlav = fs->make<TH1D>("h_CA8_nbtagged_bFlav" , ";N(b-tagged CA8 jets) ;Events" ,20 ,-0.5 ,19.5) ; 
  h_CA8_nbtagged_cFlav = fs->make<TH1D>("h_CA8_nbtagged_cFlav" , ";N(b-tagged CA8 jets) ;Events" ,20 ,-0.5 ,19.5) ; 
  h_CA8_nbtagged_lFlav = fs->make<TH1D>("h_CA8_nbtagged_lFlav" , ";N(b-tagged CA8 jets) ;Events" ,20 ,-0.5 ,19.5) ;  
  h_CA8_nbtagged_0Flav = fs->make<TH1D>("h_CA8_nbtagged_0Flav" , ";N(b-tagged CA8 jets) ;Events" ,20 ,-0.5 ,19.5) ;  

  h2_pt_AK5_CA8_bFlav  = fs->make<TH2D>("h2_pt_AK5_CA8_bFlav"  , ";p_{T} (AK5); p_{T}(CA8); Events" , 1000 ,0. ,1000. ,1000 ,0. ,1000.) ; 
  h2_pt_AK5_CA8_cFlav  = fs->make<TH2D>("h2_pt_AK5_CA8_cFlav"  , ";p_{T} (AK5); p_{T}(CA8); Events" , 1000 ,0. ,1000. ,1000 ,0. ,1000.) ; 
  h2_pt_AK5_CA8_lFlav  = fs->make<TH2D>("h2_pt_AK5_CA8_lFlav"  , ";p_{T} (AK5); p_{T}(CA8); Events" , 1000 ,0. ,1000. ,1000 ,0. ,1000.) ; 
  h2_pt_AK5_CA8_0Flav  = fs->make<TH2D>("h2_pt_AK5_CA8_0Flav"  , ";p_{T} (AK5); p_{T}(CA8); Events" , 1000 ,0. ,1000. ,1000 ,0. ,1000.) ; 
  h_DR_AK5_CA8_bFlav   = fs->make<TH1D>("h_DR_AK5_CA8_bFlav"   , ";#DeltaR(CA8, AK5) ;Events"    ,40 ,0.0  ,4.0) ; 
  h_DR_AK5_CA8_cFlav   = fs->make<TH1D>("h_DR_AK5_CA8_cFlav"   , ";#DeltaR(CA8, AK5) ;Events"    ,40 ,0.0  ,4.0) ; 
  h_DR_AK5_CA8_lFlav   = fs->make<TH1D>("h_DR_AK5_CA8_lFlav"   , ";#DeltaR(CA8, AK5) ;Events"    ,40 ,0.0  ,4.0) ; 
  h_DR_AK5_CA8_0Flav   = fs->make<TH1D>("h_DR_AK5_CA8_0Flav"   , ";#DeltaR(CA8, AK5) ;Events"    ,40 ,0.0  ,4.0) ; 


  return;  

}

double BTagEff::NearestJet( JetCollection* jets, const Jet& myJet, const Jet* nearestJet){
  double DR_min = 1000; 
  if( jets->size() > 0 ){
    for( JetCollection::const_iterator ijet = jets->begin(); ijet != jets->end(); ++ijet ){
      if( myJet.DeltaR(*ijet) < DR_min){
        DR_min = myJet.DeltaR(*ijet);
        nearestJet = &(*ijet) ;  
      }
    }
  }else{
    edm::LogError("BTagEff::NearestJet") << "JetCollection size is " << jets->size() << " !!!" ; 
  }
  return DR_min ; 
}

void BTagEff::isolateCollection( JetCollection control, JetCollection input, JetCollection* output, double dR ){
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
      output->push_back(*ii);
    } 
  }
}

// ------------ method called for each event  ------------
void BTagEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){ 
  using namespace edm;
  using namespace std;

  if(  chain_ == 0) return;

  edm::LogInfo("BTagEff::analyze") << "Starting analysis loop\n";

  for(int entry = 0; entry < maxEvents_; ++entry){ 

    if((entry%reportEvery_) == 0) edm::LogInfo("BTagEff::analyze") << " Processing event " << entry << " of " << maxEvents_ ; 

    chain_->GetEntry(entry);

    isData_   = EvtInfo.McFlag ? 0 : 1; 
    if ( isData_ ) { 
      edm::LogError("BTagEff::analyze") << " Error: Please run this module on MC only!" ; 
      return ; 
    }

    eSelector_->reset(); 
    int code = eSelector_->passCode();
    for ( int ii = 0; ii < (int)cutLevels_.size(); ++ii ) {
      if ( code >= ii ) h_cutflow -> Fill(ii);
      else break;
    }

    if( eSelector_ -> passes() ) { 

      JetCollection goodHiggsJets = eSelector_->goodHiggsJets() ; 
      JetCollection failedHiggsJets = eSelector_->nonHiggsJets() ; 

      int nJets_bFlav(0) ; 
      int nJets_cFlav(0) ; 
      int nJets_lFlav(0) ; 
      int nJets_0Flav(0) ; 
      int nbtaggedJets_bFlav(0) ; 
      int nbtaggedJets_cFlav(0) ; 
      int nbtaggedJets_lFlav(0) ; 
      int nbtaggedJets_0Flav(0) ; 
      JetSelector jetSelAK5(jetSelParams_) ; 
      JetSelector jetSelAK5BTagged(bjetSelParams_) ; 
      pat::strbitset retjetidak5 = jetSelAK5.getBitTemplate() ; 
      pat::strbitset retjetidak5btagged = jetSelAK5BTagged.getBitTemplate() ; 
      JetCollection *ak5jets_bFlav = new JetCollection(); 
      JetCollection *ak5jets_cFlav = new JetCollection(); 
      JetCollection *ak5jets_lFlav = new JetCollection(); 
      JetCollection *ak5jets_bFlav_btagged = new JetCollection(); 
      JetCollection *ak5jets_cFlav_btagged = new JetCollection(); 
      JetCollection *ak5jets_lFlav_btagged = new JetCollection(); 

      for (int ijet = 0; ijet <= JetInfo.Size; ++ijet) { 

        retjetidak5.set(false) ;
        if (jetSelAK5(JetInfo, ijet, goodHiggsJets, retjetidak5) == 0) continue ; 
        Jet thisjet(JetInfo, ijet) ; 

        if (abs(thisjet.GenFlavor()) == 5) { 
          h2_AK5_pt_eta_bFlav->Fill(thisjet.Pt(), thisjet.Eta()); 
          ak5jets_bFlav->push_back(thisjet) ; 
          ++nJets_bFlav ; 
        } 
        else if (abs(thisjet.GenFlavor()) == 4) { 
          h2_AK5_pt_eta_cFlav->Fill(thisjet.Pt(), thisjet.Eta()); 
          ak5jets_cFlav->push_back(thisjet) ;
          ++nJets_cFlav ; 
        } 
        else if (abs(thisjet.GenFlavor()) <  4 || abs(thisjet.GenFlavor()) == 21) { 
          h2_AK5_pt_eta_lFlav->Fill(thisjet.Pt(), thisjet.Eta()); 
          ak5jets_lFlav->push_back(thisjet) ; 
          ++nJets_lFlav ; 
        } 
        else {
          h2_AK5_pt_eta_0Flav->Fill(thisjet.Pt(), thisjet.Eta()); 
          ++nJets_0Flav ; 
        }

        retjetidak5btagged.set(false) ;
        if (jetSelAK5BTagged(JetInfo, ijet, goodHiggsJets, retjetidak5btagged) == 0) continue ; 

        if (abs(thisjet.GenFlavor()) == 5) { 
          h2_AK5_pt_eta_bFlav_btagged->Fill(thisjet.Pt(), thisjet.Eta()); 
          ak5jets_bFlav_btagged->push_back(thisjet) ; 
          ++nbtaggedJets_bFlav ; 
        }  
        else if (abs(thisjet.GenFlavor()) == 4) { 
          h2_AK5_pt_eta_cFlav_btagged->Fill(thisjet.Pt(), thisjet.Eta()); 
          ak5jets_cFlav_btagged->push_back(thisjet) ; 
          ++nbtaggedJets_cFlav ; 
        }  
        else if (abs(thisjet.GenFlavor()) <  4 || abs(thisjet.GenFlavor()) == 21) {
          h2_AK5_pt_eta_lFlav_btagged->Fill(thisjet.Pt(), thisjet.Eta()); 
          ak5jets_lFlav_btagged->push_back(thisjet) ; 
          ++nbtaggedJets_lFlav ; 
        }
        else {
          h2_AK5_pt_eta_0Flav_btagged->Fill(thisjet.Pt(), thisjet.Eta()); 
          ++nbtaggedJets_0Flav ; 
        }
      } //// Loop over AK5 jets 

      h_nAK5_bFlav -> Fill(nJets_bFlav); 
      h_nAK5_cFlav -> Fill(nJets_cFlav); 
      h_nAK5_lFlav -> Fill(nJets_lFlav);  
      h_nAK5_0Flav -> Fill(nJets_0Flav);  

      h_AK5_nbtagged_bFlav -> Fill(nbtaggedJets_bFlav); 
      h_AK5_nbtagged_cFlav -> Fill(nbtaggedJets_cFlav); 
      h_AK5_nbtagged_lFlav -> Fill(nbtaggedJets_lFlav);  
      h_AK5_nbtagged_0Flav -> Fill(nbtaggedJets_0Flav);  

      int nCA8_bFlav(0) ; 
      int nCA8_cFlav(0) ; 
      int nCA8_lFlav(0) ; 
      int nCA8_0Flav(0) ; 
      int nbtaggedCA8_bFlav(0) ; 
      int nbtaggedCA8_cFlav(0) ; 
      int nbtaggedCA8_lFlav(0) ; 
      int nbtaggedCA8_0Flav(0) ; 
      FatJetSelector jetSelCA8(fatJetSelParams_) ; 
      FatJetSelector jetSelCA8BTagged(fatbjetSelParams_) ; 
      pat::strbitset retjetidca8 = jetSelCA8.getBitTemplate() ; 
      pat::strbitset retjetidca8btagged = jetSelCA8BTagged.getBitTemplate() ; 
      JetCollection* ca8jets_bFlav = new JetCollection(); 
      JetCollection* ca8jets_cFlav = new JetCollection(); 
      JetCollection* ca8jets_lFlav = new JetCollection(); 
      JetCollection* ca8jets_bFlav_btagged = new JetCollection(); 
      JetCollection* ca8jets_cFlav_btagged = new JetCollection(); 
      JetCollection* ca8jets_lFlav_btagged = new JetCollection(); 

      h_nCA8->Fill(failedHiggsJets.size()) ; 

      for (JetCollection::const_iterator ica8 = failedHiggsJets.begin(); ica8 != failedHiggsJets.end(); ++ica8) { 

        if (abs(ica8->GenFlavor()) == 5) { 
          ca8jets_bFlav->push_back(*ica8) ; 
          ++nCA8_bFlav ; 
          TLorentzVector p4_ca8(ica8->Px(), ica8->Py(), ica8->Pz(), ica8->Energy()) ; 
          for (JetCollection::const_iterator iak5 = ak5jets_bFlav->begin(); iak5 != ak5jets_bFlav->end(); ++iak5) {
            TLorentzVector p4_ak5(iak5->Px(), iak5->Py(), iak5->Pz(), iak5->Energy()) ; 
            double dr = p4_ca8.DeltaR(p4_ak5) ; 
            h_DR_AK5_CA8_bFlav->Fill(dr) ; 
            if (dr < 0.05) {
              h2_CA8_pt_eta_bFlav->Fill(ica8->Pt(), ica8->Eta()); 
              h2_pt_AK5_CA8_bFlav->Fill(iak5->Pt(), ica8->Pt()) ; 
              if (iak5->CombinedSVBJetTags() > 0.679) {
                h2_CA8_pt_eta_bFlav_btagged->Fill(ica8->Pt(), ica8->Eta()) ; 
              ++nbtaggedCA8_bFlav ; 
              ca8jets_bFlav_btagged->push_back(*ica8) ; 
              }
            }
          }
          //for (JetCollection::const_iterator iak5 = ak5jets_bFlav_btagged->begin(); iak5 != ak5jets_bFlav_btagged->end(); ++iak5) {
          //  TLorentzVector p4_ak5(iak5->Px(), iak5->Py(), iak5->Pz(), iak5->Energy()) ; 
          //  double dr = p4_ca8.DeltaR(p4_ak5) ; 
          //  h_DR_AK5_CA8_bFlav->Fill(dr) ; 
          //  if (dr < 0.05) {
          //    h2_CA8_pt_eta_bFlav_btagged->Fill(ica8->Pt(), ica8->Eta()) ; 
          //    h2_pt_AK5_CA8_bFlav->Fill(iak5->Pt(), ica8->Pt()) ; 
          //    ++nbtaggedCA8_bFlav ; 
          //    ca8jets_bFlav_btagged->push_back(*ica8) ; 
          //  }
          //}
        } 
        else if (abs(ica8->GenFlavor()) == 4) { 
          h2_CA8_pt_eta_cFlav->Fill(ica8->Pt(), ica8->Eta()); 
          ca8jets_cFlav->push_back(*ica8) ;
          ++nCA8_cFlav ; 
          TLorentzVector p4_ca8(ica8->Px(), ica8->Py(), ica8->Pz(), ica8->Energy()) ; 
          for (JetCollection::const_iterator iak5 = ak5jets_cFlav_btagged->begin(); iak5 != ak5jets_cFlav_btagged->end(); ++iak5) {
            TLorentzVector p4_ak5(iak5->Px(), iak5->Py(), iak5->Pz(), iak5->Energy()) ; 
            double dr = p4_ca8.DeltaR(p4_ak5) ; 
            h_DR_AK5_CA8_cFlav->Fill(dr) ; 
            if (dr < 0.05) {
              h2_CA8_pt_eta_cFlav_btagged->Fill(ica8->Pt(), ica8->Eta()) ; 
              h2_pt_AK5_CA8_cFlav->Fill(iak5->Pt(), ica8->Pt()) ; 
              ++nbtaggedCA8_cFlav ; 
              ca8jets_cFlav_btagged->push_back(*ica8) ; 
            }
          }
        } 
        else if (abs(ica8->GenFlavor()) <  4 || abs(ica8->GenFlavor()) == 21) { 
          h2_CA8_pt_eta_lFlav->Fill(ica8->Pt(), ica8->Eta()); 
          ca8jets_lFlav->push_back(*ica8) ; 
          ++nCA8_lFlav ; 
          TLorentzVector p4_ca8(ica8->Px(), ica8->Py(), ica8->Pz(), ica8->Energy()) ; 
          for (JetCollection::const_iterator iak5 = ak5jets_lFlav_btagged->begin(); iak5 != ak5jets_lFlav_btagged->end(); ++iak5) {
            TLorentzVector p4_ak5(iak5->Px(), iak5->Py(), iak5->Pz(), iak5->Energy()) ; 
            double dr = p4_ca8.DeltaR(p4_ak5) ; 
            h_DR_AK5_CA8_lFlav->Fill(dr) ; 
            if (dr < 0.05) {
              h2_CA8_pt_eta_lFlav_btagged->Fill(ica8->Pt(), ica8->Eta()) ; 
              h2_pt_AK5_CA8_lFlav->Fill(iak5->Pt(), ica8->Pt()) ; 
              ++nbtaggedCA8_lFlav ; 
              ca8jets_lFlav_btagged->push_back(*ica8) ; 
            }
          }
        } 
        else {
          h2_CA8_pt_eta_0Flav->Fill(ica8->Pt(), ica8->Eta()); 
          ++nCA8_0Flav ; 
        }

      } //// Loop over CA8 jets 

      h_nCA8_bFlav -> Fill(nCA8_bFlav); 
      h_nCA8_cFlav -> Fill(nCA8_cFlav); 
      h_nCA8_lFlav -> Fill(nCA8_lFlav);  
      h_nCA8_0Flav -> Fill(nCA8_0Flav);  

    } /// If event passes selection 

  } //// entry loop 

}

// ------------ method called once each job just after ending the event loop  ------------
void BTagEff::endJob(){ 

  h2_AK5_btagEff_pt_eta_bFlav -> Divide(h2_AK5_pt_eta_bFlav_btagged, h2_AK5_pt_eta_bFlav); 
  h2_AK5_btagEff_pt_eta_cFlav -> Divide(h2_AK5_pt_eta_cFlav_btagged, h2_AK5_pt_eta_cFlav); 
  h2_AK5_btagEff_pt_eta_lFlav -> Divide(h2_AK5_pt_eta_lFlav_btagged, h2_AK5_pt_eta_lFlav);  
  h2_AK5_btagEff_pt_eta_0Flav -> Divide(h2_AK5_pt_eta_0Flav_btagged, h2_AK5_pt_eta_0Flav);  

  h_AK5_btagEff_pt_bFlav -> Divide(h2_AK5_pt_eta_bFlav_btagged->ProjectionX()->Rebin(50), h2_AK5_pt_eta_bFlav->ProjectionX()->Rebin(50)); 
  h_AK5_btagEff_pt_cFlav -> Divide(h2_AK5_pt_eta_cFlav_btagged->ProjectionX()->Rebin(50), h2_AK5_pt_eta_cFlav->ProjectionX()->Rebin(50)); 
  h_AK5_btagEff_pt_lFlav -> Divide(h2_AK5_pt_eta_lFlav_btagged->ProjectionX()->Rebin(50), h2_AK5_pt_eta_lFlav->ProjectionX()->Rebin(50));  
  h_AK5_btagEff_pt_0Flav -> Divide(h2_AK5_pt_eta_0Flav_btagged->ProjectionX()->Rebin(50), h2_AK5_pt_eta_0Flav->ProjectionX()->Rebin(50));  

  h_CA8_btagEff_pt_bFlav -> Divide(h2_CA8_pt_eta_bFlav_btagged->ProjectionX()->Rebin(50), h2_CA8_pt_eta_bFlav->ProjectionX()->Rebin(50)) ;
  h_CA8_btagEff_pt_cFlav -> Divide(h2_CA8_pt_eta_cFlav_btagged->ProjectionX()->Rebin(50), h2_CA8_pt_eta_cFlav->ProjectionX()->Rebin(50)) ;
  h_CA8_btagEff_pt_lFlav -> Divide(h2_CA8_pt_eta_lFlav_btagged->ProjectionX()->Rebin(50), h2_CA8_pt_eta_lFlav->ProjectionX()->Rebin(50)) ;

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void BTagEff::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BTagEff);
