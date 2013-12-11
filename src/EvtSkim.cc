// -*- C++ -*-
//
// Package:    EvtSkim
// Class:      EvtSkim
// 
/**\class EvtSkim EvtSkim.cc Bprime_kit/EvtSkim/src/EvtSkim.cc

Description: 
Analyzer class for Bprime -> b Higgs studies 
- National Taiwan University - 

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Jui-Fa Tsai 
//         Created:  
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
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
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
#include "BpbH/BprimeTobHAnalysis/interface/reRegistGen.hh"
#include "BpbH/BprimeTobHAnalysis/interface/reRegistJet.hh"
#include "BpbH/BprimeTobHAnalysis/interface/BTagSFUtil-tprime.h"
#include "BpbH/BprimeTobHAnalysis/interface/GetBTag_SF_EFF.h"

//
// class declaration
//

class EvtSkim : public edm::EDAnalyzer{
  public:
    explicit EvtSkim(const edm::ParameterSet&);
    ~EvtSkim();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();

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

    const edm::ParameterSet         bjetSelParams_ ; 
    const edm::ParameterSet         bTagSFUtilParameters_; 
    const edm::ParameterSet         evtSelParams_; 

    bool                            modifyBTags_ ; 
    int                             nJetsToModify_ ; 

    TChain*            chain_;
    TTree*		         newtree;	

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
    TH1D* 	AK5_num;
    TH1D* 	AK5_pt;
    TH1D* 	AK5_CSV;
    TH1D* 	bJet_num;
    TH1D* 	bJet_pt;
    TH1D* 	bJet_CSV;
    TH1D* 	bJetVeto_num; 
    TH1D* 	bJetVetoMatchCA8_num;

};

//
// constructors and destructor
//
EvtSkim::EvtSkim(const edm::ParameterSet& iConfig) : 
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
  bjetSelParams_(iConfig.getParameter<edm::ParameterSet>("BJetSelParams")), 
  bTagSFUtilParameters_(iConfig.getParameter<edm::ParameterSet>("BTagSFUtilParameters")),  
  evtSelParams_(iConfig.getParameter<edm::ParameterSet>("EvtSelParams")),
  modifyBTags_(iConfig.getParameter<bool>("ModifyBTags")), 
  nJetsToModify_(iConfig.getParameter<int>("NJetsToModify")), 
  isData_(0),
  evtwt_(1), 
  puweight_(1)  
{ 

  if( doPUReweighting_) LumiWeights_ = edm::LumiReWeighting(file_PUDistMC_, file_PUDistData_, hist_PUDistMC_, hist_PUDistData_);

}


EvtSkim::~EvtSkim(){ 
  delete chain_;
}

// ------------ method called once each job just before starting event loop  ------------
void EvtSkim::beginJob(){ 
  chain_  = new TChain(inputTTree_.c_str());
  newtree = fs->make<TTree>("tree", "") ; 

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

  newtree->Branch("EvtInfo.McFlag", &McFlag_, "EvtInfo.McFlag/O"); // Store weight of Evt and PU for each event
  newtree->Branch("EvtInfo.WeightEvtPU", &evtwtPu_, "EvtInfo.WeightEvtPU/D"); // Store weight of Evt and PU for each event
  newGenInfo.RegisterTree(newtree);
  newJetInfo.RegisterTree(newtree,"JetInfo");
  newFatJetInfo.RegisterTree(newtree,"FatJetInfo");
  newSubJetInfo.RegisterTree(newtree,"SubJetInfo");

  if(  maxEvents_<0 || maxEvents_>chain_->GetEntries()) maxEvents_ = chain_->GetEntries();

  eSelector_ = new EventSelector(evtSelParams_,EvtInfo,VtxInfo,JetInfo,FatJetInfo,SubJetInfo);  
  cutLevels_ = eSelector_->getCutLevels();

  h_cutflow = fs->make<TH1D>("h_cutflow" ,"Cut flow" ,(int)cutLevels_.size() ,0. ,(int)cutLevels_.size() ); 
  for (int ii = 0; ii < (int)cutLevels_.size(); ++ii) {
    h_cutflow -> GetXaxis() -> SetBinLabel(ii+1, (cutLevels_.at(ii)).c_str()) ;  
  }
  AK5_num		= fs->make<TH1D>("AK5JetInfo.Num",	"", 10,   0, 10); 
  AK5_pt 		= fs->make<TH1D>("AK5JetInfo.Pt",	"", 1500, 0, 1500);
  AK5_CSV		= fs->make<TH1D>("AK5JetInfo.CSV", 	"", 100,  0, 1.);
  bJet_num 	= fs->make<TH1D>("bJetInfo.Num",	"", 10, 0, 10); 
  bJet_pt		= fs->make<TH1D>("bJetInfo.Pt",		"", 1500, 0, 1500);
  bJet_CSV	= fs->make<TH1D>("bJetInfo.CSV", 	"", 100,  0, 1.);
  bJetVeto_num 		= fs->make<TH1D>("bJetInfo.Num.Veto",		"", 10, 0, 10); 
  bJetVetoMatchCA8_num 	= fs->make<TH1D>("bJetInfo.NumMatchToCA8.Veto",	"", 10, 0, 10);

  return;  

}

// ------------ method called for each event  ------------
void EvtSkim::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){ 
  using namespace edm;
  using namespace std;

  if(  chain_ == 0) return;

  std::string btagAlgo ;
  if (bjetSelParams_.getParameter<double>("jetCSVDiscMin") >= 0.244
      && bjetSelParams_.getParameter<double>("jetCSVDiscMin") < 0.679) btagAlgo = "CSVL" ;  
  else if (bjetSelParams_.getParameter<double>("jetCSVDiscMin") >= 0.679
      && bjetSelParams_.getParameter<double>("jetCSVDiscMin") <  0.844) btagAlgo = "CSVM" ; 
  else if (bjetSelParams_.getParameter<double>("jetCSVDiscMin") >= 0.844) btagAlgo = "CSVT" ;  

  JetSelector jetSelBTaggedAK5(bjetSelParams_) ; 
  pat::strbitset retjetidbtaggedak5 = jetSelBTaggedAK5.getBitTemplate() ; 
  BTagSFUtil* btsf = new BTagSFUtil(bTagSFUtilParameters_) ; 
  btsf -> setSeed(1) ;  

  for(int entry = 0; entry < maxEvents_; ++entry){ 

    chain_->GetEntry(entry);

    JetCollection ak5jets ; 
    std::vector< std::pair<float, float> > jet_pt_eta ; 
    if (EvtInfo.McFlag && modifyBTags_ ) {
      btsf->setSeed(entry*1e+12+EvtInfo.RunNo*1e+6+EvtInfo.EvtNo);

      for ( int ijet = 0; ijet < JetInfo.Size; ++ijet ) {
        std::pair<float, float> thisjet_pt_eta(JetInfo.Et[ijet],JetInfo.Eta[ijet]) ; 
        jet_pt_eta.push_back(thisjet_pt_eta) ; 
        Jet thisjet(JetInfo, ijet) ; 
        retjetidbtaggedak5.set(false) ;
        bool isbtagged = (bool)jetSelBTaggedAK5(JetInfo, ijet,retjetidbtaggedak5) ; 
        if (ijet < nJetsToModify_) thisjet.setIsBTagged(btagAlgo, isbtagged) ; 
        ak5jets.push_back(thisjet) ; 
      }
      btsf -> readDB(iSetup,jet_pt_eta) ; 

    }

    for (int ijet = 0; ijet < (int)ak5jets.size(); ++ijet) { 
      if ( ijet >= nJetsToModify_ ) break ; 
        edm::LogInfo("EvtSkim") << " btag algo = " << btagAlgo  ; 
      double btag_sf = btsf->getSF("MUJETSWPBTAG" + btagAlgo ,ijet) ;  
      double btag_eff = btsf->BtagEff_[0] ;  
      if ( btag_sf < 0 || btag_eff < 0 ) { 
        edm::LogWarning("EvtSkim") << " ModifyBTagsWithSF: SF < 0, something is wrong (maybe just out of range where SF measured). Doing nothing. " ; 
        continue ; 
      }
      if ( btagAlgo == "CSVL" ) ak5jets[ijet].IsBtaggedCSVL() ; 
      else if ( btagAlgo == "CSVM" ) ak5jets[ijet].IsBtaggedCSVM() ; 
      else if ( btagAlgo == "CSVT" ) ak5jets[ijet].IsBtaggedCSVT() ; 
      bool jetisbtag ; 
      btsf->modifyBTagsWithSF(jetisbtag, ak5jets[ijet].GenFlavor(), btag_sf, btag_eff);
      if ( btagAlgo == "CSVL" ) {
        if( jetisbtag != ak5jets[ijet].IsBtaggedCSVL() ) edm::LogInfo("EvtSkim") << " Jet CVSL btag been modified. " ; 
      }
      else if ( btagAlgo == "CSVM" ) {
        if( jetisbtag != ak5jets[ijet].IsBtaggedCSVM() ) edm::LogInfo("EvtSkim") << " Jet CVSM btag been modified. " ; 
      }
      else if ( btagAlgo == "CSVT" ) {
        if( jetisbtag != ak5jets[ijet].IsBtaggedCSVT() ) edm::LogInfo("EvtSkim") << " Jet CVST btag been modified. " ; 
      } 
      ak5jets[ijet].setIsBTagged(btagAlgo, jetisbtag) ; 

    }

    eSelector_->reset(); 
    int code = eSelector_->passCode();
    for ( int ii = 0; ii < (int)cutLevels_.size(); ++ii ) {
      if ( code >= ii ) h_cutflow -> Fill(ii);
      else break;
    }

    if( eSelector_ -> passes() ) { 
      if ( isData_ ) {
        McFlag_=0;
        evtwtPu_=evtwt_;
        reRegistJet(JetInfo,newJetInfo);	
        reRegistJet(FatJetInfo,newFatJetInfo);	
        reRegistJet(SubJetInfo,newSubJetInfo);	
      } else {
        McFlag_=1;
        evtwtPu_=evtwt_;
        reRegistGen(GenInfo,newGenInfo); 	
        reRegistJet(JetInfo,newJetInfo);	
        reRegistJet(FatJetInfo,newFatJetInfo);	
        reRegistJet(SubJetInfo,newSubJetInfo);	
      }
      newtree->Fill();
    }

  } //// entry loop 

}

// ------------ method called once each job just after ending the event loop  ------------
void EvtSkim::endJob(){ 
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void EvtSkim::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EvtSkim);
