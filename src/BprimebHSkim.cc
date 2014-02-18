// -*- C++ -*-
//
// Package:    BprimebHSkim
// Class:      BprimebHSkim
// 
/**\class BprimebHSkim BprimebHSkim.cc Bprime_kit/BprimebHSkim/src/BprimebHSkim.cc

Description: 
Event skimming code for Bprime->bH studies. 
Input: BprimebH Ntuples.
Output: Skimmed BprimebH Ntuples. 
- National Taiwan University - 

Implementation:
[Notes on implementation]
*/
//
// Original Author: Devdatta Majumder  
//         Created:  
// $Id$
//
//

// system include files
#include <memory>
#include <iostream>
#include <ctime>
#include <fstream>
#include <assert.h>

// Root headers 
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "BpbH/BprimeTobHAnalysis/interface/EventSelector.h"

//
// class declaration
//

class BprimebHSkim : public edm::EDAnalyzer{
  public:
    explicit BprimebHSkim(const edm::ParameterSet&);
    ~BprimebHSkim();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();

    // ----------member data ---------------------------

    //// Configurables 

    int                             maxEvents_; 
    const int                       reportEvery_; 
    const std::string               inputTTree_;
    const std::vector<std::string>  inputFiles_;

    edm::Service<TFileService> fs; 

    TChain*            chain_;
    TTree*		         newtree;	

    VertexInfoBranches VtxInfo;
    JetInfoBranches    JetInfo;
    JetInfoBranches    FatJetInfo;

    TH1D*   h_cutflow ; 

};

//
// constructors and destructor
//
BprimebHSkim::BprimebHSkim(const edm::ParameterSet& iConfig) : 
  maxEvents_(iConfig.getParameter<int>("MaxEvents")), 
  reportEvery_(iConfig.getParameter<int>("ReportEvery")),
  inputTTree_(iConfig.getParameter<std::string>("InputTTree")),
  inputFiles_(iConfig.getParameter<std::vector<std::string> >("InputFiles")), 
  newtree(0)  
{ 

}


BprimebHSkim::~BprimebHSkim(){ 
  delete chain_;
}

// ------------ method called once each job just before starting event loop  ------------
void BprimebHSkim::beginJob(){ 

  chain_  = new TChain(inputTTree_.c_str());

  h_cutflow = fs->make<TH1D>("h_cutflow" ,"Cut flow" ,4 ,0. ,4); 
  h_cutflow->GetXaxis()->SetBinLabel(1, "BeforeNtuplizer") ; 
  h_cutflow->GetXaxis()->SetBinLabel(2, "AfterNtuplizer") ; 
  h_cutflow->GetXaxis()->SetBinLabel(3, "BeforeSkim") ; 
  h_cutflow->GetXaxis()->SetBinLabel(4, "AfterSkim") ; 

  for(unsigned i=0; i<inputFiles_.size(); ++i){
    chain_->Add(inputFiles_.at(i).c_str());
    TFile *f = TFile::Open(inputFiles_.at(i).c_str(),"READ");
    TH1F* h_events = (TH1F*)f->Get("ntuple/h_events") ; 
    if ( i == 0 ) VtxInfo.Register(chain_);
    if ( i == 0 ) JetInfo.Register(chain_,"JetInfo");
    if ( i == 0 ) FatJetInfo.Register(chain_,"FatJetInfo");
    f->Close();
    if ( maxEvents_ < 0 || maxEvents_ >= chain_->GetEntries() ) {
      h_cutflow->AddBinContent(1, h_events->GetBinContent(1)) ; 
      h_cutflow->AddBinContent(2, h_events->GetBinContent(2)) ; 
    } 
  }

  if(  maxEvents_<0 || maxEvents_>chain_->GetEntries()) maxEvents_ = chain_->GetEntries();

  fs->cd() ;
  newtree = chain_->CloneTree(0) ; 
  return;  

}

// ------------ method called for each event  ------------
void BprimebHSkim::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){ 
  using namespace edm;
  using namespace std;

  if(  chain_ == 0) return;

  for(int entry = 0; entry < maxEvents_; ++entry){ 

    if((entry%reportEvery_) == 0) edm::LogInfo("Event") << entry << " of " << maxEvents_ ; 

    chain_->GetEntry(entry);

    h_cutflow -> Fill("BeforeSkim", 1) ;

    int nGoodVtxs(0) ;
    VertexSelector vtxSel(VtxInfo) ; 
    nGoodVtxs = vtxSel.NGoodVtxs(); 

    bool hasak5jet(false) ; 
    for (int iak5jet = 0; iak5jet < JetInfo.Size; ++iak5jet) { 
      if (JetInfo.Pt[iak5jet] > 30) {
        hasak5jet = true; 
        break ; 
      }
    }

    bool hasca8jet(false) ; 
    for (int ica8jet = 0; ica8jet < FatJetInfo.Size; ++ica8jet) { 
      if (FatJetInfo.Pt[ica8jet] > 200) {
        hasca8jet = true; 
        break ; 
      }
    }

    if (nGoodVtxs >= 1 && hasak5jet && hasca8jet ) {
      newtree->Fill(); 
      h_cutflow -> Fill("AfterSkim", 1) ;
    }

  } //// entry loop 

}

// ------------ method called once each job just after ending the event loop  ------------
void BprimebHSkim::endJob(){ 
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void BprimebHSkim::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BprimebHSkim);
