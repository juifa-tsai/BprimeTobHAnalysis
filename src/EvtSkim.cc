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
//#include "BpbH/BprimeTobH/interface/Njettiness.hh"
//#include "BpbH/BprimeTobH/interface/Nsubjettiness.hh"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h" 

#include "BpbH/BprimeTobHAnalysis/interface/EventSelector.h"

#include "BpbH/BprimeTobHAnalysis/interface/reRegistGen.hh"
#include "BpbH/BprimeTobHAnalysis/interface/reRegistJet.hh"

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

    const edm::ParameterSet         evtSelParams_; 

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
  evtSelParams_(iConfig.getParameter<edm::ParameterSet>("EvtSelParams")),
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

  ofstream fout("Evt_NoJets.txt"); 

  edm::LogInfo("StartingAnalysisLoop") << "Starting analysis loop\n";

  for(int entry=0; entry<maxEvents_; entry++){
    if( (entry%reportEvery_) == 0) edm::LogInfo("Event") << entry << " of " << maxEvents_; 
    chain_->GetEntry(entry);

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

  fout.close(); 

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
