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
    TH2D*   h2_ptJet_etaJet_bFlav ; 
    TH2D*   h2_ptJet_etaJet_cFlav ; 
    TH2D*   h2_ptJet_etaJet_lFlav ;  
    TH2D*   h2_ptJet_etaJet_0Flav ;  
    TH1D*   h_nJets_bFlav ; 
    TH1D*   h_nJets_cFlav ; 
    TH1D*   h_nJets_lFlav ;  
    TH1D*   h_nJets_0Flav ;  
    TH2D*   h2_ptJet_etaJet_bFlav_btagged ; 
    TH2D*   h2_ptJet_etaJet_cFlav_btagged ; 
    TH2D*   h2_ptJet_etaJet_lFlav_btagged ;  
    TH2D*   h2_ptJet_etaJet_0Flav_btagged ;  
    TH2D*   h2_btagEff_ptJet_etaJet_bFlav ; 
    TH2D*   h2_btagEff_ptJet_etaJet_cFlav ; 
    TH2D*   h2_btagEff_ptJet_etaJet_lFlav ;  
    TH2D*   h2_btagEff_ptJet_etaJet_0Flav ;  
    TH1D*   h_nbtaggedJets_bFlav ; 
    TH1D*   h_nbtaggedJets_cFlav ; 
    TH1D*   h_nbtaggedJets_lFlav ;  
    TH1D*   h_nbtaggedJets_0Flav ;  

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

  h2_ptJet_etaJet_bFlav         = fs->make<TH2D>("h2_ptJet_etaJet_bFlav"         , ";p_{T} (b flavour jet) [GeV]; #eta (jet); Events"     ,1000 ,0. ,1000. ,80 ,-4. ,4.) ; 
  h2_ptJet_etaJet_cFlav         = fs->make<TH2D>("h2_ptJet_etaJet_cFlav"         , ";p_{T} (c flavour jet) [GeV]; #eta (jet); Events"     ,1000 ,0. ,1000. ,80 ,-4. ,4.) ; 
  h2_ptJet_etaJet_lFlav         = fs->make<TH2D>("h2_ptJet_etaJet_lFlav"         , ";p_{T} (light flavour jet) [GeV]; #eta (jet); Events" ,1000 ,0. ,1000. ,80 ,-4. ,4.) ;  
  h2_ptJet_etaJet_0Flav         = fs->make<TH2D>("h2_ptJet_etaJet_0Flav"         , ";p_{T} (flavourless jet) [GeV]; #eta (jet); Events" ,1000 ,0. ,1000. ,80 ,-4. ,4.) ;  
  h_nJets_bFlav = fs->make<TH1D>("h_nJets_bFlav" , ";N(b-tagged jets) ;Events" ,20 ,-0.5 ,19.5) ; 
  h_nJets_cFlav = fs->make<TH1D>("h_nJets_cFlav" , ";N(b-tagged jets) ;Events" ,20 ,-0.5 ,19.5) ; 
  h_nJets_lFlav = fs->make<TH1D>("h_nJets_lFlav" , ";N(b-tagged jets) ;Events" ,20 ,-0.5 ,19.5) ;  
  h_nJets_0Flav = fs->make<TH1D>("h_nJets_0Flav" , ";N(b-tagged jets) ;Events" ,20 ,-0.5 ,19.5) ;  
  h2_ptJet_etaJet_bFlav_btagged = fs->make<TH2D>("h2_ptJet_etaJet_bFlav_btagged" , ";p_{T} (b flavour jet) [GeV]; #eta (jet); Events"     ,1000 ,0. ,1000. ,80 ,-4. ,4.) ; 
  h2_ptJet_etaJet_cFlav_btagged = fs->make<TH2D>("h2_ptJet_etaJet_cFlav_btagged" , ";p_{T} (c flavour jet) [GeV]; #eta (jet); Events"     ,1000 ,0. ,1000. ,80 ,-4. ,4.) ; 
  h2_ptJet_etaJet_lFlav_btagged = fs->make<TH2D>("h2_ptJet_etaJet_lFlav_btagged" , ";p_{T} (light flavour jet) [GeV]; #eta (jet); Events" ,1000 ,0. ,1000. ,80 ,-4. ,4.) ;  
  h2_ptJet_etaJet_0Flav_btagged = fs->make<TH2D>("h2_ptJet_etaJet_0Flav_btagged" , ";p_{T} (flavourless jet) [GeV]; #eta (jet); Events" ,1000 ,0. ,1000. ,80 ,-4. ,4.) ;  
  h2_btagEff_ptJet_etaJet_bFlav = fs->make<TH2D>("h2_btagEff_ptJet_etaJet_bFlav" , ";p_{T} (b flavour jet) [GeV]; #eta (jet); CSVM b-tagging efficiency"     ,1000 ,0. ,1000. ,80 ,-4. ,4.) ; 
  h2_btagEff_ptJet_etaJet_cFlav = fs->make<TH2D>("h2_btagEff_ptJet_etaJet_cFlav" , ";p_{T} (c flavour jet) [GeV]; #eta (jet); CSVM b-tagging efficiency"     ,1000 ,0. ,1000. ,80 ,-4. ,4.) ; 
  h2_btagEff_ptJet_etaJet_lFlav = fs->make<TH2D>("h2_btagEff_ptJet_etaJet_lFlav" , ";p_{T} (light flavour jet) [GeV]; #eta (jet); CSVM b-tagging efficiency" ,1000 ,0. ,1000. ,80 ,-4. ,4.) ;  
  h2_btagEff_ptJet_etaJet_0Flav = fs->make<TH2D>("h2_btagEff_ptJet_etaJet_0Flav" , ";p_{T} (flvourless jet) [GeV]; #eta (jet); CSVM b-tagging efficiency" , 1000 ,0. ,1000. ,80 ,-4. ,4.) ;  
  h_nbtaggedJets_bFlav = fs->make<TH1D>("h_nbtaggedJets_bFlav" , ";N(b-tagged jets) ;Events" ,20 ,-0.5 ,19.5) ; 
  h_nbtaggedJets_cFlav = fs->make<TH1D>("h_nbtaggedJets_cFlav" , ";N(b-tagged jets) ;Events" ,20 ,-0.5 ,19.5) ; 
  h_nbtaggedJets_lFlav = fs->make<TH1D>("h_nbtaggedJets_lFlav" , ";N(b-tagged jets) ;Events" ,20 ,-0.5 ,19.5) ;  
  h_nbtaggedJets_0Flav = fs->make<TH1D>("h_nbtaggedJets_0Flav" , ";N(b-tagged jets) ;Events" ,20 ,-0.5 ,19.5) ;  

  return;  

}

// ------------ method called for each event  ------------
void BTagEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){ 
  using namespace edm;
  using namespace std;

  if(  chain_ == 0) return;

  edm::LogInfo("BTagEff") << "Starting analysis loop\n";

  for(int entry = 0; entry < maxEvents_; ++entry){ 

    if((entry%reportEvery_) == 0) edm::LogInfo("BTagEff") << " Processing event " << entry << " of " << maxEvents_ ; 

    chain_->GetEntry(entry);

    isData_   = EvtInfo.McFlag ? 0 : 1; 
    if ( isData_ ) { 
      edm::LogError("BTagEff") << " Error: Please run this module on MC only!" ; 
      return ; 
    }

    eSelector_->reset(); 
    int code = eSelector_->passCode();
    for ( int ii = 0; ii < (int)cutLevels_.size(); ++ii ) {
      if ( code >= ii ) h_cutflow -> Fill(ii);
      else break;
    }

    if( eSelector_ -> passes() ) { 

      int nJets_bFlav(0) ; 
      int nJets_cFlav(0) ; 
      int nJets_lFlav(0) ; 
      int nJets_0Flav(0) ; 
      int nbtaggedJets_bFlav(0) ; 
      int nbtaggedJets_cFlav(0) ; 
      int nbtaggedJets_lFlav(0) ; 
      int nbtaggedJets_0Flav(0) ; 
      JetCollection goodHiggsJets = eSelector_->goodHiggsJets() ; 
      JetSelector jetSelAK5(jetSelParams_) ; 
      JetSelector jetSelAK5BTagged(bjetSelParams_) ; 
      pat::strbitset retjetidak5 = jetSelAK5.getBitTemplate() ; 
      pat::strbitset retjetidak5btagged = jetSelAK5BTagged.getBitTemplate() ; 
      JetCollection ak5jets_bFlav ; 
      JetCollection ak5jets_cFlav ; 
      JetCollection ak5jets_lFlav ; 
      JetCollection ak5jets_bFlav_btagged ; 
      JetCollection ak5jets_cFlav_btagged ; 
      JetCollection ak5jets_lFlav_btagged ; 

      for (int ijet = 0; ijet <= JetInfo.Size; ++ijet) { 

        retjetidak5.set(false) ;
        if (jetSelAK5(JetInfo, ijet, goodHiggsJets, retjetidak5) == 0) continue ; 
        Jet thisjet(JetInfo, ijet) ; 

        if (abs(thisjet.GenFlavor()) == 5) { 
          h2_ptJet_etaJet_bFlav->Fill(thisjet.Pt(), thisjet.Eta()); 
          ak5jets_bFlav.push_back(thisjet) ; 
          ++nJets_bFlav ; 
        } 
        else if (abs(thisjet.GenFlavor()) == 4) { 
          h2_ptJet_etaJet_cFlav->Fill(thisjet.Pt(), thisjet.Eta()); 
          ak5jets_cFlav.push_back(thisjet) ;
          ++nJets_cFlav ; 
        } 
        else if (abs(thisjet.GenFlavor()) <  4 || abs(thisjet.GenFlavor()) == 21) { 
          h2_ptJet_etaJet_lFlav->Fill(thisjet.Pt(), thisjet.Eta()); 
          ak5jets_lFlav.push_back(thisjet) ; 
          ++nJets_lFlav ; 
        } 
        else {
          h2_ptJet_etaJet_0Flav->Fill(thisjet.Pt(), thisjet.Eta()); 
          ++nJets_0Flav ; 
        }

        retjetidak5btagged.set(false) ;
        if (jetSelAK5BTagged(JetInfo, ijet, goodHiggsJets, retjetidak5btagged) == 0) continue ; 

        if (abs(thisjet.GenFlavor()) == 5) { 
          h2_ptJet_etaJet_bFlav_btagged->Fill(thisjet.Pt(), thisjet.Eta()); 
          ak5jets_bFlav_btagged.push_back(thisjet) ; 
          ++nbtaggedJets_bFlav ; 
        }  
        else if (abs(thisjet.GenFlavor()) == 4) { 
          h2_ptJet_etaJet_cFlav_btagged->Fill(thisjet.Pt(), thisjet.Eta()); 
          ak5jets_cFlav_btagged.push_back(thisjet) ; 
          ++nbtaggedJets_cFlav ; 
        }  
        else if (abs(thisjet.GenFlavor()) <  4 || abs(thisjet.GenFlavor()) == 21) {
          h2_ptJet_etaJet_lFlav_btagged->Fill(thisjet.Pt(), thisjet.Eta()); 
          ak5jets_lFlav_btagged.push_back(thisjet) ; 
          ++nbtaggedJets_lFlav ; 
        }
        else {
          h2_ptJet_etaJet_0Flav_btagged->Fill(thisjet.Pt(), thisjet.Eta()); 
          ++nbtaggedJets_0Flav ; 
        }
      } //// Loop over AK5 jets 

      h_nJets_bFlav -> Fill(nJets_bFlav); 
      h_nJets_cFlav -> Fill(nJets_cFlav); 
      h_nJets_lFlav -> Fill(nJets_lFlav);  
      h_nJets_0Flav -> Fill(nJets_0Flav);  

      h_nbtaggedJets_bFlav -> Fill(nbtaggedJets_bFlav); 
      h_nbtaggedJets_cFlav -> Fill(nbtaggedJets_cFlav); 
      h_nbtaggedJets_lFlav -> Fill(nbtaggedJets_lFlav);  
      h_nbtaggedJets_0Flav -> Fill(nbtaggedJets_0Flav);  

    } /// If event passes selection 

  } //// entry loop 

}

// ------------ method called once each job just after ending the event loop  ------------
void BTagEff::endJob(){ 

  h2_btagEff_ptJet_etaJet_bFlav -> Divide(h2_ptJet_etaJet_bFlav_btagged, h2_ptJet_etaJet_bFlav); 
  h2_btagEff_ptJet_etaJet_cFlav -> Divide(h2_ptJet_etaJet_cFlav_btagged, h2_ptJet_etaJet_cFlav); 
  h2_btagEff_ptJet_etaJet_lFlav -> Divide(h2_ptJet_etaJet_lFlav_btagged, h2_ptJet_etaJet_lFlav);  
  h2_btagEff_ptJet_etaJet_0Flav -> Divide(h2_ptJet_etaJet_0Flav_btagged, h2_ptJet_etaJet_0Flav);  

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
