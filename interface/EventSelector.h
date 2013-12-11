#ifndef BPRIMETOBHANALYSIS_INTERFACE_EVENTSELECTOR_H
#define BPRIMETOBHANALYSIS_INTERFACE_EVENTSELECTOR_H

#include <memory>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h" 

#include "BpbH/BprimeTobH/interface/TriggerSelector.h"
#include "BpbH/BprimeTobH/interface/VertexSelector.h"
#include "BpbH/BprimeTobH/interface/JetSelector.h"
#include "BpbH/BprimeTobH/interface/FatJetSelector.h"
#include "BpbH/BprimeTobH/interface/HTSelector.h"

class EventSelector {

  public : 

    enum SELECTYPES_t {EVTSEL, VTXSEL, JETSEL, BJETSEL, FATJETSEL, HIGGSJETSEL, HTSEL, LEPTONVETO} ; 

    EventSelector (edm::ParameterSet const& params, 
        EvtInfoBranches &evt, VertexInfoBranches &vtx, 
        JetInfoBranches &jet, JetInfoBranches &fatjet, JetInfoBranches &subjet) ; 
    ~EventSelector() ; 

    void printEventStatus();
    int passCode(); 
    void reset();
    std::vector<std::string> getCutLevels();

    bool passes();
    bool passesTrigger();
    int primaryVertex();
    int nGoodFatJets();
    int nGoodHiggsJets(); 
    int nGoodJets();
    int nGoodBJets();

    JetCollection goodJets();
    JetCollection goodBJets();
    JetCollection goodFatJets();
    JetCollection goodHiggsJets();

  private: 

    void fatjetSel() ; 
    void higgsjetSel() ; 
    void jetSel() ; 
    void bjetSel() ; 
    bool HTSel() ; 
    bool hasOverlap(Jet&, JetCollection&) ; 

    //// Configurables 
    EvtInfoBranches*    evtInfo_ ; 
    VertexInfoBranches* vtxInfo_;
    JetInfoBranches*    jetInfo_;
    JetInfoBranches*    fatjetInfo_; 
    JetInfoBranches*    subjetInfo_; 

    TriggerSelector* trigSelector_ ;
    VertexSelector*  vtxSelector_ ;
    JetSelector*     jetSelector_ ;
    JetSelector*     bjetSelector_ ;
    FatJetSelector*  fatjetSelector_ ; 
    FatJetSelector*  higgsjetSelector_ ; 
    HTSelector*      htSelector_ ; 

    std::vector<std::string> cutLevels_ ; 
    bool removeJetHiggsJetOverlap_ ; 

    int nGoodPV_ ;
    bool trigDecision_ ; 
    JetCollection goodJets_;
    JetCollection goodBJets_;
    JetCollection goodFatJets_;
    JetCollection goodHiggsJets_;
    HT HT_ ; 

    int    minNJets_ ; 
    int    maxNJets_ ; 
    int    minNBJets_ ; 
    int    maxNBJets_ ; 
    int    minNFatJets_ ; 
    int    maxNFatJets_ ; 
    int    minNHiggsJets_ ; 
    int    maxNHiggsJets_ ; 
    double HTVal_ ; 

    int passLevel_ ; 
    bool passes_ ; 

} ; 

#endif 
