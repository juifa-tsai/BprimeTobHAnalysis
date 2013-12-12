#include "BpbH/BprimeTobHAnalysis/interface/EventSelector.h"

EventSelector::EventSelector (const edm::ParameterSet& iConfig, EvtInfoBranches &evt, VertexInfoBranches &vtx,
    JetInfoBranches &jet, JetInfoBranches &fatjet, JetInfoBranches &subjet) : 
  evtInfo_(&evt),
  vtxInfo_(&vtx),
  jetInfo_(&jet),
  fatjetInfo_(&fatjet),
  subjetInfo_(&subjet),
  cutLevels_(iConfig.getParameter<std::vector<std::string>>("CutLevels")),
  removeJetHiggsJetOverlap_(iConfig.getParameter<bool>("RemoveJetHiggsJetOverlap")), 
  minNJets_(iConfig.getParameter<int>("MinNJets")),  
  maxNJets_(iConfig.getParameter<int>("MaxNJets")),  
  minNBJets_(iConfig.getParameter<int>("MinNBJets")),  
  maxNBJets_(iConfig.getParameter<int>("MaxNBJets")),  
  minNFatJets_(iConfig.getParameter<int>("MinNFatJets")),  
  maxNFatJets_(iConfig.getParameter<int>("MaxNFatJets")),  
  minNHiggsJets_(iConfig.getParameter<int>("MinNHiggsJets")),  
  maxNHiggsJets_(iConfig.getParameter<int>("MaxNHiggsJets")),  
  passes_(0) 
{

  trigSelector_    = new TriggerSelector(iConfig.getParameter<edm::ParameterSet>("TrigSelParams")) ; 
  vtxSelector_     = new VertexSelector(vtx) ;
  jetSelector_     = new JetSelector(iConfig.getParameter<edm::ParameterSet>("JetSelParams")) ; 
  bjetSelector_    = new JetSelector(iConfig.getParameter<edm::ParameterSet>("BJetSelParams")) ; 
  fatjetSelector_  = new FatJetSelector(iConfig.getParameter<edm::ParameterSet>("FatJetSelParams")) ; 
  higgsjetSelector_= new FatJetSelector(iConfig.getParameter<edm::ParameterSet>("HiggsJetSelParams")) ; 
  htSelector_      = new HTSelector(iConfig.getParameter<edm::ParameterSet>("HTSelParams")) ;  

}

EventSelector::~EventSelector () {}

bool EventSelector::passesTrigger () {
  return trigDecision_; 
}

int EventSelector::primaryVertex () {
  return nGoodPV_ ; 
}

int EventSelector::nGoodJets () {
  return goodJets_.size() ; 
}

int EventSelector::nGoodBJets () {
  return goodBJets_.size() ; 
}

int EventSelector::nGoodFatJets () {
  return goodFatJets_.size() ; 
}

int EventSelector::nGoodHiggsJets () {
  return goodHiggsJets_.size() ; 
}

JetCollection EventSelector::goodJets () {
  return goodJets_ ; 
}

JetCollection EventSelector::goodBJets () {
  return goodBJets_ ; 
}

JetCollection EventSelector::goodFatJets () {
  return goodFatJets_ ; 
}

JetCollection EventSelector::goodHiggsJets () {
  return goodHiggsJets_ ; 
}

bool EventSelector::passes() {
  bool passes = false;
  int code = passCode();
  if(code >= (int)cutLevels_.size()-1) passes=true;
  return passes;
}

int EventSelector::passCode() {
  if(cutLevels_.size()==0) return 0;

  trigDecision_ = trigSelector_->getTrigDecision(*evtInfo_) ; 
  nGoodPV_ = vtxSelector_->NGoodVtxs(); 
  fatjetSel() ; 
  higgsjetSel() ; 
  jetSel() ; 
  bjetSel() ; 
  int njets = nGoodJets() ; 
  int nbjets = nGoodBJets() ; 
  int nFatJets = nGoodFatJets() ; 
  int nHiggsJets = nGoodHiggsJets() ; 
  bool passesHT = HTSel() ; 
  int level = -1;
  for (int ilevel = 0; ilevel < (int)cutLevels_.size(); ++ilevel) {
    if      ( level < (ilevel - 1) ) break ; 
    if      ( boost::iequals(cutLevels_[ilevel], "Trigger") )      { if ( trigDecision_ )           level = ilevel ; }
    else if ( boost::iequals(cutLevels_[ilevel], "Vertex") )       { if ( nGoodPV_ >= 0 )           level = ilevel ; }
    else if ( boost::iequals(cutLevels_[ilevel], "MinFatjets") )   { if ( njets >= minNFatJets_ )   level = ilevel ; }
    else if ( boost::iequals(cutLevels_[ilevel], "MaxFatjets") )   { if ( njets <  maxNFatJets_ )   level = ilevel ; }
    else if ( boost::iequals(cutLevels_[ilevel], "MinHiggsjets") ) { if ( njets >= minNHiggsJets_ ) level = ilevel ; }
    else if ( boost::iequals(cutLevels_[ilevel], "MaxHiggsjets") ) { if ( njets <  maxNHiggsJets_ ) level = ilevel ; }
    else if ( boost::iequals(cutLevels_[ilevel], "MinJets") )      { if ( njets >= minNJets_ )      level = ilevel ; }
    else if ( boost::iequals(cutLevels_[ilevel], "MaxJets") )      { if ( njets <  maxNJets_ )      level = ilevel ; }
    else if ( boost::iequals(cutLevels_[ilevel], "MinBjets") )     { if ( njets >= minNBJets_ )     level = ilevel ; }
    else if ( boost::iequals(cutLevels_[ilevel], "MaxBjets") )     { if ( njets <  maxNBJets_ )     level = ilevel ; }
    else if ( boost::iequals(cutLevels_[ilevel], "HT") )           { if ( passesHT )                level = ilevel ; }
    else level = ilevel ; //unknown cut, let it pass
  }

  return level ; 
}

void EventSelector::jetSel () {
  pat::strbitset retjetidak5 = jetSelector_->getBitTemplate() ;
  for (int ijet = 0; ijet < jetInfo_->Size; ++ijet) {
    if (jetSelector_->operator()(*jetInfo_, ijet,retjetidak5) == 0) continue ;
    Jet thisjet(*jetInfo_, ijet) ; 
    if (removeJetHiggsJetOverlap_ && hasOverlap(thisjet, goodHiggsJets_)) continue ; 
    goodJets_.push_back(thisjet) ; 
  }

}

bool EventSelector::hasOverlap(Jet& jet, JetCollection& jetcoll) {
  bool isJetOverlap(false) ; 
  TLorentzVector jet_p4;
  jet_p4.SetPtEtaPhiM(jet.Pt(), jet.Eta(), jet.Phi(), jet.Mass()); 
    for (JetCollection::const_iterator ifatjet = jetcoll.begin(); ifatjet != jetcoll.end(); ++ifatjet) {
      TLorentzVector fatjet_p4 ; 
      fatjet_p4.SetPtEtaPhiM(ifatjet->Pt(), ifatjet->Eta(), ifatjet->Phi(), ifatjet->Mass());  
      if (jet_p4.DeltaR(fatjet_p4) < 1.2) isJetOverlap = true ; 
      break ; 
    }
  return isJetOverlap ; 
}

void EventSelector::bjetSel () {
  pat::strbitset retjetid = bjetSelector_->getBitTemplate() ;
  for (int ijet = 0; ijet < jetInfo_->Size; ++ijet) {
    if (bjetSelector_->operator()(*jetInfo_, ijet,retjetid) == 0) continue ;
    Jet thisjet(*jetInfo_, ijet) ; 
    if (removeJetHiggsJetOverlap_ && hasOverlap(thisjet, goodHiggsJets_)) continue ; 
    goodBJets_.push_back(thisjet) ; 
  }

}

void EventSelector::fatjetSel () {
  pat::strbitset retjetidca8 = fatjetSelector_->getBitTemplate() ;
  for (int ijet = 0; ijet < fatjetInfo_->Size; ++ijet) {
    if (fatjetSelector_->operator()(*fatjetInfo_, ijet, *subjetInfo_, retjetidca8) == 0) continue ;
    Jet thisfatjet(*fatjetInfo_, ijet) ; 
    goodFatJets_.push_back(thisfatjet) ; 
  }

}

void EventSelector::higgsjetSel () {
  pat::strbitset retjetidca8 = higgsjetSelector_->getBitTemplate() ;
  for (int ijet = 0; ijet < fatjetInfo_->Size; ++ijet) {
    if (higgsjetSelector_->operator()(*fatjetInfo_, ijet, *subjetInfo_, retjetidca8) == 0) continue ;
    Jet thishiggsjet(*fatjetInfo_, ijet) ; 
    goodHiggsJets_.push_back(thishiggsjet) ; 
  }

}

bool EventSelector::HTSel () {

  HT_.clearJetCollections() ; 
  HTVal_ = 0 ; 
  HT_.setJetCollection(goodJets_) ; 
  HT_.buildHT() ; 
  HTVal_ = HT_.getHT() ; 
  pat::strbitset retht = htSelector_->getBitTemplate() ; 
  retht.set(false) ; 
  return htSelector_->operator()(HT_, retht) ; 

}

void EventSelector::printEventStatus() {
  edm::LogInfo("EventSelector") << "Run " << evtInfo_->RunNo << " event " << evtInfo_->EvtNo << " status:\n";
  edm::LogInfo("EventSelector") << "\tPrimary vertex " << primaryVertex();

  edm::LogInfo("EventSelector") << "\n\t" << nGoodJets() << " good jets:";

  edm::LogInfo("EventSelector") << "\n\tPass level: " << (passLevel_ > -1 && cutLevels_.size() > 0 ? cutLevels_[passLevel_] : "none") ; 
}

void EventSelector::reset() {
  nGoodPV_ = -1;
  trigDecision_ = 0 ; 
  goodJets_.clear();
  goodBJets_.clear();
  goodFatJets_.clear();
  goodHiggsJets_.clear();
  HT_.clearJetCollections() ; 
  HTVal_ = 0 ; 

  passLevel_ = -1;
}

std::vector<std::string> EventSelector::getCutLevels() {
  return cutLevels_;
}



