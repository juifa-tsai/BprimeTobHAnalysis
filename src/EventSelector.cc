#include "BpbH/BprimeTobHAnalysis/interface/EventSelector.h"

EventSelector::EventSelector (const edm::ParameterSet& iConfig, EvtInfoBranches &evt, VertexInfoBranches &vtx,
    JetInfoBranches &jet, JetInfoBranches &fatjet, JetInfoBranches &subjet) : 
  evtInfo_(&evt),
  vtxInfo_(&vtx),
  jetInfo_(&jet),
  fatjetInfo_(&fatjet),
  subjetInfo_(&subjet),
  cutLevels_(iConfig.getParameter<std::vector<std::string>>("CutLevels")),
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
  jetSel() ; 
  bjetSel() ; 
  fatjetSel() ; 
  higgsjetSel() ; 
  int njets = nGoodJets() ; 
  int nbjets = nGoodBJets() ; 
  int nFatJets = nGoodFatJets() ; 
  int nHiggsJets = nGoodHiggsJets() ; 
  int level = -1;
  for (int ilevel = 0; ilevel < (int)cutLevels_.size(); ++ilevel) {
    if      ( level < (ilevel - 1) ) break ; 
    if      ( boost::iequals(cutLevels_[ilevel], "Trigger") )      { if ( trigDecision_ )           level = ilevel ; }
    else if ( boost::iequals(cutLevels_[ilevel], "Vertex") )       { if ( nGoodPV_ >= 0 )           level = ilevel ; }
    else if ( boost::iequals(cutLevels_[ilevel], "MinJets") )      { if ( njets >= minNJets_ )      level = ilevel ; }
    else if ( boost::iequals(cutLevels_[ilevel], "MaxJets") )      { if ( njets <  maxNJets_ )      level = ilevel ; }
    else if ( boost::iequals(cutLevels_[ilevel], "MinBjets") )     { if ( njets >= minNBJets_ )     level = ilevel ; }
    else if ( boost::iequals(cutLevels_[ilevel], "MaxBjets") )     { if ( njets <  maxNBJets_ )     level = ilevel ; }
    else if ( boost::iequals(cutLevels_[ilevel], "MinFatjets") )   { if ( njets >= minNFatJets_ )   level = ilevel ; }
    else if ( boost::iequals(cutLevels_[ilevel], "MaxFatjets") )   { if ( njets <  maxNFatJets_ )   level = ilevel ; }
    else if ( boost::iequals(cutLevels_[ilevel], "MinHiggsjets") ) { if ( njets >= minNHiggsJets_ ) level = ilevel ; }
    else if ( boost::iequals(cutLevels_[ilevel], "MaxHiggsjets") ) { if ( njets <  maxNHiggsJets_ ) level = ilevel ; }
    else level = ilevel ; //unknown cut, let it pass
  }

  return level ; 
}

void EventSelector::jetSel () {
  pat::strbitset retjetidak5 = jetSelector_->getBitTemplate() ;
  for (int ijet = 0; ijet < jetInfo_->Size; ++ijet) {
    if (jetSelector_->operator()(*jetInfo_, ijet,retjetidak5) == 0) continue ;
    goodJets_.push_back(ijet) ; 
  }

}

void EventSelector::bjetSel () {
  pat::strbitset retjetid = bjetSelector_->getBitTemplate() ;
  for (int ijet = 0; ijet < jetInfo_->Size; ++ijet) {
    if (bjetSelector_->operator()(*jetInfo_, ijet,retjetid) == 0) continue ;
    goodBJets_.push_back(ijet) ; 
  }

}

void EventSelector::fatjetSel () {
  pat::strbitset retjetidca8 = fatjetSelector_->getBitTemplate() ;
  for (int ijet = 0; ijet < fatjetInfo_->Size; ++ijet) {
    if (fatjetSelector_->operator()(*fatjetInfo_, ijet, *subjetInfo_, retjetidca8) == 0) continue ;
    goodFatJets_.push_back(ijet) ; 
  }

}

void EventSelector::higgsjetSel () {
  pat::strbitset retjetidca8 = higgsjetSelector_->getBitTemplate() ;
  for (int ijet = 0; ijet < fatjetInfo_->Size; ++ijet) {
    if (higgsjetSelector_->operator()(*fatjetInfo_, ijet, *subjetInfo_, retjetidca8) == 0) continue ;
    goodHiggsJets_.push_back(ijet) ; 
  }

}

void EventSelector::printEventStatus() {
  edm::LogInfo("EventSelector") << "Run " << evtInfo_->RunNo << " event " << evtInfo_->EvtNo << " status:\n";
  edm::LogInfo("EventSelector") << "\tPrimary vertex " << primaryVertex();

  edm::LogInfo("EventSelector") << "\n\t" << nGoodJets() << " good jets:";
  for(unsigned i=0; i< goodJets_.size(); i++)
    edm::LogInfo("EventSelector") << " " <<  goodJets_[i];

  edm::LogInfo("EventSelector") << "\n\tPass level: " << (passLevel_ > -1 && cutLevels_.size() > 0 ? cutLevels_[passLevel_] : "none") ; 
}

void EventSelector::reset() {
  nGoodPV_ = -1;
  trigDecision_ = 0 ; 
  goodJets_.clear();
  goodBJets_.clear();
  goodFatJets_.clear();
  goodHiggsJets_.clear();

  passLevel_ = -1;
}

std::vector<std::string> EventSelector::getCutLevels() {
  return cutLevels_;
}



