#include "BpbH/BprimeTobHAnalysis/interface/JMEUncertUtil.h"
#include <fstream> 
#include <boost/algorithm/string.hpp>

JMEUncertUtil::JMEUncertUtil (const edm::ParameterSet& iConfig, EvtInfoBranches &evt, JetCollection &jets, std::string stype, double jecShift) : 
  jerEta_(iConfig.getUntrackedParameter< std::vector<double> >("JEREta")), 
  jerNominal_(iConfig.getUntrackedParameter<std::vector<double> >("JERNominal")),
  jerSigmaSym_(iConfig.getUntrackedParameter<std::vector<double> >("JERSigmaSym")),
  jerSigmaNeg_(iConfig.getUntrackedParameter<std::vector<double> >("JERSigmaNeg")),
  jerSigmaPos_(iConfig.getUntrackedParameter<std::vector<double> >("JERSigmaPos")), 
  stype_(stype),
  jecShift_(jecShift),
  evt_(&evt),
  jets_(jets) 
{

  setType(stype_); 

  if (jecType_ == NOSYS) {
    edm::LogWarning("JMEUncertUtil::JMEUncertUtil") << " jecType_ = " << jecType_ << ". Returning same jet collection."; 
    for (JetCollection::const_iterator ijet = jets_.begin(); ijet != jets_.end(); ++ijet) {
      modifiedJets_.push_back(*ijet) ; 
    }  
    return ;  
  }

  if (jecType_ == JES) {
    LogDebug("JMEUncertUtil::JMEUncertUtil") << " jecType_ = " << jecType_ << ". Doing " << stype_ << " " << jecShift_ ; 
    std::string filenameJEC = iConfig.getUntrackedParameter<std::string>("FilenameJEC") ;
    ifstream fileJEC(filenameJEC.c_str());
    if ( fileJEC ) {
      jecUncert_ = new JetCorrectionUncertainty(*(new JetCorrectorParameters(filenameJEC, "Total"))); 
    }
    else {
      edm::LogError("JMEUncertUtil") << "ERROR: Couldn't open JES File, can't continue." ; 
      assert(false);
    } 
    jesUncert () ; 
  }
  else if (jecType_ == JER) {
    LogDebug("JMEUncertUtil::JMEUncertUtil") << " jecType_ = " << jecType_ << ". Doing " << stype_ << " " << jecShift_ ; 
    jerScale() ; 
  }

  return ; 

} 

JMEUncertUtil::~JMEUncertUtil () {
  modifiedJets_.clear(); 
} 

JetCollection JMEUncertUtil::GetModifiedJetColl () const {
  return modifiedJets_ ; 
}

const string JMEUncertUtil::JEC_Types [] = {"None", "JES", "JER"} ;  

void JMEUncertUtil::setType(std::string type) {
  jecType_ = NOSYS;

  for(JEC_TYPE ii = NOSYS; ii < NTYPES; ii = static_cast<JEC_TYPE>(ii+1)) { 
    if(boost::iequals(type, JEC_Types[ii])) jecType_ = ii ;  
  }
  return ; 
}

void JMEUncertUtil::jesUncert () {

  if( jecShift_ == 0.) { 
    for (JetCollection::const_iterator ijet = jets_.begin(); ijet != jets_.end(); ++ijet) {
      modifiedJets_.push_back(*ijet) ; 
    }  
    return ; 
  }

  for (JetCollection::const_iterator ijet = jets_.begin(); ijet != jets_.end(); ++ijet) {
    jecUncert_ -> setJetPt(ijet -> Pt());
    jecUncert_ -> setJetEta(ijet -> Eta()); 
    double uncert = jecUncert_ -> getUncertainty ( jecShift_ > 0. ); 
    double rescaled = ( jecShift_ > 0. ) ? ( 1. + uncert ) : ( 1. - uncert ) ;
    LogDebug("JMEUncertUtil::jesUncert") << stype_ << " " << jecShift_ << " rescaled = " << rescaled  ; 
    Jet thisjet(*ijet) ;
    thisjet.Set_Pt        ( thisjet.Pt()        * rescaled ) ; 
    thisjet.Set_Et        ( thisjet.Et()        * rescaled ) ; 
    thisjet.Set_PtCorrRaw ( thisjet.PtCorrRaw() * rescaled ) ; 
    thisjet.Set_PtCorrL2  ( thisjet.PtCorrL2 () * rescaled ) ; 
    thisjet.Set_PtCorrL7g ( thisjet.PtCorrL7g() * rescaled ) ; 
    thisjet.Set_PtCorrL7c ( thisjet.PtCorrL7c() * rescaled ) ; 
    thisjet.Set_Px        ( thisjet.Px       () * rescaled ) ; 
    thisjet.Set_Py        ( thisjet.Py       () * rescaled ) ; 
    thisjet.Set_Pz        ( thisjet.Pz       () * rescaled ) ; 
    thisjet.Set_Energy    ( thisjet.Energy   () * rescaled ) ; 
    modifiedJets_.push_back(thisjet) ; 
  }

  return ; 

}

void JMEUncertUtil::jerScale() {

  if ( !evt_ -> McFlag ) {
    for (JetCollection::const_iterator ijet = jets_.begin(); ijet != jets_.end(); ++ijet) {
      modifiedJets_.push_back(*ijet) ; 
    }  
    return ; 
  }

  if( jecShift_ == 0.) { 
    for (JetCollection::const_iterator ijet = jets_.begin(); ijet != jets_.end(); ++ijet) {
      modifiedJets_.push_back(*ijet) ; 
    }  
    return ; 
  }

  //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution

  for (JetCollection::const_iterator ijet = jets_.begin(); ijet != jets_.end(); ++ijet) {
    double rescaled(0.) ; 
    if( ijet -> Pt() <= 0.) rescaled = 1. ; // Not sure how this could happen, but just in case
    if( ijet -> GenJetPt () < 15. ) rescaled = 1.; //// Attn: hard-coded 
    int etaB = -1; 
    for(int ieta = 0; ieta < (int)jerEta_.size(); ++ieta) {
      if( etaB >= 0 ) break ;
      if( fabs( ijet -> GenJetEta() ) < jerEta_.at(ieta) ) etaB = ieta ;  
    }
    if ( etaB < 0) etaB = (int)jerEta_.size(); //must be forward 
    double sf = jerNominal_[etaB]; 
    if( jecShift_ < 0.) sf -= sqrt( pow(jerSigmaSym_[etaB], 2) + pow (jerSigmaNeg_[etaB], 2) ) ; 
    else if ( jecShift_ > 0. ) sf += sqrt( pow(jerSigmaSym_[etaB], 2)+pow(jerSigmaPos_[etaB], 2) ) ; 
    LogDebug("JMEUncertUtil::jerScale") << stype_ << " " << jecShift_ << " sf = " << sf ;
    double deltaPt = ( ijet -> Pt() - ijet ->GenJetPt() ) * sf ; 
    rescaled = std::max(0.0, ( ijet ->GenJetPt() + deltaPt ))/ijet -> Pt() ; 
    if ( rescaled <= 0. ) rescaled =  1. ; 
    Jet thisjet(*ijet) ;
    thisjet.Set_Pt        ( thisjet.Pt()        * rescaled ) ; 
    thisjet.Set_Et        ( thisjet.Et()        * rescaled ) ; 
    thisjet.Set_PtCorrRaw ( thisjet.PtCorrRaw() * rescaled ) ; 
    thisjet.Set_PtCorrL2  ( thisjet.PtCorrL2 () * rescaled ) ; 
    thisjet.Set_PtCorrL7g ( thisjet.PtCorrL7g() * rescaled ) ; 
    thisjet.Set_PtCorrL7c ( thisjet.PtCorrL7c() * rescaled ) ; 
    thisjet.Set_Px        ( thisjet.Px       () * rescaled ) ; 
    thisjet.Set_Py        ( thisjet.Py       () * rescaled ) ; 
    thisjet.Set_Pz        ( thisjet.Pz       () * rescaled ) ; 
    thisjet.Set_Energy    ( thisjet.Energy   () * rescaled ) ; 
    modifiedJets_.push_back(thisjet) ; 
  }

  return ; 

}
