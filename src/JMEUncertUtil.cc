#include "BpbH/BprimeTobHAnalysis/interface/JMEUncertUtil.h"
#include <fstream> 
#include <boost/algorithm/string.hpp>

JMEUncertUtil::JMEUncertUtil (const edm::ParameterSet& iConfig, JetCollection &jets, JetCollection &genjets, std::string stype, double jecShift) : 
  jerEta_(iConfig.getUntrackedParameter< std::vector<double> >("JEREta")), 
  jerNominal_(iConfig.getUntrackedParameter<std::vector<double> >("JERNominal")),
  jerSigmaSym_(iConfig.getUntrackedParameter<std::vector<double> >("JERSigmaSym")),
  jerSigmaNeg_(iConfig.getUntrackedParameter<std::vector<double> >("JERSigmaNeg")),
  jerSigmaPos_(iConfig.getUntrackedParameter<std::vector<double> >("JERSigmaPos")), 
  stype_(stype),
  jecShift_(jecShift),
  jets_(jets), 
  genjets_(genjets)
{

  modifiedJets_.clear(); 

  setType(stype_); 

  if (jecType_ == NOSYS) {
    edm::LogWarning("JMEUncertUtil::JMEUncertUtil") << " jecType_ = " << jecType_ << " not correct. Returning same jet collection."; 
    return ;  
  }
  else if (jecType_ == JESAK5MC || jecType_ == JESCA8MC) {
    LogDebug("JMEUncertUtil::JMEUncertUtil") << " jecType_ = " << jecType_ << ". Doing " << stype_ << " " << jecShift_ ; 
    std::string filenameJEC ; 
    if (jecType_ == JESAK5MC) filenameJEC = iConfig.getUntrackedParameter<std::string>("FilenameJECAK5MC") ;
    else filenameJEC = iConfig.getUntrackedParameter<std::string>("FilenameJECCA8MC") ;
    ifstream fileJEC(filenameJEC.c_str());
    if ( fileJEC ) {
      jecUncert_ = new JetCorrectionUncertainty(*(new JetCorrectorParameters(filenameJEC, "Total"))); 
      jesUncert () ; 
    }
    else {
      edm::LogError("JMEUncertUtil") << "ERROR: Couldn't open JES File, can't continue." ; 
      assert(false);
    } 
  }
  else if (jecType_ == JERAK5MC || jecType_ == JERCA8MC) {
    LogDebug("JMEUncertUtil::JMEUncertUtil") << " jecType_ = " << jecType_ << ". Doing " << stype_ << " " << jecShift_ ; 
    jerScale() ; 
  }
  else {
    edm::LogWarning("JMEUncertUtil::JMEUncertUtil") << " jecType_ = " << jecType_ << " not recognised. Returning same jet collection."; 
    return ;  
  }

  return ; 

} 

JMEUncertUtil::~JMEUncertUtil () {
  modifiedJets_.clear(); 
} 

JetCollection JMEUncertUtil::GetModifiedJetColl () const {
  return modifiedJets_ ; 
}

const string JMEUncertUtil::JEC_Types [] = {"NONE", "JESAK5MC", "JESCA8MC", "JERAK5MC", "JERCA8MC"} ; 

void JMEUncertUtil::setType(std::string type) {
  jecType_ = NOSYS;

  for(JEC_TYPE ii = NOSYS; ii < NTYPES; ii = static_cast<JEC_TYPE>(ii+1)) { 
    if(boost::iequals(type, JEC_Types[ii])) jecType_ = ii ;  
  }
  return ; 
}

void JMEUncertUtil::jesUncert () {

  for (JetCollection::const_iterator ijet = jets_.begin(); ijet != jets_.end(); ++ijet) {
    jecUncert_ -> setJetPt(ijet -> Pt());
    jecUncert_ -> setJetEta(ijet -> Eta()); 
    double uncert = jecUncert_ -> getUncertainty ( jecShift_ > 0. ); 
    double rescaled(1.) ;
    if ( jecShift_ > 0. ) rescaled = 1. + uncert ;
    else if ( jecShift_ < 0. ) rescaled = 1. - uncert  ; 
    else rescaled = 1. ; 
    LogDebug("JMEUncertUtil::jesUncert") << stype_ << " " << jecShift_ << " rescaled = " << rescaled  ; 
    Jet thisjet(*ijet) ;
    LogDebug("JMEUncertUtil::jesUncert") << " JES:jet pt before = " << thisjet.Pt() << " after = " << thisjet.Pt()*rescaled << " jet eta = " << thisjet.Eta() << std::endl ; 
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

  //https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution

  edm::LogInfo("JMEUncertUtil::jerScale")  << " reco jets size = " << jets_.size() << " gen jets size = " << genjets_.size() ; 

  if (jecType_ == JERCA8MC) {
    for (JetCollection::const_iterator ijet = jets_.begin(); ijet != jets_.end(); ++ijet) {
      Jet thisjet(*ijet) ;
      modifiedJets_.push_back(thisjet) ; 
    }
  }
  else {
    for (JetCollection::const_iterator ijet = jets_.begin(); ijet != jets_.end(); ++ijet) {
      //// Gen GenJet matchted to this jet 
      double rescaled(0.) ; 
      const Jet* matchedGenJet = NULL ; 
      for (JetCollection::const_iterator igenjet = genjets_.begin(); igenjet != genjets_.end(); ++ijet) {
        if ( ijet->DeltaR(*igenjet) < 0.5 ) { matchedGenJet = &(*igenjet) ; } 
        else edm::LogInfo("JMEUncertUtil::jerScale") << " gen jet pt " << igenjet->Pt() << " reco jet pt " << ijet->Pt() << " DR = " << ijet->DeltaR(*igenjet); 
      }
      if( matchedGenJet && matchedGenJet -> Pt() < 1E-6) rescaled = 1. ; 
      if( matchedGenJet && matchedGenJet -> Pt () < 15. ) rescaled = 1.; //// Attn: hard-coded 
      int etaB = -1; 
      for(int ieta = 0; ieta < (int)jerEta_.size(); ++ieta) {
        if( etaB >= 0 ) break ;
        if( fabs( ijet -> Eta() ) < jerEta_.at(ieta) ) etaB = ieta ;  
      }
      if ( etaB < 0) etaB = (int)jerEta_.size(); //must be forward 
      double sf = jerNominal_[etaB]; 
      if( jecShift_ < 0.) sf -= sqrt( pow(jerSigmaSym_[etaB], 2) + pow (jerSigmaNeg_[etaB], 2) ) ; 
      else if ( jecShift_ > 0. ) sf += sqrt( pow(jerSigmaSym_[etaB], 2)+pow(jerSigmaPos_[etaB], 2) ) ; 
      edm::LogInfo("JMEUncertUtil::jerScale") << stype_ << " " << jecShift_ << " sf = " << sf ;
      double deltaPt(0) ;
      if ( matchedGenJet ) deltaPt = ( ijet -> Pt() - matchedGenJet -> Pt() ) * sf ; 
      if ( matchedGenJet ) rescaled = std::max(0.0, ( matchedGenJet -> Pt() + deltaPt ))/ijet -> Pt() ; 
      if ( rescaled <= 0. ) rescaled =  1. ; 
      Jet thisjet(*ijet) ;
      if ( matchedGenJet ) edm::LogInfo("JMEUncertUtil::jerScale") << " Gen jet pt = " << matchedGenJet -> Pt() << " reco jet pT = " << ijet -> Pt() << " sf = " << sf << " deltaPt = " << deltaPt ; 
      else edm::LogInfo("JMEUncertUtil::jerScale") << " matched gen jet not found. " << " reco jet pT = " << ijet -> Pt() << " sf = " << sf << " deltaPt = " << deltaPt ;
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
  }

  return ; 

}
