#include "BpbH/BprimeTobHAnalysis/interface/ApplyBTagSF.h"
#include "BpbH/BprimeTobHAnalysis/interface/SFb-pt_WITHttbar_payload_EPS13.h"
#include "BpbH/BprimeTobHAnalysis/interface/SFlightFuncs_EPS2013.h"

ApplyBTagSF::ApplyBTagSF(JetCollection& jets, double bDisc, std::string algo_, std::string mode_) : 
  jets_(jets),
  bDisc_(bDisc),
  algo_(algo_),
  mode_(mode_) 
{
  btaggedJetsWithSF_.clear() ; 
}

ApplyBTagSF::~ApplyBTagSF() {
  btaggedJetsWithSF_.clear() ; 
  algo_.clear();
  mode_.clear(); 
}

JetCollection ApplyBTagSF::getBtaggedJetsWithSF () { 

  int    jet_flavor;
  double jet_et;       
  double jet_phi;        
  double jet_eta;
  bool   isTagged = false;

  double BTagSF = 0.972;
  double BTageff = 0.705;
  double LightJetSF = 1.04;
  double LightJeteff = 0.01705;

  for (JetCollection::const_iterator ijet = jets_.begin(); ijet != jets_.end(); ++ijet) {

    jet_flavor = ijet->GenFlavor() ; 
    jet_et     = ijet->Et()        ; 
    jet_phi    = ijet->Eta()       ; 
    jet_eta    = ijet->Phi()       ; 
    isTagged   = ijet->CombinedSVBJetTags() > bDisc_ ; 

    int ptbin(-1) ; 
    for (int ii = 0; ii < int(sizeof(ptmin)/sizeof(ptmin[0])); ++ii) {  
      if ( jet_et >= ptmin[ii] && jet_et < ptmax[ii] ) {
        ptbin = ii ; 
      }
    }

    double errscale(1.);
    if (abs(jet_flavor) == 4) errscale *= 2. ; 

    if (strcmp(algo_.c_str(), "CSVL") == 0) {
      if (strcmp(mode_.c_str(), "Mean") == 0)            BTagSF = SFb_CSVL->Eval(jet_et) ; 
      else if (strcmp(mode_.c_str(), "1sigmaUp") == 0)   BTagSF = SFb_CSVL->Eval(jet_et)*( 1 - (errscale*SFb_CSVL_error[ptbin]) ) ; 
      else if (strcmp(mode_.c_str(), "1sigmaDown") == 0) BTagSF = SFb_CSVL->Eval(jet_et)*( 1 + (errscale*SFb_CSVL_error[ptbin]) ) ; 
      else {
        edm::LogError("ApplyBTagSF") << " Wrong b-tagging mode_ chosen: " << mode_  << ". Choose between Mean, 1sigmaUp, or 1sigmaDown." ; 
        return btaggedJetsWithSF_ ; 
      }
    }
    else if (strcmp(algo_.c_str(), "CSVM") == 0) {
      if (strcmp(mode_.c_str(), "Mean") == 0)            BTagSF = SFb_CSVL->Eval(jet_et) ; 
      else if (strcmp(mode_.c_str(), "1sigmaUp") == 0)   BTagSF = SFb_CSVL->Eval(jet_et)*( 1 - (errscale*SFb_CSVM_error[ptbin]) ) ; 
      else if (strcmp(mode_.c_str(), "1sigmaDown") == 0) BTagSF = SFb_CSVL->Eval(jet_et)*( 1 + (errscale*SFb_CSVM_error[ptbin]) ) ; 
      else {
        edm::LogError("ApplyBTagSF") << " Wrong b-tagging mode_ chosen: " << mode_  << ". Choose between Mean, 1sigmaUp, or 1sigmaDown." ; 
        return btaggedJetsWithSF_ ; 
      }
    }
    else if (strcmp(algo_.c_str(), "CSVT") == 0) {
      if (strcmp(mode_.c_str(), "Mean") == 0)            BTagSF = SFb_CSVL->Eval(jet_et) ; 
      else if (strcmp(mode_.c_str(), "1sigmaUp") == 0)   BTagSF = SFb_CSVL->Eval(jet_et)*( 1 - (errscale*SFb_CSVT_error[ptbin]) ) ; 
      else if (strcmp(mode_.c_str(), "1sigmaDown") == 0) BTagSF = SFb_CSVL->Eval(jet_et)*( 1 + (errscale*SFb_CSVT_error[ptbin]) ) ; 
      else {
        edm::LogError("ApplyBTagSF") << " Wrong b-tagging mode_ chosen: " << mode_  << ". Choose between Mean, 1sigmaUp, or 1sigmaDown." ; 
        return btaggedJetsWithSF_ ; 
      }
    }
    else {
      edm::LogError("ApplyBTagSF") << " Wrong b-tagging algo_ chosen: " << algo_ << ". Choose between CSVL, CSVM, or CSVT." ; 
      return btaggedJetsWithSF_ ; 
    }

    if (strcmp(algo_.c_str(), "CSVL") == 0) {
      if (strcmp(mode_.c_str(), "Mean") == 0) {
        for (int ii = 0; ii < 4; ++ii) {
          if ( jet_eta >= SFlight_CSVL_etamin[ii] && jet_eta < SFlight_CSVL_etamax[ii] ) {
            LightJetSF = (GetSFlmean("CSV","L",SFlight_CSVL_etamin[ii], SFlight_CSVL_etamax[ii], "ABCD"))->Eval(jet_et) ; 
          }
        }
      }
      else if (strcmp(mode_.c_str(), "1sigmaUp") == 0) {
        for (int ii = 0; ii < 4; ++ii) {
          if ( jet_eta >= SFlight_CSVL_etamin[ii] && jet_eta < SFlight_CSVL_etamax[ii] ) {
            LightJetSF = (GetSFlmax("CSV","L",SFlight_CSVL_etamin[ii], SFlight_CSVL_etamax[ii], "ABCD"))->Eval(jet_et) ; 
          }
        }
      }
      else if (strcmp(mode_.c_str(), "1sigmaDown") == 0) {
        for (int ii = 0; ii < 4; ++ii) {
          if ( jet_eta >= SFlight_CSVL_etamin[ii] && jet_eta < SFlight_CSVL_etamax[ii] ) {
            LightJetSF = (GetSFlmin("CSV","L",SFlight_CSVL_etamin[ii], SFlight_CSVL_etamax[ii], "ABCD"))->Eval(jet_et) ; 
          }
        }
      }
      else {
        edm::LogError("ApplyBTagSF") << " Wrong b-tagging mode_ chosen: " << mode_  << ". Choose between Mean, 1sigmaUp, or 1sigmaDown." ; 
        return btaggedJetsWithSF_ ; 
      }
    }
    else if (strcmp(algo_.c_str(), "CSVM") == 0) {
      if (strcmp(mode_.c_str(), "Mean") == 0) {
        for (int ii = 0; ii < 4; ++ii) {
          if ( jet_eta >= SFlight_CSVM_etamin[ii] && jet_eta < SFlight_CSVM_etamax[ii] ) {
            LightJetSF = (GetSFlmean("CSV","M",SFlight_CSVM_etamin[ii], SFlight_CSVM_etamax[ii], "ABCD"))->Eval(jet_et) ; 
          }
        }
      }
      else if (strcmp(mode_.c_str(), "1sigmaUp") == 0) {
        for (int ii = 0; ii < 4; ++ii) {
          if ( jet_eta >= SFlight_CSVM_etamin[ii] && jet_eta < SFlight_CSVM_etamax[ii] ) {
            LightJetSF = (GetSFlmax("CSV","M",SFlight_CSVM_etamin[ii], SFlight_CSVM_etamax[ii], "ABCD"))->Eval(jet_et) ; 
          }
        }
      }
      else if (strcmp(mode_.c_str(), "1sigmaDown") == 0) {
        for (int ii = 0; ii < 4; ++ii) {
          if ( jet_eta >= SFlight_CSVM_etamin[ii] && jet_eta < SFlight_CSVM_etamax[ii] ) {
            LightJetSF = (GetSFlmin("CSV","M",SFlight_CSVM_etamin[ii], SFlight_CSVM_etamax[ii], "ABCD"))->Eval(jet_et) ; 
          }
        }
      }
      else {
        edm::LogError("ApplyBTagSF") << " Wrong b-tagging mode_ chosen: " << mode_  << ". Choose between Mean, 1sigmaUp, or 1sigmaDown." ; 
        return btaggedJetsWithSF_ ; 
      }
    }
    else if (algo_== "CSVT") {
      if (strcmp(mode_.c_str(), "Mean") == 0) {
        LightJetSF = (GetSFlmean("CSV","T", 0., 2.4, "ABCD"))->Eval(jet_et) ; 
      }
      else if (strcmp(mode_.c_str(), "1sigmaUp") == 0) {
        LightJetSF = (GetSFlmax("CSV","T", 0., 2.4, "ABCD"))->Eval(jet_et) ; 
      }
      else if (strcmp(mode_.c_str(), "1sigmaDown") == 0) {
        LightJetSF = (GetSFlmin("CSV","T", 0., 2.4, "ABCD"))->Eval(jet_et) ; 
      }
      else {
        edm::LogError("ApplyBTagSF") << " Wrong b-tagging mode_ chosen: " << mode_  << ". Choose between Mean, 1sigmaUp, or 1sigmaDown." ; 
        return btaggedJetsWithSF_ ; 
      }
    }
    else {
      edm::LogError("ApplyBTagSF") << " Wrong b-tagging algo_ chosen: " << algo_ << ". Choose between CSVL, CSVM, or CSVT." ; 
      return btaggedJetsWithSF_ ; 
    }

    double phi = jet_phi;
    double sin_phi = sin(phi*1000000);
    double seed = abs(static_cast<int>(sin_phi*100000));

    //Initialize class
    BTagSFUtil* btsfutil = new BTagSFUtil( seed );

    btsfutil->modifyBTagsWithSF(isTagged, jet_flavor, BTagSF, BTageff, LightJetSF, LightJeteff);

    if (isTagged) btaggedJetsWithSF_.push_back(*ijet) ; 

    delete btsfutil ; 

  }

  return btaggedJetsWithSF_ ; 
}

