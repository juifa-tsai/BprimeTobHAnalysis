#ifndef BPRIMETOBHANALYSIS_INTERFACE_APPLYBTAGSF_H
#define BPRIMETOBHANALYSIS_INTERFACE_APPLYBTAGSF_H

#include "BpbH/BprimeTobHAnalysis/interface/BTagSFUtil.h"
#include "BpbH/BprimeTobH/interface/JetCollection.h" 

#include <string>

class ApplyBTagSF {

  public:
  ApplyBTagSF(JetCollection&, double, std::string, std::string) ;
  ~ApplyBTagSF() ; 

  JetCollection getBtaggedJetsWithSF () ; 

  private:

  JetCollection jets_ ; 
  double bDisc_ ; 
  JetCollection btaggedJetsWithSF_ ; 

  std::string algo_ ;
  std::string mode_ ; 

};

#endif 
