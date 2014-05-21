#ifndef BPRIMETOBHANALYSIS_INTERFACE_APPLYBTAGSF_H
#define BPRIMETOBHANALYSIS_INTERFACE_APPLYBTAGSF_H

#include "BpbH/BprimeTobHAnalysis/interface/BTagSFUtil.h"
#include "BpbH/BprimeTobH/interface/JetCollection.h" 

#include <string>

class ApplyBTagSF {

  public:
    //ApplyBTagSF(JetCollection&, double, std::string, std::string) ;
    ApplyBTagSF(JetCollection&, double, double, double) ; 
    ~ApplyBTagSF() ; 

    JetCollection getBtaggedJetsWithSF () ; 

  private:

    JetCollection jets_ ; 
    double bDisc_ ; 
    std::string algo_ ; 
    JetCollection btaggedJetsWithSF_ ; 

    double SFbShift_ ; 
    double SFlShift_ ; 

    static std::map<double, std::string> btagOP_ ; 

};

#endif 
