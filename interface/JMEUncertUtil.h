#ifndef BPRIMETOBHANALYSIS_INTERFACE_JMEUNCERTUTIL_H
#define BPRIMETOBHANALYSIS_INTERFACE_JMEUNCERTUTIL_H

#include <string>
#include <vector>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "BpbH/BprimeTobH/interface/JetCollection.h"

class JetInfoBranches;

using std::string;

class JMEUncertUtil {

  public: 
    JMEUncertUtil (const edm::ParameterSet& iConfig, JetCollection &jets, std::string stype, double jecShift) ; 
    ~JMEUncertUtil () ; 

    JetCollection GetModifiedJetColl() const ;  

    static const string JEC_Types[];
    enum JEC_TYPE {NOSYS, JESAK5MC, JESCA8MC, JERAK5MC, JERCA8MC, NTYPES}  ; 

  private:

    void   setType (std::string); 
    void jesUncert () ; 
    void jerScale () ; 

    std::vector<double> jerEta_ ; 
    std::vector<double> jerNominal_ ; // Measured data/MC ratio
    std::vector<double> jerSigmaSym_ ; 
    std::vector<double> jerSigmaNeg_ ; 
    std::vector<double> jerSigmaPos_ ; 
    std::string stype_ ; 
    double jecShift_ ; 
    const JetCollection jets_ ; 

    JEC_TYPE jecType_ ; 

    JetCollection modifiedJets_ ; 
    JetCorrectionUncertainty* jecUncert_ ; 

}; 

#endif 
