#ifndef BPRIMETOBHANALYSIS_INTERFACE_APPLYHIGGSTAGSF_H
#define BPRIMETOBHANALYSIS_INTERFACE_APPLYHIGGSTAGSF_H

#include "BpbH/BprimeTobH/interface/Jet.h" 

class ApplyHiggsTagSF {
  public: 
    ApplyHiggsTagSF (double ptSubjet1, double ptSubjet2, double etaSubjet1, double etaSubjet2, int flavSubjet1, int flavSubjet2,  double csvSubjet1, double csvSubjet2) ; 
    ~ApplyHiggsTagSF() { higgsTagSF_ = 1. ; } 
    double GetHiggsTagSF () { return higgsTagSF_ ; } 

  private:
    double ptSubjet1_   ; 
    double etaSubjet2_  ; 
    double etaSubjet1_  ; 
    double ptSubjet2_   ; 
    int    flavSubjet1_ ; 
    int    flavSubjet2_ ; 
    double csvSubjet1_  ;
    double csvSubjet2_  ; 

    double higgsTagSF_ ; 

    static double higgsBBSf         ; 
    static double higgsTauTauSf     ; 
    static double higgsMuMuSf       ; 
    static double higgsCCSf         ; 
    static double higgsSSSf         ; 
    static double higgsTTSf         ; 

    static double higgsGGSf         ; 
    static double higgsGammaGammaSf ; 
    static double higgsZGammaSf     ; 
    static double higgsWWSf         ; 
    static double higgsZZSf         ; 

} ; 

#endif 
