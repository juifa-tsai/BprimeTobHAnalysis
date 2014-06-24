#ifndef BPRIMETOBHANALYSIS_INTERFACE_APPLYHIGGSTAGSF_H
#define BPRIMETOBHANALYSIS_INTERFACE_APPLYHIGGSTAGSF_H

#include "BpbH/BprimeTobH/interface/Jet.h" 

class ApplyHiggsTagSF {
  public: 
    ApplyHiggsTagSF (double ptFatJet, double ptSubjet1, double ptSubjet2, double etaFatJet, double etaSubjet1, double etaSubjet2, double phiSubjet1, double phiSubjet2, int flavSubjet1, int flavSubjet2,  double csvSubjet1, double csvSubjet2, double SFbShift, double SFlShift) ; 
    ~ApplyHiggsTagSF() { higgsTagSF_ = 1. ; } 
    double GetHiggsTagSF () { return higgsTagSF_ ; } 

  private:
    double ptFatJet_    ; 
    double ptSubjet1_   ; 
    double ptSubjet2_   ; 
    double etaFatJet_   ; 
    double etaSubjet1_  ; 
    double etaSubjet2_  ; 
    double phiSubjet1_  ; 
    double phiSubjet2_  ; 
    int    flavSubjet1_ ; 
    int    flavSubjet2_ ; 
    double csvSubjet1_  ;
    double csvSubjet2_  ; 

    double SFbShift_ ; 
    double SFlShift_ ; 

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
