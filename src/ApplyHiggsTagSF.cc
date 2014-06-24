#include "BpbH/BprimeTobHAnalysis/interface/ApplyHiggsTagSF.h" 
#include "BpbH/BprimeTobHAnalysis/interface/SFb-pt_WITHttbar_payload_EPS13.h"
#include "BpbH/BprimeTobHAnalysis/interface/SFlightFuncs_EPS2013.h"

ApplyHiggsTagSF::ApplyHiggsTagSF (double ptFatJet, double ptSubjet1, double ptSubjet2, double etaFatJet, double etaSubjet1, double etaSubjet2, double phiSubjet1, double phiSubjet2, int flavSubjet1, int flavSubjet2,  double csvSubjet1, double csvSubjet2, double SFbShift, double SFlShift) :
  ptFatJet_(ptFatJet), 
  ptSubjet1_(ptSubjet1),
  ptSubjet2_(ptSubjet2),
  etaFatJet_(etaFatJet), 
  etaSubjet1_(etaSubjet1),
  etaSubjet2_(etaSubjet2),
  phiSubjet1_(phiSubjet1),
  phiSubjet2_(phiSubjet2),
  flavSubjet1_(flavSubjet1), 
  flavSubjet2_(flavSubjet2), 
  csvSubjet1_(csvSubjet1),
  csvSubjet2_(csvSubjet2), 
  SFbShift_(SFbShift),
  SFlShift_(SFlShift), 
  higgsTagSF_(1.) 
{

  if (ptSubjet1 < 20.)  ptSubjet1 = 20. ; 
  if (ptSubjet1 > 800.) ptSubjet1 = 799.9; 
  if (ptSubjet2 < 20.)  ptSubjet2 = 20. ; 
  if (ptSubjet2 > 800.) ptSubjet2 = 799.9; 

  double drSubjets = sqrt( pow((etaSubjet1_ - etaSubjet2_),2.) + pow((phiSubjet1_ - phiSubjet2_),2.)) ;
  if (drSubjets > 0.3 && drSubjets <= 0.4) {
    SFbShift_ *= 2. ; 
    SFlShift_ *= 2. ; 
  }

  double sf1(1.), sf2(1.)  ; 

  flavSubjet1_ = abs(flavSubjet1_) ; 
  flavSubjet2_ = abs(flavSubjet2_) ; 

  if (flavSubjet1_ == 4 || flavSubjet1_ == 5) {
    double errscale(1.);
    if ( flavSubjet1_ == 4) errscale *= 2. ; 
    int ptbin(-1) ; 
    for (int ii = 0; ii < int(sizeof(ptmin)/sizeof(ptmin[0])); ++ii) {  
      if ( ptSubjet1 >= ptmin[ii] && ptSubjet1 < ptmax[ii] ) {
        ptbin = ii ; 
      }
    }
    sf1 = SFb_CSVM->Eval(ptSubjet1)*( 1 + (SFbShift_*SFb_CSVM_error[ptbin]) ) ; 
  }
  else if (flavSubjet1_ < 4 || flavSubjet1_ == 21) {
    for (int ii = 0; ii < 4; ++ii) {
      if ( etaSubjet1_ >= SFlight_CSVM_etamin[ii] && etaSubjet1_ < SFlight_CSVM_etamax[ii] ) {
        sf1 = GetSFlmean("CSV","M",SFlight_CSVM_etamin[ii], SFlight_CSVM_etamax[ii], "ABCD")->Eval(ptSubjet1) ; 
        if (SFlShift_ > 0) sf1 += SFlShift_*( (GetSFlmax("CSV","M",SFlight_CSVM_etamin[ii], SFlight_CSVM_etamax[ii], "ABCD"))->Eval(ptSubjet1) - sf1 ) ; 
        else if (SFlShift_ < 0) sf1 += abs(SFlShift_)*( (GetSFlmin("CSV","M",SFlight_CSVM_etamin[ii], SFlight_CSVM_etamax[ii], "ABCD"))->Eval(ptSubjet1) - sf1 ) ; 
        break ;
      }
    }
  }

  if (flavSubjet2_ == 4 || flavSubjet2_ == 5) {
    double errscale(1.);
    if ( flavSubjet2_ == 4) errscale *= 2. ; 
    int ptbin(-1) ; 
    for (int ii = 0; ii < int(sizeof(ptmin)/sizeof(ptmin[0])); ++ii) {  
      if ( ptSubjet2 >= ptmin[ii] && ptSubjet2 < ptmax[ii] ) {
        ptbin = ii ; 
      }
    }
    sf2 = SFb_CSVM->Eval(ptSubjet2)*( 1 + (SFbShift_*SFb_CSVM_error[ptbin]) ) ; 
  }
  else if (flavSubjet2_ < 4 || flavSubjet2_ == 21) {
    for (int ii = 0; ii < 4; ++ii) {
      if ( etaSubjet2_ >= SFlight_CSVM_etamin[ii] && etaSubjet2_ < SFlight_CSVM_etamax[ii] ) {
        sf2 = GetSFlmean("CSV","M",SFlight_CSVM_etamin[ii], SFlight_CSVM_etamax[ii], "ABCD")->Eval(ptSubjet2) ; 
        if (SFlShift_ > 0) sf2 += SFlShift_*( (GetSFlmax("CSV","M",SFlight_CSVM_etamin[ii], SFlight_CSVM_etamax[ii], "ABCD"))->Eval(ptSubjet2) - sf2 ) ; 
        else if (SFlShift_ < 0) sf2 += abs(SFlShift_)*( (GetSFlmin("CSV","M",SFlight_CSVM_etamin[ii], SFlight_CSVM_etamax[ii], "ABCD"))->Eval(ptSubjet2) - sf2 ) ; 
        break ;
      }
    }
  }

  higgsTagSF_ = sf1*sf2 ; 

  double SF_mass_partonshower = 1.007 ; 
  double SF_mass_frag = 0.955 ; 
  double SF_nsubjettiness_etaLE1 = 0.967 ; 
  double SF_nsubjettiness_etaGT1 = 0.967 ; 

  double err_SF_mass_partonshower = 0.004 ; 
  double err_SF_mass_frag = 0.011 ; 
  double err_SF_nsubjettiness_etaLE1 = .024*0.967 ; 
  double err_SF_nsubjettiness_etaGT1 = .056*0.967 ; 

  if ( abs(etaFatJet) <= 1. ) higgsTagSF_ *= SF_mass_partonshower*SF_mass_frag*SF_nsubjettiness_etaLE1 ; 
  else if ( abs(etaFatJet) > 1. && abs(etaFatJet) < 2.4 ) higgsTagSF_ *= SF_mass_partonshower*SF_mass_frag*SF_nsubjettiness_etaGT1 ; 

}

