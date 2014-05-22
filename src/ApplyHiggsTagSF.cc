#include "BpbH/BprimeTobHAnalysis/interface/ApplyHiggsTagSF.h" 
#include "BpbH/BprimeTobHAnalysis/interface/SFb-pt_WITHttbar_payload_EPS13.h"
#include "BpbH/BprimeTobHAnalysis/interface/SFlightFuncs_EPS2013.h"

ApplyHiggsTagSF::ApplyHiggsTagSF (double ptSubjet1, double ptSubjet2, double etaSubjet1, double etaSubjet2, int flavSubjet1, int flavSubjet2,  double csvSubjet1, double csvSubjet2, double SFbShift, double SFlShift) :
  ptSubjet1_(ptSubjet1),
  ptSubjet2_(ptSubjet2),
  etaSubjet1_(etaSubjet1),
  etaSubjet2_(etaSubjet2),
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

  double sf1(1.), sf2(1.)  ; 

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

}

