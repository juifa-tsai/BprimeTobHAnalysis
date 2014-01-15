#include "BpbH/BprimeTobHAnalysis/interface/ApplyHiggsTagSF.h" 
#include "BpbH/BprimeTobHAnalysis/interface/SFb-pt_WITHttbar_payload_EPS13.h"
#include "BpbH/BprimeTobHAnalysis/interface/SFlightFuncs_EPS2013.h"

ApplyHiggsTagSF::ApplyHiggsTagSF (double ptSubjet1, double ptSubjet2, double etaSubjet1, double etaSubjet2, int flavSubjet1, int flavSubjet2,  double csvSubjet1, double csvSubjet2) :
  ptSubjet1_(ptSubjet1),
  ptSubjet2_(ptSubjet2),
  etaSubjet1_(etaSubjet1),
  etaSubjet2_(etaSubjet2),
  flavSubjet1_(flavSubjet1), 
  flavSubjet2_(flavSubjet2), 
  csvSubjet1_(csvSubjet1),
  csvSubjet2_(csvSubjet2), 
  higgsTagSF_(1.) 
{

  if (ptSubjet1 < 20.)  ptSubjet1 = 20. ; 
  if (ptSubjet1 > 800.) ptSubjet1 = 800. ; 
  if (ptSubjet2 < 20.)  ptSubjet2 = 20. ; 
  if (ptSubjet2 > 800.) ptSubjet2 = 800. ; 

  double sf1(1.), sf2(1.)  ; 

  if (flavSubjet1_ == 4 || flavSubjet1_ == 5) sf1 = SFb_CSVL->Eval(ptSubjet1) ; 
  else if (flavSubjet1_ < 4 || flavSubjet1_ == 21) {
    for (int ii = 0; ii < 4; ++ii) {
      if ( etaSubjet1_ >= SFlight_CSVL_etamin[ii] && etaSubjet1_ < SFlight_CSVL_etamax[ii] ) {
        sf1 = GetSFlmean("CSV","L",SFlight_CSVL_etamin[ii], SFlight_CSVL_etamax[ii], "ABCD")->Eval(ptSubjet1) ; 
      }
    }
  }

  if (flavSubjet2_ == 4 || flavSubjet2_ == 5) sf2 = SFb_CSVL->Eval(ptSubjet2) ; 
  else if (flavSubjet2_ < 4 || flavSubjet2_ == 21) {
    for (int ii = 0; ii < 4; ++ii) {
      if ( etaSubjet2_ >= SFlight_CSVL_etamin[ii] && etaSubjet2_ < SFlight_CSVL_etamax[ii] ) {
        sf2 = GetSFlmean("CSV","L",SFlight_CSVL_etamin[ii], SFlight_CSVL_etamax[ii], "ABCD")->Eval(ptSubjet2) ; 
      }
    }
  }

  higgsTagSF_ = sf1*sf2 ; 

}

