// -*- C++ -*-
//
// Package:    PdfSystematicsAnalyzer
// Class:      PdfSystematicsAnalyzer
// 
/**\class PdfSystematicsAnalyzer PdfSystematicsAnalyzer.cc Bprime_kit/PdfSystematicsAnalyzer/src/PdfSystematicsAnalyzer.cc

Description: 
Analyzer class for Bprime -> b Higgs studies 
- National Taiwan University - 

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Eleni Petrakou,27 2-020,+41227674870,
//         Created:  Tue Jul 16 19:48:47 CEST 2013
// Second Author:    Devdatta Majumder 
// $Id$
//
//

// system include files
#include <memory>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <assert.h>
#include <vector>
#include <map>
#include <string>

// Root headers 
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "BpbH/BprimeTobHAnalysis/interface/PDFTree.h"

namespace LHAPDF {
  void initPDFSet(int nset, const std::string& filename, int member=0);
  int numberPDF(int nset);
  void usePDFMember(int nset, int member);
  double xfx(int nset, double x, double Q, int fl);
  double getXmin(int nset, int member);
  double getXmax(int nset, int member);
  double getQ2min(int nset, int member);
  double getQ2max(int nset, int member);
  void extrapolate(bool extrapolate=true);
}

//
// class declaration
//

class PdfSystematicsAnalyzer : public edm::EDAnalyzer {
  public:
    explicit PdfSystematicsAnalyzer(const edm::ParameterSet&);
    ~PdfSystematicsAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    // ----------member data ---------------------------

    //// Configurables 

    int                             maxEvents_; 
    const int                       reportEvery_; 
    const std::string               inputTTree_;
    const std::vector<std::string>  inputFiles_;
    std::vector<std::string>       pdfSetNames_;

    std::vector<double> pdf_weight0_sum_ ; 
    std::vector<std::vector<double> >pdf_weights_sum_ ; 
    std::vector<double> passed_pdf_weight0_sum_ ; 
    std::vector<std::vector<double> >passed_pdf_weights_sum_ ; 

    TChain* chain_;

    PDFTree pdfTree_ ; 

    edm::Service<TFileService> fs; 

};

//
// constructors and destructor
//
PdfSystematicsAnalyzer::PdfSystematicsAnalyzer(const edm::ParameterSet& iConfig) : 
  maxEvents_(iConfig.getParameter<int>("MaxEvents")), 
  reportEvery_(iConfig.getParameter<int>("ReportEvery")),
  inputTTree_(iConfig.getParameter<std::string>("InputTTree")),
  inputFiles_(iConfig.getParameter<std::vector<std::string> >("InputFiles")),
  pdfSetNames_(iConfig.getUntrackedParameter<std::vector<std::string> > ("PdfSetNames"))  
{ 


  if (pdfSetNames_.size()>3) {
    edm::LogWarning("PdfSystematicsAnalyzer") << pdfSetNames_.size() << " PDF sets requested on input. Using only the first 3 sets and ignoring the rest!!";
    pdfSetNames_.erase(pdfSetNames_.begin()+3,pdfSetNames_.end());
  }

}


PdfSystematicsAnalyzer::~PdfSystematicsAnalyzer() { 
  delete chain_;
}

// ------------ method called once each job just before starting event loop  ------------
void PdfSystematicsAnalyzer::beginJob() { 

  chain_ = new TChain(inputTTree_.c_str());

  for(unsigned i=0; i<inputFiles_.size(); ++i) {
    chain_->Add(inputFiles_.at(i).c_str());
    TFile *f = TFile::Open(inputFiles_.at(i).c_str(),"READ");
    f->Close();
  }

  pdfTree_.Register(chain_);

  if(maxEvents_<0 || maxEvents_>chain_->GetEntries()) maxEvents_ = chain_->GetEntries();

  pdf_weight0_sum_.reserve(pdfSetNames_.size()) ; 
  passed_pdf_weight0_sum_.reserve(pdfSetNames_.size()) ; 
  pdf_weights_sum_.reserve(pdfSetNames_.size()) ; 
  passed_pdf_weights_sum_.reserve(pdfSetNames_.size()) ; 

  for (unsigned int ik = 0; ik < pdfSetNames_.size(); ++ik) {
    LHAPDF::initPDFSet(ik+1,pdfSetNames_[ik]) ;  
    unsigned int nweights = 1 ;
    if (LHAPDF::numberPDF(ik+1) > 1) { 
      nweights += LHAPDF::numberPDF(ik+1) ;
      std::vector<double> weights ; 
      weights.reserve(nweights) ; 
      for (unsigned int ii = 0; ii < nweights; ++ii) {
        weights.push_back(0) ; 
      }
      pdf_weights_sum_.push_back(weights) ; 
      passed_pdf_weights_sum_.push_back(weights) ; 
    }
  }

  return ;  

}

// ------------ method called for each event  ------------
void PdfSystematicsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) { 
  using namespace edm;
  using namespace std;

  if(chain_ == 0) return;

  edm::LogInfo("StartingAnalysisLoop") << "Starting analysis loop\n";


  for(int entry=0; entry<maxEvents_; entry++) {
    if((entry%reportEvery_) == 0) edm::LogInfo("Event") << entry << " of " << maxEvents_ ; 

    chain_->GetEntry(entry);

    for (unsigned int ik = 0; ik < pdfSetNames_.size(); ++ik) {

      for (unsigned int ipdf=1; ipdf <= pdf_weights_sum_.size(); ++ipdf) {
        LHAPDF::usePDFMember(ik+1, ipdf); 
        double xpdf1 = LHAPDF::xfx(ik+1, pdfTree_.PDFx1, pdfTree_.qScale, pdfTree_.PDFid1);
        double xpdf2 = LHAPDF::xfx(ik+1, pdfTree_.PDFx2, pdfTree_.qScale, pdfTree_.PDFid2);
        double weight = xpdf1 * xpdf2 ;
        (pdf_weights_sum_.at(ik)).at(ipdf-1) += weight ; 
        if (pdfTree_.evtwt > 0) (passed_pdf_weights_sum_.at(ik)).at(ipdf-1) += weight ; 
      }

    }

  } //// entry loop 

}

// ------------ method called once each job just after ending the event loop  ------------
void PdfSystematicsAnalyzer::endJob() { 

  for (unsigned int ik = 0; ik < pdfSetNames_.size(); ++ik) {
    double eff0 = (passed_pdf_weights_sum_.at(ik)).at(0)/(pdf_weights_sum_.at(ik)).at(0) ; 
    double deleffP(0), deleffM(0) ;
    edm::LogInfo("PdfSystematicsAnalyzer") << pdfSetNames_.at(ik) << " eff0 = " << eff0 ; 
    for (unsigned int ii = 1; ii < passed_pdf_weights_sum_.size(); ++ii) {
      double eff = (passed_pdf_weights_sum_.at(ik)).at(ii)/(pdf_weights_sum_.at(ik)).at(ii) ;
      if (ii%2 == 0) deleffP += pow((eff - eff0)/eff0, 2.) ;
      else deleffM += pow((eff - eff0)/eff0, 2.) ;
    }
    deleffP = sqrt(deleffP) ; deleffM = sqrt(deleffM) ; 
    edm::LogInfo("PdfSystematicsAnalyzer") << pdfSetNames_.at(ik) << " pdf uncertainy  +1sigma = " << 1.00+deleffP << " -1sigma = " << 1-deleffM << "%"  ;
  }

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void PdfSystematicsAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PdfSystematicsAnalyzer);
