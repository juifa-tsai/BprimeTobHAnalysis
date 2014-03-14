// system include files
#include <memory>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <assert.h>
#include <vector>
#include <map>

// Root headers 
#include <TLorentzVector.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "BpbH/BprimeTobH/interface/format.h"

#include "BpbH/BprimeTobHAnalysis/interface/HiggsBRscaleFactors.h" 

//
// class declaration
//

class BprimeGenLevelSelector : public edm::EDAnalyzer {
  public:
    // ----------member data ---------------------------
    explicit BprimeGenLevelSelector(const edm::ParameterSet&);
    ~BprimeGenLevelSelector();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    void CreateHistos(const TString&) ; 
    void AddHisto(const TString&, const TString&, const TString&, const int&, const double&, const double&) ; 
    template <class Type>
      void FillHisto(const TString& name, const Type value, const double weight);

    TChain*            chain_;

    EvtInfoBranches    EvtInfo;
    VertexInfoBranches VtxInfo;
    GenInfoBranches    GenInfo;
    JetInfoBranches    GenJetInfo;

    edm::Service<TFileService> fs; 
    TH1D *h_bprimeDecayModes, *h_cutflow ; 
    std::map<TString, TH1D*> hmap_1d ;  

    bool isData_ ; 
    double evtwt_ ; 

    //// Configurables 

    int                             maxEvents_; 
    const int                       reportEvery_; 
    const std::string               inputTTree_;
    const std::vector<std::string>  inputFiles_;
};

BprimeGenLevelSelector::BprimeGenLevelSelector(const edm::ParameterSet& iConfig) :
  maxEvents_(iConfig.getParameter<int>("MaxEvents")), 
  reportEvery_(iConfig.getParameter<int>("ReportEvery")),
  inputTTree_(iConfig.getParameter<std::string>("InputTTree")),
  inputFiles_(iConfig.getParameter<std::vector<std::string> >("InputFiles")) 
{
}

BprimeGenLevelSelector::~BprimeGenLevelSelector() { 
  delete chain_;
}

void BprimeGenLevelSelector::beginJob() { 
  chain_ = new TChain(inputTTree_.c_str());

  for(unsigned i=0; i<inputFiles_.size(); ++i) {
    chain_->Add(inputFiles_.at(i).c_str());

    TFile *f = TFile::Open(inputFiles_.at(i).c_str(),"READ");
    f->Close();
  }

  EvtInfo.Register(chain_);
  VtxInfo.Register(chain_);
  GenInfo.Register(chain_);
  GenJetInfo.Register(chain_,"GenJetInfo");

  if(maxEvents_<0 || maxEvents_>chain_->GetEntries()) maxEvents_ = chain_->GetEntries();

  TString cutname = "SingleBp";
  CreateHistos(cutname) ; 

}

void BprimeGenLevelSelector::CreateHistos(const TString& cutname) {

  AddHisto(cutname, "_ptgen_bp"    ,";p_{T} (b');"       ,200 ,0. ,1000.) ;
  AddHisto(cutname, "_ptgen_bpbar" ,";p_{T} (#bar{b'});" ,200 ,0. ,1000.) ;
  AddHisto(cutname, "_ygen_bp"     ,";y (b');"           ,100 ,-4. ,4.) ;
  AddHisto(cutname, "_ygen_bpbar"  ,";y (#bar{b'});"     ,100 ,-4. ,4.) ;

  AddHisto(cutname, "_ptgen_H_bp"    ,";p_{T} (H from b');"       ,200 ,0. ,1000.) ;
  AddHisto(cutname, "_ptgen_H_bpbar" ,";p_{T} (H from #bar{b'});" ,200 ,0. ,1000.) ;
  AddHisto(cutname, "_ygen_H_bp"     ,";y (H from b');"           ,100 ,-5. ,5.) ;
  AddHisto(cutname, "_ygen_H_bpbar"  ,";y (H from #bar{b'});"     ,100 ,-5. ,5.) ;

  AddHisto(cutname, "_ptgen_b_bp"    ,";p_{T} (b from b');"             ,200 ,0. ,1000.) ;
  AddHisto(cutname, "_ptgen_b_bpbar" ,";p_{T} (#bar{b} from #bar{b'});" ,200 ,0. ,1000.) ;
  AddHisto(cutname, "_etagen_b_bp"   ,";#eta (b from b');"              ,100 ,-5. ,5.) ;
  AddHisto(cutname, "_etagen_b_bpbar",";#eta (#bar{b} from #bar{b'});"  ,100 ,-5. ,5.) ;

  AddHisto(cutname, "_DRHbgen_bpbar",";#DeltaR(H,b);" ,50 ,0. ,5.) ;
  AddHisto(cutname, "_DRHbgen_bp"   ,";#DeltaR(H,b);" ,50 ,0. ,5.) ;
  AddHisto(cutname, "_DRbbgen_bpbar",";#DeltaR(b,b);" ,50 ,0. ,5.) ;
  AddHisto(cutname, "_DRbbgen_bp"   ,";#DeltaR(b,b);" ,50 ,0. ,5.) ;

  h_bprimeDecayModes = fs->make<TH1D>("h_bprimeDecayModes", "b' decay modes", 11, 0, 11) ; 
  h_bprimeDecayModes -> GetXaxis() -> SetBinLabel(1  ,"N(b')") ; 
  h_bprimeDecayModes -> GetXaxis() -> SetBinLabel(2  ,"b'b' -> bHbH") ; 
  h_bprimeDecayModes -> GetXaxis() -> SetBinLabel(3  ,"b'b' -> bHbZ") ; 
  h_bprimeDecayModes -> GetXaxis() -> SetBinLabel(4  ,"b'b' -> bHtW") ; 
  h_bprimeDecayModes -> GetXaxis() -> SetBinLabel(5  ,"b'b' -> bZbZ") ; 
  h_bprimeDecayModes -> GetXaxis() -> SetBinLabel(6  ,"b'b' -> bZtW") ; 
  h_bprimeDecayModes -> GetXaxis() -> SetBinLabel(7  ,"b'b' -> tWtW") ; 
  h_bprimeDecayModes -> GetXaxis() -> SetBinLabel(8  ,"b'   -> bH")   ; 
  h_bprimeDecayModes -> GetXaxis() -> SetBinLabel(9  ,"b'   -> bZ")   ; 
  h_bprimeDecayModes -> GetXaxis() -> SetBinLabel(10 ,"b'   -> tW")   ; 
  h_bprimeDecayModes -> GetXaxis() -> SetBinLabel(11 ,"Others") ; 

  return ; 

}

void BprimeGenLevelSelector::AddHisto(const TString& cutname, const TString& histname, const TString& histtitle, const int& nbins, const double& min, const double& max) { 
  TH1D* h1d ; 
  h1d = fs->make<TH1D>(cutname+histname, cutname+histtitle, nbins, min, max);  
  h1d -> Sumw2() ; 
  hmap_1d[cutname+histname] = h1d ; 

  return ; 
}

template <class Type>
void BprimeGenLevelSelector::FillHisto(const TString& name, const Type value, const double weight){
  hmap_1d[name]->Fill(double(value),weight);
  return ; 
}

void BprimeGenLevelSelector::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) { 
  using namespace edm;
  using namespace std;

  if(chain_ == 0) return;

  edm::LogInfo("StartingAnalysisLoop") << "Starting analysis loop\n";

  for(int entry=0; entry<maxEvents_; entry++) {

    if((entry%reportEvery_) == 0) edm::LogInfo("Event") << entry << " of " << maxEvents_ ; 

    chain_->GetEntry(entry);

    //std::cout << EvtInfo.RunNo << ":" << EvtInfo.LumiNo << ":" << EvtInfo.EvtNo << std::endl ; 

    if ( !isData_ ) { //// Gen info 
      std::vector<int> gen_bp_indices, gen_higgs_indices ; 
      for (int igen = 0; igen < GenInfo.Size; ++igen) { 
        if ( abs(GenInfo.PdgID[igen]) == 7 ) { //// bprime found 
          gen_bp_indices.push_back(igen) ; 
        }
        if ( GenInfo.Status[igen] == 3 && abs(GenInfo.PdgID[igen]) == 25 && GenInfo.nDa[igen] >= 2  ) gen_higgs_indices.push_back(igen) ; 
      }
      int nbprimeBH(0), nbprimeBZ(0), nbprimeTW(0) ; 
      for (std::vector<int>::const_iterator ibp = gen_bp_indices.begin(); ibp != gen_bp_indices.end(); ++ibp) {
        int bprimeDau0(abs(GenInfo.Da0PdgID[*ibp])), bprimeDau1(abs(GenInfo.Da1PdgID[*ibp])) ;
        if ( (bprimeDau0 == 5 && bprimeDau1 == 25) || (bprimeDau0 == 25 && bprimeDau1 == 5) ) { //// b' -> bH  
          ++nbprimeBH ; 
          TLorentzVector p4_bp, p4_Higgs, p4_b ;
          p4_bp.SetPtEtaPhiM(GenInfo.Pt[*ibp], GenInfo.Eta[*ibp], GenInfo.Phi[*ibp], GenInfo.Mass[*ibp]) ; 
          if (bprimeDau0 == 25 && bprimeDau1 == 5) {
            p4_Higgs.SetPtEtaPhiM(GenInfo.Da0Pt[*ibp], GenInfo.Da0Eta[*ibp], GenInfo.Da0Phi[*ibp], GenInfo.Da0Mass[*ibp]) ; 
            p4_b.SetPtEtaPhiM(GenInfo.Da1Pt[*ibp], GenInfo.Da1Eta[*ibp], GenInfo.Da1Phi[*ibp], GenInfo.Da1Mass[*ibp]) ; 
          }
          else if (bprimeDau0 == 5 && bprimeDau1 == 25) {
            p4_b.SetPtEtaPhiM(GenInfo.Da0Pt[*ibp], GenInfo.Da0Eta[*ibp], GenInfo.Da0Phi[*ibp], GenInfo.Da0Mass[*ibp]) ; 
            p4_Higgs.SetPtEtaPhiM(GenInfo.Da1Pt[*ibp], GenInfo.Da1Eta[*ibp], GenInfo.Da1Phi[*ibp], GenInfo.Da1Mass[*ibp]) ; 
          }
          if ( GenInfo.PdgID[*ibp] == 7 ) {
            FillHisto("SingleBp_ptgen_bp",   p4_bp.Pt(), 1) ; 
            FillHisto("SingleBp_ygen_bp",    p4_bp.Rapidity(), 1) ; 
            FillHisto("SingleBp_ptgen_H_bp", p4_Higgs.Pt(), 1) ; 
            FillHisto("SingleBp_ygen_H_bp",  p4_Higgs.Rapidity(), 1) ; 
            FillHisto("SingleBp_ptgen_b_bp", p4_b.Pt(), 1) ; 
            FillHisto("SingleBp_etagen_b_bp",p4_b.Eta(), 1) ; 
            FillHisto("SingleBp_DRHbgen_bp" ,p4_b.DeltaR(p4_Higgs), 1) ; 

          }
          else if (GenInfo.PdgID[*ibp] == -7) {
            FillHisto("SingleBp_ptgen_bpbar",   p4_bp.Pt(), 1) ; 
            FillHisto("SingleBp_ygen_bpbar",    p4_bp.Rapidity(), 1) ; 
            FillHisto("SingleBp_ptgen_H_bpbar", p4_Higgs.Pt(), 1) ; 
            FillHisto("SingleBp_ygen_H_bpbar",  p4_Higgs.Rapidity(), 1) ; 
            FillHisto("SingleBp_ptgen_b_bpbar", p4_b.Pt(), 1) ; 
            FillHisto("SingleBp_etagen_b_bpbar",p4_b.Eta(), 1) ; 
            FillHisto("SingleBp_DRHbgen_bpbar" ,p4_b.DeltaR(p4_Higgs), 1) ; 
          }
        }
        else if ( (bprimeDau0 == 5 && bprimeDau1 == 23) || (bprimeDau0 == 23 && bprimeDau1 == 5) ) ++nbprimeBZ ; //// b' -> bZ 
        else if ( (bprimeDau0 == 6 && bprimeDau1 == 24) || (bprimeDau0 == 24 && bprimeDau1 == 6) ) ++nbprimeTW ; //// b' -> tW 
      }
      h_bprimeDecayModes -> AddBinContent(gen_bp_indices.size()) ; 
      if      ( nbprimeBH == 2 )                                    h_bprimeDecayModes -> AddBinContent(2)  ; 
      else if ( nbprimeBH == 1 && nbprimeBZ == 1 )                  h_bprimeDecayModes -> AddBinContent(3)  ; 
      else if ( nbprimeBH == 1 && nbprimeTW == 1 )                  h_bprimeDecayModes -> AddBinContent(4)  ; 
      else if ( nbprimeBZ == 2 )                                    h_bprimeDecayModes -> AddBinContent(5)  ; 
      else if ( nbprimeBZ == 1 && nbprimeTW == 1 )                  h_bprimeDecayModes -> AddBinContent(6)  ; 
      else if ( nbprimeTW == 2 )                                    h_bprimeDecayModes -> AddBinContent(7)  ; 
      else if ( nbprimeBH == 1 && nbprimeBZ == 0 && nbprimeTW == 0) h_bprimeDecayModes -> AddBinContent(8)  ; 
      else if ( nbprimeBH == 0 && nbprimeBZ == 1 && nbprimeTW == 0) h_bprimeDecayModes -> AddBinContent(9)  ; 
      else if ( nbprimeBH == 0 && nbprimeBZ == 0 && nbprimeTW == 1) h_bprimeDecayModes -> AddBinContent(10) ; 
      else                                                          h_bprimeDecayModes -> AddBinContent(11) ; 
      //// Higgs BR reweighting
      for (std::vector<int>::const_iterator ihig = gen_higgs_indices.begin(); ihig != gen_higgs_indices.end(); ++ihig) {
        int higgsDau0(abs(GenInfo.Da0PdgID[*ihig])), higgsDau1(abs(GenInfo.Da1PdgID[*ihig])) ;
        int higgsDau = higgsDau0+higgsDau1 ; 
        if(higgsDau==10) evtwt_ *= HiggsBRscaleFactors::higgsBBSf;
        else if(higgsDau==30) evtwt_ *= HiggsBRscaleFactors::higgsTauTauSf;
        else if(higgsDau==26) evtwt_ *= HiggsBRscaleFactors::higgsMuMuSf;
        else if(higgsDau==8)  evtwt_ *= HiggsBRscaleFactors::higgsCCSf;
        else if(higgsDau==6)  evtwt_ *= HiggsBRscaleFactors::higgsSSSf;
        else if(higgsDau==16) evtwt_ *= HiggsBRscaleFactors::higgsTTSf;
        else if(higgsDau==42) evtwt_ *= HiggsBRscaleFactors::higgsGGSf;
        else if(higgsDau==44) evtwt_ *= HiggsBRscaleFactors::higgsGammaGammaSf;
        else if(higgsDau==45) evtwt_ *= HiggsBRscaleFactors::higgsZGammaSf;
        else if(higgsDau==48) evtwt_ *= HiggsBRscaleFactors::higgsWWSf;
        else if(higgsDau==46) evtwt_ *= HiggsBRscaleFactors::higgsZZSf; 
        else evtwt_ *= 1 ; 
      } //// Higgs BR reweighting  
    } //// Gen info

  } //// entry loop 

}

// ------------ method called once each job just after ending the event loop  ------------
void BprimeGenLevelSelector::endJob() { 
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void BprimeGenLevelSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BprimeGenLevelSelector);
