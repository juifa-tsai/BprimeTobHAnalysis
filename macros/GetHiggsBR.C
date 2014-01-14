#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TChain.h>
#include <TH1D.h>
#include <TMath.h>

#include <iostream>

#include "../../../BpbH/BprimeTobH/interface/format.h"

EvtInfoBranches EvtInfo ; 
GenInfoBranches GenInfo ; 

TString fnames [10] = {
  "root://eoscms//eos/cms/store/user/devdatta/NtuplesBprimeTobH_v1/BprimeBprimeToBHBHinc_M-1000_TuneZ2star_8TeV-madgraph/BprimeTobH_v1_10_1_Xgr.root" ,
  "root://eoscms//eos/cms/store/user/devdatta/NtuplesBprimeTobH_v1/BprimeBprimeToBHBHinc_M-1000_TuneZ2star_8TeV-madgraph/BprimeTobH_v1_1_1_zu3.root" ,
  "root://eoscms//eos/cms/store/user/devdatta/NtuplesBprimeTobH_v1/BprimeBprimeToBHBHinc_M-1000_TuneZ2star_8TeV-madgraph/BprimeTobH_v1_2_1_MLQ.root" ,
  "root://eoscms//eos/cms/store/user/devdatta/NtuplesBprimeTobH_v1/BprimeBprimeToBHBHinc_M-1000_TuneZ2star_8TeV-madgraph/BprimeTobH_v1_3_1_0K3.root" ,
  "root://eoscms//eos/cms/store/user/devdatta/NtuplesBprimeTobH_v1/BprimeBprimeToBHBHinc_M-1000_TuneZ2star_8TeV-madgraph/BprimeTobH_v1_4_1_n3y.root" ,
  "root://eoscms//eos/cms/store/user/devdatta/NtuplesBprimeTobH_v1/BprimeBprimeToBHBHinc_M-1000_TuneZ2star_8TeV-madgraph/BprimeTobH_v1_5_1_8aI.root" ,
  "root://eoscms//eos/cms/store/user/devdatta/NtuplesBprimeTobH_v1/BprimeBprimeToBHBHinc_M-1000_TuneZ2star_8TeV-madgraph/BprimeTobH_v1_6_1_vHy.root" ,
  "root://eoscms//eos/cms/store/user/devdatta/NtuplesBprimeTobH_v1/BprimeBprimeToBHBHinc_M-1000_TuneZ2star_8TeV-madgraph/BprimeTobH_v1_7_1_nOK.root" ,
  "root://eoscms//eos/cms/store/user/devdatta/NtuplesBprimeTobH_v1/BprimeBprimeToBHBHinc_M-1000_TuneZ2star_8TeV-madgraph/BprimeTobH_v1_8_1_sF1.root" ,
  "root://eoscms//eos/cms/store/user/devdatta/NtuplesBprimeTobH_v1/BprimeBprimeToBHBHinc_M-1000_TuneZ2star_8TeV-madgraph/BprimeTobH_v1_9_1_agx.root"  
} ; 

TString inputTree = "ntuple/tree" ; 

/// H(125) BR from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR#Branching_Ratios
////  H->bb, H->#tau#tau, H->#mu#mu, H->cc, H->ss, H->tt, H->gg, H->#gamma#gamma, H->Z#gamma, H->WW, H->ZZ, H->everything ; 
double higgsBR[12] = { 5.77E-01, 6.37E-02, 2.21E-04, 2.67E-02, 4.39E-04, 0.00E+00, 8.55E-02, 2.29E-03, 1.55E-03, 2.16E-01, 2.66E-02, 1.00 } ; 

void GetHiggsBR (int maxEvents) {

  TChain* chain = new TChain(inputTree);

  EvtInfo.Register(chain);
  GenInfo.Register(chain);

  for (int ii = 0; ii < 10; ++ii) {
    chain -> Add(fnames[ii]);
    TFile *fin = TFile::Open(fnames[ii],"READ");
    fin->Close();
  }

  TFile* fout = new TFile("HiggsBR.root", "RECREATE") ; 
  fout->cd() ; 
  TH1D* h_NHiggs = new TH1D("h_NHiggs", "Higgs boson decay channels;Decay channel;Events;", 12, 0., 12.) ; 
  h_NHiggs -> GetXaxis() -> SetBinLabel(1 , "H->bb") ; 
  h_NHiggs -> GetXaxis() -> SetBinLabel(2 , "H->#tau#tau") ; 
  h_NHiggs -> GetXaxis() -> SetBinLabel(3 , "H->#mu#mu") ; 
  h_NHiggs -> GetXaxis() -> SetBinLabel(4 , "H->cc") ; 
  h_NHiggs -> GetXaxis() -> SetBinLabel(5 , "H->ss") ; 
  h_NHiggs -> GetXaxis() -> SetBinLabel(6 , "H->tt") ; 
  h_NHiggs -> GetXaxis() -> SetBinLabel(7 , "H->gg") ; 
  h_NHiggs -> GetXaxis() -> SetBinLabel(8 , "H->#gamma#gamma") ; 
  h_NHiggs -> GetXaxis() -> SetBinLabel(9 , "H->Z#gamma") ; 
  h_NHiggs -> GetXaxis() -> SetBinLabel(10, "H->WW") ; 
  h_NHiggs -> GetXaxis() -> SetBinLabel(11, "H->ZZ") ; 
  h_NHiggs -> GetXaxis() -> SetBinLabel(12, "H->everything") ; 
  TH1D* h_BRHiggs = new TH1D("h_BRHiggs", "Higgs boson BRs;Decay channels;BR", 12, 0., 12.) ; 
  h_BRHiggs -> GetXaxis() -> SetBinLabel(1 , "H->bb") ; 
  h_BRHiggs -> GetXaxis() -> SetBinLabel(2 , "H->#tau#tau") ; 
  h_BRHiggs -> GetXaxis() -> SetBinLabel(3 , "H->#mu#mu") ; 
  h_BRHiggs -> GetXaxis() -> SetBinLabel(4 , "H->cc") ; 
  h_BRHiggs -> GetXaxis() -> SetBinLabel(5 , "H->ss") ; 
  h_BRHiggs -> GetXaxis() -> SetBinLabel(6 , "H->tt") ; 
  h_BRHiggs -> GetXaxis() -> SetBinLabel(7 , "H->gg") ; 
  h_BRHiggs -> GetXaxis() -> SetBinLabel(8 , "H->#gamma#gamma") ; 
  h_BRHiggs -> GetXaxis() -> SetBinLabel(9 , "H->Z#gamma") ; 
  h_BRHiggs -> GetXaxis() -> SetBinLabel(10, "H->WW") ; 
  h_BRHiggs -> GetXaxis() -> SetBinLabel(11, "H->ZZ") ; 
  h_BRHiggs -> GetXaxis() -> SetBinLabel(12, "H->everything") ; 
  TH1D* h_BRSFH125 = new TH1D("h_BRSFH125", "Higgs boson SF(BR);Decay channels;SF(BR)", 12, 0., 12.) ; 
  h_BRSFH125 -> GetXaxis() -> SetBinLabel(1 , "H->bb") ; 
  h_BRSFH125 -> GetXaxis() -> SetBinLabel(2 , "H->#tau#tau") ; 
  h_BRSFH125 -> GetXaxis() -> SetBinLabel(3 , "H->#mu#mu") ; 
  h_BRSFH125 -> GetXaxis() -> SetBinLabel(4 , "H->cc") ; 
  h_BRSFH125 -> GetXaxis() -> SetBinLabel(5 , "H->ss") ; 
  h_BRSFH125 -> GetXaxis() -> SetBinLabel(6 , "H->tt") ; 
  h_BRSFH125 -> GetXaxis() -> SetBinLabel(7 , "H->gg") ; 
  h_BRSFH125 -> GetXaxis() -> SetBinLabel(8 , "H->#gamma#gamma") ; 
  h_BRSFH125 -> GetXaxis() -> SetBinLabel(9 , "H->Z#gamma") ; 
  h_BRSFH125 -> GetXaxis() -> SetBinLabel(10, "H->WW") ; 
  h_BRSFH125 -> GetXaxis() -> SetBinLabel(11, "H->ZZ") ; 
  h_BRSFH125 -> GetXaxis() -> SetBinLabel(12, "H->everything") ; 

  if(maxEvents<0 || maxEvents>chain->GetEntries()) maxEvents = chain->GetEntries();

  for(int entry = 0; entry < maxEvents; entry++) { 
    if( (entry%1000) == 0 ) std::cout << " entry = " << entry << " of " << maxEvents << std::endl ;
    chain -> GetEntry(entry);

    for (int igen=0; igen < GenInfo.Size; ++igen) {

      if ( GenInfo.Status[igen] == 3 && 
          TMath::Abs(GenInfo.PdgID[igen])==25 && 
          GenInfo.nDa[igen] >= 2 ) { 

        int higgsDau = TMath::Abs(GenInfo.Da0PdgID[igen]) ; 
        if( TMath::Abs(GenInfo.Da0PdgID[igen]) == TMath::Abs(GenInfo.Da1PdgID[igen]) ) { 

          h_BRHiggs -> Fill ("H->everything" ,1.) ; h_NHiggs -> Fill ("H->everything" ,1.) ; 
          switch (higgsDau) {
            case 5 :
              h_BRHiggs -> Fill ("H->bb" ,1.) ; h_NHiggs -> Fill ("H->bb" ,1.) ; 
              break ; 
            case 15 : 
              h_BRHiggs -> Fill ("H->#tau#tau" ,1.) ; h_NHiggs -> Fill ("H->#tau#tau" ,1.) ; 
              break ; 
            case 13 :
              h_BRHiggs -> Fill ("H->#mu#mu" ,1.) ; h_NHiggs -> Fill ("H->#mu#mu" ,1.) ; 
              break ;
            case 4 :
              h_BRHiggs -> Fill ("H->cc" ,1.) ; h_NHiggs -> Fill ("H->cc" ,1.) ; 
              break ;
            case 3 :
              h_BRHiggs -> Fill ("H->ss" ,1.) ; h_NHiggs -> Fill ("H->ss" ,1.) ; 
              break ;
            case 6 :
              h_BRHiggs -> Fill ("H->tt" ,1.) ; h_NHiggs -> Fill ("H->tt" ,1.) ; 
              break ;
            case 21 :
              h_BRHiggs -> Fill ("H->gg" ,1.) ; h_NHiggs -> Fill ("H->gg" ,1.) ; 
              break ;
            case 22 : 
              h_BRHiggs -> Fill ("H->#gamma#gamma" ,1.) ; h_NHiggs -> Fill ("H->#gamma#gamma" ,1.) ; 
              break ;
            case 24 :
              h_BRHiggs -> Fill ("H->WW" ,1.) ; h_NHiggs -> Fill ("H->WW" ,1.) ; 
              break ;
            case 23 :
              h_BRHiggs -> Fill ("H->ZZ" ,1.) ; h_NHiggs -> Fill ("H->ZZ" ,1.) ; 
              break ;
          }
        }
        else if( (TMath::Abs(GenInfo.Da0PdgID[igen]) == 23 && TMath::Abs(GenInfo.Da1PdgID[igen]) == 22) || 
            (TMath::Abs(GenInfo.Da0PdgID[igen]) == 22 && TMath::Abs(GenInfo.Da1PdgID[igen]) == 23) 
            ) {
          h_BRHiggs -> Fill ("H->everything" ,1.) ; h_NHiggs -> Fill ("H->everything" ,1.) ; 
          h_BRHiggs -> Fill ("H->Z#gamma" ,1.) ; h_NHiggs -> Fill ("H->Z#gamma" ,1.) ; 
        }
        else {
          h_BRHiggs -> Fill ("H->everything" ,1.) ; h_NHiggs -> Fill ("H->everything" ,1.) ; 
          std::cout << " dau0 = " << GenInfo.Da0PdgID[igen] << " dau 1 = " << GenInfo.Da1PdgID[igen] << std::endl ; 
        }

      } //// If status 3 Higgs bosons with two daughters 

    } //// Loop over GenParticles

  } //// Loop over entries  

  h_BRHiggs -> Scale(1./h_BRHiggs -> GetBinContent(12)) ; 

  for (int ibr = 1; ibr <= h_BRHiggs->GetNbinsX(); ++ibr) { 
    std::cout << " BR: " << h_BRHiggs->GetXaxis()->GetBinLabel(ibr) << ": " << h_BRHiggs->GetBinContent(ibr) << std::endl ; 
  }

  for (int ibr = 1; ibr <= h_BRHiggs->GetNbinsX(); ++ibr) { 
    double brsf = h_BRHiggs->GetBinContent(ibr) > 0. ? higgsBR[ibr-1]/h_BRHiggs->GetBinContent(ibr) : 1.; 
    h_BRSFH125->Fill(h_BRHiggs->GetXaxis()->GetBinLabel(ibr), brsf) ; 
    std::cout << " H(125): SF( BR(" << h_BRHiggs->GetXaxis()->GetBinLabel(ibr) << ") : " << brsf << std::endl ; 
  }
  fout->Write() ; 
  fout->Close() ; 

  return ; 

}
