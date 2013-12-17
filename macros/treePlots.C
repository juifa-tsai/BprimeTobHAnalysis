#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TString.h>
#include <TCanvas.h> 
#include <TLegend.h> 

#include <iostream>

void treePlots (TString file="/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_13_patch2_Bpbh/src/BpbH/BprimeTobHAnalysis/test/evtSkim.root") {

  TFile* fin = TFile::Open(file) ;
  TTree* tree = (TTree*)fin->Get("EvtSkim/tree") ; 

  tree->Project("jetpt", "JetInfo.Pt") ;
  tree->Project("jetpt_JESUp", "JetInfo_JESUp.Pt") ;
  tree->Project("jetpt_JESDown", "JetInfo_JESDown.Pt") ;

  tree->Project("jetpt_JER", "JetInfo_JER.Pt") ;
  tree->Project("jetpt_JERUp", "JetInfo_JERUp.Pt") ;
  tree->Project("jetpt_JERDown", "JetInfo_JERDown.Pt") ;

  tree->Project("jetpt_BJetInfo", "BJetInfo.Pt") ;
  tree->Project("jetpt_BJetInfo_SFbUp", "BJetInfo_SFbUp.Pt") ;
  tree->Project("jetpt_BJetInfo_SFbDown", "BJetInfo_SFbDown.Pt") ;
  tree->Project("jetpt_BJetInfo_SFlUp", "BJetInfo_SFlUp.Pt") ;
  tree->Project("jetpt_BJetInfo_SFlDown", "BJetInfo_SFlDown.Pt") ;

  TH1F* jetpt         = (TH1F*)gDirectory->Get("jetpt") ; 
  TH1F* jetpt_JESUp   = (TH1F*)gDirectory->Get("jetpt_JESUp") ;
  TH1F* jetpt_JESDown = (TH1F*)gDirectory->Get("jetpt_JESDown") ;

  TH1F* jetpt_JER     = (TH1F*)gDirectory->Get("jetpt_JER") ;
  TH1F* jetpt_JERUp   = (TH1F*)gDirectory->Get("jetpt_JERUp") ;
  TH1F* jetpt_JERDown = (TH1F*)gDirectory->Get("jetpt_JERDown") ;

  TH1F* jetpt_BJetInfo         = (TH1F*)gDirectory->Get("jetpt_BJetInfo") ;
  TH1F* jetpt_BJetInfo_SFbUp   = (TH1F*)gDirectory->Get("jetpt_BJetInfo_SFbUp") ;
  TH1F* jetpt_BJetInfo_SFbDown = (TH1F*)gDirectory->Get("jetpt_BJetInfo_SFbDown") ;
  TH1F* jetpt_BJetInfo_SFlUp   = (TH1F*)gDirectory->Get("jetpt_BJetInfo_SFlUp") ;
  TH1F* jetpt_BJetInfo_SFlDown = (TH1F*)gDirectory->Get("jetpt_BJetInfo_SFlDown") ;

  jetpt         -> SetLineWidth(2) ; 

  jetpt_JESUp   -> SetLineWidth(2) ; 
  jetpt_JESDown -> SetLineWidth(2) ; 
  jetpt_JESUp   -> SetLineWidth(2) ; 
  jetpt_JESDown -> SetLineWidth(2) ; 

  jetpt_JER     -> SetLineWidth(2) ; 
  jetpt_JERUp   -> SetLineWidth(2) ; 
  jetpt_JERDown -> SetLineWidth(2) ; 
  jetpt_JERUp   -> SetLineWidth(2) ; 
  jetpt_JERDown -> SetLineWidth(2) ; 

  jetpt_BJetInfo         ->SetLineWidth(2) ; 
  jetpt_BJetInfo_SFbUp   ->SetLineWidth(2) ; 
  jetpt_BJetInfo_SFbDown ->SetLineWidth(2) ; 
  jetpt_BJetInfo_SFbUp   ->SetLineWidth(2) ; 
  jetpt_BJetInfo_SFbDown ->SetLineWidth(2) ; 
  jetpt_BJetInfo_SFlUp   ->SetLineWidth(2) ; 
  jetpt_BJetInfo_SFlDown ->SetLineWidth(2) ; 
  jetpt_BJetInfo_SFlUp   ->SetLineWidth(2) ; 
  jetpt_BJetInfo_SFlDown ->SetLineWidth(2) ; 

  jetpt         ->SetLineColor(kBlack) ; 

  jetpt_JESUp   ->SetLineColor(kBlue) ; 
  jetpt_JESDown ->SetLineColor(kBlue) ; 
  jetpt_JESUp   ->SetLineStyle(2) ; 
  jetpt_JESDown ->SetLineStyle(3) ; 

  jetpt_JER     ->SetLineColor(kRed) ; 
  jetpt_JERUp   ->SetLineColor(kRed-2) ; 
  jetpt_JERDown ->SetLineColor(kRed+2) ; 
  jetpt_JERUp   ->SetLineStyle(2) ; 
  jetpt_JERDown ->SetLineStyle(3) ; 

  jetpt_BJetInfo         ->SetLineColor(kMagenta) ; 
  jetpt_BJetInfo_SFbUp   ->SetLineColor(kMagenta-2) ; 
  jetpt_BJetInfo_SFbDown ->SetLineColor(kMagenta+2) ; 
  jetpt_BJetInfo_SFbUp   ->SetLineStyle(2) ; 
  jetpt_BJetInfo_SFbDown ->SetLineStyle(3) ; 
  jetpt_BJetInfo_SFlUp   ->SetLineColor(kGreen) ; 
  jetpt_BJetInfo_SFlDown ->SetLineColor(kGreen) ; 
  jetpt_BJetInfo_SFlUp   ->SetLineStyle(2) ; 
  jetpt_BJetInfo_SFlDown ->SetLineStyle(3) ; 

  TCanvas *c0 = new TCanvas("c0_JetInfo") ; 
  c0->cd(); 
  jetpt          -> Draw() ; 
  jetpt_JESUp    ->Draw("SAMES") ; 
  jetpt_JESDown  ->Draw("SAMES") ;
  jetpt_JER      ->Draw("SAMES") ;
  jetpt_JERUp    ->Draw("SAMES") ;
  jetpt_JERDown  ->Draw("SAMES") ;

  TLegend* tleg0 = new TLegend(0.6,0.6, 0.9, 0.9, NULL, "brNDC") ;
  tleg0->SetBorderSize(1);
  tleg0->SetTextFont(132);
  tleg0->SetLineColor(1);
  tleg0->SetLineStyle(1);
  tleg0->SetLineWidth(1);
  tleg0->SetFillColor(0);
  tleg0->SetFillStyle(1001);
  tleg0->SetBorderSize(0);
  tleg0->SetTextSize(0.06);
  tleg0 -> AddEntry(jetpt         ,"jetpt        " , "l") ;  
  tleg0 -> AddEntry(jetpt_JESUp   ,"jetpt_JESUp  " , "l") ;  
  tleg0 -> AddEntry(jetpt_JESDown ,"jetpt_JESDown" , "l") ;  
  tleg0 -> AddEntry(jetpt_JER     ,"jetpt_JER    " , "l") ;  
  tleg0 -> AddEntry(jetpt_JERUp   ,"jetpt_JERUp  " , "l") ;  
  tleg0 -> AddEntry(jetpt_JERDown ,"jetpt_JERDown" , "l") ;  
  tleg0->Draw();

  std::cout << " jetpt         Entries = " <<  jetpt          -> GetEntries() << " | integral = "  << jetpt          -> Integral() << std::endl ; 
  std::cout << " jetpt_JESUp   Entries = " <<  jetpt_JESUp    -> GetEntries() << " | integral = "  << jetpt_JESUp    -> Integral() << std::endl ; 
  std::cout << " jetpt_JESDown Entries = " <<  jetpt_JESDown  -> GetEntries() << " | integral = "  << jetpt_JESDown  -> Integral() << std::endl ; 
  std::cout << " jetpt_JER     Entries = " <<  jetpt_JER      -> GetEntries() << " | integral = "  << jetpt_JER      -> Integral() << std::endl ; 
  std::cout << " jetpt_JERUp   Entries = " <<  jetpt_JERUp    -> GetEntries() << " | integral = "  << jetpt_JERUp    -> Integral() << std::endl ; 
  std::cout << " jetpt_JERDown Entries = " <<  jetpt_JERDown  -> GetEntries() << " | integral = "  << jetpt_JERDown  -> Integral() << std::endl ; 

  TCanvas *c1 = new TCanvas("c1_BJetInfo") ; 
  c1->cd(); 
  jetpt_BJetInfo         -> Draw () ; 
  jetpt_BJetInfo_SFbUp   -> Draw ("SAMES") ; 
  jetpt_BJetInfo_SFbDown -> Draw ("SAMES") ; 
  jetpt_BJetInfo_SFbUp   -> Draw ("SAMES") ; 
  jetpt_BJetInfo_SFbDown -> Draw ("SAMES") ; 
  jetpt_BJetInfo_SFlUp   -> Draw ("SAMES") ; 
  jetpt_BJetInfo_SFlDown -> Draw ("SAMES") ; 

  TLegend* tleg1 = new TLegend(0.6,0.6, 0.9, 0.9, NULL, "brNDC") ;
  tleg1->SetBorderSize(1);
  tleg1->SetTextFont(132);
  tleg1->SetLineColor(1);
  tleg1->SetLineStyle(1);
  tleg1->SetLineWidth(1);
  tleg1->SetFillColor(0);
  tleg1->SetFillStyle(1001);
  tleg1->SetBorderSize(0);
  tleg1->SetTextSize(0.06);
  tleg1 -> AddEntry(jetpt_BJetInfo         ,"jetpt_BJetInfo        " , "l") ;  
  tleg1 -> AddEntry(jetpt_BJetInfo_SFbUp   ,"jetpt_BJetInfo_SFbUp  " , "l") ;  
  tleg1 -> AddEntry(jetpt_BJetInfo_SFbDown ,"jetpt_BJetInfo_SFbDown" , "l") ;  
  tleg1 -> AddEntry(jetpt_BJetInfo_SFbUp   ,"jetpt_BJetInfo_SFbUp  " , "l") ;  
  tleg1 -> AddEntry(jetpt_BJetInfo_SFbDown ,"jetpt_BJetInfo_SFbDown" , "l") ;  
  tleg1 -> AddEntry(jetpt_BJetInfo_SFlUp   ,"jetpt_BJetInfo_SFlUp  " , "l") ;  
  tleg1 -> AddEntry(jetpt_BJetInfo_SFlDown ,"jetpt_BJetInfo_SFlDown" , "l") ;  
  tleg1->Draw();

  return ; 

}
