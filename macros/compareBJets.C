#include <cmath> 
#include <algorithm>

void compare (TString& hist) {

  TFile* f0 = TFile::Open("/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_13_patch3_Bpbh/src/BpbH/BprimeTobHAnalysis/test/bprimeTobHSkimmed_BpBp1000_bHiggsOverlapRemoved.root","READ") ; 
  TFile* f1 = TFile::Open("/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_13_patch3_Bpbh/src/BpbH/BprimeTobHAnalysis/test/bprimeTobHSkimmed_BpBp1000_bHiggsHighMassOverlapRemoved.root","READ") ; 
  TFile* f2 = TFile::Open("/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_13_patch3_Bpbh/src/BpbH/BprimeTobHAnalysis/test/bprimeTobHSkimmed_BpBp1000_bHiggsAntiHiggsOverlapRemoved.root","READ") ; 

  TH1D* h0 = (TH1D*)f0->Get("BprimebH/"+hist) ; 
  TH1D* h1 = (TH1D*)f1->Get("BprimebH/"+hist) ; 
  TH1D* h2 = (TH1D*)f2->Get("BprimebH/"+hist) ; 

  //gStyle->SetOptStat(1111) ; 
  gStyle->SetOptTitle(0) ; 

  h0->SetLineColor(kRed) ; 
  h1->SetLineColor(kBlue) ; 
  h2->SetLineColor(kGreen) ; 
  h1->SetLineWidth(2) ; 
  h1->SetLineStyle(2) ; 
  h2->SetLineWidth(2) ; 
  h2->SetLineStyle(3) ; 

  h0->SetMaximum(std::max(std::max(h0->GetMaximum(), h1->GetMaximum()), h2->GetMaximum())*1.5) ; 

  TCanvas* c1000 = new TCanvas("c1000"+hist, "c1000"+hist, 800, 600) ;
  c1000->cd();
  h0->Draw("HIST") ;
  c1000->Modified() ; c1000->Update(); 
  TPaveStats *st0 = (TPaveStats*)c1000->GetPrimitive("stats") ; 
  st0->SetName("st0") ; 
  h1->Draw("HISTSAMES") ;
  c1000->Modified() ; c1000->Update(); 
  TPaveStats *st1 = (TPaveStats*)c1000->GetPrimitive("stats") ; 
  st1->SetName("st1") ; 
  h2->Draw("HISTSAMES") ;
  c1000->Modified() ; c1000->Update(); 
  TPaveStats *st2 = (TPaveStats*)c1000->GetPrimitive("stats") ; 
  st2->SetName("st2") ; 

  st1->SetY2NDC(st0->GetY1NDC()) ; 
  st1->SetY1NDC(st1->GetY2NDC() - (st0->GetY2NDC() - st0->GetY1NDC()) ) ; 
  st2->SetY2NDC(st1->GetY1NDC()) ; 
  st2->SetY1NDC(st2->GetY2NDC() - (st1->GetY2NDC() - st1->GetY1NDC()) ) ; 
  st0->SetTextColor(kRed) ;
  st1->SetTextColor(kBlue) ; 
  st2->SetTextColor(kGreen) ; 

  TLegend * leg = new TLegend(0.12,0.75,0.78,0.88,"M(b') = 1000 GeV","brNDC") ; 
  leg->SetBorderSize(0) ;
  leg->SetFillColor(0);
  leg->AddEntry(h0, hist+": No H M(H) = [0,140] GeV overlap", "l") ;
  leg->AddEntry(h1, hist+": No H M(H) = [90, 140] GeV overlap", "l") ;
  leg->AddEntry(h2, hist+": No H M(H) = [0,75] GeV or anti-H overlap", "l") ;
  leg->Draw() ; 

  c1000->SaveAs(TString("Overlap_H_AntiH_BJets/")+c1000->GetName()+TString(".png")) ; 
  c1000->SaveAs(TString("Overlap_H_AntiH_BJets/")+c1000->GetName()+TString(".pdf")) ; 

}

void compareBJets () {

  compare("FatJetSel_BJet_Pt") ; 
  compare("FatJetSel_BJet_Eta") ; 
  compare("FatJetSel_nBJets") ; 

  compare("BJetsSel_BJet_Pt") ; 
  compare("BJetsSel_BJet_Eta") ; 
  compare("BJetsSel_HTAK5") ; 

}

