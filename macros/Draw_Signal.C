{
  gROOT->ProcessLine(".x LXBATCH_SkimmedJobs_ABCD_25Mar2014/TriggerSel_FatJets_MassPruned_Log.C") ; 
  BpBpbHbH_M500__TriggerSel_FatJets_MassPruned->GetXaxis()->SetRangeUser(0,300) ; 
  BpBpbHbH_M500__TriggerSel_FatJets_MassPruned->GetXaxis()->SetTitle("Pruned mass [GeV]") ; 
  TCanvas* c0 = new TCanvas("BpBpBHBH_PrunedMass", "Fat jet pruned mass", 800, 600) ; 
  c0->cd() ; 
  BpBpbHbH_M500__TriggerSel_FatJets_MassPruned->Draw("HIST") ; 
  BpBpbHbH_M800__TriggerSel_FatJets_MassPruned->Draw("SAMEHIST") ; 
  BpBpbHbH_M1000__TriggerSel_FatJets_MassPruned->Draw("SAMEHIST") ; 
  TLegend* leg0 = new TLegend(0.6,0.6,0.8,0.88,"CA8 jet pruned mass","brNDC") ; 
  leg0->SetBorderSize(0) ; 
  leg0->SetFillColor(0) ; 
  leg0->AddEntry(BpBpbHbH_M500__TriggerSel_FatJets_MassPruned, "M(b') = 500 GeV", "l") ; 
  leg0->AddEntry(BpBpbHbH_M800__TriggerSel_FatJets_MassPruned, "M(b') = 800 GeV", "l") ; 
  leg0->AddEntry(BpBpbHbH_M1000__TriggerSel_FatJets_MassPruned, "M(b') = 1000 GeV", "l") ; 
  leg0->Draw() ; 
  c0->SaveAs(TString(c0->GetName())+".png") ;

  gROOT->ProcessLine(".x LXBATCH_SkimmedJobs_ABCD_25Mar2014/TriggerSel_FatJets_Mass_Log.C") ; 
  BpBpbHbH_M500__TriggerSel_FatJets_Mass->GetXaxis()->SetRangeUser(0,300) ; 
  BpBpbHbH_M500__TriggerSel_FatJets_Mass->GetXaxis()->SetTitle("Mass [GeV]") ; 
  TCanvas* c1 = new TCanvas("BpBpBHBH_Mass", "Fat jet mass", 800, 600) ; 
  c1->cd() ; 
  BpBpbHbH_M500__TriggerSel_FatJets_Mass->Draw("HIST") ; 
  BpBpbHbH_M800__TriggerSel_FatJets_Mass->Draw("SAMEHIST") ; 
  BpBpbHbH_M1000__TriggerSel_FatJets_Mass->Draw("SAMEHIST") ; 
  TLegend* leg1 = new TLegend(0.6,0.6,0.8,0.88,"CA8 jet mass","brNDC") ; 
  leg1->SetBorderSize(0) ; 
  leg1->SetFillColor(0) ; 
  leg1->AddEntry(BpBpbHbH_M500__TriggerSel_FatJets_Mass, "M(b') = 500 GeV", "l") ; 
  leg1->AddEntry(BpBpbHbH_M800__TriggerSel_FatJets_Mass, "M(b') = 800 GeV", "l") ; 
  leg1->AddEntry(BpBpbHbH_M1000__TriggerSel_FatJets_Mass, "M(b') = 1000 GeV", "l") ; 
  leg1->Draw() ; 
  c1->SaveAs(TString(c1->GetName())+".png") ; 

}

