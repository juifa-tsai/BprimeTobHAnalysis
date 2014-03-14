{
  TFile *_file0 = TFile::Open("SingleBprime500_8TeV.root") ; 
  _file0->cd() ; 

  gROOT->SetStyle("Plain") ; 
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TCanvas * c0 = new TCanvas("c0", "", 800, 600) ; 
  c0->cd();
  c0->SetLogy() ; 
  h_tch_pt_bpbar->Draw("");
  h_tch_pt_bp->SetLineColor(kRed);
  h_tch_pt_bp->Draw("sames");
  h_tch_pt_bpbar->GetXaxis()->SetTitle("b' p_{t}") ; 

  TLegend* leg0 = new TLegend(0.6,0.7,0.88,0.88) ; 
  leg0->SetBorderSize(0) ;
  leg0->SetFillColor(0) ;
  leg0->AddEntry(h_tch_pt_bp, "t-channel: Single b'", "l") ;
  leg0->AddEntry(h_tch_pt_bpbar, "t-channel: Single #bar{b'}", "l") ;
  leg0->Draw() ; 

  c0->SaveAs("singleBp_pt_tch_lhe.png") ; 
  
  TCanvas * c1 = new TCanvas("c1", "", 800, 600) ; 
  c1->cd();
  h_tch_y_bpbar->Draw("");
  h_tch_y_bp->SetLineColor(kRed);
  h_tch_y_bp->Draw("sames"); 
  h_tch_y_bpbar->GetXaxis()->SetTitle("b' y") ; 

  TLegend* leg1 = new TLegend(0.6,0.7,0.88,0.88) ; 
  leg1->SetBorderSize(0) ;
  leg1->SetFillColor(0) ;
  leg1->AddEntry(h_tch_pt_bp, "t-channel: Single b'", "l") ;
  leg1->AddEntry(h_tch_pt_bpbar, "t-channel: Single #bar{b'}", "l") ;
  leg1->Draw() ; 

  c1->SaveAs("singleBp_y_tch_lhe.png") ; 

  TCanvas * c2 = new TCanvas("c2", "", 800, 600) ; 
  c2->cd();
  c2->SetLogy() ; 
  h_sch_pt_bpbar->Draw("");
  h_sch_pt_bp->SetLineColor(kRed);
  h_sch_pt_bp->Draw("sames");
  h_sch_pt_bpbar->GetXaxis()->SetTitle("b' p_{t}") ; 

  TLegend* leg2 = new TLegend(0.6,0.7,0.88,0.88) ; 
  leg2->SetBorderSize(0) ;
  leg2->SetFillColor(0) ;
  leg2->AddEntry(h_sch_pt_bp, "s-channel: Single b'", "l") ;
  leg2->AddEntry(h_sch_pt_bpbar, "s-channel: Single #bar{b'}", "l") ;
  leg2->Draw() ; 

  c2->SaveAs("singleBp_pt_sch_lhe.png") ; 
  
  TCanvas * c3 = new TCanvas("c3", "", 800, 600) ; 
  c3->cd();
  h_sch_y_bpbar->Draw("");
  h_sch_y_bp->SetLineColor(kRed);
  h_sch_y_bp->Draw("sames"); 
  h_sch_y_bpbar->GetXaxis()->SetTitle("b' y") ; 

  TLegend* leg3 = new TLegend(0.6,0.7,0.88,0.88) ; 
  leg3->SetBorderSize(0) ;
  leg3->SetFillColor(0) ;
  leg3->AddEntry(h_sch_y_bp, "s-channel: Single b'", "l") ;
  leg3->AddEntry(h_sch_y_bpbar, "s-channel: Single #bar{b'}", "l") ;
  leg3->Draw() ; 

  c3->SaveAs("singleBp_y_sch_lhe.png") ; 
}
