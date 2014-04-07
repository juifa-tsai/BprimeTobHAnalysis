#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2.h>
#include <THStack.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TString.h>
#include <TStyle.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TGraphAsymmErrors.h>
#include <TCutG.h>

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <stdlib.h>

#include "CMSstyle.C"
#include "help.C"

using namespace std;

TString filename = "/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_13_patch3_Bpbh/src/BpbH/BprimeTobHAnalysis/test/OnLxplus/ABCD_DRbH_SkimmedTrees_31Mar2014/Final_histograms_BprimebH.root" ; 
TString name = "BJetVeto_Htag_HTAK5" ; 
bool doData(true) ; 
int nRebin(4) ; 

Double_t Lint = 19700.0 ; 

TString dir4plots = "ABCD_SkimmedTrees_31Mar2014/BJetVeto" ; 

TString formata = ".pdf";
TString formatb = ".png";
TString formatc = ".C";

void abcdVal () {

  TH2D* hist2_bkg    ; 
  TH2D* hist2_ttjets ;
  TH2D* hist2_qcd    ;
  TH2D* hist2_sig0   ;
  TH2D* hist2_sig1   ;
  TH2D* hist2_sig2   ;
  TH2D* hist2_data   ;

  TFile *myFile  = TFile::Open(filename,"READ") ;
  myFile->cd();

  hist2_qcd    = (TH2D*)myFile->Get("QCD__"+name);
  hist2_ttjets = (TH2D*)myFile->Get("TTJets__"+name);
: hist2_sig0   = (TH2D*)myFile->Get("BprimeBprimeToBHBHinc_M-500__"+name);
  hist2_sig1   = (TH2D*)myFile->Get("BprimeBprimeToBHBHinc_M-800__"+name);
  hist2_sig2   = (TH2D*)myFile->Get("BprimeBprimeToBHBHinc_M-1000__"+name);
  if (doData) hist2_data = (TH2D*)myFile->Get("DATA__"+name);
  hist2_bkg = (TH2D*) hist2_ttjets->Clone("hist2_bkg");
  hist2_bkg->Add(hist2_qcd) ; 

  hist2_bkg->SetMarkerStyle(7) ;
  hist2_bkg->SetMarkerSize(1) ; 
  hist2_bkg->SetMarkerColor(kYellow-8) ; 
  hist2_sig0->SetMarkerStyle(24) ;
  hist2_sig0->SetMarkerSize(1) ; 
  hist2_sig0->SetMarkerColor(kRed) ; 
  hist2_sig1->SetMarkerStyle(25) ;
  hist2_sig1->SetMarkerSize(1) ; 
  hist2_sig1->SetMarkerColor(kGreen) ; 
  hist2_sig2->SetMarkerStyle(26) ;
  hist2_sig2->SetMarkerSize(1) ; 
  hist2_sig2->SetMarkerColor(kBlue) ; 

  TH1D* histxlow_bkg   = hist2_bkg ->ProjectionX("BKG__"+name+"_HT_antiH", 1, 1, "e") ; 
  TH1D* histxlow_sig0  = hist2_sig0->ProjectionX("BpBpbHbH_M500__"+name+"_HT_antiH", 1, 1, "e") ;
  TH1D* histxlow_sig1  = hist2_sig1->ProjectionX("BpBpbHbH_M800__"+name+"_HT_antiH", 1, 1, "e") ;
  TH1D* histxlow_sig2  = hist2_sig2->ProjectionX("BpBpbHbH_M1000__"+name+"_HT_antiH", 1, 1, "e") ;
  TH1D* histxlow_data  = hist2_data->ProjectionX("DATA__"+name+"_HT_antiH", 1, 1, "e") ; 

  TH1D* histxhigh_bkg  = hist2_bkg ->ProjectionX("BKG__"+name+"_HT_H", 2, -1, "e") ; 
  TH1D* histxhigh_sig0 = hist2_sig0->ProjectionX("BpBpbHbH_M500__"+name+"_HT_H", 2, 2, "e") ;
  TH1D* histxhigh_sig1 = hist2_sig1->ProjectionX("BpBpbHbH_M800__"+name+"_HT_H", 2, 2, "e") ;
  TH1D* histxhigh_sig2 = hist2_sig2->ProjectionX("BpBpbHbH_M1000__"+name+"_HT_H", 2, 2, "e") ;
  TH1D* histxhigh_data = hist2_data->ProjectionX("DATA__"+name+"_HT_H", 2, -1, "e") ; 

  double NAdata = histxhigh_data->Integral(76,90) ; 
  double NBdata = histxhigh_data->Integral(91,200) ; 
  double NAbyNBdata = NAdata/NBdata ; 
  double errNAbyNBdata = (2/(NAdata+NBdata))*sqrt((NAdata*NBdata)/(NAdata+NBdata)) ; 

  double NCdata = histxlow_data->Integral(76,90) ; 
  double NDdata = histxlow_data->Integral(91,200) ; 
  double NCbyNDdata = NCdata/NDdata ; 
  double errNCbyNDdata = (2/(NCdata+NDdata))*sqrt((NCdata*NDdata)/(NCdata+NDdata)) ; 

  std::cout << " DATA NA = " << NAdata << std::endl ; 
  std::cout << " DATA NB = " << NBdata << std::endl ; 
  std::cout << " DATA NC = " << NCdata << std::endl ; 
  std::cout << " DATA ND = " << NDdata << std::endl ; 
  std::cout << " DATA NA/NB = " << NAbyNBdata << " +/- " << errNAbyNBdata << std::endl ; 
  std::cout << " DATA NC/ND = " << NCbyNDdata << " +/- " << errNCbyNDdata << std::endl ; 

  double NAbkg = histxhigh_bkg->Integral(76,90) ; 
  double NBbkg = histxhigh_bkg->Integral(91,200) ; 
  double NAbyNBbkg = NAbkg/NBbkg ; 
  double errNAbyNBbkg = (2/histxlow_bkg->GetEntries())*sqrt((NBbkg/histxlow_bkg->Integral())*(1-(NBbkg/histxlow_bkg->Integral()))*histxlow_bkg->GetEntries()) ; 

  double NCbkg = histxlow_bkg->Integral(76,90) ; 
  double NDbkg = histxlow_bkg->Integral(91,200) ; 
  double NCbyNDbkg = NCbkg/NDbkg ; 
  double errNCbyNDbkg = (2/histxlow_bkg->GetEntries())*sqrt((NDbkg/histxlow_bkg->Integral())*(1-(NDbkg/histxlow_bkg->Integral()))*histxlow_bkg->GetEntries()) ; 

  std::cout << " BKG NA = " << NAbkg << std::endl ; 
  std::cout << " BKG NB = " << NBbkg << std::endl ; 
  std::cout << " BKG NC = " << NCbkg << std::endl ; 
  std::cout << " BKG ND = " << NDbkg << std::endl ; 
  std::cout << " BKG NA/NB = " << NAbyNBbkg << " +/- " << errNAbyNBbkg << std::endl ; 
  std::cout << " BKG NC/ND = " << NCbyNDbkg << " +/- " << errNCbyNDbkg << std::endl ; 

  std::cout << " M(b')500 NA = " << histxhigh_sig0->Integral(76,90)  << std::endl ; 
  std::cout << " M(b')500 NB = " << histxhigh_sig0->Integral(91,200) << std::endl ; 
  std::cout << " M(b')500 NC = " << histxlow_sig0->Integral(76,90)  << std::endl ; 
  std::cout << " M(b')500 ND = " << histxlow_sig0->Integral(91,200) << std::endl ; 
  std::cout << " Sig contamination M(b')500 = " << (histxhigh_sig0->Integral(76,90)+histxlow_sig0->Integral(76,90)+histxlow_sig0->Integral(91,200))/(histxhigh_bkg->Integral(76,90)+histxlow_bkg->Integral(76,90)+histxlow_bkg->Integral(91,200)) << std::endl ;  
  std::cout << " Sig leakage M(b')500 = " << (histxhigh_sig0->Integral(76,90)+histxlow_sig0->Integral(76,90)+histxlow_sig0->Integral(91,200))/(histxhigh_sig0->Integral(76,90)+histxhigh_sig0->Integral(91,200)+histxlow_sig0->Integral(76,90)+histxlow_sig0->Integral(91,200)) << std::endl ;  

  TH1D* h_cont_sig0 = new TH1D("h_cont_sig0", "M(b')=500 GeV ;HT(AK5) [GeV]; Signal contamination in sideband (%);", 200, 0., 2000.)  ;
  TH1D* h_leak_sig0 = new TH1D("h_leak_sig0", "M(b')=500 GeV ;HT(AK5) [GeV]; Signal leakage in sideband (%);", 200, 0., 2000.)  ;
  TH1D* h_cont_sig1 = new TH1D("h_cont_sig1", "M(b')=800 GeV ;HT(AK5) [GeV]; Signal contamination in sideband (%);", 200, 0., 2000.)  ;
  TH1D* h_leak_sig1 = new TH1D("h_leak_sig1", "M(b')=800 GeV ;HT(AK5) [GeV]; Signal leakage in sideband (%);", 200, 0., 2000.)  ;
  TH1D* h_cont_sig2 = new TH1D("h_cont_sig2", "M(b')=1000 GeV ;HT(AK5) [GeV]; Signal contamination in sideband (%);", 200, 0., 2000.)  ;
  TH1D* h_leak_sig2 = new TH1D("h_leak_sig2", "M(b')=1000 GeV ;HT(AK5) [GeV]; Signal leakage in sideband (%);", 200, 0., 2000.)  ;
  for (int ii = 0; ii < 200; ++ii) {
    double iht = ii*10. ;
    TCutG* cuthiggs = new TCutG("cuthiggs",4);
    cuthiggs->SetVarX("x") ; 
    cuthiggs->SetVarY("y") ; 
    cuthiggs -> SetPoint(0, 0   , 0.5) ; 
    cuthiggs -> SetPoint(1, 2000., 0.5) ; 
    cuthiggs -> SetPoint(2, 2000., 1.5) ; 
    cuthiggs -> SetPoint(3, 0   , 1.5) ;
    cuthiggs->SetLineColor(kBlue) ;
    cuthiggs->SetLineWidth(2) ; 
    TCutG* cutsig = new TCutG("cutsig",5);
    cutsig->SetVarX("x") ; 
    cutsig->SetVarY("y") ; 
    cutsig -> SetPoint(0, iht  , 0.5) ; 
    cutsig -> SetPoint(1, 2000., 0.5) ; 
    cutsig -> SetPoint(2, 2000., 1.5) ; 
    cutsig -> SetPoint(3, iht  , 1.5) ;
    cutsig -> SetPoint(4, iht  , 0.5) ; 
    cutsig->SetLineColor(kMagenta) ;
    cutsig->SetLineWidth(4) ; 
    cutsig->SetLineStyle(2) ; 
    TH1D *h_bkg_sb  = hist2_bkg ->ProjectionX("BKG__"+name+"_HT_sb", 0, -1, "e[-cutsig]") ; 
    TH1D *h_sig0_sig = hist2_sig0->ProjectionX("BpBpbHbH_M500__"+name+"_HT_sig", 0, -1, "e[cutsig]") ; 
    TH1D *h_sig0_sb = hist2_sig0->ProjectionX("BpBpbHbH_M500__"+name+"_HT_sb",   0, -1, "e[-cutsig]") ; 
    if ( h_bkg_sb->Integral() > 0. ) h_cont_sig0 -> Fill(iht, h_sig0_sb->Integral()/h_bkg_sb->Integral()) ; 
    h_leak_sig0 -> Fill(iht, h_sig0_sb->Integral()/(h_sig0_sig->Integral()+h_sig0_sb->Integral())) ; 
    TH1D *h_sig1_sig = hist2_sig1->ProjectionX("BpBpbHbH_M500__"+name+"_HT_sig", 0, -1, "e[cutsig]") ; 
    TH1D *h_sig1_sb = hist2_sig1->ProjectionX("BpBpbHbH_M500__"+name+"_HT_sb",   0, -1, "e[-cutsig]") ; 
    if ( h_bkg_sb->Integral() > 0. ) h_cont_sig1 -> Fill(iht, h_sig1_sb->Integral()/h_bkg_sb->Integral()) ; 
    h_leak_sig1 -> Fill(iht, h_sig1_sb->Integral()/(h_sig1_sig->Integral()+h_sig1_sb->Integral())) ; 
    TH1D *h_sig2_sig = hist2_sig2->ProjectionX("BpBpbHbH_M500__"+name+"_HT_sig", 0, -1, "e[cutsig]") ; 
    TH1D *h_sig2_sb = hist2_sig2->ProjectionX("BpBpbHbH_M500__"+name+"_HT_sb",   0, -1, "e[-cutsig]") ; 
    if ( h_bkg_sb->Integral() > 0. ) h_cont_sig2 -> Fill(iht, h_sig2_sb->Integral()/h_bkg_sb->Integral()) ; 
    h_leak_sig2 -> Fill(iht, h_sig2_sb->Integral()/(h_sig2_sig->Integral()+h_sig2_sb->Integral())) ; 
  }

  cutsig -> SetPoint(0, 900. , 0.5) ; 
  cutsig -> SetPoint(1, 2000., 0.5) ; 
  cutsig -> SetPoint(2, 2000., 1.5) ; 
  cutsig -> SetPoint(3, 900. , 1.5) ;
  cutsig -> SetPoint(4, 900. , 0.5) ; 

  TCutG* cutall = new TCutG("cutall",5);
  cutall->SetVarX("x") ; 
  cutall->SetVarY("y") ; 
  cutall -> SetPoint(0, 750  ,-0.5) ; 
  cutall -> SetPoint(1, 2000.,-0.5) ; 
  cutall -> SetPoint(2, 2000., 1.5) ; 
  cutall -> SetPoint(3, 750  , 1.5) ;
  cutall -> SetPoint(4, 750  ,-0.5) ; 
  cutall->SetLineColor(kGreen) ;
  cutall->SetLineWidth(3) ; 
  cutall->SetLineStyle(1) ; 
  
  fix(histxlow_bkg )        ; 
  fix(histxlow_sig0)        ; 
  fix(histxlow_sig1)        ; 
  fix(histxlow_sig2)        ; 
  fix(histxlow_data)        ; 

  histxhigh_sig2->Scale(10.) ; 
  fix(histxhigh_bkg )        ; 
  fix(histxhigh_sig0)        ; 
  fix(histxhigh_sig1)        ; 
  fix(histxhigh_sig2)        ; 
  fix(histxhigh_data)        ; 

  histxlow_bkg   -> Rebin(nRebin) ; 
  histxlow_sig0  -> Rebin(nRebin) ; 
  histxlow_sig1  -> Rebin(nRebin) ; 
  histxlow_sig2  -> Rebin(nRebin) ; 
  histxlow_data  -> Rebin(nRebin) ; 
  histxhigh_bkg  -> Rebin(nRebin) ; 
  histxhigh_sig0 -> Rebin(nRebin) ; 
  histxhigh_sig1 -> Rebin(nRebin) ; 
  histxhigh_sig2 -> Rebin(nRebin) ; 
  histxhigh_data -> Rebin(nRebin) ; 

  //gROOT->SetBatch(kTRUE);
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  CMSstyle() ; 

  beautify(h_cont_sig0 ,kRed   ,0, 1) ;  
  beautify(h_leak_sig0 ,kRed   ,0, 1) ;  
  beautify(h_cont_sig1 ,kBlue  ,0, 1) ;  
  beautify(h_leak_sig1 ,kBlue  ,0, 1) ;  
  beautify(h_cont_sig2 ,kGreen ,0, 1) ;  
  beautify(h_leak_sig2 ,kGreen ,0, 1) ;  

  TCanvas * c0 = new TCanvas("HT_Bkg_Sig", "HT_Bkg_Sig", 800, 600) ;
  c0->cd();
  c0->SetLogy(); 
  histxhigh_sig0 -> SetLineColor(kGreen) ; 
  histxlow_sig0  -> SetLineColor(kBlack) ; 
  histxhigh_bkg  -> SetLineColor(kRed) ; 
  histxlow_bkg   -> SetLineColor(kBlue) ; 
  histxhigh_bkg  -> Draw("HIST") ; 
  histxlow_bkg   -> Draw("HISTSAME") ; 
  histxhigh_sig0 -> Draw("HISTSAME") ; 
  histxlow_sig0  -> Draw("HISTSAME") ; 
  histxhigh_bkg  -> SetMaximum(1E+3*histxhigh_bkg->GetMaximum()) ; 

  TLegend* leg0 = new TLegend(.55,0.66,0.88,0.9, "", "brNDC")  ;
  if (name.Contains("BJetVeto")) leg0->SetHeader("B jet veto") ; 
  else if (name.Contains("BJetsSel")) leg0->SetHeader("After b jet sel.") ; 
  leg0->SetBorderSize(0);
  leg0->SetFillColor(0);
  leg0->AddEntry(histxhigh_bkg,"Background: Higgs sel","l") ;
  leg0->AddEntry(histxhigh_sig0,"M(b') = 500 GeV: Higgs sel","l") ;
  leg0->AddEntry(histxlow_bkg,"Background: Anti-Higgs sel","l") ;
  leg0->AddEntry(histxlow_sig0,"M(b') = 500 GeV: Anti-Higgs sel","l") ;
  leg0->Draw(); 

  TCanvas * c4 = new TCanvas("HT_Data", "HT_Data", 800, 600) ;
  c4->cd();
  c4->SetLogy(); 
  histxhigh_data -> SetMarkerStyle(2) ; 
  histxlow_data  -> SetMarkerStyle(3) ; 
  histxhigh_data -> SetMarkerColor(kRed) ; 
  histxlow_data  -> SetMarkerColor(kBlue) ; 
  histxhigh_data -> Draw("E1") ; 
  histxlow_data  -> Draw("E1SAME") ; 
  histxhigh_data -> SetMaximum(1E+3*histxhigh_data->GetMaximum()) ; 

  TLegend* leg0data = new TLegend(.55,0.66,0.88,0.9, "", "brNDC")  ;
  if (name.Contains("BJetVeto")) leg0data->SetHeader("B jet veto") ; 
  else if (name.Contains("BJetsSel")) leg0data->SetHeader("After b jet sel.") ; 
  leg0data->SetBorderSize(0);
  leg0data->SetFillColor(0);
  leg0data->AddEntry(histxhigh_data,"Data: Higgs sel","p") ;
  leg0data->AddEntry(histxlow_data,"Data: Anti-Higgs sel","p") ;
  leg0data->Draw(); 

  TCanvas *cd = new TCanvas("Htag_HT_Data","Htag_HT_Data",800,600);
  cd->cd();
  hist2_data->Draw() ;
  cutall->Draw("same") ; 
  cutsig->Draw("same") ; 

  TLegend* leg = new TLegend(.20,0.20,0.45,0.39, "", "brNDC")  ;
  if (name.Contains("BJetVeto")) leg->SetHeader("B jet veto") ; 
  else if (name.Contains("BJetsSel")) leg->SetHeader("After b jet sel.") ; 
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(hist2_data,"Data","l") ;
  leg->AddEntry(cutall,"ABCD area","l") ;
  leg->AddEntry(cutsig,"Signal region","l") ;
  leg->Draw(); 
  leg->Draw(); 

  TCanvas *c1 = new TCanvas("Htag_HT_Sig_Bkg","Htag_HT_Sig_Bkg",800,600);
  c1->cd();
  hist2_bkg->Draw() ;
  hist2_sig0->Draw("SAME") ; 
  cutall->Draw("same") ; 
  cutsig->Draw("same") ; 

  TLegend* leg = new TLegend(.20,0.20,0.45,0.39, "", "brNDC")  ;
  if (name.Contains("BJetVeto")) leg->SetHeader("B jet veto") ; 
  else if (name.Contains("BJetsSel")) leg->SetHeader("After b jet sel.") ; 
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(hist2_bkg,"Background","l") ;
  leg->AddEntry(hist2_sig0,"M(b') 500 GeV","l") ;
  leg->AddEntry(cutall,"ABCD area","l") ;
  leg->AddEntry(cutsig,"Signal region","l") ;
  leg->Draw(); 
  leg->Draw(); 

  TCanvas *c2 = new TCanvas("SignalContamination","Signal Contamination",800,600);
  c2->cd();
  c2->SetLogy();
  h_cont_sig0->Draw("") ;
  h_cont_sig1->Draw("SAME") ;
  h_cont_sig2->Draw("SAME") ;
  h_cont_sig0->SetMinimum(0.1*h_cont_sig2->GetMinimum()) ; 
  h_cont_sig0->SetMaximum(1E+4*h_cont_sig1->GetMaximum()) ; 

  TCanvas *c3 = new TCanvas("SignalLeakage","Signal Leakage",800,600);
  c3->cd();
  h_leak_sig0->Draw("") ;
  h_leak_sig1->Draw("SAME") ;
  h_leak_sig2->Draw("SAME") ;

  TLegend* leg = new TLegend(.20,0.76,0.45,0.9, "", "brNDC")  ;
  if (name.Contains("BJetVeto")) leg->SetHeader("B jet veto") ; 
  else if (name.Contains("BJetsSel")) leg->SetHeader("After b jet sel.") ; 
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(h_cont_sig0,"M(b') = 500 GeV","l") ;
  leg->AddEntry(h_cont_sig1,"M(b') = 800 GeV","l") ;
  leg->AddEntry(h_cont_sig2,"M(b') = 1000 GeV","l") ;
  leg->Draw(); 
  leg->Draw(); 

  TString action = "mkdir -p " + dir4plots;
  system(action);

  //c0->SaveAs(dir4plots+"/"+c0->GetName()+formata); 
  //c1->SaveAs(dir4plots+"/"+c1->GetName()+formata); 
  //c2->SaveAs(dir4plots+"/"+c2->GetName()+formata); 
  //c3->SaveAs(dir4plots+"/"+c3->GetName()+formata); 
  //cd->SaveAs(dir4plots+"/"+cd->GetName()+formata); 
  //c4->SaveAs(dir4plots+"/"+c4->GetName()+formata); 

  //c0->SaveAs(dir4plots+"/"+c0->GetName()+formatb); 
  //c1->SaveAs(dir4plots+"/"+c1->GetName()+formatb); 
  //c2->SaveAs(dir4plots+"/"+c2->GetName()+formatb); 
  //c3->SaveAs(dir4plots+"/"+c3->GetName()+formatb); 
  //cd->SaveAs(dir4plots+"/"+cd->GetName()+formatb); 
  //c4->SaveAs(dir4plots+"/"+c4->GetName()+formatb); 

  //c0->SaveAs(dir4plots+"/"+c0->GetName()+formatc); 
  //c1->SaveAs(dir4plots+"/"+c1->GetName()+formatc); 
  //c2->SaveAs(dir4plots+"/"+c2->GetName()+formatc); 
  //c3->SaveAs(dir4plots+"/"+c3->GetName()+formatc); 
  //cd->SaveAs(dir4plots+"/"+cd->GetName()+formatc); 
  //c4->SaveAs(dir4plots+"/"+c4->GetName()+formatc); 

}
