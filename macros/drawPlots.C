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

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <stdlib.h>

#include "CMSstyle.C"
#include "help.C"

using namespace std;

TString filename         = "/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_13_patch3_Bpbh/src/BpbH/BprimeTobHAnalysis/test/OnLxplus/LXBATCH_Skim_27May/Final_histograms_BprimebH.root" ; 
TString filename_JESUp   = "/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_13_patch3_Bpbh/src/BpbH/BprimeTobHAnalysis/test/OnLxplus/LXBATCH_Jobs_18Dec_JetCorr_JESUp/Final_histograms_BprimebH.root" ; 
TString filename_JESDown = "/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_13_patch3_Bpbh/src/BpbH/BprimeTobHAnalysis/test/OnLxplus/LXBATCH_Jobs_18Dec_JetCorr_JESDown/Final_histograms_BprimebH.root" ; 
TString filename_JERUp   = "/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_13_patch3_Bpbh/src/BpbH/BprimeTobHAnalysis/test/OnLxplus/LXBATCH_Jobs_18Dec_JetCorr_JERUp/Final_histograms_BprimebH.root" ; 
TString filename_JERDown = "/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_13_patch3_Bpbh/src/BpbH/BprimeTobHAnalysis/test/OnLxplus/LXBATCH_Jobs_18Dec_JetCorr_JERDown/Final_histograms_BprimebH.root" ; 
TString filename_SFbUp   = "/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_13_patch3_Bpbh/src/BpbH/BprimeTobHAnalysis/test/OnLxplus/LXBATCH_Jobs_18Dec_JetCorr_SFbUp/Final_histograms_BprimebH.root" ; 
TString filename_SFbDown = "/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_13_patch3_Bpbh/src/BpbH/BprimeTobHAnalysis/test/OnLxplus/LXBATCH_Jobs_18Dec_JetCorr_SFbDown/Final_histograms_BprimebH.root" ; 
TString filename_SFlUp   = "/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_13_patch3_Bpbh/src/BpbH/BprimeTobHAnalysis/test/OnLxplus/LXBATCH_Jobs_18Dec_JetCorr_SFlUp/Final_histograms_BprimebH.root" ; 
TString filename_SFlDown = "/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_13_patch3_Bpbh/src/BpbH/BprimeTobHAnalysis/test/OnLxplus/LXBATCH_Jobs_18Dec_JetCorr_SFlDown/Final_histograms_BprimebH.root" ; 

Double_t Lint = 19700.0 ; 
TString title1 = "CMS Preliminary, 19.7/fb at #sqrt{s} = 8 TeV";
TString datacaption = "Data"; 

TString dir4plots = "Plots_27May" ;  

TString formata = ".pdf";
TString formatb = ".png";
TString formatc = ".C";

//// Common switches  
bool web = 0;
bool setSampleName = 1;
bool systUnc = 0 ; 
bool blindSigRegion = 1; 
bool dodata = 0; 

void DrawAll () ; 
void DrawStacked(TString name, TString histotitle, bool log, bool doData, bool fExtNorm=false, int nRebin=1, bool setXRange=false, double rangeXLow=0., double rangeXHigh=0.);

void drawPlots () {

  gROOT->SetBatch(kTRUE);
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  CMSstyle() ; 

  TString action = "mkdir -p " + dir4plots;
  system(action);

  DrawAll () ; 

  return ; 

}

void DrawAll () {

  DrawStacked("h_cutflow" ,"" ,1 ,dodata ,0 ,1 ,1 ,0 ,9); 

  DrawStacked("TriggerSel_nPVtx_NoPUWt" ,"N(PV), No PU weight" ,0 ,dodata ,0 ,1 ,1 ,0 ,50); 
  DrawStacked("TriggerSel_nPVtx_PUWt" ,"N(PV)" ,0 ,dodata ,0 ,1 ,1 ,0 ,50); 

  DrawStacked("VertexSel_nFatJets" , "N(CA8 jets)" ,1 ,dodata ,0 ,1 ,1 ,0 ,10); 
  DrawStacked("VertexSel_nHJets" , "N(Higgs jets)" ,1 ,dodata ,0 ,1 ,1 ,0 ,10); 
  DrawStacked("HiggsJetSel_nAK5Jets" , "N(AK5 jets)" ,1 ,dodata ,0 ,1 ,1 ,0 ,10); 
  DrawStacked("HiggsJetSel_nBJets" , "N(b jets)" ,1 ,dodata ,0 ,1 ,1 ,0 ,10); 
  DrawStacked("BJetsSel_nAK5Jets" , "N(AK5 jets)" ,1 ,dodata ,0 ,1 ,1 ,0 ,10); 
  DrawStacked("BJetsSel_HTAK5" , "H_{T} (AK5 jets) [GeV]" ,1 ,dodata ,0 ,10 ,1 ,0 ,2500); 
  DrawStacked("HTSel_HTAK5" , "H_{T} (AK5 jets) [GeV]" ,1 ,dodata ,0 ,10 ,1 ,0 ,2500); 

 // TString cuts[7] = {"TriggerSel", "HiggsJetSel", "BJetsSel", "HTSel", "HTSel_1H_1b", "HTSel_1H_2b"} ; 

 // for (int ii = 0; ii < 5; ++ii) {
 //   TString cut = cuts[ii] ; 

 //   DrawStacked(cut+"_nJets"        ,"N(AK5 jets)"                  ,1 ,dodata ,0 ,1 ,1 ,0   ,10); 

 //   DrawStacked(cut+"_nAllAK5"      ,"N(All AK5 jets)"              ,1 ,dodata ,0 ,1 ,1 ,0   ,10); 
 //   DrawStacked(cut+"_AK5Jet_1_Pt"  ,"p_{T} (leading AK5 jet)[GeV]" ,1 ,dodata ,0 ,1 ,1 ,0   ,2000); 
 //   DrawStacked(cut+"_AK5Jet_2_Pt"  ,"p_{T} (2nd AK5 jet)[GeV]"     ,1 ,dodata ,0 ,1 ,1 ,0   ,2000); 
 //   DrawStacked(cut+"_AK5Jet_3_Pt"  ,"p_{T} (3rd AK5 jet)[GeV]"     ,1 ,dodata ,0 ,1 ,1 ,0   ,2000); 
 //   DrawStacked(cut+"_AK5Jet_4_Pt"  ,"p_{T} (4th AK5 jet)[GeV]"     ,1 ,dodata ,0 ,1 ,1 ,0   ,2000); 
 //   DrawStacked(cut+"_AK5Jet_1_Eta" ,"#eta (leading AK5 jet)[GeV]"  ,1 ,dodata ,0 ,2 ,1 ,-4. ,4.); 
 //   DrawStacked(cut+"_AK5Jet_2_Eta" ,"#eta (2nd AK5 jet)[GeV]"      ,1 ,dodata ,0 ,1 ,2 ,-4. ,4.); 
 //   DrawStacked(cut+"_AK5Jet_3_Eta" ,"#eta (3rd AK5 jet)[GeV]"      ,1 ,dodata ,0 ,1 ,2 ,-4. ,4.); 
 //   DrawStacked(cut+"_AK5Jet_4_Eta" ,"#eta (4th AK5 jet)[GeV]"      ,1 ,dodata ,0 ,1 ,2 ,-4. ,4.); 

 //   DrawStacked(cut+"_nFatJets"           ,"N(CA8 jets)"                        ,1 ,dodata ,0 ,1 ,1 ,0 ,5); 
 //   DrawStacked(cut+"_FatJets_Pt"         ,"p_{T}(CA8 jets) [GeV]"              ,1 ,dodata ,0 ,4 ,1 ,0 ,1000); 
 //   DrawStacked(cut+"_FatJets_Mass"       ,"M(CA8 jets) [GeV]"                  ,1 ,dodata ,0 ,1 ,1 ,0 ,1000); 
 //   DrawStacked(cut+"_FatJets_MassPruned" ,"pruned mass(CA8 jets) [GeV]"        ,1 ,dodata ,0 ,1 ,1 ,0 ,1000); 
 //   DrawStacked(cut+"_FatJets_tau2ByTau1" ,"#tau_{2}/#tau_{1}(Pruned CA8 jets)" ,1 ,dodata ,0 ,1 ,1 ,0 ,1); 
 //   DrawStacked(cut+"_FatJets_tau3ByTau1" ,"#tau_{3}/#tau_{1}(Pruned CA8 jets)" ,1 ,dodata ,0 ,1 ,1 ,0 ,1); 
 //   DrawStacked(cut+"_FatJets_tau3ByTau2" ,"#tau_{3}/#tau_{2}(Pruned CA8 jets)" ,1 ,dodata ,0 ,1 ,1 ,0 ,1); 

 //   DrawStacked(cut+"_SubJet1_Pt"                 ,"p_{T}(leading subjet) [GeV]"          ,1 ,dodata ,0 ,1 ,1 ,0 ,1000); 
 //   DrawStacked(cut+"_SubJet2_Pt"                 ,"p_{T}(subleading subjet) [GeV]"       ,1 ,dodata ,0 ,1 ,1 ,0 ,1000); 
 //   DrawStacked(cut+"_SubJet1_Mass"               ,"M(leading subjet) [GeV]"              ,1 ,dodata ,0 ,1 ,1 ,0 ,600); 
 //   DrawStacked(cut+"_SubJet2_Mass"               ,"M(subleading subjet) [GeV]"           ,1 ,dodata ,0 ,1 ,1 ,0 ,300); 
 //   DrawStacked(cut+"_SubJet1_CombinedSVBJetTags" ,"CSV discriminator(leading subjet)"    ,1 ,dodata ,0 ,1 ,1 ,0 ,1); 
 //   DrawStacked(cut+"_SubJet2_CombinedSVBJetTags" ,"CSV discriminator(subleading subjet)" ,1 ,dodata ,0 ,1 ,1 ,0 ,1); 

 //   DrawStacked(cut+"_nBJets"   ,"N(b-tagged AK5 jets)"           ,1 ,dodata ,0 ,1 ,1 ,0  ,5); 
 //   DrawStacked(cut+"_BJet_Pt"  ,"p_{T}(b-tagged AK5 jets) [GeV]" ,1 ,dodata ,0 ,1 ,1 ,0  ,1000); 
 //   DrawStacked(cut+"_BJet_Eta" ,"#eta(b-tagged AK5 jets) [GeV]"  ,1 ,dodata ,0 ,1 ,1 ,-3 , 3); 

 //   DrawStacked(cut+"_nHJets"        ,"N(Higgs jets)"            ,1 ,dodata ,0 ,1 ,1  ,0   ,5); 
 //   DrawStacked(cut+"_HiggsJet_Pt"   ,"p_{T} (Higgs jets) [GeV]" ,1 ,dodata ,0 ,4 ,1  ,0   ,2000); 
 //   DrawStacked(cut+"_HiggsJet_Eta"  ,"#eta (Higgs jet)"         ,1 ,dodata ,0 ,4 ,10 , -4., 4.); 
 //   DrawStacked(cut+"_HiggsJet_Mass" ,"M(Higgs jets) [GeV]"      ,1 ,dodata ,0 ,1 ,1  ,0   ,200); 

 //   DrawStacked(cut+"_nHJetsNoMDyCut"             ,"N(Higgs jets)"             ,1 ,dodata ,0 ,1 ,1  ,0   ,5); 
 //   DrawStacked(cut+"_HiggsJetNoMDyCut_Pt"        ,"p_{T} (Higgs jets) [GeV]"  ,1 ,dodata ,0 ,4 ,1  ,0   ,2000); 
 //   DrawStacked(cut+"_HiggsJetNoMDyCut_Eta"       ,"#eta (Higgs jet)"          ,1 ,dodata ,0 ,4 ,10 , -4., 4.); 
 //   DrawStacked(cut+"_HiggsJetNoMDyCut_Mass"      ,"M(Higgs jets) [GeV]"       ,1 ,dodata ,0 ,1 ,1  ,0   ,200); 
 //   DrawStacked(cut+"_HiggsJetNoMDyCut_DRsubjets" ,"#DeltaR_{y,#phi}(subjets)" ,1 ,dodata ,0 ,1 ,5  ,0   ,2); 

 //   DrawStacked(cut+"_HT"                          ,"H_{T} (Higgs + b jets) [GeV]"               ,1 ,dodata ,0 ,5 ,1 ,0 ,4000); 
 //   DrawStacked(cut+"_HTAK5"                       ,"H_{T} (AK5 jets) [GeV]"                     ,1 ,dodata ,0 ,5 ,1 ,0 ,4000); 
 //   DrawStacked(cut+"_HTAK5_leading4"              ,"H_{T} (leading four AK5 jets) [GeV]"        ,1 ,dodata ,0 ,5 ,1 ,0 ,4000); 
 //   DrawStacked(cut+"_HTCA8_leading2_AK5_leading2" ,"H_{T} (leading two AK5 and CA8 jets) [GeV]" ,1 ,dodata ,0 ,5 ,1 ,0 ,4000); 
 //   DrawStacked(cut+"_HTAllAK5"                    ,"H_{T} (All AK5 jets) [GeV]"                 ,1 ,dodata ,0 ,5 ,1 ,0 ,4000); 

 // }

 // DrawStacked("HTSel_Nbprimes"   ,"N (b' candidates)" ,1 ,dodata ,0 ,5 ,1 ,-.5  ,4.5); 
 // DrawStacked("HTSel_bprimePt"   ,"b' p_{T} [GeV]"    ,1 ,dodata ,0 ,5 ,1 ,  0. ,2000); 
 // DrawStacked("HTSel_bprimeMass" ,"b' mass [GeV]"     ,1 ,dodata ,0 ,5 ,1 ,  0. ,2000); 

 // DrawStacked("HTSel_1H_1b_Nbprimes"                    ,"N (b' candidates)"                          ,1 ,dodata ,0 ,5 ,1 ,-.5 ,4.5); 
 // DrawStacked("HTSel_1H_1b_bprimePt"                    ,"b' p_{T} [GeV]"                             ,1 ,dodata ,0 ,5 ,1, 0.  ,2000); 
 // DrawStacked("HTSel_1H_1b_bprimeMass"                  ,"b' mass [GeV]"                              ,1 ,dodata ,0 ,5 ,1, 0.  ,2000); 
 // DrawStacked("HTSel_1H_1b_HT"                          ,"H_{T} (Higgs + b jets) [GeV]"               ,1 ,dodata ,0 ,5 ,1 ,0   ,4000); 
 // DrawStacked("HTSel_1H_1b_HTAK5"                       ,"H_{T} (AK5 jets) [GeV]"                     ,1 ,dodata ,0 ,5 ,1 ,0   ,4000); 
 // DrawStacked("HTSel_1H_1b_HTAK5_leading4"              ,"H_{T} (leading four AK5 jets) [GeV]"        ,1 ,dodata ,0 ,5 ,1 ,0   ,4000); 
 // DrawStacked("HTSel_1H_1b_HTCA8_leading2_AK5_leading2" ,"H_{T} (leading two AK5 and CA8 jets) [GeV]" ,1 ,dodata ,0 ,5 ,1 ,0   ,4000); 
 // DrawStacked("HTSel_1H_1b_HTAllAK5"                    ,"H_{T} (All AK5 jets) [GeV]"                 ,1 ,dodata ,0 ,5 ,1 ,0   ,4000); 


 // DrawStacked("HTSel_1H_2b_Nbprimes"                    ,"N (b' candidates)"                          ,1 ,dodata ,0 ,5 ,1 ,-.5 ,4.5); 
 // DrawStacked("HTSel_1H_2b_bprimePt"                    ,"b' p_{T} [GeV]"                             ,1 ,dodata ,0 ,5 ,1, 0.  ,2000); 
 // DrawStacked("HTSel_1H_2b_bprimeMass"                  ,"b' mass [GeV]"                              ,1 ,dodata ,0 ,5 ,1, 0.  ,2000); 
 // DrawStacked("HTSel_1H_2b_HT"                          ,"H_{T} (Higgs + b jets) [GeV]"               ,1 ,dodata ,0 ,5 ,1 ,0   ,4000); 
 // DrawStacked("HTSel_1H_2b_HTAK5"                       ,"H_{T} (AK5 jets) [GeV]"                     ,1 ,dodata ,0 ,5 ,1 ,0   ,4000); 
 // DrawStacked("HTSel_1H_2b_HTAK5_leading4"              ,"H_{T} (leading four AK5 jets) [GeV]"        ,1 ,dodata ,0 ,5 ,1 ,0   ,4000); 
 // DrawStacked("HTSel_1H_2b_HTCA8_leading2_AK5_leading2" ,"H_{T} (leading two AK5 and CA8 jets) [GeV]" ,1 ,dodata ,0 ,5 ,1 ,0   ,4000); 
 // DrawStacked("HTSel_1H_2b_HTAllAK5"                    ,"H_{T} (All AK5 jets) [GeV]"                 ,1 ,dodata ,0 ,5 ,1 ,0   ,4000); 

  return ; 

}

void DrawStacked(TString name,
    TString histotitle,
    bool log,
    bool doData,
    bool fExtNorm,
    int nRebin,
    bool setXRange,
    double rangeXLow,
    double rangeXHigh) {

  TH1D* hist_bkg    ; 
  TH1D* hist_ttjets ;
  TH1D* hist_qcd    ;
  TH1D* hist_sig0   ;
  TH1D* hist_sig1   ;
  TH1D* hist_sig2   ;
  TH1D* hist_data   ;

  TFile *myFile  = TFile::Open(filename,"READ") ;
  myFile->cd();

  hist_qcd              = (TH1D*)myFile->Get("QCD__"+name);
  hist_ttjets           = (TH1D*)myFile->Get("TTJets__"+name);
  hist_sig0             = (TH1D*)myFile->Get("BprimeBprimeToBHBHinc_M-500__"+name);
  hist_sig1             = (TH1D*)myFile->Get("BprimeBprimeToBHBHinc_M-800__"+name);
  hist_sig2             = (TH1D*)myFile->Get("BprimeBprimeToBHBHinc_M-1000__"+name);
  if (doData) hist_data = (TH1D*)myFile->Get("DATA__"+name);

  hist_sig2->Scale(10.) ; 

  fix(hist_qcd   )           ; 
  fix(hist_ttjets)           ; 
  fix(hist_sig0  )           ; 
  fix(hist_sig1  )           ; 
  fix(hist_sig2  )           ; 
  if (doData) fix(hist_data) ; 

  if (nRebin > 1) {
    hist_ttjets -> Rebin(nRebin) ;
    hist_qcd    -> Rebin(nRebin) ;
    hist_sig0   -> Rebin(nRebin) ;
    hist_sig1   -> Rebin(nRebin) ;
    hist_sig2   -> Rebin(nRebin) ;
    if (doData) hist_data -> Rebin(nRebin) ;
  }

  hist_bkg = (TH1D*) hist_ttjets->Clone("hist_bkg");
  hist_bkg->Add(hist_qcd) ; 

  if (doData &&  name.Contains("nPVtx") ) { 
    double scalef = hist_data->Integral()/hist_bkg->Integral() ; 
    hist_bkg->Scale(scalef) ; 
  }

  //// Syst. unc. ////
  TH1D *hist_bkg_JESUp   , *hist_bkg_JESDown   , *hist_bkg_JERUp   , *hist_bkg_JERDown    ;
  TH1D *hist_ttjets_JESUp, *hist_ttjets_JESDown, *hist_ttjets_JERUp, *hist_ttjets_JERDown ;
  TH1D *hist_qcd_JESUp   , *hist_qcd_JESDown   , *hist_qcd_JERUp   , *hist_qcd_JERDown    ;
  TH1D *hist_sig0_JESUp  , *hist_sig0_JESDown  , *hist_sig0_JERUp  , *hist_sig0_JERDown   ;
  TH1D *hist_sig1_JESUp  , *hist_sig1_JESDown  , *hist_sig1_JERUp  , *hist_sig1_JERDown   ;
  TH1D *hist_sig2_JESUp  , *hist_sig2_JESDown  , *hist_sig2_JERUp  , *hist_sig2_JERDown   ;
  TH1D *hist_data_JESUp  , *hist_data_JESDown  , *hist_data_JERUp  , *hist_data_JERDown   ;

  TH1D *hist_bkg_SFbUp   , *hist_bkg_SFbDown   , *hist_bkg_SFlUp   , *hist_bkg_SFlDown    ;
  TH1D *hist_ttjets_SFbUp, *hist_ttjets_SFbDown, *hist_ttjets_SFlUp, *hist_ttjets_SFlDown ;
  TH1D *hist_qcd_SFbUp   , *hist_qcd_SFbDown   , *hist_qcd_SFlUp   , *hist_qcd_SFlDown    ;
  TH1D *hist_sig0_SFbUp  , *hist_sig0_SFbDown  , *hist_sig0_SFlUp  , *hist_sig0_SFlDown   ;
  TH1D *hist_sig1_SFbUp  , *hist_sig1_SFbDown  , *hist_sig1_SFlUp  , *hist_sig1_SFlDown   ;
  TH1D *hist_sig2_SFbUp  , *hist_sig2_SFbDown  , *hist_sig2_SFlUp  , *hist_sig2_SFlDown   ;
  TH1D *hist_data_SFbUp  , *hist_data_SFbDown  , *hist_data_SFlUp  , *hist_data_SFlDown   ;

  TFile *myFile_JESUp, *myFile_JESDown, *myFile_JERUp, *myFile_JERDown, 
        *myFile_SFbUp, *myFile_SFbDown, *myFile_SFlUp, *myFile_SFlDown ; 

  const int nbins(hist_bkg->GetNbinsX()) ; 
  double binCentre[nbins] ; 
  double binLow[nbins] ; 
  double binHigh[nbins] ; 
  double bkg[nbins] ; 
  double bkgLow[nbins] ; 
  double bkgHigh[nbins] ; 
  double norm[nbins] ; 
  double bkgProportionalUncLow [nbins] ; 
  double bkgProportionalUncHigh[nbins] ; 
  TGraphAsymmErrors* gr_bkg_uncUp_uncDown ; 
  TGraphAsymmErrors* gr_uncUp_uncDown ; 

  if ( systUnc ) {
    myFile_JESUp = TFile::Open(filename_JESUp, "READ") ; 
    myFile_JESUp->cd();
    hist_ttjets_JESUp = (TH1D*)myFile_JESUp->Get("TTJets__"+name);                         
    hist_qcd_JESUp    = (TH1D*)myFile_JESUp->Get("QCD__"+name);                      
    hist_sig0_JESUp   = (TH1D*)myFile_JESUp->Get("BprimeToBHinc_M-500__"+name); 
    hist_sig1_JESUp   = (TH1D*)myFile_JESUp->Get("BprimeToBHinc_M-800__"+name); 
    hist_sig2_JESUp   = (TH1D*)myFile_JESUp->Get("BprimeToBHinc_M-1000__"+name);
    hist_data_JESUp   = (TH1D*)myFile_JESUp->Get("DATA__"+name);                        

    myFile_JERUp = TFile::Open(filename_JERUp, "READ") ; 
    myFile_JERUp->cd();
    hist_ttjets_JERUp = (TH1D*)myFile_JERUp->Get("TTJets__"+name);                         
    hist_qcd_JERUp    = (TH1D*)myFile_JERUp->Get("QCD__"+name);                      
    hist_sig0_JERUp   = (TH1D*)myFile_JERUp->Get("BprimeToBHinc_M-500__"+name); 
    hist_sig1_JERUp   = (TH1D*)myFile_JERUp->Get("BprimeToBHinc_M-800__"+name); 
    hist_sig2_JERUp   = (TH1D*)myFile_JERUp->Get("BprimeToBHinc_M-1000__"+name);
    hist_data_JERUp   = (TH1D*)myFile_JERUp->Get("DATA__"+name);                        

    myFile_SFbUp = TFile::Open(filename_SFbUp, "READ") ; 
    myFile_SFbUp->cd();
    hist_ttjets_SFbUp = (TH1D*)myFile_SFbUp->Get("TTJets__"+name);                         
    hist_qcd_SFbUp    = (TH1D*)myFile_SFbUp->Get("QCD__"+name);                      
    hist_sig0_SFbUp   = (TH1D*)myFile_SFbUp->Get("BprimeToBHinc_M-500__"+name); 
    hist_sig1_SFbUp   = (TH1D*)myFile_SFbUp->Get("BprimeToBHinc_M-800__"+name); 
    hist_sig2_SFbUp   = (TH1D*)myFile_SFbUp->Get("BprimeToBHinc_M-1000__"+name);
    hist_data_SFbUp   = (TH1D*)myFile_SFbUp->Get("DATA__"+name);                        

    myFile_SFlUp = TFile::Open(filename_SFlUp, "READ") ; 
    myFile_SFlUp->cd();
    hist_ttjets_SFlUp = (TH1D*)myFile_SFlUp->Get("TTJets__"+name);                         
    hist_qcd_SFlUp    = (TH1D*)myFile_SFlUp->Get("QCD__"+name);                      
    hist_sig0_SFlUp   = (TH1D*)myFile_SFlUp->Get("BprimeToBHinc_M-500__"+name); 
    hist_sig1_SFlUp   = (TH1D*)myFile_SFlUp->Get("BprimeToBHinc_M-800__"+name); 
    hist_sig2_SFlUp   = (TH1D*)myFile_SFlUp->Get("BprimeToBHinc_M-1000__"+name);
    hist_data_SFlUp   = (TH1D*)myFile_SFlUp->Get("DATA__"+name);                        

    myFile_JESDown = TFile::Open(filename_JESDown, "READ") ; 
    myFile_JESDown->cd();
    hist_ttjets_JESDown = (TH1D*)myFile_JESDown->Get("TTJets__"+name);                         
    hist_qcd_JESDown    = (TH1D*)myFile_JESDown->Get("QCD__"+name);                      
    hist_sig0_JESDown   = (TH1D*)myFile_JESDown->Get("BprimeToBHinc_M-500__"+name); 
    hist_sig1_JESDown   = (TH1D*)myFile_JESDown->Get("BprimeToBHinc_M-800__"+name); 
    hist_sig2_JESDown   = (TH1D*)myFile_JESDown->Get("BprimeToBHinc_M-1000__"+name);
    hist_data_JESDown   = (TH1D*)myFile_JESDown->Get("DATA__"+name);                        

    myFile_JERDown = TFile::Open(filename_JERDown, "READ") ; 
    myFile_JERDown->cd();
    hist_ttjets_JERDown = (TH1D*)myFile_JERDown->Get("TTJets__"+name);                         
    hist_qcd_JERDown    = (TH1D*)myFile_JERDown->Get("QCD__"+name);                      
    hist_sig0_JERDown   = (TH1D*)myFile_JERDown->Get("BprimeToBHinc_M-500__"+name); 
    hist_sig1_JERDown   = (TH1D*)myFile_JERDown->Get("BprimeToBHinc_M-800__"+name); 
    hist_sig2_JERDown   = (TH1D*)myFile_JERDown->Get("BprimeToBHinc_M-1000__"+name);
    hist_data_JERDown   = (TH1D*)myFile_JERDown->Get("DATA__"+name);                        

    myFile_SFbDown = TFile::Open(filename_SFbDown, "READ") ; 
    myFile_SFbDown->cd();
    hist_ttjets_SFbDown = (TH1D*)myFile_SFbDown->Get("TTJets__"+name);                         
    hist_qcd_SFbDown    = (TH1D*)myFile_SFbDown->Get("QCD__"+name);                      
    hist_sig0_SFbDown   = (TH1D*)myFile_SFbDown->Get("BprimeToBHinc_M-500__"+name); 
    hist_sig1_SFbDown   = (TH1D*)myFile_SFbDown->Get("BprimeToBHinc_M-800__"+name); 
    hist_sig2_SFbDown   = (TH1D*)myFile_SFbDown->Get("BprimeToBHinc_M-1000__"+name);
    hist_data_SFbDown   = (TH1D*)myFile_SFbDown->Get("DATA__"+name);                        

    myFile_SFlDown = TFile::Open(filename_SFlDown, "READ") ; 
    myFile_SFlDown->cd();
    hist_ttjets_SFlDown = (TH1D*)myFile_SFlDown->Get("TTJets__"+name);                         
    hist_qcd_SFlDown    = (TH1D*)myFile_SFlDown->Get("QCD__"+name);                      
    hist_sig0_SFlDown   = (TH1D*)myFile_SFlDown->Get("BprimeToBHinc_M-500__"+name); 
    hist_sig1_SFlDown   = (TH1D*)myFile_SFlDown->Get("BprimeToBHinc_M-800__"+name); 
    hist_sig2_SFlDown   = (TH1D*)myFile_SFlDown->Get("BprimeToBHinc_M-1000__"+name);
    hist_data_SFlDown   = (TH1D*)myFile_SFlDown->Get("DATA__"+name);                        

    fix(hist_ttjets_JESUp) ; fix(hist_ttjets_JESDown); fix(hist_ttjets_JERUp); fix(hist_ttjets_JERDown) ;
    fix(hist_qcd_JESUp   ) ; fix(hist_qcd_JESDown   ); fix(hist_qcd_JERUp   ); fix(hist_qcd_JERDown   ) ;
    fix(hist_sig0_JESUp  ) ; fix(hist_sig0_JESDown  ); fix(hist_sig0_JERUp  ); fix(hist_sig0_JERDown  ) ;
    fix(hist_sig1_JESUp  ) ; fix(hist_sig1_JESDown  ); fix(hist_sig1_JERUp  ); fix(hist_sig1_JERDown  ) ;
    fix(hist_sig2_JESUp  ) ; fix(hist_sig2_JESDown  ); fix(hist_sig2_JERUp  ); fix(hist_sig2_JERDown  ) ;
    fix(hist_data_JESUp  ) ; fix(hist_data_JESDown  ); fix(hist_data_JERUp  ); fix(hist_data_JERDown  ) ;

    fix(hist_ttjets_SFbUp) ; fix(hist_ttjets_SFbDown); fix(hist_ttjets_SFlUp); fix(hist_ttjets_SFlDown) ;
    fix(hist_qcd_SFbUp   ) ; fix(hist_qcd_SFbDown   ); fix(hist_qcd_SFlUp   ); fix(hist_qcd_SFlDown   ) ;
    fix(hist_sig0_SFbUp  ) ; fix(hist_sig0_SFbDown  ); fix(hist_sig0_SFlUp  ); fix(hist_sig0_SFlDown  ) ;
    fix(hist_sig1_SFbUp  ) ; fix(hist_sig1_SFbDown  ); fix(hist_sig1_SFlUp  ); fix(hist_sig1_SFlDown  ) ;
    fix(hist_sig2_SFbUp  ) ; fix(hist_sig2_SFbDown  ); fix(hist_sig2_SFlUp  ); fix(hist_sig2_SFlDown  ) ;
    fix(hist_data_SFbUp  ) ; fix(hist_data_SFbDown  ); fix(hist_data_SFlUp  ); fix(hist_data_SFlDown  ) ;

    if (nRebin > 1) {
      hist_ttjets_JESUp->Rebin(nRebin) ; hist_ttjets_JESDown->Rebin(nRebin) ; hist_ttjets_JERUp->Rebin(nRebin) ; hist_ttjets_JERDown->Rebin(nRebin) ;
      hist_qcd_JESUp   ->Rebin(nRebin) ; hist_qcd_JESDown   ->Rebin(nRebin) ; hist_qcd_JERUp   ->Rebin(nRebin) ; hist_qcd_JERDown   ->Rebin(nRebin) ;
      hist_sig0_JESUp  ->Rebin(nRebin) ; hist_sig0_JESDown  ->Rebin(nRebin) ; hist_sig0_JERUp  ->Rebin(nRebin) ; hist_sig0_JERDown  ->Rebin(nRebin) ;
      hist_sig1_JESUp  ->Rebin(nRebin) ; hist_sig1_JESDown  ->Rebin(nRebin) ; hist_sig1_JERUp  ->Rebin(nRebin) ; hist_sig1_JERDown  ->Rebin(nRebin) ;
      hist_sig2_JESUp  ->Rebin(nRebin) ; hist_sig2_JESDown  ->Rebin(nRebin) ; hist_sig2_JERUp  ->Rebin(nRebin) ; hist_sig2_JERDown  ->Rebin(nRebin) ;
      hist_data_JESUp  ->Rebin(nRebin) ; hist_data_JESDown  ->Rebin(nRebin) ; hist_data_JERUp  ->Rebin(nRebin) ; hist_data_JERDown  ->Rebin(nRebin) ;

      hist_ttjets_SFbUp->Rebin(nRebin) ; hist_ttjets_SFbDown->Rebin(nRebin) ; hist_ttjets_SFlUp->Rebin(nRebin) ; hist_ttjets_SFlDown->Rebin(nRebin) ;
      hist_qcd_SFbUp   ->Rebin(nRebin) ; hist_qcd_SFbDown   ->Rebin(nRebin) ; hist_qcd_SFlUp   ->Rebin(nRebin) ; hist_qcd_SFlDown   ->Rebin(nRebin) ;
      hist_sig0_SFbUp  ->Rebin(nRebin) ; hist_sig0_SFbDown  ->Rebin(nRebin) ; hist_sig0_SFlUp  ->Rebin(nRebin) ; hist_sig0_SFlDown  ->Rebin(nRebin) ;
      hist_sig1_SFbUp  ->Rebin(nRebin) ; hist_sig1_SFbDown  ->Rebin(nRebin) ; hist_sig1_SFlUp  ->Rebin(nRebin) ; hist_sig1_SFlDown  ->Rebin(nRebin) ;
      hist_sig2_SFbUp  ->Rebin(nRebin) ; hist_sig2_SFbDown  ->Rebin(nRebin) ; hist_sig2_SFlUp  ->Rebin(nRebin) ; hist_sig2_SFlDown  ->Rebin(nRebin) ;
      hist_data_SFbUp  ->Rebin(nRebin) ; hist_data_SFbDown  ->Rebin(nRebin) ; hist_data_SFlUp  ->Rebin(nRebin) ; hist_data_SFlDown  ->Rebin(nRebin) ;
    }

    hist_bkg_JESUp = (TH1D*) hist_ttjets_JESUp->Clone("hist_bkg_JESUp"); 
    hist_bkg_JESUp -> Add(hist_qcd_JESUp) ; 
    hist_bkg_JERUp = (TH1D*) hist_ttjets_JERUp->Clone("hist_bkg_JERUp"); 
    hist_bkg_JERUp -> Add(hist_qcd_JERUp) ; 
    hist_bkg_SFbUp = (TH1D*) hist_ttjets_SFbUp->Clone("hist_bkg_SFbUp"); 
    hist_bkg_SFbUp -> Add(hist_qcd_SFbUp) ; 
    hist_bkg_SFlUp = (TH1D*) hist_ttjets_SFlUp->Clone("hist_bkg_SFlUp"); 
    hist_bkg_SFlUp -> Add(hist_qcd_SFlUp) ; 
    hist_bkg_JESDown = (TH1D*) hist_ttjets_JESDown->Clone("hist_bkg_JESDown"); 
    hist_bkg_JESDown -> Add(hist_qcd_JESDown) ; 
    hist_bkg_JERDown = (TH1D*) hist_ttjets_JERDown->Clone("hist_bkg_JERDown"); 
    hist_bkg_JERDown -> Add(hist_qcd_JERDown) ; 
    hist_bkg_SFbDown = (TH1D*) hist_ttjets_SFbDown->Clone("hist_bkg_SFbDown"); 
    hist_bkg_SFbDown -> Add(hist_qcd_SFbDown) ; 
    hist_bkg_SFlDown = (TH1D*) hist_ttjets_SFlDown->Clone("hist_bkg_SFlDown"); 
    hist_bkg_SFlDown -> Add(hist_qcd_SFlDown) ; 

    for (int ii = 1; ii <= hist_bkg->GetNbinsX(); ++ii) {
      double statUnc = hist_bkg->GetBinError(ii) ; 
      double jesUp = hist_bkg->GetBinContent(ii) - hist_bkg_JESUp->GetBinContent(ii) ; 
      double jerUp = hist_bkg->GetBinContent(ii) - hist_bkg_JERUp->GetBinContent(ii) ; 
      double sfbUp = hist_bkg->GetBinContent(ii) - hist_bkg_SFbUp->GetBinContent(ii) ; 
      double sflUp = hist_bkg->GetBinContent(ii) - hist_bkg_SFlUp->GetBinContent(ii) ; 
      double uncUp = TMath::Sqrt( (jesUp*jesUp) + (jerUp*jerUp) + (sfbUp*sfbUp) + (sflUp*sflUp) + (statUnc*statUnc) ) ; 
      double jesDown = hist_bkg->GetBinContent(ii) - hist_bkg_JESDown->GetBinContent(ii) ; 
      double jerDown = hist_bkg->GetBinContent(ii) - hist_bkg_JERDown->GetBinContent(ii) ; 
      double sfbDown = hist_bkg->GetBinContent(ii) - hist_bkg_SFbDown->GetBinContent(ii) ; 
      double sflDown = hist_bkg->GetBinContent(ii) - hist_bkg_SFlDown->GetBinContent(ii) ; 
      double uncDown = TMath::Sqrt( (jesDown*jesDown) + (jerDown*jerDown) + (sfbDown*sfbDown) + (sflDown*sflDown) + (statUnc*statUnc) ) ; 

      binCentre[ii-1] = hist_bkg->GetBinCenter(ii) ; 
      binLow[ii-1] = hist_bkg->GetBinWidth(ii)/2.0 ; 
      binHigh[ii-1] = hist_bkg->GetBinWidth(ii)/2.0 ; 

      bkg[ii-1] = hist_bkg->GetBinContent(ii) ; 
      bkgLow[ii-1] =  uncUp ;
      bkgHigh[ii-1] = uncDown ; 

      norm[ii-1] = 1.0 ; 
      bkgProportionalUncLow [ii-1] = bkg[ii-1] > 0 ? bkgLow[ii-1]/bkg[ii-1] : 0. ; 
      bkgProportionalUncHigh[ii-1] = bkg[ii-1] > 0 ? bkgHigh[ii-1]/bkg[ii-1] : 0. ; 
    }
    gr_bkg_uncUp_uncDown = new TGraphAsymmErrors(nbins, binCentre, bkg, binLow, binHigh, bkgLow, bkgHigh) ; 
    gr_bkg_uncUp_uncDown -> SetName("gr_bkg_uncUp_uncDown ") ; 
    gr_uncUp_uncDown = new TGraphAsymmErrors(nbins, binCentre, norm, binLow, binHigh, bkgProportionalUncLow, bkgProportionalUncHigh) ; 
    gr_uncUp_uncDown -> SetName("gr_bkg_uncUp_uncDown ") ; 

    gr_bkg_uncUp_uncDown -> SetLineColor(0) ;
    gr_bkg_uncUp_uncDown -> SetFillColor(9) ;
    gr_bkg_uncUp_uncDown -> SetFillStyle(3005) ; 

    gr_uncUp_uncDown -> SetLineColor(5) ;
    gr_uncUp_uncDown -> SetFillColor(5) ;
    gr_uncUp_uncDown -> SetFillStyle(1001) ; 

    if ( name.Contains("HTSel_HT") ) {
      double jesUp = (hist_bkg->Integral() - hist_bkg_JESUp->Integral())/hist_bkg->Integral() ; 
      double jerUp = (hist_bkg->Integral() - hist_bkg_JERUp->Integral())/hist_bkg->Integral() ; 
      double sfbUp = (hist_bkg->Integral() - hist_bkg_SFbUp->Integral())/hist_bkg->Integral() ; 
      double sflUp = (hist_bkg->Integral() - hist_bkg_SFlUp->Integral())/hist_bkg->Integral() ; 

      double jesDown = (hist_bkg->Integral() - hist_bkg_JESDown->Integral())/hist_bkg->Integral() ; 
      double jerDown = (hist_bkg->Integral() - hist_bkg_JERDown->Integral())/hist_bkg->Integral() ; 
      double sfbDown = (hist_bkg->Integral() - hist_bkg_SFbDown->Integral())/hist_bkg->Integral() ; 
      double sflDown = (hist_bkg->Integral() - hist_bkg_SFlDown->Integral())/hist_bkg->Integral() ; 

      double btagUp = TMath::Sqrt( (sfbUp*sfbUp) + (sflUp*sflUp) ) ;
      double btagDown = TMath::Sqrt( (sfbDown*sfbDown) + (sflDown*sflDown) ) ;

      std::cout << "*******************************BACKGROUND*******************************\n" ; 
      std::cout << "***" << " JES uncert + = "      << fabs(jesUp) *100. << "% - = " << fabs(jesDown) *100. << "% ***\n" ; 
      std::cout << "***" << " JER uncert + = "      << fabs(jerUp) *100. << "% - = " << fabs(jerDown) *100. << "% ***\n" ; 
      std::cout << "***" << " Btagging uncert + = " << fabs(btagUp)*100. << "% - = " << fabs(btagDown)*100. << "% ***\n" ; 
      std::cout << "************************************************************************\n" ; 

      jesUp = (hist_sig0->Integral() - hist_sig0_JESUp->Integral())/hist_sig0->Integral() ; 
      jerUp = (hist_sig0->Integral() - hist_sig0_JERUp->Integral())/hist_sig0->Integral() ; 
      sfbUp = (hist_sig0->Integral() - hist_sig0_SFbUp->Integral())/hist_sig0->Integral() ; 
      sflUp = (hist_sig0->Integral() - hist_sig0_SFlUp->Integral())/hist_sig0->Integral() ; 

      jesDown = (hist_sig0->Integral() - hist_sig0_JESDown->Integral())/hist_sig0->Integral() ; 
      jerDown = (hist_sig0->Integral() - hist_sig0_JERDown->Integral())/hist_sig0->Integral() ; 
      sfbDown = (hist_sig0->Integral() - hist_sig0_SFbDown->Integral())/hist_sig0->Integral() ; 
      sflDown = (hist_sig0->Integral() - hist_sig0_SFlDown->Integral())/hist_sig0->Integral() ; 

      btagUp = TMath::Sqrt( (sfbUp*sfbUp) + (sflUp*sflUp) ) ;
      btagDown = TMath::Sqrt( (sfbDown*sfbDown) + (sflDown*sflDown) ) ;

      std::cout << "*********************************SIGNAL Mass 500 GeV************************\n" ; 
      std::cout << "***" << " JES uncert + = "      << fabs(jesUp) *100. << "% - = " << fabs(jesDown) *100. << "% ***\n" ; 
      std::cout << "***" << " JER uncert + = "      << fabs(jerUp) *100. << "% - = " << fabs(jerDown) *100. << "% ***\n" ; 
      std::cout << "***" << " Btagging uncert + = " << fabs(btagUp)*100. << "% - = " << fabs(btagDown)*100. << "% ***\n" ; 
      std::cout << "************************************************************************\n" ; 

      jesUp = (hist_sig1->Integral() - hist_sig1_JESUp->Integral())/hist_sig1->Integral() ; 
      jerUp = (hist_sig1->Integral() - hist_sig1_JERUp->Integral())/hist_sig1->Integral() ; 
      sfbUp = (hist_sig1->Integral() - hist_sig1_SFbUp->Integral())/hist_sig1->Integral() ; 
      sflUp = (hist_sig1->Integral() - hist_sig1_SFlUp->Integral())/hist_sig1->Integral() ; 

      jesDown = (hist_sig1->Integral() - hist_sig1_JESDown->Integral())/hist_sig1->Integral() ; 
      jerDown = (hist_sig1->Integral() - hist_sig1_JERDown->Integral())/hist_sig1->Integral() ; 
      sfbDown = (hist_sig1->Integral() - hist_sig1_SFbDown->Integral())/hist_sig1->Integral() ; 
      sflDown = (hist_sig1->Integral() - hist_sig1_SFlDown->Integral())/hist_sig1->Integral() ; 

      btagUp = TMath::Sqrt( (sfbUp*sfbUp) + (sflUp*sflUp) ) ;
      btagDown = TMath::Sqrt( (sfbDown*sfbDown) + (sflDown*sflDown) ) ;

      std::cout << "*********************************SIGNAL Mass 800 GeV************************\n" ; 
      std::cout << "***" << " JES uncert + = "      << fabs(jesUp) *100. << "% - = " << fabs(jesDown) *100. << "% ***\n" ; 
      std::cout << "***" << " JER uncert + = "      << fabs(jerUp) *100. << "% - = " << fabs(jerDown) *100. << "% ***\n" ; 
      std::cout << "***" << " Btagging uncert + = " << fabs(btagUp)*100. << "% - = " << fabs(btagDown)*100. << "% ***\n" ; 
      std::cout << "************************************************************************\n" ; 

      jesUp = (hist_sig2->Integral() - hist_sig2_JESUp->Integral())/hist_sig2->Integral() ; 
      jerUp = (hist_sig2->Integral() - hist_sig2_JERUp->Integral())/hist_sig2->Integral() ; 
      sfbUp = (hist_sig2->Integral() - hist_sig2_SFbUp->Integral())/hist_sig2->Integral() ; 
      sflUp = (hist_sig2->Integral() - hist_sig2_SFlUp->Integral())/hist_sig2->Integral() ; 

      jesDown = (hist_sig2->Integral() - hist_sig2_JESDown->Integral())/hist_sig2->Integral() ; 
      jerDown = (hist_sig2->Integral() - hist_sig2_JERDown->Integral())/hist_sig2->Integral() ; 
      sfbDown = (hist_sig2->Integral() - hist_sig2_SFbDown->Integral())/hist_sig2->Integral() ; 
      sflDown = (hist_sig2->Integral() - hist_sig2_SFlDown->Integral())/hist_sig2->Integral() ; 

      btagUp = TMath::Sqrt( (sfbUp*sfbUp) + (sflUp*sflUp) ) ;
      btagDown = TMath::Sqrt( (sfbDown*sfbDown) + (sflDown*sflDown) ) ;

      std::cout << "*********************************SIGNAL Mass 1000 GeV***********************\n" ; 
      std::cout << "***" << " JES uncert + = "      << fabs(jesUp) *100. << "% - = " << fabs(jesDown) *100. << "% ***\n" ; 
      std::cout << "***" << " JER uncert + = "      << fabs(jerUp) *100. << "% - = " << fabs(jerDown) *100. << "% ***\n" ; 
      std::cout << "***" << " Btagging uncert + = " << fabs(btagUp)*100. << "% - = " << fabs(btagDown)*100. << "% ***\n" ; 
      std::cout << "************************************************************************\n" ; 
    }

  }
  //// Syst. unc. ////

  beautify(hist_qcd   ,42 ,1001 ,1 ,2) ; 
  beautify(hist_ttjets,38 ,1001 ,1 ,2) ; 
  beautify(hist_bkg   ,0  ,0    ,0 ,2) ; 
  beautify(hist_sig0  ,4  ,0    ,2 ,2) ; 
  beautify(hist_sig1  ,6  ,0    ,4 ,4) ; 
  beautify(hist_sig2  ,2  ,0    ,1 ,3) ; 
  if ( name.Contains("nPVtx") || name.Contains("h_cutflow") ) {
    hist_bkg->SetFillStyle(1001);
    hist_bkg->SetFillColor(25);
  }
  else {
    hist_bkg->SetFillStyle(3254);
    hist_bkg->SetFillColor(12);
  }
  if (doData) {
    beautify(hist_data ,1 ,0 ,1) ; 
    hist_data->SetMarkerStyle(20);
    hist_data->SetMarkerSize(0.75);
    hist_data->SetLineWidth(2);
  }

  THStack *stack = new THStack("stack","");
  stack->Add(hist_ttjets) ;
  stack->Add(hist_qcd) ; 

  TH1D *hist_mcUnc, *hist_ratio, *hist_ratioSigBlind ; 
  hist_mcUnc = (TH1D*)hist_bkg->Clone("hist_mcUnc") ; 
  hist_mcUnc->Sumw2() ; 
  hist_mcUnc->SetTitle("Stat. uncertainty on MC yield");
  hist_mcUnc->Divide(hist_bkg);

  if (doData) {
    hist_ratio = (TH1D*) hist_data->Clone("hist_ratio");
    hist_ratio->Sumw2() ; 
    hist_ratio->SetTitle("Data/MC");
    hist_ratio->Divide(hist_bkg);

    if ( blindSigRegion && name.Contains("h_cutflow")) {
      hist_ratioSigBlind = (TH1D*)hist_ratio->Clone("hist_ratioSigBlind") ; 
      hist_ratioSigBlind->SetTitle("Data/MC");
      hist_ratioSigBlind->Reset() ; 
      for (int ii = 1; ii <= hist_data->GetNbinsX(); ++ii) {
        if ( !TString(hist_ratio->GetXaxis()->GetBinLabel(ii)).EqualTo("HTSel") ) { 
          hist_ratioSigBlind->SetBinContent(ii, hist_ratio->GetBinContent(ii)) ; 
          hist_ratioSigBlind->SetBinError(ii, hist_ratio->GetBinError(ii)) ; 
        } 
      }
    }
  }

  TCanvas* c1 = new TCanvas();
  c1->cd();
  c1->SetBorderMode(0);
  c1->SetFrameBorderMode(0);

  TPad* pad0 = new TPad("pad0", "",0,0.30,1,.92) ; 
  pad0->Draw();
  pad0->cd();
  pad0->SetBorderMode(0);
  pad0->SetFrameBorderMode(0);
  beautifyTopPad(pad0) ; 
  pad0->SetLogy(log) ; 

  if (!log) {
    hist_bkg->SetMaximum( doData ? hist_data->GetMaximum()*2.5 : hist_bkg->GetMaximum()*2.5) ;
    hist_bkg->SetMinimum(0.) ; 
  }
  else {
    if (name.Contains("tau") || name.Contains("nFatJets") 
        || name.Contains("nJets") || name.Contains("TriggerSel_FatJets_Pt")
        || name.Contains("TriggerSel_SubJet1_Pt")) 
      hist_bkg->SetMaximum( doData ? hist_data->GetMaximum()*50000000 : hist_bkg->GetMaximum()*50000000) ;
    else 
      hist_bkg->SetMaximum( doData ? hist_data->GetMaximum()*5000 : hist_bkg->GetMaximum()*5000) ;
    hist_bkg->SetMinimum(0.1) ; 
  }

  TAxis* ax = hist_bkg->GetXaxis() ; 
  TAxis* ay = hist_bkg->GetYaxis() ; 
  beautifyAxis(ax) ; 
  beautifyAxis(ay) ; 
  ax->SetTitle(histotitle);
  if (doData) {
    ax->SetLabelSize( 0.0 );
    ax->SetTitleSize( 0.0 );
  }
  ay->SetTitle("Entries");
  hist_bkg->SetTitleOffset(0.83,"Y");

  if (setXRange) {
    if (rangeXLow == rangeXHigh) std::cout << "Error: X-axis low and high ranges have same value\n" ;
    else {
      hist_bkg->GetXaxis()->SetRangeUser(rangeXLow, rangeXHigh) ;
    }
  }

  if ( name.Contains("nPVtx") || name.Contains("h_cutflow") ) {
    hist_bkg->Draw("HIST");
    TH1D* hist_dataSigBlind ; 
    if (doData) {
      if (blindSigRegion)  {
        hist_dataSigBlind = (TH1D*)hist_data->Clone("hist_dataSigBlind") ; 
        hist_dataSigBlind->Reset() ; 
        for (int ii = 1; ii <= hist_data->GetNbinsX(); ++ii) {
          if ( !TString(hist_data->GetXaxis()->GetBinLabel(ii)).EqualTo("HTSel") ) { 
            hist_dataSigBlind->SetBinContent(ii, hist_data->GetBinContent(ii)) ; 
            hist_dataSigBlind->SetBinError(ii, hist_data->GetBinError(ii)) ; 
          } 
        }
        hist_dataSigBlind->Draw("SAMEE"); 
      }
      else  hist_data->Draw("SAMEE"); 
    }
    if ( name.Contains("h_cutflow") ) {
      hist_sig0->Draw("HISTSAME") ; 
      hist_sig1->Draw("HISTSAME") ; 
      hist_sig2->Draw("HISTSAME") ; 
      for (int ii = 1; ii <= hist_bkg->GetNbinsX(); ++ii) {
        if (doData) std::cout << " Data: " << hist_data->GetXaxis()->GetBinLabel(ii) << " Event: " << hist_data->GetBinContent(ii) << " +/- " << hist_data->GetBinError(ii) << std::endl ; 
        std::cout << " Bkg: " << hist_bkg->GetXaxis()->GetBinLabel(ii) << " Event: " << hist_bkg->GetBinContent(ii) << " +/- " << hist_bkg->GetBinError(ii) << std::endl ; 
        std::cout << " Sig0: " << hist_sig0->GetXaxis()->GetBinLabel(ii) << " Event: " << hist_sig0->GetBinContent(ii) << " +/- " << hist_sig0->GetBinError(ii) << std::endl ; 
        std::cout << " Sig1: " << hist_sig1->GetXaxis()->GetBinLabel(ii) << " Event: " << hist_sig1->GetBinContent(ii) << " +/- " << hist_sig1->GetBinError(ii) << std::endl ; 
        std::cout << " Sig2: " << hist_sig2->GetXaxis()->GetBinLabel(ii) << " Event: " << hist_sig2->GetBinContent(ii) << " +/- " << hist_sig2->GetBinError(ii) << std::endl ; 
      }
    }
  }
  else {
    hist_bkg->Draw("hist");
    stack->Draw("histSAME");
    if ( systUnc ) gr_bkg_uncUp_uncDown -> Draw("2Z") ; 
    else hist_bkg->Draw("samee2");
    if (doData) hist_data->Draw("SAMEE1");
    hist_sig0->Draw("HISTSAME") ; 
    hist_sig1->Draw("HISTSAME") ; 
    hist_sig2->Draw("HISTSAME") ; 
  }

  pad0->RedrawAxis();

  TPad *overlay ;
  if (name.Contains("nPVtx")) {
    overlay = new TPad("overlay","",0,0,1,1);
    overlay->SetFillStyle(4000);
    overlay->SetFillColor(0);
    overlay->SetFrameFillStyle(4000);
    overlay->Draw();
    overlay->cd();
    overlay->SetLogy();
    overlay->RedrawAxis() ;
  }

  int move_legend=0;
  if (name.Contains("HiggsJet_Mass")) move_legend = 1 ;
  TLegend *leg ;
  if (move_legend==1) {
    leg =  new TLegend(0.1,0.53,0.40,.92,NULL,"brNDC");
  }
  else {
    if (name.Contains("nJets")) {
      leg = new TLegend(0.20,0.72,0.895,0.93,NULL,"brNDC");
      leg->SetNColumns(2) ; 
    }
    else if (name.Contains("HTSel")) 
      leg = new TLegend(0.56,0.47,0.895,0.93,NULL,"brNDC");
    else 
      leg = new TLegend(0.56,0.53,0.895,0.93,NULL,"brNDC");
  }
  leg->SetBorderSize(1);
  leg->SetTextFont(132);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.06);

  if ( name.Contains("nPVtx") || name.Contains("h_cutflow") ){
    leg->AddEntry(hist_bkg      , "All backgrounds"        , "f") ; 
  }
  //if ( name.Contains("h_cutflow") || name.Contains("HTSel") ) {
  leg -> AddEntry(hist_sig0 ,"b'(500 GeV)"             ,"l") ; 
  leg -> AddEntry(hist_sig1 ,"b'(800 GeV)"             ,"l") ; 
  leg -> AddEntry(hist_sig2 ,"b'(1000 GeV)#times 10"   ,"l") ; 
  //}
  if ( !name.Contains("nPVtx") && !name.Contains("h_cutflow") ) {
    leg->AddEntry(hist_ttjets   , "t#bar{t}+jets"          ,"f");
    leg->AddEntry(hist_qcd      , "Non-t#bar{t} multijets" ,"f");
    if ( systUnc ) leg->AddEntry(gr_bkg_uncUp_uncDown, "Total bkg. error", "f") ; 
    else leg->AddEntry(hist_bkg      , "Bkg. error (stat.)"     ,"f");
  }
  if (doData) leg->AddEntry(hist_data, datacaption              ,"pl");

  leg->Draw();

  pad0->Modified();

  c1->cd();

  TPad* pad1 = new TPad("pad1", "",0,0,1,.25) ; 
  pad1->Draw() ;  
  pad1->cd() ;  
  pad1->SetGridx() ; 
  pad1->SetGridy() ; 

  if (dodata) {
    hist_ratio->SetMarkerStyle(20);
    hist_ratio->SetMarkerSize(0.75);
    hist_ratio->SetLineWidth(2);

    if ( blindSigRegion && name.Contains("h_cutflow")) {
      hist_ratioSigBlind->SetMarkerStyle(20);
      hist_ratioSigBlind->SetMarkerSize(0.75);
      hist_ratioSigBlind->SetLineWidth(2);
    }
  }

  hist_mcUnc->SetMarkerStyle(0);
  hist_mcUnc->SetMarkerSize(0);
  hist_mcUnc->SetLineWidth(0);
  hist_mcUnc->SetFillStyle(1001);
  hist_mcUnc->SetFillColor(kYellow);

  if (dodata) hist_mcUnc->GetYaxis()->SetTitle("Data/MC");
  else hist_mcUnc->GetYaxis()->SetTitle("MC uncert.");
  hist_mcUnc->SetTitleOffset(0.9,"X");
  hist_mcUnc->SetTitleOffset(0.31,"Y");
  hist_mcUnc->GetXaxis()->SetTitle(histotitle);
  hist_mcUnc->GetYaxis()->SetNdivisions( 505 );

  TAxis* ax1 = hist_mcUnc->GetXaxis();
  TAxis* ay1 = hist_mcUnc->GetYaxis();

  beautifyBottomPad(pad1,ax1,ay1) ; 

  if (setXRange) {
    if (rangeXLow == rangeXHigh) std::cout << "Error: X-axis low and high ranges have same value\n" ;
    else {
      hist_mcUnc->GetXaxis()->SetRangeUser(rangeXLow, rangeXHigh) ;
    }
  }

  hist_mcUnc->SetMinimum(0.0);
  hist_mcUnc->SetMaximum(2.6);
  hist_mcUnc->Draw("E2");
  if ( systUnc ) gr_uncUp_uncDown->Draw("2Z") ; 
  if (dodata) {
    if ( blindSigRegion && name.Contains("h_cutflow") ) hist_ratioSigBlind->Draw("SAMEE1");
    else hist_ratio->Draw("SAMEE1"); 
  } 

  pad1->Modified();

  c1->cd();

  char temp[100];
  TPaveText *plotlabel0 = new TPaveText(0.56,0.925,0.90,.95,"NDC");
  plotlabel0->SetTextColor(kBlack);
  plotlabel0->SetFillColor(kWhite);
  plotlabel0->SetBorderSize(0);
  plotlabel0->SetTextAlign(12);
  plotlabel0->SetTextSize(0.045);
  sprintf(temp, "%.1f", Lint/1000);
  plotlabel0->AddText( (string("L = ") + temp + string("/fb, ") + string("#sqrt{s} = 8 TeV")).c_str() ) ; 
  plotlabel0->Draw() ; 

  c1->Modified();
  c1->cd();
  c1->SetSelected(c1) ;

  TString name_plot=name+"_Linear"+formata;
  if(log) name_plot=name+"_Log"+formata;
  c1->SaveAs(dir4plots+"/"+name_plot);
  name_plot=name+"_Linear"+formatb;
  if(log) name_plot=name+"_Log"+formatb;
  c1->SaveAs(dir4plots+"/"+name_plot);
  name_plot=name+"_Linear"+formatc;
  if(log) name_plot=name+"_Log"+formatc;
  c1->SaveAs(dir4plots+"/"+name_plot);

  if (log && web) {  
    pad0 ->cd();
    pad0->SetLogy(false);
    c1->cd();
    c1->SaveAs(dir4plots+"/"+name+"_Linear"+formata);
  }

  return ; 

}

