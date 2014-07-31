#include <cmath>
#include <iostream>
#include <TFile.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <cstring>
#include <sstream>

#include "BpbH/BprimeTobHAnalysis/interface/LHEF.h"

#define pow2(a) pow(a,2.)

using namespace std;

void readerExample(std::string& fname) {

  std::ifstream ifs1(fname.c_str());
  LHEF::Reader reader1(ifs1);
  std::cout << reader1.headerBlock;
  std::cout << reader1.initComments;
  std::cout << "Beam A: " << reader1.heprup.EBMUP.first << " GeV, Beam B: "
    << reader1.heprup.EBMUP.second << " GeV." << std::endl;
  long ieve1(0) ;  
  unsigned pos = fname.find("_unweighted_events.lhe");  
  std::string fnamecut = fname.substr (0, pos);
  std::stringstream fout ;
  fout << fnamecut << "_8TeV.root"  ; 
  TFile *rootoutput = new TFile((fout.str()).c_str(),"RECREATE"); 
  rootoutput->cd();
  fout.str(std::string());

  fout.clear() ;
  fout << fnamecut << ";p_{T} (b');Events;"         ; 
  TH1D* h_pt_bp      = new TH1D("h_pt_bp"      , (fout.str()).c_str(), 200, 0. , 1000) ; 
  fout.clear() ;
  fout << fnamecut << ";y (b');Events;"             ; 
  TH1D* h_y_bp       = new TH1D("h_y_bp"       , (fout.str()).c_str(), 100, -4., 4.) ; 
  fout.clear() ;
  fout << fnamecut << ";p_{T} (#bar{b'});Events;"   ; 
  TH1D* h_pt_bpbar   = new TH1D("h_pt_bpbar"   , (fout.str()).c_str(), 200, 0. , 1000) ; 
  fout.clear() ;
  fout << fnamecut << ";y (#bar{b'});Events;"       ; 
  TH1D* h_y_bpbar    = new TH1D("h_y_bpbar"    , (fout.str()).c_str(), 100, -4., 4.) ; 
  fout.clear() ;
  fout << fnamecut << ";p_{T} (b'+#bar{b'});Events;"; 
  TH1D* h_pt_bpbpbar = new TH1D("h_pt_bpbpbar" , (fout.str()).c_str(), 200, 0. , 1000) ; 
  fout.clear() ;
  fout << fnamecut << ";y (b'+#bar{b'});Events;"    ; 
  TH1D* h_y_bpbpbar  = new TH1D("h_y_bpbpbar"  , (fout.str()).c_str(), 100, -4., 4.) ; 
  fout.clear() ;
  fout << fnamecut << ";p_{T} (b from b');Events;"  ; 
  TH1D* h_pt_b       = new TH1D("h_pt_b"       , (fout.str()).c_str(), 200, 0. , 1000) ; 
  fout.clear() ;
  fout << fnamecut << ";y (b from b');Events;"      ; 
  TH1D* h_y_b        = new TH1D("h_y_b"        , (fout.str()).c_str(), 100, -4., 4.) ; 
  fout.clear() ;
  fout << fnamecut << ";p_{T} (H from b');Events;"  ; 
  TH1D* h_pt_h       = new TH1D("h_pt_h"       , (fout.str()).c_str(), 200, 0. , 1000) ; 
  fout.clear() ;
  fout << fnamecut << ";y (H from b');Events;"      ; 
  TH1D* h_y_h        = new TH1D("h_y_h"        , (fout.str()).c_str(), 100, -4., 4.) ; 
  fout.clear() ;
  fout << fnamecut << "; #DeltaR(b,#bar{b});Events;"; 
  TH1D* h_dr_bb      = new TH1D("h_dr_bb"      , (fout.str()).c_str(), 50 , 0. , 4.) ; 
  fout.clear() ;
  fout << fnamecut << "; H_{T} [GeV];Events;"       ; 
  TH1D* h_htall      = new TH1D("h_htall"      , (fout.str()).c_str(), 25 , 0. , 2500.) ; 
  fout.clear() ;
  fout << fnamecut << "; H_{T} [GeV];Events;"       ; 
  TH1D* h_htsel      = new TH1D("h_htsel"      , (fout.str()).c_str(), 25 , 0. , 2500.) ; 
  fout.clear() ;
  fout << fnamecut << "; H_{T} [GeV];Efficiency;"       ; 
  TH1D* h_eff_htsel  = new TH1D("h_eff_htsel"      , (fout.str()).c_str(), 25 , 0. , 2500.) ;

  while ( reader1.readEvent() ) {

    ++ieve1;

    double mbp, ptbp, ybp, ptb, yb, pth, yh ; 
    TLorentzVector p4_bp, p4_bpbar, p4_b, p4_h, p4_bFromH, p4_bbarFromH ; 
    for(int ii=0;ii< reader1.hepeup.NUP;ii++){
      if ( reader1.hepeup.IDUP[ii] == 6000007 ) {
        p4_bp.SetPxPyPzE(reader1.hepeup.PUP[ii][0], reader1.hepeup.PUP[ii][1], reader1.hepeup.PUP[ii][2], reader1.hepeup.PUP[ii][3]) ; 
        if ( abs(p4_bp.Mag() - reader1.hepeup.PUP[ii][4])/p4_bp.Mag() > 0.05 ) std::cout << "SingleBprime: WARNING: 4-vector mass " << p4_bp.Mag() << " != on-shell mass " << reader1.hepeup.PUP[ii][4] << std::endl ; 
        h_pt_bp->Fill(p4_bp.Pt()) ; 
        h_y_bp -> Fill(p4_bp.Rapidity()) ; 
        h_pt_bpbpbar->Fill(p4_bp.Pt()) ; 
        h_y_bpbpbar -> Fill(p4_bp.Rapidity()) ; 
      }
      if ( reader1.hepeup.IDUP[ii] == -6000007 ) {
        p4_bpbar.SetPxPyPzE(reader1.hepeup.PUP[ii][0], reader1.hepeup.PUP[ii][1], reader1.hepeup.PUP[ii][2], reader1.hepeup.PUP[ii][3]) ; 
        if ( abs(p4_bpbar.Mag() - reader1.hepeup.PUP[ii][4])/p4_bpbar.Mag() > 0.05 ) std::cout << "SingleBprime: WARNING: 4-vector mass " << p4_bpbar.Mag() << " != on-shell mass " << reader1.hepeup.PUP[ii][4] << std::endl ; 
        h_pt_bpbar->Fill(p4_bpbar.Pt()) ; 
        h_y_bpbar -> Fill(p4_bpbar.Rapidity()) ; 
        h_pt_bpbpbar->Fill(p4_bpbar.Pt()) ; 
        h_y_bpbpbar -> Fill(p4_bpbar.Rapidity()) ; 
      }
      if ( abs(reader1.hepeup.IDUP[ii]) == 25 ) {
        p4_h.SetPxPyPzE(reader1.hepeup.PUP[ii][0], reader1.hepeup.PUP[ii][1], reader1.hepeup.PUP[ii][2], reader1.hepeup.PUP[ii][3]) ; 
        h_pt_h -> Fill(p4_h.Pt()) ; 
        h_y_h -> Fill(p4_h.Rapidity()) ; 
      }
      if ( abs(reader1.hepeup.IDUP[ii]) ==  5 ) {
        if ( abs(reader1.hepeup.IDUP[(reader1.hepeup.MOTHUP[ii].first)]) == 6000007 ) {
          p4_b.SetPxPyPzE(reader1.hepeup.PUP[ii][0], reader1.hepeup.PUP[ii][1], reader1.hepeup.PUP[ii][2], reader1.hepeup.PUP[ii][3]) ;
          h_pt_b -> Fill(p4_b.Pt()) ; 
          h_y_b -> Fill(p4_b.Rapidity()) ; 
        }
        else if ( abs(reader1.hepeup.IDUP[(reader1.hepeup.MOTHUP[ii].first)]) == 25 ) {
          if (reader1.hepeup.IDUP[ii] == 5) 
            p4_bFromH.SetPxPyPzE(reader1.hepeup.PUP[ii][0], reader1.hepeup.PUP[ii][1], reader1.hepeup.PUP[ii][2], reader1.hepeup.PUP[ii][3]) ;
          else if ( reader1.hepeup.IDUP[ii] == -5 )
            p4_bbarFromH.SetPxPyPzE(reader1.hepeup.PUP[ii][0], reader1.hepeup.PUP[ii][1], reader1.hepeup.PUP[ii][2], reader1.hepeup.PUP[ii][3]) ;
        }
      }
    }

    double HTall(-1), HTsel(-1) ; 
    HTall = p4_h.Pt() + p4_b.Pt() ; 
    if ( p4_h.Pt() > 300 && abs(p4_h.Rapidity()) < 2.4 
        && p4_b.Pt() > 80 && abs(p4_b.Rapidity()) < 2.4 
        && abs(p4_bFromH.Rapidity()) < 2.4 && abs(p4_bbarFromH.Rapidity()) < 2.4
        && p4_bFromH.DeltaR(p4_bbarFromH) < 0.8 ) HTsel = p4_h.Pt() + p4_b.Pt() ; 
    h_htall -> Fill(HTall) ; 
    h_htsel -> Fill(HTsel) ; 
    h_dr_bb -> Fill(p4_bFromH.DeltaR(p4_bbarFromH)) ; 

  }

  h_eff_htsel -> Divide(h_htsel, h_htall, 1, 1, "b") ; 

  rootoutput->Write();
  rootoutput->Close();

  // Now we are done.

  return ; 
}

int main(int argc, char *argv[]) { 

  if ( argc != 2 ) // argc should be 2 for correct execution
    // We print argv[0] assuming it is the program name
    std::cout << "usage: " << argv[0] << " <filename>\n" ; 
  else {
    ifstream infile ( argv[1] );
    if ( !infile.is_open() ) {
      std::cout<<"Could not open file\n"; 
    }
    else {
      std::string fname;
      while ( std::getline ( infile, fname ) && fname[0] != '#' ) { 
        readerExample(fname); 
      } 
    }
  }

  return 0;

}
