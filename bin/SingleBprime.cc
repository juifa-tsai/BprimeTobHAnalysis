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
  fout.str(std::string());
  fout.clear() ;
  fout << fnamecut << ";p_{T} (b') [GeV];" ;
  TH1D* h_pt_bp = new TH1D("h_pt_bp" ,(fout.str()).c_str() ,200 ,0.  ,1000) ; 
  fout.str(std::string());
  fout.clear() ;
  fout << fnamecut << ";y (b');" ;
  TH1D* h_y_bp  = new TH1D("h_y_bp"  ,(fout.str()).c_str() ,100 ,-4. ,4.) ; 
  fout.str(std::string());
  fout.clear() ;

  while ( reader1.readEvent() ) {
    ++ieve1;
    double mbp(0), ptbp(0), ybp(-100) ; 
    TLorentzVector p4_bp ; 
    for(int ii=0;ii< reader1.hepeup.NUP;ii++){
      if ( abs(reader1.hepeup.IDUP[ii]) == 6000007 ) {
        p4_bp.SetPxPyPzE(reader1.hepeup.PUP[ii][0], reader1.hepeup.PUP[ii][1], reader1.hepeup.PUP[ii][2], reader1.hepeup.PUP[ii][3]) ; 
        h_pt_bp->Fill(p4_bp.Pt()) ; 
        h_y_bp -> Fill(p4_bp.Rapidity()) ; 
      }
    }
  }

  rootoutput->Write();
  rootoutput->Close();

  return ; 

}

void readerExample() {
  // Open a stream connected to an event file:
  std::ifstream ifs1("singletBp_ppTobpj_test_MBp600_RL0p5_unweighted_events.lhe");

  std::ifstream ifs2("singletBp_ppTobpt_test_MBp600_RL0p5_unweighted_events.lhe");

  // Create the Reader object:
  LHEF::Reader reader1(ifs1);
  LHEF::Reader reader2(ifs2);

  // Print out the header information:
  std::cout << reader1.headerBlock;

  // Print out the addinional comments in the init block:
  std::cout << reader1.initComments;

  // Print out the beam energies:
  std::cout << "Beam A: " << reader1.heprup.EBMUP.first << " GeV, Beam B: "
    << reader1.heprup.EBMUP.second << " GeV." << std::endl;

  // Now loop over all events:
  long ieve1(0), ieve2(0);

  TFile *rootoutput = new TFile("SingleBprime500_8TeV.root","RECREATE");
  rootoutput->cd();
  TH1D* h_sch_pt_bp    = new TH1D("h_sch_pt_bp", ";p_{T} (b');", 200, 0., 1000) ; 
  TH1D* h_sch_y_bp     = new TH1D("h_sch_y_bp", ";y (b');", 100, -4., 4.) ; 
  TH1D* h_sch_pt_bpbar = new TH1D("h_sch_pt_bpbar", ";p_{T} (b');", 200, 0., 1000) ; 
  TH1D* h_sch_y_bpbar  = new TH1D("h_sch_y_bpbar", ";y (b');", 100, -4., 4.) ; 
  TH1D* h_tch_pt_bp    = new TH1D("h_tch_pt_bp", ";p_{T} (b');", 200, 0., 1000) ; 
  TH1D* h_tch_y_bp     = new TH1D("h_tch_y_bp", ";y (b');", 100, -4., 4.) ; 
  TH1D* h_tch_pt_bpbar = new TH1D("h_tch_pt_bpbar", ";p_{T} (b');", 200, 0., 1000) ; 
  TH1D* h_tch_y_bpbar  = new TH1D("h_tch_y_bpbar", ";y (b');", 100, -4., 4.) ; 

  while ( reader1.readEvent() ) {

    ++ieve1;

    double mbp, ptbp, ybp ; 
    TLorentzVector p4_bp, p4_bpbar ; 
    for(int ii=0;ii< reader1.hepeup.NUP;ii++){
      if ( reader1.hepeup.IDUP[ii] == 7 ) {
        p4_bp.SetPxPyPzE(reader1.hepeup.PUP[ii][0], reader1.hepeup.PUP[ii][1], reader1.hepeup.PUP[ii][2], reader1.hepeup.PUP[ii][3]) ; 
        if ( abs(p4_bp.Mag() - reader1.hepeup.PUP[ii][4])/p4_bp.Mag() > 0.05 ) std::cout << "SingleBprime: WARNING: 4-vector mass " << p4_bp.Mag() << " != on-shell mass " << reader1.hepeup.PUP[ii][4] << std::endl ; 
        h_sch_pt_bp->Fill(p4_bp.Pt()) ; 
        h_sch_y_bp -> Fill(p4_bp.Rapidity()) ; 
      }
      if ( reader1.hepeup.IDUP[ii] == -7 ) {
        p4_bpbar.SetPxPyPzE(reader1.hepeup.PUP[ii][0], reader1.hepeup.PUP[ii][1], reader1.hepeup.PUP[ii][2], reader1.hepeup.PUP[ii][3]) ; 
        if ( abs(p4_bpbar.Mag() - reader1.hepeup.PUP[ii][4])/p4_bpbar.Mag() > 0.05 ) std::cout << "SingleBprime: WARNING: 4-vector mass " << p4_bpbar.Mag() << " != on-shell mass " << reader1.hepeup.PUP[ii][4] << std::endl ; 
        h_sch_pt_bpbar->Fill(p4_bpbar.Pt()) ; 
        h_sch_y_bpbar -> Fill(p4_bpbar.Rapidity()) ; 
      }
    }

  }

  while ( reader2.readEvent() ) {

    ++ieve2;

    double mbp, ptbp, ybp ; 
    TLorentzVector p4_bp, p4_bpbar ; 
    for(int ii=0;ii< reader2.hepeup.NUP;ii++){
      if ( reader2.hepeup.IDUP[ii] == 7 ) {
        p4_bp.SetPxPyPzE(reader2.hepeup.PUP[ii][0], reader2.hepeup.PUP[ii][1], reader2.hepeup.PUP[ii][2], reader2.hepeup.PUP[ii][3]) ; 
        if ( abs(p4_bp.Mag() - reader2.hepeup.PUP[ii][4])/p4_bp.Mag() > 0.05 ) std::cout << "SingleBprime: WARNING: 4-vector mass " << p4_bp.Mag() << " != on-shell mass " << reader2.hepeup.PUP[ii][4] << std::endl ; 
        h_tch_pt_bp->Fill(p4_bp.Pt()) ; 
        h_tch_y_bp -> Fill(p4_bp.Rapidity()) ; 
      }
      if ( reader2.hepeup.IDUP[ii] == -7 ) {
        p4_bpbar.SetPxPyPzE(reader2.hepeup.PUP[ii][0], reader2.hepeup.PUP[ii][1], reader2.hepeup.PUP[ii][2], reader2.hepeup.PUP[ii][3]) ; 
        if ( abs(p4_bpbar.Mag() - reader2.hepeup.PUP[ii][4])/p4_bpbar.Mag() > 0.05 ) std::cout << "SingleBprime: WARNING: 4-vector mass " << p4_bpbar.Mag() << " != on-shell mass " << reader2.hepeup.PUP[ii][4] << std::endl ; 
        h_tch_pt_bpbar->Fill(p4_bpbar.Pt()) ; 
        h_tch_y_bpbar -> Fill(p4_bpbar.Rapidity()) ; 
      }
    }

  }

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


