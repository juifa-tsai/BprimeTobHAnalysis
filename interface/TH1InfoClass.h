#ifndef TH1INFO_H
#define TH1INFO_H

#include <map>
#include <string>
#include <TFile.h>

using namespace std;

enum TH1List {
  ABCDana_CutFlow_, 		
  ABCDval_CutFlow_, 		

  ABCDana_CutRegion_, 	
  ABCDana_CutRegion_1b, 
  ABCDana_CutRegion_2b, 
  ABCDval_CutRegion_, 	
  ABCDval_CutRegion_0ak5, 	
  ABCDval_CutRegion_1ak5, 	
  ABCDval_CutRegion_2ak5, 	

  ABCDana_Sumw2_A,
  ABCDana_Sumw2_B,
  ABCDana_Sumw2_C,
  ABCDana_Sumw2_D,
  ABCDana_Sumw2_A_1b,
  ABCDana_Sumw2_B_1b,
  ABCDana_Sumw2_C_1b,
  ABCDana_Sumw2_D_1b,
  ABCDana_Sumw2_A_2b,
  ABCDana_Sumw2_B_2b,
  ABCDana_Sumw2_C_2b,
  ABCDana_Sumw2_D_2b,
  ABCDval_Sumw2_A,
  ABCDval_Sumw2_B,
  ABCDval_Sumw2_C,
  ABCDval_Sumw2_D,

  ABCDana_NAK5_BeforeHTAK5,
  ABCDana_NAK5,
  ABCDana_NAK5_ABCD,
  ABCDana_NAK5_A,
  ABCDana_NAK5_B,
  ABCDana_NAK5_C,
  ABCDana_NAK5_D,
  ABCDana_Numbjet,
  ABCDana_Numbjet_ABCD,
  ABCDana_Numbjet_A,
  ABCDana_Numbjet_B,
  ABCDana_Numbjet_C,
  ABCDana_Numbjet_D,
  ABCDana_NumCA8,
  ABCDana_NumCA8_ABCD,
  ABCDana_NumCA8_A,
  ABCDana_NumCA8_B,
  ABCDana_NumCA8_C,
  ABCDana_NumCA8_D,
  ABCDana_bpMass_A,
  ABCDana_bpMass_B,
  ABCDana_bpMass_C,
  ABCDana_bpMass_D,
  ABCDana_bpPt_A,
  ABCDana_bpPt_B,
  ABCDana_bpPt_C,
  ABCDana_bpPt_D,
  ABCDana_bpEta_A,
  ABCDana_bpEta_B,
  ABCDana_bpEta_C,
  ABCDana_bpEta_D,
  ABCDana_HT_,  		
  ABCDana_HT_A,  		
  ABCDana_HT_B,  		
  ABCDana_HT_C,  		
  ABCDana_HT_D,  		
  ABCDana_HiggsMass_A, 
  ABCDana_HiggsMass_B, 
  ABCDana_HiggsMass_C, 
  ABCDana_HiggsMass_D, 
  ABCDana_CA8Pt_A, 
  ABCDana_CA8Pt_B, 
  ABCDana_CA8Pt_C, 
  ABCDana_CA8Pt_D, 
  ABCDana_dRSubJets_A, 
  ABCDana_dRSubJets_B, 
  ABCDana_dRSubJets_C, 
  ABCDana_dRSubJets_D, 
  ABCDana_Tau2ByTau1_A, 
  ABCDana_Tau2ByTau1_B, 
  ABCDana_Tau2ByTau1_C, 
  ABCDana_Tau2ByTau1_D, 
  ABCDana_Sub1Pt_A, 
  ABCDana_Sub1Pt_B, 
  ABCDana_Sub1Pt_C, 
  ABCDana_Sub1Pt_D, 
  ABCDana_Sub2Pt_A, 
  ABCDana_Sub2Pt_B, 
  ABCDana_Sub2Pt_C, 
  ABCDana_Sub2Pt_D, 
  ABCDana_Sub1Mass_A, 
  ABCDana_Sub1Mass_B, 
  ABCDana_Sub1Mass_C, 
  ABCDana_Sub1Mass_D, 
  ABCDana_Sub2Mass_A, 
  ABCDana_Sub2Mass_B, 
  ABCDana_Sub2Mass_C, 
  ABCDana_Sub2Mass_D, 
  ABCDana_Sub1CSV_A, 	
  ABCDana_Sub1CSV_B, 	
  ABCDana_Sub1CSV_C, 	
  ABCDana_Sub1CSV_D, 	
  ABCDana_Sub2CSV_A, 	
  ABCDana_Sub2CSV_B, 	
  ABCDana_Sub2CSV_C, 	
  ABCDana_Sub2CSV_D, 	
  ABCDana_JetFlavor_A, 
  ABCDana_JetFlavor_B, 
  ABCDana_JetFlavor_C, 
  ABCDana_JetFlavor_D, 
  ABCDana_HiggsMass_,  
  ABCDana_JetFlavor_,  
  ABCDana_LeadCA8_Pt_, 

  ABCDana_1b_bpMass_A,
  ABCDana_1b_bpMass_B,
  ABCDana_1b_bpMass_C,
  ABCDana_1b_bpMass_D,
  ABCDana_1b_bpPt_A,
  ABCDana_1b_bpPt_B,
  ABCDana_1b_bpPt_C,
  ABCDana_1b_bpPt_D,
  ABCDana_1b_bpEta_A,
  ABCDana_1b_bpEta_B,
  ABCDana_1b_bpEta_C,
  ABCDana_1b_bpEta_D,
  ABCDana_1b_HT_A, 
  ABCDana_1b_HT_B, 
  ABCDana_1b_HT_C, 
  ABCDana_1b_HT_D, 
  ABCDana_2b_bpMass_A,
  ABCDana_2b_bpMass_B,
  ABCDana_2b_bpMass_C,
  ABCDana_2b_bpMass_D,
  ABCDana_2b_bpPt_A,
  ABCDana_2b_bpPt_B,
  ABCDana_2b_bpPt_C,
  ABCDana_2b_bpPt_D,
  ABCDana_2b_bpEta_A,
  ABCDana_2b_bpEta_B,
  ABCDana_2b_bpEta_C,
  ABCDana_2b_bpEta_D,
  ABCDana_2b_HT_A, 
  ABCDana_2b_HT_B, 
  ABCDana_2b_HT_C, 
  ABCDana_2b_HT_D, 

  ABCDval_NumCA8_A,
  ABCDval_NumCA8_B,
  ABCDval_NumCA8_C,
  ABCDval_NumCA8_D,
  ABCDval_HT_,  
  ABCDval_HT_A, 
  ABCDval_HT_B, 
  ABCDval_HT_C, 
  ABCDval_HT_D, 
  ABCDval_HiggsMass_A, 
  ABCDval_HiggsMass_B, 
  ABCDval_HiggsMass_C, 
  ABCDval_HiggsMass_D, 
  ABCDval_CA8Pt_A, 
  ABCDval_CA8Pt_B, 
  ABCDval_CA8Pt_C, 
  ABCDval_CA8Pt_D, 
  ABCDval_dRSubJets_A, 
  ABCDval_dRSubJets_B, 
  ABCDval_dRSubJets_C, 
  ABCDval_dRSubJets_D, 
  ABCDval_Tau2ByTau1_A, 
  ABCDval_Tau2ByTau1_B, 
  ABCDval_Tau2ByTau1_C, 
  ABCDval_Tau2ByTau1_D, 
  ABCDval_Sub1Pt_A, 
  ABCDval_Sub1Pt_B, 
  ABCDval_Sub1Pt_C, 
  ABCDval_Sub1Pt_D, 
  ABCDval_Sub2Pt_A, 
  ABCDval_Sub2Pt_B, 
  ABCDval_Sub2Pt_C, 
  ABCDval_Sub2Pt_D, 
  ABCDval_Sub1Mass_A, 
  ABCDval_Sub1Mass_B, 
  ABCDval_Sub1Mass_C, 
  ABCDval_Sub1Mass_D, 
  ABCDval_Sub2Mass_A, 
  ABCDval_Sub2Mass_B, 
  ABCDval_Sub2Mass_C, 
  ABCDval_Sub2Mass_D, 
  ABCDval_Sub1CSV_A, 
  ABCDval_Sub1CSV_B, 
  ABCDval_Sub1CSV_C, 
  ABCDval_Sub1CSV_D, 
  ABCDval_Sub2CSV_A, 
  ABCDval_Sub2CSV_B, 
  ABCDval_Sub2CSV_C, 
  ABCDval_Sub2CSV_D, 
  ABCDval_JetFlavor_A, 
  ABCDval_JetFlavor_B, 
  ABCDval_JetFlavor_C, 
  ABCDval_JetFlavor_D, 
  ABCDval_HiggsMass_,  
  ABCDval_JetFlavor_,  
  ABCDval_LeadCA8_Pt_,
 	
  ABCDval_0ak5_HT_A, 
  ABCDval_0ak5_HT_B, 
  ABCDval_0ak5_HT_C, 
  ABCDval_0ak5_HT_D,
  ABCDval_1ak5_HT_A, 
  ABCDval_1ak5_HT_B, 
  ABCDval_1ak5_HT_C, 
  ABCDval_1ak5_HT_D,
  ABCDval_2ak5_HT_A, 
  ABCDval_2ak5_HT_B, 
  ABCDval_2ak5_HT_C, 
  ABCDval_2ak5_HT_D,
 
  TH1_Size_
};

class TH1Info{
  public:
    TH1Info (
        bool 	output, 
        std::string	name, 
        std::string	title, 
        std::string 	xtitle, 
        std::string	ytitle, 
        std::string	unit, 
        int  	bin, 
        double  min, 
        double  max  
        ) : 
      _output(output), 
      _name  (name), 
      _title (title), 
      _xtitle(xtitle), 
      _ytitle(ytitle), 
      _unit  (unit), 
      _bin   (bin), 
      _min   (min), 
      _max   (max) { } 
    bool 	      _output;
    std::string	_name;
    std::string	_title;
    std::string	_xtitle;
    std::string	_ytitle;
    std::string	_unit;
    int  	      _bin;
    double      _min;
    double      _max;
};
#endif 

#ifdef TH1INFO_H
#ifndef TH1INFOCLASS_H
#define TH1INFOCLASS_H

template<typename TH1_type> 
class TH1InfoClass{
  public:
    TH1InfoClass() {};
    void CreateTH1();
    void CreateTH1(edm::Service<TFileService> f);
    void CreateTH1(TFile* f, std::string dir_name ); 
    void SetTitles(); 
    void Sumw2();
    TH1_type* GetTH1(std::string _name_);
    TH1Info GetVar(std::string _name_);
    TH1Info GetVar(int index);

  private:
    map<std::string, TH1_type*> mapTH1;
    map<std::string, int> indexTH1;
    static TH1Info Var[TH1_Size_];
};

template<typename TH1_type>
TH1Info TH1InfoClass<TH1_type>::Var[TH1_Size_] = { 
  TH1Info( 0,	"ABCDana_CutFlow",		"",	"Cut flow", 			"Yields", "", 		5, 0, 5 ) , 	
  TH1Info( 0,	"ABCDval_CutFlow",		"",	"Cut flow", 			"Yields", "", 		5, 0, 5 ) , 	
  TH1Info( 0,	"ABCDana_CutRegion",		"",	"Cut Region", 			"Yields", "", 		4, 0, 4 ) , 	
  TH1Info( 0,	"ABCDana_CutRegion_1b",		"",	"Cut Region, 1bjet", 		"Yields", "", 		4, 0, 4 ) , 	
  TH1Info( 0,	"ABCDana_CutRegion_2b",		"",	"Cut Region, over 1bjet", 	"Yields", "", 		4, 0, 4 ) , 	
  TH1Info( 0,	"ABCDval_CutRegion",		"",	"Cut Region", 			"Yields", "", 		4, 0, 4 ) , 	
  TH1Info( 0,	"ABCDval_CutRegion_0ak5",	"",	"Cut Region, 0 ak5 jet",	"Yields", "", 		4, 0, 4 ) , 	
  TH1Info( 0,	"ABCDval_CutRegion_1ak5",	"",	"Cut Region, 1 ak5 jet",	"Yields", "", 		4, 0, 4 ) , 	
  TH1Info( 0,	"ABCDval_CutRegion_2ak5",	"",	"Cut Region, 2 ak5 jet",	"Yields", "", 		4, 0, 4 ) , 	
  TH1Info( 0,	"ABCDana_Sumw2_A",		"",	"SumW2 for A", 			"Yields", "", 		1, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Sumw2_B",		"",	"SumW2 for B", 			"Yields", "", 		1, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Sumw2_C",		"",	"SumW2 for C", 			"Yields", "", 		1, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Sumw2_D",		"",	"SumW2 for D", 			"Yields", "", 		1, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Sumw2_1b_A",		"",	"SumW2 for A, 1bjet",		"Yields", "", 		1, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Sumw2_1b_B",		"",	"SumW2 for B, 1bjet", 		"Yields", "", 		1, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Sumw2_1b_C",		"",	"SumW2 for C, 1bjet", 		"Yields", "", 		1, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Sumw2_1b_D",		"",	"SumW2 for D, 1bjet", 		"Yields", "", 		1, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Sumw2_2b_A",		"",	"SumW2 for A, over 1bjet",	"Yields", "", 		1, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Sumw2_2b_B",		"",	"SumW2 for B, over 1bjet", 	"Yields", "", 		1, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Sumw2_2b_C",		"",	"SumW2 for C, over 1bjet", 	"Yields", "", 		1, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Sumw2_2b_D",		"",	"SumW2 for D, over 1bjet", 	"Yields", "", 		1, 0, 1 ) , 	
  TH1Info( 0,	"ABCDval_Sumw2_A",		"",	"SumW2 for A", 			"Yields", "", 		1, 0, 1 ) , 	
  TH1Info( 0,	"ABCDval_Sumw2_B",		"",	"SumW2 for B", 			"Yields", "", 		1, 0, 1 ) , 	
  TH1Info( 0,	"ABCDval_Sumw2_C",		"",	"SumW2 for C", 			"Yields", "", 		1, 0, 1 ) , 	
  TH1Info( 0,	"ABCDval_Sumw2_D",		"",	"SumW2 for D", 			"Yields", "", 		1, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_NAK5_BeforeHTAK5",		      "",	"Num of AK5 before HTAK5;N(AK5);","Events", "", 	20, -0.5, 19.5) , 	
  TH1Info( 0,	"ABCDana_NAK5",		      "",	"Num of AK5, before ABCD;N(AK5);","Events", "", 	20, -0.5, 19.5) , 	
  TH1Info( 0,	"ABCDana_NAK5_ABCD",		"",	"Num of AK5, ABCD region;N(AK5);","Events", "", 	20, -0.5, 19.5) , 	
  TH1Info( 0,	"ABCDana_NAK5_A",	    	"",	"Num of AK5, A region;N(AK%);", 	"Events", "", 	20, -0.5, 19.5) , 	
  TH1Info( 0,	"ABCDana_NAK5_B",		    "",	"Num of AK5, B region;N(AK%);", 	"Events", "", 	20, -0.5, 19.5) , 	
  TH1Info( 0,	"ABCDana_NAK5_C",		    "",	"Num of AK5, C region;N(AK%);", 	"Events", "", 	20, -0.5, 19.5) , 	
  TH1Info( 0,	"ABCDana_NAK5_D",  		  "",	"Num of CA8, D region;N(AK%);", 	"Events", "", 	20, -0.5, 19.5) , 	
  TH1Info( 0,	"ABCDana_Numbjet",		"",	"Num of bjet, before ABCD", 	"Events", "", 		4, 0, 4 ) , 	
  TH1Info( 0,	"ABCDana_Numbjet_ABCD",		"",	"Num of bjet, ABCD region", 	"Events", "", 		4, 0, 4 ) , 	
  TH1Info( 0,	"ABCDana_Numbjet_A",		"",	"Num of bjet, A region", 	"Events", "", 		4, 0, 4 ) , 	
  TH1Info( 0,	"ABCDana_Numbjet_B",		"",	"Num of bjet, B region", 	"Events", "", 		4, 0, 4 ) , 	
  TH1Info( 0,	"ABCDana_Numbjet_C",		"",	"Num of bjet, C region", 	"Events", "", 		4, 0, 4 ) , 	
  TH1Info( 0,	"ABCDana_Numbjet_D",		"",	"Num of bjet, D region;N(AK5);", 	"Events", "", 	4, 0, 4 ) , 	
  TH1Info( 0,	"ABCDana_NumCA8",		"",	"Num of CA8, before ABCD", 	"Events", "", 		4, 0, 4 ) , 	
  TH1Info( 0,	"ABCDana_NumCA8_ABCD",		"",	"Num of CA8, ABCD region", 	"Events", "", 		4, 0, 4 ) , 	
  TH1Info( 0,	"ABCDana_NumCA8_A",		"",	"Num of CA8, A region", 	"Events", "", 		4, 0, 4 ) , 	
  TH1Info( 0,	"ABCDana_NumCA8_B",		"",	"Num of CA8, B region", 	"Events", "", 		4, 0, 4 ) , 	
  TH1Info( 0,	"ABCDana_NumCA8_C",		"",	"Num of CA8, C region", 	"Events", "", 		4, 0, 4 ) , 	
  TH1Info( 0,	"ABCDana_NumCA8_D",		"",	"Num of CA8, D region", 	"Events", "", 		4, 0, 4 ) , 	
  TH1Info( 0,	"ABCDana_bpMass_A",		"",	"M(b'), A region",		"Yields", "GeV/c^{2 ) ", 	2000, 0, 2000 ) , 	
  TH1Info( 0,	"ABCDana_bpMass_B",		"",	"M(b'), B region",		"Yields", "GeV/c^{2 ) ", 	2000, 0, 2000 ) , 	
  TH1Info( 0,	"ABCDana_bpMass_C",		"",	"M(b'), C region",		"Yields", "GeV/c^{2 ) ", 	2000, 0, 2000 ) , 	
  TH1Info( 0,	"ABCDana_bpMass_D",		"",	"M(b'), D region",		"Yields", "GeV/c^{2 ) ", 	2000, 0, 2000 ) , 	
  TH1Info( 0,	"ABCDana_bpPt_A",		"",	"p_{T ) (b'), A region",	"Yields", "GeV/c", 	1500, 0, 1500 ) , 	
  TH1Info( 0,	"ABCDana_bpPt_B",		"",	"p_{T ) (b'), B region",	"Yields", "GeV/c", 	1500, 0, 1500 ) , 	
  TH1Info( 0,	"ABCDana_bpPt_C",		"",	"p_{T ) (b'), C region",	"Yields", "GeV/c", 	1500, 0, 1500 ) , 	
  TH1Info( 0,	"ABCDana_bpPt_D",		"",	"p_{T ) (b'), D region",	"Yields", "GeV/c", 	1500, 0, 1500 ) , 	
  TH1Info( 0,	"ABCDana_bpEta_A",		"",	"Eta, A region",		"Yields", "", 	800, -4., 4. ) , 	
  TH1Info( 0,	"ABCDana_bpEta_B",		"",	"Eta, B region",		"Yields", "", 	800, -4., 4. ) , 	
  TH1Info( 0,	"ABCDana_bpEta_C",		"",	"Eta, C region",		"Yields", "", 	800, -4., 4. ) , 	
  TH1Info( 0,	"ABCDana_bpEta_D",		"",	"Eta, D region",		"Yields", "", 	800, -4., 4. ) , 	
  TH1Info( 0,	"ABCDana_HT", 			"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDana_HT_A", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDana_HT_B", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDana_HT_C", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDana_HT_D", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDana_HiggsMass_A",		"",	"Pruned Mass(Higgs), A region",	"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDana_HiggsMass_B",	 	"",	"Pruned Mass(Higgs), B region",	"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDana_HiggsMass_C",		"",	"Pruned Mass(Higgs), C region",	"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDana_HiggsMass_D",	 	"",	"Pruned Mass(Higgs), D region",	"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDana_CA8Pt_A",	 	"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDana_CA8Pt_B",	 	"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDana_CA8Pt_C",	 	"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDana_CA8Pt_D",	 	"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDana_dRSubJets_A",		"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_dRSubJets_B",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_dRSubJets_C",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_dRSubJets_D",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Tau2ByTau1_A",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Tau2ByTau1_B",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Tau2ByTau1_C",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Tau2ByTau1_D",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Sub1Pt_A",	 	"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDana_Sub1Pt_B",	 	"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDana_Sub1Pt_C",	 	"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDana_Sub1Pt_D",	 	"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDana_Sub2Pt_A",	 	"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDana_Sub2Pt_B",	 	"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDana_Sub2Pt_C",	 	"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDana_Sub2Pt_D",	 	"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDana_Sub1Mass_A",	 	"",	"",				"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDana_Sub1Mass_B",	 	"",	"",				"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDana_Sub1Mass_C",	 	"",	"",				"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDana_Sub1Mass_D",	 	"",	"",				"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDana_Sub2Mass_A",	 	"",	"",				"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDana_Sub2Mass_B",	 	"",	"",				"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDana_Sub2Mass_C",	 	"",	"",				"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDana_Sub2Mass_D",	 	"",	"",				"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDana_Sub1CSV_A",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Sub1CSV_B",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Sub1CSV_C",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Sub1CSV_D",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Sub2CSV_A",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Sub2CSV_B",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Sub2CSV_C",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_Sub2CSV_D",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDana_JetFlavor_A",	 	"",	"Jet Flavor, A region",		"Yields", "", 		50, -25, 25 ) , 	
  TH1Info( 0,	"ABCDana_JetFlavor_B",	 	"",	"Jet Flavor, B region",		"Yields", "", 		50, -25, 25 ) , 	
  TH1Info( 0,	"ABCDana_JetFlavor_C",	 	"",	"Jet Flavor, C region",		"Yields", "", 		50, -25, 25 ) , 	
  TH1Info( 0,	"ABCDana_JetFlavor_D",	 	"",	"Jet Flavor, D region",		"Yields", "", 		50, -25, 25 ) , 	
  TH1Info( 0,	"ABCDana_HiggsMass", 		"",	"Pruned Mass(Higgs)", 		"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDana_JetFlavor",	 	"",	"Jet Flavor",			"Yields", "", 		50, -25, 25 ) , 	
  TH1Info( 0,	"ABCDana_CA8Pt", 		"",	"CA8 Pt", 			"Yields", "GeV/c", 	1500,  0, 1500 ) , 	
  TH1Info( 0,	"ABCDana_1b_bpMass_A",		"",	"M(b'), A region",		"Yields", "GeV/c^{2 ) ", 	2000, 0, 2000 ) , 	
  TH1Info( 0,	"ABCDana_1b_bpMass_B",	 	"",	"M(b'), B region",		"Yields", "GeV/c^{2 ) ", 	2000, 0, 2000 ) , 	
  TH1Info( 0,	"ABCDana_1b_bpMass_C",		"",	"M(b'), C region",		"Yields", "GeV/c^{2 ) ", 	2000, 0, 2000 ) , 	
  TH1Info( 0,	"ABCDana_1b_bpMass_D",	 	"",	"M(b'), D region",		"Yields", "GeV/c^{2 ) ", 	2000, 0, 2000 ) , 	
  TH1Info( 0,	"ABCDana_1b_bpPt_A",	 	"",	"p_{T ) (b'), A region",	"Yields", "GeV/c", 	1500, 0, 1500 ) , 	
  TH1Info( 0,	"ABCDana_1b_bpPt_B",	 	"",	"p_{T ) (b'), B region",	"Yields", "GeV/c", 	1500, 0, 1500 ) , 	
  TH1Info( 0,	"ABCDana_1b_bpPt_C",	 	"",	"p_{T ) (b'), C region",	"Yields", "GeV/c", 	1500, 0, 1500 ) , 	
  TH1Info( 0,	"ABCDana_1b_bpPt_D",	 	"",	"p_{T ) (b'), D region",	"Yields", "GeV/c", 	1500, 0, 1500 ) , 	
  TH1Info( 0,	"ABCDana_1b_bpEta_A",	 	"",	"Eta, A region",		"Yields", "", 	800, -4., 4. ) , 	
  TH1Info( 0,	"ABCDana_1b_bpEta_B",	 	"",	"Eta, B region",		"Yields", "", 	800, -4., 4. ) , 	
  TH1Info( 0,	"ABCDana_1b_bpEta_C",	 	"",	"Eta, C region",		"Yields", "", 	800, -4., 4. ) , 	
  TH1Info( 0,	"ABCDana_1b_bpEta_D",	 	"",	"Eta, D region",		"Yields", "", 	800, -4., 4. ) , 	
  TH1Info( 0,	"ABCDana_1b_HT_A", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDana_1b_HT_B", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDana_1b_HT_C", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDana_1b_HT_D", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDana_2b_bpMass_A",	 	"",	"M(b'), A region",		"Yields", "GeV/c^{2 ) ", 	2000, 0, 2000 ) , 	
  TH1Info( 0,	"ABCDana_2b_bpMass_B",	 	"",	"M(b'), B region",		"Yields", "GeV/c^{2 ) ", 	2000, 0, 2000 ) , 	
  TH1Info( 0,	"ABCDana_2b_bpMass_C",	 	"",	"M(b'), C region",		"Yields", "GeV/c^{2 ) ", 	2000, 0, 2000 ) , 	
  TH1Info( 0,	"ABCDana_2b_bpMass_D",		"",	"M(b'), D region",		"Yields", "GeV/c^{2 ) ", 	2000, 0, 2000 ) , 	
  TH1Info( 0,	"ABCDana_2b_bpPt_A",	 	"",	"p_{T ) (b'), A region",	"Yields", "GeV/c", 	1500, 0, 1500 ) , 	
  TH1Info( 0,	"ABCDana_2b_bpPt_B",	 	"",	"p_{T ) (b'), B region",	"Yields", "GeV/c", 	1500, 0, 1500 ) , 	
  TH1Info( 0,	"ABCDana_2b_bpPt_C",	 	"",	"p_{T ) (b'), C region",	"Yields", "GeV/c", 	1500, 0, 1500 ) , 	
  TH1Info( 0,	"ABCDana_2b_bpPt_D",	 	"",	"p_{T ) (b'), D region",	"Yields", "GeV/c", 	1500, 0, 1500 ) , 	
  TH1Info( 0,	"ABCDana_2b_bpEta_A",	 	"",	"Eta, A region",		"Yields", "", 	800, -4., 4. ) , 	
  TH1Info( 0,	"ABCDana_2b_bpEta_B",	 	"",	"Eta, B region",		"Yields", "", 	800, -4., 4. ) , 	
  TH1Info( 0,	"ABCDana_2b_bpEta_C",	 	"",	"Eta, C region",		"Yields", "", 	800, -4., 4. ) , 	
  TH1Info( 0,	"ABCDana_2b_bpEta_D",	 	"",	"Eta, D region",		"Yields", "", 	800, -4., 4. ) , 	
  TH1Info( 0,	"ABCDana_2b_HT_A", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDana_2b_HT_B", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDana_2b_HT_C", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDana_2b_HT_D", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDval_NumCA8_A",		"",	"Num of CA8, A region", 	"Events", "", 		3, 1, 4 ) , 	
  TH1Info( 0,	"ABCDval_NumCA8_B",		"",	"Num of CA8, B region", 	"Events", "", 		3, 1, 4 ) , 	
  TH1Info( 0,	"ABCDval_NumCA8_C",		"",	"Num of CA8, C region", 	"Events", "", 		3, 1, 4 ) , 	
  TH1Info( 0,	"ABCDval_NumCA8_D",		"",	"Num of CA8, D region", 	"Events", "", 		3, 1, 4 ) , 	
  TH1Info( 0,	"ABCDval_HT", 			"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDval_HT_A", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDval_HT_B", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDval_HT_C", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDval_HT_D", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDval_HiggsMass_A",	 	"",	"Pruned Mass(Higgs), A region",	"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDval_HiggsMass_B",		"",	"Pruned Mass(Higgs), B region",	"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDval_HiggsMass_C",	 	"",	"Pruned Mass(Higgs), C region",	"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDval_HiggsMass_D",		"",	"Pruned Mass(Higgs), D region",	"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDval_CA8Pt_A",		"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDval_CA8Pt_B",		"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDval_CA8Pt_C",		"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDval_CA8Pt_D",		"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDval_dRSubJets_A",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDval_dRSubJets_B",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDval_dRSubJets_C",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDval_dRSubJets_D",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDval_Tau2ByTau1_A",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDval_Tau2ByTau1_B",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDval_Tau2ByTau1_C",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDval_Tau2ByTau1_D",		"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDval_Sub1Pt_A",	 	"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDval_Sub1Pt_B",	 	"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDval_Sub1Pt_C",	 	"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDval_Sub1Pt_D",	 	"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDval_Sub2Pt_A",	 	"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDval_Sub2Pt_B",	 	"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDval_Sub2Pt_C",	 	"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDval_Sub2Pt_D",	 	"",	"",				"Yields", "GeV/c", 	700, 0, 700 ) , 	
  TH1Info( 0,	"ABCDval_Sub1Mass_A",	 	"",	"",				"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDval_Sub1Mass_B",	 	"",	"",				"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDval_Sub1Mass_C",	 	"",	"",				"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDval_Sub1Mass_D",	 	"",	"",				"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDval_Sub2Mass_A",	 	"",	"",				"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDval_Sub2Mass_B",	 	"",	"",				"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDval_Sub2Mass_C",	 	"",	"",				"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDval_Sub2Mass_D",	 	"",	"",				"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDval_Sub1CSV_A",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDval_Sub1CSV_B",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDval_Sub1CSV_C",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDval_Sub1CSV_D",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDval_Sub2CSV_A",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDval_Sub2CSV_B",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDval_Sub2CSV_C",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDval_Sub2CSV_D",	 	"",	"",				"Yields", "", 	100, 0, 1 ) , 	
  TH1Info( 0,	"ABCDval_JetFlavor_A",	 	"",	"Jet Flavor, A region",		"Yields", "", 		50, -25, 25 ) , 	
  TH1Info( 0,	"ABCDval_JetFlavor_B",	 	"",	"Jet Flavor, B region",		"Yields", "", 		50, -25, 25 ) , 	
  TH1Info( 0,	"ABCDval_JetFlavor_C",	 	"",	"Jet Flavor, C region",		"Yields", "", 		50, -25, 25 ) , 	
  TH1Info( 0,	"ABCDval_JetFlavor_D",	 	"",	"Jet Flavor, D region",		"Yields", "", 		50, -25, 25 ) , 	
  TH1Info( 0,	"ABCDval_HiggsMass", 		"",	"Pruned Mass(Higgs)", 		"Yields", "GeV/c^{2 ) ", 	200, 0, 200 ) , 	
  TH1Info( 0,	"ABCDval_JetFlavor",	 	"",	"Jet Flavor",			"Yields", "", 		50, -25, 25 ) , 	
  TH1Info( 0,	"ABCDval_CA8Pt", 		"",	"CA8 Pt", 			"Yields", "GeV/c", 	1500,  0, 1500 ) ,  
	
  TH1Info( 0,	"ABCDval_0ak5_HT_A", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDval_0ak5_HT_B", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDval_0ak5_HT_C", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDval_0ak5_HT_D", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDval_1ak5_HT_A", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDval_1ak5_HT_B", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDval_1ak5_HT_C", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDval_1ak5_HT_D", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDval_2ak5_HT_A", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDval_2ak5_HT_B", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDval_2ak5_HT_C", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
  TH1Info( 0,	"ABCDval_2ak5_HT_D", 		"",	"HT(AK5)", 			"Yields", "GeV/c", 	2500, 0, 2500 ) , 	
};

//// Create Histogram
template<typename TH1_type> 
void TH1InfoClass<TH1_type>::CreateTH1(){
  for(int i=0; i<TH1_Size_; i++){
    mapTH1[Var[i]._name] = new TH1(Var[i]._name.c_str(),"",Var[i]._bin, Var[i]._min, Var[i]._max);
  }
}

template<typename TH1_type> 
void TH1InfoClass<TH1_type>::CreateTH1(edm::Service<TFileService> f){
  for(int i=0; i<TH1_Size_; i++){
    mapTH1[Var[i]._name] =f->make<TH1_type>(Var[i]._name.c_str(),"",Var[i]._bin, Var[i]._min, Var[i]._max);
  }
}

template<typename TH1_type> 
void TH1InfoClass<TH1_type>::CreateTH1( TFile* f, std::string dir_name="" ){
  for(int i=0; i<TH1_Size_; i++){ 
    mapTH1[Var[i]._name] =(TH1_type*)f->Get( (dir_name+Var[i]._name).c_str() );
  }

}

// Set some option for Histogram
template<typename TH1_type> 
void TH1InfoClass<TH1_type>::SetTitles(){
  for(int i=0; i<TH1_Size_; i++){ 
    //mapTH1[Var[i]._name]->SetTile(Var[i]._title.c_str());
    mapTH1[Var[i]._name]->SetX_title( (Var[i]._xtitle+" ["+Var[i]._unit+"]").c_str());
    mapTH1[Var[i]._name]->SetY_title( Var[i]._ytitle.c_str() );
  }
}

template<typename TH1_type> 
void TH1InfoClass<TH1_type>::Sumw2(){
  for(int i=0; i<TH1_Size_; i++){ 
    mapTH1.find(Var[i]._name)->second->Sumw2();
  }
}

// Get Histogram
template<typename TH1_type> 
TH1_type* TH1InfoClass<TH1_type>::GetTH1(std::string _name_){
  return mapTH1.find(_name_)->second;
}

// Get Variables
template<typename TH1_type> 
TH1Info TH1InfoClass<TH1_type>::GetVar(std::string _name_){
  return Var[indexTH1.find(_name_)->second];
}
template<typename TH1_type> 
TH1Info TH1InfoClass<TH1_type>::GetVar(int index){
  return Var[index];
}

#endif
#endif
