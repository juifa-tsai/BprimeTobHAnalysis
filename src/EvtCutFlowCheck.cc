// -*- C++ -*-
//
// Package:    EvtCutFlowCheck
// Class:      EvtCutFlowCheck
// 
/**\class EvtCutFlowCheck EvtCutFlowCheck.cc Bprime_kit/EvtCutFlowCheck/src/EvtCutFlowCheck.cc

Description: 
Analyzer class for Bprime -> b Higgs studies 
- National Taiwan University - 

Implementation:
[Notes on implementation]
 */
//
// Original Author:  Jui-Fa Tsai
//         Created:  Tue Jul 16 19:48:47 CEST 2013
// Second Author:    Devdatta Majumder 
// $Id$
//
//

// system include files
#include <memory>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <assert.h>
#include <vector>
#include <map>
#include <string>

// Root headers 
#include <TLorentzVector.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH1I.h>
#include <TEfficiency.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "BpbH/BprimeTobH/interface/format.h"
#include "BpbH/BprimeTobH/interface/TriggerBooking.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h" 

//Evt selection
#include "BpbH/BprimeTobHAnalysis/interface/EventSelector.h"
#include "BpbH/BprimeTobHAnalysis/interface/HiggsBRscaleFactors.h" 

//For mini tree
#include "BpbH/BprimeTobHAnalysis/interface/reRegistJet.hh"

//Plots information
//#include "BpbH/BprimeTobHAnalysis/interface/TH2InfoClass.h"
//#include "BpbH/BprimeTobHAnalysis/interface/TH1InfoClass.h"
///// Jet Correction
//#include "BpbH/BprimeTobHAnalysis/src/JMEUncertUtil.cc"
//#include "BpbH/BprimeTobHAnalysis/src/BTagSFUtil.cc"
//#include "BpbH/BprimeTobHAnalysis/src/ApplyBTagSF.cc"
//#include "BpbH/BprimeTobHAnalysis/src/ApplyHiggsTagSF.cc"

//#include "BpbH/BprimeTobHAnalysis/interface/JMEUncertUtil.h"
//#include "BpbH/BprimeTobHAnalysis/interface/ApplyBTagSF.h"
//#include "BpbH/BprimeTobHAnalysis/interface/ApplyHiggsTagSF.h"

//
// class declaration
//

class EvtCutFlowCheck : public edm::EDAnalyzer{
	public:
		explicit EvtCutFlowCheck(const edm::ParameterSet&);
		~EvtCutFlowCheck();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob();
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob();

		Jet CloseJetIndex( JetCollection jets_, const Jet myJet_ );
		void isolateCollection( JetCollection control_, JetCollection input_, JetCollection& output_, double dR_=1.2 );
	
		edm::LumiReWeighting LumiWeights_; 

		// ----------member data ---------------------------

		//// Configurables 

		int                             maxEvents_; 
		const int                       reportEvery_; 
		const std::string               inputTTree_;
		const std::vector<std::string>  inputFiles_;
		const edm::ParameterSet         hltPaths_; 
		const bool     DoHLTSelect_;
		const bool     DoGoodVtxSelect_;
		const bool     DoPUReweighting_;
		const std::string               file_PUDistMC_;
		const std::string               file_PUDistData_;
		const std::string               hist_PUDistMC_;
		const std::string               hist_PUDistData_;
		const double     HTAK5Min_;
		const double     HTAK5Max_;

		const double     HJetPtMin_;
		const double     HJetPtMax_;
		const double     HJetAbsEtaMin_;
		const double     HJetAbsEtaMax_;
		const double     HJetMassMin_;	
		const double     HJetMassMax_;	
		const double     HJetPrunedMassMin_;	
		const double     HJetPrunedMassMax_;	
		const double     HJetTau2ByTau1Min_;
		const double     HJetTau2ByTau1Max_;
		const double     Subjet1CSVDiscMin_;
		const double     Subjet1CSVDiscMax_;
		const double     Subjet2CSVDiscMin_;
		const double     Subjet2CSVDiscMax_;
		const double     dRSubjetsMin_;
		const double     dRSubjetsMax_;

		const double     HJetSBMassMin_;	
		const double     HJetSBMassMax_;

		const double     JetPtMin_;
		const double     JetPtMax_;
		const double     JetAbsEtaMin_;
		const double     JetAbsEtaMax_;
		const double     bJetPtMin_;
		const double     bJetPtMax_;
		const double     bJetCSVDiscMin_;
		const double     bJetCSVDiscMax_;

		const double     bVetoJetPtMin_;
		const double     bVetoJetPtMax_;
		const double     bVetoJetAbsEtaMin_;
		const double     bVetoJetAbsEtaMax_;
		const double     bVetoJetCSVDiscMin_;
		const double     bVetoJetCSVDiscMax_;

		const int     numbJetMin_;
		const int     numHiggsJetMin_;

		//const double jesShift_;
		//const double jerShift_; 
		//const double SFbShift_;
		//const double SFlShift_;

		bool BuildMinTree_;

		TChain*            chain_;

		GenInfoBranches GenInfo;
		EvtInfoBranches    EvtInfo;
		VertexInfoBranches VtxInfo;
		JetInfoBranches CA8JetInfo;
		JetInfoBranches AK5JetInfo;
		JetInfoBranches SubJetInfo;

		edm::Service<TFileService> fs; 
		
		bool isData_; 
		bool evtwt_;
		double puweight_;
		double higgsTagCorr_;

		// New branch
		//int GenEvt_;
		//double xsec_;
		bool McFlagana;
		double PUana;
		double evtWtana;
		double HTak5, HThiggsbjet;

		//variables
		int	evtPass_ana;
		int 	evtPass_val;

		//TH1D:
		TH1D* CutFlow_;
		TH1D* CutFlow_UnWt_;
		TH1D* HT_;
		TH1D* HT_UnWt_;

};

//
// constructors and destructor
//
EvtCutFlowCheck::EvtCutFlowCheck(const edm::ParameterSet& iConfig) : 
	maxEvents_(iConfig.getParameter<int>("MaxEvents")), 
	reportEvery_(iConfig.getParameter<int>("ReportEvery")),
	inputTTree_(iConfig.getParameter<std::string>("InputTTree")),
	inputFiles_(iConfig.getParameter<std::vector<std::string> >("InputFiles")),

	hltPaths_(iConfig.getParameter<edm::ParameterSet>("HLTPaths")),

	DoHLTSelect_(iConfig.getParameter<bool>("DoHLTSelect")),
	DoGoodVtxSelect_(iConfig.getParameter<bool>("DoGoodVtxSelect")),
	DoPUReweighting_(iConfig.getParameter<bool>("DoPUReweighting")),
	file_PUDistMC_(iConfig.getParameter<std::string>("File_PUDistMC")),
	file_PUDistData_(iConfig.getParameter<std::string>("File_PUDistData")),
	hist_PUDistMC_(iConfig.getParameter<std::string>("Hist_PUDistMC")),
	hist_PUDistData_(iConfig.getParameter<std::string>("Hist_PUDistData")),

	HTAK5Min_(iConfig.getParameter<double>("HTAK5Min")),
	HTAK5Max_(iConfig.getParameter<double>("HTAK5Max")),

	HJetPtMin_(iConfig.getParameter<double>("HJetPtMin")),
	HJetPtMax_(iConfig.getParameter<double>("HJetPtMax")),
	HJetAbsEtaMin_(iConfig.getParameter<double>("HJetAbsEtaMin")),
	HJetAbsEtaMax_(iConfig.getParameter<double>("HJetAbsEtaMax")),
	HJetMassMin_(iConfig.getParameter<double>("HJetMassMin")),	
	HJetMassMax_(iConfig.getParameter<double>("HJetMassMax")),	
	HJetPrunedMassMin_(iConfig.getParameter<double>("HJetPrunedMassMin")),	
	HJetPrunedMassMax_(iConfig.getParameter<double>("HJetPrunedMassMax")),	
	HJetTau2ByTau1Min_(iConfig.getParameter<double>("HJetTau2ByTau1Min")),
	HJetTau2ByTau1Max_(iConfig.getParameter<double>("HJetTau2ByTau1Max")),
	Subjet1CSVDiscMin_(iConfig.getParameter<double>("Subjet1CSVDiscMin")),
	Subjet1CSVDiscMax_(iConfig.getParameter<double>("Subjet1CSVDiscMax")),
	Subjet2CSVDiscMin_(iConfig.getParameter<double>("Subjet2CSVDiscMin")),
	Subjet2CSVDiscMax_(iConfig.getParameter<double>("Subjet2CSVDiscMax")),
	dRSubjetsMin_(iConfig.getParameter<double>("dRSubjetsMin")),
	dRSubjetsMax_(iConfig.getParameter<double>("dRSubjetsMax")),

	HJetSBMassMin_(iConfig.getParameter<double>("HJetSBMassMin")),	
	HJetSBMassMax_(iConfig.getParameter<double>("HJetSBMassMax")),

	JetPtMin_(iConfig.getParameter<double>("JetPtMin")),
	JetPtMax_(iConfig.getParameter<double>("JetPtMax")),
	JetAbsEtaMin_(iConfig.getParameter<double>("JetAbsEtaMin")),
	JetAbsEtaMax_(iConfig.getParameter<double>("JetAbsEtaMax")),
	bJetPtMin_(iConfig.getParameter<double>("bJetPtMin")),
	bJetPtMax_(iConfig.getParameter<double>("bJetPtMax")),
	bJetCSVDiscMin_(iConfig.getParameter<double>("bJetCSVDiscMin")),
	bJetCSVDiscMax_(iConfig.getParameter<double>("bJetCSVDiscMax")),

	bVetoJetPtMin_(iConfig.getParameter<double>("bVetoJetPtMin")),
	bVetoJetPtMax_(iConfig.getParameter<double>("bVetoJetPtMax")),
	bVetoJetAbsEtaMin_(iConfig.getParameter<double>("bVetoJetAbsEtaMin")),
	bVetoJetAbsEtaMax_(iConfig.getParameter<double>("bVetoJetAbsEtaMax")),
	bVetoJetCSVDiscMin_(iConfig.getParameter<double>("bVetoJetCSVDiscMin")),
	bVetoJetCSVDiscMax_(iConfig.getParameter<double>("bVetoJetCSVDiscMax")),

	numbJetMin_(iConfig.getParameter<int>("numbJetMin")),
	numHiggsJetMin_(iConfig.getParameter<int>("numHiggsJetMin")),

	BuildMinTree_(iConfig.getParameter<bool>("BuildMinTree")),
/*	higgsJetSelParame_(iConfig.getParameter<edm::ParameterSet>("HiggsJetSelParams")),
	jetSelParams_(iConfig.getParameter<edm::ParameterSet>("JetSelParams")), 
	bjetSelParams_(iConfig.getParameter<edm::ParameterSet>("BJetSelParams")), 
	evtSelParams_(iConfig.getParameter<edm::ParameterSet>("EvtSelParams")),
	jmeParams_(iConfig.getParameter<edm::ParameterSet>("JMEParams")),
	jesShift_(iConfig.getParameter<double>("JESShift")),
	jerShift_(iConfig.getParameter<double>("JERShift")),
	SFbShift_(iConfig.getParameter<double>("SFbShift")),
	SFlShift_(iConfig.getParameter<double>("SFlShift")),
*/
	isData_(0),
	evtwt_(1),
	puweight_(1),
	higgsTagCorr_(1),
	evtPass_ana(0),  
	evtPass_val(0) 
{ 
	if( DoPUReweighting_) LumiWeights_ = edm::LumiReWeighting(file_PUDistMC_, file_PUDistData_, hist_PUDistMC_, hist_PUDistData_);
}


EvtCutFlowCheck::~EvtCutFlowCheck(){ 
	delete chain_;
}

// ------------ Other function -------------
Jet EvtCutFlowCheck::CloseJetIndex( JetCollection jets_, const Jet myJet_ ){
	Jet min_Jet;
	for( JetCollection::const_iterator ijet = jets_.begin(); ijet != jets_.end(); ++ijet ){
		float DR_min = 1000.;
		if( myJet_.DeltaR(*ijet) < DR_min){
			DR_min = myJet_.DeltaR(*ijet);
			min_Jet = *ijet;
		}
	}
	if( jets_.size() != 0 ){
		return min_Jet;
	}else{
		std::cout<<"ERROR: JetCollection size is 0 !!!"<<std::endl;
		Jet emptyJet; emptyJet.Set_Pt(-1);
		return emptyJet;
	}
}
void EvtCutFlowCheck::isolateCollection( JetCollection control_, JetCollection input_, JetCollection& output_, double dR_ ){
	for( JetCollection::const_iterator _i = input_.begin(); _i != input_.end(); ++_i ){
		bool isControl(false); 
		for( JetCollection::const_iterator _c = control_.begin(); _c != control_.end(); ++_c ){
			if( _c->DeltaR(*_i) < dR_ ){
				isControl = true; 
				break; 
			}else {
				isControl = false; 
			} 
		}
		if( !isControl ){
			output_.push_back(*_i);
		} 
	}
}
// ------------ method called once each job just before starting event loop  ------------
void EvtCutFlowCheck::beginJob(){ 
	chain_  = new TChain(inputTTree_.c_str());

	for(unsigned i=0; i<inputFiles_.size(); ++i){
		chain_->Add(inputFiles_.at(i).c_str());
		TFile *f = TFile::Open(inputFiles_.at(i).c_str(),"READ");
		f->Close();
	}

	EvtInfo.Register(chain_);
	VtxInfo.Register(chain_);
	GenInfo.Register(chain_);
	AK5JetInfo.Register(chain_,"JetInfo");
	CA8JetInfo.Register(chain_,"FatJetInfo");
	SubJetInfo.Register(chain_,"SubJetInfo");

	if(  maxEvents_<0 || maxEvents_>chain_->GetEntries()) maxEvents_ = chain_->GetEntries();
	CutFlow_=fs->make<TH1D>("CutFlow", "", 8, 0, 8);
	CutFlow_UnWt_=fs->make<TH1D>("CutFlow_UnWt", "", 8, 0, 8);;
	HT_=fs->make<TH1D>("HT", "", 3000, 0, 3000);
	HT_UnWt_=fs->make<TH1D>("HT_UnWt", "", 3000, 0, 3000);

	CutFlow_->GetXaxis()->SetBinLabel(1,"All_Evt");	
	CutFlow_->GetXaxis()->SetBinLabel(2,"Trigger_Sel");	
	CutFlow_->GetXaxis()->SetBinLabel(3,"Vertex_Sel");	
	CutFlow_->GetXaxis()->SetBinLabel(4,"HiggsJets_Sel");	
	CutFlow_->GetXaxis()->SetBinLabel(5,"bJets_Sel");	
	CutFlow_->GetXaxis()->SetBinLabel(6,"HT_Sel");	
	CutFlow_UnWt_->GetXaxis()->SetBinLabel(1,"All_Evt");	
	CutFlow_UnWt_->GetXaxis()->SetBinLabel(2,"Trigger_Sel");	
	CutFlow_UnWt_->GetXaxis()->SetBinLabel(3,"Vertex_Sel");	
	CutFlow_UnWt_->GetXaxis()->SetBinLabel(4,"HiggsJets_Sel");	
	CutFlow_UnWt_->GetXaxis()->SetBinLabel(5,"bJets_Sel");	
	CutFlow_UnWt_->GetXaxis()->SetBinLabel(6,"HT_Sel");	
	
	return;  
}

// ------------ method called for each event  ------------
void EvtCutFlowCheck::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){ 
	using namespace edm;
	using namespace std;

	if(  chain_ == 0) return;

	edm::LogInfo("StartingAnalysisLoop") << "Starting analysis loop\n";
	
	//// Roop events ==================================================================================================	
	for(int entry=0; entry<maxEvents_; entry++){
		if( (entry%reportEvery_) == 0) edm::LogInfo("Event") << entry << " of " << maxEvents_;
		chain_->GetEntry(entry);

		bool passHLT(false); 
		int nGoodVtxs(0);
		double HT_AK5(0);	
		double HT_HiggsBJets(0); 
		double weight_(1);

		isData_  = EvtInfo.McFlag ? 0 : 1; 
		if( !isData_ ) evtwt_ = EvtInfo.Weight; 
		if( DoPUReweighting_ && !isData_ ) puweight_ = LumiWeights_.weight(EvtInfo.TrueIT[0]); 

		//// Higgs BR reweighting, for b'b'>bHbH sample
		double br_=1.;
		if ( !isData_ ) {
			//int nbprimeBH(0), nbprimeBZ(0), nbprimeTW(0); 
			for (int igen=0; igen < GenInfo.Size; ++igen) {
				if ( GenInfo.Status[igen] == 3 && TMath::Abs(GenInfo.PdgID[igen]) == 25 && GenInfo.nDa[igen] >= 2 ) { //// H found
					int higgsDau0 = abs(GenInfo.Da0PdgID[igen]);
					int higgsDau1 = abs(GenInfo.Da1PdgID[igen]);
					int higgsDau = higgsDau0+higgsDau1; 
					if(higgsDau==10) br_ *= HiggsBRscaleFactors::higgsBBSf;
					if(higgsDau==30) br_ *= HiggsBRscaleFactors::higgsTauTauSf;
					if(higgsDau==26) br_ *= HiggsBRscaleFactors::higgsMuMuSf;
					if(higgsDau==8)  br_ *= HiggsBRscaleFactors::higgsCCSf;
					if(higgsDau==6)  br_ *= HiggsBRscaleFactors::higgsSSSf;
					if(higgsDau==16) br_ *= HiggsBRscaleFactors::higgsTTSf;
					if(higgsDau==42) br_ *= HiggsBRscaleFactors::higgsGGSf;
					if(higgsDau==44) br_ *= HiggsBRscaleFactors::higgsGammaGammaSf;
					if(higgsDau==45) br_ *= HiggsBRscaleFactors::higgsZGammaSf;
					if(higgsDau==48) br_ *= HiggsBRscaleFactors::higgsWWSf;
					if(higgsDau==46) br_ *= HiggsBRscaleFactors::higgsZZSf; 
				}
			}
		}
		weight_ = puweight_*br_*evtwt_;

		CutFlow_->Fill(0);	
		CutFlow_UnWt_->Fill(0);
	
		//// Trigger selection =====================================================================================
		TriggerSelector trigSel(hltPaths_); 
		passHLT = trigSel.getTrigDecision(EvtInfo); 
		if( !passHLT ) continue; 
		CutFlow_->Fill(double(1),weight_);	
		CutFlow_UnWt_->Fill(double(1), 1.);

		//// Vertex selection =====================================================================================
		VertexSelector vtxSel(VtxInfo); 
		nGoodVtxs = vtxSel.NGoodVtxs(); 
		if( nGoodVtxs < 1){ edm::LogInfo("NoGoodPrimaryVertex") << " No good primary vertex "; continue; }
		CutFlow_->Fill(double(2),weight_);	
		CutFlow_UnWt_->Fill(double(2), 1.);

		////  Higgs jets selection ================================================================================ 
		JetCollection HiggsJets, HiggsLikeJets, AllHiggsJets, AllHiggsJetsDRCut;
		JetCollection AntiHiggsJets, AntiHiggsLikeJets, AllAntiHiggsJets, AllAntiHiggsJetsDRCut;
		for( int i=0; i < CA8JetInfo.Size; ++i ){
			if( CA8JetInfo.NHF[i]>=0.99 || CA8JetInfo.NEF[i]>=0.99 || CA8JetInfo.NConstituents[i]<=1 || fabs(CA8JetInfo.Eta[i]) >=2.4 || CA8JetInfo.CHF[i]<=0 || CA8JetInfo.CEF[i]>=0.99 ||  CA8JetInfo.NCH[i]<=0 ) continue; //// apply loose jet ID
			if( fabs(CA8JetInfo.Eta[i]) > HJetAbsEtaMax_ ) continue; 
			if( CA8JetInfo.Pt[i] < HJetPtMin_ || CA8JetInfo.Pt[i] > HJetPtMax_ ) continue; //// apply jet pT cut

			//// Get subjets of fat jets
			int iSub1 = CA8JetInfo.Jet_SubJet1Idx[i];
			int iSub2 = CA8JetInfo.Jet_SubJet2Idx[i];
			if( iSub1 > 25 || iSub2 > 25) continue; //// skip fat jets only one subjet 
			if( SubJetInfo.Pt[iSub1]==0. || SubJetInfo.Pt[iSub2]==0. ) continue; //// skip fat jets for which one of the subjets has pT=0
			TLorentzVector Subjet1, Subjet2;
			Subjet1.SetPtEtaPhiM(SubJetInfo.Pt[iSub1], SubJetInfo.Eta[iSub1], SubJetInfo.Phi[iSub1], SubJetInfo.Mass[iSub1]);
			Subjet2.SetPtEtaPhiM(SubJetInfo.Pt[iSub2], SubJetInfo.Eta[iSub2], SubJetInfo.Phi[iSub2], SubJetInfo.Mass[iSub2]);
			double subjet_dy = Subjet1.Rapidity() - Subjet2.Rapidity();
			double subjet_dphi = Subjet1.DeltaPhi(Subjet2);
			double subjet_dyphi = sqrt( subjet_dy*subjet_dy + subjet_dphi*subjet_dphi );

			if( subjet_dyphi <(CA8JetInfo.Mass[i]/CA8JetInfo.Pt[i]) ) continue; //// skip infrared unsafe configurations
			if( CA8JetInfo.tau2[i]/CA8JetInfo.tau1[i] > HJetTau2ByTau1Max_ ) continue;
			if( SubJetInfo.CombinedSVBJetTags[iSub1]<0.244 || SubJetInfo.CombinedSVBJetTags[iSub2]<0.244 )  continue;

			if( SubJetInfo.CombinedSVBJetTags[iSub1] < Subjet1CSVDiscMin_ && SubJetInfo.CombinedSVBJetTags[iSub2]< Subjet2CSVDiscMin_){ //Anti-Higgs
				Jet CA8Jet(CA8JetInfo, i);
				AllAntiHiggsJets.push_back(CA8Jet);
				if( CA8JetInfo.MassPruned[i] > HJetPrunedMassMin_ && CA8JetInfo.MassPruned[i] < HJetPrunedMassMax_ ){  //// apply pruned jet mass cut 
					if( subjet_dyphi >= dRSubjetsMin_ && subjet_dyphi <= dRSubjetsMax_  ){
						AntiHiggsJets.push_back(CA8Jet);
						AllAntiHiggsJetsDRCut.push_back(CA8Jet);
					}
				}else{
					AntiHiggsLikeJets.push_back(CA8Jet);
					AllAntiHiggsJetsDRCut.push_back(CA8Jet);
				}
			}else if( SubJetInfo.CombinedSVBJetTags[iSub1] >= Subjet1CSVDiscMin_ && SubJetInfo.CombinedSVBJetTags[iSub2] >= Subjet2CSVDiscMin_ ){ //Higgs
				Jet CA8Jet(CA8JetInfo, i);
				AllHiggsJets.push_back(CA8Jet);
				if( CA8JetInfo.MassPruned[i] > HJetPrunedMassMin_ && CA8JetInfo.MassPruned[i] < HJetPrunedMassMax_ ){  //// apply pruned jet mass cut 
					if( subjet_dyphi >= dRSubjetsMin_ && subjet_dyphi <= dRSubjetsMax_  ){
						HiggsJets.push_back(CA8Jet);
						AllHiggsJetsDRCut.push_back(CA8Jet);
						HT_HiggsBJets = HT_HiggsBJets + CA8JetInfo.Pt[i]; 
					}
				}else{
					HiggsLikeJets.push_back(CA8Jet);
					AllHiggsJetsDRCut.push_back(CA8Jet);
				}
			}
		}//CA8

		//// bJet and bVeto selection ================================================================================
		JetCollection bJets;
		JetCollection bJets_Veto;
		for( int i=0; i<AK5JetInfo.Size; ++i){ 
			if( AK5JetInfo.NHF[i]>=0.90 || AK5JetInfo.NEF[i]>=0.90 || AK5JetInfo.NConstituents[i]<=1 || fabs(AK5JetInfo.Eta[i]) >=2.4 || AK5JetInfo.CHF[i]<=0 || AK5JetInfo.CEF[i]>=0.99 ||  AK5JetInfo.NCH[i]<=0 ) continue; //// apply loose jet ID
			if( fabs(AK5JetInfo.Eta[i]) > JetAbsEtaMax_ ) continue; 
			if( AK5JetInfo.Pt[i] < JetPtMin_ || AK5JetInfo.Pt[i] > JetPtMax_ ) continue; 

			HT_AK5 = HT_AK5 + AK5JetInfo.Pt[i];

			if( AK5JetInfo.CombinedSVBJetTags[i] <= bJetCSVDiscMin_ ) continue;
			if( AK5JetInfo.Pt[i] <= bJetPtMin_ ) continue;
			Jet bJet(AK5JetInfo, i);
			bJets.push_back(bJet);
		}//bJet end
		/*for( int i=0; i<AK5JetInfo.Size; ++i){ 
			if( AK5JetInfo.NHF[i]>=0.90 || AK5JetInfo.NEF[i]>=0.90 || AK5JetInfo.NConstituents[i]<=1 || fabs(AK5JetInfo.Eta[i]) >=2.4 || AK5JetInfo.CHF[i]<=0 || AK5JetInfo.CEF[i]>=0.99 ||  AK5JetInfo.NCH[i]<=0 ) continue; //// apply loose jet ID
			if( fabs(AK5JetInfo.Eta[i]) > bVetoJetAbsEtaMax_ ) continue; 
			if( AK5JetInfo.Pt[i] < bVetoJetPtMin_ ) continue; 
			if( AK5JetInfo.CombinedSVBJetTags[i] <= bVetoJetCSVDiscMin_ ) continue;
			Jet bJet(AK5JetInfo, i);
			bJets_Veto.push_back(bJet);
		}//bVeto's bJet selection end*/

		// dR(b, H)>1.2 to Isolate bJet and HJet //================================================================================================================
		JetCollection bJetsNotHiggs; //Only for HT(Higgs, bjet)
		JetCollection bJetsNotAllHiggs;
		JetCollection bJetsNotAllHiggsAllAntiHiggs;
		isolateCollection( HiggsJets, 		bJets, 			bJetsNotHiggs, 			1.2); //Only for HT(Higgs, bjet)
		isolateCollection( AllHiggsJets, 	bJets, 			bJetsNotAllHiggs, 		1.2);
		isolateCollection( AllAntiHiggsJets, 	bJetsNotAllHiggs, 	bJetsNotAllHiggsAllAntiHiggs, 	1.2);

		///// Fill evt and ABCD plots //==================================================================================================================================================	
		if( HiggsJets.size() < unsigned(numHiggsJetMin_) ) continue;		
		CutFlow_->Fill(double(3),weight_);	
		CutFlow_UnWt_->Fill(double(3),1.);
	
		if( bJetsNotAllHiggsAllAntiHiggs.size() < unsigned(numbJetMin_) ) continue;	
		CutFlow_->Fill(double(4),weight_);	
		CutFlow_UnWt_->Fill(double(4),1.);

		if( HT_AK5 < HTAK5Min_ ) continue;
		CutFlow_->Fill(double(5),weight_);	
		CutFlow_UnWt_->Fill(double(5),1.);
		HT_->Fill(HT_AK5,weight_);
		HT_UnWt_->Fill(HT_AK5);

		evtPass_ana++;
	} //// entry loop 
}

// ------------ method called once each job just after ending the event loop  ------------
void EvtCutFlowCheck::endJob(){
	std::cout<<std::endl;	
	std::cout<<"evtPass_ana "<<evtPass_ana<<"/"<<maxEvents_<<std::endl;	
	std::cout<<"evtPass_val "<<evtPass_val<<"/"<<maxEvents_<<std::endl;	
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void EvtCutFlowCheck::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EvtCutFlowCheck);
