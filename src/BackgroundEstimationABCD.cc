// -*- C++ -*-
//
// Package:    BackgroundEstimationABCD
// Class:      BackgroundEstimationABCD
// 
/**\class BackgroundEstimationABCD BackgroundEstimationABCD.cc Bprime_kit/BackgroundEstimationABCD/src/BackgroundEstimationABCD.cc

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
#include <TEfficiency.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "BpbH/BprimeTobH/interface/format.h"
#include "BpbH/BprimeTobH/interface/TriggerBooking.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h" 

#include "BpbH/BprimeTobHAnalysis/interface/EventSelector.h"
#include "BpbH/BprimeTobHAnalysis/interface/JMEUncertUtil.h"
#include "BpbH/BprimeTobHAnalysis/interface/ApplyBTagSF.h"
#include "BpbH/BprimeTobHAnalysis/interface/ApplyHiggsTagSF.h"
#include "BpbH/BprimeTobHAnalysis/interface/HiggsBRscaleFactors.h" 

#include "BpbH/BprimeTobHAnalysis/interface/TH2InfoClass.h"
#include "BpbH/BprimeTobHAnalysis/interface/TH1InfoClass.h"
#include "BpbH/BprimeTobHAnalysis/interface/reRegistJet.hh"

//
// class declaration
//

class BackgroundEstimationABCD : public edm::EDAnalyzer{
	public:
		explicit BackgroundEstimationABCD(const edm::ParameterSet&);
		~BackgroundEstimationABCD();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob();
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob();

		Jet NearestJet( JetCollection jets, const Jet myJet );
		void isolateCollection( JetCollection control, JetCollection input, JetCollection& output, double dR=1.2 );
		void recoBprime( JetCollection bjets, const Jet H, vector<TLorentzVector>& p4_output );
		template <typename Type>
			void setABCDcutRegion(Type* h1);

		edm::LumiReWeighting LumiWeights_; 

		// ----------member data ---------------------------

		//// Configurables 

		int                             maxEvents_; 
		const int                       reportEvery_; 
		const std::string               inputTTree_;
		const std::vector<std::string>  inputFiles_;
		const edm::ParameterSet         hltPaths_; 
		const bool                      doPUReweighting_;
		const std::string               file_PUDistMC_;
		const std::string               file_PUDistData_;
		const std::string               hist_PUDistMC_;
		const std::string               hist_PUDistData_;

		const double jetPtMin_;
		const double bJetPtMin_;
		const double bJetCSVDiscMin_;
		const double bJetCSVDiscMax_;
		const double fatJetPtMin_ ;
		const double fatJetPtMax_ ; 
		const double fatJetPrunedMassMin_ ;
		const double fatJetPrunedMassMax_ ; 
		const double dRSubjetsMin_;
		const double dRSubjetsMax_;
		const double subj1CSVDiscMin_;
		const double subj1CSVDiscMax_;
		const double subj2CSVDiscMin_;
		const double subj2CSVDiscMax_;
		const double HTAK5Min_;
		const double HTAK5Max_;

		const double     HJetPtMin_;
		const double     HJetPtMax_;
		const double     HJetAbsEtaMax_;
		const double     HJetMassMin_;	
		const double     HJetMassMax_;	
		const double     HJetPrunedMassMin_;	
		const double     HJetPrunedMassMax_;	
		const double     HJetTau2ByTau1Min_;
		const double     HJetTau2ByTau1Max_;
		const double     HJetSBMassMin_;	
		const double     HJetSBMassMax_;

		const double     bVetoJetPtMin_;
		const double     bVetoJetPtMax_;
		const double     bVetoJetAbsEtaMin_;
		const double     bVetoJetAbsEtaMax_;
		const double     bVetoJetCSVDiscMin_;
		const double     bVetoJetCSVDiscMax_;

		const int        numbJetMin_;
		const int        numHiggsJetMin_;

		const edm::ParameterSet jetSelParams_ ; 
		const edm::ParameterSet fatJetSelParams_ ; 
		const edm::ParameterSet higgsJetSelParams_ ; 
		const edm::ParameterSet HTSelParams_ ; 
		const edm::ParameterSet evtSelParams_ ; 
		const edm::ParameterSet jmeParams_; 
		const bool   applyJEC_ ; 
		const bool   applyBTagSF_ ; 
		const double jesShift_;
		const double jerShift_; 
		const double SFbShift_;
		const double SFlShift_;

		EvtInfoBranches    EvtInfo;
		VertexInfoBranches VtxInfo;
		GenInfoBranches    GenInfo;
		JetInfoBranches    GenJetInfo;
		JetInfoBranches    JetInfo;
		JetInfoBranches    FatJetInfo;
		JetInfoBranches    SubJetInfo;
		LepInfoBranches    LepInfo;

		edm::Service<TFileService> fs; 

		TH2InfoClass<TH2D> h2;
		TH1InfoClass<TH1D> h1;

		int	evtPass_ana;
		int	evtPass_val;

		bool BuildMiniTree_;
		TChain*            chain_;
		TTree*             newtree_ana;
		JetInfoBranches HiggsJetInfoAna, AntiHiggsJetInfoAna;
		JetInfoBranches bJetInfoAna, Final_bJetInfoAna;
		JetInfoBranches HiggsSubJet1InfoAna, HiggsSubJet2InfoAna, 
				AntiHiggsSubJet1InfoAna, AntiHiggsSubJet2InfoAna;
		bool   McFlagana;
		double PUana;
		double evtWtana;
		double HTak5, HThiggsbjet;

};

//
// constructors and destructor
//
BackgroundEstimationABCD::BackgroundEstimationABCD(const edm::ParameterSet& iConfig) : 
	maxEvents_(iConfig.getParameter<int>("MaxEvents")), 
	reportEvery_(iConfig.getParameter<int>("ReportEvery")),
	inputTTree_(iConfig.getParameter<std::string>("InputTTree")),
	inputFiles_(iConfig.getParameter<std::vector<std::string> >("InputFiles")),
	hltPaths_(iConfig.getParameter<edm::ParameterSet>("HLTPaths")),
	doPUReweighting_(iConfig.getParameter<bool>("DoPUReweighting")),
	file_PUDistMC_(iConfig.getParameter<std::string>("File_PUDistMC")),
	file_PUDistData_(iConfig.getParameter<std::string>("File_PUDistData")),
	hist_PUDistMC_(iConfig.getParameter<std::string>("Hist_PUDistMC")),
	hist_PUDistData_(iConfig.getParameter<std::string>("Hist_PUDistData")),
	jetPtMin_(iConfig.getParameter<double>("JetPtMin")),
	bJetPtMin_(iConfig.getParameter<double>("BJetPtMin")),
	bJetCSVDiscMin_(iConfig.getParameter<double>("BJetCSVDiscMin")),
	bJetCSVDiscMax_(iConfig.getParameter<double>("BJetCSVDiscMax")),
	fatJetPtMin_(iConfig.getParameter<double>("FatJetPtMin")),
	fatJetPtMax_(iConfig.getParameter<double>("FatJetPtMax")), 
	fatJetPrunedMassMin_(iConfig.getParameter<double>("FatJetPrunedMassMin")),
	fatJetPrunedMassMax_(iConfig.getParameter<double>("FatJetPrunedMassMax")),
	dRSubjetsMin_(iConfig.getParameter<double>("DRSubjetsMin")),
	dRSubjetsMax_(iConfig.getParameter<double>("DRSubjetsMax")),
	subj1CSVDiscMin_(iConfig.getParameter<double>("Subjet1CSVDiscMin")),
	subj1CSVDiscMax_(iConfig.getParameter<double>("Subjet1CSVDiscMax")),
	subj2CSVDiscMin_(iConfig.getParameter<double>("Subjet2CSVDiscMin")),
	subj2CSVDiscMax_(iConfig.getParameter<double>("Subjet2CSVDiscMax")),
	HTAK5Min_(iConfig.getParameter<double>("HTAK5Min")),
	HTAK5Max_(iConfig.getParameter<double>("HTAK5Max")),
	HJetPtMin_(iConfig.getParameter<double>("HJetPtMin")),
	HJetPtMax_(iConfig.getParameter<double>("HJetPtMax")),
	HJetAbsEtaMax_(iConfig.getParameter<double>("HJetAbsEtaMax")),
	HJetMassMin_(iConfig.getParameter<double>("HJetMassMin")),	
	HJetMassMax_(iConfig.getParameter<double>("HJetMassMax")),	
	HJetPrunedMassMin_(iConfig.getParameter<double>("HJetPrunedMassMin")),	
	HJetPrunedMassMax_(iConfig.getParameter<double>("HJetPrunedMassMax")),	
	HJetTau2ByTau1Min_(iConfig.getParameter<double>("HJetTau2ByTau1Min")),
	HJetTau2ByTau1Max_(iConfig.getParameter<double>("HJetTau2ByTau1Max")),
	HJetSBMassMin_(iConfig.getParameter<double>("HJetSBMassMin")),	
	HJetSBMassMax_(iConfig.getParameter<double>("HJetSBMassMax")),
	bVetoJetPtMin_(iConfig.getParameter<double>("bVetoJetPtMin")),
	bVetoJetPtMax_(iConfig.getParameter<double>("bVetoJetPtMax")),
	bVetoJetAbsEtaMin_(iConfig.getParameter<double>("bVetoJetAbsEtaMin")),
	bVetoJetAbsEtaMax_(iConfig.getParameter<double>("bVetoJetAbsEtaMax")),
	bVetoJetCSVDiscMin_(iConfig.getParameter<double>("bVetoJetCSVDiscMin")),
	bVetoJetCSVDiscMax_(iConfig.getParameter<double>("bVetoJetCSVDiscMax")),
	numbJetMin_(iConfig.getParameter<int>("numbJetMin")),
	numHiggsJetMin_(iConfig.getParameter<int>("numHiggsJetMin")),
	jetSelParams_(iConfig.getParameter<edm::ParameterSet>("JetSelParams")), 
	fatJetSelParams_(iConfig.getParameter<edm::ParameterSet>("FatJetSelParams")), 
	higgsJetSelParams_(iConfig.getParameter<edm::ParameterSet>("HiggsJetSelParams")), 
	evtSelParams_(iConfig.getParameter<edm::ParameterSet>("EvtSelParams")),
	jmeParams_(iConfig.getParameter<edm::ParameterSet>("JMEParams")),
	applyJEC_(iConfig.getParameter<bool>("ApplyJEC")),
	applyBTagSF_(iConfig.getParameter<bool>("ApplyBTagSF")),
	jesShift_(iConfig.getParameter<double>("JESShift")),
	jerShift_(iConfig.getParameter<double>("JERShift")),
	SFbShift_(iConfig.getParameter<double>("SFbShift")),
	SFlShift_(iConfig.getParameter<double>("SFlShift")),
	evtPass_ana(0),  
	evtPass_val(0), 
	BuildMiniTree_(iConfig.getParameter<bool>("BuildMinTree")) 
{ 
	if( doPUReweighting_) LumiWeights_ = edm::LumiReWeighting(file_PUDistMC_, file_PUDistData_, hist_PUDistMC_, hist_PUDistData_);
}


BackgroundEstimationABCD::~BackgroundEstimationABCD(){ 
	delete chain_;
}

// ------------ Other function -------------
template <typename Type>
void BackgroundEstimationABCD::setABCDcutRegion(Type* h1){
	h1->GetXaxis()->SetBinLabel(1,"A");
	h1->GetXaxis()->SetBinLabel(2,"B");
	h1->GetXaxis()->SetBinLabel(3,"C");
	h1->GetXaxis()->SetBinLabel(4,"D");
}

Jet BackgroundEstimationABCD::NearestJet( JetCollection jets, const Jet myJet ){
	Jet min_Jet;
	for( JetCollection::const_iterator ijet = jets.begin(); ijet != jets.end(); ++ijet ){
		float DR_min = 1000.;
		if( myJet.DeltaR(*ijet) < DR_min){
			DR_min = myJet.DeltaR(*ijet);
			min_Jet = *ijet;
		}
	}
	if( jets.size() != 0 ){
		return min_Jet;
	}else{
		edm::LogError("BackgroundEstimationABCD::NearestJet") << "JetCollection size is 0 !!!" ; 
		Jet emptyJet; emptyJet.Set_Pt(-1);
		return emptyJet;
	}
}

void BackgroundEstimationABCD::recoBprime( JetCollection bjets, const Jet H, vector<TLorentzVector>& p4_output ){
	Jet b = *(H.NearestJet( bjets )) ;
	if( b.Pt() == -1 ) return;
	TLorentzVector p4_bp, p4_H, p4_b;
	p4_H.SetPtEtaPhiM(H.Pt(), H.Eta(), H.Phi(), H.Mass());
	p4_b.SetPtEtaPhiM(b.Pt(), b.Eta(), b.Phi(), b.Mass());
	p4_bp = p4_H + p4_b;
	p4_output.push_back( p4_bp );
}

void BackgroundEstimationABCD::isolateCollection( JetCollection control, JetCollection input, JetCollection& output, double dR ){
	for( JetCollection::const_iterator ii = input.begin(); ii != input.end(); ++ii ){
		bool isControl(false); 
		for( JetCollection::const_iterator ic = control.begin(); ic != control.end(); ++ic ){
			if( ic->DeltaR(*ii) < dR ){
				isControl = true; 
				break; 
			}else {
				isControl = false; 
			} 
		}
		if( !isControl ){
			output.push_back(*ii);
		} 
	}
}

// ------------ method called once each job just before starting event loop  ------------
void BackgroundEstimationABCD::beginJob(){ 
	chain_  = new TChain(inputTTree_.c_str());

	for(unsigned i=0; i<inputFiles_.size(); ++i){
		chain_->Add(inputFiles_.at(i).c_str());
		TFile *f = TFile::Open(inputFiles_.at(i).c_str(),"READ");
		f->Close();
	}

	EvtInfo.Register(chain_);
	VtxInfo.Register(chain_);
	GenInfo.Register(chain_);
	GenJetInfo.Register(chain_,"GenJetInfo");
	JetInfo.Register(chain_,"JetInfo");
	FatJetInfo.Register(chain_,"FatJetInfo");
	SubJetInfo.Register(chain_,"SubJetInfo");
	LepInfo.Register(chain_);

	if(  maxEvents_<0 || maxEvents_>chain_->GetEntries()) maxEvents_ = chain_->GetEntries();

	if( BuildMiniTree_ ){
		newtree_ana = fs->make<TTree>("tree", "");
		newtree_ana->Branch("EvtInfo.McFlag", 		&McFlagana, 	"EvtInfo.McFlag/O"); // store weight of evt and pu for each event
		newtree_ana->Branch("EvtInfo.PU", 		&PUana, 	"EvtInfo.PU/D"); // store weight of evt and pu for each event
		newtree_ana->Branch("EvtInfo.WeightEvt",	&evtWtana, 	"EvtInfo.WeightEvt/D"); 	
		newtree_ana->Branch("EvtInfo.HT_AK5",		&HTak5, 	"EvtInfo.HT_AK5/D"); 	
		newtree_ana->Branch("EvtInfo.HT_HiggsbJets",	&HThiggsbjet, 	"EvtInfo.HT_HiggsbJets/D"); 	
		HiggsJetInfoAna.RegisterTree(newtree_ana,"HiggsJetInfo");
		AntiHiggsJetInfoAna.RegisterTree(newtree_ana,"AntiHiggsJetInfo");
		bJetInfoAna.RegisterTree(newtree_ana,"bJetInfo");
		Final_bJetInfoAna.RegisterTree(newtree_ana,"FinalbJetInfo");
		HiggsSubJet1InfoAna.RegisterTree(newtree_ana,"HiggsSubJet1Info");
		HiggsSubJet2InfoAna.RegisterTree(newtree_ana,"HiggsSubJet2Info");
		AntiHiggsSubJet1InfoAna.RegisterTree(newtree_ana,"AntiHiggsSubJet1Info");
		AntiHiggsSubJet2InfoAna.RegisterTree(newtree_ana,"AntiHiggsSubJet2Info");
	}

	h2.CreateTH2(fs); h2.Sumw2();
	h1.CreateTH1(fs); h1.Sumw2();

	h1.GetTH1("ABCDana_CutFlow")->GetXaxis()->SetBinLabel(1,"All_Evt");	
	h1.GetTH1("ABCDana_CutFlow")->GetXaxis()->SetBinLabel(2,"Trigger_Sel");	
	h1.GetTH1("ABCDana_CutFlow")->GetXaxis()->SetBinLabel(3,"Vertex_Sel");	
	h1.GetTH1("ABCDana_CutFlow")->GetXaxis()->SetBinLabel(4,"HT_Sel");	
	h1.GetTH1("ABCDana_CutFlow")->GetXaxis()->SetBinLabel(5,"ABCD_Sel");

	h1.GetTH1("ABCDval_CutFlow")->GetXaxis()->SetBinLabel(1,"All_Evt");	
	h1.GetTH1("ABCDval_CutFlow")->GetXaxis()->SetBinLabel(2,"Trigger_Sel");	
	h1.GetTH1("ABCDval_CutFlow")->GetXaxis()->SetBinLabel(3,"Vertex_Sel");	
	h1.GetTH1("ABCDval_CutFlow")->GetXaxis()->SetBinLabel(4,"HT_Sel");	
	h1.GetTH1("ABCDval_CutFlow")->GetXaxis()->SetBinLabel(5,"ABCD_Sel");	

	setABCDcutRegion(h1.GetTH1("ABCDana_CutRegion"));	
	setABCDcutRegion(h1.GetTH1("ABCDana_CutRegion_1b"));	
	setABCDcutRegion(h1.GetTH1("ABCDana_CutRegion_2b"));	
	setABCDcutRegion(h1.GetTH1("ABCDval_CutRegion"));
	setABCDcutRegion(h1.GetTH1("ABCDval_CutRegion_0ak5"));	
	setABCDcutRegion(h1.GetTH1("ABCDval_CutRegion_1ak5"));	
	setABCDcutRegion(h1.GetTH1("ABCDval_CutRegion_2ak5"));	

	return;  
}

// ------------ method called for each event  ------------
void BackgroundEstimationABCD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){ 
	using namespace edm;
	using namespace std;

	if(  chain_ == 0) return;

	FatJetSelector jetSelCA8(fatJetSelParams_) ; 
	FatJetSelector jetSelHiggs(higgsJetSelParams_) ; 
	JetSelector jetSelAK5(jetSelParams_) ; 
	pat::strbitset retjetidak5 = jetSelAK5.getBitTemplate() ; 
	pat::strbitset retjetidca8 = jetSelCA8.getBitTemplate() ; 

	edm::LogInfo("StartingAnalysisLoop") << "Starting analysis loop\n";

	for(int entry=0; entry<maxEvents_; entry++){
		if( (entry%reportEvery_) == 0) edm::LogInfo("Event") << entry << " of " << maxEvents_;

		//// Event variables 
		bool passHLT(false); 
		int nGoodVtxs(0);
		bool isdata(0);
		double evtwt(1); 
		double puwt(1); 
		double sumw2_a, sumw2_b, sumw2_c, sumw2_d; 
		double sumw2_a_1b, sumw2_b_1b, sumw2_c_1b, sumw2_d_1b; 
		double sumw2_a_2b, sumw2_b_2b, sumw2_c_2b, sumw2_d_2b; 
		double sumw2_av, sumw2_bv, sumw2_cv, sumw2_dv;
		sumw2_a = sumw2_b = sumw2_c = sumw2_d = sumw2_a_1b = sumw2_b_1b = sumw2_c_1b = sumw2_d_1b = sumw2_a_2b = sumw2_b_2b = sumw2_c_2b = sumw2_d_2b = sumw2_av= sumw2_bv= sumw2_cv= sumw2_dv = 0.;

		chain_->GetEntry(entry);

		isdata  = EvtInfo.McFlag ? 0 : 1; 
		if( !isdata ) evtwt = EvtInfo.Weight; 
		if( doPUReweighting_ && !isdata ) puwt = LumiWeights_.weight(EvtInfo.TrueIT[0]); 

		if ( !isdata ) { //// Gen info 
			std::vector<int> gen_higgs_indices ; 
			for (int igen = 0; igen < GenInfo.Size; ++igen) { 
				if ( GenInfo.Status[igen] == 3 && abs(GenInfo.PdgID[igen]) == 25 && GenInfo.nDa[igen] >= 2  ) gen_higgs_indices.push_back(igen) ; 
			}
			//// Higgs BR reweighting
			for (std::vector<int>::const_iterator ihig = gen_higgs_indices.begin(); ihig != gen_higgs_indices.end(); ++ihig) {
				int higgsDau0(abs(GenInfo.Da0PdgID[*ihig])), higgsDau1(abs(GenInfo.Da1PdgID[*ihig])) ;
				int higgsDau = higgsDau0+higgsDau1 ; 
				if(higgsDau==10) evtwt *= HiggsBRscaleFactors::higgsBBSf;
				else if(higgsDau==30) evtwt *= HiggsBRscaleFactors::higgsTauTauSf;
				else if(higgsDau==26) evtwt *= HiggsBRscaleFactors::higgsMuMuSf;
				else if(higgsDau==8)  evtwt *= HiggsBRscaleFactors::higgsCCSf;
				else if(higgsDau==6)  evtwt *= HiggsBRscaleFactors::higgsSSSf;
				else if(higgsDau==16) evtwt *= HiggsBRscaleFactors::higgsTTSf;
				else if(higgsDau==42) evtwt *= HiggsBRscaleFactors::higgsGGSf;
				else if(higgsDau==44) evtwt *= HiggsBRscaleFactors::higgsGammaGammaSf;
				else if(higgsDau==45) evtwt *= HiggsBRscaleFactors::higgsZGammaSf;
				else if(higgsDau==48) evtwt *= HiggsBRscaleFactors::higgsWWSf;
				else if(higgsDau==46) evtwt *= HiggsBRscaleFactors::higgsZZSf; 
				else evtwt *= 1 ; 
			} //// Higgs BR reweighting  
		} //// Gen info

		evtwt *= puwt ;
		double evtwt2 = evtwt*evtwt;

		h1.GetTH1("ABCDana_CutFlow")->Fill(0);	
		h1.GetTH1("ABCDval_CutFlow")->Fill(0);

		//// Trigger selection 
		TriggerSelector trigSel(hltPaths_); 
		passHLT = trigSel.getTrigDecision(EvtInfo); 
		if( !passHLT ) continue; 
		h1.GetTH1("ABCDana_CutFlow")->Fill(double(1),evtwt);	
		h1.GetTH1("ABCDval_CutFlow")->Fill(double(1),evtwt);

		//// Vertex selection 
		VertexSelector vtxSel(VtxInfo); 
		nGoodVtxs = vtxSel.NGoodVtxs(); 
		if( nGoodVtxs < 1){ edm::LogInfo("NoGoodPrimaryVertex") << " No good primary vertex "; continue; }
		h1.GetTH1("ABCDana_CutFlow")->Fill(double(2),evtwt);	
		h1.GetTH1("ABCDval_CutFlow")->Fill(double(2),evtwt);

		JetCollection fatjets ; 
		for (int ifatjet = 0; ifatjet < FatJetInfo.Size; ++ifatjet) { 
			if (jetSelCA8(FatJetInfo, ifatjet, SubJetInfo, retjetidca8) == 0) continue ; //// pt > 150 GeV, |eta| < 2.4 tau2/tau1 < 0.5 
			Jet thisjet(FatJetInfo, ifatjet) ; 
			fatjets.push_back(thisjet) ; 
		} //// Loop over fat jets 

		//// Apply JEC and b-tagging SFs for CA8 jets in MC  
		if ( !isdata && applyJEC_ ) {
			JMEUncertUtil* jmeUtil_jer = new JMEUncertUtil(jmeParams_, fatjets, "JERCA8MC", jerShift_) ; 
			JetCollection ca8jets_jer = jmeUtil_jer->GetModifiedJetColl() ; 
			delete jmeUtil_jer ; 
			fatjets.clear() ; 

			if ( abs(jesShift_) > 1E-6 ) {
				JMEUncertUtil* jmeUtil_jes = new JMEUncertUtil(jmeParams_, ca8jets_jer, "JESCA8MC", jesShift_) ; 
				JetCollection ca8jets_jes = jmeUtil_jes->GetModifiedJetColl() ; 
				delete jmeUtil_jes ; 

				for (JetCollection::const_iterator ijet = ca8jets_jes.begin(); ijet != ca8jets_jes.end(); ++ijet) {
					Jet thisjet(*ijet) ; 
					fatjets.push_back(thisjet) ; 
				}
			}
			else {
				for (JetCollection::const_iterator ijet = ca8jets_jer.begin(); ijet != ca8jets_jer.end(); ++ijet) {
					Jet thisjet(*ijet) ; 
					fatjets.push_back(thisjet) ; 
				}
			}
		} //// Apply JEC and b-tagging SFs for CA8 jets in MC 

		JetCollection HiggsJets, HiggsLikeJets, AllHiggsJets, AllHiggsJetsDRCut;
		JetCollection AntiHiggsJets, AntiHiggsLikeJets, AllAntiHiggsJets, AllAntiHiggsJetsDRCut;
		//// Fill Higgs, anti-Higgs jet collections from pre-selected fat jets  
		for (JetCollection::const_iterator ifat = fatjets.begin(); ifat != fatjets.end(); ++ifat) {
			Jet thisjet(*ifat) ; 
			if ( thisjet.Pt() < fatJetPtMin_ || thisjet.Pt() > fatJetPtMax_ ) continue ; //// Apply fat jet pT cut  
			Jet subjet1(SubJetInfo, ifat->Jet_SubJet1Idx()) ;
			Jet subjet2(SubJetInfo, ifat->Jet_SubJet2Idx()) ;
			if (subjet1.CombinedSVBJetTags() < 0.244 || subjet2.CombinedSVBJetTags() < 0.244) continue ;  
			double subjet_dyphi = subjet1.DeltaR(subjet2) ; 

			if (subjet1.CombinedSVBJetTags() <= subj1CSVDiscMin_ && subjet2.CombinedSVBJetTags() <= subj2CSVDiscMin_) { //// subjet disc
				AllAntiHiggsJets.push_back(thisjet); 
				if (thisjet.MassPruned() > fatJetPrunedMassMin_ && thisjet.MassPruned() < fatJetPrunedMassMax_ ) { //// fat jet pruned mass 
					if( subjet_dyphi >= dRSubjetsMin_ && subjet_dyphi <= dRSubjetsMax_  ){
						AntiHiggsJets.push_back(thisjet);
						AllAntiHiggsJetsDRCut.push_back(thisjet);
					}
				}
				else {
					AntiHiggsLikeJets.push_back(thisjet);
					AllAntiHiggsJetsDRCut.push_back(thisjet);
				}
			}
			else if (subjet1.CombinedSVBJetTags() > subj1CSVDiscMin_ && subjet2.CombinedSVBJetTags() > subj2CSVDiscMin_) {
				AllHiggsJets.push_back(thisjet);
				if (thisjet.MassPruned() > HJetPrunedMassMin_ && thisjet.MassPruned() < HJetPrunedMassMax_ ) { //// fat jet pruned mass 
					if( subjet_dyphi >= dRSubjetsMin_ && subjet_dyphi <= dRSubjetsMax_  ){
						if ( !isdata && applyBTagSF_ ) { //// Apply Higgs-tagging scale factor 
							ApplyHiggsTagSF* higgsTagSF = new ApplyHiggsTagSF(double(subjet1.Pt()), double(subjet2.Pt()), 
									double(subjet1.Eta()), double(subjet2.Eta()),
									subjet1.GenFlavor(), subjet2.GenFlavor(), 
									subjet1.CombinedSVBJetTags(), subjet2.CombinedSVBJetTags()) ; 
							evtwt *= higgsTagSF->GetHiggsTagSF() ;
							delete higgsTagSF ; 
						} //// Apply Higgs-tagging scale factor  
						HiggsJets.push_back(thisjet);
						AllHiggsJetsDRCut.push_back(thisjet);
					}
				}
				else {
					HiggsLikeJets.push_back(thisjet);
					AllHiggsJetsDRCut.push_back(thisjet);
				}
			}

		} //// Fill Higgs, anti-Higgs jet collections from pre-selected fat jets 

		JetCollection AllHiggsAntiHiggsJets(AllHiggsJets.begin(), AllHiggsJets.end()) ; 
		AllHiggsAntiHiggsJets.insert(AllHiggsAntiHiggsJets.end(), AllAntiHiggsJets.begin(), AllAntiHiggsJets.end()) ; 

		JetCollection allAK5Jets, cleanedAK5Jets, ak5JetsForHT;  
		for (int ijet = 0; ijet < JetInfo.Size; ++ijet) {
			retjetidak5.set(false) ;
			if (jetSelAK5(JetInfo, ijet,retjetidak5) == 0) continue ; 
			Jet thisjet(JetInfo, ijet) ;
			allAK5Jets.push_back(thisjet) ; 
		}

		//// Apply JEC and b-tagging SFs to MC 
		if ( !isdata ) {
			if ( applyJEC_ ) { //// Apply JEC for MC   

				//// All AK5 jets 
				JMEUncertUtil* jmeUtil_jer = new JMEUncertUtil(jmeParams_, allAK5Jets, "JERAK5MC", jerShift_) ; 
				JetCollection allAK5JetsJER = jmeUtil_jer->GetModifiedJetColl() ; 
				delete jmeUtil_jer ; 
				allAK5Jets.clear() ; 

				if ( abs(jesShift_) > 1E-6 ) {
					JMEUncertUtil* jmeUtil_jes = new JMEUncertUtil(jmeParams_, allAK5JetsJER, "JESAK5MC", jesShift_) ; 
					JetCollection allAK5JetsJES = jmeUtil_jes->GetModifiedJetColl() ; 
					delete jmeUtil_jes ; 

					for (JetCollection::const_iterator ijet = allAK5JetsJES.begin(); ijet != allAK5JetsJES.end(); ++ijet) {
						Jet thisjet(*ijet) ; 
						allAK5Jets.push_back(thisjet) ; 
					}
				}
				else {
					for (JetCollection::const_iterator ijet = allAK5JetsJER.begin(); ijet != allAK5JetsJER.end(); ++ijet) {
						Jet thisjet(*ijet) ; 
						allAK5Jets.push_back(thisjet) ; 
					}
				}

			} //// Apply JEC for MC 

			if ( applyBTagSF_  ) { //// Apply b-tagging SF for MC  
				ApplyBTagSF * btagsf =  new ApplyBTagSF(allAK5Jets, 0.679, "CSVM", SFbShift_, SFlShift_) ;  
				allAK5Jets.clear() ; 
				allAK5Jets =  btagsf->getBtaggedJetsWithSF () ; 
				delete btagsf ; 
			} //// Apply b-tagging SF for MC 
		} //// Apply JEC and b-tagging SFs to MC 

		for (JetCollection::const_iterator ijet = allAK5Jets.begin(); ijet != allAK5Jets.end(); ++ijet) {
			if (ijet->Pt() > jetPtMin_ ) ak5JetsForHT.push_back(*ijet) ; 
		}

		isolateCollection (AllHiggsAntiHiggsJets, allAK5Jets, cleanedAK5Jets) ; 

		JetCollection allBJets, selectedBJets, selectedAK5Jets ;  
		for (JetCollection::const_iterator ijet = cleanedAK5Jets.begin(); ijet != cleanedAK5Jets.end(); ++ijet) {
			if (ijet->Pt() > bVetoJetPtMin_ && ijet->CombinedSVBJetTags() > bVetoJetCSVDiscMin_) allBJets.push_back(*ijet) ; 
			if (ijet->Pt() > bJetPtMin_ && ijet->CombinedSVBJetTags() > bJetCSVDiscMin_ ) selectedBJets.push_back(*ijet) ; 
			if (ijet->Pt() > bJetPtMin_ ) selectedAK5Jets.push_back(*ijet) ; 
		}

		HT HTAllAK5, MyHT ; 
		HTAllAK5.setJetCollection(ak5JetsForHT) ;  
		HTAllAK5.buildHT() ; 

		MyHT.setJetCollection(HiggsJets);
		MyHT.setJetCollection(selectedBJets) ;
		MyHT.buildHT() ; 

		///// Fill evt and ABCD plots 
		int nA=0, nAv=0;
		int nB=0, nBv=0;
		int nC=0, nCv=0;
		int nD=0, nDv=0;

		if( HTAllAK5.getHT() < HTAK5Min_ ) continue;
		h1.GetTH1("ABCDana_CutFlow")->Fill(double(3),evtwt);	
		h1.GetTH1("ABCDval_CutFlow")->Fill(double(3),evtwt);

		std::vector<TLorentzVector> p4_bprimes_A, p4_bprimes_B, p4_bprimes_C, p4_bprimes_D; 
		JetCollection Final_bJets_ABCD;
		JetCollection HiggsJets_ABCD, HiggsSubJet1_ABCD, HiggsSubJet2_ABCD;
		JetCollection AntiHiggsJets_ABCD, AntiHiggsSubJet1_ABCD, AntiHiggsSubJet2_ABCD;

		//// Higgs region
		for( JetCollection::const_iterator iH = AllHiggsJetsDRCut.begin(); iH != AllHiggsJetsDRCut.end(); iH++ ){
			int iSub1, iSub2;
			if( SubJetInfo.Pt[iH->Jet_SubJet1Idx()] > SubJetInfo.Pt[iH->Jet_SubJet2Idx()] ){
				iSub1 = iH->Jet_SubJet1Idx();
				iSub2 = iH->Jet_SubJet2Idx();
			}else{
				iSub1 = iH->Jet_SubJet2Idx();
				iSub2 = iH->Jet_SubJet1Idx();
			}
			TLorentzVector Subjet1, Subjet2;
			Subjet1.SetPtEtaPhiM(SubJetInfo.Pt[iSub1], SubJetInfo.Eta[iSub1], SubJetInfo.Phi[iSub1], SubJetInfo.Mass[iSub1]);
			Subjet2.SetPtEtaPhiM(SubJetInfo.Pt[iSub2], SubJetInfo.Eta[iSub2], SubJetInfo.Phi[iSub2], SubJetInfo.Mass[iSub2]);
			double subjet_dy = Subjet1.Rapidity() - Subjet2.Rapidity();
			double subjet_dphi = Subjet1.DeltaPhi(Subjet2);
			double subjet_dyphi = sqrt( subjet_dy*subjet_dy + subjet_dphi*subjet_dphi );

			Jet SubJet1(SubJetInfo, iSub1);
			Jet SubJet2(SubJetInfo, iSub2);

			//// B region 
			if( iH->MassPruned() > HJetPrunedMassMin_ && iH->MassPruned() < HJetPrunedMassMax_ ){
				if ( selectedBJets.size() >= unsigned(numbJetMin_) ){
					nB++; 
					HiggsJets_ABCD.push_back(*iH);
					HiggsSubJet1_ABCD.push_back(SubJet1);
					HiggsSubJet2_ABCD.push_back(SubJet2);
					Final_bJets_ABCD.push_back( *(iH->NearestJet(selectedBJets)) );
					recoBprime( selectedBJets, *iH, p4_bprimes_B);
					h1.GetTH1("ABCDana_HiggsMass")->Fill( iH->MassPruned(), evtwt);
					h1.GetTH1("ABCDana_HiggsMass_B")->Fill( iH->MassPruned(), evtwt);
					h1.GetTH1("ABCDana_CA8Pt")->Fill( iH->Pt(), evtwt);
					h1.GetTH1("ABCDana_CA8Pt_B")->Fill( iH->Pt(), evtwt);
					h1.GetTH1("ABCDana_dRSubJets_B")->Fill( subjet_dyphi, evtwt);
					h1.GetTH1("ABCDana_Tau2ByTau1_B")->Fill( iH->tau2()/iH->tau1(), evtwt);
					h1.GetTH1("ABCDana_Sub1Mass_B")->Fill( SubJetInfo.Mass[iSub1], evtwt);
					h1.GetTH1("ABCDana_Sub2Mass_B")->Fill( SubJetInfo.Mass[iSub2], evtwt);
					h1.GetTH1("ABCDana_Sub1Pt_B")->Fill( SubJetInfo.Pt[iSub1], evtwt);
					h1.GetTH1("ABCDana_Sub2Pt_B")->Fill( SubJetInfo.Pt[iSub2], evtwt);
					h1.GetTH1("ABCDana_Sub1CSV_B")->Fill( SubJetInfo.CombinedSVBJetTags[iSub1], evtwt);
					h1.GetTH1("ABCDana_Sub2CSV_B")->Fill( SubJetInfo.CombinedSVBJetTags[iSub2], evtwt);
					h2.GetTH2("ABCDana_2D")->Fill( 1., iH->MassPruned(), evtwt);
				}
				if ( allBJets.size() == 0 ){ //b-Veto
					nBv++; 
					h1.GetTH1("ABCDval_HiggsMass")->Fill( iH->MassPruned(), evtwt);
					h1.GetTH1("ABCDval_HiggsMass_B")->Fill( iH->MassPruned(), evtwt);
					h1.GetTH1("ABCDval_CA8Pt")->Fill( iH->Pt(), evtwt);
					h1.GetTH1("ABCDval_CA8Pt_B")->Fill( iH->Pt(), evtwt);
					h1.GetTH1("ABCDval_dRSubJets_B")->Fill( subjet_dyphi, evtwt);
					h1.GetTH1("ABCDval_Tau2ByTau1_B")->Fill( iH->tau2()/iH->tau1(), evtwt);
					h1.GetTH1("ABCDval_Sub1Mass_B")->Fill( SubJetInfo.Mass[iSub1], evtwt);
					h1.GetTH1("ABCDval_Sub2Mass_B")->Fill( SubJetInfo.Mass[iSub2], evtwt);
					h1.GetTH1("ABCDval_Sub1Pt_B")->Fill( SubJetInfo.Pt[iSub1], evtwt);
					h1.GetTH1("ABCDval_Sub2Pt_B")->Fill( SubJetInfo.Pt[iSub2], evtwt);
					h1.GetTH1("ABCDval_Sub1CSV_B")->Fill( SubJetInfo.CombinedSVBJetTags[iSub1], evtwt);
					h1.GetTH1("ABCDval_Sub2CSV_B")->Fill( SubJetInfo.CombinedSVBJetTags[iSub2], evtwt);
					h2.GetTH2("ABCDval_2D")->Fill( 1., iH->MassPruned(), evtwt);
				}
			}
			//// D region 
			if( iH->MassPruned() <= HJetPrunedMassMin_ || iH->MassPruned() >= HJetPrunedMassMax_ ){
				if ( selectedBJets.size() >= unsigned(numbJetMin_) ){
					if( iH->MassPruned() <= HJetSBMassMax_ && iH->MassPruned() > HJetSBMassMin_ ){ 
						nD++; 
						HiggsJets_ABCD.push_back(*iH);
						HiggsSubJet1_ABCD.push_back(SubJet1);
						HiggsSubJet2_ABCD.push_back(SubJet2);
						Final_bJets_ABCD.push_back( *(iH->NearestJet(selectedBJets)) );
						recoBprime( selectedBJets, *iH, p4_bprimes_D);
					}
					h1.GetTH1("ABCDana_HiggsMass")->Fill( iH->MassPruned(), evtwt);
					h1.GetTH1("ABCDana_HiggsMass_D")->Fill( iH->MassPruned(), evtwt);
					h1.GetTH1("ABCDana_CA8Pt")->Fill( iH->Pt(), evtwt);
					h1.GetTH1("ABCDana_CA8Pt_D")->Fill( iH->Pt(), evtwt);
					h1.GetTH1("ABCDana_dRSubJets_D")->Fill( subjet_dyphi, evtwt);
					h1.GetTH1("ABCDana_Tau2ByTau1_D")->Fill( iH->tau2()/iH->tau1(), evtwt);
					h1.GetTH1("ABCDana_Sub1Mass_D")->Fill( SubJetInfo.Mass[iSub1], evtwt);
					h1.GetTH1("ABCDana_Sub2Mass_D")->Fill( SubJetInfo.Mass[iSub2], evtwt);
					h1.GetTH1("ABCDana_Sub1Pt_D")->Fill( SubJetInfo.Pt[iSub1], evtwt);
					h1.GetTH1("ABCDana_Sub2Pt_D")->Fill( SubJetInfo.Pt[iSub2], evtwt);
					h1.GetTH1("ABCDana_Sub1CSV_D")->Fill( SubJetInfo.CombinedSVBJetTags[iSub1], evtwt);
					h1.GetTH1("ABCDana_Sub2CSV_D")->Fill( SubJetInfo.CombinedSVBJetTags[iSub2], evtwt);
					h2.GetTH2("ABCDana_2D")->Fill( 1., iH->MassPruned(), evtwt);
				}
				if ( allBJets.size() == 0 ){ //b-Veto
					if( iH->MassPruned() <= HJetSBMassMax_ && iH->MassPruned() > HJetSBMassMin_ ){ 
						nDv++; 
					} 
					h1.GetTH1("ABCDval_HiggsMass")->Fill( iH->MassPruned(), evtwt);
					h1.GetTH1("ABCDval_HiggsMass_D")->Fill( iH->MassPruned(), evtwt);
					h1.GetTH1("ABCDval_CA8Pt")->Fill( iH->Pt(), evtwt);
					h1.GetTH1("ABCDval_CA8Pt_D")->Fill( iH->Pt(), evtwt);
					h1.GetTH1("ABCDval_dRSubJets_D")->Fill( subjet_dyphi, evtwt);
					h1.GetTH1("ABCDval_Tau2ByTau1_D")->Fill( iH->tau2()/iH->tau1(), evtwt);
					h1.GetTH1("ABCDval_Sub1Mass_D")->Fill( SubJetInfo.Mass[iSub1], evtwt);
					h1.GetTH1("ABCDval_Sub2Mass_D")->Fill( SubJetInfo.Mass[iSub2], evtwt);
					h1.GetTH1("ABCDval_Sub1Pt_D")->Fill( SubJetInfo.Pt[iSub1], evtwt);
					h1.GetTH1("ABCDval_Sub2Pt_D")->Fill( SubJetInfo.Pt[iSub2], evtwt);
					h1.GetTH1("ABCDval_Sub1CSV_D")->Fill( SubJetInfo.CombinedSVBJetTags[iSub1], evtwt);
					h1.GetTH1("ABCDval_Sub2CSV_D")->Fill( SubJetInfo.CombinedSVBJetTags[iSub2], evtwt);
					h2.GetTH2("ABCDval_2D")->Fill( 1., iH->MassPruned(), evtwt);
				}
			}
		}

		//// AntiHiggs region 
		for( JetCollection::const_iterator H = AllAntiHiggsJetsDRCut.begin(); H != AllAntiHiggsJetsDRCut.end(); H++ ){
			int iSub1, iSub2;
			if( SubJetInfo.Pt[H->Jet_SubJet1Idx()] > SubJetInfo.Pt[H->Jet_SubJet2Idx()] ){
				iSub1 = H->Jet_SubJet1Idx();
				iSub2 = H->Jet_SubJet2Idx();
			}else{
				iSub1 = H->Jet_SubJet2Idx();
				iSub2 = H->Jet_SubJet1Idx();
			}
			TLorentzVector Subjet1, Subjet2;
			Subjet1.SetPtEtaPhiM(SubJetInfo.Pt[iSub1], SubJetInfo.Eta[iSub1], SubJetInfo.Phi[iSub1], SubJetInfo.Mass[iSub1]);
			Subjet2.SetPtEtaPhiM(SubJetInfo.Pt[iSub2], SubJetInfo.Eta[iSub2], SubJetInfo.Phi[iSub2], SubJetInfo.Mass[iSub2]);
			double subjet_dy = Subjet1.Rapidity() - Subjet2.Rapidity();
			double subjet_dphi = Subjet1.DeltaPhi(Subjet2);
			double subjet_dyphi = sqrt( subjet_dy*subjet_dy + subjet_dphi*subjet_dphi );

			Jet SubJet1(SubJetInfo, iSub1);
			Jet SubJet2(SubJetInfo, iSub2);

			//// A region	
			if( H->MassPruned() > HJetPrunedMassMin_ && H->MassPruned() < HJetPrunedMassMax_ ){
				if ( selectedBJets.size() >= unsigned(numbJetMin_) ){
					nA++;
					AntiHiggsJets_ABCD.push_back(*H); 
					AntiHiggsSubJet1_ABCD.push_back(SubJet1);
					AntiHiggsSubJet2_ABCD.push_back(SubJet2);
					Final_bJets_ABCD.push_back( *(H->NearestJet(selectedBJets)) );
					recoBprime( selectedBJets, *H, p4_bprimes_A);
					h1.GetTH1("ABCDana_HiggsMass")->Fill( H->MassPruned(), evtwt);
					h1.GetTH1("ABCDana_HiggsMass_A")->Fill( H->MassPruned(), evtwt);
					h1.GetTH1("ABCDana_CA8Pt")->Fill( H->Pt(), evtwt);
					h1.GetTH1("ABCDana_CA8Pt_A")->Fill( H->Pt(), evtwt);
					h1.GetTH1("ABCDana_dRSubJets_A")->Fill( subjet_dyphi, evtwt);
					h1.GetTH1("ABCDana_Tau2ByTau1_A")->Fill( H->tau2()/H->tau1(), evtwt);
					h1.GetTH1("ABCDana_Sub1Mass_A")->Fill( SubJetInfo.Mass[iSub1], evtwt);
					h1.GetTH1("ABCDana_Sub2Mass_A")->Fill( SubJetInfo.Mass[iSub2], evtwt);
					h1.GetTH1("ABCDana_Sub1Pt_A")->Fill( SubJetInfo.Pt[iSub1], evtwt);
					h1.GetTH1("ABCDana_Sub2Pt_A")->Fill( SubJetInfo.Pt[iSub2], evtwt);
					h1.GetTH1("ABCDana_Sub1CSV_A")->Fill( SubJetInfo.CombinedSVBJetTags[iSub1], evtwt);
					h1.GetTH1("ABCDana_Sub2CSV_A")->Fill( SubJetInfo.CombinedSVBJetTags[iSub2], evtwt);
					h2.GetTH2("ABCDana_2D")->Fill( 0., H->MassPruned(), evtwt);
				}
				if ( allBJets.size() == 0 ){
					nAv++; 
					h1.GetTH1("ABCDval_HiggsMass")->Fill( H->MassPruned(), evtwt);
					h1.GetTH1("ABCDval_HiggsMass_A")->Fill( H->MassPruned(), evtwt);
					h1.GetTH1("ABCDval_CA8Pt")->Fill( H->Pt(), evtwt);
					h1.GetTH1("ABCDval_CA8Pt_A")->Fill( H->Pt(), evtwt);
					h1.GetTH1("ABCDval_dRSubJets_A")->Fill( subjet_dyphi, evtwt);
					h1.GetTH1("ABCDval_Tau2ByTau1_A")->Fill( H->tau2()/H->tau1(), evtwt);
					h1.GetTH1("ABCDval_Sub1Mass_A")->Fill( SubJetInfo.Mass[iSub1], evtwt);
					h1.GetTH1("ABCDval_Sub2Mass_A")->Fill( SubJetInfo.Mass[iSub2], evtwt);
					h1.GetTH1("ABCDval_Sub1Pt_A")->Fill( SubJetInfo.Pt[iSub1], evtwt);
					h1.GetTH1("ABCDval_Sub2Pt_A")->Fill( SubJetInfo.Pt[iSub2], evtwt);
					h1.GetTH1("ABCDval_Sub1CSV_A")->Fill( SubJetInfo.CombinedSVBJetTags[iSub1], evtwt);
					h1.GetTH1("ABCDval_Sub2CSV_A")->Fill( SubJetInfo.CombinedSVBJetTags[iSub2], evtwt);
					h2.GetTH2("ABCDval_2D")->Fill( 0., H->MassPruned(), evtwt);
				}
			}
			//// C region
			if( H->MassPruned() <= HJetPrunedMassMin_ || H->MassPruned() >= HJetPrunedMassMax_ ){
				if ( selectedBJets.size() >= unsigned(numbJetMin_) ){
					if( H->MassPruned() <= HJetSBMassMax_ && H->MassPruned() > HJetSBMassMin_ ){ 
						nC++; 
						AntiHiggsJets_ABCD.push_back(*H); 
						AntiHiggsSubJet1_ABCD.push_back(SubJet1);
						AntiHiggsSubJet2_ABCD.push_back(SubJet2);	
						Final_bJets_ABCD.push_back( *(H->NearestJet(selectedBJets)) );
						recoBprime( selectedBJets, *H, p4_bprimes_C);
					} 
					h1.GetTH1("ABCDana_HiggsMass")->Fill( H->MassPruned(), evtwt);
					h1.GetTH1("ABCDana_HiggsMass_C")->Fill( H->MassPruned(), evtwt);
					h1.GetTH1("ABCDana_CA8Pt")->Fill( H->Pt(), evtwt);
					h1.GetTH1("ABCDana_CA8Pt_C")->Fill( H->Pt(), evtwt);
					h1.GetTH1("ABCDana_dRSubJets_C")->Fill( subjet_dyphi, evtwt);
					h1.GetTH1("ABCDana_Tau2ByTau1_C")->Fill( H->tau2()/H->tau1(), evtwt);
					h1.GetTH1("ABCDana_Sub1Mass_C")->Fill( SubJetInfo.Mass[iSub1], evtwt);
					h1.GetTH1("ABCDana_Sub2Mass_C")->Fill( SubJetInfo.Mass[iSub2], evtwt);
					h1.GetTH1("ABCDana_Sub1Pt_C")->Fill( SubJetInfo.Pt[iSub1], evtwt);
					h1.GetTH1("ABCDana_Sub2Pt_C")->Fill( SubJetInfo.Pt[iSub2], evtwt);
					h1.GetTH1("ABCDana_Sub1CSV_C")->Fill( SubJetInfo.CombinedSVBJetTags[iSub1], evtwt);
					h1.GetTH1("ABCDana_Sub2CSV_C")->Fill( SubJetInfo.CombinedSVBJetTags[iSub2], evtwt);
					h2.GetTH2("ABCDana_2D")->Fill( 0., H->MassPruned(), evtwt);
				}
				if ( allBJets.size() == 0){ // b-Veto
					if( H->MassPruned() <= HJetSBMassMax_ && H->MassPruned() > HJetSBMassMin_ ){ 
						nCv++; 
					} 
					h1.GetTH1("ABCDval_HiggsMass")->Fill( H->MassPruned(), evtwt);
					h1.GetTH1("ABCDval_HiggsMass_C")->Fill( H->MassPruned(), evtwt);
					h1.GetTH1("ABCDval_CA8Pt")->Fill( H->Pt(), evtwt);
					h1.GetTH1("ABCDval_CA8Pt_C")->Fill( H->Pt(), evtwt);
					h1.GetTH1("ABCDval_dRSubJets_C")->Fill( subjet_dyphi, evtwt);
					h1.GetTH1("ABCDval_Tau2ByTau1_C")->Fill( H->tau2()/H->tau1(), evtwt);
					h1.GetTH1("ABCDval_Sub1Mass_C")->Fill( SubJetInfo.Mass[iSub1], evtwt);
					h1.GetTH1("ABCDval_Sub2Mass_C")->Fill( SubJetInfo.Mass[iSub2], evtwt);
					h1.GetTH1("ABCDval_Sub1Pt_C")->Fill( SubJetInfo.Pt[iSub1], evtwt);
					h1.GetTH1("ABCDval_Sub2Pt_C")->Fill( SubJetInfo.Pt[iSub2], evtwt);
					h1.GetTH1("ABCDval_Sub1CSV_C")->Fill( SubJetInfo.CombinedSVBJetTags[iSub1], evtwt);
					h1.GetTH1("ABCDval_Sub2CSV_C")->Fill( SubJetInfo.CombinedSVBJetTags[iSub2], evtwt);
					h2.GetTH2("ABCDval_2D")->Fill( 0., H->MassPruned(), evtwt);
				}
			}
		}

		h1.GetTH1("ABCDana_NumCA8")->Fill(AllHiggsJets.size()+AllAntiHiggsJets.size());
		h1.GetTH1("ABCDana_Numbjet")->Fill(selectedBJets.size());
		if( nA+nB+nC+nD > 0){ 
			h1.GetTH1("ABCDana_Numbjet_ABCD")->Fill(selectedBJets.size());
			h1.GetTH1("ABCDana_NumCA8_ABCD")->Fill(nA+nB+nC+nD);
		}

		if( nB > 0 ){
			sumw2_b += evtwt2;
			h1.GetTH1("ABCDana_CutRegion")->Fill("B", evtwt); evtPass_ana++;
			h1.GetTH1("ABCDana_HT_B")->Fill( HTAllAK5.getHT(), evtwt);
			h1.GetTH1("ABCDana_HT")->Fill( HTAllAK5.getHT(), evtwt);
			h1.GetTH1("ABCDana_NumCA8_B")->Fill(nB);
			h1.GetTH1("ABCDana_Numbjet_B")->Fill(selectedBJets.size());
			for( vector<TLorentzVector>::const_iterator bp_ = p4_bprimes_B.begin(); bp_ != p4_bprimes_B.end(); bp_++ ){
				h1.GetTH1("ABCDana_bpMass_B")->Fill( bp_->M(), evtwt);
				h1.GetTH1("ABCDana_bpPt_B")->Fill( bp_->Pt(), evtwt);
				h1.GetTH1("ABCDana_bpEta_B")->Fill( bp_->Eta(), evtwt);
			}
			if( selectedBJets.size() == 1 ){
				sumw2_b_1b += evtwt2;
				h1.GetTH1("ABCDana_CutRegion_1b")->Fill("B", evtwt); 
				h1.GetTH1("ABCDana_1b_HT_B")->Fill( HTAllAK5.getHT(), evtwt);
				for( vector<TLorentzVector>::const_iterator bp_ = p4_bprimes_B.begin(); bp_ != p4_bprimes_B.end(); bp_++ ){
					h1.GetTH1("ABCDana_1b_bpMass_B")->Fill( bp_->M(), evtwt);
					h1.GetTH1("ABCDana_1b_bpPt_B")->Fill( bp_->Pt(), evtwt);
					h1.GetTH1("ABCDana_1b_bpEta_B")->Fill( bp_->Eta(), evtwt);
				}
			}else if( selectedBJets.size()>= 2 ){
				sumw2_b_2b += evtwt2;
				h1.GetTH1("ABCDana_CutRegion_2b")->Fill("B", evtwt); 
				h1.GetTH1("ABCDana_2b_HT_B")->Fill( HTAllAK5.getHT(), evtwt);
				for( vector<TLorentzVector>::const_iterator bp_ = p4_bprimes_B.begin(); bp_ != p4_bprimes_B.end(); bp_++ ){
					h1.GetTH1("ABCDana_2b_bpMass_B")->Fill( bp_->M(), evtwt);
					h1.GetTH1("ABCDana_2b_bpPt_B")->Fill( bp_->Pt(), evtwt);
					h1.GetTH1("ABCDana_2b_bpEta_B")->Fill( bp_->Eta(), evtwt);
				}
			}	
		}else if( nD > 0 ){ 
			sumw2_d += evtwt2;
			h1.GetTH1("ABCDana_CutRegion")->Fill("D", evtwt);
			h1.GetTH1("ABCDana_HT_D")->Fill( HTAllAK5.getHT(), evtwt);
			h1.GetTH1("ABCDana_HT")->Fill( HTAllAK5.getHT(), evtwt);
			h1.GetTH1("ABCDana_NumCA8_D")->Fill(nD);
			h1.GetTH1("ABCDana_Numbjet_D")->Fill(selectedBJets.size());
			for( vector<TLorentzVector>::const_iterator bp_ = p4_bprimes_D.begin(); bp_ != p4_bprimes_D.end(); bp_++ ){
				h1.GetTH1("ABCDana_bpMass_D")->Fill( bp_->M(), evtwt);
				h1.GetTH1("ABCDana_bpPt_D")->Fill( bp_->Pt(), evtwt);
				h1.GetTH1("ABCDana_bpEta_D")->Fill( bp_->Eta(), evtwt);
			}
			if( selectedBJets.size() == 1 ){
				sumw2_d_1b += evtwt2;
				h1.GetTH1("ABCDana_CutRegion_1b")->Fill("D", evtwt); 
				h1.GetTH1("ABCDana_1b_HT_D")->Fill( HTAllAK5.getHT(), evtwt);
				for( vector<TLorentzVector>::const_iterator bp_ = p4_bprimes_D.begin(); bp_ != p4_bprimes_D.end(); bp_++ ){
					h1.GetTH1("ABCDana_1b_bpMass_D")->Fill( bp_->M(), evtwt);
					h1.GetTH1("ABCDana_1b_bpPt_D")->Fill( bp_->Pt(), evtwt);
					h1.GetTH1("ABCDana_1b_bpEta_D")->Fill( bp_->Eta(), evtwt);
				}
			}else if( selectedBJets.size()>= 2 ){
				sumw2_d_2b += evtwt2;
				h1.GetTH1("ABCDana_CutRegion_2b")->Fill("D", evtwt); 
				h1.GetTH1("ABCDana_2b_HT_D")->Fill( HTAllAK5.getHT(), evtwt);
				for( vector<TLorentzVector>::const_iterator bp_ = p4_bprimes_D.begin(); bp_ != p4_bprimes_D.end(); bp_++ ){
					h1.GetTH1("ABCDana_2b_bpMass_D")->Fill( bp_->M(), evtwt);
					h1.GetTH1("ABCDana_2b_bpPt_D")->Fill( bp_->Pt(), evtwt);
					h1.GetTH1("ABCDana_2b_bpEta_D")->Fill( bp_->Eta(), evtwt);
				}
			}	
		}else if( nA > 0 ){
			sumw2_a += evtwt2;
			h1.GetTH1("ABCDana_CutRegion")->Fill("A", evtwt);
			h1.GetTH1("ABCDana_HT_A")->Fill( HTAllAK5.getHT(), evtwt);
			h1.GetTH1("ABCDana_HT")->Fill( HTAllAK5.getHT(), evtwt);
			h1.GetTH1("ABCDana_NumCA8_A")->Fill(nA);
			h1.GetTH1("ABCDana_Numbjet_A")->Fill(selectedBJets.size());
			for( vector<TLorentzVector>::const_iterator bp_ = p4_bprimes_A.begin(); bp_ != p4_bprimes_A.end(); bp_++ ){
				h1.GetTH1("ABCDana_bpMass_A")->Fill( bp_->M(), evtwt);
				h1.GetTH1("ABCDana_bpPt_A")->Fill( bp_->Pt(), evtwt);
				h1.GetTH1("ABCDana_bpEta_A")->Fill( bp_->Eta(), evtwt);
			}
			if( selectedBJets.size() == 1 ){
				sumw2_a_1b += evtwt2;
				h1.GetTH1("ABCDana_CutRegion_1b")->Fill("A", evtwt); 
				h1.GetTH1("ABCDana_1b_HT_A")->Fill( HTAllAK5.getHT(), evtwt);
				for( vector<TLorentzVector>::const_iterator bp_ = p4_bprimes_A.begin(); bp_ != p4_bprimes_A.end(); bp_++ ){
					h1.GetTH1("ABCDana_1b_bpMass_A")->Fill( bp_->M(), evtwt);
					h1.GetTH1("ABCDana_1b_bpPt_A")->Fill( bp_->Pt(), evtwt);
					h1.GetTH1("ABCDana_1b_bpEta_A")->Fill( bp_->Eta(), evtwt);
				}
			}else if( selectedBJets.size()>= 2 ){
				sumw2_a_2b += evtwt2;
				h1.GetTH1("ABCDana_CutRegion_2b")->Fill("A", evtwt); 
				h1.GetTH1("ABCDana_2b_HT_A")->Fill( HTAllAK5.getHT(), evtwt);
				for( vector<TLorentzVector>::const_iterator bp_ = p4_bprimes_A.begin(); bp_ != p4_bprimes_A.end(); bp_++ ){
					h1.GetTH1("ABCDana_2b_bpMass_A")->Fill( bp_->M(), evtwt);
					h1.GetTH1("ABCDana_2b_bpPt_A")->Fill( bp_->Pt(), evtwt);
					h1.GetTH1("ABCDana_2b_bpEta_A")->Fill( bp_->Eta(), evtwt);
				}
			}	
		}else if( nC > 0 ){
			sumw2_c += evtwt2;
			h1.GetTH1("ABCDana_CutRegion")->Fill("C", evtwt);
			h1.GetTH1("ABCDana_HT_C")->Fill( HTAllAK5.getHT(), evtwt);
			h1.GetTH1("ABCDana_HT")->Fill( HTAllAK5.getHT(), evtwt);
			h1.GetTH1("ABCDana_NumCA8_C")->Fill(nC);
			h1.GetTH1("ABCDana_Numbjet_C")->Fill(selectedBJets.size());
			for( vector<TLorentzVector>::const_iterator bp_ = p4_bprimes_C.begin(); bp_ != p4_bprimes_C.end(); bp_++ ){
				h1.GetTH1("ABCDana_bpMass_C")->Fill( bp_->M(), evtwt);
				h1.GetTH1("ABCDana_bpPt_C")->Fill( bp_->Pt(), evtwt);
				h1.GetTH1("ABCDana_bpEta_C")->Fill( bp_->Eta(), evtwt);
			}
			if( selectedBJets.size() == 1 ){
				sumw2_c_1b += evtwt2;
				h1.GetTH1("ABCDana_CutRegion_1b")->Fill("C", evtwt); 
				h1.GetTH1("ABCDana_1b_HT_C")->Fill( HTAllAK5.getHT(), evtwt);
				for( vector<TLorentzVector>::const_iterator bp_ = p4_bprimes_C.begin(); bp_ != p4_bprimes_C.end(); bp_++ ){
					h1.GetTH1("ABCDana_1b_bpMass_C")->Fill( bp_->M(), evtwt);
					h1.GetTH1("ABCDana_1b_bpPt_C")->Fill( bp_->Pt(), evtwt);
					h1.GetTH1("ABCDana_1b_bpEta_C")->Fill( bp_->Eta(), evtwt);
				}
			}else if( selectedBJets.size()>= 2 ){
				sumw2_c_2b += evtwt2;
				h1.GetTH1("ABCDana_CutRegion_2b")->Fill("C", evtwt); 
				h1.GetTH1("ABCDana_2b_HT_C")->Fill( HTAllAK5.getHT(), evtwt);
				for( vector<TLorentzVector>::const_iterator bp_ = p4_bprimes_C.begin(); bp_ != p4_bprimes_C.end(); bp_++ ){
					h1.GetTH1("ABCDana_2b_bpMass_C")->Fill( bp_->M(), evtwt);
					h1.GetTH1("ABCDana_2b_bpPt_C")->Fill( bp_->Pt(), evtwt);
					h1.GetTH1("ABCDana_2b_bpEta_C")->Fill( bp_->Eta(), evtwt);
				}
			}	
		}

		if( nBv > 0 ){
			h1.GetTH1("ABCDval_CutRegion")->Fill("B", evtwt); evtPass_val++;
			h1.GetTH1("ABCDval_HT_B")->Fill( HTAllAK5.getHT(), evtwt);
			h1.GetTH1("ABCDval_HT")->Fill( HTAllAK5.getHT(), evtwt);
			h1.GetTH1("ABCDval_NumCA8_B")->Fill(nBv);
			sumw2_bv += evtwt2;
			if( selectedAK5Jets.size() == 0 ){
				h1.GetTH1("ABCDval_CutRegion_0ak5")->Fill("B", evtwt); 
				h1.GetTH1("ABCDval_0ak5_HT_B")->Fill( HTAllAK5.getHT(), evtwt);
			}else if( selectedAK5Jets.size() == 1 ){
				h1.GetTH1("ABCDval_CutRegion_1ak5")->Fill("B", evtwt); 
				h1.GetTH1("ABCDval_1ak5_HT_B")->Fill( HTAllAK5.getHT(), evtwt);
			}else if( selectedAK5Jets.size()>= 2 ){
				h1.GetTH1("ABCDval_CutRegion_2ak5")->Fill("B", evtwt); 
				h1.GetTH1("ABCDval_2ak5_HT_B")->Fill( HTAllAK5.getHT(), evtwt);
			}
		}else if( nDv > 0 ){ 
			h1.GetTH1("ABCDval_CutRegion")->Fill("D", evtwt);
			h1.GetTH1("ABCDval_HT_D")->Fill( HTAllAK5.getHT(), evtwt);
			h1.GetTH1("ABCDval_HT")->Fill( HTAllAK5.getHT(), evtwt);
			h1.GetTH1("ABCDval_NumCA8_D")->Fill(nDv);
			sumw2_dv += evtwt2;
			if( selectedAK5Jets.size() == 0 ){
				h1.GetTH1("ABCDval_CutRegion_0ak5")->Fill("D", evtwt); 
				h1.GetTH1("ABCDval_0ak5_HT_D")->Fill( HTAllAK5.getHT(), evtwt);
			}else if( selectedAK5Jets.size() == 1 ){
				h1.GetTH1("ABCDval_CutRegion_1ak5")->Fill("D", evtwt); 
				h1.GetTH1("ABCDval_1ak5_HT_D")->Fill( HTAllAK5.getHT(), evtwt);
			}else if( selectedAK5Jets.size()>= 2 ){
				h1.GetTH1("ABCDval_CutRegion_2ak5")->Fill("D", evtwt); 
				h1.GetTH1("ABCDval_2ak5_HT_D")->Fill( HTAllAK5.getHT(), evtwt);
			}
		}else if( nAv > 0 ){
			h1.GetTH1("ABCDval_CutRegion")->Fill("A", evtwt);
			h1.GetTH1("ABCDval_HT_A")->Fill( HTAllAK5.getHT(), evtwt);
			h1.GetTH1("ABCDval_HT")->Fill( HTAllAK5.getHT(), evtwt);
			h1.GetTH1("ABCDval_NumCA8_A")->Fill(nAv);
			sumw2_av += evtwt2;
			if( selectedAK5Jets.size() == 0 ){
				h1.GetTH1("ABCDval_CutRegion_0ak5")->Fill("A", evtwt); 
				h1.GetTH1("ABCDval_0ak5_HT_A")->Fill( HTAllAK5.getHT(), evtwt);
			}else if( selectedAK5Jets.size() == 1 ){
				h1.GetTH1("ABCDval_CutRegion_1ak5")->Fill("A", evtwt); 
				h1.GetTH1("ABCDval_1ak5_HT_A")->Fill( HTAllAK5.getHT(), evtwt);
			}else if( selectedAK5Jets.size()>= 2 ){
				h1.GetTH1("ABCDval_CutRegion_2ak5")->Fill("A", evtwt); 
				h1.GetTH1("ABCDval_2ak5_HT_A")->Fill( HTAllAK5.getHT(), evtwt);
			}
		}else if( nCv > 0 ){
			h1.GetTH1("ABCDval_CutRegion")->Fill("C", evtwt);
			h1.GetTH1("ABCDval_HT_C")->Fill( HTAllAK5.getHT(), evtwt);
			h1.GetTH1("ABCDval_HT")->Fill( HTAllAK5.getHT(), evtwt);
			h1.GetTH1("ABCDval_NumCA8_C")->Fill(nCv);
			sumw2_cv += evtwt2;
			if( selectedAK5Jets.size() == 0 ){
				h1.GetTH1("ABCDval_CutRegion_0ak5")->Fill("C", evtwt); 
				h1.GetTH1("ABCDval_0ak5_HT_C")->Fill( HTAllAK5.getHT(), evtwt);
			}else if( selectedAK5Jets.size() == 1 ){
				h1.GetTH1("ABCDval_CutRegion_1ak5")->Fill("C", evtwt); 
				h1.GetTH1("ABCDval_1ak5_HT_C")->Fill( HTAllAK5.getHT(), evtwt);
			}else if( selectedAK5Jets.size()>= 2 ){
				h1.GetTH1("ABCDval_CutRegion_2ak5")->Fill("C", evtwt); 
				h1.GetTH1("ABCDval_2ak5_HT_C")->Fill( HTAllAK5.getHT(), evtwt);
			}
		}

		if( nB  != 0 ) h1.GetTH1("ABCDana_CutFlow")->Fill(double(4), evtwt);
		if( nBv != 0 ) h1.GetTH1("ABCDval_CutFlow")->Fill(double(4), evtwt);

		h1.GetTH1("ABCDana_Sumw2_A")->Fill( 0., sumw2_a);
		h1.GetTH1("ABCDana_Sumw2_B")->Fill( 0., sumw2_b);
		h1.GetTH1("ABCDana_Sumw2_C")->Fill( 0., sumw2_c);
		h1.GetTH1("ABCDana_Sumw2_D")->Fill( 0., sumw2_d);
		h1.GetTH1("ABCDana_Sumw2_1b_A")->Fill( 0., sumw2_a_1b);
		h1.GetTH1("ABCDana_Sumw2_1b_B")->Fill( 0., sumw2_b_1b);
		h1.GetTH1("ABCDana_Sumw2_1b_C")->Fill( 0., sumw2_c_1b);
		h1.GetTH1("ABCDana_Sumw2_1b_D")->Fill( 0., sumw2_d_1b);
		h1.GetTH1("ABCDana_Sumw2_2b_A")->Fill( 0., sumw2_a_2b);
		h1.GetTH1("ABCDana_Sumw2_2b_B")->Fill( 0., sumw2_b_2b);
		h1.GetTH1("ABCDana_Sumw2_2b_C")->Fill( 0., sumw2_c_2b);
		h1.GetTH1("ABCDana_Sumw2_2b_D")->Fill( 0., sumw2_d_2b);
		h1.GetTH1("ABCDval_Sumw2_A")->Fill( 0., sumw2_av);
		h1.GetTH1("ABCDval_Sumw2_B")->Fill( 0., sumw2_bv);
		h1.GetTH1("ABCDval_Sumw2_C")->Fill( 0., sumw2_cv);
		h1.GetTH1("ABCDval_Sumw2_D")->Fill( 0., sumw2_dv);
		//// Store new tree, new branch with Jet correction  
		if( BuildMiniTree_ ){
			if( nA+nB+nC+nD > 0 ){ 	
				McFlagana = EvtInfo.McFlag;
				PUana = puwt;
				evtWtana = evtwt;
				HTak5 = HTAllAK5.getHT();
				HThiggsbjet = MyHT.getHT();
				reRegistJet(HiggsJets_ABCD, HiggsJetInfoAna);	
				reRegistJet(AntiHiggsJets_ABCD, AntiHiggsJetInfoAna);	
				reRegistJet(HiggsSubJet1_ABCD, HiggsSubJet1InfoAna);	
				reRegistJet(HiggsSubJet2_ABCD, HiggsSubJet2InfoAna);	
				reRegistJet(AntiHiggsSubJet1_ABCD, AntiHiggsSubJet1InfoAna);	
				reRegistJet(AntiHiggsSubJet2_ABCD, AntiHiggsSubJet2InfoAna);	
				reRegistJet(selectedBJets, bJetInfoAna);	
				reRegistJet(Final_bJets_ABCD, Final_bJetInfoAna);			
				newtree_ana->Fill();	
			}
		}
	} //// entry loop 
}

// ------------ method called once each job just after ending the event loop  ------------
void BackgroundEstimationABCD::endJob(){
	edm::LogInfo("BackgroundEstimationABCD::endJob")<<"evtPass_ana "<<evtPass_ana<<"/"<<maxEvents_;	
	edm::LogInfo("BackgroundEstimationABCD::endJob")<<"evtPass_val "<<evtPass_val<<"/"<<maxEvents_;	
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void BackgroundEstimationABCD::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BackgroundEstimationABCD);
