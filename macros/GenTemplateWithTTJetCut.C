#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TText.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TPaveText.h"
#include "TAxis.h"
#include "help.C"
using namespace std;
const int A = 0;
const int B = 1;
const int C = 2;
const int D = 3;
const int ABCD = 4;

const int _nomal =0;
const int _1b =1;
const int _2b =2;
const int _cat = 3;

const int _obs = 0;
const int _exp = 1;
const int _type = 2;
const int observeRegion = B;

void GetSigma( TH1D* h, string name, int shift=0 );
void SubtractTTJets( TH1D* h, TH1D* h_tt, string name );
void drawExtrapCutTTJets( TH1D* hcut, int numZoon, int denZoon, TH1D* hExt_, TH1D* hTTJetsSB, TH1D* hTTJets, string name, bool unc=0, bool getQCD=0 );
void drawExtrap( TH1D* hcut, int numZoon, int denZoon, TH1D* hExt, string name, bool unc=0 );

void GenTemplateWithTTJetCut(){
	//string loadpath = "../test/OnLxplus/";
	//string loadpath = "../test/OnLxplus/SampleAndData";
	string loadpath = "../test/OnLxplus/SampleAndData/";
	string savepath = "results/ForScanBprimeBR_17Sep/";

	//string meanworkspace="08July_Mean"; 	
	string meanworkspace="Mean"; 	
	string meanfileName="Final_histograms_ABCD.root";

	vector<string> SFworkspace, SFfileName, SFTitle;	
	vector<string> TTAddSFTitle,TTAddSFfullName;
	//Approval results used 16July_	
	SFworkspace.push_back("JESUp"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("JESUp");
	SFworkspace.push_back("JESDown");	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("JESDown");
	SFworkspace.push_back("JERUp"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("JERUp");
	SFworkspace.push_back("JERDown");	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("JERDown");
	SFworkspace.push_back("SFbUp"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFbUp");
	SFworkspace.push_back("SFbDown");	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFbDown");
	SFworkspace.push_back("SFlUp"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFlUp");
	SFworkspace.push_back("SFlDown"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFlDown");
	//SFworkspace.push_back("16July_SFHbUp"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFHbUp");
	//SFworkspace.push_back("16July_SFHbDown");	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFHbDown");
	//SFworkspace.push_back("16July_SFHlUp"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFHlUp");
	//SFworkspace.push_back("16July_SFHlDown");	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFHlDown");
	SFworkspace.push_back("CA8Up"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("CA8Up");
	SFworkspace.push_back("CA8Down");	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("CA8Down");
	SFworkspace.push_back("PUUp"); 		SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("PUUp");
	SFworkspace.push_back("PUDown");	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("PUDown");

	//string TTSFworkspace="16July_TTJets_Unc"; 	
	string TTSFworkspace="TTJetsUnc"; 	
	string TTSFfileName="Final_histograms_ABCD.root"; 
	TTAddSFTitle.push_back("TopPtReWrtUp");  	TTAddSFfullName.push_back("TTJets_TopPtReWrtUp");
	TTAddSFTitle.push_back("TopPtReWrtDown");  	TTAddSFfullName.push_back("TTJets_TopPtReWrtDown");
	TTAddSFTitle.push_back("TTJetsMatchingUp");  	TTAddSFfullName.push_back("TTJets_MatchingUp");
	TTAddSFTitle.push_back("TTJetsMatchingDown");	TTAddSFfullName.push_back("TTJets_MatchingDown");
	TTAddSFTitle.push_back("TTJetsScaleUp");	TTAddSFfullName.push_back("TTJets_ScaleUp");
	TTAddSFTitle.push_back("TTJetsScaleDown");	TTAddSFfullName.push_back("TTJets_ScaleDown");
	const int SFSize=SFworkspace.size();
	const int TTAddSFSize=TTAddSFTitle.size();
	const int TTSFSize=SFSize+TTAddSFSize;

	string cate[_cat]; cate[_nomal]=""; cate[_1b]="_1b"; cate[_2b]="_2b"; 

	vector<string> hName; vector<int> rebin;
	hName.push_back("HT");		rebin.push_back(300);

	string dataTitle = "data"; 	string datafullName = "DATA"; 
	string TTJetsTitle = "TTJets"; 	string TTJetsfullName = "TTJets"; 

	vector<string> sigTitle, 	sigfullName; 			
	sigTitle.push_back("BHBH500"); 	sigfullName.push_back("BprimeBprimeToBHBHinc_M-500"); 
	sigTitle.push_back("BHBH600"); 	sigfullName.push_back("BprimeBprimeToBHBHinc_M-600"); 
	sigTitle.push_back("BHBH700"); 	sigfullName.push_back("BprimeBprimeToBHBHinc_M-700"); 
	sigTitle.push_back("BHBH800"); 	sigfullName.push_back("BprimeBprimeToBHBHinc_M-800"); 
	sigTitle.push_back("BHBH900"); 	sigfullName.push_back("BprimeBprimeToBHBHinc_M-900"); 
	sigTitle.push_back("BHBH1000");	sigfullName.push_back("BprimeBprimeToBHBHinc_M-1000"); 
	sigTitle.push_back("BHBH1200");	sigfullName.push_back("BprimeBprimeToBHBHinc_M-1200"); 
	sigTitle.push_back("BHBZ500"); 	sigfullName.push_back("BprimeBprimeToBHBZinc_M-500"); 
	sigTitle.push_back("BHBZ600"); 	sigfullName.push_back("BprimeBprimeToBHBZinc_M-600"); 
	sigTitle.push_back("BHBZ700"); 	sigfullName.push_back("BprimeBprimeToBHBZinc_M-700"); 
	sigTitle.push_back("BHBZ800"); 	sigfullName.push_back("BprimeBprimeToBHBZinc_M-800"); 
	sigTitle.push_back("BHBZ900"); 	sigfullName.push_back("BprimeBprimeToBHBZinc_M-900"); 
	sigTitle.push_back("BHBZ1000");	sigfullName.push_back("BprimeBprimeToBHBZinc_M-1000"); 
	sigTitle.push_back("BHBZ1200");	sigfullName.push_back("BprimeBprimeToBHBZinc_M-1200"); 
	sigTitle.push_back("BHTW500"); 	sigfullName.push_back("BprimeBprimeToBHTWinc_M-500"); 
	sigTitle.push_back("BHTW600"); 	sigfullName.push_back("BprimeBprimeToBHTWinc_M-600"); 
	sigTitle.push_back("BHTW700"); 	sigfullName.push_back("BprimeBprimeToBHTWinc_M-700"); 
	sigTitle.push_back("BHTW800"); 	sigfullName.push_back("BprimeBprimeToBHTWinc_M-800"); 
	sigTitle.push_back("BHTW900"); 	sigfullName.push_back("BprimeBprimeToBHTWinc_M-900"); 
	sigTitle.push_back("BHTW1000");	sigfullName.push_back("BprimeBprimeToBHTWinc_M-1000"); 
	sigTitle.push_back("BHTW1200");	sigfullName.push_back("BprimeBprimeToBHTWinc_M-1200"); 
	sigTitle.push_back("TWTW500"); 	sigfullName.push_back("BprimeBprimeToTWTWinc_M-500"); 
	sigTitle.push_back("TWTW600"); 	sigfullName.push_back("BprimeBprimeToTWTWinc_M-600"); 
	sigTitle.push_back("TWTW700"); 	sigfullName.push_back("BprimeBprimeToTWTWinc_M-700"); 
	sigTitle.push_back("TWTW800"); 	sigfullName.push_back("BprimeBprimeToTWTWinc_M-800"); 
	sigTitle.push_back("TWTW900"); 	sigfullName.push_back("BprimeBprimeToTWTWinc_M-900"); 
	sigTitle.push_back("TWTW1000");	sigfullName.push_back("BprimeBprimeToTWTWinc_M-1000"); 
	sigTitle.push_back("TWTW1200");	sigfullName.push_back("BprimeBprimeToTWTWinc_M-1200"); 
	sigTitle.push_back("BZTW500"); 	sigfullName.push_back("BprimeBprimeToBZTWinc_M-500"); 
	sigTitle.push_back("BZTW600"); 	sigfullName.push_back("BprimeBprimeToBZTWinc_M-600"); 
	sigTitle.push_back("BZTW700"); 	sigfullName.push_back("BprimeBprimeToBZTWinc_M-700"); 
	sigTitle.push_back("BZTW800"); 	sigfullName.push_back("BprimeBprimeToBZTWinc_M-800"); 
	sigTitle.push_back("BZTW900"); 	sigfullName.push_back("BprimeBprimeToBZTWinc_M-900"); 
	sigTitle.push_back("BZTW1000");	sigfullName.push_back("BprimeBprimeToBZTWinc_M-1000"); 
	sigTitle.push_back("BZTW1200");	sigfullName.push_back("BprimeBprimeToBZTWinc_M-1200"); 
	sigTitle.push_back("BZBZ500"); 	sigfullName.push_back("BprimeBprimeToBZBZinc_M-500"); 
	sigTitle.push_back("BZBZ600"); 	sigfullName.push_back("BprimeBprimeToBZBZinc_M-600"); 
	sigTitle.push_back("BZBZ700"); 	sigfullName.push_back("BprimeBprimeToBZBZinc_M-700"); 
	sigTitle.push_back("BZBZ800"); 	sigfullName.push_back("BprimeBprimeToBZBZinc_M-800"); 
	sigTitle.push_back("BZBZ900"); 	sigfullName.push_back("BprimeBprimeToBZBZinc_M-900"); 
	sigTitle.push_back("BZBZ1000");	sigfullName.push_back("BprimeBprimeToBZBZinc_M-1000"); 
	sigTitle.push_back("BZBZ1200");	sigfullName.push_back("BprimeBprimeToBZBZinc_M-1200"); 
	const int sigSize= sigTitle.size();

	// Load files
	TFile* fileMean = new TFile((loadpath+"/"+meanworkspace+"/"+meanfileName).c_str());
	TFile* fileTTAddSF = new TFile((loadpath+"/"+TTSFworkspace+"/"+TTSFfileName).c_str());
	TFile* fileSF[SFSize];
	for( int i=0; i<SFSize; i++ ){
		string filePath = loadpath+"/"+SFworkspace[i]+"/"+SFfileName[i];
		cout<<"Load "<<filePath<<endl;
		fileSF[i] = new TFile(filePath.c_str());
	}

	const int hNameSize = hName.size();
	for( int h=0; h<hNameSize; h++){
		cout<<"Make template for "<<hName[h]<<"... "<<endl;
		TFile* fs = new TFile((savepath+"ABCDResultTemplate_CutTTJet2_"+hName[h]+".root").c_str(), "RECREATE");

		TH1D* hDataCut[_cat];
		TH1D* hData[ABCD][_cat];
		TH1D* hTTJets[ABCD][_cat];
		TH1D* hTTJetsSF[ABCD][_cat][TTSFSize];
		TH1D* hSig[ABCD][_cat][sigSize];
		TH1D* hSigSF[ABCD][_cat][sigSize][SFSize];
		
		TH1D* outhData[_cat][_type];	
		TH1D* outhTTJets[_cat];	
		TH1D* outhTTJetsSF[_cat][TTSFSize];	
		TH1D* outhSig[_cat][sigSize];
		TH1D* outhSigSF[_cat][sigSize][SFSize];

		vector<string> TTSFTitle;
		// Load plots, rebin and over flow 		
		for( int c=0; c<_cat; c++){
			//Data
			hDataCut[c] = (TH1D*)fileMean->Get((datafullName+"__ABCDana_CutRegion"+cate[c]).c_str());
			hData[A][c] = (TH1D*)fileMean->Get((datafullName+"__ABCDana"+cate[c]+"_"+hName[h]+"_A").c_str()); hData[A][c]->Rebin(rebin[h]);
			hData[B][c] = (TH1D*)fileMean->Get((datafullName+"__ABCDana"+cate[c]+"_"+hName[h]+"_B").c_str()); hData[B][c]->Rebin(rebin[h]);
			hData[C][c] = (TH1D*)fileMean->Get((datafullName+"__ABCDana"+cate[c]+"_"+hName[h]+"_C").c_str()); hData[C][c]->Rebin(rebin[h]);
			hData[D][c] = (TH1D*)fileMean->Get((datafullName+"__ABCDana"+cate[c]+"_"+hName[h]+"_D").c_str()); hData[D][c]->Rebin(rebin[h]);
			fix(hData[A][c]);
			fix(hData[B][c]);
			fix(hData[C][c]);
			fix(hData[D][c]);
			
			//TTJets
			hTTJets[A][c] = (TH1D*)fileMean->Get((TTJetsfullName+"__ABCDana"+cate[c]+"_"+hName[h]+"_A").c_str()); hTTJets[A][c]->Rebin(rebin[h]);
			hTTJets[B][c] = (TH1D*)fileMean->Get((TTJetsfullName+"__ABCDana"+cate[c]+"_"+hName[h]+"_B").c_str()); hTTJets[B][c]->Rebin(rebin[h]);
			hTTJets[C][c] = (TH1D*)fileMean->Get((TTJetsfullName+"__ABCDana"+cate[c]+"_"+hName[h]+"_C").c_str()); hTTJets[C][c]->Rebin(rebin[h]);
			hTTJets[D][c] = (TH1D*)fileMean->Get((TTJetsfullName+"__ABCDana"+cate[c]+"_"+hName[h]+"_D").c_str()); hTTJets[D][c]->Rebin(rebin[h]);
			fix(hTTJets[A][c]);
			fix(hTTJets[B][c]);
			fix(hTTJets[C][c]);
			fix(hTTJets[D][c]);
			for( int sf=0; sf<SFSize; sf++ ){
					hTTJetsSF[A][c][sf]  = (TH1D*)fileSF[sf]->Get((TTJetsfullName+"__ABCDana"+cate[c]+"_"+hName[h]+"_A").c_str()); hTTJetsSF[A][c][sf]->Rebin(rebin[h]);
					hTTJetsSF[B][c][sf]  = (TH1D*)fileSF[sf]->Get((TTJetsfullName+"__ABCDana"+cate[c]+"_"+hName[h]+"_B").c_str()); hTTJetsSF[B][c][sf]->Rebin(rebin[h]);
					hTTJetsSF[C][c][sf]  = (TH1D*)fileSF[sf]->Get((TTJetsfullName+"__ABCDana"+cate[c]+"_"+hName[h]+"_C").c_str()); hTTJetsSF[C][c][sf]->Rebin(rebin[h]);
					hTTJetsSF[D][c][sf]  = (TH1D*)fileSF[sf]->Get((TTJetsfullName+"__ABCDana"+cate[c]+"_"+hName[h]+"_D").c_str()); hTTJetsSF[D][c][sf]->Rebin(rebin[h]);
					fix(hTTJetsSF[A][c][sf]);
					fix(hTTJetsSF[B][c][sf]);
					fix(hTTJetsSF[C][c][sf]);
					fix(hTTJetsSF[D][c][sf]);
					TTSFTitle.push_back(SFTitle[sf]);
			}
			for( int sfa=0; sfa<TTAddSFSize; sfa++ ){
					int sf=sfa+SFSize;
					hTTJetsSF[A][c][sf]  = (TH1D*)fileTTAddSF->Get((TTAddSFfullName[sfa]+"__ABCDana"+cate[c]+"_"+hName[h]+"_A").c_str()); hTTJetsSF[A][c][sf]->Rebin(rebin[h]);
					hTTJetsSF[B][c][sf]  = (TH1D*)fileTTAddSF->Get((TTAddSFfullName[sfa]+"__ABCDana"+cate[c]+"_"+hName[h]+"_B").c_str()); hTTJetsSF[B][c][sf]->Rebin(rebin[h]);
					hTTJetsSF[C][c][sf]  = (TH1D*)fileTTAddSF->Get((TTAddSFfullName[sfa]+"__ABCDana"+cate[c]+"_"+hName[h]+"_C").c_str()); hTTJetsSF[C][c][sf]->Rebin(rebin[h]);
					hTTJetsSF[D][c][sf]  = (TH1D*)fileTTAddSF->Get((TTAddSFfullName[sfa]+"__ABCDana"+cate[c]+"_"+hName[h]+"_D").c_str()); hTTJetsSF[D][c][sf]->Rebin(rebin[h]);
					fix(hTTJetsSF[A][c][sf]);
					fix(hTTJetsSF[B][c][sf]);
					fix(hTTJetsSF[C][c][sf]);
					fix(hTTJetsSF[D][c][sf]);
					TTSFTitle.push_back(TTAddSFTitle[sfa]);
			}

			//Signal MC
			for( int s=0; s<sigSize; s++ ){
				hSig[A][c][s]  = (TH1D*)fileMean->Get((sigfullName[s]+"__ABCDana"+cate[c]+"_"+hName[h]+"_A").c_str()); hSig[A][c][s]->Rebin(rebin[h]);
				hSig[B][c][s]  = (TH1D*)fileMean->Get((sigfullName[s]+"__ABCDana"+cate[c]+"_"+hName[h]+"_B").c_str()); hSig[B][c][s]->Rebin(rebin[h]);
				hSig[C][c][s]  = (TH1D*)fileMean->Get((sigfullName[s]+"__ABCDana"+cate[c]+"_"+hName[h]+"_C").c_str()); hSig[C][c][s]->Rebin(rebin[h]);
				hSig[D][c][s]  = (TH1D*)fileMean->Get((sigfullName[s]+"__ABCDana"+cate[c]+"_"+hName[h]+"_D").c_str()); hSig[D][c][s]->Rebin(rebin[h]);
				fix(hSig[A][c][s]);
				fix(hSig[B][c][s]);
				fix(hSig[C][c][s]);
				fix(hSig[D][c][s]);
				for( int sf=0; sf<SFSize; sf++ ){
					hSigSF[A][c][s][sf]  = (TH1D*)fileSF[sf]->Get((sigfullName[s]+"__ABCDana"+cate[c]+"_"+hName[h]+"_A").c_str()); hSigSF[A][c][s][sf]->Rebin(rebin[h]);
					hSigSF[B][c][s][sf]  = (TH1D*)fileSF[sf]->Get((sigfullName[s]+"__ABCDana"+cate[c]+"_"+hName[h]+"_B").c_str()); hSigSF[B][c][s][sf]->Rebin(rebin[h]);
					hSigSF[C][c][s][sf]  = (TH1D*)fileSF[sf]->Get((sigfullName[s]+"__ABCDana"+cate[c]+"_"+hName[h]+"_C").c_str()); hSigSF[C][c][s][sf]->Rebin(rebin[h]);
					hSigSF[D][c][s][sf]  = (TH1D*)fileSF[sf]->Get((sigfullName[s]+"__ABCDana"+cate[c]+"_"+hName[h]+"_D").c_str()); hSigSF[D][c][s][sf]->Rebin(rebin[h]);
					fix(hSigSF[A][c][s][sf]);
					fix(hSigSF[B][c][s][sf]);
					fix(hSigSF[C][c][s][sf]);
					fix(hSigSF[D][c][s][sf]);
				}
			}
		}

		// Store plot to new file
		//// data mc
		cout<<"	Data..."<<endl;
		outhData[_nomal][_obs] = (TH1D*)hData[observeRegion][_nomal]->Clone("CatAll_data_obs");
		outhData[_1b][_obs] = (TH1D*)hData[observeRegion][_1b]->Clone("Cat1b_data_obs");
		outhData[_2b][_obs] = (TH1D*)hData[observeRegion][_2b]->Clone("Cat2b_data_obs");

		SubtractTTJets( hData[B][_nomal], hTTJets[observeRegion][_nomal], "CatAll_data_SubtractTTJets" );
		SubtractTTJets( hData[B][_1b],    hTTJets[observeRegion][_1b], "Cat1b_data_SubtractTTJets" );
		SubtractTTJets( hData[B][_2b],    hTTJets[observeRegion][_2b], "Cat2b_data_SubtractTTJets" );

		drawExtrapCutTTJets( hDataCut[_nomal], D, C, hData[A][_nomal], hTTJets[A][_nomal], hTTJets[observeRegion][_nomal], "CatAll_background", true);
		for( int sf=0; sf<TTSFSize; sf++){
			drawExtrapCutTTJets( hDataCut[_nomal], D, C, hData[A][_nomal], hTTJetsSF[A][_nomal][sf], hTTJetsSF[observeRegion][_nomal][sf], "CatAll_background_"+TTSFTitle[sf]);
		}
		drawExtrapCutTTJets( hDataCut[_1b],    D, C, hData[A][_1b],    hTTJets[A][_1b],    hTTJets[observeRegion][_1b],    "Cat1b_background", true);
		for( int sf=0; sf<TTSFSize; sf++){
			drawExtrapCutTTJets( hDataCut[_1b], D, C, hData[A][_1b], hTTJetsSF[A][_1b][sf], hTTJetsSF[observeRegion][_1b][sf], "Cat1b_background_"+TTSFTitle[sf]);
		}
		drawExtrapCutTTJets( hDataCut[_2b],    D, C, hData[A][_2b],    hTTJets[A][_2b],    hTTJets[observeRegion][_2b],    "Cat2b_background", true);
		for( int sf=0; sf<TTSFSize; sf++){
			drawExtrapCutTTJets( hDataCut[_2b], D, C, hData[A][_2b], hTTJetsSF[A][_2b][sf], hTTJetsSF[observeRegion][_2b][sf], "Cat2b_background_"+TTSFTitle[sf]);
		}

		drawExtrap( hDataCut[_nomal], D, C, hData[A][_nomal], "CatAll_background_fulldata");
		drawExtrap( hDataCut[_1b], D, C, hData[A][_1b], "Cat1b_background_fulldata");
		drawExtrap( hDataCut[_2b], D, C, hData[A][_2b], "Cat2b_background_fulldata");

		//// TTJets mc
		cout<<"	TTJets..."<<endl;
		outhTTJets[_nomal] = (TH1D*)hTTJets[observeRegion][_nomal]->Clone("CatAll_TTJets");
		GetSigma(outhTTJets[_nomal],	"CatAll_TTJets_StatUp",  	1);
		GetSigma(outhTTJets[_nomal],	"CatAll_TTJets_StatDown", 	 -1);
		for( int sf=0; sf<TTSFSize; sf++){
			outhTTJetsSF[_nomal][sf] = (TH1D*)hTTJetsSF[observeRegion][_nomal][sf]->Clone(("CatAll_TTJets_"+TTSFTitle[sf]).c_str());
		}
		outhTTJets[_1b] = (TH1D*)hTTJets[observeRegion][_1b]->Clone("Cat1b_TTJets");
		GetSigma(outhTTJets[_1b],	"Cat1b_TTJets_StatUp",  	1);
		GetSigma(outhTTJets[_1b],	"Cat1b_TTJets_StatDown", 	 -1);
		for( int sf=0; sf<TTSFSize; sf++){
			outhTTJetsSF[_1b][sf] = (TH1D*)hTTJetsSF[observeRegion][_1b][sf]->Clone(("Cat1b_TTJets_"+TTSFTitle[sf]).c_str());
		}
		outhTTJets[_2b] = (TH1D*)hTTJets[observeRegion][_2b]->Clone("Cat2b_TTJets");
		GetSigma(outhTTJets[_2b],	"Cat2b_TTJets_StatUp",  	1);
		GetSigma(outhTTJets[_2b],	"Cat2b_TTJets_StatDown", 	 -1);
		for( int sf=0; sf<TTSFSize; sf++){
			outhTTJetsSF[_2b][sf] = (TH1D*)hTTJetsSF[observeRegion][_2b][sf]->Clone(("Cat2b_TTJets_"+TTSFTitle[sf]).c_str());
		}

		drawExtrapCutTTJets( hDataCut[_nomal], D, C, hData[A][_nomal], hTTJets[A][_nomal], hTTJets[observeRegion][_nomal], "CatAll_QCD_data", false, true);
		drawExtrapCutTTJets( hDataCut[_1b], D, C, hData[A][_1b], hTTJets[A][_1b], hTTJets[observeRegion][_1b], "Cat1b_QCD_data", false, true);
		drawExtrapCutTTJets( hDataCut[_2b], D, C, hData[A][_2b], hTTJets[A][_2b], hTTJets[observeRegion][_2b], "Cat2b_QCD_data", false, true);
		//// signal mc
		cout<<"	Signal..."<<endl;
		for( int s=0; s<sigSize; s++ ){
			// Combined
			outhSig[_nomal][s] = (TH1D*)hSig[observeRegion][_nomal][s]->Clone(("CatAll_"+sigTitle[s]).c_str());
			GetSigma(outhSig[_nomal][s],	"CatAll_"+sigTitle[s]+"_StatUp", 		 1);
			GetSigma(outhSig[_nomal][s],	"CatAll_"+sigTitle[s]+"_StatDown", 		 -1);
			for( int sf=0; sf<SFSize; sf++){
				outhSigSF[_nomal][s][sf] = (TH1D*)hSigSF[observeRegion][_nomal][s][sf]->Clone(("CatAll_"+sigTitle[s]+"_"+SFTitle[sf]).c_str());
			}
 			// 1b
			outhSig[_1b][s] = (TH1D*)hSig[observeRegion][_1b][s]->Clone(("Cat1b_"+sigTitle[s]).c_str());
			GetSigma(outhSig[_1b][s],	"Cat1b_"+sigTitle[s]+"_StatUp", 		 1);
			GetSigma(outhSig[_1b][s],	"Cat1b_"+sigTitle[s]+"_StatDown", 		 -1);
			for( int sf=0; sf<SFSize; sf++){
				outhSigSF[_1b][s][sf] = (TH1D*)hSigSF[observeRegion][_1b][s][sf]->Clone(("Cat1b_"+sigTitle[s]+"_"+SFTitle[sf]).c_str()); 
			}
			//2b
			outhSig[_2b][s] = (TH1D*)hSig[observeRegion][_2b][s]->Clone(("Cat2b_"+sigTitle[s]).c_str());
			GetSigma(outhSig[_2b][s],	"Cat2b_"+sigTitle[s]+"_StatUp", 		 1);
			GetSigma(outhSig[_2b][s],	"Cat2b_"+sigTitle[s]+"_StatDown", 		 -1);
			for( int sf=0; sf<SFSize; sf++){
				outhSigSF[_2b][s][sf] = (TH1D*)hSigSF[observeRegion][_2b][s][sf]->Clone(("Cat2b_"+sigTitle[s]+"_"+SFTitle[sf]).c_str()); 
			}
		}
		fs->Write();
	}	
}
void GetSigma( TH1D* h, string name, int shift ){
	if( abs(shift)<1 || abs(shift)>1 ) return;
	TH1D* h_ = (TH1D*)h->Clone(name.c_str()); 
	int maxBin = h->GetXaxis()->GetLast();;
	int minBin = 1;
	for( int i=minBin; i<=maxBin; i++){
		double con_(0);
		if( shift == 1 ) con_ = h->GetBinContent(i) + h->GetBinError(i);
		if( shift == -1 ) con_ = h->GetBinContent(i) - h->GetBinError(i);
		h_->SetBinContent(i, con_);	
	}
}
void SubtractTTJets( TH1D* h, TH1D* h_tt, string name ){
	TH1D* h_ = (TH1D*)h->Clone(name.c_str());
	h_->Add(h_tt, -1); 

}
void drawExtrap( TH1D* hcut, int numZoon, int denZoon, TH1D* hExt_, string name, bool unc ){
	
	TH1D* hExt = (TH1D*)hExt_->Clone(name.c_str());
	
	double yield[4], error[4]; 
	yield[A] = hcut->GetBinContent(1);
	yield[B] = hcut->GetBinContent(2);
	yield[C] = hcut->GetBinContent(3);
	yield[D] = hcut->GetBinContent(4);
	error[A] = hcut->GetBinError(1);
	error[B] = hcut->GetBinError(2);
	error[C] = hcut->GetBinError(3);
	error[D] = hcut->GetBinError(4);

	double sigmaNum = error[numZoon]/yield[numZoon];
	double sigmaDen = error[denZoon]/yield[denZoon];
	double extRatio = yield[numZoon]/yield[denZoon];

	int maxBin = hExt->GetXaxis()->GetLast();;
	int minBin = 1;

	for( int i=minBin; i<=maxBin; i++){
		double extBinContent = hExt->GetBinContent(i);
		double extEr = hExt->GetBinError(i);
		double sigmaExt = extEr/extBinContent;
		double binEr;
		if( extBinContent == 0 || yield[denZoon] == 0 ){
			binEr = 0;
		}else{
			binEr = sqrt(sigmaNum*sigmaNum+sigmaExt*sigmaExt+sigmaDen*sigmaDen)*extBinContent;
		}
		hExt->SetBinError(i, binEr);	
	}
	hExt->Scale(extRatio);	
	
	if( unc ){
		GetSigma( hExt, name+"_StatUp", 1);
		GetSigma( hExt, name+"_StatDown", -1);
	}
}
void drawExtrapCutTTJets( TH1D* hcut, int numZoon, int denZoon, TH1D* hExt_, TH1D* hTTJetsSB, TH1D* hTTJets, string name, bool unc, bool getQCD ){

	TH1D* hExt = (TH1D*)hExt_->Clone(name.c_str());
	hExt->Add(hTTJetsSB, -1);
	
	double yield[4], error[4]; 
	yield[A] = hcut->GetBinContent(1);
	yield[B] = hcut->GetBinContent(2);
	yield[C] = hcut->GetBinContent(3);
	yield[D] = hcut->GetBinContent(4);
	error[A] = hcut->GetBinError(1);
	error[B] = hcut->GetBinError(2);
	error[C] = hcut->GetBinError(3);
	error[D] = hcut->GetBinError(4);

	double sigmaNum = error[numZoon]/yield[numZoon];
	double sigmaDen = error[denZoon]/yield[denZoon];
	double extRatio = yield[numZoon]/yield[denZoon];

	int maxBin = hExt->GetXaxis()->GetLast();;
	int minBin = 1;

	for( int i=minBin; i<=maxBin; i++){
		double extBinContent = hExt->GetBinContent(i);
		double extEr = hExt->GetBinError(i);
		double sigmaExt = extEr/extBinContent;
		double binEr;
		if( extBinContent == 0 || yield[denZoon] == 0 ){
			binEr = 0;
		}else{
			binEr = sqrt(sigmaNum*sigmaNum+sigmaExt*sigmaExt+sigmaDen*sigmaDen)*extBinContent;
		}
		hExt->SetBinError(i, binEr);	
	}
	hExt->Scale(extRatio);
	if ( !getQCD )	hExt->Add(hTTJets);
	
	if( unc ){
		GetSigma( hExt, name+"_StatUp", 1);
		GetSigma( hExt, name+"_StatDown", -1);
	}
}

