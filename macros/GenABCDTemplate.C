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
TH1D* drawExtrap( TH1D* hcut, int numZoon, int denZoon, TH1D* hExt );

void GenABCDTemplate(){
	string loadpath = "../test/OnLxplus/";
	string savepath = "results/";

	string meanworkspace="26May_HT950_Mean"; 	
	string meanfileName="Final_histograms_ABCD.root";

	vector<string> SFworkspace, SFfileName, SFTitle;	
	SFworkspace.push_back("26May_HT950_JESUp"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("JESUp");
	SFworkspace.push_back("26May_HT950_JESDown"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("JESDown");
	SFworkspace.push_back("26May_HT950_JERUp"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("JERUp");
	SFworkspace.push_back("26May_HT950_JERDown"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("JERDown");
	SFworkspace.push_back("26May_HT950_SFbUp"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFbUp");
	SFworkspace.push_back("26May_HT950_SFbDown"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFbDown");
	SFworkspace.push_back("26May_HT950_SFlUp"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFlUp");
	SFworkspace.push_back("26May_HT950_SFlDown"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFlDown");
	SFworkspace.push_back("26May_HT950_SFHbUp"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFHbUp");
	SFworkspace.push_back("26May_HT950_SFHbDown"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFHbDown");
	SFworkspace.push_back("26May_HT950_SFHlUp"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFHlUp");
	SFworkspace.push_back("26May_HT950_SFHlDown"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFHlDown");
	const int SFSize=SFworkspace.size();

	string cate[_cat]; cate[_nomal]=""; cate[_1b]="_1b"; cate[_2b]="_2b"; 

	vector<string> hName; vector<int> rebin;
	hName.push_back("HT");		rebin.push_back(100);
	//hName.push_back("HT");	rebin.push_back(50);
	//hName.push_back("bpMass");	rebin.push_back(10);
	//hName.push_back("bpPt");	rebin.push_back(20);

	string dataTitle = "data"; 	string datafullName = "DATA"; 

	vector<string> sigTitle, 	sigfullName; 			
	sigTitle.push_back("BHBH450"); 	sigfullName.push_back("BprimeBprimeToBHBHinc_M-450");  
	sigTitle.push_back("BHBH500"); 	sigfullName.push_back("BprimeBprimeToBHBHinc_M-500"); 
	sigTitle.push_back("BHBH550"); 	sigfullName.push_back("BprimeBprimeToBHBHinc_M-550"); 
	sigTitle.push_back("BHBH600"); 	sigfullName.push_back("BprimeBprimeToBHBHinc_M-600"); 
	sigTitle.push_back("BHBH650");	sigfullName.push_back("BprimeBprimeToBHBHinc_M-650"); 
	sigTitle.push_back("BHBH700"); 	sigfullName.push_back("BprimeBprimeToBHBHinc_M-700"); 
	sigTitle.push_back("BHBH750"); 	sigfullName.push_back("BprimeBprimeToBHBHinc_M-750"); 
	sigTitle.push_back("BHBH800"); 	sigfullName.push_back("BprimeBprimeToBHBHinc_M-800"); 
	sigTitle.push_back("BHBH1000");	sigfullName.push_back("BprimeBprimeToBHBHinc_M-1000"); 
	sigTitle.push_back("BHBH1200");	sigfullName.push_back("BprimeBprimeToBHBHinc_M-1200"); 
	sigTitle.push_back("BHBH1500");	sigfullName.push_back("BprimeBprimeToBHBHinc_M-1500");
	const int sigSize= sigTitle.size();

	// Load files
	TFile* fileMean = new TFile((loadpath+"/"+meanworkspace+"/"+meanfileName).c_str());
	TFile* fileSF[SFSize];
	for( int i=0; i<SFSize; i++ ){
		string filePath = loadpath+"/"+SFworkspace[i]+"/"+SFfileName[i];
		cout<<"Load "<<filePath<<endl;
		fileSF[i] = new TFile(filePath.c_str());
	}
 

	const int hNameSize = hName.size();
	for( int h=0; h<hNameSize; h++){
		cout<<"Make template for "<<hName[h]<<"... "<<endl;
		TFile* fs = new TFile((savepath+"ABCDResultTemplate_ABCDana_"+hName[h]+".root").c_str(), "RECREATE");

		TH1D* hDataCut[_cat];
		TH1D* hData[ABCD][_cat];
		TH1D* hSig[ABCD][_cat][sigSize];
		TH1D* hSigSF[ABCD][_cat][sigSize][SFSize];
		
		TH1D* outhData[_cat][_type];	
		TH1D* outhSig[_cat][sigSize];
		TH1D* outhSigSF[_cat][sigSize][SFSize];

		// Load plots, rebin and over flow 		
		for( int c=0; c<_cat; c++){
			hDataCut[c] = (TH1D*)fileMean->Get((datafullName+"__ABCDana_CutRegion"+cate[c]).c_str());
			hData[A][c] = (TH1D*)fileMean->Get((datafullName+"__ABCDana"+cate[c]+"_"+hName[h]+"_A").c_str()); hData[A][c]->Rebin(rebin[h]);
			hData[B][c] = (TH1D*)fileMean->Get((datafullName+"__ABCDana"+cate[c]+"_"+hName[h]+"_B").c_str()); hData[B][c]->Rebin(rebin[h]);
			hData[C][c] = (TH1D*)fileMean->Get((datafullName+"__ABCDana"+cate[c]+"_"+hName[h]+"_C").c_str()); hData[C][c]->Rebin(rebin[h]);
			hData[D][c] = (TH1D*)fileMean->Get((datafullName+"__ABCDana"+cate[c]+"_"+hName[h]+"_D").c_str()); hData[D][c]->Rebin(rebin[h]);
			fix(hData[A][c]);
			fix(hData[B][c]);
			fix(hData[C][c]);
			fix(hData[D][c]);
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
		cout<<"	Data..."<<endl;
		outhData[_nomal][_obs] = (TH1D*)hData[observeRegion][_nomal]->Clone("CatAll_data_obs");
		outhData[_1b][_obs] = (TH1D*)hData[observeRegion][_1b]->Clone("Cat1b_data_obs");
		outhData[_2b][_obs] = (TH1D*)hData[observeRegion][_2b]->Clone("Cat2b_data_obs");

		outhData[_nomal][_exp] = (TH1D*)drawExtrap( hDataCut[_nomal], D, C, hData[A][_nomal])->Clone("CatAll_background");
		GetSigma(outhData[_nomal][_exp],"CatAll_background_SystUp", 		 1);
		GetSigma(outhData[_nomal][_exp],"CatAll_background_SystDown", 		-1);

		outhData[_1b][_exp] = (TH1D*)drawExtrap( hDataCut[_1b], D, C, hData[A][_1b])->Clone("Cat1b_background");
		GetSigma(outhData[_1b][_exp], 	"Cat1b_background_SystUp", 	 1);
		GetSigma(outhData[_1b][_exp], 	"Cat1b_background_SystDown", 	-1);

		outhData[_2b][_exp] = (TH1D*)drawExtrap( hDataCut[_2b], D, C, hData[A][_2b])->Clone("Cat2b_background");
		GetSigma(outhData[_2b][_exp], 	"Cat2b_background_SystUp", 	 1);
		GetSigma(outhData[_2b][_exp], 	"Cat2b_background_SystDown",	-1);

		cout<<"	Signal..."<<endl;
		for( int s=0; s<sigSize; s++ ){
			// Combined
			outhSig[_nomal][s] = (TH1D*)hSig[observeRegion][_nomal][s]->Clone(("CatAll_"+sigTitle[s]).c_str());
			for( int sf=0; sf<SFSize; sf++){
				outhSigSF[_nomal][s][sf] = (TH1D*)hSigSF[observeRegion][_nomal][s][sf]->Clone(("CatAll_"+sigTitle[s]+"_"+SFTitle[sf]).c_str());
			}
 			// 1b
			outhSig[_1b][s] = (TH1D*)hSig[observeRegion][_1b][s]->Clone(("Cat1b_"+sigTitle[s]).c_str());
			for( int sf=0; sf<SFSize; sf++){
				outhSigSF[_1b][s][sf] = (TH1D*)hSigSF[observeRegion][_1b][s][sf]->Clone(("Cat1b_"+sigTitle[s]+"_"+SFTitle[sf]).c_str()); 
			}
			//2b
			outhSig[_2b][s] = (TH1D*)hSig[observeRegion][_2b][s]->Clone(("Cat2b_"+sigTitle[s]).c_str());
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
TH1D* drawExtrap( TH1D* hcut, int numZoon, int denZoon, TH1D* hExt ){

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
	
	return hExt;
}
