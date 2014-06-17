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
void newbinPlots( TH1D* h1, TH1D* h2, int rebin=1, bool doOverFlow=0, double overFlowValue=1000000); 
TH1D* drawExtrap( TH1D* hcut, int numZoon, int denZoon, TH1D* hExt );

void _GenABCDTemplateBin(){
	string loadpath = "../test/OnLxplus/";
	string savepath = "results/";

	string meanworkspace="07June_HT950_Mean"; 	
	string meanfileName="Final_histograms_ABCD.root";

	vector<string> SFworkspace, SFfileName, SFTitle;	
	SFworkspace.push_back("07June_HT950_JESUp"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("JESUp");
	SFworkspace.push_back("07June_HT950_JESDown"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("JESDown");
	SFworkspace.push_back("07June_HT950_JERUp"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("JERUp");
	SFworkspace.push_back("07June_HT950_JERDown"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("JERDown");
	SFworkspace.push_back("07June_HT950_SFbUp"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFbUp");
	SFworkspace.push_back("07June_HT950_SFbDown"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFbDown");
	SFworkspace.push_back("07June_HT950_SFlUp"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFlUp");
	SFworkspace.push_back("07June_HT950_SFlDown"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFlDown");
	SFworkspace.push_back("07June_HT950_SFHbUp"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFHbUp");
	SFworkspace.push_back("07June_HT950_SFHbDown");	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFHbDown");
	SFworkspace.push_back("07June_HT950_SFHlUp"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFHlUp");
	SFworkspace.push_back("07June_HT950_SFHlDown");	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("SFHlDown");
	SFworkspace.push_back("07June_HT950_PUUp"); 	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("PUUp");
	SFworkspace.push_back("07June_HT950_PUDown");	SFfileName.push_back("Final_histograms_ABCD.root"); SFTitle.push_back("PUDown");
	SFworkspace.push_back("HMassUncert_1sigUp"); 	SFfileName.push_back("Final_histograms_BprimebH.root"); SFTitle.push_back("HMUp");
	SFworkspace.push_back("HMassUncert_1sigDown");	SFfileName.push_back("Final_histograms_BprimebH.root"); SFTitle.push_back("HMDown");
	const int SFSize=SFworkspace.size();

	string cate[_cat]; cate[_nomal]=""; cate[_1b]="_1b"; cate[_2b]="_2b"; 

	vector<string> hName; vector<int> rebin, newbinsize, newbinmin, newbinmax;
	vector<bool> catalldoOverFlow; vector<double> catallOverFlowValue;
	vector<bool> cat1bdoOverFlow; vector<double> cat1bOverFlowValue;
	vector<bool> cat2bdoOverFlow; vector<double> cat2bOverFlowValue;
	//hName.push_back("HT");	rebin.push_back(400); newbinsize.push_back(2800); newbinmin.push_back(-250); newbinmax.push_back(2550);
	//hName.push_back("HT");	rebin.push_back(300); newbinsize.push_back(2700); newbinmin.push_back(-250); newbinmax.push_back(2450);
	hName.push_back("HT");		rebin.push_back(100); newbinsize.push_back(2700); newbinmin.push_back(-250); newbinmax.push_back(2450);
	catalldoOverFlow.push_back(0); catallOverFlowValue.push_back(1000000); 
	cat1bdoOverFlow.push_back(0); cat1bOverFlowValue.push_back(1000000); 
	cat2bdoOverFlow.push_back(0); cat2bOverFlowValue.push_back(1550); //Only for 300 
	
	string dataTitle = "data"; 	string datafullName = "DATA"; 

	vector<string> sigTitle, 	sigfullName; 			
	//sigTitle.push_back("BHBH450"); 	sigfullName.push_back("BprimeBprimeToBHBHinc_M-450");  
	sigTitle.push_back("BHBH500"); 	sigfullName.push_back("BprimeBprimeToBHBHinc_M-500"); 
	//sigTitle.push_back("BHBH550"); 	sigfullName.push_back("BprimeBprimeToBHBHinc_M-550"); 
	sigTitle.push_back("BHBH600"); 	sigfullName.push_back("BprimeBprimeToBHBHinc_M-600"); 
	//sigTitle.push_back("BHBH650");	sigfullName.push_back("BprimeBprimeToBHBHinc_M-650"); 
	sigTitle.push_back("BHBH700"); 	sigfullName.push_back("BprimeBprimeToBHBHinc_M-700"); 
	//sigTitle.push_back("BHBH750"); 	sigfullName.push_back("BprimeBprimeToBHBHinc_M-750"); 
	sigTitle.push_back("BHBH800"); 	sigfullName.push_back("BprimeBprimeToBHBHinc_M-800"); 
	sigTitle.push_back("BHBH900"); 	sigfullName.push_back("BprimeBprimeToBHBHinc_M-900"); 
	sigTitle.push_back("BHBH1000");	sigfullName.push_back("BprimeBprimeToBHBHinc_M-1000"); 
	sigTitle.push_back("BHBH1200");	sigfullName.push_back("BprimeBprimeToBHBHinc_M-1200"); 
	//sigTitle.push_back("BHBH1500");	sigfullName.push_back("BprimeBprimeToBHBHinc_M-1500");
	const int sigSize= sigTitle.size();

	// Load files
	TFile* fileMean = new TFile((loadpath+"/"+meanworkspace+"/"+meanfileName).c_str());
	TFile* fileSF[SFSize];
	for( int i=0; i<SFSize; i++ ){
		string filePath = loadpath+"/"+SFworkspace[i]+"/"+SFfileName[i];
		cout<<"Load "<<filePath<<endl;
		fileSF[i] = new TFile(filePath.c_str());
	}
 
	//BprimeBprimeToBHBHinc_M-800__HTSel_2b_HTAK5
	const int hNameSize = hName.size();
	for( int h=0; h<hNameSize; h++){
		cout<<"Make template for "<<hName[h]<<" with rebin "<<rebin[h]<<"..."<<endl;
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
			hData[A][c] = (TH1D*)fileMean->Get((datafullName+"__ABCDana"+cate[c]+"_"+hName[h]+"_A").c_str()); 
			hData[B][c] = (TH1D*)fileMean->Get((datafullName+"__ABCDana"+cate[c]+"_"+hName[h]+"_B").c_str()); 
			hData[C][c] = (TH1D*)fileMean->Get((datafullName+"__ABCDana"+cate[c]+"_"+hName[h]+"_C").c_str()); 
			hData[D][c] = (TH1D*)fileMean->Get((datafullName+"__ABCDana"+cate[c]+"_"+hName[h]+"_D").c_str()); 
			fix(hData[A][c]);
			fix(hData[B][c]);
			fix(hData[C][c]);
			fix(hData[D][c]);
			for( int s=0; s<sigSize; s++ ){
				hSig[A][c][s]  = (TH1D*)fileMean->Get((sigfullName[s]+"__ABCDana"+cate[c]+"_"+hName[h]+"_A").c_str()); 
				hSig[B][c][s]  = (TH1D*)fileMean->Get((sigfullName[s]+"__ABCDana"+cate[c]+"_"+hName[h]+"_B").c_str()); 
				hSig[C][c][s]  = (TH1D*)fileMean->Get((sigfullName[s]+"__ABCDana"+cate[c]+"_"+hName[h]+"_C").c_str()); 
				hSig[D][c][s]  = (TH1D*)fileMean->Get((sigfullName[s]+"__ABCDana"+cate[c]+"_"+hName[h]+"_D").c_str()); 
				fix(hSig[A][c][s]);
				fix(hSig[B][c][s]);
				fix(hSig[C][c][s]);
				fix(hSig[D][c][s]);
				for( int sf=0; sf<SFSize; sf++ ){
					if( sf < SFSize-2 ){
						// Load all region of distribution, but only use for observeRegion
						hSigSF[A][c][s][sf]  = (TH1D*)fileSF[sf]->Get((sigfullName[s]+"__ABCDana"+cate[c]+"_"+hName[h]+"_A").c_str()); 
						hSigSF[B][c][s][sf]  = (TH1D*)fileSF[sf]->Get((sigfullName[s]+"__ABCDana"+cate[c]+"_"+hName[h]+"_B").c_str()); 
						hSigSF[C][c][s][sf]  = (TH1D*)fileSF[sf]->Get((sigfullName[s]+"__ABCDana"+cate[c]+"_"+hName[h]+"_C").c_str()); 
						hSigSF[D][c][s][sf]  = (TH1D*)fileSF[sf]->Get((sigfullName[s]+"__ABCDana"+cate[c]+"_"+hName[h]+"_D").c_str());
						fix(hSigSF[A][c][s][sf]);
						fix(hSigSF[B][c][s][sf]);
						fix(hSigSF[C][c][s][sf]);
						fix(hSigSF[D][c][s][sf]);
					}else{
						//For bug in root file
						if( s == 4 ) hSigSF[observeRegion][c][s][sf]  = (TH1D*)fileSF[sf]->Get(("/"+sigfullName[s]+"__HTSel"+cate[c]+"_"+hName[h]+"AK5").c_str()); 
						else hSigSF[observeRegion][c][s][sf]  = (TH1D*)fileSF[sf]->Get((sigfullName[s]+"__HTSel"+cate[c]+"_"+hName[h]+"AK5").c_str());
						fix(hSigSF[observeRegion][c][s][sf]);
					}
				}
			}
		}

		// Store plot to new file
		cout<<"	Data..."<<endl;
		outhData[_nomal][_obs] = new TH1D("CatAll_data_obs", "", newbinsize[h], newbinmin[h], newbinmax[h]);
		outhData[_1b][_obs] = new TH1D("Cat1b_data_obs", "", newbinsize[h], newbinmin[h], newbinmax[h]);
		outhData[_2b][_obs] = new TH1D("Cat2b_data_obs", "", newbinsize[h], newbinmin[h], newbinmax[h]);
		newbinPlots( hData[observeRegion][_nomal], outhData[_nomal][_obs], rebin[h], catalldoOverFlow[h], catallOverFlowValue[h] );
		newbinPlots( hData[observeRegion][_1b], outhData[_1b][_obs], rebin[h], cat1bdoOverFlow[h], cat1bOverFlowValue[h]);
		newbinPlots( hData[observeRegion][_2b], outhData[_2b][_obs], rebin[h], cat2bdoOverFlow[h], cat2bOverFlowValue[h]);
	
		TH1D* catallBkg = new TH1D("CatAll_background", "", newbinsize[h], newbinmin[h], newbinmax[h]);
		TH1D* cat1bBkg = new TH1D("Cat1b_background", "", newbinsize[h], newbinmin[h], newbinmax[h]);
		TH1D* cat2bBkg = new TH1D("Cat2b_background", "", newbinsize[h], newbinmin[h], newbinmax[h]);
		newbinPlots( hData[A][_nomal], catallBkg, rebin[h], catalldoOverFlow[h], catallOverFlowValue[h] );
		newbinPlots( hData[A][_1b],    cat1bBkg,  rebin[h], cat1bdoOverFlow[h], cat1bOverFlowValue[h]);
		newbinPlots( hData[A][_2b],    cat2bBkg,  rebin[h], cat2bdoOverFlow[h], cat2bOverFlowValue[h]);
		
		outhData[_nomal][_exp] = (TH1D*)drawExtrap( hDataCut[_nomal], D, C, catallBkg);
		GetSigma(outhData[_nomal][_exp],"CatAll_background_SystUp", 		 1);
		GetSigma(outhData[_nomal][_exp],"CatAll_background_SystDown", 		-1);

		outhData[_1b][_exp] = (TH1D*)drawExtrap( hDataCut[_1b], D, C, cat1bBkg);
		GetSigma(outhData[_1b][_exp], 	"Cat1b_background_SystUp", 	 1);
		GetSigma(outhData[_1b][_exp], 	"Cat1b_background_SystDown", 	-1);

		outhData[_2b][_exp] = (TH1D*)drawExtrap( hDataCut[_2b], D, C, cat2bBkg);
		GetSigma(outhData[_2b][_exp], 	"Cat2b_background_SystUp", 	 1);
		GetSigma(outhData[_2b][_exp], 	"Cat2b_background_SystDown",	-1);

		for( int s=0; s<sigSize; s++ ){
			// Combined
			outhSig[_nomal][s] = new TH1D(("CatAll_"+sigTitle[s]).c_str(), "", newbinsize[h], newbinmin[h], newbinmax[h]);
			newbinPlots( hSig[observeRegion][_nomal][s], outhSig[_nomal][s], rebin[h], catalldoOverFlow[h], catallOverFlowValue[h] );
			for( int sf=0; sf<SFSize; sf++){
				outhSigSF[_nomal][s][sf] = new TH1D(("CatAll_"+sigTitle[s]+"_"+SFTitle[sf]).c_str(), "", newbinsize[h], newbinmin[h], newbinmax[h]);
				newbinPlots( hSigSF[observeRegion][_nomal][s][sf], outhSigSF[_nomal][s][sf], rebin[h], catalldoOverFlow[h], catallOverFlowValue[h]);
			}
 			// 1b
			outhSig[_1b][s] = new TH1D(("Cat1b_"+sigTitle[s]).c_str(), "", newbinsize[h], newbinmin[h], newbinmax[h]);
			newbinPlots( hSig[observeRegion][_1b][s], outhSig[_1b][s], rebin[h], cat1bdoOverFlow[h], cat1bOverFlowValue[h] );
			for( int sf=0; sf<SFSize; sf++){
				outhSigSF[_1b][s][sf] = new TH1D(("Cat1b_"+sigTitle[s]+"_"+SFTitle[sf]).c_str(), "", newbinsize[h], newbinmin[h], newbinmax[h]);
				newbinPlots( hSigSF[observeRegion][_1b][s][sf], outhSigSF[_1b][s][sf], rebin[h], cat1bdoOverFlow[h], cat1bOverFlowValue[h] );
			}
			//2b
			outhSig[_2b][s] = new TH1D(("Cat2b_"+sigTitle[s]).c_str(), "", newbinsize[h], newbinmin[h], newbinmax[h]);
			newbinPlots( hSig[observeRegion][_2b][s], outhSig[_2b][s], rebin[h], cat2bdoOverFlow[h], cat2bOverFlowValue[h] );
			for( int sf=0; sf<SFSize; sf++){
				outhSigSF[_2b][s][sf] = new TH1D(("Cat2b_"+sigTitle[s]+"_"+SFTitle[sf]).c_str(), "", newbinsize[h], newbinmin[h], newbinmax[h]);
				newbinPlots( hSigSF[observeRegion][_2b][s][sf], outhSigSF[_2b][s][sf], rebin[h], cat2bdoOverFlow[h], cat2bOverFlowValue[h] );
			}
		}
		fs->Write();
	}	
}
void GetSigma( TH1D* h, string name, int shift ){
	if( abs(shift)<1 || abs(shift)>1 ) return;
	TH1D* h_ = (TH1D*)h->Clone(name.c_str()); 
	int maxBin = h->GetXaxis()->GetLast();
	int minBin = 1;
	for( int i=minBin; i<=maxBin; i++){
		double con_(0);
		if( shift == 1 ) con_ = h->GetBinContent(i) + h->GetBinError(i);
		if( shift == -1 ) con_ = h->GetBinContent(i) - h->GetBinError(i);
		h_->SetBinContent(i, con_);	
	}
}
void newbinPlots( TH1D* h1, TH1D* h2, int rebin, bool doOverFlow, double overFlowValue ){
	if( h1->GetBinWidth(1) != h2->GetBinWidth(1) ){
		cout<<"ERROR: Two different bin width, h1 "<<h1->GetBinWidth(1)<<", h2 "<<h2->GetBinWidth(1)<<endl;
		return;
	}	
	double binwidth = h1->GetBinWidth(1);
	double maxbin1 = h1->GetXaxis()->GetLast();	
	double minbin1  = 1;	
	double minbin2  = 1;	
	double minbinVal1 = h1->GetBinLowEdge(minbin1);	
	double minbinVal2 = h2->GetBinLowEdge(minbin2);
	if( minbinVal1 < minbinVal2 ){
		cout<<"ERROR: h1's minimum value of bin has to be small then h2's"<<endl;
		return;
	} 
	int binshift = (minbinVal1-minbinVal2)/binwidth;
	for( int b1=1; b1<=maxbin1; b1++ ){
		int b2 = b1 + binshift;
		h2->SetBinContent(b2, h1->GetBinContent(b1));
		h2->SetBinError(b2, h1->GetBinError(b1));
	}
	h2->Rebin(rebin);
	
	// overflow
	if( doOverFlow ){
		int checkBin = int(((overFlowValue-minbinVal2)/h2->GetBinWidth(1))*10)%10;
		if( checkBin != 0 ){
			cout<<"ERROR: The overflow region is wrong with bin "<<h2->GetBinWidth(1)<<", from min bin "<<minbinVal2<<endl;
			return;
		} 
		int overflowBin(0);
		double overflowBinContent(0);
		double overflowBinSumw2(0);
		double maxbin2 = h2->GetXaxis()->GetLast();
		for( int b2=1; b2<=maxbin2; b2++){
			if( h2->GetBinLowEdge(b2) >= overFlowValue ){
				if( h2->GetBinLowEdge(b2) == overFlowValue ) overflowBin=b2;
				overflowBinContent += h2->GetBinContent(b2);
				overflowBinSumw2 += h2->GetBinError(b2)*h2->GetBinError(b2);
				h2->SetBinContent(b2, 0);
				h2->SetBinError(b2, 0);
			}
		}
		h2->SetBinContent(overflowBin, overflowBinContent);
		h2->SetBinError(overflowBin, sqrt(overflowBinSumw2));
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
