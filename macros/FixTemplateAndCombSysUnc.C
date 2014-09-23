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

const int _mean = 0;
const int _up = 1;
const int _down = 2;

void overFlow( TH1D* h, double xmin, double xmax );
void clonePlots( TH1D* h1, string outname, double xmin, double xmax );
void clonePlotsCombineUncMean(TH1D* h_mean, string outname, vector<TH1D*> h_uncUP, vector<TH1D*> h_uncDown, double xmin, double xmax );
void clonePlotsCombineUnc(TH1D* h_mean, string outname, vector<TH1D*> h_unc, int sigma, double xmin, double xmax );
void combinePlots( TH1D* h_mean, TH1D* h1, TH1D* h2, string outname, int sigma, double xmin, double xmax );

void FixTemplateAndCombSysUnc(){
	//string loadpath = "results/addTopReWrt_02Aug";
	//string savepath = "results/addTopReWrt_02Aug";
	string loadpath = "results/ForScanBprimeBR_17Sep";
	string savepath = "results/ForScanBprimeBR_17Sep";

	string fileName="ABCDResultTemplate_CutTTJet2_HT.root";

	TFile* f_in = new TFile((loadpath+"/"+fileName).c_str());
	TFile* f_out = new TFile((savepath+"/Template_HT.root").c_str(), "RECREATE");

	TH1D *_data[_cat], *_dataSubtractTT[_cat], *_bkg[_cat][3], *_bkgFullData[_cat], *_QCDData[_cat];
	double xMin[_cat]={950, 950, 950};
	double xMax[_cat]={1850, 1850, 1550};

	cout<<"Clone data..."<<endl;
	_data[_nomal]=(TH1D*)f_in->Get("CatAll_data_obs"); 	clonePlots(_data[_nomal], "CatAll_data_obs", xMin[_nomal], xMax[_nomal]);
	_data[_1b]   =(TH1D*)f_in->Get("Cat1b_data_obs"); 	clonePlots(_data[_1b], "Cat1b_data_obs", xMin[_1b], xMax[_1b]);
	_data[_2b]   =(TH1D*)f_in->Get("Cat2b_data_obs"); 	clonePlots(_data[_2b], "Cat2b_data_obs", xMin[_2b], xMax[_2b]);

	_dataSubtractTT[_nomal]=(TH1D*)f_in->Get("CatAll_data_SubtractTTJets"); clonePlots(_dataSubtractTT[_nomal], "CatAll_data_SubtractTTJets", xMin[_nomal], xMax[_nomal]);
	_dataSubtractTT[_1b]   =(TH1D*)f_in->Get("Cat1b_data_SubtractTTJets"); 	clonePlots(_dataSubtractTT[_1b], "Cat1b_data_SubtractTTJets", xMin[_1b], xMax[_1b]);
	_dataSubtractTT[_2b]   =(TH1D*)f_in->Get("Cat2b_data_SubtractTTJets"); 	clonePlots(_dataSubtractTT[_2b], "Cat2b_data_SubtractTTJets", xMin[_2b], xMax[_2b]);

	cout<<"Clone background..."<<endl;
	/*_bkg[_nomal][_mean]=(TH1D*)f_in->Get("CatAll_background"); 		clonePlots(_bkg[_nomal][_mean], "CatAll_background", xMin[_nomal], xMax[_nomal]);
	_bkg[_nomal][_up]  =(TH1D*)f_in->Get("CatAll_background_StatUp"); 	clonePlots(_bkg[_nomal][_up],   "CatAll_background_StatUp", xMin[_nomal], xMax[_nomal]);
	_bkg[_nomal][_down]=(TH1D*)f_in->Get("CatAll_background_StatDown"); 	clonePlots(_bkg[_nomal][_down], "CatAll_background_StatDown", xMin[_nomal], xMax[_nomal]);
	_bkg[_1b][_mean]=(TH1D*)f_in->Get("Cat1b_background"); 			clonePlots(_bkg[_1b][_mean], "Cat1b_background", xMin[_1b], xMax[_1b]);
	_bkg[_1b][_up]  =(TH1D*)f_in->Get("Cat1b_background_StatUp");	 	clonePlots(_bkg[_1b][_up],   "Cat1b_background_StatUp", xMin[_1b], xMax[_1b]);
	_bkg[_1b][_down]=(TH1D*)f_in->Get("Cat1b_background_StatDown");	 	clonePlots(_bkg[_1b][_down], "Cat1b_background_StatDown", xMin[_1b], xMax[_1b]);
	_bkg[_2b][_mean]=(TH1D*)f_in->Get("Cat2b_background"); 			clonePlots(_bkg[_2b][_mean], "Cat2b_background", xMin[_2b], xMax[_2b]);
	_bkg[_2b][_up]  =(TH1D*)f_in->Get("Cat2b_background_StatUp");	 	clonePlots(_bkg[_2b][_up],   "Cat2b_background_StatUp", xMin[_2b], xMax[_2b]);
	_bkg[_2b][_down]=(TH1D*)f_in->Get("Cat2b_background_StatDown");	 	clonePlots(_bkg[_2b][_down], "Cat2b_background_StatDown", xMin[_2b], xMax[_2b]);*/

	_QCDData[_nomal]=(TH1D*)f_in->Get("CatAll_QCD_data");	clonePlots(_QCDData[_nomal], "CatAll_QCD_data", xMin[_nomal], xMax[_nomal]);
	_QCDData[_1b]=(TH1D*)f_in->Get("Cat1b_QCD_data");	clonePlots(_QCDData[_1b], "Cat1b_QCD_data", xMin[_1b], xMax[_1b]);
	_QCDData[_2b]=(TH1D*)f_in->Get("Cat2b_QCD_data");	clonePlots(_QCDData[_2b], "Cat2b_QCD_data", xMin[_2b], xMax[_2b]);
	
	_bkgFullData[_nomal]=(TH1D*)f_in->Get("CatAll_background_fulldata");	clonePlots(_bkgFullData[_nomal], "CatAll_background_fulldata", xMin[_nomal], xMax[_nomal]);
	_bkgFullData[_1b]=(TH1D*)f_in->Get("Cat1b_background_fulldata");	clonePlots(_bkgFullData[_1b], "Cat1b_background_fulldata", xMin[_1b], xMax[_1b]);
	_bkgFullData[_2b]=(TH1D*)f_in->Get("Cat2b_background_fulldata");	clonePlots(_bkgFullData[_2b], "Cat2b_background_fulldata", xMin[_2b], xMax[_2b]);

	string cate[_cat]={"CatAll", "Cat1b", "Cat2b"};

	vector<string> sampleName;
	sampleName.push_back("background"); 	
	sampleName.push_back("TTJets"); 	
	sampleName.push_back("BHBH500"); 	
	sampleName.push_back("BHBH600"); 	
	sampleName.push_back("BHBH700"); 	
	sampleName.push_back("BHBH800"); 	
	sampleName.push_back("BHBH900"); 	
	sampleName.push_back("BHBH1000"); 	
	sampleName.push_back("BHBH1200"); 
	sampleName.push_back("BHBZ500"); 	
	sampleName.push_back("BHBZ600"); 	
	sampleName.push_back("BHBZ700"); 	
	sampleName.push_back("BHBZ800"); 	
	sampleName.push_back("BHBZ900"); 	
	sampleName.push_back("BHBZ1000"); 	
	sampleName.push_back("BHBZ1200"); 
	sampleName.push_back("BHTW500"); 	
	sampleName.push_back("BHTW600"); 	
	sampleName.push_back("BHTW700"); 	
	sampleName.push_back("BHTW800"); 	
	sampleName.push_back("BHTW900"); 	
	sampleName.push_back("BHTW1000"); 	
	sampleName.push_back("BHTW1200"); 
	sampleName.push_back("TWTW500"); 	
	sampleName.push_back("TWTW600"); 	
	sampleName.push_back("TWTW700"); 	
	sampleName.push_back("TWTW800"); 	
	sampleName.push_back("TWTW900"); 	
	sampleName.push_back("TWTW1000"); 	
	sampleName.push_back("TWTW1200"); 
	sampleName.push_back("BZTW500"); 	
	sampleName.push_back("BZTW600"); 	
	sampleName.push_back("BZTW700"); 	
	sampleName.push_back("BZTW800"); 	
	sampleName.push_back("BZTW900"); 	
	sampleName.push_back("BZTW1000"); 	
	sampleName.push_back("BZTW1200"); 
	sampleName.push_back("BZBZ500"); 	
	sampleName.push_back("BZBZ600"); 	
	sampleName.push_back("BZBZ700"); 	
	sampleName.push_back("BZBZ800"); 	
	sampleName.push_back("BZBZ900"); 	
	sampleName.push_back("BZBZ1000"); 	
	sampleName.push_back("BZBZ1200"); 
	const int sampleNameSize=sampleName.size();	
	enum UNCNAME{ Stat, JER, JES, SFb, SFl, CA8, PU, TopPtReWrt, TTMatch, TTScale, TotalUnc};

	TH1D *_Sample[sampleNameSize][_cat], *_SampleUNC[sampleNameSize][_cat][TotalUnc][2];

	for( int s=0; s<sampleNameSize; s++){
		cout<<"Clone "<<sampleName[s]<<"..."<<endl;
		for( int c=0; c<_cat; c++ ){

			vector<TH1D*> UNCup, UNCdown;
			_Sample[s][c] = (TH1D*)f_in->Get((cate[c]+"_"+sampleName[s]).c_str());
			_SampleUNC[s][c][Stat][0] = (TH1D*)f_in->Get((cate[c]+"_"+sampleName[s]+"_StatUp").c_str());		
			_SampleUNC[s][c][Stat][1] = (TH1D*)f_in->Get((cate[c]+"_"+sampleName[s]+"_StatDown").c_str());		
			_SampleUNC[s][c][JER][0] = (TH1D*)f_in->Get((cate[c]+"_"+sampleName[s]+"_JERUp").c_str());		
			_SampleUNC[s][c][JER][1] = (TH1D*)f_in->Get((cate[c]+"_"+sampleName[s]+"_JERDown").c_str());		
			_SampleUNC[s][c][JES][0] = (TH1D*)f_in->Get((cate[c]+"_"+sampleName[s]+"_JESUp").c_str());		
			_SampleUNC[s][c][JES][1] = (TH1D*)f_in->Get((cate[c]+"_"+sampleName[s]+"_JESDown").c_str());		
			_SampleUNC[s][c][SFb][0] = (TH1D*)f_in->Get((cate[c]+"_"+sampleName[s]+"_SFbUp").c_str());		
			_SampleUNC[s][c][SFb][1] = (TH1D*)f_in->Get((cate[c]+"_"+sampleName[s]+"_SFbDown").c_str());		
			_SampleUNC[s][c][SFl][0] = (TH1D*)f_in->Get((cate[c]+"_"+sampleName[s]+"_SFlUp").c_str());	
			_SampleUNC[s][c][SFl][1] = (TH1D*)f_in->Get((cate[c]+"_"+sampleName[s]+"_SFlDown").c_str());	
			_SampleUNC[s][c][CA8][0] = (TH1D*)f_in->Get((cate[c]+"_"+sampleName[s]+"_CA8Up").c_str());	
			_SampleUNC[s][c][CA8][1] = (TH1D*)f_in->Get((cate[c]+"_"+sampleName[s]+"_CA8Down").c_str());
			_SampleUNC[s][c][PU][0] = (TH1D*)f_in->Get((cate[c]+"_"+sampleName[s]+"_PUUp").c_str());	
			_SampleUNC[s][c][PU][1] = (TH1D*)f_in->Get((cate[c]+"_"+sampleName[s]+"_PUDown").c_str());
			if( s<2 ){
				_SampleUNC[s][c][TopPtReWrt][0] = (TH1D*)f_in->Get((cate[c]+"_"+sampleName[s]+"_TopPtReWrtUp").c_str());
				_SampleUNC[s][c][TopPtReWrt][1] = (TH1D*)f_in->Get((cate[c]+"_"+sampleName[s]+"_TopPtReWrtDown").c_str());	
				_SampleUNC[s][c][TTMatch][0] = (TH1D*)f_in->Get((cate[c]+"_"+sampleName[s]+"_TTJetsMatchingUp").c_str());	
				_SampleUNC[s][c][TTMatch][1] = (TH1D*)f_in->Get((cate[c]+"_"+sampleName[s]+"_TTJetsMatchingDown").c_str());	
				_SampleUNC[s][c][TTScale][0] = (TH1D*)f_in->Get((cate[c]+"_"+sampleName[s]+"_TTJetsScaleUp").c_str());
				_SampleUNC[s][c][TTScale][1] = (TH1D*)f_in->Get((cate[c]+"_"+sampleName[s]+"_TTJetsScaleDown").c_str());
				for( int unc=0; unc<TotalUnc; unc++){
					UNCup.push_back(_SampleUNC[s][c][unc][0]);
					UNCdown.push_back(_SampleUNC[s][c][unc][1]);
				}
			}else{
				for( int unc=0; unc<(TotalUnc-3); unc++){
					UNCup.push_back(_SampleUNC[s][c][unc][0]);
					UNCdown.push_back(_SampleUNC[s][c][unc][1]);
				}
			}

			clonePlots(_Sample[s][c], cate[c]+"_"+sampleName[s], xMin[c], xMax[c]);
			clonePlotsCombineUncMean(_Sample[s][c], cate[c]+"_"+sampleName[s]+"_TotalUnc", UNCup, UNCdown, xMin[c], xMax[c]);
			clonePlotsCombineUnc(_Sample[s][c], cate[c]+"_"+sampleName[s]+"_TotalUncUp", UNCup, 1, xMin[c], xMax[c]);
			clonePlotsCombineUnc(_Sample[s][c], cate[c]+"_"+sampleName[s]+"_TotalUncDown", UNCdown, -1, xMin[c], xMax[c]);
			if( s<2 ){
				clonePlots(_SampleUNC[s][c][Stat][0], cate[c]+"_"+sampleName[s]+"_StatUp", xMin[c], xMax[c]);	
				clonePlots(_SampleUNC[s][c][Stat][1], cate[c]+"_"+sampleName[s]+"_StatDown", xMin[c], xMax[c]);	
				clonePlots(_SampleUNC[s][c][TopPtReWrt][0], cate[c]+"_"+sampleName[s]+"_TopPtReWrtUp", xMin[c], xMax[c]);	
				clonePlots(_SampleUNC[s][c][TopPtReWrt][1], cate[c]+"_"+sampleName[s]+"_TopPtReWrtDown", xMin[c], xMax[c]);	
				clonePlots(_SampleUNC[s][c][TTMatch][0], cate[c]+"_"+sampleName[s]+"_TTJetsMatchingUp", xMin[c], xMax[c]);	
				clonePlots(_SampleUNC[s][c][TTMatch][1], cate[c]+"_"+sampleName[s]+"_TTJetsMatchingDown", xMin[c], xMax[c]);	
				clonePlots(_SampleUNC[s][c][TTScale][0], cate[c]+"_"+sampleName[s]+"_TTJetsScaleUp", xMin[c], xMax[c]);	
				clonePlots(_SampleUNC[s][c][TTScale][1], cate[c]+"_"+sampleName[s]+"_TTJetsScaleDown", xMin[c], xMax[c]);	
			}else{
				clonePlots(_SampleUNC[s][c][Stat][0], cate[c]+"_"+sampleName[s]+"_SigStatUp", xMin[c], xMax[c]);	
				clonePlots(_SampleUNC[s][c][Stat][1], cate[c]+"_"+sampleName[s]+"_SigStatDown", xMin[c], xMax[c]);	
			}
			clonePlots(_SampleUNC[s][c][JER][0], cate[c]+"_"+sampleName[s]+"_JERUp", xMin[c], xMax[c]);	
			clonePlots(_SampleUNC[s][c][JER][1], cate[c]+"_"+sampleName[s]+"_JERDown", xMin[c], xMax[c]);	
			clonePlots(_SampleUNC[s][c][JES][0], cate[c]+"_"+sampleName[s]+"_JESUp", xMin[c], xMax[c]);	
			clonePlots(_SampleUNC[s][c][JES][1], cate[c]+"_"+sampleName[s]+"_JESDown", xMin[c], xMax[c]);	
			clonePlots(_SampleUNC[s][c][SFb][0], cate[c]+"_"+sampleName[s]+"_SFbUp", xMin[c], xMax[c]);	
			clonePlots(_SampleUNC[s][c][SFb][1], cate[c]+"_"+sampleName[s]+"_SFbDown", xMin[c], xMax[c]);	
			clonePlots(_SampleUNC[s][c][SFl][0], cate[c]+"_"+sampleName[s]+"_SFlUp", xMin[c], xMax[c]);	
			clonePlots(_SampleUNC[s][c][SFl][1], cate[c]+"_"+sampleName[s]+"_SFlDown", xMin[c], xMax[c]);	
			clonePlots(_SampleUNC[s][c][CA8][0], cate[c]+"_"+sampleName[s]+"_CA8Up", xMin[c], xMax[c]);	
			clonePlots(_SampleUNC[s][c][CA8][1], cate[c]+"_"+sampleName[s]+"_CA8Down", xMin[c], xMax[c]);	
			clonePlots(_SampleUNC[s][c][PU][0], cate[c]+"_"+sampleName[s]+"_PUUp", xMin[c], xMax[c]);	
			clonePlots(_SampleUNC[s][c][PU][1], cate[c]+"_"+sampleName[s]+"_PUDown", xMin[c], xMax[c]);
	
			combinePlots( _Sample[s][c], _SampleUNC[s][c][JER][0], _SampleUNC[s][c][JES][0], cate[c]+"_"+sampleName[s]+"_JECUp",      1, xMin[c], xMax[c]);		
			combinePlots( _Sample[s][c], _SampleUNC[s][c][JER][1], _SampleUNC[s][c][JES][1], cate[c]+"_"+sampleName[s]+"_JECDown",   -1, xMin[c], xMax[c]);		
			combinePlots( _Sample[s][c], _SampleUNC[s][c][SFb][0], _SampleUNC[s][c][SFl][0], cate[c]+"_"+sampleName[s]+"_bTagSFUp",   1, xMin[c], xMax[c]);		
			combinePlots( _Sample[s][c], _SampleUNC[s][c][SFb][1], _SampleUNC[s][c][SFl][1], cate[c]+"_"+sampleName[s]+"_bTagSFDown",-1, xMin[c], xMax[c]);		
//			combinePlots( _Sample[s][c], _SampleUNC[s][c][SFHb][0], _SampleUNC[s][c][SFHl][0], cate[c]+"_"+sampleName[s]+"_HTagSFUp",   1);		
//			combinePlots( _Sample[s][c], _SampleUNC[s][c][SFHb][1], _SampleUNC[s][c][SFHl][1], cate[c]+"_"+sampleName[s]+"_HTagSFDown",-1);	
/*
			clonePlots(_SampleUNC[s][c][CA8][0], cate[c]+"_"+sampleName[s]+"_CA8Up", xMin[c], xMax[c]);	
			clonePlots(_SampleUNC[s][c][CA8][1], cate[c]+"_"+sampleName[s]+"_CA8Down", xMin[c], xMax[c]);	
			clonePlots(_SampleUNC[s][c][PU][0], cate[c]+"_"+sampleName[s]+"_PUUp", xMin[c], xMax[c]);	
			clonePlots(_SampleUNC[s][c][PU][1], cate[c]+"_"+sampleName[s]+"_PUDown", xMin[c], xMax[c]);	
			*/
		}	
	}
	f_out->Write();

}
void overFlow( TH1D* h, double xmin, double xmax ){
	int binMax = h->GetXaxis()->GetLast();
	int binMin = 1;
	if( h->GetXaxis()->GetBinLowEdge(binMin) < xmin ){
		double content=0;
		double sumw2=0;
		int bin;
		for( int i=1; h->GetXaxis()->GetBinLowEdge(i)<=xmin; i++){
			content += h->GetBinContent(i);
			sumw2   += h->GetBinError(i)*h->GetBinError(i);
			if( h->GetXaxis()->GetBinLowEdge(i)== xmin ) bin=i; 
		}
		h->SetBinContent(bin, content );	
		h->SetBinError(bin,sqrt(sumw2) );	
	} 
	if( h->GetXaxis()->GetBinLowEdge(binMax) > xmax ){
		double content=0;
		double sumw2=0;
		int bin;
		for( int i=binMax; h->GetXaxis()->GetBinLowEdge(i)>=xmax; i--){
			content += h->GetBinContent(i);
			sumw2   += h->GetBinError(i)*h->GetBinError(i);
			if( h->GetXaxis()->GetBinLowEdge(i)== xmax ) bin=i;
			h->SetBinContent(i, 0); 
			h->SetBinError(i, 0); 
		}
		h->SetBinContent(bin, content );	
		h->SetBinError(bin,sqrt(sumw2) );	
	}
}
void clonePlots( TH1D* h1, string outname, double xmin, double xmax ){
	TH1D* h = (TH1D*)h1->Clone(outname.c_str());
	overFlow(h, xmin, xmax);
}
void clonePlotsCombineUncMean(TH1D* h_mean, string outname, vector<TH1D*> h_uncUP, vector<TH1D*> h_uncDown, double xmin, double xmax ){
	TH1D* h = (TH1D*)h_mean->Clone(outname.c_str());
	const int maxbin = h->GetXaxis()->GetLast();
	const int uncSize = h_uncUP.size();
	for( int bin=1; bin<=maxbin; bin++){
		double newError = 0.;
		double varSumw2 = 0.;
		for( int unc=0; unc<uncSize; unc++){
			//cout<<unc<<"/"<<uncSize<<endl;
			double meanError = (abs(h_uncUP[unc]->GetBinContent(bin)-h_mean->GetBinContent(bin))+abs(h_uncDown[unc]->GetBinContent(bin)-h_mean->GetBinContent(bin)))/2; 
			varSumw2 += meanError*meanError; 
		}
		newError = sqrt(varSumw2);
		h->SetBinError(bin, newError);
	}	
	overFlow(h, xmin, xmax);
}
void clonePlotsCombineUnc(TH1D* h_mean, string outname, vector<TH1D*> h_unc, int sigma, double xmin, double xmax ){
	TH1D* h = (TH1D*)h_mean->Clone(outname.c_str());
	const int maxbin = h->GetXaxis()->GetLast();
	const int uncSize = h_unc.size();
	for( int bin=1; bin<=maxbin; bin++){
		double newContent = 0.;
		double varSumw2 = 0.;
		for( int unc=0; unc<uncSize; unc++){
			varSumw2 += (h_unc[unc]->GetBinContent(bin)-h_mean->GetBinContent(bin))*(h_unc[unc]->GetBinContent(bin)-h_mean->GetBinContent(bin)); 
		}
		if( sigma > 0 ) newContent = h_mean->GetBinContent(bin) + sqrt(varSumw2);
		else if( sigma < 0 ) newContent = h_mean->GetBinContent(bin) - sqrt(varSumw2);
		h->SetBinContent(bin, newContent);
	}	
	overFlow(h, xmin, xmax);
}
void combinePlots( TH1D* h_mean, TH1D* h1, TH1D* h2, string outname, int sigma, double xmin, double xmax ){
	if( h_mean->GetBinWidth(1) != h2->GetBinWidth(1) || h_mean->GetBinWidth(1) != h1->GetBinWidth(1) ){
		cout<<"ERROR: The different bin width, h_mean "<<h_mean->GetBinWidth(1)<<", h1 "<<h1->GetBinWidth(1)<<", h2 "<<h2->GetBinWidth(1)<<endl;
		return;
	}
	TH1D* h = (TH1D*)h_mean->Clone(outname.c_str());
	double maxbin = h->GetXaxis()->GetLast();	
	for( int b=1; b<=maxbin; b++ ){
		double var1sqr = (h1->GetBinContent(b)-h_mean->GetBinContent(b))*(h1->GetBinContent(b)-h_mean->GetBinContent(b));
		double var2sqr = (h2->GetBinContent(b)-h_mean->GetBinContent(b))*(h2->GetBinContent(b)-h_mean->GetBinContent(b));
		double meanVar = sqrt( var1sqr + var2sqr );
		double newContent(0);
		double meanContent = h_mean->GetBinContent(b);
		if( sigma > 0 ) newContent = meanContent + meanVar;
		else if( sigma < 0 ) newContent = meanContent - meanVar;
		h->SetBinContent(b, newContent);
	}
	overFlow(h, xmin, xmax);
}
