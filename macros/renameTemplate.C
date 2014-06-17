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

void clonePlots( TH1D* h1, string outname );
void combinePlots( TH1D* h_mean, TH1D* h1, TH1D* h2, string outname, int sigma );

void renameTemplate(){
	//string loadpath = "results/bin300GeV";
	//string savepath = "results/bin300GeV";
	string loadpath = "results/bin100GeV";
	string savepath = "results/bin100GeV";

	string fileName="ABCDResultTemplate_ABCDana_HT.root";

	TFile* f_in = new TFile((loadpath+"/"+fileName).c_str());
	TFile* f_out = new TFile((loadpath+"/Template_HT.root").c_str(), "RECREATE");

	TH1D *_data[_cat], *_bkg[_cat][3];

	cout<<"Clone data..."<<endl;
	_data[_nomal]=(TH1D*)f_in->Get("CatAll_data_obs"); 	clonePlots(_data[_nomal], "CatAll_data_obs");
	_data[_1b]   =(TH1D*)f_in->Get("Cat1b_data_obs"); 	clonePlots(_data[_1b], "Cat1b_data_obs");
	_data[_2b]   =(TH1D*)f_in->Get("Cat2b_data_obs"); 	clonePlots(_data[_2b], "Cat2b_data_obs");

	cout<<"Clone background..."<<endl;
	_bkg[_nomal][_mean]=(TH1D*)f_in->Get("CatAll_background"); 		clonePlots(_bkg[_nomal][_mean], "CatAll_background");
	_bkg[_nomal][_up]  =(TH1D*)f_in->Get("CatAll_background_SystUp"); 	clonePlots(_bkg[_nomal][_up],   "CatAll_background_SystUp");
	_bkg[_nomal][_down]=(TH1D*)f_in->Get("CatAll_background_SystDown"); 	clonePlots(_bkg[_nomal][_down], "CatAll_background_SystDown");
	_bkg[_1b][_mean]=(TH1D*)f_in->Get("Cat1b_background"); 			clonePlots(_bkg[_1b][_mean], "Cat1b_background");
	_bkg[_1b][_up]  =(TH1D*)f_in->Get("Cat1b_background_SystUp");	 	clonePlots(_bkg[_1b][_up], 	"Cat1b_background_SystUp");
	_bkg[_1b][_down]=(TH1D*)f_in->Get("Cat1b_background_SystDown");	 	clonePlots(_bkg[_1b][_down], "Cat1b_background_SystDown");
	_bkg[_2b][_mean]=(TH1D*)f_in->Get("Cat2b_background"); 			clonePlots(_bkg[_2b][_mean], "Cat2b_background");
	_bkg[_2b][_up]  =(TH1D*)f_in->Get("Cat2b_background_SystUp");	 	clonePlots(_bkg[_2b][_up], 	"Cat2b_background_SystUp");
	_bkg[_2b][_down]=(TH1D*)f_in->Get("Cat2b_background_SystDown");	 	clonePlots(_bkg[_2b][_down], "Cat2b_background_SystDown");
	
	string cate[_cat]={"CatAll", "Cat1b", "Cat2b"};

	vector<string> sigName;
	sigName.push_back("BHBH500"); 	
	sigName.push_back("BHBH600"); 	
	sigName.push_back("BHBH700"); 	
	sigName.push_back("BHBH800"); 	
	sigName.push_back("BHBH900"); 	
	sigName.push_back("BHBH1000"); 	
	sigName.push_back("BHBH1200"); 
	const int sigNameSize=sigName.size();	
	
	vector<string> SFTitle;  	
	SFTitle.push_back("PU");
	SFTitle.push_back("HM");
	enum SFNAME{ JER, JES, SFb, SFl, SFHb, SFHl, PU, HM, TotalSF};

	TH1D *_Sig[sigNameSize][_cat], *_SigSF[sigNameSize][_cat][TotalSF][2];

	for( int s=0; s<sigNameSize; s++){
		cout<<"Clone Signal "<<sigName[s]<<"..."<<endl;
		for( int c=0; c<_cat; c++ ){
			_Sig[s][c] = (TH1D*)f_in->Get((cate[c]+"_"+sigName[s]).c_str());
			clonePlots(_Sig[s][c], cate[c]+"_"+sigName[s]);

			_SigSF[s][c][JER][0] = (TH1D*)f_in->Get((cate[c]+"_"+sigName[s]+"_JERUp").c_str());		
			_SigSF[s][c][JER][1] = (TH1D*)f_in->Get((cate[c]+"_"+sigName[s]+"_JERDown").c_str());		
			_SigSF[s][c][JES][0] = (TH1D*)f_in->Get((cate[c]+"_"+sigName[s]+"_JESUp").c_str());		
			_SigSF[s][c][JES][1] = (TH1D*)f_in->Get((cate[c]+"_"+sigName[s]+"_JESDown").c_str());		
			_SigSF[s][c][SFb][0] = (TH1D*)f_in->Get((cate[c]+"_"+sigName[s]+"_SFbUp").c_str());		
			_SigSF[s][c][SFb][1] = (TH1D*)f_in->Get((cate[c]+"_"+sigName[s]+"_SFbDown").c_str());		
			_SigSF[s][c][SFl][0] = (TH1D*)f_in->Get((cate[c]+"_"+sigName[s]+"_SFlUp").c_str())	;	
			_SigSF[s][c][SFl][1] = (TH1D*)f_in->Get((cate[c]+"_"+sigName[s]+"_SFlDown").c_str())	;	
			_SigSF[s][c][SFHb][0] = (TH1D*)f_in->Get((cate[c]+"_"+sigName[s]+"_SFHbUp").c_str())	;	
			_SigSF[s][c][SFHb][1] = (TH1D*)f_in->Get((cate[c]+"_"+sigName[s]+"_SFHbDown").c_str())	;	
			_SigSF[s][c][SFHl][0] = (TH1D*)f_in->Get((cate[c]+"_"+sigName[s]+"_SFHlUp").c_str())	;	
			_SigSF[s][c][SFHl][1] = (TH1D*)f_in->Get((cate[c]+"_"+sigName[s]+"_SFHlDown").c_str())	;	
			_SigSF[s][c][PU][0] = (TH1D*)f_in->Get((cate[c]+"_"+sigName[s]+"_PUUp").c_str())	;	
			_SigSF[s][c][PU][1] = (TH1D*)f_in->Get((cate[c]+"_"+sigName[s]+"_PUDown").c_str())	;	
			_SigSF[s][c][HM][0] = (TH1D*)f_in->Get((cate[c]+"_"+sigName[s]+"_HMUp").c_str())	;	
			_SigSF[s][c][HM][1] = (TH1D*)f_in->Get((cate[c]+"_"+sigName[s]+"_HMDown").c_str())	;

			combinePlots( _Sig[s][c], _SigSF[s][c][JER][0], _SigSF[s][c][JES][0], cate[c]+"_"+sigName[s]+"_JECUp",      1);		
			combinePlots( _Sig[s][c], _SigSF[s][c][JER][1], _SigSF[s][c][JES][1], cate[c]+"_"+sigName[s]+"_JECDown",   -1);		
			combinePlots( _Sig[s][c], _SigSF[s][c][SFb][0], _SigSF[s][c][SFl][0], cate[c]+"_"+sigName[s]+"_bTagSFUp",   1);		
			combinePlots( _Sig[s][c], _SigSF[s][c][SFb][1], _SigSF[s][c][SFl][1], cate[c]+"_"+sigName[s]+"_bTagSFDown",-1);		
			combinePlots( _Sig[s][c], _SigSF[s][c][SFHb][0], _SigSF[s][c][SFHl][0], cate[c]+"_"+sigName[s]+"_HTagSFUp",   1);		
			combinePlots( _Sig[s][c], _SigSF[s][c][SFHb][1], _SigSF[s][c][SFHl][1], cate[c]+"_"+sigName[s]+"_HTagSFDown",-1);	

			clonePlots(_SigSF[s][c][PU][0], cate[c]+"_"+sigName[s]+"_PUUp");	
			clonePlots(_SigSF[s][c][PU][1], cate[c]+"_"+sigName[s]+"_PUDown");	
			clonePlots(_SigSF[s][c][HM][0], cate[c]+"_"+sigName[s]+"_HMUp");	
			clonePlots(_SigSF[s][c][HM][1], cate[c]+"_"+sigName[s]+"_HMDown");	
		}	
	}
	f_out->Write();

}
void clonePlots( TH1D* h1, string outname ){
	TH1D* h = (TH1D*)h1->Clone(outname.c_str());
}
void combinePlots( TH1D* h_mean, TH1D* h1, TH1D* h2, string outname, int sigma ){
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
}
