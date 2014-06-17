#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "stdio.h"
#include "TFile.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TROOT.h"
using namespace std;
int const A_=1;
int const B_=2;
int const C_=3;
int const D_=4;
int const Mean=0;
int const Down=1;
int const Up=2;
void printLetex( FILE* output, string hname, vector<string> systType, string loadpath, string filename);
void printLetexCombine( FILE* output, string hname, vector<string> systType, string loadpath, string filename, int colum, vector<string> newType);
void CountSystABCDv2(){

	string loadpath = "results/bin300GeV";
	string savepath = "results/bin300GeV";
	string fileName = "ABCDResultTemplate_ABCDana_HT.root";
	
	vector<string> systType;
	systType.push_back("JER");
	systType.push_back("JES");
	systType.push_back("SFb");
	systType.push_back("SFl");
	systType.push_back("SFHb");
	systType.push_back("SFHl");
	systType.push_back("HM");
	systType.push_back("PU");
	vector<string> newType;
	newType.push_back("JEC");
	newType.push_back("bTagSF");
	newType.push_back("HTagSF");
	enum _cate_{ CatAll, Cat1b, Cat2b, CateSize};
	string cate[CateSize]={"CatAll", "Cat1b", "Cat2b"};

	vector<string> Title; 
	//Title.push_back("BHBH450"); 	
	Title.push_back("BHBH500"); 	
	//Title.push_back("BHBH550"); 	
	Title.push_back("BHBH600"); 	
	//Title.push_back("BHBH650");	
	Title.push_back("BHBH700"); 	
	//Title.push_back("BHBH750"); 	
	Title.push_back("BHBH800"); 	
	Title.push_back("BHBH900"); 	
	Title.push_back("BHBH1000"); 	
	Title.push_back("BHBH1200"); 	
	//Title.push_back("BHBH1500"); 	
	const int titleSize = Title.size();
	
	FILE* pFile[CateSize], *pFileCombine[CateSize];

	for( int i=0; i<CateSize; i++){
		pFile[i] = fopen((savepath+"/SystUncertainties_"+cate[i]+".tex").c_str(),"w");
		cout<<"Print letex format file "<<cate[i]<<" in "<<savepath+"SystUncertainties_"+cate[i]+".tex"<<"..."<<endl;
		for( int k=0; k<titleSize; k++){
			string hname = cate[i]+"_"+Title[k];
			fprintf ( pFile[i], "=============== %s =============== \n",Title[k].c_str());	
			printLetex( pFile[i], hname, systType, loadpath, fileName);
			fprintf ( pFile[i], "\n");	
			fprintf ( pFile[i], "\n");	
		}
		fclose(pFile[i]);
		pFileCombine[i] = fopen((savepath+"/SystUncertaintiesCombine_"+cate[i]+".tex").c_str(),"w");
		cout<<"Print letex format file "<<cate[i]<<" in "<<savepath+"SystUncertaintiesCombine_"+cate[i]+".tex"<<"..."<<endl;
		for( int k=0; k<titleSize; k++){
			string hname = cate[i]+"_"+Title[k];
			fprintf ( pFileCombine[i], "=============== %s =============== \n",Title[k].c_str());	
			printLetexCombine( pFileCombine[i], hname, systType, loadpath, fileName, 3, newType);
			fprintf ( pFileCombine[i], "\n");	
			fprintf ( pFileCombine[i], "\n");	
		}
		fclose(pFileCombine[i]);
	}

}
/* $+1 \sigma$ & $-0.68$ & $+3.27$ & $+1.42$  & $+0.01$ & $+5.04$ & $+0.07$ & $-0.19$& $< 1$\\ */
void printLetexCombine( FILE* output, string hname, vector<string> systType, string loadpath, string filename, int colum, vector<string> newType){
	if( colum != newType.size() ){ 
		cout<<"ERROR: colum size is different from newType!"<<endl;
		return;
	}
	vector<string> systNewType;	
	const int sysTypeSize=systType.size();
	const int sysNewTypeSize=sysTypeSize-colum;
	double pSigma[sysTypeSize], mSigma[sysTypeSize];
	double pNewSigma[sysNewTypeSize], mNewSigma[sysNewTypeSize];
	TFile* f = new TFile( (loadpath + "/" + filename).c_str());

	fprintf ( output, "   sigma    ");
	for( int i=0; i<sysTypeSize; i++){
		TH1D*  h_mean = (TH1D*)f->Get(hname.c_str());
		TH1D*  h_down = (TH1D*)f->Get((hname+"_"+systType[i]+"Down").c_str());
		TH1D*  h_up   = (TH1D*)f->Get((hname+"_"+systType[i]+"Up").c_str());
		pSigma[i] = ( h_up->Integral()   - h_mean->Integral() )/ h_mean->Integral() * 100;
		mSigma[i] = ( h_down->Integral() - h_mean->Integral() )/ h_mean->Integral() * 100;
		//fprintf ( output, "   %s    ", systType[i].c_str());
	}
	for( int i=0; i<colum; i++){
		int t1=i*2;	
		int t2=i*2+1;	
		pNewSigma[i] = 0+sqrt(pSigma[t1]*pSigma[t1]+pSigma[t2]*pSigma[t2] );
		mNewSigma[i] = 0-sqrt(mSigma[t1]*mSigma[t1]+mSigma[t2]*mSigma[t2] );
		systNewType.push_back(newType[i]);
		fprintf ( output, "  %s    ", systNewType[i].c_str());
	}
	for( int i=0; i<sysTypeSize; i++){
		if( i < colum*2 ) continue;
		int j = i - colum;
		systNewType.push_back(systType[i]);
		pNewSigma[j]=pSigma[i];
		mNewSigma[j]=mSigma[i];
		fprintf ( output, "  %s    ", systType[i].c_str());
	}
	fprintf ( output, "  PDF \n");
	fprintf ( output, " $+1 \\sigma$ &");
	for( int i=0; i<sysNewTypeSize; i++){
		fprintf ( output, " $%+1.2f$ &", pNewSigma[i]);
	}
	fprintf ( output, " $< 1$\\\\ \n");
	fprintf ( output, " $-1 \\sigma$ &");
	for( int i=0; i<sysNewTypeSize; i++){
		fprintf ( output, " $%+1.2f$ &", mNewSigma[i]);
	}
	fprintf ( output, " $< 1$\\\\ \n");

}
void printLetex( FILE* output, string hname, vector<string> systType, string loadpath, string filename){
	const int sysTypeSize=systType.size();
	double pSigma[sysTypeSize], mSigma[sysTypeSize];
	TFile* f = new TFile( (loadpath + "/" + filename).c_str());

	fprintf ( output, "    sigma    ");
	for( int i=0; i<sysTypeSize; i++){
		TH1D*  h_mean = (TH1D*)f->Get(hname.c_str());
		TH1D*  h_down = (TH1D*)f->Get((hname+"_"+systType[i]+"Down").c_str());
		TH1D*  h_up   = (TH1D*)f->Get((hname+"_"+systType[i]+"Up").c_str());
		pSigma[i] = ( h_up->Integral()   - h_mean->Integral() )/ h_mean->Integral() * 100;
		mSigma[i] = ( h_down->Integral()   - h_mean->Integral() )/ h_mean->Integral() * 100;
		fprintf ( output, "   %s    ", systType[i].c_str());
	}
	fprintf ( output, "  PDF \n");
	fprintf ( output, " $+1 \\sigma$ &");
	for( int i=0; i<sysTypeSize; i++){
		fprintf ( output, " $%+1.2f$ &", pSigma[i]);
	}
	fprintf ( output, " $< 1$\\\\ \n");
	fprintf ( output, " $-1 \\sigma$ &");
	for( int i=0; i<sysTypeSize; i++){
		fprintf ( output, " $%+1.2f$ &", mSigma[i]);
	}
	fprintf ( output, " $< 1$\\\\ \n");
}
