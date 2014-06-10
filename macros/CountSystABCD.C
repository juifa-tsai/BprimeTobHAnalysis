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
void printLetex( FILE* output, string hname, vector<string> systType, string loadpath, string filename[][3], int region);
void countAndStore( string title, string mean_h, string mean_f, string up_h, string up_f, string down_h, string down_f, fstream& f );
void _CountSystABCD(){

	string loadpath = "../test/OnLxplus/";
	string savepath = "results/";
	
	vector<string> systType;
	systType.push_back("JER");
	systType.push_back("JES");
	systType.push_back("SFb");
	systType.push_back("SFl");
	systType.push_back("SFHb");
	systType.push_back("SFHl");
	systType.push_back("PU");
	const int systTypeSize = systType.size();

	string fileName[systTypeSize][3];	
	fileName[0][Mean] = "07June_HT950_Mean_Sig/Final_histograms_ABCD.root";
	fileName[0][Down] = "07June_HT950_JERDown/Final_histograms_ABCD.root";
	fileName[0][Up]   = "07June_HT950_JERUp/Final_histograms_ABCD.root";
	fileName[1][Mean] = "07June_HT950_Mean_Sig/Final_histograms_ABCD.root";
	fileName[1][Down] = "07June_HT950_JESDown/Final_histograms_ABCD.root";
	fileName[1][Up]   = "07June_HT950_JESUp/Final_histograms_ABCD.root";
	fileName[2][Mean] = "07June_HT950_Mean_Sig/Final_histograms_ABCD.root";
	fileName[2][Down] = "07June_HT950_SFbDown/Final_histograms_ABCD.root";
	fileName[2][Up]   = "07June_HT950_SFbUp/Final_histograms_ABCD.root";
	fileName[3][Mean] = "07June_HT950_Mean_Sig/Final_histograms_ABCD.root";
	fileName[3][Down] = "07June_HT950_SFlDown/Final_histograms_ABCD.root";
	fileName[3][Up]   = "07June_HT950_SFlUp/Final_histograms_ABCD.root";
	fileName[4][Mean] = "07June_HT950_Mean_Sig/Final_histograms_ABCD.root";
	fileName[4][Down] = "07June_HT950_SFHbDown/Final_histograms_ABCD.root";
	fileName[4][Up]   = "07June_HT950_SFHbUp/Final_histograms_ABCD.root";
	fileName[5][Mean] = "07June_HT950_Mean_Sig/Final_histograms_ABCD.root";
	fileName[5][Down] = "07June_HT950_SFHlDown/Final_histograms_ABCD.root";
	fileName[5][Up]   = "07June_HT950_SFHlUp/Final_histograms_ABCD.root";
	fileName[6][Mean] = "07June_HT950_Mean_Sig/Final_histograms_ABCD.root";
	fileName[6][Down] = "07June_HT950_PUDown/Final_histograms_ABCD.root";
	fileName[6][Up]   = "07June_HT950_PUUp/Final_histograms_ABCD.root";

	vector<string> hName;
	hName.push_back("ABCDana_CutRegion");
	hName.push_back("ABCDana_CutRegion_1b");
	hName.push_back("ABCDana_CutRegion_2b");

	vector<string> Title, fullName; 
	//Title.push_back("BHBH450"); 	fullName.push_back("BprimeBprimeToBHBHinc_M-450");
	Title.push_back("BHBH500"); 	fullName.push_back("BprimeBprimeToBHBHinc_M-500");
	//Title.push_back("BHBH550"); 	fullName.push_back("BprimeBprimeToBHBHinc_M-550");
	Title.push_back("BHBH600"); 	fullName.push_back("BprimeBprimeToBHBHinc_M-600");
	//Title.push_back("BHBH650");	fullName.push_back("BprimeBprimeToBHBHinc_M-650");
	Title.push_back("BHBH700"); 	fullName.push_back("BprimeBprimeToBHBHinc_M-700");
	//Title.push_back("BHBH750"); 	fullName.push_back("BprimeBprimeToBHBHinc_M-750");
	Title.push_back("BHBH800"); 	fullName.push_back("BprimeBprimeToBHBHinc_M-800");
	Title.push_back("BHBH1000"); 	fullName.push_back("BprimeBprimeToBHBHinc_M-1000");
	Title.push_back("BHBH1200"); 	fullName.push_back("BprimeBprimeToBHBHinc_M-1200");
	//Title.push_back("BHBH1500"); 	fullName.push_back("BprimeBprimeToBHBHinc_M-1500");

	const int hSize = hName.size();
	const int titleSize = Title.size();
	
	fstream output[systTypeSize][hSize];
	FILE* pFile[hSize];

	for( int i=0; i<hSize; i++ ){
		cout<<"Print letex format file for "<<hName[i]<<"..."<<endl;
		pFile[i] = fopen((savepath+"SystUncertainties_"+hName[i]+".tex").c_str(),"w");
		for( int k=0; k<titleSize; k++){
			string hname = fullName[k]+"__"+hName[i];
			fprintf ( pFile[i], "=============== %s =============== \n",Title[k].c_str());	
			printLetex( pFile[i], hname, systType, loadpath, fileName, B_);
			fprintf ( pFile[i], "\n");	
			fprintf ( pFile[i], "\n");	
		}
		fclose(pFile[i]);	
	}
	for( int i=0; i<systTypeSize; i++ ){
		cout<<"[Uncertainty] "<<systType[i]<<endl;
		for( int j=0; j<hSize; j++ ){
			cout<<hName[j]<<"..."<<endl;
			string save = savepath + systType[i]+"_"+hName[j]+".txt";	
			output[i][j].open(save.c_str(),ios_base::out);
			for( int k=0; k<titleSize; k++){
				string hname = fullName[k]+"__"+hName[j];
				string fmean = loadpath + fileName[i][Mean]; 
				string fdown = loadpath + fileName[i][Down]; 
				string fup   = loadpath + fileName[i][Up];
				countAndStore( Title[k], hname, fmean, hname, fup, hname, fdown, output[i][j]);
			}
		}
	}

}
/* $+1 \sigma$ & $-0.68$ & $+3.27$ & $+1.42$  & $+0.01$ & $+5.04$ & $+0.07$ & $-0.19$& $< 1$\\ */
void printLetex( FILE* output, string hname, vector<string> systType, string loadpath, string filename[][3], int region){
	const int sysTypeSize=systType.size();
	double pSigma[sysTypeSize], mSigma[sysTypeSize];
	TFile* f_up[sysTypeSize];
	TFile* f_mean[sysTypeSize];
	TFile* f_down[sysTypeSize];

	fprintf ( output, "    sigma    ");
	for( int i=0; i<sysTypeSize; i++){
		f_up[i] = new TFile(   (loadpath + filename[i][Up]).c_str() );
		f_mean[i] = new TFile( (loadpath + filename[i][Mean]).c_str() );
		f_down[i] = new TFile( (loadpath + filename[i][Down]).c_str() );
		TH1D*  h_up   = (TH1D*)f_up[i]->Get(hname.c_str());
		TH1D*  h_mean = (TH1D*)f_mean[i]->Get(hname.c_str());
		TH1D*  h_down = (TH1D*)f_down[i]->Get(hname.c_str());
		pSigma[i] = ( h_up->GetBinContent(region)   - h_mean->GetBinContent(region) )/ h_mean->GetBinContent(region) * 100;
		mSigma[i] = ( h_down->GetBinContent(region) - h_mean->GetBinContent(region) )/ h_mean->GetBinContent(region) * 100;
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

void countAndStore( string title, string mean_h, string mean_f, string up_h, string up_f, string down_h, string down_f, fstream& f ){

	TFile* f_up   = new TFile(up_f.c_str());
	TFile* f_mean = new TFile(mean_f.c_str());
	TFile* f_down = new TFile(down_f.c_str());

	TH1D*  h_up   = (TH1D*)f_up->Get(up_h.c_str());
	TH1D*  h_mean = (TH1D*)f_mean->Get(mean_h.c_str());
	TH1D*  h_down = (TH1D*)f_down->Get(down_h.c_str());

	float A, B, C, D, Au, Bu, Cu, Du, Ad, Bd, Cd, Dd;
	A = h_mean->GetBinContent(1);
	B = h_mean->GetBinContent(2);
	C = h_mean->GetBinContent(3);
	D = h_mean->GetBinContent(4);
	Au = h_up->GetBinContent(1);
	Bu = h_up->GetBinContent(2);
	Cu = h_up->GetBinContent(3);
	Du = h_up->GetBinContent(4);
	Ad = h_down->GetBinContent(1);
	Bd = h_down->GetBinContent(2);
	Cd = h_down->GetBinContent(3);
	Dd = h_down->GetBinContent(4);

	
	float ASystUp = 100*(Au-A)/A; 
	float BSystUp = 100*(Bu-B)/B; 
	float CSystUp = 100*(Cu-C)/C; 
	float DSystUp = 100*(Du-D)/D; 
	float ASystDown = 100*(Ad-A)/A; 
	float BSystDown = 100*(Bd-B)/B; 
	float CSystDown = 100*(Cd-C)/C; 
	float DSystDown = 100*(Dd-D)/D; 
	
	f<<"****************************************************"<<endl;
	f<<"****************** "<<title<<" ******************"<<endl;
	f<<"****************************************************"<<endl;
	f<<"Region\t"<<"Num\t\t"<<"Up[%]\t\t"<<"Down[%]\t\t"<<endl;
	f<<"A\t"<<A<<"\t\t"<<ASystUp<<"\t\t"<<ASystDown<<endl;
	f<<"B\t"<<B<<"\t\t"<<BSystUp<<"\t\t"<<BSystDown<<endl;
	f<<"C\t"<<C<<"\t\t"<<CSystUp<<"\t\t"<<CSystDown<<endl;
	f<<"D\t"<<D<<"\t\t"<<DSystUp<<"\t\t"<<DSystDown<<endl;
	f<<endl;

}
