#include <iostream>
#include <fstream>
#include <sstream>
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
const int _up = 0;
const int _down = 1;

string double2str(double i){
	string s;
	stringstream ss(s);
	ss << i;
	return ss.str();
}
void overFlow( TH1D* h, double xmin, double xmax );
void clonePlotsOverflow( TH1D* h1, string outname, double xmin, double xmax );
void clonePlots( TH1D* h1, string outname );
void combineAllPlots( vector<TH1D*> inh, vector<double> weight, vector<string> name, string outname, TCanvas* c1, string canvName, bool showPlots=0, int index=0 );

void bprimeCombTemplate(){
	string loadpath = "results/ForScanBprimeBR_17Sep";
	string savepath = "results/ForScanBprimeBR_17Sep";

	string fileName="Template_HT.root";

	TFile* f_in = new TFile((loadpath+"/"+fileName).c_str());
	TFile* f_out = new TFile((savepath+"/Template_HT_bprimeComb.root").c_str(), "RECREATE");

	//enum CATECORY{ _nomal, _1b, _2b, _cat};
	enum CATECORY{ _1b, _2b, _cat};
	enum UNCNAME{ Stat, JER, JES, SFb, SFl, CA8, PU, TopPtReWrt, TTMatch, TTScale, TotalUnc};
	enum MASSPOINT{ M500, M600, M700, M800, M900, M1000, M1200, TotalMass};
	enum SAMPLE{  bHbH, bZbZ, tWtW, tWbH, tWbZ, bZbH, TotalSample};
	const int signalUnc = TotalUnc - 3;
	//string cate[_cat]={"CatAll", "Cat1b", "Cat2b"};
	string cate[_cat]={"Cat1b", "Cat2b"};
	string mass[TotalMass]={"500", "600", "700", "800", "900", "1000", "1200"};
	string sample[TotalSample]={"BHBH", "BZBZ", "TWTW", "BHTW", "BZTW", "BHBZ"};
	string uncsig[TotalUnc]={"SigStat", "JER", "JES", "SFb", "SFl", "CA8", "PU"};

	TH1D *_data[_cat], *_bkg[_cat][TotalUnc][2], *_bkgNomal[_cat], *_signal[_cat][TotalSample][TotalMass][signalUnc][2], *_signalNomal[_cat][TotalSample][TotalMass];
	cout<<"Clone data..."<<endl;
	//_data[_nomal]=(TH1D*)f_in->Get("CatAll_data_obs"); 	clonePlots(_data[_nomal], "CatAll_data_obs");
	_data[_1b]   =(TH1D*)f_in->Get("Cat1b_data_obs"); 	clonePlots(_data[_1b], "Cat1b_data_obs");
	_data[_2b]   =(TH1D*)f_in->Get("Cat2b_data_obs"); 	clonePlots(_data[_2b], "Cat2b_data_obs");

	cout<<"Clone background..."<<endl;
	for( int c=0; c<_cat; c++ ){
		_bkgNomal[c] =  (TH1D*)f_in->Get((cate[c]+"_background").c_str());
		_bkg[c][Stat][0] = (TH1D*)f_in->Get((cate[c]+"_background_StatUp").c_str());		
		_bkg[c][Stat][1] = (TH1D*)f_in->Get((cate[c]+"_background_StatDown").c_str());		
		_bkg[c][JER][0] = (TH1D*)f_in->Get((cate[c]+"_background_JERUp").c_str());		
		_bkg[c][JER][1] = (TH1D*)f_in->Get((cate[c]+"_background_JERDown").c_str());		
		_bkg[c][JES][0] = (TH1D*)f_in->Get((cate[c]+"_background_JESUp").c_str());		
		_bkg[c][JES][1] = (TH1D*)f_in->Get((cate[c]+"_background_JESDown").c_str());		
		_bkg[c][SFb][0] = (TH1D*)f_in->Get((cate[c]+"_background_SFbUp").c_str());		
		_bkg[c][SFb][1] = (TH1D*)f_in->Get((cate[c]+"_background_SFbDown").c_str());		
		_bkg[c][SFl][0] = (TH1D*)f_in->Get((cate[c]+"_background_SFlUp").c_str());	
		_bkg[c][SFl][1] = (TH1D*)f_in->Get((cate[c]+"_background_SFlDown").c_str());	
		_bkg[c][CA8][0] = (TH1D*)f_in->Get((cate[c]+"_background_CA8Up").c_str());	
		_bkg[c][CA8][1] = (TH1D*)f_in->Get((cate[c]+"_background_CA8Down").c_str());
		_bkg[c][PU][0] = (TH1D*)f_in->Get((cate[c]+"_background_PUUp").c_str());	
		_bkg[c][PU][1] = (TH1D*)f_in->Get((cate[c]+"_background_PUDown").c_str());
		_bkg[c][TopPtReWrt][0] = (TH1D*)f_in->Get((cate[c]+"_background_TopPtReWrtUp").c_str());
		_bkg[c][TopPtReWrt][1] = (TH1D*)f_in->Get((cate[c]+"_background_TopPtReWrtDown").c_str());	
		_bkg[c][TTMatch][0] = (TH1D*)f_in->Get((cate[c]+"_background_TTJetsMatchingUp").c_str());	
		_bkg[c][TTMatch][1] = (TH1D*)f_in->Get((cate[c]+"_background_TTJetsMatchingDown").c_str());	
		_bkg[c][TTScale][0] = (TH1D*)f_in->Get((cate[c]+"_background_TTJetsScaleUp").c_str());
		_bkg[c][TTScale][1] = (TH1D*)f_in->Get((cate[c]+"_background_TTJetsScaleDown").c_str());

		clonePlots(_bkgNomal[c], cate[c]+"_background");
		clonePlots(_bkg[c][Stat][0], cate[c]+"_background_StatUp");	
		clonePlots(_bkg[c][Stat][1], cate[c]+"_background_StatDown");	
		clonePlots(_bkg[c][JER][0], cate[c]+"_background_JERUp");	
		clonePlots(_bkg[c][JER][1], cate[c]+"_background_JERDown");	
		clonePlots(_bkg[c][JES][0], cate[c]+"_background_JESUp");	
		clonePlots(_bkg[c][JES][1], cate[c]+"_background_JESDown");	
		clonePlots(_bkg[c][SFb][0], cate[c]+"_background_SFbUp");	
		clonePlots(_bkg[c][SFb][1], cate[c]+"_background_SFbDown");	
		clonePlots(_bkg[c][SFl][0], cate[c]+"_background_SFlUp");	
		clonePlots(_bkg[c][SFl][1], cate[c]+"_background_SFlDown");	
		clonePlots(_bkg[c][CA8][0], cate[c]+"_background_CA8Up");	
		clonePlots(_bkg[c][CA8][1], cate[c]+"_background_CA8Down");	
		clonePlots(_bkg[c][PU][0], cate[c]+"_background_PUUp");	
		clonePlots(_bkg[c][PU][1], cate[c]+"_background_PUDown");
		clonePlots(_bkg[c][TopPtReWrt][0], cate[c]+"_background_TopPtReWrtUp");	
		clonePlots(_bkg[c][TopPtReWrt][1], cate[c]+"_background_TopPtReWrtDown");	
		clonePlots(_bkg[c][TTMatch][0], cate[c]+"_background_TTJetsMatchingUp");	
		clonePlots(_bkg[c][TTMatch][1], cate[c]+"_background_TTJetsMatchingDown");	
		clonePlots(_bkg[c][TTScale][0], cate[c]+"_background_TTJetsScaleUp");	
		clonePlots(_bkg[c][TTScale][1], cate[c]+"_background_TTJetsScaleDown");

		//Get signal MC histo
		for( int sig=0; sig<TotalSample; sig++){
			for( int m=0; m<TotalMass; m++){
				_signalNomal[c][sig][m] = (TH1D*)f_in->Get((cate[c]+"_"+sample[sig]+mass[m]).c_str());
				_signal[c][sig][m][Stat][0] = (TH1D*)f_in->Get((cate[c]+"_"+sample[sig]+mass[m]+"_SigStatUp").c_str());		
				_signal[c][sig][m][Stat][1] = (TH1D*)f_in->Get((cate[c]+"_"+sample[sig]+mass[m]+"_SigStatDown").c_str());		
				_signal[c][sig][m][JER][0] = (TH1D*)f_in->Get((cate[c]+"_"+sample[sig]+mass[m]+"_JERUp").c_str());		
				_signal[c][sig][m][JER][1] = (TH1D*)f_in->Get((cate[c]+"_"+sample[sig]+mass[m]+"_JERDown").c_str());		
				_signal[c][sig][m][JES][0] = (TH1D*)f_in->Get((cate[c]+"_"+sample[sig]+mass[m]+"_JESUp").c_str());		
				_signal[c][sig][m][JES][1] = (TH1D*)f_in->Get((cate[c]+"_"+sample[sig]+mass[m]+"_JESDown").c_str());		
				_signal[c][sig][m][SFb][0] = (TH1D*)f_in->Get((cate[c]+"_"+sample[sig]+mass[m]+"_SFbUp").c_str());		
				_signal[c][sig][m][SFb][1] = (TH1D*)f_in->Get((cate[c]+"_"+sample[sig]+mass[m]+"_SFbDown").c_str());		
				_signal[c][sig][m][SFl][0] = (TH1D*)f_in->Get((cate[c]+"_"+sample[sig]+mass[m]+"_SFlUp").c_str());	
				_signal[c][sig][m][SFl][1] = (TH1D*)f_in->Get((cate[c]+"_"+sample[sig]+mass[m]+"_SFlDown").c_str());	
				_signal[c][sig][m][CA8][0] = (TH1D*)f_in->Get((cate[c]+"_"+sample[sig]+mass[m]+"_CA8Up").c_str());	
				_signal[c][sig][m][CA8][1] = (TH1D*)f_in->Get((cate[c]+"_"+sample[sig]+mass[m]+"_CA8Down").c_str());
				_signal[c][sig][m][PU][0] = (TH1D*)f_in->Get((cate[c]+"_"+sample[sig]+mass[m]+"_PUUp").c_str());	
				_signal[c][sig][m][PU][1] = (TH1D*)f_in->Get((cate[c]+"_"+sample[sig]+mass[m]+"_PUDown").c_str());
				/*clonePlots(_signalNomal[c][sig][m], cate[c]+"_"+sample[sig]+mass[m]);
				clonePlots(_signal[c][sig][m][JER][0], cate[c]+"_"+sample[sig]+mass[m]+"_JERUp");	
				clonePlots(_signal[c][sig][m][JER][1], cate[c]+"_"+sample[sig]+mass[m]+"_JERDown");	
				clonePlots(_signal[c][sig][m][JES][0], cate[c]+"_"+sample[sig]+mass[m]+"_JESUp");	
				clonePlots(_signal[c][sig][m][JES][1], cate[c]+"_"+sample[sig]+mass[m]+"_JESDown");	
				clonePlots(_signal[c][sig][m][SFb][0], cate[c]+"_"+sample[sig]+mass[m]+"_SFbUp");	
				clonePlots(_signal[c][sig][m][SFb][1], cate[c]+"_"+sample[sig]+mass[m]+"_SFbDown");	
				clonePlots(_signal[c][sig][m][SFl][0], cate[c]+"_"+sample[sig]+mass[m]+"_SFlUp");	
				clonePlots(_signal[c][sig][m][SFl][1], cate[c]+"_"+sample[sig]+mass[m]+"_SFlDown");	
				clonePlots(_signal[c][sig][m][CA8][0], cate[c]+"_"+sample[sig]+mass[m]+"_CA8Up");	
				clonePlots(_signal[c][sig][m][CA8][1], cate[c]+"_"+sample[sig]+mass[m]+"_CA8Down");	
				clonePlots(_signal[c][sig][m][PU][0], cate[c]+"_"+sample[sig]+mass[m]+"_PUUp");	
				clonePlots(_signal[c][sig][m][PU][1], cate[c]+"_"+sample[sig]+mass[m]+"_PUDown");*/
			}
		}
	}

	TCanvas* c1 = new TCanvas("c1", "", 425, 350);
	//Caculate weight for each samples
	//| BHBH |	| 1	0	0	-0.5	0	-0.5	|| BR_bH*BR_bH |
	//| BZBZ |	| 0	1	0	0	-0.5	-0.5	|| BR_bZ*BR_bZ |
	//| TWTW |  _ 	| 0	0	1	-0.5	-0.5	0	|| BR_tW*BR_tW |
	//| TWBH |  _	| 0	0	0	2	0	0	||2*BR_tW*BR_bH|
	//| TWBZ |	| 0	0	0	0	2	0	||2*BR_tW*BR_bZ|
	//| BZBH |	| 0	0	0	0	0	2	||2*BR_bZ*BR_bH|
	int totalWeight(1);
	cout<<"\tB(bH)\tB(bZ)\tB(tW)\t|\tbHbH\tbZbZ\ttWtW\ttWbH\ttWbZ\tbZbH"<<endl;
	for( double BR_bH=0.; BR_bH<=1.; BR_bH+=0.1 ){
		double BR_bZ(0.);
		while( BR_bZ <= 1.-BR_bH){
			double BR_tW = 1. - BR_bZ - BR_bH;  if(BR_tW<0.00001) BR_tW=0.; 
			double matrix[6][6] = {{1, 0, 0, -0.5, 0, -0.5}, {0,1,0,0,-0.5, -0.5}, {0,0,1,-0.5,-0.5,0}, {0,0,0,2,0,0},{0,0,0,0,2,0}, {0,0,0,0,0,2}};
			double prob[TotalSample]={BR_bH*BR_bH, BR_bZ*BR_bZ, BR_tW*BR_tW, 2*BR_tW*BR_bH, 2*BR_tW*BR_bZ, 2*BR_bZ*BR_bH};
			double weight[TotalSample]={0., 0., 0., 0., 0., 0.};
			vector<double> weightV;
			vector<string> nameV;
			for( int out=0; out<6; out++){
				for( int s=0; s<TotalSample; s++){
					weight[out]+=matrix[out][s]*prob[s]; if(abs(weight[out]) < 0.00001) weight[out]=0.;
				}
				weightV.push_back(weight[out]);
				nameV.push_back(sample[out]);
			}
			cout<<totalWeight<<"\t"<<BR_bH<<"\t"<<BR_bZ<<"\t"<<BR_tW<<"\t|\t"<<weight[bHbH]<<"\t"<<weight[bZbZ]<<"\t"<<weight[tWtW]<<"\t"<<weight[tWbH]<<"\t"<<weight[tWbZ]<<"\t"<<weight[bZbH]<<endl;
			for( int c=0; c<_cat; c++){
				vector<TH1D*> hV;
				for( int m=0; m<TotalMass; m++){
					vector<TH1D*> hV;
					hV.push_back(_signalNomal[c][bHbH][m]); 
					hV.push_back(_signalNomal[c][bZbZ][m]); 
					hV.push_back(_signalNomal[c][tWtW][m]); 
					hV.push_back(_signalNomal[c][tWbH][m]); 
					hV.push_back(_signalNomal[c][tWbZ][m]); 
					hV.push_back(_signalNomal[c][bZbH][m]);
					string BR_bHs, BR_bZs, BR_tWs;
					if( BR_bH == 0. ) BR_bHs="00"; else if( BR_bH <= 0.91 ) BR_bHs="0"+double2str(BR_bH*10); else BR_bHs=double2str(BR_bH*10); 
					if( BR_bZ == 0. ) BR_bZs="00"; else if( BR_bZ <= 0.91 ) BR_bZs="0"+double2str(BR_bZ*10); else BR_bZs=double2str(BR_bZ*10); 
					if( BR_tW == 0. ) BR_tWs="00"; else if( BR_tW <= 0.91 ) BR_tWs="0"+double2str(BR_tW*10); else BR_tWs=double2str(BR_tW*10); 
					combineAllPlots(hV, weightV, nameV, cate[c]+"_BpM"+mass[m]+"_bH"+BR_bHs+"_bZ"+BR_bZs+"_tW"+BR_tWs, c1, "temp.png", 0);
					for( int sys=0; sys<signalUnc; sys++){
						vector<TH1D*> hVup, hVdown;
						hVup.push_back(_signal[c][bHbH][m][sys][_up]); hVdown.push_back(_signal[c][bHbH][m][sys][_down]);
						hVup.push_back(_signal[c][bZbZ][m][sys][_up]); hVdown.push_back(_signal[c][bZbZ][m][sys][_down]);
						hVup.push_back(_signal[c][tWtW][m][sys][_up]); hVdown.push_back(_signal[c][tWtW][m][sys][_down]);
						hVup.push_back(_signal[c][tWbH][m][sys][_up]); hVdown.push_back(_signal[c][tWbH][m][sys][_down]);
						hVup.push_back(_signal[c][tWbZ][m][sys][_up]); hVdown.push_back(_signal[c][tWbZ][m][sys][_down]);
						hVup.push_back(_signal[c][bZbH][m][sys][_up]); hVdown.push_back(_signal[c][bZbH][m][sys][_down]);
						combineAllPlots(hVup, weightV, nameV, cate[c]+"_BpM"+mass[m]+"_bH"+BR_bHs+"_bZ"+BR_bZs+"_tW"+BR_tWs+"_"+uncsig[sys]+"Up", c1, "temp.png", 0, 1);
						combineAllPlots(hVdown, weightV, nameV, cate[c]+"_BpM"+mass[m]+"_bH"+BR_bHs+"_bZ"+BR_bZs+"_tW"+BR_tWs+"_"+uncsig[sys]+"Down", c1, "temp.png", 0, 2);
					}
				}
			}
			BR_bZ+=0.1;
			totalWeight++;
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
void clonePlotsOverflow( TH1D* h1, string outname, double xmin, double xmax ){
	TH1D* h = (TH1D*)h1->Clone(outname.c_str());
	overFlow(h, xmin, xmax);
	h->Integral();
}
void clonePlots( TH1D* h1, string outname){
	TH1D* h = (TH1D*)h1->Clone(outname.c_str());
	h->Integral();
}
void combineAllPlots( vector<TH1D*> inh, vector<double> weight, vector<string> name, string outname, TCanvas* c1, string canvName, bool showPlots, int index ){
	if( inh.size() != weight.size() ){
		cout<<"WARNING: Size of h must the same with weight!";
		return;
	}
	const int size=inh.size();
	double max(0);
	TH1D* h[size], *hSum, *hSumtmp;
	TLegend* l = new TLegend(0.53, 0.77, 0.92, 0.92);
	l->SetFillStyle(0);
	l->SetBorderSize(0);
	for( int i=0; i<size; i++){ 
		char text[100];
		sprintf(text, "index_%d%d", i, index);
		h[i]=(TH1D*)inh[i]->Clone(text);
		if( h[i]->GetMaximum() > max ) max=h[i]->GetMaximum();
		h[i]->Scale(weight[i]);
		if( i==0 ) hSum=(TH1D*)h[i]->Clone(outname.c_str());
		else hSum->Add(h[i]);
		h[i]->SetLineColor(2+i);
		h[i]->SetLineWidth(2);
		l->AddEntry(h[i], (name[i]+" x "+double2str(weight[i])).c_str(), "l");
	}
	//TCanvas* c1 = new TCanvas("c1", "", 425, 350);
	char text[100];
	sprintf(text, "tmp_%d", index);
	hSumtmp = (TH1D*)hSum->Clone(text);
	hSumtmp->SetMaximum(max+max*10);
	hSumtmp->SetMinimum(-(max+max*10));
	hSumtmp->SetTitle(outname.c_str());
	hSumtmp->SetFillColor(16);
	l->AddEntry(hSumtmp, "Sum all", "f");
	if( showPlots ){
		hSumtmp->Draw();
		for( int i=0; i<size; i++) h[i]->Draw("SAME");
		c1->SaveAs(canvName.c_str());
	}
	for( int i=0; i<size; i++ ) delete h[i];
	//delete *h;
	delete hSumtmp;
}
