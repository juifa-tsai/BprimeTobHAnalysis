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
#include "CMSstyle.C"
using namespace std;
const int A = 0;
const int B = 1;
const int C = 2;
const int D = 3;
const int _cate_ = 3;
const int all = 0;
const int b1 = 1;
const int b2 = 2;
const bool blind = 1;
const bool unblind = 0;

void setXYTitle(TH1D* h, string xt, string yt, string xu="");
void draw( string fs, TCanvas* c1, string save, int cate, string hName, string xTitle, string yTitle, string xUnit, double xmin, double xmax, bool isBlind=1);
void FinalShapePlots(){

	CMSstyle();

	string Path="results/bin400GeV";

	vector<string> hName, xTitle, yTitle, xUnit; vector<double> xmin, xmax;
	hName.push_back("HT"); xTitle.push_back("H_{T}"); yTitle.push_back("Events"); xUnit.push_back("GeV"); xmin.push_back(800); xmax.push_back(1800);

	TCanvas* c1 = new TCanvas("c1", "", 425, 350);
	const int hNameSize = hName.size();
	for( int h=0; h<hNameSize; h++){
		cout<<"Ploting "<<hName[h]<<endl;
		string f = Path+"/ABCDResultTemplate_ABCDana_"+hName[h]+".root";
		draw( f, c1, Path, all, hName[h], xTitle[h], yTitle[h], xUnit[h], xmin[h], xmax[h], unblind );	
		draw( f, c1, Path, b1,  hName[h], xTitle[h], yTitle[h], xUnit[h], xmin[h], xmax[h], unblind );	
		draw( f, c1, Path, b2,  hName[h], xTitle[h], yTitle[h], xUnit[h], xmin[h], xmax[h], unblind );
	}
}
//CatAll_data_obs
void draw( string fs, TCanvas* c1, string save, int cate, string hName, string xTitle, string yTitle, string xUnit, double xmin, double xmax, bool isBlind ){
	c1->Clear();
	TFile* f = new TFile(fs.c_str());
	TH1D *h_data, *h_bkg, *hEr_bkg, *hsig500, *hsig800;
	string savename;
	if( cate == all ){
		h_data = (TH1D*)f->Get("CatAll_data_obs");
		h_bkg  = (TH1D*)f->Get("CatAll_background");
		hsig500 = (TH1D*)f->Get("CatAll_BHBH500");
		hsig800 = (TH1D*)f->Get("CatAll_BHBH800");
		savename = hName;
	}else if( cate == b1 ){
		h_data = (TH1D*)f->Get("Cat1b_data_obs");
		h_bkg  = (TH1D*)f->Get("Cat1b_background");
		hsig500 = (TH1D*)f->Get("Cat1b_BHBH500");
		hsig800 = (TH1D*)f->Get("Cat1b_BHBH800");
		savename = hName+"_1bjet";
	}else if( cate == b2 ){
		h_data = (TH1D*)f->Get("Cat2b_data_obs");
		h_bkg  = (TH1D*)f->Get("Cat2b_background");
		hsig500 = (TH1D*)f->Get("Cat2b_BHBH500");
		hsig800 = (TH1D*)f->Get("Cat2b_BHBH800");
		savename = hName+"_over1bjet";
	}else{
		cout<<"Error: Here is not exist this category: "<<cate<<endl;
		return;
	}
	hEr_bkg = (TH1D*)h_bkg->Clone("error");

	double max(0);	
	if( h_bkg->GetMaximum() >= hsig500->GetMaximum() ){
		if( h_bkg->GetMaximum() >= hsig800->GetMaximum()){	
			if( h_bkg->GetMaximum() >= h_data->GetMaximum())
				max = h_bkg->GetMaximum() + h_bkg->GetBinError(h_bkg->GetMaximumBin());
			else
				max = h_data->GetMaximum() + h_data->GetBinError(h_data->GetMaximumBin());
		}else{
			if( hsig800->GetMaximum() >= h_data->GetMaximum())
				max = hsig800->GetMaximum();
			else
				max = h_data->GetMaximum() + h_data->GetBinError(h_data->GetMaximumBin());
		}
	}else{
		if( hsig500->GetMaximum() >= hsig800->GetMaximum()){	
			if( hsig500->GetMaximum() >= h_data->GetMaximum())
				max = hsig500->GetMaximum();
			else
				max = h_data->GetMaximum() + h_data->GetBinError(h_data->GetMaximumBin());
		}else{
			if( hsig800->GetMaximum() >= h_data->GetMaximum())
				max = hsig800->GetMaximum();
			else
				max = h_data->GetMaximum() + h_data->GetBinError(h_data->GetMaximumBin());
		}
	}
	hEr_bkg->GetXaxis()->SetRangeUser(xmin, xmax);
	setXYTitle(hEr_bkg, xTitle, yTitle, xUnit);
	TAxis* ax = hEr_bkg->GetXaxis() ; 
	TAxis* ay = hEr_bkg->GetYaxis() ; 
	beautifyAxis(ax); 
	beautifyAxis(ay);
	ax->SetTitleOffset(0.8); 
	h_bkg->SetLineStyle(1);
	h_bkg->SetLineWidth(3);
	h_bkg->SetLineColor(9);

	h_data->SetMarkerStyle(8);
	h_data->SetMarkerSize(1);
	h_data->SetMarkerColor(1);
	h_data->SetLineColor(1);
	h_data->SetLineWidth(3);

	hsig500->SetLineStyle(5);
	hsig500->SetLineWidth(4);
	hsig500->SetLineColor(2);
	hsig800->SetLineStyle(5);
	hsig800->SetLineWidth(4);
	hsig800->SetLineColor(3);

	hEr_bkg->SetLineColor(46);
	hEr_bkg->SetFillColor(46);
	hEr_bkg->SetFillStyle(3345);	
	hEr_bkg->SetMarkerStyle(0);
	hEr_bkg->SetMarkerSize(0);

	TLegend* leg = new TLegend(0.52, 0.75, 0.93, 0.92);
	leg->SetBorderSize(1);
	leg->SetTextFont(132);
	leg->SetLineColor(1);
	leg->SetLineStyle(1);
	leg->SetLineWidth(1);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.045);

	TPaveText *tlumi_; 
	tlumi_ = new TPaveText(0.2,0.97,0.93,0.985,"NDC");
	tlumi_->AddText("CMS Preliminary 2012, L = 19.7/fb, #sqrt{s} = 8TeV");
	tlumi_->SetTextColor(kBlack);
	tlumi_->SetFillColor(kWhite);
	tlumi_->SetFillStyle(0);
	tlumi_->SetBorderSize(0);
	tlumi_->SetTextAlign(12);
	tlumi_->SetTextSize(0.05);

	TPaveText *cate_ = new TPaveText(0.2,0.25,0.4, 0.35,"NDC"); ;
	if( cate == all ) cate_->AddText("1b and #geq 2b categories");
	else if( cate == b1 ) cate_->AddText("1b category");
	else if( cate == b2 ) cate_->AddText("#geq 2b category");
	cate_->SetTextColor(kBlack);
	cate_->SetFillColor(kWhite);
	cate_->SetFillStyle(0);
	cate_->SetBorderSize(0);
	cate_->SetTextAlign(12);
	cate_->SetTextSize(0.05);

 
	if( !isBlind ) leg->AddEntry(h_data, "Data" ,"le");
	leg->AddEntry(h_bkg, "Backgrond" ,"l"); 
	leg->AddEntry(hEr_bkg, "Bkg. estimation error", "f");
	leg->AddEntry(hsig500, "b'(500)", "l");
	leg->AddEntry(hsig800, "b'(800)", "l");

	c1->SetLogy(0);
	hEr_bkg->SetMaximum(max+max/10);
	hEr_bkg->Draw("E2");
	h_bkg->Draw("HIST SAME");
	hsig500->Draw("HIST SAME");
	hsig800->Draw("HIST SAME");
	if( !isBlind )h_data->Draw("E1 SAME");
	tlumi_->Draw();
	cate_->Draw();
	leg->Draw();
	if( isBlind ) c1->SaveAs((save+"/"+savename+"_BkgFullSyst_NoData_Linear.pdf").c_str());
	else c1->SaveAs((save+"/"+savename+"_BkgFullSyst_withData_Linear.pdf").c_str());

	c1->SetLogy(1);
	hEr_bkg->SetMaximum(max*15);
	hEr_bkg->Draw("E2");
	h_bkg->Draw("HIST SAME");
	hsig500->Draw("HIST SAME");
	hsig800->Draw("HIST SAME");
	if( !isBlind )h_data->Draw("E1 SAME");
	tlumi_->Draw();
	cate_->Draw();
	leg->Draw();
	if( isBlind ) c1->SaveAs((save+"/"+savename+"_BkgFullSyst_NoData_LOG.pdf").c_str());
	else c1->SaveAs((save+"/"+savename+"_BkgFullSyst_withData_LOG.pdf").c_str());
}

void setXYTitle(TH1D* h, string xt, string yt, string xu){
	float binwidth = 0;
	char  xtitle[100], ytitle[100];
	binwidth = h->GetBinWidth(1);
	if( xu.size()!=0 ){
		sprintf( xtitle, "%s [%s]", xt.c_str(), xu.c_str());
		if (binwidth<1){
			sprintf( ytitle, "%s / %0.2f%s", yt.c_str(), binwidth, xu.c_str());
		}else{
			sprintf( ytitle, "%s / %1.0f%s", yt.c_str(), binwidth, xu.c_str());
		}
	}else{
		sprintf( xtitle, "%s", xt.c_str());
		if (binwidth<1){
			sprintf( ytitle, "%s / %0.2f", yt.c_str(), binwidth);
		}else{
			sprintf( ytitle, "%s / %1.0f", yt.c_str(), binwidth);
		}
	}
	h->SetXTitle(xtitle);
	h->SetYTitle(ytitle);
}
