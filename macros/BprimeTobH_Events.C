#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TTree.h"
using namespace std;
void BprimeTobH_Events(){

	string loadpath = "../test/OnLxplus/31Aug_Data_StoreTree/Final_histograms_ABCD.root";
	string savepath = "results/bprimeTobH_Events/";

	fstream output[3];
	string saveName[3] = {"BprimeTobH_Events", "BprimeTobH_Events_1b", "BprimeTobH_Events_2b"};	

	enum cate{ _all, _1b, _2b };

	TFile* f = new TFile(loadpath.c_str());
	TTree* t =(TTree*)f->Get("ABCD/tree");
	int RunNo, LumiNo, EvtCat;
	Long64_t EvtNo;
	t->SetBranchAddress("EvtInfo.RunNo", &RunNo);
	t->SetBranchAddress("EvtInfo.EvtNo", &EvtNo);
	t->SetBranchAddress("EvtInfo.LumiNo", &LumiNo);
	t->SetBranchAddress("EvtInfo.EvtCat", &EvtCat);

	output[_all].open((savepath+saveName[_all]+".txt").c_str(),ios_base::out);
	output[_1b].open((savepath+saveName[_1b]+".txt").c_str(),ios_base::out);
	output[_2b].open((savepath+saveName[_2b]+".txt").c_str(),ios_base::out);
	output[_all]<<"#   RunNo\tLumiNo\tEvtNo"<<endl;
	output[_1b]<<"#   RunNo\tLumiNo\tEvtNo"<<endl;
	output[_2b]<<"#   RunNo\tLumiNo\tEvtNo"<<endl;

	for( int evt=0; evt<t->GetEntries(); evt++){
		t->GetEntry(evt);
		output[_all]<<"    "<<RunNo<<"\t"<<LumiNo<<"\t"<<EvtNo<<endl;
		if( EvtCat == 1 ) output[_1b]<<"    "<<RunNo<<"\t"<<LumiNo<<"\t"<<EvtNo<<endl;
		else if( EvtCat == 2 ) output[_2b]<<"    "<<RunNo<<"\t"<<LumiNo<<"\t"<<EvtNo<<endl;
	}
}
