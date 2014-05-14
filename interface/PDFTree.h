#ifndef PDFTree_h
#define PDFTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class PDFTree {
  public :
    Int_t           McFlag;
    Int_t           PDFid1;
    Int_t           PDFid2;
    Float_t         PDFx1;
    Float_t         PDFx2;
    Float_t         PDFv1;
    Float_t         PDFv2;
    Float_t         qScale;
    Float_t         alphaQCD;
    Float_t         alphaQED;
    Double_t        evtwt;
    Double_t        puwt;

    void RegisterTree(TTree *tree) {
      tree->Branch("McFlag",   &McFlag,    "McFlag/I"   ) ; 
      tree->Branch("PDFid1",   &PDFid1,    "PDFid1/I"   ) ; 
      tree->Branch("PDFid2",   &PDFid2,    "PDFid2/I"   ) ; 
      tree->Branch("PDFx1",    &PDFx1,     "PDFx1/F"    ) ; 
      tree->Branch("PDFx2",    &PDFx2,     "PDFx2/F"    ) ; 
      tree->Branch("PDFv1",    &PDFv1,     "PDFv1/F"    ) ; 
      tree->Branch("PDFv2",    &PDFv2,     "PDFv2/F"    ) ; 
      tree->Branch("qScale",   &qScale,    "qScale/F"   ) ; 
      tree->Branch("alphaQCD", &alphaQCD,  "alphaQCD/F" ) ; 
      tree->Branch("alphaQED", &alphaQED,  "alphaQED/F" ) ; 
      tree->Branch("evtwt",    &evtwt,     "evtwt/D"    ) ; 
      tree->Branch("puwt",     &puwt,      "puwt/D"     ) ; 
    }

    void Register(TTree *tree) {
      tree->SetBranchAddress("McFlag", &McFlag, &b_EvtInfo_McFlag);
      tree->SetBranchAddress("PDFid1", &PDFid1, &b_EvtInfo_PDFid1);
      tree->SetBranchAddress("PDFid2", &PDFid2, &b_EvtInfo_PDFid2);
      tree->SetBranchAddress("PDFx1", &PDFx1, &b_EvtInfo_PDFx1);
      tree->SetBranchAddress("PDFx2", &PDFx2, &b_EvtInfo_PDFx2);
      tree->SetBranchAddress("PDFv1", &PDFv1, &b_EvtInfo_PDFv1);
      tree->SetBranchAddress("PDFv2", &PDFv2, &b_EvtInfo_PDFv2);
      tree->SetBranchAddress("qScale", &qScale, &b_EvtInfo_qScale);
      tree->SetBranchAddress("alphaQCD", &alphaQCD, &b_EvtInfo_alphaQCD);
      tree->SetBranchAddress("alphaQED", &alphaQED, &b_EvtInfo_alphaQED);
      tree->SetBranchAddress("evtwt", &evtwt, &b_evtwt);
      tree->SetBranchAddress("puwt", &puwt, &b_puwt);
    }

    TBranch        *b_EvtInfo_McFlag;   
    TBranch        *b_EvtInfo_PDFid1;   
    TBranch        *b_EvtInfo_PDFid2;   
    TBranch        *b_EvtInfo_PDFx1;    
    TBranch        *b_EvtInfo_PDFx2;    
    TBranch        *b_EvtInfo_PDFv1;    
    TBranch        *b_EvtInfo_PDFv2;    
    TBranch        *b_EvtInfo_qScale;   
    TBranch        *b_EvtInfo_alphaQCD; 
    TBranch        *b_EvtInfo_alphaQED; 
    TBranch        *b_evtwt;  
    TBranch        *b_puwt;   
};

#endif

