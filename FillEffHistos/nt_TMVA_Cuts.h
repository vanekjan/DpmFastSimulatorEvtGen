//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep  1 14:09:42 2016 by ROOT version 5.34/14
// from TTree nt/
// found on file: Dpm.toyMc.100M.root
//////////////////////////////////////////////////////////

#ifndef nt_TMVA_Cuts_h
#define nt_TMVA_Cuts_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>

using namespace std;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class nt_TMVA_Cuts {
  public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // Declaration of leaf types
    Float_t         cent;
    Float_t         vx;
    Float_t         vy;
    Float_t         vz;
    Float_t         vzIdx;
    Float_t         pid;
    Float_t         w;
    Float_t         m;
    Float_t         pt;
    Float_t         eta;
    Float_t         y;
    Float_t         phi;
    Float_t         v0x;
    Float_t         v0y;
    Float_t         v0z;
    Float_t         rM;
    Float_t         rPt;
    Float_t         rEta;
    Float_t         rY;
    Float_t         rPhi;
    Float_t         rV0x;
    Float_t         rV0y;
    Float_t         rV0z;
    Float_t         dcaDaughters;
    Float_t         dcaDaughters_helix;
    Float_t         decayLength;
    Float_t         decayLength_helix;
    Float_t         MCdecayLength;
    Float_t         dcaDpmToPv;
    Float_t         dcaDpmToPv_helix;
    Float_t         cosTheta;
    Float_t         cosTheta_helix;
    //		Float_t         angle12;
    //		Float_t         cosThetaStar;
    Float_t         kM;
    Float_t         kPt;
    Float_t         kEta;
    Float_t         kY;
    Float_t         kPhi;
    Float_t         kDca;
    Float_t         kRM;
    Float_t         kRPt;
    Float_t         kREta;
    Float_t         kRY;
    Float_t         kRPhi;
    Float_t         kRVx;
    Float_t         kRVy;
    Float_t         kRVz;
    Float_t         kRDca;
    Float_t         kRDca_helix;
    //		Float_t         kRSDca;
    Float_t         kRDcaXY;
    Float_t         kRDcaZ;
    Float_t         kRDcaXYtoSV;
    Float_t         kRDcaZtoSV;
    //		Float_t         kEtaIdx;
    //		Float_t         kPtIdx;
    Float_t         kTpc;
    Float_t         p1M;
    Float_t         p1Pt;
    Float_t         p1Eta;
    Float_t         p1Y;
    Float_t         p1Phi;
    Float_t         p1Dca;
    Float_t         p1RM;
    Float_t         p1RPt;
    Float_t         p1REta;
    Float_t         p1RY;
    Float_t         p1RPhi;
    Float_t         p1RVx;
    Float_t         p1RVy;
    Float_t         p1RVz;
    Float_t         p1RDca;
    Float_t         p1RDca_helix;
    //		Float_t         p1RSDca;
    Float_t         p1RDcaXY;
    Float_t         p1RDcaZ;
    Float_t         p1RDcaXYtoSV;
    Float_t         p1RDcaZtoSV;
    //		Float_t         p1EtaIdx;
    //		Float_t         p1PtIdx;
    Float_t         p1Tpc;
    Float_t         p2M;
    Float_t         p2Pt;
    Float_t         p2Eta;
    Float_t         p2Y;
    Float_t         p2Phi;
    Float_t         p2Dca;
    Float_t         p2RM;
    Float_t         p2RPt;
    Float_t         p2REta;
    Float_t         p2RY;
    Float_t         p2RPhi;
    Float_t         p2RVx;
    Float_t         p2RVy;
    Float_t         p2RVz;
    Float_t         p2RDca;
    Float_t         p2RDca_helix;
    //		Float_t         p2RSDca;
    Float_t         p2RDcaXY;
    Float_t         p2RDcaZ;
    Float_t         p2RDcaXYtoSV;
    Float_t         p2RDcaZtoSV;
    //		Float_t         p2EtaIdx;
    //		Float_t         p2PtIdx;
    Float_t         p2Tpc;
    Float_t         kHft;
    Float_t         p1Hft;
    Float_t         p2Hft;
    //		Float_t         kTof;
    //		Float_t         p1Tof;
    //		Float_t         p2Tof;
    //		Float_t         sA;
    //		Float_t         sB;
    //		Float_t         dcakp1;
    //		Float_t         dcap1p2;
    //		Float_t         dcap2k;
    Float_t         mdV0Max;
    Float_t         mdV0Max_helix;

    // List of branches
    TBranch        *b_cent;   //!
    TBranch        *b_vx;   //!
    TBranch        *b_vy;   //!
    TBranch        *b_vz;   //!
    TBranch        *b_vzIdx;   //!
    TBranch        *b_pid;   //!
    TBranch        *b_w;   //!
    TBranch        *b_m;   //!
    TBranch        *b_pt;   //!
    TBranch        *b_eta;   //!
    TBranch        *b_y;   //!
    TBranch        *b_phi;   //!
    TBranch        *b_v0x;   //!
    TBranch        *b_v0y;   //!
    TBranch        *b_v0z;   //!
    TBranch        *b_rM;   //!
    TBranch        *b_rPt;   //!
    TBranch        *b_rEta;   //!
    TBranch        *b_rY;   //!
    TBranch        *b_rPhi;   //!
    TBranch        *b_rV0x;   //!
    TBranch        *b_rV0y;   //!
    TBranch        *b_rV0z;   //!
    TBranch        *b_dcaDaughters;   //!
    TBranch        *b_dcaDaughters_helix;   //!
    TBranch        *b_decayLength;   //!
    TBranch        *b_decayLength_helix;   //!
    TBranch        *b_MCdecayLength;   //!
    TBranch        *b_dcaDpmToPv;   //!
    TBranch        *b_dcaDpmToPv_helix;   //!
    TBranch        *b_cosTheta;   //!
    TBranch        *b_cosTheta_helix;   //!
    //		TBranch        *b_angle12;   //!
    //		TBranch        *b_cosThetaStar;   //!
    TBranch        *b_kM;   //!
    TBranch        *b_kPt;   //!
    TBranch        *b_kEta;   //!
    TBranch        *b_kY;   //!
    TBranch        *b_kPhi;   //!
    TBranch        *b_kDca;   //!
    TBranch        *b_kRM;   //!
    TBranch        *b_kRPt;   //!
    TBranch        *b_kREta;   //!
    TBranch        *b_kRY;   //!
    TBranch        *b_kRPhi;   //!
    TBranch        *b_kRVx;   //!
    TBranch        *b_kRVy;   //!
    TBranch        *b_kRVz;   //!
    TBranch        *b_kRDca;   //!
    TBranch        *b_kRDca_helix;   //!
    //		TBranch        *b_kRSDca;   //!
    TBranch        *b_kRDcaXY;   //!
    TBranch        *b_kRDcaZ;   //!
    TBranch        *b_kRDcaXYtoSV;   //! 
    TBranch        *b_kRDcaZtoSV;   //!
    //		TBranch        *b_kEtaIdx;   //!
    //		TBranch        *b_kPtIdx;   //!
    TBranch        *b_kTpc;   //!
    TBranch        *b_p1M;   //!
    TBranch        *b_p1Pt;   //!
    TBranch        *b_p1Eta;   //!
    TBranch        *b_p1Y;   //!
    TBranch        *b_p1Phi;   //!
    TBranch        *b_p1Dca;   //!
    TBranch        *b_p1RM;   //!
    TBranch        *b_p1RPt;   //!
    TBranch        *b_p1REta;   //!
    TBranch        *b_p1RY;   //!
    TBranch        *b_p1RPhi;   //!
    TBranch        *b_p1RVx;   //!
    TBranch        *b_p1RVy;   //!
    TBranch        *b_p1RVz;   //!
    TBranch        *b_p1RDca;   //!
    TBranch        *b_p1RDca_helix;   //!
    //		TBranch        *b_p1RSDca;   //!
    TBranch        *b_p1RDcaXY;   //!
    TBranch        *b_p1RDcaZ;   //!
    TBranch        *b_p1RDcaXYtoSV;   //! 
    TBranch        *b_p1RDcaZtoSV;   //!
    //		TBranch        *b_p1EtaIdx;   //!
    //		TBranch        *b_p1PtIdx;   //!
    TBranch        *b_p1Tpc;   //!
    TBranch        *b_p2M;   //!
    TBranch        *b_p2Pt;   //!
    TBranch        *b_p2Eta;   //!
    TBranch        *b_p2Y;   //!
    TBranch        *b_p2Phi;   //!
    TBranch        *b_p2Dca;   //!
    TBranch        *b_p2RM;   //!
    TBranch        *b_p2RPt;   //!
    TBranch        *b_p2REta;   //!
    TBranch        *b_p2RY;   //!
    TBranch        *b_p2RPhi;   //!
    TBranch        *b_p2RVx;   //!
    TBranch        *b_p2RVy;   //!
    TBranch        *b_p2RVz;   //!
    TBranch        *b_p2RDca;   //!
    TBranch        *b_p2RDca_helix;   //!
    //		TBranch        *b_p2RSDca;   //!
    TBranch        *b_p2RDcaXY;   //!
    TBranch        *b_p2RDcaZ;   //!
    TBranch        *b_p2RDcaXYtoSV;   //! 
    TBranch        *b_p2RDcaZtoSV;   //!
    //		TBranch        *b_p2EtaIdx;   //!
    //		TBranch        *b_p2PtIdx;   //!
    TBranch        *b_p2Tpc;   //!
    TBranch        *b_kHft;   //!
    TBranch        *b_p1Hft;   //!
    TBranch        *b_p2Hft;   //!
    //		TBranch        *b_kTof;   //!
    //		TBranch        *b_p1Tof;   //!
    //		TBranch        *b_p2Tof;   //!
    //		TBranch        *b_sA;   //!
    //		TBranch        *b_sB;   //!
    //		TBranch        *b_dcakp1;   //!
    //		TBranch        *b_dcap1p2;   //!
    //		TBranch        *b_dcap2k;   //!
    TBranch        *b_mdV0Max;   //!
    TBranch        *b_mdV0Max_helix;   //!

    // conditions
    int con_cent;
    int con_cent_up;
    int con_cent_down;
    float con_v0z;
    float con_cosTheta;
    float con_dcaDaughters;
    float con_decayLength;
    float con_kDca;
    float con_pDca;
    float con_kRPt;
    float con_pRPt;
    float con_dca_daughters;
    float con_dV0Max;

    int ReadMode; //ReadMode: 0 - pre-cuts, 1 - analysis cust
    int pTspectrum; //pTspectrum: 0 - flat (use weights when fillign histos), 1 - Levy (don't use weights - input spectrum realistic)


    nt_TMVA_Cuts(TTree *tree=0);
    virtual ~nt_TMVA_Cuts();
    virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual void     Loop();
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);

    TString out_file_name;
    virtual void Set_out_file_name(TString c) {out_file_name = c;};

    // setting conditions
    virtual void Set_con_cent(int c) {con_cent = c;};
    virtual void Set_con_cent_up(int c) {con_cent_up = c;};
    virtual void Set_con_cent_down(int c) {con_cent_down = c;};
    virtual void Set_con_v0z(int c) {con_v0z = c;};
    virtual void Set_con_cosTheta(float c) {con_cosTheta = c;};
    virtual void Set_con_dcaDaughters(float c) {con_dcaDaughters = c;};
    virtual void Set_con_decayLength(float c) {con_decayLength = c;};
    virtual void Set_con_kDca(float c) {con_kDca = c;};
    virtual void Set_con_pDca(float c) {con_pDca = c;};
    virtual void Set_con_kRPt(float c) {con_kRPt = c;};
    virtual void Set_con_pRPt(float c) {con_pRPt = c;};
    virtual void Set_con_dV0Max(float c) {con_dV0Max = c;};
    //virtual void Set_con_();

    //set ReadMode
    virtual void Set_ReadMode(int Mode) {ReadMode = Mode;};
    virtual void Set_pTspectrum(int spectrum) {pTspectrum = spectrum;};
};

#endif

#ifdef nt_TMVA_Cuts_cxx
nt_TMVA_Cuts::nt_TMVA_Cuts(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Dpm.toyMc.100M.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("Dpm.toyMc.100M.root");
    }
    f->GetObject("nt",tree);

  }
  Init(tree);

  out_file_name = "test_new.root";

  con_cent = -9;
  con_cent_up = -9;
  con_cent_down = -9;
  con_v0z = -9;
  con_cosTheta = -9.;
  con_dcaDaughters = -9.;
  con_decayLength = -9.;
  con_kDca = -9.;
  con_pDca = -9.;
  con_kRPt = -9.;
  con_pRPt = -9.;
  con_dca_daughters = -9.;
  con_dV0Max = -9.;
}

nt_TMVA_Cuts::~nt_TMVA_Cuts()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t nt_TMVA_Cuts::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t nt_TMVA_Cuts::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void nt_TMVA_Cuts::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("cent", &cent, &b_cent);
  fChain->SetBranchAddress("vx", &vx, &b_vx);
  fChain->SetBranchAddress("vy", &vy, &b_vy);
  fChain->SetBranchAddress("vz", &vz, &b_vz);
  fChain->SetBranchAddress("vzIdx", &vzIdx, &b_vzIdx);
  fChain->SetBranchAddress("pid", &pid, &b_pid);
  fChain->SetBranchAddress("w", &w, &b_w);
  fChain->SetBranchAddress("m", &m, &b_m);
  fChain->SetBranchAddress("pt", &pt, &b_pt);
  fChain->SetBranchAddress("eta", &eta, &b_eta);
  fChain->SetBranchAddress("y", &y, &b_y);
  fChain->SetBranchAddress("phi", &phi, &b_phi);
  fChain->SetBranchAddress("v0x", &v0x, &b_v0x);
  fChain->SetBranchAddress("v0y", &v0y, &b_v0y);
  fChain->SetBranchAddress("v0z", &v0z, &b_v0z);
  fChain->SetBranchAddress("rM", &rM, &b_rM);
  fChain->SetBranchAddress("rPt", &rPt, &b_rPt);
  fChain->SetBranchAddress("rEta", &rEta, &b_rEta);
  fChain->SetBranchAddress("rY", &rY, &b_rY);
  fChain->SetBranchAddress("rPhi", &rPhi, &b_rPhi);
  fChain->SetBranchAddress("rV0x", &rV0x, &b_rV0x);
  fChain->SetBranchAddress("rV0y", &rV0y, &b_rV0y);
  fChain->SetBranchAddress("rV0z", &rV0z, &b_rV0z);
  fChain->SetBranchAddress("dcaDaughters", &dcaDaughters, &b_dcaDaughters);
  fChain->SetBranchAddress("dcaDaughters_helix", &dcaDaughters_helix, &b_dcaDaughters_helix);
  fChain->SetBranchAddress("decayLength", &decayLength, &b_decayLength);
  fChain->SetBranchAddress("decayLength_helix", &decayLength_helix, &b_decayLength_helix);
  fChain->SetBranchAddress("MCdecayLength", &MCdecayLength, &b_MCdecayLength);
  fChain->SetBranchAddress("dcaDpmToPv", &dcaDpmToPv, &b_dcaDpmToPv);
  fChain->SetBranchAddress("dcaDpmToPv_helix", &dcaDpmToPv_helix, &b_dcaDpmToPv_helix);
  fChain->SetBranchAddress("cosTheta", &cosTheta, &b_cosTheta);
  fChain->SetBranchAddress("cosTheta_helix", &cosTheta_helix, &b_cosTheta_helix);
  //	fChain->SetBranchAddress("angle12", &angle12, &b_angle12);
  //	fChain->SetBranchAddress("cosThetaStar", &cosThetaStar, &b_cosThetaStar);
  fChain->SetBranchAddress("kM", &kM, &b_kM);
  fChain->SetBranchAddress("kPt", &kPt, &b_kPt);
  fChain->SetBranchAddress("kEta", &kEta, &b_kEta);
  fChain->SetBranchAddress("kY", &kY, &b_kY);
  fChain->SetBranchAddress("kPhi", &kPhi, &b_kPhi);
  fChain->SetBranchAddress("kDca", &kDca, &b_kDca);
  fChain->SetBranchAddress("kRM", &kRM, &b_kRM);
  fChain->SetBranchAddress("kRPt", &kRPt, &b_kRPt);
  fChain->SetBranchAddress("kREta", &kREta, &b_kREta);
  fChain->SetBranchAddress("kRY", &kRY, &b_kRY);
  fChain->SetBranchAddress("kRPhi", &kRPhi, &b_kRPhi);
  fChain->SetBranchAddress("kRVx", &kRVx, &b_kRVx);
  fChain->SetBranchAddress("kRVy", &kRVy, &b_kRVy);
  fChain->SetBranchAddress("kRVz", &kRVz, &b_kRVz);
  fChain->SetBranchAddress("kRDca", &kRDca, &b_kRDca);
  fChain->SetBranchAddress("kRDca_helix", &kRDca_helix, &b_kRDca_helix);
  //	fChain->SetBranchAddress("kRSDca", &kRSDca, &b_kRSDca);
  //	fChain->SetBranchAddress("kRDcaXY", &kRDcaXY, &b_kRDcaXY);
  //	fChain->SetBranchAddress("kRDcaZ", &kRDcaZ, &b_kRDcaZ);
  //	fChain->SetBranchAddress("kEtaIdx", &kEtaIdx, &b_kEtaIdx);
  //	fChain->SetBranchAddress("kPtIdx", &kPtIdx, &b_kPtIdx);
  fChain->SetBranchAddress("kTpc", &kTpc, &b_kTpc);
  fChain->SetBranchAddress("p1M", &p1M, &b_p1M);
  fChain->SetBranchAddress("p1Pt", &p1Pt, &b_p1Pt);
  fChain->SetBranchAddress("p1Eta", &p1Eta, &b_p1Eta);
  fChain->SetBranchAddress("p1Y", &p1Y, &b_p1Y);
  fChain->SetBranchAddress("p1Phi", &p1Phi, &b_p1Phi);
  fChain->SetBranchAddress("p1Dca", &p1Dca, &b_p1Dca);
  fChain->SetBranchAddress("p1RM", &p1RM, &b_p1RM);
  fChain->SetBranchAddress("p1RPt", &p1RPt, &b_p1RPt);
  fChain->SetBranchAddress("p1REta", &p1REta, &b_p1REta);
  fChain->SetBranchAddress("p1RY", &p1RY, &b_p1RY);
  fChain->SetBranchAddress("p1RPhi", &p1RPhi, &b_p1RPhi);
  fChain->SetBranchAddress("p1RVx", &p1RVx, &b_p1RVx);
  fChain->SetBranchAddress("p1RVy", &p1RVy, &b_p1RVy);
  fChain->SetBranchAddress("p1RVz", &p1RVz, &b_p1RVz);
  fChain->SetBranchAddress("p1RDca", &p1RDca, &b_p1RDca);
  fChain->SetBranchAddress("p1RDca_helix", &p1RDca_helix, &b_p1RDca_helix);
  //	fChain->SetBranchAddress("p1RSDca", &p1RSDca, &b_p1RSDca);
  fChain->SetBranchAddress("p1RDcaXY", &p1RDcaXY, &b_p1RDcaXY);
  fChain->SetBranchAddress("p1RDcaZ", &p1RDcaZ, &b_p1RDcaZ);
  fChain->SetBranchAddress("p1RDcaXYtoSV", &p1RDcaXYtoSV, &b_p1RDcaXYtoSV);
  fChain->SetBranchAddress("p1RDcaZtoSV", &p1RDcaZtoSV, &b_p1RDcaZtoSV);
  //	fChain->SetBranchAddress("p1EtaIdx", &p1EtaIdx, &b_p1EtaIdx);
  //	fChain->SetBranchAddress("p1PtIdx", &p1PtIdx, &b_p1PtIdx);
  fChain->SetBranchAddress("p1Tpc", &p1Tpc, &b_p1Tpc);
  fChain->SetBranchAddress("p2M", &p2M, &b_p2M);
  fChain->SetBranchAddress("p2Pt", &p2Pt, &b_p2Pt);
  fChain->SetBranchAddress("p2Eta", &p2Eta, &b_p2Eta);
  fChain->SetBranchAddress("p2Y", &p2Y, &b_p2Y);
  fChain->SetBranchAddress("p2Phi", &p2Phi, &b_p2Phi);
  fChain->SetBranchAddress("p2Dca", &p2Dca, &b_p2Dca);
  fChain->SetBranchAddress("p2RM", &p2RM, &b_p2RM);
  fChain->SetBranchAddress("p2RPt", &p2RPt, &b_p2RPt);
  fChain->SetBranchAddress("p2REta", &p2REta, &b_p2REta);
  fChain->SetBranchAddress("p2RY", &p2RY, &b_p2RY);
  fChain->SetBranchAddress("p2RPhi", &p2RPhi, &b_p2RPhi);
  fChain->SetBranchAddress("p2RVx", &p2RVx, &b_p2RVx);
  fChain->SetBranchAddress("p2RVy", &p2RVy, &b_p2RVy);
  fChain->SetBranchAddress("p2RVz", &p2RVz, &b_p2RVz);
  fChain->SetBranchAddress("p2RDca", &p2RDca, &b_p2RDca);
  fChain->SetBranchAddress("p2RDca_helix", &p2RDca_helix, &b_p2RDca_helix);
  //	fChain->SetBranchAddress("p2RSDca", &p2RSDca, &b_p2RSDca);
  fChain->SetBranchAddress("p2RDcaXY", &p2RDcaXY, &b_p2RDcaXY);
  fChain->SetBranchAddress("p2RDcaZ", &p2RDcaZ, &b_p2RDcaZ);
  fChain->SetBranchAddress("p2RDcaXYtoSV", &p2RDcaXYtoSV, &b_p2RDcaXYtoSV);
  fChain->SetBranchAddress("p2RDcaZtoSV", &p2RDcaZtoSV, &b_p2RDcaZtoSV);
  //	fChain->SetBranchAddress("p2EtaIdx", &p2EtaIdx, &b_p2EtaIdx);
  //	fChain->SetBranchAddress("p2PtIdx", &p2PtIdx, &b_p2PtIdx);
  fChain->SetBranchAddress("p2Tpc", &p2Tpc, &b_p2Tpc);
  fChain->SetBranchAddress("kHft", &kHft, &b_kHft);
  fChain->SetBranchAddress("p1Hft", &p1Hft, &b_p1Hft);
  fChain->SetBranchAddress("p2Hft", &p2Hft, &b_p2Hft);
  //	fChain->SetBranchAddress("kTof", &kTof, &b_kTof);
  //	fChain->SetBranchAddress("p1Tof", &p1Tof, &b_p1Tof);
  //	fChain->SetBranchAddress("p2Tof", &p2Tof, &b_p2Tof);
  //	fChain->SetBranchAddress("sA", &sA, &b_sA);
  //	fChain->SetBranchAddress("sB", &sB, &b_sB);
  //	fChain->SetBranchAddress("dcakp1", &dcakp1, &b_dcakp1);
  //	fChain->SetBranchAddress("dcap1p2", &dcap1p2, &b_dcap1p2);
  //	fChain->SetBranchAddress("dcap2k", &dcap2k, &b_dcap2k);
  fChain->SetBranchAddress("mdV0Max", &mdV0Max, &b_mdV0Max);
  fChain->SetBranchAddress("mdV0Max_helix", &mdV0Max_helix, &b_mdV0Max_helix);
  Notify();
}

Bool_t nt_TMVA_Cuts::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void nt_TMVA_Cuts::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t nt_TMVA_Cuts::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef nt_cxx
