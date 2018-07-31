#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>

using namespace std;

void SaveFastSimHistos()
{
  TFile *inFile = new TFile("/star/u/vanekjan/500GBStorage/vanekjan/myDpmEvtGenFastSim/myOutput/2018-07-31_13-54/merge/output.root", "read");    

    TTree          *fChain = (TTree*)inFile->Get("nt");   //!pointer to the analyzed TTree or TChain
		
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
		Float_t         decayLength;
		Float_t         dcaDpmToPv;
		Float_t         cosTheta;
		Float_t         angle12;
		Float_t         cosThetaStar;
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
		Float_t         kRSDca;
		Float_t         kRDcaXY;
		Float_t         kRDcaZ;
		Float_t         kEtaIdx;
		Float_t         kPtIdx;
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
		Float_t         p1RSDca;
		Float_t         p1RDcaXY;
		Float_t         p1RDcaZ;
		Float_t         p1EtaIdx;
		Float_t         p1PtIdx;
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
		Float_t         p2RSDca;
		Float_t         p2RDcaXY;
		Float_t         p2RDcaZ;
		Float_t         p2EtaIdx;
		Float_t         p2PtIdx;
		Float_t         p2Tpc;
		Float_t         kHft;
		Float_t         p1Hft;
		Float_t         p2Hft;
		Float_t         kTof;
		Float_t         p1Tof;
		Float_t         p2Tof;
		Float_t         sA;
		Float_t         sB;
		Float_t         dcakp1;
		Float_t         dcap1p2;
		Float_t         dcap2k;
		Float_t         mdV0Max;

  fChain->SetBranchAddress("cent", &cent);
	fChain->SetBranchAddress("vx", &vx);
	fChain->SetBranchAddress("vy", &vy);
	fChain->SetBranchAddress("vz", &vz);
	fChain->SetBranchAddress("vzIdx", &vzIdx);
	fChain->SetBranchAddress("pid", &pid);
	fChain->SetBranchAddress("w", &w);
	fChain->SetBranchAddress("m", &m);
	fChain->SetBranchAddress("pt", &pt);
	fChain->SetBranchAddress("eta", &eta);
	fChain->SetBranchAddress("y", &y);
	fChain->SetBranchAddress("phi", &phi);
	fChain->SetBranchAddress("v0x", &v0x);
	fChain->SetBranchAddress("v0y", &v0y);
	fChain->SetBranchAddress("v0z", &v0z);
	fChain->SetBranchAddress("rM", &rM);
	fChain->SetBranchAddress("rPt", &rPt);
	fChain->SetBranchAddress("rEta", &rEta);
	fChain->SetBranchAddress("rY", &rY);
	fChain->SetBranchAddress("rPhi", &rPhi);
	fChain->SetBranchAddress("rV0x", &rV0x);
	fChain->SetBranchAddress("rV0y", &rV0y);
	fChain->SetBranchAddress("rV0z", &rV0z);
	fChain->SetBranchAddress("dcaDaughters", &dcaDaughters);
	fChain->SetBranchAddress("decayLength", &decayLength);
	fChain->SetBranchAddress("dcaDpmToPv", &dcaDpmToPv);
	fChain->SetBranchAddress("cosTheta", &cosTheta);
	fChain->SetBranchAddress("angle12", &angle12);
	fChain->SetBranchAddress("cosThetaStar", &cosThetaStar);
	fChain->SetBranchAddress("kM", &kM);
	fChain->SetBranchAddress("kPt", &kPt);
	fChain->SetBranchAddress("kEta", &kEta);
	fChain->SetBranchAddress("kY", &kY);
	fChain->SetBranchAddress("kPhi", &kPhi);
	fChain->SetBranchAddress("kDca", &kDca);
	fChain->SetBranchAddress("kRM", &kRM);
	fChain->SetBranchAddress("kRPt", &kRPt);
	fChain->SetBranchAddress("kREta", &kREta);
	fChain->SetBranchAddress("kRY", &kRY);
	fChain->SetBranchAddress("kRPhi", &kRPhi);
	fChain->SetBranchAddress("kRVx", &kRVx);
	fChain->SetBranchAddress("kRVy", &kRVy);
	fChain->SetBranchAddress("kRVz", &kRVz);
	fChain->SetBranchAddress("kRDca", &kRDca);
	fChain->SetBranchAddress("kRSDca", &kRSDca);
	fChain->SetBranchAddress("kRDcaXY", &kRDcaXY);
	fChain->SetBranchAddress("kRDcaZ", &kRDcaZ);
	fChain->SetBranchAddress("kEtaIdx", &kEtaIdx);
	fChain->SetBranchAddress("kPtIdx", &kPtIdx);
	fChain->SetBranchAddress("kTpc", &kTpc);
	fChain->SetBranchAddress("p1M", &p1M);
	fChain->SetBranchAddress("p1Pt", &p1Pt);
	fChain->SetBranchAddress("p1Eta", &p1Eta);
	fChain->SetBranchAddress("p1Y", &p1Y);
	fChain->SetBranchAddress("p1Phi", &p1Phi);
	fChain->SetBranchAddress("p1Dca", &p1Dca);
	fChain->SetBranchAddress("p1RM", &p1RM);
	fChain->SetBranchAddress("p1RPt", &p1RPt);
	fChain->SetBranchAddress("p1REta", &p1REta);
	fChain->SetBranchAddress("p1RY", &p1RY);
	fChain->SetBranchAddress("p1RPhi", &p1RPhi);
	fChain->SetBranchAddress("p1RVx", &p1RVx);
	fChain->SetBranchAddress("p1RVy", &p1RVy);
	fChain->SetBranchAddress("p1RVz", &p1RVz);
	fChain->SetBranchAddress("p1RDca", &p1RDca);
	fChain->SetBranchAddress("p1RSDca", &p1RSDca);
	fChain->SetBranchAddress("p1RDcaXY", &p1RDcaXY);
	fChain->SetBranchAddress("p1RDcaZ", &p1RDcaZ);
	fChain->SetBranchAddress("p1EtaIdx", &p1EtaIdx);
	fChain->SetBranchAddress("p1PtIdx", &p1PtIdx);
	fChain->SetBranchAddress("p1Tpc", &p1Tpc);
	fChain->SetBranchAddress("p2M", &p2M);
	fChain->SetBranchAddress("p2Pt", &p2Pt);
	fChain->SetBranchAddress("p2Eta", &p2Eta);
	fChain->SetBranchAddress("p2Y", &p2Y);
	fChain->SetBranchAddress("p2Phi", &p2Phi);
	fChain->SetBranchAddress("p2Dca", &p2Dca);
	fChain->SetBranchAddress("p2RM", &p2RM);
	fChain->SetBranchAddress("p2RPt", &p2RPt);
	fChain->SetBranchAddress("p2REta", &p2REta);
	fChain->SetBranchAddress("p2RY", &p2RY);
	fChain->SetBranchAddress("p2RPhi", &p2RPhi);
	fChain->SetBranchAddress("p2RVx", &p2RVx);
	fChain->SetBranchAddress("p2RVy", &p2RVy);
	fChain->SetBranchAddress("p2RVz", &p2RVz);
	fChain->SetBranchAddress("p2RDca", &p2RDca);
	fChain->SetBranchAddress("p2RSDca", &p2RSDca);
	fChain->SetBranchAddress("p2RDcaXY", &p2RDcaXY);
	fChain->SetBranchAddress("p2RDcaZ", &p2RDcaZ);
	fChain->SetBranchAddress("p2EtaIdx", &p2EtaIdx);
	fChain->SetBranchAddress("p2PtIdx", &p2PtIdx);
	fChain->SetBranchAddress("p2Tpc", &p2Tpc);
	fChain->SetBranchAddress("kHft", &kHft);
	fChain->SetBranchAddress("p1Hft", &p1Hft);
	fChain->SetBranchAddress("p2Hft", &p2Hft);
	fChain->SetBranchAddress("kTof", &kTof);
	fChain->SetBranchAddress("p1Tof", &p1Tof);
	fChain->SetBranchAddress("p2Tof", &p2Tof);
	fChain->SetBranchAddress("sA", &sA);
	fChain->SetBranchAddress("sB", &sB);
	fChain->SetBranchAddress("dcakp1", &dcakp1);
	fChain->SetBranchAddress("dcap1p2", &dcap1p2);
	fChain->SetBranchAddress("dcap2k", &dcap2k);
	fChain->SetBranchAddress("mdV0Max", &mdV0Max);
  //______________________________________________________________________________________________________________________

	if (fChain == 0) return;

  const int nPtBins = 12;  
  double pT_bins[nPtBins+1] = { 1., 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10. }; //pT binning - for analysis

  TFile *outFile = new TFile("InputHistosEvtGen_new_02.root", "recreate");

  //const int nPtBins = 6;
  //float const pT_bins[nPtBins+1] = {0., 1., 2., 3., 5., 7., 10.}; //my TMVA pT bins, for D0 nSignal estimation

	
	TH1F *h_dpm_phi = new TH1F("h_dpm_phi","Generated_Dpm_phi",64, -3.2, 3.2);
  TH1F *h_dpm_eta = new TH1F("h_dpm_eta","Generated_Dpm_eta", 140, -1.2, 1.2);
	TH2F *h_dpm_eta_phi = new TH2F("h_dpm_eta_phi","Generated_Dpm_eta_phi", 140, -1.2, 1.2, 64, -3.2, 3.2);

	h_dpm_phi->Sumw2();
	h_dpm_eta_phi->Sumw2();


  TH1F *h_k_DCA = new TH1F("h_k_DCA", "h_k_DCA", 100, 0, 2000);
  TH1F *h_pi1_DCA = new TH1F("h_pi1_DCA", "h_pi1_DCA", 100, 0, 2000);
  TH1F *h_pi2_DCA = new TH1F("h_pi2_DCA", "h_pi2_DCA", 100, 0, 2000);

  TH1F *h_k_rDCA = new TH1F("h_k_rDCA", "h_k_rDCA", 100, 0, 2000);
  TH1F *h_pi1_rDCA = new TH1F("h_pi1_rDCA", "h_pi1_rDCA", 100, 0, 2000);
  TH1F *h_pi2_rDCA = new TH1F("h_pi2_rDCA", "h_pi2_rDCA", 100, 0, 2000);

  TH1F *h_DecayLength = new TH1F("h_DecayLength", "h_DecayLength", 100, 0, 2000);

  TH1F *h_PairDCAmax = new TH1F("h_PairDCAmax", "h_PairDCAmax", 200, 0, 200);

  TH1F *h_cosTheta = new TH1F("h_cosTheta", "h_cosTheta", 200, -1, 1);

  TH1F *h_dV0Max = new TH1F("h_dV0Max", "h_dV0Max", 100, 0, 250);

  TH1F *h_k_rPt = new TH1F("h_k_rPt", "h_k_rPt", 100, 0, 10);
  TH1F *h_pi1_rPt = new TH1F("h_pi1_rPt", "h_pi1_rPt", 100, 0, 10);
  TH1F *h_pi2_rPt = new TH1F("h_pi2_rPt", "h_pi2_rPt", 100, 0, 10);

  TH1F *h_k_rEta = new TH1F("h_k_rEta", "h_k_rEta", 100, -1.1, 1.1);
  TH1F *h_pi1_rEta = new TH1F("h_pi1_rEta", "h_pi1_rEta", 100, -1.1, 1.1);
  TH1F *h_pi2_rEta = new TH1F("h_pi2_rEta", "h_pi2_rEta", 100, -1.1, 1.1);

  TH1F *h_v0z = new TH1F("h_v0z", "h_v0z", 100, -60000, 60000);

  TH1D *h_TPC_match = new TH1D("h_TPC_match", "h_TPC_match", 2, 0, 2);
  TH1D *h_HFT_match = new TH1D("h_HFT_match", "h_HFT_match", 2, 0, 2);

  TH1D *h_HFT_match_kaon = new TH1D("h_HFT_match_kaon", "h_HFT_match_kaon", 2, 0, 2);



	Long64_t nentries = fChain->GetEntries();
	cout << "#events: " << nentries << endl;
  
  outFile->cd();

  Long64_t TPC[2] = {0,0}; 

	for (Long64_t jentry=0; jentry<nentries;jentry++) 
  {
		Long64_t ientry = fChain->GetEntry(jentry);
		if (ientry < 0) break;
    if(jentry%(nentries/10)==0) std::cout<<((jentry+10)*100/nentries)<<" % done..."<<std::endl;

    fChain->GetEntry(jentry);

      

		h_dpm_phi->Fill(phi);
    h_dpm_eta->Fill(eta);
		h_dpm_eta_phi->Fill(eta, phi);

		h_v0z->Fill(v0z);

    h_k_rEta->Fill(kREta);
    h_pi1_rEta->Fill(p1REta);
    h_pi2_rEta->Fill(p2REta);

    
		h_k_rPt->Fill(kRPt);
    h_pi1_rPt->Fill(p1RPt);
    h_pi2_rPt->Fill(p2RPt);

    if(kHft != 1)
    {
      h_HFT_match_kaon->Fill(0.5);
    }
    else
    {
      h_HFT_match_kaon->Fill(1.5);
    }
    
    if(kHft != 1 || p1Hft != 1 || p2Hft != 1)
    {
      h_HFT_match->Fill(0.5);
    }
    else
    {
      h_HFT_match->Fill(1.5);
    }
    
		if(kTpc != 1 || p1Tpc != 1 || p2Tpc != 1)
    {
      h_TPC_match->Fill(0.5);
      TPC[0] += 1;
    }
    else
    {
      h_TPC_match->Fill(1.5);
      TPC[1] += 1;
    }

		
		//if (kTof != 1 || p1Tof != 1 || p2Tof != 1) continue;
		
		h_cosTheta->Fill(cosTheta);
		
		h_PairDCAmax->Fill(dcaDaughters);
		
		h_DecayLength->Fill(decayLength);
		
		h_k_DCA->Fill(kDca);
    h_pi1_DCA->Fill(p1Dca);
    h_pi2_DCA->Fill(p2Dca);

    h_k_rDCA->Fill(kRDca);
    h_pi1_rDCA->Fill(p1RDca);
    h_pi2_rDCA->Fill(p2RDca);
		

		
		h_dV0Max->Fill(mdV0Max);
		    
		
	} //end entries loop

  cout<<"No TPC match: "<<TPC[0]<<endl;
  cout<<"Good TPC match: "<<TPC[1]<<endl;

  h_dpm_phi->Write();
    h_dpm_eta->Write();
		h_dpm_eta_phi->Write();

		h_v0z->Write();

    h_k_rEta->Write();
    h_pi1_rEta->Write();
    h_pi2_rEta->Write();

		h_k_rPt->Write();
    h_pi1_rPt->Write();
    h_pi2_rPt->Write();

    h_HFT_match->Write();
    h_HFT_match_kaon->Write();		
    h_TPC_match->Write();
		
		h_cosTheta->Write();		
		h_PairDCAmax->Write();		
		h_DecayLength->Write();
    h_dV0Max->Write();
		
		h_k_DCA->Write();
    h_pi1_DCA->Write();
    h_pi2_DCA->Write();

    h_k_rDCA->Write();
    h_pi1_rDCA->Write();
    h_pi2_rDCA->Write();


	outFile->Close();
  inFile->Close();
}
