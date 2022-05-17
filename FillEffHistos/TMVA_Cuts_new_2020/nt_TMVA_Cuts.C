#define nt_TMVA_Cuts_cxx
#include "nt_TMVA_Cuts.h"
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>

#include "./TMVA_cuts/tmvaCuts.h" //contains analysis pre-cuts, set ReadMode = 1 to overwrite them with TMVA cuts

using namespace tmvaCuts;
using namespace std;

void nt_TMVA_Cuts::Loop()
{
  //   In a ROOT session, you can do:
  //      Root > .L nt.C
  //      Root > nt t
  //      Root > t.GetEntry(12); // Fill t data members with entry number 12
  //      Root > t.Show();       // Show values of entry 12
  //      Root > t.Show(16);     // Read and show values of entry 16
  //      Root > t.Loop();       // Loop on all entries
  //

  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch

  if (fChain == 0) return;

  // generated Dpm histograms
  TH1F *h_dpm_pt[9];
  for(unsigned int i = 0; i<9; i++)
  {
    h_dpm_pt[i] = new TH1F(Form("h_dpm_pt%i", i),Form("Generated_Dpm_cent_%i", i),nPtBins, pT_bins);
    h_dpm_pt[i]->Sumw2();

  }
  TH1F *h_dpm_phi = new TH1F("h_dpm_phi","Generated_Dpm_phi",64, -3.2, 3.2);
  TH2F *h_dpm_eta_phi = new TH2F("h_dpm_eta_phi","Generated_Dpm_eta_phi", 140, -1.2, 1.2, 64, -3.2, 3.2);

  h_dpm_phi->Sumw2();
  h_dpm_eta_phi->Sumw2();
  //---------------------------------------------------------------------------------------------
  // reconstructed Dpm histograms without TOF matching
  TH1F *h_r_dpm_pt[9];
  for(unsigned int i = 0; i<9; i++)
  {
    h_r_dpm_pt[i] = new TH1F(Form("h_r_dpm_pt%i", i),Form("Reconstructed_Dpm_no_TOF_cent_%i", i),nPtBins, pT_bins);
    h_r_dpm_pt[i]->Sumw2();

  }
  TH1F *h_r_dpm_phi = new TH1F("h_r_dpm_phi","Reconstructed_Dpm_no_TOF_phi",64, -3.2, 3.2);
  TH2F *h_r_dpm_eta_phi = new TH2F("h_r_dpm_eta_phi","Reconstructed_Dpm_no_TOF_eta_phi", 140, -1.2, 1.2, 64, -3.2, 3.2);

  h_r_dpm_phi->Sumw2();
  h_r_dpm_eta_phi->Sumw2();
  //----------------------------------------------------------------------------------------------------

  // reconstructed Dpm only with cuts, but without TPC or HFT matching
  TH1F *h_r_dpm_pt_no_match[9];
  for(unsigned int i = 0; i<9; i++)
  {
    h_r_dpm_pt_no_match[i] = new TH1F(Form("h_r_dpm_pt_no_match%i", i),Form("Reconstructed_Dpm_no_match_no_TOF_cent_%i", i),nPtBins, pT_bins);
    h_r_dpm_pt_no_match[i]->Sumw2();

  }

  //----------------------------------------------------------------------------------------------------
  
  TH1F *h_r_dpm_pt_geometric[9];
  for(unsigned int i = 0; i<9; i++)
  {
    h_r_dpm_pt_geometric[i] = new TH1F(Form("h_r_dpm_pt_geometric%i", i),Form("Reconstructed_Dpm_geometric_cuts_%i", i),nPtBins, pT_bins);
    h_r_dpm_pt_geometric[i]->Sumw2();
  }
       
    
  // reconstructed Dpm with TPC matching without HFT matching
  TH1F *h_r_dpm_pt_TPC_match[9];
  for(unsigned int i = 0; i<9; i++)
  {
    h_r_dpm_pt_TPC_match[i] = new TH1F(Form("h_r_dpm_pt_TPC_match%i", i),Form("Reconstructed_Dpm_TPC_match_no_TOF_cent_%i", i),nPtBins, pT_bins);
    h_r_dpm_pt_TPC_match[i]->Sumw2();

  }

  //----------------------------------------------------------------------------------------------------

  // reconstructed Dpm with HFT matching without TPC matching
  TH1F *h_r_dpm_pt_HFT_match[9];
  for(unsigned int i = 0; i<9; i++)
  {
    h_r_dpm_pt_HFT_match[i] = new TH1F(Form("h_r_dpm_pt_HFT_match%i", i),Form("Reconstructed_Dpm_HFT_match_no_TOF_cent_%i", i),nPtBins, pT_bins);
    h_r_dpm_pt_HFT_match[i]->Sumw2();

  }

  //----------------------------------------------------------------------------------------------------

  // reconstructed Dpm with TPC matching without HFT matching, without topo cuts
  TH1F *h_r_dpm_pt_TPC_match_no_topo_cuts[9];
  for(unsigned int i = 0; i<9; i++)
  {
    h_r_dpm_pt_TPC_match_no_topo_cuts[i] = new TH1F(Form("h_r_dpm_pt_TPC_match_no_topo_cuts%i", i),Form("Reconstructed_Dpm_TPC_match_no_TOF_no_topo_cent_%i", i),nPtBins, pT_bins);
    h_r_dpm_pt_TPC_match_no_topo_cuts[i]->Sumw2();

  }

  //----------------------------------------------------------------------------------------------------

  // reconstructed Dpm with HFT matching without TPC matching, without topo cuts
  TH1F *h_r_dpm_pt_HFT_match_no_topo_cuts[9];
  for(unsigned int i = 0; i<9; i++)
  {
    h_r_dpm_pt_HFT_match_no_topo_cuts[i] = new TH1F(Form("h_r_dpm_pt_HFT_match_no_topo_cuts%i", i),Form("Reconstructed_Dpm_HFT_match_no_TOF_no_topo_cent_%i", i),nPtBins, pT_bins);
    h_r_dpm_pt_HFT_match_no_topo_cuts[i]->Sumw2();

  }

  //----------------------------------------------------------------------------------------------------

  //invariant mass hitograms
  TH1F *h_D_inv_mass_base = new TH1F("D_inv_mass", "D_inv_mass", 40, 1.7, 2.1); //base histogram for D+- inv. mass

  TH1F *h_D_inv_mass_array[nPtBins][nCentBins]; //inv. mass histograms - same as in analysis

  for(unsigned int i = 0; i < nCentBins; i++)
  {
    for(unsigned int j = 0; j < nPtBins; j++)
    {
      //cout<<i<<" "<<j<<endl;
      h_D_inv_mass_array[j][i] = (TH1F*)h_D_inv_mass_base->Clone();
      h_D_inv_mass_array[j][i]->SetNameTitle(Form("D_inv_mass_cent_%i_pT_%i", i, j), Form("D_inv_mass_cent_%i_pT_%i", i, j));
    }
  }  
  //-------------------------------------------------------------------------------------------------------
  //load input files with cuts-----------------------------------------------------------------------------
  if(ReadMode != 0) //ReadMode == 0 is for analysis pre-cuts stored in tmvaCuts.h
  {
    ifstream CutsFile[nCentBins-1];

    for(unsigned int j = 0; j < nCentBins-1; j++)
    {
      ostringstream filename;

      // load cuts (0 - pre-cuts, 1 - ana, 2 - loose, 3 - tight)
      if(ReadMode == 1)
      {
        filename << "./TMVA_cuts/005_new_training_2020/analysis/cuts_cent_" << j << ".txt";
      }
      if(ReadMode == 2)
      {
        filename << "./TMVA_cuts/005_new_training_2020/loose/50_percent/cuts_loose_cent_" << j << ".txt";
      }
      if(ReadMode == 3)
      {
        filename << "./TMVA_cuts/005_new_training_2020/tight/50_percent/cuts_tight_cent_" << j << ".txt";
      }


      CutsFile[j].open(filename.str().c_str());
      if(!CutsFile[j].is_open())
      {
        cout<<"Failed to open file with cuts for centrality "<<j<<"!"<<endl;
        return;
      }
    }

    string CutsLine;

    for(unsigned int cent = 0; cent < nCentBins-1; cent++)
    {
      int lineNo = 1;

      while(getline(CutsFile[cent], CutsLine))
      {

        stringstream CutsLineStream(CutsLine);
        for(unsigned int pTbin = 0; pTbin < nPtBins; pTbin++)
        {
          if(lineNo==1)
          {
            CutsLineStream>>k_dca_cut[cent][pTbin];
          }
          if(lineNo==2)
          {
            CutsLineStream>>pi1_dca_cut[cent][pTbin];
            pi2_dca_cut[cent][pTbin] = pi1_dca_cut[cent][pTbin];
          }
          if(lineNo==3)
          {
            CutsLineStream>>mdcaMax_cut[cent][pTbin];
          }
          if(lineNo==4)
          {
            CutsLineStream>>D_decayL_min_cut[cent][pTbin];
          }
          if(lineNo==5)
          {
            CutsLineStream>>D_cos_theta_cut[cent][pTbin];
          }

        }//end for (pT)

        lineNo++;

      }//end while
    }//end for (centrality)
  }//end if for ReadMode
  //----------------------------------------------------------------------------------------



  // now cuts
  TH1D *n_cuts = new TH1D("n_cuts","n_cuts",100,0,100);

  Long64_t nentries = fChain->GetEntries();
  cout << "#events: " << nentries << endl;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    if(jentry%(nentries/10)==0) std::cout<<((jentry+10)*100/nentries)<<" % done..."<<std::endl;		
    nb = fChain->GetEntry(jentry);
    nbytes += nb;

    if(pTspectrum == 1) w = 1; //pTspectrum == 1 -> Levy input pt spectrum - already realistic, no need to use weights, set in run_nt_TMVA_Cuts.C

    if(cent > -0.5 and cent < 0.5) h_dpm_pt[0]->Fill(rPt,w);
    if(cent > 0.5 and cent < 1.5) h_dpm_pt[1]->Fill(rPt,w);
    if(cent > 1.5 and cent < 2.5) h_dpm_pt[2]->Fill(rPt,w);
    if(cent > 2.5 and cent < 3.5) h_dpm_pt[3]->Fill(rPt,w);
    if(cent > 3.5 and cent < 4.5) h_dpm_pt[4]->Fill(rPt,w);
    if(cent > 4.5 and cent < 5.5) h_dpm_pt[5]->Fill(rPt,w);
    if(cent > 5.5 and cent < 6.5) h_dpm_pt[6]->Fill(rPt,w);
    if(cent > 6.5 and cent < 7.5) h_dpm_pt[7]->Fill(rPt,w);
    if(cent > 7.5 and cent < 8.5) h_dpm_pt[8]->Fill(rPt,w);


    int centBin_9 = (int)cent; //set centrality bin)



    h_dpm_phi->Fill(phi);
    h_dpm_eta_phi->Fill(eta, phi);

    n_cuts->Fill(1);

    if (cent > con_cent_up) continue;
    if (cent < con_cent_down) continue;
    n_cuts->Fill(2);

    if ( !(fabs(rY)<1) ) continue; //cut on reco rapidity of the D meson

    if (fabs(kREta) > 1 || fabs(p1REta) > 1 || fabs(p2REta) > 1) continue;


    // kRPt
    if (kRPt < con_kRPt) continue;
    n_cuts->Fill(11);
    // p1RPt
    if (p1RPt < con_pRPt) continue;
    n_cuts->Fill(12);
    // p2RPt
    if (p2RPt < con_pRPt) continue;
    n_cuts->Fill(13);


    h_r_dpm_pt_geometric[centBin_9]->Fill(rPt, w);

    // HFT
    int DpmHFT = 1;

    if (kHft != 1 || p1Hft != 1 || p2Hft != 1) DpmHFT = 0; //continue;
    n_cuts->Fill(3);
    // TPC
    int DpmTPC = 1;

    if (kTpc != 1 || p1Tpc != 1 || p2Tpc != 1) DpmTPC = 0;//continue;

    n_cuts->Fill(4);



    int ptBin = -1; //find pT bin (as in real data, not TMVA!!!)

    for(int j = 0; j < nPtBins; j++) //loop over pT bins
    {
      if(rPt > pT_bins[j] && rPt <= pT_bins[j+1])
      {
        ptBin = j;
      }

    }




    int centBin = -1;

    if( centBin_9 == 7 || centBin_9 == 8 ) centBin = 0; //0-10%
    if( centBin_9 == 4 || centBin_9 == 5 || centBin_9 == 6  ) centBin = 1; //10-40%
    if( centBin_9 == 0 || centBin_9 == 1 || centBin_9 ==  2 || centBin_9 ==  3 ) centBin = 2; //40-80%


    if(ptBin < 0 || centBin < 0 ) continue; //reject Dpm with pT and centrality out of range


    if(DpmHFT == 1)
    {
      h_r_dpm_pt_HFT_match_no_topo_cuts[centBin_9]->Fill(rPt, w);
    }

    if(DpmTPC == 1)
    {
      h_r_dpm_pt_TPC_match_no_topo_cuts[centBin_9]->Fill(rPt, w);
    }

    //for current version of results
    if(  p1RDca > 2000 || p2RDca > 2000 || kRDca > 2000  ) continue;

    if((dcaDaughters < mdcaMax_cut[centBin][ptBin]*10000) &&
        (decayLength > D_decayL_min_cut[centBin][ptBin]*10000 ) &&
        (cosTheta > D_cos_theta_cut[centBin][ptBin] ) &&
        (mdV0Max <  250) && //manually tuned cut on Delta_max
        ( p1RDca > pi1_dca_cut[centBin][ptBin]*10000 && p2RDca > pi2_dca_cut[centBin][ptBin]*10000 && kRDca > k_dca_cut[centBin][ptBin]*10000 ) )
    {

      h_r_dpm_pt_no_match[centBin_9]->Fill(rPt, w);

      if(DpmHFT == 1)
      {
        h_r_dpm_pt_HFT_match[centBin_9]->Fill(rPt, w);
      }

      if(DpmTPC == 1)
      {
        h_r_dpm_pt_TPC_match[centBin_9]->Fill(rPt, w);
      }

      if(DpmHFT == 1 && DpmTPC == 1)
      {
        h_r_dpm_pt[centBin_9]->Fill(rPt,w);

        h_r_dpm_phi->Fill(phi);
        h_r_dpm_eta_phi->Fill(eta, phi);

        h_D_inv_mass_array[ptBin][centBin]->Fill(rM);
      }


    }


  } //end entries loop
  //_____________________________________________________________________________

  //write histograms to a new file
  TFile *f = new TFile(out_file_name, "recreate");

  //generated Dpm
  for(unsigned int j = 0; j<9; j++ )
  {
    h_dpm_pt[j]->Write();
  }
  h_dpm_phi->Write();
  h_dpm_eta_phi->Write();

  //reconstructed Dpm, no TOF matching
  for(unsigned int j = 0; j<9; j++ )
  {
    h_r_dpm_pt[j]->Write();
  }
  h_r_dpm_phi->Write();
  h_r_dpm_eta_phi->Write();

  //reconstructed Dpm with only geometrical cuts
  for(unsigned int j = 0; j<9; j++ )
  {
    h_r_dpm_pt_geometric[j]->Write();
  }

  
  //reconstructed Dpm, only cuts no matching  
  for(unsigned int j = 0; j<9; j++ )
  {
    h_r_dpm_pt_no_match[j]->Write();
  }

  //reconstructed Dpm, cuts and TPC match  
  //no totpo cuts, with eta and pT
  for(unsigned int j = 0; j<9; j++ )
  {
    h_r_dpm_pt_TPC_match_no_topo_cuts[j]->Write();
  }
  //pure TPC matching
  for(unsigned int j = 0; j<9; j++ )
  {
    h_r_dpm_pt_TPC_match[j]->Write();
  }

  
  //reconstructed Dpm, cuts and HFT match
  //no topo cuts,  with eta and pT
  for(unsigned int j = 0; j<9; j++ )
  {
    h_r_dpm_pt_HFT_match_no_topo_cuts[j]->Write();
  }
  //pure HFT matching
  for(unsigned int j = 0; j<9; j++ )
  {
    h_r_dpm_pt_HFT_match[j]->Write();
  }

  

  for(unsigned int j = 0; j < nCentBins; j++)
  {
    for(unsigned int k = 0; k < nPtBins; k++)
    {
      h_D_inv_mass_array[k][j]->Write();
    }
  }


  n_cuts->Write();


  f->Close();
}
