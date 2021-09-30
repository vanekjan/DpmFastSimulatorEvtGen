#include <TString.h>
//#include <TFile.h>
//#include <TTree.h>
//#include "nt.h"
//#include <iostream>
//#include"TROOT.h"

// /star/u/vanekjan/500GBStorage/vanekjan/myDpmEvtGenFastSim/myOutput/2018-07-31_13-54_EvtGen_final/merge/output.root //physics stream
// /star/u/vanekjan/500GBStorage/vanekjan/myDpmEvtGenFastSim/myOutput/2019-03-21_04-15_EvtGen_new_pT_resolution_test/merge/output.root

// /star/u/vanekjan/500GBStorage/vanekjan/myDpmEvtGenFastSim/myOutput/2019-09-10_12-37/merge/output.root //physics stream - new (corrected Delta_max)

// /star/u/vanekjan/500GBStorage/vanekjan/myDpmEvtGenFastSim/myOutput/2019-08-10_12-26/merge //sst+nosst streams
//
//
// /star/u/vanekjan/500GBStorage/vanekjan/myDpmEvtGenFastSim/myOutput/2019-08-22_05-58/merge/output.root //real p_T 02 (high p_T)
//
// 2019-10-17_10-15_new_DCAxy //new Run16 DCAxy distributions
//
//
// 2020-05-12_04-44_new_inputs_from_data_full_stat_sst_nosst/
//
// 2020-05-12_03-58_new_inputs_from_data_full_stat_physics/
//
// 2021-09-15_04-29_physics_new_inputs_Vz_weight/
//
// 2021-09-15_05-36_sst_new_inputs_Vz_weight/

void run_nt_TMVA_Cuts(TString infile = "/star/u/vanekjan/pwg/vanekjan/myDpmEvtGenFastSim/myOutput/2021-09-15_04-29_physics_new_inputs_Vz_weight/merge/output.root", 
                      TString outfile = "/star/u/vanekjan/pwg/vanekjan/myDpmEvtGenFastSim/myOutput/Histo_output/TMVA/physics/TMVA_cuts_2020/Delta_max_cut_250/Dpm.out_eff_new_11_pT_bins_TMVA_ana_cuts_physics_new_prod_no_decL_up_cut_DCA_up_cut_loose_Delta_max_Vz_refmult_weight_new_hist.toyMc.root") 
  {
  std::cout << "start " << std::endl;
	gROOT->ProcessLine(".L nt_TMVA_Cuts.C++");
	//gROOT->ProcessLine(".L nt.C");
	TFile *f_input = new TFile(infile, "open");
	TTree *tree_nt = (TTree*)f_input->Get("nt");
	nt_TMVA_Cuts n(tree_nt);
	n.Set_out_file_name(outfile);
	//n.Set_con_cent(8);
	n.Set_con_cent_up(9);
	n.Set_con_cent_down(-1);
	n.Set_con_v0z(60000.);
	
//analysis pT cut
  n.Set_con_kRPt(0.3);
  n.Set_con_pRPt(0.3);

//for pT cut systematic error
//  n.Set_con_kRPt(0.5);
//  n.Set_con_pRPt(0.5);

  n.Set_ReadMode(1);//ReadMode = 0 - pre-cuts, ReadMode = 1 - TMVA analysis cuts, ReadMode = 2 - TMVA loose cuts, ReadMode = 3 - TMVA tight cuts, ReadMode = 4 - old topo. cuts

  //set topological cuts in tmvaCuts.h
	//_____________________________________________________________________________
	
	n.Loop();
	std::cout << "end " << std::endl;
	f_input->Close();
}
