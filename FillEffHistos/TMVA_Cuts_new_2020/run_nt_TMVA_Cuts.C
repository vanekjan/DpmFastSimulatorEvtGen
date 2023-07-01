#include <TString.h>

// folders with output from fast-sim

// 2021-11-09_11-01_physics_new_inputs_Vz_weight_new_embedd/
// 2021-11-15_03-02_physics_new_inputs_Vz_weight_new_embedd_Levy_pT/
//
// 2021-09-15_05-36_sst_new_inputs_Vz_weight/
// 2021-10-22_09-44_sst_new_inputs_Vz_weight_Levy_pT/

void run_nt_TMVA_Cuts(TString infile = "/star/u/vanekjan/pwg/vanekjan/myDpmEvtGenFastSim/myOutput/2023-07-01_08-24/merge/output.root", 
                      TString outfile = "/star/u/vanekjan/pwg/vanekjan/myDpmEvtGenFastSim/myOutput/Histo_output/Code_QA/Dpm.out_eff.toyMc.root") 
  {
  std::cout << "start " << std::endl;
	gROOT->ProcessLine(".L nt_TMVA_Cuts.C++");

	TFile *f_input = new TFile(infile, "open");
	TTree *tree_nt = (TTree*)f_input->Get("nt");
	nt_TMVA_Cuts n(tree_nt);
	n.Set_out_file_name(outfile);

//centrality range in bins defined same as in StRefMultCorr
	n.Set_con_cent_up(9);
	n.Set_con_cent_down(-1);
	
//analysis pT cut
  n.Set_con_kRPt(0.3);
  n.Set_con_pRPt(0.3);

//for pT cut systematic error
//  n.Set_con_kRPt(0.5);
//  n.Set_con_pRPt(0.5);

  n.Set_ReadMode(2);//ReadMode = 0 - pre-cuts, ReadMode = 1 - TMVA analysis cuts, ReadMode = 2 - TMVA loose cuts, ReadMode = 3 - TMVA tight cuts, ReadMode = 4 - old topo. cuts
  n.Set_pTspectrum(0);// pTspectrum = 0 - flat input pT spectrum (use weights), 1 - Levy input spectrum (do not use weights - already realistic pT shape)


//topological cuts set in TMVA_cuts
//_____________________________________________________________________________
	
	n.Loop();
	std::cout << "end " << std::endl;
	f_input->Close();
}
