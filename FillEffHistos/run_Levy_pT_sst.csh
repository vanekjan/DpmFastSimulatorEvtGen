#!/bin/csh

starver SL19e

#use with sst stream only

#Levy pT, analysis TMVA cuts, analysis daugter pT cut
root4star -b -q run_nt_TMVA_Cuts.C\(\"$1\",\"./output/sst/Dpm.out_eff_TMVA_ana_cuts_sst_Levy_pT.toyMc.root\",1,1,0\)

#Levy pT, loose TMVA cuts, analysis daugter pT cut
root4star -b -q run_nt_TMVA_Cuts.C\(\"$1\",\"./output/sst/Dpm.out_eff_TMVA_loose_cuts_sst_Levy_pT.toyMc.root\",1,2,0\)

#Levy pT, tight TMVA cuts, analysis daugter pT cut
root4star -b -q run_nt_TMVA_Cuts.C\(\"$1\",\"./output/sst/Dpm.out_eff_TMVA_tight_cuts_sst_Levy_pT.toyMc.root\",1,3,0\)

#Levy pT,  analysis TMVA cuts, varied daugter pT cut
root4star -b -q run_nt_TMVA_Cuts.C\(\"$1\",\"./output/sst/Dpm.out_eff_TMVA_ana_cuts_sst_pT_cut_var_500_Levy_pT.toyMc.root\",1,0,1\)

