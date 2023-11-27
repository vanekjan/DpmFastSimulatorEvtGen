#!/bin/csh

starver SL19e

#use with physics stream only

#flat pT, analysis TMVA cuts, analysis daugter pT cut
root4star -b -q run_nt_TMVA_Cuts.C\(\"$1\",\"./output/physics/Dpm.out_eff_TMVA_ana_cuts_physics.toyMc.root\",0,1,0\)

#flat pT, loose TMVA cuts, analysis daugter pT cut
root4star -b -q run_nt_TMVA_Cuts.C\(\"$1\",\"./output/physics/Dpm.out_eff_TMVA_loose_cuts_physics.toyMc.root\",0,2,0\)

#flat pT, tight TMVA cuts, analysis daugter pT cut
root4star -b -q run_nt_TMVA_Cuts.C\(\"$1\",\"./output/physics/Dpm.out_eff_TMVA_tight_cuts_physics.toyMc.root\",0,3,0\)

#flat pT,  analysis TMVA cuts, varied daugter pT cut
root4star -b -q run_nt_TMVA_Cuts.C\(\"$1\",\"./output/physics/Dpm.out_eff_TMVA_ana_cuts_physics_pT_cut_var_500.toyMc.root\",0,0,1\)

