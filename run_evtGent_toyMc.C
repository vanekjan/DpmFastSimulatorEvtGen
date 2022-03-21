void run_evtGent_toyMc(int nEvents = 1e3)
{

  //gROOT->ProcessLine(".L bfc.C");
  /*  {
      TString simple = "y2014 geant gstar usexgeom agml ";
      bfc(0, simple );
      }
      */
  //delete chain;

  // gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  // loadSharedLibraries();

  gSystem->Load("libStTableUtilities.so");
  gSystem->Load("libTable");
  gSystem->Load("libPhysics");
  gSystem->Load("St_base");
  gSystem->Load("StChain");
  gSystem->Load("St_Tables");
  gSystem->Load("StUtilities");        // new addition 22jul99
  //gSystem->Load("StTreeMaker");
  //gSystem->Load("StIOMaker");
  // gSystem->Load("StarClassLibrary");
  //gSystem->Load("StTriggerDataMaker"); // new starting from April 2003
  //gSystem->Load("StBichsel");
  //gSystem->Load("StEvent");
  //gSystem->Load("StEventUtilities");
  //gSystem->Load("StDbLib");
  //gSystem->Load("StEmcUtil");
  //gSystem->Load("StTofUtil");
  //gSystem->Load("StPmdUtil");
  //gSystem->Load("StPreEclMaker");
  //gSystem->Load("StStrangeMuDstMaker");
  //gSystem->Load("StMuDSTMaker");
  //


  gSystem->Load( "libVMC.so" );
  gSystem->Load( "libSt_g2t.so" );
  gSystem->Load( "libSt_geant_Maker.so" );
  gSystem->Load( "StarGeneratorUtil.so" );
  gSystem->Load( "StarGeneratorEvent.so" );
  gSystem->Load( "StarGeneratorBase.so" );
  // gSystem->Load( "libMathMore.so"   );  
  //gSystem->Load( "libHijing1_383.so");
  //gSystem->Load( "libKinematics.so");
  //gSystem->Load( "xgeometry.so"     );

  gSystem->Load("libHepMC2_06_09.so");
  gSystem->Load("libPythia8_1_86.so");
  gSystem->Load("libPhotos3_61.so");
  gSystem->Load("libTauola1_1_5.so");
  gSystem->Load(".sl73_gcc485/lib/libEvtGen1_06_00.so");

  gInterpreter->AddIncludePath("StRoot/StarGenerator/EvtGen1_06_00");
  gInterpreter->AddIncludePath("$STAR/StRoot/StarGenerator/HepMC2_06_09");
  //gInterpreter->AddIncludePath("/common/star/star64/packages/SL16k/StRoot/StarGenerator/HepMC2_06_09"); //probably just for PDSF

  //cout<<"Processing evtGen_toyMc with..."<<endl;
  //
  //StChain *chain;
  //chain = new StChain();
  //
  StChain chain("myChain");

  gROOT->ProcessLine(Form(".x evtGen_toyMc.C+(%i , %i, %i, %i)",nEvents, 40, 80, 1)); // added possibility to set a custom centrality range and shape of input pT spectrum (0 - flat, 1 - Levy)

  //chain->Finish();
  chain.Finish();
  //delete chain;
  //


  return;
}
