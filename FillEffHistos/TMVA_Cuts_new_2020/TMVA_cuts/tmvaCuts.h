#ifndef TMVACUTS_H
#define TMVACUTS_H

//analysis pre-cuts

namespace tmvaCuts
{
    //int const totalNumberOfEvents = 90.e6; // Run14 dataset - probably not used
    int const nCentBins = 4;
    int   const nPtBins = 11;
    //float const pT_bins[nPtBins+1] = { 1., 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10. }; //my old pT bins, 12 bins
    //float const pT_bins[nPtBins+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5., 6., 7., 10. }; //my pT bins, new (wrong)
    float const pT_bins[nPtBins+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5., 6., 8., 10. }; //my pT bins, new (same as Guannan, 11 bins)

    //int   const nPtBins_TMVA = 12;
    //float const pT_bins_TMVA[nPtBins_TMVA+1] = { 1., 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10. }; //my pT bins, for analysis and some versions of TMVA

    const int nPtBins_TMVA = 11;
    //float const pT_bins_TMVA[nPtBins_TMVA+1] = {0., 1., 2., 3., 5., 7., 10.}; //my TMVA pT bins
    //float const pT_bins_TMVA[nPtBins_TMVA+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5., 6., 7., 10. }; //my TMVA pT bins (old, wrong)
    float const pT_bins_TMVA[nPtBins_TMVA+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5., 6., 8., 10. }; //my TMVA pT bins (new, same as Guannan in analysis)

    //float const PtBins[nPtBins+1] = {0., 1., 2., 3., 5., 7., 10.}; //original Xinyue, 6 pT bins
    //float const BDTmin[nPtBins] = {0.2123, 0.1872, 0.1752, 0.1666, 0.2061, 0.1709, 0.2381};

    //-------------deifne cuts--------------------------------------------------

     //Dpm
    Float_t D_y_cut = 1.;

    //TPC track quality cuts
    //Int_t pi1_nHitFit_cut, pi2_nHitFit_cut, k_nHitFit_cut;
    Int_t pi1_nHitFit_cut = 20;
    Int_t pi2_nHitFit_cut = 20;
    Int_t k_nHitFit_cut = 20;

    //Float_t pi1_nHitsMax_cut, pi2_nHitsMax_cut, k_nHitsMax_cut;
    Float_t pi1_nHitsMax_cut = 0.52;
    Float_t pi2_nHitsMax_cut = 0.52;
    Float_t k_nHitsMax_cut = 0.52;

    //PID cuts
    //Float_t k_nSigma_cut; //pi_nSigma set in production code
    Float_t k_nSigma_cut = 2;
    //pi_nSigma already set in production

    //Float_t pi1_TOFinvbeta_cut, pi2_TOFinvbeta_cut, k_TOFinvbeta_cut;
    Float_t pi1_TOFinvbeta_cut = 0.03;
    Float_t pi2_TOFinvbeta_cut = 0.03;
    Float_t k_TOFinvbeta_cut = 0.03;


    //topological production pre-cuts
    //will be used when ReadMode = 0

    //DCA to primary vertex of daughters
     Float_t k_dca_cut[nCentBins][nPtBins] = {
      {0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006}, //0-10%
      {0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006}, //10-40%
      {0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006}, //40-80%
      {0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006}  //0-80%
    };
    Float_t pi1_dca_cut[nCentBins][nPtBins] = {
      {0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006}, //0-10%
      {0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006}, //10-40%
      {0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006}, //40-80%
      {0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006}  //0-80%
    };
    Float_t pi2_dca_cut[nCentBins][nPtBins] = {
      {0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006}, //0-10%
      {0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006}, //10-40%
      {0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006}, //40-80%
      {0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006}  //0-80%
    };

    Float_t mdcaMax_cut[nCentBins][nPtBins] = {
      {0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011}, //0-10%
      {0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011}, //10-40%
      {0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011}, //40-80%
      {0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011}  //0-80%
    }; //DCA between K-pi (pi-pi) pairs
    Float_t D_decayL_min_cut[nCentBins][nPtBins] = {
      {0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011}, //0-10%
      {0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011}, //10-40%
      {0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011}, //40-80%
      {0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011}  //0-80%
    }; //D decay length
    Float_t D_decayL_max_cut = 100; //maximum is fixed

    Float_t D_cos_theta_cut[nCentBins][nPtBins] = {
      {0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995}, //0-10%
      {0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995}, //10-40%
      {0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995}, //40-80%
      {0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995, 0.995}  //0-80%
    }; //pointing angle
    Float_t D_dV0Max_cut[nCentBins][nPtBins] = {
      {100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100}, //0-10%
      {100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100}, //10-40%
      {100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100}, //40-80%
      {100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100}  //0-80%
    }; //Max distance between recostructed pairs vertecies

}
#endif
