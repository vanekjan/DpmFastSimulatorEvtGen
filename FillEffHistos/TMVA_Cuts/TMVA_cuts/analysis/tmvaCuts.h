#ifndef TMVACUTS_H
#define TMVACUTS_H

namespace tmvaCuts
{
    //int const totalNumberOfEvents = 90.e6; // Run14 dataset - probably not used
    int const nCentBins = 4;
    int   const nPtBins = 12;
    float const pT_bins[nPtBins+1] = { 1., 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10. }; //my pT bins

    //int   const nPtBins_TMVA = 12;
    //float const pT_bins_TMVA[nPtBins_TMVA+1] = { 1., 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10. }; //my pT bins, for analysis and some versions of TMVA

    const int nPtBins_TMVA = 6;
    float const pT_bins_TMVA[nPtBins_TMVA+1] = {0., 1., 2., 3., 5., 7., 10.}; //my TMVA pT bins, for D0 nSignal estimation

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


    //topological cuts from TMVA
    //see .../Analysis/Root/TMVA/Dpm/Cuts/"centrality"/weights/"file_name".xml

    //DCA to primary vertex of daughters (testing sample of cuts - optimized on central bin)
     Float_t k_dca_cut[nCentBins][nPtBins] = {
      {0.019, 0.012, 0.012, 0.0097, 0.0089, 0.0084, 0.0072, 0.0072, 0.0073, 0.0073, 0.0073, 0.0073}, //0-10%
      {0.019, 0.012, 0.012, 0.0097, 0.0089, 0.0084, 0.0072, 0.0072, 0.0073, 0.0073, 0.0073, 0.0073}, //10-40%
      {0.019, 0.012, 0.012, 0.0097, 0.0089, 0.0084, 0.0072, 0.0072, 0.0073, 0.0073, 0.0073, 0.0073}, //40-80%
      {0.019, 0.012, 0.012, 0.0097, 0.0089, 0.0084, 0.0072, 0.0072, 0.0073, 0.0073, 0.0073, 0.0073}  //0-80%
    };
    Float_t pi1_dca_cut[nCentBins][nPtBins] = {
      {0.018, 0.012, 0.012, 0.011, 0.0098, 0.0092, 0.0091, 0.0090, 0.0091, 0.0091, 0.0091, 0.0091}, //0-10%
      {0.018, 0.012, 0.012, 0.011, 0.0098, 0.0092, 0.0091, 0.0090, 0.0091, 0.0091, 0.0091, 0.0091}, //10-40%
      {0.018, 0.012, 0.012, 0.011, 0.0098, 0.0092, 0.0091, 0.0090, 0.0091, 0.0091, 0.0091, 0.0091}, //40-80%
      {0.018, 0.012, 0.012, 0.011, 0.0098, 0.0092, 0.0091, 0.0090, 0.0091, 0.0091, 0.0091, 0.0091}  //0-80%
    };
    Float_t pi2_dca_cut[nCentBins][nPtBins] = {
      {0.019, 0.013, 0.015, 0.011, 0.0098, 0.0093, 0.0090, 0.0092, 0.0093, 0.0093, 0.0093, 0.0093}, //0-10%
      {0.019, 0.013, 0.015, 0.011, 0.0098, 0.0093, 0.0090, 0.0092, 0.0093, 0.0093, 0.0093, 0.0093}, //10-40%
      {0.019, 0.013, 0.015, 0.011, 0.0098, 0.0093, 0.0090, 0.0092, 0.0093, 0.0093, 0.0093, 0.0093}, //40-80%
      {0.019, 0.013, 0.015, 0.011, 0.0098, 0.0093, 0.0090, 0.0092, 0.0093, 0.0093, 0.0093, 0.0093}  //0-80%
    };

    Float_t mdcaMax_cut[nCentBins][nPtBins] = {
      {0.0087, 0.0088, 0.0087, 0.0089, 0.0088, 0.0089, 0.0089, 0.0089, 0.0088, 0.0088, 0.0088, 0.0088}, //0-10%
      {0.0087, 0.0088, 0.0087, 0.0089, 0.0088, 0.0089, 0.0089, 0.0089, 0.0088, 0.0088, 0.0088, 0.0088}, //10-40%
      {0.0087, 0.0088, 0.0087, 0.0089, 0.0088, 0.0089, 0.0089, 0.0089, 0.0088, 0.0088, 0.0088, 0.0088}, //40-80%
      {0.0087, 0.0088, 0.0087, 0.0089, 0.0088, 0.0089, 0.0089, 0.0089, 0.0088, 0.0088, 0.0088, 0.0088}  //0-80%
    }; //DCA between K-pi (pi-pi) pairs
    Float_t D_decayL_min_cut[nCentBins][nPtBins] = {
      {0.032, 0.046, 0.041, 0.034, 0.042, 0.037, 0.042, 0.035, 0.034, 0.034, 0.034, 0.034}, //0-10%
      {0.032, 0.046, 0.041, 0.034, 0.042, 0.037, 0.042, 0.035, 0.034, 0.034, 0.034, 0.034}, //10-40%
      {0.032, 0.046, 0.041, 0.034, 0.042, 0.037, 0.042, 0.035, 0.034, 0.034, 0.034, 0.034}, //40-80%
      {0.032, 0.046, 0.041, 0.034, 0.042, 0.037, 0.042, 0.035, 0.034, 0.034, 0.034, 0.034}  //0-80%
    }; //D decay length
    Float_t D_decayL_max_cut = 0.2; //maximum is fixed

    Float_t D_cos_theta_cut[nCentBins][nPtBins] = {
      {0.9977, 0.9977, 0.9983, 0.9983, 0.998, 0.9986, 0.9986, 0.9983, 0.9976, 0.9976, 0.9976, 0.9976}, //0-10%
      {0.9977, 0.9977, 0.9983, 0.9983, 0.998, 0.9986, 0.9986, 0.9983, 0.9976, 0.9976, 0.9976, 0.9976}, //10-40%
      {0.9977, 0.9977, 0.9983, 0.9983, 0.998, 0.9986, 0.9986, 0.9983, 0.9976, 0.9976, 0.9976, 0.9976}, //40-80%
      {0.9977, 0.9977, 0.9983, 0.9983, 0.998, 0.9986, 0.9986, 0.9983, 0.9976, 0.9976, 0.9976, 0.9976}  //0-80%
    }; //pointing angle
    Float_t D_dV0Max_cut[nCentBins][nPtBins] = {
      {0.0057, 0.0074, 0.0105, 0.0119, 0.0137, 0.0163, 0.0203, 0.0210, 0.0212, 0.0212, 0.0212, 0.0212}, //0-10%
      {0.0057, 0.0074, 0.0105, 0.0119, 0.0137, 0.0163, 0.0203, 0.0210, 0.0212, 0.0212, 0.0212, 0.0212}, //10-40%
      {0.0057, 0.0074, 0.0105, 0.0119, 0.0137, 0.0163, 0.0203, 0.0210, 0.0212, 0.0212, 0.0212, 0.0212}, //40-80%
      {0.0057, 0.0074, 0.0105, 0.0119, 0.0137, 0.0163, 0.0203, 0.0210, 0.0212, 0.0212, 0.0212, 0.0212}  //0-80%
    }; //Max distance between recostructed pairs vertecies

}
#endif
