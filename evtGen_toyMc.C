/* *********************************************************************
 *  ROOT macro - Toy Monte Carlo Simulation for D0 decay
 *  Includes Momentum Resolution, DCA, hft ration, TPC efficiency ...
 *  Example for phi --> K+K-
 *
 *  Authors:
 *            Guannan Xie (guannanxie@lbl.gov)
 *            **Mustafa Mustafa (mmustafa@lbl.gov)
 *            Hao Qiu (hqiu@lbl.gov)
 *
 *  ** Code Maintainer
 *
 * *********************************************************************
 */

#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TF1.h"
#include "TClonesArray.h"
#include "TPythia6.h"
#include "TPythia6Decayer.h"
#include "TRandom3.h"
#include "TParticle.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "phys_constants.h"
#include "SystemOfUnits.h"
#include "TParticlePDG.h"

#include "EvtGen/EvtGen.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtHepMCEvent.hh"
#include "EvtGenBase/EvtSimpleRandomEngine.hh"
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenExternal/EvtExternalGenList.hh"
#include "StarEvtGenDecayer.h"

//for dca from helix
#include "StThreeVectorF.hh"
#include "StPhysicalHelixD.hh"

using namespace std;

//----------------------FROM PYTHIA FAST-SIM-------------------------------------------------------------------------------

void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& p1Mom, TLorentzVector const& p2Mom, TVector3 v00, int centLow, int centUp);

TLorentzVector smearMom(TLorentzVector const& b, TF1 const * const fMomResolution);
TVector3 smearPos(TLorentzVector const& mom, TLorentzVector const& rMom, TVector3 const& pos);
TVector3 smearPosData(int iParticleIndex, double vz, int cent, TLorentzVector const& rMom, TVector3 const& pos);
float dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaHelix(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex, float charge);
float dcaSigned(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaXY(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaZ(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dca1To2(TVector3 const& p1, TVector3 const& pos1, TVector3 const& p2, TVector3 const& pos2, TVector3& v0);
float dca1To2_helix(TVector3 const& p1, TVector3 const& pos1, float charge1, TVector3 const& p2, TVector3 const& pos2, float charge2, TVector3 const& PrimVertex, TVector3& v0);
TVector3 getVertex(int centrality);
bool matchHft(int iParticleIndex, double vz, int cent, TLorentzVector const& mom);
bool tpcReconstructed(int iParticleIndex, float charge, int cent, TLorentzVector const& mom);

int getPtIndexDca(double);
int getEtaIndexDca(double);
int getVzIndexDca(double);
int getPhiIndexDca(double);

int getPtIndexHftRatio(double);
int getEtaIndexHftRatio(double);
int getVzIndexHftRatio(double);
int getPhiIndexHftRatio(double);

TNtuple* nt;
TFile* result;

TF1* fKaonPlusMomResolution = NULL;
TF1* fKaonMinusMomResolution = NULL;
TF1* fPionPlusMomResolution = NULL;
TF1* fPionMinusMomResolution = NULL;

TF1* fWeightFunction = NULL;
const Int_t nParticles = 2;
const Int_t nCentHftRatio = 9;
int const centDef[nCentHftRatio+1] = {80, 70, 60, 50, 40, 30, 20, 10, 5, 0}; //centrality definition in %

// HFT ratio binning
const Int_t nEtasHftRatio = 10;
const Int_t nVzsHftRatio = 6;
const Int_t nPtBinsHftRatio = 36;
const Double_t EtaEdgeHftRatio[nEtasHftRatio + 1] =
{
  -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0
};
const Double_t VzEdgeHftRatio[nVzsHftRatio + 1] = {-6.0e4, -4.0e4, -2.0e4, 0.0, 2.0e4, 4.0e4, 6.0e4}; //vertex position in micrometers in simulation!!!
/*{

  -6., -4., -2.,  0., 2.,  4.,  6.

  };
  */
const Double_t ptEdgeHftRatio[nPtBinsHftRatio + 1] =
{
  0.3, 0.4, 0.5, 0.6 , 0.7 , 0.8 , 0.9 ,
  1. , 1.1 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.7 , 1.8 , 1.9 ,
  2. , 2.2 , 2.4 , 2.6 , 2.8 , 3.0 ,
  3.4 , 3.8 , 4.2 , 4.6 , 5.0 ,  5.5 ,
  6. , 6.5 , 7.0 , 8.0 , 9.0 , 10. , 11.,  12.0
};

const Int_t nPhisHftRatio = 11;
const Double_t PhiEdgeHftRatio[nPhisHftRatio + 1] =
{
  -3.14159 , -2.80359 , -2.17527 , -1.54696 , -0.918637 , -0.290319 , 0.338 , 0.966319 , 1.59464 , 2.22296 , 2.85127 , 3.14159 //Sector by Sector  // sector number 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3
};

const Int_t nPhisDca = 11;
const Double_t PhiEdgeDca[nPhisDca + 1] =
{
  -3.14159 , -2.80359 , -2.17527 , -1.54696 , -0.918637 , -0.290319 , 0.338 , 0.966319 , 1.59464 , 2.22296 , 2.85127 , 3.14159 //Sector by Sector  // sector number 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3
};

// DCA binning
int const nVzsDca = 4;
float const VzEdgeDca[nVzsDca + 1] = {   -6.e4, -3.e4, 0, 3.e4, 6.e4}; //vertex position in micrometers in simulation!!!

int const nEtasDca = 5;
float const EtaEdgeDca[nEtasDca + 1] = { -1.0, -0.6, -0.2, 0.2, 0.6, 1.0};

const Int_t nPtBinsDca = 19;
const Double_t ptEdgeDca[nPtBinsDca + 1] =
{
  0.3, 0.4, 0.5,
  0.6,  0.7 , 0.8 , 0.9 ,
  1. , 1.25 , 1.5 , 1.75 , 2.  , 2.25 , 2.5 , 2.75 , 3.0 , 3.5,
  4.  , 6. , 12.0
};

TH1D* h1Vz[nCentHftRatio];

TH1D* hHftRatio1[nParticles][nEtasHftRatio][nVzsHftRatio][nPhisHftRatio][nCentHftRatio];
TH1D* hHftRatioCorrect[nParticles][nEtasHftRatio][nVzsHftRatio][nPhisHftRatio][1];
int const nCentDca = 9;
TH2D* h2Dca[nParticles][nEtasDca][nVzsDca][nCentDca][nPtBinsDca];

TH1D* hTpcPiPlus[nCentHftRatio];
TH1D* hTpcPiMinus[nCentHftRatio];
TH1D* hTpcKPlus[nCentHftRatio];
TH1D* hTpcKMinus[nCentHftRatio];

//TF1 *InvEffWeight;

string outFileName = "Dpm.toyMc.root"; //default output file name (do not change - see submit XML)

std::pair<float, float> const momentumRange(0, 11); //momentum range
//std::pair<float, float> const momentumRange(0, 3); //for low-pT to increase statistics

float const gVzCut = 6.0e4;
float const acceptanceRapidity = 1.0;
float const sigmaPos0 = 15.2;
float const pxlLayer1Thickness = 0.00486;
float const sigmaVertexCent[nCentHftRatio] = {31., 18.1, 12.8, 9.3, 7.2, 5.9, 5., 4.6, 4.};

TH1D* h_pi_tof_eff;
TH1D* h_k_tof_eff;
//-------------------------------------------------------------------------------------------------------------------------

void initEvtGen();
void decayAndFill(TString name, TLorentzVector* b, double const weight, TClonesArray& daughters); //Decay particle - name
void decayAndFill(int PDG_id, TLorentzVector* b, double const weight, TClonesArray& daughters, int centLow, int centUp); //Decay particle - PDG ID
void getKinematics(TLorentzVector& b, double const mass, const int pTspectrum);
void bookObjects(int centLow, int centUp);
void write();

StarEvtGenDecayer* starEvtGenDecayer = NULL;

//============== main  program ==================

//centrality range in %, default range is 0% to 80%, pTspectrum = 0 - flat pT, 1 - Levy (weight function)
//specifically set in run_evtGent_toyMc.C
void evtGen_toyMc(int npart = 1000, int centLow = 0, int centUp = 80, int pTspectrum = 0) 
{
  cout<<"Starting EvtGen"<<endl;
  initEvtGen();
  cout<<"initEvtGen() done..."<<endl;

  gRandom->SetSeed();
  bookObjects(centLow, centUp);

  TClonesArray ptl("TParticle", 10); //array of 10 TParticles

  for (int ipart = 0; ipart < npart; ipart++)
  {
    if (!(ipart % 100000))
      cout << "____________ ipart = " << ipart / static_cast<float>(npart) << " ________________" << endl;

    TLorentzVector* b_d = new TLorentzVector;
    TLorentzVector* b_d2 = new TLorentzVector;

    getKinematics(*b_d, M_D_PLUS, pTspectrum);//D+ (PDG - 411) - generates kinematics of D+
    getKinematics(*b_d2, M_D_PLUS, pTspectrum);//D+ (PDG - 411) - generates kinematics of D-

    decayAndFill(411, b_d, fWeightFunction->Eval(b_d->Perp()), ptl, centLow, centUp);//D+
    decayAndFill(-411, b_d2, fWeightFunction->Eval(b_d2->Perp()), ptl, centLow, centUp);//D-

  }


  cout<<"Write!"<<endl;
  write();
  cout<<"Written!"<<endl;



  return;

}

//not used in this version of fast-sim, see the next decayAndFill(...)
void decayAndFill(TString name, TLorentzVector* b, double const weight, TClonesArray& daughters)//decay for lambda(1520)
{
  TLorentzVector kMom;
  TLorentzVector piMom;
  TLorentzVector pMom;

  cout << "check point here A " << endl;

  starEvtGenDecayer->Decay(name, b); //Decay works also with PDG ID of the particle

  cout << "check point here B " << endl;

  starEvtGenDecayer->ImportParticles(&daughters); //get daughters from decay

  cout << "check point here C " << endl;

  int nTrk = daughters.GetEntriesFast(); //get number of generated tracks (all daughters)
  cout << "nTrk = " << nTrk << endl;

  for (int iTrk = 0; iTrk < nTrk; ++iTrk) //go through all tracks
  {
    TParticle* ptl0 = (TParticle*)daughters.At(iTrk);

    cout << "iTrk = " << iTrk << " , " << ptl0->GetPdgCode() << endl;

    if (abs(ptl0->GetPdgCode()) == 211) //save momentum of pion, kaon and proton
    {
      // case 211://pion
      ptl0->Momentum(piMom);
    }
    if (abs(ptl0->GetPdgCode()) == 321)
    {
      // case 321://kaon
      ptl0->Momentum(kMom);
    }
    if (abs(ptl0->GetPdgCode()) == 2212)
    {
      // case 2212://proton
      ptl0->Momentum(pMom);
    }
  }

  daughters.Clear();

  //fill(kf, b, weight, kMom, p1Mom, p2Mom, v00);

}

void decayAndFill(int PDG_id, TLorentzVector* b, double const weight, TClonesArray& daughters, int centLow, int centUp)//decay for D+
{

  starEvtGenDecayer->Decay(PDG_id, b); //Decay particle

  starEvtGenDecayer->ImportParticles(&daughters); //get daughters from decay

  TLorentzVector kMom;
  TLorentzVector p1Mom;
  TLorentzVector p2Mom;
  TVector3 v00;

  int nTrk = daughters.GetEntriesFast();
  for (int iTrk = 0; iTrk < nTrk; ++iTrk)
  {
    TParticle* ptl0 = (TParticle*)daughters.At(iTrk);

    switch (abs(ptl0->GetPdgCode()))
    {
      case 321:
        ptl0->Momentum(kMom);
        v00.SetXYZ(ptl0->Vx() * 1000., ptl0->Vy() * 1000., ptl0->Vz() * 1000.); // converted to μm
        break;
      case 211:
        if (!p1Mom.P()) ptl0->Momentum(p1Mom);
        else ptl0->Momentum(p2Mom);
        break;
      default:
        break;
    }
  }

  daughters.Clear();

  fill(PDG_id, b, weight, kMom, p1Mom, p2Mom, v00, centLow, centUp);

}

void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& p1Mom, TLorentzVector const& p2Mom, TVector3 v00, int centLow, int centUp)
{
  int const centrality_rndm = gRandom->Uniform(centLow, centUp); //generate random centrality in % - new version with possibility to set any centrality range

  int centrality = -1; //defaul value

  for(unsigned int i = 0; i < nCentHftRatio; i++) //determine centrality bin from % (different centrality bin width for most central collisions)
  {
    if(centrality_rndm < centDef[i] && centrality_rndm >= centDef[i+1])
    {
      centrality = i;
    }
  }

  float const MCdecayLength = v00.Mag(); //generated decay length, from generator, before smearing

  TVector3 const vertex = getVertex(centrality);

  v00 += vertex; //shift MC secondary vertex along the beam by v_z

  // smear momentum
  TLorentzVector kRMom;
  TLorentzVector p1RMom;
  TLorentzVector p2RMom;

  if(kf > 0) //D+
  {
    kRMom  = smearMom(kMom,  fKaonMinusMomResolution);
    p1RMom = smearMom(p1Mom, fPionPlusMomResolution);
    p2RMom = smearMom(p2Mom, fPionPlusMomResolution);
  }
  else //D-
  {
    kRMom  = smearMom(kMom,  fKaonPlusMomResolution);
    p1RMom = smearMom(p1Mom, fPionMinusMomResolution);
    p2RMom = smearMom(p2Mom, fPionMinusMomResolution);
  }


  // smear position
  TVector3 const kRPos  = smearPosData(1, vertex.z(), centrality, kRMom,  v00);
  TVector3 const p1RPos = smearPosData(0, vertex.z(), centrality, p1RMom, v00);
  TVector3 const p2RPos = smearPosData(0, vertex.z(), centrality, p2RMom, v00);


  // reconstruct
  int const charge = kf > 0 ? 1 : -1;

  TLorentzVector const rMom = kRMom + p1RMom + p2RMom;
  float const kDca = dca(kMom.Vect(), v00, vertex);
  float const p1Dca = dca(p1Mom.Vect(), v00, vertex);
  float const p2Dca = dca(p2Mom.Vect(), v00, vertex);

  float const kRDca = dca(kRMom.Vect(), kRPos, vertex);
  float const kRDca_helix = dcaHelix(kRMom.Vect(), kRPos, vertex, -1*charge); //dca from helix, kaon has opposite charge than D meson
  float const kRSDca = dcaSigned(kRMom.Vect(), kRPos, vertex);
  float const kRDcaXY = dcaXY(kRMom.Vect(), kRPos, vertex);
  float const kRDcaZ = dcaZ(kRMom.Vect(), kRPos, vertex);
  float const kRDcaXYtoSV = dcaXY(kRMom.Vect(), kRPos, v00);
  float const kRDcaZtoSV = dcaZ(kRMom.Vect(), kRPos, v00);

  float const p1RDca = dca(p1RMom.Vect(), p1RPos, vertex);
  float const p1RDca_helix = dcaHelix(p1RMom.Vect(), p1RPos, vertex, charge); //dca from helix, pion has the same charge as D meson
  float const p1RSDca = dcaSigned(p1RMom.Vect(), p1RPos, vertex);
  float const p1RDcaXY = dcaXY(p1RMom.Vect(), p1RPos, vertex);
  float const p1RDcaZ = dcaZ(p1RMom.Vect(), p1RPos, vertex);
  float const p1RDcaXYtoSV = dcaXY(p1RMom.Vect(), p1RPos, v00);
  float const p1RDcaZtoSV = dcaZ(p1RMom.Vect(), p1RPos, v00);

  float const p2RDca = dca(p2RMom.Vect(), p2RPos, vertex);
  float const p2RDca_helix = dcaHelix(p2RMom.Vect(), p2RPos, vertex, charge); //dca from helix, pion has the same charge as D meson
  float const p2RSDca = dcaSigned(p2RMom.Vect(), p2RPos, vertex);
  float const p2RDcaXY = dcaXY(p2RMom.Vect(), p2RPos, vertex);
  float const p2RDcaZ = dcaZ(p2RMom.Vect(), p2RPos, vertex);
  float const p2RDcaXYtoSV = dcaXY(p2RMom.Vect(), p2RPos, v00);
  float const p2RDcaZtoSV = dcaZ(p2RMom.Vect(), p2RPos, v00);

  TVector3 v12, v23, v31, v0;
  float const dca12 = dca1To2(kRMom.Vect(), kRPos, p1RMom.Vect(), p1RPos, v12);
  float const dca23 = dca1To2(p1RMom.Vect(), p1RPos, p2RMom.Vect(), p2RPos, v23);
  float const dca31 = dca1To2(p2RMom.Vect(), p2RPos, kRMom.Vect(), kRPos, v31);

  float dcaDaughters = dca12 > dca31 ? dca12 : dca31;
  dcaDaughters = dcaDaughters > dca23 ? dcaDaughters : dca23;
  v0 = (v12 + v23 + v31) * 0.333333333;

  float const decayLength = (v0 - vertex).Mag();
  float const dcaDpmToPv = dca(rMom.Vect(), v0, vertex);
  float const cosTheta = (v0 - vertex).Unit().Dot(rMom.Vect().Unit());
  float const angle12 = kRMom.Vect().Angle(p1RMom.Vect());

  // Distance between v12 and v23
  float const v12_23 = (v12 - v23).Mag();
  // Distance between v23 and v31
  float const v23_31 = (v23 - v31).Mag();
  // Distance between v31 and v12
  float const v31_12 = (v31 - v12).Mag();

  //maximum dist between reo v0's to be averaging
  float const max12 =  v12_23 > v23_31 ? v12_23 : v23_31 ;
  float const mdV0Max = max12 > v31_12 ? max12 : v31_12;

  //calculate also from helix
  TVector3 v12_helix, v23_helix, v31_helix, v0_helix;
  float const dca12_helix = dca1To2_helix(kRMom.Vect(), kRPos, -1*charge, p1RMom.Vect(), p1RPos, charge, vertex, v12_helix);
  float const dca23_helix = dca1To2_helix(p1RMom.Vect(), p1RPos, charge, p2RMom.Vect(), p2RPos, charge, vertex, v23_helix);
  float const dca31_helix = dca1To2_helix(p2RMom.Vect(), p2RPos, charge, kRMom.Vect(), kRPos, -1*charge, vertex, v31_helix);

  float dcaDaughters_helix = dca12_helix > dca31_helix ? dca12_helix : dca31_helix;
  dcaDaughters_helix = dcaDaughters_helix > dca23_helix ? dcaDaughters_helix : dca23_helix;
  v0_helix = (v12_helix + v23_helix + v31_helix) * 0.333333333;

  float const decayLength_helix = (v0_helix - vertex).Mag();
  float const dcaDpmToPv_helix = dcaHelix(rMom.Vect(), v0_helix, vertex, charge);
  float const cosTheta_helix = (v0_helix - vertex).Unit().Dot(rMom.Vect().Unit());

  // Distance between v12 and v23
  float const v12_23_helix = (v12_helix - v23_helix).Mag();
  // Distance between v23 and v31
  float const v23_31_helix = (v23_helix - v31_helix).Mag();
  // Distance between v31 and v12
  float const v31_12_helix = (v31_helix - v12_helix).Mag();

  //maximum dist between reo v0's to be averaging
  float const max12_helix =  v12_23_helix > v23_31_helix ? v12_23_helix : v23_31_helix ;
  float const mdV0Max_helix = max12_helix > v31_12_helix ? max12_helix : v31_12_helix;


  TLorentzVector kRMomRest = kRMom;
  TVector3 beta;
  beta.SetMagThetaPhi(rMom.Beta(), rMom.Theta(), rMom.Phi());
  kRMomRest.Boost(-beta);
  float const cosThetaStar = rMom.Vect().Unit().Dot(kRMomRest.Vect().Unit());


  // save
  float arr[110];
  int iArr = 0;
  arr[iArr++] = centrality;
  arr[iArr++] = vertex.X();
  arr[iArr++] = vertex.Y();
  arr[iArr++] = vertex.Z();
  arr[iArr++] = getVzIndexDca(vertex.Z());

  arr[iArr++] = kf;
  arr[iArr++] = weight;
  arr[iArr++] = b->M();
  arr[iArr++] = b->Perp();
  arr[iArr++] = b->PseudoRapidity();
  arr[iArr++] = b->Rapidity();
  arr[iArr++] = b->Phi();
  arr[iArr++] = v00.X();
  arr[iArr++] = v00.Y();
  arr[iArr++] = v00.Z();

  arr[iArr++] = rMom.M();
  arr[iArr++] = rMom.Perp();
  arr[iArr++] = rMom.PseudoRapidity();
  arr[iArr++] = rMom.Rapidity();
  arr[iArr++] = rMom.Phi();
  arr[iArr++] = v0.X();
  arr[iArr++] = v0.Y();
  arr[iArr++] = v0.Z();

  arr[iArr++] = dcaDaughters;
  arr[iArr++] = dcaDaughters_helix;
  arr[iArr++] = decayLength; //reconstrucetd decay length
  arr[iArr++] = decayLength_helix;
  arr[iArr++] = MCdecayLength; //MC decay length
  arr[iArr++] = dcaDpmToPv;
  arr[iArr++] = dcaDpmToPv_helix;
  arr[iArr++] = cosTheta;
  arr[iArr++] = cosTheta_helix;

  arr[iArr++] = kMom.M();
  arr[iArr++] = kMom.Perp();
  arr[iArr++] = kMom.PseudoRapidity();
  arr[iArr++] = kMom.Rapidity();
  arr[iArr++] = kMom.Phi();
  arr[iArr++] = kDca;

  arr[iArr++] = kRMom.M();
  arr[iArr++] = kRMom.Perp();
  arr[iArr++] = kRMom.PseudoRapidity();
  arr[iArr++] = kRMom.Rapidity();
  arr[iArr++] = kRMom.Phi();
  arr[iArr++] = kRPos.X();
  arr[iArr++] = kRPos.Y();
  arr[iArr++] = kRPos.Z();
  arr[iArr++] = kRDca;
  arr[iArr++] = kRDca_helix;
  arr[iArr++] = kRDcaXY;
  arr[iArr++] = kRDcaZ;
  arr[iArr++] = kRDcaXYtoSV;
  arr[iArr++] = kRDcaZtoSV;
  arr[iArr++] = tpcReconstructed(1, -1 * charge, centrality, kRMom);

  arr[iArr++] = p1Mom.M();
  arr[iArr++] = p1Mom.Perp();
  arr[iArr++] = p1Mom.PseudoRapidity();
  arr[iArr++] = p1Mom.Rapidity();
  arr[iArr++] = p1Mom.Phi();
  arr[iArr++] = p1Dca;

  arr[iArr++] = p1RMom.M();
  arr[iArr++] = p1RMom.Perp();
  arr[iArr++] = p1RMom.PseudoRapidity();
  arr[iArr++] = p1RMom.Rapidity();
  arr[iArr++] = p1RMom.Phi();
  arr[iArr++] = p1RPos.X();
  arr[iArr++] = p1RPos.Y();
  arr[iArr++] = p1RPos.Z();
  arr[iArr++] = p1RDca;
  arr[iArr++] = p1RDca_helix;
  arr[iArr++] = p1RDcaXY;
  arr[iArr++] = p1RDcaZ;
  arr[iArr++] = p1RDcaXYtoSV;
  arr[iArr++] = p1RDcaZtoSV;
  arr[iArr++] = tpcReconstructed(0, charge, centrality, p1RMom);

  arr[iArr++] = p2Mom.M();
  arr[iArr++] = p2Mom.Perp();
  arr[iArr++] = p2Mom.PseudoRapidity();
  arr[iArr++] = p2Mom.Rapidity();
  arr[iArr++] = p2Mom.Phi();
  arr[iArr++] = p2Dca;

  arr[iArr++] = p2RMom.M();
  arr[iArr++] = p2RMom.Perp();
  arr[iArr++] = p2RMom.PseudoRapidity();
  arr[iArr++] = p2RMom.Rapidity();
  arr[iArr++] = p2RMom.Phi();
  arr[iArr++] = p2RPos.X();
  arr[iArr++] = p2RPos.Y();
  arr[iArr++] = p2RPos.Z();
  arr[iArr++] = p2RDca;
  arr[iArr++] = p2RDca_helix;
  arr[iArr++] = p2RDcaXY;
  arr[iArr++] = p2RDcaZ;
  arr[iArr++] = p2RDcaXYtoSV;
  arr[iArr++] = p2RDcaZtoSV;
  arr[iArr++] = tpcReconstructed(0, charge, centrality, p2RMom);

  arr[iArr++] = matchHft(1, vertex.z(), centrality, kRMom);
  arr[iArr++] = matchHft(0, vertex.z(), centrality, p1RMom);
  arr[iArr++] = matchHft(0, vertex.z(), centrality, p2RMom);

  arr[iArr++] = mdV0Max;
  arr[iArr++] = mdV0Max_helix;

  nt->Fill(arr);
}


void getKinematics(TLorentzVector& b, double const mass, const int pTspectrum)
{
  float pt = 0;
  if(pTspectrum == 0) pt = gRandom->Uniform(momentumRange.first, momentumRange.second); //flat pT distribution
  if(pTspectrum == 1) pt = fWeightFunction->GetRandom(momentumRange.first, momentumRange.second); //realistic pT distribution

  float const y = gRandom->Uniform(-acceptanceRapidity, acceptanceRapidity);
  float const phi = TMath::TwoPi() * gRandom->Rndm();

  float const mT = sqrt(mass * mass + pt * pt);
  float const pz = mT * sinh(y);
  float const E = mT * cosh(y);

  b.SetPxPyPzE(pt * cos(phi), pt * sin(phi) , pz, E);
}

//___________
void bookObjects(int centLow, int centUp)
{

  int centrality9_low = -1; //defaul value
  int centrality9_up = -1; //defaul value

  //determine centrality bin from % (different centrality bin width for most central collisions) 
  //the int centrality bins go from peripheral to most central - the up/low are "reversed"
  for(unsigned int i = 0; i < nCentHftRatio+1; i++)
  {
    if(centLow == centDef[i] )
    {
      centrality9_up = i;
      cout<<centrality9_up<<endl;
    }

    if(centUp == centDef[i] )
    {
      centrality9_low = i;
      cout<<centrality9_low<<endl;
    }
  }



  cout << "Loading input files ..." << endl;
  cout<<endl;

  cout << "Loading input momentum resolution ..." << endl;
  TFile f("./input/Momentum_resolution_Run16.root");
  fPionPlusMomResolution = (TF1*)f.Get("PiPlusMomResFit")->Clone("PiPlusMomResFit");
  fPionMinusMomResolution = (TF1*)f.Get("PiMinusMomResFit")->Clone("PiMinusMomResFit");
  fKaonPlusMomResolution = (TF1*)f.Get("KPlusMomResFit")->Clone("KPlusMomResFit");
  fKaonMinusMomResolution = (TF1*)f.Get("KMinusMomResFit")->Clone("KMinusMomResFit");
  f.Close();


  cout << "Loading input spectra ..." << endl;
  TFile fPP("./input/pp200_spectra.root");
  fWeightFunction = (TF1*)fPP.Get("run12/f1Levy")->Clone("f1Levy");
  fPP.Close();

  TFile fVertex("./input/Vz_Cent_Run16.root");

  for(int ii = centrality9_low; ii < centrality9_up; ii++)
  {
    h1Vz[ii] = (TH1D*)(fVertex.Get(Form("mh1VzWg_%i", ii)));
    h1Vz[ii]->SetDirectory(0);
  }

  fVertex.Close();

  cout << "Loading input HFT ratios and DCA ..." << endl;


  TFile *fHftRatio1 = new TFile("./input/HFT_Matching_Ratios_Run16.root", "read");
  TFile *fDca1 = new TFile("./input/DCA_Resolution_Run16.root", "read");


  for (int iParticle = 0; iParticle < nParticles; ++iParticle)
  {
    for(int iCent = centrality9_low; iCent < centrality9_up; iCent++)
    {
      // HFT ratio
      for (int iEta = 0; iEta < nEtasHftRatio; ++iEta)
      {
        for (int iVz = 0; iVz < nVzsHftRatio; ++iVz)
        {
          for (int iPhi = 0; iPhi < nPhisHftRatio; ++iPhi)
          {
            hHftRatio1[iParticle][iEta][iVz][iPhi][iCent] = (TH1D*)fHftRatio1->Get(Form("mh1HFT1PtCentPartEtaVzPhi_%i_%i_%i_%i_%i", iParticle, iEta, iVz, iPhi, iCent)); //Run16 names
            //hHftRatio1[iParticle][iEta][iVz][iPhi][iCent] = (TH1D*)fHftRatio1->Get(Form("mh1HFT1PtCentPartEtaVzPhiRatio_%i_%i_%i_%i_%i", iParticle, iEta, iVz, iPhi, iCent)); //Run14 names
            hHftRatio1[iParticle][iEta][iVz][iPhi][iCent]->SetDirectory(0);
          }
        }
      }
    }
    cout << "Finished loading HFT Ratio: " <<  endl;


    for(int iCent = centrality9_low; iCent < centrality9_up; iCent++)
    {
      // DCA
      for (int iEta = 0; iEta < nEtasDca; ++iEta)
      {
        for (int iVz = 0; iVz < nVzsDca; ++iVz)
        {
          for (int iPt = 0; iPt < nPtBinsDca; ++iPt)
          {
            h2Dca[iParticle][iEta][iVz][iCent][iPt] = (TH2D*)((fDca1->Get(Form("mh2DcaPtCentPartEtaVzPhi_%i_%i_%i_%i_%i", iParticle, iEta, iVz, iCent, iPt))));
            h2Dca[iParticle][iEta][iVz][iCent][iPt]->SetDirectory(0);
          }
        }
      }
    }    
  }
  cout << "Finished loading Dca: " <<  endl;


  std::cout << "Loading HFT ratios Correction factor for inclusive and primary..." << std::endl;
  TFile *fHftRatioCorrect = new TFile("./input/HFT_Ratio_Correction_Hijing_QM.root", "read");
  for (int iParticle = 0; iParticle < nParticles; ++iParticle)
  {
    for(int iCent = 0; iCent < 1; iCent++)
    {

      for (int iEta = 0; iEta < nEtasHftRatio; ++iEta)
      {
        for (int iVz = 0; iVz < nVzsHftRatio; ++iVz)
        {
          for (int iPhi = 0; iPhi < nPhisHftRatio; ++iPhi)
          {
            hHftRatioCorrect[iParticle][iEta][iVz][iPhi][iCent]  = (TH1D*)fHftRatioCorrect->Get(Form("mhHFTRatio_%i_%i_%i_%i_%i", iParticle, iEta, iVz, iPhi, iCent));
            hHftRatioCorrect[iParticle][iEta][iVz][iPhi][iCent]->SetDirectory(0);
          }
        }
      }
    }
  }
    
  fHftRatioCorrect->Close();

  fHftRatio1->Close();
  fDca1->Close();

  cout << " Loading TPC tracking efficiencies " << endl;

  TFile fTpcPiPlus("./input/Eff_PionPlus_embedding.root"); //all new for Run16, SL16j
  TFile fTpcPiMinus("./input/Eff_PionMinus_embedding.root");
  TFile fTpcKPlus("./input/Eff_KaonPlus_embedding.root");
  TFile fTpcKMinus("./input/Eff_KaonMinus_embedding.root");


  for(int iCent = centrality9_low; iCent < centrality9_up; iCent++)
  {
    hTpcPiPlus[iCent] = (TH1D*)fTpcPiPlus.Get(Form("TrackEffCent%i", iCent));
    hTpcPiPlus[iCent]->SetDirectory(0);
    hTpcPiMinus[iCent] = (TH1D*)fTpcPiMinus.Get(Form("TrackEffCent%i", iCent));
    hTpcPiMinus[iCent] ->SetDirectory(0);
    hTpcKPlus[iCent] = (TH1D*)fTpcKPlus.Get(Form("TrackEffCent%i", iCent));
    hTpcKPlus[iCent]->SetDirectory(0);
    hTpcKMinus[iCent] = (TH1D*)fTpcKMinus.Get(Form("TrackEffCent%i", iCent));
    hTpcKMinus[iCent]->SetDirectory(0);
  }

  fTpcPiPlus.Close();
  fTpcPiMinus.Close();
  fTpcKPlus.Close();
  fTpcKMinus.Close();

  result = new TFile(outFileName.c_str(), "recreate");
  result->SetCompressionLevel(1);
  result->cd();

  int BufSize = (int)pow(2., 16.);

  nt = new TNtuple("nt", "", "cent:vx:vy:vz:vzIdx:"
      "pid:w:m:pt:eta:y:phi:v0x:v0y:v0z:" // MC Dpm
      "rM:rPt:rEta:rY:rPhi:rV0x:rV0y:rV0z:" // Rc Dpm
      "dcaDaughters:dcaDaughters_helix:decayLength:decayLength_helix:MCdecayLength:dcaDpmToPv:dcaDpmToPv_helix:cosTheta:cosTheta_helix:" // Rc pair
      "kM:kPt:kEta:kY:kPhi:kDca:" // MC Kaon
      "kRM:kRPt:kREta:kRY:kRPhi:kRVx:kRVy:kRVz:kRDca:kRDca_helix:kRDcaXY:kRDcaZ:kRDcaXYtoSV:kRDcaZtoSV:kTpc:" // Rc Kaon
      "p1M:p1Pt:p1Eta:p1Y:p1Phi:p1Dca:" // MC Pion1
      "p1RM:p1RPt:p1REta:p1RY:p1RPhi:p1RVx:p1RVy:p1RVz:p1RDca:p1RDca_helix:p1RDcaXY:p1RDcaZ:p1RDcaXYtoSV:p1RDcaZtoSV:p1Tpc:" // Rc Pion1
      "p2M:p2Pt:p2Eta:p2Y:p2Phi:p2Dca:" // MC Pion2
      "p2RM:p2RPt:p2REta:p2RY:p2RPhi:p2RVx:p2RVy:p2RVz:p2RDca:p2RDca_helix:p2RDcaXY:p2RDcaZ:p2RDcaXYtoSV:p2RDcaZtoSV:p2Tpc:" // Rc Pion2
      "kHft:p1Hft:p2Hft:mdV0Max:mdV0Max_helix");

  cout << "Done with loading all files ..." << endl;
}
//_______________________________________________________________________________________________________________________


float dca(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
  TVector3 posDiff = pos - vertex;
  return fabs(p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff));
}

float dcaHelix(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex, float charge) //DCA to PV - calcualted same way as in StPicoTrack
{
  //all position vectors used for calcualtion with helix transformed form mum to cm

  StThreeVectorF p_work(p.X(), p.Y(), p.Z()); //create StThreeVectorF of momentum
  StThreeVectorF pos_work(pos.X()/1e4, pos.Y()/1e4, pos.Z()/1e4); //create StThreeVectorF of position

  StThreeVectorF vertex_work(vertex.X()/1e4, vertex.Y()/1e4, vertex.Z()/1e4); //create StThreeVectorF of PV


  StPhysicalHelixD helix_work(p_work, pos_work, -5*kilogauss , charge); //true helix with origin near SV

  helix_work.moveOrigin(helix_work.pathLength(vertex_work)); // move origin to DCA to PV

  StThreeVectorF origin = helix_work.origin(); //new helix origin at DCA


  return (origin - vertex_work).mag()*1e4;

}

float dcaSigned(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
  TVector3 posDiff = pos - vertex;
  float sign = posDiff.x() * p.y() - posDiff.y() * p.x() > 0 ? +1 : -1;

  return sign * p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff);
}

float dcaXY(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
  TVector3 newPos(pos);
  newPos.SetZ(0);

  TVector3 newP(p);
  newP.SetZ(0);

  TVector3 newVertex(vertex);
  newVertex.SetZ(0);

  TVector3 posDiff = newPos - newVertex;
  float sign = posDiff.x() * p.y() - posDiff.y() * p.x() > 0 ? +1 : -1;
  return sign * newP.Cross(posDiff.Cross(newP)).Unit().Dot(posDiff);
}

float dcaZ(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex)
{
  TVector3 posDiff = pos - vertex;
  if (sin(p.Theta()) == 0) return 0;
  else return (-posDiff.x() * cos(p.Phi()) - posDiff.y() * sin(p.Phi())) * cos(p.Theta()) / sin(p.Theta()) + posDiff.z();
}

float dca1To2(TVector3 const& p1, TVector3 const& pos1, TVector3 const& p2, TVector3 const& pos2, TVector3& v0)
{
  TVector3 posDiff = pos2 - pos1;
  TVector3 pu1 = p1.Unit();
  TVector3 pu2 = p2.Unit();
  double pu1Pu2 = pu1.Dot(pu2);
  double g = posDiff.Dot(pu1);
  double k = posDiff.Dot(pu2);
  double s2 = (k - pu1Pu2 * g) / (pu1Pu2 * pu1Pu2 - 1.);
  double s1 = g + s2 * pu1Pu2;
  TVector3 posDca1 = pos1 + pu1 * s1;
  TVector3 posDca2 = pos2 + pu2 * s2;
  v0 = 0.5 * (posDca1 + posDca2);
  return (posDca1 - posDca2).Mag();
}

float dca1To2_helix(TVector3 const& p1, TVector3 const& pos1, float charge1, TVector3 const& p2, TVector3 const& pos2, float charge2, TVector3 const& PrimVertex, TVector3& v0)
{
  //all position vectors used for calcualtion with helix transformed form mum to cm

  StThreeVectorF PrimVertex_work(PrimVertex.X()/1e4, PrimVertex.Y()/1e4, PrimVertex.Z()/1e4); //primary vertex needed to move origin of helixes

  //particle 1  
  StThreeVectorF p1_work(p1.X(), p1.Y(), p1.Z()); //create StThreeVectorF of momentum
  StThreeVectorF pos1_work(pos1.X()/1e4, pos1.Y()/1e4, pos1.Z()/1e4); //create StThreeVectorF of position

  StPhysicalHelixD helix1_work(p1_work, pos1_work, -5*kilogauss , charge1); //true helix with ture origin near SV
  helix1_work.moveOrigin(helix1_work.pathLength(PrimVertex_work)); //move helix origin to DCA to PV

  StThreeVectorF p1_at_DCA_work = helix1_work.momentum(-5*kilogauss); // new momentum at DCA  
  StPhysicalHelixD helix1_line_work(p1_at_DCA_work, helix1_work.origin(), 0, charge1); //use straight lines to calculate pair DCA (same as in data)
  //______________________________________________________________________________________________________________________
  //particle 2
  StThreeVectorF p2_work(p2.X(), p2.Y(), p2.Z()); //create StThreeVectorF of momentum
  StThreeVectorF pos2_work(pos2.X()/1e4, pos2.Y()/1e4, pos2.Z()/1e4); //create StThreeVectorF of position

  StPhysicalHelixD helix2_work(p2_work, pos2_work, -5*kilogauss , charge2); //true helix with ture origin near S
  helix2_work.moveOrigin(helix2_work.pathLength(PrimVertex_work)); //move helix origin to DCA to PV

  StThreeVectorF p2_at_DCA_work = helix2_work.momentum(-5*kilogauss); // new momentum at DCA
  StPhysicalHelixD helix2_line_work(p2_at_DCA_work, helix2_work.origin(), 0, charge2); //use straight lines to calculate pair DCA (same as in data)
  //______________________________________________________________________________________________________________________

  //calculate pair DCA
  pair<float,float> const ss = helix1_line_work.pathLengths(helix2_line_work);
  StThreeVectorF posDca1_work = helix1_line_work.at(ss.first);
  StThreeVectorF posDca2_work = helix2_line_work.at(ss.second); 

  TVector3 posDca1(posDca1_work.x()*1e4, posDca1_work.y()*1e4, posDca1_work.z()*1e4);
  TVector3 posDca2(posDca2_work.x()*1e4, posDca2_work.y()*1e4, posDca2_work.z()*1e4);

  v0 = 0.5 * (posDca1 + posDca2); //secondary vertex of the pair
  return (posDca1 - posDca2).Mag(); //pair DCA
}

TLorentzVector smearMom(TLorentzVector const& b, TF1 const * const fMomResolution)
{
  float const pt = b.Perp();
  float const sPt = gRandom->Gaus(pt, pt * fMomResolution->Eval(pt));

  TLorentzVector sMom;
  sMom.SetXYZM(sPt * cos(b.Phi()), sPt * sin(b.Phi()), sPt * sinh(b.PseudoRapidity()), b.M());
  return sMom;
}

TVector3 smearPos(TLorentzVector const& mom, TLorentzVector const& rMom, TVector3 const& pos)
{
  float thetaMCS = 13.6 / mom.Beta() / rMom.P() / 1000 * sqrt(pxlLayer1Thickness / fabs(sin(mom.Theta())));
  float sigmaMCS = thetaMCS * 28000 / fabs(sin(mom.Theta()));
  float sigmaPos = sqrt(pow(sigmaMCS, 2) + pow(sigmaPos0, 2));

  return TVector3(gRandom->Gaus(pos.X(), sigmaPos), gRandom->Gaus(pos.Y(), sigmaPos), gRandom->Gaus(pos.Z(), sigmaPos));
}

int getPtIndexDca(double pT)
{
  for (int i = 0; i < nPtBinsDca; i++)
  {
    if ((pT >= ptEdgeDca[i]) && (pT < ptEdgeDca[i + 1]))
      return i;
  }

  if(pT > 0 && pT < ptEdgeDca[0] ) return 0; //if pT is smaller than lower edge of first bin (i = 0) - use first bin as aporximation

  return nPtBinsDca - 1 ; //in all other cases pT is larger than upper edge of last bin - use last bin as approximation

}

int getEtaIndexDca(double Eta)
{
  for (int i = 0; i < nEtasDca; i++)
  {
    if ((Eta >= EtaEdgeDca[i]) && (Eta < EtaEdgeDca[i + 1]))
      return i;
  }
  return nEtasDca - 1 ;
}

int getVzIndexDca(double Vz)
{
  for (int i = 0; i < nVzsDca; i++)
  {
    if ((Vz >= VzEdgeDca[i]) && (Vz < VzEdgeDca[i + 1]))
      return i;
  }
  return nVzsDca - 1 ;
}

int getPhiIndexDca(double Phi)
{
  for (int i = 0; i < nPhisDca; i++)
  {
    if ((Phi >= PhiEdgeDca[i]) && (Phi < PhiEdgeDca[i + 1]))
      return i;
  }
  return nPhisDca - 1 ;
}

int getPtIndexHftRatio(double pT)
{
  for (int i = 0; i < nPtBinsHftRatio; i++)
  {
    if ((pT >= ptEdgeHftRatio[i]) && (pT < ptEdgeHftRatio[i + 1]))
      return i;
  }
  return nPtBinsHftRatio - 1 ;
}

int getEtaIndexHftRatio(double Eta)
{
  for (int i = 0; i < nEtasHftRatio; i++)
  {
    if ((Eta >= EtaEdgeHftRatio[i]) && (Eta < EtaEdgeHftRatio[i + 1]))
      return i;
  }
  return nEtasHftRatio - 1 ;
}

int getVzIndexHftRatio(double Vz)
{
  for (int i = 0; i < nVzsHftRatio; i++)
  {
    if ((Vz >= VzEdgeHftRatio[i]) && (Vz < VzEdgeHftRatio[i + 1]))
      return i;
  }
  return nVzsHftRatio - 1 ;
}

int getPhiIndexHftRatio(double Phi)
{
  for (int i = 0; i < nPhisHftRatio; i++)
  {
    if ((Phi >= PhiEdgeHftRatio[i]) && (Phi < PhiEdgeHftRatio[i + 1]))
      return i;
  }
  return nPhisHftRatio - 1 ;
}

TVector3 smearPosData(int const iParticleIndex, double const vz, int cent, TLorentzVector const& rMom, TVector3 const& pos)
{
  int const iEtaIndex = getEtaIndexDca(rMom.PseudoRapidity());
  int const iVzIndex = getVzIndexDca(vz);
  int const iPtIndex = getPtIndexDca(rMom.Perp());

  double sigmaPosZ = 0;
  double sigmaPosXY = 0;

  h2Dca[iParticleIndex][iEtaIndex][iVzIndex][cent][iPtIndex]->GetRandom2(sigmaPosXY, sigmaPosZ);
  if(sigmaPosXY==0 || sigmaPosZ==0 || sigmaPosXY != sigmaPosXY || sigmaPosZ != sigmaPosZ)
  {
    cout<<iParticleIndex<<" "<<iEtaIndex<<" "<<iVzIndex<<" "<<cent<<" "<<iPtIndex<<endl;
  }

  sigmaPosZ *= 1.e4; //convert to micrometers
  sigmaPosXY *= 1.e4;

  TVector3 newPos(pos);
  newPos.SetZ(0);
  TVector3 momPerp(-rMom.Vect().Y(), rMom.Vect().X(), 0.0); //vector perpendicular to momentum in xy plane

  newPos -= momPerp.Unit() * sigmaPosXY;

  return TVector3(newPos.X(), newPos.Y(), pos.Z() + sigmaPosZ);
}

TVector3 getVertex(int const centrality)
{
  double rdmVz;

  if (h1Vz[centrality]->GetEntries() == 0) rdmVz = 0.;
  else
  {
    do rdmVz = h1Vz[centrality]->GetRandom() * 1e4;
    while (fabs(rdmVz) > gVzCut);
  }

  return TVector3(0., 0., rdmVz);
}

bool tpcReconstructed(int iParticleIndex, float charge, int cent, TLorentzVector const& mom)
{
  TH1D* h = NULL;

  if (iParticleIndex == 0)
  {
    if (charge > 0) h = hTpcPiPlus[cent];
    else h = hTpcPiMinus[cent];
  }
  else
  {
    if (charge > 0) h = hTpcKPlus[cent];
    else h = hTpcKMinus[cent];
  }

  int const bin = h->FindBin(mom.Perp());

  return gRandom->Rndm() < h->GetBinContent(bin);
}

bool matchHft(int const iParticleIndex, double const vz, int const cent, TLorentzVector const& mom)
{
  int const iEtaIndex = getEtaIndexHftRatio(mom.PseudoRapidity());
    int const iVzIndex = getVzIndexHftRatio(vz);
    int const iPhiIndex = getPhiIndexHftRatio(mom.Phi());

    int bin = hHftRatio1[iParticleIndex][iEtaIndex][iVzIndex][iPhiIndex][cent]->FindBin(mom.Perp());

    int const bin_corr = hHftRatioCorrect[iParticleIndex][iEtaIndex][iVzIndex][iPhiIndex][0]->FindBin(mom.Perp());

  //additional correction same as for Run14
  double AdditionalCorrectFactor = hHftRatioCorrect[iParticleIndex][iEtaIndex][iVzIndex][iPhiIndex][0]->GetBinContent(bin_corr);
  if(fabs(AdditionalCorrectFactor - 0.) < 1.e-5) AdditionalCorrectFactor = 1.0;
  if(mom.Perp() > 3.0) AdditionalCorrectFactor = 1.0;
  if(cent == 0 && mom.Perp() > 4.0) 
  {    
    return gRandom->Rndm() < hHftRatio1[iParticleIndex][iEtaIndex][iVzIndex][iPhiIndex][cent+1]->GetBinContent(bin)/AdditionalCorrectFactor;
  }
  else
  { 
    return gRandom->Rndm() < hHftRatio1[iParticleIndex][iEtaIndex][iVzIndex][iPhiIndex][cent]->GetBinContent(bin)/AdditionalCorrectFactor;
  }

}
//______________________________________________________________________________________________________________________

void write()
{
  result->cd();
  nt->Write();
  result->Close();
}

void initEvtGen()
{
  cout<<"initEvtGen start..."<<endl;
  EvtRandomEngine* eng = 0;
  eng = new EvtSimpleRandomEngine();
  cout<<"setting random engine..."<<endl;
  EvtRandom::setRandomEngine((EvtRandomEngine*)eng);
  cout<<"done"<<endl;
  EvtAbsRadCorr* radCorrEngine = 0;
  std::list<EvtDecayBase*> extraModels;

  EvtExternalGenList genList;
  radCorrEngine = genList.getPhotosModel();
  extraModels = genList.getListOfModels();

  TString Decay_DEC="StRoot/StarGenerator/EvtGen1_06_00/DECAY.DEC";
  TString Evt_pdl="StRoot/StarGenerator/EvtGen1_06_00/evt.pdl";
  EvtGen *myGenerator=new EvtGen(Decay_DEC,Evt_pdl,(EvtRandomEngine*)eng,radCorrEngine, &extraModels);
  starEvtGenDecayer=new StarEvtGenDecayer(myGenerator);
  starEvtGenDecayer->SetDecayTable("Dpm.D_DALITZ.DEC");
  starEvtGenDecayer->SetDebug(0);
}
