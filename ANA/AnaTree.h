#include <sstream>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <TRandom.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TMath.h>
#include <TFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVector3.h>
#include <TSystem.h>
#include <THelixTrack.h>
#include <StHelixD.h>
#include <StPhysicalHelixD.h>
#include <SystemOfUnits.h>
//#include "LHAPDF/LHAPDF.h"
//#include "LHAPDF/Version.h"
//#include "/gpfs01/star/pwg/mkelsey/EIC/PYTHIA/ANA/EICVertexFitter/EICVertexFitter.h"
using namespace std;

TLorentzVector * smearMomATHENA(TLorentzVector const * b, TGraph *g);
TLorentzVector * smearMomEIC(TLorentzVector const * b);
TLorentzVector * smearMomEIC2(TLorentzVector const * b, TH2F* h, TH2F* h2, int id);
TLorentzVector * smearMomEIC_SlowParticle(TLorentzVector const * b);
bool isFlag(int id, int arr[]);
TVector3 smearPosATHENA(TLorentzVector const* b, TVector3 const& pos, TGraph *g1, TGraph *g2);
TVector3 smearPosEIC(TLorentzVector const* b, TVector3 const& pos);
TVector3 smearPosEICLBNL(TLorentzVector const* b, TVector3 const& pos);
TVector3 smearPosEICLANL(TLorentzVector const* b, TVector3 const& pos);
TVector3 smearPosEIC2(TLorentzVector const* b, TVector3 const& pos,TH1F* DCARes_T_Eta1,TH1F* DCARes_T_Eta2,TH1F* DCARes_T_Eta3,TH1F* DCARes_Z_Eta1,TH1F* DCARes_Z_Eta2,TH1F* DCARes_Z_Eta3);
TVector3 smearPosEIC3(TLorentzVector const * b, TVector3 const& pos, TH1F* hp,TH1F *hm,TH1F* DCARes_Z_Eta1,TH1F* DCARes_Z_Eta2,TH1F* DCARes_Z_Eta3, int id);
float dcaSigned(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
float dcaXY(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex, int m);
float dcaZ(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex);
//bool isGoodD0(TLorentzVector const& P,TVector3 const& p1, TVector3 const& pos1,TVector3 const& p2, TVector3 const& pos2,TVector3 const& vertex);
bool isGoodD0(TLorentzVector* P,TVector3 p1, TVector3 pos1,TVector3 p2, TVector3 pos2,TVector3 vertex, int d1,int d2);
bool isGoodD0_XY(TLorentzVector* P,TVector3 p1, TVector3 pos1,TVector3 p2, TVector3 pos2,TVector3 vertex, int d1, int d2,TH2F* h2,int i, double e1, double e2);
bool isSlowPion(TVector3 p1);
int isFromB(int id);
int isFromD(int id);
bool isKPiDecay(int id1, int id2);
bool isWSKPiDecay(int id1, int id2);
bool iseeDecay(int id1, int id2);
bool isK2PiDecay(int id1, int id2, int id3);
bool passTracking(TH2F*h,double p, double eta);
bool passTrackingATHENA(double p, double eta);
int getEtaIdx(double eta);
int getPtIdx(double pt);
bool passPID(int p1, int p2, double pt1, double pt2, double eta1, double eta2);
