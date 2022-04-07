#include <sstream>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <TChain.h>
#include <TH1F.h>
#include <TF1.h>
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
using namespace std;

int getEvent(string line,TH2F *h);
int getTrack(string line,int id[], float px[], float py[], float pz[], float E[], int fDau[], int lDau[],int par[]);
void doAna(int id[], float px[], float py[], float pz[], float E[], int fDau[], int lDau[], int par[], int nTracks);
TH2F *hxQ2;
TH1I *hEvents;
//TTree *ch;
float mpx[1000]; float mpy[1000]; float mpz[1000]; float me[1000];
float mm[1000]; float mvx[1000]; float mvy[1000]; float mvz[1000];
int mid[1000]; int mpid[1000];int mfd[1000]; int mld[1000];
float my; float mx; float mq2; float mw2; int mntracks; int msp;int mflag;
float mphi; float mxtp;float mxbp; float mtp;float mbp; 
