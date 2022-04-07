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
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TLegend.h>

using namespace std;

TH1F *hCS_x1[100];
TH1F *hCS_x2[100];
TH1F *hF2_x[100];

TH1F *hCS_x1_eAu[100];
TH1F *hCS_x2_eAu[100];
TH1F *hF2_x_eAu[100];

TH1F *hCS_x1_eAu_pseudo[100];
TH1F *hCS_x2_eAu_pseudo[100];
TH1F *hF2_x_eAu_pseudo[100];

TH1F *hCS_x1_eAu_er[100];
TH1F *hCS_x2_eAu_er[100];
TH1F *hF2_x_eAu_er[100];
TH1F *hCS_x1_eAu_ernew[100];
TH1F *hCS_x2_eAu_ernew[100];
TH1F *hF2_x_eAu_ernew[100];

TH2F *hCS_x1_eAu_toy[100];
TH2F *hCS_x2_eAu_toy[100];
TH2F *hF2_x_eAu_toy[100];
TH2F *hCS_x1_eAu_toynew[100];
TH2F *hCS_x2_eAu_toynew[100];
TH2F *hF2_x_eAu_toynew[100];


void setXbins(int nbins_x, double binning_x[]){
    for(int i = 0;i<nbins_x+1;i++)
        binning_x[i] = -4. + i * 4. / nbins_x;
}
void setQ2bins(int nbins_q2, double binning_q2[]){
    for(int i = 0;i<nbins_q2+1;i++)
        binning_q2[i] = i * 2.6 / nbins_q2;
}
void norm(TH1F *h, char label[]);
void getRat(TH1F* h, TH1F* h1,int er);
void fitF2(int nbins_q2, int nbins_x,double binx[], double biny[],TH1F** hF2_x,TH1F** hCS_x1,TH1F** hCS_x2);
double getF2(double v1, double v2, double e1, double e2, double qq, double xx);
double getF2er(double v1, double v2, double e1, double e2, double qq, double xx);
void fill(char file1[100],TH1F** hCS_x1,TH1F** hCS_x1_eAu,TH1F* hq2bin,TH2F* Weights,TH2F *Weights_eAu,TH2F* N_GLUON,TH2F* vars,int vary, double ww);
void fill(char file1[100],TH1F** hCS_x1,TH1F** hCS_x1_eAu,TH1F* hq2bin,TH2F* Weights,TH2F *Weights_eAu,TH2F* N_GLUON,int vary, double ww);
double getError(TH1F* h, double c, int a);
double getError1(TH1F* h, double c, int a);
void getReplica(TH2F* n, TH2F** set,int nset,double rik[][40],int vary_iter);
void getNewReplica(TH2F* n, TH2F** set,int nset, double weights[],double rik[][40], int nvary);
void getPDFError(TH2F* n, TH2F** _set,int vary_iter);
double getNewPDFError(TH2F* n, TH2F** _set,int nset, double weights[],double rik[][40], int nvary);
