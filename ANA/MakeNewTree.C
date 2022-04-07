#include "MakeTree.h"
//=========================== Only change these
char infile[100] = "../../Charm/out_541/out_541.root";
char outfile[100] = "D0Tree_541_new.root";
//=========================== 

float Mq2;
float Mx;
float Mxtp;
float Mw;

void MakeNewTree(){
    TFile *partonf = new TFile("PartonWeights.root","READ");
    TH2F* N_GLUON = (TH2F*)partonf->Get("N_GLUON");N_GLUON->SetDirectory(0);    
    TFile *outFile = new TFile(outfile,"RECREATE");
    TTree *tree = new TTree("tree","tree");
    tree->Branch("mq2",&Mq2,"mq2/F");
    tree->Branch("mx",&Mx,"mx/F");
    tree->Branch("mw",&Mw,"mw/F");
    tree->Branch("mxtp",&Mxtp,"mxtp/F");

    doLoop(infile,tree,N_GLUON);

    tree->Write();
    outFile->Close();
}
int isFromB(int id){
    if( fabs(id) == 5 ||     //b
        fabs(id) == 511 ||     //B0
	fabs(id) == 521 ||     //B+
	fabs(id) == 531 ||     //Bs
	fabs(id) == 541 ||     //Bc
	fabs(id) == 5122 ||    //Lambda b
	fabs(id) == 5112 ||    //Sigma b-
	fabs(id) == 5212 ||    //Sigma b0
	fabs(id) == 5222 ||    //Sigam b+
	fabs(id) == 5132 ||    //Cascade b-
	fabs(id) == 5232 ||    //Cascade b0
	fabs(id) == 5332 ||    //Omega b-
	fabs(id) == 553 ||     //Upsilon 1S
	fabs(id) == 100553 ||  //Upsilon 2S
	fabs(id) == 200553)    //Upsilon 3S
        return 1;
    else return 0;
}
void doLoop(std::string loopfile,TTree *tree, TH2F* N_GLUON){
    float mpx[1000];float mpy[1000]; float mpz[1000]; float me[1000];
    float mvx[1000];float mvy[1000]; float mvz[1000];
    int mid[1000]; int mpid[1000];int mfd[1000]; int mld[1000];
    float my; float mx; float mq2; float mw2; int mntracks; int msp;
    float mphi; float mxtp;float mxbp; float mtp;float mbp;
    TChain *ch = new TChain("ch","ch");
    char _file[100];
    sprintf(_file,"%s",loopfile.c_str());
    ch->AddFile(_file);
    ch->SetBranchAddress("mntracks",&mntracks);
    ch->SetBranchAddress("mpx",mpx);
    ch->SetBranchAddress("mpy",mpy);
    ch->SetBranchAddress("mpz",mpz);
    ch->SetBranchAddress("mvx",mvx);
    ch->SetBranchAddress("mvy",mvy);
    ch->SetBranchAddress("mvz",mvz);
    ch->SetBranchAddress("me",me);
    ch->SetBranchAddress("mid",mid);
    ch->SetBranchAddress("mpid",mpid);
    ch->SetBranchAddress("mfd",mfd);
    ch->SetBranchAddress("mld",mld);
    ch->SetBranchAddress("mx",&mx);
    ch->SetBranchAddress("mq2",&mq2);
    ch->SetBranchAddress("mw2",&mw2);
    ch->SetBranchAddress("msp",&msp);
    ch->SetBranchAddress("my",&my);
    ch->SetBranchAddress("mphi",&mphi);
    ch->SetBranchAddress("mxtp",&mxtp);
    ch->SetBranchAddress("mxbp",&mxbp);
    ch->SetBranchAddress("mtp",&mtp);
    ch->SetBranchAddress("mbp",&mbp);
    int iloop = ch->GetEntries();
    cout << iloop << " entries" << endl;
    for(int i =0  ;i < iloop; i++){
        ch->GetEvent(i);

        //if(!(msp==135 || msp==136))continue; //Select Photon Gluon Fusion
	if(mtp!=21){
	    cout <<" Now printing event record" << endl;
	    for(int tr = 0; tr<mntracks ; tr++){
		cout << tr+1 << ": " << mid[tr] << " " << mtp << endl;
	    }
	    continue;
	}
	
        for(int tr = 0; tr<mntracks ; tr++){
            if(fabs(mid[tr])!=421 && fabs(mid[tr])!=411 && fabs(mid[tr])!=431 && fabs(mid[tr])!=4122)continue;
	    Mq2 = mq2;
            Mx = mx;
            Mxtp = mxtp;
	    Mw = N_GLUON->GetBinContent(N_GLUON->FindBin(log10(mxtp),log10(mq2)));
	    tree->Fill();
        }	
    }
    delete ch;
}
