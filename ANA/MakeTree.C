#include "MakeTree.h"
//=========================== Only change these
double q2_low=10;
double q2_high=10000;
char infile[100] = "541List_fulltree_q210.list";
int limit = 10000; 
char outfile[100] = "D0Tree_541_q210_full.root";
//=========================== 

float Mq2;
float Mx;
float Mxtp;
float Mw;

void MakeTree(){
    TFile *partonf = new TFile("PartonWeights.root","READ");
    TH2F* N_GLUON = (TH2F*)partonf->Get("N_GLUON");N_GLUON->SetDirectory(0);    
    TFile *outFile = new TFile(outfile,"RECREATE");
    TTree *tree = new TTree("tree","tree");
    tree->Branch("mq2",&Mq2,"mq2/F");
    tree->Branch("mx",&Mx,"mx/F");
    tree->Branch("mw",&Mw,"mw/F");
    tree->Branch("mxtp",&Mxtp,"mxtp/F");

    
    std::ifstream listOfFiles(infile);
    int iter = 0;
    if (listOfFiles.is_open())
    {
	std::string file;
        while (getline(listOfFiles, file) && iter < limit)
        {
            cout << "Looping on new file:" << file << endl;
            doLoop(file,tree,N_GLUON);
	    iter++;
        }
    }
    else
    {
        cout << "Could not open list of files. ABORT!" << endl;
        return;
    }

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
    for(int i =0  ;i < iloop; i++){
        ch->GetEvent(i);

        if(!(msp==135 || msp==136))continue; //Select Photon Gluon Fusion
	if(mtp!=21){
	    cout <<" Now printing event record" << endl;
	    for(int tr = 0; tr<mntracks ; tr++){
		cout << tr+1 << ": " << mid[tr] << " " << mtp << endl;
	    }
	    continue;
	}
	if(mq2 < q2_low || mq2> q2_high)continue;

        for(int tr = 0; tr<mntracks ; tr++){
            if(fabs(mid[tr])!=421 && fabs(mid[tr])!=411 && fabs(mid[tr])!=431 && fabs(mid[tr])!=4122)continue;
            int motherID = 0;
            if(mpid[tr]>0)motherID=mid[mpid[tr]-1];
            int gmotherLine = mpid[mpid[tr]-1];
            int gmotherID=0;
            if(gmotherLine>0)gmotherID = mid[gmotherLine-1];
            int ggmotherLine=-1;
            if(gmotherLine>0)ggmotherLine=mpid[mpid[mpid[tr]-1]-1];
            int ggmotherID=0;
            if(ggmotherLine>0)ggmotherID = mid[ggmotherLine-1];
            int gggmotherLine=-1;
            if(ggmotherLine>0)gggmotherLine=mpid[mpid[mpid[mpid[tr]-1]-1]-1];
            int gggmotherID=0;
            if(gggmotherLine>0)gggmotherID = mid[gggmotherLine-1];
            int ggggmotherLine=-1;
            if(gggmotherLine>0)ggggmotherLine=mpid[mpid[mpid[mpid[mpid[tr]-1]-1]-1]-1];
            int ggggmotherID=0;
            if(ggggmotherLine>0)ggggmotherID = mid[ggggmotherLine-1];
            int gggggmotherLine=-1;
            if(ggggmotherLine>0)gggggmotherLine=mpid[mpid[mpid[mpid[mpid[mpid[tr]-1]-1]-1]-1]-1];
            int gggggmotherID=0;
            if(gggggmotherLine>0)gggggmotherID = mid[gggggmotherLine-1];
	    int flag = 0;
            flag+=isFromB(motherID);
            flag+=isFromB(gmotherID);
            flag+=isFromB(ggmotherID);
            flag+=isFromB(gggmotherID);
            flag+=isFromB(ggggmotherID);
            flag+=isFromB(gggggmotherID);
            if(flag!=0)continue;
            Mq2 = mq2;
            Mx = mx;
            Mxtp = mxtp;
	    Mw = N_GLUON->GetBinContent(N_GLUON->FindBin(log10(mxtp),log10(mq2)));
	    tree->Fill();
        }	
    }
    delete ch;
}
