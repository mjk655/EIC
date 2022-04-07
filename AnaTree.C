#include "AnaTree.h"

TH2F *hB0_EtaP = new TH2F("hB0_EtaP","hB0_EtaP",160,-8,8,700,-2,5);
TH2F *hB0bar_EtaP = new TH2F("hB0bar_EtaP","hB0bar_EtaP",160,-8,8,700,-2,5);
TH2F *hBp_EtaP = new TH2F("hBp_EtaP","hBp_EtaP",160,-8,8,700,-2,5);
TH2F *hBm_EtaP = new TH2F("hBm_EtaP","hBm_EtaP",160,-8,8,700,-2,5);
TH2F *hBs_EtaP = new TH2F("hBs_EtaP","hBs_EtaP",160,-8,8,700,-2,5);
TH2F *hBsbar_EtaP = new TH2F("hBsbar_EtaP","hBsbar_EtaP",160,-8,8,700,-2,5);
TH2F *hLb_EtaP = new TH2F("hLb_EtaP","hLb_EtaP",160,-8,8,700,-2,5);
TH2F *hLbbar_EtaP = new TH2F("hLbbar_EtaP","hLbbar_EtaP",160,-8,8,700,-2,5);
TH2F *hUp_EtaP = new TH2F("hUp_EtaP","hUp_EtaP",160,-8,8,700,-2,5);
TH2F *hJpsi_EtaP = new TH2F("hJpsi_EtaP","hJpsi_EtaP",160,-8,8,700,-2,5);
TH2F *hD0_EtaP = new TH2F("hD0_EtaP","hD0_EtaP",160,-8,8,700,-2,5);
TH2F *hD0bar_EtaP = new TH2F("hD0bar_EtaP","hD0bar_EtaP",160,-8,8,700,-2,5);
TH2F *hUp_Muon_EtaP = new TH2F("hUp_Muon_EtaP","hUp_Muon_EtaP",160,-8,8,700,-2,5);
TH2F *hUp_Electron_EtaP = new TH2F("hUp_Electron_EtaP","hUp_Electron_EtaP",160,-8,8,700,-2,5);
TH2F *hB2D0_EtaP = new TH2F("hB2D0_EtaP","hB2D0_EtaP",160,-8,8,700,-2,5);
TH2F *hB2Dp_EtaP = new TH2F("hB2Dp_EtaP","hB2Dp_EtaP",160,-8,8,700,-2,5);
TH2F *hB2L_EtaP = new TH2F("hB2L_EtaP","hB2L_EtaP",160,-8,8,700,-2,5);
TH2F *hB2Jpsi_EtaP = new TH2F("hB2Jpsi_EtaP","hB2Jpsi_EtaP",160,-8,8,700,-2,5);
TH2F *hB2K_EtaP = new TH2F("hB2K_EtaP","hB2K_EtaP",160,-8,8,700,-2,5);
TH2F *he_EtaP = new TH2F("he_EtaP","he_EtaP",160,-8,8,700,-2,5);
TTree *ch = new TTree("ch","ch");

void AnaTree(int process=1){
   
    char filename[100];
    sprintf(filename,"Hist_%i.root",process);
    cout << "\n ->Now in ROOT macro \n ->Will save data in " << filename << " \n " << endl; 
    int loop=-1;//100000;
    TFile *f = new TFile(filename,"RECREATE");

//    TTree *ch = new TTree("ch","ch");
    ch->Branch("mntracks",&mntracks,"mntracks/I");
    ch->Branch("mpx",mpx,"mpx[mntracks]/F");
    ch->Branch("mpy",mpy,"mpy[mntracks]/F");
    ch->Branch("mpz",mpz,"mpz[mntracks]/F");
    ch->Branch("mvx",mvx,"mvx[mntracks]/F");
    ch->Branch("mvy",mvy,"mvy[mntracks]/F");
    ch->Branch("mvz",mvz,"mvz[mntracks]/F");
    ch->Branch("mm",mm,"mm[mntracks]/F");
    ch->Branch("me",me,"me[mntracks]/F");
    ch->Branch("mid",mid,"mid[mntracks]/I");
    ch->Branch("mpid",mpid,"mpid[mntracks]/I");
    ch->Branch("mfd",mfd,"mfd[mntracks]/I");
    ch->Branch("mld",mld,"mld[mntracks]/I");
    ch->Branch("mx",&mx,"mx/F");
    ch->Branch("mxtp",&mxtp,"mxtp/F");
    ch->Branch("mxbp",&mxbp,"mxbp/F");
    ch->Branch("mtp",&mtp,"mtp/F");
    ch->Branch("mbp",&mbp,"mbp/F");
    ch->Branch("mphi",&mphi,"mphi/F");
    ch->Branch("my",&my,"my/F");
    ch->Branch("mq2",&mq2,"mq2/F");
    ch->Branch("mw2",&mw2,"mw2/F");
    ch->Branch("msp",&msp,"msp/I");
    ch->Branch("mflag",&mflag,"mflag/I");    
    TH2F *hxQ2 = new TH2F("hxQ2","hxQ2",100,-5,0,100,-1,5);
    TH1I *hEvents = new TH1I("hEvents","hEvents",100000,1E6+0.5-1000,1E6+100000.5-1000);
    int nTracks=0;
    int nevents=0;
    int newEvent=0;
    int flagevts=0;
    ifstream data("ep_allQ2.20x100.small.txt");//pythia.ep.20x100.5Mevents.5.RadCor=0.Q2-0.9.txt");
    if(data.is_open()){
	string line;
	std::getline (data,line);
	if (line.find("PYTHIA EVENT FILE") != string::npos){
	    for(int i = 0;i<5;i++)std::getline (data,line);//skip past header                                                                                                                                                                                             
	    newEvent=1;
	}else{
	    cout <<"Fatal error 1 "<<endl;
	    exit;
	}
	int id[200];
	int fDau[200];
	int lDau[200];
	int par[200];
	float px[200];
	float py[200];
	float pz[200];
	float E[200];
	int flag=0;
	while(true  && nevents!=loop){
	    if(newEvent==1){
		if(nevents%10000 == 0)cout << "On event " << nevents << endl;
		std::getline (data,line);
		if(data.eof())break;
		nTracks = getEvent(line,hxQ2);
		newEvent=2;// Make sure next line is == divider
		flag=0;
	    }
	    if(newEvent==2){
		std::getline (data,line);
		if (line.find("============================================") != string::npos)newEvent=3;
		else{
		    cout <<"Fatal error 2 "<<endl;
		    cout << line << endl;
		    exit;
		}
		newEvent=3;//Go to track level lines
	    }
	    if(newEvent==3){//Fill track information
		std::getline (data,line);
		if (line.find("=============== Event finished ===============") != string::npos){
		    if(flag){
			mflag=flag;
			doAna(id,px,py,pz,E,fDau,lDau,par,nTracks);
		    }else{
			mflag=0;
			doAna(id,px,py,pz,E,fDau,lDau,par,nTracks);//COmment out here if not doing bkg sample
		    }

		    flagevts+=flag;
		    newEvent=1;
		    nevents++;
		}
		else{
		    int prox = getTrack(line,id,px,py,pz,E,fDau,lDau,par); 
		    if(prox>0)flag=1;
		}
	    }
	}
    }
    else{
	cout <<"Fatal Error: File not opened" << endl;
    }
    if(nevents>0)cout <<"\n\nNumber of flaged events " << flagevts << " out of total events " << nevents << " ratio " << float(flagevts)/float(nevents) << endl;

    hEvents->Fill(nevents);
    
    ch->Write();
    hxQ2->Write();
    hEvents->Write();
    hB0_EtaP->Write();
    hB0bar_EtaP->Write();
    hBp_EtaP->Write();
    hBm_EtaP->Write();
    hBs_EtaP->Write();
    hBsbar_EtaP->Write();
    hLb_EtaP->Write();
    hLbbar_EtaP->Write(); 
    hUp_EtaP->Write();
    hJpsi_EtaP->Write();
    hUp_Muon_EtaP->Write();
    hUp_Electron_EtaP->Write();
    hD0_EtaP->Write();
    hD0bar_EtaP->Write();
    hB2D0_EtaP->Write();
    hB2Dp_EtaP->Write();
    hB2L_EtaP->Write();
    hB2Jpsi_EtaP->Write();
    hB2K_EtaP->Write();    
}
int getEvent(string line,TH2F *h){

    vector <float> x; 
    istringstream iss(line);
    for(double val; iss >> val; )
	x.push_back(val);

    h->Fill(log10(x[12]),log10(x[11]));
    mphi = x[15];

    mtp = x[5];
    mxtp = x[6];
    mbp = x[7];
    mxbp = x[8];
    my = x[10];
    mq2 = x[11];
    mx = x[12];
    mw2 = x[13];
    mntracks = int(x[29]);
    msp = x[3];
    return x[29];
}
int getTrack(string line,int id[], float px[], float py[], float pz[], float E[], int fDau[], int lDau[],int par[]){
    
    vector <float> x;
    istringstream iss(line);
    for(double val; iss >> val; )
        x.push_back(val);

    id[int(x[0]-1)]   = x[2];//PID
    px[int(x[0]-1)]   = x[6];//Px
    py[int(x[0]-1)]   = x[7];//Py
    pz[int(x[0]-1)]   = x[8];//Pz
    E[int(x[0]-1)]    = x[9];//E
    fDau[int(x[0]-1)] = x[4];// fist daughter index
    lDau[int(x[0]-1)] = x[5];// last daughter index
    par[int(x[0]-1)]  = x[3];// parent index;
    mpx[int(x[0]-1)] = x[6];
    mpy[int(x[0]-1)] = x[7];
    mpz[int(x[0]-1)] = x[8];
    mm[int(x[0]-1)]  = x[10];
    mvx[int(x[0]-1)] = x[11];
    mvy[int(x[0]-1)] = x[12];
    mvz[int(x[0]-1)] = x[13];
    me[int(x[0]-1)]  = x[9];
    mfd[int(x[0]-1)] = int(x[4]);
    mld[int(x[0]-1)] = int(x[5]);
    mpid[int(x[0]-1)] = int(x[3]);
    mid[int(x[0]-1)] = int(x[2]);
    if(fabs(x[2]) == 443)cout <<" JPsi ISUB " << msp << endl;
    if(fabs(x[2]) == 553)cout <<" Upsilon ISUB " << msp << endl;
    if(fabs(x[2]) == 421 ||     //D0
       fabs(x[2]) == 411 ||     //D+    
       fabs(x[2]) == 431 ||     //Ds
       fabs(x[2]) == 4122 ||    //Lc          
       fabs(x[2]) == 443 ||     //Jpsi
       fabs(x[2]) == 511 ||     //B0
       fabs(x[2]) == 521 ||     //B+
       fabs(x[2]) == 531 ||     //Bs 
       fabs(x[2]) == 541 ||     //Bc
       fabs(x[2]) == 5122 ||    //Lambda b
       fabs(x[2]) == 5112 ||    //Sigma b-
       fabs(x[2]) == 5212 ||    //Sigma b0
       fabs(x[2]) == 5222 ||    //Sigam b+
       fabs(x[2]) == 5132 ||    //Cascade b-
       fabs(x[2]) == 5232 ||    //Cascade b0
       fabs(x[2]) == 5332 ||    //Omega b-
       fabs(x[2]) == 553 ||     //Upsilon 1S
       fabs(x[2]) == 100553 ||  //Upsilon 2S
       fabs(x[2]) == 200553)    //Upsilon 3S
	return 1;
    else return 0;
}
void doAna(int id[], float px[], float py[], float pz[], float E[], int fDau[], int lDau[], int par[], int nTracks){
    ch->Fill();
    for(int i =0; i < nTracks; i++){
	//cout << i << " " << nTracks << " "  << id[i] << " " << px[i] << " " << py[i] << " " << pz[i] << " " << E[i] << " " << fDau[i] << " " << lDau[i] << endl;
	double  p = sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]);
	double eta = TMath::ATanH(pz[i]/p);
	if(id[i] == 421) hD0_EtaP->Fill(eta,log10(p));
	if(id[i] == -421) hD0bar_EtaP->Fill(eta,log10(p));
	if(id[i] == 511) hB0_EtaP->Fill(eta,log10(p));
	if(id[i] == -511) hB0bar_EtaP->Fill(eta,log10(p));
	if(id[i] == 521) hBp_EtaP->Fill(eta,log10(p));
	if(id[i] == -521) hBm_EtaP->Fill(eta,log10(p));
	if(id[i] == 531) hBs_EtaP->Fill(eta,log10(p));
	if(id[i] == -531) hBsbar_EtaP->Fill(eta,log10(p));
	if(id[i] == 5122) hLb_EtaP->Fill(eta,log10(p));
	if(id[i] == -5122) hLbbar_EtaP->Fill(eta,log10(p));
        if(id[i] == 443)hJpsi_EtaP->Fill(eta,log10(p));
	if(id[i] == 553 || id[i] == 100553 || id[i] == 200553) hUp_EtaP->Fill(eta,log10(p));
	if( fabs(id[i]) == 13 && (id[par[i]-1] == 553 || id[par[i]-1] == 100553 || id[par[i]-1] == 200553))hUp_Muon_EtaP->Fill(eta,log10(p));
	if( fabs(id[i]) == 11 && (id[par[i]-1] == 553 || id[par[i]-1] == 100553 || id[par[i]-1] == 200553))hUp_Electron_EtaP->Fill(eta,log10(p));
	if( fabs(id[i]) == 421 && (fabs(id[par[i]-1]) == 521 || fabs(id[par[i]-1]) == 511))hB2D0_EtaP->Fill(eta,log10(p));
	if( fabs(id[i]) == 411 && (fabs(id[par[i]-1]) == 521 || fabs(id[par[i]-1]) == 511))hB2Dp_EtaP->Fill(eta,log10(p));
	if( (fabs(id[i]) == 11 || fabs(id[i]) == 13) && (fabs(id[par[i]-1]) == 521 || fabs(id[par[i]-1]) == 511))hB2L_EtaP->Fill(eta,log10(p));
	if( id[i] == 443 && (fabs(id[par[i]-1]) == 521 || fabs(id[par[i]-1]) == 511))hB2Jpsi_EtaP->Fill(eta,log10(p));
	if( (fabs(id[i]) == 321 || fabs(id[i]) == 311) && (id[par[i]-1] == 521 || id[par[i]-1] == 511))hB2K_EtaP->Fill(eta,log10(p));
    }
}
