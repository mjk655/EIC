#include "AnaTree.h"

TH3F * D0_DecayL = new TH3F("D0_DecayL","D0_DecayL",40,0,20,200,0,2,40,-4,4);
TH3F * D0_DCA = new TH3F("D0_DCA","D0_DCA",40,0,20,100,0,0.1,40,-4,4);
TH3F * D0_DDCA = new TH3F("D0_DDCA","D0_DDCA",40,0,20,500,0,0.5,40,-4,4);
TH3F * D0_D1DCA = new TH3F("D0_D1DCA","D0_D1DCA",40,0,20,500,0,0.5,40,-4,4);
TH3F * D0_D2DCA = new TH3F("D0_D2DCA","D0_D2DCA",40,0,20,500,0,0.5,40,-4,4);
TH3F * D0_CosTheta = new TH3F("D0_CosTheta","D0_CosTheta",40,0,20,100,0.5,1,40,-4,4);
TH3F * D0_Boost = new TH3F("D0_Boost","D0_Boost",40,0,20,200,0.,1,40,-4,4);
TH3F * D0_DecayL_WS = new TH3F("D0_DecayL_WS","D0_DecayL_WS",40,0,20,200,0,2,40,-4,4);
TH3F * D0_DCA_WS = new TH3F("D0_DCA_WS","D0_DCA_WS",40,0,20,100,0,0.1,40,-4,4);
TH3F * D0_DDCA_WS = new TH3F("D0_DDCA_WS","D0_DDCA_WS",40,0,20,500,0,0.5,40,-4,4);
TH3F * D0_D1DCA_WS = new TH3F("D0_D1DCA_WS","D0_D1DCA_WS",40,0,20,500,0,0.5,40,-4,4);
TH3F * D0_D2DCA_WS = new TH3F("D0_D2DCA_WS","D0_D2DCA_WS",40,0,20,500,0,0.5,40,-4,4);
TH3F * D0_CosTheta_WS = new TH3F("D0_CosTheta_WS","D0_CosTheta_WS",40,0,20,100,0.5,1,40,-4,4);
TH3F * D0_Boost_WS = new TH3F("D0_Boost_WS","D0_Boost_WS",40,0,20,200,0.,1,40,-4,4);

void DoAnaATHENA(int index=1){
    bool VTX_FIT = false;
    bool onlyD0 = true;
    cout <<"Now doing DoAnaATHENA.C " << endl;
    char infile[100] = "in.list";// set to in.list for job submission
    cout <<"Reading in Nuclear Modification Weights " << endl;

    TFile *partonfp = new TFile("PartonWeights_PreIC_NNPDF.root","READ");
    TH2F* N_GLUONp = (TH2F*)partonfp->Get("N_GLUON");N_GLUONp->SetName("N_GLUONp");
    TH2F* N_UQp = (TH2F*)partonfp->Get("N_UQ");N_UQp->SetName("N_UQp");
    TH2F* N_DQp = (TH2F*)partonfp->Get("N_DQ");N_DQp->SetName("N_DQp");
    TH2F* N_SQp = (TH2F*)partonfp->Get("N_SQ");N_SQp->SetName("N_SQp");
    TH2F* N_CQp = (TH2F*)partonfp->Get("N_CQ");N_CQp->SetName("N_CQp");
    TH2F* N_UBARQp = (TH2F*)partonfp->Get("N_UBARQ");N_UBARQp->SetName("N_UBARQp");
    TH2F* N_DBARQp = (TH2F*)partonfp->Get("N_DBARQ");N_DBARQp->SetName("N_DBARQp");
    TH2F* N_SBARQp = (TH2F*)partonfp->Get("N_SBARQ");N_SBARQp->SetName("N_SBARQp");
    TH2F* N_CBARQp = (TH2F*)partonfp->Get("N_CBARQ");N_CBARQp->SetName("N_CBARQp");

    TFile *partonf = new TFile("PartonWeights.root","READ");
    TH2F* N_GLUON = (TH2F*)partonf->Get("N_GLUON");
    TH2F* N_UQ = (TH2F*)partonf->Get("N_UQ");
    TH2F* N_DQ = (TH2F*)partonf->Get("N_DQ");
    TH2F* N_SQ = (TH2F*)partonf->Get("N_SQ");
    TH2F* N_CQ = (TH2F*)partonf->Get("N_CQ");
    TH2F* N_UBARQ = (TH2F*)partonf->Get("N_UBARQ");
    TH2F* N_DBARQ = (TH2F*)partonf->Get("N_DBARQ");
    TH2F* N_SBARQ = (TH2F*)partonf->Get("N_SBARQ");
    TH2F* N_CBARQ = (TH2F*)partonf->Get("N_CBARQ");

    cout <<"Reading in Vertex Resolution " << endl;
    TFile *vtxres = new TFile("VertexRes_ATHENA.root","READ");
    TH1F* VertexRes_X = (TH1F*) vtxres->Get("VertexRes_X");VertexRes_X->SetName("VertexRes_X");
    TH1F* VertexRes_Y = (TH1F*) vtxres->Get("VertexRes_Y");VertexRes_Y->SetName("VertexRes_Y");
    TH1F* VertexRes_Z = (TH1F*) vtxres->Get("VertexRes_Z");VertexRes_Z->SetName("VertexRes_Z");
    VertexRes_X->SetDirectory(0);
    VertexRes_Y->SetDirectory(0);
    VertexRes_Z->SetDirectory(0);
    cout <<"Reading in Eff and p resolution " << endl;
    TFile *trackeff = new TFile("TrackEff_EtaP_ver2_p8.root","READ");
    TH2F* TrackEff_EtaP = (TH2F*) trackeff->Get("TrackEff_EtaP");
    TrackEff_EtaP->SetDirectory(0);


    TFile *track_res = new TFile("ATHENA_Resolutions_r.root","READ");
    TH1F * Res_Handler = (TH1F*)track_res->Get("Res_Handler");
    int const npt = Res_Handler->GetNbinsX();
    TGraph *gmom_res[npt];
    TGraph *gdca_rphi_res[npt];
    TGraph *gdca_z_res[npt];
    for(int i = 0; i<npt; i++){
        gmom_res[i] = (TGraph*)track_res->Get(Form("gmom_res_%i",i));
        gdca_rphi_res[i] = (TGraph*)track_res->Get(Form("gdca_rphi_res_%i",i));
        gdca_z_res[i] = (TGraph*)track_res->Get(Form("gdca_z_res_%i",i));
    }
    


    cout << "Done with  DCA histograms " << endl;
    TChain *ch = new TChain("ch","ch");
    TH2F *_hxQ2 = new TH2F("_hxQ2","_hxQ2",100,-5,0,100,-1,5);
    TH1I *_hEvents = new TH1I("_hEvents","_hEvents",100000,1E6+0.5-1000,1E6+100000.5-1000);
    std::ifstream listOfFiles(infile);
    if (listOfFiles.is_open())
    {
	std::string file;
        while (getline(listOfFiles, file))
        {
	    char _file[100];
	    sprintf(_file,"%s",file.c_str());
	    //TFile *f = new TFile(_file,"READ");
	    //TH2F* hxQ2 = (TH2F*) f->Get("hxQ2");
	    //TH1F* hEvents = (TH1F*) f->Get("hEvents");
	    cout << "DoAna - Adding :" << file << endl;
	    ch->AddFile(_file);
		//_hxQ2->Add(hxQ2);
		//_hEvents->Add(hEvents);
	    //delete hxQ2;delete hEvents;delete f;
	}
    }
    else
    {
        cout << "DoAna - Could not open list of files. ABORT!" << endl;
        return;
    }

    float mpx[1000];float mpy[1000]; float mpz[1000]; float me[1000];
    float mvx[1000];float mvy[1000]; float mvz[1000];
    int mid[1000]; int mpid[1000];int mfd[1000]; int mld[1000];
    float my; float mx; float mq2; float mw2; int mntracks; int msp;
    float mphi; float mxtp;float mxbp; float mtp;float mbp;

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

// ======
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    TH3::SetDefaultSumw2();

//=======================// ANA related items //=======================// 
    char output[100];
    sprintf(output,"out_%i.root",index);
    cout<<"Saving output to "<< output << endl;
    TFile *savefile = new TFile(output,"RECREATE");
    savefile->cd();
    int particle_flag[3] = {421,411,443}; // Truth flag for HF particles of interest  
    TH3F *D0_Lepton_Corr_PGF = new TH3F("D0_Lepton_Corr_PGF","D0_Lepton_Corr_PGF",100,-5,0,100,0,5,100,0,1);
    TH3F *D0_Lepton_Corr_L0 = new TH3F("D0_Lepton_Corr_L0","D0_Lepton_Corr_L0",100,-5,0,100,0,5,100,0,1);
    TH3F *D0_Lepton_Corr_Other = new TH3F("D0_Lepton_Corr_Other","D0_Lepton_Corr_Other",100,-5,0,100,0,5,100,0,1);
    TH3F *Incl_RedCS = new TH3F("Incl_RedCS","Incl_RedCS",100,-5,0,100,0,5,100,1.8,1.9);
    TH3F *Incl_nRedCS = new TH3F("Incl_nRedCS","Incl_nRedCS",100,-5,0,100,0,5,100,1.8,1.9);

    TH3F *D0_RedCS_BKG = new TH3F("D0_RedCS_BKG","D0_RedCS_BKG",100,-5,0,100,0,5,100,1.8,1.9);
    TH3F *D0_RedCS_Topo_BKG = new TH3F("D0_RedCS_Topo_BKG","D0_RedCS_Topo_BKG",100,-5,0,100,0,5,100,1.8,1.9);
    
    TH3F *D0_RedCS_T = new TH3F("D0_RedCS_T","D0_RedCS_T",100,-5,0,100,0,5,100,1.8,1.9);
    TH3F *D0_RedCS_T_Acc = new TH3F("D0_RedCS_T_Acc","D0_RedCS_T_Acc",100,-5,0,100,0,5,100,1.8,1.9);
    TH3F *B2D0_RedCS_T_Acc = new TH3F("B2D0_RedCS_T_Acc","B2D0_RedCS_T_Acc",100,-5,0,100,0,5,100,1.8,1.9);

    TH3F *D0_RedCS_Topo_T = new TH3F("D0_RedCS_Topo_T","D0_RedCS_Topo_T",100,-5,0,100,0,5,100,1.8,1.9);

    TH3F *D0_nRedCS_Topo_T = new TH3F("D0_nRedCS_Topo_T","D0_nRedCS_Topo_T",100,-5,0,100,0,5,100,1.8,1.9);
    TH3F *D0_nRedCS_Topo_BKG = new TH3F("D0_nRedCS_Topo_BKG","D0_nRedCS_Topo_BKG",100,-5,0,100,0,5,100,1.8,1.9);
    TH3F *D0_nRedCS_T_Acc = new TH3F("D0_nRedCS_T_Acc","D0_nRedCS_T_Acc",100,-5,0,100,0,5,100,1.8,1.9);

    TH3F *D0L_nRedCS_Topo_T = new TH3F("D0L_nRedCS_Topo_T","D0L_nRedCS_Topo_T",100,-5,0,100,0,5,100,0,1);
    TH3F *D0L_nRedCS_Topo_T_L0 = new TH3F("D0L_nRedCS_Topo_T_L0","D0L_nRedCS_Topo_T_L0",100,-5,0,100,0,5,100,0,1);
    TH3F *D0L_nRedCS_Topo_T_PGF = new TH3F("D0L_nRedCS_Topo_T_PGF","D0L_nRedCS_Topo_T_PGF",100,-5,0,100,0,5,100,0,1);
    TH3F *D0L_nRedCS_Topo_T_Other = new TH3F("D0L_nRedCS_Topo_T_Other","D0L_nRedCS_Topo_T_Other",100,-5,0,100,0,5,100,0,1);
    TH3F *D0L_nRedCS_Topo_BKG = new TH3F("D0L_nRedCS_Topo_BKG","D0L_nRedCS_Topo_BKG",100,-5,0,100,0,5,100,0,1);
    TH3F *D0L_nRedCS_T_Acc = new TH3F("D0L_nRedCS_T_Acc","D0L_nRedCS_T_Acc",100,-5,0,100,0,5,100,0,1);
    TH3F *D0L_RedCS_Topo_T = new TH3F("D0L_RedCS_Topo_T","D0L_RedCS_Topo_T",100,-5,0,100,0,5,100,0,1);
    TH3F *D0L_RedCS_Topo_T_L0 = new TH3F("D0L_RedCS_Topo_T_L0","D0L_RedCS_Topo_T_L0",100,-5,0,100,0,5,100,0,1);
    TH3F *D0L_RedCS_Topo_T_PGF = new TH3F("D0L_RedCS_Topo_T_PGF","D0L_RedCS_Topo_T_PGF",100,-5,0,100,0,5,100,0,1);
    TH3F *D0L_RedCS_Topo_T_Other = new TH3F("D0L_RedCS_Topo_T_Other","D0L_RedCS_Topo_T_Other",100,-5,0,100,0,5,100,0,1);
    TH3F *D0L_RedCS_Topo_BKG = new TH3F("D0L_RedCS_Topo_BKG","D0L_RedCS_Topo_BKG",100,-5,0,100,0,5,100,0,1);
    TH3F *D0L_RedCS_T_Acc = new TH3F("D0L_RedCS_T_Acc","D0L_RedCS_T_Acc",100,-5,0,100,0,5,100,0,1);



    TH3F *D0_DCA_RedCS_Topo_BKG = new TH3F("D0_DCA_RedCS_Topo_BKG","D0_DCA_RedCS_Topo_BKG",100,-5,0,100,0,5,400,-2,2);
    TH3F *D0_DCA_RedCS_Topo_Charm = new TH3F("D0_DCA_RedCS_Topo_Charm","D0_DCA_RedCS_Topo_Charm",100,-5,0,100,0,5,400,-2,2);
    TH3F *D0_DCA_RedCS_Topo_Bottom = new TH3F("D0_DCA_RedCS_Topo_Bottom","D0_DCA_RedCS_Topo_Bottom",100,-5,0,100,0,5,400,-2,2);

    TH2D *hD02pi_DCATheta=new TH2D("hD02Pi_DCATheta","D^{0}#rightarrow#pi |DCA_{xy}| vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,0.1);//100.,0.,50.);          
    

    TH2D *hD0_PTheta=new TH2D("hD0_PTheta","D^{0} Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//100.,0.,50.);
    TH2D *hD0bar_PTheta=new TH2D("hD0bar_PTheta","#bar{D}^{0} Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);
    TH2D *hDp_PTheta=new TH2D("hDp_PTheta","D^{+} Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);
    TH2D *hDpbar_PTheta=new TH2D("hDpbar_PTheta","D^{-} Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);
    TH2D *hDs_PTheta=new TH2D("hDs_PTheta","D^{+}_{s} Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);
    TH2D *hDsbar_PTheta=new TH2D("hDsbar_PTheta","D^{-}_{s} Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);
    TH2D *hLc_PTheta=new TH2D("hLc_PTheta","#Lambda_{c}^{+} Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);
    TH2D *hLcbar_PTheta=new TH2D("hLcbar_PTheta","#Lambda_{c}^{-} Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);
    TH2D *hJpsi_PTheta=new TH2D("hJpsi_PTheta","J/#psi Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);
    TH2D *hJpsi2e_PTheta=new TH2D("hJpsi2e_PTheta","J/#psi#rightarrow#it{e} Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);    
    TH2D *hJpsi2mu_PTheta=new TH2D("hJpsi2mu_PTheta","J/#psi#rightarrow#mu Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);              
    TH2D *hB0_PTheta=new TH2D("hB0_PTheta","B^{0} Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);
    TH2D *hB0bar_PTheta=new TH2D("hB0bar_PTheta","#bar{B}^{0} Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);
    TH2D *hBp_PTheta=new TH2D("hBp_PTheta","B^{+} Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);
    TH2D *hBpbar_PTheta=new TH2D("hBpbar_PTheta","B^{-} Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);
    TH2D *hBs_PTheta=new TH2D("hBs_PTheta","B^{0}_{s} Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);
    TH2D *hBsbar_PTheta=new TH2D("hBsbar_PTheta","#bar{B}^{0}_{s} Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);
    TH2D *hLb_PTheta=new TH2D("hLb_PTheta","#Lambda_{b}^{0} Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);
    TH2D *hLbbar_PTheta=new TH2D("hLbbar_PTheta","#bar{#Lambda}_{b}^{0} Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);
    TH2D *hB2D0_PTheta=new TH2D("hB2D0_PTheta","B#rightarrow D^{0} Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);
    TH2D *hB2e_PTheta=new TH2D("hB2e_PTheta","B#rightarrow #it{e} Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);

    TH2D *hD02K_PTheta=new TH2D("hD02K_PTheta","D^{0}#rightarrowK Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);
    TH2D *hD02pi_PTheta=new TH2D("hD02pi_PTheta","D^{0}#rightarrow#pi Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);
    TH2D *hD02K_x_PTheta=new TH2D("hD02K_x_PTheta","D^{0}#rightarrowK x_{B}>0.1 Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);
    TH2D *hD02pi_x_PTheta=new TH2D("hD02pi_x_PTheta","D^{0}#rightarrow#pi x_{B}>0.1 Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);
    TH2D *hDst_PTheta=new TH2D("hDst_PTheta","D^{*+} Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);   
    TH2D *hDst2pi_PTheta=new TH2D("hDst2pi_PTheta","D^{*+}#rightarrow#pi^{+} Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);     
    TH2D *hDst2D0_PTheta=new TH2D("hDst2D0_PTheta","D^{*+}#rightarrowD^{0} Momentum vs Theta (10 fb^{-1} )",36,0.,TMath::Pi(),1000,0,5);//-2,3);//100.,0.,50.);  

    TH2F *hD0Acc = new TH2F("hD0Acc","hD0Acc",200,0,20,5,0,5);

    TH3F *Jpsi_M = new TH3F("Jpsi_M","Jpsi_M",600,2.7,3.3,100,0,20,100,-10,10);
    TH3F *D0_M = new TH3F("D0_M","D0_M",20,0,10,600,1.5,2.1,40,-4,4);
    TH3F *D0_M_2 = new TH3F("D0_M_2","D0_M_2",20,0,10,600,1.5,2.1,40,-4,4);
    TH3F *D0_M_3 = new TH3F("D0_M_3","D0_M_3",20,0,10,600,1.5,2.1,40,-4,4);
    TH3F *Dp_M = new TH3F("Dp_M","Dp_M",600,1.5,2.1,100,0,20,100,-10,10);


    TH2F *Jpsi = new TH2F("Jpsi","Jpsi",200,0,20,1000,-1,1);
    TH2F *B2Jpsi = new TH2F("B2Jpsi","B2Jpsi",200,0,20,1000,-1,1);
    TH2F *PromptJpsi = new TH2F("PromptJpsi","PromptJpsi",200,0,20,1000,-1,1);
    //TH2F *D0 = new TH2F("D0","D0",200,0,20,2000,-2,2);
    TH2F *D0_NV = new TH2F("D0_NV","D0_NV",200,0,20,2000,-2,2);
    TH2F *B2D0 = new TH2F("B2D0","B2D0",200,0,20,2000,-2,2);
    TH2F *B2D0_NV = new TH2F("B2D0_NV","B2D0_NV",200,0,20,2000,-2,2);
    TH2F *PromptD0 = new TH2F("PromptD0","PromptD0",200,0,20,2000,-2,2);    
    TH2F *PromptD0_NV = new TH2F("PromptD0_NV","PromptD0_NV",200,0,20,2000,-2,2);
    //TH2F *D0_Pass = new TH2F("D0_Pass","D0_Pass",200,0,20,2000,-2,2);
    TH3F *pres_pi = new TH3F("pres_pi","pres_pi",200,0,50,1000,-0.1,0.1,40,-4,4);
    TH3F *D0 = new TH3F("D0","D0",20,0,10,600,1.5,2.1,40,-4,4);
    TH3F *D0_Pass = new TH3F("D0_Pass","D0_Pass",20,0,10,600,1.5,2.1,40,-4,4);
    TH2F *BKGD0_Pass = new TH2F("BKGD0_Pass","BKGD0_Pass",200,0,20,2000,-2,2);
    TH2F *B2D0_Pass = new TH2F("B2D0_Pass","B2D0_Pass",200,0,20,2000,-2,2);
    TH2F *PromptD0_Pass = new TH2F("PromptD0_Pass","PromptD0_Pass",200,0,20,2000,-2,2);
    TH2F *Dp = new TH2F("Dp","Dp",200,0,20,1000,-2,2);
    TH2F *B2Dp = new TH2F("B2Dp","B2Dp",200,0,20,2000,-2,2);
    TH2F *PromptDp = new TH2F("PromptDp","PromptDp",200,0,20,2000,-2,2);
    TH2F *elec = new TH2F("elec","elec",200,0,20,10000,-10,10);
    TH2F *elec_NV = new TH2F("elec_NV","elec_NV",200,0,20,10000,-10,10);
    TH2F *B2e = new TH2F("B2e","B2e",200,0,20,10000,-10,10);
    TH2F *B2e_NV = new TH2F("B2e_NV","B2e_NV",200,0,20,10000,-10,10);
    TH2F *D2e = new TH2F("D2e","D2e",200,0,20,10000,-10,10);
    TH2F *D2e_NV = new TH2F("D2e_NV","D2e_NV",200,0,20,10000,-10,10);
    TH2F *hxQ2_B = new TH2F("hxQ2_B","hxQ2_B",100,-5,0,100,0,5);
    TH2F *hxQ2_D = new TH2F("hxQ2_D","hxQ2_D",100,-5,0,100,0,5);
    TH3F *hxQ2_D_Eta = new TH3F("hxQ2_D_Eta","hxQ2_D_Eta",100,-5,0,100,0,5,60,-3,3);
    TH2F *hxQ2_B_DIS = new TH2F("hxQ2_B_DIS","hxQ2_B_DIS",100,-5,0,100,0,5);
    TH2F *hxQ2_D_DIS = new TH2F("hxQ2_D_DIS","hxQ2_D_DIS",100,-5,0,100,0,5);
    TH2F *hxQ2_B_L0DIS = new TH2F("hxQ2_B_L0DIS","hxQ2_B_L0DIS",100,-5,0,100,0,5);
    TH2F *hxQ2_D_L0DIS = new TH2F("hxQ2_D_L0DIS","hxQ2_D_L0DIS",100,-5,0,100,0,5);
    TH2F *hxQ2_B_Ph = new TH2F("hxQ2_B_Ph","hxQ2_B_Ph",100,-5,0,100,0,5);
    TH2F *hxQ2_D_Ph = new TH2F("hxQ2_D_Ph","hxQ2_D_Ph",100,-5,0,100,0,5);
    TH2F *hxQ2_B_Re = new TH2F("hxQ2_B_Re","hxQ2_B_Re",100,-5,0,100,0,5);
    TH2F *hxQ2_D_Re = new TH2F("hxQ2_D_Re","hxQ2_D_Re",100,-5,0,100,0,5);
    TH2F *res = new TH2F("res","res",200,0,20,10000,-1,1);
    TH2F *dcaxy = new TH2F("dcaxy","dcaxy",200,0,20,1000,-5,5);
    TH2F *dcaz = new TH2F("dcaz","dcaz",200,0,20,1000,-5,5);
    TH2F *dcaz1 = new TH2F("dcaz1","dcaz1",200,0,20,1000,-5,5);
    TH2F *dcaz2 = new TH2F("dcaz2","dcaz2",200,0,20,1000,-5,5);
    TH2F *dcaz3 = new TH2F("dcaz3","dcaz3",200,0,20,1000,-5,5);

    TH2F *dca = new TH2F("dca","dca",200,0,20,1000,-5,5);
    TH2F *etapt = new TH2F("etapt","etapt",600,-4,4,200,0,20);

    TH2F *hVtx_Mult_XY = new TH2F("hVtx_Mult_XY","hVtx_Mult_XY",100,-0.5,100-0.5,1000,-1,1);
    TH2F *hVtx_Mult_Y = new TH2F("hVtx_Mult_Y","hVtx_Mult_Y",100,-0.5,100-0.5,1000,-1,1);
    TH2F *hVtx_Mult_X = new TH2F("hVtx_Mult_X","hVtx_Mult_X",100,-0.5,100-0.5,1000,-1,1);
    TH2F *hVtx_Mult_Z = new TH2F("hVtx_Mult_Z","hVtx_Mult_Z",100,-0.5,100-0.5,1000,-1,1);
    TH1F *hVtx_Mult = new TH1F("hVtx_Mult","hVtx_Mult",100,-0.5,100-0.5);
    TH1F *hVtx_Mult_Eta3 = new TH1F("hVtx_Mult_Eta3","hVtx_Mult_Eta3",100,-0.5,100-0.5);
    TH1F *hVtx_MultHF = new TH1F("hVtx_MultHF","hVtx_MultHF",100,-0.5,100-0.5);
    TH1F *hVtx_MultHF_Eta3 = new TH1F("hVtx_MultHF_Eta3","hVtx_MultHF_Eta3",100,-0.5,100-0.5);
    TH2F *hVtx_MultPt = new TH2F("hVtx_MultPt","hVtx_MultPt",100,-0.5,100-0.5,200,0,20);
    TH2F *hVtx_MultPt_Eta3 = new TH2F("hVtx_MultPt_Eta3","hVtx_MultPt_Eta3",100,-0.5,100-0.5,200,0,20);
    TH2F *hVtx_MultPtHF = new TH2F("hVtx_MultPtHF","hVtx_MultPtHF",100,-0.5,100-0.5,200,0,20);
    TH2F *hVtx_MultPtHF_Eta3 = new TH2F("hVtx_MultPtHF_Eta3","hVtx_MultPtHF_Eta3",100,-0.5,100-0.5,200,0,20);
    int iloop = ch->GetEntries();//10E6;//3.108E6;//ch->GetEntries();

    //if(VTX_FIT)iloop=5E6;
    cout << "Actually " << ch->GetEntries() << " events in files " << endl;
    TF1 *ff=new TF1("ff","1.13383e+01 + (2.15858e+01)/sqrt(x)+x*x*(-2.40998e-02)",0,100);
    TF1 *ffnew=new TF1("ffnew","53/sqrt(x)/1000.",0,100);
    TF1 *ffnewz=new TF1("ffnewz","53/sqrt(x)/1000.",0,100);
    TF1 *ff1=new TF1("ff1","1",-TMath::Pi(),TMath::Pi());
    cout <<"Starting loop of " << iloop << " events " << endl;
    for(int i =0  ;i < iloop; i++){
	ch->GetEvent(i);
	_hxQ2->Fill(log10(mx),log10(mq2));

	double WW = 1.;
	
	int bin_w = N_GLUON->FindBin(log10(mxtp),log10(mq2));
	if(mtp == 21) WW = N_GLUON->GetBinContent(bin_w);
	if(mtp == 2) WW = N_UQ->GetBinContent(bin_w);
	if(mtp == 1) WW = N_DQ->GetBinContent(bin_w);
	if(mtp == 3) WW = N_SQ->GetBinContent(bin_w);
	if(mtp == 4) WW = N_CQ->GetBinContent(bin_w);
	if(mtp == 5) WW = N_CQ->GetBinContent(bin_w);
	if(mtp == -2) WW = N_UBARQ->GetBinContent(bin_w);
	if(mtp == -1) WW = N_DBARQ->GetBinContent(bin_w);
	if(mtp == -3) WW = N_SBARQ->GetBinContent(bin_w);
	if(mtp == -4) WW = N_CBARQ->GetBinContent(bin_w);
	if(mtp == -5) WW = N_CBARQ->GetBinContent(bin_w);//Do same for B
	if(WW==0)WW=1;
	if(WW>100)WW=1;
	//cout <<" Weight " << WW << " " <<mxtp << " " << log10(mxtp) << " " << log10(mq2) << " " << mtp<< endl;
	double WWp = 1.;
	if(mtp == 21) WWp = N_GLUONp->GetBinContent(bin_w);
        if(mtp == 2) WWp = N_UQp->GetBinContent(bin_w);
        if(mtp == 1) WWp = N_DQp->GetBinContent(bin_w);
        if(mtp == 3) WWp = N_SQp->GetBinContent(bin_w);
        if(mtp == 4) WWp = N_CQp->GetBinContent(bin_w);
	if(mtp == 5) WWp = N_CQp->GetBinContent(bin_w);
        if(mtp == -2) WWp = N_UBARQp->GetBinContent(bin_w);
        if(mtp == -1) WWp = N_DBARQp->GetBinContent(bin_w);
        if(mtp == -3) WWp = N_SBARQp->GetBinContent(bin_w);
        if(mtp == -4) WWp = N_CBARQp->GetBinContent(bin_w);
        if(mtp == -5) WWp = N_CBARQp->GetBinContent(bin_w);//Do same for B                                                                                                                                                                                                  
        if(WWp==0)WWp=1;
	if(WWp>100)WWp=1;
	//Turing all parton weights off here
	WWp = 1;
	WW = 1;

	Incl_RedCS->Fill(log10(mx),log10(mq2),1.85);
	Incl_nRedCS->Fill(log10(mx),log10(mq2),1.85,WW);
	if(i%100000 == 0)cout << "On event " << i << " out of " << iloop << " " << float(i)/iloop*100 << "%" << endl;
	int bflag=0;
	int cflag=0;
	int pflag=-1;
	if(msp==131 || msp==132 || msp==135 || msp==136 )pflag=1;
	if(msp==99)pflag=4;
	if(msp==11 || msp==12 || msp==13 || msp==28 || msp==53 || msp==68 )pflag=2;
	if(msp==91 || msp==92 || msp==93 || msp==94 || msp==95 || msp==96 )pflag=3;
	if(pflag==-1)cout << "HERE "<< msp << endl;

	TVector3 vertex(0,0,0);
	double YY = my*my/(1+(1-my)*(1-my));

	//Track level loop
        //First doing primary vertex fit
	/*
	EICVertexFitter *pVtx = new EICVertexFitter("pVtx");
	if(VTX_FIT){
	    pVtx->setVerbose(false);
	    pVtx->setDebug(false);
	    pVtx->init(1);	
	    //pVtx->setSeedZ(gRandom->Gaus(0, 50));
	    }*/
	int mult=0;
	int mult3=0;
	int hf_flag=0;
	int hf_flag3=0;
	vector<float> pts;
	for(int tr = 0; tr<mntracks ; tr++){
	    double  p = sqrt(mpx[tr]*mpx[tr]+mpy[tr]*mpy[tr]+mpz[tr]*mpz[tr]);
	    double pt = sqrt(mpx[tr]*mpx[tr]+mpy[tr]*mpy[tr]);
	    if(isFromB(mid[tr])||isFromD(mid[tr])){
		if(fabs(TMath::ATanH(mpz[tr]/p))<1)hf_flag=1;
		if(fabs(TMath::ATanH(mpz[tr]/p))<3.5)hf_flag3=1;
	    }
	    if(fabs(mid[tr])==211 || fabs(mid[tr])==321 || fabs(mid[tr])==2212){
		
		if(fabs(TMath::ATanH(mpz[tr]/p))<1){
		    /*if(VTX_FIT){
			TLorentzVector *p1 = new TLorentzVector(mpx[tr],mpy[tr],mpz[tr],me[tr]);
			TVector3 p1_vtx(mvx[tr],mvy[tr],mvz[tr]);
			TLorentzVector *ps1 = smearMomEIC(p1);
			TVector3 pos = smearPosEIC(ps1,p1_vtx);
			TVector3 mom = ps1->Vect();
			pVtx->addParticle(mom,pos);
			pts.push_back(ps1->Perp());
			delete ps1;delete p1;
			}*/
		    if(!passTrackingATHENA(p,(fabs(TMath::ATanH(mpz[tr]/p)))))mult++;
		}
		if(fabs(TMath::ATanH(mpz[tr]/p))<3.5){
		    /*
		    if(VTX_FIT){
                        if(passTrackingATHENA(p,(fabs(TMath::ATanH(mpz[tr]/p))))){
			    TLorentzVector *p1 = new TLorentzVector(mpx[tr],mpy[tr],mpz[tr],me[tr]);
			    TVector3 p1_vtx(mvx[tr],mvy[tr],mvz[tr]);
			    int bin1 = 1;
			    if(fabs(p1->PseudoRapidity())<=3.5)bin1 = Res_Handler->FindBin(p1->PseudoRapidity());
			    TLorentzVector *ps1 = smearMomATHENA(p1,gmom_res[bin1-1]);
			    TVector3 pos = smearPosATHENA(ps1,p1_vtx,gdca_rphi_res[bin1-1],gdca_z_res[bin1-1]);
			    TVector3 mom = ps1->Vect();
			    pVtx->addParticle(mom,pos,gdca_rphi_res[bin1-1]->Eval(p),gdca_z_res[bin1-1]->Eval(p));
			    pts.push_back(ps1->Perp());
			    delete ps1;delete p1;
			}
		    }
		    */
		    if(passTrackingATHENA(p,(fabs(TMath::ATanH(mpz[tr]/p)))))
			mult3++;
		}
	    }
	}
/*	if(VTX_FIT){
	    //pVtx->doFit(0);
	    pVtx->doFit(1);
	    pVtx->doFit(2);
	    
	    vertex += pVtx->getVtx();
	    
	    double xyres = sqrt(vertex.X()*vertex.X()+vertex.Y()*vertex.Y());
	    hVtx_Mult_XY->Fill(mult3,xyres);
	    hVtx_Mult_X->Fill(mult3,vertex.X());
	    hVtx_Mult_Y->Fill(mult3,vertex.Y());
	    hVtx_Mult_Z->Fill(mult3,vertex.Z());
	    
	    hVtx_Mult->Fill(mult);
	    hVtx_Mult_Eta3->Fill(mult3);
	    if(hf_flag==1)hVtx_MultHF->Fill(mult);
	    if(hf_flag3==1)hVtx_MultHF_Eta3->Fill(mult3);
	    for(int j = 0;j<fabs(pts.size());j++){
		hVtx_MultPt->Fill(mult,pts[j]);
		hVtx_MultPt_Eta3->Fill(mult3,pts[j]);
		if(hf_flag==1)hVtx_MultPtHF->Fill(mult,pts[j]);
		if(hf_flag3==1)hVtx_MultPtHF_Eta3->Fill(mult3,pts[j]);
	    }
	}
	delete pVtx;
*/
	//For fast sim with no vtx fit do following
	double v_phi = ff1->GetRandom();
	//double v_r = gRandom->Gaus(0,ff->Eval(mult3))/1000.;
	//if(!VTX_FIT)vertex.SetXYZ(v_r*TMath::Cos(v_phi),v_r*TMath::Sin(v_phi),0);
	double rann = gRandom->Gaus(0,ffnew->Eval(mult3));
	double rann1 = gRandom->Gaus(0,ffnew->Eval(mult3));
	//if(!VTX_FIT)vertex.SetXYZ(gRandom->Gaus(0,ffnew->Eval(mult3)), gRandom->Gaus(0,ffnew->Eval(mult3)),gRandom->Gaus(0,ffnewz->Eval(mult3)));
	if(mult3<2)mult3=2;
	if(mult<2)mult=2;
	double newx = gRandom->Gaus(0,VertexRes_X->GetBinContent(VertexRes_X->FindBin(mult3+0.1)))/1000.;
	double newy = gRandom->Gaus(0,VertexRes_Y->GetBinContent(VertexRes_Y->FindBin(mult3+0.1)))/1000.;
	double newz = gRandom->Gaus(0,VertexRes_Z->GetBinContent(VertexRes_Z->FindBin(mult+0.1)))/1000.;
	if(!VTX_FIT)vertex.SetXYZ(newx,newy,newz);
//==========================================
	int check1=0;
	int check2=0;
	mphi = atan2(mpy[2],mpx[2]);
	if(mphi<0)mphi+=TMath::Pi()*2.;
	for(int tr = 0; tr<mntracks ; tr++){
            double  p = sqrt(mpx[tr]*mpx[tr]+mpy[tr]*mpy[tr]+mpz[tr]*mpz[tr]);
            double eta = TMath::ATanH(mpz[tr]/p);
            double theta = atan2(sqrt(mpx[tr]*mpx[tr]+mpy[tr]*mpy[tr]),mpz[tr]);
	    double phi =  atan2(mpy[tr],mpx[tr]);
	    if(phi<0)phi+=TMath::Pi()*2.;     
//if(theta<0)theta+=TMath::Pi()*2.;                                                                                                                                                                                                                                     //theta = theta*180./TMath::Pi();                                                                                                                                                                                                                       
            //if(fabs(mid[tr])==421)cout <<p << " " << theta << " " << eta << endl;                                                                                                                                                                                  
            //Mother particle ID                                                                                                                                                                                                                                              
            int motherID = 0;
            if(mpid[tr]>0)motherID=mid[mpid[tr]-1];//mother id                                                                                                                                                                                                                
            //Grandmother particle ID                                                                                                                                                                                                                                         
            int gmotherLine = mpid[mpid[tr]-1];//Find if grandmother particle exists                                                                                                                                                                                          
            int gmotherID=0;//Default no grandmother                                                                                                                                                                                                                          
            if(gmotherLine>0)gmotherID = mid[gmotherLine-1];//Has grandmother                                                                                                                                                                                                 
            //Great-grandmother particle ID                                                                                                                                                                                                                                   
            int ggmotherLine=-1;
            if(gmotherLine>0)ggmotherLine=mpid[mpid[mpid[tr]-1]-1];//Find if great-grandmother particle exists (to account for B->D*->D->e)                                                                                                                                   
            int ggmotherID=0;//Default no great-grandmother                                                                                                                                                                                                                   
            if(ggmotherLine>0)ggmotherID = mid[ggmotherLine-1];//Has great-grandmother                                                                                                                                                                                        
            //Great-great-grandmother particle ID                                                                                                                                                                                                                             
            int gggmotherLine=-1;
            if(ggmotherLine>0)gggmotherLine=mpid[mpid[mpid[mpid[tr]-1]-1]-1];//Find if great-great-grandmother particle exists (to account for B->D**->D*->D->e)                                                                                                              
            int gggmotherID=0;//Default no great-great-grandmother                                                                                                                                                                                                            
            if(gggmotherLine>0)gggmotherID = mid[gggmotherLine-1];//Has great-great-grandmother                                                                                                                                                                               
            //g-g-g-grandmother particle ID                                                                                                                                                                                                                                   
            int ggggmotherLine=-1;//Default no great-great-grandmother                                                                                                                                                                                                        
            if(gggmotherLine>0)ggggmotherLine=mpid[mpid[mpid[mpid[mpid[tr]-1]-1]-1]-1];//Find if great-great-grandmother particle exists (to account for higher B->D**->D*->D->e)                                                                                             
            int ggggmotherID=0;
            if(ggggmotherLine>0)ggggmotherID = mid[ggggmotherLine-1];//Has great-great-grandmother                                                                                                                                                                            
            //g-g-g-g-grandmother particle ID                                                                                                                                                                                                                                 
            int gggggmotherLine=-1;//Default no great-great-grandmother                                                                                                                                                                                                       
            if(ggggmotherLine>0)gggggmotherLine=mpid[mpid[mpid[mpid[mpid[mpid[tr]-1]-1]-1]-1]-1];//Find if great-great-grandmother particle exists (to account for higher B->D**->D*->D->e)                                                                                   
            int gggggmotherID=0;
            if(gggggmotherLine>0)gggggmotherID = mid[gggggmotherLine-1];//Has great-great-grandmother                                                                                                                                                                         
            int flag = 0;
            flag+=isFromB(motherID);
            flag+=isFromB(gmotherID);
            flag+=isFromB(ggmotherID);
            flag+=isFromB(gggmotherID);
            flag+=isFromB(ggggmotherID);
            flag+=isFromB(gggggmotherID);
	    if(onlyD0 && flag!=0)continue;

	    if(fabs(mid[tr])==211 || fabs(mid[tr])==321 || fabs(mid[tr])==2212){
                int bkg_flag=0;
		TLorentzVector *p1 = new TLorentzVector(mpx[tr],mpy[tr],mpz[tr],me[tr]);
                if(sqrt(mpx[tr]*mpx[tr]+mpy[tr]*mpy[tr])>0){
		    //TLorentzVector *ps1 = (TLorentzVector *)smearMomEIC2(p1,MomResEtaVsP,MomResEtaVsP_K,fabs(mid[tr]));
		    int bin1 = 1;
                    if(fabs(p1->PseudoRapidity())<=3.5)bin1 = Res_Handler->FindBin(p1->PseudoRapidity());
  
                    TLorentzVector *ps1 = smearMomATHENA(p1,gmom_res[bin1-1]);
		    //TLorentzVector *ps1 = smearMomEIC(p1);  
		    TVector3 p1_vtx(mvx[tr],mvy[tr],mvz[tr]);
		    
		    //TVector3 ps1_vtx = smearPosEIC3(ps1,p1_vtx,hDCATRes_Plus[getPtIdx(ps1->Pt())][getEtaIdx(ps1->PseudoRapidity())],hDCATRes_Minus[getPtIdx(ps1->Pt())][getEtaIdx(ps1->PseudoRapidity())],DCARes_Z_Eta1,DCARes_Z_Eta2,DCARes_Z_Eta3,mid[tr]);
		    //TVector3 ps1_vtx = smearPosEIC(ps1,p1_vtx);      
		    TVector3 ps1_vtx = smearPosATHENA(ps1,p1_vtx,gdca_rphi_res[bin1-1],gdca_z_res[bin1-1]);
		    for(int tr1 = tr+1; tr1<mntracks ; tr1++){
			if(!(fabs(mid[tr1])==211 || fabs(mid[tr1])==321 || fabs(mid[tr1])==2212))continue;
			if(sqrt(mpx[tr1]*mpx[tr1]+mpy[tr1]*mpy[tr1])==0)continue;
			/*if(
			  !(
			    (fabs(mid[tr])==211 && fabs(mid[tr1])==321) ||
				(fabs(mid[tr1])==211 && fabs(mid[tr])==321)
				)
				){
			    continue; //any sign k-pi pair                                                                                                                                                                   
			    }*/
			TLorentzVector *p2 = new TLorentzVector(mpx[tr1],mpy[tr1],mpz[tr1],me[tr1]);
			if(!passPID(mid[tr],mid[tr1],p1->P(),p2->P(),p1->PseudoRapidity(),p2->PseudoRapidity())){
			    delete p2;
			    continue;
			}
			if(mid[tr]*mid[tr1]>0){
			    bkg_flag=2;
			}
			else{
			    bkg_flag=1;
			}
			//TLorentzVector *ps2 = (TLorentzVector *)smearMomEIC2(p2,MomResEtaVsP,MomResEtaVsP_K,fabs(mid[tr1]));
			int bin2 = 1;
			if(fabs(p2->PseudoRapidity())<=3.5)bin2 = Res_Handler->FindBin(p2->PseudoRapidity());
			TLorentzVector *ps2 = smearMomATHENA(p2,gmom_res[bin2-1]);
			//TLorentzVector *ps2 = smearMomEIC(p2);  
			TVector3 p2_vtx(mvx[tr1],mvy[tr1],mvz[tr1]);
			//TVector3 ps2_vtx = smearPosEIC3(ps2,p2_vtx,hDCATRes_Plus[getPtIdx(ps2->Pt())][getEtaIdx(ps2->PseudoRapidity())],hDCATRes_Minus[getPtIdx(ps2->Pt())][getEtaIdx(ps2->PseudoRapidity())],DCARes_Z_Eta1,DCARes_Z_Eta2,DCARes_Z_Eta3,mid[tr1]);
			//TVector3 ps2_vtx = smearPosEIC(ps2,p2_vtx);      
			TVector3 ps2_vtx = smearPosATHENA(ps2,p2_vtx,gdca_rphi_res[bin2-1],gdca_z_res[bin2-1]);
			TLorentzVector *Ps = new TLorentzVector(*ps1+*ps2);
			double d0_m = Ps->M();
			int dofill = 0;
			TVector3 Ps_vtx = (ps1_vtx+ps2_vtx)*0.5;
			double _dca = dcaXY(Ps->Vect(), Ps_vtx, vertex,mult3);
			if(!passTrackingATHENA(p1->P(),fabs(p1->PseudoRapidity())) || !passTrackingATHENA(p2->P(),fabs(p2->PseudoRapidity()))){
			    delete p2; delete ps2; delete Ps;
			    continue;
			}
			if(bkg_flag==1 && fabs(Ps->PseudoRapidity())<3.5 && fabs(Ps->M()-1.86484)<0.02)D0_RedCS_BKG->Fill(log10(mx),log10(mq2),Ps->M());
			if(bkg_flag==2)dofill=2;
			if(bkg_flag==1)D0_M->Fill(Ps->Pt(),Ps->M(),Ps->PseudoRapidity());
			if(isGoodD0_XY(Ps,ps1->Vect(),ps1_vtx,ps2->Vect(),ps2_vtx,vertex,fabs(mid[tr]),fabs(mid[tr1]),TrackEff_EtaP,dofill,ps1->E(),ps2->E())){
			    if(bkg_flag==1 && fabs(Ps->PseudoRapidity())<3.5 && fabs(Ps->M()-1.86484)<0.02){
				check1=1;
				double dphi = Ps->Phi();
				if(dphi<0)dphi+=TMath::Pi()*2.;
				double delta_phi = fabs(dphi - mphi)/TMath::Pi();
				if(delta_phi>1) delta_phi =1 - delta_phi;
				BKGD0_Pass->Fill(Ps->Pt(),_dca);
				D0_RedCS_Topo_BKG->Fill(log10(mx),log10(mq2),Ps->M());
				D0_nRedCS_Topo_BKG->Fill(log10(mx),log10(mq2),Ps->M(),WW);
				D0L_RedCS_Topo_BKG->Fill(log10(mx),log10(mq2),delta_phi,WWp);
                                D0L_nRedCS_Topo_BKG->Fill(log10(mx),log10(mq2),delta_phi,WW*WWp);
				D0_DCA_RedCS_Topo_BKG->Fill(log10(mx),log10(mq2),_dca);
			    }
			    if(bkg_flag==1)D0_M_2->Fill(Ps->Pt(),Ps->M(),Ps->PseudoRapidity());
			    else if(bkg_flag==2)D0_M_3->Fill(Ps->Pt(),Ps->M(),Ps->PseudoRapidity());
			}
			delete p2; delete ps2; delete Ps;
		    }
		    delete ps1;
		}
		delete p1;
	    }
	    	    
//if(isFromB(gggmotherID))cout << mid[tr] << " mother " << motherID << "  " << gmotherID << " " << ggmotherID << " " << gggmotherID << " " << ggggmotherID << " " << ggggmotherID <<" " << gggggmotherID<<endl; 
	    //TVector3 vertex(0,0,0);
	    if(mid[tr]==421)hD0_PTheta->Fill(theta,log10(p*100.));
            if(mid[tr]==-421)hD0bar_PTheta->Fill(theta,log10(p*100.));
            if(mid[tr]==411)hDp_PTheta->Fill(theta,log10(p*100.));
            if(mid[tr]==-411)hDpbar_PTheta->Fill(theta,log10(p*100.));
            if(mid[tr]==431)hDs_PTheta->Fill(theta,log10(p*100.));
            if(mid[tr]==-431)hDsbar_PTheta->Fill(theta,log10(p*100.));
            if(mid[tr]==4122)hLc_PTheta->Fill(theta,log10(p*100.));
            if(mid[tr]==-4122)hLcbar_PTheta->Fill(theta,log10(p*100.));
	    if(mid[tr]==443)hJpsi_PTheta->Fill(theta,log10(p*100.));
            if(mid[tr]==511)hB0_PTheta->Fill(theta,log10(p*100.));
            if(mid[tr]==-511)hB0bar_PTheta->Fill(theta,log10(p*100.));
            if(mid[tr]==521)hBp_PTheta->Fill(theta,log10(p*100.));
            if(mid[tr]==-521)hBpbar_PTheta->Fill(theta,log10(p*100.));
            if(mid[tr]==531)hBs_PTheta->Fill(theta,log10(p*100.));
            if(mid[tr]==-531)hBsbar_PTheta->Fill(theta,log10(p*100.));
            if(mid[tr]==5122)hLb_PTheta->Fill(theta,log10(p*100.));
            if(mid[tr]==-5122)hLbbar_PTheta->Fill(theta,log10(p*100.));

	    if(fabs(mid[tr])==11 && motherID==443)hJpsi2e_PTheta->Fill(theta,log10(p*100.));
            if(fabs(mid[tr])==13 && motherID==443)hJpsi2mu_PTheta->Fill(theta,log10(p*100.));
	    if(fabs(mid[tr])==421 && flag>0)hB2D0_PTheta->Fill(theta,log10(p*100.));
	    if(fabs(mid[tr])==11 && isFromB(motherID))hB2e_PTheta->Fill(theta,log10(p*100.));
            if(fabs(mid[tr])==321 && fabs(motherID)==421)hD02K_PTheta->Fill(theta,log10(p*100.));
            if(fabs(mid[tr])==211 && fabs(motherID)==421)hD02pi_PTheta->Fill(theta,log10(p*100.));
            if(fabs(mid[tr])==321 && fabs(motherID)==421 && mx>0.1)hD02K_x_PTheta->Fill(theta,log10(p*100.));
            if(fabs(mid[tr])==211 && fabs(motherID)==421 && mx>0.1)hD02pi_x_PTheta->Fill(theta,log10(p*100.));

            if(fabs(mid[tr])==413)hDst_PTheta->Fill(theta,log10(p*100.));
            if(fabs(mid[tr])==211 && fabs(motherID)==413)hDst2pi_PTheta->Fill(theta,log10(p*100.));
            if(fabs(mid[tr])==421 && fabs(motherID)==413)hDst2D0_PTheta->Fill(theta,log10(p*100.));

	    // Doing some correaltion studies between D0 and scattered lepton
	    if(fabs(mid[tr])==421 && fabs(eta<3.5)){
		double delta_phi = fabs(phi - mphi)/TMath::Pi();
		if(delta_phi>1) delta_phi =1 - delta_phi;
		if(pflag==1)D0_Lepton_Corr_PGF->Fill(log10(mx),log10(mq2),delta_phi);
		else if(pflag==4)D0_Lepton_Corr_L0->Fill(log10(mx),log10(mq2),delta_phi);
		else D0_Lepton_Corr_Other->Fill(log10(mx),log10(mq2),delta_phi);
	    }
	    // End corr studies
	    /*
	    if(0){//fabs(eta<3) && fabs(mid[tr])==211 && (fabs(motherID)==2 || fabs(motherID)==1)){//!(isFromB(motherID) || isFromD(motherID))){
		TLorentzVector *p1 = new TLorentzVector(mpx[tr],mpy[tr],mpz[tr],me[tr]);
		TVector3 p1_vtx(mvx[tr],mvy[tr],mvz[tr]);
		TLorentzVector * ps1 = smearMomEIC2(p1,MomResEtaVsP,MomResEtaVsP_K,fabs(mid[tr]));
		TVector3 ps1_vtx = smearPosEIC2(ps1,p1_vtx,DCARes_T_Eta1,DCARes_T_Eta2,DCARes_T_Eta3,DCARes_Z_Eta1,DCARes_Z_Eta2,DCARes_Z_Eta3);
		double _dca = dcaSigned(ps1->Vect(), ps1_vtx, vertex);
                double _dcaxy = dcaXY(ps1->Vect(), ps1_vtx, vertex,mult3);
                double _dcaz = dcaZ(ps1->Vect(), ps1_vtx, vertex);
		dca->Fill(p1->Pt(),_dca);
		dcaxy->Fill(p1->Pt(),_dcaxy);
		dcaz->Fill(p1->Pt(),_dcaz);
		if(fabs(eta<1))dcaz1->Fill(p1->Pt(),_dcaz);
		else if (fabs(eta<2))dcaz2->Fill(p1->Pt(),_dcaz);
		else if (fabs(eta<3))dcaz3->Fill(p1->Pt(),_dcaz);
		res->Fill(p1->Pt(),(p1->Pt()-ps1->Pt())/p1->Pt());
		etapt->Fill(eta,p1->Pt());
		delete p1; delete ps1;
	    }
	    */
	    if(isFromB(mid[tr]))bflag=1;
	    if(isFromD(mid[tr]) && !(isFromB(motherID) || isFromB(gmotherID)|| isFromB(ggmotherID)|| isFromB(gggmotherID)|| isFromB(ggggmotherID)|| isFromB(gggggmotherID)) )cflag=1;
	    //Looking at electron channel first
	    /*
	    if(fabs(mid[tr])==11){
		if(isFromB(motherID) || isFromD(motherID)){// FROM HF DECAY
		    //cout << "here " << motherID << " " << gmotherID << " " << gmotherLine << " ggm " <<  ggmotherID << " " << ggmotherLine << endl;
		    TLorentzVector *p1 = new TLorentzVector(mpx[tr],mpy[tr],mpz[tr],me[tr]);
		    if(fabs(p1->PseudoRapidity())<5){
			TVector3 p1_vtx(mvx[tr],mvy[tr],mvz[tr]);
			TLorentzVector * ps1 = smearMomEIC2(p1,MomResEtaVsP,MomResEtaVsP_K,fabs(mid[tr]));
			TVector3 ps1_vtx = smearPosEIC2(ps1,p1_vtx,DCARes_T_Eta1,DCARes_T_Eta2,DCARes_T_Eta3,DCARes_Z_Eta1,DCARes_Z_Eta2,DCARes_Z_Eta3);
			if(fabs(ps1->PseudoRapidity())<1){
			    double _dca = dcaXY(ps1->Vect(), ps1_vtx, vertex,mult3);
			    TVector3 vertex_ns(0,0,0);
			    double _dca_ns = dcaXY(ps1->Vect(), ps1_vtx, vertex_ns,mult3);
			    
			    elec->Fill(ps1->Pt(),_dca);
			    elec_NV->Fill(ps1->Pt(),_dca_ns);
			    if(isFromB(motherID) || isFromB(gmotherID) || isFromB(ggmotherID)){
				B2e->Fill(ps1->Pt(),_dca);
				B2e_NV->Fill(ps1->Pt(),_dca_ns);
			    }
			    else if(isFromD(motherID) && !(isFromB(motherID) || isFromB(gmotherID)|| isFromB(ggmotherID))){
				D2e->Fill(ps1->Pt(),_dca);
				D2e_NV->Fill(ps1->Pt(),_dca_ns);
			    }
			}
			delete ps1;
		    }
		    delete p1;
		}
	    }
	    */
	    //Now looking for Hadron channels
	    if(isFlag(fabs(mid[tr]),particle_flag)){
		//Find num. daughter particles
		if(fabs(mid[tr])==421)hxQ2_D_Eta->Fill(log10(mx),log10(mq2),eta);
		int num = fabs(mfd[tr] - mld[tr]);
		if(num==1 && (fabs(mid[tr]) == 421||fabs(mid[tr]) == 443)){//Two-body decays
		    if(!iseeDecay(mid[mfd[tr]-1],mid[mld[tr]-1]) && !isKPiDecay(mid[mfd[tr]-1],mid[mld[tr]-1])) continue;
		    TLorentzVector *p1 = new TLorentzVector(mpx[mfd[tr]-1],mpy[mfd[tr]-1],mpz[mfd[tr]-1],me[mfd[tr]-1]);
		    TLorentzVector *p2 = new TLorentzVector(mpx[mld[tr]-1],mpy[mld[tr]-1],mpz[mld[tr]-1],me[mld[tr]-1]);
		    TVector3 p1_vtx(mvx[mfd[tr]-1],mvy[mfd[tr]-1],mvz[mfd[tr]-1]);
		    TVector3 p2_vtx(mvx[mld[tr]-1],mvy[mld[tr]-1],mvz[mld[tr]-1]);
		    TLorentzVector *P = new TLorentzVector(*p1+*p2);
		    //TLorentzVector *ps1 = smearMomEIC2(p1,MomResEtaVsP,MomResEtaVsP_K,mid[mfd[tr]-1]);
		    //TLorentzVector *ps2 = smearMomEIC2(p2,MomResEtaVsP,MomResEtaVsP_K,mid[mld[tr]-1]);

		    //Below method is done to not have to touch later code where I cut on track eta
		    int bin1 = 1;
		    if(fabs(p1->PseudoRapidity())<=3.5)bin1 = Res_Handler->FindBin(p1->PseudoRapidity());
		    int bin2 = 1;
		    if(fabs(p2->PseudoRapidity())<=3.5)bin2 = Res_Handler->FindBin(p2->PseudoRapidity());

		    TLorentzVector *ps1 = smearMomATHENA(p1,gmom_res[bin1-1]);  
		    TLorentzVector *ps2 = smearMomATHENA(p2,gmom_res[bin2-1]);
		    TLorentzVector *Ps = new TLorentzVector(*ps1+*ps2);
		    TVector3 P_vtx = (p1_vtx+p2_vtx)*0.5;
		    if(fabs(mid[tr]) == 421){// && pflag==1){
			if(fabs(P->PseudoRapidity())<1){
			    hD0Acc->Fill(P->Pt(),0.5);
			    if(p1->Pt()>0.035 && p2->Pt()>0.035)hD0Acc->Fill(P->Pt(),1.5);
			    if(p1->Pt()>0.071 && p2->Pt()>0.071)hD0Acc->Fill(P->Pt(),2.5);
			    if(p1->Pt()>0.097 && p2->Pt()>0.097)hD0Acc->Fill(P->Pt(),3.5);
			    if(p1->Pt()>0.194 && p2->Pt()>0.194)hD0Acc->Fill(P->Pt(),4.5);
			}
			double dphi = P->Phi();
			if(dphi<0)dphi+=TMath::Pi()*2.;
			double delta_phi = fabs(dphi - mphi)/TMath::Pi();
			if(delta_phi>1) delta_phi =1 - delta_phi;
			if(flag==0)D0_RedCS_T_Acc->Fill(log10(mx),log10(mq2),Ps->M());
			if(flag==0)D0_nRedCS_T_Acc->Fill(log10(mx),log10(mq2),Ps->M(),WW);
			if(flag==0)D0L_RedCS_T_Acc->Fill(log10(mx),log10(mq2),delta_phi,WWp);
                        if(flag==0)D0L_nRedCS_T_Acc->Fill(log10(mx),log10(mq2),delta_phi,WW*WWp);
			else if(flag>0)B2D0_RedCS_T_Acc->Fill(log10(mx),log10(mq2),dcaXY(P->Vect(), P_vtx, vertex,mult3));
		    }
		    
		    if(fabs(Ps->PseudoRapidity())>5 || !passTrackingATHENA(p1->P(),fabs(p1->PseudoRapidity())) || !passTrackingATHENA(p2->P(),fabs(p2->PseudoRapidity()))){
			delete p1; delete p2; delete Ps;delete ps1; delete ps2;delete P;
			continue;
		    }
		    if(!passPID(mid[mfd[tr]-1],mid[mld[tr]-1],p1->P(),p2->P(),p1->PseudoRapidity(),p2->PseudoRapidity())){
			delete p1; delete p2; delete Ps;delete ps1; delete ps2;delete P;
			continue;
		    }

		    if(fabs(mid[tr]) == 421 && fabs(Ps->M()-1.86484)<0.02){
			if(fabs(Ps->PseudoRapidity())<3.5)D0_RedCS_T->Fill(log10(mx),log10(mq2),Ps->M());
		    }

		    TVector3 ps1_vtx = smearPosATHENA(ps1,p1_vtx,gdca_rphi_res[bin1-1],gdca_z_res[bin1-1]);      
		    TVector3 ps2_vtx = smearPosATHENA(ps2,p2_vtx,gdca_rphi_res[bin2-1],gdca_z_res[bin2-1]);      
		    TVector3 Ps_vtx = (ps1_vtx+ps2_vtx)*0.5;
		    double _dca = dcaXY(Ps->Vect(), Ps_vtx, vertex,mult3);
		    double pi_dca = -1;
		    TVector3 vertex_ns(0,0,0);
		    if(fabs(mid[tr]) == 421){
			if(fabs(mid[mfd[tr]-1])==211)pi_dca = fabs(dcaXY(p1->Vect(),p1_vtx,vertex_ns,mult3));
			else pi_dca = fabs(dcaXY(p2->Vect(),p2_vtx,vertex_ns,mult3));
			hD02pi_DCATheta->Fill(theta,pi_dca);
			if(pi_dca>1)cout << "greater DCA for pion " << endl;
		    }
		    double _dca_ns = dcaXY(Ps->Vect(), Ps_vtx, vertex_ns,mult3);
		    bool D0flag = false;
		    if(fabs(mid[tr]) == 421)D0flag=isGoodD0_XY(Ps,ps1->Vect(),ps1_vtx,ps2->Vect(),ps2_vtx,vertex,fabs(mid[mfd[tr]-1]),fabs(mid[mld[tr]-1]),TrackEff_EtaP,1,ps1->E(),ps2->E());
		    if(fabs(mid[tr]) == 421){
			if(isKPiDecay(mid[mfd[tr]-1],mid[mld[tr]-1])){
			    if(mid[mfd[tr]-1]==211)pres_pi->Fill(p1->P(),(p1->P()-ps1->P())/p1->P(),p1->PseudoRapidity());
			    else pres_pi->Fill(p2->P(),(p2->P()-ps2->P())/p2->P(),p2->PseudoRapidity());
			    D0->Fill(Ps->Pt(),Ps->M(),Ps->PseudoRapidity());
			    D0_NV->Fill(Ps->Pt(),_dca_ns);
			    if(flag>0){
				B2D0->Fill(Ps->Pt(),_dca);
				B2D0_NV->Fill(P->Pt(),_dca_ns);
			    }			
			    else if(flag==0){
				PromptD0->Fill(Ps->Pt(),_dca);
				PromptD0_NV->Fill(Ps->Pt(),_dca_ns);
				
			    }
			    if(D0flag){
				if(fabs(Ps->PseudoRapidity())<3.5 && fabs(Ps->M()-1.86484)<0.02){
				    check2=1;
				    double dphi = Ps->Phi();
				    if(dphi<0)dphi+=TMath::Pi()*2.;
				    double delta_phi = fabs(dphi - mphi)/TMath::Pi();
				    if(delta_phi>1) delta_phi = 1 - delta_phi;
				    D0_RedCS_Topo_T->Fill(log10(mx),log10(mq2),Ps->M());
				    D0_nRedCS_Topo_T->Fill(log10(mx),log10(mq2),Ps->M(),WW);
				    D0L_RedCS_Topo_T->Fill(log10(mx),log10(mq2),delta_phi,WWp);
                                    D0L_nRedCS_Topo_T->Fill(log10(mx),log10(mq2),delta_phi,WW*WWp);
				    //cout << "here " << log10(mx) << " " << log10(mq2) << " " << delta_phi << " " << dphi << " "<< mphi << endl;
				    if(pflag==4)D0L_RedCS_Topo_T_L0->Fill(log10(mx),log10(mq2),delta_phi,WWp);
				    else if(pflag==1)D0L_RedCS_Topo_T_PGF->Fill(log10(mx),log10(mq2),delta_phi,WWp);
				    else D0L_RedCS_Topo_T_Other->Fill(log10(mx),log10(mq2),delta_phi,WWp);
				    if(pflag==4)D0L_nRedCS_Topo_T_L0->Fill(log10(mx),log10(mq2),delta_phi,WW*WWp);
                                    else if(pflag==1)D0L_nRedCS_Topo_T_PGF->Fill(log10(mx),log10(mq2),delta_phi,WW*WWp);
                                    else D0L_nRedCS_Topo_T_Other->Fill(log10(mx),log10(mq2),delta_phi,WW*WWp);
				    if(flag>0)D0_DCA_RedCS_Topo_Bottom->Fill(log10(mx),log10(mq2),_dca);
				    else if(flag==0)D0_DCA_RedCS_Topo_Charm->Fill(log10(mx),log10(mq2),_dca);
				    
				    if(fabs(Ps->PseudoRapidity())<10){
					D0_Pass->Fill(Ps->Pt(),Ps->M(),Ps->PseudoRapidity());
					if(flag>0)B2D0_Pass->Fill(Ps->Pt(),_dca);
					else if(flag==0)PromptD0_Pass->Fill(Ps->Pt(),_dca);
				    }
				}
			    }
			}
		    }
		    if(fabs(mid[tr]) == 443){
			if(iseeDecay(mid[mfd[tr]-1],mid[mld[tr]-1])){
			    Jpsi_M->Fill(Ps->M(),Ps->P(),Ps->PseudoRapidity());
			    Jpsi->Fill(Ps->Pt(),_dca);
			    if(flag>0)B2Jpsi->Fill(Ps->Pt(),_dca);
			    else if(flag==0)PromptJpsi->Fill(Ps->Pt(),_dca);
			}
		    }
		    delete Ps;delete P;
                    delete ps2; delete ps1;
                    delete p2; delete p1;
		}
		/*
		else if(0){//num==2 && fabs(mid[tr]) == 411){//Three-body decays
		    if(!isK2PiDecay(mid[mfd[tr]-1],mid[mfd[tr]],mid[mld[tr]-1]))continue;                                                                                                         
                    TLorentzVector *p1= new TLorentzVector(mpx[mfd[tr]-1],mpy[mfd[tr]-1],mpz[mfd[tr]-1],me[mfd[tr]-1]);
		    if(fabs(p1->PseudoRapidity())>5.25)continue;
                    TLorentzVector *p2 = new TLorentzVector(mpx[mfd[tr]],mpy[mfd[tr]],mpz[mfd[tr]],me[mfd[tr]]);
		    if(fabs(p2->PseudoRapidity())>5.25)continue;
		    TLorentzVector *p3 = new TLorentzVector(mpx[mld[tr]-1],mpy[mld[tr]-1],mpz[mld[tr]-1],me[mld[tr]-1]);
		    if(fabs(p3->PseudoRapidity())>5.25)continue;
		    TVector3 p1_vtx(mvx[mfd[tr]-1],mvy[mfd[tr]-1],mvz[mfd[tr]-1]);
                    TVector3 p2_vtx(mvx[mfd[tr]],mvy[mfd[tr]],mvz[mfd[tr]]);
		    TVector3 p3_vtx(mvx[mld[tr]-1],mvy[mld[tr]-1],mvz[mld[tr]-1]);
                    TLorentzVector *P= new TLorentzVector(*p1+*p2+*p3);
		    TLorentzVector *ps1 = smearMomEIC2(p1,MomResEtaVsP,MomResEtaVsP_K,fabs(mid[tr]));
                    TLorentzVector *ps2 = smearMomEIC2(p2,MomResEtaVsP,MomResEtaVsP_K,fabs(mid[tr]));
		    TLorentzVector *ps3 = smearMomEIC2(p3,MomResEtaVsP,MomResEtaVsP_K,fabs(mid[tr]));
		    if(fabs(ps1->PseudoRapidity())>1 || fabs(ps2->PseudoRapidity())>1 || fabs(ps3->PseudoRapidity())>1)continue;
		    TLorentzVector *Ps = new TLorentzVector(*ps1+*ps2+*ps3);
                    TVector3 ps1_vtx = smearPosEIC2(ps1,p1_vtx,DCARes_T_Eta1,DCARes_T_Eta2,DCARes_T_Eta3,DCARes_Z_Eta1,DCARes_Z_Eta2,DCARes_Z_Eta3);
                    TVector3 ps2_vtx = smearPosEIC2(ps2,p2_vtx,DCARes_T_Eta1,DCARes_T_Eta2,DCARes_T_Eta3,DCARes_Z_Eta1,DCARes_Z_Eta2,DCARes_Z_Eta3);
		    TVector3 ps3_vtx = smearPosEIC2(ps3,p3_vtx,DCARes_T_Eta1,DCARes_T_Eta2,DCARes_T_Eta3,DCARes_Z_Eta1,DCARes_Z_Eta2,DCARes_Z_Eta3);
                    TVector3 Ps_vtx = (ps1_vtx+ps2_vtx+ps3_vtx)*0.3333333;
                    double _dca = dcaXY(Ps->Vect(), Ps_vtx, vertex,mult3);
		    if(fabs(mid[tr]) == 411){
			Dp_M->Fill(Ps->M(),Ps->P(),Ps->PseudoRapidity());
			Dp->Fill(Ps->Pt(),_dca);
                        if(flag>0)B2Dp->Fill(Ps->Pt(),_dca);
                        else if(flag==0)PromptDp->Fill(Ps->Pt(),_dca);
                    }
		    delete Ps;delete P;
		    delete ps3;delete ps2; delete ps1;
		    delete p3;delete p2; delete p1;
		}
		*/
		else{
		    continue;
		}
	    }
	}//track loop
	
	if(bflag==1)hxQ2_B->Fill(log10(mx),log10(mq2));
	if(cflag==1)hxQ2_D->Fill(log10(mx),log10(mq2));
	if(pflag==1){
	    if(bflag==1)hxQ2_B_DIS->Fill(log10(mx),log10(mq2));
	    if(cflag==1)hxQ2_D_DIS->Fill(log10(mx),log10(mq2));
	}
	if(pflag==4){
            if(bflag==1)hxQ2_B_L0DIS->Fill(log10(mx),log10(mq2));
            if(cflag==1)hxQ2_D_L0DIS->Fill(log10(mx),log10(mq2));
        }
	if(pflag==3){
            if(bflag==1)hxQ2_B_Ph->Fill(log10(mx),log10(mq2));
            if(cflag==1)hxQ2_D_Ph->Fill(log10(mx),log10(mq2));
        }
        if(pflag==2){
            if(bflag==1)hxQ2_B_Re->Fill(log10(mx),log10(mq2));
            if(cflag==1)hxQ2_D_Re->Fill(log10(mx),log10(mq2));
        }
    }//event loop
    if(1){
	D0_Lepton_Corr_PGF->Write();
	D0_Lepton_Corr_L0->Write();
	D0_Lepton_Corr_Other->Write();
	hD02pi_DCATheta->Write();
	hD0Acc->Write();
	Incl_RedCS->Write();
	Incl_nRedCS->Write();
	D0_RedCS_Topo_BKG->Write();
	D0_nRedCS_Topo_BKG->Write();
	D0_RedCS_BKG->Write();
	D0_DCA_RedCS_Topo_BKG->Write();
	D0_DCA_RedCS_Topo_Charm->Write();
	D0_DCA_RedCS_Topo_Bottom->Write();
	D0_RedCS_T->Write();
	D0_RedCS_T_Acc->Write();
	D0_nRedCS_T_Acc->Write();
	B2D0_RedCS_T_Acc->Write();
	D0_RedCS_Topo_T->Write();
	D0_nRedCS_Topo_T->Write();

	D0L_nRedCS_Topo_T->Write();
	D0L_RedCS_Topo_T->Write();
	D0L_nRedCS_Topo_T_PGF->Write();
        D0L_RedCS_Topo_T_PGF->Write();

	D0L_nRedCS_Topo_T_L0->Write();
        D0L_RedCS_Topo_T_L0->Write();
	D0L_nRedCS_Topo_T_Other->Write();
        D0L_RedCS_Topo_T_Other->Write();


	D0L_RedCS_Topo_BKG->Write();
	D0L_nRedCS_Topo_BKG->Write();
	D0L_nRedCS_T_Acc->Write();
	D0L_RedCS_T_Acc->Write();

	pres_pi->Write();
	hD0_PTheta->Write();
	hD0bar_PTheta->Write();
	hDp_PTheta->Write();
        hDpbar_PTheta->Write();
	hDs_PTheta->Write();
        hDsbar_PTheta->Write();
	hLc_PTheta->Write();
        hLcbar_PTheta->Write();
        hJpsi_PTheta->Write();
	hJpsi2e_PTheta->Write();
	hJpsi2mu_PTheta->Write();
	hB0_PTheta->Write();
        hB0bar_PTheta->Write();
        hBp_PTheta->Write();
	hBpbar_PTheta->Write();
        hBs_PTheta->Write();
	hBsbar_PTheta->Write();
        hLb_PTheta->Write();
	hLbbar_PTheta->Write();
	hB2D0_PTheta->Write();
	hB2e_PTheta->Write();
	hD02K_PTheta->Write();
	hD02pi_PTheta->Write();
	hD02K_x_PTheta->Write();
	hD02pi_x_PTheta->Write();
	hDst_PTheta->Write();
        hDst2D0_PTheta->Write();
        hDst2pi_PTheta->Write();
	D0_M->Write();
	D0_M_2->Write();
	D0_M_3->Write();
	Dp_M->Write();
	Jpsi_M->Write();
	D0->Write();
	D0_NV->Write();
	B2D0->Write();
	B2D0_NV->Write();
	PromptD0->Write();
	PromptD0_NV->Write();
	D0_Pass->Write();
        BKGD0_Pass->Write();
	B2D0_Pass->Write();
        PromptD0_Pass->Write();
	Dp->Write();
	B2Dp->Write();
	PromptDp->Write();
	Jpsi->Write();
	B2Jpsi->Write();
	PromptJpsi->Write();
	elec->Write();
	B2e->Write();
	D2e->Write();
	elec_NV->Write();
        B2e_NV->Write();
        D2e_NV->Write();
	hxQ2_B->Write();
	hxQ2_D->Write();
	hxQ2_B_DIS->Write();
        hxQ2_D_DIS->Write();
	hxQ2_B_L0DIS->Write();
        hxQ2_D_L0DIS->Write();
	hxQ2_B_Ph->Write();
        hxQ2_D_Ph->Write();
	hxQ2_B_Re->Write();
        hxQ2_D_Re->Write();
	hxQ2_D_Eta->Write();
	res->Write();
	dca->Write();
	dcaxy->Write();
	dcaz->Write();
	dcaz1->Write();
	dcaz2->Write();
	dcaz3->Write();
	etapt->Write();

	D0_DecayL->Write();
	D0_DCA->Write();
	D0_DDCA->Write();
	D0_D1DCA->Write();
	D0_D2DCA->Write();
	D0_CosTheta->Write();
	D0_Boost->Write();
	D0_DecayL_WS->Write();
        D0_DCA_WS->Write();
        D0_DDCA_WS->Write();
        D0_D1DCA_WS->Write();
        D0_D2DCA_WS->Write();
        D0_CosTheta_WS->Write();
	D0_Boost_WS->Write();
	_hxQ2->Write("hxQ2");
	_hEvents->Write("hEvents");
	hVtx_Mult_XY->Write();
	hVtx_Mult_Z->Write();
	hVtx_Mult_Y->Write();
	hVtx_Mult_X->Write();
	hVtx_Mult->Write();
	hVtx_Mult_Eta3->Write();
	hVtx_MultHF->Write();
        hVtx_MultHF_Eta3->Write();
	hVtx_MultPt->Write();
        hVtx_MultPt_Eta3->Write();
        hVtx_MultPtHF->Write();
        hVtx_MultPtHF_Eta3->Write();
    }
}
TLorentzVector * smearMomEIC(TLorentzVector const * b)
{
    float const pt = b->Perp();
    float const p = b->P();
    float sPt = -1;
    float sP = -1;
    float eta = b->PseudoRapidity();
    
    if(fabs(eta)<1)sP = gRandom->Gaus(p, p*sqrt(p*p*0.0002*0.0002+0.005*0.005));
    else if(eta < -2.5) sP = gRandom->Gaus(p, p*sqrt(p*p*0.001*0.001+0.02*0.02));//don't have lower limit; set limits on track cuts  
    else if(eta > -2.5 && eta < -2.0) sP = gRandom->Gaus(p, p*sqrt(p*p*0.0002*0.0002+0.01*0.01));
    else if(eta > -2.0 && eta < -1.0) sP = gRandom->Gaus(p, p*sqrt(p*p*0.0002*0.0002+0.01*0.01));
    else if(eta > 1 && eta < 2.5) sP = gRandom->Gaus(p, p*sqrt(p*p*0.0002*0.0002+0.01*0.01));
    else if(eta > 2.5) sP = gRandom->Gaus(p, p*sqrt(p*p*0.001*0.001+0.02*0.02));//dont have upper limit; set limits on track cuts  
    else sP = gRandom->Gaus(p, p*sqrt(p*p*0.001*0.001+0.005*0.005));//catch all not really needed
    sPt = sP*TMath::Sin(b->Theta());
    TLorentzVector *sMom = new TLorentzVector();
    sMom->SetXYZM(sPt * cos(b->Phi()), sPt * sin(b->Phi()), sPt * sinh(eta), b->M());
    return sMom;
}
TLorentzVector * smearMomEIC2(TLorentzVector const * b, TH2F *h, TH2F* h2,int id)
{
    float const pt = b->Perp();
    float const p = b->P();
    float  eta = fabs(b->PseudoRapidity());
    if(eta>=3)eta=2.999;
    float sPt = -1;
    float sP = -1;
    double val = 0;
    if(id==211)val = h->GetBinContent(h->GetXaxis()->FindBin(p),h->GetYaxis()->FindBin(fabs(eta)));
    else if (id==321)val = h2->GetBinContent(h2->GetXaxis()->FindBin(p),h2->GetYaxis()->FindBin(fabs(eta)));
    if(val==0 || val<0 || val>1){
	val=sqrt(p*p*0.0005*0.0005+0.005*0.005);
	//cout << "ERROR IN MOM. SMEARING VAL = " << id << " " << val <<  " " << p << " " << eta <<endl;
    }
    sP = gRandom->Gaus(p, p*val);
    sPt = sP/cosh(b->PseudoRapidity());
    TLorentzVector *sMom = new TLorentzVector();
    sMom->SetXYZM(sPt * cos(b->Phi()), sPt * sin(b->Phi()), sPt * sinh(b->PseudoRapidity()), b->M());
    return sMom;
}
bool isFlag(int id, int arr[]){
    for(int i = 0; i < fabs(sizeof(arr)) ; i++)
    {
        if(arr[i] == id)
            return 1;
    }
    return 0;
}
TVector3 smearPosATHENA(TLorentzVector const * b, TVector3 const& pos, TGraph * g_rphi,  TGraph * g_z){
    float const pt = b->Perp();
    float const p = b->P();
    float rand_xy = 0;
    float rand_z = 0;
    double eta =  TMath::ATanH(b->Pz()/b->P());//Doing it this way to suppress TVector3 Warnings; depends on any pre-cuts           
    // resolutions are in microns, need in mm
    rand_xy=gRandom->Gaus(0,g_rphi->Eval(p))/1000.;
    rand_z=gRandom->Gaus(0,g_z->Eval(p))/1000.;

    TVector3 newPos(pos);
    newPos.SetZ(0);
    TVector3 momPerp(-b->Vect().Y(), b->Vect().X(), 0.0);
    newPos -= momPerp.Unit() * rand_xy;

    return TVector3(newPos.X(), newPos.Y(), pos.Z() + rand_z);

}
TLorentzVector * smearMomATHENA(TLorentzVector const * b,TGraph * g)
{
    float const pt = b->Perp();
    float const p = b->P();
    float sPt = -1;
    float sP = -1;
    float eta = b->PseudoRapidity();
    // resolutions are in relative fraction delta p/p, need to smear with absolute delta p 
    sP = gRandom->Gaus(p,p*g->Eval(p));
    sPt = sP*TMath::Sin(b->Theta());
    TLorentzVector *sMom = new TLorentzVector();
    sMom->SetXYZM(sPt * cos(b->Phi()), sPt * sin(b->Phi()), sPt * sinh(eta), b->M());
    return sMom;
}


TVector3 smearPosEIC(TLorentzVector const * b, TVector3 const& pos){
    
    float const pt = b->Perp();
    float rand_xy = 0;
    double eta =  TMath::ATanH(b->Pz()/b->P());//Doing it this way to suppress TVector3 Warnings; depends on any pre-cuts    
    
    if(fabs(eta) <= 1)rand_xy=gRandom->Gaus(0, sqrt(5.*5.+30.*30./pt/pt))/1000.;
    else if(eta < -2.5)rand_xy=gRandom->Gaus(0, sqrt(15.*15.+60.*60./pt/pt))/1000.;//dont have lower limit; set limits on track cuts
    else if(eta > -2.5 && eta < -2.0)rand_xy=gRandom->Gaus(0, sqrt(15.*15.+60.*60./pt/pt))/1000.;
    else if(eta > -2.0 && eta < -1.0)rand_xy=gRandom->Gaus(0, sqrt(10.*10.+40.*40./pt/pt))/1000.;
    else if(eta > 1.0 && eta < 2.0)rand_xy=gRandom->Gaus(0, sqrt(10.*10.+40.*40./pt/pt))/1000.;
    else if(eta > 2.0 && eta < 2.5)rand_xy=gRandom->Gaus(0, sqrt(15.*15.+60.*60./pt/pt))/1000.;
    else if(eta > 2.5 && eta < 3.0)rand_xy=gRandom->Gaus(0, sqrt(15.*15.+60.*60./pt/pt))/1000.;
    else if(eta > 3.0)rand_xy=gRandom->Gaus(0, sqrt(60.*60.+30.*30./pt/pt))/1000.;//dont have upper limit; set limits on track cuts  
    else rand_xy=gRandom->Gaus(0, sqrt(15.*15.+60.*60./pt/pt))/1000.;//catch all not really needed
    float rand_z;
    if(fabs(eta)<=1)rand_z = gRandom->Gaus(0, sqrt(15.*15.+60.*60./pt/pt))/1000.;
    else rand_z = rand_xy/sin(b->Theta());
    TVector3 newPos(pos);
    newPos.SetZ(0);
    TVector3 momPerp(-b->Vect().Y(), b->Vect().X(), 0.0);
    newPos -= momPerp.Unit() * rand_xy;

    return TVector3(newPos.X(), newPos.Y(), pos.Z() + rand_z);
}
TVector3 smearPosEICLBNL(TLorentzVector const * b, TVector3 const& pos){
    float const pt = b->Perp();
    double eta =  TMath::ATanH(b->Pz()/b->P());//Doing it this way to suppress TVector3 Warnings   
    float rand_xy = 0;
    if(fabs(eta)<1)rand_xy=gRandom->Gaus(0, sqrt(5.*5.+20.*20./pt/pt))/1000.;
    else if(fabs(eta)>1 && fabs(eta)<2)rand_xy=gRandom->Gaus(0, sqrt(10.*10.+25.*25./pt/pt))/1000.;
    else if(fabs(eta)>2)rand_xy=gRandom->Gaus(0, sqrt(10.*10.+30.*30./pt/pt))/1000.;
    if(rand_xy==0)cout << "Problem here " << eta << endl;
    float rand_z;
    if(fabs(eta)<=1)rand_z = gRandom->Gaus(0, sqrt(5.*5.+20.*20./pt/pt))/1000.;
    else rand_z = rand_xy/sin(b->Theta());
    TVector3 newPos(pos);
    newPos.SetZ(0);
    TVector3 momPerp(-b->Vect().Y(), b->Vect().X(), 0.0);
    newPos -= momPerp.Unit() * rand_xy;

    return TVector3(newPos.X(), newPos.Y(), pos.Z() + rand_z);
}
TVector3 smearPosEICLANL(TLorentzVector const * b, TVector3 const& pos){
    float const pt = b->Perp();
    double eta =  TMath::ATanH(b->Pz()/b->P());//Doing it this way to suppress TVector3 Warnings                                                                                                                                                                      
    float rand_xy = 0;
    if(fabs(eta)<1)rand_xy=gRandom->Gaus(0, sqrt(25.*25./pt/pt))/1000.;
    else if(fabs(eta)>1 && fabs(eta)<2)rand_xy=gRandom->Gaus(0, sqrt(20.*20.+30.*30./pt/pt))/1000.;
    else if(fabs(eta)>2)rand_xy=gRandom->Gaus(0, sqrt(40.*40.+30.*30./pt/pt))/1000.;
    float rand_z;
    if(fabs(eta)<=1)rand_z = gRandom->Gaus(0, sqrt(5.*5.+20.*20./pt/pt))/1000.;
    else rand_z = rand_xy/sin(b->Theta());
    TVector3 newPos(pos);
    newPos.SetZ(0);
    TVector3 momPerp(-b->Vect().Y(), b->Vect().X(), 0.0);
    newPos -= momPerp.Unit() * rand_xy;

    return TVector3(newPos.X(), newPos.Y(), pos.Z() + rand_z);
}
TVector3 smearPosEIC2(TLorentzVector const * b, TVector3 const& pos,TH1F* DCARes_T_Eta1,TH1F* DCARes_T_Eta2,TH1F* DCARes_T_Eta3,TH1F* DCARes_Z_Eta1,TH1F* DCARes_Z_Eta2,TH1F* DCARes_Z_Eta3){
    float const pt = b->Perp();
    float rand_xy = 0;
    double eta =  TMath::ATanH(b->Pz()/b->P());//Doing it this way to suppress TVector3 Warnings        
    if(fabs(eta)<1)rand_xy= gRandom->Gaus(0, DCARes_T_Eta1->GetBinContent(DCARes_T_Eta1->FindBin(pt)))/1000.;
    else if(fabs(eta)>1 && fabs(eta)<2)rand_xy= gRandom->Gaus(0, DCARes_T_Eta2->GetBinContent(DCARes_T_Eta2->FindBin(pt)))/1000.;
    else if(fabs(eta)>2 && fabs(eta)<5)rand_xy= gRandom->Gaus(0, DCARes_T_Eta3->GetBinContent(DCARes_T_Eta3->FindBin(pt)))/1000.;
    
    float rand_z = 0;
    if(fabs(eta)<1)rand_z= gRandom->Gaus(0, DCARes_Z_Eta1->GetBinContent(DCARes_Z_Eta1->FindBin(pt)))/1000.;
    else if(fabs(eta)>1 && fabs(eta)<2)rand_z= gRandom->Gaus(0, DCARes_Z_Eta2->GetBinContent(DCARes_Z_Eta2->FindBin(pt)))/1000.;
    else if(fabs(eta)>2 && fabs(eta)<5)rand_z= gRandom->Gaus(0, DCARes_Z_Eta3->GetBinContent(DCARes_Z_Eta3->FindBin(pt)))/1000.;

    if(rand_z==0 || rand_z>200) rand_z = rand_xy/sin(b->Theta());
    TVector3 newPos(pos);
    newPos.SetZ(0);
    TVector3 momPerp(-b->Vect().Y(), b->Vect().X(), 0.0);
    newPos -= momPerp.Unit() * rand_xy;
    return TVector3(newPos.X(), newPos.Y(), pos.Z() + rand_z);
}
TVector3 smearPosEIC3(TLorentzVector const * b, TVector3 const& pos, TH1F* hp,TH1F *hm,TH1F* DCARes_Z_Eta1,TH1F* DCARes_Z_Eta2,TH1F* DCARes_Z_Eta3, int id){
    float const pt = b->Perp();

    float rand_xy = 0;
    if(id>0)rand_xy = hp->GetRandom()*10;
    else rand_xy = hm->GetRandom()*10;
    float rand_z;
    double eta =  TMath::ATanH(b->Pz()/b->P());//Doing it this way to suppress TVector3 Warnings                                                                                                                                                                          
    if(fabs(eta)<1)rand_z= gRandom->Gaus(0, DCARes_Z_Eta1->GetBinContent(DCARes_Z_Eta1->FindBin(pt)))/1000.;
    else if(fabs(eta)>1 && fabs(eta)<2)rand_z= gRandom->Gaus(0, DCARes_Z_Eta2->GetBinContent(DCARes_Z_Eta2->FindBin(pt)))/1000.;
    else if(fabs(eta)>2 && fabs(eta)<5)rand_z= gRandom->Gaus(0, DCARes_Z_Eta3->GetBinContent(DCARes_Z_Eta3->FindBin(pt)))/1000.;
    if(rand_z==0 || rand_z>200) rand_z = rand_xy/sin(b->Theta());

    TVector3 newPos(pos);
    newPos.SetZ(0);
    TVector3 momPerp(-b->Vect().Y(), b->Vect().X(), 0.0);
    newPos -= momPerp.Unit() * rand_xy;

    return TVector3(newPos.X(), newPos.Y(), pos.Z() + rand_z);
}


float dcaSigned(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex){
    TVector3 posDiff = pos - vertex;
    float sign = posDiff.x() * p.y() - posDiff.y() * p.x() > 0 ? +1 : -1;
    return sign * p.Cross(posDiff.Cross(p)).Unit().Dot(posDiff);
}
float dcaXY(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex,int mult){
    TVector3 newPos(pos);
    newPos.SetZ(0);

    TVector3 newP(p);
    newP.SetZ(0);

    TVector3 newVertex(vertex);
    newVertex.SetZ(0);

    TVector3 posDiff = newPos - newVertex;
    float sign = posDiff.x() * p.y() - posDiff.y() * p.x() > 0 ? +1 : -1;
    return sign * (newP.Cross(posDiff.Cross(newP)).Unit().Dot(posDiff));
}
float dcaZ(TVector3 const& p, TVector3 const& pos, TVector3 const& vertex){
    TVector3 posDiff = pos - vertex;
    if (sin(p.Theta()) == 0) return 0;
    else return (-posDiff.x() * cos(p.Phi()) - posDiff.y() * sin(p.Phi())) * cos(p.Theta()) / sin(p.Theta()) + posDiff.z();
}
bool isGoodD0(TLorentzVector *P,TVector3 p1, TVector3 pos1,TVector3 p2, TVector3 pos2,TVector3 vertex, int d1, int d2){
    TVector3 dca_vec = pos1-pos2;
    double dca_kpi = dca_vec.Mag();
    double dca_p1 = -999;
    double dca_p2 = -999;
    if(d1==211){
	dca_p1 = dcaSigned(p1,pos1,vertex);
	dca_p2 = dcaSigned(p2,pos2,vertex);
    }
    else{
	dca_p1 = dcaSigned(p2,pos2,vertex);
	dca_p2 = dcaSigned(p1,pos1,vertex);
    }
    TVector3 DecayL = (pos1+pos2)*0.5 - vertex;
    double D0_vtx = DecayL.Mag();
    double D0_dca = fabs(dcaSigned(P->Vect(), DecayL, vertex));
    TVector3 D0_vec = P->Vect();
    double D0_costheta = TMath::Cos(D0_vec.Angle(DecayL));
    double eta = P->PseudoRapidity();
    D0_DecayL->Fill(P->Pt(),D0_vtx,eta);
    D0_DCA->Fill(P->Pt(),D0_dca,eta);
    D0_DDCA->Fill(P->Pt(),dca_kpi,eta);
    D0_D1DCA->Fill(P->Pt(),dca_p1,eta);
    D0_D2DCA->Fill(P->Pt(),dca_p2,eta);
    D0_CosTheta->Fill(P->Pt(),D0_costheta,eta);

    if(D0_vtx>0.1 && fabs(dca_kpi)<0.07 && fabs(D0_dca)>0.01 && fabs(dca_p1)>0.03 && fabs(dca_p2)>0.03 && D0_costheta>0.95)
	return true;
    return false;

}
bool isGoodD0_XY(TLorentzVector *P,TVector3 p1, TVector3 pos1,TVector3 p2, TVector3 pos2,TVector3 vertex, int d1, int d2, TH2F *TrackEff_EtaP, int fill, double e1, double e2){
    


    StPhysicalHelixD Helix1(p1, pos1, 0.*kilogauss, 1);
    StPhysicalHelixD Helix2(p2, pos2, 0.*kilogauss, -1);

    double Path1 = (Helix1.pathLengths(Helix2)).first;
    double Path2 = (Helix1.pathLengths(Helix2)).second;

    TVector3 newPos1 = Helix1.at(Path1);
    TVector3 newPos2 = Helix2.at(Path2);
    TVector3 newMom1 = Helix1.momentumAt(Path1, 0. * kilogauss);
    TVector3 newMom2 = Helix2.momentumAt(Path2, 0. * kilogauss);

    TVector3 newpos1(newPos1.X(),newPos1.Y(),0);
    TVector3 newp1(newMom1.X(),newMom1.Y(),0);
    TVector3 newpos2(newPos2.X(),newPos2.Y(),0);
    TVector3 newp2(newMom2.X(),newMom2.Y(),0);
    TVector3 newvertex(vertex.X(),vertex.Y(),0);

    TVector3 dca_vec = newpos1-newpos2;
    double dca_kpi = dca_vec.Mag();
    double dca_p1 = -999;
    double dca_p2 = -999;
    if(d1==211){
	dca_p1 = fabs(dcaXY(p1,pos1,newvertex,-1));
	dca_p2 = fabs(dcaXY(p2,pos2,newvertex,-1));
    }
    else{
	dca_p1 = fabs(dcaXY(p2,pos2,newvertex,-1));
        dca_p2 = fabs(dcaXY(p1,pos1,newvertex,-1));
    }
    TVector3 DecayL = (newpos1+newpos2)*0.5 - newvertex;
    double D0_vtx = DecayL.Mag();
    double D0_dca = fabs(dcaXY(P->Vect(), DecayL, newvertex,-1));
    TVector3 D0_vec(P->X(),P->Y(),0);
    double D0_costheta = TMath::Cos(D0_vec.Angle(DecayL));
    double eta = P->PseudoRapidity();
    
//if(!passTracking(TrackEff_EtaP,p1.Perp(),fabs(p1.PseudoRapidity())) || !passTracking(TrackEff_EtaP,p2.Perp(),fabs(p2.PseudoRapidity())))return false;
    if(fabs(P->M()-1.870)<0.03){
	if(fill==1){
	    D0_DecayL->Fill(P->Pt(),D0_vtx,eta);
	    D0_DCA->Fill(P->Pt(),D0_dca,eta);
	    D0_DDCA->Fill(P->Pt(),dca_kpi,eta);
	    D0_D1DCA->Fill(P->Pt(),dca_p1,eta);
	    D0_D2DCA->Fill(P->Pt(),dca_p2,eta);
	    D0_CosTheta->Fill(P->Pt(),D0_costheta,eta);
	}
	if(fill==2){
	    D0_DecayL_WS->Fill(P->Pt(),D0_vtx,eta);
	    D0_DCA_WS->Fill(P->Pt(),D0_dca,eta);
	    D0_DDCA_WS->Fill(P->Pt(),dca_kpi,eta);
	    D0_D1DCA_WS->Fill(P->Pt(),dca_p1,eta);
	    D0_D2DCA_WS->Fill(P->Pt(),dca_p2,eta);
	    D0_CosTheta_WS->Fill(P->Pt(),D0_costheta,eta);
	} 
    }
    if(D0_vtx>0.04 && fabs(dca_kpi)<0.150 &&  D0_costheta>0.98)
        return true;
    return false;

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
int isFromD(int id){
    if(fabs(id) == 411 ||     //D+                                                                                                                                                            
       fabs(id) == 421 ||     //D0                                                                                                                                                       
       fabs(id) == 431 ||     //Ds                                                                                                                                                       
       fabs(id) == 4122 ||    //Lambda c                                                                                                                                              
       fabs(id) == 443)     //Jpsi                                                                                                                                                             
       	return 1;
    else return 0;
}
bool isKPiDecay(int id1, int id2){
    if( (id1 == 211 && id2 == -321) || //pi+ K-
	(id1 == -211 && id2 == 321) || //pi- K+
	(id1 == -321 && id2 == 211) || //pi+ K-   
	(id1 == 321 && id2 == -211)) //pi- K+  
	return true;
    return false;
}
bool isWSKPiDecay(int id1, int id2){
    if( (id1 == 211 && id2 == 321) || //pi+ K-                                                                                                                                                                                                                 
        (id1 == -211 && id2 == -321) || //pi- K+                                                                                                                                                                                                              
        (id1 == -321 && id2 == -211) || //pi+ K-                                                                                                                                                                                                              
        (id1 == 321 && id2 == 211)) //pi- K+                                                                                                                                                                                                                   
        return true;
    return false;
}
bool iseeDecay(int id1, int id2){
    if( (id1 == 11 && id2 == -11) || 
        (id1 == -11 && id2 == 11)) 
        return true;
    return false;
}
bool isK2PiDecay(int id1, int id2, int id3){
    if( (id1 == 211 && id2 == 211 && id3 == -321) || //pi+ pi+ K-                                                                                                       
	(id1 == -211 && id2 == -211 && id3 == 321) || //pi- pi- K+
	(id1 == 211 && id2 == -321 && id3 == 211) || //pi+ K- pi+                                                                                                                                 
	(id1 == -211 && id2 == 321 && id3 == -211) || //pi- K+ pi-   
	(id1 == -321 && id2 == 211 && id3 == 211) || //K- pi+ pi+                                                                                                                              
        (id1 == 321 && id2 == -211 && id3 == -211)) //K+ pi- pi- 
        return true;
    return false;


}
bool passTracking(TH2F*h,double p, double eta){
    if(eta>3)return false;//eta=3;
    int bin1 = h->GetXaxis()->FindBin(p);
    int bin2 = h->GetYaxis()->FindBin(eta);
    double val = h->GetBinContent(bin1,bin2);
    double ran = gRandom->Uniform(0,1);
    if(ran<val)return true;
    return false;
}

bool passTrackingATHENA(double p, double eta){
    if(fabs(eta)>3.5)return false;
    if(p<0.5)return false;
    return true;
}
int getEtaIdx(double eta){
    if(fabs(eta)<1)return 0;
    else if(fabs(eta)>1 && fabs(eta)<2)return 1;
    else if(fabs(eta)>2 && fabs(eta)<200)return 2;
}
int getPtIdx(double pt){
    double bin[38+1]={   0,0.2,0.4,0.6,0.8,
                         1,1.2,1.4,1.6,1.8,
                         2,2.2,2.4,2.6,2.8,
                         3,3.2,3.4,3.6,3.8,
                         4,4.2,4.4,4.6,4.8,
                         5,5.4,5.8,
                         6,6.4,6.8,
                         7.5,
                         8.5,
                         9.5,
                         11,13,15,
                         25,30000};
    for(int i = 0;i<37;i++)if(pt>bin[i] && pt<bin[i+1])return i;
    return 37;
    
}
bool passPID(int p1, int p2, double pt1, double pt2, double eta1, double eta2){
    
    double val1=-1;double val2=-1;
    
    if(fabs(eta1)<=1)val1=6;
    else if(eta1<-1) val1=10;
    else if(eta1>1) val1=50;
    else val1=6;//Not really needed

    if(fabs(eta2)<=1)val2=6;
    else if(eta2<-1) val2=10;
    else if(eta2>1)val2=50;
    else val2=6;//Not really needed
    
    if(pt1<=val1 && pt2<=val2){// both have PID 
	if(isKPiDecay(p1,p2) || isWSKPiDecay(p1,p2)){
	    return true;
	}else return false;
    }else if(pt1<=val1 &&pt2>val2){// only first can have PID
	if( (fabs(p1)==211 || fabs(p1)==321) ){
	    return true;
	} else return false;
    }else if (pt1>val1 &&pt2<=val2){// only second can have PID
	if( (fabs(p2)==211 || fabs(p2)==321) ){
            return true;
        } else return false;
    }else if (pt1>val1 && pt2>val2){// no PID
	return true;
    }
    else return false;
}
