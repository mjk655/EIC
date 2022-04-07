double const nbins_q2 = 13;
double const nbins_x  = 20;
double const nbins_y  = 4;
double binning_q2[nbins_q2+1];//={0.575,0.675,0.85,1,1.5,2.075,2.175,3};
double binning_x[nbins_x+1];//  ={-5,-4.75,-4.5,-4.25,-4,-3.75,-3.5,-3.25,-3,-2.75,-2.5,-2.25,-2,-1.75,-1.5,-1.25,-1,-0.75,-0.5,-0.25,0};
double binning_y[nbins_y+1]  ={0,0.25,0.5,0.75,1};//.,2./3.,1};
double BR =0.0395;
void setXbins(){
    for(int i = 0;i<nbins_x+1;i++){
	binning_x[i] = -4 + i * 4 / nbins_x;
	cout <<"x binning " << binning_x[i] << endl;
    }
}
void setQ2bins(){
    for(int i = 0;i<nbins_q2+1;i++){
        binning_q2[i] = i * 2.6 / nbins_q2;
    }
}
int markers[100]={20,21,22,23,24,25,26,27,28,29,30,20,21,23,24,25,26,27,28,29,30,21,22,23,24,25,26,27,28,29,30};
int colors[100]={1,14,33,4,9,38,8,30,40,42,46,6,2};
void Check(){
    gROOT->ProcessLine(".x ~/myStyle.C");
    setXbins();
    setQ2bins();
    TFile *f1 = new TFile("out_10_100_fulltree_nucl.root");//out_10_100_ATHENA.root","READ");
    TFile *f2 = new TFile("out_10_100_ATHENA.root");//out_10_100_ATHENA.root","READ");
    h3d1 = (TH3F*) f1->Get("D0_RedCS_Topo_T");h3d1->SetName("h3d1");h3d1->SetDirectory(0);
    //h3d_bkg1 = (TH3F*) f1->Get("D0_RedCS_BKG");h3d_bkg1->SetName("h3d_bkg1");h3d_bkg1->SetDirectory(0);
    h3d_bkg1 = (TH3F*) f2->Get("D0_RedCS_Topo_T");h3d_bkg1->SetName("h3d_bkg1");h3d_bkg1->SetDirectory(0);
    //h3d1->GetXaxis()->SetRangeUser(binning_x[14],binning_x[15]);
    //h3d1->GetYaxis()->SetRangeUser(binning_q2[8],binning_q2[9]);
    //h3d_bkg1->GetXaxis()->SetRangeUser(binning_x[14],binning_x[15]);
    //h3d_bkg1->GetYaxis()->SetRangeUser(binning_q2[8],binning_q2[9]);
    h1 = (TH1F*)h3d1->Project3D("z");
    h2 = (TH1F*)h3d_bkg1->Project3D("z");
    h1->Draw("hist");
    h2->Draw("PE same");
    cout << h1->Integral() << " " << h2->Integral() << endl;
}
void Calc(int save=1){
    char outname[100]="ep_output_ATHENA1.root";
    char system[10]="ep";
    gROOT->ProcessLine(".x ~/myStyle.C");
    gStyle->SetPalette(51);
    double scale2[nbins_q2]={0.015,0.025,0.025,0.03,0.05,0.1,0.2,0.4,0.6,1.5,4,2,3};
    double scale[nbins_q2]={0.015,0.025,0.025,0.03,0.05,0.1,0.2,0.4,0.6,1.5,4,2,3};
    char scenario[100]=" ";//Perfect PID + Vertexing";

    double CS = 1;
    double CS1 = 1;
    double CS2 = 1;
    double CS11 = 1;
    double CS22 = 1;
    double CS1_10 = 1;
    double CS2_10 = 1;
    double CS11_10 = 1;
    double CS22_10 = 1;
//============================ 18x275 GeV
    TFile *f = new TFile("./out_18_275.root","READ");
    h3d = (TH3F*) f->Get("D0_RedCS_Topo_T");h3d->SetName("h3d");h3d->SetDirectory(0);
    h3d_bkg = (TH3F*) f->Get("D0_RedCS_Topo_BKG");h3d_bkg->SetName("h3d_bkg");h3d_bkg->SetDirectory(0);
    h3d_d = (TH3F*) f->Get("D0_RedCS_T_Acc");h3d_d->SetName("h3d_d");h3d_d->SetDirectory(0);
    TH2F* _hxQ2 = (TH2F*) f->Get("hxQ2");_hxQ2->SetName("_hxQ2");
    CS = getXsec(_hxQ2,1./0.94);
    cout << "Calculated Lumi in 18x275 is " << CS << endl;
//============================
    TFile *f1 = new TFile("./out_10_100_ATHENA.root","READ");
    //TFile *f1 = new TFile("out_10_100_fulltree_nucl.root","READ");
    h3d1 = (TH3F*) f1->Get("D0_RedCS_Topo_T");h3d1->SetName("h3d1");h3d1->SetDirectory(0);
    h3d_bkg1 = (TH3F*) f1->Get("D0_RedCS_Topo_BKG");h3d_bkg1->SetName("h3d_bkg1");h3d_bkg1->SetDirectory(0);
    h3d_d1 = (TH3F*) f1->Get("D0_RedCS_T_Acc");h3d_d1->SetName("h3d_d1");h3d_d1->SetDirectory(0);
    TH2F* _hxQ21 = (TH2F*) f1->Get("hxQ2");_hxQ21->SetName("_hxQ21");
    CS1 = getXsec(_hxQ21,1/0.65);
    //CS1 = getXsec(_hxQ21,1/0.56);
    cout << "Calculated Lumi in 10x100 is " << CS1 << endl;
//============================  
    TFile *f2 = new TFile("./out_5_41_ATHENA.root","READ");
    //TFile *f2 = new TFile("./out_5_41_fulltree_nucl.root","READ");
    h3d2 = (TH3F*) f2->Get("D0_RedCS_Topo_T");h3d2->SetName("h3d2");h3d2->SetDirectory(0);
    h3d_bkg2 = (TH3F*) f2->Get("D0_RedCS_Topo_BKG");h3d_bkg2->SetName("h3d_bkg2");h3d_bkg2->SetDirectory(0);
    h3d_d2 = (TH3F*) f2->Get("D0_RedCS_T_Acc");h3d_d2->SetName("h3d_d2");h3d_d2->SetDirectory(0);
    TH2F* _hxQ22 = (TH2F*) f2->Get("hxQ2");_hxQ22->SetName("_hxQ22");
    CS2 = getXsec(_hxQ22,1/0.42);
    cout << "Calculated Lumi in 5x41 is " << CS2 << endl;
//============================                      
    TFile *f1_10 = new TFile("./out_10_100_Q210_ATHENA.root","READ");
    //TFile *f1_10 = new TFile("out_10_100_fulltree_nucl_Q210.root","READ");
    h3d1_10 = (TH3F*) f1_10->Get("D0_RedCS_Topo_T");h3d1_10->SetName("h3d1_10");h3d1_10->SetDirectory(0);
    h3d_bkg1_10 = (TH3F*) f1_10->Get("D0_RedCS_Topo_BKG");h3d_bkg1_10->SetName("h3d_bkg1_10");h3d_bkg1_10->SetDirectory(0);
    h3d_d1_10 = (TH3F*) f1_10->Get("D0_RedCS_T_Acc");h3d_d1_10->SetName("h3d_d1_10");h3d_d1_10->SetDirectory(0);
    TH2F* _hxQ21_10 = (TH2F*) f1_10->Get("hxQ2");_hxQ21_10->SetName("_hxQ21_10");
    CS1_10 = getXsec(_hxQ21_10,1/0.042);
    //CS1_10 = getXsec(_hxQ21_10,1/0.56);
    cout << "Calculated Lumi in 10x100 q2=10 is " << CS1_10 << endl;
//============================                                            
    TFile *f2_10 = new TFile("./out_5_41_Q210_ATHENA.root","READ");
    //TFile *f2_10 = new TFile("./out_5_41_fulltree_nucl_Q210.root","READ");
    h3d2_10 = (TH3F*) f2_10->Get("D0_RedCS_Topo_T");h3d2_10->SetName("h3d2_10");h3d2_10->SetDirectory(0);
    h3d_bkg2_10 = (TH3F*) f2_10->Get("D0_RedCS_Topo_BKG");h3d_bkg2_10->SetName("h3d_bkg2_10");h3d_bkg2_10->SetDirectory(0);
    h3d_d2_10 = (TH3F*) f2_10->Get("D0_RedCS_T_Acc");h3d_d2_10->SetName("h3d_d2_10");h3d_d2_10->SetDirectory(0);
    TH2F* _hxQ22_10 = (TH2F*) f2_10->Get("hxQ2");_hxQ22_10->SetName("_hxQ22_10");
    CS2_10 = getXsec(_hxQ22_10,1/0.0204);
    cout << "Calculated Lumi in 5x41 q2=10 is " << CS2_10 << endl;
//============================                                                                                                                                                                                                                             
    TFile *f11 = new TFile("./out_10_100_PID.root","READ");
    h3d11 = (TH3F*) f11->Get("D0_RedCS_Topo_T");h3d11->SetName("h3d1");h3d11->SetDirectory(0);
    h3d_bkg11 = (TH3F*) f11->Get("D0_RedCS_Topo_BKG");h3d_bkg11->SetName("h3d_bkg1");h3d_bkg11->SetDirectory(0);
    h3d_d11 = (TH3F*) f11->Get("D0_RedCS_T_Acc");h3d_d11->SetName("h3d_d1");h3d_d11->SetDirectory(0);
    TH2F* _hxQ211 = (TH2F*) f11->Get("hxQ2");_hxQ211->SetName("_hxQ211");
    CS11 = getXsec(_hxQ211,1/0.65);
    cout << "Calculated Lumi in 10x100 is " << CS11 << endl;
//============================                                                                                                                                                                                                                        
    TFile *f22 = new TFile("./out_5_41_PID.root","READ");
    h3d22 = (TH3F*) f22->Get("D0_RedCS_Topo_T");h3d22->SetName("h3d22");h3d22->SetDirectory(0);
    h3d_bkg22 = (TH3F*) f22->Get("D0_RedCS_Topo_BKG");h3d_bkg22->SetName("h3d_bkg22");h3d_bkg22->SetDirectory(0);
    h3d_d22 = (TH3F*) f22->Get("D0_RedCS_T_Acc");h3d_d22->SetName("h3d_d22");h3d_d22->SetDirectory(0);
    TH2F* _hxQ222 = (TH2F*) f22->Get("hxQ2");_hxQ222->SetName("_hxQ222");
    CS22 = getXsec(_hxQ222,1/0.42);
    cout << "Calculated Lumi in 5x41 is " << CS22 << endl;
//============================  
    TFile *f11_10 = new TFile("./out_10_100_Q210_PID.root","READ");
    h3d11_10 = (TH3F*) f11_10->Get("D0_RedCS_Topo_T");h3d11_10->SetName("h3d1_10");h3d11_10->SetDirectory(0);
    h3d_bkg11_10 = (TH3F*) f11_10->Get("D0_RedCS_Topo_BKG");h3d_bkg11_10->SetName("h3d_bkg1_10");h3d_bkg11_10->SetDirectory(0);
    h3d_d11_10 = (TH3F*) f11_10->Get("D0_RedCS_T_Acc");h3d_d11_10->SetName("h3d_d1_10");h3d_d11_10->SetDirectory(0);
    TH2F* _hxQ211_10 = (TH2F*) f11_10->Get("hxQ2");_hxQ211_10->SetName("_hxQ211_10");
    CS11_10 = getXsec(_hxQ211_10,1/0.042);
    cout << "Calculated Lumi in 10x100 q2=10 is " << CS11_10 << endl;
//============================ 
    TFile *f22_10 = new TFile("./out_5_41_Q210_PID.root","READ");
    h3d22_10 = (TH3F*) f22_10->Get("D0_RedCS_Topo_T");h3d22_10->SetName("h3d22_10");h3d22_10->SetDirectory(0);
    h3d_bkg22_10 = (TH3F*) f22_10->Get("D0_RedCS_Topo_BKG");h3d_bkg22_10->SetName("h3d_bkg22_10");h3d_bkg22_10->SetDirectory(0);
    h3d_d22_10 = (TH3F*) f22_10->Get("D0_RedCS_T_Acc");h3d_d22_10->SetName("h3d_d22_10");h3d_d22_10->SetDirectory(0);
    TH2F* _hxQ222_10 = (TH2F*) f22_10->Get("hxQ2");_hxQ222_10->SetName("_hxQ222_10");
    CS22_10 = getXsec(_hxQ222_10,1/0.0204);
    cout << "Calculated Lumi in 5x41 q2=10 is " << CS22_10 << endl;
//============================ 
    
    setXbins();
    setQ2bins();
    TH1F *hCS_x[nbins_q2];
    TH1F *hdCS_x[nbins_q2];
    TH1F *hmean_x[nbins_q2];
    TH1F *hmean_q2[nbins_q2];
    
    TH1F *hCS_x1[nbins_q2];
    TH1F *hdCS_x1[nbins_q2];
    TH1F *hmean_x1[nbins_q2];
    TH1F *hmean_q21[nbins_q2];

    TH1F *hCS_x2[nbins_q2];
    TH1F *hdCS_x2[nbins_q2];
    TH1F *hmean_x2[nbins_q2];
    TH1F *hmean_q22[nbins_q2];

    TH1F *hCS_x11[nbins_q2];
    TH1F *hdCS_x11[nbins_q2];
    TH1F *hmean_x11[nbins_q2];
    TH1F *hmean_q211[nbins_q2];

    TH1F *hCS_x22[nbins_q2];
    TH1F *hdCS_x22[nbins_q2];
    TH1F *hmean_x22[nbins_q2];
    TH1F *hmean_q222[nbins_q2];
    
    TH1F *hF2_x[nbins_q2];
    TH1F *hdF2_x[nbins_q2];
    TH1F *hsdF2_x[nbins_q2];
    TH1F *hs1dF2_x[nbins_q2];
    TH1F *hF2_x2[nbins_q2];
    TH1F *hdF2_x2[nbins_q2];
    TH1F *hsdF2_x2[nbins_q2];
    TH1F *hs1dF2_x2[nbins_q2];
    
    for(int i = 0;i < nbins_q2; i++){
	char name[100];
	sprintf(name,"hCS_x_%i",i);
	hCS_x[i] = new TH1F(name,name,nbins_x,binning_x);
        sprintf(name,"hmean_x_%i",i);
        hmean_x[i] = new TH1F(name,name,nbins_x,binning_x);
	sprintf(name,"hmean_q2_%i",i);
        hmean_q2[i] = new TH1F(name,name,nbins_x,binning_x);
	sprintf(name,"hdCS_x_%i",i);
        hdCS_x[i] = new TH1F(name,name,nbins_x,binning_x);
	
        sprintf(name,"hCS_x1_%i",i);
        hCS_x1[i] = new TH1F(name,name,nbins_x,binning_x);
        sprintf(name,"hmean_x1_%i",i);
        hmean_x1[i] = new TH1F(name,name,nbins_x,binning_x);
        sprintf(name,"hmean_q21_%i",i);
        hmean_q21[i] = new TH1F(name,name,nbins_x,binning_x);
        sprintf(name,"hdCS_x1_%i",i);
        hdCS_x1[i] = new TH1F(name,name,nbins_x,binning_x);

	sprintf(name,"hCS_x2_%i",i);
        hCS_x2[i] = new TH1F(name,name,nbins_x,binning_x);
        sprintf(name,"hmean_x2_%i",i);
        hmean_x2[i] = new TH1F(name,name,nbins_x,binning_x);
        sprintf(name,"hmean_q22_%i",i);
        hmean_q22[i] = new TH1F(name,name,nbins_x,binning_x);
        sprintf(name,"hdCS_x2_%i",i);
        hdCS_x2[i] = new TH1F(name,name,nbins_x,binning_x);


	sprintf(name,"hCS_x11_%i",i);
        hCS_x11[i] = new TH1F(name,name,nbins_x,binning_x);
        sprintf(name,"hmean_x11_%i",i);
        hmean_x11[i] = new TH1F(name,name,nbins_x,binning_x);
        sprintf(name,"hmean_q211_%i",i);
        hmean_q211[i] = new TH1F(name,name,nbins_x,binning_x);
        sprintf(name,"hdCS_x11_%i",i);
        hdCS_x11[i] = new TH1F(name,name,nbins_x,binning_x);

        sprintf(name,"hCS_x22_%i",i);
        hCS_x22[i] = new TH1F(name,name,nbins_x,binning_x);
        sprintf(name,"hmean_x22_%i",i);
        hmean_x22[i] = new TH1F(name,name,nbins_x,binning_x);
        sprintf(name,"hmean_q222_%i",i);
        hmean_q222[i] = new TH1F(name,name,nbins_x,binning_x);
        sprintf(name,"hdCS_x22_%i",i);
        hdCS_x22[i] = new TH1F(name,name,nbins_x,binning_x);

	sprintf(name,"hF2_x_%i",i);
        hF2_x[i] = new TH1F(name,name,nbins_x,binning_x);
        sprintf(name,"hdF2_x_%i",i);
	hdF2_x[i] = new TH1F(name,name,nbins_x,binning_x);
	sprintf(name,"hsdF2_x_%i",i);
        hsdF2_x[i] = new TH1F(name,name,nbins_x,binning_x);
	sprintf(name,"hs1dF2_x_%i",i);
        hs1dF2_x[i] = new TH1F(name,name,nbins_x,binning_x);
	sprintf(name,"hF2_x2_%i",i);
        hF2_x2[i] = new TH1F(name,name,nbins_x,binning_x);
        sprintf(name,"hdF2_x2_%i",i);
        hdF2_x2[i] = new TH1F(name,name,nbins_x,binning_x);
	sprintf(name,"hsdF2_x2_%i",i);
        hsdF2_x2[i] = new TH1F(name,name,nbins_x,binning_x);
	sprintf(name,"hs1dF2_x2_%i",i);
        hs1dF2_x2[i] = new TH1F(name,name,nbins_x,binning_x);
    }
    double qs[20];
    double qs1[20];
    double qs2[20];
    double qs11[20];
    double qs22[20];
    for(int i = 0;i<nbins_q2;i++){
	qs[i]  = getRCSVsX(hmean_q2[i],hmean_x[i],hCS_x[i],hdCS_x[i],h3d,h3d_d,h3d_bkg,i+1,CS,1,10);
	if(i<5){
	    qs1[i] = getRCSVsX(hmean_q21[i],hmean_x1[i],hCS_x1[i],hdCS_x1[i],h3d1,h3d_d1,h3d_bkg1,i+1,CS1,2,10);
	    qs2[i] = getRCSVsX(hmean_q22[i],hmean_x2[i],hCS_x2[i],hdCS_x2[i],h3d2,h3d_d2,h3d_bkg2,i+1,CS2,3,10);
	    qs11[i] = getRCSVsX(hmean_q211[i],hmean_x11[i],hCS_x11[i],hdCS_x11[i],h3d11,h3d_d11,h3d_bkg11,i+1,CS11,2,10);
	    qs22[i] = getRCSVsX(hmean_q222[i],hmean_x22[i],hCS_x22[i],hdCS_x22[i],h3d22,h3d_d22,h3d_bkg22,i+1,CS22,3,10);
	}
	else{
	    qs1[i] = getRCSVsX(hmean_q21[i],hmean_x1[i],hCS_x1[i],hdCS_x1[i],h3d1_10,h3d_d1_10,h3d_bkg1_10,i+1,CS1_10,2,10);
            qs2[i] = getRCSVsX(hmean_q22[i],hmean_x2[i],hCS_x2[i],hdCS_x2[i],h3d2_10,h3d_d2_10,h3d_bkg2_10,i+1,CS2_10,3,10);
            qs11[i] = getRCSVsX(hmean_q211[i],hmean_x11[i],hCS_x11[i],hdCS_x11[i],h3d11_10,h3d_d11_10,h3d_bkg11_10,i+1,CS11_10,2,10);
            qs22[i] = getRCSVsX(hmean_q222[i],hmean_x22[i],hCS_x22[i],hdCS_x22[i],h3d22_10,h3d_d22_10,h3d_bkg22_10,i+1,CS22_10,3,10);
	}
    }
    TCanvas *cfit1 = new TCanvas("cfit1","cfit1");//,1800,1200);
    TCanvas *cfit2 = new TCanvas("cfit2","cfit2");//,1800,1200);                                                                                                                                                                                                
    //  cfit1->Divide(5,5);
    double qs_1[20];
    double qs_2[20];
    for(int i = 0;i<nbins_q2;i++){
        qs_1[i]  = getF2VsX(hF2_x[i],hdF2_x[i],hsdF2_x[i],hs1dF2_x[i],hCS_x[i],hmean_x[i],hmean_q2[i],hCS_x1[i],hmean_x1[i],hmean_q21[i],hCS_x2[i],hmean_x2[i],hmean_q22[i],i+1,cfit1);
	qs_2[i]  = getF2VsX(hF2_x2[i],hdF2_x2[i],hsdF2_x2[i],hs1dF2_x2[i],hCS_x[i],hmean_x[i],hmean_q2[i],hCS_x11[i],hmean_x11[i],hmean_q211[i],hCS_x22[i],hmean_x22[i],hmean_q222[i],i+1,cfit2);
    }
    TLegend *leg = new TLegend(0.17,0.62,0.4,0.87);
    TLegend *leg_2 = new TLegend(0.75,0.7,0.87,0.94);
    TLegend *legF2 = new TLegend(0.15,0.7,0.35,0.94);
    TLegend *legF22 = new TLegend(0.35,0.7,0.55,0.94);
    legF2->SetTextSize(0.02);
    legF22->SetTextSize(0.02);
    leg->SetTextSize(0.02);
    leg_2->SetTextSize(0.02);
    TLegend *leg1 = new TLegend(0.18,0.22,0.33,0.85);
    TLegend *leg2 = new TLegend(0.18,0.4,0.3,0.84);
    for(int i = 0;i<nbins_q2-2;i++){
	char label[100];
	sprintf(label,"Q^{2} = %1.1f GeV^{2}",qs[i]);
	//leg->AddEntry(hCS_x[i],label,"PE");
	leg1->AddEntry(hCS_x[i],label,"l");
	leg2->AddEntry(hCS_x[i],label,"l");
	sprintf(label,"Q^{2} = %1.1f GeV^{2}, C = %1.3f ",qs[i],scale[i]);
	if(i<6){
	    sprintf(label,"Q^{2} = %1.1f GeV^{2}, C = %1.3f ",qs[i],scale[i]);
	    leg->AddEntry(hCS_x[i],label,"PE");
	}        
	else {
	    sprintf(label,"Q^{2} = %1.1f GeV^{2}, C = %1.1f ",qs[i],scale[i]);
	    leg_2->AddEntry(hCS_x[i],label,"PE");
	}
	if(i<6){
	    sprintf(label,"Q^{2} = %1.1f GeV^{2}, C = %1.3f ",qs[i],scale[i]);
	    legF2->AddEntry(hCS_x[i],label,"PE");
	}
	else {
	    sprintf(label,"Q^{2} = %1.1f GeV^{2}, C = %1.1f ",qs[i],scale[i]);
	    legF22->AddEntry(hCS_x[i],label,"PE");
	}
    }







    if(save){

	TFile *outfile = new TFile(outname,"RECREATE");
	char nname[100];
	for(int i = 0;i<nbins_q2;i++){
	    sprintf(nname,"%s_10100_Q2Bin_%i",system,i);
	    hCS_x1[i]->Write(nname);
	    sprintf(nname,"%s_541_Q2Bin_%i",system,i);
	    hCS_x2[i]->Write(nname);
	    sprintf(nname,"%s_CharmF2_Q2Bin_%i",system,i);
	    hF2_x[i]->Write(nname);
	}
    }
	


    
    TLatex lat;

/*    TCanvas *cx = new TCanvas("cx","cx");
    hCS_x[0]->GetYaxis()->SetRangeUser(0,hCS_x[0]->GetMaximum()*5);
    for(int i = 0;i<nbins_q2-2;i++)hCS_x[i]->Draw("PE X0 same");
    leg->Draw("same");
    lat.DrawLatex(-3.8,hCS_x[0]->GetMaximum()*0.9,"e+p 18+275 GeV 10 fb^{-1}");
    lat.DrawLatex(-3.8,hCS_x[0]->GetMaximum()*0.8,scenario);

    TCanvas *cdx = new TCanvas("cdx","cdx");
    hdCS_x[0]->GetYaxis()->SetRangeUser(0,hdCS_x[0]->GetMaximum()*2);
    for(int i = 0;i<nbins_q2-2;i++)hdCS_x[i]->Draw("hist same");
    leg1->Draw("same");
    lat.DrawLatex(-3.8,hdCS_x[0]->GetMaximum()*0.9,"e+p 18+275 GeV 10 fb^{-1}");
    lat.DrawLatex(-3.8,hdCS_x[0]->GetMaximum()*0.8,scenario);
*/
//==================================================== CS  ==========================       
    TGraphErrors *gCS1 = new TGraphErrors();
    TGraphErrors *gCS2 = new TGraphErrors();
    TGraphErrors *gF2 = new TGraphErrors();
//    TCanvas *cx1 = new TCanvas("cx1","cx1"
    hCS_x1[0]->GetXaxis()->SetRangeUser(-3.3,0);//hF2_x[0]->GetMaximum()*5);                                                                                                                                                                          
    hCS_x1[0]->GetYaxis()->SetRangeUser(0.000001,5);//hF2_x[0]->GetMaximum()*5);          
    for(int i = 0;i<nbins_q2-2;i++){
	hCS_x1[i]->Scale(scale[i]);
	hCS_x1[i]->GetYaxis()->SetRangeUser(0.000001,5);
//	hCS_x1[i]->Draw("PE X0 same");
    }
//    legF22->Draw("same");
//    legF2->Draw("same");
//    lat.DrawLatex(-1.6,1,"e+p 10+100 GeV 10 fb^{-1}");
//    lat.DrawLatex(-1.5,0.15,scenario);
//    gPad->SetLogy(1);

/*    TCanvas *cdx1 = new TCanvas("cdx1","cdx1");
    hdCS_x1[0]->GetYaxis()->SetRangeUser(0,0.7);//hdCS_x1[0]->GetMaximum()*2);
    for(int i = 0;i<nbins_q2-2;i++)hdCS_x1[i]->Draw("hist same");
    leg1->Draw("same");
    lat.DrawLatex(-3.8,hdCS_x1[0]->GetMaximum()*0.9,"e+p 10+100 GeV 10 fb^{-1}");
    lat.DrawLatex(-3.8,hdCS_x1[0]->GetMaximum()*0.8,scenario);
*/
//    TCanvas *cx2 = new TCanvas("cx2","cx2");
    hCS_x2[0]->GetXaxis()->SetRangeUser(-3.3,0);//hF2_x[0]->GetMaximum()*5);                                                                                                                                                                         
    hCS_x2[0]->GetYaxis()->SetRangeUser(0.000001,5);//hF2_x[0]->GetMaximum()*5);        
    for(int i = 0;i<nbins_q2-2;i++){
	hCS_x2[i]->Scale(scale[i]);
        hCS_x2[i]->GetYaxis()->SetRangeUser(0.000001,5);
//	hCS_x2[i]->Draw("PE X0 same");
    }
//    legF22->Draw("same");
//    legF2->Draw("same");
//    lat.DrawLatex(-1.6,1,"e+p 5+41 GeV 10 fb^{-1}");
//    lat.DrawLatex(-1.5,0.15,scenario);
//    gPad->SetLogy(1);

/*    TCanvas *cdx2 = new TCanvas("cdx2","cdx2");
    hdCS_x2[0]->GetYaxis()->SetRangeUser(0,0.7);//hdCS_x2[0]->GetMaximum()*2);
    for(int i = 0;i<nbins_q2-2;i++)hdCS_x2[i]->Draw("hist same");
    leg1->Draw("same");
    lat.DrawLatex(-3.8,hdCS_x2[0]->GetMaximum()*0.9,"e+p 5+41 GeV 10 fb^{-1}");
    lat.DrawLatex(-3.8,hdCS_x2[0]->GetMaximum()*0.8,scenario);
*/
//=================== All in one TGraph =====
    int cnt=1;
    for(int i = 0;i<nbins_q2-2;i++){
	for(int j = 1; j<hCS_x1[i]->GetNbinsX()+1;j++){
	    double err = hCS_x1[i]->GetBinError(j);
	    double x = hCS_x1[i]->GetBinCenter(j);
	    double y = hCS_x1[i]->GetBinContent(j);
	    //if(i==5 && (j > hCS_x1[i]->GetNbinsX()-3))continue;
	    //if(i==6 && (j > hCS_x1[i]->GetNbinsX()-2))continue;
	    //if(i==7 && (j > hCS_x1[i]->GetNbinsX()-2))continue;
	    if(y>0){
		gCS1->SetPoint(cnt,x,y);
		gCS1->SetPointError(cnt,0,err);
		cnt++;
	    }
	}
    }
    cnt=1;
    for(int i = 0;i<nbins_q2-2;i++){
	for(int j = 1; j<hCS_x2[i]->GetNbinsX()+1;j++){
            double err = hCS_x2[i]->GetBinError(j);
	    double x = hCS_x2[i]->GetBinCenter(j);
            double y = hCS_x2[i]->GetBinContent(j);
            
	    if(y>0){
                gCS2->SetPoint(cnt,x+0.04,y*1.2);
		gCS2->SetPointError(cnt,0,err);
		cnt++;
            }
	}
    }
    gCS2->SetMarkerStyle(25);gCS2->SetMarkerColor(kRed);gCS2->SetLineColor(kRed);
    //gCS1->SetLineStyle(7);
    TLatex lat;
    TCanvas *c12 = new TCanvas("c12","c12");
    gCS1->GetXaxis()->SetTitle("log(x_{B})");
    gCS1->GetYaxis()->SetTitle("#sigma_{r}^{c#bar{c}}(x_{B},Q^{2}) #times C");
    gCS1->GetYaxis()->SetLimits(0.00000001,5);
    gCS1->GetYaxis()->SetRangeUser(0.000005,0.4);
    gCS1->GetXaxis()->SetLimits(-6.9,0);
    gCS1->GetXaxis()->SetRangeUser(-3.9,0);
    gCS1->Draw("APE");
    gCS2->Draw("same PE");
    
    for(int i = 0;i<nbins_q2-2;i++){
	hCS_x1[i]->SetLineWidth(0.2);
	hCS_x2[i]->SetLineWidth(0.2);
	hCS_x1[i]->SetLineStyle(7);
	hCS_x1[i]->SetLineColor(1);
	hCS_x2[i]->SetLineColor(kRed);
	int lowbin = 0;hCS_x1[i]->GetBinWithContent(hCS_x1[i]->GetMinimum(),lowbin,0,nbins_q2,0);double low = hCS_x1[i]->GetBinLowEdge(hCS_x1[i]->GetMinimumBin());
	double high = hCS_x1[i]->GetBinLowEdge(hCS_x1[i]->GetMaximumBin())+hCS_x1[i]->GetBinWidth(hCS_x1[i]->GetMaximumBin());
	cout << low << " " << high << endl;
	hCS_x1[i]->GetXaxis()->SetLimits(high,low);
	hCS_x1[i]->GetXaxis()->SetRangeUser(high,low);
	//hCS_x1[i]->Draw("E5 same");    
	hCS_x2[i]->SetLineStyle(7);
	hCS_x2[i]->Scale(1.2);
//	hCS_x2[i]->Draw("E5 same");
    }
    gPad->SetLogy(1);
    TLegend *leg_cs = new TLegend(0.17,0.2,0.4,0.4);
    leg_cs->SetTextSize(0.045);
    leg_cs->SetHeader("PYTHIA6 e+p");
    leg_cs->AddEntry(gCS1,"10#times100 GeV 10 fb^{-1}","P");
    leg_cs->AddEntry(gCS2,"5#times41 GeV 10 fb^{-1} (x1.2)","P");
    lat.SetTextSize(0.04);
    double latx[nbins_q2-2]={-3.6,-3.5,-3.4,-3.0,-2.9,-2.46,-2.07,-2.0,-2.0,-1.8,-1.6};
    double laty[nbins_q2-2]={0.0002,0.00052,0.0013,0.003,0.006,0.012,0.023,0.038,0.065,0.1,0.18};
    double lata[nbins_q2-2]={0,0,0,-3,-5,-10,-10,-10,-10,-15,-15};
    //double latx[nbins_q2-2]={-4.2,-4.1,-4,-3.7,-3.5,-3.26,-3.07,-2.83,-2.6,-2.3,-2};
    //double laty[nbins_q2-2]={0.0002,0.00052,0.0013,0.003,0.006,0.012,0.023,0.038,0.065,0.13,0.23};
    //double lata[nbins_q2-2]={0,0,0,-3,-5,-10,-10,-10,-10,-15,-15};
    for(int i = 0;i<nbins_q2-2;i++){
//	if(i>3)continue;
	cout << "Drawing " << i << endl;
	lat.SetTextAngle(lata[i]);
	//if(i<3)lat.DrawLatex(latx[i],laty[i],Form("#color[4]{Q^{2}=%1.1f GeV^{2}, C=%1.3f}",qs[i],scale[i]));
	//else if(i<5) lat.DrawLatex(latx[i],laty[i],Form("#color[4]{Q^{2}=%1.f GeV^{2}, C=%1.2f}",qs[i],scale[i]));
	//else lat.DrawLatex(latx[i],laty[i],Form("#color[4]{Q^{2}=%1.f GeV^{2}, C=%1.1f}",qs[i],scale[i]));
	if(i<3)lat.DrawLatex(latx[i],laty[i],Form("#color[4]{%1.1f, %1.3f}",qs[i],scale[i]));
        else if(i<5) lat.DrawLatex(latx[i],laty[i],Form("#color[4]{%1.f, %1.2f}",qs[i],scale[i]));
        else lat.DrawLatex(latx[i],laty[i],Form("#color[4]{%1.f, %1.1f}",qs[i],scale[i]));
    }
    leg_cs->Draw("same");
    lat.SetTextAngle(0);
    lat.DrawLatex(-0.8,0.1,"#color[4]{Q^{2} (GeV^{2}), C}");

//    TCanvas *cf2 = new TCanvas("cf2","cf2");
    hF2_x[0]->GetXaxis()->SetRangeUser(-3.,0);//hF2_x[0]->GetMaximum()*5);            
    hF2_x[0]->GetYaxis()->SetRangeUser(0.0001,0.1);//hF2_x[0]->GetMaximum()*5);
    for(int i = 0;i<nbins_q2-2;i++){
	hF2_x[i]->Scale(scale[i]);
	hF2_x[i]->GetYaxis()->SetRangeUser(0.0002,0.1);
	//hF2_x[i]->DrawClone("PE X0 same");    
	//hF2_x2[i]->DrawClone("PE X0 same");
    }
/*    legF22->Draw("same");
    legF2->Draw("same");
    lat.DrawLatex(-2.8,0.4,scenario);
    gPad->SetLogy(1);

*/
    cnt=1;
    for(int i = 0;i<nbins_q2-2;i++){
        for(int j = 1; j<hF2_x[i]->GetNbinsX()+1;j++){
            double err = hF2_x[i]->GetBinError(j);
            double x = hF2_x[i]->GetBinCenter(j);
            double y = hF2_x[i]->GetBinContent(j);
	    hF2_x[i]->SetBinError(j,y*0.03);
	    //if(i==5 && (j > hF2_x[i]->GetNbinsX()-3))continue;
            //if(i==6 && (j > hF2_x[i]->GetNbinsX()-2))continue;
            //if(i==7 && (j > hF2_x[i]->GetNbinsX()-2))continue;
	    if(y>0){
                gF2->SetPoint(cnt,x,y);
                gF2->SetPointError(cnt,0,err);
                cnt++;
            }
        }
    }
    TCanvas *cf2 = new TCanvas("cf2","cf2"); 
    gF2->GetXaxis()->SetTitle("log(x_{B})");
    gF2->GetYaxis()->SetTitle("F_{2}^{c#bar{c}}(x_{B},Q^{2}) #times C");
    gF2->GetYaxis()->SetLimits(0.0000001,5);
    gF2->GetYaxis()->SetRangeUser(0.00002,0.2);
    gF2->GetXaxis()->SetLimits(-6.9,0);
    gF2->GetXaxis()->SetRangeUser(-3.,0);
    gF2->Draw("APE");
    for(int i = 0;i<nbins_q2-2;i++){
	hF2_x[i]->SetLineColor(1);
	hF2_x[i]->SetMarkerColor(1);
	hF2_x[i]->SetMarkerStyle(8);
        hF2_x[i]->SetLineWidth(0.2);
	hF2_x[i]->SetFillColor(1);
	hF2_x[i]->Draw("E5 ][ same");
    }    
    gPad->SetLogy(1);
    TLegend *leg_f2 = new TLegend(0.15,0.85,0.4,0.9);
    leg_f2->SetTextSize(0.045);
    leg_f2->SetHeader("PYTHIA6 e+p #sqrt{s} #approx 63,29 GeV");
    leg_f2->Draw("same");
    double latx1[nbins_q2-2]={-3.45,-3.4,-3.25,-2.95,-2.75,-2.5,-2.3,-2.1,-1.88,-1.65,-1.45};
    double laty1[nbins_q2-2]={0.00012,0.00035,0.0008,0.002,0.004,0.007,0.012,0.02,0.032,0.04,0.06};
    double lata1[nbins_q2-2]={0,0,0,-8,-13,-15,-15,-17,-20,-20,-25};

    for(int i = 0;i<nbins_q2-2;i++){
	lat.SetTextAngle(lata1[i]);
        //if(i<3)lat.DrawLatex(latx1[i],laty1[i],Form("#color[4]{Q^{2}=%1.1f GeV^{2}, C=%1.3f}",qs[i],scale[i]));
        //else if(i<5) lat.DrawLatex(latx1[i],laty1[i],Form("#color[4]{Q^{2}=%1.f GeV^{2}, C=%1.2f}",qs[i],scale[i]));
        //else lat.DrawLatex(latx1[i],laty1[i],Form("#color[4]{Q^{2}=%1.f GeV^{2}, C=%1.1f}",qs[i],scale[i]));
	if(i<3)lat.DrawLatex(latx1[i],laty1[i],Form("#color[4]{%1.1f, %1.3f}",qs[i],scale[i]));
        else if(i<5) lat.DrawLatex(latx1[i],laty1[i],Form("#color[4]{%1.f, %1.2f}",qs[i],scale[i]));
        else lat.DrawLatex(latx1[i],laty1[i],Form("#color[4]{%1.f, %1.1f}",qs[i],scale[i]));
    }
    lat.SetTextAngle(0);
    lat.DrawLatex(-0.8,0.1,"#color[4]{Q^{2} (GeV^{2}), C}");
    /*TCanvas *cdf2 = new TCanvas("cdf2","cdf2",2000,1200);
    cdf2->Divide(4,3);
    for(int i = 0;i<nbins_q2-2;i++){
	cdf2->cd(i+1);
	hdF2_x[i]->GetYaxis()->SetRangeUser(0.001,1);
	hdF2_x[i]->SetLineColor(1);
	hdF2_x[i]->SetMarkerColor(1);
	hdF2_x[i]->SetMarkerStyle(8);
	hdF2_x[i]->DrawClone("PE X0 same");
	gPad->SetLogy(1);
	//gPad->SetGrid(1,1);
	//hsdF2_x[i]->DrawClone("hist same");
	//hs1dF2_x[i]->DrawClone("hist same");
	lat.SetTextSize(0.09);
	lat.SetTextAngle(0);
	lat.DrawLatex(-3.8,0.45,Form("Q^{2}=%1.f GeV^{2}",qs[i]));
	//gStyle->SetPadLeftMargin(0.2);
//leg2->Draw("same");
	}
    */
//    lat.DrawLatex(-3.8,hdF2_x[0]->GetMaximum()*0.9,scenario);
/*
    TF1 *line = new TF1("line","1",-100,100);
    line->SetLineStyle(7);

    TCanvas *cf2r = new TCanvas("cf2r","cf2r",2000,1200);                                    
    cf2r->Divide(4,3);
    for(int i = 0;i<nbins_q2-2;i++){
	cf2r->cd(i+1);
	hdF2_x[i]->GetXaxis()->SetRangeUser(-3.7,0);//hF2_x[0]->GetMaximum()*5);                                                                                                                                                                                
	hdF2_x[i]->GetYaxis()->SetRangeUser(0.8,1.3);//hF2_x[0]->GetMaximum()*5);                                                                                                                                                                             
	hdF2_x[i]->GetYaxis()->SetTitle("F_{2}^{c#bar{c}} Error Ratio");
	getRat(hdF2_x[i],hdF2_x2[i]);
	hdF2_x[i]->SetLineColor(1);
	hdF2_x[i]->DrawClone("PE X0 same");
	line->Draw("l same");
	lat.SetTextSize(0.08);
	lat.DrawLatex(-3.7,1.24,Form("Q^{2}=%1.f GeV^{2}",qs[i]));
    }
*/    
    TCanvas *cf2r2 = new TCanvas("cf2r2","cf2r2");
    hdF2_x[0]->SetMarkerColor(1);
    hdF2_x[0]->SetLineColor(1);
    hdF2_x[0]->SetMarkerStyle(8);
    hdF2_x[5]->SetMarkerColor(kRed);
    hdF2_x[5]->SetLineColor(kRed);
    hdF2_x[5]->SetMarkerStyle(22);
    hdF2_x[9]->SetMarkerColor(kBlue);
    hdF2_x[9]->SetLineColor(kBlue);
    hdF2_x[9]->SetMarkerStyle(21);
    hdF2_x[0]->GetYaxis()->SetRangeUser(0.002,1);
    hdF2_x[0]->GetXaxis()->SetRangeUser(-3,0);
    TGraph *gd1 = new TGraph();
    TGraph *gd2 = new TGraph();
    TGraph *gd3 = new TGraph();
    int cnt1=0;
    int cnt2=0;
    int cnt3=0;
    for(int i = 1;i < hdF2_x[0]->GetNbinsX();i++){
	double v1 = hdF2_x[0]->GetBinContent(i);
	double v2 = hdF2_x[5]->GetBinContent(i);
	double v3 = hdF2_x[9]->GetBinContent(i);
	if(v1<1 && v1 > 0.002 ){
	    cout << " HERE " << v1 << endl;
	    gd1->SetPoint(gd1->GetN(),hdF2_x[0]->GetBinCenter(i),v1);
	    cnt1++;
	}	
	if(v2<1 && v2 > 0.002){
	    gd2->SetPoint(gd2->GetN(),hdF2_x[0]->GetBinCenter(i),v2);
	    cnt2++;
	}
	if(v3<1 && v3 > 0.002){
	    gd3->SetPoint(gd3->GetN(),hdF2_x[0]->GetBinCenter(i),v3);
	    cnt3++;
	}
    }
    gd1->SetLineColor(1);
    gd2->SetLineColor(kRed);
    gd3->SetLineColor(kBlue);
TLegend *legg = new TLegend(0.16,0.68,0.38,0.92);
    sprintf(label,"Q^{2} = %1.f GeV^{2}",qs[0]);
    legg->AddEntry(hdF2_x[0],label,"PE");
    sprintf(label,"Q^{2} = %1.f GeV^{2}",qs[5]);
    legg->AddEntry(hdF2_x[5],label,"PE");
    sprintf(label,"Q^{2} = %1.f GeV^{2}",qs[9]);
    legg->AddEntry(hdF2_x[9],label,"PE");

    hdF2_x[0]->Draw("PE X0 same");
    gd1->Draw("l same");
    gd2->Draw("l same");
    gd3->Draw("l same");
    //line->Draw("same l");
    hdF2_x[5]->Draw("PE X0 same");
    hdF2_x[9]->Draw("PE X0 same");
    legg->Draw("same");
    gPad->SetLogy(1);        
//lat.DrawLatex(-3.6,0.85,"Detector Matrix/LBNL Vertexing");
    













}
double getRCSVsX(TH1F *means_q2, TH1F *means_x, TH1F * fill, TH1F *dfill,TH3F *h3d,TH3F *h3d_d,TH3F *h3d_bkg,int qbin, double xsec, int ss, double L){
    //cout << "==> Now in getRCSVsX " << endl;
    double width = 1;
    if(qbin!=0){
	//cout << " => Setting range from " << binning_q2[qbin-1] << " to " << binning_q2[qbin]<< endl;
        h3d->GetYaxis()->SetRangeUser(binning_q2[qbin-1],binning_q2[qbin]);
	h3d_d->GetYaxis()->SetRangeUser(binning_q2[qbin-1],binning_q2[qbin]);
	h3d_bkg->GetYaxis()->SetRangeUser(binning_q2[qbin-1],binning_q2[qbin]);
	double mean_Q2 = TMath::Power(10,h3d->Project3D("y")->GetMean());
	width = fabs(TMath::Power(10,binning_q2[qbin])-TMath::Power(10,binning_q2[qbin-1]));
    }
    //cout<<" => mean_q2 in reagion is " << mean_q2 << " and width " << width << endl;
    double up = L/xsec;
    double lumi = L * 1E+12; // in mb       
    for(int i = 1; i< nbins_x+1; i++){
	//Set Ranges in 3D histograms
        h3d->GetXaxis()->SetRangeUser(binning_x[i-1],binning_x[i]);
	h3d_d->GetXaxis()->SetRangeUser(binning_x[i-1],binning_x[i]);
	h3d_bkg->GetXaxis()->SetRangeUser(binning_x[i-1],binning_x[i]);
	double mean_x  = TMath::Power(10,h3d->Project3D("x")->GetMean());
        double mean_q2 = TMath::Power(10,h3d->Project3D("y")->GetMean());
	if(mean_x==0||mean_q2==0){
	    means_x->SetBinContent(i,0);
	    means_q2->SetBinContent(i,0);
	    if(nbins_x==20)cout << "\n\n\n\n ====HERE==== " << qbin << endl;
	    continue;
	}
	means_x->SetBinContent(i,mean_x);
        means_q2->SetBinContent(i,mean_q2);	
	char name[20];
        sprintf(name,"hh%i",i);
        temp = (TH1D*) h3d->Project3D("z");temp->SetName(name);
	temp1 = (TH1D*) h3d_d->Project3D("z");
	temp2 = (TH1D*) h3d_bkg->Project3D("z");
	if(temp->Integral()<=0 || temp1->Integral()<=0)continue;
	double eff = temp->Integral(-1,-1)/temp1->Integral(-1,-1);
	double val = temp->Integral(-1,-1);
	double val2 = temp2->Integral(-1,-1) - temp->Integral(-1,-1);
	if(0){//ss==2){
	    cout <<" HERE " << i<< " qbin " << qbin << " " << val << " " << temp2->Integral(-1,-1) << " " << val2 << " " << val/sqrt(val+val2) << endl;
	}
//cout <<"HERE" << mean_x << " " << mean_q2 <<" " << val << " " << val2 << endl;
	if(val2<0){
	    cout <<"Uh oh with the yields " << endl;
	    val2=1;
//continue;
	}
	if(val/sqrt(val+val2)<1)continue;
	double s = 0;
	if(ss==1)s = 4.*18.*275.;
	if(ss==2)s = 4.*10.*100.;
	if(ss==3)s = 4.*5.*41.;
        double mean_y = mean_q2/s/mean_x;
	double conv = 1./0.565577/0.3894;
	double factor = mean_x * TMath::Power(mean_q2,2)/(2*3.141592/137.06/137.06*(1+(1-mean_y)*(1-mean_y)));
        double scale = 0.5 * conv * 1./eff * 1./BR * 1./lumi * 1./fabs(TMath::Power(10,binning_x[i])-TMath::Power(10,binning_x[i-1])) * 1./width * factor;
	fill->SetBinContent(i,up*val*scale);
	fill->SetBinError(i,sqrt(up)*scale*sqrt(val+val2));
	dfill->SetBinContent(i,((sqrt(up)*sqrt(val+val2))/(up*val)));
	dfill->SetBinError(i,0);
        delete temp1;delete temp; 
    }
    //Next ad-hoc method to remove first data point if at edge of acceptance producing ugly dip in distribution 
    double tempp=0;
    for(int i = 1; i<fill->GetNbinsX()+1;i++){
	double val = fill->GetBinContent(i);
	if(val>0 && tempp==0)tempp=val;
	else if(tempp<val){
	    fill->SetBinContent(i-1,0);
	    dfill->SetBinContent(i-1,0);
	    break;
	}
    }

    fill->GetXaxis()->SetRangeUser(-4,0);
    fill->GetXaxis()->SetTitle("log(x_{B})");
    fill->GetYaxis()->SetTitle("#sigma_{r}^{c#overline{c}}(x_{B},Q^{2}) #times C");
    fill->GetYaxis()->SetNdivisions(505);
    dfill->GetYaxis()->SetTitle("Relative Error on #sigma_{r}^{c#overline{c}}");//S/#sqrt{S+B}");
    //dfill->GetYaxis()->SetRangeUser(0,0.1);
    dfill->GetXaxis()->SetTitle("log(x_{B})");
    dfill->GetXaxis()->SetRangeUser(-4,0);

    h3d->GetXaxis()->SetRangeUser(-100,100);
    h3d->GetYaxis()->SetRangeUser(-100,100);
    h3d_d->GetXaxis()->SetRangeUser(-100,100);
    h3d_d->GetYaxis()->SetRangeUser(-100,100);
    fill->SetLineColor(colors[qbin-1]);
    fill->SetMarkerColor(colors[qbin-1]);
    fill->SetMarkerStyle(markers[qbin-1]);
    dfill->SetLineColor(colors[qbin-1]);
    dfill->SetMarkerColor(colors[qbin-1]);
    dfill->SetMarkerStyle(markers[qbin-1]);    
    return mean_Q2;
}
double getF2VsX(TH1F * fill,TH1F * dfill, TH1F * dsfill,TH1F * ds1fill,TH1F *cs1,TH1F *mean_x1, TH1F* mean_q21,  TH1F *cs2,TH1F *mean_x2, TH1F* mean_q22,  TH1F *cs3,TH1F *mean_x3, TH1F* mean_q23,int qbin,TCanvas *cc){
    TLegend *leggg = new TLegend(0.65,0.65,0.9,0.9);
    double s1 = 4.*18.*275.;
    double s2 = 4.*10.*100.;
    double s3 = 4.*5.*41.;
    int cnt=0;
    for(int i = 1; i< nbins_x+1; i++){
        
	double mean_y1 = mean_q21->GetBinContent(i)/s1/mean_x1->GetBinContent(i);
	double mean_y2 = mean_q22->GetBinContent(i)/s2/mean_x2->GetBinContent(i);
	double mean_y3 = mean_q23->GetBinContent(i)/s3/mean_x3->GetBinContent(i);
	double yy1 = mean_y1*mean_y1 / (1+ (1-mean_y1)*(1-mean_y1));
	double yy2 = mean_y2*mean_y2 / (1+ (1-mean_y2)*(1-mean_y2));
        double yy3 = mean_y3*mean_y3 / (1+ (1-mean_y3)*(1-mean_y3));
	double val1 = cs1->GetBinContent(i);
        double val2 = cs2->GetBinContent(i);
        double val3 = cs3->GetBinContent(i);
	if(0){//qbin==1){
	    cout << "\n<q2> " << mean_q21->GetBinContent(i) << " " << mean_q22->GetBinContent(i) << " " << mean_q23->GetBinContent(i) << endl;
	    cout << "<x> " << mean_x1->GetBinContent(i) << " " << mean_x2->GetBinContent(i) << " " << mean_x3->GetBinContent(i) << endl;
	    cout << "y " << mean_y1 << " " << mean_y2 << " " << mean_y3 << endl;
	    cout << "YY " << yy1 << " " << yy2 << " " << yy3 << "\n"<< endl;
	    cout << "val " << val1 << " " << val2 << " " << val3 << "\n"<< endl;
	}
	if(val1+val2==0 || val1+val3==0 || val2+val3==0)continue;//Make sure there are at least two points to fit
	if(val2==0 || val3==0)continue;//Using only two energies
	double er1 = cs1->GetBinError(i);
	double er2 = cs2->GetBinError(i);
	double er3 = cs3->GetBinError(i);
	/*if(val1==0){
	    double axis_x[3] = {yy2,yy3};
	    double axis_ex[3] = {0,0};
	    double axis_y[3] = {val2,val3};
	    double axis_ey[3] = {er2,er3};
	    TGraphErrors gr(2,axis_x,axis_y,axis_ex,axis_ey);
	}else if(val2==0){
	*/  
	    double axis_x[2] = {yy2,yy3};
            double axis_ex[2] = {0,0};
            double axis_y[2] = {val2,val3};
            double axis_ey[2] = {er2,er3};
	    double axis_ey1[2] = {er2*sqrt(2),er3};
	    double axis_ey2[2] = {er2,er3*sqrt(10)};

            TGraphErrors gr(2,axis_x,axis_y,axis_ex,axis_ey);
	    TGraphErrors gr1(2,axis_x,axis_y,axis_ex,axis_ey1);
	    TGraphErrors gr2(2,axis_x,axis_y,axis_ex,axis_ey2);
	    /*}else if(val3==0){
	    
	    double axis_x[3] = {yy1,yy2};
            double axis_ex[3] = {0,0};
            double axis_y[3] = {val1,val2};
            double axis_ey[3] = {er1,er2};
            TGraphErrors gr(2,axis_x,axis_y,axis_ex,axis_ey);
	}else{
	    double axis_x[3] = {yy1,yy2,yy3};
            double axis_ex[3] = {0,0,0};
            double axis_y[3] = {val1,val2,val3};
            double axis_ey[3] = {er1,er2,er3};
            TGraphErrors gr(3,axis_x,axis_y,axis_ex,axis_ey);
	}
	    */
	TF1 fit("fit","[0]+[1]*x",0,1);
	fit.SetLineStyle(7);
	if(cnt==0)fit.SetLineColor(1);
	if(cnt==1)fit.SetLineColor(kBlue);
	if(cnt==2)fit.SetLineColor(kRed);
	if(cnt==3)fit.SetLineColor(kGreen-2);
	fit.SetParLimits(0,0,0.5);
	gr.Fit("fit","M EX0");
	TLatex lat;
	gr.GetYaxis()->SetNdivisions(505);
	
	gr.GetXaxis()->SetLimits(-1,1);
	gr.GetYaxis()->SetLimits(0,1);
	gr.GetXaxis()->SetRangeUser(-0.05,axis_x[1]*1.2);
	gr.GetYaxis()->SetRangeUser(0,axis_y[0]*2);
	gr.GetYaxis()->SetTitle("#sigma_{r}^{c#bar{c}}(x_{B},Q^{2})");
	gr.GetXaxis()->SetTitle("y^{2}/Y^{+}");
	if(qbin==3 && cnt>-1 && cnt<4){
	    //cc->cd(i);
	    char label[100];
	    leggg->SetHeader("Q^{2} = 3.1 GeV^{2}");
	    sprintf(label,"x = %1.3f",mean_x2->GetBinContent(i));
	    TH1F *prox = new TH1F("prox","prox",10,0,10);
	    prox->SetMarkerSize(1.5);
	    if(cnt==0){
                gr.SetMarkerStyle(8);
                gr.SetLineColor(1);
                gr.SetMarkerColor(1);
                prox->SetMarkerStyle(8);
                prox->SetLineColor(1);
                prox->SetMarkerColor(1);
            }
	    if(cnt==1){
		gr.SetMarkerStyle(21);
		gr.SetLineColor(kBlue);
		gr.SetMarkerColor(kBlue);
		prox->SetMarkerStyle(21);
		prox->SetLineColor(kBlue);
		prox->SetMarkerColor(kBlue);
	    }
	    if(cnt==2){

                gr.SetMarkerStyle(22);
		gr.SetLineColor(kRed);
		gr.SetMarkerColor(kRed);
		prox->SetMarkerStyle(22);
		prox->SetLineColor(kRed);
                prox->SetMarkerColor(kRed);

            }
	    if(cnt==3){

                gr.SetMarkerStyle(23);
		gr.SetLineColor(kGreen-2);
		gr.SetMarkerColor(kGreen-2);
		prox->SetMarkerStyle(23);
		prox->SetLineColor(kGreen-2);
                prox->SetMarkerColor(kGreen-2);

            }
	    gr.SetMarkerSize(1.5);
	    leggg->AddEntry(prox,label,"p");    
	    if(cnt==0){
		cc->cd();
		gr.GetYaxis()->SetLimits(0,1);
		gr.GetXaxis()->SetLimits(-1,1);
		gr.GetYaxis()->SetRangeUser(0,0.05/0.35);
		gr.GetXaxis()->SetRangeUser(-0.05,0.6);
		gr.DrawClone("APE");
	    }
	    else gr.DrawClone("PE same");
	    //lat.DrawLatex(axis_x[1]*0.5,axis_y[0]*1.4,Form("Q2=%1.3f,x=%1.3f",mean_q22->GetBinContent(i),mean_x2->GetBinContent(i)));
	    leggg->Draw("same");
	}
	cnt++;
	if(fit.GetParameter(0)>0){
	    fill->SetBinContent(i,fit.GetParameter(0));
	    fill->SetBinError(i,fit.GetParError(0));
	    dfill->SetBinContent(i,fit.GetParError(0 )/fit.GetParameter(0));
	    dfill->SetBinError(i,0);
	}
	gr1.Fit("fit","M EX0");
	if(fit.GetParameter(0)>0){
	    dsfill->SetBinContent(i,fit.GetParError(0)/fit.GetParameter(0) );
            dsfill->SetBinError(i,0);
        }
        gr2.Fit("fit","M EX0");
        if(fit.GetParameter(0)>0){
            ds1fill->SetBinContent(i,fit.GetParError(0)/fit.GetParameter(0));
            ds1fill->SetBinError(i,0);
        }
    }
//    leggg->Draw("same");
    fill->SetLineColor(colors[qbin-1]);
    fill->SetMarkerColor(colors[qbin-1]);
    fill->SetMarkerStyle(markers[qbin-1]);
    dfill->SetLineColor(colors[qbin-1]);
    dfill->SetMarkerColor(colors[qbin-1]);
    dfill->SetMarkerStyle(markers[qbin-1]);
    dsfill->SetLineColor(colors[qbin-1]);
    dsfill->SetMarkerColor(colors[qbin-1]);
    dsfill->SetMarkerStyle(markers[qbin-1]);
    dsfill->SetLineStyle(7);
    ds1fill->SetLineColor(colors[qbin-1]);
    ds1fill->SetMarkerColor(colors[qbin-1]);
    ds1fill->SetMarkerStyle(markers[qbin-1]);
    ds1fill->SetLineStyle(8);

    fill->GetXaxis()->SetRangeUser(-4,0);
    fill->GetXaxis()->SetTitle("log(x_{B})");
    fill->GetYaxis()->SetTitle("F_{2}^{c#bar{c}}(x_{B},Q^{2}) #times C");
    fill->GetYaxis()->SetNdivisions(505);
    dfill->GetXaxis()->SetRangeUser(-4,0);
    dfill->GetXaxis()->SetTitle("log(x_{B})");
    dfill->GetYaxis()->SetTitle("#sigma(F_{2}^{c#bar{c}})/F_{2}^{c#bar{c}} (Stat.)");
    dfill->GetYaxis()->SetNdivisions(505);
    //dfill->GetXaxis()->SetTitleSize(0.08);
    //dfill->GetYaxis()->SetTitleSize(0.08);

    return 1;
}
void norm(TH1F *h, char label[]){
    double norm = h->Integral();
    h->GetXaxis()->SetTitle(label);
    h->GetYaxis()->SetTitle("Arb. Units");
    
    for(int i = 1; i < h->GetNbinsX()+1;i++){
	double v = h->GetBinContent(i);
	double e = h->GetBinError(i);
	if(v>0){
	    h->SetBinContent(i,v/norm);
	    h->SetBinError(i,e/norm);
	}
	else{
	    h->SetBinContent(i,0);
            h->SetBinError(i,0);
	}
    }
    h->GetYaxis()->SetRangeUser(0.0,h->GetMaximum()*1.2);
}

double getXsec(TH1F *h){

    double sum=0;
    for(int i = 1;i<h->GetNbinsX()+1;i++){
	sum+= h->GetBinCenter(i)*h->GetBinContent(i);
    }
    cout << "Total events " << sum << endl;
    sum *= 1/0.94 * 1E-9;//0.000001;
    cout << "Int. lumi " << sum << endl;

    return sum;
}
double getXsec(TH2F *h, double conv){

    double sum=h->Integral();

    cout << "Total events " << sum << endl;
    sum *= conv * 1E-9;//0.000001;                                                                                                                                                                                                                                      
    cout << "Int. lumi " << sum << endl;

    return sum;
}
void getRat(TH1F* h, TH1F* h1){
    for(int i = 1;i<h->GetNbinsX()+1;i++){
	double val1 = h->GetBinContent(i);
	double er1 = h->GetBinError(i);
	double val2 = h1->GetBinContent(i);
	double er2 = h1->GetBinError(i);
	if(val1>0 && val2>0){
	    if(val1/val2<1.5){
		h->SetBinContent(i,val1/val2);
		h->SetBinError(i,val1/val2*sqrt(er1*er1/val1/val1+er2*er2/val2/val2));
	    }
	}
    }
}
