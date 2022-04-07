#include "fitdata1D.C"
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
void ICCalc(int save=0){
    char outname[100]="ep_output.root";
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
//============================
    TFile *f1 = new TFile("out_10_100_IC_Q210.root","READ");
    h3d1 = (TH3F*) f1->Get("D0L_RedCS_Topo_T");h3d1->SetName("h3d1");h3d1->SetDirectory(0);
    h3d_bkg1 = (TH3F*) f1->Get("D0L_RedCS_Topo_BKG");h3d_bkg1->SetName("h3d_bkg1");h3d_bkg1->SetDirectory(0);
    h3d_d1 = (TH3F*) f1->Get("D0L_RedCS_T_Acc");h3d_d1->SetName("h3d_d1");h3d_d1->SetDirectory(0);
    
    
    TH2F* _hxQ21 = (TH2F*) f1->Get("hxQ2");_hxQ21->SetName("_hxQ21");
    CS1 = getXsec(_hxQ21,1/0.044);
    //CS1 = getXsec(_hxQ21,1/0.56);
    cout << "Calculated Lumi in 10x100 is " << CS1 << endl;
//============================  
    TFile *f2 = new TFile("./out_10_100_IC_Q210.root","READ");
    h3d2 = (TH3F*) f2->Get("D0L_nRedCS_Topo_T");h3d2->SetName("h3d2");h3d2->SetDirectory(0);
    h3d_bkg2 = (TH3F*) f2->Get("D0L_nRedCS_Topo_BKG");h3d_bkg2->SetName("h3d_bkg2");h3d_bkg2->SetDirectory(0);
    h3d_d2 = (TH3F*) f2->Get("D0L_nRedCS_T_Acc");h3d_d2->SetName("h3d_d2");h3d_d2->SetDirectory(0);
    TH2F* _hxQ22 = (TH2F*) f2->Get("hxQ2");_hxQ22->SetName("_hxQ22");
    CS2 = getXsec(_hxQ22,1/0.044);
    cout << "Calculated Lumi in IC is " << CS2 << endl;
//============================                      
    TFile *f1_10 = new TFile("out_10_100_IC_Q210.root","READ");
    h3d1_10 = (TH3F*) f1_10->Get("D0L_RedCS_Topo_T");h3d1_10->SetName("h3d1_10");h3d1_10->SetDirectory(0);
    h3d_bkg1_10 = (TH3F*) f1_10->Get("D0L_RedCS_Topo_BKG");h3d_bkg1_10->SetName("h3d_bkg1_10");h3d_bkg1_10->SetDirectory(0);
    h3d_d1_10 = (TH3F*) f1_10->Get("D0L_RedCS_T_Acc");h3d_d1_10->SetName("h3d_d1_10");h3d_d1_10->SetDirectory(0);
    TH2F* _hxQ21_10 = (TH2F*) f1_10->Get("hxQ2");_hxQ21_10->SetName("_hxQ21_10");
    CS1_10 = getXsec(_hxQ21_10,1/0.044);
    //CS1_10 = getXsec(_hxQ21_10,1/0.56);
    cout << "Calculated Lumi in 10x100 q2=10 is " << CS1_10 << endl;
//============================                                            
    TFile *f2_10 = new TFile("./out_10_100_IC_Q210.root","READ");
    h3d2_10 = (TH3F*) f2_10->Get("D0L_nRedCS_Topo_T");h3d2_10->SetName("h3d2_10");h3d2_10->SetDirectory(0);
    h3d_bkg2_10 = (TH3F*) f2_10->Get("D0L_nRedCS_Topo_BKG");h3d_bkg2_10->SetName("h3d_bkg2_10");h3d_bkg2_10->SetDirectory(0);
    h3d_d2_10 = (TH3F*) f2_10->Get("D0L_nRedCS_T_Acc");h3d_d2_10->SetName("h3d_d2_10");h3d_d2_10->SetDirectory(0);
    TH2F* _hxQ22_10 = (TH2F*) f2_10->Get("hxQ2");_hxQ22_10->SetName("_hxQ22_10");
    CS2_10 = getXsec(_hxQ22_10,1/0.044);
    cout << "Calculated Lumi in IC q2=10 is " << CS2_10 << endl;
//============================                                                                                                                                                                                                                             
    
    
    setXbins();
    setQ2bins();
    TH1F *hCS_x1[nbins_q2];        
    TH1F *hCSL0_x1[nbins_q2];
    TH1F *hdCS_x1[nbins_q2];
    TH1F *hmean_x1[nbins_q2];
    TH1F *hmean_q21[nbins_q2];


    TH1F *hCSL0_x2[nbins_q2];

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
        
    TH1F *hR_x[nbins_q2];

    for(int i = 0;i < nbins_q2; i++){
	char name[100];
        sprintf(name,"hCS_x1_%i",i);
        hCS_x1[i] = new TH1F(name,name,nbins_x,binning_x);
	sprintf(name,"hCSL0_x1_%i",i);
        hCSL0_x1[i] = new TH1F(name,name,nbins_x,binning_x);
	sprintf(name,"hmean_x1_%i",i);
        hmean_x1[i] = new TH1F(name,name,nbins_x,binning_x);
        sprintf(name,"hmean_q21_%i",i);
        hmean_q21[i] = new TH1F(name,name,nbins_x,binning_x);
        sprintf(name,"hdCS_x1_%i",i);
        hdCS_x1[i] = new TH1F(name,name,nbins_x,binning_x);

	sprintf(name,"hCS_x2_%i",i);
        hCS_x2[i] = new TH1F(name,name,nbins_x,binning_x);
	sprintf(name,"hCSL0_x2_%i",i);
        hCSL0_x2[i] = new TH1F(name,name,nbins_x,binning_x);
	sprintf(name,"hmean_x2_%i",i);
        hmean_x2[i] = new TH1F(name,name,nbins_x,binning_x);
        sprintf(name,"hmean_q22_%i",i);
        hmean_q22[i] = new TH1F(name,name,nbins_x,binning_x);
        sprintf(name,"hdCS_x2_%i",i);
        hdCS_x2[i] = new TH1F(name,name,nbins_x,binning_x);

	sprintf(name,"hR_x_%i",i);
	hR_x[i] = new TH1F(name,name,nbins_x,binning_x);


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
    }
    double qs[20];
    double qs1[20];
    double qs2[20];
    double qs11[20];
    double qs22[20];
    for(int i = 5;i<nbins_q2;i++){
	if(i<5){
	    qs1[i] = getRCSVsX(hmean_q21[i],hmean_x1[i],hCS_x1[i],hCSL0_x1[i],hdCS_x1[i],h3d1,h3d_d1,h3d_bkg1,i+1,CS1,2,10);
	    qs2[i] = getRCSVsX(hmean_q22[i],hmean_x2[i],hCS_x2[i],hCSL0_x2[i],hdCS_x2[i],h3d2,h3d_d2,h3d_bkg2,i+1,CS2,3,10);
	    //qs11[i] = getRCSVsX(hmean_q211[i],hmean_x11[i],hCS_x11[i],hdCS_x11[i],h3d11,h3d_d11,h3d_bkg11,i+1,CS11,2,10);
	    //qs22[i] = getRCSVsX(hmean_q222[i],hmean_x22[i],hCS_x22[i],hdCS_x22[i],h3d22,h3d_d22,h3d_bkg22,i+1,CS22,3,10);
	}
	else{
	    qs1[i] = getRCSVsX(hmean_q21[i],hmean_x1[i],hCS_x1[i],hCSL0_x1[i],hdCS_x1[i],h3d1_10,h3d_d1_10,h3d_bkg1_10,i+1,CS1_10,2,10);
            qs2[i] = getRCSVsX(hmean_q22[i],hmean_x2[i],hCS_x2[i],hCSL0_x2[i],hdCS_x2[i],h3d2_10,h3d_d2_10,h3d_bkg2_10,i+1,CS2_10,3,10);
            //qs11[i] = getRCSVsX(hmean_q211[i],hmean_x11[i],hCS_x11[i],hdCS_x11[i],h3d11_10,h3d_d11_10,h3d_bkg11_10,i+1,CS11_10,2,10);
            //qs22[i] = getRCSVsX(hmean_q222[i],hmean_x22[i],hCS_x22[i],hdCS_x22[i],h3d22_10,h3d_d22_10,h3d_bkg22_10,i+1,CS22_10,3,10);
	}
    }
    double qs_1[20];
    
    TLatex lat;
    TCanvas *c11 = new TCanvas("c11","c11");
    hCS_x2[5]->SetMarkerStyle(24);
    hCS_x2[5]->SetLineColor(kRed);
    hCS_x2[5]->SetMarkerColor(kRed);
    hCS_x1[5]->Draw("PE");
    hCS_x2[5]->Draw("PE same");


//==================================================== CS  ==========================       
    TGraphErrors *gCS1 = new TGraphErrors();
    TGraphErrors *gCS2 = new TGraphErrors();
    TGraphErrors *gF2 = new TGraphErrors();
    TGraphErrors *gR = new TGraphErrors();

    hCS_x1[0]->GetXaxis()->SetRangeUser(-3.3,0);//hF2_x[0]->GetMaximum()*5);                                                                                                                                                                          
    hCS_x1[0]->GetYaxis()->SetRangeUser(0.000001,5);//hF2_x[0]->GetMaximum()*5);          
    for(int i = 5;i<nbins_q2-2;i++){
	hCS_x1[i]->Scale(scale[i]);
	hCS_x1[i]->GetYaxis()->SetRangeUser(0.000001,5);
    }
    hCS_x2[0]->GetXaxis()->SetRangeUser(-3.3,0);//hF2_x[0]->GetMaximum()*5);                                                                                                                                                                         
    hCS_x2[0]->GetYaxis()->SetRangeUser(0.000001,5);//hF2_x[0]->GetMaximum()*5);        
    for(int i = 5;i<nbins_q2-2;i++){
	hCS_x2[i]->Scale(scale[i]);
        hCS_x2[i]->GetYaxis()->SetRangeUser(0.000001,5);
    }
//=================== All in one TGraph =====
    int cnt=1;
    for(int i = 5;i<nbins_q2-2;i++){
	for(int j = 1; j<hCS_x1[i]->GetNbinsX()+1;j++){
	    double err = hCS_x1[i]->GetBinError(j);
	    double x = hCS_x1[i]->GetBinCenter(j);
	    double y = hCS_x1[i]->GetBinContent(j);
	    if(y>0){
		gCS1->SetPoint(cnt,x,y);
		gCS1->SetPointError(cnt,0,err);
		cnt++;
	    }
	}
    }
    cnt=1;
    for(int i = 5;i<nbins_q2-2;i++){
	for(int j = 1; j<hCS_x2[i]->GetNbinsX()+1;j++){
            double err = hCS_x2[i]->GetBinError(j);
	    double x = hCS_x2[i]->GetBinCenter(j);
            double y = hCS_x2[i]->GetBinContent(j);
            
	    if(y>0){
                gCS2->SetPoint(cnt,x+0.04,y);
		gCS2->SetPointError(cnt,0,err);
		cnt++;
            }
	}
    }
    gCS2->SetMarkerStyle(25);gCS2->SetMarkerColor(kRed);gCS2->SetLineColor(kRed);
    cnt=1;
    for(int i = 5;i<nbins_q2-2;i++){
        for(int j = 1; j<hCS_x2[i]->GetNbinsX()+1;j++){
            double err = hCS_x1[i]->GetBinError(j);
            double x = hCS_x1[i]->GetBinCenter(j);
            double y = hCS_x1[i]->GetBinContent(j);

	    double err1 = hCS_x2[i]->GetBinError(j);
	    double y1 = hCS_x2[i]->GetBinContent(j);
	    if(y>0){
		y = y1/y;
		//err = sqrt(err*err/y/y+err1*err1/y1/y1);
		err = sqrt(err*err/y/y);
		hR_x[i]->SetBinContent(j,y);
		hR_x[i]->SetBinError(j,err);

                cnt++;
            }
        }
    }


    TLatex lat;
    
    TCanvas *c12 = new TCanvas("c12","c12");
    gCS1->GetXaxis()->SetTitle("log(x_{B})");
    gCS1->GetYaxis()->SetTitle("#sigma_{r}^{c#bar{c}}(x_{B},Q^{2}) #times C");
    gCS1->GetYaxis()->SetLimits(0.00000001,5);
    gCS1->GetYaxis()->SetRangeUser(0.00005,0.4);
    gCS1->GetXaxis()->SetLimits(-6.9,0);
    gCS1->GetXaxis()->SetRangeUser(-2.6,0);
    gCS1->Draw("APE");
    gCS2->Draw("same PE");
    
    for(int i = 5;i<nbins_q2-2;i++){
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
	hCS_x2[i]->SetLineStyle(7);
	//hCS_x2[i]->Scale(1.2);
    }
    gPad->SetLogy(1);
    TLegend *leg_cs = new TLegend(0.17,0.2,0.4,0.4);
    leg_cs->SetTextSize(0.045);
    leg_cs->SetHeader("PYTHIA6 e+p 10#times100 10 fb^{-1}");
    leg_cs->AddEntry(gCS1,"CT14NNLO","P");
    leg_cs->AddEntry(gCS2,"CT14NNLO w/ IC","P");
    lat.SetTextSize(0.04);
    double latx[nbins_q2-2]={-3.6,-3.5,-3.4,-3.0,-2.9,-2.46,-2.07,-2.0,-2.0,-1.8,-1.6};
    double laty[nbins_q2-2]={0.0002,0.00052,0.0013,0.003,0.006,0.012,0.023,0.038,0.065,0.1,0.18};
    double lata[nbins_q2-2]={0,0,0,-3,-5,-10,-10,-10,-10,-15,-15};
    for(int i = 5;i<nbins_q2-2;i++){
	cout << "Drawing " << i << endl;
	lat.SetTextAngle(lata[i]);
	if(i<3)lat.DrawLatex(latx[i],laty[i],Form("#color[4]{%1.1f, %1.3f}",qs1[i],scale[i]));
        else if(i<5) lat.DrawLatex(latx[i],laty[i],Form("#color[4]{%1.f, %1.2f}",qs1[i],scale[i]));
        else lat.DrawLatex(latx[i],laty[i],Form("#color[4]{%1.f, %1.1f}",qs1[i],scale[i]));
    }
    leg_cs->Draw("same");
    lat.SetTextAngle(0);
    lat.DrawLatex(-0.8,0.1,"#color[4]{Q^{2} (GeV^{2}), C}");

    TF1 *line = new TF1("line","1",-100,100);
    
    TCanvas *c13 = new TCanvas("c13","err",1500,700);
    // Canvas setup                                                                                                                                                                                                                                    
    const Int_t Nx1 = 3;
    const Int_t Ny1 = 2;
    CanvasPartition(c13,Nx1,Ny1,0.07,0.05,0.1,0.05);
    TPad *pad[Nx1][Ny1];
    for (Int_t i=0;i<Nx1;i++) {
        for (Int_t j=0;j<Ny1;j++) {
            c13->cd(0);
            // Get the pads previously created.                                                                                                                                                                                                         
            char pname[16];
            sprintf(pname,"pad_%i_%i",i,j);
            pad[i][j] = (TPad*) gROOT->FindObject(pname);
            pad[i][j]->Draw();
            /// pad[i][j]->SetFillStyle(4000);                                                                                                                                                                                                         
            pad[i][j]->SetFrameFillStyle(4000);
            pad[i][j]->cd();
        }
    }
    line->SetLineStyle(7);
    for(int i = 5;i<nbins_q2-2;i++){
	if(i<3+5)pad[i-5][1]->cd();
	else pad[i-5-3][0]->cd();
	if(i==2+5 || i ==5+5)hR_x[i]->GetXaxis()->SetRangeUser(-2.2,0);
	else hR_x[i]->GetXaxis()->SetRangeUser(-2.2,-0.3);
	hR_x[i]->GetYaxis()->SetRangeUser(0.5,4.5);
	hR_x[i]->GetYaxis()->SetNdivisions(5);
        hR_x[i]->GetXaxis()->SetTitle("log(x_{B})");
        hR_x[i]->GetYaxis()->SetTitle("R(#sigma^{c#bar{c}}_{r}) With IC/ No IC");
	hR_x[i]->GetXaxis()->CenterTitle(1);
	hR_x[i]->GetYaxis()->CenterTitle(1);
	
	hR_x[i]->DrawClone("PE X0 same");
	line->Draw("same");
	lat.DrawLatex(-3.2,2.7,"e+p 10#times100 GeV");
        lat.DrawLatex(-3.2,2.5,Form("Q^{2} = %1.f GeV^{2}",qs1[i]));

    }
}
double getRCSVsX(TH1F *means_q2, TH1F *means_x, TH1F * fill,TH1F * fill1, TH1F *dfill,TH3F *h3d,TH3F *h3d_d,TH3F *h3d_bkg,int qbin, double xsec, int ss, double L){
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
        TH1F* temp = (TH1F*) h3d->Project3D("z");temp->SetName(name);
	TH1F* temp1 = (TH1F*) h3d_d->Project3D("z");
	TH1F* temp2 = (TH1F*) h3d_bkg->Project3D("z");
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
	fill1->SetBinContent(i,up* val*scale);                                                                                                                                                                                                   
	fill1->SetBinError(i,((sqrt(up)*sqrt(val+val2))/(up*val)));

	char pdffile[100];
        sprintf(pdffile,"fits/%1.2f_x_%1.2f_Q2.pdf",mean_x,mean_q2);
	for(int j = 1;j<temp->GetNbinsX()+1;j++){
	    double errr = temp->GetBinError(j);
	    temp->SetBinError(j,errr*sqrt(val+val2)/sqrt(val));
	}
	//double fit_v[2];
	//fit_v = 0;//fitdata1D(temp,temp,temp,0,pdffile,log10(mean_x),mean_q2);
	//double vv1 = *(fit_v+0);
	//double vv2 = *(fit_v+1);
	//cout << "HEERE " << up* vv1*scale << " " << *(fit_v+1) << endl;
	//if(vv1>0 && vv1<10e6){
	    //fill1->SetBinContent(i,up* vv1*scale);
	    //fill1->SetBinError(i,sqrt(up)*scale* vv2);
	//}
	
	delete temp2;delete temp1;delete temp; 
    }
    //Next ad-hoc method to remove first data point if at edge of acceptance producing ugly dip in distribution 
    /*   double tempp=0;
    for(int i = 1; i<fill->GetNbinsX()+1;i++){
	double val = fill->GetBinContent(i);
	if(val>0 && tempp==0)tempp=val;
	else if(tempp<val){
	    fill->SetBinContent(i-1,0);
	    dfill->SetBinContent(i-1,0);
	    break;
	}
	}*/
//    fix(fill1,dfill);
//    fix(fill1,dfill);
/*    fix(fill1,dfill);
    fix(fill1,dfill);
    fix(fill1,dfill);
    fix(fill1,dfill);
    fix(fill1,dfill);
    fix(fill1,dfill);
    fix(fill1,dfill);
    fix(fill1,dfill);
*/
    fix(fill,dfill);
    fix(fill,dfill);
    fix(fill,dfill);
    fill->GetXaxis()->SetRangeUser(-4,0);
    fill->GetXaxis()->SetTitle("log(x_{B})");
    fill->GetYaxis()->SetTitle("#sigma_{r}^{c#bar{c}}(x_{B},Q^{2}) #times C");
    fill->GetYaxis()->SetNdivisions(505);
    dfill->GetYaxis()->SetTitle("Relative Error on #sigma_{r}^{c#bar{c}}");//S/#sqrt{S+B}");
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
void fix(TH1F* fill, TH1F* dfill){
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
void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin)
{
    if (!C) return;

    // Setup Pad layout:                                                                                                                                                                                                                                
    Float_t vSpacing = 0.0;
    Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;

    Float_t hSpacing = 0.0;
    Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;

    Float_t vposd,vposu,vmard,vmaru,vfactor;
    Float_t hposl,hposr,hmarl,hmarr,hfactor;

    for (Int_t i=0;i<Nx;i++) {

        if (i==0) {
            hposl = 0.0;
            hposr = lMargin + hStep;
            hfactor = hposr-hposl;
            hmarl = lMargin / hfactor;
            hmarr = 0.0;
        } else if (i == Nx-1) {
            hposl = hposr + hSpacing;
            hposr = hposl + hStep + rMargin;
            hfactor = hposr-hposl;
            hmarl = 0.0;
            hmarr = rMargin / (hposr-hposl);
        } else {
            hposl = hposr + hSpacing;
            hposr = hposl + hStep;
            hfactor = hposr-hposl;
            hmarl = 0.0;
            hmarr = 0.0;
        }

        for (Int_t j=0;j<Ny;j++) {
	    if (j==0) {
                vposd = 0.0;
                vposu = bMargin + vStep;
                vfactor = vposu-vposd;
                vmard = bMargin / vfactor;
                vmaru = 0.0;
            } else if (j == Ny-1) {
                vposd = vposu + vSpacing;
                vposu = vposd + vStep + tMargin;
                vfactor = vposu-vposd;
                vmard = 0.0;
                vmaru = tMargin / (vposu-vposd);
            } else {
                vposd = vposu + vSpacing;
                vposu = vposd + vStep;
                vfactor = vposu-vposd;
                vmard = 0.0;
                vmaru = 0.0;
            }

            C->cd(0);

            char name[16];
            sprintf(name,"pad_%i_%i",i,j);
            TPad *pad = (TPad*) gROOT->FindObject(name);
            if (pad) delete pad;
            pad = new TPad(name,"",hposl,vposd,hposr,vposu);
            pad->SetLeftMargin(hmarl);
            pad->SetRightMargin(hmarr);
            pad->SetBottomMargin(vmard);
            pad->SetTopMargin(vmaru);

            pad->SetFrameBorderMode(0);
            pad->SetBorderMode(0);
            pad->SetBorderSize(0);

            pad->Draw();
        }
    }
}
