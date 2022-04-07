void CanvasPartition(TCanvas *C,const Int_t Nx = 2,const Int_t Ny = 2,
                     Float_t lMargin = 0.15, Float_t rMargin = 0.05,
                     Float_t bMargin = 0.15, Float_t tMargin = 0.05);

void setXbins(int nbins_x, double binning_x[]){
    for(int i = 0;i<nbins_x+1;i++)
        binning_x[i] = -4. + i * 4. / nbins_x;
}
void setQ2bins(int nbins_q2, double binning_q2[]){
    for(int i = 0;i<nbins_q2+1;i++)
        binning_q2[i] = i * 2.6 / nbins_q2;
}
void Plot_Alt(int save=0){
    gROOT->ProcessLine(".x ~/myStyle.C");
    int const nbins_q2 = 13;
    int const nbins_x  = 20;
    double binning_q2[nbins_q2+1];
    double binning_x[nbins_x+1];
    setXbins(nbins_x,binning_x);
    setQ2bins(nbins_q2,binning_q2);
    double qq[20] = {1.3,2,3.1,5,8,12,20,31,52,76,120};
    double scale[nbins_q2]={0.002,0.008,0.01,0.008,0.5,1.5,3,5,10,20,55,100,400};
    TH1F *hCS_x1[100];
    TH1F *hCS_x2[100];
    TH1F *hF2_x[100];
    TH1F *hCS_x1_er[100];
    TH1F *hCS_x2_er[100];
    TH1F *hF2_x_er[100];


    TH1F *hF2_x_IC1[100];
    TH1F *hF2_x_IC2[100];


    TFile *fin1 = new TFile("ep_output_ATHENA.root","READ");
    for(int i = 0;i<nbins_q2-2;i++){
	hF2_x_IC1[i]= (TH1F*)fin1->Get(Form("ep_541_Q2Bin_%i",i));
	hF2_x_IC1[i]->SetName(Form("hF2_x_IC1_%i",i));
	if(i==9)hF2_x_IC1[i]->SetBinContent(20,0);
    }


    TFile *fin3 = new TFile("ep_output_ATHENA.root","READ");
    for(int i = 0;i<nbins_q2-2;i++){//nodifying the index                                                                                                                                                                                                         
	hF2_x[i]= (TH1F*)fin3->Get(Form("ep_10100_Q2Bin_%i",i));
    }

    TF1 *line = new TF1("line","1",-100,100);
    line->SetLineStyle(1);
    TLatex lat;
    TH1F *hF2_xx[100];
    TH1F *hF2_xx_er[100];

    TH1F *hF2_x_IC11[100];
    TH1F *hF2_x_IC22[100];
    for(int i = 0;i <100;i++){
	hF2_x_IC11[i] = new TH1F(Form("hF2_x_IC11_%i",i),Form("hF2_x_IC11_%i",i),nbins_q2,binning_q2);
	hF2_x_IC22[i] = new TH1F(Form("hF2_x_IC22_%i",i),Form("hF2_x_IC22_%i",i),nbins_q2,binning_q2);
	hF2_xx[i] = new TH1F(Form("hF2_xx_%i",i),Form("hF2_xx_%i",i),nbins_q2,binning_q2);
	hF2_xx_er[i] = new TH1F(Form("hF2_xx_er_%i",i),Form("hF2_xx_er_%i",i),nbins_q2,binning_q2);
    }
    
    for(int i = 0;i < nbins_q2-2 ; i++){
	for(int j = 0;j < nbins_x ; j++){
	    hF2_xx[j]->SetBinContent(i+1, hF2_x[i]->GetBinContent(j+1));
	    hF2_xx[j]->SetBinError(i+1, hF2_x[i]->GetBinError(j+1));
	    hF2_x_IC11[j]->SetBinContent(i+1, hF2_x_IC1[i]->GetBinContent(j+1));
	    hF2_x_IC11[j]->SetBinError(i+1, hF2_x_IC1[i]->GetBinError(j+1));
	    
	}
    }



    for(int i = 0;i < nbins_x ; i++){ 
	double fac = TMath::Power(4,(nbins_x-i));
	hF2_xx[i]->Scale(fac);//scale[i]);
	hF2_x_IC11[i]->Scale(fac);//scale[i]);
	hF2_x_IC11[i]->SetLineColor(kRed);
	hF2_x_IC11[i]->SetLineWidth(2);
	hF2_xx_er[i]->Scale(fac);//scale[i]);
    }
    TGraphErrors *gF2 = new TGraphErrors();
    TGraphErrors *gF2_IC1 = new TGraphErrors();
    
    TGraphAsymmErrors *gF2_er = new TGraphAsymmErrors();
    int cnt=1;
    for(int i = 0;i<nbins_x;i++){
        for(int j = 1; j<hF2_xx[i]->GetNbinsX()+1;j++){
	    double fac = 1;//TMath::Power(2,(i+1));
            double err = hF2_xx[i]->GetBinError(j)*fac;
	    double err2 = hF2_xx_er[i]->GetBinError(j)*fac;
	    double x = hF2_xx[i]->GetBinCenter(j);
            double y = hF2_xx[i]->GetBinContent(j)*fac;
	    double y1 = hF2_x_IC11[i]->GetBinContent(j)*fac;
	    if(y>0){
                gF2->SetPoint(cnt,x,y);
                gF2->SetPointError(cnt,0,err);
		gF2_er->SetPoint(cnt,x,y);
                if(j!=19)gF2_er->SetPointError(cnt,0.1,0.1,err2,err2);
		else if(j == 19 && i%2== 0 )gF2_er->SetPointError(cnt,0.1,0.08,err2,err2);
		else gF2_er->SetPointError(cnt,0.08,0.1,err2,err2);
		gF2_IC1->SetPoint(cnt,x,y1);
		cnt++;
            }
        }
    }
    double latx[nbins_q2-2]={-2.79,-2.8,-2.6,-2.3,-2.1,-1.95,-1.75,-1.55,-1.4,-1.2,-1.1};
    double laty[nbins_q2-2]={0.0007,0.015,0.2,2,13,100,500,2000,10000,50000,200000};
    gF2_er->SetLineColor(34);
    gF2_er->SetFillColor(33);
    gF2_er->SetFillStyle(1000);
    gF2_IC1->SetLineColor(kRed);
    TCanvas *cf2 = new TCanvas("cf2","cf2",1000,800);
    gPad->SetLeftMargin(0.2);
    gF2->GetXaxis()->SetTitle("log(Q^{2}/GeV^{2})");
    gF2->GetYaxis()->SetTitle("F_{2}^{c#bar{c}}(x_{B},Q^{2}) #times C");
    gF2->GetYaxis()->SetLimits(0.0000001,5);
    gF2->GetYaxis()->SetTitleOffset(1.3);
    gF2->GetYaxis()->SetRangeUser(0.0002,80000000);
    gF2->GetXaxis()->SetLimits(0,2.4);
    gF2->GetXaxis()->SetRangeUser(-3.,0);
    gF2->GetYaxis()->SetNdivisions(4);
    gF2->GetXaxis()->CenterTitle(1);
    gF2->GetYaxis()->CenterTitle(1);
    gF2->Draw("APE");
    //gF2_er->Draw("E5 same");
    gF2->Draw("PE same");
    for(int i = 6;i<nbins_x-1;i++){
	hF2_x_IC11[i]->DrawClone("hist ][ same");
	
	lat.SetTextSize(0.03);
	/*if(i<3)lat.DrawLatex(latx[i]-0.1,laty[i]*100,Form("#color[4]{%1.1f[%i]}",qq[i],i+1));//,scale[i]));
	else if(i<5) lat.DrawLatex(latx[i]-0.1,laty[i]*100,Form("#color[4]{%1.f[%i]}",qq[i],i+1));//,scale[i]));
	else lat.DrawLatex(latx[i]-0.1,laty[i]*100,Form("#color[4]{%1.f[%i]}",qq[i],i+1));//,scale[i]));
	*/
    }
    TLegend *leg = new TLegend(0.65,0.65,0.85,0.92);
    leg->SetTextSize(0.02);
    
    leg->SetHeader("PYTHIA+CT14 NNLO e+p #sqrt{s}=63,29 GeV");
    leg->AddEntry(gF2,"Projection Without IC","PE");
    //leg->AddEntry(gF2_er,"PDF Uncertainties","f");
    leg->AddEntry(hF2_x_IC11[1],"With IC: BHPS(1)","l");
    leg->Draw("same");
    lat.SetTextSize(0.04);
    lat.DrawLatex(-1,0.002*100,"#color[4]{Q^{2} (GeV/#it{c}^{2})[i]}");
    lat.DrawLatex(-1,0.0006*100,"#color[4]{C=7^{i}}");

    gPad->SetLogy();
    


    TLegend *leg1 = new TLegend(0.65,0.3,0.9,0.7);
    leg1->SetTextSize(0.04);
    leg1->SetHeader("PYTHIA e+p #sqrt{s}=63,29 GeV @ 10 fb^{-1}");
    leg1->AddEntry(gF2,"Projection With No IC","PE");
    leg1->AddEntry(gF2_er,"PDF Uncertainties","f");
    leg1->AddEntry(hF2_x_IC11[1],"With IC: BHPS","l");
    /*
    TCanvas *c11 = new TCanvas("c111","f2",800,800);
    // Canvas setup                                                                                                                                                                                                                                    
    const Int_t Nx1 = 1;
    const Int_t Ny1 = 3;
    CanvasPartition(c11,Nx1,Ny1,0.1,0.09,0.1,0.05);
    TPad *pad1[Nx1][Ny1];
    for (Int_t i=0;i<Nx1;i++) {
        for (Int_t j=0;j<Ny1;j++) {
            c11->cd(0);
            // Get the pads previously created.                                                                                                                                                                                                         
            char pname[16];
            sprintf(pname,"pad_%i_%i",i,j);
            pad1[i][j] = (TPad*) gROOT->FindObject(pname);
            pad1[i][j]->Draw();
            /// pad[i][j]->SetFillStyle(4000);                                                                                                                                                                                                         
            pad1[i][j]->SetFrameFillStyle(4000);
            pad1[i][j]->cd();
        }
    }

    for(int i = 0; i< 3; i++){
	int bin;
	if(i==0)bin=6;
	if(i==1)bin=8;
	if(i==2)bin=10;
	pad1[0][2-i]->cd();
	for(int ii = 1; ii< hF2_x[bin]->GetNbinsX()+1;ii++){
	    double val = hF2_x[bin]->GetBinContent(ii);
	    double er = hF2_x[bin]->GetBinError(ii);
	    double val1 = hF2_x_IC1[bin]->GetBinContent(ii);
	    double val2 = hF2_x_IC2[bin]->GetBinContent(ii);
	    double er2 =hF2_x_er[bin]->GetBinError(ii);
	    hF2_x_er[bin]->SetLineColor(34);
	    hF2_x_er[bin]->SetFillColor(33);
	    hF2_x_er[bin]->SetFillStyle(1000);

	    if(val>0){
		hF2_x[bin]->SetBinContent(ii,val/val);
		hF2_x[bin]->SetBinError(ii,val/val*er/val);
		hF2_x_er[bin]->SetBinContent(ii,val/val);
		hF2_x_er[bin]->SetBinError(ii,val/val*er2/val);
		hF2_x_IC1[bin]->SetBinContent(ii,val1/val);
		hF2_x_IC1[bin]->SetBinError(ii,0);
		hF2_x_IC2[bin]->SetBinContent(ii,val2/val);
		hF2_x_IC2[bin]->SetBinError(ii,0);
	    }
	    else {
		hF2_x_er[bin]->SetBinContent(ii,0);
                hF2_x_er[bin]->SetBinError(ii,0);

		hF2_x_IC1[bin]->SetBinContent(ii,0);
		hF2_x_IC1[bin]->SetBinError(ii,0);
		hF2_x_IC2[bin]->SetBinContent(ii,0);
                hF2_x_IC2[bin]->SetBinError(ii,0);
	    }
	    hF2_x[bin]->GetXaxis()->SetTitle("log(x_{B})");
	    hF2_x[bin]->GetYaxis()->SetTitle("Ratio to No IC Scenario");
	    hF2_x[bin]->GetXaxis()->CenterTitle(1);
	    hF2_x[bin]->GetYaxis()->CenterTitle(1);
	    hF2_x[bin]->GetYaxis()->SetTitleOffset(0.4);
	    
	    hF2_x[bin]->GetYaxis()->SetRangeUser(0.6,7.2);
	    hF2_x[bin]->GetXaxis()->SetRangeUser(-3.,0);
	    hF2_x[bin]->GetYaxis()->SetNdivisions(4);

	}
	hF2_x[bin]->DrawClone("PE X0");
	hF2_x_er[bin]->DrawClone("E2 same");
	hF2_x_IC1[bin]->DrawClone("hist same ][");
	hF2_x_IC2[bin]->DrawClone("hist same ][");
	line->Draw("same");
	if(i==0)leg1->Draw("same");
	lat.SetTextSize(0.1);

	if(bin<3)lat.DrawLatex(-1.,2.7,Form("#color[1]{%1.1f}",qq[bin]));
	else lat.DrawLatex(-1.,2.7,Form("#color[1]{%1.f}",qq[bin]));
	gPad->RedrawAxis();
    }
    
*/ 
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
