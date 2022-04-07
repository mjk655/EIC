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
void PlotCS(int save=0){
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
    

    TH1F *hF2_x_IC1[100];
    TH1F *hF2_x_IC2[100];


    TFile *fin1 = new TFile("ep_output_ATHENA.root","READ");
    for(int i = 0;i<nbins_q2-2;i++){
	hF2_x_IC1[i]= (TH1F*)fin1->Get(Form("ep_541_Q2Bin_%i",i));
	hF2_x_IC1[i]->SetName(Form("hF2_x_IC1_%i",i));
    }

    TFile *fin3 = new TFile("ep_output_ATHENA.root","READ");
    for(int i = 0;i<nbins_q2-2;i++){//nodifying the index                                                                                                                                                                                                         
	hF2_x[i]= (TH1F*)fin3->Get(Form("ep_10100_Q2Bin_%i",i));
    }

    TF1 *line = new TF1("line","1",-100,100);
    line->SetLineStyle(1);
    TLatex lat;
    
    

    for(int i = 0;i < nbins_q2-2 ; i++){ 
	double fac = TMath::Power(7,(i+1));
	hF2_x[i]->Scale(fac);//scale[i]);
	hF2_x_IC1[i]->Scale(fac);//scale[i]);
	hF2_x_IC1[i]->SetLineColor(kRed);
	hF2_x_IC1[i]->SetMarkerColor(kRed);
	hF2_x_IC1[i]->SetMarkerStyle(25);
	hF2_x_IC1[i]->SetMarkerSize(1.2);
	hF2_x[i]->GetYaxis()->SetRangeUser(0.0002,100000);
    }
    TGraphErrors *gF2 = new TGraphErrors();
    TGraphErrors *gF2_IC1 = new TGraphErrors();
    
    int cnt=1;
    for(int i = 0;i<nbins_q2-2;i++){
        for(int j = 1; j<hF2_x[i]->GetNbinsX()+1;j++){
	    double fac = 1;//TMath::Power(7,(i+1))*0.01;
            double err = hF2_x[i]->GetBinError(j)*fac;
	    double x = hF2_x[i]->GetBinCenter(j);
            double y = hF2_x[i]->GetBinContent(j)*fac;
	    double y1 = hF2_x_IC1[i]->GetBinContent(j)*fac;
	    double err1 = hF2_x_IC1[i]->GetBinError(j)*fac;
	    if(y>0){
                gF2->SetPoint(cnt,x,y);
                gF2->SetPointError(cnt,0,err);
		gF2_IC1->SetPoint(cnt,x,y1);
		gF2_IC1->SetPoint(cnt,x,err1);
                cnt++;
            }
        }
//	hF2_x_IC1[i]->Scale(1.2);
    }
    double latx[nbins_q2-2]={-2.79,-2.85,-2.65,-2.35,-2.15,-2.0,-1.8,-1.6,-1.45,-1.3,-1.2};
    double laty[nbins_q2-2]={0.0005,0.01,0.15,1,10,80,400,2000,17000,70000,300000};
    TCanvas *cf2 = new TCanvas("cf2","cf2",800,800);
    gPad->SetLeftMargin(0.2);
    gF2_IC1->SetLineColor(kRed);
    gF2_IC1->SetMarkerColor(kRed);
    gF2_IC1->SetMarkerStyle(25);
//    gF2_IC1->SetLineStyle(7);
    gF2->GetXaxis()->SetTitle("log_{10}(x_{B})");
    gF2->GetYaxis()->SetTitle("#sigma_{r}^{c#bar{c}}(x_{B},Q^{2}) #times C");
    gF2->GetYaxis()->SetLimits(0.0000001,5);
    gF2->GetYaxis()->SetTitleOffset(1.3);
    gF2->GetYaxis()->SetRangeUser(0.001,200000000);
    gF2->GetXaxis()->SetLimits(-6.9,0);
    gF2->GetXaxis()->SetRangeUser(-3.6,0);
    gF2->GetYaxis()->SetNdivisions(4);
    gF2->GetXaxis()->CenterTitle(1);
    gF2->GetYaxis()->CenterTitle(1);
    gF2->Draw("APE");
//    gF2_IC1->Draw("PE same");
    for(int i = 0;i<nbins_q2-2;i++){
	hF2_x[i]->DrawClone("hist ][  same");
	hF2_x_IC1[i]->DrawClone("PE X0 same");
	hF2_x_IC1[i]->SetLineStyle(7);
	hF2_x_IC1[i]->DrawClone("hist ][ same");
	scale[i] = TMath::Power(7,(i+1))/100.;
	lat.SetTextSize(0.03);
	if(i<3)lat.DrawLatex(latx[i]-0.7,laty[i]*200,Form("#color[4]{%1.1f[%i]}",qq[i],i+1));//,scale[i]));
	else if(i<5) lat.DrawLatex(latx[i]-0.7,laty[i]*200,Form("#color[4]{%1.f[%i]}",qq[i],i+1));//,scale[i]));
	else lat.DrawLatex(latx[i]-0.7,laty[i]*200,Form("#color[4]{%1.f[%i]}",qq[i],i+1));//,scale[i]));
	
    }
    TLegend *leg = new TLegend(0.23,0.75,0.4,0.92);
    leg->SetTextSize(0.032);
    
    leg->SetHeader("PYTHIA e+p, #it{L}=10 fb^{-1}");
    leg->AddEntry(gF2,"#sqrt{s}=63 GeV","PEL");
    leg->AddEntry(gF2_IC1,"#sqrt{s}=29 GeV","PEL");
    leg->Draw("same");
    lat.SetTextSize(0.032);
    lat.DrawLatex(-1.,0.03,"#color[4]{Q^{2} (GeV/#it{c}^{2})[i]}");
    lat.DrawLatex(-1.,0.008,"#color[4]{C=7^{i}}");

    gPad->SetLogy();
    

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
