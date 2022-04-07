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
void GetErrorRatio(int save=0){
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


    TFile *fin1 = new TFile("ep_output_new.root","READ");
    for(int i = 0;i<nbins_q2-2;i++){
	hF2_x_IC1[i]= (TH1F*)fin1->Get(Form("ep_CharmF2_Q2Bin_%i",i));
	hF2_x_IC1[i]->SetName(Form("hF2_x_IC1_%i",i));
    }
    TFile *fin3 = new TFile("ep_output_ATHENA1.root","READ");
    for(int i = 0;i<nbins_q2-2;i++){//nodifying the index                                                                                                                                                                                            
	hF2_x[i]= (TH1F*)fin3->Get(Form("ep_CharmF2_Q2Bin_%i",i));
    }
    cout<< "HERE " << endl;
    TH2F *hR = new TH2F("hR","hR",nbins_x,binning_x,nbins_q2,binning_q2);
    
    for(int i = 1; i<hR->GetNbinsX()+1;i++){
	for(int j =1; j<nbins_q2-2;j++){
	    cout << i << " " << j << endl;
	    double er1 = hF2_x[j-1]->GetBinError(i);
	    double er2 = hF2_x_IC1[j-1]->GetBinError(i);
	    if(er1>0 && er2>0) hR->SetBinContent(i,j,er1/er2);
	    else hR->SetBinContent(i,j,-1);

	}   
    }
//    hR->GetZaxis()->SetTitle("Uncer. Ratio (ATHENA/Det. Mat.)");
    hR->GetZaxis()->SetRangeUser(0.2,1.2);
    hR->GetZaxis()->SetTitleOffset(0.25);
    hR->GetXaxis()->SetRangeUser(-3,0);
    hR->GetYaxis()->SetRangeUser(0,2.6);
    hR->GetXaxis()->SetTitle("x_{B}");
    hR->GetYaxis()->SetTitle("log_{10}(Q^{2} [GeV^{2}])");

    TCanvas *c1 = new TCanvas("c1","c1");
    hR->Draw("COLZ");
    gStyle->SetPadRightMargin(0.5);


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
