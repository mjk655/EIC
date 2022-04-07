void setXbins(int nbins_x, double binning_x[]){
    for(int i = 0;i<nbins_x+1;i++)
        binning_x[i] = -4. + i * 4. / nbins_x;
}
void setQ2bins(int nbins_q2, double binning_q2[]){
    for(int i = 0;i<nbins_q2+1;i++)
        binning_q2[i] = i * 2.6 / nbins_q2;
}
void Plot_EPPS_F2(int save=1){
    gROOT->ProcessLine(".x ~/myStyle.C");
    int const nvary =10000;//                                                                                                                                                                                                                                              
    int const nbins_q2 = 13;
    int const nbins_x  = 20;
    double binning_q2[nbins_q2+1];
    double binning_x[nbins_x+1];
    setXbins(nbins_x,binning_x);
    setQ2bins(nbins_q2,binning_q2);
    TH2F* hF2_Replica[nbins_q2][nvary];    
    TH1F *hCS_x1[10000];
    TH1F *hCS_x2[10000];
    TH1F *hF2_x[10000];

    TH1F *hCS_x1_eAu[10000];
    TH1F *hCS_x2_eAu[10000];
    TH1F *hF2_x_eAu[10000];

    TH1F *hCS_x1_eAu_pseudo[10000];
    TH1F *hCS_x2_eAu_pseudo[10000];
    TH1F *hF2_x_eAu_pseudo[10000];

    TH1F *hCS_x1_eAu_er[10000];
    TH1F *hCS_x2_eAu_er[10000];
    TH1F *hF2_x_eAu_er[10000];
    TH1F *hCS_x1_eAu_ernew[10000];
    TH1F *hCS_x2_eAu_ernew[10000];
    TH1F *hF2_x_eAu_ernew[10000];
    TH2F* N_GLUON[10000];
    TH2F* N_GLUON_Weight[10000];
    TH2F* N_GLUON_Replica[10000];
    TH2F* N_GLUON_Comb[10000];
    TFile *fin = new TFile("results11_NewWeights_highstat.root","READ");
    TH1F *REP_WEIGHTS = (TH1F*)fin->Get("REP_WEIGHTS");
    for(int i = 1;i<nbins_q2-2;i++){
	hCS_x2_eAu_er[i]= (TH1F*)fin->Get(Form("hCS_x2_eAu_er_%i",i));
	hCS_x2_eAu_ernew[i]= (TH1F*)fin->Get(Form("hCS_x2_eAu_ernew_%i",i));
	hCS_x2_eAu[i]= (TH1F*)fin->Get(Form("hCS_x2_eAu_%i",i));
	hCS_x1_eAu_er[i]= (TH1F*)fin->Get(Form("hCS_x1_eAu_er_%i",i));
	hCS_x1_eAu_ernew[i]= (TH1F*)fin->Get(Form("hCS_x1_eAu_ernew_%i",i));
	hCS_x1_eAu[i]= (TH1F*)fin->Get(Form("hCS_x1_eAu_%i",i));
	hF2_x_eAu_er[i]= (TH1F*)fin->Get(Form("hF2_x_eAu_er_%i",i));
	hF2_x_eAu_ernew[i]= (TH1F*)fin->Get(Form("hF2_x_eAu_ernew_%i",i));
	hF2_x_eAu[i]= (TH1F*)fin->Get(Form("hF2_x_eAu_%i",i));
	hF2_x_eAu_pseudo[i]= (TH1F*)fin->Get(Form("hF2_x_eAu_%i_pseudo",i));
	hCS_x1_eAu_pseudo[i]= (TH1F*)fin->Get(Form("hCS_x1_eAu_%i_pseudo",i));
	hCS_x2_eAu_pseudo[i]= (TH1F*)fin->Get(Form("hCS_x2_eAu_%i_pseudo",i));
    }
    for(int i = 0;i<nvary;i++){
	N_GLUON_Weight[i]= (TH2F*)fin->Get(Form("N_GLUON_Weight_%i",i));
	N_GLUON_Replica[i]= (TH2F*)fin->Get(Form("N_GLUON_Replica_%i",i));
	N_GLUON[i]= (TH2F*)fin->Get(Form("N_GLUON_%i",i));
	N_GLUON_Comb[i]= (TH2F*)fin->Get(Form("N_GLUON_Comb_%i",i));
	for(int ii = 1;ii<nbins_q2-2;ii++)hF2_Replica[ii][i] = (TH2F*)fin->Get(Form("hF2_Replica_%i_%i",ii,i));
    }
    TH1F* ERRS[30];
    TH1F* ERRS_W[30];
    TH1F * TEMP = (TH1F*)N_GLUON[0]->ProjectionX();
    for(int i = 0;i<30;i++){
        ERRS[i] = (TH1F*)TEMP->Clone(Form("ERRS_%i",i));
        ERRS_W[i] = (TH1F*)TEMP->Clone(Form("ERRS_W_%i",i));
    }
    for(int i = 2;i<nbins_q2-2;i++){
	cout << "QBIN " << i << endl;
	for(int j = 1;j<nbins_x+1;j++){
	    double cen = hF2_x_eAu[i-1]->GetBinContent(j);
	    double sum = 0; 
	    double sum1 = 0;
	    TH1F *temp = new TH1F("temp","temp",10000,-5,5);
	    for(int r = 0;r<nvary;r++){
		double val = hF2_Replica[i-1][r]->GetBinContent(j);
		temp->Fill(cen-val);
		sum += REP_WEIGHTS->GetBinContent(r) * (cen-val)*(cen-val);
	       
		sum1 += (cen-val)*(cen-val);
		if(i==2 && j==11) cout << "DEBUG cen " << cen << " new " << val << " sum " << sum << " " << sum1 <<   " " << REP_WEIGHTS->GetBinContent(r) << endl; 
	    }
	    if(i==2 && j==11) cout << "DEBUG var " << sqrt(sum/nvary) << " " << sqrt(sum1/nvary) << endl;
	    hF2_x_eAu_ernew[i-1]->SetBinError(j,sqrt(sum/nvary));
	    delete temp;
	}
    }


    TF1 *line = new TF1("line","1",-100,100);
    line->SetLineStyle(7);
    TLatex lat;
    TLegend *leg = new TLegend(0.17,0.62,0.4,0.87);

    double qq[20] = {1.3,2,3.1,5,8,12,20,31,52,76,120};

    TCanvas *c11 = new TCanvas("c11","f2",1200,900);
    // Canvas setup                                                                                                                                                                                                                                                      
    const Int_t Nx1 = 3;
    const Int_t Ny1 = 3;
    CanvasPartition(c11,Nx1,Ny1,0.07,0.05,0.1,0.05);
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
    for(int i = 1;i < nbins_q2-3 ; i++){
        if(i<4)pad1[i-1][2]->cd();
        else if(i<7)pad1[i-4][1]->cd();
        else pad1[i-7][0]->cd();
        hF2_x_eAu_er[i]->SetFillColor(1);
        hF2_x_eAu_er[i]->SetFillStyle(3004);
        hF2_x_eAu_er[i]->SetMarkerSize(0);
        hF2_x_eAu_ernew[i]->SetFillColor(38);
        hF2_x_eAu_ernew[i]->SetLineColor(38);
        hF2_x_eAu_ernew[i]->SetFillStyle(1000);
        hF2_x_eAu_ernew[i]->SetMarkerSize(0);
        hF2_x_eAu_pseudo[i]->GetYaxis()->SetRangeUser(0.4,1.65);
        hF2_x_eAu_pseudo[i]->GetXaxis()->SetRangeUser(-3.4,-0.3);
        hF2_x_eAu_pseudo[i]->GetYaxis()->SetNdivisions(5);
        hF2_x_eAu_pseudo[i]->GetXaxis()->SetTitle("log_{10}(x_{B})");
        hF2_x_eAu_pseudo[i]->GetYaxis()->SetTitle("R_{eA}(F_{2}^{c#bar{c}})");
        hF2_x_eAu_pseudo[i]->GetXaxis()->CenterTitle(1);
        hF2_x_eAu_pseudo[i]->GetYaxis()->CenterTitle(1);
        hF2_x_eAu_pseudo[i]->GetYaxis()->SetTitleSize(0.1);
        hF2_x_eAu_pseudo[i]->GetXaxis()->SetTitleSize(0.1);
        if(i==3 || i==6 || i== 9)hF2_x_eAu_pseudo[i]->GetXaxis()->SetRangeUser(-3.4,0);
        if(i==1 || i==4){
            if(i==1)hF2_x_eAu_pseudo[i]->GetYaxis()->SetTitleSize(0.117);
            else hF2_x_eAu_pseudo[i]->GetYaxis()->SetTitleSize(0.14);
            if(i==1)hF2_x_eAu_pseudo[i]->GetYaxis()->SetTitleOffset(0.7);
            else hF2_x_eAu_pseudo[i]->GetYaxis()->SetTitleOffset(0.62);
            if(i==1)hF2_x_eAu_pseudo[i]->GetYaxis()->SetLabelSize(0.07);
            else hF2_x_eAu_pseudo[i]->GetYaxis()->SetLabelSize(0.083);

	}
        else{
            hF2_x_eAu_pseudo[i]->GetXaxis()->SetLabelSize(0.08);

            hF2_x_eAu_pseudo[i]->GetXaxis()->SetTitleOffset(0.75);
            hF2_x_eAu_pseudo[i]->GetYaxis()->SetTitleOffset(0.75);
        }
        if(i==1){
            hF2_x_eAu_pseudo[i]->SetBinContent(13,-100);
            hF2_x_eAu_er[i]->SetBinContent(13,-100);
            hF2_x_eAu_ernew[i]->SetBinContent(13,-100);
        }
        if(i==2){
            hF2_x_eAu_pseudo[i]->SetBinContent(15,-100);
            hF2_x_eAu_er[i]->SetBinContent(15,-100);
            hF2_x_eAu_ernew[i]->SetBinContent(15,-100);
        }
        if(i==3){
            hF2_x_eAu_pseudo[i]->SetBinContent(16,-100);
            hF2_x_eAu_er[i]->SetBinContent(16,-100);
            hF2_x_eAu_ernew[i]->SetBinContent(16,-100);
        }
        if(i==4){
            hF2_x_eAu_pseudo[i]->SetBinContent(17,-100);
            hF2_x_eAu_er[i]->SetBinContent(17,-100);
            hF2_x_eAu_ernew[i]->SetBinContent(17,-100);
        }
        if(i==6){
            hF2_x_eAu_er[i]->SetBinContent(18,-100);
            hF2_x_eAu_ernew[i]->SetBinContent(18,-100);
        }
	if(i==7){
            hF2_x_eAu_pseudo[i]->SetBinContent(19,-100);
            hF2_x_eAu_er[i]->SetBinContent(19,-100);
            hF2_x_eAu_ernew[i]->SetBinContent(19,-100);
        }
        if(i==8){
            hF2_x_eAu_pseudo[i]->SetBinContent(19,-100);
            hF2_x_eAu_er[i]->SetBinContent(19,-100);
            hF2_x_eAu_ernew[i]->SetBinContent(19,-100);
        }
        if(i==9){
            hF2_x_eAu_pseudo[i]->SetBinContent(18,-100);
            hF2_x_eAu_er[i]->SetBinContent(18,-100);
            hF2_x_eAu_ernew[i]->SetBinContent(18,-100);
            hF2_x_eAu_pseudo[i]->SetBinContent(19,-100);
            hF2_x_eAu_er[i]->SetBinContent(19,-100);
            hF2_x_eAu_ernew[i]->SetBinContent(19,-100);
        }
        if(i==10){
            hF2_x_eAu_pseudo[i]->SetBinContent(19,-100);
            hF2_x_eAu_er[i]->SetBinContent(19,-100);
            hF2_x_eAu_ernew[i]->SetBinContent(19,-100);
        }
        if(i==2){
            TLegend *canv = new TLegend(0.08,0.54,0.4,0.8);
            canv->SetTextSize(0.05);
            canv->AddEntry(hF2_x_eAu_pseudo[i],"EIC Pseudodata","Pl");
            canv->AddEntry(hF2_x_eAu_er[i],"EPPS16 Baseline","f");
            canv->AddEntry(hF2_x_eAu_ernew[i],"EPPS16 With F_{2}^{c#bar{c}}","f");
        }

        hF2_x_eAu_pseudo[i]->Draw("PE X0");
        hF2_x_eAu_er[i]->Draw("E2 same");
        hF2_x_eAu_ernew[i]->Draw("E2 same");
        hF2_x_eAu_pseudo[i]->Draw("PE X0 same");
        line->Draw("same L");
        lat.SetTextSize(0.08);
        if(i==1) lat.DrawLatex(-3.2,1.53,"PYTHIA e+p/Au #sqrt{s}=63,29 GeV" );
        //if(i==1) lat.DrawLatex(-3.2,1.4,"#it{L} = 1 fb^{-1}/nucleon" );                                                                                                                                                                                                
        lat.SetTextSize(0.067);
        if(i==1|| i==2 || i==3) lat.SetTextSize(0.075);
        if(i==4|| i==5 || i==6) lat.SetTextSize(0.09);

	lat.DrawLatex(-3,0.6,Form("Q^{2} = %1.f GeV^{2}",qq[i]));

	gPad->RedrawAxis();

	if(i==2)canv->Draw("same");
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
