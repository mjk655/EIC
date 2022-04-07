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
void PlotSepGluon_Alt(int save=0){
    gROOT->ProcessLine(".x ~/myStyle.C");
    
    int const nbins_q2 = 13;
    int const nbins_x  = 20;
    double binning_q2[nbins_q2+1];
    double binning_x[nbins_x+1];
    setXbins(nbins_x,binning_x);
    setQ2bins(nbins_q2,binning_q2);
    double qq[20] = {1.3,2,3.1,5,8,12,20,31,52,76,120};

    TH1D* Base_EPPS[30];
    TH1D* New_EPPS[30];
    TH1D* ERRS_EPPS[30];
    TH1D* ERRS_New_EPPS[30];    
    TFile *fin = new TFile("EPPS_Gluon_NewWeights_noT.root","READ");
    for(int i = 1;i<nbins_q2-2;i++){
	ERRS_New_EPPS[i] = (TH1D*)fin->Get(Form("ERRS_New_EPPS_%i",i));//naming error
	ERRS_EPPS[i] =(TH1D*)fin->Get(Form("ERRS_EPPS_%i",i));
	New_EPPS[i]=(TH1D*)fin->Get(Form("New_EPPS_%i",i));
	Base_EPPS[i]=(TH1D*)fin->Get(Form("Base_EPPS_%i",i));
    }


    TH1D* Base_CTEQ[30];
    TH1D* New_CTEQ[30];
    TH1D* ERRS_CTEQ[30];
    TH1D* ERRS_New_CTEQ[30];
    TFile *fin1 = new TFile("CTEQ_Gluon_NewWeights_noT.root","READ");
    for(int i = 1;i<nbins_q2-2;i++){
        ERRS_New_CTEQ[i] = (TH1D*)fin1->Get(Form("ERRS_New_CTEQ_%i",i));
        ERRS_CTEQ[i] =(TH1D*)fin1->Get(Form("ERRS_CTEQ_%i",i));
        New_CTEQ[i]=(TH1D*)fin1->Get(Form("New_CTEQ_%i",i));
        Base_CTEQ[i]=(TH1D*)fin1->Get(Form("Base_CTEQ_%i",i));
    }

    TH1D* Base_NNPDF[30];
    TH1D* New_NNPDF[30];
    TH1D* ERRS_NNPDF[30];
    TH1D* ERRS_New_NNPDF[30];
    TFile *fin2 = new TFile("NNPDF_Gluon_NewWeights.root","READ");
    for(int i = 1;i<nbins_q2-2;i++){
        ERRS_New_NNPDF[i] = (TH1D*)fin2->Get(Form("ERRS_New_NNPDF_%i",i));
        ERRS_NNPDF[i] =(TH1D*)fin2->Get(Form("ERRS_NNPDF_%i",i));
        New_NNPDF[i]=(TH1D*)fin2->Get(Form("New_NNPDF_%i",i));
        Base_NNPDF[i]=(TH1D*)fin2->Get(Form("Base_NNPDF_%i",i));
    }

    TCanvas *c111 = new TCanvas("c111","gluon ratios",1300,600);
    // Canvas setup
    const Int_t Nx = 3;
    const Int_t Ny = 2;
    CanvasPartition(c111,Nx,Ny,0.09,0.05,0.15,0.05);
    TPad *pad[Nx][Ny];
    for (Int_t i=0;i<Nx;i++) {
	for (Int_t j=0;j<Ny;j++) {
	    c111->cd(0);
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
    TF1 *line = new TF1("line","1",-100,100);
    line->SetLineStyle(1);
    //line->SetLineWidth(0.5);
    
    TF1 *line1 = new TF1("line1","5",-100,100);
    line1->SetLineStyle(1);
    TF1 *line11 = new TF1("line11","10",-100,100);
    line11->SetLineStyle(1);
    TF1 *line2 = new TF1("line2","0.01",-100,100);
    line2->SetLineStyle(1);
    line2->SetLineWidth(0.5);




    TLatex lat;
    for(int i = 3;i > 0; i+=-1){
	c111->cd(i);
	pad[i-1][1]->cd();
	int bin = 11-1;//2 / 7 / or 11
	if(i==1)  Base_EPPS[bin]->GetXaxis()->SetRangeUser(-4,-0.01);
	if(i==2)  Base_EPPS[bin]->GetXaxis()->SetRangeUser(-3.99,-0.01);
	if(i==3)  Base_EPPS[bin]->GetXaxis()->SetRangeUser(-3.99,0);
	Base_EPPS[bin]->GetYaxis()->SetRangeUser(0.1,2.1);
	Base_EPPS[bin]->GetYaxis()->SetTitle("R_{g}^{Au} (x,Q^{2})");
	Base_EPPS[bin]->GetXaxis()->SetTitle("log_{10}(x_{g})");
	Base_EPPS[bin]->GetYaxis()->SetTitleSize(0.13);
	Base_EPPS[bin]->GetYaxis()->SetTitleOffset(0.7);
	Base_EPPS[bin]->GetYaxis()->SetLabelSize(0.09);
	Base_EPPS[bin]->GetYaxis()->CenterTitle(1);
        Base_EPPS[bin]->GetXaxis()->CenterTitle(1);
	Base_EPPS[bin]->GetYaxis()->SetNdivisions(5);
	
	New_EPPS[bin]->SetFillColor(15);
	Base_EPPS[bin]->SetLineColor(15);
	New_EPPS[bin]->SetFillStyle(1000);
	New_EPPS[bin]->SetLineStyle(1);
	Base_EPPS[bin]->SetMarkerSize(0);
	New_EPPS[bin]->SetMarkerSize(0);
	New_EPPS[bin]->SetLineColor(12);
	TH1F* temp1 = (TH1F*)Base_EPPS[bin]->Clone("temp1");
        temp1->SetFillColor(1);
        temp1->SetFillStyle(3006);
        temp1->SetMarkerSize(0);
	TH1F* temp11 = (TH1F*)New_EPPS[bin]->Clone("temp1");
        temp11->SetFillColor(15);
        temp11->SetFillStyle(1000);
        temp11->SetMarkerSize(0);

	
	New_CTEQ[bin]->SetFillColor(38);
	Base_CTEQ[bin]->SetLineColor(38);
        New_CTEQ[bin]->SetFillStyle(1000);
	New_CTEQ[bin]->SetLineStyle(7);
        Base_CTEQ[bin]->SetMarkerSize(0);
        New_CTEQ[bin]->SetMarkerSize(0);
        New_CTEQ[bin]->SetLineColor(kBlue+1);
        TH1F* temp2 = (TH1F*)Base_CTEQ[bin]->Clone("temp2");
        temp2->SetFillColor(kBlue);
        temp2->SetFillStyle(3013);
        temp2->SetMarkerSize(0);
	
        New_NNPDF[bin]->SetFillColor(45);
        Base_NNPDF[bin]->SetLineColor(46);
	New_NNPDF[bin]->SetFillStyle(1000);
	New_NNPDF[bin]->SetLineStyle(5);
        Base_NNPDF[bin]->SetMarkerSize(0);
        New_NNPDF[bin]->SetMarkerSize(0);
        New_NNPDF[bin]->SetLineColor(kRed+1);
        TH1F* temp3 = (TH1F*)Base_NNPDF[bin]->Clone("temp3");
        temp3->SetFillColor(kRed);
        temp3->SetFillStyle(3007);
        temp3->SetMarkerSize(0);

	
	temp1->SetLineColor(0);
	temp2->SetLineColor(0);
	temp3->SetLineColor(0);
	if(i==1){
	    TLegend *canv = new TLegend(0.31,0.57,0.54,0.87);
	    canv->SetTextSize(0.055);
	    canv->SetHeader("EPPS16");
	    canv->AddEntry(temp1,"Baseline","f");
	    canv->AddEntry(temp11,"With F_{2}^{c#bar{c}}","fl");
	}
	if(i==2){
            TLegend *canv = new TLegend(0.1,0.57,0.39,0.87);
            canv->SetTextSize(0.055);
            canv->SetHeader("nCTEQ15");
            canv->AddEntry(temp2,"Baseline","f");
            canv->AddEntry(New_CTEQ[bin],"With F_{2}^{c#bar{c}}","fl");
        }
        if(i==3){
            TLegend *canv = new TLegend(0.1,0.57,0.39,0.87);
            canv->SetTextSize(0.055);
            canv->SetHeader("nNNPDF2.0");
            canv->AddEntry(temp3,"Baseline","f");
            canv->AddEntry(New_NNPDF[bin],"With F_{2}^{c#bar{c}}","fl");
        }
	
	if(i!=1){
	    Base_EPPS[bin]->SetLineColor(0);
	    Base_EPPS[bin]->SetFillColor(0);
	}
	Base_EPPS[bin]->DrawClone("hist same ][");
	if(i==1)temp1->DrawClone("E2 same");
	if(i==2)temp2->DrawClone("E2 same");
	if(i==3)temp3->DrawClone("E2 same");
	
	if(i==1)Base_EPPS[bin]->DrawClone("hist same ][");
	if(i==1)New_EPPS[bin]->DrawClone("E2 same");
	if(i==1)New_EPPS[bin]->SetFillColor(0);
	if(i==1)New_EPPS[bin]->DrawClone("hist same ][");
	if(i==2)Base_CTEQ[bin]->DrawClone("hist same ][");
        if(i==2)New_CTEQ[bin]->DrawClone("E2 same");
	if(i==2)New_CTEQ[bin]->SetFillColor(0);
	if(i==2)New_CTEQ[bin]->DrawClone("Hist same ][");
	if(i==3)Base_NNPDF[bin]->DrawClone("hist same ][");
        if(i==3)New_NNPDF[bin]->DrawClone("E2 same");
	if(i==3)New_NNPDF[bin]->SetFillColor(0);
	if(i==3)New_NNPDF[bin]->DrawClone("hist same ][");
	if(i==3)New_EPPS[bin]->SetFillColor(15);
	if(i==3)New_CTEQ[bin]->SetFillColor(38);
	if(i==3)New_NNPDF[bin]->SetFillColor(45);


	line->Draw("same L");
	canv->Draw("same");
	lat.SetTextSize(0.07);
	lat.DrawLatex(-1.6,1.9,Form("Q^{2} = %1.f GeV^{2}",qq[bin]));
	gPad->RedrawAxis();
	pad[i-1][0]->cd();
	if(i==1)ERRS_EPPS[bin]->GetXaxis()->SetRangeUser(-4,-0.01);
	if(i==2)ERRS_EPPS[bin]->GetXaxis()->SetRangeUser(-3.99,-0.01);
	if(i==3)ERRS_EPPS[bin]->GetXaxis()->SetRangeUser(-3.99,0);
	ERRS_EPPS[bin]->GetYaxis()->SetTitle("Red. Fact.");                                                                                                                                                                           
        ERRS_EPPS[bin]->GetYaxis()->SetRangeUser(0,13);
	ERRS_EPPS[bin]->GetYaxis()->SetNdivisions(5);
        ERRS_EPPS[bin]->GetXaxis()->SetTitle("log_{10}(x_{g})");
	ERRS_EPPS[bin]->GetYaxis()->SetLabelSize(0.13);
	ERRS_EPPS[bin]->GetXaxis()->SetLabelSize(0.15);
	ERRS_EPPS[bin]->GetXaxis()->SetTitleSize(0.2);
	ERRS_EPPS[bin]->GetYaxis()->SetTitleSize(0.15);
	ERRS_EPPS[bin]->GetXaxis()->CenterTitle(1);
	ERRS_EPPS[bin]->GetYaxis()->CenterTitle(1);
	
	ERRS_EPPS[bin]->GetXaxis()->SetTitleOffset(0.8);
	ERRS_EPPS[bin]->GetYaxis()->SetTitleOffset(0.6);

	ERRS_New_EPPS[bin]->SetMarkerColor(15);                                                                                                                                                                                                       
        ERRS_New_EPPS[bin]->SetMarkerStyle(24);
	ERRS_EPPS[bin]->SetLineColor(15);
        ERRS_EPPS[bin]->SetLineStyle(1);

	ERRS_New_CTEQ[bin]->SetMarkerColor(38);
	ERRS_New_CTEQ[bin]->SetMarkerStyle(27);
	ERRS_CTEQ[bin]->SetLineColor(38);
	ERRS_CTEQ[bin]->SetLineStyle(7);

	ERRS_New_NNPDF[bin]->SetMarkerColor(46);
        ERRS_New_NNPDF[bin]->SetMarkerStyle(30);
        ERRS_NNPDF[bin]->SetLineColor(46);
        ERRS_NNPDF[bin]->SetLineStyle(5);

	ERRS_New_EPPS[bin]->SetMarkerSize(0.6);
	ERRS_New_CTEQ[bin]->SetMarkerSize(0.6);
	ERRS_New_NNPDF[bin]->SetMarkerSize(0.6);
			
	if(i!=1)ERRS_EPPS[bin]->SetLineColor(0);
	ERRS_EPPS[bin]->DrawClone("hist ][");
        ERRS_EPPS[bin]->SetLineWidth(3);
	ERRS_CTEQ[bin]->SetLineWidth(3);
	ERRS_NNPDF[bin]->SetLineWidth(3);
	if(i==1){
	    ERRS_EPPS[bin]->Divide(ERRS_New_EPPS[bin]);
	    ERRS_EPPS[bin]->DrawClone("hist ][");
	}
	if(i==2){
	    ERRS_CTEQ[bin]->Divide(ERRS_New_CTEQ[bin]);
	    ERRS_CTEQ[bin]->DrawClone("same hist ][");
        }
	if(i==3){
	    ERRS_NNPDF[bin]->Divide(ERRS_New_NNPDF[bin]);
	    ERRS_NNPDF[bin]->DrawClone("same hist ][");
        }
	ERRS_New_EPPS[bin]->SetMarkerSize(1.5);
        ERRS_New_CTEQ[bin]->SetMarkerSize(1.5);
        ERRS_New_NNPDF[bin]->SetMarkerSize(1.5);
	if(i==1){
            if(bin==1)TLegend *canv1 = new TLegend(0.3,0.33,0.6,0.58);
	    else TLegend *canv1 = new TLegend(0.58,0.7,0.88,0.95);
            canv1->SetTextSize(0.055);
            canv1->SetHeader("EPPS16");
            canv1->AddEntry(ERRS_EPPS[bin],"Baseline","l");
            canv1->AddEntry(ERRS_New_EPPS[bin],"With F_{2}^{c#bar{c}}","p");
       }
	if(i==2){
            if(bin!=1)TLegend *canv1 = new TLegend(0.45,0.7,0.7,0.95);
            else TLegend *canv1 = new TLegend(0.09,0.33,0.34,0.58);
	    canv1->SetTextSize(0.055);
	    canv1->SetHeader("nCTEQ15");
	    canv1->AddEntry(ERRS_CTEQ[bin],"Baseline","l");
            canv1->AddEntry(ERRS_New_CTEQ[bin],"With F_{2}^{c#bar{c}}","p");
        }
        if(i==3){
            TLegend *canv1 = new TLegend(0.45,0.7,0.7,0.95);
            canv1->SetTextSize(0.055);
	    canv1->SetHeader("nNNPDF2.0");
            canv1->AddEntry(ERRS_NNPDF[bin],"Baseline","l");
            canv1->AddEntry(ERRS_New_NNPDF[bin],"With F_{2}^{c#bar{c}}","p");
        }
//	canv1->Draw("same");
	line1->Draw("same");
	line11->Draw("same");
//	line2->Draw("same");
//	gPad->SetLogy();
	gPad->RedrawAxis();
    }
    if(save){
        c111->SaveAs("paper/GluonRatios.pdf");
        c111->SaveAs("paper/GluonRatios.C");
        c111->SaveAs("paper/GluonRatios_2GeV.ps");
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
		vposu = 0.35;
		vfactor = vposu-vposd;
		vmard = bMargin / vfactor;
		vmaru = 0.0;
	    } else if (j == Ny-1) {
		vposd = 0.35;
		vposu = 1 - tMargin;
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
