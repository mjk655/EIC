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
void Plot(int save=0){
    gROOT->ProcessLine(".x ~/myStyle.C");
    int const nvary =50;//                                                                                                                                                                                                                                              
    int const nbins_q2 = 13;
    int const nbins_x  = 20;
    double binning_q2[nbins_q2+1];
    double binning_x[nbins_x+1];
    setXbins(nbins_x,binning_x);
    setQ2bins(nbins_q2,binning_q2);
    
    TH1F *hCS_x1[100];
    TH1F *hCS_x2[100];
    TH1F *hF2_x[100];
    TH1F *hCS_x1_eAu[100];
    TH1F *hCS_x2_eAu[100];
    TH1F *hF2_x_eAu[100];
    TH1F *hCS_x1_eAu_pseudo[100];
    TH1F *hCS_x2_eAu_pseudo[100];
    TH1F *hF2_x_eAu_pseudo[100];
    TH1F *hCS_x1_eAu_er[100];
    TH1F *hCS_x2_eAu_er[100];
    TH1F *hF2_x_eAu_er[100];
    TH1F *hCS_x1_eAu_ernew[100];
    TH1F *hCS_x2_eAu_ernew[100];
    TH1F *hF2_x_eAu_ernew[100];
    TH2F* N_GLUON[152];
    TH2F* N_GLUON_Weight[97];
    TH2F* N_GLUON_Comb[97];
    TFile *fin = new TFile("results1.root","READ");
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
    for(int i = 0;i<40/2;i++){
	N_GLUON_Weight[i]= (TH2F*)fin->Get(Form("N_GLUON_Weight_%i",i));
	N_GLUON[i]= (TH2F*)fin->Get(Form("N_GLUON_%i",i));
	N_GLUON_Comb[i]= (TH2F*)fin->Get(Form("N_GLUON_Comb_%i",i));
    }

    TH1F* ERRS[30];
    TH1F* ERRS_W[30];
    TH1F * TEMP = (TH1F*)N_GLUON[0]->ProjectionX();
    TH1F * TEMP1 = (TH1F*)N_GLUON_Comb[0]->ProjectionX();
    for(int i = 0;i<30;i++){
        ERRS[i] = (TH1F*)TEMP->Clone(Form("ERRS_%i",i));
        ERRS_W[i] = (TH1F*)TEMP1->Clone(Form("ERRS_W_%i",i));
    }


    double qq[20] = {1.3,2,3.1,5,8,12,20,31,52,76,120};
    TF1 *line = new TF1("line","1",-100,100);
    line->SetLineStyle(7);
    TLatex lat;
    TLegend *leg = new TLegend(0.17,0.62,0.4,0.87);
    
/*           
	     TCanvas *c1 = new TCanvas("c1","10-100 Ratios",2000,800);                                                                                                                                                                                                        
	     c1->Divide(5,2);                                                                                                                                                                                                                                   

    Double_t W = 0.2; // Pad WidthX                                                                                                                                                                                                                     
    Int_t Nx = 5; // Number of pads along X                                                                                                                                                                                                             
    Double_t Xm = (1-(Nx*W))/2; // X Margin                                                                                                                                                                                                             
    Double_t dw = (W*0.1)/4;

    Double_t Wy = 0.5; // Pad Width Y                                                                                                                                                                                                                   
    Int_t Ny = 2; // Number of pads along Y                                                                                                                                                                                                             
    Double_t Ym = (1-(Ny*Wy))/2; // y Margin                                                                                                                                                                                                           
    Double_t dwy = (Wy*0.1)/4;
    TCanvas *c111 = new TCanvas("c111","gluon ratios",1500,600);
    c111->Divide(5,2);
    TPad *p1 = new TPad("p1", "p1", Xm,Ym+Wy+dwy, Xm+W+dw, Ym+2*Wy-dwy, Ym+2*Wy-dwy, 0, 0);   p1->SetRightMargin(0);p1->SetBottomMargin(0);p1->Draw();
    TPad *p2 = new TPad("p2", "p2", Xm+W+dw, Ym+Wy+dwy, Xm+2*W-dw,  Ym+2*Wy-dwy, 0, 0, 0);  p2->SetRightMargin(0);p2->SetBottomMargin(0);p2->SetLeftMargin(0);p2->Draw();
    TPad *p3 = new TPad("p3", "p3", Xm+2*W+dw, Ym+Wy+dwy, Xm+3*W-dw,  Ym+2*Wy-dwy, 0, 0, 0);p3->SetRightMargin(0);p3->SetBottomMargin(0);p3->SetLeftMargin(0);p3->Draw();
    TPad *p4 = new TPad("p4", "p4", Xm+3*W+dw, Ym+Wy+dwy, Xm+4*W-dw,  Ym+2*Wy-dwy, 0, 0, 0);p4->SetRightMargin(0);p4->SetBottomMargin(0);p4->SetLeftMargin(0);p4->Draw();
    TPad *p5 = new TPad("p5", "p5", Xm+4*W-dw, Ym+Wy+dwy, Xm+5*W,  Ym+2*Wy-dwy, 0, 0, 0);   p5->SetLeftMargin(0);p5->SetBottomMargin(0);p5->Draw();

    TPad *p6 = new TPad("p6", "p6", Xm,Ym, Xm+W+dw, Ym+Wy+dwy, Ym+Wy+dwy, 0, 0);   p6->SetRightMargin(0);p6->SetTopMargin(0);p6->Draw();
    TPad *p7 = new TPad("p7", "p7", Xm+W+dw, Ym, Xm+2*W-dw,  Ym+Wy+dwy, 0, 0, 0);  p7->SetRightMargin(0);p7->SetTopMargin(0);p7->SetLeftMargin(0);p7->Draw();
    TPad *p8 = new TPad("p8", "p8", Xm+2*W+dw, Ym, Xm+3*W-dw,  Ym+Wy+dwy, 0, 0, 0);p8->SetRightMargin(0);p8->SetTopMargin(0);p8->SetLeftMargin(0);p8->Draw();
    TPad *p9 = new TPad("p9", "p9", Xm+3*W+dw, Ym, Xm+4*W-dw,  Ym+Wy+dwy, 0, 0, 0);p9->SetRightMargin(0);p9->SetTopMargin(0);p9->SetLeftMargin(0);p9->Draw();
    TPad *p10 = new TPad("p10", "p10", Xm+4*W-dw, Ym, Xm+5*W,  Ym+Wy+dwy, 0, 0, 0);p10->SetLeftMargin(0);p10->SetTopMargin(0);p10->Draw();
              
    for(int i = 1;i<nbins_q2-2;i++){                                                                                                                                                                                                                                 
        c1->cd(i);                                                                                                                                                                                                                                                   
        hCS_x1_eAu_er[i]->SetFillColor(1);
        hCS_x1_eAu_er[i]->SetFillStyle(3004);
        hCS_x1_eAu_er[i]->SetMarkerSize(0);
        hCS_x1_eAu_ernew[i]->SetFillColor(38);
        hCS_x1_eAu_ernew[i]->SetFillStyle(1000);
        hCS_x1_eAu_ernew[i]->SetMarkerSize(0);                                                                                                                                                                                                    
        hCS_x1_eAu_pseudo[i]->GetXaxis()->SetRangeUser(-3.5,0);                                                                                                                                                                                 
        hCS_x1_eAu_pseudo[i]->GetYaxis()->SetRangeUser(0.5,2);                                                                                                                                                                                      
        hCS_x1_eAu_pseudo[i]->GetYaxis()->SetNdivisions(5);                                                                                                                                                                                         
        hCS_x1_eAu_pseudo[i]->GetXaxis()->SetTitle("log(x_{B})");                                                                                                                                                                                  
        hCS_x1_eAu_pseudo[i]->GetYaxis()->SetTitle("#sigma^{c#bar{c}}_{r}(e+Au)/#sigma^{c#bar{c}}_{r}(e+p)");                                                                                                                                  
        hCS_x1_eAu_pseudo[i]->Draw("PE X0 same");                                                                                                                                                                                                  
        hCS_x1_eAu_er[i]->Draw("E2 same");                                                                                                                                                                                                         
        hCS_x1_eAu_ernew[i]->Draw("E2 same");                                                                                                                                                                                                   
        hCS_x1_eAu_pseudo[i]->Draw("PE X0 same");                                                                                                                                                                                                    
        lat.DrawLatex(-3.2,1.7,"e+p/Au 10#times100 GeV");                                                                                                                                                                                           
        lat.DrawLatex(-3.2,1.5,Form("Q^{2} = %1.f GeV^{2}",qq[i]));                                                                                                                                                                               
        line->Draw("same");                                                                                                                                                                                                                                          
    }                                                                                                                                                                                                                                                                
    
    TCanvas *c2 = new TCanvas("c2","5-41 Ratios",2000,800);                                                                                                                                                                                                          
    c2->Divide(5,2);                                                                                                                                                                                                                                                 
    for(int i = 1;i<nbins_q2-2;i++){                                                                                                                                                                                                                                 
        c2->cd(i);                                                                                                                                                                                                                                                   
        hCS_x2_eAu_er[i]->SetFillColor(1);                                                                                                                                                                                                                           
        hCS_x2_eAu_er[i]->SetFillStyle(3004);                                                                                                                                                                                                                        
        hCS_x2_eAu_er[i]->SetMarkerSize(0);                                                                                                                                                                                                                          
        hCS_x2_eAu_ernew[i]->SetFillColor(38);                                                                                                                                                                                                                       
        hCS_x2_eAu_ernew[i]->SetFillStyle(1000);                                                                                                                                                                                                                     
        hCS_x2_eAu_ernew[i]->SetMarkerSize(0);                                                                                                                                                                                                                       
        hCS_x2_eAu_pseudo[i]->GetXaxis()->SetRangeUser(-3.5,0);                                                                                                                                                                                                      
        hCS_x2_eAu_pseudo[i]->GetYaxis()->SetRangeUser(0.5,2);                                                                                                                                                                                                     
        hCS_x2_eAu_pseudo[i]->GetYaxis()->SetNdivisions(5);                                                                                                                                                                                                          
        hCS_x2_eAu_pseudo[i]->GetXaxis()->SetTitle("log(x_{B})");                                                                                                                                                                                                    
        hCS_x2_eAu_pseudo[i]->GetYaxis()->SetTitle("#sigma^{c#bar{c}}_{r}(e+Au)/#sigma^{c#bar{c}}_{r}(e+p)");                                                                                                                                                        
        hCS_x2_eAu_pseudo[i]->Draw("PE X0 same");                                                                                                                                                                                                                    
        hCS_x2_eAu_er[i]->Draw("E2 same");
        hCS_x2_eAu_ernew[i]->Draw("E2 same");                                                                                                                                                                                                                        
        hCS_x2_eAu_pseudo[i]->Draw("PE X0 same");                                                                                                                                                                                                                    
        lat.DrawLatex(-3.2,1.7,"e+p/Au 5#times41 GeV");                                                                                                                                                                                                              
        lat.DrawLatex(-3.2,1.5,Form("Q^{2} = %1.f GeV^{2}",qq[i]));                                                                                                                                                                                                  
        line->Draw("same");                                                                                                                                                                                                                     
	}                                                                                                                                                                                                                                                    */            

/*


    TCanvas *c11 = new TCanvas("c111","f2",1500,700);
    // Canvas setup                                                                                                                                                                                                                                    
    const Int_t Nx1 = 5;
    const Int_t Ny1 = 2;
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


    for(int i = 1;i < nbins_q2-2 ; i++){                                                                                                                                                                                                                     if(i<6)pad1[i-1][1]->cd();
	else pad1[i-1-5][0]->cd();
	
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
        hF2_x_eAu_pseudo[i]->GetXaxis()->SetTitle("log(x_{B})");
        hF2_x_eAu_pseudo[i]->GetYaxis()->SetTitle("R_{eA}(F_{2}^{c#bar{c}})");
	hF2_x_eAu_pseudo[i]->GetXaxis()->CenterTitle(1);
	hF2_x_eAu_pseudo[i]->GetYaxis()->CenterTitle(1);
	hF2_x_eAu_pseudo[i]->GetYaxis()->SetTitleSize(0.1);
	hF2_x_eAu_pseudo[i]->GetXaxis()->SetTitleSize(0.1);
	
	if(i==6 || i==10){
	    hF2_x_eAu_pseudo[i]->GetYaxis()->SetTitleSize(0.087);
	    hF2_x_eAu_pseudo[i]->GetXaxis()->SetTitleSize(0.08);

	}
	else{
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
	if(i==10){                                                                                                                                                                                                                             
            TLegend *canv = new TLegend(0.1,0.65,0.4,0.9);                                                                                                                                                                                         
            canv->SetTextSize(0.05);                                                                                                                                                                       
	    canv->AddEntry(hF2_x_eAu_pseudo[i],"Pseudo-data","Pl");
            canv->AddEntry(hF2_x_eAu_er[i],"Baseline","f");                                                                                                                                                                                   
            canv->AddEntry(hF2_x_eAu_ernew[i],"With F_{2}^{c#bar{c}}","f");                                                                                                                                                                     
        }    

        hF2_x_eAu_pseudo[i]->Draw("PE X0");
        hF2_x_eAu_er[i]->Draw("E2 same");
        hF2_x_eAu_ernew[i]->Draw("E2 same");
        hF2_x_eAu_pseudo[i]->Draw("PE X0 same");
        line->Draw("same L");
	lat.SetTextSize(0.06);
	if(i==1) lat.DrawLatex(-3.2,1.53,"PYTHIA e+p/Au #sqrt{s}=63,29 GeV" );
	
	lat.SetTextSize(0.08);
	if(i==6|| i==10 || i==1 || i==5) lat.SetTextSize(0.065);
	lat.DrawLatex(-3,0.6,Form("Q^{2} = %1.f GeV^{2}",qq[i]));
	if(i==10)canv->Draw("same");	
	gPad->RedrawAxis();     
    }

*/

    TCanvas *c111 = new TCanvas("c111","gluon ratios",1300,600);
    // Canvas setup
    const Int_t Nx = 3;
    const Int_t Ny = 2;
    CanvasPartition(c111,Nx,Ny,0.09,0.05,0.12,0.05);
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

    

    for(int i = 3;i > 0; i+=-1){
	c111->cd(i);
	pad[i-1][1]->cd();
	int bin = 1;
	if(i==1)bin = 2;
	if(i==2)bin = 7;
	if(i==3)bin = 11;
	
	TH1F* True = (TH1F*)N_GLUON[0]->ProjectionX("True",bin,bin);
	TH1F* Mod = (TH1F*)N_GLUON_Weight[0]->ProjectionX("Mod",bin,bin);
	TH1F* Ini = (TH1F*)N_GLUON_Comb[0]->ProjectionX("Ini",bin,bin);

	//True->SetLineStyle(7);
	True->GetYaxis()->SetRangeUser(0.1,2);
	if(i==1)  True->GetXaxis()->SetRangeUser(-4,-0.01);
	if(i==2)  True->GetXaxis()->SetRangeUser(-3.99,-0.01);
	if(i==3)  True->GetXaxis()->SetRangeUser(-3.99,0);
	True->GetYaxis()->SetTitle("R_{g}^{Au} (x,Q^{2})");
	True->GetXaxis()->SetTitle("log_{10}(x)");
	True->GetYaxis()->SetTitleSize(0.14);
	True->GetYaxis()->SetTitleOffset(0.7);
	True->GetYaxis()->SetLabelSize(0.08);
	True->GetYaxis()->CenterTitle(1);
        True->GetXaxis()->CenterTitle(1);
	True->GetYaxis()->SetNdivisions(5);
	Ini->SetFillColor(1);
	Ini->SetFillStyle(3004);
	Mod->SetFillColor(38);
	Mod->SetFillStyle(1000);
	Ini->SetMarkerSize(0);
	Mod->SetMarkerSize(0);
	Mod->SetLineColor(38);

	if(i==1){
	    TLegend *canv = new TLegend(0.3,0.57,0.5,0.87);
	    canv->SetTextSize(0.07);
	    canv->AddEntry(Ini,"Baseline","f");
	    canv->AddEntry(Mod,"With F_{2}^{c#bar{c}}","f");
	}
	
	for(int ii = 1; ii< Ini->GetNbinsX()+1; ii++){
	    double sum1 = True->GetBinContent(ii);
	    double sum2 = True->GetBinContent(ii);
	    double sum21 = 0;
	    double sum22 = 0;
	    for(int j = 0; j< 40/2; j++){
		TH1F *_temp2 = (TH1F*)N_GLUON_Weight[j]->ProjectionX("_temp2",bin,bin);
		TH1F *_temp1 = (TH1F*)N_GLUON_Comb[j]->ProjectionX("_temp1",bin,bin);
		sum22+=(sum2 - _temp2->GetBinContent(ii))*(sum2 - _temp2->GetBinContent(ii));
		sum21+=(sum1 - _temp1->GetBinContent(ii))*(sum1 - _temp1->GetBinContent(ii));
		delete _temp1; delete _temp2;
	    }
	    Mod->SetBinContent(ii,sum2);
	    Ini->SetBinContent(ii,sum1);
	    Mod->SetBinError(ii,sqrt(sum22));
	    Ini->SetBinError(ii,sqrt(sum21));
	    if(sum2>0){
		ERRS_W[bin]->SetBinContent(ii,sqrt(sum22)/sum2);
		ERRS[bin]->SetBinContent(ii,sqrt(sum21)/sum1);
	    }else{
		ERRS_W[bin]->SetBinContent(ii,0);
                ERRS[bin]->SetBinContent(ii,0);
	    }
	}
	True->DrawClone("hist same ][");
	Ini->DrawClone("E2 same");
	Mod->DrawClone("E2 same");
	True->DrawClone("hist same ][");
	line->Draw("same L");
	if(i==1)canv->Draw("same");
	lat.SetTextSize(0.08);
	lat.DrawLatex(-1.7,1.7,Form("Q^{2} = %1.f GeV^{2}",qq[bin-1]));
	gPad->RedrawAxis();
	pad[i-1][0]->cd();
	if(i==1)ERRS[bin]->GetXaxis()->SetRangeUser(-4,-0.01);
	if(i==2)ERRS[bin]->GetXaxis()->SetRangeUser(-3.99,-0.01);
	if(i==3)ERRS[bin]->GetXaxis()->SetRangeUser(-3.99,0);
	ERRS[bin]->GetYaxis()->SetTitle("Rel. Uncertainty");                                                                                                                                                                           
        ERRS[bin]->GetYaxis()->SetRangeUser(0.0,1);                                                                                                                                                                                                
        ERRS[bin]->GetXaxis()->SetTitle("log_{10}(#it{x})");
	ERRS[bin]->GetYaxis()->SetLabelSize(0.08);
	ERRS[bin]->GetXaxis()->SetLabelSize(0.08);
	ERRS[bin]->GetXaxis()->SetTitleSize(0.1);
	ERRS[bin]->GetYaxis()->SetTitleSize(0.1);
	ERRS[bin]->GetXaxis()->CenterTitle(1);
	ERRS[bin]->GetYaxis()->CenterTitle(1);
	ERRS_W[bin]->SetLineColor(38);                                                                                                                                                                                                          
        ERRS[bin]->SetLineStyle(7);
	ERRS[bin]->DrawClone("hist ][");
        ERRS_W[bin]->DrawClone("hist same ][");
	if(i==3){
            TLegend *canv = new TLegend(0.45,0.69,0.65,0.94);
            canv->SetTextSize(0.07);
	    canv->AddEntry(ERRS[bin],"Baseline","l");
            canv->AddEntry(Mod,"With F_{2}^{c#bar{c}}","l");
        }
	if(i==3)canv->Draw("same");
	//gPad->SetLogy();
	gPad->RedrawAxis();
    }
  
/*
    TCanvas *c1111 = new TCanvas("c1111","gluon ratios errors",2000,800);
    c1111->Divide(5,2);
    for(int i = 1;i < nbins_q2-2; i++){
	c1111->cd(i);
	ERRS[i]->GetYaxis()->SetTitle("Relative Gluon nPDF Uncertainty");
	ERRS[i]->GetYaxis()->SetRangeUser(0.01,10);
	ERRS[i]->GetXaxis()->SetTitle("log(x)");
	ERRS_W[i]->SetLineColor(38);
	ERRS[i]->SetLineStyle(7);
	ERRS[i]->DrawClone("hist");
	ERRS_W[i]->DrawClone("hist same");
	lat.DrawLatex(-3.7,1.7,Form("Q^{2} = %1.f GeV^{2}",qq[i])); 
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
