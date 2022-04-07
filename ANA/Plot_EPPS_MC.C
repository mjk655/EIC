void setXbins(int nbins_x, double binning_x[]){
    for(int i = 0;i<nbins_x+1;i++)
        binning_x[i] = -4. + i * 4. / nbins_x;
}
void setQ2bins(int nbins_q2, double binning_q2[]){
    for(int i = 0;i<nbins_q2+1;i++)
        binning_q2[i] = i * 2.6 / nbins_q2;
}
void Plot_EPPS_MC(int save=0){
    gROOT->ProcessLine(".x ~/myStyle.C");
    int const nvary =1000;//                                                                                                                                                                                                                                              
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
    TFile *fin = new TFile("results11_NewWeights_Sys_Chi2.root","READ");
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
	for(int j = 1;j<nbins_x+1;j++){
	    double cen = hF2_x_eAu[i-1]->GetBinContent(j);
	    double sum = 0; 
	    double sum1 = 0;
	    TH1F *temp = new TH1F("temp","temp",10000,-5,5);
	    for(int r = 0;r<nvary;r++){
		double val = hF2_Replica[i-1][r]->GetBinContent(j);
		temp->Fill(cen-val);
		sum += REP_WEIGHTS->GetBinContent(r+1) * (cen-val)*(cen-val);
	       
		sum1 += (cen-val)*(cen-val);
		if(i==2 && j==11) cout << "DEBUG cen " << cen << " new " << val << " sum " << sum << " " << sum1 <<   " " << REP_WEIGHTS->GetBinContent(r) << endl; 
	    }
	    if(i==2 && j==11) cout << "DEBUG var " << sqrt(sum/nvary) << " " << sqrt(sum1/nvary) << endl;
	    hF2_x_eAu_ernew[i-1]->SetBinError(j,sqrt(sum/nvary));
	    delete temp;
	}
    }




    double qq[20] = {1.3,2,3.1,5,8,12,20,31,52,76,120};
    TF1 *line = new TF1("line","1",-100,100);
    line->SetLineStyle(7);
    TLatex lat;
    TLegend *leg = new TLegend(0.17,0.62,0.4,0.87);

    /*
    TCanvas *c1 = new TCanvas("c1","10-100 Ratios",2000,800);                                                                                                                                                                                                        
    c1->Divide(5,2);                                                                                                                                                                                                                                                 
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
    }                                                                                                                                                                                                                                                                
    */

  TCanvas *c11 = new TCanvas("c11","f2 ratios",2000,800);
    c11->Divide(5,2);

    for(int i = 1;i < nbins_q2-2; i++){                                                                                                                                                                                                                     
        hF2_x_eAu_er[i]->SetFillColor(1);
        hF2_x_eAu_er[i]->SetFillStyle(3004);
        hF2_x_eAu_er[i]->SetMarkerSize(0);
        hF2_x_eAu_ernew[i]->SetFillColor(38);
        hF2_x_eAu_ernew[i]->SetFillStyle(1000);
        hF2_x_eAu_ernew[i]->SetMarkerSize(0);
        hF2_x_eAu_pseudo[i]->GetYaxis()->SetRangeUser(0.5,2);
        hF2_x_eAu_pseudo[i]->GetXaxis()->SetRangeUser(-3.5,0);
        hF2_x_eAu_pseudo[i]->GetYaxis()->SetNdivisions(5);
        hF2_x_eAu_pseudo[i]->GetXaxis()->SetTitle("log(x_{B})");
        hF2_x_eAu_pseudo[i]->GetYaxis()->SetTitle("F_{2}^{c#bar{c}}(x,Q^{2}) R_{eA}");
        c11->cd(i);
        hF2_x_eAu_pseudo[i]->Draw("PE X0 same");
        hF2_x_eAu_er[i]->Draw("E2 same");
        hF2_x_eAu_ernew[i]->Draw("E2 same");
        hF2_x_eAu_pseudo[i]->Draw("PE X0 same");
        line->Draw("same L");
        lat.DrawLatex(-3.2,1.7,Form("Q^{2} = %1.f GeV^{2}",qq[i]));
    }
    
    TCanvas *c111 = new TCanvas("c111","gluon ratios",2000,800);
    c111->Divide(5,2);
    if(save)TFile *fout = new TFile("EPPS_Gluon_NewWeights_1k.root","RECREATE");




    for(int i = 1;i < 2; i++){
//	for(int i = 1;i < nbins_q2-2; i++){
	c111->cd(i);
	TH1F* True = (TH1F*)N_GLUON[0]->ProjectionX("True",i+1,i+1);
	TH1F* Mod = (TH1F*)N_GLUON[0]->ProjectionX("Mod",i+1,i+1);
	TH1F* Ini = (TH1F*)N_GLUON[0]->ProjectionX("Ini",i+1,i+1);
	
	//True->SetLineStyle(7);
	True->GetYaxis()->SetRangeUser(0,2);
	True->GetYaxis()->SetTitle("Nuclear Gluon Ratio");
	True->GetXaxis()->SetTitle("log(x)");
	Ini->SetFillColor(1);
	Ini->SetFillStyle(3004);
	Mod->SetFillColor(38);
	Mod->SetFillStyle(1000);
	Ini->SetMarkerSize(0);
	Mod->SetMarkerSize(0);
	Mod->SetLineColor(38);
	cout <<"test " << i << endl;
	for(int ii = 1; ii< Ini->GetNbinsX()+1; ii++){
	    double sum1 = True->GetBinContent(ii);
	    double sum2 = True->GetBinContent(ii);
	    double sum21 = 0;
	    double sum22 = 0;

	    for(int j = 0; j< 40/2; j++){
                TH1F *_temp2 = (TH1F*)N_GLUON_Weight[j]->ProjectionX("_temp2",i+1,i+1);
                TH1F *_temp1 = (TH1F*)N_GLUON_Comb[j]->ProjectionX("_temp1",i+1,i+1);
		sum21+=(sum1 - _temp1->GetBinContent(ii))*(sum1 - _temp1->GetBinContent(ii));
                delete _temp1; delete _temp2;
            }

	    for(int j = 0; j< nvary; j++){
		double _val = N_GLUON_Replica[j]->GetBinContent(ii,i+1);
		if(_val == 0) cout << i << " " << j << " " << endl;
		if(_val<0) _val = 0 ;
		sum22+=REP_WEIGHTS->GetBinContent(j+1)*(sum2 - _val)*(sum2 - _val);
		//sum22+=(sum2 - _val)*(sum2 - _val);
		if(i==1 && ii==11) {
		    cout <<" >PDF DEBUG " << _val << " " << sum22 <<  " " << REP_WEIGHTS->GetBinContent(j+1) << endl;
		}
	    }
	    Mod->SetBinContent(ii,sum2);
	    Ini->SetBinContent(ii,sum1);
	    Mod->SetBinError(ii,sqrt(sum22/nvary));
	    Ini->SetBinError(ii,sqrt(sum21));
	    if(sum2>0){
		ERRS_W[i]->SetBinContent(ii,sqrt(sum22/nvary)/sum2);
                ERRS[i]->SetBinContent(ii,sqrt(sum21)/sum1);

            }else{
                ERRS_W[i]->SetBinContent(ii,0);
                ERRS[i]->SetBinContent(ii,0);
            }
	}

	True->DrawClone("hist same");
	Ini->DrawClone("E2 same");
	Mod->DrawClone("E2 same");
	True->DrawClone("hist same");
	line->Draw("same L");
	lat.DrawLatex(-3.7,1.7,Form("Q^{2} = %1.f GeV^{2}",qq[i]));

	    if(save){
	Ini->Write(Form("Base_EPPS_%i",i));
        Mod->Write(Form("New_EPPS_%i",i));
        ERRS_W[i]->Write(Form("ERRS_New_EPPS_%i",i));
        ERRS[i]->Write(Form("ERRS_EPPS_%i",i));
	}
    }
/*
    TCanvas *c1111 = new TCanvas("c1111","gluon ratios errors",2000,800);
    c1111->Divide(5,2);
    for(int i = 1;i < nbins_q2-2; i++){
        c1111->cd(i);
        ERRS[i]->GetYaxis()->SetTitle("Rel. Gluon nPDF Uncertainty");
	ERRS[i]->GetYaxis()->SetRangeUser(0,1.5);
        ERRS[i]->GetXaxis()->SetTitle("log(x)");
        ERRS_W[i]->SetLineColor(38);
        ERRS[i]->SetLineStyle(7);
        ERRS[i]->DrawClone("hist");
        ERRS_W[i]->DrawClone("hist same");
	lat.DrawLatex(-3.7,1.7,Form("Q^{2} = %1.f GeV^{2}",qq[i]));
    }
*/

}
