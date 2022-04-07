#include "PlotErrors.h"

void PlotErrors(){
    gROOT->ProcessLine(".x ~/myStyle.C");
    gStyle->SetPalette(51);
    int const nvary = 20;
    int const nbins_q2 = 13;
    int const nbins_x  = 20;
    double binning_q2[nbins_q2+1];
    double binning_x[nbins_x+1];
    setXbins(nbins_x,binning_x);
    setQ2bins(nbins_q2,binning_q2);
    TFile *partonf = new TFile("PartonWeights.root","READ");
    TH2F* N_GLUON = (TH2F*)partonf->Get("N_GLUON");N_GLUON->SetDirectory(0);
    
    TFile *f1 = new TFile("./Weight_File.root","READ");
    TH2F *Weights_10100    = (TH2F*)f1->Get("Weights_10100");
    TH2F *Weights_541      = (TH2F*)f1->Get("Weights_541");
    TH2F *Errors_10100     = (TH2F*)f1->Get("Errors_10100");
    TH2F *Errors_541       = (TH2F*)f1->Get("Errors_541");
    TH2F *Weights_10100_eAu= (TH2F*)f1->Get("Weights_10100_eAu");
    TH2F *Weights_541_eAu  = (TH2F*)f1->Get("Weights_541_eAu");
    TH2F *Errors_10100_eAu = (TH2F*)f1->Get("Errors_10100_eAu");
    TH2F *Errors_541_eAu   = (TH2F*)f1->Get("Errors_541_eAu");

    TH1F *hq2bin = new TH1F("hq2bin","hq2bin",nbins_q2,binning_q2);

    for(int i = 0;i < nbins_q2; i++){
	hCS_x1[i] = new TH1F(Form("hCS_x1_%i",i),Form("hCS_x1_%i",i),nbins_x,binning_x);
        hCS_x2[i] = new TH1F(Form("hCS_x2_%i",i),Form("hCS_x2_%i",i),nbins_x,binning_x);
	hF2_x[i] = new TH1F(Form("hF2_x_%i",i),Form("hF2_x_%i",i),nbins_x,binning_x);
	hCS_x1_eAu[i] = new TH1F(Form("hCS_x1_%i_eAu",i),Form("hCS_x1_%i_eAu",i),nbins_x,binning_x);
	hCS_x2_eAu[i] = new TH1F(Form("hCS_x2_%i_eAu",i),Form("hCS_x2_%i_eAu",i),nbins_x,binning_x);
	hF2_x_eAu[i] = new TH1F(Form("hF2_x_%i_eAu",i),Form("hF2_x_%i_eAu",i),nbins_x,binning_x);
	hCS_x1_eAu_er[i] = new TH1F(Form("hCS_x1_%i_eAu_er",i),Form("hCS_x1_%i_eAu_er",i),nbins_x,binning_x);
	hCS_x2_eAu_er[i] = new TH1F(Form("hCS_x2_%i_eAu_er",i),Form("hCS_x2_%i_eAu_er",i),nbins_x,binning_x);
	hF2_x_eAu_er[i] = new TH1F(Form("hF2_x_%i_eAu_er",i),Form("hF2_x_%i_eAu_er",i),nbins_x,binning_x);
	hCS_x1_eAu_toy[i] = new TH2F(Form("hCS_x1_%i_eAu_toy",i),Form("hCS_x1_%i_eAu_toy",i),nbins_x,binning_x,200,0,2);
	hCS_x2_eAu_toy[i] = new TH2F(Form("hCS_x2_%i_eAu_toy",i),Form("hCS_x2_%i_eAu_toy",i),nbins_x,binning_x,200,0,2);
	hF2_x_eAu_toy[i] = new TH2F(Form("hF2_x_%i_eAu_toy",i),Form("hF2_x_%i_eAu_toy",i),nbins_x,binning_x,200,0,2);
    }
//===================================================== Central values =================================================
    cout << "Filling central values " << endl;
    cout << "Filling 10x100 histograms " << endl;
    fill("D0Tree_10100_full.root",hCS_x1,hCS_x1_eAu,hq2bin,Weights_10100,Weights_10100_eAu,N_GLUON,N_GLUON,0,0.75);
    fill("D0Tree_10100_q210_full.root",hCS_x1,hCS_x1_eAu,hq2bin,Weights_10100,Weights_10100_eAu,N_GLUON,N_GLUON,0,0.75);
    cout << "Filling 5x41 histograms " << endl;
    fill("D0Tree_541_full.root",hCS_x2,hCS_x2_eAu,hq2bin,Weights_541,Weights_541_eAu,N_GLUON,N_GLUON,0,1);
    fill("D0Tree_541_q210_full.root",hCS_x2,hCS_x2_eAu,hq2bin,Weights_541,Weights_541_eAu,N_GLUON,N_GLUON,0,1);
    cout << "Setting proper errors " << endl;
    for(int i =1; i < nbins_q2-2;i++){
	for(int j = 1; j< nbins_x+1;j++){
	    hCS_x1[i]->SetBinError(j,Errors_10100->GetBinContent(j,i));
	    hCS_x2[i]->SetBinError(j,Errors_541->GetBinContent(j,i));
	    hCS_x1_eAu[i]->SetBinError(j,Errors_10100_eAu->GetBinContent(j,i));
            hCS_x2_eAu[i]->SetBinError(j,Errors_541_eAu->GetBinContent(j,i));
	}
    }
    cout << "Fitting F2 values " << endl;
    fitF2(nbins_q2,nbins_x,binning_q2,binning_x,hF2_x,hCS_x1,hCS_x2);
    fitF2(nbins_q2,nbins_x,binning_q2,binning_x,hF2_x_eAu,hCS_x1_eAu,hCS_x2_eAu);    
//===================================================== nPDF Errors ====================================================
    cout << ">> Starting PDF variation " << endl;
    for(int i = 0; i < nvary ; i++){
	// Generate sigma variation histogram
	TH2F* sigma_hist = (TH2F*)N_GLUON->Clone();
	for(int j=1;j<sigma_hist->GetNbinsY()+1;j++){
	    double sig = gRandom->Gaus(0,1);
	    for(int k=1;k<sigma_hist->GetNbinsX()+1;k++)
		sigma_hist->SetBinContent(k,j,sig);
	}
	TH1F *hCS_x1_eAu_temp[nbins_q2];
	TH1F *hCS_x2_eAu_temp[nbins_q2];
	TH1F *hF2_x_eAu_temp[nbins_q2];
	
	for(int ii = 0;ii < nbins_q2; ii++){
	    hCS_x1_eAu_temp[ii] = new TH1F(Form("hCS_x1_%i_eAu_temp",ii),Form("hCS_x1_%i_eAu_temp",ii),nbins_x,binning_x);
	    hCS_x2_eAu_temp[ii] = new TH1F(Form("hCS_x2_%i_eAu_temp",ii),Form("hCS_x2_%i_eAu_temp",ii),nbins_x,binning_x);
	    hF2_x_eAu_temp[ii] = new TH1F(Form("hF2_x_%i_eAu_temp",ii),Form("hF2_x_%i_eAu_temp",ii),nbins_x,binning_x);
	}
	fill("D0Tree_10100_full.root",hCS_x1,hCS_x1_eAu_temp,hq2bin,Weights_10100,Weights_10100_eAu,N_GLUON,sigma_hist,1,0.75);
	fill("D0Tree_10100_q210_full.root",hCS_x1,hCS_x1_eAu_temp,hq2bin,Weights_10100,Weights_10100_eAu,N_GLUON,sigma_hist,1,0.75);
	fill("D0Tree_541_full.root",hCS_x2,hCS_x2_eAu_temp,hq2bin,Weights_541,Weights_541_eAu,N_GLUON,sigma_hist,1,1);
	fill("D0Tree_541_q210_full.root",hCS_x2,hCS_x2_eAu_temp,hq2bin,Weights_541,Weights_541_eAu,N_GLUON,sigma_hist,1,1);
	fitF2(nbins_q2,nbins_x,binning_q2,binning_x,hF2_x_eAu_temp,hCS_x1_eAu_temp,hCS_x2_eAu_temp);
	
	for(int ii = 0;ii < nbins_q2-2; ii++){
	    getRat(hF2_x_eAu_temp[ii],hF2_x[ii],1);
	    getRat(hCS_x1_eAu_temp[ii],hCS_x1[ii],1);
	    getRat(hCS_x2_eAu_temp[ii],hCS_x2[ii],1);
	    for(int jj = 1;jj < nbins_x+1; jj++){
		double val1 = hCS_x1_eAu_temp[ii]->GetBinContent(jj);
		double val2 = hCS_x2_eAu_temp[ii]->GetBinContent(jj);
		double val3 = hF2_x_eAu_temp[ii]->GetBinContent(jj);
		if(val1>0)hCS_x1_eAu_toy[ii]->Fill(binning_x[jj-1]+0.00000001,val1);
		if(val2>0)hCS_x2_eAu_toy[ii]->Fill(binning_x[jj-1]+0.00000001,val2);
		if(val3>0)hF2_x_eAu_toy[ii]->Fill(binning_x[jj-1]+0.00000001,val3);
	    }
	}
	delete sigma_hist;
	for(int ii = 0;ii < nbins_q2; ii++){
	    delete hCS_x1_eAu_temp[ii];
	    delete hCS_x2_eAu_temp[ii]; 
	    delete hF2_x_eAu_temp[ii];
	}
    }
    cout << " Done with PDF varying " << endl;
    for(int ii = 0;ii < nbins_q2-2; ii++){
	for(int jj = 1;jj < nbins_x+1; jj++){
	    TH1F* temp1 = (TH1F*)hCS_x1_eAu_toy[ii]->ProjectionY("temp1",jj,jj);
	    TH1F* temp2 = (TH1F*)hCS_x2_eAu_toy[ii]->ProjectionY("temp2",jj,jj);
	    TH1F* temp3 = (TH1F*)hF2_x_eAu_toy[ii]->ProjectionY("temp3",jj,jj);
	    if(temp1->Integral()>0)hCS_x1_eAu_er[ii]->SetBinError(jj,temp1->GetRMS());
	    if(temp2->Integral()>0)hCS_x2_eAu_er[ii]->SetBinError(jj,temp2->GetRMS());
	    if(temp3->Integral()>0)hF2_x_eAu_er[ii]->SetBinError(jj,temp3->GetRMS());
	    delete temp1; delete temp2; delete temp3;
	}
    }
//============================ This is the drawing part =======================================================
    double qq[20] = {1.3,2,3.1,5,8,12,20,31,52,76,120};
    TF1 *line = new TF1("line","1",-100,100);
    line->SetLineStyle(7);
    TLatex lat;
    TLegend *leg = new TLegend(0.17,0.62,0.4,0.87);
    
    TCanvas *c11 = new TCanvas("c11","f2 ratios",2000,800);
    c11->Divide(5,2);
    
    for(int i = 1;i < nbins_q2-2; i++){
	getRat(hF2_x_eAu[i],hF2_x[i],0);
	for(int j = 1; j< nbins_x+1;j++)hF2_x_eAu_er[i]->SetBinContent(j,hF2_x_eAu[i]->GetBinContent(j));
	    
      	hF2_x_eAu_er[i]->SetFillColor(38);
	hF2_x_eAu_er[i]->SetFillStyle(1000);
	hF2_x_eAu_er[i]->SetMarkerSize(0);
	hF2_x_eAu[i]->GetYaxis()->SetRangeUser(0.5,2);
	hF2_x_eAu[i]->GetXaxis()->SetRangeUser(-3.5,0);
	hF2_x_eAu[i]->GetYaxis()->SetNdivisions(5);
	hF2_x_eAu[i]->GetXaxis()->SetTitle("log(x_{B})");
	hF2_x_eAu[i]->GetYaxis()->SetTitle("F_{2}^{c#bar{c}}(x,Q^{2}) R_{eA}");
	c11->cd(i);
	hF2_x_eAu[i]->Draw("PE X0 same");
	hF2_x_eAu_er[i]->Draw("E2 same");
        hF2_x_eAu[i]->Draw("PE X0 same");
	line->Draw("same L");
	lat.DrawLatex(-3.2,1.5,Form("Q^{2} = %1.f GeV^{2}",qq[i]));
    }

    TCanvas *c1 = new TCanvas("c1","10-100 Ratios",2000,800);
    c1->Divide(5,2);
    for(int i = 1;i<nbins_q2-2;i++){
	c1->cd(i);
	getRat(hCS_x1_eAu[i],hCS_x1[i],0);
	for(int j = 1; j< nbins_x+1;j++)hCS_x1_eAu_er[i]->SetBinContent(j,hCS_x1_eAu[i]->GetBinContent(j));
	hCS_x1_eAu_er[i]->SetFillColor(38);
	hCS_x1_eAu_er[i]->SetFillStyle(1000);
        hCS_x1_eAu_er[i]->SetMarkerSize(0);
	hCS_x1_eAu[i]->GetXaxis()->SetRangeUser(-3.5,0);
	hCS_x1_eAu[i]->GetYaxis()->SetRangeUser(0.5,2);
	hCS_x1_eAu[i]->GetYaxis()->SetNdivisions(5);
	hCS_x1_eAu[i]->GetXaxis()->SetTitle("log(x_{B})");
	hCS_x1_eAu[i]->GetYaxis()->SetTitle("#sigma^{c#bar{c}}_{r}(e+Au)/#sigma^{c#bar{c}}_{r}(e+p)");
	hCS_x1_eAu[i]->Draw("PE X0 same");
	hCS_x1_eAu_er[i]->Draw("E2 same");
        hCS_x1_eAu[i]->Draw("PE X0 same");
	lat.DrawLatex(-3.2,1.7,"e+p/Au 10#times100 GeV");
	lat.DrawLatex(-3.2,1.5,Form("Q^{2} = %1.f GeV^{2}",qq[i]));
	line->Draw("same");
    }

    TCanvas *c2 = new TCanvas("c2","5-41 Ratios",2000,800);
    c2->Divide(5,2);
    for(int i = 1;i<nbins_q2-2;i++){
	c2->cd(i);
        getRat(hCS_x2_eAu[i],hCS_x2[i],0);
	for(int j = 1; j< nbins_x+1;j++)hCS_x2_eAu_er[i]->SetBinContent(j,hCS_x2_eAu[i]->GetBinContent(j));		
	hCS_x2_eAu_er[i]->SetFillColor(38);
        hCS_x2_eAu_er[i]->SetFillStyle(1000);
        hCS_x2_eAu_er[i]->SetMarkerSize(0);
	hCS_x2_eAu[i]->GetXaxis()->SetRangeUser(-3.5,0);
        hCS_x2_eAu[i]->GetYaxis()->SetRangeUser(0.5,2);
	hCS_x2_eAu[i]->GetYaxis()->SetNdivisions(5);
	hCS_x2_eAu[i]->GetXaxis()->SetTitle("log(x_{B})");
	hCS_x2_eAu[i]->GetYaxis()->SetTitle("#sigma^{c#bar{c}}_{r}(e+Au)/#sigma^{c#bar{c}}_{r}(e+p)");
        hCS_x2_eAu[i]->Draw("PE X0 same");
	hCS_x2_eAu_er[i]->Draw("E2 same");
	hCS_x2_eAu[i]->Draw("PE X0 same");
	lat.DrawLatex(-3.2,1.7,"e+p/Au 5#times41 GeV");
	lat.DrawLatex(-3.2,1.5,Form("Q^{2} = %1.f GeV^{2}",qq[i]));
	line->Draw("same");
    }

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

void getRat(TH1F* h, TH1F* h1,int er){
    for(int i = 1;i<h->GetNbinsX()+1;i++){
	double val1 = h->GetBinContent(i);
	double er1 = h->GetBinError(i);
	double val2 = h1->GetBinContent(i);
	double er2 = h1->GetBinError(i);
	if(val1>0 && val2>0){
	    if(er==0){
		h->SetBinContent(i,val1/val2);
		h->SetBinError(i,val1/val2*sqrt(er1*er1/val1/val1+er2*er2/val2/val2));
	    }
	    else{
		h->SetBinContent(i,val1/val2);
                h->SetBinError(i,val1/val2*sqrt(er1*er1/val1/val1));
	    }
	}
	else{
	    h->SetBinContent(i,0);
	    h->SetBinError(i,0);
	}
    }
}
void fitF2(int nbins_q2, int nbins_x, double binning_q2[], double binning_x[],TH1F** hF2_x,TH1F** hCS_x1,TH1F** hCS_x2){
    for(int i = 0; i < nbins_q2-2;i++){
        for(int j = 1; j< nbins_x+1;j++){
            double x = TMath::Power(10,(binning_q2[i+1]-binning_q2[i])/2);
            double q = TMath::Power(10,(binning_x[j]-binning_x[j-1])/2);
            if(hCS_x1[i]->GetBinContent(j)>0 && hCS_x2[i]->GetBinContent(j)>0){
                double f2 = getF2(hCS_x1[i]->GetBinContent(j),hCS_x2[i]->GetBinContent(j),hCS_x1[i]->GetBinError(j),hCS_x2[i]->GetBinError(j),x,q);
                double f2er = getF2er(hCS_x1[i]->GetBinContent(j),hCS_x2[i]->GetBinContent(j),hCS_x1[i]->GetBinError(j),hCS_x2[i]->GetBinError(j),x,q);
                if(f2 !=1){
                    hF2_x[i]->SetBinContent(j,f2);
                    hF2_x[i]->SetBinError(j,f2er);
                }
            }
	}
    }
}

double getF2(double v1, double v2, double e1, double e2, double qq, double xx){
    
    double s1 = 4.*10.*100.;
    double s2 = 4.*5.*41.;
    double mean_y1 = qq/s1/xx;
    double mean_y2 = qq/s2/xx;
    double yy1 = mean_y1*mean_y1 / (1+ (1-mean_y1)*(1-mean_y1));
    double yy2 = mean_y2*mean_y2 / (1+ (1-mean_y2)*(1-mean_y2));
    double axis_x[2] = {yy1,yy2};
    double axis_ex[2] = {0,0};
    double axis_y[2] = {v1,v2};
    double axis_ey[2] = {e1,e2};
    TGraphErrors gr(2,axis_x,axis_y,axis_ex,axis_ey);
    TF1 fit("fit","[0]+[1]*x",0,1);
    fit.SetParLimits(0,0,0.5);
    gr.Fit("fit","M EX0");
    double error = fit.GetParError(0);
    double val=  fit.GetParameter(0);
    return val;
}
double getF2er(double v1, double v2, double e1, double e2, double qq, double xx){
    double s1 = 4.*10.*100.;
    double s2 = 4.*5.*41.;
    double mean_y1 = qq/s1/xx;
    double mean_y2 = qq/s2/xx;
    double yy1 = mean_y1*mean_y1 / (1+ (1-mean_y1)*(1-mean_y1));
    double yy2 = mean_y2*mean_y2 / (1+ (1-mean_y2)*(1-mean_y2));
    double axis_x[2] = {yy1,yy2};
    double axis_ex[2] = {0,0};
    double axis_y[2] = {v1,v2};
    double axis_ey[2] = {e1,e2};
    TGraphErrors gr(2,axis_x,axis_y,axis_ex,axis_ey);
    TF1 fit("fit","[0]+[1]*x",0,1);
    fit.SetParLimits(0,0,0.5);
    gr.Fit("fit","M EX0");
    double error = fit.GetParError(0);
    double val=  fit.GetParameter(0);
    return error;
}
void fill(char file1[100],TH1F** hCS_x1,TH1F** hCS_x1_eAu,TH1F* hq2bin,TH2F* Weights,TH2F *Weights_eAu,TH2F* N_GLUON,TH2F* vars,int vary,double ww){
    
    float mx; float mq2; float mxtp;float mw;
    TChain *tree = new TChain("tree","tree");
    tree->AddFile(file1);
    tree->SetBranchAddress("mq2",&mq2);
    tree->SetBranchAddress("mx",&mx);
    tree->SetBranchAddress("mw",&mw);
    tree->SetBranchAddress("mxtp",&mxtp);
    
    Long64_t iloop = tree->GetEntries()*ww;
    
    for(Long64_t i =0  ;i < iloop; i++){
        if(i%1000000 == 0) cout << "> On " << i << " out of " << iloop << endl;
	tree->GetEntry(i);

        int bin = hq2bin->FindBin(log10(mq2));
	
        double ww1 = Weights->GetBinContent(Weights->FindBin(log10(mx),log10(mq2)));
	double ww2 = Weights_eAu->GetBinContent(Weights_eAu->FindBin(log10(mx),log10(mq2)));
	double wwn = 1;
	
	int bin_g = N_GLUON->FindBin(log10(mxtp),log10(mq2));
	if(vary==0)
	    wwn = N_GLUON->GetBinContent(bin_g);
	else{
	    wwn = N_GLUON->GetBinContent(bin_g) * (1 + N_GLUON->GetBinError(bin_g)*vars->GetBinContent(bin_g));
	}
	
	if(vary==0)hCS_x1[bin]->Fill(log10(mx),ww1/ww);
	hCS_x1_eAu[bin]->Fill(log10(mx),ww2*wwn/ww);
    }
    delete tree;
}

