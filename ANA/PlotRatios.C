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
void PlotRatios(int save=1){
    gROOT->ProcessLine(".x ~/myStyle.C");
    gStyle->SetPalette(51);
    double scale2[nbins_q2]={0.015,0.025,0.025,0.03,0.05,0.1,0.2,0.4,0.6,1.5,4,2,3};
    double scale[nbins_q2]={0.015,0.025,0.025,0.03,0.05,0.1,0.2,0.4,0.6,1.5,4,2,3};
    setXbins();
    setQ2bins();
    TH1F *hCS_x1[nbins_q2];
    TH1F *hCS_x2[nbins_q2];
    TH1F *hCS_x1eAu[nbins_q2];
    TH1F *hCS_x2eAu[nbins_q2];
    TH1F *hF2_x[nbins_q2];
    TH1F *hF2_xeAu[nbins_q2];
    
    char name[100];char nname[100];
    //============================     
    TFile *f1 = new TFile("./ep_output.root","READ");
    for(int i = 0;i < nbins_q2; i++){
	sprintf(nname,"ep_10100_Q2Bin_%i",i);
	sprintf(name,"hCS_x1_%i",i);
	hCS_x1[i] = (TH1F*)f1->Get(nname);hCS_x1[i]->SetName(name);
	sprintf(nname,"ep_541_Q2Bin_%i",i);
	sprintf(name,"hCS_x2_%i",i);
        hCS_x2[i] = (TH1F*)f1->Get(nname);hCS_x2[i]->SetName(name);
	sprintf(nname,"ep_CharmF2_Q2Bin_%i",i);
	sprintf(name,"hF2_x_%i",i);
	hF2_x[i] = (TH1F*) f1->Get(nname);hF2_x[i]->SetName(name);
    }
//============================  
    TFile *f2 = new TFile("./eAu_output.root","READ");
    for(int i = 0;i < nbins_q2; i++){
        sprintf(nname,"eAu_10100_Q2Bin_%i",i);
        sprintf(name,"hCS_xeAu_%i",i);
        hCS_x1eAu[i] = (TH1F*)f2->Get(nname);hCS_x1eAu[i]->SetName(name);
        sprintf(nname,"eAu_541_Q2Bin_%i",i);
        sprintf(name,"hCS_x2eAu_%i",i);
        hCS_x2eAu[i] = (TH1F*)f2->Get(nname);hCS_x2eAu[i]->SetName(name);
        sprintf(nname,"eAu_CharmF2_Q2Bin_%i",i);
        sprintf(name,"hF2_xeAu_%i",i);
        hF2_xeAu[i] = (TH1F*) f2->Get(nname);hF2_xeAu[i]->SetName(name);
    }
//============================ 
    double qq[20] = {1.3,2,3.1,5,8,12,20,31,52,76,120};
    TF1 *line = new TF1("line","1",-100,100);
    line->SetLineStyle(7);
    TLatex lat;
    TLegend *leg = new TLegend(0.17,0.62,0.4,0.87);

    TCanvas *c11 = new TCanvas("c11","f2 ratios",2000,800);
    c11->Divide(5,2);
    TF1 *line = new TF1("line","1",-100,100);
    for(int i = 1;i < nbins_q2-2; i++){
	getRat(hF2_xeAu[i],hF2_x[i]);
	hF2_xeAu[i]->GetYaxis()->SetRangeUser(0.5,2);
	hF2_xeAu[i]->GetXaxis()->SetRangeUser(-3.5,0);
	hF2_xeAu[i]->GetYaxis()->SetNdivisions(5);
	hF2_xeAu[i]->GetXaxis()->SetTitle("log(x_{B})");
	hF2_xeAu[i]->GetYaxis()->SetTitle("F_{2}^{c#bar{c}} R_{eA}(x,Q^{2})");
	c11->cd(i);
	hF2_xeAu[i]->Draw("PE X0 same");
	line->Draw("same L");
	lat.DrawLatex(-3.2,1.5,Form("Q^{2} = %1.f GeV^{2}",qq[i]));
    }

    TCanvas *c1 = new TCanvas("c1","10-100 Ratios",2000,800);
    c1->Divide(5,2);
    for(int i = 1;i<nbins_q2-2;i++){
	c1->cd(i);
	getRat(hCS_x1eAu[i],hCS_x1[i]);
	hCS_x1eAu[i]->GetXaxis()->SetRangeUser(-3.5,0);
	hCS_x1eAu[i]->GetYaxis()->SetRangeUser(0.5,2);
	hCS_x1eAu[i]->GetYaxis()->SetNdivisions(5);
	hCS_x1eAu[i]->GetXaxis()->SetTitle("log(x_{B})");
	hCS_x1eAu[i]->GetYaxis()->SetTitle("#sigma^{c#bar{c}}_{r}(e+Au)/#sigma^{c#bar{c}}_{r}(e+p)");
	hCS_x1eAu[i]->Draw("PE X0 same");
	lat.DrawLatex(-3.2,1.7,"e+p/Au 10+100 GeV");
	lat.DrawLatex(-3.2,1.5,Form("Q^{2} = %1.f GeV^{2}",qq[i]));
	line->Draw("same");
    }

    TCanvas *c2 = new TCanvas("c2","5-41 Ratios",2000,800);
    c2->Divide(5,2);
    for(int i = 1;i<nbins_q2-2;i++){
	c2->cd(i);
        getRat(hCS_x2eAu[i],hCS_x2[i]);
        hCS_x2eAu[i]->GetXaxis()->SetRangeUser(-3.5,0);
        hCS_x2eAu[i]->GetYaxis()->SetRangeUser(0.5,2);
	hCS_x2eAu[i]->GetYaxis()->SetNdivisions(5);
	hCS_x2eAu[i]->GetXaxis()->SetTitle("log(x_{B})");
	hCS_x2eAu[i]->GetYaxis()->SetTitle("#sigma^{c#bar{c}}_{r}(e+Au)/#sigma^{c#bar{c}}_{r}(e+p)");
        hCS_x2eAu[i]->Draw("PE X0 same");
	lat.DrawLatex(-3.2,1.7,"e+p/Au 5+41 GeV");
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

void getRat(TH1F* h, TH1F* h1){
    for(int i = 1;i<h->GetNbinsX()+1;i++){
	double val1 = h->GetBinContent(i);
	double er1 = h->GetBinError(i);
	double val2 = h1->GetBinContent(i);
	double er2 = h1->GetBinError(i);
	if(val1>0 && val2>0){
	    if(1){
		h->SetBinContent(i,val1/val2);
		h->SetBinError(i,val1/val2*sqrt(er1*er1/val1/val1+er2*er2/val2/val2));
	    }
	}
    }
}
