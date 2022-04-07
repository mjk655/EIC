void PlotLeptonCorrelation(int save=1){
    gROOT->ProcessLine(".x ~/myStyle.C");
    gStyle->SetPalette(51);
    TFile *f = new TFile("./out_10_100_IC_Q210.root");
    TH3F* D0_Lepton_Corr_PGF = (TH3F*)f->Get("D0L_nRedCS_T_Acc");
    TH3F* D0_Lepton_Corr_L0 = (TH3F*)f->Get("D0L_nRedCS_Topo_T_L0");
    //TH3F* D0_Lepton_Corr_PGF = (TH3F*)f->Get("D0_Lepton_Corr_PGF");
    //TH3F* D0_Lepton_Corr_L0 = (TH3F*)f->Get("D0_Lepton_Corr_L0");
    TH1F* PGF = (TH1F*)D0_Lepton_Corr_PGF->Project3D("z");
    TH1F* L0 = (TH1F*)D0_Lepton_Corr_L0->Project3D("z");

    PGF->SetLineColor(kBlue);
    L0->SetLineColor(kRed);
    L0->SetLineStyle(7);
    TLegend *leg = new TLegend(0.2,0.6,0.5,0.9);
    leg->SetHeader("PYTHIA e+p 10+100 GeV; D^{0} |#eta|<3");
    leg->AddEntry(PGF,"Photon-Gluon Fusion","l");
    leg->AddEntry(L0,"Extrinsic #gamma*+#it{c}#rightarrow#it{c}","l");
    leg->SetTextSize(0.05);
    TCanvas *c1 = new TCanvas("c1","c1");
    PGF->GetYaxis()->SetTitle("Arb. Units");
    PGF->GetYaxis()->SetRangeUser(0,L0->GetMaximum()*1.3);
    PGF->GetXaxis()->SetTitle("|#phi(D^{0}) - #phi(e')|/#pi");
    norm(PGF);
    norm(L0);
    PGF->GetYaxis()->SetRangeUser(0,L0->GetMaximum()*1.3);
    PGF->Draw("hist");
    L0->Draw("hist same");
    leg->Draw("same");


}
void norm(TH1F *h){
    double norm = h->Integral();
        
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
