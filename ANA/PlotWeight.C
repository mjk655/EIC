void PlotWeight(){


    gROOT->ProcessLine(".x ~/myStyle.C");
    double csum1 = 0;
    double csum2 = 0;

    TFile *fin2 = new TFile("results2_NewWeights_Sys.root");//results11_NewWeights_Sys_highstat.root");//results11_NewWeights.root","READ");                                                                                                                                                            
                                                                                                                                                                                                                                                                
    TH1F *REP_WEIGHTS2 = (TH1F*)fin2->Get("REP_WEIGHTS");
    REP_WEIGHTS2->SetName("REP_WEIGHTS1");
    TH1F *chi22 = (TH1F*)fin2->Get("chi2");
    chi22->SetName("chi22");


    TFile *fin1 = new TFile("results2_NewWeights_Sys.root");//results11_NewWeights_Sys_Chi2.root");//results11_NewWeights.root","READ");                                                                                                                                                                  
    TH1F *REP_WEIGHTS1 = (TH1F*)fin1->Get("REP_WEIGHTS");
    REP_WEIGHTS1->SetName("REP_WEIGHTS1");
    TH1F *chi21 = (TH1F*)fin1->Get("chi2");
    chi21->SetName("chi21");


    TFile *fin = new TFile("results2_NewWeights_Sys.root");//results11_NewWeights.root","READ");                                                                                                                                                                   
    TH1F *REP_WEIGHTS = (TH1F*)fin->Get("REP_WEIGHTS");
    TH1F *chi2 = (TH1F*)fin->Get("chi2");

    getEfc(REP_WEIGHTS,"Sys");
    getEfc(REP_WEIGHTS1,"noT");
    getEfc(REP_WEIGHTS2,"chi2");
    for(int i = 0;i<100;i++){
	double th = 0 + i*0.1;
	getW(chi21,th);
    }
    double ave = 0;
    
    TGraph *g2 = new TGraph(REP_WEIGHTS2->GetNbinsX());
    TGraph *g1 = new TGraph(REP_WEIGHTS1->GetNbinsX());
    TGraph *g = new TGraph(REP_WEIGHTS->GetNbinsX());

    TH1F *hh = new TH1F("hh","hh",300,-20,10);
    TH1F *hh1 = new TH1F("hh1","hh1",300,-20,10);
    TH1F *hh2 = new TH1F("hh2","hh2",300,-20,10);
    TH1F *hp = new TH1F("hp","hp",400,0,20);
    TH1F *hp1 = new TH1F("hp1","hp1",400,0,20);
    TH1F *hp2 = new TH1F("hp2","hp2",400,0,20);
    TH1F *hc = new TH1F("hc","hc",400,0,200);
    TH1F *hcw = new TH1F("hcw","hcw",400,0,200);
    TH1F *hcw1 = new TH1F("hcw1","hcw1",400,0,200);
    TH1F *hcw2 = new TH1F("hcw2","hcw2",400,0,200);
    TH1F *hc1 = new TH1F("hc1","hc1",400,0,200);
    TH1F *hc2 = new TH1F("hc2","hc2",400,0,200);
    for(int i = 1;i<REP_WEIGHTS->GetNbinsX()+1;i++){
	hc->Fill(chi2->GetBinContent(i)/43.);
	hcw->Fill(chi2->GetBinContent(i)/43.,REP_WEIGHTS->GetBinContent(i));
	//hp->Fill(1./REP_WEIGHTS->GetNbinsX() * REP_WEIGHTS->GetBinContent(i));
	ave +=chi2->GetBinContent(i)/43./REP_WEIGHTS->GetNbinsX(); 
	if(REP_WEIGHTS->GetBinContent(i)>0){
	    hh->Fill(log10(REP_WEIGHTS->GetBinContent(i)));
	    g->SetPoint(i,chi2->GetBinContent(i)/43.,log10(REP_WEIGHTS->GetBinContent(i)));
	    if(chi2->GetBinContent(i)/43.<1.5)csum1+=1;//REP_WEIGHTS->GetBinContent(i);
	    else csum2+=1;//REP_WEIGHTS->GetBinContent(i);
	}
    }
    for(int i = 1;i<REP_WEIGHTS1->GetNbinsX()+1;i++){
	hc1->Fill(chi21->GetBinContent(i)/43.);
	hcw1->Fill(chi21->GetBinContent(i)/43.,REP_WEIGHTS1->GetBinContent(i));
	if(REP_WEIGHTS1->GetBinContent(i)>0){
	    hh1->Fill(log10(REP_WEIGHTS1->GetBinContent(i)));
	    g1->SetPoint(i,chi21->GetBinContent(i)/43.,log10(REP_WEIGHTS1->GetBinContent(i)));
	    //if(chi21->GetBinContent(i)/43.<1.5)csum1+=REP_WEIGHTS1->GetBinContent(i);
            //else csum2+=REP_WEIGHTS1->GetBinContent(i);
	}
    }
    for(int i = 1;i<REP_WEIGHTS2->GetNbinsX()+1;i++){
	hc2->Fill(chi22->GetBinContent(i)/43.);
	hcw2->Fill(chi22->GetBinContent(i)/43.,REP_WEIGHTS2->GetBinContent(i));
	//hp2->Fill(1./REP_WEIGHTS2->GetNbinsX() * REP_WEIGHTS2->GetBinContent(i));
	if(REP_WEIGHTS2->GetBinContent(i)>0){
            hh2->Fill(log10(REP_WEIGHTS2->GetBinContent(i)));
	    g2->SetPoint(i,chi22->GetBinContent(i)/43.,log10(REP_WEIGHTS2->GetBinContent(i)));
	    //if(chi22->GetBinContent(i)/43.<=1.5)csum1+=1;//REP_WEIGHTS2->GetBinContent(i);
            //else csum2+=1;//REP_WEIGHTS2->GetBinContent(i);
	}
    }

    for(int i = 1;i<hp->GetNbinsX()+1;i++){
	double low = hp->GetBinLowEdge(i);
	double high= hp->GetBinLowEdge(i)+hp->GetBinWidth(i);
	double sum = 0;

	for(int i = 1;i<REP_WEIGHTS->GetNbinsX()+1;i++){
	    double val = REP_WEIGHTS->GetBinContent(i);
	    double cc2 = chi2->GetBinContent(i)/43.;
	
	    if(cc2>= low && cc2<high)sum+=val;
	}
	cout <<low << " " << high << " " << sum << endl;;
	hp->SetBinContent(i,sum);
    }

    TCanvas *c13 = new TCanvas("c13","c13");
    hp2->SetLineColor(kRed);
    hp->Draw("hist");
    hp2->Draw("hist same");



    TCanvas *c1 = new TCanvas("c1","c1");
    hh->GetXaxis()->CenterTitle(1);
    hh->GetYaxis()->CenterTitle(1);
    hh->SetLineColor(kBlue);
    hh2->SetLineColor(kRed);
    hh->GetXaxis()->SetTitle("log_{10}(w_{k})");
    hh->GetYaxis()->SetTitle("Normalized Counts");
    hh->DrawNormalized("hist");
    hh1->DrawNormalized("same hist");
    //hh->Draw("hist");
    //hh1->Draw("same hist");
    // hc2->SetLineColor(kRed);
    hh2->DrawNormalized("same hist");
    TCanvas *c2 = new TCanvas("c2","c2");

    hcw2->SetLineColor(kRed);
    hc->GetXaxis()->SetTitle("#chi^{2}/n.d.f.");
    hc->GetYaxis()->SetTitle("Arb. Units");
    hc->GetXaxis()->CenterTitle(1);
    hc->GetYaxis()->CenterTitle(1);
    hc2->Draw("hist");                                                                                                                                                                                                                                        
    hcw2->Draw("same hist");                                                                                                                                                                                                                                  
    //hc2->DrawNormalized("hist");
    //gPad->SetLogy();
    //gPad->SetLogx();
    //hc1->Draw("same hist");

    TCanvas *c22 = new TCanvas("c22","c22");
    g->GetXaxis()->SetTitle("#chi^{2}/n.d.f.");
    g->GetYaxis()->SetTitle("log_{10}(w_{k})");
    g->GetXaxis()->CenterTitle(1);
    g->GetYaxis()->CenterTitle(1);
    g->GetXaxis()->SetLimits(-100,1000);
    g->GetYaxis()->SetRangeUser(-3,3);
    g->GetXaxis()->SetRangeUser(0.5,11);
    g->SetMarkerColor(kBlue);
    g2->SetMarkerColor(kRed);
    g1->SetMarkerStyle(2);
    g2->SetMarkerStyle(30);
    g->SetMarkerStyle(24);
    TLegend *leg = new TLegend(0.5,0.65,0.8,0.8);
    leg->SetTextSize(0.045);
    leg->AddEntry(g,"Giele-Keller","p");
    leg->AddEntry(g2,"Scaled Giele-Keller","p");
    leg->AddEntry(g1,"#chi^{2}","p");
    //g->GetYaxis()->SetRangerUser(-3,2.5);
    //g->GetXaxis()->SetRangerUser(0.5,10);
    g->Draw("AP");
    g1->Draw("sameP");
    g2->Draw("sameP");
    leg->Draw("same");
    gPad->SetLogx();
    gPad->RedrawAxis();

    cout << "chi2 ave " << hc->GetMean() << " " << ave <<endl;
    cout << "weighter chi2 ave " << hcw->GetMean() << endl;
    
    cout << "chi2 rms " << hc->GetRMS() << endl;
    cout << "weighter chi2 ave " << hcw2->GetMean() << endl;

    cout << "chi2 ratio <1.5 " << csum1 << " outside " << csum2 << endl;;

}
void getEfc(TH1F *h, char mess[100]){
    double sum = 0;
    double n = h->GetNbinsX();
    for(int i = 1;i<h->GetNbinsX()+1;i++){
	double val  = h->GetBinContent(i);
	if(val>0)sum += 1/n * val * log(n/val);
    }

    cout << mess << " " << n << " " << exp(sum)  << endl;

}
void getW(TH1F *h,double t){

    double d = 0 ;
    
    for(int i = 1; i< h->GetNbinsX()+1;i++){
	double _chi2 = h->GetBinContent(i);
	if(_chi2>t*43.)continue;
	d += TMath::Exp(-_chi2/2.)/h->GetNbinsX();
    }
    double sum = 0;

    for(int i =1; i< h->GetNbinsX()+1;i++){
        double _chi2 = h->GetBinContent(i);
        if(_chi2>t*43.)continue;
	sum += TMath::Exp(-_chi2/2.)/d;
    }
    
    cout << "  " << t << " SUM " <<  d << endl;

} 
