#include "PlotErrors.h"

void DoAna_NNPDF(){
    bool save = true;
    gROOT->ProcessLine(".x ~/myStyle.C");
    gStyle->SetPalette(51);

    double sys = gRandom->Gaus(0,0.04);    
    int const nvary =900;//
    int const nbins_q2 = 13;
    int const nbins_x  = 20;
    TH1F* chi2 = new TH1F("chi2","chi2",nvary,0,nvary);//chi2 values for re-weighting
    TH1F *REP_WEIGHTS = new TH1F("REP_WEIGHTS","REP_WEIGHTS",nvary,0,nvary);
    double wk[nvary];
    double rik[nvary][nvary];//
    double binning_q2[nbins_q2+1];
    double binning_x[nbins_x+1];
    setXbins(nbins_x,binning_x);
    setQ2bins(nbins_q2,binning_q2);
    TFile *partonf = new TFile("GluonWeights_NNPDF.root","READ");
    
    TH2F* N_GLUON[nvary+1];
    TH2F* N_GLUON_Replica[nvary];
    TH2F* N_GLUON_Weight[nvary+1];
    TH2F* N_GLUON_Comb[nvary+1];
    

    for(int i = 0;i <nvary+1;i++){// only need 0 - 97
	N_GLUON[i] = (TH2F*)partonf->Get(Form("N_GLUON_%i",i));
	N_GLUON[i]->SetDirectory(0);
    }
    partonf->Close();
    TFile *f1 = new TFile("./Weight_File.root","READ");
    TH2F *Weights_10100    = (TH2F*)f1->Get("Weights_10100");Weights_10100->SetDirectory(0);
    TH2F *Weights_541      = (TH2F*)f1->Get("Weights_541");Weights_541->SetDirectory(0);
    TH2F *Errors_10100     = (TH2F*)f1->Get("Errors_10100");Errors_10100->SetDirectory(0);
    TH2F *Errors_541       = (TH2F*)f1->Get("Errors_541");Errors_541->SetDirectory(0);
    TH2F *Weights_10100_eAu= (TH2F*)f1->Get("Weights_10100_eAu");Weights_10100_eAu->SetDirectory(0);
    TH2F *Weights_541_eAu  = (TH2F*)f1->Get("Weights_541_eAu");Weights_541_eAu->SetDirectory(0);
    TH2F *Errors_10100_eAu = (TH2F*)f1->Get("Errors_10100_eAu");Errors_10100_eAu->SetDirectory(0);
    TH2F *Errors_541_eAu   = (TH2F*)f1->Get("Errors_541_eAu");Errors_541_eAu->SetDirectory(0);
    f1->Close();
    TH1F *hq2bin = new TH1F("hq2bin","hq2bin",nbins_q2,binning_q2);
    for(int i = 0;i < nbins_q2+5; i++){
	hCS_x1[i] = new TH1F(Form("hCS_x1_%i",i),Form("hCS_x1_%i",i),nbins_x,binning_x);
        hCS_x2[i] = new TH1F(Form("hCS_x2_%i",i),Form("hCS_x2_%i",i),nbins_x,binning_x);
	hF2_x[i] = new TH1F(Form("hF2_x_%i",i),Form("hF2_x_%i",i),nbins_x,binning_x);
	hCS_x1_eAu[i] = new TH1F(Form("hCS_x1_%i_eAu",i),Form("hCS_x1_%i_eAu",i),nbins_x,binning_x);
	hCS_x2_eAu[i] = new TH1F(Form("hCS_x2_%i_eAu",i),Form("hCS_x2_%i_eAu",i),nbins_x,binning_x);
	hF2_x_eAu[i] = new TH1F(Form("hF2_x_%i_eAu",i),Form("hF2_x_%i_eAu",i),nbins_x,binning_x);
	hCS_x1_eAu_pseudo[i] = new TH1F(Form("hCS_x1_%i_eAu_pseudo",i),Form("hCS_x1_%i_eAu_pseudo",i),nbins_x,binning_x);
        hCS_x2_eAu_pseudo[i] = new TH1F(Form("hCS_x2_%i_eAu_pseudo",i),Form("hCS_x2_%i_eAu_pseudo",i),nbins_x,binning_x);
        hF2_x_eAu_pseudo[i] = new TH1F(Form("hF2_x_%i_eAu_pseudo",i),Form("hF2_x_%i_eAu_pseudo",i),nbins_x,binning_x);

	hCS_x1_eAu_er[i] = new TH1F(Form("hCS_x1_%i_eAu_er",i),Form("hCS_x1_%i_eAu_er",i),nbins_x,binning_x);
	hCS_x2_eAu_er[i] = new TH1F(Form("hCS_x2_%i_eAu_er",i),Form("hCS_x2_%i_eAu_er",i),nbins_x,binning_x);
	hF2_x_eAu_er[i] = new TH1F(Form("hF2_x_%i_eAu_er",i),Form("hF2_x_%i_eAu_er",i),nbins_x,binning_x);
	hCS_x1_eAu_ernew[i] = new TH1F(Form("hCS_x1_%i_eAu_ernew",i),Form("hCS_x1_%i_eAu_ernew",i),nbins_x,binning_x);
        hCS_x2_eAu_ernew[i] = new TH1F(Form("hCS_x2_%i_eAu_ernew",i),Form("hCS_x2_%i_eAu_ernew",i),nbins_x,binning_x);
        hF2_x_eAu_ernew[i] = new TH1F(Form("hF2_x_%i_eAu_ernew",i),Form("hF2_x_%i_eAu_ernew",i),nbins_x,binning_x);
	hCS_x1_eAu_toy[i] = new TH2F(Form("hCS_x1_%i_eAu_toy",i),Form("hCS_x1_%i_eAu_toy",i),nbins_x,binning_x,4000,0,4);
	hCS_x2_eAu_toy[i] = new TH2F(Form("hCS_x2_%i_eAu_toy",i),Form("hCS_x2_%i_eAu_toy",i),nbins_x,binning_x,4000,0,4);
	hF2_x_eAu_toy[i] = new TH2F(Form("hF2_x_%i_eAu_toy",i),Form("hF2_x_%i_eAu_toy",i),nbins_x,binning_x,4000,0,4);
	hCS_x1_eAu_toynew[i] = new TH2F(Form("hCS_x1_%i_eAu_toynew",i),Form("hCS_x1_%i_eAu_toynew",i),nbins_x,binning_x,40000,-.01,.01);
        hCS_x2_eAu_toynew[i] = new TH2F(Form("hCS_x2_%i_eAu_toynew",i),Form("hCS_x2_%i_eAu_toynew",i),nbins_x,binning_x,40000,-.01,.01);
        hF2_x_eAu_toynew[i] = new TH2F(Form("hF2_x_%i_eAu_toynew",i),Form("hF2_x_%i_eAu_toynew",i),nbins_x,binning_x,40000,-.01,.01);
    }
    cout << "===================================================================" << endl;
//===================================================== Central values =================================================
    cout << "===================================================================" << endl;
    cout << "Filling central values " << endl;
    cout << "Filling 10x100 histograms " << endl;
    fill("D0Tree_10100_new.root",hCS_x1,hCS_x1_eAu,hq2bin,Weights_10100,Weights_10100_eAu,N_GLUON[0],0,1);
    cout << "Filling 5x41 histograms " << endl;
    fill("D0Tree_541_new.root",hCS_x2,hCS_x2_eAu,hq2bin,Weights_541,Weights_541_eAu,N_GLUON[0],0,1);
    cout << "Setting proper errors " << endl;
    for(int i =1; i < nbins_q2-2;i++){
	for(int j = 1; j< nbins_x+1;j++){
	    // Scaling here bc errors are too small right now
	    double scale = sqrt(10); // scaling to 5/1 fb-1 for e+p and e+Au
	    hCS_x1[i-1]->SetBinError(j,scale*Errors_10100->GetBinContent(j,i));
	    hCS_x2[i-1]->SetBinError(j,scale*Errors_541->GetBinContent(j,i));
	    hCS_x1_eAu[i-1]->SetBinError(j,scale*Errors_10100_eAu->GetBinContent(j,i));
            hCS_x2_eAu[i-1]->SetBinError(j,scale*Errors_541_eAu->GetBinContent(j,i));
	}
    }
    cout << "Fitting F2 values " << endl;
    fitF2(nbins_q2,nbins_x,binning_q2,binning_x,hF2_x,hCS_x1,hCS_x2);
    fitF2(nbins_q2,nbins_x,binning_q2,binning_x,hF2_x_eAu,hCS_x1_eAu,hCS_x2_eAu);  
    cout << "calculating N points to re-weight with" << endl;
    int Npts = 0; 
    for(int ii = 1;ii < nbins_q2-2; ii++){
	getRat(hF2_x_eAu[ii],hF2_x[ii],0);
        getRat(hCS_x1_eAu[ii],hCS_x1[ii],0);
        getRat(hCS_x2_eAu[ii],hCS_x2[ii],0);
	for(int jj = 1;jj < nbins_x+1; jj++){
	    if(hF2_x[ii]->GetBinContent(jj)>0 && hF2_x_eAu[ii]->GetBinContent(jj)>0){
		if(hF2_x_eAu[ii]->GetBinError(jj)/hF2_x_eAu[ii]->GetBinContent(jj)<0.4 && hF2_x_eAu[ii]->GetBinContent(jj)<2){
		    Npts++;
		}
	    }
	}
	
	for(int jj = 1;jj < nbins_x+1; jj++){
	    double val = hF2_x_eAu[ii]->GetBinContent(jj);
	    double er = hF2_x_eAu[ii]->GetBinError(jj);
	    double newval = gRandom->Gaus(val,er);
	    newval += sys*val;
	    hF2_x_eAu_pseudo[ii]->SetBinContent(jj,newval);
	    hF2_x_eAu_pseudo[ii]->SetBinError(jj,hF2_x_eAu[ii]->GetBinError(jj));
	    val = hCS_x1_eAu[ii]->GetBinContent(jj);
            er = hCS_x1_eAu[ii]->GetBinError(jj);
            newval = gRandom->Gaus(val,er);
	    newval += sys*val;
	    hCS_x1_eAu_pseudo[ii]->SetBinContent(jj,newval);
	    hCS_x1_eAu_pseudo[ii]->SetBinError(jj,hCS_x1_eAu[ii]->GetBinError(jj));
	    val = hCS_x2_eAu[ii]->GetBinContent(jj);
            er = hCS_x2_eAu[ii]->GetBinError(jj);
            newval = gRandom->Gaus(val,er);
	    newval += sys*val;
	    hCS_x2_eAu_pseudo[ii]->SetBinContent(jj,newval);
	    hCS_x2_eAu_pseudo[ii]->SetBinError(jj,hCS_x2_eAu[ii]->GetBinError(jj));
	    hF2_x_eAu_er[ii]->SetBinContent(jj,hF2_x_eAu[ii]->GetBinContent(jj));
	    hCS_x1_eAu_er[ii]->SetBinContent(jj,hCS_x1_eAu[ii]->GetBinContent(jj));
	    hCS_x2_eAu_er[ii]->SetBinContent(jj,hCS_x2_eAu[ii]->GetBinContent(jj));
	    hF2_x_eAu_ernew[ii]->SetBinContent(jj,hF2_x_eAu[ii]->GetBinContent(jj));
            hCS_x1_eAu_ernew[ii]->SetBinContent(jj,hCS_x1_eAu[ii]->GetBinContent(jj));
            hCS_x2_eAu_ernew[ii]->SetBinContent(jj,hCS_x2_eAu[ii]->GetBinContent(jj));
	}
    }
    cout << "> N points for re-weight " << Npts << endl;
    cout << "===================================================================" << endl;
//===================================================== nPDF error on x-sec ====================================================
    cout << "===================================================================" << endl;
    cout << ">> Starting Eigen PDF variation " << endl;
    int iter = 0;
    for(int i = 1; i <nvary+1 ; i++){
	cout <<">>> On " << i << " out of " << nvary << " variations "  << endl;

	TH1F *hCS_x1_eAu_temp[nbins_q2];
        TH1F *hCS_x2_eAu_temp[nbins_q2];
        TH1F *hF2_x_eAu_temp[nbins_q2];
	
        for(int ii = 0;ii < nbins_q2; ii++){
            hCS_x1_eAu_temp[ii] = new TH1F(Form("hCS_x1_%i_eAu_temp",ii),Form("hCS_x1_%i_eAu_temp",ii),nbins_x,binning_x);
            hCS_x2_eAu_temp[ii] = new TH1F(Form("hCS_x2_%i_eAu_temp",ii),Form("hCS_x2_%i_eAu_temp",ii),nbins_x,binning_x);
            hF2_x_eAu_temp[ii] = new TH1F(Form("hF2_x_%i_eAu_temp",ii),Form("hF2_x_%i_eAu_temp",ii),nbins_x,binning_x);
        }
	TH2F* New_N_GLUON = (TH2F*) N_GLUON[i]->Clone("New_N_GLUON");
	//getPDFError(New_N_GLUON,N_GLUON,i);
	N_GLUON_Comb[iter] = (TH2F*) New_N_GLUON->Clone(Form("N_GLUON_Comb_%i",iter));
    
        fill("D0Tree_10100_new.root",hCS_x1,hCS_x1_eAu_temp,hq2bin,Weights_10100,Weights_10100_eAu,New_N_GLUON,1,1);
        fill("D0Tree_541_new.root",hCS_x2,hCS_x2_eAu_temp,hq2bin,Weights_541,Weights_541_eAu,New_N_GLUON,1,1);
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
        for(int ii = 0;ii < nbins_q2; ii++){
            delete hCS_x1_eAu_temp[ii];
            delete hCS_x2_eAu_temp[ii];
            delete hF2_x_eAu_temp[ii];
        }
	iter++;
    }
    cout << " Done with PDF varying \n" << endl;

    for(int ii = 0;ii < nbins_q2-2; ii++){
        for(int jj = 1;jj < nbins_x+1; jj++){
            TH1F* temp1 = (TH1F*)hCS_x1_eAu_toy[ii]->ProjectionY("temp1",jj,jj);
            TH1F* temp2 = (TH1F*)hCS_x2_eAu_toy[ii]->ProjectionY("temp2",jj,jj);
            TH1F* temp3 = (TH1F*)hF2_x_eAu_toy[ii]->ProjectionY("temp3",jj,jj);
            double c1 = hCS_x1_eAu_er[ii]->GetBinContent(jj);
            double c2 = hCS_x2_eAu_er[ii]->GetBinContent(jj);
            double c3 = hF2_x_eAu_er[ii]->GetBinContent(jj);
            double er1 = getError(temp1,c1,0);
            double er2 = getError(temp2,c2,0);
            double er3 = getError(temp3,c3,0);
            if(temp1->Integral()>0)hCS_x1_eAu_er[ii]->SetBinError(jj,er1);
            if(temp2->Integral()>0)hCS_x2_eAu_er[ii]->SetBinError(jj,er2);
            if(temp3->Integral()>0)hF2_x_eAu_er[ii]->SetBinError(jj,er3);
            delete temp1; delete temp2; delete temp3;
        }
    }
//===================================================== nPDF Replicas ====================================================
    cout << "===================================================================" << endl;
    cout << ">> Starting Replica PDF variation " << endl;
    for(int i = 0; i < nvary ; i++){
	cout <<">>> On " << i+1 << " out of " << nvary << " variations " << endl;
	double _chi2=0;
	TH1F *hCS_x1_eAu_temp[nbins_q2];
	TH1F *hCS_x2_eAu_temp[nbins_q2];
	TH1F *hF2_x_eAu_temp[nbins_q2];
	
	for(int ii = 0;ii < nbins_q2; ii++){
	    hCS_x1_eAu_temp[ii] = new TH1F(Form("hCS_x1_%i_eAu_temp",ii),Form("hCS_x1_%i_eAu_temp",ii),nbins_x,binning_x);
	    hCS_x2_eAu_temp[ii] = new TH1F(Form("hCS_x2_%i_eAu_temp",ii),Form("hCS_x2_%i_eAu_temp",ii),nbins_x,binning_x);
	    hF2_x_eAu_temp[ii] = new TH1F(Form("hF2_x_%i_eAu_temp",ii),Form("hF2_x_%i_eAu_temp",ii),nbins_x,binning_x);
	}
	
	TH2F* New_N_GLUON = (TH2F*) N_GLUON[i+1]->Clone("New_N_GLUON"); 
	//getReplica(New_N_GLUON,N_GLUON,40,rik,i);
	N_GLUON_Replica[i] = (TH2F*) New_N_GLUON->Clone(Form("N_GLUON_Replica_%i",i));
	
	fill("D0Tree_10100_new.root",hCS_x1,hCS_x1_eAu_temp,hq2bin,Weights_10100,Weights_10100_eAu,New_N_GLUON,1,1);
	fill("D0Tree_541_new.root",hCS_x2,hCS_x2_eAu_temp,hq2bin,Weights_541,Weights_541_eAu,New_N_GLUON,1,1);
	fitF2(nbins_q2,nbins_x,binning_q2,binning_x,hF2_x_eAu_temp,hCS_x1_eAu_temp,hCS_x2_eAu_temp);
	
	for(int ii = 1;ii < nbins_q2-2; ii++){
	    getRat(hF2_x_eAu_temp[ii],hF2_x[ii],1);
	    getRat(hCS_x1_eAu_temp[ii],hCS_x1[ii],1);
	    getRat(hCS_x2_eAu_temp[ii],hCS_x2[ii],1);
	    for(int jj = 1;jj < nbins_x+1; jj++){
		double val1 = hCS_x1_eAu_temp[ii]->GetBinContent(jj);
		double val2 = hCS_x2_eAu_temp[ii]->GetBinContent(jj);
		double val3 = hF2_x_eAu_temp[ii]->GetBinContent(jj);		
		// Doing chi2 calculcation here
		double cen3 = hF2_x_eAu_pseudo[ii]->GetBinContent(jj);
		double err3 = hF2_x_eAu_pseudo[ii]->GetBinError(jj);
		if(val3>0 && cen3>0 && err3>0 && cen3<2){
		    if(err3/val3<0.4){//HERE
			//if(i==1)cout << "debug chi2 " << val3 << " " << cen3 << " " << err3 << " " << (val3-cen3)*(val3-cen3) /err3 /err3 << endl;
			_chi2 += (val3-cen3)*(val3-cen3) /err3 /err3;   
		    }
		}
	    }
	}
	for(int ii = 0;ii < nbins_q2; ii++){
	    delete hCS_x1_eAu_temp[ii];
	    delete hCS_x2_eAu_temp[ii]; 
	    delete hF2_x_eAu_temp[ii];
	}
	chi2->SetBinContent(i+1,_chi2);
    }
    cout << " Done with PDF varying \n" << endl;
    
    double deno = 0;
    for(int i = 0; i < nvary ; i++){
	double _chi2 = chi2->GetBinContent(i+1);
	deno += TMath::Power(_chi2,0.5*(Npts-1)) * TMath::Exp(-_chi2/2.) / nvary;//chi2
	//deno += TMath::Exp(-_chi2/2.) / nvary;//GK
	//cout <<"contribution to weight sum for " << i << " replica: " << TMath::Power(_chi2,0.5*(Npts-1)) * TMath::Exp(-_chi2/2.) << " chi2 " <<  _chi2 << endl;
	cout <<"contribution to weight sum for " << i << " replica: " << TMath::Exp(-_chi2/2.) << " chi2 " <<  _chi2 << " sum " << deno << endl;
    }
    cout <<"calculated sum of weights (deno): " << deno << endl; 
    for(int i = 0; i < nvary ; i++){
	double _chi2 = chi2->GetBinContent(i+1);
	wk[i] = TMath::Power(_chi2,0.5*(Npts-1)) * TMath::Exp(-_chi2/2.) / deno;//chi2
	//wk[i] = TMath::Exp(-_chi2/2.) / deno;//GK
	REP_WEIGHTS->SetBinContent(i+1,wk[i]);
	cout <<"calculated weights for " <<i <<"-th replica: " <<  wk[i] << endl;
    }
    cout << "===================================================================" << endl;
//===================================================== New nPDF Errors on x-sec ====================================================  
    cout << "===================================================================" << endl;
    cout << ">> Starting New PDF variation " << endl;
    iter = 0;
    double eigen_weights[1000];
    for(int i = 1; i < nvary+1 ; i++){
	cout <<">>> On " << i << " out of " << nvary << " variations " << endl;
        double _chi2=0;
        TH1F *hCS_x1_eAu_temp[nbins_q2];
        TH1F *hCS_x2_eAu_temp[nbins_q2];
        TH1F *hF2_x_eAu_temp[nbins_q2];

        for(int ii = 0;ii < nbins_q2; ii++){
	    hCS_x1_eAu_temp[ii] = new TH1F(Form("hCS_x1_%i_eAu_temp",ii),Form("hCS_x1_%i_eAu_temp",ii),nbins_x,binning_x);
            hCS_x2_eAu_temp[ii] = new TH1F(Form("hCS_x2_%i_eAu_temp",ii),Form("hCS_x2_%i_eAu_temp",ii),nbins_x,binning_x);
            hF2_x_eAu_temp[ii] = new TH1F(Form("hF2_x_%i_eAu_temp",ii),Form("hF2_x_%i_eAu_temp",ii),nbins_x,binning_x);
        }
	TH2F* New_N_GLUON = (TH2F*) N_GLUON[i]->Clone("New_N_GLUON");
	eigen_weights[i-1] = wk[i-1];//getNewPDFError(New_N_GLUON,N_GLUON,i,wk,rik,nvary);
	N_GLUON_Weight[iter] = (TH2F*) New_N_GLUON->Clone(Form("N_GLUON_Weight_%i",iter));
	fill("D0Tree_10100_new.root",hCS_x1,hCS_x1_eAu_temp,hq2bin,Weights_10100,Weights_10100_eAu,New_N_GLUON,1,1);
	fill("D0Tree_541_new.root",hCS_x2,hCS_x2_eAu_temp,hq2bin,Weights_541,Weights_541_eAu,New_N_GLUON,1,1);
	fitF2(nbins_q2,nbins_x,binning_q2,binning_x,hF2_x_eAu_temp,hCS_x1_eAu_temp,hCS_x2_eAu_temp);
	for(int ii = 0;ii < nbins_q2-2; ii++){
            getRat(hF2_x_eAu_temp[ii],hF2_x[ii],1);
            getRat(hCS_x1_eAu_temp[ii],hCS_x1[ii],1);
            getRat(hCS_x2_eAu_temp[ii],hCS_x2[ii],1);
            for(int jj = 1;jj < nbins_x+1; jj++){
                double val1 = hCS_x1_eAu_temp[ii]->GetBinContent(jj)-hCS_x1_eAu_ernew[ii]->GetBinContent(jj);
                double val2 = hCS_x2_eAu_temp[ii]->GetBinContent(jj)-hCS_x2_eAu_ernew[ii]->GetBinContent(jj);
                double val3 = hF2_x_eAu_temp[ii]->GetBinContent(jj)-hF2_x_eAu_ernew[ii]->GetBinContent(jj);
		if(val1>0)hCS_x1_eAu_toynew[ii]->Fill(binning_x[jj-1]+0.00000001,wk[i-1]*val1*val1/nvary);
                if(val2>0)hCS_x2_eAu_toynew[ii]->Fill(binning_x[jj-1]+0.00000001,wk[i-1]*val2*val2/nvary);
                if(val3>0)hF2_x_eAu_toynew[ii]->Fill(binning_x[jj-1]+0.00000001,wk[i-1]*val3*val3/nvary);
	    }
        }
	for(int ii = 0;ii < nbins_q2; ii++){
            delete hCS_x1_eAu_temp[ii];
            delete hCS_x2_eAu_temp[ii];
            delete hF2_x_eAu_temp[ii];
        }
	iter++;
    }
    cout << " Done with new PDF varying " << endl;
    for(int ii = 0;ii < nbins_q2-2; ii++){
        for(int jj = 1;jj < nbins_x+1; jj++){
            TH1F* temp1 = (TH1F*)hCS_x1_eAu_toynew[ii]->ProjectionY("temp1",jj,jj);
            TH1F* temp2 = (TH1F*)hCS_x2_eAu_toynew[ii]->ProjectionY("temp2",jj,jj);
            TH1F* temp3 = (TH1F*)hF2_x_eAu_toynew[ii]->ProjectionY("temp3",jj,jj);
            double c1 = hCS_x1_eAu_ernew[ii]->GetBinContent(jj);
            double c2 = hCS_x2_eAu_ernew[ii]->GetBinContent(jj);
            double c3 = hF2_x_eAu_ernew[ii]->GetBinContent(jj);
            double er1 = getError1(temp1,c1,0);
            double er2 = getError1(temp2,c2,0);
            double er3 = getError1(temp3,c3,0);
	    if(temp1->Integral()>0)hCS_x1_eAu_ernew[ii]->SetBinError(jj,er1);
            if(temp2->Integral()>0)hCS_x2_eAu_ernew[ii]->SetBinError(jj,er2);
            if(temp3->Integral()>0)hF2_x_eAu_ernew[ii]->SetBinError(jj,er3);
            delete temp1; delete temp2; delete temp3;
        }
    }
    cout << "===================================================================" << endl;
    for(int i = 1; i < nvary+1 ; i++){
	cout << " >>>>> Eigen weights " << i << " " << eigen_weights[i-1] << endl;
	REP_WEIGHTS->SetBinContent(i,eigen_weights[i-1]);
    }
//============================ This is the drawing part =======================================================
    cout << "===================================================================" << endl;
    cout << "Now entering drawing portion " << endl;
    TH1F* True = (TH1F*)N_GLUON[0]->ProjectionX("True",6,6);
    TH1F* Mod = (TH1F*)N_GLUON_Weight[0]->ProjectionX("Mod",6,6); 
    TH1F* Ini = (TH1F*)N_GLUON_Comb[0]->ProjectionX("Ini",6,6);
    cout << "> Doing gluon part " << endl;
    True->GetYaxis()->SetRangeUser(0,2);
    Ini->SetFillColor(14);
    Ini->SetFillStyle(1000);
    Mod->SetFillColor(kBlue);
    Mod->SetFillStyle(3004);
    Ini->SetMarkerSize(0);
    Mod->SetMarkerSize(0);
    Mod->SetLineColor(kBlue);
    for(int i = 1; i< Ini->GetNbinsX()+1; i++){
	double sum1 = True->GetBinContent(i);
	double sum2 = True->GetBinContent(i);
	double sum21 = 0;
        double sum22 = 0;
	for(int j = 1; j< nvary+1; j++){
            TH1F *_temp2 = (TH1F*)N_GLUON[j]->ProjectionX("_temp2",6,6);
            TH1F *_temp1 = (TH1F*)N_GLUON[j]->ProjectionX("_temp1",6,6);
            sum22+= REP_WEIGHTS->GetBinContent(j)*(sum2 - _temp2->GetBinContent(i))*(sum2 - _temp2->GetBinContent(i));
            sum21+=(sum1 - _temp1->GetBinContent(i))*(sum1 - _temp1->GetBinContent(i));
	    delete _temp1; delete _temp2;
	}
	Mod->SetBinContent(i,sum2);
	Ini->SetBinContent(i,sum1);
	Mod->SetBinError(i,sqrt(sum22/nvary));
        Ini->SetBinError(i,sqrt(sum21/nvary));
    }
    TCanvas *c111 = new TCanvas("c111","Gluon ratio");
    True->Draw("hist same");
    Ini->Draw("E2 same");
    Mod->Draw("E2 same");
    True->Draw("hist same");
    cout << "Done with gluon part" << endl;
    double qq[20] = {1.3,2,3.1,5,8,12,20,31,52,76,120};
    TF1 *line = new TF1("line","1",-100,100);
    line->SetLineStyle(7);
    TLatex lat;
    TLegend *leg = new TLegend(0.17,0.62,0.4,0.87);
    
    TCanvas *c11 = new TCanvas("c11","f2 ratios",2000,800);
    c11->Divide(5,2);
    
    for(int i = 1;i < 4;i++){//nbins_q2-2; i++){
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
	lat.DrawLatex(-3.2,1.5,Form("Q^{2} = %1.f GeV^{2}",qq[i]));
    }
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
    cout << "Done with drawing, now saving if flagged" << endl;
    if(save){
	TFile *fout  = new TFile("results2_NewWeights_Sys.root","RECREATE");
	REP_WEIGHTS->Write("REP_WEIGHTS");
	chi2->Write("chi2");
	for(int i = 0;i<nbins_q2-2;i++){
	    hCS_x2_eAu_er[i]->Write(Form("hCS_x2_eAu_er_%i",i));
	    hCS_x2_eAu_ernew[i]->Write(Form("hCS_x2_eAu_ernew_%i",i));
	    hCS_x2_eAu[i]->Write(Form("hCS_x2_eAu_%i",i));
	    hCS_x1_eAu_er[i]->Write(Form("hCS_x1_eAu_er_%i",i));
	    hCS_x1_eAu_ernew[i]->Write(Form("hCS_x1_eAu_ernew_%i",i));
	    hCS_x1_eAu[i]->Write(Form("hCS_x1_eAu_%i",i));
	    hF2_x_eAu_er[i]->Write(Form("hF2_x_eAu_er_%i",i));
	    hF2_x_eAu_ernew[i]->Write(Form("hF2_x_eAu_ernew_%i",i));
	    hF2_x_eAu[i]->Write(Form("hF2_x_eAu_%i",i));
	    hF2_x_eAu_pseudo[i]->Write(Form("hF2_x_eAu_%i_pseudo",i));
	    hCS_x1_eAu_pseudo[i]->Write(Form("hCS_x1_eAu_%i_pseudo",i));
	    hCS_x2_eAu_pseudo[i]->Write(Form("hCS_x2_eAu_%i_pseudo",i));
	}
	for(int i = 0;i<1000;i++){
	    N_GLUON_Weight[i]->Write(Form("N_GLUON_Weight_%i",i));
	    N_GLUON[i]->Write(Form("N_GLUON_%i",i));
	    N_GLUON_Comb[i]->Write(Form("N_GLUON_Comb_%i",i));
	}
	for(int i = 0;i<nvary;i++){
	    N_GLUON_Replica[i]->Write(Form("N_GLUON_Replica_%i",i));
	}
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
    gr.Fit("fit","Q");
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
    gr.Fit("fit","Q");
    double error = fit.GetParError(0);
    double val=  fit.GetParameter(0);
    return error;
}
void fill(char file1[100],TH1F** hCS_x1,TH1F** hCS_x1_eAu,TH1F* hq2bin,TH2F* Weights,TH2F *Weights_eAu,TH2F* N_GLUON,int vary,double ww){
    
    float mx; float mq2; float mxtp;float mw;
    TChain *tree = new TChain("tree","tree");
    tree->AddFile(file1);
    tree->SetBranchAddress("mq2",&mq2);
    tree->SetBranchAddress("mx",&mx);
    tree->SetBranchAddress("mw",&mw);
    tree->SetBranchAddress("mxtp",&mxtp);
    
    Long64_t iloop = tree->GetEntries()*ww;
    
    for(Long64_t i =0  ;i < iloop; i++){
        //if(i%1000000 == 0) cout << "> On " << i << " out of " << iloop << endl;
	tree->GetEntry(i);
	
        int bin = hq2bin->FindBin(log10(mq2));
	if(bin>12)continue;
        double ww1 = Weights->GetBinContent(Weights->FindBin(log10(mx),log10(mq2)));
	double ww2 = Weights_eAu->GetBinContent(Weights_eAu->FindBin(log10(mx),log10(mq2)));
	double wwn = 1;
	int bin_g = N_GLUON->FindBin(log10(mxtp),log10(mq2));
	wwn = N_GLUON->GetBinContent(bin_g);
	if(wwn<0)wwn=0;
	if(vary==0)
	    hCS_x1[bin-1]->Fill(log10(mx),ww1/ww);
	hCS_x1_eAu[bin-1]->Fill(log10(mx),wwn*ww2*ww);
    }
    delete tree;
}
double getError(TH1F* h, double c, int debug){
    double val=0;int n=0;
    if(c==0)return 0;
    for(int i = 1;i<h->GetNbinsX();i++){
	if(h->GetBinContent(i)==0)continue;
	double vv = h->GetBinCenter(i);
	val += h->GetBinContent(i)*(vv-c)*(vv-c);
	n+=h->GetBinContent(i);
	if(debug==1) cout <<"Err calc " << c << " " << vv << " " << n << endl;
	//cout << " error here " << sqrt(val/n) << " " << h->GetBinCenter(i) << " " << h->GetBinContent(i) << " " << c << endl;
    }

    return sqrt(val/n);
}
double getError1(TH1F* h, double c, int debug){
    double val=0;int n=0;
    if(c==0)return 0;
    for(int i = 1;i<h->GetNbinsX();i++){
        if(h->GetBinContent(i)==0)continue;
        double vv = h->GetBinCenter(i);
        val += h->GetBinContent(i)*(vv);
        n+=h->GetBinContent(i);
        if(debug==1) cout <<"Err calc " << c << " " << vv << " " << n << endl;
        //cout << " error here " << sqrt(val/n) << " " << h->GetBinCenter(i) << " " << h->GetBinContent(i) << " " << c << endl;                                                                                                                      
    }

    return sqrt(val);
}
void getPDFError(TH2F* n, TH2F** _set,int nset){
    // nset should be 97 - 1 for EPPS16 (0 is central value,1-40 for EPS16 erros and 40-97 for CT14NLO erros)                                                                                                                                                            
    int xbins = n->GetNbinsX();
    int ybins = n->GetNbinsY();
    for(int x = 1; x< xbins+1;x++){
        for(int y =1; y< ybins+1-8;y++){
	    double cc = _set[0]->GetBinContent(x,y);
	    // Even nset = S{+} and odd = S{-} and nset should always be odd
	    double val = cc + ((_set[nset+1]->GetBinContent(x,y)-cc)-(_set[nset]->GetBinContent(x,y)-cc))/2.;
	    n->SetBinContent(x,y,val);
        }
    }
}
void getReplica(TH2F* n, TH2F** _set,int nset, double rik[][40], int vary_iter){
    // nset should be 97 - 1 for EPPS16 (0 is central value,1-40 for EPS16 erros and 40-97 for CT14NLO erros)
    int xbins = n->GetNbinsX();
    int ybins = n->GetNbinsY();
    for(int x = 1; x< xbins+1;x++){
       	for(int y =1; y< ybins+1-8;y++){
	    double cc = _set[0]->GetBinContent(x,y);
	    double val = cc;//initialize with central value
	    for(int i = 1; i < nset ; i = i+2){// odd = S{-1} even = S{+1}
		double ran = gRandom->Gaus(0,1);
		rik[vary_iter][i-1] = ran; 
		double _val = 0;
		_val = ((_set[nset+1]->GetBinContent(x,y)-cc)-(_set[nset]->GetBinContent(x,y) - cc))/2. * ran;
		val+=_val;
	    }
	    n->SetBinContent(x,y,val);
	}
    }
}
double getNewPDFError(TH2F* n, TH2F** _set,int nset, double weights[],double rik[][40], int nvary){
    // nset should be 97 - 1 for EPPS16 (0 is central value,1-40 for EPS16 erros and 40-97 for CT14NLO erros)                                                                                                                                             
    int xbins = n->GetNbinsX();
    int ybins = n->GetNbinsY();

    double sum=0;
    for(int rep = 0; rep<nvary; rep++){
	double ran = rik[rep][nset-1];
	double ww = weights[rep];
	if(ran==0 || ww ==0) cout <<"MAJOR ERROR IN NEW PDF-- CHECK!!!" << endl;
	sum+= ran * ww / nvary;
    }
    
    for(int x = 1; x< xbins+1;x++){
        for(int y =1; y< ybins+1-8;y++){
            double cc = _set[0]->GetBinContent(x,y);
            double val = cc;//initialize with central value                                                                                                                                                                                        
	    double _val = ((_set[nset+1]->GetBinContent(x,y)-cc)-(_set[nset]->GetBinContent(x,y) - cc))/2.;
	    _val*=sum;
	    val+=_val;
	    n->SetBinContent(x,y,val);
	}
    }
    return sum;
}
