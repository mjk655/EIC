#include "prettycanvasAll.C"
gSystem->Load("libRooFit.so");

using namespace RooFit;

double * fitdata1D(TH1F *data,TH1F* data_l0,TH1F* data_pgf,int save,char dir1[],double xx, double qq)
{
    double fit_v[2];
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
/////Minuit Options
    bool hesse = true;
    bool minos = false;

    const double _start =0;  
    const double _end = 1;
    RooRealVar *DCA = new RooRealVar("DCA","DCA",_start,_end);
    RooDataHist *dataHist = new RooDataHist("dataHist","dataHist",RooArgList(*DCA),data);
    RooDataHist *dataHist_l0 = new RooDataHist("dataHist_l0","dataHist_l0",RooArgList(*DCA),data_l0);
    RooDataHist *dataHist_pgf = new RooDataHist("dataHist_pgf","dataHist_pgf",RooArgList(*DCA),data_pgf);
    
//Signals
    //RooHistPdf g1("g1","g1",RooArgSet(*DCA),*dataHist_l0);
    //RooHistPdf g2("g2","g2",RooArgSet(*DCA),*dataHist_pgf);

    RooRealVar mean1("mean1","mean1",1,0,2);
    RooRealVar sigma1("sigma1","sigma1",0.1,0,0.5);
    RooGaussian g1("g11_1","g1",*DCA,mean1,sigma1);

    RooRealVar a1("sigma1","a1",0.1,-1,1);
    RooPolynomial g2("g2","g2",*DCA,RooArgList(a1),1);


    RooProdPdf l0_signal("l0_signal","l0_signal",g1);
    RooProdPdf pgf_signal("pgf_signal","pgf_signal",g2); 
    
    RooRealVar l0_yield("l0_yield","l0_yield",10000,0,4e9);
    RooRealVar pgf_yield("pgf_yield","pgf_yield",10000,0,4e9);

    RooArgList shapeList;
    RooArgList yieldList;

    shapeList.add(l0_signal);
    shapeList.add(pgf_signal);
        
    yieldList.add(l0_yield);
    yieldList.add(pgf_yield);
//All together now
    RooAddPdf completePDF("completePDF","completePDF",shapeList,yieldList);
   
//Now fit All Data
    RooAbsReal *nllData = completePDF.createNLL(*dataHist,Extended(kTRUE));
    RooAddition * nll = new RooAddition("nll","nll",RooArgSet(*nllData));  
    RooMinuit m (*nll);
    m.setVerbose(kFALSE);
    m.migrad();
    if(hesse)m.hesse();
    if(minos)m.minos();

    RooPlot *DCAFrame = DCA->frame();
    dataHist->plotOn(DCAFrame,Name("HistM"));
    completePDF.plotOn(DCAFrame,Name("DCA_Curve"));
    completePDF.plotOn(DCAFrame,Components(RooArgSet(l0_signal)),LineStyle(8),LineColor(kRed),Name("pi"));
    completePDF.plotOn(DCAFrame,Components(RooArgSet(pgf_signal)),LineStyle(7),LineColor(kGreen-2),FillStyle(1000),FillColor(kGreen-1),Name("b"));
    
    RooArgSet *freeparam = completePDF.getParameters(dataHist);
    int numfreeparam = (freeparam->selectByAttrib("Constant",kFALSE))->getSize();
    TPaveText* txt = new TPaveText(0.4,0.8,0.6,0.9,"BRNDC");
    txt->SetTextSize(0.08);
    txt->SetTextColor(1);
    txt->SetBorderSize(0.1);
    txt->SetFillColor(0);
    txt->SetFillStyle(0);
   TPaveText* txt1 = new TPaveText(0.56,0.55,0.9,0.71,"BRNDC");
    txt1->SetTextSize(0.07);
    txt1->SetTextColor(1);
    txt1->SetBorderSize(0.1);
    txt1->SetFillColor(0);
    txt1->SetFillStyle(0);
    char chi2[50];
    char entr[50];
    char ptbin[50];    
    sprintf(chi2,"#chi^{2}/NDF=%1.2f",DCAFrame->chiSquare("DCA_Curve","HistM",numfreeparam));
    //sprintf(entr,"Signal Yield: %i #pm %i",(int)e_signal_yield.getVal(),(int)TMath::Max(fabs(e_signal_yield.getErrorLo()),e_signal_yield.getErrorHi()));
    if(1){  
        sprintf(ptbin,"log_{10}(x_{B})-Q^{2} (%1.1f,%1.1f)",xx,qq);
	txt->AddText(ptbin);
	DCAFrame->addObject(txt);
    }    
//Objects for saving yields and so forth
    double yields[5];
    double yields_errors[5];
    yields[0]=(int)l0_yield.getVal();
    yields[1]=(int)pgf_yield.getVal();
    yields_errors[0]=l0_yield.getError();
    yields_errors[1]=pgf_yield.getError();
    cout << "====================" << endl;   
    cout << "Chi2 of Mass Fit: " << chi2 << endl;
    cout << "====================" << endl;   
    cout<< "Yields: \n";
    cout << "====================" << endl;   
    cout<<"L0: "<< yields[0] <<" +/- " << yields_errors[0] << " " << endl;
    cout<<"PGF: "    << yields[1] <<" +/- " << yields_errors[1] << " " << endl;
    fit_v[0] = yields[0];
    fit_v[1] = yields_errors[0];

    if(save==1){
	TCanvas *c1 = new TCanvas("c1","full fit");
	prettycanvasAll(c1,DCAFrame,data->GetMaximum());
	c1->SaveAs(Form(dir1));
    }
    return fit_v;
}

