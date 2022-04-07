void prettycanvasAllleg(TCanvas *c1,RooPlot *Frame, RooPlot *Pull, int log, double max)
{
    TPaveText *Name = new TPaveText(gStyle->GetPadLeftMargin() + 0.06,
                            0.83 - gStyle->GetPadTopMargin(),
                            gStyle->GetPadLeftMargin() + 0.3,
                            0.99 - gStyle->GetPadTopMargin(),
                                        "BRNDC");

//    Name->AddText("        #it{STAR}");
    Name->AddText("#color[2]{#it{STAR}}");// Au+Au #sqrt{s_{NN}} = 200 GeV");
    //Name->AddText("STAR Preliminary");
    Name->SetTextSize(0.08);
    Name->SetFillColorAlpha(0, 0);
    Name->SetTextAlign(12);
    Name->SetBorderSize(0); 

    TPaveText *Name3 = new TPaveText(0.17,0.8,0.2,0.89,"BRNDC");
    Name3->AddText("(b)");
    Name3->SetFillColorAlpha(0, 0);
    Name3->SetTextAlign(12);
    Name3->SetBorderSize(0);
    Name3->SetTextSize(0.095);

    TPaveText *Name1 = new TPaveText(gStyle->GetPadLeftMargin() + 0.27,
				     0.72,
				    gStyle->GetPadLeftMargin() + 0.79,
				     0.87,
				    "BRNDC");
    TPaveText *Name2 = new TPaveText(gStyle->GetPadLeftMargin() + 0.55,//0.27,                                                                                                                                                                                     
                                     0.68,
				     gStyle->GetPadLeftMargin() + 0.69,
                                     0.77,
				     "BRNDC");

//    Name->AddText("        #it{STAR}");                                                                                                                                            
    Name1->AddText("            Au+Au #sqrt{s_{NN}} = 200 GeV");
    Name1->AddText(" ");
    Name1->AddText("                                         0-80%");
    Name1->SetTextSize(0.08);
    Name1->SetFillColorAlpha(0, 0);
    Name1->SetTextAlign(12);
    Name1->SetBorderSize(0);
    Name2->SetTextSize(0.08);
    Name2->SetFillColorAlpha(0, 0);
    Name2->SetTextAlign(12);
    Name2->SetBorderSize(0);
    TPad *pad1 = new TPad("pad1","pad1",0.0,0.25,1.0,1.0);
    TPad *pad2 = new TPad("pad2","pad2",0.0,0.0,1.0,0.255);
    TGaxis::SetMaxDigits(3);
    pad1->SetBottomMargin(0.025);
    //pad1->SetTopMargin(0.08);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.5);  
    pad2->SetBorderMode(0);
    pad1->SetBorderMode(0);
    c1->SetBorderMode(0);  
    pad1->Draw();
    pad2->Draw(); 
    pad1->cd(); 
    Frame->GetYaxis()->SetTitleSize(0.1); 
    Frame->GetYaxis()->SetTitleOffset(0.6);
    Frame->SetTitle(""); 
    Frame->GetXaxis()->SetTitle("log(DCA/cm)");
    Frame->GetYaxis()->SetTitle("Candidates / (0.04)");
      Frame->GetYaxis()->CenterTitle();
    Frame->GetYaxis()->SetRangeUser(0.2,max*10);
    Pull->GetXaxis()->SetTitle("log(DCA/cm)"); 
    //Frame->GetXaxis()->SetTitleSize(0.3);
    //Frame->GetXaxis()->SetLabelSize(0.1);  
   
    //Frame->GetYaxis()->SetNdivisions(4,0); 
    Frame->Draw();
    leg = new TLegend(0.17,0.2,0.5,0.8);
    leg->SetFillStyle(0);
    TH1F *temp1 = new TH1F("temp1","temp1",2,0,1);
    TH1F *temp11 = new TH1F("temp11","temp11",2,0,1);
    TH1F *temp2= new TH1F("temp2","temp2",2,0,1);
    TH1F *temp3=new TH1F("temp3","temp3",2,0,1);
    TH1F *temp4=new TH1F("temp4","temp4",2,0,1);
    TH1F *temp5=new TH1F("temp5","temp5",2,0,1);
    temp1->SetLineColor(1);//kOrange-3);
    
    temp11->SetLineColor(kCyan);
    temp2->SetLineColor(kRed);
    temp3->SetLineColor(kMagenta);
    temp4->SetLineColor(kGreen-2);
    temp5->SetLineColor(kBlue);
    temp11->SetLineStyle(4);
    temp1->SetLineStyle(2);
    temp2->SetLineStyle(8);
    temp3->SetLineStyle(9);
    temp4->SetLineStyle(7);
    temp5->SetLineStyle(1);



    leg->AddEntry(temp5,"Combined template","l");
    leg->AddEntry(temp1,"Photonic #it{e}","l");
    //leg->AddEntry(temp11,"Photonic #it{e} mis-matched","l");
    leg->AddEntry(temp2,"Hadrons","l");
    leg->AddEntry(temp3,"#it{c}#rightarrow#it{e}","l");
    // leg->AddEntry(temp5,"#it{D^{+}}#rightarrow#it{e}","l");
    leg->AddEntry(temp4,"#it{b}#rightarrow#it{e}","l");
//    leg->AddEntry(temp5,"Combined template","l");
//    Name->Draw("SAME");
    Name1->Draw("SAME");
    Name3->Draw("SAME");
    leg->Draw("same");

    gPad->RedrawAxis();
    if(log)pad1->SetLogy();
    pad2->cd();
    Pull->SetTitle(0);
    Pull->GetYaxis()->SetTitle("Pull");
    Pull->GetYaxis()->SetTitleSize(0.2);
    Pull->GetYaxis()->SetLabelSize(0.15);
    Pull->GetYaxis()->CenterTitle();
    Pull->GetXaxis()->CenterTitle();
    Pull->GetYaxis()->SetTitleOffset(0.3);
    Pull->GetYaxis()->SetNdivisions(505);
    Pull->GetYaxis()->SetRangeUser(-5.9,5.9);
    Pull->GetXaxis()->SetTitleSize(0.3);
    Pull->GetXaxis()->SetLabelSize(0.2);  
    Pull->GetXaxis()->SetTitleOffset(0.7);
    Pull->Draw();
    //temp1->Delete();
    // temp2->Delete();
    //temp3->Delete();
    //temp4->Delete();
}void prettycanvasAll(TCanvas *c1,RooPlot *Frame, double max)
{
    TH1F *temp1 = new TH1F("temp1","temp1",2,0,1);
    TH1F *temp2= new TH1F("temp2","temp2",2,0,1);
    temp1->SetLineColor(kRed);
    temp2->SetLineColor(kGreen-2);
    temp1->SetLineStyle(8);
    temp2->SetLineStyle(7);
    leg = new TLegend(0.2,0.6,0.5,0.75);
    leg->AddEntry(temp1,"#gamma*+#it{c}#rightarrow#it{c}","L");
    leg->AddEntry(temp2,"PGF","L");

    Frame->GetXaxis()->SetTitle("|#phi(D^{0}) - #phi(e')|/#pi");
    Frame->GetYaxis()->SetTitle("dN/d#Delta#phi");
    Frame->GetYaxis()->SetRangeUser(0,max*1.3);
    Frame->Draw();
    leg->Draw("same");
}
void prettycanvasAll1(TCanvas *c1,RooPlot *Frame, RooPlot *Pull, int log, double max)
{
    TPaveText *Name = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
                            0.75 - gStyle->GetPadTopMargin(),
                            gStyle->GetPadLeftMargin() + 0.31,
                            0.95 - gStyle->GetPadTopMargin(),
                                        "BRNDC");

    Name->AddText("");
    //Name->SetFillColorAlpha(0, 0);
    Name->SetTextAlign(12);
    Name->SetBorderSize(0); 
    TPad *pad1 = new TPad("pad1","pad1",0.0,0.,1.0,1.0);
    TPad *pad2 = new TPad("pad2","pad2",0.0,0.0,1.0,0.255);
    TGaxis::SetMaxDigits(3);
    //pad1->SetBottomMargin(0.02);
    pad1->SetTopMargin(0.08);
    pad2->SetTopMargin(0.015);
    pad2->SetBottomMargin(0.5);  
    pad2->SetBorderMode(0);
    pad1->SetBorderMode(0);
    c1->SetBorderMode(0);  
    pad1->Draw();
    //pad2->Draw(); 
    pad1->cd(); 
    Frame->GetYaxis()->SetTitleSize(0.1); 
    Frame->GetYaxis()->SetTitleOffset(0.6);
    Frame->SetTitle(""); 
    Frame->GetXaxis()->SetTitle("log_{10}(|DCA|/cm)");
    Pull->GetXaxis()->SetTitle("log_{10}(|DCA|/cm)"); 
    //Frame->GetXaxis()->SetTitleSize(0.3);
    //Frame->GetXaxis()->SetLabelSize(0.1);  
   
    //Frame->GetYaxis()->SetNdivisions(4,0); 
    Frame->Draw();
    //Name->Draw("SAME");
    if(log)pad1->SetLogy();
    //pad2->cd();
    Pull->SetTitle(0);
    Pull->GetYaxis()->SetTitle("Pull");
    Pull->GetYaxis()->SetTitleSize(0.15);
    Pull->GetYaxis()->SetLabelSize(0.1);
    Pull->GetYaxis()->CenterTitle();
    Pull->GetYaxis()->SetTitleOffset(0.3);
    Pull->GetYaxis()->SetNdivisions(505);
    Pull->GetYaxis()->SetRangeUser(-9.9,9.9);
    Pull->GetXaxis()->SetTitleSize(0.3);
    Pull->GetXaxis()->SetLabelSize(0.15);  
    Pull->GetXaxis()->SetTitleOffset(0.63);
    ///Pull->Draw();
}
void prettycanvasAll2(TCanvas *c1,RooPlot *Frame, RooPlot *Pull, int log)
{
    TPaveText *Name = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
                            0.75 - gStyle->GetPadTopMargin(),
                            gStyle->GetPadLeftMargin() + 0.31,
                            0.95 - gStyle->GetPadTopMargin(),
                                        "BRNDC");

    Name->AddText("");
    //Name->SetFillColorAlpha(0, 0);
    Name->SetTextAlign(12);
    Name->SetBorderSize(0); 
    TPad *pad1 = new TPad("pad1","pad1",0.0,0.25,1.0,1.0);
    TPad *pad2 = new TPad("pad2","pad2",0.0,0.0,1.0,0.255);
    TGaxis::SetMaxDigits(3);
    pad1->SetBottomMargin(0.02);
    pad1->SetTopMargin(0.08);
    pad2->SetTopMargin(0.015);
    pad2->SetBottomMargin(0.5);  
    pad2->SetBorderMode(0);
    pad1->SetBorderMode(0);
    c1->SetBorderMode(0);  
    pad1->Draw();
    pad2->Draw(); 
    pad1->cd(); 
    Frame->GetYaxis()->SetTitleSize(0.1); 
    Frame->GetYaxis()->SetTitleOffset(0.6);
    Frame->SetTitle(""); 
    Frame->GetXaxis()->SetTitle("n#sigma_{e}");
    Pull->GetXaxis()->SetTitle("n#sigma_{e}"); 
    //Frame->GetXaxis()->SetTitleSize(0.3);
    //Frame->GetXaxis()->SetLabelSize(0.1);  
   
    //Frame->GetYaxis()->SetNdivisions(4,0); 
    Frame->Draw();
    //Name->Draw("SAME");
    if(log)pad1->SetLogy();
    pad2->cd();
    Pull->SetTitle(0);
    Pull->GetYaxis()->SetTitle("Pull");
    Pull->GetYaxis()->SetTitleSize(0.15);
    Pull->GetYaxis()->SetLabelSize(0.1);
    Pull->GetYaxis()->CenterTitle();
    Pull->GetYaxis()->SetTitleOffset(0.3);
    Pull->GetYaxis()->SetNdivisions(505);
    Pull->GetYaxis()->SetRangeUser(-9.9,9.9);
    Pull->GetXaxis()->SetTitleSize(0.3);
    Pull->GetXaxis()->SetLabelSize(0.15);  
    Pull->GetXaxis()->SetTitleOffset(0.63);
    Pull->Draw();
}
