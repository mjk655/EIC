void test(){

    double deno= 0;
    int nvary = 100;
    for (int i = 1;i<nvary+1;i++){
	double _chi2 = i;
	deno += TMath::Exp(-_chi2/2./13.) / nvary;//GK    
    }
    
    for(int i = 1;i<nvary+1;i++){
        double _chi2 = i;
        double ww= TMath::Exp(-_chi2/2./13.) / deno;//GK                                                                                                                                                                                                                             
	cout << "weight " << ww << endl;
    }



}
