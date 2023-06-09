    #include "Analyzer.hpp"

double generpoiss(double *x, double *p)
{
    // Fit function from Pavel Parygin (CERN)
    double xx = x[0];
    double u0 = p[4];
    double s0 = p[5];
    double g = p[1];
    double sG = p[2];
    double u = p[3];
    double lambda = p[0];
    double N = p[6];
    double npeaks = p[7];
    double fitval = 0.;

    for (int k=0; k<npeaks+1; k++) {
        double uk = u + k*lambda;
        double G = u * pow(uk, k-1) * exp(-uk) / TMath::Factorial(k);
        double sk2 = s0*s0 + k*sG*sG;
        double dx = xx - u0 - k*g;
        fitval += N * G * exp(-0.5*dx*dx/sk2) / sqrt(2*TMath::Pi()) / sqrt(sk2);
    }

    return fitval;
}

// Analyzer::Analyzer(std::string runkey, std::string runlist_filename, int runstart, int runstop, bool stdout_flag):_runkey(runkey),_stdout_flag(stdout_flag),_runstart(runstart),_runstop(runstop)
Analyzer::Analyzer(int argv, char *argc[])
{

    InputArguments input(argv, argc);
    parameters = input.Get();

    runlist_file.open(parameters.runlist_filename);
        
    if (!runlist_file.is_open())
    {
        std::cout << "Error occured when trying to open runlist_file! " << parameters.runlist_filename << std::endl;
        exit(0);
    }            
}

Analyzer::~Analyzer()
{
    runlist_file.close();
}

void Analyzer::Processing()
{
    ParameterParser parser(parameters.argKeys);
    switch (parser.ParseThis(parameters.runkey))
    {
        case 3:
            ProcessingAll();
            break;
    }
    // if (parameters.runkey == "ALL") ProcessingAll();
}

void Analyzer::ProcessingAll()
{
    std::string line, runtype_line;
    std::size_t found;

    while (std::getline(runlist_file, line))
    {
        found = line.find("run");
        if (found != std::string::npos) std::cout << std::stoi(line.substr(4, 4)) << std::endl;
        // {
            // runtype_line = line.substr(found + 5, line.size() - (found + 4));
            // switch (RunTypes[runtype_line])
            // {
            //     case 0:
            //         std::cout << "BVscan processing" << std::endl;
            //         break;
            //     case 1:
            //         std::cout << "LED processing" << std::endl;
            //         break;
            //     case 2:
            //         std::cout << "All batch processing" << std::endl;
            //         break;
            //     default:
            //         std::cout << "Wrong runkey" << std::endl;
            //         break;
            // }

        // }
        else std::cout << "Just skipping this ->" << line << std::endl;
    }
        // 
        // std::size_t found = line.find(key);
        // if (found != std::string::npos) std::cout << std::stoi(line.substr(4, 4)) << std::endl;

}

    // void Analyzer::fitBVscan(int runstart, int runstop)
    // {
    // 	int range = runstop - runstart;
    // 	int k = 0;
    // 	float BVs1[range+1], eBVs1[range+1], Gains1[range+1], eGains1[range+1], BVs2[range+1], eBVs2[range+1], Gains2[range+1], eGains2[range+1], dT1[range+1], dT2[range+1];
    // 	float T1 = 0, T2 = 0;
    // 	//float testBV[10] = {40.50, 40.75, 41., 41.25, 41.5, 41.75, 42.0, 42.25, 42.5, 42.75};

    // 	FinalData output;

    // 	TString str;
    // 	str = "BVscan_run";
    // 	str += runstart;
    // 	str += "_run";
    // 	str += runstop;

    //     std::string str_result = "output_results.txt";
    //     FILE* ptr_res = fopen(str_result.c_str(), "a+");
    //     if (ptr_res == NULL) {
    //         printf("cannot open output file\n");
    //         exit(0);
    //     }

    // 	//TCanvas *c1 = new TCanvas(str+"_ch1",str+"_ch1",800,600);

    // 	for (int i = runstart; i <= runstop; i++)
    // 	{	
    // 		TString filename;
    // 		if (i<1000) filename += "output00000";
    //         	else filename += "output0000";
    // 		filename += i;
    // 		filename += ".root";		
    // 		output = fitBVstep(filename);

    // 		BVs1[k] = output.BV1;
    // 		eBVs1[k] = 0.01;
    // 		Gains1[k] = output.Gain1;
    // 		eGains1[k] = output.eGain1;
    // 		T1 += output.T1;
    //         dT1[k] = output.T1;

    // 		BVs2[k] = output.BV2;
    // 		eBVs2[k] = 0.01;
    // 		Gains2[k] = output.Gain2;
    // 		eGains2[k] = output.eGain2;
    // 		T2 += output.T2;
    //         dT2[k] = output.T2;

    // 		k++;
    // 	}

    // 	TGraphErrors *gr1 = new TGraphErrors(k,BVs1,Gains1,eBVs1,eGains1);
    // 	TF1 *f1 = new TF1("linear_ch1", "[1]*(x-[0])");
    // 	gr1->Fit(f1,"QM0+");
    //     gr1->GetFunction("linear_ch1")->ResetBit(TF1::kNotDraw);


    // 	TGraphErrors *gr2 = new TGraphErrors(k,BVs2,Gains2,eBVs2,eGains2);
    // 	TF1 *f2 = new TF1("linear_ch2", "[1]*(x-[0])");
    // 	gr2->Fit(f2,"QM0+");
    //     gr2->GetFunction("linear_ch2")->ResetBit(TF1::kNotDraw);


    // 	//std::cout << f1->GetParameter(0) << "," << f1->GetParError(0) << "," << f2->GetParameter(0) << "," << f2->GetParError(0) << "," << f1->GetParameter(1) << "," << f1->GetParError(1) << "," << f2->GetParameter(1) << "," << f2->GetParError(1) << std::endl;

    // 	fprintf(ptr_res,"%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",runstart,runstop,f1->GetParameter(0),f1->GetParError(0),f2->GetParameter(0),f2->GetParError(0),f1->GetParameter(1),f1->GetParError(1),f2->GetParameter(1),f2->GetParError(1), 1.0*T1/k, 1.0*T2/k, dT1[0]-dT1[range], dT2[0]-dT2[range]);

    // 	delete gr1;
    // 	delete gr2;
    // 	delete f1;
    // 	delete f2;

    //     fclose(ptr_res);

    // }


FinalData Analyzer::fitBVstep(const char* fileName)
{
    TString fileName_full = "/home/qfl/online/midas_digi/daq/";
    fileName_full += fileName;
    //cout << "Analyzing file " << fileName << endl;sssss
    TFile *file = TFile::Open(fileName_full);
    TH1F *htemp1, *htemp2, *hBV1, *hBV2, *hI1, *hI2, *hT1, *hT2, *hB1, *hB2, *hL3, *hLb3;;// *hQ4
    TH2F *hmean_rms1, *hmean_rms2;
    TProfile *hmean1, *hmean2, *hrms1, *hrms2;

    gErrorIgnoreLevel = 6001;

    file->GetObject("ch1/charge_spectrum_cut_ch1",htemp1);//add cut!!!
    file->GetObject("ch1/waveformMEAN_ch1",hmean1);
    file->GetObject("ch1/waveformRMS_ch1",hrms1);
    file->GetObject("ch1/Mean_vs_RMS_ch1",hmean_rms1);
    file->GetObject("ch1/baseline_spectrum_ch1",hB1);
    file->GetObject("gpib/gpib_0",hI1);
    file->GetObject("gpib/gpib_2",hBV1);
    file->GetObject("temperatures/temperature_0",hT1);
    TF1 *ch1 = GPFit(htemp1,0);

    file->GetObject("ch2/charge_spectrum_cut_ch2",htemp2);//add cut!!!
    file->GetObject("ch2/waveformMEAN_ch2",hmean2);
    file->GetObject("ch2/waveformRMS_ch2",hrms2);
    file->GetObject("ch2/Mean_vs_RMS_ch2",hmean_rms2);
    file->GetObject("ch2/baseline_spectrum_ch2",hB2);
    file->GetObject("gpib/gpib_1",hI2);
    file->GetObject("gpib/gpib_3",hBV2);
    file->GetObject("temperatures/temperature_1",hT2);
    TF1 *ch2 = GPFit(htemp2,0);

    //   file->GetObject("ch4/charge_spectrum_ch4",hQ4);
    //   TF1 *ch4 = GFit(hQ4, 0);

    FinalData output;

    file->GetObject("ch3/charge_spectrum_ch3",hL3);
    TF1 *ch3 = GFit(hL3, 0, 256, 1, 3);
    output.LYSO_Yield = ch3->GetParameter(1);

    file->GetObject("ch3/charge_spectrum_cut_ch3",hLb3);
    if (hLb3->GetEntries()!=0) {
        TF1 *ch3b = DGFit(hLb3, 0, 4);
        output.LYSO_Gain = ch3b->GetParameter(1) - ch3b->GetParameter(0);
    }
    else {
        output.LYSO_Gain = 1;
    }  

    output.Gain1 = ch1->GetParameter(1);
    output.eGain1 = ch1->GetParError(1);

    output.Gain2 = ch2->GetParameter(1);
    output.eGain2 = ch2->GetParError(1);

    output.BV1 = hBV1->Integral()/hBV1->GetEntries();
    output.BV2 = hBV2->Integral()/hBV2->GetEntries();

    output.current1 = hI1->Integral()/hI1->GetEntries();
    output.current2 = hI2->Integral()/hI2->GetEntries();

    output.T1 = hT1->Integral()/hT1->GetEntries();
    output.T2 = hT2->Integral()/hT2->GetEntries();

    //cout << output.current1 << "\t" << output.current2 << endl;

    output.mean1 = hmean1->GetSum()/hmean1->GetEntries();
    output.mean2 = hmean2->GetSum()/hmean2->GetEntries();

    output.rms1 = hrms1->GetSum()/hrms1->GetEntries();
    output.rms2 = hrms2->GetSum()/hrms2->GetEntries();

    //hmean_rms1->GetXaxis()->SetRangeUser(-100,10);
    //output.baseline1 = hmean_rms1->GetMean(2);
    //   output.rmsb1 = hmean_rms1->GetRMS(2);
    TF1 *ch1b = GFit(hB1, 0, 1, 3, 1);
    //output.baseline1 = hB1->GetXaxis()->GetBinCenter(hB1->GetMaximumBin());
    output.baseline1 = ch1b->GetParameter(1);
    output.rmsb1 = hB1->GetRMS();

    //hmean_rms2->GetXaxis()->SetRangeUser(-100,10);
    //output.baseline2 = hmean_rms2->GetMean(2);
    //output.rmsb2 = hmean_rms2->GetRMS(2);
    //output.baseline2 = hB2->GetMean();
    TF1 *ch2b = GFit(hB2, 0, 1, 3, 1);
    //output.baseline2 = hB2->GetXaxis()->GetBinCenter(hB2->GetMaximumBin());
    output.baseline2 = ch2b->GetParameter(1);
    output.rmsb2 = hB2->GetRMS();
    
    delete htemp1;
    delete htemp2;
    delete hBV1;
    delete hBV2;
    delete hI1;
    delete hI2;
    delete hrms1;
    delete hrms2;
    delete hmean1;
    delete hmean2;
    delete hT1;
    delete hT2;
    delete ch1;
    delete ch2;
    delete hmean_rms1;
    delete hmean_rms2;

    //delete hQ4;

    file->Close();

    delete file;

    return output;

}

void Analyzer::fitGain(int run)
{
    std::string str_result = "output_results.txt";
    FILE* ptr_res = fopen(str_result.c_str(), "a+");
    if (ptr_res == NULL) {
        printf("cannot open output file\n");
        exit(0);
    }

	FinalData o;
	TString filename;
	if (run<1000) filename += "output00000";
	else filename += "output0000";
	filename += run;
	filename += ".root";		
	o = fitBVstep(filename);

	fprintf(ptr_res,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", run,
                                                                            o.Gain1, o.Gain2,
                                                                            o.eGain1, o.eGain2,
                                                                            o.mean1, o.mean2,
                                                                            o.baseline1, o.baseline2,
                                                                            o.rms1, o.rms2,
                                                                            o.rmsb1, o.rmsb2,
                                                                            o.current1*1e9, o.current2*1e9,
                                                                            o.T1, o.T2,
                                                                            o.LYSO_Yield, o.LYSO_Gain);

    fclose(ptr_res);
}


TF1 *Analyzer::GPFit(TH1F* hist2proc, bool verbose, int xlow, int xhigh){

    // TCanvas c1("c1", "c1")
    double npe_peaks, mu_hist, par[6];
    double binwidth;
    int ntotal, n0;

    hist2proc->Rebin(20);
    //hist2proc->GetXaxis()->UnZoom();
    binwidth=hist2proc->GetBinWidth(1);
    // uncomment this for day 7
    xlow = hist2proc->GetBinCenter(hist2proc->FindFirstBinAbove(0,1));
    xhigh = hist2proc->GetBinCenter(hist2proc->FindLastBinAbove(0,1)); // in case of only GP in histogram
    hist2proc -> GetXaxis() -> SetRangeUser(xlow, xhigh);
    ntotal=hist2proc->Integral();
    if (verbose) std::cout<<"Spectrum created"<<std::endl;

    // Find all peaks and write their positions to vector
    TSpectrum s(15);
    int nfound = s.Search(hist2proc, 3, "goff", 0.05);
    double *_xpeaks = s.GetPositionX();
    std::vector<float> xpeaks;
    xpeaks.clear();
    if (_xpeaks==nullptr) {
        std::cout << "_xpeaks is nullptr" << std::endl;
    } else {
        for (int p = 0; p < nfound; p++) xpeaks.push_back(_xpeaks[p]);
    }
    std::sort(xpeaks.begin(), xpeaks.end());

    if (verbose) std::cout << "Gauss fiting" << std::endl;
    // Fit first two peaks with Gauss function. Fit range is 1/3 of distance between two peaks
    TF1 g1("peak0","gaus",xpeaks[0]-(xpeaks[1]-xpeaks[0])/3,xpeaks[0]+(xpeaks[1]-xpeaks[0])/3);
    g1.SetLineColor(kCyan);
    hist2proc -> Fit(&g1,"RQ0");
    hist2proc->GetFunction("peak0")->ResetBit(TF1::kNotDraw);

    TF1 g2("peak1","gaus",xpeaks[1]-(xpeaks[1]-xpeaks[0])/3,xpeaks[1]+(xpeaks[1]-xpeaks[0])/3);
    g2.SetLineColor(kCyan);
    hist2proc -> Fit(&g2,"RQ0+");
    hist2proc->GetFunction("peak1")->ResetBit(TF1::kNotDraw);

    fp->g1Chi2 = g2.GetChisquare();
    fp->g0NDF = g1.GetNDF();
    fp->g1NDF = g2.GetNDF();
    fp->g0Integral = g1.Integral(xpeaks[0]-(xpeaks[1]-xpeaks[0])/3,xpeaks[0]+(xpeaks[1]-xpeaks[0])/3);
    fp->g1Integral = g2.Integral(xpeaks[1]-(xpeaks[1]-xpeaks[0])/3,xpeaks[1]+(xpeaks[1]-xpeaks[0])/3);


    g1.SetLineColor(kGreen+2);
    g2.SetLineColor(kGreen+2);
    g1.GetParameters(&par[0]);
    g2.GetParameters(&par[3]);

    // Get mean and sigma of two Gaussians
    fp->gN0 = g1.GetParameter(0); // g1params[0]
    fp->gm0 = g1.GetParameter(1);
    fp->gs0 = g1.GetParameter(2);
    fp->gN1 = g2.GetParameter(0);
    fp->gm1 = g2.GetParameter(1);
    fp->gs1 = g2.GetParameter(2);
    fp->gain_spe_01 = fp->gm1 - fp->gm0;

    // Set new range for histogram (Mu1-Mu0)/2 around Mu0
    hist2proc->GetXaxis()->SetRangeUser(fp->gm0-(fp->gm1-fp->gm0)/2,fp->gm0+(fp->gm1-fp->gm0)/2);
    n0 = hist2proc -> Integral();
    hist2proc->GetXaxis()->SetRangeUser(xlow, xhigh);
    mu_hist= -TMath::Log(1.*n0/ntotal);
    //npe_peaks = 3;
    npe_peaks = TMath::Nint((xhigh-fp->gm0)/(fp->gm1-fp->gm0));
    //cout << npe_peaks << endl;

    // Construct GP function in range [Mu0-(Mu1-Mu0)/2, xhigh]
    TF1 *func = new TF1("generpoiss", generpoiss, fp->gm0-(fp->gm1-fp->gm0)/2, xhigh, 8);
    // Set initial parameters for GP and limits for them
    func -> SetParameters( 0.1, fp->gm1-fp->gm0, sqrt(abs(fp->gs1*fp->gs1 - fp->gs0*fp->gs0)), mu_hist, fp->gm0, fp->gs0, binwidth*ntotal, npe_peaks); //initial parameters
    func -> SetParLimits(1, 0, 1.5*(fp->gm1-fp->gm0));
    //func -> FixParameter(0, 0);
    //func -> FixParameter(7, npe_peaks);
    func -> FixParameter(7, 5);
    func -> SetParNames("lambda", "gain_SPE", "sigma_gain", "mu_avg", "mu0", "sigma0", "NxW", "npe_peaks");
    func -> SetNpx(1000);
    if (verbose) std::cout<<"GP fiting"<<std::endl;
    hist2proc -> Fit(func, "RQ0+");
    hist2proc->GetFunction("generpoiss")->ResetBit(TF1::kNotDraw);
    //hist2proc->GetXaxis()->UnZoom();

    // put parameters of GP to fp structure
    fp->gain = func->GetParameter(1);
    fp->u0 = func->GetParameter(4);
    fp->s0 = func->GetParameter(5);
    fp->sG = func->GetParameter(2);
    fp->lambda = func->GetParameter(0);
    fp->navg = func->GetParameter(3);
    fp->npeaks = func->GetParameter(7);

    fp->gain_err = func->GetParError(1);

    double Chi2 = func->GetChisquare();
    double NDF = func->GetNDF();

    fp->chi2_ndf = Chi2/NDF;

    if (verbose) {

        std::cout << std::endl;
        std::cout << "GenPoiss Fit is Done. Parameters:" << std::endl;
        std::cout << "Baseline Sigma0 = " << fp->s0 << std::endl;
        std::cout << "SPE Gain = " << fp->gain << " +- " << fp->gain_err << std::endl;
        std::cout << "SPE Gain Width = " << fp->sG << std::endl;
        std::cout << "SPE Lambda = " << fp->lambda << std::endl;
        std::cout << "SPE Mu = " << fp->navg << std::endl;
    
        std::cout << "Chi^2/NDF is " << Chi2/NDF << std::endl;
        std::cout << std::endl;
    }

    // return GP function
    return (func);
}

TF1 *Analyzer::GFit(TH1F *hist, bool verbose, int rebin, int multA, int multB) {

    hist->Rebin(rebin);

    double lower_bound = hist->GetXaxis()->GetBinLowEdge(hist->GetMaximumBin()-5);
    double upper_bound = hist->GetXaxis()->GetBinUpEdge(hist->GetMaximumBin()+5);
    int gnpar = 3, lgnpar = 4;
    double gparameters[gnpar], parameters[lgnpar];

    if (verbose) {
        std::cout << "Starting 1 Gaussian Fit...." << std::endl;
        std::cout << std::endl;
    }

    TF1 *G_func = new TF1("1 fit","gaus",lower_bound, upper_bound);
    hist->Fit(G_func,"RQ0+");
    hist->GetFunction("1 fit")->ResetBit(TF1::kNotDraw);
    G_func->GetParameters(gparameters);

    double Chi2 = G_func->GetChisquare();
    double NDF = G_func->GetNDF();

    parameters[0] = gparameters[0];
    parameters[1] = gparameters[1];
    parameters[2] = gparameters[2];

    lower_bound = parameters[1]-parameters[2]/3*multA;
    upper_bound = parameters[1]+parameters[2]/3*multB;

    if (verbose) {
        std::cout << "Starting 2 Gaussian Fit...." << std::endl;
        std::cout << std::endl;
    }

    TF1 *LG_func = new TF1("2 fit","gaus",lower_bound, upper_bound);
    LG_func->SetParameters(parameters);
    LG_func->SetLineColor(kGreen);
    hist->Fit(LG_func,"RQ0+");
    hist->GetFunction("2 fit")->ResetBit(TF1::kNotDraw);
    LG_func->GetParameters(parameters);

    if (verbose) {
        std::cout << "Gauss Fit is Done. Parameters:" << std::endl;
        std::cout << "Maximum Gauss = " << parameters[0] << std::endl;
        std::cout << "Mean Gauss = " << parameters[1] << std::endl;
        std::cout << "Sigma Gauss = " << parameters[2] << std::endl;
    
        //Chi2 = LG_func->GetChisquare();
        //NDF = LG_func->GetNDF();

        //std::cout << std::endl;
        //std::cout << "Chi^2/NDF is " << Chi2 << " / " << NDF << std::endl;
        std::cout << "Chi^2/NDF is " << Chi2/NDF << std::endl;
        std::cout << std::endl;

    }

    return LG_func;
}

TF1 *Analyzer::DGFit(TH1F *hist, bool verbose, int rebin) {

    hist->Rebin(rebin);

    int gnpar = 3;
    double gparameters1[gnpar], gparameters2[gnpar];

    TSpectrum s(2);
    int nfound = s.Search(hist, 3, "goff", 0.05);
    double *_xpeaks = s.GetPositionX();
    std::vector<float> xpeaks;
    xpeaks.clear();
    if (_xpeaks==nullptr) {
        std::cout << "_xpeaks is nullptr" << std::endl;
    } else {
        for (int p = 0; p < nfound; p++) xpeaks.push_back(_xpeaks[p]);
    }
    std::sort(xpeaks.begin(), xpeaks.end());

    if (verbose) std::cout << "Double Gauss fiting" << std::endl;

    TF1 *g1 = new TF1("peak11","gaus",xpeaks[0]-(xpeaks[1]-xpeaks[0])/3,xpeaks[0]+(xpeaks[1]-xpeaks[0])/3);
    g1->SetLineColor(kCyan);
    g1->SetParameter(1,xpeaks[0]);
    hist->Fit(g1,"RQ0");
    g1->GetParameters(gparameters1);
    hist->GetFunction("peak11")->ResetBit(TF1::kNotDraw);

    TF1 *g2 = new TF1("peak12","gaus",xpeaks[1]-(xpeaks[1]-xpeaks[0])/3,xpeaks[1]+(xpeaks[1]-xpeaks[0])/3);
    g2->SetParameter(1,xpeaks[1]);
    g2->SetLineColor(kCyan);
    hist->Fit(g2,"RQ0+");
    g2->GetParameters(gparameters2);
    hist->GetFunction("peak12")->ResetBit(TF1::kNotDraw);


    if (verbose) {

        double Chi2 = g1->GetChisquare();
        double NDF = g1->GetNDF();

        std::cout << "1 st Gauss Peak Fit is Done. Parameters:" << std::endl;
        std::cout << "Maximum 1st Gauss = " << gparameters1[0] << std::endl;
        std::cout << "Mean 1st Gauss = " << gparameters1[1] << std::endl;
        std::cout << "Sigma 1st Gauss = " << gparameters1[2] << std::endl;
        std::cout << "Chi^2/NDF is " << Chi2/NDF << std::endl;
        std::cout << std::endl;


        Chi2 = g2->GetChisquare();
        NDF = g2->GetNDF();

        std::cout << "2 st Gauss Peak Fit is Done. Parameters:" << std::endl;
        std::cout << "Maximum 2st Gauss = " << gparameters2[0] << std::endl;
        std::cout << "Mean 2st Gauss = " << gparameters2[1] << std::endl;
        std::cout << "Sigma 2st Gauss = " << gparameters2[2] << std::endl;
        std::cout << "Chi^2/NDF is " << Chi2/NDF << std::endl;
        std::cout << std::endl;

    }

    g2->SetParameter(0,gparameters1[1]);
    g2->SetParError(0,g1->GetParError(1)); 

    return g2;

}