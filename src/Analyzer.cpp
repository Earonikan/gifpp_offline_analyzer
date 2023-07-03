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

Analyzer::Analyzer(int argv, char *argc[])
{

    InputArguments input(argv, argc);
    parameters = input.Get();

    runlist_file.open(parameters.runlist_filename);

    if (parameters.root_tree) {root_file = new TFile("GIFpp_Data.root","recreate");}
        
    if (!runlist_file.is_open())
    {
        std::cout << "Error occured when trying to open runlist_file! " << parameters.runlist_filename << std::endl;
        exit(0);
    }            
}

Analyzer::~Analyzer()
{
	if (parameters.verbose > 0) std::cout << "Closing Analyzer!" << std::endl;
    runlist_file.close();
    if (parameters.root_tree)
    {
        tree->Write();
        tree_bv->Write();
        root_file->Write();
        root_file->Close();
        // std::cout << "HERE" << std::endl;
    }
}

void Analyzer::Processing()
{
    ParameterParser parser(parameters.runTypes);
    //std::cout << parser.ParseThis(parameters.runkey) << std::endl;
    if (parameters.root_tree)
    {
        //monitoring runs tree
        tree = new TTree("GIFpp_DATA", "GIFpp data tree");
		
		tree->Branch("run", &output.run);
        tree->Branch("runtype", &output.runtype);
        tree->Branch("time_stamp", &output.timestamp);
        tree->Branch("Gain1", &output.Gain1);
        tree->Branch("Gain2", &output.Gain2);
        tree->Branch("eGain1", &output.eGain1);
        tree->Branch("eGain2", &output.eGain2);
        tree->Branch("mean1", &output.mean1);
        tree->Branch("mean2", &output.mean2);
        tree->Branch("baseline1", &output.baseline1);
        tree->Branch("baseline2", &output.baseline2);
        tree->Branch("rms1", &output.rms1);
        tree->Branch("rms2", &output.rms2);
        tree->Branch("rmsb1", &output.rmsb1);
        tree->Branch("rmsb2", &output.rmsb2);
        tree->Branch("current1", &output.current1);
        tree->Branch("current2", &output.current2);
        tree->Branch("BV1", &output.BV1);
        tree->Branch("BV2", &output.BV2);
        tree->Branch("T1", &output.T1);
        tree->Branch("T2", &output.T2);
        tree->Branch("LYSO_Yield", &output.LYSO_Yield);
        tree->Branch("LYSO_Gain", &output.LYSO_Gain);
        tree->Branch("dip_15402", &output.dip_15402);
        tree->Branch("dip_filter", &output.dip_filter);


        //BV tree
        tree_bv = new TTree("GIFpp_DATA_BVSCANs", "GIFpp data BV tree");
		
		tree_bv->Branch("runstart", &output_bv.runstart);
        tree_bv->Branch("runstop", &output_bv.runstop);
        tree_bv->Branch("time_stamp", &output_bv.timestamp);
        tree_bv->Branch("BDV1", &output_bv.BDV1);
        tree_bv->Branch("BDV2", &output_bv.BDV2);
        tree_bv->Branch("eBDV1", &output_bv.eBDV1);
        tree_bv->Branch("eBDV2", &output_bv.eBDV2);
        tree_bv->Branch("Slope1", &output_bv.Slope1);
        tree_bv->Branch("Slope2", &output_bv.Slope2);
        tree_bv->Branch("eSlope1", &output_bv.eSlope1);
        tree_bv->Branch("eSlope2", &output_bv.eSlope2);
        tree_bv->Branch("T1", &output_bv.T1);
        tree_bv->Branch("T2", &output_bv.T2);
        tree_bv->Branch("dip_15402", &output_bv.dip_15402);
        tree_bv->Branch("dip_filter", &output_bv.dip_filter);

        //Source tree
        tree_ann = new TTree("GIFpp_DATA_ANNEALING", "GIFpp data from annealing");
       
        tree_ann->Branch("run", &output_ann.run);
        tree_ann->Branch("time_stamp", &output.timestamp);
        tree_ann->Branch("MIP1", &output_ann.MIP1);
        tree_ann->Branch("eMIP1", &output_ann.eMIP1); 
        tree_ann->Branch("MIP2", &output_ann.MIP2);
        tree_ann->Branch("eMIP2", &output_ann.eMIP2);
        tree_ann->Branch("Gain1", &output.Gain2);
        tree_ann->Branch("eGain1", &output.eGain2);  
        tree_ann->Branch("Gain2", &output.Gain2);
        tree_ann->Branch("eGain2", &output.eGain2); 
        tree_ann->Branch("BV1", &output.BV1);
        tree_ann->Branch("T1", &output.T1);
        tree_ann->Branch("BV2", &output.BV2);
        tree_ann->Branch("T2", &output.T2);

    }
    switch (parser.ParseThis(parameters.runkey))
    {
        case 1:
            ProcessingBVSCAN();
            break;
        case 2:
            ProcessingLED();
            break;
        case 3:
            ProcessingIrradiation();
            break;
        case 4:
            ProcessingAnnealing();
            break;
            
    }
    //if (parameters.runkey == "ALL") ProcessingAll();
}

void Analyzer::ProcessingLED()
{
    std::string line;
    std::size_t found;
    int runnum;
    //std::cout << "HERE"<< std::endl;
    while (std::getline(runlist_file, line))
    {
        found = line.find("run");
        if (found != std::string::npos)
        {   
            runnum = std::stoi(line.substr(4, 4));
            if ((runnum >= parameters.runstart) && (runnum <= parameters.runstop)) //std::cout << std::stoi(line.substr(4, 4)) << std::endl;
            {
                found = line.find("type");
                runtype_ = line.substr(found + 5, line.size() - (found + 4));
                if (runtype_ == "LED")
                {
                    if (parameters.verbose > 0) std::cout << "Processing run " << runnum << " runtype " << runtype_ << std::endl;
                	AnalyzeRun(runnum);
                }
            }
        }
    }
}

void Analyzer::ProcessingBVSCAN()
{
    std::string line;
    std::size_t found;
    int runnum, prevrun = 0;
    //std::cout << "HERE"<< std::endl;
    while (std::getline(runlist_file, line))
    {
        found = line.find("run");
        if (found != std::string::npos)
        {   
            runnum = std::stoi(line.substr(4, 4));
            if ((runnum >= parameters.runstart) && (runnum <= parameters.runstop) && (runnum > prevrun)) //std::cout << std::stoi(line.substr(4, 4)) << std::endl;
            {
                found = line.find("type");
                runtype_ = line.substr(found + 5, line.size() - (found + 4));
                if (runtype_ == "BVSCAN")
                {
                    if (parameters.verbose > 0) std::cout << "Processing runs " << runnum << " " << runnum + 8 << " runtype " << runtype_ << std::endl;
                    AnalyzeBVSCAN(runnum, runnum + 8);
                    prevrun = runnum + 8;
                }
            }
        }
    }
}

void Analyzer::ProcessingIrradiation()
{
    std::string line;
    std::size_t found;
    int runnum, prevrun = 0;

    while (std::getline(runlist_file, line))
    {
        found = line.find("run");
        if (found != std::string::npos)
        {   
            runnum = std::stoi(line.substr(4, 4));
            if ((runnum >= parameters.runstart) && (runnum <= parameters.runstop) && (runnum > prevrun)) //std::cout << std::stoi(line.substr(4, 4)) << std::endl;
            {
                found = line.find("type");
				runtype_ = line.substr(found + 5, line.size() - (found + 4));
				if (runtype_ == "LED") {
                // std::cout << line.substr(found + 5, line.size() - (found + 4)) << std::endl;
                	if (parameters.verbose > 0) std::cout << "Processing run " << runnum << " runtype " << runtype_ << std::endl;
                	AnalyzeRun(runnum);
                    prevrun = runnum;
				}
                if (runtype_ == "BVSCAN") {
                    if (parameters.verbose > 0) std::cout << "Processing runs " << runnum << " " << runnum + 8 << " runtype " << runtype_ << std::endl;
                    AnalyzeBVSCAN(runnum, runnum + 8);
                    prevrun = runnum + 8;
                }

            }
        }
        else std::cout << "Just skipping this ->" << line << std::endl;
    }
}

void Analyzer::ProcessingAnnealing()
{
    std::string line;
    std::size_t found;
    int runnum, prevrun = 0;

    while (std::getline(runlist_file, line))
    {
        found = line.find("run");
        if (found != std::string::npos)
        {   
            runnum = std::stoi(line.substr(4, 4));
            if ((runnum >= parameters.runstart) && (runnum <= parameters.runstop) && (runnum > prevrun)) //std::cout << std::stoi(line.substr(4, 4)) << std::endl;
            {
                found = line.find("type");
				runtype_ = line.substr(found + 5, line.size() - (found + 4));
				if (runtype_ == "LED") {
                // std::cout << line.substr(found + 5, line.size() - (found + 4)) << std::endl;
                	if (parameters.verbose > 0) std::cout << "Processing run " << runnum << " runtype " << runtype_ << std::endl;
                	AnalyzeRun(runnum);
                    prevrun = runnum;
				}
                if (runtype_ == "SOURCE") {
                // std::cout << line.substr(found + 5, line.size() - (found + 4)) << std::endl;
                	if (parameters.verbose > 0) std::cout << "Processing run " << runnum << " runtype " << runtype_ << std::endl;
                	AnalyzeSource(runnum);
                    prevrun = runnum;
				}
                if (runtype_ == "BVSCAN") {
                    if (parameters.verbose > 0) std::cout << "Processing runs " << runnum << " " << runnum + 8 << " runtype " << runtype_ << std::endl;
                    AnalyzeBVSCAN(runnum, runnum + 8);
                    prevrun = runnum + 8;
                }
                if (runtype_ == "COSMIC") {
                // std::cout << line.substr(found + 5, line.size() - (found + 4)) << std::endl;
                	if (parameters.verbose > 0) std::cout << "Processing run " << runnum << " runtype " << runtype_ << std::endl;
                	AnalyzeCOSMIC(runnum);
                    prevrun = runnum;
				}
            }
        }
        else std::cout << "Just skipping this ->" << line << std::endl;
    }
}

void Analyzer::AnalyzeSource(int run) {

    output.run = run;
    output.runtype = runtype_;
	std::string filename;
	filename +=	pathname;
	if (run<1000) filename += "output00000";
	else filename += "output0000";
	filename +=  std::to_string(run);
	filename += ".root";

	TFile *file = TFile::Open(filename.c_str());
    if (!file)
	{
		if (parameters.verbose > 0) std::cout << "Run " << run << " file doesn't exist." << std::endl;
		return;
	}

	TH1F *htemp;
	file->GetObject("ch2/charge_spectrum_ch2", htemp);

	float peak;
	htemp->Rebin(100);
	TSpectrum s(2);
	s.Search(htemp, 3, "goff", 0.001);
	double *_xpeaks = s.GetPositionX();
	if (_xpeaks[0] > _xpeaks[1]) peak = _xpeaks[0];
	else peak = _xpeaks[1];

	TF1 *g = new TF1("MIP", "gaus", peak-peak/3, peak+peak/5);
	g->SetLineColor(kCyan);
	g->SetParameter(1, peak);

	htemp->Fit(g,"RQ0+");
	htemp->GetFunction("MIP")->ResetBit(TF1::kNotDraw);
    output_ann.MIP2 = g->GetParameter(1);
    output_ann.eMIP2 = g->GetParError(1);
    output_ann.run = run;
    
    if (!parameters.stdout_flag) {
        std::string str_result = "MIP_results.txt";
        FILE* ptr_res = fopen(str_result.c_str(), "a+");
        if (ptr_res == NULL)
        {
            printf("cannot open output file\n");
                exit(0);
        }
        fprintf(ptr_res,"%d\t%ld\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", output_ann.run, output.timestamp,
                                                                        output_ann.MIP2, output_ann.eMIP2,
                                                                        output.Gain2, output.eGain2,
                                                                        output.mean2, output.baseline2,
                                                                        output.rms2, output.rmsb2,
                                                                        output.current2, output.BV2,
                                                                        output.T2);
        fclose(ptr_res);
    }
    else {
        printf("%d\t%ld\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", output_ann.run, output.timestamp,
                                                                        output_ann.MIP2, output_ann.eMIP2,
                                                                        output.Gain2, output.eGain2,
                                                                        output.mean2, output.baseline2,
                                                                        output.rms2, output.rmsb2,
                                                                        output.current2, output.BV2,
                                                                        output.T2);
    }

    if (parameters.root_tree)
    {        
        tree_ann->Fill();
    }
	
    delete htemp;
    delete g;
    file->Close();
    delete file;

}

void Analyzer::AnalyzeCOSMIC(int run) {

    output_ann.run = run;
    output.runtype = runtype_;
	std::string filename;
	filename +=	pathname;
	if (run<1000) filename += "output00000";
	else filename += "output0000";
	filename +=  std::to_string(run);
	filename += ".root";

	TFile *file = TFile::Open(filename.c_str());
    if (!file)
	{
		if (parameters.verbose > 0) std::cout << "Run " << run << " file doesn't exist." << std::endl;
		return;
	}

	TH1F *htemp1, *htemp2;
	file->GetObject("ch1/charge_spectrum_ch1", htemp1);
    file->GetObject("ch2/charge_spectrum_ch2", htemp2);

	float peak;
	htemp1->Rebin(500);
	TSpectrum s(2);
	s.Search(htemp1, 3, "goff", 0.001);
	double *_xpeaks = s.GetPositionX();
	if (_xpeaks[0] > _xpeaks[1]) peak = _xpeaks[0];
	else peak = _xpeaks[1];

	TF1 *g1 = new TF1("MIP1", "gaus", peak-peak/3, peak+peak/5);
	g1->SetLineColor(kCyan);
	g1->SetParameter(1, peak);

	htemp1->Fit(g1,"RQ0+");
	htemp1->GetFunction("MIP1")->ResetBit(TF1::kNotDraw);
    output_ann.MIP1 = g1->GetParameter(1);
    output_ann.eMIP1 = g1->GetParError(1);

    htemp2->Rebin(500);
	s.Search(htemp2, 3, "goff", 0.001);
	_xpeaks = s.GetPositionX();
	if (_xpeaks[0] > _xpeaks[1]) peak = _xpeaks[0];
	else peak = _xpeaks[1];

	TF1 *g2 = new TF1("MIP2", "gaus", peak-peak/3, peak+peak/5);
	g2->SetLineColor(kCyan);
	g2->SetParameter(1, peak);

	htemp2->Fit(g2,"RQ0+");
	htemp2->GetFunction("MIP2")->ResetBit(TF1::kNotDraw);
    output_ann.MIP2 = g2->GetParameter(1);
    output_ann.eMIP2 = g2->GetParError(1);

    
    if (!parameters.stdout_flag) {
        std::string str_result = "MIP_results.txt";
        FILE* ptr_res = fopen(str_result.c_str(), "a+");
        if (ptr_res == NULL)
        {
            printf("cannot open output file\n");
                exit(0);
        }
        fprintf(ptr_res,"%d\t%ld\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", output_ann.run, output.timestamp,
                                                                                    output_ann.MIP1, output_ann.eMIP1,
                                                                                    output_ann.MIP2, output_ann.eMIP2,
                                                                                    output.Gain1, output.eGain1,
                                                                                    output.Gain2, output.eGain2,
                                                                                    output.BV1, output.T1,
                                                                                    output.BV2, output.T2);
        fclose(ptr_res);
    }
    else {
        printf("%d\t%ld\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", output_ann.run, output.timestamp,
                                                                            output_ann.MIP1, output_ann.eMIP1,
                                                                            output_ann.MIP2, output_ann.eMIP2,
                                                                            output.Gain1, output.eGain1,
                                                                            output.Gain2, output.eGain2,
                                                                            output.BV1, output.T1,
                                                                            output.BV2, output.T2);
    }

    if (parameters.root_tree)
    {        
        tree_ann->Fill();
    }
	
    delete htemp1;
    delete g1;
    delete htemp2;
    delete g2;
    file->Close();
    delete file;

}

void Analyzer::AnalyzeBVSCAN(int runstart, int runstop)
{
    	int range = runstop - runstart;
    	int k = 0;
    	float BVs1[range+1], eBVs1[range+1], Gains1[range+1], eGains1[range+1], BVs2[range+1], eBVs2[range+1], Gains2[range+1], eGains2[range+1], dT1[range+1], dT2[range+1], time_stamps[range+1];
    	float T1 = 0, T2 = 0, dip_15402 = 0, dip_filter = 0;

    	for (int i = runstart; i <= runstop; i++)
    	{	
            AnalyzeRun(i);

    		BVs1[k] = output.BV1;
    		eBVs1[k] = 0.01;
    		Gains1[k] = output.Gain1;
    		eGains1[k] = output.eGain1;
    		T1 += output.T1;
            dT1[k] = output.T1;

    		BVs2[k] = output.BV2;
    		eBVs2[k] = 0.01;
    		Gains2[k] = output.Gain2;
    		eGains2[k] = output.eGain2;
    		T2 += output.T2;
            dT2[k] = output.T2;

            time_stamps[k] = output.timestamp;
            dip_15402 += output.dip_15402;
            dip_filter += output.dip_filter;

    		k++;           
    	}

    	TGraphErrors *gr1 = new TGraphErrors(k, BVs1, Gains1, eBVs1, eGains1);
    	TF1 *f1 = new TF1("linear_ch1", "[1]*(x-[0])", 20, 50);
        f1->SetParameter(0, 37.8);
        f1->SetParameter(1, 110);
    	gr1->Fit(f1,"RQM0+");
        gr1->GetFunction("linear_ch1")->ResetBit(TF1::kNotDraw);

    	TGraphErrors *gr2 = new TGraphErrors(k, BVs2, Gains2, eBVs2, eGains2);
    	TF1 *f2 = new TF1("linear_ch2", "[1]*(x-[0])", 20, 50);
        f2->SetParameter(0, 37.9);
        f2->SetParameter(1, 105);
    	gr2->Fit(f2,"RQM0+");
        gr2->GetFunction("linear_ch2")->ResetBit(TF1::kNotDraw);

        output_bv.runstart = runstart;
        output_bv.runstop = runstop;
        output_bv.timestamp = time_stamps[4];

        output_bv.BDV1 = f1->GetParameter(0);
        output_bv.BDV2 = f2->GetParameter(0);
        output_bv.eBDV1 = f1->GetParError(0);
        output_bv.eBDV2 = f2->GetParError(0);
        output_bv.Slope1 = f1->GetParameter(1);
        output_bv.Slope2 = f2->GetParameter(1);
        output_bv.eSlope1 = f1->GetParError(1);
        output_bv.eSlope2 = f2->GetParError(1);
        output_bv.T1 = 1.0*T1/k;
        output_bv.T2 = 1.0*T2/k;
        output_bv.eT1 = dT1[0]-dT1[range];
        output_bv.eT2 = dT2[0]-dT2[range];

        output_bv.dip_15402 = 1.0*dip_15402/k;
        output_bv.dip_filter = 1.0*dip_filter/k;

        if (!parameters.stdout_flag) {
            std::string str_result = "BVSCAN_results.txt";
            FILE* ptr_res = fopen(str_result.c_str(), "a+");
            if (ptr_res == NULL)
            {
                printf("cannot open output file\n");
                exit(0);
            }
            fprintf(ptr_res,"%d\t%d\t%ld\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", output_bv.runstart, output_bv.runstop, output_bv.timestamp,
                                                                                        output_bv.BDV1, output_bv.BDV2,
                                                                                        output_bv.eBDV1, output_bv.eBDV2,
                                                                                        output_bv.Slope1, output_bv.Slope2,
                                                                                        output_bv.eSlope1, output_bv.eSlope2,
                                                                                        output_bv.T1, output_bv.T2,
                                                                                        output_bv.eT1, output_bv.eT2,
                                                                                        output_bv.dip_15402, output_bv.dip_filter);
            fclose(ptr_res);
        }
        else {
            printf("%d\t%d\t%ld\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", output_bv.runstart, output_bv.runstop, output_bv.timestamp,
                                                                                        output_bv.BDV1, output_bv.BDV2,
                                                                                        output_bv.eBDV1, output_bv.eBDV2,
                                                                                        output_bv.Slope1, output_bv.Slope2,
                                                                                        output_bv.eSlope1, output_bv.eSlope2,
                                                                                        output_bv.T1, output_bv.T2,
                                                                                        output_bv.eT1, output_bv.eT2,
                                                                                        output_bv.dip_15402, output_bv.dip_filter);
        }

        if (parameters.root_tree)
        {        
            tree_bv->Fill();
        }

    	delete gr1;
    	delete gr2;
    	delete f1;
    	delete f2;
}

void Analyzer::AnalyzeRun(int run)
{
    output.run = run;
    output.runtype = runtype_;
	std::string filename;
	filename +=	pathname;
	if (run<1000) filename += "output00000";
	else filename += "output0000";
	filename +=  std::to_string(run);
	filename += ".root";

	TFile *file = TFile::Open(filename.c_str());	
	if (!file)
	{
		if (parameters.verbose > 0) std::cout << "Run " << run << " file doesn't exist." << std::endl;
		return;
	}
    else fitRun(file, run);

	if (!parameters.stdout_flag) {
		std::string str_result = "irradiation_results.txt";
		FILE* ptr_res = fopen(str_result.c_str(), "a+");
		if (ptr_res == NULL)
		{
		    printf("cannot open output file\n");
		    exit(0);
		}
		// std::cout << TimeConverterToSec(parameters.data_start, o.timestart) << std::endl;
		fprintf(ptr_res,"%d\t%s\t%ld\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", output.run,
												                                    output.runtype.c_str(), output.timestamp,
				                                                                    output.Gain1, output.Gain2,
				                                                                    output.eGain1, output.eGain2,
				                                                                    output.mean1, output.mean2,
				                                                                    output.baseline1, output.baseline2,
				                                                                    output.rms1, output.rms2,
				                                                                    output.rmsb1, output.rmsb2,
				                                                                    output.current1, output.current2,
                                                                                    output.BV1, output.BV2,
				                                                                    output.T1, output.T2,
				                                                                    output.LYSO_Yield, output.LYSO_Gain,
												                                    output.dip_15402, output.dip_filter);

    	fclose(ptr_res);

	}
	else {
		printf("%d\t%s\t%ld\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", output.run,
												                                    output.runtype.c_str(), output.timestamp,
				                                                                    output.Gain1, output.Gain2,
				                                                                    output.eGain1, output.eGain2,
				                                                                    output.mean1, output.mean2,
				                                                                    output.baseline1, output.baseline2,
				                                                                    output.rms1, output.rms2,
				                                                                    output.rmsb1, output.rmsb2,
				                                                                    output.current1, output.current2,
                                                                                    output.BV1, output.BV2,
				                                                                    output.T1, output.T2,
				                                                                    output.LYSO_Yield, output.LYSO_Gain,
												                                    output.dip_15402, output.dip_filter);
	}
    if (parameters.root_tree)
    {        
        tree->Fill();
        
        std::string path = "run_" + std::to_string(run);
        root_file->mkdir(path.c_str());
        hcoll.h_ch1->SetDirectory(root_file->GetDirectory(path.c_str()));
        hcoll.h_ch1->Write();
        hcoll.h_ch2->SetDirectory(root_file->GetDirectory(path.c_str()));
        hcoll.h_ch2->Write();
        hcoll.h_ch1b->SetDirectory(root_file->GetDirectory(path.c_str()));
        hcoll.h_ch1b->Write();
        hcoll.h_ch2b->SetDirectory(root_file->GetDirectory(path.c_str()));
        hcoll.h_ch2b->Write();
        if (parameters.runkey != "ANNEALING") {
            hcoll.h_ch3->SetDirectory(root_file->GetDirectory(path.c_str()));
            hcoll.h_ch3->Write();
            hcoll.h_ch3b->SetDirectory(root_file->GetDirectory(path.c_str()));
            hcoll.h_ch3b->Write();
        }
        // tree->Write();
        //std::cout << output.runtype.c_str() << " " << output.timestamp << " " << output.Gain1 << std::endl;
    }

    file->Close();
    delete file;    
}

void Analyzer::fitRun(TFile *file, int run)
{
    // std::cout << "Analyzing file " << file->GetName() << std::endl;
    // TFile *file = TFile::Open(filename.c_str());
    TH1F *htemp1, *htemp2, *hBV1, *hBV2, *hI1, *hI2, *hT1, *hT2, *hB1, *hB2, *hL3, *hLb3, *hDip1, *hDip2, *hDipF;// *hQ4
    TH2F *hmean_rms1, *hmean_rms2;
    TProfile *hmean1, *hmean2, *hrms1, *hrms2;
	TTree *tr;

    gErrorIgnoreLevel = 6001;

    file->GetObject("ch1/charge_spectrum_cut_ch1", htemp1);//add cut!!!
    file->GetObject("ch1/waveformMEAN_ch1",hmean1);
    file->GetObject("ch1/waveformRMS_ch1",hrms1);
    file->GetObject("ch1/Mean_vs_RMS_ch1",hmean_rms1);
    file->GetObject("ch1/baseline_spectrum_ch1",hB1);
    file->GetObject("gpib/gpib_0",hI1);
    file->GetObject("gpib/gpib_2",hBV1);
    file->GetObject("temperatures/temperature_0",hT1);
    // std::cout << "HERE" << std::endl;
    TF1 *ch1 = GPFit(htemp1);

    file->GetObject("ch2/charge_spectrum_cut_ch2",htemp2);//add cut!!!
    file->GetObject("ch2/waveformMEAN_ch2",hmean2);
    file->GetObject("ch2/waveformRMS_ch2",hrms2);
    file->GetObject("ch2/Mean_vs_RMS_ch2",hmean_rms2);
    file->GetObject("ch2/baseline_spectrum_ch2",hB2);
    file->GetObject("gpib/gpib_1",hI2);
    file->GetObject("gpib/gpib_3",hBV2);
    file->GetObject("temperatures/temperature_1",hT2);
    TF1 *ch2 = GPFit(htemp2);

    file->GetObject("ch3/charge_spectrum_ch3",hL3);
    if (hL3) {
        TF1 *ch3 = GFit(hL3, 256, 1, 3);
        output.LYSO_Yield = ch3->GetParameter(1);
        delete ch3;
    }
    else {
        output.LYSO_Yield = 0;
    }

    file->GetObject("ch3/charge_spectrum_cut_ch3",hLb3);
	
    if (hLb3) {
        if (hLb3->GetEntries()!=0) {
            TF1 *ch3b = DGFit(hLb3, 4);
            output.LYSO_Gain = ch3b->GetParameter(1) - ch3b->GetParameter(0);
            delete ch3b;
        }
        else {
            output.LYSO_Gain = 1;
        } 
    }

    output.Gain1 = ch1->GetParameter(1);
    output.eGain1 = ch1->GetParError(1);

    output.Gain2 = ch2->GetParameter(1);
    output.eGain2 = ch2->GetParError(1);


    output.BV1 = hBV1->Integral()/hBV1->GetEntries();
    output.BV2 = hBV2->Integral()/hBV2->GetEntries();

    output.current1 = hI1->Integral()/hI1->GetEntries()*1e9;
    output.current2 = hI2->Integral()/hI2->GetEntries()*1e9;

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
    TF1 *ch1b = GFit(hB1, 1, 3, 1);
    //output.baseline1 = hB1->GetXaxis()->GetBinCenter(hB1->GetMaximumBin());
    output.baseline1 = ch1b->GetParameter(1);
    output.rmsb1 = hB1->GetRMS();

    //hmean_rms2->GetXaxis()->SetRangeUser(-100,10);
    //output.baseline2 = hmean_rms2->GetMean(2);
    //output.rmsb2 = hmean_rms2->GetRMS(2);
    //output.baseline2 = hB2->GetMean();
    TF1 *ch2b = GFit(hB2, 1, 3, 1);
    //output.baseline2 = hB2->GetXaxis()->GetBinCenter(hB2->GetMaximumBin());
    output.baseline2 = ch2b->GetParameter(1);
    output.rmsb2 = hB2->GetRMS();

	file->GetObject("DIP/DIP_0",hDip1);
	file->GetObject("DIP/DIP_1",hDip2);
	file->GetObject("DIP/DIP_2",hDipF);

	output.dip_15402 = hDip1->Integral()/hDip1->GetEntries();
	output.dip_15403 = hDip2->Integral()/hDip2->GetEntries();

    float bsum = 0;
    int nsum = 0, nlast = hDipF->FindLastBinAbove(0);
    for (int i = 0; i < nlast; ++i) {
        if (hDipF->GetBinContent(i) > 0) {
            bsum += hDipF->GetBinContent(i);
            nsum++;
        }
    }
	output.dip_filter = 1.0*bsum/nsum;

	file->GetObject("Run_Data", tr);
	TString *cstr;
	tr->SetBranchAddress("time_start", &cstr); 
	tr->GetEntry(0);
	std::string sample_string(cstr->Data());
	sample_string = sample_string.substr(0,24);

    output.timestamp = TimeConverterToSec(parameters.data_start, sample_string);

    std::string ch1_name = "run" + std::to_string(run) + "_ch1_charge_cut";
    std::string ch2_name = "run" + std::to_string(run) + "_ch2_charge_cut";
    std::string ch3_name = "run" + std::to_string(run) + "_ch3_charge_cut";
    std::string ch3_name2 = "run" + std::to_string(run) + "_ch3_charge";
    std::string ch1b_name = "run" + std::to_string(run) + "_ch1_baseline";
    std::string ch2b_name = "run" + std::to_string(run) + "_ch2_baseline";

    hcoll.h_ch1 = (TH1F*)htemp1->Clone(ch1_name.c_str());
    hcoll.h_ch1->SetDirectory(0);
    hcoll.h_ch2 = (TH1F*)htemp2->Clone(ch2_name.c_str());
    hcoll.h_ch2->SetDirectory(0);
    if (hLb3) {
        hcoll.h_ch3b = (TH1F*)hLb3->Clone(ch3_name.c_str());
        hcoll.h_ch3b->SetDirectory(0);
    }
    if (hLb3) {
        hcoll.h_ch3 = (TH1F*)hL3->Clone(ch3_name2.c_str());
        hcoll.h_ch3->SetDirectory(0);
    }
    hcoll.h_ch1b = (TH1F*)hB1->Clone(ch1b_name.c_str());
    hcoll.h_ch1b->SetDirectory(0);
    hcoll.h_ch2b = (TH1F*)hB2->Clone(ch2b_name.c_str());
    hcoll.h_ch2b->SetDirectory(0);


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
    delete hL3;
    delete hLb3;
	delete ch1b;
	delete ch2b;
	delete hB2;
	delete hB1;
	delete hDip1;
	delete hDip2;
	delete hDipF;
    delete tr;
	//delete cstr;

}


TF1 *Analyzer::GPFit(TH1F* hist2proc, int xlow, int xhigh){

    // TCanvas c1("c1", "c1")
    double npe_peaks, mu_hist, par[6];
    double binwidth;
    int ntotal, n0;

    hist2proc->Rebin(10);
    //hist2proc->GetXaxis()->UnZoom();
    binwidth=hist2proc->GetBinWidth(1);
    // uncomment this for day 7
    xlow = hist2proc->GetBinCenter(hist2proc->FindFirstBinAbove(0,1));
    xhigh = hist2proc->GetBinCenter(hist2proc->FindLastBinAbove(0,1)); // in case of only GP in histogram
    hist2proc -> GetXaxis() -> SetRangeUser(xlow, xhigh);
    ntotal=hist2proc->Integral();
    if (parameters.verbose>1) std::cout<<"Spectrum created"<<std::endl;

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

    if (parameters.verbose>1) std::cout << "Gauss fiting 1" << std::endl;
    // Fit first two peaks with Gauss function. Fit range is 1/3 of distance between two peaks
    TF1 g1("peak0","gaus",xpeaks[0]-(xpeaks[1]-xpeaks[0])/3,xpeaks[0]+(xpeaks[1]-xpeaks[0])/3);
    g1.SetLineColor(kCyan);
    hist2proc -> Fit(&g1,"RQ0");
    hist2proc->GetFunction("peak0")->ResetBit(TF1::kNotDraw);

    if (parameters.verbose>1) std::cout << "Gauss fiting 2" << std::endl;
    TF1 g2("peak1","gaus",xpeaks[1]-(xpeaks[1]-xpeaks[0])/3,xpeaks[1]+(xpeaks[1]-xpeaks[0])/3);
    g2.SetLineColor(kCyan);
    hist2proc -> Fit(&g2,"RQ0");
    hist2proc->GetFunction("peak1")->ResetBit(TF1::kNotDraw);
    fp.g1Chi2 = g2.GetChisquare();
    fp.g0NDF = g1.GetNDF();
    fp.g1NDF = g2.GetNDF();

    fp.g0Integral = g1.Integral(xpeaks[0]-(xpeaks[1]-xpeaks[0])/3,xpeaks[0]+(xpeaks[1]-xpeaks[0])/3);
    fp.g1Integral = g2.Integral(xpeaks[1]-(xpeaks[1]-xpeaks[0])/3,xpeaks[1]+(xpeaks[1]-xpeaks[0])/3);

    g1.SetLineColor(kGreen+2);
    g2.SetLineColor(kGreen+2);
    g1.GetParameters(&par[0]);
    g2.GetParameters(&par[3]);

    // Get mean and sigma of two Gaussians
    fp.gN0 = g1.GetParameter(0); // g1params[0]
    fp.gm0 = g1.GetParameter(1);
    fp.gs0 = g1.GetParameter(2);
    fp.gN1 = g2.GetParameter(0);
    fp.gm1 = g2.GetParameter(1);
    fp.gs1 = g2.GetParameter(2);
    fp.gain_spe_01 = fp.gm1 - fp.gm0;


    // Set new range for histogram (Mu1-Mu0)/2 around Mu0
    hist2proc->GetXaxis()->SetRangeUser(fp.gm0-(fp.gm1-fp.gm0)/2,fp.gm0+(fp.gm1-fp.gm0)/2);
    n0 = hist2proc -> Integral();
    hist2proc->GetXaxis()->SetRangeUser(xlow, xhigh);
    mu_hist= -TMath::Log(1.*n0/ntotal);
    //npe_peaks = 3;
    npe_peaks = TMath::Nint((xhigh-fp.gm0)/(fp.gm1-fp.gm0));
    //cout << npe_peaks << endl;
    

    // Construct GP function in range [Mu0-(Mu1-Mu0)/2, xhigh]
    TF1 *func = new TF1("generpoiss", generpoiss, fp.gm0-(fp.gm1-fp.gm0)/2, xhigh, 8);
    // Set initial parameters for GP and limits for them
    func -> SetParameters( 0.1, fp.gm1-fp.gm0, sqrt(abs(fp.gs1*fp.gs1 - fp.gs0*fp.gs0)), mu_hist, fp.gm0, fp.gs0, binwidth*ntotal, npe_peaks); //initial parameters
    func -> SetParLimits(1, 0, 1.5*(fp.gm1-fp.gm0));
    //func -> FixParameter(0, 0);
    //func -> FixParameter(7, npe_peaks);
    func -> FixParameter(7, 5);
    func -> SetParNames("lambda", "gain_SPE", "sigma_gain", "mu_avg", "mu0", "sigma0", "NxW", "npe_peaks");
    func -> SetNpx(1000);
    if (parameters.verbose>1) std::cout<<"GP fiting"<<std::endl;
    hist2proc -> Fit(func, "RQ0+");
    hist2proc->GetFunction("generpoiss")->ResetBit(TF1::kNotDraw);
    //hist2proc->GetXaxis()->UnZoom();

    // put parameters of GP to fp structure
    fp.gain = func->GetParameter(1);
    fp.u0 = func->GetParameter(4);
    fp.s0 = func->GetParameter(5);
    fp.sG = func->GetParameter(2);
    fp.lambda = func->GetParameter(0);
    fp.navg = func->GetParameter(3);
    fp.npeaks = func->GetParameter(7);

    fp.gain_err = func->GetParError(1);

    double Chi2 = func->GetChisquare();
    double NDF = func->GetNDF();

    fp.chi2_ndf = Chi2/NDF;

    if (parameters.verbose>1) {

        std::cout << std::endl;
        std::cout << "GenPoiss Fit is Done. Parameters:" << std::endl;
        std::cout << "Baseline Sigma0 = " << fp.s0 << std::endl;
        std::cout << "SPE Gain = " << fp.gain << " +- " << fp.gain_err << std::endl;
        std::cout << "SPE Gain Width = " << fp.sG << std::endl;
        std::cout << "SPE Lambda = " << fp.lambda << std::endl;
        std::cout << "SPE Mu = " << fp.navg << std::endl;
    
        std::cout << "Chi^2/NDF is " << Chi2/NDF << std::endl;
        std::cout << std::endl;
    }
	
	//delete func;

    // return GP function
    return (func);
}

TF1 *Analyzer::GFit(TH1F *hist, int rebin, int multA, int multB) {

    hist->Rebin(rebin);

    double lower_bound = hist->GetXaxis()->GetBinLowEdge(hist->GetMaximumBin()-5);
    double upper_bound = hist->GetXaxis()->GetBinUpEdge(hist->GetMaximumBin()+5);
    int gnpar = 3;
    double gparameters[gnpar];

    if (parameters.verbose>1) {
        std::cout << "1 iteration Gaussian Fit...." << std::endl;
        std::cout << std::endl;
    }

    TF1 *G_func = new TF1("1 fit","gaus",lower_bound, upper_bound);
    hist->Fit(G_func,"RQ0+");
    hist->GetFunction("1 fit")->ResetBit(TF1::kNotDraw);
    G_func->GetParameters(gparameters);

    double Chi2 = G_func->GetChisquare();
    double NDF = G_func->GetNDF();

    lower_bound = gparameters[1]-gparameters[2]/3*multA;
    upper_bound = gparameters[1]+gparameters[2]/3*multB;

    if (parameters.verbose > 1) {
        std::cout << "2 iteration Gaussian Fit...." << std::endl;
        std::cout << std::endl;
    }

    TF1 *LG_func = new TF1("2 fit","gaus",lower_bound, upper_bound);
    LG_func->SetParameters(gparameters);
    LG_func->SetLineColor(kGreen);
    hist->Fit(LG_func,"RQ0+");
    hist->GetFunction("2 fit")->ResetBit(TF1::kNotDraw);
    LG_func->GetParameters(gparameters);

    if (parameters.verbose>1) {
        std::cout << "Gaussian Fit is Done. Parameters:" << std::endl;
        std::cout << "Maximum Gauss = " << gparameters[0] << std::endl;
        std::cout << "Mean Gauss = " << gparameters[1] << std::endl;
        std::cout << "Sigma Gauss = " << gparameters[2] << std::endl;
    
        //Chi2 = LG_func->GetChisquare();
        //NDF = LG_func->GetNDF();

        //std::cout << std::endl;
        //std::cout << "Chi^2/NDF is " << Chi2 << " / " << NDF << std::endl;
        std::cout << "Chi^2/NDF is " << Chi2/NDF << std::endl;
        std::cout << std::endl;

    }

	delete G_func;

    return LG_func;
}

TF1 *Analyzer::DGFit(TH1F *hist, int rebin) {

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

    if (parameters.verbose>1) std::cout << "Double Gauss fiting, peak 1" << std::endl;

    TF1 *g1 = new TF1("peak11","gaus",xpeaks[0]-(xpeaks[1]-xpeaks[0])/3,xpeaks[0]+(xpeaks[1]-xpeaks[0])/3);
    g1->SetLineColor(kCyan);
    g1->SetParameter(1,xpeaks[0]);
    hist->Fit(g1,"RQ0");
    g1->GetParameters(gparameters1);
    hist->GetFunction("peak11")->ResetBit(TF1::kNotDraw);

	if (parameters.verbose>1) std::cout << "Double Gauss fiting, peak 2" << std::endl;

    TF1 *g2 = new TF1("peak12","gaus",xpeaks[1]-(xpeaks[1]-xpeaks[0])/3,xpeaks[1]+(xpeaks[1]-xpeaks[0])/3);
    g2->SetParameter(1,xpeaks[1]);
    g2->SetLineColor(kCyan);
    hist->Fit(g2,"RQ0+");
    g2->GetParameters(gparameters2);
    hist->GetFunction("peak12")->ResetBit(TF1::kNotDraw);


    if (parameters.verbose>1) {

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

	//delete g1;

    return g2;

}


long int Analyzer::TimeConverterToSec(std::string date1, std::string date2)
{
	int year1, mon1, day1, hour1, min1, sec1;
	int year2, mon2, day2, hour2, min2, sec2;
	long int t1, t2;

	year1 = (std::stoi(date1.substr(20,4))-2023)*365*24*3600;
	mon1 = dt.month[date1.substr(4,3)]*24*3600;
	day1 = std::stoi(date1.substr(8,2))*24*3600;
	hour1 = std::stoi(date1.substr(11,2))*3600;
	min1 = std::stoi(date1.substr(14,2))*60;
	sec1 = std::stoi(date1.substr(17,2));
	t1 = sec1 + min1 + hour1 + day1 + mon1 + year1;

	year2 = (std::stoi(date2.substr(20,4))-2023)*365*24*3600;
	mon2 = dt.month[date2.substr(4,3)]*24*3600;
	day2 = std::stoi(date2.substr(8,2))*24*3600;
	hour2 = std::stoi(date2.substr(11,2))*3600;
	min2 = std::stoi(date2.substr(14,2))*60;
	sec2 = std::stoi(date2.substr(17,2));
	t2 = sec2 + min2 + hour2 + day2 + mon2 + year2;
	//std::cout << t1.tm_year << t1.tm_mon << t1.tm_mday << t1.tm_hour << t1.tm_min << t1.tm_sec << std::endl;

	return t2-t1;//static_cast<long int>(timeSinceEpoch);

}