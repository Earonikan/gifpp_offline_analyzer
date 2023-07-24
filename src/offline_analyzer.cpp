#include <iostream>
#include <TSystemDirectory.h>
#include <TROOT.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH2F.h>
#include <TFile.h>
#include <TProfile.h>

#include "Fitter.h"

struct data
{
	float BV1;
	float Gain1;
	float eGain1;
	float mean1;
    float baseline1;
	float rms1;
    float rmsb1;
	float current1;
	float T1;

	float BV2;
	float Gain2;
	float eGain2;
	float mean2;
    float baseline2;
	float rms2;
    float rmsb2;
	float current2;
	float T2;

	float LYSO_Yield;
	float LYSO_eYield;
	float LYSO_Gain;
	float LYSO_eGain;
};

struct gp_data
{
	GenPoisParams ch1;
	GenPoisParams ch2;
};

void fitBVscan(int num1, int num2);
gp_data fitBVstep_spe(const char* fileName);
void fitGain(int run);
void fitSPE(int run);
data fitBVstep(const char* fileName);

int main(int argv, char *argc[])
{
    const char *ext = ".root";
    // const char *ext = "ch1.png";
    TSystemDirectory dir(".", "/home/qfl/online/midas_digi/daq");
    TList *files = dir.GetListOfFiles();
    
    int num1, runstart, runstop;
    runstart = atoi(argc[1]);
    runstop = atoi(argc[2]);
    if (files) {

        TSystemFile *file;
        TString fname;
        TIter next(files);

        while ((file=(TSystemFile*)next())) {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(ext)) {
                // std::cout << fname.Data() << std::endl;
                // sscanf(fname.Data(), "BVscan_run%d_run%d_ch1.png", &num1, &num2);
                // if ((num1>1576))
                // {
                    // std::cout << "Processing BV scan, started from run " << num1+1 << " finished with run " << num2;// << std::endl;
                    // auto start = std::chrono::steady_clock::now();
                    // fitBVscan(num1+1, num2);
                    // auto end = std::chrono::steady_clock::now();
                    // std::cout << ", elapsed time in seconds: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " sec." << std::endl;
                // }
                sscanf(fname.Data(), "output0000%d.root", &num1);
                if ((num1 >= runstart)&&(num1 <= runstop))
                {
                    std::cout << "Processing run " << num1 << std::endl;
                    fitSPE(num1);


                }
            }
        }

        delete file;
    }

    delete files;

    return 0;
}

void fitBVscan(int runstart, int runstop)
{
	int range = runstop - runstart;
	int k = 0;
	float BVs1[range+1], eBVs1[range+1], Gains1[range+1], eGains1[range+1], BVs2[range+1], eBVs2[range+1], Gains2[range+1], eGains2[range+1], dT1[range+1], dT2[range+1];
	float T1 = 0, T2 = 0;
	//float testBV[10] = {40.50, 40.75, 41., 41.25, 41.5, 41.75, 42.0, 42.25, 42.5, 42.75};

	data output;

	TString str;
	str = "BVscan_run";
	str += runstart;
	str += "_run";
	str += runstop;

    std::string str_result = "output_results_BVscans.txt";
    FILE* ptr_res = fopen(str_result.c_str(), "a+");
    if (ptr_res == NULL) {
        printf("cannot open output file\n");
        exit(0);
    }

	//TCanvas *c1 = new TCanvas(str+"_ch1",str+"_ch1",800,600);

	for (int i = runstart; i <= runstop; i++)
	{	
		TString filename;
		if (i<1000) filename += "output00000";
        	else filename += "output0000";
		filename += i;
		filename += ".root";		
		output = fitBVstep(filename);

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

		k++;
	}

	TGraphErrors *gr1 = new TGraphErrors(k,BVs1,Gains1,eBVs1,eGains1);
	TF1 *f1 = new TF1("linear_ch1", "[1]*(x-[0])");
	gr1->Fit(f1,"QM0+");
    gr1->GetFunction("linear_ch1")->ResetBit(TF1::kNotDraw);


	TGraphErrors *gr2 = new TGraphErrors(k,BVs2,Gains2,eBVs2,eGains2);
	TF1 *f2 = new TF1("linear_ch2", "[1]*(x-[0])");
	gr2->Fit(f2,"QM0+");
    gr2->GetFunction("linear_ch2")->ResetBit(TF1::kNotDraw);


	//std::cout << f1->GetParameter(0) << "," << f1->GetParError(0) << "," << f2->GetParameter(0) << "," << f2->GetParError(0) << "," << f1->GetParameter(1) << "," << f1->GetParError(1) << "," << f2->GetParameter(1) << "," << f2->GetParError(1) << std::endl;

	fprintf(ptr_res,"%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",runstart,runstop,f1->GetParameter(0),f1->GetParError(0),f2->GetParameter(0),f2->GetParError(0),f1->GetParameter(1),f1->GetParError(1),f2->GetParameter(1),f2->GetParError(1), 1.0*T1/k, 1.0*T2/k, dT1[0]-dT1[range], dT2[0]-dT2[range]);

	delete gr1;
	delete gr2;
	delete f1;
	delete f2;

    fclose(ptr_res);

}


data fitBVstep(const char* fileName)
{
    TString fileName_full = "/home/qfl/online/midas_digi/daq/";
    fileName_full += fileName;
    //cout << "Analyzing file " << fileName << endl;sssss
    TFile *file = TFile::Open(fileName_full);
    TH1F *htemp1, *htemp2, *hBV1, *hBV2, *hI1, *hI2, *hT1, *hT2, *hB1, *hB2, *hL3, *hLb3;
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

    data output;

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

void fitGain(int run)
{
    std::string str_result = "output_results10.txt";
    FILE* ptr_res = fopen(str_result.c_str(), "a+");
    if (ptr_res == NULL) {
        printf("cannot open output file\n");
        exit(0);
    }

	data o;
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

void fitSPE(int run)
{
    std::string str_result = "output_results11.txt";
    FILE* ptr_res = fopen(str_result.c_str(), "a+");
    if (ptr_res == NULL) {
        printf("cannot open output file\n");
        exit(0);
    }

	gp_data o;
	TString filename;
	if (run<1000) filename += "output00000";
	else filename += "output0000";
	filename += run;
	filename += ".root";		
	o = fitBVstep_spe(filename);

	fprintf(ptr_res,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", run,
                                                                            o.ch1.gain, o.ch2.gain,
                                                                            o.ch1.gain_err, o.ch2.gain_err,
                                                                            o.ch1.u0, o.ch2.u0,
                                                                            o.ch1.s0, o.ch2.s0,
                                                                            o.ch1.sG, o.ch2.sG,
                                                                            o.ch1.lambda, o.ch2.lambda,
                                                                            o.ch1.navg, o.ch2.navg,
                                                                            o.ch1.hmean, o.ch2.hmean,
                                                                            o.ch1.hrms, o.ch2.hrms);

    fclose(ptr_res);
}


gp_data fitBVstep_spe(const char* fileName)
{
    TString fileName_full = "/home/qfl/online/midas_digi/daq/";
    fileName_full += fileName;
    //cout << "Analyzing file " << fileName << endl;sssss
    TFile *file = TFile::Open(fileName_full);
    TH1F *htemp1, *htemp2;

    gErrorIgnoreLevel = 6001;

    gp_data output;
 

    file->GetObject("ch1/charge_spectrum_cut_ch1",htemp1);
    output.ch1 = GPFit_par(htemp1,0);//,-500,2000);

    file->GetObject("ch2/charge_spectrum_cut_ch2",htemp2);
    output.ch2 = GPFit_par(htemp2,0);//,-500,2000);


    delete htemp1;
    delete htemp2;

    file->Close();

    delete file;

    return output;

}

