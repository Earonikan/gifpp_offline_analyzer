#pragma once

#include "stdafx.hpp"
#include "DataStructs.hpp"
#include "InputArguments.hpp"

struct FinalData;
struct BVData;
struct GenPoisParams;
struct Parameters;

double generpoiss(double *x, double *p);

class Analyzer
{

private:
    void fitRun(TFile *file, int run);
    void AnalyzeRun(int run);
    void AnalyzeBVSCAN(int num1, int num2);
    void AnalyzeSource(int run);
    void AnalyzeCOSMIC(int run);

    TF1 *GPFit(TH1F* hist2proc, int xlow = -5000, int xhigh = 50000);
    TF1 *GFit(TH1F *hist, int rebin, int multA, int multB);
    TF1 *DGFit(TH1F *hist, int rebin);

    std::ifstream runlist_file;
    Parameters parameters;
    GenPoisParams fp;
    HistoCollection hcoll;
    std::string runtype_;
    Date dt;
    TFile *root_file;
    TTree *tree;
    TTree *tree_bv;
    TTree *tree_ann;
    FinalData output;
    BVData output_bv;
    ANNData output_ann;

    void ProcessingMonitoring();
    void ProcessingLED();
    void ProcessingBVSCAN();
    void ProcessingAll();
    void ProcessingMIPs();
    long int TimeConverterToSec(std::string date1, std::string date2);
public:
    Analyzer(int argv, char *argc[]);
    // Analyzer(std::string key, std::string runlist_filename, int runstart, int runstop, bool stdout_flag);
    void Processing();
    ~Analyzer();


};
