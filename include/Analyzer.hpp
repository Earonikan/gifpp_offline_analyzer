#pragma once

#include "stdafx.hpp"
#include "DataStructs.hpp"
#include "InputArguments.hpp"

struct FinalData;
struct GenPoisParams;
struct Parameters;

double generpoiss(double *x, double *p);

class Analyzer
{

private:
    // void fitBVscan(int num1, int num2);
    FinalData fitBVstep(std::string pathname);
    void fitGain(int run);

    TF1 *GPFit(TH1F* hist2proc, int xlow = -5000, int xhigh = 50000);
    TF1 *GFit(TH1F *hist, int rebin, int multA, int multB);
    TF1 *DGFit(TH1F *hist, int rebin);

    std::ifstream runlist_file;
    Parameters parameters;
    GenPoisParams fp;
    std::string runtype_;
    std::string pathname = "/home/qfl/online/midas_digi/daq/";
    Date dt;

    void ProcessingMonitoring();
    void ProcessingLED();
    long int TimeConverterToSec(std::string date1, std::string date2);
public:
    Analyzer(int argv, char *argc[]);
    // Analyzer(std::string key, std::string runlist_filename, int runstart, int runstop, bool stdout_flag);
    void Processing();
    ~Analyzer();


};
