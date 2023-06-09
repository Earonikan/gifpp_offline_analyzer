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
    FinalData fitBVstep(const char* fileName);
    void fitGain(int run);

    TF1 *GPFit(TH1F* hist2proc, bool verbose, int xlow = -5000, int xhigh = 50000);
    TF1 *GFit(TH1F *hist, bool verbose, int rebin, int multA, int multB);
    TF1 *DGFit(TH1F *hist, bool verbose, int rebin);

    std::ifstream runlist_file;
    Parameters parameters;

    // enum RunTypes {};

    std::map<std::string, int> RunTypes = {{"BVSCAN", 1}, {"LED", 2}, {"ALL", 3}, {"COINCIDENCE", 4}, {"COSMIC", 5}, {"LYSO_SELF", 6}};

    GenPoisParams *fp;

    void ProcessingAll();
public:
    Analyzer(int argv, char *argc[]);
    // Analyzer(std::string key, std::string runlist_filename, int runstart, int runstop, bool stdout_flag);
    void Processing();
    ~Analyzer();


};