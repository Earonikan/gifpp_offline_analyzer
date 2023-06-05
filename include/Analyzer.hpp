#pragma once

#include "stdafx.hpp"

struct data;
struct GenPoisParams;

double generpoiss(double *x, double *p);

class Analyzer
{

private:
    // void fitBVscan(int num1, int num2);
    data fitBVstep(const char* fileName);
    void fitGain(int run);

    TF1 *GPFit(TH1F* hist2proc, bool verbose, int xlow=-5000, int xhigh=50000);
    TF1 *GFit(TH1F *hist);

    std::ifstream runlist_file;

    GenPoisParams *fp;
public:
    Analyzer(std::string key, std::string runlist_filename, int runstart, int runstop);
    ~Analyzer() {runlist_file.close();}


};