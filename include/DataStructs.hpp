#pragma once

#include "stdafx.hpp"

struct Parameters
{
    std::string runlist_filename = "/home/qfl/online/midas_digi/daq/runs_history.txt";
    std::string pathname = "/home/qfl/online/midas_digi/daq/";
    int runstart = 0;
    int runstop = 999999;
    bool stdout_flag = false;
    bool root_tree = false;
    bool annealing = false;
    int verbose = 0;
    std::vector<std::string> argKeys = {"--runtype", "--runlist", "--runnumber", "--runrange", "--stdout", "--verbose", "--root", "--pathname", "--annealing"};
    std::vector<std::string> runTypes = {"BVSCAN", "LED", "MONITORING", "MIP", "ALL"};
    std::string runkey = runTypes[2];
    std::string data_start = "Wed Jun  7 22:39:50 2023";
    std::string data_annealing = "Wed Jun 28 08:10:54 2023";
};

struct Date
{
	std::map<std::string, int> month {{"Jan", 0}, {"Feb", 31}, {"Mar", 59}, {"Apr", 90}, {"May", 120}, {"Jun", 151}, {"Jul", 181}, {"Aug", 212}, {"Sep", 243}, {"Oct", 273}, {"Nov", 304}, {"Dec", 334}};

};

struct FinalData
{
    int run;
    std::string runtype;

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

    float dip_15402;
    float dip_15403;
    float dip_filter;
	
	long int timestamp;
};

struct BVData {
    int runstart;
    int runstop;

    long int timestamp;

    float BDV1;
    float BDV2;
    float eBDV1;
    float eBDV2;
    float Slope1;
    float Slope2;
    float eSlope1;
    float eSlope2;
    float T1;
    float T2;
    float eT1;
    float eT2;

    float dip_15402;
    float dip_filter;

};

struct ANNData {
    int run;

    float MIP1;
    float eMIP1;
    float MIP2;
    float eMIP2;

};

struct HistoCollection
{
    TH1F *h_ch1;
    TH1F *h_ch2;
    TH1F *h_ch3;
    TH1F *h_ch3b;
    TH1F *h_ch1b;
    TH1F *h_ch2b;
};

struct GenPoisParams {

    /* GP fit */
    double gain, u0, s0, sG, lambda, navg, npeaks, chi2_ndf;
    /* GP fit errors*/
    double gain_err;

    /* Gaussian fits */
    double gain_spe_01, gs0, gs1, gm0, gm1;
    double gN0, gN1;
    double g0Chi2, g1Chi2, g0NDF, g1NDF;
    /* Integrals of gaussian fits */
    double g0Integral, g1Integral;

};
