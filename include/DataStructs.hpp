#pragma once

#include "stdafx.hpp"

struct Parameters
{
    std::string runlist_filename = "/home/qfl/online/midas_digi/daq/runs_history.txt";
    int runstart = 0;
    int runstop = 0;
    bool stdout_flag = 0;
    std::vector<std::string> argKeys = {"--runtype", "--filename", "--runnumber", "--runrange", "--stdout"};
    std::vector<std::string> runTypes = {"BVSCAN", "LED", "ALL", "COINCIDENCE", "COSMIC", "LYSO_SELF"};
    std::string runkey = runTypes[2];
};

struct FinalData
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
