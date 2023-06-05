#pragma once

#include <TH1.h>
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TPaveStats.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TPolyMarker.h"
#include "TSpectrum.h"
#include "TTree.h"
#include "TString.h"
#include "string.h"

#include <fstream>
#include <iostream>
#include <dirent.h>
#include <sstream>


// using namespace std;

// n0 vas not initialized
int n0;

double generpoiss(double *x, double *p){
  /* Function description is here - https://indico.cern.ch/event/725114/contributions/2983038/attachments/1641357/2621413/he_sipm_dpg2704_v2.pdf#page=6 */
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

struct GenPoisParams{
  GenPoisParams():gain(0),u0(0),s0(0),sG(0),lambda(0),navg(0),npeaks(0),chi2_ndf(0),gain_err(0),gain_spe_01(0),gs0(0),gs1(0),gm0(0),gm1(0),gN0(0),gN1(0),g0Chi2(0), g1Chi2(0),g0NDF(0),g1NDF(0),g0Integral(0),g1Integral(0){}

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

  /* Info from histogram */

} fp;

TF1 *GPFit(TH1F* hist2proc, bool verbose, int xlow=-5000, int xhigh=50000){

    // TCanvas c1("c1", "c1")
  double npe_peaks, mu_hist, par[6], lastbin_x;
  double binwidth;
  int ntotal;

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
  hist2proc -> Fit(&g1,"RQM0");
  hist2proc->GetFunction("peak0")->ResetBit(TF1::kNotDraw);

  TF1 g2("peak1","gaus",xpeaks[1]-(xpeaks[1]-xpeaks[0])/3,xpeaks[1]+(xpeaks[1]-xpeaks[0])/3);
  g2.SetLineColor(kCyan);
  hist2proc -> Fit(&g2,"RQM0+");
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
  func -> FixParameter(7, npe_peaks);
  func -> SetParNames("lambda", "gain_SPE", "sigma_gain", "mu_avg", "mu0", "sigma0", "NxW", "npe_peaks");
  func -> SetNpx(1000);
  if (verbose) std::cout<<"GP fiting"<<std::endl;
  hist2proc -> Fit(func, "RQM0+");
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

  if (verbose) {

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

  // return GP function
  return (func);
}

TF1* GFit(TH1F *hist, bool verbose) {

  double max = hist->GetMaximum();
  double lower_bound = hist->GetXaxis()->GetBinLowEdge(hist->FindFirstBinAbove(max/10,1,5,-1));
  double upper_bound = hist->GetXaxis()->GetBinUpEdge(hist->FindLastBinAbove(2*max/3));
  int gnpar = 3, lgnpar = 4;
  double gparameters[gnpar],parameters[lgnpar];

  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function

  /*std::cout << "Starting 1 Gaussian Fit...." << std::endl;
  std::cout << std::endl;*/

  TF1 *G_func = new TF1("1 fit","gaus",lower_bound, upper_bound);
  hist->Fit(G_func,"RQM0+");
  G_func->GetParameters(gparameters);

  /*std::cout << "Gaussian Fit 1 is Done. Parameters:" << std::endl;
  std::cout << "Maximum Gauss = " << gparameters[0] << std::endl;
  std::cout << "Mean Gauss = " << gparameters[1] << std::endl;
  std::cout << "Sigma Gauss = " << gparameters[2] << std::endl;
  std::cout << std::endl;*/

  double Chi2 = G_func->GetChisquare();
  double NDF = G_func->GetNDF();

  /*std::cout << "Chi^2/NDF is " << Chi2 << " / " << NDF << std::endl;
  std::cout << std::endl;*/

  parameters[0] = gparameters[0];
  parameters[1] = gparameters[1];
  parameters[2] = gparameters[2];

  lower_bound = parameters[1]-2*parameters[2];
  upper_bound = parameters[1]+parameters[2];

  //std::cout << "Starting 2 Gaussian Fit...." << std::endl;
  //std::cout << std::endl;  

  TF1 *LG_func = new TF1("2 fit","gaus",lower_bound, upper_bound);
  LG_func->SetParameters(parameters);
  //LG_function->ReleaseParameter(0);
  //LG_func->SetParLimits(0, 10, 5000);
  //LG_func->SetParLimits(3, 0, 500000);
  //LG_func->FixParameter(3, TMath::Sqrt(1049.09 + 27.53166 * gparameters[1]));
  LG_func->SetLineColor(kGreen);
  hist->Fit(LG_func,"RQM0+");
  LG_func->GetParameters(parameters);
  //hist->GetXaxis()->UnZoom();

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

