#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TGraphErrors.h"

using std::cin;
using std::cout;
using std::endl;

void graphNbkg_eb()
{

	const int n = 7;
	float eb[n] = {0.4, 0.48641, 0.5, 0.6, 0.605622, 0.7, 0.8};
	float error_x[n] = {0, 1.10e-02,0, 0,0, 0, 0};

	//Signal Region (2btag)
	float NQCD_2[n]   = {3.0167e+03, 2.503e+03, 2.4147e+03, 2.21e+03, 2.2031e+03 , 1.8895e+03, 0.708e+03};
	float error_y2[n] = {1.30e+02,   1.41e+02, 1.32e+02, 1.34e+02, 1.36e+02 , 0.80e+02, 0.72e+02};
	
	//1btag region
	float NQCD_1[n] = {3.511e+04, 3.6754e+04, 3.7109e+04, 3.9728e+04, 3.9865e+04, 4.177e+04, 4.2947e+04};
	float error_y1[n] = {406, 456, 331, 208, 287, 222, 221};

	//Control Region (0btag)
	float NQCD_0[n] = {1.574e+05, 1.612e+05, 1.617e+05, 1.6429e+05, 1.6439e+05, 1.656e+05, 1.6621e+05};
	float error_y0[n] = {487, 588, 433, 418, 418, 413, 413};
	
	//for each region 1 graph
	TGraphErrors *gr2 = new TGraphErrors(n,eb,NQCD_2,error_x,error_y2);
	gr2->SetTitle("N_{QCD}^{(2)} vs e_{b}");
	gr2->SetMarkerColor(2);
	gr2->SetMarkerStyle(20);
	gr2->GetXaxis()->SetTitle("e_{b}");
	gr2->GetYaxis()->SetTitle("N_{QCD}^{(2)}");
	gr2->GetYaxis()->SetTitleOffset(1.4);
	TCanvas *c2 = new TCanvas("2btag", "2btag", 800, 600);
	gr2->Draw("AP");


	TGraphErrors *gr1 = new TGraphErrors(n, eb, NQCD_1, error_x, error_y1);
	gr1->SetTitle("N_{QCD}^{(1)} vs e_{b}");
	gr1->SetMarkerColor(2);
	gr1->SetMarkerStyle(20);
	gr1->GetXaxis()->SetTitle("e_{b}");
	gr1->GetYaxis()->SetTitle("N_{QCD}^{(1)}");
	gr1->GetYaxis()->SetTitleOffset(1.5);
	TCanvas *c1 = new TCanvas("1btag", "1btag", 800,600);
	gr1->Draw("AP");


	TGraphErrors *gr0 = new TGraphErrors(n, eb, NQCD_0, error_x, error_y0);
	gr0->SetTitle("N_{QCD}^{(0)} vs e_{b}");
	gr0->SetMarkerColor(2);
	gr0->SetMarkerStyle(20);
	gr0->GetXaxis()->SetTitle("e_{b}");
	gr0->GetYaxis()->SetTitle("N_{QCD}^{(0)}");
	gr0->GetYaxis()->SetTitleOffset(1.3);
	TCanvas *c0 = new TCanvas("0btag", "0btag", 800,600);
	gr0->Draw("AP");
}