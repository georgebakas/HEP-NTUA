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

	const int n = 6;
	float eb[n] = {0.4, 0.5, 0.57371, 0.6, 0.633934, 0.7};
	float error_x[n] = {0, 0,9.8e-03, 0,0, 0};

	//Signal Region (2btag)
	float NQCD_2[n]   = {5.2073e+03, 4.6677e+03, 4.4306e+03, 4.3853e+03, 4.3472e+03 , 4.3681e+03};
	float error_y2[n] = {169,450,177,175,178,236};
	
	//1btag region
	float NQCD_1[n] = {4.2122e+04, 4.6380e+04, 4.9321e+04, 5.0315e+04, 5.1537e+04, 5.3287e+04};
	float error_y1[n] = {485,467,536,347,347,430};

	//Control Region (0btag)
	float NQCD_0[n] = {1.5924e+05, 1.6543e+05, 1.6842e+05, 1.6920e+05, 1.7002e+05, 1.7117e+05};
	float error_y0[n] = {562,479,564,440,437,432};
	


	//for each region 1 graph
	TGraphErrors *gr2 = new TGraphErrors(n,eb,NQCD_2,error_x,error_y2);
	gr2->SetTitle("N_{QCD}^{(2)} vs e_{b}");
	gr2->SetMarkerColor(8);
	gr2->SetMarkerStyle(23);
	gr2->GetXaxis()->SetTitle("e_{b}");
	gr2->GetYaxis()->SetTitle("N_{QCD}^{(2)}");
	gr2->GetYaxis()->SetTitleOffset(1.4);
	TCanvas *c2 = new TCanvas("2btag", "2btag", 800, 600);
	gr2->Draw("AP");


	TGraphErrors *gr1 = new TGraphErrors(n, eb, NQCD_1, error_x, error_y1);
	gr1->SetTitle("N_{QCD}^{(1)} vs e_{b}");
	gr1->SetMarkerColor(8);
	gr1->SetMarkerStyle(23);
	gr1->GetXaxis()->SetTitle("e_{b}");
	gr1->GetYaxis()->SetTitle("N_{QCD}^{(1)}");
	gr1->GetYaxis()->SetTitleOffset(1.5);
	TCanvas *c1 = new TCanvas("1btag", "1btag", 800,600);
	gr1->Draw("AP");


	TGraphErrors *gr0 = new TGraphErrors(n, eb, NQCD_0, error_x, error_y0);
	gr0->SetTitle("N_{QCD}^{(0)} vs e_{b}");
	gr0->SetMarkerColor(8);
	gr0->SetMarkerStyle(23);
	gr0->GetXaxis()->SetTitle("e_{b}");
	gr0->GetYaxis()->SetTitle("N_{QCD}^{(0)}");
	gr0->GetYaxis()->SetTitleOffset(1.3);
	TCanvas *c0 = new TCanvas("0btag", "0btag", 800,600);
	gr0->Draw("AP");
}