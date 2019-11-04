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
	float eb[n] = {0.4, 0.5, 0.56029, 0.6, 0.629909, 0.7, 0.8};
	float error_x[n] = {0,0, 1.17e-02, 0,0, 0, 0};

	//Signal Region (2btag)
	float NQCD_2[n] = {4.1645e+03,3.2985e+03, 2.9980e+03,2.8784e+03,2.84e+03 ,2.8091e+03, 2.9322e+03};
	float error_y2[n] = {1.36e+02, 1.37e+02, 1.43e+02, 1.38e+02,1.55e+02 ,1.41e+02, 1.45e+02};
	
	//1btag region
	float NQCD_1[n] = {2.6618e+04, 2.8082e+04, 2.8973e+04, 2.9865e+04,3.0542e+04 ,3.2121e+04, 3.4154e+04};
	float error_y1[n] = {3.80e+02, 3.20e+02, 3.93e+02, 2.76e+02,2.62e+02 ,2.36e+02, 2.13e+02};

	//Control Region (0btag)
	float NQCD_0[n] = {8.1973e+04, 8.5459e+04, 8.7019e+04, 8.7826e+04,8.832e+04 ,8.9197e+04, 9.9892e+04};
	float error_y0[n] = {4.08e+02, 3.36e+02, 4.15e+02, 3.15e+02, 3.13e+02,3.10e+02, 3.08e+02};
	
	//for each region 1 graph
	TGraphErrors *gr2 = new TGraphErrors(n,eb,NQCD_2,error_x,error_y2);
	gr2->SetTitle("N_{QCD}^{(2)} vs e_{b}");
	gr2->SetMarkerColor(4);
	gr2->SetMarkerStyle(21);
	gr2->GetXaxis()->SetTitle("e_{b}");
	gr2->GetYaxis()->SetTitle("N_{QCD}^{(2)}");
	gr2->GetYaxis()->SetTitleOffset(1.4);
	TCanvas *c2 = new TCanvas("2btag", "2btag", 800, 600);
	gr2->Draw("AP");


	TGraphErrors *gr1 = new TGraphErrors(n, eb, NQCD_1, error_x, error_y1);
	gr1->SetTitle("N_{QCD}^{(1)} vs e_{b}");
	gr1->SetMarkerColor(4);
	gr1->SetMarkerStyle(21);
	gr1->GetXaxis()->SetTitle("e_{b}");
	gr1->GetYaxis()->SetTitle("N_{QCD}^{(1)}");
	gr1->GetYaxis()->SetTitleOffset(1.5);
	TCanvas *c1 = new TCanvas("1btag", "1btag", 800,600);
	gr1->Draw("AP");


	TGraphErrors *gr0 = new TGraphErrors(n, eb, NQCD_0, error_x, error_y0);
	gr0->SetTitle("N_{QCD}^{(0)} vs e_{b}");
	gr0->SetMarkerColor(4);
	gr0->SetMarkerStyle(21);
	gr0->GetXaxis()->SetTitle("e_{b}");
	gr0->GetYaxis()->SetTitle("N_{QCD}^{(0)}");
	gr0->GetYaxis()->SetTitleOffset(1.3);
	TCanvas *c0 = new TCanvas("0btag", "0btag", 800,600);
	gr0->Draw("AP");
}