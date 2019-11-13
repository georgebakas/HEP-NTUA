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
	float eb[n] = {0.4, 0.5, 0.6, 0.60786, 0.629909, 0.7};
	float error_x[n] = {0,0, 1.17e-02, 0,0, 0};

	//Signal Region (2btag)
	float NQCD_2[n] = {3.9683e+03,3.2907e+03, 3.0841e+03,3.0802e+03,3.0552e+03 ,3.0336e+03};
	float error_y2[n] = {84,136,140,145,142,87};
	
	//1btag region
	float NQCD_1[n] = {2.6065e+04, 2.8069e+04, 3.059e+04, 3.0796e+04,3.1318e+04 ,3.2638e+04};
	float error_y1[n] = {203,316,271,461,275,232};

	//Control Region (0btag)
	float NQCD_0[n] = {8.2355e+04, 8.5484e+04, 8.7999e+04, 8.8132e+04,8.8492e+04 ,8.9367e+04};
	float error_y0[n] = {309,372,323,448,370,308};
	
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