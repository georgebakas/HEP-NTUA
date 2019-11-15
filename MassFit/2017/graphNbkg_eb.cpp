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
	float eb[n] = {0.4, 0.5,0.54097, 0.6, 0.605622, 0.7};
	float error_x[n] = {0, 0,1.08e-02, 0,0, 0};

	//Signal Region (2btag)
	float NQCD_2[n]   = {2.9028e+03, 2.5498e+03, 2.4652e+03, 2.3919e+03, 2.3906e+03 , 1.59e+03};
	float error_y2[n] = {83,127,132,136,145,768};
	
	//1btag region
	float NQCD_1[n] = {3.3287e+04, 3.6013e+04, 3.7067e+04, 3.8502e+04, 3.8634e+04, 3.9360e+04};
	float error_y1[n] = {247,325,457,526,226,247};

	//Control Region (0btag)
	float NQCD_0[n] = {1.4992e+05, 1.5369e+05, 1.5483e+05, 1.5603e+05, 1.5612e+05, 1.5704e+05};
	float error_y0[n] = {500, 439,541, 419,378,413};
	


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