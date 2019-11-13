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

#include "TemplateConstants.h"

using std::cin;
using std::cout;
using std::endl;

void graphNbkg_eb(TString year)
{
	initNqcdVsEb();
	//for each region 1 graph
	TGraphErrors *gr2 = new TGraphErrors(n,eb,NQCD[2],error_x,error_y[2]);
	gr2->SetTitle("N_{QCD}^{(2)} vs e_{b}");
	gr2->SetMarkerColor(2);
	gr2->SetMarkerStyle(20);
	gr2->GetXaxis()->SetTitle("e_{b}");
	gr2->GetYaxis()->SetTitle("N_{QCD}^{(2)}");
	gr2->GetYaxis()->SetTitleOffset(1.4);
	TCanvas *c2 = new TCanvas("2btag", "2btag", 800, 600);
	gr2->Draw("AP");


	TGraphErrors *gr1 = new TGraphErrors(n, eb, NQCD[1], error_x, error_y[1]);
	gr1->SetTitle("N_{QCD}^{(1)} vs e_{b}");
	gr1->SetMarkerColor(2);
	gr1->SetMarkerStyle(20);
	gr1->GetXaxis()->SetTitle("e_{b}");
	gr1->GetYaxis()->SetTitle("N_{QCD}^{(1)}");
	gr1->GetYaxis()->SetTitleOffset(1.5);
	TCanvas *c1 = new TCanvas("1btag", "1btag", 800,600);
	gr1->Draw("AP");


	TGraphErrors *gr0 = new TGraphErrors(n, eb, NQCD[0], error_x, error[0]);
	gr0->SetTitle("N_{QCD}^{(0)} vs e_{b}");
	gr0->SetMarkerColor(2);
	gr0->SetMarkerStyle(20);
	gr0->GetXaxis()->SetTitle("e_{b}");
	gr0->GetYaxis()->SetTitle("N_{QCD}^{(0)}");
	gr0->GetYaxis()->SetTitleOffset(1.3);
	TCanvas *c0 = new TCanvas("0btag", "0btag", 800,600);
	gr0->Draw("AP");
}