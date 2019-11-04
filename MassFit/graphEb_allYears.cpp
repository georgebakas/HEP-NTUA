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

void graphEb_allYears()
{
	const int n = 7;
	//2016
	/*
	//Signal Region (2btag)
	float NQCD_2[n] = {4.1645e+03,3.2985e+03, 2.9980e+03,2.8784e+03,2.84e+03 ,2.8091e+03, 2.9322e+03};
	float error_y2[n] = {1.36e+02, 1.37e+02, 1.43e+02, 1.38e+02,1.55e+02 ,1.41e+02, 1.45e+02};
	
	//1btag region
	float NQCD_1[n] = {2.6618e+04, 2.8082e+04, 2.8973e+04, 2.9865e+04,3.0542e+04 ,3.2121e+04, 3.4154e+04};
	float error_y1[n] = {3.80e+02, 3.20e+02, 3.93e+02, 2.76e+02,2.62e+02 ,2.36e+02, 2.13e+02};

	//Control Region (0btag)
	float NQCD_0[n] = {8.1973e+04, 8.5459e+04, 8.7019e+04, 8.7826e+04,8.832e+04 ,8.9197e+04, 9.9892e+04};
	float error_y0[n] = {4.08e+02, 3.36e+02, 4.15e+02, 3.15e+02, 3.13e+02,3.10e+02, 3.08e+02};
	*/

	float eb16[n] = {0.4, 0.5, 0.56029, 0.6, 0.629909, 0.7, 0.8};
	float error_x16[n] = {0,0, 1.17e-02, 0,0, 0, 0};
	
	float NQCD16[3][n] = {
		{8.1973e+04, 8.5459e+04, 8.7019e+04, 8.7826e+04,8.832e+04 ,8.9197e+04, 9.9892e+04}, //0btag
		{2.6618e+04, 2.8082e+04, 2.8973e+04, 2.9865e+04,3.0542e+04 ,3.2121e+04, 3.4154e+04}, //1btag
		{4.1645e+03,3.2985e+03, 2.9980e+03,2.8784e+03,2.84e+03 ,2.8091e+03, 2.9322e+03} //2btag
	};

	float errory16[3][n]= {
		{4.08e+02, 3.36e+02, 4.15e+02, 3.15e+02, 3.13e+02,3.10e+02, 3.08e+02}, //0btag
		{3.80e+02, 3.20e+02, 3.93e+02, 2.76e+02,2.62e+02 ,2.36e+02, 2.13e+02}, //1btag
		{1.36e+02, 1.37e+02, 1.43e+02, 1.38e+02,1.55e+02 ,1.41e+02, 1.45e+02} //2btag
	};

	//2017
	
	/*
	//Signal Region (2btag)
	NQCD17[2]   = {3.0167e+03, 2.503e+03, 2.4147e+03, 2.21e+03, 2.2031e+03 , 1.8895e+03, 0.708e+03};
	errory17[2]  = {1.30e+02,   1.41e+02, 1.32e+02, 1.34e+02, 1.36e+02 , 0.80e+02, 0.72e+02};
	
	//1btag region
	NQCD17[1]  = {3.511e+04, 3.6754e+04, 3.7109e+04, 3.9728e+04, 3.9865e+04, 4.177e+04, 4.2947e+04};
	errory17[1] = {406, 456, 331, 208, 287, 222, 221};

	//Control Region (0btag)
	NQCD17[0]   = {1.574e+05, 1.612e+05, 1.617e+05, 1.6429e+05, 1.6439e+05, 1.656e+05, 1.6621e+05};
	errory17[0] = {487, 588, 433, 418, 418, 413, 413};
	*/
	float eb17[n] = {0.4, 0.48641, 0.5, 0.6, 0.605622, 0.7, 0.8};
	float error_x17[n] = {0, 1.10e-02,0, 0,0, 0, 0};

	float NQCD17[3][n] = {
		{1.574e+05, 1.612e+05, 1.617e+05, 1.6429e+05, 1.6439e+05, 1.656e+05, 1.6621e+05}, //0btag
		{3.511e+04, 3.6754e+04, 3.7109e+04, 3.9728e+04, 3.9865e+04, 4.177e+04, 4.2947e+04}, //1btag
		{3.0167e+03, 2.503e+03, 2.4147e+03, 2.21e+03, 2.2031e+03 , 1.8895e+03, 0.708e+03} //2btag

	};

	float errory17[3][n]= {
		{487, 588, 433, 418, 418, 413, 413}, //0btag
		{406, 456, 331, 208, 287, 222, 221}, //1btag
		{1.30e+02,   1.41e+02, 1.32e+02, 1.34e+02, 1.36e+02 , 0.80e+02, 0.72e+02} //2btag
	};

	TMultiGraph *mg[3];
	TGraphErrors *gr16[3], *gr17[3];
	TCanvas *can[3];
	TLegend *leg[3];
	//scale in relation to the 1st element?? could be ok 
	for(int i=0; i<3; i++)
	{	
		float norm16 = NQCD16[i][0];
		float norm17 = NQCD17[i][0];
		float errNorm16 = errory16[i][0];
		float errNorm17 = errory17[i][0];
		for(int j=0; j<n; j++)
		{
			NQCD16[i][j] = NQCD16[i][j]/norm16;
			NQCD17[i][j] = NQCD17[i][j]/norm17;
			errory16[i][j] = 0;
			errory17[i][j] = 0;
		}
	}
	for(int i=0; i<3; i++)
	{	

		mg[i] = new TMultiGraph();
		gr16[i] = new TGraphErrors(n,eb16,NQCD16[i],error_x16,errory16[i]);
		gr17[i] = new TGraphErrors(n,eb17,NQCD17[i],error_x17,errory17[i]);
		gr16[i]->SetTitle(TString::Format("N_{QCD}^{(%d)} vs e_{b}", i));
		gr16[i]->SetMarkerColor(4);
		gr16[i]->SetLineColor(4);
		gr16[i]->SetMarkerStyle(21);
		gr16[i]->GetXaxis()->SetTitle("e_{b}");
		gr16[i]->GetYaxis()->SetTitle(TString::Format("N_{QCD}^{(%d)}", i));
		gr16[i]->GetYaxis()->SetTitleOffset(1.4);
		gr17[i]->SetMarkerColor(2);
		gr17[i]->SetLineColor(2);
		gr17[i]->SetMarkerStyle(20);
		gr17[i]->GetXaxis()->SetTitle("e_{b}");
		gr17[i]->GetYaxis()->SetTitle(TString::Format("N_{QCD}^{(%d)}", i));
		gr17[i]->GetYaxis()->SetTitleOffset(1.4);
		mg[i]->Add(gr17[i]);
		mg[i]->Add(gr16[i]);
		can[i] = new TCanvas(TString::Format("can_%d", i), TString::Format("can_%d", i), 800,600);
		leg[i] = new TLegend(0.5, 0.6,0.7,0.8);
		leg[i]->AddEntry(gr16[i], "2016", "lp");
		leg[i]->AddEntry(gr17[i], "2017", "lp");
		can[i]->cd();
		mg[i]->Draw("AP");
		leg[i]->Draw();
	}







	
}