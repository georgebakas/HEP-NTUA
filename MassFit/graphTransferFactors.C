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

void graphTransferFactors()
{
	TFile *inf[3];
	TH1F *htf[3];
	std::vector<int> col ={4,2,8};
	std::vector<int> marker = {21,20,23};
	TCanvas *can[3];
	for(int i =0; i<3; i++)
	{
		int year = 2016+i;
		inf[i] = TFile::Open(TString::Format("%d/TransferFactor.root", year));
		htf[i] = (TH1F*)inf[i]->Get("TransferFactor_hist");
		//htf->SetName("R_{yield} transfer factor");
		htf[i]->SetTitle(TString::Format("R_{yield} transfer factor %d", year));
		htf[i]->GetYaxis()->SetTitle("#frac{N_{Region}}{N_{Ext.Region}}");
		htf[i]->SetMarkerStyle(marker[i]);
		htf[i]->SetMarkerColor(col[i]);
		htf[i]->SetLineColor(col[i]);
		can[i] = new TCanvas(TString::Format("can_%d", year),TString::Format("can_%d", year),800,600);
		htf[i]->Draw("hist");
	}
}