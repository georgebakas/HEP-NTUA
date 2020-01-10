#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"

using std::cin;
using std::cout;
using std::endl;

std::vector<TString> histoNames;

void initHistoNames()
{
	histoNames.push_back("b_quarks"); //0
	histoNames.push_back("c_quarks"); //1
	histoNames.push_back("uds_quarks"); //2
	histoNames.push_back("gluons"); //3
}

void plotBTagEfficiencies()
{
	gStyle->SetOptStat(0);
	gStyle->SetPaintTextFormat("0.2f");
	initHistoNames();
	TFile *inf = TFile::Open("BTaggingEfficiency_massCutLoose.root");

	TH2F *hRecoParton[histoNames.size()];
	TH2F *hRecoParton_clones[histoNames.size()];
	TH2F *hParton[histoNames.size()];
	TH2F *hAll;
	TH2F *hbParton;

	hbParton = (TH2F*)inf->Get("hTagParton_Signal_histo_Mtt-700-1000_b_quarks");
	for(int ih =0; ih<histoNames.size(); ih++)
	{
		hRecoParton[ih] = (TH2F*)inf->Get(TString::Format("hTagRecoParton_Signal_histo_Mtt-700-1000_%s", histoNames[ih].Data()));
		hRecoParton_clones[ih] = (TH2F*)hRecoParton[ih]->Clone(TString::Format("hTagRecoParton_%s_Clone",histoNames[ih].Data()));

		hParton[ih] = (TH2F*)inf->Get(TString::Format("hTagParton_Signal_histo_Mtt-700-1000_%s", histoNames[ih].Data()));
		if(ih ==0) hAll = (TH2F*)hRecoParton[ih]->Clone("hAll");
		else hAll->Add(hRecoParton[ih]);
	}

	TCanvas *canDiv[histoNames.size()];
	//now work with the TH2 that need to be divided
	for(int i =0; i<histoNames.size(); i++)
	{
		canDiv[i] = new TCanvas(TString::Format("canvas_%d", i+1),TString::Format("canvas_%d", i+1), 800,600);
		canDiv[i]->cd();
		hRecoParton_clones[i]->Divide(hParton[i]);
		hRecoParton_clones[i]->Draw("text colz");
	}


	//for acceptance and efficiency
	float eb = hRecoParton[0]->Integral() / hbParton->Integral();
	float acceptance = hRecoParton[0]->Integral() / hAll->Integral();
	cout<<"btagging efficiency = "<<eb<<endl; 
	cout<<"btagging acceptance = "<<acceptance<<endl;

	TH1F *hAcceptance = new TH1F("hPurity", "hPurity", 4,0,5);

	float acc[histoNames.size()];
	for(int i =0; i<histoNames.size(); i++)
	{
		acc[i] = hRecoParton[i]->Integral() /hAll->Integral();
		hAcceptance->SetBinContent(i+1, acc[i]);
		hAcceptance->GetXaxis()->SetBinLabel(i+1,histoNames[i].Data());
	}
	TCanvas *canAcc = new TCanvas("canAcc", "canAcc", 800, 600);
	hAcceptance->SetFillColor(kRed);
	hAcceptance->SetFillStyle(3003);
	hAcceptance->Draw();

}
