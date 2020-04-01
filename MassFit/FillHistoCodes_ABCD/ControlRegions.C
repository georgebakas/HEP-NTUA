#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"

#include "TemplateConstants.h"
using std::cin;
using std::cout;
using std::endl;

/*
In this code I will try to test the ABCD method
We will have 4 Control Regions:
For both jets:
A. 2jets top tagged (Medium WP) and 0jets b-tagged (Loose WP)
B. 0jets top tagged (Loose WP)  and 0jets b-tagged (Loose WP)
C. 2jets top tagged (Medium WP) and 2jets b-tagged (Medium WP)
D. 0jets top tagged (Loose WP)  and 2jets b-tagged (Loose WP)

What we do:

Find for all files:
1. Data: Remove ttbar and subdominant bkg (both taken from MC) 

*/



void ControlRegions(TString year = "2016")
{

	TFile *ttFileNominal = TFile::Open(TString::Format("%s/Histo_TT_NominalMC_100_reduced.root", year.Data()));
	//TFile *ttFileMtt = TFile::Open(TString::Format("%s/Histo_TT_Mtt-700toInf_100_reduced.root", year.Data())); 

	TFile *infData = TFile::Open(TString::Format("%s/Histo_Data_%s_100_reduced.root", year.Data(), year.Data()));	
	TFile *infQCD = TFile::Open(TString::Format("%s/Histo_QCD_HT300toInf_100_reduced.root", year.Data()));
	TFile *infBkg = TFile::Open(TString::Format("%s/Histo_SubdominantBkgs_100_reduced.root", year.Data()));


	//get all the control regions and find the ttbar yield in every region
	const int controlRegions = 4;
	TString regions[controlRegions] = {"A", "B", "C", "D"};

	float totalYield[controlRegions];
	float ttFraction[controlRegions];
	float totalData[controlRegions];

	TH1F *hData[controlRegions], *hTTMtt[controlRegions], *hQCD[controlRegions], *hBkg[controlRegions], *hTTNominal[controlRegions];

	for(int ireg=0; ireg<controlRegions; ireg++)
	{
		hData[ireg] = (TH1F*)infData->Get(TString::Format("hWt_mTop_CR%s_expYield",regions[ireg].Data()));
		//hTTMtt[ireg] = (TH1F*)ttFileMtt->Get(TString::Format("hWt_mTop_CR%s_expYield",regions[ireg].Data()));
		hTTNominal[ireg] = (TH1F*)ttFileNominal->Get(TString::Format("hWt_mTop_CR%s_expYield",regions[ireg].Data()));
		hQCD[ireg]= (TH1F*)infQCD->Get(TString::Format("hWt_mTop_CR%s_expYield",regions[ireg].Data()));
		hBkg[ireg] = (TH1F*)infBkg->Get(TString::Format("hWt_mTop_CR%s_expYield",regions[ireg].Data()));

		totalYield[ireg] = hTTNominal[ireg]->Integral() + hQCD[ireg]->Integral() + hBkg[ireg]->Integral();
		totalData[ireg]  = hData[ireg]->Integral();
		ttFraction[ireg] = hTTNominal[ireg]->Integral() / totalData[ireg];

		
		cout<<"--------------"<<endl;
		cout<<"MC total: "<<totalYield[ireg]<<" Data: "<<totalData[ireg]<<endl;
		cout<<"ttbar: "<<hTTNominal[ireg]->Integral()<<endl;
		cout<<"The ttbar fraction in CR_"<<regions[ireg]<<" is: "<<ttFraction[ireg]<<endl;
		cout<<"MC/data: "<<totalYield[ireg]/totalData[ireg]<<endl;
	}
	cout<<"----------"<<endl;
	//now apply abcd method
	float qcdExp[controlRegions];
	for(int ireg=0; ireg<controlRegions; ireg++)
	{
		qcdExp[ireg] = hData[ireg]->Integral() - hTTNominal[ireg]->Integral() - hBkg[ireg]->Integral();
		cout<<"qcd from data region "<<regions[ireg]<<": "<<qcdExp[ireg]<<endl;
	}
	//is it consistent with abcd method??
	cout<<"-----------"<<endl;

	//expectation 
	cout<<"Expected qcd in signal from Data -subdominant bkg - ttbar  is: "<<qcdExp[2]<<endl;
	float CR_S = qcdExp[0] * qcdExp[3] / qcdExp[1];
	cout<<"From ABCD method: "<<CR_S<<endl;
	cout<<"From MC qcd expectation: "<<hQCD[2]->Integral()<<endl;



}
