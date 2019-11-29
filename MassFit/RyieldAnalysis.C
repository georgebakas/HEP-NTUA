/*
	This is a piece of code dedicated for Finding Ryield
*/

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

using namespace std;
#include "TemplateConstants.h"

using std::cin;
using std::cout;
using std::endl;	


void RyieldAnalysis(TString year = "2016", bool free_eb = true)
{
	initFilesMapping(free_eb);
	//we need both reduced and extended files for Data and ttbar MC
	TFile *dataFiles[2], *mcFiles[2];
	TString red = "";
	for(int i =0; i<sizeof(dataFiles)/sizeof(dataFiles[0]); i++)
	{
		if(i==1) red = "_reduced";
		//cout<<TString::Format("%s/Histo_Data_%s_100%s.root",year.Data(), red.Data())<<endl;
		dataFiles[i] = TFile::Open(TString::Format("%s/Histo_Data_%s_100%s.root",year.Data(),year.Data(), red.Data()));
		mcFiles[i] = TFile::Open(TString::Format("%s/Histo_TT_Mtt-700toInf_100%s.root"	,year.Data(), red.Data()));
	} 
	//all of the th1's are in the CR
	TH1F *hDExt, *hDRed;
	TH1F *hMCExt, *hMCRed;

	hDExt = (TH1F*)dataFiles[0]->Get("hWt_mTop_0btag_expYield");//hWt_mTop_0btag_expYield
	hDRed = (TH1F*)dataFiles[1]->Get("hWt_mTop_0btag_expYield");

	hMCExt = (TH1F*)mcFiles[0]->Get("hWt_mTop_0btag_expYield");
	hMCRed = (TH1F*)mcFiles[1]->Get("hWt_mTop_0btag_expYield");

	//Find R0 and R1 where R0 = DRed / DExt and R1 = (DRed - MCRed) /(DExt - MCExt)

	float R[2];

	R[0] = hDRed->Integral() / hDExt->Integral();
	R[1] = (hDRed->Integral() - hMCRed->Integral() )/ (hDExt->Integral() - hMCExt->Integral());

	cout<<"(R[0]-R[1])/R[0] = "<< (R[0]-R[1])/R[0]<<endl;
	cout<<"----"<<year<<"----"<<endl;
	cout<<"R0 (just data) = "<<R[0]<<endl;
	cout<<"R1 (taking MC into account) = "<<R[1]<<endl;
	//now 2nd step is to find Nqcd, reduced in SR
	float NQCD_reduced_CR = hDRed->Integral() - hMCRed->Integral();

	float NQCD_reduced_SR = (Nbkg2Constants[TString::Format("Nbkg%s", year.Data())]/Nbkg0Constants[TString::Format("Nbkg%s", year.Data())]) * NQCD_reduced_CR ;

	cout<<"NQCD in Reduced (SR): "<<NQCD_reduced_SR<<endl;


}