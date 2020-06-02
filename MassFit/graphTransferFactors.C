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

void graphTransferFactors(TString year = "2016", bool bEnriched = false)
{
  gStyle->SetOptStat(0);
  TH1F *hfData;
  TH1F *hfQCD;
  TFile *outF, *outFQCD;
  TString bEnr = "_HT300toInf_100";
  if(bEnriched) bEnr = "_bEnriched_HT200toInf_100";


	TFile *infData, *infDataReduced, *infQCD, *infQCDReduced;
	TString str=year;
	infData = TFile::Open(TString::Format("%s/Histo_Data_%s_100.root", str.Data(),str.Data()));
	infDataReduced = TFile::Open(TString::Format("%s/Histo_Data_%s_100_reduced.root", str.Data(),str.Data()));

	infQCD = TFile::Open(TString::Format("%s/Histo_QCD%s.root", str.Data(),bEnr.Data()));
	infQCDReduced = TFile::Open(TString::Format("%s/Histo_QCD%s_reduced.root",str.Data(),bEnr.Data()));

	TH1F *hData[2], *hDataReduced[2];
	TH1F *hQCD[2], *hQCDReduced[2];
	hfData = new TH1F(TString::Format("hTransf_%s",year.Data()), TString::Format("hTransf_%s",year.Data()),2,0,2);
    hfQCD  = new TH1F(TString::Format("hTransfClosure_%s",year.Data()), TString::Format("hTransfClosure_%s",year.Data()),2,0,2);

	//TString histoNames[3] = {"SR", "1Btag", "CR"};
	float tFactorData[2], tFactorQCD[2], tFactorQCDError[2], tFactorDataError[2];
	std::vector<TString> names = {"0btag", "2btag"};

    float x[2] = {0,2};

	hData[0] = (TH1F*)infData->Get(TString::Format("hWt_mTop_%dbtag_expYield",0));
	hDataReduced[0] = (TH1F*)infDataReduced->Get(TString::Format("hWt_mTop_%dbtag_expYield",0));
	cout<<"hData_0btag_reduced: "<<hDataReduced[0]->Integral()<<endl;
	cout<<"hData_0btag_extended: "<<hData[0]->Integral()<<endl;

	hQCD[0] = (TH1F*)infQCD->Get(TString::Format("hWt_mTop_%dbtag_expYield",0));
	hQCDReduced[0] = (TH1F*)infQCDReduced->Get(TString::Format("hWt_mTop_%dbtag_expYield",0));
	cout<<"mcQCD_0btag_reduced: "<<hQCDReduced[0]->Integral()<<endl;
	cout<<"mcQCD_0btag_extended: "<<hQCD[0]->Integral()<<endl;

	tFactorData[0] = (hDataReduced[0]->GetEntries() / hData[0]->GetEntries());
	tFactorDataError[0] = TMath::Sqrt((hDataReduced[0]->GetEntries()*(hData[0]->GetEntries() + hDataReduced[0]->GetEntries()))/ TMath::Power(hData[0]->GetEntries(),3));
	tFactorQCD[0] = (hQCDReduced[0]->Integral() / hQCD[0]->Integral());
	//cout<<tFactorData[i]<<endl;
	Double_t intError, intErrorReduced;
	float integral = hQCD[0]->IntegralAndError(1, hQCD[0]->GetNbinsX(), intError);
	float integralReduced = hQCDReduced[0]->IntegralAndError(1, hQCDReduced[0]->GetNbinsX(), intErrorReduced);

	tFactorQCDError[0] = TMath::Sqrt(TMath::Power(intErrorReduced,2)/TMath::Power(hQCD[0]->Integral(),2) + TMath::Power(hQCDReduced[0]->Integral(),2)
			*TMath::Power(intError,2)/TMath::Power(hQCD[0]->Integral(),4)  );


	for(int i =1; i<2; i++)
	{
		hData[i] = (TH1F*)infData->Get(TString::Format("hWt_mTop_%dbtag_expYield",i+1));
		hDataReduced[i] = (TH1F*)infDataReduced->Get(TString::Format("hWt_mTop_%dbtag_expYield",i+1));
		//cout<<"hData_2btag_reduced: "<<hDataReduced[i]->GetEntries()<<endl;
		//cout<<"hData_2btag_extended: "<<hData[i]->GetEntries()<<endl;

		hQCD[i] = (TH1F*)infQCD->Get(TString::Format("hWt_mTop_%dbtag_expYield",i+1));
		hQCDReduced[i] = (TH1F*)infQCDReduced->Get(TString::Format("hWt_mTop_%dbtag_expYield",i+1));
		cout<<"mcQCD_2btag_reduced: "<<hQCDReduced[i]->Integral()<<endl;
		cout<<"mcQCD_2btag_extended: "<<hQCD[i]->Integral()<<endl;

		tFactorData[i] = (hDataReduced[i]->GetEntries() / hData[i]->GetEntries());
		tFactorDataError[i] = TMath::Sqrt((hDataReduced[i]->GetEntries()*(hData[i]->GetEntries() + hDataReduced[i]->GetEntries()))/ TMath::Power(hData[i]->GetEntries(),3));
		tFactorQCD[i] = (hQCDReduced[i]->Integral() / hQCD[i]->Integral());
		//cout<<tFactorData[i]<<endl;
		Double_t intError, intErrorReduced;
		float integral = hQCD[i]->IntegralAndError(1, hQCD[i]->GetNbinsX(), intError);
		float integralReduced = hQCDReduced[i]->IntegralAndError(1, hQCDReduced[i]->GetNbinsX(), intErrorReduced);

		tFactorQCDError[i] = TMath::Sqrt(TMath::Power(intErrorReduced,2)/TMath::Power(hQCD[i]->Integral(),2) + TMath::Power(hQCDReduced[i]->Integral(),2)
				*TMath::Power(intError,2)/TMath::Power(hQCD[i]->Integral(),4)  );
	}
	for(int i =0; i<2; i++)
	{
		hfData->SetBinContent(i+1, tFactorData[i]);
     	hfData->SetBinError(i+1, tFactorDataError[i]);
  	  	hfData->GetXaxis()->SetBinLabel(i+1,names[i].Data());

  	  	hfQCD->SetBinContent(i+1, tFactorQCD[i]);
     	hfQCD->SetBinError(i+1, tFactorQCDError[i]);
  	  	hfQCD->GetXaxis()->SetBinLabel(i+1,names[i].Data());
	}

   //GetXaxis()->SetBinLabel(i+1,histoNames[i].Data());

    std::vector<int> col ={4,2,8};
	std::vector<int> marker = {21,20,23};
	TCanvas *can, *canQCD;


  	int i(0);
  	if(year.EqualTo("2017")) i = 1;
  	else if(year.EqualTo("2018")) i = 2;

  	  //int year = 2016+i;
  	  //TString year = "2016";
  	  //now plot them Data
  	  outF = new TFile(TString::Format("%s/TransferFactor%s.root",year.Data(), bEnr.Data()), "RECREATE");
  	  hfData->SetTitle(TString::Format("R_{yield} transfer factor %s", year.Data()));
	  hfData->GetYaxis()->SetTitle("#frac{N_{Region}}{N_{Ext.Region}}");
      hfData->SetMarkerStyle(marker[i]);
	  hfData->SetMarkerColor(col[i]);
	  hfData->SetLineColor(col[i]);
	  hfData->GetYaxis()->SetTitleOffset(1.25);
      hfData->GetYaxis()->SetRangeUser(0.1,0.8);
	  can = new TCanvas(TString::Format("can_%s", year.Data()),TString::Format("can_%s", year.Data()),800,600);
	  hfData->Draw("hist e");
	  can->Print(TString::Format("%s/Ryield/TransferFactor%s.pdf",year.Data(), bEnr.Data()),"pdf");

	  hfData->Write("dataTransferFactor");
	  //now plot them QCD
  	  hfQCD->SetTitle(TString::Format("R_{yield} transfer factor %s (Closure Test)", year.Data()));
	  hfQCD->GetYaxis()->SetTitle("#frac{N_{Region}}{N_{Ext.Region}}");
      hfQCD->SetMarkerStyle(marker[i]);
	  hfQCD->SetMarkerColor(col[i]);
	  hfQCD->SetLineColor(col[i]);
	  hfQCD->GetYaxis()->SetTitleOffset(1.25);
      hfQCD->GetYaxis()->SetRangeUser(0.1,0.5);
	  canQCD = new TCanvas(TString::Format("canQCD_%s", year.Data()),TString::Format("canQCD_%s", year.Data()),800,600);
	  hfQCD->Draw("hist e");
	  canQCD->Print(TString::Format("%s/Ryield/TransferFactor_ClosureIntegral%s.pdf",year.Data(), bEnr.Data()),"pdf");

	  hfQCD->Write("ClosureTest_TransferFactor");



}//eof
