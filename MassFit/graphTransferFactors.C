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

void graphTransferFactors(TString year = "2016")
{
  gStyle->SetOptStat(0);
  TH1F *hfData; 
  TH1F *hfQCD; 
  TFile *outF, *outFQCD;
	
	TFile *infData, *infDataReduced, *infQCD, *infQCDReduced;
	TString str=year;
	infData = TFile::Open(TString::Format("%s/Histo_Data_%s_100.root", str.Data(),str.Data()));
	infDataReduced = TFile::Open(TString::Format("%s/Histo_Data_%s_100_reduced.root", str.Data(),str.Data()));

	infQCD = TFile::Open(TString::Format("%s/Histo_QCD_HT300toInf_100.root", str.Data()));
	infQCDReduced = TFile::Open(TString::Format("%s/Histo_QCD_HT300toInf_100_reduced.root",str.Data()));

	TH1F *hData[3], *hDataReduced[3];
	TH1F *hQCD[3], *hQCDReduced[3];
	hfData = new TH1F(TString::Format("hTransf_%s",year.Data()), TString::Format("hTransf_%s",year.Data()),3,0,3);	
    hfQCD  = new TH1F(TString::Format("hTransfClosure_%s",year.Data()), TString::Format("hTransfClosure_%s",year.Data()),3,0,3);
     	
	//TString histoNames[3] = {"SR", "1Btag", "CR"};
	float tFactorData[3], tFactorQCD[3], tFactorQCDError[3], tFactorDataError[3];
	std::vector<TString> names = {"0btag", "1btag", "2btag"};
    
    float x[3] = {0,1,2};

	for(int i =0; i<3; i++)
	{
		hData[i] = (TH1F*)infData->Get(TString::Format("hWt_mTop_%dbtag_expYield",i));
		hDataReduced[i] = (TH1F*)infDataReduced->Get(TString::Format("hWt_mTop_%dbtag_expYield",i));

		hQCD[i] = (TH1F*)infQCD->Get(TString::Format("hWt_mTop_%dbtag_expYield",i));
		hQCDReduced[i] = (TH1F*)infQCDReduced->Get(TString::Format("hWt_mTop_%dbtag_expYield",i));

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
	for(int i =0; i<3; i++)
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
  	if(year.EqualTo("2017") || year.EqualTo("2017_Loose")) i = 1;
  	else if(year.EqualTo("2018") || year.EqualTo("2018_Loose")) i = 2;
  	
  	  //int year = 2016+i;
  	  //TString year = "2016";
  	  //now plot them Data
  	  outF = new TFile(TString::Format("%s/TransferFactor.root",year.Data()), "RECREATE");	  
  	  hfData->SetTitle(TString::Format("R_{yield} transfer factor %s", year.Data()));
	  hfData->GetYaxis()->SetTitle("#frac{N_{Region}}{N_{Ext.Region}}");
      hfData->SetMarkerStyle(marker[i]);
	  hfData->SetMarkerColor(col[i]);
	  hfData->SetLineColor(col[i]);
	  hfData->GetYaxis()->SetTitleOffset(1.25);
      hfData->GetYaxis()->SetRangeUser(0.1,0.8);
	  can = new TCanvas(TString::Format("can_%s", year.Data()),TString::Format("can_%s", year.Data()),800,600);
	  hfData->Draw("hist e");
	  can->Print(TString::Format("%s/Ryield/TransferFactor.pdf",year.Data()),"pdf");

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
	  canQCD->Print(TString::Format("%s/Ryield/TransferFactor_ClosureIntegral.pdf",year.Data()),"pdf");

	  hfQCD->Write("ClosureTest_TransferFactor");
  

	

   

		



}//eof

