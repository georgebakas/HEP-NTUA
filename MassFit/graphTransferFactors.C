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
  TH1F *hfData[3]; 
  TH1F *hfQCD[3]; 
  TFile *outF[3], *outFQCD[3];
  for(int y = 2016; y<2019; y++)
  {		
	TFile *infData, *infDataReduced, *infQCD, *infQCDReduced;

	infData = TFile::Open(TString::Format("%d/Histo_Data_%d_100.root", y,y));
	infDataReduced = TFile::Open(TString::Format("%d/Histo_Data_%d_100_reduced.root", y,y));
	infQCD = TFile::Open(TString::Format("%d/Histo_QCD_HT300toInf_100.root", y));
	infQCDReduced = TFile::Open(TString::Format("%d/Histo_QCD_HT300toInf_100_reduced.root",y));

	TH1F *hData[3], *hDataReduced[3];
	TH1F *hQCD[3], *hQCDReduced[3];
	hfData[y-2016] = new TH1F(TString::Format("hTransf_%d",y), TString::Format("hTransf_%d",y),3,0,3);	
    hfQCD[y-2016]  = new TH1F(TString::Format("hTransfClosure_%d",y), TString::Format("hTransfClosure_%d",y),3,0,3);
     	
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
		
		hfData[y-2016]->SetBinContent(i+1, tFactorData[i]);
     	hfData[y-2016]->SetBinError(i+1, tFactorDataError[i]);
  	  	hfData[y-2016]->GetXaxis()->SetBinLabel(i+1,names[i].Data());
  	  	cout<<hfData[y-2016]->GetBinError(i+1)<<endl;
  	  	hfQCD[y-2016]->SetBinContent(i+1, tFactorQCD[i]);
     	hfQCD[y-2016]->SetBinError(i+1, tFactorQCDError[i]);
  	  	hfQCD[y-2016]->GetXaxis()->SetBinLabel(i+1,names[i].Data());
	}
		
   //GetXaxis()->SetBinLabel(i+1,histoNames[i].Data());
    
	} //end of all years
    std::vector<int> col ={4,2,8};
	std::vector<int> marker = {21,20,23};
	TCanvas *can[3], *canQCD[3];
  	

  	for(int i =0; i<3; i++)
  	{
  	  int year = 2016+i;
  	  //now plot them Data
  	  outF[i] = new TFile(TString::Format("%d/TransferFactor.root",year), "RECREATE");	  
  	  hfData[i]->SetTitle(TString::Format("R_{yield} transfer factor %d", year));
	  hfData[i]->GetYaxis()->SetTitle("#frac{N_{Region}}{N_{Ext.Region}}");
      hfData[i]->SetMarkerStyle(marker[i]);
	  hfData[i]->SetMarkerColor(col[i]);
	  hfData[i]->SetLineColor(col[i]);
	  hfData[i]->GetYaxis()->SetTitleOffset(1.25);
      hfData[i]->GetYaxis()->SetRangeUser(0.1,0.8);
	  can[i] = new TCanvas(TString::Format("can_%d", year),TString::Format("can_%d", year),800,600);
	  hfData[i]->Draw("hist e");
	  can[i]->Print(TString::Format("%d/Ryield/TransferFactor.pdf",year),"pdf");

	  hfData[i]->Write("dataTransferFactor");
	  //now plot them QCD
  	  hfQCD[i]->SetTitle(TString::Format("R_{yield} transfer factor %d (Closure Test)", year));
	  hfQCD[i]->GetYaxis()->SetTitle("#frac{N_{Region}}{N_{Ext.Region}}");
      hfQCD[i]->SetMarkerStyle(marker[i]);
	  hfQCD[i]->SetMarkerColor(col[i]);
	  hfQCD[i]->SetLineColor(col[i]);
	  hfQCD[i]->GetYaxis()->SetTitleOffset(1.25);
      hfQCD[i]->GetYaxis()->SetRangeUser(0.1,0.5);
	  canQCD[i] = new TCanvas(TString::Format("canQCD_%d", year),TString::Format("canQCD_%d", year),800,600);
	  hfQCD[i]->Draw("hist e");
	  canQCD[i]->Print(TString::Format("%d/Ryield/TransferFactor_ClosureIntegral.pdf",year),"pdf");

	  hfQCD[i]->Write("ClosureTest_TransferFactor");
    }

	

   

		



}//eof

