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

void graphTransferFactorsSpecific(TString year = "2016", TString variable = "jetPt0", int mJJCut = 1200);

void graphTransferFactors_2TeV(TString year = "2016", int mJJCut = 1200)
{
  const int NVAR = 4;
  TString vars[] = {"chi", "cosTheta_0", "cosTheta_1", "mJJ"};
  for(int ivar =0; ivar<sizeof(vars)/sizeof(vars[0]); ivar++)
  {
    graphTransferFactorsSpecific(year, vars[ivar], mJJCut);
    //break;
  }
}


void graphTransferFactorsSpecific(TString year = "2016", TString variable = "chi", int mJJCut = 1200)
{
  gStyle->SetOptStat(0);
  TH1F *hfData;
  TH1F *hfQCD;
  TFile *outF, *outFQCD;


	TFile *infData, *infDataReduced, *infQCD, *infQCDReduced;
	TString str=year;
  //_reduced_ is for mTop cut --> [120,220]GeV!!!!!!
   //extended is mJJ > 1000 GeV
	infData = TFile::Open(TString::Format("%s/Histo_Data_%s_reduced_1000.root", str.Data(),str.Data()));
  //my Signal region, is the reduced with mJJ > %d GeV
	infDataReduced = TFile::Open(TString::Format("%s/Histo_Data_%s_reduced_%d.root", str.Data(),str.Data(), mJJCut));

	infQCD = TFile::Open(TString::Format("%s/Histo_QCD_HT300toInf_reduced_1000.root", str.Data())); //this is now my SR
	infQCDReduced = TFile::Open(TString::Format("%s/Histo_QCD_HT300toInf_reduced_%d.root",str.Data(), mJJCut)); //this is SR extension

	TH1F *hData[2], *hDataReduced[2];
	TH1F *hQCD[2], *hQCDReduced[2];
	hfData = new TH1F(TString::Format("hTransf_%s",year.Data()), TString::Format("hTransf_%s",year.Data()),2,0,2);
  hfQCD  = new TH1F(TString::Format("hTransfClosure_%s",year.Data()), TString::Format("hTransfClosure_%s",year.Data()),2,0,2);

	//TString histoNames[3] = {"SR", "1Btag", "CR"};
	float tFactorData[2], tFactorQCD[2], tFactorQCDError[2], tFactorDataError[2];
	std::vector<TString> names = {"0btag", "2btag"};

    float x[2] = {0,2};
/*
	hData[0] = (TH1F*)infData->Get(TString::Format("hWt_%s_%dbtag_expYield",variable.Data(),0));
	hDataReduced[0] = (TH1F*)infDataReduced->Get(TString::Format("hWt_%s_%dbtag_expYield",variable.Data(),0));
  double temp;
  cout<<"hData_0btag_reduced: "<<hDataReduced[0]->Integral()<<endl;
	cout<<"hData_0btag_reduced integral and error: "<<hDataReduced[0]->IntegralAndError(1, hData[0]->GetNbinsX(),temp)<<endl;
	cout<<"hData_0btag_extended: "<<hData[0]->Integral()<<endl;

	hQCD[0] = (TH1F*)infQCD->Get(TString::Format("hWt_%s_%dbtag_expYield",variable.Data(),0));
	hQCDReduced[0] = (TH1F*)infQCDReduced->Get(TString::Format("hWt_%s_%dbtag_expYield",variable.Data(),0));
	cout<<"mcQCD_0btag_reduced: "<<hQCDReduced[0]->Integral()<<endl;
	cout<<"mcQCD_0btag_extended: "<<hQCD[0]->Integral()<<endl;

	tFactorData[0] = (hDataReduced[0]->Integral() / hData[0]->Integral());
	tFactorDataError[0] = TMath::Sqrt((hDataReduced[0]->GetEntries()*(hData[0]->GetEntries() + hDataReduced[0]->GetEntries()))/ TMath::Power(hData[0]->GetEntries(),3));
	tFactorQCD[0] = (hQCDReduced[0]->Integral() / hQCD[0]->Integral());
	//cout<<tFactorData[i]<<endl;
	Double_t intError, intErrorReduced;
	float integral = hQCD[0]->IntegralAndError(1, hQCD[0]->GetNbinsX(), intError);
	float integralReduced = hQCDReduced[0]->IntegralAndError(1, hQCDReduced[0]->GetNbinsX(), intErrorReduced);

	tFactorQCDError[0] = TMath::Sqrt(TMath::Power(intErrorReduced,2)/TMath::Power(hQCD[0]->Integral(),2) + TMath::Power(hQCDReduced[0]->Integral(),2)
			*TMath::Power(intError,2)/TMath::Power(hQCD[0]->Integral(),4)  );

*/
  cout<<variable<<endl;
	for(int i =0; i<2; i++)
	{
		if(i==0)
    {
      hData[i] = (TH1F*)infData->Get(TString::Format("hWt_%s_%dbtag_expYield",variable.Data(),i));
		  hDataReduced[i] = (TH1F*)infDataReduced->Get(TString::Format("hWt_%s_%dbtag_expYield",variable.Data(),i));
      hQCD[i] = (TH1F*)infQCD->Get(TString::Format("hWt_%s_%dbtag_expYield",variable.Data(),i));
  		hQCDReduced[i] = (TH1F*)infQCDReduced->Get(TString::Format("hWt_%s_%dbtag_expYield",variable.Data(),i));

    }
    else
    {
      hData[i] = (TH1F*)infData->Get(TString::Format("hWt_%s_%dbtag_expYield",variable.Data(),i+1));
      hDataReduced[i] = (TH1F*)infDataReduced->Get(TString::Format("hWt_%s_%dbtag_expYield",variable.Data(),i+1));
      hQCD[i] = (TH1F*)infQCD->Get(TString::Format("hWt_%s_%dbtag_expYield",variable.Data(),i+1));
  		hQCDReduced[i] = (TH1F*)infQCDReduced->Get(TString::Format("hWt_%s_%dbtag_expYield",variable.Data(),i+1));

    }

		tFactorData[i] = (hDataReduced[i]->Integral() / hData[i]->Integral());
    Double_t intErrorData, intErrorDataReduced;
    Double_t integralData = hData[i]->IntegralAndError(1, hData[i]->GetNbinsX(), intErrorData);
		Double_t integralDataReduced = hDataReduced[i]->IntegralAndError(1, hDataReduced[i]->GetNbinsX(), intErrorDataReduced);

    tFactorDataError[i] = TMath::Sqrt(TMath::Power(intErrorDataReduced/integralData,2) + TMath::Power(integralDataReduced*intErrorData/TMath::Power(integralData,2),2) );

    tFactorQCD[i] = (hQCDReduced[i]->Integral() / hQCD[i]->Integral());
		cout<<tFactorData[i]<<endl;
    cout<<tFactorQCD[i]<<endl;
		Double_t intError, intErrorReduced;
		Double_t integral = hQCD[i]->IntegralAndError(1, hQCD[i]->GetNbinsX(), intError);
		Double_t integralReduced = hQCDReduced[i]->IntegralAndError(1, hQCDReduced[i]->GetNbinsX(), intErrorReduced);

    tFactorQCDError[i] = TMath::Sqrt(TMath::Power(intErrorReduced/integral,2) + TMath::Power(integralReduced*intError/TMath::Power(integral,2),2) );


    if(i==0)
    {
      cout<<"data_0btag_reduced: "<<integralDataReduced<<" ± "<<intErrorDataReduced<<endl;
		  cout<<"data_0btag_extended: "<<integralData<<" ± "<<intErrorData<<endl;
      cout<<"mcQCD_0btag_reduced: "<<integralReduced<<" ± "<<intErrorReduced<<endl;
		  cout<<"mcQCD_0btag_extended: "<<integral<<" ± "<<intError<<endl;
    }
    else
    {
      cout<<"mcQCD_2btag_reduced: "<<integralReduced<<" ± "<<intErrorReduced<<endl;
		  cout<<"mcQCD_2btag_extended: "<<integral<<" ± "<<intError<<endl;
    }
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
  	  outF = new TFile(TString::Format("%s/Ryield/TransferFactor_%s_%d.root",year.Data(), variable.Data(), mJJCut), "RECREATE");
  	  hfData->SetTitle(TString::Format("R_{yield} transfer factor %s %s", year.Data(),variable.Data()));
	     hfData->GetYaxis()->SetTitle("#frac{N_{Region}}{N_{Ext.Region}}");
      hfData->SetMarkerStyle(marker[i]);
	  hfData->SetMarkerColor(col[i]);
	  hfData->SetLineColor(col[i]);
	  hfData->GetYaxis()->SetTitleOffset(1.25);
      hfData->GetYaxis()->SetRangeUser(0,1);
	  can = new TCanvas(TString::Format("can_%s_%s", year.Data(),variable.Data()),TString::Format("can_%s_%s", year.Data(),variable.Data()),800,600);
	  hfData->Draw("hist e");
	  can->Print(TString::Format("%s/Ryield/TransferFactor_%s_%d.pdf",year.Data(),variable.Data(), mJJCut),"pdf");

	  hfData->Write("dataTransferFactor");
	  //now plot them QCD
  	  hfQCD->SetTitle(TString::Format("R_{yield} transfer factor %s %s(Closure Test)", year.Data(),variable.Data()));
	  hfQCD->GetYaxis()->SetTitle("#frac{N_{Region}}{N_{Ext.Region}}");
      hfQCD->SetMarkerStyle(marker[i]);
	  hfQCD->SetMarkerColor(col[i]);
	  hfQCD->SetLineColor(col[i]);
	  hfQCD->GetYaxis()->SetTitleOffset(1.25);
      hfQCD->GetYaxis()->SetRangeUser(0,1);
	  canQCD = new TCanvas(TString::Format("canQCD_%s_%s", year.Data(),variable.Data()),TString::Format("canQCD_%s_%s", year.Data(),variable.Data()),800,600);
	  hfQCD->Draw("hist e");
	  canQCD->Print(TString::Format("%s/Ryield/TransferFactor_ClosureIntegral_%s_%d.pdf",year.Data(), variable.Data(), mJJCut),"pdf");

	  hfQCD->Write("ClosureTest_TransferFactor");



}//eof
