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

#include "../CMS_plots/CMS_lumi.C"
#include "../CMS_plots/tdrstyle.C"

using std::cin;
using std::cout;
using std::endl;

void graphTransferFactorsSpecific(TString year = "2016_preVFP", bool bEnriched = false, TString variable = "jetPt0");

void graphTransferFactors(TString year = "2016_preVFP", bool bEnriched = false)
{
  const int NVAR = 11;
  TString vars[] = {"mJJ", "ptJJ", "yJJ", "jetPt0", "jetPt1", "jetY0", "jetY1", "chi", "cosTheta_0", "cosTheta_1", "mTop"};
  for(int ivar =0; ivar<sizeof(vars)/sizeof(vars[0]); ivar++)
  {
    graphTransferFactorsSpecific(year, bEnriched, vars[ivar]);
  }
}


void graphTransferFactorsSpecific(TString year = "2016_preVFP", bool bEnriched = false, TString variable = "jetPt0")
{
  gStyle->SetOptStat(0);
  setTDRStyle();
  TH1F *hfData;
  TH1F *hfQCD;
  TFile *outF, *outFQCD;
  TString bEnr = "_HT300toInf_100";
  if(bEnriched) bEnr = "_bEnriched_HT200toInf_100";


	TFile *infData, *infDataReduced, *infQCD, *infQCDReduced;
	TString str=year;
	infData = TFile::Open(TString::Format("../MassFit/%s/Histo_Data_%s_100.root", str.Data(),str.Data()));
	infDataReduced = TFile::Open(TString::Format("../MassFit/%s/Histo_Data_%s_100_reduced_UnequalBinning.root", str.Data(),str.Data()));

	infQCD = TFile::Open(TString::Format("../MassFit/%s/Histo_QCD%s.root", str.Data(),bEnr.Data()));
	infQCDReduced = TFile::Open(TString::Format("../MassFit/%s/Histo_QCD%s_reduced_UnequalBinning.root",str.Data(),bEnr.Data()));

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
  cout<<"hintegralDataReduced[0]: "<<hDataReduced[0]->Integral()<<endl;
	cout<<"hintegralDataReduced[0] integral and error: "<<hDataReduced[0]->IntegralAndError(1, hData[0]->GetNbinsX(),temp)<<endl;
	cout<<"hintegralData[0]: "<<hData[0]->Integral()<<endl;

	hQCD[0] = (TH1F*)infQCD->Get(TString::Format("hWt_%s_%dbtag_expYield",variable.Data(),0));
	hQCDReduced[0] = (TH1F*)infQCDReduced->Get(TString::Format("hWt_%s_%dbtag_expYield",variable.Data(),0));
	cout<<"mcintegralReduced[0]: "<<hQCDReduced[0]->Integral()<<endl;
	cout<<"mcintegral[0]: "<<hQCD[0]->Integral()<<endl;

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
	
	hData[0] = (TH1F*)infData->Get(TString::Format("hWt_%s_%dbtag_expYield",variable.Data(),0));
	hDataReduced[0] = (TH1F*)infDataReduced->Get(TString::Format("hWt_%s_%dbtag_expYield",variable.Data(),0));
	hQCD[0] = (TH1F*)infQCD->Get(TString::Format("hWt_%s_%dbtag_expYield",variable.Data(),0));
	hQCDReduced[0] = (TH1F*)infQCDReduced->Get(TString::Format("hWt_%s_%dbtag_expYield",variable.Data(),0));

	hData[1] = (TH1F*)infData->Get(TString::Format("hWt_%s_%dbtag_expYield",variable.Data(),2));
	hDataReduced[1] = (TH1F*)infDataReduced->Get(TString::Format("hWt_%s_%dbtag_expYield",variable.Data(),2));
	hQCD[1] = (TH1F*)infQCD->Get(TString::Format("hWt_%s_%dbtag_expYield",variable.Data(),2));
	hQCDReduced[1] = (TH1F*)infQCDReduced->Get(TString::Format("hWt_%s_%dbtag_expYield",variable.Data(),2));


	tFactorData[0] = (hDataReduced[0]->Integral() / hData[0]->Integral());
	tFactorData[1] = (hDataReduced[1]->Integral() / hData[1]->Integral());
    
	Double_t intErrorData[2], intErrorDataReduced[2];
    
	Double_t integralData[2];
	integralData[0] = hData[0]->IntegralAndError(1, hData[0]->GetNbinsX(), intErrorData[0]);
	integralData[1] = hData[1]->IntegralAndError(1, hData[1]->GetNbinsX(), intErrorData[1]);
	
	Double_t integralDataReduced[2];
	integralDataReduced[0] = hDataReduced[0]->IntegralAndError(1, hDataReduced[0]->GetNbinsX(), intErrorDataReduced[0]);
	integralDataReduced[1] = hDataReduced[1]->IntegralAndError(1, hDataReduced[1]->GetNbinsX(), intErrorDataReduced[1]);

    tFactorDataError[0] = TMath::Sqrt(TMath::Power(intErrorDataReduced[0]/integralData[0],2) + TMath::Power(integralDataReduced[0]*intErrorData[0]/TMath::Power(integralData[0],2),2) );
	tFactorDataError[1] = TMath::Sqrt(TMath::Power(intErrorDataReduced[1]/integralData[1],2) + TMath::Power(integralDataReduced[1]*intErrorData[1]/TMath::Power(integralData[1],2),2) );

    tFactorQCD[0] = (hQCDReduced[0]->Integral() / hQCD[0]->Integral());
	tFactorQCD[1] = (hQCDReduced[1]->Integral() / hQCD[1]->Integral());
	cout<<tFactorData[0]<<"±"<<tFactorDataError[0]<<endl;;
    cout<<tFactorQCD[0]<<endl;
	cout<<tFactorData[1]<<"±"<<tFactorDataError[1]<<endl;;
    cout<<tFactorQCD[1]<<endl;
	
	Double_t intError[2], intErrorReduced[2];
	Double_t integral[2], integralReduced[2];
	integral[0] = hQCD[0]->IntegralAndError(1, hQCD[0]->GetNbinsX(), intError[0]);
	integral[1] = hQCD[1]->IntegralAndError(1, hQCD[1]->GetNbinsX(), intError[1]);
	
	integralReduced[0] = hQCDReduced[0]->IntegralAndError(1, hQCDReduced[0]->GetNbinsX(), intErrorReduced[0]);
	integralReduced[1] = hQCDReduced[1]->IntegralAndError(1, hQCDReduced[1]->GetNbinsX(), intErrorReduced[1]);

    tFactorQCDError[0] = TMath::Sqrt(TMath::Power(intErrorReduced[0]/integral[0],2) + TMath::Power(integralReduced[0]*intError[0]/TMath::Power(integral[0],2),2) );
	tFactorQCDError[1] = TMath::Sqrt(TMath::Power(intErrorReduced[1]/integral[1],2) + TMath::Power(integralReduced[1]*intError[1]/TMath::Power(integral[1],2),2) );
	

	////////////////////////////////////////////////
	// rYield = (integralReduced[1] / integral[1]) / (integralReduced[0] / integral[0]);
  	Double_t rYield = tFactorQCD[1] / tFactorQCD[0];
	  
	// rYield_correction = (integralDataReduced[0] / integralData[0]);

	Double_t rYield_correction = tFactorData[0];

  	Double_t rYield_error = TMath::Sqrt(TMath::Power(integral[0] / (integralReduced[0] * integral[1]) * intErrorReduced[1], 2) +
                             TMath::Power(integralReduced[1] / (integralReduced[0] * integral[1]) * intError[0], 2) +
                             TMath::Power((integralReduced[1] * integral[0]) / (TMath::Power(integralReduced[0], 2) * integral[1]) * intErrorReduced[0], 2) +
                             TMath::Power((integralReduced[1] * integral[0]) / (integralReduced[0] * TMath::Power(integral[1], 2)) * intError[1], 2));

 	Double_t rYield_correction_error = TMath::Sqrt(TMath::Power(1 / integralData[0] * intErrorDataReduced[0], 2) +
                                                 TMath::Power(integralDataReduced[0] / TMath::Power(integralData[0], 2) * intErrorData[0], 2));
  	Double_t error = TMath::Sqrt(TMath::Power(rYield_correction * rYield_error, 2) + TMath::Power(rYield * rYield_correction_error, 2));
 	rYield = rYield * rYield_correction;
  	std::cout << variable << ": " << rYield << " " << error << std::endl;


	
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
  	  outF = new TFile(TString::Format("../MassFit/%s/Ryield/TransferFactor%s_%s.root",year.Data(), bEnr.Data(), variable.Data()), "RECREATE");
  	  hfData->SetTitle(TString::Format("R_{yield} %s %s", year.Data(),variable.Data()));
	  hfData->GetYaxis()->SetTitle("#frac{N_{Region}}{N_{Ext.Region}}");
      hfData->SetMarkerStyle(marker[i]);
	  hfData->SetMarkerColor(col[i]);
	  hfData->SetLineColor(col[i]);
	  hfData->GetYaxis()->SetTitleOffset(1.25);
      hfData->GetYaxis()->SetRangeUser(0.1,0.8);
	  can = new TCanvas(TString::Format("can_%s_%s", year.Data(),variable.Data()),TString::Format("can_%s_%s", year.Data(),variable.Data()),800,600);
	  hfData->Draw("hist e");
	  int iPeriod = 13;
      int iPos = 0;
      writeExtraText=true;
      CMS_lumi(can, iPeriod, iPos);
	  can->Print(TString::Format("../MassFit/%s/Ryield/TransferFactor%s_%s.pdf",year.Data(), 
	  							bEnr.Data(),variable.Data()),"pdf");

	  hfData->Write("dataTransferFactor");
	  //now plot them QCD
  	  hfQCD->SetTitle(TString::Format("R_{yield} %s %s (Closure Test)", year.Data(),variable.Data()));
	  hfQCD->GetYaxis()->SetTitle("#frac{N_{Region}}{N_{Ext.Region}}");
      hfQCD->SetMarkerStyle(marker[i]);
	  hfQCD->SetMarkerColor(col[i]);
	  hfQCD->SetLineColor(col[i]);
	  hfQCD->GetYaxis()->SetTitleOffset(1.25);
      hfQCD->GetYaxis()->SetRangeUser(0.1,0.5);
	  canQCD = new TCanvas(TString::Format("canQCD_%s_%s", year.Data(),variable.Data()),TString::Format("canQCD_%s_%s", year.Data(),variable.Data()),800,600);
	  hfQCD->Draw("hist e");
	  
	  hfQCD->Write("ClosureTest_TransferFactor");
      writeExtraText=true;
      CMS_lumi(canQCD, iPeriod, iPos);
	  canQCD->Print(TString::Format("../MassFit/%s/Ryield/TransferFactor_ClosureIntegral%s_%s.pdf",
	  				year.Data(), bEnr.Data(),variable.Data()),"pdf");


	  TH1F *hCorrectedRyield = new TH1F("CorrectedRYield", "CorrectedRYield", 1, 0, 1);
	  hCorrectedRyield->SetBinContent(1, rYield);
      hCorrectedRyield->SetBinError(1, error);

	  hCorrectedRyield->Write("CorrectedRYield");

}//eof
