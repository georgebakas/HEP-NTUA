/*
	This is a piece of code dedicated for Signal Extraction:
	For every variable x_reco we need to extract the signal using the formula:

	S(x) = D(x) - Ryield * Nbkg * QCDShapeCorrection * Q(x) - B(x)
	or
	S(x) = D(x) - (Ryield * Nbkg) * QCDShapeCorrection * Q(x) - B(x)
	where:
	1.S(x) is the extracted signal
	2.Ryield is the NSR / NSRA taken from the file--> year/TransferFactor.root
	3.Nbkg is the number of QCD events in the SR taken from the simultaneous fit leaving the eb free

	or 2*3 = N_q(reduced) = (Nqcd_2_extended from fit / Nqcd_0_extended from fit ) * Nqcd_reduced_0,
			 where Nqcd_red_0 = Data_reduced_0 - TTMC_red_0
	4.QCDShapeCorrection is taken from the corrected shape (only for 2017 and 2018 for specific variables) from files: FitOutput.root
	5.Q(x) is taken from the CR from data
	6.B(x) is the signal region of the Subdominant bkg (exp yield) taken from MC
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

using std::cin;
using std::cout;
using std::endl;

#include "TemplateConstants.h"

TH1F *getRebinned(TH1F *h, float BND[], int N);
void SignalExtractionSpecific(TString year = "2016", TString variable = "jetPt0", bool free_eb = true);

TH1F *getRebinned(TH1F *h, float BND[], int N)
{
  TString name = TString(h->GetName())+"_Rebinned";
  TH1F *hNew = new TH1F(name, name, N, BND);

  for(int i=0;i<h->GetNbinsX();i++) {
    float x = h->GetBinCenter(i+1);
    float y = h->GetBinContent(i+1);
    float e = h->GetBinError(i+1);
    for(int j=0;j<hNew->GetNbinsX();j++) {
      float x1 = hNew->GetBinLowEdge(j+1);
      float x2 = x1+hNew->GetBinWidth(j+1);
      if ((x>x1) && (x<x2)) {
        float yNew = hNew->GetBinContent(j+1);
        float eNew = hNew->GetBinError(j+1);
        hNew->SetBinContent(j+1,yNew+y);
        hNew->SetBinError(j+1,sqrt(e*e+eNew*eNew));
        break;
      }
    }
  }
  return hNew;

}

void SignalExtraction(TString year)
{
	TString vars[] = {"mJJ", "ptJJ", "yJJ", "jetPt0", "jetPt1"};
	for(int i =0; i<sizeof(vars)/sizeof(vars[0]); i++)
	{
		SignalExtractionSpecific(year, vars[i], true);
	}
}


void SignalExtractionSpecific(TString year = "2016", TString variable = "jetPt0", bool free_eb = true)
{
	initFilesMapping(free_eb);
	gStyle->SetOptStat(0);
	//open the signal file: get D(x) and Q(x) for every variable
	TFile *infData = TFile::Open(TString::Format("%s/Histo_Data_%s_100.root", year.Data(), year.Data()));

	TH1F *hD = (TH1F*)infData->Get(TString::Format("hWt_%s_2btag", variable.Data()));
	TH1F *hQ = (TH1F*)infData->Get(TString::Format("hWt_%s_0btag", variable.Data()));

	//open the file to get the Ryield
	TFile *infRyield = TFile::Open(TString::Format("%s/TransferFactor.root",year.Data()));
	TH1F *hRyield = (TH1F*)infRyield->Get("ClosureTest_TransferFactor");
	float Ryield = hRyield->GetBinContent();
	float Ryield_error = hRyield->GetBinError();

	cout<<hRyield->GetBinContent(1)<<endl;
	cout<<hRyield->GetBinContent(2)<<endl;
	cout<<hRyield->GetBinContent(3)<<endl;

	//open the file to get the Nbkg
	float NQCD = Nbkg2Constants[TString::Format("Nbkg%s",year.Data())];
	float NQCD_error = Nbkg2ConstantsErrors[TString::Format("Nbkg%s_error",year.Data())];

	//Subdominant bkgs files
	TFile *infSub = TFile::Open(TString::Format("%s/Histo_SubdominantBkgs_100.root", year.Data()));
	TH1F *hSub = (TH1F*)infSub->Get(TString::Format("hWt_%s_2btag_expYield", variable.Data()));
	//here I will import correction factors for QCD if needed...

	//now i need the rebin function for aaaaaall my hists so that I am conistent
	//I will include binning in the TemplateConstants.h
	
	std::vector< std::vector <Float_t> > const BND = {{1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, 4000, 5000}, //mjj
											 {0,60,150,300,450,600,750,950,1100,1300}, //ptjj
											 {-2.4,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.4}, //yjj
											 {400,450,500,570,650,750,850,950,1100,1300,1500}, //jetPt0	
											 {400,450,500,570,650,750,850,950,1100,1300,1500}, //jetPt1
											 {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}, //jetY0
											 {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}}; //jetY1
	int nBins[BND.size()];
	for (int i = 0; i<BND.size(); i++) nBins[i] = BND[i].size()-1;

	//get the respected integer from the mapping in the TemplateConstants.h
	int selVar = variableConstant[variable.Data()];
	cout<<selVar<<endl;
	cout<<NQCD<<endl;
	

	//use this template for the initialization of the BND source given as input for the rebinned histos
	float tempBND[nBins[selVar]+1];
	std::copy(BND[selVar].begin(), BND[selVar].end(), tempBND);	

	//rebin all histograms with the same method
	TH1F *hD_rebinned, *hQ_rebinned, *hSub_rebinned;
	hD_rebinned = getRebinned(hD, tempBND, nBins[selVar]);
	hQ_rebinned = getRebinned(hQ, tempBND, nBins[selVar]);
	hSub_rebinned = getRebinned(hSub, tempBND, nBins[selVar]);

	TH1F *hSignal = (TH1F*)hD_rebinned->Clone(TString::Format("hSignal_%s",variable.Data()));

	//work on the elements for the QCD so that we have the right Q(x)
	//Ryield * Nbkg  * Q(x)
	cout<<"Ryield: "<<Ryield<<endl;
	cout<<variable<<endl;
	
	float SF[hQ_rebinned->GetNbinsX()];
	//QCD correction factor for shape
	if(year.EqualTo("2017") || year.EqualTo("2018"))
	{
		if(variable.EqualTo("jetPt0") || variable.EqualTo("jetPt1") || variable.EqualTo("mJJ") || (year.EqualTo("2017") && variable.EqualTo("ptJJ"))) 		
		{
			TFile *fitFile =  TFile::Open(TString::Format("../QCD_ClosureTests_All/QCD_Closure_%s/FitOutput.root",year.Data()));
	 	 	TF1 *fitResult = (TF1*)fitFile->Get(TString::Format("func_%s",variable.Data()));
	  		//scale now all over the bins
	  		//cout<<"beta [0]: "<<fitResult->GetParameter("beta")<<" ± "<<fitResult->GetParError(0)<<endl;
	  		//cout<<"delta [1]: "<<fitResult->GetParameter("delta")<<" ± "<<fitResult->GetParError(1)<<endl;
	  		//cout<<"alpha [2]: "<<fitResult->GetParameter("alpha")<<" ± "<<fitResult->GetParError(2)<<endl;
	  		for(int i=0; i<hQ_rebinned->GetNbinsX(); i++)
	  		{
	  			float chi = hQ_rebinned->GetBinCenter(i+1);
	    		SF[i] = fitResult->Eval(chi);
	    	}
    	}
    	else for(int i=0; i<hQ_rebinned->GetNbinsX(); i++) SF[i] = 1;
  	}
  	else
  	{
  		for(int i=0; i<hQ_rebinned->GetNbinsX(); i++) SF[i] = 1;
	}


	hQ_rebinned->Scale(1./hQ_rebinned->Integral());
	for(int i =0; i<hQ_rebinned->GetNbinsX(); i++)
	{	
		//cout<<SF[i]<<endl;
		float oldContent = hQ_rebinned->GetBinContent(i+1);
		float oldError = hQ_rebinned->GetBinError(i+1);
		float newContent = Ryield * NQCD * oldContent * SF[i];
		//cout<<Ryield * NQCD * oldContent * SF[i]<<endl;
		//cout<<NQCD2_reduced[year.Data()] * oldContent *SF[i]<<endl;
		//float newContent = NQCD2_reduced[year.Data()] * oldContent *SF[i];
		float newError   = TMath::Sqrt(TMath::Power(NQCD*oldContent*Ryield_error,2) + TMath::Power(NQCD*oldError*Ryield,2)+
										TMath::Power(NQCD_error*oldContent*Ryield,2));
		hQ_rebinned->SetBinContent(i+1, newContent);
		hQ_rebinned->SetBinError(i+1, newError);
		//now setThe content for the hSignal
	}
 	hSignal->Add(hQ_rebinned,-1);
 	hSignal->Add(hSub_rebinned,-1);

 	//for reviewing get the MC signal
 	TFile *infSignalMC = TFile::Open(TString::Format("%s/Histo_TT_Mtt-700toInf_100.root", year.Data()));
 	TH1F *hSMC = (TH1F*)infSignalMC->Get(TString::Format("hWt_%s_2btag_expYield", variable.Data()));

 	//cout<<"Entries: "<<hSMC->GetEntries()<<endl;
 	//cout<<"Integral: "<<hSMC->Integral()<<endl;
 	cout<<"-------"<<endl;
 	hSignal->SetLineColor(kBlack);
 	hSMC->SetLineColor(kBlue);
 	hSignal->SetMarkerStyle(20);
 	hSignal->SetMarkerColor(kBlack);

 	TLegend *leg = new TLegend(0.65,0.7,0.9,0.9);
 	leg->AddEntry(hSignal, "Data", "lpe");
 	leg->AddEntry(hSMC, "MC", "lp");

 	TCanvas *can = new TCanvas(TString::Format("can_%s",variable.Data()),TString::Format("can_%s",variable.Data()) , 800,600);
 	hSignal->Scale(1,"width");
 	hSMC->Scale(1, "width");
 	hSMC->GetYaxis()->SetTitle("#frac{d#sigma}{d#chi} [pb]");
 	hSMC->SetTitle(TString::Format("Data vs MC %s for %s ",year.Data(), variable.Data()));
 		if(!variable.EqualTo("yJJ")) gPad->SetLogy();

 	hSMC->Draw("e");
 	hSignal->Draw("same e0");
 	leg->Draw();

 	TString path;
 	if(free_eb) path = TString::Format("%s/FiducialMeasurement/free_eb/fiducial_%s.pdf",year.Data(),variable.Data());
 	else path = TString::Format("%s/FiducialMeasurement/fixed_eb/fiducial_%s.pdf",year.Data(),variable.Data());

	//can->Print(path,"pdf");

}