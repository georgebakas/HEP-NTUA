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
We will have 16 Control Regions:
For each jet: 
1. Leading jet:
	a. 1top1b
	b. 1top0b
	c. 0top1b
	d. 0top0b
2. Sub-Leading jet:
	a. 1top1b
	b. 1top0b
	c. 0top1b
	d. 0top0b

The combination of each of the subcategories will give 16 in total

Control Regions:

Region LeadingJet SubleadingJet
A      0t0b       0t0b
B      0t0b       0t1b
C      1t0b       0t0b
D      1t0b       0t1b

E      0t0b       1t0b
F      1t0b       1t0b
G      0t1b       1t0b
H      0t1b       0t1b

I      0t1b       0t0b
J      0t0b       1t1b
K      1t0b       1t1b
L      0t1b       1t1b

M      1t1b       1t0b
N      1t1b       0t1b
O      1t1b       0t0b
Signal 1t1b       1t1b

*/



void ControlRegions(TString year = "2016")
{

	TFile *ttFileNominal = TFile::Open(TString::Format("%s/Histo_TT_NominalMC_100_reduced.root", year.Data()));
	TFile *ttFileMtt = TFile::Open(TString::Format("%s/Histo_TT_Mtt-700toInf_100_reduced.root", year.Data())); 

	TFile *infData = TFile::Open(TString::Format("%s/Histo_Data_%s_100_reduced.root", year.Data(), year.Data()));	
	TFile *infQCD = TFile::Open(TString::Format("%s/Histo_QCD_HT300toInf_100_reduced.root", year.Data()));
	TFile *infBkg = TFile::Open(TString::Format("%s/Histo_SubdominantBkgs_100_reduced.root", year.Data()));


	//get all the control regions and find the ttbar yield in every region
	const int controlRegions = 16;
	TString regions[controlRegions] = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "Signal"};

	float totalYield[controlRegions];
	float ttFraction[controlRegions];

	TH1F *hData[controlRegions], *hTTMtt[controlRegions], *hQCD[controlRegions], *hBkg[controlRegions], *hTTNominal[controlRegions];

	for(int ireg=0; ireg<controlRegions; ireg++)
	{
		hData[ireg] = (TH1F*)infData->Get(TString::Format("hWt_mTop_CR%s_expYield",regions[ireg].Data()));
		hTTMtt[ireg] = (TH1F*)ttFileMtt->Get(TString::Format("hWt_mTop_CR%s_expYield",regions[ireg].Data()));
		hTTNominal[ireg] = (TH1F*)ttFileNominal->Get(TString::Format("hWt_mTop_CR%s_expYield",regions[ireg].Data()));
		hQCD[ireg]= (TH1F*)infQCD->Get(TString::Format("hWt_mTop_CR%s_expYield",regions[ireg].Data()));
		hBkg[ireg] = (TH1F*)infBkg->Get(TString::Format("hWt_mTop_CR%s_expYield",regions[ireg].Data()));

		//totalYield[ireg] = hTTNominal[ireg]->Integral() + hQCD[ireg]->Integral() + hBkg[ireg]->Integral();
		totalYield[ireg] = hData[ireg]->Integral();
		ttFraction[ireg] = hTTNominal[ireg]->Integral() / totalYield[ireg];
		
		//cout<<totalYield[ireg]<<endl;
		cout<<"The ttbar fraction in CR_"<<regions[ireg]<<" is: "<<ttFraction[ireg]<<endl;

	}



}
