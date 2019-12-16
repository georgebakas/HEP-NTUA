#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TUnfold.h"


#include "TemplateConstants.h"

void plotAllYears(TString variable = "jetPt0")
{
	initFilesMapping();
	TFile *inf16 = TFile::Open("2016/FiducialMeasurement/EqualBinning/ABCDMethod/free_eb/SignalHistograms_ABCDMethod_freeEb.root");
	TFile *inf17 = TFile::Open("2017/FiducialMeasurement/EqualBinning/ABCDMethod/free_eb/SignalHistograms_ABCDMethod_freeEb.root");
	TFile *inf18 = TFile::Open("2018/FiducialMeasurement/EqualBinning/ABCDMethod/free_eb/SignalHistograms_ABCDMethod_freeEb.root");



	TH1F *h16 = (TH1F*)inf16->Get(TString::Format("hSignal_%s", variable.Data()));
	TH1F *h17 = (TH1F*)inf17->Get(TString::Format("hSignal_%s", variable.Data()));
	TH1F *h18 = (TH1F*)inf18->Get(TString::Format("hSignal_%s", variable.Data()));

	TH1F *hMC16 = (TH1F*)inf16->Get(TString::Format("hSMC_%s", variable.Data()));
	TH1F *hMC17 = (TH1F*)inf17->Get(TString::Format("hSMC_%s", variable.Data()));
	TH1F *hMC18 = (TH1F*)inf18->Get(TString::Format("hSMC_%s", variable.Data()));


	h16->SetMarkerColor(kBlue);
	h17->SetMarkerColor(kRed);
	h18->SetMarkerColor(kGreen);

	h16->SetLineColor(kBlue);
	h17->SetLineColor(kRed);
	h18->SetLineColor(kGreen);

	h16->SetMarkerStyle(20);
	h17->SetMarkerStyle(21);
	h18->SetMarkerStyle(22);

	hMC16->SetMarkerColor(kBlue);
	hMC17->SetMarkerColor(kRed);
	hMC18->SetMarkerColor(kGreen);

	hMC16->SetLineColor(kBlue);
	hMC17->SetLineColor(kRed);
	hMC18->SetLineColor(kGreen);

	hMC16->SetMarkerStyle(20);
	hMC17->SetMarkerStyle(21);
	hMC18->SetMarkerStyle(22);

	h16->Scale(1/luminosity["2016"], "width");
	h17->Scale(1/luminosity["2017"], "width");
	h18->Scale(1/luminosity["2018"], "width");
	
	hMC16->Scale(1/luminosity["2016"], "width");
	hMC17->Scale(1/luminosity["2017"], "width");
	hMC18->Scale(1/luminosity["2018"], "width");

	TCanvas *canData = new TCanvas ("canData", "canData", 800,600);
	h18->SetTitle("Fid. Measurements All years combined");
	h18->GetXaxis()->SetTitle("jetPt0 (GeV)");
	h18->GetYaxis()->SetTitle("#frac{d#sigma}{d#chi} [pb]");
	TLegend *leg = new TLegend(0.65,0.7,0.9,0.9);
    leg->AddEntry(h16, "Data '16", "lpe");
    leg->AddEntry(h17, "Data '17", "lpe");
    leg->AddEntry(h18, "Data '18", "lpe");
    
	h18->Draw();
	h17->Draw("same");
	h16->Draw("same");
	gPad->SetLogy();
	leg->Draw();

	TCanvas *canMC = new TCanvas ("canMC", "canMC", 800,600);
	hMC18->SetTitle("Fid. Measurements All years combined");
	hMC18->GetXaxis()->SetTitle("jetPt0 (GeV)");
	hMC18->GetYaxis()->SetTitle("#frac{d#sigma}{d#chi} [pb]");
	TLegend *legMC = new TLegend(0.65,0.7,0.9,0.9);
    legMC->AddEntry(hMC16, "MC '16", "lpe");
    legMC->AddEntry(hMC17, "MC '17", "lpe");
    legMC->AddEntry(hMC18, "MC '18", "lpe");

	hMC18->Draw();
	hMC17->Draw("same");
	hMC16->Draw("same");
	gPad->SetLogy();
	legMC->Draw();

	canData->Print(TString::Format("allYearsPlots/FiducialMeasurements/dataFiducial_%s.pdf",variable.Data()),"pdf");
	canMC->Print(TString::Format("allYearsPlots/FiducialMeasurements/mcFiducial_%s.pdf",variable.Data()),"pdf");

}