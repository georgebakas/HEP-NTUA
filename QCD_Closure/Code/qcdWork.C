#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"

using std::cin;
using std::cout;
using std::endl;

std::vector<TString> listOfFiles;
std::vector<float> XSEC;
std::vector<TString> histoNames;
float LUMI = 35900;
TString eosPath;

void initFileNames()
{
  eosPath = "/eos/cms/store/user/ipapakri/ttbar/MC/";
  listOfFiles.push_back("Signal/TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8_legacy2016_deepAK8.root");
  listOfFiles.push_back("Bkg/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
  listOfFiles.push_back("Bkg/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
  listOfFiles.push_back("Bkg/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
  listOfFiles.push_back("Bkg/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
  listOfFiles.push_back("Bkg/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
  listOfFiles.push_back("Bkg/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
}

void initXsections()
{
  XSEC.push_back(16.74); //this is the TT xsec
  XSEC.push_back(3.67e+5);
  XSEC.push_back(2.94e+4);
  XSEC.push_back(6.524e+03);
  XSEC.push_back(1.064e+03);
  XSEC.push_back(121.5);
  XSEC.push_back(2.542e+01);
}


void initGlobals()
{
  initFileNames();
  initXsections();
}
void qcdWork(float selMvaCut = 0.1)
{
  initGlobals();
  std::vector<float> weights(0);
  for(int f=0; f<listOfFiles.size(); f++)
  //[0] is the mtt and from [1] up to [6] its the Bkg
  //the histograms from QCD are already scaled with LUMI
  {
	TFile *file =TFile::Open(eosPath+listOfFiles[f]);
	float norm = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
    float weight = XSEC[f]/norm;
    weights.push_back(weight);  
  }
  
  
  
  //files to get histograms
  TFile *infQCD = TFile::Open(TString::Format("Closure_QCDBkg_Chi_%0.1f.root", selMvaCut));
  TFile *infTT  = TFile::Open(TString::Format("Output_TT_QCD_Reco_Chi_%0.1f.root", selMvaCut));
  
  //TH1F used for the TRatioPlot
  TH1F *hChi_QCD_CR = (TH1F*)infQCD->Get("hChi_QCD_CR");
  TH1F *hCos_QCD_CR = (TH1F*)infQCD->Get("hCos_QCD_CR");
  
  TH1F *hChi_QCD_SR = (TH1F*)infQCD->Get("hChi_QCD_SR");
  TH1F *hCos_QCD_SR = (TH1F*)infQCD->Get("hCos_QCD_SR");  
  
  //TH1F used for the Top contamination differences
  TH1F *hChiCR_TT = (TH1F*)infTT->Get("hChiCR_TT");
  TH1F *hCosCR_TT = (TH1F*)infTT->Get("hCosCR_TT"); 
  hChiCR_TT->Scale(weights[0]*LUMI);
  hCosCR_TT->Scale(weights[0]*LUMI);
  hCosCR_TT->GetXaxis()->SetTitle("|cos(#theta^{*})|");
  hChiCR_TT->SetLineColor(kBlue);
  hCosCR_TT->SetLineColor(kBlue);
  
  TH1F *hChi_QCD_CR_Clone=(TH1F*)hChi_QCD_CR->Clone("hChi_QCD_CR_Clone");
  TH1F *hCos_QCD_CR_Clone=(TH1F*)hCos_QCD_CR->Clone("hCos_QCD_CR_Clone");
  hCos_QCD_CR_Clone->GetXaxis()->SetTitle("|cos(#theta^{*})|");
  
  TCanvas *canDifferences = new TCanvas("can chi", "can chi", 700, 600);
  TLegend *leg = new TLegend(0.5,0.6,0.7,0.8);
  leg->AddEntry(hChiCR_TT, "TT in CR", "l");
  leg->AddEntry(hChi_QCD_CR_Clone, "QCD in CR", "l");
  hChi_QCD_CR_Clone->SetLineColor(kRed);
  hChi_QCD_CR_Clone->Draw();
  hChiCR_TT->Draw("same");
  leg->Draw();
   
  TCanvas *canDifferencesCos = new TCanvas("can cos", "can cos", 700, 600);
  TLegend *legCos= new TLegend(0.5,0.6,0.7,0.8);
  legCos->AddEntry(hCosCR_TT, "TT in CR", "l");
  legCos->AddEntry(hCos_QCD_CR_Clone, "QCD in CR", "l");
  hCos_QCD_CR_Clone->SetLineColor(kRed);
  hCos_QCD_CR_Clone->Draw();
  hCosCR_TT->Draw("same");
  legCos->Draw();


  auto c0 = new TCanvas("#chi TT contamination", "#chi TT contamination", 700,600);
  auto rp_chiContamination = new TRatioPlot(hChiCR_TT, hChi_QCD_CR_Clone);
  c0->SetTicks(0,1);
  rp_chiContamination->Draw();
  c0->Update();

  auto c1 = new TCanvas("cos TT contamination", "cos TT contamination", 700,600);
  auto rp_cosContamination = new TRatioPlot(hCosCR_TT, hCos_QCD_CR_Clone);
  c1->SetTicks(0,1);
  rp_cosContamination->Draw();
  c1->Update();    


  auto c3 = new TCanvas("#chi", "#chi", 700,600);
  TLegend *leg_chi = new TLegend(0.6,0.7,0.8,0.9);
  hChi_QCD_SR->Scale(1./hChi_QCD_SR->Integral());
  hChi_QCD_CR->Scale(1./hChi_QCD_CR->Integral());
  leg_chi->AddEntry(hChi_QCD_SR, "#chi SR 2 btag", "l"); 
  leg_chi->AddEntry(hChi_QCD_CR, "#chi CR 0 btag", "l"); 
  auto rp_chi = new TRatioPlot(hChi_QCD_SR, hChi_QCD_CR);
  c3->SetTicks(0,1);
  rp_chi->Draw();
  leg_chi->Draw();
  c3->Update();  
  

  auto c2 = new TCanvas("Cos(#theta)", "Cos(#theta)", 700,600);
  TLegend *leg_cos = new TLegend(0.6,0.7,0.8,0.9);
  hCos_QCD_SR->Scale(1./hCos_QCD_SR->Integral());
  hCos_QCD_CR->Scale(1./hCos_QCD_CR->Integral());
  leg_cos->AddEntry(hCos_QCD_SR, "cos(#theta^{*}) SR 2 btag", "l"); 
  leg_cos->AddEntry(hCos_QCD_CR, "cos(#theta^{*}) CR 0 btag", "l"); 
  auto rp_cos = new TRatioPlot(hCos_QCD_SR, hCos_QCD_CR);
  c2->SetTicks(0,1);
  rp_cos->Draw();
  c2->Update();
  leg_cos->Draw();
}
