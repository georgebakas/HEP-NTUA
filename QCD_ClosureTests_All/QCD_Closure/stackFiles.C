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

std::vector<TString> histoNames;

void initHistoNames()
{
	histoNames.push_back("QCD_300_500");
	histoNames.push_back("QCD_500_700");
	histoNames.push_back("QCD_700_1000");
	histoNames.push_back("QCD_1000_1500");
	histoNames.push_back("QCD_1500_2000");
	histoNames.push_back("QCD_2000_Inf");
}


void stackFiles()
{
  initHistoNames();
  TFile *inf = TFile::Open("htBkgFile.root");
  THStack *hStack = new THStack("ht hStack", "ht hStack");
  //hStack->GetXaxis()->SetTitle("hT (GeV)");
  //hStack->GetYaxis()->SetTitle("Expected Yield");
  std::vector<Color_t> colors= {kBlue, kBlack, kRed, kMagenta,kGreen, kCyan};

  TLegend *leg = new TLegend(0.5, 0.6, 0.7, 0.8);
  TH1F *h_hT[histoNames.size()];
  for(int f=0; f<histoNames.size(); f++)
  {
  	h_hT[f] = (TH1F*)inf->Get(histoNames[f].Data());
  	h_hT[f]->SetLineColor(colors[f]);
  	h_hT[f]->SetFillColor(colors[f]);
    h_hT[f]->SetOption("hist");
  	h_hT[f]->GetXaxis()->SetTitle("hT (GeV)");
  	h_hT[f]->GetYaxis()->SetTitle("Expected Yield");
    h_hT[f]->Draw();
  	
  }

  for(int f=0; f<histoNames.size(); f++)
  {
    hStack->Add(h_hT[f]);
    leg->AddEntry(h_hT[f],histoNames[f].Data(), "f");
  }


  hStack->Draw("hist");
  leg->Draw();
}