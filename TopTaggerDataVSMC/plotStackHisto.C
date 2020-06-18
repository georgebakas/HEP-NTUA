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

void plotStackHisto(TString year = "2016")
{
  //get the files from the directory
  //data file
  TFile *infData = TFile::Open(TString::Format("%s/TopTaggerHisto_Data_%s_100.root", year.Data(), year.Data()));

  //tt nominal file:
  TFile *infTT = TFile::Open(TString::Format("%s/TopTaggerHisto_TT_NominalMC_100.root", year.Data()));
  //qcd mc file
  TFile *infQCD = TFile::Open(TString::Format("%s/TopTaggerHisto_QCD_HT300toInf_100.root", year.Data()));
  //subdominant file:
  TFile *infSub = TFile::Open(TString::Format("%s/TopTaggerHisto_SubdominantBkgs_100.root", year.Data()));


  //now get the histograms
  TH1F *hData, *hTT, *hQCD, *hSub;

  hData = (TH1F*)infData->Get("hWt_mva_2btag_expYield");
  hTT = (TH1F*)infTT->Get("hWt_mva_2btag_expYield");
  //hQCD = (TH1F*)infData->Get("hWt_mva_0btag_expYield");
  hQCD = (TH1F*)infQCD->Get("hWt_mva_2btag_expYield");
  hSub = (TH1F*)infSub->Get("hWt_mva_2btag_expYield");

  //make them pretty :D
  hTT->SetLineColor(kRed-9);
  hTT->SetMarkerColor(kRed-9);
  hTT->SetFillColor(kRed-9);

  hQCD->SetLineColor(kBlue-6);
  hQCD->SetMarkerColor(kBlue-6);
  hQCD->SetFillColor(kBlue-6);

  hSub->SetLineColor(kTeal-7);
  hSub->SetMarkerColor(kTeal-7);
  hSub->SetFillColor(kTeal-7);

  hData->SetLineColor(kBlack);
  hData->SetMarkerStyle(20);
  hData->SetMarkerColor(kBlack);

  hData->Rebin(2);
  hQCD->Rebin(2);
  hTT->Rebin(2);
  hSub->Rebin(2);

  /*hData->Scale(1/hData->Integral());
  hTT->Scale(1/hTT->Integral());
  hQCD->Scale(1/hQCD->Integral());
  hSub->Scale(1/hSub->Integral()); */

  THStack *hs = new THStack("Data vs MC", "Data vs MC;TopTagger Output;Number of Events");
  hs->Add(hQCD);
  hs->Add(hSub);
  hs->Add(hTT);

  TCanvas *can = new TCanvas("can_", "can_", 800, 600);
  TLegend *leg = new TLegend(0.65,0.7,0.9,0.9);
  can->cd();
  TPad *closure_pad2 = new TPad("cp2","cp2",0.,0.,1.,0.3);
  closure_pad2->Draw();
  closure_pad2->SetTopMargin(0.05);
  closure_pad2->SetBottomMargin(0.25);
  closure_pad2->SetGrid();

  TPad *closure_pad1 = new TPad("cp1","cp1",0.,0.3,1.,1.);
  closure_pad1->Draw();
  closure_pad1->SetBottomMargin(0.01);
  closure_pad1->cd();

  hData->GetYaxis()->SetTitleSize(20);
  hData->GetYaxis()->SetTitleFont(43);
  hData->GetYaxis()->SetTitleOffset(1.4);


  leg->AddEntry(hData, "Data", "lep");
  leg->AddEntry(hTT, "TTbar", "f");
  leg->AddEntry(hQCD, "QCD", "f");
  leg->AddEntry(hSub, "Subdominant", "f");
  leg->Draw();

  hs->Draw("hist");
  hData->Draw("same E");
  hs->GetYaxis()->SetTitle("Number of Events");

  closure_pad2->cd();
  TH1F *hDenom = (TH1F*)hQCD->Clone("hDenom");
  hDenom->Add(hSub);
  hDenom->Add(hTT);
  TH1F *hNum = (TH1F*)hData->Clone("hNum");
  hNum->Divide(hDenom);
  hNum->SetTitle("");
  hNum->GetYaxis()->SetTitle("#frac{Data}{MC}");
  hNum->GetXaxis()->SetTitle("TopTagger Output");
  hNum->GetYaxis()->SetTitleSize(20);
  hNum->GetYaxis()->SetTitleFont(43);
  hNum->GetYaxis()->SetTitleOffset(1.3);
  hNum->GetYaxis()->SetLabelFont(43);
  hNum->GetYaxis()->SetLabelSize(15);
  hNum->GetXaxis()->SetTitleSize(0.1);
  hNum->GetXaxis()->SetLabelFont(43);
  hNum->GetXaxis()->SetLabelSize(13);

  hNum->Draw();
  can->Print(TString::Format("%s/TopTaggerDatavsMC.pdf",year.Data()),"pdf");

}
