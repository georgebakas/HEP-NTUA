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
#include "TemplateConstants.h"
bool useDataQCD;

void plotStackHisto_Variable(TString year, TFile *infData, TFile *infTT, TFile *infQCD, TFile *infSub, TString variable, TString leadingStr);

void plotStackHisto(TString year, bool useDataQCD1)
{
  useDataQCD = useDataQCD1;
  //get the files from the directory
  //data file
  TFile *infData = TFile::Open(TString::Format("%s/TopTaggerHisto_Data_%s_100.root", year.Data(), year.Data()));
  //tt nominal file:
  TFile *infTT = TFile::Open(TString::Format("%s/TopTaggerHisto_TT_NominalMC_100.root", year.Data()));
  //qcd mc file
  TFile *infQCD = TFile::Open(TString::Format("%s/TopTaggerHisto_QCD_HT300toInf_100.root", year.Data()));
  //subdominant file:
  TFile *infSub = TFile::Open(TString::Format("%s/TopTaggerHisto_SubdominantBkgs_100.root", year.Data()));

  const int NVAR =11;
  TString leadStr[] = {"leading", "subleading"};
  TString varReco[NVAR]   = {"topTagger","jetTau3", "jetTau2", "jetTau1","jetMassSub0","jetMassSub1",
                              "ecfB1N2", "ecfB1N3","ecfB2N2", "ecfB2N3", "JetPtOverSumPt"};
  for(int ivar = 0; ivar< NVAR; ivar++)
  {
    for(int ilead = 0; ilead<2; ilead++)
    {
      plotStackHisto_Variable(year, infData, infTT, infQCD, infSub, varReco[ivar], leadStr[ilead]);
    }
  }
}


void plotStackHisto_Variable(TString year, TFile *infData, TFile *infTT, TFile *infQCD, TFile *infSub, TString variable, TString leadingStr)
{
  initFilesMapping(false);
  //now get the histograms
  TH1F *hData, *hTT, *hQCD, *hSub;

  hData = (TH1F*)infData->Get(TString::Format("hWt_%s_2btag_expYield_%s", variable.Data(),leadingStr.Data()));
  hTT = (TH1F*)infTT->Get(TString::Format("hWt_%s_2btag_expYield_%s", variable.Data(),leadingStr.Data()));
  //if use data, uncomment
  TString qcdStr = "qcdMC";
  if(useDataQCD)
  {
    qcdStr = "qcdData";
    hQCD = (TH1F*)infData->Get(TString::Format("hWt_%s_0btag_expYield_%s", variable.Data(),leadingStr.Data()));
    hQCD->Scale(1./hQCD->Integral());
    hQCD->Scale(Nbkg2Constants[year]);
  }
  else
    hQCD = (TH1F*)infQCD->Get(TString::Format("hWt_%s_2btag_expYield_%s", variable.Data(),leadingStr.Data()));

  hSub = (TH1F*)infSub->Get(TString::Format("hWt_%s_2btag_expYield_%s", variable.Data(),leadingStr.Data()));

  //scale ttbar with signal strength
  hTT->Scale(signalStrenth[year]);

  //scale qcd with Data
  //we use a k-factor
  TH1F *hQCD_tempFromData = (TH1F*)hData->Clone("hQCD_tempFromData");
  hQCD_tempFromData->Add(hTT, -1);
  hQCD_tempFromData->Add(hSub, -1);

  float qcdScaleFactor = hQCD_tempFromData->Integral()/hQCD->Integral();
  cout<<"qcdScaleFactor: "<<qcdScaleFactor<<endl;
  hQCD->Scale(qcdScaleFactor);

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
  hSub->Rebin(2);
  hQCD->Rebin(2);
  hTT->Rebin(2);
  THStack *hs = new THStack("Data vs MC", "Data vs MC;TopTagger Output;Number of Events");
  hs->Add(hSub);
  hs->Add(hQCD);
  hs->Add(hTT);

  TCanvas *can = new TCanvas(TString::Format("can_%s_%s",variable.Data(),leadingStr.Data()), TString::Format("can_%s_%s",variable.Data(),leadingStr.Data()), 800, 600);
  TLegend *leg;
  if(!variable.Contains("topTagger") && !variable.EqualTo("ecfB1N2"))
    leg = new TLegend(0.7,0.7,0.9,0.9);
  else
    leg = new TLegend(0.10,0.7,0.25,0.9);
  can->cd();
  TPad *closure_pad2 = new TPad(TString::Format("cp2_%s_%s",variable.Data(),leadingStr.Data()),TString::Format("cp2_%s_%s",variable.Data(),leadingStr.Data()),0.,0.,1.,0.3);
  closure_pad2->Draw();
  closure_pad2->SetTopMargin(0.05);
  closure_pad2->SetBottomMargin(0.25);
  closure_pad2->SetGrid();

  TPad *closure_pad1 = new TPad(TString::Format("cp1_%s_%s",variable.Data(),leadingStr.Data()),TString::Format("cp2_%s_%s",variable.Data(),leadingStr.Data()),0.,0.3,1.,1.);
  closure_pad1->Draw();
  closure_pad1->SetBottomMargin(0.01);
  closure_pad1->cd();

  hData->GetYaxis()->SetTitleSize(20);
  hData->GetYaxis()->SetTitleFont(43);
  hData->GetYaxis()->SetTitleOffset(1.4);
  hData->GetYaxis()->SetRangeUser(0, hData->GetMaximum() * 1.2);

  leg->AddEntry(hData, "Data", "lep");
  leg->AddEntry(hTT, "TTbar", "f");
  leg->AddEntry(hQCD, "QCD", "f");
  leg->AddEntry(hSub, "Subdominant", "f");

  hs->Draw("hist");
  hData->Draw("same E");
  hs->GetYaxis()->SetTitle("Number of Events");
  leg->Draw();


  closure_pad2->cd();
  TH1F *hDenom = (TH1F*)hQCD->Clone("hDenom");
  hDenom->Add(hSub);
  hDenom->Add(hTT);
  TH1F *hNum = (TH1F*)hData->Clone("hNum");
  hNum->Divide(hDenom);
  hNum->SetTitle("");
  hNum->GetYaxis()->SetRangeUser(0,3);
  hNum->GetYaxis()->SetTitle("#frac{Data}{MC}");
  hNum->GetXaxis()->SetTitle(variable);
  hNum->GetYaxis()->SetTitleSize(20);
  hNum->GetYaxis()->SetTitleFont(43);
  hNum->GetYaxis()->SetTitleOffset(1.3);
  hNum->GetYaxis()->SetLabelFont(43);
  hNum->GetYaxis()->SetLabelSize(15);
  hNum->GetXaxis()->SetTitleSize(0.1);
  hNum->GetXaxis()->SetLabelFont(43);
  hNum->GetXaxis()->SetLabelSize(13);

  hNum->Draw();
  can->Print(TString::Format("%s/plots/%s/TopTaggerDatavsMC_%s_%s.pdf",year.Data(), qcdStr.Data(), variable.Data(), leadingStr.Data()),"pdf");

}
