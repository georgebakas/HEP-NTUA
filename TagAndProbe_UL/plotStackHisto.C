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

void plotStackHisto_Variable(TString year, TFile *infData, TFile *infTT, TFile *infQCD, TFile *infSub, TString variable, TString regions);

void plotStackHisto(TString year)
{

  //get the files from the directory
  //data file
  TFile *infData = TFile::Open(TString::Format("../MassFit/%s/TagAndProbeHisto_Data.root", year.Data()));
  //tt nominal file:
  TFile *infTT = TFile::Open(TString::Format("%s/Nominal/combined/TagAndProbeHisto_1000_TT_Nominal.root", year.Data()));
  //qcd mc file
  TFile *infQCD = TFile::Open(TString::Format("../MassFit/%s/TagAndProbeHisto_QCD_HT300toInf.root",year.Data()));
  //subdominant file:
  TFile *infSub = TFile::Open(TString::Format("../MassFit/%s/TagAndProbeHisto_SubdominantBkgs.root",year.Data()));
  const int NVAR =8;
  TString regions[2] = {"hSRBTightAndSR_", "hSRBTightAndProbe_"};
  TString varReco[NVAR]   = {"mJJ", "ptJJ", "yJJ", "jetPt0", "jetPt1", "jetY0", "jetY1", "mTop_Leading"};
  for(int ivar = 0; ivar< NVAR; ivar++)
  {
    for(int i = 0; i<2; i++)
    {
      plotStackHisto_Variable(year, infData, infTT, infQCD, infSub, varReco[ivar],regions[i]);
    }
  }
}


void plotStackHisto_Variable(TString year, TFile *infData, TFile *infTT, TFile *infQCD, TFile *infSub, TString variable, TString regions)
{
  initFilesMapping();
  //now get the histograms
  TH1F *hData, *hTT, *hQCD, *hSub;
  hData = (TH1F*)infData->Get(TString::Format("%s%s_expYield", regions.Data(),variable.Data()));
  hTT = (TH1F*)infTT->Get(TString::Format("%s%s_expYield", regions.Data(),variable.Data()));
  //if use data, uncomment
  hQCD = (TH1F*)infQCD->Get(TString::Format("%s%s_expYield", regions.Data(),variable.Data()));

  hSub = (TH1F*)infSub->Get(TString::Format("%s%s_expYield", regions.Data(),variable.Data()));

  //scale ttbar with signal strength
  if(regions.EqualTo("hSRBTightAndSR_")) hTT->Scale(ttbarSigStrength_TagNSR[year]);
  else hTT->Scale(ttbarSigStrength_TagNProbe[year]);

  // get the NQCD 
  TFile *masFitResultsFile = TFile::Open(TString::Format("%s/Nominal/MassFitResults_%sSignalTemplates_.root",
                                                        year.Data(), 
                                                        regions.Data()));
  RooWorkspace *w = (RooWorkspace*)masFitResultsFile->Get("w");
  RooRealVar *value = (RooRealVar *)w->var("nFitQCD_2b");

  float val = value->getValV();
  float error = value->getError();
  masFitResultsFile->Close();
  // scale qcd with data driven method
  hQCD->Scale(val/hQCD->Integral());
  //scale qcd with Data
  //we use a k-factor
  TH1F *hQCD_tempFromData = (TH1F*)hData->Clone("hQCD_tempFromData");
  // hQCD_tempFromData scale this with the NQCD
  hQCD_tempFromData->Add(hTT,-1);
  hQCD_tempFromData->Add(hSub,-1);

  float qcdScaleFactor = hQCD_tempFromData->Integral()/hQCD->Integral();
  cout<<"qcdScaleFactor: "<<qcdScaleFactor<<endl;
  if(qcdScaleFactor > 0) hQCD->Scale(qcdScaleFactor);
  // cout<<"--------"<<endl;
  //cout<<regions<<endl;
  

  if(variable.Contains("mTop"))
  {
    //cout<<"hQCD_tempFromData->Integral(): "<<hQCD_tempFromData->Integral()<<endl;
    // cout<<"hData->Integral(): "<<hData->Integral()<<endl;
    // cout<<"hTT->Integral(): "<<hTT->Integral()<<endl;
    // cout<<"hSub->Integral(): "<<hSub->Integral()<<endl;
    // cout<<"hQCD->Integral(): "<<hQCD->Integral()<<endl;
    hData->Rebin(4);
    hTT->Rebin(4);
    hQCD->Rebin(4);
    hSub->Rebin(4);
  }

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


  THStack *hs = new THStack("Data vs MC", "Data vs MC;TagAndProbe;Number of Events");
  hs->Add(hSub);
  hs->Add(hQCD);
  hs->Add(hTT);

  TCanvas *can = new TCanvas(TString::Format("can_%s_%s",variable.Data(),regions.Data()), TString::Format("can_%s_%s",variable.Data(),regions.Data()), 800, 600);
  TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
  can->cd();
  TPad *closure_pad2 = new TPad(TString::Format("cp2_%s_%s",variable.Data(),regions.Data()),TString::Format("cp2_%s_%s",variable.Data(),regions.Data()),0.,0.,1.,0.3);
  closure_pad2->Draw();
  closure_pad2->SetTopMargin(0.05);
  closure_pad2->SetBottomMargin(0.25);
  closure_pad2->SetGrid();

  TPad *closure_pad1 = new TPad(TString::Format("cp1_%s_%s",variable.Data(),regions.Data()),TString::Format("cp2_%s_%s",variable.Data(),regions.Data()),0.,0.3,1.,1.);
  closure_pad1->Draw();
  closure_pad1->SetBottomMargin(0.01);
  closure_pad1->cd();

  hData->GetYaxis()->SetTitleSize(20);
  hData->GetYaxis()->SetTitleFont(43);
  hData->GetYaxis()->SetTitleOffset(1.4);
  hData->GetYaxis()->SetRangeUser(0, hData->GetMaximum() * 2);

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
  can->Print(TString::Format("%s/Nominal/plots/stacks/TagAndProbe_%s%s.pdf",year.Data(),  regions.Data(), variable.Data()),"pdf");

}
