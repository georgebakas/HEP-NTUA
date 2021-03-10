#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "../CMS_plots/tdrstyle.C"
#include "../CMS_plots/CMS_lumi.C"

using std::cin;
using std::cout;
using std::endl;
#include "TemplateConstants_FillHistograms.h"
int mass, width;

void plotStackHisto_Variable(TString year, TFile *infData, TFile *infTT, TFile *infQCD, TFile *infSub, TString variable, int mJJCut);

void plotStackHisto(TString year, int mJJCut = 2000, int selMass= 2500, int selWidth=25)
{
  setTDRStyle();
  initFilesMapping();
  mass = selMass;
  width = selWidth;
  //get the files from the directory
  //data file
  TFile *infData = TFile::Open(TString::Format("%s/Histo_Data_%s_reduced_%d.root",year.Data(), year.Data(), mJJCut));
  //tt nominal file: or should I take it from signal extractions??
  //ttnominal because we believe that the signal strength will not be modified with a further mass cut
  //this is nominal
  TFile *infTT = TFile::Open(TString::Format("%s/Histo_TT_NominalMC_reduced_%d.root",year.Data(), mJJCut));
  //TFile *infTT = TFile::Open(TString::Format("%s/FiducialMeasurement_2TeV/EqualBinning/SignalHistograms_%s.root",year.Data(),"chi"));

  //qcd mc file
  //thake qcd from mc and scale it accordingly
  TFile *infQCD = TFile::Open(TString::Format("%s/Histo_QCD_HT300toInf_reduced_%d.root",year.Data(), mJJCut));

  //subdominant file:
  TFile *infSub = TFile::Open(TString::Format("%s/Histo_SubdominantBkgs_reduced_%d.root",year.Data(), mJJCut));

  const int NVAR =4;
  TString varReco[NVAR]   = {"chi","cosTheta_0", "cosTheta_1", "mJJ"};

  for(int ivar = 0; ivar< NVAR; ivar++)
  {
    plotStackHisto_Variable(year, infData, infTT, infQCD, infSub, varReco[ivar], mJJCut);
  }
}


void plotStackHisto_Variable(TString year, TFile *infData, TFile *infTT, TFile *infQCD, TFile *infSub, TString variable, int mJJCut)
{
  //initFilesMapping();
  //now get the histograms
  TH1F *hData, *hTT, *hQCD, *hSub;

  hData = (TH1F*)infData->Get(TString::Format("hWt_%s_2btag_expYield", variable.Data()));

  hTT = (TH1F*)infTT->Get(TString::Format("hWt_%s_2btag_expYield", variable.Data()));
  //hTT = (TH1F*)infTT->Get(TString::Format("hSignal_%s", variable.Data())); //to use if SignalExtraction
  //if use data, uncomment
  hQCD = (TH1F*)infQCD->Get(TString::Format("hWt_%s_2btag_expYield", variable.Data()));

  hSub = (TH1F*)infSub->Get(TString::Format("hWt_%s_2btag_expYield", variable.Data()));


  //scale ttbar with signal strength when using MC
  hTT->Scale(ttbarSigStrength[year]);

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

  THStack *hs = new THStack("Data vs MC", "Data vs MC;Top Angular Dists;Number of Events");
  hs->Add(hSub);
  hs->Add(hQCD);
  hs->Add(hTT);

  TCanvas *can = new TCanvas(TString::Format("can_%s",variable.Data()), TString::Format("can_%s",variable.Data()), 800, 600);
  TLegend *leg;
  if(variable.Contains("cos"))leg = new TLegend(0.40,0.7,0.65,0.9);
  else leg = new TLegend(0.7,0.7,0.9,0.9);
  can->cd();
  TPad *closure_pad2 = new TPad(TString::Format("cp2_%s",variable.Data()),TString::Format("cp2_%s",variable.Data()),0.,0.,1.,0.3);
  closure_pad2->Draw();
  closure_pad2->SetTopMargin(0.05);
  closure_pad2->SetBottomMargin(0.25);
  closure_pad2->SetGrid();

  TPad *closure_pad1 = new TPad(TString::Format("cp1_%s",variable.Data()),TString::Format("cp2_%s",variable.Data()),0.,0.3,1.,1.);
  closure_pad1->Draw();
  closure_pad1->SetBottomMargin(0.01);
  if(variable.EqualTo("mJJ")) closure_pad1->SetLogy();
  closure_pad1->cd();

  hData->GetYaxis()->SetTitleSize(20);
  hData->GetYaxis()->SetTitleFont(43);
  hData->GetYaxis()->SetTitleOffset(1.4);
  hData->GetYaxis()->SetRangeUser(0.01, hData->GetMaximum() * 1.2);

  //add the Zprime contribution
  TFile *infZprime;
  if(year.EqualTo("2016")) infZprime = TFile::Open(TString::Format("../Zprime/%s/HistoMassWindows_ZprimeToTT_M-%d_W-%d_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root", year.Data(), mass, width));
  else if (year.EqualTo("2017")) infZprime = TFile::Open(TString::Format("../Zprime/%s/HistoMassWindows_ZprimeToTT_M%d_W%d_TuneCP2_PSweights_13TeV-madgraph-pythiaMLM-pythia8_19UL.root", year.Data(), mass, width));  
  else infZprime = TFile::Open(TString::Format("../Zprime/%s/HistoMassWindows_ZprimeToTT_M%d_W%d_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8.root", year.Data(), mass, width));
  //TFile *infZprime = TFile::Open(TString::Format("../Zprime/%s/HistoMassWindows_ZprimeToTT_M-%d_W-%d_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root", year.Data(), mass, width));
  TH1F *hZ = (TH1F*)infZprime->Get(TString::Format("hReco_%s_%d", variable.Data(), mJJCut));
  hZ->SetLineColor(kGray);
  //hZ->SetMarkerColor(kCyan);
  //hZ->SetMarkerStyle(5);
  hZ->SetFillColor(kGray);
  hZ->SetFillStyle(3021);


  leg->AddEntry(hData, "Data", "lep");
  leg->AddEntry(hTT, "TTbar", "f");
  leg->AddEntry(hQCD, "QCD", "f");
  leg->AddEntry(hSub, "Subdominant", "f");
  leg->AddEntry(hZ, "Zprime", "lep");

  cout<<"hZ integral: "<<hZ->Integral()<<endl;

  hs->Draw("hist");
  hZ->Draw("hist same E");
  hData->Draw("hist same E");
  hs->GetYaxis()->SetTitle("Number of Events");
  hs->SetMaximum(hs->GetMaximum()* 2);
  leg->Draw();


  closure_pad2->cd();
  TH1F *hDenom = (TH1F*)hQCD->Clone("hDenom");
  hDenom->Add(hSub);
  hDenom->Add(hTT);
  TH1F *hNum = (TH1F*)hData->Clone("hNum");
  hNum->Divide(hDenom);
  hNum->SetTitle("");
  hNum->GetYaxis()->SetRangeUser(0.01,3);
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

  TString lumi_str = TString::Format("%0.1f", luminosity[year]/1000);
  lumi_13TeV = lumi_str+" fb^{-1}";

  //lumi_sqrtS = "13 TeV";
  int iPeriod = 4;
  int iPos = 0;
  //writeExtraText=true;
  CMS_lumi(closure_pad1, iPeriod, iPos);
  can->Print(TString::Format("%s/StackPlots/DatavsMC_%s_%d_M%dW%d.pdf",year.Data(), variable.Data(), mJJCut, mass, width),"pdf");

}
