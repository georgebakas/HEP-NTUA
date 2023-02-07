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
#include "TemplateConstants.h"
int mass, width;
void plotStackHisto_Variable(TString year, TFile *infData, TFile *infTT, TFile *infQCD, TFile *infSub, TString variable, int btagFlag, bool topTaggerSF);

void plotStackHisto_TopTaggerErrors(TString year, int mJJCut = 1000, TString region = "SR", bool topTaggerSF = true)
{
  initFilesMapping();
  //setTDRStyle();
  //get the files from the directory
  //data file
  TFile *infData = TFile::Open(TString::Format("../MassFit/%s/Histo_Data_%s_100_reduced_UnequalBinning.root",year.Data(), year.Data()));
  //tt nominal file: or should I take it from signal extractions??
  //ttnominal because we believe that the signal strength will not be modified with a further mass cut
  //this is nominal
  TFile *infTT = TFile::Open(TString::Format("%s/Nominal/combined/HistoReduced_1000_TT.root",year.Data()));
  //TFile *infTT = TFile::Open(TString::Format("%s/FiducialMeasurement_2TeV/EqualBinning/SignalHistograms_%s.root",year.Data(),"chi"));

  //qcd mc file
  //thake qcd from mc and scale it accordingly
  TFile *infQCD = TFile::Open(TString::Format("../MassFit/%s/Histo_QCD_HT300toInf_100_reduced_UnequalBinning.root",year.Data()));

  //subdominant file:
  TFile *infSub = TFile::Open(TString::Format("../MassFit/%s/Histo_SubdominantBkgs_100_reduced_UnequalBinning.root",year.Data()));

  //define btag region:
  int btagFlag = 2;
  if (region.EqualTo("CR")) 
    btagFlag = 0;

  const int NVAR =12;
  TString varReco[NVAR] = {"jetPt0", "jetPt1", "jetY0", "jetY1", "ptJJ", "yJJ", "chi","cosTheta_0", "cosTheta_1", 
                          "mJJ", "mTop", "jetMassSoftDrop"};

  for(int ivar = 0; ivar< NVAR; ivar++)
  {
    plotStackHisto_Variable(year, infData, infTT, infQCD, infSub, varReco[ivar], btagFlag, topTaggerSF);
  }
}


void plotStackHisto_Variable(TString year, TFile *infData, TFile *infTT, TFile *infQCD, TFile *infSub, TString variable, int btagFlag, bool topTaggerSF)
{
  //initFilesMapping();
  //now get the histograms
  TH1F *hData, *hTT, *hQCD, *hSub;

  hData = (TH1F*)infData->Get(TString::Format("hWt_%s_%dbtag_expYield", variable.Data(), btagFlag));

  TString ttvar;
  if (variable.EqualTo("cosTheta_0"))
    ttvar = "cosThjetEta0";
  else if (variable.EqualTo("cosTheta_1"))
    ttvar = "cosThjetEta1";
  else if (variable.EqualTo("mTop"))
    ttvar = "mTop_Leading";
  else if (variable.EqualTo("jetMassSoftDrop"))
    ttvar = "mTop_Subleading";
  else 
    ttvar = variable;
  
  hTT = (TH1F*)infTT->Get(TString::Format("hWt_%s_%dbtag", ttvar.Data(), btagFlag));
  //hTT = (TH1F*)infTT->Get(TString::Format("hSignal_%s", variable.Data())); //to use if SignalExtraction
  //if use data, uncomment
  hQCD = (TH1F*)infQCD->Get(TString::Format("hWt_%s_%dbtag_expYield", variable.Data(), btagFlag));

  hSub = (TH1F*)infSub->Get(TString::Format("hWt_%s_%dbtag_expYield", variable.Data(), btagFlag));


  //scale ttbar with signal strength when using MC
  hTT->Scale(ttbarSigStrength[year]);

  //scale qcd with Data
  //we use a k-factor
  TH1F *hQCD_tempFromData = (TH1F*)hData->Clone("hQCD_tempFromData");
  hQCD_tempFromData->Add(hTT, -1);
  hQCD_tempFromData->Add(hSub, -1);

  // this is needed only if I want to plot the stuff related with the scale factors
  if (topTaggerSF)
  {
    float top_tagger_sf = (topTaggerSF_data[year.Data()]/topTaggerSF_sim[year.Data()]);
    cout<<top_tagger_sf<<endl;
    cout<< topTaggerSF_data[year.Data()]<<endl;
    cout<< topTaggerSF_sim[year.Data()]<<endl;
    cout<<"----------------------"<<endl;
    float top_tagger_sf_error = TMath::Sqrt( 
                                TMath::Power(topTaggerSF_data_error[year.Data()]/topTaggerSF_sim[year.Data()],2) +
                                TMath::Power(topTaggerSF_data[year.Data()]*topTaggerSF_sim_error[year.Data()]/TMath::Power(topTaggerSF_sim[year.Data()],2),2));
    //now apply these errors accordingly
    for(int ibin=1; ibin<=hTT->GetNbinsX(); ibin++)
    { 
        float new_value = hTT->GetBinContent(ibin) * top_tagger_sf;
        float new_error = TMath::Sqrt(TMath::Power(hTT->GetBinError(ibin) ,2)+ 
                                      TMath::Power(top_tagger_sf_error ,2));
        hTT->SetBinError(ibin, new_error);
        hTT->SetBinContent(ibin, new_value);
    }
  }

  float qcdScaleFactor = hQCD_tempFromData->Integral()/hQCD->Integral();
  cout<<"qcdScaleFactor: "<<qcdScaleFactor<<endl;
  hQCD->Scale(qcdScaleFactor);
  
  TH1F *hSum = (TH1F*)hSub->Clone("sum");
  hSum->Add(hQCD);
  hSum->Add(hTT);
  hSum->SetLineColor(kMagenta);
  hSum->SetMarkerStyle(21);
  hSum->SetMarkerColor(kMagenta);

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
  

  THStack *hs = new THStack("", ";Top Angular Dists;Number of Events");
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
  //if(variable.EqualTo("mJJ")) closure_pad1->SetLogy();
  closure_pad1->cd();

  hData->GetYaxis()->SetTitleSize(20);
  hData->GetYaxis()->SetTitleFont(43);
  hData->GetYaxis()->SetTitleOffset(1.4);
  hData->GetYaxis()->SetRangeUser(0.01, hData->GetMaximum() * 1.2);



  leg->AddEntry(hData, "Data", "lep");
  leg->AddEntry(hTT, "TTbar", "f");
  leg->AddEntry(hQCD, "QCD", "f");
  leg->AddEntry(hSub, "Subdominant", "f");


  hs->Draw("hist");
  hData->Draw("hist same E");
  hSum->Draw("SAME");
  hs->GetYaxis()->SetTitle("Number of Events");
  
  if (variable.EqualTo("chi")) hs->SetMaximum(hs->GetMaximum());
  else hs->SetMaximum(hs->GetMaximum()* 2);
  hs->SetMinimum(0.001);
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
  if (variable.EqualTo("chi"))
    hNum->GetXaxis()->SetTitle("#chi");
  else if (variable.Contains("cosTheta"))
    hNum->GetXaxis()->SetTitle("|cos(#theta^{*})|");
  else if (variable.Contains("jetY0"))
    hNum->GetXaxis()->SetTitle("Leading Jet abs Y");
  else if (variable.Contains("jetY0"))
    hNum->GetXaxis()->SetTitle("Leading Jet abs Y");
  else if (variable.Contains("jetY1"))
    hNum->GetXaxis()->SetTitle("Second Leading Jet abs Y");
  else if (variable.Contains("jetPt0"))
    hNum->GetXaxis()->SetTitle("Leading Jet p_{T}");
  else if (variable.Contains("jetPt1"))
    hNum->GetXaxis()->SetTitle("Second Leading Jet p_{T}");
  else if (variable.Contains("yJJ"))
    hNum->GetXaxis()->SetTitle(variable);
  else 
    hNum->GetXaxis()->SetTitle(TString::Format("%s (GeV)",variable.Data()));
  hNum->GetYaxis()->SetTitleSize(20);
  hNum->GetYaxis()->SetTitleFont(43);
  hNum->GetYaxis()->SetTitleOffset(1.3);
  hNum->GetYaxis()->SetLabelFont(43);
  hNum->GetYaxis()->SetLabelSize(15);
  hNum->GetXaxis()->SetTitleSize(0.1);
  hNum->GetXaxis()->SetLabelFont(43);
  hNum->GetXaxis()->SetLabelSize(13);
  hNum->Draw();

  TString lumi_str = TString::Format("%0.1f", luminosity["luminosity"+year]/1000);
  lumi_13TeV = lumi_str+" fb^{-1}";

  float extraTextFactor = 0.14;
  int iPeriod = 13;
  int iPos = 0;
  writeExtraText=true;
  CMS_lumi(closure_pad1, year, iPos);

  if (topTaggerSF)
  {
    can->Print(TString::Format("%s/StackPlots_TopTagSF/DatavsMC_%s_%dbTag_E2.pdf",year.Data(), variable.Data(), btagFlag),"pdf");
  }
  else
    can->Print(TString::Format("%s/StackPlots/DatavsMC_%s_%dbTag_E2.pdf",year.Data(), variable.Data(), btagFlag),"pdf");

}
