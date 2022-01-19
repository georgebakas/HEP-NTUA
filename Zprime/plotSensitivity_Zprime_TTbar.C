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

void plotStackSensitivity_Variable(TString year, TFile *infTT, TString variable, int mJJCut);

void plotSensitivity_Zprime_TTbar(TString year, int mJJCut = 2000, int selMass= 2500, int selWidth=25)
{
  setTDRStyle();
  initFilesMapping(false);
  mass = selMass;
  width = selWidth;

  //tt nominal file: or should I take it from signal extractions??
  //ttnominal because we believe that the signal strength will not be modified with a further mass cut
  //this is nominal
  TFile *infTT = TFile::Open(TString::Format("../VariationHandling/%s/Nominal/combined/HistoReduced_1000_TT.root",year.Data()));

  const int NVAR =4;
  TString varReco[NVAR]   = {"chi","cosTheta_0", "cosTheta_1", "mJJ"};

  for(int ivar = 0; ivar< NVAR; ivar++)
  {
    plotStackSensitivity_Variable(year, infTT, varReco[ivar], mJJCut);
  }
}


void plotStackSensitivity_Variable(TString year, TFile *infTT, TString variable, int mJJCut)
{
  initFilesMapping(false);
  //now get the histograms
  TH1F *hData, *hTT, *hQCD, *hSub;

  hTT = (TH1F*)infTT->Get(TString::Format("hWt_%s_2btag_expYield", variable.Data()));

  //scale ttbar with signal strength when using MC
  hTT->Scale(ttbarSigStrength[year]);


  //make them pretty :D
  hTT->SetLineColor(kRed-9);
  hTT->SetMarkerColor(kRed-9);
  hTT->SetFillColor(kRed-9);
  hTT->Scale(1./hTT->Integral());

  TCanvas *can = new TCanvas(TString::Format("can_%s",variable.Data()), TString::Format("can_%s",variable.Data()), 800, 600);
  TLegend *leg;
  if(variable.Contains("cos"))leg = new TLegend(0.40,0.7,0.65,0.9);
  else leg = new TLegend(0.7,0.7,0.9,0.9);
  can->cd();

  // Zprime contribution
  TFile *infZprime;
  if(year.EqualTo("2016_preVFP")) infZprime = TFile::Open(TString::Format("%s/HistoMassWindows_ZprimeToTT_M%d_W%d_TuneCP2_PSweights_13TeV-madgraph-pythiaMLM-pythia8_20UL.root", year.Data(), mass, width));
  else if (year.EqualTo("2017")) infZprime = TFile::Open(TString::Format("%s/HistoMassWindows_ZprimeToTT_M%d_W%d_TuneCP2_PSweights_13TeV-madgraph-pythiaMLM-pythia8_20UL.root", year.Data(), mass, width));  
  else infZprime = TFile::Open(TString::Format("%s/HistoMassWindows_ZprimeToTT_M%d_W%d_TuneCP2_PSweights_13TeV-madgraph-pythiaMLM-pythia8_20UL.root", year.Data(), mass, width));
  //TFile *infZprime = TFile::Open(TString::Format("../Zprime/%s/HistoMassWindows_ZprimeToTT_M-%d_W-%d_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root", year.Data(), mass, width));
  TH1F *hZ = (TH1F*)infZprime->Get(TString::Format("hReco_%s_%d", variable.Data(), mJJCut));
  hZ->SetLineColor(kBlack);
  hZ->SetMarkerStyle(20);
  hZ->SetMarkerColor(kBlack);
  //scale to get the shape
  hZ->Scale(1./hZ->Integral());
  if (variable.EqualTo("mJJ")) hZ->GetXaxis()->SetTitle(variable+" (GeV)");
  else hZ->GetXaxis()->SetTitle(variable);

  leg->AddEntry(hTT, "TTbar", "f");
  leg->AddEntry(hZ, "Zprime", "lep");

  hTT->Draw("hist");
  hZ->Draw("hist same E");
  hTT->SetMaximum(hTT->GetMaximum()* 2);
  leg->Draw();

  TString lumi_str = TString::Format("%0.1f", luminosity["luminosity"+year]/1000);
  lumi_13TeV = lumi_str+" fb^{-1}";

  //lumi_sqrtS = "13 TeV";
  int iPeriod = 4;
  int iPos = 0;
  //writeExtraText=true;
  CMS_lumi(can, iPeriod, iPos);
  can->Print(TString::Format("%s/SensitivityPlots/TTbar_Vs_Zprime_%s_%d_M%dW%d.pdf",year.Data(), variable.Data(), mJJCut, mass, width),"pdf");


  //chi2 test for comparison of 2 histograms
  Double_t chi2_result = hTT->Chi2Test(hZ, "WWP");
  cout<<chi2_result<<endl;
}
