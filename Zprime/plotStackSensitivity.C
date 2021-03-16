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

void plotStackSensitivity_Variable(TString year, TFile *infData, TFile *infTT, TFile *infQCD, TFile *infSub, TString variable, int mJJCut);

void plotStackSensitivity(TString year, int mJJCut = 2000, int selMass= 2500, int selWidth=25)
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
    plotStackSensitivity_Variable(year, infData, infTT, infQCD, infSub, varReco[ivar], mJJCut);
  }
}


void plotStackSensitivity_Variable(TString year, TFile *infData, TFile *infTT, TFile *infQCD, TFile *infSub, TString variable, int mJJCut)
{
  initFilesMapping();
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
  hTT->Scale(1./hTT->Integral());

  hQCD->SetLineColor(kBlue-6);
  hQCD->SetMarkerColor(kBlue-6);
  hQCD->SetFillColor(kBlue-6);
  hQCD->Scale(1./hQCD->Integral());

  hSub->SetLineColor(kTeal-7);
  hSub->SetMarkerColor(kTeal-7);
  hSub->SetFillColor(kTeal-7);
  hSub->Scale(1./hSub->Integral());


  THStack *hs;
  if(variable.EqualTo("mJJ")) hs = new THStack("Data vs MC", TString::Format("Data vs MC;%s (GeV);Shape", variable.Data()));
  else hs = new THStack("Data vs MC", TString::Format("Data vs MC;%s;Shape", variable.Data()));
  hs->Add(hSub);
  hs->Add(hQCD);
  hs->Add(hTT);

  TCanvas *can = new TCanvas(TString::Format("can_%s",variable.Data()), TString::Format("can_%s",variable.Data()), 800, 600);
  TLegend *leg;
  if(variable.Contains("cos"))leg = new TLegend(0.40,0.7,0.65,0.9);
  else leg = new TLegend(0.7,0.7,0.9,0.9);
  can->cd();

  //add the Zprime contribution
  TFile *infZprime;
  if(year.EqualTo("2016")) infZprime = TFile::Open(TString::Format("../Zprime/%s/HistoMassWindows_ZprimeToTT_M-%d_W-%d_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root", year.Data(), mass, width));
  else if (year.EqualTo("2017")) infZprime = TFile::Open(TString::Format("../Zprime/%s/HistoMassWindows_ZprimeToTT_M%d_W%d_TuneCP2_PSweights_13TeV-madgraph-pythiaMLM-pythia8_19UL.root", year.Data(), mass, width));  
  else infZprime = TFile::Open(TString::Format("../Zprime/%s/HistoMassWindows_ZprimeToTT_M%d_W%d_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8.root", year.Data(), mass, width));
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
  leg->AddEntry(hQCD, "QCD", "f");
  leg->AddEntry(hSub, "Subdominant", "f");
  leg->AddEntry(hZ, "Zprime", "lep");

  hs->Draw("hist");
  can->Update();
  hZ->Draw("hist same E");
  hs->SetMaximum(hs->GetMaximum()* 2);
  leg->Draw();

  TString lumi_str = TString::Format("%0.1f", luminosity[year]/1000);
  lumi_13TeV = lumi_str+" fb^{-1}";

  //lumi_sqrtS = "13 TeV";
  int iPeriod = 4;
  int iPos = 0;
  //writeExtraText=true;
  CMS_lumi(can, iPeriod, iPos);
  can->Print(TString::Format("%s/SensitivityPlots/DatavsMC_%s_%d_M%dW%d.pdf",year.Data(), variable.Data(), mJJCut, mass, width),"pdf");

}
