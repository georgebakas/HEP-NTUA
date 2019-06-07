#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TRatioPlot.h"

using std::cin;
using std::cout;
using std::endl;

void plot_QCDClosure(TString recoVar = "mJJ")
{
  //TT file 
  TFile *infTT = TFile::Open("SignalOutput_AllRegions.root");  
  //QCD file
  TFile *infBkg = TFile::Open("BkgOutput_AllRegions.root");
   
  TH1F *hBkg_CR[2], *hBkg_SR[2], *hBkg_CRExpYield[2];  
  TH1F *hSig_CR[2];    
  
  //TH1F used for the TRatioPlot
  hBkg_CR[0] = (TH1F*)infBkg->Get(TString::Format("CR_tTagger_%s",recoVar.Data()));
  hBkg_CR[1] = (TH1F*)infBkg->Get(TString::Format("CR_deepAK8_0.6_%s",recoVar.Data()));
  hBkg_SR[0] = (TH1F*)infBkg->Get(TString::Format("SR_tTagger_%s",recoVar.Data()));
  hBkg_SR[1] = (TH1F*)infBkg->Get(TString::Format("SR_deepAK8_0.6_%s",recoVar.Data()));
  
  cout<<hBkg_CR[0]->GetEntries()<<endl;
  cout<<hBkg_CR[1]->GetEntries()<<endl;
  
  hBkg_CRExpYield[0] = (TH1F*)infBkg->Get(TString::Format("CR_tTagger_%s_expYield",recoVar.Data()));
  hBkg_CRExpYield[1] = (TH1F*)infBkg->Get(TString::Format("CR_deepAK8_0.6_%s_expYield",recoVar.Data()));
  hBkg_CRExpYield[0]->SetLineColor(kRed);
  hBkg_CRExpYield[1]->SetLineColor(kRed);
  
  hSig_CR[0] = (TH1F*)infTT->Get(TString::Format("CR_tTagger_%s_expYield",recoVar.Data()));
  hSig_CR[1] = (TH1F*)infTT->Get(TString::Format("CR_deepAK8_0.6_%s_expYield",recoVar.Data()));
  
  TRatioPlot *trQCD[2];
  
  //trQCD[0] = (TRatioPlot*)infBkg->Get(TString::Format("TRatioPlot_tTagger_%s",recoVar.Data()));
  //trQCD[1] = (TRatioPlot*)infBkg->Get(TString::Format("TRatioPlot_deepAK8_0.6_%s",recoVar.Data()));
  
  auto c1 = new TCanvas("QCD closure Test tTagger", "QCD closure Test tTagger", 700,600);
  trQCD[0] = new TRatioPlot(hBkg_SR[0],hBkg_CR[0]);
  c1->SetTicks(0,1);
  trQCD[0]->Draw();
  c1->Update();
  trQCD[0]->GetLowerRefYaxis()->SetTitle("ratio");
  trQCD[0]->GetUpperRefYaxis()->SetTitle("Exp. Yield");
  
  auto c2 = new TCanvas("QCD closure Test deepAK8", "QCD closure Test deepAK8", 700,600);
  c2->SetTicks(0,1);
  trQCD[1] = new TRatioPlot(hBkg_SR[1],hBkg_CR[1]);
  trQCD[1]->Draw();
  c2->Update();
  trQCD[1]->GetLowerRefYaxis()->SetTitle("ratio");
  trQCD[1]->GetUpperRefYaxis()->SetTitle("Exp. Yield");
   

  auto c3 = new TCanvas("TT contamination tTagger", "TT contamination tTagger", 700,600);
  auto rp_tTagger = new TRatioPlot(hSig_CR[0],hBkg_CRExpYield[0] );
  c3->SetTicks(0,1);
  rp_tTagger->Draw();
  c3->Update();

  auto c4 = new TCanvas("TT contamination deepAK8", "TT contamination deepAK8", 700,600);
  auto rp_deepAK8 = new TRatioPlot(hSig_CR[1],hBkg_CRExpYield[1]);
  c4->SetTicks(0,1);
  rp_deepAK8->Draw();
  c4->Update();

}