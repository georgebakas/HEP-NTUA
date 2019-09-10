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
  TFile *infTT = TFile::Open("SignalOutput_AllRegions_0.10_deepCSV.root");  
  //QCD file
  TFile *infBkg = TFile::Open("BkgOutput_AllRegions_0.10_deepCSV.root");
   
  TH1F *hBkg_CR[2], *hBkg_SR[2], *hBkg_CRExpYield[2];  
  TH1F *hSig_CR[2];    
  TH1F *hSig_1Btag[2], *hBkg_1Btag[2];
  
  //TH1F used for the TRatioPlot
  hBkg_CR[0] = (TH1F*)infBkg->Get(TString::Format("CR_tTagger_%s",recoVar.Data()));
  hBkg_CR[1] = (TH1F*)infBkg->Get(TString::Format("CR_deepAK8_0.60_%s",recoVar.Data()));
  hBkg_SR[0] = (TH1F*)infBkg->Get(TString::Format("SR_tTagger_%s",recoVar.Data()));
  hBkg_SR[1] = (TH1F*)infBkg->Get(TString::Format("SR_deepAK8_0.60_%s",recoVar.Data()));
  
  hBkg_1Btag[0] = (TH1F*)infBkg->Get(TString::Format("h1Btag_tTagger_%s",recoVar.Data()));
  hBkg_1Btag[1] = (TH1F*)infBkg->Get(TString::Format("h1Btag_deepAK8_%s",recoVar.Data()));
  
  
  //These are bkg Signal region and bkg region
  hBkg_CR[0]->SetTitle("QCD Closure tTagger");
  hBkg_SR[0]->SetTitle("QCD Closure tTagger");
  hBkg_CR[1]->SetTitle("QCD Closure deepAK8");
  hBkg_SR[1]->SetTitle("QCD Closure deepAK8");
  
  
  hBkg_CRExpYield[0] = (TH1F*)infBkg->Get(TString::Format("CR_tTagger_%s_expYield",recoVar.Data()));
  hBkg_CRExpYield[1] = (TH1F*)infBkg->Get(TString::Format("CR_deepAK8_0.60_%s_expYield",recoVar.Data()));
  hBkg_CRExpYield[0]->SetLineColor(kRed);
  hBkg_CRExpYield[1]->SetLineColor(kRed);
  
  hBkg_CRExpYield[0]->SetTitle("TT Contamination tTagger");
  hBkg_CRExpYield[1]->SetTitle("TT Contamination deepAK8");
  
  hSig_CR[0] = (TH1F*)infTT->Get(TString::Format("CR_tTagger_%s_expYield",recoVar.Data()));
  hSig_CR[0] -> SetLineColor(kBlue);
  hSig_CR[0] -> SetTitle("TT Contamination tTagger");
  hSig_CR[1] = (TH1F*)infTT->Get(TString::Format("CR_deepAK8_0.60_%s_expYield",recoVar.Data()));
  hSig_CR[1] -> SetLineColor(kBlue);
  hSig_CR[1] -> SetTitle("TT Contamination deepAK8");
  
  TLegend *closureLegend = new TLegend(0.5,0.6,0.7,0.8);
  closureLegend->AddEntry(hBkg_SR[0],"Signal Region (2btag)", "l");
  closureLegend->AddEntry(hBkg_CR[0],"Control Region (0btag)", "l");
  closureLegend->AddEntry(hBkg_1Btag[0],"1btag Region", "l");
  
  auto c1 = new TCanvas("QCD closure Test tTagger", "QCD closure Test tTagger", 700,600);
  auto *closure_pad2 = new TPad("closure_pad2","closure_pad2",0.,0.,1.,0.3); 
  closure_pad2->Draw();
  closure_pad2->SetTopMargin(0.05);
  closure_pad2->SetBottomMargin(0.1);
  closure_pad2->SetGrid();

  auto *closure_pad1 = new TPad("closure_pad1","closure_pad1",0.,0.3,1.,1.);  
  closure_pad1->Draw();
  closure_pad1->SetBottomMargin(0.001);
  closure_pad1->cd();
  //closure_pad1->SetGrid();
  hBkg_SR[0]->GetYaxis()->SetTitleSize(20);
  hBkg_SR[0]->GetYaxis()->SetTitleFont(43);
  hBkg_SR[0]->GetYaxis()->SetTitleOffset(1.55);  
  // h2 settings
  
  //hBkg_CR[0]->ResetAttFill();
  hBkg_CR[0]->ResetAttLine();
  hBkg_CR[0]->ResetAttMarker();
  hBkg_CR[0]->SetLineColor(33); 
  hBkg_CR[0]->SetFillStyle(3001);
  hBkg_CR[0]->SetFillColor(33);
  
  hBkg_SR[0]->Draw();
  hBkg_CR[0]->Draw("Hist same");
  hBkg_1Btag[0]->Draw("same");
  closureLegend->Draw();
  
  
  
  TH1F *hClosure[2];
  closure_pad2->cd();
  hClosure[0] = (TH1F*)hBkg_SR[0]->Clone("hClosure_0"); 
  hClosure[0]->SetTitle("");
  hClosure[0]->GetYaxis()->SetTitle("ratio SR/CR");
  hClosure[0]->GetYaxis()->SetTitleSize(14);
  hClosure[0]->GetYaxis()->SetTitleFont(43);
  hClosure[0]->GetYaxis()->SetTitleOffset(1.55);
  hClosure[0]->GetYaxis()->SetLabelFont(43);
  hClosure[0]->GetYaxis()->SetLabelSize(15);
  hClosure[0]->GetXaxis()->SetTitleSize(1);
  hClosure[0]->GetXaxis()->SetTitleFont(43);
  hClosure[0]->GetXaxis()->SetTitleOffset(4.);
  hClosure[0]->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  hClosure[0]->GetXaxis()->SetLabelSize(15);
  hClosure[1] = (TH1F*)hBkg_1Btag[0]->Clone("hClosure_1"); 
  hClosure[0]->Divide(hBkg_CR[0]);
  hClosure[1]->Divide(hBkg_CR[0]);
  hClosure[0]->SetLineColor(kBlack);
  hClosure[1]->SetLineColor(kRed);
  hClosure[0]->Draw();
  hClosure[1]->Draw("same");
  
  /*
  
  auto c1 = new TCanvas("QCD closure Test tTagger", "QCD closure Test tTagger", 700,600);
  trQCD[0] = new TRatioPlot(hBkg_SR[0],hBkg_CR[0]);
  c1->SetTicks(0,1);
  trQCD[0]->Draw();
  c1->Update();
  trQCD[0]->GetLowerRefYaxis()->SetTitle("ratio");
  
  
  auto c2 = new TCanvas("QCD closure Test deepAK8", "QCD closure Test deepAK8", 700,600);
  c2->SetTicks(0,1);
  trQCD[1] = new TRatioPlot(hBkg_SR[1],hBkg_CR[1]);
  trQCD[1]->Draw();
  c2->Update();
  trQCD[1]->GetLowerRefYaxis()->SetTitle("ratio");
  
  auto c4 = new TCanvas("TT contamination deepAK8", "TT contamination deepAK8", 700,600);
  auto rp_deepAK8 = new TRatioPlot(hSig_CR[1],hBkg_CRExpYield[1]);
  c4->SetTicks(0,1);
  rp_deepAK8->Draw();
  rp_deepAK8->GetUpperRefYaxis()->SetTitle("Exp. Yield");
  rp_deepAK8->GetLowerRefYaxis()->SetTitle("ratio");
  c4->Update();
 */  

  auto c3 = new TCanvas("TT contamination tTagger", "TT contamination tTagger", 700,600);
  auto rp_tTagger = new TRatioPlot(hSig_CR[0],hBkg_CRExpYield[0]);
  c3->SetTicks(0,1);
  rp_tTagger->Draw();
  rp_tTagger->GetUpperRefYaxis()->SetTitle("Exp. Yield");
  rp_tTagger->GetLowerRefYaxis()->SetTitle("ratio");
  c3->Update();

  
}
