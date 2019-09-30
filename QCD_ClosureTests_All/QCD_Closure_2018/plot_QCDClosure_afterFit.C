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

void plot_QCDClosure_afterFit(TString recoVar = "jetPt0")
{
  //TT file 
  TFile *infTT = TFile::Open("SignalOutput_AllRegions_0.10_deepCSV.root");  
  //TFile *infTT = TFile::Open("SignalOutput_AllRegions_0.00_CSVv2.root");  
  //QCD file
  TFile *infBkg = TFile::Open("BkgOutput_AllRegions_0.10_deepCSV.root");
  //TFile *infBkg = TFile::Open("BkgOutput_AllRegions_0.00_CSVv2.root");
  TH1F *hBkg_CR, *hBkg_SR, *hBkg_CRExpYield;  
  TH1F *hSig_CR;    
  TH1F *hSig_1Btag, *hBkg_1Btag;
  
  //TH1F used for the TRatioPlot
  hBkg_CR = (TH1F*)infBkg->Get(TString::Format("CR_tTagger_%s",recoVar.Data()));
  hBkg_SR = (TH1F*)infBkg->Get(TString::Format("SR_tTagger_%s",recoVar.Data()));
  
  hBkg_1Btag = (TH1F*)infBkg->Get(TString::Format("h1Btag_tTagger_%s",recoVar.Data()));
  

  //get the fit result
  TFile *fitFile =  TFile::Open("FitOutput.root");
  TF1 *fitResult = (TF1*)fitFile->Get(TString::Format("func_%s",recoVar.Data()));
  //scale now all over the bins

  
  cout<<"beta [0]: "<<fitResult->GetParameter("beta")<<" ± "<<fitResult->GetParError(0)<<endl;
  cout<<"delta [1]: "<<fitResult->GetParameter("delta")<<" ± "<<fitResult->GetParError(1)<<endl;
  cout<<"alpha [2]: "<<fitResult->GetParameter("alpha")<<" ± "<<fitResult->GetParError(2)<<endl;


  int NBINS = hBkg_CR->GetNbinsX();

  for(int ibin=1; ibin<= NBINS; ibin++)
  {
    float binContent = hBkg_CR->GetBinContent(ibin);
    float chi = hBkg_CR->GetBinCenter(ibin);
    float SF = fitResult->Eval(chi);

    hBkg_CR->SetBinContent(ibin, binContent * SF);
  }
  
  //These are bkg Signal region and bkg region
  hBkg_CR->SetTitle("QCD Closure tTagger");
  hBkg_SR->SetTitle("QCD Closure tTagger");
  
  
  hBkg_CRExpYield = (TH1F*)infBkg->Get(TString::Format("CR_tTagger_%s_expYield",recoVar.Data()));
  hBkg_CRExpYield->SetLineColor(kRed);
  
  hBkg_CRExpYield->SetTitle("TT Contamination tTagger");
  
  hSig_CR = (TH1F*)infTT->Get(TString::Format("CR_tTagger_%s_expYield",recoVar.Data()));
  hSig_CR -> SetLineColor(kBlue);
  hSig_CR -> SetTitle("TT Contamination tTagger");
/*
  if(recoVar.EqualTo("jetMassSoftDrop"))
  {
    hBkg_CR[0]->Rebin(4);
    hBkg_CR[1]->Rebin(4);

    hBkg_SR[0]->Rebin(4);
    hBkg_SR[1]->Rebin(4);

    hBkg_1Btag[0]->Rebin(4);
    //hBkg_1Btag[1]->Rebin(4);

    hBkg_CRExpYield[0]->Rebin(4);
    hBkg_CRExpYield[1]->Rebin(4);

    hSig_CR[0]->Rebin(4);
    hSig_CR[1]->Rebin(4);

  }
  */
  TLegend *closureLegend = new TLegend(0.5,0.6,0.7,0.8);
  closureLegend->AddEntry(hBkg_SR,"Signal Region (2btag)", "l");
  closureLegend->AddEntry(hBkg_CR,"Control Region (0btag)", "f");
  closureLegend->AddEntry(hBkg_1Btag,"1btag Region", "l");
  
  auto c1 = new TCanvas("QCD closure Test tTagger", "QCD closure Test tTagger", 700,600);
  auto *closure_pad2 = new TPad("closure_pad2","closure_pad2",0.,0.,1.,0.3); 
  closure_pad2->Draw();
  closure_pad2->SetTopMargin(0.05);
  closure_pad2->SetBottomMargin(0.3);
  closure_pad2->SetGrid();

  auto *closure_pad1 = new TPad("closure_pad1","closure_pad1",0.,0.3,1.,1.);  
  closure_pad1->Draw();
  closure_pad1->SetBottomMargin(0.005);
  closure_pad1->cd();
  //closure_pad1->SetGrid();
  hBkg_SR->GetYaxis()->SetTitleSize(20);
  hBkg_SR->GetYaxis()->SetTitleFont(43);
  hBkg_SR->GetYaxis()->SetTitleOffset(1.4);  
  //hBkg_SR[0]->GetXaxis()->SetTitleOffset(1.5);
  hBkg_SR->GetXaxis()->SetTitle(recoVar+" (GeV)");
  // h2 settings
  
  //hBkg_CR[0]->ResetAttFill();
  hBkg_CR->ResetAttLine();
  hBkg_CR->ResetAttMarker();
  hBkg_CR->SetLineColor(33); 
  hBkg_CR->SetFillStyle(3001);
  hBkg_CR->SetFillColor(33);
  
  hBkg_1Btag->SetLineColor(kBlue);
  hBkg_SR->Draw();
  hBkg_CR->Draw("Hist same");
  hBkg_1Btag->Draw("same");
  closureLegend->Draw();
  
  
  
  TH1F *hClosure[2];
  closure_pad2->cd();
  hClosure[0] = (TH1F*)hBkg_SR->Clone("hClosure_0"); 
  hClosure[0]->SetTitle("");
  hClosure[0]->GetYaxis()->SetTitle("ratio SR/CR");
  hClosure[0]->GetYaxis()->SetTitleSize(14);
  hClosure[0]->GetYaxis()->SetTitleFont(43);
  hClosure[0]->GetYaxis()->SetTitleOffset(1.55);
  hClosure[0]->GetYaxis()->SetLabelFont(43);
  hClosure[0]->GetYaxis()->SetLabelSize(15);
  hClosure[0]->GetXaxis()->SetTitleSize(0.09);
  //hClosure[0]->GetXaxis()->SetTitleFont(43);
  //hClosure[0]->GetXaxis()->SetTitleOffset(4.);
  //hClosure[0]->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  hClosure[0]->GetXaxis()->SetLabelSize(0.09);
  hClosure[1] = (TH1F*)hBkg_1Btag->Clone("hClosure_1"); 
  hClosure[0]->Divide(hBkg_CR);
  hClosure[1]->Divide(hBkg_CR);
  hClosure[0]->SetLineColor(kRed);
  hClosure[1]->SetLineColor(kBlue);
  hClosure[0]->GetXaxis()->SetTitle(recoVar+" (GeV)");
  //hClosure[0]->GetXaxis()->SetTitleOffset(1);
  hClosure[0]->Draw();
  hClosure[1]->Draw("same");
  


  auto c3 = new TCanvas("TT contamination tTagger", "TT contamination tTagger", 700,600);
  auto rp_tTagger = new TRatioPlot(hSig_CR,hBkg_CRExpYield);
  c3->SetTicks(0,1);
  rp_tTagger->Draw();
  rp_tTagger->GetUpperRefYaxis()->SetTitle("Exp. Yield");
  rp_tTagger->GetLowerRefYaxis()->SetTitle("ratio");
  c3->Update();

}
