#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRatioPlot.h"

using std::cin;
using std::cout;
using std::endl;

#include "TemplateConstants.h"
void plotYearVar(TString year, TString recoVar = "jetPt0", bool useScaleFactor= true, TString fitRecoVar= "leadingJetPt");
void ratioPlot(TString year, TH1F *hNum ,TH1F *hDenom, TString recoVar, TString reason, bool useScaleFactor, TF1 *fitResult);

void plot_QCDClosure_afterFit_combined(TString year = "2016", bool useScaleFactor = true)
{
  const int NVAR = 9;
  TString recoVar[NVAR] = {"jetPt0", "mJJ", "ptJJ", "yJJ", "jetPt1", "mTop", "jetMassSoftDrop", "jetY0", "jetY1"};
  TString fitRecoVar[NVAR] = {"leadingJetPt", "mJJ", "ptJJ", "yJJ", "subleadingJetPt", "leadingJetMassSoftDrop", "subleadingJetMassSoftDrop",
                              "leadingJetY", "subleadingJetY"};

  for(int ivar = 0; ivar<NVAR; ivar++)
  {
    plotYearVar(year,recoVar[ivar], useScaleFactor, fitRecoVar[ivar]);
    //break;
  }
}
void plotYearVar(TString year, TString recoVar = "jetPt0", bool useScaleFactor= true, TString fitRecoVar= "leadingJetPt")
{

  //I need 2 files for each:
  //1. Signal and Control Region files which contain everything with medium WP's
  TString yearSignalRegion = year;
  TString yearControlRegion = year;

  //TT files
  TFile *infTT = TFile::Open(TString::Format("../MassFit/%s/Histo_TT_NominalMC_100_reduced_UnequalBinning.root",year.Data()));

  //QCD files
  TFile *infBkg = TFile::Open(TString::Format("../MassFit/%s/Histo_QCD_HT300toInf_100_reduced_UnequalBinning.root",year.Data()));
  TH1F *hBkg_CR, *hBkg_SR;
  TH1F *hSig_CRExpYield, *hBkg_CRExpYield;

  //for closure test
  hBkg_CR = (TH1F*)infBkg->Get(TString::Format("hWt_%s_0btag_expYield",recoVar.Data()));
  hBkg_SR = (TH1F*)infBkg->Get(TString::Format("hWt_%s_2btag_expYield",recoVar.Data()));

  //These are bkg Signal region and bkg region
  hBkg_CR->SetTitle("QCD_{CR} 0btag");
  hBkg_SR->SetTitle("QCD_{SR} 2btag");

  hBkg_CR->Scale(1./hBkg_CR->Integral());
  hBkg_SR->Scale(1./hBkg_SR->Integral());
  cout<<"hBkg_CR Entries: "<<hBkg_CR->GetEntries()<<endl;
  cout<<"hBkg_SR Entries: "<<hBkg_SR->GetEntries()<<endl;


  if(recoVar.EqualTo("mTop") || recoVar.EqualTo("jetMassSoftDrop"))
  {
    hBkg_CR->Rebin(4);
    hBkg_SR->Rebin(4);
  }


  //get the fit result

  TFile *fitFile;
  if(recoVar.EqualTo("mTop") || recoVar.EqualTo("jetMassSoftDrop")) fitFile = TFile::Open(TString::Format("closureTest_fitResults_%s_extended.root",year.Data()));
  else fitFile = TFile::Open(TString::Format("closureTest_fitResults_%s_reduced.root",year.Data()));
  TF1 *fitResult = (TF1*)fitFile->Get(TString::Format("FitFunction_%s",fitRecoVar.Data()));

  TString reas = "QCDClosure";
  TF1 *tempFR2;

  if(useScaleFactor)
  {
    int NBINS = hBkg_CR->GetNbinsX();
    float SF;
    if(((recoVar.Contains("jetPt") || recoVar.EqualTo("mJJ")) && !year.EqualTo("2016")) || (recoVar.EqualTo("jetPt0") && year.EqualTo("2016")))
    {
      for(int ibin=1; ibin<= NBINS; ibin++)
      {
        cout<<"----------"<<endl;
        cout<<"before fit: "<<hBkg_SR->GetBinContent(ibin)/hBkg_CR->GetBinContent(ibin)<<endl;
        float binContent = hBkg_CR->GetBinContent(ibin);
        float chi = hBkg_CR->GetBinCenter(ibin);
        SF = fitResult->Eval(chi);
        //cout<<SF<<endl;
        hBkg_CR->SetBinContent(ibin, binContent * SF);
        cout<<"after fit: "<<hBkg_SR->GetBinContent(ibin)/hBkg_CR->GetBinContent(ibin)<<endl;
      }
    }

    ratioPlot(year, hBkg_SR, hBkg_CR,recoVar, reas, useScaleFactor, fitResult);
  }
  else
  {
    if(((recoVar.Contains("jetPt") || recoVar.EqualTo("mJJ")) && !year.EqualTo("2016")) || (recoVar.EqualTo("jetPt0") && year.EqualTo("2016"))) ratioPlot(year, hBkg_SR, hBkg_CR,recoVar, reas, false, fitResult);
  }


}



void ratioPlot(TString year, TH1F *hNum ,TH1F *hDenom, TString recoVar, TString reason, bool useScaleFactor, TF1 *fitResult)
{
  cout<<"after fitINSIDE: "<<hNum->GetBinContent(9)/hDenom->GetBinContent(9)<<endl;
  TString titleNum = hNum->GetTitle();
  TString titleDenom = hDenom->GetTitle();
  TLegend *closureLegend = new TLegend(0.65,0.65,0.9,0.9);
  closureLegend->AddEntry(hNum,titleNum, "lep");
  closureLegend->AddEntry(hDenom,titleDenom, "lep");

  hNum->SetTitle("QCD Closure");
  hDenom->SetTitle("QCD Closure");


  auto c1 = new TCanvas(reason+recoVar, reason+recoVar, 800,700);
  auto *closure_pad2 = new TPad("closure_pad2","closure_pad2",0.,0.,1.,0.4);
  closure_pad2->Draw();
  closure_pad2->SetTopMargin(0.05);
  closure_pad2->SetBottomMargin(0.3);
  closure_pad2->SetGrid();

  auto *closure_pad1 = new TPad("closure_pad1","closure_pad1",0.,0.4,1.,1.);
  closure_pad1->Draw();
  closure_pad1->SetBottomMargin(0.005);
  closure_pad1->cd();
  //closure_pad1->SetGrid();
  hNum->GetYaxis()->SetTitleSize(20);
  hNum->GetYaxis()->SetTitleFont(43);
  hNum->GetYaxis()->SetTitleOffset(1.4);
  hNum->ResetAttLine();
  hNum->ResetAttMarker();
  hNum->SetLineColor(kRed-7);
  hNum->SetMarkerStyle(22);
  hNum->SetMarkerColor(kRed-7);
  //hBkg_SR[0]->GetXaxis()->SetTitleOffset(1.5);
  if(!recoVar.EqualTo("yJJ")) hNum->GetXaxis()->SetTitle(recoVar+" (GeV)");
  else hNum->GetXaxis()->SetTitle(recoVar);
  // h2 settings

  //hBkg_CR[0]->ResetAttFill();
  hDenom->ResetAttLine();
  hDenom->ResetAttMarker();
  hDenom->SetLineColor(kBlue-2);
  hDenom->SetMarkerStyle(21);
  hDenom->SetMarkerColor(kBlue-2);


  cout<<"CLOSURE"<<endl;
  hNum->Draw();
  hDenom->Draw("same");

  closureLegend->Draw();
  closure_pad1->SetLogy();


  TH1F *hRatio;
  closure_pad2->cd();
  hRatio = (TH1F*)hNum->Clone("hNumerator_0");
  hRatio->SetTitle("");
  hRatio->ResetAttMarker();
  //hRatio->GetYaxis()->SetTitle(TString::Format("ratio %s/%s",hNum->GetTitle(),hDenom->GetTitle()));
  hRatio->GetYaxis()->SetTitle("ratio Red/Blue");
  hRatio->GetYaxis()->SetTitleSize(20);
  hRatio->GetYaxis()->SetTitleFont(43);
  hRatio->GetYaxis()->SetTitleOffset(1.55);
  hRatio->GetYaxis()->SetLabelFont(43);
  hRatio->GetYaxis()->SetLabelSize(15);
  hRatio->GetXaxis()->SetTitleSize(0.06);


  hRatio->GetYaxis()->SetRangeUser(0,2);
  hRatio->GetXaxis()->SetLabelSize(0.06);
  hRatio->Divide(hDenom);
  hRatio->SetLineColor(kBlack);
  if(!recoVar.EqualTo("yJJ")) hRatio->GetXaxis()->SetTitle(recoVar+" (GeV)");
  else hRatio->GetXaxis()->SetTitle(recoVar);
  //hRatio->GetXaxis()->SetTitleOffset(1);

  if(useScaleFactor)
  {
    hRatio->Draw();
    recoVar = recoVar+"_ScaledWithFitResults";
  }
  else
  {
    hRatio->Fit(fitResult, "R");
    hRatio->Draw();
    fitResult->Draw("same");
  }
  c1->Print(TString::Format("./plotsMedium/%s/qcdClosure_%s.pdf",year.Data(),recoVar.Data()),"pdf");
}
