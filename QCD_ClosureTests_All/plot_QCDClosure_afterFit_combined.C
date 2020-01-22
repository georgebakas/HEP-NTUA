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

#include "TemplateConstants.h"
void plotYearVar(TString year, TString recoVar = "jetPt0");
void ratioPlot(TString year, TH1F *hNum ,TH1F *hDenom, TString recoVar, TString reason, bool isClosure);

void plot_QCDClosure_afterFit_combined(TString year = "2016")
{
	const int NVAR = 7;
	TString recoVar[NVAR] = {"jetPt0", "mJJ", "ptJJ", "yJJ", "jetPt1", "jetMassSoftDrop0", "jetMassSoftDrop1"};

	for(int ivar = 0; ivar<NVAR; ivar++)
	{
		plotYearVar(year,recoVar[ivar]);
	}
}
void plotYearVar(TString year, TString recoVar = "jetPt0")
{
  initFilesMapping();
   
  //I need 2 files for each:
  //1. Signal Region files which contain everything with medium WP's
  //2. Control Region files which contain everything with loose WP's
  TString yearSignalRegion = year;
  TString yearControlRegion = year + "_Loose";
  
  //TT files
  TFile *infTT = TFile::Open(filesSignal[yearSignalRegion.Data()]);  
  TFile *infTTLoose = TFile::Open(filesSignal[yearControlRegion.Data()]);  
  
  //QCD files 
  TFile *infBkg = TFile::Open(filesBkg[yearSignalRegion.Data()]);
  TFile *infBkgLoose = TFile::Open(filesBkg[yearControlRegion.Data()]);
  TH1F *hBkg_CR, *hBkg_SR;  
  TH1F *hSig_CRExpYield, *hBkg_CRExpYield;    
  
  //for closure test
  hBkg_CR = (TH1F*)infBkgLoose->Get(TString::Format("CR_tTagger_%s_expYield",recoVar.Data())); //from file with loose WP's
  hBkg_SR = (TH1F*)infBkg->Get(TString::Format("SR_tTagger_%s_expYield",recoVar.Data())); //from file with medium WP's

  //These are bkg Signal region and bkg region
  hBkg_CR->SetTitle("QCD_{CR} 0btag (Loose)");
  hBkg_SR->SetTitle("QCD_{SR} 2btag (Medium)");
  
  hBkg_CR->Scale(1./hBkg_CR->Integral());
  hBkg_SR->Scale(1./hBkg_SR->Integral());
  cout<<"hBkg_CR Loose Entries: "<<hBkg_CR->GetEntries()<<endl;
  cout<<"hBkg_SR Medium Entries: "<<hBkg_SR->GetEntries()<<endl;

  
  //now for contamination
  hBkg_CRExpYield = (TH1F*)infBkgLoose->Get(TString::Format("CR_tTagger_%s_expYield",recoVar.Data()));
  hBkg_CRExpYield->SetLineColor(kRed);
  
  hBkg_CRExpYield->SetTitle("QCD_{CR} Exp.Yield");
  
  hSig_CRExpYield = (TH1F*)infTTLoose->Get(TString::Format("CR_tTagger_%s_expYield",recoVar.Data()));
  hSig_CRExpYield->SetLineColor(kBlue);
  hSig_CRExpYield->SetTitle("TT_{CR} Exp.Yield");

  if(recoVar.EqualTo("jetMassSoftDrop0") || recoVar.EqualTo("jetMassSoftDrop1"))
  {
    hBkg_CR->Rebin(5);
    hBkg_SR->Rebin(5);
    hBkg_CRExpYield->Rebin(5);
    hSig_CRExpYield->Rebin(5);
  }
  
  cout<<"hBkg_CR_expYield Loose Entries: "<<hBkg_CRExpYield->GetEntries()<<endl;
  cout<<"hSig_CR_expYield Loose Entries: "<<hSig_CRExpYield->GetEntries()<<endl;
  //this is ttbar contamination
  TString reas = "TTbar Contamination";
  ratioPlot(year, hSig_CRExpYield,hBkg_CRExpYield, recoVar, reas, false);

  //this is closure test
  reas = "QCD Closure";
  ratioPlot(year, hBkg_SR, hBkg_CR,recoVar, reas, true);

}



void ratioPlot(TString year, TH1F *hNum ,TH1F *hDenom, TString recoVar, TString reason, bool isClosure)
{
  TString titleNum = hNum->GetTitle();
  TString titleDenom = hDenom->GetTitle();
  TLegend *closureLegend = new TLegend(0.65,0.65,0.9,0.9);
  closureLegend->AddEntry(hNum,titleNum, "lep");
  closureLegend->AddEntry(hDenom,titleDenom, "lep");

  if(isClosure)
  {
  	hNum->SetTitle("QCD Closure");
  	hDenom->SetTitle("QCD Closure");
  }
  else
  {
  	hNum->SetTitle("TTbar Contamination");
  	hDenom->SetTitle("TTbar Contamination");
  }
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
  
  if(isClosure)
  {	
  	//this is for qcd closure
  	cout<<"CLOSURE"<<endl;
 	hNum->Draw();
  	hDenom->Draw("same");
  }
  else
  {
  	cout<<"ttbar contamination"<<endl;
  	hDenom->Draw();
  	hNum->Draw("same");
  	hDenom->GetYaxis()->SetRangeUser(10E-2,hDenom->GetMaximum()+200);
  }
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


  if(isClosure) hRatio->GetYaxis()->SetRangeUser(0,2);
  hRatio->GetXaxis()->SetLabelSize(0.06);
  hRatio->Divide(hDenom);
  hRatio->SetLineColor(kBlack);
  if(!recoVar.EqualTo("yJJ")) hRatio->GetXaxis()->SetTitle(recoVar+" (GeV)");
  else hRatio->GetXaxis()->SetTitle(recoVar);
  //hRatio->GetXaxis()->SetTitleOffset(1);
  hRatio->Draw();
  if(isClosure) c1->Print(TString::Format("../TopTaggerEfficiencies/plotsCombined_LooseCR_MediumSR/%s/qcdClosure_%s.pdf",year.Data(),recoVar.Data()),"pdf");
  else  c1->Print(TString::Format("../TopTaggerEfficiencies/plotsCombined_LooseCR_MediumSR/%s/ttContamination_%s.pdf",year.Data(),recoVar.Data()),"pdf");
}
//needed for fit
/*
  //get the fit result
  TFile *fitFile; 
  TF1 *fitResult; 
  bool useSF=false;
  
  if((recoVar.EqualTo("mJJ") || recoVar.EqualTo("jetPt1") || recoVar.EqualTo("jetPt0"))
      && !year.EqualTo("2016_Loose") && !year.EqualTo("2016"))
  {
    cout<<"in"<<endl;
    useSF = true; 
    fitFile =  TFile::Open(TString::Format("QCD_Closure_%s/FitOutput.root",year.Data()));
    fitResult = (TF1*)fitFile->Get(TString::Format("func_%s",recoVar.Data()));

  }
  if(year.EqualTo("2017") && recoVar.EqualTo("ptJJ"))
  {
    useSF = true; 
    fitFile =  TFile::Open(TString::Format("QCD_Closure_%s/FitOutput.root",year.Data()));
    fitResult = (TF1*)fitFile->Get(TString::Format("func_%s",recoVar.Data()));
  } 
   
  TH1F *hBkg_CR_uncorrected = (TH1F*)hBkg_CR->Clone("hBkg_CR_uncorrected");
  int NBINS = hBkg_CR->GetNbinsX();
  float SF =1;
  for(int ibin=1; ibin<= NBINS; ibin++)
  {
    float binContent = hBkg_CR->GetBinContent(ibin);
    float chi = hBkg_CR->GetBinCenter(ibin);
    if(useSF)  SF = fitResult->Eval(chi);

    hBkg_CR->SetBinContent(ibin, binContent * SF);
  }
  */
