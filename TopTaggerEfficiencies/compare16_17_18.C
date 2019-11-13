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

void compare16_17_18(TString recoVar = "jetPt0",TString partonVar = "partonPt0", bool isDeepCSV = true)
{
   TFile *eff[3], *signal[3], *bkg[3];
   
   //check what year you want to compare 
   //if compare2016 is true then you compare just 2016
   std::vector<Color_t> colors = {kBlue,kRed, kGreen,kBlack, kMagenta};
    cout<<"It's deepCSV"<<endl;
	   //for closure tests [0] is 16, [1] is 17 and [2] is 18
	   signal[0] = TFile::Open("./QCD_Closure/SignalOutput_AllRegions_0.2_deepCSV_1File.root"); 
	   signal[1] = TFile::Open("./QCD_Closure_2017/SignalOutput_AllRegions_0.2_deepCSV_1File.root");
	   signal[2] = TFile::Open("./QCD_Closure_2018/SignalOutput_AllRegions_0.2_deepCSV_1File.root");

	   bkg[0] = TFile::Open("./QCD_Closure/BkgOutput_AllRegions_0.2_deepCSV_1File.root");
	   bkg[1] = TFile::Open("./QCD_Closure_2017/BkgOutput_AllRegions_0.2_deepCSV_1File.root");
	   bkg[2] = TFile::Open("./QCD_Closure_2018/BkgOutput_AllRegions_0.2_deepCSV_1File.root");
	   
 	   //for efficiencies
     eff[0] = TFile::Open("./TopTagger_Efficiencies/Efficiencies_allVars_tTagger_0.20_deepCSV.root"); 
	   eff[1] = TFile::Open("./TopTagger_Efficiencies_2017/Efficiencies_allVars_tTagger_0.00_deepCSV.root");
	   eff[2] = TFile::Open("./TopTagger_Efficiencies_2018/Efficiencies_allVars_tTagger_0.10_deepCSV.root");

     TFile *oldInf = TFile::Open("PartonEfficiencyAll_July19.root");
   
   
   
   //things we need:
   TEfficiency *eff16[2], *acc16[2];   
   TEfficiency *eff17[2], *acc17[2];
   TEfficiency *eff18[2], *acc18[2];

   TEfficiency *effOld16, *accOld16;
   effOld16 = (TEfficiency*)oldInf->Get("Eff_jetPt0_Parton_Nominal");
   accOld16 = (TEfficiency*)oldInf->Get("Eff_jetPt0_Parton_Nominal_common");
   
   //1. Efficiency for 2016 and acceptance for same year
   eff16[0] = (TEfficiency*)eff[0]->Get(TString::Format("Sig_Parton_tTagger_0.20_%s", partonVar.Data()));
   acc16[0] = (TEfficiency*)eff[0]->Get(TString::Format("Sig_Reco_tTagger_0.20_%s", recoVar.Data()));
   
   eff16[1] = (TEfficiency*)eff[0]->Get(TString::Format("Sig_Parton_oldMva_0.8_%s", partonVar.Data()));
   acc16[1] = (TEfficiency*)eff[0]->Get(TString::Format("Sig_Reco_oldMva_0.8_%s", recoVar.Data()));
   
   //2. Efficiency for 2017 and acceptance for same year both for tTagger and eventTagger
   eff17[0] = (TEfficiency*)eff[1]->Get(TString::Format("Sig_Parton_tTagger_0.00_%s", partonVar.Data()));
   acc17[0] = (TEfficiency*)eff[1]->Get(TString::Format("Sig_Reco_tTagger_0.00_%s", recoVar.Data()));
   
   eff17[1] = (TEfficiency*)eff[1]->Get(TString::Format("Sig_Parton_oldMva_0.80_%s", partonVar.Data()));
   acc17[1] = (TEfficiency*)eff[1]->Get(TString::Format("Sig_Reco_oldMva_0.80_%s", recoVar.Data()));
   
   //3. Efficiency for 2018 and acceptance for same year
   eff18[0] = (TEfficiency*)eff[2]->Get(TString::Format("Sig_Parton_tTagger_0.10_%s", partonVar.Data()));
   acc18[0] = (TEfficiency*)eff[2]->Get(TString::Format("Sig_Reco_tTagger_0.10_%s", recoVar.Data()));
   
   eff18[1] = (TEfficiency*)eff[2]->Get(TString::Format("Sig_Parton_oldMva_0.80_%s", partonVar.Data()));
   acc18[1] = (TEfficiency*)eff[2]->Get(TString::Format("Sig_Reco_oldMva_0.80_%s", recoVar.Data()));
  
   
   eff16[0]->SetMarkerStyle(21);
   eff17[0]->SetMarkerStyle(22);
   eff18[0]->SetMarkerStyle(23);
   effOld16->SetMarkerStyle(20);
   eff16[0]->SetMarkerColor(colors[0]);
   eff17[0]->SetMarkerColor(colors[1]);
   eff18[0]->SetMarkerColor(colors[2]);
   effOld16->SetMarkerColor(colors[3]);
   eff16[0]->SetLineColor(colors[0]);
   eff17[0]->SetLineColor(colors[1]);
   eff18[0]->SetLineColor(colors[2]);
   effOld16->SetLineColor(colors[3]);
   
   acc16[0]->SetMarkerStyle(21);
   acc17[0]->SetMarkerStyle(22);
   acc18[0]->SetMarkerStyle(23);
   accOld16->SetMarkerStyle(20);
   acc16[0]->SetLineColor(colors[0]);
   acc17[0]->SetLineColor(colors[1]);
   acc18[0]->SetLineColor(colors[2]);
   acc16[0]->SetMarkerColor(colors[0]);
   acc17[0]->SetMarkerColor(colors[1]);
   acc18[0]->SetMarkerColor(colors[2]);
   accOld16->SetMarkerColor(colors[3]);
   accOld16->SetLineColor(colors[3]);
   
   TLegend *effLeg = new TLegend(0.5,0.6,0.7,0.8);
   effLeg->AddEntry(eff16[0], "tTagger '16", "lp");
   effLeg->AddEntry(eff17[0], "tTagger '17", "lp");
   effLeg->AddEntry(eff18[0], "tTagger '18", "lp");
   effLeg->AddEntry(effOld16, "2016 analysis", "lp");

   TCanvas *can_eff = new TCanvas("Efficiency can", "Efficiency can", 700, 600);
   eff18[0]->SetTitle(TString::Format("Efficiency '16,'17,'18;%s (GeV);Efficiency",recoVar.Data())); 
   eff18[0]->Draw();
   eff17[0]->Draw("same");
   eff16[0]->Draw("same");
   effOld16->Draw("same");
   effLeg->Draw();
   
   
   TCanvas *can_acc = new TCanvas("Acceptance can", "Acceptance can", 700, 600);
   acc18[0]->SetTitle(TString::Format("Acceptance '16,'17,'18;%s (GeV);Acceptance",recoVar.Data()));  
   acc18[0]->Draw();
   acc17[0]->Draw("same");
   acc16[0]->Draw("same");
   accOld16->Draw("same");
   effLeg->Draw();
   

   //now we compare the closure tests for the 2016 and 2017 MC 
  TH1F *hBkg_CR[3], *hBkg_SR[3], *hBkg_CRExpYield[3];  
  TH1F *hSig_CR[3];    
  TH1F *hSig_1Btag[3], *hBkg_1Btag[3];
   //TH1F used for the TRatioPlot [0] is 16 and [1] is 2017 and [2] is 2018
  hBkg_CR[0] = (TH1F*)bkg[0]->Get(TString::Format("CR_tTagger_%s",recoVar.Data()));
  hBkg_CR[1] = (TH1F*)bkg[1]->Get(TString::Format("CR_tTagger_%s",recoVar.Data()));
  hBkg_CR[2] = (TH1F*)bkg[2]->Get(TString::Format("CR_tTagger_%s",recoVar.Data()));
  hBkg_SR[0] = (TH1F*)bkg[0]->Get(TString::Format("SR_tTagger_%s",recoVar.Data()));
  hBkg_SR[1] = (TH1F*)bkg[1]->Get(TString::Format("SR_tTagger_%s",recoVar.Data()));
  hBkg_SR[2] = (TH1F*)bkg[2]->Get(TString::Format("SR_tTagger_%s",recoVar.Data()));
  
  
  hBkg_1Btag[0] = (TH1F*)bkg[0]->Get(TString::Format("h1Btag_tTagger_%s",recoVar.Data()));
  hBkg_1Btag[1] = (TH1F*)bkg[1]->Get(TString::Format("h1Btag_tTagger_%s",recoVar.Data()));  
  hBkg_1Btag[2] = (TH1F*)bkg[2]->Get(TString::Format("h1Btag_tTagger_%s",recoVar.Data()));  
   
  hBkg_CRExpYield[0] = (TH1F*)bkg[0]->Get(TString::Format("CR_tTagger_%s_expYield",recoVar.Data()));
  hBkg_CRExpYield[1] = (TH1F*)bkg[1]->Get(TString::Format("CR_tTagger_%s_expYield",recoVar.Data()));
  hBkg_CRExpYield[2] = (TH1F*)bkg[2]->Get(TString::Format("CR_tTagger_%s_expYield",recoVar.Data()));
  hBkg_CRExpYield[0]->SetLineColor(kRed);
  hBkg_CRExpYield[1]->SetLineColor(kRed);
  hBkg_CRExpYield[2]->SetLineColor(kRed);
  
  hBkg_CRExpYield[0]->SetTitle("TT Contamination tTagger '16");
  hBkg_CRExpYield[1]->SetTitle("TT Contamination tTagger '17");
  hBkg_CRExpYield[2]->SetTitle("TT Contamination tTagger '18");
  
  hSig_CR[0] = (TH1F*)signal[0]->Get(TString::Format("CR_tTagger_%s_expYield",recoVar.Data()));
  hSig_CR[0] -> SetLineColor(kBlue);
  hSig_CR[0] -> SetTitle("TT Contamination tTagger '16");
  hSig_CR[1] = (TH1F*)signal[1]->Get(TString::Format("CR_tTagger_%s_expYield",recoVar.Data()));
  hSig_CR[1] -> SetLineColor(kBlue);
  hSig_CR[1] -> SetTitle("TT Contamination tTagger '17");
  hSig_CR[2] = (TH1F*)signal[2]->Get(TString::Format("CR_tTagger_%s_expYield",recoVar.Data()));
  hSig_CR[2] -> SetLineColor(kBlue);
  hSig_CR[2] -> SetTitle("TT Contamination tTagger '18");
  
  
  TLegend *closureLegend = new TLegend(0.5,0.6,0.7,0.8);
  closureLegend->AddEntry(hBkg_SR[0],"Signal Region (2btag)", "l");
  closureLegend->AddEntry(hBkg_CR[0],"Control Region (0btag)", "l");
  closureLegend->AddEntry(hBkg_1Btag[0],"1btag Region", "l");
   
  //----------------------------------FOR 2018---------------------------------------------------------------------
  auto c18 = new TCanvas("QCD closure Test tTagger '18", "QCD closure Test tTagger '18", 700,600);
  auto *closure_pad182 = new TPad("closure_pad18","closure_pad18",0.,0.,1.,0.3); 
  closure_pad182->Draw();
  closure_pad182->SetTopMargin(0.05);
  closure_pad182->SetBottomMargin(0.2);
  closure_pad182->SetGrid();

  auto *closure_pad181 = new TPad("closure_pad182","closure_pad182",0.,0.3,1.,1.);  
  closure_pad181->Draw();
  closure_pad181->SetBottomMargin(0.001);
  closure_pad181->cd();
  //closure_pad1->SetGrid();
  hBkg_SR[2]->GetYaxis()->SetTitleSize(20);
  hBkg_SR[2]->GetYaxis()->SetTitleFont(43);
  hBkg_SR[2]->GetYaxis()->SetTitleOffset(1.55);  
  // h2 settings
  hBkg_SR[2]->Draw();
  hBkg_CR[2]->Draw("same");
  hBkg_1Btag[2]->Draw("same");
  closureLegend->Draw();

  TH1F *hClosure18[2];
  closure_pad182->cd();
  hClosure18[0] = (TH1F*)hBkg_SR[2]->Clone("hClosure18_0"); 
  hClosure18[1] = (TH1F*)hBkg_1Btag[2]->Clone("hClosure18_1"); 
  hClosure18[0]->Divide(hBkg_CR[2]);
  hClosure18[1]->Divide(hBkg_CR[2]);
  hClosure18[0]->Draw();
  hClosure18[0]->SetLineColor(kRed);
  hClosure18[1]->SetLineColor(kBlack);
  hClosure18[1]->Draw("same");  
   
  hClosure18[0]->SetTitle(""); 
  hClosure18[0]->GetYaxis()->SetTitle("ratio SR/CR");
  hClosure18[0]->GetYaxis()->SetTitleSize(14);
  hClosure18[0]->GetYaxis()->SetTitleFont(43);
  hClosure18[0]->GetYaxis()->SetTitleOffset(1.55);
  hClosure18[0]->GetYaxis()->SetLabelFont(43);
  hClosure18[0]->GetYaxis()->SetLabelSize(15);
  
  hClosure18[0]->GetXaxis()->SetTitleSize(15);
  hClosure18[0]->GetXaxis()->SetTitleFont(43);
  hClosure18[0]->GetXaxis()->SetTitleOffset(3.2);
  hClosure18[0]->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  hClosure18[0]->GetXaxis()->SetLabelSize(15); 
  //----------------------------------END 2018---------------------------------------------------------------------
   
   
  //----------------------------------FOR 2017---------------------------------------------------------------------   
  auto c1 = new TCanvas("QCD closure Test tTagger '17", "QCD closure Test tTagger '17", 700,600);
  auto *closure_pad2 = new TPad("closure_pad2","closure_pad2",0.,0.,1.,0.3); 
  closure_pad2->Draw();
  closure_pad2->SetTopMargin(0.05);
  closure_pad2->SetBottomMargin(0.2);
  closure_pad2->SetGrid();

  auto *closure_pad1 = new TPad("closure_pad1","closure_pad1",0.,0.3,1.,1.);  
  closure_pad1->Draw();
  closure_pad1->SetBottomMargin(0.001);
  closure_pad1->cd();
  //closure_pad1->SetGrid();
  hBkg_SR[1]->GetYaxis()->SetTitleSize(20);
  hBkg_SR[1]->GetYaxis()->SetTitleFont(43);
  hBkg_SR[1]->GetYaxis()->SetTitleOffset(1.55);  
  // h2 settings
  hBkg_SR[1]->Draw();
  hBkg_CR[1]->Draw("same");
  hBkg_1Btag[1]->Draw("same");
  closureLegend->Draw();

  TH1F *hClosure17[2];
  closure_pad2->cd();
  hClosure17[0] = (TH1F*)hBkg_SR[1]->Clone("hClosure17_0"); 
  hClosure17[1] = (TH1F*)hBkg_1Btag[1]->Clone("hClosure17_1"); 
  hClosure17[0]->Divide(hBkg_CR[1]);
  hClosure17[1]->Divide(hBkg_CR[1]);
  hClosure17[0]->Draw();
  hClosure17[0]->SetLineColor(kRed);
  hClosure17[1]->SetLineColor(kBlack);
  hClosure17[1]->Draw("same");  
   
  hClosure17[0]->SetTitle(""); 
  hClosure17[0]->GetYaxis()->SetTitle("ratio SR/CR");
  hClosure17[0]->GetYaxis()->SetTitleSize(14);
  hClosure17[0]->GetYaxis()->SetTitleFont(43);
  hClosure17[0]->GetYaxis()->SetTitleOffset(1.55);
  hClosure17[0]->GetYaxis()->SetLabelFont(43);
  hClosure17[0]->GetYaxis()->SetLabelSize(15);
  
  hClosure17[0]->GetXaxis()->SetTitleSize(15);
  hClosure17[0]->GetXaxis()->SetTitleFont(43);
  hClosure17[0]->GetXaxis()->SetTitleOffset(3.2);
  hClosure17[0]->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  hClosure17[0]->GetXaxis()->SetLabelSize(15); 
  //----------------------------------END 2017--------------------------------------------------------------------- 
  
  
  //----------------------------------FOR 2016---------------------------------------------------------------------
  auto c16 = new TCanvas("QCD closure Test tTagger '16", "QCD closure Test tTagger '16", 700,600);
  auto *closure_pad2_16 = new TPad("closure_pad2_16","closure_pad2_16",0.,0.,1.,0.3); 
  closure_pad2_16->Draw();
  closure_pad2_16->SetTopMargin(0.05);
  closure_pad2_16->SetBottomMargin(0.2);
  closure_pad2_16->SetGrid();

  auto *closure_pad1_16 = new TPad("closure_pad1_16","closure_pad1_16",0.,0.3,1.,1.);  
  closure_pad1_16->Draw();
  closure_pad1_16->SetBottomMargin(0.001);
  closure_pad1_16->cd();
  //closure_pad1->SetGrid();
  hBkg_SR[0]->GetYaxis()->SetTitleSize(20);
  hBkg_SR[0]->GetYaxis()->SetTitleFont(43);
  hBkg_SR[0]->GetYaxis()->SetTitleOffset(1.55);  
  // h2 settings
  hBkg_SR[0]->Draw();
  hBkg_CR[0]->Draw("same");
  hBkg_1Btag[0]->Draw("same");
  closureLegend->Draw();

  TH1F *hClosure16[2];
  closure_pad2_16->cd();
  hClosure16[0] = (TH1F*)hBkg_SR[0]->Clone("hClosure16_0"); 
  hClosure16[1] = (TH1F*)hBkg_1Btag[0]->Clone("hClosure16_1"); 
  hClosure16[0]->Divide(hBkg_CR[0]);
  hClosure16[1]->Divide(hBkg_CR[0]);
  hClosure16[0]->Draw();
  hClosure16[0]->SetLineColor(kRed);
  hClosure16[1]->SetLineColor(kBlack);
  hClosure16[1]->Draw("same"); 
  
  hClosure16[0]->SetTitle("");
  hClosure16[0]->GetYaxis()->SetTitle("ratio SR/CR");
  hClosure16[0]->GetYaxis()->SetTitleSize(14);
  hClosure16[0]->GetYaxis()->SetTitleFont(43);
  hClosure16[0]->GetYaxis()->SetTitleOffset(1.55);
  hClosure16[0]->GetYaxis()->SetLabelFont(43);
  hClosure16[0]->GetYaxis()->SetLabelSize(15);
  
  hClosure16[0]->GetXaxis()->SetTitleSize(15);
  hClosure16[0]->GetXaxis()->SetTitleFont(43);
  hClosure16[0]->GetXaxis()->SetTitleOffset(3.2);
  hClosure16[0]->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  hClosure16[0]->GetXaxis()->SetLabelSize(15);
  //----------------------------------END 2016--------------------------------------------------------------------- 
   

  //TT contamination for 2016
  auto cTT16 = new TCanvas("TT contamination tTagger '16", "TT contamination tTagger '16", 700,600);
  auto rp_tTagger16 = new TRatioPlot(hSig_CR[0],hBkg_CRExpYield[0]);
  cTT16->SetTicks(0,1);
  rp_tTagger16->Draw();
  rp_tTagger16->GetUpperRefYaxis()->SetTitle("Exp. Yield");
  rp_tTagger16->GetLowerRefYaxis()->SetTitle("ratio");
  cTT16->Update();  
  
  //TT contamination for 2017 
  auto cTT17 = new TCanvas("TT contamination tTagger '17", "TT contamination tTagger '17", 700,600);
  auto rp_tTagger17 = new TRatioPlot(hSig_CR[1],hBkg_CRExpYield[1]);
  cTT17->SetTicks(0,1);
  rp_tTagger17->Draw();
  rp_tTagger17->GetUpperRefYaxis()->SetTitle("Exp. Yield");
  rp_tTagger17->GetLowerRefYaxis()->SetTitle("ratio");
  cTT17->Update();
  
  //TT contamination for 2018
  auto  cTT18 = new TCanvas("TT contamination tTagger '18", "TT contamination tTagger '18", 700,600);
  auto rp_tTagger18 = new TRatioPlot(hSig_CR[2],hBkg_CRExpYield[2]);
  cTT18->SetTicks(0,1);
  rp_tTagger18->Draw();
  rp_tTagger18->GetUpperRefYaxis()->SetTitle("Exp. Yield");
  rp_tTagger18->GetLowerRefYaxis()->SetTitle("ratio");
  cTT18->Update();
  
  
  
  
  
  
}
