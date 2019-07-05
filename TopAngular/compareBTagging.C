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

std::vector<TString> listOfFiles;
std::vector<float> XSEC;
std::vector<TString> histoNames;
float LUMI = 41530;
TString eosPath;

void initFileNames()
{
	
  //listOfFiles.push_back("/eos/cms/store/user/ipapakri/ttbar/MC/Signal/2016/TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8_legacy2016_deepAK8.root");
  listOfFiles.push_back("/eos/cms/store/user/ipapakri/ttbar/MC/Signal/2017/TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root");
  listOfFiles.push_back("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/2017/QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8.root");
  listOfFiles.push_back("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/2017/QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8.root");
  listOfFiles.push_back("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/2017/QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8.root");
  listOfFiles.push_back("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/2017/QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8.root");
  listOfFiles.push_back("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/2017/QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8.root");
  listOfFiles.push_back("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/2017/QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8.root");
}

void initXsections()
{
  XSEC.push_back(16.74); //this is the TT xsec
  XSEC.push_back(3.67e+5);
  XSEC.push_back(2.94e+4);
  XSEC.push_back(6.524e+03);
  XSEC.push_back(1.064e+03);
  XSEC.push_back(121.5);
  XSEC.push_back(2.542e+01);
}


void initGlobals()
{
  initFileNames();
  initXsections();
}
void compareBTagging(TString recoVar = "chi",TString closureVar = "Chi", bool is2016 = false)
{
  initGlobals();
  std::vector<float> weights(0);
  for(int f=0; f<listOfFiles.size(); f++)
  //[0] is the mtt and from [1] up to [6] its the Bkg
  //the histograms from QCD are already scaled with LUMI
  {
	TFile *file =TFile::Open(listOfFiles[f]);
	float norm = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
    float weight = XSEC[f]/norm;
    weights.push_back(weight);  
  }
   TFile *eff[2], *signal[2], *bkg[2];
   
   //check what year you want to compare 
   //if compare2016 is true then you compare just 2016
   std::vector<Color_t> colors = {kBlue,kRed, kBlack,  kMagenta};
   if(is2016)
   {
	   //for closure tests [0] is CSVv2  and [1] is deepCSV
	   signal[0] = TFile::Open("./Code_2016/SignalOutput_AllRegions_0.2.root");
	   signal[1] = TFile::Open("./Code_2016/SignalOutput_AllRegions_0.2.root"); //NEEDS CHANGE

	   bkg[1] = TFile::Open("./Code_2016/BkgOutput_AllRegions_0.2.root");
	   bkg[1] = TFile::Open("./Code_2016/BkgOutput_AllRegions_0.2.root"); //NEEDS CHAGE
	   
 	   //for efficiencies
       eff[1] = TFile::Open("./Code_2016/deepAK8_efficiencies_allVars_tTagger_0.2.root");
       eff[1] = TFile::Open("./Code_2016/deepAK8_efficiencies_allVars_tTagger_0.2.root"); //NEEDS CHANGE
   }
   else
   {
	   //for closure tests [0] is CSVv2  and [1] is deepCSV
	   signal[0] = TFile::Open("./Code_2017/Output_TT_QCD_Reco_Chi_0.1.root");
	   signal[1] = TFile::Open("./Code_2017/Output_TT_QCD_Reco_Chi_0.1_deepCSV.root");

	   bkg[0] = TFile::Open("./Code_2017/Closure_QCDBkg_Chi_0.1.root");
	   bkg[1] = TFile::Open("./Code_2017/Closure_QCDBkg_Chi_0.1_deepCSV.root");
	   
 	   //for efficiencies/ acceptance purity and stability
	   eff[0] = TFile::Open("./Code_2017/ResponseMatricesChiCos_0.1.root");
	   eff[1] = TFile::Open("./Code_2017/ResponseMatricesChiCos_0.1_deepCSV.root");

   }
   
   //things we need:
   //1. Efficiency for 2017 and acceptance for same year both for tTagger and eventTagger
   TEfficiency *effCSVv2[2], *accCSVv2[2];
   effCSVv2[0] = (TEfficiency*)eff[0]->Get(TString::Format("%sEfficiency", recoVar.Data()));
   accCSVv2[0] = (TEfficiency*)eff[0]->Get(TString::Format("%sAcceptance", recoVar.Data()));
   
   //2. Efficiency for 2016 and acceptance for same year
   TEfficiency *effDeepCSV[2], *accDeepCSV[2];
   effDeepCSV[0] = (TEfficiency*)eff[1]->Get(TString::Format("%sEfficiency", recoVar.Data()));
   accDeepCSV[0] = (TEfficiency*)eff[1]->Get(TString::Format("%sAcceptance", recoVar.Data()));

   
  
   //efficiency coloring
   effCSVv2[0]->SetMarkerStyle(21);
   effDeepCSV[0]->SetMarkerStyle(22);
   effCSVv2[0]->SetMarkerColor(colors[0]);
   effDeepCSV[0]->SetMarkerColor(colors[1]);
   effCSVv2[0]->SetLineColor(colors[0]);
   effDeepCSV[0]->SetLineColor(colors[1]);

   //acceptance coloring
   accCSVv2[0]->SetMarkerStyle(21);
   accDeepCSV[0]->SetMarkerStyle(22);
   accCSVv2[0]->SetLineColor(colors[0]);
   accDeepCSV[0]->SetLineColor(colors[1]);
   accCSVv2[0]->SetMarkerColor(colors[0]);
   accDeepCSV[0]->SetMarkerColor(colors[1]);

   
   TLegend *effLeg = new TLegend(0.5,0.6,0.7,0.8);
   effLeg->AddEntry(effCSVv2[0], "tTagger CSVv2", "lp");
   effLeg->AddEntry(effDeepCSV[0], "tTagger deepCSV", "lp");
   TCanvas *can_eff = new TCanvas("Efficiency can", "Efficiency can", 700, 600);
   effCSVv2[0]->SetTitle(TString::Format("Efficiency CSVv2, DeepCSV;%s (GeV);Efficiency",recoVar.Data())); 
   effCSVv2[0]->Draw();
   effDeepCSV[0]->Draw("same");
   effLeg->Draw();
   
  
   TCanvas *can_acc = new TCanvas("Acceptance can", "Acceptance can", 700, 600);
   accCSVv2[0]->Draw();
   accCSVv2[0]->SetTitle(TString::Format("Acceptance CSVv2, DeepCSV;%s (GeV);Acceptance",recoVar.Data())); 
   accDeepCSV[0]->Draw("same");
   effLeg->Draw();
   
  
  //-----------------------------------------------------------------------------------------------------------------
  //now we compare the closure tests for the CSVv2 and deepCSV  
  TH1F *hBkg_CR[2], *hBkg_SR[2], *hBkg_CRExpYield[2];  
  TH1F *hSig_CR[2];    
  TH1F *hSig_1Btag[2], *hBkg_1Btag[2];
   //TH1F used for the TRatioPlot [0] is CSVv2 and [1] is DeepCSV
  hBkg_CR[0] = (TH1F*)bkg[0]->Get(TString::Format("h%s_QCD_CR",closureVar.Data()));
  hBkg_CR[1] = (TH1F*)bkg[1]->Get(TString::Format("h%s_QCD_CR",closureVar.Data()));
  hBkg_SR[0] = (TH1F*)bkg[0]->Get(TString::Format("h%s_QCD_SR",closureVar.Data()));
  hBkg_SR[1] = (TH1F*)bkg[1]->Get(TString::Format("h%s_QCD_SR",closureVar.Data()));
  hBkg_1Btag[0] = (TH1F*)bkg[0]->Get(TString::Format("h%s_QCD_1btag",closureVar.Data()));
  hBkg_1Btag[1] = (TH1F*)bkg[1]->Get(TString::Format("h%s_QCD_1btag",closureVar.Data()));  

  hBkg_CRExpYield[0] = (TH1F*)hBkg_CR[0]->Clone(TString::Format("CR_tTagger_%s_expYield_0",closureVar.Data()));
  hBkg_CRExpYield[1] = (TH1F*)hBkg_CR[1]->Clone(TString::Format("CR_tTagger_%s_expYield_1",closureVar.Data()));
  hBkg_CRExpYield[0]->SetLineColor(kRed);
  hBkg_CRExpYield[1]->SetLineColor(kRed);
  cout<<"HERE"<<endl;
  //now you have to scale these to 1.
  hBkg_CR[0] ->Scale(1./hBkg_CR[0]->Integral());
  hBkg_CR[1] ->Scale(1./hBkg_CR[1]->Integral());
  hBkg_SR[0] ->Scale(1./hBkg_SR[0]->Integral());
  hBkg_SR[1] ->Scale(1./hBkg_SR[1]->Integral());
  hBkg_1Btag[0] ->Scale(1./hBkg_1Btag[0]->Integral());
  hBkg_1Btag[1] ->Scale(1./hBkg_1Btag[1]->Integral());
  
  hBkg_CRExpYield[0]->SetTitle("TT Contamination tTagger CSVv2");
  hBkg_CRExpYield[1]->SetTitle("TT Contamination tTagger deepCSV");
  
  hSig_CR[0] = (TH1F*)signal[0]->Get(TString::Format("h%sCR_TT",closureVar.Data()));
  hSig_CR[0] -> SetLineColor(kBlue);
  hSig_CR[0] -> Scale(weights[0]*LUMI);
  hSig_CR[0] -> SetTitle("TT Contamination tTagger CSVv2");
  hSig_CR[1] = (TH1F*)signal[1]->Get(TString::Format("h%sCR_TT",closureVar.Data()));
  hSig_CR[1] -> SetLineColor(kBlue);
  hSig_CR[1] -> Scale(weights[0]*LUMI);
  hSig_CR[1] -> SetTitle("TT Contamination tTagger deepCSV");
  
  
  TLegend *closureLegend = new TLegend(0.5,0.6,0.7,0.8);
  closureLegend->AddEntry(hBkg_SR[0],"Signal Region (2btag)", "l");
  closureLegend->AddEntry(hBkg_CR[0],"Control Region (0btag)", "l");
  closureLegend->AddEntry(hBkg_1Btag[0],"1btag Region", "l");
   
   
   //FOR CSVv2
  auto c1 = new TCanvas("QCD closure Test tTagger CSVv2", "QCD closure Test tTagger CSVv2", 700,600);
  

  auto *closure_pad1 = new TPad("closure_pad1","closure_pad1",0.,0.3,1.,1.);  
  closure_pad1->Draw();
  closure_pad1->SetBottomMargin(0.001);
  closure_pad1->cd();
  //closure_pad1->SetGrid();
  hBkg_SR[0]->GetYaxis()->SetTitleSize(20);
  hBkg_SR[0]->GetYaxis()->SetTitleFont(43);
  hBkg_SR[0]->GetYaxis()->SetTitleOffset(1.55);  
  // h2 settings
  hBkg_SR[0]->Draw();
  hBkg_CR[0]->Draw("same");
  hBkg_1Btag[0]->Draw("same");
  closureLegend->Draw();
  
  c1->cd();
  auto *closure_pad2 = new TPad("closure_pad2","closure_pad2",0.,0.,1.,0.3); 
  closure_pad2->Draw();
  closure_pad2->SetTopMargin(0.05);
  closure_pad2->SetBottomMargin(0.2);
  closure_pad2->SetGrid();

  TH1F *hClosureCSVv2[2];
  closure_pad2->cd();
  hClosureCSVv2[0] = (TH1F*)hBkg_SR[0]->Clone("hClosureCSVv2_0"); 
  hClosureCSVv2[0]->SetTitle("");
  
  hClosureCSVv2[1] = (TH1F*)hBkg_1Btag[0]->Clone("hClosureCSVv2_1"); 
  hClosureCSVv2[0]->Divide(hBkg_CR[0]);
  hClosureCSVv2[1]->Divide(hBkg_CR[0]);
  hClosureCSVv2[0]->Draw();
  hClosureCSVv2[0]->SetLineColor(kRed);
  hClosureCSVv2[1]->SetLineColor(kBlack);
  hClosureCSVv2[1]->Draw("same");  
   
  hClosureCSVv2[0]->GetYaxis()->SetTitle("ratio SR/CR");
  hClosureCSVv2[0]->GetYaxis()->SetTitleSize(14);
  hClosureCSVv2[0]->GetYaxis()->SetTitleFont(43);
  hClosureCSVv2[0]->GetYaxis()->SetTitleOffset(1.55);
  hClosureCSVv2[0]->GetYaxis()->SetLabelFont(43);
  hClosureCSVv2[0]->GetYaxis()->SetLabelSize(15);
  
  hClosureCSVv2[0]->GetXaxis()->SetTitleSize(15);
  hClosureCSVv2[0]->GetXaxis()->SetTitleFont(43);
  hClosureCSVv2[0]->GetXaxis()->SetTitleOffset(3.2);
  hClosureCSVv2[0]->GetXaxis()->SetLabelSize(15);
  hClosureCSVv2[0]->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   
  
  
  
  //FOR DeepCSV
  auto c16 = new TCanvas("QCD closure Test tTagger deepCSV", "QCD closure Test tTagger deepCSV", 700,600);
 

  auto *closure_pad1_16 = new TPad("closure_pad1_deepCSV","closure_pad1_deepCSV",0.,0.3,1.,1.);  
  closure_pad1_16->Draw();
  closure_pad1_16->SetBottomMargin(0.001);
  closure_pad1_16->cd();
  //closure_pad1->SetGrid();
  hBkg_SR[1]->GetYaxis()->SetTitleSize(20);
  hBkg_SR[1]->GetYaxis()->SetTitleFont(43);
  hBkg_SR[1]->GetYaxis()->SetTitleOffset(1.55);  
  // h2 settings
  hBkg_SR[1]->Draw();
  hBkg_CR[1]->Draw("same");
  hBkg_1Btag[1]->Draw("same");
  closureLegend->Draw();

  c16->cd();
   auto *closure_pad2_16 = new TPad("closure_pad2_deepCSV","closure_pad2_deepCSV",0.,0.,1.,0.28	); 
  closure_pad2_16->Draw();
  closure_pad2_16->SetTopMargin(0.05);
  closure_pad2_16->SetBottomMargin(0.2);
  closure_pad2_16->SetGrid();
  
  TH1F *hClosureDeepCSV[2];
  closure_pad2_16->cd();
  hClosureDeepCSV[0] = (TH1F*)hBkg_SR[1]->Clone("hClosureDeepCSV_0"); 
  hClosureDeepCSV[0] = (TH1F*)hBkg_SR[1]->Clone("hClosureDeepCSV_0"); 
  hClosureDeepCSV[1] = (TH1F*)hBkg_1Btag[1]->Clone("hClosureDeepCSV_1"); 
  hClosureDeepCSV[0]->Divide(hBkg_CR[0]);
  hClosureDeepCSV[1]->Divide(hBkg_CR[0]);
 
  hClosureDeepCSV[0]->Draw();
  hClosureDeepCSV[0]->SetLineColor(kRed);
  hClosureDeepCSV[1]->SetLineColor(kBlack);
  hClosureDeepCSV[1]->Draw("same"); 
  
  hClosureDeepCSV[0]->SetTitle("");
  
  hClosureDeepCSV[0]->GetYaxis()->SetTitle("ratio SR/CR");
  hClosureDeepCSV[0]->GetYaxis()->SetTitleSize(14);
  hClosureDeepCSV[0]->GetYaxis()->SetTitleFont(43);
  hClosureDeepCSV[0]->GetYaxis()->SetTitleOffset(1.55);
  hClosureDeepCSV[0]->GetYaxis()->SetLabelFont(43);
  hClosureDeepCSV[0]->GetYaxis()->SetLabelSize(15);
  
  hClosureDeepCSV[0]->GetXaxis()->SetTitleSize(15);
  hClosureDeepCSV[0]->GetXaxis()->SetTitleFont(43);
  hClosureDeepCSV[0]->GetXaxis()->SetTitleOffset(3.2);
  hClosureDeepCSV[0]->GetXaxis()->SetLabelSize(15);
  hClosureDeepCSV[0]->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  
  
  auto c3 = new TCanvas("TT contamination tTagger CSVv2", "TT contamination tTagger CSVv2", 700,600);
  auto rp_tTaggerCSVv2 = new TRatioPlot(hSig_CR[0],hBkg_CRExpYield[0]);
  c3->SetTicks(0,1);
  rp_tTaggerCSVv2->Draw();
  rp_tTaggerCSVv2->GetUpperRefYaxis()->SetTitle("Exp. Yield");
  rp_tTaggerCSVv2->GetLowerRefYaxis()->SetTitle("ratio");
  c3->Update();
  
  
  auto c4 = new TCanvas("TT contamination tTagger deepCSV", "TT contamination tTagger deepCSV", 700,600);
  auto rp_tTaggerDeepCSV = new TRatioPlot(hSig_CR[1],hBkg_CRExpYield[1]);
  c4->SetTicks(0,1);
  rp_tTaggerDeepCSV->Draw();
  rp_tTaggerDeepCSV->GetUpperRefYaxis()->SetTitle("Exp. Yield");
  rp_tTaggerDeepCSV->GetLowerRefYaxis()->SetTitle("ratio");
  c3->Update();
	 
}