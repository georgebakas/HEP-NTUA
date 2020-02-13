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
void ratioPlot(TH1F *hNum ,TH1F *hDenom, TString recoVar, TString reason, TString level);
bool globalIsNominalMC;
TString globalYear1;

void plotTT_yearly(TString year = "2016")
{
	globalYear1 = year;
	globalYear1.Remove(TString::kBoth,'2');
	globalYear1.Remove(TString::kBoth,'0');

	const int NVAR = 9;
	TString recoVar[NVAR] = {"jetPt0", "mJJ", "ptJJ", "yJJ", "jetPt1","jetY0", "jetY1"
							 ,"mTop", "jetMassSoftDrop"};

	for(int ivar = 0; ivar<NVAR; ivar++)
	{
		plotYearVar(year,recoVar[ivar]);
		//break;
	}
}
void plotYearVar(TString year, TString recoVar = "jetPt0")
{
  initFilesMapping();

  std::vector<TString> histoNamesTTMtt;
  std::vector<TString> histoNamesTTNominal;
  histoNamesTTMtt.clear();
  histoNamesTTNominal.clear();

  if(!globalIsNominalMC)
  {
  	histoNamesTTMtt.push_back("Signal_histo_Mtt_700_1000"); 
  	histoNamesTTMtt.push_back("Signal_histo_Mtt_1000_Inf");

  	if(!year.EqualTo("2016"))
  	{
  		histoNamesTTNominal.push_back("Signal_histo_TTHadronic");
  		histoNamesTTNominal.push_back("Signal_histo_TTSemiLeptonic");
  		histoNamesTTNominal.push_back("Signal_histo_TTTo2L2Nu");
  	}
  	else 
  		histoNamesTTNominal.push_back("Signal_histo_Nominal");
  }
  
  //0 is 2017 and //1 is 2018
  //TT files
  TFile *infTT[2];
  infTT[0] = TFile::Open(filesttbar[TString::Format("%s",year.Data())]); 
  infTT[1] = TFile::Open(filesttbar[TString::Format("%sNominal",year.Data())]);  

  TH1F *hSigAll_Nominal, *hSigAll_Mtt;
  

  /*
  all
  hScaledXSEC_mJJ
  hPartonScaledXSEC_mJJ
  hParticleScaledXSEC_mJJ

  per slice
  h_Signal_histo_Mtt_700_1000_mJJ
  hParton_Signal_histo_Mtt_700_1000_mJJ
  hParticle_Signal_histo_Mtt_700_1000_mJJ
  */

  //get the histograms 
  //tt slices
  //TString level = ""; //reco 
  //TString level = "Parton";
  //TString level = "Particle";
 
  //scaled TT to xsec and LUMI
  hSigAll_Mtt = (TH1F*)infTT[0]->Get(TString::Format("h%sScaledXSEC_%s",level.Data(), recoVar.Data()));
  hSigAll_Nominal = (TH1F*)infTT[1]->Get(TString::Format("h%sScaledXSEC_%s",level.Data(), recoVar.Data()));
  hSigAll_Mtt->SetTitle(TString::Format("TT %s %s Mtt",level.Data(), year.Data()));
  hSigAll_Nominal->SetTitle(TString::Format("TT %s %s Nominal",level.Data(), year.Data()));
  hSigAll_Mtt->SetName(TString::Format("TT %s %s Mtt",level.Data(), year.Data()));
  hSigAll_Nominal->SetName(TString::Format("TT %s %s Nominal",level.Data(), year.Data()));
  
  TString reason;  
  reason = "TT MC Overlay all slices";
  ratioPlot(hSigAll_Mtt, hSigAll_Nominal,recoVar,reason, level);
  
  histoNamesTTNominal.clear();
  histoNamesTTMtt.clear();
}



void ratioPlot(TH1F *hNum ,TH1F *hDenom, TString recoVar, TString reason, TString level)
{
  TString titleNum = hNum->GetTitle();
  TString titleDenom = hDenom->GetTitle();
  TLegend *closureLegend = new TLegend(0.65,0.73,0.9,0.9);
  closureLegend->AddEntry(hNum,titleNum, "lep");
  closureLegend->AddEntry(hDenom,titleDenom, "lep");

  hNum->Scale(1./hNum->Integral());
  hDenom->Scale(1./hDenom->Integral());
  //cout<<"---------------"<<endl;
  //cout<<reason<<endl;
  //cout<<"hNum name: "<<hNum->GetTitle()<<endl;
  //cout<<"hDenom name: "<<hDenom->GetTitle()<<endl;

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
  if(!recoVar.EqualTo("yJJ") || !recoVar.EqualTo("jetY0") || !recoVar.EqualTo("jetY1")) hNum->GetXaxis()->SetTitle(recoVar+" (GeV)");
  else hNum->GetXaxis()->SetTitle(recoVar);
  // h2 settings
  
  //hBkg_CR[0]->ResetAttFill();
  hDenom->ResetAttLine();
  hDenom->ResetAttMarker();
  hDenom->SetLineColor(kBlue-2); 
  hDenom->SetMarkerStyle(21);
  hDenom->SetMarkerColor(kBlue-2);
  
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


  //hRatio->GetYaxis()->SetRangeUser(0,2);
  hRatio->GetXaxis()->SetLabelSize(0.06);
  hRatio->Divide(hDenom);
  hRatio->SetLineColor(kBlack);
  if(!recoVar.EqualTo("yJJ") || !recoVar.EqualTo("jetY0") || !recoVar.EqualTo("jetY1")) hRatio->GetXaxis()->SetTitle(recoVar+" (GeV)");
  else hRatio->GetXaxis()->SetTitle(recoVar);
  //hRatio->GetXaxis()->SetTitleOffset(1);
  hRatio->Draw();
  if(level.EqualTo("")) level = "Reco";
  c1->Print(TString::Format("./YearlyComparison/%s/%s/%s/comparison_mc_%s.pdf",globalYear1.Data(),
  							level.Data(),recoVar.Data(),reason.Data()),"pdf");
}