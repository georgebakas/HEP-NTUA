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
void plotYearVar(TString year, TString year2, TString recoVar = "jetPt0");
void ratioPlot(TH1F *hNum ,TH1F *hDenom, TString recoVar, TString reason);
bool globalIsNominalMC;
TString globalYear1, globalYear2;

void plotTT_QCD(TString year1 = "2017", TString year2 = "2018", bool isNominalMC = false)
{
	globalIsNominalMC = isNominalMC;
	globalYear1 = year1;
	globalYear2 = year2;
	globalYear1.Remove(TString::kBoth,'2');
	globalYear2.Remove(TString::kBoth,'2');
	globalYear1.Remove(TString::kBoth,'0');
	globalYear2.Remove(TString::kBoth,'0');

	const int NVAR = 11;
	TString recoVar[NVAR]   = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1","jetPhi0","jetPhi1","mTop0", "mTop1"};

	for(int ivar = 0; ivar<NVAR; ivar++)
	{
		plotYearVar(year1, year2,recoVar[ivar]);
		//break;
	}
}
void plotYearVar(TString year1, TString year2 ,TString recoVar = "jetPt0")
{

  initFilesMapping();

  std::vector<TString> histoNamesTT;
  std::vector<TString> histoNamesQCD;
  histoNamesTT.clear();
  histoNamesQCD.clear();

  if(!globalIsNominalMC)
  {
  	histoNamesTT.push_back("Signal_histo_Mtt_700_1000");
  	histoNamesTT.push_back("Signal_histo_Mtt_1000_Inf");
  }
  else
  {
  	histoNamesTT.push_back("Signal_histo_TTHadronic");
  	histoNamesTT.push_back("Signal_histo_TTSemiLeptonic");
  	histoNamesTT.push_back("Signal_histo_TTTo2L2Nu");
  }
  histoNamesQCD.push_back("QCD_histo_300_500");
  histoNamesQCD.push_back("QCD_histo_500_700");
  histoNamesQCD.push_back("QCD_histo_700_1000");
  histoNamesQCD.push_back("QCD_histo_1000_1500");
  histoNamesQCD.push_back("QCD_histo_1500_2000");
  histoNamesQCD.push_back("QCD_histo_2000_Inf");

  //0 is 2017 and //1 is 2018
  //TT files
  TFile *infTT[2];
  cout<<histoNamesTT.size()<<endl;

  TString temp = "";
  if(globalIsNominalMC) temp = "Nominal";
  infTT[0] = TFile::Open(filesttbar[TString::Format("%s%s",year1.Data(),temp.Data())]);
  infTT[1] = TFile::Open(filesttbar[TString::Format("%s%s",year2.Data(),temp.Data())]);

  //QCD files
  TFile *infBkg[2];
  infBkg[0] = TFile::Open(filesqcd[year1.Data()]);
  infBkg[1] = TFile::Open(filesqcd[year2.Data()]);

  TH1F *hQCDAll[2], *hSigAll[2];

  TH1F *hQCD_slice[2][6], *hSig_slice[2][histoNamesTT.size()];


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
  TString level = "Parton";
  //TString level = "Particle";
  for(int i =0; i<histoNamesTT.size(); i++)
  {
		cout<<TString::Format("h%s_%s_%s",level.Data(),histoNamesTT[i].Data(),recoVar.Data())<<endl;
  	hSig_slice[0][i] = (TH1F*)infTT[0]->Get(TString::Format("h%s_%s_%s",level.Data(),histoNamesTT[i].Data(),recoVar.Data()));
  	hSig_slice[1][i] = (TH1F*)infTT[1]->Get(TString::Format("h%s_%s_%s",level.Data(),histoNamesTT[i].Data(),recoVar.Data()));
  	hSig_slice[0][i]->SetTitle(TString::Format("Slice %s: %s %s",level.Data(), histoNamesTT[i].Data(), year1.Data()));
  	hSig_slice[1][i]->SetTitle(TString::Format("Slice %s: %s %s",level.Data(), histoNamesTT[i].Data(), year2.Data()));
  	hSig_slice[0][i]->SetName(TString::Format("Slice %s: %s %s",level.Data(), histoNamesTT[i].Data(), year1.Data()));
  	hSig_slice[1][i]->SetName(TString::Format("Slice %s: %s %s",level.Data(), histoNamesTT[i].Data(), year2.Data()));
  }
  /*
  //qcd slices
  for(int i =0; i<histoNamesQCD.size(); i++)
  {
  	hQCD_slice[0][i] = (TH1F*)infBkg[0]->Get(TString::Format("h_%s_%s",histoNamesQCD[i].Data(),recoVar.Data()));
  	hQCD_slice[1][i] = (TH1F*)infBkg[1]->Get(TString::Format("h_%s_%s",histoNamesQCD[i].Data(),recoVar.Data()));

  	hQCD_slice[0][i]->SetTitle(TString::Format("Slice: %s %s", histoNamesQCD[i].Data(), year1.Data()));
  	hQCD_slice[1][i]->SetTitle(TString::Format("Slice: %s %s", histoNamesQCD[i].Data(), year2.Data()));
  	hQCD_slice[0][i]->SetName(TString::Format("Slice: %s %s", histoNamesQCD[i].Data(), year1.Data()));
  	hQCD_slice[1][i]->SetName(TString::Format("Slice: %s %s", histoNamesQCD[i].Data(), year2.Data()));
  }*/
  //scaled TT to xsec and LUMI
  hSigAll[0] = (TH1F*)infTT[0]->Get(TString::Format("h%sScaledXSEC_%s",level.Data(), recoVar.Data()));
  hSigAll[1] = (TH1F*)infTT[1]->Get(TString::Format("h%sScaledXSEC_%s",level.Data(), recoVar.Data()));
  hSigAll[0]->SetTitle(TString::Format("TT %s All Slices %s",level.Data(), year1.Data()));
  hSigAll[1]->SetTitle(TString::Format("TT %s All Slices %s",level.Data(), year2.Data()));
  hSigAll[0]->SetName(TString::Format("TT %s All Slices %s",level.Data(), year1.Data()));
  hSigAll[1]->SetName(TString::Format("TT %s All Slices %s",level.Data(), year2.Data()));

  //scaled qcd to xsec and LUMI
  /*
  hQCDAll[0] = (TH1F*)infBkg[0]->Get(TString::Format("hScaledXSEC_%s", recoVar.Data()));
  hQCDAll[1] = (TH1F*)infBkg[1]->Get(TString::Format("hScaledXSEC_%s", recoVar.Data()));
  hQCDAll[0]->SetTitle(TString::Format("QCD All Slices %s", year1.Data()));
  hQCDAll[1]->SetTitle(TString::Format("QCD All Slices %s", year2.Data()));
  hQCDAll[0]->SetName(TString::Format("QCD All Slices %s", year1.Data()));
  hQCDAll[1]->SetName(TString::Format("QCD All Slices %s", year2.Data()));
  */
  TString reason;
  //we do 17/18 so numerator is 17 --> [0] and denominator is 18 --> [1]

  for(int i =0; i<histoNamesTT.size(); i++)
  {
  	reason = "TT MC Overlay " + histoNamesTT[i];
  	ratioPlot(hSig_slice[0][i],hSig_slice[1][i], recoVar, reason);
  }

  /*
  for(int i =0; i<histoNamesQCD.size(); i++)
  {
  	reason = "QCD MC Overlay "+ histoNamesQCD[i];
  	ratioPlot(hQCD_slice[0][i],hQCD_slice[1][i], recoVar, reason);
  }*/


  reason = "TT MC Overlay all slices";
  ratioPlot(hSigAll[0], hSigAll[1],recoVar,reason);

  //reason = "QCD MC Overlay all slices";
  //ratioPlot(hQCDAll[0], hQCDAll[1],recoVar,reason);

  histoNamesQCD.clear();
  histoNamesTT.clear();
}



void ratioPlot(TH1F *hNum ,TH1F *hDenom, TString recoVar, TString reason)
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
  TString temp;
  if(globalIsNominalMC) temp = "NominalMC";
  else temp = "TT_Mtt";
  c1->Print(TString::Format("./ComparisonScaledIntegral/%s/%s/comparison_mc%s_%s_%sParton.pdf",temp.Data(),recoVar.Data(),
  							 globalYear1.Data(), globalYear2.Data(), reason.Data()),"pdf");
}
