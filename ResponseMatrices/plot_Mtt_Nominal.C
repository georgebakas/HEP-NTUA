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
void plotYearVar(TString year ,TString recoVar, TString varParton, TString varParticle );
void ratioPlot(TEfficiency *effNum, TEfficiency *effDenom,TH1F *hNum ,TH1F *hDenom, TString recoVar, TString reason);

void plot_Mtt_Nominal(TString year="2017")
{
	const int NVAR = 7;
	TString recoVar[NVAR] = {"mJJ", "ptJJ", "yJJ","jetPt0", "jetPt1","jetY0", "jetY1"};
	TString varParton[NVAR] = {"mTTbarParton", "ptTTbarParton", "yTTbarParton","partonPt0", "partonPt1", "partonY0", "partonY1"};
    TString varParticle[NVAR] = {"mJJGen", "ptJJGen", "yJJGen","genjetPt0", "genjetPt1", "genjetY0", "genjetY1"};
	
	for(int ivar = 0; ivar<NVAR-2; ivar++)
	{
		plotYearVar(year,recoVar[ivar],varParton[ivar], varParticle[ivar]);
		//break;
	}
}
void plotYearVar(TString year ,TString recoVar, TString varParton, TString varParticle)
{
  gStyle->SetOptStat(0);
  std::vector<TString> histoNames;
  histoNames.clear();

  histoNames.push_back("Signal_histo_Mtt"); 
  histoNames.push_back("Signal_histo_Nominal");
  
  //0 is 2017 and //1 is 2018
  //TT files

  TFile *infTT[2];
  infTT[0] = TFile::Open(TString::Format("%s/EqualBins/ResponsesEfficiency_%s.root", year.Data(), year.Data()));  
  infTT[1] = TFile::Open(TString::Format("%s/EqualBins/ResponsesEfficiencyNominalMC_%s.root", year.Data(), year.Data())); 

  TEfficiency *effParton[2], *effParticle[2];
  TEfficiency *accParton[2], *accParticle[2];

  TEfficiency *effParton_[2], *effParticle_[2];
  TEfficiency *accParton_[2], *accParticle_[2];

  //get the efficiencies and acceptances
  
  //get total and passed histo
  TH1F *hTotalPartonEff[2], *hPassedPartonEff[2];
  TH1F *hTotalPartonAcc[2], *hPassedPartonAcc[2];
  TH1F *hTotalParticleEff[2], *hPassedParticleEff[2];
  TH1F *hTotalParticleAcc[2], *hPassedParticleAcc[2];

  for(int i =0; i<2; i++)
  {
  	//get parton efficiency and make it into a th1
  	effParton[i] = (TEfficiency*)infTT[i]->Get(TString::Format("EfficiencyParton_%s",varParton.Data()));
  	effParton_[i] = (TEfficiency*)infTT[i]->Get(TString::Format("EfficiencyParton_%s",varParton.Data()));
  	hTotalPartonEff[i] = (TH1F*)effParton[i]->GetTotalHistogram();
    hPassedPartonEff[i] = (TH1F*)effParton[i]->GetPassedHistogram();
    hPassedPartonEff[i]->Divide(hTotalPartonEff[i]);

	//get particle efficiency and make it into a th1
	effParticle[i] = (TEfficiency*)infTT[i]->Get(TString::Format("EfficiencyParticle_%s",varParticle.Data()));
	effParticle_[i] = (TEfficiency*)infTT[i]->Get(TString::Format("EfficiencyParticle_%s",varParticle.Data()));
    hTotalParticleEff[i] = (TH1F*)effParticle[i]->GetTotalHistogram();
    hPassedParticleEff[i] = (TH1F*)effParticle[i]->GetPassedHistogram();
    hPassedParticleEff[i]->Divide(hTotalParticleEff[i]);

    //get parton acceptance and make it into a th1
    accParton[i] = (TEfficiency*)infTT[i]->Get(TString::Format("AcceptanceParton_%s",recoVar.Data()));
    accParton_[i] = (TEfficiency*)infTT[i]->Get(TString::Format("AcceptanceParton_%s",recoVar.Data()));
    hTotalPartonAcc[i] = (TH1F*)accParton[i]->GetTotalHistogram();
    hPassedPartonAcc[i] = (TH1F*)accParton[i]->GetPassedHistogram();
    hPassedPartonAcc[i]->Divide(hTotalPartonAcc[i]);

    //get particle acceptance and make it into a th1
    accParticle[i] = (TEfficiency*)infTT[i]->Get(TString::Format("AcceptanceParticle_%s",recoVar.Data()));
    accParticle_[i] = (TEfficiency*)infTT[i]->Get(TString::Format("AcceptanceParticle_%s",recoVar.Data()));
    hTotalParticleAcc[i] = (TH1F*)accParticle[i]->GetTotalHistogram();
    hPassedParticleAcc[i] = (TH1F*)accParticle[i]->GetPassedHistogram();
    hPassedParticleAcc[i]->Divide(hTotalParticleAcc[i]);

  }

  hPassedPartonEff[0]->SetTitle("Mtt 700-Inf");
  hPassedPartonEff[1]->SetTitle("TT Nominal");
  hPassedPartonEff[0]->SetName(TString::Format("Mtt 700-Inf Parton Eff %s", year.Data()));
  hPassedPartonEff[1]->SetName(TString::Format("TT Nominal MC Parton Eff %s", year.Data()));

  
  hPassedParticleEff[0]->SetTitle("Mtt 700-Inf");
  hPassedParticleEff[1]->SetTitle("TT Nominal");
  hPassedParticleEff[0]->SetName(TString::Format("Mtt 700-Inf Particle Eff %s", year.Data()));
  hPassedParticleEff[1]->SetName(TString::Format("TT Nominal MC Particle Eff %s", year.Data()));

  hPassedPartonAcc[0]->SetTitle("Mtt 700-Inf");
  hPassedPartonAcc[1]->SetTitle("TT Nominal");
  hPassedPartonAcc[0]->SetName(TString::Format("Mtt 700-Inf Parton Acc %s", year.Data()));
  hPassedPartonAcc[1]->SetName(TString::Format("TT Nominal MC Parton Acc %s", year.Data()));

  hPassedParticleAcc[0]->SetTitle("Mtt 700-Inf");
  hPassedParticleAcc[1]->SetTitle("TT Nominal");
  hPassedParticleAcc[0]->SetName(TString::Format("Mtt 700-Inf Particle Acc %s", year.Data()));
  hPassedParticleAcc[1]->SetName(TString::Format("TT Nominal MC Particle Acc %s", year.Data()));
  

  
  TString reason;
  
  reason = "Parton Efficiency " + year;
  ratioPlot(effParton_[0], effParton_[1],hPassedPartonEff[0], hPassedPartonEff[1],recoVar,reason);

  reason = "Particle Efficiency " + year;
  ratioPlot(effParticle_[0], effParticle_[1],hPassedParticleEff[0], hPassedParticleEff[1],recoVar,reason);

  reason = "Parton Acceptance " + year;
  ratioPlot(accParton_[0],accParton_[1],hPassedPartonAcc[0], hPassedPartonAcc[1],recoVar,reason);

  reason = "Particle Acceptance " + year;
  ratioPlot(accParticle_[0],accParticle_[1],hPassedParticleAcc[0], hPassedParticleAcc[1],recoVar,reason);
  
  histoNames.clear();

  
}



void ratioPlot(TEfficiency *effNum, TEfficiency *effDenom, TH1F *hNum ,TH1F *hDenom, TString recoVar, TString reason)
{

  TString titleNum = hNum->GetTitle();
  TString titleDenom = hDenom->GetTitle();
  TLegend *closureLegend = new TLegend(0.65,0.73,0.9,0.9);
  closureLegend->AddEntry(effNum,titleNum, "lep");
  closureLegend->AddEntry(effDenom,titleDenom, "lep");

  hNum->SetTitle(reason);
  hDenom->SetTitle(reason);
  //cout<<"---------------"<<endl;
  //cout<<reason<<endl;
  //cout<<"hNum name: "<<hNum->GetTitle()<<endl;
  //cout<<"hDenom name: "<<hDenom->GetTitle()<<endl;

  auto c1 = new TCanvas(reason+recoVar, reason+recoVar, 800,700);
  auto *closure_pad2 = new TPad("closure_pad2","closure_pad2",0.,0.,1.,0.45); 
  closure_pad2->Draw();
  closure_pad2->SetTopMargin(0.05);
  closure_pad2->SetBottomMargin(0.3);
  closure_pad2->SetGrid();

  auto *closure_pad1 = new TPad("closure_pad1","closure_pad1",0.,0.45,1.,1.);  
  closure_pad1->Draw();
  closure_pad1->SetBottomMargin(0.005);
  closure_pad1->cd();
  //closure_pad1->SetGrid();
  TH1F *hEfficiencyNum, *hEfficiencyDenom;
/*
  if(reason.Contains("Efficiency"))
  { 
  	if(reason.Contains("Parton")) hNum->GetYaxis()->SetRangeUser(0,0.1);
  	else hNum->GetYaxis()->SetRangeUser(0,0.4);
  }
  else hNum->GetYaxis()->SetRangeUser(hNum->GetMaximum()-0.4, hNum->GetMaximum()+0.2);
  */
  //effNum->GetYaxis()->SetTitleSize(20);
  //effNum->GetYaxis()->SetTitleFont(43);
  //effNum->GetYaxis()->SetTitleOffset(1.4); 
  //effNum->ResetAttLine();
  //effNum->ResetAttMarker();
  effNum->SetLineColor(kRed-7);
  effNum->SetMarkerStyle(22);
  effNum->SetMarkerColor(kRed-7);

  //effDenom->ResetAttLine();
  //effDenom->ResetAttMarker();
  effDenom->SetLineColor(kBlue-2); 
  effDenom->SetMarkerStyle(21);
  effDenom->SetMarkerColor(kBlue-2);
  
  effNum->Draw();
  effDenom->Draw("same"); 
  //hBkg_SR[0]->GetXaxis()->SetTitleOffset(1.5);
  //if(!recoVar.EqualTo("yJJ") && !recoVar.EqualTo("jetY0") && !recoVar.EqualTo("jetY1")) hNum->GetXaxis()->SetTitle(recoVar+" (GeV)");
  //else hNum->GetXaxis()->SetTitle(recoVar);
  // h2 settings
  
  //hBkg_CR[0]->ResetAttFill();
  
  
  closureLegend->Draw();
  //closure_pad1->SetLogy();
  
  
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

  reason.ReplaceAll(" ","");
  cout<<reason<<endl;
  hRatio->GetXaxis()->SetLabelSize(0.06);
  hRatio->GetYaxis()->SetRangeUser(0.3,1.6);
  hRatio->Divide(hDenom);
  hRatio->SetLineColor(kBlack);
  if(!recoVar.EqualTo("yJJ") && !recoVar.EqualTo("jetY0") && !recoVar.EqualTo("jetY1")) hRatio->GetXaxis()->SetTitle(recoVar+" (GeV)");
  else hRatio->GetXaxis()->SetTitle(recoVar);
  hRatio->Draw();
  c1->Print(TString::Format("./nominalHighMtt_comparison/%s/comparison_%s.pdf",recoVar.Data(), reason.Data()),"pdf");
}
