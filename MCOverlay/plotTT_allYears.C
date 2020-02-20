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
void plotVar(TString recoVar = "jetPt0");
void ratioPlot(TH1F *hNum ,TH1F *hDenom, TString recoVar, TString reason, TString level);
bool globalIsNominalMC;
TString globalYear1;

void plotTT_allYears(bool isNominalMC = false)
{
	globalIsNominalMC = isNominalMC;
	const int NVAR = 11;
	TString recoVar[NVAR]   = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1","jetPhi0","jetPhi1","mTop0", "mTop1"};

	for(int ivar = 0; ivar<NVAR; ivar++)
	{
		plotVar(recoVar[ivar]);
		//break;
	}
}
void plotVar(TString recoVar = "jetPt0")
{
  initFilesMapping();
  //0 is 2017 and //1 is 2018
  //TT files
  TString nominal = "";
  if(globalIsNominalMC) nominal = "Nominal";

  TFile *infTT[3];
  infTT[0] = TFile::Open(filesttbar[TString::Format("2016%s",nominal.Data())]); 
  infTT[1] = TFile::Open(filesttbar[TString::Format("2017%s",nominal.Data())]);  
  infTT[2] = TFile::Open(filesttbar[TString::Format("2018%s",nominal.Data())]);
  TH1F *hReco[3], *hParton[3], *hParticle[3];
  

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
  std::vector<Color_t> col = {kBlue, kRed, kGreen-6}; //this refers to each year
  TString year[] = {"'16", "'17", "'18"};
  TString level[] = {"Reco", "Parton", "Particle"};
  TCanvas *can = new TCanvas(TString::Format("canParton %s", recoVar.Data()), 
                   TString::Format("canParton %s", recoVar.Data()), 800, 700);
  TLegend *leg = new TLegend(0.65,0.73,0.9,0.9);
  //scaled TT to xsec and LUMI
  if(!globalIsNominalMC) nominal = "Mtt";
  for(int i = 0; i < 3; i++) //loop on all years...
  {
	  hParton[i] = (TH1F*)infTT[i]->Get(TString::Format("hPartonScaledXSEC_%s", recoVar.Data()));
	  hParton[i]->SetTitle(TString::Format("TT Parton %s MC", nominal.Data()));
	  hParton[i]->SetName(TString::Format("TT Parton %s MC", nominal.Data()));
	 
	  hParton[i]->SetLineColor(col[i]);
	  hParton[i]->Scale(1./hParton[i]->Integral());
  }	

  for(int y=0; y<3; y++)
  {
  	leg->AddEntry(hParton[y],year[y],"lep");
  }

  can->cd();
  hParton[0]->Draw();
  gPad->SetLogy();
  hParton[1]->Draw("same");
  hParton[2]->Draw("same");
  leg->Draw();

	can->Print(TString::Format("./allYears/%s/%s/allyears_%sParton.pdf",nominal.Data(),
									recoVar.Data(), recoVar.Data()),"pdf");
  
  
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