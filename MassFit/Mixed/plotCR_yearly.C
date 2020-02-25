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

/*
	Purpose of this code is to plot the CR between data and MC for each year
*/
using std::cin;
using std::cout;
using std::endl;

#include "TemplateConstants.h"
void plotYearVar(TString year, TString recoVar = "jetPt0");
void ratioPlot(TH1F *hNum ,TH1F *hDenom_0, TH1F *hDenom_1, TString recoVar, TString reason);
bool globalIsNominalMC;
TString globalYear1;

void plotCR_yearly(TString year = "2016")
{
	globalYear1 = year;
	const int NVAR = 9;
	TString recoVar[NVAR]   = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1", "mTop", "jetMassSoftDrop"};
               
	for(int ivar = 0; ivar<NVAR; ivar++)
	{
		plotYearVar(year,recoVar[ivar]);
		//break;
	}
}
void plotYearVar(TString year, TString recoVar = "jetPt0")
{
  //TT files 
  TFile *infTT[3];
  infTT[0] = TFile::Open(TString::Format("%s/Histo_Data_%s_100_reduced_Loose.root",year.Data(),year.Data())); //data
  infTT[1] = TFile::Open(TString::Format("%s/Histo_TT_NominalMC_100_reduced_Loose.root",year.Data())); //nominal mc
  infTT[2] = TFile::Open(TString::Format("%s/Histo_TT_Mtt-700toInf_100_reduced_Loose.root",year.Data()));  //mtt mc

  TH1F *hCR_Data, *hCR_Nominal, *hCR_Mtt;
 
  //scaled TT to xsec and LUMI
  hCR_Data = (TH1F*)infTT[0]->Get(TString::Format("hWt_%s_0btag", recoVar.Data()));
  hCR_Nominal = (TH1F*)infTT[1]->Get(TString::Format("hWt_%s_0btag_expYield", recoVar.Data()));
  hCR_Mtt = (TH1F*)infTT[2]->Get(TString::Format("hWt_%s_0btag_expYield", recoVar.Data())); //hWt_jetPt0_0btag_expYield
  
  hCR_Mtt->SetTitle(TString::Format("TT %s Mtt",year.Data()));
  hCR_Nominal->SetTitle(TString::Format("TT %s Nominal", year.Data()));
  hCR_Data->SetTitle(TString::Format("TT %s Data", year.Data()));
  
  hCR_Mtt->SetName(TString::Format("TT %s Mtt", year.Data()));
  hCR_Nominal->SetName(TString::Format("TT %s Nominal", year.Data()));
  hCR_Data->SetName(TString::Format("TT %s Data", year.Data()));
  
  TString reason;  
  reason = "CR_ShapeComparison";
  ratioPlot(hCR_Data, hCR_Nominal, hCR_Mtt, recoVar, reason);
  
}



void ratioPlot(TH1F *hNum ,TH1F *hDenom_0, TH1F *hDenom_1, TString recoVar, TString reason)
{
  TString titleNum = hNum->GetTitle();
  TString titleDenom_0 = hDenom_0->GetTitle();
  TString titleDenom_1 = hDenom_1->GetTitle();
  TLegend *closureLegend = new TLegend(0.65,0.73,0.9,0.9);
  closureLegend->AddEntry(hNum,titleNum, "lep");
  closureLegend->AddEntry(hDenom_0,titleDenom_0, "lep");
  closureLegend->AddEntry(hDenom_1,titleDenom_1, "lep");

  hNum->Scale(1./hNum->Integral());
  hDenom_0->Scale(1./hDenom_0->Integral());
  hDenom_1->Scale(1./hDenom_1->Integral());
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
  hNum->SetLineColor(kBlue-2);
  hNum->SetMarkerStyle(22);
  hNum->SetMarkerColor(kBlue-2); 
  //hBkg_SR[0]->GetXaxis()->SetTitleOffset(1.5);
  if(!recoVar.EqualTo("yJJ") || !recoVar.EqualTo("jetY0") || !recoVar.EqualTo("jetY1")) hNum->GetXaxis()->SetTitle(recoVar+" (GeV)");
  else hNum->GetXaxis()->SetTitle(recoVar);
  // h2 settings
  
  //hBkg_CR[0]->ResetAttFill();
  hDenom_0->ResetAttLine();
  hDenom_0->ResetAttMarker();
  hDenom_0->SetLineColor(kRed-7); 
  hDenom_0->SetMarkerStyle(21);
  hDenom_0->SetMarkerColor(kRed-7);

  hDenom_1->ResetAttLine();
  hDenom_1->ResetAttMarker();
  hDenom_1->SetLineColor(kGreen+3); 
  hDenom_1->SetMarkerStyle(21);
  hDenom_1->SetMarkerColor(kGreen+3);
  
  hDenom_0->Draw();
  hDenom_1->Draw("same");
  hNum->Draw("same");
  
  closureLegend->Draw();
  closure_pad1->SetLogy();
  
  TH1F *hRatio[2];
  closure_pad2->cd();
  for(int ir = 0; ir<2; ir++)
  {
  	hRatio[ir] = (TH1F*)hNum->Clone(TString::Format("hNumerator_%d",ir)); 
  	hRatio[ir]->SetTitle("");
  	hRatio[ir]->ResetAttMarker();
  	hRatio[ir]->GetYaxis()->SetTitle("#frac{Data}{MC}");
  	hRatio[ir]->GetYaxis()->SetTitleSize(20);
  	hRatio[ir]->GetYaxis()->SetTitleFont(43);
  	hRatio[ir]->GetYaxis()->SetTitleOffset(1.55);
 	hRatio[ir]->GetYaxis()->SetLabelFont(43);
  	hRatio[ir]->GetYaxis()->SetLabelSize(15);
  	hRatio[ir]->GetXaxis()->SetTitleSize(0.06);
  	hRatio[ir]->GetXaxis()->SetLabelSize(0.06);
  	if(!recoVar.EqualTo("yJJ") || !recoVar.EqualTo("jetY0") || !recoVar.EqualTo("jetY1")) hRatio[ir]->GetXaxis()->SetTitle(recoVar+" (GeV)");
    else hRatio[ir]->GetXaxis()->SetTitle(recoVar);
  }
  
  hRatio[0]->Divide(hDenom_0);
  hRatio[1]->Divide(hDenom_1);
  hRatio[0]->SetLineColor(kRed-7);
  hRatio[1]->SetLineColor(kGreen+3);


  hRatio[0]->Draw();
  hRatio[1]->Draw("same");
  c1->Print(TString::Format("./YearlyCRShapeComparison/%s/CRCShapeComparison_%s.pdf",globalYear1.Data(),recoVar.Data()),"pdf");
}