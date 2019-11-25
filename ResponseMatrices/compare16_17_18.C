#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLatex.h"

using std::cin;
using std::cout;
using std::endl;

void compare16_17_18(TString recoVar = "jetPt0",TString partonVar = "partonPt0", TString particleVar = "genjetPt0")
{
   gStyle->SetOptStat(0);
   gStyle->SetPaintTextFormat("4.4f");
   std::vector<Color_t> colors = {kBlue,kRed, kGreen,kBlack, kMagenta};
	
   bool isAngular = false;
   if(recoVar.EqualTo("yJJ") || recoVar.EqualTo("jetY0") || recoVar.EqualTo("jetY1")) isAngular = true;	   
 	
   TFile *eff[3];
   //for efficiencies and acceptance
   eff[0] = TFile::Open("./2016/UnequalBins/ResponsesEfficiency_2016.root"); 
   eff[1] = TFile::Open("./2017/UnequalBins/ResponsesEfficiency_2017.root");
   eff[2] = TFile::Open("./2018/UnequalBins/ResponsesEfficiency_2018.root");

   TFile *oldInf[2];
   oldInf[0] = TFile::Open("PartonEfficiencyAll_July19.root");
   oldInf[1] = TFile::Open("ParticleEfficiencyAll_July19.root");
   
   //things we need:
   TEfficiency *eff16[2], *acc16[2];   
   TEfficiency *eff17[2], *acc17[2];
   TEfficiency *eff18[2], *acc18[2];

   TEfficiency *effOld16[2], *accOld16[2];

   effOld16[0] = (TEfficiency*)oldInf[0]->Get(TString::Format("Eff_%s_Parton_Nominal",recoVar.Data()));
   accOld16[0] = (TEfficiency*)oldInf[0]->Get(TString::Format("Eff_%s_Parton_Nominal_common",recoVar.Data()));

   effOld16[1] = (TEfficiency*)oldInf[1]->Get(TString::Format("Eff_%s_Gen_Nominal",recoVar.Data()));
   accOld16[1] = (TEfficiency*)oldInf[1]->Get(TString::Format("Eff_%s_Gen_Nominal_common",recoVar.Data()));

   
   //1. Efficiency for 2016 and acceptance for same year for the parton level
   eff16[0] = (TEfficiency*)eff[0]->Get(TString::Format("EfficiencyParton_%s", partonVar.Data()));
   acc16[0] = (TEfficiency*)eff[0]->Get(TString::Format("AcceptanceParton_%s", recoVar.Data()));
   
   //2. Efficiency for 2017 and acceptance for same year for parton level
   eff17[0] = (TEfficiency*)eff[1]->Get(TString::Format("EfficiencyParton_%s", partonVar.Data()));
   acc17[0] = (TEfficiency*)eff[1]->Get(TString::Format("AcceptanceParton_%s", recoVar.Data()));
   
   //3. Efficiency for 2018 and acceptance for same year for parton level
   eff18[0] = (TEfficiency*)eff[2]->Get(TString::Format("EfficiencyParton_%s", partonVar.Data()));
   acc18[0] = (TEfficiency*)eff[2]->Get(TString::Format("AcceptanceParton_%s", recoVar.Data()));
  


   //1. Efficiency for 2016 and acceptance for same year for the particle level
   eff16[1] = (TEfficiency*)eff[0]->Get(TString::Format("EfficiencyParticle_%s", particleVar.Data()));
   acc16[1] = (TEfficiency*)eff[0]->Get(TString::Format("AcceptanceParticle_%s", recoVar.Data()));
   
   //2. Efficiency for 2017 and acceptance for same year for particle level
   eff17[1] = (TEfficiency*)eff[1]->Get(TString::Format("EfficiencyParticle_%s", particleVar.Data()));
   acc17[1] = (TEfficiency*)eff[1]->Get(TString::Format("AcceptanceParticle_%s", recoVar.Data()));
   
   //3. Efficiency for 2018 and acceptance for same year for particle level
   eff18[1] = (TEfficiency*)eff[2]->Get(TString::Format("EfficiencyParticle_%s", particleVar.Data()));
   acc18[1] = (TEfficiency*)eff[2]->Get(TString::Format("AcceptanceParticle_%s", recoVar.Data()));


   for(int i =0; i<sizeof(eff16)/sizeof(eff16[0]); i++)
   {
	   eff16[i]->SetMarkerStyle(21);
	   eff17[i]->SetMarkerStyle(22);
	   eff18[i]->SetMarkerStyle(23);
	   effOld16[i]->SetMarkerStyle(20);
	   eff16[i]->SetMarkerColor(colors[0]);
	   eff17[i]->SetMarkerColor(colors[1]);
	   eff18[i]->SetMarkerColor(colors[2]);
	   effOld16[i]->SetMarkerColor(colors[3]);
	   eff16[i]->SetLineColor(colors[0]);
	   eff17[i]->SetLineColor(colors[1]);
	   eff18[i]->SetLineColor(colors[2]);
	   effOld16[i]->SetLineColor(colors[3]);
	   
	   acc16[i]->SetMarkerStyle(21);
	   acc17[i]->SetMarkerStyle(22);
	   acc18[i]->SetMarkerStyle(23);
	   accOld16[i]->SetMarkerStyle(20);
	   acc16[i]->SetLineColor(colors[0]);
	   acc17[i]->SetLineColor(colors[1]);
	   acc18[i]->SetLineColor(colors[2]);
	   acc16[i]->SetMarkerColor(colors[0]);
	   acc17[i]->SetMarkerColor(colors[1]);
	   acc18[i]->SetMarkerColor(colors[2]);
	   accOld16[i]->SetMarkerColor(colors[3]);
	   accOld16[i]->SetLineColor(colors[3]);
   }

   TLegend *effLeg = new TLegend(0.65,0.7,0.9,0.9);
   effLeg->AddEntry(eff16[0], "tTagger '16", "lp");
   effLeg->AddEntry(eff17[0], "tTagger '17", "lp");
   effLeg->AddEntry(eff18[0], "tTagger '18", "lp");
   effLeg->AddEntry(effOld16[0], "2016 analysis", "lp");

   TCanvas *can_eff[2], *can_acc[2];
   TString phaseSpace[] = {"Parton", "Particle"};
   TString years[] = {"2016", "2017", "2018"};

   for(int i =0; i<sizeof(eff16)/sizeof(eff16[0]); i++)
   {  
	   can_eff[i] = new TCanvas(TString::Format("Efficiency can_%s",phaseSpace[i].Data()), TString::Format("Efficiency can_%s",phaseSpace[i].Data()), 700, 600);
	   eff18[i]->SetTitle(TString::Format("%s Efficiency '16,'17,'18;%s (GeV);Efficiency",phaseSpace[i].Data(),recoVar.Data())); 
	   eff18[i]->Draw();
	   eff17[i]->Draw("same");
	   eff16[i]->Draw("same");
	   effOld16[i]->Draw("same");
	   effLeg->Draw();
      double xmin, xmax, ymin, ymax;
      can_eff[i]->cd();
      gPad->GetRangeAxis(xmin,ymin,xmax,ymax);
      if(i==0)gPad->Range(xmin,0,xmax,0.1);
      else gPad->Range(xmin,0,xmax,0.4);
	   can_eff[i]->Print(TString::Format("plots/%s/Efficiency%s_%s.pdf",recoVar.Data(),phaseSpace[i].Data(),recoVar.Data()),"pdf");
	   
	   can_acc[i] = new TCanvas(TString::Format("Acceptance can_%s",phaseSpace[i].Data()), TString::Format("Acceptance can_%s",phaseSpace[i].Data()), 700, 600);
	   acc18[i]->SetTitle(TString::Format("%s Acceptance '16,'17,'18;%s (GeV);Acceptance",phaseSpace[i].Data(),recoVar.Data()));  
	   acc18[i]->Draw();
	   acc17[i]->Draw("same");
	   acc16[i]->Draw("same");
	   accOld16[i]->Draw("same");
      effLeg->Draw(); 
      can_acc[i]->cd();
      gPad->GetRangeAxis(xmin,ymin,xmax,ymax);
      gPad->Range(xmin,0.4,xmax,1.1);
	   can_eff[i]->Print(TString::Format("plots/%s/Acceptance%s_%s.pdf",recoVar.Data() ,phaseSpace[i].Data(),recoVar.Data()),"pdf");
   }

   
   TH2F *hResponses[3][2]; //3 is for years , 16-0, 17-1, 18-2 and 2 is for parton particle: parton-0, particle 1
   TCanvas *canResponse[3][2];
   for(int iy = 0; iy<sizeof(years)/sizeof(years[0]); iy++)
   {
   	for(int i =0; i<sizeof(eff16)/sizeof(eff16[0]); i++)
   	{
   		canResponse[iy][i] = new TCanvas(TString::Format("Response Reco-%s %s %s",phaseSpace[i].Data(), recoVar.Data(), years[iy].Data()),
   		 								 TString::Format("Response Reco-%s %s %s",phaseSpace[i].Data(), recoVar.Data(), years[iy].Data()),800,600);
   		hResponses[iy][i] = (TH2F*)eff[iy]->Get(TString::Format("h%sResponse_%s",phaseSpace[i].Data() ,recoVar.Data()));
   		hResponses[iy][i]->Scale(1./hResponses[iy][i]->Integral());
   		hResponses[iy][i]->SetTitle(TString::Format("Response Reco-%s %s %s",phaseSpace[i].Data(), recoVar.Data(), years[iy].Data()));
   		if(isAngular) 
   		{	
   			hResponses[iy][i]->GetYaxis()->SetTitle(recoVar.Data());
   			if(i==0) hResponses[iy][i]->GetXaxis()->SetTitle(partonVar.Data());
   			else hResponses[iy][i]->GetXaxis()->SetTitle(particleVar.Data());
   		}
   		else
   		{

   			hResponses[iy][i]->GetYaxis()->SetTitle(TString::Format("%s (GeV)",recoVar.Data()));
   			if(i==0) hResponses[iy][i]->GetXaxis()->SetTitle(TString::Format("%s (GeV)",partonVar.Data()));
   			else hResponses[iy][i]->GetXaxis()->SetTitle(TString::Format("%s (GeV)",particleVar.Data()));
  		}
   		hResponses[iy][i]->Draw("colz text");
   		canResponse[iy][i]->Print(TString::Format("%s/UnequalBins/%s/%sResponseMatrix_%s.pdf",years[iy].Data(),recoVar.Data() ,phaseSpace[i].Data(),recoVar.Data()),"pdf");
   	}
   }


  
}
