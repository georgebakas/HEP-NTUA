#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLatex.h"

#include "TemplateConstants_Response.h"
using std::cin;
using std::cout;
using std::endl;

void plotEfficiencyResponse(TString recoVar = "jetPt0",TString partonVar = "partonPt0", TString particleVar = "genjetPt0",
					 bool isEqual = true, bool isNominal = true);

void compare16_17_18(bool isEqual, bool isNominal)
{
  const int NVAR =3;
  TString varReco[NVAR]   = {"chi", "cosTheta_0", "cosTheta_1"};
  TString varParton[NVAR] = {"chiParton", "cosThetaParton_0", "cosThetaParton_1"};
  TString varParticle[NVAR] = {"chiParticle", "cosThetaParticle_0", "cosThetaParticle_1"};

  for(int ivar = 0; ivar<NVAR; ivar++)
  {
  	plotEfficiencyResponse(varReco[ivar], varParton[ivar], varParticle[ivar], isEqual, isNominal);
  }
}

void plotEfficiencyResponse(TString recoVar = "jetPt0",TString partonVar = "partonPt0", TString particleVar = "genjetPt0",
					 bool isEqual = true, bool isNominal = true)
{
   initFilesMapping();
   gStyle->SetOptStat(0);
   gStyle->SetPaintTextFormat("2.2f");
   std::vector<Color_t> colors = {kBlue,kRed, kGreen,kBlack, kMagenta};

   bool isAngular = true;

   TFile *eff[3];
   //for efficiencies and acceptance
   TString binning = "EqualBins";
   if(!isEqual) binning = "UnequalBins";

   TString nominal = "";
   if(isNominal) nominal = "NominalMC";
   //TEfficiency *eff16[2];

   eff[0] = TFile::Open(TString::Format("./2016/%s/ResponsesEfficiency%s_2016.root",binning.Data(), nominal.Data()));
   eff[1] = TFile::Open(TString::Format("./2017/%s/ResponsesEfficiency%s_2017.root",binning.Data(), nominal.Data()));
   eff[2] = TFile::Open(TString::Format("./2018/%s/ResponsesEfficiency%s_2018.root",binning.Data(), nominal.Data()));


   //things we need:
   TEfficiency *eff16[2], *acc16[2];
   TEfficiency *eff17[2], *acc17[2];
   TEfficiency *eff18[2], *acc18[2];



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
	   eff16[i]->SetMarkerColor(colors[0]);
	   eff17[i]->SetMarkerColor(colors[1]);
	   eff18[i]->SetMarkerColor(colors[2]);
	   eff16[i]->SetLineColor(colors[0]);
	   eff17[i]->SetLineColor(colors[1]);
	   eff18[i]->SetLineColor(colors[2]);

	   acc16[i]->SetMarkerStyle(21);
	   acc17[i]->SetMarkerStyle(22);
	   acc18[i]->SetMarkerStyle(23);
	   acc16[i]->SetLineColor(colors[0]);
	   acc17[i]->SetLineColor(colors[1]);
	   acc18[i]->SetLineColor(colors[2]);
	   acc16[i]->SetMarkerColor(colors[0]);
	   acc17[i]->SetMarkerColor(colors[1]);
	   acc18[i]->SetMarkerColor(colors[2]);
   }

   TLegend *effLeg = new TLegend(0.65,0.73,0.9,0.9);
   effLeg->AddEntry(eff16[0], "tTagger '16", "lp");
   effLeg->AddEntry(eff17[0], "tTagger '17", "lp");
   effLeg->AddEntry(eff18[0], "tTagger '18", "lp");

   TCanvas *can_eff[2], *can_acc[2];
   TString phaseSpace[] = {"Parton", "Particle"};
   TString years[] = {"2016", "2017", "2018"};

   for(int i =0; i<sizeof(eff16)/sizeof(eff16[0]); i++)
   {
	   can_eff[i] = new TCanvas(TString::Format("Efficiency can_%s%s",recoVar.Data(),phaseSpace[i].Data()), TString::Format("Efficiency can_%s%s",recoVar.Data(),phaseSpace[i].Data()), 700, 600);
	   eff18[i]->SetTitle(TString::Format("%s Efficiency '16,'17,'18 %s;%s;Efficiency",phaseSpace[i].Data(), nominal.Data(),recoVar.Data()));
	   eff18[i]->Draw();
	   eff17[i]->Draw("same");
	   eff16[i]->Draw("same");
	   effLeg->Draw();

		gPad->Update();
		auto graph = eff18[i]->GetPaintedGraph();
		graph->GetXaxis()->SetRangeUser(BNDmin[recoVar],BNDmax[recoVar]);
		if(recoVar.Contains("cosTheta"))graph->GetXaxis()->SetTitle("|"+recoVar+"|");
		else graph->GetXaxis()->SetTitle(recoVar.Data());
		if(i==0)
		{
			graph->SetMinimum(0.);
			graph->SetMaximum(0.1);
		}
		else
		{
			graph->SetMinimum(0.0);
			graph->SetMaximum(0.4);
		}
		gPad->Update();

      //if(i==0)gPad->Range(xmin,0,xmax,0.2);
      //else gPad->Range(xmin,0,xmax,0.4);
	   can_eff[i]->Print(TString::Format("eff_acc/%s/%sEfficiency%s_%s.pdf",binning.Data(),phaseSpace[i].Data(),nominal.Data(),recoVar.Data()),"pdf");

	   can_acc[i] = new TCanvas(TString::Format("Acceptance can_%s%s",recoVar.Data(),phaseSpace[i].Data()), TString::Format("Acceptance can_%s%s",recoVar.Data(),phaseSpace[i].Data()), 700, 600);
	   acc18[i]->SetTitle(TString::Format("%s Acceptance '16,'17,'18 %s;%s;Acceptance",phaseSpace[i].Data(), nominal.Data(),recoVar.Data()));
	   acc18[i]->Draw();
	   acc17[i]->Draw("same");
	   acc16[i]->Draw("same");
     effLeg->Draw();
   	 gPad->Update();
		auto graphAcc = acc18[i]->GetPaintedGraph();
	 	graphAcc->SetMinimum(0.5);
		graphAcc->SetMaximum(1.2);
		graphAcc->GetXaxis()->SetRangeUser(BNDmin[recoVar],BNDmax[recoVar]);
		if(recoVar.Contains("cosTheta"))graphAcc->GetXaxis()->SetTitle("|"+recoVar+"|");
		else graphAcc->GetXaxis()->SetTitle(recoVar.Data());
		gPad->Update();

		 can_acc[i]->Print(TString::Format("eff_acc/%s/%sAcceptance%s_%s.pdf",binning.Data(),phaseSpace[i].Data(),nominal.Data(),recoVar.Data()),"pdf");
   }


   //TString years[] = {"2016", "2017", "2018"};
   //TString phaseSpace[] = {"Parton", "Particle"};
   TH2F *hResponses[3][2]; //3 is for years , 16-0, 17-1, 18-2 and 2 is for parton particle: parton-0, particle 1
   TCanvas *canResponse[3][2];
   for(int iy = 0; iy<sizeof(years)/sizeof(years[0]); iy++)
   {
   	for(int i =0; i<sizeof(eff16)/sizeof(eff16[0]); i++)
   	{
   		canResponse[iy][i] = new TCanvas(TString::Format("Response Reco-%s %s %s",phaseSpace[i].Data(), recoVar.Data(), years[iy].Data()),
   		 								 TString::Format("Response Reco-%s %s %s",phaseSpace[i].Data(), recoVar.Data(), years[iy].Data()),800,600);
   		hResponses[iy][i] = (TH2F*)eff[iy]->Get(TString::Format("h%sResponse_%s",phaseSpace[i].Data() ,recoVar.Data()));
   		//hResponses[iy][i]->Scale(1./hResponses[iy][i]->Integral());
      if(recoVar.EqualTo("jetPt0")) cout<<hResponses[iy][i]->Integral()<<endl;
   		hResponses[iy][i]->SetTitle(TString::Format("Response Reco-%s %s %s %s",phaseSpace[i].Data(), recoVar.Data(), years[iy].Data(), nominal.Data()));
   		if(recoVar.Contains("cosTheta"))
   		{
   			hResponses[iy][i]->GetYaxis()->SetTitle("|"+recoVar+"|");
   			if(i==0) hResponses[iy][i]->GetXaxis()->SetTitle("|"+partonVar+"|");
   			else hResponses[iy][i]->GetXaxis()->SetTitle("|"+particleVar+"|");
   		}
   		else
   		{

   			hResponses[iy][i]->GetYaxis()->SetTitle(recoVar.Data());
   			if(i==0) hResponses[iy][i]->GetXaxis()->SetTitle(partonVar.Data());
   			else hResponses[iy][i]->GetXaxis()->SetTitle(particleVar.Data());
  		}
   		hResponses[iy][i]->Draw("colz text");
   		canResponse[iy][i]->Print(TString::Format("%s/%s/%sResponseMatrix%s_%s.pdf",years[iy].Data(),binning.Data(),phaseSpace[i].Data(),nominal.Data(),recoVar.Data()),"pdf");
   	}
   }

}
