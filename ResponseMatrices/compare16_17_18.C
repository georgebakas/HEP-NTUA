#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLatex.h"

#include "TemplateConstants.h"
using std::cin;
using std::cout;
using std::endl;

void plotEfficiencyResponse(TString recoVar = "jetPt0",TString partonVar = "partonPt0", TString particleVar = "genjetPt0",
					 bool isEqual = true, bool isNominal = true);

void compare16_17_18(bool isEqual = false, bool isNominal = true)
{
  const int NVAR =7;
  TString varReco[NVAR]   = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1"};
  TString varParton[NVAR] = {"mTTbarParton", "ptTTbarParton", "yTTbarParton","partonPt0", "partonPt1", "partonY0", "partonY1"};
  TString varParticle[NVAR] = {"mJJGen", "ptJJGen", "yJJGen","genjetPt0", "genjetPt1", "genjetY0", "genjetY1"};

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
   std::vector<Color_t> colors = {kBlue, kMagenta, kRed, kGreen, kBlack};

   bool isAngular = false;
   if(recoVar.EqualTo("yJJ") || recoVar.EqualTo("jetY0") || recoVar.EqualTo("jetY1")) isAngular = true;

   TFile *eff[4];
   //for efficiencies and acceptance
   TString binning = "EqualBins";
   if(!isEqual) binning = "UnequalBins";

   TString nominal = "";
   if(isNominal) nominal = "NominalMC";

   eff[0] = TFile::Open(TString::Format("./UL2016_preVFP/%s/ResponsesEfficiency%s_2016_preVFP.root",binning.Data(), nominal.Data()));
   eff[1] = TFile::Open(TString::Format("./UL2016_postVFP/%s/ResponsesEfficiency%s_2016_postVFP.root",binning.Data(), nominal.Data()));
   eff[2] = TFile::Open(TString::Format("./UL2017/%s/ResponsesEfficiency%s_2017.root",binning.Data(), nominal.Data()));
   eff[3] = TFile::Open(TString::Format("./UL2018/%s/ResponsesEfficiency%s_2018.root",binning.Data(), nominal.Data()));

   TFile *oldInf[2];
   oldInf[0] = TFile::Open("PartonEfficiencyAll_July19.root");
   oldInf[1] = TFile::Open("ParticleEfficiencyAll_July19.root");

   //things we need:
   TEfficiency *eff16_pre[2], *acc16_pre[2];
   TEfficiency *eff16_post[2], *acc16_post[2];
   TEfficiency *eff17[2], *acc17[2];
   TEfficiency *eff18[2], *acc18[2];

   TEfficiency *effOld16[2], *accOld16[2];

   effOld16[0] = (TEfficiency*)oldInf[0]->Get(TString::Format("Eff_%s_Parton_Nominal",recoVar.Data()));
   accOld16[0] = (TEfficiency*)oldInf[0]->Get(TString::Format("Eff_%s_Parton_Nominal_common",recoVar.Data()));

   effOld16[1] = (TEfficiency*)oldInf[1]->Get(TString::Format("Eff_%s_Gen_Nominal",recoVar.Data()));
   accOld16[1] = (TEfficiency*)oldInf[1]->Get(TString::Format("Eff_%s_Gen_Nominal_common",recoVar.Data()));


   //1. Efficiency for 2016_preVFP and acceptance for same year for the parton level
   eff16_pre[0] = (TEfficiency*)eff[0]->Get(TString::Format("EfficiencyParton_%s", partonVar.Data()));
   acc16_pre[0] = (TEfficiency*)eff[0]->Get(TString::Format("AcceptanceParton_%s", recoVar.Data()));

   //2. Efficiency for 2016_postVFP and acceptance for same year for the parton level
   eff16_post[0] = (TEfficiency*)eff[1]->Get(TString::Format("EfficiencyParton_%s", partonVar.Data()));
   acc16_post[0] = (TEfficiency*)eff[1]->Get(TString::Format("AcceptanceParton_%s", recoVar.Data()));

   //3. Efficiency for 2017 and acceptance for same year for parton level
   eff17[0] = (TEfficiency*)eff[2]->Get(TString::Format("EfficiencyParton_%s", partonVar.Data()));
   acc17[0] = (TEfficiency*)eff[2]->Get(TString::Format("AcceptanceParton_%s", recoVar.Data()));

   //4. Efficiency for 2018 and acceptance for same year for parton level
   eff18[0] = (TEfficiency*)eff[3]->Get(TString::Format("EfficiencyParton_%s", partonVar.Data()));
   acc18[0] = (TEfficiency*)eff[3]->Get(TString::Format("AcceptanceParton_%s", recoVar.Data()));



   //1. Efficiency for 2016_preVFP and acceptance for same year for the particle level
   eff16_pre[1] = (TEfficiency*)eff[0]->Get(TString::Format("EfficiencyParticle_%s", particleVar.Data()));
   acc16_pre[1] = (TEfficiency*)eff[0]->Get(TString::Format("AcceptanceParticle_%s", recoVar.Data()));

   //1. Efficiency for 2016_postVFP and acceptance for same year for the particle level
   eff16_post[1] = (TEfficiency*)eff[1]->Get(TString::Format("EfficiencyParticle_%s", particleVar.Data()));
   acc16_post[1] = (TEfficiency*)eff[1]->Get(TString::Format("AcceptanceParticle_%s", recoVar.Data()));

   //2. Efficiency for 2017 and acceptance for same year for particle level
   eff17[1] = (TEfficiency*)eff[2]->Get(TString::Format("EfficiencyParticle_%s", particleVar.Data()));
   acc17[1] = (TEfficiency*)eff[2]->Get(TString::Format("AcceptanceParticle_%s", recoVar.Data()));

   //3. Efficiency for 2018 and acceptance for same year for particle level
   eff18[1] = (TEfficiency*)eff[3]->Get(TString::Format("EfficiencyParticle_%s", particleVar.Data()));
   acc18[1] = (TEfficiency*)eff[3]->Get(TString::Format("AcceptanceParticle_%s", recoVar.Data()));


   for(int i =0; i<sizeof(eff16)/sizeof(eff16[0]); i++)
   {
	   eff16_pre[i]->SetMarkerStyle(21);
	   eff16_post[i]->SetMarkerStyle(22);
	   eff17[i]->SetMarkerStyle(23);
	   eff18[i]->SetMarkerStyle(24);
	   effOld16[i]->SetMarkerStyle(20);
	   
	   eff16_pre[i]->SetMarkerColor(colors[0]);
	   eff16_post[i]->SetMarkerColor(colors[0]);
	   eff17[i]->SetMarkerColor(colors[2]);
	   eff18[i]->SetMarkerColor(colors[3]);
	   effOld16[i]->SetMarkerColor(colors[4]);
	   
	   eff16_pre[i]->SetLineColor(colors[0]);
	   eff16_post[i]->SetLineColor(colors[1]);
	   eff17[i]->SetLineColor(colors[2]);
	   eff18[i]->SetLineColor(colors[3]);
	   effOld16[i]->SetLineColor(colors[4]);

	   acc16_pre[i]->SetMarkerStyle(21);
	   acc16_post[i]->SetMarkerStyle(22);
	   acc17[i]->SetMarkerStyle(23);
	   acc18[i]->SetMarkerStyle(24);
	   accOld16[i]->SetMarkerStyle(20);
	   
	   acc16_pre[i]->SetLineColor(colors[0]);
	   acc16_post[i]->SetLineColor(colors[1]);
	   acc17[i]->SetLineColor(colors[2]);
	   acc18[i]->SetLineColor(colors[3]);
	   accOld16[i]->SetLineColor(colors[4]);

	   acc16_pre[i]->SetMarkerColor(colors[0]);
	   acc16_post[i]->SetMarkerColor(colors[1]);
	   acc17[i]->SetMarkerColor(colors[2]);
	   acc18[i]->SetMarkerColor(colors[3]);
	   accOld16[i]->SetMarkerColor(colors[4]);
	   
   }

   TLegend *effLeg = new TLegend(0.65,0.73,0.9,0.9);
   effLeg->AddEntry(eff16_pre[0], "tTagger '16 preVFP", "lp");
   effLeg->AddEntry(eff16_post[0], "tTagger '16 postVFP", "lp");
   effLeg->AddEntry(eff17[0], "tTagger '17", "lp");
   effLeg->AddEntry(eff18[0], "tTagger '18", "lp");
   effLeg->AddEntry(effOld16[0], "2016 analysis", "lp");

   TCanvas *can_eff[2], *can_acc[2];
   TString phaseSpace[] = {"Parton", "Particle"};
   TString years[] = {"2016 preVFP", "2016 postVFP", "2017", "2018"};

   for(int i =0; i<sizeof(eff16)/sizeof(eff16[0]); i++)
   {
	   can_eff[i] = new TCanvas(TString::Format("Efficiency can_%s%s",recoVar.Data(),phaseSpace[i].Data()), TString::Format("Efficiency can_%s%s",recoVar.Data(),phaseSpace[i].Data()), 700, 600);
	   eff18[i]->SetTitle(TString::Format("%s Efficiency '16,'17,'18 %s;%s (GeV);Efficiency",phaseSpace[i].Data(), nominal.Data(),recoVar.Data()));
	   eff18[i]->Draw();
	   eff17[i]->Draw("same");
	   eff16_pre[i]->Draw("same");
	   eff16_post[i]->Draw("same");
	   effOld16[i]->Draw("same");
	   effLeg->Draw();

		gPad->Update();
		auto graph = eff18[i]->GetPaintedGraph();
		graph->GetXaxis()->SetRangeUser(BNDmin[recoVar],BNDmax[recoVar]);
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
<<<<<<< HEAD
	   can_eff[i]->Print(TString::Format("plots%s/%s/%s/Efficiency_BTaggingSF%s_%s.pdf",nominal.Data(),binning.Data(), recoVar.Data() ,phaseSpace[i].Data(),recoVar.Data()),"pdf");
=======
	   can_eff[i]->Print(TString::Format("plots%s_UL/%s/%s/Efficiency_BTaggingSF%s_%s.pdf",nominal.Data(),binning.Data(), recoVar.Data() ,phaseSpace[i].Data(),recoVar.Data()),"pdf");
>>>>>>> 9f561e5a5dcd08bb10ed9c7b77ceb119c18806f2

	   can_acc[i] = new TCanvas(TString::Format("Acceptance can_%s%s",recoVar.Data(),phaseSpace[i].Data()), TString::Format("Acceptance can_%s%s",recoVar.Data(),phaseSpace[i].Data()), 700, 600);
	   acc18[i]->SetTitle(TString::Format("%s Acceptance '16,'17,'18 %s;%s (GeV);Acceptance",phaseSpace[i].Data(), nominal.Data(),recoVar.Data()));
	   acc18[i]->Draw();
	   acc17[i]->Draw("same");
	   acc16_pre[i]->Draw("same");
	   acc16_post[i]->Draw("same");
	   accOld16[i]->Draw("same");
       effLeg->Draw();
       gPad->Update();
	   auto graphAcc = acc18[i]->GetPaintedGraph();
	   graphAcc->SetMinimum(0.5);
	   graphAcc->SetMaximum(1.2);
	   graphAcc->GetXaxis()->SetRangeUser(BNDmin[recoVar],BNDmax[recoVar]);
	   gPad->Update();

<<<<<<< HEAD
	   can_acc[i]->Print(TString::Format("plots%s/%s/%s/Acceptance_BTaggingSF%s_%s.pdf",nominal.Data(),binning.Data(), recoVar.Data() ,phaseSpace[i].Data(),recoVar.Data()),"pdf");
=======
	   can_acc[i]->Print(TString::Format("plots%s_UL/%s/%s/Acceptance_BTaggingSF%s_%s.pdf",nominal.Data(),binning.Data(), recoVar.Data() ,phaseSpace[i].Data(),recoVar.Data()),"pdf");
>>>>>>> 9f561e5a5dcd08bb10ed9c7b77ceb119c18806f2
   }

   TH2F *hResponses[4][2]; //4 is for years , 16_preVFP-0, 16_postVFP-1, 17-2, 18-3 and 2 is for parton particle: parton-0, particle 1
   TCanvas *canResponse[4][2];
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
<<<<<<< HEAD
   		canResponse[iy][i]->Print(TString::Format("%s/%s/%s/%sResponseMatrixBtaggingSF%s_%s.pdf",years[iy].Data(),binning.Data(), recoVar.Data() ,phaseSpace[i].Data(),nominal.Data(),recoVar.Data()),"pdf");
=======
   		canResponse[iy][i]->Print(TString::Format("UL%s/%s/%s/%sResponseMatrixBtaggingSF%s_%s.pdf",years[iy].Data(),binning.Data(), recoVar.Data() ,phaseSpace[i].Data(),nominal.Data(),recoVar.Data()),"pdf");
>>>>>>> 9f561e5a5dcd08bb10ed9c7b77ceb119c18806f2
   	}
   }



}
