#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"

using std::cin;
using std::cout;
using std::endl;

void compareAngular(int ZprimeMass = 2000, bool is2016=false, bool isParton=false)
{
  gStyle->SetOptStat(0);
  TString eosPath;
  if(is2016) eosPath = "/eos/cms/store/user/gbakas/ZprimeToTT/mc/2016/";
  else eosPath = "/eos/cms/store/user/gbakas/ZprimeToTT/mc/2017/";
  
  const int sizeBins = 4;
  int sizeWidth;
  //float BND[sizeBins+1] = {1000, 2000, 3000, 4000, 5000};
  std::vector<TString> width;
  TString tempMass;
  std::vector<float> XSEC;
  std::vector<float> WEIGHT;
  float NORM;
  float LUMI;
  LUMI = 41530;
  
  if(ZprimeMass == 2000)
  {
	  sizeWidth = 3;
	  tempMass = "2TeV";
	  width = {"20", "200", "600"};
	  XSEC  = {0.1451, 0.01592, 0.005714};
  }
  else if(ZprimeMass == 3000)
  {
	  tempMass = "3TeV";
	  if(is2016)
	  {
	    sizeWidth = 2;
		width = {"30", "300"};
		XSEC  = {0.01309, 0.001809};
	  }
	  else
	  {
		sizeWidth = 3;
		width = {"30", "300", "900"};
	  }
  }
  else if(ZprimeMass == 4000)
  {
	  sizeWidth = 3;
	  tempMass = "4TeV";
	  width = {"40", "400", "1200"};
	  XSEC  = {0.001568, 0.0003247, 0.000193};
  }
  else if(ZprimeMass == 5000)
  {
	  sizeWidth = 2;
	  tempMass = "5TeV";
	  width = {"50", "500"};
	  XSEC  = {0.0002212, 8.66E-05};
  }
  else if(ZprimeMass == 2500)
  {
	  sizeWidth = 2;
	  tempMass = "2.5TeV";
	  width = {"25", "250"};
	  XSEC  = {0.04169,0.005061};
	  
  }
  TFile *inf_Zprime, *inf_TT;
  if(isParton) 
  {
	  inf_Zprime = TFile::Open(TString::Format("Output_M%s_%s.root", tempMass.Data(), "Parton"));
	  inf_TT = TFile::Open(TString::Format("Output_TT_QCD_%s.root", "Parton"));
  }
  else 
  {
	  inf_Zprime = TFile::Open(TString::Format("Output_M%s_%s.root", tempMass.Data(), "Reco"));
	  inf_TT = TFile::Open(TString::Format("Output_TT_QCD_%s.root", "Reco"));
  }
  

  float NORM_TT(0), WEIGHT_TT(0);
  float XSEC_TT = 16.74;

  TH1F *h_mTTbarParton[sizeWidth], *h_mTTbarParton_TT;
  
 
  TH1F *hAngularDistZ[sizeBins][sizeWidth], *hChiZ[sizeBins][sizeWidth];
  TH1F *hAngularDist[sizeBins], *hChi[sizeBins];
  std::vector<TString> massLimits = {"1000-2500","2500-3500", "3500-5000", "5000-Inf"};
  //std::vector<TString> legends    = {""}
  std::vector<Color_t> cols = {kRed, kBlack, kGreen};
  for(int i=0; i<sizeBins; i++)
  {
	  TString temp = massLimits[i];
	  
	  //ttbar from QCD
	  h_mTTbarParton_TT = (TH1F*)inf_TT->Get("h_mTTbarParton_TT");
	  hAngularDist[i]   = (TH1F*)inf_TT->Get(TString::Format("abs_cosTheta_%s_TT",temp.Data()));
	  hChi[i] 			= (TH1F*)inf_TT->Get(TString::Format("chi_%s_TT",temp.Data()));
	  if(i ==0)
	  {	
		  TFile *ttSelFile;
		  if(is2016) ttSelFile = TFile::Open("/eos/cms/store/user/gbakas/ttbar/topTagger/April19/TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8_Copy.root");
		  else       ttSelFile = TFile::Open("/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2017/Signal/TT_Mtt-700to1000_TuneCP5_13TeV-powheg-pythia8_Copy.root");
		  
		  NORM_TT   = ((TH1F*)ttSelFile->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
		  WEIGHT_TT = (XSEC_TT/NORM_TT)*LUMI;
	  }
	  //ttbar from Zprime for a specific mass and different widths
	  for(int wid = 0; wid<sizeWidth; wid++)
	  { 
		h_mTTbarParton[wid]   = (TH1F*)inf_Zprime->Get(TString::Format("h_mTTbarParton_M%s_W%s",tempMass.Data(), width[wid].Data() ));
		hAngularDistZ[i][wid] = (TH1F*)inf_Zprime->Get(TString::Format("abs_cosTheta_%s_M%s_W%s",temp.Data(), tempMass.Data(), width[wid].Data()));
		hChiZ[i][wid]		  = (TH1F*)inf_Zprime->Get(TString::Format("chi_%s_M%s_W%s",temp.Data(), tempMass.Data(), width[wid].Data()));
	  }
  }
    
  TFile *selFile;	
  for(int wid=0; wid<sizeWidth; wid++)
  {
  	if(is2016) selFile = TFile::Open(TString::Format(eosPath+"ZprimeToTT_M-%d_W-%s_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Copy.root", ZprimeMass, width[wid].Data())); 
	else 	   selFile = TFile::Open(TString::Format(eosPath+"ZprimeToTT_M%d_W%s_TuneCP2_13TeV-madgraphMLM-pythia8_Copy.root", ZprimeMass, width[wid].Data())); 
	NORM = ((TH1F*)selFile->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
	WEIGHT.push_back((XSEC[wid]/NORM)*LUMI);
  }
  
  TCanvas *can_theta[sizeBins], *can_chi[sizeBins], *can_mTTbarParton;
  TLegend *leg_mTTbarParton, *leg_theta[sizeBins], *leg_chi[sizeBins];
  
  can_mTTbarParton = new TCanvas("mTTbarParton", "mTTbarParton", 900, 600);
  leg_mTTbarParton = new TLegend(0.5, 0.5, 0.7, 0.7);
  h_mTTbarParton_TT->SetLineColor(kBlue);
  //h_mTTbarParton_TT->Scale(WEIGHT_TT);
  h_mTTbarParton_TT->Scale(1./h_mTTbarParton_TT->Integral());
  h_mTTbarParton_TT->Draw();
  
  for(int wid = 0; wid<sizeWidth; wid++)
  {
	//h_mTTbarParton[wid]->Scale(WEIGHT[wid]);
	h_mTTbarParton[wid]->Scale(1./h_mTTbarParton[wid]->Integral());
	h_mTTbarParton[wid]->SetLineColor(cols[wid]);
	h_mTTbarParton[wid]->Draw("same");
	leg_mTTbarParton->AddEntry(h_mTTbarParton[wid], TString::Format("Zprime, M=%s, W=%s", tempMass.Data(), width[wid].Data()), "l");

  }
  leg_mTTbarParton->AddEntry(h_mTTbarParton_TT, "TT sample", "l");
  leg_mTTbarParton->Draw();
  
  
  for(int i =0; i<sizeBins; i++)
  {
	TString temp = massLimits[i];
	can_theta[i] = new TCanvas(TString::Format("can theta %d", i),TString::Format("can theta %d", i),900, 600); 
	can_theta[i]->cd();
	leg_theta[i] = new TLegend(0.5, 0.5, 0.7, 0.7); 
    //ttbar from QCD
	hAngularDist[i]->SetLineColor(kBlue);
	//hAngularDist[i]->Scale(WEIGHT_TT);
	hAngularDist[i]->Scale(1./hAngularDist[i]->Integral());
    hAngularDist[i]->Draw();
	leg_theta[i]->AddEntry(hAngularDist[i], TString::Format("#||{cos(#theta)} %s_TT",temp.Data() ), "l");
	//ttbar from Zprime for several widths at a certain mass cos distributions
	for(int wid=0; wid<sizeWidth; wid++)
	{
		hAngularDistZ[i][wid]->SetLineColor(cols[wid]);
		//hAngularDistZ[i][wid]->Scale(WEIGHT[wid]);
		hAngularDistZ[i][wid]->Scale(1./hAngularDistZ[i][wid]->Integral());
		hAngularDistZ[i][wid]->Draw("same");
		leg_theta[i]->AddEntry(hAngularDistZ[i][wid], TString::Format("#||{cos(#theta)} %s_M%s_W%s",temp.Data(), tempMass.Data(), width[wid].Data()), "l");
	}
	leg_theta[i]->Draw();
	
	//----------------------------------------------------------------------------------------------------------------------
	
	can_chi[i] = new TCanvas(TString::Format("can #chi %d", i),TString::Format("can #chi %d", i),900, 600); 
	can_chi[i]->cd();
	leg_chi[i] = new TLegend(0.5, 0.5, 0.7, 0.7); 
	//ttbar from QCD
	hChi[i]->SetLineColor(kBlue);
	//hChi[i]->Scale(WEIGHT_TT);
	hChi[i]->Scale(1./hChi[i]->Integral());
    hChi[i]->Draw();
	leg_chi[i]->AddEntry(hAngularDist[i], TString::Format("#chi %s_TT",temp.Data() ), "l");
	//ttbar from Zprime for several widths at a certain mass for chi distributions
	for(int wid=0; wid<sizeWidth; wid++)
	{
		//hChiZ[i][wid]->Scale(WEIGHT[wid]);
		hChiZ[i][wid]->Scale(1./hChiZ[i][wid]->Integral());
		hChiZ[i][wid]->SetLineColor(cols[wid]); 
		hChiZ[i][wid]->Draw("same");
		leg_chi[i]->AddEntry(hAngularDistZ[i][wid], TString::Format("#chi %s_M%s_W%s",temp.Data(), tempMass.Data(), width[wid].Data()), "l");
	}
	leg_chi[i]->Draw(); 
  }

	
}
