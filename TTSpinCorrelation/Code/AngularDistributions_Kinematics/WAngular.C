#include <TTree.h>
#include <TFile.h>
#include <vector>
#include <string>
#include <iostream>
#include <iterator>
#include <TH1F.h>
#include <TString.h>
#include <TFile.h>
#include "TLorentzVector.h"
#include <TLatex.h>
#include <TCanvas.h>
#include "TStyle.h"
#include <TLegend.h>


TVector3 getBoostVector(TLorentzVector p4_1, TLorentzVector p4_2, TLorentzVector &p4CombinedVector);

void WAngular()
{
	

  gStyle->SetOptStat(0);
  //reco angular distributions
  TString inputfile = "/afs/cern.ch/work/g/gbakas/public/NewMCSamples/flatTree_TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root";
  
  TFile *fileData2016 = new TFile(inputfile, "READ"); 
  TTree *treeData2016 = (TTree*)fileData2016->Get("boosted/events");
cout<<"Read File"<<endl;

  int category = 0, nJets = 0;
  int nvtx = 0;
  float mva,  ht;
  std::vector<float> *jetMassSub0 = 0, *jetMassSub1 = 0;
  std::vector<float> *jetEtaSub0 = 0, *jetEtaSub1 = 0;
  std::vector<float> *jetPtSub0 = 0, *jetPtSub1 = 0;
  std::vector<float> *jetPhiSub0 = 0, *jetPhiSub1 = 0;

  std::vector<bool> *triggerBits = 0;
  std::vector<float> *jetPt = 0;
  std::vector<float> *jetEta = 0;
  std::vector<float> *jetPhi = 0;
  std::vector<float> *jetMass = 0;
  std::vector<float> *jetBtagSub0 = 0, *jetBtagSub1 = 0;

  std::vector<float> *jetMassSoftDrop = 0;
  int triggerBit = 2;

  float binValuesJetMassSub[] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100,
                                  105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200,
                                  205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 300};
  int nbinJetMassSub = sizeof(binValuesJetMassSub)/sizeof(float) - 1; 
  std::cout<<"nbinJetMassSub: "<<nbinJetMassSub<<std::endl;
  
  treeData2016->SetBranchAddress("triggerBit", &triggerBits);
  treeData2016->SetBranchAddress("category", &category);
  treeData2016->SetBranchAddress("nJets", &nJets);
  treeData2016->SetBranchAddress("ht", &ht);
  treeData2016->SetBranchAddress("mva", &mva);
  treeData2016->SetBranchAddress("jetBtagSub0", &jetBtagSub0);
  treeData2016->SetBranchAddress("jetBtagSub1", &jetBtagSub1);
  treeData2016->SetBranchAddress("jetMassSoftDrop", &jetMassSoftDrop);
  treeData2016->SetBranchAddress("nvtx", &nvtx);
  
  treeData2016->SetBranchAddress("jetPt", &jetPt);
  treeData2016->SetBranchAddress("jetEta", &jetEta);
  treeData2016->SetBranchAddress("jetPhi", &jetPhi);
  treeData2016->SetBranchAddress("jetMass", &jetMass);
  
  treeData2016->SetBranchAddress("jetMassSub0", &jetMassSub0);
  treeData2016->SetBranchAddress("jetMassSub1", &jetMassSub1);
  treeData2016->SetBranchAddress("jetPtSub0", &jetPtSub0);
  treeData2016->SetBranchAddress("jetPtSub1", &jetPtSub1);
  treeData2016->SetBranchAddress("jetPhiSub1", &jetPhiSub1);
  treeData2016->SetBranchAddress("jetPhiSub1", &jetPhiSub1);
  treeData2016->SetBranchAddress("jetEtaSub0", &jetEtaSub0);
  treeData2016->SetBranchAddress("jetEtaSub1", &jetEtaSub1);

  

  TH1F* dPhi_Reco = new TH1F("W candidate #Delta#phi Reco" , "W candidate #Delta#phi Reco", 50, -3, 3); 
  TH1F* dPhi_Parton = new TH1F("W #Delta#phi Parton" , "W #Delta#phi Parton", 50, -3, 3); 

  TH1F* h_PtParton = new TH1F("P_{T} parton","P_{T} parton", 50, 350,1000);
  TH1F* h_PhiParton = new TH1F("#phi parton","#phi parton", 50, -3,3);
  TH1F* h_EtaParton = new TH1F("#eta parton","#eta parton", 50, 0,3);
  TH1F* h_MassParton = new TH1F("Mass parton","Mass parton", 50, 0,500);
  
  
  TH1F* h_PtReco = new TH1F("P_{T} reco","P_{T} reco", 50, 350,1000);
  TH1F* h_PhiReco = new TH1F("#phi reco","#phi reco", 50, -3,3);
  TH1F* h_EtaReco = new TH1F("#eta reco","#eta reco", 50, 0,3);
  TH1F* h_MassReco = new TH1F("Mass reco","Mass reco", 50, 0,500);
  
  TH1F* h_massSoftDrop = new TH1F("Jet Mass Soft Drop", "Jet Mass Soft Drop", 50, 0, 1500);
  TH1F* h_mass = new TH1F("Jet Mass", "Jet Mass", 50, 0, 1500);
										
  TLorentzVector p4T_0(0,0,0,0), p4T_1(0,0,0,0), p4W_0(0,0,0,0), p4W_1(0,0,0,0);
  TLorentzVector p4TTbar(0,0,0,0);
  cout<<"starting..."<<endl;
  for(int i=0; i<treeData2016->GetEntries(); i++)
  {
    treeData2016->GetEntry(i);
	
	for(int k=0; k<=1; k++)
	{
		h_mass->Fill((*jetMass)[k]);
		h_massSoftDrop->Fill((*jetMassSoftDrop)[k]);
	}
    if((*triggerBits)[triggerBit] && category == 2  && nJets > 1 && mva > 0.8 && nvtx >= 15 && nvtx <= 25)
    {
		
		//here get the top variables eta, phi, mass, pT to reconstruct the top 4-Vector
		p4T_0.SetPtEtaPhiM((*jetPt)[0], (*jetEta)[0], (*jetPhi)[0], (*jetMass)[0]);
		p4T_1.SetPtEtaPhiM((*jetPt)[1], (*jetEta)[1], (*jetPhi)[1], (*jetMass)[1]);
		if ((*jetPt)[0] > 350)
		{
			h_PtReco->Fill((*jetPt)[0]);
			h_PhiReco->Fill((*jetPhi)[0]);
			h_EtaReco->Fill((*jetEta)[0]);
			h_MassReco->Fill((*jetMassSoftDrop)[0]);
		}
		if ((*jetPt)[1] > 350) 
		{
			h_PtReco->Fill((*jetPt)[1]);
			h_PhiReco->Fill((*jetPhi)[1]);
			h_EtaReco->Fill((*jetEta)[1]);
			h_MassReco->Fill((*jetMassSoftDrop)[1]);
			
		
			//jetMassSubHisto_2016->Fill((*jetMass)[0]);
			//jetMassSubHisto_2016->Fill((*jetMass)[1]);
			
			//I believe that this is going to need identification meaning that you must check that this is actually a Top quark so I need access to the parton quantity
			//
			//
		  TVector3 TTbar_boostVector = getBoostVector( p4T_0, p4T_1, p4TTbar);
		  p4TTbar.Print();
		  int bothW = 0;
		  
		  for(std::vector<float>::iterator ptIter = jetPt->begin(); ptIter!=jetPt->end(); ++ptIter)
		  {
				 bothW = 0;
				int index = ptIter - jetPt->begin();
				if((*jetMassSoftDrop)[index] >= 150 && (*jetMassSoftDrop)[index] <= 210 && *ptIter > 450)
				{
				  if((*jetPtSub0)[index] > 400)
				  {
					if((*jetBtagSub0)[index] < 0.5426)
					{
						p4W_0.SetPtEtaPhiM((*jetPt)[index], (*jetEta)[index], (*jetPhi)[index], (*jetMassSub0)[index]);
						bothW++;
					}
				  }
				  if((*jetPtSub1)[index] > 400)
				  {
					 if((*jetBtagSub1)[index] < 0.5426)
					{   
						p4W_1.SetPtEtaPhiM((*jetPt)[index], (*jetEta)[index], (*jetPhi)[index], (*jetMassSub0)[index]);  
						bothW++;
					}
				  }
					//this means if found both W's 
					//if(bothW ==2)
					//{
						dPhi_Reco->Fill(p4W_0.Phi()-p4W_1.Phi());
					//}
				}
				
				
			
			}
	  
		}
    }
  }
  
 
 
 //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 //start with the parton W angular diistributions
  TFile *f = TFile::Open("/afs/cern.ch/work/g/gbakas/private/TTSpinCorrelation/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");
  //TFile *f = TFile::Open("/afs/cern.ch/user/g/gbakas/CRAB3-jobs2016/CMSSW_8_0_6/src/UserCode/TopAnalysis/test/flatTree_10K-01-12.root");
  TTree *t = (TTree*)f->Get("ttbar/events");
  vector <float>*particleMass(0);
  vector <float>*particlePt(0);
  vector <float>*particleEta(0);
  vector <float>*particlePhi(0);
  //this is particle id
  vector <int>  *particleId(0);

  //declare the branches
  TBranch *b_eta, *b_pt, *b_phi, *b_mass, *b_id;
  t->SetBranchAddress("particleId",&particleId, &b_id);
  t->SetBranchAddress("particleMass",&particleMass, &b_mass);
  t->SetBranchAddress("particlePt",&particlePt, &b_pt);
  t->SetBranchAddress("particleEta",&particleEta, &b_eta);
  t->SetBranchAddress("particlePhi",&particlePhi, &b_phi);
  long numOfEntries = t->GetEntries();

  TLorentzVector p4T_parton(0,0,0,0), p4Tbar_parton(0,0,0,0), p4Wminus_parton(0,0,0,0), p4Wplus_parton(0,0,0,0);
  bool isTopPtOk = 0;
  cout<<"Starting with parton level...."<<endl;
  for(int k =0; k<=numOfEntries; k++)
  {
	  
	  t->GetEntry(k);
	  for(int l=0; l<=5; l++)
	  {
	  //verify the particle Id
			if( (*particleId)[l] == 6 || (*particleId)[l] == -6)
			{
				
				if((*particlePt) [l] > 350)
				{
					isTopPtOk = 1;
					h_PtParton->Fill((*particlePt)[l]);
					h_PhiParton->Fill((*particlePhi)[l]);
					h_EtaParton->Fill((*particleEta)[l]);
					h_MassParton->Fill((*particleMass)[l]);
					//cout<<(*particleMass)[l]<<endl;

					if((*particleId)[l] == 6)
					{
						p4T_parton.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
						
					}
					else 
					{
						p4Tbar_parton.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
					}
				}
				
			}
			//this is W
			if( ((*particleId)[l] == 24 || (*particleId)[l] == -24) && isTopPtOk)
			{
				if((*particleMass)[l] <=210 && (*particleMass)[l] > 150 && (*particlePt)[l] > 400)
				{
					if((*particleId)[l] == 24) 
					{	
						//get the boost from top
						p4Wplus_parton.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
					
					}
					else 
					{
						//get the boost from tbar
						p4Wminus_parton.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);

					}	
				}
			}
	   } 
		
		
		

		TLorentzVector p4TTbar_parton(0,0,0,0);
		TVector3 partonTTbarBoostVector	= getBoostVector(p4T_parton, p4Tbar_parton, p4TTbar_parton);
		//p4TTbar_parton.Print();
		//p4Wminus_parton.Boost(partonTTbarBoostVector);
		//p4Wplus_parton.Boost(partonTTbarBoostVector);
		
		dPhi_Parton->Fill(fabs(p4Wminus_parton.Phi()- p4Wplus_parton.Phi()));
		
  }
  
  //TCanvas *can2 = new TCanvas("test", "test", 900, 600);
 // dPhi_Reco->Draw();
  
 double norm =1;
 Double_t scale_1 = norm/(h_PtParton->Integral());
 h_PtParton->Scale(scale_1);
 Double_t scale_2 = norm/(h_PtReco->Integral());
 h_PtReco->Scale(scale_2);


 h_PtParton->GetXaxis()->SetTitle("P_{T} (GeV)");
 h_PtParton->GetYaxis()->SetTitle("Events");
 h_PtReco->GetXaxis()->SetTitle("P_{T} (GeV)");
 h_PtReco->GetYaxis()->SetTitle("Events");
 
 
  TCanvas *can3 = new TCanvas("Pt", "Pt", 900,600);
  TLegend *leg_pt = new TLegend(0.86,0.7,0.99,0.9);
  h_PtParton->SetLineColor(kRed);
  h_PtReco->SetLineColor(kBlack);
  leg_pt->AddEntry(h_PtParton, "h_PtParton", "l");
  leg_pt->AddEntry(h_PtReco, "h_PtReco", "l"); 
  h_PtParton->Draw();
  h_PtReco->Draw("same");
  leg_pt->Draw();
  
  
  //just checking difference between jetMass and jetMassSoftDrop
  TCanvas *canMass1 = new TCanvas("mass", "mass", 900,600);
  TLegend *leg = new TLegend(0.86,0.7,0.99,0.9); 
  h_mass->SetLineColor(kRed);
  h_massSoftDrop->SetLineColor(kBlue);
  h_mass->GetXaxis()->SetTitle("Mass (GeV)");
  h_massSoftDrop->GetXaxis()->SetTitle("Mass (GeV)");
  h_mass->GetYaxis()->SetTitle("Events");
  h_massSoftDrop->GetYaxis()->SetTitle("Events");
  leg->AddEntry(h_mass,"Mass","l");
  leg->AddEntry(h_massSoftDrop,"Mass SoftDrop","l"); 
  h_mass->Draw();
  h_massSoftDrop->Draw("same");
  leg->Draw();
  
  
  //phi distributions for reco and parton
  Double_t scale_phi1 = norm/(h_PhiParton->Integral());
  h_PhiParton->Scale(scale_phi1);
  Double_t scale_phi2 = norm/(h_PhiReco->Integral());
  h_PhiReco->Scale(scale_phi2);
  
  
  h_PhiParton->SetLineColor(kRed);
  h_PhiReco->SetLineColor(kBlue);
  h_PhiParton->GetXaxis()->SetTitle("#phi");
  h_PhiReco->GetXaxis()->SetTitle("#phi");
  
  TCanvas *canPhi = new TCanvas("#phi can", "#phi can", 900, 600);
  TLegend *legPhi = new TLegend(0.86,0.7,0.99,0.9);
  h_PhiParton->Draw();
  h_PhiReco->Draw("same");
  legPhi->AddEntry(h_PhiParton,"#phi parton","l");
  legPhi->AddEntry(h_PhiReco,"#phi reco","l"); 
  legPhi->Draw();
  
  //eta distributions for reco and parton
  Double_t scale_eta1 = norm/(h_EtaParton->Integral());
  h_EtaParton->Scale(scale_eta1);
  Double_t scale_eta2 = norm/(h_EtaReco->Integral());
  h_EtaReco->Scale(scale_eta2);
  
  h_EtaParton->SetLineColor(kRed);
  h_EtaReco->SetLineColor(kBlue);
  h_EtaParton->GetXaxis()->SetTitle("#eta");
  h_EtaReco->GetXaxis()->SetTitle("#eta");
  
  TCanvas *canEta = new TCanvas("#eta can", "#eta can", 900, 600);
  TLegend *legEta = new TLegend(0.86,0.7,0.99,0.9);
  h_EtaParton->Draw();
  h_EtaReco->Draw("same");
  legEta->AddEntry(h_EtaParton,"#eta parton","l");
  legEta->AddEntry(h_EtaReco,"#eta reco","l"); 
  legEta->Draw();
  
  //mass distributions for reco and parton
  Double_t scale_mass1 = norm/(h_MassParton->Integral());
  h_MassParton->Scale(scale_mass1);
  Double_t scale_mass2 = norm/(h_MassReco->Integral());
  h_MassReco->Scale(scale_mass2);
  
  h_MassParton->SetLineColor(kRed);
  h_MassReco->SetLineColor(kBlue);
  h_MassParton->GetXaxis()->SetTitle("Mass (GeV)");
  h_MassReco->GetXaxis()->SetTitle("Mass (GeV)");
  
  
  TCanvas *canMass = new TCanvas("Mass can", "Mass can", 900, 600);
  h_MassParton->Draw();
  h_MassReco->Draw("same");
  TLegend *legMass = new TLegend(0.86,0.7,0.99,0.9);
  h_MassParton->Draw();
  h_MassReco->Draw("same");
  legMass->AddEntry(h_MassParton,"Mass parton","l");
  legMass->AddEntry(h_MassReco,"Mass reco","l"); 
  legMass->Draw();
 
  
  
  //for(auto & lep: myLeptons) if( deltaR(lep->eta(),lep->phi(),ijet->eta(),ijet->phi()) < DRmax ) isLeptonMatched = true;
  //this is lepton matching and it is something you don't want
  //you need to match the top, maybe use b-tagging technique..?? like taking dR of W candidate and b candidate?
  
  /*
 double norm =1;
 Double_t scale_1 = norm/(dPhi_Parton->Integral());
 dPhi_Parton->Scale(scale_1);
 Double_t scale_2 = norm/(dPhi_Reco->Integral());
 dPhi_Reco->Scale(scale_2);
 
 dPhi_Reco->GetXaxis()->SetTitle("#Delta#Phi");
 dPhi_Parton->GetXaxis()->SetTitle("#Delta#Phi");
 dPhi_Reco->GetYaxis()->SetTitle("Norm Events");
 dPhi_Parton->GetYaxis()->SetTitle("Norm Events");
 
 
 
 TCanvas *can1 = new TCanvas("#DeltaPhi Canvas", "#DeltaPhi Canvas", 900,600);
 dPhi_Reco->SetLineColor(kBlack);
 dPhi_Parton->SetLineColor(kRed);
 dPhi_Reco->Draw();
 dPhi_Parton->Draw("same");
 */
 
 
}
TVector3 getBoostVector(TLorentzVector p4_1, TLorentzVector p4_2, TLorentzVector &p4CombinedVector)
{
	//define the combined Lorentz vector of ttbar 
	//TLorentzVector p4CombinedVector;
	p4CombinedVector.SetPxPyPzE(p4_1.Px()+p4_2.Px(),p4_1.Py()+p4_2.Py(), p4_1.Pz()+p4_2.Pz(), p4_1.Energy()+p4_2.Energy()); 
	//get boost from this vector
	TVector3 TTbar_boostVector = p4CombinedVector.BoostVector();
	p4CombinedVector.Boost(-TTbar_boostVector);
	return -TTbar_boostVector;
}

