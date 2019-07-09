#include <iostream>
#include <TString.h>
#include <TH1F.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include "TMath.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//This code is made to make some checks to the spin correlation test file.
//The file includes only gen level particles (at parton level) from MC and we receive only their Id's
//We know that top has id of +/- 6  and the decay products are W+ 24, W- -24, b +5, etc
//Here I will check the number of quark pairs as decay products of W bosons vs the number of lepton and neutrino pairs coming from W bosons
//Also I will plot the Particle Mass for every Particle Id that exists in the root file.
//
//
//Written by: Georgios Bakas
//19/11/2018
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

using namespace std;
void boostLorentzVectorTTbar(TLorentzVector p4_1, TLorentzVector p4_2, TLorentzVector &p4New_1, TLorentzVector &p4New_2);
void boostLorentzVector(TLorentzVector p4, TLorentzVector &p4_New, TVector3 boostVec);

void LorentzVector_checks()
{
  gStyle->SetOptStat(0);
  //TFile *f = TFile::Open("/afs/cern.ch/user/g/gbakas/CRAB3-jobs2016/CMSSW_8_0_6/src/UserCode/TopAnalysis/test/flatTree_10KNew.root");
  TFile *f = TFile::Open("/afs/cern.ch/work/g/gbakas/private/TTSpinCorrelation/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");
  TTree *t = (TTree*)f->Get("ttbar/events");
  cout<<"opened tree"<<endl;
  const char* momentum[4] = {"P_{x}", "P_{y}", "P_{z}", "E"};
  
  TH1F *h_PR_CM = new TH1F("#Delta(Pseudorapity_{top} - Pseudorapity_{tbar}) in CM frame", "#Delta(Pseudorapity_{top} - Pseudorapity_{tbar}) in CM frame", 50,-5,5);
  TH1F *h_PR_Lab = new TH1F("#Delta(Pseudorapity_{top} - Pseudorapity_{tbar}) in Lab frame", "#Delta(Pseudorapity_{top} - Pseudorapity_{tbar}) in Lab frame", 50,-5,5);
  
  TH1F *h_PR_CM_Wb = new TH1F("#Delta(Pseudorapity_{W^{+}} - Pseudorapity_{b}) in CM frame", "#Delta(Pseudorapity_{W^{+}} - Pseudorapity_{b}) in top rest frame", 50,-5,5);
  TH1F *h_PR_Lab_Wb = new TH1F("#Delta(Pseudorapity_{W^{+}} - Pseudorapity_{b}) in Lab frame", "#Delta(Pseudorapity_{W^{+}} - Pseudorapity_{b}) in Lab Frame", 50,-5,5);
  
  TH1F *h_PR_CM_Wbbar = new TH1F("#Delta(Pseudorapity_{W^{-}} - Pseudorapity_{#bar{b}}) in CM frame", "#Delta(Pseudorapity_{W^{-}} - Pseudorapity_{#bar{b}}) in tbar rest frame", 50,-5,5);
  TH1F *h_PR_Lab_Wbbar = new TH1F("#Delta(Pseudorapity_{W^{-}} - Pseudorapity_{#bar{b}}) in Lab frame", "#Delta(Pseudorapity_{W^{-}} - Pseudorapity_{#bar{b}}) in Lab frame", 50,-5,5);
  
  TH1F *h_4momentum[4];
  h_4momentum[0] = new TH1F("#Delta p_{x} Top quark - p_{x} (W^{+} + b)", "#Delta p_{x} Top quark - p_{x} (W^{+} + b)", 50, -0.05,0.05);
  h_4momentum[1] = new TH1F("#Delta p_{y} Top quark - p_{y} (W^{+} + b)", "#Delta p_{y} Top quark - p_{y} (W^{+} + b)", 50, -0.05,0.05);
  h_4momentum[2] = new TH1F("#Delta p_{z} Top quark - p_{z} (W^{+} + b)", "#Delta p_{z} Top quark - p_{z} (W^{+} + b)", 50, -0.05,0.05);
  h_4momentum[3] = new TH1F("#Delta E_{Top quark} - E_{W^{+} + b}", "#Delta E_{Top quark} - E_{W^{+} + b}", 50, -0.05,0.05);
 //get Branches and their values, here there are all vector floats
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
  //cout<<numOfEntries<<endl;

  //define the Lorentz vectors
  TLorentzVector p4T(0,0,0,0), p4Tbar(0,0,0,0),p4W(0,0,0,0), p4Wplus(0,0,0,0), p4Wminus(0,0,0,0), p4b(0,0,0,0), p4bbar(0,0,0,0);
  TLorentzVector p4T_TTbarBoost(0,0,0,0), p4Tbar_TTbarBoost(0,0,0,0), p4Wminus_topBoost(0,0,0,0), p4Wplus_topBoost(0,0,0,0), p4b_topBoost(0,0,0,0), p4bbar_topBoost(0,0,0,0);

  //loop over all entries
  for(int k =0; k<=numOfEntries; k++)
  {
    
    t->GetEntry(k);
	for(int l=0; l<=5; l++)
	{
			//verify the particle Id
			if( (*particleId)[l] == 6 || (*particleId)[l] == -6)
			{
				if((*particleId)[l] == 6)
				{
					p4T.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
				}
				else 
				{
					p4Tbar.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
				}
				
				
				//Printf("(P,eta,phi,E)=(%f,%f,%f,%f)",  p4T.P(),p4T.Eta(),p4T.Phi(),p4T.E());
			}
			//this is W
			if( (*particleId)[l] == 24 || (*particleId)[l] == -24)
			{
				
				p4W.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
				if((*particleId)[l] == 24) //w+
				{	
					//get the boost from top
					p4Wplus.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
					boostLorentzVector(p4Wplus, p4Wplus_topBoost, p4T.BoostVector());
				}
				else //w-
				{
					//get the boost from tbar
					p4Wminus.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
					boostLorentzVector(p4Wminus, p4Wminus_topBoost, p4Tbar.BoostVector());

				}
				
				
			}
			//this is b, b bar from top
			else if (fabs((*particleId)[l]) == 5)
			{
				//p4b.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
				
				if((*particleId)[l]==5) //b coming from top
				{
					p4b.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
					boostLorentzVector(p4b, p4b_topBoost, p4T.BoostVector());
				}
				else //bbar coming from tbar
				{
					p4bbar.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
					boostLorentzVector(p4bbar, p4bbar_topBoost, p4Tbar.BoostVector());
				}
			}
			
			
			
		 
		}
		
		boostLorentzVectorTTbar(p4T, p4Tbar, p4T_TTbarBoost, p4Tbar_TTbarBoost);
		h_PR_CM->Fill(p4T_TTbarBoost.PseudoRapidity() - p4Tbar_TTbarBoost.PseudoRapidity());
		h_PR_Lab->Fill(p4T.PseudoRapidity() - p4Tbar.PseudoRapidity());

	
		//these are the histograms that show the PseudoRapidity of the W and b in the top system and the antitop system respectivelly
		h_PR_CM_Wb->Fill(p4Wplus_topBoost.PseudoRapidity()-p4b_topBoost.PseudoRapidity());
		h_PR_Lab_Wb->Fill(p4Wplus.PseudoRapidity()-p4b.PseudoRapidity());
			
		h_PR_CM_Wbbar->Fill(p4Wminus_topBoost.PseudoRapidity()-p4bbar_topBoost.PseudoRapidity());
		h_PR_Lab_Wbbar->Fill(p4Wminus.PseudoRapidity()-p4bbar.PseudoRapidity());
		
	
	//check that the p4T = p4W + p4b
	TLorentzVector tempVector;
	tempVector.SetPxPyPzE(p4T.Px()-p4Wplus.Px()-p4b.Px(),p4T.Py()-p4Wplus.Py()-p4b.Py(), p4T.Pz()-p4Wplus.Pz()-p4b.Pz(), p4T.Energy()-p4Wplus.Energy()-p4b.Energy());
	//tempVector.Print();
	h_4momentum[0]->Fill(p4T.Px()-p4Wplus.Px()-p4b.Px());
	h_4momentum[1]->Fill(p4T.Py()-p4Wplus.Py()-p4b.Py());
	h_4momentum[2]->Fill(p4T.Pz()-p4Wplus.Pz()-p4b.Pz());	
	h_4momentum[3]->Fill(p4T.Energy()-p4Wplus.Energy()-p4b.Energy());
  }


  TFile *outputFile =  TFile::Open("LorentzChecks.root", "RECREATE");
  TLegend *leg1 = new TLegend(0.9,0.7,0.6,0.6);
  outputFile->cd();
  TCanvas *can1 =  new TCanvas("canvas 1","canvas 1", 900,600);
  h_PR_CM->Draw();
  h_PR_CM->GetXaxis()->SetTitle("PseudoRapidity");
  h_PR_CM->GetYaxis()->SetTitle("Events");
  h_PR_Lab->GetXaxis()->SetTitle("PseudoRapidity");
  h_PR_Lab->GetYaxis()->SetTitle("Events");
  h_PR_CM->SetLineColor(kRed);
  leg1->AddEntry(h_PR_CM, "PseudoRapidity in ttbar CM system");
  leg1->AddEntry(h_PR_Lab, "PseudoRapidity in ttbar Lab system");
  h_PR_Lab->Draw("same");
  leg1->Draw();
  
  TCanvas *can2 = new TCanvas("canvas 2","canvas 2", 900,600);
  TLegend *leg2 = new TLegend(0.9,0.7,0.6,0.6);
  h_PR_CM_Wb->GetXaxis()->SetTitle("PseudoRapidity");
  h_PR_CM_Wb->GetYaxis()->SetTitle("Events");
  h_PR_Lab_Wb->GetXaxis()->SetTitle("PseudoRapidity");
  h_PR_Lab_Wb->GetYaxis()->SetTitle("Events");
  h_PR_CM_Wb->Draw();
  h_PR_CM_Wb->SetLineColor(kRed);
  h_PR_Lab_Wb->Draw("same");
  leg2->AddEntry(h_PR_CM_Wb, "PseudoRapidity in Wb CM system");
  leg2->AddEntry(h_PR_Lab_Wb, "PseudoRapidity in Wb Lab system");
  leg2->Draw();
  
  TCanvas *can3 = new TCanvas("canvas 3","canvas 3", 900,600);
  TLegend *leg3 = new TLegend(0.9,0.7,0.6,0.6);
  h_PR_CM_Wbbar->GetXaxis()->SetTitle("PseudoRapidity");
  h_PR_CM_Wbbar->GetYaxis()->SetTitle("Events");
  h_PR_Lab_Wbbar->GetXaxis()->SetTitle("PseudoRapidity");
  h_PR_Lab_Wbbar->GetYaxis()->SetTitle("Events");
  h_PR_CM_Wbbar->Draw();
  h_PR_CM_Wbbar->SetLineColor(kRed);
  h_PR_Lab_Wbbar->Draw("same");
  leg3->AddEntry(h_PR_CM_Wb, "PseudoRapidity in Wbbar CM system");
  leg3->AddEntry(h_PR_Lab_Wb, "PseudoRapidity in Wbbar Lab system");
  leg3->Draw();
  
  //outputFile->Close();
  //f->Close();
  TString tempString;
  for(int i =0; i<=3; i++)
  {
	  tempString = momentum[i];
	  h_4momentum[i]->GetXaxis()->SetTitle(tempString+" [GeV]");
	  h_4momentum[i]->GetYaxis()->SetTitle("Events");
	  h_4momentum[i]->Write();
  }
  h_PR_CM_Wb->Write();
  h_PR_CM_Wbbar->Write();
  h_PR_Lab_Wb->Write();
  h_PR_Lab_Wbbar->Write();
  h_PR_CM->Write();
  h_PR_Lab->Write();
  

}


void boostLorentzVectorTTbar(TLorentzVector p4_1, TLorentzVector p4_2, TLorentzVector &p4New_1, TLorentzVector &p4New_2)
{
	
	//define the combined Lorentz vector of ttbar 
	TLorentzVector p4CombinedVector;
	p4CombinedVector.SetPxPyPzE(p4_1.Px()+p4_2.Px(),p4_1.Py()+p4_2.Py(), p4_1.Pz()+p4_2.Pz(), p4_1.Energy()+p4_2.Energy()); 
	//get boost from this vector
	TVector3 TTbar_boostVector = p4CombinedVector.BoostVector();
	p4New_1.Boost(-TTbar_boostVector);
	p4New_2.Boost(-TTbar_boostVector);
	//p4CombinedVector.Print();
	
	//check if it is correct by boosting the combined also and check if PT is 0 there!
	TLorentzVector p4CombinedVector_NewSys = p4CombinedVector;
	p4CombinedVector_NewSys.Boost(-TTbar_boostVector);
	cout<<"after boost"<<endl;
	p4CombinedVector_NewSys.Print();

}

void boostLorentzVector(TLorentzVector p4, TLorentzVector &p4_New, TVector3 boostVec)
{
	p4_New = p4;
	p4_New.Boost(boostVec);
	
}
