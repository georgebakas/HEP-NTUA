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
#include "THStack.h"
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//This code will check the angular distributions of the W decay products
//It will check the angular distributions for the semileptonic decay as well as the fully hadronic decay
//
//
//Written by: Georgios Bakas
//29/11/2018
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

using namespace std;
void boostLorentzVectorTTbar(TLorentzVector p4_1, TLorentzVector &p4_2, TVector3 boostVec);
TVector3 getBoostVectorTTbar(TLorentzVector p4_1, TLorentzVector p4_2);
Float_t helAngle(TLorentzVector particle, TLorentzVector parent, TLorentzVector grandparent);

void Semileptonic_AngularDistributions(int bins = 20)
{
  //gStyle->SetOptStat(0);
  
  //TFile *f = TFile::Open("/afs/cern.ch/work/g/gbakas/private/TTSpinCorrelation/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");
  TFile *f = TFile::Open("/afs/cern.ch/work/g/gbakas/private/TTSpinCorrelation/boostedTTBar_SpCorr/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_BoostedPart.root");
  TTree *t = (TTree*)f->Get("boosted/events");
  cout<<"opened tree"<<endl;
  
  
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
   
  
  
  //this is semileptonic cos(theta) where theta is the angle between the d-quark and the charged lepton in the ttbar CM frame
  TH1F *h_semilepton = new TH1F ("Semileptonic cos#theta Distribution", "Semileptonic cos#theta Distribution", 20, -1,1); 
  TH1F *h_semileptonPhi[2], *h_HelicityAngle[2];
  h_semileptonPhi[0] = new TH1F ("Semileptonic #Delta#phi Distribution", "Semileptonic #Delta#phi Distribution", bins, -3,3); 
  h_semileptonPhi[1] = new TH1F ("Semileptonic #Sigma#phi Distribution", "Semileptonic #Sigma#phi Distribution", bins, -3,3); 
  
  h_HelicityAngle[0] = new TH1F("Angle distribtution cos(#theta_{l,d})", "Angle distribtution cos(#theta_{l,d})", bins, -1,1); 
  h_HelicityAngle[1] = new TH1F("Helicity angle distribtution cos(#theta^{*}_{l})cos(#theta_{d}))", "Helicity angle distribtution cos(#theta^{*}_{l})cos(#theta_{d})",bins, -1,1); 
  
  TLorentzVector p4T(0,0,0,0), p4Tbar(0,0,0,0),p4W(0,0,0,0), p4Wplus(0,0,0,0), p4Wminus(0,0,0,0), p4b(0,0,0,0), p4bbar(0,0,0,0), p4Down(0,0,0,0), p4DownBar(0,0,0,0);
  TLorentzVector p4T_TTbarBoost(0,0,0,0), p4Tbar_TTbarBoost(0,0,0,0), p4Wminus_boosted(0,0,0,0), p4Wplus_boosted(0,0,0,0), p4b_boosted(0,0,0,0), p4bbar_boosted(0,0,0,0);
  TLorentzVector p4LeptonSemi(0,0,0,0), p4LeptonSemi_boosted(0,0,0,0), p4Down_boosted(0,0,0,0), p4DownBar_boosted(0,0,0,0); 
 
  
  int *foundQuark, *foundLepton, *foundLepton_2, *foundLepton_3, *foundQuark_2, *foundQuark_3;
  int quarkIds[] = {1,2,3,4,5};
  int leptonIds[] = {11,12,13,14,15,16};
  int semileptonic = 0;
  int fullyHadronic = 0;
  int dilepton = 0;
  int dileptonMassCut = 0;
  long numOfEntries = t->GetEntries();
  
  Float_t angle_quark, angle_lepton, angle_1;
  //cout<<numOfEntries<<endl;
  int semileptonic_quarkFirst = 0;
  int semileptonic_leptonFirst = 0;
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
				
			}
			//this is W
			if( (*particleId)[l] == 24 || (*particleId)[l] == -24)
			{
				
				if((*particleId)[l] == 24) 
				{	
					//get the boost from top
					p4Wplus.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
					
				}
				else 
				{
					//get the boost from tbar
					p4Wminus.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);

				}

			}
			//this is b, b bar from top
			else if (fabs((*particleId)[l]) == 5)
			{
				
				if((*particleId)[l]==5) // coming from top
				{
					p4b.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
				}
				else //bbar coming from tbar
				{
					p4bbar.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
				}
			}
		
		 
	}
	
	//get the boost vector and boost everything with this vector
		TVector3 ttbarBoostVector = getBoostVectorTTbar(p4T,p4Tbar);
		boostLorentzVectorTTbar(p4T, p4T_TTbarBoost,ttbarBoostVector );
		boostLorentzVectorTTbar(p4Tbar, p4Tbar_TTbarBoost,ttbarBoostVector );

		boostLorentzVectorTTbar(p4Wplus,p4Wplus_boosted, ttbarBoostVector);
		boostLorentzVectorTTbar(p4Wminus,p4Wminus_boosted, ttbarBoostVector);
		boostLorentzVectorTTbar(p4b,p4b_boosted, ttbarBoostVector);
		boostLorentzVectorTTbar(p4bbar,p4bbar_boosted, ttbarBoostVector);

		foundQuark_2 = std::find(quarkIds,quarkIds+5, fabs((*particleId)[6])); 
		foundLepton_2 = std::find(leptonIds,leptonIds+6, fabs((*particleId)[6]));
		
		TVector3 tempVector;
		bool foundDownQuark;
		if(foundQuark_2 != quarkIds+5) //at least semilepton
		{
			int choiceDown=0;
			for(int l=6; l<=7; l++)
			{
				if(fabs((*particleId)[l] == 1))
				{
					foundDownQuark = 1;
					if ((*particleId)[l] == 1)p4Down.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
					else p4DownBar.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
					choiceDown = l;
				}
			}
	
		    for(int r=8; r<=9; r++)
			{
				
				foundLepton_3 = std::find(leptonIds,leptonIds+6, fabs((*particleId)[r]));
				if(foundLepton_3 !=leptonIds+6 && foundDownQuark == 1 && (fabs((*particleId)[r]) == 11 || fabs((*particleId)[r]) == 13 || fabs((*particleId)[r]) == 15))
				{
					    //semileptonic
						semileptonic_quarkFirst++;
						semileptonic++;
						p4LeptonSemi.SetPtEtaPhiM((*particlePt)[r], (*particleEta)[r], (*particlePhi)[r], (*particleMass)[r]);
						boostLorentzVectorTTbar(p4LeptonSemi, p4LeptonSemi_boosted, ttbarBoostVector);
	
						if ((*particleId)[r] > 0) //if you get l- you need the angle with the d quark from the other top -> t
						{
							angle_quark  = helAngle(p4DownBar, p4Wplus, p4T);
							angle_lepton = helAngle(p4LeptonSemi, p4Wminus, p4Tbar);
							
							boostLorentzVectorTTbar(p4DownBar, p4DownBar_boosted, ttbarBoostVector);
							angle_1 = p4DownBar_boosted.Angle(p4LeptonSemi_boosted.Vect());
						}
						else if ((*particleId)[r] < 0) //if you have l+ this means you need the angle with the tbar quark coming from the other top -> tbar
						{
							angle_quark  = helAngle(p4Down, p4Wminus, p4Tbar);
							angle_lepton = helAngle(p4LeptonSemi, p4Wplus, p4T);
							
							boostLorentzVectorTTbar(p4Down, p4Down_boosted, ttbarBoostVector);
							angle_1 = p4Down_boosted.Angle(p4LeptonSemi_boosted.Vect());
						}
						
						
						
						h_HelicityAngle[0]->Fill(TMath::Cos(angle_1));
						h_HelicityAngle[1]->Fill(TMath::Cos(angle_quark)*TMath::Cos(angle_lepton));
						
				}
				
					
			}
							
				
		}
		else if(foundLepton_2 != leptonIds+6)
		{
			    int choice = 0;
				for(int r=6; r<=7; r++)
				{
					if(fabs((*particleId)[r]) == 11 || fabs((*particleId)[r]) == 13 || fabs((*particleId)[r]) == 15)
					{
						choice = r;
						p4LeptonSemi.SetPtEtaPhiM((*particlePt)[r], (*particleEta)[r], (*particlePhi)[r], (*particleMass)[r]);
						boostLorentzVectorTTbar(p4LeptonSemi, p4LeptonSemi_boosted, ttbarBoostVector);
					}
				}
				
				
				foundQuark_3 = std::find(quarkIds,quarkIds+5, fabs((*particleId)[8])); 
				
				if(foundQuark_3 != quarkIds+5)
				{
					//this is semileptonic decay
					
					semileptonic_leptonFirst++;
					//this is only for d quark that has high A coefficient 
					for(int r=8; r<=9; r++)
					{
						if(fabs((*particleId)[r]) == 1)
						{
							if ((*particleId)[r] == 1)
							{
								p4Down.SetPtEtaPhiM((*particlePt)[r], (*particleEta)[r], (*particlePhi)[r], (*particleMass)[r]);
								boostLorentzVectorTTbar(p4Down, p4Down_boosted, ttbarBoostVector);

							}
							else 
							{
								p4DownBar.SetPtEtaPhiM((*particlePt)[r], (*particleEta)[r], (*particlePhi)[r], (*particleMass)[r]);
								boostLorentzVectorTTbar(p4DownBar, p4DownBar_boosted, ttbarBoostVector);
							}
						    
						}
						
					}
					if ((*particleId)[choice] > 0) //if you get l- (coming from tbar) you need the angle with the b quark from the other top -> b
					{
						angle_quark  = helAngle(p4DownBar, p4Wplus, p4T);
						angle_lepton = helAngle(p4LeptonSemi, p4Wminus, p4Tbar);

						boostLorentzVectorTTbar(p4LeptonSemi, p4LeptonSemi_boosted, ttbarBoostVector);
						angle_1 = p4DownBar_boosted.Angle(p4LeptonSemi_boosted.Vect());
						
					}
					else if ((*particleId)[choice] < 0) //if you have l+ (coming from top) this means you need the angle with the bbar quark coming from the other top -> bbar
					{
						angle_quark  = helAngle(p4Down, p4Wminus, p4Tbar);
						angle_lepton = helAngle(p4LeptonSemi, p4Wplus, p4T);

						boostLorentzVectorTTbar(p4LeptonSemi, p4LeptonSemi_boosted, ttbarBoostVector);
						angle_1 = p4Down_boosted.Angle(p4LeptonSemi_boosted.Vect());
					}
					semileptonic++;
					h_HelicityAngle[0]->Fill(TMath::Cos(angle_1));
					h_HelicityAngle[1]->Fill(TMath::Cos(angle_quark)*TMath::Cos(angle_lepton));
				}
					
					
				}
				
			
	
	
		
		
  }
  
  
  cout<<"Semileptonic: "<<semileptonic<<endl;
  cout<<"Semi quark first: "<<semileptonic_quarkFirst<<endl;
  cout<<"Semi lepton first: "<<semileptonic_leptonFirst<<endl;
  
  Double_t norm =1;
  for(int i=0; i<=h_HelicityAngle[1]->GetNbinsX(); i++)
  {

	 for(int j=0; j<=1; j++)
	 {		
		 float width_h_HelicityAngle = h_HelicityAngle[j]->GetXaxis()->GetBinWidth(i+1);
		 float tempValue_h_HelicityAngle = h_HelicityAngle[j]->GetBinContent(i+1);
 		 h_HelicityAngle[j]->SetBinContent(i+1, tempValue_h_HelicityAngle/width_h_HelicityAngle); 
	 }
  }
 
  
  Double_t scale_2 = norm/(h_HelicityAngle[1]->Integral());
  h_HelicityAngle[1]->Scale(scale_2);
  h_HelicityAngle[1]->GetXaxis()->SetTitle("cos(#theta^{*}_{l})cos(#theta^{*}_{d})");
  h_HelicityAngle[1]->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{d(cos(#theta^{*}_{l})cos(#theta^{*}_{d}))}");
  TCanvas *can4 = new TCanvas("can4", "can4", 900,600);
  //h_HelicityAngle[1]->SetFillColor(kRed);
  h_HelicityAngle[1]->Draw();
  
  
  
  
  Double_t scale_1 = norm/(h_HelicityAngle[0]->Integral());
  h_HelicityAngle[0]->Scale(scale_1);
  //h_HelicityAngle[0]->SetFillColor(kBlue);
  h_HelicityAngle[0]->GetXaxis()->SetTitle("cos(#theta_{l,d})");
  h_HelicityAngle[0]->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{d(cos(#theta_{l,d}))}");
  TCanvas *can1 = new TCanvas("can1", "can1", 900,600);
  h_HelicityAngle[0]->Draw();
  
  
  

  TFile *newFile = new TFile("Semileptonic_AngularDistributions.root", "RECREATE");
  h_HelicityAngle[1]->Write();
  h_HelicityAngle[0]->Write();
 
}


TVector3 getBoostVectorTTbar(TLorentzVector p4_1, TLorentzVector p4_2)
{
	//define the combined Lorentz vector of ttbar 
	TLorentzVector p4CombinedVector;
	p4CombinedVector.SetPxPyPzE(p4_1.Px()+p4_2.Px(),p4_1.Py()+p4_2.Py(), p4_1.Pz()+p4_2.Pz(), p4_1.Energy()+p4_2.Energy()); 
	//get boost from this vector
	TVector3 TTbar_boostVector = p4CombinedVector.BoostVector();
	
	return -TTbar_boostVector;
}

void boostLorentzVectorTTbar(TLorentzVector p4, TLorentzVector &p4_New, TVector3 boostVec)
{
	p4_New = p4;
	p4_New.Boost(boostVec);
	
}

//this function returns the helicity angle = angle between product of parent and the grandparent particle in the parent rest frame
Float_t helAngle(TLorentzVector particle, TLorentzVector parent, TLorentzVector grandparent)
{
	TVector3 parentBoost = -(parent.BoostVector());
	
	particle.Boost(parentBoost);
	grandparent.Boost(parentBoost);
	
	Float_t angle = grandparent.Angle(particle.Vect());
	
	return angle;	
	
}
