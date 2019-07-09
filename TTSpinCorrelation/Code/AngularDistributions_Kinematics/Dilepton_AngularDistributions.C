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
//This code will check the angular distributions of the W decay products ONLY IN DILEPTONIC FINAL STATES
//
//
//Written by: Georgios Bakas
//5/12/2018
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

using namespace std;
void boostLorentzVectorTTbar(TLorentzVector p4_1, TLorentzVector &p4_2, TVector3 boostVec);
TVector3 getBoostVectorTTbar(TLorentzVector p4_1, TLorentzVector p4_2);
Float_t helAngle(TLorentzVector particle, TLorentzVector parent, TLorentzVector grandparent);

void Dilepton_AngularDistributions(int bins = 10)
{
  //gStyle->SetOptStat(0);
  //TFile *f = TFile::Open("/afs/cern.ch/user/g/gbakas/CRAB3-jobs2016/CMSSW_8_0_6/src/UserCode/TopAnalysis/test/flatTree_10KNew.root");
  //TFile *f = TFile::Open("/afs/cern.ch/user/g/gbakas/CRAB3-jobs2016/CMSSW_8_0_6/src/UserCode/TopAnalysis/test/flatTree_10K-01-12.root");
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
  
  
  TH1F *h_dilepton = new TH1F ("Dilepton #Delta#phi Distribution", "Dilepton #Delta#phi Distribution", bins, 0,TMath::Pi()); //this is in the CM frame Delta Phi for delepton decay (lepton1, 2)
  //this is semileptonic cos(theta) where theta is the angle between the d-quark and the charged lepton in the ttbar CM frame
  TH1F *h_HelicityAngle[2], *h_HelicityAngle2[2];
  h_HelicityAngle[0] = new TH1F("Helicity angle distribtution cos(#theta_{l^{+}l^{-}})", "Helicity angle distribtution cos(#theta_{l^{+}l^{-}})", bins, -1,1); 
  h_HelicityAngle2[0] = new TH1F("Helicity angle distribtution cos(#theta_{l^{+}l^{-}})", "Helicity angle distribtution cos(#theta_{l^{+}l^{-}})", bins, -1,1); 
  h_HelicityAngle[1] = new TH1F("Helicity angle distribtution cos(#theta_{l^{+}})cos(#theta_{l^{-}}))", "Helicity angle distribtution cos(#theta_{l^{+}})cos(#theta_{l^{-}})",bins, -1,1); 
  h_HelicityAngle2[1] = new TH1F("Helicity angle distribtution cos(#theta_{l^{+}})cos(#theta_{l^{-}}))", "Helicity angle distribtution cos(#theta_{l^{+}})cos(#theta_{l^{-}})",bins, -1,1); 
  
  TLorentzVector p4T(0,0,0,0), p4Tbar(0,0,0,0),p4W(0,0,0,0), p4Wplus(0,0,0,0), p4Wminus(0,0,0,0), p4b(0,0,0,0), p4bbar(0,0,0,0);
  TLorentzVector p4T_TTbarBoost(0,0,0,0), p4Tbar_TTbarBoost(0,0,0,0), p4Wminus_boosted(0,0,0,0), p4Wplus_boosted(0,0,0,0), p4b_boosted(0,0,0,0), p4bbar_boosted(0,0,0,0);
  TLorentzVector p4Lepton2(0,0,0,0), p4Down(0,0,0,0), p4LeptonSemi(0,0,0,0), p4Lepton2_boosted(0,0,0,0), p4LeptonSemi_boosted(0,0,0,0), p4Down_boosted(0,0,0,0); 
  
    
  
  int *foundLepton, *foundLepton_2, *foundLepton_3;
  int leptonIds[] = {11,12,13,14,15,16};
  int dilepton = 0;
  long numOfEntries = t->GetEntries();
  
  Float_t angle_minus, angle_plus, angle_plus2, angle_minus2;
  //cout<<numOfEntries<<endl;
  
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
		

		foundLepton_2 = std::find(leptonIds,leptonIds+6, fabs((*particleId)[6]));
		Double_t angle;
		TVector3 tempVector;
		bool foundDownQuark;
		
		if(foundLepton_2 != leptonIds+6)
		{
			    int choice = 0;
				if(fabs((*particleId)[6]) == 11 || fabs((*particleId)[6]) == 13 || fabs((*particleId)[6]) == 15)
				{
					choice = 6;
					p4LeptonSemi.SetPtEtaPhiM((*particlePt)[6], (*particleEta)[6], (*particlePhi)[6], (*particleMass)[6]);
					boostLorentzVectorTTbar(p4LeptonSemi, p4LeptonSemi_boosted, ttbarBoostVector);
				}
				else if(fabs((*particleId)[7]) == 11 || fabs((*particleId)[7]) == 13 || fabs((*particleId)[7]) == 15)
				{
					choice = 7;
					p4LeptonSemi.SetPtEtaPhiM((*particlePt)[7], (*particleEta)[7], (*particlePhi)[7], (*particleMass)[7]);
					boostLorentzVectorTTbar(p4LeptonSemi, p4LeptonSemi_boosted, ttbarBoostVector);
				}
				
				foundLepton_3 = std::find(leptonIds,leptonIds+6, fabs((*particleId)[8]));
				if(foundLepton_3 !=leptonIds+6)
				{
						//dileptonic
						int choice = 0;
						for(int r=8; r<=9; r++)
						{
							if(fabs((*particleId)[r]) == 11 || fabs((*particleId)[r] == 13) || fabs((*particleId)[r]) == 15)
							{
								p4Lepton2.SetPtEtaPhiM((*particlePt)[r], (*particleEta)[r], (*particlePhi)[r], (*particleMass)[r]);		
								//this is delta phi of the charged leptons in the lab frame
								h_dilepton->Fill(fabs(p4LeptonSemi.Phi()-p4Lepton2.Phi()));
								choice = r;
								dilepton++;
								break;
							}
						
						}
						//now I need to fill distributions for cos(theta_hel_l+l-) and cos(theta_hel_l+)*cos(theta_hel_l-)
						if((*particleId)[choice] > 0) //this means its a lepton - so it comes from W- so from a anti-top
						{						
							angle_minus = helAngle(p4Lepton2,p4Wminus, p4Tbar); 
							angle_plus = helAngle(p4LeptonSemi, p4Wplus, p4T);
							
							//will also try to do it another way... there is explanation that helicity angle is calculated by the angle between the lepton in the top rest frame and the top momentum 
							//of the top in the ttbar rest frame
							//this is for lepton- --> lepton2
							TLorentzVector p4Lepton2_tBoost= p4Lepton2;		
							p4Lepton2_tBoost.Boost(-(p4Tbar.BoostVector()));
							angle_minus2 = p4Lepton2_tBoost.Angle(p4Tbar_TTbarBoost.Vect());
							
							TLorentzVector p4Lepton1_tBoost= p4LeptonSemi;		
							p4Lepton1_tBoost.Boost(-(p4T.BoostVector()));
							angle_plus2 = p4Lepton1_tBoost.Angle(p4T_TTbarBoost.Vect());
						}
						else
						{
							angle_minus = helAngle(p4LeptonSemi,p4Wminus, p4Tbar); 
							angle_plus = helAngle(p4Lepton2, p4Wplus, p4T);		

							TLorentzVector p4Lepton2_tBoost= p4Lepton2;		
							p4Lepton2_tBoost.Boost(-(p4T.BoostVector()));
							angle_plus2 = p4Lepton2_tBoost.Angle(p4T_TTbarBoost.Vect());
							
							TLorentzVector p4Lepton1_tBoost= p4LeptonSemi;		
							p4Lepton1_tBoost.Boost(-(p4Tbar.BoostVector()));
							angle_minus2 = p4Lepton1_tBoost.Angle(p4Tbar_TTbarBoost.Vect());							
						}
						h_HelicityAngle[0]->Fill(TMath::Cos(angle_minus-angle_plus));
						h_HelicityAngle[1]->Fill(TMath::Cos(angle_minus)*TMath::Cos(angle_plus));
						
						h_HelicityAngle2[0]->Fill(TMath::Cos(angle_minus2-angle_plus2));
						h_HelicityAngle2[1]->Fill(TMath::Cos(angle_minus2)*TMath::Cos(angle_plus2));		
						
				}
			//cout<<dilepton<<endl;	
		}	
	
	
		
		
  }
  
  

  
  h_dilepton->Scale(1./h_dilepton->Integral(), "width");
  h_dilepton->GetXaxis()->SetTitle("#Delta#phi");
  h_dilepton->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{d(#Delta#phi)}");
  TCanvas *can1 = new TCanvas("can1", "can1", 900,600);
  h_dilepton->Draw();

  
  h_HelicityAngle[0]->Scale(1./h_HelicityAngle[0]->Integral(),"width");
  h_HelicityAngle2[0]->Scale(1./h_HelicityAngle2[0]->Integral(),"width");
  h_HelicityAngle[0]->GetXaxis()->SetTitle("cos(#theta_{l^{+}l^{-}})");
  h_HelicityAngle[0]->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{d(cos(#theta_{l^{+}l^{-}}))}");
  TCanvas *can3 = new TCanvas("can3", "can3", 900,600);
  h_HelicityAngle[0]->Draw();
  h_HelicityAngle2[0]->SetLineColor(kRed);
  h_HelicityAngle2[0]->Draw("same");
  
  h_HelicityAngle[1]->Scale(1./h_HelicityAngle[1]->Integral(),"width");
  h_HelicityAngle2[1]->Scale(1./h_HelicityAngle2[1]->Integral(),"width");
  h_HelicityAngle[1]->GetXaxis()->SetTitle("cos(#theta_{l^{+}})cos(#theta_{l^{-}})");
  h_HelicityAngle[1]->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{d(cos(#theta_{l^{+}})cos(#theta_{l^{-}}))}");
  TCanvas *can4 = new TCanvas("can4", "can4", 900,600);
  h_HelicityAngle[1]->Draw();
  h_HelicityAngle2[1]->SetLineColor(kRed);
  h_HelicityAngle2[1]->Draw("same");
  
  


  TFile *newFile = new TFile("Dilepton_AngularDistributions.root", "RECREATE");
  h_dilepton->Write("Dilepton #Delta#Phi");  
  h_HelicityAngle[0]->Write();
  h_HelicityAngle2[0]->Write();
  h_HelicityAngle[1]->Write();
  h_HelicityAngle2[1]->Write();
 
 
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
