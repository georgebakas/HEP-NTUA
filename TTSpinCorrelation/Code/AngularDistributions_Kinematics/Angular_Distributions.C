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

void Angular_Distributions(int bins = 10, bool scSample = true)
{
  //gStyle->SetOptStat(0);
  TString isSC;
  if(scSample) isSC = "SC";
  else isSC = "noSC";
  TFile *f;
  if(scSample) f = TFile::Open("/eos/cms/store/user/gbakas/spinCorrelation/mc-2018/TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root");
  else f = TFile::Open("/eos/cms/store/user/gbakas/spinCorrelation/mc-2018/TT_noSC_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root");
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
  TH1F *h_dileptonMassCut = new TH1F ("Dilepton #Delta#phi Distribution with mass cut", "Dilepton #Delta#phi Distribution with m_{t#bar{t}} >450", bins, 0,TMath::Pi()); //this is in the CM frame Delta Phi for delepton decay (lepton1, 2)
  //this is semileptonic cos(theta) where theta is the angle between the d-quark and the charged lepton in the ttbar CM frame
  TH1F *h_semilepton = new TH1F ("Semileptonic cos#theta Distribution", "Semileptonic cos#theta Distribution", 20, -1,1); 
  TH1F *h_semileptonPhi[2];
  h_semileptonPhi[0] = new TH1F ("Semileptonic #Delta#phi Distribution", "Semileptonic #Delta#phi Distribution", bins, -3,3); 
  h_semileptonPhi[1] = new TH1F ("Semileptonic #Sigma#phi Distribution", "Semileptonic #Sigma#phi Distribution", bins, -3,3); 
  
  TH1F *h_HelicityAngle[2];
  h_HelicityAngle[0] = new TH1F("Helicity angle distribtution cos(#theta_{l^{+}l^{-}})", "Helicity angle distribtution cos(#theta_{l^{+}l^{-}})", bins, -1,1); 
  h_HelicityAngle[1] = new TH1F("Helicity angle distribtution cos(#theta_{l^{+}})cos(#theta_{l^{-}}))", "Helicity angle distribtution cos(#theta_{l^{+}})cos(#theta_{l^{-}})",bins, -1,1); 
  
  TH1F *h_coscos   = new TH1F("cos(lepton1)*cos(lepton2)", "cos(lepton1)*cos(lepton2)", 20, -1,1);
  TH1F *h_deltaEta = new TH1F("Delta Eta of dilepton", "Delta Eta of dilepton", 20, 0,3);
  TH1F *h_cosPhiDilepton = new TH1F("cos(Phi_{lepton1,2}) in dilepton", "cos(Phi_{lepton1,2}) in dilepton", 20, -1,1);
  
  h_coscos->GetXaxis()->SetTitle("coscos");
  h_deltaEta->GetXaxis()->SetTitle("#Delta#eta");
  h_cosPhiDilepton->GetXaxis()->SetTitle("cos(#phi_{lep 1,2})");
  
  TLorentzVector p4T(0,0,0,0), p4Tbar(0,0,0,0),p4W(0,0,0,0), p4Wplus(0,0,0,0), p4Wminus(0,0,0,0), p4b(0,0,0,0), p4bbar(0,0,0,0);
  TLorentzVector p4T_TTbarBoost(0,0,0,0), p4Tbar_TTbarBoost(0,0,0,0), p4Wminus_boosted(0,0,0,0), p4Wplus_boosted(0,0,0,0), p4b_boosted(0,0,0,0), p4bbar_boosted(0,0,0,0);
  TLorentzVector p4Lepton2(0,0,0,0), p4Down(0,0,0,0), p4LeptonSemi(0,0,0,0), p4Lepton2_boosted(0,0,0,0), p4LeptonSemi_boosted(0,0,0,0), p4Down_boosted(0,0,0,0); 
  
   
  /*
  used for semileptonic delta phi and sigma phi distributions
  for these angular distribtutions
  1. Boost in ttbar rest frame
  2. Boost in top products in top rest frame and antitop producs in the respective rest frame
  3. Get phi and calculate deltaPhi and sigmaPhi of the decay products
  Comment: Here I will try with the b , bbar quarks and maybe with Down , Downbar quarks
  */
  TLorentzVector p4Lepton_semiBoosted(0,0,0,0), p4Down_semiBoosted(0,0,0,0), p4DownBar_semiBoosted(0,0,0,0), p4b_semiBoosted(0,0,0,0), p4bbar_semiBoosted; 
  
  
  int *foundQuark, *foundLepton, *foundLepton_2, *foundLepton_3, *foundQuark_2, *foundQuark_3;
  int quarkIds[] = {1,2,3,4,5};
  int leptonIds[] = {11,12,13,14,15,16};
  int semileptonic = 0;
  int fullyHadronic = 0;
  int dilepton = 0;
  int dileptonMassCut = 0;
  long numOfEntries = t->GetEntries();
  
  Float_t angle_minus, angle_plus;
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
		//cout<<"entry: "<<k<<endl;
		//p4T.Print();
		//p4Tbar.Print();
		boostLorentzVectorTTbar(p4T, p4T_TTbarBoost,ttbarBoostVector );
		boostLorentzVectorTTbar(p4Tbar, p4Tbar_TTbarBoost,ttbarBoostVector );
		//cout<<"after boost"<<endl;
		//p4T_TTbarBoost.Print();
		//p4Tbar_TTbarBoost.Print();
		//TLorentzVector p4CombinedVector(0,0,0,0),p4CombinedBoostedVector(0,0,0,0),p4CombinedBoostedVector_2(0,0,0,0) ;
		
		
		/*
		p4CombinedVector.SetPxPyPzE(p4T.Px()+p4Tbar.Px(),p4T.Py()+p4Tbar.Py(), p4T.Pz()+p4Tbar.Pz(), p4T.Energy()+p4Tbar.Energy()); 
		boostLorentzVectorTTbar(p4CombinedVector, p4CombinedBoostedVector, ttbarBoostVector);
		p4CombinedBoostedVector_2.SetPxPyPzE(p4T_TTbarBoost.Px()+p4Tbar_TTbarBoost.Px(),p4T_TTbarBoost.Py()+p4Tbar_TTbarBoost.Py(), p4T_TTbarBoost.Pz()+p4Tbar_TTbarBoost.Pz(), p4T_TTbarBoost.Energy()+p4Tbar_TTbarBoost.Energy()); 
		
		std::cout<<"before boost "<<endl;
		p4CombinedVector.Print();
		std::cout<<"after boost "<< endl;
		p4CombinedBoostedVector.Print();
		p4CombinedBoostedVector_2.Print();
		*/
		
		boostLorentzVectorTTbar(p4Wplus,p4Wplus_boosted, ttbarBoostVector);
		boostLorentzVectorTTbar(p4Wminus,p4Wminus_boosted, ttbarBoostVector);
		boostLorentzVectorTTbar(p4b,p4b_boosted, ttbarBoostVector);
		boostLorentzVectorTTbar(p4bbar,p4bbar_boosted, ttbarBoostVector);
		
	//boost W decay products	
	
	//this is q, qbar from W and there are in positions 3,4,7,8
//if(p4T_TTbarBoost.M() + p4Tbar_TTbarBoost.M() < 400)
	//{
		foundQuark_2 = std::find(quarkIds,quarkIds+5, fabs((*particleId)[6])); 
		foundLepton_2 = std::find(leptonIds,leptonIds+6, fabs((*particleId)[6]));
		Double_t angle;
		TVector3 tempVector;
		bool foundDownQuark;
		if(foundQuark_2 != quarkIds+5) //at least semilepton
		{
			if(fabs((*particleId)[6] == 1))
			{
				foundDownQuark = 1;
				p4Down.SetPtEtaPhiM((*particlePt)[6], (*particleEta)[6], (*particlePhi)[6], (*particleMass)[6]);
				boostLorentzVectorTTbar(p4Down,p4Down_boosted, ttbarBoostVector);
				
			}
			else if(fabs((*particleId)[7] == 1))
			{
				foundDownQuark = 1;
				p4Down.SetPtEtaPhiM((*particlePt)[7], (*particleEta)[7], (*particlePhi)[7], (*particleMass)[7]);
				boostLorentzVectorTTbar(p4Down,p4Down_boosted, ttbarBoostVector);		
			}
			//cout<<"-----------------------------------------------------------------"<<endl;
			//cout<<"At 6 we have: "<<(*particleId)[6]<<endl;
			//cout<<"At 7 we have: "<<(*particleId)[7]<<endl;
			
			foundQuark_3 = std::find(quarkIds,quarkIds+5, fabs((*particleId)[8])); 
			foundLepton_3 = std::find(leptonIds,leptonIds+6, fabs((*particleId)[8]));
			
			if(foundQuark_3 != quarkIds+5)
			{
				//fully hadronic
				fullyHadronic++;
			}
			else if(foundLepton_3 !=leptonIds+6)
			{
				//semileptonic
				semileptonic++;
				if(foundDownQuark == 1)
				{
					if(fabs((*particleId)[8]) == 11 || fabs((*particleId)[8]) == 13 || fabs((*particleId)[8]) == 15)
					{
						p4LeptonSemi.SetPtEtaPhiM((*particlePt)[8], (*particleEta)[8], (*particlePhi)[8], (*particleMass)[8]);
						boostLorentzVectorTTbar(p4LeptonSemi, p4LeptonSemi_boosted, ttbarBoostVector);
						//angle = p4LeptonSemi_boosted.Angle(p4Down_boosted.Vect());
						//h_semilepton->Fill(TMath::Cos(angle));
						
						if ((*particleId)[8] > 0) //if you get l- you need the angle with the b quark from the other top -> b
						{
							boostLorentzVectorTTbar(p4LeptonSemi_boosted, p4Lepton_semiBoosted, -(p4Tbar_TTbarBoost.BoostVector()));
							boostLorentzVectorTTbar(p4b_boosted, p4b_semiBoosted, -(p4T_TTbarBoost.BoostVector()));
							h_semileptonPhi[0]->Fill(p4Lepton_semiBoosted.Phi() - p4b_boosted.Phi()); // φ - φBar --> top decay prod - tbar decay prod
							h_semileptonPhi[1]->Fill(p4Lepton_semiBoosted.Phi() + p4b_boosted.Phi());
							angle = p4b_boosted.Angle(p4Lepton_semiBoosted.Vect());
							h_semilepton->Fill(TMath::Cos(angle));
						}
						else if ((*particleId)[8] < 0) //if you have l+ this means you need the angle with the bbar quark coming from the other top -> bbar
						{
							boostLorentzVectorTTbar(p4LeptonSemi_boosted, p4Lepton_semiBoosted, -(p4T_TTbarBoost.BoostVector()));
							boostLorentzVectorTTbar(p4bbar_boosted, p4bbar_semiBoosted, -(p4Tbar_TTbarBoost.BoostVector()));
							h_semileptonPhi[0]->Fill(p4Lepton_semiBoosted.Phi() - p4bbar_semiBoosted.Phi()); // φ - φBar --> top decay prod - tbar decay prod
							h_semileptonPhi[1]->Fill(p4Lepton_semiBoosted.Phi() + p4bbar_semiBoosted.Phi()); // φ - φBar --> top decay prod + tbar decay prod
							angle = p4bbar_semiBoosted.Angle(p4Lepton_semiBoosted.Vect());
							h_semilepton->Fill(TMath::Cos(angle));
						}
						
					}
					else if(fabs((*particleId)[9]) == 11 || fabs((*particleId)[9]) == 13 || fabs((*particleId)[9]) == 15)
					{
						p4LeptonSemi.SetPtEtaPhiM((*particlePt)[9], (*particleEta)[9], (*particlePhi)[9], (*particleMass)[9]);
						boostLorentzVectorTTbar(p4LeptonSemi, p4LeptonSemi_boosted, ttbarBoostVector);
								
						
						if ((*particleId)[9] > 0) //if you get l- (coming from tbar) you need the angle with the b quark from the other top -> b
						{
							boostLorentzVectorTTbar(p4LeptonSemi_boosted, p4Lepton_semiBoosted, -(p4Tbar_TTbarBoost.BoostVector()));
							boostLorentzVectorTTbar(p4b_boosted, p4b_semiBoosted, -(p4T_TTbarBoost.BoostVector()));
							h_semileptonPhi[0]->Fill(p4Lepton_semiBoosted.Phi() - p4b_boosted.Phi()); // φ - φBar --> top decay prod - tbar decay prod
							h_semileptonPhi[1]->Fill(p4Lepton_semiBoosted.Phi() + p4b_boosted.Phi());
							angle = p4b_boosted.Angle(p4Lepton_semiBoosted.Vect());
							h_semilepton->Fill(TMath::Cos(angle));
						}
						else if ((*particleId)[9] < 0)//if you have l+ (coming from top) this means you need the angle with the bbar quark coming from the other top -> bbar
						{
							boostLorentzVectorTTbar(p4LeptonSemi_boosted, p4Lepton_semiBoosted, -(p4T_TTbarBoost.BoostVector()));
							boostLorentzVectorTTbar(p4bbar_boosted, p4bbar_semiBoosted, -(p4Tbar_TTbarBoost.BoostVector()));
							h_semileptonPhi[0]->Fill(p4Lepton_semiBoosted.Phi() - p4bbar_semiBoosted.Phi()); // φ - φBar --> top decay prod - tbar decay prod
							h_semileptonPhi[1]->Fill(p4Lepton_semiBoosted.Phi() + p4bbar_semiBoosted.Phi()); // φ - φBar --> top decay prod + tbar decay prod
							angle = p4bbar_semiBoosted.Angle(p4Lepton_semiBoosted.Vect());
							h_semilepton->Fill(TMath::Cos(angle));
						}
					}
				}
				
			}
					
				
		}
		else if(foundLepton_2 != leptonIds+6)
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
				
				
				foundQuark_3 = std::find(quarkIds,quarkIds+5, fabs((*particleId)[8])); 
				foundLepton_3 = std::find(leptonIds,leptonIds+6, fabs((*particleId)[8]));
				
				if(foundQuark_3 != quarkIds+5)
				{
					//this is semileptonic decay
					semileptonic++;
					//this is only for d quark that has high A coefficient 
					if(fabs((*particleId)[8]) == 1)
					{
						p4Down.SetPtEtaPhiM((*particlePt)[8], (*particleEta)[8], (*particlePhi)[8], (*particleMass)[8]);
						boostLorentzVectorTTbar(p4Down,p4Down_boosted, ttbarBoostVector);
						
		

						if ((*particleId)[choice] > 0) //if you get l- (coming from tbar) you need the angle with the b quark from the other top -> b
						{
							boostLorentzVectorTTbar(p4LeptonSemi_boosted, p4Lepton_semiBoosted, -(p4Tbar_TTbarBoost.BoostVector()));
							boostLorentzVectorTTbar(p4b_boosted, p4b_semiBoosted, -(p4T_TTbarBoost.BoostVector()));
							h_semileptonPhi[0]->Fill(p4Lepton_semiBoosted.Phi() - p4b_boosted.Phi()); // φ - φBar --> top decay prod - tbar decay prod
							h_semileptonPhi[1]->Fill(p4Lepton_semiBoosted.Phi() + p4b_boosted.Phi());
							angle = p4b_boosted.Angle(p4Lepton_semiBoosted.Vect());
							h_semilepton->Fill(TMath::Cos(angle));
						}
						else if ((*particleId)[choice] < 0) //if you have l+ (coming from top) this means you need the angle with the bbar quark coming from the other top -> bbar
						{
							boostLorentzVectorTTbar(p4LeptonSemi_boosted, p4Lepton_semiBoosted, -(p4T_TTbarBoost.BoostVector()));
							boostLorentzVectorTTbar(p4bbar_boosted, p4bbar_semiBoosted, -(p4Tbar_TTbarBoost.BoostVector()));
							h_semileptonPhi[0]->Fill(p4Lepton_semiBoosted.Phi() - p4bbar_semiBoosted.Phi()); // φ - φBar --> top decay prod - tbar decay prod
							h_semileptonPhi[1]->Fill(p4Lepton_semiBoosted.Phi() + p4bbar_semiBoosted.Phi()); // φ - φBar --> top decay prod + tbar decay prod
							angle = p4bbar_semiBoosted.Angle(p4Lepton_semiBoosted.Vect());
							h_semilepton->Fill(TMath::Cos(angle));
						}
					}
					else if(fabs((*particleId)[9]) == 1)
					{
						p4Down.SetPtEtaPhiM((*particlePt)[9], (*particleEta)[9], (*particlePhi)[9], (*particleMass)[9]);
						boostLorentzVectorTTbar(p4Down,p4Down_boosted, ttbarBoostVector);
						
						
						if ((*particleId)[choice] > 0) //if you get l- (coming from tbar) you need the angle with the b quark from the other top -> b
						{
							boostLorentzVectorTTbar(p4LeptonSemi_boosted, p4Lepton_semiBoosted, -(p4Tbar_TTbarBoost.BoostVector()));
							boostLorentzVectorTTbar(p4b_boosted, p4b_semiBoosted, -(p4T_TTbarBoost.BoostVector()));
							h_semileptonPhi[0]->Fill(p4Lepton_semiBoosted.Phi() - p4b_boosted.Phi()); // φ - φBar --> top decay prod - tbar decay prod
							h_semileptonPhi[1]->Fill(p4Lepton_semiBoosted.Phi() + p4b_boosted.Phi());
							angle = p4b_boosted.Angle(p4Lepton_semiBoosted.Vect());
							h_semilepton->Fill(TMath::Cos(angle));
						}
						else if ((*particleId)[choice] < 0) //if you have l+ (coming from top) this means you need the angle with the bbar quark coming from the other top -> bbar
						{
							boostLorentzVectorTTbar(p4LeptonSemi_boosted, p4Lepton_semiBoosted, -(p4T_TTbarBoost.BoostVector()));
							boostLorentzVectorTTbar(p4bbar_boosted, p4bbar_semiBoosted, -(p4Tbar_TTbarBoost.BoostVector()));
							h_semileptonPhi[0]->Fill(p4Lepton_semiBoosted.Phi() - p4bbar_semiBoosted.Phi()); // φ - φBar --> top decay prod - tbar decay prod
							h_semileptonPhi[1]->Fill(p4Lepton_semiBoosted.Phi() + p4bbar_semiBoosted.Phi()); // φ - φBar --> top decay prod + tbar decay prod
							angle = p4bbar_semiBoosted.Angle(p4Lepton_semiBoosted.Vect());
							h_semilepton->Fill(TMath::Cos(angle));
						}
						

					}
				}
				else if(foundLepton_3 !=leptonIds+6)
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
								h_coscos->Fill(TMath::Cos(p4LeptonSemi.Phi())*TMath::Cos(p4Lepton2.Phi()));
								h_deltaEta->Fill(fabs(p4LeptonSemi.Eta()-p4Lepton2.Eta()));
								h_cosPhiDilepton->Fill(TMath::Cos(p4LeptonSemi.Phi()));
								h_cosPhiDilepton->Fill(TMath::Cos(p4Lepton2.Phi()));
								
								choice = r;
								dilepton++;
								if(p4T.M() + p4Tbar.M() > 450 && p4T.M() )
								{
									dileptonMassCut++;
									h_dileptonMassCut->Fill(fabs(p4LeptonSemi.Phi()-p4Lepton2.Phi()));
								}
								break;
							}
						
						}
						//now I need to fill distributions for cos(theta_hel_l+l-) and cos(theta_hel_l+)*cos(theta_hel_l-)
						if((*particleId)[choice] > 0) //this means its a lepton - so it comes from W-
						{						
							angle_minus = helAngle(p4Lepton2,p4Wminus, p4Tbar); 
							angle_plus = helAngle(p4LeptonSemi, p4Wplus, p4T);
						}
						else
						{
							angle_minus = helAngle(p4LeptonSemi,p4Wminus, p4Tbar); 
							angle_plus = helAngle(p4Lepton2, p4Wplus, p4T);							
						}
						h_HelicityAngle[0]->Fill(TMath::Cos(angle_minus-angle_plus));
						h_HelicityAngle[1]->Fill(TMath::Cos(angle_minus)*TMath::Cos(angle_plus));

							
						
				}
//		cout<<dilepton<<endl;	
		}	
	
	//}
		
		
  }
  
  
  cout<<"Fully Hadronic: "<<fullyHadronic<<endl;
  cout<<"Semileptonic: "<<semileptonic<<endl;
  cout<<"Dileptonic: "<<dilepton<<endl;
  cout<<"Dileptonic with mass cut: "<<dileptonMassCut<<endl;
  
  
  Double_t norm =1;
  Double_t scale = norm /(h_dilepton->Integral());
  for(int i=0; i<=h_dilepton->GetNbinsX(); i++)
  {
	 float width1 = h_dilepton->GetXaxis()->GetBinWidth(i+1);
	 float widMassCut = h_dileptonMassCut->GetXaxis()->GetBinWidth(i+1);
	 float width2 = h_semilepton->GetXaxis()->GetBinWidth(i+1);
	 float tempValue = h_dilepton->GetBinContent(i+1);
	 float tempVal_massCut = h_dileptonMassCut->GetBinContent(i+1);
	 h_dilepton->SetBinContent(i+1,tempValue/width1); 
	 h_dileptonMassCut->SetBinContent(i+1, tempVal_massCut/widMassCut);
	 float tempValue2 = h_semilepton->GetBinContent(i+1);
	 h_semilepton->SetBinContent(i+1, tempValue2/width2);

		


	 for(int j=0; j<=1; j++)
	 {
		 float wid = h_HelicityAngle[j]->GetXaxis()->GetBinWidth(i+1);
		 float temp = h_HelicityAngle[j]->GetBinContent(i+1);
		 h_HelicityAngle[j]->SetBinContent(i+1, temp/wid);
		 
		 float width_semileptonPhi = h_semileptonPhi[j]->GetXaxis()->GetBinWidth(i+1);
		 float tempValue_semileptonPhi = h_semileptonPhi[j]->GetBinContent(i+1);
 		 h_semileptonPhi[j]->SetBinContent(i+1, tempValue_semileptonPhi/width_semileptonPhi);
		 
	 }
  }
  h_dilepton->Scale(scale);
  h_dilepton->GetXaxis()->SetTitle("#Delta#phi");
  h_dilepton->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{d(#Delta#phi)}");
  TCanvas *can1 = new TCanvas("can1", "can1", 900,600);
  h_dilepton->Draw();
  //h_dileptonMassCut->SetLineColor(kRed);
  //h_dileptonMassCut->Scale(scale);
  //h_dileptonMassCut->Draw("same");
  
  Double_t scale_semilepton = norm/(h_semilepton->Integral());
  h_semilepton->Scale(scale_semilepton);
  h_semilepton->GetXaxis()->SetTitle("cos#theta");
  h_semilepton->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{d(cos#theta)}");
  TCanvas *can2 = new TCanvas("can2", "can2", 900,600);
  h_semilepton->Draw();
  
  Double_t scale_1 = norm/(h_HelicityAngle[0]->Integral());
  h_HelicityAngle[0]->Scale(scale_1);
  h_HelicityAngle[0]->GetXaxis()->SetTitle("cos(#theta_{l^{+}l^{-}})");
  h_HelicityAngle[0]->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{d(cos(#theta_{l^{+}l^{-}}))}");
  TCanvas *can3 = new TCanvas("can3", "can3", 900,600);
  h_HelicityAngle[0]->Draw();
  
  Double_t scale_2 = norm/(h_HelicityAngle[1]->Integral());
  h_HelicityAngle[1]->Scale(scale_1);
  h_HelicityAngle[1]->GetXaxis()->SetTitle("cos(#theta_{l^{+}})cos(#theta_{l^{-}})");
  h_HelicityAngle[1]->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{d(cos(#theta_{l^{+}})cos(#theta_{l^{-}}))}");
  TCanvas *can4 = new TCanvas("can4", "can4", 900,600);
  h_HelicityAngle[1]->Draw();
	
  Double_t scale_semileptonPhi_1 = norm/(h_semileptonPhi[0]->Integral());
  h_semileptonPhi[0]->Scale(scale_semileptonPhi_1);
  h_semileptonPhi[0]->GetXaxis()->SetTitle("#Delta#phi");
  h_semileptonPhi[0]->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{d(#Delta#phi)}");
  TCanvas *can5 = new TCanvas("can5", "can5", 900,600);
  h_semileptonPhi[0]->Draw();
  
  Double_t scale_semileptonPhi_2 = norm/(h_semileptonPhi[1]->Integral());
  h_semileptonPhi[1]->Scale(scale_semileptonPhi_2);
  h_semileptonPhi[1]->GetXaxis()->SetTitle("#Sigma#phi");
  h_semileptonPhi[1]->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{d(#Sigma#phi)}");
  TCanvas *can6 = new TCanvas("can6", "can6", 900,600);
  h_semileptonPhi[1]->Draw();
  
  
  TCanvas *can7 = new TCanvas("can7", "can7", 900,600);
  h_coscos->Draw();

  TCanvas *can8 = new TCanvas("can8", "can8", 900,600);
  h_deltaEta->Draw();
  
  TCanvas *can9 = new TCanvas("can9", "can9", 900,600);
  h_cosPhiDilepton->Draw();
  

  TFile *newFile = new TFile("Angular_Distributions.root", "UPDATE");
  h_dilepton->Write(TString::Format("Dilepton_deltaPhi_%s",isSC.Data()));
  h_dileptonMassCut->Write(TString::Format("h_dileptonMassCut_%s",isSC.Data()));
  h_semilepton->Write(TString::Format("cosTheta_%s", isSC.Data()));
  h_HelicityAngle[0]->Write(TString::Format("h_HelicityAngle_0_%s",isSC.Data()));
  h_HelicityAngle[1]->Write(TString::Format("h_HelicityAngle_1_%s",isSC.Data()));
  h_semileptonPhi[0]->Write(TString::Format("h_semileptonPhi_0_%s",isSC.Data()));
  h_semileptonPhi[1]->Write(TString::Format("h_semileptonPhi_1_%s",isSC.Data()));

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
