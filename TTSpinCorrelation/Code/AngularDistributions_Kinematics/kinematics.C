#include <iostream>
#include <TString.h>
#include <TH1F.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include "TMath.h"
#include "TLorentzVector.h"
#include "TLatex.h"
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

void kinematics(bool scSample = true)
{
  TString isSC;
  if(scSample) isSC = "SC";
  else isSC = "noSC";
  TFile *f;
  if(scSample) f = TFile::Open("/eos/cms/store/user/gbakas/spinCorrelation/mc-2018/TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root");
  else f = TFile::Open("/eos/cms/store/user/gbakas/spinCorrelation/mc-2018/TT_noSC_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root");
  TTree *t = (TTree*)f->Get("boosted/events");
  cout<<"opened tree"<<endl;
  
  //this is a histogram that contains the mass of every particle Id
  TH1F *h_mass[10], *h_pt[10], *h_E[10], *h_eta[10], *h_phi[10]; 
  const char* particleNames[10] = {"top", "W", "b quark", "u quark", "d quark", "c quark", "s quark", "e", "muon", "tau"};
  double_t upLimits[10] = {220, 110, 10, 0.02, 0.02, 3, 0.3, 0.001, 0.3, 3};
  double_t downLimits[10] = {150, 50, 0, 0 ,0 ,0 ,0 ,0 ,0 ,0};
  for(int i =1; i<=10; i++)
  {
	  TString tempLimUp, tempLimDown, tempName;
	  //tempLimUp = upLimits[i];
	  //tempLimDown = downLimits[i];
	  tempName = particleNames[i-1];
	  h_mass[i-1] = new TH1F(tempName +" Mass",tempName +" Mass", 50, downLimits[i-1], upLimits[i-1]);
	  h_pt[i-1] = new TH1F(tempName +" P_{T}",tempName +" P_{T}", 100, 0, 1500);
	  h_E[i-1] = new TH1F(tempName +" Energy",tempName +" Energy", 100, 0, 1500);
	  h_eta[i-1] = new TH1F(tempName +" Eta",tempName +" Eta", 50, -4, 4);
	  h_phi[i-1] = new TH1F(tempName +" Phi",tempName +" Phi", 50, -4, 4);
  }
  TH1F *h_ids = new TH1F("ids", "ids", 1000, -30,30);
  TH1F *h_dPhi_W = new TH1F("#Delta Phi W", "#Delta Phi", 50, -4,4);
  TH1F *h_dEta_W = new TH1F("#Delta Eta W", "#Delta Eta", 50, -4,4);
  TH1F *h_dPhi_L = new TH1F("#Delta Phi Leptons", "#Delta Eta Leptons", 50, -4,4);
  TH1F *h_dEta_L = new TH1F("#Delta Eta Leptons", "#Delta Eta Leptons", 50, -4,4);
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
  cout<<numOfEntries<<endl;
  int quarkIds[] = {1,2,3,4,5};
  int leptonIds[] = {11,12,13,14,15,16};
  int WProdPositions[] = {3,4,7,8}; //these are the positions in the vector float in which the products of W exist
  float nleptons=0;
  float nquarks = 0;
  int l =0;
  int *foundQuark, *foundLepton, *foundLepton_2, *foundLepton_3, *foundQuark_2, *foundQuark_3;
  int semileptonic = 0;
  int fullyHadronic = 0;
  int dilepton = 0;
  int semileptonic_quarkFirst = 0;
  int semileptonic_leptonFirst = 0;
  int numberOfTops(0);
  //define the Lorentz vectors	
  TLorentzVector p4T(0,0,0,0), p4Tbar(0,0,0,0),p4W(0,0,0,0), p4Wplus(0,0,0,0), p4Wminus(0,0,0,0), p4b(0,0,0,0), p4bbar(0,0,0,0),
				 p4Lepton1(0,0,0,0), p4Lepton2(0,0,0,0);
  
  
  //loop over all entries
  for(int k =0; k<=numOfEntries; k++)
  {
    t->GetEntry(k);
    /*
	cout<<"-----------------------------------------------------------------------"<<endl;
	cout<<(*particleId)[0]<<endl;
	cout<<(*particleId)[1]<<endl;
	cout<<(*particleId)[2]<<endl;
	cout<<(*particleId)[3]<<endl;
	cout<<(*particleId)[4]<<endl;
	cout<<(*particleId)[5]<<endl;
	cout<<(*particleId)[6]<<endl;
	cout<<(*particleId)[7]<<endl;
	cout<<(*particleId)[8]<<endl;
	cout<<(*particleId)[9]<<endl;
		*/		
		
		for(int l=0; l<=5; l++)
		{
				//verify the particle Id 
				if( (*particleId)[l] == 6 || (*particleId)[l] == -6)
				{
					//cout<<(*particleMass)[1]<<endl;
					h_mass[0] ->Fill((*particleMass)[l]);	
					p4T.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
					h_pt[0] ->Fill(p4T.P());
					h_E[0] ->Fill(p4T.E());
					h_eta[0] ->Fill(p4T.Eta());
					h_phi[0] ->Fill(p4T.Phi());
					numberOfTops++;
					//Printf("(P,eta,phi,E)=(%f,%f,%f,%f)",  p4T.P(),p4T.Eta(),p4T.Phi(),p4T.E());
				}
			
				//this is W
				if( (*particleId)[l] == 24 || (*particleId)[l] == -24)
				{
					p4W.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
					h_mass[1] ->Fill((*particleMass)[l]);
					h_pt[1] ->Fill(p4W.P());
					h_E[1] ->Fill(p4W.E());
					h_eta[1] ->Fill(p4W.Eta());
					h_phi[1] ->Fill(p4W.Phi());
					
					if((*particleId)[l] == 24 ) p4Wplus.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
					else p4Wminus.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
					h_dPhi_W->Fill(p4Wminus.Phi()-p4Wplus.Phi());
					h_dEta_W->Fill(p4Wminus.Eta()-p4Wplus.Eta());
				}
				//this is b, b bar from top
				else if (fabs((*particleId)[l]) == 5)
				{
					p4b.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
					h_mass[2] ->Fill((*particleMass)[l]);
					h_pt[2] ->Fill(p4b.P());
					h_E[2] ->Fill(p4b.E());
					h_eta[2] ->Fill(p4b.Eta());
					h_phi[2] ->Fill(p4b.Phi());
				}
			}
			
			
			//this is q, qbar from W and there are in positions 3,4,7,8
			for(int j=6; j<=9; j++)
			{
				foundQuark = std::find(quarkIds,quarkIds+5, fabs((*particleId)[j])); 
				foundLepton = std::find(leptonIds,leptonIds+6, fabs((*particleId)[j])); 
				int tempInt = j;
				//cout<<tempInt<<endl;
				//cout<<k<<endl;
				if (foundQuark !=quarkIds+5 && fabs((*particleId)[tempInt]) != 5)
				{
					 nquarks++;
					//now check exactly what quark is this 
					if(fabs((*particleId)[tempInt]) == 1) //down quark with id = 1
					{
						//std::cout<<"This is down "<<endl;
						h_mass[4]  ->Fill((*particleMass)[tempInt]);
					}
					if(fabs((*particleId)[tempInt]) == 2) //up quark with id =2
					{	
						//std::cout<<"This is up "<<endl;
						h_mass[3]  ->Fill((*particleMass)[tempInt]);
					}
					if(fabs((*particleId)[tempInt]) == 3) //s quark with id = 3
					{				
						//std::cout<<"This is strange "<<endl;
						h_mass[6]  ->Fill((*particleMass)[tempInt]);
					}
					if(fabs((*particleId)[tempInt]) == 4) //c quark with id = 4
					{
					//std::cout<<"This is charm "<<endl;
						h_mass[5]  ->Fill((*particleMass)[tempInt]);
					}
					/* if(fabs((*particleId)[j]) == 5) //b quark (rare!) with id =5
					{
						//std::cout<<"This is bottom "<<endl;
						h_mass[2]  ->Fill((*particleMass)[j]);
					}*/
				}	

				if(foundLepton != leptonIds+6)
				{	
					nleptons++;
					//now check exactly what lepton is this 		
					if(fabs((*particleId)[tempInt]) == 11) //electron +/- with abs_id = 11
					{
						h_mass[7]  ->Fill((*particleMass)[tempInt]);
					}
					if(fabs((*particleId)[tempInt]) == 13)  //muon +/- with abs_id = 13
					{
						h_mass[8]  ->Fill((*particleMass)[tempInt]);
					}
					if(fabs((*particleId)[tempInt]) == 15)  //tau +/- with abs_id = 15
					{
						h_mass[9]  ->Fill((*particleMass)[tempInt]);
					}
	
					
				}
				
			}
			
			//in this part I check how the W decays, so I compute how many fully hadronic, semileptonic, dileptonic
			//we know that we have decay products from [6]-[9] and they are by 2
			//so [6][7] are W1 products and [8][9] are W2 products
			//check with respect to [3]
			foundQuark_2 = std::find(quarkIds,quarkIds+5, fabs((*particleId)[6])); 
			foundLepton_2 = std::find(leptonIds,leptonIds+6, fabs((*particleId)[6]));
			if(foundQuark_2 != quarkIds+5) //at least semilepton
			{
				//cout<<"-----------------------------------------------------------------"<<endl;
				//cout<<"At 6 we have: "<<(*particleId)[6]<<endl;
				//cout<<"At 7 we have: "<<(*particleId)[7]<<endl;
				
				foundQuark_3 = std::find(quarkIds,quarkIds+5, fabs((*particleId)[8])); 
				foundLepton_3 = std::find(leptonIds,leptonIds+6, fabs((*particleId)[8]));
				
				if(foundQuark_3 != quarkIds+5)
				{
					//cout<<"hadronic"<<endl;
					//cout<<"At 8 we have: "<<(*particleId)[8]<<endl;
					//cout<<"At 9 we have: "<<(*particleId)[9]<<endl;
					fullyHadronic++; //this is fully hadronic decay
				}
				else if(foundLepton_3 !=leptonIds+6)
				{
					//cout<<"semileptonic: "<<endl;
					//cout<<"At 8 we have: "<<(*particleId)[8]<<endl;
					//cout<<"At 9 we have: "<<(*particleId)[9]<<endl;
					semileptonic++; //this is semileptonic decay
					semileptonic_quarkFirst++;
				}
				
				
			}
			else if(foundLepton_2 != leptonIds+6)
			{
				p4Lepton1.SetPtEtaPhiM((*particlePt)[6], (*particleEta)[6], (*particlePhi)[6], (*particleMass)[6]);
				foundQuark_3 = std::find(quarkIds,quarkIds+5, fabs((*particleId)[8])); 
				foundLepton_3 = std::find(leptonIds,leptonIds+6, fabs((*particleId)[8]));
				
				if(foundQuark_3 != quarkIds+5)
				{
					semileptonic++; //this is semileptonic decay
					semileptonic_leptonFirst++;
				}
				else if(foundLepton_3 !=leptonIds+6)
				{
					p4Lepton2.SetPtEtaPhiM((*particlePt)[8], (*particleEta)[8], (*particlePhi)[8], (*particleMass)[8]);
					dilepton++; //this is dileptonic decay
					h_dPhi_L->Fill(p4Wminus.Phi()-p4Wplus.Phi());
					h_dEta_L->Fill(p4Wminus.Eta()-p4Wplus.Eta());
				}
				
			}
	
	
	  
  }
  
  cout<<"Number of tops: "<<numberOfTops<<endl;
  cout<<"Fully Hadronic: "<<fullyHadronic<<endl;
  cout<<"Semileptonic: "<<semileptonic<<endl;
  cout<<"Dileptonic: "<<dilepton<<endl;
  cout<<"Semi quark first: "<<semileptonic_quarkFirst<<endl;
  cout<<"Semi lepton first: "<<semileptonic_leptonFirst<<endl;
  
  TFile *massesFile =  TFile::Open("Kinematics.root", "UPDATE");
  massesFile->cd();
  
  
  
  for(int i =1; i<=3; i++)
  {
	
	h_mass[i-1]->GetXaxis()->SetTitle("mass [GeV]");
	h_mass[i-1]->GetYaxis()->SetTitle("Events");
	h_mass[i-1]->Write(TString::Format("h_mass_%d_%s",(i-1), isSC.Data()));
	h_pt[i-1]->Write(TString::Format("h_pt_%d_%s",(i-1), isSC.Data()));
	h_E[i-1]->Write(TString::Format("h_E_%d_%s",(i-1), isSC.Data()));
	h_eta[i-1]->Write(TString::Format("h_eta_%d_%s",(i-1), isSC.Data()));
	h_phi[i-1]->Write(TString::Format("h_phi_%d_%s",(i-1), isSC.Data()));
  }
  h_dEta_W->Write(TString::Format("h_dEta_W_%s",isSC.Data()));
  h_dPhi_W->Write(TString::Format("h_dPhi_W_%s",isSC.Data()));
  h_dEta_L->Write(TString::Format("h_dEta_L_%s",isSC.Data()));
  h_dPhi_L->Write(TString::Format("h_dPhi_L_%s",isSC.Data()));
  cout<<h_mass[0]->GetEntries()<<endl;
  cout<<h_mass[1]->GetEntries()<<endl;
  cout<<h_mass[2]->GetEntries()<<endl;
  
  
  cout<<"Number of Quarks: "<<nquarks<<endl;
  cout<<"Number of Leptons: "<<nleptons<<endl;
  cout<<"Fraction: "<<nquarks/nleptons <<endl;
  
  //p4T.Print();
  //p4Tbar.Print();
  
  massesFile->Close();
  f->Close();
  
}
