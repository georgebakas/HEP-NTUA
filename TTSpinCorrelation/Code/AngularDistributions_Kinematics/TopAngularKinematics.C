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
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

using namespace std;

TVector3 getBoostVector(TLorentzVector p4_1, TLorentzVector p4_2, TLorentzVector &p4CombinedVector);


void TopAngularKinematics()
{
  //TFile *f = TFile::Open("/afs/cern.ch/user/g/gbakas/CRAB3-jobs2016/CMSSW_8_0_6/src/UserCode/TopAnalysis/prod/ttbar/flatTree_100K-13-2-2019.root");

  //TFile *f = TFile::Open("/afs/cern.ch/work/g/gbakas/private/TTSpinCorrelation/boostedTTBar_SpCorr/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_BoostedPart.root");
  TFile *f = TFile::Open("/eos/cms/store/user/gbakas/spinCorrelation/mc-2018/TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root");
  //TFile *f = TFile::Open("/eos/cms/store/user/gbakas/spinCorrelation/mc-2018/TT_noSC_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root");
  TTree *t = (TTree*)f->Get("boosted/events");
  cout<<"opened tree"<<endl;
  int size = 50;
  //this is a histogram that contains the mass of every particle Id
  const char* particleNames = "top";
  TString tempName = particleNames;
  TH1F *h_mass = new TH1F(tempName +" Mass",tempName +" Mass", 20, 165, 185);
  TH1F *h_pt = new TH1F(tempName +" P_{T}",tempName +" P_{T}", 100, 0, 1500);
  TH1F *h_E = new TH1F(tempName +" Energy",tempName +" Energy", 100, 0, 1500);
  TH1F *h_eta = new TH1F(tempName +" Eta",tempName +" Eta", 50, -4, 4);
  TH1F *h_phi = new TH1F(tempName +" Phi",tempName +" Phi", 50, -4, 4);
  
  TH1F *h_coscos = new TH1F("cos#phi cos#phi", "cos#phi cos#phi", size, -1, 1);
  TH1F *h_cos = new TH1F("cos#phi", "cos#phi ", size, -1, 1);
  TH1F *h_deltaEta = new TH1F("#Delta#eta ","#Delta#eta ", size,0, 1);
  TH1F *h_deltaPhi = new TH1F("#Delta#phi ","#Delta#phi ", size, 0, 1);
  TH1F *h_etaTop = new TH1F("#eta top", "#eta top", size,-3,3);
  TH1F *h_phiTop = new TH1F("#phi top", "#phi top", size,-3,3);
  

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


  int leptonIds[] = {11,13,15};
  int *foundLepton_1, *foundLepton_2;
  int flag(0);
  int numberOfTops(0), dilepton(0);
  //define the Lorentz vectors	
  TLorentzVector p4T(0,0,0,0), p4Tbar(0,0,0,0), p4leptonPlus(0,0,0,0), p4leptonMinus(0,0,0,0); //all of these are in Lab frame
  TLorentzVector p4T_ttbarZMF(0,0,0,0), p4Tbar_ttbarZMF(0,0,0,0);
  TLorentzVector leptonPlus_topZMF(0,0,0,0), leptonMinus_tbarZMF(0,0,0,0); //leptons in top and antitop Rest frame
  TLorentzVector p4leptonPlus_ttbarZMF(0,0,0,0), p4leptonMinus_ttbarZMF(0,0,0,0); //leptons in ttbar rest frame --> Zero Momentum Frame
  TLorentzVector p4TTbar(0,0,0,0);
  int decade(0);
  float _cosphi, _coscos;
  bool isDilepton;
  
  //loop over all entries
  for(int k =0; k<=numOfEntries; k++)
  {
	double progress = 10.0*k/(1.0*numOfEntries);
    int i = TMath::FloorNint(progress); 
    if (i > decade) 
      cout<<10*i<<" %"<<endl;
    decade = i;
    t->GetEntry(k);
    flag =0;
	    for(int l=0; l<=9; l++)
		{
				//verify the particle Id 
				if( (*particleId)[l] == 6)
				{
					//cout<<(*particleMass)[1]<<endl;
					h_mass ->Fill((*particleMass)[l]);	
					p4T.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
					h_pt ->Fill(p4T.Pt());
					h_E ->Fill(p4T.E());
					h_eta ->Fill(p4T.Eta());
					h_phi ->Fill(p4T.Phi());
					
					//Printf("(P,eta,phi,E)=(%f,%f,%f,%f)",  p4T.Pt(),p4T.Eta(),p4T.Phi(),p4T.E());
				}
				if((*particleId)[l] == -6)
				{
					p4Tbar.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);

					h_eta->Fill(p4Tbar.Eta());
					h_mass ->Fill((*particleMass)[l]);	
					h_pt ->Fill(p4Tbar.Pt());
					h_phi->Fill(p4Tbar.Phi());

				}	
				
		

			if((*particleId)[l] == 11 || (*particleId)[l]==13 ) //this is electron or muon-
			{
				flag++;
				//cout<<"Minus Sign lepton with pdgId: "<<(*particleId)[l]<<endl;	
				p4leptonMinus.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);	
			}
			else if((*particleId)[l] ==-11 || (*particleId)[l]==-13 )//this is positron on muon+
			{
				flag++;
				//cout<<"Plus Sign lepton with pdgId: "<<(*particleId)[l]<<endl;
				p4leptonPlus.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);	
			}


		}
				
		if (flag==2)
		{
		//	cout<<"--------------------------------------------------------------------"<<endl;
			isDilepton = true;
			dilepton++;
		//	cout<<"Dilepton Found!"<<endl;
		}
		else isDilepton =false;
		
		TVector3 ttbarBoostVector = getBoostVector(p4T, p4Tbar, p4TTbar);
		
		//clone all the lorentz vectors that I will neeed
		//top /antitop in ZMF
		p4T_ttbarZMF.SetPtEtaPhiM(p4T.Pt(), p4T.Eta(), p4T.Phi(), p4T.M());
		p4Tbar_ttbarZMF.SetPtEtaPhiM(p4Tbar.Pt(), p4Tbar.Eta(), p4Tbar.Phi(), p4Tbar.M());
		p4T_ttbarZMF.Boost(ttbarBoostVector);
		p4Tbar_ttbarZMF.Boost(ttbarBoostVector);
		if(isDilepton)
		{
			//leptons in ZMF
			p4leptonMinus_ttbarZMF.SetPtEtaPhiM(p4leptonMinus.Pt(), p4leptonMinus.Eta(), p4leptonMinus.Phi(), p4leptonMinus.M());
			p4leptonPlus_ttbarZMF.SetPtEtaPhiM(p4leptonPlus.Pt(), p4leptonPlus.Eta(), p4leptonPlus.Phi(), p4leptonPlus.M());
			p4leptonMinus_ttbarZMF.Boost(ttbarBoostVector);
			p4leptonPlus_ttbarZMF.Boost(ttbarBoostVector);
			
			//leptons in top / antitop rest frame
			leptonMinus_tbarZMF.SetPtEtaPhiM(p4leptonMinus.Pt(), p4leptonMinus.Eta(), p4leptonMinus.Phi(), p4leptonMinus.M()); 
			leptonPlus_topZMF.SetPtEtaPhiM(p4leptonPlus.Pt(), p4leptonPlus.Eta(), p4leptonPlus.Phi(), p4leptonPlus.M());
			
			TVector3 topBoostVector  = p4T.BoostVector();
			TVector3 tbarBoostVector = p4Tbar.BoostVector();
			leptonMinus_tbarZMF.Boost(-tbarBoostVector);
			leptonPlus_topZMF.Boost(-topBoostVector);
			
			//now that I am in lab frame I will plot deltaPhi, deltaEta
			float pi = 3.14;
			h_deltaPhi->Fill(p4leptonMinus.DeltaPhi(p4leptonPlus)/TMath::Pi());
			h_deltaEta->Fill(fabs(TMath::TanH((p4leptonMinus.Eta()-p4leptonPlus.Eta())/2.)));
			//_deltaEta->Fill(fabs(p4leptonMinus.Eta()- p4leptonPlus.Phi()));
			
			
			_coscos = TMath::Cos(leptonMinus_tbarZMF.Angle(p4Tbar_ttbarZMF.Vect())) * TMath::Cos(leptonPlus_topZMF.Angle(p4T_ttbarZMF.Vect()));
			h_coscos->Fill(_coscos);
			
			_cosphi = TMath::Cos(leptonMinus_tbarZMF.Angle(leptonPlus_topZMF.Vect()));
			h_cos->Fill(_cosphi);
		}
	  
  }
 cout<<"number of Top pairs: "<<numOfEntries<<endl;
 cout<<"Number of lepton pairs: "<<dilepton<<endl; 

  TCanvas *can1 = new TCanvas("can1","can1", 900, 600);
  h_coscos->Draw();

  TCanvas *can2 = new TCanvas("can2","can2", 900, 600);
  h_cos->Draw();

  TCanvas *can3 = new TCanvas("can3","can3", 900, 600);
  h_deltaEta->Draw();

  TCanvas *can4 = new TCanvas("can4","can4", 900, 600);
  h_deltaPhi->Draw();
  
  TFile *outf = new TFile("TopAngularKinematicsFile.root", "UPDATE");
  outf->cd();
  h_deltaEta->Write("h_deltaEta_SC");
  h_deltaPhi->Write("h_deltaPhi_SC");
  h_cos->Write("h_cos_SC");
  h_coscos->Write("h_coscos_SC");
  /*
  TCanvas *can5 = new TCanvas("can5","can5", 900, 600);
  h_eta->Draw();
  
  TCanvas *can6 = new TCanvas("can6","can6", 900, 600);
  h_pt->Draw();
      
  TCanvas *can7 = new TCanvas("can7","can7", 900, 600);
  h_mass->Draw();
  
  TCanvas *can8 = new TCanvas("can8","can8", 900, 600);
  h_phi->Draw();
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
