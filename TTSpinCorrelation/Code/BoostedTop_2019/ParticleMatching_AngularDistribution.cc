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
#include <TEfficiency.h>


#include "ParticleMatching_AngularDistribution.h"
using namespace std;

TVector3 getBoostVector(TLorentzVector p4_1, TLorentzVector p4_2, TLorentzVector &p4CombinedVector);


void ParticleMatching_AngularDistribution()
{
	

  gStyle->SetOptStat(0);
  //reco angular distributions
  //TString inputfile = "/afs/cern.ch/work/g/gbakas/private/TTSpinCorrelation/boostedTTBar_SpCorr/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_BoostedPart.root";
  TString inputfile = "/eos/cms/store/user/gbakas/spinCorrelation/mc-2018/TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root";
  int size = 20;
  TFile *inf = new TFile(inputfile, "READ"); 
  TTree *trIn = (TTree*)inf->Get("boosted/events");
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
  
  trIn->SetBranchAddress("triggerBit", &triggerBits);
  trIn->SetBranchAddress("category", &category);
  trIn->SetBranchAddress("nJets", &nJets);
  trIn->SetBranchAddress("ht", &ht);
  trIn->SetBranchAddress("mva", &mva);
  trIn->SetBranchAddress("jetBtagSub0", &jetBtagSub0);
  trIn->SetBranchAddress("jetBtagSub1", &jetBtagSub1);
  trIn->SetBranchAddress("jetMassSoftDrop", &jetMassSoftDrop);
  trIn->SetBranchAddress("nvtx", &nvtx);
  
  trIn->SetBranchAddress("jetPt", &jetPt);
  trIn->SetBranchAddress("jetEta", &jetEta);
  trIn->SetBranchAddress("jetPhi", &jetPhi);
  trIn->SetBranchAddress("jetMass", &jetMass);
  
  trIn->SetBranchAddress("jetMassSub0", &jetMassSub0);
  trIn->SetBranchAddress("jetMassSub1", &jetMassSub1);
  trIn->SetBranchAddress("jetPtSub0", &jetPtSub0);
  trIn->SetBranchAddress("jetPtSub1", &jetPtSub1);
  trIn->SetBranchAddress("jetPhiSub0", &jetPhiSub0);
  trIn->SetBranchAddress("jetPhiSub1", &jetPhiSub1);
  trIn->SetBranchAddress("jetEtaSub0", &jetEtaSub0);
  trIn->SetBranchAddress("jetEtaSub1", &jetEtaSub1);
	
  vector <float>*particleMass(0);
  vector <float>*particlePt(0);
  vector <float>*particleEta(0);
  vector <float>*particlePhi(0);
  //this is particle id
  vector <int>  *particleId(0);

  //declare the branches
  TBranch *b_eta, *b_pt, *b_phi, *b_mass, *b_id;
  trIn->SetBranchAddress("particleId",&particleId, &b_id);
  trIn->SetBranchAddress("particleMass",&particleMass, &b_mass);
  trIn->SetBranchAddress("particlePt",&particlePt, &b_pt);
  trIn->SetBranchAddress("particleEta",&particleEta, &b_eta);
  trIn->SetBranchAddress("particlePhi",&particlePhi, &b_phi);
  

  
										
  TLorentzVector p4T_0(0,0,0,0), p4T_1(0,0,0,0), p4W_0(0,0,0,0), p4W_1(0,0,0,0);
  TLorentzVector p4TTbar(0,0,0,0);
  cout<<"starting..."<<endl;
  
  TLorentzVector p4T_parton(0,0,0,0), p4Tbar_parton(0,0,0,0), p4Wminus_parton(0,0,0,0), p4Wplus_parton(0,0,0,0);
  bool isTopPtOk = 0;
  cout<<"Starting with parton level...."<<endl;
  
  
  pt_Jet    = new std::vector<float>;  
  eta_Jet   = new std::vector<float>;  
  phi_Jet   = new std::vector<float>;  
  mass_Jet  = new std::vector<float>; 

  pt_SubJet0   = new std::vector<float>;  
  pt_SubJet1   = new std::vector<float>;
  eta_SubJet0  = new std::vector<float>;  
  eta_SubJet1  = new std::vector<float>;
  phi_SubJet0  = new std::vector<float>;  
  phi_SubJet1  = new std::vector<float>;
  mass_SubJet0 = new std::vector<float>;  
  mass_SubJet1 = new std::vector<float>;
  
  matched_jetBtagSub0 = new std::vector<float>;
  matched_jetBtagSub1  = new std::vector<float>;
  
  
  partonMatchIdx = new std::vector<int>;  
  partonMatchDR  = new std::vector<float>;  
  
  int numOfJets(0);
  float dRmin_all(0.4);
  
//------------------------------------------------------HISTOGRAMS------------------------------------------------------------------------------------------------------------------
  TH1F* dPhi_Reco = new TH1F("W candidate #Delta#phi Reco" , "W candidate #Delta#phi Reco", size, -3, 3); 
  TH1F* dPhi_Parton = new TH1F("W #Delta#phi Parton" , "W #Delta#phi Parton", size, -3, 3); 

  TH1F* h_PtParton[2];
  TH1F* h_PhiParton[2];
  TH1F* h_EtaParton[2];
  TH1F* h_MassParton[2];
  h_PtParton[0]   = new TH1F("P_{T} Top parton","P_{T} Top parton", size, 400,1000);
  h_PtParton[1]   = new TH1F("P_{T} W candidate parton","P_{T} W candidate parton", size, 400,1000);
  h_PhiParton[0]  = new TH1F("#phi Top parton","#phi Top parton", size, -3,3);
  h_PhiParton[1]  = new TH1F("#phi W candidate parton","#phi W candidate parton", size, -3,3);
  h_EtaParton[0]  = new TH1F("#eta Top  parton","#eta Top parton", size, -3,3);
  h_EtaParton[1]  = new TH1F("#eta W candidate parton","#eta W candidate parton", size, -3,3);
  h_MassParton[0] = new TH1F("Mass Top parton","Mass Top parton", size, 0,500);
  h_MassParton[1] = new TH1F("Mass parton W candidate","Mass  W candidate parton", size, 0,500);
  
  
  TH1F* h_PtReco[2];
  TH1F* h_EtaReco[2];
  TH1F* h_PhiReco[2];
  TH1F* h_MassReco[2];
  h_PtReco[0]     = new TH1F("P_{T} Top reco","P_{T} Top reco", size, 400,1000);
  h_PtReco[1]     = new TH1F("P_{T} W candidate reco","P_{T} W candidate reco", size, 400,1000);
  h_PhiReco[0] = new TH1F("#phi Top reco","#phi Top reco", size, -3,3);
  h_PhiReco[1] = new TH1F("#phi W candidate reco","#phi W candidate reco", size, -3,3);
  h_EtaReco[0] = new TH1F("#eta Top reco","#eta Top reco", size, -3,3);
  h_EtaReco[1] = new TH1F("#eta W candidate reco","#eta W candidate reco", size, -3,3);
  h_MassReco[0] = new TH1F("Mass Top reco","Mass Top reco", size, 0,500);
  h_MassReco[1] = new TH1F("Mass W candidate reco","Mass W candidate reco", size, 0,500);
  
  TH1F* h_massSoftDrop = new TH1F("Jet Mass Soft Drop", "Jet Mass Soft Drop", size, 0, 1500);
  TH1F* h_mass = new TH1F("Jet Mass", "Jet Mass", size, 0, 1500);
//------------------------------------------------------END OF HISTOGRAMS------------------------------------------------------------------------------------------------------------

    cout<<"ok"<<endl;
for(int i=0; i<trIn->GetEntries(); i++)
{
	 
	trIn->GetEntry(i);    
	if((*triggerBits)[triggerBit] && category == 2  && nJets > 1 && mva > 0.8 && nvtx >= 15 && nvtx <= 25 && (*jetPt)[1] > 350)
    {
	
	    for(std::vector<float>::iterator ptIter = jetPt->begin(); ptIter!=jetPt->end(); ++ptIter)
		{
				int index = ptIter - jetPt->begin();
				//if(index >2 ) cout<<index<<endl;
				if((*jetMassSoftDrop)[index] >= 150 && (*jetMassSoftDrop)[index] <= 210 && (*jetPt)[index] > 350)
				{ 
				  numOfJets++;
				  pt_Jet  ->push_back((*jetPt)[index]);
				  eta_Jet ->push_back((*jetEta)[index]);
				  phi_Jet ->push_back((*jetPhi)[index]);
				  mass_Jet->push_back((*jetMassSoftDrop)[index]);
				   //I also need the subjets pushed back....
				  pt_SubJet0   ->push_back((*jetPtSub0)[index]);
				  pt_SubJet1   ->push_back((*jetPtSub1)[index]);
				  eta_SubJet0  ->push_back((*jetEtaSub0)[index]);
				  eta_SubJet1  ->push_back((*jetEtaSub1)[index]);
				  phi_SubJet0  ->push_back((*jetPhiSub0)[index]);
				  phi_SubJet1  ->push_back((*jetPhiSub1)[index]);
				  mass_SubJet0 ->push_back((*jetMassSub0)[index]);
				  mass_SubJet1 ->push_back((*jetMassSub1)[index]);
				  matched_jetBtagSub0 ->push_back((*jetBtagSub0)[index]);
				  matched_jetBtagSub1 ->push_back((*jetBtagSub1)[index]);
				  
				}				
				//cout<<numOfJets<<endl;
		} 
    
		//cout<<numOfJets<<endl;
		//---------------------------------------------------------------- parton level ------------------------------------------------------------------------------------------------
		//try the particle matching....
		for(int l=0; l<=5; l++)
		{
				
			if (fabs((*particleId)[l]) == 6) {
					
						int bothW = 0;
			//----- match partons with jets ------------
				float dRmin(1000);
				int imatch(-1);
				for(int j=0;j<numOfJets;j++) {
				  
				  float dR = TMath::Sqrt( TMath::Power((*particleEta)[l]-(*eta_Jet)[j],2) +TMath::Power((*particlePhi)[l]-(*phi_Jet)[j],2));
				  if (dR < dRmin) {
					imatch = j;
					dRmin = dR;
				  }
				}
				
			
				if(dRmin < dRmin_all)
				{
					
					//partonMatchIdx->push_back(imatch);
					//partonMatchDR->push_back(dRmin);
					//now I have the matched particles so I will fill the histograms
					h_PtParton[0]  ->Fill((*particlePt)[l]);
					h_PtReco[0]    ->Fill((*pt_Jet)[imatch]);
					
					h_EtaParton[0] ->Fill((*particleEta)[l]);
					h_EtaReco[0]   ->Fill((*eta_Jet)[imatch]);
					
					h_PhiParton[0] ->Fill((*particlePhi)[l]);
					h_PhiReco[0]   ->Fill((*phi_Jet)[imatch]);
					//out<<"ok"<<endl;
					 
					//cout<<particleId->size()<<endl;
					for(int k = 2; k<=5; k++)
					{	 
						
						
						if(fabs((*particleId)[k])== 24 && (*particlePt)[k] > 400)
						{
							
							//W matching
							
							//1st pass selection cuts and then check the dR to match
							if((*pt_SubJet0)[imatch] > 400 && (*matched_jetBtagSub0)[imatch] < 0.5426)
							{
								float dR_SubJet0 = TMath::Sqrt( TMath::Power((*particleEta)[k]-(*eta_SubJet0)[imatch],2) +TMath::Power((*particlePhi)[k]-(*phi_SubJet0)[imatch],2));
								//if(dR_SubJet0 < 0.4)
								//{
									h_PtParton[1]  ->Fill((*particlePt)[k]);
									h_PtReco[1]    ->Fill((*pt_SubJet0)[imatch]);
				
									h_EtaParton[1] ->Fill((*particleEta)[k]);
									h_EtaReco[1]   ->Fill((*eta_SubJet0)[imatch]);
							
									h_PhiParton[1] ->Fill((*particlePhi)[k]);
									h_PhiReco[1]   ->Fill((*phi_SubJet0)[imatch]);
								//}
							}
							if((*pt_SubJet1)[imatch] > 400 && (*matched_jetBtagSub1)[imatch] < 0.5426)
							{
								float dR_SubJet1 = TMath::Sqrt( TMath::Power((*particleEta)[k]-(*eta_SubJet1)[imatch],2) +TMath::Power((*particlePhi)[k]-(*phi_SubJet1)[imatch],2));
								//if(dR_SubJet1 < 0.4)
								//{
									h_PtParton[1]  ->Fill((*particlePt)[k]);
									h_PtReco[1]    ->Fill((*pt_SubJet1)[imatch]);
										
									h_EtaParton[1] ->Fill((*particleEta)[k]);
									h_EtaReco[1]   ->Fill((*eta_SubJet1)[imatch]); 
										
									h_PhiParton[1] ->Fill((*particlePhi)[k]);
									h_PhiReco[1]   ->Fill((*phi_SubJet1)[imatch]);
								//}
							}
							
						} 
						
						
					}
				}
				
			  } 
			  
			  //now that the reco/parton is matched, I need to see their distributions
			  
			  
		}
	}
	partonMatchIdx->clear();
	partonMatchDR ->clear();
	pt_Jet   ->clear();  
	eta_Jet  ->clear();
	phi_Jet  ->clear();
	mass_Jet ->clear(); 
	numOfJets = 0;
	pt_SubJet0   ->clear();
	pt_SubJet1   ->clear();
	eta_SubJet0  ->clear();
	eta_SubJet1  ->clear();
    phi_SubJet0  ->clear();
	phi_SubJet1  ->clear();
	mass_SubJet0 ->clear();
	mass_SubJet1 ->clear();
	matched_jetBtagSub0 ->clear();
	matched_jetBtagSub1 ->clear();
}
  
	  
		
 for(int i =0; i<2; i++)
 {	 
	 double norm =1;
	 //pt distributions
	 Double_t scale_1 = norm/(h_PtParton[i]->Integral());
	 h_PtParton[i]->Scale(scale_1);
	 Double_t scale_2 = norm/(h_PtReco[i]->Integral());
	 h_PtReco[i]->Scale(scale_2);
	 
	 h_PtParton[i]->GetXaxis()->SetTitle("P_{T} (GeV)");
	 //h_PtParton[i]->GetYaxis()->SetTitle("Events");
	 h_PtReco[i]->GetXaxis()->SetTitle("P_{T} (GeV)");
	 //h_PtReco[i]->GetYaxis()->SetTitle("Events");
	 h_PtParton[i]->SetLineColor(kRed);
	 h_PtReco[i]->SetLineColor(kBlue);
	 
	 //phi distributions
	 Double_t scale_phi1 = norm/(h_PhiParton[i]->Integral());
	 h_PhiParton[i]->Scale(scale_phi1);
	 Double_t scale_phi2 = norm/(h_PhiReco[i]->Integral());
	 h_PhiReco[i]->Scale(scale_phi2);
	 
	 h_PhiParton[i]->SetLineColor(kRed);
	 h_PhiReco[i]->SetLineColor(kBlue);
     h_PhiParton[i]->GetXaxis()->SetTitle("#phi");
	 h_PhiReco[i]->GetXaxis()->SetTitle("#phi");
	 
	 //eta distributions 
	 Double_t scale_eta1 = norm/(h_EtaParton[i]->Integral());
	 h_EtaParton[i]->Scale(scale_eta1);
	 Double_t scale_eta2 = norm/(h_EtaReco[i]->Integral());
	 h_EtaReco[i]->Scale(scale_eta2);
	 
	 h_EtaParton[i]->SetLineColor(kRed);
	 h_EtaReco[i]->SetLineColor(kBlue);
	 h_EtaParton[i]->GetXaxis()->SetTitle("#eta");
	 h_EtaReco[i]->GetXaxis()->SetTitle("#eta");
 }

 
  TCanvas *canPtTop = new TCanvas("Pt Top", "Pt Top", 900,600);
  TLegend *leg_pt = new TLegend(0.86,0.7,0.99,0.9);
  leg_pt->AddEntry(h_PtParton[0], "h_PtParton", "l");
  leg_pt->AddEntry(h_PtReco[0], "h_PtReco", "l"); 
  h_PtParton[0]->Draw();
  h_PtReco[0]->Draw("same");
  leg_pt->Draw();
  
  
  TCanvas *canPtWCandidate = new TCanvas("Pt WCandidate", "Pt WCandidate", 900,600);
  h_PtParton[1]->Draw();
  h_PtReco[1]->Draw("same");
  leg_pt->Draw();

  //phi distributions for reco and parton
  TCanvas *canPhiTop = new TCanvas("#phi Top can", "#phi Top can", 900, 600);
  TLegend *legPhi = new TLegend(0.86,0.7,0.99,0.9);
  h_PhiParton[0]->Draw();
  h_PhiReco[0]->Draw("same");
  legPhi->AddEntry(h_PhiParton[0],"#phi parton","l");
  legPhi->AddEntry(h_PhiReco[0],"#phi reco","l"); 
  legPhi->Draw();
  
  TCanvas *canPhiW = new TCanvas("#phi WCandidate can", "#phi WCandidate can", 900, 600);
  h_PhiParton[1]->Draw();
  h_PhiReco[1]->Draw("same");
  legPhi->Draw();
  
  //eta distributions  
  TCanvas *canEtaTop = new TCanvas("#eta Top can", "#eta Top can", 900, 600);
  TLegend *legEta = new TLegend(0.86,0.7,0.99,0.9);
  h_EtaParton[0]->Draw();
  h_EtaReco[0]->Draw("same");
  legEta->AddEntry(h_EtaParton[0],"#eta parton","l");
  legEta->AddEntry(h_EtaReco[0],"#eta reco","l"); 
  legEta->Draw();
  
  TCanvas *canEtaWCandidate = new TCanvas("#eta WCandidate can", "#eta WCandidate can", 900, 600);
  h_EtaParton[1]->Draw();
  h_EtaReco[1]->Draw("same");
  legEta->Draw();
  
  
//  TEfficiency *eff_topPt= 0, *eff_topPhi = 0, *eff_topEta = 0;
  
 // if(TEfficiency::CheckConsistency(*h_PtParton[0],*h_PtReco[0]))
//{
 // TH1F* eff_topPt = new TH1F("pt div", "pt div", size, 400,1000);
 // TH1F* eff_topPhi = new TH1F("#phi div","#phi div", size, -3,3);
 // TH1F* eff_topEta = new TH1F("#eta div", "eta div", size, -3,3);

  TH1F* eff_topPt  = (TH1F*)h_PtReco[0]->Clone("eff_topPt");
  TH1F* eff_topPhi = (TH1F*)h_PhiReco[0]->Clone("eff_topPhi");
  TH1F* eff_topEta = (TH1F*)h_EtaReco[0]->Clone("eff_topEta");

  eff_topPt ->Divide(h_PtParton[0]);
  eff_topPhi->Divide(h_PhiParton[0]);
  eff_topEta->Divide(h_EtaParton[0]);
  // this will write the TEfficiency object to "myfile.root"
  // AND pEff will be attached to the current directory

  
  TCanvas *can_effPt  = new TCanvas("Pt Efficiency can", "Pt efficienct", 900, 600);
  eff_topPt->Draw();
  
  TCanvas *can_effPhi = new TCanvas("#phi eff can", "#phi eff can", 900, 600);
  eff_topPhi->Draw();
  
  TCanvas *can_effEta = new TCanvas("#eta eff can", "#eta eff can", 900, 600);
  eff_topEta->Draw();
//}
  /*
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
 
  
  */
 
  delete partonMatchIdx;
  delete partonMatchDR;
  delete pt_Jet;  
  delete eta_Jet;
  delete phi_Jet;
  delete mass_Jet; 
  delete pt_SubJet0;
  delete pt_SubJet1;
  delete eta_SubJet0;
  delete eta_SubJet1;
  delete phi_SubJet0;
  delete phi_SubJet1;
  delete mass_SubJet0;
  delete mass_SubJet1;
  delete matched_jetBtagSub0;
  delete matched_jetBtagSub1;
 
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


