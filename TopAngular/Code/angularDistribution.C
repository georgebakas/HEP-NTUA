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

TVector3 getBoostVector(TLorentzVector p4_1, TLorentzVector p4_2, TLorentzVector &p4CombinedVector);


void angularDistribution(TString fileName = "/eos/cms/store/user/gbakas/ZprimeToTT/ZprimeToTT_M-2000_W-200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_2016_Copy.root", 
						float selMvaCut=0.3, float floatBTag = 0.8838)
{
  //gStyle->SetOptStat(0);
  TFile *inf     = TFile::Open(fileName);
  TTree *trIN    = (TTree*)inf->Get("events");
  
  float XSEC = 832.;
  float NORM = ((TH1F*)inf->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
  float weight = XSEC/NORM;
  int nJets,nLeptons;
  float genEvtWeight;
  vector<float> *jetPt(0),*tau3(0),*tau2(0),*tau1(0);
  vector<float> *jetMassSub0(0), *jetMassSub1(0);
  vector<float> *jetMassSoftDrop(0);

  float mva(0);
  vector<float> *jetTtag(0);
  vector<bool> *bit(0);
  float mTTbarParton(0),mJJ(0), yTTbarParton(0), ptTTbarParton(0);
  int  category(0);
  //matching info
  vector<float> *jetPhi(0), *jetEta(0);
  vector<int> *partonId(0), *partonMatchIdx(0);
  
  vector<float> *partonEta(0), *partonPhi(0), *partonMatchDR(0),  *partonPt(0), *partonE(0), *partonMass(0);
  std::vector<int> *addedIndexes = new std::vector<int>(0);
  std::vector<float> *jetBtagSub0(0), *jetBtagSub1(0);
          
  //------- input tree --------------
  trIN->SetBranchAddress("nJets"          ,&nJets);
  trIN->SetBranchAddress("nLeptons"       ,&nLeptons);
  trIN->SetBranchAddress("jetPt"          ,&jetPt);
  trIN->SetBranchAddress("jetEta"         ,&jetEta);
  trIN->SetBranchAddress("jetPhi"         ,&jetPhi);
  trIN->SetBranchAddress("jetTau3"        ,&tau3);
  trIN->SetBranchAddress("jetTau2"        ,&tau2);
  trIN->SetBranchAddress("jetTau1"        ,&tau1);
  trIN->SetBranchAddress("triggerBit"     ,&bit);
  trIN->SetBranchAddress("genEvtWeight"   ,&genEvtWeight);
  trIN->SetBranchAddress("jetMassSub0"    ,&jetMassSub0);
  trIN->SetBranchAddress("jetMassSub1"    ,&jetMassSub1);
  trIN->SetBranchAddress("jetBtagSub0"	  ,&jetBtagSub0);
  trIN->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1);
  trIN->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
  trIN->SetBranchAddress("mJJ"   		  ,&mJJ);
  trIN->SetBranchAddress("mva"	  		  ,&mva);
  trIN->SetBranchAddress("category"	  	  ,&category);
  trIN->SetBranchAddress("partonMatchIdx" ,&partonMatchIdx);
  trIN->SetBranchAddress("partonMatchDR"  ,&partonMatchDR);
  trIN->SetBranchAddress("jetTtagCategory",&jetTtag);
  //parton variables
  trIN->SetBranchAddress("mTTbarParton"   ,&mTTbarParton);
  trIN->SetBranchAddress("yTTbarParton"   ,&yTTbarParton);
  trIN->SetBranchAddress("ptTTbarParton"  ,&ptTTbarParton);
  trIN->SetBranchAddress("partonPt"		  ,&partonPt);
  trIN->SetBranchAddress("partonEta"	  ,&partonEta);
  trIN->SetBranchAddress("partonPhi" 	  ,&partonPhi);
  trIN->SetBranchAddress("partonE"	 	  ,&partonE);
  trIN->SetBranchAddress("partonMass"	  ,&partonMass);

  
  
  
  TLorentzVector p4T[2], p4TTbar, p4T_ZMF[2];

  int decade(0);
  int NN = trIN->GetEntries();

  //int NN = 10000;
  const int sizeBins = 6;
  float BND[sizeBins+1] = {800,1000, 1200, 1500, 2000, 2500, 3000};
  int counter =0;
  
  TH1F *h_mTTbarParton = new TH1F("mTTbarParton", "mTTbarParton histogram", 50, 800,3000);
  
  TH1F *hAngularDist[sizeBins];
  std::vector<TString> massLimits = {"800-100","1000-1200","1200-1500", "1500-2000", "2000-2500", "2500-3000"};
  
  for(int i=0; i<=(sizeBins-1); i++)
  {
	  TString temp = massLimits[i];
	  hAngularDist[i] = new TH1F(TString::Format("cos(#theta) Distribution for mass limit: %s (GeV)", temp.Data()),TString::Format("cos(#theta) Distribution for mass limit: %s (GeV)", temp.Data()), 20, 0,1);
  }
  
  
  cout<<"Reading "<<NN<<" entries"<<endl;
  for(int iev=0;iev<NN;iev++) 
  {
    double progress = 10.0*iev/(1.0*NN);
    int k = TMath::FloorNint(progress); 
    if (k > decade) 
      cout<<10*k<<" %"<<endl;
    decade = k;
    trIN->GetEntry(iev);
			
	
	bool recoCuts   = nJets > 1  && fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1] <2.4) &&(*jetPt)[0] > 400 && (*jetPt)[1] > 400;
	bool partonCuts = (*partonPt)[0] > 400 && (*partonPt)[1] > 400 && fabs((*partonEta)[0]) < 2.4 && fabs((*partonEta)[1]) < 2.4 &&  mTTbarParton > 800;
	bool btagging   = ((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag) && ((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag);
	bool topTagger  = (*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut;
	//cout<<"ok"<<endl;
	if(recoCuts && partonCuts && topTagger && btagging)
	{
		h_mTTbarParton->Fill(mTTbarParton);
		//split the mttbar phase space into regions and for each region I will calculate the thetas dists
		p4T[0].SetPtEtaPhiM((*partonPt)[0], (*partonEta)[0], (*partonPhi)[0], (*partonMass)[0]);
		p4T[1].SetPtEtaPhiM((*partonPt)[1], (*partonEta)[1], (*partonPhi)[1], (*partonMass)[1]);
		
		TVector3 ttbarBoostVector = getBoostVector(p4T[0], p4T[1], p4TTbar);
		
		p4T_ZMF[0].SetPtEtaPhiM(p4T[0].Pt(), p4T[0].Eta(), p4T[0].Phi(), p4T[0].M());
		p4T_ZMF[1].SetPtEtaPhiM(p4T[1].Pt(), p4T[1].Eta(), p4T[1].Phi(), p4T[1].M());
		p4T_ZMF[0].Boost(ttbarBoostVector);
		p4T_ZMF[1].Boost(ttbarBoostVector);
				
		//cout<<"-------------------------------------"<<endl;		
		//cout<< p4T_ZMF[0].Pt()<<endl;
		//cout<< p4T_ZMF[1].Pt()<<endl;
		
		if(mTTbarParton >= 800 && mTTbarParton < 1000 )
		{
			hAngularDist[0]->Fill(TMath::Cos(p4T_ZMF[0].Theta()));
			hAngularDist[0]->Fill(TMath::Cos(p4T_ZMF[1].Theta()));
		}
		else if(mTTbarParton >= 1000 && mTTbarParton < 1200)
		{
			hAngularDist[1]->Fill(TMath::Cos(p4T_ZMF[0].Theta()));
			hAngularDist[1]->Fill(TMath::Cos(p4T_ZMF[1].Theta()));
		}
		else if(mTTbarParton >= 1200 && mTTbarParton < 1500)
		{
			hAngularDist[2]->Fill(TMath::Cos(p4T_ZMF[0].Theta()));
			hAngularDist[2]->Fill(TMath::Cos(p4T_ZMF[1].Theta()));
		}
		else if(mTTbarParton >= 1500 && mTTbarParton < 2000)
		{
			hAngularDist[3]->Fill(TMath::Cos(p4T_ZMF[0].Theta()));
			hAngularDist[3]->Fill(TMath::Cos(p4T_ZMF[1].Theta()));
		}
		else if(mTTbarParton >= 2000 && mTTbarParton < 2500)
		{
			hAngularDist[4]->Fill(TMath::Cos(p4T_ZMF[0].Theta()));
			hAngularDist[4]->Fill(TMath::Cos(p4T_ZMF[1].Theta()));
		}
		else if(mTTbarParton >= 2500 && mTTbarParton < 3000)
		{
			hAngularDist[5]->Fill(TMath::Cos(p4T_ZMF[0].Theta()));
			hAngularDist[5]->Fill(TMath::Cos(p4T_ZMF[1].Theta()));
		}
	}
  }	
  
  TCanvas *can[sizeBins];
  for(int i =0; i<=(sizeBins-1); i++)
  {
	can[i] = new TCanvas(TString::Format("can %d", (i+1)), TString::Format("can %d", (i+1)), 900, 600);
	hAngularDist[i]->GetXaxis()->SetTitle("cos(#theta)");
	hAngularDist[i]->Draw();
  }
  
  TCanvas *can_mTTbarParton = new TCanvas("mTTbarParton can", "mTTbarParton can", 900, 600);
  h_mTTbarParton->GetXaxis()->SetTitle("mTTbarParton (GeV)");
  h_mTTbarParton->Draw();
  
  
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

