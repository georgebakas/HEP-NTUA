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


void angularDistribution(TString file = "/eos/cms/store/user/gbakas/ZprimeToTT/mc/2017/ZprimeToTT_M3000_W300_TuneCP2_13TeV-madgraphMLM-pythia8_Copy.root", 
						float selMvaCut=0.3, float floatBTag = 0.8838, bool isZprime= true,bool isParton=true, int ZprimeMass = 2000, TString width = "200" )
{
	
//TString TTbarFile = "/eos/cms/store/user/gbakas/ttbar/topTagger/April19/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_Copy.root"	
  //gStyle->SetOptStat(0);
  TFile *inf     = TFile::Open(file);
  TTree *trIN    = (TTree*)inf->Get("events");
  //cout<<"here"<<endl;
  
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
  const int sizeBins = 1;
  //float BND[sizeBins+1] = {1000, 2000, 3000, 4000, 5000};
  //float BND[sizeBins+1] = {1000, 2500, 3500, 5000, 6000};
  float BND[sizeBins+1] = {1000, 6000};
  int counter =0;
  
  TH1F *h_mTTbarParton;
  if(isZprime)h_mTTbarParton  = new TH1F("mTTbarParton", "mTTbarParton histogram", 60, 1000,ZprimeMass+1000);
  else h_mTTbarParton = new TH1F("mTTbarParton", "mTTbarParton histogram", 60, 1000,6000);  
  
  TH1F *hAngularDist[sizeBins], *hChi[sizeBins];
  std::vector<TString> massLimits = {"1000-2500","2500-3500", "3500-5000", "5000-Inf"};
  
  for(int i=0; i<=(sizeBins-1); i++)
  {
	  TString temp = massLimits[i];
	  hAngularDist[i] = new TH1F(TString::Format("#||{cos(#theta)} Distribution for mass limit: %s (GeV)", temp.Data()),TString::Format("#||{cos(#theta)} Distribution for mass limit: %s (GeV)", temp.Data()), 20, 0,1);
	  hChi[i] = new TH1F(TString::Format("#chi Distribution for mass limit: %s (GeV)", temp.Data()),TString::Format("#chi Distribution for mass limit: %s (GeV)", temp.Data()), 20, 0,16);
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
	bool recoCuts, partonCuts, btagging, topTagger, massCut;	
	if (nJets >1)
	{
		recoCuts   = fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1] <2.4) &&(*jetPt)[0] > 400 && (*jetPt)[1] > 400 && nLeptons==0;
		partonCuts = (*partonPt)[0] > 400 && (*partonPt)[1] > 400 && fabs((*partonEta)[0]) < 2.4 && fabs((*partonEta)[1]) < 2.4 &&  mTTbarParton > 1000;
		btagging   = ((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag) && ((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag);
		topTagger  = (*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut;
		massCut    = (*jetMassSoftDrop)[0] > 120 && (*jetMassSoftDrop)[0] < 220 && (*jetMassSoftDrop)[1] > 120 && (*jetMassSoftDrop)[1] < 220;
		//cout<<"ok"<<endl;
		if(recoCuts && partonCuts && topTagger && btagging)
		{
			h_mTTbarParton->Fill(mTTbarParton);
			//split the mttbar phase space into regions and for each region I will calculate the thetas dists
			if(isParton)
			{
				p4T[0].SetPtEtaPhiM((*partonPt)[0], (*partonEta)[0], (*partonPhi)[0], (*partonMass)[0]);
				p4T[1].SetPtEtaPhiM((*partonPt)[1], (*partonEta)[1], (*partonPhi)[1], (*partonMass)[1]);
			}
			else
			{
				p4T[0].SetPtEtaPhiM((*jetPt)[0], (*jetEta)[0], (*jetPhi)[0], (*jetMassSoftDrop)[0]);
				p4T[1].SetPtEtaPhiM((*jetPt)[1], (*jetEta)[1], (*jetPhi)[1], (*jetMassSoftDrop)[1]);
			}
			
			TVector3 ttbarBoostVector = getBoostVector(p4T[0], p4T[1], p4TTbar);
			
			p4T_ZMF[0].SetPtEtaPhiM(p4T[0].Pt(), p4T[0].Eta(), p4T[0].Phi(), p4T[0].M());
			p4T_ZMF[1].SetPtEtaPhiM(p4T[1].Pt(), p4T[1].Eta(), p4T[1].Phi(), p4T[1].M());
			p4T_ZMF[0].Boost(ttbarBoostVector);
			p4T_ZMF[1].Boost(ttbarBoostVector);
					
			//cout<<"-------------------------------------"<<endl;		
			//cout<< p4T_ZMF[0].Pt()<<endl;
			//cout<< p4T_ZMF[1].Pt()<<endl;
			float chi0(0), chi1(0);
			chi0 = (1 + fabs(TMath::Cos(p4T_ZMF[0].Theta()))) / ( 1 - fabs(TMath::Cos(p4T_ZMF[0].Theta())));
			chi1 = (1 + fabs(TMath::Cos(p4T_ZMF[1].Theta()))) / ( 1 - fabs(TMath::Cos(p4T_ZMF[1].Theta())));
			//chi0 = TMath::Exp(p4T_ZMF[0].Rapidity() - p4T_ZMF[1].Rapidity());
			
			
			if(mTTbarParton >= BND[0] && mTTbarParton < BND[1] )
			{
				hChi[0] ->Fill(chi0);
				hChi[0] ->Fill(chi1);
				hAngularDist[0]->Fill(fabs(TMath::Cos(p4T_ZMF[0].Theta())));
				hAngularDist[0]->Fill(fabs(TMath::Cos(p4T_ZMF[1].Theta())));
			}
			/*
			else if(mTTbarParton >= BND[1] && mTTbarParton < BND[2])
			{
				hChi[1] ->Fill(chi0);
				hChi[1] ->Fill(chi1);
				hAngularDist[1]->Fill(fabs(TMath::Cos(p4T_ZMF[0].Theta())));
				hAngularDist[1]->Fill(fabs(TMath::Cos(p4T_ZMF[1].Theta())));
			}
			else if(mTTbarParton >= BND[2] && mTTbarParton < BND[3])
			{
				hChi[2] ->Fill(chi0);
				hChi[2] ->Fill(chi1);
				hAngularDist[2]->Fill(fabs(TMath::Cos(p4T_ZMF[0].Theta())));
				hAngularDist[2]->Fill(fabs(TMath::Cos(p4T_ZMF[1].Theta())));
			}
			else if(mTTbarParton >= BND[3] )
			{
				hChi[3] ->Fill(chi0);
				hChi[3] ->Fill(chi1);
				hAngularDist[3]->Fill(fabs(TMath::Cos(p4T_ZMF[0].Theta())));
				hAngularDist[3]->Fill(fabs(TMath::Cos(p4T_ZMF[1].Theta())));
			}*/
			
		}
	}
  }	
  
  TCanvas *can[sizeBins];
  for(int i =0; i<sizeBins; i++)
  {
	can[i] = new TCanvas(TString::Format("can %d", (i+1)), TString::Format("can %d", (i+1)), 900, 600);
	if(isParton) hAngularDist[i]->GetXaxis()->SetTitle("cos(#theta)");
	else hAngularDist[i]->GetXaxis()->SetTitle(TString::Format("cos(#theta) %s","reco"));
	hAngularDist[i]->Draw();
  }
  
  TCanvas *canChi[sizeBins];
  for(int i =0; i<sizeBins; i++)
  {
	canChi[i] = new TCanvas(TString::Format("can chi %d", (i+1)), TString::Format("can chi %d", (i+1)), 900, 600);
	if(isParton) hChi[i]->GetXaxis()->SetTitle("#chi");
	else hChi[i]->GetXaxis()->SetTitle(TString::Format("#chi %s", "reco"));
	hChi[i]->Draw();
  }
  
  TCanvas *can_mTTbarParton = new TCanvas("mTTbarParton can", "mTTbarParton can", 900, 600);
  h_mTTbarParton->GetXaxis()->SetTitle("mTTbarParton (GeV)");
  h_mTTbarParton->Draw();
  
  
  TString tempMass;
  if(ZprimeMass == 2000) tempMass = "2TeV";
  else if(ZprimeMass == 3000) tempMass = "3TeV";
  else if(ZprimeMass == 4000) tempMass = "4TeV";
  else if(ZprimeMass == 2500) tempMass = "2.5TeV";
  else if(ZprimeMass == 5000) tempMass = "5TeV";
    
  TFile *outf;
  TString recoParton ="";
  if(isParton) recoParton = "Parton";
  else recoParton = "Reco";
  if(isZprime)  outf = new TFile(TString::Format("Output_M%s_%s_test.root", tempMass.Data(), recoParton.Data()), "UPDATE");
  else outf = new TFile(TString::Format("Output_TT_QCD_%s_test.root", recoParton.Data()), "UPDATE");
  
  if (isZprime) h_mTTbarParton->Write(TString::Format("h_mTTbarParton_M%s_W%s",tempMass.Data(), width.Data() ));
  else h_mTTbarParton->Write("h_mTTbarParton_TT"); 
  
  
  
  for(int i =0; i<sizeBins; i++)
  {
	  TString temp = massLimits[i];
	  if(isZprime) 
	  {
		  hAngularDist[i]->Write(TString::Format("abs_cosTheta_%s_M%s_W%s",temp.Data(), tempMass.Data(), width.Data()));
		  hChi[i]->Write(TString::Format("chi_%s_M%s_W%s",temp.Data(), tempMass.Data(), width.Data()));
	  }
	  else 
	  {
		  hAngularDist[i]->Write(TString::Format("abs_cosTheta_%s_TT",temp.Data() ));
		  hChi[i]->Write(TString::Format("chi_%s_TT", temp.Data()));
	  }
	  
  }	  
  
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

