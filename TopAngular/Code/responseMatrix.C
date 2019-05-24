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


void responseMatrix(TString file = "/eos/cms/store/user/gbakas/ttbar/topTagger/April19/TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8_Copy.root", 
						float selMvaCut=0.3, float floatBTag = 0.8838, bool isZprime=false, int ZprimeMass = 2000, TString width = "200" )
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

  
  
  
  TLorentzVector p4T_parton[2], p4TTbar_parton, p4T_parton_ZMF[2];
  TLorentzVector p4T_reco[2], p4TTbar_reco, p4T_reco_ZMF[2];

  int decade(0);
  int NN = trIN->GetEntries();

  //int NN = 10000;
  const int sizeBins = 5;
  const int chiSize =9;
  //float BND[sizeBins+1] = {1000, 2000, 3000, 4000, 5000};
  //float BND[sizeBins+1] = {1000, 2400, 3000, 3600,4200,4800,5400, 6000, 13000};
  float BND[sizeBins+1] = {1000,1600,2200,3000,3600,6000};
  float BND_chi[chiSize+1] = {1,2,3,4,5,6,8,10,13,16};
  //std::vector<TString> massLimits = {"1000-2400","2400-3000", "3000-3600","3600-4200","4200-4800",
	//								 "4200-4800","4800-5400", "5400-6000","6000-Inf"};
  std::vector<TString> massLimits = {"1000-1600","1600-2200", "2200-3000","3000-3600", "3600-6000"};									 

  int counter =0;
  std::vector< std::vector <Float_t> > const Resp_BND = {{1,2,3,4,5,6,7,8,9,10,12,14,16},
														 {1,2,3,4,5,6,7,8,9,10,12,14,16},
														 {1,2,3,4,5,6,7,8,9,10,12,14,16},
														 {1,2,3,4,5,6,7,8,9,10,12,14,16},
														 {1,2,3,4,5,6,7,8,9,10,12,14,16}};
														 //{1,2,3,4,5,6,7,8,9,10,12,14,16},
														 //{1,2,3,4,5,6,7,8,9,10,12,14,16},
														 //{1,2,3,4,5,6,7,8,9,10,12,14,16}};
  
 int NBINS[sizeBins] = {12,12,12,12,12}; 
 
 
  TH2F *responseMatrix[sizeBins]; 
  TH2F *responseMatrix_mTTTbar = new TH2F("mTTbarParton response", "mTTbarParton response", sizeBins, BND, sizeBins, BND); 
  TH2F *responseMatrix_chi = new TH2F("#chi response", "#chi response", chiSize, BND_chi, chiSize, BND_chi); 
 
  
  TH1F *hDivReco   = new TH1F("hDivReco","hDivReco", chiSize, BND_chi);
  TH1F *hDivParton = new TH1F("hDivParton","hDivParton", chiSize, BND_chi);
 // std::vector<TString> massLimits = {"1000-2500","2500-3500", "3500-5000", "5000-Inf"};

  for(int i=0; i<sizeBins; i++)
  {
	  float tempBND[NBINS[i]+1];
      std::copy(Resp_BND[i].begin(), Resp_BND[i].end(), tempBND);
	  
	  TString temp = massLimits[i];
	  responseMatrix[i] = new TH2F(TString::Format("#chi response matrix for mass limit: %s (GeV)", temp.Data()),TString::Format("#chi response matrix for mass limit: %s (GeV)", temp.Data()),NBINS[i], tempBND,NBINS[i], tempBND);
  }
  //for matching
  std::vector<int> *jetMatchedIndexes = new std::vector<int>(0);
  std::vector<float> *jetMatchedDr = new std::vector<float>(0);
  float jetDr_(0), eta_(0), phi_(0), mass_(0), pt_(0);  
  
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
		bool isMatched = true;
       //----------------------MATCHING------------------------------------------------------
       jetMatchedIndexes->clear();
       jetMatchedDr->clear();
       std::vector<int>::iterator it = std::find(partonMatchIdx->begin(), partonMatchIdx->end(), iev);
       //get all entries that match our jet.
       while(it != partonMatchIdx->end())
       {
           int index = it - partonMatchIdx->begin();
		   jetMatchedIndexes->push_back(index); //has the positions where I found the jet i in partonMatchedIdx
           jetMatchedDr->push_back((*partonMatchDR)[index]); //same here for the DR: DR that correspond to the jet i

           //cout<<"jetFound at: "<<index<<endl;
           ++it;
           it = std::find(it, partonMatchIdx->end(), iev);
	   }
       //if we actually selected something
       if(jetMatchedIndexes->size() > 0)
       {
        
			float dRmin = (*jetMatchedDr)[0];
			int indexMin = (*jetMatchedIndexes)[0];
			
			//cout<<"dRmin[0]: "<<dRmin<<endl;
			for(int k=1; k<jetMatchedIndexes->size(); k++)
			{
				//cout<<"jetMatchedIndexes at k = "<<k<<" is: "<<(*jetMatchedIndexes)[k]<<endl;
				//cout<<"jetMatchedDr at k =  "<<k<<" is: "<<(*jetMatchedDr)[k]<<endl;
				if((*jetMatchedDr)[k] < dRmin)
				{
					dRmin = (*jetMatchedDr)[k];
					indexMin = (*jetMatchedIndexes)[k];
				}
			//cout<<"dRmin is: "<<dRmin<<endl;
			}
			if(dRmin < 0.4)
			{
				isMatched = true;
				jetDr_ = dRmin;
				pt_    = (*jetPt)[(*partonMatchIdx)[indexMin]];
				mass_  = (*jetMassSoftDrop)[(*partonMatchIdx)[indexMin]];
				eta_   = (*jetEta)[(*partonMatchIdx)[indexMin]];
				phi_   = (*jetPhi)[(*partonMatchIdx)[indexMin]];
			}
	    }
		//---------------------------end of MATCHING---------------------------------------------------------
		// Do anything ONLY if matching is ok
		if(isMatched)
		{
			//split the mttbar phase space into regions and for each region I will calculate the thetas dists
			responseMatrix_mTTTbar->Fill(mTTbarParton, mJJ);
			
			p4T_parton[0].SetPtEtaPhiM((*partonPt)[0], (*partonEta)[0], (*partonPhi)[0], (*partonMass)[0]);
			p4T_parton[1].SetPtEtaPhiM((*partonPt)[1], (*partonEta)[1], (*partonPhi)[1], (*partonMass)[1]);
			
			p4T_reco[0].SetPtEtaPhiM((*jetPt)[0], (*jetEta)[0], (*jetPhi)[0], (*jetMassSoftDrop)[0]);
			p4T_reco[1].SetPtEtaPhiM((*jetPt)[1], (*jetEta)[1], (*jetPhi)[1], (*jetMassSoftDrop)[1]);
			
			TVector3 ttbarBoostVector_parton = getBoostVector(p4T_parton[0], p4T_parton[1], p4TTbar_parton);
			
			p4T_parton_ZMF[0].SetPtEtaPhiM(p4T_parton[0].Pt(), p4T_parton[0].Eta(), p4T_parton[0].Phi(), p4T_parton[0].M());
			p4T_parton_ZMF[1].SetPtEtaPhiM(p4T_parton[1].Pt(), p4T_parton[1].Eta(), p4T_parton[1].Phi(), p4T_parton[1].M());
			p4T_parton_ZMF[0].Boost(ttbarBoostVector_parton);
			p4T_parton_ZMF[1].Boost(ttbarBoostVector_parton);
			
			TVector3 ttbarBoostVector_reco = getBoostVector(p4T_reco[0], p4T_reco[1], p4TTbar_reco);
			
			p4T_reco_ZMF[0].SetPtEtaPhiM(p4T_reco[0].Pt(), p4T_reco[0].Eta(), p4T_reco[0].Phi(), p4T_reco[0].M());
			p4T_reco_ZMF[1].SetPtEtaPhiM(p4T_reco[1].Pt(), p4T_reco[1].Eta(), p4T_reco[1].Phi(), p4T_reco[1].M());
			p4T_reco_ZMF[0].Boost(ttbarBoostVector_reco);
			p4T_reco_ZMF[1].Boost(ttbarBoostVector_reco);
					
			//cout<<"-------------------------------------"<<endl;		
			//cout<< p4T_ZMF[0].Pt()<<endl;
			//cout<< p4T_ZMF[1].Pt()<<endl;
			float chi0_parton(0), chi1_parton(0);
			float chi0_reco(0), chi1_reco(0);
			chi0_parton = (1 + fabs(TMath::Cos(p4T_parton_ZMF[0].Theta()))) / ( 1 - fabs(TMath::Cos(p4T_parton_ZMF[0].Theta())));
			chi1_parton = (1 + fabs(TMath::Cos(p4T_parton_ZMF[1].Theta()))) / ( 1 - fabs(TMath::Cos(p4T_parton_ZMF[1].Theta())));
			//chi0 = TMath::Exp(p4T_ZMF[0].Rapidity() - p4T_ZMF[1].Rapidity());
			
			chi0_reco = (1 + fabs(TMath::Cos(p4T_reco_ZMF[0].Theta()))) / ( 1 - fabs(TMath::Cos(p4T_reco_ZMF[0].Theta())));
			chi1_reco = (1 + fabs(TMath::Cos(p4T_reco_ZMF[1].Theta()))) / ( 1 - fabs(TMath::Cos(p4T_reco_ZMF[1].Theta())));
			
			/*
			cout<<"-------------------------------------------------------------"<<endl;
			cout<<chi0_parton<<endl;
			cout<<chi0_reco<<endl;
			cout<<chi1_parton<<endl;
			cout<<chi1_reco<<endl;
			*/
			if(recoCuts && partonCuts && topTagger && btagging)
			{
				responseMatrix_chi->Fill(chi0_parton, chi0_reco);
				responseMatrix_chi->Fill(chi1_parton, chi1_reco);
				bool found =false;
				int imass = 0;
				while(!found && imass<sizeBins)
				{
					//cout<<imass<<endl;
					//cout<<"------------------------"<<endl;
					if(mTTbarParton >= BND[imass] && mTTbarParton < BND[imass+1] )
					{
						found =true;
						responseMatrix[imass]->Fill(chi0_parton, chi0_reco);
						responseMatrix[imass]->Fill(chi1_parton, chi1_reco);
					}
					imass++;
				}
				
			}
			if(recoCuts && topTagger && btagging)
			{
			  hDivReco->Fill(chi0_reco);
			  hDivReco->Fill(chi1_reco);
			}
			if(partonCuts)
			{
			  hDivParton->Fill(chi0_parton);
			  hDivParton->Fill(chi1_parton);
			}
		
	  }//-----end of isMatched loop 
	  
	}//---end of nJets >1 loop
	
 }//----end of iev loop	
  
  TFile *outf = new TFile("resp.root", "RECREATE");
  
  TCanvas *can_resp[sizeBins];	
  for(int imass=0; imass<sizeBins; imass++)
  {
	can_resp[imass] = new TCanvas(TString::Format("responseMatrix %d",imass+1),TString::Format("responseMatrix %d", imass+1), 900, 600); 
	can_resp[imass]->cd();
	responseMatrix[imass]->Draw("textcolz");
	responseMatrix[imass]->Write();
  }
  
  
  TCanvas *can_mTTbarParton = new TCanvas("responseMatrix_chi can", "responseMatrix_chi can", 900, 600);
  responseMatrix_chi->GetXaxis()->SetTitle("#chi parton");
  responseMatrix_chi->GetYaxis()->SetTitle("#chi reco");  
  responseMatrix_chi->Draw("textColz");
  responseMatrix_chi->Write();
  
  
  
  
  //this is now for purity and stability for chi
  //first find sums
  float bins[chiSize+1];
  for(int i=0; i<chiSize+1; i++)
  {
	  bins[i]=i+1;
  }
  TH1F *purity = new TH1F ("purity", "purity", chiSize, bins);
  TH1F *stability = new TH1F ("stability", "stability", chiSize, bins);
  float sumOfRows[chiSize], sumOfCols[chiSize];
  for(int i=0; i<chiSize+1; i++)
  {
	 sumOfCols[i] = ((TH1D*)responseMatrix_chi->ProjectionX())->GetBinContent(i);
	 sumOfRows[i] = ((TH1D*)responseMatrix_chi->ProjectionY())->GetBinContent(i);
	 
	 for(int j=0; j<chiSize+1; j++)
	 {
		if(i==j)
		{
			float initContent = responseMatrix_chi->GetBinContent(i,j);
			purity->SetBinContent(i,initContent/sumOfCols[i]);
			stability->SetBinContent(i,initContent/sumOfRows[i]);
		}		
	 }
  }
  //purity: sum all over the columns and find binContent(i,j)/SumOfColumn(j) for all i 
  //stability: sum all over the lines and find binContent/SumOfLine(i) for all jetBtagSub0
  TCanvas *can_pur = new TCanvas("test", "test", 900, 600);
  TLegend *leg_purityStability = new TLegend(0.5,0.6,0.7,0.9);
  leg_purityStability->AddEntry(purity, "Purity", "l");
  leg_purityStability->AddEntry(stability, "Stability", "l");
  purity->GetXaxis()->SetTitle("Bin Number");
  stability->GetXaxis()->SetTitle("Bin Number");
  purity->SetLineColor(kRed);
  stability->SetLineColor(kBlue);
  purity->Draw();
  stability->Draw("same");
  leg_purityStability->Draw();
  
  
  TH1D *efficiency, *acceptance;
  
  efficiency = responseMatrix_chi->ProjectionX();
  acceptance = responseMatrix_chi->ProjectionY();
  
  acceptance->Divide(hDivReco);
  efficiency->Divide(hDivParton);
  TCanvas *can_test = new TCanvas("acceptance, efficiency can", "acceptance, efficiency can", 900, 600);
  TLegend *leg_accEff = new TLegend(0.5,0.6,0.7,0.9);
  efficiency->GetXaxis()->SetTitle("#chi");
  acceptance->GetXaxis()->SetTitle("#chi");
  efficiency->SetLineColor(kRed);
  acceptance->SetLineColor(kBlue);
  efficiency->Draw();
  acceptance->Draw("same");
  leg_accEff->AddEntry(acceptance, "Acceptance", "l");
  leg_accEff->AddEntry(efficiency, "Efficiency", "l");
  leg_accEff->Draw();

  
  
  acceptance->Write("chiAcceptance");
  efficiency->Write("chiEfficiency");
  stability->Write("chiStability");
  purity->Write("chiPurity");
  
  //outf->Close();
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

