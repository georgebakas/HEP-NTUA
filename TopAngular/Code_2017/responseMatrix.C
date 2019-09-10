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


void responseMatrix(TString file = "/eos/cms/store/user/gbakas/ttbar/topTagger_old/mc-2017/Signal/TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root", 
						float selMvaCut=0.1, float floatBTag = 0.8838, bool isDeepCSV = true, bool isZprime=false, int ZprimeMass = 2000, TString width = "200" )
{
	
//TString TTbarFile = "/eos/cms/store/user/gbakas/ttbar/topTagger/April19/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_Copy.root"	
  gStyle->SetOptStat(0);
  TFile *inf     = TFile::Open(file);
  TTree *trIN    = (TTree*)inf->Get("boosted/events");
  //cout<<"here"<<endl;
  float deepCSVFloat = 0.4941;
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
  vector<bool> *bit = new vector<bool>;
  float mTTbarParton(0),mJJ(0), yTTbarParton(0), ptTTbarParton(0);
  int  category(0);
  //matching info
  vector<float> *jetPhi(0), *jetEta(0);
  vector<int> *partonId(0), *partonMatchIdx(0);
  
  vector<float> *partonEta(0), *partonPhi(0), *partonMatchDR(0),  *partonPt(0), *partonE(0), *partonMass(0);
  std::vector<int> *addedIndexes = new std::vector<int>(0);
  std::vector<float> *jetBtagSub0(0), *jetBtagSub1(0);
  std::vector<float> *jetBtagSub0DCSVbb(0), *jetBtagSub1DCSVbb(0);
  std::vector<float> *jetBtagSub0DCSVbbb(0), *jetBtagSub1DCSVbbb(0);      
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
   //deepCSV
  trIN->SetBranchAddress("jetBtagSub0DCSVbb" ,&jetBtagSub0DCSVbb);
  trIN->SetBranchAddress("jetBtagSub1DCSVbb" ,&jetBtagSub1DCSVbb);
  trIN->SetBranchAddress("jetBtagSub0DCSVbbb",&jetBtagSub0DCSVbbb);
  trIN->SetBranchAddress("jetBtagSub1DCSVbbb",&jetBtagSub1DCSVbbb);
  
  
  TLorentzVector p4T_parton[2], p4TTbar_parton, p4T_parton_ZMF[2];
  TLorentzVector p4T_reco[2], p4TTbar_reco, p4T_reco_ZMF[2];

  int decade(0); 
  int NN = trIN->GetEntries();

  //int NN = 10000;
  const int sizeBins = 4;
  const int chiSize =11;
  //float BND[sizeBins+1] = {1000, 2000, 3000, 4000, 5000};
  //float BND[sizeBins+1] = {1000, 2400, 3000, 3600,4200,4800,5400, 6000, 13000};
  float BND[sizeBins+1] = {1000,1600,2200,3200,6000};
  float BND_chi[chiSize+1] = {1,2,3,4,5,6,7,8,9,10,13,16};
  //std::vector<TString> massLimits = {"1000-2400","2400-3000", "3000-3600","3600-4200","4200-4800",
	//								 "4200-4800","4800-5400", "5400-6000","6000-Inf"};
  std::vector<TString> massLimits = {"1000-1600","1600-2200", "2200-3200","3200-6000"};									 

  int counter =0;
  std::vector< std::vector <Float_t> > const Resp_BND = {{1,2,3,4,5,6,7,8,9,10,13,16},
														 {1,2,3,4,5,6,7,8,9,10,13,16},
														 {1,2,3,4,5,6,7,8,9,10,13,16},
														 {1,2,3,4,5,6,7,8,9,10,13,16}};
														 //{1,2,3,4,5,6,8,10,13,16}};
														 //{1,2,3,4,5,6,7,8,9,10,12,14,16},
														 //{1,2,3,4,5,6,7,8,9,10,12,14,16},
														 //{1,2,3,4,5,6,7,8,9,10,12,14,16}};
  
 int NBINS[sizeBins] = {chiSize,chiSize,chiSize,chiSize}; 
 const int cosSize = 10;
 //float BND_cos[cosSize+1] = {0,0.1,0.2,0.3,0.4,0.6,0.8,1};
 
  TH2F *responseMatrix[sizeBins]; 
  TH2F *responseMatrix_mTTTbar = new TH2F("mTTbarParton response", "mTTbarParton response", sizeBins, BND, sizeBins, BND); 
  TH2F *responseMatrix_chi = new TH2F("#chi response", "#chi response", chiSize, BND_chi, chiSize, BND_chi); 
  TH2F *responseMatrix_cos = new TH2F("|cos(#theta)| response", "|cos(#theta)| response",cosSize ,0, 1, cosSize, 0,1); 
 
  TH1F *hDivReco      = new TH1F("hDivReco","hDivReco", chiSize, BND_chi);
  TH1F *hDivRecoCos   = new TH1F("hDivRecoCos","hDivRecoCos", cosSize, 0,1);
  TH1F *hDivParton    = new TH1F("hDivParton","hDivParton", chiSize, BND_chi);
  TH1F *hDivPartonCos = new TH1F("hDivPartonCos","hDivPartonCos", cosSize, 0,1);
 // std::vector<TString> massLimits = {"1000-2500","2500-3500", "3500-5000", "5000-Inf"};
  TH1F *efficiency    = new TH1F("efficiency #chi", "efficiency #chi", chiSize, BND_chi); 
  TH1F *efficiencyCos = new TH1F("efficiency cos(#theta)", "efficiency cos(#theta)", cosSize, 0,1); 
  TH1F *acceptance    = new TH1F("acceptance #chi", "acceptance #chi", chiSize, BND_chi); 
  TH1F *acceptanceCos = new TH1F("acceptance cos(#theta)", "acceptance cos(#theta)", cosSize, 0,1); 
 
 //for testing
 TH1F *hChiPartonTest[3], *hChiRecoTest[3];
 hChiPartonTest[0]   = new TH1F("hChiPartonTest with exp ZMF","hChiPartonTest with exp ZMF", chiSize, BND_chi);
 hChiPartonTest[1]   = new TH1F("hChiPartonTest with cos","hChiPartonTest with cos", chiSize, BND_chi);
 hChiPartonTest[2]   = new TH1F("hChiPartonTest with #DeltaY LB","hChiPartonTest #DeltaY LB", chiSize, BND_chi);
 
 hChiRecoTest[0]   = new TH1F("hChiRecoTest with exp","hChiRecoTest with exp", chiSize, BND_chi);
 hChiRecoTest[1]   = new TH1F("hChiRecoTest with cos","hChiRecoTest with cos", chiSize, BND_chi);
 hChiRecoTest[2]   = new TH1F("hChiRecoTest with #DeltaY LB","hChiRecoTest #DeltaY LB", chiSize, BND_chi);


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
  std::vector<float> *eta_ = new std::vector<float>(0);
  std::vector<float> *phi_ = new std::vector<float>(0);
  std::vector<float> *mass_ = new std::vector<float>(0);
  std::vector<float> *pt_ = new std::vector<float>(0);
  std::vector<float> *jetBtagSub0_ = new std::vector<float>(0);
  std::vector<float> *jetBtagSub1_ = new std::vector<float>(0);
  std::vector<float> *jetTtag_ = new std::vector<float>(0);
  
  std::vector<float> *partonPt_ = new std::vector<float>(0);
  std::vector<float> *partonEta_ = new std::vector<float>(0);
  std::vector<float> *partonMass_ = new std::vector<float>(0);
  std::vector<float> *partonPhi_ = new std::vector<float>(0);
  
  std::vector<float> *jetBtagSub0DCSVbb_ = new std::vector<float>(0);
  std::vector<float> *jetBtagSub1DCSVbb_ = new std::vector<float>(0);
  std::vector<float> *jetBtagSub0DCSVbbb_ = new std::vector<float>(0);
  std::vector<float> *jetBtagSub1DCSVbbb_ = new std::vector<float>(0);
  
  float jetDr_(0);// eta_(0), phi_(0), mass_(0), pt_(0);  
  
  cout<<"Reading "<<NN<<" entries"<<endl;
  for(int iev=0;iev<NN;iev++) 
  {
    double progress = 10.0*iev/(1.0*NN);
    int k = TMath::FloorNint(progress); 
    if (k > decade) 
      cout<<10*k<<" %"<<endl;
    decade = k;
    trIN->GetEntry(iev);
	int isMatched = 0;
	bool recoCuts, partonCuts, btagging, topTagger, massCut;	
	eta_->clear();
	mass_->clear();
	pt_->clear();
	phi_->clear();
	jetBtagSub0_->clear();
	jetBtagSub1_->clear();
	jetTtag_->clear();
	
	
	partonPt_->clear();
	partonMass_->clear();
	partonPhi_->clear();
	partonEta_->clear();
	
	jetBtagSub0DCSVbb_->clear();
    jetBtagSub1DCSVbb_->clear();
    jetBtagSub0DCSVbbb_->clear();
    jetBtagSub1DCSVbbb_->clear();
	
	if (nJets >1)
	{	

		//----------------------MATCHING------------------------------------------------------
		
		for(int ijet =0; ijet<nJets; ijet++)
		{
			//cout<<"ok"<<endl;
		   jetMatchedIndexes->clear();
		   jetMatchedDr->clear();
		   std::vector<int>::iterator it = std::find(partonMatchIdx->begin(), partonMatchIdx->end(), ijet);
		   //get all entries that match our jet.
		   while(it != partonMatchIdx->end())
		   {
			   int index = it - partonMatchIdx->begin();
			   jetMatchedIndexes->push_back(index); //has the positions where I found the jet i in partonMatchedIdx
			   jetMatchedDr->push_back((*partonMatchDR)[index]); //same here for the DR: DR that correspond to the jet i

			   //cout<<"jetFound at: "<<index<<endl;
			   ++it;
			   it = std::find(it, partonMatchIdx->end(), ijet);
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
					isMatched++;
					//cout<<"int isMatched: "<<isMatched<<endl;
					jetDr_ = dRmin;
					pt_->push_back((*jetPt)[(*partonMatchIdx)[indexMin]]);
					mass_->push_back((*jetMassSoftDrop)[(*partonMatchIdx)[indexMin]]);
					eta_->push_back((*jetEta)[(*partonMatchIdx)[indexMin]]);
					phi_->push_back( (*jetPhi)[(*partonMatchIdx)[indexMin]]);
					jetBtagSub0_->push_back( (*jetBtagSub0)[(*partonMatchIdx)[indexMin]]);
					jetBtagSub1_->push_back( (*jetBtagSub1)[(*partonMatchIdx)[indexMin]]);
					jetTtag_->push_back( (*jetTtag)[(*partonMatchIdx)[indexMin]]);
					
					jetBtagSub0DCSVbb_->push_back((*jetBtagSub0DCSVbb)[(*partonMatchIdx)[indexMin]]);
					jetBtagSub1DCSVbb_->push_back((*jetBtagSub1DCSVbb)[(*partonMatchIdx)[indexMin]]);
					jetBtagSub0DCSVbbb_->push_back((*jetBtagSub0DCSVbbb)[(*partonMatchIdx)[indexMin]]);
					jetBtagSub1DCSVbbb_->push_back((*jetBtagSub1DCSVbbb)[(*partonMatchIdx)[indexMin]]);
					
					partonPt_->push_back((*partonPt)[indexMin]);
					partonMass_->push_back((*partonMass)[indexMin]);
					partonPhi_->push_back((*partonPhi)[indexMin]);
					partonEta_->push_back((*partonEta)[indexMin]);
				}
					
		   }
		   
		 }//---------------------------end of MATCHING---------------------------------------------------------
	     if(isMatched >1)
		 {
			    float dCSVScoreSub0[2], dCSVScoreSub1[2];
				dCSVScoreSub0[0] = (*jetBtagSub0DCSVbb_)[0] + (*jetBtagSub0DCSVbbb_)[0];
				dCSVScoreSub0[1] = (*jetBtagSub0DCSVbb_)[1] + (*jetBtagSub0DCSVbbb_)[1];
				dCSVScoreSub1[0] = (*jetBtagSub1DCSVbb_)[0] + (*jetBtagSub1DCSVbbb_)[0];
				dCSVScoreSub1[1] = (*jetBtagSub1DCSVbb_)[1] + (*jetBtagSub1DCSVbbb_)[1];
				bool CSVv2Cut, deepCSVCut;
				//cout<<"------------------------------------------------------"<<endl;
				//cout<<"pt[0]: "<<(*pt_)[0]<<endl;
				//cout<<"pt[1]: "<<(*pt_)[1]<<endl;
				// Do anything ONLY if matching is ok
				recoCuts   = fabs((*eta_)[0]) < 2.4 && fabs((*eta_)[1]) <2.4 && (*pt_)[0] > 400 && (*pt_)[1] > 400 && nLeptons==0 && mJJ > 1000 && (*bit)[5];
				partonCuts = (*partonPt_)[0] > 400 && (*partonPt_)[1] > 400 && fabs((*partonEta_)[0]) < 2.4 && fabs((*partonEta_)[1]) < 2.4 &&  mTTbarParton > 1000 && (*bit)[5];
				CSVv2Cut   = ((*jetBtagSub0_)[0] > floatBTag || (*jetBtagSub1_)[0] > floatBTag) && ((*jetBtagSub0_)[1] > floatBTag || (*jetBtagSub1_)[1] > floatBTag);
				deepCSVCut = (dCSVScoreSub0[0] > deepCSVFloat || dCSVScoreSub1[0] > deepCSVFloat) && (dCSVScoreSub0[1] > deepCSVFloat || dCSVScoreSub1[1] > deepCSVFloat);
				topTagger  = (*jetTtag_)[0] > selMvaCut && (*jetTtag_)[1] > selMvaCut;
				massCut    = (*mass_)[0] > 120 && (*mass_)[0] < 220 && (*mass_)[1] > 120 && (*mass_)[1] < 220;
			
				if(isDeepCSV) btagging = deepCSVCut;
				else btagging = CSVv2Cut;
				
				//split the mttbar phase space into regions and for each region I will calculate the thetas dists
				responseMatrix_mTTTbar->Fill(mTTbarParton, mJJ);
			
				p4T_parton[0].SetPtEtaPhiM((*partonPt_)[0], (*partonEta_)[0], (*partonPhi_)[0], (*partonMass_)[0]);
				p4T_parton[1].SetPtEtaPhiM((*partonPt_)[1], (*partonEta_)[1], (*partonPhi_)[1], (*partonMass_)[1]);
			
				p4T_reco[0].SetPtEtaPhiM((*pt_)[0], (*eta_)[0], (*phi_)[0], (*mass_)[0]);
				p4T_reco[1].SetPtEtaPhiM((*pt_)[1], (*eta_)[1], (*phi_)[1], (*mass_)[1]);
				
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
						
				
				float chi_parton(0);
				float chi_reco(0);
				
				float chi_parton_TEST= (1 + fabs(TMath::Cos(p4T_parton_ZMF[0].Theta()))) / ( 1 - fabs(TMath::Cos(p4T_parton_ZMF[0].Theta())));
				float chi_reco_TEST  = (1 + fabs(TMath::Cos(p4T_reco_ZMF[0].Theta()))) / ( 1 - fabs(TMath::Cos(p4T_reco_ZMF[0].Theta())));
				
				chi_parton = TMath::Exp(fabs(p4T_parton_ZMF[0].Rapidity() - p4T_parton_ZMF[1].Rapidity()));
				chi_reco   = TMath::Exp(fabs(p4T_reco_ZMF[0].Rapidity() - p4T_reco_ZMF[1].Rapidity()));
				
				float chi_partonLF = TMath::Exp(fabs(p4T_parton[0].Rapidity() - p4T_parton[1].Rapidity()));
				float chi_recoLF   = TMath::Exp(fabs(p4T_reco[0].Rapidity() - p4T_reco[1].Rapidity()));
				
				//float chi0_partonT = TMath::Exp(fabs(p4T_parton[0].Rapidity() - p4T_parton[1].Rapidity()));
				//float chi0_recoT   = TMath::Exp(fabs(p4T_reco[0].Rapidity() - p4T_reco[1].Rapidity()));
				//cout<<TMath::Exp(fabs(p4T_reco_ZMF[0].Rapidity() - p4T_reco_ZMF[1].Rapidity()))<<endl;
				hChiPartonTest[0]->Fill(chi_parton,genEvtWeight);
				hChiPartonTest[1]->Fill(chi_parton_TEST,genEvtWeight);
				hChiPartonTest[2]->Fill(chi_partonLF,genEvtWeight);
				
				hChiRecoTest[0]->Fill(chi_reco, genEvtWeight);
				hChiRecoTest[1]->Fill(chi_reco_TEST, genEvtWeight);
				hChiRecoTest[2]->Fill(chi_recoLF, genEvtWeight);
				
				
				//cout<<"-------------------------------------------------------------"<<endl;
				//cout<<chi_parton<<endl;
				//cout<<chi0_partonT<<endl;
				
				if(recoCuts && partonCuts && topTagger && btagging && massCut)
				{
					responseMatrix_chi->Fill(chi_parton, chi_reco);
					responseMatrix_cos->Fill(fabs(TMath::Cos(p4T_parton_ZMF[0].Theta())), fabs(TMath::Cos(p4T_reco_ZMF[0].Theta())));
					
					efficiency->Fill(chi_parton,genEvtWeight);
					acceptance->Fill(chi_reco,genEvtWeight);
					
					efficiencyCos->Fill(fabs(TMath::Cos(p4T_parton_ZMF[0].Theta())),genEvtWeight);
					acceptanceCos->Fill(fabs(TMath::Cos(p4T_reco_ZMF[0].Theta())),genEvtWeight);
					
					//responseMatrix_chi->Fill(chi1_parton, chi1_reco);
					bool found =false;
					int imass = 0;
					while(!found && imass<sizeBins)
					{
						//cout<<imass<<endl;
						//cout<<"------------------------"<<endl;
						if(mTTbarParton >= BND[imass] && mTTbarParton < BND[imass+1] )
						{
							found =true;
							responseMatrix[imass]->Fill(chi_parton, chi_reco);
							//responseMatrix[imass]->Fill(chi1_parton, chi1_reco);
						}
						imass++;
					}
					
				}
				if(recoCuts && topTagger && btagging && massCut)
				{
				  hDivReco->Fill(chi_reco,genEvtWeight);
				  hDivRecoCos->Fill(fabs(TMath::Cos(p4T_reco_ZMF[0].Theta())),genEvtWeight);
				}

				
		 }//----end of isMatched
	  
	}//---end of nJets >1 loop
	
 }//----end of iev loop	
 
 hChiPartonTest[0]->GetXaxis()->SetTitle("#chi");
 hChiPartonTest[1]->GetXaxis()->SetTitle("#chi");
 hChiPartonTest[2]->GetXaxis()->SetTitle("#chi");
 TCanvas *c1 = new TCanvas ("c1", "c1", 700, 600);
 TLegend *leg1 = new TLegend(0.5,0.6,0.7,0.8);
 leg1->AddEntry(hChiPartonTest[0],"Chi parton #DeltaY ZMF", "l");
 leg1->AddEntry(hChiPartonTest[1],"Chi parton with Cos", "l");
 leg1->AddEntry(hChiPartonTest[2],"Chi parton #DeltaY Lab Frame", "l"); 
 hChiPartonTest[0]->SetLineColor(kBlue);
 hChiPartonTest[0]->Draw();
 hChiPartonTest[1]->SetLineColor(kRed);
 hChiPartonTest[1]->Draw("same");
 hChiPartonTest[2]->SetLineColor(kMagenta);
 hChiPartonTest[2]->Draw("same");
 leg1->Draw();
 
 
 TCanvas *c2 = new TCanvas ("c2", "c2", 700, 600);
 TLegend *leg2 = new TLegend(0.5,0.6,0.7,0.8);
 leg2->AddEntry(hChiRecoTest[0],"Chi reco Exp #DeltaY ZMF", "l");
 leg2->AddEntry(hChiRecoTest[1],"Chi reco with Cos", "l");
 leg2->AddEntry(hChiRecoTest[2],"Chi reco #DeltaY Lab Frame", "l"); 

 hChiRecoTest[0]->GetXaxis()->SetTitle("#chi");
 hChiRecoTest[1]->GetXaxis()->SetTitle("#chi");
 hChiRecoTest[2]->GetXaxis()->SetTitle("#chi");
 hChiRecoTest[0]->SetLineColor(kBlue);
 hChiRecoTest[0]->Draw();
 hChiRecoTest[1]->SetLineColor(kRed);
 hChiRecoTest[1]->Draw("same");
 hChiRecoTest[2]->SetLineColor(kMagenta);
 hChiRecoTest[2]->Draw("same");
 leg2->Draw();
  //--------------------------------------------START OF EVENT COUNTER LOOP -------------------------------------------------------------------

  //now another for that fills the denominators for the parton efficiencies
  //loop over other tree -> eventCounter
  TFile *infCnt = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Signal/2017/TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8_EventCounter.root");
  TTree *trCnt = (TTree*)infCnt->Get("eventCounter/events");
  float ptTTbarPartonCnt(0), mTTbarPartonCnt(0), yTTbarPartonCnt(0), genEvtWeightCnt(0);
  float partonPtCnt[2], partonEtaCnt[2],partonYCnt[2], partonMassCnt[2], phiTopParton[2];
  //tree for eventCounter
  trCnt->SetBranchAddress("ptTopParton"    ,&partonPtCnt);
  trCnt->SetBranchAddress("etaTopParton"   ,&partonEtaCnt);
  trCnt->SetBranchAddress("yTopParton"     ,&partonYCnt);
  trCnt->SetBranchAddress("mTopParton"     ,&partonMassCnt);
  trCnt->SetBranchAddress("phiTopParton"   ,&phiTopParton);
  trCnt->SetBranchAddress("mTTbarParton"   ,&mTTbarPartonCnt);
  trCnt->SetBranchAddress("yTTbarParton"   ,&yTTbarPartonCnt);
  trCnt->SetBranchAddress("ptTTbarParton"  ,&ptTTbarPartonCnt);
  trCnt->SetBranchAddress("genEvtWeight"   ,&genEvtWeightCnt);

  int NNCnt = trCnt->GetEntries();
  TH1F *hChiEventCounter[3];
  hChiEventCounter[0]   = new TH1F("hChiEventCounter exp ZMF","hChiEventCounter exp ZMF", chiSize, BND_chi);
  hChiEventCounter[1]   = new TH1F("hChiEventCounter with COS","hChiEventCounter with cos", chiSize, BND_chi);
  hChiEventCounter[2]   = new TH1F("hChiEventCounter exp y1-y2 lab frame","hChiEventCounter exp y1-y2 lab frame", chiSize, BND_chi);
  
  TLorentzVector p4T_EventCounter[2],p4T_EventCounter_ZMF[2] ,p4TTbar_EventCounter;
  for(int iev = 0; iev < NNCnt; iev++)
  {
	trCnt->GetEntry(iev);
	
    bool partonCuts = fabs(partonEtaCnt[0]) < 2.4 && fabs(partonEtaCnt[1]) <2.4 && partonPtCnt[0] > 400 && partonPtCnt[1] > 400 && mTTbarPartonCnt > 1000;
	if(partonCuts)
	{
	  p4T_EventCounter[0].SetPtEtaPhiM(partonPtCnt[0], partonEtaCnt[0], phiTopParton[0], partonMassCnt[0]);	
	  p4T_EventCounter[1].SetPtEtaPhiM(partonPtCnt[1], partonEtaCnt[1], phiTopParton[1], partonMassCnt[1]);	
	
	  TVector3 ttbarBoostVector_partonEvent = getBoostVector(p4T_EventCounter[0], p4T_EventCounter[1], p4TTbar_EventCounter);	
	  p4T_EventCounter_ZMF[0].SetPtEtaPhiM(p4T_EventCounter[0].Pt(), p4T_EventCounter[0].Eta(), p4T_EventCounter[0].Phi(), p4T_EventCounter[0].M());
	  p4T_EventCounter_ZMF[1].SetPtEtaPhiM(p4T_EventCounter[1].Pt(), p4T_EventCounter[1].Eta(), p4T_EventCounter[1].Phi(), p4T_EventCounter[1].M());
	  p4T_EventCounter_ZMF[0].Boost(ttbarBoostVector_partonEvent);
      p4T_EventCounter_ZMF[1].Boost(ttbarBoostVector_partonEvent);
   	
      float chi_parton = TMath::Exp(fabs(p4T_EventCounter_ZMF[0].Rapidity() - p4T_EventCounter_ZMF[1].Rapidity()));	  
      float chi_partonEvent = TMath::Exp(fabs(partonYCnt[0] - partonYCnt[1]));	  
      float chi_partonEventTEST = TMath::Exp(fabs(p4T_EventCounter[0].Rapidity() - p4T_EventCounter[1].Rapidity()));	  
	  float chi_partonCos = (1 + fabs(TMath::Cos(p4T_EventCounter_ZMF[0].Theta()))) / ( 1 - fabs(TMath::Cos(p4T_EventCounter_ZMF[0].Theta())));
	    
	  	/*
	  cout<<"--------------Combined--------------------"<<endl;
	  cout<<p4TTbar_EventCounter.Pt()<<endl;
	  cout<<p4TTbar_EventCounter.M()<<endl;
	  cout<<p4TTbar_EventCounter.Eta()<<endl;
	  cout<<p4TTbar_EventCounter.Rapidity()<<endl;
	  cout<<"--------------Y--------------------"<<endl;
	  cout<<chi_parton<<endl;
	  cout<<chi_partonEvent<<endl;
	  cout<<chi_partonEventTEST<<endl;
	  cout<<chi_partonCos<<endl;
		*/
	  hChiEventCounter[0]->Fill(chi_parton,genEvtWeightCnt);
	  hChiEventCounter[1]->Fill(chi_partonCos,genEvtWeightCnt);
	  hChiEventCounter[2]->Fill(chi_partonEvent,genEvtWeightCnt);
	  hDivParton->Fill(chi_parton,genEvtWeightCnt);
	  hDivPartonCos->Fill(fabs(TMath::Cos(p4T_EventCounter_ZMF[0].Theta())),genEvtWeightCnt);
	}
         
  }
  cout<<hDivParton->GetBinContent(1)<<endl;
  //--------------------------------------------END OF EVENT COUNTER LOOP -------------------------------------------------------------------
  TCanvas *c3 = new TCanvas("c3", "c3", 700, 600);
  TLegend *leg3 = new TLegend(0.5,0.6,0.7,0.8);
  leg3->AddEntry(hChiEventCounter[0],"Chi EvntCnt #DeltaY ZMF", "l");
  leg3->AddEntry(hChiEventCounter[1],"Chi EvntCnt with Cos", "l"); 
  leg3->AddEntry(hChiEventCounter[2],"Chi EvntCnt #DeltaY Lab Frame", "l"); 
  hChiEventCounter[0]->SetLineColor(kBlue);
  hChiEventCounter[0]->GetXaxis()->SetTitle("#chi");
  hChiEventCounter[1]->GetXaxis()->SetTitle("#chi");
  hChiEventCounter[2]->GetXaxis()->SetTitle("#chi");
  hChiEventCounter[0]->Draw();
  hChiEventCounter[1]->SetLineColor(kRed);
  hChiEventCounter[1]->Draw("same");
  hChiEventCounter[2]->SetLineColor(kMagenta);
  hChiEventCounter[2]->Draw("same");
  leg3->Draw();
  
  TFile *outf;
  if(isDeepCSV) outf = new TFile(TString::Format("ResponseMatricesChiCos_%0.1f_deepCSV.root",selMvaCut), "RECREATE");
  else outf = new TFile(TString::Format("ResponseMatricesChiCos_%0.1f.root",selMvaCut), "RECREATE");
  outf->cd();
  hDivParton->Write("testParton");
  TCanvas *can_resp[sizeBins];	
  for(int imass=0; imass<sizeBins; imass++)
  {
	can_resp[imass] = new TCanvas(TString::Format("responseMatrix_%d",imass+1),TString::Format("responseMatrix_%d", imass+1), 900, 600); 
	can_resp[imass]->cd();
	responseMatrix[imass]->Draw("textcolz");
	responseMatrix[imass]->Write();
	can_resp[imass]->Write();
  }
  
  
  TCanvas *can_chiResp = new TCanvas("responseMatrix_chi_can", "responseMatrix_chi_can", 900, 600);
  responseMatrix_chi->GetXaxis()->SetTitle("#chi parton");
  responseMatrix_chi->GetYaxis()->SetTitle("#chi reco");  
  responseMatrix_chi->Draw("textColz");
  responseMatrix_chi->Write("chiResponse");
  can_chiResp->Write();
  
  TCanvas *can_cosThetaResp = new TCanvas("responseMatrix_cosTheta_can", "responseMatrix_cosTheta_can", 900, 600);
  responseMatrix_cos->GetXaxis()->SetTitle("|cos(#theta)| parton");
  responseMatrix_cos->GetYaxis()->SetTitle("|cos(#theta)| reco");  
  responseMatrix_cos->Draw("textColz");
  responseMatrix_cos->Write("cosResponse");
  can_cosThetaResp->Write();
  
  
  
  
  //this is now for purity and stability for chi
  //first find sums
  float bins[chiSize+1];
  for(int i=0; i<chiSize+1; i++)
  {
	  bins[i]=i+1;
  }
  TH1F *purity = new TH1F ("purity_Chi", "purity_Chi", chiSize, bins);
  TH1F *stability = new TH1F ("stability_Chi", "stability_Chi", chiSize, bins);
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
  
  //this is now for purity and stability for cos(theta)
  //first find sums
  float binsCos[cosSize+1];
  for(int i=0; i<cosSize; i++)
  {
	  bins[i]=i+1;
  }
  TH1F *purityCos = new TH1F ("purity_CosTheta", "purity_CosTheta", cosSize, bins);
  TH1F *stabilityCos = new TH1F ("stability_CosTheta", "stability_CosTheta", cosSize, bins);
  float sumOfRowsCos[cosSize], sumOfColsCos[cosSize];
  for(int i=0; i<cosSize+1; i++)
  {
	 sumOfColsCos[i] = ((TH1D*)responseMatrix_cos->ProjectionX())->GetBinContent(i);
	 sumOfRowsCos[i] = ((TH1D*)responseMatrix_cos->ProjectionY())->GetBinContent(i);
	 
	 for(int j=0; j<chiSize+1; j++)
	 {
		if(i==j)
		{
			float initContent = responseMatrix_cos->GetBinContent(i,j);
			purityCos->SetBinContent(i,initContent/sumOfColsCos[i]);
			stabilityCos->SetBinContent(i,initContent/sumOfRowsCos[i]);
		}		
	 }
  }
  //purity: sum all over the columns and find binContent(i,j)/SumOfColumn(j) for all i 
  //stability: sum all over the lines and find binContent/SumOfLine(i) for all jetBtagSub0
  TCanvas *can_pur = new TCanvas("purityStability_Chi", "purityStability_Chi", 900, 600);
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
  
  
  TCanvas *can_purCos = new TCanvas("purityStability_Cos", "purityStability_Cos", 900, 600);
  TLegend *leg_purityStabilityCos = new TLegend(0.5,0.6,0.7,0.9);
  leg_purityStabilityCos->AddEntry(purityCos, "Purity", "l");
  leg_purityStabilityCos->AddEntry(stabilityCos, "Stability", "l");
  purityCos->GetXaxis()->SetTitle("Bin Number");
  stabilityCos->GetXaxis()->SetTitle("Bin Number");
  purityCos->SetLineColor(kRed);
  stabilityCos->SetLineColor(kBlue);
  purityCos->Draw();
  stabilityCos->Draw("same");
  leg_purityStabilityCos->Draw();
  
  
  /* 
  efficiency = responseMatrix_chi->ProjectionX();
  acceptance = responseMatrix_chi->ProjectionY();
  
  efficiencyCos = responseMatrix_cos->ProjectionX();
  acceptanceCos = responseMatrix_cos->ProjectionY();
	*/
  TEfficiency *acc = new TEfficiency(*acceptance,*hDivReco);
  TEfficiency *eff = new TEfficiency(*efficiency,*hDivParton);
  TCanvas *can_effAcc = new TCanvas("acceptance, efficiency can", "acceptance, efficiency can", 900, 600);
  TLegend *leg_accEff = new TLegend(0.5,0.6,0.7,0.9);
  eff->SetTitle("Efficiency;#chi;Efficiency");
  acc->SetTitle("Acceptance;#chi;Acceptance");
  eff->SetLineColor(kRed);
  acc->SetLineColor(kBlue);
  eff->Draw();
  acc->Draw("same");
  leg_accEff->AddEntry(acc, "Acceptance", "l");
  leg_accEff->AddEntry(eff, "Efficiency", "l");
  leg_accEff->Draw();
  
  TEfficiency *accCos = new TEfficiency(*acceptanceCos,*hDivRecoCos);
  TEfficiency *effCos = new TEfficiency(*efficiencyCos,*hDivPartonCos);
  TCanvas *can_effAccCos = new TCanvas("acceptance, efficiency can for Cos", "acceptance, efficiency can for Cos", 900, 600);
  TLegend *leg_accEffCos = new TLegend(0.5,0.6,0.7,0.9);
  effCos->SetTitle("Efficiency;|cos(#theta)|;Efficiency");
  accCos->SetTitle("Acceptance;|cos(#theta)|;Acceptance");
  effCos->SetLineColor(kRed);
  accCos->SetLineColor(kBlue);
  effCos->Draw();
  accCos->Draw("same");
  leg_accEffCos->AddEntry(acc, "Acceptance", "l");
  leg_accEffCos->AddEntry(eff, "Efficiency", "l");
  leg_accEffCos->Draw();

  
  
  acc->Write("chiAcceptance");
  eff->Write("chiEfficiency");
  stability->Write("chiStability");
  purity->Write("chiPurity");
  accCos->Write("cosAcceptance");
  effCos->Write("cosEfficiency");
  stabilityCos->Write("cosStability");
  purityCos->Write("cosPurity");
  can_purCos->Write("PurityStability_Cos_can");
  can_pur->Write("PurityStability_Chi_can");
  can_effAccCos->Write("EfficiencyAcceptance_Cos_can");
  can_effAcc->Write("EfficiencyAcceptance_Chi_can");
  
  
  //outf->Close();
}

TVector3 getBoostVector(TLorentzVector p4_1, TLorentzVector p4_2, TLorentzVector &p4CombinedVector)
{
	//define the combined Lorentz vector of ttbar 
	//TLorentzVector p4CombinedVector;
	p4CombinedVector.SetPxPyPzE((p4_1+p4_2).Px(),(p4_1+p4_2).Py(), (p4_1+p4_2).Pz(), (p4_1+p4_2).Energy()); 
	//get boost from this vector
	TVector3 TTbar_boostVector = p4CombinedVector.BoostVector();
	p4CombinedVector.Boost(-TTbar_boostVector);
	return -TTbar_boostVector;
}

