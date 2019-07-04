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


void atlasComparison(TString file = "/eos/cms/store/user/ipapakri/ttbar/MC/Signal/TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8_legacy2016_deepAK8.root", 
						float selMvaCut=0.1, float floatBTag = 0.8838, bool isDeepCSV = true, bool isZprime= false,bool isParton=false, int ZprimeMass = 2000, TString width = "200" )
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
  vector<bool> *bit(0);
  float mTTbarParton(0),mJJ(0), yTTbarParton(0), ptTTbarParton(0);
  int  category(0);
  //matching info 
  vector<float> *jetPhi(0), *jetEta(0);
  vector<int> *partonId(0), *partonMatchIdx(0);
  float yTopParton[2];
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
  trIN->SetBranchAddress("yTopParton"	  ,&yTopParton);
  //deepCSV
  trIN->SetBranchAddress("jetBtagSub0DCSVbb" ,&jetBtagSub0DCSVbb);
  trIN->SetBranchAddress("jetBtagSub1DCSVbb" ,&jetBtagSub1DCSVbb);
  trIN->SetBranchAddress("jetBtagSub0DCSVbbb",&jetBtagSub0DCSVbbb);
  trIN->SetBranchAddress("jetBtagSub1DCSVbbb",&jetBtagSub1DCSVbbb);
  
  
  
  TLorentzVector p4T[2], p4TTbar, p4T_ZMF[2];

  int decade(0);
  int NN = trIN->GetEntries();
  
  const int chiSize =9;
  float BND_chi[chiSize+1] = {1,1.5,2,3,4,5,7,10,13,16};
  const int cosSize = 6;
  float BND_cos[cosSize+1] = {0,0.2,0.4,0.6,0.7,0.8,1};
  TH1F *h_Chi_all = new TH1F("#chi dist", "#chi dist", chiSize, BND_chi);
  TH1F *h_Cos_all = new TH1F("#||{cos(#theta)} dist", "#||{cos(#theta)} dist", cosSize, BND_cos);
    
  TH1F *hChi_ATLAS = new TH1F("#chiAtlas_dist", "#chiAtlas_dist", chiSize, BND_chi);
  TH1F *hCos_ATLAS = new TH1F("#||{cos(#theta)}_Atlas_dist", "#||{cos(#theta)}_Atlas_dist", cosSize, BND_cos);
  

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
  float jetDr_(0);
  
  cout<<"Reading "<<NN<<" entries"<<endl;
  for(int iev=0;iev<NN;iev++) 
  {
    double progress = 10.0*iev/(1.0*NN);
    int k = TMath::FloorNint(progress); 
    if (k > decade) 
      cout<<10*k<<" %"<<endl;
    decade = k;
    trIN->GetEntry(iev);
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
	
	int isMatched=0;
	bool recoCuts, partonCuts, btagging, topTagger, massCut, btaggingReverted;	
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
					
					partonPt_->push_back( (*partonPt)[indexMin]);
					partonMass_->push_back( (*partonMass)[indexMin]);
					partonPhi_->push_back( (*partonPhi)[indexMin]);
					partonEta_->push_back( (*partonEta)[indexMin]);
					
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
				
				recoCuts   = fabs((*eta_)[0]) < 2.4 && fabs((*eta_)[1]) <2.4 && (*pt_)[0] > 500 && (*pt_)[1] > 500 && nLeptons==0 && mJJ > 1000;
				partonCuts = (*partonPt_)[0] > 500 && (*partonPt_)[1] > 500 && fabs((*partonEta_)[0]) < 2.4 && fabs((*partonEta_)[1]) < 2.4 &&  mTTbarParton > 1000;
				CSVv2Cut   = ((*jetBtagSub0_)[0] > floatBTag || (*jetBtagSub1_)[0] > floatBTag) && ((*jetBtagSub0_)[1] > floatBTag || (*jetBtagSub1_)[1] > floatBTag);
				deepCSVCut = (dCSVScoreSub0[0] > deepCSVFloat || dCSVScoreSub1[0] > deepCSVFloat) && (dCSVScoreSub1[1] > deepCSVFloat || dCSVScoreSub1[1] > deepCSVFloat);
				topTagger  = (*jetTtag_)[0] > selMvaCut && (*jetTtag_)[1] > selMvaCut;
				massCut    = (*mass_)[0] > 120 && (*mass_)[0] < 220 && (*mass_)[1] > 120 && (*mass_)[1] < 220;
				//btaggingReverted = category==0;
				//btaggingReverted = ((*jetBtagSub0_)[0] < floatBTag && (*jetBtagSub1_)[0] < floatBTag) && ((*jetBtagSub0_)[1] < floatBTag && (*jetBtagSub1_)[1] < floatBTag);
				if(isDeepCSV) btagging = deepCSVCut;
				else btagging = CSVv2Cut;
				//cout<<"ok"<<endl;
					
					//split the mttbar phase space into regions and for each region I will calculate the thetas dists
					bool currentSelection;
					if(isParton)
					{
						p4T[0].SetPtEtaPhiM((*partonPt_)[0], (*partonEta_)[0], (*partonPhi_)[0], (*partonMass_)[0]);
						p4T[1].SetPtEtaPhiM((*partonPt_)[1], (*partonEta_)[1], (*partonPhi_)[1], (*partonMass_)[1]);
						currentSelection = partonCuts;
					}
					else
					{
						p4T[0].SetPtEtaPhiM((*pt_)[0], (*eta_)[0], (*phi_)[0], (*mass_)[0]);
						p4T[1].SetPtEtaPhiM((*pt_)[1], (*eta_)[1], (*phi_)[1], (*mass_)[1]);
						currentSelection = recoCuts && topTagger && btagging;
					}
					
					TVector3 ttbarBoostVector = getBoostVector(p4T[0], p4T[1], p4TTbar);
					
					p4T_ZMF[0].SetPtEtaPhiM(p4T[0].Pt(), p4T[0].Eta(), p4T[0].Phi(), p4T[0].M());
					p4T_ZMF[1].SetPtEtaPhiM(p4T[1].Pt(), p4T[1].Eta(), p4T[1].Phi(), p4T[1].M());
					p4T_ZMF[0].Boost(ttbarBoostVector);
					p4T_ZMF[1].Boost(ttbarBoostVector);

					float chi0(0);
					chi0 = (1 + fabs(TMath::Cos(p4T_ZMF[0].Theta()))) / ( 1 - fabs(TMath::Cos(p4T_ZMF[0].Theta())));
					float yStarExp  = TMath::Exp(fabs(p4T_ZMF[0].Rapidity() - p4T_ZMF[1].Rapidity()));
					
					if(currentSelection)
					{
						h_Chi_all->Fill(yStarExp,genEvtWeight);	
						h_Cos_all->Fill(fabs(TMath::Cos(p4T_ZMF[0].Theta())),genEvtWeight);
				
						
						if(chi0-yStarExp < 0) cout<<"Found less that 0"<<endl;
	
				    }//----end of cuts

				   
				
		}//----end of isMatched	
	}//----end of nJets
  }	//---end of event loop
  
  float chiAtlas[chiSize] = {0.606788, 0.39109, 0.230487, 0.121592, 0.0689784, 0.033708, 0.011961, 0,0};
  float cosAtlas[cosSize] = {1.57669,1.47934,1.26476,0.93499,0.611159,0.140864};
  
  float chiAtlasErrors[chiSize] = {0.0027431,0.00167096,0.000785063,0.000495496,0.000329766,0.000168031,7.31112E-5,0,0};
  float cosAtlasErrors[cosSize] = {0.00703382,0.00499111,0.0035704,0.00380381,0.00277803,0.00101854};
  //set Bin Content for atlas distributions and errors
  for(int i =0; i<chiSize; i++)
  {
	  hChi_ATLAS->SetBinContent(i+1, chiAtlas[i]);
	  hChi_ATLAS->SetBinError(i+1, chiAtlasErrors[i]);
  }
  for(int i =0; i<cosSize; i++)
  {
	  hCos_ATLAS->SetBinContent(i+1, cosAtlas[i]);
	  hCos_ATLAS->SetBinError(i+1, cosAtlasErrors[i]);
  }
  
  TCanvas *c1 = new TCanvas("hChi over hChiATLAS", "hChi over hChiATLAS", 900,600);
  TLegend *leg_rp = new TLegend(0.6,0.7,0.8,0.9);
  leg_rp->AddEntry(h_Chi_all,"#chi CMS", "l");
  leg_rp->AddEntry(hChi_ATLAS ,"#chi ATLAS", "l");
  h_Chi_all->Scale(1./h_Chi_all->Integral(),"width");
  h_Chi_all->GetXaxis()->SetTitle("#chi");
  hChi_ATLAS->GetXaxis()->SetTitle("#chi");
  hChi_ATLAS->SetLineColor(kRed);
  
  TRatioPlot *rp = new TRatioPlot(hChi_ATLAS, h_Chi_all);
  c1->SetTicks(0,1);
  rp->Draw();
  leg_rp->Draw();
  c1->Update();
  
  TCanvas *c2 = new TCanvas("hCosTheta over hCostThetaATLAS", "hCosTheta over hCostThetaATLAS", 900,600);
  TLegend *leg_rp2 = new TLegend(0.6,0.7,0.8,0.9);
  leg_rp2->AddEntry(h_Cos_all,"#chi CMS", "l");
  leg_rp2->AddEntry(hCos_ATLAS ,"#chi ATLAS", "l");
  h_Cos_all->Scale(1./h_Cos_all->Integral(),"width");
  h_Cos_all->GetXaxis()->SetTitle("|cos(#theta^{*})|");
  hCos_ATLAS->GetXaxis()->SetTitle("|cos(#theta^{*})|");
  hCos_ATLAS->SetLineColor(kRed);
  
  TRatioPlot *rp2 = new TRatioPlot(hCos_ATLAS, h_Cos_all);
  c2->SetTicks(0,1);
  rp2->Draw();
  leg_rp2->Draw();
  c2->Update();

/*
  TString recoParton ="";
  if(isParton) recoParton = "Parton";
  else recoParton = "Reco";
  TFile *outf = new TFile(TString::Format("Output_TT_QCD_%s_Chi_%0.1f.root", recoParton.Data(),selMvaCut), "RECREATE");
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

