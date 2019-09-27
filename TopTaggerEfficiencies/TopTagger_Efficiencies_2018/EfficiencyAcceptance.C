#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"


using std::cin;
using std::cout; 
using std::endl;


std::vector<TString> listOfFiles;
std::vector<float> XSEC;
std::vector<TString> histoNames;
float LUMI = 59740;
TString eosPath; 
float deepCSVFloat = 0.4184;

void initFileNames()
{
  
  eosPath = "/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2018/Signal/";  
  listOfFiles.push_back("TT_Mtt-700to1000_TuneCP5_PSweights_13TeV-powheg-pythia8.root");
  listOfFiles.push_back("TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root");
}

void initXsections()
{
  XSEC.push_back(69.64);
  XSEC.push_back(16.74);
}

void initHistoNames()
{
  histoNames.push_back("Signal_histo_Mtt_700_1000");
  histoNames.push_back("Signal_histo_Mtt_1000_Inf");
}

void initGlobals()
{
  initFileNames();
  initXsections();
  initHistoNames();
}

void EfficiencyAcceptance(TString varName,TString varNameReco, int varNum , int sizeBins,float BND[] , float selMvaCut=0.3, float floatBTag = 0.8838, bool saveTtagger= true, float deepAK8CutValue = 0.6, float oldSelMVA = 0.8, bool isDeepCSV=true)
{
  initGlobals();
  
  TH1F *h_out_parton[3][listOfFiles.size()];
  TH1F *h_out_reco[3][listOfFiles.size()];
  
  TH1F* h_denominator_parton[listOfFiles.size()]; //one histogram for parton denominator because its the same for all types (tTagger, deepAK8, oldMva)
  TH1F* h_denominator_reco[3][listOfFiles.size()];
  std::vector<float> weights;
  
  //int sizeBins = 10;
  //float BND[11] = {1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, 4000, 5000};
  //int MIN = 1000; 
  //int MAX = 5000;
  
  for(int f=0; f<listOfFiles.size(); f++)
  {
    h_out_parton[0][f] = new TH1F(TString::Format("%s_%s_parton", histoNames[f].Data(), "tTagger"), TString::Format("%s_%s_parton", histoNames[f].Data(), "tTagger"), sizeBins, BND);
    h_out_parton[1][f] = new TH1F(TString::Format("%s_%s_parton", histoNames[f].Data(), "deepAK8"), TString::Format("%s_%s_parton", histoNames[f].Data(), "deepAK8"), sizeBins, BND);
	h_out_parton[2][f] = new TH1F(TString::Format("%s_%s_parton", histoNames[f].Data(), "oldMva"), TString::Format("%s_%s_parton", histoNames[f].Data(), "oldMva"), sizeBins, BND);
    
    h_out_reco[0][f] = new TH1F(TString::Format("%s_%s_reco", histoNames[f].Data(), "tTagger"), TString::Format("%s_%s_reco", histoNames[f].Data(), "tTagger"), sizeBins, BND);
    h_out_reco[1][f] = new TH1F(TString::Format("%s_%s_reco", histoNames[f].Data(), "deepAK8"), TString::Format("%s_%s_reco", histoNames[f].Data(), "deepAK8"), sizeBins, BND);
	h_out_reco[2][f] = new TH1F(TString::Format("%s_%s_reco", histoNames[f].Data(), "oldMva"), TString::Format("%s_%s_reco", histoNames[f].Data(), "oldMva"), sizeBins, BND);
	
	
    
    h_denominator_parton[f] = new TH1F(TString::Format("denom_%s_parton",histoNames[f].Data()), TString::Format("denom_%s_parton", histoNames[f].Data()), sizeBins, BND);
	h_denominator_parton[f]->Sumw2();
    
    h_denominator_reco[0][f] = new TH1F(TString::Format("denom_%s_%s_reco","tTagger",histoNames[f].Data()), TString::Format("denom_%s_%s_reco","tTagger", histoNames[f].Data()), sizeBins, BND);
    h_denominator_reco[1][f] = new TH1F(TString::Format("denom_%s_%s_reco","deepAK8", histoNames[f].Data()), TString::Format("denom_%s_%s_reco","deepAK8",histoNames[f].Data()), sizeBins, BND);
    h_denominator_reco[2][f] = new TH1F(TString::Format("denom_%s_%s_reco","oldMva",histoNames[f].Data()), TString::Format("denom_%s_%s_reco","oldMva", histoNames[f].Data()), sizeBins, BND);
	h_denominator_reco[0][f]->Sumw2();
	h_denominator_reco[1][f]->Sumw2();
	h_denominator_reco[2][f]->Sumw2();
    
    int nJets,nLeptons;
    float genEvtWeight;
    vector<bool>  *bit(0);
    vector<bool>  *matchedJet(0);
    vector<float> *jetPt(0),*tau3(0),*tau2(0),*tau1(0);
    vector<float> *jetMassSub0(0), *jetMassSub1(0);
    vector<float> *jetMassSoftDrop(0), *partonEta(0), *partonPt(0), *partonPhi(0), *partonMass(0),*partonMatchDR(0);

    float mva(0), yTTbarParton(0), ptTTbarParton(0);
    vector<float> *jetTtag(0);
    float mJJ(0), yJJ(0), ptJJ(0);
    int  category(0);
	std::vector<float> *deepAK8(0);
    std::vector<float> *jetPhi(0), *jetEta(0);
    std::vector<int> *addedIndexes = new std::vector<int>(0);
    std::vector<float> *jetBtagSub0(0), *jetBtagSub1(0);
    vector<int> *partonId(0), *partonMatchIdx(0);
    float mTTbarParton(0);
    std::vector<float> *jetBtagSub0DCSVbb(0), *jetBtagSub1DCSVbb(0);
    std::vector<float> *jetBtagSub0DCSVbbb(0), *jetBtagSub1DCSVbbb(0);
	
    std::cout<<"Working in file: "<<listOfFiles[f]<<std::endl;
    TFile *file = TFile::Open(eosPath+listOfFiles[f]);
    TTree *trIN = (TTree*)file->Get("boosted/events");
    	
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
    trIN->SetBranchAddress("mJJ"   			,&mJJ);
    trIN->SetBranchAddress("yJJ"   			,&yJJ);
    trIN->SetBranchAddress("ptJJ"   		,&ptJJ);
	
	trIN->SetBranchAddress("mTTbarParton"	,&mTTbarParton);
	trIN->SetBranchAddress("yTTbarParton"	,&yTTbarParton);
	trIN->SetBranchAddress("ptTTbarParton"	,&ptTTbarParton);
	trIN->SetBranchAddress("partonPt"	    ,&partonPt);
	trIN->SetBranchAddress("partonEta"		,&partonEta);
	trIN->SetBranchAddress("partonMass"     ,&partonMass);
	trIN->SetBranchAddress("partonMatchDR"  ,&partonMatchDR);
	trIN->SetBranchAddress("partonMatchIdx" ,&partonMatchIdx);
	trIN->SetBranchAddress("partonPhi"      ,&partonPhi);
    
	trIN->SetBranchAddress("jetBtagSub0"	,&jetBtagSub0);
    trIN->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1);
    trIN->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
  
    trIN->SetBranchAddress("mva"	  		,&mva);
    trIN->SetBranchAddress("category"	  	,&category);
	trIN->SetBranchAddress("deepAK8Tagger"  ,&deepAK8);
    trIN->SetBranchAddress("jetTtagCategory",&jetTtag);
	//deepCSV
    trIN->SetBranchAddress("jetBtagSub0DCSVbb" ,&jetBtagSub0DCSVbb);
    trIN->SetBranchAddress("jetBtagSub1DCSVbb" ,&jetBtagSub1DCSVbb);
    trIN->SetBranchAddress("jetBtagSub0DCSVbbb",&jetBtagSub0DCSVbbb);
    trIN->SetBranchAddress("jetBtagSub1DCSVbbb",&jetBtagSub1DCSVbbb);
	
	
    float norm = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
	float weight = XSEC[f]/norm;
	weights.push_back(weight);
    
	//for matching
	std::vector<int> *jetMatchedIndexes = new std::vector<int>(0);
	std::vector<float> *jetMatchedDr = new std::vector<float>(0);
	std::vector<float> *eta_ = new std::vector<float>(0);
	std::vector<float> *phi_ = new std::vector<float>(0);
	std::vector<float> *mass_ = new std::vector<float>(0);
	std::vector<float> *pt_ = new std::vector<float>(0);
	std::vector<float> *deepAK8_ = new std::vector<float>(0);
	std::vector<float> *jetTtag_ = new std::vector<float>(0);
	std::vector<float> *jetBtagSub0_ = new std::vector<float>(0);
	std::vector<float> *jetBtagSub1_ = new std::vector<float>(0);

	std::vector<float> *partonPt_ = new std::vector<float>(0);
	std::vector<float> *partonEta_ = new std::vector<float>(0);
	std::vector<float> *partonMass_ = new std::vector<float>(0);
	std::vector<float> *partonPhi_ = new std::vector<float>(0);
	
	std::vector<float> *jetBtagSub0DCSVbb_ = new std::vector<float>(0);
	std::vector<float> *jetBtagSub1DCSVbb_ = new std::vector<float>(0);
	std::vector<float> *jetBtagSub0DCSVbbb_ = new std::vector<float>(0);
	std::vector<float> *jetBtagSub1DCSVbbb_ = new std::vector<float>(0);
	
	float jetDr_(0);
	
	
    int decade(0);
    int NN = trIN->GetEntries();
    //NN = 100000;
    std::cout<<"Entries: "<<NN<<std::endl;
	
	bool isJetVar;
	if(varNum > 3)
	{
		isJetVar = true;
	}
	float xParton(0), xReco(0);
	float xJetParton[2], xJetReco[2];
	 
	
    for(int iev=0;iev<NN;iev++) 
    {
		double progress = 10.0*iev/(1.0*NN);
      int k = TMath::FloorNint(progress); 
      if (k > decade) 
        cout<<10*k<<" %"<<endl;
      decade = k;
      trIN->GetEntry(iev);
		int isMatched =0;
		eta_->clear();
        mass_->clear();
        pt_->clear();
        phi_->clear();
        jetBtagSub0_->clear();
        jetBtagSub1_->clear();
        jetTtag_->clear();
        deepAK8_->clear();
		
       partonPt_->clear();
       partonMass_->clear();
       partonPhi_->clear();
       partonEta_->clear();
	  
	   jetBtagSub0DCSVbb_->clear();
	   jetBtagSub1DCSVbb_->clear();
	   jetBtagSub0DCSVbbb_->clear();
	   jetBtagSub1DCSVbbb_->clear();
	   
		if (varNum ==1)
		{
			xParton = mTTbarParton;
			xReco   = mJJ;
		}
		else if (varNum ==2)
		{
			xParton = ptTTbarParton;
			xReco   = ptJJ;
		}
		else if (varNum ==3) 
		{
			xParton = yTTbarParton;
			xReco   = yJJ;
		}
		
      
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
					//cout<<"dRdRmin: "<<dRmin<<endl;
					//cout<<indexMin<<endl;
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
						deepAK8_->push_back( (*deepAK8)[(*partonMatchIdx)[indexMin]]);
						
						jetBtagSub0DCSVbb_->push_back((*jetBtagSub0DCSVbb)[(*partonMatchIdx)[indexMin]]);
						jetBtagSub1DCSVbb_->push_back((*jetBtagSub1DCSVbb)[(*partonMatchIdx)[indexMin]]);
						jetBtagSub0DCSVbbb_->push_back((*jetBtagSub0DCSVbbb)[(*partonMatchIdx)[indexMin]]);
						jetBtagSub1DCSVbbb_->push_back((*jetBtagSub1DCSVbbb)[(*partonMatchIdx)[indexMin]]);

						partonPt_->push_back( (*partonPt)[indexMin]);
						partonMass_->push_back( (*partonMass)[indexMin]);
						partonPhi_->push_back( (*partonPhi)[indexMin]);
						partonEta_->push_back( (*partonEta)[indexMin]);
						
						
						//cout<<(*partonMatchIdx)[indexMin]<<endl;
						//cout<<"------------"<<endl;
					}
					  
			   }
			   
			}
			if(isMatched >1)
			{
				
				if (varNum ==4)
				{
					xJetParton[0] = (*partonPt_)[0];
					xJetParton[1] = (*partonPt_)[1];
					
					xJetReco[0] = (*pt_)[0];
					xJetReco[1] = (*pt_)[1];
				}
				else if (varNum ==5)
				{					
					xJetParton[0] = (*partonEta_)[0];
					xJetParton[1] = (*partonEta_)[1];
				
					xJetReco[0] = (*eta_)[0];
					xJetReco[1] = (*eta_)[1];
				}
			}
		
	  //---------------------------end of MATCHING---------------------------------------------------------
	  bool recoCuts, partonCuts; 
	  bool massCut = (*mass_)[0] > 120 && (*mass_)[0] < 220 && (*mass_)[1] > 120 && (*mass_)[1] < 220;
	  bool CSVv2Cut = ((*jetBtagSub0_)[0] > floatBTag || (*jetBtagSub1_)[0] > floatBTag) && ((*jetBtagSub0_)[1] > floatBTag || (*jetBtagSub1_)[1] > floatBTag);
	  bool tTaggerCut = (*jetTtag_)[0] > selMvaCut && (*jetTtag_)[1] > selMvaCut;
	  bool deepAK8Cut = (*deepAK8_)[0] > deepAK8CutValue && (*deepAK8_)[1] > deepAK8CutValue;
	  recoCuts = nJets > 1 && fabs((*eta_)[0]) < 2.4 && fabs((*eta_)[1] <2.4) && (*pt_)[0] > 400 && (*pt_)[1] > 400 && mJJ > 1000 && massCut && nLeptons==0;
	  partonCuts = fabs((*partonEta_)[0]) < 2.4 && fabs((*partonEta_)[1] <2.4) && (*partonPt_)[0] > 400 && (*partonPt_)[1] > 400 && mTTbarParton > 1000;
	  bool deepCSV = (((*jetBtagSub0DCSVbb_)[0] + (*jetBtagSub0DCSVbbb_)[0])> deepCSVFloat || ((*jetBtagSub1DCSVbb_)[0] + (*jetBtagSub1DCSVbbb_)[0])> deepCSVFloat) && 
					 (((*jetBtagSub0DCSVbb_)[1] + (*jetBtagSub0DCSVbbb_)[1])> deepCSVFloat || ((*jetBtagSub1DCSVbb_)[1] + (*jetBtagSub1DCSVbbb_)[1])> deepCSVFloat);
			
      bool btagCut;			
	  if(isDeepCSV) btagCut = deepCSV;
	  else btagCut = CSVv2Cut;
	  
	  /*
	  cout<<"-------------------------------------------------------"<<endl;
	  cout<<"jet 0 sub 0: "<<((*jetBtagSub0DCSVbb_)[0] + (*jetBtagSub0DCSVbbb_)[0])<<endl;
	  cout<<"jet 1 sub 0: "<<((*jetBtagSub0DCSVbb_)[1] + (*jetBtagSub0DCSVbbb_)[1])<<endl;
	  cout<<"jet 0 sub 1: "<<((*jetBtagSub1DCSVbb_)[0] + (*jetBtagSub1DCSVbbb_)[0])<<endl;
	  cout<<"jet 1 sub 1: "<<((*jetBtagSub1DCSVbb_)[1] + (*jetBtagSub1DCSVbbb_)[1])<<endl;
	  
     if( isMatched > 1 )
	 {
	  
	  cout<<"------------------------------"<<endl;
	  cout<<"(*pt_)[0]: "<< (*pt_)[0]<<endl;
	  cout<<"(*pt_)[1]: "<< (*pt_)[1]<<endl;
	 
	  for(int ijet = 0; ijet< nJets; ijet++)
	  {
		 cout<<"(*jetPt)["<<ijet<<"]: "<< (*jetPt)[ijet]<<endl;
	  }
	  cout<<"(*partonPt_)[0]: "<<(*partonPt_)[0]<<endl;
	  cout<<"(*partonPt_)[1]: "<<(*partonPt_)[1]<<endl;
	  cout<<"partonPt[0]: "<<(*partonPt)[0]<<endl;
	  cout<<"partonPt[1]: "<<(*partonPt)[1]<<endl;
	  cout<<"(*partonMatchIdx)[0]: "<<(*partonMatchIdx)[0]<<endl;
	  cout<<"(*partonMatchIdx)[1]: "<<(*partonMatchIdx)[1]<<endl;
	  }*/
	  
	  
	  	 
	  //fill the denominators
	  //1. denominator passing only reco cuts for topTagger
	  if(recoCuts && btagCut && tTaggerCut)
	  {
	    if(!isJetVar)
		{
			h_denominator_reco[0][f]->Fill(xReco, genEvtWeight);
		}
	    else
		{	   
			  h_denominator_reco[0][f]->Fill(xJetReco[0], genEvtWeight);
			  h_denominator_reco[0][f]->Fill(xJetReco[1], genEvtWeight);
		}
	  }
	  // denominator passing only reco cuts for deepAK8
	  if(recoCuts && btagCut && deepAK8Cut)
	  {
	    if(!isJetVar)
		{
			h_denominator_reco[1][f]->Fill(xReco, genEvtWeight);
		}
	    else
		{	   
			h_denominator_reco[1][f]->Fill(xJetReco[0], genEvtWeight);
			h_denominator_reco[1][f]->Fill(xJetReco[1], genEvtWeight);
		}
	  }
	  // denominator passing only reco cuts for old mva
	  if(recoCuts && category==2 && mva >0.8)
	  {
	    if(!isJetVar)
		{
			h_denominator_reco[2][f]->Fill(xReco, genEvtWeight);
		}
	    else
		{	   
			h_denominator_reco[2][f]->Fill(xJetReco[0], genEvtWeight);
			h_denominator_reco[2][f]->Fill(xJetReco[1], genEvtWeight);
		}
	  }
	  //2. fill the histograms pass reco and parton cuts numerators for efficiencies and acceptance
      if(partonCuts && recoCuts)
	  {
        //1. category where both are fully top tagged
         //top tagger efficiency
	      if(tTaggerCut && btagCut)
		  {
			
			if(!isJetVar) 
			{
			   h_out_parton[0][f]->Fill(xParton, genEvtWeight);
			   h_out_reco[0][f]->Fill(xReco, genEvtWeight);
			}
	     	else
			{
			  if(isMatched > 1)
			  {
				h_out_parton[0][f]->Fill(xJetParton[0], genEvtWeight);
				h_out_parton[0][f]->Fill(xJetParton[1], genEvtWeight);
				h_out_reco[0][f]->Fill(xJetReco[0], genEvtWeight);
	            h_out_reco[0][f]->Fill(xJetReco[1], genEvtWeight); 
			  }   
			}
			
		  }//---- end of tTagger
		  
          //deep AK8 tagger
          if(deepAK8Cut && btagCut)
          {
			if(!isJetVar) 
			{
			   h_out_parton[1][f]->Fill(xParton, genEvtWeight);
			   h_out_reco[1][f]->Fill(xReco, genEvtWeight);
			}
	     	else
			{
			  if(isMatched > 1)
			  {
				h_out_parton[1][f]->Fill(xJetParton[0], genEvtWeight);
				h_out_parton[1][f]->Fill(xJetParton[1], genEvtWeight);
				h_out_reco[1][f]->Fill(xJetReco[0], genEvtWeight);
	            h_out_reco[1][f]->Fill(xJetReco[1], genEvtWeight); 
			  }   
			}
          }//---- end of deepAK8
		  
		  //old tagger kkousour
		  if(mva > oldSelMVA && category==2)
          {
            if(!isJetVar) 
			{
			   h_out_parton[2][f]->Fill(xParton, genEvtWeight);
			   h_out_reco[2][f]->Fill(xReco, genEvtWeight);
			}
	     	else
			{
			  if(isMatched > 1)
			  {
				h_out_parton[2][f]->Fill(xJetParton[0], genEvtWeight);
				h_out_parton[2][f]->Fill(xJetParton[1], genEvtWeight);
				h_out_reco[2][f]->Fill(xJetReco[0], genEvtWeight);
	            h_out_reco[2][f]->Fill(xJetReco[1], genEvtWeight); 
			  }   
			}
          }//----end of old mva
		  	
      }//----- end of selection cuts parton and reco 
	  
	  
    }//---end the event for
	
	//--------------------------------------------START OF EVENT COUNTER LOOP -------------------------------------------------------------------

  //now another for that fills the denominators for the parton efficiencies 
  //loop over other tree -> eventCounter
  TTree *trCnt = (TTree*)file->Get("eventCounter/events");
  float ptTTbarPartonCnt(0), mTTbarPartonCnt(0), yTTbarPartonCnt(0);
  float partonPtCnt[2], partonEtaCnt[2],partonYCnt[2],genEvtWeightCnt; 	
  //tree for eventCounter		
  trCnt->SetBranchAddress("ptTopParton"    ,&partonPtCnt);
  trCnt->SetBranchAddress("etaTopParton"   ,&partonEtaCnt);
  trCnt->SetBranchAddress("yTopParton"     ,&partonYCnt);
  trCnt->SetBranchAddress("mTTbarParton"   ,&mTTbarPartonCnt);
  trCnt->SetBranchAddress("yTTbarParton"   ,&yTTbarPartonCnt);
  trCnt->SetBranchAddress("ptTTbarParton"  ,&ptTTbarPartonCnt);
  trCnt->SetBranchAddress("genEvtWeight"   ,&genEvtWeightCnt);   
  int NNCnt = trCnt->GetEntries();
  
  for(int iev = 0; iev < NNCnt; iev++)
  {
	  trCnt->GetEntry(iev);
	  if (varNum ==1) xParton = mTTbarPartonCnt;
	  else if (varNum ==2) xParton = ptTTbarPartonCnt;
	  else if (varNum ==3) xParton = yTTbarPartonCnt;
	  else if (varNum ==4)
	  {
	   xJetParton[0] = partonPtCnt[0];
	   xJetParton[1] = partonPtCnt[1];
	  }
	  else if (varNum ==5)
	  {					
 	   xJetParton[0] = partonEtaCnt[0];
	   xJetParton[1] = partonEtaCnt[1];
	  }
	  bool partonCuts = fabs(partonEtaCnt[0]) < 2.4 && fabs(partonEtaCnt[1]) <2.4 && partonPtCnt[0] > 400 && partonPtCnt[1] > 400 && mTTbarPartonCnt > 1000;
		
	  if(partonCuts)
	  {
		if(!isJetVar) h_denominator_parton[f]->Fill(xParton, genEvtWeightCnt);
		else 
		{
	    	h_denominator_parton[f]->Fill(xJetParton[0], genEvtWeightCnt);
			h_denominator_parton[f]->Fill(xJetParton[1], genEvtWeightCnt);
		}
			
	  }  
  }

  //--------------------------------------------END OF EVENT COUNTER LOOP -------------------------------------------------------------------
}
  
  
  
  
  for(int i=0; i<sizeof(h_out_reco)/sizeof(*h_out_reco); i++)
  {
	//for every slice
	for(int j=0; j<sizeof(*h_out_reco)/sizeof(*h_out_reco[0]); j++)
    {
      h_out_parton[i][j]->Scale(weights[j]);
      h_out_reco[i][j]->Scale(weights[j]);
    }
    
    for(int j=1; j<sizeof(*h_out_reco)/sizeof(*h_out_reco[0]); j++)
	{
	  //Add them to get the whole phase space
	  h_out_parton[i][0]->Add(h_out_parton[i][j]);
      std::cout<<h_out_reco[i][0]->GetName()<<" "<<h_out_reco[i][j]->GetName()<<std::endl;
      h_out_reco[i][0]->Add(h_out_reco[i][j]);
	}
  }
  
  
  for(int i=0; i<sizeof(h_denominator_reco)/sizeof(*h_denominator_reco); i++)
  {
	//for every slice
	for(int j=0; j<sizeof(*h_denominator_reco)/sizeof(*h_denominator_reco[0]); j++)
    {
      h_denominator_reco[i][j]->Scale(weights[j]);
    }
    
    for(int j=1; j<sizeof(*h_denominator_reco)/sizeof(*h_denominator_reco[0]); j++)
	{
	  //Add them to get the whole phase space
      h_denominator_reco[i][0]->Add(h_out_reco[i][j]);
	}
  }
  
  h_denominator_parton[0]->Scale(weights[0]);
  for(int i=1; i<sizeof(h_denominator_parton)/sizeof(*h_denominator_parton); i++) 
  {
    h_denominator_parton[i]->Scale(weights[i]);
	h_denominator_parton[0]->Add(h_denominator_parton[i]);
  }
  
  TEfficiency *efficiency_parton[3], *efficiency_reco[3];
 
  //efficiency for parton quantity and for topTagger (new)
  efficiency_parton[0]  = new TEfficiency(*h_out_parton[0][0], *h_denominator_parton[0]);
  if (varNum ==1 || varNum ==2 || varNum ==4) efficiency_parton[0]->SetTitle(TString::Format("Parton_tTagger;%s (GeV);Efficiency",varName.Data()));
  else efficiency_parton[0]->SetTitle(TString::Format("Parton_tTagger;%s ;Efficiency",varName.Data()));
  efficiency_parton[0]->SetStatisticOption(TEfficiency::kFNormal);
  efficiency_parton[0]->SetUseWeightedEvents();
  efficiency_parton[0]->SetLineColor(kRed);
	  
  //efficiency for parton quantity and for deepAK8
  efficiency_parton[1]  = new TEfficiency(*h_out_parton[1][0], *h_denominator_parton[0]);
  if (varNum ==1 || varNum ==2 || varNum ==4) efficiency_parton[1]->SetTitle(TString::Format("Parton_deepAK8; %s (GeV);Efficiency",varName.Data()));
  else efficiency_parton[1]->SetTitle(TString::Format("Parton_deepAK8; %s;Efficiency",varName.Data()));
  efficiency_parton[1]->SetStatisticOption(TEfficiency::kFNormal);
  efficiency_parton[1]->SetUseWeightedEvents();
  efficiency_parton[1]->SetLineColor(kBlue);
 
  //efficiency for parton quantity and for old top tagger mva (old)
  efficiency_parton[2]  = new TEfficiency(*h_out_parton[2][0], *h_denominator_parton[0]);
  if (varNum ==1 || varNum ==2 || varNum ==4) efficiency_parton[2]->SetTitle(TString::Format("Parton_oldMva; %s (GeV);Efficiency",varName.Data()));
  else efficiency_parton[2]->SetTitle(TString::Format("Parton_oldMva; %s;Efficiency",varName.Data()));
  efficiency_parton[2]->SetStatisticOption(TEfficiency::kFNormal);
  efficiency_parton[2]->SetUseWeightedEvents();
  efficiency_parton[2]->SetLineColor(kMagenta);

  //efficiency for reco quantity and for topTagger (new)
  efficiency_reco[0]  = new TEfficiency(*h_out_reco[0][0], *h_denominator_reco[0][0]);
  if (varNum ==1 || varNum ==2 || varNum ==4)  efficiency_reco[0]->SetTitle(TString::Format("Reco_tTagger; %s(GeV);Efficiency",varNameReco.Data()));
  else efficiency_reco[0]->SetTitle(TString::Format("Reco_tTagger; %s;Efficiency", varNameReco.Data()));
  efficiency_reco[0]->SetStatisticOption(TEfficiency::kFNormal);
  efficiency_reco[0]->SetUseWeightedEvents();
  efficiency_reco[0]->SetLineColor(kRed);

  //efficiency for reco quantity and for deepAK8  
  efficiency_reco[1]  = new TEfficiency(*h_out_reco[1][0], *h_denominator_reco[1][0]);
  if (varNum ==1 || varNum ==2 || varNum ==4) efficiency_reco[1]->SetTitle(TString::Format("Reco_deepAK8;%s (GeV);Efficiency",varNameReco.Data()));
  else efficiency_reco[1]->SetTitle(TString::Format("Reco_deepAK8; %s;Efficiency",varNameReco.Data()));
  efficiency_reco[1]->SetStatisticOption(TEfficiency::kFNormal);
  efficiency_reco[1]->SetUseWeightedEvents();
  efficiency_reco[1]->SetLineColor(kBlue);  
  
  //efficiency for parton quantity and for old top tagger mva (old)	
  efficiency_reco[2]  = new TEfficiency(*h_out_reco[2][0], *h_denominator_reco[2][0]);
  if (varNum ==1 || varNum ==2 || varNum ==4) efficiency_reco[2]->SetTitle(TString::Format("Reco_oldMva;%s (GeV);Efficiency",varNameReco.Data()));
  else efficiency_reco[2]->SetTitle(TString::Format("Reco_deepAK8; %s;Efficiency",varNameReco.Data()));
  efficiency_reco[2]->SetStatisticOption(TEfficiency::kFNormal);
  efficiency_reco[2]->SetUseWeightedEvents();
  efficiency_reco[2]->SetLineColor(kMagenta);  
  

  TCanvas *c1 = new TCanvas("partonCanvas", "partonCanvas", 600, 500);
  c1->cd();
  efficiency_parton[0]->Draw();
  efficiency_parton[1]->Draw("SAME");
  efficiency_parton[2]->Draw("SAME");

  
  TCanvas *c2 = new TCanvas("recoCanvas", "recoCanvas", 600, 500);
  c2->cd();
  efficiency_reco[0]->Draw();
  efficiency_reco[1]->Draw("SAME");
  efficiency_reco[2]->Draw("SAME");
  
  TFile *outFile;
  if(isDeepCSV) outFile = TFile::Open(TString::Format("deepAK8_efficiencies_allVars_tTagger_%0.2f_deepCSV.root", selMvaCut), "UPDATE");
  else outFile = TFile::Open(TString::Format("deepAK8_efficiencies_allVars_tTagger_%0.2f.root", selMvaCut), "UPDATE");
  
  
  if(saveTtagger)
  {
	//save the topTagger new efficiency
    efficiency_parton[0]->Write(TString::Format("Sig_Parton_tTagger_%0.2f_%s", selMvaCut,  varName.Data()));
    efficiency_reco[0]->Write(TString::Format("Sig_Reco_tTagger_%0.2f_%s", selMvaCut,  varNameReco.Data()));
	
	//save the topTagger mva old efficiency
	efficiency_parton[2]->Write(TString::Format("Sig_Parton_oldMva_%0.2f_%s", oldSelMVA,  varName.Data()));
    efficiency_reco[2]->Write(TString::Format("Sig_Reco_oldMva_%0.2f_%s", oldSelMVA,  varNameReco.Data()));
  }
  efficiency_parton[1]->Write(TString::Format("Sig_Parton_deepAK8_%0.2f_%s", deepAK8CutValue,  varName.Data()));
  efficiency_reco[1]->Write(TString::Format("Sig_Reco_deepAK8_%0.2f_%s", deepAK8CutValue, varNameReco.Data()));
  outFile->Close();


  listOfFiles.clear();
  histoNames.clear();
  XSEC.clear();
}
