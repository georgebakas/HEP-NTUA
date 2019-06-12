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


std::vector<TString> listOfFiles;
std::vector<float> XSEC;
std::vector<TString> histoNames;
float LUMI = 35900;
TString eosPath;
bool isSignal;

void initFileNames()
{
  if(isSignal)
  {
	eosPath = "/eos/cms/store/user/ipapakri/ttbar/MC/Signal/";  
	listOfFiles.push_back("TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8_legacy2016_deepAK8.root");
	listOfFiles.push_back("TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8_legacy2016_deepAK8.root");
  }
  else
  {
	eosPath = "/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/";
	listOfFiles.push_back("QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
	listOfFiles.push_back("QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
	listOfFiles.push_back("QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
	listOfFiles.push_back("QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
	listOfFiles.push_back("QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
	listOfFiles.push_back("QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
  }
}

void initXsections()
{
  if(isSignal)
  {
	XSEC.push_back(69.64);
	XSEC.push_back(16.74);
  }
  else
  {
	XSEC.push_back(3.67e+5);
	XSEC.push_back(2.94e+4);
	XSEC.push_back(6.524e+03);
	XSEC.push_back(1.064e+03);
	XSEC.push_back(121.5);
	XSEC.push_back(2.542e+01);
  }
}

void initHistoNames()
{
  if(isSignal)
  {
	histoNames.push_back("Signal_histo_Mtt_700_1000"); 
	histoNames.push_back("Signal_histo_Mtt_1000_Inf");
  }
  else
  {
	histoNames.push_back("QCD_histo_Mtt_300_500");
	histoNames.push_back("QCD_histo_Mtt_500_700");
	histoNames.push_back("QCD_histo_Mtt_700_1000");
	histoNames.push_back("QCD_histo_Mtt_1000_1500");
	histoNames.push_back("QCD_histo_Mtt_1500_2000");
	histoNames.push_back("QCD_histo_Mtt_2000_Inf");
  }
}

void initGlobals()
{
  initFileNames();
  initXsections();
  initHistoNames();
}
 
void qcdClosure(TString varName,TString varNameReco, int varNum,  int sizeBins,float BND[] ,  bool saveTTagger= true, float selMvaCut=0.2, float floatBTag = 0.8838, bool 				isSig = false, float deepAK8CutValue= 0.6)
{
  isSignal = isSig;
  initGlobals();	
//TString TTbarFile = "/eos/cms/store/user/gbakas/ttbar/topTagger/April19/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_Copy.root"	
  gStyle->SetOptStat(0);
 
  
  int fileSize = listOfFiles.size();
  TH1F *h_Chi_all[fileSize],*h_Cos_all[fileSize], *hChiRevertBtag[fileSize], *hCosRevertBtag[fileSize]; 
  TFile *inf;
  vector<float> weights(0);
  
 //number of checks: one for Top Tagger and one for DAK8
 const int nChecks = 2;
 //initialize the required histograms 
 TH1F *hCR[3][listOfFiles.size()];
 TH1F *hSR[3][listOfFiles.size()];
 TH1F *hTestBtag[nChecks][listOfFiles.size()];
 TH1F *hTest1Btag[nChecks][listOfFiles.size()];
 
 for(int f=0; f<listOfFiles.size(); f++)
 {
  cout<<"Entering "<<listOfFiles[f]<<endl;
  inf = TFile::Open(eosPath+listOfFiles[f]);	 
  TTree *trIN    = (TTree*)inf->Get("boosted/events");  
  
  float NORM = ((TH1F*)inf->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
  weights.push_back(XSEC[f]/NORM);
  
  int decade(0);
  int NN = trIN->GetEntries();
  
  int nJets,nLeptons;
  float genEvtWeight;
  vector<float> *jetPt(0),*tau3(0),*tau2(0),*tau1(0);
  vector<float> *jetMassSub0(0), *jetMassSub1(0);
  vector<float> *jetMassSoftDrop(0);

  float mva(0);
  vector<float> *jetTtag(0);
  vector<bool> *bit(0);
  float mTTbarParton(0),mJJ(0), yJJ(0), ptJJ(0), yTTbarParton(0), ptTTbarParton(0);
  int  category(0);
  //matching info 
  vector<float> *jetPhi(0), *jetEta(0);
  vector<int> *partonId(0), *partonMatchIdx(0);
  
  vector<float> *partonEta(0), *partonPhi(0), *partonMatchDR(0),  *partonPt(0), *partonE(0), *partonMass(0), *deepAK8(0);
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
  trIN->SetBranchAddress("mJJ"   		  ,&mJJ);
  trIN->SetBranchAddress("yJJ"   		  ,&yJJ);
  trIN->SetBranchAddress("ptJJ"   		  ,&ptJJ);
  trIN->SetBranchAddress("jetBtagSub0"	  ,&jetBtagSub0);
  trIN->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1);
  trIN->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
  trIN->SetBranchAddress("mva"	  		  ,&mva);
  trIN->SetBranchAddress("category"	  	  ,&category);
  trIN->SetBranchAddress("deepAk8"        ,&deepAK8);
  trIN->SetBranchAddress("jetTtagCategory",&jetTtag);
  
  if(isSignal)
  {
	trIN->SetBranchAddress("mTTbarParton"	,&mTTbarParton);
	trIN->SetBranchAddress("yTTbarParton"	,&yTTbarParton);
	trIN->SetBranchAddress("ptTTbarParton"	,&ptTTbarParton);
	trIN->SetBranchAddress("partonPt"	    ,&partonPt);
	trIN->SetBranchAddress("partonEta"		,&partonEta);
	trIN->SetBranchAddress("partonMass"     ,&partonMass);
	trIN->SetBranchAddress("partonMatchDR"  ,&partonMatchDR);
	trIN->SetBranchAddress("partonMatchIdx" ,&partonMatchIdx);
	trIN->SetBranchAddress("partonPhi"      ,&partonPhi);
  }
  bool isJetVar;
  if(varNum > 3)
  {
	isJetVar = true;
  }
  float xParton(0), xReco(0);
  float xJetParton[2], xJetReco[2];
  
  //book the histograms
  //histograms for Signal/QCD in CR 
  hCR[0][f] = new TH1F(TString::Format("hCR_%s_%s", "tTagger",histoNames[f].Data()), TString::Format("hCR_%s_%s","tTagger",histoNames[f].Data()), sizeBins, BND);
  hCR[1][f] = new TH1F(TString::Format("hCR_%s_%s", "deepAK8",histoNames[f].Data()), TString::Format("hCR_%s_%s","deepAK8",histoNames[f].Data()), sizeBins, BND);
  hCR[2][f] = new TH1F(TString::Format("hCR_%s_%s", "oldMva" ,histoNames[f].Data()), TString::Format("hCR_%s_%s","oldMva" ,histoNames[f].Data()), sizeBins, BND);
  
  //histograms for Signal/QCD in SR
  hSR[0][f] = new TH1F(TString::Format("hSR_%s_%s", "tTagger",histoNames[f].Data()), TString::Format("hSR_%s_%s","tTagger",histoNames[f].Data()), sizeBins, BND);
  hSR[1][f] = new TH1F(TString::Format("hSR_%s_%s", "deepAK8",histoNames[f].Data()), TString::Format("hSR_%s_%s","deepAK8",histoNames[f].Data()), sizeBins, BND);
  hSR[2][f] = new TH1F(TString::Format("hSR_%s_%s", "oldMva" ,histoNames[f].Data()), TString::Format("hSR_%s_%s","oldMva" ,histoNames[f].Data()), sizeBins, BND);
  
  hTestBtag[0][f] = new TH1F(TString::Format("hSR_%s_%s", "noBtag_tTagger",histoNames[f].Data()), TString::Format("hSR_%s_%s","noBtag_tTagger",histoNames[f].Data()), sizeBins, BND);
  hTestBtag[1][f] = new TH1F(TString::Format("hSR_%s_%s", "noBtag_deepAK8",histoNames[f].Data()), TString::Format("hSR_%s_%s","noBtag_deepAK8",histoNames[f].Data()), sizeBins, BND);
  hTest1Btag[0][f] = new TH1F(TString::Format("hSR_%s_%s", "1Btag_tTagger",histoNames[f].Data()), TString::Format("hSR_%s_%s","1Btag_tTagger",histoNames[f].Data()), sizeBins, BND);
  hTest1Btag[1][f] = new TH1F(TString::Format("hSR_%s_%s", "1Btag_deepAK8",histoNames[f].Data()), TString::Format("hSR_%s_%s","1Btag_deepAK8",histoNames[f].Data()), sizeBins, BND);
  
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

	if (varNum ==1)
	{
	 if(isSignal) xParton = mTTbarParton;
	 xReco   = mJJ;
	}
	else if (varNum ==2)
	{
	 if(isSignal) xParton = ptTTbarParton;
	 xReco   = ptJJ;
	}
	else if (varNum ==3) 
	{
	 if(isSignal) xParton = yTTbarParton;
	 xReco   = yJJ;
	}
	bool partonCuts, recoCuts, btagging, massCut, revertBtag,tTaggerCut, deepAK8Cut, btag1,oldMva;	
	if (nJets >1)
	{				
		//matching only if we have Signal/
		if(isSignal)
		{
		//----------------------MATCHING------------------------------------------------------
		
		for(int ijet =0; ijet<nJets; ijet++)
		{
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
			
			for(int k=1; k<jetMatchedIndexes->size(); k++)
			{
				if((*jetMatchedDr)[k] < dRmin)
				{
					dRmin = (*jetMatchedDr)[k];
					indexMin = (*jetMatchedIndexes)[k];
				}
			
			}
			
			if(dRmin < 0.4)
			{
				isMatched++;
				jetDr_ = dRmin;
				//RECO MATCHED
				pt_->push_back((*jetPt)[(*partonMatchIdx)[indexMin]]);
				mass_->push_back((*jetMassSoftDrop)[(*partonMatchIdx)[indexMin]]);
				eta_->push_back((*jetEta)[(*partonMatchIdx)[indexMin]]);
				phi_->push_back( (*jetPhi)[(*partonMatchIdx)[indexMin]]);
				jetBtagSub0_->push_back( (*jetBtagSub0)[(*partonMatchIdx)[indexMin]]);
				jetBtagSub1_->push_back( (*jetBtagSub1)[(*partonMatchIdx)[indexMin]]);
				jetTtag_->push_back((*jetTtag)[(*partonMatchIdx)[indexMin]]);
				deepAK8_->push_back((*deepAK8)[(*partonMatchIdx)[indexMin]]);
				//PARTON MATCHED
				partonPt_->push_back( (*partonPt)[indexMin]);
				partonMass_->push_back( (*partonMass)[indexMin]);
				partonPhi_->push_back( (*partonPhi)[indexMin]);
				partonEta_->push_back( (*partonEta)[indexMin]);
			}		  
		 }//----end of if jetMatchedIndexes > 0
		}//----end of for on all jets for mathching
	   
	  //---------------------------------------END OF MATCHING------------------------------------------------------
	  recoCuts   = fabs((*eta_)[0]) < 2.4 && fabs((*eta_)[1]) <2.4 && (*pt_)[0] > 400 && (*pt_)[1] > 400 && nLeptons==0 && mJJ > 1000 && nJets > 1;
	  partonCuts = fabs((*partonEta_)[0]) < 2.4 && fabs((*partonEta_)[1] <2.4) && (*partonPt_)[0] > 400 && (*partonPt_)[1] > 400 && mTTbarParton > 1000;
	  btagging   = ((*jetBtagSub0_)[0] > floatBTag || (*jetBtagSub1_)[0] > floatBTag) && ((*jetBtagSub0_)[1] > floatBTag || (*jetBtagSub1_)[1] > floatBTag);
      massCut    = (*mass_)[0] > 120 && (*mass_)[0] < 220 && (*mass_)[1] > 120 && (*mass_)[1] < 220;
	  oldMva 	 =  mva >0.8;
	  revertBtag = ((*jetBtagSub0_)[0] < floatBTag && (*jetBtagSub1_)[0] < floatBTag) && ((*jetBtagSub0_)[1] < floatBTag && (*jetBtagSub1_)[1] < floatBTag);
	  btag1 	 = (((*jetBtagSub0_)[0] > floatBTag || (*jetBtagSub1_)[0] > floatBTag) && ((*jetBtagSub0_)[1] < floatBTag || (*jetBtagSub1_)[1] < floatBTag)) ||
					(((*jetBtagSub0_)[0] < floatBTag || (*jetBtagSub1_)[0] < floatBTag) && ((*jetBtagSub0_)[1] > floatBTag || (*jetBtagSub1_)[1] > floatBTag));

	  tTaggerCut = (*jetTtag_)[0] > selMvaCut && (*jetTtag_)[1] > selMvaCut;
	  deepAK8Cut = (*deepAK8_)[0] > deepAK8CutValue && (*deepAK8_)[1] > deepAK8CutValue;
	  
	  if(isMatched > 1)
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
	  

	 
    }//----end of isSignal so that we do this only when we deal with signal	   
	
	else //we are in QCD samples
	{
		
	  recoCuts   = fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1]) <2.4 && (*jetPt)[0] > 400 && (*jetPt)[1] > 400 &&  mJJ > 1000; //nLeptons==0 &&
	  btagging   = ((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag) && ((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag);
      massCut    = (*jetMassSoftDrop)[0] > 120 && (*jetMassSoftDrop)[0] < 220 && (*jetMassSoftDrop)[1] > 120 && (*jetMassSoftDrop)[1] < 220;
	  revertBtag = ((*jetBtagSub0)[0] < floatBTag && (*jetBtagSub1)[0] < floatBTag) && ((*jetBtagSub0)[1] < floatBTag && (*jetBtagSub1)[1] < floatBTag);
	  tTaggerCut = (*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut;
	  oldMva 	 =  mva >0.8;
	  deepAK8Cut = (*deepAK8)[0] > deepAK8CutValue && (*deepAK8)[1] > deepAK8CutValue;
	  btag1 	 = (((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag) && ((*jetBtagSub0)[1] < floatBTag || (*jetBtagSub1)[1] < floatBTag)) ||
					(((*jetBtagSub0)[0] < floatBTag || (*jetBtagSub1)[0] < floatBTag) && ((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag));
	  
	  if (varNum ==4)
	  {
		xJetReco[0] = (*jetPt)[0];
		xJetReco[1] = (*jetPt)[1];
	  }
	  else if (varNum ==5)
	  {					
	    xJetReco[0] = (*jetEta)[0];
		xJetReco[1] = (*jetEta)[1];
	  }
	  
	  
	}//---end of else of isSignal
	/*
	 cout<<"---------------------"<<endl;
	 cout<<(*jetEta)[0]<<endl;
	 cout<<(*jetEta)[1]<<endl;
	 cout<<(*jetPt)[0]<<endl;
	 cout<<(*jetPt)[1]<<endl;
	 cout<<mJJ<<endl;
	 cout<<"jet 0: "<<(*jetBtagSub0)[0]<<"  "<<(*jetBtagSub1)[0]<<endl;
	 cout<<"jet 1: "<<(*jetBtagSub0)[1]<<"  "<<(*jetBtagSub1)[1]<<endl;
	 cout<<(*jetMassSoftDrop)[1]<<"   "<<(*jetMassSoftDrop)[0]<<endl;
	 cout<<(*jetTtag)[0]<<"   "<<(*jetTtag)[0]<<endl;
	 cout<<(*deepAK8)[0]<<"   "<<(*deepAK8)[0]<<endl;
	 cout<<"booleans"<<endl;
	 cout<<"SR tTagger: "<<(recoCuts && btagging && massCut && tTaggerCut)<<endl;
	 cout<<"SR deepAK8: "<<(recoCuts && btagging && massCut && deepAK8Cut)<<endl;
	 cout<<"CR tTagger: "<<(recoCuts && revertBtag && massCut && tTaggerCut)<<endl;
	 cout<<"CR deepAK8: "<<(recoCuts && revertBtag && massCut && deepAK8Cut)<<endl;
	 cout<<"1 btag tTagger: "<<(recoCuts && btag1 && massCut && tTaggerCut)<<endl;
	 cout<<"1 btag deepAK8: "<<(recoCuts && btag1 && massCut && deepAK8Cut)<<endl;
	 cout<<"No b tag requirement tTagger: "<<(recoCuts && massCut && tTaggerCut)<<endl;
	 cout<<"No b tag requirement deepAK8: "<<(recoCuts && massCut && deepAK8Cut)<<endl;
	*/
	 //Signal Region with tTagger
	 if(recoCuts && btagging && massCut && tTaggerCut)
	 {
		if(!isJetVar) hSR[0][f]->Fill(xReco);
	    else 
		{
			hSR[0][f]->Fill(xJetReco[0]);
	        hSR[0][f]->Fill(xJetReco[1]);
		}
	  }
	  //Control Region with tTagger
	 if(recoCuts && revertBtag && massCut && tTaggerCut)
	 {
		if(!isJetVar) hCR[0][f]->Fill(xReco);
	    else 
		{
			hCR[0][f]->Fill(xJetReco[0]);
	        hCR[0][f]->Fill(xJetReco[1]);
		}
	 }
	  
	 //Signal Region with DeepAK8
	 if(recoCuts && btagging && massCut && deepAK8Cut)
	 {
		if(!isJetVar) hSR[1][f]->Fill(xReco);
	    else 
		{
			hSR[1][f]->Fill(xJetReco[0]);
	        hSR[1][f]->Fill(xJetReco[1]);
		}
	 }
	 //Control Region with DeepAK8
	 if(recoCuts && revertBtag && massCut && deepAK8Cut)
	 {
		if(!isJetVar) hCR[1][f]->Fill(xReco);
	    else 
		{
			hCR[1][f]->Fill(xJetReco[0]);
	        hCR[1][f]->Fill(xJetReco[1]);
		}
	 }
	 //Signal Region with oldMva
	 if(recoCuts && btagging && massCut && oldMva)
	 {
		if(!isJetVar) hSR[2][f]->Fill(xReco);
	    else 
		{
			hSR[2][f]->Fill(xJetReco[0]);
	        hSR[2][f]->Fill(xJetReco[1]);
		}
	 }
	 //Control Region with oldMva
	 if(recoCuts && revertBtag && massCut && oldMva)
	 {
		if(!isJetVar) hCR[2][f]->Fill(xReco);
	    else 
		{
			hCR[2][f]->Fill(xJetReco[0]);
	        hCR[2][f]->Fill(xJetReco[1]);
		}
	 }
	 
	 //Control Region with DeepAK8
	 //this is just for testing
	 //I Want to see how the tagger reacts even without the btagging cut
	 if(recoCuts && massCut && tTaggerCut)
	 {
		if(!isJetVar) hTestBtag[0][f]->Fill(xReco);
	    else 
		{
			hTestBtag[0][f]->Fill(xJetReco[0]);
	        hTestBtag[0][f]->Fill(xJetReco[1]);
		}
	 }
	 if(recoCuts && massCut && deepAK8Cut)
	 {
		if(!isJetVar) hTestBtag[1][f]->Fill(xReco);
	    else 
		{
			hTestBtag[1][f]->Fill(xJetReco[0]);
	        hTestBtag[1][f]->Fill(xJetReco[1]);
		}
	 }
	 //1 btag region
	 if(recoCuts && massCut && tTaggerCut && btag1)
	 {
		if(!isJetVar) hTest1Btag[0][f]->Fill(xReco);
	    else 
		{
			hTest1Btag[0][f]->Fill(xJetReco[0]);
	        hTest1Btag[0][f]->Fill(xJetReco[1]);
		}
	 }
	 if(recoCuts && massCut && deepAK8Cut && btag1)
	 {
		if(!isJetVar) hTest1Btag[1][f]->Fill(xReco);
	    else 
		{
			hTest1Btag[1][f]->Fill(xJetReco[0]);
	        hTest1Btag[1][f]->Fill(xJetReco[1]);
		}
	 }
	
					

	}//----end of nJets
  }	//---end of event loop
  }//----end of fileSize loop 
  
  TH1F *hCR_Clone[2][listOfFiles.size()];
  TH1F *hSR_Clone[2][listOfFiles.size()];
   for(int i=0; i<2; i++)
  {
	//for every slice
	
	for(int j=0; j<listOfFiles.size(); j++)
    {
	  if(i==0)
	  {
		  hCR_Clone[i][j]=(TH1F*)hCR[i][j]->Clone(TString::Format("hCR_%s_%s_Clone","tTagger",histoNames[j].Data()));	
		  hSR_Clone[i][j]=(TH1F*)hSR[i][j]->Clone(TString::Format("hSR_%s_%s_Clone","tTagger",histoNames[j].Data()));	
	  }
	  else 
	  {
		  hCR_Clone[i][j]=(TH1F*)hCR[i][j]->Clone(TString::Format("hCR_%s_%s_Clone","deepAK8",histoNames[j].Data()));	
		  hSR_Clone[i][j]=(TH1F*)hSR[i][j]->Clone(TString::Format("hSR_%s_%s_Clone","deepAK8",histoNames[j].Data()));	
	  }
	  
      hCR_Clone[i][j]->Scale(weights[j]*LUMI); //this is 0 btagged (CR)
      hSR_Clone[i][j]->Scale(weights[j]*LUMI); //this is 2 btagged (SR)
	
	  hCR[i][j]->Scale(weights[j]);
	  hSR[i][j]->Scale(weights[j]);
	  
	  hTest1Btag[i][j]->Scale(weights[j]*LUMI); //this is 1 btag
	  
	  hTestBtag[i][j]->Scale(weights[j]*LUMI); //this is no b tag
	 
    }
    
    for(int j=1; j<listOfFiles.size(); j++)
	{
	  //Add them to get the whole phase space
	  hCR[i][0]->Add(hCR[i][j]);
      hSR[i][0]->Add(hSR[i][j]);
      hCR_Clone[i][0]->Add(hCR_Clone[i][j]);
      hSR_Clone[i][0]->Add(hSR_Clone[i][j]);
      hTestBtag[i][0]->Add(hTestBtag[i][j]);	  
      hTest1Btag[i][0]->Add(hTest1Btag[i][j]);	  
	}
  }
  //for old mva
  for(int i=2; i<3; i++)
  {
	  for(int j =0;j <listOfFiles.size(); j++)
	  {
		  hSR[i][j]->Scale(weights[j]);
		  hCR[i][j]->Scale(weights[j]);
	  }
	  for(int j =1; j<listOfFiles.size(); j++)
	  {
		  hCR[i][0]->Add(hCR[i][j]);
		  hSR[i][0]->Add(hSR[i][j]);
	  }
  }	  
  
   
  TCanvas *test1 = new TCanvas("test1", "test1", 700, 600);  
  TLegend *tLeg = new TLegend(0.6,0.7,0.8,0.9);
  hCR_Clone[1][0]->SetLineColor(kRed);
  hTestBtag[1][0]->SetLineColor(kMagenta);
  hSR_Clone[1][0]->SetLineColor(kBlue);
  hTest1Btag[1][0]->SetLineColor(kGreen);
  tLeg->AddEntry(hSR_Clone[1][0], "Signal Region (2btag)", "l");
  tLeg->AddEntry(hCR_Clone[1][0], "Control Region (0btag)", "l");
  tLeg->AddEntry(hTestBtag[1][0], "Region No b tagging requirement", "l");
  tLeg->AddEntry(hTest1Btag[1][0], "Region 1b tag requirement", "l");
  hSR_Clone[1][0]->Draw();
  hCR_Clone[1][0]->Draw("same");
  hTestBtag[1][0]->Draw("same");
  hTest1Btag[1][0]->Draw("same");
  tLeg->Draw();
  
  test1->Write("testCan_deepAK8");
  
  TCanvas *test2 = new TCanvas("test2 tTagger", "test2 tTagger", 700, 600);  
  TLegend *tLeg2 = new TLegend(0.6,0.7,0.8,0.9);
  hCR_Clone[0][0]->SetLineColor(kRed);
  hTestBtag[0][0]->SetLineColor(kMagenta);
  hSR_Clone[0][0]->SetLineColor(kBlue);
  hTest1Btag[0][0]->SetLineColor(kGreen);
  tLeg2->AddEntry(hSR_Clone[0][0], "Signal Region (2btag)", "l");
  tLeg2->AddEntry(hCR_Clone[0][0], "Control Region (0btag)", "l");
  tLeg2->AddEntry(hTestBtag[0][0], "Region No b tagging requirement", "l");
  tLeg2->AddEntry(hTest1Btag[0][0], "Region 1btag requirement", "l");
  hSR_Clone[0][0]->Draw();
  hCR_Clone[0][0]->Draw("same");
  hTestBtag[0][0]->Draw("same");
  hTest1Btag[0][0]->Draw("same");
  tLeg2->Draw();
  
  
  cout<<"no btag tTagger:"<<hTestBtag[0][0]->GetEntries()<<endl;
  cout<<"no btag deepAk8:"<<hTestBtag[1][0]->GetEntries()<<endl;
  cout<<"0btag tTagger:"<<hCR_Clone[0][0]->GetEntries()<<endl;
  cout<<"0btag deepAk8:"<<hCR_Clone[1][0]->GetEntries()<<endl;
  cout<<"1btag tTagger:"<<hTest1Btag[0][0]->GetEntries()<<endl;
  cout<<"1btag deepAk8:"<<hTest1Btag[1][0]->GetEntries()<<endl;
  cout<<"2btag tTagger:"<<hSR_Clone[0][0]->GetEntries()<<endl;
  cout<<"2btag deepAk8:"<<hSR_Clone[1][0]->GetEntries()<<endl;
  
  
   
  hCR[0][0]->SetTitle("CR tTagger");
  hCR[1][0]->SetTitle("CR deepAK8");
  hSR[0][0]->SetTitle("SR tTagger");
  hSR[1][0]->SetTitle("SR deepAK8");

  hCR[0][0]->Scale(1./hCR[0][0]->Integral());
  hCR[1][0]->Scale(1./hCR[1][0]->Integral());
  hSR[0][0]->Scale(1./hSR[0][0]->Integral());
  hSR[1][0]->Scale(1./hSR[1][0]->Integral());
  hCR[2][0]->Scale(1./hCR[2][0]->Integral());
  hSR[2][0]->Scale(1./hSR[2][0]->Integral());
  
  hCR[0][0]->SetLineColor(kBlue);
  hCR[1][0]->SetLineColor(kBlue);
  hCR[2][0]->SetLineColor(kBlue);
  hSR[0][0]->SetLineColor(kRed);
  hSR[1][0]->SetLineColor(kRed);
  hSR[2][0]->SetLineColor(kRed);
  
  for(int i=0; i<sizeof(hCR)/sizeof(*hCR); i++)
  {
    if(varNum ==1 || varNum ==2 || varNum == 4)
	{
		hCR[i][0]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameReco.Data()));
		hSR[i][0]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameReco.Data()));
		hCR_Clone[i][0]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameReco.Data()));
		hSR_Clone[i][0]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameReco.Data()));
	}
	else
	{
		hCR[i][0]->GetXaxis()->SetTitle(TString::Format("%s", varNameReco.Data()));
		hSR[i][0]->GetXaxis()->SetTitle(TString::Format("%s", varNameReco.Data()));
		hCR_Clone[i][0]->GetXaxis()->SetTitle(TString::Format("%s", varNameReco.Data()));
		hSR_Clone[i][0]->GetXaxis()->SetTitle(TString::Format("%s", varNameReco.Data()));
	}
  }
  TCanvas *can_CR_tTagger = new TCanvas("can_CR tTagger", "can_CR tTagger", 700, 600);
  TLegend *leg_tTagger = new TLegend(0.6,0.7,0.8,0.9);
  leg_tTagger->AddEntry(hCR[0][0],"Control Region", "l");
  leg_tTagger->AddEntry(hSR[0][0],"Signal Region", "l");
  hCR[0][0]->Draw();
  hSR[0][0]->Draw("same");
  leg_tTagger->Draw();
  
  TCanvas *can_CR_deepAK8 = new TCanvas("can_CR deepAK8", "can_CR deepAK8", 700, 600);
  hCR[1][0]->Draw();
  
  hSR[1][0]->Draw("same");
  TLegend *leg_deepAK8 = new TLegend(0.6,0.7,0.8,0.9);
  leg_deepAK8->AddEntry(hCR[1][0],"Control Region", "l");
  leg_deepAK8->AddEntry(hSR[1][0],"Signal Region", "l");
  leg_deepAK8->Draw();
  
  
  TCanvas *can_CR_oldMva = new TCanvas("can_CR oldMva", "can_CR oldMva", 700, 600);
  hCR[2][0]->Draw();
  
  hSR[2][0]->Draw("same");
  TLegend *leg_oldMva = new TLegend(0.6,0.7,0.8,0.9);
  leg_oldMva->AddEntry(hCR[2][0],"Control Region", "l");
  leg_oldMva->AddEntry(hSR[2][0],"Signal Region", "l");
  leg_oldMva->Draw();
  
  
  TFile *outFile;
  if(isSignal) outFile = new TFile("SignalOutput_AllRegions.root", "UPDATE");
  else outFile = new TFile("BkgOutput_AllRegions_1.root", "UPDATE");
  
  outFile->cd();
  if(saveTTagger) 
  {
	  hSR[0][0]->Write(TString::Format("SR_tTagger_%s", varNameReco.Data()));
	  hCR[0][0]->Write(TString::Format("CR_tTagger_%s", varNameReco.Data()));
	  hSR_Clone[0][0]->Write(TString::Format("SR_tTagger_%s_expYield", varNameReco.Data()));
	  hCR_Clone[0][0]->Write(TString::Format("CR_tTagger_%s_expYield", varNameReco.Data()));
  }
  hSR[1][0]->Write(TString::Format("SR_deepAK8_%0.1f_%s",deepAK8CutValue,varNameReco.Data()));
  hCR[1][0]->Write(TString::Format("CR_deepAK8_%0.1f_%s",deepAK8CutValue, varNameReco.Data()));
  hSR_Clone[1][0]->Write(TString::Format("SR_deepAK8_%0.1f_%s_expYield",deepAK8CutValue, varNameReco.Data()));
  hCR_Clone[1][0]->Write(TString::Format("CR_deepAK8_%0.1f_%s_expYield",deepAK8CutValue, varNameReco.Data()));
  
  test1->Write("testCan_deepAK8");
  test2->Write("testCan_tTagger");
  
  
  if(!isSignal)
  {
	auto c1 = new TCanvas("Control Region tTagger", "Control Region tTagger", 900,600);
	auto rp = new TRatioPlot(hSR[0][0],hCR[0][0]);
	c1->SetTicks(0,1);
	rp->Draw();
	//leg_deepAK8->Draw();
	c1->Update();
  
	auto c2 = new TCanvas("Control Region deepAK8", "Control Region deepAK8", 900,600);
	auto rp2 = new TRatioPlot(hSR[1][0],hCR[1][0]);
	c2->SetTicks(0,1);
	rp2->Draw();
	//leg_deepAK8->Draw();
	c2->Update();
	
	auto c3 = new TCanvas("Control Region old mva", "Control Region oldMva", 900,600);
	auto rp3 = new TRatioPlot(hSR[2][0],hCR[2][0]);
	c3->SetTicks(0,1);
	rp3->Draw();
	//leg_deepAK8->Draw();
	c3->Update();
	
	if(saveTTagger) c1->Write(TString::Format("TRatioPlot_tTagger_%s", varNameReco.Data()));
    c2->Write(TString::Format("TRatioPlot_deepAK8_%0.1f_%s",deepAK8CutValue,varNameReco.Data()));
  }
  
  outFile->Close();
  listOfFiles.clear();
  XSEC.clear();
  histoNames.clear();
  
}

  