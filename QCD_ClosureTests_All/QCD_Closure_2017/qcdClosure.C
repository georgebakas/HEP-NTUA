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
float LUMI = 41530;
TString eosPath;
bool isSignal;
float deepCSVFloat = 0.4941;

void initFileNames()
{
  if(isSignal)
  {
	eosPath = "/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2017/Signal/";  
	listOfFiles.push_back("TT_Mtt-700to1000_TuneCP5_13TeV-powheg-pythia8.root");
	listOfFiles.push_back("TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root");
  }
  else
  {
	eosPath = "/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2017/Bkg/";
	listOfFiles.push_back("QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8.root");
	listOfFiles.push_back("QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8.root");
	listOfFiles.push_back("QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8.root");
	listOfFiles.push_back("QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8.root");
	listOfFiles.push_back("QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8.root");
	listOfFiles.push_back("QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8.root");
	
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
 
void qcdClosure(TString varName,TString varNameReco, int varNum,  int sizeBins,float BND[] ,  bool saveTTagger= true, float selMvaCut=0.2, float floatBTag = 0.8838, bool 				isSig = false, float deepAK8CutValue= 0.6, bool isDeepCSV= true)
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
 TH1F *h1Btag[3][listOfFiles.size()];
 
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
  vector<bool> *bit = new vector<bool>;
  float mTTbarParton(0),mJJ(0), yJJ(0), ptJJ(0), yTTbarParton(0), ptTTbarParton(0);
  int  category(0);
  //matching info 
  vector<float> *jetPhi(0), *jetEta(0), *jetY(0);
  vector<int> *partonId(0), *partonMatchIdx(0);
  
  vector<float> *partonEta(0), *partonPhi(0), *partonMatchDR(0),  *partonPt(0), *partonE(0), *partonMass(0), *deepAK8(0);
  std::vector<int> *addedIndexes = new std::vector<int>(0);
  std::vector<float> *jetBtagSub0(0), *jetBtagSub1(0);
  std::vector<float> *jetBtagSub0DCSVbb(0), *jetBtagSub1DCSVbb(0);
  std::vector<float> *jetBtagSub0DCSVbbb(0), *jetBtagSub1DCSVbbb(0);
          
  //------- input tree --------------
  trIN->SetBranchAddress("nJets"          ,&nJets);
  trIN->SetBranchAddress("nLeptons"       ,&nLeptons);
  trIN->SetBranchAddress("jetPt"          ,&jetPt);
  trIN->SetBranchAddress("jetEta"         ,&jetEta);
  trIN->SetBranchAddress("jetY"           ,&jetY);
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
  trIN->SetBranchAddress("deepAK8Tagger"  ,&deepAK8);
  trIN->SetBranchAddress("jetTtagCategory",&jetTtag);
  
  //deepCSV
  trIN->SetBranchAddress("jetBtagSub0DCSVbb" ,&jetBtagSub0DCSVbb);
  trIN->SetBranchAddress("jetBtagSub1DCSVbb" ,&jetBtagSub1DCSVbb);
  trIN->SetBranchAddress("jetBtagSub0DCSVbbb",&jetBtagSub0DCSVbbb);
  trIN->SetBranchAddress("jetBtagSub1DCSVbbb",&jetBtagSub1DCSVbbb);
  
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
  
  h1Btag[0][f] = new TH1F(TString::Format("h%s_%s", "1Btag_tTagger",histoNames[f].Data()), TString::Format("h%s_%s","1Btag_tTagger",histoNames[f].Data()), sizeBins, BND);
  h1Btag[1][f] = new TH1F(TString::Format("h%s_%s", "1Btag_deepAK8",histoNames[f].Data()), TString::Format("h%s_%s","1Btag_deepAK8",histoNames[f].Data()), sizeBins, BND);
  h1Btag[2][f] = new TH1F(TString::Format("h%s_%s", "1Btag_oldMva",histoNames[f].Data()), TString::Format("h%s_%s","1Btag_oldMva",histoNames[f].Data()), sizeBins, BND);
  
  //for matching
  std::vector<int> *jetMatchedIndexes = new std::vector<int>(0);
  std::vector<float> *jetMatchedDr = new std::vector<float>(0);
  std::vector<float> *eta_ = new std::vector<float>(0);
  std::vector<float> *y_ = new std::vector<float>(0);
  std::vector<float> *phi_ = new std::vector<float>(0);
  std::vector<float> *mass_ = new std::vector<float>(0);
  std::vector<float> *pt_ = new std::vector<float>(0);
  std::vector<float> *deepAK8_ = new std::vector<float>(0);
  std::vector<float> *jetTtag_ = new std::vector<float>(0);
  std::vector<float> *jetBtagSub0_ = new std::vector<float>(0);
  std::vector<float> *jetBtagSub1_ = new std::vector<float>(0);

  std::vector<float> *partonPt_ = new std::vector<float>(0);
  std::vector<float> *partonEta_ = new std::vector<float>(0);
  std::vector<float> *partonY_ = new std::vector<float>(0);
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
	
	int isMatched =0;
	eta_->clear();
	y_->clear();
    mass_->clear();
    pt_->clear();
    phi_->clear();
    jetBtagSub0_->clear();
    jetBtagSub1_->clear();
    jetTtag_->clear();
    deepAK8_->clear();

    partonPt_->clear();
    partonY_->clear();
    partonMass_->clear();
    partonPhi_->clear();
    partonEta_->clear();
	
	jetBtagSub0DCSVbb_->clear();
    jetBtagSub1DCSVbb_->clear();
    jetBtagSub0DCSVbbb_->clear();
    jetBtagSub1DCSVbbb_->clear();

	
	bool partonCuts, recoCuts, massCut, tTaggerCut, deepAK8Cut, oldMva;
    bool deepCSV, btag1DeepCSV, revertBtagDeepCSV;	
	bool CSVv2Cut, revertBtagCSVv2, btag1CSVv2; 
	bool btagCut, revertBtag, btag1;
	
	if (nJets >1)
	{		
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
				y_->push_back((*jetY)[(*partonMatchIdx)[indexMin]]);
				phi_->push_back( (*jetPhi)[(*partonMatchIdx)[indexMin]]);
				jetBtagSub0_->push_back( (*jetBtagSub0)[(*partonMatchIdx)[indexMin]]);
				jetBtagSub1_->push_back( (*jetBtagSub1)[(*partonMatchIdx)[indexMin]]);
				jetTtag_->push_back((*jetTtag)[(*partonMatchIdx)[indexMin]]);
				deepAK8_->push_back((*deepAK8)[(*partonMatchIdx)[indexMin]]);
				
				jetBtagSub0DCSVbb_->push_back((*jetBtagSub0DCSVbb)[(*partonMatchIdx)[indexMin]]);
				jetBtagSub1DCSVbb_->push_back((*jetBtagSub1DCSVbb)[(*partonMatchIdx)[indexMin]]);
				jetBtagSub0DCSVbbb_->push_back((*jetBtagSub0DCSVbbb)[(*partonMatchIdx)[indexMin]]);
				jetBtagSub1DCSVbbb_->push_back((*jetBtagSub1DCSVbbb)[(*partonMatchIdx)[indexMin]]);
	
				//PARTON MATCHED
				partonPt_->push_back( (*partonPt)[indexMin]);
				partonMass_->push_back( (*partonMass)[indexMin]);
				partonPhi_->push_back( (*partonPhi)[indexMin]);
				partonEta_->push_back( (*partonEta)[indexMin]);
				//here misssing partonY
			}		  
		 }//----end of if jetMatchedIndexes > 0
		}//----end of for on all jets for mathching
	   
	  //---------------------------------------END OF MATCHING------------------------------------------------------
	  float dCSVScoreSub0[2], dCSVScoreSub1[2];
	  dCSVScoreSub0[0] = (*jetBtagSub0DCSVbb_)[0] + (*jetBtagSub0DCSVbbb_)[0];
	  dCSVScoreSub0[1] = (*jetBtagSub0DCSVbb_)[1] + (*jetBtagSub0DCSVbbb_)[1];
	  dCSVScoreSub1[0] = (*jetBtagSub1DCSVbb_)[0] + (*jetBtagSub1DCSVbbb_)[0];
	  dCSVScoreSub1[1] = (*jetBtagSub1DCSVbb_)[1] + (*jetBtagSub1DCSVbbb_)[1];
	  
	  recoCuts   = fabs((*eta_)[0]) < 2.4 && fabs((*eta_)[1]) <2.4 && (*pt_)[0] > 400 && (*pt_)[1] > 400 && nLeptons==0 && mJJ > 1000 && nJets > 1 && (*bit)[5];
	  partonCuts = fabs((*partonEta_)[0]) < 2.4 && fabs((*partonEta_)[1] <2.4) && (*partonPt_)[0] > 400 && (*partonPt_)[1] > 400 && mTTbarParton > 1000;
      massCut    = (*mass_)[0] > 120 && (*mass_)[0] < 220 && (*mass_)[1] > 120 && (*mass_)[1] < 220;
	  oldMva 	 =  mva >0.8;
	  tTaggerCut = (*jetTtag_)[0] > selMvaCut && (*jetTtag_)[1] > selMvaCut;
	  deepAK8Cut = (*deepAK8_)[0] > deepAK8CutValue && (*deepAK8_)[1] > deepAK8CutValue;
	  //2 btag category with csvv2 and deepCSV
	  CSVv2Cut   = ((*jetBtagSub0_)[0] > floatBTag || (*jetBtagSub1_)[0] > floatBTag) && ((*jetBtagSub0_)[1] > floatBTag || (*jetBtagSub1_)[1] > floatBTag);

	  deepCSV    = (((*jetBtagSub0DCSVbb_)[0] + (*jetBtagSub0DCSVbbb_)[0])> deepCSVFloat || ((*jetBtagSub1DCSVbb_)[0] + (*jetBtagSub1DCSVbbb_)[0])> deepCSVFloat) && 
					 (((*jetBtagSub0DCSVbb_)[1] + (*jetBtagSub0DCSVbbb_)[1])> deepCSVFloat || ((*jetBtagSub1DCSVbb_)[1] + (*jetBtagSub1DCSVbbb_)[1])> deepCSVFloat);
	  //1 btag category with csvv2 and deepCSV			 
	  btag1CSVv2 = (((*jetBtagSub0_)[0] > floatBTag || (*jetBtagSub1_)[0] > floatBTag) && ((*jetBtagSub0_)[1] < floatBTag || (*jetBtagSub1_)[1] < floatBTag)) ||
					(((*jetBtagSub0_)[0] < floatBTag || (*jetBtagSub1_)[0] < floatBTag) && ((*jetBtagSub0_)[1] > floatBTag || (*jetBtagSub1_)[1] > floatBTag));
					
	  btag1DeepCSV	= ((dCSVScoreSub0[0] > deepCSVFloat || dCSVScoreSub1[0] > deepCSVFloat) && (dCSVScoreSub0[1] < deepCSVFloat && dCSVScoreSub1[1] < deepCSVFloat)) ||
					  ((dCSVScoreSub0[0] < deepCSVFloat && dCSVScoreSub1[0] < deepCSVFloat) && (dCSVScoreSub0[1] > deepCSVFloat || dCSVScoreSub1[1] > deepCSVFloat));
	 
	 //0 btag category with csvv2 and deepCSV
      revertBtagCSVv2 = ((*jetBtagSub0)[0] < floatBTag && (*jetBtagSub1)[0] < floatBTag) && ((*jetBtagSub0)[1] < floatBTag && (*jetBtagSub1)[1] < floatBTag);
      revertBtagDeepCSV = (dCSVScoreSub0[0] < deepCSVFloat &&  dCSVScoreSub1[0] < deepCSVFloat) && (dCSVScoreSub0[1] < deepCSVFloat && dCSVScoreSub1[1] < deepCSVFloat);
	  
	  
	  
	  
	  if(isMatched > 1)
	  {
		if (varNum ==4)
		{
		   xParton = (*partonPt_)[0];
		   xReco   = (*pt_)[0];
		}
		else if (varNum ==5)
		{					
		   xParton = (*partonPt_)[1];
		   xReco   = (*pt_)[1];
		}
		else if(varNum ==6)
			xReco = fabs((*y_)[0]);
		else if (varNum == 7)
			xReco = fabs((*y_)[1]);
		  
	  }
	  

	 
    }//----end of isSignal so that we do this only when we deal with signal	   
	
	else //we are in QCD samples
	{
      float dCSVScoreSub0[2], dCSVScoreSub1[2];
	  dCSVScoreSub0[0] = (*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0];
	  dCSVScoreSub0[1] = (*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1];
	  dCSVScoreSub1[0] = (*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0];
	  dCSVScoreSub1[1] = (*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1];
	  
	  recoCuts   = fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1]) <2.4 && (*jetPt)[0] > 400 && (*jetPt)[1] > 400 &&  mJJ > 1000 && (*bit)[5] && nLeptons==0;
      massCut    = (*jetMassSoftDrop)[0] > 120 && (*jetMassSoftDrop)[0] < 220 && (*jetMassSoftDrop)[1] > 120 && (*jetMassSoftDrop)[1] < 220;
	  tTaggerCut = (*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut;
	  oldMva 	 =  mva >0.8;
	  deepAK8Cut = (*deepAK8)[0] > deepAK8CutValue && (*deepAK8)[1] > deepAK8CutValue;
	  //2 btag category with csvv2 and deepCSV
	  CSVv2Cut   = ((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag) && ((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag);

	  deepCSV    = (((*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0])> deepCSVFloat) && 
					 (((*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1])> deepCSVFloat);
	  //1 btag category with csvv2 and deepCSV			 
	  btag1CSVv2 = (((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag) && ((*jetBtagSub0)[1] < floatBTag || (*jetBtagSub1)[1] < floatBTag)) ||
					(((*jetBtagSub0)[0] < floatBTag || (*jetBtagSub1)[0] < floatBTag) && ((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag));
					
	 btag1DeepCSV	= ((dCSVScoreSub0[0] > deepCSVFloat || dCSVScoreSub1[0] > deepCSVFloat) && (dCSVScoreSub0[1] < deepCSVFloat && dCSVScoreSub1[1] < deepCSVFloat)) ||
					  ((dCSVScoreSub0[0] < deepCSVFloat && dCSVScoreSub1[0] < deepCSVFloat) && (dCSVScoreSub0[1] > deepCSVFloat || dCSVScoreSub1[1] > deepCSVFloat));
	 
	  //0 btag category with csvv2 and deepCSV
      revertBtagCSVv2 = ((*jetBtagSub0)[0] < floatBTag && (*jetBtagSub1)[0] < floatBTag) && ((*jetBtagSub0)[1] < floatBTag && (*jetBtagSub1)[1] < floatBTag);
      revertBtagDeepCSV = (dCSVScoreSub0[0] < deepCSVFloat &&  dCSVScoreSub1[0] < deepCSVFloat) && (dCSVScoreSub0[1] < deepCSVFloat && dCSVScoreSub1[1] < deepCSVFloat);
	
	  
	  if (varNum ==4)
	  {
		xParton = (*jetPt)[0];
		xReco   = (*jetPt)[0];
   	  }
	  else if (varNum ==5)
	  {					
	    xParton = (*jetPt)[1];
	    xReco   = (*jetPt)[1];
	  }
	  else if(varNum ==6)
		xReco = fabs((*jetY)[0]);
	  else if (varNum == 7)
		xReco = fabs((*jetY)[1]);
	  
	  
	}//---end of else of isSignal
	
	if(isDeepCSV) 
    {
	  btagCut = deepCSV;
	  revertBtag = revertBtagDeepCSV;
	  btag1 = btag1DeepCSV;
	}
	else 
	{
	  btagCut = CSVv2Cut;
	  revertBtag = revertBtagCSVv2;
	  btag1 = btag1CSVv2;
    }
	
	 //Signal Region with tTagger
	 if(recoCuts && btagCut && massCut && tTaggerCut)
		hSR[0][f]->Fill(xReco,genEvtWeight);
	  //Control Region with tTagger
	 if(recoCuts && revertBtag && massCut && tTaggerCut)
	    hCR[0][f]->Fill(xReco,genEvtWeight);
	  
	 //Signal Region with DeepAK8
	 if(recoCuts && btagCut && massCut && deepAK8Cut)
		hSR[1][f]->Fill(xReco,genEvtWeight);

	 //Control Region with DeepAK8
	 if(recoCuts && revertBtag && massCut && deepAK8Cut)
		hCR[1][f]->Fill(xReco,genEvtWeight);


	 //Signal Region with oldMva
	 if(recoCuts && btagCut && massCut && oldMva)
		hSR[2][f]->Fill(xReco,genEvtWeight);
	    
	 //Control Region with oldMva
	 if(recoCuts && revertBtag && massCut && oldMva)
		hCR[2][f]->Fill(xReco,genEvtWeight);
	    
	 
	 
	 
	 //1 btag region with tTagger
	 if(recoCuts && massCut && tTaggerCut && btag1)
		h1Btag[0][f]->Fill(xReco,genEvtWeight);
	    
	 //1 btag region with deepAK8
	 if(recoCuts && massCut && deepAK8Cut && btag1)
		h1Btag[1][f]->Fill(xReco,genEvtWeight);
	    
	 //1 btag region with oldMva
	 if(recoCuts && massCut && oldMva && btag1)
		 h1Btag[2][f]->Fill(xReco,genEvtWeight);
	    
	 
	

	}//----end of nJets
  }	//---end of event loop
  }//----end of fileSize loop 
  
  TH1F *hCR_Clone[3][listOfFiles.size()];
  TH1F *hSR_Clone[3][listOfFiles.size()];
   for(int i=0; i<3; i++)
  {
	//for every slice
	
	for(int j=0; j<listOfFiles.size(); j++)
    {
	  if(i==0)
	  {
		  hCR_Clone[i][j]=(TH1F*)hCR[i][j]->Clone(TString::Format("hCR_%s_%s_Clone","tTagger",histoNames[j].Data()));	
		  hSR_Clone[i][j]=(TH1F*)hSR[i][j]->Clone(TString::Format("hSR_%s_%s_Clone","tTagger",histoNames[j].Data()));	
	  }
	  else if(i==1)
	  {
		  hCR_Clone[i][j]=(TH1F*)hCR[i][j]->Clone(TString::Format("hCR_%s_%s_Clone","deepAK8",histoNames[j].Data()));	
		  hSR_Clone[i][j]=(TH1F*)hSR[i][j]->Clone(TString::Format("hSR_%s_%s_Clone","deepAK8",histoNames[j].Data()));	
	  }
	  else if(i==2)
	  {
		  hCR_Clone[i][j]=(TH1F*)hCR[i][j]->Clone(TString::Format("hCR_%s_%s_Clone","oldMva",histoNames[j].Data()));	
		  hSR_Clone[i][j]=(TH1F*)hSR[i][j]->Clone(TString::Format("hSR_%s_%s_Clone","oldMva",histoNames[j].Data()));	
	  }
	  
      hCR_Clone[i][j]->Scale(weights[j]*LUMI); //this is 0 btagged (CR)
      hSR_Clone[i][j]->Scale(weights[j]*LUMI); //this is 2 btagged (SR)
	
	  hCR[i][j]->Scale(weights[j]); //this is CR
	  hSR[i][j]->Scale(weights[j]); //this is Signal region
	  h1Btag[i][j]->Scale(weights[j]); //this is 1 btag
	  	 
    }
    
    for(int j=1; j<listOfFiles.size(); j++)
	{
	  //Add them to get the whole phase space
	  hCR[i][0]->Add(hCR[i][j]);
      hSR[i][0]->Add(hSR[i][j]);
      hCR_Clone[i][0]->Add(hCR_Clone[i][j]);
      hSR_Clone[i][0]->Add(hSR_Clone[i][j]);
      h1Btag[i][0]->Add(h1Btag[i][j]);	  
	}
  }
  
  
   
  hCR[0][0]->SetTitle("CR tTagger");
  hCR[1][0]->SetTitle("CR deepAK8");
  hSR[0][0]->SetTitle("SR tTagger");
  hSR[1][0]->SetTitle("SR deepAK8");

 //scale Control region to 1
  hCR[0][0]->Scale(1./hCR[0][0]->Integral(),"width");
  hCR[1][0]->Scale(1./hCR[1][0]->Integral(),"width");
  hCR[2][0]->Scale(1./hCR[2][0]->Integral(),"width");
  
  //scale Signal Region to 1
  hSR[0][0]->Scale(1./hSR[0][0]->Integral(),"width");
  hSR[1][0]->Scale(1./hSR[1][0]->Integral(),"width");
  hSR[2][0]->Scale(1./hSR[2][0]->Integral(),"width");
  
  //scale 1btag region to 1
  h1Btag[0][0]->Scale(1./h1Btag[0][0]->Integral(),"width");
  h1Btag[1][0]->Scale(1./h1Btag[1][0]->Integral(),"width");
  h1Btag[2][0]->Scale(1./h1Btag[2][0]->Integral(),"width");
  
  //set line Colors CR
  hCR[0][0]->SetLineColor(kBlue);
  hCR[1][0]->SetLineColor(kBlue);
  hCR[2][0]->SetLineColor(kBlue);
  //set line Colors SR
  hSR[0][0]->SetLineColor(kRed);
  hSR[1][0]->SetLineColor(kRed);
  hSR[2][0]->SetLineColor(kRed);
  
  //set line Colors 1Btag region
  h1Btag[0][0]->SetLineColor(kMagenta);
  h1Btag[1][0]->SetLineColor(kMagenta);
  h1Btag[2][0]->SetLineColor(kMagenta);
  
  for(int i=0; i<sizeof(hCR)/sizeof(*hCR); i++)
  {
    if(varNum ==1 || varNum ==2 || varNum == 4 || varNum == 5)
	{
		hCR[i][0]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameReco.Data()));
		hSR[i][0]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameReco.Data()));
		h1Btag[i][0]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameReco.Data()));
		hCR_Clone[i][0]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameReco.Data()));
		hSR_Clone[i][0]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameReco.Data()));
	}
	else
	{
		hCR[i][0]->GetXaxis()->SetTitle(TString::Format("%s", varNameReco.Data()));
		hSR[i][0]->GetXaxis()->SetTitle(TString::Format("%s", varNameReco.Data()));
		h1Btag[i][0]->GetXaxis()->SetTitle(TString::Format("%s", varNameReco.Data()));
		hCR_Clone[i][0]->GetXaxis()->SetTitle(TString::Format("%s", varNameReco.Data()));
		hSR_Clone[i][0]->GetXaxis()->SetTitle(TString::Format("%s", varNameReco.Data()));
	}
  }
  /*
  TCanvas *can_CR_tTagger = new TCanvas("can_CR tTagger", "can_CR tTagger", 700, 600);
  TLegend *leg_tTagger = new TLegend(0.6,0.7,0.8,0.9);
  leg_tTagger->AddEntry(hCR[0][0],"Control Region", "l");
  leg_tTagger->AddEntry(hSR[0][0],"Signal Region", "l");
  leg_tTagger->AddEntry(h1Btag[0][0],"1Btag Region", "l");
  hCR[0][0]->Draw();
  hSR[0][0]->Draw("same");
  h1Btag[0][0]->Draw("same");
  leg_tTagger->Draw();
  
  TCanvas *can_CR_deepAK8 = new TCanvas("can_CR deepAK8", "can_CR deepAK8", 700, 600);
  hCR[1][0]->Draw();
  hSR[1][0]->Draw("same");
  h1Btag[1][0]->Draw("same");
  TLegend *leg_deepAK8 = new TLegend(0.6,0.7,0.8,0.9);
  leg_deepAK8->AddEntry(hCR[1][0],"Control Region", "l");
  leg_deepAK8->AddEntry(hSR[1][0],"Signal Region", "l");
  leg_deepAK8->AddEntry(h1Btag[1][0],"1Btag Region", "l");
  leg_deepAK8->Draw();
  
  
  TCanvas *can_CR_oldMva = new TCanvas("can_CR oldMva", "can_CR oldMva", 700, 600);
  hCR[2][0]->Draw();
  hSR[2][0]->Draw("same");
  h1Btag[2][0]->Draw("same");
  TLegend *leg_oldMva = new TLegend(0.6,0.7,0.8,0.9);
  leg_oldMva->AddEntry(hCR[2][0],"Control Region", "l");
  leg_oldMva->AddEntry(hSR[2][0],"Signal Region", "l");
  leg_oldMva->AddEntry(h1Btag[2][0],"1Btag Region", "l");
  leg_oldMva->Draw();
  */
  
  TFile *outFile;
  if(isSignal) 
  {
	  if(isDeepCSV) outFile = new TFile(TString::Format("SignalOutput_AllRegions_%0.2f_deepCSV.root", selMvaCut), "UPDATE");
	  else outFile = new TFile(TString::Format("SignalOutput_AllRegions_%0.2f.root", selMvaCut), "UPDATE");
  }
  else 
  {
	  if(isDeepCSV) outFile = new TFile(TString::Format("BkgOutput_AllRegions_%0.2f_deepCSV.root",selMvaCut), "UPDATE");
	  else outFile = new TFile(TString::Format("BkgOutput_AllRegions_%0.2f.root",selMvaCut), "UPDATE");
  }
  
  outFile->cd();
  if(saveTTagger) 
  {
	  hSR[0][0]->Write(TString::Format("SR_tTagger_%s", varNameReco.Data()));
	  hCR[0][0]->Write(TString::Format("CR_tTagger_%s", varNameReco.Data()));
	  h1Btag[0][0]->Write(TString::Format("h1Btag_tTagger_%s", varNameReco.Data()));
	  hSR_Clone[0][0]->Write(TString::Format("SR_tTagger_%s_expYield", varNameReco.Data()));
	  hCR_Clone[0][0]->Write(TString::Format("CR_tTagger_%s_expYield", varNameReco.Data()));
	  
	  hSR[2][0]->Write(TString::Format("SR_oldMva_%s", varNameReco.Data()));
	  hCR[2][0]->Write(TString::Format("CR_oldMva_%s", varNameReco.Data()));
	  h1Btag[2][0]->Write(TString::Format("h1Btag_oldMva_%s", varNameReco.Data()));
	  hSR_Clone[2][0]->Write(TString::Format("SR_oldMva_%s_expYield", varNameReco.Data()));
	  hCR_Clone[2][0]->Write(TString::Format("CR_oldMva_%s_expYield", varNameReco.Data()));
  }
  hSR[1][0]->Write(TString::Format("SR_deepAK8_%0.2f_%s",deepAK8CutValue,varNameReco.Data()));
  hCR[1][0]->Write(TString::Format("CR_deepAK8_%0.2f_%s",deepAK8CutValue, varNameReco.Data()));
  h1Btag[1][0]->Write(TString::Format("h1Btag_deepAK8_%0.2f_%s",deepAK8CutValue, varNameReco.Data()));
  hSR_Clone[1][0]->Write(TString::Format("SR_deepAK8_%0.2f_%s_expYield",deepAK8CutValue, varNameReco.Data()));
  hCR_Clone[1][0]->Write(TString::Format("CR_deepAK8_%0.2f_%s_expYield",deepAK8CutValue, varNameReco.Data()));
  
  /*
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
	
	if(saveTTagger)
	{
	   c1->Write(TString::Format("TRatioPlot_tTagger_%s", varNameReco.Data()));
	   c3->Write(TString::Format("TRatioPlot_oldMva_%s", varNameReco.Data()));
	}
    c2->Write(TString::Format("TRatioPlot_deepAK8_%0.2f_%s",deepAK8CutValue,varNameReco.Data()));
  }
  */
  outFile->Close();
  listOfFiles.clear();
  XSEC.clear();
  histoNames.clear();
  
}