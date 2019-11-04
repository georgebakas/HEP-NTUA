#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"

using std::cin;
using std::cout;
using std::endl;


std::vector<TString> listOfFiles;
std::vector<float> XSEC;
std::vector<TString> histoNames;
std::vector<TString> fileHistoNames;
float LUMI = 35922;
TString eosPath;
int selection;
float deepCSVFloat = 0.6321;
float selMvaCut = 0.2;

void initFiles()
{
	eosPath = "/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/Signal/";
	listOfFiles.push_back("TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");
    listOfFiles.push_back("TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");
}
void initXSec()
{
	XSEC.push_back(69.64);
	XSEC.push_back(16.74);	
}

void initFileHistoNames()
{
	fileHistoNames.push_back("Signal_histo_Mtt-700-1000");
	fileHistoNames.push_back("Signal_histo_Mtt-1000-Inf");
}

void initHistoNames()
{
	histoNames.push_back("b_quarks"); //0
	histoNames.push_back("c_quarks"); //1
	histoNames.push_back("uds_quarks"); //2
	histoNames.push_back("gluons"); //3
}

void init()
{
	initFiles();
	initXSec();
	initFileHistoNames();
	initHistoNames();
}
void CalculateBTaggingEfficiency()
{

  const int N_ETA = 10;
  const int N_PT = 10;
  float BND_PT[N_PT+1]= {};
  float BND_ETA[N_ETA+1] = {}; 
  
  gStyle->SetOptStat(0);
  init();
  TH2F *hBTagParton[listOfFiles.size()][4], *hBTagRecoParton[listOfFiles.size()][4]; 
  TH2F *hBTagReco[listOfFiles.size()];
 vector<float> weights(0);

 std::vector<int> events(0);
 for (int f =0; f<listOfFiles.size(); f++)
 {

	TFile *inf = TFile::Open(eosPath+listOfFiles[f]);

	//book b-tagging efficiency histograms
	for(int ihist = 0; ihist<4; ihist++)
	{
		hBTagParton[f][ihist] = new TH2F(TString::Format("hTagParton_%s_%s", fileHistoNames[f].Data(),histoNames[ihist].Data()), 
			TString::Format("hTagParton_%s_%s", fileHistoNames[f].Data(),histoNames[ihist].Data()), N_PT,0,1500,N_ETA,0,2.4 );
    hBTagRecoParton[f][ihist] = new TH2F(TString::Format("hTagRecoParton_%s_%s", fileHistoNames[f].Data(),histoNames[ihist].Data()), 
      TString::Format("hTagRecoParton_%s_%s", fileHistoNames[f].Data(),histoNames[ihist].Data()), N_PT,0,1500,N_ETA,0,2.4 );
	}

  hBTagReco[f] = new TH2F(TString::Format("hBTagReco_%s", fileHistoNames[f].Data()), 
			TString::Format("hBTagReco_%s", fileHistoNames[f].Data()), N_PT, 0,1500, N_ETA, 0,2.4);

  TTree *trIN    = (TTree*)inf->Get("boosted/events");
  float NORM = ((TH1F*)inf->Get("eventCounter/GenEventWeight"))->GetSumOfWeights(); 
  weights.push_back((XSEC[f]/NORM) * LUMI );  
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

  std::vector<float> *jetFlavorSub0(0), *jetFlavorSub1(0);
  std::vector<float> *jetPtSub0(0), *jetPtSub1(0);
  std::vector<float> *jetEtaSub0(0), *jetEtaSub1(0);
  //------- input tree --------------
  int evtNo;
  trIN->SetBranchAddress("evtNo"          ,&evtNo);
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
  trIN->SetBranchAddress("mJJ"        ,&mJJ);
  trIN->SetBranchAddress("yJJ"        ,&yJJ);
  trIN->SetBranchAddress("ptJJ"         ,&ptJJ);
  trIN->SetBranchAddress("jetBtagSub0"    ,&jetBtagSub0);
  trIN->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1);
  trIN->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
  trIN->SetBranchAddress("mva"          ,&mva);
  trIN->SetBranchAddress("category"       ,&category);
  trIN->SetBranchAddress("deepAK8Tagger"  ,&deepAK8);
  trIN->SetBranchAddress("jetTtagCategory",&jetTtag);

  trIN->SetBranchAddress("mTTbarParton" ,&mTTbarParton);
  trIN->SetBranchAddress("yTTbarParton" ,&yTTbarParton);
  trIN->SetBranchAddress("ptTTbarParton"  ,&ptTTbarParton);
  trIN->SetBranchAddress("partonPt"     ,&partonPt);
  trIN->SetBranchAddress("partonEta"    ,&partonEta);
  trIN->SetBranchAddress("partonMass"     ,&partonMass);
  trIN->SetBranchAddress("partonMatchDR"  ,&partonMatchDR);
  trIN->SetBranchAddress("partonMatchIdx" ,&partonMatchIdx);
  trIN->SetBranchAddress("partonPhi"      ,&partonPhi);
  trIN->SetBranchAddress("jetFlavorSub0"  ,&jetFlavorSub0);
  trIN->SetBranchAddress("jetFlavorSub1"  ,&jetFlavorSub1);
  trIN->SetBranchAddress("jetPtSub0"  ,&jetPtSub0);
  trIN->SetBranchAddress("jetPtSub1"  ,&jetPtSub1);
  trIN->SetBranchAddress("jetEtaSub0"  ,&jetEtaSub0);
  trIN->SetBranchAddress("jetEtaSub1"  ,&jetEtaSub1);
  
  //deepCSV
  trIN->SetBranchAddress("jetBtagSub0DCSVbb" ,&jetBtagSub0DCSVbb);
  trIN->SetBranchAddress("jetBtagSub1DCSVbb" ,&jetBtagSub1DCSVbb);
  trIN->SetBranchAddress("jetBtagSub0DCSVbbb",&jetBtagSub0DCSVbbb);
  trIN->SetBranchAddress("jetBtagSub1DCSVbbb",&jetBtagSub1DCSVbbb);

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
  std::vector<float> *jetPtSub0_ = new std::vector<float>(0);
  std::vector<float> *jetPtSub1_ = new std::vector<float>(0);
  std::vector<float> *jetEtaSub0_ = new std::vector<float>(0);
  std::vector<float> *jetEtaSub1_ = new std::vector<float>(0);

  std::vector<float> *partonPt_ = new std::vector<float>(0);
  std::vector<float> *partonEta_ = new std::vector<float>(0);
  std::vector<float> *partonY_ = new std::vector<float>(0);
  std::vector<float> *partonMass_ = new std::vector<float>(0);
  std::vector<float> *partonPhi_ = new std::vector<float>(0);
  std::vector<float> *jetFlavorSub0_ = new std::vector<float>(0);
  std::vector<float> *jetFlavorSub1_ = new std::vector<float>(0);
  
  std::vector<float> *jetBtagSub0DCSVbb_ = new std::vector<float>(0);
  std::vector<float> *jetBtagSub1DCSVbb_ = new std::vector<float>(0);
  std::vector<float> *jetBtagSub0DCSVbbb_ = new std::vector<float>(0);
  std::vector<float> *jetBtagSub1DCSVbbb_ = new std::vector<float>(0);



  cout<<"Reading "<<NN<<" entries"<<endl;
  for(int iev=0;iev<NN;iev++) 
  {
    double progress = 10.0*iev/(1.0*NN);
    int k = TMath::FloorNint(progress); 
    if (k > decade) 
      cout<<10*k<<" %"<<endl;
    decade = k;
    trIN->GetEntry(iev);
    //trIN->GetEntry(277816);

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
	jetFlavorSub1_->clear();
	jetFlavorSub0_->clear();
	jetPtSub0_->clear();
	jetPtSub1_->clear();
	jetEtaSub0_->clear();
	jetEtaSub1_->clear();
  
    jetBtagSub0DCSVbb_->clear();
    jetBtagSub1DCSVbb_->clear();
    jetBtagSub0DCSVbbb_->clear();
    jetBtagSub1DCSVbbb_->clear();


  if (nJets >1)
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

	        jetPtSub0_->push_back((*jetPtSub0)[(*partonMatchIdx)[indexMin]]);
	        jetPtSub1_->push_back((*jetPtSub1)[(*partonMatchIdx)[indexMin]]);

	        jetEtaSub0_->push_back((*jetEtaSub0)[(*partonMatchIdx)[indexMin]]);
	        jetEtaSub1_->push_back((*jetEtaSub1)[(*partonMatchIdx)[indexMin]]);
	  		  jetFlavorSub0_->push_back((*jetFlavorSub0)[(*partonMatchIdx)[indexMin]]);
	        jetFlavorSub1_->push_back((*jetFlavorSub1)[(*partonMatchIdx)[indexMin]]);

	        jetBtagSub0DCSVbb_->push_back((*jetBtagSub0DCSVbb)[(*partonMatchIdx)[indexMin]]);
	        jetBtagSub1DCSVbb_->push_back((*jetBtagSub1DCSVbb)[(*partonMatchIdx)[indexMin]]);
	        jetBtagSub0DCSVbbb_->push_back((*jetBtagSub0DCSVbbb)[(*partonMatchIdx)[indexMin]]);
	        jetBtagSub1DCSVbbb_->push_back((*jetBtagSub1DCSVbbb)[(*partonMatchIdx)[indexMin]]);
	  
	        //PARTON MATCHED
	        partonPt_->push_back((*partonPt)[indexMin]);
	        partonMass_->push_back((*partonMass)[indexMin]);
	        partonPhi_->push_back((*partonPhi)[indexMin]);
	        partonEta_->push_back((*partonEta)[indexMin]);
	        
	        //here misssing partonY
	      }     
	     }//----end of if jetMatchedIndexes > 0
	    }//----end of for on all jets for mathching


	bool recoCuts, partonCuts, massCut, tTaggerCut, deepCSV, btag1DeepCSV, revertBtagDeepCSV, tightMassCut;
	float dCSVScoreSub0[2], dCSVScoreSub1[2];
    dCSVScoreSub0[0] = (*jetBtagSub0DCSVbb_)[0] + (*jetBtagSub0DCSVbbb_)[0];
    dCSVScoreSub0[1] = (*jetBtagSub0DCSVbb_)[1] + (*jetBtagSub0DCSVbbb_)[1];
    dCSVScoreSub1[0] = (*jetBtagSub1DCSVbb_)[0] + (*jetBtagSub1DCSVbbb_)[0];
    dCSVScoreSub1[1] = (*jetBtagSub1DCSVbb_)[1] + (*jetBtagSub1DCSVbbb_)[1];
    
    recoCuts   = fabs((*eta_)[0]) < 2.4 && fabs((*eta_)[1]) <2.4 && (*pt_)[0] > 400 && (*pt_)[1] > 400 && nLeptons==0 && mJJ > 1000 && nJets > 1 && (*bit)[2];
    partonCuts = fabs((*partonEta_)[0]) < 2.4 && fabs((*partonEta_)[1]) <2.4 && (*partonPt_)[0] > 400 && (*partonPt_)[1] > 400 && mTTbarParton > 1000;
    massCut    = (*mass_)[0] > 50 && (*mass_)[0] < 300 && (*mass_)[1] > 50 && (*mass_)[1] < 300;
    tightMassCut    = (*mass_)[0] > 120 && (*mass_)[0] < 220 && (*mass_)[1] > 120 && (*mass_)[1] < 220;
    tTaggerCut = (*jetTtag_)[0] > selMvaCut && (*jetTtag_)[1] > selMvaCut;
    //2 btag category with deepCSV
    deepCSV    = (((*jetBtagSub0DCSVbb_)[0] + (*jetBtagSub0DCSVbbb_)[0])> deepCSVFloat || ((*jetBtagSub1DCSVbb_)[0] + (*jetBtagSub1DCSVbbb_)[0])> deepCSVFloat) && 
           (((*jetBtagSub0DCSVbb_)[1] + (*jetBtagSub0DCSVbbb_)[1])> deepCSVFloat || ((*jetBtagSub1DCSVbb_)[1] + (*jetBtagSub1DCSVbbb_)[1])> deepCSVFloat);
    //1 btag category with deepCSV                
   	btag1DeepCSV  = ((dCSVScoreSub0[0] > deepCSVFloat || dCSVScoreSub1[0] > deepCSVFloat) && (dCSVScoreSub0[1] < deepCSVFloat && dCSVScoreSub1[1] < deepCSVFloat)) ||
            ((dCSVScoreSub0[0] < deepCSVFloat && dCSVScoreSub1[0] < deepCSVFloat) && (dCSVScoreSub0[1] > deepCSVFloat || dCSVScoreSub1[1] > deepCSVFloat));
   
    //0 btag category with deepCSV
    revertBtagDeepCSV = (dCSVScoreSub0[0] < deepCSVFloat &&  dCSVScoreSub1[0] < deepCSVFloat) && (dCSVScoreSub0[1] < deepCSVFloat && dCSVScoreSub1[1] < deepCSVFloat);
/*
    cout<<"---------------------"<<endl;
    cout<<"dCSVScoreSub0 jet0: "<<dCSVScoreSub0[0]<<endl;
    cout<<"dCSVScoreSub0 jet1: "<<dCSVScoreSub0[1]<<endl;
    cout<<"dCSVScoreSub1 jet0: "<<dCSVScoreSub1[0]<<endl;
    cout<<"dCSVScoreSub1 jet1: "<<dCSVScoreSub1[1]<<endl;

    cout<<"abs flavour Sub0 parton0: "<<fabs((*jetFlavorSub0_)[0])<<endl;
    cout<<"abs flavour Sub0 parton1: "<<fabs((*jetFlavorSub0_)[1])<<endl;
    cout<<"abs flavour Sub1 parton0: "<<fabs((*jetFlavorSub1_)[0])<<endl;
    cout<<"abs flavour Sub1 parton1: "<<fabs((*jetFlavorSub1_)[1])<<endl;

    cout<<"topTagger: "<<tTaggerCut<<endl;
    cout<<"massCut: "<<massCut<<endl;
    cout<<"partonCuts: "<<partonCuts<<endl;;
    cout<<"recoCuts: "<<recoCuts<<endl;

    cout<<"partonEta0: "<<fabs((*partonEta_)[0])<<endl;
    cout<<"partonEta1: "<<fabs((*partonEta_)[1])<<endl;
    cout<<"partonPt0: "<<(*partonPt_)[0]<<endl;
    cout<<"partonPt1: "<<(*partonPt_)[1]<<endl;
    cout<<"mTTbarParton: "<<mTTbarParton<<endl;
  */
    if(isMatched > 1)
    {
      if(recoCuts && tTaggerCut &&  tightMassCut && partonCuts)
      {
        //3 categories
        //1: pass reco cuts (meaning only the btagging deepCSV requirement)
        //2: pass reco and parton cuts (meaning that they pass the deepCSV and then categorize them in 4 regions depending on the parton flavour)
        //3: pass only parton cuts (meaning that we categorize 4 regions depending on the parton flavour)

        //you have to do this on the 2 leading jets on their 0 and 1 subjets
        //we do it per subjets:
        for(int ijet =0; ijet<2; ijet++)
        {         
          if(dCSVScoreSub0[ijet] > deepCSVFloat) //this means that the first subjet is btagged
          {
            //fill the region 1
            hBTagReco[f]->Fill((*jetPtSub0_)[ijet], fabs((*jetEtaSub0_)[ijet]));

            //now for category 2 for the first subjet
            if(fabs((*jetFlavorSub0_)[ijet])  ==5)
              hBTagRecoParton[f][0]->Fill((*jetPtSub0_)[ijet], fabs((*jetEtaSub0_)[ijet]), genEvtWeight);
            else if(fabs((*jetFlavorSub0_)[ijet]) == 4)
              hBTagRecoParton[f][1]->Fill((*jetPtSub0_)[ijet], fabs((*jetEtaSub0_)[ijet]), genEvtWeight);
            else if(fabs((*jetFlavorSub0_)[ijet]) == 1 || fabs((*jetFlavorSub0_)[ijet]) == 2 || fabs((*jetFlavorSub0_)[ijet]) == 3)
              hBTagRecoParton[f][2]->Fill((*jetPtSub0_)[ijet], fabs((*jetEtaSub0_)[ijet]), genEvtWeight);
            else if(fabs((*jetFlavorSub0_)[ijet]) == 21)
              hBTagRecoParton[f][3]->Fill((*jetPtSub0_)[ijet], fabs((*jetEtaSub0_)[ijet]), genEvtWeight);

          }
          if(dCSVScoreSub1[ijet] > deepCSVFloat) //this means that the second subjet is btagged
          {
            //fill the region 1
            hBTagReco[f]->Fill((*jetPtSub1_)[ijet], fabs((*jetEtaSub1_)[ijet]));

            //now for category 2 for the second subjet
            if(fabs((*jetFlavorSub1_)[ijet]) ==5)
              hBTagRecoParton[f][0]->Fill((*jetPtSub1_)[ijet], fabs((*jetEtaSub1_)[ijet]), genEvtWeight);
            else if(fabs((*jetFlavorSub1_)[ijet]) == 4)
              hBTagRecoParton[f][1]->Fill((*jetPtSub1_)[ijet], fabs((*jetEtaSub1_)[ijet]), genEvtWeight);
            else if(fabs((*jetFlavorSub1_)[ijet]) == 1 || fabs((*jetFlavorSub1_)[ijet]) == 2 || fabs((*jetFlavorSub1_)[ijet]) == 3)
              hBTagRecoParton[f][2]->Fill((*jetPtSub1_)[ijet], fabs((*jetEtaSub1_)[ijet]), genEvtWeight);
            else if(fabs((*jetFlavorSub1_)[ijet]) == 21)
              hBTagRecoParton[f][3]->Fill((*jetPtSub1_)[ijet], fabs((*jetEtaSub1_)[ijet]), genEvtWeight);
          }

          //category 3 only flavours for 1st subjet
          if(fabs((*jetFlavorSub0_)[ijet]) ==5)
            hBTagParton[f][0]->Fill((*jetPtSub0_)[ijet], fabs((*jetEtaSub0_)[ijet]), genEvtWeight);
          else if(fabs((*jetFlavorSub0_)[ijet]) == 4)
            hBTagParton[f][1]->Fill((*jetPtSub0_)[ijet], fabs((*jetEtaSub0_)[ijet]), genEvtWeight);
          else if(fabs((*jetFlavorSub0_)[ijet]) == 1 || fabs((*jetFlavorSub0_)[ijet] == 2) || fabs((*jetFlavorSub0_)[ijet]) == 3)
            hBTagParton[f][2]->Fill((*jetPtSub0_)[ijet], fabs((*jetEtaSub0_)[ijet]), genEvtWeight);
          else if(fabs((*jetFlavorSub0_)[ijet]) == 21)
            hBTagParton[f][3]->Fill((*jetPtSub0_)[ijet], fabs((*jetEtaSub0_)[ijet]), genEvtWeight);

          //category 3 only flavours for 2nd subjet
          if(fabs((*jetFlavorSub1_)[ijet]) ==5)
            hBTagParton[f][0]->Fill((*jetPtSub1_)[ijet], fabs((*jetEtaSub1_)[ijet]), genEvtWeight);      
          else if(fabs((*jetFlavorSub1_)[ijet]) == 4)
            hBTagParton[f][1]->Fill((*jetPtSub1_)[ijet], fabs((*jetEtaSub1_)[ijet]), genEvtWeight);
          else if(fabs((*jetFlavorSub1_)[ijet]) == 1 || fabs((*jetFlavorSub1_)[ijet]) == 2 || fabs((*jetFlavorSub1_)[ijet]) == 3)
            hBTagParton[f][2]->Fill((*jetPtSub1_)[ijet], fabs((*jetEtaSub1_)[ijet]), genEvtWeight);
          else if(fabs((*jetFlavorSub1_)[ijet]) == 21)
            hBTagParton[f][3]->Fill((*jetPtSub1_)[ijet], fabs((*jetEtaSub1_)[ijet]), genEvtWeight);
        

        }//end of jet loop

      	

      } //end of selection
    }//is matched


	}//----if nJets>1
  }//----end of NN loop

 }//for listoffiles.size()


 //scale them with XSEC * LUMI / GENEVENTS
 for(int f = 0; f<listOfFiles.size(); f++)
 {
 	hBTagReco[f]->GetXaxis()->SetTitle("subJet P_{T} (GeV)");
 	hBTagReco[f]->GetYaxis()->SetTitle("subJet |#eta |");
 	hBTagReco[f]->Scale(weights[f]);
 	for(int i=0; i<histoNames.size(); i++)
 	{
 		hBTagParton[f][i]->GetXaxis()->SetTitle("subJet P_{T} (GeV)");
 		hBTagParton[f][i]->GetYaxis()->SetTitle("subJet |#eta|");
 		hBTagParton[f][i]->Scale(weights[f]);

    hBTagRecoParton[f][i]->GetXaxis()->SetTitle("subJet P_{T} (GeV)");
    hBTagRecoParton[f][i]->GetYaxis()->SetTitle("subJet |#eta|");
    hBTagRecoParton[f][i]->Scale(weights[f]);
 	}
 }

 for(int f = 1; f<listOfFiles.size(); f++)
 {
 	hBTagReco[0]->Add(hBTagReco[f]);
 	for(int i=0; i<histoNames.size(); i++)
 	{
 		hBTagParton[0][i]->Add(hBTagParton[f][i]);
    hBTagRecoParton[0][i]->Add(hBTagRecoParton[f][i]);
 	}
 }

/*
 //now plot them
 TCanvas *can[5];
 for(int i =0; i<histoNames.size(); i++)
 {
 	can[i] = new TCanvas(TString::Format("can_%d",i),TString::Format("can_%d",i), 800, 600 );
 	can[i]->cd();
 	hBTagParton[0][i]->Draw("text colz");
 }
 
 can[4] = new TCanvas(TString::Format("can_%d",4),TString::Format("can_%d",4), 800, 600 );
 can[4]->cd();
 hBTagReco[0]->Draw("text colz");
 */

 TFile *outf = new TFile("BTaggingEfficiency_massCutLoose_tightMassCut.root", "RECREATE");
 outf->cd();
 //this is to save the histograms
 hBTagReco[0]->Write();
 for(int i =0; i<histoNames.size(); i++)
 {
 	hBTagParton[0][i]->Write();
  hBTagRecoParton[0][i]->Write();
 }

 cout<<"btagging efficiency eb = "<<hBTagRecoParton[0][0]->Integral() / hBTagParton[0][0]->Integral()<<endl;
 cout<<"btagging acceptance = "<<hBTagRecoParton[0][0]->Integral() / (hBTagRecoParton[0][0]->Integral() + hBTagRecoParton[0][1]->Integral()+ hBTagRecoParton[0][2]->Integral()+ hBTagRecoParton[0][3]->Integral())<<endl;

}
