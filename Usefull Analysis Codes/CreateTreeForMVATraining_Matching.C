#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TMVA/Tools.h"

using std::cin;
using std::cout;
using std::endl;

//float deltaRho


void EfficiencyPurity_topMatching(TString TREENAME, bool matching = false)
{
  TFile *inf     = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Signal/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_legacy2016_copy.root");
  TTree *trIN    = (TTree*)inf->Get("events");

  int nJets,nLeptons;
  float genEvtWeight;
  vector<bool>  *bit(0);
  vector<bool>  *matchedJet(0);
  vector<float> *pt(0),*tau3(0),*tau2(0),*tau1(0),*mass(0);
  vector<float> *jetMassSub0(0), *jetMassSub1(0);

  float pt_(-1),tau3_(-1),tau2_(-1),tau1_(-1),mass_(-1), jetMassSub0_(-1), jetMassSub1_(-1);
  float jetDr_(-1);

  //matching info
  vector<float> *jetPhi(0), *jetEta(0);
  vector<int> *partonId(0), *partonMatchIdx(0);
  vector<float> *partonEta(0), *partonPhi(0), *partonMatchDR(0), *partonPt(0);
  std::vector<int> *addedIndexes = new std::vector<int>(0);
                
  //------- input tree --------------
  trIN->SetBranchAddress("nJets"          ,&nJets);
  trIN->SetBranchAddress("nLeptons"       ,&nLeptons);
  trIN->SetBranchAddress("jetPt"          ,&pt);
  trIN->SetBranchAddress("jetEta"         ,&jetEta);
  trIN->SetBranchAddress("jetPhi"         ,&jetPhi);
  trIN->SetBranchAddress("jetMassSoftDrop",&mass);
  trIN->SetBranchAddress("jetTau3"        ,&tau3);
  trIN->SetBranchAddress("jetTau2"        ,&tau2);
  trIN->SetBranchAddress("jetTau1"        ,&tau1);
  trIN->SetBranchAddress("triggerBit"     ,&bit);
  trIN->SetBranchAddress("genEvtWeight"   ,&genEvtWeight);
  trIN->SetBranchAddress("jetMassSub0"    ,&jetMassSub0);
  trIN->SetBranchAddress("jetMassSub1"    ,&jetMassSub1);
  
  
  trIn->SetBranchAddress("partonPt"		  ,&partonPt);
  trIn->SetBranchAddress("partonEta"	  ,&partonEta);
  trIn->SetBranchAddress("partonPhi"	  ,&partonPhi);
  
  trIN->SetBranchAddress("partonMatchIdx" ,&partonMatchIdx);
  trIN->SetBranchAddress("partonMatchDR"  ,&partonMatchDR);
 

  int decade(0);
  int NN = trIN->GetEntries();
  std::vector<int> *jetMatchedIndexes = new std::vector<int>(0);
  std::vector<float> *jetMatchedDr = new std::vector<float>(0);
  //int NN = 10000;
  
  TH1F *h_Eff[2], *h_Purity[2];
  h_Eff[0] = new TH1F("partons matched with a reco jet passing the selection criteria", "partons matched with a reco jet passing the selection criteria", 50, 200, 2500);
  h_Eff[1] = new TH1F("all partons", "all partons", 50, 200, 2500);
  
  h_Purity[0] = new TH1F("number of selected reco jets that are matched with a parton", "number of selected reco jets that are matched with a parton", 50, 200, 2500);
  h_Purity[1] = new TH1F("number of selected reco jets", "number of selected reco jets", 50, 200, 2500);
  
  
  cout<<"Reading "<<NN<<" entries"<<endl;
for(int iev=0;iev<NN;iev++) {
    double progress = 10.0*iev/(1.0*NN);
    int k = TMath::FloorNint(progress); 
    if (k > decade) 
      cout<<10*k<<" %"<<endl;
    decade = k;
    trIN->GetEntry(iev);
   
   //for every event fill the denominator of the efficiency (partonPt)
   int partonSize = partonPt->size();
    for(int k= 0; k<=partonSize; k++)
	{
		h_Eff[1] ->Fill((*partonPt)[k]);
	}
	
    //for every jet do the matching or store its properties
    for(int i=0; i<nJets; i++)
    {
      pt_          =-1;
      mass_        =-1;
      tau3_        =-1;
      tau2_        =-1;
      tau1_        =-1;
      jetMassSub0_ =-1;
      jetMassSub1_ =-1;
      
	  
	  //here I need to put selection cuts for jets
	  if(mva > 0.8 && fabs((*jetEta)[i]) && jetPt[i] > 400)
	  {
		    h_Purity[1] ->Fill((*jetPt)[i]);
			
			//fill the denominator of the purity (jetPt, reco variable in general);
			jetMatchedIndexes->clear();
			jetMatchedDr->clear();
			std::vector<int>::iterator it = std::find(partonMatchIdx->begin(), partonMatchIdx->end(), i);
			//get all entries that match our jet.
			while(it != partonMatchIdx->end())
			{
			  int index = it - partonMatchIdx->begin();
			  jetMatchedIndexes->push_back(index); //has the positions where I found the jet i in partonMatchedIdx
			  jetMatchedDr->push_back((*partonMatchDR)[index]); //same here foer the DR: DR that correspond to the jet i
			  ++it;
			  it = std::find(it, partonMatchIdx->end(), i);
			}
			//if we actually selected something
			if(jetMatchedIndexes->size() > 0)
			{
			  float dRmin = (*jetMatchedDr)[0];
			  float indexMin = (*jetMatchedIndexes)[0];
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
				jetDr_ = dRmin;
				pt_    = (*pt)[(*partonMatchIdx)[indexMin]];
				mass_  = (*mass)[(*partonMatchIdx)[indexMin]];
				tau3_  = (*tau3)[(*partonMatchIdx)[indexMin]];
				tau2_  = (*tau2)[(*partonMatchIdx)[indexMin]];
				tau1_  = (*tau1)[(*partonMatchIdx)[indexMin]];
				jetMassSub0_ = (*jetMassSub0)[(*partonMatchIdx)[indexMin]];
				jetMassSub1_ = (*jetMassSub1)[(*partonMatchIdx)[indexMin]];
				
				//fill (*partonPt)[indexMin] numerator efficiency
				h_Eff[0]->Fill((*partonPt)[indexMin]);
				//fill purity (*jetPt)[(*partonMatchIdx)[indexMin]] purity
				h_Purity[0]->Fill((*jetPt)[(*partonMatchIdx)[indexMin]];
			  }
			}
		  
	  }//----if selection 
    }//for --- an all jets
} // for---- all events

  cout<<"Problematic events: "<<jetsWithProblem<<endl;
  cout<<"Total top jets: "<<totalTopJets<<std::endl;
  cout<<"==== Found "<<trOUT->GetEntries()<<" events ===="<<endl; 
  outf->cd();

  
  TEfficiency *efficiency = new TEfficiency(h_Eff[0], h_Eff[1]);
  TEfficiency *purity	  = new TEfficiency(h_Purity[0], h_Purity[1]);
  
  TCanvas *can_eff = new TCanvas("efficiency canvas","efficiency canvas", 900,600);
  efficiency->Draw();

  TCanvas *can_purity = new TCanvas("purity canvas", "purity canvas", 900,600);
  purity->Draw();
  
  //inf->Close();
}






















