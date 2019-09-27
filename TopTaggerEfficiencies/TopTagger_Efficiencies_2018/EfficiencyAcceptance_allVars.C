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
  listOfFiles.push_back(eosPath+"TT_Mtt-700to1000_TuneCP5_PSweights_13TeV-powheg-pythia8.root");
  listOfFiles.push_back(eosPath+"TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root");
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

void EfficiencyAcceptance_allVars(float selMvaCut=-0.1, bool saveTtagger= true)
{
  initGlobals();
  const int NVAR =7;
  const int N_MJJ = 10;
  const int N_PTJJ = 9;
  const int N_YJJ = 8;
  const int N_PT = 10;
  const int N_JETY = 12;
  const int N_JETMASS = 40;
  



  int NBINS[NVAR] = {N_MJJ, N_PTJJ, N_YJJ, N_PT, N_PT ,N_JETY, N_JETY};
  std::vector< std::vector <Float_t> > const BND = {{1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, 4000, 5000}, //mjj
											 {0,60,150,300,450,600,750,950,1100,1300}, //ptjj
											 {-2.4,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.4}, //yjj
											 {400,450,500,570,650,750,850,950,1100,1300,1500}, //jetPt0	
											 {400,450,500,570,650,750,850,950,1100,1300,1500}, //jetPt1
											 {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}, //jetY0
											 {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}}; //jetY1

  TString varReco[NVAR]   = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1"}; 
  TString varParton[NVAR] = {"mTTbarParton", "ptTTbarParton", "yTTbarParton","partonPt0", "partonPt1", "partonY0", "partonY1"}; 
  TH1F *h_out_parton[listOfFiles.size()][NVAR];
  TH1F *h_out_reco[listOfFiles.size()][NVAR];
  
  TH1F* h_denominator_parton[listOfFiles.size()][NVAR]; //one histogram for parton denominator because its the same for all types (tTagger, deepAK8, oldMva)
  TH1F* h_denominator_reco[listOfFiles.size()][NVAR];
  std::vector<float> weights;
  
  //int sizeBins = 10;
  //float BND[11] = {1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, 4000, 5000};
  //int MIN = 1000; 
  //int MAX = 5000;
  
  for(int f=0; f<listOfFiles.size(); f++)
  {
  	for(int ivar =0; ivar<NVAR; ivar++)
  	{
	  	int sizeBins = NBINS[ivar];
	    float tempBND[NBINS[ivar]+1];
	    std::copy(BND[ivar].begin(), BND[ivar].end(), tempBND);	


		h_out_parton[f][ivar] = new TH1F(TString::Format("%s_%s_parton_%s", histoNames[f].Data(), "tTagger",varParton[ivar].Data()), TString::Format("%s_%s_parton_%s", histoNames[f].Data(), "tTagger",varParton[ivar].Data()), sizeBins, tempBND);
		h_out_reco[f][ivar] = new TH1F(TString::Format("%s_%s_reco_%s", histoNames[f].Data(), "tTagger",varReco[ivar].Data()), TString::Format("%s_%s_reco_%s", histoNames[f].Data(), "tTagger",varReco[ivar].Data()), sizeBins, tempBND);
	    
	    h_denominator_parton[f][ivar] = new TH1F(TString::Format("denom_%s_parton_%s",histoNames[f].Data(), varParton[ivar].Data()), TString::Format("denom_%s_parton_%s", histoNames[f].Data(),varParton[ivar].Data()), sizeBins, tempBND);
		h_denominator_parton[f][ivar]->Sumw2();
		    
		h_denominator_reco[f][ivar] = new TH1F(TString::Format("denom_%s_%s_reco_%s","tTagger",histoNames[f].Data(),varReco[ivar].Data()), TString::Format("denom_%s_%s_reco_%s","tTagger", histoNames[f].Data(),varReco[ivar].Data()), sizeBins, tempBND);
		h_denominator_reco[f][ivar]->Sumw2();

    }
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
    //TFile *file = TFile::Open(eosPath+listOfFiles[f]);
    TFile *file = TFile::Open(listOfFiles[f]);
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
	int evtNo;
	trIN->SetBranchAddress("evtNo", &evtNo);
	
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
	
	std::vector<int> indexes(0);
	std::vector<int> events(0);
    int decade(0);
    int NN = trIN->GetEntries();
    //NN = 100000;
    std::cout<<"Entries: "<<NN<<std::endl;
	std::vector<float> xRecoAll(0);
	std::vector<float> xPartonAll(0);

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
	   
	   xPartonAll.clear();
	   xRecoAll.clear();
		

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
	if(isMatched > 1)
    {
   	   	int leadingPt =0;
      	int subleadingPt = 1;

  		if((*pt_)[0] < (*pt_)[1])
  		{
   		    subleadingPt =0;
   		    leadingPt = 1;
   		}
    	xRecoAll.push_back(mJJ);
		xRecoAll.push_back(ptJJ);
		xRecoAll.push_back(yJJ);
		xRecoAll.push_back((*pt_)[leadingPt]);
		xRecoAll.push_back((*pt_)[subleadingPt]);
		//xRecoAll.push_back(fabs((*y_)[0]));
		//xRecoAll.push_back(fabs((*y_)[1])); 
				
		xPartonAll.push_back(mTTbarParton);
		xPartonAll.push_back(ptTTbarParton);
		xPartonAll.push_back(yTTbarParton);
		xPartonAll.push_back((*partonPt_)[leadingPt]);
		xPartonAll.push_back((*partonPt_)[subleadingPt]);
		//xPartonAll.push_back(fabs((*partonY)[0]));
		//xPartonAll.push_back(fabs((*partonY)[1])); 
	  		

		
	  //---------------------------end of MATCHING---------------------------------------------------------
	  bool recoCuts, partonCuts; 
	  bool massCut = (*mass_)[0] > 120 && (*mass_)[0] < 220 && (*mass_)[1] > 120 && (*mass_)[1] < 220;
	  bool tTaggerCut = (*jetTtag_)[0] > selMvaCut && (*jetTtag_)[1] > selMvaCut;
	  recoCuts = nJets > 1 && fabs((*eta_)[0]) < 2.4 && fabs((*eta_)[1]) <2.4 && (*pt_)[0] > 400 && (*pt_)[1] > 400 && mJJ > 1000 && massCut && nLeptons==0 && (*bit)[5];
	  partonCuts = fabs((*partonEta_)[0]) < 2.4 && fabs((*partonEta_)[1]) <2.4 && (*partonPt_)[0] > 400 && (*partonPt_)[1] > 400 && mTTbarParton > 1000;
	  bool deepCSV = (((*jetBtagSub0DCSVbb_)[0] + (*jetBtagSub0DCSVbbb_)[0])> deepCSVFloat || ((*jetBtagSub1DCSVbb_)[0] + (*jetBtagSub1DCSVbbb_)[0])> deepCSVFloat) && 
					 (((*jetBtagSub0DCSVbb_)[1] + (*jetBtagSub0DCSVbbb_)[1])> deepCSVFloat || ((*jetBtagSub1DCSVbb_)[1] + (*jetBtagSub1DCSVbbb_)[1])> deepCSVFloat);
			
      bool btagCut;			
	  btagCut = deepCSV;
	   
	  //qcout<<"----------"<<endl;
	  
		  //fill the denominators
		  //1. denominator passing only reco cuts for topTagger
		  if(recoCuts && btagCut && tTaggerCut)
		  {
		  	for(int ivar = 0; ivar < NVAR-2; ivar++)
	  		{
				h_denominator_reco[f][ivar]->Fill(xRecoAll[ivar], genEvtWeight);
			}
		  } 
		  //2. fill the histograms pass reco and parton cuts numerators for efficiencies and acceptance
	      if(partonCuts && recoCuts && tTaggerCut && btagCut)
		  {
			  	for(int ivar = 0; ivar < NVAR-2; ivar++)
	  			{
			  	   //if(ivar == 0)
			  	   //{
			  	   //	indexes.push_back(iev);
			  	   //	events.push_back(evtNo);
			  	   //}
				   h_out_parton[f][ivar]->Fill(xPartonAll[ivar], genEvtWeight);
				   h_out_reco[f][ivar]->Fill(xRecoAll[ivar], genEvtWeight);
				}//---- end of the ivar loop
			  	
	      }//----- end of selection cuts parton and reco 
	  
	 }//----- end of is matched 	
    }//---end the event for

	
	//--------------------------------------------START OF EVENT COUNTER LOOP -------------------------------------------------------------------

  //now another for that fills the denominators for the parton efficiencies 
  //loop over other tree -> eventCounter
  TTree *trCnt = (TTree*)file->Get("eventCounter/events");
  float ptTTbarPartonCnt(0), mTTbarPartonCnt(0), yTTbarPartonCnt(0);
  float partonPtCnt[2], partonEtaCnt[2],partonYCnt[2]; 	
  float genEvtWeightCnt;
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
	  
	  std::vector<float> xPartonAllCnt(0);
	  xPartonAllCnt.clear();
	  bool partonCuts = fabs(partonEtaCnt[0]) < 2.4 && fabs(partonEtaCnt[1]) <2.4 && partonPtCnt[0] > 400 && partonPtCnt[1] > 400 && mTTbarPartonCnt > 1000;
	  
	  xPartonAllCnt.push_back(mTTbarPartonCnt);
	  xPartonAllCnt.push_back(ptTTbarPartonCnt);
	  xPartonAllCnt.push_back(yTTbarPartonCnt);

	  if(partonPtCnt[0] >= partonPtCnt[1])
	  {
	    xPartonAllCnt.push_back(partonPtCnt[0]);
	    xPartonAllCnt.push_back(partonPtCnt[1]);
	  }
	  else
	  {
	  	xPartonAllCnt.push_back(partonPtCnt[1]);
	    xPartonAllCnt.push_back(partonPtCnt[0]);
	  }
	  //xPartonAll.push_back(partonYCnt[0]);
	  //xPartonAll.push_back(partonYCnt[1]);

	  for(int ivar = 0; ivar < NVAR-2; ivar++)
	  {
	  	if(partonCuts)
			h_denominator_parton[f][ivar]->Fill(xPartonAllCnt[ivar], genEvtWeightCnt);
	  }
  }

  //--------------------------------------------END OF EVENT COUNTER LOOP -------------------------------------------------------------------

/*
std::ofstream myfile;
myfile.open (TString::Format("entries_%s.txt",histoNames[f].Data()), std::ios_base::app);
for(int i=0; i<events.size(); i++)
{ 
    myfile << "i: " <<indexes[i]<<" event: " << events[i];
    myfile << std::endl;
}
indexes.clear();
events.clear();
myfile.close();
*/


}//----end of file loop
  
  
  for(int ivar =0; ivar<NVAR-2; ivar++)
  {

	//for every slice
	for(int j=0; j<listOfFiles.size(); j++)
    {
      h_out_parton[j][ivar]->Scale(weights[j]*LUMI);
      h_out_reco[j][ivar]->Scale(weights[j]*LUMI);
      h_denominator_reco[j][ivar]->Scale(weights[j]*LUMI);
    }
    
    for(int j=1; j<listOfFiles.size(); j++)
	{
	  //Add them to get the whole phase space
	  h_out_parton[0][ivar]->Add(h_out_parton[j][ivar]);
      //std::cout<<h_out_reco[i][0]->GetName()<<" "<<h_out_reco[i][j]->GetName()<<std::endl;
      h_out_reco[0][ivar]->Add(h_out_reco[j][ivar]);
      h_denominator_reco[0][ivar]->Add(h_denominator_reco[j][ivar]);
	}

  
  
  
  h_denominator_parton[0][ivar]->Scale(weights[0]*LUMI);
  for(int f=1; f<listOfFiles.size(); f++) 
  {
    h_denominator_parton[f][ivar]->Scale(weights[f]*LUMI);
	h_denominator_parton[0][ivar]->Add(h_denominator_parton[f][ivar]);
  }
  }
  TEfficiency *efficiency_parton[NVAR], *efficiency_reco[NVAR];
 
  //efficiency for parton quantity and for topTagger (new)

  for(int ivar = 0; ivar< NVAR-2; ivar++)
  {
  	if(h_denominator_parton[0][ivar]->GetBinContent(0) > 0)
  		h_denominator_parton[0][ivar]->SetBinContent(0,0.0);
  	if(h_out_parton[0][ivar]->GetBinContent(0) > 0)
  		h_out_parton[0][ivar]->SetBinContent(0,0.0);

  	if(h_denominator_reco[0][ivar]->GetBinContent(0) > 0)
  		h_denominator_reco[0][ivar]->SetBinContent(0,0.0);
  	if(h_out_reco[0][ivar]->GetBinContent(0) > 0)
  		h_out_reco[0][ivar]->SetBinContent(0,0.0);
  }


  for(int ivar = 0; ivar< NVAR-2; ivar++)
  {




  cout<<"--------"<<endl;
  cout<<"parton "<<ivar<<endl;
  efficiency_parton[ivar]  = new TEfficiency(*h_out_parton[0][ivar], *h_denominator_parton[0][ivar]);
  if (ivar ==0 || ivar ==1 || ivar ==3 || ivar ==4) efficiency_parton[ivar]->SetTitle(TString::Format("Parton_tTagger_%s;%s (GeV);Efficiency",varParton[ivar].Data(),varParton[ivar].Data()));
  else efficiency_parton[ivar]->SetTitle(TString::Format("Parton_tTagger_%s;%s ;Efficiency",varParton[ivar].Data(),varParton[ivar].Data()));
  efficiency_parton[ivar]->SetStatisticOption(TEfficiency::kFNormal);
  efficiency_parton[ivar]->SetUseWeightedEvents();
  efficiency_parton[ivar]->SetLineColor(kRed);
	  
  cout<<"reco"<<endl;
  efficiency_reco[ivar]  = new TEfficiency(*h_out_reco[0][ivar], *h_denominator_reco[0][ivar]);
  if (ivar ==0 ||ivar ==1 || ivar ==3 || ivar ==4)  efficiency_reco[ivar]->SetTitle(TString::Format("Reco_tTagger_%s; %s(GeV);Efficiency",varReco[ivar].Data(),varReco[ivar].Data()));
  else efficiency_reco[ivar]->SetTitle(TString::Format("Reco_tTagger_%s; %s;Efficiency",varReco[ivar].Data(), varReco[ivar].Data()));
  efficiency_reco[ivar]->SetStatisticOption(TEfficiency::kFNormal);
  efficiency_reco[ivar]->SetUseWeightedEvents();
  efficiency_reco[ivar]->SetLineColor(kRed);

/*
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
  */
  }
  TFile *outFile;
  outFile = TFile::Open(TString::Format("Efficiencies_allVars_tTagger_%0.2f_deepCSV.root", selMvaCut), "UPDATE");
  
  for(int ivar = 0; ivar<NVAR-2; ivar++)
  {
  	

  	h_out_reco[0][ivar]->Write();
  	h_denominator_reco[0][ivar]->Write();
  	h_out_parton[0][ivar]->Write();
  	h_denominator_parton[0][ivar]->Write();
    efficiency_parton[ivar]->Write(TString::Format("Sig_Parton_tTagger_%0.2f_%s", selMvaCut,  varParton[ivar].Data()));
    efficiency_reco[ivar]->Write(TString::Format("Sig_Reco_tTagger_%0.2f_%s", selMvaCut,  varReco[ivar].Data()));
  }//end of ivar

  outFile->Close();
  listOfFiles.clear();
  histoNames.clear();
  XSEC.clear();
 }
