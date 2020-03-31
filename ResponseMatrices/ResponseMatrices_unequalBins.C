#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1.h"
#include "TCanvas.h"


using std::cin;
using std::cout;
using std::endl;
#include <fstream>

std::vector<float> XSEC;
std::vector<TString> histoNames;
std::vector<TString> fileNames;

#include "TemplateConstants.h"
TString globalYear;
bool globalIsNominalMC;

void initXsections()
{
  if(!globalIsNominalMC)
  {
  	if(globalYear.EqualTo("2016"))
  	{
  		//XSEC.push_back(80.73);
  		//XSEC.push_back(21.35);
  		XSEC.push_back(80.78);
        XSEC.push_back(21.43);
  	}
  	else if(globalYear.EqualTo("2017"))
  	{
  		//XSEC.push_back(231.01);
  		//XSEC.push_back(61.08);
  	 	XSEC.push_back(0.);
       	XSEC.push_back(20.09);
  	} 
  	else if(globalYear.EqualTo("2018"))
  	{
  		//XSEC.push_back(114.66);
  		//XSEC.push_back(30.32);
  		XSEC.push_back(0.);
        XSEC.push_back(20.09);
  	}
  }
  else
  {
  	//this is for nominal files
 	if(globalYear.EqualTo("2016")) XSEC.push_back(832.);
 	else
 	{
 		XSEC.push_back(377.96);
 		XSEC.push_back(365.34);
 		XSEC.push_back(88.29);
 	}
  }
}

void initHistoNames()
{
  if(!globalIsNominalMC)
  {
	  histoNames.push_back("Signal_histo_Mtt_700_1000"); 
	  histoNames.push_back("Signal_histo_Mtt_1000_Inf");

	  fileNames.push_back("700-1000");
	  fileNames.push_back("1000-Inf");
  }
  else
  {
	  if(globalYear.EqualTo("2016"))
	  {
	  	fileNames.push_back("TTNominal");
	  	histoNames.push_back("Signal_histo_Nominal");
	  }
	  else
	  {
	  	fileNames.push_back("TTHadronic_0");
	  	fileNames.push_back("TTSemiLeptonic_0");
	  	fileNames.push_back("TTTo2L2Nu_0");

	  	histoNames.push_back("Signal_histo_TTHadronic_0");
	  	histoNames.push_back("Signal_histo_TTSemiLeptonic_0");
	  	histoNames.push_back("Signal_histo_TTTo2L2Nu_0");

	  }
  }
}

void ResponseMatrices_unequalBins(TString year = "2016", bool isNominalMC=false)
{
  globalIsNominalMC = isNominalMC;
  globalYear = year;
  initFilesMapping();
  initHistoNames();
  initXsections();

  float deepCSVFloat = floatConstants[TString::Format("btagWP%s",year.Data())];
  float selMvaCut = topTaggerConstants[TString::Format("topTagger%s",year.Data())];
  float LUMI = luminosity[TString::Format("luminosity%s", year.Data())];
  cout<<LUMI<<endl;
 
  std::vector< std::vector <Float_t> > const BND_gen = {{1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, 4000, 5000}, //mjj
													{0,60,150,300,450,600,750,950,1100,1300}, //ptjj
													{-2.4,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.4}, //yjj
		   	                                        {400,450,500,570,650,750,850,950,1100,1300,1500}, //jetPt0     
													{400,450,500,570,650,750,850,950,1100,1300,1500}, //jetPt1
													{0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}, //jetY0
                                                    {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}}; //jetY1
  
  std::vector< std::vector <Float_t> > const BND_reco = {{1000, 1100,1200,1300, 1400,1500, 1600,1700, 1800,1900, 2000,2200, 2400,2600, 2800,3000, 3200,3600, 4000,4500, 5000}, //mjj 21
													{0,30,60,105,150,225,300,375,450,525,600,675,750,850,950,1025,1100,1200,1300}, //ptjj 19
													{-2.4,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.4}, //yjj
		   	                                        {400,425,450,475,500,535,570,610,650,700,750,800,850,900,950,1025,1100,1200,1300,1400,1500}, //jetPt0 21   
													{400,425,450,475,500,535,570,610,650,700,750,800,850,900,950,1025,1100,1200,1300,1400,1500}, //jetPt1 21   
													{0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4}, //jetY0
                                                    {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4}}; //jetY0 25

  
  float triggerFloat;
  if(year.EqualTo("2016")) triggerFloat = 2;
  else triggerFloat = 5;

  int NBINS[BND_reco.size()];
  const int NVAR = 7;
  for (int i = 0; i<BND_reco.size(); i++) NBINS[i] = BND_reco[i].size()-1;

  int NBINSParton[BND_reco.size()];	
  for(int i =0; i<BND_reco.size(); i++) NBINSParton[i] = BND_gen[i].size()-1;
  TString varReco[NVAR]   = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1"}; 
  TString varParton[NVAR] = {"mTTbarParton", "ptTTbarParton", "yTTbarParton","partonPt0", "partonPt1", "partonY0", "partonY1"}; 
  TString varParticle[NVAR] = {"mJJGen", "ptJJGen", "yJJGen","genjetPt0", "genjetPt1", "genjetY0", "genjetY1"}; 

  std::vector<float> weights;
  TH1F *hParton[fileNames.size()][NVAR], *hParticle[fileNames.size()][NVAR], *hReco[fileNames.size()][NVAR];
  TH1F *hRecoParton[fileNames.size()][NVAR], *hPartonReco[fileNames.size()][NVAR];
  TH1F *hRecoParticle[fileNames.size()][NVAR], *hParticleReco[fileNames.size()][NVAR];
  TH2F *hPartonResponse[fileNames.size()][NVAR], *hParticleResponse[fileNames.size()][NVAR];
  cout<<fileNames.size()<<endl;

  ofstream myFile;
  myFile.open("matching.txt");

for(int f=0; f<fileNames.size(); f++)
{
    cout<<XSEC[f]<<endl;
  	//declare the histograms
  	for(int ivar =0; ivar<NVAR; ivar++)
  	{	

  		//for the hists used for acceptance I need binning as reco m
  		//for the hists used for efficiency I need binning as parton/particle n = 2m
  		 int sizeBins = NBINS[ivar];
         float tempBND[NBINS[ivar]+1];
         std::copy(BND_reco[ivar].begin(), BND_reco[ivar].end(), tempBND);

         int sizeBinsPartonParticle = NBINSParton[ivar];
         float tempBNDPartonParticle[NBINSParton[ivar]+1];
         std::copy(BND_gen[ivar].begin(), BND_gen[ivar].end(), tempBNDPartonParticle);
  		 
  		 //denominators for parton efficiency (hParton vs parton), particle eff (hParticle vs particle) and acceptance (same for both hReco vs reco)
  		 //this is events pass parton cuts vs parton quantity --> BND_gen
  		 hParton[f][ivar] = new TH1F(TString::Format("hParton_%s_%s", histoNames[f].Data(),varParton[ivar].Data()), TString::Format("hParton_%s_%s", histoNames[f].Data(),varParton[ivar].Data()), sizeBinsPartonParticle, tempBNDPartonParticle);
         //this is events pass reco cuts used as denominator vs reco quantity --> BND_reco
         hReco[f][ivar] = new TH1F(TString::Format("hReco_%s_%s", histoNames[f].Data(),varReco[ivar].Data()), TString::Format("hReco_%s_%s", histoNames[f].Data(),varReco[ivar].Data()), sizeBins, tempBND);
         //this is events pass particle cuts only used as denominator for particle eff vs particle -->BND
         hParticle[f][ivar] = new TH1F(TString::Format("hParticle_%s_%s", histoNames[f].Data(),varParticle[ivar].Data()), TString::Format("hParticle_%s_%s", histoNames[f].Data(),varParticle[ivar].Data()), sizeBinsPartonParticle, tempBNDPartonParticle);
         
         hParton[f][ivar]->Sumw2();
         hReco[f][ivar]->Sumw2();
         hParticle[f][ivar]->Sumw2();

         //numerator for parton efficiency (hPartonReco vs parton) and acceptance (hRecoParton vs reco)
         //this is RecoParton vs reco --> BND_reco
         hRecoParton[f][ivar] = new TH1F(TString::Format("hRecoParton_%s_%s", histoNames[f].Data(),varReco[ivar].Data()), TString::Format("hRecoParton_%s_%s", histoNames[f].Data(),varReco[ivar].Data()), sizeBins, tempBND);
         hPartonReco[f][ivar] = new TH1F(TString::Format("hPartonReco_%s_%s", histoNames[f].Data(),varParton[ivar].Data()), TString::Format("hPartonReco_%s_%s", histoNames[f].Data(),varParton[ivar].Data()), sizeBinsPartonParticle, tempBNDPartonParticle);

         //numerator for particle efficiency (hParticleReco vs particle) and acceptance (hRecoParticle vs reco)
         hRecoParticle[f][ivar] = new TH1F(TString::Format("hRecoParticle_%s_%s", histoNames[f].Data(),varReco[ivar].Data()), TString::Format("hRecoParticle_%s_%s", histoNames[f].Data(),varReco[ivar].Data()), sizeBins, tempBND);
         hParticleReco[f][ivar] = new TH1F(TString::Format("hParticleReco_%s_%s", histoNames[f].Data(),varParticle[ivar].Data()), TString::Format("hParticleReco_%s_%s", histoNames[f].Data(),varParticle[ivar].Data()), sizeBinsPartonParticle, tempBNDPartonParticle);

         //response matrices
         //x-axis: parton or particle and y-axis: reco (detector level)
         //x-axis will have fewer bins than y-axis
         hPartonResponse[f][ivar] = new TH2F(TString::Format("hPartonResponse%s_%s", histoNames[f].Data(),varReco[ivar].Data()), TString::Format("hPartonResponse%s_%s", histoNames[f].Data(),varReco[ivar].Data()), 
         	sizeBinsPartonParticle, tempBNDPartonParticle,sizeBins, tempBND);
         hParticleResponse[f][ivar] = new TH2F(TString::Format("hParticleResponse%s_%s", histoNames[f].Data(),varReco[ivar].Data()), TString::Format("hParticleResponse%s_%s", histoNames[f].Data(),varReco[ivar].Data()), 
         	sizeBinsPartonParticle, tempBNDPartonParticle,sizeBins, tempBND);

    }
    int nJets,nLeptons, category(0);
    vector<bool>  *bit(0),*matchedJet(0);
    //reco vars:
    std::vector<float> *jetPt(0), *jetY(0), *jetEta(0), *jetPhi(0), *jetTtag(0);
    float genEvtWeight(0);
    float mJJ(0), ptJJ(0), yJJ(0),mva(0);
    vector<float> *tau3(0),*tau2(0),*tau1(0);
	  vector<float> *jetMassSub0(0), *jetMassSub1(0), *jetBtagSub0(0), *jetBtagSub1(0);
    vector<float> *jetMassSoftDrop(0), *partonEta(0);

    //parton
    std::vector<float> *partonPt(0), *partonPhi(0), *partonMass(0),*partonMatchDR(0), *partonY(0);
    float yTTbarParton(0), ptTTbarParton(0), mTTbarParton(0);
    
    //float yTopParton[2], ptTopParton[2]; 
    std::vector<int> *addedIndexes = new std::vector<int>(0);
    vector<int> *partonId(0), *partonMatchIdx(0);

    std::vector<float> *jetBtagSub0DCSVbb(0), *jetBtagSub1DCSVbb(0);
    std::vector<float> *jetBtagSub0DCSVbbb(0), *jetBtagSub1DCSVbbb(0);
	
    //particle
    std::vector<float> *genjetPt(0), *genjetY(0), *genjetEta(0), *genSoftDropMass(0), *genjetMassSoftDrop(0);
    std::vector<int> *nSubGenJets(0);
    int nJetsGen(0);
    float mJJGen(0), ptJJGen(0), yJJGen(0);

    cout<<"Working in file: "<<files[year.Data()][fileNames[f].Data()]<<endl;
    TFile *file = TFile::Open(eospath[year.Data()]+files[year.Data()][fileNames[f].Data()]);
    
    TTree *trIN = (TTree*)file->Get("boosted/events"); 
    	
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
    //trIN->SetBranchAddress("jetMassSub1"    ,&jetMassSub1);
    trIN->SetBranchAddress("mJJ"   			,&mJJ);
    trIN->SetBranchAddress("yJJ"   			,&yJJ);
    trIN->SetBranchAddress("ptJJ"   		,&ptJJ);
	trIN->SetBranchAddress("jetBtagSub0"	,&jetBtagSub0);
    trIN->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1);
    trIN->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
  	trIN->SetBranchAddress("jetTtagCategory",&jetTtag);
  	//deepCSV
    trIN->SetBranchAddress("jetBtagSub0DCSVbb" ,&jetBtagSub0DCSVbb);
    trIN->SetBranchAddress("jetBtagSub1DCSVbb" ,&jetBtagSub1DCSVbb);
    trIN->SetBranchAddress("jetBtagSub0DCSVbbb",&jetBtagSub0DCSVbbb);
    trIN->SetBranchAddress("jetBtagSub1DCSVbbb",&jetBtagSub1DCSVbbb);
    //general
    int evtNo;
	trIN->SetBranchAddress("evtNo", &evtNo);

    //parton
	trIN->SetBranchAddress("mTTbarParton"	,&mTTbarParton);
	trIN->SetBranchAddress("yTTbarParton"	,&yTTbarParton);
	trIN->SetBranchAddress("ptTTbarParton"	,&ptTTbarParton);
	trIN->SetBranchAddress("partonPt"	    ,&partonPt);
	trIN->SetBranchAddress("partonEta"		,&partonEta);
	trIN->SetBranchAddress("partonMass"     ,&partonMass);
	trIN->SetBranchAddress("partonMatchDR"  ,&partonMatchDR);
	trIN->SetBranchAddress("partonMatchIdx" ,&partonMatchIdx);
	trIN->SetBranchAddress("partonPhi"      ,&partonPhi);
	trIN->SetBranchAddress("partonY"        ,&partonY);
    
    //particle
    trIN->SetBranchAddress("mJJGen"			,&mJJGen);
    trIN->SetBranchAddress("ptJJGen"		,&ptJJGen);
    trIN->SetBranchAddress("yJJGen"			,&yJJGen);
    trIN->SetBranchAddress("genjetPt"		,&genjetPt);
    trIN->SetBranchAddress("genjetY"		,&genjetY);
    trIN->SetBranchAddress("genjetEta"		,&genjetEta);
    trIN->SetBranchAddress("nJetsGen"		,&nJetsGen);
    trIN->SetBranchAddress("genjetMassSoftDrop", &genjetMassSoftDrop);

	
    float norm = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
	float weight = XSEC[f]/norm;
	weights.push_back(weight);
    
	//for parton matching
	std::vector<int> *jetMatchedIndexes = new std::vector<int>(0);
	std::vector<float> *jetMatchedDr = new std::vector<float>(0);
	std::vector<float> *eta_ = new std::vector<float>(0);
	std::vector<float> *phi_ = new std::vector<float>(0);
	std::vector<float> *y_   = new std::vector<float>(0);
	std::vector<float> *mass_ = new std::vector<float>(0);
	std::vector<float> *pt_ = new std::vector<float>(0);
	std::vector<float> *jetTtag_ = new std::vector<float>(0);
	std::vector<float> *jetBtagSub0_ = new std::vector<float>(0);
	std::vector<float> *jetBtagSub1_ = new std::vector<float>(0);

	std::vector<float> *partonPt_ = new std::vector<float>(0);
	std::vector<float> *partonEta_ = new std::vector<float>(0);
	std::vector<float> *partonMass_ = new std::vector<float>(0);
	std::vector<float> *partonPhi_ = new std::vector<float>(0);
	std::vector<float> *partonY_ = new std::vector<float>(0);
	
	std::vector<float> *jetBtagSub0DCSVbb_ = new std::vector<float>(0);
	std::vector<float> *jetBtagSub1DCSVbb_ = new std::vector<float>(0);
	std::vector<float> *jetBtagSub0DCSVbbb_ = new std::vector<float>(0);
	std::vector<float> *jetBtagSub1DCSVbbb_ = new std::vector<float>(0);
	
	float jetDr_(0);

    int decade(0);
    int NN = trIN->GetEntries();
    //NN = 100000;
    std::cout<<"Entries: "<<NN<<std::endl;
	std::vector<float> xRecoAll(0);
	std::vector<float> xPartonAll(0);
	std::vector<float> xParticleAll(0);

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
        y_->clear();
        jetBtagSub0_->clear();
        jetBtagSub1_->clear();
        jetTtag_->clear();
		
       partonPt_->clear();
       partonMass_->clear();
       partonPhi_->clear();
       partonEta_->clear();
       partonY_->clear();
	  
	   jetBtagSub0DCSVbb_->clear();
	   jetBtagSub1DCSVbb_->clear();
	   jetBtagSub0DCSVbbb_->clear();
	   jetBtagSub1DCSVbbb_->clear();
	   
	   xPartonAll.clear();
	   xRecoAll.clear();
	   xParticleAll.clear();

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
						phi_->push_back((*jetPhi)[(*partonMatchIdx)[indexMin]]);
						y_->push_back((*jetY)[(*partonMatchIdx)[indexMin]]);
						jetBtagSub0_->push_back( (*jetBtagSub0)[(*partonMatchIdx)[indexMin]]);
						//jetBtagSub1_->push_back( (*jetBtagSub1)[(*partonMatchIdx)[indexMin]]);
						jetTtag_->push_back( (*jetTtag)[(*partonMatchIdx)[indexMin]]);
						
						jetBtagSub0DCSVbb_->push_back((*jetBtagSub0DCSVbb)[(*partonMatchIdx)[indexMin]]);
						jetBtagSub1DCSVbb_->push_back((*jetBtagSub1DCSVbb)[(*partonMatchIdx)[indexMin]]);
						jetBtagSub0DCSVbbb_->push_back((*jetBtagSub0DCSVbbb)[(*partonMatchIdx)[indexMin]]);
						jetBtagSub1DCSVbbb_->push_back((*jetBtagSub1DCSVbbb)[(*partonMatchIdx)[indexMin]]);

						partonPt_->push_back((*partonPt)[indexMin]);
						partonMass_->push_back((*partonMass)[indexMin]);
						partonPhi_->push_back((*partonPhi)[indexMin]);
						partonEta_->push_back((*partonEta)[indexMin]);
						partonY_->push_back((*partonY)[indexMin]);
						
						
						//cout<<(*partonMatchIdx)[indexMin]<<endl;
						//cout<<"------------"<<endl;
					}
					  
			   }
			   
			}
	if(isMatched > 1)
    {
      if((*pt_)[0] > 1300 && (*pt_)[1] < 1400)
      {
        myFile<<"Before Parton: ";  
        for (int i = 0; i < partonPt->size(); i++)
        {
          myFile<<(*partonPt)[i];
        }
        myFile<<"Before Reco: ";
        for (int i = 0; i < jetPt->size(); i++)
        {
          myFile<<(*jetPt)[i];
          myFile<<" ";
          myFile<<(*partonMatchIdx)[i];
          myFile<<" ";
          myFile<<(*partonMatchDR)[i];
          myFile<<" ";
        }
        
        myFile<<"\nAfter Parton: ";
        for (unsigned int i = 0; i < partonPt_->size(); i++)
        {
          myFile<<(*partonPt_)[i]<<" ";
        }
        myFile<<"\nAfter Reco: ";
        for (unsigned int i = 0; i < pt_->size(); i++)
        {
          myFile<<(*pt_)[i]<<" ";
        }
        
        myFile<<"\n";
      }
      //continue;

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
		xRecoAll.push_back(fabs((*y_)[leadingPt]));
		xRecoAll.push_back(fabs((*y_)[subleadingPt])); 
				
		xPartonAll.push_back(mTTbarParton);
		xPartonAll.push_back(ptTTbarParton);
		xPartonAll.push_back(yTTbarParton);
		xPartonAll.push_back((*partonPt_)[leadingPt]);
		xPartonAll.push_back((*partonPt_)[subleadingPt]);
		xPartonAll.push_back(fabs((*partonY_)[leadingPt]));
		xPartonAll.push_back(fabs((*partonY_)[subleadingPt])); 

		xParticleAll.push_back(mJJGen);
		xParticleAll.push_back(ptJJGen);
		xParticleAll.push_back(yJJGen);
		xParticleAll.push_back((*genjetPt)[leadingPt]);
		xParticleAll.push_back((*genjetPt)[subleadingPt]);	
		xParticleAll.push_back(fabs((*genjetY)[leadingPt]));
		xParticleAll.push_back(fabs((*genjetY)[subleadingPt]));	
      	
	  		

		
	  //---------------------------end of MATCHING---------------------------------------------------------
	  bool recoCuts, partonCuts, particleCuts; 
	  bool massCut = (*mass_)[0] > 120 && (*mass_)[0] < 220 && (*mass_)[1] > 120 && (*mass_)[1] < 220;
	  bool tTaggerCut = (*jetTtag_)[0] > selMvaCut && (*jetTtag_)[1] > selMvaCut;
	  recoCuts = nJets > 1 && fabs((*eta_)[0]) < 2.4 && fabs((*eta_)[1]) <2.4 && (*pt_)[0] > 400 && (*pt_)[1] > 400 && mJJ > 1000 && massCut && nLeptons==0 && (*bit)[triggerFloat];
	  partonCuts = fabs((*partonEta_)[0]) < 2.4 && fabs((*partonEta_)[1]) <2.4 && (*partonPt_)[0] > 400 && (*partonPt_)[1] > 400 && mTTbarParton > 1000;
	  particleCuts = fabs((*genjetEta)[0]) < 2.4 && fabs((*genjetEta)[1]) <2.4 && (*genjetPt)[0] > 400 && (*genjetPt)[1] > 400 && mJJGen > 1000 && nJetsGen >1 &&
	  				 (*genjetMassSoftDrop)[0] > 120 && (*genjetMassSoftDrop)[0] < 220 && (*genjetMassSoftDrop)[1] > 120 && (*genjetMassSoftDrop)[1] < 220;
	  bool deepCSV = (((*jetBtagSub0DCSVbb_)[0] + (*jetBtagSub0DCSVbbb_)[0])> deepCSVFloat || ((*jetBtagSub1DCSVbb_)[0] + (*jetBtagSub1DCSVbbb_)[0])> deepCSVFloat) && 
					 (((*jetBtagSub0DCSVbb_)[1] + (*jetBtagSub0DCSVbbb_)[1])> deepCSVFloat || ((*jetBtagSub1DCSVbb_)[1] + (*jetBtagSub1DCSVbbb_)[1])> deepCSVFloat);
			
      bool btagCut;			
	  btagCut = deepCSV;
	   
	  //qcout<<"----------"<<endl;
	  
		  //fill the denominators
		  //1. denominator passing only reco cuts for topTagger (same for parton and particle)
		  if(recoCuts && btagCut && tTaggerCut) 
		  {
		  	for(int ivar = 0; ivar < NVAR; ivar++)
	  		{
				hReco[f][ivar]->Fill(xRecoAll[ivar], genEvtWeight);				
			}
		  } 
		  //2. fill the histograms pass reco and parton cuts numerators for efficiencies and acceptance
		  //fill also response matrix 
	      if(partonCuts && recoCuts && tTaggerCut && btagCut)
		  {
			  	for(int ivar = 0; ivar < NVAR; ivar++)
	  			{			  	  
				   hPartonReco[f][ivar]->Fill(xPartonAll[ivar], genEvtWeight);
				   hRecoParton[f][ivar]->Fill(xRecoAll[ivar], genEvtWeight);
           //myfile << evtNo << "\n";
				   hPartonResponse[f][ivar]->Fill(xPartonAll[ivar] ,xRecoAll[ivar], genEvtWeight *weights[f]*LUMI);
				}//---- end of the ivar loop
			  	
	      }//----- end of selection cuts parton and reco 

	      //3. Events pass reco and particle 
	      if(particleCuts && recoCuts && tTaggerCut && btagCut)
	      {
	      	for(int ivar = 0; ivar < NVAR; ivar++)
	  		{
	      		hParticleReco[f][ivar]->Fill(xParticleAll[ivar], genEvtWeight);
	      		hRecoParticle[f][ivar]->Fill(xRecoAll[ivar], genEvtWeight);

	      		hParticleResponse[f][ivar]->Fill(xParticleAll[ivar], xRecoAll[ivar], genEvtWeight*weights[f]*LUMI);
	      	}
	      }
	      if(particleCuts) 
	      {
	      	for(int ivar = 0; ivar < NVAR; ivar++)
	  		{
	      		hParticle[f][ivar]->Fill(xParticleAll[ivar], genEvtWeight);
	      	}
	      }
	  
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

	    xPartonAll.push_back(fabs(partonYCnt[0]));
	    xPartonAll.push_back(fabs(partonYCnt[1]));
	  }
	  else
	  {
	  	xPartonAllCnt.push_back(partonPtCnt[1]);
	    xPartonAllCnt.push_back(partonPtCnt[0]);

	    xPartonAll.push_back(fabs(partonYCnt[1]));
	    xPartonAll.push_back(fabs(partonYCnt[0]));
	  }


	  for(int ivar = 0; ivar < NVAR; ivar++)
	  {
	  	if(partonCuts)
			hParton[f][ivar]->Fill(xPartonAllCnt[ivar], genEvtWeightCnt);
	  }
  }

  //--------------------------------------------END OF EVENT COUNTER LOOP ------------------------------------------------------------------
}//----end of file loop
  
  
  for(int ivar =0; ivar<NVAR; ivar++)
  {
	//for every slice
	for(int j=0; j<fileNames.size(); j++)
    {
      hReco[j][ivar]->Scale(weights[j]*LUMI);
      hParticle[j][ivar]->Scale(weights[j]*LUMI);

      hRecoParticle[j][ivar]->Scale(weights[j]*LUMI);
      hParticleReco[j][ivar]->Scale(weights[j]*LUMI);

      hRecoParton[j][ivar]->Scale(weights[j]*LUMI);
      hPartonReco[j][ivar]->Scale(weights[j]*LUMI);
    }
    
    for(int j=1; j<fileNames.size(); j++)
	{
	  //Add them to get the whole phase space
      hReco[0][ivar]->Add(hReco[j][ivar]);
      hParticle[0][ivar]->Add(hParticle[j][ivar]);

      hRecoParticle[0][ivar]->Add(hRecoParticle[j][ivar]);
      hParticleReco[0][ivar]->Add(hParticleReco[j][ivar]);

      hRecoParton[0][ivar]->Add(hRecoParton[j][ivar]);
      hPartonReco[0][ivar]->Add(hPartonReco[j][ivar]);

      hPartonResponse[0][ivar]->Add(hPartonResponse[j][ivar]);
      hParticleResponse[0][ivar]->Add(hParticleResponse[j][ivar]);

	}
  
  hParton[0][ivar]->Scale(weights[0]*LUMI);
  for(int f=1; f<fileNames.size(); f++) 
  {
    hParton[f][ivar]->Scale(weights[f]*LUMI);
	  hParton[0][ivar]->Add(hParton[f][ivar]);
  }

  }
  TEfficiency *efficiency_parton[NVAR], *acceptance_parton[NVAR];
  TEfficiency *efficiency_particle[NVAR], *acceptance_particle[NVAR];
 
  //efficiency for parton quantity and for topTagger (new)


  for(int ivar = 0; ivar< NVAR; ivar++)
  {

  cout<<"--------"<<endl;
  cout<<"parton "<<ivar<<endl;
  efficiency_parton[ivar]  = new TEfficiency(*hPartonReco[0][ivar], *hParton[0][ivar]);
  if (ivar ==0 || ivar ==1 || ivar ==3 || ivar ==4) efficiency_parton[ivar]->SetTitle(TString::Format("EfficiencyParton_%s;%s (GeV);Efficiency Parton",varParton[ivar].Data(),varParton[ivar].Data()));
  else efficiency_parton[ivar]->SetTitle(TString::Format("EfficiencyParton_%s;%s ;Efficiency Parton",varParton[ivar].Data(),varParton[ivar].Data()));
  efficiency_parton[ivar]->SetStatisticOption(TEfficiency::kFNormal);
  efficiency_parton[ivar]->SetUseWeightedEvents();
  efficiency_parton[ivar]->SetLineColor(kRed);
	  
  acceptance_parton[ivar]  = new TEfficiency(*hRecoParton[0][ivar], *hReco[0][ivar]);
  if (ivar ==0 ||ivar ==1 || ivar ==3 || ivar ==4)  acceptance_parton[ivar]->SetTitle(TString::Format("AcceptanceParton_%s; %s(GeV);Acceptance Parton",varReco[ivar].Data(),varReco[ivar].Data()));
  else acceptance_parton[ivar]->SetTitle(TString::Format("AcceptanceParton_%s; %s;Acceptance Parton",varReco[ivar].Data(), varReco[ivar].Data()));
  acceptance_parton[ivar]->SetStatisticOption(TEfficiency::kFNormal);
  acceptance_parton[ivar]->SetUseWeightedEvents();
  acceptance_parton[ivar]->SetLineColor(kRed);


  cout<<"-------"<<endl;
  cout<<"particle "<<ivar<<endl;
  efficiency_particle[ivar]  = new TEfficiency(*hParticleReco[0][ivar], *hParticle[0][ivar]);
  if (ivar ==0 || ivar ==1 || ivar ==3 || ivar ==4) efficiency_particle[ivar]->SetTitle(TString::Format("EfficiencyParticle_%s;%s (GeV);Efficiency Particle",varParticle[ivar].Data(),varParticle[ivar].Data()));
  else efficiency_particle[ivar]->SetTitle(TString::Format("EfficiencyParticle_%s;%s ;Efficiency Particle",varParticle[ivar].Data(),varParticle[ivar].Data()));
  efficiency_particle[ivar]->SetStatisticOption(TEfficiency::kFNormal);
  efficiency_particle[ivar]->SetUseWeightedEvents();
  efficiency_particle[ivar]->SetLineColor(kRed);
	  
  acceptance_particle[ivar]  = new TEfficiency(*hRecoParticle[0][ivar], *hReco[0][ivar]);
  if (ivar ==0 ||ivar ==1 || ivar ==3 || ivar ==4)  acceptance_particle[ivar]->SetTitle(TString::Format("AcceptanceParticle_%s; %s(GeV);Acceptance particle",varReco[ivar].Data(),varReco[ivar].Data()));
  else acceptance_particle[ivar]->SetTitle(TString::Format("AcceptanceParticle_%s; %s;Acceptance particle",varReco[ivar].Data(), varReco[ivar].Data()));
  acceptance_particle[ivar]->SetStatisticOption(TEfficiency::kFNormal);
  acceptance_particle[ivar]->SetUseWeightedEvents();
  acceptance_particle[ivar]->SetLineColor(kRed);


  }
  //WARNING !!!
  //for purity and stability: when we have unequal binning we cannot find purity and stability but only when we have a
  //square matrx 


  TFile *outFile;
  TString nominal ="";
  if(isNominalMC) nominal = "NominalMC";
  outFile = TFile::Open(TString::Format("%s/UnequalBins/ResponsesEfficiency%s_%s.root", year.Data(),nominal.Data(),year.Data()), "RECREATE");
  //outFile->cd();
  //write them to file
  for(int ivar = 0; ivar<NVAR; ivar++)
  {
  	efficiency_parton[ivar]->Write(TString::Format("EfficiencyParton_%s",varParton[ivar].Data()));
  	efficiency_particle[ivar]->Write(TString::Format("EfficiencyParticle_%s",varParticle[ivar].Data()));

  	acceptance_parton[ivar]->Write(TString::Format("AcceptanceParton_%s",varReco[ivar].Data()));
  	acceptance_particle[ivar]->Write(TString::Format("AcceptanceParticle_%s",varReco[ivar].Data()));

  	hPartonResponse[0][ivar]->Write(TString::Format("hPartonResponse_%s", varReco[ivar].Data()));
  	hParticleResponse[0][ivar]->Write(TString::Format("hParticleResponse_%s", varReco[ivar].Data()));

  }//end of ivar
  myFile.close();
  outFile->Close();
  fileNames.clear();
  histoNames.clear();
  XSEC.clear();
 }
