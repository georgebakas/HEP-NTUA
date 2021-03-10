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


std::vector<TString> histoNames;
std::vector<TString> fileNames;

#include "TemplateConstants.h"
TVector3 getBoostVector(TLorentzVector p4_1, TLorentzVector p4_2, TLorentzVector &p4CombinedVector);

TString globalYear;

void FillHistograms_massWindow_DY(TString file_name, TString mass_name, TString year = "2016")
{
  globalYear = year;
  initFilesMapping(false);
  cout<<"YEAR: "<<year<<endl;
  cout<<"file_name: "<<file_name<<endl;
  cout<<"mass_name: "<<mass_name<<endl;
  int triggerFloat;
  if(year.Contains("2016")) triggerFloat = 2;
  else triggerFloat = 5;

  float deepCSVFloat = floatConstants[TString::Format("btagWP%s",year.Data())];
  float selMvaCut = topTaggerConstants[TString::Format("topTagger%s",year.Data())];
  float LUMI = luminosity[TString::Format("luminosity%s", year.Data())];
  float XSEC = XSECAll[year.Data()][mass_name.Data()];
  cout<<"XSEC: "<<XSEC<<endl;
  std::vector< std::vector <Float_t> > const BND = {{1,2,3,4,5,6,7,8,9,10,13,16}, //chi
                                                    {-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1}, //|cosTheta*| leading
                                                    {-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1},//|cosTheta*| subleading
                                                    {1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,1700,1800,1900,2000,2200,2600,3000,5000}}; //mJJ

  int NBINS[BND.size()];
  const int NVAR = 4;
  for (int i = 0; i<BND.size(); i++) NBINS[i] = BND[i].size()-1;
  TString varReco[NVAR]   = {"chi", "cosTheta_0", "cosTheta_1", "mJJ"};
  TString varParton[NVAR] = {"chiParton", "cosThetaParton_0", "cosThetaParton_1", "mJJGen"};
  TString varParticle[NVAR] = {"chiParticle", "cosThetaParticle_0", "cosThetaParticle_1", "mTTbarParton"};

  const int NWINDOWS=6;
  int massWindows[NWINDOWS] = {1000,1200,1400,1600,1800,2000};

  float weights;
  TH1F *hParton_Lo1[NVAR][NWINDOWS], *hParticle_Lo1[NVAR][NWINDOWS], *hReco_Lo1[NVAR][NWINDOWS];
  TH1F *hParton_Gr1[NVAR][NWINDOWS], *hParticle_Gr1[NVAR][NWINDOWS], *hReco_Gr1[NVAR][NWINDOWS];

  	//declare the histograms
  for(int ivar =0; ivar<NVAR; ivar++)
  {
    for(int iwind =0; iwind<NWINDOWS; iwind++)
    {
      int sizeBins = NBINS[ivar];
      float tempBND[NBINS[ivar]+1];
      std::copy(BND[ivar].begin(), BND[ivar].end(), tempBND);
  		hParton_Lo1[ivar][iwind] = new TH1F(TString::Format("hParton_Lo1_%s_%d", varParton[ivar].Data(), massWindows[iwind]), TString::Format("hParton_%s_%d",varParton[ivar].Data(),massWindows[iwind]), sizeBins, tempBND);
      hReco_Lo1[ivar][iwind] = new TH1F(TString::Format("hReco_Lo1_%s_%d", varReco[ivar].Data(),massWindows[iwind]), TString::Format("hReco_%s_%d",varReco[ivar].Data(),massWindows[iwind]), sizeBins, tempBND);
      hParticle_Lo1[ivar][iwind] = new TH1F(TString::Format("hParticle_Lo1_%s_%d", varParticle[ivar].Data(),massWindows[iwind]), TString::Format("hParticle_%s_%d", varParticle[ivar].Data(),massWindows[iwind]), sizeBins, tempBND);

      hParton_Lo1[ivar][iwind]->Sumw2();
      hReco_Lo1[ivar][iwind]->Sumw2();
      hParticle_Lo1[ivar][iwind]->Sumw2();

      hParton_Gr1[ivar][iwind] = new TH1F(TString::Format("hParton_Gr1_%s_%d", varParton[ivar].Data(), massWindows[iwind]), TString::Format("hParton_%s_%d",varParton[ivar].Data(),massWindows[iwind]), sizeBins, tempBND);
      hReco_Gr1[ivar][iwind] = new TH1F(TString::Format("hReco_Gr1_%s_%d", varReco[ivar].Data(),massWindows[iwind]), TString::Format("hReco_%s_%d",varReco[ivar].Data(),massWindows[iwind]), sizeBins, tempBND);
      hParticle_Gr1[ivar][iwind] = new TH1F(TString::Format("hParticle_Gr1_%s_%d", varParticle[ivar].Data(),massWindows[iwind]), TString::Format("hParticle_%s_%d", varParticle[ivar].Data(),massWindows[iwind]), sizeBins, tempBND);

      hParton_Gr1[ivar][iwind]->Sumw2();
      hReco_Gr1[ivar][iwind]->Sumw2();
      hParticle_Gr1[ivar][iwind]->Sumw2();
    }

  }//end of ivar loop
    int nJets,nLeptons, category(0);
    vector<bool>  *bit(0),*matchedJet(0);
    //reco vars:
    std::vector<float> *jetPt(0), *jetY(0), *jetEta(0), *jetPhi(0), *jetTtag(0);
    float genEvtWeight(0);
    double  bTagEvntWeight(0);
    float mJJ(0), ptJJ(0), yJJ(0),mva(0);
    vector<float> *tau3(0),*tau2(0),*tau1(0);
	  vector<float> *jetMassSub0(0), *jetMassSub1(0), *jetBtagSub0(0), *jetBtagSub1(0);
    vector<float> *jetMassSoftDrop(0), *partonEta(0);

    //parton
    std::vector<float> *partonPt(0), *partonPhi(0), *partonMass(0),*partonMatchDR(0);
    float yTTbarParton(0), ptTTbarParton(0), mTTbarParton(0);

    //float yTopParton[2], ptTopParton[2];
    std::vector<int> *addedIndexes = new std::vector<int>(0);
    vector<int> *partonId(0), *partonMatchIdx(0);

    std::vector<float> *jetBtagSub0DCSVbb(0), *jetBtagSub1DCSVbb(0);
    std::vector<float> *jetBtagSub0DCSVbbb(0), *jetBtagSub1DCSVbbb(0);

    //particle
    std::vector<float> *genjetPt(0), *genjetY(0), *genjetEta(0), *genSoftDropMass(0), *genjetMassSoftDrop(0),*genjetPhi(0);
    std::vector<int> *nSubGenJets(0);
    int nJetsGen(0);
    float mJJGen(0), ptJJGen(0), yJJGen(0);

    TFile *file = TFile::Open(eospath[year.Data()]+file_name.Data());

    TTree *trIN = (TTree*)file->Get("boosted/events");

	//------- input tree --------------
    trIN->SetBranchAddress("nJets"          ,&nJets);
    trIN->SetBranchAddress("nLeptons"       ,&nLeptons);
    trIN->SetBranchAddress("jetPt"          ,&jetPt);
    trIN->SetBranchAddress("jetY"           ,&jetY);
    trIN->SetBranchAddress("jetEta"         ,&jetEta);
    trIN->SetBranchAddress("jetPhi"         ,&jetPhi);
    trIN->SetBranchAddress("jetTau3"        ,&tau3);
    trIN->SetBranchAddress("jetTau2"        ,&tau2);
    trIN->SetBranchAddress("jetTau1"        ,&tau1);
    trIN->SetBranchAddress("triggerBit"     ,&bit);
    trIN->SetBranchAddress("genEvtWeight"   ,&genEvtWeight);
    trIN->SetBranchAddress("bTagEvntWeight"   ,&bTagEvntWeight);
    trIN->SetBranchAddress("jetMassSub0"    ,&jetMassSub0);
    trIN->SetBranchAddress("jetMassSub1"    ,&jetMassSub1);
    trIN->SetBranchAddress("mJJ"   			,&mJJ);
    trIN->SetBranchAddress("yJJ"   			,&yJJ);
    trIN->SetBranchAddress("ptJJ"   		,&ptJJ);
	trIN->SetBranchAddress("jetBtagSub0"	,&jetBtagSub0);
    //trIN->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1);
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

    //particle
    trIN->SetBranchAddress("mJJGen"			,&mJJGen);
    trIN->SetBranchAddress("ptJJGen"		,&ptJJGen);
    trIN->SetBranchAddress("yJJGen"			,&yJJGen);
    trIN->SetBranchAddress("genjetPt"		,&genjetPt);
    trIN->SetBranchAddress("genjetY"		,&genjetY);
    trIN->SetBranchAddress("genjetPhi"		,&genjetPhi);
    trIN->SetBranchAddress("genjetEta"		,&genjetEta);
    trIN->SetBranchAddress("nJetsGen"		,&nJetsGen);
    trIN->SetBranchAddress("genjetMassSoftDrop", &genjetMassSoftDrop);


    float norm = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
	float weight = XSEC/norm;
	weights = weight;

	//for parton matching
	std::vector<int> *jetMatchedIndexes = new std::vector<int>(0);
	std::vector<float> *jetMatchedDr = new std::vector<float>(0);
	std::vector<float> *eta_ = new std::vector<float>(0);
	std::vector<float> *phi_ = new std::vector<float>(0);
	std::vector<float> *mass_ = new std::vector<float>(0);
	std::vector<float> *pt_ = new std::vector<float>(0);
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
      /*
		int isMatched =0;
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
    */
	   xPartonAll.clear();
	   xRecoAll.clear();
	   xParticleAll.clear();

      //----------------------MATCHING------------------------------------------------------
    /*
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
						//jetBtagSub1_->push_back( (*jetBtagSub1)[(*partonMatchIdx)[indexMin]]);
						jetTtag_->push_back( (*jetTtag)[(*partonMatchIdx)[indexMin]]);

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

			}*/
	if(nJets > 1)
    {
   	   	int leadingPt =0;
      	int subleadingPt = 1;

  		if((*jetPt)[0] < (*jetPt)[1])
  		{
   		    subleadingPt =0;
   		    leadingPt = 1;
   		}

		TLorentzVector p4T[2], p4T_ZMF[2], p4TTbar;
    	p4T[leadingPt].SetPtEtaPhiM((*jetPt)[leadingPt], (*jetEta)[leadingPt], (*jetPhi)[leadingPt], (*jetMassSoftDrop)[leadingPt]);
   		p4T[subleadingPt].SetPtEtaPhiM((*jetPt)[subleadingPt], (*jetEta)[subleadingPt], (*jetPhi)[subleadingPt], (*jetMassSoftDrop)[subleadingPt]);

  	    TVector3 ttbarBoostVector = getBoostVector(p4T[leadingPt], p4T[subleadingPt], p4TTbar);

    	p4T_ZMF[0].SetPtEtaPhiM(p4T[leadingPt].Pt(), p4T[leadingPt].Eta(), p4T[leadingPt].Phi(), p4T[leadingPt].M());
    	p4T_ZMF[1].SetPtEtaPhiM(p4T[subleadingPt].Pt(), p4T[subleadingPt].Eta(), p4T[subleadingPt].Phi(), p4T[subleadingPt].M());
     	p4T_ZMF[0].Boost(ttbarBoostVector);
     	p4T_ZMF[1].Boost(ttbarBoostVector);

	    float yStarExp  = TMath::Exp(fabs(p4T_ZMF[0].Rapidity() - p4T_ZMF[1].Rapidity())); //this is chi = e^(|y*|) , y* = 1/2(y1-y0)

     	xRecoAll.push_back(yStarExp); //this is chi
     	xRecoAll.push_back(TMath::Cos(p4T_ZMF[0].Theta())); //this is |cos(theta*)| leading
     	xRecoAll.push_back(TMath::Cos(p4T_ZMF[1].Theta())); //this is |cos(theta*)| subleading
      xRecoAll.push_back(mJJ); //this is mTTbarParton

		//now parton
		TLorentzVector p4TParton[2], p4T_ZMFParton[2], p4TTbarParton;
    	p4TParton[leadingPt].SetPtEtaPhiM((*partonPt)[leadingPt], (*partonEta)[leadingPt], (*partonPhi)[leadingPt], (*partonMass)[leadingPt]);
   		p4TParton[subleadingPt].SetPtEtaPhiM((*partonPt)[subleadingPt], (*partonEta)[subleadingPt], (*partonPhi)[subleadingPt], (*partonMass)[subleadingPt]);

  	    TVector3 ttbarBoostVectorParton = getBoostVector(p4TParton[leadingPt], p4TParton[subleadingPt], p4TTbarParton);

    	p4T_ZMFParton[0].SetPtEtaPhiM(p4TParton[leadingPt].Pt(), p4TParton[leadingPt].Eta(), p4TParton[leadingPt].Phi(), p4TParton[leadingPt].M());
    	p4T_ZMFParton[1].SetPtEtaPhiM(p4TParton[subleadingPt].Pt(), p4TParton[subleadingPt].Eta(), p4TParton[subleadingPt].Phi(), p4TParton[subleadingPt].M());
     	p4T_ZMFParton[0].Boost(ttbarBoostVectorParton);
     	p4T_ZMFParton[1].Boost(ttbarBoostVectorParton);

	    float yStarExpParton  = TMath::Exp(fabs(p4T_ZMFParton[0].Rapidity() - p4T_ZMFParton[1].Rapidity())); //this is chi = e^(|y*|) , y* = 1/2(y1-y0)

		xPartonAll.push_back(yStarExpParton);
		xPartonAll.push_back(TMath::Cos(p4T_ZMFParton[0].Theta())); //this is |cos(theta*)| leading
		xPartonAll.push_back(TMath::Cos(p4T_ZMFParton[1].Theta())); //this is |cos(theta*)| subleading
    xPartonAll.push_back(mTTbarParton); //this is mTTbarParton

		//now particle
		TLorentzVector p4TParticle[2], p4T_ZMFParticle[2], p4TTbarParticle;
    	p4TParticle[leadingPt].SetPtEtaPhiM((*genjetPt)[leadingPt], (*genjetEta)[leadingPt], (*genjetPhi)[leadingPt], (*genjetMassSoftDrop)[leadingPt]);
   		p4TParticle[subleadingPt].SetPtEtaPhiM((*genjetPt)[subleadingPt], (*genjetEta)[subleadingPt], (*genjetPhi)[subleadingPt], (*genjetMassSoftDrop)[subleadingPt]);

  	    TVector3 ttbarBoostVectorParticle = getBoostVector(p4TParticle[leadingPt], p4TParticle[subleadingPt], p4TTbarParticle);

    	p4T_ZMFParticle[0].SetPtEtaPhiM(p4TParticle[leadingPt].Pt(), p4TParticle[leadingPt].Eta(), p4TParticle[leadingPt].Phi(), p4TParticle[leadingPt].M());
    	p4T_ZMFParticle[1].SetPtEtaPhiM(p4TParticle[subleadingPt].Pt(), p4TParticle[subleadingPt].Eta(), p4TParticle[subleadingPt].Phi(), p4TParticle[subleadingPt].M());
     	p4T_ZMFParticle[0].Boost(ttbarBoostVectorParticle);
     	p4T_ZMFParticle[1].Boost(ttbarBoostVectorParticle);

	    float yStarExpParticle  = TMath::Exp(fabs(p4T_ZMFParticle[0].Rapidity() - p4T_ZMFParticle[1].Rapidity())); //this is chi = e^(|y*|) , y* = 1/2(y1-y0)

		xParticleAll.push_back(yStarExpParticle);
		xParticleAll.push_back(TMath::Cos(p4T_ZMFParticle[0].Theta())); //this is |cos(theta*)| leading
		xParticleAll.push_back(TMath::Cos(p4T_ZMFParticle[1].Theta())); //this is |cos(theta*)| subleading
    xParticleAll.push_back(mJJGen); //this is particle mJJ

	  //---------------------------end of MATCHING---------------------------------------------------------
	  bool recoCuts, partonCuts, particleCuts;
	  bool massCut = (*jetMassSoftDrop)[0] > 120 && (*jetMassSoftDrop)[0] < 220 && (*jetMassSoftDrop)[1] > 120 && (*jetMassSoftDrop)[1] < 220;
	  bool tTaggerCut = (*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut;
	  recoCuts = nJets > 1 && fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1]) <2.4 && (*jetPt)[0] > 400 && (*jetPt)[1] > 400 && mJJ > 1000 && massCut && nLeptons==0 && (*bit)[triggerFloat];
	  partonCuts = fabs((*partonEta)[0]) < 2.4 && fabs((*partonEta)[1]) <2.4 && (*partonPt)[0] > 400 && (*partonPt)[1] > 400 && mTTbarParton > 1000;
	  particleCuts = fabs((*genjetEta)[0]) < 2.4 && fabs((*genjetEta)[1]) && (*genjetPt)[0] > 400 && (*genjetPt)[1] > 400 && mJJGen > 1000 && nJetsGen >1 &&
	  				 (*genjetMassSoftDrop)[0] > 120 && (*genjetMassSoftDrop)[0] < 220 && (*genjetMassSoftDrop)[1] > 120 && (*genjetMassSoftDrop)[1] < 220;
	  bool deepCSV = (((*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0])> deepCSVFloat) &&
					 (((*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1])> deepCSVFloat);

      bool btagCut;
	  btagCut = deepCSV;
	  //qcout<<"----------"<<endl;
    bool massWindowCut;
    float delta_y = fabs((*jetY)[0] - (*jetY)[1]);
		  //fill the denominators
		  //1. denominator passing only reco cuts for topTagger (same for parton and particle)
		  if(recoCuts && btagCut && tTaggerCut)
		  {
        for(int iwind =0; iwind<NWINDOWS; iwind++)
        {
          massWindowCut = mJJ > massWindows[iwind];
          if(massWindowCut)
          {
            if (delta_y > 1)
            {
      		  	for(int ivar = 0; ivar < NVAR; ivar++)
      	  		{
                hReco_Gr1[ivar][iwind]->Fill(xRecoAll[ivar], genEvtWeight*bTagEvntWeight);
      			  }
            }
            else
            {
      		  	for(int ivar = 0; ivar < NVAR; ivar++)
      	  		{
                hReco_Lo1[ivar][iwind]->Fill(xRecoAll[ivar], genEvtWeight*bTagEvntWeight);
      			  }
            }
          }
        }
		  }
      //2. events pass particle cuts
	    if(particleCuts)
	    {
        for(int iwind =0; iwind<NWINDOWS; iwind++)
        {
          massWindowCut = mJJGen > massWindows[iwind];
          if(massWindowCut)
          {
            if (delta_y > 1)
            {
      	      for(int ivar = 0; ivar < NVAR; ivar++)
      	  	  {
      	      	hParticle_Gr1[ivar][iwind]->Fill(xParticleAll[ivar], genEvtWeight*bTagEvntWeight);
      	      }
            }
            else
            {
              for(int ivar = 0; ivar < NVAR; ivar++)
      	  	  {
      	      	hParticle_Lo1[ivar][iwind]->Fill(xParticleAll[ivar], genEvtWeight*bTagEvntWeight);
      	      }
            }
          }
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
  float partonPhiCnt[2], partonMCnt[2];
  //tree for eventCounter
  trCnt->SetBranchAddress("ptTopParton"    ,&partonPtCnt);
  trCnt->SetBranchAddress("phiTopParton"   ,&partonPhiCnt);
  trCnt->SetBranchAddress("mTopParton"     ,&partonMCnt );
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
	  int leadingPt = 0;
	  int subleadingPt = 1;
	  if(partonPtCnt[0] < partonPtCnt[1])
	  {
	  	leadingPt = 1;
	  	subleadingPt = 0;
	  }
	  TLorentzVector p4TParton[2], p4T_ZMFPartonCnt[2], p4TTbarParton;
      p4TParton[leadingPt].SetPtEtaPhiM(partonPtCnt[leadingPt], partonEtaCnt[leadingPt], partonPhiCnt[leadingPt], partonMCnt[leadingPt]);
   	  p4TParton[subleadingPt].SetPtEtaPhiM(partonPtCnt[subleadingPt], partonEtaCnt[subleadingPt], partonPhiCnt[subleadingPt], partonMCnt[subleadingPt]);

  	  TVector3 ttbarBoostVectorParton = getBoostVector(p4TParton[leadingPt], p4TParton[subleadingPt], p4TTbarParton);

      p4T_ZMFPartonCnt[0].SetPtEtaPhiM(p4TParton[leadingPt].Pt(), p4TParton[leadingPt].Eta(), p4TParton[leadingPt].Phi(), p4TParton[leadingPt].M());
      p4T_ZMFPartonCnt[1].SetPtEtaPhiM(p4TParton[subleadingPt].Pt(), p4TParton[subleadingPt].Eta(), p4TParton[subleadingPt].Phi(), p4TParton[subleadingPt].M());
      p4T_ZMFPartonCnt[0].Boost(ttbarBoostVectorParton);
      p4T_ZMFPartonCnt[1].Boost(ttbarBoostVectorParton);

	  float yStarExpPartonCnt  = TMath::Exp(fabs(p4T_ZMFPartonCnt[0].Rapidity() - p4T_ZMFPartonCnt[1].Rapidity())); //this is chi = e^(|y*|) , y* = 1/2(y1-y0)

	  xPartonAllCnt.push_back(yStarExpPartonCnt);
	  xPartonAllCnt.push_back(TMath::Cos(p4T_ZMFPartonCnt[0].Theta())); //this is |cos(theta*)| leading
	  xPartonAllCnt.push_back(TMath::Cos(p4T_ZMFPartonCnt[1].Theta())); //this is |cos(theta*)| subleading
    xPartonAllCnt.push_back(mTTbarPartonCnt); //mTTbarPartonCnt

    float delta_y_cnt = fabs(partonPhiCnt[leadingPt] - partonPhiCnt[subleadingPt]);
    bool massWindowCut;
	  for(int ivar = 0; ivar < NVAR; ivar++)
	  {
      for(int iwind =0; iwind<NWINDOWS; iwind++)
      {
        massWindowCut = mTTbarPartonCnt > massWindows[iwind];
        if(massWindowCut && partonCuts)
        {
          if (delta_y_cnt > 1)
            hParton_Gr1[ivar][iwind]->Fill(xPartonAllCnt[ivar], genEvtWeightCnt);
          else
            hParton_Lo1[ivar][iwind]->Fill(xPartonAllCnt[ivar], genEvtWeightCnt);
        }

      }
	  }
  }

  //--------------------------------------------END OF EVENT COUNTER LOOP ------------------------------------------------------------------


  for(int ivar =0; ivar<NVAR; ivar++)
  {
    for(int iwind= 0; iwind<NWINDOWS; iwind++)
    {
      hReco_Lo1[ivar][iwind]->Scale(weights*LUMI);
      hParticle_Lo1[ivar][iwind]->Scale(weights*LUMI);
      hParton_Lo1[ivar][iwind]->Scale(weights*LUMI);

      hReco_Gr1[ivar][iwind]->Scale(weights*LUMI);
      hParticle_Gr1[ivar][iwind]->Scale(weights*LUMI);
      hParton_Gr1[ivar][iwind]->Scale(weights*LUMI);
    }
  }//end of loop on all vars


  TFile *outFile;
  outFile = TFile::Open(TString::Format("%s/HistoMassWindows_%s", year.Data(),file_name.Data()), "RECREATE");
  //outFile->cd();
  //write them to file
  for(int ivar = 0; ivar<NVAR; ivar++)
  {
    for(int iwind= 0; iwind<NWINDOWS; iwind++)
    {
      hParton_Lo1[ivar][iwind]->Write(TString::Format("hParton_%s_%d", varReco[ivar].Data(), massWindows[iwind]));
      hParticle_Lo1[ivar][iwind]->Write(TString::Format("hParticle_%s_%d", varReco[ivar].Data(),massWindows[iwind]));
      hReco_Lo1[ivar][iwind]->Write(TString::Format("hReco_%s_%d", varReco[ivar].Data(),massWindows[iwind]));

      hParton_Gr1[ivar][iwind]->Write(TString::Format("hParton_%s_%d", varReco[ivar].Data(), massWindows[iwind]));
      hParticle_Gr1[ivar][iwind]->Write(TString::Format("hParticle_%s_%d", varReco[ivar].Data(),massWindows[iwind]));
      hReco_Gr1[ivar][iwind]->Write(TString::Format("hReco_%s_%d", varReco[ivar].Data(),massWindows[iwind]));
    }
  }//end of ivar

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
