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

void FillHistograms_Reduced_bTag(TString file_name, TString ttbar_process, TString year = "2016", TString variation = "up")
{
  globalYear = year;
  initFilesMapping();
  cout<<"YEAR: "<<year<<endl;
  cout<<"file_name: "<<file_name<<endl;
  cout<<"ttbar_process: "<<ttbar_process<<endl;
  float triggerFloat;
  if(year.EqualTo("2016")) triggerFloat = 2;
  else triggerFloat = 5;

  float deepCSVFloat = floatConstants[TString::Format("btagWP%s",year.Data())];
  float selMvaCut = topTaggerConstants[TString::Format("topTagger%s",year.Data())];
  float LUMI = luminosity[TString::Format("luminosity%s", year.Data())];
  float LUMI_CR = luminosityCR[TString::Format("luminosity%s", year.Data())];
  float XSEC = XSECAll[year.Data()][ttbar_process.Data()];
  cout<<"XSEC: "<<XSEC<<endl;
  const int NVAR = 12;
  std::vector< std::vector <Float_t> > const BND = {{1000, 1200, 1400, 1600, 1800, 2000, 2400, 3000, 5000},//mJJ
                       {0, 60, 150, 300, 450, 600, 850, 1100, 1300}, //ptJJ
                       {-2.4, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.4}, //yjj
                       {400, 450, 500, 570, 650, 800, 1000, 1250, 1500}, //jetPt0
                       {400, 450, 500, 570, 650, 800, 1000, 1250, 1500}, //jetPt1
                       {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.1, 2.4}, //jetY0
                       {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.1, 2.4}, //jetY1
                       {1,2,3,4,5,6,7,8,9,10,13,16}, //chi
                       {-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1}, //|cosTheta*| leading
                       {-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1}}; //|cosTheta*| subleading

  int NBINS[BND.size()];
  for (int i = 0; i<BND.size()+2; i++)
  {
    if(i<10) NBINS[i] = BND[i].size()-1;
    else NBINS[i] = 100;
  }
  TString varReco[NVAR]   = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1","chi","cosTheta_0", "cosTheta_1","mTop_Leading", "mTop_Subleading"};

  float weights;
  TH1F *hReco[NVAR], *hRecoCR[NVAR];
  int maxJetMass = 220;
  int minJetMass = 120;
  	//declare the histograms
  	for(int ivar =0; ivar<NVAR; ivar++)
  	{
      int sizeBins = NBINS[ivar];
      if(ivar>=10)
      {
        hReco[ivar] = new TH1F(TString::Format("hReco_%s", varReco[ivar].Data()), TString::Format("hReco_%s",varReco[ivar].Data()), sizeBins, minJetMass,maxJetMass);
        hRecoCR[ivar] = new TH1F(TString::Format("hRecoCR_%s", varReco[ivar].Data()), TString::Format("hRecoCR_%s",varReco[ivar].Data()), sizeBins, minJetMass,maxJetMass);
      }
      else
      {
        int sizeBins = NBINS[ivar];
        float tempBND[NBINS[ivar]+1];
        std::copy(BND[ivar].begin(), BND[ivar].end(), tempBND);
        hReco[ivar] = new TH1F(TString::Format("hReco_%s", varReco[ivar].Data()), TString::Format("hReco_%s",varReco[ivar].Data()), sizeBins, tempBND);
        hRecoCR[ivar] = new TH1F(TString::Format("hRecoCR_%s", varReco[ivar].Data()), TString::Format("hRecoCR_%s",varReco[ivar].Data()), sizeBins, tempBND);
      }
      hReco[ivar]->Sumw2();
      hRecoCR[ivar]->Sumw2();
    }
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


    TFile *file = TFile::Open(file_name.Data());

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
    if(variation.EqualTo("up")) trIN->SetBranchAddress("bTagEvntWeight_up", &bTagEvntWeight);
    else if(variation.EqualTo("down")) trIN->SetBranchAddress("bTagEvntWeight_down", &bTagEvntWeight);
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



    float norm = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
	float weight = XSEC/norm;
	weights = weight;

	//for parton matching
	std::vector<int> *jetMatchedIndexes = new std::vector<int>(0);
	std::vector<float> *jetMatchedDr = new std::vector<float>(0);
	std::vector<float> *eta_ = new std::vector<float>(0);
  std::vector<float> *y_ = new std::vector<float>(0);
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
        y_->clear();
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
            y_->push_back((*jetY)[(*partonMatchIdx)[indexMin]]);
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

		TLorentzVector p4T[2], p4T_ZMF[2], p4TTbar;
    	p4T[leadingPt].SetPtEtaPhiM((*pt_)[leadingPt], (*eta_)[leadingPt], (*phi_)[leadingPt], (*mass_)[leadingPt]);
   		p4T[subleadingPt].SetPtEtaPhiM((*pt_)[subleadingPt], (*eta_)[subleadingPt], (*phi_)[subleadingPt], (*mass_)[subleadingPt]);

  	    TVector3 ttbarBoostVector = getBoostVector(p4T[leadingPt], p4T[subleadingPt], p4TTbar);

    	p4T_ZMF[0].SetPtEtaPhiM(p4T[leadingPt].Pt(), p4T[leadingPt].Eta(), p4T[leadingPt].Phi(), p4T[leadingPt].M());
    	p4T_ZMF[1].SetPtEtaPhiM(p4T[subleadingPt].Pt(), p4T[subleadingPt].Eta(), p4T[subleadingPt].Phi(), p4T[subleadingPt].M());
     	p4T_ZMF[0].Boost(ttbarBoostVector);
     	p4T_ZMF[1].Boost(ttbarBoostVector);

	    float yStarExp  = TMath::Exp(fabs(p4T_ZMF[0].Rapidity() - p4T_ZMF[1].Rapidity())); //this is chi = e^(|y*|) , y* = 1/2(y1-y0)

      xRecoAll.push_back(mJJ);
      xRecoAll.push_back(ptJJ);
      xRecoAll.push_back(yJJ);
      xRecoAll.push_back((*pt_)[leadingPt]);
      xRecoAll.push_back((*pt_)[subleadingPt]);
      xRecoAll.push_back(fabs((*y_)[leadingPt]));
      xRecoAll.push_back(fabs((*y_)[subleadingPt]));
      xRecoAll.push_back(yStarExp); //this is chi
     	xRecoAll.push_back(TMath::Cos(p4T_ZMF[0].Theta())); //this is |cos(theta*)| leading
     	xRecoAll.push_back(TMath::Cos(p4T_ZMF[1].Theta())); //this is |cos(theta*)| subleading
      xRecoAll.push_back((*mass_)[leadingPt]);
      xRecoAll.push_back((*mass_)[subleadingPt]);


	  //---------------------------end of MATCHING---------------------------------------------------------
	  bool recoCuts;
	  bool massCut = (*mass_)[0] > 120 && (*mass_)[0] < 220 && (*mass_)[1] > 120 && (*mass_)[1] < 220;
	  bool tTaggerCut = (*jetTtag_)[0] > selMvaCut && (*jetTtag_)[1] > selMvaCut;
	  recoCuts = nJets > 1 && fabs((*eta_)[0]) < 2.4 && fabs((*eta_)[1]) <2.4 && (*pt_)[0] > 400 && (*pt_)[1] > 400 && mJJ > 1000 && massCut && nLeptons==0 && (*bit)[triggerFloat];
	  bool deepCSV = (((*jetBtagSub0DCSVbb_)[0] + (*jetBtagSub0DCSVbbb_)[0])> deepCSVFloat || ((*jetBtagSub1DCSVbb_)[0] + (*jetBtagSub1DCSVbbb_)[0])> deepCSVFloat) &&
					 (((*jetBtagSub0DCSVbb_)[1] + (*jetBtagSub0DCSVbbb_)[1])> deepCSVFloat || ((*jetBtagSub1DCSVbb_)[1] + (*jetBtagSub1DCSVbbb_)[1])> deepCSVFloat);

      bool btagCut;
	  btagCut = deepCSV;
    //Signal Region 2btags
		if(recoCuts && btagCut && tTaggerCut)
		{
		  for(int ivar = 0; ivar < NVAR; ivar++)
	  	{
		    hReco[ivar]->Fill(xRecoAll[ivar], genEvtWeight*bTagEvntWeight);
		  }
	  }
    //Control Region 0btag
    if(recoCuts && !btagCut && tTaggerCut)
	  {
	  	for(int ivar = 0; ivar < NVAR; ivar++)
  		{
		    hRecoCR[ivar]->Fill(xRecoAll[ivar], genEvtWeight*bTagEvntWeight);
		  }
	  }

	 }//----- end of is matched
  }//---end the event for


    for(int ivar =0; ivar<NVAR; ivar++)
    {
      hReco[ivar]->Scale(weights*LUMI);
      hRecoCR[ivar]->Scale(weights*LUMI_CR);
    }//end of loop on all vars


    TFile *outFile;
    outFile = TFile::Open(TString::Format("%s/bTagVariation/HistoReduced_%s_%s.root", year.Data(), ttbar_process.Data(), variation.Data()), "RECREATE");
    //outFile->cd();
    //write them to file
    for(int ivar = 0; ivar<NVAR; ivar++)
    {
      hReco[ivar]->Write(TString::Format("hWt_%s_2btag", varReco[ivar].Data()));
      hRecoCR[ivar]->Write(TString::Format("hWt_%s_0btag", varReco[ivar].Data()));
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
