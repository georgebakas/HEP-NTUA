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


#include "TemplateConstants.h"
TVector3 getBoostVector(TLorentzVector p4_1, TLorentzVector p4_2, TLorentzVector &p4CombinedVector);

TString globalYear;
bool globalIsNominalMC;


void GetHistograms_bTag(TString file_name, TString ttbar_process, TString year = "2016",  TString variation = "up" ,float mJJCut=1000)
{
  globalYear = year;
  initFilesMapping();
  cout<<"YEAR: "<<year<<endl;
  cout<<"file_name: "<<file_name<<endl;
  cout<<"ttbar_process: "<<ttbar_process<<endl;
  float triggerFloat;
  if(year.Contains("2016")) triggerFloat = 2;
  else triggerFloat = 5;

  float deepCSVFloat = floatConstants[TString::Format("btagWP%s",year.Data())];
  float selMvaCut = topTaggerConstants[TString::Format("topTagger%s",year.Data())];
  float LUMI = luminosity[TString::Format("luminosity%s", year.Data())];
  float LUMI_CR = luminosityCR[TString::Format("luminosity%s", year.Data())];
  float XSEC = XSECAll[year.Data()][ttbar_process.Data()];
  cout<<"XSEC: "<<XSEC<<endl;
  const int NVAR = 10;
  std::vector< std::vector <Float_t> > const BND_gen = {{1000, 1200, 1400, 1600, 1800, 2000, 2400, 3000, 5000}, //mJJ
                                                        {0, 60, 150, 300, 450, 850, 1300}, //ptjj
                                                        {-2.4, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.4}, //yjj
                                                        {400, 450, 500, 570, 650, 800, 1100, 1500}, //jetpt0
                                                        {400, 450, 500, 570, 650, 800, 1100, 1500}, //jetpt1
                                                        {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.4}, //jetY0
                                                        {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.4}, //jetY1
                                                        {1,2,3,4,5,6,7,8,9,10,13,16}, //chi
                                                        {-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1}, //|cosTheta*| leading
                                                        {-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1}}; //|cosTheta*| subleading

  std::vector< std::vector <Float_t> > const BND_reco = {{1000, 1200, 1400, 1600, 1800, 2000, 2400, 3000, 5000}, //mJJ
                                                {0, 60, 150, 300, 450, 850, 1300}, //ptJJ
                                                {-2.4, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.4}, //yjj
                                                {400, 450, 500, 570, 650, 800, 1100, 1500}, //jetPt0
                                                {400, 450, 500, 570, 650, 800, 1100, 1500}, //jetPt1
                                                {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.4}, //jetY0
                                                {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.4}, //jetY1
                                                {1,2,3,4,5,6,7,8,9,10,13,16}, //chi
                                                {-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1}, //|cosTheta*| leading
                                                {-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1}}; //|cosTheta*| subleading

  int NBINS[BND_reco.size()];
  for (int i = 0; i<BND_reco.size(); i++) NBINS[i] = BND_reco[i].size()-1;

  int NBINSParton[BND_reco.size()];
  for(int i =0; i<BND_reco.size(); i++) NBINSParton[i] = BND_gen[i].size()-1;

  TString varReco[NVAR]   = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1", "chi", "cosTheta_0", "cosTheta_1"};
  TString varParton[NVAR] = {"mTTbarParton", "ptTTbarParton", "yTTbarParton","partonPt0", "partonPt1", "partonY0", "partonY1","chiParton", "cosThetaParton_0", "cosThetaParton_1"};
  TString varParticle[NVAR] = {"mJJGen", "ptJJGen", "yJJGen","genjetPt0", "genjetPt1", "genjetY0", "genjetY1", "chiParticle", "cosThetaParticle_0", "cosThetaParticle_1"};

  TH1F *hParton[NVAR], *hParticle[NVAR], *hReco[NVAR];
  TH1F *hRecoParton[NVAR], *hPartonReco[NVAR];
  TH1F *hRecoParticle[NVAR], *hParticleReco[NVAR];
  	//declare the histograms
  	for(int ivar =0; ivar<NVAR; ivar++)
  	{
      int sizeBins = NBINS[ivar];
        float tempBND[NBINS[ivar]+1];
        std::copy(BND_reco[ivar].begin(), BND_reco[ivar].end(), tempBND);

        int sizeBinsPartonParticle = NBINSParton[ivar];
        float tempBNDPartonParticle[NBINSParton[ivar]+1];
        std::copy(BND_gen[ivar].begin(), BND_gen[ivar].end(), tempBNDPartonParticle);
  		 //denominators for parton efficiency (hParton vs parton), particle eff (hParticle vs particle) and acceptance (same for both hReco vs reco)
       //denominators for parton efficiency (hParton vs parton), particle eff (hParticle vs particle) and acceptance (same for both hReco vs reco)
       //this is events pass parton cuts vs parton quantity --> BND_gen
       hParton[ivar] = new TH1F(TString::Format("hParton_%s", varParton[ivar].Data()), TString::Format("hParton_%s", varParton[ivar].Data()), sizeBinsPartonParticle, tempBNDPartonParticle);
      //this is events pass reco cuts used as denominator vs reco quantity --> BND_reco
       hReco[ivar] = new TH1F(TString::Format("hReco_%s", varReco[ivar].Data()), TString::Format("hReco_%s", varReco[ivar].Data()), sizeBins, tempBND);
      //this is events pass particle cuts only used as denominator for particle eff vs particle -->BND
      hParticle[ivar] = new TH1F(TString::Format("hParticle_%s", varParticle[ivar].Data()), TString::Format("hParticle_%s", varParticle[ivar].Data()), sizeBinsPartonParticle, tempBNDPartonParticle);

         //numerator for parton efficiency (hPartonReco vs parton) and acceptance (hRecoParton vs reco)
         //this is RecoParton vs reco --> BND_reco
         hRecoParton[ivar] = new TH1F(TString::Format("hRecoParton_%s", varReco[ivar].Data()), TString::Format("hRecoParton_%s", varReco[ivar].Data()), sizeBins, tempBND);
         hPartonReco[ivar] = new TH1F(TString::Format("hPartonReco_%s", varParton[ivar].Data()), TString::Format("hPartonReco_%s", varParton[ivar].Data()), sizeBinsPartonParticle, tempBNDPartonParticle);

         //numerator for particle efficiency (hParticleReco vs particle) and acceptance (hRecoParticle vs reco)
         hRecoParticle[ivar] = new TH1F(TString::Format("hRecoParticle_%s", varReco[ivar].Data()), TString::Format("hRecoParticle_%s", varReco[ivar].Data()), sizeBins, tempBND);
         hParticleReco[ivar] = new TH1F(TString::Format("hParticleReco_%s", varParticle[ivar].Data()), TString::Format("hParticleReco_%s", varParticle[ivar].Data()), sizeBinsPartonParticle, tempBNDPartonParticle);
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
    vector<float> *jetMassSoftDrop(0), *partonEta(0), *partonY(0);

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
  trIN->SetBranchAddress("partonY"		  ,&partonY);
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
	float weights = XSEC/norm;

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
  std::vector<float> *partonY_ = new std::vector<float>(0);
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
 
	   xPartonAll.clear();
	   xRecoAll.clear();
	   xParticleAll.clear();

	if(nJets > 1)
    {

   	   	int leadingPt =0;
      	int subleadingPt = 1;

  		if((*jetPt)[0] < (*jetPt)[1])
  		{
   		    subleadingPt =0;
   		    leadingPt = 1;
   		}

      int leadingPt_parton =0;
      int subleadingPt_parton = 1;
      if((*partonPt)[0] < (*partonPt)[1])
  		{
   		  subleadingPt_parton =0;
   		  leadingPt_parton = 1;
   		}

      int leadingPt_particle =0;
      int subleadingPt_particle = 1;
      if((*genjetPt)[0] < (*genjetPt)[1])
  		{

   		  subleadingPt_particle =0;
   		  leadingPt_particle = 1;
   		}

      //reco for angular
    TLorentzVector p4T[2], p4T_ZMF[2], p4TTbar;
    p4T[leadingPt].SetPtEtaPhiM((*jetPt)[leadingPt], (*jetEta)[leadingPt], (*jetPhi)[leadingPt], (*jetMassSoftDrop)[leadingPt]);
   	p4T[subleadingPt].SetPtEtaPhiM((*jetPt)[subleadingPt], (*jetEta)[subleadingPt], (*jetPhi)[subleadingPt], (*jetMassSoftDrop)[subleadingPt]);

  	TVector3 ttbarBoostVector = getBoostVector(p4T[leadingPt], p4T[subleadingPt], p4TTbar);

    p4T_ZMF[0].SetPtEtaPhiM(p4T[leadingPt].Pt(), p4T[leadingPt].Eta(), p4T[leadingPt].Phi(), p4T[leadingPt].M());
    p4T_ZMF[1].SetPtEtaPhiM(p4T[subleadingPt].Pt(), p4T[subleadingPt].Eta(), p4T[subleadingPt].Phi(), p4T[subleadingPt].M());
    p4T_ZMF[0].Boost(ttbarBoostVector);
    p4T_ZMF[1].Boost(ttbarBoostVector);
	    float yStarExp  = TMath::Exp(fabs(p4T_ZMF[0].Rapidity() - p4T_ZMF[1].Rapidity())); //this is chi = e^(|y*|) , y* = 1/2(y1-y0)


    xRecoAll.push_back(mJJ);
		xRecoAll.push_back(ptJJ);
		xRecoAll.push_back(yJJ);
		xRecoAll.push_back((*jetPt)[leadingPt]);
		xRecoAll.push_back((*jetPt)[subleadingPt]);
		xRecoAll.push_back(fabs((*jetY)[leadingPt]));
		xRecoAll.push_back(fabs((*jetY)[subleadingPt]));
    xRecoAll.push_back(yStarExp); //this is chi
    xRecoAll.push_back(TMath::Cos(p4T_ZMF[0].Theta())); //this is |cos(theta*)| leading
    xRecoAll.push_back(TMath::Cos(p4T_ZMF[1].Theta())); //this is |cos(theta*)| subleading

	  //now parton
		TLorentzVector p4TParton[2], p4T_ZMFParton[2], p4TTbarParton;
    p4TParton[leadingPt_parton].SetPtEtaPhiM((*partonPt)[leadingPt_parton], (*partonEta)[leadingPt_parton], (*partonPhi)[leadingPt_parton], (*partonMass)[leadingPt_parton]);
   	p4TParton[subleadingPt_parton].SetPtEtaPhiM((*partonPt)[subleadingPt_parton], (*partonEta)[subleadingPt_parton], (*partonPhi)[subleadingPt_parton], (*partonMass)[subleadingPt_parton]);

  	TVector3 ttbarBoostVectorParton = getBoostVector(p4TParton[leadingPt], p4TParton[subleadingPt], p4TTbarParton);

    p4T_ZMFParton[0].SetPtEtaPhiM(p4TParton[leadingPt_parton].Pt(), p4TParton[leadingPt_parton].Eta(), p4TParton[leadingPt_parton].Phi(), p4TParton[leadingPt_parton].M());
    p4T_ZMFParton[1].SetPtEtaPhiM(p4TParton[subleadingPt_parton].Pt(), p4TParton[subleadingPt_parton].Eta(), p4TParton[subleadingPt_parton].Phi(), p4TParton[subleadingPt_parton].M());
    p4T_ZMFParton[0].Boost(ttbarBoostVectorParton);
    p4T_ZMFParton[1].Boost(ttbarBoostVectorParton);

	  float yStarExpParton  = TMath::Exp(fabs(p4T_ZMFParton[0].Rapidity() - p4T_ZMFParton[1].Rapidity())); //this is chi = e^(|y*|) , y* = 1/2(y1-y0)

		xPartonAll.push_back(mTTbarParton);
		xPartonAll.push_back(ptTTbarParton);
		xPartonAll.push_back(yTTbarParton);
		xPartonAll.push_back((*partonPt)[leadingPt_parton]);
		xPartonAll.push_back((*partonPt)[subleadingPt_parton]);
		xPartonAll.push_back(fabs((*partonY)[leadingPt_parton]));
		xPartonAll.push_back(fabs((*partonY)[subleadingPt_parton]));
    xPartonAll.push_back(yStarExpParton);
		xPartonAll.push_back(TMath::Cos(p4T_ZMFParton[0].Theta())); //this is |cos(theta*)| leading
		xPartonAll.push_back(TMath::Cos(p4T_ZMFParton[1].Theta())); //this is |cos(theta*)| subleading

    //now particle
		TLorentzVector p4TParticle[2], p4T_ZMFParticle[2], p4TTbarParticle;
    p4TParticle[leadingPt_particle].SetPtEtaPhiM((*genjetPt)[leadingPt_particle], (*genjetEta)[leadingPt_particle], (*genjetPhi)[leadingPt_particle], (*genjetMassSoftDrop)[leadingPt_particle]);
   	p4TParticle[subleadingPt_particle].SetPtEtaPhiM((*genjetPt)[subleadingPt_particle], (*genjetEta)[subleadingPt_particle], (*genjetPhi)[subleadingPt_particle], (*genjetMassSoftDrop)[subleadingPt_particle]);

  	TVector3 ttbarBoostVectorParticle = getBoostVector(p4TParticle[leadingPt_particle], p4TParticle[subleadingPt_particle], p4TTbarParticle);

    p4T_ZMFParticle[0].SetPtEtaPhiM(p4TParticle[leadingPt_particle].Pt(), p4TParticle[leadingPt_particle].Eta(), p4TParticle[leadingPt_particle].Phi(), p4TParticle[leadingPt_particle].M());
    p4T_ZMFParticle[1].SetPtEtaPhiM(p4TParticle[subleadingPt_particle].Pt(), p4TParticle[subleadingPt_particle].Eta(), p4TParticle[subleadingPt_particle].Phi(), p4TParticle[subleadingPt_particle].M());
    p4T_ZMFParticle[0].Boost(ttbarBoostVectorParticle);
    p4T_ZMFParticle[1].Boost(ttbarBoostVectorParticle);

	  float yStarExpParticle  = TMath::Exp(fabs(p4T_ZMFParticle[0].Rapidity() - p4T_ZMFParticle[1].Rapidity())); //this is chi = e^(|y*|) , y* = 1/2(y1-y0)

		xParticleAll.push_back(mJJGen);
		xParticleAll.push_back(ptJJGen);
		xParticleAll.push_back(yJJGen);
		xParticleAll.push_back((*genjetPt)[leadingPt_particle]);
		xParticleAll.push_back((*genjetPt)[subleadingPt_particle]);
		xParticleAll.push_back(fabs((*genjetY)[leadingPt_particle]));
		xParticleAll.push_back(fabs((*genjetY)[subleadingPt_particle]));
    xParticleAll.push_back(yStarExpParticle);
		xParticleAll.push_back(TMath::Cos(p4T_ZMFParticle[0].Theta())); //this is |cos(theta*)| leading
		xParticleAll.push_back(TMath::Cos(p4T_ZMFParticle[1].Theta())); //this is |cos(theta*)| subleading

	  //---------------------------end of MATCHING---------------------------------------------------------
    bool recoCuts, partonCuts, particleCuts;
	  bool massCut = (*jetMassSoftDrop)[0] > 120 && (*jetMassSoftDrop)[0] < 220 && (*jetMassSoftDrop)[1] > 120 && (*jetMassSoftDrop)[1] < 220;
	  bool tTaggerCut = (*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut;
	  recoCuts = nJets > 1 && fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1]) <2.4 && (*jetPt)[0] > 400 && (*jetPt)[1] > 400 && mJJ > 1000 && massCut && nLeptons==0 && (*bit)[triggerFloat];
	  partonCuts = fabs((*partonEta)[0]) < 2.4 && fabs((*partonEta)[1]) <2.4 && (*partonPt)[0] > 400 && (*partonPt)[1] > 400 && mTTbarParton > 1000;
	  particleCuts = fabs((*genjetEta)[0]) < 2.4 && fabs((*genjetEta)[1]) <2.4 && (*genjetPt)[0] > 400 && (*genjetPt)[1] > 400 && mJJGen > 1000 && nJetsGen >1 &&
	  				 (*genjetMassSoftDrop)[0] > 120 && (*genjetMassSoftDrop)[0] < 220 && (*genjetMassSoftDrop)[1] > 120 && (*genjetMassSoftDrop)[1] < 220;
	  bool deepCSV = (((*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0])> deepCSVFloat) &&
					 (((*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1])> deepCSVFloat);

      bool btagCut;
	  btagCut = deepCSV;

	  //qcout<<"----------"<<endl;

		  //fill the denominators
		  //1. denominator passing only reco cuts for topTagger (same for parton and particle)
		  if(recoCuts && btagCut && tTaggerCut)
		  {
		  	for(int ivar = 0; ivar < NVAR; ivar++)
	  		{
				hReco[ivar]->Fill(xRecoAll[ivar], genEvtWeight*bTagEvntWeight);
			}
		  }
		  //2. fill the histograms pass reco and parton cuts numerators for efficiencies and acceptance
		  //fill also response matrix
	      if(partonCuts && recoCuts && tTaggerCut && btagCut)
		  {
			  	for(int ivar = 0; ivar < NVAR; ivar++)
	  			{
				   hPartonReco[ivar]->Fill(xPartonAll[ivar], genEvtWeight*bTagEvntWeight);
				   hRecoParton[ivar]->Fill(xRecoAll[ivar], genEvtWeight*bTagEvntWeight);
  				}//---- end of the ivar loop

	      }//----- end of selection cuts parton and reco

	      //3. Events pass reco and particle
	      if(particleCuts && recoCuts && tTaggerCut && btagCut)
	      {
	      	for(int ivar = 0; ivar < NVAR; ivar++)
	  		  {
	      		hParticleReco[ivar]->Fill(xParticleAll[ivar], genEvtWeight*bTagEvntWeight);
	      		hRecoParticle[ivar]->Fill(xRecoAll[ivar], genEvtWeight*bTagEvntWeight);
	      	}
	      }
	      if(particleCuts)
	      {
	      	for(int ivar = 0; ivar < NVAR; ivar++)
	  		{
	      		hParticle[ivar]->Fill(xParticleAll[ivar], genEvtWeight*bTagEvntWeight);
	      	}
	      }

	 }//----- end of nJets > 1
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

  /*
  float norm = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
	float weights = XSEC/norm;
  */
  for(int iev = 0; iev < NNCnt; iev++)
  {
	  trCnt->GetEntry(iev);

	  std::vector<float> xPartonAllCnt(0);
	  xPartonAllCnt.clear();
	  bool partonCuts = fabs(partonEtaCnt[0]) < 2.4 && fabs(partonEtaCnt[1]) <2.4 && partonPtCnt[0] > 400 && partonPtCnt[1] > 400 && mTTbarPartonCnt > mJJCut;
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
    xPartonAllCnt.push_back(mTTbarPartonCnt);
	  xPartonAllCnt.push_back(ptTTbarPartonCnt);
	  xPartonAllCnt.push_back(yTTbarPartonCnt);
    xPartonAllCnt.push_back(partonPtCnt[leadingPt]);
    xPartonAllCnt.push_back(partonPtCnt[subleadingPt]);
    xPartonAllCnt.push_back(fabs(partonYCnt[leadingPt]));
    xPartonAllCnt.push_back(fabs(partonYCnt[subleadingPt]));
	  xPartonAllCnt.push_back(yStarExpPartonCnt);
	  xPartonAllCnt.push_back(TMath::Cos(p4T_ZMFPartonCnt[0].Theta())); //this is |cos(theta*)| leading
	  xPartonAllCnt.push_back(TMath::Cos(p4T_ZMFPartonCnt[1].Theta())); //this is |cos(theta*)| subleading

	  for(int ivar = 0; ivar < NVAR; ivar++)
	  {
	  	if(partonCuts)
			hParton[ivar]->Fill(xPartonAllCnt[ivar], genEvtWeightCnt);
	  }
  }


  for(int ivar =0; ivar<NVAR; ivar++)
  {
    hReco[ivar]->Scale(weights*LUMI);
    hParticle[ivar]->Scale(weights*LUMI);

    hRecoParticle[ivar]->Scale(weights*LUMI);
    hParticleReco[ivar]->Scale(weights*LUMI);

    hRecoParton[ivar]->Scale(weights*LUMI);
    hPartonReco[ivar]->Scale(weights*LUMI);

    hParton[ivar]->Scale(weights*LUMI);
  }
  

  TFile *outFile = TFile::Open(TString::Format("%s/bTagVariation/Histograms_%s_%s.root", year.Data(), ttbar_process.Data(),variation.Data()), "RECREATE");
  //outFile->cd();
  //write them to file
  for(int ivar = 0; ivar<NVAR; ivar++)
  {
  	hReco[ivar]->Write();
    hParticle[ivar]->Write();

    hRecoParticle[ivar]->Write();
    hParticleReco[ivar]->Write();

    hRecoParton[ivar]->Write();
    hPartonReco[ivar]->Write();

    hParton[ivar]->Write();
  }//end of ivar

  outFile->Close();
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
