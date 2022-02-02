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

void FillHistograms_Reduced_JES(TString file_name, TString ttbar_process, TString jes_variation,  TString year = "2016", float mJJCut = 1000)
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
  const int NVAR = 12;
  std::vector< std::vector <Float_t> > const BND = {{1000, 1200, 1400, 1600, 1800, 2000, 2400, 3000, 5000},
                       {0, 60, 150, 300, 450, 850, 1300}, //ptJJ
                       {-2.4, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.4}, //yjj
                       {450, 500, 570, 650, 800, 1100, 1500}, //jetPt0
                       {400, 450, 500, 570, 650, 800, 1100, 1500}, //jetPt1
                       {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.4}, //jetY0
                       {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.4}, //jetY1
                       {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 16}, //chi
                       {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1}, //|cosTheta*| leading
                       {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1}}; //|cosTheta*| subleading


  int NBINS[BND.size()];
  for (int i = 0; i<BND.size()+2; i++)
  {
    if(i<10) NBINS[i] = BND[i].size()-1;
    else NBINS[i] = 100;
  }
  TString varReco[NVAR]   = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1","chi","cosThjetEta0", "cosThjetEta1","mTop_Leading", "mTop_Subleading"};

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

    TTree *trIN = (TTree*)file->Get(TString::Format("%s/events", jes_variation.Data()));

	//------- input tree --------------
    trIN->SetBranchAddress("nJets"          ,&nJets);
    trIN->SetBranchAddress("nLeptons"       ,&nLeptons);
    trIN->SetBranchAddress("jetPt"          ,&jetPt);
    trIN->SetBranchAddress("jetEta"         ,&jetEta);
    trIN->SetBranchAddress("jetY"           ,&jetY);
    trIN->SetBranchAddress("jetPhi"         ,&jetPhi);
    trIN->SetBranchAddress("triggerBit"     ,&bit);
    trIN->SetBranchAddress("genEvtWeight"   ,&genEvtWeight);
    trIN->SetBranchAddress("bTagEvntWeight", &bTagEvntWeight);
    trIN->SetBranchAddress("mJJ"   			,&mJJ);
    trIN->SetBranchAddress("yJJ"   			,&yJJ);
    trIN->SetBranchAddress("ptJJ"   		,&ptJJ);
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

  cout<<ttbar_process<<endl;
  /*
  cout<<eospath_nom[year]<<ttNominalFiles[year][ttbar_process]<<endl;
  TFile *inf_GenEventWeight = TFile::Open(eospath_nom[year]+ttNominalFiles[year][ttbar_process]);

  float norm = ((TH1F*)inf_GenEventWeight->Get("eventCounter/GenEventWeight"))->GetSumOfWeights(); */

  float norm = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
	float weight = XSEC/norm;
	weights = weight;
  cout<<"norm: "<<norm<<endl;

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
	   xRecoAll.clear();
    if(nJets>1)
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

      xRecoAll.push_back(mJJ);
      xRecoAll.push_back(ptJJ);
      xRecoAll.push_back(yJJ);
      xRecoAll.push_back((*jetPt)[leadingPt]);
      xRecoAll.push_back((*jetPt)[subleadingPt]);
      xRecoAll.push_back(fabs((*jetY)[leadingPt]));
      xRecoAll.push_back(fabs((*jetY)[subleadingPt]));
      xRecoAll.push_back(yStarExp); //this is chi
     	xRecoAll.push_back(fabs(TMath::Cos(p4T_ZMF[0].Theta()))); //this is |cos(theta*)| leading
     	xRecoAll.push_back(fabs(TMath::Cos(p4T_ZMF[1].Theta()))); //this is |cos(theta*)| subleading
      xRecoAll.push_back((*jetMassSoftDrop)[leadingPt]);
      xRecoAll.push_back((*jetMassSoftDrop)[subleadingPt]);


	  //---------------------------end of MATCHING---------------------------------------------------------
    float dCSVScoreSub0[2], dCSVScoreSub1[2];
    dCSVScoreSub0[0] = (*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0];
    dCSVScoreSub0[1] = (*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1];
    dCSVScoreSub1[0] = (*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0];
    dCSVScoreSub1[1] = (*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1];

	  bool recoCuts;
	  bool massCut = (*jetMassSoftDrop)[0] > 120 && (*jetMassSoftDrop)[0] < 220 && (*jetMassSoftDrop)[1] > 120 && (*jetMassSoftDrop)[1] < 220;
	  bool tTaggerCut = (*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut;
	  recoCuts = nJets > 1 && fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1]) <2.4 && (*jetPt)[0] > 450 && (*jetPt)[1] > 400 && mJJ> mJJCut && massCut && nLeptons==0 && (*bit)[triggerFloat];
	  bool deepCSV = (((*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0])> deepCSVFloat) &&
					 (((*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1])> deepCSVFloat);

    bool revertBtag = (dCSVScoreSub0[0] < deepCSVFloat &&  dCSVScoreSub1[0] < deepCSVFloat) && (dCSVScoreSub0[1] < deepCSVFloat && dCSVScoreSub1[1] < deepCSVFloat);

    bool btagCut;
	  btagCut = deepCSV;
    //Signal Region 2btags
		if(recoCuts && btagCut && tTaggerCut)
		{
		  for(int ivar = 0; ivar < NVAR; ivar++)
	  	{
        double weights_temp = genEvtWeight * bTagEvntWeight;
		    hReco[ivar]->Fill(xRecoAll[ivar], genEvtWeight*bTagEvntWeight);
        //cout<<genEvtWeight*bTagEvntWeight<<endl;
		  }
	  }
    //Control Region 0btag
    if(recoCuts && revertBtag && tTaggerCut)
	  {
	  	for(int ivar = 0; ivar < NVAR; ivar++)
  		{
		    hRecoCR[ivar]->Fill(xRecoAll[ivar], genEvtWeight*bTagEvntWeight);
		  }
	  }
  }//---end nJets>1
  }//---end the event for


    for(int ivar =0; ivar<NVAR; ivar++)
    {
      //cout<<ivar<<endl;
      //cout<<"Integral before scaling: " <<hReco[ivar]->Integral()<<endl;
      hReco[ivar]->Scale(weights*LUMI);
      hRecoCR[ivar]->Scale(weights*LUMI_CR);
    }//end of loop on all vars


    TFile *outFile;
    //outFile = TFile::Open(TString::Format("%s/JES/HistoReduced_%d_%s_%s.root", year.Data(),(int)mJJCut, ttbar_process.Data(), jes_variation.Data()), "RECREATE");
    outFile = TFile::Open(TString::Format("HistoReduced_%d_%s_%s.root", (int)mJJCut, ttbar_process.Data(), jes_variation.Data()), "RECREATE");
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
