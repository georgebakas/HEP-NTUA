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

void TagAndProbe_extra_data(TString file_name, TString year = "2016", float mJJCut = 1000)
{
  globalYear = year;
  initFilesMapping();
  cout<<"YEAR: "<<year<<endl;
  cout<<"file_name: "<<file_name<<endl;
  float triggerFloat;
  float triggerFloatCR;
  if(year.Contains("2016"))
  {
    triggerFloatCR = 4;
    triggerFloat = 2;
  } 
  else
  {
    triggerFloat = 5;
    triggerFloatCR = 5;
  } 

  float deepCSVFloat = floatConstants[TString::Format("btagWP%s",year.Data())];
  float selMvaCut = topTaggerConstants[TString::Format("topTagger%s",year.Data())];
  float LUMI = luminosity[TString::Format("luminosity%s", year.Data())];
  float LUMI_CR = luminosityCR[TString::Format("luminosity%s", year.Data())];
  float tightTopTaggerCut = 0.8;
  const int NVAR = 3;
  int NBINS = 20;
  TString varReco[NVAR]   = {"jetPt", "jetEta", "topTagger"};

  // I need the following histograms 
  // 1. Eta, JetPt and Top Tagging for MC filled with probe jets from events passing tight + probe 
  // 2. Eta, JetPt and Top Tagging double filled with events passing SR + SR --> double fill 
  // 3. Eta, JetPt and Top Tagging filled random jet which events passing SR + SR --> choose randomly 
  // 2D eta, topTagger --> 1 
  // 2D eta, topTagger --> 2
  // 2D pT, topTagger --> 3
  TH1F *h_Probe[NVAR], *h_SR_random[NVAR], *h_SR_double[NVAR];
  TH2F *hPtTopTagger_probe, *hPtTopTagger_random, *hPtTopTagger_double;
  TH2F *hEtaTopTagger_probe, *hEtaTopTagger_random, *hEtaTopTagger_double;
  int maxJetMass = 220;
  int minJetMass = 120;
  	//declare the histograms
  	for(int ivar =0; ivar<NVAR; ivar++)
  	{
      int sizeBins = NBINS;
      if (ivar==0)
      {
        h_Probe[ivar] = new TH1F(TString::Format("hProbe%s", varReco[ivar].Data()), TString::Format("hProbe%s",varReco[ivar].Data()), sizeBins, 400, 1500);
        h_SR_random[ivar] = new TH1F(TString::Format("h_SR_random%s", varReco[ivar].Data()), TString::Format("h_SR_random%s",varReco[ivar].Data()), sizeBins, 400, 1500);
        h_SR_double[ivar] = new TH1F(TString::Format("h_SR_double%s", varReco[ivar].Data()), TString::Format("h_SR_double%s",varReco[ivar].Data()), sizeBins, 400, 1500);
      }
      else if (ivar==1)
      {
        h_Probe[ivar] = new TH1F(TString::Format("hProbe%s", varReco[ivar].Data()), TString::Format("hProbe%s",varReco[ivar].Data()), sizeBins, 0, 2.5);
        h_SR_random[ivar] = new TH1F(TString::Format("h_SR_random%s", varReco[ivar].Data()), TString::Format("h_SR_random%s",varReco[ivar].Data()), sizeBins, 0, 2.5);
        h_SR_double[ivar] = new TH1F(TString::Format("h_SR_double%s", varReco[ivar].Data()), TString::Format("h_SR_double%s",varReco[ivar].Data()), sizeBins, 0, 2.5);
      }
      else if (ivar==2)
      {
        h_Probe[ivar] = new TH1F(TString::Format("hProbe%s", varReco[ivar].Data()), TString::Format("hProbe%s",varReco[ivar].Data()), sizeBins, -1, 1);
        h_SR_random[ivar] = new TH1F(TString::Format("h_SR_random%s", varReco[ivar].Data()), TString::Format("h_SR_random%s",varReco[ivar].Data()), sizeBins, -1, 1);
        h_SR_double[ivar] = new TH1F(TString::Format("h_SR_double%s", varReco[ivar].Data()), TString::Format("h_SR_double%s",varReco[ivar].Data()), sizeBins, -1, 1);
      }
      h_Probe[ivar]->Sumw2();
      h_SR_random[ivar]->Sumw2();
      h_SR_double[ivar]->Sumw2();
    }
    // x-axis: pT, y-axis: topTagger
    hPtTopTagger_probe = new TH2F("hPtTopTagger_probe", "hPtTopTagger_probe", NBINS, 400, 1500, NBINS, -1, 1);
    hPtTopTagger_random = new TH2F("hPtTopTagger_random", "hPtTopTagger_random", NBINS, 400, 1500, NBINS, -1, 1);
    hPtTopTagger_double = new TH2F("hPtTopTagger_double", "hPtTopTagger_double", NBINS, 400, 1500, NBINS, -1, 1);

    // x-axis: eta, y-axis: topTagger
    hEtaTopTagger_probe = new TH2F("hEtaTopTagger_probe", "hEtaTopTagger_probe", NBINS, 0, 2.5, NBINS, -1, 1);
    hEtaTopTagger_random = new TH2F("hEtaTopTagger_random", "hEtaTopTagger_random", NBINS, 0, 2.5, NBINS, -1, 1);
    hEtaTopTagger_double = new TH2F("hEtaTopTagger_double", "hEtaTopTagger_double", NBINS, 0, 2.5, NBINS, -1, 1);


    int nJets,nLeptons, category(0);
    vector<bool>  *bit(0),*matchedJet(0);
    //reco vars:
    std::vector<float> *jetPt(0), *jetY(0), *jetEta(0), *jetPhi(0), *jetTtag(0), *deepAK8(0);;
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
    //trIN->SetBranchAddress("jetTau3"        ,&tau3);
    //trIN->SetBranchAddress("jetTau2"        ,&tau2);
    //trIN->SetBranchAddress("jetTau1"        ,&tau1);
    trIN->SetBranchAddress("triggerBit"     ,&bit);
    trIN->SetBranchAddress("bTagEvntWeight", &bTagEvntWeight);
    //trIN->SetBranchAddress("jetMassSub0"    ,&jetMassSub0);
    //trIN->SetBranchAddress("jetMassSub1"    ,&jetMassSub1);
    trIN->SetBranchAddress("mJJ"   			,&mJJ);
    trIN->SetBranchAddress("yJJ"   			,&yJJ);
    trIN->SetBranchAddress("ptJJ"   		,&ptJJ);
	//trIN->SetBranchAddress("jetBtagSub0"	,&jetBtagSub0);
    //trIN->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1);
    trIN->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
    //trIN->SetBranchAddress("deepAK8Tagger"  ,&deepAK8);
  	trIN->SetBranchAddress("jetTtagCategory",&jetTtag);
  	//deepCSV
    trIN->SetBranchAddress("jetBtagSub0DCSVbb" ,&jetBtagSub0DCSVbb);
    trIN->SetBranchAddress("jetBtagSub1DCSVbb" ,&jetBtagSub1DCSVbb);
    trIN->SetBranchAddress("jetBtagSub0DCSVbbb",&jetBtagSub0DCSVbbb);
    trIN->SetBranchAddress("jetBtagSub1DCSVbbb",&jetBtagSub1DCSVbbb);
    //general
    int evtNo;
	trIN->SetBranchAddress("evtNo", &evtNo);

    float tTaggerTight = 0;
    float tTaggerOther = 0;

    int decade(0);
    int NN = trIN->GetEntries();
    //NN = 100000;
    std::cout<<"Entries: "<<NN<<std::endl;
	  std::vector<float> xRecoAll(0);
    std::vector<float> xRecoAll_random(0);
    std::vector<float> xRecoAllSR_double(0);

    for(int iev=0;iev<NN;iev++)
    {
        tTaggerTight = 0;
        tTaggerOther = 0;
        double progress = 10.0*iev/(1.0*NN);
        int k = TMath::FloorNint(progress);
        if (k > decade)
            cout<<10*k<<" %"<<endl;
        decade = k;
        trIN->GetEntry(iev);
        xRecoAll.clear();
        xRecoAll_random.clear();
        xRecoAllSR_double.clear();

        if(nJets>1)
        {
            int leadingPt = 0;
            int subleadingPt = 1;

            TRandom2 *randJet = new TRandom2();
            int tightJet=0;
            int otherJet=0;
            if (randJet->Rndm() > 0.5)
            {
                tightJet = 1;
                otherJet = 0;
            }
            else
            {
                tightJet = 0;
                otherJet = 1;
            }

            tTaggerTight = (*jetTtag)[tightJet];
            tTaggerOther = (*jetTtag)[otherJet];
            //for probe jets
            xRecoAll.push_back((*jetPt)[otherJet]);
            xRecoAll.push_back(fabs((*jetEta)[otherJet]));
            xRecoAll.push_back((*jetTtag)[otherJet]);

            //for random jet in SR
            TRandom2 *ran = new TRandom2();
            int rndm_jet = ran->Integer(2);
            xRecoAll_random.push_back((*jetPt)[rndm_jet]);
            xRecoAll_random.push_back(fabs((*jetEta)[rndm_jet]));
            xRecoAll_random.push_back((*jetTtag)[rndm_jet]);

            // for double fill 
            xRecoAllSR_double.push_back((*jetPt)[0]);
            xRecoAllSR_double.push_back(fabs((*jetEta)[0]));
            xRecoAllSR_double.push_back((*jetTtag)[0]);
            xRecoAllSR_double.push_back((*jetPt)[1]);
            xRecoAllSR_double.push_back(fabs((*jetEta)[1]));
            xRecoAllSR_double.push_back((*jetTtag)[1]);
            

            float dCSVScoreSub0[2], dCSVScoreSub1[2];
            dCSVScoreSub0[0] = (*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0];
            dCSVScoreSub0[1] = (*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1];
            dCSVScoreSub1[0] = (*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0];
            dCSVScoreSub1[1] = (*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1];

            bool recoCuts;
            bool massCut = (*jetMassSoftDrop)[0] > 120 && (*jetMassSoftDrop)[0] < 220 && (*jetMassSoftDrop)[1] > 120 && (*jetMassSoftDrop)[1] < 220;
            bool tTaggerCut = (*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut;
            recoCuts = nJets > 1 && fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1]) <2.4 && (*jetPt)[leadingPt] > 450 && (*jetPt)[subleadingPt] > 400 && mJJ> mJJCut && massCut && nLeptons==0;
            bool deepCSV = (((*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0])> deepCSVFloat) &&
                            (((*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1])> deepCSVFloat);

            bool revertBtag = (dCSVScoreSub0[0] < deepCSVFloat &&  dCSVScoreSub1[0] < deepCSVFloat) && (dCSVScoreSub0[1] < deepCSVFloat && dCSVScoreSub1[1] < deepCSVFloat);

            bool btagCut;
            btagCut = deepCSV;
            for(int ivar = 0; ivar < xRecoAll.size(); ivar++)
            {
              //Probe Region 2btags
              if(recoCuts && btagCut && (*bit)[triggerFloat])
              { 
                  // tight and probe region 
                  if(tTaggerTight > tightTopTaggerCut)
                  {
                      h_Probe[ivar]->Fill(xRecoAll[ivar]);
                  }
                  // SR and SR region 
                  if(tTaggerTight > selMvaCut && tTaggerOther > selMvaCut)
                  {
                      h_SR_random[ivar]->Fill(xRecoAll_random[ivar]);
                      //double fill 
                      h_SR_double[ivar]->Fill(xRecoAllSR_double[ivar]);
                      h_SR_double[ivar]->Fill(xRecoAllSR_double[ivar+3]);
                  }
              }
            }

            // 2D plots here 
            //TH2F *hPtTopTagger_probe, *hPtTopTagger_random, *hPtTopTagger_double;
            // TH2F *hEtaTopTagger_probe, *hEtaTopTagger_random, *hEtaTopTagger_double;
            if(recoCuts && btagCut && (*bit)[triggerFloat])
            { 
                // tight and probe region 
                if(tTaggerTight > tightTopTaggerCut)
                {
                    hPtTopTagger_probe->Fill(xRecoAll[0], xRecoAll[2]);
                    hEtaTopTagger_probe->Fill(xRecoAll[1], xRecoAll[2]);
                }
                // SR and SR region 
                if(tTaggerTight > selMvaCut && tTaggerOther > selMvaCut)
                {
                    hPtTopTagger_random->Fill(xRecoAll_random[0], xRecoAll_random[2]);
                    hEtaTopTagger_random->Fill(xRecoAll_random[1], xRecoAll_random[2]);
                    //double fill 
                    hPtTopTagger_double->Fill(xRecoAllSR_double[0], xRecoAllSR_double[2]);
                    hPtTopTagger_double->Fill(xRecoAllSR_double[3], xRecoAllSR_double[5]);
                    hEtaTopTagger_double->Fill(xRecoAllSR_double[1], xRecoAllSR_double[2]);
                    hEtaTopTagger_double->Fill(xRecoAllSR_double[4], xRecoAllSR_double[5]);
                }
            }

            // fill tagNprobe 

        }//---end nJets>1
    }//---end the event for


    TFile *outFile;
    outFile = TFile::Open(TString::Format("%s/TTbarExtraAnalysis_data.root", year.Data()), "RECREATE");

    for(int ivar =0; ivar<NVAR; ivar++)
    {
        h_Probe[ivar]->Write(TString::Format("hProbe%s", varReco[ivar].Data()));
        h_SR_random[ivar]->Write(TString::Format("h_SR_random%s", varReco[ivar].Data()));
        h_SR_double[ivar]->Write(TString::Format("h_SR_double%s", varReco[ivar].Data()));
    }
    hPtTopTagger_probe->Write("hPtTopTagger_probe");
    hPtTopTagger_random->Write("hPtTopTagger_random");
    hPtTopTagger_double->Write("hPtTopTagger_double");
    hEtaTopTagger_probe->Write("hEtaTopTagger_probe");
    hEtaTopTagger_random->Write("hEtaTopTagger_random");
    hEtaTopTagger_double->Write("hEtaTopTagger_double");

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
