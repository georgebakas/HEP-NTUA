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

void FillLeptonsTight(TString file_name, TString ttbar_process, TString year = "2018")
{
  globalYear = year;
  initFilesMapping();
  cout<<"YEAR: "<<year<<endl;
  cout<<"file_name: "<<file_name<<endl;
  cout<<"ttbar_process: "<<ttbar_process<<endl;
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
  float XSEC = XSECAll[year.Data()][ttbar_process.Data()];
  float tightTopTaggerCut = 0.8;
  cout<<"XSEC: "<<XSEC<<endl;
  const int NVAR = 1;
  std::vector< std::vector <Float_t> > const BND = {{450, 500, 570, 650, 800, 1100, 1500}};//jetPt0


  int NBINS[BND.size()];
  for (int i = 0; i<BND.size()+4; i++)
  {
    if(i<10) NBINS[i] = BND[i].size()-1;
    else NBINS[i] = 100;
  }
  TString varReco[NVAR]   = {"jetPt0"};

  float weights;
  TH1F *hLeptons[NVAR], *hHadronic[NVAR];
  int maxJetMass = 220;
  int minJetMass = 120;
  	//declare the histograms
  	for(int ivar =0; ivar<NVAR; ivar++)
  	{
      int sizeBins = NBINS[ivar];
      float tempBND[NBINS[ivar]+1];
      std::copy(BND[ivar].begin(), BND[ivar].end(), tempBND);
      hLeptons[ivar] = new TH1F(TString::Format("hLeptons%s", varReco[ivar].Data()), TString::Format("hLeptons%s",varReco[ivar].Data()), sizeBins, tempBND);
      hLeptons[ivar]->Sumw2();

      hHadronic[ivar] = new TH1F(TString::Format("hHadronic%s", varReco[ivar].Data()), TString::Format("hHadronic%s",varReco[ivar].Data()), sizeBins, tempBND);
      hHadronic[ivar]->Sumw2();
    }
    int nJets,nLeptons, category(0);
    vector<bool>  *bit(0),*matchedJet(0);
    //reco vars:
    std::vector<float> *jetPt(0), *jetY(0), *jetEta(0), *jetPhi(0), *jetTtag(0), *deepAK8(0);;
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

    // leptons 
    std::vector<int> *lepId(0);
    std::vector<float> *lepPt(0), *lepEta(0), *lepPhi(0), *lepE(0);


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
    trIN->SetBranchAddress("genEvtWeight"   ,&genEvtWeight);
    trIN->SetBranchAddress("bTagEvntWeight", &bTagEvntWeight);
    //trIN->SetBranchAddress("jetMassSub0"    ,&jetMassSub0);
    //trIN->SetBranchAddress("jetMassSub1"    ,&jetMassSub1);

    //lepton vars 
    trIN->SetBranchAddress("lepId" ,&lepId);
    trIN->SetBranchAddress("lepPt" ,&lepPt);
    trIN->SetBranchAddress("lepEta",&lepEta);
    trIN->SetBranchAddress("lepPhi", &lepPhi);
    trIN->SetBranchAddress("lepE", &lepE);


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

    float tTaggerTight = 0;
    float tTaggerOther = 0;

    float norm = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
	  float weight = XSEC/norm;
	  weights = weight;
    cout<<"norm: "<<norm<<endl;

    int decade(0);
    int NN = trIN->GetEntries();
    //NN = 100000;
    std::cout<<"Entries: "<<NN<<std::endl;
	  std::vector<float> xRecoAll(0), xRecoAll_hadronic(0);
    TRandom2 *randJet = new TRandom2();

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
        xRecoAll_hadronic.clear();
        if(nJets>0)
        {
            int tightJet=0;
            if (nJets >1)
            { 
              if (randJet->Rndm() > 0.5)
                  tightJet = 1;
              else
                  tightJet = 0;
            }
            else 
              tightJet = 0;

            tTaggerTight = (*jetTtag)[tightJet];

            // fill only the 1 jet
            xRecoAll.push_back((*jetPt)[tightJet]);

            float dCSVScoreSub0[2], dCSVScoreSub1[2];
            dCSVScoreSub0[0] = (*jetBtagSub0DCSVbb)[tightJet] + (*jetBtagSub0DCSVbbb)[tightJet];
            dCSVScoreSub1[0] = (*jetBtagSub1DCSVbb)[tightJet] + (*jetBtagSub1DCSVbbb)[tightJet];

            // needed for leptons + jets
            bool massCut_leptons = (*jetMassSoftDrop)[tightJet] > 120 && (*jetMassSoftDrop)[tightJet] < 220;
            // for (int l=0; l<lepId->size(); l++)
            // {
            //   cout<<"----------------"<<endl;
            //   cout<<"lepId: "<<(*lepId)[l]<<endl;
            //   cout<<"lepPt: "<<(*lepPt)[l]<<endl;
            // }
            bool recoCuts_leptons = nJets > 0 && (*jetPt)[tightJet] > 450 && nLeptons==1 && fabs((*lepId)[0])==13 && (*lepPt)[0] > 35;
            
            // common requirements 
            bool btagCut = (((*jetBtagSub0DCSVbb)[tightJet] + (*jetBtagSub0DCSVbbb)[tightJet])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[tightJet] + (*jetBtagSub1DCSVbbb)[tightJet])> deepCSVFloat);            

            for(int ivar = 0; ivar < xRecoAll.size(); ivar++)
            {
              //Probe Region 2btags
              if(recoCuts_leptons && btagCut)
              {   
                  // tight 
                  if(tTaggerTight > tightTopTaggerCut)
                  {
                      double weights_temp = genEvtWeight * bTagEvntWeight;
                      hLeptons[ivar]->Fill(xRecoAll[ivar], genEvtWeight*bTagEvntWeight);
                  }
              }
            } 

        }//---end nJets>0


        // now check hadronic channel tight
        if(nJets>1)
        {
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

            xRecoAll_hadronic.push_back((*jetPt)[tightJet]);


            float dCSVScoreSub0[2], dCSVScoreSub1[2];
            dCSVScoreSub0[0] = (*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0];
            dCSVScoreSub0[1] = (*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1];
            dCSVScoreSub1[0] = (*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0];
            dCSVScoreSub1[1] = (*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1];

            bool massCut = (*jetMassSoftDrop)[0] > 120 && (*jetMassSoftDrop)[0] < 220 && (*jetMassSoftDrop)[1] > 120 && (*jetMassSoftDrop)[1] < 220;
            bool tTaggerCut = (*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut;
            bool recoCuts = nJets > 1 && fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1]) <2.4 && (*jetPt)[0] > 450 && (*jetPt)[1] > 400 && mJJ> 1000 && massCut && nLeptons==0;
            bool btagCut = (((*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0])> deepCSVFloat) &&
                            (((*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1])> deepCSVFloat);


            for(int ivar = 0; ivar < xRecoAll_hadronic.size(); ivar++)
            {
              //Probe Region 2btags
              if(recoCuts && btagCut && (*bit)[triggerFloat])
              {
                  if(tTaggerTight > tightTopTaggerCut && tTaggerOther > tightTopTaggerCut)
                  {
                      double weights_temp = genEvtWeight * bTagEvntWeight;
                      hHadronic[ivar]->Fill(xRecoAll_hadronic[ivar], genEvtWeight*bTagEvntWeight);
                  }

              }
            }
        }//---end nJets>1



        

    }//---end the event for


    for(int ivar =0; ivar<NVAR; ivar++)
    {
        hLeptons[ivar]->Scale(weights*LUMI);
        hHadronic[ivar]->Scale(weights*LUMI);
    }


    TFile *outFile;
    outFile = TFile::Open(TString::Format("%s/Nominal/LeptonAnalysis_%s.root", year.Data(), ttbar_process.Data()), "RECREATE");

    for(int ivar = 0; ivar<NVAR; ivar++)
    {
        hLeptons[ivar]->Write(TString::Format("hLeptons%s_expYield", varReco[ivar].Data()));
        hHadronic[ivar]->Write(TString::Format("hHadronic%s_expYield", varReco[ivar].Data()));
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
