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

void TagAndProbe_nominal(TString file_name, TString ttbar_process, TString year = "2016", float mJJCut = 1000)
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
  const int NVAR = 14;
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
  for (int i = 0; i<BND.size()+4; i++)
  {
    if(i<10) NBINS[i] = BND[i].size()-1;
    else NBINS[i] = 100;
  }
  TString varReco[NVAR]   = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1","chi","cosThetEta0", "cosThetEta1",\
                            "mTop_Leading", "mTop_Subleading", "topTagger_leading", "topTagger_Subleading"};

  float weights;
  TH1F *h_Numerator[NVAR], *h_Denominator[NVAR];
  int maxJetMass = 220;
  int minJetMass = 120;
  	//declare the histograms
  	for(int ivar =0; ivar<NVAR; ivar++)
  	{
      int sizeBins = NBINS[ivar];
      if ((ivar==10) || (ivar==11))
      {
        h_Numerator[ivar] = new TH1F(TString::Format("hSRBTightAndSR%s", varReco[ivar].Data()), TString::Format("hSRBTightAndSR%s",varReco[ivar].Data()), sizeBins, minJetMass,maxJetMass);
        h_Denominator[ivar] = new TH1F(TString::Format("hSRBTightAndProbe_%s", varReco[ivar].Data()), TString::Format("hSRBTightAndProbe_%s",varReco[ivar].Data()), sizeBins, minJetMass,maxJetMass);
      }
      else if ((ivar==12)  || (ivar==13))//top tagger [-1,1]
      {
        h_Numerator[ivar] = new TH1F(TString::Format("hSRBTightAndSR%s", varReco[ivar].Data()), TString::Format("hSRBTightAndSR%s",varReco[ivar].Data()), sizeBins, -1, 1);
        h_Denominator[ivar] = new TH1F(TString::Format("hSRBTightAndProbe_%s", varReco[ivar].Data()), TString::Format("hSRBTightAndProbe_%s",varReco[ivar].Data()), sizeBins, -1, 1);
      }
      else if ((ivar==14) || (ivar==15))//deep ak8 [0,1]
      {
        h_Numerator[ivar] = new TH1F(TString::Format("hSRBTightAndSR%s", varReco[ivar].Data()), TString::Format("hSRBTightAndSR%s",varReco[ivar].Data()), sizeBins, 0,1);
        h_Denominator[ivar] = new TH1F(TString::Format("hSRBTightAndProbe_%s", varReco[ivar].Data()), TString::Format("hSRBTightAndProbe_%s",varReco[ivar].Data()), sizeBins, 0,1);
      }
      else
      {
        int sizeBins = NBINS[ivar];
        float tempBND[NBINS[ivar]+1];
        std::copy(BND[ivar].begin(), BND[ivar].end(), tempBND);
        h_Numerator[ivar] = new TH1F(TString::Format("hSRBTightAndSR%s", varReco[ivar].Data()), TString::Format("hSRBTightAndSR%s",varReco[ivar].Data()), sizeBins, tempBND);
        h_Denominator[ivar] = new TH1F(TString::Format("hSRBTightAndProbe_%s", varReco[ivar].Data()), TString::Format("hSRBTightAndProbe_%s",varReco[ivar].Data()), sizeBins, tempBND);
      }
      h_Numerator[ivar]->Sumw2();
      h_Denominator[ivar]->Sumw2();
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
	std::vector<float> xRecoAll(0);

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
            // top tagger leading and subleading 
            xRecoAll.push_back((*jetTtag)[leadingPt]);
            xRecoAll.push_back((*jetTtag)[subleadingPt]);
            // deep AK8 leading and subleading
            //xRecoAll.push_back((*deepAK8)[leadingPt]);
            //xRecoAll.push_back((*deepAK8)[subleadingPt]);


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
                  if(tTaggerTight > tightTopTaggerCut)
                  {
                      double weights_temp = genEvtWeight * bTagEvntWeight;
                      h_Denominator[ivar]->Fill(xRecoAll[ivar], genEvtWeight*bTagEvntWeight);
                  }
                  if(tTaggerTight > tightTopTaggerCut && tTaggerOther > selMvaCut)
                  {
                      double weights_temp = genEvtWeight * bTagEvntWeight;
                      h_Numerator[ivar]->Fill(xRecoAll[ivar], genEvtWeight*bTagEvntWeight);
                  }
              }
            }

            // fill tagNprobe 

        }//---end nJets>1
    }//---end the event for


    for(int ivar =0; ivar<NVAR; ivar++)
    {
        h_Denominator[ivar]->Scale(weights*LUMI);
        h_Numerator[ivar]->Scale(weights*LUMI);
    }


    TFile *outFile;
    outFile = TFile::Open(TString::Format("%s/Nominal/TagAndProbeHisto_TT_%d_%s.root", year.Data(),(int)mJJCut, ttbar_process.Data()), "RECREATE");

    for(int ivar = 0; ivar<NVAR; ivar++)
    {
        h_Numerator[ivar]->Write(TString::Format("hSRBTightAndSR_%s_expYield", varReco[ivar].Data()));
        h_Denominator[ivar]->Write(TString::Format("hSRBTightAndProbe_%s_expYield", varReco[ivar].Data()));
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
