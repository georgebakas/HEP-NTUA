#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TMVA/Tools.h"

using std::cin;
using std::cout;
using std::endl;

//float deltaRho

TString formOutFileName(TString str)
{
  TString outStr;
  TString outFileName;
  Ssiz_t from = 0;
  while(str.Tokenize(outStr, from, "/"))
  {
    if(outStr.Contains(".root"))
    {
      outFileName = outStr;
      outFileName.ReplaceAll(".root", "");
      break;
    }
  }
  return outFileName;
}

void CreateTreeForMVATraining_Matching(TString TREENAME, bool matching = false)
{
  TFile *inf     = TFile::Open(TREENAME+".root");
  TTree *trIN    = (TTree*)inf->Get("boosted/events");
  //When opening a file from eos
  //TString outFileName = formOutFileName(TREENAME+".root");
  TFile *outf    = TFile::Open(TREENAME+"_Train.root","RECREATE");

  int nJets,nLeptons;
  float genEvtWeight;
  vector<bool>  *bit(0);
  vector<bool>  *matchedJet(0);
  vector<float> *pt(0), *tau4(0), *tau3(0),*tau2(0),*tau1(0),*mass(0);
  vector<float> *ecfB1N2(0), *ecfB1N3(0), *ecfB2N2(0), *ecfB2N3(0);
  vector<float> *jetMassSub0(0), *jetMassSub1(0);

  float pt_(-1),tau4_(-1), tau3_(-1),tau2_(-1),tau1_(-1),mass_(-1), jetMassSub0_(-1), jetMassSub1_(-1);
  float ecfB1N2_(-1), ecfB1N3_(-1), ecfB2N2_(-1), ecfB2N3_(-1);
  float jetDr_(-1);
  float JetPtOverSumPt_(-1);

  //matching info
  vector<float> *jetPhi(0), *jetEta(0);
  vector<int> *partonId(0), *partonMatchIdx(0);
  vector<float> *partonEta(0), *partonPhi(0), *partonMatchDR(0);
  std::vector<int> *addedIndexes = new std::vector<int>(0);

  //------- input tree --------------
  trIN->SetBranchAddress("nJets"          ,&nJets);
  trIN->SetBranchAddress("nLeptons"       ,&nLeptons);
  trIN->SetBranchAddress("jetPt"          ,&pt);
  trIN->SetBranchAddress("jetMassSoftDrop",&mass);
  trIN->SetBranchAddress("jetTau4"        ,&tau4);
  trIN->SetBranchAddress("jetTau3"        ,&tau3);
  trIN->SetBranchAddress("jetTau2"        ,&tau2);
  trIN->SetBranchAddress("jetTau1"        ,&tau1);
  trIN->SetBranchAddress("ecfB1N2"        ,&ecfB1N2);
  trIN->SetBranchAddress("ecfB1N3"        ,&ecfB1N3);
  trIN->SetBranchAddress("ecfB2N2"        ,&ecfB2N2);
  trIN->SetBranchAddress("ecfB2N3"        ,&ecfB2N3);

  trIN->SetBranchAddress("triggerBit"     ,&bit);
  trIN->SetBranchAddress("genEvtWeight"   ,&genEvtWeight);
  trIN->SetBranchAddress("jetMassSub0"    ,&jetMassSub0);
  trIN->SetBranchAddress("jetMassSub1"    ,&jetMassSub1);

  if(matching)
  {
    trIN->SetBranchAddress("partonMatchIdx" ,&partonMatchIdx);
    trIN->SetBranchAddress("partonMatchDR"  ,&partonMatchDR);
  }

  //------- output tree --------------
  TTree *trOUT = new TTree("events","events");

  trOUT->Branch("jetPt"           ,&pt_          ,"pt_/F");
  trOUT->Branch("jetMassSoftDrop" ,&mass_        ,"mass_/F");
  trOUT->Branch("jetTau4"         ,&tau4_        ,"tau4_/F");
  trOUT->Branch("jetTau3"         ,&tau3_        ,"tau3_/F");
  trOUT->Branch("jetTau2"         ,&tau2_        ,"tau2_/F");
  trOUT->Branch("jetTau1"         ,&tau1_        ,"tau1_/F");
  trOUT->Branch("ecfB1N2"         ,&ecfB1N2_     ,"ecfB1N2_/F");
  trOUT->Branch("ecfB1N3"         ,&ecfB1N3_     ,"ecfB1N3_/F");
  trOUT->Branch("ecfB2N2"         ,&ecfB2N2_     ,"ecfB2N2_/F");
  trOUT->Branch("ecfB2N3"         ,&ecfB2N3_     ,"ecfB2N3_/F");
  trOUT->Branch("JetPtOverSumPt"  ,&JetPtOverSumPt_, "JetPtOverSumPt_/F"); //jetPt_i / Sum of all jet pt's

  trOUT->Branch("jetMassSub0"     ,&jetMassSub0_ ,"jetMassSub0_/F");
  trOUT->Branch("jetMassSub1"     ,&jetMassSub1_ ,"jetMassSub1_/F");
  trOUT->Branch("jetDr"     ,&jetDr_ ,"jetDr_/F");
  int jetsWithProblem = 0;
  int totalTopJets = 0;
  int decade(0);
  int NN = trIN->GetEntries();
  std::vector<int> *jetMatchedIndexes = new std::vector<int>(2);
  std::vector<float> *jetMatchedDr = new std::vector<float>(0);
  //int NN = 10000;
  cout<<"Reading "<<NN<<" entries"<<endl;
  for(int iev=0;iev<NN;iev++) {
    double progress = 10.0*iev/(1.0*NN);
    int k = TMath::FloorNint(progress);
    if (k > decade)
      cout<<10*k<<" %"<<endl;
    decade = k;
    trIN->GetEntry(iev);
    //if (!(*bit)[2]) continue;
    //basic kinematic cuts for the event
    if (nJets < 2) continue;
    if (nLeptons > 0) continue;
    if ((*pt)[1] < 400) continue;
    if ((*mass[0]) < 50 && (*mass[0]) > 300 && (*mass[1]) < 50 && (*mass[2]) > 300) continue;

    float sumPt =0;
    for(int i =0; i<nJets; i++)
    {
      pt_ = -1;
      sumPt = sumPt + (*pt)[i];
    }
    //for every jet do the matching or store its properties
    for(int i=0; i<nJets; i++)
    {

      pt_          =-1;
      mass_        =-1;
      tau4_        =-1;
      tau3_        =-1;
      tau2_        =-1;
      tau1_        =-1;
      jetMassSub0_ =-1;
      jetMassSub1_ =-1;
      ecfB1N2_	   =-1;
      ecfB1N3_	   =-1;
      ecfB2N2_	   =-1;
      ecfB2N3_	   =-1;

      if(matching)
      {
        jetMatchedIndexes->clear();
        jetMatchedDr->clear();
        std::vector<int>::iterator it = std::find(partonMatchIdx->begin(), partonMatchIdx->end(), i);
        //get all entries that match our jet.
        while(it != partonMatchIdx->end())
        {
          int index = it - partonMatchIdx->begin();
          jetMatchedIndexes->push_back(index);
          jetMatchedDr->push_back((*partonMatchDR)[index]);
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

            jetDr_ 	 = dRmin;
            pt_    	 = (*pt)[(*partonMatchIdx)[indexMin]];
            mass_ 	 = (*mass)[(*partonMatchIdx)[indexMin]];
            tau4_    = (*tau4)[(*partonMatchIdx)[indexMin]];
            tau3_    = (*tau3)[(*partonMatchIdx)[indexMin]];
            tau2_ 	 = (*tau2)[(*partonMatchIdx)[indexMin]];
            tau1_ 	 = (*tau1)[(*partonMatchIdx)[indexMin]];
            ecfB1N2_ = (*ecfB1N2)[(*partonMatchIdx)[indexMin]];
            ecfB1N3_ = (*ecfB1N3)[(*partonMatchIdx)[indexMin]];
            ecfB2N2_ = (*ecfB2N2)[(*partonMatchIdx)[indexMin]];
            ecfB2N3_ = (*ecfB2N3)[(*partonMatchIdx)[indexMin]];

            JetPtOverSumPt_ = ((*pt)[(*partonMatchIdx)[indexMin]])/sumPt;

            jetMassSub0_ = (*jetMassSub0)[(*partonMatchIdx)[indexMin]];
            jetMassSub1_ = (*jetMassSub1)[(*partonMatchIdx)[indexMin]];
            trOUT->Fill();
          }
        }
      }
      else
      {
        pt_          =(*pt)[i];
        mass_        =(*mass)[i];
        tau4_        =(*tau4)[i];
        tau3_        =(*tau3)[i];
        tau2_        =(*tau2)[i];
        tau1_        =(*tau1)[i];
        ecfB1N2_     =(*ecfB1N2)[i];
        ecfB1N3_     =(*ecfB1N3)[i];
        ecfB2N2_     =(*ecfB2N2)[i];
        ecfB2N3_     =(*ecfB2N3)[i];
        jetMassSub0_ =(*jetMassSub0)[i];
        jetMassSub1_ =(*jetMassSub1)[i];

        JetPtOverSumPt_ = (*pt)[i]/sumPt;
        trOUT->Fill();
      }
    }
  }

  cout<<"Problematic events: "<<jetsWithProblem<<endl;
  cout<<"Total top jets: "<<totalTopJets<<std::endl;
  cout<<"==== Found "<<trOUT->GetEntries()<<" events ===="<<endl;
  outf->cd();
  trOUT->Write();
  outf->Close();
  inf->Close();
}
