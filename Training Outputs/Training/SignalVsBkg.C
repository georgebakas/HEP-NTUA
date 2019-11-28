#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/Tools.h"

#include "TFile.h"
#include "TTree.h"
#include "TCut.h" 
#include "TH1F.h"


using namespace TMVA;
using namespace TMath;

void SignalVsBkg()
{
  char name[1000];
  //float XSEC[6] = {3.67e+5,2.94e+4,6.524e+03,1.064e+03,121.5,2.542e+01};
  float XSEC[5] = {3.67e+5,2.94e+4,1.064e+03,121.5,2.542e+01};
  float NORM[5];
  TFile *bkgSrcNorm[5];//*bkgSrcNorm[6];

  //bkgSrc[0] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1_legacy2016_Train.root");
  //bkgSrc[1] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1_legacy2016_Train.root");
  //bkgSrc[2] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1_legacy2016_Train.root");
  //bkgSrc[3] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1_legacy2016_Train.root");
  //bkgSrc[4] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1_legacy2016_Train.root");
  //bkgSrc[5] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1_legacy2016_Train.root");

  bkgSrcNorm[0] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  bkgSrcNorm[1] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  //bkgSrcNorm[2] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  bkgSrcNorm[2] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  bkgSrcNorm[3] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  bkgSrcNorm[4] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");

  TFile *sigSrc = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Signal/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_legacy2016_Copy.root");
  //TFile *sigSrc = TFile::Open("/eos/cms/store/user/gbakas/ttbar/topTagger/April19/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_Train.root.root");
  TTree *sigTree = (TTree*)sigSrc->Get("events");

  vector<float> *jetPt(0),*tau3(0),*tau2(0),*tau1(0);
  vector<float> *jetMassSub0(0), *jetMassSub1(0);
  vector<float> *jetMassSoftDrop(0);

  float mva(0);
  float mTTbarParton(0),mJJ(0);
  int  category(0), nJets(0);
  //matching info
  vector<float> *jetPhi(0), *jetEta(0);
  vector<float> *partonEta(0), *partonPhi(0),  *partonPt(0);
  std::vector<float> *jetBtagSub0(0), *jetBtagSub1(0);

    //------- input tree --------------
  sigTree->SetBranchAddress("nJets"          ,&nJets);
  //sigTree->SetBranchAddress("nLeptons"       ,&nLeptons);
  sigTree->SetBranchAddress("jetPt"          ,&jetPt);
  sigTree->SetBranchAddress("jetEta"         ,&jetEta);
  sigTree->SetBranchAddress("jetPhi"         ,&jetPhi);
  sigTree->SetBranchAddress("jetTau3"        ,&tau3);
  sigTree->SetBranchAddress("jetTau2"        ,&tau2);
  sigTree->SetBranchAddress("jetTau1"        ,&tau1);
  //sigTree->SetBranchAddress("triggerBit"     ,&bit);
  //sigTree->SetBranchAddress("genEvtWeight"   ,&genEvtWeight);
  sigTree->SetBranchAddress("jetMassSub0"    ,&jetMassSub0);
  sigTree->SetBranchAddress("jetMassSub1"    ,&jetMassSub1);
  sigTree->SetBranchAddress("mTTbarParton"   ,&mTTbarParton);
  sigTree->SetBranchAddress("jetBtagSub0"    ,&jetBtagSub0);
  sigTree->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1);
  sigTree->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
  sigTree->SetBranchAddress("mJJ"            ,&mJJ);
  sigTree->SetBranchAddress("partonPt"       ,&partonPt);
  sigTree->SetBranchAddress("partonEta"      ,&partonEta);
  sigTree->SetBranchAddress("partonPhi"      ,&partonPhi);
  sigTree->SetBranchAddress("mva"            ,&mva);
  sigTree->SetBranchAddress("category"       ,&category);

  
  
  //TFile *outf = new TFile("boosted_MVA.root","RECREATE");

  int sigEvents = sigTree->GetEntries();
  TH1F *hSig_0 = new TH1F("signal btagging 0", "signal btagging 0", 50, 0 , 2);
  TH1F *hSig_1 = new TH1F("signal btagging 1", "signal btagging 2", 50, 0 , 2);

  for(int iev =0; iev <=sigEvents; iev++)
  {
    sigTree->GetEntry(iev);

    //without any selection info
    for(int nj =0; nj <nJets; nj++)
    {
      hSig_0->Fill((*jetBtagSub0)[nj]);
      hSig_1->Fill((*jetBtagSub1)[nj]);  
    }
  

  }

  

cout<<"entering bkg"<<endl;
//bkg----------------------------------------------------------
  vector<float> *jetPtBkg(0),*tau3Bkg(0),*tau2Bkg(0),*tau1Bkg(0);
  vector<float> *jetMassSub0Bkg(0), *jetMassSub1Bkg(0);
  vector<float> *jetMassSoftDropBkg(0);

  float mTTbarPartonBkg(0),mJJBkg(0);
  int  categoryBkg(0), nJetsBkg(0);
  //matching info
  vector<float> *jetPhiBkg(0), *jetEtaBkg(0);
  vector<float> *partonEtaBkg(0), *partonPhiBkg(0),  *partonPtBkg(0);
  std::vector<float> *jetBtagSub0Bkg(0), *jetBtagSub1Bkg(0);

 TTree *bkgTree[5];
 TH1F  *hBkg_0[5], *hBkg_1[5];
 for(int nBkg = 0; nBkg <5; nBkg++)
 {

  cout<<"crashed here"<<endl;
    bkgTree[nBkg] = (TTree*)bkgSrcNorm[nBkg]->Get("boosted/events"); 
      
    bkgTree[nBkg]->SetBranchAddress("nJets"          ,&nJetsBkg);
    //bkgTree[nBkg]->SetBranchAddress("jetPt"          ,&jetPtBkg);
    //bkgTree[nBkg]->SetBranchAddress("jetEta"         ,&jetEtaBkg);
    //bkgTree[nBkg]->SetBranchAddress("jetPhi"         ,&jetPhiBkg);
    //bkgTree[nBkg]->SetBranchAddress("jetTau3"        ,&tau3Bkg);
    //bkgTree[nBkg]->SetBranchAddress("jetTau2"        ,&tau2Bkg);
    //bkgTree[nBkg]->SetBranchAddress("jetTau1"        ,&tau1Bkg);
    //bkgTree[nBkg]->SetBranchAddress("jetMassSub0"    ,&jetMassSub0Bkg);
    //bkgTree[nBkg]->SetBranchAddress("jetMassSub1"    ,&jetMassSub1Bkg);
    //bkgTree[nBkg]->SetBranchAddress("mTTbarParton"   ,&mTTbarPartonBkg);
    bkgTree[nBkg]->SetBranchAddress("jetBtagSub0"    ,&jetBtagSub0Bkg);
    bkgTree[nBkg]->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1Bkg);
    //bkgTree[nBkg]->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
    //bkgTree[nBkg]->SetBranchAddress("mJJ"            ,&mJJBkg);
   // bkgTree[nBkg]->SetBranchAddress("partonPt"       ,&partonPt);
    //bkgTree[nBkg]->SetBranchAddress("partonEta"      ,&partonEta);
    //bkgTree[nBkg]->SetBranchAddress("partonPhi"      ,&partonPhi);
    //bkgTree[nBkg]->SetBranchAddress("mva"            ,&mva);
   // bkgTree[nBkg]->SetBranchAddress("category"       ,&categoryBkg);
  cout<<"crashed here"<<endl;

    hBkg_0[nBkg] = new TH1F (TString::Format("bkg %d 0", nBkg+1),TString::Format("bkg %d 0", nBkg+1), 50,0,2);
    hBkg_1[nBkg] = new TH1F (TString::Format("bkg %d 1", nBkg+1),TString::Format("bkg %d 1", nBkg+1), 50,0,2);

    for(int bev=0; bev<=bkgTree[nBkg]->GetEntries(); bev++)
    {
      bkgTree[nBkg]->GetEntry(bev);

      //without any selection info
      for(int nj =0; nj <nJetsBkg; nj++)
      {
        hBkg_0[nBkg]->Fill((*jetBtagSub0Bkg)[nj]);
        hBkg_1[nBkg]->Fill((*jetBtagSub1Bkg)[nj]);  
      }
    }

 }
  float LUMI = 35500.;
  float NORMsig = ((TH1F*)sigSrc->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
  //normalize signal
  //hSig_1->Scale((832./NORMsig) *LUMI);
  //hSig_0->Scale((832./NORMsig) *LUMI);

  //normalize bkg
  for(int k=0;k<5;k++) {
    NORM[k] = ((TH1F*)bkgSrcNorm[k]->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
    //hBkg_0[k]->Scale((XSEC[k]/NORM[k])*LUMI);
    //hBkg_1[k]->Scale((XSEC[k]/NORM[k])*LUMI);
    if(k !=0)
    {
      hBkg_1[0]->Add(hBkg_1[k]);
      hBkg_0[0]->Add(hBkg_0[k]);
    }
  }

  hBkg_0[0]->Scale(1./hBkg_0[0]->Integral());
  hBkg_1[0]->Scale(1./hBkg_1[0]->Integral());
  hSig_1->Scale(1./hSig_1->Integral());
  hSig_0->Scale(1./hSig_0->Integral());

  hSig_0->GetXaxis()->SetTitle("jetBtagSub0");
  hSig_1->GetXaxis()->SetTitle("jetBtagSub1");
  hBkg_0[0]->GetXaxis()->SetTitle("jetBtagSub0");
  hBkg_1[0]->GetXaxis()->SetTitle("jetBtagSub1");


  TCanvas *can0 = new TCanvas("can0", "can0", 900, 600);
  
  hBkg_0[0]->SetLineColor(kRed);
  hBkg_0[0]->Draw();
  hSig_0->Draw("same");

  TCanvas *can1 = new TCanvas("can1", "can1", 900, 600);
  hBkg_1[0]->SetLineColor(kRed);
  hBkg_1[0]->Draw();
  hSig_1->Draw("same");




  //outf->Close();
}
