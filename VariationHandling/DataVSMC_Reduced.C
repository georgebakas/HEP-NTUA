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
#include "TemplateConstantsDataVSMC.h"
TVector3 getBoostVector(TLorentzVector p4_1, TLorentzVector p4_2, TLorentzVector &p4CombinedVector);


void DataVSMC_Reduced(TString eos, TString process, TString data, TString year = "2016")
{
  cout<<data<<endl;
  initFilesMapping();
  initFilesXSEC();
  cout<<"YEAR: "<<year<<endl;
  TString file_name;
  float XSEC(1);
  if(data.EqualTo("true"))
  {
    file_name = dataFiles[year.Data()];
  }  
  else 
  {
    file_name = mcFiles[year.Data()][process.Data()];
    XSEC = XSECall[year.Data()][process.Data()];
    cout<<"XSEC: "<<XSEC<<endl;
  }
  cout<<"file_name: "<<file_name<<endl;
  cout<<"process: "<<process<<endl;
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
  

  gStyle->SetOptStat(0);
  const int NVAR =13;
  const int N_MVA = 50;

  cout<<"topTagger: "<<selMvaCut<<endl;
  cout<<"deepCSVFloat: "<<deepCSVFloat<<endl;

  //jetMassSub0_[ijet]

  TString varReco[NVAR]   = {"mTop","topTagger","jetTau3", "jetTau2", "jetTau1","jetMassSub0","jetMassSub1",
                              "ecfB1N2", "ecfB1N3","ecfB2N2", "ecfB2N3", "JetPtOverSumPt", "deltaPhi"};

  TFile *inf;
  float weights(0);

 //number of checks: one for Top Tagger and one for DAK8
 //initialize the required histograms
 TH1F *hCR_Leading[NVAR], *hCR_SubLeading[NVAR];
 TH1F *hSR_Leading[NVAR], *hSR_SubLeading[NVAR];

  int counter(0);
  TString temp = "";
  if (data.EqualTo("false"))
  {
    if (file_name.Contains("TTTo"))
      temp = "Signal/";
    else
      temp = "Bkg/";
  }
  
  cout<<eos<<temp<<file_name<<endl;
  TFile *file = TFile::Open(eos+temp+file_name); 
  TTree *trIN = (TTree*)file->Get("boosted/events");
  if(data.EqualTo("false"))
  {
    float norm = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
	  weights = XSEC/norm;
  }

  cout<<"LUMI: "<<LUMI<<endl;
  cout<<"LUMI_CR: "<<LUMI_CR<<endl;
  int decade(0);
  int NN = trIN->GetEntries();

  int nJets,nLeptons;
  float genEvtWeight(0);
  double bTagEvntWeight(0);
  vector<float> *jetPt(0),*jetTau3(0),*jetTau2(0),*jetTau1(0);
  vector<float> *jetMassSub0(0), *jetMassSub1(0);
  vector<float> *jetMassSoftDrop(0);

  vector<float> *jetTtag(0);
  vector<float> *ecfB1N2(0), *ecfB1N3(0), *ecfB2N2(0), *ecfB2N3(0);
  vector<bool> *bit = new vector<bool>;
  float mTTbarParton(0),mJJ(0), yJJ(0), ptJJ(0), yTTbarParton(0), ptTTbarParton(0);
  int  category(0);
  //matching info
  vector<float> *jetPhi(0), *jetEta(0), *jetY(0);
  vector<int> *partonId(0), *partonMatchIdx(0);

  vector<float> *partonEta(0), *partonPhi(0), *partonMatchDR(0),  *partonPt(0), *partonE(0), *partonMass(0), *deepAK8(0);
  std::vector<int> *addedIndexes = new std::vector<int>(0);
  std::vector<float> *jetBtagSub0(0), *jetBtagSub1(0);
  std::vector<float> *jetBtagSub0DCSVbb(0), *jetBtagSub1DCSVbb(0);
  std::vector<float> *jetBtagSub0DCSVbbb(0), *jetBtagSub1DCSVbbb(0);
  //------- input tree --------------
  trIN->SetBranchAddress("nJets"          ,&nJets);
  trIN->SetBranchAddress("nLeptons"       ,&nLeptons);
  trIN->SetBranchAddress("jetPt"          ,&jetPt);
  trIN->SetBranchAddress("jetEta"         ,&jetEta);
  trIN->SetBranchAddress("jetY"           ,&jetY);
  trIN->SetBranchAddress("jetPhi"         ,&jetPhi);
  trIN->SetBranchAddress("jetTau3"        ,&jetTau3);
  trIN->SetBranchAddress("jetTau2"        ,&jetTau2);
  trIN->SetBranchAddress("jetTau1"        ,&jetTau1);
  trIN->SetBranchAddress("triggerBit"     ,&bit);
  trIN->SetBranchAddress("genEvtWeight"   ,&genEvtWeight);
  trIN->SetBranchAddress("bTagEvntWeight"  ,&bTagEvntWeight);
  trIN->SetBranchAddress("jetMassSub0"    ,&jetMassSub0);
  trIN->SetBranchAddress("jetMassSub1"    ,&jetMassSub1);
  trIN->SetBranchAddress("mJJ"        ,&mJJ);
  trIN->SetBranchAddress("yJJ"        ,&yJJ);
  trIN->SetBranchAddress("ptJJ"         ,&ptJJ);
  trIN->SetBranchAddress("jetBtagSub0"    ,&jetBtagSub0);
  //trIN->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1);
  trIN->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
  trIN->SetBranchAddress("category"       ,&category);
  trIN->SetBranchAddress("deepAK8Tagger"  ,&deepAK8);
  trIN->SetBranchAddress("jetTtagCategory",&jetTtag);
  trIN->SetBranchAddress("ecfB1N2",&ecfB1N2);
  trIN->SetBranchAddress("ecfB1N3",&ecfB1N3);
  trIN->SetBranchAddress("ecfB2N2",&ecfB2N2);
  trIN->SetBranchAddress("ecfB2N3",&ecfB2N3);

  //deepCSV
  trIN->SetBranchAddress("jetBtagSub0DCSVbb" ,&jetBtagSub0DCSVbb);
  trIN->SetBranchAddress("jetBtagSub1DCSVbb" ,&jetBtagSub1DCSVbb);
  trIN->SetBranchAddress("jetBtagSub0DCSVbbb",&jetBtagSub0DCSVbbb);
  trIN->SetBranchAddress("jetBtagSub1DCSVbbb",&jetBtagSub1DCSVbbb);

  if(file_name.Contains("TTTo"))
  {
    trIN->SetBranchAddress("mTTbarParton" ,&mTTbarParton);
    trIN->SetBranchAddress("yTTbarParton" ,&yTTbarParton);
    trIN->SetBranchAddress("ptTTbarParton"  ,&ptTTbarParton);
    trIN->SetBranchAddress("partonPt"     ,&partonPt);
    trIN->SetBranchAddress("partonEta"    ,&partonEta);
    trIN->SetBranchAddress("partonMass"     ,&partonMass);
    trIN->SetBranchAddress("partonMatchDR"  ,&partonMatchDR);
    trIN->SetBranchAddress("partonMatchIdx" ,&partonMatchIdx);
    trIN->SetBranchAddress("partonPhi"      ,&partonPhi);
  }

  float xParton(0);
  float xReco_leading(0), xReco_subleading(0);
  std::vector<float> xRecoAll_Leading(0);
  std::vector<float> xRecoAll_SubLeading(0);
  //book the histograms
  //histograms for Signal/QCD in CR
  float xMin[NVAR] = {120,-1,0,0,0, 0,0,0,0, 0,0,0,-3};
  float xMax[NVAR] = {220,1,0.3,0.4, 0.6,220, 120,0.5,5, 0.5,5,1,3};

  for(int ivar =0; ivar< NVAR; ivar++)
  {
    hCR_Leading[ivar] = new TH1F(TString::Format("hCR_Leading_%s_%s", "tTagger", varReco[ivar].Data()), 
                            TString::Format("hCR_Leading_%s_%s","tTagger", varReco[ivar].Data()), N_MVA, xMin[ivar],xMax[ivar]);
    hCR_SubLeading[ivar] = new TH1F(TString::Format("hCR_SubLeading_%s_%s", "tTagger", varReco[ivar].Data()), 
                            TString::Format("hCR_SubLeading_%s_%s","tTagger", varReco[ivar].Data()), N_MVA, xMin[ivar],xMax[ivar]);
    hSR_Leading[ivar] = new TH1F(TString::Format("hSR_Leading_%s_%s", "tTagger", varReco[ivar].Data()), 
                            TString::Format("hSR_Leading_%s_%s","tTagger", varReco[ivar].Data()), N_MVA, xMin[ivar],xMax[ivar]);
    hSR_SubLeading[ivar] = new TH1F(TString::Format("hSR_SubLeading_%s_%s", "tTagger", varReco[ivar].Data()), 
                            TString::Format("hSR_SubLeading%s_%s","tTagger", varReco[ivar].Data()), N_MVA, xMin[ivar],xMax[ivar]);
  }


  cout<<"Reading "<<NN<<" entries"<<endl;
  for(int iev=4;iev<NN;iev++)
  {
    double progress = 10.0*iev/(1.0*NN);
    int k = TMath::FloorNint(progress);
    if (k > decade)
      cout<<10*k<<" %"<<endl;
    decade = k;
    trIN->GetEntry(iev);

  xRecoAll_Leading.clear();
  xRecoAll_SubLeading.clear();
  bool partonCuts, recoCuts, massCut, tTaggerCut, triggerCR, triggerSR;
  bool deepCSV, revertBtagDeepCSV;
  bool btagCut, revertBtag;

  if (nJets >1)
  {
    float dCSVScoreSub0[2], dCSVScoreSub1[2];
    dCSVScoreSub0[0] = (*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0];
    dCSVScoreSub0[1] = (*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1];
    dCSVScoreSub1[0] = (*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0];
    dCSVScoreSub1[1] = (*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1];

    recoCuts   = fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1]) <2.4 && (*jetPt)[0] > 450 && (*jetPt)[1] > 400 &&  mJJ > 1000 && nLeptons==0;
    triggerSR  = (*bit)[triggerFloat];
    triggerCR  = (*bit)[triggerFloatCR];
    massCut    = (*jetMassSoftDrop)[0] > 120 && (*jetMassSoftDrop)[0] < 220 && (*jetMassSoftDrop)[1] > 120 && (*jetMassSoftDrop)[1] < 220;
    //tTaggerCut = (*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut;
    tTaggerCut = true;
    //2 btag category with csvv2 and deepCSV
    deepCSV    = (((*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0])> deepCSVFloat) &&
           (((*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1])> deepCSVFloat);

    //0 btag category with deepCSV
    revertBtagDeepCSV = (dCSVScoreSub0[0] < deepCSVFloat &&  dCSVScoreSub1[0] < deepCSVFloat) && (dCSVScoreSub0[1] < deepCSVFloat && dCSVScoreSub1[1] < deepCSVFloat);

    float JetPtOverSumPt[2];
    float sumPt = 0;
    for(int ijet = 0; ijet<nJets; ijet++)
    {
      sumPt = sumPt + (*jetPt)[ijet];
    }

    JetPtOverSumPt[0] = (*jetPt)[0]/sumPt;
    JetPtOverSumPt[1] = (*jetPt)[1]/sumPt;
    //cout<<"ok"<<endl;
    float deltaPhi = fabs((*jetPhi)[0] - (*jetPhi)[1]);
    //cout<<deltaPhi<<endl;
    xRecoAll_Leading.push_back((*jetMassSoftDrop)[0]);
    xRecoAll_Leading.push_back((*jetTtag)[0]);
    //xRecoAll_Leading.push_back((*jetTau3)[0]);
    //xRecoAll_Leading.push_back((*jetTau2)[0]);
    //xRecoAll_Leading.push_back((*jetTau1)[0]);
    //xRecoAll_Leading.push_back((*jetMassSub0)[0]);
    //xRecoAll_Leading.push_back((*jetMassSub1)[0]);
    //xRecoAll_Leading.push_back((*ecfB1N2)[0]);
    //xRecoAll_Leading.push_back((*ecfB1N3)[0]);
    //xRecoAll_Leading.push_back((*ecfB2N2)[0]);
    //xRecoAll_Leading.push_back((*ecfB2N3)[0]);
    xRecoAll_Leading.push_back(JetPtOverSumPt[0]);
    xRecoAll_Leading.push_back(deltaPhi);

    xRecoAll_SubLeading.push_back((*jetMassSoftDrop)[1]);
    xRecoAll_SubLeading.push_back((*jetTtag)[1]);
    //xRecoAll_SubLeading.push_back((*jetTau3)[1]);
    //xRecoAll_SubLeading.push_back((*jetTau2)[1]);
    //xRecoAll_SubLeading.push_back((*jetTau1)[1]);
    //xRecoAll_SubLeading.push_back((*jetMassSub0)[1]);
    //xRecoAll_SubLeading.push_back((*jetMassSub1)[1]);
    //xRecoAll_SubLeading.push_back((*ecfB1N2)[1]);
    //xRecoAll_SubLeading.push_back((*ecfB1N3)[1]);
    //xRecoAll_SubLeading.push_back((*ecfB2N2)[1]);
    //xRecoAll_SubLeading.push_back((*ecfB2N3)[1]);
    xRecoAll_SubLeading.push_back(JetPtOverSumPt[1]);
    xRecoAll_SubLeading.push_back(deltaPhi);

  btagCut = deepCSV;
  revertBtag = revertBtagDeepCSV;

   //cout<<"------"<<endl;
   for(int ivar = 0; ivar <xRecoAll_Leading.size(); ivar ++)
   {

     xReco_leading = xRecoAll_Leading[ivar];
     xReco_subleading = xRecoAll_SubLeading[ivar];
     if(data.EqualTo("true")){
       genEvtWeight =1;
       bTagEvntWeight = 1;
     }

     //Signal Region with tTagger
     if(recoCuts && btagCut && massCut && tTaggerCut && triggerSR)
     {
      hSR_Leading[ivar]->Fill(xReco_leading,genEvtWeight*bTagEvntWeight);
      hSR_SubLeading[ivar]->Fill(xReco_subleading,genEvtWeight*bTagEvntWeight);
     }
     //Control Region with tTagger
     if(recoCuts && revertBtag && massCut && tTaggerCut && triggerCR)
     {
      hCR_Leading[ivar]->Fill(xReco_leading,genEvtWeight*bTagEvntWeight);
      hCR_SubLeading[ivar]->Fill(xReco_subleading,genEvtWeight*bTagEvntWeight);
     }

  }//----end of nJets
  } //---end of event loop

  }//----end of fileSize loop

  TH1F *hCR_Leading_Clone[NVAR];
  TH1F *hCR_SubLeading_Clone[NVAR];
  TH1F *hSR_Leading_Clone[NVAR];
  TH1F *hSR_SubLeading_Clone[NVAR];

  // loop on variables 
  for(int ivar= 0; ivar<NVAR; ivar++)
  {
    TString varNameReco = varReco[ivar];
    
    hCR_Leading_Clone[ivar]=(TH1F*)hCR_Leading[ivar]->Clone(TString::Format("hCR_%s_Clone","tTagger"));
    hSR_Leading_Clone[ivar]=(TH1F*)hSR_Leading[ivar]->Clone(TString::Format("hSR_%s_Clone","tTagger"));
    hCR_SubLeading_Clone[ivar]=(TH1F*)hCR_SubLeading[ivar]->Clone(TString::Format("hCR_%s_Clone","tTagger"));
    hSR_SubLeading_Clone[ivar]=(TH1F*)hSR_SubLeading[ivar]->Clone(TString::Format("hSR_%s_Clone","tTagger"));

    if(data.EqualTo("false"))
    {
      hCR_Leading_Clone[ivar]->Scale(weights*LUMI_CR); //this is 0 btagged (CR)
      hCR_SubLeading_Clone[ivar]->Scale(weights*LUMI_CR); //this is 0 btagged (CR)
      hSR_Leading_Clone[ivar]->Scale(weights*LUMI); //this is 2 btagged (SR)
      hSR_SubLeading_Clone[ivar]->Scale(weights*LUMI); //this is 2 btagged (SR)

      hCR_Leading[ivar]->Scale(weights); //this is CR
      hCR_SubLeading[ivar]->Scale(weights); //this is CR
      hSR_Leading[ivar]->Scale(weights); //this is Signal region
      hSR_SubLeading[ivar]->Scale(weights); //this is Signal region
    }

  }

  TFile *outFile = new TFile(TString::Format("%s/Nominal/DataVSMC/%s.root",
                              year.Data(), process.Data()), "RECREATE"); 

  for(int ivar = 0; ivar<NVAR; ivar++)
  {
    TString varNameReco = varReco[ivar];
    outFile->cd();
    hSR_Leading[ivar]->Write(TString::Format("hWt_%s_2btag_leading", varNameReco.Data()));
    hSR_SubLeading[ivar]->Write(TString::Format("hWt_%s_2btag_subleading", varNameReco.Data()));

    hCR_Leading[ivar]->Write(TString::Format("hWt_%s_0btag_leading", varNameReco.Data()));
    hCR_SubLeading[ivar]->Write(TString::Format("hWt_%s_0btag_subleading", varNameReco.Data()));


    hSR_Leading_Clone[ivar]->Write(TString::Format("hWt_%s_2btag_expYield_leading", varNameReco.Data()));
    hSR_SubLeading_Clone[ivar]->Write(TString::Format("hWt_%s_2btag_expYield_subleading", varNameReco.Data()));

    hCR_Leading_Clone[ivar]->Write(TString::Format("hWt_%s_0btag_expYield_leading", varNameReco.Data()));
    hCR_SubLeading_Clone[ivar]->Write(TString::Format("hWt_%s_0btag_expYield_subleading", varNameReco.Data()));
  }
}
