#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"

#include "TemplateConstants.h"
using std::cin;
using std::cout;
using std::endl;


std::vector<TString> listOfFiles;
std::vector<float> XSEC;
std::vector<TString> histoNames;
float LUMI;
float LUMI_CR;
float deepCSVFloat;
TString eosPath;
TString year;
int selection;

void initFileNames()
{

  if(selection ==0) //data
  {
    eosPath = TString::Format("%s%s/",eosDataPath.Data(), year.Data());
    listOfFiles.push_back(dataFiles[year.Data()]);
  }
  else if(selection ==1) //signal ttbar mc mtt
  {
    eosPath = TString::Format("%s%s/Signal/",eosPathMC.Data(), year.Data());
    cout<<eosPath<<endl;
    cout<<mttFiles[year.Data()]["700-1000"]<<endl;
    listOfFiles.push_back(mttFiles[year.Data()]["700-1000"]);
    listOfFiles.push_back(mttFiles[year.Data()]["1000-Inf"]);
  }
  else if(selection ==2) //bkg mc
  {
    eosPath = TString::Format("%s%s/Bkg/",eosPathMC.Data(), year.Data());
    listOfFiles.push_back(qcdBkgFiles[year.Data()]["300-500"]);
    listOfFiles.push_back(qcdBkgFiles[year.Data()]["500-700"]);
    listOfFiles.push_back(qcdBkgFiles[year.Data()]["700-1000"]);
    listOfFiles.push_back(qcdBkgFiles[year.Data()]["1000-1500"]);
    listOfFiles.push_back(qcdBkgFiles[year.Data()]["1500-2000"]);
    listOfFiles.push_back(qcdBkgFiles[year.Data()]["2000-Inf"]);
  }
  else if(selection ==3) //subdominant bkgs
  {
    eosPath = TString::Format("%s%s/Bkg/",eosPathMC.Data(), year.Data());
    if(year.EqualTo("2016"))
    {
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["DY"]);
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["WJetsToQQ_HT180"]);
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_tW_top_5f"]);
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_tW_antitop_5f"]);
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_t-channel_top_4f"]);
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_t-channel_antitop_4f"]);
    }
    else
    {
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["DY"]);
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["WJetsToQQ_HT400"]);
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["WJetsToQQ_HT600"]);
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_tW_top_5f"]);
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_tW_antitop_5f"]);
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_t-channel_top_4f"]);
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_t-channel_antitop_4f"]);
      if(year.EqualTo("2018"))
      {
        listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_t-channel_top_5f"]);
        listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_t-channel_antitop_5f"]);
      }
    }
  }
  else if(selection ==4) //signal ttbar mc nominal
  {
    cout<<"nominal!!!"<<endl;
    eosPath = TString::Format("%s%s/Signal/",eosPathMC.Data(), year.Data());
    cout<<eosPath<<endl;
    if(year.EqualTo("2016")) listOfFiles.push_back(ttNominalFiles[year.Data()]["TTNominal"]);
    else
    {
      listOfFiles.push_back(ttNominalFiles[year.Data()]["TTHadronic_0"]);
      listOfFiles.push_back(ttNominalFiles[year.Data()]["TTSemiLeptonic_0"]);
      listOfFiles.push_back(ttNominalFiles[year.Data()]["TTTo2L2Nu_0"]);
    }
  }
}

void initXsections()
{
  if(selection ==1) //signal ttbar mc
  {
    XSEC.push_back(mttXSEC[year.Data()]["700-1000"]);
    XSEC.push_back(mttXSEC[year.Data()]["1000-Inf"]);
  }
  else if(selection ==2) //bkg mc
  {
    XSEC.push_back(qcdBkgXSEC[year.Data()]["300-500"]);
    XSEC.push_back(qcdBkgXSEC[year.Data()]["500-700"]);
    XSEC.push_back(qcdBkgXSEC[year.Data()]["700-1000"]);
    XSEC.push_back(qcdBkgXSEC[year.Data()]["1000-1500"]);
    XSEC.push_back(qcdBkgXSEC[year.Data()]["1500-2000"]);
    XSEC.push_back(qcdBkgXSEC[year.Data()]["2000-Inf"]);
  }
  else if(selection ==3) //subdominant bkgs
  {
    if(year.EqualTo("2016"))
    {
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["DY"]);
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["WJetsToQQ_HT180"]);
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_tW_top_5f"]);
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_tW_antitop_5f"]);
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_t-channel_top_4f"]);
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_t-channel_antitop_4f"]);
    }
    else
    {
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["DY"]);
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["WJetsToQQ_HT400"]);
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["WJetsToQQ_HT600"]);
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_tW_top_5f"]);
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_tW_antitop_5f"]);
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_t-channel_top_4f"]);
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_t-channel_antitop_4f"]);
      if(year.EqualTo("2018"))
      {
        XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_t-channel_top_5f"]);
        XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_t-channel_antitop_5f"]);
      }
    }
  }
  else if(selection ==4)
  {
    if(year.EqualTo("2016")) XSEC.push_back(ttNominalXSEC[year.Data()]["TTNominal"]);
    else
    {
      XSEC.push_back(ttNominalXSEC[year.Data()]["TTHadronic_0"]);
      XSEC.push_back(ttNominalXSEC[year.Data()]["TTSemiLeptonic_0"]);
      XSEC.push_back(ttNominalXSEC[year.Data()]["TTTo2L2Nu_0"]);
    }
  }

}

void initHistoNames()
{

  if(selection ==0) histoNames.push_back(TString::Format("Data_%s", year.Data()));
  else if(selection ==1)
  {
    histoNames.push_back("Signal_histo_Mtt_700_1000");
    histoNames.push_back("Signal_histo_Mtt_1000_Inf");
  }
  else if (selection ==2)
  {
    histoNames.push_back("QCD_histo_300_500");
    histoNames.push_back("QCD_histo_500_700");
    histoNames.push_back("QCD_histo_700_1000");
    histoNames.push_back("QCD_histo_1000_1500");
    histoNames.push_back("QCD_histo_1500_2000");
    histoNames.push_back("QCD_histo_2000_Inf");
  }
  else if(selection ==3)
  {

    if(year.EqualTo("2016"))
    {
       histoNames.push_back("DYJetsToQQ_HT180");
       histoNames.push_back("WJetsToQQ_HT180");
       histoNames.push_back("ST_tW_top_5f_inclusiveDecays");
       histoNames.push_back("ST_tW_antitop_5f_inclusiveDecays");
       histoNames.push_back("ST_t-channel_top_4f_inclusiveDecays");
       histoNames.push_back("ST_t-channel_antitop_4f_inclusiveDecays");
    }
    else
    {
       histoNames.push_back("DYJetsToQQ_HT180");
       histoNames.push_back("WJetsToQQ_HT400");
       histoNames.push_back("WJetsToQQ_HT600");
       histoNames.push_back("ST_tW_top_5f_inclusiveDecays");
       histoNames.push_back("ST_tW_antitop_5f_inclusiveDecays");
       histoNames.push_back("ST_t-channel_top_4f_inclusiveDecays");
       histoNames.push_back("ST_t-channel_antitop_4f_inclusiveDecays");
      if(year.EqualTo("2018"))
      {
        histoNames.push_back("ST_t-channel_top_5f");
        histoNames.push_back("ST_t-channel_antitop_5f");
      }
    }
  }
  else if(selection ==4)
  {
    if(year.EqualTo("2016")) histoNames.push_back("TTNominal");
    else
    {
      histoNames.push_back("TTHadronic");
      histoNames.push_back("TTSemiLeptonic");
      histoNames.push_back("TTTo2L2Nu");
    }
  }

}

void initGlobals()
{
  initFileNames();
  initXsections();
  initHistoNames();
}

void BTagSFHistograms(TString y="2016", int sel = 0, bool isLoose=false)
{
  year =y;
  initFilesMapping(isLoose);
  deepCSVFloat = deepCSVFloatMap[year.Data()];
  selection = sel;
  LUMI = luminosity[year.Data()];
  LUMI_CR = luminosityCR[year.Data()];
  initGlobals();
  gStyle->SetOptStat(0);
  const int NVAR =1;
  const int N_MVA = 100;

  float selMvaCut=topTaggerCuts[year];

  cout<<"triggerSRConst[year.Data()]]: "<<triggerSRConst[year.Data()]<<endl;
  cout<<"triggerCRConst[year.Data()]]: "<<triggerCRConst[year.Data()]<<endl;
  cout<<"topTagger: "<<selMvaCut<<endl;
  cout<<"deepCSVFloat: "<<deepCSVFloat<<endl;

  //jetMassSub0_[ijet]

  TString varReco[NVAR] = {"bTagEvntWeight"};

  int fileSize = listOfFiles.size();
  TFile *inf;
  vector<float> weights(0);

 //number of checks: one for Top Tagger and one for DAK8
 //initialize the required histograms
 TH1F *hCR[listOfFiles.size()][NVAR];
 TH1F *hSR[listOfFiles.size()][NVAR];

 for(int f=0; f<listOfFiles.size(); f++)
 {
  int counter(0);
  cout<<"Entering "<<eosPath+listOfFiles[f]<<endl;
  inf = TFile::Open(eosPath+listOfFiles[f]);
  TTree *trIN    = (TTree*)inf->Get("boosted/events");
  //cout<<"XSEC: "<<XSEC[f]<<endl;
  if(selection != 0)
  {
    float NORM = ((TH1F*)inf->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
    weights.push_back(XSEC[f]/NORM);
  }
  cout<<"file: "<<eosPath+listOfFiles[f]<<endl;
  //cout<<"weight: "<<weights[f]<<endl;
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

  if(selection == 1 || selection ==4)
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
  float xReco(0);
  std::vector<float> xRecoAll(0);
  //book the histograms
  //histograms for Signal/QCD in CR
  float xMin[NVAR] = {0.0};
  float xMax[NVAR] = {1.2};

  for(int ivar =0; ivar< NVAR; ivar++)
  {
    hCR[f][ivar] = new TH1F(TString::Format("hCR_%s_%s",histoNames[f].Data(),varReco[ivar].Data()), TString::Format("hCR_%s_%s",histoNames[f].Data(),varReco[ivar].Data()), N_MVA, xMin[ivar],xMax[ivar]);
    hSR[f][ivar] = new TH1F(TString::Format("hSR_%s_%s",histoNames[f].Data(),varReco[ivar].Data()), TString::Format("hSR_%s_%s",histoNames[f].Data(),varReco[ivar].Data()), N_MVA, xMin[ivar],xMax[ivar]);

  }
  //for matching

  std::vector<int> *jetMatchedIndexes = new std::vector<int>(0);
  std::vector<float> *jetMatchedDr = new std::vector<float>(0);
  std::vector<float> *eta_ = new std::vector<float>(0);
  std::vector<float> *y_ = new std::vector<float>(0);
  std::vector<float> *phi_ = new std::vector<float>(0);
  std::vector<float> *mass_ = new std::vector<float>(0);
  std::vector<float> *pt_ = new std::vector<float>(0);

  std::vector<float> *jetTtag_ = new std::vector<float>(0);
  std::vector<float> *jetTau3_ = new std::vector<float>(0);
  std::vector<float> *jetTau2_ = new std::vector<float>(0);
  std::vector<float> *jetTau1_ = new std::vector<float>(0);
  std::vector<float> *jetMassSub0_ = new std::vector<float>(0);
  std::vector<float> *jetMassSub1_ = new std::vector<float>(0);
  std::vector<float> *jetBtagSub0_ = new std::vector<float>(0);
  std::vector<float> *jetBtagSub1_ = new std::vector<float>(0);

  std::vector<float> *ecfB1N2_ = new std::vector<float>(0);
  std::vector<float> *ecfB1N3_ = new std::vector<float>(0);
  std::vector<float> *ecfB2N2_ = new std::vector<float>(0);
  std::vector<float> *ecfB2N3_ = new std::vector<float>(0);

  std::vector<float> *partonPt_ = new std::vector<float>(0);
  std::vector<float> *partonY_ = new std::vector<float>(0);
  std::vector<float> *partonMass_ = new std::vector<float>(0);
  std::vector<float> *partonPhi_ = new std::vector<float>(0);
  std::vector<float> *partonEta_ = new std::vector<float>(0);

  std::vector<float> *jetBtagSub0DCSVbb_ = new std::vector<float>(0);
  std::vector<float> *jetBtagSub1DCSVbb_ = new std::vector<float>(0);
  std::vector<float> *jetBtagSub0DCSVbbb_ = new std::vector<float>(0);
  std::vector<float> *jetBtagSub1DCSVbbb_ = new std::vector<float>(0);
  float jetDr_(0);

  cout<<"Reading "<<NN<<" entries"<<endl;
  for(int iev=4;iev<NN;iev++)
  {
    double progress = 10.0*iev/(1.0*NN);
    int k = TMath::FloorNint(progress);
    if (k > decade)
      cout<<10*k<<" %"<<endl;
    decade = k;
    trIN->GetEntry(iev);
  int isMatched =0;
  eta_->clear();
  y_->clear();
    mass_->clear();
    pt_->clear();
    phi_->clear();
    jetBtagSub0_->clear();
    //jetBtagSub1_->clear();

    jetTtag_->clear();
    ecfB1N2_->clear();
    ecfB1N3_->clear();
    ecfB2N2_->clear();
    ecfB2N3_->clear();

    jetTau3_->clear();
    jetTau2_->clear();
    jetTau1_->clear();
    jetMassSub0_->clear();
    jetMassSub1_->clear();

    partonPt_->clear();
    partonY_->clear();
    partonMass_->clear();
    partonPhi_->clear();
    partonEta_->clear();

    jetBtagSub0DCSVbb_->clear();
    jetBtagSub1DCSVbb_->clear();
    jetBtagSub0DCSVbbb_->clear();
    jetBtagSub1DCSVbbb_->clear();

  xRecoAll.clear();
  bool partonCuts, recoCuts, massCut, tTaggerCut, triggerCR, triggerSR;
  bool deepCSV, revertBtagDeepCSV;
  bool btagCut, revertBtag;

  if (nJets >1)
  {
    //matching only if we have Signal ttbar MC
    if(selection == 1 || selection ==4)
    {
    //----------------------MATCHING------------------------------------------------------

    for(int ijet =0; ijet<nJets; ijet++)
    {
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
        isMatched++;
        jetDr_ = dRmin;
        //RECO MATCHED
        pt_->push_back((*jetPt)[(*partonMatchIdx)[indexMin]]);
        mass_->push_back((*jetMassSoftDrop)[(*partonMatchIdx)[indexMin]]);
        eta_->push_back((*jetEta)[(*partonMatchIdx)[indexMin]]);
        y_->push_back((*jetY)[(*partonMatchIdx)[indexMin]]);
        phi_->push_back( (*jetPhi)[(*partonMatchIdx)[indexMin]]);
        jetBtagSub0_->push_back( (*jetBtagSub0)[(*partonMatchIdx)[indexMin]]);
   //     jetBtagSub1_->push_back( (*jetBtagSub1)[(*partonMatchIdx)[indexMin]]);
        jetTtag_->push_back((*jetTtag)[(*partonMatchIdx)[indexMin]]);
        jetTau1_->push_back((*jetTau1)[(*partonMatchIdx)[indexMin]]);
        jetTau2_->push_back((*jetTau2)[(*partonMatchIdx)[indexMin]]);
        jetTau3_->push_back((*jetTau3)[(*partonMatchIdx)[indexMin]]);

        ecfB1N2_->push_back((*ecfB1N2)[(*partonMatchIdx)[indexMin]]);
        ecfB1N3_->push_back((*ecfB1N3)[(*partonMatchIdx)[indexMin]]);
        ecfB2N2_->push_back((*ecfB2N2)[(*partonMatchIdx)[indexMin]]);
        ecfB2N3_->push_back((*ecfB2N3)[(*partonMatchIdx)[indexMin]]);

        jetMassSub0_->push_back((*jetMassSub0)[(*partonMatchIdx)[indexMin]]);
        jetMassSub1_->push_back((*jetMassSub1)[(*partonMatchIdx)[indexMin]]);

        jetBtagSub0DCSVbb_->push_back((*jetBtagSub0DCSVbb)[(*partonMatchIdx)[indexMin]]);
        jetBtagSub1DCSVbb_->push_back((*jetBtagSub1DCSVbb)[(*partonMatchIdx)[indexMin]]);
        jetBtagSub0DCSVbbb_->push_back((*jetBtagSub0DCSVbbb)[(*partonMatchIdx)[indexMin]]);
        jetBtagSub1DCSVbbb_->push_back((*jetBtagSub1DCSVbbb)[(*partonMatchIdx)[indexMin]]);

        //PARTON MATCHED
        partonPt_->push_back((*partonPt)[indexMin]);
        partonMass_->push_back((*partonMass)[indexMin]);
        partonPhi_->push_back((*partonPhi)[indexMin]);
        partonEta_->push_back((*partonEta)[indexMin]);
        //here misssing partonY
      }
     }//----end of if jetMatchedIndexes > 0
    }//----end of for on all jets for mathching

    //---------------------------------------END OF MATCHING------------------------------------------------------
    float dCSVScoreSub0[2], dCSVScoreSub1[2];
    dCSVScoreSub0[0] = (*jetBtagSub0DCSVbb_)[0] + (*jetBtagSub0DCSVbbb_)[0];
    dCSVScoreSub0[1] = (*jetBtagSub0DCSVbb_)[1] + (*jetBtagSub0DCSVbbb_)[1];
    dCSVScoreSub1[0] = (*jetBtagSub1DCSVbb_)[0] + (*jetBtagSub1DCSVbbb_)[0];
    dCSVScoreSub1[1] = (*jetBtagSub1DCSVbb_)[1] + (*jetBtagSub1DCSVbbb_)[1];
   // cout<<"right before the cut init"<<endl;
    recoCuts   = fabs((*eta_)[0]) < 2.4 && fabs((*eta_)[1]) <2.4 && (*pt_)[0] > 400 && (*pt_)[1] > 400 && nLeptons==0 && mJJ > 1000 && nJets > 1;
    triggerSR  = (*bit)[triggerSRConst[year.Data()]];
    triggerCR  = (*bit)[triggerCRConst[year.Data()]];
    partonCuts = fabs((*partonEta_)[0]) < 2.4 && fabs((*partonEta_)[1]) <2.4 && (*partonPt_)[0] > 400 && (*partonPt_)[1] > 400 && mTTbarParton > 1000;
    massCut    = (*mass_)[0] > 120 && (*mass_)[0] < 220 && (*mass_)[1] > 120 && (*mass_)[1] < 220;
    tTaggerCut = true;
    //2 btag category with deepCSV
    deepCSV    = (((*jetBtagSub0DCSVbb_)[0] + (*jetBtagSub0DCSVbbb_)[0])> deepCSVFloat || ((*jetBtagSub1DCSVbb_)[0] + (*jetBtagSub1DCSVbbb_)[0])> deepCSVFloat) &&
           (((*jetBtagSub0DCSVbb_)[1] + (*jetBtagSub0DCSVbbb_)[1])> deepCSVFloat || ((*jetBtagSub1DCSVbb_)[1] + (*jetBtagSub1DCSVbbb_)[1])> deepCSVFloat);
   //0 btag category with deepCSV
    revertBtagDeepCSV = (dCSVScoreSub0[0] < deepCSVFloat &&  dCSVScoreSub1[0] < deepCSVFloat) && (dCSVScoreSub0[1] < deepCSVFloat && dCSVScoreSub1[1] < deepCSVFloat);
    if(isMatched > 1)
    {
    //cout<<"ok"<<endl;
    xRecoAll.push_back(bTagEvntWeight);
    }//if is matched
  }//----end of selection ==1 so that we do this only when we deal with signal MC

  else //we are in QCD samples or Subdominant BKG or Data sample
  {
    float dCSVScoreSub0[2], dCSVScoreSub1[2];
    dCSVScoreSub0[0] = (*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0];
    dCSVScoreSub0[1] = (*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1];
    dCSVScoreSub1[0] = (*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0];
    dCSVScoreSub1[1] = (*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1];

    recoCuts   = fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1]) <2.4 && (*jetPt)[0] > 400 && (*jetPt)[1] > 400 &&  mJJ > 1000 && nLeptons==0;
    triggerSR  = (*bit)[triggerSRConst[year.Data()]];
    triggerCR  = (*bit)[triggerCRConst[year.Data()]];
    massCut    = (*jetMassSoftDrop)[0] > 120 && (*jetMassSoftDrop)[0] < 220 && (*jetMassSoftDrop)[1] > 120 && (*jetMassSoftDrop)[1] < 220;
    //tTaggerCut = (*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut;
    tTaggerCut = true;
    //2 btag category with csvv2 and deepCSV
    deepCSV    = (((*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0])> deepCSVFloat) &&
           (((*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1])> deepCSVFloat);

    //0 btag category with deepCSV
    revertBtagDeepCSV = (dCSVScoreSub0[0] < deepCSVFloat &&  dCSVScoreSub1[0] < deepCSVFloat) && (dCSVScoreSub0[1] < deepCSVFloat && dCSVScoreSub1[1] < deepCSVFloat);

    xRecoAll.push_back(bTagEvntWeight);

  }//---end of else of isSignal

  btagCut = deepCSV;
  revertBtag = revertBtagDeepCSV;

   //cout<<"------"<<endl;
   for(int ivar = 0; ivar <xRecoAll.size(); ivar ++)
   {

     xReco = xRecoAll[ivar];
     if(selection == 0){
       genEvtWeight =1;
       bTagEvntWeight = 1;
     }

     //Signal Region with tTagger
     if(recoCuts && btagCut && massCut && tTaggerCut && triggerSR)
      hSR[f][ivar]->Fill(xReco,genEvtWeight);

     //Control Region with tTagger
     if(recoCuts && revertBtag && massCut && tTaggerCut && triggerCR)
      hCR[f][ivar]->Fill(xReco,genEvtWeight);

   }

  }//----end of nJets
  } //---end of event loop

  }//----end of fileSize loop

  TH1F *hCR_Clone[listOfFiles.size()][NVAR];
  TH1F *hSR_Clone[listOfFiles.size()][NVAR];

  TH1F *hCRTemp[NVAR], *hSRTemp[NVAR];
  TH1F *hCRTemp_Clone[NVAR], *hSRTemp_Clone[NVAR];

  for(int ivar= 0; ivar<NVAR; ivar++)
  {
    TString varNameReco = varReco[ivar];
    //for every slice
    for(int j=0; j<listOfFiles.size(); j++)
    {
      hCR_Clone[j][ivar]=(TH1F*)hCR[j][ivar]->Clone(TString::Format("hCR_%s_Clone",histoNames[j].Data()));
      hSR_Clone[j][ivar]=(TH1F*)hSR[j][ivar]->Clone(TString::Format("hSR_%s_Clone",histoNames[j].Data()));

      if(j==0)
      {
        hCRTemp[ivar]=(TH1F*)hCR[j][ivar]->Clone(TString::Format("hCRTemp_%s",histoNames[j].Data()));
        hSRTemp[ivar]=(TH1F*)hSR[j][ivar]->Clone(TString::Format("hSRTemp_%s",histoNames[j].Data()));

        hCRTemp_Clone[ivar]=(TH1F*)hCR[j][ivar]->Clone(TString::Format("hCRTemp_%s_Clone",histoNames[j].Data()));
        hSRTemp_Clone[ivar]=(TH1F*)hSR[j][ivar]->Clone(TString::Format("hSRTemp_%s_Clone",histoNames[j].Data()));
      }
        if(selection !=0)
        {
         hCR_Clone[j][ivar]->Scale(weights[j]*LUMI_CR); //this is 0 btagged (CR)
         hSR_Clone[j][ivar]->Scale(weights[j]*LUMI); //this is 2 btagged (SR)

         hCR[j][ivar]->Scale(weights[j]); //this is CR
         hSR[j][ivar]->Scale(weights[j]); //this is Signal region
        }
    }


    for(int j=1; j<listOfFiles.size(); j++)
    {
      cout<<"inside the loop!"<<endl;
      //Add them to get the whole phase space
      hCRTemp[ivar]->Add(hCR[j][ivar]);
      hSRTemp[ivar]->Add(hSR[j][ivar]);

      hCRTemp_Clone[ivar]->Add(hCR_Clone[j][ivar]);
      hSRTemp_Clone[ivar]->Add(hSR_Clone[j][ivar]);
    }


  }
  TFile *outFile;

  if(selection ==1)
    outFile = new TFile(TString::Format("%s/TopTaggerHisto_TT_Mtt-700toInf_100.root",year.Data()), "RECREATE");
  else if(selection ==2)
    outFile = new TFile(TString::Format("%s/TopTaggerHisto_QCD_HT300toInf_100.root",year.Data()), "RECREATE");
  else if(selection ==3)
    outFile = new TFile(TString::Format("%s/TopTaggerHisto_SubdominantBkgs_100.root",year.Data()), "RECREATE");
  else if(selection ==4)
    outFile = new TFile(TString::Format("%s/TopTaggerHisto_TT_NominalMC_100.root",year.Data()), "RECREATE");


  for(int ivar = 0; ivar<NVAR; ivar++)
  {
  TString varNameReco = varReco[ivar];
   outFile->cd();
   hSRTemp[ivar]->Write(TString::Format("hWt_%s_2btag", varNameReco.Data()));
   hCRTemp[ivar]->Write(TString::Format("hWt_%s_0btag", varNameReco.Data()));
   hSRTemp_Clone[ivar]->Write(TString::Format("hWt_%s_2btag_expYield", varNameReco.Data()));
   hCRTemp_Clone[ivar]->Write(TString::Format("hWt_%s_0btag_expYield", varNameReco.Data()));

   for(int j=1; j<listOfFiles.size(); j++)
   {
     hSR[j][ivar]->Write(TString::Format("hWt_%s_2btag_%s", varNameReco.Data(), histoNames[j].Data()));
     hCR[j][ivar]->Write(TString::Format("hWt_%s_0btag_%s", varNameReco.Data(), histoNames[j].Data()));
     hSR_Clone[j][ivar]->Write(TString::Format("hWt_%s_2btag_%s_expYield", varNameReco.Data(), histoNames[j].Data()));
     hCR_Clone[j][ivar]->Write(TString::Format("hWt_%s_0btag_%s_expYield", varNameReco.Data(), histoNames[j].Data()));
   }


 }
  //outFile->Close();
  listOfFiles.clear();
  XSEC.clear();
  histoNames.clear();

}
