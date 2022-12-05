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

TVector3 getBoostVector(TLorentzVector p4_1, TLorentzVector p4_2, TLorentzVector &p4CombinedVector);

void initFileNames(TString year)
{

  if(selection ==0) //data
  {
    eosPath = TString::Format("%s%s/",eosDataPath.Data(), year.Data());
    listOfFiles.push_back(dataFiles[year.Data()]);
  }
  else if(selection ==1) //signal ttbar mc mtt
  {
    eosPath = TString::Format("%s%s/Signal/",eosPathMC.Data(), year.Data());
    // cout<<eosPath<<endl;
    // cout<<mttFiles[year.Data()]["700-1000"]<<endl;
    //listOfFiles.push_back(mttFiles[year.Data()][""]);
    //listOfFiles.push_back(mttFiles[year.Data()][""]);
    //listOfFiles.push_back(mttFiles[year.Data()][""]);
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
    listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_tW_top_5f"]);
    listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_tW_antitop_5f"]);
    listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_t-channel_top_4f"]);
    if (!year.EqualTo("2017")) listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_t-channel_antitop_4f"]);
    listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_t-channel_top_5f"]);
    listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_t-channel_antitop_5f"]);
    listOfFiles.push_back(subdominantBkgFiles[year.Data()]["WJetsToQQ_HT-200to400"]);
    listOfFiles.push_back(subdominantBkgFiles[year.Data()]["WJetsToQQ_HT-400to600"]);
    listOfFiles.push_back(subdominantBkgFiles[year.Data()]["WJetsToQQ_HT-600to800"]);
    listOfFiles.push_back(subdominantBkgFiles[year.Data()]["WJetsToQQ_HT-800toInf"]);
  }
  else if(selection ==4) //signal ttbar mc nominal
  {
    cout<<"nominal!!!"<<endl;
    eosPath = TString::Format("%s%s/Signal/",eosPathMC.Data(), year.Data());
    cout<<eosPath<<endl;
    listOfFiles.push_back(ttNominalFiles[year.Data()]["TTHadronic"]);
    listOfFiles.push_back(ttNominalFiles[year.Data()]["TTSemiLeptonic"]);
    listOfFiles.push_back(ttNominalFiles[year.Data()]["TTTo2L2Nu"]);
  }
}

void initXsections(TString year)
{
  if(selection ==1) //signal ttbar mc
  {
    //XSEC.push_back(XSECAll[year.Data()]["TTToHadronic"]);
    //XSEC.push_back(XSECAll[year.Data()]["TTToSemiLeptonic"]);
    //XSEC.push_back(XSECAll[year.Data()]["TTTo2L2Nu"]);
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
    XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_tW_top_5f"]);
    XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_tW_antitop_5f"]);
    XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_t-channel_top_4f"]);
    if (!year.EqualTo("2017")) XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_t-channel_antitop_4f"]);
    XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_t-channel_top_5f"]);
    XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_t-channel_antitop_5f"]);
    XSEC.push_back(subdominantBkgXSEC[year.Data()]["WJetsToQQ_HT-200to400"]);
    XSEC.push_back(subdominantBkgXSEC[year.Data()]["WJetsToQQ_HT-400to600"]);
    XSEC.push_back(subdominantBkgXSEC[year.Data()]["WJetsToQQ_HT-600to800"]);
    XSEC.push_back(subdominantBkgXSEC[year.Data()]["WJetsToQQ_HT-800toInf"]);
  }
  else if(selection ==4)
  {
    XSEC.push_back(ttNominalXSEC[year.Data()]["TTHadronic"]);
    XSEC.push_back(ttNominalXSEC[year.Data()]["TTSemiLeptonic"]);
    XSEC.push_back(ttNominalXSEC[year.Data()]["TTTo2L2Nu"]);

  }

}

void initHistoNames(TString year)
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
    histoNames.push_back("ST_tW_top_5f_inclusiveDecays");
    histoNames.push_back("ST_tW_antitop_5f_inclusiveDecays");
    histoNames.push_back("ST_t-channel_top_4f_inclusiveDecays");
    if (!year.EqualTo("2017")) histoNames.push_back("ST_t-channel_antitop_4f_inclusiveDecays");
    histoNames.push_back("ST_t-channel_top_5f_inclusiveDecays");
    histoNames.push_back("ST_t-channel_antitop_5f_inclusiveDecays");
    histoNames.push_back("WJetsToQQ_HT-200to400");
    histoNames.push_back("WJetsToQQ_HT-400to600");
    histoNames.push_back("WJetsToQQ_HT-600to800");
    histoNames.push_back("WJetsToQQ_HT-800toInf");
  }
  else if(selection ==4)
  {
    histoNames.push_back("TTHadronic");
    histoNames.push_back("TTSemiLeptonic");
    histoNames.push_back("TTTo2L2Nu");
  }

}

void initGlobals(TString year)
{
  initFileNames(year);
  initXsections(year);
  initHistoNames(year);
}

void TagAndProbe(TString y="2016_preVFP", int sel = 0)
{
  year =y;
  initFilesMapping(false);
  deepCSVFloat = deepCSVFloatMap[year.Data()];
  selection = sel;
  LUMI = luminosity[year.Data()];
  LUMI_CR = luminosityCR[year.Data()];
  initGlobals(year);
  gStyle->SetOptStat(0);
  const int NVAR =14;

  float selMvaCut=topTaggerCuts[year];

  cout<<"triggerSRConst[year.Data()]]: "<<triggerSRConst[year.Data()]<<endl;
  cout<<"triggerCRConst[year.Data()]]: "<<triggerCRConst[year.Data()]<<endl;
  cout<<"topTagger: "<<selMvaCut<<endl;
  cout<<"deepCSVFloat: "<<deepCSVFloat<<endl;

  TString varReco[NVAR]   = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1","chi","cosThjetEta0", "cosThjetEta1",\
                            "mTop_Leading", "mTop_Subleading", "topTagger_leading", "topTagger_Subleading"};


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

  TFile *inf;
  int fileSize = listOfFiles.size();
  vector<float> weights(0);

  int NBINS[BND.size()];
  for (int i = 0; i<BND.size()+6; i++)
  {
    if(i<10) NBINS[i] = BND[i].size()-1;
    else NBINS[i] = 100;
  }

 //initialize the required histograms
  TH1F *h_Numerator[listOfFiles.size()][NVAR];
  TH1F *h_Denominator[listOfFiles.size()][NVAR];
  int maxJetMass = 220;
  int minJetMass = 120;

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
  vector<float> *jetPt(0),*tau3(0),*tau2(0),*tau1(0);
  vector<float> *jetMassSub0(0), *jetMassSub1(0);
  vector<float> *jetMassSoftDrop(0);

  float mva(0);
  vector<float> *jetTtag(0);
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
  trIN->SetBranchAddress("jetTau3"        ,&tau3);
  trIN->SetBranchAddress("jetTau2"        ,&tau2);
  trIN->SetBranchAddress("jetTau1"        ,&tau1);
  trIN->SetBranchAddress("triggerBit"     ,&bit);
  trIN->SetBranchAddress("genEvtWeight"   ,&genEvtWeight);
  trIN->SetBranchAddress("bTagEvntWeight"  ,&bTagEvntWeight);
  trIN->SetBranchAddress("jetMassSub0"    ,&jetMassSub0);
  trIN->SetBranchAddress("jetMassSub1"    ,&jetMassSub1);
  trIN->SetBranchAddress("mJJ"        ,&mJJ);
  trIN->SetBranchAddress("yJJ"        ,&yJJ);
  trIN->SetBranchAddress("ptJJ"         ,&ptJJ);
  //trIN->SetBranchAddress("jetBtagSub0"    ,&jetBtagSub0);
  //trIN->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1);
  trIN->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
  trIN->SetBranchAddress("mva"          ,&mva);
  trIN->SetBranchAddress("category"       ,&category);
  trIN->SetBranchAddress("deepAK8Tagger"  ,&deepAK8);
  trIN->SetBranchAddress("jetTtagCategory",&jetTtag);

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
  float tTaggerTight = 0;
  float tTaggerOther = 0;
  float tightTopTaggerCut = 0.8;

  float xParton(0), xReco(0);
  std::vector<float> xRecoAll(0);
  //book the histograms
  //histograms for Signal/QCD in CR
  for(int ivar =0; ivar< NVAR; ivar++)
  {
    int sizeBins = NBINS[ivar];
    if ((ivar==10) || (ivar==11))
    {
      h_Numerator[f][ivar] = new TH1F(TString::Format("hSRBTightAndSR%s_%s", varReco[ivar].Data(), histoNames[f].Data()), TString::Format("hSRBTightAndSR%s_%s",varReco[ivar].Data(), histoNames[f].Data()), sizeBins, minJetMass,maxJetMass);
      h_Denominator[f][ivar] = new TH1F(TString::Format("hSRBTightAndProbe_%s_%s", varReco[ivar].Data(), histoNames[f].Data()), TString::Format("hSRBTightAndProbe_%s_%s",varReco[ivar].Data(), histoNames[f].Data()), sizeBins, minJetMass,maxJetMass);
    }
    else if ((ivar==12)  || (ivar==13))//top tagger [-1,1]
    {
      h_Numerator[f][ivar] = new TH1F(TString::Format("hSRBTightAndSR%s_%s", varReco[ivar].Data(), histoNames[f].Data()), TString::Format("hSRBTightAndSR%s_%s",varReco[ivar].Data(), histoNames[f].Data()), sizeBins, -1, 1);
      h_Denominator[f][ivar] = new TH1F(TString::Format("hSRBTightAndProbe_%s_%s", varReco[ivar].Data(), histoNames[f].Data()), TString::Format("hSRBTightAndProbe_%s_%s",varReco[ivar].Data(), histoNames[f].Data()), sizeBins, -1, 1);
    }
    else if ((ivar==14) || (ivar==15))//deep ak8 [0,1]
    {
      h_Numerator[f][ivar] = new TH1F(TString::Format("hSRBTightAndSR%s_%s", varReco[ivar].Data(), histoNames[f].Data()), TString::Format("hSRBTightAndSR%s_%s",varReco[ivar].Data(), histoNames[f].Data()), sizeBins, 0,1);
      h_Denominator[f][ivar] = new TH1F(TString::Format("hSRBTightAndProbe_%s_%s", varReco[ivar].Data(), histoNames[f].Data()), TString::Format("hSRBTightAndProbe_%s_%s",varReco[ivar].Data(), histoNames[f].Data()), sizeBins, 0,1);
    }
    else
    {
      int sizeBins = NBINS[ivar];
      float tempBND[NBINS[ivar]+1];
      std::copy(BND[ivar].begin(), BND[ivar].end(), tempBND);
      h_Numerator[f][ivar] = new TH1F(TString::Format("hSRBTightAndSR%s_%s", varReco[ivar].Data(), histoNames[f].Data()), TString::Format("hSRBTightAndSR%s_%s",varReco[ivar].Data(), histoNames[f].Data()), sizeBins, tempBND);
      h_Denominator[f][ivar] = new TH1F(TString::Format("hSRBTightAndProbe_%s_%s", varReco[ivar].Data(), histoNames[f].Data()), TString::Format("hSRBTightAndProbe_%s_%s",varReco[ivar].Data(), histoNames[f].Data()), sizeBins, tempBND);
    }
    h_Numerator[f][ivar]->Sumw2();
    h_Denominator[f][ivar]->Sumw2();

  }


  cout<<"Reading "<<NN<<" entries"<<endl;
  for(int iev=0;iev<NN;iev++)
  {
    double progress = 10.0*iev/(1.0*NN);
    int k = TMath::FloorNint(progress);
    if (k > decade)
      cout<<10*k<<" %"<<endl;
    decade = k;
    trIN->GetEntry(iev);

    xRecoAll.clear();
    bool partonCuts, recoCuts, massCut, tTaggerCut, triggerCR, triggerSR;
    bool deepCSV, btag1DeepCSV, revertBtagDeepCSV;
    bool btagCut, revertBtag, btag1;
    if (nJets >1)
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
  
    float dCSVScoreSub0[2], dCSVScoreSub1[2];
    dCSVScoreSub0[0] = (*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0];
    dCSVScoreSub0[1] = (*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1];
    dCSVScoreSub1[0] = (*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0];
    dCSVScoreSub1[1] = (*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1];

    recoCuts   = fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1]) <2.4 && (*jetPt)[0] > 450 && (*jetPt)[1] > 400 &&  mJJ > 1000 && nLeptons==0;
    triggerSR  = (*bit)[triggerSRConst[year.Data()]];
    triggerCR  = (*bit)[triggerCRConst[year.Data()]];
    massCut    = (*jetMassSoftDrop)[0] > 50 && (*jetMassSoftDrop)[0] < 300 && (*jetMassSoftDrop)[1] > 50 && (*jetMassSoftDrop)[1] < 300;
    tTaggerCut = (*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut;
    //2 btag category with csvv2 and deepCSV
    deepCSV    = (((*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0])> deepCSVFloat) &&
          (((*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1])> deepCSVFloat);
    //1 btag category with  deepCSV
    btag1DeepCSV  = ((dCSVScoreSub0[0] > deepCSVFloat || dCSVScoreSub1[0] > deepCSVFloat) && (dCSVScoreSub0[1] < deepCSVFloat && dCSVScoreSub1[1] < deepCSVFloat)) ||
            ((dCSVScoreSub0[0] < deepCSVFloat && dCSVScoreSub1[0] < deepCSVFloat) && (dCSVScoreSub0[1] > deepCSVFloat || dCSVScoreSub1[1] > deepCSVFloat));

    //0 btag category with deepCSV
      revertBtagDeepCSV = (dCSVScoreSub0[0] < deepCSVFloat &&  dCSVScoreSub1[0] < deepCSVFloat) && (dCSVScoreSub0[1] < deepCSVFloat && dCSVScoreSub1[1] < deepCSVFloat);

    TLorentzVector p4T[2], p4T_ZMF[2], p4TTbar;
    p4T[0].SetPtEtaPhiM((*jetPt)[0], (*jetEta)[0], (*jetPhi)[0], (*jetMassSoftDrop)[0]);
    p4T[1].SetPtEtaPhiM((*jetPt)[1], (*jetEta)[1], (*jetPhi)[1], (*jetMassSoftDrop)[1]);

    TVector3 ttbarBoostVector = getBoostVector(p4T[0], p4T[1], p4TTbar);

    p4T_ZMF[0].SetPtEtaPhiM(p4T[0].Pt(), p4T[0].Eta(), p4T[0].Phi(), p4T[0].M());
    p4T_ZMF[1].SetPtEtaPhiM(p4T[1].Pt(), p4T[1].Eta(), p4T[1].Phi(), p4T[1].M());
    p4T_ZMF[0].Boost(ttbarBoostVector);
    p4T_ZMF[1].Boost(ttbarBoostVector);

    float chi0(0), chi1(0);
     //chi0 = (1 + fabs(TMath::Cos(p4T_ZMF[0].Theta()))) / ( 1 - fabs(TMath::Cos(p4T_ZMF[0].Theta())));
     //chi1 = (1 + fabs(TMath::Cos(p4T_ZMF[1].Theta()))) / ( 1 - fabs(TMath::Cos(p4T_ZMF[1].Theta())));
    float yStarExp  = TMath::Exp(fabs(p4T_ZMF[0].Rapidity() - p4T_ZMF[1].Rapidity())); //this is chi = e^(y*) , y* = 1/2(y1-y0)

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

  btagCut = deepCSV;
  revertBtag = revertBtagDeepCSV;
  btag1 = btag1DeepCSV;

   //Signal Region with tTagger
    for(int ivar = 0; ivar <xRecoAll.size(); ivar ++)
    {
      //Probe Region 2btags
      if(recoCuts && btagCut && triggerSR)
      {
          if (selection == 0)
          {
            genEvtWeight=1;
            bTagEvntWeight=1;
          }
          if(tTaggerTight > tightTopTaggerCut)
          {
            double weights_temp = genEvtWeight * bTagEvntWeight;
            h_Denominator[f][ivar]->Fill(xRecoAll[ivar], genEvtWeight*bTagEvntWeight);
            //cout<<"h_Denominator "<<xRecoAll[ivar]<<endl;
              
          }
          if(tTaggerTight > tightTopTaggerCut && tTaggerOther > selMvaCut)
          {
            double weights_temp = genEvtWeight * bTagEvntWeight;
            h_Numerator[f][ivar]->Fill(xRecoAll[ivar], genEvtWeight*bTagEvntWeight);
            //cout<<"h_Numerator "<<xRecoAll[ivar]<<endl;
          } 
      }
    }

    }//----end of nJets
  } //---end of event loop

  cout<<"counter: "<<counter<<endl;
  }//----end of fileSize loop

  for(int ivar= 0; ivar<NVAR; ivar++)
  {
    //for every slice
    for(int j=0; j<listOfFiles.size(); j++)
    {
        if(selection !=0)
        {
          h_Numerator[j][ivar]->Scale(weights[j]*LUMI);  
          h_Denominator[j][ivar]->Scale(weights[j]*LUMI); 
        }
    }


    for(int j=1; j<listOfFiles.size(); j++)
    {
      cout<<"inside the loop!"<<endl;
      //Add them to get the whole phase space
      h_Denominator[0][ivar]->Add(h_Denominator[j][ivar]);
      h_Numerator[0][ivar]->Add(h_Numerator[j][ivar]);
    }


  }
  TFile *outFile;
  if(selection == 0)
    outFile = new TFile(TString::Format("../%s/TagAndProbeHisto_Data.root",year.Data()), "RECREATE");
  else if(selection == 2)
    outFile = new TFile(TString::Format("../%s/TagAndProbeHisto_QCD_HT300toInf.root",year.Data()), "RECREATE");
  else if(selection == 3)
    outFile = new TFile(TString::Format("../%s/TagAndProbeHisto_SubdominantBkgs.root",year.Data()), "RECREATE");


  for(int ivar = 0; ivar<NVAR; ivar++)
  {
    outFile->cd();
    h_Numerator[0][ivar]->Write(TString::Format("hSRBTightAndSR_%s_expYield", varReco[ivar].Data()));
    h_Denominator[0][ivar]->Write(TString::Format("hSRBTightAndProbe_%s_expYield", varReco[ivar].Data()));
  }
  //outFile->Close();

  listOfFiles.clear();
  XSEC.clear();
  histoNames.clear();

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
