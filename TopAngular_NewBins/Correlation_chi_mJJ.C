#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"

#include "TemplateConstants_FillHistograms.h"
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
    listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_tW_top_5f"]);
    listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_tW_antitop_5f"]);
    listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_t-channel_top_4f"]);
    if (!year.EqualTo("2017")) listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_t-channel_antitop_4f"]);
    listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_t-channel_top_5f"]);
    listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_t-channel_antitop_5f"]);
  }
  else if(selection ==4) //signal ttbar mc nominal
  {
    cout<<"nominal!!!"<<endl;
    eosPath = TString::Format("%s%s/Signal/",eosPathMC.Data(), year.Data());
    cout<<eosPath<<endl;
    listOfFiles.push_back(ttNominalFiles[year.Data()]["TTHadronic_0"]);
    listOfFiles.push_back(ttNominalFiles[year.Data()]["TTSemiLeptonic_0"]);
    listOfFiles.push_back(ttNominalFiles[year.Data()]["TTTo2L2Nu_0"]);
  }
}

void initXsections(TString year)
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
    XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_tW_top_5f"]);
    XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_tW_antitop_5f"]);
    XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_t-channel_top_4f"]);
    if (!year.EqualTo("2017")) XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_t-channel_antitop_4f"]);
    XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_t-channel_top_5f"]);
    XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_t-channel_antitop_5f"]);
  }
  else if(selection ==4)
  {
    XSEC.push_back(ttNominalXSEC[year.Data()]["TTHadronic_0"]);
    XSEC.push_back(ttNominalXSEC[year.Data()]["TTSemiLeptonic_0"]);
    XSEC.push_back(ttNominalXSEC[year.Data()]["TTTo2L2Nu_0"]);

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
  }
  else if(selection ==4)
  {
    histoNames.push_back("TTHadronic_0");
    histoNames.push_back("TTSemiLeptonic_0");
    histoNames.push_back("TTTo2L2Nu_0");
  }

}

void initGlobals(TString year)
{
  initFileNames(year);
  initXsections(year);
  initHistoNames(year);
}

void Correlation_chi_mJJ(TString y="2016", int sel = 0, int massWindow=1000)
{
  year =y;
  initFilesMapping();
  deepCSVFloat = deepCSVFloatMap[year.Data()];
  selection = sel;
  LUMI = luminosity[year.Data()];
  LUMI_CR = luminosityCR[year.Data()];
  initGlobals(year);
  gStyle->SetOptStat(0);

  const int NVAR = 4;
  const int chiSize =11;
  const int cosSize = 10;
  const int mJJSize = 7;

  float selMvaCut=topTaggerCuts[year];

  cout<<"triggerSRConst[year.Data()]]: "<<triggerSRConst[year.Data()]<<endl;
  cout<<"triggerCRConst[year.Data()]]: "<<triggerCRConst[year.Data()]<<endl;
  cout<<"topTagger: "<<selMvaCut<<endl;
  cout<<"deepCSVFloat: "<<deepCSVFloat<<endl;

  int NBINS[NVAR] = {chiSize, cosSize, cosSize, mJJSize};
  std::vector< std::vector <Float_t> > const BND = {{1,2,3,4,5,6,7,8,9,10,13,16}, //chi
                                                    {-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1}, //|cosTheta*| leading
                                                    {-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1}, //|cosTheta*| subleading
                                                    {1000, 1200, 1400, 1600, 1800, 2000, 2400, 5000}}; //mJJ
  TString varReco[NVAR]   = {"chi", "cosTheta_0", "cosTheta_1", "mJJ"};

  int fileSize = listOfFiles.size();
  TFile *inf;
  vector<float> weights(0);

 //number of checks: one for Top Tagger and one for DAK8
 const int nChecks = 2;
 //initialize the required histograms
 TH2F *hCor[listOfFiles.size()];

 for(int f=0; f<listOfFiles.size(); f++)
 {
  int counter(0);
  cout<<"Entering "<<eosPath+listOfFiles[f]<<endl;
  inf = TFile::Open(eosPath+listOfFiles[f]);
  TTree *trIN    = (TTree*)inf->Get("boosted/events");

  if(selection != 0)
  {
    float NORM = ((TH1F*)inf->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
    weights.push_back(XSEC[f]/NORM);
  }
  cout<<"file: "<<eosPath+listOfFiles[f]<<endl;
 // cout<<"weight: "<<weights[f]<<endl;
  cout<<"LUMI: "<<LUMI<<endl;
  cout<<"LUMI_CR: "<<LUMI_CR<<endl;
  int decade(0);
  int NN = trIN->GetEntries();

  int nJets,nLeptons;
  float genEvtWeight(0);
  double  bTagEvntWeight(0);
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
  trIN->SetBranchAddress("jetBtagSub0"    ,&jetBtagSub0);
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

  float xParton(0), xReco(0);
  std::vector<float> xRecoAll(0);
  //book the TH2 histogram 
  int sizeBins_chi = NBINS[0];
  int sizeBins_mJJ = NBINS[3];
  float tempBND_chi[NBINS[0]+1];
  float tempBND_mJJ[NBINS[3]+1];
  std::copy(BND[0].begin(), BND[0].end(), tempBND_chi);
  std::copy(BND[3].begin(), BND[3].end(), tempBND_mJJ);
  hCor[f] = new TH2F(TString::Format("hCor_%s",histoNames[f].Data()), TString::Format("hCor_%s",histoNames[f].Data()), 
                    sizeBins_mJJ, tempBND_mJJ, sizeBins_chi, tempBND_chi);

  //for matching
  std::vector<int> *jetMatchedIndexes = new std::vector<int>(0);
  std::vector<float> *jetMatchedDr = new std::vector<float>(0);
  std::vector<float> *eta_ = new std::vector<float>(0);
  std::vector<float> *y_ = new std::vector<float>(0);
  std::vector<float> *phi_ = new std::vector<float>(0);
  std::vector<float> *mass_ = new std::vector<float>(0);
  std::vector<float> *pt_ = new std::vector<float>(0);
  std::vector<float> *deepAK8_ = new std::vector<float>(0);
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

  cout<<"Reading "<<NN<<" entries"<<endl;
  for(int iev=0;iev<NN;iev++)
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
    jetBtagSub1_->clear();
    jetTtag_->clear();
    deepAK8_->clear();

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
  bool deepCSV, btag1DeepCSV, revertBtagDeepCSV;
  bool btagCut, revertBtag, btag1;

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
        deepAK8_->push_back((*deepAK8)[(*partonMatchIdx)[indexMin]]);

        jetBtagSub0DCSVbb_->push_back((*jetBtagSub0DCSVbb)[(*partonMatchIdx)[indexMin]]);
        jetBtagSub1DCSVbb_->push_back((*jetBtagSub1DCSVbb)[(*partonMatchIdx)[indexMin]]);
        jetBtagSub0DCSVbbb_->push_back((*jetBtagSub0DCSVbbb)[(*partonMatchIdx)[indexMin]]);
        jetBtagSub1DCSVbbb_->push_back((*jetBtagSub1DCSVbbb)[(*partonMatchIdx)[indexMin]]);

        //PARTON MATCHED
        partonPt_->push_back( (*partonPt)[indexMin]);
        partonMass_->push_back( (*partonMass)[indexMin]);
        partonPhi_->push_back( (*partonPhi)[indexMin]);
        partonEta_->push_back( (*partonEta)[indexMin]);
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

    recoCuts   = fabs((*eta_)[0]) < 2.4 && fabs((*eta_)[1]) <2.4 && (*pt_)[0] > 400 && (*pt_)[1] > 400 && nLeptons==0 && mJJ > massWindow && nJets > 1;
    triggerSR  = (*bit)[triggerSRConst[year.Data()]];
    triggerCR  = (*bit)[triggerCRConst[year.Data()]];
    partonCuts = fabs((*partonEta_)[0]) < 2.4 && fabs((*partonEta_)[1]) <2.4 && (*partonPt_)[0] > 400 && (*partonPt_)[1] > 400 && mTTbarParton > 1000;
    massCut    = (*mass_)[0] > 120 && (*mass_)[0] < 220 && (*mass_)[1] > 120 && (*mass_)[1] < 220;
    tTaggerCut = (*jetTtag_)[0] > selMvaCut && (*jetTtag_)[1] > selMvaCut;
    //2 btag category with deepCSV
    deepCSV    = (((*jetBtagSub0DCSVbb_)[0] + (*jetBtagSub0DCSVbbb_)[0])> deepCSVFloat || ((*jetBtagSub1DCSVbb_)[0] + (*jetBtagSub1DCSVbbb_)[0])> deepCSVFloat) &&
           (((*jetBtagSub0DCSVbb_)[1] + (*jetBtagSub0DCSVbbb_)[1])> deepCSVFloat || ((*jetBtagSub1DCSVbb_)[1] + (*jetBtagSub1DCSVbbb_)[1])> deepCSVFloat);
    //1 btag category with deepCSV
    btag1DeepCSV  = ((dCSVScoreSub0[0] > deepCSVFloat || dCSVScoreSub1[0] > deepCSVFloat) && (dCSVScoreSub0[1] < deepCSVFloat && dCSVScoreSub1[1] < deepCSVFloat)) ||
            ((dCSVScoreSub0[0] < deepCSVFloat && dCSVScoreSub1[0] < deepCSVFloat) && (dCSVScoreSub0[1] > deepCSVFloat || dCSVScoreSub1[1] > deepCSVFloat));

   //0 btag category with deepCSV
    revertBtagDeepCSV = (dCSVScoreSub0[0] < deepCSVFloat &&  dCSVScoreSub1[0] < deepCSVFloat) && (dCSVScoreSub0[1] < deepCSVFloat && dCSVScoreSub1[1] < deepCSVFloat);

    if(isMatched > 1)
    {
      int leadingPt = 0;
      int subleadingPt = 1;
      if ((*pt_)[0] < (*pt_)[1])
      {
        leadingPt = 1;
        subleadingPt = 0;
      }
    
     TLorentzVector p4T[2], p4T_ZMF[2], p4TTbar;
     p4T[leadingPt].SetPtEtaPhiM((*pt_)[leadingPt], (*eta_)[leadingPt], (*phi_)[leadingPt], (*mass_)[leadingPt]);
     p4T[subleadingPt].SetPtEtaPhiM((*pt_)[subleadingPt], (*eta_)[subleadingPt], (*phi_)[subleadingPt], (*mass_)[subleadingPt]);

     TVector3 ttbarBoostVector = getBoostVector(p4T[leadingPt], p4T[subleadingPt], p4TTbar);

     p4T_ZMF[0].SetPtEtaPhiM(p4T[leadingPt].Pt(), p4T[leadingPt].Eta(), p4T[leadingPt].Phi(), p4T[leadingPt].M());
     p4T_ZMF[1].SetPtEtaPhiM(p4T[subleadingPt].Pt(), p4T[subleadingPt].Eta(), p4T[subleadingPt].Phi(), p4T[subleadingPt].M());
     p4T_ZMF[0].Boost(ttbarBoostVector);
     p4T_ZMF[1].Boost(ttbarBoostVector);


     float chi0(0), chi1(0);
     //chi0 = (1 + fabs(TMath::Cos(p4T_ZMF[0].Theta()))) / ( 1 - fabs(TMath::Cos(p4T_ZMF[0].Theta())));
     //chi1 = (1 + fabs(TMath::Cos(p4T_ZMF[1].Theta()))) / ( 1 - fabs(TMath::Cos(p4T_ZMF[1].Theta())));
     float yStarExp  = TMath::Exp(fabs(p4T_ZMF[0].Rapidity() - p4T_ZMF[1].Rapidity())); //this is chi = e^(|y*|) , y* = 1/2(y1-y0)

     xRecoAll.push_back(yStarExp); //this is chi
     xRecoAll.push_back(mJJ); //this is dijet mass

    }
    else continue;



    }//----end of selection ==1 so that we do this only when we deal with signal MC

  else //we are in QCD samples or Subdominant BKG or Data sample
  {
    if (jetBtagSub0DCSVbb->size() < 2) continue; 
    float dCSVScoreSub0[2], dCSVScoreSub1[2];
    dCSVScoreSub0[0] = (*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0];
    dCSVScoreSub0[1] = (*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1];
    dCSVScoreSub1[0] = (*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0];
    dCSVScoreSub1[1] = (*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1];

    int leadingPt = 0;
    int subleadingPt = 1;
    if ((*jetPt)[0] < (*jetPt)[1])
    {
      leadingPt = 1;
      subleadingPt = 0;
    }

    recoCuts   = fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1]) <2.4 && (*jetPt)[0] > 400 && (*jetPt)[1] > 400 &&  mJJ > massWindow && nLeptons==0;
    triggerSR  = (*bit)[triggerSRConst[year.Data()]];
    triggerCR  = (*bit)[triggerCRConst[year.Data()]];
    massCut    = (*jetMassSoftDrop)[0] > 120 && (*jetMassSoftDrop)[0] < 220 && (*jetMassSoftDrop)[1] > 120 && (*jetMassSoftDrop)[1] < 220;
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
     p4T[leadingPt].SetPtEtaPhiM((*jetPt)[leadingPt], (*jetEta)[leadingPt], (*jetPhi)[leadingPt], (*jetMassSoftDrop)[leadingPt]);
     p4T[subleadingPt].SetPtEtaPhiM((*jetPt)[subleadingPt], (*jetEta)[subleadingPt], (*jetPhi)[subleadingPt], (*jetMassSoftDrop)[subleadingPt]);

     TVector3 ttbarBoostVector = getBoostVector(p4T[leadingPt], p4T[subleadingPt], p4TTbar);

     p4T_ZMF[0].SetPtEtaPhiM(p4T[leadingPt].Pt(), p4T[leadingPt].Eta(), p4T[leadingPt].Phi(), p4T[leadingPt].M());
     p4T_ZMF[1].SetPtEtaPhiM(p4T[subleadingPt].Pt(), p4T[subleadingPt].Eta(), p4T[subleadingPt].Phi(), p4T[subleadingPt].M());
     p4T_ZMF[0].Boost(ttbarBoostVector);
     p4T_ZMF[1].Boost(ttbarBoostVector);


     float chi0(0), chi1(0);
     //chi0 = (1 + fabs(TMath::Cos(p4T_ZMF[0].Theta()))) / ( 1 - fabs(TMath::Cos(p4T_ZMF[0].Theta())));
     //chi1 = (1 + fabs(TMath::Cos(p4T_ZMF[1].Theta()))) / ( 1 - fabs(TMath::Cos(p4T_ZMF[1].Theta())));
     float yStarExp  = TMath::Exp(fabs(p4T_ZMF[0].Rapidity() - p4T_ZMF[1].Rapidity())); //this is chi = e^(y*) , y* = 1/2(y1-y0)

     xRecoAll.push_back(yStarExp); //this is chi
     xRecoAll.push_back(mJJ); //this is dijet mass

  }//---end of else of isSignal

  btagCut = deepCSV;
  revertBtag = revertBtagDeepCSV;
  btag1 = btag1DeepCSV;


     if(selection == 0){
       genEvtWeight =1;
       bTagEvntWeight = 1;
     }

    //Signal Region with tTagger
    if(recoCuts && btagCut && massCut && tTaggerCut && triggerSR)
      hCor[f]->Fill(xRecoAll[1], xRecoAll[0], genEvtWeight*bTagEvntWeight);

  }//----end of nJets
  } //---end of event loop

  }//----end of fileSize loop

  TH2F *hCor_Clone[listOfFiles.size()];

  //for every slice
  for(int j=0; j<listOfFiles.size(); j++)
  {
    hCor_Clone[j]=(TH2F*)hCor[j]->Clone(TString::Format("hCor_%s_Clone",histoNames[j].Data()));

    if(selection !=0)
    {
        hCor_Clone[j]->Scale(weights[j]*LUMI); //this is 0 btagged (CR)
        hCor[j]->Scale(weights[j]); //this is Signal region
    }
  }


  for(int j=1; j<listOfFiles.size(); j++)
  {
    //Add them to get the whole phase space
    hCor[0]->Add(hCor[j]);
    hCor_Clone[0]->Add(hCor_Clone[j]);
  }

  TFile *outFile;
  if(selection ==0)
    outFile = new TFile(TString::Format("%s/Correlation/Histo_Data_%s_reduced_%d.root",year.Data(),year.Data(),massWindow), "RECREATE");
  else if(selection ==1)
    outFile = new TFile(TString::Format("%s/Correlation/Histo_TT_Mtt-700toInf_reduced_%d.root",year.Data(),massWindow), "RECREATE");
  else if(selection ==2)
    outFile = new TFile(TString::Format("%s/Correlation/Histo_QCD_HT300toInf_reduced_%d.root",year.Data(),massWindow), "RECREATE");
  else if(selection ==3)
    outFile = new TFile(TString::Format("%s/Correlation/Histo_SubdominantBkgs_reduced_%d.root",year.Data(),massWindow), "RECREATE");
  else if(selection ==4)
    outFile = new TFile(TString::Format("%s/Correlation/Histo_TT_NominalMC_reduced_%d.root",year.Data(), massWindow), "RECREATE");


  hCor[0]->GetXaxis()->SetTitle("mJJ (GeV)");
  hCor_Clone[0]->GetXaxis()->SetTitle("mJJ (GeV)");

  hCor[0]->GetYaxis()->SetTitle("chi (#chi)");
  hCor_Clone[0]->GetYaxis()->SetTitle("chi (#chi)");

  outFile->cd();
  hCor[0]->Write("hCor_mJJ_chi");
  hCor_Clone[0]->Write("hCor_mJJ_chi_expYield");


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
