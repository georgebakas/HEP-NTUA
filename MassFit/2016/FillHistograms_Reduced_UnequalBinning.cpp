#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"

using std::cin;
using std::cout;
using std::endl;


std::vector<TString> listOfFiles;
std::vector<float> XSEC;
std::vector<TString> histoNames;
float LUMI = 35922;
float LUMI_CR = 1670;
TString eosPath;
int selection;
float deepCSVFloat = 0.6321; //medium WP
//float deepCSVFloat = 0.2217; //Loose WP

void initFileNames()
{
  
  if(selection ==0) //data
  {
  eosPath = "/eos/cms/store/user/gbakas/ttbar/JetHT/2016/";  
  listOfFiles.push_back("JetHT_Run2016-17Jul2018.root");
  }
  else if(selection ==1) //signal ttbar mc
  {
  eosPath = "/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/Signal/";  
    listOfFiles.push_back("TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");
    listOfFiles.push_back("TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");
  }
  else if(selection ==2) //bkg mc
  {
  eosPath = "/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/Bkg/";
  listOfFiles.push_back("QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  listOfFiles.push_back("QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  listOfFiles.push_back("QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  listOfFiles.push_back("QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  listOfFiles.push_back("QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  listOfFiles.push_back("QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  }
  else if(selection ==3) //subdominant bkgs
  {
    eosPath = "/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/Bkg/";
  listOfFiles.push_back("DYJetsToQQ_HT180_13TeV-madgraphMLM-pythia8.root");
  listOfFiles.push_back("WJetsToQQ_HT180_13TeV-madgraphMLM-pythia8.root");
  listOfFiles.push_back("ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
  listOfFiles.push_back("ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
  listOfFiles.push_back("ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root");
  listOfFiles.push_back("ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root");
  }
  else if(selection ==4)
  {
    eosPath = "/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/Signal/";
    listOfFiles.push_back("TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");
  }
}

void initXsections()
{
  if(selection ==1)
  {
  XSEC.push_back(69.64);
  XSEC.push_back(16.74);
  }
  else if(selection ==2)
  {
  XSEC.push_back(3.67e+5);
  XSEC.push_back(2.94e+4);
  XSEC.push_back(6.524e+03);
  XSEC.push_back(1.064e+03);
  XSEC.push_back(121.5);
  XSEC.push_back(2.542e+01);
  }
  else if(selection ==3)
  {
  XSEC.push_back(1208);
  XSEC.push_back(3105);
  XSEC.push_back(38.09);
  XSEC.push_back(38.06);
  XSEC.push_back(35.6);
  XSEC.push_back(35.6);
  }
  else if(selection ==4)
  {
  XSEC.push_back(832);  
  }
}

void initHistoNames()
{
  
  if(selection ==0) histoNames.push_back("Data_2016");
  else if(selection ==1)
  {
  histoNames.push_back("Signal_histo_Mtt_700_1000"); 
  histoNames.push_back("Signal_histo_Mtt_1000_Inf");
  }
  else if (selection ==2)
  {
  histoNames.push_back("QCD_histo_Mtt_300_500");
  histoNames.push_back("QCD_histo_Mtt_500_700");
  histoNames.push_back("QCD_histo_Mtt_700_1000");
  histoNames.push_back("QCD_histo_Mtt_1000_1500");
  histoNames.push_back("QCD_histo_Mtt_1500_2000");
  histoNames.push_back("QCD_histo_Mtt_2000_Inf");
  }
  else if(selection ==3)
  {
    histoNames.push_back("DYJetsToQQ_HT180");
    histoNames.push_back("WJetsToQQ_HT180");
    histoNames.push_back("ST_tW_top_5f_inclusiveDecays");
    histoNames.push_back("ST_tW_antitop_5f_inclusiveDecays");
    histoNames.push_back("ST_t-channel_top_4f_inclusiveDecays");
    histoNames.push_back("ST_t-channel_antitop_4f_inclusiveDecays");
  }
  else if(selection ==4)
  {
    histoNames.push_back("Signal_histo_NominalMC");
  }

}


void initGlobals()
{
  initFileNames();
  initXsections();
  initHistoNames();
}
 
void FillHistograms_Reduced_UnequalBinning(int sel = 0)
{
  selection = sel;
  initGlobals();  
  gStyle->SetOptStat(0);
  const int NVAR =11;
  const int N_MJJ = 20;
  const int N_PTJJ = 18;
  const int N_YJJ = 8;
  const int N_PT = 20;
  const int N_JETY = 24;
  const int N_JETMASS = 100;
  const int N_MVA = 100;

  float selMvaCut=0.2;
  float floatBTag = 0.8838;
  
  int NBINS[NVAR] = {N_MJJ, N_PTJJ, N_YJJ, N_PT, N_PT ,N_JETY, N_JETY,N_MVA, N_MVA ,N_JETMASS, N_JETMASS};
  std::vector< std::vector <Float_t> > const BND = {{1000, 1100,1200,1300, 1400,1500, 1600,1700, 1800,1900, 2000,2200, 2400,2600, 2800,3000, 3200,3600, 4000,4500, 5000}, //mjj 21
                                                   {0,30,60,105,150,225,300,375,450,525,600,675,750,850,950,1025,1100,1200,1300}, //ptjj 19
                                                   {-2.4,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.4}, //yjj
                                                   {400,425,450,475,500,535,570,610,650,700,750,800,850,900,950,1025,1100,1200,1300,1400,1500}, //jetPt0 21   
                                                   {400,425,450,475,500,535,570,610,650,700,750,800,850,900,950,1025,1100,1200,1300,1400,1500}, //jetPt1 21   
                                                   {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4}, //jetY0
                                                   {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4}}; //jetY0 25


  
  TString varReco[NVAR]   = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1",
               "mva", "topTagger1", "mTop", "jetMassSoftDrop"};  
  
  int fileSize = listOfFiles.size();
  TFile *inf;
  vector<float> weights(0);
  
 //number of checks: one for Top Tagger and one for DAK8
 const int nChecks = 2;
 //initialize the required histograms 
 TH1F *hCR[listOfFiles.size()][NVAR];
 TH1F *hSR[listOfFiles.size()][NVAR];
 TH1F *h1Btag[listOfFiles.size()][NVAR];
 
 for(int f=0; f<listOfFiles.size(); f++)
 {
  int counter;
  cout<<"Entering "<<listOfFiles[f]<<endl;
  inf = TFile::Open(eosPath+listOfFiles[f]);   
  TTree *trIN    = (TTree*)inf->Get("boosted/events");  
  
  if(selection != 0)
  {
    float NORM = ((TH1F*)inf->Get("eventCounter/GenEventWeight"))->GetSumOfWeights(); 
    weights.push_back(XSEC[f]/NORM);  
  }
  
  int decade(0);
  int NN = trIN->GetEntries();
  
  int nJets,nLeptons;
  float genEvtWeight;
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
  trIN->SetBranchAddress("jetMassSub0"    ,&jetMassSub0);
  trIN->SetBranchAddress("jetMassSub1"    ,&jetMassSub1);
  trIN->SetBranchAddress("mJJ"        ,&mJJ);
  trIN->SetBranchAddress("yJJ"        ,&yJJ);
  trIN->SetBranchAddress("ptJJ"         ,&ptJJ);
  trIN->SetBranchAddress("jetBtagSub0"    ,&jetBtagSub0);
 // trIN->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1);
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
  //book the histograms
  //histograms for Signal/QCD in CR 
  for(int ivar =0; ivar< NVAR; ivar++)
  {
    int sizeBins = NBINS[ivar];
    if(ivar < 7)
    {
      float tempBND[NBINS[ivar]+1];
      std::copy(BND[ivar].begin(), BND[ivar].end(), tempBND); 
      hCR[f][ivar] = new TH1F(TString::Format("hCR_%s_%s_%s", "tTagger",histoNames[f].Data(),varReco[ivar].Data()), TString::Format("hCR_%s_%s_%s","tTagger",histoNames[f].Data(),varReco[ivar].Data()), sizeBins, tempBND);
      hSR[f][ivar] = new TH1F(TString::Format("hSR_%s_%s_%s", "tTagger",histoNames[f].Data(),varReco[ivar].Data()), TString::Format("hSR_%s_%s_%s","tTagger",histoNames[f].Data(),varReco[ivar].Data()), sizeBins, tempBND);
      h1Btag[f][ivar] = new TH1F(TString::Format("h%s_%s_%s", "1Btag_tTagger",histoNames[f].Data(),varReco[ivar].Data()), TString::Format("h%s_%s_%s","1Btag_tTagger",histoNames[f].Data(),varReco[ivar].Data()), sizeBins, tempBND);
    }
    else if (ivar == 9 || ivar == 10)
    {
      hCR[f][ivar] = new TH1F(TString::Format("hCR_%s_%s_%s", "tTagger",histoNames[f].Data(),varReco[ivar].Data()), TString::Format("hCR_%s_%s_%s","tTagger",histoNames[f].Data(),varReco[ivar].Data()), sizeBins, 120,220);
      hSR[f][ivar] = new TH1F(TString::Format("hSR_%s_%s_%s", "tTagger",histoNames[f].Data(),varReco[ivar].Data()), TString::Format("hSR_%s_%s_%s","tTagger",histoNames[f].Data(),varReco[ivar].Data()), sizeBins, 120,220);
    h1Btag[f][ivar] = new TH1F(TString::Format("h%s_%s_%s", "1Btag_tTagger",histoNames[f].Data(),varReco[ivar].Data()), TString::Format("h%s_%s_%s","1Btag_tTagger",histoNames[f].Data(),varReco[ivar].Data()), sizeBins, 120,220);
    }
    else if (ivar == 7 || ivar == 8)
    {
      hCR[f][ivar] = new TH1F(TString::Format("hCR_%s_%s_%s", "tTagger",histoNames[f].Data(),varReco[ivar].Data()), TString::Format("hCR_%s_%s_%s","tTagger",histoNames[f].Data(),varReco[ivar].Data()), sizeBins, -1,1);
      hSR[f][ivar] = new TH1F(TString::Format("hSR_%s_%s_%s", "tTagger",histoNames[f].Data(),varReco[ivar].Data()), TString::Format("hSR_%s_%s_%s","tTagger",histoNames[f].Data(),varReco[ivar].Data()), sizeBins, -1,1);
    h1Btag[f][ivar] = new TH1F(TString::Format("h%s_%s_%s", "1Btag_tTagger",histoNames[f].Data(),varReco[ivar].Data()), TString::Format("h%s_%s_%s","1Btag_tTagger",histoNames[f].Data(),varReco[ivar].Data()), sizeBins, -1,1);
    }

  }
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
    
    recoCuts   = fabs((*eta_)[0]) < 2.4 && fabs((*eta_)[1]) <2.4 && (*pt_)[0] > 400 && (*pt_)[1] > 400 && nLeptons==0 && mJJ > 1000 && nJets > 1;
    triggerSR  = (*bit)[2];
    triggerCR  = (*bit)[4];
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
        
    xRecoAll.push_back(mJJ);
    xRecoAll.push_back(ptJJ);
    xRecoAll.push_back(yJJ);
    xRecoAll.push_back((*pt_)[leadingPt]);
    xRecoAll.push_back((*pt_)[subleadingPt]);
    xRecoAll.push_back(fabs((*y_)[leadingPt]));
    xRecoAll.push_back(fabs((*y_)[subleadingPt])); 
    xRecoAll.push_back((*jetTtag_)[0]);
    xRecoAll.push_back((*jetTtag_)[1]);
    xRecoAll.push_back((*mass_)[leadingPt]);    
    for(int ijet=1; ijet<nJets; ijet++)
      xRecoAll.push_back((*mass_)[ijet]);
    }
    else continue;


   
    }//----end of selection ==1 so that we do this only when we deal with signal MC 
  
  else //we are in QCD samples or Subdominant BKG or Data sample
  {
    float dCSVScoreSub0[2], dCSVScoreSub1[2];
    dCSVScoreSub0[0] = (*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0];
    dCSVScoreSub0[1] = (*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1];
    dCSVScoreSub1[0] = (*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0];
    dCSVScoreSub1[1] = (*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1];
    
    recoCuts   = fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1]) <2.4 && (*jetPt)[0] > 400 && (*jetPt)[1] > 400 &&  mJJ > 1000 && nLeptons==0;
    triggerSR  = (*bit)[2];
    triggerCR  = (*bit)[4];
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
  
   
   xRecoAll.push_back(mJJ);
   xRecoAll.push_back(ptJJ);
   xRecoAll.push_back(yJJ);
   xRecoAll.push_back((*jetPt)[0]);
   xRecoAll.push_back((*jetPt)[1]);
   if(selection ==0)
   {
    xRecoAll.push_back(1);
    xRecoAll.push_back(1);
   }
   else
   {
    xRecoAll.push_back(fabs((*jetY)[0]));
    xRecoAll.push_back(fabs((*jetY)[1]));
   }
   xRecoAll.push_back((*jetTtag)[0]);
   xRecoAll.push_back((*jetTtag)[1]);
   xRecoAll.push_back((*jetMassSoftDrop)[0]);
   for(int ijet=1; ijet<nJets; ijet++)
    xRecoAll.push_back((*jetMassSoftDrop)[ijet]);
  
  
  
    
  }//---end of else of isSignal
  
  btagCut = deepCSV;
  revertBtag = revertBtagDeepCSV;
  btag1 = btag1DeepCSV;

   //Signal Region with tTagger
   for(int ivar = 0; ivar <xRecoAll.size(); ivar ++)
   {
     //cout<<"enter loop"<<endl;
     xReco = xRecoAll[ivar];
     //genEventWeight is set probably to a value and this is why the histos have so many entries
     if(selection == 0) genEvtWeight =1;
     //for the jetMassSoftDrop just keep it simple from 50 to 300 GeV
     if(ivar < 10)
     {
       if(recoCuts && btagCut && massCut && tTaggerCut && triggerSR)
       {
        counter++;
        hSR[f][ivar]->Fill(xReco,genEvtWeight);
       }
        //Control Region with tTagger
       if(recoCuts && revertBtag && massCut && tTaggerCut && triggerCR)
        hCR[f][ivar]->Fill(xReco,genEvtWeight);
       //1 btag region with tTagger
       if(recoCuts && massCut && tTaggerCut && btag1 && triggerSR)
        h1Btag[f][ivar]->Fill(xReco,genEvtWeight);  
    }
    else
    {
      //Signal Region with tTagger
      if(recoCuts && btagCut && massCut && tTaggerCut && triggerSR)
      {
        hSR[f][10]->Fill(xReco,genEvtWeight);
      }
        
      //Control Region with tTagger
      if(recoCuts && revertBtag && massCut && tTaggerCut && triggerCR)
      {
        hCR[f][10]->Fill(xReco,genEvtWeight);  
      }
        
      //1 btag region with tTagger
      if(recoCuts && massCut && tTaggerCut && btag1 && triggerSR)
      {
        h1Btag[f][10]->Fill(xReco,genEvtWeight);  
      } 
    }
   }
   
  

  }//----end of nJets
  } //---end of event loop

  cout<<"counter: "<<counter<<endl;
  }//----end of fileSize loop 
  
  TH1F *hCR_Clone[listOfFiles.size()][NVAR];
  TH1F *hSR_Clone[listOfFiles.size()][NVAR];
  TH1F *h1Btag_Clone[listOfFiles.size()][NVAR];
  
  cout<<hSR[0][9]->GetEntries()<<endl;
  hSR[0][9]->Draw();

  for(int ivar= 0; ivar<NVAR; ivar++)
  {
    //for every slice
    for(int j=0; j<listOfFiles.size(); j++)
    {
      hCR_Clone[j][ivar]=(TH1F*)hCR[j][ivar]->Clone(TString::Format("hCR_%s_%s_Clone","tTagger",histoNames[j].Data())); 
      hSR_Clone[j][ivar]=(TH1F*)hSR[j][ivar]->Clone(TString::Format("hSR_%s_%s_Clone","tTagger",histoNames[j].Data()));
      h1Btag_Clone[j][ivar]=(TH1F*)h1Btag[j][ivar]->Clone(TString::Format("h1Btag_%s_%s_Clone","tTagger",histoNames[j].Data())); 
      
        if(selection !=0)
        {
         hCR_Clone[j][ivar]->Scale(weights[j]*LUMI_CR); //this is 0 btagged (CR)
         hSR_Clone[j][ivar]->Scale(weights[j]*LUMI); //this is 2 btagged (SR)
         h1Btag_Clone[j][ivar]->Scale(weights[j]*LUMI); //this is 1 btagged 
    
         hCR[j][ivar]->Scale(weights[j]); //this is CR
         hSR[j][ivar]->Scale(weights[j]); //this is Signal region
         h1Btag[j][ivar]->Scale(weights[j]); //this is 1 btag
        }
    }
    
  
    for(int j=1; j<listOfFiles.size(); j++)
    {
      cout<<"inside the loop!"<<endl;
      //Add them to get the whole phase space
      hCR[0][ivar]->Add(hCR[j][ivar]);
      hSR[0][ivar]->Add(hSR[j][ivar]);
      h1Btag[0][ivar]->Add(h1Btag[j][ivar]);    
      hCR_Clone[0][ivar]->Add(hCR_Clone[j][ivar]);
      hSR_Clone[0][ivar]->Add(hSR_Clone[j][ivar]);
      h1Btag_Clone[0][ivar]->Add(h1Btag_Clone[j][ivar]);
      
    } 
  
    
  }

  TFile *outFile;
  if(selection ==0)
    outFile = new TFile("Histo_Data_2016_100_reduced_UnequalBinning.root", "RECREATE");
  else if(selection ==1)
    outFile = new TFile("Histo_TT_Mtt-700toInf_100_reduced_UnequalBinning.root", "RECREATE");
  else if(selection ==2)
    outFile = new TFile("Histo_QCD_HT300toInf_100_reduced_UnequalBinning.root", "RECREATE");
  else if(selection ==3)
    outFile = new TFile("Histo_SubdominantBkgs_100_reduced_UnequalBinning.root", "RECREATE");
  else if(selection ==4)
    outFile = new TFile("Histo_TT_NominalMC_100_reduced_UnequalBinning.root", "RECREATE");

  for(int ivar = 0; ivar<NVAR; ivar++)
  {

  TString varNameReco = varReco[ivar];
  cout<<varNameReco<<endl;
    if(ivar ==0 || ivar ==1 || ivar == 3 || ivar == 4 || ivar == 9 || ivar ==10 )
  {
    hCR[0][ivar]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameReco.Data()));
    hSR[0][ivar]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameReco.Data()));
    h1Btag[0][ivar]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameReco.Data()));

    hCR_Clone[0][ivar]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameReco.Data()));
    hSR_Clone[0][ivar]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameReco.Data()));
    h1Btag_Clone[0][ivar]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameReco.Data()));
  }
  else
  {
    hCR[0][ivar]->GetXaxis()->SetTitle(TString::Format("%s", varNameReco.Data()));
    hSR[0][ivar]->GetXaxis()->SetTitle(TString::Format("%s", varNameReco.Data()));
    h1Btag[0][ivar]->GetXaxis()->SetTitle(TString::Format("%s", varNameReco.Data()));
    hCR_Clone[0][ivar]->GetXaxis()->SetTitle(TString::Format("%s", varNameReco.Data()));
    hSR_Clone[0][ivar]->GetXaxis()->SetTitle(TString::Format("%s", varNameReco.Data()));
    h1Btag_Clone[0][ivar]->GetXaxis()->SetTitle(TString::Format("%s", varNameReco.Data()));
  }
  
  
  
  
  
  if(ivar == 8)
  {
    hSR[0][ivar-1]->Add(hSR[0][ivar]);
    hCR[0][ivar-1]->Add(hCR[0][ivar]);
    h1Btag[0][ivar-1]->Add(h1Btag[0][ivar]);
    hSR_Clone[0][ivar-1]->Add(hSR_Clone[0][ivar]);
    hCR_Clone[0][ivar-1]->Add(hCR_Clone[0][ivar]);
    h1Btag_Clone[0][ivar-1]->Add(h1Btag_Clone[0][ivar]); 
  }
  
   outFile->cd();
   hSR[0][ivar]->Write(TString::Format("hWt_%s_2btag", varNameReco.Data()));
   hCR[0][ivar]->Write(TString::Format("hWt_%s_0btag", varNameReco.Data()));
   h1Btag[0][ivar]->Write(TString::Format("hWt_%s_1btag", varNameReco.Data()));
   hSR_Clone[0][ivar]->Write(TString::Format("hWt_%s_2btag_expYield", varNameReco.Data()));
   hCR_Clone[0][ivar]->Write(TString::Format("hWt_%s_0btag_expYield", varNameReco.Data()));
   h1Btag_Clone[0][ivar]->Write(TString::Format("hWt_%s_1btag_expYield", varNameReco.Data()));
  
  
 }
  //outFile->Close();
  
  listOfFiles.clear();
  XSEC.clear();
  histoNames.clear();
  
}
