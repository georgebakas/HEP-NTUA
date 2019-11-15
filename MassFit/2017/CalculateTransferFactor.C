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
float LUMI = 41530;
TString eosPath;
float deepCSVFloat = 0.4941;

void initFileNames()
{
  eosPath = "/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2017/Bkg/";
  listOfFiles.push_back("QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8.root");
  listOfFiles.push_back("QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8.root");
  listOfFiles.push_back("QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8.root");
  listOfFiles.push_back("QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8.root");
  listOfFiles.push_back("QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8.root");
  listOfFiles.push_back("QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8.root");
}

void initXsections()
{
  XSEC.push_back(3.67e+5);
  XSEC.push_back(2.94e+4); 
  XSEC.push_back(6.524e+03);
  XSEC.push_back(1.064e+03);
  XSEC.push_back(121.5);
  XSEC.push_back(2.542e+01);
}

void initHistoNames()
{
  histoNames.push_back("QCD_histo_Mtt_300_500");
  histoNames.push_back("QCD_histo_Mtt_500_700");
  histoNames.push_back("QCD_histo_Mtt_700_1000");
  histoNames.push_back("QCD_histo_Mtt_1000_1500");
  histoNames.push_back("QCD_histo_Mtt_1500_2000");
  histoNames.push_back("QCD_histo_Mtt_2000_Inf");
}

void initGlobals()
{
  initFileNames();
  initXsections();
  initHistoNames();
}
 
void CalculateTransferFactor()
{

  initGlobals();  
  gStyle->SetOptStat(0);
  const int N_PT = 10;
  float selMvaCut=0.2;
  float floatBTag = 0.8838;
  const float BND[N_PT+1] = {400,450,500,570,650,750,850,950,1100,1300,1500}; //jetPt0 
  int fileSize = listOfFiles.size();
  TFile *inf;
  vector<float> weights(0);
  
 //initialize the required histograms 
 TH1F *hCR[listOfFiles.size()];
 TH1F *hSR[listOfFiles.size()];
 TH1F *h1Btag[listOfFiles.size()];

 TH1F *hCRA[listOfFiles.size()];
 TH1F *hSRA[listOfFiles.size()];
 TH1F *h1BtagA[listOfFiles.size()];
 
 for(int f=0; f<listOfFiles.size(); f++)
 {
  int counter;
  cout<<"Entering "<<listOfFiles[f]<<endl;
  inf = TFile::Open(eosPath+listOfFiles[f]);   
  TTree *trIN    = (TTree*)inf->Get("boosted/events");  
  
  float NORM = ((TH1F*)inf->Get("eventCounter/GenEventWeight"))->GetSumOfWeights(); 
  weights.push_back(XSEC[f]/NORM);  
  
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
  trIN->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1);
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

  //book the histograms
  //histograms for Signal/QCD in CR 

  int sizeBins = N_PT; 
  hCR[f] = new TH1F(TString::Format("hCR_%s_%s_%s", "tTagger",histoNames[f].Data(),"jetPt0"), TString::Format("hCR_%s_%s_%s","tTagger",histoNames[f].Data(),"jetPt0"), sizeBins, BND);
  hSR[f] = new TH1F(TString::Format("hSR_%s_%s_%s", "tTagger",histoNames[f].Data(),"jetPt0"), TString::Format("hSR_%s_%s_%s","tTagger",histoNames[f].Data(),"jetPt0"), sizeBins, BND);
  h1Btag[f] = new TH1F(TString::Format("h%s_%s_%s", "1Btag_tTagger",histoNames[f].Data(),"jetPt0"), TString::Format("h%s_%s_%s","1Btag_tTagger",histoNames[f].Data(),"jetPt0"), sizeBins, BND);

  hCRA[f] = new TH1F(TString::Format("hCRA_%s_%s_%s", "tTagger",histoNames[f].Data(),"jetPt0"), TString::Format("hCRA_%s_%s_%s","tTagger",histoNames[f].Data(),"jetPt0"), sizeBins, BND);
  hSRA[f] = new TH1F(TString::Format("hSRA_%s_%s_%s", "tTagger",histoNames[f].Data(),"jetPt0"), TString::Format("hSRA_%s_%s_%s","tTagger",histoNames[f].Data(),"jetPt0"), sizeBins, BND);
  h1BtagA[f] = new TH1F(TString::Format("h%s_%s_%s", "1BtagA_tTagger",histoNames[f].Data(),"jetPt0"), TString::Format("h%s_%s_%s","1BtagA_tTagger",histoNames[f].Data(),"jetPt0"), sizeBins, BND);
  
  cout<<"Reading "<<NN<<" entries"<<endl;
  for(int iev=0;iev<NN;iev++) 
  {
    double progress = 10.0*iev/(1.0*NN);
    int k = TMath::FloorNint(progress); 
    if (k > decade) 
      cout<<10*k<<" %"<<endl;
    decade = k;
    trIN->GetEntry(iev);
 
  bool partonCuts, recoCuts, massCut, tTaggerCut, massCutTight;
  bool deepCSV, btag1DeepCSV, revertBtagDeepCSV;  
  bool btagCut, revertBtag, btag1;
  
  if (nJets >1)
  { 
  

    float dCSVScoreSub0[2], dCSVScoreSub1[2];
    dCSVScoreSub0[0] = (*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0];
    dCSVScoreSub0[1] = (*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1];
    dCSVScoreSub1[0] = (*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0];
    dCSVScoreSub1[1] = (*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1];
    
    recoCuts   = fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1]) <2.4 && (*jetPt)[0] > 400 && (*jetPt)[1] > 400 &&  mJJ > 1000 && (*bit)[2] && nLeptons==0;
    massCut    = (*jetMassSoftDrop)[0] > 50 && (*jetMassSoftDrop)[0] < 300 && (*jetMassSoftDrop)[1] > 50 && (*jetMassSoftDrop)[1] < 300;
    massCutTight = (*jetMassSoftDrop)[0] > 120 && (*jetMassSoftDrop)[0] < 220 && (*jetMassSoftDrop)[1] > 120 && (*jetMassSoftDrop)[1] < 220;
    tTaggerCut = (*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut;
    //2 btag category with csvv2 and deepCSV
    deepCSV    = (((*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0])> deepCSVFloat) && 
           (((*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1])> deepCSVFloat);
    //1 btag category with  deepCSV   
    btag1DeepCSV  = ((dCSVScoreSub0[0] > deepCSVFloat || dCSVScoreSub1[0] > deepCSVFloat) && (dCSVScoreSub0[1] < deepCSVFloat && dCSVScoreSub1[1] < deepCSVFloat)) ||
            ((dCSVScoreSub0[0] < deepCSVFloat && dCSVScoreSub1[0] < deepCSVFloat) && (dCSVScoreSub0[1] > deepCSVFloat || dCSVScoreSub1[1] > deepCSVFloat));
   
    //0 btag category with deepCSV
      revertBtagDeepCSV = (dCSVScoreSub0[0] < deepCSVFloat &&  dCSVScoreSub1[0] < deepCSVFloat) && (dCSVScoreSub0[1] < deepCSVFloat && dCSVScoreSub1[1] < deepCSVFloat);  
  
      
  btagCut = deepCSV;
  revertBtag = revertBtagDeepCSV;
  btag1 = btag1DeepCSV;
  float xReco;
   //Signal Region with tTagger
     //cout<<"enter loop"<<endl;
     xReco = (*jetPt)[0];
     //for the jetMassSoftDrop just keep it simple from 50 to 300 GeV

     if(recoCuts && btagCut && massCutTight && tTaggerCut)       
      hSR[f]->Fill(xReco,genEvtWeight);
        //Control Region with tTagger
       if(recoCuts && revertBtag && massCutTight && tTaggerCut)
        hCR[f]->Fill(xReco,genEvtWeight);
       //1 btag region with tTagger
       if(recoCuts && massCutTight && tTaggerCut && btag1)
        h1Btag[f]->Fill(xReco,genEvtWeight); 
   	

   	//Extended regions (SRA, CRA, 1BTAG A)
       if(recoCuts && btagCut && massCut && tTaggerCut)
        hSRA[f]->Fill(xReco,genEvtWeight);
        //Control Region with tTagger
       if(recoCuts && revertBtag && massCut && tTaggerCut)
        hCRA[f]->Fill(xReco,genEvtWeight);
       //1 btag region with tTagger
       if(recoCuts && massCut && tTaggerCut && btag1)
        h1BtagA[f]->Fill(xReco,genEvtWeight);  
  

  }//----end of nJets
  } //---end of event loop

  }//----end of fileSize loop 
  


    //for every slice
    for(int j=0; j<listOfFiles.size(); j++)
    { 
         hCR[j]->Scale(weights[j]); //this is CR
         hSR[j]->Scale(weights[j]); //this is Signal region
         h1Btag[j]->Scale(weights[j]); //this is 1 btag

         hCRA[j]->Scale(weights[j]); //this is CRA
         hSRA[j]->Scale(weights[j]); //this is Signal region A
         h1BtagA[j]->Scale(weights[j]); //this is 1 btag A
    }    
  
    for(int j=1; j<listOfFiles.size(); j++)
    {
      //Add them to get the whole phase space
      hCR[0]->Add(hCR[j]);
      hSR[0]->Add(hSR[j]);
      h1Btag[0]->Add(h1Btag[j]); 

      hCRA[0]->Add(hCRA[j]);
      hSRA[0]->Add(hSRA[j]);
      h1BtagA[0]->Add(h1BtagA[j]);    
    } 
  
    
  
  TFile *outFile =  new TFile("TransferFactor.root", "RECREATE");

  hCR[0]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", "jetPt0"));
  hSR[0]->GetXaxis()->SetTitle(TString::Format("%s (GeV)","jetPt0"));
  h1Btag[0]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", "jetPt0"));
  
  hCRA[0]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", "jetPt0"));
  hSRA[0]->GetXaxis()->SetTitle(TString::Format("%s (GeV)","jetPt0"));
  h1BtagA[0]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", "jetPt0"));
 
  
  outFile->cd();
  hSR[0]->Write(TString::Format("hWt_%s_2btag", "jetPt0"));
  hCR[0]->Write(TString::Format("hWt_%s_0btag", "jetPt0"));
  h1Btag[0]->Write(TString::Format("hWt_%s_1btag", "jetPt0"));

  hSRA[0]->Write(TString::Format("hWt_%s_2btagA", "jetPt0"));
  hCRA[0]->Write(TString::Format("hWt_%s_0btagA", "jetPt0"));
  h1BtagA[0]->Write(TString::Format("hWt_%s_1btagA", "jetPt0"));


  float tFactor[3];

  tFactor[0] = (hCR[0]->GetEntries() / hCRA[0]->GetEntries());
  tFactor[1] = (h1Btag[0]->GetEntries() / h1BtagA[0]->GetEntries());
  tFactor[2] = (hSR[0]->GetEntries() / hSRA[0]->GetEntries());

  float x[3] = {0,1,2};
  //GetXaxis()->SetBinLabel(i+1,histoNames[i].Data());
  std::vector<TString> names = {"0btag", "1btag", "2btag"};
  TH1F *hf = new TH1F("hTransf", "hTransf",3,0,3);
  for(int i =0; i<3; i++ )
  {
    hf->SetBinContent(i+1, tFactor[i]);
  	hf->GetXaxis()->SetBinLabel(i+1,names[i].Data());
  }
  hf->Write("TransferFactor_hist");
  //outFile->Close();
  listOfFiles.clear();
  XSEC.clear();
  histoNames.clear();
  
}
