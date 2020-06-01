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
TString eosPath;
bool isSignal;
float deepCSVFloat = 0.6321;

void initFileNames()
{
	eosPath = "/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/Bkg/";
	listOfFiles.push_back("QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
	listOfFiles.push_back("QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
	listOfFiles.push_back("QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
	listOfFiles.push_back("QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
	listOfFiles.push_back("QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
	listOfFiles.push_back("QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
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
	histoNames.push_back("QCD_300_500");
	histoNames.push_back("QCD_500_700");
	histoNames.push_back("QCD_700_1000");
	histoNames.push_back("QCD_1000_1500");
	histoNames.push_back("QCD_1500_2000");
	histoNames.push_back("QCD_2000_Inf");
}

void initGlobals()
{
  initFileNames();
  initXsections();
  initHistoNames();
}
 
void htBkgs()
{
  initGlobals();	
//TString TTbarFile = "/eos/cms/store/user/gbakas/ttbar/topTagger/April19/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_Copy.root"	
  gStyle->SetOptStat(0);
 
  
  int fileSize = listOfFiles.size();
  TH1F *h_Chi_all[fileSize],*h_Cos_all[fileSize], *hChiRevertBtag[fileSize], *hCosRevertBtag[fileSize]; 
  TFile *inf;
  vector<float> weights(0);
  
 //initialize the required histograms 
 TH1F *h_hT[listOfFiles.size()];

 
 for(int f=0; f<listOfFiles.size(); f++)
 {

  h_hT[f] = new TH1F(histoNames[f].Data(), histoNames[f].Data(), 50, 300, 3000);
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
  float ht(0); 

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
  trIN->SetBranchAddress("mJJ"   		  ,&mJJ);
  trIN->SetBranchAddress("yJJ"   		  ,&yJJ);
  trIN->SetBranchAddress("ptJJ"   		  ,&ptJJ);
  trIN->SetBranchAddress("jetBtagSub0"	  ,&jetBtagSub0);
  trIN->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1);
  trIN->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
  trIN->SetBranchAddress("mva"	  		  ,&mva);
  trIN->SetBranchAddress("category"	  	  ,&category);
  trIN->SetBranchAddress("deepAK8Tagger"  ,&deepAK8);
  trIN->SetBranchAddress("jetTtagCategory",&jetTtag);	
  trIN->SetBranchAddress("ht"			  ,&ht);
  //deepCSV
  trIN->SetBranchAddress("jetBtagSub0DCSVbb" ,&jetBtagSub0DCSVbb);
  trIN->SetBranchAddress("jetBtagSub1DCSVbb" ,&jetBtagSub1DCSVbb);
  trIN->SetBranchAddress("jetBtagSub0DCSVbbb",&jetBtagSub0DCSVbbb);
  trIN->SetBranchAddress("jetBtagSub1DCSVbbb",&jetBtagSub1DCSVbbb);


  
 
  cout<<"Reading "<<NN<<" entries"<<endl;
  for(int iev=0;iev<NN;iev++) 
  {
    double progress = 10.0*iev/(1.0*NN);
    int k = TMath::FloorNint(progress); 
    if (k > decade) 
      cout<<10*k<<" %"<<endl;
    decade = k;
    trIN->GetEntry(iev);
    h_hT[f]->Fill(ht, genEvtWeight);
	
	
  }	//---end of event loop
  }//----end of fileSize loop 
  
  for(int f=0; f<listOfFiles.size(); f++)
  {
  	h_hT[f]->Scale(weights[f]*LUMI);
  }


  TFile *outFile = new TFile("htBkgFile.root", "RECREATE");

  outFile->cd();
  for(int f=0; f<listOfFiles.size(); f++)
  {
  	h_hT[f]->Write();
  }
}

  
