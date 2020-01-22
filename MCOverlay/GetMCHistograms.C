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
TString eosPath;
TString year;
int selection;

void initFileNames()
{
  if(selection ==1) //signal ttbar mc
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

}

void initHistoNames()
{
  
  if(selection ==1)
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

}

void initGlobals()
{
  initFileNames();
  initXsections();
  initHistoNames();
}
 
void GetMCHistograms(TString y="2016", int sel = 0)
{
  year =y;
  initFilesMapping();
  selection = sel;
  initGlobals();  
  gStyle->SetOptStat(0);
  const int NVAR =11;
  TString varReco[NVAR]   = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1",
               "mva", "topTagger1", "mTop", "jetMassSoftDrop"};  
  
  int fileSize = listOfFiles.size();
  TFile *inf;
  vector<float> weights(0);
  
 
 //initialize the required histograms 
 TH1F *h[listOfFiles.size()][NVAR];
 for(int f=0; f<listOfFiles.size(); f++)
 {
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
  trIN->SetBranchAddress("jetPt"          ,&jetPt);
  trIN->SetBranchAddress("jetEta"         ,&jetEta);
  trIN->SetBranchAddress("jetY"           ,&jetY);
  trIN->SetBranchAddress("jetPhi"         ,&jetPhi);
  trIN->SetBranchAddress("genEvtWeight"   ,&genEvtWeight);
  trIN->SetBranchAddress("mJJ"        ,&mJJ);
  trIN->SetBranchAddress("yJJ"        ,&yJJ);
  trIN->SetBranchAddress("ptJJ"         ,&ptJJ);
  trIN->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
  trIN->SetBranchAddress("mva"          ,&mva);

  
  if(selection == 1)
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

  float xReco(0);
  std::vector<float> xRecoAll(0);
  
  int NBINS = 100;
  float x_min[NVAR] = {1000,0,-2.4, 400, 400, 0.0, 0.0};
  float x_max[NVAR] = {5000, 1300, 2.4, 1500, 1500, 2.4, 2.4};
  //book the histograms
  //histograms for Signal/QCD in CR 
  for(int ivar =0; ivar< NVAR; ivar++)
  {
    if(ivar < 7)
      h[f][ivar] = new TH1F(TString::Format("h_%s_%s",histoNames[f].Data(),varReco[ivar].Data()), TString::Format("h_%s_%s",histoNames[f].Data(),varReco[ivar].Data()), NBINS, x_min[ivar], x_max[ivar]);
    else if (ivar == 9 || ivar == 10)
      h[f][ivar] = new TH1F(TString::Format("h_%s_%s",histoNames[f].Data(),varReco[ivar].Data()), TString::Format("h_%s_%s",histoNames[f].Data(),varReco[ivar].Data()), NBINS, 50,300);
    else if (ivar == 7 || ivar == 8)
      h[f][ivar] = new TH1F(TString::Format("h_%s_%s",histoNames[f].Data(),varReco[ivar].Data()), TString::Format("h_%s_%s_",histoNames[f].Data(),varReco[ivar].Data()), NBINS, -1,1);

  }  


  //event loop
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
  xRecoAll.push_back(mJJ);
  xRecoAll.push_back(ptJJ);
  xRecoAll.push_back(yJJ);
  xRecoAll.push_back((*jetPt)[0]);
  xRecoAll.push_back((*jetPt)[1]);

  xRecoAll.push_back(fabs((*jetY)[0]));
  xRecoAll.push_back(fabs((*jetY)[1]));
  xRecoAll.push_back((*jetTtag)[0]);
  xRecoAll.push_back((*jetTtag)[1]);
  xRecoAll.push_back((*jetMassSoftDrop)[0]);
  for(int ijet=1; ijet<nJets; ijet++)
   xRecoAll.push_back((*jetMassSoftDrop)[ijet]);

   for(int ivar = 0; ivar <xRecoAll.size(); ivar ++)
   {
     //cout<<"enter loop"<<endl;
     xReco = xRecoAll[ivar];
     //for the jetMassSoftDrop just keep it simple from 50 to 300 GeV
     if(ivar < 10)
       h[f][ivar]->Fill(xReco, genEvtWeight);
     else
       h[f][10]->Fill(xReco,genEvtWeight);
        
   }


  } //---end of event loop

  }//----end of fileSize loop 
  
  //this will be used for the combined scale to XSEC histogram
  TH1F *h_Clone[listOfFiles.size()][NVAR];

  for(int ivar= 0; ivar<NVAR; ivar++)
  {
    //for every slice
    for(int j=0; j<listOfFiles.size(); j++)
    {
        h_Clone[j][ivar]=(TH1F*)h[j][ivar]->Clone(TString::Format("h_%s_Clone",histoNames[j].Data()));       
        h_Clone[j][ivar]->Scale(weights[j]); //this is the new to be added afterwards
    }
    
    for(int j=1; j<listOfFiles.size(); j++)
    {      
      //Add them to get the whole phase space
      h_Clone[0][ivar]->Add(h_Clone[j][ivar]);
    } 
  
    
  }//end of ivar loop

  TFile *outFile;
  if(selection ==1)
    outFile = new TFile(TString::Format("HistoMC_TT_Mtt-700toInf_100bins_%s.root",year.Data()), "RECREATE");
  else if(selection ==2)
    outFile = new TFile(TString::Format("HistoMC_QCD_HT300toInf_100bins_%s.root",year.Data()), "RECREATE");



  for(int ivar = 0; ivar<NVAR; ivar++)
  {

    TString varNameReco = varReco[ivar];
    cout<<varNameReco<<endl;
    for(int f=0; f<listOfFiles.size(); f++)
    {
      if(ivar ==0 || ivar ==1 || ivar == 3 || ivar == 4 || ivar == 9 || ivar ==10 )
      {
        h[f][ivar]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameReco.Data()));
        h_Clone[0][ivar]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameReco.Data()));
      }
      else
      {
        h[f][ivar]->GetXaxis()->SetTitle(TString::Format("%s", varNameReco.Data()));
        h_Clone[0][ivar]->GetXaxis()->SetTitle(TString::Format("%s", varNameReco.Data()));
      }
    }
    
    if(ivar == 8)
    {
      for(int f=0; f<listOfFiles.size(); f++)
      {
        h[f][ivar-1]->Add(h[f][ivar]);
      }
      h_Clone[0][ivar-1]->Add(h_Clone[0][ivar]);
    }

  }//end of ivar loop
  

  //write histograms to file:
  for(int ivar=0; ivar<NVAR; ivar++)
  {
    for(int f=0; f<listOfFiles.size(); f++)
    {
      h[f][ivar]->Write(TString::Format("h_%s_%s",histoNames[f].Data(),varReco[ivar].Data()));
    }
    h_Clone[0][ivar]->Write(TString::Format("hScaledXSEC_%s", varReco[ivar].Data()));
  }
  listOfFiles.clear();
  XSEC.clear();
  histoNames.clear();
  
}
