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
float LUMI;
bool globalIsNominalMC;

void initFileNames()
{
  eosPath = TString::Format("%s%s/Signal/",eosPathMC.Data(), year.Data());  
  cout<<eosPath<<endl;
  if(!globalIsNominalMC)
  {
      listOfFiles.push_back(ttFiles[year.Data()]["700-1000"]);
      listOfFiles.push_back(ttFiles[year.Data()]["1000-Inf"]);
  }
  else
  {
      if(year.EqualTo("2016"))
          listOfFiles.push_back(ttFiles[year.Data()]["TTNominal"]);
      else
      {
          listOfFiles.push_back(ttFiles[year.Data()]["TTHadronic"]);
          listOfFiles.push_back(ttFiles[year.Data()]["TTSemiLeptonic"]);
          listOfFiles.push_back(ttFiles[year.Data()]["TTTo2L2Nu"]);
    }
  }
}

void initXsections()
{
    if(!globalIsNominalMC)
  {
       XSEC.push_back(ttXSEC[year.Data()]["700-1000"]);
       XSEC.push_back(ttXSEC[year.Data()]["1000-Inf"]);
  }
  else
  {
      if(year.EqualTo("2016"))
          XSEC.push_back(ttXSEC[year.Data()]["TTNominal"]);
      else
      {
          XSEC.push_back(ttXSEC[year.Data()]["TTHadronic"]);
          XSEC.push_back(ttXSEC[year.Data()]["TTSemiLeptonic"]);
          XSEC.push_back(ttXSEC[year.Data()]["TTTo2L2Nu"]);
    }
  }
}

void initHistoNames()
{
  
  if(!globalIsNominalMC)
  {
      histoNames.push_back("Signal_histo_Mtt_700_1000"); 
      histoNames.push_back("Signal_histo_Mtt_1000_Inf");
  }
  else
  {
      if(year.EqualTo("2016"))
          histoNames.push_back("Signal_histo_Nominal");
      else
      {
      histoNames.push_back("Signal_histo_TTHadronic");
      histoNames.push_back("Signal_histo_TTSemiLeptonic");
      histoNames.push_back("Signal_histo_TTTo2L2Nu");
    }
  }
}

void initGlobals()
{
  initFileNames();
  initXsections();
  initHistoNames();
}
 
void GetPartonMCHistograms(TString y="2017", bool isNominalMC= false)
{
  globalIsNominalMC = isNominalMC;  
  year =y;
  initFilesMapping();
  initGlobals();  
  gStyle->SetOptStat(0);
  const int NVAR =11;
  TString varReco[NVAR]   = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1","jetPhi0","jetPhi1","mTop0", "mTop1"};  
  TString varParton[NVAR] = {"mTTbarParton", "ptTTbarParton", "yTTbarParton","Leading PartonPt","Subleading PartonPt", 
                              "Leading jetY", "Subleading JetY","Leading jetPhi", "Subleading JetPhi","Leading mTop", "Subleading mTop"};  
  
  int fileSize = listOfFiles.size();
  TFile *inf;
  vector<float> weights(0);
  
 LUMI = luminosity[year];
 //initialize the required histograms 
 TH1F *h[listOfFiles.size()][NVAR], *hParton[listOfFiles.size()][NVAR], *hParticle[listOfFiles.size()][NVAR];
 for(int f=0; f<listOfFiles.size(); f++)
 {
  cout<<"Entering "<<eosPath<<listOfFiles[f]<<endl;
  cout<<"XSEC of slice: "<<XSEC[f]<<endl;
  cout<<"LUMI: "<<LUMI<<endl;
  inf = TFile::Open(eosPath+listOfFiles[f]);   
  TTree *trINCnt = (TTree*)inf->Get("eventCounter/events"); 
  
  float NORM = ((TH1F*)inf->Get("eventCounter/GenEventWeight"))->GetSumOfWeights(); 
  weights.push_back(XSEC[f]/NORM);  
  
  float genEvtWeight;
  float mTTbarParton(0),yTTbarParton(0), ptTTbarParton(0);
  
  float partonEta[2], partonPhi[2], partonY[2], partonPt[2], partonE[2], partonMass[2];
  //------- input tree --------------
  
  trINCnt->SetBranchAddress("genEvtWeight", &genEvtWeight);
  trINCnt->SetBranchAddress("mTTbarParton" ,&mTTbarParton);
  trINCnt->SetBranchAddress("yTTbarParton" ,&yTTbarParton);
  trINCnt->SetBranchAddress("ptTTbarParton"  ,&ptTTbarParton);
  trINCnt->SetBranchAddress("ptTopParton"    ,&partonPt);
  trINCnt->SetBranchAddress("etaTopParton"   ,&partonEta);
  trINCnt->SetBranchAddress("yTopParton"     ,&partonY);
  trINCnt->SetBranchAddress("mTopParton"     ,&partonMass);
  trINCnt->SetBranchAddress("phiTopParton"   ,&partonPhi);
  //trINCnt->SetBranchAddress("partonMatchDR"  ,&partonMatchDR);
  //trINCnt->SetBranchAddress("partonMatchIdx" ,&partonMatchIdx);
  float xParton(0);
  std::vector<float> xPartonAll(0);

  int NBINS = 100;
  float x_min[NVAR] = {0,0,-3.5, 0,0, 0.0, 0.0, -3.5,-3.5,0, 0};
  float x_max[NVAR] = {5000, 1300, 3.5, 1500, 1500, 3.5, 3.5, 3.5,3.5,300, 300};
  //book the histograms
  //histograms for Signal/QCD in CR 
  for(int ivar =0; ivar< NVAR; ivar++)
  {
    hParton[f][ivar] = new TH1F(TString::Format("hParton_%s_%s",histoNames[f].Data(),varReco[ivar].Data()), TString::Format("hParton_%s_%s",histoNames[f].Data(),varReco[ivar].Data()), NBINS, x_min[ivar], x_max[ivar]);
  }  

  //now for parton
  
    int NNCnt = trINCnt->GetEntries();
    for(int iev =0; iev<NNCnt; iev++)
    {
      xPartonAll.clear();
      trINCnt->GetEntry(iev);

        //cout<<mTTbarParton<<endl;
      xPartonAll.push_back(mTTbarParton);
      xPartonAll.push_back(ptTTbarParton);
      xPartonAll.push_back(yTTbarParton);
      xPartonAll.push_back(partonPt[0]);
      xPartonAll.push_back(partonPt[1]);
      xPartonAll.push_back(fabs(partonY[0]));
      xPartonAll.push_back(fabs(partonY[1]));
      xPartonAll.push_back(partonPhi[0]);
      xPartonAll.push_back(partonPhi[1]);
      xPartonAll.push_back(partonMass[0]);
      xPartonAll.push_back(partonMass[1]);

      for(int ivar = 0; ivar <xPartonAll.size(); ivar ++)
        {
         xParton = xPartonAll[ivar];
         hParton[f][ivar]->Fill(xParton, genEvtWeight);
      }
    }

  }//----end of fileSize loop 
  
  //this will be used for the combined scale to XSEC histogram
  TH1F *hParton_Clone[listOfFiles.size()][NVAR];
  for(int ivar= 0; ivar<NVAR; ivar++)
  {
    //for every slice
    for(int j=0; j<listOfFiles.size(); j++)
    {
      hParton_Clone[j][ivar]=(TH1F*)hParton[j][ivar]->Clone(TString::Format("hParton_%s_Clone",histoNames[j].Data()));       
      hParton_Clone[j][ivar]->Scale(weights[j]*LUMI); //this is the new to be added afterwards
      hParton[j][ivar]->Scale(weights[j]*LUMI);
    }
    for(int j=1; j<listOfFiles.size(); j++)
    {      
      hParton_Clone[0][ivar]->Add(hParton_Clone[j][ivar]);
    }  
  }//end of ivar loop

  TFile *outFile;
  if(!isNominalMC)
    outFile = new TFile(TString::Format("HistoMC_TT_Mtt-700toInf_100bins_%sParton.root",year.Data()), "RECREATE");
  else
    outFile = new TFile(TString::Format("HistoMC_TT_NominalMC_100bins_%sParton.root",year.Data()), "RECREATE");
  
  for(int ivar = 0; ivar<NVAR; ivar++)
  {
    TString varNameParton = varParton[ivar];
   cout<<varNameParton<<endl;
    for(int f=0; f<listOfFiles.size(); f++)
    {
      if(ivar ==0 || ivar ==1 || ivar == 3 || ivar == 4 || ivar == 9 || ivar ==10 )
      {
        hParton[f][ivar]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameParton.Data()));
        hParton_Clone[0][ivar]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameParton.Data()));
      }
      else
      {
        hParton[f][ivar]->GetXaxis()->SetTitle(TString::Format("%s", varNameParton.Data()));
        hParton_Clone[0][ivar]->GetXaxis()->SetTitle(TString::Format("%s", varNameParton.Data()));
      }
    }
    
  }//end of ivar loop
  

  //write histograms to file:
  for(int ivar=0; ivar<NVAR; ivar++)
  {
    for(int f=0; f<listOfFiles.size(); f++)
    {
       hParton[f][ivar]->Write(TString::Format("hParton_%s_%s",histoNames[f].Data(),varReco[ivar].Data()));      
    }
    hParton_Clone[0][ivar]->Write(TString::Format("hPartonScaledXSEC_%s", varReco[ivar].Data()));
    
  }
  listOfFiles.clear();
  XSEC.clear();
  histoNames.clear();
  
}
