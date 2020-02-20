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
  if(selection ==1) //signal ttbar mc
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
 
void GetMCHistograms(TString y="2017", int sel = 1, bool isNominalMC= false)
{
  globalIsNominalMC = isNominalMC;	
  year =y;
  initFilesMapping();
  selection = sel;
  initGlobals();  
  gStyle->SetOptStat(0);
  const int NVAR =9;
  TString varReco[NVAR]   = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1","mTop", "jetMassSoftDrop"};  
  
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
  TTree *trIN    = (TTree*)inf->Get("boosted/events");
  TTree *trINCnt = (TTree*)inf->Get("eventCounter/events"); 
  
  float NORM = ((TH1F*)inf->Get("eventCounter/GenEventWeight"))->GetSumOfWeights(); 
  weights.push_back(XSEC[f]/NORM);  
  
  int decade(0);
  int NN = trIN->GetEntries();
  
  int nJets,nLeptons;
  float genEvtWeight, genEvtWeightCnt;
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
  
  float partonEta[2], partonPhi[2], partonY[2], partonPt[2], partonE[2], partonMass[2];
  std::vector<float>  *deepAK8(0);
  std::vector<int> *addedIndexes = new std::vector<int>(0);
  std::vector<float> *jetBtagSub0(0), *jetBtagSub1(0);
  std::vector<float> *jetBtagSub0DCSVbb(0), *jetBtagSub1DCSVbb(0);
  std::vector<float> *jetBtagSub0DCSVbbb(0), *jetBtagSub1DCSVbbb(0);

  //particle
  std::vector<float> *genjetPt(0), *genjetY(0), *genjetEta(0), *genSoftDropMass(0), *genjetMassSoftDrop(0);
  std::vector<int> *nSubGenJets(0);
  int nJetsGen(0);
  float mJJGen(0), ptJJGen(0), yJJGen(0);
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
  trIN->SetBranchAddress("jetTtagCategory",&jetTtag);
  
  if(selection == 1)
  {
  trINCnt->SetBranchAddress("genEvtWeight" ,&genEvtWeightCnt);  
  trINCnt->SetBranchAddress("mTTbarParton" ,&mTTbarParton);
  trINCnt->SetBranchAddress("yTTbarParton" ,&yTTbarParton);
  trINCnt->SetBranchAddress("ptTTbarParton"  ,&ptTTbarParton);
  trINCnt->SetBranchAddress("ptTopParton"    ,&partonPt);
  trINCnt->SetBranchAddress("etaTopParton"   ,&partonEta);
  trINCnt->SetBranchAddress("yTopParton"     ,&partonY);
  trINCnt->SetBranchAddress("mTopParton"     ,&partonMass);
  //trINCnt->SetBranchAddress("partonMatchDR"  ,&partonMatchDR);
  //trINCnt->SetBranchAddress("partonMatchIdx" ,&partonMatchIdx);
  trINCnt->SetBranchAddress("phiTopParton"      ,&partonPhi);

  trIN->SetBranchAddress("mJJGen"			,&mJJGen);
  trIN->SetBranchAddress("ptJJGen"		,&ptJJGen);
  trIN->SetBranchAddress("yJJGen"			,&yJJGen);
  trIN->SetBranchAddress("genjetPt"		,&genjetPt);
  trIN->SetBranchAddress("genjetY"		,&genjetY);
  trIN->SetBranchAddress("genjetEta"		,&genjetEta);
  trIN->SetBranchAddress("nJetsGen"		,&nJetsGen);
  trIN->SetBranchAddress("genjetMassSoftDrop", &genjetMassSoftDrop);

  }

  float xReco(0), xParton(0), xParticle(0);
  std::vector<float> xRecoAll(0);
  std::vector<float> xPartonAll(0);
  std::vector<float> xParticleAll(0);

  int NBINS = 100;
  float x_min[NVAR] = {700,0,-3, 200, 200, 0.0, 0.0, 0, 0};
  float x_max[NVAR] = {5000, 1300, 3, 1500, 1500, 3, 3, 300,300};
  //book the histograms
  //histograms for Signal/QCD in CR 
  for(int ivar =0; ivar< NVAR; ivar++)
  {
    h[f][ivar] = new TH1F(TString::Format("h_%s_%s",histoNames[f].Data(),varReco[ivar].Data()), TString::Format("h_%s_%s",histoNames[f].Data(),varReco[ivar].Data()), NBINS, x_min[ivar], x_max[ivar]);
    hParton[f][ivar] = new TH1F(TString::Format("hParton_%s_%s",histoNames[f].Data(),varReco[ivar].Data()), TString::Format("hParton_%s_%s",histoNames[f].Data(),varReco[ivar].Data()), NBINS, x_min[ivar], x_max[ivar]);
    hParticle[f][ivar] = new TH1F(TString::Format("hParticle_%s_%s",histoNames[f].Data(),varReco[ivar].Data()), TString::Format("hParticle_%s_%s",histoNames[f].Data(),varReco[ivar].Data()), NBINS, x_min[ivar], x_max[ivar]);
  }  


  //event loop
  //åcout<<"Reading "<<NN<<" entries"<<endl;
  
  for(int iev=0;iev<NN;iev++) 
  {
    double progress = 10.0*iev/(1.0*NN);
    int k = TMath::FloorNint(progress); 
    if (k > decade) 
      cout<<10*k<<" %"<<endl;
    decade = k;
   trIN->GetEntry(iev);
  	if(nJets >1)
 	{
	  xRecoAll.clear();
   	xParticleAll.clear();
	  
    xRecoAll.push_back(mJJ);
	  xRecoAll.push_back(ptJJ);
	  xRecoAll.push_back(yJJ);
	  xRecoAll.push_back((*jetPt)[0]);
	  xRecoAll.push_back((*jetPt)[1]);

	  xRecoAll.push_back(fabs((*jetY)[0]));
	  xRecoAll.push_back(fabs((*jetY)[1]));
	  xRecoAll.push_back((*jetMassSoftDrop)[0]);

	  for(int ijet=1; ijet<nJets; ijet++)
	   xRecoAll.push_back((*jetMassSoftDrop)[ijet]);
	  
	   for(int ivar = 0; ivar <xRecoAll.size(); ivar ++)
	   {
	     xReco = xRecoAll[ivar];
	     if(ivar < 8)
	       h[f][ivar]->Fill(xReco, genEvtWeight);
	     else
	       h[f][8]->Fill(xReco,genEvtWeight);
	        
	   }

	   if(selection == 1)
	   {  
		 //cout<<"here"<<endl;
		  xParticleAll.push_back(mJJGen);
		  xParticleAll.push_back(ptJJGen);
		  xParticleAll.push_back(yJJGen);
		  xParticleAll.push_back((*genjetPt)[0]);
		  xParticleAll.push_back((*genjetPt)[1]);
		  xParticleAll.push_back(fabs((*genjetY)[0]));
		  xParticleAll.push_back(fabs((*genjetY)[1]));
		  xParticleAll.push_back((*genjetMassSoftDrop)[0]);

		  for(int ivar = 0; ivar <xParticleAll.size(); ivar ++)
	   	  {
		     xParticle = xParticleAll[ivar];  
		     hParticle[f][ivar]->Fill(xParticle, genEvtWeight);
	   	  }	
	   }

  	}//end of nJets
  } //---end of event loop 

  //now for parton
  
  if(selection ==1)
  {
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
		  xPartonAll.push_back(partonMass[0]);

		  for(int ivar = 0; ivar <xParticleAll.size(); ivar ++)
	   	  {
		     xParton = xPartonAll[ivar];
		     hParton[f][ivar]->Fill(xParton, genEvtWeightCnt);
		  }
	  }
  }

  }//----end of fileSize loop 
  
  //this will be used for the combined scale to XSEC histogram
  TH1F *h_Clone[listOfFiles.size()][NVAR], *hParton_Clone[listOfFiles.size()][NVAR], *hParticle_Clone[listOfFiles.size()][NVAR];

  for(int ivar= 0; ivar<NVAR; ivar++)
  {
    //for every slice
    for(int j=0; j<listOfFiles.size(); j++)
    {
    	//reco
        h_Clone[j][ivar]=(TH1F*)h[j][ivar]->Clone(TString::Format("h_%s_Clone",histoNames[j].Data()));       
        h_Clone[j][ivar]->Scale(weights[j]*LUMI); //this is the new to be added afterwards
        h[j][ivar]->Scale(weights[j]*LUMI);
      	//parton
      	if(selection ==1)
      	{
	        hParton_Clone[j][ivar]=(TH1F*)hParton[j][ivar]->Clone(TString::Format("hParton_%s_Clone",histoNames[j].Data()));       
	        hParton_Clone[j][ivar]->Scale(weights[j]*LUMI); //this is the new to be added afterwards
	        hParton[j][ivar]->Scale(weights[j]*LUMI);
	      	//particle
	        hParticle_Clone[j][ivar]=(TH1F*)hParticle[j][ivar]->Clone(TString::Format("hParticle_%s_Clone",histoNames[j].Data()));       
	        hParticle_Clone[j][ivar]->Scale(weights[j]*LUMI); //this is the new to be added afterwards
	        hParticle[j][ivar]->Scale(weights[j]*LUMI);
    	}
    }
    
    for(int j=1; j<listOfFiles.size(); j++)
    {      
      //Add them to get the whole phase space
      h_Clone[0][ivar]->Add(h_Clone[j][ivar]);
      if(selection ==1)
      {
      	hParton_Clone[0][ivar]->Add(hParton_Clone[j][ivar]);
      	hParticle_Clone[0][ivar]->Add(hParticle_Clone[j][ivar]);
  	  }
    } 
  
    
  }//end of ivar loop

  TFile *outFile;
  if(selection ==1)
  	if(!isNominalMC)
    	outFile = new TFile(TString::Format("HistoMC_TT_Mtt-700toInf_100bins_%s.root",year.Data()), "RECREATE");
    else
    	outFile = new TFile(TString::Format("HistoMC_TT_NominalMC_100bins_%s.root",year.Data()), "RECREATE");
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
    
  }//end of ivar loop
  

  //write histograms to file:
  for(int ivar=0; ivar<NVAR; ivar++)
  {
    for(int f=0; f<listOfFiles.size(); f++)
    {
      h[f][ivar]->Write(TString::Format("h_%s_%s",histoNames[f].Data(),varReco[ivar].Data()));
      if(selection ==1)
      {
      	hParton[f][ivar]->Write(TString::Format("hParton_%s_%s",histoNames[f].Data(),varReco[ivar].Data()));
      	hParticle[f][ivar]->Write(TString::Format("hParticle_%s_%s",histoNames[f].Data(),varReco[ivar].Data()));
      }
    }
    h_Clone[0][ivar]->Write(TString::Format("hScaledXSEC_%s", varReco[ivar].Data()));
    if(selection ==1)
    {
	  hParton_Clone[0][ivar]->Write(TString::Format("hPartonScaledXSEC_%s", varReco[ivar].Data()));
      hParticle_Clone[0][ivar]->Write(TString::Format("hParticleScaledXSEC_%s", varReco[ivar].Data())); 
    }
  }
  listOfFiles.clear();
  XSEC.clear();
  histoNames.clear();
  
}
