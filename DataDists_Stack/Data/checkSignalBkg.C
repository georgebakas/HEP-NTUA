#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"


using std::cin;
using std::cout; 
using std::endl;


std::vector<TString> listOfFiles;
std::vector<float> XSEC;
std::vector<TString> histoNames;
float LUMI = 35900;
TString eosPath; 
float deepCSVFloat = 0.6321;
int selection;

void initFileNames()
{
  
  if(selection ==1) //data
  {
	eosPath = "/eos/cms/store/user/gbakas/ttbar/JetHT/2016/";  
	listOfFiles.push_back("JetHT_Run2016-17Jul2018.root");
  }
  else if(selection ==2) //signal ttbar mc
  {
	eosPath = "/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/Signal/";  
    listOfFiles.push_back("TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");
    listOfFiles.push_back("TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");
  }
  else if(selection ==3) //bkg mc
  {
	eosPath = "/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/Bkg/";
	listOfFiles.push_back("QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
	listOfFiles.push_back("QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
	listOfFiles.push_back("QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
	listOfFiles.push_back("QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
	listOfFiles.push_back("QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
	listOfFiles.push_back("QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  }
  else if(selection ==4) //subdominant bkgs
  {
  	eosPath = "/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/Bkg/";
	listOfFiles.push_back("DYJetsToQQ_HT180_13TeV-madgraphMLM-pythia8.root");
	listOfFiles.push_back("WJetsToQQ_HT180_13TeV-madgraphMLM-pythia8.root");
	listOfFiles.push_back("ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
	listOfFiles.push_back("ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
	listOfFiles.push_back("ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root");
	listOfFiles.push_back("ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root");
  }
}

void initXsections()
{
  if(selection ==2)
  {
	XSEC.push_back(69.64);
	XSEC.push_back(16.74);
  }
  else if(selection ==3)
  {
	XSEC.push_back(3.67e+5);
	XSEC.push_back(2.94e+4);
	XSEC.push_back(6.524e+03);
	XSEC.push_back(1.064e+03);
	XSEC.push_back(121.5);
	XSEC.push_back(2.542e+01);
  }
  else if(selection ==4)
  {
	XSEC.push_back(1208);
	XSEC.push_back(3105);
	XSEC.push_back(38.09);
	XSEC.push_back(38.06);
	XSEC.push_back(35.6);
	XSEC.push_back(35.6);
  }
}

void initHistoNames()
{
	
  if(selection ==1) histoNames.push_back("Data_2016");
  else if(selection ==2)
  {
	histoNames.push_back("Signal_histo_Mtt_700_1000"); 
	histoNames.push_back("Signal_histo_Mtt_1000_Inf");
  }
  else if (selection ==3)
  {
	histoNames.push_back("QCD_histo_Mtt_300_500");
	histoNames.push_back("QCD_histo_Mtt_500_700");
	histoNames.push_back("QCD_histo_Mtt_700_1000");
	histoNames.push_back("QCD_histo_Mtt_1000_1500");
	histoNames.push_back("QCD_histo_Mtt_1500_2000");
	histoNames.push_back("QCD_histo_Mtt_2000_Inf");
  }
  else
  {
  	histoNames.push_back("DYJetsToQQ_HT180");
  	histoNames.push_back("WJetsToQQ_HT180");
  	histoNames.push_back("ST_tW_top_5f_inclusiveDecays");
  	histoNames.push_back("ST_tW_antitop_5f_inclusiveDecays");
  	histoNames.push_back("ST_t-channel_top_4f_inclusiveDecays");
  	histoNames.push_back("ST_t-channel_antitop_4f_inclusiveDecays");
  }
}

void initGlobals()
{
  initFileNames();
  initXsections();
  initHistoNames();
}

void checkSignalBkg(int tempSelection=1)
{
  
  selection = tempSelection;
  initGlobals();
  int fileSize = listOfFiles.size();
  vector<float> weights(0);
 
  
  const int NVAR =7;
  const int N_MJJ = 10;
  const int N_PTJJ = 9;
  const int N_YJJ = 8;
  const int N_PT = 10;
  const int N_JETY = 12;
  const int N_JETMASS = 40;

  TH1F *hSignal[2][fileSize][NVAR], *hBkg[2][fileSize][NVAR];

  bool saveTTagger= true;
  float selMvaCut=0.2;
  float floatBTag = 0.8838;
  bool isDeepCSV= true;
  
  int NBINS[NVAR] = {N_MJJ, N_PTJJ, N_YJJ, N_PT, N_PT ,/*N_JETY, N_JETY,*/ N_JETMASS, N_JETMASS};
  std::vector< std::vector <Float_t> > const BND = {{1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, 4000, 5000}, //mjj
											 {0,60,150,300,450,600,750,950,1100,1300}, //ptjj
											 {-2.4,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.4}, //yjj
											 {400,450,500,570,650,750,850,950,1100,1300,1500}, //jetPt0	
											 {400,450,500,570,650,750,850,950,1100,1300,1500}}; //jetPt1
											 //{0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}, //jetY0
											 //{0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}}; //jetY1
	
  TString varNameReco[NVAR]   = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", /*"jetY0", "jetY1",*/ "jetMassSoftDrop0", "jetMassSoftDrop1"}; 

  for(int f=0; f<fileSize; f++)
  {
	for(int ivar = 0; ivar<NVAR; ivar++)
	{
		int sizeBins = NBINS[ivar];
		if(ivar < NVAR-2)
		{
			float tempBND[NBINS[ivar]+1];
	  	    std::copy(BND[ivar].begin(), BND[ivar].end(), tempBND);

			hSignal[0][f][ivar] = new TH1F(TString::Format("%s_%s_SignalData_%s", varNameReco[ivar].Data(),"tTagger", histoNames[f].Data()), TString::Format("%s_%s_SignalData_%s" ,varNameReco[ivar].Data(),"tTagger", histoNames[f].Data()), sizeBins, tempBND);
			hSignal[1][f][ivar] = new TH1F(TString::Format("%s_%s_SignalData_%s", varNameReco[ivar].Data(),"oldMva", histoNames[f].Data()) , TString::Format("%s_%s_SignalData_%s" ,varNameReco[ivar].Data(),"oldMva", histoNames[f].Data()) , sizeBins, tempBND);
			
			hBkg[0][f][ivar] = new TH1F(TString::Format("%s_%s_BkgData_%s",varNameReco[ivar].Data(), "tTagger", histoNames[f].Data()), TString::Format("%s_%s_BkgData_%s" ,varNameReco[ivar].Data(),"tTagger", histoNames[f].Data()), sizeBins, tempBND);
			hBkg[1][f][ivar] = new TH1F(TString::Format("%s_%s_BkgData_%s",varNameReco[ivar].Data(), "oldMva", histoNames[f].Data()) , TString::Format("%s_%s_BkgData_%s" ,varNameReco[ivar].Data(),"oldMva", histoNames[f].Data()) , sizeBins, tempBND); 
		}
		else
		{ 
			hSignal[0][f][ivar] = new TH1F(TString::Format("%s_%s_SignalData_%s", varNameReco[ivar].Data(),"tTagger", histoNames[f].Data()), TString::Format("%s_%s_SignalData_%s" ,varNameReco[ivar].Data(),"tTagger", histoNames[f].Data()), 50,50,300);
			hSignal[1][f][ivar] = new TH1F(TString::Format("%s_%s_SignalData_%s", varNameReco[ivar].Data(),"oldMva", histoNames[f].Data()) , TString::Format("%s_%s_SignalData_%s" ,varNameReco[ivar].Data(),"oldMva", histoNames[f].Data()) , 50,50,300);
			
			hBkg[0][f][ivar] = new TH1F(TString::Format("%s_%s_BkgData_%s",varNameReco[ivar].Data(), "tTagger", histoNames[f].Data()), TString::Format("%s_%s_BkgData_%s" ,varNameReco[ivar].Data(),"tTagger", histoNames[f].Data()),  50,50,300);
			hBkg[1][f][ivar] = new TH1F(TString::Format("%s_%s_BkgData_%s",varNameReco[ivar].Data(), "oldMva", histoNames[f].Data()) , TString::Format("%s_%s_BkgData_%s" ,varNameReco[ivar].Data(),"oldMva", histoNames[f].Data()) , 50,50,300); 
			
		}
		hSignal[0][f][ivar]->GetXaxis()->SetTitle(varNameReco[ivar].Data());
		hSignal[1][f][ivar]->GetXaxis()->SetTitle(varNameReco[ivar].Data());
		hBkg[0][f][ivar]->GetXaxis()->SetTitle(varNameReco[ivar].Data());
		hBkg[1][f][ivar]->GetXaxis()->SetTitle(varNameReco[ivar].Data());
	}
    int nJets,nLeptons;
    float genEvtWeight;
    vector<bool>  *triggerBit = new vector<bool>;
    vector<bool>  *matchedJet(0);
    vector<float> *jetPt(0),*tau3(0),*tau2(0),*tau1(0);
    vector<float> *jetMassSub0(0), *jetMassSub1(0);
    vector<float> *jetMassSoftDrop(0);

    vector<float> *jetTtag(0);
    float mJJ(0), yJJ(0), ptJJ(0), mva(0);
    int  category(0);
    std::vector<float> *jetPhi(0), *jetEta(0), *jetY(0);
    std::vector<float> *jetBtagSub0(0), *jetBtagSub1(0);
    std::vector<float> *jetBtagSub0DCSVbb(0), *jetBtagSub1DCSVbb(0);
    std::vector<float> *jetBtagSub0DCSVbbb(0), *jetBtagSub1DCSVbbb(0);
	
    std::cout<<"Working in file: "<<listOfFiles[f]<<std::endl;
    TFile *file = TFile::Open(eosPath+listOfFiles[f]);
    TTree *trIN = (TTree*)file->Get("boosted/events");
    if(selection != 1)
    {
		float NORM = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights(); 
		weights.push_back(XSEC[f]/NORM); 	
	}
	//------- input tree --------------
    trIN->SetBranchAddress("nJets"          ,&nJets);
    trIN->SetBranchAddress("nLeptons"       ,&nLeptons);
    trIN->SetBranchAddress("jetPt"          ,&jetPt);
    trIN->SetBranchAddress("jetY"           ,&jetY);
    trIN->SetBranchAddress("jetEta"         ,&jetEta);
    trIN->SetBranchAddress("jetPhi"         ,&jetPhi);
    trIN->SetBranchAddress("jetTau3"        ,&tau3);
    trIN->SetBranchAddress("jetTau2"        ,&tau2);
    trIN->SetBranchAddress("jetTau1"        ,&tau1);
    trIN->SetBranchAddress("triggerBit"     ,&triggerBit);
    trIN->SetBranchAddress("jetMassSub0"    ,&jetMassSub0);
    trIN->SetBranchAddress("jetMassSub1"    ,&jetMassSub1);
    trIN->SetBranchAddress("genEvtWeight"   ,&genEvtWeight);
    trIN->SetBranchAddress("mJJ"   			,&mJJ);
    trIN->SetBranchAddress("yJJ"   			,&yJJ);
    trIN->SetBranchAddress("ptJJ"   		,&ptJJ);
    
	trIN->SetBranchAddress("jetBtagSub0"	,&jetBtagSub0);
    trIN->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1);
    trIN->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
  
    trIN->SetBranchAddress("mva"	  		,&mva);
    trIN->SetBranchAddress("category"	  	,&category);
    trIN->SetBranchAddress("jetTtagCategory",&jetTtag);
	//deepCSV
    trIN->SetBranchAddress("jetBtagSub0DCSVbb" ,&jetBtagSub0DCSVbb);
    trIN->SetBranchAddress("jetBtagSub1DCSVbb" ,&jetBtagSub1DCSVbb);
    trIN->SetBranchAddress("jetBtagSub0DCSVbbb",&jetBtagSub0DCSVbbb);
    trIN->SetBranchAddress("jetBtagSub1DCSVbbb",&jetBtagSub1DCSVbbb);
		
    int decade(0);
    int NN = trIN->GetEntries();
    //NN = 100000;
    std::cout<<"Entries: "<<NN<<std::endl;
	
	 
	std::vector<float> xReco(0);

    for(int iev=0;iev<NN;iev++)
    { 
		double progress = 10.0*iev/(1.0*NN);
      int k = TMath::FloorNint(progress); 
      if (k > decade) 
        cout<<10*k<<" %"<<endl;
      decade = k;
      trIN->GetEntry(iev);

      xReco.clear();
		
	     
	   //cout<<"--------------------"<<endl;
	   //cout<<xReco<<endl;

	  if(nJets >1)
	  {
      xReco.push_back(mJJ);
      xReco.push_back(ptJJ);
      xReco.push_back(yJJ);
      xReco.push_back((*jetPt)[0]);
      xReco.push_back((*jetPt)[1]);
      //xReco.push_back((*jetY)[0]);
      //xReco.push_back((*jetY)[0]);
      xReco.push_back((*jetMassSoftDrop)[0]);
      xReco.push_back((*jetMassSoftDrop)[1]);
		
		//cout<<"-----------------"<<endl;
		//cout<< xReco<<endl;
	  bool recoCuts; 
	  bool massCut = (*jetMassSoftDrop)[0] > 120 && (*jetMassSoftDrop)[0] < 220 && (*jetMassSoftDrop)[1] > 120 && (*jetMassSoftDrop)[1] < 220;
	  bool CSVv2Cut = ((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag) && ((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag);
	  bool tTaggerCut = (*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut;
	  recoCuts = nJets > 1 && fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1] <2.4) && (*jetPt)[0] > 400 && (*jetPt)[1] > 400 && mJJ > 1000  && nLeptons==0;
	  
	  bool triggerCut = (*triggerBit)[2];
	  
	  bool deepCSV = (((*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0])> deepCSVFloat) && 
					 (((*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1])> deepCSVFloat);
	  
      bool btagCut = deepCSV;
					  
	  float dCSVScoreSub0[2], dCSVScoreSub1[2];
	  dCSVScoreSub0[0] = (*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0];
	  dCSVScoreSub0[1] = (*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1];
	  dCSVScoreSub1[0] = (*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0];
	  dCSVScoreSub1[1] = (*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1];
	  bool revertBtagDeepCSV = (dCSVScoreSub0[0] < deepCSVFloat &&  dCSVScoreSub1[0] < deepCSVFloat) && (dCSVScoreSub0[1] < deepCSVFloat && dCSVScoreSub1[1] < deepCSVFloat);
	  


	  for(int ivar = 0; ivar < NVAR; ivar++)
	  {

	  	  if(ivar == 5 || ivar ==6) 
		 	massCut=(*jetMassSoftDrop)[0] > 50 && (*jetMassSoftDrop)[0] < 300 && (*jetMassSoftDrop)[1] > 50 && (*jetMassSoftDrop)[1] < 300;
		  //fill the denominators
		  //1. fill signal for tTagger
		  if(recoCuts && massCut && btagCut && tTaggerCut && triggerCut)
				hSignal[0][f][ivar]->Fill(xReco[ivar]);
		  //2. Fill signal for oldMva
		  if(recoCuts && massCut && category==2 && mva >0.8)
				hSignal[1][f][ivar]->Fill(xReco[ivar]);
		  //3. fill bkg for tTagger
		  if(recoCuts && massCut && revertBtagDeepCSV && tTaggerCut)
			  	hBkg[0][f][ivar]->Fill(xReco[ivar]);

		  //4. fill bkg for oldMva
		  if(recoCuts && massCut && category==0 && mva >0.8)
				hBkg[1][f][ivar]->Fill(xReco[ivar]);
	  }
	  
	  }//--- end of nJets if
    }//---end the event for
  cout<<"event loop ended"<<endl;
  

  for(int ivar = 0; ivar< NVAR; ivar++)
  {
	  hBkg[0][f][ivar]->SetLineColor(kRed);
	  hBkg[0][f][ivar]->SetMarkerColor(kRed);
	  hBkg[0][f][ivar]->SetMarkerStyle(22);
	  
	  hBkg[1][f][ivar]->SetLineColor(kRed);
	  hBkg[1][f][ivar]->SetMarkerColor(kRed);
	  hBkg[1][f][ivar]->SetMarkerStyle(23);
	  
	  hSignal[0][f][ivar]->SetLineColor(kBlue);
	  hSignal[0][f][ivar]->SetMarkerColor(kBlue);
	  hSignal[0][f][ivar]->SetMarkerStyle(24);
	  
	  hSignal[1][f][ivar]->SetLineColor(kBlue);
	  hSignal[1][f][ivar]->SetMarkerColor(kBlue);
	  hSignal[1][f][ivar]->SetMarkerStyle(25);
  }
  
  
} //-----end of file loop
  
  //for selection 2 and 3 we need to 
  //for tTagger and oldMva
  TFile *outfSignal = new TFile("SignalRegion_Data_MC.root","UPDATE");
  TFile *outfBkg = new TFile("ControlRegion_Data_MC.root","UPDATE");

  for(int ivar =0; ivar < NVAR; ivar++)
  {
	  if(selection !=1)
	  {
		  for(int i=0; i<2; i++)
		  {
			  for(int j =0;j <listOfFiles.size(); j++)
			  {
				  hSignal[i][j][ivar]->Scale(weights[j]*LUMI);
				  hBkg[i][j][ivar]->Scale(weights[j]*LUMI);
			  }
			  for(int j =1; j<listOfFiles.size(); j++)
			  {
				  hSignal[i][0][ivar]->Add(hSignal[i][j][ivar]);
				  hBkg[i][0][ivar]->Add(hBkg[i][j][ivar]);
			  }
		  }	
	  }
	  else
	  {
		  for(int j =0;j <listOfFiles.size(); j++)
		  {
			  hSignal[0][j][ivar]->Scale(1);
			  hSignal[1][j][ivar]->Scale(1);
			  hBkg[0][j][ivar]->Scale(1);
			  hBkg[1][j][ivar]->Scale(1);
		  }
	  }
  
  
  TString dataMC;
  if(selection ==1) dataMC = "data";
  else if(selection ==2) dataMC ="ttbarMC";
  else if(selection ==3) dataMC ="bkgMC";
  else if(selection ==4) dataMC ="subdominantBkgMC";
  
  outfSignal->cd();
  hSignal[0][0][ivar]->Write(TString::Format("%s_%sSignal_tTagger",varNameReco[ivar].Data(), dataMC.Data()));
  hSignal[1][0][ivar]->Write(TString::Format("%s_%sSignal_oldMva",varNameReco[ivar].Data(), dataMC.Data()));
  
  outfBkg->cd();
  hBkg[0][0][ivar]->Write(TString::Format("%s_%sBkg_tTagger",varNameReco[ivar].Data(), dataMC.Data()));
  hBkg[1][0][ivar]->Write(TString::Format("%s_%sBkg_oldMva",varNameReco[ivar].Data(), dataMC.Data()));
  
  }
  outfBkg->Close();
  outfSignal->Close();
  listOfFiles.clear();
  histoNames.clear();
  XSEC.clear();
}
