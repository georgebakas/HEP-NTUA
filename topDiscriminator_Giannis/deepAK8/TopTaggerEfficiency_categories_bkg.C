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


void TopTaggerEfficiency_categories_bkg(TString varName,TString varNameReco, int varNum, int sizeBins,float BND[], float selMvaCut=0.6, float floatBTag = 0.5803, bool computeMassSoftDrop = true)
{
  
	//for categories look comments below
	//for slices look the listOfFiles variable
  TH1F *h_Bkg[5][6];
  TH1F *h_JetMassSoftDrop[5][6];
cout<<"hello"<<endl;
  
  TString eosPath = "/eos/cms/store/user/gbakas/ttbar/topTagger/April19/";
  std::vector<TString> listOfFiles = {
	  "QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Copy.root",
	  "QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Copy.root",
	  "QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Copy.root",
	  "QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Copy.root",
	  "QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Copy.root",
	  "QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Copy.root"
  };
  
  float XSEC[6] = {3.67e+5,2.94e+4,6.524e+03,1.064e+03,121.5,2.542e+01};
 
  //for histograms
  // for var binning 
  //int sizeBins = 16; //for regular binning
  //const int sizeBins = 10;
  //float BND[11] = {1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, 4000, 5000};
  
  
  //const int sizeBins = 8;
  //double BND[sizeBins+1] = {1000, 1200, 1400, 1600, 1800, 2000, 2200, 2500, 3000};

  std::vector<float> weights;
  
  for(int f=0; f<listOfFiles.size(); f++)
  {
		//initialize histograms
		for(int i=0; i<sizeof(h_Bkg)/sizeof(*h_Bkg); i++)
		{
			h_Bkg[i][f] = new TH1F(TString::Format("Category_%i_Slice_%i vs mJJ", (i+1), (f+1)), TString::Format("Category_%i_Slice_%i vs mJJ", (i+1), (f+1)), sizeBins, BND);
			h_JetMassSoftDrop[i][f] = new TH1F(TString::Format("Category_%i_Slice_%i vs jetMassSoftDrop", (i+1), (f+1)), TString::Format("Category_%i_Slice_%i vs jetMassSoftDrop", (i+1), (f+1)),20,50,300);
			
		}

	  TFile *file = TFile::Open(eosPath+listOfFiles[f]);
	  TTree *trIN = (TTree*)file->Get("events");
	
	  float norm = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
	  float weight = XSEC[f]/norm;
	  weights.push_back(weight);

    int nJets,nLeptons;
    float genEvtWeight;
    vector<bool>  *bit(0);
    vector<bool>  *matchedJet(0);
    vector<float> *jetPt(0),*tau3(0),*tau2(0),*tau1(0);
    vector<float> *jetMassSub0(0), *jetMassSub1(0);
    vector<float> *jetMassSoftDrop(0);

    float pt_(-1),tau3_(-1),tau2_(-1),tau1_(-1),mass_(-1), jetMassSub0_(-1), jetMassSub1_(-1);
    float jetDr_(-1), mva(0);
    vector<float> *jetTtag(0);
    float mTTbarParton(0),mJJ(0), yJJ(0), ptJJ(0);
    int  category(0);
    //matching info
    vector<float> *jetPhi(0), *jetEta(0);
    std::vector<int> *addedIndexes = new std::vector<int>(0);
    std::vector<float> *jetBtagSub0(0), *jetBtagSub1(0);
    std::vector <bool>*triggerBit(0);       
    //------- input tree --------------
    trIN->SetBranchAddress("nJets"          ,&nJets);
    trIN->SetBranchAddress("nLeptons"       ,&nLeptons);
    trIN->SetBranchAddress("jetPt"          ,&jetPt);
    trIN->SetBranchAddress("jetEta"         ,&jetEta);
    trIN->SetBranchAddress("jetPhi"         ,&jetPhi);
    trIN->SetBranchAddress("jetTau3"        ,&tau3);
    trIN->SetBranchAddress("jetTau2"        ,&tau2);
    trIN->SetBranchAddress("jetTau1"        ,&tau1);
    trIN->SetBranchAddress("triggerBit"     ,&bit);
    trIN->SetBranchAddress("genEvtWeight"   ,&genEvtWeight);
    trIN->SetBranchAddress("jetMassSub0"    ,&jetMassSub0);
    trIN->SetBranchAddress("jetMassSub1"    ,&jetMassSub1);
    trIN->SetBranchAddress("mJJ"   			,&mJJ);
    trIN->SetBranchAddress("yJJ"   			,&yJJ);
    trIN->SetBranchAddress("ptJJ"   		,&ptJJ);
    trIN->SetBranchAddress("jetBtagSub0"	,&jetBtagSub0);
    trIN->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1);
    trIN->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
    trIN->SetBranchAddress("mTTbarParton"   ,&mTTbarParton);
	trIN->SetBranchAddress("triggerBit"		,&triggerBit);


    trIN->SetBranchAddress("mva"	  		,&mva);
    trIN->SetBranchAddress("category"	  	,&category);
  
    trIN->SetBranchAddress("jetTtagCategory",&jetTtag);
 

   int decade(0);
   int NN = trIN->GetEntries();
   std::vector<int> *jetMatchedIndexes = new std::vector<int>(0);
   std::vector<float> *jetMatchedDr = new std::vector<float>(0);
   
   bool isJetVar;
   if(varNum > 3)
   {
	   isJetVar = true;
   }
   //only reco in the bkg files
   float xReco(0);
   float xJetReco[2];
   
   cout<<"here"<<endl;
   cout<<"Reading "<<NN<<" entries"<<endl;
   for(int iev=0;iev<NN;iev++) 
   {
     double progress = 10.0*iev/(1.0*NN);
     int k = TMath::FloorNint(progress); 
     if (k > decade) 
       cout<<10*k<<" %"<<endl;
     decade = k;
     trIN->GetEntry(iev);
    //for every event fill the denominator of the efficiency (partonPt)
		if (varNum ==1)
		{
			xReco   = mJJ;
		}
		else if (varNum ==2)
		{
			xReco   = ptJJ;
		}
		else if (varNum ==3) 
		{
			xReco   = yJJ;
		}
		else if (varNum ==4)
		{
			if(nJets >1)
			{
				xJetReco[0] = (*jetPt)[0];
				xJetReco[1] = (*jetPt)[1];	
			}
			//cout<<"in"<<endl;
		}
		else if (varNum ==5)
		{
			if(nJets)
			{
				xJetReco[0] = (*jetEta)[0];
				xJetReco[1] = (*jetEta)[1];
			}
		}
	 //cout<<"ok"<<endl;
	 if(nJets > 1  && fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1] <2.4)
	   &&(*jetPt)[0] > 400 && (*jetPt)[1] > 400 && mJJ > 1000 && (*triggerBit)[2])
	 {	
	  //split them in 5 categories
	  //1. category where both are fully top tagged
	  if((*jetMassSoftDrop)[0] > 120 && (*jetMassSoftDrop)[0] < 220 && (*jetMassSoftDrop)[1] > 120 && (*jetMassSoftDrop)[1] < 220)
	  {
	    if((*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut && ((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag)
		   && ((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag))
		{
			 if(!isJetVar) 
			 {
				 h_Bkg[0][f]->Fill(xReco);
			 }
			 else
			 {
				h_Bkg[0][f]->Fill(xJetReco[0]);
				h_Bkg[0][f]->Fill(xJetReco[1]);
			 }
		}
		  //------ category 2 and 3 for new jetTtager---------------------------------------------------------
		if((*jetTtag)[0] > selMvaCut  
		   && ((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag))
		//meaning at least 1 fully tagged
		{
			//2. This category has 1 jet full top tagged and 1 mva from jet [1]
		  if((*jetTtag)[1] > selMvaCut) 
		  {
			 if(!isJetVar) 
			 {
				 h_Bkg[1][f]->Fill(xReco);
			 }
			 else
			 {
				h_Bkg[1][f]->Fill(xJetReco[0]);
				h_Bkg[1][f]->Fill(xJetReco[1]);
			 }
		  }
		  //or just b tagging
		  if(((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag)) 
		  {
  			 if(!isJetVar) 
			 {
				 h_Bkg[2][f]->Fill(xReco);
			 }
			 else
			 {
				h_Bkg[2][f]->Fill(xJetReco[0]);
				h_Bkg[2][f]->Fill(xJetReco[1]);
			 }
		  }
		}
		else if((*jetTtag)[1] > selMvaCut  && ((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag )) //the other at least is fully tagged
		{
		  //2. This category has 1 jet full top tagged and 1 mva from jet [1]
		  if( (*jetTtag)[0] > selMvaCut) 
		  {
		  	 if(!isJetVar) 
			 {
				 h_Bkg[1][f]->Fill(xReco);
			 }
			 else
			 {
				h_Bkg[1][f]->Fill(xJetReco[0]);
				h_Bkg[1][f]->Fill(xJetReco[1]);
			 }
		  }
		  //or just the btagging
		  if(((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag)) 
		  {
		     if(!isJetVar) 
			 {
				 h_Bkg[2][f]->Fill(xReco);
			 }
			 else
			 {
				h_Bkg[2][f]->Fill(xJetReco[0]);
				h_Bkg[2][f]->Fill(xJetReco[1]);
			 }
		  }
		}

		//------------------------------------------------------------ 4th category for jetTtager -----------------------------------------------------------------
		//4th category is both jets passing only mva
		if((*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut)
		{
		 	 if(!isJetVar) 
			 {
				 h_Bkg[3][f]->Fill(xReco);
			 }
			 else
			 {
				h_Bkg[3][f]->Fill(xJetReco[0]);
				h_Bkg[3][f]->Fill(xJetReco[1]);
			 }
		}
			
		//------------------------------------------------------------ 5th category for jetTtager -----------------------------------------------------------------
		//5th category is both jets passing only btagging
		if(((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag) && ((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag))
		{
		 	if(!isJetVar) 
			 {
				 h_Bkg[4][f]->Fill(xReco);
			 }
			 else
			 {
				h_Bkg[4][f]->Fill(xJetReco[0]);
				h_Bkg[4][f]->Fill(xJetReco[1]);
			 }
		}

	  }
	  
	  //-----------just for Jet mass soft drop because we need the full spectrum-----------------------------------------------------------------------------
	  if(computeMassSoftDrop)
	  {
		  if((*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut && ((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag)
			   && ((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag))
			{
				h_JetMassSoftDrop[0][f]->Fill((*jetMassSoftDrop)[0]);
				h_JetMassSoftDrop[0][f]->Fill((*jetMassSoftDrop)[1]);
			}
			  //------ category 2 and 3 for new jetTtager---------------------------------------------------------
			if((*jetTtag)[0] > selMvaCut  
			   && ((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag))
			//meaning at least 1 fully tagged
			{
				//2. This category has 1 jet full top tagged and 1 mva from jet [1]
			  if((*jetTtag)[1] > selMvaCut) 
			  {
				 h_JetMassSoftDrop[1][f]->Fill((*jetMassSoftDrop)[0]);
				 h_JetMassSoftDrop[1][f]->Fill((*jetMassSoftDrop)[1]);
			  }
			  //or just b tagging
			  if(((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag)) 
			  {
				h_JetMassSoftDrop[2][f]->Fill((*jetMassSoftDrop)[0]);
				h_JetMassSoftDrop[2][f]->Fill((*jetMassSoftDrop)[1]);
				
			  }
			}
			else if((*jetTtag)[1] > selMvaCut  && ((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag )) //the other at least is fully tagged
			{
			  //2. This category has 1 jet full top tagged and 1 mva from jet [1]
			  if( (*jetTtag)[0] > selMvaCut) 
			  {
				h_JetMassSoftDrop[1][f]->Fill((*jetMassSoftDrop)[0]);
				h_JetMassSoftDrop[1][f]->Fill((*jetMassSoftDrop)[1]);
			  }
			  //or just the btagging
			  if(((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag)) 
			  {
				h_JetMassSoftDrop[2][f]->Fill((*jetMassSoftDrop)[0]);
				h_JetMassSoftDrop[2][f]->Fill((*jetMassSoftDrop)[1]);
			  }
			}

			//------------------------------------------------------------ 4th category for jetTtager -----------------------------------------------------------------
			//4th category is both jets passing only mva
			if((*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut)
			{
				h_JetMassSoftDrop[3][f]->Fill((*jetMassSoftDrop)[0]);
				h_JetMassSoftDrop[3][f]->Fill((*jetMassSoftDrop)[1]);
			}
				
			//------------------------------------------------------------ 5th category for jetTtager -----------------------------------------------------------------
			//5th category is both jets passing only btagging
			if(((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag) && ((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag))
			{
				h_JetMassSoftDrop[4][f]->Fill((*jetMassSoftDrop)[0]);
				h_JetMassSoftDrop[4][f]->Fill((*jetMassSoftDrop)[1]);
			}
	  }
	  
	}
   }
      // for---- all events 
   file->Close();
   
  }//for all filfes

	//Scale the histograms to get the yield
	//loop over categories
	//Do the same for signal over Bkg 
	
	for(int i=0; i<sizeof(h_Bkg)/sizeof(*h_Bkg); i++)
	{
		//for every slice
		for(int j=0; j<sizeof(*h_Bkg)/sizeof(*h_Bkg[0]); j++)
		{
      //scale it to get the yield
			std::cout<<"Sacling "<<i<<" "<<j<<std::endl;
			if(h_Bkg[i][j]->GetSumw2N() == 0)
			{
				h_Bkg[i][j]->Sumw2(kTRUE);
			}
			h_Bkg[i][j]->Scale(weights[j]);
		}
    
		for(int j=1; j<sizeof(*h_Bkg)/sizeof(*h_Bkg[0]); j++)
		{
			//Add them to get the whole phase space
			std::cout<<"Appending "<<i<<" "<<j<<std::endl;
			h_Bkg[i][0]->Add(h_Bkg[i][j]);
		}
	}
	
	//scale the jetMassSoftDrop and find the yield 
	if(computeMassSoftDrop)
	{
		for(int i=0; i<sizeof(h_JetMassSoftDrop)/sizeof(*h_JetMassSoftDrop); i++)
		{
			//for every slice
			for(int j=0; j<sizeof(*h_JetMassSoftDrop)/sizeof(*h_JetMassSoftDrop[0]); j++)
			{
		  //scale it to get the yield
				std::cout<<"Sacling "<<i<<" "<<j<<std::endl;
				if(h_JetMassSoftDrop[i][j]->GetSumw2N() == 0)
				{
					h_JetMassSoftDrop[i][j]->Sumw2(kTRUE);
				}
				h_JetMassSoftDrop[i][j]->Scale(weights[j]*35900);
			}
		
			for(int j=1; j<sizeof(*h_JetMassSoftDrop)/sizeof(*h_JetMassSoftDrop[0]); j++)
			{
				//Add them to get the whole phase space
				std::cout<<"Appending "<<i<<" "<<j<<std::endl;
				h_JetMassSoftDrop[i][0]->Add(h_JetMassSoftDrop[i][j]);
			}
		}
	}
	//Now the first entry of every row h_Eff[i][0] contains the yield for the whole space.

	

  std::vector<Color_t> cols = {kBlue, kBlack, kRed, kMagenta, kGreen};
  std::vector<TString> categoryTitles = {"Both b and top tagged",
	                                     "One b and top tagged and one top tagged",
										 "One b and top tagged and one b tagged",
										 "Both top tagged",
										 "Both b tagged"};
  
  TEfficiency *efficiency[5];
  TString stringBTag;
  if (floatBTag < 0.6) stringBTag = "Loose";
  else  stringBTag = "Medium";
  TFile *outf = TFile::Open(TString::Format("Efficiencies_mva%0.1f_btag%s.root",selMvaCut, stringBTag.Data()), "UPDATE");
  

  outf->cd();
  if(varNum ==1 || varNum ==2 ||varNum ==4) 
  {
	  h_Bkg[0][0]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameReco.Data()));
  }
  else 
  {
	  h_Bkg[0][0]->GetXaxis()->SetTitle(TString::Format("%s", varNameReco.Data()));
  }
  h_Bkg[0][0]->Write(TString::Format("hBkg_Category_1_%s", varNameReco.Data()));
  
  if(computeMassSoftDrop)
  {
	h_JetMassSoftDrop[0][0]->GetXaxis()->SetTitle("jetMassSoftDrop (GeV)");
	h_JetMassSoftDrop[0][0]->Write("hBkg_JetMassSoftDrop_Category_1");  
  }

  for(int i=1; i<sizeof(h_Bkg)/sizeof(*h_Bkg); i++)
  {
	//efficiency[i]  = new TEfficiency(*h_Eff[i][0], *h_denominator[0]);
	//efficiency[i]->SetStatisticOption(TEfficiency::kFNormal);
	//efficiency[i]->SetUseWeightedEvents();
	//std::cout<<h_Eff[i][0]->GetEntries()<<std::endl;
	//can_eff[i] = new TCanvas(TString::Format("Efficiency canvas for category: %d",i+1),TString::Format("Efficiency canvas for category: %d",i+1), 900,600);
	//efficiency[i]->SetTitle(TString::Format("Efficiency vs mTTbarParton for mva > %.1f for category: %d ;mTTbarParton (GeV);Efficiency", selMvaCut, i+1));
	//efficiency[i]->SetLineColor(cols[i]);
	//efficiency[i]->Draw("SAME");
	//h_Eff[i][0]->Draw("peSAME");
	outf->cd();
	//efficiency[i]->Write(TString::Format("Bkg_Category_%i", (i+1)));
	if(varNum ==1 || varNum ==2 ||varNum ==4) 
	{
	  h_Bkg[i][0]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameReco.Data()));
	}
	else 
	{
	  h_Bkg[i][0]->GetXaxis()->SetTitle(TString::Format("%s", varNameReco.Data()));
	}
	
	if(computeMassSoftDrop)
	{
		h_JetMassSoftDrop[i][0]->SetLineColor(cols[i]);
		h_JetMassSoftDrop[i][0]->SetMarkerColor(cols[i]);
		h_JetMassSoftDrop[i][0]->GetXaxis()->SetTitle("jetMassSoftDrop (GeV)");
		h_JetMassSoftDrop[i][0]->Write(TString::Format("hBkg_JetMassSoftDrop_Category_%i", (i+1)));
	}
	
	h_Bkg[i][0]->Write(TString::Format("hBkg_Category_%i_%s", (i+1), varNameReco.Data()));
	//c2->cd();
	//h_Eff[i][0]->SetLineColor(cols[i]);
	//h_Eff[i][0]->SetMarkerColor(cols[i]);
	
	//h_Eff[i][0]->Draw("pesame");
	//leg1->AddEntry(efficiency[i], categoryTitles[i], "l");
  }


	//h_denominator[0]->SetMarkerColor(kRed);
	//h_denominator[0]->Draw("pesame");
  
  outf->Close();
	
}





