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


void TopTaggerEfficiency_categories_signal_Mtt(TString varName,TString varNameReco, int varNum, int sizeBins,float BND[], float selMvaCut=0.6, float floatBTag = 0.5803, bool computeMassSoftDrop = true)
{
  TH1F *h_Eff[5][2];
  TH1F *h_denominator[2];
  TH1F *h_Mtt[5][2];
  TH1F *h_JetMassSoftDrop[5][2];
  
  TString eosPath = "/eos/cms/store/user/gbakas/ttbar/topTagger/April19/";
  std::vector<TString> listOfFiles = {
	  "TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8_Copy.root",
	  "TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8_Copy.root"
  };
  
  float XSEC[2] = {69.64 ,16.74};
  std::vector<float> weights;
  float LUMI = 35900;
//  const int sizeBins = 10;
//  float BND[11] = {1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, 4000, 5000};
  
  for(int f=0; f<listOfFiles.size(); f++)
  {
		for(int i=0; i<sizeof(h_Eff)/sizeof(*h_Eff); i++)
		{
			h_Eff[i][f] = new TH1F(TString::Format("Category_%i_Slice_%i", (i+1), (f+1)), TString::Format("Category_%i_Slice_%i", (i+1), (f+1)), sizeBins, BND);
			h_Mtt[i][f] = new TH1F(TString::Format("hSigMtt_%i_Slice_%i", (i+1), (f+1)), TString::Format("hSigMtt_%i_Slice_%i", (i+1), (f+1)), sizeBins, BND);
			h_Eff[i][f]->Sumw2();
			h_JetMassSoftDrop[i][f] = new TH1F(TString::Format("JetMassSoftDrop Category_%i_Slice_%i", (i+1), (f+1)), TString::Format("JetMassSoftDrop Category_%i_Slice_%i", (i+1), (f+1)), 20, 50,300);
		}
		h_denominator[f] = new TH1F(TString::Format("Denom_Slice_%i", (f+1)), TString::Format("Denom_Slice_%i", (f+1)), sizeBins, BND);
		h_denominator[f]->Sumw2();

	  TFile *file = TFile::Open(eosPath+listOfFiles[f]);
    TTree *trIN = (TTree*)file->Get("events");
	
	  float norm = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
	  float weight = XSEC[f]/norm;
	  weights.push_back(weight);
	
    int nJets,nLeptons;
    float genEvtWeight;
    vector<bool>  *bit(0);
    vector<bool>  *matchedJet(0);
    vector<float> *jetPt(0),*tau3(0),*tau2(0),*tau1(0),*mass(0);
    vector<float> *jetMassSub0(0), *jetMassSub1(0);
    vector<float> *jetMassSoftDrop(0), *partonEta(0), *partonPt(0);

    float pt_(-1),tau3_(-1),tau2_(-1),tau1_(-1),mass_(-1), jetMassSub0_(-1), jetMassSub1_(-1);
    float jetDr_(-1), mva(0), mTTbarParton(0), yTTbarParton(0), ptTTbarParton(0);
    vector<float> *jetTtag(0);
    float mJJ(0), yJJ(0), ptJJ(0);
    int  category(0);
    vector<float> *jetPhi(0), *jetEta(0);
    std::vector<int> *addedIndexes = new std::vector<int>(0);
    std::vector<float> *jetBtagSub0(0), *jetBtagSub1(0);
          
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
	trIN->SetBranchAddress("mTTbarParton"	,&mTTbarParton);
	trIN->SetBranchAddress("yTTbarParton"	,&yTTbarParton);
	trIN->SetBranchAddress("ptTTbarParton"	,&ptTTbarParton);
	trIN->SetBranchAddress("partonPt"	    ,&partonPt);
	trIN->SetBranchAddress("partonEta"		,&partonEta);
	
    trIN->SetBranchAddress("jetBtagSub0"	,&jetBtagSub0);
    trIN->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1);
    trIN->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
  
    trIN->SetBranchAddress("mva"	  		,&mva);
    trIN->SetBranchAddress("category"	  	,&category);
  
    trIN->SetBranchAddress("jetTtagCategory",&jetTtag);

	
   int decade(0);
   int NN = trIN->GetEntries();
   std::vector<int> *jetMatchedIndexes = new std::vector<int>(0);
   std::vector<float> *jetMatchedDr = new std::vector<float>(0);
   //NN = 500000;
   int counter =0;
   
   bool isJetVar;
   if(varNum > 3)
   {
	   isJetVar = true;
   }
   float xParton(0), xReco(0);
   float xJetParton[2], xJetReco[2];
   
    cout<<"Reading "<<NN<<" entries"<<endl;
    for(int iev=0;iev<NN;iev++) 
    {
      double progress = 10.0*iev/(1.0*NN);
      int k = TMath::FloorNint(progress); 
      if (k > decade) 
        cout<<10*k<<" %"<<endl;
      decade = k;
      trIN->GetEntry(iev);
	  
	  //cout<< (*jetEta)[0]<<endl;
	  //cout<< (*jetEta)[1]<<endl;
		
		if (varNum ==1)
		{
			xParton = mTTbarParton;
			xReco   = mJJ;
		}
		else if (varNum ==2)
		{
			xParton = ptTTbarParton;
			xReco   = ptJJ;
		}
		else if (varNum ==3) 
		{
			xParton = yTTbarParton;
			xReco   = yJJ;
		}
		else if (varNum ==4)
		{
			
			//xJetParton.push_back((*partonPt)[0]);
			//xJetParton.push_back((*partonPt)[1]);
			
			//xJetReco.push_back((*jetPt)[0]);
			//xJetReco.push_back((*jetPt)[1]);
			xJetParton[0] = (*partonPt)[0];
			xJetParton[1] = (*partonPt)[1];
			
			xJetReco[0] = (*jetPt)[0];
			xJetReco[1] = (*jetPt)[1];
			
		}
		else if (varNum ==5)
		{
			
			//xJetParton.push_back((*partonEta)[0]);
			//xJetParton.push_back((*partonEta)[1]);
			
			//xJetReco.push_back((*jetEta)[0]);
			//xJetReco.push_back((*jetEta)[1]);
			
			xJetParton[0] = (*partonEta)[0];
			xJetParton[1] = (*partonEta)[1];
			
			xJetReco[0] = (*jetEta)[0];
			xJetReco[1] = (*jetEta)[1];
		}
	
	  //cout<<"ok"<<endl;
	  if(nJets > 1  && fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1] <2.4)
		&& (*jetPt)[0] > 400 && (*jetPt)[1] > 400 && mTTbarParton > 1000)
	  {	
		//cout<<"ok"<<endl;
		 
		 if(!isJetVar) h_denominator[f]->Fill(xParton);
		 else
		 {
			 h_denominator[f]->Fill(xJetParton[0]);
			 h_denominator[f]->Fill(xJetParton[1]);
		 }
		
		//split them in 5 categories
		//1. category where both are fully top tagged
		if((*jetMassSoftDrop)[0] > 120 && (*jetMassSoftDrop)[0] < 220 && (*jetMassSoftDrop)[1] > 120 && (*jetMassSoftDrop)[1] < 220)
		{
	      if((*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut && ((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag)
		    && ((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag))
		  {
	
			 
			 if(!isJetVar) 
			 {
				 h_Eff[0][f]->Fill(xParton);
				 h_Mtt[0][f]->Fill(xReco);
			 }
			 else
			 {
				h_Eff[0][f]->Fill(xJetParton[0]);
				h_Eff[0][f]->Fill(xJetParton[1]);
				
				h_Mtt[0][f]->Fill(xJetReco[0]);
				h_Mtt[0][f]->Fill(xJetReco[1]);
			 }
		    
		  }
    
		    //------ category 2 and 3 for new jetTtager---------------------------------------------------------
		  if((*jetTtag)[0] > selMvaCut  &&
		     ((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag))
			 //meaning at least 1 fully tagged
		  {
		  	//2. This category has 1 jet full top tagged and 1 mva from jet [1]
		  	if((*jetTtag)[1] > selMvaCut) 
		  	{
			
  			    if(!isJetVar) 
				{
					h_Eff[1][f]->Fill(xParton);
					h_Mtt[1][f]->Fill(xReco);
				}
				else
				{
					h_Eff[1][f]->Fill(xJetParton[0]);
					h_Eff[1][f]->Fill(xJetParton[1]);
					
					h_Mtt[1][f]->Fill(xJetReco[0]);
					h_Mtt[1][f]->Fill(xJetReco[1]);
				}
		  
		  	}
		  	//or just b tagging
		  	if(((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag)) 
		  	{
			 
  			    if(!isJetVar) 
				{
					h_Eff[2][f]->Fill(xParton);
					h_Mtt[2][f]->Fill(xReco);
				}
				else
				{
					h_Eff[2][f]->Fill(xJetParton[0]);
					h_Eff[2][f]->Fill(xJetParton[1]);
					
					h_Mtt[2][f]->Fill(xJetReco[0]);
					h_Mtt[2][f]->Fill(xJetReco[1]);
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
					h_Eff[1][f]->Fill(xParton);
					h_Mtt[1][f]->Fill(xReco);
				}
				else
				{
					h_Eff[1][f]->Fill(xJetParton[0]);
					h_Eff[1][f]->Fill(xJetParton[1]);
					
					h_Mtt[1][f]->Fill(xJetReco[0]);
					h_Mtt[1][f]->Fill(xJetReco[1]);
				}
		  	}
		  	//or just the btagging
		  	if(((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag)) 
		  	{
			  
		  		if(!isJetVar) 
				{
					h_Eff[2][f]->Fill(xParton);
					h_Mtt[2][f]->Fill(xReco);
				}
				else
				{
					h_Eff[2][f]->Fill(xJetParton[0]);
					h_Eff[2][f]->Fill(xJetParton[1]);
					
					h_Mtt[2][f]->Fill(xJetReco[0]);
					h_Mtt[2][f]->Fill(xJetReco[1]);
				}
			  
		  	}
		  }

			//------------------------------------------------------------ 4th category for jetTtager -----------------------------------------------------------------
			//4th category is both jets passing only mva
			if((*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut)
			{
			  
		  		if(!isJetVar) 
				{
					h_Eff[3][f]->Fill(xParton);
					h_Mtt[3][f]->Fill(xReco);
				}
				else
				{
					h_Eff[3][f]->Fill(xJetParton[0]);
					h_Eff[3][f]->Fill(xJetParton[1]);
					
					h_Mtt[3][f]->Fill((xJetReco)[0]);
					h_Mtt[3][f]->Fill((xJetReco)[1]);
				}
			}
			
			//------------------------------------------------------------ 5th category for jetTtager -----------------------------------------------------------------
			//5th category is both jets passing only btagging
			if(((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag) && ((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag))
			{
		  		if(!isJetVar) 
				{
					h_Eff[4][f]->Fill(xParton);
					h_Mtt[4][f]->Fill(xReco);
				}
				else
				{
					h_Eff[4][f]->Fill(xJetParton[0]);
					h_Eff[4][f]->Fill(xJetParton[1]);
					
					h_Mtt[4][f]->Fill(xJetReco[0]);
					h_Mtt[4][f]->Fill(xJetReco[1]);
				}
			}
		  } //soft drop mass window 
		  
		  
		  //--------------------------------jet Mass Soft Drop----------------------------------------------------------------------------------------------
		  
		  if(computeMassSoftDrop)
		  {
			  if((*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut && ((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag)
				&& ((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag))
			  {
		
				  h_JetMassSoftDrop[0][f]->Fill((*jetMassSoftDrop)[0]);
				  h_JetMassSoftDrop[0][f]->Fill((*jetMassSoftDrop)[1]);
				
			  }
		
				//------ category 2 and 3 for new jetTtager---------------------------------------------------------
			  if((*jetTtag)[0] > selMvaCut  &&
				 ((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag))
				 //meaning at least 1 fully tagged
			  {
				//2. This category has 1 jet full top tagged and 1 mva from jet [1]
				if((*jetTtag)[1] > selMvaCut) 
				{
				
					h_JetMassSoftDrop[1][f]->Fill((*jetMassSoftDrop)[0]);
					h_JetMassSoftDrop[1][f]->Fill((*jetMassSoftDrop)[1]);
			  
				  //cat2++;
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
				//cat2++;
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
  for(int i=0; i<sizeof(h_Eff)/sizeof(*h_Eff); i++)
	{
		//for every slice
		for(int j=0; j<sizeof(*h_Eff)/sizeof(*h_Eff[0]); j++)
		{
			h_Eff[i][j]->Scale(weights[j]);
		}
    
		for(int j=1; j<sizeof(*h_Eff)/sizeof(*h_Eff[0]); j++)
		{
			//Add them to get the whole phase space
			std::cout<<"Appending "<<i<<" "<<j<<std::endl;
			h_Eff[i][0]->Add(h_Eff[i][j]);
		}

	}
	
	//do same for input to be used later for signal mtt over bkg
	for(int i=0; i<sizeof(h_Mtt)/sizeof(*h_Mtt); i++)
	{
		//for every slice
		for(int j=0; j<sizeof(*h_Mtt)/sizeof(*h_Mtt[0]); j++)
		{
      //scale it to get the yield
			//std::cout<<"Sacling "<<i<<" "<<j<<std::endl;
			h_Mtt[i][j]->Scale(weights[j]);
		}
    
		for(int j=1; j<sizeof(*h_Mtt)/sizeof(*h_Mtt[0]); j++)
		{
			//Add them to get the whole phase space
			std::cout<<"Appending "<<i<<" "<<j<<std::endl;
			h_Mtt[i][0]->Add(h_Mtt[i][j]);
		}

	}
	if(computeMassSoftDrop)
	{
		//do same for input to be used later for jetMassSoftDrop
		for(int i=0; i<sizeof(h_JetMassSoftDrop)/sizeof(*h_JetMassSoftDrop); i++)
		{
			//for every slice
			for(int j=0; j<sizeof(*h_JetMassSoftDrop)/sizeof(*h_JetMassSoftDrop[0]); j++)
			{
		  //scale it to get the yield
				//std::cout<<"Sacling "<<i<<" "<<j<<std::endl;
				h_JetMassSoftDrop[i][j]->Scale(weights[j]*LUMI);
			}
		
			for(int j=1; j<sizeof(*h_JetMassSoftDrop)/sizeof(*h_JetMassSoftDrop[0]); j++)
			{
				//Add them to get the whole phase space
				std::cout<<"Appending "<<i<<" "<<j<<std::endl;
				h_JetMassSoftDrop[i][0]->Add(h_JetMassSoftDrop[i][j]);
				h_JetMassSoftDrop[i][0]->GetXaxis()->SetTitle("jetMassSoftDrop (GeV)");
			}
		
		}
	}

	h_denominator[0]->Scale(weights[0]);
	for(int i=1; i<sizeof(h_denominator)/sizeof(*h_denominator); i++) 
	{
		h_denominator[i]->Scale(weights[i]);
		h_denominator[0]->Add(h_denominator[i]);
	}


  std::vector<Color_t> cols = {kBlue, kBlack, kRed, kMagenta, kGreen};
  
  std::vector<TString> categoryTitles = {"Both b and top tagged",
	                                     "One b and top tagged and one top tagged",
										 "One b and top tagged and one b tagged",
										 "Both top tagged",
										 "Both b tagged"};
  
  TEfficiency *efficiency[5];
  
  TString stringBTag;
  if (floatBTag == (float)0.5803) stringBTag = "Loose";
  else if(floatBTag == (float)0.8838) stringBTag = "Medium";
  TFile *outf = TFile::Open(TString::Format("Efficiencies_mva%0.1f_btag%s.root",selMvaCut, stringBTag.Data()), "UPDATE");
  
  efficiency[0]  = new TEfficiency(*h_Eff[0][0], *h_denominator[0]);
  efficiency[0]->SetStatisticOption(TEfficiency::kFNormal);
  efficiency[0]->SetUseWeightedEvents();
  TCanvas *c1 = new TCanvas(TString::Format("Efficiency canvas for category: %d",1),TString::Format("Efficiency canvas for category: %d",1), 900,600);
  if(computeMassSoftDrop)
  {
	h_JetMassSoftDrop[0][0]->GetXaxis()->SetTitle("jetMassSoftDrop (GeV)");
	h_JetMassSoftDrop[0][0]->Write(TString::Format("hMtt_JetMassSoftDrop_Sig_Category_%i", 1));
  }

  
  //for the vars mass, pt of each top and the combined ptTTbar
	
  efficiency[0]->SetLineColor(cols[0]);
  efficiency[0]->Draw();
  
  if (varNum ==1 || varNum ==2 || varNum ==4) 
  {
	  efficiency[0]->SetTitle(TString::Format("Efficiency vs %s for mva > %.1f;%s (GeV);Efficiency",varName.Data(),selMvaCut, varName.Data()));
	  h_Mtt[0][0]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameReco.Data()));
	  h_denominator[0]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varName.Data()));
  }
  else
  {
	  efficiency[0]->SetTitle(TString::Format("Efficiency vs %s for mva > %.1f;%s;Efficiency",varName.Data(),selMvaCut, varName.Data()));
	  h_Mtt[0][0]->GetXaxis()->SetTitle(TString::Format("%s", varNameReco.Data()));
	  h_denominator[0]->GetXaxis()->SetTitle(TString::Format("%s", varName.Data()));	  

  }	  
	  
	  
  efficiency[0]->Write(TString::Format("Mtt_Sig_Category_%i_%s", 1, varName.Data()));
  h_Mtt[0][0]->Write(TString::Format("hMtt_Sig_Category_%i_%s", 1, varNameReco.Data()));
	
  TLegend *leg1 = new TLegend(0.5, 0.5, 0.7, 0.7);
  leg1->AddEntry(efficiency[0], categoryTitles[0], "l");

  for(int i=1; i<sizeof(h_Eff)/sizeof(*h_Eff); i++)
  {
	  if (varNum ==1 || varNum ==2 || varNum ==4) h_Eff[i][0]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varName.Data()));
	  else h_Eff[i][0]->GetXaxis()->SetTitle(TString::Format("%s", varName.Data()));
	  
	  efficiency[i]  = new TEfficiency(*h_Eff[i][0], *h_denominator[0]);
	  efficiency[i]->SetStatisticOption(TEfficiency::kFNormal);
	  efficiency[i]->SetUseWeightedEvents();
	  if(computeMassSoftDrop)
	  {
        h_JetMassSoftDrop[i][0]->SetLineColor(cols[i]); 
		h_JetMassSoftDrop[i][0]->SetMarkerColor(cols[i]);
		h_JetMassSoftDrop[i][0]->GetXaxis()->SetTitle("jetMassSoftDrop (GeV)");
		h_JetMassSoftDrop[i][0]->Write(TString::Format("hMtt_JetMassSoftDrop_Sig_Category_%i", (i+1)));
	  }

	  if (varNum ==1 || varNum ==2 || varNum ==4) 
	  {
		  efficiency[i]->SetTitle(TString::Format("Efficiency for %s for mva > %.1f for category: %d ;%s (GeV);Efficiency", varName.Data(), selMvaCut, i+1, varName.Data()));
		  h_Mtt[i][0]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", varNameReco.Data()));
	  }
	  else 
	  {
		  efficiency[i]->SetTitle(TString::Format("Efficiency for %s for mva > %.1f for category: %d ;%s;Efficiency", varNameReco.Data(), selMvaCut, i+1, varName.Data()));
		  h_Mtt[i][0]->GetXaxis()->SetTitle(TString::Format("%s", varNameReco.Data()));
	  }
	  efficiency[i]->SetLineColor(cols[i]);
	  efficiency[i]->Draw("SAME");
	  
	  efficiency[i]->Write(TString::Format("Mtt_Sig_Category_%i_%s", (i+1), varName.Data()));
	  h_Mtt[i][0]->Write(TString::Format("hMtt_Sig_Category_%i_%s", (i+1), varNameReco.Data()));
	  leg1->AddEntry(efficiency[i], categoryTitles[i], "l");
  }

  leg1->Draw();

  outf->Close();
}





