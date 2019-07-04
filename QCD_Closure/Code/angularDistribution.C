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

TVector3 getBoostVector(TLorentzVector p4_1, TLorentzVector p4_2, TLorentzVector &p4CombinedVector);


void angularDistribution(TString file = "/eos/cms/store/user/ipapakri/ttbar/MC/Signal/TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8_legacy2016_deepAK8.root", 
						float selMvaCut=0.1, float floatBTag = 0.8838, bool isZprime= false,bool isParton=false, int ZprimeMass = 2000, TString width = "200" )
{
	
//TString TTbarFile = "/eos/cms/store/user/gbakas/ttbar/topTagger/April19/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_Copy.root"	
  gStyle->SetOptStat(0);
  TFile *inf     = TFile::Open(file);
  TTree *trIN    = (TTree*)inf->Get("boosted/events");
  //cout<<"here"<<endl;
  
  float XSEC = 832.;
  float NORM = ((TH1F*)inf->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
  float weight = XSEC/NORM;
  int nJets,nLeptons;
  float genEvtWeight;
  vector<float> *jetPt(0),*tau3(0),*tau2(0),*tau1(0);
  vector<float> *jetMassSub0(0), *jetMassSub1(0);
  vector<float> *jetMassSoftDrop(0);

  float mva(0);
  vector<float> *jetTtag(0);
  vector<bool> *bit(0);
  float mTTbarParton(0),mJJ(0), yTTbarParton(0), ptTTbarParton(0);
  int  category(0);
  //matching info 
  vector<float> *jetPhi(0), *jetEta(0);
  vector<int> *partonId(0), *partonMatchIdx(0);
  float yTopParton[2];
  vector<float> *partonEta(0), *partonPhi(0), *partonMatchDR(0),  *partonPt(0), *partonE(0), *partonMass(0);
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
  trIN->SetBranchAddress("jetBtagSub0"	  ,&jetBtagSub0);
  trIN->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1);
  trIN->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
  trIN->SetBranchAddress("mJJ"   		  ,&mJJ);
  trIN->SetBranchAddress("mva"	  		  ,&mva);
  trIN->SetBranchAddress("category"	  	  ,&category);
  trIN->SetBranchAddress("partonMatchIdx" ,&partonMatchIdx);
  trIN->SetBranchAddress("partonMatchDR"  ,&partonMatchDR);
  trIN->SetBranchAddress("jetTtagCategory",&jetTtag);
  //parton variables
  trIN->SetBranchAddress("mTTbarParton"   ,&mTTbarParton);
  trIN->SetBranchAddress("yTTbarParton"   ,&yTTbarParton);
  trIN->SetBranchAddress("ptTTbarParton"  ,&ptTTbarParton);
  trIN->SetBranchAddress("partonPt"		  ,&partonPt);
  trIN->SetBranchAddress("partonEta"	  ,&partonEta);
  trIN->SetBranchAddress("partonPhi" 	  ,&partonPhi);
  trIN->SetBranchAddress("partonE"	 	  ,&partonE);
  trIN->SetBranchAddress("partonMass"	  ,&partonMass);
  trIN->SetBranchAddress("yTopParton"	  ,&yTopParton);

  
  
  
  TLorentzVector p4T[2], p4TTbar, p4T_ZMF[2];

  int decade(0);
  int NN = trIN->GetEntries();

  //int NN = 10000;
  const int sizeBins = 1;
  //float BND[sizeBins+1] = {1000, 2000, 3000, 4000, 5000};
  //float BND[sizeBins+1] = {1000, 2500, 3500, 5000, 6000};
  float BND[sizeBins+1] = {1000, 6000};
  int counter =0;
  
  TH1F *h_mTTbarParton;
  if(isZprime)h_mTTbarParton  = new TH1F("mTTbarParton", "mTTbarParton histogram", 60, 1000,ZprimeMass+1000);
  else h_mTTbarParton = new TH1F("mTTbarParton", "mTTbarParton histogram", 60, 1000,6000);  
  
  const int chiSize =11;
  float BND_chi[chiSize+1] = {1,2,3,4,5,6,7,8,9,10,13,16};
  const int cosSize = 10;
  float BND_cos[cosSize+1] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
  TH1F *h_Chi_all = new TH1F("#chi dist", "#chi dist", chiSize, BND_chi);
  TH1F *h_Cos_all = new TH1F("#||{cos(#theta)} dist", "#||{cos(#theta)} dist", cosSize, BND_cos);
  
  //relation between x from costheta and x from exp(2y*)
  TH1F *hChiExp = new TH1F("#chi from e^{y^{*}} dist", "#chi from e^{y^{*}} dist", chiSize, BND_chi); //actual distribution for exp(2y*)
  TH1F *hChi_Diff = new TH1F("#chi_{cos(#theta)} - #chi_{e^{|2y^{*}|}} dist", "#chi_{cos(#theta)} - #chi_{e^{|2y^{*}|}} dist", 10, 0,1); //this is the x_fromCosTheta - x_fromExp
  
  TH1F *hChiBoost = new TH1F("#chi boost dist", "#chi boost dist", chiSize, BND_chi);
  TH1F *hCosBoost = new TH1F("#||{cos(#theta)} boost dist", "#||{cos(#theta)} boost dist", cosSize, BND_cos);
  
  TH1F *hChiCR = new TH1F("#chi CR dist", "#chi CR dist", chiSize, BND_chi);
  TH1F *hCosCR = new TH1F("#||{cos(#theta)} CR dist", "#||{cos(#theta)} CR dist", cosSize, BND_cos);
  
  TH1F *hAngularDist[sizeBins], *hChi[sizeBins];
  std::vector<TString> massLimits = {"1000-2500","2500-3500", "3500-5000", "5000-Inf"};
  
  for(int i=0; i<=(sizeBins-1); i++)
  {
	  TString temp = massLimits[i];
	  hAngularDist[i] = new TH1F(TString::Format("#||{cos(#theta)} Distribution for mass limit: %s (GeV)", temp.Data()),TString::Format("#||{cos(#theta)} Distribution for mass limit: %s (GeV)", temp.Data()), 20, 0,1);
	  hChi[i] = new TH1F(TString::Format("#chi Distribution for mass limit: %s (GeV)", temp.Data()),TString::Format("#chi Distribution for mass limit: %s (GeV)", temp.Data()), 20, 0,16);
  }
  //for matching
  std::vector<int> *jetMatchedIndexes = new std::vector<int>(0);
  std::vector<float> *jetMatchedDr = new std::vector<float>(0);
  std::vector<float> *eta_ = new std::vector<float>(0);
  std::vector<float> *phi_ = new std::vector<float>(0);
  std::vector<float> *mass_ = new std::vector<float>(0);
  std::vector<float> *pt_ = new std::vector<float>(0);
  std::vector<float> *jetBtagSub0_ = new std::vector<float>(0);
  std::vector<float> *jetBtagSub1_ = new std::vector<float>(0);
  std::vector<float> *jetTtag_ = new std::vector<float>(0);
  
  std::vector<float> *partonPt_ = new std::vector<float>(0);
  std::vector<float> *partonEta_ = new std::vector<float>(0);
  std::vector<float> *partonMass_ = new std::vector<float>(0);
  std::vector<float> *partonPhi_ = new std::vector<float>(0);
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
	eta_->clear();
	mass_->clear();
	pt_->clear();
	phi_->clear();
	jetBtagSub0_->clear();
	jetBtagSub1_->clear();
	jetTtag_->clear();
	
	partonPt_->clear();
	partonMass_->clear();
	partonPhi_->clear();
	partonEta_->clear();
	
	int isMatched=0;
	bool recoCuts, partonCuts, btagging, topTagger, massCut, btaggingReverted;	
	if (nJets >1)
	{
		
		//----------------------MATCHING------------------------------------------------------
		
		for(int ijet =0; ijet<nJets; ijet++)
		{
			//cout<<"ok"<<endl;
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
				
				//cout<<"dRmin[0]: "<<dRmin<<endl;
				for(int k=1; k<jetMatchedIndexes->size(); k++)
				{
					//cout<<"jetMatchedIndexes at k = "<<k<<" is: "<<(*jetMatchedIndexes)[k]<<endl;
					//cout<<"jetMatchedDr at k =  "<<k<<" is: "<<(*jetMatchedDr)[k]<<endl;
					if((*jetMatchedDr)[k] < dRmin)
					{
						dRmin = (*jetMatchedDr)[k];
						indexMin = (*jetMatchedIndexes)[k];
					}
				//cout<<"dRmin is: "<<dRmin<<endl;
				}
				if(dRmin < 0.4)
				{
					isMatched++;
					//cout<<"int isMatched: "<<isMatched<<endl;
					jetDr_ = dRmin;
					pt_->push_back((*jetPt)[(*partonMatchIdx)[indexMin]]);
					mass_->push_back((*jetMassSoftDrop)[(*partonMatchIdx)[indexMin]]);
					eta_->push_back((*jetEta)[(*partonMatchIdx)[indexMin]]);
					phi_->push_back( (*jetPhi)[(*partonMatchIdx)[indexMin]]);
					jetBtagSub0_->push_back( (*jetBtagSub0)[(*partonMatchIdx)[indexMin]]);
					jetBtagSub1_->push_back( (*jetBtagSub1)[(*partonMatchIdx)[indexMin]]);
					jetTtag_->push_back( (*jetTtag)[(*partonMatchIdx)[indexMin]]);
					
					partonPt_->push_back( (*partonPt)[indexMin]);
					partonMass_->push_back( (*partonMass)[indexMin]);
					partonPhi_->push_back( (*partonPhi)[indexMin]);
					partonEta_->push_back( (*partonEta)[indexMin]);
					
				}
					
		   }
		   
		 }//---------------------------end of MATCHING---------------------------------------------------------
		 
		if(isMatched >1)	 
		{
				recoCuts   = fabs((*eta_)[0]) < 2.4 && fabs((*eta_)[1]) <2.4 && (*pt_)[0] > 400 && (*pt_)[1] > 400 && nLeptons==0 && mJJ > 1000;
				partonCuts = (*partonPt_)[0] > 400 && (*partonPt_)[1] > 400 && fabs((*partonEta_)[0]) < 2.4 && fabs((*partonEta_)[1]) < 2.4 &&  mTTbarParton > 1000;
				btagging   = ((*jetBtagSub0_)[0] > floatBTag || (*jetBtagSub1_)[0] > floatBTag) && ((*jetBtagSub0_)[1] > floatBTag || (*jetBtagSub1_)[1] > floatBTag);
				topTagger  = (*jetTtag_)[0] > selMvaCut && (*jetTtag_)[1] > selMvaCut;
				massCut    = (*mass_)[0] > 120 && (*mass_)[0] < 220 && (*mass_)[1] > 120 && (*mass_)[1] < 220;
				//btaggingReverted = category==0;
				btaggingReverted = ((*jetBtagSub0_)[0] < floatBTag && (*jetBtagSub1_)[0] < floatBTag) && ((*jetBtagSub0_)[1] < floatBTag && (*jetBtagSub1_)[1] < floatBTag);

				//cout<<"ok"<<endl;
					
					//split the mttbar phase space into regions and for each region I will calculate the thetas dists
					bool currentSelection;
					if(isParton)
					{
						p4T[0].SetPtEtaPhiM((*partonPt_)[0], (*partonEta_)[0], (*partonPhi_)[0], (*partonMass_)[0]);
						p4T[1].SetPtEtaPhiM((*partonPt_)[1], (*partonEta_)[1], (*partonPhi_)[1], (*partonMass_)[1]);
						currentSelection = partonCuts;
					}
					else
					{
						p4T[0].SetPtEtaPhiM((*pt_)[0], (*eta_)[0], (*phi_)[0], (*mass_)[0]);
						p4T[1].SetPtEtaPhiM((*pt_)[1], (*eta_)[1], (*phi_)[1], (*mass_)[1]);
						currentSelection = recoCuts && topTagger && btagging;
					}
					
					TVector3 ttbarBoostVector = getBoostVector(p4T[0], p4T[1], p4TTbar);
					
					p4T_ZMF[0].SetPtEtaPhiM(p4T[0].Pt(), p4T[0].Eta(), p4T[0].Phi(), p4T[0].M());
					p4T_ZMF[1].SetPtEtaPhiM(p4T[1].Pt(), p4T[1].Eta(), p4T[1].Phi(), p4T[1].M());
					p4T_ZMF[0].Boost(ttbarBoostVector);
					p4T_ZMF[1].Boost(ttbarBoostVector);
							
					//cout<<"-------------------------------------"<<endl;		
					//cout<< p4T_ZMF[0].Pt()<<endl;
					//cout<< p4T_ZMF[1].Pt()<<endl;
					float chi0(0), chi1(0);
					chi0 = (1 + fabs(TMath::Cos(p4T_ZMF[0].Theta()))) / ( 1 - fabs(TMath::Cos(p4T_ZMF[0].Theta())));
					//chi1 = (1 + fabs(TMath::Cos(p4T_ZMF[1].Theta()))) / ( 1 - fabs(TMath::Cos(p4T_ZMF[1].Theta())));
					//chi0 = TMath::Exp(fabs(p4T_ZMF[0].Rapidity() - p4T_ZMF[1].Rapidity()));
					float yStarExp  = TMath::Exp(fabs(p4T_ZMF[0].Rapidity() - p4T_ZMF[1].Rapidity()));
					
					if(currentSelection)
					{
						//h_Chi_all->Fill(chi0);
						h_Chi_all->Fill(yStarExp);
	
						//float testYstar = TMath::Exp(fabs(p4T[0].Rapidity() - p4T[1].Rapidity()));
						hChiExp->Fill(yStarExp);

						h_Cos_all->Fill(fabs(TMath::Cos(p4T_ZMF[0].Theta())));
						//h_Cos_all->Fill(fabs(TMath::Cos(p4T_ZMF[1].Theta())));
						
						hChi_Diff->Fill(chi0-yStarExp);
						
						if(chi0-yStarExp < 0) cout<<"Found less that 0"<<endl;
						
						if(fabs((0.5)*(p4T_ZMF[0].Rapidity() + p4T_ZMF[1].Rapidity())) < 1.1 )
						{
							hChiBoost->Fill(chi0);
							hCosBoost->Fill(fabs(TMath::Cos(p4T_ZMF[0].Theta())));
						}	
				    }//----end of cuts
				
				   if(recoCuts && topTagger && btaggingReverted)
				   {
				 	  hChiCR->Fill(yStarExp);
				  	  hCosCR->Fill(fabs(TMath::Cos(p4T_ZMF[0].Theta())));
				   }
				   if(recoCuts && partonCuts && topTagger && btagging)
				   {
					 h_mTTbarParton->Fill(mTTbarParton);
					 bool found =false;
					 int imass = 0;
					 while(!found && imass<sizeBins)
					 {
					 	//cout<<imass<<endl;
						//cout<<"------------------------"<<endl;
					 	 if(mTTbarParton >= BND[imass] && mTTbarParton < BND[imass+1] )
						 {
							 found =true;
							 hAngularDist[imass]->Fill(fabs(TMath::Cos(p4T_ZMF[0].Theta())));
							 //hAngularDist[imass]->Fill(fabs(TMath::Cos(p4T_ZMF[1].Theta())));
							 hChi[imass]->Fill(yStarExp);
							 //hChi[imass]->Fill(chi9);
							 //hChi[imass]->Fill(chi1);
						 }
						 imass++;
					 }
				   }
				
		}//----end of isMatched	
	}//----end of nJets
  }	//---end of event loop
  
  TCanvas *can[sizeBins];
  for(int i =0; i<sizeBins; i++)
  {
	can[i] = new TCanvas(TString::Format("can %d", (i+1)), TString::Format("can %d", (i+1)), 900, 600);
	if(isParton) hAngularDist[i]->GetXaxis()->SetTitle("cos(#theta)");
	else hAngularDist[i]->GetXaxis()->SetTitle(TString::Format("cos(#theta) %s","reco"));
	hAngularDist[i]->Draw();
  }
  
  TCanvas *canChi[sizeBins];
  for(int i =0; i<sizeBins; i++)
  {
	canChi[i] = new TCanvas(TString::Format("can chi %d", (i+1)), TString::Format("can chi %d", (i+1)), 900, 600);
	if(isParton) hChi[i]->GetXaxis()->SetTitle("#chi");
	else hChi[i]->GetXaxis()->SetTitle(TString::Format("#chi %s", "reco"));
	hChi[i]->Draw();
  }
  
  TCanvas *can_mTTbarParton = new TCanvas("mTTbarParton can", "mTTbarParton can", 900, 600);
  h_mTTbarParton->GetXaxis()->SetTitle("mTTbarParton (GeV)");
  h_mTTbarParton->Draw();
  
  TCanvas *can_chiDiff = new TCanvas("chi Diff","chi Diff" ,900, 600);
  hChi_Diff->GetXaxis()->SetTitle("#chi");
  hChi_Diff->Draw();
  
  TCanvas *can_cosAll = new TCanvas("|cos(#theta)|can", "#|cos(#theta)| exp can", 900, 600);
  h_Cos_all->GetXaxis()->SetTitle("|cos(#theta^{*})|");
  h_Cos_all->Scale(1./h_Cos_all->Integral(), "width");
  h_Cos_all->Draw();
  
  TCanvas *can_chiAll = new TCanvas("#chi All can", "#chi all can", 900, 600);
  h_Chi_all->GetXaxis()->SetTitle("#chi");
  h_Chi_all->Scale(1./h_Chi_all->Integral(), "width");
  h_Chi_all->Draw();
  //can_chiAll->BuildLegend();

  TCanvas *can_chiCR = new TCanvas("#chi hChiCR can", "#chi hChiCR can", 900, 600);
  hChiCR->SetTitle("#chi top contamination in CR ");
  hChiCR->GetXaxis()->SetTitle("#chi");
  //hChiCR->Scale(1./hChiCR->Integral(), "width");
  hChiCR->Draw();
  
  TCanvas *can_cosCR = new TCanvas("#chi =hCosCR can", "#chi hCosCR can", 900, 600);
  hCosCR->SetTitle("cos(#theta) top contamination in CR ");
  hCosCR->GetXaxis()->SetTitle("cos(#theta)");
  //hCosCR->Scale(1./hCosCR->Integral(), "width");
  hCosCR->Draw();

 
  auto c3 = new TCanvas("#chi No yBoost/yBoost", "#chi No yBoost/yBoost", 900,600);
  TLegend *leg_chi = new TLegend(0.6,0.7,0.8,0.9);
  leg_chi->AddEntry(h_Cos_all, "#chi no |y_{boost}| cut", "l"); 
  leg_chi->AddEntry(hCosBoost, "#chi with |y_{boost}| < 1.19", "l");  
  //h_Chi_all->GetXaxis()->SetTitle("#chi");
  //h_Chi_all->Scale(1./h_Chi_all->Integral());
  hChiBoost->Scale(1./hChiBoost->Integral(), "width");
  hChiBoost->SetLineColor(kRed);
  auto rp_chi = new TRatioPlot(h_Chi_all, hChiBoost);
  c3->SetTicks(0,1);
  rp_chi->Draw();
  leg_chi->Draw();
  c3->Update();  
  
  auto c2 = new TCanvas("Cos No yBoost/yBoost", "Cos No yBoost/yBoost", 900,600);
  TLegend *leg_cos = new TLegend(0.6,0.7,0.8,0.9);
  leg_cos->AddEntry(h_Cos_all, "cos(#theta^{*}) no |y_{boost}| cut", "l"); 
  leg_cos->AddEntry(hCosBoost, "cos(#theta^{*}) with |y_{boost}| < 1.19", "l"); 
  //h_Cos_all->GetXaxis()->SetTitle("#||{cos(#theta)}");
  //h_Cos_all->Scale(1./h_Cos_all->Integral());
  hCosBoost->Scale(1./hCosBoost->Integral(), "width");
  hCosBoost->SetLineColor(kRed);
  auto rp_cos = new TRatioPlot(h_Cos_all, hCosBoost);
  c2->SetTicks(0,1);
  rp_cos->Draw();
  c2->Update();
  leg_cos->Draw();

  
  auto c1 = new TCanvas("hChi over hChiExp", "hChi over hChiExp", 900,600);
  TLegend *leg_rp = new TLegend(0.6,0.7,0.8,0.9);
  leg_rp->AddEntry(h_Chi_all,"#chi distibution using angle #theta", "l");
  leg_rp->AddEntry(hChiExp  ,"#chi distibution using e^{|2y^{*}|}", "l");
  //h_Chi_all->GetXaxis()->SetTitle("#chi");
  //h_Chi_all->Scale(1./h_Chi_all->Integral());
  hChiExp->Scale(1./hChiExp->Integral(),"width");
  hChiExp->SetLineColor(kRed);
  auto rp = new TRatioPlot(h_Chi_all, hChiExp);
  c1->SetTicks(0,1);
  rp->Draw();
  leg_rp->Draw();
  c1->Update();
    
  TString tempMass;
  if(ZprimeMass == 2000) tempMass = "2TeV";
  else if(ZprimeMass == 3000) tempMass = "3TeV";
  else if(ZprimeMass == 4000) tempMass = "4TeV";
  else if(ZprimeMass == 2500) tempMass = "2.5TeV";
  else if(ZprimeMass == 5000) tempMass = "5TeV";
    
  TFile *outf;
  TString recoParton ="";
  if(isParton) recoParton = "Parton";
  else recoParton = "Reco";
  if(isZprime)  outf = new TFile(TString::Format("Output_M%s_%s_Chi.root", tempMass.Data(), recoParton.Data()), "UPDATE");
  else outf = new TFile(TString::Format("Output_TT_QCD_%s_Chi_%0.1f.root", recoParton.Data(),selMvaCut), "RECREATE");
  
  if (isZprime) 
  {
	  h_mTTbarParton->Write(TString::Format("h_mTTbarParton_M%s_W%s",tempMass.Data(), width.Data()));
	  h_Cos_all->Write(TString::Format("hCos_M%s_W%s",tempMass.Data(), width.Data()));
	  h_Chi_all->Write(TString::Format("hChi_M%s_W%s",tempMass.Data(), width.Data()));
	  hChiBoost->Write(TString::Format("hChiBoost_M%s_W%s",tempMass.Data(), width.Data()));
	  hCosBoost->Write(TString::Format("hCosBoost_M%s_W%s",tempMass.Data(), width.Data()));
	  hChiExp->Write(TString::Format("hChiExp_M%s_W%s",tempMass.Data(), width.Data()));
  }
  else 
  {
	  h_mTTbarParton->Write("h_mTTbarParton_TT"); 
	  h_Cos_all->Write("hCos_TT");
	  h_Chi_all->Write("hChi_TT");
	  hChiBoost->Write("hChiBoost_TT");
	  hCosBoost->Write("hCosBoost_TT");
	  hChiExp->Write("hChiExp_TT");
	  hCosCR->Write("hCosCR_TT");
	  hChiCR->Write("hChiCR_TT");
  }
  
  
  
  for(int i =0; i<sizeBins; i++)
  {
	  TString temp = massLimits[i];
	  if(isZprime) 
	  {
		  hAngularDist[i]->Write(TString::Format("abs_cosTheta_%s_M%s_W%s",temp.Data(), tempMass.Data(), width.Data()));
		  hChi[i]->Write(TString::Format("chi_%s_M%s_W%s",temp.Data(), tempMass.Data(), width.Data()));
	  }
	  else 
	  {
		  hAngularDist[i]->Write(TString::Format("abs_cosTheta_%s_TT",temp.Data() ));
		  hChi[i]->Write(TString::Format("chi_%s_TT", temp.Data()));
	  }
	  
  }
  
  
  
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

