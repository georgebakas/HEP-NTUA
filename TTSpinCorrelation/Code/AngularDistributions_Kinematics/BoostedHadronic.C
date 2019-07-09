#include <iostream>
#include <TString.h>
#include <TH1F.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include "TMath.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TCanvas.h"
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

using namespace std;

TVector3 getBoostVector(TLorentzVector p4_1, TLorentzVector p4_2, TLorentzVector &p4CombinedVector);
std::vector<TString> listOfFiles;
std::vector<TString> histoNames;
std::vector<float> XSEC;
float LUMI = 59740;

void initFiles()
{
   //listOfFiles.push_back("/afs/cern.ch/user/g/gbakas/CMSSW_10_2_6/src/KKousour/TopAnalysis/prod/ttbar/flatTreeDCSV_Signal.root");
   listOfFiles.push_back("/eos/cms/store/user/gbakas/spinCorrelation/mc-2018/TT_noSC_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root");
   listOfFiles.push_back("/eos/cms/store/user/gbakas/spinCorrelation/mc-2018/TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root");
}

void initHistNames()
{
	histoNames.push_back("SC_sample");
	histoNames.push_back("noSC_sample");
}

void initXsections()
{
	XSEC.push_back(16.74);
	XSEC.push_back(16.74);
}
void initGlobals()
{
	initFiles();
	initHistNames();
	initXsections();
}


void BoostedHadronic()
{
  gStyle->SetOptStat(0);
  initGlobals();
  TFile *f;
  std::vector<int> quarkIds;
  std::vector<int> leptonIds;
  leptonIds.push_back(-11);
  leptonIds.push_back(11);
  leptonIds.push_back(-13);
  leptonIds.push_back(13); 
  
  for(int i=1; i<=4; i++) 
  {
	  quarkIds.push_back(-i);
	  quarkIds.push_back(i);
  }	

  std::vector<int>::iterator it;
  std::vector<int>::iterator it_lep;
  
  //Declare histograms
  TH1F *hW_angle[2];  
  TH1F *hTopW_angle[2];  
  TH1F *hDeltaPhi[2];  
  TH1F *hDeltaPhiW[2];  
  
  const int binsAngle=32;
  const int binsPhi = 30;
  for(int f=0; f<listOfFiles.size(); f++)
  {
	  
	hW_angle[f] = new TH1F(TString::Format("W+/- angle in ZMF %s", histoNames[f].Data()),TString::Format("W+/- angle in ZMF %s", histoNames[f].Data()), binsAngle, 0.8,1);  
	hTopW_angle[f] = new TH1F(TString::Format("top/W angle in ZMF %s", histoNames[f].Data()),TString::Format("top/W angle in ZMF %s", histoNames[f].Data()), binsAngle, 0.8,1);  
	hDeltaPhi[f] = new TH1F(TString::Format("dilepton #Delta#phi %s", histoNames[f].Data()),TString::Format("dilepton #Delta#phi %s", histoNames[f].Data()), binsPhi, 0,3);  
	hDeltaPhiW[f] = new TH1F(TString::Format("dilepton #Delta#phi W+/- %s", histoNames[f].Data()),TString::Format("dilepton #Delta#phi W+/- %s", histoNames[f].Data()), binsPhi, 0,3);  
	



	TFile *inf = TFile::Open(listOfFiles[f]);
	TTree *trInf = (TTree*)inf->Get("boosted/events");
    vector <float>*particleMass(0);
    vector <float>*particlePt(0);
    vector <float>*particleEta(0);
    vector <float>*particlePhi(0);
    //this is particle id
    vector <int>  *particleId(0);
    cout<<listOfFiles[f]<<endl;
    float genEvtWeight(0);
	
	
	//declare the branches
    trInf->SetBranchAddress("particleId"   ,&particleId); 
    trInf->SetBranchAddress("particleMass" ,&particleMass); 
    trInf->SetBranchAddress("particlePt"   ,&particlePt); 
    trInf->SetBranchAddress("particleEta"  ,&particleEta); 
	trInf->SetBranchAddress("genEvtWeight" ,&genEvtWeight);
    trInf->SetBranchAddress("particlePhi"  ,&particlePhi); 
	
	long NN = trInf->GetEntries();
	cout<<NN<<endl;
	
	
	//These are all in Lab frame
	TLorentzVector p4T(0,0,0,0), p4Tbar(0,0,0,0); //declaration of top 4vectors (top goes to Wplus and tbar goes to Wminus)
	TLorentzVector p4Wplus(0,0,0,0), p4Wminus(0,0,0,0); //declaration of W 4vectors
	TLorentzVector p4b(0,0,0,0), p4bbar(0,0,0,0); //declaration of b 4vectors (top goes to b and tbar goes to bbar)
	TLorentzVector p4Lepton[2]; //declaration of lepton 4vectors 
	
	//What comes out of W
	TLorentzVector p4Wplus_quark[2];
	TLorentzVector p4Wminus_quark[2];
	
	TLorentzVector p4T_ttbarZMF(0,0,0,0),p4Tbar_ttbarZMF(0,0,0,0); //these are tops boosted into ZMF
	
	int hadronic(0), semileptonic(0), dileptonic(0), topPairs(0);;
	int decade(0);
	for(int i =0; i<NN; i++)
	{
	  int allHadrons(0);
	  int allLeptons(0);
	  double progress = 10.0*i/(1.0*NN);
      int k = TMath::FloorNint(progress); 
      if (k > decade) 
       cout<<10*k<<" %"<<endl;
      decade = k;
	  
      trInf->GetEntry(i);
      //cout<<particleId->size()<<endl;
	//  cout<<"------------------------"<<endl;
	//  cout<<"size of particleId: "<<particleId->size()<<endl;
	  for(int l=0; l<particleId->size(); l++)
	  {
		  
		  //verify the particle Id 
		  if( (*particleId)[l] == 6)
		  {
			topPairs++;
			p4T.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);			
		  }
		  if((*particleId)[l] == -6)
		  {
			p4Tbar.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
		  }	
          if((*particleId)[l] == 24)
		  {
			p4Wplus.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);  
	      }
		  else if((*particleId)[l] == -24)
		  {
			p4Wminus.SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
		  }
		  
		  it = std::find(quarkIds.begin(),quarkIds.end(), fabs((*particleId)[l])); 
		  if( it != quarkIds.end()) 
		  {
			  //cout<<(*particleId)[l]<<endl;
			  allHadrons++; //means you found hadrons
			  
		  }
		  
		 //find dilepton events
      	 it_lep = std::find(leptonIds.begin(),leptonIds.end(), fabs((*particleId)[l])); 
         if(it_lep !=leptonIds.end())
		 {
			 allLeptons++;
			 if(allLeptons==1) p4Lepton[0].SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
			 else if(allLeptons==2) p4Lepton[1].SetPtEtaPhiM((*particlePt)[l], (*particleEta)[l], (*particlePhi)[l], (*particleMass)[l]);
		 }
		  
		  
		  //now find the angles 
		  
	  }
	  
	  //check dilepton
	  if(allLeptons ==2)
	  {
		  dileptonic++;
		  float deltaPhi = fabs(p4Lepton[0].Phi()-p4Lepton[1].Phi());
		  hDeltaPhi[f]->Fill(deltaPhi);
	  }
	  else if(allLeptons ==1)
	  {
		  semileptonic++;
	  }
	  
	  //after iteration work only with the hadronic events
	  if(allHadrons == 4) 
	  {
		  hadronic++;
		  //get boost vector
		  bool partonCuts = p4T.Pt() > 400 &&  p4Tbar.Pt() > 400 && fabs(p4T.Eta()) < 2.4 && fabs(p4Tbar.Eta()) <2.4;
		  if(partonCuts)
		  {
			  TLorentzVector p4TTbar;
			  TVector3 ttbarBoostVector = getBoostVector(p4T, p4Tbar, p4TTbar);
			  p4T_ttbarZMF = p4T;
			  p4Tbar_ttbarZMF = p4Tbar;
			  
			  //boost in ZMF
			  p4T_ttbarZMF.Boost(ttbarBoostVector);
	    	  p4Tbar_ttbarZMF.Boost(ttbarBoostVector);
			  
			  //boost the W's also
			  TLorentzVector p4Wminus_ZMF(0,0,0,0), p4Wplus_ZMF(0,0,0,0);
			  p4Wminus_ZMF = p4Wminus;
			  p4Wplus_ZMF  = p4Wplus;
			  p4Wminus_ZMF.Boost(ttbarBoostVector);
		  	  p4Wplus_ZMF.Boost(ttbarBoostVector);
			  
			  //find angles between Wplus and Wminus
			  Double_t thetaW = p4Wminus_ZMF.Angle(p4Wplus_ZMF.Vect());
			  float cosThetaW = fabs(TMath::Cos(thetaW));
			  hW_angle[f]->Fill(cosThetaW, genEvtWeight);
			  
			  //Angle between top and W in ZMF
			  Double_t thetaTopW[2];
			  thetaTopW[0] = fabs(TMath::Cos(p4T_ttbarZMF.Angle(p4Wplus_ZMF.Vect())));
			  thetaTopW[1] = fabs(TMath::Cos(p4Tbar_ttbarZMF.Angle(p4Wminus_ZMF.Vect())));
			  
			  hTopW_angle[f] -> Fill(thetaTopW[0]);
			  hTopW_angle[f] -> Fill(thetaTopW[1]);
			  
			  float deltaPhiW = fabs(p4Wminus.Phi()-p4Wplus.Phi());
			  hDeltaPhiW[f]->Fill(deltaPhiW);
		  }

	  }
		
		
	} 

	particleMass->clear();
	particlePt->clear();
	particlePhi->clear();
	particleEta->clear();
	particleId->clear();


	cout<<"number of Top pairs events: "<<topPairs<<endl;
    cout<<"Number of hadronic events: "<<hadronic<<endl;
    cout<<"Number of semileptonic events: "<<semileptonic<<endl;
    cout<<"Number of dileptonic events: "<<dileptonic<<endl;
    cout<<"Number of sum of all processes: "<<hadronic+semileptonic+dileptonic<<endl;
   } //end of files loop
	
	
	TLegend *leg = new TLegend(0.5,0.6,0.7,0.8);
	
	TCanvas *can_Wangle = new TCanvas("can_Wangle", "can_Wangle", 700,600);
	for(int f=0; f<listOfFiles.size(); f++)
	{
	   hW_angle[f]->Scale(1./hW_angle[f]->Integral(), "width");
	   leg->AddEntry(hTopW_angle[f], histoNames[f], "l");
	   if(f==0) hW_angle[f]->Draw("E");
	   else 
	   	{
	   		hW_angle[f]->SetLineColor(kRed);
	   		hW_angle[f]->Draw("Esame"); 
	   	}
	}
	leg->Draw();
	
	TCanvas *can_WTopAngle = new TCanvas("can_WTopAngle", "can_WTopAngle", 700,600);
	for(int f=0; f<listOfFiles.size(); f++)
	{
	   hTopW_angle[f]->Scale(1./hTopW_angle[f]->Integral(), "width");
	   if(f==0)hTopW_angle[f]->Draw("E");
	   else
	   {
	   	 hTopW_angle[f]->SetLineColor(kRed);
	   	 hTopW_angle[f]->Draw("Esame"); 
	   }
	}
	leg->Draw();
	
	TCanvas *can_deltaPhi = new TCanvas("can_deltaPhi", "can_deltaPhi", 700,600);
	for(int f=0; f<listOfFiles.size(); f++)
	{
	   hDeltaPhi[f]->Scale(1./hDeltaPhi[f]->Integral(), "width");
	   if(f==0)hDeltaPhi[f]->Draw("E");
	   else
	   {
	   	 hDeltaPhi[f]->SetLineColor(kRed);
	   	 hDeltaPhi[f]->Draw("Esame"); 
	   }
	}
	leg->Draw();
	
	TCanvas *can_deltaPhiW = new TCanvas("can_deltaPhiW", "can_deltaPhiW", 700,600);
	for(int f=0; f<listOfFiles.size(); f++)
	{
	   hDeltaPhiW[f]->Scale(1./hDeltaPhiW[f]->Integral(), "width");
	   if(f==0)hDeltaPhiW[f]->Draw("E");
	   else
	   {
	   	 hDeltaPhiW[f]->SetLineColor(kRed);
	   	 hDeltaPhiW[f]->Draw("Esame"); 
	   }
	}
	leg->Draw();
	
	//measure sensitivity between SC and noSC with chi-square method
	Double_t res[2][binsAngle];
	Double_t resDeltaPhi[2][binsPhi];
	cout<<"W+/- cos(Angle) chi2 value"<<endl;
	hW_angle[0]->Chi2Test(hW_angle[1],"WW P",res[0]);
	cout<<"W/top cos(Angle) chi2 value"<<endl;
	hTopW_angle[0]->Chi2Test(hTopW_angle[1],"WW P",res[1]);
	cout<<"W deltaPhi chi2 value"<<endl;
	hDeltaPhiW[0]->Chi2Test(hDeltaPhiW[1], "WW P CHI2/NDF", resDeltaPhi[1]);
    Double_t x[binsAngle];
	//Graph for Residuals
    for (Int_t i=1; i<binsAngle; i++) x[i]= 0+ i*(0.00625);
    TGraph *resgr = new TGraph(binsAngle,x,res[0]);
    resgr->GetXaxis()->SetRangeUser(4,16);
    resgr->GetYaxis()->SetRangeUser(-3.5,3.5);
    resgr->GetYaxis()->SetTitle("Normalized Residuals");
    resgr->SetMarkerStyle(21);
    resgr->SetMarkerColor(2);
    resgr->SetMarkerSize(.9);
    resgr->SetTitle("Normalized Residuals for hW_angle");
	
    TCanvas *canRes = new TCanvas("canRes", "canRes", 700, 600);
    resgr->Draw("APL");

    TGraph *resgr2 = new TGraph(binsAngle,x,res[1]);
    resgr2->GetXaxis()->SetRangeUser(4,16);
    resgr2->GetYaxis()->SetRangeUser(-3.5,3.5);
    resgr2->GetYaxis()->SetTitle("Normalized Residuals");
    resgr2->SetMarkerStyle(21);
    resgr2->SetMarkerColor(2);
    resgr2->SetMarkerSize(.9);
    resgr2->SetTitle("Normalized Residuals for hTopW_angle");
    TCanvas *canRes2 = new TCanvas("canRes2", "canRes2", 700, 600);
    resgr2->Draw("APL");
	
	Double_t xPhi[binsPhi];
	//residuals for deltaPhi of W
	for (Int_t i=1; i<binsPhi; i++) xPhi[i]= 0+ i*(0.1);
    TGraph *resPhi = new TGraph(binsPhi,xPhi,resDeltaPhi[0]);
    resPhi->GetXaxis()->SetRangeUser(4,16);
    resPhi->GetYaxis()->SetRangeUser(-3.5,3.5);
    resPhi->GetYaxis()->SetTitle("Normalized Residuals");
    resPhi->SetMarkerStyle(21);
    resPhi->SetMarkerColor(2);
    resPhi->SetMarkerSize(.9);
    resPhi->SetTitle("Normalized Residuals for deltaPhiW");
	
    TCanvas *canResPhi = new TCanvas("canResPhi", "canResPhi", 700, 600);
    resPhi->Draw("APL");
	
	
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
