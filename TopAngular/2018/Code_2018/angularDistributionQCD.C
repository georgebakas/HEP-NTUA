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

std::vector<TString> listOfFiles;
std::vector<float> XSEC;
std::vector<TString> histoNames;
float LUMI = 59740;
TString eosPath;
float deepCSVFloat = 0.4184;
void initFileNames()
{
  eosPath = "/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2018/Bkg/";
  listOfFiles.push_back("QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8.root");
  listOfFiles.push_back("QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8.root");
  listOfFiles.push_back("QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8.root");
  listOfFiles.push_back("QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8.root");
  listOfFiles.push_back("QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8.root");
  listOfFiles.push_back("QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8.root");

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

void angularDistributionQCD(float selMvaCut=0.1, float floatBTag = 0.8838, float isDeepCSV = true)
{

  initGlobals();	
//TString TTbarFile = "/eos/cms/store/user/gbakas/ttbar/topTagger/April19/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_Copy.root"	
  gStyle->SetOptStat(0);
  
  
  //cout<<"here"<<endl;
  
  //parton variables not existing in QCD samples
  //int NN = 10000;
  const int sizeBins = 4;
  //float BND[sizeBins+1] = {1000, 2000, 3000, 4000, 5000};
  //float BND[sizeBins+1] = {1000, 2500, 3500, 5000, 6000};
  float BND[sizeBins+1] = {1000,1600,2200,3200,6000};
  int counter =0;
    
  const int chiSize =11;
  float BND_chi[chiSize+1] = {1,2,3,4,5,6,7,8,9,10,13,16};
  const int cosSize = 10;
  float BND_cos[cosSize+1] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
  
  int fileSize = listOfFiles.size();
  TH1F *h_Chi_all[fileSize],*h_Cos_all[fileSize], *hChiRevertBtag[fileSize], *hCosRevertBtag[fileSize], *hChi1btag[fileSize], *hCos1btag[fileSize]; 
  TFile *inf;
  vector<float> weights(0);
  
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
  float mTTbarParton(0),mJJ(0), yTTbarParton(0), ptTTbarParton(0);
  int  category(0);
  //matching info 
  vector<float> *jetPhi(0), *jetEta(0);
  vector<int> *partonId(0), *partonMatchIdx(0);
  
  vector<float> *partonEta(0), *partonPhi(0), *partonMatchDR(0),  *partonPt(0), *partonE(0), *partonMass(0);
  std::vector<int> *addedIndexes = new std::vector<int>(0);
  std::vector<float> *jetBtagSub0(0), *jetBtagSub1(0);
  std::vector<float> *jetBtagSub0DCSVbb(0), *jetBtagSub1DCSVbb(0);
  std::vector<float> *jetBtagSub0DCSVbbb(0), *jetBtagSub1DCSVbbb(0);       
  
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
  trIN->SetBranchAddress("jetTtagCategory",&jetTtag);
   //deepCSV
  trIN->SetBranchAddress("jetBtagSub0DCSVbb" ,&jetBtagSub0DCSVbb);
  trIN->SetBranchAddress("jetBtagSub1DCSVbb" ,&jetBtagSub1DCSVbb);
  trIN->SetBranchAddress("jetBtagSub0DCSVbbb",&jetBtagSub0DCSVbbb);
  trIN->SetBranchAddress("jetBtagSub1DCSVbbb",&jetBtagSub1DCSVbbb);
  
  TLorentzVector p4T[2], p4TTbar, p4T_ZMF[2];
  
  h_Chi_all[f] = new TH1F(TString::Format("%s_#chi_SR", histoNames[f].Data()), TString::Format("%s_#chi_SR", histoNames[f].Data()), chiSize, BND_chi);
  h_Cos_all[f] = new TH1F(TString::Format("%s_#||{cos(#theta)}_SR", histoNames[f].Data()),TString::Format("%s_#||{cos(#theta)}_SR", histoNames[f].Data()), cosSize, BND_cos);
  

  hChiRevertBtag[f] = new TH1F(TString::Format("%s_#chi_CR", histoNames[f].Data()), TString::Format("%s_#chi_CR", histoNames[f].Data()), chiSize, BND_chi);
  hCosRevertBtag[f] = new TH1F(TString::Format("%s_#||{cos(#theta)}_CR", histoNames[f].Data()),TString::Format("%s_#||{cos(#theta)}_CR",histoNames[f].Data()),cosSize, BND_cos);
  
  hChi1btag[f] = new TH1F(TString::Format("%s_#chi_1btag", histoNames[f].Data()), TString::Format("%s_#chi_1btag", histoNames[f].Data()), chiSize, BND_chi);
  hCos1btag[f] = new TH1F(TString::Format("%s_#||{cos(#theta)}_1btag", histoNames[f].Data()),TString::Format("%s_#||{cos(#theta)}_1btag",histoNames[f].Data()),cosSize, BND_cos);
  
  //TH1F *hAngularDist[sizeBins], *hChi[sizeBins];
  //std::vector<TString> massLimits = {"1000-2500","2500-3500", "3500-5000", "5000-Inf"};
  std::vector<TString> massLimits = {"1000-1600","1600-2200", "2200-3200","3200-6000"};
  /*
  for(int i=0; i<=(sizeBins-1); i++)
  {
	  TString temp = massLimits[i];
	  hAngularDist[i] = new TH1F(TString::Format("#||{cos(#theta)} Distribution for mass limit: %s (GeV)", temp.Data()),TString::Format("#||{cos(#theta)} Distribution for mass limit: %s (GeV)", temp.Data()), 20, 0,1);
	  hChi[i] = new TH1F(TString::Format("#chi Distribution for mass limit: %s (GeV)", temp.Data()),TString::Format("#chi Distribution for mass limit: %s (GeV)", temp.Data()), 20, 0,16);
  }*/
  //for matching
  std::vector<int> *jetMatchedIndexes = new std::vector<int>(0);
  std::vector<float> *jetMatchedDr = new std::vector<float>(0);
  std::vector<float> *eta_ = new std::vector<float>(0);
  std::vector<float> *phi_ = new std::vector<float>(0);
  std::vector<float> *mass_ = new std::vector<float>(0);
  std::vector<float> *pt_ = new std::vector<float>(0);
  
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

	int isMatched=0;
	bool recoCuts, btagging, topTagger, massCut, revertBtag, btag1;	
	bool CSVv2Cut, deepCSV, btag1CSVv2, btag1DeepCSV, revertBtagCSVv2, revertBtagDeepCSV;

	if (nJets >1)
	{		
		  float dCSVScoreSub0[2], dCSVScoreSub1[2];
          dCSVScoreSub0[0] = (*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0];
          dCSVScoreSub0[1] = (*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1];
          dCSVScoreSub1[0] = (*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0];
          dCSVScoreSub1[1] = (*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1];
		  bool CSVv2Cut, deepCSV, btag1CSVv2, btag1DeepCSV, revertBtagCSVv2, revertBtagDeepCSV;
		  massCut    = (*jetMassSoftDrop)[0] > 120 && (*jetMassSoftDrop)[0] < 220 && (*jetMassSoftDrop)[1] > 120 && (*jetMassSoftDrop)[1] < 220;		  
          recoCuts   = fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1]) <2.4 && (*jetPt)[0] > 400 && (*jetPt)[1] > 400 &&  mJJ > 1000 && nLeptons==0 && (*bit)[5] && massCut;
          topTagger = (*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut;
          //2 btag category with csvv2 and deepCSV
          CSVv2Cut   = ((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag) && ((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag);

          deepCSV    = (((*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0])> deepCSVFloat) &&
                                         (((*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1])> deepCSVFloat);
          //1 btag category with csvv2 and deepCSV
          btag1CSVv2 = (((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag) && ((*jetBtagSub0)[1] < floatBTag || (*jetBtagSub1)[1] < floatBTag)) ||
                                        (((*jetBtagSub0)[0] < floatBTag || (*jetBtagSub1)[0] < floatBTag) && ((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag));

          btag1DeepCSV	= ((dCSVScoreSub0[0] > deepCSVFloat || dCSVScoreSub1[0] > deepCSVFloat) && (dCSVScoreSub0[1] < deepCSVFloat && dCSVScoreSub1[1] < deepCSVFloat)) ||
					  ((dCSVScoreSub0[0] < deepCSVFloat && dCSVScoreSub1[0] < deepCSVFloat) && (dCSVScoreSub0[1] > deepCSVFloat || dCSVScoreSub1[1] > deepCSVFloat));

          //0 btag category with csvv2 and deepCSV
		 revertBtagCSVv2 = ((*jetBtagSub0)[0] < floatBTag && (*jetBtagSub1)[0] < floatBTag) && ((*jetBtagSub0)[1] < floatBTag && (*jetBtagSub1)[1] < floatBTag);
         revertBtagDeepCSV = (dCSVScoreSub0[0] < deepCSVFloat &&  dCSVScoreSub1[0] < deepCSVFloat) && (dCSVScoreSub0[1] < deepCSVFloat && dCSVScoreSub1[1] < deepCSVFloat);
		
		if(isDeepCSV)
		{
          btagging = deepCSV;
          revertBtag = revertBtagDeepCSV;
          btag1 = btag1DeepCSV;
        }
        else
        {
          btagging = CSVv2Cut;
          revertBtag = revertBtagCSVv2;
          btag1 = btag1CSVv2;
		}

		
		p4T[0].SetPtEtaPhiM((*jetPt)[0], (*jetEta)[0], (*jetPhi)[0], (*jetMassSoftDrop)[0]);
		p4T[1].SetPtEtaPhiM((*jetPt)[1], (*jetEta)[1], (*jetPhi)[1], (*jetMassSoftDrop)[1]);
		
		TVector3 ttbarBoostVector = getBoostVector(p4T[0], p4T[1], p4TTbar);		
		p4T_ZMF[0].SetPtEtaPhiM(p4T[0].Pt(), p4T[0].Eta(), p4T[0].Phi(), p4T[0].M());
		p4T_ZMF[1].SetPtEtaPhiM(p4T[1].Pt(), p4T[1].Eta(), p4T[1].Phi(), p4T[1].M());
		p4T_ZMF[0].Boost(ttbarBoostVector);
		p4T_ZMF[1].Boost(ttbarBoostVector);
							
		//cout<<"-------------------------------------"<<endl;		
		//cout<< p4T_ZMF[0].Pt()<<endl;
		//cout<< p4T_ZMF[1].Pt()<<endl;
		float chi0(0), chi1(0);
		//chi0 = (1 + fabs(TMath::Cos(p4T_ZMF[0].Theta()))) / ( 1 - fabs(TMath::Cos(p4T_ZMF[0].Theta())));
		//chi1 = (1 + fabs(TMath::Cos(p4T_ZMF[1].Theta()))) / ( 1 - fabs(TMath::Cos(p4T_ZMF[1].Theta())));	
		chi0 = TMath::Exp(fabs(p4T_ZMF[0].Rapidity() - p4T_ZMF[1].Rapidity()));		//take it from exp	
		if(recoCuts && topTagger && btagging)
		{
		  //cout<<"SR chi: "<<chi0<<endl;
		  h_Chi_all[f]->Fill(chi0,genEvtWeight);
		  h_Cos_all[f]->Fill(fabs(TMath::Cos(p4T_ZMF[0].Theta())),genEvtWeight);
		}//----end of Signal Region cuts
			
		if(recoCuts && topTagger && revertBtag)
		{
		 //cout<<"CR chi: "<<chi0<<endl;
		  hChiRevertBtag[f]->Fill(chi0,genEvtWeight);
		  hCosRevertBtag[f]->Fill(fabs(TMath::Cos(p4T_ZMF[0].Theta())),genEvtWeight);
		} //---end of control region cuts
		
		if(recoCuts && topTagger && btag1)
		{
		 //cout<<"CR chi: "<<chi0<<endl;
		  hChi1btag[f]->Fill(chi0,genEvtWeight);
		  hCos1btag[f]->Fill(fabs(TMath::Cos(p4T_ZMF[0].Theta())),genEvtWeight);
		} //1 btag region
					

	}//----end of nJets
  }	//---end of event loop
  }//----end of fileSize loop 


    
  for(int i=0; i< fileSize; i++)
  {
	h_Cos_all[i]->GetXaxis()->SetTitle("|cos(#theta^{*})|");  
	h_Chi_all[i]->GetXaxis()->SetTitle("#chi");
	
	hChiRevertBtag[i]->GetXaxis()->SetTitle("#chi");
    hChiRevertBtag[i]->SetLineColor(kRed);
	
	hCosRevertBtag[i]->GetXaxis()->SetTitle("|cos(#theta^{*})|");
    hCosRevertBtag[i]->SetLineColor(kRed);
	
	hChi1btag[i]->GetXaxis()->SetTitle("#chi");
    hChi1btag[i]->SetLineColor(kMagenta);
	
	hCos1btag[i]->GetXaxis()->SetTitle("|cos(#theta^{*})|");
    hCos1btag[i]->SetLineColor(kMagenta);
	
	
	h_Chi_all[i]->Scale(weights[i]*LUMI);
    hChiRevertBtag[i]->Scale(weights[i]*LUMI);
   	hChi1btag[i]->Scale(weights[i]*LUMI);

	h_Cos_all[i]->Scale(weights[i]*LUMI);
    hCosRevertBtag[i]->Scale(weights[i]*LUMI);
    hCos1btag[i]->Scale(weights[i]*LUMI);
  
	if(i !=0) 
	{	
     h_Chi_all[0]->Add(h_Chi_all[i]);
	 h_Cos_all[0]->Add(h_Cos_all[i]);
	 
	 hChiRevertBtag[0]->Add(hChiRevertBtag[i]);
	 hCosRevertBtag[0]->Add(hCosRevertBtag[i]);
	 
	 hCos1btag[0]->Add(hCos1btag[i]);
	 hChi1btag[0]->Add(hChi1btag[i]);
	}
  }

  TFile *outf;
  if(isDeepCSV) outf = new TFile(TString::Format("Closure_QCDBkg_Chi_%0.1f_deepCSV.root",selMvaCut), "RECREATE");
  else outf = new TFile(TString::Format("Closure_QCDBkg_Chi_%0.1f.root",selMvaCut), "RECREATE");
  outf->cd();
  h_Cos_all[0]->Write("hCos_QCD_SR");
  h_Chi_all[0]->Write("hChi_QCD_SR");
  hChiRevertBtag[0]->Write("hChi_QCD_CR");
  hCosRevertBtag[0]->Write("hCos_QCD_CR");
  hCos1btag[0]->Write("hCos_QCD_1btag");
  hChi1btag[0]->Write("hChi_QCD_1btag");
  
  listOfFiles.clear();
  XSEC.clear();
  histoNames.clear();
  
  
  /*
  TCanvas *canChiBoth = new TCanvas("#chi qcd bkg closure", "#chi qcd bkg closure", 700, 600);
  
  h_Chi_all[0]->Draw();
  hChiRevertBtag[0]->Draw("same");
  canChiBoth->Write("ChiCanvas");
  
  TCanvas *canCosBoth = new TCanvas("cos qcd bkg closure", "cos qcd bkg closure", 700, 600);
  
  h_Cos_all[0]->Draw();
  hCosRevertBtag[0]->Draw("same");
  canCosBoth->Write("CosCanvas");
  
	TString tempMass;
	if(ZprimeMass == 2000) tempMass = "2TeV";
	else if(ZprimeMass == 3000) tempMass = "3TeV";
	else if(ZprimeMass == 4000) tempMass = "4TeV";
	else if(ZprimeMass == 2500) tempMass = "2.5TeV";
	else if(ZprimeMass == 5000) tempMass = "5TeV";
  */
  
 


  /*
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
  */
  
  
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

