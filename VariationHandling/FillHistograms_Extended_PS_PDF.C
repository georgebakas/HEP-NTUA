#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1.h"
#include "TCanvas.h"
//for the pdf and ps weights I need to loop all over the weights!!
//for each weight I get a ps or pdf value and fill the histogram with it


using std::cin;
using std::cout;
using std::endl;


std::vector<TString> histoNames;
std::vector<TString> fileNames;

#include "TemplateConstants.h"
TVector3 getBoostVector(TLorentzVector p4_1, TLorentzVector p4_2, TLorentzVector &p4CombinedVector);

TString globalYear;

void FillHistograms_Extended_PS_PDF(TString file_name, TString ttbar_process, TString year = "2016", TString weightType="PSWeights", float mJJCut = 1000)
{
  globalYear = year;
  initFilesMapping();
  cout<<"YEAR: "<<year<<endl;
  cout<<"file_name: "<<file_name<<endl;
  cout<<"ttbar_process: "<<ttbar_process<<endl;
  float triggerFloat;
  if(year.Contains("2016")) triggerFloat = 2;
  else triggerFloat = 5;

  float deepCSVFloat = floatConstants[TString::Format("btagWP%s",year.Data())];
  float selMvaCut = topTaggerConstants[TString::Format("topTagger%s",year.Data())];
  float LUMI = luminosity[TString::Format("luminosity%s", year.Data())];
  float LUMI_CR = luminosityCR[TString::Format("luminosity%s", year.Data())];
  float XSEC = XSECAll[year.Data()][ttbar_process.Data()];
  cout<<"XSEC: "<<XSEC<<endl;
  const int NVAR = 12;
  std::vector< std::vector <Float_t> > const BND = {{1000, 1200, 1400, 1600, 1800, 2000, 2400, 3000, 5000},
                       {0, 60, 150, 300, 450, 850, 1300}, //ptJJ
                       {-2.4, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.4}, //yjj
                       {400, 450, 500, 570, 650, 800, 1100, 1500}, //jetPt0
                       {400, 450, 500, 570, 650, 800, 1100, 1500}, //jetPt1
                       {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.4}, //jetY0
                       {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.4}, //jetY1
                       {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 16}, //chi
                       {-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1}, //|cosTheta*| leading
                       {-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1}}; //|cosTheta*| subleading; //jetY1

  int NBINS[BND.size()+2];
  for (int i = 0; i<BND.size()+2; i++)
  {
    if(i<10) NBINS[i] = BND[i].size()-1;
    else NBINS[i] = 100;
  }
  TString varReco[NVAR]   = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1","chi","cosThjetEta0", "cosThjetEta1","mTop_Leading", "mTop_Subleading"};

  float weights;
  int maxJetMass = 300;
  int minJetMass = 50;

    int nJets,nLeptons, category(0);
    vector<bool>  *bit(0),*matchedJet(0);
    //reco vars:
    std::vector<float> *jetPt(0), *jetY(0), *jetEta(0), *jetPhi(0), *jetTtag(0);
    float genEvtWeight(0);
    double  bTagEvntWeight(0);
    float mJJ(0), ptJJ(0), yJJ(0),mva(0);
    vector<float> *tau3(0),*tau2(0),*tau1(0);
	  vector<float> *jetMassSub0(0), *jetMassSub1(0), *jetBtagSub0(0), *jetBtagSub1(0);
    vector<float> *jetMassSoftDrop(0), *partonEta(0);

    //parton
    std::vector<float> *partonPt(0), *partonPhi(0), *partonMass(0),*partonMatchDR(0);
    float yTTbarParton(0), ptTTbarParton(0), mTTbarParton(0);

    //float yTopParton[2], ptTopParton[2];
    std::vector<int> *addedIndexes = new std::vector<int>(0);
    vector<int> *partonId(0), *partonMatchIdx(0);

    std::vector<float> *jetBtagSub0DCSVbb(0), *jetBtagSub1DCSVbb(0);
    std::vector<float> *jetBtagSub0DCSVbbb(0), *jetBtagSub1DCSVbbb(0);

    //ps_weights
    std::vector<float> *psWeights(0),*pdfWeights(0), *scaleWeights(0);

    int nJetsGen(0);
    float mJJGen(0), ptJJGen(0), yJJGen(0);
    TFile *file = TFile::Open(file_name.Data());

    TTree *trIN = (TTree*)file->Get("boosted/events");

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
    trIN->SetBranchAddress("psWeights"      ,&psWeights);
    trIN->SetBranchAddress("pdfWeights"     ,&pdfWeights);
    trIN->SetBranchAddress("scaleWeights"   ,&scaleWeights);
    trIN->SetBranchAddress("bTagEvntWeight" ,&bTagEvntWeight);
    trIN->SetBranchAddress("jetMassSub0"    ,&jetMassSub0);
    trIN->SetBranchAddress("jetMassSub1"    ,&jetMassSub1);
    trIN->SetBranchAddress("mJJ"   			,&mJJ);
    trIN->SetBranchAddress("yJJ"   			,&yJJ);
    trIN->SetBranchAddress("ptJJ"   		,&ptJJ);
	trIN->SetBranchAddress("jetBtagSub0"	,&jetBtagSub0);
    //trIN->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1);
    trIN->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
  	trIN->SetBranchAddress("jetTtagCategory",&jetTtag);
  	//deepCSV
    trIN->SetBranchAddress("jetBtagSub0DCSVbb" ,&jetBtagSub0DCSVbb);
    trIN->SetBranchAddress("jetBtagSub1DCSVbb" ,&jetBtagSub1DCSVbb);
    trIN->SetBranchAddress("jetBtagSub0DCSVbbb",&jetBtagSub0DCSVbbb);
    trIN->SetBranchAddress("jetBtagSub1DCSVbbb",&jetBtagSub1DCSVbbb);
    //general
    int evtNo;
	trIN->SetBranchAddress("evtNo", &evtNo);

    //parton
	trIN->SetBranchAddress("mTTbarParton"	,&mTTbarParton);
	trIN->SetBranchAddress("yTTbarParton"	,&yTTbarParton);
	trIN->SetBranchAddress("ptTTbarParton"	,&ptTTbarParton);
	trIN->SetBranchAddress("partonPt"	    ,&partonPt);
	trIN->SetBranchAddress("partonEta"		,&partonEta);
	trIN->SetBranchAddress("partonMass"     ,&partonMass);
	trIN->SetBranchAddress("partonMatchDR"  ,&partonMatchDR);
	trIN->SetBranchAddress("partonMatchIdx" ,&partonMatchIdx);
	trIN->SetBranchAddress("partonPhi"      ,&partonPhi);

  float norm = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
	float weight = XSEC/norm;
	weights = weight;


  //declare the histograms
  /*
  trIN->GetEntry(1);
  if(weightType.EqualTo("PSWeights")) weightsSize = weightsSize;
  else weightsSize = pdfWeights->size(); */
  int weightsSize(1);
  if(weightType.EqualTo("PSWeights")) weightsSize = ps_weights.size();
  else if(weightType.EqualTo("PDFWeights")) weightsSize = pdf_weights.size();
  else weightsSize = scale_weights.size();

  TH1F *hReco[weightsSize][NVAR], *hRecoCR[weightsSize][NVAR];
  for(int ivar =0; ivar<BND.size()+2; ivar++)
  {
    for(int iweight=0; iweight<weightsSize; iweight++)
    {
      TString weightName;
      if(weightType.EqualTo("PSWeights")) weightName = ps_weights[iweight];
      else if(weightType.EqualTo("PDFWeights")) weightName = pdf_weights[iweight];
      else weightName = scale_weights[iweight];

      int sizeBins = NBINS[ivar];
      if(ivar>=10)
      {
        hReco[iweight][ivar] = new TH1F(TString::Format("hReco_%s_%s", varReco[ivar].Data(), weightName.Data()), TString::Format("hReco_%s_%s",varReco[ivar].Data(), weightName.Data()), sizeBins, minJetMass,maxJetMass);
        hRecoCR[iweight][ivar] = new TH1F(TString::Format("hRecoCR_%s_%s", varReco[ivar].Data(), weightName.Data()), TString::Format("hRecoCR_%s_%s",varReco[ivar].Data(), weightName.Data()), sizeBins, minJetMass,maxJetMass);
      }
      else
      {
        int sizeBins = NBINS[ivar];
        float tempBND[NBINS[ivar]+1];
        std::copy(BND[ivar].begin(), BND[ivar].end(), tempBND);
        hReco[iweight][ivar] = new TH1F(TString::Format("hReco_%s_%s", varReco[ivar].Data(), weightName.Data()), TString::Format("hReco_%s_%s",varReco[ivar].Data(), weightName.Data()), sizeBins, tempBND);
        hRecoCR[iweight][ivar] = new TH1F(TString::Format("hRecoCR_%s_%s", varReco[ivar].Data(), weightName.Data()), TString::Format("hRecoCR_%s_%s",varReco[ivar].Data(), weightName.Data()), sizeBins, tempBND);
      }
      hReco[iweight][ivar]->Sumw2();
      hRecoCR[iweight][ivar]->Sumw2();
    }
  }


    int decade(0);
    int NN = trIN->GetEntries();
    //NN = 100000;
    std::cout<<"Entries: "<<NN<<std::endl;
	std::vector<float> xRecoAll(0);

    for(int iev=0;iev<NN;iev++)
    {
		double progress = 10.0*iev/(1.0*NN);
      int k = TMath::FloorNint(progress);
      if (k > decade)
        cout<<10*k<<" %"<<endl;
      decade = k;
      trIN->GetEntry(iev);

	   xRecoAll.clear();

	if(nJets > 1)
    {
   	   	int leadingPt =0;
      	int subleadingPt = 1;

  		if((*jetPt)[0] < (*jetPt)[1])
  		{
   		    subleadingPt =0;
   		    leadingPt = 1;
   		}

		  TLorentzVector p4T[2], p4T_ZMF[2], p4TTbar;
    	p4T[leadingPt].SetPtEtaPhiM((*jetPt)[leadingPt], (*jetEta)[leadingPt], (*jetPhi)[leadingPt], (*jetMassSoftDrop)[leadingPt]);
   		p4T[subleadingPt].SetPtEtaPhiM((*jetPt)[subleadingPt], (*jetEta)[subleadingPt], (*jetPhi)[subleadingPt], (*jetMassSoftDrop)[subleadingPt]);
  	  TVector3 ttbarBoostVector = getBoostVector(p4T[leadingPt], p4T[subleadingPt], p4TTbar);

    	p4T_ZMF[0].SetPtEtaPhiM(p4T[leadingPt].Pt(), p4T[leadingPt].Eta(), p4T[leadingPt].Phi(), p4T[leadingPt].M());
    	p4T_ZMF[1].SetPtEtaPhiM(p4T[subleadingPt].Pt(), p4T[subleadingPt].Eta(), p4T[subleadingPt].Phi(), p4T[subleadingPt].M());
     	p4T_ZMF[0].Boost(ttbarBoostVector);
     	p4T_ZMF[1].Boost(ttbarBoostVector);
	    float yStarExp  = TMath::Exp(fabs(p4T_ZMF[0].Rapidity() - p4T_ZMF[1].Rapidity())); //this is chi = e^(|y*|) , y* = 1/2(y1-y0)

      xRecoAll.push_back(mJJ);
      xRecoAll.push_back(ptJJ);
      xRecoAll.push_back(yJJ);
      xRecoAll.push_back((*jetPt)[leadingPt]);
      xRecoAll.push_back((*jetPt)[subleadingPt]);
      xRecoAll.push_back(fabs((*jetY)[leadingPt]));
      xRecoAll.push_back(fabs((*jetY)[subleadingPt]));
     	xRecoAll.push_back(yStarExp); //this is chi
     	xRecoAll.push_back(TMath::Cos(p4T_ZMF[0].Theta())); //this is |cos(theta*)| leading
     	xRecoAll.push_back(TMath::Cos(p4T_ZMF[1].Theta())); //this is |cos(theta*)| subleading
      xRecoAll.push_back((*jetMassSoftDrop)[leadingPt]);
      xRecoAll.push_back((*jetMassSoftDrop)[subleadingPt]);

	  //---------------------------end of MATCHING---------------------------------------------------------
    float dCSVScoreSub0[2], dCSVScoreSub1[2];
    dCSVScoreSub0[0] = (*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0];
    dCSVScoreSub0[1] = (*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1];
    dCSVScoreSub1[0] = (*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0];
    dCSVScoreSub1[1] = (*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1];

	  bool recoCuts;
	  bool massCut = (*jetMassSoftDrop)[0] > 50 && (*jetMassSoftDrop)[0] < 300 && (*jetMassSoftDrop)[1] > 50 && (*jetMassSoftDrop)[1] < 300;
	  bool tTaggerCut = (*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut;
	  recoCuts = nJets > 1 && fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1]) <2.4 && (*jetPt)[0] > 400 && (*jetPt)[1] > 400 && mJJ> mJJCut && massCut && nLeptons==0 && (*bit)[triggerFloat];
	  bool deepCSV = (((*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0])> deepCSVFloat) &&
					 (((*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1])> deepCSVFloat || ((*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1])> deepCSVFloat);
    bool revertBtag = (dCSVScoreSub0[0] < deepCSVFloat &&  dCSVScoreSub1[0] < deepCSVFloat) && (dCSVScoreSub0[1] < deepCSVFloat && dCSVScoreSub1[1] < deepCSVFloat);
      bool btagCut;
	  btagCut = deepCSV;
    float extra_weight(1);
    //Signal Region 2btags
		if(recoCuts && btagCut && tTaggerCut)
		{
		  for(int ivar = 0; ivar < NVAR; ivar++)
	  	{
        for(int iweight=0; iweight<weightsSize; iweight++)
        {
          if(weightType.EqualTo("PSWeights")) extra_weight = (*psWeights)[iweight];
          else if(weightType.EqualTo("PDFWeights")) extra_weight = (*pdfWeights)[iweight];
          else extra_weight = (*scaleWeights)[iweight];
		      hReco[iweight][ivar]->Fill(xRecoAll[ivar], genEvtWeight*bTagEvntWeight*extra_weight);
        }
		  }
	  }
    //Control Region 0btag
    if(recoCuts && revertBtag && tTaggerCut)
	  {
	  	for(int ivar = 0; ivar < NVAR; ivar++)
  		{
        for(int iweight=0; iweight<weightsSize; iweight++)
        {
          if(weightType.EqualTo("PSWeights")) extra_weight = (*psWeights)[iweight];
          else if(weightType.EqualTo("PDFWeights")) extra_weight = (*pdfWeights)[iweight];
          else extra_weight = (*scaleWeights)[iweight];
		      hRecoCR[iweight][ivar]->Fill(xRecoAll[ivar], genEvtWeight*bTagEvntWeight*extra_weight);
        }
		  }
	  }

  }//----- end of nJets>1
}//---end the event for


  for(int ivar =0; ivar<NVAR; ivar++)
  {
    for(int iweight=0; iweight<weightsSize; iweight++)
    {
      hReco[iweight][ivar]->Scale(weights*LUMI);
      hRecoCR[iweight][ivar]->Scale(weights*LUMI_CR);
    }
  }//end of loop on all vars


  TFile *outFile[weightsSize];
  for(int iweight=0; iweight<weightsSize; iweight++ )
  {
    TString weightName;
    if(weightType.EqualTo("PSWeights")) weightName = ps_weights[iweight];
    else if(weightType.EqualTo("PDFWeights")) weightName = pdf_weights[iweight];
    else weightName = scale_weights[iweight];

    outFile[iweight] = TFile::Open(TString::Format("%s/%s/Histo_%d_%s_%s.root", year.Data(), weightType.Data(), (int)mJJCut, ttbar_process.Data(), weightName.Data()), "RECREATE");
    outFile[iweight]->cd();
    //write them to file
    for(int ivar = 0; ivar<NVAR; ivar++)
    {
      hReco[iweight][ivar]->Write(TString::Format("hWt_%s_2btag", varReco[ivar].Data()));
      hRecoCR[iweight][ivar]->Write(TString::Format("hWt_%s_0btag", varReco[ivar].Data()));
    }//end of ivar loop
    outFile[iweight]->Close();
  }//end of iweight loop

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
