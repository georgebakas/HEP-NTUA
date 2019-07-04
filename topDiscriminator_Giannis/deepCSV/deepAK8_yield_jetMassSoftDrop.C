#include "BoostedAnalysisInitializer.h"
R__LOAD_LIBRARY(BoostedAnalysisInitializer_cxx)

void savePlot(TH1F* histogram, TString name, TFile* file)
{
  if(!file->GetListOfKeys()->Contains(name))
  {
    histogram->SetName(name);
    histogram->SetTitle(name);
    histogram->Write(name);
  }
}

void deepAK8_yield_jetMassSoftDrop(bool isSignal, float deepAK8Cut = 0.6, float tTaggerCut = 0.2, float mvaCut = 0.8)
{
  BoostedAnalysisInitializer utils("2018", isSignal);
  
  TH1F *h_out[3][utils.listOfFiles.size()];
  
  TH1F* h_denominator_parton[utils.listOfFiles.size()];
  TH1F* h_denominator_reco[utils.listOfFiles.size()];
  std::vector<float> weights;
  
  int sizeBins = 20;
  float BND[11] = {1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, 4000, 5000};
  int MIN = 50; 
  int MAX = 250;
  
  BoostedTopCategoryDiscriminator_tTagger *tagger = new BoostedTopCategoryDiscriminator_tTagger("/afs/cern.ch/user/g/gbakas/public/ForGiannis/Training_2018/BoostedMVA/weights/boosted_MVA_BDTCat.weights.xml");
  
  int passedEvents = 0;
  for(int f=0; f<utils.listOfFiles.size(); f++)
  {
    h_out[0][f] = new TH1F(TString::Format("%s_%s_reco", utils.histoNames[f].Data(), "tTagger"), TString::Format("%s_%s_reco", utils.histoNames[f].Data(), "tTagger"), sizeBins, MIN, MAX);
    h_out[1][f] = new TH1F(TString::Format("%s_%s_reco", utils.histoNames[f].Data(), "deepAK8"), TString::Format("%s_%s_reco", utils.histoNames[f].Data(), "deepAK8"), sizeBins, MIN, MAX);
    h_out[2][f] = new TH1F(TString::Format("%s_%s_reco", utils.histoNames[f].Data(), "eventMVA"), TString::Format("%s_%s_reco", utils.histoNames[f].Data(), "eventMVA"), sizeBins, MIN, MAX);
    
    
    int nJets,nLeptons;
    std::vector<float> *jetPt(0), *jetEta(0), *jetPhi(0), *deepAK8(0);
    std::vector<float> *mass(0), *tau3(0), *tau2(0), *tau1(0);
    std::vector<float> *jetMassSub0(0), *jetMassSub1(0);
    std::vector<float> /* *jetTtag(0),*/ *jetMassSoftDrop(0);
    std::vector<float> *jetBtagSub0DCSVbb(0), *jetBtagSub0DCSVbbb(0);
    std::vector<float> *jetBtagSub1DCSVbb(0), *jetBtagSub1DCSVbbb(0);
    std::vector<float> *ecfB1N2(0), *ecfB1N3(0), *ecfB2N2(0), *ecfB2N3(0);
    
    float mTTbarParton(0), mJJ(0), mva(0), eventWeight(0), ht(0);
    
    std::vector<float> *jetTtag = new std::vector<float>();
    
    std::cout<<"Working in file: "<<utils.eosPath+utils.listOfFiles[f]<<std::endl;
    TFile *file = TFile::Open(utils.eosPath+utils.listOfFiles[f]);
    TTree *trIN = (TTree*)file->Get("boosted/events");
    
    trIN->SetBranchAddress("nJets"          ,&nJets);
    trIN->SetBranchAddress("nLeptons"       ,&nLeptons);
    trIN->SetBranchAddress("jetPt"          ,&jetPt);
    trIN->SetBranchAddress("deepAK8Tagger"  ,&deepAK8);
    trIN->SetBranchAddress("mva"            ,&mva);
    
    trIN->SetBranchAddress("jetEta"         ,&jetEta);
    trIN->SetBranchAddress("jetPhi"         ,&jetPhi);
    trIN->SetBranchAddress("jetMassSoftDrop",&mass);
    trIN->SetBranchAddress("jetTau3"        ,&tau3);
    trIN->SetBranchAddress("jetTau2"        ,&tau2);
    trIN->SetBranchAddress("jetTau1"        ,&tau1);
    trIN->SetBranchAddress("ecfB1N2"        ,&ecfB1N2);
    trIN->SetBranchAddress("ecfB1N3"        ,&ecfB1N3);
    trIN->SetBranchAddress("ecfB2N2"        ,&ecfB2N2);
    trIN->SetBranchAddress("ecfB2N3"        ,&ecfB2N3);
    trIN->SetBranchAddress("jetMassSub0"    ,&jetMassSub0);
    trIN->SetBranchAddress("jetMassSub1"    ,&jetMassSub1);
	  trIN->SetBranchAddress("mTTbarParton"   ,&mTTbarParton);
    trIN->SetBranchAddress("mJJ"   			    ,&mJJ);
    trIN->SetBranchAddress("jetBtagSub0DCSVbb", &jetBtagSub0DCSVbb);
    trIN->SetBranchAddress("jetBtagSub0DCSVbbb", &jetBtagSub0DCSVbbb);
    trIN->SetBranchAddress("jetBtagSub1DCSVbb", &jetBtagSub1DCSVbb);
    trIN->SetBranchAddress("jetBtagSub1DCSVbbb", &jetBtagSub1DCSVbbb);
    trIN->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
    //trIN->SetBranchAddress("jetTtagCategory",&jetTtag);
    trIN->SetBranchAddress("ht", &ht);
    trIN->SetBranchAddress("genEvtWeight", &eventWeight);
    
    float norm = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
	  float weight = utils.XSEC[f]/norm;
	  weights.push_back(weight);
    
    int decade(0);
    int NN = trIN->GetEntries();
    //NN = 10000;
    float floatBTag = 0.4184;
    std::cout<<"Entries: "<<NN<<std::endl;
    
    for(int iev=0; iev<NN; iev++)
    {
      double progress = 10.0*iev/(1.0*NN);
      int k = TMath::FloorNint(progress); 
      if (k > decade) 
        cout<<10*k<<" %"<<endl;
      decade = k;
      trIN->GetEntry(iev);
      
      if(nJets > 1 && fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1] <2.4) && (*jetPt)[0] > 400 && (*jetPt)[1] > 400 && mJJ > 1000 && nLeptons == 0/*&& (*triggerBit)[2]*/ )
      {
        jetTtag->clear();
        for(int i=0; i<nJets; i++)
        {
          if((*jetPt)[i] >= 400)
          {
            float jetPtOverSumPt = (*jetPt)[i]/ht;
            float score = tagger->eval((*tau3)[i], (*tau2)[i], (*tau1)[i], (*jetMassSub0)[i], (*jetMassSub1)[i], (*ecfB1N2)[i], (*ecfB1N3)[i], (*ecfB2N2)[i], (*ecfB2N3)[i], jetPtOverSumPt, (*jetPt)[i]);
          
            jetTtag->push_back(score);
          }
        }
        
        if((*jetTtag)[0] > tTaggerCut && (*jetTtag)[1] > tTaggerCut && 
           (((*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0])> floatBTag || 
            ((*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0])> floatBTag) && 
           (((*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1])> floatBTag || 
            ((*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1])> floatBTag))
        {
          h_out[0][f]->Fill((*jetMassSoftDrop)[0], eventWeight);
          h_out[0][f]->Fill((*jetMassSoftDrop)[1], eventWeight);
        }
        
        if((*deepAK8)[0] > deepAK8Cut && (*deepAK8)[1] > deepAK8Cut && 
           (((*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0])> floatBTag || 
            ((*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0])> floatBTag) && 
           (((*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1])> floatBTag || 
            ((*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1])> floatBTag))
        {
          h_out[1][f]->Fill((*jetMassSoftDrop)[0], eventWeight);
          h_out[1][f]->Fill((*jetMassSoftDrop)[1], eventWeight);
        }
        
        if(mva > mvaCut && 
           (((*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0])> floatBTag || 
            ((*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0])> floatBTag) && 
           (((*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1])> floatBTag || 
            ((*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1])> floatBTag))
        {
          //std::cout<<"Passed: "<<((*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0])<<" | "<<((*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0]);
          //std::cout<<" | "<<((*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1])<<" | "<<((*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1])<<std::endl;
          h_out[2][f]->Fill((*jetMassSoftDrop)[0], eventWeight);
          h_out[2][f]->Fill((*jetMassSoftDrop)[1], eventWeight);
          passedEvents++;
        }
      }
    }
    file->Close();
  }
  
  std::cout<<"LUMI: "<<utils.LUMI<<std::endl;
  for(int i=0; i<sizeof(h_out)/sizeof(*h_out); i++)
  {
	//for every slice
  	for(int j=0; j<sizeof(*h_out)/sizeof(*h_out[0]); j++)
    {
      h_out[i][j]->Scale(weights[j]*utils.LUMI);
    }
    
    for(int j=1; j<sizeof(*h_out)/sizeof(*h_out[0]); j++)
	  {
	    h_out[i][0]->Add(h_out[i][j]);
	  }
  }
  
  std::cout<<"Passed events: "<<passedEvents<<std::endl;
  
  TFile *outFile = TFile::Open("deepAK8_yields_jetMassSoftDrop.root", "UPDATE");
  TString prefix = (isSignal?"Sig":"Bkg");
  
  savePlot(h_out[0][0], TString::Format("%s_yield_jetMassSoftDrop_tTagger_%.1f", prefix.Data(), tTaggerCut), outFile);
  savePlot(h_out[2][0], TString::Format("%s_yield_jetMassSoftDrop_eventMVA_%.1f", prefix.Data(), mvaCut), outFile);
  savePlot(h_out[1][0], TString::Format("%s_yield_jetMassSoftDrop_deepAK8_%.1f", prefix.Data(), deepAK8Cut), outFile);
  
  outFile->Close();
}