#include "BoostedAnalysisInitializer.h"
R__LOAD_LIBRARY(BoostedAnalysisInitializer_cxx)

void deepAK8_efficienciesOnlyRecoCuts(float deepAK8Cut = 0.6, float tTaggerCut = 0.2, float mvaCut = 0.8)
{
  BoostedAnalysisInitializer utils("2018", true);
  
  TH1F *h_out_reco[3][utils.listOfFiles.size()];
  
  TH1F* h_denominator_reco[utils.listOfFiles.size()];
  std::vector<float> weights;
  
  int sizeBins = 10;
  float BND[11] = {1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, 4000, 5000};
  int MIN = 1000; 
  int MAX = 5000;
  
  BoostedTopCategoryDiscriminator_tTagger *tagger = new BoostedTopCategoryDiscriminator_tTagger("/afs/cern.ch/work/i/ipapakri/private/analysis/topDiscriminator/deepAK8/boosted_MVA_BDTCat.weights.xml");
  
  for(int f=0; f<utils.listOfFiles.size(); f++)
  {
    h_out_reco[0][f] = new TH1F(TString::Format("%s_%s_reco", utils.histoNames[f].Data(), "tTagger"), TString::Format("%s_%s_reco", utils.histoNames[f].Data(), "tTagger"), sizeBins, BND);
    h_out_reco[1][f] = new TH1F(TString::Format("%s_%s_reco", utils.histoNames[f].Data(), "deepAK8"), TString::Format("%s_%s_reco", utils.histoNames[f].Data(), "deepAK8"), sizeBins, BND);
    h_out_reco[2][f] = new TH1F(TString::Format("%s_%s_reco", utils.histoNames[f].Data(), "eventMVA"), TString::Format("%s_%s_reco", utils.histoNames[f].Data(), "eventMVA"), sizeBins, BND);
    
    h_denominator_reco[f] = new TH1F(TString::Format("denom_%s_reco", utils.histoNames[f].Data()), TString::Format("denom_%s_reco", utils.histoNames[f].Data()), sizeBins, BND);
	  h_denominator_reco[f]->Sumw2();
    
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
    TTree *eventTree = (TTree*)file->Get("boosted/events");
    
    eventTree->SetBranchAddress("nJets"          ,&nJets);
    eventTree->SetBranchAddress("nLeptons"       ,&nLeptons);
    eventTree->SetBranchAddress("jetPt"          ,&jetPt);
    eventTree->SetBranchAddress("deepAK8Tagger"        ,&deepAK8);
    eventTree->SetBranchAddress("mva"            ,&mva);
    
    eventTree->SetBranchAddress("jetEta"         ,&jetEta);
    eventTree->SetBranchAddress("jetPhi"         ,&jetPhi);
    eventTree->SetBranchAddress("jetMassSoftDrop",&mass);
    eventTree->SetBranchAddress("jetTau3"        ,&tau3);
    eventTree->SetBranchAddress("jetTau2"        ,&tau2);
    eventTree->SetBranchAddress("jetTau1"        ,&tau1);
    eventTree->SetBranchAddress("ecfB1N2"        ,&ecfB1N2);
    eventTree->SetBranchAddress("ecfB1N3"        ,&ecfB1N3);
    eventTree->SetBranchAddress("ecfB2N2"        ,&ecfB2N2);
    eventTree->SetBranchAddress("ecfB2N3"        ,&ecfB2N3);
    eventTree->SetBranchAddress("jetMassSub0"    ,&jetMassSub0);
    eventTree->SetBranchAddress("jetMassSub1"    ,&jetMassSub1);
	  eventTree->SetBranchAddress("mTTbarParton"   ,&mTTbarParton);
    eventTree->SetBranchAddress("mJJ"   			   ,&mJJ);
    eventTree->SetBranchAddress("jetBtagSub0DCSVbb", &jetBtagSub0DCSVbb);
    eventTree->SetBranchAddress("jetBtagSub0DCSVbbb", &jetBtagSub0DCSVbbb);
    eventTree->SetBranchAddress("jetBtagSub1DCSVbb", &jetBtagSub1DCSVbb);
    eventTree->SetBranchAddress("jetBtagSub1DCSVbbb", &jetBtagSub1DCSVbbb);
    eventTree->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
    eventTree->SetBranchAddress("ht", &ht);
    /*eventTree->SetBranchAddress("jetTtagCategory",&jetTtag);*/
    eventTree->SetBranchAddress("genEvtWeight", &eventWeight);
    
    float norm = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
	  float weight = utils.XSEC[f]/norm;
	  weights.push_back(weight);
    
    int decade(0);
    int NN = eventTree->GetEntries();
    //NN = 100000;
    float floatBTag = 0.4184;
    std::cout<<"Entries: "<<NN<<std::endl;
    for(int iev=0;iev<NN;iev++) 
    {
      double progress = 10.0*iev/(1.0*NN);
      int k = TMath::FloorNint(progress); 
      if (k > decade) 
        cout<<10*k<<" %"<<endl;
      decade = k;
      eventTree->GetEntry(iev);
      
      if(nJets > 1 && fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1] <2.4) && (*jetPt)[0] > 400 && (*jetPt)[1] > 400 && mJJ > 1000 && nLeptons == 0/*&& (*triggerBit)[2]*/ )
	    {
        h_denominator_reco[f]->Fill(mJJ, eventWeight);
        //1. category where both are fully top tagged
          
		    if((*jetMassSoftDrop)[0] > 120 && (*jetMassSoftDrop)[0] < 220 && (*jetMassSoftDrop)[1] > 120 && (*jetMassSoftDrop)[1] < 220)
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
            h_out_reco[0][f]->Fill(mJJ, eventWeight);
		      }
          
          if((*deepAK8)[0] > deepAK8Cut && (*deepAK8)[1] > deepAK8Cut && 
             (((*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0])> floatBTag || 
              ((*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0])> floatBTag) && 
             (((*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1])> floatBTag || 
              ((*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1])> floatBTag))
          {
            h_out_reco[1][f]->Fill(mJJ, eventWeight);
          }
          
          if(mva > mvaCut && 
             (((*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0])> floatBTag || 
              ((*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0])> floatBTag) && 
             (((*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1])> floatBTag || 
              ((*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1])> floatBTag))
          {
            h_out_reco[2][f]->Fill(mJJ, eventWeight);
          }
        }
      }
    }
    file->Close();
  }
  
  for(int i=0; i<sizeof(h_out_reco)/sizeof(*h_out_reco); i++)
  {
	//for every slice
	  for(int j=0; j<sizeof(*h_out_reco)/sizeof(*h_out_reco[0]); j++)
    {
      h_out_reco[i][j]->Scale(weights[j]);
    }
    
    for(int j=1; j<sizeof(*h_out_reco)/sizeof(*h_out_reco[0]); j++)
	  {
	  //Add them to get the whole phase space
      std::cout<<h_out_reco[i][0]->GetName()<<" "<<h_out_reco[i][j]->GetName()<<std::endl;
      h_out_reco[i][0]->Add(h_out_reco[i][j]);
	  }
  }
  
  h_denominator_reco[0]->Scale(weights[0]);
  
  for(int i=1; i<sizeof(h_denominator_reco)/sizeof(*h_denominator_reco); i++) 
  {
    h_denominator_reco[i]->Scale(weights[i]);
	  h_denominator_reco[0]->Add(h_denominator_reco[i]);
  }
  
  TEfficiency *efficiency_reco[3];
  
  efficiency_reco[0]  = new TEfficiency(*h_out_reco[0][0], *h_denominator_reco[0]);
  efficiency_reco[0]->SetTitle("reco_tTagger;mJJ (GeV);Efficiency");
  efficiency_reco[0]->SetStatisticOption(TEfficiency::kFNormal);
  efficiency_reco[0]->SetUseWeightedEvents();
  efficiency_reco[0]->SetLineColor(kRed);
  
  efficiency_reco[1]  = new TEfficiency(*h_out_reco[1][0], *h_denominator_reco[0]);
  efficiency_reco[1]->SetTitle("reco_deepAK8;mJJ (GeV);Efficiency");
  efficiency_reco[1]->SetStatisticOption(TEfficiency::kFNormal);
  efficiency_reco[1]->SetUseWeightedEvents();
  efficiency_reco[1]->SetLineColor(kBlue);
  
  efficiency_reco[2]  = new TEfficiency(*h_out_reco[2][0], *h_denominator_reco[0]);
  efficiency_reco[2]->SetTitle("reco_eventMVA;mJJ (GeV);Efficiency");
  efficiency_reco[2]->SetStatisticOption(TEfficiency::kFNormal);
  efficiency_reco[2]->SetUseWeightedEvents();
  efficiency_reco[2]->SetLineColor(kBlack); 
  
  /*TCanvas *c2 = new TCanvas("recoCanvas", "recoCanvas", 600, 500);
  c2->cd();
  efficiency_reco[0]->Draw();
  efficiency_reco[1]->Draw("SAME");
  efficiency_reco[2]->Draw("SAME");*/
  
  TFile *outFile = TFile::Open("deepAK8_efficienciesOnlyRecoCuts.root", "UPDATE");
  
  if(!outFile->GetListOfKeys()->Contains(TString::Format("Sig_Reco_tTagger_%.1f", tTaggerCut)))
  {
    efficiency_reco[0]->Write(TString::Format("Sig_Reco_tTagger_%.1f", tTaggerCut));
  }
  
  if(!outFile->GetListOfKeys()->Contains(TString::Format("Sig_Reco_eventMVA_%.1f", mvaCut)))
  {
    efficiency_reco[2]->Write(TString::Format("Sig_Reco_eventMVA_%.1f", mvaCut));
  }
  
  if(!outFile->GetListOfKeys()->Contains(TString::Format("Sig_Reco_deepAK8_%.1f", deepAK8Cut)))
  {
    efficiency_reco[1]->Write(TString::Format("Sig_Reco_deepAK8_%.1f", deepAK8Cut));
  }
  
}