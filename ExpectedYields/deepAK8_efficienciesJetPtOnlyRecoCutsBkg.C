//#include "BoostedAnalysisInitializer.h"
//R__LOAD_LIBRARY(BoostedAnalysisInitializer_cxx)
#include "BoostedAnalysisInitializer.h"

void deepAK8_efficienciesJetPtOnlyRecoCutsBkg(TString inYear = "2017",float selMvaCut = 0.0, float deepAK8Cut = 0.6,  float mvaCut = 0.8)
{
	TString year = inYear;
  BoostedAnalysisInitializer utils(year.Data(), false);
  
  TH1F *h_out_reco[3][utils.listOfFiles.size()];
  
  TH1F* h_denominator_reco[utils.listOfFiles.size()];
  std::vector<float> weights;
  
  const int N_PT = 10;
	float BND_PT[N_PT+1] = {400,450,500,570,650,750,850,950,1100,1300,1500};
  /*
  BoostedTopCategoryDiscriminator_tTagger *tagger = new BoostedTopCategoryDiscriminator_tTagger("/afs/cern.ch/work/i/ipapakri/private/analysis/topDiscriminator/deepAK8/boosted_MVA_BDTCat.weights.xml");
  */
  int deepAK8events = 0;
  int tTaggerEvents = 0;
  float floatBTag;
  
  if (year.EqualTo("2018"))
	floatBTag = 0.4184;
  else if(year.EqualTo("2017"))
	floatBTag = 0.4941;
  else if(year.EqualTo("2016"))
	floatBTag = 0.6321;

  cout<<selMvaCut<<endl;
  cout<<floatBTag<<endl;
  cout<<year<<endl;
  for(int f=0; f<utils.listOfFiles.size(); f++)
  {
    h_out_reco[0][f] = new TH1F(TString::Format("%s_%s_reco_jetPt", utils.histoNames[f].Data(), "tTagger"), TString::Format("%s_%s_reco_jetPt", utils.histoNames[f].Data(), "tTagger"), N_PT, BND_PT);
    h_out_reco[1][f] = new TH1F(TString::Format("%s_%s_reco_jetPt", utils.histoNames[f].Data(), "deepAK8"), TString::Format("%s_%s_reco_jetPt", utils.histoNames[f].Data(), "deepAK8"), N_PT, BND_PT);
    h_out_reco[2][f] = new TH1F(TString::Format("%s_%s_reco_jetPt", utils.histoNames[f].Data(), "eventMVA"), TString::Format("%s_%s_reco_jetPt", utils.histoNames[f].Data(), "eventMVA"), N_PT, BND_PT);
    
    h_denominator_reco[f] = new TH1F(TString::Format("denom_%s_reco_jetPt", utils.histoNames[f].Data()), TString::Format("denom_%s_reco_jetPt", utils.histoNames[f].Data()), N_PT, BND_PT);
	  h_denominator_reco[f]->Sumw2();
    
    int nJets,nLeptons;
    std::vector<float> *jetPt(0), *jetEta(0), *jetPhi(0), *deepAK8(0);
    std::vector<float> *mass(0), *tau3(0), *tau2(0), *tau1(0);
    std::vector<float> *jetMassSub0(0), *jetMassSub1(0);
    std::vector<float> /* *jetTtag(0),*/ *jetMassSoftDrop(0);
    std::vector<float> *jetBtagSub0DCSVbb(0), *jetBtagSub0DCSVbbb(0);
    std::vector<float> *jetBtagSub1DCSVbb(0), *jetBtagSub1DCSVbbb(0);
    std::vector<float> *ecfB1N2(0), *ecfB1N3(0), *ecfB2N2(0), *ecfB2N3(0);
    
    float mTTbarParton(0), mJJ(0), mva(0), ht(0), eventWeight(0);
    
    std::vector<float> *jetTtag = new std::vector<float>();
    
    std::cout<<"Working in file: "<<utils.eosPath+utils.listOfFiles[f]<<std::endl;
    TFile *file = TFile::Open(utils.eosPath+utils.listOfFiles[f]);
    TTree *trIN;
	if(year.EqualTo("2017") || year.EqualTo("2016")) trIN = (TTree*)file->Get("boosted/events");
	else if(year.EqualTo("2018")) trIN = (TTree*)file->Get("boosted/events");
    
    trIN->SetBranchAddress("nJets"          ,&nJets);
    trIN->SetBranchAddress("nLeptons"       ,&nLeptons);
    trIN->SetBranchAddress("jetPt"          ,&jetPt);
    trIN->SetBranchAddress("deepAK8Tagger"        ,&deepAK8);
    //trIN->SetBranchAddress("deepAk8"        ,&deepAK8);
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
    
    trIN->SetBranchAddress("mJJ"   			    ,&mJJ);
    trIN->SetBranchAddress("jetBtagSub0DCSVbb", &jetBtagSub0DCSVbb);
    trIN->SetBranchAddress("jetBtagSub0DCSVbbb", &jetBtagSub0DCSVbbb);
    trIN->SetBranchAddress("jetBtagSub1DCSVbb", &jetBtagSub1DCSVbb);
    trIN->SetBranchAddress("jetBtagSub1DCSVbbb", &jetBtagSub1DCSVbbb);
    trIN->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
    trIN->SetBranchAddress("ht", &ht);
    trIN->SetBranchAddress("jetTtagCategory",&jetTtag);
    trIN->SetBranchAddress("genEvtWeight", &eventWeight);
    
    float norm = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
	  float weight = utils.XSEC[f]/norm;
	  weights.push_back(weight);
    
    int decade(0);
    int NN = trIN->GetEntries();
    //NN = 100000;
    
    std::cout<<"Entries: "<<NN<<std::endl;
    for(int iev=0;iev<NN;iev++) 
    {
      double progress = 10.0*iev/(1.0*NN);
      int k = TMath::FloorNint(progress); 
      if (k > decade) 
        cout<<10*k<<" %"<<endl;
      decade = k;
      trIN->GetEntry(iev);
      
      if(nJets > 1 && fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1] <2.4) && (*jetPt)[0] > 400 && (*jetPt)[1] > 400 && mJJ > 1000 && nLeptons == 0/*&& (*triggerBit)[2]*/ )
	    {        
        h_denominator_reco[f]->Fill((*jetPt)[0], eventWeight);
        h_denominator_reco[f]->Fill((*jetPt)[1], eventWeight);
        
        //1. category where both are fully top tagged
          
		    if((*jetMassSoftDrop)[0] > 120 && (*jetMassSoftDrop)[0] < 220 && (*jetMassSoftDrop)[1] > 120 && (*jetMassSoftDrop)[1] < 220)
		    {
				/*
          jetTtag->clear();
          for(int i=0; i<nJets; i++)
          {
            if((*jetPt)[i] >= 400)
            {
              float jetPtOverSumPt = (*jetPt)[i]/ht;
              float score = tagger->eval((*tau3)[i], (*tau2)[i], (*tau1)[i], (*jetMassSub0)[i], (*jetMassSub1)[i], (*ecfB1N2)[i], (*ecfB1N3)[i], (*ecfB2N2)[i], (*ecfB2N3)[i], jetPtOverSumPt, (*jetPt)[i]);
            
              jetTtag->push_back(score);
            }
          }*/
	        if((*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut && 
             (((*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0])> floatBTag || 
              ((*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0])> floatBTag) && 
             (((*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1])> floatBTag || 
              ((*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1])> floatBTag))
		      {
            tTaggerEvents++;
            h_out_reco[0][f]->Fill((*jetPt)[0], eventWeight);
            h_out_reco[0][f]->Fill((*jetPt)[1], eventWeight);
		      }
          
          if((*deepAK8)[0] > deepAK8Cut && (*deepAK8)[1] > deepAK8Cut && 
             (((*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0])> floatBTag || 
              ((*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0])> floatBTag) && 
             (((*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1])> floatBTag || 
              ((*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1])> floatBTag))
          {
            deepAK8events++;
            h_out_reco[1][f]->Fill((*jetPt)[0], eventWeight);
            h_out_reco[1][f]->Fill((*jetPt)[1], eventWeight);
          }
          
          if(mva > mvaCut && 
             (((*jetBtagSub0DCSVbb)[0] + (*jetBtagSub0DCSVbbb)[0])> floatBTag || 
              ((*jetBtagSub1DCSVbb)[0] + (*jetBtagSub1DCSVbbb)[0])> floatBTag) && 
             (((*jetBtagSub0DCSVbb)[1] + (*jetBtagSub0DCSVbbb)[1])> floatBTag || 
              ((*jetBtagSub1DCSVbb)[1] + (*jetBtagSub1DCSVbbb)[1])> floatBTag))
          {
            h_out_reco[2][f]->Fill((*jetPt)[0], eventWeight);
            h_out_reco[2][f]->Fill((*jetPt)[1], eventWeight);
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
  std::cout<<h_out_reco[0][0]->GetEntries()<<std::endl;
  std::cout<<h_denominator_reco[0]->GetEntries()<<std::endl;
  std::cout<<h_out_reco[1][0]->GetEntries()<<std::endl;
  std::cout<<h_denominator_reco[0]->GetEntries()<<std::endl;
  std::cout<<h_out_reco[2][0]->GetEntries()<<std::endl;
  std::cout<<h_denominator_reco[0]->GetEntries()<<std::endl;
  efficiency_reco[0]  = new TEfficiency(*h_out_reco[0][0], *h_denominator_reco[0]);
  efficiency_reco[0]->SetTitle("reco_tTagger;jetPt (GeV);Efficiency");
  efficiency_reco[0]->SetStatisticOption(TEfficiency::kFNormal);
  efficiency_reco[0]->SetUseWeightedEvents();
  efficiency_reco[0]->SetLineColor(kRed);
  efficiency_reco[1]  = new TEfficiency(*h_out_reco[1][0], *h_denominator_reco[0]);
  efficiency_reco[1]->SetTitle("reco_deepAK8;jetPt (GeV);Efficiency");
  efficiency_reco[1]->SetStatisticOption(TEfficiency::kFNormal);
  efficiency_reco[1]->SetUseWeightedEvents();
  efficiency_reco[1]->SetLineColor(kBlue); 
  efficiency_reco[2]  = new TEfficiency(*h_out_reco[2][0], *h_denominator_reco[0]);
  efficiency_reco[2]->SetTitle("reco_eventMVA;jetPt (GeV);Efficiency");
  efficiency_reco[2]->SetStatisticOption(TEfficiency::kFNormal);
  efficiency_reco[2]->SetUseWeightedEvents();
  efficiency_reco[2]->SetLineColor(kBlack); 
  
  /*TCanvas *c2 = new TCanvas("recoCanvas", "recoCanvas", 600, 500);
  c2->cd();
  //efficiency_reco[0]->Draw();
  //efficiency_reco[1]->Draw("SAME");
  //efficiency_reco[2]->Draw("SAME");
  std::cout<<efficiency_reco[0]->GetLineColor()<<std::endl;
  std::cout<<efficiency_reco[1]->GetLineColor()<<std::endl;
  std::cout<<efficiency_reco[2]->GetLineColor()<<std::endl;
  */
  std::cout<<"Test"<<std::endl;
  TFile *outFile = TFile::Open("deepAK8_efficienciesJetPtOnlyRecoCutsBkg.root", "UPDATE");
  std::cout<<"Test"<<std::endl;
  if(!outFile->GetListOfKeys()->Contains(TString::Format("Bkg_Reco_jetPt_tTagger_%.2f", selMvaCut)))
  {
    efficiency_reco[0]->Write(TString::Format("Bkg_Reco_jetPt_tTagger_%.2f", selMvaCut));
  }
  std::cout<<"Test"<<std::endl;
  if(!outFile->GetListOfKeys()->Contains(TString::Format("Bkg_Reco_jetPt_eventMVA_%.2f", mvaCut)))
  {
    efficiency_reco[2]->Write(TString::Format("Bkg_Reco_jetPt_eventMVA_%.2f", mvaCut));
  }
  std::cout<<"Test"<<std::endl;
  if(!outFile->GetListOfKeys()->Contains(TString::Format("Bkg_Reco_jetPt_deepAK8_%.2f", deepAK8Cut)))
  {
    efficiency_reco[1]->Write(TString::Format("Bkg_Reco_jetPt_deepAK8_%.2f", deepAK8Cut));
  }
  std::cout<<"Test"<<std::endl;
}