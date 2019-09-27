std::vector<TString> listOfFiles;
std::vector<float> XSEC;
std::vector<TString> histoNames;
//float LUMI = 35900;
float LUMI = 41530;
TString eosPath; 

void initFileNames()
{
  eosPath = "/eos/cms/store/user/ipapakri/ttbar/MC/Signal/2017/";  
  listOfFiles.push_back("TT_Mtt-700to1000_TuneCP5_13TeV-powheg-pythia8.root");
  listOfFiles.push_back("TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root");
}

void initXsections()
{
  XSEC.push_back(69.64);
  XSEC.push_back(16.74);
}

void initHistoNames()
{
  histoNames.push_back("Signal_histo_Mtt_700_1000");
  histoNames.push_back("Signal_histo_Mtt_1000_Inf");
}

void initGlobals()
{
  initFileNames();
  initXsections();
  initHistoNames();
}

void deepAK8_efficiencies(bool saveTtagger = false, float deepAK8Cut = 0.6, float selMvaCut = 0.3)
{
  initGlobals();
  
  TH1F *h_out_parton[2][listOfFiles.size()];
  TH1F *h_out_reco[2][listOfFiles.size()];
  
  TH1F* h_denominator_parton[listOfFiles.size()];
  TH1F* h_denominator_reco[listOfFiles.size()];
  std::vector<float> weights;
  
  int sizeBins = 10;
  float BND[11] = {1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, 4000, 5000};
  int MIN = 1000; 
  int MAX = 5000;
  
  for(int f=0; f<listOfFiles.size(); f++)
  {
    h_out_parton[0][f] = new TH1F(TString::Format("%s_%s_parton", histoNames[f].Data(), "tTagger"), TString::Format("%s_%s_parton", histoNames[f].Data(), "tTagger"), sizeBins, BND);
    h_out_parton[1][f] = new TH1F(TString::Format("%s_%s_parton", histoNames[f].Data(), "deepAK8"), TString::Format("%s_%s_parton", histoNames[f].Data(), "deepAK8"), sizeBins, BND);
    
    h_out_reco[0][f] = new TH1F(TString::Format("%s_%s_reco", histoNames[f].Data(), "tTagger"), TString::Format("%s_%s_reco", histoNames[f].Data(), "tTagger"), sizeBins, BND);
    h_out_reco[1][f] = new TH1F(TString::Format("%s_%s_reco", histoNames[f].Data(), "deepAK8"), TString::Format("%s_%s_reco", histoNames[f].Data(), "deepAK8"), sizeBins, BND);
    
    h_denominator_parton[f] = new TH1F(TString::Format("denom_%s_parton", histoNames[f].Data()), TString::Format("denom_%s_parton", histoNames[f].Data()), sizeBins, BND);
	  h_denominator_parton[f]->Sumw2();
    
    h_denominator_reco[f] = new TH1F(TString::Format("denom_%s_reco", histoNames[f].Data()), TString::Format("denom_%s_reco", histoNames[f].Data()), sizeBins, BND);
	  h_denominator_reco[f]->Sumw2();
    
    int nJets,nLeptons;
    std::vector<float> *jetPt(0), *jetEta(0), *jetPhi(0), *deepAK8(0);
    std::vector<float> *mass(0), *tau3(0), *tau2(0), *tau1(0);
    std::vector<float> *jetMassSub0(0), *jetMassSub1(0), *jetBtagSub0(0);
    std::vector<float> *jetBtagSub1(0), *jetTtag(0), *jetMassSoftDrop(0);
    
    float mTTbarParton(0), mJJ(0);
    
    std::cout<<"Working in file: "<<listOfFiles[f]<<std::endl;
    TFile *file = TFile::Open(eosPath+listOfFiles[f]);
    TTree *trIN = (TTree*)file->Get("boosted/events");
    
    trIN->SetBranchAddress("nJets"          ,&nJets);
    trIN->SetBranchAddress("nLeptons"       ,&nLeptons);
    trIN->SetBranchAddress("jetPt"          ,&jetPt);
    trIN->SetBranchAddress("deepAk8"        ,&deepAK8);
    
    trIN->SetBranchAddress("jetEta"         ,&jetEta);
    trIN->SetBranchAddress("jetPhi"         ,&jetPhi);
    trIN->SetBranchAddress("jetMassSoftDrop",&mass);
    trIN->SetBranchAddress("jetTau3"        ,&tau3);
    trIN->SetBranchAddress("jetTau2"        ,&tau2);
    trIN->SetBranchAddress("jetTau1"        ,&tau1);
    trIN->SetBranchAddress("jetMassSub0"    ,&jetMassSub0);
    trIN->SetBranchAddress("jetMassSub1"    ,&jetMassSub1);
	  trIN->SetBranchAddress("mTTbarParton"   ,&mTTbarParton);
    trIN->SetBranchAddress("mJJ"   			,&mJJ);
    trIN->SetBranchAddress("jetBtagSub0"	,&jetBtagSub0);
    trIN->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1);
    trIN->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
    trIN->SetBranchAddress("jetTtagCategory",&jetTtag);
    
    float norm = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
	float weight = XSEC[f]/norm;
	weights.push_back(weight);
    
    int decade(0);
    int NN = trIN->GetEntries();
    //NN = 100000;
    float floatBTag = 0.8838;
    std::cout<<"Entries: "<<NN<<std::endl;
    for(int iev=0;iev<NN;iev++) 
    {
      double progress = 10.0*iev/(1.0*NN);
      int k = TMath::FloorNint(progress); 
      if (k > decade) 
        cout<<10*k<<" %"<<endl;
      decade = k;
      trIN->GetEntry(iev);
      
      if(nJets > 1 && fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1] <2.4) && (*jetPt)[0] > 400 && (*jetPt)[1] > 400 && mTTbarParton > 1000 /*&& (*triggerBit)[2]*/ )
	  {
          
        h_denominator_parton[f]->Fill(mTTbarParton);
        h_denominator_reco[f]->Fill(mJJ);
        //1. category where both are fully top tagged
          
		if((*jetMassSoftDrop)[0] > 120 && (*jetMassSoftDrop)[0] < 220 && (*jetMassSoftDrop)[1] > 120 && (*jetMassSoftDrop)[1] < 220)
		{
	      if((*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut && ((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag)
		    && ((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag))
		  {
			h_out_parton[0][f]->Fill(mTTbarParton);
            h_out_reco[0][f]->Fill(mJJ);
		  }
          
          if((*deepAK8)[0] > deepAK8Cut && (*deepAK8)[1] > deepAK8Cut &&  ((*jetBtagSub0)[0] > floatBTag || (*jetBtagSub1)[0] > floatBTag)
            && ((*jetBtagSub0)[1] > floatBTag || (*jetBtagSub1)[1] > floatBTag))
          {
            h_out_parton[1][f]->Fill(mTTbarParton);
            h_out_reco[1][f]->Fill(mJJ);
          }
        }
      }
    }
    file->Close();
  }
  
  for(int i=0; i<sizeof(h_out_parton)/sizeof(*h_out_parton); i++)
  {
	//for every slice
	for(int j=0; j<sizeof(*h_out_parton)/sizeof(*h_out_parton[0]); j++)
    {
      h_out_parton[i][j]->Scale(weights[j]);
      h_out_reco[i][j]->Scale(weights[j]);
    }
    
    for(int j=1; j<sizeof(*h_out_parton)/sizeof(*h_out_parton[0]); j++)
	{
	  //Add them to get the whole phase space
	  h_out_parton[i][0]->Add(h_out_parton[i][j]);
      std::cout<<h_out_reco[i][0]->GetName()<<" "<<h_out_reco[i][j]->GetName()<<std::endl;
      h_out_reco[i][0]->Add(h_out_reco[i][j]);
	}
  }
  
  h_denominator_parton[0]->Scale(weights[0]);
  h_denominator_reco[0]->Scale(weights[0]);
  
  for(int i=1; i<sizeof(h_denominator_parton)/sizeof(*h_denominator_parton); i++) 
  {
	h_denominator_parton[i]->Scale(weights[i]);
	h_denominator_parton[0]->Add(h_denominator_parton[i]);
    h_denominator_reco[i]->Scale(weights[i]);
	h_denominator_reco[0]->Add(h_denominator_reco[i]);
  }
  
  TEfficiency *efficiency_parton[2], *efficiency_reco[2];
  
  efficiency_parton[0]  = new TEfficiency(*h_out_parton[0][0], *h_denominator_parton[0]);
  efficiency_parton[0]->SetTitle("Parton_tTagger;mTTbarParton (GeV);Efficiency");
  efficiency_parton[0]->SetStatisticOption(TEfficiency::kFNormal);
  efficiency_parton[0]->SetUseWeightedEvents();
  efficiency_parton[0]->SetLineColor(kRed);
  
  efficiency_parton[1]  = new TEfficiency(*h_out_parton[1][0], *h_denominator_parton[0]);
  efficiency_parton[1]->SetTitle("Parton_deepAK8;mTTbarParton (GeV);Efficiency");
  efficiency_parton[1]->SetStatisticOption(TEfficiency::kFNormal);
  efficiency_parton[1]->SetUseWeightedEvents();
  efficiency_parton[1]->SetLineColor(kBlue);
  
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
  
  TCanvas *c1 = new TCanvas("partonCanvas", "partonCanvas", 600, 500);
  c1->cd();
  efficiency_parton[0]->Draw();
  efficiency_parton[1]->Draw("SAME");
  
  TCanvas *c2 = new TCanvas("recoCanvas", "recoCanvas", 600, 500);
  c2->cd();
  efficiency_reco[0]->Draw();
  efficiency_reco[1]->Draw("SAME");
  
  TFile *outFile = TFile::Open("deepAK8_efficiencies.root", "UPDATE");
  
  if(saveTtagger)
  {
    efficiency_parton[0]->Write(TString::Format("Sig_Parton_tTagger_%.1f", selMvaCut));
    efficiency_reco[0]->Write(TString::Format("Sig_Reco_tTagger_%.1f", selMvaCut));
  }
  efficiency_parton[1]->Write(TString::Format("Sig_Parton_deepAK8_%.1f", deepAK8Cut));
  efficiency_reco[1]->Write(TString::Format("Sig_Reco_deepAK8_%.1f", deepAK8Cut));
  
}