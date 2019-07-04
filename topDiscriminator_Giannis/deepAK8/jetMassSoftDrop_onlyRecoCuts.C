std::vector<TString> listOfFiles;
std::vector<float> XSEC;
std::vector<TString> histoNames;
float LUMI = 35900;
TString eosPath; 

void initFileNames(bool isSignal)
{
  if(isSignal)
  {
    eosPath = "/eos/cms/store/user/ipapakri/ttbar/MC/Signal/";  
    listOfFiles.push_back("TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8_legacy2016_deepAK8.root");
    listOfFiles.push_back("TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8_legacy2016_deepAK8.root");
  }
  else
  {
    eosPath = "/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/";
    listOfFiles.push_back("QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
    listOfFiles.push_back("QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
    listOfFiles.push_back("QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
    listOfFiles.push_back("QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
    listOfFiles.push_back("QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
    listOfFiles.push_back("QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
  }
}


void initXsections(bool isSignal)
{
  if(isSignal)
  {
    XSEC.push_back(69.64);
    XSEC.push_back(16.74);
  }
  else
  {
    XSEC.push_back(3.67e+5);
    XSEC.push_back(2.94e+4);
    XSEC.push_back(6.524e+03);
    XSEC.push_back(1.064e+03);
    XSEC.push_back(121.5);
    XSEC.push_back(2.542e+01);
  }
}

void initHistoNames(bool isSignal)
{
  if(isSignal)
  {
    histoNames.push_back("Signal_histo_Mtt_700_1000");
    histoNames.push_back("Signal_histo_Mtt_1000_Inf");
  }
  else
  {
    histoNames.push_back("QCD_histo_Mtt_300_500");
    histoNames.push_back("QCD_histo_Mtt_500_700");
    histoNames.push_back("QCD_histo_Mtt_700_1000");
    histoNames.push_back("QCD_histo_Mtt_1000_1500");
    histoNames.push_back("QCD_histo_Mtt_1500_2000");
    histoNames.push_back("QCD_histo_Mtt_2000_Inf");
  }
}

void initGlobals(bool isSignal)
{
  initFileNames(isSignal);
  initXsections(isSignal);
  initHistoNames(isSignal);
}

void jetMassSoftDrop_onlyRecoCuts(bool isSignal)
{
  gStyle->SetOptStat(0);
  initGlobals(isSignal);
  TString suffix = (isSignal?"Sig":"Bkg");
  //TH1F *h_out_reco[2][listOfFiles.size()];
  TH1F *h_out_reco[listOfFiles.size()];
  std::vector<float> weights;
  
  for(int f=0; f<listOfFiles.size(); f++)
  {
    h_out_reco[f] = new TH1F(TString::Format("jetMassSoftDrop_%s", suffix.Data()), TString::Format("jetMassSoftDrop_%s", suffix.Data()), 60, 50, 450);
    int nJets,nLeptons;
    std::vector<float> *jetPt(0), *jetEta(0), *jetPhi(0), *deepAK8(0);
    std::vector<float> *mass(0), *tau3(0), *tau2(0), *tau1(0);
    std::vector<float> *jetMassSub0(0), *jetMassSub1(0), *jetBtagSub0(0);
    std::vector<float> *jetBtagSub1(0), *jetTtag(0), *jetMassSoftDrop(0);
    
    float mTTbarParton(0), mJJ(0), mva(0);
    
    std::cout<<"Working in file: "<<listOfFiles[f]<<std::endl;
    TFile *file = TFile::Open(eosPath+listOfFiles[f]);
    TTree *eventTree = (TTree*)file->Get("boosted/events");
    
    eventTree->SetBranchAddress("nJets"          ,&nJets);
    eventTree->SetBranchAddress("nLeptons"       ,&nLeptons);
    eventTree->SetBranchAddress("jetPt"          ,&jetPt);
    eventTree->SetBranchAddress("deepAk8"        ,&deepAK8);
    eventTree->SetBranchAddress("mva"            ,&mva);
    
    eventTree->SetBranchAddress("jetEta"         ,&jetEta);
    eventTree->SetBranchAddress("jetPhi"         ,&jetPhi);
    eventTree->SetBranchAddress("jetMassSoftDrop",&mass);
    eventTree->SetBranchAddress("jetTau3"        ,&tau3);
    eventTree->SetBranchAddress("jetTau2"        ,&tau2);
    eventTree->SetBranchAddress("jetTau1"        ,&tau1);
    eventTree->SetBranchAddress("jetMassSub0"    ,&jetMassSub0);
    eventTree->SetBranchAddress("jetMassSub1"    ,&jetMassSub1);
	  eventTree->SetBranchAddress("mTTbarParton"   ,&mTTbarParton);
    eventTree->SetBranchAddress("mJJ"   			   ,&mJJ);
    eventTree->SetBranchAddress("jetBtagSub0"	   ,&jetBtagSub0);
    eventTree->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1);
    eventTree->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
    eventTree->SetBranchAddress("jetTtagCategory",&jetTtag);
    
    float norm = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
	  float weight = XSEC[f]/norm;
    weights.push_back(weight);
    
    int decade(0);
    int NN = eventTree->GetEntries();
    //NN = 100000;
    std::cout<<"Entries: "<<NN<<std::endl;
    
    for(int iev=0;iev<NN;iev++) 
    {
      double progress = 10.0*iev/(1.0*NN);
      int k = TMath::FloorNint(progress); 
      if (k > decade) 
        cout<<10*k<<" %"<<endl;
      decade = k;
      eventTree->GetEntry(iev);
      
      if(nJets > 1 && fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1] <2.4) && (*jetPt)[0] > 400 && (*jetPt)[1] > 400 && mJJ > 1000 /*&& (*triggerBit)[2]*/ )
      {
        h_out_reco[f]->Fill((*jetMassSoftDrop)[0]);
        h_out_reco[f]->Fill((*jetMassSoftDrop)[1]);
      }
    }
  }
  
  for(int i=0; i<sizeof(h_out_reco)/sizeof(*h_out_reco); i++)
  {
	//for every slice
    h_out_reco[i]->Scale(weights[i]);
  }
  
  for(int i=1; i<sizeof(h_out_reco)/sizeof(*h_out_reco); i++)
	{
	//Add them to get the whole phase space
    h_out_reco[0]->Add(h_out_reco[i]);
	}
  
  h_out_reco[0]->Scale(LUMI);
  
  TFile* f = new TFile("jetMassSoftDrop.root", "UPDATE");
  
  TCanvas *c1 = new TCanvas("c1", "c1", 600, 500);
  c1->cd();
  h_out_reco[0]->Draw();
  
  f->cd();
  h_out_reco[0]->Write();
  f->Close();
}