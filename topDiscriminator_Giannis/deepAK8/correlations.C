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

void correlations(bool isSignal)
{
  gStyle->SetOptStat(0);
  initGlobals(isSignal);
  TString suffix = (isSignal?"Sig":"Bkg");
  //TH1F *h_out_reco[2][listOfFiles.size()];
  TH2F *h_out_reco[2];
  h_out_reco[0] = new TH2F(TString::Format("Correlation_plot_tTagger_%s", suffix.Data()), TString::Format("Correlation_plot_tTagger_%s", suffix.Data()), 20, -1, 1, 10, 0, 1);
  h_out_reco[1] = new TH2F(TString::Format("Correlation_plot_deepAK8_%s", suffix.Data()), TString::Format("Correlation_plot_deepAK8_%s", suffix.Data()), 10, 0, 1, 10, 0, 1);
  std::vector<float> weights;
  
  for(int f=0; f<listOfFiles.size(); f++)
  {
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
        if((*jetMassSoftDrop)[0] > 120 && (*jetMassSoftDrop)[0] < 220 && (*jetMassSoftDrop)[1] > 120 && (*jetMassSoftDrop)[1] < 220)
		    {
          for(int jet = 0; jet<2; jet++)
          {
            if((*jetBtagSub0)[jet] > (*jetBtagSub1)[jet])
            {
              h_out_reco[0]->Fill((*jetTtag)[jet], (*jetBtagSub0)[jet], weight);
              h_out_reco[1]->Fill((*deepAK8)[jet], (*jetBtagSub0)[jet], weight);
            }
            else
            {
              h_out_reco[0]->Fill((*jetTtag)[jet], (*jetBtagSub1)[jet], weight);
              h_out_reco[1]->Fill((*deepAK8)[jet], (*jetBtagSub1)[jet], weight);
            }
          }
        }
      }
    }
  }
  
  h_out_reco[0]->Scale(LUMI);
  h_out_reco[1]->Scale(LUMI);
  
  TCanvas *c1 = new TCanvas("tTagger vs CSVV2", "tTagger vs CSVV2", 600, 500);
  c1->cd();
  h_out_reco[0]->Scale(1./h_out_reco[0]->Integral());
  h_out_reco[0]->Draw("COLZ");
  
  TCanvas *c2 = new TCanvas("deepAK8 vs CSVV2", "deepAK8 vs CSVV2", 600, 500);
  c2->cd();
  h_out_reco[1]->Scale(1./h_out_reco[1]->Integral());
  h_out_reco[1]->Draw("COLZ");
  
  TFile *outFile = new TFile("correlationPlots.root", "UPDATE");
  outFile->cd();
  h_out_reco[0]->Write();
  h_out_reco[1]->Write();
  outFile->Close();
}