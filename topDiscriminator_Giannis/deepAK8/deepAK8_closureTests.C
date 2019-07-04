std::vector<TString> listOfFiles;
std::vector<float> XSEC;
std::vector<TString> histoNames;
float LUMI = 35900;
TString eosPath; 

void initFileNames()
{
  eosPath = "/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/";
  listOfFiles.push_back("QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
  listOfFiles.push_back("QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
  listOfFiles.push_back("QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
  listOfFiles.push_back("QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
  listOfFiles.push_back("QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
  listOfFiles.push_back("QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_deepAK8.root");
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

/*
bool taggerCuts(float mass, std:vector<float> topTaggerScores, float topTaggerCut, std::vector<float> bTaggerScores, float bTaggerCut)
{
  if(mass > 1000)
  {
    if()
  }
  else
  {
    return false;
  }
}
*/

void deepAK8_closureTests(bool saveTtagger = false, float deepAK8Cut = 0.6, float selMvaCut = 0.3, float mvaCut = 0.8)
{
  initGlobals();
  
  TH1F *h_out_reco[3][listOfFiles.size()];
  TH1F *h_out_recoCR[3][listOfFiles.size()];
  
  TH1F* h_denominator_reco[listOfFiles.size()];
  std::vector<float> weights;
  
  int sizeBins = 10;
  float BND[11] = {1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, 4000, 5000};
  int MIN = 1000; 
  int MAX = 5000;
  int deepAK8events = 0;
  int tTaggerEvents = 0;
  for(int f=0; f<listOfFiles.size(); f++)
  {
    /*
    h_out_reco[0][f] = new TH1F(TString::Format("%s_%s_reco", histoNames[f].Data(), "tTagger"), TString::Format("%s_%s_reco", histoNames[f].Data(), "tTagger"), sizeBins, BND);
    h_out_reco[1][f] = new TH1F(TString::Format("%s_%s_reco", histoNames[f].Data(), "deepAK8"), TString::Format("%s_%s_reco", histoNames[f].Data(), "deepAK8"), sizeBins, BND);
    h_out_reco[2][f] = new TH1F(TString::Format("%s_%s_reco", histoNames[f].Data(), "eventMVA"), TString::Format("%s_%s_reco", histoNames[f].Data(), "eventMVA"), sizeBins, BND);
    */
    h_out_recoCR[0][f] = new TH1F(TString::Format("%s_%s_recoCR", histoNames[f].Data(), "tTagger"), TString::Format("%s_%s_recoCR", histoNames[f].Data(), "tTagger"), sizeBins, BND);
    h_out_recoCR[1][f] = new TH1F(TString::Format("%s_%s_recoCR", histoNames[f].Data(), "deepAK8"), TString::Format("%s_%s_recoCR", histoNames[f].Data(), "deepAK8"), sizeBins, BND);
    h_out_recoCR[2][f] = new TH1F(TString::Format("%s_%s_recoCR", histoNames[f].Data(), "eventMVA"), TString::Format("%s_%s_recoCR", histoNames[f].Data(), "eventMVA"), sizeBins, BND);

    h_denominator_reco[f] = new TH1F(TString::Format("denom_%s_reco", histoNames[f].Data()), TString::Format("denom_%s_reco", histoNames[f].Data()), sizeBins, BND);
	  h_denominator_reco[f]->Sumw2();
    
    int nJets,nLeptons;
    std::vector<float> *jetPt(0), *jetEta(0), *jetPhi(0), *deepAK8(0);
    std::vector<float> *mass(0), *tau3(0), *tau2(0), *tau1(0);
    std::vector<float> *jetMassSub0(0), *jetMassSub1(0), *jetBtagSub0(0);
    std::vector<float> *jetBtagSub1(0), *jetTtag(0), *jetMassSoftDrop(0);
    
    float mTTbarParton(0), mJJ(0), mva(0);
    
    std::cout<<"Working in file: "<<listOfFiles[f]<<std::endl;
    TFile *file = TFile::Open(eosPath+listOfFiles[f]);
    TTree *trIN = (TTree*)file->Get("boosted/events");
    
    trIN->SetBranchAddress("nJets"          ,&nJets);
    trIN->SetBranchAddress("nLeptons"       ,&nLeptons);
    trIN->SetBranchAddress("jetPt"          ,&jetPt);
    trIN->SetBranchAddress("deepAk8"        ,&deepAK8);
    trIN->SetBranchAddress("mva"            ,&mva);
    
    trIN->SetBranchAddress("jetEta"         ,&jetEta);
    trIN->SetBranchAddress("jetPhi"         ,&jetPhi);
    trIN->SetBranchAddress("jetMassSoftDrop",&mass);
    trIN->SetBranchAddress("jetTau3"        ,&tau3);
    trIN->SetBranchAddress("jetTau2"        ,&tau2);
    trIN->SetBranchAddress("jetTau1"        ,&tau1);
    trIN->SetBranchAddress("jetMassSub0"    ,&jetMassSub0);
    trIN->SetBranchAddress("jetMassSub1"    ,&jetMassSub1);
	//trIN->SetBranchAddress("mTTbarParton"   ,&mTTbarParton);
    trIN->SetBranchAddress("mJJ"   			    ,&mJJ);
    trIN->SetBranchAddress("jetBtagSub0"	  ,&jetBtagSub0);
    trIN->SetBranchAddress("jetBtagSub1"    ,&jetBtagSub1);
    trIN->SetBranchAddress("jetMassSoftDrop",&jetMassSoftDrop);
    trIN->SetBranchAddress("jetTtagCategory",&jetTtag);
    
    float norm = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
	  float weight = XSEC[f]/norm;
	  weights.push_back(weight);
    
    int decade(0);
    int NN = trIN->GetEntries();
    //NN = 10000;
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
      
      if(nJets > 1 && fabs((*jetEta)[0]) < 2.4 && fabs((*jetEta)[1] <2.4) && (*jetPt)[0] > 400 && (*jetPt)[1] > 400 && mJJ > 1000 /*&& (*triggerBit)[2]*/ )
	    {
        //h_denominator_reco[f]->Fill(mJJ);
        
        //1. category where both are fully top tagged
          
		    if((*jetMassSoftDrop)[0] > 120 && (*jetMassSoftDrop)[0] < 220 && (*jetMassSoftDrop)[1] > 120 && (*jetMassSoftDrop)[1] < 220)
		    {
	        if((*jetTtag)[0] > selMvaCut && (*jetTtag)[1] > selMvaCut && ((*jetBtagSub0)[0] < floatBTag && (*jetBtagSub1)[0] < floatBTag)
            && ((*jetBtagSub0)[1] < floatBTag && (*jetBtagSub1)[1] < floatBTag))
		      {
            tTaggerEvents++;
            h_out_recoCR[0][f]->Fill(mJJ);
		      }
          
          if((*deepAK8)[0] > deepAK8Cut && (*deepAK8)[1] > deepAK8Cut && ((*jetBtagSub0)[0] < floatBTag && (*jetBtagSub1)[0] < floatBTag)
            && ((*jetBtagSub0)[1] < floatBTag && (*jetBtagSub1)[1] < floatBTag))
          {
            deepAK8events++;
            h_out_recoCR[1][f]->Fill(mJJ);
          }
          
          if(mva > mvaCut && ((*jetBtagSub0)[0] < floatBTag && (*jetBtagSub1)[0] < floatBTag)
            && ((*jetBtagSub0)[1] < floatBTag && (*jetBtagSub1)[1] < floatBTag))
          {
            h_out_recoCR[2][f]->Fill(mJJ);
          }
        }
      }
    }
    file->Close();
  }
  
  std::cout<<"Events deepAK8: "<<deepAK8events<<std::endl;
  std::cout<<"Events topTagger: "<<tTaggerEvents<<std::endl;
  
  for(int i=0; i<sizeof(h_out_recoCR)/sizeof(*h_out_recoCR); i++)
  {
	//for every slice
  	for(int j=0; j<sizeof(*h_out_recoCR)/sizeof(*h_out_recoCR[0]); j++)
    {
      h_out_recoCR[i][j]->Scale(weights[j]);
    }
    
    for(int j=1; j<sizeof(*h_out_recoCR)/sizeof(*h_out_recoCR[0]); j++)
	  {
	    h_out_recoCR[i][0]->Add(h_out_recoCR[i][j]);
	  }
  }
  
  std::cout<<"Test"<<std::endl;
  
  TCanvas *c1 = new TCanvas("c1", "c1", 600, 500);
  h_out_recoCR[0][0]->SetMarkerStyle(20);
  h_out_recoCR[0][0]->SetMarkerColor(kRed);
  h_out_recoCR[0][0]->SetLineColor(kRed);
  h_out_recoCR[0][0]->Draw("P");
  h_out_recoCR[1][0]->SetMarkerStyle(21);
  h_out_recoCR[1][0]->SetMarkerColor(kBlue);
  h_out_recoCR[1][0]->SetLineColor(kBlue);
  h_out_recoCR[1][0]->Draw("PSAME");
  h_out_recoCR[2][0]->SetMarkerStyle(22);
  h_out_recoCR[2][0]->SetMarkerColor(kMagenta);
  h_out_recoCR[2][0]->SetLineColor(kMagenta);
  h_out_recoCR[2][0]->Draw("PSAME");

  TFile *outFile = TFile::Open("deepAK8_closureTestsOnlyRecoCutsBkg.root", "UPDATE");
  
  if(saveTtagger)
  {
    h_out_recoCR[0][0]->Write(TString::Format("Bkg_RecoCR_tTagger_%.1f", selMvaCut));
    h_out_recoCR[2][0]->Write(TString::Format("Bkg_RecoCR_eventMVA_%.1f", mvaCut));
  }
  h_out_recoCR[1][0]->Write(TString::Format("Bkg_RecoCR_deepAK8_%.1f", deepAK8Cut));
  
}