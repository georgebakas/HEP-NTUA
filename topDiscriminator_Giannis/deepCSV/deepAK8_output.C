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

std::vector<TString> getFileNames()
{
  std::vector<TString> fileNames;
  fileNames.push_back("/eos/cms/store/user/ipapakri/ttbar/MC/Signal/TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8_legacy2016_deepAK8.root");
  fileNames.push_back("/eos/cms/store/user/ipapakri/ttbar/MC/Signal/TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8_legacy2016_deepAK8.root");
  
  return fileNames;
}

void deepAK8_output(bool isSignal = true)
{
  initGlobals(isSignal);
  
  for(std::vector<TString>::iterator it=listOfFiles.begin(); it!=listOfFiles.end(); ++it)
  {
    std::cout<<*it<<std::endl;
  }
  
  for(std::vector<TString>::iterator it=histoNames.begin(); it!=histoNames.end(); ++it)
  {
    std::cout<<*it<<std::endl;
  }
  
  for(std::vector<float>::iterator it=XSEC.begin(); it!=XSEC.end(); ++it)
  {
    std::cout<<*it<<std::endl;
  }
  
  TH1F *h_out[2][listOfFiles.size()];
  std::vector<float> weights;
  
  for(int f=0; f<listOfFiles.size(); f++)
  {
      
    h_out[0][f] = new TH1F(TString::Format("%s_tTagger", histoNames[f].Data()), TString::Format("%s_tTagger", histoNames[f].Data()), 100, -1, 1);
    h_out[1][f] = new TH1F(TString::Format("%s_deepAK8", histoNames[f].Data()), TString::Format("%s_deepAK8", histoNames[f].Data()), 50, 0, 1);
    
    int nJets,nLeptons;
    vector<float> *jetPt(0), *deepAK8(0), *jetTtag(0);
    
    std::cout<<"Working in file: "<<listOfFiles[f]<<std::endl;
    TFile *file = TFile::Open(eosPath+listOfFiles[f]);
    TTree *trIN = (TTree*)file->Get("boosted/events");
    
    trIN->SetBranchAddress("nJets"          ,&nJets);
    trIN->SetBranchAddress("nLeptons"       ,&nLeptons);
    trIN->SetBranchAddress("jetPt"          ,&jetPt);
    trIN->SetBranchAddress("deepAk8"        ,&deepAK8);
    trIN->SetBranchAddress("jetTtagCategory",&jetTtag);
    
    float norm = ((TH1F*)file->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
	float weight = XSEC[f]/norm;
	weights.push_back(weight);
    
    int decade(0);
    int entries = trIN->GetEntries();
    //entries = 10000;
    for(int i=0; i<entries; i++)
    {
      double progress = 10.0*i/(1.0*entries);
      int k = TMath::FloorNint(progress); 
      if (k > decade) 
        cout<<10*k<<" %"<<endl;
      decade = k;
      
      trIN->GetEntry(i);
      if (nJets < 2) continue;
      if (nLeptons > 0) continue;
      if ((*jetPt)[1] < 400) continue;
      
      for(int j=0; j<nJets; j++)
      {
        h_out[0][f]->Fill((*jetTtag)[j]);
        h_out[1][f]->Fill((*deepAK8)[j]);
      }
    }
  }
  
  TFile *outFile = TFile::Open("deepAK8_outPut.root", "UPDATE");
  std::cout<<"Size: "<<sizeof(h_out)/sizeof(*h_out)<<std::endl;
  for(int i=0; i<sizeof(h_out)/sizeof(*h_out); i++)
  {
    h_out[0][i]->Scale(weights[i]);
    h_out[1][i]->Scale(weights[i]);
  }
  
  for(int i=1; i<sizeof(*h_out)/sizeof(*h_out[0]); i++)
  {
    h_out[0][0]->Add(h_out[0][i]);
    h_out[1][0]->Add(h_out[1][i]);
  }
  
  if(isSignal)
  {
    h_out[0][0]->Scale(1./h_out[0][0]->Integral(), "width");
    h_out[0][0]->Write("tTagger_output_Signal");
    h_out[1][0]->Scale(1./h_out[1][0]->Integral(), "width");
    h_out[1][0]->Write("deepAK8_output_Signal");
  }
  else
  {
    h_out[0][0]->Scale(1./h_out[0][0]->Integral(), "width");
    h_out[0][0]->Write("tTagger_output_Bkg");
    h_out[1][0]->Scale(1./h_out[1][0]->Integral(), "width");
    h_out[1][0]->Write("tTagger_output_Bkg");
  }
  
  TCanvas *c1 = new TCanvas("c1", "c1", 600, 500);
  
  h_out[0][0]->Draw();
  h_out[1][0]->SetLineColor(kRed);
  h_out[1][0]->Draw("SAME");
  outFile->Close();
}

/*
void deepAK8_output_experimental()
{
  std::vector<TString> fileNames = getFileNames();
  std::vector<TH1F*> histograms;
  std::cout<<"Size: "<<fileNames.size()<<std::endl;
  int i=1;
  for(std::vector<TString>::iterator fileName = fileNames.begin(); fileName != fileNames.end(); ++fileName)
  {
    std::cout<<"FileName: "<<*fileName<<std::endl;
    TFile *inputFile = TFile::Open(*fileName);
    std::cout<<"Name: "<<inputFile->GetName()<<std::endl;
    
    int nJets, nLeptons;
    std::vector<float> *jetPt, *deepAK8;
  
    trIN->SetBranchAddress("nJets"          ,&nJets);
    trIN->SetBranchAddress("nLeptons"       ,&nLeptons);
    trIN->SetBranchAddress("jetPt"          ,&jetPt);
    trIN->SetBranchAddress("deepAK8"         ,&deepAK8);
    
    TH1F* histo = new TH1F(TString::Format("Histo_File_%i", i), TString::Format("Histo_File_%i", i), 100, 0, 1);
    histo->Fill(0.5);
    histograms.push_back(histo);
    
    inputFile->Close();
    i++;
  }
  
  std::cout<<"Size: "<<histograms.size()<<std::endl;
  for(std::vector<TH1F*>::iterator histogram = histograms.begin(); histogram != histograms.end(); ++histogram)
  {
    std::cout<<typeid(*histogram).name()<<std::endl;
    std::cout<<"Name: "<<(*histogram).GetName()<<" entries: "<<(*histogram).GetEntries()<<std::endl;
  }
}
*/