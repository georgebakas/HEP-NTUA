void plotVariable_distributions()
{
    TFile *f = TFile::Open("/afs/cern.ch/work/i/ipapakri/private/analysis/train/topJetDiscriminator/boosted_MVA.root");

    TTree *trainTree = (TTree*)f->Get("BoostedMVA/TestTree");

    std::cout<<trainTree->GetEntries()<<std::endl;
    //trainTree->Print();
    int id;
    int signalId = 0, bkgId = 1;
    float jetMassSub0;
    float jetPt;

    int *classId;
    //jetMassSub0
    TH1F* jetMassSub0_bkg = new TH1F("jetMassSub0_bkg", "jetMassSub0_bkg", 50, 0, 300);
    TH1F* jetMassSub0_Signal = new TH1F("jetMassSub0_Signal", "jetMassSub0_Signal", 30, 0, 300);

    //jetMassSub1
    TH1F* jetMassSub1_bkg = new TH1F("jetMassSub1_bkg", "jetMassSub1_bkg", 50, 0, 300);
    TH1F* jetMassSub1_Signal = new TH1F("jetMassSub1_Signal", "jetMassSub1_Signal", 30, 0, 300);

    //tau3
    TH1F* tau3_bkg = new TH1F("tau3_bkg", "tau3_bkg", 50, 0, 300);
    TH1F* tau3_Signal = new TH1F("tau3_Signal", "tau3_Signal", 30, 0, 300);

    //tau2
    TH1F* tau2_bkg = new TH1F("tau2_bkg", "tau2_bkg", 50, 0, 300);
    TH1F* tau2_Signal = new TH1F("tau2_Signal", "tau2_Signal", 30, 0, 300);

    //tau1
    TH1F* tau1_bkg = new TH1F("tau1_bkg", "tau1_bkg", 50, 0, 300);
    TH1F* tau1_Signal = new TH1F("tau1_Signal", "tau1_Signal", 30, 0, 300);

    trainTree->SetBranchAddress("classID", &id);
    trainTree->SetBranchAddress("jetMassSub0", &jetMassSub0);
    trainTree->SetBranchAddress("jetMassSub1", &jetMassSub1);
    trainTree->SetBranchAddress("tau3", &tau3);
    trainTree->SetBranchAddress("tau2", &tau2);
    trainTree->SetBranchAddress("tau1", &tau1);
    trainTree->SetBranchAddress("jetPt", &jetPt);

    for(int i=0; i<trainTree->GetEntries(); i++)
    {
      trainTree->GetEntry(i);
      if(jetPt>=1200)
      {
        if(id == signalId)
        {
            jetMassSub0_Signal->Fill(jetMassSub0);
        }
        else if(id == bkgId)
        {
          jetMassSub0_bkg->Fill(jetMassSub0);
        }
      }
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 500);
    jetMassSub0_bkg->SetLineColor(kRed);
    jetMassSub0_bkg->Scale(1./jetMassSub0_bkg->Integral());
    jetMassSub0_Signal->SetLineColor(kBlue);
    jetMassSub0_Signal->Scale(1./jetMassSub0_Signal->Integral());
    jetMassSub0_Signal->Draw("HIST");
    
    jetMassSub0_bkg->Draw("HISTSAME");
    std::cout<<"Signal entries: "<<jetMassSub0_Signal->GetEntries()<<std::endl;
    std::cout<<"Background entries: "<<jetMassSub0_bkg->GetEntries()<<std::endl;
}