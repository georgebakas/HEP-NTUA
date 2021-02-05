void UpdateTopTagger(TString oldFileName)
{

    TString weightsFile = "/afs/cern.ch/work/g/gbakas/public/TrainingOutputs/UL_Training/Training_2018/BoostedMVA/weights/boosted_MVA_BDTCat.weights.xml";
    
    std::cout<<"Openning: "<<oldFileName<<std::endl;
    TFile *oldFile = TFile::Open(oldFileName+".root");

    TH1F* triggerNames = (TH1F*) oldFile->Get("boosted/TriggerNames");
    TH1F* cutFlow = (TH1F*) oldFile->Get("boosted/CutFlow");
    TH1F* triggerPass = (TH1F*) oldFile->Get("boosted/TriggerPass");
    TTree* tr = (TTree*) oldFile->Get("boosted/events");

    tr->SetBranchStatus("jetTtagCategory", 0);
    /*TBranch *b_temp = tr->GetBranch("jetTtagCategory");
    tr->GetListOfBranches()->Remove(b_temp);
    TLeaf *l = tr->GetLeaf("jetTtagCategory");
    tr->GetListOfLeaves()->Remove(l);
    cout<<"deleted previous branch!"<<endl; */

    TString newFileName = oldFileName + "_new.root";
    std::cout<<"Output file is: "<<newFileName<<std::endl;
    TFile *newFile = TFile::Open(newFileName, "RECREATE");
    if(oldFile->cd("eventCounter"))
    {
      TH1F* GenEventWeight = (TH1F*) oldFile->Get("eventCounter/GenEventWeight");
      TH1F* pileup = (TH1F*) oldFile->Get("eventCounter/pileup");
      TTree* events = (TTree*) oldFile->Get("eventCounter/events");

      newFile->mkdir("eventCounter");
      newFile->cd("eventCounter");

      GenEventWeight->Write("GenEventWeight");
      pileup->Write("pileup");

      TTree *newTree = events->CloneTree(-1, "fast");
      newTree->Write();
      newFile->cd();
    }
    newFile->mkdir("boosted");
    newFile->cd("boosted");

    triggerNames->Write("TriggerNames");
    cutFlow->Write("CutFlow");
    triggerPass->Write("TriggerPass");

    //
    TTree *outputTree = tr->CloneTree(-1, "fast");

    std::vector<float> *jetPt = new std::vector<float>();
    std::vector<float> *tau1  = new std::vector<float>();
    std::vector<float> *tau2  = new std::vector<float>();
    std::vector<float> *tau3  = new std::vector<float>();
    std::vector<float> *jetMassSub1  = new std::vector<float>();
    std::vector<float> *jetMassSub0  = new std::vector<float>();
    std::vector<float> *ecfB1N2  = new std::vector<float>();
    std::vector<float> *ecfB1N3  = new std::vector<float>();
    std::vector<float> *ecfB2N2  = new std::vector<float>();
    std::vector<float> *ecfB2N3  = new std::vector<float>();
    std::vector<float> *jetTtagCategory = new std::vector<float>();
    float ht;
    int nJets;
    //TTree *tr = (TTree*) f->Get("boosted/events");
    outputTree->SetBranchAddress("jetPt", &jetPt);
    outputTree->SetBranchAddress("jetTau1", &tau1);
    outputTree->SetBranchAddress("jetTau2", &tau2);
    outputTree->SetBranchAddress("jetTau3", &tau3);
    outputTree->SetBranchAddress("jetTau3", &tau3);
    outputTree->SetBranchAddress("jetMassSub0", &jetMassSub0);
    outputTree->SetBranchAddress("jetMassSub1", &jetMassSub1);
    outputTree->SetBranchAddress("ecfB1N2", &ecfB1N2);
    outputTree->SetBranchAddress("ecfB1N3", &ecfB1N3);
    outputTree->SetBranchAddress("ecfB2N2", &ecfB2N2);
    outputTree->SetBranchAddress("ecfB2N3", &ecfB2N3);
    outputTree->SetBranchAddress("ht", &ht);
    outputTree->SetBranchAddress("nJets", &nJets);


    TBranch *bpt = outputTree->Branch("jetTtagCategory", "vector<float>", &jetTtagCategory);

    TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");
    float var[11];

    reader->AddVariable("jetTau3",&var[0]);
    reader->AddVariable("jetTau2",&var[1]);
    reader->AddVariable("jetTau1",&var[2]);
    reader->AddVariable("jetMassSub0",&var[3]);
    reader->AddVariable("jetMassSub1",&var[4]);
    reader->AddVariable("ecfB1N2",&var[5]);
    reader->AddVariable("ecfB1N3",&var[6]);
    reader->AddVariable("ecfB2N2",&var[7]);
    reader->AddVariable("ecfB2N3",&var[8]);
    reader->AddVariable("JetPtOverSumPt",&var[9]);

    reader->AddSpectator("jetPt",&var[10]);

    reader->BookMVA("BDTCat", weightsFile);

    for(int evt = 0; evt< tr->GetEntries(); evt++)
    {
      outputTree->GetEntry(evt);

      for(int ijet = 0; ijet<nJets; ijet++)
      {
        if((*jetPt)[ijet] > 400)
        {
          var[0] = (*tau3)[ijet];
          var[1] = (*tau2)[ijet];
          var[2] = (*tau1)[ijet];
          var[3] = (*jetMassSub0)[ijet];
          var[4] = (*jetMassSub1)[ijet];
          var[5] = (*ecfB1N2)[ijet];
          var[6] = (*ecfB1N3)[ijet];
          var[7] = (*ecfB2N2)[ijet];
          var[8] = (*ecfB2N3)[ijet];
          var[9] = (*jetPt)[ijet]/ht;
          var[10] = (*jetPt)[ijet];

          jetTtagCategory->push_back(reader->EvaluateMVA("BDTCat"));
        }
        else
        {
          jetTtagCategory->push_back(-10);
        }
      }

      bpt->Fill();
      jetTtagCategory->clear();
    }
    newFile->cd("boosted");
    //outputTree->Write();

    newFile->Write();
    oldFile->Close();
    newFile->Close();
}
