R__LOAD_LIBRARY(BoostedTopCategoryDiscriminator_cc)
#include "BoostedTopCategoryDiscriminator.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>

void evaluateTopMVA_category()
{
  
  TString inputFile = "/eos/cms/store/user/ipapakri/ttbar/MC/Signal/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_legacy2016_Copy.root";
  
  //TString weightsCategory = "/afs/cern.ch/work/i/ipapakri/private/analysis/train/topJetDiscriminator/BoostedMVA/weights/boosted_MVA_BDTCat.weights.xml";
  TString weightsCategory = "/afs/cern.ch/user/g/gbakas/public/ForGiannis/Training/BoostedMVA/weights/boosted_MVA_BDTCat.weights.xml";
  TFile *f = new TFile(inputFile, "update");
  BoostedTopCategoryDiscriminator *topDiscriminator = new BoostedTopCategoryDiscriminator(weightsCategory);
  TTree* signalTree = (TTree*) f->Get("boosted/events");

  std::vector<float> *jetPt = 0, *jetTau1 = 0, *jetTau2 = 0, *jetTau3 = 0;
  std::vector<float> *jetMassSub0 = 0, *jetMassSub1 = 0;
  std::vector<float> *jetTtag = new std::vector<float>;
  std::vector<bool> *matchedJet = 0;
  int nJets;
  
  signalTree->SetBranchAddress("jetPt", &jetPt);
  signalTree->SetBranchAddress("jetTau1", &jetTau1);
  signalTree->SetBranchAddress("jetTau2", &jetTau2);
  signalTree->SetBranchAddress("jetTau3", &jetTau3);
  signalTree->SetBranchAddress("jetMassSub0", &jetMassSub0);
  signalTree->SetBranchAddress("jetMassSub1", &jetMassSub1);
  signalTree->SetBranchAddress("nJets", &nJets);

  TBranch *jetTtagBranch = signalTree->Branch("jetTtagCategory", "vector<float>", &jetTtag);

  for(int i=0; i<signalTree->GetEntries(); i++)
  {
    jetTtag->clear();
    signalTree->GetEntry(i);

    for(int j=0; j<nJets; j++)
    {
      float score;
      if((*jetPt)[j]>=400)
      {
        score = topDiscriminator->eval((*jetTau3)[j], (*jetTau2)[j], (*jetTau1)[j], (*jetMassSub0)[j], (*jetMassSub1)[j], (*jetPt)[j]);
      }
      else
      {
        score = -10;
      }
      
      jetTtag->push_back(score);
    }

    int returnCode = jetTtagBranch->Fill();
  }
  
  std::cout<<"Loop ended"<<std::endl;
  
  signalTree->Print();
  signalTree->Write();
  f->Close();

  delete f;
  delete jetTtag;
  delete topDiscriminator;
}
