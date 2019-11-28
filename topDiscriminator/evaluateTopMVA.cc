#include "BoostedTopDiscriminator.h"
#include "BoostedTopCategoryDiscriminator.h"
R__LOAD_LIBRARY(BoostedTopDiscriminator.h)
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>

/*
namespace {
  int loadMyLibraryTriggerFunc() {
    gSystem->Load("BoostedTopDiscriminator.h");
    return 0;
  }
  int loadMyLibraryTrigger = loadMyLibraryTriggerFunc();
}*/

void evaluateTopMVA()
{
  TString inputFile = "/eos/cms/store/user/ipapakri/ttbar/MC/Signal/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_legacy2016_copy.root";
  
  TString weights = "/afs/cern.ch/work/i/ipapakri/private/analysis/train/topJetDiscriminator/BoostedMVA_whole/weights/boosted_MVA_MLP.weights.xml";
  TString weightsCategory = "/afs/cern.ch/work/i/ipapakri/private/analysis/train/topJetDiscriminator/BoostedMVA/weights/boosted_MVA_BDTCat.weights.xml";

  TFile *f = new TFile(inputFile, "update");
  BoostedTopDiscriminator *topDiscriminator = new BoostedTopDiscriminator(weights);
  BoostedTopCategoryDiscriminator *topDiscriminatorCategory = new BoostedTopCategoryDiscriminator(weightsCategory);
  TTree* signalTree = (TTree*) f->Get("boosted/events");

  std::vector<float> *jetPt = 0, *jetTau1 = 0, *jetTau2 = 0, *jetTau3 = 0;
  std::vector<float> *jetMassSub0 = 0, *jetMassSub1 = 0;
  std::vector<float> *jetTtag = new std::vector<float>;
  std::vector<float> *jetTtagCategory = new std::vector<float>;
  std::vector<bool> *matchedJet = 0;
  int nJets;
  
  signalTree->SetBranchAddress("jetPt", &jetPt);
  signalTree->SetBranchAddress("jetTau1", &jetTau1);
  signalTree->SetBranchAddress("jetTau2", &jetTau2);
  signalTree->SetBranchAddress("jetTau3", &jetTau3);
  signalTree->SetBranchAddress("jetMassSub0", &jetMassSub0);
  signalTree->SetBranchAddress("jetMassSub1", &jetMassSub1);
  signalTree->SetBranchAddress("nJets", &nJets);
  TBranch *jetTtagBranch = signalTree->Branch("jetTtag", "vector<float>", &jetTtag);
  TBranch *jetTtagCategoryBranch = signalTree->Branch("jetTtagCategory", "vector<float>", &jetTtagCategory);

  for(int i=0; i<signalTree->GetEntries(); i++)
  {
    jetTtag->clear();
    jetTtagCategory->clear();
    signalTree->GetEntry(i);
    //std::cout<<"nJets: "<<nJets<<std::endl;
    for(int j=0; j<nJets; j++)
    {
        float score = topDiscriminator->eval((*jetTau3)[j], (*jetTau2)[j], (*jetTau1)[j], (*jetMassSub0)[j], (*jetMassSub1)[j]);
        //std::cout<<"Score: "<<score<<" matched: "<<(*matchedJet)[j]<<std::endl;
        jetTtag->push_back(score);
        //if(score > 0.5) std::cout<<"Good score"<<std::endl;
        float scoreCategory;
        if((*jetPt)[j] >= 400)
        {
          scoreCategory = topDiscriminatorCategory->eval((*jetTau3)[j], (*jetTau2)[j], (*jetTau1)[j], (*jetMassSub0)[j], (*jetMassSub1)[j], (*jetPt)[j]);
        }
        else
        {
          scoreCategory = -10;
        }

        jetTtagCategory->push_back(scoreCategory);
    }
    //std::cout<<"Size: "<<jetTtag->size()<<std::endl;
    int returnCode = jetTtagBranch->Fill();
    int returnCodeCategory = jetTtagCategoryBranch->Fill();
    //std::cout<<returnCode<<std::endl;
  }
  
  std::cout<<"Loop ended"<<std::endl;
  signalTree->Print();
  signalTree->Write();
  f->Close();

  delete f;
  delete jetTtag;
  delete topDiscriminator;
}
