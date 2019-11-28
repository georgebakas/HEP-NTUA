#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h" 
#include "TH1F.h"
using namespace TMVA;
using namespace TMath;

void trainMVABoosted();

void trainMVABoosted()
{
  char name[1000];
  float XSEC[6] = {3.67e+5,2.94e+4,6.524e+03,1.064e+03,121.5,2.542e+01};
  float NORM[6];
  TFile *bkgSrc[6],*bkgSrcNorm[6];

  bkgSrc[0] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1_legacy2016_Train.root");
  bkgSrc[1] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1_legacy2016_Train.root");
  bkgSrc[2] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1_legacy2016_Train.root");
  bkgSrc[3] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1_legacy2016_Train.root");
  bkgSrc[4] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1_legacy2016_Train.root");
  bkgSrc[5] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1_legacy2016_Train.root");

  bkgSrcNorm[0] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1_legacy2016.root");
  bkgSrcNorm[1] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1_legacy2016.root");
  bkgSrcNorm[2] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1_legacy2016.root");
  bkgSrcNorm[3] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1_legacy2016.root");
  bkgSrcNorm[4] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1_legacy2016.root");
  bkgSrcNorm[5] = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1_legacy2016.root");

  TFile *sigSrc = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Signal/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_legacy2016_Train.root");
  TTree *sigTree = (TTree*)sigSrc->Get("events");
  
  TTree *bkgTree[6];
  
  TFile *outf = new TFile("boosted_MVA.root","RECREATE");
  TMVA::Factory* factory = new TMVA::Factory("boosted_MVA",outf,"!V:!Silent:Color:DrawProgressBar:Transformations=I;G:AnalysisType=Classification");

  TMVA::DataLoader *dataloader = new TMVA::DataLoader("BoostedMVA");

  dataloader->AddSignalTree(sigTree);
 
  for(int k=0;k<6;k++) {
    NORM[k] = ((TH1F*)bkgSrcNorm[k]->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
    bkgTree[k] = (TTree*)bkgSrc[k]->Get("events");
    dataloader->AddBackgroundTree(bkgTree[k],XSEC[k]/NORM[k]);
  }
  
  const int NVAR = 5;
  TString VAR[NVAR] = {
    "jetTau3",
    "jetTau2",
    "jetTau1",
    "jetMassSub0",
    "jetMassSub1"
  };
  char TYPE[NVAR] = {
    'F','F','F','F','F'
  };

  for(int i=0;i<NVAR;i++) {
    dataloader->AddVariable(VAR[i],TYPE[i]);
  }

  //dataloader->AddSpectator("jetPt",'F');

  sprintf(name,"nTrain_Signal=%d:nTrain_Background=%d:nTest_Signal=%d:nTest_Background=%d",80000,120000,80000,120000);
  dataloader->PrepareTrainingAndTestTree("",name);

  factory->BookMethod(dataloader, TMVA::Types::kFisher,"Fisher");
  factory->BookMethod(dataloader, TMVA::Types::kBDT,"BDT","!H:!V:NTrees=500:BoostType=Grad:Shrinkage=0.1");
  factory->BookMethod(dataloader, TMVA::Types::kMLP,"MLP","NCycles=1000:HiddenLayers=N+10,N-2:TrainingMethod=BFGS:VarTransform=Norm");
  //factory->BookMethod(dataloader, TMVA::Types::kMLP,"MLP","NCycles=800:HiddenLayers=15,10,5,2:TrainingMethod=BP:VarTransform=Norm");
  
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods(); 
  outf->Close();
}
