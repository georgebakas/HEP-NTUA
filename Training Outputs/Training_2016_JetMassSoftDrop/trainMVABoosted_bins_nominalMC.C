#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/Tools.h"

#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TH1F.h"


using namespace TMVA;
using namespace TMath;

void trainMVABoosted_bins();

void trainMVABoosted_bins()
{
  char name[1000];
  float XSEC[6] = {351400,32260,6830,1207,119.1,25.16};
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

  //TFile *sigSrc = TFile::Open("/eos/cms/store/user/ipapakri/ttbar/MC/Signal/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_legacy2016_Copy.root");
  TFile *sigSrc = TFile::Open("/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/Signal/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_Train.root");
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


  const int NVAR = 10;
  TString VAR[NVAR] = {
    "jetTau3",
    "jetTau2",
    "jetTau1",
    "jetMassSub0",
    "jetMassSub1",
    "ecfB1N2",
    "ecfB1N3",
    "ecfB2N2",
    "ecfB2N3",
    "JetPtOverSumPt"

  };
  char TYPE[NVAR] = {
    'F','F','F','F','F','F','F','F','F','F'
  };

  for(int i=0;i<NVAR;i++) {
    dataloader->AddVariable(VAR[i],TYPE[i]);
  }

  dataloader->AddSpectator("jetPt",'F');
  dataloader->AddSpectator("jetMassSoftDrop", 'F');

  /*sprintf(name,"nTrain_Signal=%d:nTrain_Background=%d:nTest_Signal=%d:nTest_Background=%d",80000,120000,80000,120000);
  dataloader->PrepareTrainingAndTestTree("",name);*/
  TCut selectionCut = "jetPt >= 400";
  sprintf(name,"nTrain_Signal=%d:nTrain_Background=%d:nTest_Signal=%d:nTest_Background=%d:NormMode=EqualNumEvents",0,0,0,0);//:VerboseLevel=Debug
  dataloader->PrepareTrainingAndTestTree(selectionCut, name);

  //factory->BookMethod(dataloader, TMVA::Types::kFisher,"Fisher");
  //factory->BookMethod(dataloader, TMVA::Types::kBDT,"BDT","!H:!V:NTrees=500:BoostType=Grad:Shrinkage=0.1");
  //factory->BookMethod(dataloader, TMVA::Types::kMLP,"MLP","NCycles=1000:HiddenLayers=N+10,N-2:TrainingMethod=BFGS:VarTransform=Norm");
  //factory->BookMethod(dataloader, TMVA::Types::kMLP,"MLP","NCycles=800:HiddenLayers=15,10,5,2:TrainingMethod=BP:VarTransform=Norm");
  TString myVars  = VAR[0];
  for(int i=1; i<NVAR; i++)
  {
    myVars += ":";
    myVars += VAR[i];
  }
  std::cout<<"myVars: "<<myVars<<std::endl;
  TMVA::MethodCategory *mCat = 0;

  TMVA::MethodBase*     fiCat = factory->BookMethod(dataloader, TMVA::Types::kCategory, "FisherCat", "");
  mCat = dynamic_cast<TMVA::MethodCategory*>(fiCat);
  mCat->AddMethod("jetPt>=400 && jetPt<600", myVars, TMVA::Types::kFisher, "Fisher_400_600", "");
  mCat->AddMethod("jetPt>=600 && jetPt<800", myVars, TMVA::Types::kFisher, "Fisher_600_800", "");
  mCat->AddMethod("jetPt>=800 && jetPt<1200", myVars, TMVA::Types::kFisher, "Fisher_800_1200", "");
  mCat->AddMethod("jetPt>=1200", myVars, TMVA::Types::kFisher, "Fisher_1200_Inf", "");


  TMVA::MethodBase* bdtCat = factory->BookMethod(dataloader, TMVA::Types::kCategory, "BDTCat");
  mCat = dynamic_cast<TMVA::MethodCategory*>(bdtCat);
  mCat->AddMethod("jetPt>=400 && jetPt<600", myVars, TMVA::Types::kBDT, "BDT_400_600", "!H:!V:NTrees=500:BoostType=Grad:Shrinkage=0.1");
  mCat->AddMethod("jetPt>=600 && jetPt<800", myVars, TMVA::Types::kBDT, "BDT_600_800", "!H:!V:NTrees=500:BoostType=Grad:Shrinkage=0.1");
  mCat->AddMethod("jetPt>=800 && jetPt<1200", myVars, TMVA::Types::kBDT, "BDT_800_1200", "!H:!V:NTrees=500:BoostType=Grad:Shrinkage=0.1");
  mCat->AddMethod("jetPt>=1200", myVars, TMVA::Types::kBDT, "BDT_1200_Inf", "!H:!V:NTrees=500:BoostType=Grad:Shrinkage=0.1");

	/*
  TMVA::MethodBase* mplCat = factory->BookMethod(dataloader, TMVA::Types::kCategory, "MLPCat");
  mCat = dynamic_cast<TMVA::MethodCategory*>(mplCat);
  mCat->AddMethod("jetPt>=400 && jetPt<600", myVars, TMVA::Types::kMLP, "MLP_400_600", "NCycles=1000:HiddenLayers=N+10,N-2:TrainingMethod=BFGS:VarTransform=Norm");
  mCat->AddMethod("jetPt>=600 && jetPt<800", myVars, TMVA::Types::kMLP, "MLP_600_800", "NCycles=1000:HiddenLayers=N+10,N-2:TrainingMethod=BFGS:VarTransform=Norm");
  mCat->AddMethod("jetPt>=800 && jetPt<1200", myVars, TMVA::Types::kMLP, "MLP_800_1200", "NCycles=1000:HiddenLayers=N+10,N-2:TrainingMethod=BFGS:VarTransform=Norm");
  mCat->AddMethod("jetPt>=1200", myVars, TMVA::Types::kMLP, "MLP_1200_Inf", "NCycles=1000:HiddenLayers=N+10,N-2:TrainingMethod=BFGS:VarTransform=Norm");
	*/
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  outf->Close();
}
