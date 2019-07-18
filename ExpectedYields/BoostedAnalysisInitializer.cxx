#define BoostedAnalysisInitializer_cxx
#include "BoostedAnalysisInitializer.h"

BoostedAnalysisInitializer::BoostedAnalysisInitializer(TString year, bool isSignal)
{
  this->year = year;
  this->isSignal = isSignal;
  if(isSignal)
  {
    eosPath = "/eos/cms/store/user/gbakas/ttbar/topTagger/mc-"+year+"/Signal/";
  }
  else
  {
    eosPath = "";//"/eos/cms/store/user/gbakas/ttbar/topTagger/mc-"+year+"/Bkg/";
  }
  
  if(year.EqualTo("2016"))
  {
    LUMI = 35900;
  }
  else if(year.EqualTo("2017"))
  {
    LUMI = 41350;
  }
  else if(year.EqualTo("2018"))
  {
    LUMI = 59740;
  }
  initialize();
}

void BoostedAnalysisInitializer::initialize()
{
  initXsections();
  initHistoNames();
  if(isSignal)
  {
    initSignalFiles();
  }
  else
  {
    initBkgFiles();
  }
}

void BoostedAnalysisInitializer::initXsections()
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

void BoostedAnalysisInitializer::initHistoNames()
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

void BoostedAnalysisInitializer::initSignalFiles()
{
  if(year.EqualTo("2016"))
  {
    listOfFiles.push_back("TT_Mtt-700to1000_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");
    listOfFiles.push_back("TT_Mtt-1000toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8.root");
  }
  else if(year.EqualTo("2017"))
  {
    listOfFiles.push_back("TT_Mtt-700to1000_TuneCP5_13TeV-powheg-pythia8.root");
    listOfFiles.push_back("TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root");
  }
  else if(year.EqualTo("2018"))
  {
    listOfFiles.push_back("TT_Mtt-700to1000_TuneCP5_PSweights_13TeV-powheg-pythia8.root");
    listOfFiles.push_back("TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8.root");
  }
}

void BoostedAnalysisInitializer::initBkgFiles()
{
  if(year.EqualTo("2016"))
  {
    listOfFiles.push_back("/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/Bkg/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
    listOfFiles.push_back("/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/Bkg/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
    listOfFiles.push_back("/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/Bkg/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
    listOfFiles.push_back("/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/Bkg/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
    listOfFiles.push_back("/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/Bkg/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
    listOfFiles.push_back("/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2016/Bkg/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  }
  else if(year.EqualTo("2017"))
  {
    listOfFiles.push_back("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/2017/QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8.root");
    listOfFiles.push_back("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/2017/QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8.root");
    listOfFiles.push_back("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/2017/QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8.root");
    listOfFiles.push_back("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/2017/QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8.root");
    listOfFiles.push_back("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/2017/QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8.root");
    listOfFiles.push_back("/eos/cms/store/user/ipapakri/ttbar/MC/Bkg/2017/QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8.root");
  }
  else if(year.EqualTo("2018"))
  {
    listOfFiles.push_back("/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2018/Bkg/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8.root");
    listOfFiles.push_back("/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2018/Bkg/QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8.root");
    listOfFiles.push_back("/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2018/Bkg/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8.root");
    listOfFiles.push_back("/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2018/Bkg/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8.root");
    listOfFiles.push_back("/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2018/Bkg/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8.root");
    listOfFiles.push_back("/eos/cms/store/user/gbakas/ttbar/topTagger/mc-2018/Bkg/QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8.root");
  }
}

#undef BoostedAnalysisInitializer_cxx