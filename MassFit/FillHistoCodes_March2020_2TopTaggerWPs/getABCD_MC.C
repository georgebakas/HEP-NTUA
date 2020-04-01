#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRatioPlot.h"

using std::cin;
using std::cout;
using std::endl;
/*
	In this code we try and get the ABCD method by defining 4 regions
	A: 2jets top tagged (Medium WP) and 0 jets b-tagged (Loose WP)
	B: 0jets top tagged (Loose WP)  and 0 jets 0-tagged (Loose WP)
	C: 2jets top tagged (Medium WP) and 2 jets b-tagged (Medium WP)
	D: 0jets top tagged (Loose WP) and 2 jets-btagged (Medium WP)

To get expected yield in SR qcd:
For data:
Fill all A,B,C,D regions 
Extract the ttbar and other bkg contribution
The remainder is the expected QCD

Then I will find the QCD from the MC in region C


*/

void initFileNames()
{
  
  if(selection ==0) //data
  {
    eosPath = TString::Format("%s%s/",eosDataPath.Data(), year.Data());  
    listOfFiles.push_back(dataFiles[year.Data()]);
  }
  else if(selection ==1) //signal ttbar mc mtt
  {
    eosPath = TString::Format("%s%s/Signal/",eosPathMC.Data(), year.Data());  
    cout<<eosPath<<endl;
    cout<<mttFiles[year.Data()]["700-1000"]<<endl;
    listOfFiles.push_back(mttFiles[year.Data()]["700-1000"]);
    listOfFiles.push_back(mttFiles[year.Data()]["1000-Inf"]);
  }
  else if(selection ==2) //bkg mc
  {
    eosPath = TString::Format("%s%s/Bkg/",eosPathMC.Data(), year.Data());  
    listOfFiles.push_back(qcdBkgFiles[year.Data()]["300-500"]);
    listOfFiles.push_back(qcdBkgFiles[year.Data()]["500-700"]);
    listOfFiles.push_back(qcdBkgFiles[year.Data()]["700-1000"]);
    listOfFiles.push_back(qcdBkgFiles[year.Data()]["1000-1500"]);
    listOfFiles.push_back(qcdBkgFiles[year.Data()]["1500-2000"]);
    listOfFiles.push_back(qcdBkgFiles[year.Data()]["2000-Inf"]);
  }
  else if(selection ==3) //subdominant bkgs
  {
    eosPath = TString::Format("%s%s/Bkg/",eosPathMC.Data(), year.Data());  
    if(year.EqualTo("2016"))
    {
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["DY"]);
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["WJetsToQQ_HT180"]);
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_tW_top_5f"]);
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_tW_antitop_5f"]);
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_t-channel_top_4f"]);
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_t-channel_antitop_4f"]);
    }
    else
    {
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["DY"]);
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["WJetsToQQ_HT400"]);
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["WJetsToQQ_HT600"]);
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_tW_top_5f"]);
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_tW_antitop_5f"]);
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_t-channel_top_4f"]);
      listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_t-channel_antitop_4f"]);
      if(year.EqualTo("2018"))
      {
        listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_t-channel_top_5f"]);
        listOfFiles.push_back(subdominantBkgFiles[year.Data()]["ST_t-channel_antitop_5f"]);
      }
    }
  }
  else if(selection ==4) //signal ttbar mc nominal
  { 
    cout<<"nominal!!!"<<endl;
    eosPath = TString::Format("%s%s/Signal/",eosPathMC.Data(), year.Data());  
    cout<<eosPath<<endl;
    if(year.EqualTo("2016")) listOfFiles.push_back(ttNominalFiles[year.Data()]["TTNominal"]);
    else
    {
      listOfFiles.push_back(ttNominalFiles[year.Data()]["TTHadronic_0"]);
      listOfFiles.push_back(ttNominalFiles[year.Data()]["TTSemiLeptonic_0"]);
      listOfFiles.push_back(ttNominalFiles[year.Data()]["TTTo2L2Nu_0"]);
    }
  }
}

void initXsections()
{
  if(selection ==1) //signal ttbar mc
  {
    XSEC.push_back(mttXSEC[year.Data()]["700-1000"]);
    XSEC.push_back(mttXSEC[year.Data()]["1000-Inf"]);
  }
  else if(selection ==2) //bkg mc
  {
    XSEC.push_back(qcdBkgXSEC[year.Data()]["300-500"]);
    XSEC.push_back(qcdBkgXSEC[year.Data()]["500-700"]);
    XSEC.push_back(qcdBkgXSEC[year.Data()]["700-1000"]);
    XSEC.push_back(qcdBkgXSEC[year.Data()]["1000-1500"]);
    XSEC.push_back(qcdBkgXSEC[year.Data()]["1500-2000"]);
    XSEC.push_back(qcdBkgXSEC[year.Data()]["2000-Inf"]);
  }
  else if(selection ==3) //subdominant bkgs
  {
    if(year.EqualTo("2016"))
    {
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["DY"]);
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["WJetsToQQ_HT180"]);
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_tW_top_5f"]);
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_tW_antitop_5f"]);
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_t-channel_top_4f"]);
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_t-channel_antitop_4f"]);
    }
    else
    {
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["DY"]);
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["WJetsToQQ_HT400"]);
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["WJetsToQQ_HT600"]);
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_tW_top_5f"]);
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_tW_antitop_5f"]);
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_t-channel_top_4f"]);
      XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_t-channel_antitop_4f"]);
      if(year.EqualTo("2018"))
      {
        XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_t-channel_top_5f"]);
        XSEC.push_back(subdominantBkgXSEC[year.Data()]["ST_t-channel_antitop_5f"]);
      }
    }
  }
  else if(selection ==4)
  {
    if(year.EqualTo("2016")) XSEC.push_back(ttNominalXSEC[year.Data()]["TTNominal"]);
    else
    {
      XSEC.push_back(ttNominalXSEC[year.Data()]["TTHadronic_0"]);
      XSEC.push_back(ttNominalXSEC[year.Data()]["TTSemiLeptonic_0"]);
      XSEC.push_back(ttNominalXSEC[year.Data()]["TTTo2L2Nu_0"]);
    }
  }

}

void initHistoNames()
{
  
  if(selection ==0) histoNames.push_back(TString::Format("Data_%s", year.Data()));
  else if(selection ==1)
  {
    histoNames.push_back("Signal_histo_Mtt_700_1000"); 
    histoNames.push_back("Signal_histo_Mtt_1000_Inf");
  }
  else if (selection ==2)
  {
    histoNames.push_back("QCD_histo_300_500");
    histoNames.push_back("QCD_histo_500_700");
    histoNames.push_back("QCD_histo_700_1000");
    histoNames.push_back("QCD_histo_1000_1500");
    histoNames.push_back("QCD_histo_1500_2000");
    histoNames.push_back("QCD_histo_2000_Inf");
  }
  else if(selection ==3)
  {
   
    if(year.EqualTo("2016"))
    {
       histoNames.push_back("DYJetsToQQ_HT180");
       histoNames.push_back("WJetsToQQ_HT180");
       histoNames.push_back("ST_tW_top_5f_inclusiveDecays");
       histoNames.push_back("ST_tW_antitop_5f_inclusiveDecays");
       histoNames.push_back("ST_t-channel_top_4f_inclusiveDecays");
       histoNames.push_back("ST_t-channel_antitop_4f_inclusiveDecays");
    }
    else
    {
       histoNames.push_back("DYJetsToQQ_HT180");
       histoNames.push_back("WJetsToQQ_HT400");
       histoNames.push_back("WJetsToQQ_HT600");
       histoNames.push_back("ST_tW_top_5f_inclusiveDecays");
       histoNames.push_back("ST_tW_antitop_5f_inclusiveDecays");
       histoNames.push_back("ST_t-channel_top_4f_inclusiveDecays");
       histoNames.push_back("ST_t-channel_antitop_4f_inclusiveDecays");
      if(year.EqualTo("2018"))
      {
        histoNames.push_back("ST_t-channel_top_5f");
        histoNames.push_back("ST_t-channel_antitop_5f");
      }
    }
  }
  else if(selection ==4)
  {
    if(year.EqualTo("2016")) histoNames.push_back("TTNominal");
    else
    {
      histoNames.push_back("TTHadronic_0");
      histoNames.push_back("TTSemiLeptonic_0");
      histoNames.push_back("TTTo2L2Nu_0");
    }
  }

}

void initGlobals()
{
  initFileNames();
  initXsections();
  initHistoNames();
}



void getABCD_MC(TString year = "2016")