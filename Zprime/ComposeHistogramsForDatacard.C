#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1.h"
#include "TCanvas.h"

#include "TemplateConstants.h"

using std::cin;
using std::cout;
using std::endl;



void ComposeHistogramsForDatacard(TString year="2017", int mJJCut = 1000)
{
  //initFilesMapping(false);
  //read the data file and get the histogram for our new SR
  TString histoName = "hWt_chi_2btag_expYield";
  TFile *infDataFile = TFile::Open(TString::Format("../MassFit/%s/Histo_Data_%s_100_reduced_UnequalBinning.root",year.Data(), year.Data()));
  //our new SR is: SR (old) + mJJ > mJJCut
  TH1F *hData = (TH1F*)infDataFile->Get(histoName);

  //move on to ttbar
  //TFile *infTTfile = TFile::Open(TString::Format("%s/Histo_TT_NominalMC_reduced_%d.root",year.Data(), mJJCut));
  //select the ttbar from extracted signal
  TFile *infTTfile = TFile::Open(TString::Format("../MassFit/%s/FiducialMeasurement/UnequalBinning/SignalHistograms_chi.root",year.Data()));
  //our new SR is: SR (old) + mJJ > mJJCut
  TH1F *hTT = (TH1F*)infTTfile->Get("hSignal_chi");
  //hTT->Scale(ttbarSigStrength[year]); //only to be used when looking at ttbar from mc

  //get qcd
  //mc file
  //TFile *infQCDfile = TFile::Open(TString::Format("%s/Histo_QCD_HT300ToInf_reduced_%d.root",year.Data(), mJJCut));
  //TH1F *hQCD = (TH1F*)infQCDfile->Get(histoName);
  //I have written also the extracted qcd signal in the signal extraction file so take it from there
  //our new SR is: SR (old) + mJJ > mJJCut
  TH1F *hQCD = (TH1F*)infTTfile->Get("hQCD_chi");

  //move on to subdominant
  TFile *infSubFile = TFile::Open(TString::Format("../MassFit/%s/Histo_SubdominantBkgs_100_reduced_UnequalBinning.root",year.Data()));
  //our new SR is: SR (old) + mJJ > mJJCut
  TH1F *hSub = (TH1F*)infSubFile->Get(histoName);

  //output files
  TFile *outf_data = new TFile(TString::Format("%s/Datacard/DataFile_%d.root",year.Data(), mJJCut), "RECREATE");
  outf_data->cd();
  hData->Write("h_Data");

  TFile *outf_processes = new TFile(TString::Format("%s/Datacard/ProcessesFile_%d.root", year.Data(),mJJCut), "RECREATE");
  outf_processes->cd();
  hTT->Write("h_chi_ttbar");
  hQCD->Write("h_chi_qcd");
  hSub->Write("h_chi_Subdominant");


  //get the Zprime files:
  const int number_of_masses = 6;
  int masses[] = {2000,2500,3000,3500,4000, 4500};
  float widths[] = {.01, .1, .3};
  TFile *infZprime;
  TFile *outf_Zprime;
  TH1F *hZprime;

  for (int imass = 0; imass<number_of_masses; imass++)
  {
    if (masses[imass] == 3500 && year.Contains("2016")) continue;
    for(int iw = 0; iw<1; iw++)
    {

      if(imass==2 && iw ==2) continue;
      float width = masses[imass]* widths[iw];
      cout<<"mass: "<<masses[imass]<<" width:"<<(int)width<<endl;

      if(year.EqualTo("2016")) infZprime = TFile::Open(TString::Format("%s/HistoMassWindows_ZprimeToTT_M-%d_W-%d_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root", year.Data(), masses[imass],(int)width));
      else if (year.EqualTo("2017")) infZprime = TFile::Open(TString::Format("%s/HistoMassWindows_ZprimeToTT_M%d_W%d_TuneCP2_PSweights_13TeV-madgraph-pythiaMLM-pythia8_20UL.root", year.Data(), masses[imass],(int)width));
      else infZprime = TFile::Open(TString::Format("%s/HistoMassWindows_ZprimeToTT_M%d_W%d_TuneCP2_PSweights_13TeV-madgraph-pythiaMLM-pythia8_20UL.root", year.Data(), masses[imass],(int)width));

      infZprime->cd();
      hZprime = (TH1F*)infZprime->Get(TString::Format("hReco_chi_%d", mJJCut));
      outf_Zprime = new TFile(TString::Format("%s/Datacard/ZprimeFile_%d_%d_massCut%d.root",
                              year.Data(), masses[imass], (int)width, mJJCut), "RECREATE");
      outf_Zprime->cd();
      hZprime->Write("h_chi_Zprime");
    }
  }

}
