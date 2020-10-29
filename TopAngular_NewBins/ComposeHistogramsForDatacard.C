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



void ComposeHistogramsForDatacard(TString year="2016", int mJJCut = 1000)
{
  initFilesMapping();
  //read the data file and get the histogram for our new SR
  TString histoName = "hWt_chi_2btag_expYield";
  TFile *infDataFile = TFile::Open(TString::Format("%s/Histo_Data_%s_reduced_%d.root",year.Data(), year.Data(), mJJCut));
  //our new SR is: SR (old) + mJJ > mJJCut
  TH1F *hData = (TH1F*)infDataFile->Get(histoName);

  //move on to ttbar
  //TFile *infTTfile = TFile::Open(TString::Format("%s/Histo_TT_NominalMC_reduced_%d.root",year.Data(), mJJCut));
  //select the ttbar from extracted signal
  TFile *infTTfile = TFile::Open(TString::Format("%s/FiducialMeasurement_2TeV/EqualBinning/SignalHistograms_chi.root",year.Data()));
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
  TFile *infSubFile = TFile::Open(TString::Format("%s/Histo_SubdominantBkgs_reduced_%d.root",year.Data(), mJJCut));
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
  const int number_of_masses = 5;
  int masses[] = {1000,2000,2500,3000,4000};
  float widths[] = {.01, .1, .3};
  TFile *infZprime;
  TFile *outf_Zprime;
  TH1F *hZprime;

  for (int imass = 0; imass<number_of_masses; imass++)
  {
    for(int iw = 0; iw<3; iw++)
    {
      if(imass==2 && iw ==2) continue;
      float width = masses[imass]* widths[iw];
      infZprime = TFile::Open(TString::Format("../Zprime/%s/HistoMassWindows_ZprimeToTT_M-%d_W-%d_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",
                          year.Data(),masses[imass],(int)width));
      infZprime->cd();
      hZprime = (TH1F*)infZprime->Get(TString::Format("hReco_chi_%d", mJJCut));
      outf_Zprime = new TFile(TString::Format("%s/Datacard/ZprimeFile_%d_%d_massCut%d.root",
                              year.Data(), masses[imass], (int)width, mJJCut), "RECREATE");
      outf_Zprime->cd();
      hZprime->Write("h_chi_Zprime");
    }
  }

}
