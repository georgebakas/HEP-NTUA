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
void Significance_ZprimeMass_Width(TString year="2016", int z_mass=2000, float z_width=0.01);


void Significance(TString year="2016")
{
  const int number_of_masses = 10;
  int masses[] = {1200, 1400, 1600, 1800, 2000, 2500, 3000, 3500, 4000, 4500};
  float widths[] = {.01};//, .1, .3};
  TFile *infZprime;
  TFile *outf_Zprime;
  TH1F *hZprime;

  for (int imass = 0; imass<number_of_masses; imass++)
  {
    //if (masses[imass] == 1400 && year.Contains("2016")) continue;
    for(int iw = 0; iw<1; iw++)
    {
      //if(imass==1 && iw ==2) continue;
      Significance_ZprimeMass_Width(year, masses[imass], widths[iw]);
    }
  }
}

void Significance_ZprimeMass_Width(TString year="2016", int z_mass=2000, float z_width=0.01)
{
  initFilesMapping(false);
  int n = 6;
  float significance_values[6];
  float mJJCuts[] = {1000,1200,1400,1600,1800,2000};
  float width = z_mass*z_width;

  for (int icut=0; icut<n; icut++)
  {
    int mJJCut = (int)mJJCuts[icut];
    cout<<mJJCut<<endl;
    //---------------------------------- START OF DATA --------------------------------------------
    //read the data file and get the histogram for our new SR
    TString histoName = "hWt_chi_2btag_expYield";
    TFile *infDataFile = TFile::Open(TString::Format("../TopAngular_NewBins/%s/Histo_Data_%s_reduced_%d.root",year.Data(), year.Data(), mJJCut));
    //our new SR is: SR (old) + mJJ > mJJCut
    TH1F *hData = (TH1F*)infDataFile->Get(histoName);

    //---------------------------------- END OF DATA --------------------------------------------

    //---------------------------------- START OF TTBAR --------------------------------------------
    //move on to ttbar
    TFile *infTTfile = TFile::Open(TString::Format("../TopAngular_NewBins/%s/Histo_TT_NominalMC_reduced_%d.root",year.Data(), mJJCut));
    //select the ttbar from extracted signal
    //TFile *infTTfile = TFile::Open(TString::Format("../MassFit/%s/FiducialMeasurement/UnequalBinning/SignalHistograms_chi.root",year.Data()));
    //our new SR is: SR (old) + mJJ > mJJCut
    TH1F *hTT = (TH1F*)infTTfile->Get(histoName);
    hTT->Scale(ttbarSigStrength[year]); //only to be used when looking at ttbar from mc

    //---------------------------------- END OF TTBAR --------------------------------------------

    //---------------------------------- START OF SUBDOMINANT --------------------------------------------
    //move on to subdominant
    TFile *infSubFile = TFile::Open(TString::Format("../TopAngular_NewBins/%s/Histo_SubdominantBkgs_reduced_%d.root",year.Data(), mJJCut));
    //our new SR is: SR (old) + mJJ > mJJCut
    TH1F *hSub = (TH1F*)infSubFile->Get(histoName);
    //---------------------------------- END OF SUBDOMINANT --------------------------------------------

    //---------------------------------- START OF QCD --------------------------------------------

    //mc file
    TFile *infQCDfile = TFile::Open(TString::Format("../TopAngular_NewBins/%s/Histo_QCD_HT300ToInf_reduced_%d.root",year.Data(), mJJCut));
    TH1F *hQCD = (TH1F*)infQCDfile->Get("hWt_chi_2btag_expYield");
    //I have written also the extracted qcd signal in the signal extraction file so take it from there
    //our new SR is: SR (old) + mJJ > mJJCut
    //TH1F *hQCD = (TH1F*)infTTfile->Get("hQCD_chi");

    //scale qcd with Data
    //we use a k-factor
    TH1F *hQCD_tempFromData = (TH1F*)hData->Clone("hQCD_tempFromData");
    hQCD_tempFromData->Add(hTT, -1);
    hQCD_tempFromData->Add(hSub, -1);

    float qcdScaleFactor = hQCD_tempFromData->Integral()/hQCD->Integral();
    cout<<"qcdScaleFactor: "<<qcdScaleFactor<<endl;
    hQCD->Scale(qcdScaleFactor);


    //---------------------------------- END OF QCD --------------------------------------------


    // get integral from ttbar, qcd, and subdomint -->bkgs
    float total_bkg_integral = hTT->Integral();
    total_bkg_integral += hQCD->Integral();
    total_bkg_integral += hSub->Integral();
    cout<<"mJJ cut: "<<mJJCut<<endl;
    cout<<"total bkg integral: "<<total_bkg_integral<<endl;

    //get the Zprime files:
    TFile *infZprime;
    TH1F *hZprime;

    //cout<<"mass: "<<masses[imass]<<" width:"<<(int)width<<endl;
    if(year.EqualTo("2016_preVFP")) infZprime = TFile::Open(TString::Format("%s/HistoMassWindows_ZprimeToTT_M%d_W%d_TuneCP2_PSweights_13TeV-madgraph-pythiaMLM-pythia8_20UL.root", year.Data(), z_mass, (int)width));
    else if (year.EqualTo("2017")) infZprime = TFile::Open(TString::Format("%s/HistoMassWindows_ZprimeToTT_M%d_W%d_TuneCP2_PSweights_13TeV-madgraph-pythiaMLM-pythia8_20UL.root", year.Data(), z_mass, (int)width));
    else infZprime = TFile::Open(TString::Format("%s/HistoMassWindows_ZprimeToTT_M%d_W%d_TuneCP2_PSweights_13TeV-madgraph-pythiaMLM-pythia8_20UL.root", year.Data(), z_mass, (int)width));

    //infZprime->cd();
    hZprime = (TH1F*)infZprime->Get(TString::Format("hReco_chi_%d", mJJCut));

    float signal = hZprime->Integral();
    float significance = signal / TMath::Sqrt(signal+total_bkg_integral);
    cout<<"Signal integral: "<<signal<<endl;
    cout<<"Mass: "<<z_mass<<", width: "<<width<<" with significance: "<<significance<<endl;

    significance_values[icut] = significance;
  }

  //now draw a tgraph

  TCanvas *c1 = new TCanvas(TString::Format("Significance_%d_%d", z_mass, (int)width),TString::Format("Significance_%d_%d", z_mass, (int)width),800,600);
  TGraph* gr = new TGraph(n,mJJCuts,significance_values);
  gr->SetTitle(TString::Format("Significance_M%d_W%d", z_mass, (int)width));
  gr->GetYaxis()->SetTitleOffset(1.56);
  gr->GetXaxis()->SetTitle("mJJCut (GeV)");
  gr->GetYaxis()->SetTitle("Significance");
  cout<<"Maximum at: "<<gr->GetMaximum()<<endl;
  gr->Draw("A*");
  c1->Print(TString::Format("%s/SignificancePlots/Significance_%d_%d.pdf", year.Data(), z_mass, (int)width), "pdf");

}
