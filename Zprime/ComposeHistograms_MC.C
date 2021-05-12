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

//----ONLY TO CHECK THE FILES!!!!-------------------//
/*
  Here I use the QCD and Signal from signal extraction:
  1. Signal is the MCD ttbar in the SR scaled to ttbarSigStrength
  2. QCD is the MC in the SR scaled to Data
*/

void ComposeHistograms_MC(TString year="2017", int mJJCut = 1000)
{
  initFilesMapping(false);
  //---------------------------------- START OF DATA --------------------------------------------
  //read the data file and get the histogram for our new SR
  TString histoName = "hWt_chi_2btag_expYield";
  TFile *infDataFile = TFile::Open(TString::Format("../TopAngular_NewBins/%s/Histo_Data_%s_reduced_%d.root",year.Data(), year.Data(), mJJCut));
  //our new SR is: SR (old) + mJJ > mJJCut
  TH1F *hData = (TH1F*)infDataFile->Get(histoName);

  //---------------------------------- END OF DATA --------------------------------------------

  //---------------------------------- START OF TTBAR --------------------------------------------
  //move on to ttbar
  TFile *infTTfile = TFile::Open(TString::Format("../VariationHandling/%s/Nominal/combined/HistoReduced_%d_TT.root",year.Data(), mJJCut));
  //select the ttbar from extracted signal
  //TFile *infTTfile = TFile::Open(TString::Format("../MassFit/%s/FiducialMeasurement/UnequalBinning/SignalHistograms_chi.root",year.Data()));
  //our new SR is: SR (old) + mJJ > mJJCut
  TH1F *hTT = (TH1F*)infTTfile->Get("hWt_chi_2btag");
  hTT->Scale(ttbarSigStrength[year]); //only to be used when looking at ttbar from mc

  //---------------------------------- END OF TTBAR --------------------------------------------

  //---------------------------------- START OF TTBAR variations --------------------------------------------
  //Scale Weights
  //we only want ISRDef Lo and Hi
  TFile *infTTfile_isr_lo = TFile::Open(TString::Format("../VariationHandling/%s/PSWeights/combined/HistoReduced_%d_TT_PSWeights_isrDefLo.root",year.Data(), mJJCut));
  TH1F *h_TT_isr_def_lo = (TH1F*)infTTfile_isr_lo->Get("hWt_chi_2btag");

  TFile *infTTfile_isr_hi = TFile::Open(TString::Format("../VariationHandling/%s/PSWeights/combined/HistoReduced_%d_TT_PSWeights_isrDefHi.root",year.Data(), mJJCut));
  TH1F *h_TT_isr_def_hi = (TH1F*)infTTfile_isr_hi->Get("hWt_chi_2btag");

  TFile *infTTfile_fsr_lo = TFile::Open(TString::Format("../VariationHandling/%s/PSWeights/combined/HistoReduced_%d_TT_PSWeights_fsrDefLo.root",year.Data(), mJJCut));
  TH1F *h_TT_fsr_def_lo = (TH1F*)infTTfile_fsr_lo->Get("hWt_chi_2btag");

  TFile *infTTfile_fsr_hi = TFile::Open(TString::Format("../VariationHandling/%s/PSWeights/combined/HistoReduced_%d_TT_PSWeights_fsrDefHi.root",year.Data(), mJJCut));
  TH1F *h_TT_fsr_def_hi = (TH1F*)infTTfile_fsr_hi->Get("hWt_chi_2btag");
  //here we get the FSR Def Lo and Hi


  //get the JES (only nominal Smeared Up, Down and Shifted Up, Down)
  TFile *infTTfile_smeared_up = TFile::Open(TString::Format("../VariationHandling/%s/JES/combined/HistoReduced_%d_TT_JES_boostedSmearedUp.root",year.Data(), mJJCut));
  TH1F *h_TT_smeared_up = (TH1F*)infTTfile_smeared_up->Get("hWt_chi_2btag");

  TFile *infTTfile_smeared_down = TFile::Open(TString::Format("../VariationHandling/%s/JES/combined/HistoReduced_%d_TT_JES_boostedSmearedDown.root",year.Data(), mJJCut));
  TH1F *h_TT_smeared_down = (TH1F*)infTTfile_smeared_up->Get("hWt_chi_2btag");

  TFile *infTTfile_shifted_up = TFile::Open(TString::Format("../VariationHandling/%s/JES/combined/HistoReduced_%d_TT_JES_boostedShiftedUp.root",year.Data(), mJJCut));
  TH1F *h_TT_shifted_up = (TH1F*)infTTfile_smeared_up->Get("hWt_chi_2btag");

  TFile *infTTfile_shifted_down = TFile::Open(TString::Format("../VariationHandling/%s/JES/combined/HistoReduced_%d_TT_JES_boostedShiftedUp.root",year.Data(), mJJCut));
  TH1F *h_TT_shifted_down = (TH1F*)infTTfile_smeared_up->Get("hWt_chi_2btag");

  // Scale Weights 2,3,4,5,7,9
  const int N_scale=6;
  int scale_weight_nums[N_scale] = {2,3,4,5,7,9};
  TFile *infTTfile_scale;
  int var_scale[hTT->GetNbinsX()];
  for (int i=0; i<N_scale; i++)
  {
    infTTfile_scale = TFile::Open(TString::Format("../VariationHandling/%s/ScaleWeights/combined/HistoReduced_%d_TT_scale_%d.root",year.Data(), mJJCut, scale_weight_nums[i]));
    TH1F *h_ = (TH1F*)infTTfile_scale->Get("hWt_chi_2btag");
    cout<<"i var integral: "<<h_->Integral()<<endl;

    for (int ibin=1; ibin<=h_->GetNbinsX(); ibin++)
    {
      if (i==0) var_scale[ibin] =0;
      var_scale[ibin] += TMath::Power(hTT->GetBinContent(ibin)-h_->GetBinContent(ibin), 2);
    }
  }
  for (int ibin=1; ibin<=hTT->GetNbinsX(); ibin++)
  {
    var_scale[ibin] = TMath::Sqrt(var_scale[ibin]/(N_scale));
    cout<<ibin<<" "<<var_scale[ibin]<<" "<<TMath::Sqrt(var_scale[ibin]/(N_scale))<<endl;
  }
  //get TT th1 and add or remove RMS per bin
  TH1F *hTT_scale_up, *hTT_scale_low;
  hTT_scale_up = (TH1F*)hTT->Clone("hTT_SPS_up");
  hTT_scale_low = (TH1F*)hTT->Clone("hTT_PS_low");
  for (int ibin=1; ibin<=hTT->GetNbinsX(); ibin++)
  {
    hTT_scale_up->SetBinContent(ibin, hTT_scale_up->GetBinContent(ibin) + var_scale[ibin]);
    hTT_scale_low->SetBinContent(ibin, hTT_scale_low->GetBinContent(ibin) - var_scale[ibin]);
  }


  //PDF Weights 0-100
  int NPDF = 100;
  TFile *infTTfile_pdf;
  int var[hTT->GetNbinsX()];
  for (int i=0; i<=NPDF; i++)
  {
    infTTfile_pdf = TFile::Open(TString::Format("../VariationHandling/%s/PDFWeights/combined/HistoReduced_%d_TT_pdf_%d.root",year.Data(), mJJCut, i));
    TH1F *h_ = (TH1F*)infTTfile_pdf->Get("hWt_chi_2btag");

    for (int ibin=1; ibin<=h_->GetNbinsX(); ibin++)
    {
      if (i==0) var[ibin] =0;
      var[ibin] += TMath::Power(hTT->GetBinContent(ibin)-h_->GetBinContent(ibin), 2);
    }
  }
  for (int ibin=1; ibin<=hTT->GetNbinsX(); ibin++)
    var[ibin] = TMath::Sqrt(var[ibin]/(NPDF+1));
  //get TT th1 and add or remove RMS per bin
  TH1F *hTT_PDF_up, *hTT_PDF_low;
  hTT_PDF_up = (TH1F*)hTT->Clone("hTT_PDF_up");
  hTT_PDF_low = (TH1F*)hTT->Clone("hTT_PDF_low");
  for (int ibin=1; ibin<=hTT->GetNbinsX(); ibin++)
  {
    hTT_PDF_up->SetBinContent(ibin, hTT_PDF_up->GetBinContent(ibin) + var[ibin]);
    hTT_PDF_low->SetBinContent(ibin, hTT_PDF_low->GetBinContent(ibin) - var[ibin]);
  }

  //bTag Variations
  TFile *infTTfile_btagUp = TFile::Open(TString::Format("../VariationHandling/%s/bTagVariation/combined/HistoReduced_%d_TT_bTagVariation_up.root",year.Data(), mJJCut));
  TFile *infTTfile_btagDown = TFile::Open(TString::Format("../VariationHandling/%s/bTagVariation/combined/HistoReduced_%d_TT_bTagVariation_down.root",year.Data(), mJJCut));

  TH1F *h_TT_btag_up = (TH1F*)infTTfile_btagUp->Get("hWt_chi_2btag");
  TH1F *h_TT_btag_down = (TH1F*)infTTfile_btagDown->Get("hWt_chi_2btag");


  //---------------------------------- END OF TTBAR variations--------------------------------------------



  //---------------------------------- START OF SUBDOMINANT --------------------------------------------
  //move on to subdominant
  TFile *infSubFile = TFile::Open(TString::Format("../TopAngular_NewBins/%s/Histo_SubdominantBkgs_reduced_%d.root",year.Data(), mJJCut));
  //our new SR is: SR (old) + mJJ > mJJCut
  TH1F *hSub = (TH1F*)infSubFile->Get(histoName);
  //---------------------------------- END OF SUBDOMINANT --------------------------------------------

  //---------------------------------- START OF QCD --------------------------------------------

  //mc file
  TFile *infQCDfile = TFile::Open(TString::Format("../TopAngular_NewBins/%s/Histo_QCD_HT300toInf_reduced_%d.root",year.Data(), mJJCut));
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

  //---------------------------------- Start handling of histograms --------------------------------------------
  //output files
  TFile *outf_data = new TFile(TString::Format("%s/Datacard_MC/DataFile_%d.root",year.Data(), mJJCut), "RECREATE");
  outf_data->cd();
  hData->Write("h_Data");

  TFile *outf_processes = new TFile(TString::Format("%s/Datacard_MC/ProcessesFile_%d.root", year.Data(),mJJCut), "RECREATE");
  outf_processes->cd();
  hTT->Write("h_chi_ttbar");
  hQCD->Write("h_chi_qcd");
  hSub->Write("h_chi_Subdominant");

  //scale variations
  hTT_PDF_up->Scale(ttbarSigStrength[year]);
  hTT_PDF_low->Scale(ttbarSigStrength[year]);
  //pdf rms
  hTT_PDF_up->Write("h_chi_ttbar_pdfUp");
  hTT_PDF_low->Write("h_chi_ttbar_pdfDown");

  //scale scale-var
  hTT_scale_up->Scale(ttbarSigStrength[year]);
  hTT_scale_low->Scale(ttbarSigStrength[year]);
  //scale var rms
  hTT_scale_up->Write("h_chi_ttbar_scaleUp");
  hTT_scale_low->Write("h_chi_ttbar_scaleDown");

  //scale ps Weights
  h_TT_isr_def_lo->Scale(ttbarSigStrength[year]);
  h_TT_isr_def_hi->Scale(ttbarSigStrength[year]);
  h_TT_fsr_def_lo->Scale(ttbarSigStrength[year]);
  h_TT_fsr_def_hi->Scale(ttbarSigStrength[year]);
  //ps weights
  h_TT_isr_def_lo->Write("h_chi_ttbar_isrUp");
  h_TT_isr_def_hi->Write("h_chi_ttbar_isrDown");
  h_TT_fsr_def_lo->Write("h_chi_ttbar_fsrUp");
  h_TT_fsr_def_hi->Write("h_chi_ttbar_fsrDown");

  //jes weights
  h_TT_shifted_down->Write("h_chi_ttbar_shiftedDown");
  h_TT_shifted_up->Write("h_chi_ttbar_shiftedUp");
  h_TT_smeared_down->Write("h_chi_ttbar_smearedDown");
  h_TT_smeared_up->Write("h_chi_ttbar_smearedUp");

  //scale btag variations
  h_TT_btag_up->Scale(ttbarSigStrength[year]);
  h_TT_btag_down->Scale(ttbarSigStrength[year]);
  //btag variations
  h_TT_btag_up->Write("h_chi_ttbar_btagUp");
  h_TT_btag_down->Write("h_chi_ttbar_btagDown");


  //get the Zprime files:
  const int number_of_masses = 10;
  int masses[] = {1200,1400,1600,1800,2000,2500,3000,3500,4000, 4500};
  float widths[] = {.01, .1, .3};
  TFile *infZprime;
  TFile *outf_Zprime;
  TH1F *hZprime;

  for (int imass = 0; imass<number_of_masses; imass++)
  {
    if (masses[imass] == 1400 && year.Contains("2016")) continue;
    for(int iw = 0; iw<1; iw++)
    {

      if(imass==2 && iw ==2) continue;
      float width = masses[imass]* widths[iw];
      cout<<"mass: "<<masses[imass]<<" width:"<<(int)width<<endl;

      if(year.EqualTo("2016_preVFP")) infZprime = TFile::Open(TString::Format("%s/HistoMassWindows_ZprimeToTT_M%d_W%d_TuneCP2_PSweights_13TeV-madgraph-pythiaMLM-pythia8_20UL.root", year.Data(), masses[imass],(int)width));
      else if (year.EqualTo("2017")) infZprime = TFile::Open(TString::Format("%s/HistoMassWindows_ZprimeToTT_M%d_W%d_TuneCP2_PSweights_13TeV-madgraph-pythiaMLM-pythia8_20UL.root", year.Data(), masses[imass],(int)width));
      else infZprime = TFile::Open(TString::Format("%s/HistoMassWindows_ZprimeToTT_M%d_W%d_TuneCP2_PSweights_13TeV-madgraph-pythiaMLM-pythia8_20UL.root", year.Data(), masses[imass],(int)width));

      infZprime->cd();
      hZprime = (TH1F*)infZprime->Get(TString::Format("hReco_chi_%d", mJJCut));
      outf_Zprime = new TFile(TString::Format("%s/Datacard_MC/ZprimeFile_%d_%d_massCut%d.root",
                              year.Data(), masses[imass], (int)width, mJJCut), "RECREATE");
      outf_Zprime->cd();
      hZprime->Write("h_chi_Zprime");
    }
  }

}
