#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1.h"
#include "TCanvas.h"

using std::cin;
using std::cout;
using std::endl;
#include "TemplateConstants.h"

void PlotVariablesMC(TString file_names_Zprime, TString level, TString year, TString histoNames, TString mJJCut)
{
  initFilesMapping();
  float LUMI = luminosity["luminosity"+year];
  const int NVAR = 4;
  TString vars[NVAR] = {"chi", "cosTheta_0", "cosTheta_1", "mJJ"};
  int cols[] = {kRed, kRed+1, kRed+3,kBlue, kBlue+1, kBlue+3,kMagenta, kMagenta+1, kMagenta+3, kGreen, kGreen+1, kGreen+3, kYellow+2,kYellow+3, kYellow+4,
                kTeal, kTeal-3, kTeal-6};
  //open the mc ttbar file:
  TFile *inf_fid[NVAR], *inf_par;
  TH1F *hSig[NVAR];
  if(level.EqualTo("Reco"))
  {
    /*
    inf_fid[0] = TFile::Open(TString::Format("../TopAngular_NewBins/%s/FiducialMeasurement/EqualBinning/SignalHistograms_chi.root", year.Data()));
    inf_fid[1] = TFile::Open(TString::Format("../TopAngular_NewBins/%s/FiducialMeasurement/EqualBinning/SignalHistograms_cosTheta_0.root", year.Data()));
    inf_fid[2] = TFile::Open(TString::Format("../TopAngular_NewBins/%s/FiducialMeasurement/EqualBinning/SignalHistograms_cosTheta_1.root", year.Data()));
    */
    inf_fid[0] = TFile::Open(TString::Format("../TopAngular_NewBins/%s/Histo_TT_NominalMC_reduced_%s.root", year.Data(), mJJCut.Data()));
    for(int ivar = 0; ivar<NVAR; ivar++)
    {
      //hSig[ivar] = (TH1F*)inf_fid[ivar]->Get(TString::Format("hSignal_%s", vars[ivar].Data()));
      hSig[ivar] = (TH1F*)inf_fid[0]->Get(TString::Format("hWt_%s_2btag_expYield", vars[ivar].Data()));
      float integral = hSig[ivar]->Integral();
      hSig[ivar]->Scale(1/LUMI, "width");
      hSig[ivar]->Scale(1/integral);
    }

  }
  else if(level.EqualTo("Parton") || level.EqualTo("Particle"))
  {
    inf_par = TFile::Open(TString::Format("../TopAngular_NewBins/%s/%sMeasurements/Data/OutputFile.root", year.Data(), level.Data()));
    for(int ivar = 0; ivar<NVAR; ivar++)
    {
      //hSig[ivar] = (TH1F*)inf_par->Get(TString::Format("hUnfoldFinal_%s", vars[ivar].Data()));
      hSig[ivar] = (TH1F*)inf_par->Get(TString::Format("hTheoryFinal_%s", vars[ivar].Data()));
      float integral = hSig[ivar]->Integral();
      hSig[ivar]->Scale(1/LUMI, "width");
      hSig[ivar]->Scale(1/integral);
    }
  }

  file_names_Zprime.ReplaceAll(']', "");
  file_names_Zprime.ReplaceAll('[', "");
  file_names_Zprime.ReplaceAll(' ', "");

  histoNames.ReplaceAll(']', "");
  histoNames.ReplaceAll('[', "");
  histoNames.ReplaceAll(' ', "");
  TObjArray *tx = file_names_Zprime.Tokenize(",");
  TObjArray *t_names = histoNames.Tokenize(",");
  //tx->Print();
  //t_names->Print();
  TH1F *hZ[NVAR];
  TCanvas *can[NVAR];
  for(int ivar= 0; ivar<NVAR; ivar++)
  {
    //cout<<year<<"/"<<((TObjString *)(tx->At(i)))->String()<<endl;
    can[ivar] = new TCanvas(TString::Format("can_%s", vars[ivar].Data()), TString::Format("can_%s", vars[ivar].Data()), 800,600);
    hSig[ivar]->SetLineColor(kBlack);
    hSig[ivar]->SetMarkerStyle(20);
    hSig[ivar]->SetMarkerColor(kBlack);
    hSig[ivar]->Draw();
    for (Int_t i = 0; i < tx->GetEntries(); i++)
    {
      cout<<year+"/"+((TObjString *)(tx->At(i)))->String()<<endl;
      TFile *infZprime = TFile::Open(year+"/"+((TObjString *)(tx->At(i)))->String());
      hZ[ivar] = (TH1F*)infZprime->Get(TString::Format("h%s_%s_%s", level.Data(), vars[ivar].Data(), mJJCut.Data()));
      hZ[ivar]->SetName(((TObjString *)(t_names->At(i)))->String());
      hZ[ivar]->SetTitle(((TObjString *)(t_names->At(i)))->String());
      hZ[ivar]->SetLineColor(cols[i]);
      float integral = hZ[ivar]->Integral();
      hZ[ivar]->Scale(1/LUMI,"width");
      hZ[ivar]->Scale(1/integral);
      hZ[ivar]->Draw("same");
    }


  }




}
