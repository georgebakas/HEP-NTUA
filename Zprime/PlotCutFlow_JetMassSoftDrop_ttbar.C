#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1.h"
#include "TCanvas.h"
#include "../CMS_plots/tdrstyle.C"
#include "../CMS_plots/CMS_lumi.C"


using std::cin;
using std::cout;
using std::endl;
#include "TemplateConstants.h"
//this code will be used for stacks
void PlotCutFlow_JetMassSoftDrop_ttbar(TString year)
{
  TString histo_names[3] = {"TTToHadronic", "TTToSemiLeptonic", "TTTo2L2Nu"};
  bool isTTbar = true;
  initFilesMapping(isTTbar);
  cout<<year<<endl;
  gStyle->SetOptStat(0);
  setTDRStyle();
  float LUMI = luminosity["luminosity"+year];
  const int NVAR = 2;
  TString vars[NVAR]   = {"jetMass_0", "jetMass_1"};
  int cols[] = {kRed,kBlue, kMagenta,kGreen+3, kTeal, kBlack};
  TString cuts[] = {"Baseline", "jetPt", "mJJ", "btag", "topTagger"};

  lumi_13TeV = TString::Format("%0.1f fb^{-1}", luminosity["luminosity"+year]/1000);
  //lumi_sqrtS = "13 TeV";
  int iPeriod = 4;
  int iPos = 0;
  writeExtraText=true;
  //baseline --> 0,
  //cut_jetPt -->1
  //cut_mJJ -->2 (1500)
  //btagCut-->3
  //tTaggerCut-->4
  const int NCUTS = 5;
  TH1F *hCut[3][NVAR][NCUTS];
  TCanvas *can[NVAR];
  TLegend *leg[NVAR];
  for(int ifile; ifile<3; ifile++)
  {
    TFile *inf = TFile::Open(files[year][histo_names[ifile]]);
    for(int ivar = 0; ivar<NVAR; ivar++)
    {
      for(int icut=0; icut<NCUTS; icut++)
      {
        hCut[ifile][ivar][icut] =(TH1F*)inf->Get(TString::Format("hReco_%s_%d",vars[ivar].Data(), icut));
        if(ifile >0) hCut[0][ivar][icut]->Add(hCut[ifile][ivar][icut]);
      }//end of icut loop
    }//end of ivar loop
  }//end of ifle loop

  for(int ivar = 0; ivar<NVAR; ivar++)
  {
    can[ivar] = new TCanvas("can_"+vars[ivar], "can_"+vars[ivar], 800,600);
    leg[ivar] = new TLegend(0.7, 0.7,0.9, 0.9);
    for(int icut=0; icut<NCUTS; icut++)
    {
      hCut[0][ivar][icut]->SetLineColor(cols[icut]);
      hCut[0][ivar][icut]->SetMarkerStyle(20+icut);
      hCut[0][ivar][icut]->SetMarkerColor(cols[icut]);
      hCut[0][ivar][icut]->GetXaxis()->SetTitle(vars[ivar] + " (GeV)");
      hCut[0][ivar][icut]->GetYaxis()->SetTitle("Expected Yield");
      if(icut==0) hCut[0][ivar][icut]->Draw();
      else hCut[0][ivar][icut]->Draw("same");
      leg[ivar]->AddEntry(hCut[0][ivar][icut],cuts[icut],"ELP");
    }//end of icut loop
    leg[ivar]->Draw();
    CMS_lumi(can[ivar], iPeriod, iPos);
    can[ivar]->Print(TString::Format("%s/CutFlowPlots_JetMassSoftDrop/CutFlow_%s_ttbar.pdf", year.Data(), vars[ivar].Data()), "pdf");

  } //end of ivar loop


}
