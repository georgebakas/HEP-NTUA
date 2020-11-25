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
void PlotCutFlow_JetMassSoftDrop(TString year, TString infStr, TString histo_names)
{
  bool isTTbar = true;
  initFilesMapping(isTTbar);
  cout<<year<<endl;
  cout<<infStr<<endl;
  cout<<histo_names<<endl;
  gStyle->SetOptStat(0);
  setTDRStyle();
  float LUMI = luminosity["luminosity"+year];
  const int NVAR = 2;
  TString vars[NVAR]   = {"jetMass_0", "jetMass_1"};
  int cols[] = {kRed,kBlue, kMagenta,kGreen+3, kTeal, kBlack};
  TString cuts[] = {"Baseline", "jetPt", "mJJ", "btag", "topTagger"};

  TFile *inf = TFile::Open(year+"/"+infStr);

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
  TH1F *hCut[NVAR][NCUTS];
  TCanvas *can[NVAR];
  TLegend *leg[NVAR];
  for(int ivar = 0; ivar<NVAR; ivar++)
  {
    can[ivar] = new TCanvas("can_"+vars[ivar], "can_"+vars[ivar], 800,600);
    leg[ivar] = new TLegend(0.7, 0.7,0.9, 0.9);
    for(int icut=0; icut<NCUTS; icut++)
    {
      hCut[ivar][icut] =(TH1F*)inf->Get(TString::Format("hReco_%s_%d",vars[ivar].Data(), icut));
      hCut[ivar][icut]->Rebin(2);
      hCut[ivar][icut]->SetLineColor(cols[icut]);
      hCut[ivar][icut]->SetMarkerStyle(20+icut);
      hCut[ivar][icut]->SetMarkerColor(cols[icut]);
      hCut[ivar][icut]->GetXaxis()->SetTitle(vars[ivar] + " (GeV)");

      hCut[ivar][icut]->GetYaxis()->SetTitle("Expected Yield");
      if(icut==0)hCut[ivar][icut]->Draw();
      else hCut[ivar][icut]->Draw("same");
      leg[ivar]->AddEntry(hCut[ivar][icut],cuts[icut],"ELP");
    }
    leg[ivar]->Draw();
    CMS_lumi(can[ivar], iPeriod, iPos);
    can[ivar]->Print(TString::Format("%s/CutFlowPlots_JetMassSoftDrop/CutFlow_%s_%s.pdf", year.Data(), vars[ivar].Data(), histo_names.Data()), "pdf");
  }




}
