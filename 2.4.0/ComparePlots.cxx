#include "BASE.h"
#include "Settings.h"

#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TFile.h>
#include <TStyle.h>


//#include <vector>

void draw(bool data, bool normalized)
{
  TString baseInputDir = "/afs/cern.ch/work/g/gbakas/public/HEP-NTUA/";
  TString partonParticleStr = "Parton";
  float extraTextFactor = 0;
  std::cout<< "skata "<<std::endl;

  std::vector<Color_t> colors = {kRed, kGreen, kBlue, kMagenta};
  for (unsigned int var = 0; var < AnalysisConstants::unfoldingVariables.size(); var++)
  {
    TCanvas *c1 = new TCanvas(TString::Format("canvas_%s%s",
                                              AnalysisConstants::unfoldingVariables[var].Data(),
                                              (normalized ? "_normalized" : "")),
                              TString::Format("canvas_%s%s",
                                              AnalysisConstants::unfoldingVariables[var].Data(),
                                              (normalized ? "_normalized" : "")),
                              600, 600);
    TPad *upperPad = new TPad("upperPad", "upperPad", 0, 0.4, 1, 1.0);
    if (!data)
    {
      upperPad->SetBottomMargin(0.05);
      upperPad->Draw();
      upperPad->cd();
      upperPad->SetLogy();
    }

    c1->SetLogy();
    TFile *file = TFile::Open("outFile.root");
    TLegend *leg;
    if (AnalysisConstants::unfoldingVariables[var].EqualTo("yJJ") ||
        AnalysisConstants::unfoldingVariables[var].EqualTo("subleadingJetY") ||
        AnalysisConstants::unfoldingVariables[var].EqualTo("leadingJetY"))
    {
      leg = new TLegend(0.35, 0.1, 0.65, 0.4);
    }
    else
    {
      leg = new TLegend(0.6, 0.6, 0.9, 0.9);
    }
    TH1F *comb = (TH1F*)file->Get(TString::Format("combined_%s%s",
                                                   AnalysisConstants::unfoldingVariables[var].Data(),
                                                   (normalized ? "_normalized" : "")));

    comb->SetTitle("");
    comb->SetMarkerStyle(20);
    comb->SetMarkerColor(kBlack);
    if (normalized)
    {
      comb->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dx}");
    }
    else
    {
      comb->GetYaxis()->SetTitle("#frac{d#sigma}{dx}");
    }

    if (data)
    {
      comb->GetYaxis()->SetLabelFont(42);
      comb->GetYaxis()->SetLabelSize(0.035);
      comb->GetYaxis()->SetTitleFont(42);
      comb->GetYaxis()->SetTitleSize(0.035);

      comb->GetXaxis()->SetLabelFont(42);
      comb->GetXaxis()->SetLabelSize(0.035);
      comb->GetXaxis()->SetTitleFont(42);
      comb->GetXaxis()->SetTitleSize(0.035);
      comb->GetXaxis()->SetTitle(AnalysisConstants::axisTitles[var]);
    }
    else
    {
      comb->GetYaxis()->SetLabelFont(42);
      comb->GetYaxis()->SetLabelSize(0.05);
      comb->GetYaxis()->SetTitleFont(42);
      comb->GetYaxis()->SetTitleSize(0.05);
      comb->GetYaxis()->SetTitleOffset(0.5);

      comb->GetXaxis()->SetLabelSize(0);
    }
    comb->SetLineColor(kBlack);
    comb->Draw("pehist");
    std::cout<< "ok here malaka"<<std::endl;
    leg->AddEntry(comb, "comb", "ple");
    TPad *lowerPad = new TPad("lowerPad", "lowerPad", 0, 0.05, 1, 0.4);
    if (!data)
    {
      lowerPad->SetBottomMargin(0.22);
      lowerPad->SetTopMargin(0);
      c1->cd();
      lowerPad->Draw();
    }
    //file->Close();
    for (unsigned int y = 0; y < AnalysisConstants::years.size(); y++)
    {
      file = TFile::Open(TString::Format("%s/Unfolding/%s/%sMeasurements/Data%s/OutputFile.root",
                                          baseInputDir.Data(),
                                          AnalysisConstants::years[y].Data(),
                                          partonParticleStr.Data(),
                                          (normalized ? "_Norm" : "")));
      TH1F *h;
      if (data)
      {
        std::cout<<"data: "<< AnalysisConstants::years[y]<<std::endl;
        h = (TH1F *)file->Get(TString::Format("hUnfold_%s",AnalysisConstants::unfoldingVariables[var].Data()));

        h->SetLineColor(colors[y]);
        h->Draw("lsame");
        leg->AddEntry(h, AnalysisConstants::years[y], "l");
      }
      else
      {
        c1->cd();
        upperPad->cd();
        h = (TH1F *)file->Get(TString::Format("hTheory_%s", AnalysisConstants::unfoldingVariables[var].Data()));

        h->SetLineColor(colors[y]);
        h->Draw("lsame");
        leg->AddEntry(h, AnalysisConstants::years[y], "l");

        TH1F *ratio = (TH1F *)comb->Clone(TString::Format("ratio_%s%s",
                                                          AnalysisConstants::years[y].Data(),
                                                          (normalized ? "_normalized" : "")));
        ratio->SetLineColor(colors[y]);
        ratio->SetMarkerStyle(1);

        ratio->GetXaxis()->SetLabelFont(42);
        ratio->GetXaxis()->SetTitleFont(42);
        ratio->GetXaxis()->SetTitle(AnalysisConstants::axisTitles[var]);
        ratio->GetXaxis()->SetLabelSize(0.08);
        ratio->GetXaxis()->SetLabelOffset(0.015);
        ratio->GetXaxis()->SetNdivisions(410);
        ratio->GetYaxis()->SetNdivisions(505);
        ratio->GetYaxis()->SetLabelSize(0.08);
        ratio->GetYaxis()->SetTitle("#frac{Data}{MC}");
        ratio->GetYaxis()->SetTitleSize(0.08);
        ratio->GetYaxis()->SetTitleOffset(0.43);
        ratio->GetYaxis()->CenterTitle(true);

        ratio->GetXaxis()->SetTitleSize(0.115);
        ratio->GetXaxis()->SetTitleOffset(1.);

        ratio->Divide(h);
        c1->cd();
        lowerPad->cd();
        ratio->Draw("lsame");
      }
    }

    if (!data)
    {
      c1->cd();
      upperPad->cd();
    }
    leg->Draw();
    if (data)
    {
      lumiTextSize = 0.5;
      cmsTextSize = 0.5;
      extraTextFactor = 0.14;
      writeExtraText = true;
      CMS_lumi(c1, 13, 0);
    }
    else
    {
      /*lumiTextSize = 0.5;
      cmsTextSize = 0.5;
      extraTextFactor = 0.14;
      writeExtraText = true;*/
      extraTextFactor = 0.12;
      writeExtraText = true;
      CMS_lumi(upperPad, 13, 0);
    }

    c1->SaveAs(TString::Format("results/%s%s%s.png",
                               AnalysisConstants::unfoldingVariables[var].Data(),
                               (data ? "" : "_MC"),
                               (normalized ? "_normalized" : "")));
  }
}

void ComparePlots(bool data = true, bool normalised = false)
{
  gStyle->SetOptStat(0);
  AnalysisConstants::initConstants();
  std::cout<< "edw"<<std::endl;
  draw(data, normalised);
}