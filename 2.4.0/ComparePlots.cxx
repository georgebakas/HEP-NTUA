#include "../BASE.h"
#include "Settings.h"

#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TFile.h>
#include <TStyle.h>

//#include <vector>

void draw()
{
}

void ComparePlots(bool data = true)
{
  gStyle->SetOptStat(0);
  AnalysisConstants::initConstants();
  AnalysisConstants::debug = true;
  std::vector<Color_t> colors = {kRed, kGreen, kBlue, kMagenta};
  for (unsigned int var = 0; var < AnalysisConstants::unfoldingVariables.size(); var++)
  {
    TCanvas *c1 = new TCanvas(TString::Format("canvas_%s", AnalysisConstants::unfoldingVariables[var].Data()),
                              TString::Format("canvas_%s", AnalysisConstants::unfoldingVariables[var].Data()),
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
    TFile *file = TFile::Open(TString::Format("%s/Blue/2.4.0/outFile.root",
                                              AnalysisConstants::baseDir.Data()));
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
    TH1F *comb = (TH1F *)file->Get(TString::Format("combined_%s",
                                                   AnalysisConstants::unfoldingVariables[var].Data()));
    std::cout << comb->GetName() << " " << comb->GetYaxis()->GetLabelFont() << std::endl;
    std::cout << comb->GetYaxis()->GetLabelOffset() << std::endl;
    std::cout << comb->GetYaxis()->GetLabelSize() << std::endl;
    comb->SetTitle("");
    comb->SetMarkerStyle(20);
    comb->SetMarkerColor(kBlack);
    comb->GetYaxis()->SetTitle("#frac{d#sigma}{dx}");

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
      std::cout << "offset: " << comb->GetYaxis()->GetTitleOffset() << std::endl;
      comb->GetYaxis()->SetTitleOffset(0.5);

      comb->GetXaxis()->SetLabelSize(0);
    }
    comb->SetLineColor(kBlack);
    comb->Draw("pehist");

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
      file = TFile::Open(TString::Format("%s/Unfolding/results/%s/Nominal%s%s/UnfoldingResults_%s.root",
                                         AnalysisConstants::baseDir.Data(),
                                         AnalysisConstants::years[y].Data(),
                                         (AnalysisConstants::isUL ? "/UL" : ""),
                                         AnalysisConstants::currentlyWorkingDirectory[AnalysisConstants::years[y]].Data(),
                                         AnalysisConstants::years[y].Data()));
      TH1F *h;
      if (data)
      {
        h = (TH1F *)file->Get(TString::Format("unfoldedHistogram_%s",
                                              AnalysisConstants::unfoldingVariables[var].Data()));

        h->SetLineColor(colors[y]);
        h->Draw("lsame");
        leg->AddEntry(h, AnalysisConstants::years[y], "l");
      }
      else
      {
        c1->cd();
        upperPad->cd();
        h = (TH1F *)file->Get(TString::Format("theory_%s",
                                              AnalysisConstants::unfoldingVariables[var].Data()));
        h->SetLineColor(colors[y]);
        h->Draw("lsame");
        leg->AddEntry(h, AnalysisConstants::years[y], "l");

        TH1F *ratio = (TH1F *)comb->Clone(TString::Format("ratio_%s", AnalysisConstants::years[y].Data()));
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

    c1->SaveAs(TString::Format("results/%s%s.png",
                               AnalysisConstants::unfoldingVariables[var].Data(),
                               (data ? "" : "_MC")));
  }
}