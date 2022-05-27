#include "../BASE.h"

#include "FinalResultsConstants.h"

#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMatrixD.h"

void DrawWithRatio(TCanvas *can, std::vector<TH1F *> histograms, int index, double chiSquare, double chiSquareAmc, bool normalized)
{
  std::vector<Color_t> colors = {kBlack, kBlack, kRed, kRed, kBlue, kBlue};
  std::vector<Color_t> fillColors = {kGray, kGray, kRed, kRed, kBlue, kBlue};
  std::vector<Int_t> markerStyle = {20, 1, 1, 1, 1, 1};
  std::vector<int> fillStyles = {1001, 1001, 3675, 3675, 3257, 3257};
  std::vector<TString> drawOptions = {"E2", "SAME EP", "SAME E2", "SAME", "SAME E2", "SAME"};
  std::vector<bool> drawRatio = {true, true, true, true, true, true};
  std::vector<bool> putInLegend = {true, true, true, true, true, true};
  std::vector<TString> legendTitles = {"Data", "Total unc.",
                                       TString::Format("amc@NLO+Pythia8 %.2f", chiSquareAmc),
                                       "amc@NLO+Pythia8 unc",
                                       TString::Format("Powheg+Pythia8 %.2f", chiSquare),
                                       "Powheg+Pythia8 unc"};
  std::vector<TString> legendDrawOption = {"EP", "f", "EL", "f", "EL", "f"};

  TPad *upperPad = new TPad("upperPad", "upperPad", 0, 0.3, 1, 1.0);
  upperPad->SetBottomMargin(0.05);
  upperPad->Draw();
  upperPad->cd();
  if (AnalysisConstants::axisInLogScale[index])
  {
    upperPad->SetLogy();
  }
  can->cd();
  TPad *lowerPad = new TPad("lowerPad", "lowerPad", 0, 0.05, 1, 0.3);
  lowerPad->SetBottomMargin(0.22);
  lowerPad->SetTopMargin(0.03);
  lowerPad->Draw();
  lowerPad->SetGridy();
  lowerPad->cd();

  TLegend *leg1 = new TLegend(AnalysisConstants::FinalResultsConstans::legendPositions[index][0],
                              AnalysisConstants::FinalResultsConstans::legendPositions[index][1],
                              AnalysisConstants::FinalResultsConstans::legendPositions[index][2],
                              AnalysisConstants::FinalResultsConstans::legendPositions[index][3]);

  TH1F *dataHist = (TH1F *)histograms[1]->Clone("DataHist");

  for (unsigned int i = 0; i < histograms.size(); i++)
  {
    can->cd();
    upperPad->cd();
    if (AnalysisConstants::axisInLogScale[index])
    {
      upperPad->SetLogy();
    }
    TH1F *hist = histograms[i];

    hist->SetMarkerColor(colors[i]);
    hist->SetLineColor(colors[i]);
    hist->SetFillColor(fillColors[i]);
    hist->SetMarkerStyle(markerStyle[i]);
    hist->SetLineWidth(2);
    hist->SetFillStyle(fillStyles[i]);
    if (normalized)
    {
      hist->GetYaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::partonYAxisValuesNormalized[index][0],
                                     AnalysisConstants::FinalResultsConstans::partonYAxisValuesNormalized[index][1]);
    }
    else
    {
      hist->GetYaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::partonYAxisValues[index][0],
                                     AnalysisConstants::FinalResultsConstans::partonYAxisValues[index][1]);
    }

    TH1F *ratio = (TH1F *)hist->Clone("ratio");
    ratio->Divide(dataHist);
    for (int bin = 1; bin <= hist->GetNbinsX(); bin++)
    {
      ratio->SetBinContent(bin, ratio->GetBinContent(bin) - 1);
    }
    hist->GetXaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::partonXAxisValues[index][0],
                                   AnalysisConstants::FinalResultsConstans::partonXAxisValues[index][1]);
    std::cout << hist->GetTitle() << " " << hist->GetLineColor() << " " << hist->GetFillColor() << std::endl;
    hist->Draw(drawOptions[i]);

    can->cd();
    lowerPad->cd();

    if (drawRatio[i])
    {
      ratio->GetYaxis()->SetRangeUser(-2, 2);
      ratio->GetYaxis()->SetNdivisions(505);
      ratio->GetYaxis()->SetLabelSize(0.1);
      ratio->GetYaxis()->SetTitle("The./data - 1");
      ratio->GetYaxis()->SetTitleSize(0.115);
      ratio->GetYaxis()->SetTitleOffset(0.35);
      ratio->GetYaxis()->CenterTitle(true);
      ratio->GetXaxis()->SetTitleSize(0.115);
      ratio->GetXaxis()->SetTitleOffset(1.);
      ratio->GetXaxis()->SetLabelSize(0.12);
      ratio->GetXaxis()->SetLabelOffset(0.015);
      ratio->GetXaxis()->SetTitle(AnalysisConstants::partonAxisTitles[index]);
      ratio->GetXaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::partonXAxisValues[index][0],
                                      AnalysisConstants::FinalResultsConstans::partonXAxisValues[index][1]);
      ratio->Draw(drawOptions[i]);
    }

    if (putInLegend[i])
    {
      leg1->AddEntry(hist, legendTitles[i], legendDrawOption[i]);
    }
  }

  can->cd();
  upperPad->cd();
  leg1->SetBorderSize(0);
  leg1->Draw();

  /*TList *primitives = upperPad->GetListOfPrimitives();

  for (int i = 0; i < primitives->GetSize(); i++)
  {
    std::cout << primitives->At(i)->GetName() << std::endl;
  }*/

  extraTextFactor = 0.14;
  writeExtraText = true;
  CMS_lumi(upperPad, "combined", 0);
}

void DrawWithChiSquare(TString subDir = "")
{
  AnalysisConstants::subDir = subDir;
  AnalysisConstants::initConstants();

  TFile *resultsFile = TFile::Open(TString::Format("%s/FinalResults/results/%soutputFile.root",
                                                   AnalysisConstants::baseDir.Data(),
                                                   AnalysisConstants::subDir.Data()));
  TFile *chiSquareFiles = TFile::Open(TString::Format("%s/ChiSquare/results/%schiSquareResults.root",
                                                      AnalysisConstants::baseDir.Data(),
                                                      AnalysisConstants::subDir.Data()));
  chiSquareFiles->ls();

  for (int i = 0; i < AnalysisConstants::unfoldingVariables.size(); i++)
  {
    TString variable = AnalysisConstants::unfoldingVariables[i];
    TCanvas *c = new TCanvas(TString::Format("c_%s",
                                             variable.Data()),
                             TString::Format("c_%s",
                                             variable.Data()),
                             600, 600);
    std::vector<TH1F *> histograms;

    TH1F *f = (TH1F *)resultsFile->Get(TString::Format("FinalResult_%s",
                                                       variable.Data()));
    histograms.push_back(f);
    f = (TH1F *)resultsFile->Get(TString::Format("unfoldedHistogram_%s",
                                                 variable.Data()));
    histograms.push_back(f);
    f = (TH1F *)resultsFile->Get(TString::Format("FinalTheoryAmcAtNlo_%s",
                                                 variable.Data()));
    histograms.push_back(f);
    f = (TH1F *)resultsFile->Get(TString::Format("theoryAmcAtNloHistogramValue_%s",
                                                 variable.Data()));
    histograms.push_back(f);
    f = (TH1F *)resultsFile->Get(TString::Format("FinalTheory_%s",
                                                 variable.Data()));
    histograms.push_back(f);
    f = (TH1F *)resultsFile->Get(TString::Format("theoryHistogramValue_%s",
                                                 variable.Data()));
    histograms.push_back(f);

    TMatrixD *chiSquare = (TMatrixD *)chiSquareFiles->Get(TString::Format("chiSquare_%s",
                                                                          variable.Data()));
    TMatrixD *chiSquare_amc = (TMatrixD *)chiSquareFiles->Get(TString::Format("chiSquare_%s_amc",
                                                                              variable.Data()));
    DrawWithRatio(c, histograms, i, chiSquare->operator()(0, 0), chiSquare_amc->operator()(0, 0), false);

    c->SaveAs(TString::Format("%s/FinalResults/results/%sFinalResult_%s_chiSquare.pdf",
                              AnalysisConstants::baseDir.Data(),
                              AnalysisConstants::subDir.Data(),
                              variable.Data()),
              "pdf");
    c->SaveAs(TString::Format("%s/FinalResults/results/%sFinalResult_%s_chiSquare.png",
                              AnalysisConstants::baseDir.Data(),
                              AnalysisConstants::subDir.Data(),
                              variable.Data()),
              "png");
  }
}