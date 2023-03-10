#include "../BASE.h"

#include "TMath.h"
#include "TString.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "ChiSquareConstants.h"

void CalculateChiSquare(TMatrixD covariance, TH1F *nominalHistogram, TH1F *theoryHistogram,
                        TString title, TLegend *leg)
{
  std::cout << "Skata" << std::endl;
  covariance.Invert();
  TMatrixD *chiSquare = new TMatrixD(1, 1);
  TMatrixD *bins = new TMatrixD(nominalHistogram->GetNbinsX(), 1);
  TMatrixD *binsT = new TMatrixD(nominalHistogram->GetNbinsX(), 1);
  for (int bin = 1; bin <= nominalHistogram->GetNbinsX(); bin++)
  {
    bins->operator()(bin - 1, 0) = TMath::Abs(nominalHistogram->GetBinContent(bin) -
                                              theoryHistogram->GetBinContent(bin));
    binsT->operator()(bin - 1, 0) = bins->operator()(bin - 1, 0);
  }
  binsT->T();

  TMatrixD *intermediate = new TMatrixD(1, nominalHistogram->GetNbinsX());
  intermediate->Mult(*binsT, covariance);
  chiSquare->Mult(*intermediate, *bins);
  chiSquare->Print();

  title += TString::Format(" #chi^{2}: %.2f", chiSquare->operator()(0, 0));
  // title += TString::Format(" p-value: %f", TMath::Prob(chiSquare->operator()(0, 0), nominalHistogram->GetNbinsX() - 1));
  leg->AddEntry(theoryHistogram, title, "l");
}

void ChiSquare_2()
{
  TString theoryYear = "2018";
  AnalysisConstants::initConstants();
  TString baseInputDir = AnalysisConstants::baseDir;
  baseInputDir = TString::Format("%s/Unfolding/results/combined%s",
                                 baseInputDir.Data(),
                                 (AnalysisConstants::isUL ? "/UL" : ""));
  TString outputDir = TString::Format("%s/ChiSquare/results",
                                      AnalysisConstants::baseDir.Data());
  CheckAndCreateDirectory(outputDir);

  TFile *nominalFile = TFile::Open(TString::Format("%s/Nominal/UnfoldingResults.root",
                                                   baseInputDir.Data()));
  TFile *finalResultFile = TFile::Open(TString::Format("%s/FinalResults/results/outputFile.root",
                                                       AnalysisConstants::baseDir.Data()));
  TFile *theoryFile = TFile::Open(TString::Format("%s/EfficiencyAcceptance_Nominal_%s.root",
                                                  GetInputFilesPath(theoryYear).Data(),
                                                  theoryYear.Data()));
  TFile *theoryFileAmcAtNlo = TFile::Open(TString::Format("%s/EfficiencyAcceptance_amcatnlo_Nominal_%s.root",
                                                          GetInputFilesPath(theoryYear).Data(),
                                                          theoryYear.Data()));

  TFile *theoryFileHerwigAtNlo = TFile::Open(TString::Format("%s/EfficiencyAcceptance_amcatnlo_Nominal_%s.root",
                                                          GetInputFilesPath(theoryYear).Data(),
                                                          theoryYear.Data()));

  for (int var = 0; var < AnalysisConstants::unfoldingVariables.size(); var++)
  {
    TString variable = AnalysisConstants::unfoldingVariables[var];
    std::cout << variable << std::endl;

    TH1F *nominalHistogram = (TH1F *)nominalFile->Get(TString::Format("unfoldedHistogram_%s",
                                                                      variable.Data()));
    /*TH1F *theoryHistogram = (TH1F *)theoryFile->Get(TString::Format("FinalTheory_%s",
                                                                    variable.Data()));
    TH1F *theoryAmcAtNloHistogram = (TH1F *)theoryFile->Get(TString::Format("FinalTheoryAmcAtNlo_%s",
                                                                            variable.Data()));*/
    TH1F *theoryHistogram = MergeHistograms<TH1F>(theoryFile,
                                                  theoryYear,
                                                  TString::Format("EfficiencyDenom_%s",
                                                                  AnalysisConstants::partonVariables[var].Data()),
                                                  AnalysisConstants::luminositiesSR[theoryYear]);
    theoryHistogram->Scale(1. / AnalysisConstants::luminositiesSR[theoryYear], "width");
    TH1F *theoryAmcAtNloHistogram = MergeHistograms<TH1F>(theoryFileAmcAtNlo,
                                                          theoryYear,
                                                          TString::Format("EfficiencyDenom_%s",
                                                                          AnalysisConstants::partonVariables[var].Data()),
                                                          AnalysisConstants::luminositiesSR[theoryYear]);
    theoryAmcAtNloHistogram->Scale(1. / AnalysisConstants::luminositiesSR[theoryYear], "width");

    if (variable.Contains("AbsY"))
    {
      nominalHistogram->Scale(1 / 2.);
      theoryHistogram->Scale(1. / 2.);
      theoryAmcAtNloHistogram->Scale(1. / 2.);
    }

    TMatrixD *covariance = new TMatrixD(nominalHistogram->GetNbinsX(), nominalHistogram->GetNbinsX());
    for (int bin = 1; bin <= nominalHistogram->GetNbinsX(); bin++)
    {
      covariance->operator()(bin - 1, bin - 1) = TMath::Power(nominalHistogram->GetBinError(bin), 2);
    }

    for (unsigned int var = 0; var < AnalysisConstants::variations.size(); var++)
    {
      TString variation = AnalysisConstants::variations[var];
      if (variation.Contains("Up"))
      {
        std::cout << variation << std::endl;
        TString variationDown = variation;
        TFile *variationFileUp = TFile::Open(TString::Format("%s/%s/UnfoldingResults.root",
                                                             baseInputDir.Data(),
                                                             variation.Data()));
        TH1F *variationHistogramUp = (TH1F *)variationFileUp->Get(TString::Format("unfoldedHistogram_%s",
                                                                                  variable.Data()));
        TFile *variationFileDown = TFile::Open(TString::Format("%s/%s/UnfoldingResults.root",
                                                               baseInputDir.Data(),
                                                               variationDown.Data()));
        TH1F *variationHistogramDown = (TH1F *)variationFileDown->Get(TString::Format("unfoldedHistogram_%s",
                                                                                      variable.Data()));
        if (variable.Contains("AbsY"))
        {
          variationHistogramUp->Scale(1. / 2.);
          variationHistogramDown->Scale(1. / 2.);
        }
        for (int binI = 1; binI <= variationHistogramUp->GetNbinsX(); binI++)
        {
          float error_I_up = TMath::Abs(nominalHistogram->GetBinContent(binI) -
                                        variationHistogramUp->GetBinContent(binI));
          float error_I_down = TMath::Abs(nominalHistogram->GetBinContent(binI) -
                                          variationHistogramDown->GetBinContent(binI));
          for (int binJ = binI; binJ <= variationHistogramUp->GetNbinsX(); binJ++)
          {
            float error_J_up = TMath::Abs(nominalHistogram->GetBinContent(binJ) -
                                          variationHistogramUp->GetBinContent(binJ));
            float error_J_down = TMath::Abs(nominalHistogram->GetBinContent(binJ) -
                                            variationHistogramDown->GetBinContent(binJ));

            covariance->operator()(binI - 1, binJ - 1) += 1. / 2. * (error_I_up * error_J_up + error_I_down * error_J_down);

            covariance->operator()(binJ - 1, binI - 1) = covariance->operator()(binI - 1, binJ - 1);
          }
        }
        variationFileUp->Close();
        variationFileDown->Close();
      }
      else if (variation.Contains("Down") || variation.Contains("DOWN") || variation.Contains("mtop"))
      {
        continue;
      }
      else
      {
        std::cout << variation << std::endl;
        TFile *variationFile = TFile::Open(TString::Format("%s/%s/UnfoldingResults.root",
                                                           baseInputDir.Data(),
                                                           variation.Data()));
        TH1F *variationHistogram = (TH1F *)variationFile->Get(TString::Format("unfoldedHistogram_%s",
                                                                              variable.Data()));
        if (variable.Contains("AbsY"))
        {
          variationHistogram->Scale(1. / 2.);
        }
        for (int binI = 1; binI <= variationHistogram->GetNbinsX(); binI++)
        {
          float error_I = TMath::Abs(nominalHistogram->GetBinContent(binI) -
                                     variationHistogram->GetBinContent(binI));
          for (int binJ = binI; binJ <= variationHistogram->GetNbinsX(); binJ++)
          {
            float error_J = TMath::Abs(nominalHistogram->GetBinContent(binJ) -
                                       variationHistogram->GetBinContent(binJ));

            covariance->operator()(binI - 1, binJ - 1) += (error_I * error_J);
            covariance->operator()(binJ - 1, binI - 1) = covariance->operator()(binI - 1, binJ - 1);
          }
        }
        variationFile->Close();
      }
    }
    covariance->Print();

    TCanvas *c1 = new TCanvas(TString::Format("c_%s",
                                              variable.Data()),
                              TString::Format("c_%s",
                                              variable.Data()),
                              600, 600);
    TH1F *finalResult = (TH1F *)finalResultFile->Get(TString::Format("FinalResult_%s",
                                                                     variable.Data()));
    TLegend *leg = new TLegend(AnalysisConstants::ChiSquareConstants::legendPositions[var][0],
                               AnalysisConstants::ChiSquareConstants::legendPositions[var][1],
                               AnalysisConstants::ChiSquareConstants::legendPositions[var][2],
                               AnalysisConstants::ChiSquareConstants::legendPositions[var][3]);
    if (!variable.Contains("AbsY"))
    {
      c1->SetLogy();
    }

    finalResult->Draw("E2");
    nominalHistogram->Draw("SAME EP");
    nominalHistogram->SetMarkerStyle(20);
    nominalHistogram->SetLineColor(kBlack);
    finalResult->GetXaxis()->SetTitle(AnalysisConstants::partonAxisTitles[var]);
    finalResult->GetXaxis()->SetLabelSize(0.035);
    finalResult->GetXaxis()->SetLabelOffset(0.015);
    leg->AddEntry(nominalHistogram, "Data", "E P");
    leg->AddEntry(finalResult, "Total unc.", "f");
    CalculateChiSquare(*covariance, nominalHistogram, theoryHistogram, "POW + PY8", leg);
    CalculateChiSquare(*covariance, nominalHistogram, theoryAmcAtNloHistogram, "AMC + PY8", leg);
    theoryHistogram->SetLineColor(kRed);
    theoryHistogram->Draw("SAME");
    theoryAmcAtNloHistogram->SetLineColor(kGreen);
    theoryAmcAtNloHistogram->Draw("SAME");
    leg->Draw();

    cmsTextSize = 0.5;
    lumiTextSize = 0.5;
    extraTextFactor = 0.13;
    writeExtraText = true;
    CMS_lumi(c1, "combined", 0);

    c1->SaveAs(TString::Format("%s/%s.png",
                               outputDir.Data(),
                               variable.Data()),
               "png");
    c1->SaveAs(TString::Format("%s/%s.pdf",
                               outputDir.Data(),
                               variable.Data()),
               "pdf");
    return;
  }
}