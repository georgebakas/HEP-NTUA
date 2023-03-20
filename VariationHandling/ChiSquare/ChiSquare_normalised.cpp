#include "../BASE.h"

#include "TMath.h"
#include "TString.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "ChiSquareConstants.h"
#include "FinalResultsConstants.h"
#include "../../CMS_plots/CMS_lumi.C"
#include "../../CMS_plots/tdrstyle.C"
#include "TemplateConstants.h"

TMatrixD *CalculateChiSquare(TMatrixD covariance, TH1F *nominalHistogram, TH1F *theoryHistogram,
                             TString title, TLegend *leg)
{
  std::cout << "Skata" << std::endl;
  covariance.Print();
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

  return chiSquare;
}

void ChiSquare_normalised(TString subDir = "", TString partonParticle = "Particle")
{
  TString theoryYear = "2018";

  AnalysisConstants::subDir = subDir;
  AnalysisConstants::initConstants();
  bool normalized = true;
  TString baseDir = "/Users/georgebakas/Documents/HEP-NTUA_ul/VariationHandling";
  TString outputDir = TString::Format("%s/ChiSquare/results/",
                                      baseDir.Data());
  
  CheckAndCreateDirectory(outputDir);
  TFile *outputFile = TFile::Open(TString::Format("%s/%s_chiSquareResults_normalised.root",
                                                  outputDir.Data(),
                                                  partonParticle.Data()),
                                                  "recreate");

  TFile *nominalFile = TFile::Open(TString::Format("%s/UnfoldedCombined/Nominal/OutputFile%s.root",
                                                    baseDir.Data(),
                                                    partonParticle.Data()));

  TFile *nominalFileTopSF = TFile::Open(TString::Format("%s/UnfoldedCombined/NominalTopSF/OutputFile%s.root",
                                                    baseDir.Data(),
                                                    partonParticle.Data()));

  TFile *finalResultFile = TFile::Open(TString::Format("%s/FinalResults/results/%s_outputFileNorm.root",
                                                       baseDir.Data(),
                                                       partonParticle.Data())); 
                                                  
  TFile *theoryFileAmcAtNlo = TFile::Open(TString::Format(
                            "../../VariationHandling_Theory_amc@NLO/%s/Nominal/Histograms_TTJets.root", theoryYear.Data()));


  TFile *theoryFileHerwig = TFile::Open(TString::Format(
                            "../../VariationHandling_Theory_Herwig/%s/Nominal/Histograms_TT.root", theoryYear.Data()));

  for (int var = 0; var < AnalysisConstants::unfoldingVariables.size(); var++)
  {
    TString variable = AnalysisConstants::unfoldingVariables[var];

    std::cout << variable << std::endl;

    TH1F *nominalHistogram = (TH1F *)nominalFile->Get(TString::Format("hUnfold%s_%s",
                                                                    (normalized ? "Norm" : "Final"), 
                                                                    variable.Data()));

    TH1F *nominalHistogramTopSF = (TH1F *)nominalFileTopSF->Get(TString::Format("hUnfold%s_%s",
                                                                        (normalized ? "Norm" : "Final"), 
                                                                        variable.Data()));

    TH1F *theoryHistogram = (TH1F *)nominalFile->Get(TString::Format("hTheory%s_%s",
                                                                    (normalized ? "Norm" : ""),
                                                                    variable.Data()));

    TH1F *theoryAmcAtNloHistogram;
    if (partonParticle.EqualTo("Parton")) theoryAmcAtNloHistogram = (TH1F *)theoryFileAmcAtNlo->Get(TString::Format("hParton_%s", 
                                                            AnalysisConstants::partonVariables[var].Data()));
    else theoryAmcAtNloHistogram = (TH1F *)theoryFileAmcAtNlo->Get(TString::Format("hParticle_%s", 
                                    AnalysisConstants::particleVariables[var].Data()));
    
    theoryAmcAtNloHistogram->Scale(1. / AnalysisConstants::luminositiesSR["2018"], "width");
    float_t theoryAmcAtNloYield = theoryAmcAtNloHistogram->Integral();
    if (normalized)
    {
        theoryAmcAtNloHistogram->Scale(1 / theoryAmcAtNloYield);
    }

    //herwig
    TH1F *theoryHerwigHistogram;
    if (partonParticle.EqualTo("Parton")) theoryHerwigHistogram = (TH1F *)theoryFileHerwig->Get(TString::Format("hParton_%s", 
                                                            AnalysisConstants::partonVariables[var].Data()));
    else theoryHerwigHistogram = (TH1F *)theoryFileHerwig->Get(TString::Format("hParticle_%s", 
                                    AnalysisConstants::particleVariables[var].Data()));
    
    theoryHerwigHistogram->Scale(1. / AnalysisConstants::luminositiesSR["2018"], "width");
    float_t theoryHerwigYield = theoryHerwigHistogram->Integral();
    if (normalized)
    {
        theoryHerwigHistogram->Scale(1 / theoryHerwigYield);
    }

    TMatrixD *covariance = new TMatrixD(nominalHistogram->GetNbinsX(), nominalHistogram->GetNbinsX());
    TMatrixD *covarianceTheory = new TMatrixD(nominalHistogram->GetNbinsX(), nominalHistogram->GetNbinsX());
    TMatrixD *covarianceTheoryAmcAtNlo = new TMatrixD(nominalHistogram->GetNbinsX(), nominalHistogram->GetNbinsX());
    TMatrixD *covarianceTheoryHerwig = new TMatrixD(nominalHistogram->GetNbinsX(), nominalHistogram->GetNbinsX());
    for (int bin = 1; bin <= nominalHistogram->GetNbinsX(); bin++)
    {
      covariance->operator()(bin - 1, bin - 1) = TMath::Power(nominalHistogram->GetBinError(bin), 2);
      covarianceTheory->operator()(bin - 1, bin - 1) = TMath::Power(theoryHistogram->GetBinError(bin), 2);
      covarianceTheoryAmcAtNlo->operator()(bin - 1, bin - 1) = TMath::Power(theoryAmcAtNloHistogram->GetBinError(bin), 2);
      covarianceTheoryHerwig->operator()(bin - 1, bin - 1) = TMath::Power(theoryHerwigHistogram->GetBinError(bin), 2);
    }

    for (unsigned int v = 0; v < AnalysisConstants::variations.size(); v++)
    {

      TString variation = AnalysisConstants::variations[v];
      //cout<<var<<" "<<variation<<endl;
      TString tempVariation = "";
      if (variation.Contains("isr") || variation.Contains("fsr"))
          tempVariation = "PSWeights";
      else if (variation.EqualTo("bTag_up") || variation.EqualTo("bTag_down"))
          tempVariation = "bTagVariation";
      else if (variation.EqualTo("up") || variation.EqualTo("down"))
          tempVariation = "topTaggingVariation";
      else if (variation.Contains("pdf"))
          tempVariation = "PDFWeights";
      else if (variation.Contains("scale"))
          tempVariation = "ScaleWeights";
      else tempVariation = "JES";

      cout<<"----"<<endl;
      cout<<variation.Data()<<endl;
      cout<<tempVariation.Data()<<endl;

      if (tempVariation.Contains("topTagging"))
            nominalHistogram = nominalHistogramTopSF;


      if (variation.EqualTo("topTaggingup")) continue;
      if (variation.EqualTo("topTaggingdown")) continue;
      
      if (tempVariation.Contains("bTag") || variation.Contains("JES"))
      {
        if (variation.Contains("bTag")) variation = variation.ReplaceAll("bTag_", "");
        
        TString variationDown = variation;
        if (variation.Contains("Up"))
        {
          variationDown.ReplaceAll("Up", "Down");
        }
        else if (variation.Contains("up"))
        {
          variationDown.ReplaceAll("up", "down");
        }
        else 
        {
          variationDown.ReplaceAll("UP", "DOWN");
        }

        TFile *variationFileUp = TFile::Open(TString::Format("%s/UnfoldedCombined/%s/OutputFile%s_%s.root",
                                                                baseDir.Data(),
                                                                tempVariation.Data(),
                                                                partonParticle.Data(),
                                                                variation.Data()));

        TH1F *variationHistogramUp = (TH1F *)variationFileUp->Get(TString::Format("hUnfold%s_%s",
                                                                            (normalized ? "Norm" : "Final"), 
                                                                            variable.Data()));
        TFile *variationFileDown = TFile::Open(TString::Format("%s/UnfoldedCombined/%s/OutputFile%s_%s.root",
                                                                baseDir.Data(),
                                                                tempVariation.Data(),
                                                                partonParticle.Data(),
                                                                variation.Data()));

        TH1F *variationHistogramDown = (TH1F *)variationFileUp->Get(TString::Format("hUnfold%s_%s",
                                                                            (normalized ? "Norm" : "Final"), 
                                                                            variable.Data()));

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
      else
      {
        if (variation.Contains("pdf_99")) continue;
        if (variation.Contains("pdf_98")) continue;
        if (variation.Contains("pdf_100")) continue;
        if (variation.Contains("pdf_1")) continue;
        if (variation.Contains("pdf_35")) continue;
        if (variation.Contains("pdf_40")) continue;
        if (variation.Contains("pdf_56")) continue;
        if (variation.Contains("pdf_57")) continue;
        if (variation.Contains("pdf_59")) continue;
        if (variation.Contains("pdf_62")) continue;
        if (variation.Contains("pdf_83")) continue;
        if (variation.Contains("pdf_86")) continue;
        if (variation.Contains("pdf_85")) continue;
        if (variation.Contains("pdf_76")) continue;
        if (variation.Contains("pdf_42")) continue;
        if (variation.Contains("pdf_58")) continue;
        if (variation.Contains("pdf_69")) continue;
        if (variation.Contains("pdf_93")) continue; 

        /*cout<<TString::Format("%s/UnfoldedCombined/%s/OutputFile%s_%s.root",
                                                                baseDir.Data(),
                                                                tempVariation.Data(),
                                                                partonParticle.Data(),
                                                                variation.Data())<<endl; */
        TFile *variationFile = TFile::Open(TString::Format("%s/UnfoldedCombined/%s/OutputFile%s_%s.root",
                                                                baseDir.Data(),
                                                                tempVariation.Data(),
                                                                partonParticle.Data(),
                                                                variation.Data()));
        TH1F *variationHistogram = (TH1F *)variationFile->Get(TString::Format("hUnfold%s_%s",
                                                                            (normalized ? "Norm" : "Final"), 
                                                                            variable.Data()));

        if (variable.Contains("jetY"))
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


    covarianceTheory->operator+=(*covariance);
    covarianceTheoryAmcAtNlo->operator+=(*covariance);
    covarianceTheoryHerwig->operator+=(*covariance);

    finalResult->Draw("E2");
    nominalHistogram->Draw("SAME EP");
    nominalHistogram->SetMarkerStyle(20);
    nominalHistogram->SetLineColor(kBlack);
    finalResult->GetXaxis()->SetTitle(AnalysisConstants::partonAxisTitles[var]);
    finalResult->GetXaxis()->SetLabelSize(0.035);
    finalResult->GetXaxis()->SetLabelOffset(0.015);
    leg->AddEntry(nominalHistogram, "Data", "E P");
    leg->AddEntry(finalResult, "Total unc.", "f");
    TMatrixD *chiSquare = CalculateChiSquare(*covariance, nominalHistogram, theoryHistogram, "POW + PY8", leg);
    TMatrixD *chiSquareAmc = CalculateChiSquare(*covariance, nominalHistogram, theoryAmcAtNloHistogram, "AMC + PY8", leg);
    TMatrixD *chiSquareHerwig = CalculateChiSquare(*covarianceTheoryHerwig, nominalHistogram, theoryHerwigHistogram, "POW + Hewig", leg);
    theoryHistogram->SetLineColor(kRed);
    theoryHistogram->Draw("SAME");
    theoryAmcAtNloHistogram->SetLineColor(kGreen);
    theoryAmcAtNloHistogram->Draw("SAME");
    theoryHerwigHistogram->SetLineColor(kMagenta);
    theoryHerwigHistogram->Draw("SAME");
    leg->Draw();

    cmsTextSize = 0.5;
    lumiTextSize = 0.5;
    extraTextFactor = 0.13;
    writeExtraText = true;
    CMS_lumi(c1, "combined", 0);

    c1->SaveAs(TString::Format("%s/%s_%s_norm.pdf",
                               outputDir.Data(),
                               partonParticle.Data(),
                               variable.Data()),
               "png");
    c1->SaveAs(TString::Format("%s/%s_%s_norm.pdf",
                               outputDir.Data(),
                               partonParticle.Data(),
                               variable.Data()),
               "pdf");

    outputFile->cd();
    chiSquare->Write(TString::Format("chiSquare_%s_normalized",
                                     variable.Data()));
    chiSquareAmc->Write(TString::Format("chiSquare_%s_amc_normalized",
                                        variable.Data()));
    chiSquareHerwig->Write(TString::Format("chiSquare_%s_herwig_normalized",
                                        variable.Data()));
  }
  outputFile->Close();
}