#include "../BASE.h"

#include <TFile.h>
#include <TH1F.h>
#include <TString.h>
#include <TF1.h>
#include <RooFitResult.h>
#include <RooRealVar.h>
#include <TStyle.h>
#include <TLegend.h>

TFile *dataFile, *subFile;
TH1F *signalHistogram, *bkgHistogram, *subHistogram;
TString year, variation;
TString massFitResultsFileName;

/*
 * Applies the corrections calculated from the closure tests to the
 * corresponding histograms.
 * @param variable: The variable we are extracting the signal for
 */
void ApplyClosureTestScaleFactor(TString variable)
{
  //std::cout << AnalysisConstants::closureScaleFactors[year] << std::endl;
  if (!gSystem->AccessPathName(AnalysisConstants::closureScaleFactors[year].Data()))
  {
    TFile *sfFile = TFile::Open(AnalysisConstants::closureScaleFactors[year].Data());
    if (sfFile->GetListOfKeys()->Contains(TString::Format("FitFunction_%s", variable.Data())))
    {
      std::cout << "Correcting with closure test results for variable: " << variable << std::endl;
      TF1 *func = (TF1 *)sfFile->Get(TString::Format("FitFunction_%s", variable.Data()));
      for (int i = 0; i <= bkgHistogram->GetNbinsX(); i++)
      {
        float sf = func->Eval(bkgHistogram->GetBinCenter(i));
        bkgHistogram->SetBinContent(i, bkgHistogram->GetBinContent(i) * sf);
        bkgHistogram->SetBinError(i, bkgHistogram->GetBinError(i) * sf);
      }
    }
    sfFile->Close();
  }
}

/*
 * Returns that integral for a spesific histogram
 * @param f: the file that the histogram is in
 * @param histoName: the name of the histogram for which we want
 *                   the integral.
 */
float GetHistoIntegral(TString fileName, TString pattern, float lumi)
{
  TFile *f = TFile::Open(fileName);
  TH1F *tempHisto = MergeHistograms<TH1F>(f, year, pattern, lumi);
  //std::cout << pattern << " " << tempHisto->GetEntries() << std::endl;
  f->Close();
  return tempHisto->Integral();
}

Double_t GetHistoIntegralAndError(TString fileName, TString pattern, float lumi, Double_t &error)
{
  TFile *f = TFile::Open(fileName);
  TH1F *tempHisto = MergeHistograms<TH1F>(f, year, pattern, lumi);

  Double_t integral = tempHisto->IntegralAndError(1, tempHisto->GetNbinsX(), error);
  Double_t returnIntegral = tempHisto->Integral();
  f->Close();
  return returnIntegral;
}

//deprecated
/*float ABCDMethod(TString var)
{	
	float N_data_0btag_reduced = GetHistoIntegral(dataFile, TString::Format("0btag_%s", var.Data()));
	float N_mcSig_0btag_reduced = GetHistoIntegral(sigFile, TString::Format("0btag_%s", var.Data()));
	float N_qcd_2btag_extended = AnalysisConstants::fitConstants[TString::Format("%s_2btag", year.Data())];
	float N_qcd_0btag_extended = AnalysisConstants::fitConstants[TString::Format("%s_0btag", year.Data())];
	float N_qcd_2btag_reduced = N_qcd_2btag_extended * (N_data_0btag_reduced - N_mcSig_0btag_reduced) / N_qcd_0btag_extended;
	std::cout<<"ABCD: "<<N_qcd_2btag_reduced<<std::endl;
	return N_qcd_2btag_reduced;
}*/

/*
 * Calculates the Ryield transfer factor that transfers
 * the qcd shape from the extended to the reduced region
 */
Double_t CalculateRyield(TString variable, Double_t &error)
{
  Double_t rYield, rYield_error;
  //reduced
  Double_t mcSig_2btag_reduced_error, mcSig_0btag_reduced_error, data_0btag_reduced_error;

  Double_t mcSig_2btag_reduced = GetHistoIntegralAndError(TString::Format("%s/output_%s_mcBkg_reduced.root",
                                                                          GetInputFilesPath(year).Data(),
                                                                          year.Data()),
                                                          TString::Format("2btag_%s",
                                                                          variable.Data()),
                                                          AnalysisConstants::luminositiesSR[year],
                                                          mcSig_2btag_reduced_error);
  Double_t mcSig_0btag_reduced = GetHistoIntegralAndError(TString::Format("%s/output_%s_mcBkg_reduced.root",
                                                                          GetInputFilesPath(year).Data(),
                                                                          year.Data()),
                                                          TString::Format("0btag_%s",
                                                                          variable.Data()),
                                                          AnalysisConstants::luminositiesCR[year],
                                                          mcSig_0btag_reduced_error);
  Double_t data_0btag_reduced = GetHistoIntegralAndError(TString::Format("%s/output_%s_Data_reduced.root",
                                                                         GetInputFilesPath(year).Data(),
                                                                         year.Data()),
                                                         TString::Format("0btag_%s",
                                                                         variable.Data()),
                                                         AnalysisConstants::luminositiesCR[year],
                                                         data_0btag_reduced_error);

  //extended
  Double_t mcSig_2btag_extended_error, mcSig_0btag_extended_error, data_0btag_extended_error;
  Double_t mcSig_2btag_extended = GetHistoIntegralAndError(TString::Format("%s/output_%s_mcBkg_extended.root",
                                                                           GetInputFilesPath(year).Data(),
                                                                           year.Data()),
                                                           TString::Format("2btag_%s",
                                                                           variable.Data()),
                                                           AnalysisConstants::luminositiesSR[year],
                                                           mcSig_2btag_extended_error);
  Double_t mcSig_0btag_extended = GetHistoIntegralAndError(TString::Format("%s/output_%s_mcBkg_extended.root",
                                                                           GetInputFilesPath(year).Data(),
                                                                           year.Data()),
                                                           TString::Format("0btag_%s",
                                                                           variable.Data()),
                                                           AnalysisConstants::luminositiesCR[year],
                                                           mcSig_0btag_extended_error);
  Double_t data_0btag_extended = GetHistoIntegralAndError(TString::Format("%s/output_%s_Data_extended.root",
                                                                          GetInputFilesPath(year).Data(),
                                                                          year.Data()),
                                                          TString::Format("0btag_%s",
                                                                          variable.Data()),
                                                          AnalysisConstants::luminositiesCR[year],
                                                          data_0btag_extended_error);

  rYield = (mcSig_2btag_reduced / mcSig_2btag_extended) / (mcSig_0btag_reduced / mcSig_0btag_extended);
  Double_t rYield_correction = (data_0btag_reduced / data_0btag_extended);

  rYield_error = TMath::Sqrt(TMath::Power(mcSig_0btag_extended / (mcSig_0btag_reduced * mcSig_2btag_extended) * mcSig_2btag_reduced_error, 2) +
                             TMath::Power(mcSig_2btag_reduced / (mcSig_0btag_reduced * mcSig_2btag_extended) * mcSig_0btag_extended_error, 2) +
                             TMath::Power((mcSig_2btag_reduced * mcSig_0btag_extended) / (TMath::Power(mcSig_0btag_reduced, 2) * mcSig_2btag_extended) * mcSig_0btag_reduced_error, 2) +
                             TMath::Power((mcSig_2btag_reduced * mcSig_0btag_extended) / (mcSig_0btag_reduced * TMath::Power(mcSig_2btag_extended, 2)) * mcSig_2btag_extended_error, 2));

  Double_t rYield_correction_error = TMath::Sqrt(TMath::Power(1 / data_0btag_extended * data_0btag_reduced_error, 2) +
                                                 TMath::Power(data_0btag_reduced / TMath::Power(data_0btag_extended, 2) * data_0btag_extended_error, 2));
  error = TMath::Sqrt(TMath::Power(rYield_correction * rYield_error, 2) + TMath::Power(rYield * rYield_correction_error, 2));
  rYield = rYield * rYield_correction;
  std::cout << variable << ": " << rYield << " " << error << std::endl;

  return rYield;
}

/*
 * Returns the number of qcd events calculated in the fit
 */
void GetNqcd(Double_t &val, Double_t &error)
{
  TFile *masFitResultsFile = TFile::Open(TString::Format("%s/fitResults_%s.root",
                                                         massFitResultsFileName.Data(),
                                                         year.Data()));
  masFitResultsFile->GetName();
  RooFitResult *fitResult = (RooFitResult *)masFitResultsFile->Get(TString::Format("fitResults_%s", year.Data()));
  RooRealVar *value = (RooRealVar *)fitResult->floatParsFinal().find("nFitQCD_2b");

  val = value->getValV();
  error = value->getError();
  masFitResultsFile->Close();
}

/*
 * Get the merged histograms that are to be used in the rest of the code
 * 
 * @param variable: the variable we are examining
 */
void GetHistograms(TString variable)
{
  std::cout << dataFile->GetName() << std::endl;
  std::cout << subFile->GetName() << std::endl;
  signalHistogram = MergeHistograms<TH1F>(dataFile, year, TString::Format("2btag_%s", variable.Data()), AnalysisConstants::luminositiesSR[year]);
  bkgHistogram = MergeHistograms<TH1F>(dataFile, year, TString::Format("0btag_%s", variable.Data()), AnalysisConstants::luminositiesCR[year], true);
  subHistogram = MergeHistograms<TH1F>(subFile, year, TString::Format("2btag_%s", variable.Data()), AnalysisConstants::luminositiesSR[year]);
}

/*
 * Does the signal extraction by using the following method:
 * S(x) = D(x) - Ryield * Nqcd * Q(x) - B(x)
 * Where:
 * D(x) is the signal region (2btag) from the data sample
 * Q(x) is the shape of the Contol Region (0btag) taken from the data
 * Nqcd is the number of QCD event calculated from the mass fit
 * Ryield is the transfer factor from the extended to the reduced region
 *        since we do the fit in a different region than the one were we
 *        calculate the cross section
 * B(x) is the signal region taken from the subdominant background samples.
 *
 * @param variable: the variable for we we extract the signal
 */
TH1F *ExtractSignal(TString variable)
{
  ApplyClosureTestScaleFactor(variable);
  Double_t rYield_error;
  Double_t rYield = CalculateRyield(variable, rYield_error);
  Double_t Nqcd, Nqcd_error;
  GetNqcd(Nqcd, Nqcd_error); //get NQCD  from the fit result.
  std::cout << Nqcd << " " << Nqcd_error << std::endl;

  for (int i = 1; i <= bkgHistogram->GetNbinsX(); i++)
  {
    Double_t binError = TMath::Sqrt(TMath::Power(bkgHistogram->GetBinContent(i) * Nqcd * rYield_error, 2) +
                                    TMath::Power(bkgHistogram->GetBinError(i) * Nqcd * rYield, 2) +
                                    TMath::Power(bkgHistogram->GetBinContent(i) * Nqcd_error * rYield, 2));

    Double_t binContent = bkgHistogram->GetBinContent(i) * Nqcd * rYield;
    bkgHistogram->SetBinError(i, binError);
    bkgHistogram->SetBinContent(i, binContent);
  }
  TH1F *oriSignal = (TH1F *)signalHistogram->Clone("oriSignal");
  signalHistogram->Add(bkgHistogram, -1);

  TLegend *leg = new TLegend(0.6, 0.6, 0.9, 0.9);
  leg->AddEntry(signalHistogram, "sig", "l");
  leg->AddEntry(oriSignal, "ori", "l");
  leg->AddEntry(bkgHistogram, "bkg", "l");
  signalHistogram->Draw();
  oriSignal->SetLineColor(kRed);
  bkgHistogram->SetLineColor(kGreen);
  oriSignal->Draw("same");
  bkgHistogram->Draw("same");
  leg->Draw();
  return signalHistogram;
}

void ExtractSignal(TFile *outFile, TString variable)
{
  if (variable.EqualTo(""))
  {
    for (int i = 0; i < AnalysisConstants::variables.size(); i++)
    {
      GetHistograms(AnalysisConstants::variables[i]);
      TH1F *extractedSignal = ExtractSignal(AnalysisConstants::variables[i]);
      outFile->cd();
      extractedSignal->Write(TString::Format("Signal_%s", AnalysisConstants::variables[i].Data()));
    }
  }
  else
  {
    GetHistograms(variable);
    TH1F *extractedSignal = ExtractSignal(variable);

    outFile->cd();
    extractedSignal->Write(TString::Format("Signal_%s", variable.Data()));
  }
}

void SignalExtraction(TString y, bool variations = false, TString variable = "")
{
  year = y;
  gStyle->SetOptStat(0);
  AnalysisConstants::subDir = "";
  AnalysisConstants::initConstants();
  if (variations)
  {
    for (int i = 0; i < AnalysisConstants::variations.size(); i++)
    {
      dataFile = TFile::Open(TString::Format("%s/output_%s_Data_reduced.root",
                                             GetInputFilesPath(year).Data(),
                                             year.Data()));
      subFile = TFile::Open(TString::Format("%s/output_%s_mcSub_reduced.root",
                                            GetInputFilesPath(year).Data(),
                                            year.Data()));
      variation = AnalysisConstants::variations[i];

      TString outputDir = TString::Format("%s/SignalExtraction/results/%s%s/%s%s",
                                          AnalysisConstants::baseDir.Data(),
                                          year.Data(),
                                          (AnalysisConstants::isUL ? "/UL" : ""),
                                          variation.Data(),
                                          AnalysisConstants::currentlyWorkingDirectory[year].Data());
      CheckAndCreateDirectory(outputDir);

      TFile *outFile = TFile::Open(TString::Format("%s/ExtractedSignal_%s.root",
                                                   outputDir.Data(), year.Data()),
                                   "RECREATE");
      massFitResultsFileName = TString::Format("%s/MassFit/results/%s%s/%s%s/extended",
                                               AnalysisConstants::baseDir.Data(),
                                               year.Data(),
                                               (AnalysisConstants::isUL ? "/UL" : ""),
                                               variation.Data(),
                                               AnalysisConstants::currentlyWorkingDirectory[year].Data());
      //std::cout << massFitResultsFileName << std::endl;

      ExtractSignal(outFile, variable);
      outFile->Close();
      dataFile->Close();
      subFile->Close();
    }
  }
  else
  {
    TString variation = "Nominal";
    dataFile = TFile::Open(TString::Format("%s/output_%s_Data_reduced.root",
                                           GetInputFilesPath(year).Data(),
                                           year.Data()));
    subFile = TFile::Open(TString::Format("%s/output_%s_mcSub_reduced.root",
                                          GetInputFilesPath(year).Data(),
                                          year.Data()));
    TString outputDir = TString::Format("%s/SignalExtraction/results/%s%s/%s%s",
                                        AnalysisConstants::baseDir.Data(),
                                        year.Data(),
                                        (AnalysisConstants::isUL ? "/UL" : ""),
                                        variation.Data(),
                                        AnalysisConstants::currentlyWorkingDirectory[year].Data());
    CheckAndCreateDirectory(outputDir);

    TFile *outFile = TFile::Open(TString::Format("%s/ExtractedSignal_%s.root",
                                                 outputDir.Data(),
                                                 year.Data()),
                                 "RECREATE");
    massFitResultsFileName = TString::Format("%s/MassFit/results/%s%s/%s%s/extended",
                                             AnalysisConstants::baseDir.Data(),
                                             year.Data(),
                                             (AnalysisConstants::isUL ? "/UL" : ""),
                                             variation.Data(),
                                             AnalysisConstants::currentlyWorkingDirectory[year]
                                                 .Data());

    ExtractSignal(outFile, variable);
    //outFile->Close();

    //dataFile->Close();
    //subFile->Close();
  }
}