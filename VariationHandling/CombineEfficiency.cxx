#include "BASE.h"

#include <iostream>

#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TEfficiency.h"
#include "TGraph.h"

std::vector<TString> listFiles(const char *dirname="", const char *var="", const char *ext=".root")
{
  std::vector<TString> list_of_files;
   TSystemDirectory dir(dirname, dirname);
   TList *files = dir.GetListOfFiles();
   if (files) {
     TSystemFile *file;
     TString fname;
     TIter next(files);
     while ((file=(TSystemFile*)next())) {
       fname = file->GetName();
       if (!file->IsDirectory() && fname.EndsWith(ext) && fname.Contains(var)) {
         //cout << fname.Data() << endl;
         list_of_files.push_back(fname.Data());
       }
     }
   }
   return list_of_files;
}


void CombineEfficiency(TString variation="Nominal", TString varParton="Parton")
{ 

  AnalysisConstants::initConstants();
  TString baseInputDir = "/afs/cern.ch/work/g/gbakas/public/HEP-NTUA/";
  //TString baseInputDir = "/Users/georgebakas/Documents/HEP-NTUA_ul/";
  baseInputDir = TString::Format("%s/VariationHandling/", baseInputDir.Data());

  TString outFileDir = TString::Format("%sEfficiencyCombined/%s", baseInputDir.Data(), variation.Data());
  // Define formats for Figures and Latex file
  const TString ForVal = "%1.6f";
  const TString ForUnc = "%1.6f";
  const TString ForWei = "%1.3f";
  const TString ForRho = "%1.2f";
  const TString ForPul = ForRho;
  const TString ForUni = "pb";
  //std::vector<TString> variation_dirs = {"Nominal", "JES", "bTagVariation", "SystematicsFiles", "PSWeights", "PDFWeights", "ScaleWeights"};

  static const Int_t NumEst = AnalysisConstants::years.size();
  TString NamEst[NumEst];
  for (unsigned int i = 0; i < AnalysisConstants::years.size(); i++)
  {
    NamEst[i] = AnalysisConstants::years[i];
  }
  static const Int_t NumUnc = 1; //AnalysisConstants::variations.size() + 1;
  TString NamUnc[NumUnc];
  NamUnc[0] = "Stat";
  static const Int_t NumObs = 1;

  //loop on all variables
  for (unsigned int var = 0; var < AnalysisConstants::unfoldingVariables.size(); var++)
  {
    TString NamObs[NumObs];
    static const Int_t LenXEst = NumEst * (NumUnc + 1);
    Double_t XEst[LenXEst];

    Int_t IWhichObs[] = {0, 0, 0, 0};
    NamObs[0] = AnalysisConstants::unfoldingVariables[var];


    TString variable = AnalysisConstants::unfoldingVariables[var];
    // use 2018 as reference
    std::vector<TString> variationFiles_Hadronic = listFiles(TString::Format("%s/2016_postVFP/Responses%s/",
                                            baseInputDir.Data(), variation.Data()), "TTToHadronic");
    
    std::vector<TString> variationFiles_SemiLeptonic = listFiles(TString::Format("%s/2016_postVFP/Responses%s/",
                                            baseInputDir.Data(), variation.Data()), "TTToSemiLeptonic");
    
    std::vector<TString> variationFiles_Dilepton = listFiles(TString::Format("%s/2016_postVFP/Responses%s/",
                                            baseInputDir.Data(), variation.Data()), "TTTo2L2Nu");

    std::cout << "variable: " << variable << std::endl;
    std::cout << variationFiles_Hadronic.size() << std::endl;
    std::cout << variationFiles_SemiLeptonic.size() << std::endl;
    std::cout << variationFiles_Dilepton.size() << std::endl;

    // loop on all files (same for each year)
    for (unsigned int jvar=0; jvar<variationFiles_Hadronic.size(); jvar++)
    {

      std::vector<TH1F *> originalHistograms_numerator;
      std::vector<TH1F *> originalHistograms_denominator;
      //std::vector<TH1F *> uncHistograms;
      std::vector<TGraph *> weightGraphs;
      TFile *outFile = TFile::Open(TString::Format("%s/CombEfficiency%s_%s_%s",
                                    outFileDir.Data(),
                                    varParton.Data(),
                                    variable.Data(), variationFiles_Hadronic[jvar].Data()), "RECREATE");
      
      //find each file per year
      for (unsigned int y = 0; y < AnalysisConstants::years.size(); y++)
      {
        NamEst[y] = AnalysisConstants::years[y];
        TFile *inf_had = TFile::Open(TString::Format("%s/%s/Responses%s/%s",
                                  baseInputDir.Data(),
                                  AnalysisConstants::years[y].Data(),
                                  variation.Data(),
                                  variationFiles_Hadronic[jvar].Data()));

        TFile *inf_sem = TFile::Open(TString::Format("%s/%s/Responses%s/%s",
                                  baseInputDir.Data(),
                                  AnalysisConstants::years[y].Data(),
                                  variation.Data(),
                                  variationFiles_SemiLeptonic[jvar].Data()));

        TFile *inf_dil = TFile::Open(TString::Format("%s/%s/Responses%s/%s",
                                  baseInputDir.Data(),
                                  AnalysisConstants::years[y].Data(),
                                  variation.Data(),
                                  variationFiles_Dilepton[jvar].Data()));

        TString tempVar;
        if (varParton.EqualTo("Parton")) tempVar = AnalysisConstants::partonVariables[var];
        else tempVar = AnalysisConstants::particleVariables[var];

        TEfficiency *eff_had = (TEfficiency*)inf_had->Get(TString::Format("Efficiency%s_%s",varParton.Data(), tempVar.Data()));
        TEfficiency *eff_sem = (TEfficiency*)inf_sem->Get(TString::Format("Efficiency%s_%s",varParton.Data(), tempVar.Data()));
        TEfficiency *eff_dil = (TEfficiency*)inf_dil->Get(TString::Format("Efficiency%s_%s",varParton.Data(), tempVar.Data()));

        TEfficiency *efficiency = (TEfficiency*)eff_had->Clone();
        *efficiency += (*eff_sem);
        *efficiency += (*eff_dil); 

        //Get numerator and denominator
        TH1F *numerator = (TH1F*)efficiency->GetPassedHistogram();
        numerator->Sumw2();
        TH1F *denominator = (TH1F*)efficiency->GetTotalHistogram();
        denominator->Sumw2();
        
        //Get binning
        float bins[numerator->GetNbinsX() + 1];
        for (int bin = 0; bin < numerator->GetNbinsX(); bin++)
        {
          bins[bin] = numerator->GetBinLowEdge(bin + 1);
        }
        bins[numerator->GetNbinsX()] = numerator->GetBinLowEdge(numerator->GetNbinsX() + 1);
        //f->Divide(numerator, denominator, 1., 1., "B");

        numerator->SetDirectory(0);
        numerator->SetName(TString::Format("%s_%s",
                                  numerator->GetName(),
                                  AnalysisConstants::years[y].Data()));

        denominator->SetName(TString::Format("%s_%s",
                                  denominator->GetName(),
                                  AnalysisConstants::years[y].Data()));

        originalHistograms_numerator.push_back(numerator);
        originalHistograms_denominator.push_back(denominator);

        inf_had->Close();
        inf_sem->Close();
        inf_dil->Close();
      }

    Float_t *bins = GetHistogramBins(originalHistograms_numerator[0]);

    TH1F *resultsHisto = new TH1F(TString::Format("combined_%s",
                                                  AnalysisConstants::unfoldingVariables[var].Data()),
                                  TString::Format("combined_%s",
                                                  AnalysisConstants::unfoldingVariables[var].Data()),
                                  originalHistograms_numerator[0]->GetNbinsX(),
                                  bins);
    
    cout<< "--------------------------------" << endl;
    TH1F *numerator_all = new TH1F("numerator_all", "numerator_all", originalHistograms_numerator[0]->GetNbinsX(), bins);
    TH1F *denominator_all = new TH1F("denominator_all", "denominator_all", originalHistograms_numerator[0]->GetNbinsX(), bins);

    for (int ibin = 1; ibin <=originalHistograms_numerator[0]->GetNbinsX(); ibin++)
    {
    
    float numerator_content = originalHistograms_numerator[0]->GetBinContent(ibin) + 
                                originalHistograms_numerator[1]->GetBinContent(ibin) +
                                originalHistograms_numerator[2]->GetBinContent(ibin) +
                                originalHistograms_numerator[3]->GetBinContent(ibin);
    float denominator_content = originalHistograms_denominator[0]->GetBinContent(ibin) + 
                                originalHistograms_denominator[1]->GetBinContent(ibin) +
                                originalHistograms_denominator[2]->GetBinContent(ibin) +
                                originalHistograms_denominator[3]->GetBinContent(ibin);

    float numerator_error = TMath::Sqrt(TMath::Power(originalHistograms_numerator[0]->GetBinError(ibin), 2) + 
                                TMath::Power(originalHistograms_numerator[1]->GetBinError(ibin), 2) +
                                TMath::Power(originalHistograms_numerator[2]->GetBinError(ibin), 2) +
                                TMath::Power(originalHistograms_numerator[3]->GetBinError(ibin), 2)) + 
                                  // I have to input here the correlation coefficients 
                                2*(AnalysisConstants::correlations[variation]).correlations[1]*originalHistograms_numerator[0]->GetBinError(ibin)*originalHistograms_numerator[1]->GetBinError(ibin) + 
                                2*(AnalysisConstants::correlations[variation]).correlations[2]*originalHistograms_numerator[0]->GetBinError(ibin)*originalHistograms_numerator[2]->GetBinError(ibin) + 
                                2*(AnalysisConstants::correlations[variation]).correlations[3]*originalHistograms_numerator[0]->GetBinError(ibin)*originalHistograms_numerator[3]->GetBinError(ibin) + 
                                2*(AnalysisConstants::correlations[variation]).correlations[6]*originalHistograms_numerator[1]->GetBinError(ibin)*originalHistograms_numerator[2]->GetBinError(ibin) + 
                                2*(AnalysisConstants::correlations[variation]).correlations[7]*originalHistograms_numerator[1]->GetBinError(ibin)*originalHistograms_numerator[3]->GetBinError(ibin) + 
                                2*(AnalysisConstants::correlations[variation]).correlations[11]*originalHistograms_numerator[2]->GetBinError(ibin)*originalHistograms_numerator[3]->GetBinError(ibin);
    
    float denominator_error = TMath::Sqrt(TMath::Power(originalHistograms_denominator[0]->GetBinError(ibin), 2) + 
                                TMath::Power(originalHistograms_denominator[1]->GetBinError(ibin), 2) +
                                TMath::Power(originalHistograms_denominator[2]->GetBinError(ibin), 2) +
                                TMath::Power(originalHistograms_denominator[3]->GetBinError(ibin), 2)) + 
                                  // I have to input here the correlation coefficients 
                                2*(AnalysisConstants::correlations[variation]).correlations[1]*originalHistograms_denominator[0]->GetBinError(ibin)*originalHistograms_denominator[1]->GetBinError(ibin) + 
                                2*(AnalysisConstants::correlations[variation]).correlations[2]*originalHistograms_denominator[0]->GetBinError(ibin)*originalHistograms_denominator[2]->GetBinError(ibin) + 
                                2*(AnalysisConstants::correlations[variation]).correlations[3]*originalHistograms_denominator[0]->GetBinError(ibin)*originalHistograms_denominator[3]->GetBinError(ibin) + 
                                2*(AnalysisConstants::correlations[variation]).correlations[6]*originalHistograms_denominator[1]->GetBinError(ibin)*originalHistograms_denominator[2]->GetBinError(ibin) + 
                                2*(AnalysisConstants::correlations[variation]).correlations[7]*originalHistograms_denominator[1]->GetBinError(ibin)*originalHistograms_denominator[3]->GetBinError(ibin) + 
                                2*(AnalysisConstants::correlations[variation]).correlations[11]*originalHistograms_denominator[2]->GetBinError(ibin)*originalHistograms_denominator[3]->GetBinError(ibin);
    
    numerator_all->SetBinContent(ibin, numerator_content);
    denominator_all->SetBinContent(ibin, denominator_content);
    numerator_all->SetBinError(ibin, numerator_error);
    denominator_all->SetBinError(ibin, denominator_error);
    
    }
    numerator_all->Divide(denominator_all);

    outFile->cd();
    numerator_all->Write("efficiency");

    outFile->Close();
    delete bins;
    delete resultsHisto;
    }
  }
}