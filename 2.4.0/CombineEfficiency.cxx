#include "BASE.h"

#include <iostream>

#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"

#include "Blue.h"

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
  //gSystem->Load("libBlue.so");
  TString baseInputDir = "/afs/cern.ch/work/g/gbakas/public/HEP-NTUA/";
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

      std::vector<TH1F *> originalHistograms;
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
        //Efficiency
        TH1F *f = new TH1F("Efficiency", "Efficiency", numerator->GetNbinsX(), bins);
        f->Divide(numerator, denominator, 1., 1., "B");

        f->SetDirectory(0);
        f->SetName(TString::Format("%s_%s",
                                  f->GetName(),
                                  AnalysisConstants::years[y].Data()));
        originalHistograms.push_back(f);

        TGraph *g = new TGraph(f->GetNbinsX());
        g->SetName(TString::Format("weights_Efficiency_%s_%s",
                                  variable.Data(),
                                  AnalysisConstants::years[y].Data()));
        weightGraphs.push_back(g);

        inf_had->Close();
        inf_sem->Close();
        inf_dil->Close();
       
      }


    Float_t *bins = GetHistogramBins(originalHistograms[0]);

    TH1F *resultsHisto = new TH1F(TString::Format("combined_%s",
                                                  AnalysisConstants::unfoldingVariables[var].Data()),
                                  TString::Format("combined_%s",
                                                  AnalysisConstants::unfoldingVariables[var].Data()),
                                  originalHistograms[0]->GetNbinsX(),
                                  bins);
    for (int bin = 1; bin <= originalHistograms[0]->GetNbinsX(); bin++)
    {
      std::cout << "Bin: " << bin << std::endl;
      for (unsigned int h = 0; h < originalHistograms.size(); h++)
      {
        XEst[(NumUnc + 1) * h] = originalHistograms[h]->GetBinContent(bin);
        XEst[(NumUnc + 1) * h + 1] = originalHistograms[h]->GetBinError(bin);
        std::cout << "Value: " << (NumUnc + 1) * h << " " << XEst[(NumUnc + 1) * h] << std::endl;
        std::cout << "Statistical uncertainty: " << (NumUnc + 1) * h + 1 << " " << XEst[(NumUnc + 1) * h + 1] << std::endl;
        /*
        for (unsigned int i = 0; i < AnalysisConstants::variations.size(); i++)
        {
          XEst[(NumUnc + 1) * h + (i + 2)] = TMath::Abs((uncHistograms[(h * AnalysisConstants::variations.size()) + i]->GetBinContent(bin)) - originalHistograms[h]->GetBinContent(bin));
          std::cout << AnalysisConstants::variations[i] << " " << (NumUnc + 1) * h + (i + 2) << " " << XEst[(NumUnc + 1) * h + (i + 2)] << std::endl;
        }
        */

      }

      Blue *myBlue = new Blue(NumEst, NumUnc, NumObs, &IWhichObs[0]);
      myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForUni);
      myBlue->FillNamObs(&NamObs[0]);
      myBlue->FillNamEst(&NamEst[0]);
      myBlue->FillNamUnc(&NamUnc[0]);

      Int_t ind = 0;
      for (Int_t i = 0; i < NumEst; i++)
      {
        myBlue->FillEst(i, &XEst[ind]);
        ind = ind + NumUnc + 1;
      }
      /*
      for (int k = 0; k < NumUnc; k++)
      {
        if (k == 0)
        {
          myBlue->FillCor(k, 0.0);
        }
        else
        {
          myBlue->FillCor(k, &(AnalysisConstants::correlations[AnalysisConstants::variations[k - 1]].correlations[0]));
        }
      } */

      myBlue->FixInp();
      //myBlue->PrintCov();
      //myBlue->PrintCov();
      myBlue->Solve();

      TMatrixD *result = new TMatrixD(NumObs, NumUnc + 1);
      TMatrixD *unc = new TMatrixD(NumObs, 1);
      TMatrixD *cov = new TMatrixD(NumEst, NumEst);
      TMatrixD *rho = new TMatrixD(NumEst, NumEst);
      myBlue->PrintResult();
      myBlue->GetResult(result);
      myBlue->GetUncert(unc);
      myBlue->GetCov(cov);
      myBlue->GetRho(rho);
      //myBlue->SetPrintLevel(1);
      //myBlue->InspectResult();

      //rho->Print();
      //std::cout << "Covariance" << std::endl;
      //cov->Print();

      resultsHisto->SetBinContent(bin, result->operator()(0, 0));
      resultsHisto->SetBinError(bin, unc->operator()(0, 0));
      std::cout << "results:" << std::endl;

      TMatrixD *weights = new TMatrixD(NumEst, NumObs);
      myBlue->GetWeight(weights);
      //weights->Print();

      for (unsigned int i = 0; i < weightGraphs.size(); i++)
      {
        weightGraphs[i]->SetPoint(bin - 1, resultsHisto->GetBinCenter(bin), weights->operator()(i, 0));
      }
      delete result;
      delete unc;
      delete myBlue;
      delete rho;
    }

    outFile->cd();
    resultsHisto->Write();
    for (unsigned int i = 0; i < weightGraphs.size(); i++)
    {
      weightGraphs[i]->Write();
    }

    delete bins;
    delete resultsHisto;
    }
  }
}
