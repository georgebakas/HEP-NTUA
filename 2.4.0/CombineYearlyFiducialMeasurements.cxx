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


void CombineYearlyFiducialMeasurements(TString variation="Nominal")
{
  AnalysisConstants::initConstants();
  //gSystem->Load("libBlue.so");
  TString baseInputDir = "/afs/cern.ch/work/g/gbakas/public/HEP-NTUA/";
  baseInputDir = TString::Format("%s/VariationHandling/", baseInputDir.Data());

  TString outFileDir = TString::Format("%s/FiducialCombined/%s", baseInputDir.Data(), variation.Data());
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
  TString NamObs[NumObs];
  static const Int_t LenXEst = NumEst * (NumUnc + 1);
  Double_t XEst[LenXEst];

  Int_t IWhichObs[] = {0, 0, 0, 0};

  
  //loop on all variables
  for (unsigned int var = 0; var < AnalysisConstants::unfoldingVariables.size(); var++)
  {
    NamObs[0] = AnalysisConstants::unfoldingVariables[var];
    std::vector<TH1F *> originalHistograms;
    //std::vector<TH1F *> uncHistograms;
    std::vector<TGraph *> weightGraphs;

    TString variable = AnalysisConstants::unfoldingVariables[var];
    std::cout << "variable: " << variable << std::endl;

    // use 2018 as reference 
    std::vector<TString> variationFiles = listFiles(TString::Format("%s/2018/%s/FiducialMeasurement/", 
                                            baseInputDir.Data(), variation.Data()), variable.Data());

    
    // loop on all files (same for each year)
    for (unsigned int jvar=0; jvar<variationFiles.size(); jvar++)
    {
      TFile *outFile = TFile::Open(TString::Format("%s/%s/Comb_%s_%s.root", 
                                    outFileDir.Data(), 
                                    variation.Data(), 
                                    variationFiles[jvar].Data(),
                                    variable.Data()), "RECREATE");
      std::cout<<variationFiles[jvar]<<endl;
      //find each file per year
      for (unsigned int y = 0; y < AnalysisConstants::years.size(); y++)
      {
        NamEst[y] = AnalysisConstants::years[y];
        
        TFile *file = TFile::Open(TString::Format("%s/%s/%s/FiducialMeasurement/%s",
                                  baseInputDir.Data(),
                                  AnalysisConstants::years[y].Data(),
                                  variation.Data(),
                                  variationFiles[jvar].Data()));
        

        TH1F *f = (TH1F *)file->Get(TString::Format("Signal_%s",
                                                    variable.Data()));

        //convert to differential cross section
        f->Scale(1. / AnalysisConstants::luminositiesSR[AnalysisConstants::years[y]], "width");

        f->SetDirectory(0);
        f->SetName(TString::Format("%s_%s",
                                  f->GetName(),
                                  AnalysisConstants::years[y].Data()));
        originalHistograms.push_back(f);

        TGraph *g = new TGraph(f->GetNbinsX());
        g->SetName(TString::Format("weights_%s_%s",
                                  AnalysisConstants::unfoldingVariables[var].Data(),
                                  AnalysisConstants::years[y].Data()));

        weightGraphs.push_back(g);

        file->Close();
        /*
        for (unsigned int i = 0; i < variation_dirs.size(); i++)
        {
          std::vector<TString> variationFiles = listFiles(TString::Format("%s/%s/%s/FiducialMeasurement/", baseInputDir.Data(), AnalysisConstants::years[y].Data(), variation_dirs[i].Data()), variable.Data());
          for (unsigned int jvar=0; jvar<variationFiles.size(); jvar++)
          {
            NamUnc[i + 1] = AnalysisConstants::variations[i];
            TFile *file = TFile::Open(TString::Format("%s/%s/%s/FiducialMeasurement/%s",
                                                    baseInputDir.Data(),
                                                    AnalysisConstants::years[y].Data(),
                                                    variation_dirs[i].Data(),
                                                    variationFiles[jvar].Data()));
            f = (TH1F *)file->Get(TString::Format("Signal_%s",
                                                  AnalysisConstants::unfoldingVariables[var].Data()));
            f->Scale(1. / AnalysisConstants::luminositiesSR[AnalysisConstants::years[y]], "width");

            f->SetDirectory(0);
            f->SetName(TString::Format("%s_%s_%s",
                                      f->GetName(),
                                      AnalysisConstants::variations[i].Data(),
                                      AnalysisConstants::years[y].Data()));
            uncHistograms.push_back(f);
            file->Close();
          }
        } 
        */
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
        std::cout << std::endl;
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
      }

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