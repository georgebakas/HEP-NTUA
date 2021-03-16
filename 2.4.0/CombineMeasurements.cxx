#include "BASE.h"

#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TGraph.h>

#include "Blue.h"

void CombineMeasurements(TFile *outFile, bool normalized = false)
{
  gSystem->Load("libBlue.so");
  cout<<"here"<<endl;
  AnalysisConstants::initConstants();
  TString baseInputDir = "/afs/cern.ch/work/g/gbakas/public/HEP-NTUA/";
  baseInputDir = TString::Format("%s/Unfolding", baseInputDir.Data());
  std::cout<<"years SIZE: "<<  AnalysisConstants::years.size() << std::endl;

  // Define formats for Figures and Latex file
  const TString ForVal = "%1.6f";
  const TString ForUnc = "%1.6f";
  const TString ForWei = "%1.3f";
  const TString ForRho = "%1.2f";
  const TString ForPul = ForRho;
  const TString ForUni = "pb";

  std::cout << "Skata" << std::endl;
  static const Int_t NumEst = AnalysisConstants::years.size();
  
  TString NamEst[NumEst];
  for (unsigned int i = 0; i < AnalysisConstants::years.size(); i++)
  {
    NamEst[i] = AnalysisConstants::years[i];
  }
  static const Int_t NumUnc = AnalysisConstants::variations.size() + 1;
  TString NamUnc[NumUnc];
  NamUnc[0] = "Stat";
  static const Int_t NumObs = 1;
  TString NamObs[NumObs];
  static const Int_t LenXEst = NumEst * (NumUnc + 1);
  Double_t XEst[LenXEst];

  Int_t IWhichObs[NumEst] = {0, 0, 0, 0};

  for (unsigned int var = 0; var < AnalysisConstants::unfoldingVariables.size(); var++)
  {
    NamObs[0] = AnalysisConstants::unfoldingVariables[var];
    std::cout<< "size: "<< AnalysisConstants::unfoldingVariables.size() << endl;
    std::cout<< "name: "<<NamObs[0]<<std::endl;

    std::vector<TH1F *> originalHistograms;
    std::vector<TH1F *> uncHistograms;
    std::vector<TGraph *> weightGraphs;
    TString variable = AnalysisConstants::unfoldingVariables[var];
    std::cout << "variable: " << variable.Data() << std::endl;

    TString partonParticleStr = "Parton";
    for (unsigned int y = 0; y < AnalysisConstants::years.size(); y++)
    {
      NamEst[y] = AnalysisConstants::years[y];
      TString fileName = TString::Format("%s/%s/%sMeasurements/Data%s/OutputFile.root",
                                                baseInputDir.Data(),
                                                AnalysisConstants::years[y].Data(),
                                                partonParticleStr.Data(),
                                                (normalized ? "_Norm" : ""));
      std::cout << "year: " << NamEst[y] << endl;
      std::cout << "fileName: "<< fileName << std::endl;
      TFile *file = TFile::Open(TString::Format("%s/%s/%sMeasurements/Data%s/OutputFile.root",
                                                baseInputDir.Data(),
                                                AnalysisConstants::years[y].Data(),
                                                partonParticleStr.Data(),
                                                (normalized ? "_Norm" : "")));
      std::cout << file->GetName() <<std::endl;
      std::cout << TString::Format("hUnfold_%s",
                                   variable.Data() )<< std::endl;
                                   //(normalized ? "_normalized" : ""))
      TH1F *f = (TH1F *)file->Get(TString::Format("hUnfold_%s",
                                                  variable.Data()));
      f->SetDirectory(0);
      f->SetName(TString::Format("%s_%s",
                                 f->GetName(),
                                 AnalysisConstants::years[y].Data()));
      
      originalHistograms.push_back(f);
      TGraph *g = new TGraph(f->GetNbinsX());
      g->SetName(TString::Format("weights_%s%s_%s",
                                 AnalysisConstants::unfoldingVariables[var].Data(),
                                 (normalized ? "_normalized" : ""),
                                 AnalysisConstants::years[y].Data()));
      weightGraphs.push_back(g);

      file->Close();

      for (unsigned int i = 0; i < AnalysisConstants::variations.size(); i++)
      {
        NamUnc[i + 1] = AnalysisConstants::variations[i];
        file = TFile::Open(TString::Format("%s/%s/%s%s/UnfoldingResults_%s.root",
                                           baseInputDir.Data(),
                                           AnalysisConstants::years[y].Data(),
                                           AnalysisConstants::variations[i].Data(),
                                           AnalysisConstants::currentlyWorkingDirectory[AnalysisConstants::years[y]].Data(),
                                           AnalysisConstants::years[y].Data()));
        f = (TH1F *)file->Get(TString::Format("unfoldedHistogram_%s%s",
                                              AnalysisConstants::unfoldingVariables[var].Data(),
                                              (normalized ? "_normalized" : "")));
        std::cout << TString::Format("unfoldedHistogram_%s%s",
                                     AnalysisConstants::unfoldingVariables[var].Data(),
                                     (normalized ? "_normalized" : ""))
                  << std::endl;
        f->SetDirectory(0);
        f->SetName(TString::Format("%s_%s_%s%s",
                                   f->GetName(),
                                   AnalysisConstants::variations[i].Data(),
                                   AnalysisConstants::years[y].Data(),
                                   (normalized ? "_normalized" : "")));
        std::cout << TString::Format("%s_%s_%s%s",
                                     f->GetName(),
                                     AnalysisConstants::variations[i].Data(),
                                     AnalysisConstants::years[y].Data(),
                                     (normalized ? "_normalized" : ""))
                  << std::endl;
        uncHistograms.push_back(f);
        file->Close();
      }
    }

    Float_t *bins = GetHistogramBins(originalHistograms[0]);

    TH1F *resultsHisto = new TH1F(TString::Format("combined_%s%s",
                                                  AnalysisConstants::unfoldingVariables[var].Data(),
                                                  (normalized ? "_normalized" : "")),
                                  TString::Format("combined_%s%s",
                                                  AnalysisConstants::unfoldingVariables[var].Data(),
                                                  (normalized ? "_normalized" : "")),
                                  originalHistograms[0]->GetNbinsX(),
                                  bins);

    std::cout << resultsHisto->GetName() << std::endl;

    for (int bin = 1; bin <= originalHistograms[0]->GetNbinsX(); bin++)
    {
      for (unsigned int h = 0; h < originalHistograms.size(); h++)
      {
        XEst[(NumUnc + 1) * h] = originalHistograms[h]->GetBinContent(bin);
        XEst[(NumUnc + 1) * h + 1] = originalHistograms[h]->GetBinError(bin);
        std::cout << "Value: " << (NumUnc + 1) * h << " " << XEst[(NumUnc + 1) * h] << std::endl;
        std::cout << "Statistical uncertainty: " << (NumUnc + 1) * h + 1 << " " << XEst[(NumUnc + 1) * h + 1] << std::endl;
        for (unsigned int i = 0; i < AnalysisConstants::variations.size(); i++)
        {
          XEst[(NumUnc + 1) * h + (i + 2)] = TMath::Abs((uncHistograms[(h * AnalysisConstants::variations.size()) + i]->GetBinContent(bin) - originalHistograms[h]->GetBinContent(bin)));
          std::cout << AnalysisConstants::variations[i] << " " << (NumUnc + 1) * h + (i + 2) << " " << XEst[(NumUnc + 1) * h + (i + 2)] << std::endl;
        }
        std::cout << std::endl;
      }

      // Construct Object

      Blue *myBlue = new Blue(NumEst, NumUnc, NumObs, &IWhichObs[0]);
      myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForUni);
      myBlue->FillNamObs(&NamObs[0]);
      myBlue->FillNamEst(&NamEst[0]);
      myBlue->FillNamUnc(&NamUnc[0]);
      // Fill estimates
      Int_t ind = 0;
      for (Int_t i = 0; i < NumEst; i++)
      {
        myBlue->FillEst(i, &XEst[ind]);
        ind = ind + NumUnc + 1;
      }

      // Fill correlations
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
      //myBlue->PrintCor();
      myBlue->PrintCov();
      myBlue->PrintCov();
      //myBlue->PrintStatus();
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
      //result->Print();
      //unc->Print();
      rho->Print();
      cov->Print();
      //myBlue->PrintResult();
      resultsHisto->SetBinContent(bin, result->operator()(0, 0));
      resultsHisto->SetBinError(bin, unc->operator()(0, 0));
      std::cout << "results:" << std::endl;

      TMatrixD *weights = new TMatrixD(NumEst, NumObs);
      myBlue->GetWeight(weights);
      weights->Print();

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

    std::cout<< "S K A T A RE MALAKA!!!!!"<<std::endl;
  }

  
}