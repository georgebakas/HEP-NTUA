#include "BASE.h"

#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TGraph.h>

//Combine cross section ta the fiducial level without using the Blue method 

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


void CombineFiducialMeasurements_massVariations()
{

  AnalysisConstants::initConstants();

  TString baseInputDir = "/afs/cern.ch/work/g/gbakas/public/HEP-NTUA/";
  baseInputDir = TString::Format("%s/VariationHandling/", baseInputDir.Data());

  // Define formats for Figures and Latex file
  const TString ForVal = "%1.6f";
  const TString ForUnc = "%1.6f";
  const TString ForWei = "%1.3f";
  const TString ForRho = "%1.2f";
  const TString ForPul = ForRho;
  const TString ForUni = "pb";
  std::vector<TString> variation_dirs = {"SystematicsFiles"};

  std::vector<TString> systematic_vars = {"mtop166p5", "mtop169p5", "mtop171p5", "mtop173p5", "mtop175p5", "mtop178p5"};


  const Int_t NumEst = 4;//AnalysisConstants::years.size();
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
  
  TFile *outFile = new TFile("testFile_massVariations.root", "RECREATE");
  for (unsigned int var = 0; var < AnalysisConstants::unfoldingVariables.size(); var++)
  {
    for (unsigned int isys = 0; isys < systematic_vars.size(); isys++)
    {
        NamObs[0] = AnalysisConstants::unfoldingVariables[var];
        std::vector<TH1F *> originalHistograms;
        std::vector<TH1F *> uncHistograms;
        std::vector<TGraph *> weightGraphs;
        TString variable = AnalysisConstants::unfoldingVariables[var];

        std::cout << "variable: " << variable << std::endl;
        for (unsigned int y = 0; y < AnalysisConstants::years.size(); y++)
        {
            NamEst[y] = AnalysisConstants::years[y];        
            TFile *file = TFile::Open(TString::Format("%s/%s/SystematicsFiles/FiducialMeasurement/SignalHistograms_%s_MassFitResults_SignalTemplates_%s.root",
                                                    baseInputDir.Data(),
                                                    AnalysisConstants::years[y].Data(),
                                                    variable.Data(),
                                                    systematic_vars[isys].Data()));

            TH1F *f = (TH1F *)file->Get(TString::Format("hSignal_%s",
                                                    variable.Data()));
            cout << AnalysisConstants::years[y]<< endl;
            for (int ibin=1; ibin<=f->GetNbinsX(); ibin++ )
            {
                cout<< f->GetBinContent(ibin) << " err: "<< f->GetBinError(ibin) <<endl;
            }
            cout<<"integral: "<< f->Integral()<<endl;
            //convert to differential cross section
            //f->Scale(1. / AnalysisConstants::luminositiesSR[AnalysisConstants::years[y]], "width");

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
            
        }
    

    Float_t *bins = GetHistogramBins(originalHistograms[0]);

    TH1F *resultsHisto = new TH1F(TString::Format("combined_%s_%s",
                                                AnalysisConstants::unfoldingVariables[var].Data(),
                                                systematic_vars[isys].Data()),
                                TString::Format("combined_%s_%s",
                                                AnalysisConstants::unfoldingVariables[var].Data(),
                                                systematic_vars[isys].Data()),
                                originalHistograms[0]->GetNbinsX(),
                                bins);

    cout<< "--------------------------------" << endl;
    for (int ibin = 1; ibin <=originalHistograms[0]->GetNbinsX(); ibin++)
    {
    resultsHisto -> SetBinContent(ibin, originalHistograms[0]->GetBinContent(ibin) + 
                                originalHistograms[1]->GetBinContent(ibin) +
                                originalHistograms[2]->GetBinContent(ibin) +
                                originalHistograms[3]->GetBinContent(ibin));

    resultsHisto -> SetBinError(ibin, TMath::Sqrt(TMath::Power(originalHistograms[0]->GetBinError(ibin), 2) + 
                                TMath::Power(originalHistograms[1]->GetBinError(ibin), 2) +
                                TMath::Power(originalHistograms[2]->GetBinError(ibin), 2) +
                                TMath::Power(originalHistograms[3]->GetBinError(ibin), 2)) + 
                                  // I have to input here the correlation coefficients 
                                2*(AnalysisConstants::correlations["Nominal"]).correlations[1]*originalHistograms[0]->GetBinError(ibin)*originalHistograms[1]->GetBinError(ibin) + 
                                2*(AnalysisConstants::correlations["Nominal"]).correlations[2]*originalHistograms[0]->GetBinError(ibin)*originalHistograms[2]->GetBinError(ibin) + 
                                2*(AnalysisConstants::correlations["Nominal"]).correlations[3]*originalHistograms[0]->GetBinError(ibin)*originalHistograms[3]->GetBinError(ibin) + 
                                2*(AnalysisConstants::correlations["Nominal"]).correlations[6]*originalHistograms[1]->GetBinError(ibin)*originalHistograms[2]->GetBinError(ibin) + 
                                2*(AnalysisConstants::correlations["Nominal"]).correlations[7]*originalHistograms[1]->GetBinError(ibin)*originalHistograms[3]->GetBinError(ibin) + 
                                2*(AnalysisConstants::correlations["Nominal"]).correlations[11]*originalHistograms[2]->GetBinError(ibin)*originalHistograms[3]->GetBinError(ibin));
    
    cout<< resultsHisto ->GetBinContent(ibin) << " with error "<<resultsHisto->GetBinError(ibin)<<endl;
    
    }
    
    outFile->cd();
    resultsHisto->Write();

    delete bins;
    delete resultsHisto;
    }
  }
}
