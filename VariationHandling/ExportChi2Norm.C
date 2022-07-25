#include "BASE.h"

#include "FinalResultsConstants.h"

#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "../CMS_plots/CMS_lumi.C"
#include "../CMS_plots/tdrstyle.C"

void ExportChi2Norm(TString partonParticle = "Parton")
{
    AnalysisConstants::initConstants();
    TString baseDir = "/Users/georgebakas/Documents/HEP-NTUA_ul/VariationHandling";
    //TString baseDir = "/afs/cern.ch/work/g/gbakas/public/HEP-NTUA/VariationHandling";
    
    TFile *resultsFile = TFile::Open(TString::Format("%s/FinalResults/results/%s_outputFileNorm.root",
                                                baseDir.Data(),
                                                partonParticle.Data()));
    
    TFile *chiSquareFiles = TFile::Open(TString::Format("%s/ChiSquare/results/%s_chiSquareResults_normalised.root",
                                                    baseDir.Data(),
                                                    partonParticle.Data()));

    //write to a txt file the Kolmogorov tests
    ofstream myfile;
    myfile.open (TString::Format("FinalResults/resultsChi2/%sAnalytical_resultsNorm.txt", partonParticle.Data()));
    myfile<<"------------"<< "\n";    
    myfile<<"Variable & "<<"NDOF & chi2 Theory & chi2 Theory AMC@NLO & P-value Theory & P-value AMC@NLO \\\\" <<"\n";

    for (int i = 0; i < AnalysisConstants::unfoldingVariables.size(); i++)
    {
        TString variable = AnalysisConstants::unfoldingVariables[i];
        TMatrixD *chiSquare = (TMatrixD *)chiSquareFiles->Get(TString::Format("chiSquare_%s_normalized",
                                                                        variable.Data()));
    TMatrixD *chiSquare_amc = (TMatrixD *)chiSquareFiles->Get(TString::Format("chiSquare_%s_amc_normalized",
                                                                            variable.Data()));

        TH1F *f = (TH1F *)resultsFile->Get(TString::Format("FinalResult_%s",
                                                    variable.Data()));

        Double_t chi2Theory = chiSquare->operator()(0, 0);
        Double_t chi2Theory_amcNlo = chiSquare_amc->operator()(0, 0);
        Int_t ndof = f->GetNbinsX()-1;

        Double_t p_valTheory = TMath::Prob(chi2Theory, ndof-1);
        Double_t p_valTheory_amcNlo = TMath::Prob(chi2Theory_amcNlo, ndof-1);

        myfile<<variable<<" & "<<ndof <<" & "<<chi2Theory<<" & "<<chi2Theory_amcNlo<<" & "<<p_valTheory <<" & "<< p_valTheory_amcNlo <<" \\\\ " <<"\n";
        myfile<<"\\hline"<< "\n";
    }
    myfile.close();
}