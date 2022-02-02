#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include <TParameter.h>

using namespace std;

using std::cin;
using std::cout;
using std::endl;

#include "TemplateConstants.h"
using namespace RooFit;

void GetDifference(TString year, TString weightType, TString inputFitFile)
{
    initFilesMapping();
    RooMsgService::instance().setSilentMode(kTRUE);
    RooMsgService::instance().setStreamStatus(0,kFALSE);
    RooMsgService::instance().setStreamStatus(1,kFALSE);
    gStyle->SetOptStat(0);

    //open Nominal mass file to get the Nbkg
    TFile *fitFile_nom = TFile::Open(TString::Format("%s/Nominal/MassFitResults_SignalTemplates_.root",year.Data()));
    RooFitResult  *fitResult_nom = (RooFitResult*)fitFile_nom->Get(TString::Format("fitResults_%s", year.Data()));
    RooWorkspace  *w_nom = (RooWorkspace*)fitFile_nom->Get("w");
    //float NQCD = Nbkg2Constants[TString::Format("Nbkg%s",year.Data())];
    //float NQCD_error = Nbkg2ConstantsErrors[TString::Format("Nbkg%s_error",year.Data())];
    RooRealVar *value_nom = (RooRealVar*)w_nom->var("nFitQCD_2b");
    float NQCD_nominal = value_nom->getVal();
    float NQCD_nominal_error = value_nom->getError();

    //open the file to get the Nbkg
    cout<<TString::Format("%s/%s/%s",year.Data(), weightType.Data(), inputFitFile.Data())<<endl;
    TFile *fitFile = TFile::Open(TString::Format("%s/%s/%s",year.Data(), weightType.Data(), inputFitFile.Data()));
    RooFitResult  *fitResult = (RooFitResult*)fitFile->Get(TString::Format("fitResults_%s", year.Data()));
    RooWorkspace  *w = (RooWorkspace*)fitFile->Get("w");
    //float NQCD = Nbkg2Constants[TString::Format("Nbkg%s",year.Data())];
    //float NQCD_error = Nbkg2ConstantsErrors[TString::Format("Nbkg%s_error",year.Data())];
    RooRealVar *value = (RooRealVar*)w->var("nFitQCD_2b");
    float NQCD = value->getVal();
    float NQCD_error = value->getError();

    float total = NQCD_nominal/NQCD;
    float total_error = TMath::Sqrt(TMath::Power(1/NQCD * NQCD_nominal_error,2) + 
                                    TMath::Power((NQCD_nominal/TMath::Power(NQCD,2))*NQCD_error,2));
    cout<<"----- weightType: "<<weightType<< " inputFitFile: "<< inputFitFile<< " -----"<<endl;
    cout<<"NQCD/NQCD_nominal: "<<  total <<endl;
    cout<<"NQCD/NQCD_nominal error "<< total_error <<endl;


    TFile *output_file = TFile::Open(TString::Format("%s/MassDifferences.root",year.Data()), "UPDATE");
    output_file->cd();
    TParameter<float> *p = new TParameter<float>(inputFitFile, total);
    TParameter<float> *p_e = new TParameter<float>(TString::Format("%s_error",inputFitFile.Data()), total_error);
    p->Write(inputFitFile);
    p_e->Write(TString::Format("%s_error",inputFitFile.Data()));
    output_file->Close();
}
