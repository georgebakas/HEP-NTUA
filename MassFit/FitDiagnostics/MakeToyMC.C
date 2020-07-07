#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include <iostream>
#include <fstream>
using namespace RooFit;

#include "../TemplateConstants.h"

TString globYear;

void MakeToyMC(TString year="2016")
{
  initFilesMapping();
  globYear = year;

  // Name of variable
	RooRealVar mTop("mTop", "mTop", 50., 300.);
  Float_t sig_factor = ttbarSigStrength[globYear];

  //this is the output file
  TFile *resultsFile = new TFile(TString::Format("%s/ToyMC_Data.root", globYear.Data()), "RECREATE");


  //here i need to take the initial integrals for all my processes
  //So 1st thing is to open the root files

  //the processes in the ttbar are:
  //ttbar
  TFile *infTT = TFile::Open(TString::Format("../%s/Histo_TT_NominalMC_100.root",globYear.Data()));
  //data
  TFile *infData = TFile::Open(TString::Format("../%s/Histo_Data_%s_100.root",globYear.Data(),globYear.Data()));
  //qcd
  TFile *infQCD = TFile::Open(TString::Format("../%s/Histo_QCD_HT300toInf_100.root",globYear.Data()));
  //subdominant
  TFile *infSub = TFile::Open(TString::Format("../%s/Histo_SubdominantBkgs_100.root",globYear.Data()));

  float initialSignalStrength = 0.5;
  float step = 0.1;

  for (int i = 0; i <= 10; i++)
  {

  float signalStrength = initialSignalStrength + i * step;
  TH1F *hTT_in = (TH1F*)infTT->Get("hWt_mTop_2btag");
  TH1F *hQCD_in = (TH1F*)infQCD->Get("hWt_mTop_2btag");
  TH1F *hData_in = (TH1F*)infData->Get("hWt_mTop_2btag");
  TH1F *hSub_in = (TH1F*)infSub->Get("hWt_mTop_2btag");

  //get the initial yields:
  Float_t nTT_in = hTT_in->Integral();
  Float_t nData_in = hData_in->Integral();
  Float_t nSub_in = hSub_in->Integral();

  //scale qcd with data (k-factor)
  //Float_t nQCD_in = hQCD_in->Integral();
  Float_t nQCD_in = nData_in - nTT_in - nSub_in;
  hQCD_in -> Scale(nQCD_in/(hQCD_in -> Integral()));

  //Create the fit model from the templates
  TFile *fTemplatesBkg = TFile::Open(TString::Format("../%s/templates_Bkg_100.root", year.Data()));
  TFile *fTemplatesSig = TFile::Open(TString::Format("../%s/templates_Sig_100.root", year.Data()));
  RooWorkspace *wTemplatesBkg = (RooWorkspace *)fTemplatesBkg->Get("w");
  RooWorkspace *wTemplatesSig = (RooWorkspace *)fTemplatesSig->Get("w");
  RooRealVar *x = (RooRealVar*)wTemplatesSig->var("mTop");

  RooAbsPdf *pdf_bkg_2b = (RooAbsPdf *)wTemplatesBkg->pdf("bkg_pdf_2btag");
  RooAbsPdf *pdf_qcd_2b = (RooAbsPdf *)wTemplatesBkg->pdf("qcd_pdf");
  RooAbsPdf *pdf_signal_2b = (RooAbsPdf *)wTemplatesSig->pdf("ttbar_pdf_2btag");

  RooRealVar *kMassScale = (RooRealVar*)wTemplatesSig->var("kMassScale");
  RooRealVar *kMassResol = (RooRealVar*)wTemplatesSig->var("kMassResol");
  kMassScale->setConstant(false);
  kMassResol->setConstant(false);

  float sP = 0.0025;
  RooRealVar *kQCD2b_0 = new RooRealVar("kQCD_2b","kQCD_2b", sP);//, min, max);
  kQCD2b_0->setConstant(false);
  RooFormulaVar qcdCor_2b("qcdCor_2b","(1+@0*@1)",RooArgList(*x,*kQCD2b_0));
  //corrected QCD
  RooEffProd pdf_qcdCor_2b("qcdCor_pdf_2b","qcdCor_pdf_2b",*pdf_qcd_2b,qcdCor_2b);

  RooRealVar *nFitBkg2b = new RooRealVar("nFitBkg_2b", "nFitBkg_2b", nSub_in, 0, 10e+3);
  RooRealVar *nFitQCD2b = new RooRealVar("nFitQCD_2b", "nFitQCD_2b", nQCD_in, 0, 10e+4);
  RooRealVar *nFitSig2b = new RooRealVar("nFitSig2b", "nFitSig2b", signalStrength * nTT_in, 100, 10e+4);

  RooAddPdf model("model_2b", "model_2b",
                      RooArgList(*pdf_signal_2b, pdf_qcdCor_2b, *pdf_bkg_2b),
                      RooArgList(*nFitSig2b, *nFitQCD2b, *nFitBkg2b));

  // generated yields for bkgs
	Float_t nTT_gen,nQCD_gen,nSub_gen,nData_gen;

	TRandom3 ran;
  // ----- Histograms to save the results ------------

	//tt histograms
  TH1F *ttbarHisto = new TH1F(TString::Format("ttbarHisto_%.2f_%s", signalStrength, year.Data()),
                              TString::Format("ttbarHisto_%.2f_%s", signalStrength, year.Data()),
                              75, -0.1 * signalStrength * nTT_in, 0.1 * signalStrength * nTT_in);

  TH1F *ttbarPullHisto = new TH1F(TString::Format("ttbarPullHisto_%.2f_%s", signalStrength, year.Data()),
                                  TString::Format("ttbarPullHisto_%.2f_%s", signalStrength, year.Data()),
                                  100, -10, 10);

  TH1F *ttbarErrorHisto = new TH1F(TString::Format("ttbarErrorHisto_%.2f_%s", signalStrength, year.Data()),
                                   TString::Format("ttbarErrorHisto_%.2f_%s", signalStrength, year.Data()),
                                   100, 80, 180);


  for (int i = 0; i < 1; i++)
    {
      std::cout << "runNo: " << i << std::endl;
      ran.SetSeed(0);

      nSub_gen = ran.Poisson(nSub_in);
      nQCD_gen = ran.Poisson(nQCD_in);

      nData_gen = nSub_gen + nQCD_gen + signalStrength * nTT_in;

      RooAddPdf pseudodata("pseudodata", "pseudodata",
                           RooArgList(*pdf_signal_2b, pdf_qcdCor_2b, *pdf_bkg_2b),
                           RooArgList(RooFit::RooConst(signalStrength * nTT_in),
                                      RooFit::RooConst(nQCD_gen),
                                      RooFit::RooConst(nSub_gen)));

      RooDataSet *data_gen = pseudodata.generate(mTop, nData_gen);

      RooDataHist *data = data_gen->binnedClone();

      RooFitResult *fit_result = model.fitTo(*data);//, RooFit::PrintLevel(-1)); //, RooFit::Save(), RooFit::Extended(kTRUE));
      fit_result->Print();

      /*RooPlot *frame = mTop.frame();
      data_gen->plotOn(frame);
      pseudodata.plotOn(frame, RooFit::Components("ttbar_pdf"), RooFit::LineColor(kRed));
      pseudodata.plotOn(frame, RooFit::Components("qcd_pdf"), RooFit::LineColor(kGreen));
      pseudodata.plotOn(frame, RooFit::Components("sub_pdf"), RooFit::LineColor(kMagenta));
      frame->Draw();
      dataHisto->Draw("SAME");*/

      //TFile *f = new TFile("output.root", "RECREATE");

      //data->Write();
      //f->Close();
      ttbarHisto->Fill(nFitSig2b->getVal() - signalStrength * nTT_in);
      ttbarPullHisto->Fill((nFitSig2b->getVal() - signalStrength * nTT_in) / nFitSig2b->getError());
      ttbarErrorHisto->Fill(nFitSig2b->getError());

      nFitBkg2b->setVal(nSub_gen);
      nFitQCD2b->setVal(nQCD_gen);
      nFitSig2b->setVal(signalStrength * nTT_in);
    }
    resultsFile->cd();
    ttbarHisto->Write();
    ttbarPullHisto->Write();
    ttbarErrorHisto->Write();
  }
  resultsFile->Close();

}
