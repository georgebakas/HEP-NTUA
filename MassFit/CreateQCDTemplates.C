#include "QCDBkgConstants.h"

void CreateQCDTemplates(TString year, TString CUT = "")
{
  initQCDParams();
  gROOT->ForceStyle();

  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);

  RooRealVar *kMassScale = new RooRealVar("kMassScale","kMassScale",1.0,0.5,1.5);
  RooRealVar *kMassResol = new RooRealVar("kMassResol","kMassResol",1.0,0.5,1.5);
  kMassScale->setConstant(kTRUE);
  kMassResol->setConstant(kTRUE);

  TString VAR,TAG;
  float XMIN,XMAX;

  VAR = "mTop";
  XMIN = 50.;
  XMAX = 300.;

  RooWorkspace *w = new RooWorkspace("w","workspace");

  //---- define observable ------------------------
  RooRealVar *x = new RooRealVar("mTop","mTop",XMIN,XMAX);
  w->import(*x);
  //---- first do the data template ---------------
  TFile *infData = TFile::Open(TString::Format("%s/Histo_Data_%s_100.root", year.Data(),year.Data()));
  TH1F *hData = (TH1F*)infData->Get("hWt_mTop_0btag_expYield");
  cout<<"CR Entries from data: "<<hData->GetEntries()<<endl;
  //because of contamination we need to substract the ttbar from the CR.
  //we do that by extracting the 0btag th1 using nominal or mtt mc
  TFile *infTTMC = TFile::Open(TString::Format("%s/Histo_TT_NominalMC_100.root",year.Data()));  //nominal
  //TFile *infTTMC = TFile::Open(TString::Format("%s/Histo_TT_Mtt-700toInf_100.root",year.Data())); //mtt
  TH1F *hCR_MC = (TH1F*)infTTMC->Get("hWt_mTop_0btag_expYield");

  TFile *infBkg = TFile::Open(TString::Format("%s/Histo_SubdominantBkgs_100.root",year.Data()));
  TH1F *hCR_MCSubdominant = (TH1F*)infBkg->Get("hWt_mTop_0btag_expYield");
  //TH1F *hDataBefore = (TH1F*)hData->Clone("test");

  //hDataBefore->SetLineColor(kRed);
  hData->Add(hCR_MCSubdominant,-1);
  hData->Add(hCR_MC,-1);

  hData->Rebin(2);

  cout<<"After Substr: "<<hData->GetEntries()<<endl;
  /*
  hData->Draw();
  hDataBefore->Draw("same");
  hCR_MCSubdominant->SetLineColor(kGreen);
  hCR_MC->SetLineColor(kMagenta);
  hCR_MC->Draw("same");
  hCR_MCSubdominant->Draw("same");
  return;
  */
  RooDataHist *roohData = new RooDataHist("roohistData","roohistData",RooArgList(*x),hData);
  TH1F *hDataJet = (TH1F*)infData->Get("hWt_jetMassSoftDrop_0btag");
  TH1F *hDataEvt = (TH1F*)infData->Get("hWt_mJJ_0btag");
  RooRealVar *fBkgJet = new RooRealVar("fBkgJet","fBkgJet",hDataJet->Integral()/hData->Integral());
  RooRealVar *fBkgEvt = new RooRealVar("fBkgEvt","fBkgEvt",hDataEvt->Integral()/hData->Integral());
  //---- QCD -----------------------------------
  RooRealVar bQCD0("qcd_b0","qcd_b0",qcdParams[year]["qcd_b0"],0,2.);
  RooRealVar bQCD1("qcd_b1","qcd_b1",qcdParams[year]["qcd_b1"],0,2.);
  RooRealVar bQCD2("qcd_b2","qcd_b2",qcdParams[year]["qcd_b2"],0,2.);
  RooRealVar bQCD3("qcd_b3","qcd_b3",qcdParams[year]["qcd_b3"],0,2.);
  RooRealVar bQCD4("qcd_b4","qcd_b4",qcdParams[year]["qcd_b4"],0,2.);
  RooBernstein qcd1("qcd_brn","qcd_brn",*x,RooArgList(bQCD0,bQCD1,bQCD2,bQCD3,bQCD4));


  RooRealVar mQCD("qcd_mean" ,"qcd_mean",qcdParams[year]["qcd_mean"],130,300);
  RooRealVar sQCD("qcd_sigma","qcd_sigma",qcdParams[year]["qcd_sigma"],10,200);
  RooGaussian qcd2("qcd_gaus" ,"qcd_gaus",*x,mQCD,sQCD);


  RooRealVar mWqcd("qcd_meanW", "qcd_meanW", 77, 70, 90);
  RooRealVar sWqcd("qcd_sigmaW", "qcd_sigmaW", 20, 0, 40);
  RooGaussian qcd3("qcd_gausW", "qcd_gausW", *x, mWqcd, sWqcd);

  RooRealVar fqcd1("qcd_f1","qcd_f1",qcdParams[year]["qcd_f1"],0,1);
  RooRealVar fqcd2("qcd_f2","qcd_f2",0.5,0,1);

  RooAddPdf *qcd = new RooAddPdf("qcd_pdf","qcd_pdf",RooArgList(qcd1,qcd2), RooArgList(fqcd1));
  qcd->Print();

  //---- plots ---------------------------------------------------
  TCanvas *canQCD = new TCanvas("Template_QCD_"+CUT,"Template_QCD_"+CUT,900,600);

  RooFitResult *res = qcd->fitTo(*roohData,RooFit::Save());
  res->Print();
  RooPlot *frameQCD = x->frame();
  roohData->plotOn(frameQCD);
  qcd->plotOn(frameQCD);
  qcd->plotOn(frameQCD,RooFit::Components("qcd_brn"),RooFit::LineColor(kRed),RooFit::LineWidth(2),RooFit::LineStyle(2));
  qcd->plotOn(frameQCD,RooFit::Components("qcd_gaus"),RooFit::LineColor(kGreen+1),RooFit::LineWidth(2),RooFit::LineStyle(2));
  //qcd->plotOn(frameQCD,RooFit::Components("qcd_gausW"),RooFit::LineColor(kOrange+1),RooFit::LineWidth(2),RooFit::LineStyle(2));
  frameQCD->GetXaxis()->SetTitle("m_{t} (GeV)");
  frameQCD->Draw();
  gPad->Update();
  canQCD->Print(TString::Format("%s/plots/templateResults/"+TString(canQCD->GetName())+".pdf", year.Data()));

  RooArgSet *parsQCD = (RooArgSet*)qcd->getParameters(roohData);
  parsQCD->setAttribAll("Constant",true);

  //we do that by extracting the 0b
  w->import(*qcd);
  w->writeToFile(TString::Format("%s/templates_QCD_"+CUT+"100.root", year.Data()));

}
