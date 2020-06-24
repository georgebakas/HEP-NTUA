#include "QCDBkgConstants.h"

void CreateBkgTemplates(TString year, TString CUT = "")
{
  initQCDParams();
  gROOT->ForceStyle();

  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);

  TFile *infBkg;
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

  TFile *infBkg_temp = TFile::Open(TString::Format("%s/Histo_SubdominantBkgs_100.root",year.Data()));
  TH1F *hCR_MCSubdominant = (TH1F*)infBkg_temp->Get("hWt_mTop_0btag_expYield");
  TH1F *hDataBefore = (TH1F*)hData->Clone("test");

  //hDataBefore->SetLineColor(kRed);
  hData->Add(hCR_MC,-1);
  hData->Add(hCR_MCSubdominant,-1);
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

  w->import(*qcd);
  w->import(*fBkgJet);
  w->import(*fBkgEvt);
  //TH1F *hMeanTop[2],*hSigmaTop[2],*hMeanW[2],*hSigmaW[2];

  float LUMI(0);
  if (year.EqualTo("2016")) LUMI = 35920;
  else if (year.EqualTo("2017")) LUMI = 41530;
  else if (year.EqualTo("2018")) LUMI = 59740;

  for(int icat=0;icat<3;icat++) {
    if (icat==1) continue;
    infBkg = TFile::Open(TString::Format("%s/Histo_SubdominantBkgs_100.root",year.Data()));


    TString CAT = TString::Format("%dbtag",icat);
    TAG = CUT+"_"+CAT;

    //---- do the bkg templates -------------
    TH1F *hBkg = (TH1F*)infBkg->Get("hWt_"+VAR+TAG+"_expYield");
    cout<<"icat "<<icat<<": "<<hBkg->Integral()<<endl;
    RooDataHist *roohBkg = new RooDataHist("roohistBkg","roohistBkg",RooArgList(*x),hBkg);

    RooRealVar mW("bkg_meanW_"+CAT,"meanW_"+CAT,80,70,90);
    RooRealVar sW("bkg_sigmaW_"+CAT,"sigmaW_"+CAT,5,0,15);
    RooFormulaVar mWShift("bkg_meanWShifted_"+CAT,"@0*@1",RooArgList(mW,*(kMassScale)));
    RooFormulaVar sWShift("bkg_sigmaWShifted_"+CAT,"@0*@1",RooArgList(sW,*(kMassResol)));

    RooGaussian pdfW("bkg_pdfW_"+CAT,"bkg_pdfW_"+CAT,*x,mWShift,sWShift);


    RooRealVar mBkgTop("bkg_meanTop_"+CAT,"bkg_meanTop_"+CAT,172,150,180);
    RooRealVar sBkgTop("bkg_sigmaTop_"+CAT,"bkg_sigmaTop_"+CAT,15,5,30);

    RooFormulaVar mBkgTopShift("bkg_meanTopShifted_"+CAT,"@0*@1",RooArgList(mBkgTop,*(kMassScale)));
    RooFormulaVar sBkgTopShift("bkg_sigmaTopShifted_"+CAT,"@0*@1",RooArgList(sBkgTop,*(kMassResol)));

    RooGaussian bkgTop("bkg_pdfTop_"+CAT,"bkg_pdfTop_"+CAT,*x,mBkgTopShift,sBkgTopShift);

    RooRealVar bBkg0("bkg_b0_"+CAT,"bkg_b0_"+CAT,0.5,0,1);
    RooRealVar bBkg1("bkg_b1_"+CAT,"bkg_b1_"+CAT,0.5,0,1);
    RooRealVar bBkg2("bkg_b2_"+CAT,"bkg_b2_"+CAT,0.5,0,1);
    RooRealVar bBkg3("bkg_b3_"+CAT,"bkg_b3_"+CAT,0.5,0,1);
    RooRealVar bBkg4("bkg_b4_"+CAT,"bkg_b4_"+CAT,0.5,0,1);
    RooRealVar bBkg5("bkg_b5_"+CAT,"bkg_b5_"+CAT,0.5,0,1);
    RooRealVar bBkg6("bkg_b6_"+CAT,"bkg_b6_"+CAT,0.5,0,1);
    RooRealVar bBkg7("bkg_b7_"+CAT,"bkg_b7_"+CAT,0.5,0,1);
    RooRealVar bBkg8("bkg_b8_"+CAT,"bkg_b8_"+CAT,0.5,0,1);

    RooBernstein bkgComb("bkg_pdfComb_"+CAT,"bkg_pdfComb_"+CAT,*x,RooArgList(bBkg0,bBkg1,bBkg2));

    RooRealVar fbkg1("bkg_f1_"+CAT,"bkg_f1_"+CAT,0.9,0.01,1);
    RooRealVar fbkg2("bkg_f2_"+CAT,"bkg_f2_"+CAT,0.1,0.01,1);

    RooAddPdf *bkg = new RooAddPdf("bkg_pdf_"+CAT,"bkg_pdf_"+CAT,RooArgList(bkgTop,pdfW,bkgComb),RooArgList(fbkg1,fbkg2));
    res = bkg->fitTo(*roohBkg,RooFit::Save());
    res->Print();

    TCanvas *canBkg = new TCanvas("Template_Bkg_"+CUT+"_"+CAT,"Template_Bkg_"+CUT+"_"+CAT,900,600);
    RooPlot *frameBkg = x->frame();
    roohBkg->plotOn(frameBkg);
    bkg->plotOn(frameBkg);
    bkg->plotOn(frameBkg,RooFit::Components("bkg_pdfComb_"+CAT),RooFit::LineColor(kRed),RooFit::LineWidth(2),RooFit::LineStyle(2));
    bkg->plotOn(frameBkg,RooFit::Components("bkg_pdfTop_"+CAT),RooFit::LineColor(kGreen+1),RooFit::LineWidth(2),RooFit::LineStyle(2));
    bkg->plotOn(frameBkg,RooFit::Components("bkg_pdfW_"+CAT),RooFit::LineColor(kOrange+1),RooFit::LineWidth(2),RooFit::LineStyle(2));
    frameBkg->GetXaxis()->SetTitle("m_{t} (GeV)");
    frameBkg->Draw();
    gPad->Update();
    canBkg->Print(TString::Format("%s/plots/templateResults/"+TString(canBkg->GetName())+".pdf", year.Data()));

    RooArgSet *parsBkg = (RooArgSet*)bkg->getParameters(roohData);
    parsBkg->setAttribAll("Constant",true);

    w->import(*bkg);
  }

  w->writeToFile(TString::Format("%s/templates_Bkg_"+CUT+"100.root", year.Data()));




}
