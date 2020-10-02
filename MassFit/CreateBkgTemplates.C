  void CreateBkgTemplates(TString year, int icat, TString CUT = "")
{
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

  TFile *infBkg = TFile::Open(TString::Format("%s/Histo_SubdominantBkgs_100.root",year.Data()));

  float LUMI(0);
  if (year.EqualTo("2016")) LUMI = 35920;
  else if (year.EqualTo("2017")) LUMI = 41530;
  else if (year.EqualTo("2018")) LUMI = 59740;
  cout<<icat<<endl;
  if (icat == 1)
    return;
    //infBkg = TFile::Open(TString::Format("%s/Histo_SubdominantBkgs_100.root", year.Data()));

    TString CAT = TString::Format("%dbtag", icat);
    TAG = CUT + "_" + CAT;
    TH1F *hBkg;
    //---- do the bkg templates -------------
    if (icat == 0)
    {
      hBkg = (TH1F *)infBkg->Get("hWt_mTop_0btag_expYield");
    }
    else if (icat == 2)
    {
      hBkg = (TH1F *)infBkg ->Get("hWt_mTop_2btag_expYield");
    }
    hBkg->Rebin(2);

    RooDataHist *roohBkg = new RooDataHist("roohistBkg", "roohistBkg", RooArgList(*x), hBkg);

    RooRealVar mW("bkg_meanW_" + CAT, "meanW_" + CAT, 80, 70, 90);
    RooRealVar sW("bkg_sigmaW_" + CAT, "sigmaW_" + CAT, 5, 0, 15);
    RooGaussian pdfW("bkg_pdfW_" + CAT, "bkg_pdfW_" + CAT, *x, mW, sW);

    RooRealVar mBkgTop("bkg_meanTop_" + CAT, "bkg_meanTop_" + CAT, 172, 150, 180);
    RooRealVar sBkgTop("bkg_sigmaTop_" + CAT, "bkg_sigmaTop_" + CAT, 15, 5, 30);
    RooGaussian bkgTop("bkg_pdfTop_" + CAT, "bkg_pdfTop_" + CAT, *x, mBkgTop, sBkgTop);

    RooRealVar bBkg0("bkg_b0_" + CAT, "bkg_b0_" + CAT, 0.5, 0, 1);
    RooRealVar bBkg1("bkg_b1_" + CAT, "bkg_b1_" + CAT, 0.5, 0, 1);
    RooRealVar bBkg2("bkg_b2_" + CAT, "bkg_b2_" + CAT, 0.5, 0, 1);
    RooRealVar bBkg3("bkg_b3_" + CAT, "bkg_b3_" + CAT, 0.5, 0, 1);
    RooRealVar bBkg4("bkg_b4_" + CAT, "bkg_b4_" + CAT, 0.5, 0, 1);
    RooRealVar bBkg5("bkg_b5_" + CAT, "bkg_b5_" + CAT, 0.5, 0, 1);
    RooRealVar bBkg6("bkg_b6_" + CAT, "bkg_b6_" + CAT, 0.5, 0, 1);
    RooRealVar bBkg7("bkg_b7_" + CAT, "bkg_b7_" + CAT, 0.5, 0, 1);
    RooRealVar bBkg8("bkg_b8_" + CAT, "bkg_b8_" + CAT, 0.5, 0, 1);

    RooBernstein bkgComb("bkg_pdfComb_" + CAT, "bkg_pdfComb_" + CAT, *x, RooArgList(bBkg0, bBkg1, bBkg2));

    RooRealVar fbkg1("bkg_f1_" + CAT, "bkg_f1_" + CAT, 0.9, 0.01, 1);
    RooRealVar fbkg2("bkg_f2_" + CAT, "bkg_f2_" + CAT, 0.1, 0.01, 1);

    RooAddPdf *bkg = new RooAddPdf("bkg_pdf_" + CAT, "bkg_pdf_" + CAT, RooArgList(bkgTop, pdfW, bkgComb), RooArgList(fbkg1, fbkg2));
    RooFitResult *res = bkg->fitTo(*roohBkg, RooFit::Save());
    res->Print();

    TCanvas *canBkg = new TCanvas("Template_Bkg_" + CUT + "_" + CAT, "Template_Bkg_" + CUT + "_" + CAT, 900, 600);
    RooPlot *frameBkg = x->frame();
    roohBkg->plotOn(frameBkg);
    bkg->plotOn(frameBkg);
    bkg->plotOn(frameBkg, RooFit::Components("bkg_pdfComb_" + CAT), RooFit::LineColor(kRed), RooFit::LineWidth(2), RooFit::LineStyle(2));
    bkg->plotOn(frameBkg, RooFit::Components("bkg_pdfTop_" + CAT), RooFit::LineColor(kGreen + 1), RooFit::LineWidth(2), RooFit::LineStyle(2));
    bkg->plotOn(frameBkg, RooFit::Components("bkg_pdfW_" + CAT), RooFit::LineColor(kOrange + 1), RooFit::LineWidth(2), RooFit::LineStyle(2));
    frameBkg->GetXaxis()->SetTitle("m_{t} (GeV)");
    frameBkg->Draw();
    gPad->Update();
    canBkg->Print(TString::Format("%s/plots/templateResults/" + TString(canBkg->GetName()) + ".pdf", year.Data()));

    RooArgSet *parsBkg = (RooArgSet*)bkg->getParameters(roohBkg);
    parsBkg->setAttribAll("Constant", true);

    w->import(*bkg);


  w->writeToFile(TString::Format("%s/templates_Bkg_"+CUT+"100_"+icat+".root", year.Data()));




}
