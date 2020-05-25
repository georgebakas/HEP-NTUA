void CreateSignalTemplates(TString year, TString CUT = "")
{
  gROOT->ForceStyle();

  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);

  //TFile *infMC = TFile::Open(TString::Format("%s/Histo_TT_Mtt-700toInf_100.root",year.Data()));
  TFile *infMC = TFile::Open(TString::Format("%s/Histo_TT_NominalMC_100.root",year.Data()));
  float normMC;
  float XSEC(832.);
  float LUMI(0);
  if (year.EqualTo("2016")) LUMI = 35920;
  else if (year.EqualTo("2017")) LUMI = 41530;
  else if (year.EqualTo("2018")) LUMI = 59740;


  std::vector<float> XSECMTT;
  XSECMTT.push_back(69.64);
  XSECMTT.push_back(16.74);

  RooRealVar *kMassScale = new RooRealVar("kMassScale","kMassScale",1.0,0.5,1.5);
  RooRealVar *kMassResol = new RooRealVar("kMassResol","kMassResol",1.0,0.5,1.5);
  kMassScale->setConstant(kTRUE);
  kMassResol->setConstant(kTRUE);

  RooRealVar *YieldTT;
  RooRealVar *AccTT;
  TH1F *hMC;
  RooDataHist *roohMC;
  RooAddPdf *signal;

  TString VAR,TAG;
  float XMIN,XMAX;

  VAR = "mTop";
  XMIN = 50.;
  XMAX = 300.;

  RooWorkspace *w = new RooWorkspace("w","workspace");

  //---- define observable ------------------------
  RooRealVar *x = new RooRealVar("mTop","mTop",XMIN,XMAX);
  w->import(*x);

  RooRealVar mW("meanW","meanW",80,70,90);
  RooRealVar sW("sigmaW","sigmaW",5,0,15);
  RooFormulaVar mWShift("meanWShifted","@0*@1",RooArgList(mW,*(kMassScale)));
  RooFormulaVar sWShift("sigmaWShifted","@0*@1",RooArgList(sW,*(kMassResol)));

  RooGaussian pdfW("pdfW","pdfW",*x,mWShift,sWShift);

  for(int icat=0;icat<3;icat++) {
    TString CAT = TString::Format("%dbtag",icat);
    TAG = CUT+"_"+CAT;
    if (icat==0) {
      mW.setConstant(false);
      sW.setConstant(false);
    }
    else {
      mW.setConstant(false);
      sW.setConstant(false);
    }
    cout<<TAG<<endl;
    //---- then do the signal templates -------------
    hMC = (TH1F*)infMC->Get("hWt_"+VAR+TAG);
    TH1F *hMC_yield = (TH1F*)infMC->Get("hWt_"+VAR+TAG+"_expYield");

    TCanvas *canS;
    RooRealVar *fSigJet,*fSigEvt;
    TH1F *hSig  = (TH1F*)infMC->Get("hWt_mTop"+TAG+"_expYield");
    TH1F *hSigJet = (TH1F*)infMC->Get("hWt_jetMassSoftDrop"+TAG+"_expYield");
    TH1F *hSigEvt = (TH1F*)infMC->Get("hWt_mJJ"+TAG+"_expYield");
    fSigJet = new RooRealVar("fSigJet_"+CAT,"fSigJet_"+CAT,hSigJet->Integral()/hSig->Integral());
    fSigEvt = new RooRealVar("fSigEvt_"+CAT,"fSigEvt_"+CAT,hSigEvt->Integral()/hSig->Integral());

    double error(0.0);
    float signal_yield, signal_error;
    signal_yield =  hMC_yield->IntegralAndError(1,hMC_yield->GetNbinsX(),error);
    signal_error =  error;

    YieldTT = new RooRealVar("YieldTT_"+CAT,"YieldTT_"+CAT,signal_yield);
    YieldTT->setError(signal_error);

    AccTT = new RooRealVar("AccTT_"+CAT,"AccTT_"+CAT,signal_yield/(XSEC*LUMI));
    AccTT->setError(signal_error/(XSEC*LUMI));

    cout<<"Yield "<<CAT<<": "<<signal_yield<<" +/- "<<signal_error<<endl;

    roohMC = new RooDataHist("roohistTT_"+CAT,"roohistTT_"+CAT,RooArgList(*x),hMC_yield);
    RooRealVar mTop("ttbar_meanTop_"+CAT,"ttbar_meanTop_"+CAT,172,150,180);
    RooRealVar sTop("ttbar_sigmaTop_"+CAT,"ttbar_sigmaTop_"+CAT,20,5,30);

    RooFormulaVar mTopShift("ttbar_meanTopShifted_"+CAT,"@0*@1",RooArgList(mTop,*(kMassScale)));
    RooFormulaVar sTopShift("ttbar_sigmaTopShifted_"+CAT,"@0*@1",RooArgList(sTop,*(kMassResol)));

    RooGaussian sigTop("ttbar_pdfTop_"+CAT,"ttbar_pdfTop_"+CAT,*x,mTopShift,sTopShift);

    RooRealVar mW("ttbar_meanW_"+CAT,"ttbar_meanW_"+CAT,90,70,100);
    RooRealVar sW("ttbar_sigmaW_"+CAT,"ttbar_sigmaW_"+CAT,5,5,10);

    //RooFormulaVar mWShift("ttbar_meanWShifted_"+CAT,"@0*@1",RooArgList(mW,*(kMassScale)));
    //RooFormulaVar sWShift("ttbar_sigmaWShifted_"+CAT,"@0*@1",RooArgList(sW,*(kMassResol)));

    RooGaussian sigW("ttbar_pdfW_"+CAT,"ttbar_pdfW_"+CAT,*x,mW,sW   );

    RooRealVar bSig0("ttbar_b0_"+CAT,"ttbar_b0_"+CAT,0.5,0,1);
    RooRealVar bSig1("ttbar_b1_"+CAT,"ttbar_b1_"+CAT,0.5,0,1);
    RooRealVar bSig2("ttbar_b2_"+CAT,"ttbar_b2_"+CAT,0.5,0,1);
    RooRealVar bSig3("ttbar_b3_"+CAT,"ttbar_b3_"+CAT,0.5,0,1);
    RooRealVar bSig4("ttbar_b4_"+CAT,"ttbar_b4_"+CAT,0.5,0,1);
    RooRealVar bSig5("ttbar_b5_"+CAT,"ttbar_b5_"+CAT,0.5,0,1);
    RooRealVar bSig6("ttbar_b6_"+CAT,"ttbar_b6_"+CAT,0.5,0,1);
    RooRealVar bSig7("ttbar_b7_"+CAT,"ttbar_b7_"+CAT,0.5,0,1);
    RooRealVar bSig8("ttbar_b8_"+CAT,"ttbar_b8_"+CAT,0.5,0,1);

    RooBernstein sigComb("ttbar_bkg_"+CAT,"ttbar_bkg_"+CAT,*x,RooArgList(bSig0,bSig1,bSig2,bSig3));

    RooRealVar m1("m_"+CAT, "m1"+CAT, 130,120,200);
    RooRealVar s1("s_"+CAT, "s1"+CAT, 50, 10, 100);
    RooGaussian g1("gaus_"+CAT, "g1"+CAT, *x, m1,s1);

    RooRealVar fsig1("ttbar_f1_"+CAT,"ttbar_f1_"+CAT,0,0,1);
    RooRealVar fsig2("ttbar_f2_"+CAT,"ttbar_f2_"+CAT,0.1,0.01,1);
    RooRealVar fsig3("ttbar_f3_"+CAT,"ttbar_f3_"+CAT,0.1,0.01,1);

    RooAddPdf *signal = new RooAddPdf("ttbar_pdf_"+CAT,"ttbar_pdf_"+CAT,RooArgList(sigTop,sigW,sigComb,g1),RooArgList(fsig1,fsig2,fsig3));
    cout<<"ttbar_pdf_"+CAT<<endl;

    canS = new TCanvas("Template_TT_"+CAT+"_"+CUT,"Template_TT_"+CAT+"_"+CUT,900,600);

    RooFitResult *res = signal->fitTo(*roohMC,RooFit::Save());

    res->Print();
    RooPlot *frameS = x->frame();
    roohMC->plotOn(frameS);
    signal->plotOn(frameS);
    signal->plotOn(frameS,RooFit::Components("ttbar_pdfTop_"+CAT),RooFit::LineColor(kRed),RooFit::LineWidth(2),RooFit::LineStyle(2));
    signal->plotOn(frameS,RooFit::Components("ttbar_pdfW_"+CAT),RooFit::LineColor(kOrange),RooFit::LineWidth(2),RooFit::LineStyle(2));
    signal->plotOn(frameS,RooFit::Components("ttbar_bkg_"+CAT),RooFit::LineColor(kGreen+1),RooFit::LineWidth(2),RooFit::LineStyle(2));
    signal->plotOn(frameS,RooFit::Components("gaus_"+CAT),RooFit::LineColor(kMagenta+1),RooFit::LineWidth(2),RooFit::LineStyle(2));
    frameS->GetXaxis()->SetTitle("m_{t} (GeV)");
    frameS->Draw();
    gPad->Update();
    canS->Print(TString::Format("%s/plots/templateResults/"+TString(canS->GetName())+".pdf",year.Data()));

    RooArgSet *parsSig = (RooArgSet*)signal->getParameters(roohMC);
    parsSig->setAttribAll("Constant",true);

    w->import(*signal);
    w->import(*YieldTT);
    w->import(*AccTT);
    w->import(*fSigJet);
    w->import(*fSigEvt);
  }

  w->writeToFile(TString::Format("%s/templates_Sig_"+CUT+"100.root",year.Data()));
}
