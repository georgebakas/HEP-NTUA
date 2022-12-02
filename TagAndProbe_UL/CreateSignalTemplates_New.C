TH1F *addHistograms(TH1F *h1, TH1F *h2, TH1F *h3)
{
  TH1F *hComb = (TH1F*)h1->Clone(h1->GetName());
  hComb->Add(h2);
  hComb->Add(h3);
  return hComb;
}

void CreateSignalTemplates_New(TString year, TString dir, TString inputFile, TString selection="probe", TString CUT = "")
{
  TString selectedRegion;
  if(selection.EqualTo("probe")) selectedRegion = "hSRBTightAndProbe";
  else if (selection.EqualTo("SR")) selectedRegion = "hSRBTightAndSR"; 
  gROOT->ForceStyle();
  gROOT->ForceStyle();
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);
  TString VAR,TAG;
  float XMIN,XMAX;
  XMIN = 120;
  XMAX = 220;
  TString tempFileName;
  if(dir.EqualTo("PDFWeights")) tempFileName = "_pdf_"+inputFile;
  else if (dir.EqualTo("ScaleWeights")) tempFileName = "_scale_"+inputFile;
  else if (dir.EqualTo("Nominal")) tempFileName = "Nominal";
  else if (dir.EqualTo("bTagVariation")) tempFileName = "_bTagVariation_"+inputFile;
  else if (dir.EqualTo("PSWeights")) tempFileName = "_PSWeights_"+inputFile;
  else if (dir.EqualTo("JES")) tempFileName = "_JES_"+inputFile;
  else tempFileName = "_"+inputFile;


  TFile *inf_had, *inf_sem, *inf_dil;
  cout<< TString::Format("%s/%s/combined/TagAndProbeHisto_1000_TT_%s.root", year.Data(), dir.Data(), tempFileName.Data()) << endl;
  inf_had = TFile::Open(TString::Format("%s/%s/combined/TagAndProbeHisto_1000_TT_%s.root", year.Data(), dir.Data(), tempFileName.Data())); //nominal Hadronic
  //inf_sem = TFile::Open(TString::Format("%s/%s/Histo_1000_TTToSemiLeptonic%s.root", year.Data(), dir.Data(), tempFileName.Data())); //nominal SemiLeptonic
  //inf_dil = TFile::Open(TString::Format("%s/%s/Histo_1000_TTTo2L2Nu%s.root", year.Data(), dir.Data(), tempFileName.Data())); //nominal 2L2Nu

  //histograms to be added
  TH1F *h_had, *h_sem, *h_dil;

  RooRealVar *kMassScale = new RooRealVar("kMassScale","kMassScale",1.0,0.5,1.5);
  RooRealVar *kMassResol = new RooRealVar("kMassResol","kMassResol",1.0,0.5,1.5);
  kMassScale->setConstant(kTRUE);
  kMassResol->setConstant(kTRUE);
  
  RooRealVar *YieldTT;
  RooDataHist *roohMC;



  RooWorkspace *w = new RooWorkspace("w","workspace");

  //---- define observable ------------------------
  RooRealVar *x = new RooRealVar("mTop","mTop",XMIN,XMAX);
  w->import(*x);


  for(int icat=2;icat<3;icat++) 
  {
    if(icat ==1) continue;
    //TH1F *hMC_yield = (TH1F*)inf_had->Get(TString::Format("hWt_mTop_Leading_%dbtag", icat));
    TH1F *hMC_yield = (TH1F*)inf_had->Get(TString::Format("%s_mTop_Leading_expYield",selectedRegion.Data()));

    //TH1F *hMC_yield = addHistograms(h_had, h_sem, h_dil);
    TString CAT = TString::Format("%dbtag",icat);
    TAG = CUT+"_"+CAT;
    cout<<TAG<<endl;
    //---- then do the signal templates -------------
    hMC_yield->Rebin(2);
    double error(0.0);
    float signal_yield, signal_error;
    signal_yield =  hMC_yield->IntegralAndError(1,hMC_yield->GetNbinsX(),error);
    signal_error =  error;
    cout<<"Yield "<<CAT<<": "<<signal_yield<<" +/- "<<signal_error<<endl;
    cout<<"Entries: "<<hMC_yield->GetEntries()<<endl;

    TCanvas *canS;

    YieldTT = new RooRealVar("YieldTT_"+CAT,"YieldTT_"+CAT,signal_yield);
    YieldTT->setError(signal_error);

    roohMC = new RooDataHist("roohistTT_"+CAT,"roohistTT_"+CAT,RooArgList(*x),hMC_yield);
    RooRealVar mTop("ttbar_meanTop_"+CAT,"ttbar_meanTop_"+CAT,172,150,180);
    RooRealVar sTop("ttbar_sigmaTop_"+CAT,"ttbar_sigmaTop_"+CAT,20,5,30);

    RooFormulaVar mTopShift("ttbar_meanTopShifted_" + CAT, "@0*@1", RooArgList(mTop, (*kMassScale)));
    RooFormulaVar sTopShift("ttbar_sigmaTopShifted_" + CAT, "@0*@1", RooArgList(sTop, (*kMassResol)));

    RooGaussian sigTop("ttbar_pdfTop_" + CAT, "ttbar_pdfTop_" + CAT, *x, mTopShift, sTopShift);

    RooRealVar mW("ttbar_meanW_"+CAT,"ttbar_meanW_"+CAT,90,70,100);
    RooRealVar sW("ttbar_sigmaW_"+CAT,"ttbar_sigmaW_"+CAT,5,5,10);
    mW.setConstant(false);
    sW.setConstant(false);

    RooFormulaVar mWShift("ttbar_meanWShifted_"+CAT,"@0*@1",RooArgList(mW,*(kMassScale)));
    RooFormulaVar sWShift("ttbar_sigmaWShifted_"+CAT,"@0*@1",RooArgList(sW,*(kMassResol)));

    RooGaussian sigW("ttbar_pdfW_"+CAT,"ttbar_pdfW_"+CAT,*x,mWShift,sWShift);

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
    RooRealVar fsig3("ttbar_f3_"+CAT,"ttbar_f3_"+CAT,0.1,0.,1);

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
    canS->Print(TString::Format("%s/%s/plots/templateResults/"+TString(canS->GetName())+".pdf",year.Data(), dir.Data()));

    RooArgSet *parsSig = (RooArgSet*)signal->getParameters(roohMC);
    parsSig->setAttribAll("Constant",true);

    w->import(*signal);
    w->import(*YieldTT);
  }

  //cout<<TString::Format("%s/%s/SignalTemplates_%s.root",year.Data(), dir.Data(), inputFile.Data())<<endl;
  w->writeToFile(TString::Format("%s/%s/SignalTemplates_%s.root",year.Data(), dir.Data(), inputFile.Data()));
}
