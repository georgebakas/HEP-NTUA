void CreateTemplates(TString CUT = "")
{
  gROOT->ForceStyle();
  
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);

  const int NMC = 1;

  TString SAMPLE[NMC] = {
    "TT_Mtt-700toInf_TuneCUETP8M2T4_13TeV-powheg-pythia8_New"
    //"TT_TuneCUETP8M2T4_13TeV-powheg-pythia8",
    //"TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8",  
    //"TT_TuneCUETP8M2T4_13TeV-powheg-fsrdown-pythia8",
    //"TT_TuneCUETP8M2T4_13TeV-powheg-fsrup-pythia8",
    //"TT_TuneCUETP8M2T4_13TeV-powheg-isrdown-pythia8",
    //"TT_TuneCUETP8M2T4_13TeV-powheg-isrdown-pythia8",
    //"TT_TuneCUETP8M2T4down_13TeV-powheg-pythia8",
    //"TT_TuneCUETP8M2T4up_13TeV-powheg-pythia8",
    //"TT_TuneEE5C_13TeV-powheg-herwigpp",
  };
  TFile *infMC[NMC];
  float normMC[NMC];
  for(int k=0;k<NMC;k++) {
    infMC[k] = TFile::Open("Histo_"+SAMPLE[k]+".root");
  }
  
  TFile *infData = TFile::Open("Histo_JetHT_Run2016-17Jul2018_New.root");
  TFile *infBkg  = TFile::Open("Histo_SubdominantBkgs_New.root");
  /*
  TFile *infW  = TFile::Open("Histo_WJetsToQQ_HT180_13TeV-madgraphMLM-pythia8.root");
  TFile *infDY = TFile::Open("Histo_DYJetsToQQ_HT180_13TeV-madgraphMLM-pythia8.root");
  TFile *infST_tW_top = TFile::Open("Histo_ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
  TFile *infST_tW_antitop = TFile::Open("Histo_ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
  TFile *infST_t_top = TFile::Open("Histo_ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root");
  TFile *infST_t_antitop = TFile::Open("Histo_ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root");
  */
  float XSEC(832.);
  float LUMI(37000.);
  
  const int NSRC = 1;
  TString ALIAS[NSRC] = {
    "Nominal"
    /*
    ,"FSRdown","FSRup","ISRdown","ISRup","TUNEdown","TUNEup","HERWIGpp"
    "NoWeight","BtagUp","BtagDown","TrigUp","TrigDown","JER","JERUp","JERDown","JESUp","JESDown",
    "aMC@NLO","Madgraph","ScaleDown","ScaleUp","mtop1665","mtop1695","mtop1715","mtop1735","mtop1755","mtop1785",
    "mpiOFF","herwigpp"
    */
  };
  RooRealVar *kMassScale = new RooRealVar("kMassScale","kMassScale",1.0,0.5,1.5);
  RooRealVar *kMassResol = new RooRealVar("kMassResol","kMassResol",1.0,0.5,1.5);
  kMassScale->setConstant(kTRUE);
  kMassResol->setConstant(kTRUE);

  RooRealVar *YieldTT[NSRC];
  RooRealVar *AccTT[NSRC];
  TH1F *hMC[NSRC];
  TH1F *hData;
  RooDataHist *roohMC[NSRC];
  RooDataHist *roohData,*roohDataCor,*roohDataCorUp,*roohDataCorDown;
  RooAddPdf *qcd,*qcdCor,*qcdCorUp,*qcdCorDown;
  RooAddPdf *signal[NSRC];
  
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
  
  hData = (TH1F*)infData->Get("hWt_mTop_"+CUT+"0btag");
  hData->Rebin(2);
  roohData = new RooDataHist("roohistData","roohistData",RooArgList(*x),hData);  
  TH1F *hDataJet = (TH1F*)infData->Get("hWt_jetMassSoftDrop_"+CUT+"0btag");
  TH1F *hDataEvt = (TH1F*)infData->Get("hWt_mva_"+CUT+"0btag"); 
  RooRealVar *fBkgJet = new RooRealVar("fBkgJet","fBkgJet",hDataJet->Integral()/hData->Integral());
  RooRealVar *fBkgEvt = new RooRealVar("fBkgEvt","fBkgEvt",hDataEvt->Integral()/hData->Integral());
  //---- QCD -----------------------------------
  RooRealVar bQCD0("qcd_b0","qcd_b0",0.5,0,1);
  RooRealVar bQCD1("qcd_b1","qcd_b1",0.5,0,1);
  RooRealVar bQCD2("qcd_b2","qcd_b2",0.5,0,1);
  RooRealVar bQCD3("qcd_b3","qcd_b3",0.5,0,1);
  RooBernstein qcd1("qcd_brn","qcd_brn",*x,RooArgList(bQCD0,bQCD1,bQCD2,bQCD3));

  RooRealVar mQCD("qcd_mean" ,"qcd_mean",140,130,300);
  RooRealVar sQCD("qcd_sigma","qcd_sigma",50,10,200);
  RooGaussian qcd2("qcd_gaus" ,"qcd_gaus",*x,mQCD,sQCD);

  RooRealVar fqcd("qcd_f","qcd_f",0.5,0,1);

  qcd = new RooAddPdf("qcd_pdf","qcd_pdf",qcd1,qcd2,fqcd);
  
  //---- plots ---------------------------------------------------
  TCanvas *canQCD = new TCanvas("Template_QCD_"+CUT,"Template_QCD_"+CUT,900,600);

  RooFitResult *res = qcd->fitTo(*roohData,RooFit::Save());
  //res->Print();
  RooPlot *frameQCD = x->frame();
  roohData->plotOn(frameQCD);
  qcd->plotOn(frameQCD);
  qcd->plotOn(frameQCD,RooFit::Components("qcd_brn"),RooFit::LineColor(kRed),RooFit::LineWidth(2),RooFit::LineStyle(2));
  qcd->plotOn(frameQCD,RooFit::Components("qcd_gaus"),RooFit::LineColor(kGreen+1),RooFit::LineWidth(2),RooFit::LineStyle(2));
  frameQCD->GetXaxis()->SetTitle("m_{t} (GeV)");
  frameQCD->Draw();
  gPad->Update();
  canQCD->Print("plots/"+TString(canQCD->GetName())+".pdf");

  RooArgSet *parsQCD = (RooArgSet*)qcd->getParameters(roohData);
  parsQCD->setAttribAll("Constant",true);

  w->import(*qcd);
  w->import(*fBkgJet);
  w->import(*fBkgEvt);

  RooRealVar mW("meanW","meanW",80,70,90);
  RooRealVar sW("sigmaW","sigmaW",5,0,15);
  RooFormulaVar mWShift("meanWShifted","@0*@1",RooArgList(mW,*(kMassScale)));
  RooFormulaVar sWShift("sigmaWShifted","@0*@1",RooArgList(sW,*(kMassResol)));

  RooGaussian pdfW("pdfW","pdfW",*x,mWShift,sWShift);

  TH1F *hMeanTop[2],*hSigmaTop[2],*hMeanW[2],*hSigmaW[2];

  for(int icat=0;icat<2;icat++) {
    TString CAT = TString::Format("%dbtag",icat+1);
    TAG = CUT+"_"+CAT;
    if (icat==0) {
      mW.setConstant(false);
      sW.setConstant(false);
    }
    else {
      mW.setConstant(false);
      sW.setConstant(false);
    }
    //---- do the bkg templates -------------
    cout<<"hWt_"+VAR+"_"+TAG+"_expYield"<<endl;
    TH1F *hBkg = (TH1F*)infBkg->Get("hWt_"+VAR+TAG+"_expYield"); 
    /*
    TH1F *hW = (TH1F*)infW->Get("boosted/hWt_"+VAR+"_"+TAG);
    TH1F *hDY = (TH1F*)infW->Get("boosted/hWt_"+VAR+"_"+TAG);
    TH1F *hST_tW_top = (TH1F*)infST_tW_top->Get("boosted/hWt_"+VAR+"_"+TAG);
    TH1F *hST_tW_antitop = (TH1F*)infST_tW_antitop->Get("boosted/hWt_"+VAR+"_"+TAG);
    TH1F *hST_t_top = (TH1F*)infST_t_top->Get("boosted/hWt_"+VAR+"_"+TAG);
    TH1F *hST_t_antitop = (TH1F*)infST_t_antitop->Get("boosted/hWt_"+VAR+"_"+TAG);
    float normW = ((TH1F*)infW->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
    float normDY = ((TH1F*)infDY->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
    float normST_tW_top = ((TH1F*)infST_tW_top->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
    float normST_tW_antitop = ((TH1F*)infST_tW_antitop->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
    float normST_t_top = ((TH1F*)infST_t_top->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
    float normST_t_antitop = ((TH1F*)infST_t_antitop->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
    hW->Sumw2();
    hDY->Sumw2();
    hST_tW_top->Sumw2();
    hST_tW_antitop->Sumw2();
    hST_t_top->Sumw2();
    hST_t_antitop->Sumw2();
    hW->Scale(3539*LUMI/normW);
    hDY->Scale(1460.*LUMI/normDY);
    hST_tW_top->Scale(35.6*LUMI/normST_tW_top);
    hST_tW_antitop->Scale(35.6*LUMI/normST_tW_antitop);
    hST_t_top->Scale(136.02*LUMI/normST_t_top);
    hST_t_antitop->Scale(80.95*LUMI/normST_t_antitop);
    TH1F *hBkg = (TH1F*)hST_tW_top->Clone("hBkg");
    //hBkg->Add(hW);
    //hBkg->Add(hDY);
    //hBkg->Add(hST_tW_top);
    hBkg->Add(hST_tW_antitop);
    hBkg->Add(hST_t_top);
    hBkg->Add(hST_t_antitop);
    */
    hBkg->Rebin(2);
    RooDataHist *roohBkg = new RooDataHist("roohistBkg","roohistBkg",RooArgList(*x),hBkg); 
    RooRealVar mBkg1("bkg_mean1_"+CAT,"bkg_mean1_"+CAT,172,150,180);
    RooRealVar sBkg1("bkg_sigma1_"+CAT,"bkg_sigma1_"+CAT,20,5,30);

    RooFormulaVar mBkg1Shift("bkg_mean1Shifted_"+CAT,"@0*@1",RooArgList(mBkg1,*(kMassScale)));
    RooFormulaVar sBkg1Shift("bkg_sigma1Shifted_"+CAT,"@0*@1",RooArgList(sBkg1,*(kMassResol)));

    //RooGaussian bkg1("bkg_gaus1_"+CAT,"bkg_gaus1_"+CAT,*x,mBkg1,sBkg1);
    RooGaussian bkg1("bkg_gaus1_"+CAT,"bkg_gaus1_"+CAT,*x,mBkg1Shift,sBkg1Shift);

    RooRealVar mBkg2("bkg_mean2_"+CAT,"bkg_mean2_"+CAT,80,70,90);
    RooRealVar sBkg2("bkg_sigma2_"+CAT,"bkg_sigma2_"+CAT,5,5,20);

    RooFormulaVar mBkg2Shift("bkg_mean2Shifted_"+CAT,"@0*@1",RooArgList(mBkg2,*(kMassScale)));
    RooFormulaVar sBkg2Shift("bkg_sigma2Shifted_"+CAT,"@0*@1",RooArgList(sBkg2,*(kMassResol)));
    
    RooGaussian bkg2("bkg_gaus2_"+CAT,"bkg_gaus2_"+CAT,*x,mBkg2Shift,sBkg2Shift);
   
    RooRealVar bBkg0("bkg_b0_"+CAT,"bkg_b0_"+CAT,0.5,0,1);
    RooRealVar bBkg1("bkg_b1_"+CAT,"bkg_b1_"+CAT,0.5,0,1);
    RooRealVar bBkg2("bkg_b2_"+CAT,"bkg_b2_"+CAT,0.5,0,1); 
    RooRealVar bBkg3("bkg_b3_"+CAT,"bkg_b3_"+CAT,0.5,0,1);
    RooRealVar bBkg4("bkg_b4_"+CAT,"bkg_b4_"+CAT,0.5,0,1);
    RooRealVar bBkg5("bkg_b5_"+CAT,"bkg_b5_"+CAT,0.5,0,1); 
    RooRealVar bBkg6("bkg_b6_"+CAT,"bkg_b6_"+CAT,0.5,0,1);
    RooRealVar bBkg7("bkg_b7_"+CAT,"bkg_b7_"+CAT,0.5,0,1);
    RooRealVar bBkg8("bkg_b8_"+CAT,"bkg_b8_"+CAT,0.5,0,1);

    RooBernstein bkg3("bkg_bkg_"+CAT,"bkg_bkg_"+CAT,*x,RooArgList(bBkg0,bBkg1,bBkg2)); 

    RooRealVar fbkg1("bkg_f1_"+CAT,"bkg_f1_"+CAT,0.9,0.01,1);
    RooRealVar fbkg2("bkg_f2_"+CAT,"bkg_f2_"+CAT,0.1,0.01,1);

    RooAddPdf *bkg = new RooAddPdf("bkg_pdf_"+CAT,"bkg_pdf_"+CAT,RooArgList(bkg1,bkg2,bkg3),RooArgList(fbkg1,fbkg2));
    res = bkg->fitTo(*roohBkg,RooFit::Save());  
    //res->Print();

    TCanvas *canBkg = new TCanvas("Template_Bkg_"+CUT+"_"+CAT,"Template_Bkg_"+CUT+"_"+CAT,900,600);
    RooPlot *frameBkg = x->frame();
    roohBkg->plotOn(frameBkg);
    bkg->plotOn(frameBkg);
    bkg->plotOn(frameBkg,RooFit::Components("bkg_bkg_"+CAT),RooFit::LineColor(kRed),RooFit::LineWidth(2),RooFit::LineStyle(2));
    bkg->plotOn(frameBkg,RooFit::Components("bkg_gaus1_"+CAT),RooFit::LineColor(kGreen+1),RooFit::LineWidth(2),RooFit::LineStyle(2));
    bkg->plotOn(frameBkg,RooFit::Components("bkg_gaus2_"+CAT),RooFit::LineColor(kOrange+1),RooFit::LineWidth(2),RooFit::LineStyle(2));
    frameBkg->GetXaxis()->SetTitle("m_{t} (GeV)");
    frameBkg->Draw();
    gPad->Update();
    canBkg->Print("plots/"+TString(canBkg->GetName())+".pdf");

    RooArgSet *parsBkg = (RooArgSet*)bkg->getParameters(roohData);
    parsBkg->setAttribAll("Constant",true);

    w->import(*bkg);
  
    //---- then do the signal templates -------------  
    /*
    hMC[1]  = (TH1F*)infMC[0]->Get(TYPE+"/h_"+VAR+"_"+TAG);
    hMC[2]  = (TH1F*)infMC[0]->Get(TYPE+"/hWtBtagUp_"+VAR+"_"+TAG);
    hMC[4]  = (TH1F*)infMC[0]->Get(TYPE+"/hWtTrigUp_"+VAR+"_"+TAG);
    if (TYPE == "resolved") {
      hMC[3]  = (TH1F*)infMC[0]->Get(TYPE+"/hWtBtagDo_"+VAR+"_"+TAG);
      hMC[5]  = (TH1F*)infMC[0]->Get(TYPE+"/hWtTrigDo_"+VAR+"_"+TAG);
    }
    else {
      hMC[3]  = (TH1F*)infMC[0]->Get(TYPE+"/hWtBtagDown_"+VAR+"_"+TAG);
      hMC[5]  = (TH1F*)infMC[0]->Get(TYPE+"/hWtTrigDown_"+VAR+"_"+TAG);
    }
    hMC[6]  = (TH1F*)infMC[0]->Get(TYPE+"Smeared/hWt_"+VAR+"_"+TAG);
    hMC[7]  = (TH1F*)infMC[0]->Get(TYPE+"SmearedUp/hWt_"+VAR+"_"+TAG);
    hMC[8]  = (TH1F*)infMC[0]->Get(TYPE+"SmearedDown/hWt_"+VAR+"_"+TAG);
    hMC[9]  = (TH1F*)infMC[0]->Get(TYPE+"ShiftedUp/hWt_"+VAR+"_"+TAG);
    hMC[10] = (TH1F*)infMC[0]->Get(TYPE+"ShiftedDown/hWt_"+VAR+"_"+TAG);
    hMC[11] = (TH1F*)infMC[1]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
    hMC[12] = (TH1F*)infMC[2]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
    hMC[13] = (TH1F*)infMC[3]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
    hMC[14] = (TH1F*)infMC[4]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
    hMC[15] = (TH1F*)infMC[5]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
    hMC[16] = (TH1F*)infMC[6]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
    hMC[17] = (TH1F*)infMC[7]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
    hMC[18] = (TH1F*)infMC[8]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
    hMC[19] = (TH1F*)infMC[9]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
    hMC[20] = (TH1F*)infMC[10]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
    hMC[21] = (TH1F*)infMC[11]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
    hMC[22] = (TH1F*)infMC[12]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
    */
    TCanvas *canS[NSRC];
    RooRealVar *fSigJet[NSRC],*fSigEvt[NSRC];
    hMeanTop[icat]  = new TH1F("hmTopMean_"+CAT,"hmTopMean_"+CAT,NSRC,0,NSRC);
    hSigmaTop[icat] = new TH1F("hmTopSigma_"+CAT,"hmTopSigma_"+CAT,NSRC,0,NSRC);
    hMeanW[icat]    = new TH1F("hmWMean_"+CAT,"hmWMean_"+CAT,NSRC,0,NSRC);
    hSigmaW[icat]   = new TH1F("hmWSigma_"+CAT,"hmWSigma_"+CAT,NSRC,0,NSRC);
    for(int k=0;k<1;k++) {
      hMC[k]  = (TH1F*)infMC[k]->Get("hWt_"+VAR+TAG+"_expYield"); 
      if (k < 11) {
        //normMC[k] = ((TH1F*)infMC[k]->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
        //cout<<ALIAS[k]<<" "<<normMC[k]<<" "<<((TH1F*)infMC[0]->Get("eventCounter/GenEventWeight"))->GetEntries()<<endl;
      }
      else {
        //normMC[k] = ((TH1F*)infMC[k-10]->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
        //cout<<ALIAS[k]<<" "<<normMC[k]<<" "<<((TH1F*)infMC[k-10]->Get("eventCounter/GenEventWeight"))->GetEntries()<<endl;
      }
      //hMC[k]->Rebin(2);
      cout<<"hWt_jetMassSoftDrop_"+TAG<<endl;
      TH1F *hSig  = (TH1F*)infMC[k]->Get("hWt_mTop"+TAG);
      TH1F *hSigJet = (TH1F*)infMC[k]->Get("hWt_jetMassSoftDrop"+TAG);
      TH1F *hSigEvt = (TH1F*)infMC[k]->Get("hWt_mva"+TAG);
      fSigJet[k] = new RooRealVar("fSigJet_"+ALIAS[k]+"_"+CAT,"fSigJet_"+ALIAS[k]+"_"+CAT,hSigJet->Integral()/hSig->Integral());
      fSigEvt[k] = new RooRealVar("fSigEvt_"+ALIAS[k]+"_"+CAT,"fSigEvt_"+ALIAS[k]+"_"+CAT,hSigEvt->Integral()/hSig->Integral());

      double error(0.0);
      float signal_yield = hMC[k]->IntegralAndError(1,hMC[k]->GetNbinsX(),error);
      float signal_error = error;
    
      YieldTT[k] = new RooRealVar("YieldTT_"+ALIAS[k]+"_"+CAT,"YieldTT_"+ALIAS[k]+"_"+CAT,signal_yield);
      YieldTT[k]->setError(signal_error);

      AccTT[k] = new RooRealVar("AccTT_"+ALIAS[k]+"_"+CAT,"AccTT_"+ALIAS[k]+"_"+CAT,signal_yield/(XSEC*LUMI));
      AccTT[k]->setError(signal_error/(XSEC*LUMI));

      cout<<"Yield "<<ALIAS[k]+"_"+CAT<<": "<<signal_yield<<" +/- "<<signal_error<<endl;
    
      roohMC[k] = new RooDataHist("roohistTT_"+ALIAS[k]+"_"+CAT,"roohistTT_"+ALIAS[k]+"_"+CAT,RooArgList(*x),hMC[k]);    
      RooRealVar m1("ttbar_mean1_"+ALIAS[k]+"_"+CAT,"ttbar_mean1_"+ALIAS[k]+"_"+CAT,172,150,180);
      RooRealVar s1("ttbar_sigma1_"+ALIAS[k]+"_"+CAT,"ttbar_sigma1_"+ALIAS[k]+"_"+CAT,20,5,30);

      RooFormulaVar m1Shift("ttbar_mean1Shifted_"+ALIAS[k]+"_"+CAT,"@0*@1",RooArgList(m1,*(kMassScale)));
      RooFormulaVar s1Shift("ttbar_sigma1Shifted_"+ALIAS[k]+"_"+CAT,"@0*@1",RooArgList(s1,*(kMassResol)));

      RooGaussian sig1("ttbar_gaus1_"+ALIAS[k]+"_"+CAT,"ttbar_gaus2_"+ALIAS[k]+"_"+CAT,*x,m1Shift,s1Shift);

      RooRealVar m2("ttbar_mean2_"+ALIAS[k]+"_"+CAT,"ttbar_mean2_"+ALIAS[k]+"_"+CAT,80,70,90);
      RooRealVar s2("ttbar_sigma2_"+ALIAS[k]+"_"+CAT,"ttbar_sigma2_"+ALIAS[k]+"_"+CAT,5,5,20);

      RooFormulaVar m2Shift("ttbar_mean2Shifted_"+ALIAS[k]+"_"+CAT,"@0*@1",RooArgList(m2,*(kMassScale)));
      RooFormulaVar s2Shift("ttbar_sigma2Shifted_"+ALIAS[k]+"_"+CAT,"@0*@1",RooArgList(s2,*(kMassResol)));

      //RooFormulaVar m2Shift("ttbar_mean2Shifted_"+ALIAS[k]+"_"+CAT,"@0*@1",RooArgList(mW,*(kMassScale)));
      //RooFormulaVar s2Shift("ttbar_sigma2Shifted_"+ALIAS[k]+"_"+CAT,"@0*@1",RooArgList(sW,*(kMassResol)));
   
      RooGaussian sig2("ttbar_gaus2_"+ALIAS[k]+"_"+CAT,"ttbar_gaus2_"+ALIAS[k]+"_"+CAT,*x,m2Shift,s2Shift);

      RooRealVar bSig0("ttbar_b0_"+ALIAS[k]+"_"+CAT,"ttbar_b0_"+ALIAS[k]+"_"+CAT,0.5,0,1);
      RooRealVar bSig1("ttbar_b1_"+ALIAS[k]+"_"+CAT,"ttbar_b1_"+ALIAS[k]+"_"+CAT,0.5,0,1);
      RooRealVar bSig2("ttbar_b2_"+ALIAS[k]+"_"+CAT,"ttbar_b2_"+ALIAS[k]+"_"+CAT,0.5,0,1); 
      RooRealVar bSig3("ttbar_b3_"+ALIAS[k]+"_"+CAT,"ttbar_b3_"+ALIAS[k]+"_"+CAT,0.5,0,1);
      RooRealVar bSig4("ttbar_b4_"+ALIAS[k]+"_"+CAT,"ttbar_b4_"+ALIAS[k]+"_"+CAT,0.5,0,1);
      RooRealVar bSig5("ttbar_b5_"+ALIAS[k]+"_"+CAT,"ttbar_b5_"+ALIAS[k]+"_"+CAT,0.5,0,1); 
      RooRealVar bSig6("ttbar_b6_"+ALIAS[k]+"_"+CAT,"ttbar_b6_"+ALIAS[k]+"_"+CAT,0.5,0,1);
      RooRealVar bSig7("ttbar_b7_"+ALIAS[k]+"_"+CAT,"ttbar_b7_"+ALIAS[k]+"_"+CAT,0.5,0,1);
      RooRealVar bSig8("ttbar_b8_"+ALIAS[k]+"_"+CAT,"ttbar_b8_"+ALIAS[k]+"_"+CAT,0.5,0,1);

      RooBernstein sig3("ttbar_bkg_"+ALIAS[k]+"_"+CAT,"ttbar_bkg_"+ALIAS[k]+"_"+CAT,*x,RooArgList(bSig0,bSig1,bSig2,bSig3,bSig4,bSig5,bSig6,bSig7)); 

      //RooRealVar mL("ttbar_mean3_"+ALIAS[k]+"_"+CAT,"ttbar_mean3_"+ALIAS[k]+"_"+CAT,90,0,200);
      //RooRealVar sL("ttbar_sigma3_"+ALIAS[k]+"_"+CAT,"ttbar_sigma3_"+ALIAS[k]+"_"+CAT,10,0,100);
      //RooLandau sig3("ttbar_bkg_"+ALIAS[k],"ttbar_bkg_"+ALIAS[k],*x,mL,sL);

      RooRealVar fsig1("ttbar_f1_"+ALIAS[k]+"_"+CAT,"ttbar_f1_"+ALIAS[k]+"_"+CAT,0.5,0,1);
      RooRealVar fsig2("ttbar_f2_"+ALIAS[k]+"_"+CAT,"ttbar_f2_"+ALIAS[k]+"_"+CAT,0.1,0.01,1);

      RooAddPdf *signal = new RooAddPdf("ttbar_pdf_"+ALIAS[k]+"_"+CAT,"ttbar_pdf_"+ALIAS[k]+"_"+CAT,RooArgList(sig1,sig2,sig3),RooArgList(fsig1,fsig2));

      canS[k] = new TCanvas("Template_TT_"+CAT+"_"+CUT+"_"+ALIAS[k],"Template_TT_"+CAT+"_"+CUT+"_"+ALIAS[k],900,600);

      RooFitResult *res = signal->fitTo(*roohMC[k],RooFit::Save());
 
      cout<<ALIAS[k]<<endl;
      cout<<"mean1 = "<<m1.getVal()<<" +/- "<<m1.getError()<<", sigma1 = "<<s1.getVal()<<" +/- "<<s1.getError()<<endl;
      cout<<"mean2 = "<<m2.getVal()<<" +/- "<<m2.getError()<<", sigma2 = "<<s2.getVal()<<" +/- "<<s2.getError()<<endl;
      hMeanTop[icat]->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
      hSigmaTop[icat]->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
      hMeanTop[icat]->SetBinContent(k+1,m1.getVal());
      hMeanTop[icat]->SetBinError(k+1,m1.getError());
      hSigmaTop[icat]->SetBinContent(k+1,s1.getVal());
      hSigmaTop[icat]->SetBinError(k+1,s1.getError());

      hMeanW[icat]->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
      hSigmaW[icat]->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
      hMeanW[icat]->SetBinContent(k+1,m2.getVal());
      hMeanW[icat]->SetBinError(k+1,m2.getError());
      hSigmaW[icat]->SetBinContent(k+1,s2.getVal());
      hSigmaW[icat]->SetBinError(k+1,s2.getError());

      //res->Print();
      RooPlot *frameS = x->frame();
      roohMC[k]->plotOn(frameS);
      signal->plotOn(frameS);
      signal->plotOn(frameS,RooFit::Components("ttbar_gaus1_"+ALIAS[k]+"_"+CAT),RooFit::LineColor(kRed),RooFit::LineWidth(2),RooFit::LineStyle(2));
      signal->plotOn(frameS,RooFit::Components("ttbar_gaus2_"+ALIAS[k]+"_"+CAT),RooFit::LineColor(kOrange),RooFit::LineWidth(2),RooFit::LineStyle(2));
      signal->plotOn(frameS,RooFit::Components("ttbar_bkg_"+ALIAS[k]+"_"+CAT),RooFit::LineColor(kGreen+1),RooFit::LineWidth(2),RooFit::LineStyle(2));
      frameS->GetXaxis()->SetTitle("m_{t} (GeV)");
      frameS->Draw();
      gPad->Update();
      canS[k]->Print("plots/"+TString(canS[k]->GetName())+".pdf");

      RooArgSet *parsSig = (RooArgSet*)signal->getParameters(roohMC[k]);
      parsSig->setAttribAll("Constant",true);

      w->import(*signal);
      w->import(*YieldTT[k]);
      w->import(*AccTT[k]);
      w->import(*fSigJet[k]);
      w->import(*fSigEvt[k]);
    }
  }

  TLegend *leg = new TLegend(0.6,0.75,0.9,0.9);
  leg->AddEntry(hMeanTop[0],"1btag","P");
  leg->AddEntry(hMeanTop[1],"2btag","P");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
/*
  TCanvas *can_mean_top = new TCanvas("mTopMean","mTopMean",900,600);
  hMeanTop[1]->SetLineColor(kRed);
  hMeanTop[1]->SetMarkerColor(kRed);
  hMeanTop[1]->SetMarkerStyle(21);
  hMeanTop[0]->GetYaxis()->SetRangeUser(170,180);
  hMeanTop[0]->GetYaxis()->SetTitle("Top mass mean (GeV)");
  hMeanTop[0]->Draw();
  hMeanTop[1]->Draw("same");
  leg->Draw();

  TCanvas *can_mean_W = new TCanvas("mWMean","mWMean",900,600);
  hMeanW[1]->SetLineColor(kRed);
  hMeanW[1]->SetMarkerColor(kRed);
  hMeanW[1]->SetMarkerStyle(21);
  hMeanW[0]->GetYaxis()->SetRangeUser(70,90);
  hMeanW[0]->GetYaxis()->SetTitle("W mass mean (GeV)");
  hMeanW[0]->Draw();
  hMeanW[1]->Draw("same");
  leg->Draw();

  TCanvas *can_sigma_top = new TCanvas("mTopSigma","mTopSigma",900,600);
  hSigmaTop[1]->SetLineColor(kRed);
  hSigmaTop[1]->SetMarkerColor(kRed);
  hSigmaTop[1]->SetMarkerStyle(21);
  hSigmaTop[0]->GetYaxis()->SetRangeUser(10,15);
  hSigmaTop[0]->GetYaxis()->SetTitle("Top mass sigma (GeV)");
  hSigmaTop[0]->Draw();
  hSigmaTop[1]->Draw("same");
  leg->Draw();

  TCanvas *can_sigma_W = new TCanvas("mWSigma","mWSigma",900,600);
  hSigmaW[1]->SetLineColor(kRed);
  hSigmaW[1]->SetMarkerColor(kRed);
  hSigmaW[1]->SetMarkerStyle(21);
  hSigmaW[0]->GetYaxis()->SetRangeUser(0,15);
  hSigmaW[0]->GetYaxis()->SetTitle("W mass sigma (GeV)");
  hSigmaW[0]->Draw();
  hSigmaW[1]->Draw("same");
  leg->Draw();
  //w->Print();
*/
  w->writeToFile("templates_"+CUT+"_workspace.root");
}                            

