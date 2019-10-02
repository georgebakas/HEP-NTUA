void CreateBkgTemplates(TString CUT)
{
  gROOT->ForceStyle();
  
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);

  TFile *infData = TFile::Open("../Histo_JetHT.root");
  TFile *infW  = TFile::Open("../Histo_WJetsToQQ_HT180_13TeV-madgraphMLM-pythia8.root");
  TFile *infDY = TFile::Open("../Histo_DYJetsToQQ_HT180_13TeV-madgraphMLM-pythia8.root");
  TFile *infST_tW_top = TFile::Open("../Histo_ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
  TFile *infST_tW_antitop = TFile::Open("../Histo_ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
  TFile *infST_t_top = TFile::Open("../Histo_ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root");
  TFile *infST_t_antitop = TFile::Open("../Histo_ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root");

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
  
  TH1F *hData = (TH1F*)infData->Get("boosted/h_mTop_"+CUT+"_0btag");
  hData->Rebin(2);
  RooDataHist *roohData = new RooDataHist("roohistData","roohistData",RooArgList(*x),hData);  
  TH1F *hDataJet = (TH1F*)infData->Get("boosted/h_jetMassSoftDrop_"+CUT+"_0btag");
  TH1F *hDataEvt = (TH1F*)infData->Get("boosted/h_mva_"+CUT+"_0btag"); 
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

  RooAddPdf *qcd = new RooAddPdf("qcd_pdf","qcd_pdf",qcd1,qcd2,fqcd);
  
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

  //TH1F *hMeanTop[2],*hSigmaTop[2],*hMeanW[2],*hSigmaW[2];

  float LUMI(37000);

  for(int icat=0;icat<2;icat++) {
    TString CAT = TString::Format("%dbtag",icat+1);
    TAG = CUT+"_"+CAT;
    
    //---- do the bkg templates -------------
    TH1F *hW = (TH1F*)infW->Get("boosted/hWt_"+VAR+"_"+TAG);
    TH1F *hDY = (TH1F*)infDY->Get("boosted/hWt_"+VAR+"_"+TAG);
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
    hW->Scale(3539.*LUMI/normW);
    hDY->Scale(1460.*LUMI/normDY);
    hST_tW_top->Scale(35.6*LUMI/normST_tW_top);
    hST_tW_antitop->Scale(35.6*LUMI/normST_tW_antitop);
    hST_t_top->Scale(136.02*LUMI/normST_t_top);
    hST_t_antitop->Scale(80.95*LUMI/normST_t_antitop);
    TH1F *hBkg = (TH1F*)hW->Clone("hBkg");
    //hBkg->Add(hW);
    hBkg->Add(hDY);
    hBkg->Add(hST_tW_top);
    hBkg->Add(hST_tW_antitop);
    hBkg->Add(hST_t_top);
    hBkg->Add(hST_t_antitop);

    cout<<"Expected number of events: "<<endl;
    cout<<"WJets:         "<<hW->Integral()<<endl; 
    cout<<"DYJets:        "<<hDY->Integral()<<endl;
    cout<<"ST_tW_top:     "<<hST_tW_top->Integral()<<endl;
    cout<<"ST_tW_antitop: "<<hST_tW_antitop->Integral()<<endl;
    cout<<"ST_t_top:      "<<hST_t_top->Integral()<<endl;
    cout<<"ST_t_antitop:  "<<hST_t_antitop->Integral()<<endl;
  
    hBkg->Rebin(5);
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
    //res->Print();

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
    canBkg->Print("plots/"+TString(canBkg->GetName())+".pdf");

    RooArgSet *parsBkg = (RooArgSet*)bkg->getParameters(roohData);
    parsBkg->setAttribAll("Constant",true);

    w->import(*bkg);
  }  

  w->writeToFile("templates_Bkg_"+CUT+".root");
}                            

