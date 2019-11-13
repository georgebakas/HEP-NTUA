void CreateQCDTemplates_MC(TString CUT = "")
{
  gROOT->ForceStyle();
  
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);

  TFile *infData = TFile::Open("Histo_QCD_HT300toInf_TuneCP5_13TeV-madgraph-pythia8_New_100.root"); 

  TFile *infBkg = TFile::Open("Histo_SubdominantBkgs_New_100.root");
  RooRealVar *kMassScale = new RooRealVar("kMassScale","kMassScale",1.0,0.5,1.5);
  RooRealVar *kMassResol = new RooRealVar("kMassResol","kMassResol",1.0,0.5,1.5);
  kMassScale->setConstant(kTRUE);
  kMassResol->setConstant(kTRUE);
  
  TString VAR,TAG;
  float XMIN,XMAX;

  VAR = "mTop";
  XMIN = 0.;
  XMAX = 300.; 

  RooWorkspace *w = new RooWorkspace("w","workspace");

  //---- define observable ------------------------
  RooRealVar *x = new RooRealVar("mTop","mTop",XMIN,XMAX);
  w->import(*x);
  //---- first do the data template ---------------
  
  TH1F *hData = (TH1F*)infData->Get("hWt_mTop_0btag_expYield");
  hData->Rebin(2);
  RooDataHist *roohData = new RooDataHist("roohistData","roohistData",RooArgList(*x),hData);  
  TH1F *hDataJet = (TH1F*)infData->Get("hWt_jetMassSoftDrop_0btag");
  TH1F *hDataEvt = (TH1F*)infData->Get("hWt_mJJ_0btag"); 
  RooRealVar *fBkgJet = new RooRealVar("fBkgJet","fBkgJet",hDataJet->Integral()/hData->Integral());
  RooRealVar *fBkgEvt = new RooRealVar("fBkgEvt","fBkgEvt",hDataEvt->Integral()/hData->Integral());
  //---- QCD -----------------------------------
  RooRealVar bQCD0("qcd_b0","qcd_b0",0.5,0,1.6);
  RooRealVar bQCD1("qcd_b1","qcd_b1",0.5,0,1.6);
  RooRealVar bQCD2("qcd_b2","qcd_b2",0.5,0,1.6);
  RooRealVar bQCD3("qcd_b3","qcd_b3",0.5,0,1.6);
  RooRealVar bQCD4("qcd_b4","qcd_b4",0.5,0,1.6);
  RooRealVar bQCD5("qcd_b5","qcd_b5",0.5,0,1.6);
  RooRealVar bQCD6("qcd_b6","qcd_b6",0.5,0,1.6);
  RooRealVar bQCD7("qcd_b7","qcd_b7",0.5,0,1.6);
  RooRealVar bQCD8("qcd_b8","qcd_b8",0.5,0,1.6);
  RooBernstein qcd1("qcd_brn","qcd_brn",*x,RooArgList(bQCD0,bQCD1,bQCD2,bQCD3, bQCD4));

  RooRealVar mQCD("qcd_mean" ,"qcd_mean",150,130,300);
  RooRealVar sQCD("qcd_sigma","qcd_sigma",50,10,200);
  RooGaussian qcd2("qcd_gaus" ,"qcd_gaus",*x,mQCD,sQCD);

  RooRealVar fqcd1("qcd_f1","qcd_f1",0.5,0,1);

  RooAddPdf *qcd = new RooAddPdf("qcd_pdf","qcd_pdf",RooArgList(qcd1,qcd2), RooArgList(fqcd1));
  
  //---- plots ---------------------------------------------------
  TCanvas *canQCD = new TCanvas("Template_QCD_"+CUT,"Template_QCD_"+CUT,900,600);

  RooFitResult *res = qcd->fitTo(*roohData,RooFit::Save()); 
  res->Print();
  RooPlot *frameQCD = x->frame();
  roohData->plotOn(frameQCD);
  qcd->plotOn(frameQCD);
  qcd->plotOn(frameQCD,RooFit::Components("qcd_brn"),RooFit::LineColor(kRed),RooFit::LineWidth(2),RooFit::LineStyle(2));
  qcd->plotOn(frameQCD,RooFit::Components("qcd_gaussian2"),RooFit::LineColor(kRed),RooFit::LineWidth(2),RooFit::LineStyle(2));
  qcd->plotOn(frameQCD,RooFit::Components("qcd_gaus"),RooFit::LineColor(kGreen+1),RooFit::LineWidth(2),RooFit::LineStyle(2));
  frameQCD->GetXaxis()->SetTitle("m_{t} (GeV)");
  frameQCD->Draw();
  gPad->Update();
  canQCD->Print("plots/"+TString(canQCD->GetName())+".pdf");

  //RooArgSet *parsQCD = (RooArgSet*)qcd->getParameters(roohData);
  //parsQCD->setAttribAll("Constant",true);

  w->import(*qcd);
  w->import(*fBkgJet);
  w->import(*fBkgEvt);
}