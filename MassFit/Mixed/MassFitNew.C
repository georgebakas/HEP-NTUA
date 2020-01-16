//#include "../../CMS_lumi.C"

using namespace RooFit;
void MassFitNew(TString year = "2016", TString ALIAS="",TString CUT="", int REBIN= 5)
{
  gROOT->ForceStyle();
  
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);
  
  //Take SR data from Medium WP
  TFile *inf = TFile::Open(TString::Format("%s/Histo_Data_%s_100.root", year.Data(), year.Data()));
  TH1F *h2b  = (TH1F*)inf->Get("hWt_mTop"+CUT+"_2btag");
  //h2b->Rebin(REBIN);
  // -----------------------------------------
  const float LUMI = 35922;
  
  TFile *fTemplatesBkg = TFile::Open(TString::Format("%s/templates_Bkg_"+CUT+"100.root", year.Data()));
  TFile *fTemplatesSig = TFile::Open(TString::Format("%s/templates_Sig_"+CUT+"100.root",year.Data()));
  RooWorkspace *wTemplatesBkg = (RooWorkspace*)fTemplatesBkg->Get("w");
  RooWorkspace *wTemplatesSig = (RooWorkspace*)fTemplatesSig->Get("w");
  
  RooRealVar *x = (RooRealVar*)wTemplatesSig->var("mTop");
  RooRealVar *yieldTT = (RooRealVar*)wTemplatesSig->var("YieldTT_2btag");
  
  RooRealVar *kMassScale = (RooRealVar*)wTemplatesSig->var("kMassScale");
  RooRealVar *kMassResol = (RooRealVar*)wTemplatesSig->var("kMassResol");
  kMassScale->setConstant(false);
  kMassResol->setConstant(false);
    
  RooDataHist *roohist_data_2b = new RooDataHist("roohist_data_2b","roohist_data_2b",*x,h2b);

  RooCategory sample("sample","sample");
  sample.defineType("2btag");
  
  RooDataHist combData("combData","combData",*x,Index(sample),Import("2btag",*h2b));
    
  
  RooAbsPdf *pdf_bkg_2b = (RooAbsPdf*)wTemplatesBkg->pdf("bkg_pdf_2btag");
  
  RooAbsPdf *pdf_qcd_2b = (RooAbsPdf*)wTemplatesBkg->pdf("qcd_pdf");
  
  //---- QCD correction factor ---------------------------
  
  RooRealVar kQCD2b("kQCD_2b","kQCD_2b",1e-3,-1,1);
  
  kQCD2b.setConstant(false);
  
  RooFormulaVar qcdCor_2b("qcdCor_2b","1+@0*@1",RooArgList(*x,kQCD2b));
  //---- corrected QCD -----------------------------------
  
  RooEffProd pdf_qcdCor_2b("qcdCor_pdf_2b","qcdCor_pdf_2b",*pdf_qcd_2b,qcdCor_2b);
  
  RooRealVar *nFitBkg2b = new RooRealVar("nFitBkg_2b","nFitBkg_2b",400,0,1e+4);
 
  RooRealVar *nFitQCD2b = new RooRealVar("nFitQCD_2b","nFitQCD_2b",10000,0,1e+5);  

  
  
  RooRealVar *nFitSig2b = new RooRealVar("nFitSig2b","nFitSig2b",2000,100,1e+5);
  
  RooAbsPdf *pdf_signal_2b = (RooAbsPdf*)wTemplatesSig->pdf("ttbar_pdf_2btag"); 

  RooAddPdf *model_2b = new RooAddPdf("model_2b","model_2b",RooArgList(*pdf_signal_2b,pdf_qcdCor_2b,*pdf_bkg_2b),RooArgList(*nFitSig2b,*nFitQCD2b,*nFitBkg2b));


  RooSimultaneous simPdf("simPdf","simPdf",sample);
  simPdf.addPdf(*model_2b,"2btag");
    
  RooFitResult *res = simPdf.fitTo(combData,RooFit::Save(),RooFit::Extended(kTRUE));

  cout<<"Signal strength: r = "<<nFitSig2b->getVal()/yieldTT->getVal()<<endl;

  res->Print();

  //cout<<"correlation = "<<res->correlation(*nFitSig,*btagEff)<<endl;
  //cout<<"correlation = "<<res->correlation(*nFitQCD2b,*btagEff)<<endl;
  /*
  RooAbsReal *nll = simPdf.createNLL(combData,NumCPU(2));
  RooMinuit(*nll).migrad();
  RooAbsReal *pll_Ntt = nll->createProfile(*nFitSig);
  RooPlot *frameNtt = nFitSig->frame(Bins(20),Range(nFitSig->getVal()-2*nFitSig->getError(),nFitSig->getVal()+2*nFitSig->getError()),Title(""));
  nll->plotOn(frameNtt,ShiftToZero(),LineStyle(2));
  pll_Ntt->plotOn(frameNtt,ShiftToZero());
  frameNtt->SetMinimum(0);
  frameNtt->SetMaximum(2);
  frameNtt->Draw();
  */
  
  RooPlot *frame2b = x->frame();
  combData.plotOn(frame2b,Cut("sample==sample::2btag"),DrawOption("EP")); 
  //simPdf.plotOn(frame2b,Slice(sample,"2btag"),ProjWData(sample,combData),VisualizeError(*res,1),FillColor(kOrange),MoveToBack());
  simPdf.plotOn(frame2b,Slice(sample,"2btag"),ProjWData(sample,combData),LineColor(kBlue),MoveToBack());
  RooHist *pull2b = frame2b->pullHist();
  simPdf.plotOn(frame2b,Slice(sample,"2btag"),Components("qcd_pdf"),ProjWData(sample,combData),LineColor(kGreen+2),LineWidth(3),LineStyle(7));
  simPdf.plotOn(frame2b,Slice(sample,"2btag"),Components("ttbar_pdf_2btag"),ProjWData(sample,combData),DrawOption("FL"),LineColor(kRed),LineWidth(0),FillColor(kRed-10),MoveToBack());
  simPdf.plotOn(frame2b,Slice(sample,"2btag"),Components("bkg_pdf_2btag"),ProjWData(sample,combData),LineColor(kBlack),LineWidth(3),LineStyle(5)); 

  RooPlot *frame2bPull = x->frame();
  frame2bPull->addPlotable(pull2b,"p");

  TCanvas *can2b = new TCanvas("MassFit_2btag_"+ALIAS+"_"+CUT,"MassFit_2btag_"+ALIAS+"_"+CUT,800,704);
  can2b->cd(1)->SetBottomMargin(0.32);
  frame2b->GetYaxis()->SetTitle("Events");
  frame2b->GetXaxis()->SetTitle("");
  frame2b->GetXaxis()->SetLabelSize(0.0);
  frame2b->Draw();
  
  //can2b->GetListOfPrimitives()->Print();
  TLegend *leg = new TLegend(0.6,0.68,0.9,0.9);
  leg->SetHeader("Hadronic t#bar{t} decay");
  leg->AddEntry(can2b->GetPrimitive("h_combData_Cut[sample==sample::2btag]"),"Data","E1LP");
  leg->AddEntry(can2b->GetPrimitive("model_2b_Norm[mTop]"),"Fit model","L");
  leg->AddEntry(can2b->GetPrimitive("model_2b_Norm[mTop]_Comp[ttbar_pdf_2btag]"),"t#bar{t}","F");
  leg->AddEntry(can2b->GetPrimitive("model_2b_Norm[mTop]_Comp[qcd_pdf]"),"QCD multijets","L");
  leg->AddEntry(can2b->GetPrimitive("model_2b_Norm[mTop]_Comp[bkg_pdf_2btag]"),"Other backgrounds","L");
  leg->Draw();

  TPad *pad2b = new TPad("pad2b","pad2b",0.,0.,1.,1.);
  pad2b->SetTopMargin(0.7);
  pad2b->SetFillColor(0);
  pad2b->SetFillStyle(0);
  pad2b->Draw();
  pad2b->cd(0);
  pad2b->SetGridy();
  frame2bPull->SetMinimum(-5);
  frame2bPull->SetMaximum(5);
  frame2bPull->GetYaxis()->SetNdivisions(505);
  frame2bPull->GetXaxis()->SetTitleOffset(0.9);
  frame2bPull->GetYaxis()->SetTitleOffset(1.2);
  frame2bPull->GetYaxis()->SetTickLength(0.06);
  frame2bPull->GetYaxis()->SetTitleSize(0.05);
  frame2bPull->GetYaxis()->SetTitleSize(0.03);
  frame2bPull->GetYaxis()->SetLabelSize(0.03);
  frame2bPull->GetYaxis()->SetTitle("(Data-Fit)/Error");
  frame2bPull->GetXaxis()->SetTitle("m^{t} (GeV)");
  frame2bPull->Draw();
  
  gPad->RedrawAxis();

  //CMS_lumi(can2b,4,0);

  can2b->Print(TString::Format("%s/plots/SimpleMassFit/", year.Data())+TString(can2b->GetName())+".pdf");
  

  RooWorkspace *wOut = new RooWorkspace("w","workspace");
  wOut->import(*nFitQCD2b);
  wOut->import(*nFitSig2b);
  wOut->import(*yieldTT);
  wOut->writeToFile(TString::Format("%s/MassFitResults_",year.Data())+ALIAS+"_"+CUT+".root");

}

