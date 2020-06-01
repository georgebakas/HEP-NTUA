using namespace RooFit;

void shapeComparison_CR(TString year)
{
  gROOT->ForceStyle();

  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);

  //comparison of 3 things:
  //template without any corrections

  TFile *fTemplatesBkg = TFile::Open(TString::Format("%s/templates_Bkg_100.root", year.Data()));
  TFile *fTemplatesSig = TFile::Open(TString::Format("%s/templates_Sig_100.root", year.Data()));
  RooWorkspace *wTemplatesBkg = (RooWorkspace*)fTemplatesBkg->Get("w");
  RooWorkspace *wTemplatesSig = (RooWorkspace*)fTemplatesSig->Get("w");


  RooRealVar *x = (RooRealVar*)wTemplatesBkg->var("mTop");
  RooAbsPdf *pdf_qcd_2b = (RooAbsPdf*)wTemplatesBkg->pdf("qcd_pdf");
  //now the pdf of the corrected qcd

  TFile *fCorrectedFit = TFile::Open(TString::Format("%s/MassFitResults__.root", year.Data()));
  RooWorkspace *wTemplatesBkgCorr = (RooWorkspace*)fCorrectedFit->Get("w");
  RooAbsPdf *pdf_qcdCor_2b = (RooAbsPdf*)wTemplatesBkgCorr->pdf("qcdCor_pdf_2b");

  TCanvas *can2b = new TCanvas(TString::Format("shapeComparison_CR_%s",year.Data()),
                               TString::Format("shapeComparison_CR_%s",year.Data()),800,704);
  TLegend *leg = new TLegend(0.6,0.68,0.9,0.9);

  RooPlot *testFrame = x->frame();
  pdf_qcd_2b->plotOn(testFrame, RooFit::Name("noCorr"), LineColor(kGreen+2));
  pdf_qcdCor_2b->plotOn(testFrame, RooFit::Name("Corrected"), LineColor(kOrange+1));
  testFrame->Draw();

  can2b->Print(TString::Format("%s/plots/shapeComparisonBeforeAfterCorrection/"+TString(can2b->GetName())+"_data.pdf", year.Data()));
}
