//#include "../../CMS_lumi.C"

using namespace RooFit;

std::vector<TCanvas*> canvases;
std::vector<TCanvas*> correlationCanvases;
void draw(TString CAT, TString year, RooRealVar *x, RooDataHist &combData, RooSimultaneous &simPdf, RooCategory &sample)
{
	TCanvas *can = new TCanvas(TString::Format("MassFit_s_%s_%s", CAT.Data(), year.Data()),
	                           TString::Format("MassFit_s_%s_%s", CAT.Data(), year.Data()), 800, 704);
	can->cd()->SetBottomMargin(0.32);

	RooPlot* frame2b = x->frame();
  frame2b->GetYaxis()->SetTitle("Events");
	frame2b->GetYaxis()->SetTitleOffset(1.5);
  frame2b->GetXaxis()->SetTitle("");
  frame2b->GetXaxis()->SetLabelSize(0.0);
	combData.plotOn(frame2b, Cut(TString::Format("sample==sample::%s", CAT.Data())));
	simPdf.plotOn(frame2b, Slice(sample, CAT), ProjWData(sample, combData));
	RooHist *pull2b = frame2b->pullHist();
	simPdf.plotOn(frame2b, Slice(sample, CAT), Components("qcd_pdf"), ProjWData(sample, combData), LineColor(kGreen+2), LineWidth(3), LineStyle(7));
	simPdf.plotOn(frame2b, Slice(sample, CAT), Components(TString::Format("ttbar_pdf_%s", CAT.Data())), ProjWData(sample, combData), DrawOption("FL"), LineColor(kRed), LineWidth(0), FillColor(kRed-10), MoveToBack());
	simPdf.plotOn(frame2b, Slice(sample, CAT), Components(TString::Format("bkg_pdf_%s", CAT.Data())), ProjWData(sample, combData), LineColor(kBlack), LineWidth(3), LineStyle(5));
	frame2b->Draw();

	RooPlot* frame2bPull = x->frame();
	frame2bPull->addPlotable(pull2b, "p");

	TLegend *leg = new TLegend(0.6, 0.68, 0.9, 0.8999);
	leg->SetHeader("Hadronic t#bar{t} decay");
  leg->AddEntry(can->GetPrimitive(TString::Format("h_combData_Cut[sample==sample::%s]", CAT.Data())),"Data","E1LP");
  leg->AddEntry(can->GetPrimitive(TString::Format("model_%s_Norm[mTop]", CAT.Data())),"Fit model","L");
  leg->AddEntry(can->GetPrimitive(TString::Format("model_%s_Norm[mTop]_Comp[ttbar_pdf_%s]", CAT.Data(), CAT.Data())),"t#bar{t}","F");
  leg->AddEntry(can->GetPrimitive(TString::Format("model_%s_Norm[mTop]_Comp[qcd_pdf]", CAT.Data())),"QCD multijets","L");
  leg->AddEntry(can->GetPrimitive(TString::Format("model_%s_Norm[mTop]_Comp[bkg_pdf_%s]", CAT.Data(), CAT.Data())),"Other backgrounds","L");
	leg->SetMargin(0.3);
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
	can->cd(1);
	canvases.push_back(can);
}

void correlation(RooRealVar *x, RooRealVar *y, RooFitResult *res, TString titleX, TString titleY)
{
	TCanvas *can = new TCanvas(TString::Format("Correlation_%s_vs_%s", x->GetName(), y->GetName()),
	                           TString::Format("Correlation_%s_vs_%s", x->GetName(), y->GetName()),
														 900, 600);
  RooPlot *frame = new RooPlot(*x, *y, x->getVal()-1.5*x->getError(),x->getVal()+1.5*x->getError(),y->getVal()-1.5*y->getError(),y->getVal()+1.5*y->getError());
  res->plotOn(frame,*x,*y,"ME12ABHV");
  frame->GetXaxis()->SetTitle(titleX);
  frame->GetYaxis()->SetTitle(titleY);
	frame->GetYaxis()->SetTitleOffset(1.5);
  frame->Draw();

	correlationCanvases.push_back(can);
}

void MassFitNew(TString year = "2016", TString selection = "probe",TString ALIAS="",TString CUT="", int REBIN= 5)
{
  TString selectedRegion;
  if(selection.EqualTo("probe")) selectedRegion = "hSRBTightAndProbe";
  else if (selection.EqualTo("SR")) selectedRegion = "hSRBTightAndSR"; 
  gROOT->ForceStyle();

  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);

  //Take SR data from Medium WP
  TFile *inf = TFile::Open(TString::Format("%s/TagAndProbeHisto_Data_%s_100_reduced_UnequalBinning.root", year.Data(), year.Data()));
  TH1F *h2b  = (TH1F*)inf->Get(TString::Format("%s_mTop_expYield",selectedRegion.Data()));
  //h2b->Rebin(2);
  // -----------------------------------------
  TFile *fTemplatesBkg = TFile::Open(TString::Format("%s/templates_Bkg_%s_100.root", year.Data(),selectedRegion.Data()));
  TFile *fTemplatesSig = TFile::Open(TString::Format("%s/templates_Sig_%s_100.root",year.Data(), selectedRegion.Data()));
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

  RooAbsPdf *pdf_signal_2b = (RooAbsPdf*)wTemplatesSig->pdf("ttbar_pdf_2btag");
  //RooAbsPdf *pdf_signal_0b = (RooAbsPdf*)wTemplatesSig->pdf("ttbar_pdf_0btag");

  //---- QCD correction factor ---------------------------

  float sP, min, max;
  //sP = 1e-3;
  min = -1;
  max = 10;
  sP = 0.0025;
  RooRealVar *kQCD2b_0 = new RooRealVar("kQCD_2b","kQCD_2b", sP);//, min, max);
  RooRealVar *mBar = new RooRealVar("mBar", "mBar", 175, 50, 300);
  mBar->setConstant(true);
  kQCD2b_0->setConstant(false);

  //RooFormulaVar qcdCor_2b("qcdCor","(1+@0*@1)/(1+@1*@2)",RooArgList(*x,*kQCD2b_0,*mBar));
  RooFormulaVar qcdCor_2b("qcdCor_2b","(1+@0*@1)",RooArgList(*x,*kQCD2b_0));
  //---- corrected QCD -----------------------------------

  RooEffProd pdf_qcdCor_2b("qcdCor_pdf_2b","qcdCor_pdf_2b",*pdf_qcd_2b,qcdCor_2b);

  RooRealVar *nFitBkg2b = new RooRealVar("nFitBkg_","nFitBkg_",400,0,10e+3);

  RooRealVar *nFitQCD2b = new RooRealVar("nFitQCD_","nFitQCD_",10000,0,10e+4);

  RooRealVar *nFitSig2b = new RooRealVar("nFitSig","nFitSig",2000,100,10e+4);

  RooAddPdf *model_2b = new RooAddPdf("model_2b","model_2b",RooArgList(*pdf_signal_2b,*pdf_qcd_2b,*pdf_bkg_2b),RooArgList(*nFitSig2b,*nFitQCD2b,*nFitBkg2b));

  RooSimultaneous simPdf("simPdf","simPdf",sample);
  simPdf.addPdf(*model_2b,"2btag");

  RooFitResult *res = simPdf.fitTo(combData,RooFit::Save(),RooFit::Extended(kTRUE));

	//Signal strengh and error propagation:
	float sigError = TMath::Sqrt(TMath::Power(nFitSig2b->getError()/yieldTT->getVal(),2) + TMath::Power(nFitSig2b->getVal()*nFitSig2b->getError()/TMath::Power(yieldTT->getVal(),2),2));
  cout<<"Signal strength: r = "<<nFitSig2b->getVal()/yieldTT->getVal()<<" ± "<< sigError<<endl;

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

  can2b->Print(TString::Format("%s/plots/SimpleMassFit/", year.Data())+TString(can2b->GetName())+selectedRegion+".pdf");


  RooWorkspace *wOut = new RooWorkspace("w","workspace");
  //wOut->import(*pdf_qcd_2b);
  wOut->import(pdf_qcdCor_2b);
  wOut->import(*nFitQCD2b);
  wOut->import(*nFitSig2b);
  wOut->import(*yieldTT);
  wOut->writeToFile(TString::Format("%s/MassFitResults_%s",year.Data(),selectedRegion.Data())+ALIAS+"_"+CUT+".root");
  //wOut->writeToFile(TString::Format("%s/MassFitResultsCorrectedFit_",year.Data())+ALIAS+"_"+CUT+".root");
  //wOut->writeToFile(TString::Format("%s/MassFitResultsNoCorrection_",year.Data())+ALIAS+"_"+CUT+".root");

  /*
  correlation(kQCD2b_0, nFitBkg2b, res, "kQCD2b ", "nFitBkg2b");
  correlation(kQCD2b_0, nFitQCD2b, res, "kQCD2b", "nFitQCD2b");
  correlation(kQCD2b_0, nFitSig2b, res, "kQCD2b", "nFitSig2b");
  correlation(kQCD2b_0, kMassResol, res, "kQCD2b", "kMassResol");
  correlation(kQCD2b_0, kMassScale, res, "kQCD2b", "kMassScale");

  for(int i=0; i<correlationCanvases.size(); i++)
  {
	  correlationCanvases[i]->Print(TString::Format("%s/correlationPlots/%s.pdf", year.Data(), correlationCanvases[i]->GetName()), "pdf");
  }	*/
}
