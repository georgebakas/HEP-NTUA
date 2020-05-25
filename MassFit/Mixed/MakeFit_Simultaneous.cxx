	#include "TemplateConstants.h"

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

void MakeFit_Simultaneous(TString year = "2016", bool setConstant = false)
{
	gROOT->ForceStyle();

  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);

	initFilesMapping();

	TFile *file2 = TFile::Open(TString::Format("%s/Histo_Data_%s_100.root", year.Data(), year.Data()));
	TFile *file0 = TFile::Open(TString::Format("%s/Histo_Data_%s_100_Loose.root", year.Data(), year.Data()));
	TH1F* h0b = (TH1F*) file0->Get("hWt_mTop_0btag_expYield");
	TH1F* h2b = (TH1F*) file2->Get("hWt_mTop_2btag_expYield");
	//h0b->Rebin(2);
	//h2b->Rebin(2);

	//TFile *fileMC0 = TFile::Open(TString::Format("%s/Histo_TT_Mtt-700toInf_100_Loose.root", year.Data()));
	//TFile *fileMC2 = TFile::Open(TString::Format("%s/Histo_TT_Mtt-700toInf_100.root", year.Data()));
	TFile *fileMC0 = TFile::Open(TString::Format("%s/Histo_TT_NominalMC_100_Loose.root", year.Data()));
	TFile *fileMC2 = TFile::Open(TString::Format("%s/Histo_TT_NominalMC_100.root", year.Data()));
	TH1F* h0b_TT = (TH1F*)fileMC0->Get("hWt_mTop_0btag_expYield");
	TH1F* h2b_TT = (TH1F*)fileMC2->Get("hWt_mTop_2btag_expYield");
	float Ntt_expected = h2b_TT->Integral() + h0b_TT->Integral();

	TFile *fileSub0 = TFile::Open(TString::Format("%s/Histo_SubdominantBkgs_100_Loose.root", year.Data()));
	TFile *fileSub2 = TFile::Open(TString::Format("%s/Histo_SubdominantBkgs_100.root", year.Data()));
	TH1F* h0b_Bkg = (TH1F*)fileSub0->Get("hWt_mTop_0btag_expYield");
	TH1F* h2b_Bkg = (TH1F*)fileSub2->Get("hWt_mTop_2btag_expYield");

	TFile *fTemplatesBkg = TFile::Open(TString::Format("%s/templates_Bkg_100.root", year.Data()));
	TFile *fTemplatesSig = TFile::Open(TString::Format("%s/templates_Sig_100.root", year.Data()));
	RooWorkspace *wTemplatesBkg = (RooWorkspace*)fTemplatesBkg->Get("w");
	RooWorkspace *wTemplatesSig = (RooWorkspace*)fTemplatesSig->Get("w");

	RooRealVar *x = (RooRealVar*)wTemplatesSig->var("mTop");

	RooRealVar *kMassScale = (RooRealVar*)wTemplatesSig->var("kMassScale");
	RooRealVar *kMassResol = (RooRealVar*)wTemplatesSig->var("kMassResol");
	kMassScale->setConstant(false);
  kMassResol->setConstant(false);

	RooDataHist *roohist_data_0b = new RooDataHist("roohist_data_0b", "roohist_data_0b", *x, h0b);
	RooDataHist *roohist_data_2b = new RooDataHist("roohist_data_2b", "roohist_data_2b", *x, h2b);

	RooCategory sample("sample","sample");
	sample.defineType("0btag");
  sample.defineType("2btag");

  RooDataHist combData("combData","combData",*x, Index(sample),Import("0btag",*h0b),Import("2btag",*h2b));

	//Subdominant bkgs
	RooAbsPdf *pdf_bkg_0b = (RooAbsPdf*)wTemplatesBkg->pdf("bkg_pdf_0btag");
	RooAbsPdf *pdf_bkg_2b = (RooAbsPdf*)wTemplatesBkg->pdf("bkg_pdf_2btag");

	//QCD
	RooAbsPdf *pdf_qcd_0b = (RooAbsPdf*)wTemplatesBkg->pdf("qcd_pdf");
	RooAbsPdf *pdf_qcd_2b = (RooAbsPdf*)wTemplatesBkg->pdf("qcd_pdf");

	//QCD correction factor
	RooRealVar *kQCD2b = new RooRealVar("kQCD_2b", "kQCD_2b", 10e-4, -1, 1);
	kQCD2b->setConstant(false);
	RooFormulaVar qcdCor_2b("qcdCor_2b", "1+@0*@1", RooArgList(*x, *kQCD2b));

	//corected QCD
	RooEffProd pdf_qcdCor_2b("qcdCor_pdf_2b", "qcdCor_pdf_2b", *pdf_qcd_2b, qcdCor_2b);

	RooRealVar *nFitBkg0b = new RooRealVar("nFitBkg_0b", "nFitBkg_0b", h0b_Bkg->Integral(),0.9*h0b_Bkg->Integral(),1.1*h0b_Bkg->Integral());
	RooRealVar *nFitBkg2b = new RooRealVar("nFitBkg_2b", "nFitBkg_2b", h2b_Bkg->Integral(),0.9*h2b_Bkg->Integral(),1.1*h2b_Bkg->Integral());

	RooRealVar *nFitQCD0b = new RooRealVar("nFitQCD_0b", "nFitQCD_0b", 90000, 0, 1.2e+6);
	RooRealVar *nFitQCD2b = new RooRealVar("nFitQCD_2b", "nFitQCD_2b", 3000, 0, 1e+5);

	RooRealVar *nFitSig = new RooRealVar("nFitSig", "nFitSig", Ntt_expected, 0.5*Ntt_expected, 1.5*Ntt_expected);
	RooRealVar *nFitSig0b = new RooRealVar("nFitSig0b", "nFitSig0b", h0b_TT->Integral(), 0.6*h0b_TT->Integral(), 1.4*h0b_TT->Integral());
	RooRealVar *nFitSig2b = new RooRealVar("nFitSig2b", "nFitSig2b", h2b_TT->Integral(), 0.6*h2b_TT->Integral(), 1.4*h2b_TT->Integral());
	//RooRealVar *btagEff   = new RooRealVar("btagEff", "btagEff", floatConstants[TString::Format("bTagEff%s",year.Data())],0.1,1);
  RooRealVar *btagEff_0   = new RooRealVar("btagEff_0", "btagEff_0",floatConstants[TString::Format("bTagEff%s",year.Data())],0.4,1);
	RooRealVar *btagEff_2   = new RooRealVar("btagEff_2", "btagEff_2",floatConstants[TString::Format("bTagEff%s",year.Data())],0.4,1);
	if(setConstant)
	{
		btagEff_0->setConstant(false);
		btagEff_2->setConstant(true);
	}


	RooFormulaVar nSig0b("nSig_0b", "(1-@0)*(1-@0)*@1", RooArgList(*btagEff_0, *nFitSig));
	RooFormulaVar nSig2b("nSig_2b", "@0*@0*@1", RooArgList(*btagEff_2, *nFitSig));

	RooAbsPdf *pdf_signal_0b = (RooAbsPdf*)wTemplatesSig->pdf("ttbar_pdf_0btag");
	RooAbsPdf *pdf_signal_2b = (RooAbsPdf*)wTemplatesSig->pdf("ttbar_pdf_2btag");

	RooAddPdf *model_0b = new RooAddPdf("model_0btag", "model_0btag", RooArgList(*pdf_signal_0b, *pdf_qcd_0b, *pdf_bkg_0b), RooArgList(nSig0b, *nFitQCD0b, *nFitBkg0b));
	RooAddPdf *model_2b = new RooAddPdf("model_2btag", "model_2btag", RooArgList(*pdf_signal_2b, pdf_qcdCor_2b, *pdf_bkg_2b), RooArgList(nSig2b, *nFitQCD2b, *nFitBkg2b));

	RooSimultaneous simPdf("simPdf", "simPdf", sample);
	simPdf.addPdf(*model_0b, "0btag");
	simPdf.addPdf(*model_2b, "2btag");

	RooFitResult *res = simPdf.fitTo(combData, Save(), Extended(kTRUE));
	res->Print();

	std::cout<<"N0_observed = "<<nSig0b.getVal()<<", N2_observed = "<<nSig2b.getVal()<<std::endl;
	float Ntt_observed = nSig0b.getVal() + nSig2b.getVal();

	std::cout<<"Ntt expected: "<<Ntt_expected<<std::endl;
	std::cout<<"Ntt observed: "<<Ntt_observed<<std::endl;
	std::cout<<"Signal strength r: "<<Ntt_observed/Ntt_expected<<std::endl;
	std::cout<<"Singal strength r in 2btag: "<<nSig2b.getVal()/h2b_TT->Integral()<<endl;
	std::cout<<"Singal strength r in 0btag: "<<nSig0b.getVal()/h0b_TT->Integral()<<endl;

	TString CAT = "2btag";
	draw("2btag", year, x, combData, simPdf, sample);
	draw("0btag", year, x, combData, simPdf, sample);

	for(int i=0; i< canvases.size(); i++)
	{
		canvases[i]->Print(TString::Format("%s/plots/SimultaneousFit_2regions/%s.pdf", year.Data(), canvases[i]->GetName()), "pdf");
	}
	/*
	correlation(nFitSig, btagEff, res, "t#bar{t} events", "btag efficiency");
	correlation(kMassResol, btagEff, res, "kMassResol", "btag efficiency");
	correlation(kMassScale, btagEff, res, "kMassScale", "btag efficiency");

	correlation(nFitQCD2b, btagEff, res, "QCD events (2b category)", "btag efficiency");
	correlation(nFitQCD2b, kMassResol, res, "QCD events (2b category)", "kMassResol");
	correlation(nFitQCD2b, kMassScale, res, "QCD events (2b category)", "kMassScale");
	correlation(nFitQCD2b, kQCD1b, res, "QCD events (2b category)", "QCD_1b correction factor");
	correlation(nFitQCD2b, kQCD2b, res, "QCD events (2b category)", "QCD_2b correction factor");
	correlation(nFitQCD2b, nFitBkg0b, res, "QCD events (2b category)", "subdominant bkg events (0b category)");
	correlation(nFitQCD2b, nFitBkg1b, res, "QCD events (2b category)", "subdominant bkg events (1b category)");
	correlation(nFitQCD2b, nFitBkg2b, res, "QCD events (2b category)", "subdominant bkg events (2b category)");
	correlation(nFitQCD2b, nFitQCD0b, res, "QCD events (2b category)", "QCD events (0b category)");
	correlation(nFitQCD2b, nFitQCD1b, res, "QCD events (2b category)", "QCD events (1b category)");
	correlation(nFitQCD2b, nFitSig, res, "QCD events (2b category)", "ttbar events (all categories)");

	for(int i=0; i<correlationCanvases.size(); i++)
	{
	  correlationCanvases[i]->Print(TString::Format("%s/correlationPlots/%s.pdf", year.Data(), correlationCanvases[i]->GetName()), "pdf");
	}
	*/
}
