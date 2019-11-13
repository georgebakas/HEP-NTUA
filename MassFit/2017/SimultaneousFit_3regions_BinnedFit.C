using namespace RooFit;
void SimultaneousFit_3regions_new(int REBIN =2)
{
  gROOT->ForceStyle();
  TString CUT = "Cut5";
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);
  
  TFile *inf = TFile::Open("Histo_JetHT_Run2017-31Mar2018_New_100.root");
  TH1F *h0b  = (TH1F*)inf->Get("hWt_mTop_0btag");
  TH1F *h1b  = (TH1F*)inf->Get("hWt_mTop_1btag");
  TH1F *h2b  = (TH1F*)inf->Get("hWt_mTop_2btag");
  h0b->Rebin(REBIN);
  h2b->Rebin(REBIN);
  h1b->Rebin(REBIN);
  // -----------------------------------------
  const float LUMI = 41530;

  TFile *infTT = TFile::Open("Histo_TT_Mtt-700toInf_TuneCP5_13TeV-powheg-pythia8_New_100.root");
  TH1F *h2b_TT = (TH1F*)infTT->Get("hWt_mTop_2btag_expYield");
  TH1F *h1b_TT = (TH1F*)infTT->Get("hWt_mTop_1btag_expYield");
  TH1F *h0b_TT = (TH1F*)infTT->Get("hWt_mTop_0btag_expYield");
  float Ntt_expected = h2b_TT->Integral() + h1b_TT->Integral() + h0b_TT->Integral();

  TFile *infSubBkg = TFile::Open("Histo_SubdominantBkgs_New_100.root");
  TH1F *h2b_Bkg = (TH1F*)infSubBkg->Get("hWt_mTop_2btag_expYield");
  TH1F *h1b_Bkg = (TH1F*)infSubBkg->Get("hWt_mTop_1btag_expYield");
  TH1F *h0b_Bkg = (TH1F*)infSubBkg->Get("hWt_mTop_0btag_expYield");

  TFile *infQCDBkg = TFile::Open("Histo_QCD_HT300toInf_TuneCP5_13TeV-madgraph-pythia8_New_100.root");
  TH1F *h2b_QCDBkg = (TH1F*)infQCDBkg->Get("hWt_mTop_2btag_expYield");
  TH1F *h1b_QCDBkg = (TH1F*)infQCDBkg->Get("hWt_mTop_1btag_expYield");
  TH1F *h0b_QCDBkg = (TH1F*)infQCDBkg->Get("hWt_mTop_0btag_expYield");

  
  TFile *fTemplatesBkg = TFile::Open("templates_Bkg_100_NewTest.root");
  TFile *fTemplatesSig = TFile::Open("templates_Sig_100.root");
  RooWorkspace *wTemplatesBkg = (RooWorkspace*)fTemplatesBkg->Get("w");
  RooWorkspace *wTemplatesSig = (RooWorkspace*)fTemplatesSig->Get("w");
  
  RooRealVar *x = (RooRealVar*)wTemplatesSig->var("mTop");
  
  RooRealVar *kMassScale = (RooRealVar*)wTemplatesSig->var("kMassScale");
  RooRealVar *kMassResol = (RooRealVar*)wTemplatesSig->var("kMassResol");
  kMassScale->setConstant(false);
  kMassResol->setConstant(false);
  
  RooDataHist *roohist_data_0b = new RooDataHist("roohist_data_0b","roohist_data_0b",*x,h0b);
  RooDataHist *roohist_data_2b = new RooDataHist("roohist_data_2b","roohist_data_2b",*x,h2b);
  RooDataHist *roohist_data_1b = new RooDataHist("roohist_data_1b","roohist_data_1b",*x,h1b);

  RooPlot *fram1 = x->frame();
  roohist_data_1b->plotOn(fram1);
  fram1->Draw();

  RooCategory sample("sample","sample");
  sample.defineType("0btag");
  sample.defineType("1btag");
  sample.defineType("2btag");
  
  RooDataHist combData("combData","combData",*x,Index(sample),Import("0btag",*h0b),Import("1btag",*h1b),Import("2btag",*h2b));
  //subdominant bkg tempates 
  RooAbsPdf *pdf_bkg_0b = (RooAbsPdf*)wTemplatesBkg->pdf("bkg_pdf_0btag");
  RooAbsPdf *pdf_bkg_2b = (RooAbsPdf*)wTemplatesBkg->pdf("bkg_pdf_2btag");
  RooAbsPdf *pdf_bkg_1b = (RooAbsPdf*)wTemplatesBkg->pdf("bkg_pdf_1btag");
  //CR templates taken from data for QCD bkg -->This is the difference, I will give as input an RooHistPdf
  //RooAbsPdf *pdf_qcd_0b = (RooAbsPdf*)wTemplatesBkg->pdf("qcd_pdf");
  //RooAbsPdf *pdf_qcd_1b = (RooAbsPdf*)wTemplatesBkg->pdf("qcd_pdf");
  //RooAbsPdf *pdf_qcd_2b = (RooAbsPdf*)wTemplatesBkg->pdf("qcd_pdf");

  RooHistPdf *pdf_qcd_0b = (RooHistPdf*)wTemplatesBkg->pdf("QCDhist_pdf");
  RooHistPdf *pdf_qcd_1b = (RooHistPdf*)wTemplatesBkg->pdf("QCDhist_pdf");
  RooHistPdf *pdf_qcd_2b = (RooHistPdf*)wTemplatesBkg->pdf("QCDhist_pdf");
  
  //---- QCD correction factor ---------------------------
  RooRealVar kQCD1b("kQCD_1b","kQCD_1b",1e-3,-1,1);
  RooRealVar kQCD2b("kQCD_2b","kQCD_2b",1e-3,-1,1);
  kQCD1b.setConstant(false);
  kQCD2b.setConstant(false);
  RooFormulaVar qcdCor_1b("qcdCor_1b","1+@0*@1",RooArgList(*x,kQCD1b)); 
  RooFormulaVar qcdCor_2b("qcdCor_2b","1+@0*@1",RooArgList(*x,kQCD2b));

  //---- corrected QCD -----------------------------------
  RooEffProd pdf_qcdCor_1b("qcdCor_pdf_1b","qcdCor_pdf_1b",*pdf_qcd_1b,qcdCor_1b);
  RooEffProd pdf_qcdCor_2b("qcdCor_pdf_2b","qcdCor_pdf_2b",*pdf_qcd_2b,qcdCor_2b);
  
  RooRealVar *nFitBkg0b = new RooRealVar("nFitBkg_0b","nFitBkg_0b",h0b_Bkg->Integral(),0.7*h0b_Bkg->Integral() ,1.3*h0b_Bkg->Integral()); //oti dinei to MC ± 30%
  RooRealVar *nFitBkg1b = new RooRealVar("nFitBkg_1b","nFitBkg_1b",h1b_Bkg->Integral(),0.7*h1b_Bkg->Integral() ,1.3*h1b_Bkg->Integral());
  RooRealVar *nFitBkg2b = new RooRealVar("nFitBkg_2b","nFitBkg_2b",h2b_Bkg->Integral(),0.7*h2b_Bkg->Integral() ,1.3*h2b_Bkg->Integral());

  RooRealVar *nFitQCD0b = new RooRealVar("nFitQCD_0b","nFitQCD_0b",90000,0,1.2e+6);
  RooRealVar *nFitQCD1b = new RooRealVar("nFitQCD_1b","nFitQCD_1b",35000,0,1e+5);
  RooRealVar *nFitQCD2b = new RooRealVar("nFitQCD_2b","nFitQCD_2b",3000,0,1e+4);  
/*
  cout<<"h0b_QCDBkg: "<<h0b_QCDBkg->Integral()<<endl;
  cout<<"h1b_QCDBkg: "<<h1b_QCDBkg->Integral()<<endl;
  cout<<"h2b_QCDBkg: "<<h2b_QCDBkg->Integral()<<endl;


  RooRealVar *nFitQCD0b = new RooRealVar("nFitQCD_0b","nFitQCD_0b",h0b_QCDBkg->Integral(),0.5*h0b_QCDBkg->Integral(),1.5*h0b_QCDBkg->Integral());
  RooRealVar *nFitQCD1b = new RooRealVar("nFitQCD_1b","nFitQCD_1b",h1b_QCDBkg->Integral(),0.5*h1b_QCDBkg->Integral(),1.5*h1b_QCDBkg->Integral());
  RooRealVar *nFitQCD2b = new RooRealVar("nFitQCD_2b","nFitQCD_2b",h2b_QCDBkg->Integral(),0.5*h2b_QCDBkg->Integral(),1.5*h2b_QCDBkg->Integral()); 
*/

  RooRealVar *nFitSig   = new RooRealVar("nFitSig","nFitSig",Ntt_expected,0.5*Ntt_expected,1.5*Ntt_expected);
  RooRealVar *nFitSig0b = new RooRealVar("nFitSig0b","nFitSig0b",h0b_TT->Integral(),0.6* h0b_TT->Integral(),1.4*h0b_TT->Integral());
  RooRealVar *nFitSig1b = new RooRealVar("nFitSig1b","nFitSig1b",h1b_TT->Integral(),0.6* h1b_TT->Integral(),1.4*h1b_TT->Integral());
  RooRealVar *nFitSig2b = new RooRealVar("nFitSig2b","nFitSig2b",h2b_TT->Integral(),0.6* h2b_TT->Integral(),1.4*h2b_TT->Integral());
  RooRealVar *btagEff   = new RooRealVar("btagEff","btagEff",0.605622,0.4,0.8);
  //RooRealVar *btagEff = new RooRealVar("btagEff", "btagEff", 0.605622);
  btagEff->setConstant(true);

  RooFormulaVar nSig0b("nSig_0b","(1-@0)*(1-@0)*@1",RooArgList(*btagEff,*nFitSig)); 
  RooFormulaVar nSig2b("nSig_2b","@0*@0*@1",RooArgList(*btagEff,*nFitSig));
  RooFormulaVar nSig1b("nSig_1b","2*(1-@0)*@0*@1", RooArgList(*btagEff, *nFitSig));

  RooAbsPdf *pdf_signal_0b = (RooAbsPdf*)wTemplatesSig->pdf("ttbar_pdf_0btag");
  RooAbsPdf *pdf_signal_1b = (RooAbsPdf*)wTemplatesSig->pdf("ttbar_pdf_1btag");
  RooAbsPdf *pdf_signal_2b = (RooAbsPdf*)wTemplatesSig->pdf("ttbar_pdf_2btag");

  RooAddPdf *model_0b = new RooAddPdf("model_0b","model_0b",RooArgList(*pdf_signal_0b,*pdf_qcd_0b,*pdf_bkg_0b),RooArgList(nSig0b,*nFitQCD0b,*nFitBkg0b)); 
  //RooAddPdf *model_0b = new RooAddPdf("model_0b","model_0b",RooArgList(*pdf_signal_0b,*pdf_qcd_0b),RooArgList(nSig0b,*nFitQCD0b)); 
  RooAddPdf *model_1b = new RooAddPdf("model_1b","model_1b",RooArgList(*pdf_signal_1b,pdf_qcdCor_1b,*pdf_bkg_1b),RooArgList(nSig1b,*nFitQCD1b,*nFitBkg1b));
  RooAddPdf *model_2b = new RooAddPdf("model_2b","model_2b",RooArgList(*pdf_signal_2b,pdf_qcdCor_2b,*pdf_bkg_2b),RooArgList(nSig2b,*nFitQCD2b,*nFitBkg2b));

  RooSimultaneous simPdf("simPdf","simPdf",sample);
  simPdf.addPdf(*model_0b,"0btag");
  simPdf.addPdf(*model_1b,"1btag");
  simPdf.addPdf(*model_2b,"2btag");
    
  RooFitResult *res = simPdf.fitTo(combData,RooFit::Save(),RooFit::Extended(kTRUE));
  res->Print();


  cout<<"N0_observed = "<<nSig0b.getVal()<<", N2_observed = "<<nSig2b.getVal()<<", N1_observed = "<<nSig1b.getVal()<<endl;
  float Ntt_observed = nSig0b.getVal()+ nSig1b.getVal() + nSig2b.getVal();
  //cout<<"N0 = "<<nFitSig0b->getVal()<<", N2 = "<<nFitSig2b->getVal()<<endl;



  //expected Ntt is the sum of all 0btag and 1 btag and 2btag
  //observed you have to find using nSig2b and nSig0b and nSig1b
  cout<<"Ntt expected = "<<Ntt_expected<<endl;
  cout<<"Ntt observed = "<<Ntt_observed<<endl;
  cout<<"r = "<< Ntt_observed/Ntt_expected<<endl;
  

  
  TCanvas *canEll1 = new TCanvas("Correlation_BTagEff_vs_NTT_"+CUT,"Correlation_BTagEff_vs_NTT_"+CUT,900,600);
  RooPlot *frameEll1 = new RooPlot(*nFitSig,*btagEff,nFitSig->getVal()-1.5*nFitSig->getError(),nFitSig->getVal()+1.5*nFitSig->getError(),btagEff->getVal()-1.5*btagEff->getError(),btagEff->getVal()+1.5*btagEff->getError());
  res->plotOn(frameEll1,*nFitSig,*btagEff,"ME12ABHV");
  frameEll1->GetXaxis()->SetTitle("t#bar{t} events");
  frameEll1->GetYaxis()->SetTitle("btag efficiency");
  frameEll1->Draw();

  TCanvas *canEll3 = new TCanvas("Correlation_BTagEff_vs_kMassResol_"+CUT,"Correlation_BTagEff_vs_kMassResol_"+CUT,900,600);
  RooPlot *frameEll3 = new RooPlot(*kMassResol,*btagEff,kMassResol->getVal()-1.5*kMassResol->getError(),kMassResol->getVal()+1.5*kMassResol->getError(),btagEff->getVal()-1.5*btagEff->getError(),btagEff->getVal()+1.5*btagEff->getError());
  res->plotOn(frameEll3,*kMassResol,*btagEff,"ME12ABHV");
  frameEll3->GetXaxis()->SetTitle("kMassResol");
  frameEll3->GetYaxis()->SetTitle("btag efficiency");
  frameEll3->Draw();

  TCanvas *canEll4 = new TCanvas("Correlation_BTagEff_vs_kMassScale_"+CUT,"Correlation_BTagEff_vs_kMassScale_"+CUT,900,600);
  RooPlot *frameEll4 = new RooPlot(*kMassScale,*btagEff,kMassScale->getVal()-1.5*kMassScale->getError(),kMassScale->getVal()+1.5*kMassScale->getError(),btagEff->getVal()-1.5*btagEff->getError(),btagEff->getVal()+1.5*btagEff->getError());
  res->plotOn(frameEll4,*kMassScale,*btagEff,"ME12ABHV");
  frameEll4->GetXaxis()->SetNdivisions(505);
  frameEll4->GetXaxis()->SetTitle("kMassScale");
  frameEll4->GetYaxis()->SetTitle("btag efficiency");
  frameEll4->Draw();

  /*
      Floating Parameter
  --------------------  
               btagEff  if this runs as a free parameter in the fit 
            kMassResol    
            kMassScale    
               kQCD_1b    
               kQCD_2b    
            nFitBkg_0b    
            nFitBkg_1b    
            nFitBkg_2b    
            nFitQCD_0b    
            nFitQCD_1b    
            nFitQCD_2b    
               nFitSig    
  */

  //correlation plots for NQCD2b vs all the fit parameters mentioned above in the array
  TCanvas *canEll2 = new TCanvas("Correlation_NQCD2b_vs_BTagEff"+CUT,"Correlation_BTagEff_vs_NQCD2b_"+CUT,800,600);
  RooPlot *frameEll2 = new RooPlot(*nFitQCD2b,*btagEff,nFitQCD2b->getVal()-1.5*nFitQCD2b->getError(),nFitQCD2b->getVal()+1.5*nFitQCD2b->getError(),btagEff->getVal()-1.5*btagEff->getError(),btagEff->getVal()+1.5*btagEff->getError());
  res->plotOn(frameEll2,*nFitQCD2b,*btagEff,"ME12ABHV");
  frameEll2->GetXaxis()->SetTitle("QCD events (2b category)");
  frameEll2->GetYaxis()->SetTitle("btag efficiency");
  frameEll2->Draw();

  TCanvas *canEll5 = new TCanvas("Correlation_NQCD2b_vs_kMassResol", "Correlation_NQCD2b_vs_kMassResol", 800,600);
  RooPlot *frameEll5 = new RooPlot(*nFitQCD2b,*kMassResol,nFitQCD2b->getVal()-1.5*nFitQCD2b->getError(),nFitQCD2b->getVal()+1.5*nFitQCD2b->getError(),kMassResol->getVal()-1.5*kMassResol->getError(),kMassResol->getVal()+1.5*kMassResol->getError());
  res->plotOn(frameEll5,*nFitQCD2b,*kMassResol,"ME12ABHV");
  frameEll5->GetXaxis()->SetTitle("QCD events (2b category)");
  frameEll5->GetYaxis()->SetTitle("kMassResol");
  frameEll5->Draw();

  TCanvas *canEll6 = new TCanvas("Correlation_NQCD2b_vs_kMassScale", "Correlation_NQCD2b_vs_kMassScale", 800,600);
  RooPlot *frameEll6 = new RooPlot(*nFitQCD2b,*kMassScale,nFitQCD2b->getVal()-1.5*nFitQCD2b->getError(),nFitQCD2b->getVal()+1.5*nFitQCD2b->getError(),kMassScale->getVal()-1.5*kMassScale->getError(),kMassScale->getVal()+1.5*kMassScale->getError());
  res->plotOn(frameEll6,*nFitQCD2b,*kMassScale,"ME12ABHV");
  frameEll6->GetXaxis()->SetTitle("QCD events (2b category)");
  frameEll6->GetYaxis()->SetTitle("kMassScale");
  frameEll6->Draw();

  TCanvas *canEll7 = new TCanvas("Correlation_NQCD2b_vs_kQCD_1b", "Correlation_NQCD2b_vs_kQCD_1b", 800,600);
  RooPlot *frameEll7 = new RooPlot(*nFitQCD2b,kQCD1b,nFitQCD2b->getVal()-1.5*nFitQCD2b->getError(),nFitQCD2b->getVal()+1.5*nFitQCD2b->getError(),kQCD1b.getVal()-1.5*kQCD1b.getError(),kQCD1b.getVal()+1.5*kQCD1b.getError());
  res->plotOn(frameEll7,*nFitQCD2b,kQCD1b,"ME12ABHV");
  frameEll7->GetXaxis()->SetTitle("QCD events (2b category)");
  frameEll7->GetYaxis()->SetTitle("QCD_1b correction factor");
  frameEll7->Draw();

  TCanvas *canEll8 = new TCanvas("Correlation_NQCD2b_vs_kQCD_2b", "Correlation_NQCD2b_vs_kQCD_1b", 800,600);
  RooPlot *frameEll8 = new RooPlot(*nFitQCD2b,kQCD2b,nFitQCD2b->getVal()-1.5*nFitQCD2b->getError(),nFitQCD2b->getVal()+1.5*nFitQCD2b->getError(),kQCD2b.getVal()-1.5*kQCD2b.getError(),kQCD2b.getVal()+1.5*kQCD2b.getError());
  res->plotOn(frameEll8,*nFitQCD2b,kQCD2b,"ME12ABHV");
  frameEll8->GetXaxis()->SetTitle("QCD events (2b category)");
  frameEll8->GetYaxis()->SetTitle("QCD_2b correction factor");
  frameEll8->Draw();

  TCanvas *canEll9 = new TCanvas("Correlation_NQCD2b_vs_nFitBkg_0b", "Correlation_NQCD2b_vs_nFitBkg_0b", 800,600);
  RooPlot *frameEll9 = new RooPlot(*nFitQCD2b,*nFitBkg0b,nFitQCD2b->getVal()-1.5*nFitQCD2b->getError(),nFitQCD2b->getVal()+1.5*nFitQCD2b->getError(),nFitBkg0b->getVal()-1.5*nFitBkg0b->getError(),nFitBkg0b->getVal()+1.5*nFitBkg0b->getError());
  res->plotOn(frameEll9,*nFitQCD2b,*nFitBkg0b,"ME12ABHV");
  frameEll9->GetXaxis()->SetTitle("QCD events (2b category)");
  frameEll9->GetYaxis()->SetTitle("subdominant bkg events (0b category)");
  frameEll9->Draw();

  TCanvas *canEll10 = new TCanvas("Correlation_NQCD2b_vs_nFitBkg_1b", "Correlation_NQCD2b_vs_nFitBkg_1b", 800,600);
  RooPlot *frameEll10 = new RooPlot(*nFitQCD2b,*nFitBkg1b,nFitQCD2b->getVal()-1.5*nFitQCD2b->getError(),nFitQCD2b->getVal()+1.5*nFitQCD2b->getError(),nFitBkg1b->getVal()-1.5*nFitBkg1b->getError(),nFitBkg1b->getVal()+1.5*nFitBkg1b->getError());
  res->plotOn(frameEll10,*nFitQCD2b,*nFitBkg1b,"ME12ABHV");
  frameEll10->GetXaxis()->SetTitle("QCD events (2b category)");
  frameEll10->GetYaxis()->SetTitle("subdominant bkg events (1b category)");
  frameEll10->Draw();

  TCanvas *canEll11 = new TCanvas("Correlation_NQCD2b_vs_nFitBkg_2b", "Correlation_NQCD2b_vs_nFitBkg_2b", 800,600);
  RooPlot *frameEll11 = new RooPlot(*nFitQCD2b,*nFitBkg2b,nFitQCD2b->getVal()-1.5*nFitQCD2b->getError(),nFitQCD2b->getVal()+1.5*nFitQCD2b->getError(),nFitBkg2b->getVal()-1.5*nFitBkg2b->getError(),nFitBkg2b->getVal()+1.5*nFitBkg2b->getError());
  res->plotOn(frameEll11,*nFitQCD2b,*nFitBkg2b,"ME12ABHV");
  frameEll11->GetXaxis()->SetTitle("QCD events (2b category)");
  frameEll11->GetYaxis()->SetTitle("subdominant bkg events (2b category)");
  frameEll11->Draw();

  TCanvas *canEll12 = new TCanvas("Correlation_NQCD2b_vs_nFitQCD0b", "Correlation_NQCD2b_vs_nFitQCD0b", 800,600);
  RooPlot *frameEll12 = new RooPlot(*nFitQCD2b,*nFitQCD0b,nFitQCD2b->getVal()-1.5*nFitQCD2b->getError(),nFitQCD2b->getVal()+1.5*nFitQCD2b->getError(),nFitQCD0b->getVal()-1.5*nFitQCD0b->getError(),nFitQCD0b->getVal()+1.5*nFitQCD0b->getError());
  res->plotOn(frameEll12,*nFitQCD2b,*nFitQCD0b,"ME12ABHV");
  frameEll12->GetXaxis()->SetTitle("QCD events (2b category)");
  frameEll12->GetYaxis()->SetTitle("QCD events (0b category)");
  frameEll12->Draw();

  TCanvas *canEll13 = new TCanvas("Correlation_NQCD2b_vs_nFitQCD1b", "Correlation_NQCD2b_vs_nFitQCD1b", 800,600);
  RooPlot *frameEll13 = new RooPlot(*nFitQCD2b,*nFitQCD1b,nFitQCD2b->getVal()-1.5*nFitQCD2b->getError(),nFitQCD2b->getVal()+1.5*nFitQCD2b->getError(),nFitQCD1b->getVal()-1.5*nFitQCD1b->getError(),nFitQCD1b->getVal()+1.5*nFitQCD1b->getError());
  res->plotOn(frameEll13,*nFitQCD2b,*nFitQCD1b,"ME12ABHV");
  frameEll13->GetXaxis()->SetTitle("QCD events (2b category)");
  frameEll13->GetYaxis()->SetTitle("QCD events (1b category)");
  frameEll13->Draw();

  TCanvas *canEll14 = new TCanvas("Correlation_NQCD2b_vs_nFitSig", "Correlation_NQCD2b_vs_nFitSig", 800,600);
  RooPlot *frameEll14 = new RooPlot(*nFitQCD2b,*nFitSig,nFitQCD2b->getVal()-1.5*nFitQCD2b->getError(),nFitQCD2b->getVal()+1.5*nFitQCD2b->getError(),nFitSig->getVal()-1.5*nFitSig->getError(),nFitSig->getVal()+1.5*nFitSig->getError());
  res->plotOn(frameEll14,*nFitQCD2b,*nFitSig,"ME12ABHV");
  frameEll14->GetXaxis()->SetTitle("QCD events (2b category)");
  frameEll14->GetYaxis()->SetTitle("ttbar events (all categories)");
  frameEll14->Draw();
  /*
  //cout<<"correlation = "<<res->correlation(*nFitSig,*btagEff)<<endl;
  //cout<<"correlation = "<<res->correlation(*nFitQCD2b,*btagEff)<<endl;
  
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
  //-------------------- 0 bTag fit results -------------------------------
  
  RooPlot *frame0b = x->frame();
  combData.plotOn(frame0b,Cut("sample==sample::0btag")); 
  simPdf.plotOn(frame0b,Slice(sample,"0btag"),ProjWData(sample,combData));
  RooHist *pull0b = frame0b->pullHist();
  simPdf.plotOn(frame0b,Slice(sample,"0btag"),Components("qcd_pdf"),ProjWData(sample,combData),LineColor(kGreen),LineWidth(2),LineStyle(2));
  simPdf.plotOn(frame0b,Slice(sample,"0btag"),Components("ttbar_pdf_0btag"),ProjWData(sample,combData),LineColor(kRed),LineWidth(2),LineStyle(1));
  simPdf.plotOn(frame0b,Slice(sample,"0btag"),Components("bkg_pdf_0btag"),ProjWData(sample,combData),LineColor(kOrange+3),LineWidth(2),LineStyle(5)); 

  RooPlot *frame0bPull = x->frame();
  frame0bPull->addPlotable(pull0b,"p");

  TCanvas *can0b = new TCanvas("SimFit_0btag_"+CUT,"SimFit_0btag_"+CUT,900,600);
  can0b->cd(1)->SetBottomMargin(0.3);
  frame0b->GetXaxis()->SetTitle("");
  frame0b->GetXaxis()->SetLabelSize(0.0);
  frame0b->Draw();
  
  TPad *pad0b = new TPad("pad0b","pad0b",0.,0.,1.,1.);
  pad0b->SetTopMargin(0.7);
  pad0b->SetFillColor(0);
  pad0b->SetFillStyle(0);
  pad0b->Draw();
  pad0b->cd(0);
  pad0b->SetGridy();
  frame0bPull->SetMinimum(-5);
  frame0bPull->SetMaximum(5);
  frame0bPull->GetYaxis()->SetNdivisions(505);
  frame0bPull->GetXaxis()->SetTitleOffset(0.9);
  frame0bPull->GetYaxis()->SetTitleOffset(0.8);
  frame0bPull->GetYaxis()->SetTickLength(0.06);
  frame0bPull->GetYaxis()->SetTitleSize(0.05);
  frame0bPull->GetYaxis()->SetTitleSize(0.03);
  frame0bPull->GetYaxis()->SetLabelSize(0.03);
  frame0bPull->GetYaxis()->SetTitle("(Data-Fit)/Error");
  frame0bPull->GetXaxis()->SetTitle("m_{t} (GeV)");
  frame0bPull->Draw();
  

  //-------1 Btag Fit Results-----------------------------
  RooPlot *frame1b = x->frame();
  combData.plotOn(frame1b,Cut("sample==sample::1btag")); 
  simPdf.plotOn(frame1b,Slice(sample,"1btag"),ProjWData(sample,combData));
  RooHist *pull1b = frame1b->pullHist();
  simPdf.plotOn(frame1b,Slice(sample,"1btag"),Components("qcd_pdf"),ProjWData(sample,combData),LineColor(kGreen+1),LineWidth(2),LineStyle(2));
  simPdf.plotOn(frame1b,Slice(sample,"1btag"),Components("ttbar_pdf_1btag"),ProjWData(sample,combData),LineColor(kRed),LineWidth(2),LineStyle(1));
  simPdf.plotOn(frame1b,Slice(sample,"1btag"),Components("bkg_pdf_1btag"),ProjWData(sample,combData),LineColor(kOrange+3),LineWidth(2),LineStyle(5)); 

  RooPlot *frame1bPull = x->frame();
  frame1bPull->addPlotable(pull1b,"p");

  TCanvas *can1b = new TCanvas("SimFit_1btag_"+CUT,"SimFit_1btag_"+CUT,900,600);
  can1b->cd(1)->SetBottomMargin(0.3);
  frame1b->GetXaxis()->SetTitle("");
  frame1b->GetXaxis()->SetLabelSize(0.0);
  frame1b->Draw();
  
  TPad *pad1b = new TPad("pad1b","pad1b",0.,0.,1.,1.);
  pad1b->SetTopMargin(0.7);
  pad1b->SetFillColor(0);
  pad1b->SetFillStyle(0);
  pad1b->Draw();
  pad1b->cd(0);
  pad1b->SetGridy();
  frame1bPull->SetMinimum(-5);
  frame1bPull->SetMaximum(5);
  frame1bPull->GetYaxis()->SetNdivisions(505);
  frame1bPull->GetXaxis()->SetTitleOffset(0.9);
  frame1bPull->GetYaxis()->SetTitleOffset(0.8);
  frame1bPull->GetYaxis()->SetTickLength(0.06);
  frame1bPull->GetYaxis()->SetTitleSize(0.05);
  frame1bPull->GetYaxis()->SetTitleSize(0.03);
  frame1bPull->GetYaxis()->SetLabelSize(0.03);
  frame1bPull->GetYaxis()->SetTitle("(Data-Fit)/Error");
  frame1bPull->GetXaxis()->SetTitle("m_{t} (GeV)");
  frame1bPull->Draw();

  //-------------------- 2 bTag fit results -------------------------------
  RooPlot *frame2b = x->frame();
  combData.plotOn(frame2b,Cut("sample==sample::2btag")); 
  simPdf.plotOn(frame2b,Slice(sample,"2btag"),ProjWData(sample,combData));
  RooHist *pull2b = frame2b->pullHist();
  simPdf.plotOn(frame2b,Slice(sample,"2btag"),Components("qcd_pdf"),ProjWData(sample,combData),LineColor(kGreen+1),LineWidth(2),LineStyle(2));
  simPdf.plotOn(frame2b,Slice(sample,"2btag"),Components("ttbar_pdf_2btag"),ProjWData(sample,combData),LineColor(kRed),LineWidth(2),LineStyle(1));
  simPdf.plotOn(frame2b,Slice(sample,"2btag"),Components("bkg_pdf_2btag"),ProjWData(sample,combData),LineColor(kOrange+3),LineWidth(2),LineStyle(5)); 

  RooPlot *frame2bPull = x->frame();
  frame2bPull->addPlotable(pull2b,"p");

  TCanvas *can2b = new TCanvas("SimFit_2btag_"+CUT,"SimFit_2btag_"+CUT,900,600);
  can2b->cd(1)->SetBottomMargin(0.3);
  frame2b->GetXaxis()->SetTitle("");
  frame2b->GetXaxis()->SetLabelSize(0.0);
  frame2b->Draw();
  
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
  frame2bPull->GetYaxis()->SetTitleOffset(0.8);
  frame2bPull->GetYaxis()->SetTickLength(0.06);
  frame2bPull->GetYaxis()->SetTitleSize(0.05);
  frame2bPull->GetYaxis()->SetTitleSize(0.03);
  frame2bPull->GetYaxis()->SetLabelSize(0.03);
  frame2bPull->GetYaxis()->SetTitle("(Data-Fit)/Error");
  frame2bPull->GetXaxis()->SetTitle("m_{t} (GeV)");
  frame2bPull->Draw();
  

  
/*

  //other processes
  can0b->Print(TString(can0b->GetName())+".pdf");
  can2b->Print(TString(can2b->GetName())+".pdf");
  canEll1->Print(TString(canEll1->GetName())+".pdf");
  canEll2->Print(TString(canEll2->GetName())+".pdf");
  canEll3->Print(TString(canEll3->GetName())+".pdf");
  canEll4->Print(TString(canEll4->GetName())+".pdf");

  RooWorkspace *wOut = new RooWorkspace("w","workspace");
  wOut->import(*nFitQCD0b);
  wOut->import(*nFitQCD2b);
  wOut->import(*nFitSig);
  wOut->import(*btagEff);
  wOut->writeToFile("SimFitResults_"+CUT+".root");

  
    model->plotOn(frame);
    
    float chi2 = frame->chiSquare(2);
    float xsec = fsig->getVal()*nFitSig->getVal()/LUMI;
    float xsec_err = fsig->getVal()*nFitSig->getError()/LUMI;
    float xsec_exp = fsig->getVal()*exp_sig->getVal()/LUMI;
    float xsec_exp_err = fsig->getVal()*exp_sig->getError()/LUMI;

    cout<<"Signal          = "<<nFitSig->getVal()<<" +/- "<<nFitSig->getError()<<endl;
    cout<<"QCD             = "<<nFitQCD->getVal()<<" +/- "<<nFitQCD->getError()<<endl;
    cout<<"Bkg             = "<<nFitBkg->getVal()<<" +/- "<<nFitBkg->getError()<<endl;
    cout<<"Exp. Signal     = "<<exp_sig->getVal()<<" +/- "<<exp_sig->getError()<<endl;
    cout<<"r               = "<<nFitSig->getVal()/exp_sig->getVal()<<endl;
    cout<<"chi2/ndof       = "<<chi2<<endl;
    cout<<"Signal fraction = "<<fsig->getVal()<<endl;
    cout<<"QCD fraction    = "<<fqcd->getVal()<<endl;
    cout<<"correlation     = "<<res->correlation(*nFitSig,*nFitQCD)<<endl;
    
    hYieldSig->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
    hYieldSig->SetBinContent(k+1,nFitSig->getVal());
    hYieldSig->SetBinError(k+1,nFitSig->getError());
    hYieldBkg->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
    hYieldBkg->SetBinContent(k+1,nFitQCD->getVal());
    hYieldBkg->SetBinError(k+1,nFitQCD->getError());
    hYieldCor->SetBinContent(k+1,res->correlation(*nFitSig,*nFitQCD));
    
    hYieldSigExp->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
    hYieldSigExp->SetBinContent(k+1,exp_sig->getVal());
    hYieldSigExp->SetBinError(k+1,exp_sig->getError());

    hAccSigExp->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
    hAccSigExp->SetBinContent(k+1,exp_acc->getVal());
    hAccSigExp->SetBinError(k+1,exp_acc->getError());

    hXsecFid->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
    hXsecFid->SetBinContent(k+1,xsec);
    hXsecFid->SetBinError(k+1,xsec_err);

    hXsecFidExp->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
    hXsecFidExp->SetBinContent(k+1,xsec_exp);
    hXsecFidExp->SetBinError(k+1,xsec_exp_err);
    
    hXsecExtr->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
    hXsecExtr->SetBinContent(k+1,xsec/(fsig->getVal()*exp_acc->getVal()));
    hXsecExtr->SetBinError(k+1,(xsec/(fsig->getVal()*exp_acc->getVal()))*sqrt(pow(xsec_err/xsec,2)+pow(exp_acc->getError()/exp_acc->getVal(),2)));
    
    hChi2->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
    hChi2->SetBinContent(k+1,chi2);
    
    hSigFrac->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
    hSigFrac->SetBinContent(k+1,fsig->getVal());
    
    hBkgFrac->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
    hBkgFrac->SetBinContent(k+1,fqcd->getVal());
    
    RooHist *hpull = frame->pullHist();
    //model->plotOn(frame,RooFit::VisualizeError(*res,1,kFALSE),RooFit::FillColor(kGray),RooFit::MoveToBack());
    model->plotOn(frame,RooFit::Components("qcdCor_pdf"),RooFit::LineColor(kGreen+1),RooFit::LineWidth(2),RooFit::LineStyle(2));
    model->plotOn(frame,RooFit::Components("bkg_pdf"),RooFit::LineColor(kBlack),RooFit::LineWidth(2),RooFit::LineStyle(4));
    model->plotOn(frame,RooFit::Components("ttbar_pdf_"+TTBARNAME),RooFit::LineColor(kRed),RooFit::LineWidth(2),RooFit::LineStyle(1));

    can[k] = new TCanvas("Fitted_"+CUT+"_"+CAT+"_mTop_"+ALIAS[k],"Fitted_"+CUT+"_"+CAT+"_mTop_"+ALIAS[k],900,600);
    can[k]->cd(1)->SetBottomMargin(0.3);
    frame->SetMinimum(0.1);
    frame->GetYaxis()->SetNdivisions(505);
    frame->GetXaxis()->SetTitle("");
    frame->GetXaxis()->SetLabelSize(0.0);
    frame->Draw(); 

    RooHist  *hist_data = (RooHist*)frame->findObject("h_roohist_data",RooHist::Class());
    RooCurve *curve_all = (RooCurve*)frame->findObject("model_Norm[mTop]",RooCurve::Class());
    RooCurve *curve_sig = (RooCurve*)frame->findObject("model_Norm[mTop]_Comp[ttbar_pdf_"+TTBARNAME+"]",RooCurve::Class());
    RooCurve *curve_qcd = (RooCurve*)frame->findObject("model_Norm[mTop]_Comp[qcdCor_pdf]",RooCurve::Class());
    RooCurve *curve_bkg = (RooCurve*)frame->findObject("model_Norm[mTop]_Comp[bkg_pdf]",RooCurve::Class());

    TLegend *leg = new TLegend(0.75,0.7,0.92,0.9);
    leg->SetHeader(ALIAS[k]);
    leg->AddEntry(hist_data,"data","P");
    leg->AddEntry(curve_all,"total","L");
    leg->AddEntry(curve_sig,"ttbar","L");
    leg->AddEntry(curve_qcd,"qcd","L");
    leg->AddEntry(curve_bkg,"other bkg","L"); 
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    leg->Draw();
        
    RooPlot *frame2 = x->frame();
    frame2->addPlotable(hpull,"p");
  
    TPad *pad = new TPad("pad","pad",0.,0.,1.,1.);
    pad->SetTopMargin(0.7);
    //pad->SetGridy();
    pad->SetFillColor(0);
    pad->SetFillStyle(0);
    pad->Draw();
    pad->cd(0);
    pad->SetGridy();
    frame2->SetMinimum(-5);
    frame2->SetMaximum(5);
    frame2->GetYaxis()->SetNdivisions(505);
    frame2->GetXaxis()->SetTitleOffset(0.9);
    frame2->GetYaxis()->SetTitleOffset(0.8);
    frame2->GetYaxis()->SetTickLength(0.06);
    frame2->GetYaxis()->SetTitleSize(0.05);
    frame2->GetYaxis()->SetTitleSize(0.03);
    frame2->GetYaxis()->SetLabelSize(0.03);
    frame2->GetYaxis()->SetTitle("(Data-Fit)/Error");
    frame2->GetXaxis()->SetTitle("m_{t} (GeV)");
    frame2->Draw();   

    can[k]->Print("plots/"+TString(can[k]->GetName())+".pdf");
    can[k]->Print("plots/"+TString(can[k]->GetName())+".png");
    
    RooAbsReal *nll = model->createNLL(*roohist_data,RooFit::NumCPU(2));
    RooMinuit(*nll).migrad();
    RooAbsReal *pll_sig = nll->createProfile(*nFitSig);
    //RooAbsReal *pll_bkg = nll->createProfile(RooArgSet(*nFitSig,kBkg));
    RooPlot *frame_nll = nFitSig->frame(RooFit::Bins(10),RooFit::Range(nFitSig->getVal()-2*nFitSig->getError(),nFitSig->getVal()+2*nFitSig->getError()));
    pll_sig->plotOn(frame_nll,RooFit::ShiftToZero());
    //pll_bkg->plotOn(frame_nll,RooFit::LineStyle(3),RooFit::ShiftToZero(),RooFit::ShiftToZero());
    nll->plotOn(frame_nll,RooFit::LineStyle(kDashed),RooFit::ShiftToZero());
    TCanvas *c = new TCanvas("c","c",900,600);
    frame_nll->SetMinimum(0.0);
    frame_nll->SetMaximum(3.0); 
    frame_nll->Draw();
    */


  /* 
  Withing SRA
  TH2

  */
}

