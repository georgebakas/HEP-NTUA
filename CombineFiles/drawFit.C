void drawFit(TString year, TString fit, float XMIN=0, float XMAX=16)/*, TString channel, int kbins, float XMIN, float XMAX)*/
{
  //TString bins = to_string(kbins);
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);

  TFile *inf = TFile::Open(year+"/fitDiagnostics."+fit+".root", "READ");


  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);

  TGraphAsymmErrors *gData = (TGraphAsymmErrors*)inf->Get("shapes_fit_s/SR_C/data");
  TH1D *hOverall    = (TH1D*)inf->Get("shapes_fit_s/SR_C/total");
  TH1D *hSignal     = (TH1D*)inf->Get("shapes_fit_s/SR_C/total_signal");


  THStack *h_stack = new THStack("h_stack", "hstack");
  const int NBKG = 3;
  TString bkg[NBKG]={"ttbar", "Subdominant", "qcd"};
  TH1D* h[NBKG];
  for(int i=0; i<NBKG; i++)
  {
	  h[i]=(TH1D*)inf->Get("shapes_fit_s/SR_C/"+bkg[i]);
	  if(i==2)
	  {
		h[i]->SetFillColor(30);
		h[i]->SetLineColor(30);
	  }
	  else
	  {
		h[i]->SetFillColor(37+i);
		h[i]->SetLineColor(37+i);
	  }

	  h[i]->SetLineWidth(2);
	  h[i]->SetFillStyle(3004);
  }
  for(int i=NBKG-1; i>=0;i--)
	  h_stack->Add(h[i]);

  h[2]->SetFillColor(30);//qcd
  h[2]->SetLineColor(30);

  h[0]->SetLineColor(kBlue-3);//TTbar
  h[0]->SetFillColor(kBlue-3);
  h[1]->SetLineColor(kMagenta);//Subdominant
  h[1]->SetFillColor(kMagenta);


  TH1D *hOverallAux = (TH1D*)hOverall->Clone("aux");

  TGraphAsymmErrors *gResiduals = (TGraphAsymmErrors*)gData->Clone("Residuals");
  TH1D *hFitError = (TH1D*)hOverall->Clone("FitError");

  double x,y,ey,eyl,eyh,yFit,eFit;
  for(int i=0;i<gResiduals->GetN();i++) {
    gResiduals->GetPoint(i,x,y);
    eyl  = gResiduals->GetErrorYlow(i);
    eyh  = gResiduals->GetErrorYhigh(i);
    ey   = gResiduals->GetErrorY(i);
    yFit = hOverall->GetBinContent(i+1);
    eFit = hOverall->GetBinError(i+1);


	gResiduals->SetPoint(i,x,(y-yFit)/ey);
    //gResiduals->SetPointEYlow(i,sqrt(eyl*eyl+eFit*eFit));
    //gResiduals->SetPointEYhigh(i,sqrt(eyh*eyh+eFit*eFit));
    gResiduals->SetPointEYlow(i,eyl/ey);
    gResiduals->SetPointEYhigh(i,eyh/ey);
    hFitError->SetBinContent(i+1,0);
    hFitError->SetBinError(i+1,eFit/ey);
    //cout<<i<<" "<<x<<" "<<y<<" "<<yFit<<endl;
  }
  hFitError->SetMarkerSize(0);
  hFitError->SetLineStyle(3);
  hFitError->SetFillColor(kOrange);

  hOverall->SetFillColor(kOrange);
  hOverall->SetMarkerSize(0);
  hOverall->SetLineWidth(3);
  hOverall->SetLineColor(kBlue);
  hOverallAux->SetLineWidth(3);
  hOverallAux->SetLineColor(kBlue);

  gData->SetMarkerStyle(20);
  gData->SetMarkerColor(kBlack);
  gData->SetLineColor(kBlack);
  gData->SetLineWidth(2);
  gResiduals->SetLineWidth(2);

  hSignal->SetLineWidth(2);
  hSignal->SetLineColor(kRed);
  hSignal->SetFillColor(kRed);
  hSignal->SetFillStyle(3005);


  TCanvas *can_fit = new TCanvas("fit_"+fit+"_SR_C","fit_"+fit+"_SR_C",800,600);

  can_fit->SetBottomMargin(0.35);

  hOverall->GetXaxis()->SetLabelOffset(999);
  hOverall->GetXaxis()->SetLabelSize(0);
  hOverall->GetXaxis()->SetTickSize(0);
  hOverall->GetYaxis()->SetTitle("Events");
  hOverall->GetYaxis()->SetTitleFont(62);
  hOverall->SetMaximum(1.5*hOverall->GetBinContent(hOverall->GetMaximumBin()));
  hOverall->Draw("E2");
  hOverallAux->Draw("same hist");
  h_stack->Add(hSignal);
  h_stack->Draw("same hist");


  gData->Draw("same EPZ");

  TLegend *leg = new TLegend(0.75,0.75,0.95,0.9);
  leg->AddEntry(gData,"data","ELP");
  leg->AddEntry(hOverall,"total fit","FL");
  leg->AddEntry(hSignal,"signal","F");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  for(int i=0; i<NBKG; i++)
  {
	  leg->AddEntry(h[i], bkg[i], "F");
  }
  leg->Draw();

  gPad->Update();

  TGaxis *axis = new TGaxis(gPad->GetUxmin(),0,gPad->GetUxmax(),0,XMIN,XMAX,510,"");
  axis->SetLabelSize(0);
  axis->SetTitleSize(0);
  axis->SetName("axis");
  axis->Draw("same");

  TPad *pad = new TPad("pad","pad",0.,0.,1.,1.);
  pad->SetTopMargin(0.68);
  pad->SetFillColor(0);
  pad->SetFillStyle(0);
  pad->Draw();
  pad->cd(0);
  pad->SetGridy();

  hFitError->GetXaxis()->SetLabelOffset(999);
  hFitError->GetXaxis()->SetLabelSize(0);
  hFitError->GetXaxis()->SetTickSize(0);
  hFitError->GetYaxis()->SetRangeUser(-3,3);
  hFitError->GetYaxis()->SetLabelSize(0.03);
  hFitError->GetYaxis()->SetTitleOffset(0.99);
  hFitError->GetXaxis()->SetTitleOffset(0.99);
  hFitError->GetYaxis()->SetTitleFont(62);
  hFitError->GetYaxis()->SetTitle("Pulls");
  hFitError->GetXaxis()->SetTitle("#chi");
  hFitError->Draw("histE2");
  gResiduals->Draw("samePZ");

  gPad->Update();

  TGaxis *axisPad = new TGaxis(gPad->GetUxmin(),-3,gPad->GetUxmax(),-3,XMIN,XMAX,510,"");
  axisPad->SetLabelFont(42);
  axisPad->SetLabelSize(0.05);
  axisPad->SetTitleSize(0.05);

  axisPad->SetName("axisPad");
  axisPad->Draw("same");

  can_fit->Print(year+"/plots/postFit_"+fit+".pdf", "pdf");
}
