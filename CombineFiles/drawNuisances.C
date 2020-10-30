void ComputeFraction(TH1D *h,TH2D *COV);

void ComputeFraction(TH1D *h,TH2D *COV)
{
  //--- compute fraction ---------
  float sum(0.0),sum1(0.0);
  int k = h->GetNbinsX()-5;

  for(int i=0;i<h->GetNbinsX();i++) {
    float y = h->GetBinContent(i+1);
    sum += y;
    if (i >= k) {
      sum1 += y;
    }
  }
  float sumV(0.0),sumV1(0.0),sumV2(0.0);
  for(int i=0;i<h->GetNbinsX();i++) {
    for(int j=0;j<h->GetNbinsX();j++) {
      float v = COV->GetBinContent(i+1,j+1);
      //cout<<i<<" "<<j<<" "<<v<<endl;
      sumV += v;
      if (i >= k) {
        sumV1 += v;
        if (j >= k) {
          sumV2 += v;
        }
      }
    }
  }
  float f  = sum1/sum;
  float ef = (1./sum)*sqrt(sumV2-2*(sum1/sum)*sumV1+(sum1*sum1/(sum*sum))*sumV);
  cout<<"f = "<<f<<" +/- "<<ef<<", "<<ef/f<<endl;
}
//void drawNuisances(TString channel, int kbins, float XMIN, float XMAX, TString TITLE)
void drawNuisances(TString channel="", float XMIN=0, float XMAX=12)/*, TString channel, int kbins, float XMIN, float XMAX)*/
{
  channel = "SR_C";
  TString TITLE = "#chi";
  //TString bins = to_string(kbins);
  TString bins = "11";
  gROOT->ForceStyle();

  TFile *inf = TFile::Open("fitDiagnostics_fit_result.root", "READ");

  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);

  TH2D *hCorBin_p = (TH2D*)inf->Get("shapes_prefit/overall_total_covar");
  TH2D *hCorBin_s = (TH2D*)inf->Get("shapes_fit_s/overall_total_covar");
  TH1D *hSignal_p = (TH1D*)inf->Get("shapes_prefit/total_signal");
  TH1D *hSignal_s = (TH1D*)inf->Get("shapes_fit_s/total_signal");

  //ComputeFraction(hSignal_p,hCorBin_p);
  // ComputeFraction(hSignal_s,hCorBin_s);


  //hSignal_p->Sumw2();
  //hSignal_s->Sumw2();
  hSignal_p->Scale(1./hSignal_p->Integral());
  hSignal_s->Scale(1./hSignal_s->Integral());

  RooFitResult *fit_s  = (RooFitResult*)inf->Get("fit_s");
  RooFitResult *fit_b  = (RooFitResult*)inf->Get("fit_b");
  RooArgSet    *prefit = (RooArgSet*)inf->Get("nuisances_prefit");

  RooArgList fpf_b = fit_b->floatParsFinal();
  RooArgList fpf_s = fit_s->floatParsFinal();
  //RooArgList fpf_p = prefit->floatParsFinal();

  TH1F *hPulls_p = new TH1F("hPulls_p","hPulls_p",fpf_s.getSize()-1,0,fpf_s.getSize()-1);
  TH1F *hPulls_s = new TH1F("hPulls_s","hPulls_s",fpf_s.getSize()-1,0,fpf_s.getSize()-1);
  int N(0);
  //cout<<"fpf_s.getSize() = "<<fpf_s.getSize()<<endl;
  for(int i=0;i<fpf_s.getSize();i++) {
	  //cout<<"--i-"<<endl;
    RooAbsArg *nuis_s(fpf_s.at(i));
    TString name(nuis_s->GetName());
    RooRealVar *par_s = (RooRealVar*)fpf_s.find(name);
    //RooRealVar *par_b = (RooRealVar*)fpf_s.find(name); OLD
    RooRealVar *par_b = (RooRealVar*)fpf_b.find(name);
    RooRealVar *par_p = (RooRealVar*)prefit->find(name);
    //cout<<name<<" "<<par_s->getVal()<<" +/- "<<par_s->getError()<<endl;
	//cout<<"par_s : "<<par_s->getVal()<<endl;
	//cout<<"name = "<<name<<endl;

	//takes only nuis of datacard not the extra
     if (par_p) {
		if(!name.Contains("bin")){

      //cout<<name<<" "<<par_s->getVal()<<" +/- "<<par_s->getError()<<endl;
      //cout<<name<<" "<<par_p->getVal()<<" +/- "<<par_p->getError()<<endl;
      hPulls_s->SetBinContent(N+1,par_s->getVal());
      hPulls_s->SetBinError(N+1,par_s->getError());
      hPulls_p->SetBinContent(N+1,0);
      hPulls_p->SetBinError(N+1,1);
      hPulls_p->GetXaxis()->SetBinLabel(N+1,name);
      N++;
		}
    }
	//cout<<"----"<<endl;
  }

  TCanvas *can_signal = new TCanvas("SignalShape_"+channel,"SignalShape_"+channel,900,600);
  hSignal_p->SetFillColor(kRed-10);
  hSignal_p->SetMarkerColor(kRed);
  hSignal_p->SetMarkerStyle(21);
  hSignal_p->SetLineColor(kBlack);
  hSignal_s->SetMarkerColor(kBlue);
  hSignal_s->SetLineColor(kBlue);
  hSignal_s->SetLineWidth(3);
  hSignal_p->Draw("E2");
  hSignal_s->Draw("Esame");
  hSignal_p->GetYaxis()->SetTitle("Probability Density");
  hSignal_p->GetXaxis()->SetLabelOffset(999);
  hSignal_p->GetXaxis()->SetLabelSize(0);
  hSignal_p->GetXaxis()->SetTickSize(0);

  TLegend *leg = new TLegend(0.2,0.7,0.5,0.9);
  leg->AddEntry(hSignal_p,"Prefit","PF");
  leg->AddEntry(hSignal_s,"Postfit","PLE");
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->Draw();

  gPad->Update();

  TGaxis *axis = new TGaxis(gPad->GetUxmin(),gPad->GetUymin(),gPad->GetUxmax(),gPad->GetUymin(),XMIN,XMAX,510,"");
  axis->SetLabelFont(42);
  axis->SetLabelSize(0.05);
  axis->SetTitleSize(0.05);
  axis->SetTitle(TITLE);
  axis->SetName("axisPad");
  axis->Draw("same");

  TCanvas *can_pulls = new TCanvas("Pulls_"+channel,"Pulls_"+channel,900,600);
  hPulls_s->SetMarkerColor(kBlue);
  hPulls_s->SetMarkerStyle(20);
  hPulls_s->SetLineColor(kBlue);
  hPulls_s->SetLineWidth(3);
  hPulls_p->SetMarkerSize(0);
  hPulls_p->SetLineStyle(2);
  hPulls_p->SetFillColor(kRed-10);
  hPulls_p->GetYaxis()->SetRangeUser(-3,3);
  hPulls_p->GetYaxis()->SetTitle("Pull");
  hPulls_p->GetXaxis()->SetTitle("Nuisance parameter");
  hPulls_p->GetXaxis()->SetLabelFont(62);
  hPulls_p->Draw("histE2");
  hPulls_s->Draw("Esame");
}
