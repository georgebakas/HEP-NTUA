#include "TemplateConstants.h"

void SignalStrengthBefore_AfterBTagSFComparison()
{
  initFilesMapping();

  float sigStrength[3] = {
    ttbarSigStrength["2016"],
    ttbarSigStrength["2017"],
    ttbarSigStrength["2018"]
  };
  //get errors
  float sigStrengthError[3] = {
    ttbarSigStrengthError["2016"],
    ttbarSigStrengthError["2017"],
    ttbarSigStrengthError["2018"]
  };

  float sigStrength_noBTagSF[3] = {
    ttbarSigStrength_noBTagSF["2016"],
    ttbarSigStrength_noBTagSF["2017"],
    ttbarSigStrength_noBTagSF["2018"]
  };

  float sigStrengthError_noBTagSF[3] = {
    ttbarSigStrengthError_noBTagSF["2016"],
    ttbarSigStrengthError_noBTagSF["2017"],
    ttbarSigStrengthError_noBTagSF["2018"]
  };

  //ttbarSigStrengthError
  const int n = 3;
  float e_x[n] = {0,0,0};
  float x[n] = {0.5,1.5,2.5};
  float x_noBTag[n] = {0.7,1.7,2.7};

  TGraphErrors *gr = new TGraphErrors(n,x, sigStrength, e_x, sigStrengthError);
  TGraphErrors *gr_noBTagSF = new TGraphErrors(n,x_noBTag, sigStrength_noBTagSF,e_x, sigStrengthError_noBTagSF);

  gr->SetMarkerColor(kRed);
  gr->SetMarkerStyle(22);
  gr_noBTagSF->SetMarkerColor(kBlue);
  gr_noBTagSF->SetMarkerStyle(20);

  //use TMultiGraph to draw the two graphs together
  TCanvas *can = new TCanvas("can", "can", 800, 600);
  TLegend *leg = new TLegend(0.4,0.75,0.6,0.9);
  leg->AddEntry(gr, "bTag SFs", "lep");
  leg->AddEntry(gr_noBTagSF, "NO bTag SFs", "lep");

  TString Xname[n] = {"2016","2017","2018"};
  //gr->GetHistogram()->GetXaxis()->Set(3, 0, 3);
  //gr_noBTagSF->GetHistogram()->GetXaxis()->Set(3, 0, 3);
  //for(int i=0; i<n; i++){
  //  gr->GetHistogram()->GetXaxis()->SetBinLabel(i+1, Xname[i]);
  //  gr_noBTagSF->GetHistogram()->GetXaxis()->SetBinLabel(i+1, Xname[i]);
  //}


  TMultiGraph  *mg  = new TMultiGraph();
  mg->Add(gr);
  mg->Add(gr_noBTagSF);
  mg->GetHistogram()->GetXaxis()->Set(3, 0, 3);
  for(int i=0; i<n; i++){
    mg->GetHistogram()->GetXaxis()->SetBinLabel(i+1, Xname[i]);
  }
  mg->GetYaxis()->SetTitle("t#bar{t} Sig. Strength");
  mg->GetYaxis()->SetRangeUser(0.5, 0.8);

  mg->Draw("AP");
  leg->Draw();

  can->Print("SignalStrengthBefore_AfterBTagSFComparison.pdf", "pdf");

}
