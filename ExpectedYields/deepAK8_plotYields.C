void deepAK8_plotYields()
{
  gStyle->SetOptStat(0);
  TFile *f16 = TFile::Open("2016/rootFiles/deepAK8_yields.root");
  TFile *f = TFile::Open("2018/rootFiles/deepAK8_yields.root");
  TCanvas *c1 = new TCanvas("c1", "c1", 600, 500);
  
  TH1F *sig_tTagger_16 = (TH1F*)f16->Get("Yield_Reco_tTagger_0.20");
  sig_tTagger_16->SetMarkerColor(kGreen);
  sig_tTagger_16->SetMarkerStyle(20);
  sig_tTagger_16->SetLineColor(kGreen);
  sig_tTagger_16->GetYaxis()->SetTitle("Events");
  sig_tTagger_16->SetTitle("Yields mJJ");
  sig_tTagger_16->Draw("P");
  
  TH1F *bkg_tTagger_16 = (TH1F*) f16->Get("Bkg_Reco_tTagger_0.20");
  bkg_tTagger_16->SetMarkerColor(kOrange);
  bkg_tTagger_16->SetMarkerStyle(21);
  bkg_tTagger_16->SetLineColor(kOrange);
  bkg_tTagger_16->Draw("PSAME");
  
  TH1F *sig_tTagger = (TH1F*) f->Get("Yield_Reco_tTagger_0.1");
  sig_tTagger->SetMarkerColor(kBlack);
  sig_tTagger->SetMarkerStyle(22);
  sig_tTagger->SetLineColor(kBlack);
  sig_tTagger->Draw("PSAME");
  sig_tTagger->GetYaxis()->SetTitle("Events");
  sig_tTagger->SetTitle("Yields mJJ");
  
  TH1F *bkg_tTagger = (TH1F*) f->Get("Bkg_Reco_tTagger_0.1");
  bkg_tTagger->SetMarkerColor(kMagenta);
  bkg_tTagger->SetMarkerStyle(23);
  bkg_tTagger->SetLineColor(kMagenta);
  bkg_tTagger->Draw("PSAME");
  
  TH1F *sig_tTagger_2 = (TH1F*) f->Get("Yield_Reco_tTagger_0.15");
  sig_tTagger_2->SetMarkerColor(kRed);
  sig_tTagger_2->SetMarkerStyle(24);
  sig_tTagger_2->SetLineColor(kRed);
  sig_tTagger_2->Draw("PSAME");
  
  TH1F *bkg_tTagger_2 = (TH1F*) f->Get("Bkg_Reco_tTagger_0.15");
  bkg_tTagger_2->SetMarkerColor(kBlue);
  bkg_tTagger_2->SetMarkerStyle(25);
  bkg_tTagger_2->SetLineColor(kBlue);
  bkg_tTagger_2->Draw("PSAME");
  
  TLegend *leg = new TLegend(0.4, 0.4, 0.8, 0.8);
  leg->AddEntry(sig_tTagger, "Signal tTagger '18(0.1)", "lp");
  leg->AddEntry(sig_tTagger_2, "Signal tTagger '18 (0.15)", "lp");
  leg->AddEntry(sig_tTagger_16, "Signal tTagger '16(0.2)", "lp");
  leg->AddEntry(bkg_tTagger, "Bkg tTagger '18(0.1)", "lp");
  leg->AddEntry(bkg_tTagger_2, "Bkg tTagger '18(0.15)", "lp");
  leg->AddEntry(bkg_tTagger_16, "Bkg tTagger '16(0.2)", "lp");
  leg->Draw();
}