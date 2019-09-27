void deepAK8_plotYields(float workingPoint)
{
  gStyle->SetOptStat(0);
  TFile *f = TFile::Open("deepAK8_yields.root");
  TCanvas *c1 = new TCanvas("c1", "c1", 600, 500);
  
  TH1F *sig_deepAK8 = (TH1F*) f->Get(TString::Format("Yield_Reco_deepAK8_%.1f", workingPoint));
  sig_deepAK8->SetMarkerColor(kRed);
  sig_deepAK8->SetLineColor(kRed);
  sig_deepAK8->SetTitle(TString::Format("Yield deepAK8 (%.1f) vs tTagger", workingPoint));
  sig_deepAK8->GetYaxis()->SetTitle("Events");
  sig_deepAK8->Draw();
  
  TH1F *bkg_deepAK8 = (TH1F*) f->Get(TString::Format("Bkg_Reco_deepAK8_%.1f", workingPoint));
  bkg_deepAK8->SetMarkerColor(kBlack);
  bkg_deepAK8->SetLineColor(kBlack);
  bkg_deepAK8->Draw("SAME");
  
  TH1F *sig_tTagger = (TH1F*) f->Get("Yield_Reco_tTagger_0.3");
  sig_tTagger->SetMarkerColor(kGreen);
  sig_tTagger->SetLineColor(kGreen);
  sig_tTagger->Draw("SAME");
  
  TH1F *bkg_tTagger = (TH1F*) f->Get("Bkg_Reco_tTagger_0.3");
  bkg_tTagger->SetMarkerColor(kMagenta);
  bkg_tTagger->SetLineColor(kMagenta);
  bkg_tTagger->Draw("SAME");
  
  TLegend *leg = new TLegend(0.4, 0.4, 0.8, 0.8);
  leg->AddEntry(sig_deepAK8, "Signal deepAK8", "lp");
  leg->AddEntry(bkg_deepAK8, "Bkg deepAK8", "lp");
  leg->AddEntry(sig_tTagger, "Signal tTagger", "lp");
  leg->AddEntry(bkg_tTagger, "Bkg tTagger", "lp");
  leg->Draw();
}