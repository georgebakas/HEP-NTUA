void deepAK8_plotYields_jetPt(float workingPoint)
{
  gStyle->SetOptStat(0);
  TFile *f = TFile::Open("2018/rootFiles/deepAK8_yields_jetPt.root");
  TCanvas *c1 = new TCanvas("c1", "c1", 600, 500);
  
  TH1F *sig_deepAK8 = (TH1F*) f->Get(TString::Format("Yield_Reco_jetPt_deepAK8_%.1f", workingPoint));
  sig_deepAK8->SetMarkerColor(kRed);
  sig_deepAK8->SetMarkerStyle(20);
  sig_deepAK8->SetLineColor(kRed);
  //sig_deepAK8->SetTitle(TString::Format("Yield deepAK8 (%.1f) vs tTagger", workingPoint));
  sig_deepAK8->SetTitle("Yields");
  sig_deepAK8->GetYaxis()->SetTitle("Events");
  sig_deepAK8->Draw("P");
  
  TH1F *bkg_deepAK8 = (TH1F*) f->Get(TString::Format("Bkg_Reco_jetPt_deepAK8_%.1f", workingPoint));
  bkg_deepAK8->SetMarkerColor(kBlack);
  bkg_deepAK8->SetMarkerStyle(21);
  bkg_deepAK8->SetLineColor(kBlack);
  bkg_deepAK8->Draw("PSAME");
  
  TH1F *sig_tTagger = (TH1F*) f->Get("Yield_Reco_jetPt_tTagger_0.2");
  sig_tTagger->SetMarkerColor(kGreen);
  sig_tTagger->SetMarkerStyle(22);
  sig_tTagger->SetLineColor(kGreen);
  sig_tTagger->Draw("PSAME");
  
  TH1F *bkg_tTagger = (TH1F*) f->Get("Bkg_Reco_jetPt_tTagger_0.2");
  bkg_tTagger->SetMarkerColor(kMagenta);
  bkg_tTagger->SetMarkerStyle(23);
  bkg_tTagger->SetLineColor(kMagenta);
  bkg_tTagger->Draw("PSAME");
  
  TH1F *sig_mva = (TH1F*) f->Get("Yield_Reco_jetPt_eventMVA_0.8");
  sig_mva->SetMarkerColor(kOrange);
  sig_mva->SetMarkerStyle(24);
  sig_mva->SetLineColor(kOrange);
  sig_mva->Draw("PSAME");
  
  TH1F *bkg_mva = (TH1F*) f->Get("Bkg_Reco_jetPt_eventMVA_0.8");
  bkg_mva->SetMarkerColor(kBlue);
  bkg_mva->SetMarkerStyle(25);
  bkg_mva->SetLineColor(kBlue);
  bkg_mva->Draw("PSAME");
  
  TLegend *leg = new TLegend(0.4, 0.4, 0.8, 0.8);
  leg->AddEntry(sig_deepAK8, "Signal deepAK8", "lp");
  leg->AddEntry(bkg_deepAK8, "Bkg deepAK8", "lp");
  leg->AddEntry(sig_tTagger, "Signal tTagger", "lp");
  leg->AddEntry(bkg_tTagger, "Bkg tTagger", "lp");
  leg->AddEntry(sig_mva, "Signal event mva", "lp");
  leg->AddEntry(bkg_mva, "Bkg event mva", "lp");
  leg->Draw();
}