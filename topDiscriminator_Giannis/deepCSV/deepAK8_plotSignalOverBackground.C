void deepAK8_plotSignalOverBackground()
{
  gStyle->SetOptStat(0);
  std::vector<float> workingPoints = {0.6/*, 0.5, 0.4, 0.3, 0.2*/};
  std::vector<Color_t> colors = {kRed, kGreen, kBlack, kCyan, kOrange};
  TFile *f = TFile::Open("2018/rootFiles/signal_over_background.root");
  
  TCanvas *reco_Canvas = new TCanvas("reco_Canvas", "reco_Canvas", 600, 500);
  reco_Canvas->SetLogy();
  
  TLegend *reco_legend = new TLegend(0.5, 0.5, 0.8, 0.8);
  
  reco_Canvas->cd();
  TH1F* h_reco_tTagger = (TH1F*) f->Get("signal_over_bkg_tTagger_0.2");
  h_reco_tTagger->SetMarkerColor(kBlue);
  h_reco_tTagger->SetMarkerStyle(20);
  h_reco_tTagger->SetLineColor(kBlue);
  h_reco_tTagger->SetTitle("Signal Over Background");
  h_reco_tTagger->Draw("P");
  reco_legend->AddEntry(h_reco_tTagger, "topTagger (0.2)", "lp");
  
  TH1F* h_reco_mva = (TH1F*) f->Get("signal_over_bkg_eventMVA_0.8");
  h_reco_mva->SetMarkerColor(kMagenta);
  h_reco_mva->SetMarkerStyle(21);
  h_reco_mva->SetLineColor(kMagenta);
  h_reco_mva->SetTitle("Signal Over Background");
  h_reco_mva->Draw("PSAME");
  reco_legend->AddEntry(h_reco_mva, "event mva (0.8)", "lp");
  
  for(int i=0; i<workingPoints.size(); i++)
  {
    reco_Canvas->cd();
    TH1F* h_reco = (TH1F*) f->Get(TString::Format("signal_over_bkg_deepAK8_%.1f", workingPoints[i]));
    h_reco->SetLineColor(colors[i]);
    h_reco->SetMarkerColor(colors[i]);
    h_reco->SetMarkerStyle(22+i);
    h_reco->Draw("PSAME");
    reco_legend->AddEntry(h_reco, TString::Format("deepAK8 (%.1f)", workingPoints[i]), "lp");
  }
  
  reco_Canvas->cd();
  reco_legend->Draw();
}