void deepAK8_plotSignalOverBackground()
{
  gStyle->SetOptStat(0);
  std::vector<float> workingPoints = {0.6/*, 0.5, 0.4, 0.3, 0.2*/};
  std::vector<Color_t> colors = {kRed, kGreen, kBlack, kCyan, kOrange};
  TFile *f1 = TFile::Open("2016/rootFiles/signal_over_background.root");
  TFile *f = TFile::Open("2017/rootFiles/signal_over_background.root");

  
  TCanvas *reco_Canvas = new TCanvas("reco_Canvas", "reco_Canvas", 600, 500);
  reco_Canvas->SetLogy();
  
  TLegend *reco_legend = new TLegend(0.5, 0.5, 0.8, 0.8);
  
  reco_Canvas->cd();
  TH1F* h_reco_tTagger = (TH1F*) f->Get("signal_over_bkg_tTagger_0.0");
  h_reco_tTagger->SetMarkerColor(kBlue);
  h_reco_tTagger->SetMarkerStyle(20);
  h_reco_tTagger->SetLineColor(kBlue);
  h_reco_tTagger->SetTitle("Signal Over Background");
  h_reco_tTagger->Draw("P");
  reco_legend->AddEntry(h_reco_tTagger, "topTagger '17 (0.0)", "lp");
  
  
  TH1F* h_reco_tTagger_2 = (TH1F*) f->Get("signal_over_bkg_tTagger_0.1");
  h_reco_tTagger_2->SetMarkerColor(kBlack);
  h_reco_tTagger_2->SetMarkerStyle(21);
  h_reco_tTagger_2->SetLineColor(kBlack);
  h_reco_tTagger_2->SetTitle("Signal Over Background");
  h_reco_tTagger_2->Draw("PSAME");
  reco_legend->AddEntry(h_reco_tTagger_2, "topTagger '17 (0.1)", "lp");
  
    
  TH1F* h_reco_tTagger_16 = (TH1F*) f1->Get("signal_over_bkg_tTagger_0.2");
  h_reco_tTagger_16->SetMarkerColor(kRed);
  h_reco_tTagger_16->SetMarkerStyle(22);
  h_reco_tTagger_16->SetLineColor(kRed);
  h_reco_tTagger_16->SetTitle("Signal Over Background");
  h_reco_tTagger_16->Draw("PSAME");
  reco_legend->AddEntry(h_reco_tTagger_16, "topTagger '16 (0.2)", "lp");
  
  /*
  for(int i=0; i<workingPoints.size(); i++)
  {
    reco_Canvas->cd();
    TH1F* h_reco = (TH1F*) f->Get(TString::Format("signal_over_bkg_deepAK8_%.1f", workingPoints[i]));
    h_reco->SetLineColor(colors[i]);
    h_reco->SetMarkerColor(colors[i]);
    h_reco->SetMarkerStyle(22+i);
    h_reco->Draw("PSAME");
    reco_legend->AddEntry(h_reco, TString::Format("deepAK8 (%.1f)", workingPoints[i]), "lp");
  } */
  
  reco_Canvas->cd();
  reco_legend->Draw();
}