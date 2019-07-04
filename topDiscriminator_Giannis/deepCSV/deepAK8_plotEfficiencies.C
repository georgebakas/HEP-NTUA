void deepAK8_plotEfficiencies()
{
  std::vector<float> workingPoints = {0.6, 0.5, 0.4, 0.3, 0.2};
  std::vector<Color_t> colors = {kRed, kGreen, kBlack, kMagenta, kOrange};
  TFile *f = TFile::Open("deepAK8_efficiencies.root");
  
  TCanvas *reco_Canvas = new TCanvas("reco_Canvas", "reco_Canvas", 600, 500);
  TCanvas *parton_Canvas = new TCanvas("parton_Canvas", "parton_Canvas", 600, 500);
  
  TLegend *reco_legend = new TLegend(0.5, 0.5, 0.8, 0.8);
  TLegend *parton_legend = new TLegend(0.5, 0.5, 0.8, 0.8);
  
  reco_Canvas->cd();
  TH1F* h_reco_tTagger = (TH1F*) f->Get("Sig_Reco_tTagger_0.3");
  h_reco_tTagger->SetMarkerColor(kBlue);
  h_reco_tTagger->SetLineColor(kBlue);
  h_reco_tTagger->SetTitle("Efficiencies Reco Level");
  h_reco_tTagger->Draw();
  //h_reco_tTagger->GetYaxis()->SetRangeUser(0, 0.15);
  reco_legend->AddEntry(h_reco_tTagger, "topTagger (0.3)", "l");
  
  parton_Canvas->cd();
  TH1F* h_parton_tTagger = (TH1F*) f->Get("Sig_Parton_tTagger_0.3");
  h_parton_tTagger->SetMarkerColor(kBlue);
  h_parton_tTagger->SetLineColor(kBlue);
  h_parton_tTagger->SetTitle("Efficiencies Parton Level");
  h_parton_tTagger->Draw();
  //h_parton_tTagger->GetYaxis()->SetRangeUser(0, 0.15);
  parton_legend->AddEntry(h_parton_tTagger, "topTagger (0.3)", "l");
  
  for(int i=0; i<workingPoints.size(); i++)
  {
    reco_Canvas->cd();
    TH1F* h_reco = (TH1F*) f->Get(TString::Format("Sig_Reco_deepAK8_%.1f", workingPoints[i]));
    h_reco->SetLineColor(colors[i]);
    h_reco->SetMarkerColor(colors[i]);
    h_reco->Draw("SAME");
    reco_legend->AddEntry(h_reco, TString::Format("deepAK8 (%.1f)", workingPoints[i]), "l");
    
    parton_Canvas->cd();
    TH1F* h_parton = (TH1F*) f->Get(TString::Format("Sig_Parton_deepAK8_%.1f", workingPoints[i]));
    h_parton->SetLineColor(colors[i]);
    h_parton->SetMarkerColor(colors[i]);
    h_parton->Draw("SAME");
    parton_legend->AddEntry(h_parton, TString::Format("deepAK8 (%.1f)", workingPoints[i]), "l");
  }
  
  reco_Canvas->cd();
  reco_legend->Draw();
  
  parton_Canvas->cd();
  parton_legend->Draw();
  
}
