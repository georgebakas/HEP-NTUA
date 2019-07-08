void plotCombined()
{
  TFile *_file0 = TFile::Open("/afs/cern.ch/work/i/ipapakri/private/analysis/topDiscriminator/deepCSV/2016/rootFiles/signal_over_background.root");
  TH1F* sig_2016 = (TH1F*)_file0->Get("signal_over_bkg_tTagger_0.2");
  //TH1F* bkg_2016 = (TH1F*)_file0->Get("Bkg_Reco_tTagger_0.2");

  TFile *_file1 = TFile::Open("/afs/cern.ch/work/i/ipapakri/private/analysis/topDiscriminator/deepAK8/2016/rootFiles/signal_over_background.root");
  TH1F* sig_2017 = (TH1F*)_file1->Get("signal_over_bkg_tTagger_0.2"); 
  //TH1F* bkg_2017 = (TH1F*)_file1->Get("Bkg_Reco_tTagger_0.2");

  //TFile *_file2 = TFile::Open("2018/rootFiles/signal_over_background.root");
  //TH1F* sig_2018 = (TH1F*)_file2->Get("signal_over_bkg_tTagger_0.2");
  //TH1F* bkg_2018 = (TH1F*)_file2->Get("Bkg_Reco_tTagger_0.2");
  
  gStyle->SetOptStat(0);
  sig_2016->SetMarkerStyle(20);
  sig_2016->SetMarkerColor(kRed);
  sig_2016->SetLineColor(kRed);
  sig_2017->SetLineColor(kBlue);
  sig_2017->SetMarkerColor(kBlue);
  sig_2017->SetMarkerStyle(21);
  //sig_2018->SetMarkerStyle(22);
  //sig_2018->SetMarkerColor(kGreen);
  //sig_2018->SetLineColor(kGreen);

/*  bkg_2016->SetLineColor(kMagenta);
  bkg_2016->SetMarkerColor(kMagenta);
  bkg_2016->SetMarkerStyle(23);
  bkg_2017->SetMarkerStyle(24);
  bkg_2017->SetMarkerColor(kBlack);
  bkg_2017->SetLineColor(kBlack);
  bkg_2018->SetLineColor(kCyan);
  bkg_2018->SetMarkerColor(kCyan);
  bkg_2018->SetMarkerStyle(25);
*/
  sig_2016->SetTitle("Signal over Background");
  sig_2016->GetYaxis()->SetTitle("S/B");
  sig_2016->GetXaxis()->SetTitle("jetPt (GeV)");

  TLegend *leg = new TLegend(0.15, 0.5, 0.35, 0.7);
  leg->AddEntry(sig_2016, "deepCSV 2016", "pl");
  leg->AddEntry(sig_2017, "CSVv2 2016", "pl");
  //leg->AddEntry(sig_2018, "2018", "pl");
  //leg->AddEntry(bkg_2016, "bkg 2016", "pl");
  //leg->AddEntry(bkg_2017, "bkg 2017", "pl");
  //leg->AddEntry(bkg_2018, "bkg_2018", "pl");

  sig_2016->Draw();
  sig_2017->Draw("PSAME");
  //sig_2018->Draw("PSAME");
  //bkg_2016->Draw("PSAME");
  //bkg_2017->Draw("PSAME");
  //bkg_2018->Draw("PSAME");

  leg->Draw();
}
