void deepAK8_yields_jetPt()
{
  BoostedAnalysisInitializer utils("2018", true);
  std::cout<<utils.LUMI<<std::endl;
  
  std::vector<float> workingPoints = {0.6/*, 0.5, 0.4, 0.3, 0.2*/};
  
  TFile *sigFile = TFile::Open("2018/rootFiles/deepAK8_efficienciesJetPtOnlyRecoCuts.root");
  TFile *bkgFile = TFile::Open("2018/rootFiles/deepAK8_efficienciesJetPtOnlyRecoCutsBkg.root");
  
  TFile *yieldFile = TFile::Open("deepAK8_yields_jetPt.root", "UPDATE");
  for(std::vector<float>::iterator workingPoint=workingPoints.begin(); workingPoint!=workingPoints.end(); ++workingPoint)
  {
    TEfficiency *sig_Eff_deepAK8 = (TEfficiency*) sigFile->Get(TString::Format("Sig_Reco_jetPt_deepAK8_%.1f", *workingPoint));
    TEfficiency *bkg_Eff_deepAK8 = (TEfficiency*) bkgFile->Get(TString::Format("Bkg_Reco_jetPt_deepAK8_%.1f", *workingPoint));
    
    yieldFile->cd();
    TH1F *sig_yield_deepAK8 = (TH1F*) sig_Eff_deepAK8->GetPassedHistogram();
    sig_yield_deepAK8->Scale(utils.LUMI);
    sig_yield_deepAK8->Write(TString::Format("Yield_Reco_jetPt_deepAK8_%.1f", *workingPoint));
    TH1F *bkg_yield_deepAK8 = (TH1F*) bkg_Eff_deepAK8->GetPassedHistogram();
    bkg_yield_deepAK8->Scale(utils.LUMI);
    bkg_yield_deepAK8->Write(TString::Format("Bkg_Reco_jetPt_deepAK8_%.1f", *workingPoint));
  }
  
  TEfficiency *sig_Eff_tTagger = (TEfficiency*) sigFile->Get("Sig_Reco_jetPt_tTagger_0.2");
  TEfficiency *bkg_Eff_tTagger = (TEfficiency*) bkgFile->Get("Bkg_Reco_jetPt_tTagger_0.2");
  
  TH1F *sig_yield_tTagger = (TH1F*) sig_Eff_tTagger->GetPassedHistogram();
  sig_yield_tTagger->Scale(utils.LUMI);
  sig_yield_tTagger->Write("Yield_Reco_jetPt_tTagger_0.2");
  TH1F *bkg_yield_tTagger = (TH1F*) bkg_Eff_tTagger->GetPassedHistogram();
  bkg_yield_tTagger->Scale(utils.LUMI);
  bkg_yield_tTagger->Write("Bkg_Reco_jetPt_tTagger_0.2");
  
  TEfficiency *sig_Eff_mva = (TEfficiency*) sigFile->Get("Sig_Reco_jetPt_eventMVA_0.8");
  TEfficiency *bkg_Eff_mva = (TEfficiency*) bkgFile->Get("Bkg_Reco_jetPt_eventMVA_0.8");
  
  TH1F *sig_yield_mva = (TH1F*) sig_Eff_mva->GetPassedHistogram();
  sig_yield_mva->Scale(utils.LUMI);
  sig_yield_mva->Write("Yield_Reco_jetPt_eventMVA_0.8");
  TH1F *bkg_yield_mva = (TH1F*) bkg_Eff_mva->GetPassedHistogram();
  bkg_yield_mva->Scale(utils.LUMI);
  bkg_yield_mva->Write("Bkg_Reco_jetPt_eventMVA_0.8");
}
