void sigOverBkg()
{
    TFile *signalFile = TFile::Open("2016/rootFiles/deepAK8_efficienciesOnlyRecoCuts.root");
    TFile *bkgFile = TFile::Open("2016/rootFiles/deepAK8_efficienciesOnlyRecoCutsBkg.root");
    
    std::vector<float> workingPoints = {0.6/*, 0.5, 0.4, 0.3, 0.2*/};
    TFile *outFile = TFile::Open("2016/rootFiles/signal_over_background.root", "UPDATE");
    
    for(std::vector<float>::iterator workingPoint=workingPoints.begin(); workingPoint!=workingPoints.end(); ++workingPoint)
    {
      TEfficiency* eff_sig_deepAK8 = (TEfficiency*) signalFile->Get(TString::Format("Sig_Reco_deepAK8_%.1f", *workingPoint));
      TEfficiency* eff_bkg_deepAK8 = (TEfficiency*) bkgFile->Get(TString::Format("Bkg_Reco_deepAK8_%.2f", *workingPoint));
      
      TH1F* events_sig_deepAK8 = (TH1F*) eff_sig_deepAK8->GetPassedHistogram();
      TH1F* events_bkg_deepAK8 = (TH1F*) eff_bkg_deepAK8->GetPassedHistogram();
      events_sig_deepAK8->Divide(events_bkg_deepAK8);
      
      events_sig_deepAK8->Write(TString::Format("signal_over_bkg_deepAK8_%.1f", *workingPoint));
    }
    
    TEfficiency* eff_sig_tTagger = (TEfficiency*) signalFile->Get("Sig_Reco_tTagger_0.20");
    TEfficiency* eff_bkg_tTagger = (TEfficiency*) bkgFile->Get("Bkg_Reco_tTagger_0.20");
    
    TH1F* events_sig_tTagger = (TH1F*) eff_sig_tTagger->GetPassedHistogram();
    TH1F* events_bkg_tTagger = (TH1F*) eff_bkg_tTagger->GetPassedHistogram();
    
    events_sig_tTagger->Divide(events_bkg_tTagger);
    events_sig_tTagger->Write("signal_over_bkg_tTagger_0.2");
    
    TEfficiency* eff_sig_mva = (TEfficiency*) signalFile->Get("Sig_Reco_eventMVA_0.8");
    TEfficiency* eff_bkg_mva = (TEfficiency*) bkgFile->Get("Bkg_Reco_eventMVA_0.80");
    
    TH1F* events_sig_mva = (TH1F*) eff_sig_mva->GetPassedHistogram();
    TH1F* events_bkg_mva = (TH1F*) eff_bkg_mva->GetPassedHistogram();
    
    events_sig_mva->Divide(events_bkg_mva);
    events_sig_mva->Write("signal_over_bkg_eventMVA_0.8");
    
    outFile->Close();
}
