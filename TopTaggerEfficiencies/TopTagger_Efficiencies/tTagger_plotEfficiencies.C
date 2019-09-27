void tTagger_plotEfficiencies()
{
  std::vector<float> workingPoints = {0.1,0.2, 0.3, 0.4, 0.5, 0.6};
  std::vector<Color_t> colors = {kBlue, kBlack, kRed, kGreen, kCyan, kOrange};
  
  
  
  TLegend *reco_legend = new TLegend(0.5, 0.5, 0.8, 0.8);
  TLegend *parton_legend = new TLegend(0.5, 0.5, 0.8, 0.8);
  
  TString varReco[5]   = {"mJJ","jetPt", "ptJJ", "yJJ",  "jetEta"}; 
  TString varParton[5] = {"mTTbarParton", "partonPt", "ptTTbarParton", "yTTbarParton", "partonEta"}; 
  TCanvas *reco_Canvas[5], *parton_Canvas[5];
  TEfficiency *h_reco_tTagger[5], *h_reco_oldMva;
  TEfficiency *h_parton_tTagger[5],*h_parton_oldMva;
  int initMarker = 22;
  
  for(int ivar= 0; ivar< 2; ivar++)
  {
   reco_Canvas[ivar] = new TCanvas(TString::Format("reco_Canvas_%s",varReco[ivar].Data()),TString::Format("reco_Canvas_%s", varReco[ivar].Data()), 700, 600);
   parton_Canvas[ivar] = new TCanvas(TString::Format("parton_Canvas_%s",varParton[ivar].Data()), TString::Format("parton_Canvas_%s", varParton[ivar].Data()), 700, 600);
   
   
   for(int iwp = 0; iwp<workingPoints.size(); iwp++)
   {
   
    TFile *f = TFile::Open(TString::Format("deepAK8_efficiencies_allVars_tTagger_%0.1f.root", workingPoints[iwp]));
	parton_Canvas[ivar]->cd();
	if(iwp ==0)
	{
		h_parton_oldMva = (TEfficiency*) f->Get(TString::Format("Sig_Parton_oldMva_0.8_%s", varParton[ivar].Data()));
		h_parton_oldMva->SetMarkerColor(kMagenta);
		h_parton_oldMva->SetLineColor(kMagenta);
		h_parton_oldMva->SetMarkerStyle(21);
		h_parton_oldMva->SetTitle("Efficiency Comparison");
		h_parton_oldMva->Draw();
		if(ivar ==0) parton_legend->AddEntry(h_parton_oldMva, "event MVA (0.8)", "lp");
	}	
	h_parton_tTagger[iwp] = (TEfficiency*) f->Get(TString::Format("Sig_Parton_tTagger_%0.1f_%s",workingPoints[iwp], varParton[ivar].Data()));  
	h_parton_tTagger[iwp]->SetMarkerColor(colors[iwp]);
	h_parton_tTagger[iwp]->SetLineColor(colors[iwp]);
	h_parton_tTagger[iwp]->SetMarkerStyle(initMarker);
	if(ivar == 0 || ivar ==1 || ivar ==3) h_parton_tTagger[iwp]->SetTitle(TString::Format("Efficiency Parton Level;%s (GeV);Efficiency",varParton[ivar].Data()) );
	else h_parton_tTagger[iwp]->SetTitle(TString::Format("Efficiency Comparison;%s;Efficiency",varParton[ivar].Data()) );
	h_parton_tTagger[iwp]->Draw("same");
    if(ivar ==0) parton_legend->AddEntry(h_parton_tTagger[iwp],TString::Format("topTagger (%0.1f)",workingPoints[iwp]) , "lp");
	
	
	reco_Canvas[ivar]->cd();
    if(iwp ==0)
	{
		h_reco_oldMva = (TEfficiency*) f->Get(TString::Format("Sig_Reco_oldMva_0.8_%s",varReco[ivar].Data() ));
		h_reco_oldMva->SetMarkerColor(kMagenta);
		h_reco_oldMva->SetLineColor(kMagenta);
		h_reco_oldMva->SetMarkerStyle(21);
		if(ivar == 0 || ivar ==1 || ivar ==3) h_reco_oldMva->SetTitle(TString::Format("Acceptance Comparison;%s (GeV);Acceptance",varReco[ivar].Data()) );
		else h_reco_oldMva->SetTitle(TString::Format("Acceptance Comparison;%s;Acceptance",varReco[ivar].Data()) );
		h_reco_oldMva->Draw();
		if(ivar ==0) reco_legend->AddEntry(h_reco_oldMva, "event MVA (0.8)", "lp");
	}

	h_reco_tTagger[iwp] = (TEfficiency*)f->Get(TString::Format("Sig_Reco_tTagger_%0.1f_%s",workingPoints[iwp],varReco[ivar].Data() ));
	h_reco_tTagger[iwp]->SetMarkerColor(colors[iwp]);
	h_reco_tTagger[iwp]->SetLineColor(colors[iwp]);
	h_reco_tTagger[iwp]->SetMarkerStyle(initMarker);
	if(ivar == 0 || ivar ==1 || ivar ==3) h_reco_tTagger[iwp]->SetTitle(TString::Format("Acceptance Comparison;%s (GeV);Acceptance",varReco[ivar].Data()) );
	else h_reco_tTagger[iwp]->SetTitle(TString::Format("Acceptance Comparison;%s;Acceptance",varReco[ivar].Data()) );
	h_reco_tTagger[iwp]->Draw("same");
	if(ivar ==0) reco_legend->AddEntry(h_reco_tTagger[iwp], TString::Format("topTagger (%0.1f)",workingPoints[iwp]), "lp");
	  
   }  
	reco_Canvas[ivar]->cd();
	reco_legend->Draw();
	  
	parton_Canvas[ivar]->cd();
	parton_legend->Draw();

    
  }
  
}
