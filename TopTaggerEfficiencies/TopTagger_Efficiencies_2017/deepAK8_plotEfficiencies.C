void deepAK8_plotEfficiencies()
{
  std::vector<float> workingPoints = {0.6, 0.5, 0.4, 0.3, 0.2};
  std::vector<Color_t> colors = {kRed, kGreen, kBlack, kCyan, kOrange};
  TFile *f = TFile::Open("deepAK8_efficiencies_allVars_tTagger_0.2.root");
  TFile *f2 = TFile::Open("deepAK8_efficiencies_allVars_tTagger_0.2.root");
  
  
  TLegend *reco_legend = new TLegend(0.5, 0.5, 0.8, 0.8);
  TLegend *parton_legend = new TLegend(0.5, 0.5, 0.8, 0.8);
  
  TString varReco[5]   = {"mJJ", "jetPt","ptJJ", "yJJ",  "jetEta"}; 
  TString varParton[5] = {"mTTbarParton","partonPt", "ptTTbarParton", "yTTbarParton",  "partonEta"}; 
  TCanvas *reco_Canvas[5], *parton_Canvas[5];
  
  for(int ivar= 0; ivar< 5; ivar++)
  {
  
	reco_Canvas[ivar] = new TCanvas(TString::Format("reco_Canvas_%s",varReco[ivar].Data()),TString::Format("reco_Canvas_%s", varReco[ivar].Data()), 700, 600);
	parton_Canvas[ivar] = new TCanvas(TString::Format("parton_Canvas_%s",varParton[ivar].Data()), TString::Format("parton_Canvas_%s", varParton[ivar].Data()), 700, 600);
  
	  reco_Canvas[ivar]->cd();
	  TEfficiency* h_reco_tTagger = (TEfficiency*)f->Get(TString::Format("Sig_Reco_tTagger_0.2_%s",varReco[ivar].Data() ));
	  h_reco_tTagger->SetMarkerColor(kBlue);
	  h_reco_tTagger->SetLineColor(kBlue);
	  h_reco_tTagger->SetMarkerStyle(20);
	  if(ivar == 0 || ivar ==1 || ivar ==3) h_reco_tTagger->SetTitle(TString::Format("Acceptance Comparison;%s (GeV);Acceptance",varReco[ivar].Data()) );
	  else h_reco_tTagger->SetTitle(TString::Format("Acceptance Comparison;%s;Acceptance",varReco[ivar].Data()) );
	  h_reco_tTagger->Draw();
	  if(ivar ==0) reco_legend->AddEntry(h_reco_tTagger, "topTagger (0.3)", "lp");
	  
	  TEfficiency* h_reco_oldMva = (TEfficiency*) f->Get(TString::Format("Sig_Reco_oldMva_0.8_%s",varReco[ivar].Data() ));
	  h_reco_oldMva->SetMarkerColor(kMagenta);
	  h_reco_oldMva->SetLineColor(kMagenta);
	  h_reco_oldMva->SetMarkerStyle(21);
	  h_reco_oldMva->Draw("same");
	  //h_reco_tTagger->GetYaxis()->SetRangeUser(0, 0.15);
	  if(ivar ==0) reco_legend->AddEntry(h_reco_oldMva, "event MVA (0.8)", "lp");
	  
	  parton_Canvas[ivar]->cd();
	  TEfficiency* h_parton_tTagger = (TEfficiency*) f->Get(TString::Format("Sig_Parton_tTagger_0.2_%s", varParton[ivar].Data()));
	  h_parton_tTagger->SetMarkerColor(kBlue);
	  h_parton_tTagger->SetLineColor(kBlue);
	  h_parton_tTagger->SetMarkerStyle(20);
	  h_parton_tTagger->SetTitle("Efficiency Comparison");
	  h_parton_tTagger->Draw();
	  //h_parton_tTagger->GetYaxis()->SetRangeUser(0, 0.15);
	  if(ivar ==0) parton_legend->AddEntry(h_parton_tTagger, "topTagger (0.3)", "lp");
	  
	  TEfficiency* h_parton_oldMva = (TEfficiency*) f->Get(TString::Format("Sig_Parton_oldMva_0.8_%s", varParton[ivar].Data()));
	  h_parton_oldMva->SetMarkerColor(kMagenta);
	  h_parton_oldMva->SetLineColor(kMagenta);
	  h_parton_oldMva->SetMarkerStyle(21);
	  h_parton_oldMva->Draw("same");
	  //h_parton_tTagger->GetYaxis()->SetRangeUser(0, 0.15);
	  if(ivar ==0) parton_legend->AddEntry(h_parton_oldMva, "event MVA (0.8)", "lp");
	  
	  for(int i=0; i<workingPoints.size(); i++)
	  {
		reco_Canvas[ivar]->cd();
		TEfficiency* h_reco = (TEfficiency*) f2->Get(TString::Format("Sig_Reco_deepAK8_%.1f_%s", workingPoints[i], varReco[ivar].Data()));
		h_reco->SetLineColor(colors[i]);
		h_reco->SetMarkerColor(colors[i]);
	    h_reco->SetMarkerStyle(22);		
		h_reco->Draw("SAME");
		if(ivar ==0) reco_legend->AddEntry(h_reco, TString::Format("deepAK8 (%.1f)", workingPoints[i]), "lp");
		
		parton_Canvas[ivar]->cd();
		TEfficiency* h_parton = (TEfficiency*) f2->Get(TString::Format("Sig_Parton_deepAK8_%.1f_%s", workingPoints[i], varParton[ivar].Data()));
		h_parton->SetLineColor(colors[i]);
		h_parton->SetMarkerColor(colors[i]);
		h_parton->SetMarkerStyle(22);
		h_parton->Draw("SAME");
		if(ivar ==0) parton_legend->AddEntry(h_parton, TString::Format("deepAK8 (%.1f)", workingPoints[i]), "lp");
	  }
	  
	  reco_Canvas[ivar]->cd();
	  reco_legend->Draw();
	  
	  parton_Canvas[ivar]->cd();
	  parton_legend->Draw();

  }
  
}
