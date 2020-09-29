
void Contamination(TString year)
{
  gStyle->SetOptStat(0);

  //data template
  TFile *infData = TFile::Open(TString::Format("%s/Histo_Data_%s_100.root", year.Data(),year.Data()));
  //qcd template
  TFile *infQCDMC = TFile::Open(TString::Format("%s/Histo_QCD_HT300toInf_100.root", year.Data()));
  //ttbar
  TFile *infTTMC = TFile::Open(TString::Format("%s/Histo_TT_NominalMC_100.root",year.Data()));  //nominal
  //TFile *infTTMC = TFile::Open(TString::Format("%s/Histo_TT_Mtt-700toInf_100.root",year.Data())); //mtt
  //subdominant bkgs
  TFile *infSubBkg = TFile::Open(TString::Format("%s/Histo_SubdominantBkgs_100.root",year.Data()));  //nominal


  //get the histograms for all the variables:
  const int NVAR = 9;
  TString recoVarAll[NVAR] = {"jetPt0", "mJJ", "ptJJ", "yJJ", "jetPt1", "jetY0","jetY1", "mTop", "jetMassSoftDrop"};


  TH1F *hData[NVAR], *hCR_tt[NVAR], *hCR_Subdominant[NVAR], *hDataBefore[NVAR], *hCR_QCD[NVAR];
  TH1F *hRatio[NVAR][3];
  TCanvas *c1[NVAR];
  TPad *closure_pad2[NVAR], *closure_pad1[NVAR];
  TLegend *legend[NVAR];
  for(int ivar = 0; ivar<NVAR; ivar++)
  {
    TString recoVar = recoVarAll[ivar];
    hData[ivar] = (TH1F*)infData->Get(TString::Format("hWt_%s_0btag_expYield",recoVar.Data()));
    hCR_tt[ivar] = (TH1F*)infTTMC->Get(TString::Format("hWt_%s_0btag_expYield",recoVar.Data()));
    hCR_Subdominant[ivar] = (TH1F*)infSubBkg->Get(TString::Format("hWt_%s_0btag_expYield",recoVar.Data()));
    hCR_QCD[ivar] = (TH1F*)infQCDMC->Get(TString::Format("hWt_%s_0btag_expYield",recoVar.Data()));

    if(recoVar.EqualTo("mTop"))
    {
      cout<<"QCD entries: "<<hCR_QCD[ivar]->GetEntries()<<endl;
      cout<<"Subdominant entries: "<<hCR_Subdominant[ivar]->GetEntries()<<endl;
      cout<<"TT entries: "<<hCR_tt[ivar]->GetEntries()<<endl;
      cout<<"QCD Integral: "<<hCR_QCD[ivar]->Integral()<<endl;
      cout<<"Subdominant Integral: "<<hCR_Subdominant[ivar]->Integral()<<endl;
      cout<<"TT Integral: "<<hCR_tt[ivar]->Integral()<<endl;
      cout<<"-------------"<<endl;
    }
    //hData[ivar]->Add(hCR_tt,-1);
    //hData[ivar]->Add(hCR_Subdominant,-1);
    //hData->Draw();

    c1[ivar] = new TCanvas(TString::Format("can_%s",recoVar.Data()), TString::Format("can_%s",recoVar.Data()), 800,700);
    closure_pad2[ivar] = new TPad(TString::Format("cp2_%s",recoVar.Data()),TString::Format("cp2_%s",recoVar.Data()),0.,0.,1.,0.6);
    closure_pad2[ivar]->Draw();
    closure_pad2[ivar]->SetTopMargin(0.05);
    closure_pad2[ivar]->SetBottomMargin(0.2);
    closure_pad2[ivar]->SetGrid();

    closure_pad1[ivar] = new TPad(TString::Format("cp1_%s",recoVar.Data()),TString::Format("cp1_%s",recoVar.Data()),0.,0.6,1.,1.);
    closure_pad1[ivar]->Draw();
    closure_pad1[ivar]->SetBottomMargin(0.01);
    closure_pad1[ivar]->cd();

    hCR_QCD[ivar]->GetYaxis()->SetTitleSize(20);
    hCR_QCD[ivar]->GetYaxis()->SetTitleFont(43);
    hCR_QCD[ivar]->GetYaxis()->SetTitleOffset(1.4);


    if(recoVar.EqualTo("yJJ"))
    {
      hCR_QCD[ivar]->GetXaxis()->SetTitle(recoVar);
      hCR_tt[ivar]->GetXaxis()->SetTitle(recoVar);
      hCR_Subdominant[ivar]->GetXaxis()->SetTitle(recoVar);
    }
    else if(recoVar.EqualTo("jetY0") || recoVar.EqualTo("jetY1"))
    {
      hCR_QCD[ivar]->GetXaxis()->SetTitle(TString::Format("#||{%s}",recoVar.Data()));
      hCR_Subdominant[ivar]->GetXaxis()->SetTitle(TString::Format("#||{%s}",recoVar.Data()));
      hCR_tt[ivar]->GetXaxis()->SetTitle(TString::Format("#||{%s}",recoVar.Data()));
    }
    else
    {
      hCR_QCD[ivar]->GetXaxis()->SetTitle(recoVar + " (GeV)");
      hCR_tt[ivar]->GetXaxis()->SetTitle(recoVar + " (GeV)");
      hCR_Subdominant[ivar]->GetXaxis()->SetTitle(recoVar + " (GeV)");
    }

    if(recoVar.EqualTo("mTop") || recoVar.EqualTo("jetMassSoftDrop"))
    {
      hCR_QCD[ivar]->Rebin(2);
      hCR_tt[ivar]->Rebin(2);
      hCR_Subdominant[ivar]->Rebin(2);
    }

    hCR_QCD[ivar]->GetYaxis()->SetRangeUser(0, 1.1*hCR_QCD[ivar]->GetMaximum());

    hCR_QCD[ivar]->ResetAttLine();
    hCR_QCD[ivar]->ResetAttMarker();
    hCR_QCD[ivar]->SetLineColor(kRed);
    hCR_QCD[ivar]->SetMarkerStyle(21);
    hCR_QCD[ivar]->SetMarkerColor(kRed);

    hCR_tt[ivar]->ResetAttLine();
    hCR_tt[ivar]->ResetAttMarker();
    hCR_tt[ivar]->SetLineColor(kBlue);
    hCR_tt[ivar]->SetMarkerStyle(22);
    hCR_tt[ivar]->SetMarkerColor(kBlue);

    hCR_Subdominant[ivar]->ResetAttLine();
    hCR_Subdominant[ivar]->ResetAttMarker();
    hCR_Subdominant[ivar]->SetLineColor(kMagenta);
    hCR_Subdominant[ivar]->SetMarkerStyle(23);
    hCR_Subdominant[ivar]->SetMarkerColor(kMagenta);

    legend[ivar] = new TLegend(0.65,0.65,0.9,0.9);
    legend[ivar]->AddEntry(hCR_QCD[ivar], "QCD", "lep");
    legend[ivar]->AddEntry(hCR_tt[ivar], "ttbar", "lep");
    legend[ivar]->AddEntry(hCR_Subdominant[ivar], "Sub. Bkgs", "lep");


    hCR_QCD[ivar]->Draw();
    hCR_tt[ivar]->Draw("same");
    hCR_Subdominant[ivar]->Draw("same");
    legend[ivar]->Draw();

    //now the ratios:

    hRatio[ivar][0] = (TH1F*)hCR_QCD[ivar]->Clone(TString::Format("hRatio_QCD_%s",recoVar.Data()));
    hRatio[ivar][1] = (TH1F*)hCR_tt[ivar]->Clone(TString::Format("hRatio_Sub_%s",recoVar.Data()));
    hRatio[ivar][2] = (TH1F*)hCR_Subdominant[ivar]->Clone(TString::Format("hRatio_tt_%s",recoVar.Data()));


    //we want tt / qcd and sub / qcd
    hRatio[ivar][0]->Divide(hCR_QCD[ivar]);
    hRatio[ivar][1]->Divide(hCR_QCD[ivar]);
    hRatio[ivar][2]->Divide(hCR_QCD[ivar]);


    closure_pad2[ivar]->cd();
    hRatio[ivar][0]->ResetAttMarker();
    hRatio[ivar][1]->ResetAttMarker();
    hRatio[ivar][2]->ResetAttMarker();

    //hRatio->GetYaxis()->SetTitle(TString::Format("ratio %s/%s",hNum->GetTitle(),hDenom->GetTitle()));
    hRatio[ivar][0]->GetYaxis()->SetTitle("ratio Red/Blue");
    for(int i =0; i<3; i++)
    {
      hRatio[ivar][i]->SetTitle("");
      hRatio[ivar][i]->GetYaxis()->SetTitleSize(20);
      hRatio[ivar][i]->GetYaxis()->SetTitleFont(43);
      hRatio[ivar][i]->GetYaxis()->SetTitleOffset(1.55);
      hRatio[ivar][i]->GetYaxis()->SetLabelFont(43);
      hRatio[ivar][i]->GetYaxis()->SetLabelSize(15);
      hRatio[ivar][i]->GetXaxis()->SetTitleSize(0.06);
      hRatio[ivar][i]->GetYaxis()->SetRangeUser(-0.001,0.2);
      hRatio[ivar][i]->GetXaxis()->SetLabelSize(0.04);
    }

    //hRatio[ivar][0]->Draw();

    hRatio[ivar][1]->Draw("hist");
    hRatio[ivar][2]->Draw("hist same");
    //break;

    //c1[ivar]->Print(TString::Format("./ContaminationPlots/%s/Medium/contamination_%s.pdf",year.Data(),recoVar.Data()),"pdf");

  }

 }
