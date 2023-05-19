
#include <stdio.h>
#include "TemplateConstants.h"

void plotLeptons(TString year = "2018")
{
    gStyle->SetOptStat(0);
    TFile *_file0 = TFile::Open(TString::Format("%s/Nominal/combined/LeptonAnalysis_Nominal.root", year.Data()));
    
    // first get hadrons
    // this is tight and signal regions
    TH1F *hhTR = (TH1F*)_file0->Get("hHadronic_TightRegjetPt0_expYield");
    // this is tight and probe
    TH1F *hhTP = (TH1F*)_file0->Get("hHadronic_TightProbejetPt0_expYield");
    // this is both tight
    TH1F *hh = (TH1F*)_file0->Get("hHadronicjetPt0_expYield");

    // now leptons
    // this is tight jet and muon
    TH1F *hL = (TH1F*)_file0->Get("hLeptonsjetPt0_expYield");

    // scale all on their integral
    hhTR->Scale(1./hhTR->Integral());
    hhTP->Scale(1./hhTP->Integral());
    hh->Scale(1./hh->Integral());
    hL->Scale(1./hL->Integral());


    // plot 
    hL->SetLineColor(kRed);
    hh->SetLineColor(kBlue);
    hhTR->SetLineColor(kBlack);
    hhTP->SetLineColor(kMagenta);


    TCanvas *can = new TCanvas("effCanPt", "effCanPt", 800, 600);
    TLegend *leg;
    //if (year.Contains("2016")) leg = new TLegend(0.5, 0.75, 0.7, 0.9);
    //else 
    leg = new TLegend(0.7, 0.5, 0.9, 0.75);
    leg->AddEntry(hL,"Tight Jet + Lepton", "l");
    leg->AddEntry(hhTP,"1 Tight Jet + Probe", "l");
    leg->AddEntry(hhTR,"1 Tight Jet + 1 SR", "l");
    leg->AddEntry(hh,"2 Tight Jets", "l");


    hL->GetYaxis()->SetRangeUser(0, 0.4);
    hL->GetXaxis()->SetTitle("Leading Jet pT");
    hL->SetTitle("ttbar Channel Analysis");
    hL->SetName("ttbar Channel Analysis");
    hL->Draw();

    hhTP->Draw("same");
    hhTR->Draw("same");
    hh->Draw("same");
    
    leg->Draw();
    can->Print("leptonAnalysis.pdf", "pdf");
}
