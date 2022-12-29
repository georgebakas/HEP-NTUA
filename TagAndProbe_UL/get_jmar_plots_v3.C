#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"

using std::cin;
using std::cout;
using std::endl;
#include "TemplateConstants.h"


void plot_diff_region(TH1F *hProbe, TH1F *hSR, TString variable, TString data_mc="Data", TString random_type="random", TString year="2018")
{
    //plot the normalized plots
    TCanvas *can = new TCanvas(TString::Format("can_%s_%s_%s", data_mc.Data(), random_type.Data(), variable.Data()), 
                                TString::Format("can_%s_%s_%s",data_mc.Data(), random_type.Data(), variable.Data()), 800, 600);
    TLegend *leg = new TLegend(0.4,0.7,0.6,0.9);
    TPad *closure_pad2 = new TPad(TString::Format("cp2_%s",variable.Data()),TString::Format("cp2_%s",variable.Data()),0.,0.,1.,0.3);
    closure_pad2->Draw();
    closure_pad2->SetTopMargin(0.05);
    closure_pad2->SetBottomMargin(0.25);
    closure_pad2->SetGrid();

    TPad *closure_pad1 = new TPad(TString::Format("cp1_%s",variable.Data()),TString::Format("cp2_%s",variable.Data()),0.,0.3,1.,1.);
    closure_pad1->Draw();
    closure_pad1->SetBottomMargin(0.01);
    closure_pad1->cd();
    hProbe->Draw();
    hSR->Draw("same");
    leg->AddEntry(hProbe, TString::Format("%s Probe", data_mc.Data()), "lep");
    leg->AddEntry(hSR, TString::Format("%s SR (%s jet)", data_mc.Data(),random_type.Data()), "lep");
    leg->Draw();

    closure_pad2->cd();
    TH1F *hDenom = (TH1F*)hSR->Clone("hDenom");
    TH1F *hNum = (TH1F*)hProbe->Clone("hNum");
    hNum->Divide(hDenom);
    hNum->SetTitle("");
    hNum->GetYaxis()->SetRangeUser(0,2);
    hNum->GetYaxis()->SetTitle("#frac{Data Probe}{Data SR}");
    hNum->GetXaxis()->SetTitle(variable);
    hNum->GetXaxis()->SetTitle(variable);
    hNum->GetYaxis()->SetTitleSize(20);
    hNum->GetYaxis()->SetTitleFont(43);
    hNum->GetYaxis()->SetTitleOffset(1.3);
    hNum->GetYaxis()->SetLabelFont(43);
    hNum->GetYaxis()->SetLabelSize(15);
    hNum->GetXaxis()->SetTitleSize(0.1);
    hNum->GetXaxis()->SetLabelFont(43);
    hNum->GetXaxis()->SetLabelSize(13);
    hNum->Draw();

    can->Print(TString::Format("%s/extra_jmar_analysis/%s_probe_%s_%s.pdf", year.Data(), data_mc.Data(), random_type.Data(), variable.Data()), "pdf");
}


void plotStackHisto_Variable(TString year, TFile *infData, TFile *infTT, TString variable)
{
  gStyle->SetOptStat(0);
  initFilesMapping();
  int NREGIONS = 3;
  //now get the histograms
  TH1F *hData[NREGIONS];
  hData[0] = (TH1F*)infData->Get(TString::Format("hProbe%s%s", (variable.Contains("top") ? "" : "jet"), variable.Data()));
  hData[1] = (TH1F*)infData->Get(TString::Format("h_SR_random%s%s", (variable.Contains("top") ? "" : "jet"), variable.Data()));
  hData[2] = (TH1F*)infData->Get(TString::Format("h_SR_double%s%s", (variable.Contains("top") ? "" : "jet"), variable.Data()));

  TH1F *hTT[NREGIONS];
  hTT[0] = (TH1F*)infTT->Get(TString::Format("hProbe%s%s", (variable.Contains("top") ? "" : "jet"), variable.Data()));
  hTT[1] = (TH1F*)infTT->Get(TString::Format("h_SR_random%s%s", (variable.Contains("top") ? "" : "jet"), variable.Data()));
  hTT[2] = (TH1F*)infTT->Get(TString::Format("h_SR_double%s%s", (variable.Contains("top") ? "" : "jet"), variable.Data()));

  for(int ir = 0; ir < NREGIONS; ir++)
  {
    hTT[ir]->Scale(1./hTT[ir]->Integral());
    hData[ir]->Scale(1./hData[ir]->Integral());

    hTT[ir]->GetYaxis()->SetTitle("Normalized Distribution");
    hTT[ir]->GetXaxis()->SetTitle(TString::Format("%s", variable.Data()));
    hData[ir]->GetYaxis()->SetTitle("Normalized Distribution");
    hData[ir]->GetXaxis()->SetTitle(TString::Format("%s", variable.Data()));
    
    if (ir ==0)
    {
        hTT[ir]->SetLineColor(kBlack);
        hTT[ir]->SetMarkerStyle(20);
        hTT[ir]->SetMarkerColor(kBlack);

        hData[ir]->SetLineColor(kBlack);
        hData[ir]->SetMarkerStyle(20);
        hData[ir]->SetMarkerColor(kBlack);
    }
    else
    {
        hTT[ir]->SetLineColor(kRed-9);
        hTT[ir]->SetMarkerColor(kRed-9);
        hTT[ir]->SetFillColor(kRed-9);

        hData[ir]->SetLineColor(kRed-9);
        hData[ir]->SetMarkerColor(kRed-9);
        hData[ir]->SetFillColor(kRed-9);
    }
  }
    // plot normalized distributions 
    plot_diff_region(hData[0], hData[1], variable, "Data", "random", year);
    plot_diff_region(hData[0], hData[2], variable, "Data", "double", year);

    plot_diff_region(hTT[0], hTT[1], variable, "MC", "random", year);
    plot_diff_region(hTT[0], hTT[2], variable, "MC", "double", year);

    // now plot the 2D plots 
    if (!variable.EqualTo("topTagger"))
    {
        TH2F *hData2D[3], *hMC2D[3];
        TString regions[3] = {"probe", "random", "double"};

        for (int ir=0; ir<NREGIONS; ir++)
        {
            hData2D[ir] = (TH2F*)infData->Get(TString::Format("h%sTopTagger_%s", variable.Data(), regions[ir].Data()));
            hMC2D[ir] = (TH2F*)infTT->Get(TString::Format("h%sTopTagger_%s", variable.Data(), regions[ir].Data()));
            hData2D[ir]->GetXaxis()->SetTitle(variable);
            hMC2D[ir]->GetXaxis()->SetTitle(variable);
            hData2D[ir]->GetYaxis()->SetTitle("Top Tagger Score");
            hMC2D[ir]->GetYaxis()->SetTitle("Top Tagger Score");
            hData2D[ir]->SetTitle(TString::Format("Data: %s - TopTagger Distribution (%s Region)", variable.Data(), (regions[ir].Contains("probe") ? "Probe" : "SR")));
            hMC2D[ir]->SetTitle(TString::Format("MC: %s - TopTagger Distribution (%s Region)", variable.Data(), (regions[ir].Contains("probe") ? "Probe" : "SR")));
        
            TCanvas *can_data = new TCanvas(TString::Format("can_%d_%s", ir, variable.Data()), 
                                TString::Format("can_%d_%s", ir, variable.Data()), 800, 600);
            hData2D[ir]->Draw("text colz");
            can_data->Print(TString::Format("%s/extra_jmar_analysis/TH2_Data_%s_topTagger.pdf", year.Data(),  variable.Data()),"pdf");

            TCanvas *can_mc = new TCanvas(TString::Format("can_mc_%d_%s", ir, variable.Data()), 
                                TString::Format("can_mc_%d_%s", ir, variable.Data()), 800, 600);
            hMC2D[ir]->Draw("text colz");
            can_mc->Print(TString::Format("%s/extra_jmar_analysis/TH2_MC_%s_topTagger.pdf", year.Data(),  variable.Data()),"pdf");

        }

        
        

    }
}


void get_jmar_plots_v3(TString year)
{

    //get the files from the directory
    //data file
    TFile *infData = TFile::Open(TString::Format("%s/TTbarExtraAnalysis_data.root", year.Data()));
    //tt nominal file:
    TFile *infTT = TFile::Open(TString::Format("%s/Nominal/combined/TTbarExtraAnalysis_TT_Nominal.root", year.Data()));

    //now add the # entries, etc into text files 
    //write to a txt file the Kolmogorov tests
    ofstream myfile;
    myfile.open (TString::Format("%s/TagNProbe_AnalysisResults.txt", year.Data()));
    myfile<<"------------"<< "\n";    
    myfile<<"Data || MC & Region & Events || Integral"<<"\n";

    float data_probe_entries = ((TH1F*)infData->Get("hProbejetPt"))->GetEntries();
    float data_sr_entries = ((TH1F*)infData->Get("h_SR_randomjetPt"))->GetEntries();

    myfile<<"Data & Tag N Probe & Entries: "<<data_probe_entries<<"\n";
    myfile<<"Data & SR and SR & Entries: "<<data_sr_entries<<"\n";

    float mc_probe_entries = ((TH1F*)infTT->Get("hProbejetPt"))->GetEntries();
    float mc_probe_integral = ((TH1F*)infTT->Get("hProbejetPt"))->Integral();
    
    float mc_sr_entries = ((TH1F*)infTT->Get("h_SR_randomjetPt"))->GetEntries();
    float mc_sr_integral = ((TH1F*)infTT->Get("h_SR_randomjetPt"))->Integral();

    myfile<<"MC & Tag N Probe & Entries: "<<mc_probe_entries<<"\n";
    myfile<<"MC & SR and SR & Entries: "<<mc_sr_entries<<"\n";

    myfile<<"MC & Tag N Probe & Integral: "<<mc_probe_integral<<"\n";
    myfile<<"MC & SR and SR & Integral: "<<mc_sr_integral<<"\n";

    TFile *inf_tagNProbe_tt = TFile::Open(TString::Format("%s/Nominal/combined/TagAndProbeHisto_1000_TT_Nominal.root", year.Data()));
    float tag_sr_mc_entries = ((TH1F*)inf_tagNProbe_tt->Get("hSRBTightAndSR_jetPt0_expYield"))->GetEntries();
    float tag_sr_mc_integral = ((TH1F*)inf_tagNProbe_tt->Get("hSRBTightAndSR_jetPt0_expYield"))->Integral();

    TFile *inf_tagNProbe_data = TFile::Open(TString::Format("../MassFit/%s/TagAndProbeHisto_Data.root", year.Data()));
    float tag_sr_data_entries = ((TH1F*)inf_tagNProbe_data->Get("hSRBTightAndSR_jetPt0_expYield"))->GetEntries();

    myfile<<"MC & Tag and SR & Entries: "<<tag_sr_mc_entries<<"\n";
    myfile<<"MC & Tag and SR & Integral: "<<tag_sr_mc_integral<<"\n";

    myfile<<"Data & Tag and SR & Entries: "<<tag_sr_data_entries<<"\n";

    myfile.close();

    const int NVAR =3;
    TString varReco[NVAR]   = {"Pt", "Eta", "topTagger"};
    for(int ivar = 0; ivar< NVAR; ivar++)
    {
        plotStackHisto_Variable(year, infData, infTT, varReco[ivar]);
    }
}