#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"


using std::cin;
using std::cout; 
using std::endl;

void plotSignalBkg(TString varNameReco = "jetPt0")
{ 
  gStyle->SetPalette(kOcean);
  TFile *infSignal = TFile::Open("SignalRegion_Data_MC.root");
  TFile *infBkg    = TFile::Open("ControlRegion_Data_MC.root");
  
  TH1F *hSignal_tTagger[4], *hBkg_tTagger[4];
  TH1F *hSignal_oldMva[4], *hBkg_oldMva[4];
  
  //0 is from Data, 1 is from ttbarMC and 2 is from QCD MC
  
  //1st is the signal region for the tTagger and the oldMva
  hSignal_tTagger[0] = (TH1F*)infSignal->Get(TString::Format("%s_%sSignal_tTagger",varNameReco.Data(), "data"));
  hSignal_tTagger[1] = (TH1F*)infSignal->Get(TString::Format("%s_%sSignal_tTagger",varNameReco.Data(), "ttbarMC"));
  hSignal_tTagger[2] = (TH1F*)infSignal->Get(TString::Format("%s_%sSignal_tTagger",varNameReco.Data(), "bkgMC"));
  hSignal_tTagger[3] = (TH1F*)infSignal->Get(TString::Format("%s_%sSignal_tTagger",varNameReco.Data(), "subdominantBkgMC"));
  
  hSignal_oldMva[0] = (TH1F*)infSignal->Get(TString::Format("%s_%sSignal_oldMva",varNameReco.Data(), "data"));
  hSignal_oldMva[1] = (TH1F*)infSignal->Get(TString::Format("%s_%sSignal_oldMva",varNameReco.Data(), "ttbarMC"));
  hSignal_oldMva[2] = (TH1F*)infSignal->Get(TString::Format("%s_%sSignal_oldMva",varNameReco.Data(), "bkgMC"));
  hSignal_oldMva[3] = (TH1F*)infSignal->Get(TString::Format("%s_%sSignal_oldMva",varNameReco.Data(), "subdominantBkgMC"));
  
  //2nd is the bkg (control) region for the tTagger and the oldMva
  hBkg_tTagger[0] = (TH1F*)infBkg->Get(TString::Format("%s_%sBkg_tTagger",varNameReco.Data(), "data"));
  hBkg_tTagger[1] = (TH1F*)infBkg->Get(TString::Format("%s_%sBkg_tTagger",varNameReco.Data(), "ttbarMC"));
  hBkg_tTagger[2] = (TH1F*)infBkg->Get(TString::Format("%s_%sBkg_tTagger",varNameReco.Data(), "bkgMC"));
  hBkg_tTagger[3] = (TH1F*)infBkg->Get(TString::Format("%s_%sBkg_tTagger",varNameReco.Data(), "subdominantBkgMC"));
  
  hBkg_oldMva[0] = (TH1F*)infBkg->Get(TString::Format("%s_%sBkg_oldMva",varNameReco.Data(), "data"));
  hBkg_oldMva[1] = (TH1F*)infBkg->Get(TString::Format("%s_%sBkg_oldMva",varNameReco.Data(), "ttbarMC"));
  hBkg_oldMva[2] = (TH1F*)infBkg->Get(TString::Format("%s_%sBkg_oldMva",varNameReco.Data(), "bkgMC"));
  hBkg_oldMva[3] = (TH1F*)infBkg->Get(TString::Format("%s_%sBkg_oldMva",varNameReco.Data(), "subdominantBkgMC"));
  
  
  //now we need to plot
  //Stack for the signal region
  THStack *hStack_tTagger = new THStack("tTagger bkg+signal","tTagger bkg+signal");
  THStack *hStack_oldMva  = new THStack("oldMva bkg+signal","oldMva bkg+signal" );
  
  //let's make them pretty :) (Signal Region)
  hSignal_tTagger[0]->SetLineColor(kBlack);
  hSignal_tTagger[0]->SetMarkerStyle(2);
  hSignal_tTagger[0]->SetMarkerColor(kBlack);
  //hSignal_tTagger[0]->SetLineWidth(2);

  hSignal_tTagger[1]->SetLineColor(38);
  hSignal_tTagger[1]->SetFillColor(38);
  //hSignal_tTagger[1]->SetMarkerStyle(21); 
  hSignal_tTagger[1]->SetMarkerColor(38);
  
  hSignal_tTagger[2]->SetLineColor(46);
  hSignal_tTagger[2]->SetFillColor(46);
  //hSignal_tTagger[2]->SetMarkerStyle(21);
  hSignal_tTagger[2]->SetMarkerColor(46);
  
  hSignal_tTagger[3]->SetLineColor(29);
  hSignal_tTagger[3]->SetFillColor(29);
  hSignal_tTagger[3]->SetMarkerColor(29);
  
  hSignal_oldMva[0]->SetLineColor(kBlack);
  hSignal_oldMva[0]->SetMarkerColor(kBlack);
  hSignal_oldMva[0]->SetMarkerStyle(2);

  hSignal_oldMva[1]->SetLineColor(38);
  hSignal_oldMva[1]->SetFillColor(38);
  hSignal_oldMva[1]->SetMarkerColor(38);
  //hSignal_oldMva[1]->SetMarkerStyle(21);

  hSignal_oldMva[2]->SetLineColor(46);
  hSignal_oldMva[2]->SetFillColor(46);
  hSignal_oldMva[2]->SetMarkerColor(46);
  //hSignal_oldMva[2]->SetMarkerStyle(21); 

  hSignal_oldMva[3]->SetLineColor(29);
  hSignal_oldMva[3]->SetFillColor(29);
  hSignal_oldMva[3]->SetMarkerColor(29);

  hStack_tTagger->Add(hSignal_tTagger[3]);
  hStack_tTagger->Add(hSignal_tTagger[2]);
  hStack_tTagger->Add(hSignal_tTagger[1]);  
  
  hStack_oldMva->Add(hSignal_oldMva[3]);
  hStack_oldMva->Add(hSignal_oldMva[2]);
  hStack_oldMva->Add(hSignal_oldMva[1]);
  
  
  //now plot them
  
  TCanvas *can_tTagger = new TCanvas("tTagger canvas", "tTagger canvas", 800,600);
  TLegend *leg_tTagger = new TLegend(0.5, 0.6,0.7, 0.8);
  hStack_tTagger->Draw("hist");
  hSignal_tTagger[0]->Draw("same");
  //hSignal_tTagger[1]->Draw("same");
  //hSignal_tTagger[2]->Draw("same");
  leg_tTagger->AddEntry(hSignal_tTagger[0], "Data", "lep");
  leg_tTagger->AddEntry(hSignal_tTagger[1], "ttbar", "f");
  leg_tTagger->AddEntry(hSignal_tTagger[2], "QCD multijets", "f");
  leg_tTagger->AddEntry(hSignal_tTagger[3], "Subdominant Bkg", "f");
  leg_tTagger->Draw();
  
  TCanvas *can_oldMva = new TCanvas("oldMva canvas", "oldMva canvas", 800,600);
  TLegend *leg_oldMva = new TLegend(0.5, 0.6,0.7, 0.8);
  hStack_oldMva->Draw("hist");
  hSignal_oldMva[0]->Draw("same");
  //hSignal[1]->Draw("same");
  //hBkg[1]->Draw("same");
  leg_oldMva->AddEntry(hSignal_oldMva[0], "Data", "lep");
  leg_oldMva->AddEntry(hSignal_oldMva[1], "ttbar", "f");
  leg_oldMva->AddEntry(hSignal_oldMva[2], "QCD multijets", "f");
  leg_oldMva->AddEntry(hSignal_oldMva[3], "Subdominant Bkg", "f");
  leg_oldMva->Draw();
  
  
  TCanvas *can_Ratio = new TCanvas("ratio" , "ratio", 800, 600);
  TH1F *hSignal_DataClone = (TH1F*)hSignal_tTagger[0]->Clone("hSignal_DataClone"); 
  TH1F *hSignal_ttbarMC   = (TH1F*)hSignal_tTagger[1]->Clone("hSignal_ttbarMC"); 
  TH1F *hSignal_bkgMC     = (TH1F*)hSignal_tTagger[2]->Clone("hSignal_bkgMC"); 
  
  hSignal_ttbarMC->Add(hSignal_bkgMC);
  hSignal_DataClone->Divide(hSignal_ttbarMC);
  
  hSignal_DataClone->Draw();
  
  //Control Region
  
  THStack *hStackBkg_tTagger = new THStack("tTagger bkg+signal control region","tTagger bkg+signal control region");
  THStack *hStackBkg_oldMva  = new THStack("oldMva bkg+signal control region","oldMva bkg+signal control region");
  
  hBkg_tTagger[0]->SetLineColor(kBlack);
  hBkg_tTagger[0]->SetMarkerColor(kBlack);
  hBkg_tTagger[0]->SetMarkerStyle(2);

  hBkg_tTagger[1]->SetLineColor(38);
  hBkg_tTagger[1]->SetFillColor(38);
  //hBkg_tTagger[1]->SetMarkerStyle(21); 
  hBkg_tTagger[1]->SetMarkerColor(38);
  hStackBkg_tTagger->Add(hBkg_tTagger[1]);
  
  hBkg_tTagger[2]->SetLineColor(46);
  hBkg_tTagger[2]->SetFillColor(46);
  //hBkg_tTagger[2]->SetMarkerStyle(21);
  hBkg_tTagger[2]->SetMarkerColor(46);

  hBkg_tTagger[3]->SetLineColor(29);
  hBkg_tTagger[3]->SetFillColor(29);
  hBkg_tTagger[3]->SetMarkerColor(29);

  hStackBkg_tTagger->Add(hBkg_tTagger[2]);
  hStackBkg_tTagger->Add(hBkg_tTagger[3]);
  
  hBkg_oldMva[0]->SetLineColor(kBlack);
  hBkg_oldMva[0]->SetMarkerColor(kBlack);
  hBkg_oldMva[0]->SetMarkerStyle(2);

  hBkg_oldMva[1]->SetLineColor(38);
  hBkg_oldMva[1]->SetFillColor(38);
  hBkg_oldMva[1]->SetMarkerColor(38);
  //hBkg_oldMva[1]->SetMarkerStyle(21);

  hBkg_oldMva[2]->SetLineColor(46);
  hBkg_oldMva[2]->SetFillColor(46);
  hBkg_oldMva[2]->SetMarkerColor(46);
  //hBkg_oldMva[2]->SetMarkerStyle(21); 
  
  hBkg_oldMva[3]->SetLineColor(29);
  hBkg_oldMva[3]->SetFillColor(29);
  hBkg_oldMva[3]->SetMarkerColor(29);

  hStackBkg_oldMva->Add(hBkg_oldMva[1]);
  hStackBkg_oldMva->Add(hBkg_oldMva[2]);
  hStackBkg_oldMva->Add(hBkg_oldMva[3]);
  
  //now plot them
  
  TCanvas *can_tTagger2 = new TCanvas("tTagger canvas2", "tTagger canvas2", 800,600);
  TLegend *leg_tTagger2 = new TLegend(0.5, 0.6,0.7, 0.8);
  hStackBkg_tTagger->Draw("hist");
  hBkg_tTagger[0]->Draw("same");
  //hSignal_tTagger[1]->Draw("same");
  //hSignal_tTagger[2]->Draw("same");
  leg_tTagger2->AddEntry(hBkg_tTagger[0], "Data", "lep");
  leg_tTagger2->AddEntry(hBkg_tTagger[1], "ttbar", "f");
  leg_tTagger2->AddEntry(hBkg_tTagger[2], "QCD multijets", "f");
  leg_tTagger2->AddEntry(hBkg_tTagger[3], "Subdominant Bkg", "f");
  leg_tTagger2->Draw();
  
  TCanvas *can_oldMva2 = new TCanvas("oldMva canvas2", "oldMva canvas2", 800,600);
  TLegend *leg_oldMva2= new TLegend(0.5, 0.6,0.7, 0.8);
  hStackBkg_oldMva->Draw("hist");
  hBkg_oldMva[0]->Draw("same");
  //hSignal[1]->Draw("same");
  //hBkg[1]->Draw("same");
  leg_oldMva2->AddEntry(hBkg_oldMva[0], "Data", "lep");
  leg_oldMva2->AddEntry(hBkg_oldMva[1], "ttbar", "f");
  leg_oldMva2->AddEntry(hBkg_oldMva[2], "QCD multijets", "f");
  leg_oldMva2->AddEntry(hBkg_oldMva[3], "Subdominant Bkg", "f");
  leg_oldMva2->Draw();
  
}