/*
  In this code I will calculate the tag and probe efficiency of the top tagger
  N(1) probe is: #events pass SRb with top tight cut and SRb without any top tagger CUT
  N(2) probe is: #events pass SRb with top tight cut and SR (standard top tagger cut)

  N(1) serves as the denominator
  N(2) serves as the numerator

  as efficiency:= N(2) / N(1) for ttbar files, and data - QCD - subdominant files
*/

#include <stdio.h>
#include "TemplateConstants.h"

void CalculateEfficiency(TString year = "2016")
{
  initFilesMapping();
  //get the files from the directory
  //data file
  TFile *infData = TFile::Open(TString::Format("../MassFit/%s/TagAndProbeHisto_Data.root", year.Data()));
  //tt nominal file:
  TFile *infTT = TFile::Open(TString::Format("%s/Nominal/combined/TagAndProbeHisto_1000_TT_Nominal.root", year.Data()));
  //TFile *infTT = TFile::Open(TString::Format("../MassFit/%s/TagAndProbeHisto_1000_TT_Nominal.root", year.Data()));
  //qcd mc file
  TFile *infQCD = TFile::Open(TString::Format("../MassFit/%s/TagAndProbeHisto_QCD_HT300toInf.root",year.Data()));
  //subdominant file:
  TFile *infSub = TFile::Open(TString::Format("../MassFit/%s/TagAndProbeHisto_SubdominantBkgs.root",year.Data()));

  TString regions[2] = {"hSRBTightAndProbe_", "hSRBTightAndSR_"};
  //get the histograms for the jetPt0 variable
  //[0] is denominator and [1] is numerator
  TH1F *hData[2], *hQCD[2], *hSub[2], *hTT[2];
  hData[0] = (TH1F*)infData->Get("hSRBTightAndProbe_jetPt0_expYield");
  hData[1] = (TH1F*)infData->Get("hSRBTightAndSR_jetPt0_expYield");

  hTT[0] = (TH1F*)infTT->Get("hSRBTightAndProbe_jetPt0_expYield");
  hTT[1] = (TH1F*)infTT->Get("hSRBTightAndSR_jetPt0_expYield");

  hSub[0] = (TH1F*)infSub->Get("hSRBTightAndProbe_jetPt0_expYield");
  hSub[1] = (TH1F*)infSub->Get("hSRBTightAndSR_jetPt0_expYield");

  hQCD[0] = (TH1F*)infQCD->Get("hSRBTightAndProbe_jetPt0_expYield");
  hQCD[1] = (TH1F*)infQCD->Get("hSRBTightAndSR_jetPt0_expYield");

  cout<<"hData Numerator: "<<hData[1]->Integral()<<endl;
  cout<<"hData denominator: "<<hData[0]->Integral()<<endl;

  cout<<"hSub Numerator: "<<hSub[1]->Integral()<<endl;
  cout<<"hSub denominator: "<<hSub[0]->Integral()<<endl;

  cout<<"hQCD Numerator: "<<hQCD[1]->Integral()<<endl;
  cout<<"hQCD denominator: "<<hQCD[0]->Integral()<<endl;

  cout<<"hTT Numerator: "<<hTT[1]->Integral()<<endl;
  cout<<"hTT denominator: "<<hTT[0]->Integral()<<endl;
  
  for(int i =0; i<2; i++)
  { 
    // scale QCD to its shape
    // get the NQCD 
    TFile *masFitResultsFile = TFile::Open(TString::Format("%s/Nominal/MassFitResults_%sSignalTemplates_.root",
                                                          year.Data(),  
                                                          regions[i].Data()));
    RooWorkspace *w = (RooWorkspace*)masFitResultsFile->Get("w");
    RooRealVar *value = (RooRealVar *)w->var("nFitQCD_2b");

    float val = value->getValV();
    float nqcd_error = value->getError();
    masFitResultsFile->Close();
    
    
    Double_t integral_error;
    float integral_qcd = hQCD[i]->IntegralAndError(1, hQCD[i]->GetNbinsX(), integral_error);
    for(int ibin=1; ibin<=hQCD[i]->GetNbinsX(); ibin++)
    {
      float error = hQCD[i]->GetBinError(ibin);
      float bin_value = hQCD[i]->GetBinContent(ibin);
      float new_error = TMath::Sqrt(TMath::Power((error*val)/integral_qcd,2) + 
                              TMath::Power((nqcd_error*bin_value)/integral_qcd, 2) + 
                              TMath::Power(bin_value*val*integral_error/TMath::Power(integral_qcd,2),2));
      hQCD[i]->SetBinError(ibin, new_error);
      hQCD[i]->SetBinContent(ibin, val*hQCD[i]->GetBinContent(ibin)/integral_qcd);
    }
    //hQCD[i]->Scale(val/integral_qcd);
  
  } 
  for(int i =0; i<2; i++)
  { 

    for(int ibin=1; ibin<=hData[i]->GetNbinsX(); ibin++)
    { 
        float new_error = TMath::Sqrt(TMath::Power(hQCD[i]->GetBinError(ibin) ,2)+ 
                                TMath::Power(hSub[i]->GetBinError(ibin) ,2) +
                                TMath::Power(hData[i]->GetBinError(ibin), 2));
        hData[i]->SetBinError(ibin, new_error);
        hData[i]->SetBinContent(ibin, hData[i]->GetBinContent(ibin) - hSub[i]->GetBinContent(ibin) - hQCD[i]->GetBinContent(ibin));
    }
    
    //hData[i]->Add(hSub[i],-1);
    //hData[i]->Add(hQCD[i],-1);
  }


  
  //scale the ttbar with its signal strength
  //hTT[0]->Scale(ttbarSigStrength_TagNProbe[year.Data()]);
  //hTT[1]->Scale(ttbarSigStrength_TagNSR[year.Data()]);

  //now measure the efficiency for data and ttbar files
  float eff_data = hData[1]->Integral() / hData[0]->Integral();
  float eff_tt   = hTT[1]->Integral() / hTT[0]->Integral();

  Double_t error_data[2], error_tt[2];
  Double_t intData[2], intTT[2];

  for(int i=0; i<2; i++)
  {
    intData[i] = hData[i]->IntegralAndError(1,hData[i]->GetNbinsX(),error_data[i]);
    intTT[i] = hTT[i]->IntegralAndError(1,hTT[i]->GetNbinsX(),error_tt[i]);
  }

  //calculate errors:
  float eff_data_error = TMath::Sqrt( TMath::Power(error_data[1]/intData[0],2) + TMath::Power(error_data[0] * intData[1]/ TMath::Power(intData[0],2),2));
  float eff_tt_error = TMath::Sqrt( TMath::Power(error_tt[1]/intTT[0],2) + TMath::Power(error_tt[0] * intTT[1]/ TMath::Power(intTT[0],2),2));

  //now calculate it per pT region:
  //jetPt0 binning: {400,425,450,475,500,535,570,610,650,700,750,800,850,900,950,1025,1100,1200,1300,1400,1500}
  //[i][0] is the ptregion i, denominator and [i][1] is pt region i, numerator
  Double_t integralPtRegionsData[3][2],integralPtRegionsData_error[3][2];
  Double_t integralPtRegionsTT[3][2], integralPtRegionsTT_error[3][2];

  float eff_PtRegionsData[3];
  float eff_PtRegionsTT[3];
  //calculate errors in pt regions:
  float eff_data_errorPtRegions[3];
  float eff_tt_errorPtRegions[3];
  float eff_tt_errorPtRegions_systematic[3];


  float start[3] = {1,3,5};
  float end[3] = {2,4,6};
  // {450, 500, 570, 650, 800, 1100, 1500}, //jetPt0
  // loop on all regions
  for(int ipt=0; ipt<3;ipt++)
  {
    integralPtRegionsData[ipt][0] = hData[0]->IntegralAndError(start[ipt],end[ipt],integralPtRegionsData_error[ipt][0]);
    integralPtRegionsData[ipt][1] = hData[1]->IntegralAndError(start[ipt],end[ipt],integralPtRegionsData_error[ipt][1]);

    integralPtRegionsTT[ipt][0] = hTT[0]->IntegralAndError(start[ipt],end[ipt],integralPtRegionsTT_error[ipt][0]);
    integralPtRegionsTT[ipt][1] = hTT[1]->IntegralAndError(start[ipt],end[ipt],integralPtRegionsTT_error[ipt][1]);


    eff_PtRegionsData[ipt] = integralPtRegionsData[ipt][1]/integralPtRegionsData[ipt][0];
    //eff_PtRegionsTT[ipt] = integralPtRegionsTT[ipt][1]/integralPtRegionsTT[ipt][0];
    eff_PtRegionsTT[ipt] = ttbar_mc_tagnprobe_eff[year.Data()][start[ipt]];
    
    eff_data_errorPtRegions[ipt] = TMath::Sqrt( 
                      TMath::Power(integralPtRegionsData_error[ipt][1]/integralPtRegionsData[ipt][0],2) +  // 
                      TMath::Power(integralPtRegionsData_error[ipt][0] * integralPtRegionsData[ipt][1]/ TMath::Power(integralPtRegionsData[ipt][0],2),2)
                      
                      );
    // this is statistical error
    eff_tt_errorPtRegions[ipt] = TMath::Sqrt( TMath::Power(integralPtRegionsTT_error[ipt][1]/integralPtRegionsTT[ipt][0],2) + TMath::Power(integralPtRegionsTT_error[ipt][0] * integralPtRegionsTT[ipt][1]/ TMath::Power(integralPtRegionsTT[ipt][0],2),2));
    // this is systematic error
    eff_tt_errorPtRegions_systematic[ipt] = ttbar_mc_tagnprobe_eff_error[year.Data()][start[ipt]];
  }


  FILE *fp;
  TString str = TString::Format("%s/Output_%s.txt",year.Data(), year.Data());
  fp = fopen(str.Data(),"w");
  fprintf(fp, "Efficiency--\n" );
  fprintf(fp, "eff data: %f ± %f\n",eff_data, eff_data_error);
  fprintf(fp, "eff ttbar: %f ± (stat) %f ± (systematic) %f\n",ttbar_mc_tagnprobe_eff[year.Data()][0], eff_tt_error, ttbar_mc_tagnprobe_eff_error[year.Data()][0]);
  fprintf(fp, "-----------\n" );
  fprintf(fp, "Efficiency per Pt region\n");
  fprintf(fp, "eff data pT[400-600]: %f ± %f\n",eff_PtRegionsData[0], eff_data_errorPtRegions[0]);
  fprintf(fp, "eff ttbar pT[400-600]: %f ± (stat) %f ± (systematic) %f\n",ttbar_mc_tagnprobe_eff[year.Data()][1], eff_tt_errorPtRegions[0], eff_tt_errorPtRegions_systematic[0]);
  fprintf(fp, "-----------\n" );
  fprintf(fp, "eff data pT[600-800]: %f ± %f\n",eff_PtRegionsData[1], eff_data_errorPtRegions[1]);
  fprintf(fp, "eff ttbar pT[600-800]: %f ± (stat) %f ± (systematic) %f\n",ttbar_mc_tagnprobe_eff[year.Data()][3], eff_tt_errorPtRegions[1], eff_tt_errorPtRegions_systematic[1]);
  fprintf(fp, "-----------\n" );
  fprintf(fp, "eff data pT[800-Inf]: %f ± %f\n",eff_PtRegionsData[2], eff_data_errorPtRegions[2]);
  fprintf(fp, "eff ttbar pT[800-Inf]: %f ± (stat) %f ± (systematic) %f\n",ttbar_mc_tagnprobe_eff[year.Data()][5], eff_tt_errorPtRegions[2], eff_tt_errorPtRegions_systematic[2]);

  cout<<"eff data: "<<eff_data<<" ± "<<eff_data_error<<endl;
  cout<<"eff ttbar: "<<eff_tt<<" ± "<<eff_tt_error<<endl;
  fclose(fp);



  TCanvas *can = new TCanvas("effCanPt", "effCanPt", 800, 600);
  TLegend *leg;
  //if (year.Contains("2016")) leg = new TLegend(0.5, 0.75, 0.7, 0.9);
  //else 
  leg = new TLegend(0.5, 0.15, 0.7, 0.3);

  int n = 3;
  float xData[] = {1,2,3};
  float xTT[] = {1.2,2.2,3.2};
  float xTT_systematic[] = {1.4,2.4,3.4};
  std::vector<TString> x_labels = {"pT[400-600]", "pT[600-800]", "pT[800-Inf]"};
  float ex[]= {0,0,0};
  TGraphErrors *grData = new TGraphErrors(n,xData,eff_PtRegionsData,ex,eff_data_errorPtRegions);
  TGraphErrors *grTT = new TGraphErrors(n,xTT,eff_PtRegionsTT,ex,eff_tt_errorPtRegions);
  TGraphErrors *grTT_systematic = new TGraphErrors(n,xTT_systematic,eff_PtRegionsTT,ex,eff_tt_errorPtRegions_systematic);

  leg->AddEntry(grData,"Data", "lep");
  leg->AddEntry(grTT,"TT Statistical", "lep");
  leg->AddEntry(grTT_systematic,"TT (Sys) + (Stat)", "lep");
  TMultiGraph *mg = new TMultiGraph();
  grData->SetTitle("TGraphErrors Example");
  grData->SetMarkerColor(kBlue);
  grData->SetMarkerStyle(21);
  grTT->SetMarkerColor(kRed);
  grTT->SetMarkerStyle(20);
  grTT_systematic->SetMarkerColor(kMagenta);
  grTT_systematic->SetMarkerStyle(22);

  // data central line 
  Double_t plotX_data[3] = {1, 2, 4};
  Double_t plotY_data[3] = {eff_data, eff_data, eff_data};
  Double_t plotX_data_error[3] = {0, 0, 0};
  Double_t plotY_data_error[3] = {eff_data_error, eff_data_error, eff_data_error};
  TGraphErrors *data_central_err = new TGraphErrors(3,plotX_data,plotY_data, plotX_data_error, plotY_data_error);
  data_central_err->SetFillColor(kBlue);
  data_central_err->SetLineColor(kBlue);
  data_central_err->SetLineWidth(2.0);
  data_central_err->SetFillStyle(3004);

  // MC central line 
  Double_t plotX_mc[3] = {1, 2, 4};
  Double_t plotY_mc[3] = {ttbar_mc_tagnprobe_eff[year.Data()][0], ttbar_mc_tagnprobe_eff[year.Data()][0], ttbar_mc_tagnprobe_eff[year.Data()][0]};
  Double_t plotX_mc_error[3] = {0, 0, 0};
  Double_t plotY_mc_error[3] = {ttbar_mc_tagnprobe_eff_error[year.Data()][0], ttbar_mc_tagnprobe_eff_error[year.Data()][0], ttbar_mc_tagnprobe_eff_error[year.Data()][0]};
  TGraphErrors *mc_central_err = new TGraphErrors(3,plotX_mc,plotY_mc, plotX_mc_error, plotY_mc_error);
  mc_central_err->SetFillColor(kMagenta);
  mc_central_err->SetLineColor(kMagenta);
  mc_central_err->SetLineWidth(2.0);
  mc_central_err->SetFillStyle(3005);

  leg->AddEntry(data_central_err,"Total Eff. Data", "f");
  leg->AddEntry(mc_central_err,"Total Eff. MC", "f");

  mg->Add(grData, "AP");
  mg->Add(grTT, "AP");
  mg->Add(grTT_systematic, "AP");
  mg->Add(data_central_err, "A3L");
  mg->Add(mc_central_err, "same A3L");
  
  mg->GetYaxis()->SetRangeUser(0.4, 1.05);
  mg->Draw("a");
  //data_central_err->Draw("same a3");
  
  for (int i=0;i<n;i++)
  {
    //cout<<mg->GetXaxis()->GetBinLabel(i+1)<<endl;
    //mg->GetXaxis()->SetBinLabel(i+1,x_labels[i].Data());
    //cout<<mg->GetXaxis()->GetBinLabel(i+1)<<endl;
  } 
  leg->Draw();
  can->Update();
  can->Print(TString::Format("%s/Nominal/plots/Efficiency_perPtRegion.pdf",year.Data()),"pdf");

}
