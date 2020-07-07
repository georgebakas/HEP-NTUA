/*
  In this code I will calculate the tag and probe efficiency of the top tagger
  N(1) probe is: #events pass SRb with top tight cut and SRb without any top tagger CUT
  N(2) probe is: #events pass SRb with top tight cut and SR (standard top tagger cut)

  N(1) serves as the denominator
  N(2) serves as the numerator

  as efficiency:= N(2) / N(1) for ttbar files, and data - QCD - subdominant files
*/

#include <stdio.h>

void CalculateEfficiency(TString year = "2016")
{
  //get the files:
  TFile *infData = TFile::Open(TString::Format("%s/TagAndProbeHisto_Data_%s_100_reduced_UnequalBinning.root", year.Data(), year.Data()));

  TFile *infQCD = TFile::Open(TString::Format("%s/TagAndProbeHisto_QCD_HT300toInf_100_reduced_UnequalBinning.root", year.Data())); //qcd file
  TFile *infTT = TFile::Open(TString::Format("%s/TagAndProbeHisto_TT_NominalMC_100_reduced_UnequalBinning.root", year.Data())); //ttbar file
  TFile *infSub = TFile::Open(TString::Format("%s/TagAndProbeHisto_SubdominantBkgs_100_reduced_UnequalBinning.root", year.Data())); //subdominant file


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

  //remove bkg contributions from data
  for(int i =0; i<2; i++)
  {
    hData[i]->Add(hSub[i],-1);
    hData[i]->Add(hQCD[i],-1);
  }

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

  FILE *fp;
  TString str = TString::Format("%s/Output_%s.txt",year.Data(), year.Data());
  fp = fopen(str.Data(),"w");
  fprintf(fp, "eff data: %f ± %f\n",eff_data, eff_data_error);
  fprintf(fp, "eff ttbar: %f ± %f\n",eff_tt, eff_tt_error);


  cout<<"eff data: "<<eff_data<<" ± "<<eff_data_error<<endl;
  cout<<"eff ttbar: "<<eff_tt<<" ± "<<eff_tt_error<<endl;
  fclose(fp);

}
