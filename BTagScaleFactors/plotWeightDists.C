void plotWeightDists(TString year)
{
  //open the ttbar nominal file
  TFile *inf = TFile::Open(TString::Format("%s/TopTaggerHisto_TT_NominalMC_100.root",year.Data()));
  int col[3] = {kAzure, kRed, kGreen+3};
  //get the histogram for combined files
  TH1F *hBTag = (TH1F*)inf->Get("hWt_bTagEvntWeight_2btag");
  hBTag->GetXaxis()->SetRangeUser(0.7,1.1);
  //hBTag->SetName("bTagEvntWeight Distribution in SR");
  hBTag->SetTitle(TString::Format("%s bTagEvntWeight Distribution in SR", year.Data()));

  int colorToUse;
  if(year.EqualTo("2016")) colorToUse = col[0];
  else if(year.EqualTo("2017")) colorToUse = col[1];
  else if(year.EqualTo("2018")) colorToUse = col[2];
  hBTag->SetLineColor(colorToUse);
  hBTag->SetFillColor(colorToUse);

  TCanvas *can = new TCanvas("bTagEvntWeight dist","bTagEvntWeight dist", 800,600);
  hBTag->Draw("hist");
  can->Print(TString::Format("%s/bTagEvntWeight_dist.pdf",year.Data()),"pdf");

}
