#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TRatioPlot.h"

using std::cin;
using std::cout;
using std::endl;

/*
We are using this file in order to fit the ratio of the BKG SR over the Signal SR
*/


void FitQCDClosure(TString recoVar = "jetPt0")
{
  // {400,450,500,570,650,750,850,950,1100,1300,1500}, //jetPt0 
  //QCD file
  TFile *infBkg = TFile::Open("BkgOutput_AllRegions_0.10_deepCSV.root");
  //TFile *infBkg = TFile::Open("BkgOutput_AllRegions_0.00_CSVv2.root");
  TH1F *hBkg_CR, *hBkg_SR;  
  
  //TH1F used for the TRatioPlot
  hBkg_CR = (TH1F*)infBkg->Get(TString::Format("CR_tTagger_%s",recoVar.Data()));
  hBkg_SR = (TH1F*)infBkg->Get(TString::Format("SR_tTagger_%s",recoVar.Data()));

  //These are bkg Signal region and bkg region
  hBkg_CR->SetTitle("QCD Closure tTagger");
  hBkg_SR->SetTitle("QCD Closure tTagger");
    
  TH1F *hClosure;
  hClosure = (TH1F*)hBkg_SR->Clone("hClosure"); 
  hClosure->SetTitle("");
  hClosure->GetYaxis()->SetTitle("ratio SR/CR");
  hClosure->Divide(hBkg_CR);

  TF1 *f1;

  //Now fit the hClosure which is the Bkg Signal Region over the Bkg Control Region   
  if(recoVar.EqualTo("jetPt0") || recoVar.EqualTo("jetPt1") || recoVar.EqualTo("mJJ"))
  {
    f1 = new TF1(TString::Format("func_%s", recoVar.Data()),"(0.8+[0]*pow(x,[2]))/(1+[1]*pow(x,[2]+2))",400,1500);
    f1->SetParameters(1.2,1.5, 1);
    f1->SetParNames("beta","delta", "alpha");
    //tested, ok
  }
  else if(recoVar.EqualTo("ptJJ")) //not yet ok
  {
    f1 = new TF1(TString::Format("func_%s", recoVar.Data()),"(1-[0]*pow(x,[2]))/(1+[1]*pow(x,[2]+2))",0,1300);
    f1->SetParameters(1.2,1.2,1);
    f1->SetParNames("beta","delta", "alpha");
    //seems ok needs more testing, check with Kostas
  }


  hClosure->Fit(TString::Format("func_%s",recoVar.Data()));

  TFitResult fitResult;
  f1->SetFitResult(fitResult);

  TFile *outFile = new TFile("FitOutput.root", "UPDATE");
  fitResult.Write();
  f1->Write();
  //outFile->Close();
  //infBkg->Close();

}



