#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TUnfold.h"
#include "TUnfoldDensity.h"

using namespace std;

using std::cin;
using std::cout;
using std::endl;

void CompareResults(TString year = "2016")
{
  gStyle->SetOptStat(0);
  const int NVAR = 7;
  TString variable[NVAR] = {"mJJ", "ptJJ", "yJJ", "jetPt0", "jetPt1","jetY0", "jetY1"};

  //Our file with the unfolding results:
  TFile *inf;

  //top-18-013 files
  TFile *top18013file, *top18013file_theory;

  //TH1F here for the unfolded
  TH1F *hFiducial[NVAR], *hFiducialNorm[NVAR];
  //for theory 16,17,18:
  TH1F *hTheory[NVAR], *hTheoryNorm[NVAR];

  //TH1F here for the top18013
  TH1F *hTOP18013[NVAR], *hTOP18013Norm[NVAR];
  //th1f for the top18013 mc files
  TH1F *hTOP18013_theory[NVAR], *hTOP18013Norm_theory[NVAR];


  TCanvas *can[NVAR], *canNorm[NVAR];
  TLegend *leg[NVAR], *legNorm[NVAR];

  for(int ivar =0; ivar<NVAR; ivar++)
  {
    inf = TFile::Open(TString::Format("%s/FiducialMeasurement/UnequalBinning/SignalHistograms_%s.root", year.Data(), variable[ivar].Data()));
    cout<<variable[ivar]<<endl;
    //get our plots:
    hFiducial[ivar] = (TH1F*)inf->Get(TString::Format("hSignal_%s", variable[ivar].Data()));
    hTheory[ivar] = (TH1F*)inf->Get(TString::Format("hSMC_%s", variable[ivar].Data()));
    hFiducial[ivar]->SetLineColor(kBlue);
    hFiducial[ivar]->SetMarkerColor(kBlue);
    hFiducial[ivar]->SetMarkerStyle(22);

    hTheory[ivar]->SetLineColor(kRed);
    hTheory[ivar]->SetMarkerColor(kRed);
    hTheory[ivar]->SetMarkerStyle(20);
    if(variable[ivar].Contains("jetY"))
      hFiducial[ivar]->GetXaxis()->SetTitle("|"+variable[ivar]+"|");
    else if(variable[ivar].EqualTo("yJJ"))
      hFiducial[ivar]->GetXaxis()->SetTitle(variable[ivar]);
    else
      hFiducial[ivar]->GetXaxis()->SetTitle(variable[ivar]+" (GeV)");

    //normalized 16,17,18:
    hFiducialNorm[ivar] = (TH1F*)hFiducial[ivar]->Clone(TString::Format("hSignalNormalised_%s", variable[ivar].Data()));
    hFiducialNorm[ivar]->Scale(1./hFiducial[ivar]->Integral());
    hTheoryNorm[ivar] = (TH1F*)hTheory[ivar]->Clone(TString::Format("hTheoryNormalised_%s", variable[ivar].Data()));
    hTheoryNorm[ivar]->Scale(1./hTheory[ivar]->Integral());


    //top-18-013 plots:
    top18013file = TFile::Open(TString::Format("../Unfolding/Results-TOP18013/CrossSection_Parton_%s.root",variable[ivar].Data()));
    hTOP18013[ivar] = (TH1F*)top18013file->Get("CrossSection_Fiducial_Nominal");
    hTOP18013Norm[ivar] = (TH1F*)top18013file->Get("NormCrossSection_Fiducial_Nominal");

    hTOP18013_theory[ivar] = (TH1F*)top18013file->Get("Fiducial_CrossSection_PowhegPythia8");
    hTOP18013Norm_theory[ivar] = (TH1F*)top18013file->Get("Fiducial_NormCrossSection_PowhegPythia8");

    //now for 16,17,18
    hFiducial[ivar]->Divide(hTheory[ivar]);
    hFiducialNorm[ivar]->Divide(hTheoryNorm[ivar]);
    //TOP18013:
    hTOP18013[ivar]->Divide(hTOP18013_theory[ivar]);
    hTOP18013Norm[ivar]->Divide(hTOP18013Norm_theory[ivar]);

    //draw the unfolded and extrapolated with the mc result
    can[ivar] = new TCanvas(TString::Format("can_%s",variable[ivar].Data()),TString::Format("can_%s",variable[ivar].Data()) , 800,600);
    can[ivar]->cd();
    leg[ivar] = new TLegend(0.65,0.7,0.9,0.9);

    hFiducial[ivar]->SetTitle(TString::Format("Fiducial DataOverMC ratio (%s, TOP18013)",year.Data()));
    hTOP18013[ivar]->SetTitle(TString::Format("Fiducial DataOverMC ratio (%s, TOP18013}",year.Data()));

    hTOP18013[ivar]->SetLineColor(kRed);
    hTOP18013[ivar]->SetMarkerColor(kRed);
    hTOP18013[ivar]->SetMarkerStyle(23);

    hFiducial[ivar]->GetYaxis()->SetTitle("#frac{Data}{Theory}");
    hTOP18013[ivar]->GetYaxis()->SetTitle("#frac{Data}{Theory}");
    hFiducial[ivar]->GetYaxis()->SetRangeUser(0,3);
    hTOP18013[ivar]->GetYaxis()->SetRangeUser(0,3);
    leg[ivar]->AddEntry(hFiducial[ivar], "New Measurement", "lpe");
    leg[ivar]->AddEntry(hTOP18013[ivar], "TOP-18-013", "lpe");

    hFiducial[ivar]->Draw("same");
    hTOP18013[ivar]->Draw("same");
  	leg[ivar]->Draw();

    //normalized plots:
    canNorm[ivar] = new TCanvas(TString::Format("canNorm_%s",variable[ivar].Data()),TString::Format("canNorm_%s",variable[ivar].Data()) , 800,600);
    canNorm[ivar]->cd();
    legNorm[ivar] = new TLegend(0.65,0.7,0.9,0.9);

    hFiducialNorm[ivar]->SetTitle(TString::Format("Normalised Fiducial DataOverMC ratio (%s, TOP18013}",year.Data()));
    hTOP18013Norm[ivar]->SetTitle(TString::Format("Normalised Fiducial DataOverMC ratio (%s, TOP18013}",year.Data()));

    hTOP18013Norm[ivar]->SetLineColor(kRed);
    hTOP18013Norm[ivar]->SetMarkerColor(kRed);
    hTOP18013Norm[ivar]->SetMarkerStyle(23);

    hFiducialNorm[ivar]->GetYaxis()->SetTitle("Norm #frac{Data}{Theory}");
    hTOP18013Norm[ivar]->GetYaxis()->SetTitle("Norm #frac{Data}{Theory}");
    hFiducialNorm[ivar]->GetYaxis()->SetRangeUser(0,3);
    hTOP18013Norm[ivar]->GetYaxis()->SetRangeUser(0,3);
    legNorm[ivar]->AddEntry(hFiducialNorm[ivar], "New Measurement ", "lpe");
    legNorm[ivar]->AddEntry(hTOP18013Norm[ivar], "TOP-18-013", "lpe");

    hFiducialNorm[ivar]->Draw("same");
    hTOP18013Norm[ivar]->Draw("same");
  	legNorm[ivar]->Draw();

    can[ivar]->Print(TString::Format("ComparisonWithTOP-18-013/%s/Fiducial_%s.pdf", year.Data(), variable[ivar].Data()), "pdf");
    canNorm[ivar]->Print(TString::Format("ComparisonWithTOP-18-013/%s/NormFiducial_%s.pdf", year.Data(), variable[ivar].Data()), "pdf");
  }


}
