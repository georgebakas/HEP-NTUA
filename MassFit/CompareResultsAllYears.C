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

void CompareResultsAllYears()
{
  gStyle->SetOptStat(0);

  const int NVAR = 7;
  TString variable[NVAR] = {"mJJ", "ptJJ", "yJJ", "jetPt0", "jetPt1","jetY0", "jetY1"};

  //Our file with the unfolding results:
  TFile *inf[3];

  //TH1F here for the unfolded
  TH1F *hFiducial[3][NVAR], *hFiducialNorm[3][NVAR];
  //clones for double ratio
  TH1F *hFiducial_Clone[3][NVAR], *hFiducialNorm_Clone[3][NVAR];
  //for theory 16,17,18:
  TH1F *hTheory[3][NVAR], *hTheoryNorm[3][NVAR];

  TString years[] = {"2016", "2017", "2018"};
  Int_t colors[] = {kBlue, kRed, kGreen+2};

  TCanvas *can[NVAR], *canNorm[NVAR];
  TLegend *leg[NVAR], *legNorm[NVAR];

  TPad *closure_padRatio[NVAR], *closure_padRatio_norm[NVAR];
  TPad *closure_pad1[NVAR], *closure_pad1_norm[NVAR];
  TH1F *hTemp, *hTemp_Norm;
  for(int ivar =0; ivar<NVAR; ivar++)
  {
    inf[0] = TFile::Open(TString::Format("2016/FiducialMeasurement/EqualBinning/SignalHistograms_%s.root", variable[ivar].Data()));
    inf[1] = TFile::Open(TString::Format("2017/FiducialMeasurement/EqualBinning/SignalHistograms_%s.root", variable[ivar].Data()));
    inf[2] = TFile::Open(TString::Format("2018/FiducialMeasurement/EqualBinning/SignalHistograms_%s.root", variable[ivar].Data()));

    cout<<variable[ivar]<<endl;
    can[ivar] = new TCanvas(TString::Format("can_%s",variable[ivar].Data()),TString::Format("can_%s",variable[ivar].Data()) , 800,600);
    can[ivar]->cd();
    closure_pad1[ivar] =  new TPad(TString::Format("closure_pad1%d",ivar),TString::Format("closure_pad1%d", ivar),0.,0.3,1.,1.);
    closure_pad1[ivar]->Draw();
    closure_padRatio[ivar] = new TPad(TString::Format("closure_pad2%d",ivar),TString::Format("closure_pad2%d", ivar),0.,0.,1.,0.3);
    closure_padRatio[ivar]->Draw();

    canNorm[ivar] = new TCanvas(TString::Format("canNorm_%s",variable[ivar].Data()),TString::Format("canNorm_%s",variable[ivar].Data()) , 800,600);
    canNorm[ivar]->cd();
    closure_padRatio_norm[ivar] = new TPad(TString::Format("norm_closure_pad2%d",ivar),TString::Format("closure_pad2%d", ivar),0.,0.,1.,0.3);
    closure_padRatio_norm[ivar]->Draw();
    closure_pad1_norm[ivar] =  new TPad(TString::Format("norm_closure_pad1%d",ivar),TString::Format("closure_pad1%d", ivar),0.,0.3,1.,1.);
    closure_pad1_norm[ivar]->Draw();

    leg[ivar] = new TLegend(0.65,0.7,0.9,0.9);
    legNorm[ivar] = new TLegend(0.65,0.7,0.9,0.9);

    for(int iy = 0; iy<3; iy++)
    {
      //get our plots:
      hFiducial[iy][ivar] = (TH1F*)inf[iy]->Get(TString::Format("hSignal_%s", variable[ivar].Data()));
      hTheory[iy][ivar] = (TH1F*)inf[iy]->Get(TString::Format("hSMC_%s", variable[ivar].Data()));
      hFiducial[iy][ivar]->SetLineColor(colors[iy]);
      hFiducial[iy][ivar]->SetMarkerColor(colors[iy]);
      hFiducial[iy][ivar]->SetMarkerStyle(20+iy);

      if(variable[ivar].Contains("jetY"))
        hTheory[iy][ivar]->GetXaxis()->SetTitle("|"+variable[ivar]+"|");

      //normalized 16,17,18:
      hFiducialNorm[iy][ivar] = (TH1F*)hFiducial[iy][ivar]->Clone(TString::Format("hSignalNormalised_%s", variable[ivar].Data()));
      hFiducialNorm[iy][ivar]->Scale(1./hFiducial[iy][ivar]->Integral());
      hTheoryNorm[iy][ivar] = (TH1F*)hTheory[iy][ivar]->Clone(TString::Format("hSMCNorm_%s", variable[ivar].Data()));
      hTheoryNorm[iy][ivar]->Scale(1./hTheory[iy][ivar]->Integral());
      //now for 16,17,18
      hFiducial[iy][ivar]->Divide(hTheory[iy][ivar]);
      hFiducialNorm[iy][ivar]->Divide(hTheoryNorm[iy][ivar]);
      //get the clones of these because I need to do their ratios also (double ratio)
      hFiducial_Clone[iy][ivar] = (TH1F*)hFiducial[iy][ivar]->Clone(TString::Format("%s_Clone", hFiducial[iy][ivar]->GetName()));
      hFiducialNorm_Clone[iy][ivar] = (TH1F*)hFiducialNorm[iy][ivar]->Clone(TString::Format("%s_Clone", hFiducialNorm[iy][ivar]->GetName()));

      //draw the unfolded and extrapolated with the mc result
      can[ivar]->cd();
      //if(iy==0)closure_pad1[ivar]->Draw();
      closure_pad1[ivar]->SetBottomMargin(0.005);
      closure_pad1[ivar]->cd();

      hFiducial[iy][ivar]->SetTitle("Ficucial Ratio Comparison ('16,'17,'18)");
      hFiducial[iy][ivar]->GetYaxis()->SetTitle("#frac{Data}{Theory}");
      hFiducial[iy][ivar]->GetYaxis()->SetRangeUser(0,2);
      leg[ivar]->AddEntry(hFiducial[iy][ivar], TString::Format("%s",years[iy].Data()), "lpe");
      //for the double ratio:
      //divide with 2016 as reference
      if(iy==0)
      {
          hTemp = (TH1F*)hFiducial_Clone[0][ivar]->Clone("hTemp");
          hTemp_Norm= (TH1F*)hFiducialNorm_Clone[0][ivar]->Clone("hTemp_Norm");
      }
      hFiducial_Clone[iy][ivar]->Divide(hTemp);
      hFiducial[iy][ivar]->Draw("same");
    	leg[ivar]->Draw();

      closure_padRatio[ivar]->SetTopMargin(0.05);
      closure_padRatio[ivar]->SetBottomMargin(0.3);
      //closure_padRatio[ivar]->SetGrid();
      closure_padRatio[ivar]->cd();

      hFiducial_Clone[iy][ivar]->SetTitle("");
      hFiducial_Clone[iy][ivar]->GetYaxis()->SetTitle("#frac{Other year}{2016}");
      hFiducial_Clone[iy][ivar]->GetYaxis()->SetTitleSize(14);
      hFiducial_Clone[iy][ivar]->GetYaxis()->SetTitleFont(43);
      hFiducial_Clone[iy][ivar]->GetYaxis()->SetTitleOffset(1.55);
      hFiducial_Clone[iy][ivar]->GetYaxis()->SetLabelFont(43);
      hFiducial_Clone[iy][ivar]->GetYaxis()->SetLabelSize(15);
      hFiducial_Clone[iy][ivar]->GetXaxis()->SetLabelSize(0.09);
      hFiducial_Clone[iy][ivar]->GetXaxis()->SetTitleSize(0.09);
      hFiducial_Clone[iy][ivar]->GetYaxis()->SetRangeUser(0,2);

      hFiducial_Clone[iy][ivar]->Draw("hist same");

      //normalized plots:

      canNorm[ivar]->cd();

      hFiducialNorm[iy][ivar]->SetTitle("Normalized Fiducial Ratio Comparison ('16,'17,'18)");
      hFiducialNorm[iy][ivar]->GetYaxis()->SetTitle("Normalised #frac{Data}{Theory}");
      hFiducialNorm[iy][ivar]->GetYaxis()->SetRangeUser(0,2);
      legNorm[ivar]->AddEntry(hFiducial[iy][ivar], TString::Format("%s",years[iy].Data()), "lpe");
      hFiducialNorm_Clone[iy][ivar]->Divide(hTemp_Norm);
      closure_padRatio_norm[ivar]->SetTopMargin(0.05);
      closure_padRatio_norm[ivar]->SetBottomMargin(0.3);
      //closure_padRatio_norm[ivar]->SetGrid();
      closure_padRatio_norm[ivar]->cd();
      hFiducialNorm_Clone[iy][ivar]->Draw("hist same");
      hFiducialNorm_Clone[iy][ivar]->SetTitle("");
      hFiducialNorm_Clone[iy][ivar]->GetYaxis()->SetTitle("Norm #frac{Other year}{2016}");
      hFiducialNorm_Clone[iy][ivar]->GetYaxis()->SetTitleSize(14);
      hFiducialNorm_Clone[iy][ivar]->GetYaxis()->SetTitleFont(43);
      hFiducialNorm_Clone[iy][ivar]->GetYaxis()->SetTitleOffset(1.55);
      hFiducialNorm_Clone[iy][ivar]->GetYaxis()->SetLabelFont(43);
      hFiducialNorm_Clone[iy][ivar]->GetYaxis()->SetLabelSize(15);
      hFiducialNorm_Clone[iy][ivar]->GetXaxis()->SetLabelSize(0.09);
      hFiducialNorm_Clone[iy][ivar]->GetXaxis()->SetTitleSize(0.09);
      hFiducialNorm_Clone[iy][ivar]->GetYaxis()->SetRangeUser(0,2);
      closure_pad1_norm[ivar]->cd();
      closure_pad1_norm[ivar]->SetBottomMargin(0.005);

      hFiducialNorm[iy][ivar]->Draw("same");
      legNorm[ivar]->Draw();

    }//------end of loop on years

    can[ivar]->Print(TString::Format("ComparisonAllYears/FiducialCrossSection_%s.pdf", variable[ivar].Data()), "pdf");
    canNorm[ivar]->Print(TString::Format("ComparisonAllYears/NormFiducialCrossSection_%s.pdf", variable[ivar].Data()), "pdf");
    //break;
  }//------ end of loop on nvars


}
