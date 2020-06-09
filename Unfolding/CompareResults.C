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
#include "TemplateConstantsUnfold.h"

TH1F *getRebinned(TH1F *h, float BND[], int N);


TH1F *getRebinned(TH1F *h, float BND[], int N)
{
  TString name = TString(h->GetName())+"_Rebinned";
  TH1F *hNew = new TH1F(name, name, N, BND);

  for(int i=0;i<h->GetNbinsX();i++) {
    float x = h->GetBinCenter(i+1);
    float y = h->GetBinContent(i+1);
    float e = h->GetBinError(i+1);
    for(int j=0;j<hNew->GetNbinsX();j++) {
      float x1 = hNew->GetBinLowEdge(j+1);
      float x2 = x1+hNew->GetBinWidth(j+1);
      if ((x>x1) && (x<x2)) {
        float yNew = hNew->GetBinContent(j+1);
        float eNew = hNew->GetBinError(j+1);
        hNew->SetBinContent(j+1,yNew+y);
        hNew->SetBinError(j+1,sqrt(e*e+eNew*eNew));
        break;
      }
    }
  }
  return hNew;

}

void CompareResults(TString year = "2016", bool isParton = true, int unfoldMethod = 1)
{
  gStyle->SetOptStat(0);
  initFilesMapping();

  TString varParton = "Parton";
  TString varParton_top = varParton;
  if(!isParton)
  {
    varParton = "Particle";
    varParton_top = "Gen";

  }

  TString unfMethodStr = "";
  if(unfoldMethod ==2) unfMethodStr = "_LCurveMethod";
  else if(unfoldMethod ==3) unfMethodStr = "_RhoMethod";

  const int NVAR = 7;
  TString variable[NVAR] = {"mJJ", "ptJJ", "yJJ", "jetPt0", "jetPt1","jetY0", "jetY1"};
  TString variableGen[NVAR] = {"mJJGen", "ptJJGen", "yJJGen", "genjetPt0", "genjetPt1","genjetY0", "genjetY1"};
  TString variableParton[NVAR] = {"mTTbarParton", "ptTTbarParton", "yTTbarParton", "partonPt0", "partonPt1","partonY0", "partonY1"};

  //Our file with the unfolding results:
  TFile *inf = TFile::Open(TString::Format("%s/%sMeasurements/Data/OutputFile%s.root", year.Data(), varParton.Data(), unfMethodStr.Data()));

  //top-18-013 files
  TFile *top18013file;

  //TH1F here for the unfolded
  TH1F *hUnfolded[NVAR], *hTOP18013[NVAR];
  //TH1F here for the normalized
  TH1F *hUnfNorm[NVAR], *hTOP18013Norm[NVAR];
  //for theory:
  TH1F *hTheory[NVAR], *hTheoryNorm[NVAR];

  TCanvas *can[NVAR], *canNorm[NVAR];
  TLegend *leg[NVAR], *legNorm[NVAR];

  for(int ivar =0; ivar<NVAR; ivar++)
  {
    cout<<variable[ivar]<<endl;
    //get our plots:
    hUnfolded[ivar] = (TH1F*)inf->Get(TString::Format("hUnfold_%s", variable[ivar].Data()));
    hTheory[ivar] = (TH1F*)inf->Get(TString::Format("hTheory_%s", variable[ivar].Data()));
    hUnfolded[ivar]->SetLineColor(kBlue);
    hUnfolded[ivar]->SetMarkerColor(kBlue);
    hUnfolded[ivar]->SetMarkerStyle(22);

    hTheory[ivar]->SetLineColor(kRed);
    hTheory[ivar]->SetMarkerColor(kRed);
    hTheory[ivar]->SetMarkerStyle(20);
    if(variable[ivar].Contains("jetY"))
    {
      if(isParton)
        hTheory[ivar]->GetXaxis()->SetTitle("|"+variableParton[ivar]+"|");
      else
        hTheory[ivar]->GetXaxis()->SetTitle("|"+variableGen[ivar]+"|");
    }

    //top-18-013 plots:
    top18013file = TFile::Open(TString::Format("Results-TOP18013/CrossSection_%s_%s.root",varParton_top.Data(), variable[ivar].Data()));
    hTOP18013[ivar] = (TH1F*)top18013file->Get(TString::Format("CrossSection_%s_Nominal",varParton_top.Data()));
    hTOP18013Norm[ivar] = (TH1F*)top18013file->Get(TString::Format("NormCrossSection_%s_Nominal",varParton_top.Data()));

    hTOP18013[ivar]->SetLineColor(kGreen+4);
    hTOP18013[ivar]->SetMarkerColor(kGreen+4);
    hTOP18013[ivar]->SetMarkerStyle(23);

    hTOP18013Norm[ivar]->SetLineColor(kGreen+4);
    hTOP18013Norm[ivar]->SetMarkerColor(kGreen+4);
    hTOP18013Norm[ivar]->SetMarkerStyle(23);

    //normalized:
    hUnfNorm[ivar] = (TH1F*)hUnfolded[ivar]->Clone(TString::Format("hUnfoldNormalised_%s", variable[ivar].Data()));
    hUnfNorm[ivar]->Scale(luminosity[year]/((TH1F*)inf->Get(TString::Format("hUnfoldFinal_%s", variable[ivar].Data())))->Integral());
    hTheoryNorm[ivar] = (TH1F*)hTheory[ivar]->Clone(TString::Format("hTheoryNormalised_%s", variable[ivar].Data()));
    hTheoryNorm[ivar]->Scale(luminosity[year]/((TH1F*)inf->Get(TString::Format("hTheoryFinal_%s", variable[ivar].Data())))->Integral());

    //draw the unfolded and extrapolated with the mc result
    can[ivar] = new TCanvas(TString::Format("can_%s",variable[ivar].Data()),TString::Format("can_%s",variable[ivar].Data()) , 800,600);
    can[ivar]->cd();
    leg[ivar] = new TLegend(0.65,0.7,0.9,0.9);

    hUnfolded[ivar]->SetTitle(TString::Format("%s Unfolded vs TOP18013 %s",varParton.Data(),year.Data()));
    hTheory[ivar]->SetTitle(TString::Format("%s Unfolded vs TOP18013 %s",varParton.Data(),year.Data()));
    hTOP18013[ivar]->SetTitle(TString::Format("%s Unfolded vs TOP18013 %s",varParton.Data(), year.Data()));

    hUnfolded[ivar]->GetYaxis()->SetTitle(TString::Format("#frac{d#sigma}{d(%s)} [pb]",variable[ivar].Data()));
    hTOP18013[ivar]->GetYaxis()->SetTitle(TString::Format("#frac{d#sigma}{d(%s)} [pb]",variable[ivar].Data()));
    hTheory[ivar]->GetYaxis()->SetTitle(TString::Format("#frac{d#sigma}{d(%s)} [pb]",variable[ivar].Data()));
    leg[ivar]->AddEntry(hTheory[ivar], "Theory", "lpe");
    leg[ivar]->AddEntry(hUnfolded[ivar], "Unfolded", "lpe");
    leg[ivar]->AddEntry(hTOP18013[ivar], "TOP-18- 013", "lpe");

    hTheory[ivar]->Draw("same");
    hUnfolded[ivar]->Draw("same");
    hTOP18013[ivar]->Draw("same");
  	leg[ivar]->Draw();
  	if(!variable[ivar].Contains("jetY")) gPad->SetLogy();

    //normalized plots:
    canNorm[ivar] = new TCanvas(TString::Format("canNorm_%s",variable[ivar].Data()),TString::Format("canNorm_%s",variable[ivar].Data()) , 800,600);
    canNorm[ivar]->cd();
    legNorm[ivar] = new TLegend(0.65,0.7,0.9,0.9);

    hUnfNorm[ivar]->SetTitle(TString::Format("Normalised %s Unfolded vs TOP18013 %s",varParton.Data(),year.Data()));
    hTheoryNorm[ivar]->SetTitle(TString::Format("Normalised %s Unfolded vs TOP18013 %s",varParton.Data(),year.Data()));
    hTOP18013Norm[ivar]->SetTitle(TString::Format("Normalised %s Unfolded vs TOP18013 %s",varParton.Data(),year.Data()));

    hUnfNorm[ivar]->GetYaxis()->SetTitle(TString::Format("#frac{1}{#sigma} #frac{d#sigma}{d(%s)}",variable[ivar].Data()));
    hTOP18013Norm[ivar]->GetYaxis()->SetTitle(TString::Format("#frac{1}{#sigma} #frac{d#sigma}{d(%s)}",variable[ivar].Data()));
    hTheoryNorm[ivar]->GetYaxis()->SetTitle(TString::Format("#frac{1}{#sigma} #frac{d#sigma}{d(%s)}",variable[ivar].Data()));
    legNorm[ivar]->AddEntry(hTheoryNorm[ivar], "Theory", "lpe");
    legNorm[ivar]->AddEntry(hUnfNorm[ivar], "Unfolded", "lpe");
    legNorm[ivar]->AddEntry(hTOP18013Norm[ivar], "TOP-18- 013", "lpe");

    hTheoryNorm[ivar]->Draw("same");
    hUnfNorm[ivar]->Draw("same");
    hTOP18013Norm[ivar]->Draw("same");
  	legNorm[ivar]->Draw();
    if(!variable[ivar].Contains("jetY")) gPad->SetLogy();


    can[ivar]->Print(TString::Format("ComparisonWithTOP-18-013/%s/%s/DiffCrossSection_%s%s.pdf", varParton.Data(),year.Data(), variable[ivar].Data(), unfMethodStr.Data()), "pdf");
    canNorm[ivar]->Print(TString::Format("ComparisonWithTOP-18-013/%s/%s/NormDiffCrossSection_%s%s.pdf",varParton.Data(), year.Data(), variable[ivar].Data(), unfMethodStr.Data()), "pdf");
  }


}

/* THESE ARE NEEDED ONLY FOR THE RATIOS
if(!variable[ivar].Contains("jetY"))
{
  hTOP18013[ivar] = getRebinned(hTOP18013[ivar],tempBNDGen, NBINS_GEN[ivar]);
  hTOP18013Norm[ivar] = getRebinned(hTOP18013Norm[ivar],tempBNDGen, NBINS_GEN[ivar]);
}

auto *closure_padRatio = new TPad("closure_pad2","closure_pad2",0.,0.,1.,0.3);
closure_padRatio->Draw();
closure_padRatio->SetTopMargin(0.05);
closure_padRatio->SetBottomMargin(0.3);
closure_padRatio->SetGrid();

auto *closure_pad1 = new TPad("closure_pad1","closure_pad1",0.,0.3,1.,1.);
closure_pad1->Draw();
closure_pad1->SetBottomMargin(0.005);
closure_pad1->cd();
closure_padRatio->cd();

hUnfTemp[ivar] = (TH1F*)hUnfolded[ivar]->Clone(TString::Format("hUnfTemp_%s", variable[ivar].Data()));
hTOP18013Temp[ivar] = (TH1F*)hTOP18013[ivar]->Clone(TString::Format("hTOP18013Temp_%s", variable[ivar].Data()));
hUnfTemp[ivar]->Divide(hTOP18013Temp[ivar]);
hUnfTemp[ivar]->Draw();

hUnfTemp[ivar]->SetTitle("");
hUnfTemp[ivar]->GetYaxis()->SetTitle("#frac{Unfolded}{TOP-18-013}");
hUnfTemp[ivar]->GetYaxis()->SetTitleSize(14);
hUnfTemp[ivar]->GetYaxis()->SetTitleFont(43);
hUnfTemp[ivar]->GetYaxis()->SetTitleOffset(1.55);
hUnfTemp[ivar]->GetYaxis()->SetLabelFont(43);
hUnfTemp[ivar]->GetYaxis()->SetLabelSize(15);
hUnfTemp[ivar]->GetXaxis()->SetTitleSize(0.09);
hUnfTemp[ivar]->GetYaxis()->SetRangeUser(0,2);

hUnfTemp[ivar]->SetLineColor(kRed);
hUnfTemp[ivar]->SetMarkerStyle(20);
hUnfTemp[ivar]->SetMarkerColor(kRed);
hUnfTemp[ivar]->Draw();
hUnfTemp[ivar]->GetXaxis()->SetLabelSize(0.09);
*/
