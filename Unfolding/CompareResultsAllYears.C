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


void CompareResultsAllYears(bool isParton = true, int unfoldMethod = 1)
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
  TFile *inf[3];
  inf[0] = TFile::Open(TString::Format("2016/%sMeasurements/Data/OutputFile.root", varParton.Data()));
  inf[1] = TFile::Open(TString::Format("2017/%sMeasurements/Data/OutputFile.root", varParton.Data()));
  inf[2] = TFile::Open(TString::Format("2018/%sMeasurements/Data/OutputFile.root", varParton.Data()));

  //TH1F here for the unfolded
  TH1F *hUnfolded[3][NVAR], *hUnfNorm[3][NVAR];
  //clones for double ratio
  TH1F *hUnfolded_Clone[3][NVAR], *hUnfNorm_Clone[3][NVAR];
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
      hUnfolded[iy][ivar] = (TH1F*)inf[iy]->Get(TString::Format("hUnfold_%s", variable[ivar].Data()));
      hTheory[iy][ivar] = (TH1F*)inf[iy]->Get(TString::Format("hTheory_%s", variable[ivar].Data()));
      hUnfolded[iy][ivar]->SetLineColor(colors[iy]);
      hUnfolded[iy][ivar]->SetMarkerColor(colors[iy]);
      hUnfolded[iy][ivar]->SetMarkerStyle(20+iy);

      if(variable[ivar].Contains("jetY"))
      {
        if(isParton)
          hTheory[iy][ivar]->GetXaxis()->SetTitle("|"+variableParton[ivar]+"|");
        else
          hTheory[iy][ivar]->GetXaxis()->SetTitle("|"+variableGen[ivar]+"|");
      }

      //normalized 16,17,18:
      hUnfNorm[iy][ivar] = (TH1F*)hUnfolded[iy][ivar]->Clone(TString::Format("hUnfoldNormalised_%s", variable[ivar].Data()));
      hUnfNorm[iy][ivar]->Scale(luminosity[years[iy]]/((TH1F*)inf[iy]->Get(TString::Format("hUnfoldFinal_%s", variable[ivar].Data())))->Integral());
      hTheoryNorm[iy][ivar] = (TH1F*)hTheory[iy][ivar]->Clone(TString::Format("hTheoryNormalised_%s", variable[ivar].Data()));
      hTheoryNorm[iy][ivar]->Scale(luminosity[years[iy]]/((TH1F*)inf[iy]->Get(TString::Format("hTheoryFinal_%s", variable[ivar].Data())))->Integral());

      //now for 16,17,18
      hUnfolded[iy][ivar]->Divide(hTheory[iy][ivar]);
      hUnfNorm[iy][ivar]->Divide(hTheoryNorm[iy][ivar]);

      //get the clones of these because I need to do their ratios also (double ratio)
      hUnfolded_Clone[iy][ivar] = (TH1F*)hUnfolded[iy][ivar]->Clone(TString::Format("%s_Clone", hUnfolded[iy][ivar]->GetName()));
      hUnfNorm_Clone[iy][ivar] = (TH1F*)hUnfNorm[iy][ivar]->Clone(TString::Format("%s_Clone", hUnfNorm[iy][ivar]->GetName()));

      //draw the unfolded and extrapolated with the mc result
      can[ivar]->cd();
      //if(iy==0)closure_pad1[ivar]->Draw();
      closure_pad1[ivar]->SetBottomMargin(0.005);
      closure_pad1[ivar]->cd();

      hUnfolded[iy][ivar]->SetTitle(TString::Format("%s Unfolded Ratio Comparison ('16,'17,'18)",varParton.Data()));
      hUnfolded[iy][ivar]->GetYaxis()->SetTitle("#frac{Data}{Theory}");
      hUnfolded[iy][ivar]->GetYaxis()->SetRangeUser(0,2);
      leg[ivar]->AddEntry(hUnfolded[iy][ivar], TString::Format("%s",years[iy].Data()), "lpe");

      //for the double ratio:
      //divide with 2016 as reference
      if(iy==0)
      {
          hTemp = (TH1F*)hUnfolded_Clone[0][ivar]->Clone("hTemp");
          hTemp_Norm= (TH1F*)hUnfNorm_Clone[0][ivar]->Clone("hTemp_Norm");
      }
      hUnfolded_Clone[iy][ivar]->Divide(hTemp);

      hUnfolded[iy][ivar]->Draw("same");
    	leg[ivar]->Draw();

      closure_padRatio[ivar]->SetTopMargin(0.05);
      closure_padRatio[ivar]->SetBottomMargin(0.3);
      //closure_padRatio[ivar]->SetGrid();
      closure_padRatio[ivar]->cd();

      hUnfolded_Clone[iy][ivar]->SetTitle("");
      hUnfolded_Clone[iy][ivar]->GetYaxis()->SetTitle("#frac{Other year}{2016}");
      hUnfolded_Clone[iy][ivar]->GetYaxis()->SetTitleSize(14);
      hUnfolded_Clone[iy][ivar]->GetYaxis()->SetTitleFont(43);
      hUnfolded_Clone[iy][ivar]->GetYaxis()->SetTitleOffset(1.55);
      hUnfolded_Clone[iy][ivar]->GetYaxis()->SetLabelFont(43);
      hUnfolded_Clone[iy][ivar]->GetYaxis()->SetLabelSize(15);
      hUnfolded_Clone[iy][ivar]->GetXaxis()->SetLabelSize(0.09);
      hUnfolded_Clone[iy][ivar]->GetXaxis()->SetTitleSize(0.09);
      hUnfolded_Clone[iy][ivar]->GetYaxis()->SetRangeUser(0,2);

      hUnfolded_Clone[iy][ivar]->Draw("hist same");

      //normalized plots:

      canNorm[ivar]->cd();

      hUnfNorm[iy][ivar]->SetTitle(TString::Format("Normalized %s Unfolded Ratio Comparison ('16,'17,'18)",varParton.Data()));
      hUnfNorm[iy][ivar]->GetYaxis()->SetTitle("Normalised #frac{Data}{Theory}");
      hUnfNorm[iy][ivar]->GetYaxis()->SetRangeUser(0,2);
      legNorm[ivar]->AddEntry(hUnfolded[iy][ivar], TString::Format("%s",years[iy].Data()), "lpe");

      hUnfNorm_Clone[iy][ivar]->Divide(hTemp_Norm);

      closure_padRatio_norm[ivar]->SetTopMargin(0.05);
      closure_padRatio_norm[ivar]->SetBottomMargin(0.3);
      //closure_padRatio_norm[ivar]->SetGrid();
      closure_padRatio_norm[ivar]->cd();
      hUnfNorm_Clone[iy][ivar]->Draw("hist same");

      hUnfNorm_Clone[iy][ivar]->SetTitle("");
      hUnfNorm_Clone[iy][ivar]->GetYaxis()->SetTitle("Norm #frac{Other year}{2016}");
      hUnfNorm_Clone[iy][ivar]->GetYaxis()->SetTitleSize(14);
      hUnfNorm_Clone[iy][ivar]->GetYaxis()->SetTitleFont(43);
      hUnfNorm_Clone[iy][ivar]->GetYaxis()->SetTitleOffset(1.55);
      hUnfNorm_Clone[iy][ivar]->GetYaxis()->SetLabelFont(43);
      hUnfNorm_Clone[iy][ivar]->GetYaxis()->SetLabelSize(15);
      hUnfNorm_Clone[iy][ivar]->GetXaxis()->SetLabelSize(0.09);
      hUnfNorm_Clone[iy][ivar]->GetXaxis()->SetTitleSize(0.09);
      hUnfNorm_Clone[iy][ivar]->GetYaxis()->SetRangeUser(0,2);

      closure_pad1_norm[ivar]->cd();
      closure_pad1_norm[ivar]->SetBottomMargin(0.005);

      hUnfNorm[iy][ivar]->Draw("same");
    	legNorm[ivar]->Draw();


    }//------end of loop on years

    can[ivar]->Print(TString::Format("ComparisonAllYears/%s/DiffCrossSection_%s%s.pdf", varParton.Data(), variable[ivar].Data(), unfMethodStr.Data()), "pdf");
    canNorm[ivar]->Print(TString::Format("ComparisonAllYears/%s/NormDiffCrossSection_%s%s.pdf",varParton.Data(), variable[ivar].Data(), unfMethodStr.Data()), "pdf");
    //break;
  }//------ end of loop on nvars


}
