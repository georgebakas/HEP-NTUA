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

#include "FinalResultsConstants.h"


#include "../CMS_plots/CMS_lumi.C"
#include "../CMS_plots/tdrstyle.C"
#include "TemplateConstants.h"



void CompareResultsAllYearsToTOP18013(bool isParton = true, bool normalized = false)
{
  gStyle->SetOptStat(0);
  initFilesMapping();

  TString partonParticle = "Parton";
  TString varParton_top = partonParticle;
  if(!isParton)
  {
    partonParticle = "Particle";
    varParton_top = "Gen";

  }

  TString baseDir = "/Users/georgebakas/Documents/HEP-NTUA_ul/VariationHandling";
  AnalysisConstants::initConstants();

  const int NVAR = 7;
  TString variables[NVAR] = {"mJJ", "ptJJ", "yJJ", "jetPt0", "jetPt1","jetY0", "jetY1"};
  TString variableGen[NVAR] = {"mJJGen", "ptJJGen", "yJJGen", "genjetPt0", "genjetPt1","genjetY0", "genjetY1"};
  TString variableParton[NVAR] = {"mTTbarParton", "ptTTbarParton", "yTTbarParton", "partonPt0", "partonPt1","partonY0", "partonY1"};

  //Our file with the unfolding results:
  TFile *inf;

  TFile *nominalFile = TFile::Open(TString::Format("%s/UnfoldedCombined/Nominal/OutputFile%s.root",
                                                    baseDir.Data(),
                                                    partonParticle.Data()));

  TString outputDir = TString::Format("%s/ComparisonWithTOP18013/",
                                        baseDir.Data());

  for (unsigned int v = 0; v < NVAR; v++)
  {  
    TString variable = AnalysisConstants::unfoldingVariables[v];
    // Giannis and George results 
    TH1F *nominalHistogram = (TH1F *)nominalFile->Get(TString::Format("hUnfold%s_%s",
                                                                        (normalized ? "Norm" : "Final"), 
                                                                        variable.Data()));

    TH1F *finalResult = (TH1F *)nominalHistogram->Clone(TString::Format("FinalResult%s_%s",
                                                                        (normalized ? "Norm" : "Final"), 
                                                                        variable.Data()));

    TH1F *theoryHistogram = (TH1F *)nominalFile->Get(TString::Format("hTheory%s_%s",
                                                                        (normalized ? "Norm" : ""),
                                                                        variable.Data()));

    TH1F *finalTheory = (TH1F *)theoryHistogram->Clone(TString::Format("FinalTheory%s_%s",
                                                                        (normalized ? "Norm" : ""),
                                                                        variable.Data()));

    //top-18-013 files and histograms
    TFile *top18013file = TFile::Open(TString::Format("../Unfolding/Results-TOP18013/CrossSection_%s_%s.root",varParton_top.Data(), variable.Data()));
    //top-18-013 data
    TH1F *hTOP18013 = (TH1F*)top18013file->Get(TString::Format("%sCrossSection_%s_Nominal",
                                                              (normalized ? "Norm" : ""),
                                                              varParton_top.Data()));

    TH1F *hTOP18013_Clone = (TH1F*)hTOP18013->Clone(TString::Format("Top-18-013Result%s_%s",
                                                              (normalized ? "Norm" : ""),
                                                              variable.Data()));

    //top-18-013 theory plots:
    TH1F *hTOP18013_theory = (TH1F*)top18013file->Get(TString::Format("%sCrossSection_PowhegPythia8",
                                                              (normalized ? "Norm" : "")));

    TH1F *hTOP18013_theoryClone = (TH1F*)hTOP18013->Clone(TString::Format("Top-18-013Theory%s_%s",
                                                              (normalized ? "Norm" : ""),
                                                              variable.Data()));

    TCanvas *can;
    TLegend *leg;
    TPad *closure_padRatio, *closure_pad1;
    TH1F *hTemp, *hTemp_run2;

    cout<<variable<<endl;
    can = new TCanvas(TString::Format("can_%s",variable.Data()),TString::Format("can_%s",variable.Data()) , 800,600);
    can->cd();
    //closure_pad1 =  new TPad(TString::Format("closure_pad1%d",v),TString::Format("closure_pad1%d", v),0.,0.,1.,1.);
    //closure_pad1->Draw();
    //closure_padRatio = new TPad(TString::Format("closure_pad2%d",v),TString::Format("closure_pad2%d", v),0.,0.,1.,0.3);
    //closure_padRatio->Draw();
    leg = new TLegend(0.65,0.7,0.9,0.9);

    // kkousour plots
    hTOP18013->SetLineColor(kBlack);
    hTOP18013->SetMarkerColor(kBlack);
    hTOP18013->SetMarkerStyle(20);

    //get our plots:
    nominalHistogram->SetLineColor(kRed);
    nominalHistogram->SetMarkerColor(kRed);
    nominalHistogram->SetMarkerStyle(20);

    if(variables[v].Contains("jetY"))
    {
      if(isParton)
      {
        nominalHistogram->GetXaxis()->SetTitle("|"+variableParton[v]+"|");
        hTOP18013->GetXaxis()->SetTitle("|"+variableParton[v]+"|");
      }
      else
      {
        hTOP18013->GetXaxis()->SetTitle("|"+variableGen[v]+"|");
        nominalHistogram->GetXaxis()->SetTitle("|"+variableGen[v]+"|");
      }
        
    }

    //now divide with theory 
    nominalHistogram->Divide(theoryHistogram);

    //divide top-18-013 histograms
    hTOP18013->Divide(hTOP18013_theory);
    

    //draw the unfolded and extrapolated with the mc result
    can->cd();
    //closure_pad1->SetBottomMargin(0.005);
    //closure_pad1->cd();

    hTOP18013->SetTitle("Unfolded Ratio Comparison (TOP-18-013)");
    hTOP18013->GetYaxis()->SetTitle("#frac{Data}{Nom. Theory}");
    hTOP18013->GetYaxis()->SetRangeUser(0.2,1.8);
    if (partonParticle.EqualTo("Parton"))
        hTOP18013->GetXaxis()->SetTitle(AnalysisConstants::partonAxisTitles[v]);
    else
        hTOP18013->GetXaxis()->SetTitle(AnalysisConstants::particleAxisTitles[v]);
    leg->AddEntry(hTOP18013, "TOP-18-013 2016", "lpe");
    leg->AddEntry(nominalHistogram, "Combined Full Run-II", "lpe");

    hTOP18013->Draw();
    nominalHistogram->Draw("same");
    leg->Draw();

    //for the double ratio:
    //divide with 2016 top-18-013 as reference
    /*
    hTemp = (TH1F*)hTOP18013->Clone("hTemp");
    hTemp_run2 = (TH1F*)nominalHistogram->Clone("hTempFull2");
    hTemp_run2 ->Divide(hTemp_run2);

    closure_padRatio->SetTopMargin(0.05);
    closure_padRatio->SetBottomMargin(0.3);
    //closure_padRatio->SetGrid();
    closure_padRatio->cd();

    hTemp_run2->SetTitle("");
    hTemp_run2->GetYaxis()->SetTitle("#frac{Other year}{2016}");
    hTemp_run2->GetYaxis()->SetTitleSize(14);
    hTemp_run2->GetYaxis()->SetTitleFont(43);
    hTemp_run2->GetYaxis()->SetTitleOffset(1.55);
    hTemp_run2->GetYaxis()->SetLabelFont(43);
    hTemp_run2->GetYaxis()->SetLabelSize(15);
    hTemp_run2->GetXaxis()->SetLabelSize(0.09);
    hTemp_run2->GetXaxis()->SetTitleSize(0.09);
    hTemp_run2->GetYaxis()->SetRangeUser(0,2);

    hTemp_run2->Draw("hist same");
    */
  
    can->SaveAs(TString::Format("%s/%s%s.pdf",
                                outputDir.Data(),
                                variable.Data(),
                                (normalized ? "_normalized" : "")),
                                "pdf");
  }


}
