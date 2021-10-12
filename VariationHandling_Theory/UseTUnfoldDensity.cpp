#include "TH1F.h"
#include "TH2F.h"

#include "TUnfoldBinning.h"
#include "TUnfoldDensity.h"

/*
 * Perform unfolding with TUnfoldDensity
*/
TH1 *UnfoldDensity(TH2F *responseMatrix, TH1F *recoHistogram, TString varReco, TString varParton)
{
  //Construct reco binning
  const int leadingRecoBins = responseMatrix->ProjectionY()->GetNbinsX();
  Double_t BND_reco[leadingRecoBins + 1];
  for (int i = 1; i <= leadingRecoBins + 1; i++)
  {
    BND_reco[i - 1] = responseMatrix->ProjectionY()->GetBinLowEdge(i);
  }
  TUnfoldBinning *reco = new TUnfoldBinning("reco");
  reco->AddAxis(varReco, leadingRecoBins, BND_reco,
                false, true);

  //Construct parton binning
  const int leadingPartonBins = responseMatrix->ProjectionX()->GetNbinsX();
  Double_t BND_parton[leadingPartonBins + 1];
  for (int i = 1; i <= leadingPartonBins + 1; i++)
  {
    BND_parton[i - 1] = responseMatrix->ProjectionX()->GetBinLowEdge(i);
  }
  TUnfoldBinning *parton = new TUnfoldBinning("parton");
  parton->AddAxis(varParton, leadingPartonBins, BND_parton,
                  false, true);

  //Create bkg th1 with overflow bin content
  TH1 *recoBkg = parton->CreateHistogram("recoBkg", true, 0, "recoBkg", "");
  Float_t fakes = 0.;
  for (int i = 1; i <= responseMatrix->GetNbinsY() + 1; i++)
  {
    fakes += responseMatrix->GetBinContent(responseMatrix->GetNbinsX() + 1, i);
    responseMatrix->SetBinContent(responseMatrix->GetNbinsX() + 1, i, 0);
  }

  for (int i = 1; i <= responseMatrix->GetNbinsX() + 1; i++)
  {
    fakes += responseMatrix->GetBinContent(i, responseMatrix->GetNbinsY() + 1);
    responseMatrix->SetBinContent(i, responseMatrix->GetNbinsX() + 1, 0);
  }
  recoBkg->SetBinContent(responseMatrix->GetNbinsX() + 1, fakes);

  TUnfoldDensity unfold(responseMatrix, TUnfold::kHistMapOutputHoriz,
                        TUnfold::kRegModeCurvature, TUnfold::kEConstraintNone, // TUnfold::kEConstraintArea,
                        TUnfoldDensity::kDensityModeNone,
                        parton, reco);

  std::cout << unfold.SetInput(recoHistogram) << std::endl;

  unfold.SubtractBackground(recoBkg, "BGR", 1., 0.);
  //Int_t iBest = unfold.ScanLcurve(nScan, 0., 0., &lCurve, &logTauX, &logTauY);
  std::cout << "Unfold result: " << unfold.DoUnfold(0) << std::endl;

  TH1 *unfoldedHistogram = unfold.GetOutput("unfoldedHistogram", 0, 0, 0, kTRUE);
  unfoldedHistogram->SetTitle("unfoldedHistogram");
  return unfoldedHistogram;
}

void UseTUnfoldDensity()
{
  TH2F *responseMatrix;
  TH1F *signal;
  //TH1F *unfoldedHistogram = (TH1F *)UnfoldDensity(responseMatrix, signal, variable, variableParton);
}
