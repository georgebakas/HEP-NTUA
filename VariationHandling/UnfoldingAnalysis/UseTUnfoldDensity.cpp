#include "TH1F.h"
#include "TH2F.h"

#include "TUnfoldBinning.h"
#include "TUnfoldDensity.h"

/*

In TUnfoldDensity, two methods are implemented to determine tau**2

1. ScanLcurve() locate the tau where the L-curve plot has a "kink" this function is implemented in the TUnfold class
2. ScanTau() finds the solution such that some variable (e.g. global correlation coefficient) is minimized. This function is implemented in the TUnfoldDensity class

Each of the algorithms has its own advantages and disadvantages. The algorithm (1) does not work if the input data are too similar to the MC prediction. 
Typical ** no-go ** cases of the L-curve scan are:

* the number of measurements is too small (e.g. ny=nx)
* the input data have no statistical fluctuations [identical MC events are used to fill the matrix of migrations and the vector y for a "closure test"]

The algorithm (2) only works if the variable does have a real minimum as a function of tau. 
If global correlations are minimized, the situation is as follows: 
* The matrix of migration typically introduces negative correlations. 
* The area constraint introduces some positive correlation. 
* Regularisation on the "size" introduces no correlation. 
* Regularisation on 1st or 2nd derivatives adds positive correlations.
* For these reasons, "size" regularisation does not work well with the tau-scan: the higher tau, the smaller rho, but there is no minimum. 
* As a result, large values of tau (too strong regularisation) are found. 
* In contrast, the tau-scan is expected to work better with 1st or 2nd derivative regularisation, because at some point the negative correlations from migrations are approximately cancelled by the positive correlations from the regularisation conditions.

whichever algorithm is used, the output has to be checked:

* The L-curve should have approximate L-shape and the final choice of tau should not be at the very edge of the scanned region
* The scan result should have a well-defined minimum and the final choice of tau should sit right in the minimum

Comments:
EDensityMode:
  - choice of regularisation scale factors to cinstruct the matrix L --> curvature matrix 
  kDensityModeNone 	
  no scale factors, matrix L is similar to unity matrix

  kDensityModeBinWidth 	
  scale factors from multidimensional bin width

  kDensityModeUser 	
  scale factors from user function in TUnfoldBinning

  kDensityModeBinWidthAndUser 	
  scale factors from multidimensional bin width and user function

EScanTauMode
  - scan mode for correlation scan

  kEScanTauRhoAvg 	
  average global correlation coefficient (from TUnfold::GetRhoI())

  kEScanTauRhoMax 	
  maximum global correlation coefficient (from TUnfold::GetRhoI())

  kEScanTauRhoAvgSys 	
  average global correlation coefficient (from TUnfoldSys::GetRhoItotal())

  kEScanTauRhoMaxSys 	
  maximum global correlation coefficient (from TUnfoldSys::GetRhoItotal())

  kEScanTauRhoSquareAvg 	
  average global correlation coefficient squared (from TUnfold::GetRhoI())

  kEScanTauRhoSquareAvgSys 	
  average global correlation coefficient squared (from TUnfoldSys::GetRhoItotal())





*/



/*
 * Perform unfolding with TUnfoldDensity
*/
TH1 *UnfoldDensity_complex(TH2F *responseMatrix, TH1F *recoHistogram, TString varReco, TString varParton)
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

  /*
  kRegModeNone: no regularisation, or defined later by RegularizeXXX() methods
  kRegModeSize: regularise the amplitude of the output distribution
  kRegModeDerivative: regularize the 1st derivative of the output distribution
  kRegModeCurvature: regularize the 2nd derivative of the output distribution
  kRegModeMixed: mixed regularisation pattern
  */
  // preserve the area
  TUnfold::EConstraint constraintMode= TUnfold::kEConstraintArea; // TUnfold::kEConstraintNone

  // basic choice of regularisation scheme:
  //    curvature (second derivative)
  TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature; // TUnfold::kRegModeCurvature

  // density flags
  TUnfoldDensity::EDensityMode densityFlags= TUnfoldDensity::kDensityModeBinWidth; //TUnfoldDensity::kDensityModeNone
  // detailed steering for regularisation
  const char *REGULARISATION_DISTRIBUTION=0;
  const char *REGULARISATION_AXISSTEERING=0;

  /*
    axisSteering is a string with several tokens, separated by a semicolon: "axisName[options];axisName[options];...".
    axisName: the name of an axis. The special name * matches all. So the argument distribution selects one (or all) distributions. Within the selected distribution(s), steering options may be specified for each axis (or for all axes) to define the regularisation conditions.
    options one or several character as follows:
    u : exclude underflow bin from derivatives along this axis
    o : exclude overflow bin from derivatives along this axis
    U : exclude underflow bin
    O : exclude overflow bin
    b : use bin width for derivative calculation
    B : same as 'b', in addition normalize to average bin width
    N : completely exclude derivatives along this axis
    p : axis is periodic (e.g. azimuthal angle), so include derivatives built from combinations involving bins at both ends of the axis "wrap around"
    example: axisSteering="*[UOB]" uses bin widths to calculate derivatives but underflow/overflow bins are not regularized
  */
  TUnfoldDensity unfold(responseMatrix, TUnfold::kHistMapOutputHoriz,
                        regMode, constraintMode, densityFlags,
                        parton, reco, 
                        REGULARISATION_DISTRIBUTION, REGULARISATION_AXISSTEERING);

  std::cout << unfold.SetInput(recoHistogram) << std::endl;
  unfold.SubtractBackground(recoBkg, "BGR", 1., 0.);

  // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  TH2 *histL= unfold.GetL("L");
  for(Int_t j=1;j<=histL->GetNbinsY();j++) {
      cout<<"L["<<unfold.GetLBinning()->GetBinName(j)<<"]";
      for(Int_t i=1;i<=histL->GetNbinsX();i++) {
          Double_t c=histL->GetBinContent(i,j);
          if(c!=0.0) cout<<" ["<<i<<"]="<<c;
      }
      cout<<"\n";
  }

  // run the unfolding
  // here, tau is determined by scanning the global correlation coefficients
  Int_t nScan=30;
  TSpline *rhoLogTau=0;
  TGraph *lCurve=0;
  // for determining tau, scan the correlation coefficients
  // correlation coefficients may be probed for all distributions
  // or only for selected distributions
  // underflow/overflow bins may be included/excluded
  //
  const char *SCAN_DISTRIBUTION=0;
  // example: axisSteering="*[UOB]" uses bin widths to calculate derivatives but underflow/overflow bins are not regularized
  const char *SCAN_AXISSTEERING=0;
  cout<<"trying to fetch iBest.. "<<endl;
  Int_t iBest=unfold.ScanTau(nScan,0.,1.,&rhoLogTau,
                            TUnfoldDensity::kEScanTauRhoAvgSys,
                            SCAN_DISTRIBUTION,SCAN_AXISSTEERING,
                            &lCurve);
  cout<<"iBEST: "<< iBest<<endl;
  // create graphs with one point to visualize best choice of tau
  Double_t t[1],rho[1],x[1],y[1];
  rhoLogTau->GetKnot(iBest,t[0],rho[0]);
  lCurve->GetPoint(iBest,x[0],y[0]);
  TGraph *bestRhoLogTau=new TGraph(1,t,rho);
  TGraph *bestLCurve=new TGraph(1,x,y);
  Double_t *tAll=new Double_t[nScan],*rhoAll=new Double_t[nScan];
  for(Int_t i=0;i<nScan;i++) {
    rhoLogTau->GetKnot(i,tAll[i],rhoAll[i]);
  }
  TGraph *knots=new TGraph(nScan,tAll,rhoAll);

  cout<<"chi**2="<<unfold.GetChi2A()<<"+"<<unfold.GetChi2L()
      <<" / "<<unfold.GetNdf()<<"\n";

  TCanvas *c = new TCanvas("c1", "c1", 800, 600);
  rhoLogTau->Draw();
  // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Int_t iBest = unfold.ScanLcurve(nScan, 0., 0., &lCurve, &logTauX, &logTauY);
  // std::cout << "Unfold result: " << unfold.DoUnfold(iBest) << std::endl;
  
  //=======================================================================================================================================
  // Step 4: retrieve and plot unfolding results
 
  TH1 *unfoldedHistogram = unfold.GetOutput("unfoldedHistogram_complex", 0, 0, 0, kTRUE);
  unfoldedHistogram->SetTitle("unfoldedHistogram_complex");

  // get matrix of probabilities
  TH2 *histProbability=unfold.GetProbabilityMatrix("histProbability");
  // get global correlation coefficients
  /*TH1 *histGlobalCorr = unfold.GetRhoItotal("histGlobalCorr",0,0,0,kFALSE);*/
  TH1 *histGlobalCorrScan = unfold.GetRhoItotal("histGlobalCorrScan",0,SCAN_DISTRIBUTION,SCAN_AXISSTEERING,kFALSE);
  TH2 *histCorrCoeff = unfold.GetRhoIJtotal("histCorrCoeff",0,0,0,kFALSE);

  // get covariance matrix 
  TH2 *covMatrix = new TH2D("cov matrix", "cov matrix", responseMatrix->GetNbinsX(), 0, responseMatrix->GetNbinsX(), responseMatrix->GetNbinsX(), 0, responseMatrix->GetNbinsX());
  unfold.GetEmatrix(covMatrix);
  
  //=======================================================================================================================================
  

  TCanvas canvas;
  canvas.Print("complex_unfolding.ps[");
 
  //========== page 1 ============
  // tau-scan, global correlations, correlation coefficients
  canvas.Clear();
  canvas.Divide(3,2);
 
  // (1) matrix of probabilities
  canvas.cd(1);
  // gStyle->SetPalette(kOcean);
  histProbability->Draw("colz 1");

  // (2) scan of correlation vs tau
  canvas.cd(2);
  rhoLogTau->Draw();
  knots->Draw("*");
  bestRhoLogTau->SetMarkerColor(kRed);
  bestRhoLogTau->Draw("*");
  
  // (3) global correlation coefficients for the distributions
  //     used during the scan
  canvas.cd(3);
  // gStyle->SetPalette(kOcean);
  histCorrCoeff->Draw("colz1");


  // (4) global correlation coefficients for the distributions
  //     used during the scan
  canvas.cd(4);
  // gStyle->SetPalette(kOcean);
  histGlobalCorrScan->Draw("HIST");
  
  // (5) L-curve
  canvas.cd(5);
  lCurve->Draw("AL");
  bestLCurve->SetMarkerColor(kRed);
  bestLCurve->Draw("*");

  // (6) covariance matrix 
  canvas.cd(6);
  // canvas.SetLogz();
  // gStyle->SetLogz();
  covMatrix->Draw("COLZ 1");

  canvas.Print("complex_unfolding.ps");

  canvas.Print("complex_unfolding.ps]");
  

  TCanvas *can = new TCanvas("test", "test", 800, 600);
  can->SetLogz();
  covMatrix->Draw("COLZ 1");


  return unfoldedHistogram;
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


TH1 *UnfoldDensity_regular(TH2F *responseMatrix, TH1F *recoHistogram, TString varReco, TString varParton)
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
  
  //=======================================================================================================================================
  // Step 4: retrieve and plot unfolding results

  // get matrix of probabilities
  TH2 *histProbability=unfold.GetProbabilityMatrix("histProbability");
  // get global correlation coefficients
  const char *SCAN_DISTRIBUTION=0;
  // example: axisSteering="*[UOB]" uses bin widths to calculate derivatives but underflow/overflow bins are not regularized
  const char *SCAN_AXISSTEERING="*[UOB]";
  TH1 *histGlobalCorrScan = unfold.GetRhoItotal("histGlobalCorrScan",0,SCAN_DISTRIBUTION,SCAN_AXISSTEERING,kFALSE);
  TH2 *histCorrCoeff = unfold.GetRhoIJtotal("histCorrCoeff",0,0,0,kFALSE);

  // get covariance matrix 
  TH2 *covMatrix = new TH2D("cov matrix", "cov matrix", responseMatrix->GetNbinsX(), 0, responseMatrix->GetNbinsX(), responseMatrix->GetNbinsX(), 0, responseMatrix->GetNbinsX());
  unfold.GetEmatrix(covMatrix);
  
  //=======================================================================================================================================
  

  TCanvas canvas;
  canvas.Print("regular_unfolding.ps[");
 
  //========== page 1 ============
  // tau-scan, global correlations, correlation coefficients
  canvas.Clear();
  canvas.Divide(2,2);

  // (1) matrix of probabilities
  canvas.cd(1);
  // gStyle->SetPalette(kOcean);
  histProbability->Draw("colz 1");

  
  // (3) global correlation coefficients for the distributions
  //     used during the scan
  canvas.cd(2);
  // gStyle->SetPalette(kOcean);
  histCorrCoeff->Draw("colz1");


  // (4) global correlation coefficients for the distributions
  //     used during the scan
  canvas.cd(3);
  // gStyle->SetPalette(kOcean);
  histGlobalCorrScan->Draw("HIST");

  // (6) covariance matrix 
  canvas.cd(4);
  // canvas.SetLogz();
  covMatrix->Draw("COLZ 1");

  canvas.Print("regular_unfolding.ps");

  canvas.Print("regular_unfolding.ps]");
  

  TCanvas *can = new TCanvas("test regular", "test regular", 800, 600);
  can->SetLogz();
  covMatrix->Draw("COLZ 1");


  // return histogram 
  TH1 *unfoldedHistogram = unfold.GetOutput("unfoldedHistogram_regular", 0, 0, 0, kTRUE);
  unfoldedHistogram->SetTitle("unfoldedHistogram_regular");

  return unfoldedHistogram;
}

// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


void UseTUnfoldDensity()
{
  TH2F *responseMatrix;
  TH1F *signal;
  //TH1F *unfoldedHistogram = (TH1F *)UnfoldDensity(responseMatrix, signal, variable, variableParton);
}


// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





TH1* UnfoldDensity_lcurve(TH2F *responseMatrix, TH1F *recoHistogram, TString varReco, TString varParton)
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

  /*
  kRegModeNone: no regularisation, or defined later by RegularizeXXX() methods
  kRegModeSize: regularise the amplitude of the output distribution
  kRegModeDerivative: regularize the 1st derivative of the output distribution
  kRegModeCurvature: regularize the 2nd derivative of the output distribution
  kRegModeMixed: mixed regularisation pattern
  */
  // preserve the area
  TUnfold::EConstraint constraintMode = TUnfold::kEConstraintNone; // TUnfold::kEConstraintArea

  // basic choice of regularisation scheme:
  //    curvature (second derivative)
  // curvature only sets results 
  TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature; // TUnfold::kRegModeCurvature

  // density flags
  TUnfoldDensity::EDensityMode densityFlags= TUnfoldDensity::kDensityModeBinWidth; //TUnfoldDensity::kDensityModeNone
  // detailed steering for regularisation
  const char *REGULARISATION_DISTRIBUTION=0;
  const char *REGULARISATION_AXISSTEERING=0;


  TUnfoldDensity unfold(responseMatrix, TUnfold::kHistMapOutputHoriz,
                        regMode, constraintMode, densityFlags,
                        parton, reco, 
                        REGULARISATION_DISTRIBUTION, REGULARISATION_AXISSTEERING);

  std::cout << unfold.SetInput(recoHistogram) << std::endl;
  unfold.SubtractBackground(recoBkg, "BGR", 1., 0.);

  // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  TH2 *histL= unfold.GetL("L");
  for(Int_t j=1;j<=histL->GetNbinsY();j++) {
      cout<<"L["<<unfold.GetLBinning()->GetBinName(j)<<"]";
      for(Int_t i=1;i<=histL->GetNbinsX();i++) {
          Double_t c=histL->GetBinContent(i,j);
          if(c!=0.0) cout<<" ["<<i<<"]="<<c;
      }
      cout<<"\n";
  }

  // run the unfolding
  const char *SCAN_DISTRIBUTION=0;
  // example: axisSteering="*[UOB]" uses bin widths to calculate derivatives but underflow/overflow bins are not regularized
  const char *SCAN_AXISSTEERING=0;
  cout<<"trying to fetch iBest.. "<<endl;
  
  // here, tau is determined by scanning the global correlation coefficients
  Int_t nScan=100;
  TGraph *lCurve=0;
  TSpline *logTauX, *logTauY;
  TSpline *logTauCurvature;
  // return value: the coordinate number in the logTauX,logTauY graphs corresponding to the "final" choice of tau

  // Recommendation: always check logTauCurvature, it should be a peaked function (similar to a Gaussian), the maximum corresponding to the final choice of tau. 
  // Also, check the lCurve it should be approximately L-shaped. If in doubt, adjust tauMin and tauMax until the results are satisfactory.
  
  Int_t iBest=unfold.ScanLcurve(nScan, 0, 1, &lCurve, &logTauX, &logTauY, &logTauCurvature);
  
  cout<<"iBEST: "<< iBest<<endl;
  
  // create graphs with one point to visualize best choice of tau
  TCanvas *cL = new TCanvas("lCurve", "lCurve", 800, 600);
  lCurve->Draw();

  // create graphs with one point to visualize best choice of tau
  TCanvas *cX = new TCanvas("lcurve cx", "lcurve cx", 800, 600);
  logTauX->Draw();

  TCanvas *cy = new TCanvas("lcurve cy", "lcurve cy", 800, 600);
  logTauY->Draw();

  TCanvas *curve = new TCanvas("lcurve curvature", "lcurve curvature", 800, 600);
  logTauCurvature->Draw();  

  // create graphs for some info 
  Double_t x[1], y[1];
  lCurve->GetPoint(iBest,x[0],y[0]);
  TGraph *bestLCurve=new TGraph(1,x,y);

  Double_t xC[1], yC[1];
  logTauCurvature->GetKnot(iBest,xC[0],yC[0]);
  TGraph *bestCurvature=new TGraph(1,xC,yC);
  
  Double_t xt, yt;
  std::cout << "Unfold L-curve result: " << lCurve->GetPoint(iBest, xt, yt)<< std::endl;
  std::cout << "logTauX result: " << xt << std::endl;
  std::cout << "logTauY result: " << yt << std::endl;
  logTauX->GetKnot(iBest, xt, yt);
  std::cout << "Log tau result: " << xt << std::endl;

  logTauY->GetKnot(iBest, xt, yt);
  std::cout << "Log Tau result: " << xt << std::endl;
  std::cout << "bestTau result: " << unfold.GetTau() << std::endl;
  
  // unfold with best tau from l-curve
  // unfold.DoUnfold(unfold.GetTau());
  
  //=======================================================================================================================================
  // Step 4: retrieve and plot unfolding results
 
  TH1 *unfoldedHistogram = unfold.GetOutput("unfoldedHistogram_lcurve", 0, 0, 0, kTRUE);
  unfoldedHistogram->SetTitle("unfoldedHistogram_lcurve");

  // get matrix of probabilities
  TH2 *histProbability=unfold.GetProbabilityMatrix("histProbability");
  // get global correlation coefficients
  TH1 *histGlobalCorrScan = unfold.GetRhoItotal("histGlobalCorrScan",0,SCAN_DISTRIBUTION,SCAN_AXISSTEERING,kFALSE);
  TH2 *histCorrCoeff = unfold.GetRhoIJtotal("histCorrCoeff",0,0,0,kFALSE);

  // get covariance matrix 
  TH2 *covMatrix = new TH2D("cov matrix", "cov matrix", responseMatrix->GetNbinsX(), 0, responseMatrix->GetNbinsX(), responseMatrix->GetNbinsX(), 0, responseMatrix->GetNbinsX());
  unfold.GetEmatrix(covMatrix);
  
  //=======================================================================================================================================
  
  TCanvas canvas;
  canvas.Print("lcurve_unfolding.ps[");

  //========== page 1 ============
  // tau-scan, global correlations, correlation coefficients
  canvas.Clear();
  canvas.Divide(3,2);
 
  // (1) matrix of probabilities
  canvas.cd(1);
  histProbability->Draw("colz 1");

  // (2) scan of L-curve
  canvas.cd(2);
  lCurve->Draw("AL");
  bestLCurve->SetMarkerColor(kRed);
  bestLCurve->Draw("same *");
  
  // (3) global correlation coefficients for the distributions
  //     used during the scan
  canvas.cd(3);
  histCorrCoeff->Draw("colz1");


  // (4) global correlation coefficients for the distributions
  //     used during the scan
  canvas.cd(4);
  histGlobalCorrScan->Draw("HIST");
  
  // (5) logTauCurvature
  canvas.cd(5);
  logTauCurvature->Draw();
  bestCurvature->SetMarkerColor(kRed);
  bestCurvature->Draw("same *");

  // (6) covariance matrix 
  canvas.cd(6);
  covMatrix->Draw("COLZ 1");

  canvas.Print("lcurve_unfolding.ps");

  canvas.Print("lcurve_unfolding.ps]");
  

  TCanvas *can = new TCanvas("test lcurve_unfolding", "test lcurve_unfolding", 800, 600);
  can->SetLogz();
  covMatrix->Draw("COLZ 1");


  return unfoldedHistogram; 
}