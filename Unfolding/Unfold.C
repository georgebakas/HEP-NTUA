/*
  For Unfolding we need the files:
  1. Signal Extracted S_i(xReco) from SignalExtraction.cpp in MassFit folder 
    MassFit/year/FiducialMeasurements/ABCDMethod/free_eb/SignalHistograms_ABCDMethod_freeEb.root
    we need to rebin the histogram taken from the S_i because it has bins same as parton
  2. Acceptance from files in ../ResponseMatrices/year/UnequalBins/ResponsesEfficiency_year.root more bins for acc
  3. Response matrix from ../ResponseMatrices/year/UnequalBins/ResponsesEfficiency_year.root 
  4. Efficiency from same files but with fewer bins
  5. Lumi taken from this file for every year

*/

#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TGraphErrors.h"

using namespace std;

using std::cin;
using std::cout;
using std::endl;

#include "TemplateConstantsUnfold.h"

TH1F *getRebinned(TH1F *h, float BND[], int N);
void SignalExtractionSpecific(TString year = "2016", TString variable = "jetPt0", bool free_eb = true, bool ABCDMethod = true);

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


void Unfold(TString year = "2016")
{

  std::vector< std::vector <Float_t> > const BND_reco = {{1000, 1100,1200,1300, 1400,1500, 1600,1700, 1800,1900, 2000,2200, 2400,2600, 2800,3000, 3200,3600, 4000,4500, 5000}, //mjj 21
                                                        {0,30,60,105,150,225,300,375,450,525,600,675,750,850,950,1025,1100,1200,1300}, //ptjj 19
                                                        {-2.4,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.4}, //yjj
                                                        {400,425,450,475,500,535,570,610,650,700,750,800,850,900,950,1025,1100,1200,1300,1400,1500}, //jetPt0 21
                                                        {400,425,450,475,500,535,570,610,650,700,750,800,850,900,950,1025,1100,1200,1300,1400,1500}}; //jetPt1 $
                                                        //{0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4}, //jetY0 25
                                                        //{0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4}}; //jetY1 25
  

  TFile *signalFile = TFile::Open(TString::Format("../MassFit/%s/FiducialMeasurement/UnequalBinning/ABCDMethod/free_eb/SignalHistograms_ABCDMethod_freeEb.root", 
                                  year.Data()));    
  int NBINS[BND_reco.size()];
  const int NVAR = 7;
  for (int i = 0; i<BND_reco.size(); i++) NBINS[i] = BND_reco[i].size()-1;

  TH1F *inSig[BND_reco.size()], *hSig[BND_reco.size()];
  TString variable[NVAR] = {"mJJ", "ptJJ", "yJJ", "jetPt0", "jetPt1","jetY0", "jetY1"};

  TCanvas *can[BND_reco.size()];
  for(int ivar =0; ivar<BND_reco.size(); ivar++)
  {
    int sizeBins = NBINS[ivar];
    float tempBND[NBINS[ivar]+1];
    std::copy(BND_reco[ivar].begin(), BND_reco[ivar].end(), tempBND);
    
    //from signal file get the initial S_i and then rebin it
    hSig[ivar] = (TH1F*)signalFile->Get(TString::Format("hSignal_%s",variable[ivar].Data()));  
    //now rebin it
    can[ivar] = new TCanvas(TString::Format("can_%d",ivar), TString::Format("can_%d",ivar),800,600);
    //inSig[ivar]->Draw();
    hSig[ivar]->Draw();
  }
  

  
  









}