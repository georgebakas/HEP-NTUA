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
#include "TUnfold.h"

using namespace std;

using std::cin;
using std::cout;
using std::endl;

#include "TemplateConstantsUnfold.h"

TH1F *getRebinned(TH1F *h, float BND[], int N);
TH1 *unfoldedOutput(TH2F *hResponse_, TH1F *hReco, float BND[], int sizeBins);

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


void Unfold(TString year = "2016", bool isParton = true)
{
  initFilesMapping();
  std::vector< std::vector <Float_t> > const BND_reco = {{1000, 1100,1200,1300, 1400,1500, 1600,1700, 1800,1900, 2000,2200, 2400,2600, 2800,3000, 3200,3600, 4000,4500, 5000}, //mjj 21
                                                        {0,30,60,105,150,225,300,375,450,525,600,675,750,850,950,1025,1100,1200,1300}, //ptjj 19
                                                        {-2.4,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.4}, //yjj
                                                        {400,425,450,475,500,535,570,610,650,700,750,800,850,900,950,1025,1100,1200,1300,1400,1500}, //jetPt0 21
                                                        {400,425,450,475,500,535,570,610,650,700,750,800,850,900,950,1025,1100,1200,1300,1400,1500}}; //jetPt1 $
                                                        //{0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4}, //jetY0 25
                                                        //{0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4}}; //jetY1 25
  std::vector< std::vector <Float_t> > const BND_gen = {{1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, 4000, 5000}, //mjj
                                                        {0,60,150,300,450,600,750,950,1100,1300}, //ptjj
                                                        {-2.4,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.4}, //yjj
                                                        {400,450,500,570,650,750,850,950,1100,1300,1500}, //jetPt0     
                                                        {400,450,500,570,650,750,850,950,1100,1300,1500}}; //jetPt1
                                                        //{0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}, //jetY0
                                                        //{0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}}; //jetY1

  float LUMI = luminosity[year];
  TFile *signalFile = TFile::Open(TString::Format("../MassFit/%s/FiducialMeasurement/UnequalBinning/ABCDMethod/free_eb/SignalHistograms_ABCDMethod_freeEb.root", 
                                  year.Data()));    

  TFile *effAccInf = TFile::Open(TString::Format("../ResponseMatrices/%s/UnequalBins/ResponsesEfficiency_%s.root", year.Data(), year.Data()));

  //whether parton or particle
  TString varParton = "Parton";
  if(!isParton) varParton = "Particle";

  int NBINS[BND_reco.size()];
  int NBINS_GEN[BND_gen.size()];
  const int NVAR = 7;
  for (int i = 0; i<BND_reco.size(); i++) NBINS[i] = BND_reco[i].size()-1;
  for (int i = 0; i<BND_gen.size(); i++)  NBINS_GEN[i] = BND_gen[i].size()-1;

  TH1F *inSig[BND_reco.size()], *hSig[BND_reco.size()];
  TString variable[NVAR] = {"mJJ", "ptJJ", "yJJ", "jetPt0", "jetPt1","jetY0", "jetY1"};
  TString variableGen[NVAR] = {"mJJGen", "ptJJGen", "yJJGen", "genjetPt0", "genjetPt1","genjetYt0", "genjetY1"};
  TString variableParton[NVAR] = {"mTTbarParton", "ptTTbarParton", "yTTbarParton", "partonPt0", "partonPt1","partonY0", "partonY1"};    

  TH2F *hResponse[BND_reco.size()];
  TUnfold *unf[BND_reco.size()];
  TH1 *hUnf[BND_reco.size()];
  TH1F *hUnf_Clone[BND_reco.size()];

  TCanvas *can[BND_reco.size()];
  for(int ivar =0; ivar<BND_reco.size(); ivar++)
  {
    int sizeBins = NBINS[ivar];
    float tempBND[NBINS[ivar]+1];
    std::copy(BND_reco[ivar].begin(), BND_reco[ivar].end(), tempBND);

    float tempBNDGen[NBINS_GEN[ivar]+1];
    std::copy(BND_gen[ivar].begin(), BND_gen[ivar].end(), tempBNDGen);
    cout<<"------"<<endl;
    cout<<sizeBins<<endl;
    cout<<NBINS_GEN[ivar]<<endl; 
    //from signal file get the initial S_j with j bins ~ 2* parton bins (i)
    hSig[ivar] = (TH1F*)signalFile->Get(TString::Format("hSignal_%s",variable[ivar].Data()));  
    can[ivar] = new TCanvas(TString::Format("can_%d",ivar), TString::Format("can_%d",ivar),800,600);
    //inSig[ivar]->Draw();
    //hSig[ivar]->Draw();
    
    //set the new content and get acceptance
    TEfficiency *acceptance =  (TEfficiency*)effAccInf->Get(TString::Format("Acceptance%s_%s",varParton.Data(), variable[ivar].Data()));
    for(int j =1; j<=hSig[ivar]->GetNbinsX(); j++)
    {
      float acc = acceptance->GetEfficiency(j);
      hSig[ivar]->SetBinContent(j, acc*hSig[ivar]->GetBinContent(j));
      //handle errors as well--> asymmetric error from acceptance
      float accError = (acceptance->GetEfficiencyErrorLow(j) + acceptance->GetEfficiencyErrorUp(j))/2;
      hSig[ivar]->SetBinError(j, accError*hSig[ivar]->GetBinError(j));
    }

    TString tempVar;
    if(isParton) 
      tempVar = variableParton[ivar];
    else
      tempVar = variableGen[ivar];

    //get response matrix
    hResponse[ivar] = (TH2F*)effAccInf->Get(TString::Format("h%sResponse_%s",varParton.Data(), variable[ivar].Data()));
   
    //this will be used to unfold result
    cout<<"entering unfolding method!"<<endl;
    hUnf[ivar] = new TH1F(TString::Format("hUnf_%s", variable[ivar].Data()), TString::Format("hUnf_%s", variable[ivar].Data()), NBINS_GEN[ivar], tempBNDGen);
    hUnf[ivar] = unfoldedOutput(hResponse[ivar], hSig[ivar], tempBNDGen, NBINS_GEN[ivar]);
    TString axisTitle = variable[ivar];
    if(variable[ivar].EqualTo("yJJ")) 
    	hUnf[ivar]->GetXaxis()->SetTitle(variable[ivar]);
    else 
    	hUnf[ivar]->GetXaxis()->SetTitle(TString::Format("%s [GeV]", variable[ivar].Data()));
    hUnf[ivar]->GetYaxis()->SetTitle(TString::Format("#frac{d#sigma}{d#chi} %s", varParton.Data()));
    hUnf[ivar]->GetYaxis()->SetTitleOffset(1.4);
    hUnf_Clone[ivar] = (TH1F*)hUnf[ivar]->Clone(TString::Format("hUnf_Clone %s",variable[ivar].Data()));
    can[ivar] ->SetLogy();
    cout<<"hUnf received!!"<<endl;
    cout<<variable[ivar]<<endl;
    TEfficiency *efficiency =  (TEfficiency*)effAccInf->Get(TString::Format("Efficiency%s_%s",varParton.Data(), tempVar.Data()));
    for(int i =1; i<=hUnf[ivar]->GetNbinsX(); i++)
    {
      float eff = efficiency->GetEfficiency(i);
      cout<<"i= "<<i<<" eff="<<eff<<endl;
      cout<<hUnf_Clone[ivar]->GetBinContent(i)<<endl;
      hUnf[ivar]->SetBinContent(i, hUnf[ivar]->GetBinContent(i)/eff);
      cout<<hUnf[ivar]->GetBinContent(i)<<endl;
      //handle errors as well--> asymmetric error from acceptance
      float effError = (efficiency->GetEfficiencyErrorLow(i) + efficiency->GetEfficiencyErrorUp(i))/2;
      hUnf[ivar]->SetBinError(i, hUnf[ivar]->GetBinError(i)/effError);
    }
    hUnf[ivar]->Scale(1/LUMI, "width");
    hUnf_Clone[ivar]->Scale(1/LUMI, "width");
    hUnf_Clone[ivar]->SetLineColor(kRed);
	
    //check shape
    //hUnf[ivar]->Scale(1./hUnf[ivar]->Integral());
    //hUnf_Clone[ivar]->Scale(1./hUnf_Clone[ivar]->Integral());

	hUnf[ivar]->Draw();
	hUnf_Clone[ivar]->Draw("same");
  }
  
}

TH1 *unfoldedOutput(TH2F *hResponse_, TH1F *hReco, float BND[], int sizeBins)
{	
	//TCanvas *can_response = new TCanvas("can_response", "can_response", 800,600);
	//hResponse_->Draw("BOX");
	/*
    //binning control needed here!!
    TH1F *hUnf_ = new TH1F("hUnf_", "hUnf",sizeBins, BND);
    TUnfold *unf = new TUnfold(hResponse_,TUnfold::kHistMapOutputVert, TUnfold::kRegModeSize);

    int size = hUnf_->GetSize()-2;
    int *binMap=new Int_t[size+2]; 
    for(int i=1;i<=size;i++) binMap[i]=i;
    binMap[0]=-1;
    binMap[size+1]=-1;

    float biasScale = 0.0;
    int nScan=30; //this is number of scans Double t tauMin=1.0E−9;
    double tauMax=1.0;
    double tauMin=10E-8;
    int iBest;
    TSpline *logTauX,*logTauY;
    TGraph *lCurve ;

    iBest=unf->ScanLcurve(nScan ,tauMin ,tauMax, &lCurve);
    cout<<"Best tau at: "<<iBest;
    cout<<" with tau="<<unf->GetTau()<<endl;

    unf->DoUnfold(unf->GetTau(),hReco, biasScale);
    //unf->DoUnfold(tauMin,hReco, biasScale);
    unf->GetOutput(hUnf_, binMap);

    cout<<"returning hUnf_ ..."<<endl;
    return hUnf_;/*/

    TUnfoldDensity unfold(hResponse_,TUnfold::kHistMapOutputVert);
    unfold.SetInput(hReco);
      //========================================================================
	  // the unfolding is done here
	  //
	  // scan L curve and find best point
	  Int_t nScan=30;
	  // use automatic L-curve scan: start with taumin=taumax=0.0
	  //Double_t tauMin=0.0;
	  //Double_t tauMax=0.0;
	  Double_t tauMax=1.0;
      Double_t tauMin=10E-8;
	  Int_t iBest;
	  TSpline *logTauX,*logTauY;
	  TGraph *lCurve;
	  // if required, report Info messages (for debugging the L-curve scan)
	//#ifdef VERBOSE_LCURVE_SCAN
	  //Int_t oldinfo=gErrorIgnoreLevel;
	  //gErrorIgnoreLevel=kInfo;
	//#endif
	  // this method scans the parameter tau and finds the kink in the L curve
	  // finally, the unfolding is done for the best choice of tau
	  iBest=unfold.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
	  // if required, switch to previous log-level
	//#ifdef VERBOSE_LCURVE_SCAN
	 // gErrorIgnoreLevel=oldinfo;
	//#endif
	  //==========================================================================
	  //==========================================================================
	  // print some results
	  //
	  std::cout<<"tau="<<unfold.GetTau()<<"\n";
	  std::cout<<"chi**2="<<unfold.GetChi2A()<<"+"<<unfold.GetChi2L()
	           <<" / "<<unfold.GetNdf()<<"\n";
	  std::cout<<"chi**2(sys)="<<unfold.GetChi2Sys()<<"\n";
	  //==========================================================================
	  // create graphs with one point to visualize the best choice of tau
	  //
	  Double_t t[1],x[1],y[1];
	  logTauX->GetKnot(iBest,t[0],x[0]);
	  logTauY->GetKnot(iBest,t[0],y[0]);
	  TGraph *bestLcurve=new TGraph(1,x,y);
	  TGraph *bestLogTauLogChi2=new TGraph(1,t,x);
	  //==========================================================================  
	  // retreive results into histograms
	  // get unfolded distribution
	  TH1 *histMunfold=unfold.GetOutput("Unfolded");

	  // get error matrix (input distribution [stat] errors only)
	  // TH2D *histEmatData=unfold.GetEmatrix("EmatData");
	  // get total error matrix:
	  //   migration matrix uncorrelated and correlated systematic errors
	  //   added in quadrature to the data statistical errors
	  TH2 *histEmatTotal=unfold.GetEmatrixTotal("EmatTotal");
	  // create data histogram with the total errors
	  /*TH1D *histTotalError=
	     new TH1D("TotalError",";mass(gen)",nGen,xminGen,xmaxGen);
	  for(Int_t bin=1;bin<=sizeBins;bin++) {
	    histTotalError->SetBinContent(bin,histMunfold->GetBinContent(bin));
	    histTotalError->SetBinError
	       (bin,TMath::Sqrt(histEmatTotal->GetBinContent(bin,bin)));
	  }*/

	  // get global correlation coefficients
	  // for this calculation one has to specify whether the
	  // underflow/overflow bins are included or not
	  // default: include all binsv
	  // here: exclude underflow and overflow bins
	  TH2 *gHistInvEMatrix;
	  TH1 *histRhoi=unfold.GetRhoItotal("rho_I",
	                                    0, // use default title
	                                    0, // all distributions
	                                    "*[UO]", // discard underflow and overflow bins on all axes
	                                    kTRUE, // use original binning
	                                    &gHistInvEMatrix // store inverse of error matrix
	                                    );
	  //TCanvas *can_test = new TCanvas("can_test", "can_test", 800,600);
	  //histRhoi->Draw();

	  return histMunfold;
}
  


