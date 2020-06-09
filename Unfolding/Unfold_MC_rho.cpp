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
#include "TUnfoldDensity.h"

using namespace std;

using std::cin;
using std::cout;
using std::endl;

#include "TemplateConstantsUnfold.h"

TH1F *getRebinned(TH1F *h, float BND[], int N);
TH1 *unfoldedOutput(TH2F *hResponse_, TH1F *hReco, float BND[], int sizeBins, TString variable);
TH1 *unfoldedOutputRho(TH2F *hResponse_, TH1F *hReco, float BND[], int sizeBins, TString variable);
TH1 *unfoldedOutput_LCurve(TH2F *hResponse_, TH1F *hReco, float BND[], int sizeBins, TString variable);
TString year;
TString varParton;


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


void Unfold_MC_rho(TString inYear = "2016", bool isParton = true, int unfoldMethod=3)
{
  year = inYear;
  initFilesMapping();
  gStyle->SetOptStat(0);

  std::vector< std::vector <Float_t> > const BND_reco = {{1000, 1100,1200,1300, 1400,1500, 1600,1700, 1800,1900, 2000,2200, 2400,2600,2800,3000,3500,4200,5000}, //mjj 19
                                                        {0,30,60,105,150,225,300,375,450,525,600,675,750,850,950,1100,1300}, //ptjj 17
                                                        {-2.4,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.4}, //yjj
                                                        {400,425,450,475,500,535,570,610,650,700,750,800,850,900,950,1025,1100,1200,1300,1400,1500}, //jetPt0 21
                                                        {400,425,450,475,500,535,570,610,650,700,750,800,850,900,950,1100,1300,1500}, //jetPt1 18
                                                        {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4}, //jetY0 25
                                                        {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4}}; //jetY1 25

   std::vector< std::vector <Float_t> > const BND_gen = {{1000, 1200, 1400, 1600, 1800, 2000, 2400, 2800, 3200, 4000, 5000}, //mjj
                                                         {0,60,150,300,450,600,750,1000,1300}, //ptjj
                                                         {-2.4,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.4}, //yjj
                                                         {400,450,500,570,650,750,850,950,1100,1300,1500}, //jetPt0
                                                         {400,450,500,570,650,750,850,1000,1200,1500}, //jetPt1
                                                         {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}, //jetY0
                                                         {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4}}; //jetY1

  float LUMI = luminosity[year];
  //get the files:
  //1. the signal file has the fiducial measurements that are going to be used as input
  TFile *signalFile;

  //TFile *signalFile = TFile::Open(TString::Format("../MassFit/Mixed/%s/Histo_TT_Mtt-700toInf_100_reduced_UnequalBinning.root",
  //								  year.Data()));
  //i also have the files in the same files with the data:
  //TFile *signalFile = TFile::Open(TString::Format("../MassFit/%s/FiducialMeasurement/UnequalBinning/SignalHistograms.root",
                    //year.Data()));
  //2. This file has the response matrices as well as the efficiency and acceptance for the signal procedure
  TFile *effAccInf = TFile::Open(TString::Format("../ResponseMatrices/%s/UnequalBins/ResponsesEfficiencyNominalMC_%s.root", year.Data(), year.Data()));
  //TFile *effAccInf = TFile::Open(TString::Format("../ResponseMatrices/%s/UnequalBins/ResponsesEfficiency_%s.root", year.Data(), year.Data()));
  //3. This file is the theoretical parton/particle file that we use for comparison
  TFile *infTheory = TFile::Open(TString::Format("%s/TheoryTemplatesNominalMC.root", year.Data()));


  //whether parton or particle, from the choice of the user
  varParton = "Parton";
  if(!isParton) varParton = "Particle";

  //get the number of bins for each
  int NBINS[BND_reco.size()];
  int NBINS_GEN[BND_gen.size()];
  const int NVAR = 7;
  for (int i = 0; i<BND_reco.size(); i++)
  	NBINS[i] = BND_reco[i].size()-1;
  for (int i = 0; i<BND_gen.size(); i++)
  	NBINS_GEN[i] = BND_gen[i].size()-1;

  TH1F *inSig[BND_reco.size()], *hSig[BND_reco.size()], *hSig_Init[BND_reco.size()];
  TString variable[NVAR] = {"mJJ", "ptJJ", "yJJ", "jetPt0", "jetPt1","jetY0", "jetY1"};
  TString variableGen[NVAR] = {"mJJGen", "ptJJGen", "yJJGen", "genjetPt0", "genjetPt1","genjetY0", "genjetY1"};
  TString variableParton[NVAR] = {"mTTbarParton", "ptTTbarParton", "yTTbarParton", "partonPt0", "partonPt1","partonY0", "partonY1"};

  TH2F *hResponse[BND_reco.size()];
  TUnfold *unf[BND_reco.size()];
  TH1 *hUnf[BND_reco.size()];

  //theoretical Histogram
  TH1F *hTheory[BND_reco.size()], *hTheory_2[BND_reco.size()];
  TH1F *hErrorBefore[BND_reco.size()], *hErrorAfter[BND_reco.size()];


  TCanvas *can[BND_reco.size()], *can_rho[BND_reco.size()], *canError[BND_reco.size()];
  TLegend *leg[BND_reco.size()];
  TH1F *hUnfTemp[BND_reco.size()], *hTheoryTemp[BND_reco.size()];
  TH1F *hUnfFinal[BND_reco.size()], *hTheoryFinal[BND_reco.size()];


  TString unfMethodStr = "";
  if(unfoldMethod ==2) unfMethodStr = "_LCurveMethod";
  else if(unfoldMethod ==3) unfMethodStr = "_RhoMethod";
  //TFile *gfile = TFile::Open("2016/output_2016_mcSig_nom_reduced.root");
  TFile *outf = TFile::Open(TString::Format("%s/%sMeasurements/MC/OutputFile%s.root", year.Data(), varParton.Data(), unfMethodStr.Data()),"RECREATE");


  for(int ivar = 0; ivar<BND_reco.size(); ivar++)
  {

    signalFile = TFile::Open(TString::Format("../MassFit/%s/Histo_TT_NominalMC_100_reduced_UnequalBinning.root",
                     year.Data(), variable[ivar].Data()));

    int sizeBins = NBINS[ivar];
    float tempBND[NBINS[ivar]+1];
    std::copy(BND_reco[ivar].begin(), BND_reco[ivar].end(), tempBND);

    float tempBNDGen[NBINS_GEN[ivar]+1];
    std::copy(BND_gen[ivar].begin(), BND_gen[ivar].end(), tempBNDGen);
    //from signal file get the initial S_j with j bins ~ 2* parton bins (i)
    hSig_Init[ivar] = (TH1F*)signalFile->Get(TString::Format("hWt_%s_2btag_expYield",variable[ivar].Data()));
    //could also use this from data file because I also store the MC output:
    //hSig_Init[ivar] = (TH1F*)signalFile->Get(TString::Format("hSMC_%s",variable[ivar].Data()));

    for(int j=1; j<hSig_Init[ivar]->GetNbinsX(); j++)
    {
    	cout<<"hSig_Init "<<j<<": "<<hSig_Init[ivar]->GetBinContent(j)<< " hSig_Init Error: "<<hSig_Init[ivar]->GetBinError(j)<<endl;
    }
    //hSig[ivar] = getRebinned(hSig_Init[ivar], tempBND, NBINS[ivar]);
    hSig[ivar] = (TH1F*)hSig_Init[ivar]->Clone(TString::Format("hSig_%s", variable[ivar].Data()));

    //hSig[ivar] = (TH1F*)gfile->Get("2btag_mJJ_nominal");
    //hSig[ivar]->Scale(832 * 35920);
    leg[ivar] = new TLegend(0.65,0.7,0.9,0.9);

    //set the new content and get acceptance
    TEfficiency *acceptance =  (TEfficiency*)effAccInf->Get(TString::Format("Acceptance%s_%s",varParton.Data(), variable[ivar].Data()));

    for(int j =1; j<hSig[ivar]->GetNbinsX()+1; j++)
    {
      float acc = acceptance->GetEfficiency(j);
      hSig[ivar]->SetBinContent(j, acc*hSig[ivar]->GetBinContent(j));
      //handle errors as well--> asymmetric error from acceptance
      //cout<<"hSig "<<j<<": "<<hSig[ivar]->GetBinContent(j)<< "hSig Error: "<<hSig[ivar]->GetBinError(j)<<endl;
      //cout<<"hSig Error "<<j<<": "<<hSig[ivar]->GetBinError(j)<<endl;

      float accError = (acceptance->GetEfficiencyErrorLow(j) + acceptance->GetEfficiencyErrorUp(j))/2;
      cout<<"acceptance: "<<acc<<endl;
      if(accError > 0)
      	  //error propagation
      	hSig[ivar]->SetBinError(j,TMath::Sqrt(TMath::Power(accError*hSig[ivar]->GetBinContent(j),2) + TMath::Power(hSig[ivar]->GetBinError(j)*acc,2)));
	  else
	  	hSig[ivar]->SetBinError(j, hSig[ivar]->GetBinError(j));
    }
    //break;

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
    if(unfoldMethod ==1)
    	hUnf[ivar] = unfoldedOutput(hResponse[ivar], hSig[ivar], tempBNDGen, NBINS_GEN[ivar], variable[ivar]);
    else if(unfoldMethod ==2)
    	hUnf[ivar] = unfoldedOutput_LCurve(hResponse[ivar], hSig[ivar], tempBNDGen, NBINS_GEN[ivar], variable[ivar]);
    else
    	hUnf[ivar] = unfoldedOutputRho(hResponse[ivar], hSig[ivar], tempBNDGen, NBINS_GEN[ivar], variable[ivar]);
    //continue;
    if(variable[ivar].EqualTo("yJJ"))
      hUnf[ivar]->GetXaxis()->SetTitle(variable[ivar]);
    else if(variable[ivar].Contains("jetY"))
    	hUnf[ivar]->GetXaxis()->SetTitle("|"+variable[ivar]+"|");
    else
    	hUnf[ivar]->GetXaxis()->SetTitle(TString::Format("%s (GeV)", variable[ivar].Data()));

    hUnf[ivar]->GetYaxis()->SetTitle(TString::Format("#frac{d#sigma}{d#chi} %s", varParton.Data()));
    hUnf[ivar]->GetYaxis()->SetTitleOffset(1.4);

     //here get the errors:


    hErrorBefore[ivar] = (TH1F*)hSig[ivar]->Clone(TString::Format("hErrorBefore_%s", variable[ivar].Data()));
    hErrorAfter[ivar] = (TH1F*)hUnf[ivar]->Clone(TString::Format("hErrorAfter_%s", variable[ivar].Data()));
    hErrorBefore[ivar]->Clear();
    hErrorAfter[ivar]->Clear();

    for (int i =1; i<hUnf[ivar]->GetNbinsX()+1; i++)
    {
      hErrorAfter[ivar]->SetBinContent(i,hUnf[ivar]->GetBinError(i));
      cout<<"hUnf after: "<<hUnf[ivar]->GetBinContent(i)<<endl;
      //cout<<"Error after for bin "<<i<<": "<<hErrorAfter[ivar]->GetBinError(i)<<endl;
    }
    for(int j=1; j<hSig[ivar]->GetNbinsX()+1; j++)
    {
      hErrorBefore[ivar]->SetBinContent(j,hSig[ivar]->GetBinError(j));
      //cout<<"Error before for bin "<<j<<": "<<hErrorBefore[ivar]->GetBinError(j)<<endl;
    }

    //break;
    TEfficiency *efficiency =  (TEfficiency*)effAccInf->Get(TString::Format("Efficiency%s_%s",varParton.Data(), tempVar.Data()));

    if(!variable[ivar].EqualTo("yJJ"))hErrorBefore[ivar]->Rebin(2);
    hErrorBefore[ivar]->SetLineColor(kBlue);
    hErrorBefore[ivar]->SetMarkerColor(kBlue);
    hErrorBefore[ivar]->SetMarkerStyle(20);
    hErrorAfter[ivar]->SetLineColor(kRed);
    hErrorAfter[ivar]->SetMarkerColor(kRed);
    hErrorAfter[ivar]->SetMarkerStyle(23);

    canError[ivar] = new TCanvas(TString::Format("canError_%s",variable[ivar].Data()),TString::Format("canError_%s",variable[ivar].Data()) , 800,600);
    canError[ivar]->cd();
    hErrorAfter[ivar]->Draw("hist");
    hErrorBefore[ivar]->Draw("hist same");


    for(int i =1; i<hUnf[ivar]->GetNbinsX()+1; i++)
    {
      float eff = efficiency->GetEfficiency(i);
      if(eff >0)
      {
	      cout<<"i= "<<i<<" eff="<<eff<<endl;
        float oldContent = hUnf[ivar]->GetBinContent(i);
	      float newContent = hUnf[ivar]->GetBinContent(i)/eff;
	      hUnf[ivar]->SetBinContent(i, newContent);
          //cout<<"old: "<<oldContent<<endl;
	      //cout<<"new: "<<hUnf[ivar]->GetBinContent(i)<<endl;
	      //handle errors as well--> asymmetric error from efficiency
	      float effError = (efficiency->GetEfficiencyErrorLow(i) + efficiency->GetEfficiencyErrorUp(i))/2;
	      hUnf[ivar]->SetBinError(i,TMath::Sqrt(TMath::Power(effError*hUnf[ivar]->GetBinContent(i)/TMath::Power(eff,2),2) + TMath::Power(hUnf[ivar]->GetBinError(i)/eff,2)));
  	  }
    }

    can[ivar] = new TCanvas(TString::Format("can_%s",variable[ivar].Data()),TString::Format("can_%s",variable[ivar].Data()) , 800,600);
    can[ivar]->cd();
    auto *closure_padRatio = new TPad("closure_pad2","closure_pad2",0.,0.,1.,0.3);
    closure_padRatio->Draw();
    closure_padRatio->SetTopMargin(0.05);
    closure_padRatio->SetBottomMargin(0.3);
    closure_padRatio->SetGrid();

    auto *closure_pad1 = new TPad("closure_pad1","closure_pad1",0.,0.3,1.,1.);
    closure_pad1->Draw();
    closure_pad1->SetBottomMargin(0.005);
    closure_pad1->cd();

    //theory histogram:
    hTheory_2[ivar] = (TH1F*)infTheory->Get(TString::Format("h%s_%s", varParton.Data(), variable[ivar].Data()));
    hTheory[ivar] = (TH1F*)efficiency->GetCopyTotalHisto();


    //hTheory_2[ivar]->SetLineColor(kGreen+2);
    hUnfFinal[ivar]= (TH1F*)hUnf[ivar]->Clone(TString::Format("hUnfFinal_%s", variable[ivar].Data()));
    hTheoryFinal[ivar]= (TH1F*)hTheory[ivar]->Clone(TString::Format("hTheoryFinal_%s", variable[ivar].Data()));
    hTheory[ivar]->Scale(1/luminosity[year], "width");
    hUnf[ivar]->Scale(1/luminosity[year], "width");

    hUnf[ivar]->SetLineColor(kBlue);
    hTheory[ivar]->SetLineColor(kRed);
    hUnf[ivar]->SetMarkerStyle(20);
    hUnf[ivar]->SetMarkerColor(kBlue);
    hUnf[ivar]->SetTitle(TString::Format("%s Unfolded vs Theory %s %s",varParton.Data(),variable[ivar].Data(),year.Data()));
    hTheory[ivar]->SetTitle(TString::Format("%s Unfolded vs Theory %s %s",varParton.Data(),variable[ivar].Data(),year.Data()));

    hUnf[ivar]->GetYaxis()->SetTitle("#frac{d#sigma}{d#chi}");
    hTheory[ivar]->GetYaxis()->SetTitle("#frac{d#sigma}{d#chi}");
    hTheory[ivar]->SetMarkerStyle(23);
    hTheory[ivar]->SetMarkerColor(kRed);

    leg[ivar]->AddEntry(hUnf[ivar], "Unfolded", "lpe");
    leg[ivar]->AddEntry(hTheory[ivar], "Theory", "lpe");


  	hTheory[ivar]->Draw();
  	hUnf[ivar]->Draw("same");
  	leg[ivar]->Draw();

  	if(!variable[ivar].EqualTo("jetY0") && !variable[ivar].EqualTo("jetY1")) gPad->SetLogy();

  	closure_padRatio->cd();
  	hUnfTemp[ivar] = (TH1F*)hUnf[ivar]->Clone(TString::Format("hUnf_%s", variable[ivar].Data()));
  	hTheoryTemp[ivar] = (TH1F*)hTheory[ivar]->Clone(TString::Format("hTheory_%s", variable[ivar].Data()));

  	hUnfTemp[ivar]->Divide(hTheoryTemp[ivar]);
  	hUnfTemp[ivar]->Draw();

  	hUnfTemp[ivar]->SetTitle("");
    hUnfTemp[ivar]->GetYaxis()->SetTitle("#frac{Unfolded}{Theory}");
    hUnfTemp[ivar]->GetYaxis()->SetTitleSize(14);
    hUnfTemp[ivar]->GetYaxis()->SetTitleFont(43);
    hUnfTemp[ivar]->GetYaxis()->SetTitleOffset(1.55);
    hUnfTemp[ivar]->GetYaxis()->SetLabelFont(43);
    hUnfTemp[ivar]->GetYaxis()->SetLabelSize(15);
    hUnfTemp[ivar]->GetXaxis()->SetTitleSize(0.09);
    hUnfTemp[ivar]->GetYaxis()->SetRangeUser(0,3);

    hUnfTemp[ivar]->SetLineColor(kRed);
    hUnfTemp[ivar]->SetMarkerStyle(20);
    hUnfTemp[ivar]->SetMarkerColor(kRed);
    hUnfTemp[ivar]->Draw();
    hUnfTemp[ivar]->GetXaxis()->SetLabelSize(0.09);


    //hTheory_2[ivar]->Draw("same");
   // cout<<"------"<<endl;
   // cout<<"theory: "<<hTheory_2[ivar]->Integral()<<endl;
    //cout<<"theory from eff: "<<hTheory[ivar]->Integral()<<endl;

    outf->cd();
    hTheory[ivar]->Write(TString::Format("hTheory_%s", variable[ivar].Data()));
    hTheoryFinal[ivar]->Write(TString::Format("hTheoryFinal%s", variable[ivar].Data()));
  	hUnf[ivar]->Write(TString::Format("hUnfold_%s", variable[ivar].Data()));
    hUnfFinal[ivar]->Write(TString::Format("hUnfoldFinal_%s", variable[ivar].Data()));
  	hErrorAfter[ivar]->Write(TString::Format("hErrorAfter_%s", variable[ivar].Data()));
    hErrorBefore[ivar]->Write(TString::Format("hErrorBefore_%s", variable[ivar].Data()));
    can[ivar]->Print(TString::Format("%s/%sMeasurements/MC/Unfold_%s%s.pdf",year.Data(),varParton.Data(),variable[ivar].Data(), unfMethodStr.Data()), "pdf");
 //   break;
  }

}

TH1 *unfoldedOutput(TH2F *hResponse_, TH1F *hReco, float BND[], int sizeBins, TString variable)
{
  cout<<"Using Standard Unfolding method method"<<endl;
	//TCanvas *can_response = new TCanvas("can_response", "can_response", 800,600);
	//hResponse_->Draw("BOX");
  cout<<variable<<endl;
  TUnfold unfold(hResponse_,TUnfold::kHistMapOutputHoriz, TUnfold::kRegModeNone, TUnfold::kEConstraintArea);


  /*
    Return value: nError1+10000*nError2:

    -->nError1: number of bins where the uncertainty is zero.
                these bins either are not used for the unfolding (if oneOverZeroError==0) or 1/uncertainty is set to oneOverZeroError.
    -->nError2: return values>10000 are fatal errors, because the unfolding can not be done.
                The number nError2 corresponds to the number of truth bins which are not constrained by data points
  */
  if(unfold.SetInput(hReco)>=10000) {
    std::cout<<"Unfolding result may be wrong\n";
  }
      //========================================================================
	  // the unfolding is done here
  //========================================================================
  // the unfolding is done here

  unfold.DoUnfold(10E-9);
  //TH1 *histMunfold=unfold.GetOutput("Unfolded");

  //set up a bin map, excluding underflow and overflow bins
  // the binMap relates the the output of the unfolding to the final
  // histogram bins
  Int_t *binMap=new Int_t[sizeBins+2];
  for(Int_t i=1;i<=sizeBins;i++) binMap[i]=i;
  binMap[0]=-1;
  binMap[sizeBins+1]=-1;

  TH1F *histMunfold = new TH1F (TString::Format("UnfoldedOutput_%s",variable.Data()),TString::Format("UnfoldedOutput_%s",variable.Data()),
                               sizeBins, BND);
  unfold.GetOutput(histMunfold, binMap);
  return histMunfold;

}

TH1 *unfoldedOutputRho(TH2F *hResponse_, TH1F *hReco, float BND[], int sizeBins, TString variable)
{
  cout<<"Using Rho method"<<endl;

  TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
  TUnfold::EConstraint constraintMode = TUnfold::kEConstraintArea;
  TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;

  const TArrayD *arrParton = hResponse_->GetXaxis()->GetXbins();
  Double_t partonBins[arrParton->GetSize()];
  for (int i = 0; i < arrParton->GetSize(); i++)
  {
    partonBins[i] = arrParton->At(i);
  }

  TUnfoldBinning *partonBinning = new TUnfoldBinning("partonBinning");
  TUnfoldBinning *partonDistribution = partonBinning->AddBinning("partonDistribution");
  partonDistribution->AddAxis("partonPt", arrParton->GetSize() - 1, partonBins, false, false);
  //partonBinning->AddAxis("partonPt", arrParton->GetSize() - 1, partonBins, false, false);
  //std::cout << "partonBinning: " << (partonBinning->GetDistributionNumberOfBins() > 0) << std::endl;
  //std::cout << "partonDistribution: " << partonDistribution->GetDistributionNumberOfBins() << std::endl;

  const TArrayD *arrReco = hResponse_->GetYaxis()->GetXbins();
  Double_t recoBins[arrReco->GetSize()];
  for (int i = 0; i < arrReco->GetSize(); i++)
  {
    recoBins[i] = arrReco->At(i);
  }
  TUnfoldBinning *recoBinning = new TUnfoldBinning("recoBinning");
  TUnfoldBinning *recoDistribution = recoBinning->AddBinning("recoDistribution");
  recoDistribution->AddAxis("recoPt", arrReco->GetSize() - 1, recoBins, false, false);
  //recoBinning->AddAxis("recoPt", arrReco->GetSize() - 1, recoBins, false, false);
  //std::cout << "partonBinning: " << (recoBinning->GetDistributionNumberOfBins() > 0) << std::endl;
  //std::cout << "partonDistribution: " << recoDistribution->GetDistributionNumberOfBins() << std::endl;

  std::cout << "Bins: " << arrReco->GetSize() << " " << arrParton->GetSize() << std::endl;
  const char *REGULARISATION_DISTRIBUTION = 0;
  const char *REGULARISATION_AXISSTEERING = "*[B]";


  TUnfoldDensity unfold(hResponse_, TUnfold::kHistMapOutputHoriz,
                        regMode, constraintMode, densityFlags,
                        partonBinning, recoBinning,
                        REGULARISATION_DISTRIBUTION,
                        REGULARISATION_AXISSTEERING);


   unfold.SetInput(hReco);

  Int_t nScan = 30;
  TSpline *rhoLogTau = 0;
  TGraph *lCurve = 0;

  // for determining tau, scan the correlation coefficients
  // correlation coefficients may be probed for all distributions
  // or only for selected distributions
  // underflow/overflow bins may be included/excluded
  //
  const char *SCAN_DISTRIBUTION = "partonBinning";
  const char *SCAN_AXISSTEERING = 0;

  Int_t iBest = unfold.ScanTau(nScan, 0., 2., &rhoLogTau,
                               TUnfoldDensity::kEScanTauRhoMax,
                               //TUnfoldDensity::kEScanTauRhoAvg,
                               SCAN_DISTRIBUTION, SCAN_AXISSTEERING,
                               &lCurve);

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

  std::cout <<variable <<" Tau: " << unfold.GetTau() << std::endl;

  TCanvas *c1 = new TCanvas(TString::Format("Canvas_%s",variable.Data()), TString::Format("Canvas_%s",variable.Data()), 600, 500);

  rhoLogTau->Draw();
  knots->Draw("*");
  bestRhoLogTau->SetMarkerColor(kRed);
  bestRhoLogTau->Draw("*");

  c1->Print(TString::Format("%s/%sMeasurements/MC/GlobalCorrelationGraph_%s.pdf",year.Data(), varParton.Data(), variable.Data()),"pdf");

  TH1F *histMunfold = (TH1F*)unfold.GetOutput(TString::Format("UnfoldedOutput_%s",variable.Data())); //, 0, 0, 0, kFALSE);
  return histMunfold;


}


TH1 *unfoldedOutput_LCurve(TH2F *hResponse_, TH1F *hReco, float BND[], int sizeBins, TString variable)
{
  cout<<"Using Standard Unfolding method method"<<endl;
  //TCanvas *can_response = new TCanvas("can_response", "can_response", 800,600);
  cout<<variable<<endl;
  TUnfold unfold(hResponse_,TUnfold::kHistMapOutputHoriz, TUnfold::kRegModeSize, TUnfold::kEConstraintArea);

  if(unfold.SetInput(hReco)>=10000) {
    std::cout<<"Unfolding result may be wrong\n";
  }
  //========================================================================
  // the unfolding is done here
  //========================================================================
  // the unfolding is done here
  //
  // scan L curve and find best point
  Int_t nScan=30;
  // use automatic L-curve scan: start with taumin=taumax=0.0
  Double_t tauMin=0.0;
  Double_t tauMax=0.0;
  Int_t iBest;
  TSpline *logTauX,*logTauY;
  TGraph *lCurve;

  // this method scans the parameter tau and finds the kink in the L curve
  // finally, the unfolding is done for the best choice of tau
  iBest=unfold.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
  std::cout<<"tau="<<unfold.GetTau()<<endl;
  //unfold.DoUnfold(unfold.GetTau());
  //TH1 *histMunfold=unfold.GetOutput("Unfolded");

  //set up a bin map, excluding underflow and overflow bins
  // the binMap relates the the output of the unfolding to the final
  // histogram bins
  Int_t *binMap=new Int_t[sizeBins+2];
  for(Int_t i=1;i<=sizeBins;i++) binMap[i]=i;
  binMap[0]=-1;
  binMap[sizeBins+1]=-1;

  TH1F *histMunfold = new TH1F (TString::Format("UnfoldedOutput_%s",variable.Data()),TString::Format("UnfoldedOutput_%s",variable.Data()),
                               sizeBins, BND);
  unfold.GetOutput(histMunfold, binMap);
  return histMunfold;

}
