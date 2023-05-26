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
#include "../../CMS_plots/CMS_lumi.C"
#include "../../CMS_plots/tdrstyle.C"

using namespace std;

#include "UseTUnfoldDensity.cpp"

using std::cin;
using std::cout;
using std::endl;

#include "TemplateConstants.h"


void PlotHistogramErrors(TH1F* hist1, TH1F *hist2, TH1F *hist3);

void Unfold_Combined(TString dir, TString inputFile, TString isParton = "true", int unfoldmethod=1)
{
  bool isNorm = false;
  initFilesMapping();
  gStyle->SetOptStat(0);
  TString year = "2018";
  TString tempFileName;
  if(dir.EqualTo("PDFWeights")) tempFileName = "_pdf_"+inputFile;
  else if (dir.EqualTo("ScaleWeights")) tempFileName = "_scale_"+inputFile;
  else if (dir.EqualTo("Nominal")) tempFileName = "";
  else if (dir.EqualTo("NominalTopSF")) tempFileName = "";
  else tempFileName = "_"+inputFile;
  std::vector< std::vector <Float_t> > const BND_gen = {//{1000, 1200, 1400, 1600, 1800, 2000, 2400, 3000, 5000}, //mJJ
                                                        //{0, 60, 150, 300, 450, 850, 1300}, //ptjj
                                                        //{-2.4, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.4}, //yjj
                                                        // {450, 500, 570, 650, 750, 850, 950, 1100, 1300, 1500, 2000}}; //jetPt0
                                                        // {400, 450, 500, 570, 650, 800, 1100, 1600}}; //jetPt1
                                                        {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.4}}; //jetY0
                                                        //{0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.4}}; //jetY1
                                                        // {1,2,3,4,5,6,7,8,9,10,13,16}, //chi
                                                        // {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1}, //|cosTheta*| leading
                                                        // {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1}}; //|cosTheta*| subleading

  // std::vector< std::vector <Float_t> > const BND_reco ={{1000, 1200, 1400, 1600, 1800, 2000, 2400, 3000, 5000}, //mJJ
  //                                                       {0, 60, 150, 300, 450, 850, 1300}, //ptJJ
  //                                                       {-2.4, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.4}, //yjj
  //                                                       {450, 500, 570, 650, 800, 1100, 1500}, //jetPt0
  //                                                       {400, 450, 500, 570, 650, 800, 1100, 1500}, //jetPt1
  //                                                       {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.4}, //jetY0
  //                                                       {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.4}, //jetY1
  //                                                       {1,2,3,4,5,6,7,8,9,10,13,16}, //chi
  //                                                       {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1}, //|cosTheta*| leading
  //                                                       {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1}}; //|cosTheta*| subleading
   // will be used for reco
  std::vector< std::vector <Float_t> > const BND_reco = { //{1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2700, 3000, 4000, 5000}, //mJJ
                       //{0, 30, 60, 100, 150, 225, 300, 375, 450, 600, 850, 1000, 1300}, //ptJJ
                       //{-2.4, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.4}, //yjj
                      //  {450, 475, 500, 535, 570, 610, 650, 700, 750, 800, 850, 900, 950, 1000, 1100, 1200, 1300, 1500, 2000}}; //jetPt0
                      //  {400, 425, 450, 475, 500, 535, 570, 610, 650, 700, 800, 1000, 1100, 1300, 1600}}; //jetPt1
                       {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0, 2.4}}; //jetY0
                       //{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0, 2.4}}; //jetY1
                      //  {1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 11, 13, 14, 16}, //chi
                      //  {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1}, //|cosTheta*| leading
                      //  {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1}}; //|cosTheta*| subleading


  //get the files:
  //1. the signal file has the fiducial measurements that are going to be used as input
  TFile *signalFile;

  //whether parton or particle, from the choice of the user
  TString varParton = "Parton";
  if(isParton.EqualTo("false")) varParton = "Particle";

  //get the number of bins for each
  int NBINS[BND_reco.size()];
  int NBINS_GEN[BND_gen.size()];
  // const int NVAR = 7;
  const int NVAR = 1;
  for (int i = 0; i<BND_reco.size(); i++)
  	NBINS[i] = BND_reco[i].size()-1;
  for (int i = 0; i<BND_gen.size(); i++)
  	NBINS_GEN[i] = BND_gen[i].size()-1;

  TH1F *hSig[BND_reco.size()], *hSig_orthogonal[BND_reco.size()];
  // TString variable[NVAR]   = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1"};
  // TString variableParton[NVAR] = {"mTTbarParton", "ptTTbarParton", "yTTbarParton","partonPt0", "partonPt1", "partonY0", "partonY1"};
  // TString variableGen[NVAR] = {"mJJGen", "ptJJGen", "yJJGen","genjetPt0", "genjetPt1", "genjetY0", "genjetY1"};
  // TString variable[NVAR]   = {"jetPt0"};
  // TString variableParton[NVAR] = {"partonPt0"};
  // TString variableGen[NVAR] = {"genjetPt0"};
  TString variable[NVAR]   = {"jetY0"};
  TString variableParton[NVAR] = {"partonY0"};
  TString variableGen[NVAR] = {"genjetY0"};
  //2. This file has the response matrices as well as the efficiency and acceptance for the signal procedure
  //handle responses efficiency, acceptance
  TH2F *hResponse[BND_reco.size()], *hResponse_orthogonal[BND_reco.size()];
  // TH1F *efficiency[BND_reco.size()], *acceptance[BND_reco.size()];
  TH1F *efficiency_denom[BND_reco.size()], *efficiency_denom_orthogonal[BND_reco.size()];
  TH1F *fiducial_theory[BND_reco.size()];
  TH1F *hFidTheory[BND_reco.size()], *hFidTheory_norm[BND_reco.size()];
  TH1F *hTheory[BND_reco.size()], *hTheory_norm[BND_reco.size()];
  TH1F *hTheory_orthogonal[BND_reco.size()], *hTheory_norm_orthogonal[BND_reco.size()];

  TH2F *hResponse_had[BND_reco.size()], *hResponse_sem[BND_reco.size()], *hResponse_dil[BND_reco.size()];
  TH2F *hResponse_had_orthogonal[BND_reco.size()], *hResponse_sem_orthogonal[BND_reco.size()], *hResponse_dil_orthogonal[BND_reco.size()];
  TEfficiency *efficiency[BND_reco.size()], *acceptance[BND_reco.size()];
  TEfficiency *efficiency_orthogonal[BND_reco.size()], *acceptance_orthogonal[BND_reco.size()];

  // get responses and add files 
  cout<<TString::Format("%s/Responses%s/ResponsesEfficiency_TTToHadronic%s.root", year.Data(), dir.Data(), tempFileName.Data())<<endl;
  TFile *inf_had = TFile::Open(TString::Format("%s/Responses%s/ResponsesEfficiency_TTToHadronic%s.root", year.Data(), dir.Data(), tempFileName.Data()));
  TFile *inf_sem = TFile::Open(TString::Format("%s/Responses%s/ResponsesEfficiency_TTToSemiLeptonic%s.root", year.Data(), dir.Data(), tempFileName.Data()));
  TFile *inf_dil = TFile::Open(TString::Format("%s/Responses%s/ResponsesEfficiency_TTTo2L2Nu%s.root", year.Data(), dir.Data(), tempFileName.Data()));

  TFile *inf_had_orthogonal = TFile::Open(TString::Format("%s/Responses%s/ResponsesEfficiency_TTToHadronic%s_orthogonal.root", year.Data(), dir.Data(), tempFileName.Data()));
  TFile *inf_sem_orthogonal = TFile::Open(TString::Format("%s/Responses%s/ResponsesEfficiency_TTToSemiLeptonic%s_orthogonal.root", year.Data(), dir.Data(), tempFileName.Data()));
  TFile *inf_dil_orthogonal = TFile::Open(TString::Format("%s/Responses%s/ResponsesEfficiency_TTTo2L2Nu%s_orthogonal.root", year.Data(), dir.Data(), tempFileName.Data()));

  float LUMI = luminosity[TString::Format("luminosity%s", year.Data())];
  //get response matrix
  for(int ivar = 0; ivar<BND_reco.size(); ivar++)
  {
      TString tempVar;
      if(isParton)
          tempVar = variableParton[ivar];
      else
          tempVar = variableGen[ivar];

      cout<< TString::Format("h%sResponse_%s",varParton.Data(), variable[ivar].Data())<<endl;

      hResponse_had[ivar] = (TH2F*)inf_had->Get(TString::Format("h%sResponse_%s",varParton.Data(), variable[ivar].Data()));
      hResponse_sem[ivar] = (TH2F*)inf_sem->Get(TString::Format("h%sResponse_%s",varParton.Data(), variable[ivar].Data()));
      hResponse_dil[ivar] = (TH2F*)inf_dil->Get(TString::Format("h%sResponse_%s",varParton.Data(), variable[ivar].Data()));

      hResponse_had_orthogonal[ivar] = (TH2F*)inf_had_orthogonal->Get(TString::Format("h%sResponse_%s",varParton.Data(), variable[ivar].Data()));
      hResponse_sem_orthogonal[ivar] = (TH2F*)inf_sem_orthogonal->Get(TString::Format("h%sResponse_%s",varParton.Data(), variable[ivar].Data()));
      hResponse_dil_orthogonal[ivar] = (TH2F*)inf_dil_orthogonal->Get(TString::Format("h%sResponse_%s",varParton.Data(), variable[ivar].Data()));

      // add all (regardless of ttbar process) to the hResponse per variable
      // responses are already scaled to XSEC and lumi
      hResponse[ivar] = (TH2F*)hResponse_had[ivar]->Clone();
      hResponse[ivar]->Add(hResponse_sem[ivar]);
      hResponse[ivar]->Add(hResponse_dil[ivar]);

      hResponse_orthogonal[ivar] = (TH2F*)hResponse_had_orthogonal[ivar]->Clone();
      hResponse_orthogonal[ivar]->Add(hResponse_sem_orthogonal[ivar]);
      hResponse_orthogonal[ivar]->Add(hResponse_dil_orthogonal[ivar]);

      /*
      cout<<"----"<<variable[ivar]<<"----"<<endl;
      cout<<"Integral responses: "<<hResponse[ivar]->Integral()<<endl;
      cout<<"Entries responses: "<<hResponse[ivar]->GetEntries()<<endl;
      */

    // read efficiency and acceptance from responses files directories
    TEfficiency *eff_had = (TEfficiency*)inf_had->Get(TString::Format("Efficiency%s_%s",varParton.Data(), tempVar.Data()));
    TEfficiency *eff_sem = (TEfficiency*)inf_sem->Get(TString::Format("Efficiency%s_%s",varParton.Data(), tempVar.Data()));
    TEfficiency *eff_dil = (TEfficiency*)inf_dil->Get(TString::Format("Efficiency%s_%s",varParton.Data(), tempVar.Data()));

    TEfficiency *eff_had_orthogonal = (TEfficiency*)inf_had_orthogonal->Get(TString::Format("Efficiency%s_%s",varParton.Data(), tempVar.Data()));
    TEfficiency *eff_sem_orthogonal = (TEfficiency*)inf_sem_orthogonal->Get(TString::Format("Efficiency%s_%s",varParton.Data(), tempVar.Data()));
    TEfficiency *eff_dil_orthogonal = (TEfficiency*)inf_dil_orthogonal->Get(TString::Format("Efficiency%s_%s",varParton.Data(), tempVar.Data()));

    TEfficiency *acc_had = (TEfficiency*)inf_had->Get(TString::Format("Acceptance%s_%s",varParton.Data(), variable[ivar].Data()));
    TEfficiency *acc_sem = (TEfficiency*)inf_sem->Get(TString::Format("Acceptance%s_%s",varParton.Data(), variable[ivar].Data()));
    TEfficiency *acc_dil = (TEfficiency*)inf_dil->Get(TString::Format("Acceptance%s_%s",varParton.Data(), variable[ivar].Data()));

    TEfficiency *acc_had_orthogonal = (TEfficiency*)inf_had_orthogonal->Get(TString::Format("Acceptance%s_%s",varParton.Data(), variable[ivar].Data()));
    TEfficiency *acc_sem_orthogonal = (TEfficiency*)inf_sem_orthogonal->Get(TString::Format("Acceptance%s_%s",varParton.Data(), variable[ivar].Data()));
    TEfficiency *acc_dil_orthogonal = (TEfficiency*)inf_dil_orthogonal->Get(TString::Format("Acceptance%s_%s",varParton.Data(), variable[ivar].Data()));

    efficiency[ivar] = (TEfficiency*)eff_had->Clone();
    *efficiency[ivar] += (*eff_sem);
    *efficiency[ivar] += (*eff_dil); 

    acceptance[ivar] = (TEfficiency*)acc_had->Clone();
    *acceptance[ivar] += (*acc_sem);
    *acceptance[ivar] += (*acc_dil); 

    // orthogonal eff, acc
    efficiency_orthogonal[ivar] = (TEfficiency*)eff_had_orthogonal->Clone();
    *efficiency_orthogonal[ivar] += (*eff_sem_orthogonal);
    *efficiency_orthogonal[ivar] += (*eff_dil_orthogonal); 

    acceptance_orthogonal[ivar] = (TEfficiency*)acc_had_orthogonal->Clone();
    *acceptance_orthogonal[ivar] += (*acc_sem_orthogonal);
    *acceptance_orthogonal[ivar] += (*acc_dil_orthogonal); 

    efficiency_denom[ivar] = (TH1F*)efficiency[ivar]->GetTotalHistogram();
    efficiency_denom_orthogonal[ivar] = (TH1F*)efficiency_orthogonal[ivar]->GetTotalHistogram();

  }
  
  TUnfold *unf[BND_reco.size()];

  TCanvas *can[BND_reco.size()], *can_rho[BND_reco.size()], *canError[BND_reco.size()];
  TH1 *hUnf[BND_reco.size()], *hUnf_complex[BND_reco.size()], *hUnf_lcurve[BND_reco.size()];
  TH1F *hUnfTemp[BND_reco.size()], *hUnfTemp_complex[BND_reco.size()], *hUnfTemp_lcurve[BND_reco.size()];
  TH1F *hUnfFinal[BND_reco.size()], *hUnfFinal_complex[BND_reco.size()], *hUnfFinal_lcurve[BND_reco.size()];

  // this is the file that contains all combined fiducial measurements
  signalFile = TFile::Open(TString::Format("%s/Nominal/combined/HistoReduced_1000_TT.root", year.Data()));

  TFile *outf = TFile::Open(TString::Format("UnfoldedCombined/%s/OutputFile%s%s.root", year.Data(), varParton.Data(), tempFileName.Data()),"RECREATE");

  // read the files from our analysis
  // TFile *inf_ttX = TFile::Open("UnfoldingResults_ttX.root");
  TFile *signalFile_orthogonal = TFile::Open(TString::Format("../%s/Nominal/combined/HistoReduced_1000_TT.root", year.Data()));

  // get the file after acceptance and unfold 
  // TH1F *httX = (TH1F*)inf_ttX->Get("Common_leadingJetPt");
  // get also histogram after unfolding
  // TH1F *httX = (TH1F*)inf_ttX->Get("hWt_jetPt0_2btag");
  // TH1F *httX_CS = (TH1F*)inf_ttX->Get("unfoldedHistogram_leadingJetPt");
  // TH1F *httX_theory = (TH1F*)inf_ttX->Get("theory_leadingJetPt");

  // loop on each variable
  for(int ivar = 0; ivar<BND_reco.size(); ivar++)
  {
    int sizeBins = NBINS[ivar];
    float tempBND[NBINS[ivar]+1];
    std::copy(BND_reco[ivar].begin(), BND_reco[ivar].end(), tempBND);

    float tempBNDGen[NBINS_GEN[ivar]+1];
    std::copy(BND_gen[ivar].begin(), BND_gen[ivar].end(), tempBNDGen);

    //from signal file get the initial S_j with j bins ~ 2* parton bins (i)
    cout<<TString::Format("hWt_%s_2btag",
                                  variable[ivar].Data())<<endl;
    
    hSig_orthogonal[ivar] = (TH1F*)signalFile_orthogonal->Get(TString::Format("hWt_%s_2btag",
                                  variable[ivar].Data()));
    
    hSig[ivar] = (TH1F*)signalFile->Get(TString::Format("hWt_%s_2btag",
                                  variable[ivar].Data()));
    
    TH1F *hSignal_init = (TH1F*)hSig[ivar]->Clone(TString::Format("hSig_init_%s", variable[ivar].Data()));
    TH1F *hSignal_orthogonal_init = (TH1F*)hSig_orthogonal[ivar]->Clone(TString::Format("hSig_init_%s", variable[ivar].Data()));
    //set the new content and get acceptance
    cout<<"The variable is: "<<variable[ivar]<<endl;

    for(int j =1; j<hSig[ivar]->GetNbinsX()+1; j++)
    {
      float oldContent = hSig[ivar]->GetBinContent(j);
      float oldContentError = hSig[ivar]->GetBinError(j);
      float acc = acceptance[ivar]->GetEfficiency(j);
      float accError = (acceptance[ivar]->GetEfficiencyErrorLow(j) + acceptance[ivar]->GetEfficiencyErrorUp(j))/2;

      hSig[ivar]->SetBinContent(j, acc*oldContent);
      
      //error propagation
      if(accError > 0) hSig[ivar]->SetBinError(j,TMath::Sqrt(TMath::Power(accError*oldContent,2) + TMath::Power(oldContentError*acc,2)));
      else hSig[ivar]->SetBinError(j, hSig[ivar]->GetBinError(j));

    }

    for(int j =1; j<hSig_orthogonal[ivar]->GetNbinsX()+1; j++)
    {
      float oldContent_orthogonal = hSig_orthogonal[ivar]->GetBinContent(j);
      float oldContentError_orthogonal = hSig_orthogonal[ivar]->GetBinError(j);
      float acc_orthogonal = acceptance_orthogonal[ivar]->GetEfficiency(j);
      float accError_orthogonal = (acceptance_orthogonal[ivar]->GetEfficiencyErrorLow(j) + acceptance_orthogonal[ivar]->GetEfficiencyErrorUp(j))/2;
      
      hSig_orthogonal[ivar]->SetBinContent(j, acc_orthogonal*oldContent_orthogonal);

      //error propagation
      if(accError_orthogonal > 0) hSig_orthogonal[ivar]->SetBinError(j,
                                TMath::Sqrt(TMath::Power(accError_orthogonal*oldContent_orthogonal,2) + TMath::Power(oldContentError_orthogonal*acc_orthogonal,2)));
      else hSig_orthogonal[ivar]->SetBinError(j, hSig_orthogonal[ivar]->GetBinError(j));
    }
    
    TString tempVar;
    if(isParton)
      tempVar = variableParton[ivar];
    else
      tempVar = variableGen[ivar];


    hUnf[ivar] = new TH1F(TString::Format("hUnf_%s", variable[ivar].Data()), TString::Format("hUnf_%s", variable[ivar].Data()), NBINS_GEN[ivar], tempBNDGen);
    hUnf_complex[ivar] = new TH1F(TString::Format("hUnf_%s", variable[ivar].Data()), TString::Format("hUnf_%s", variable[ivar].Data()), NBINS_GEN[ivar], tempBNDGen);
    
    //unfold with methods:

    hUnf[ivar] = UnfoldDensity_regular(hResponse[ivar], hSig_orthogonal[ivar], variable[ivar], tempVar);
    
    hUnf[ivar]->GetYaxis()->SetTitle(TString::Format("#frac{d#sigma}{d#chi} %s", varParton.Data()));
    hUnf[ivar]->GetYaxis()->SetTitleOffset(1.4);

    // complex global correlation 
    hUnf_complex[ivar] = UnfoldDensity_complex(hResponse[ivar], hSig[ivar], variable[ivar], tempVar);
    
    hUnf_complex[ivar]->GetYaxis()->SetTitle(TString::Format("#frac{d#sigma}{d#chi} %s", varParton.Data()));
    hUnf_complex[ivar]->GetYaxis()->SetTitleOffset(1.4);

    // lcurve method
    hUnf_lcurve[ivar] = UnfoldDensity_lcurve(hResponse[ivar], hSig[ivar], variable[ivar], tempVar);
    
    hUnf_lcurve[ivar]->GetYaxis()->SetTitle(TString::Format("#frac{d#sigma}{d#chi} %s", varParton.Data()));
    hUnf_lcurve[ivar]->GetYaxis()->SetTitleOffset(1.4);

    // this is just after unfolding
    TH1F *hUnf_init = (TH1F*)hUnf[ivar]->Clone(TString::Format("hUnf_init_%s", variable[ivar].Data()));
    TH1F *hUnf_init_complex = (TH1F*)hUnf_complex[ivar]->Clone(TString::Format("hUnf_init_complex_%s)", variable[ivar].Data()));
    TH1F *hUnf_init_lcurve = (TH1F*)hUnf_lcurve[ivar]->Clone(TString::Format("hUnf_init_lcurve_%s)", variable[ivar].Data()));

    PlotHistogramErrors(hUnf_init, hUnf_init_complex, hUnf_init_lcurve);

    // handle efficiency
    for(int i=1; i<hUnf[ivar]->GetNbinsX()+1; i++)
    {
      
      float eff = efficiency[ivar]->GetEfficiency(i);
      // cout<<"eff "<<i<<" "<<eff<<endl;
      /* If eff --> TEfficiency object: GetEfficiency(i); */
      if(eff >0)
      {
        float oldContent = hUnf[ivar]->GetBinContent(i);
        float newContent = hUnf[ivar]->GetBinContent(i)/eff;

        float oldContent_complex = hUnf_complex[ivar]->GetBinContent(i);
        float newContent_complex = hUnf_complex[ivar]->GetBinContent(i)/eff;

        float oldContent_lcurve = hUnf_lcurve[ivar]->GetBinContent(i);
        float newContent_lcurve = hUnf_lcurve[ivar]->GetBinContent(i)/eff;

        float effError = (efficiency[ivar]->GetEfficiencyErrorLow(i) + efficiency[ivar]->GetEfficiencyErrorUp(i))/2;
        
        // regular unfolding without regularization
        float sqrt1 = TMath::Power((1/eff)*hUnf[ivar]->GetBinError(i),2);
        float sqrt2 = TMath::Power((hUnf[ivar]->GetBinContent(i)*effError),2)/TMath::Power(eff,4);
        hUnf[ivar]->SetBinContent(i, newContent);
        hUnf[ivar]->SetBinError(i, TMath::Sqrt(sqrt1+sqrt2));

        // complex unfolding
        float sqrt1_c = TMath::Power((1/eff)*hUnf_complex[ivar]->GetBinError(i),2);
        float sqrt2_c = TMath::Power((hUnf_complex[ivar]->GetBinContent(i)*effError),2)/TMath::Power(eff,4);
        hUnf_complex[ivar]->SetBinContent(i, newContent_complex);
        hUnf_complex[ivar]->SetBinError(i, TMath::Sqrt(sqrt1_c+sqrt2_c));
        

        // lcurve unfolding
        float sqrt1_l = TMath::Power((1/eff)*hUnf_lcurve[ivar]->GetBinError(i),2);
        float sqrt2_l = TMath::Power((hUnf_lcurve[ivar]->GetBinContent(i)*effError),2)/TMath::Power(eff,4);
        hUnf_lcurve[ivar]->SetBinContent(i, newContent_lcurve);
        hUnf_lcurve[ivar]->SetBinError(i, TMath::Sqrt(sqrt1_l+sqrt2_l));
        cout<< " -------- "<<endl;
        cout<< " -------- "<<endl;
        cout<<"regular value ± errror: "<< newContent << " ± "<<TMath::Sqrt(sqrt1+sqrt2)<<endl;
        cout<<"complex value ± error: "<< newContent_complex << " ± "<<TMath::Sqrt(sqrt1_c+sqrt2_c)<<endl;
        cout<<"lcurve value ± error: "<< newContent_lcurve << " ± "<<TMath::Sqrt(sqrt1_l+sqrt2_l)<<endl;

      }
    }

    //draw the unfolded and extrapolated with the mc result
    can[ivar] = new TCanvas(TString::Format("can_%s",variable[ivar].Data()),TString::Format("can_%s",variable[ivar].Data()) , 800,600);
    can[ivar]->cd();

    cout<<"LUMI "<<LUMI<<endl;

    auto *closure_padRatio = new TPad("closure_pad2","closure_pad2",0.,0.,1.,0.3);
    closure_padRatio->Draw();
    closure_padRatio->SetTopMargin(0.05);
    closure_padRatio->SetBottomMargin(0.3);
    closure_padRatio->SetGrid();

    auto *closure_pad1 = new TPad("closure_pad1","closure_pad1",0.,0.3,1.,1.);
    closure_pad1->Draw();
    closure_pad1->SetBottomMargin(0.005);
    closure_pad1->cd();
    cout<< "ok"<<endl;
    
    // get theory from efficiency
    hTheory[ivar] = (TH1F*)efficiency_denom[ivar]->Clone();
    hTheory_orthogonal[ivar] = (TH1F*)efficiency_denom_orthogonal[ivar]->Clone();
    TH1F *hTheory_notScaled = (TH1F*)hTheory[ivar]->Clone(TString::Format("hTheory_%s_notScaled", variable[ivar].Data()));
    hTheory_notScaled->SetTitle(TString::Format("hTheory_%s_notScaled", variable[ivar].Data()));

    // this is diff cross section 
    hTheory[ivar]->Scale(1/LUMI, "width");
    hTheory[ivar]->SetTitle(TString::Format("hTheory_%s", variable[ivar].Data()));

    hTheory_orthogonal[ivar]->Scale(1/LUMI, "width");
    hTheory_orthogonal[ivar]->SetTitle(TString::Format("hTheoryOrthogonal_%s", variable[ivar].Data()));
    
    //this is differential cross section dsigma / dX  = S_i / L * dXi
    hUnf[ivar]->Scale(1/LUMI, "width");
    hUnf_complex[ivar]->Scale(1/LUMI, "width");
    hUnf_lcurve[ivar]->Scale(1/LUMI, "width");

    hUnfFinal[ivar]= (TH1F*)hUnf[ivar]->Clone(TString::Format("hUnfFinal_%s", variable[ivar].Data()));
    hUnfFinal_complex[ivar]= (TH1F*)hUnf_complex[ivar]->Clone(TString::Format("hUnfFinal_%s_complex", variable[ivar].Data()));
    hUnfFinal_lcurve[ivar]= (TH1F*)hUnf_lcurve[ivar]->Clone(TString::Format("hUnfFinal_%s_lcurve", variable[ivar].Data()));
    

    // make them pretty
    hUnfFinal[ivar]->GetYaxis()->SetTitle("#frac{d#sigma}{d#chi}");
    hUnfFinal_complex[ivar]->GetYaxis()->SetTitle("#frac{d#sigma}{d#chi}");
    hUnfFinal_lcurve[ivar]->GetYaxis()->SetTitle("#frac{d#sigma}{d#chi}");
    // httX_CS->GetYaxis()->SetTitle("#frac{d#sigma}{d#chi}");
    hUnfFinal[ivar]->SetLineColor(kBlack);
    hUnfFinal_complex[ivar]->SetLineColor(kMagenta);
    hUnfFinal_lcurve[ivar]->SetLineColor(kRed);
    // hTheory_orthogonal[ivar]->SetLineColor(kGreen+2);

    hUnfFinal[ivar]->SetMarkerStyle(20);
    hUnfFinal[ivar]->SetMarkerColor(kBlack);
    // httX_CS->SetLineColor(kBlack);
    
    hTheory[ivar]->Draw();
    // hTheory_orthogonal[ivar]->Draw("same");
    hUnfFinal[ivar]->Draw("P same");
    hUnfFinal_complex[ivar]->Draw("same");
    hUnfFinal_lcurve[ivar]->Draw("same");



    TLegend* legend = new TLegend(0.5, 0.5, 0.9, 0.9);
    // legend->AddEntry(httX_CS, "Orthogonal, No Regul.", "lep");
    legend->AddEntry(hTheory[ivar], "Theory", "lep");
    legend->AddEntry(hUnfFinal[ivar], "Orthogonal Bins, No Regul.", "lep");
    legend->AddEntry(hUnfFinal_complex[ivar], "Double Bins GlobalCorrelation", "lep");
    legend->AddEntry(hUnfFinal_lcurve[ivar], "Double Bins LCurve", "lep");
    legend->Draw();
 
    closure_padRatio->cd();
    hUnfTemp[ivar] = (TH1F*)hUnfFinal[ivar]->Clone(TString::Format("hUnf_%s", variable[ivar].Data()));
    hUnfTemp_complex[ivar] = (TH1F*)hUnfFinal_complex[ivar]->Clone(TString::Format("hUnf_complex_%s", variable[ivar].Data()));
    hUnfTemp_lcurve[ivar] = (TH1F*)hUnfFinal_lcurve[ivar]->Clone(TString::Format("hUnf_lcurve_%s", variable[ivar].Data()));
    
    
    hUnfTemp[ivar]->Divide(hTheory_orthogonal[ivar]);
    hUnfTemp_complex[ivar]->Divide(hTheory[ivar]);
    hUnfTemp_lcurve[ivar]->Divide(hTheory[ivar]);
    // hUnfTemp_ttX->Divide(httX_theory);

    hUnfTemp[ivar]->SetTitle("");
    hUnfTemp[ivar]->GetYaxis()->SetTitle("#frac{Unfolded}{Theory}");
    hUnfTemp[ivar]->GetXaxis()->SetTitle(variable[ivar].Data());
    hUnfTemp[ivar]->GetYaxis()->SetTitleSize(14);
    hUnfTemp[ivar]->GetYaxis()->SetTitleFont(43);
    hUnfTemp[ivar]->GetYaxis()->SetTitleOffset(1.55);
    hUnfTemp[ivar]->GetYaxis()->SetLabelFont(43);
    hUnfTemp[ivar]->GetYaxis()->SetLabelSize(15);
    hUnfTemp[ivar]->GetXaxis()->SetTitleSize(0.09);
    hUnfTemp[ivar]->GetYaxis()->SetRangeUser(0,2);


    hUnfTemp_lcurve[ivar]->SetTitle("");
    hUnfTemp_lcurve[ivar]->GetYaxis()->SetTitle("#frac{Unfolded}{Theory}");
    hUnfTemp_lcurve[ivar]->GetXaxis()->SetTitle(variable[ivar].Data());
    hUnfTemp_lcurve[ivar]->GetYaxis()->SetTitleSize(14);
    hUnfTemp_lcurve[ivar]->GetYaxis()->SetTitleFont(43);
    hUnfTemp_lcurve[ivar]->GetYaxis()->SetTitleOffset(1.55);
    hUnfTemp_lcurve[ivar]->GetYaxis()->SetLabelFont(43);
    hUnfTemp_lcurve[ivar]->GetYaxis()->SetLabelSize(15);
    hUnfTemp_lcurve[ivar]->GetXaxis()->SetTitleSize(0.09);
    hUnfTemp_lcurve[ivar]->GetYaxis()->SetRangeUser(0,2);
    hUnfTemp_lcurve[ivar]->SetLineColor(kBlack);
    hUnfTemp_lcurve[ivar]->SetMarkerStyle(23);
    hUnfTemp_lcurve[ivar]->SetMarkerColor(kBlack);


    hUnfTemp_complex[ivar]->SetTitle("");
    hUnfTemp_complex[ivar]->GetYaxis()->SetTitle("#frac{Unfolded}{Theory}");
    hUnfTemp_complex[ivar]->GetXaxis()->SetTitle(variable[ivar].Data());
    hUnfTemp_complex[ivar]->GetYaxis()->SetTitleSize(14);
    hUnfTemp_complex[ivar]->GetYaxis()->SetTitleFont(43);
    hUnfTemp_complex[ivar]->GetYaxis()->SetTitleOffset(1.55);
    hUnfTemp_complex[ivar]->GetYaxis()->SetLabelFont(43);
    hUnfTemp_complex[ivar]->GetYaxis()->SetLabelSize(15);
    hUnfTemp_complex[ivar]->GetXaxis()->SetTitleSize(0.09);
    hUnfTemp_complex[ivar]->GetYaxis()->SetRangeUser(0,2);
    hUnfTemp_complex[ivar]->SetLineColor(kMagenta);
    hUnfTemp_complex[ivar]->SetMarkerStyle(21);
    hUnfTemp_complex[ivar]->SetMarkerColor(kMagenta);
    
    hUnfTemp[ivar]->SetLineColor(kBlack);
    hUnfTemp[ivar]->SetMarkerStyle(20);
    hUnfTemp[ivar]->SetMarkerColor(kBlack);
    
    hUnfTemp[ivar]->Draw();
    hUnfTemp_complex[ivar]->Draw("same");
    hUnfTemp_lcurve[ivar]->Draw("same");
    hUnfTemp[ivar]->GetXaxis()->SetLabelSize(0.09);
    hUnfTemp_complex[ivar]->GetXaxis()->SetLabelSize(0.09);
    hUnfTemp_lcurve[ivar]->GetXaxis()->SetLabelSize(0.09);

    lumi_13TeV = TString::Format("%0.1f fb^{-1}", luminosity["luminosityAll"]/1000);
    int iPeriod = 13;
    int iPos = 0;
    extraTextFactor = 0.14;
    writeExtraText=true;
    CMS_lumi(closure_pad1, "combined", iPos);
  
  }
  
  // signalFile->Close();
  // signalFile_orthogonal->Close();

}



void PlotHistogramErrors(TH1F* hist1, TH1F *hist2, TH1F *hist3)
{
    // Create a canvas
    TCanvas* canvas = new TCanvas("canvas", "Histogram Errors", 800, 600);

    // Set the plotting style
    gStyle->SetOptStat(0);

    // Set histogram attributes
    hist1->SetMarkerStyle(20);
    hist1->SetMarkerSize(0.8);
    hist1->SetLineColor(kRed);
    // hist1->SetFillColorAlpha(kRed, 0.2); // Set histogram fill color with transparency

    hist2->SetMarkerStyle(21);
    hist2->SetMarkerSize(0.8);
    hist2->SetLineColor(kMagenta);
    // hist2->SetFillColorAlpha(kMagenta, 0.2);

    hist3->SetMarkerStyle(23);
    hist3->SetMarkerSize(0.8);
    hist3->SetLineColor(kBlack);
    // hist3->SetFillColorAlpha(kBlack, 0.2);

    // Create a histogram with only the error bars
    TH1F* histErrors1 = (TH1F*)hist1->Clone();
    histErrors1->Reset();
    histErrors1->GetYaxis()->SetTitle("Error");
    histErrors1->SetMarkerSize(0);
    histErrors1->SetFillColor(0);
    histErrors1->SetFillStyle(0);

    TH1F* histErrors2 = (TH1F*)hist2->Clone();
    histErrors2->Reset();
    histErrors2->GetYaxis()->SetTitle("Error");
    histErrors2->SetMarkerSize(0);
    histErrors2->SetFillColor(0);
    histErrors2->SetFillStyle(0);

    TH1F* histErrors3 = (TH1F*)hist3->Clone();
    histErrors3->Reset();
    histErrors3->GetYaxis()->SetTitle("Error");
    histErrors3->SetMarkerSize(0);
    histErrors3->SetFillColor(0);
    histErrors3->SetFillStyle(0);

    // Set only the errors
    for (int i = 1; i <= hist1->GetNbinsX(); i++) {
        double binError = hist1->GetBinError(i);
        histErrors1->SetBinContent(i, binError);

        binError = hist2->GetBinError(i);
        histErrors2->SetBinContent(i, binError);

        binError = hist3->GetBinError(i);
        histErrors3->SetBinContent(i, binError);

    }

    // Draw the histograms with errors
    histErrors1->Draw();
    histErrors2->Draw("SAME");
    histErrors3->Draw("SAME");

    // Create a legend
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(histErrors1, "Regular", "lep");
    legend->AddEntry(histErrors2, "GlobalCorrelation", "lep");
    legend->AddEntry(histErrors3, "LCurve", "lep");
    legend->Draw();

    // Update the canvas
    canvas->Update();
    canvas->Print("error_propagation_unfolding.pdf", "pdf");
}