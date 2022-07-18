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
#include "../CMS_plots/CMS_lumi.C"
#include "../CMS_plots/tdrstyle.C"

using namespace std;

#include "UseTUnfoldDensity.cpp"

using std::cin;
using std::cout;
using std::endl;

#include "TemplateConstants.h"

TString varParton;

void Unfold_Combined(TString dir, TString inputFile, TString isParton = "true")
{
  bool isNorm = false;
  initFilesMapping();
  //setTDRStyle();
  gStyle->SetOptStat(0);

  TString tempFileName;
  if(dir.EqualTo("PDFWeights")) tempFileName = "_pdf_"+inputFile;
  else if (dir.EqualTo("ScaleWeights")) tempFileName = "_scale_"+inputFile;
  else if (dir.EqualTo("Nominal")) tempFileName = "";
  else tempFileName = "_"+inputFile;

  std::vector< std::vector <Float_t> > const BND_gen = {{1000, 1200, 1400, 1600, 1800, 2000, 2400, 3000, 5000}, //mJJ
                                                        {0, 60, 150, 300, 450, 850, 1300}, //ptjj
                                                        {-2.4, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.4}, //yjj
                                                        {450, 500, 570, 650, 800, 1100, 1500}, //jetpt0
                                                        {400, 450, 500, 570, 650, 800, 1100, 1500}, //jetpt1
                                                        {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.4}, //jetY0
                                                        {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.4}, //jetY1
                                                        {1,2,3,4,5,6,7,8,9,10,13,16}, //chi
                                                        {0,0.2,0.4,0.6,0.7,0.8,1}, //|cosTheta*| leading
                                                        {0,0.2,0.4,0.6,0.7,0.8,1}}; //|cosTheta*| subleading

  std::vector< std::vector <Float_t> > const BND_reco ={{1000, 1200, 1400, 1600, 1800, 2000, 2400, 3000, 5000}, //mJJ
                                                        {0, 60, 150, 300, 450, 850, 1300}, //ptJJ
                                                        {-2.4, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.4}, //yjj
                                                        {450, 500, 570, 650, 800, 1100, 1500}, //jetPt0
                                                        {400, 450, 500, 570, 650, 800, 1100, 1500}, //jetPt1
                                                        {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.4}, //jetY0
                                                        {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.4}, //jetY1
                                                        {1,2,3,4,5,6,7,8,9,10,13,16}, //chi
                                                        {0,0.2,0.4,0.6,0.7,0.8,1}, //|cosTheta*| leading
                                                        {0,0.2,0.4,0.6,0.7,0.8,1}}; //|cosTheta*| subleading


  //get the files:
  //1. the signal file has the fiducial measurements that are going to be used as input
  TFile *signalFile;

  //whether parton or particle, from the choice of the user
  varParton = "Parton";
  if(isParton.EqualTo("false")) varParton = "Particle";

  //get the number of bins for each
  int NBINS[BND_reco.size()];
  int NBINS_GEN[BND_gen.size()];
  const int NVAR = 10;
  for (int i = 0; i<BND_reco.size(); i++)
  	NBINS[i] = BND_reco[i].size()-1;
  for (int i = 0; i<BND_gen.size(); i++)
  	NBINS_GEN[i] = BND_gen[i].size()-1;

  TH1F *inSig[BND_reco.size()], *hSig[BND_reco.size()];
  TString variable[NVAR]   = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1", "chi", "cosTheta_0", "cosTheta_1"};
  TString variableParton[NVAR] = {"mTTbarParton", "ptTTbarParton", "yTTbarParton","partonPt0", "partonPt1", "partonY0", "partonY1","chiParton", "cosThetaParton_0", "cosThetaParton_1"};
  TString variableGen[NVAR] = {"mJJGen", "ptJJGen", "yJJGen","genjetPt0", "genjetPt1", "genjetY0", "genjetY1", "chiParticle", "cosThetaParticle_0", "cosThetaParticle_1"};
  //2. This file has the response matrices as well as the efficiency and acceptance for the signal procedure
  //handle responses efficiency, acceptance
  TH2F *hResponse[BND_reco.size()];
  TH1F *efficiency[BND_reco.size()], *acceptance[BND_reco.size()];
  TH1F *efficiency_denom[BND_reco.size()];
  TH1F *fiducial_theory[BND_reco.size()];
  TH1F *hFidTheory[BND_reco.size()], *hFidTheory_norm[BND_reco.size()];
  TH1F *hTheory[BND_reco.size()], *hTheory_norm[BND_reco.size()];
  float LUMI = 0;

  TEfficiency *acc_had[BND_reco.size()], *acc_sem[BND_reco.size()], *acc_dil[BND_reco.size()];
  TEfficiency *eff_had[BND_reco.size()], *eff_sem[BND_reco.size()], *eff_dil[BND_reco.size()];

  TH2F *hResponse_had[4][BND_reco.size()], *hResponse_sem[4][BND_reco.size()], *hResponse_dil[4][BND_reco.size()];

  TString years[4] = {"2016_preVFP", "2016_postVFP", "2017", "2018"};
  for (int iy=0; iy<4; iy++)
  {
    cout<<TString::Format("%s/Responses%s/ResponsesEfficiency_TTToHadronic%s.root", years[iy].Data(), dir.Data(), tempFileName.Data())<<endl;
    TFile *inf_had = TFile::Open(TString::Format("%s/Responses%s/ResponsesEfficiency_TTToHadronic%s.root", years[iy].Data(), dir.Data(), tempFileName.Data()));
    TFile *inf_sem = TFile::Open(TString::Format("%s/Responses%s/ResponsesEfficiency_TTToSemiLeptonic%s.root", years[iy].Data(), dir.Data(), tempFileName.Data()));
    TFile *inf_dil = TFile::Open(TString::Format("%s/Responses%s/ResponsesEfficiency_TTTo2L2Nu%s.root", years[iy].Data(), dir.Data(), tempFileName.Data()));

    LUMI += luminosity["luminosity"+years[iy]];
    //get response matrix
    for(int ivar = 0; ivar<BND_reco.size(); ivar++)
    {
        TString tempVar;
        if(isParton)
            tempVar = variableParton[ivar];
        else
            tempVar = variableGen[ivar];

        hResponse_had[iy][ivar] = (TH2F*)inf_had->Get(TString::Format("h%sResponse_%s",varParton.Data(), variable[ivar].Data()));
        hResponse_sem[iy][ivar] = (TH2F*)inf_sem->Get(TString::Format("h%sResponse_%s",varParton.Data(), variable[ivar].Data()));
        hResponse_dil[iy][ivar] = (TH2F*)inf_dil->Get(TString::Format("h%sResponse_%s",varParton.Data(), variable[ivar].Data()));

        // add all (regardless of ttbar process) to the hResponse per variable
        // responses are already scaled to XSEC and lumi
        if (iy == 0) hResponse[ivar] = (TH2F*)hResponse_had[iy][ivar]->Clone();
        else hResponse[ivar] ->Add(hResponse_had[iy][ivar]);
        hResponse[ivar]->Add(hResponse_sem[iy][ivar]);
        hResponse[ivar]->Add(hResponse_dil[iy][ivar]);

        /*
        cout<<"----"<<variable[ivar]<<"----"<<endl;
        cout<<"Integral responses: "<<hResponse[ivar]->Integral()<<endl;
        cout<<"Entries responses: "<<hResponse[ivar]->GetEntries()<<endl;
        */

       // read efficiency and acceptance from combined directories
      TFile *acc_inf, *eff_inf;
      if (dir.EqualTo("Nominal"))
      {
        acc_inf = TFile::Open(TString::Format("AcceptanceCombined/%s/CombAcceptance%s_%s_ResponsesEfficiency_TTToHadronic.root",
                                    dir.Data(), varParton.Data(), variable[ivar].Data()));

        eff_inf = TFile::Open(TString::Format("EfficiencyCombined/%s/CombEfficiency%s_%s_ResponsesEfficiency_TTToHadronic.root",
                                    dir.Data(), varParton.Data(), variable[ivar].Data()));
      }
      else
      {
        //cout<< TString::Format("AcceptanceCombined/%s/CombAcceptance%s_%s_ResponsesEfficiency_TTToHadronic%s.root",
        //                            dir.Data(), varParton.Data(), variable[ivar].Data(), tempFileName.Data()) <<endl;

        acc_inf = TFile::Open(TString::Format("AcceptanceCombined/%s/CombAcceptance%s_%s_ResponsesEfficiency_TTToHadronic%s.root",
                                    dir.Data(), varParton.Data(), variable[ivar].Data(), tempFileName.Data()));

        eff_inf = TFile::Open(TString::Format("EfficiencyCombined/%s/CombEfficiency%s_%s_ResponsesEfficiency_TTToHadronic%s.root",
                                    dir.Data(), varParton.Data(), variable[ivar].Data(), tempFileName.Data()));
      }

      
      acceptance[ivar] = (TH1F*)acc_inf->Get("acceptance");
      efficiency[ivar] = (TH1F*)eff_inf->Get("efficiency");
      efficiency_denom[ivar] = (TH1F*)eff_inf->Get("denominator_efficiency");
      fiducial_theory[ivar] = (TH1F*)acc_inf->Get("fiducial_theory");

    }
  }


  TUnfold *unf[BND_reco.size()];

  TCanvas *can[BND_reco.size()], *can_rho[BND_reco.size()], *canError[BND_reco.size()];
  TH1 *hUnf[BND_reco.size()];
  TH1F *hUnfTemp[BND_reco.size()];
  TH1F *hUnfFinal[BND_reco.size()];

  //if (dir.EqualTo("Nominal")) inputFile = "";
  TFile *outf = TFile::Open(TString::Format("UnfoldedCombined/%s/OutputFile%s%s.root", dir.Data(), varParton.Data(), tempFileName.Data()),"RECREATE");
  

  for(int ivar = 0; ivar<BND_reco.size(); ivar++)
  {
    // this is the file that contains all combined fiducial measurements
    signalFile = TFile::Open("testFile_.root");

    int sizeBins = NBINS[ivar];
    float tempBND[NBINS[ivar]+1];
    std::copy(BND_reco[ivar].begin(), BND_reco[ivar].end(), tempBND);

    float tempBNDGen[NBINS_GEN[ivar]+1];
    std::copy(BND_gen[ivar].begin(), BND_gen[ivar].end(), tempBNDGen);

    //from signal file get the initial S_j with j bins ~ 2* parton bins (i)

    TString tempFileName_signal;
    if(dir.EqualTo("PDFWeights")) tempFileName_signal = "_pdfVariation"+inputFile;
    else if (dir.EqualTo("ScaleWeights")) tempFileName_signal = "_scaleWeight"+inputFile;
    else if (dir.EqualTo("Nominal")) tempFileName_signal = "";
    else if (dir.EqualTo("bTagVariation"))
    {
      if (inputFile.Contains("up")) inputFile = "up";
      else if (inputFile.Contains("down")) inputFile = "down";
      tempFileName_signal = "_bTag"+inputFile;
    } 
    else tempFileName_signal = "_"+inputFile;

    if (dir.EqualTo("Nominal"))
    {
      hSig[ivar] = (TH1F*)signalFile->Get(TString::Format("combined_%s",
                                    variable[ivar].Data()));
    }
    else 
    {
      cout << TString::Format("combined%s_%s",
                                    tempFileName_signal.Data(),
                                    variable[ivar].Data()) << endl; 

      hSig[ivar] = (TH1F*)signalFile->Get(TString::Format("combined%s_%s",
                                    tempFileName_signal.Data(),
                                    variable[ivar].Data()));
    }
    
    TH1F *hSignal_init = (TH1F*)hSig[ivar]->Clone(TString::Format("hSig_init_%s", variable[ivar].Data()));
    cout <<"before "<< hSignal_init->Integral() <<endl;
    //set the new content and get acceptance
    cout<<"The variable is: "<<variable[ivar]<<endl;

    for(int j =1; j<hSig[ivar]->GetNbinsX()+1; j++)
    {
      float oldContent = hSig[ivar]->GetBinContent(j);
      float oldContentError = hSig[ivar]->GetBinError(j);
      float acc = acceptance[ivar]->GetBinContent(j);
      float accError = acceptance[ivar]->GetBinError(j);
      /*here we have code if we have acceptance as TEfficiency object
      (acceptance[ivar]->GetEfficiencyErrorLow(j) + acceptance[ivar]->GetEfficiencyErrorUp(j))/2;*/

      hSig[ivar]->SetBinContent(j, acc*oldContent);
      //error propagation
      if(accError > 0) hSig[ivar]->SetBinError(j,TMath::Sqrt(TMath::Power(accError*oldContent,2) + TMath::Power(oldContentError*acc,2)));
      else hSig[ivar]->SetBinError(j, hSig[ivar]->GetBinError(j));

    }
    TH1F *hSignal_acceptance = (TH1F*)hSig[ivar]->Clone(TString::Format("hSig_acceptance_%s", variable[ivar].Data()));
    
    TString tempVar;
    if(isParton)
      tempVar = variableParton[ivar];
    else
      tempVar = variableGen[ivar];


    hUnf[ivar] = new TH1F(TString::Format("hUnf_%s", variable[ivar].Data()), TString::Format("hUnf_%s", variable[ivar].Data()), NBINS_GEN[ivar], tempBNDGen);
    hUnf[ivar] = UnfoldDensity(hResponse[ivar], hSig[ivar], variable[ivar], tempVar);
    hUnf[ivar]->GetYaxis()->SetTitle(TString::Format("#frac{d#sigma}{d#chi} %s", varParton.Data()));
    hUnf[ivar]->GetYaxis()->SetTitleOffset(1.4);

    TH1F *hUnf_init = (TH1F*)hUnf[ivar]->Clone(TString::Format("hUnf_init_%s", variable[ivar].Data()));

    for(int i=1; i<hUnf[ivar]->GetNbinsX()+1; i++)
    {
      float eff = efficiency[ivar]->GetBinContent(i);
      /* If eff --> TEfficiency object: GetEfficiency(i); */
      if(eff >0)
      {
        float oldContent = hUnf[ivar]->GetBinContent(i);
	    float newContent = hUnf[ivar]->GetBinContent(i)/eff;

	    float effError = efficiency[ivar]->GetBinError(i);
        /*If eff --> TEfficiency object:
        (efficiency[ivar]->GetEfficiencyErrorLow(i) + efficiency[ivar]->GetEfficiencyErrorUp(i))/2; */

        float sqrt1 = TMath::Power((1/eff)*hUnf[ivar]->GetBinError(i),2);
        float sqrt2 = TMath::Power((hUnf[ivar]->GetBinContent(i)*effError),2)/TMath::Power(eff,4);
        hUnf[ivar]->SetBinContent(i, newContent);
        hUnf[ivar]->SetBinError(i, TMath::Sqrt(sqrt1+sqrt2));
  	  }
    }

    //draw the unfolded and extrapolated with the mc result
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
    cout<< "ok"<<endl;
    // get theory from efficiency
    hTheory[ivar] = (TH1F*)efficiency_denom[ivar]->Clone();
    TH1F *hTheory_notScaled = (TH1F*)hTheory[ivar]->Clone(TString::Format("hTheory_%s_notScaled", variable[ivar].Data()));
    hTheory_notScaled->SetTitle(TString::Format("hTheory_%s_notScaled", variable[ivar].Data()));
    
    //this is fiducial diff cross section 
    hFidTheory[ivar] = (TH1F*)fiducial_theory[ivar]->Clone();
    hFidTheory[ivar]->Scale(1/LUMI, "width");
    hFidTheory[ivar]->SetTitle(TString::Format("hFidTheory_%s", variable[ivar].Data()));

    // this is Fiducial normalized cross section 
    hFidTheory_norm[ivar] = (TH1F*)hFidTheory[ivar]->Clone(TString::Format("hFidTheoryNorm_%s", variable[ivar].Data()));
    hFidTheory_norm[ivar]->Scale(1/hFidTheory[ivar]->Integral());
    hFidTheory_norm[ivar]->SetTitle(TString::Format("hFidTheoryNorm_%s", variable[ivar].Data()));

    // this is diff cross section 
    hTheory[ivar]->Scale(1/LUMI, "width");
    hTheory[ivar]->SetTitle(TString::Format("hTheory_%s", variable[ivar].Data()));
    // this is normalized cross section 
    hTheory_norm[ivar] = (TH1F*)hTheory[ivar]->Clone(TString::Format("hTheoryNorm_%s", variable[ivar].Data()));
    hTheory_norm[ivar]->Scale(1/hTheory_norm[ivar]->Integral());
    hTheory_norm[ivar]->SetTitle(TString::Format("hTheoryNorm_%s", variable[ivar].Data()));
    
    TH1F *hUnfFinal_notScaled= (TH1F*)hUnf[ivar]->Clone(TString::Format("hUnfFinal_%s_notScaled", variable[ivar].Data()));
    //this is differential cross section dsigma / dX  = S_i / L * dXi
    hUnf[ivar]->Scale(1/LUMI, "width");
    hUnfFinal[ivar]= (TH1F*)hUnf[ivar]->Clone(TString::Format("hUnfFinal_%s", variable[ivar].Data()));
    

    hUnf[ivar]->Scale(LUMI/hUnfFinal[ivar]->Integral());
    //hUnf[ivar]->SetTitle(TString::Format("%s Unfolded vs Theory %s",varParton.Data(),variable[ivar].Data()));

    TH1F *hUnf_norm = (TH1F*)hUnf[ivar]->Clone(TString::Format("hUnfNorm_%s", variable[ivar].Data()));
    hUnf_norm->Scale(1./hUnf[ivar]->Integral());
    /*
    cout<<"-----Integral------"<<endl;
    Double_t error, integral;
    integral = hUnf[ivar]->IntegralAndError(1,hUnf[ivar]->GetNbinsX(), error);
    cout<<variable[ivar].Data()<<" hUnfolded: "<<integral<<" ± "<<error<<endl; */

    hUnf[ivar]->GetYaxis()->SetTitle("#frac{d#sigma}{d#chi}");
    hUnf_norm->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d#chi}");

  	hUnf[ivar]->Draw("same");

  	closure_padRatio->cd();
  	hUnfTemp[ivar] = (TH1F*)hUnf[ivar]->Clone(TString::Format("hUnf_%s", variable[ivar].Data()));

  	hUnfTemp[ivar]->Draw();

  	hUnfTemp[ivar]->SetTitle("");
    hUnfTemp[ivar]->GetYaxis()->SetTitle("#frac{Unfolded}{Theory}");
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

    lumi_13TeV = TString::Format("%0.1f fb^{-1}", LUMI/1000);
    //lumi_sqrtS = "13 TeV";
    int iPeriod = 4;
    int iPos = 0;
    writeExtraText=true;
    CMS_lumi(closure_pad1, iPeriod, iPos);

    cout <<"after "<< hSignal_init->Integral() <<endl;
    outf->cd();
    // fiducial signal
    hSignal_init->Write(TString::Format("hSigInit_%s", variable[ivar].Data()));

    // fiducial times acceptance 
    hSignal_acceptance->Write(TString::Format("hSigAcceptance_%s", variable[ivar].Data()));

    // unfolded before efficiency 
    hUnf_init->Write(TString::Format("hUnfInit_%s", variable[ivar].Data()));

    //unfolded after efficiency (not scaled)
    hUnfFinal_notScaled->Write(TString::Format("hUnfoldFinal_NoScaled_%s", variable[ivar].Data()));

    //unfolded after efficiency (scaled to LUMI, width)
    hUnfFinal[ivar]->Write(TString::Format("hUnfoldFinal_%s", variable[ivar].Data()));
    
    // normalised result
    hUnf_norm->Write(TString::Format("hUnfoldNorm_%s", variable[ivar].Data()));
    
    // fiducial theory bulk
    fiducial_theory[ivar]->Write(TString::Format("hFidTheory_NoScaled_%s", variable[ivar].Data()));

    // Fiducial theory diff xsec
    hFidTheory[ivar]->Write(TString::Format("hFidTheory_%s", variable[ivar].Data()));

    // Fiducial theory norm diff xsec 
    hFidTheory_norm[ivar]->Write(TString::Format("hFidTheoryNorm_%s", variable[ivar].Data())); 

    // theory bulk 
    hTheory_notScaled->Write(TString::Format("hTheory_NoScaled_%s", variable[ivar].Data()));

    // theory diff xsec
    hTheory[ivar]->Write(TString::Format("hTheory_%s", variable[ivar].Data()));

    // theory norm diff xsec 
    hTheory_norm[ivar]->Write(TString::Format("hTheoryNorm_%s", variable[ivar].Data())); 

    signalFile->Close();
  }

}