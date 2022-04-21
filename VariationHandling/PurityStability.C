#include "../CMS_plots/tdrstyle.C"
#include "../CMS_plots/CMS_lumi.C"

using namespace std;
using std::cin;
using std::cout;
using std::endl;

#include "TemplateConstants.h"

TString varParton;

void PurityStability(TString isParton = "true")
{
  initFilesMapping();
  //setTDRStyle();
  gStyle->SetOptStat(0);

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
  const int NVAR = 10;

  TH1F *inSig[BND_reco.size()], *hSig[BND_reco.size()];
  TString variable[NVAR]   = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1", "chi", "cosTheta_0", "cosTheta_1"};
  TString variableParton[NVAR] = {"mTTbarParton", "ptTTbarParton", "yTTbarParton","partonPt0", "partonPt1", "partonY0", "partonY1","chiParton", "cosThetaParton_0", "cosThetaParton_1"};
  TString variableGen[NVAR] = {"mJJGen", "ptJJGen", "yJJGen","genjetPt0", "genjetPt1", "genjetY0", "genjetY1", "chiParticle", "cosThetaParticle_0", "cosThetaParticle_1"};

  float LUMI = 0;


  TH1F *purity_had[4][BND_reco.size()], *purity_semi[4][BND_reco.size()], *purity_dil[4][BND_reco.size()];
  TH1F *stability_had[4][BND_reco.size()], *stability_semi[4][BND_reco.size()], *stability_dil[4][BND_reco.size()];

  TH1F *purity[4], *stability[4];

  TString years[4] = {"2016_preVFP", "2016_postVFP", "2017", "2018"};
  for (int iy=0; iy<4; iy++)
  { 
    TFile *inf_had = TFile::Open(TString::Format("%s/ResponsesNominal/ResponsesEfficiency_TTToHadronic.root", years[iy].Data()));
    TFile *inf_sem = TFile::Open(TString::Format("%s/ResponsesNominal/ResponsesEfficiency_TTToSemiLeptonic.root", years[iy].Data()));
    TFile *inf_dil = TFile::Open(TString::Format("%s/ResponsesNominal/ResponsesEfficiency_TTTo2L2Nu.root", years[iy].Data()));

    LUMI += luminosity["luminosity"+years[iy]];
    //get response matrix
    for(int ivar = 0; ivar<NVAR; ivar++)
    {
        TString tempVar;
        if(isParton)
            tempVar = variableParton[ivar];
        else
            tempVar = variableGen[ivar];
        
        stability_had[iy][ivar] = (TH1F*)inf_had->Get(TString::Format("Stability%s_%s",varParton.Data(), variable[ivar].Data()));
        stability_semi[iy][ivar] = (TH1F*)inf_sem->Get(TString::Format("Stability%s_%s",varParton.Data(), variable[ivar].Data()));
        stability_dil[iy][ivar] = (TH1F*)inf_dil->Get(TString::Format("Stability%s_%s",varParton.Data(), variable[ivar].Data()));
        
        purity_had[iy][ivar] = (TH1F*)inf_had->Get(TString::Format("Purity%s_%s",varParton.Data(), variable[ivar].Data()));
        purity_semi[iy][ivar] = (TH1F*)inf_sem->Get(TString::Format("Purity%s_%s",varParton.Data(), variable[ivar].Data()));
        purity_dil[iy][ivar] = (TH1F*)inf_dil->Get(TString::Format("Purity%s_%s",varParton.Data(), variable[ivar].Data()));

        // stabilities are already scaled to XSEC and lumi
        if (iy == 0) {
            stability[ivar] = (TH1F*)stability_had[iy][ivar]->Clone();
            purity[ivar] = (TH1F*)purity_had[iy][ivar]->Clone();
            cout<<"stab0 "<<stability[ivar]->GetName()<<endl;
            cout<<"purity 0 "<< purity[ivar]->GetName()<<endl;
        }
        else {
            stability[ivar] ->Add(stability_had[iy][ivar]);
            //purity[ivar] ->Add(purity_had[iy][ivar]);

        }
        stability[ivar]->Add(stability_semi[iy][ivar]);
        stability[ivar]->Add(stability_dil[iy][ivar]); 
        //purity[ivar]->Add(purity_semi[iy][ivar]);
        //purity[ivar]->Add(purity_dil[iy][ivar]);

        //cout<< "year "<<years[iy]<< ", ivar "<<ivar<< " "<<variable[ivar]<<endl;
        //cout<<"stabilities: " <<stability_had[iy][ivar]->GetName()<< endl;
        //cout<<"purity: " <<purity_had[iy][ivar]->GetName()<< endl;
    }
  }
    
  lumi_13TeV = TString::Format("%0.1f fb^{-1}", luminosity["luminosityAll"]/1000);
  int iPeriod = 13;
  int iPos = 0;
  //extraTextFactor = 0.14;
  writeExtraText=true;
  
  for(int ivar = 0; ivar<BND_reco.size(); ivar++)
  {

    // plot purity and stability in 1 plot for the combined
    TCanvas *can_purStab = new TCanvas(TString::Format("canEff_%s",variable[ivar].Data()),
                                    TString::Format("canEff_%s",variable[ivar].Data()) , 800,600);
    
    can_purStab->cd();
    TLegend *leg = new TLegend(0.7, 0.75, 0.9, 0.9);
    leg->AddEntry(stability[ivar], "Stability", "lep");
    leg->AddEntry(purity[ivar], "Purity", "lep");
    stability[ivar]->Draw();
    purity[ivar]->Draw("same");

    purity[ivar]->SetLineColor(kBlack);
    stability[ivar]->SetLineColor(kBlue);
    stability[ivar]->GetYaxis()->SetRangeUser(0.0, 1);
    gPad->Update(); 
    
    CMS_lumi(can_purStab, "combined", iPos);
    leg->Draw();
    can_purStab->SaveAs(TString::Format("Comparison_EffAccResponses/combined/PurityStability%s_%s.pdf", 
                    varParton.Data(), variable[ivar].Data()), "pdf");
  }
}