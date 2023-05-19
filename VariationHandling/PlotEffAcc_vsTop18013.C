/*
Plot efficiencies and responses for the 2016preVFP, 2016postVFP, 2017 and 2018
For the Response matrices something like chi2 testing is needed.
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
#include "../CMS_plots/tdrstyle.C"
#include "../CMS_plots/CMS_lumi.C"
#include "FinalResultsConstants.h"
using namespace std;


using std::cin;
using std::cout;
using std::endl;

#include "TemplateConstants.h"


void PlotEffAcc_vsTop18013(bool isParton = false)
{   
    /* init variables and values */
    initFilesMapping();
    AnalysisConstants::initConstants();
    //setTDRStyle();
    gStyle->SetOptStat(0);
    const int NVAR = 10;
    TString variable[NVAR]   = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1", "chi", "cosTheta_0", "cosTheta_1"};
    TString variableParton[NVAR] = {"mTTbarParton", "ptTTbarParton", "yTTbarParton","partonPt0", "partonPt1", "partonY0", "partonY1","chiParton", "cosThetaParton_0", "cosThetaParton_1"};
    TString variableGen[NVAR] = {"mJJGen", "ptJJGen", "yJJGen","genjetPt0", "genjetPt1", "genjetY0", "genjetY1", "chiParticle", "cosThetaParticle_0", "cosThetaParticle_1"};
    TString expr[NVAR] = {"(GeV)", "(GeV)", "", "(GeV)", "(GeV)", "", "", "", "", ""};
    
    std::vector<int> color = {kRed, kGreen, kBlue, kBlack};
    TString varParton = "Parton";
    if(!isParton) varParton = "Particle";
    /* Loop on all variables and create canvases */

    TString years[4] = {"2016_preVFP", "2016_postVFP", "2017", "2018"};
    
    TString variations[6] = {"bTagVariation", "topTaggingVariation", "JES", "PDFWeights", "PSWeights", "ScaleWeights"};
    TH1D *purityParton[NVAR], *stabilityParton[10];

    for (int ivar=0; ivar<NVAR-3; ivar++)
    {
        // define efficiencies for each process for each year 
        TEfficiency *acc_had[4], *acc_sem[4], *acc_dil[4];
        TEfficiency *eff_had[4], *eff_sem[4], *eff_dil[4];
        // define responses per process for each year
        TH2F *hResponse_had[4], *hResponse_sem[4], *hResponse_dil[4];

        // define efficiencies total per year
        TEfficiency *acceptance[4], *efficiency[4];
        // define th2 response matrices per year
        TH2F *hResponse[4];
        
        // legend 
        TLegend *leg_acc = new TLegend(0.7, 0.25, 0.9, 0.45);
        TLegend *leg = new TLegend(0.7, 0.75, 0.9, 0.9);


        TFile *f_acc = TFile::Open(TString::Format("AcceptanceCombined/Nominal/CombAcceptance%s_%s_ResponsesEfficiency_TTToHadronic.root", 
                                varParton.Data(), variable[ivar].Data()));
        TFile *f_eff = TFile::Open(TString::Format("EfficiencyCombined/Nominal/CombEfficiency%s_%s_ResponsesEfficiency_TTToHadronic.root", 
                                varParton.Data(), variable[ivar].Data()));
        TH1F *eff_combined = (TH1F*)f_eff->Get("efficiency");
        TH1F *acc_combined = (TH1F*)f_acc->Get("acceptance");

        // now read the Top-18-013 analysis plots
        // ../ResponseMatrices/PartonEfficiencyAll_July19.root
        TString oldVarParton = "Parton";
        if (varParton.EqualTo("Particle")) oldVarParton = "Gen";
        TFile *old_eff_acc = TFile::Open(TString::Format("../ResponseMatrices/%sEfficiencyAll_July19.root", varParton.Data()));
        TEfficiency *eff = (TEfficiency*)old_eff_acc->Get(TString::Format("Eff_%s_%s_Nominal", variable[ivar].Data(), oldVarParton.Data()));
        TEfficiency *acc = (TEfficiency*)old_eff_acc->Get(TString::Format("Eff_%s_%s_Nominal_common", variable[ivar].Data(), oldVarParton.Data()));


        eff_combined->SetMarkerStyle(28);
        eff_combined->SetMarkerSize(1.2);
        eff_combined->SetMarkerColor(kMagenta);
        eff_combined->SetLineColor(kMagenta);
        eff_combined->SetTitle(TString::Format("Efficiency for %s", variable[ivar].Data()));
        eff_combined->SetName(TString::Format("Efficiency for %s", variable[ivar].Data()));

        eff->SetTitle(TString::Format("Efficiency; %s %s; Efficiency", variable[ivar].Data(), expr[ivar].Data()));
        acc->SetTitle(TString::Format("Acceptance; %s %s; Acceptance", variable[ivar].Data(), expr[ivar].Data()));


        acc_combined->SetMarkerStyle(28);
        acc_combined->SetMarkerSize(1.2);
        acc_combined->SetMarkerColor(kMagenta);
        acc_combined->SetLineColor(kMagenta);
        acc_combined->SetTitle(TString::Format("Acceptance for %s", variable[ivar].Data()));
        acc_combined->SetName(TString::Format("Acceptance for %s", variable[ivar].Data()));


        leg->AddEntry(eff_combined, "Combined", "lep");
        leg->AddEntry(eff, "TOP-18-013", "lep");
        leg_acc->AddEntry(acc_combined, "Combined", "lep");
        leg_acc->AddEntry(eff, "TOP-18-013", "lep");

        lumi_13TeV = TString::Format("%0.1f fb^{-1}", luminosity["luminosityAll"]/1000);
        int iPeriod = 13;
        int iPos = 0;
        //extraTextFactor = 0.14;
        writeExtraText=true;
        
        // plot acceptance 

        TCanvas *can_acc = new TCanvas(TString::Format("canAcc_%s",variable[ivar].Data()),
                                        TString::Format("canAcc_%s",variable[ivar].Data()) , 800,600);
        can_acc->cd();
        leg->Draw();
        acc->Draw();
        acc_combined->Draw("same");
        

        gPad->Update(); 
        auto graph = acc->GetPaintedGraph(); 
        graph->SetMinimum(0.4);
        graph->SetMaximum(1.); 
        gPad->Update(); 

        // CMS_lumi(can_acc, "combined", iPos);
        leg_acc->Draw();    

        can_acc->SaveAs(TString::Format("Comparison_EffAccResponses_TOP18013/Acceptance%s_%s.pdf", 
                        varParton.Data(), variable[ivar].Data()), "pdf");
        

        // plot efficiency 
        TCanvas *can_eff = new TCanvas(TString::Format("canEff_%s",variable[ivar].Data()),
                                        TString::Format("canEff_%s",variable[ivar].Data()) , 800,600);
        
        can_eff->cd();
        leg->Draw();
        eff->Draw();
        eff_combined->Draw("same");
        
        gPad->Update(); 

        auto graph_eff = eff->GetPaintedGraph(); 
        graph_eff->SetMinimum(0);
        if (isParton)
            graph_eff->SetMaximum(0.1);
        else 
            graph_eff->SetMaximum(0.3);
        gPad->Update(); 
        
        // CMS_lumi(can_eff, "combined", iPos);
        leg->Draw();
        
        can_eff->SaveAs(TString::Format("Comparison_EffAccResponses_TOP18013/Efficiency%s_%s.pdf", 
                        varParton.Data(), variable[ivar].Data()), "pdf");

        
        // CMS_lumi(can_eff_comb, "combined", iPos);
        
        //break;
    } // end of variables loop 

}