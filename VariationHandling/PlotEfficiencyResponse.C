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
#include "TUnfold.h"
#include "TUnfoldDensity.h"
#include "../CMS_plots/tdrstyle.C"
#include "../CMS_plots/CMS_lumi.C"

using namespace std;

#include "UseTUnfoldDensity.cpp"

using std::cin;
using std::cout;
using std::endl;

#include "TemplateConstants.h"


void PlotEfficiencyResponse(bool isParton = true)
{   
    /* init variables and values */
    initFilesMapping();
    setTDRStyle();
    gStyle->SetOptStat(0);
    const int NVAR = 10;
    TString variable[NVAR]   = {"mJJ", "ptJJ", "yJJ","jetPt0","jetPt1", "jetY0", "jetY1", "chi", "cosTheta_0", "cosTheta_1"};
    TString variableParton[NVAR] = {"mTTbarParton", "ptTTbarParton", "yTTbarParton","partonPt0", "partonPt1", "partonY0", "partonY1","chiParton", "cosThetaParton_0", "cosThetaParton_1"};
    TString variableGen[NVAR] = {"mJJGen", "ptJJGen", "yJJGen","genjetPt0", "genjetPt1", "genjetY0", "genjetY1", "chiParticle", "cosThetaParticle_0", "cosThetaParticle_1"};

    std::vector<int> color = {kRed, kGreen, kBlue, kBlack};
    TString varParton = "Parton";
    if(!isParton) varParton = "Particle";
    /* Loop on all variables and create canvases */

    TString years[4] = {"2016_preVFP", "2016_postVFP", "2017", "2018"};
    
    TString variations[5] = {"bTagVariation", "JES", "PDFWeights", "PSWeights", "ScaleWeights"};

    //write to a txt file the Kolmogorov tests
    ofstream myfile;
    myfile.open ("Comparison_EffAccResponses/KS_results.txt");

    for (int ivar=0; ivar<NVAR; ivar++)
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
        
        eff_combined->SetMarkerStyle(28);
        eff_combined->SetMarkerSize(1.2);
        eff_combined->SetMarkerColor(kMagenta);
        eff_combined->SetLineColor(kMagenta);

        acc_combined->SetMarkerStyle(28);
        acc_combined->SetMarkerSize(1.2);
        acc_combined->SetMarkerColor(kMagenta);
        acc_combined->SetLineColor(kMagenta);

        // get the error propagation


        //get response matrix
        for(int iy = 0; iy<4; iy++)
        {   
            TString tempVar;
            if(isParton)
                tempVar = variableParton[ivar];
            else
                tempVar = variableGen[ivar];

            TFile *inf_had = TFile::Open(TString::Format("%s/ResponsesNominal/ResponsesEfficiency_TTToHadronic.root", years[iy].Data()));
            TFile *inf_sem = TFile::Open(TString::Format("%s/ResponsesNominal/ResponsesEfficiency_TTToSemiLeptonic.root", years[iy].Data()));
            TFile *inf_dil = TFile::Open(TString::Format("%s/ResponsesNominal/ResponsesEfficiency_TTTo2L2Nu.root", years[iy].Data()));

            
            cout << years[iy]<< "\n";
            // get responses 
            hResponse_had[iy] = (TH2F*)inf_had->Get(TString::Format("h%sResponse_%s",varParton.Data(), variable[ivar].Data()));
            hResponse_sem[iy] = (TH2F*)inf_sem->Get(TString::Format("h%sResponse_%s",varParton.Data(), variable[ivar].Data()));
            hResponse_dil[iy] = (TH2F*)inf_dil->Get(TString::Format("h%sResponse_%s",varParton.Data(), variable[ivar].Data()));

            //get efficiencies and acceptances 
            TEfficiency *eff_had = (TEfficiency*)inf_had->Get(TString::Format("Efficiency%s_%s",varParton.Data(), tempVar.Data()));
            TEfficiency *eff_sem = (TEfficiency*)inf_sem->Get(TString::Format("Efficiency%s_%s",varParton.Data(), tempVar.Data()));
            TEfficiency *eff_dil = (TEfficiency*)inf_dil->Get(TString::Format("Efficiency%s_%s",varParton.Data(), tempVar.Data()));

            //get efficiencies and acceptances 
            TEfficiency *acc_had = (TEfficiency*)inf_had->Get(TString::Format("Acceptance%s_%s",varParton.Data(), variable[ivar].Data()));
            TEfficiency *acc_sem = (TEfficiency*)inf_sem->Get(TString::Format("Acceptance%s_%s",varParton.Data(), variable[ivar].Data()));
            TEfficiency *acc_dil = (TEfficiency*)inf_dil->Get(TString::Format("Acceptance%s_%s",varParton.Data(), variable[ivar].Data()));
            
            // add all (regardless of ttbar process) to the hResponse per variable
            // responses are already scaled to XSEC and lumi
            hResponse[iy] = (TH2F*)hResponse_had[iy]->Clone();
            hResponse[iy]->Add(hResponse_sem[iy]);
            hResponse[iy]->Add(hResponse_dil[iy]);

            efficiency[iy] = (TEfficiency*)eff_had->Clone();
            *efficiency[iy] += (*eff_sem);
            *efficiency[iy] += (*eff_dil); 

            acceptance[iy] = (TEfficiency*)acc_had->Clone();
            *acceptance[iy] += (*acc_sem);
            *acceptance[iy] += (*acc_dil); 
            
            efficiency[iy]->SetMarkerStyle(20+iy);
            efficiency[iy]->SetMarkerSize(1.2);
            efficiency[iy]->SetMarkerColor(color[iy]);
            efficiency[iy]->SetLineColor(color[iy]);

            acceptance[iy]->SetMarkerStyle(20+iy);
            acceptance[iy]->SetMarkerSize(1.2);
            acceptance[iy]->SetMarkerColor(color[iy]);
            acceptance[iy]->SetLineColor(color[iy]);
            
            leg->AddEntry(efficiency[iy], years[iy], "lep");
            leg_acc->AddEntry(efficiency[iy], years[iy], "lep");
        } // end of years loop

        leg->AddEntry(eff_combined, "Combined", "lep");
        leg_acc->AddEntry(acc_combined, "Combined", "lep");

        lumi_13TeV = TString::Format("%0.1f fb^{-1}", luminosity["luminosityAll"]/1000);
        int iPeriod = 4;
        int iPos = 1;
        writeExtraText=true;
        
        // plot acceptance 

        TCanvas *can_acc = new TCanvas(TString::Format("canAcc_%s",variable[ivar].Data()),
                                        TString::Format("canAcc_%s",variable[ivar].Data()) , 800,600);
        can_acc->cd();
        leg->Draw();
        acceptance[0]->Draw();
        acceptance[1]->Draw("same");
        acceptance[2]->Draw("same");
        acceptance[3]->Draw("same");
        acc_combined->Draw("same");

        gPad->Update(); 
        auto graph = acceptance[0]->GetPaintedGraph(); 
        graph->SetMinimum(0.4);
        graph->SetMaximum(1.); 
        gPad->Update(); 

        CMS_lumi(can_acc, iPeriod, iPos);
        leg_acc->Draw();    
        can_acc->SaveAs(TString::Format("Comparison_EffAccResponses/Acceptance%s_%s.png", 
                        varParton.Data(), variable[ivar].Data()));
        

        // plot efficiency 
        TCanvas *can_eff = new TCanvas(TString::Format("canEff_%s",variable[ivar].Data()),
                                        TString::Format("canEff_%s",variable[ivar].Data()) , 800,600);
        
        can_eff->cd();
        leg->Draw();
        efficiency[0]->Draw();
        efficiency[1]->Draw("same");
        efficiency[2]->Draw("same");
        efficiency[3]->Draw("same");
        eff_combined->Draw("same");

        gPad->Update(); 
        auto graph_eff = efficiency[0]->GetPaintedGraph(); 
        graph_eff->SetMinimum(0);
        if (isParton)
            graph_eff->SetMaximum(0.1);
        else 
            graph_eff->SetMaximum(0.3);
        gPad->Update(); 
        
        CMS_lumi(can_eff, iPeriod, iPos);
        leg->Draw();

        
        can_eff->SaveAs(TString::Format("Comparison_EffAccResponses/Efficiency%s_%s.png", 
                        varParton.Data(), variable[ivar].Data()));

        // now deal with the responses
        // I will get the results from the Kolmogorov test 
        // The returned function value is the probability of test (much less than one means NOT compatible)
        
        myfile<<"------------"<< "\n";
        myfile<< variable[ivar]<< "\n";
        myfile<<"16pre-16post: " <<hResponse[0]->KolmogorovTest(hResponse[1])<< "\n";
        myfile<<"16pre-17: " <<hResponse[0]->KolmogorovTest(hResponse[2])<< "\n";
        myfile<<"16pre-18: " <<hResponse[0]->KolmogorovTest(hResponse[3])<< "\n";

        myfile<<"16post-16pre: " <<hResponse[1]->KolmogorovTest(hResponse[0])<< "\n";
        myfile<<"16post-17: " <<hResponse[1]->KolmogorovTest(hResponse[2])<< "\n";
        myfile<<"16post-18: " <<hResponse[1]->KolmogorovTest(hResponse[3])<< "\n";

        myfile<<"17-16pre: " <<hResponse[2]->KolmogorovTest(hResponse[0])<< "\n";
        myfile<<"17-16post: " <<hResponse[2]->KolmogorovTest(hResponse[1])<< "\n";
        myfile<<"17-18: " <<hResponse[2]->KolmogorovTest(hResponse[3])<< "\n";

        myfile<<"18-16pre: " <<hResponse[3]->KolmogorovTest(hResponse[0])<< "\n";
        myfile<<"18-16post: " <<hResponse[3]->KolmogorovTest(hResponse[1])<< "\n";
        myfile<<"18-17: " <<hResponse[3]->KolmogorovTest(hResponse[2])<< "\n";
        
        myfile<<"------------"<< "\n";
        
        
    } // end of variables loop 

    myfile.close();
}