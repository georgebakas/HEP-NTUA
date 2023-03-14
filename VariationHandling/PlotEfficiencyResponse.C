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


void PlotEfficiencyResponse(bool isParton = false)
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
        int iPeriod = 13;
        int iPos = 0;
        //extraTextFactor = 0.14;
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

        CMS_lumi(can_acc, "combined", iPos);
        leg_acc->Draw();    

        can_acc->SaveAs(TString::Format("Comparison_EffAccResponses/combined/Acceptance%s_%s.pdf", 
                        varParton.Data(), variable[ivar].Data()), "pdf");
        

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
        
        CMS_lumi(can_eff, "combined", iPos);
        leg->Draw();
        
        can_eff->SaveAs(TString::Format("Comparison_EffAccResponses/combined/Efficiency%s_%s.pdf", 
                        varParton.Data(), variable[ivar].Data()), "pdf");

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
        
        // plot the response matrices (these are scaled to their integrals!!)
        // x-axis is particle or parton while y is reco
        TH2F *hCombined = (TH2F*)hResponse[0]->Clone(TString::Format(
                                    "CombinedResponseMatrix_%s", variable[ivar].Data()));
        TString tempVar;
        if(isParton)
            tempVar = variableParton[ivar];
        else
            tempVar = variableGen[ivar];
        for (int iy=0; iy<4; iy++)
        {   
            if (iy!=0) hCombined->Add(hResponse[iy]);
            TCanvas *can_response_yearly = new TCanvas(
                            TString::Format("canRespYearly%s_%s", variable[ivar].Data(),years[iy].Data()),
                            TString::Format("canRespYearly%s_%s", variable[ivar].Data(),years[iy].Data()), 800, 600);
            can_response_yearly->cd();
            //hResponse[iy]->Scale(1/hResponse[iy]->Integral());
            gStyle->SetPaintTextFormat("4.1f");
            hResponse[iy]->GetXaxis()->SetTitle(tempVar+" "+expr[ivar]);
            hResponse[iy]->GetYaxis()->SetTitle(variable[ivar]+" "+expr[ivar]);
            hResponse[iy]->SetTitle("");
            hResponse[iy]->Draw("text colz0");
            CMS_lumi(can_response_yearly, "combined", iPos);
            can_response_yearly->SaveAs(TString::Format("Comparison_EffAccResponses/%s/Response%s_%s.pdf", 
                        years[iy].Data(), varParton.Data(), variable[ivar].Data()), "pdf");
        }
        
        TCanvas *can_response_comb = new TCanvas(
                            TString::Format("canRespComb%s", variable[ivar].Data()),
                            TString::Format("canRespComb%s", variable[ivar].Data()), 800, 600);
        can_response_comb->cd();
        gStyle->SetPaintTextFormat("4.1f");
        hCombined->GetXaxis()->SetTitle(tempVar+" "+expr[ivar]);
        hCombined->GetYaxis()->SetTitle(variable[ivar]+" "+expr[ivar]);
        //hCombined->Scale(1/hCombined->Integral());
        hCombined->Draw("textcolz0");
        CMS_lumi(can_response_comb, "combined", iPos);
        can_response_comb->SaveAs(TString::Format("Comparison_EffAccResponses/combined/Response%s_%s.pdf", 
                        varParton.Data(), variable[ivar].Data()), "pdf");

        int sizeBins = hCombined->GetXaxis()->GetNbins();
        purityParton[ivar] 	  = new TH1D(TString::Format("PurityParton_%s", variable[ivar].Data()),
                                            TString::Format("PurityParton_%s", variable[ivar].Data()), sizeBins, 0, sizeBins+1);
  	    stabilityParton[ivar] = new TH1D(TString::Format("StabilityParton_%s", variable[ivar].Data()),
                                            TString::Format("StabilityParton_%s", variable[ivar].Data()), sizeBins, 0, sizeBins+1);

        //get the stability and purity for the combined:
        float sumOfRowsParton[sizeBins], sumOfColsParton[sizeBins];
        float sumOfRowsParticle[sizeBins], sumOfColsParticle[sizeBins];

        for(int i=1; i<=sizeBins; i++)
        {
            sumOfColsParton[i] = ((TH1D*)hCombined->ProjectionX())->GetBinContent(i);
            sumOfRowsParton[i] = ((TH1D*)hCombined->ProjectionY())->GetBinContent(i);

            for(int j=1; j<=sizeBins; j++)
            {
                if(i==j)
                {
                    float initContentParton = hCombined->GetBinContent(i,j);
                    purityParton[ivar]->SetBinContent(i,initContentParton/sumOfColsParton[i]);
                    stabilityParton[ivar]->SetBinContent(i,initContentParton/sumOfRowsParton[i]);
                }
            }
        }
        purityParton[ivar]->SetLineColor(kBlack);
        stabilityParton[ivar]->SetLineColor(kBlue);

        // plot purity and stability 

        TCanvas *can_purStab = new TCanvas(TString::Format("can_purStab%s",variable[ivar].Data()),
                                        TString::Format("can_purStab%s",variable[ivar].Data()) , 800,600);
        TLegend *leg_purstab = new TLegend(0.3, 0.2, 0.5, 0.4);
        can_purStab->cd();
        purityParton[ivar]->GetYaxis()->SetRangeUser(0,1);
        purityParton[ivar]->GetXaxis()->SetTitle("Bin Number");
        purityParton[ivar]->Draw();
        purityParton[ivar]->SetName("");
        purityParton[ivar]->SetTitle("");
        stabilityParton[ivar]->Draw("same");

        leg_purstab->SetHeader(variable[ivar]);
        leg_purstab->AddEntry(purityParton[ivar], "Purity" ,"lp");
        leg_purstab->AddEntry(stabilityParton[ivar], "Stability", "lp");
        
        CMS_lumi(can_purStab, "combined", iPos);
        leg_purstab->Draw();
        
        can_purStab->SaveAs(TString::Format("Comparison_EffAccResponses/combined/PurityStability%s_%s.pdf", 
                        varParton.Data(), variable[ivar].Data()), "pdf");
        
        //plot efficiency and acceptance per year and combined in the same plot 
        for (int iy=0; iy<4; iy++)
        {   
            TCanvas *can_eff_yearly = new TCanvas(
                            TString::Format("canEffYearly%s_%s", variable[ivar].Data(),years[iy].Data()),
                            TString::Format("canEffYearly%s_%s", variable[ivar].Data(),years[iy].Data()), 800, 600);
            can_eff_yearly->cd();
            TEfficiency *effCombYearly = (TEfficiency*)efficiency[iy]->Clone("effCombYearly");
            TEfficiency *accCombYearly = (TEfficiency*)acceptance[iy]->Clone("accCombYearly");
            accCombYearly->SetMarkerStyle(20);
            accCombYearly->SetMarkerSize(0.7);
            accCombYearly->SetMarkerColor(kBlue);
            accCombYearly->SetLineColor(kBlue);
            accCombYearly->SetFillColorAlpha(kBlue, 0.4);
            accCombYearly->SetFillStyle(1001);
            
            effCombYearly->SetMarkerStyle(20);
            effCombYearly->SetMarkerSize(0.7);
            effCombYearly->SetMarkerColor(kRed);
            effCombYearly->SetLineColor(kRed);
            effCombYearly->SetFillColorAlpha(kRed, 0.4); 
            effCombYearly->SetFillStyle(1001);
            accCombYearly->SetTitle("");
            accCombYearly->Draw();
            effCombYearly->Draw("same");     
            TLegend *legEffAcc = new TLegend(0.5, 0.5, 0.6, 0.6);
            legEffAcc->AddEntry(effCombYearly, "f1", "l");
            legEffAcc->AddEntry(accCombYearly, "f2", "l");
            legEffAcc->Draw();
            gPad->Update(); 
            auto graph_eff_yearly = accCombYearly->GetPaintedGraph(); 
            graph_eff_yearly->SetMinimum(0);
            graph_eff_yearly->SetMaximum(1);
            graph_eff_yearly->GetYaxis()->SetTitle("Fractions");
            CMS_lumi(can_eff_yearly, "combined", iPos);
            can_eff_yearly->SaveAs(TString::Format("Comparison_EffAccResponses/%s/EfficiencyAcceptance%s_%s.pdf", 
                        years[iy].Data(), varParton.Data(), variable[ivar].Data()), "pdf");
        }

        TCanvas *can_eff_comb = new TCanvas(
                        TString::Format("canEffYearly%s", variable[ivar].Data()),
                        TString::Format("canEffYearly%s", variable[ivar].Data()), 800, 600);
        
        TLegend *legEffAcc = new TLegend(0.5, 0.5, 0.6, 0.6);
        
        can_eff_comb->cd();
        TH1F *effCombClone = (TH1F*)eff_combined->Clone("effCombClone");
        TH1F *accCombClone = (TH1F*)acc_combined->Clone("accCombClone");
        effCombClone->SetTitle("");
        effCombClone->SetMarkerColor(kRed);
        effCombClone->SetLineColor(kRed);
        accCombClone->SetMarkerColor(kBlue);
        accCombClone->SetLineColor(kBlue);
        effCombClone->SetFillStyle(1001);
        accCombClone->SetFillStyle(1001);
        effCombClone->SetFillColorAlpha(kRed, 0.4);
        accCombClone->SetFillColorAlpha(kBlue, 0.4);
        effCombClone->GetXaxis()->SetTitle(variable[ivar]+" "+expr[ivar]);
        accCombClone->GetXaxis()->SetTitle(variable[ivar]+" "+expr[ivar]);
        accCombClone->SetTitle("");
        effCombClone->Draw("E2");
        accCombClone->Draw("E2 same");
        effCombClone->SetMinimum(0);
        effCombClone->SetMaximum(1);
        legEffAcc->AddEntry(effCombClone, "f1", "l");
        legEffAcc->AddEntry(accCombClone, "f2", "l");
        legEffAcc->SetBorderSize(0);
        legEffAcc->Draw();

        CMS_lumi(can_eff_comb, "combined", iPos);
        can_eff_comb->SaveAs(TString::Format("Comparison_EffAccResponses/combined/EfficiencyAcceptance%s_%s.pdf", 
                        varParton.Data(), variable[ivar].Data()), "pdf");
        
        //break;
    } // end of variables loop 

    myfile.close();
}