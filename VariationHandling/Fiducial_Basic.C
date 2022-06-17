#include "BASE.h"

#include "FinalResultsConstants.h"

#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "../CMS_plots/CMS_lumi.C"
#include "../CMS_plots/tdrstyle.C"
#include "TemplateConstants.h"

void DrawWithRatio(TCanvas *can, std::vector<TH1F *> histograms, int index, bool isNormalized, TString partonParticle)
{
    std::vector<Color_t> colors = {kRed, kBlue};
    std::vector<bool> putInLegend = {true, true};
    std::vector<TString> legendTitles = {"MC", "Data"};
    std::vector<TString> legendDrawOption = {"l", "l"};
    std::vector<TString> drawOptions = {"", "SAME EP"};
    std::vector<bool> drawRatio = {true, true};



    TPad *upperPad = new TPad("upperPad", "upperPad", 0, 0.3, 1, 1.0);
    upperPad->SetBottomMargin(0.05);
    upperPad->Draw();
    upperPad->cd();

    //histogramsToDraw.push_back(finalTheory);
    //histogramsToDraw.push_back(finalResult);
    
    if (AnalysisConstants::axisInLogScale[index])
    {
        upperPad->SetLogy();
    }
    can->cd();
    TPad *lowerPad = new TPad("lowerPad", "lowerPad", 0, 0.05, 1, 0.3);
    lowerPad->SetBottomMargin(0.22);
    lowerPad->SetTopMargin(0.03);
    lowerPad->Draw();
    lowerPad->SetGridy();
    lowerPad->cd();

    TLegend *leg1 = new TLegend(0.69, 0.7, 0.9, 0.9);
    TH1F *dataHist = (TH1F *)histograms[1]->Clone("DataHist");

    for (unsigned int i = 0; i < histograms.size(); i++)
    {
        
        can->cd();
        upperPad->cd();
        if (AnalysisConstants::axisInLogScale[index])
        {
            upperPad->SetLogy();
        }
        TH1F *hist = histograms[i];
        hist->SetTitle("");
        hist->SetLineColor(colors[i]);
        if (!isNormalized)
            hist->GetYaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::fiducialYAxisValues[index][0],
                                        AnalysisConstants::FinalResultsConstans::fiducialYAxisValues[index][1]);
        else 
            hist->GetYaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::partonYAxisValuesNormalized[index][0],
                                        AnalysisConstants::FinalResultsConstans::partonYAxisValuesNormalized[index][1]);

        TH1F *ratio = (TH1F *)hist->Clone("ratio");
        ratio->Divide(dataHist);
        if (partonParticle.EqualTo("Parton")){
        hist->GetXaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::partonXAxisValues[index][0],
                                    AnalysisConstants::FinalResultsConstans::partonXAxisValues[index][1]);
        }
        else{
            hist->GetXaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::particleXAxisValues[index][0],
                                    AnalysisConstants::FinalResultsConstans::particleXAxisValues[index][1]);
        }
        hist->Draw(drawOptions[i]);

        can->cd();
        lowerPad->cd();

        if (drawRatio[i])
        {
        ratio->GetYaxis()->SetRangeUser(0, 2);
        ratio->GetYaxis()->SetNdivisions(505);
        ratio->GetYaxis()->SetLabelSize(0.1);
        ratio->GetYaxis()->SetTitle("#frac{Data}{MC}");
        ratio->GetYaxis()->SetTitleSize(0.115);
        ratio->GetYaxis()->SetTitleOffset(0.35);
        ratio->GetYaxis()->CenterTitle(true);
        ratio->GetXaxis()->SetTitleSize(0.115);
        ratio->GetXaxis()->SetTitleOffset(1.);
        ratio->GetXaxis()->SetLabelSize(0.12);
        ratio->GetXaxis()->SetLabelOffset(0.015);
        if (partonParticle.EqualTo("Parton")){
            ratio->GetXaxis()->SetTitle(AnalysisConstants::fiducialAxisTitles[index]);
            ratio->GetXaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::partonXAxisValues[index][0],
                                            AnalysisConstants::FinalResultsConstans::partonXAxisValues[index][1]);
        }

        ratio->Draw(drawOptions[i]);
        }

        if (putInLegend[i])
        {
        leg1->AddEntry(hist, legendTitles[i], legendDrawOption[i]);
        }
    }

    can->cd();
    upperPad->cd();
    leg1->Draw();

    TList *primitives = upperPad->GetListOfPrimitives();

    for (int i = 0; i < primitives->GetSize(); i++)
    {
        std::cout << primitives->At(i)->GetName() << std::endl;
    }

    lumi_13TeV = TString::Format("%0.1f fb^{-1}", luminosity["luminosityAll"]/1000);
    int iPeriod = 13;
    int iPos = 0;
    extraTextFactor = 0.14;
    writeExtraText=true;

    CMS_lumi(upperPad, "combined", iPos);
    }
    

void Fiducial_Basic(bool normalized = true)
{
    initFilesMapping();
    gStyle->SetOptStat(0);
    TString partonParticle = "Parton";
    TString baseDir = "/Users/georgebakas/Documents/HEP-NTUA_ul/VariationHandling";
    AnalysisConstants::initConstants();
    TFile *nominalFile = TFile::Open("testFile_.root");
    TFile *nominalFileTheory = TFile::Open(TString::Format("%s/UnfoldedCombined/Nominal/OutputFileParton.root",
                                                                baseDir.Data()));
    TString outputDir = TString::Format("%s/FiducialMeasurementCombined/results",
                                        baseDir.Data());
    CheckAndCreateDirectory("results/FiducialMeasurement");

    TString theoryYear = "2018";

    TFile *nominalAmcAtNloFile = TFile::Open(TString::Format(
                            "../VariationHandling_Theory_amc@NLO/%s/Nominal/Histograms_TTJets.root", theoryYear.Data()));

    for (unsigned int v = 0; v < AnalysisConstants::variables.size(); v++)
    {   
        TString variable = AnalysisConstants::variables[v];

        TH1F *nominalHistogram = (TH1F *)nominalFile->Get(TString::Format("combined_%s",
                                                                        variable.Data()));

        nominalHistogram->Scale(1. / AnalysisConstants::luminositiesSR["comb"], "width");
        float_t nominalHistogramYield = nominalHistogram->Integral();
        if (normalized) nominalHistogram->Scale(1 / nominalHistogramYield);

        TH1F *finalResult = (TH1F *)nominalHistogram->Clone(TString::Format("FinalResult%s_%s",
                                                                        (normalized ? "Norm" : "Final"), 
                                                                        variable.Data()));
                                                                

        TH1F *theoryHistogram = (TH1F *)nominalFileTheory->Get(TString::Format("hFidTheory%s_%s",
                                                                        (normalized ? "Norm" : ""),
                                                                        variable.Data()));

        TH1F *finalTheory = (TH1F *)theoryHistogram->Clone(TString::Format("FinalTheory%s_%s",
                                                                        (normalized ? "Norm" : ""),
                                                                        variable.Data()));
        

        TCanvas *c1 = new TCanvas(TString::Format("FinalResult_%s",
                                                variable.Data()),
                                TString::Format("FinalResult_%s",
                                                variable.Data()),
                                600, 600);

        TH1F *theoryHistogramValue = (TH1F *)theoryHistogram->Clone(TString::Format("theoryHistogramValue%s",
                                                                                    (normalized ? "_normalized" : "")));
        std::vector<TH1F *> histogramsToDraw;
        
        histogramsToDraw.push_back(finalTheory);
        histogramsToDraw.push_back(finalResult);

        DrawWithRatio(c1, histogramsToDraw, v, normalized, partonParticle);

        TString draw_name;
        if (partonParticle.EqualTo("Particle")) draw_name = AnalysisConstants::particleVariables[v].Data();
        else draw_name = AnalysisConstants::partonVariables[v].Data();
        
        c1->SaveAs(TString::Format("%s/FinalResult_%s%s.pdf",
                                outputDir.Data(),
                                draw_name.Data(),
                                (normalized ? "_normalized" : "")),
                "pdf");
    }
}
