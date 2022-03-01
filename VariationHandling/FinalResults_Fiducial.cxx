#include "BASE.h"

#include "FinalResultsConstants.h"

#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "../CMS_plots/CMS_lumi.C"
#include "../CMS_plots/tdrstyle.C"

void DrawWithRatio(TCanvas *can, std::vector<TH1F *> histograms, int index, bool isNormalized, TString partonParticle)
{
    std::vector<Color_t> colors = {kBlack, kBlack, kRed, kRed, kBlue, kBlue};
    std::vector<Color_t> fillColors = {kGray, kGray, kRed, kRed, kBlue, kBlue};
    std::vector<Int_t> markerStyle = {20, 1, 1, 1, 1, 1};
    std::vector<int> fillStyles = {1001, 1001, 3675, 3675, 3257, 3257};
    std::vector<TString> drawOptions = {"E2", "SAME EP", "SAME E2", "SAME", "SAME E2", "SAME"};
    std::vector<bool> drawRatio = {true, true, true, true, true, true};
    std::vector<bool> putInLegend = {true, true, true, true, true, true};
    std::vector<TString> legendTitles = {"Data", "Total unc.", "amc@NLO+Pythia8", "amc@NLO+Pythia8 unc","Powheg+Pythia8", "Powheg+Pythia8 unc"};
    std::vector<TString> legendDrawOption = {"EP", "f", "EL", "f", "EL", "f"};

    TPad *upperPad = new TPad("upperPad", "upperPad", 0, 0.3, 1, 1.0);
    upperPad->SetBottomMargin(0.05);
    upperPad->Draw();
    upperPad->cd();
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

  TLegend *leg1 = new TLegend(AnalysisConstants::FinalResultsConstans::legendPositions[index][0],
                                AnalysisConstants::FinalResultsConstans::legendPositions[index][1],
                                AnalysisConstants::FinalResultsConstans::legendPositions[index][2],
                                AnalysisConstants::FinalResultsConstans::legendPositions[index][3]);
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
        //cout<<hist->GetTitle()<<endl;
        hist->SetTitle("");
        hist->SetMarkerColor(colors[i]);
        hist->SetLineColor(colors[i]);
        hist->SetFillColor(fillColors[i]);
        hist->SetMarkerStyle(markerStyle[i]);
        hist->SetLineWidth(2);
        hist->SetFillStyle(fillStyles[i]);
        if (partonParticle.EqualTo("Parton")){
            if (!isNormalized)
                hist->GetYaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::partonYAxisValues[index][0],
                                            AnalysisConstants::FinalResultsConstans::partonYAxisValues[index][1]);
            else 
                hist->GetYaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::partonYAxisValuesNormalized[index][0],
                                            AnalysisConstants::FinalResultsConstans::partonYAxisValuesNormalized[index][1]);
        }
        else{
            if (!isNormalized)
                hist->GetYaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::particleYAxisValues[index][0],
                                            AnalysisConstants::FinalResultsConstans::particleYAxisValues[index][1]);
            else 
                hist->GetYaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::particleYAxisValuesNormalized[index][0],
                                            AnalysisConstants::FinalResultsConstans::particleYAxisValuesNormalized[index][1]);
        }

        TH1F *ratio = (TH1F *)hist->Clone("ratio");
        ratio->Divide(dataHist);
        for (int bin = 1; bin <= hist->GetNbinsX(); bin++)
        {
        ratio->SetBinContent(bin, ratio->GetBinContent(bin) - 1);
        }
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
        ratio->GetYaxis()->SetRangeUser(-2, 2);
        ratio->GetYaxis()->SetNdivisions(505);
        ratio->GetYaxis()->SetLabelSize(0.1);
        ratio->GetYaxis()->SetTitle("The./data - 1");
        ratio->GetYaxis()->SetTitleSize(0.115);
        ratio->GetYaxis()->SetTitleOffset(0.35);
        ratio->GetYaxis()->CenterTitle(true);
        ratio->GetXaxis()->SetTitleSize(0.115);
        ratio->GetXaxis()->SetTitleOffset(1.);
        ratio->GetXaxis()->SetLabelSize(0.12);
        ratio->GetXaxis()->SetLabelOffset(0.015);
        if (partonParticle.EqualTo("Parton")){
            ratio->GetXaxis()->SetTitle(AnalysisConstants::partonAxisTitles[index]);
            ratio->GetXaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::partonXAxisValues[index][0],
                                            AnalysisConstants::FinalResultsConstans::partonXAxisValues[index][1]);
        }
        else{
            ratio->GetXaxis()->SetTitle(AnalysisConstants::particleAxisTitles[index]);
            ratio->GetXaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::particleXAxisValues[index][0],
                                            AnalysisConstants::FinalResultsConstans::particleXAxisValues[index][1]);
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
    leg1->SetBorderSize(0);
    leg1->Draw();

    TList *primitives = upperPad->GetListOfPrimitives();

    for (int i = 0; i < primitives->GetSize(); i++)
    {
        std::cout << primitives->At(i)->GetName() << std::endl;
    }

    //float extraTextFactor = 0.14;
    int iPeriod = 13;
    int iPos = 0;
    writeExtraText=true;
    CMS_lumi(upperPad, iPeriod, iPos);
    }

    void AddSystematicToErrorBar(TH1F *nominal, TH1F *systematic)
    {
    for (int bin = 1; bin <= nominal->GetNbinsX(); bin++)
    {
        float nominalValue = nominal->GetBinContent(bin);
        float variationValue = systematic->GetBinContent(bin);
        float difference = TMath::Abs(nominalValue - variationValue);
        /*if ((difference / nominalValue) > 0.3)
        {
            cout<<systematic->GetName()<<endl;
            break;
        } */
        float error = TMath::Sqrt(TMath::Power(nominal->GetBinError(bin), 2) +
                                TMath::Power(difference, 2));
        nominal->SetBinError(bin, error);
    }
    }

void FinalResults_Fiducial(bool normalized = true)
{
    //setTDRStyle();
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

        TH1F *theoryAmcAtNloHistogram = (TH1F *)nominalAmcAtNloFile->Get(TString::Format("hReco_%s", 
                                        variable.Data()));
        
        theoryAmcAtNloHistogram->Scale(1. / AnalysisConstants::luminositiesSR["2018"], "width");
        float_t theoryAmcAtNloYield = theoryAmcAtNloHistogram->Integral();
        if (normalized)
            theoryAmcAtNloHistogram->Scale(1 / theoryAmcAtNloYield);

        TH1F *finalTheoryAmcAtNlo = (TH1F *)theoryAmcAtNloHistogram->Clone(TString::Format("FinalTheoryAmcAtNlo%s_%s",
                                                                            (normalized ? "Norm" : ""),
                                                                            variable.Data()));

        TCanvas *c1 = new TCanvas(TString::Format("FinalResult_%s",
                                                variable.Data()),
                                TString::Format("FinalResult_%s",
                                                variable.Data()),
                                600, 600);
        

        for (unsigned int var = 0; var < AnalysisConstants::variationsFiducial.size(); var++)
        {
            TString variation = AnalysisConstants::variationsFiducial[var];
            //cout<<var<<" "<<variation<<endl;
            TString tempVariation = "";
            if (variation.Contains("isr") || variation.Contains("fsr"))
                tempVariation = "PSWeights";
            else if (variation.Contains("up") || variation.Contains("down"))
                tempVariation = "bTagVariation";
            else if (variation.Contains("pdf"))
                tempVariation = "PDFWeights";
            else if (variation.Contains("scale"))
                tempVariation = "ScaleWeights";
            else tempVariation = "JES";

            //if (tempVariation.Contains("JES")) continue;

            if (variation.Contains("pdfVariation99")) continue;
            if (variation.Contains("pdfVariation98")) continue;
            if (variation.Contains("pdfVariation100")) continue;

            /* cout<<TString::Format("%s/UnfoldedCombined/%s/OutputFile%s_%s.root",
                                                                baseDir.Data(),
                                                                tempVariation.Data(),
                                                                partonParticle.Data(),
                                                                variation.Data())<<endl; */
            
            cout<<TString::Format("combined_%s_%s", variation.Data(), variable.Data())<<endl;
            
            TH1F *variationHistogram = (TH1F*)nominalFile->Get(TString::Format("combined_%s_%s", 
                                                                            variation.Data(),
                                                                            variable.Data()));
            
            variationHistogram->Scale(1. / AnalysisConstants::luminositiesSR["comb"], "width");
            float_t variationHistogramYield = variationHistogram->Integral();
            if (normalized) variationHistogram->Scale(1 / variationHistogramYield);

            AddSystematicToErrorBar(finalResult, variationHistogram);
            /*
            if (tempVariation.Contains("JES")) continue;
            if (tempVariation.Contains("bTagVariation")) continue;
            if (tempVariation.Contains("PDFWeights")) continue;
            if (tempVariation.EqualTo("PSWeights"))
                cout<<variation<<endl;*/
            // get it from the unfolded file 
            TString tempVariationForTheory = AnalysisConstants::variations[var];
            TFile *variationFileTheory = TFile::Open(TString::Format("%s/UnfoldedCombined/%s/OutputFileParton_%s.root",
                                                                baseDir.Data(),
                                                                tempVariation.Data(),
                                                                tempVariationForTheory.Data()));
            
            TH1F *variationTheory = (TH1F *)variationFileTheory->Get(TString::Format("hFidTheory%s_%s",
                                                                            (normalized ? "Norm" : ""),
                                                                            variable.Data()));

            AddSystematicToErrorBar(finalTheory, variationTheory);
            variationFileTheory->Close();
            // check AMC@NLO pdf and scale 
            if (variation.Contains("pdf") ||
                variation.Contains("scale")
            )
            {   
                if (tempVariation.EqualTo("JES")) continue;
                if (tempVariation.EqualTo("PSWeights")) continue;
                
                TString tempVariationForTheory = AnalysisConstants::variations[var];
                TFile *variationTheoryAmcAtNloFile = TFile::Open(TString::Format(
                                        "../VariationHandling_Theory_amc@NLO/%s/%s/Histograms_TTJets_%s.root",
                                        theoryYear.Data(),
                                        tempVariation.Data(),
                                        tempVariationForTheory.Data()));
                // this histogram has the expected yield!!
                // thus, you need to scale to get xsec
                TH1F *variationTheoryAmcAtNlo = (TH1F *)variationTheoryAmcAtNloFile->Get(TString::Format("hReco_%s_%s", 
                                                                        variable.Data(),
                                                                        tempVariationForTheory.Data()));

                
                if (variation.Contains("pdf") ||
                    variation.Contains("scale"))
                {
                    variationTheoryAmcAtNlo->Scale(2.);
                } 
                if((/*variation.Contains("scale_3") || variation.Contains("scale_5") || */variation.Contains("scale_9")) 
                    && (AnalysisConstants::partonVariables[v].Contains("chi") || AnalysisConstants::partonVariables[v].Contains("cos")))
                    continue;
                
                variationTheoryAmcAtNlo->Scale(1. / AnalysisConstants::luminositiesSR["2018"], "width");
                
                float_t variationTheoryAmcAtNloYield = variationTheoryAmcAtNlo->Integral();
                // cout<<variationTheoryAmcAtNlo->GetBinContent(1)<<" "<<finalTheoryAmcAtNlo->GetBinContent(1)<<endl;
                // gb has already scaled
                if (normalized)
                {
                    variationTheoryAmcAtNlo->Scale(1 / variationTheoryAmcAtNloYield);
                }
                AddSystematicToErrorBar(finalTheoryAmcAtNlo, variationTheoryAmcAtNlo);
                variationTheoryAmcAtNloFile->Close();
            }
            
        }
        if (variable.Contains("jetY") && !normalized)
        {
            nominalHistogram->Scale(0.5);
            finalResult->Scale(0.5);
            finalTheory->Scale(0.5);
            theoryHistogram->Scale(0.5);
            finalTheoryAmcAtNlo->Scale(0.5);
            theoryAmcAtNloHistogram->Scale(0.5);
        }

        TH1F *theoryHistogramValue = (TH1F *)theoryHistogram->Clone(TString::Format("theoryHistogramValue%s",
                                                                                    (normalized ? "_normalized" : "")));
        TH1F *theoryAmcAtNloHistogramValue = (TH1F *)theoryAmcAtNloHistogram->Clone(TString::Format("theoryAmcAtNloHistogramValue%s",
                                                                                                    (normalized ? "_normalized" : "")));
        std::vector<TH1F *> histogramsToDraw;
        histogramsToDraw.push_back(finalResult);
        histogramsToDraw.push_back(nominalHistogram);
        histogramsToDraw.push_back(finalTheoryAmcAtNlo);
        histogramsToDraw.push_back(theoryAmcAtNloHistogramValue);
        histogramsToDraw.push_back(finalTheory);
        histogramsToDraw.push_back(theoryHistogramValue);

        DrawWithRatio(c1, histogramsToDraw, v, normalized, partonParticle);


        TString draw_name;
        if (partonParticle.EqualTo("Particle")) draw_name = AnalysisConstants::particleVariables[v].Data();
        else draw_name = AnalysisConstants::partonVariables[v].Data();
            

        /* c1->SaveAs(TString::Format("%s/FinalResult_%s%s.png",
                                outputDir.Data(),
                                draw_name.Data(),
                                (normalized ? "_normalized" : "")),
                "png");

        c1->SaveAs(TString::Format("%s/FinalResult_%s%s.svg",
                                outputDir.Data(),
                                variable.Data(),
                                (normalized ? "_normalized" : "")),
                "svg"); */
        
        c1->SaveAs(TString::Format("%s/FinalResult_%s%s.pdf",
                                outputDir.Data(),
                                draw_name.Data(),
                                (normalized ? "_normalized" : "")),
                "pdf");
    }
}
