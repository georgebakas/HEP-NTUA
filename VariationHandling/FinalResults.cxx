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
    std::vector<Color_t> colors = {kBlack, kBlack, kRed, kRed, kGreen, kGreen, kBlue, kBlue};
    std::vector<Color_t> fillColors = {kGray, kGray, kRed, kRed, kGreen, kGreen, kBlue, kBlue};
    std::vector<Int_t> markerStyle = {20, 1, 1, 1, 1, 1, 1, 1};
    std::vector<int> fillStyles = {1001, 1001, 3475, 3475, 3495, 3495, 3457, 3457};
    std::vector<TString> drawOptions = {"E2", "SAME EP", "SAME EP", "SAME", "SAME EP", "SAME", "SAME EP", "SAME"};
    std::vector<bool> drawRatio = {true, true, true, false, true, false, true, false};
    std::vector<bool> putInLegend = {true, true, true, false, true, false, true, false};
    std::vector<TString> legendTitles = {"Data", "Total unc.", "amc@NLO+Pythia8", "amc@NLO+Pythia8 unc", "Powheg+Herwig", "Powheg+Herwig unc", "Powheg+Pythia8", "Powheg+Pythia8 unc"};
    std::vector<TString> legendDrawOption = {"EP", "f", "EL", "f", "EL", "f", "EL", "f"};
    
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
            {
                hist->GetYaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::partonYAxisValues[index][0],
                                            AnalysisConstants::FinalResultsConstans::partonYAxisValues[index][1]);
                hist->GetYaxis()->SetTitle("#frac{d#sigma}{dX}");
            }
            else 
            {
                hist->GetYaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::partonYAxisValuesNormalized[index][0],
                                            AnalysisConstants::FinalResultsConstans::partonYAxisValuesNormalized[index][1]);
                hist->GetYaxis()->SetTitle("#frac{1}{#sigma}#frac{d#sigma}{dX}");
            }
        }
        else{
            if (!isNormalized)
            {
                hist->GetYaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::particleYAxisValues[index][0],
                                            AnalysisConstants::FinalResultsConstans::particleYAxisValues[index][1]);
                hist->GetYaxis()->SetTitle("#frac{d#sigma}{dX}");
            }
            else 
            {
                hist->GetYaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::particleYAxisValuesNormalized[index][0],
                                            AnalysisConstants::FinalResultsConstans::particleYAxisValuesNormalized[index][1]);
                hist->GetYaxis()->SetTitle("#frac{1}{#sigma}#frac{d#sigma}{dX}");
            }
        }
        hist->GetYaxis()->SetLabelSize(0.035);
        hist->GetYaxis()->SetTitleSize(0.042);
        hist->GetYaxis()->SetTitleOffset(1);

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

        if (i < 2) {
            hist->Draw(drawOptions[i]);
        }
        else {
            if (i % 2 == 0)
            {
                TGraphAsymmErrors *gr = new TGraphAsymmErrors(hist->GetNbinsX());

                for (int j = 1; j <= hist->GetNbinsX(); j++)
                {
                float binWidth = hist->GetBinWidth(j);
                float x = hist->GetBinLowEdge(j) + binWidth / 3. * (i / 2);
                if (x == hist->GetBinLowEdge(j + 1))
                {
                    x -= binWidth / 10.;
                }
                gr->SetPoint(j - 1, x, hist->GetBinContent(j));
                gr->SetPointEXhigh(j - 1, 0);
                gr->SetPointEXlow(j - 1, 0);
                gr->SetPointEYhigh(j - 1, hist->GetBinError(j));
                gr->SetPointEYlow(j - 1, hist->GetBinError(j));
                }
                gr->SetMarkerStyle(hist->GetMarkerStyle());
                gr->SetMarkerColor(hist->GetMarkerColor());
                gr->SetLineColor(hist->GetLineColor());
                gr->Draw(drawOptions[i]);
            }
        }
        // hist->Draw(drawOptions[i]);

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
            
            if (i < 2) {
                ratio->Draw(drawOptions[i]);
            }
            else {
                TGraphAsymmErrors *gr = new TGraphAsymmErrors(ratio->GetNbinsX());

                for (int j = 1; j <= ratio->GetNbinsX(); j++)
                {
                float binWidth = ratio->GetBinWidth(j);
                float x = ratio->GetBinLowEdge(j) + binWidth / 3. * (i / 2);
                if (x == ratio->GetBinLowEdge(j + 1))
                {
                    x -= binWidth / 10.;
                }
                gr->SetPoint(j - 1, x, ratio->GetBinContent(j));
                gr->SetPointEXhigh(j - 1, 0);
                gr->SetPointEXlow(j - 1, 0);
                gr->SetPointEYhigh(j - 1, ratio->GetBinError(j));
                gr->SetPointEYlow(j - 1, ratio->GetBinError(j));
                }
                gr->SetMarkerStyle(ratio->GetMarkerStyle());
                gr->SetMarkerColor(ratio->GetMarkerColor());
                gr->SetLineColor(ratio->GetLineColor());
                gr->Draw(drawOptions[i]);
            }
            if (i > 2 && (i % 2 != 0))
            {
            continue;
            }
            // ratio->Draw(drawOptions[i]);
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

    lumi_13TeV = TString::Format("%0.1f fb^{-1}", luminosity["luminosityAll"]/1000);
    int iPeriod = 13;
    int iPos = 0;
    extraTextFactor = 0.14;
    writeExtraText=true;
        
    CMS_lumi(upperPad, "combined", iPos);
    }

void CalcualteSystematicUncertainty(TH1F* nominal, TH1F *systematic, TH1F *systematicUp, TH1F* systematicDown)
{
    for (int bin = 0; bin < nominal->GetNbinsX(); bin++)
    {
        double nominalValue = nominal->GetBinContent(bin + 1);
        double nominalError = nominal->GetBinError(bin + 1);
        double variationValue = systematic->GetBinContent(bin + 1);
        double variationErrorUp = systematicUp->GetBinContent(bin + 1);
        double variationErrorDown = systematicDown->GetBinContent(bin + 1);
        double variation = variationValue - nominalValue;

        if(variation != 0)
        {
            if(variation > 0)
            {
                systematicUp->SetBinContent(bin + 1, variationErrorUp + (variation * variation));
            }
            else
            {
                systematicDown->SetBinContent(bin + 1, variationErrorUp + (variation * variation));
            }
        }
    }
}

void InitSystematicHistogram(TH1F* histogram, TH1F* systematicHistogram)
{
    systematicHistogram->Reset();
    for (int bin = 0; bin < histogram->GetNbinsX(); bin++)
    {
        systematicHistogram->SetBinContent(bin + 1, histogram->GetBinError(bin + 1) * histogram->GetBinError(bin + 1));
    }
}

void AddSystematicToErrorBar(TH1F *nominal, TH1F *systematic)
{
    for (int bin = 1; bin <= nominal->GetNbinsX(); bin++)
    {
        float nominalValue = nominal->GetBinContent(bin);
        float variationValue = systematic->GetBinContent(bin);
        float difference = TMath::Abs(nominalValue - variationValue);
        float error = TMath::Sqrt(TMath::Power(nominal->GetBinError(bin), 2) +
                                TMath::Power(difference, 2));
        nominal->SetBinError(bin, error);
    }
}

void CalculateTotalSystematicUncertainty(TH1F* histogram, TH1F* varationUp, TH1F* variationDown)
{
    for (int bin = 0; bin < histogram->GetNbinsX(); bin++)
    {
        double errorUp = TMath::Sqrt(varationUp->GetBinContent(bin + 1));
        double errorDown = TMath::Sqrt(variationDown->GetBinContent(bin + 1));
        histogram->SetBinError(bin + 1, 0.5 * (errorUp + errorDown));
    }
}


void FinalResults(TString partonParticle = "Parton", bool normalized = false)
{
    initFilesMapping();
    //setTDRStyle();
    gStyle->SetOptStat(0);
    //TString partonParticle = "Parton";
    TString baseDir = "/Users/georgebakas/Documents/HEP-NTUA_ul/VariationHandling";
    //TString baseDir = "/afs/cern.ch/work/g/gbakas/public/HEP-NTUA/VariationHandling";
    AnalysisConstants::initConstants();

    TFile *finalResultFile = TFile::Open(TString::Format("%s/FinalResults/results/%s_outputFile%s.root",
                                                    baseDir.Data(),
                                                    partonParticle.Data(),
                                                    (normalized ? "Norm" : "")), "RECREATE");

    TFile *nominalFile = TFile::Open(TString::Format("%s/UnfoldedCombined/Nominal/OutputFile%s.root",
                                                    baseDir.Data(),
                                                    partonParticle.Data()));
    TString outputDir = TString::Format("%s/FinalResults/results",
                                        baseDir.Data());
    CheckAndCreateDirectory("results");

    TString theoryYear = "2018";
    TFile *nominalAmcAtNloFile = TFile::Open(TString::Format(
                            "../VariationHandling_Theory_amc@NLO/%s/Nominal/Histograms_TTJets.root", theoryYear.Data()));
    TFile *nominalHerwigFile = TFile::Open(TString::Format(
                            "../VariationHandling_Theory_herwig/%s/Nominal/Histograms_TT.root", theoryYear.Data()));

    for (unsigned int v = 0; v < AnalysisConstants::unfoldingVariables.size(); v++)
    {   
        TString variable = AnalysisConstants::unfoldingVariables[v];

        TH1F *nominalHistogram = (TH1F *)nominalFile->Get(TString::Format("hUnfold%s_%s",
                                                                        (normalized ? "Norm" : "Final"), 
                                                                        variable.Data()));

        TH1F *finalResult = (TH1F *)nominalHistogram->Clone(TString::Format("FinalResult%s_%s",
                                                                        (normalized ? "Norm" : "Final"), 
                                                                        variable.Data()));

        TH1F *finalResultUp = (TH1F *)nominalHistogram->Clone(TString::Format("FinalResultErrorUp_%s%s",
                                                                        variable.Data(),
                                                                        (normalized ? "_normalized" : "")));
        InitSystematicHistogram(finalResult, finalResultUp);
        TH1F *finalResultDown = (TH1F *)nominalHistogram->Clone(TString::Format("FinalResultErrorDown_%s%s",
                                                                            variable.Data(),
                                                                            (normalized ? "_normalized" : "")));
        InitSystematicHistogram(finalResult, finalResultDown);

        TH1F *theoryHistogram = (TH1F *)nominalFile->Get(TString::Format("hTheory%s_%s",
                                                                        (normalized ? "Norm" : ""),
                                                                        variable.Data()));
        
        TH1F *finalTheory = (TH1F *)theoryHistogram->Clone(TString::Format("FinalTheory%s_%s",
                                                                        (normalized ? "Norm" : ""),
                                                                        variable.Data()));

        TH1F* theoryHistogrampUp = (TH1F *)theoryHistogram->Clone(TString::Format("FinalTheoryErrorUp_%s%s",
                                                                        variable.Data(),
                                                                        (normalized ? "_normalized" : "")));
        
        InitSystematicHistogram(theoryHistogram, theoryHistogrampUp);

        TH1F* theoryHistogrampDown = (TH1F *)theoryHistogram->Clone(TString::Format("FinalTheoryErrorDown_%s%s",
                                                                        variable.Data(),
                                                                        (normalized ? "_normalized" : "")));
        InitSystematicHistogram(theoryHistogram, theoryHistogrampDown);

        // AMC@NLO 
        TH1F *theoryAmcAtNloHistogram, *theoryAmcAtNloHistogramUp, *theoryAmcAtNloHistogramDown;
        if (partonParticle.EqualTo("Parton")) theoryAmcAtNloHistogram = (TH1F *)nominalAmcAtNloFile->Get(TString::Format("hParton_%s", 
                                                                AnalysisConstants::partonVariables[v].Data()));
        else theoryAmcAtNloHistogram = (TH1F *)nominalAmcAtNloFile->Get(TString::Format("hParticle_%s", 
                                        AnalysisConstants::particleVariables[v].Data()));
        
        theoryAmcAtNloHistogram->Scale(1. / AnalysisConstants::luminositiesSR["2018"], "width");

        // Herwig samples
        TH1F *theoryHerwigHistogram, *theoryHerwigHistogramDown, *theoryHerwigHistogramUp;
        if (partonParticle.EqualTo("Parton")) theoryHerwigHistogram = (TH1F *)nominalHerwigFile->Get(TString::Format("hParton_%s", 
                                                                AnalysisConstants::partonVariables[v].Data()));
        else theoryHerwigHistogram = (TH1F *)nominalHerwigFile->Get(TString::Format("hParticle_%s", 
                                        AnalysisConstants::particleVariables[v].Data()));
        
        theoryHerwigHistogram->Scale(1. / AnalysisConstants::luminositiesSR["2018"], "width");

        float_t theoryAmcAtNloYield = theoryAmcAtNloHistogram->Integral();
        float_t herwigYield = theoryHerwigHistogram->Integral();
        if (normalized)
        {
            theoryAmcAtNloHistogram->Scale(1. / theoryAmcAtNloYield);
            theoryHerwigHistogram->Scale(1. / herwigYield);
        }

        TH1F *finalTheoryAmcAtNlo = (TH1F *)theoryAmcAtNloHistogram->Clone(TString::Format("FinalTheoryAmcAtNlo%s_%s",
                                                                            (normalized ? "Norm" : ""),
                                                                            variable.Data()));

        theoryAmcAtNloHistogramUp = (TH1F *)theoryAmcAtNloHistogram->Clone(TString::Format("FinalTheoryAmcAtNloErrorUp_%s%s",
                                                                                            variable.Data(),
                                                                                            (normalized ? "_normalized" : "")));
        InitSystematicHistogram(theoryAmcAtNloHistogram, theoryAmcAtNloHistogramUp);

        theoryAmcAtNloHistogramDown = (TH1F *)theoryAmcAtNloHistogram->Clone(TString::Format("FinalTheoryAmcAtNloErrorDown_%s%s",
                                                                                                variable.Data(),
                                                                                                (normalized ? "_normalized" : "")));
        InitSystematicHistogram(theoryAmcAtNloHistogram, theoryAmcAtNloHistogramDown);
        
        //hewrig histogram 
        TH1F *finalTheoryHerwig = (TH1F *)theoryHerwigHistogram->Clone(TString::Format("FinalTheoryHerwig%s_%s",
                                                                            (normalized ? "Norm" : ""),
                                                                            variable.Data()));

        theoryHerwigHistogramUp = (TH1F *)theoryHerwigHistogram->Clone(TString::Format("FinalTheoryHerwigErrorUp_%s%s",
                                                                                            variable.Data(),
                                                                                            (normalized ? "_normalized" : "")));
        InitSystematicHistogram(theoryHerwigHistogram, theoryHerwigHistogramUp);
        theoryHerwigHistogramDown = (TH1F *)theoryHerwigHistogram->Clone(TString::Format("FinalTheoryHerwigErrorDown_%s%s",
                                                                                                variable.Data(),
                                                                                                (normalized ? "_normalized" : "")));
        InitSystematicHistogram(theoryHerwigHistogram, theoryHerwigHistogramDown);

        TCanvas *c1 = new TCanvas(TString::Format("FinalResult_%s",
                                                variable.Data()),
                                TString::Format("FinalResult_%s",
                                                variable.Data()),
                                600, 600);
        
        for (unsigned int var = 0; var < AnalysisConstants::variations.size(); var++)
        {
            TString variation = AnalysisConstants::variations[var];
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

            if (variation.Contains("pdf_99")) continue;
            if (variation.Contains("pdf_98")) continue;
            if (variation.Contains("pdf_100")) continue;
            if (variation.Contains("pdf_1")) continue;
            if (variation.Contains("pdf_35")) continue;
            if (variation.Contains("pdf_40")) continue;
            if (variation.Contains("pdf_56")) continue;
            if (variation.Contains("pdf_57")) continue;
            if (variation.Contains("pdf_59")) continue;
            if (variation.Contains("pdf_62")) continue;
            if (variation.Contains("pdf_83")) continue;
            if (variation.Contains("pdf_86")) continue;
            if (variation.Contains("pdf_85")) continue;
            if (variation.Contains("pdf_76")) continue;
            if (variation.Contains("pdf_42")) continue;
            if (variation.Contains("pdf_58")) continue;
            if (variation.Contains("pdf_69")) continue;
            if (variation.Contains("pdf_93")) continue; 

            cout<<TString::Format("%s/UnfoldedCombined/%s/OutputFile%s_%s.root",
                                                                baseDir.Data(),
                                                                tempVariation.Data(),
                                                                partonParticle.Data(),
                                                                variation.Data())<<endl;
            TFile *variationFile = TFile::Open(TString::Format("%s/UnfoldedCombined/%s/OutputFile%s_%s.root",
                                                                baseDir.Data(),
                                                                tempVariation.Data(),
                                                                partonParticle.Data(),
                                                                variation.Data()));

            TH1F *variationHistogram = (TH1F *)variationFile->Get(TString::Format("hUnfold%s_%s",
                                                                            (normalized ? "Norm" : "Final"), 
                                                                            variable.Data()));
            // AddSystematicToErrorBar(finalResult, variationHistogram);
            CalcualteSystematicUncertainty(finalResult, variationHistogram, finalResultUp, finalResultDown);

            // check AMC@NLO pdf and scale 
            if ( variation.Contains("pdf") ||
                variation.Contains("scale") || variation.Contains("ps"))
            {   
                // if (tempVariation.EqualTo("PSWeights")) continue;
                TH1F *variationTheory = (TH1F *)variationFile->Get(TString::Format("hTheory%s_%s",
                                                                            (normalized ? "Norm" : ""),
                                                                            variable.Data()));

                // AddSystematicToErrorBar(finalTheory, variationTheory);
                CalcualteSystematicUncertainty(finalTheory, variationTheory, theoryHistogrampUp, theoryHistogrampDown);
                variationFile->Close();


                TFile *variationTheoryAmcAtNloFile = TFile::Open(TString::Format(
                                        "../VariationHandling_Theory_amc@NLO/%s/%s/Histograms_TTJets_%s.root",
                                        theoryYear.Data(),
                                        tempVariation.Data(),
                                        variation.Data()));

                TH1F *variationTheoryAmcAtNlo;
                if (partonParticle.EqualTo("Parton")) variationTheoryAmcAtNlo = (TH1F *)variationTheoryAmcAtNloFile->Get(TString::Format("hParton_%s_%s", 
                                                                        AnalysisConstants::partonVariables[v].Data(),
                                                                        variation.Data()));
                else variationTheoryAmcAtNlo = (TH1F *)variationTheoryAmcAtNloFile->Get(TString::Format("hParticle_%s_%s", 
                                                AnalysisConstants::particleVariables[v].Data(),
                                                variation.Data()));

                //if((/*variation.Contains("scale_3") || variation.Contains("scale_5") || */variation.Contains("scale_9")) 
                //    && (AnalysisConstants::partonVariables[v].Contains("chi") || AnalysisConstants::partonVariables[v].Contains("cos")))
                //    continue; 
                
                variationTheoryAmcAtNlo->Scale(1. / AnalysisConstants::luminositiesSR["2018"], "width");
                
                float_t variationTheoryAmcAtNloYield = variationTheoryAmcAtNlo->Integral();
                cout<<variationTheoryAmcAtNlo->GetBinContent(1)<<" "<<finalTheoryAmcAtNlo->GetBinContent(1)<<endl;
                // gb has already scaled
                if (normalized)
                {
                    variationTheoryAmcAtNlo->Scale(1. / variationTheoryAmcAtNloYield);
                }
                // AddSystematicToErrorBar(finalTheoryAmcAtNlo, variationTheoryAmcAtNlo);
                CalcualteSystematicUncertainty(finalTheoryAmcAtNlo, variationTheoryAmcAtNlo, theoryAmcAtNloHistogramUp, theoryAmcAtNloHistogramDown);
                variationTheoryAmcAtNloFile->Close();
                
                // --------------------------------------------------------------------------------------------------------------------------------------------------                
                //now the same for the herwig samples
                // --------------------------------------------------------------------------------------------------------------------------------------------------                
                TFile *variationTheoryHerwigFile = TFile::Open(TString::Format(
                                        "../VariationHandling_Theory_Herwig/%s/%s/Histograms_TT_%s.root",
                                        theoryYear.Data(),
                                        tempVariation.Data(),
                                        variation.Data()));
                TH1F *variationTheoryHerwig;
                if (partonParticle.EqualTo("Parton")) variationTheoryHerwig = (TH1F *)variationTheoryHerwigFile->Get(TString::Format("hParton_%s_%s", 
                                                                        AnalysisConstants::partonVariables[v].Data(),
                                                                        variation.Data()));
                else variationTheoryHerwig = (TH1F *)variationTheoryHerwigFile->Get(TString::Format("hParticle_%s_%s", 
                                                AnalysisConstants::particleVariables[v].Data(),
                                                variation.Data()));
                
                variationTheoryHerwig->Scale(1. / AnalysisConstants::luminositiesSR["2018"], "width");
                
                float_t variationTheoryHerwigYield = variationTheoryHerwig->Integral();
                // gb has already scaled
                if (normalized)
                {
                    variationTheoryHerwig->Scale(1. / variationTheoryHerwigYield);
                }
                // AddSystematicToErrorBar(finalTheoryHerwig, variationTheoryHerwig);
                CalcualteSystematicUncertainty(finalTheoryHerwig, variationTheoryHerwig, theoryHerwigHistogramUp, theoryHerwigHistogramDown);
                variationTheoryHerwigFile->Close();
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
            finalTheoryHerwig->Scale(0.5);
            theoryHerwigHistogram->Scale(0.5);
        }
        c1->cd();
        TH1F *theoryHistogramValue = (TH1F *)theoryHistogram->Clone(TString::Format("theoryHistogramValue%s",
                                                                                    (normalized ? "_normalized" : "")));
        TH1F *theoryAmcAtNloHistogramValue = (TH1F *)theoryAmcAtNloHistogram->Clone(TString::Format("theoryAmcAtNloHistogramValue%s",
                                                                                                    (normalized ? "_normalized" : "")));
        TH1F *theoryHerwigHistogramValue = (TH1F *)theoryHerwigHistogram->Clone(TString::Format("theoryHerwigHistogramValue%s",
                                                                                                    (normalized ? "_normalized" : "")));
        std::vector<TH1F *> histogramsToDraw;
        CalculateTotalSystematicUncertainty(finalResult, finalResultUp, finalResultDown);
        histogramsToDraw.push_back(finalResult);
        histogramsToDraw.push_back(nominalHistogram);

        // now amc@nlo 
        CalculateTotalSystematicUncertainty(finalTheoryAmcAtNlo, theoryAmcAtNloHistogramUp, theoryAmcAtNloHistogramDown);
        histogramsToDraw.push_back(finalTheoryAmcAtNlo);
        histogramsToDraw.push_back(theoryAmcAtNloHistogramValue);

        // now herwig 
        CalculateTotalSystematicUncertainty(finalTheoryHerwig, theoryHerwigHistogramUp, theoryHerwigHistogramDown);
        histogramsToDraw.push_back(finalTheoryHerwig);
        histogramsToDraw.push_back(theoryHerwigHistogramValue);

        // now theory 
        CalculateTotalSystematicUncertainty(finalTheory, theoryHistogrampUp, theoryHistogrampDown);
        histogramsToDraw.push_back(finalTheory);
        histogramsToDraw.push_back(theoryHistogramValue);
        
        DrawWithRatio(c1, histogramsToDraw, v, normalized, partonParticle);
        
        TString draw_name;
        if (partonParticle.EqualTo("Particle")) draw_name = AnalysisConstants::particleVariables[v].Data();
        else draw_name = AnalysisConstants::partonVariables[v].Data();
        
        c1->SaveAs(TString::Format("%s/FinalResult_%s%s.pdf",
                                outputDir.Data(),
                                draw_name.Data(),
                                (normalized ? "_normalized" : "")),
                "pdf");

        finalResultFile->cd();
        finalResult->Write((TString::Format("FinalResult_%s", variable.Data())));
        nominalHistogram->Write((TString::Format("nominalHistogram%s", variable.Data())));
        finalTheoryAmcAtNlo->Write((TString::Format("finalTheoryAmcAtNlo%s", variable.Data())));
        finalTheoryHerwig->Write((TString::Format("finalTheoryHerwig%s", variable.Data())));
        theoryAmcAtNloHistogramValue->Write((TString::Format("theoryAmcAtNloHistogramValue%s", variable.Data())));
        theoryHerwigHistogramValue->Write((TString::Format("theoryHerwigHistogramValue%s", variable.Data())));
        finalTheory->Write((TString::Format("finalTheory%s", variable.Data())));
        theoryHistogramValue->Write((TString::Format("theoryHistogramValue%s", variable.Data())));
    }
    finalResultFile->Close();
}
