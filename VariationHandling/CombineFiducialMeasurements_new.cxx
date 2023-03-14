//#include "BASE.h"

#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TCanvas.h> 
#include <TLegend.h>
#include <TMath.h> 

#include "FinalResultsConstants.h"

#include "../CMS_plots/CMS_lumi.C"
#include "../CMS_plots/tdrstyle.C"

//Combine cross section ta the fiducial level without using the Blue method 

/*
    Combination Strategy:

    1. Get Nominal files and add the bulkly
    2. Create TH1F vectors for the uncertainties -> 1 vector for each year
    3. list all files of a current year (2016_postVFP) and get THESE as templates 
        This way we will ensure that all files are added correctly.
    
*/
void DrawWithRatio(TCanvas *can, std::vector<TH1F *> histograms, int index, bool isNormalized)
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
        if (!isNormalized)
            hist->GetYaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::partonYAxisValues[index][0],
                                        AnalysisConstants::FinalResultsConstans::partonYAxisValues[index][1]);
        else 
            hist->GetYaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::partonYAxisValuesNormalized[index][0],
                                        AnalysisConstants::FinalResultsConstans::partonYAxisValuesNormalized[index][1]);

        TH1F *ratio = (TH1F *)hist->Clone("ratio");
        ratio->Divide(dataHist);
        for (int bin = 1; bin <= hist->GetNbinsX(); bin++)
        {
        ratio->SetBinContent(bin, ratio->GetBinContent(bin) - 1);
        }
        hist->GetXaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::partonXAxisValues[index][0],
                                    AnalysisConstants::FinalResultsConstans::partonXAxisValues[index][1]);
        
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
        ratio->GetXaxis()->SetTitle(AnalysisConstants::partonAxisTitles[index]);
        ratio->GetXaxis()->SetRangeUser(AnalysisConstants::FinalResultsConstans::partonXAxisValues[index][0],
                                    AnalysisConstants::FinalResultsConstans::partonXAxisValues[index][1]);
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
    CMS_lumi(upperPad, "combined", iPos);
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


std::vector<TString> listFiles(const char *dirname="", const char *var="", const char *ext=".root")
{
    std::vector<TString> list_of_files;
    TSystemDirectory dir(dirname, dirname);
    TList *files = dir.GetListOfFiles();
    if (files) {
        TSystemFile *file;
        TString fname;
        TIter next(files);
        while ((file=(TSystemFile*)next())) {
        fname = file->GetName();
        if (!file->IsDirectory() && fname.EndsWith(ext) && fname.Contains(var)) {
         //cout << fname.Data() << endl;
            list_of_files.push_back(fname.Data());
        }
    }
    }
    return list_of_files;
}


void CombineFiducialMeasurements_new()
{
    TFile *outFile = new TFile("testFile_.root", "RECREATE");
    AnalysisConstants::initConstants();

    TString baseInputDir = "/afs/cern.ch/work/g/gbakas/public/HEP-NTUA/";
    //TString baseInputDir = "/Users/georgebakas/Documents/HEP-NTUA_ul/";
    baseInputDir = TString::Format("%s/VariationHandling/", baseInputDir.Data());

    // Define formats for Figures and Latex file
    const TString ForVal = "%1.6f";
    const TString ForUnc = "%1.6f";
    const TString ForWei = "%1.3f";
    const TString ForRho = "%1.2f";
    const TString ForPul = ForRho;
    const TString ForUni = "pb";
    std::vector<TString> variation_dirs = {"Nominal", "NominalTopSF", "PSWeights", "bTagVariation", "topTaggingVariation", "JES",  "ScaleWeights", "PDFWeights"/*"SystematicsFiles",*/ };

    const int num_years = 4; 
    TString NamEst[num_years];
    for (unsigned int i = 0; i < AnalysisConstants::years.size(); i++)
    {
        NamEst[i] = AnalysisConstants::years[i];
    }
    
    for (int ivar = 0; ivar<=AnalysisConstants::variations.size(); ivar++)
    {
        cout<<AnalysisConstants::variations[ivar]<<endl;
    }
    // loop on all variables 
    for (unsigned int var = 0; var < AnalysisConstants::unfoldingVariables.size(); var++)
    {
        std::vector<TH1F *> originalHistograms;
    
        std::vector<TH1F *> uncHistograms16pre;
        std::vector<TH1F *> uncHistograms16post;
        std::vector<TH1F *> uncHistograms17;
        std::vector<TH1F *> uncHistograms18;

        TString variable = AnalysisConstants::unfoldingVariables[var];
        cout<<"variable: "<<variable << endl;

        TCanvas *c1 = new TCanvas(TString::Format("FinalFiducialResult_%s",
                                                variable.Data()),
                                TString::Format("FinalFiducialResult_%s",
                                                variable.Data()),
                                600, 600);

        //loop on each year 
        for (unsigned int y = 0; y < AnalysisConstants::years.size(); y++)
        {
            NamEst[y] = AnalysisConstants::years[y];  

            // Step 1 is to handle the nominal histograms       
            TFile *file = TFile::Open(TString::Format("%s/%s/Nominal/FiducialMeasurement/SignalHistograms_%s_MassFitResults_SignalTemplates_.root",
                                                baseInputDir.Data(),
                                                AnalysisConstants::years[y].Data(),
                                                variable.Data()));
            TH1F *f = (TH1F *)file->Get(TString::Format("hSignal_%s",
                                                variable.Data()));
            f->SetDirectory(0);
            f->SetName(TString::Format("%s_%s",
                                        variable.Data(),
                                        AnalysisConstants::years[y].Data()));
            originalHistograms.push_back(f);
            file->Close();

            //loop on all variation directories 
            for (unsigned int i = 1; i < variation_dirs.size(); i++)
            {
                //search on a STATIC year for the vaiation name!!
                TString testYear = "2016_postVFP";
                std::vector<TString> variationFiles = listFiles(TString::Format("%s/%s/%s/FiducialMeasurement/", 
                                                    baseInputDir.Data(), 
                                                    testYear.Data(),
                                                    variation_dirs[i].Data()), 
                                                    variable.Data());
                
                for (int jvar=0; jvar<variationFiles.size(); jvar++)
                {
                    //cout<<variationFiles[jvar] <<endl;
                    if (variationFiles[jvar].Contains("nom0") || variationFiles[jvar].Contains("nom1")) continue;
                    TFile *file = TFile::Open(TString::Format("%s/%s/%s/FiducialMeasurement/%s",
                                                            baseInputDir.Data(),
                                                            AnalysisConstants::years[y].Data(),
                                                            variation_dirs[i].Data(),
                                                            variationFiles[jvar].Data()));

                    f = (TH1F *)file->Get(TString::Format("hSignal_%s",
                                                AnalysisConstants::unfoldingVariables[var].Data())); 
                    
                    //do this to get coherent name
                    TString hname = variationFiles[jvar];
                    hname.ReplaceAll(TString::Format("SignalHistograms_%s_MassFitResults_SignalTemplates_", variable.Data()), "");
                    hname.ReplaceAll(".root", "");
                    
                    if (variation_dirs[i].Contains("PDF")) hname = TString::Format("pdfVariation%s", hname.Data());
                    if (variation_dirs[i].Contains("Scale")) hname = TString::Format("scaleWeight%s", hname.Data());
                    if (variation_dirs[i].Contains("bTag")) hname = TString::Format("bTag%s", hname.Data());

                    f->SetDirectory(0);
                    f->SetName(TString::Format("%s_%s_%s",
                                                hname.Data(),
                                                AnalysisConstants::unfoldingVariables[var].Data(),
                                                AnalysisConstants::years[y].Data()));

                    if (y==0)uncHistograms16pre.push_back(f);
                    else if (y==1)uncHistograms16post.push_back(f);
                    else if (y==2)uncHistograms17.push_back(f);
                    else if (y==3)uncHistograms18.push_back(f);
                    file->Close();                       
                } //loop on variation directory files 

            } // loop on variation directories
        } //loop on years

        cout<< uncHistograms16pre.size()<<endl;
        cout<< uncHistograms16post.size()<<endl;
        cout<< uncHistograms17.size()<<endl;
        cout<< uncHistograms18.size()<<endl;
        // get the bins
        Float_t *bins = GetHistogramBins(originalHistograms[0]);
        // create a combined results histogram for each variable 
        TH1F *resultsHisto = new TH1F(TString::Format("combined_%s",
                                                    AnalysisConstants::unfoldingVariables[var].Data()),
                                    TString::Format("combined_%s",
                                                    AnalysisConstants::unfoldingVariables[var].Data()),
                                    originalHistograms[0]->GetNbinsX(),
                                    bins);

        cout<< "--------------------------------" << endl;
        for (int ibin = 1; ibin <=originalHistograms[0]->GetNbinsX(); ibin++)
        {
            resultsHisto -> SetBinContent(ibin, originalHistograms[0]->GetBinContent(ibin) + 
                                    originalHistograms[1]->GetBinContent(ibin) +
                                    originalHistograms[2]->GetBinContent(ibin) +
                                    originalHistograms[3]->GetBinContent(ibin));
            resultsHisto -> SetBinError(ibin, TMath::Sqrt(TMath::Power(originalHistograms[0]->GetBinError(ibin), 2) + 
                                    TMath::Power(originalHistograms[1]->GetBinError(ibin), 2) +
                                    TMath::Power(originalHistograms[2]->GetBinError(ibin), 2) +
                                    TMath::Power(originalHistograms[3]->GetBinError(ibin), 2)) + 
                                    // I have to input here the correlation coefficients 
                                    2*(AnalysisConstants::correlations["Nominal"]).correlations[1]*originalHistograms[0]->GetBinError(ibin)*originalHistograms[1]->GetBinError(ibin) + 
                                    2*(AnalysisConstants::correlations["Nominal"]).correlations[2]*originalHistograms[0]->GetBinError(ibin)*originalHistograms[2]->GetBinError(ibin) + 
                                    2*(AnalysisConstants::correlations["Nominal"]).correlations[3]*originalHistograms[0]->GetBinError(ibin)*originalHistograms[3]->GetBinError(ibin) + 
                                    2*(AnalysisConstants::correlations["Nominal"]).correlations[6]*originalHistograms[1]->GetBinError(ibin)*originalHistograms[2]->GetBinError(ibin) + 
                                    2*(AnalysisConstants::correlations["Nominal"]).correlations[7]*originalHistograms[1]->GetBinError(ibin)*originalHistograms[3]->GetBinError(ibin) + 
                                    2*(AnalysisConstants::correlations["Nominal"]).correlations[11]*originalHistograms[2]->GetBinError(ibin)*originalHistograms[3]->GetBinError(ibin));        
            cout<< resultsHisto ->GetBinContent(ibin) << " with error "<<resultsHisto->GetBinError(ibin)<<endl;
        }
        
        //now deal with the variations 

        TString testYear = "2016_postVFP";
        std::vector<TH1F *> uncHistograms_total;
        for (int unc = 0; unc<uncHistograms16post.size(); unc++)
        { 
            TString hname = uncHistograms16post[unc]->GetName();
            hname.ReplaceAll("_2016_postVFP", "");
            
            Float_t *bins = GetHistogramBins(originalHistograms[0]);
            TH1F *f = new TH1F(TString::Format("combined_%s",
                                        hname.Data()),
                                TString::Format("combined_%s",
                                    hname.Data()),
                                    originalHistograms[0]->GetNbinsX(),
                                    bins);
            uncHistograms_total.push_back(f);
        }

        // per variation 
        //cout<< AnalysisConstants::variations.size() << endl;
        for (unsigned int i = 0; i < AnalysisConstants::variations.size(); i++)
        {
            //cout<< "------------------" <<endl;
            //cout<< "16pre: "<< uncHistograms16pre[i]->GetName() << endl;
            //cout<< "16post: "<< uncHistograms16post[i]->GetName() << endl;
            //cout<< "17: "<< uncHistograms17[i]->GetName() << endl;
            //cout<< "18: "<< uncHistograms18[i]->GetName() << endl;
            // per bin 
            for (int bin = 1; bin <= originalHistograms[0]->GetNbinsX(); bin++)
            {
                //this is to get it per year while looping on years
                uncHistograms_total[i] -> SetBinContent(bin, 
                                            uncHistograms16pre[i]->GetBinContent(bin) +
                                            uncHistograms16post[i]->GetBinContent(bin) +
                                            uncHistograms17[i]->GetBinContent(bin) +
                                            uncHistograms18[i]->GetBinContent(bin));

                uncHistograms_total[i] -> SetBinError(bin, TMath::Sqrt(TMath::Power(uncHistograms16pre[i]->GetBinError(bin), 2) + 
                                            TMath::Power(uncHistograms16post[i]->GetBinError(bin), 2) +
                                            TMath::Power(uncHistograms17[i]->GetBinError(bin), 2) +
                                            TMath::Power(uncHistograms18[i]->GetBinError(bin), 2)) + 
                                            // I have to input here the correlation coefficients 
                                            2*(AnalysisConstants::correlations[AnalysisConstants::variations[i]]).correlations[1]*uncHistograms16pre[i]->GetBinError(bin)*uncHistograms16post[i]->GetBinError(bin) + 
                                            2*(AnalysisConstants::correlations[AnalysisConstants::variations[i]]).correlations[2]*uncHistograms17[i]->GetBinError(bin)*uncHistograms17[i]->GetBinError(bin) + 
                                            2*(AnalysisConstants::correlations[AnalysisConstants::variations[i]]).correlations[3]*uncHistograms16pre[i]->GetBinError(bin)*uncHistograms18[i]->GetBinError(bin) + 
                                            2*(AnalysisConstants::correlations[AnalysisConstants::variations[i]]).correlations[6]*uncHistograms16post[i]->GetBinError(bin)*uncHistograms17[i]->GetBinError(bin) + 
                                            2*(AnalysisConstants::correlations[AnalysisConstants::variations[i]]).correlations[7]*uncHistograms16post[i]->GetBinError(bin)*uncHistograms18[i]->GetBinError(bin) + 
                                            2*(AnalysisConstants::correlations[AnalysisConstants::variations[i]]).correlations[11]*uncHistograms17[i]->GetBinError(bin)*uncHistograms18[i]->GetBinError(bin));
            }
        } 

        outFile->cd();
        resultsHisto->Write();
        for (int unc = 0; unc<AnalysisConstants::variations.size(); unc++)
        {
            //cout<<"NAME "<<uncHistograms_total[unc]->GetName()<<endl;
            uncHistograms_total[unc]->Write();
        }

        delete bins;
        delete resultsHisto;
        
        uncHistograms16pre.clear();
        uncHistograms16post.clear();
        uncHistograms17.clear();
        uncHistograms18.clear();
        originalHistograms.clear();

    } // loop on all variables    


    //DrawWithRatio(c1, histogramsToDraw, v, normalized);

    /*1->SaveAs(TString::Format("%s/FinalResult_%s%s.pdf",
                                outputDir.Data(),
                                draw_name.Data(),
                                (normalized ? "_normalized" : "")),
                "pdf");*/

    outFile->Close();
} //oef 

