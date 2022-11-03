#include "BASE.h"

#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TMath.h>
//Combine cross section ta the fiducial level without using the Blue method 

/*
    Combination Strategy:

    1. Get Nominal files and add the bulkly
    2. Create TH1F vectors for the uncertainties -> 1 vector for each year
    3. list all files of a current year (2016_postVFP) and get THESE as templates 
        This way we will ensure that all files are added correctly.
    
*/

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


void CombineMeasurements(bool isParton = true)
{   
    TString fname = "Particle";
    if (isParton)
        fname = "Parton";
    
    TFile *outFile = new TFile(TString::Format("testFile_Theory%s.root", fname.Data()), "RECREATE");
    AnalysisConstants::initConstants();

     //TString baseInputDir = "/afs/cern.ch/work/g/gbakas/public/HEP-NTUA/";
    TString baseInputDir = "/Users/georgebakas/Documents/HEP-NTUA_ul/";
    baseInputDir = TString::Format("%s/VariationHandling_Theory_amc@NLO/", baseInputDir.Data());
    const int NVAR = 10;
    TString variable[NVAR]   = {"jetPt0","jetPt1","yJJ", "ptJJ", "mJJ", "jetY0", "jetY1", "chi", "cosTheta_0", "cosTheta_1"};
    TString variableParton[NVAR] = {"partonPt0", "partonPt1", "yTTbarParton", "ptTTbarParton", "mTTbarParton", "partonY0", "partonY1","chiParton", "cosThetaParton_0", "cosThetaParton_1"};
    TString variableGen[NVAR] = {"genjetPt0", "genjetPt1", "yJJGen", "ptJJGen", "mJJGen", "genjetY0", "genjetY1", "chiParticle", "cosThetaParticle_0", "cosThetaParticle_1"};

    // Define formats for Figures and Latex file
    const TString ForVal = "%1.6f";
    const TString ForUnc = "%1.6f";
    const TString ForWei = "%1.3f";
    const TString ForRho = "%1.2f";
    const TString ForPul = ForRho;
    const TString ForUni = "pb";
    std::vector<TString> variation_dirs = {"Nominal", "bTagVariation",  "ScaleWeights", "PDFWeights"};

    const int num_years = 4; 
    TString NamEst[num_years];
    float LUMI = 0;
    for (unsigned int i = 0; i < AnalysisConstants::years.size(); i++)
    {
        NamEst[i] = AnalysisConstants::years[i];
        LUMI += AnalysisConstants::luminositiesSR[AnalysisConstants::years[i]];
    }
    
    for (int ivar = 0; ivar<=AnalysisConstants::variations.size(); ivar++)
    {
        cout<<ivar << " " <<AnalysisConstants::variations[ivar]<<endl;
    }

    
    // loop on all variables 
    for (unsigned int var = 0; var < AnalysisConstants::unfoldingVariables.size(); var++)
    {
        std::vector<TH1F *> originalHistograms;
        std::vector<TH1F *> uncHistograms16pre;
        std::vector<TH1F *> uncHistograms16post;
        std::vector<TH1F *> uncHistograms17;
        std::vector<TH1F *> uncHistograms18;
        TString tempVar;
        TString histoName;
        if(isParton)
        {   
            histoName = "hParton_";
            tempVar = variableParton[var];
        }
        else
        {
            histoName = "hParticle_";
            tempVar = variableGen[var];
        }
        TString variable = tempVar;
        cout<<"variable: "<<variable << endl;

        //loop on each year 
        for (unsigned int y = 0; y < -1; y++)
        {
            NamEst[y] = "2018";  

            // Step 1 is to handle the nominal histograms       
            TFile *file = TFile::Open(TString::Format("%s/%s/Nominal/Histograms_TTJets.root",
                                                baseInputDir.Data(),
                                                AnalysisConstants::years[y].Data()));
            TH1F *f = (TH1F *)file->Get(TString::Format("%s%s",
                                                histoName.Data(), 
                                                tempVar.Data()));

            cout<< TString::Format("%s/%s/Nominal/Histograms_TTJets.root",
                                                baseInputDir.Data(),
                                                AnalysisConstants::years[y].Data()) <<endl;
            cout<< TString::Format("%s%s",
                                    histoName.Data(), tempVar.Data()) << endl;
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
                TString testYear = "2018";
                std::vector<TString> variationFiles = listFiles(TString::Format("%s/%s/%s/", 
                                                    baseInputDir.Data(), 
                                                    testYear.Data(),
                                                    variation_dirs[i].Data()), ".root");
                cout<<variation_dirs[i]<< " "<<AnalysisConstants::years[y]<<endl;
                for (int jvar=0; jvar<variationFiles.size(); jvar++)
                {
                    //cout<<variationFiles[jvar] <<endl;
                    if (variationFiles[jvar].Contains("nom0") || variationFiles[jvar].Contains("nom1")) continue;
                    TFile *file = TFile::Open(TString::Format("%s/%s/%s/%s",
                                                            baseInputDir.Data(),
                                                            AnalysisConstants::years[y].Data(),
                                                            variation_dirs[i].Data(),
                                                            variationFiles[jvar].Data()));
                    cout<< TString::Format("%s/%s/%s/%s",
                                                            baseInputDir.Data(),
                                                            AnalysisConstants::years[y].Data(),
                                                            variation_dirs[i].Data(),
                                                            variationFiles[jvar].Data())<<endl;

                    //do this to get coherent name
                    TString hname = variationFiles[jvar];
                    hname.ReplaceAll("Histograms_TTJets_", "");
                    hname.ReplaceAll(".root", "");
                    cout<<"here"<<endl;
                    if (variation_dirs[i].Contains("bTag"))
                    {   
                        f = (TH1F *)file->Get(TString::Format("%s%s",
                                                    histoName.Data(),
                                                    tempVar.Data()));
                        cout<<TString::Format("%s%s",histoName.Data(), 
                                                    tempVar.Data())<<endl;
                    }
                    else
                    {
                        f = (TH1F *)file->Get(TString::Format("%s%s_%s",
                                                    histoName.Data(),
                                                    tempVar.Data(), 
                                                    hname.Data()));
                        cout<<TString::Format("%s%s_%s",histoName.Data(), 
                                                    tempVar.Data(), 
                                                    hname.Data())<<endl;
                    }
                    if (variation_dirs[i].Contains("PDF"))
                    {   
                        hname.ReplaceAll("pdf_", "");
                        hname = TString::Format("pdfVariation%s", hname.Data());
                    }
                    if (variation_dirs[i].Contains("Scale"))
                    {
                        hname.ReplaceAll("scale_", "");
                        hname = TString::Format("scaleWeight%s", hname.Data());
                    } 
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
            //per bin
            for (int bin = 1; bin <= originalHistograms[0]->GetNbinsX(); bin++)
            {   
                //error divided by 2 due to amc@nlo bug 
                uncHistograms16pre[i]->SetBinError(bin , uncHistograms16pre[i]->GetBinError(bin)/2); 
                uncHistograms16post[i]->SetBinError(bin, uncHistograms16post[i]->GetBinError(bin)/2);
                uncHistograms17[i]->SetBinError(bin,  uncHistograms17[i]->GetBinError(bin)/2);
                uncHistograms18[i]->SetBinError(bin, uncHistograms18[i]->GetBinError(bin)/2);
                
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

        // now you need to scale with 1/LUMI and width 
        cout<< LUMI << endl;
        resultsHisto->Scale(1/LUMI, "width");

        //TH1F *resultsHisto_norm = (TH1F*)resultsHisto->Clone(TString::Format("%s_norm", 
        //                    resultsHisto->GetName()));
        //resultsHisto_norm->Scale(LUMI/resultsHisto->Integral());  

        // scale the variations
        //TH1F *uncHistograms_total_norm[AnalysisConstants::variations.size()];

        for (int unc = 0; unc<AnalysisConstants::variations.size(); unc++)
        {
            uncHistograms_total[unc]->Scale(1/LUMI, "width");
            //cout<<"NAME "<<uncHistograms_total[unc]->GetName()<<endl;
            //uncHistograms_total_norm[unc] = (TH1F*)uncHistograms_total[unc]->Clone(TString::Format("%s_norm", 
            //                uncHistograms_total[unc]->GetName()));
            //uncHistograms_total_norm[unc]->Scale(LUMI/uncHistograms_total[unc]->Integral());   
        }

        outFile->cd();
        resultsHisto->Write();
        //resultsHisto_norm->Write(TString::Format("%s_norm", resultsHisto->GetName()));
        for (int unc = 0; unc<AnalysisConstants::variations.size(); unc++)
        {
            uncHistograms_total[unc]->Write();
            //uncHistograms_total_norm[unc]->Write(TString::Format("%s_norm", uncHistograms_total[unc]->GetName()));
        }

        delete bins;
        delete resultsHisto;
        
        uncHistograms16pre.clear();
        uncHistograms16post.clear();
        uncHistograms17.clear();
        uncHistograms18.clear();
        originalHistograms.clear();

    } // loop on all variables    

    outFile->Close();
} //oef 

