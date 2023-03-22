/*
    In this code I will calculate the tag and probe efficiency of the top tagger
    N(1) probe is: #events pass SRb with top tight cut and SRb without any top tagger CUT
    N(2) probe is: #events pass SRb with top tight cut and SR (standard top tagger cut)

    N(1) serves as the denominator
    N(2) serves as the numerator

    as efficiency:= N(2) / N(1) for ttbar files, and data - QCD - subdominant files
*/

#include <stdio.h>
#include "TemplateConstants.h"

void CalculateSF_FitProcess(TString year = "2017")
{
    initFilesMapping();
    //get the files from the directory
    //data file
    TFile *infData = TFile::Open(TString::Format("../MassFit/%s/TagAndProbeHisto_Data_newBins.root", year.Data()));
    //tt nominal file:
    TFile *infTT = TFile::Open(TString::Format("%s/Nominal/combined/TagAndProbeHisto_1000_TT_Nominal_newBins.root", year.Data()));
    //qcd mc file
    TFile *infQCD = TFile::Open(TString::Format("../MassFit/%s/TagAndProbeHisto_QCD_HT300toInf_newBins.root",year.Data()));
    //subdominant file:
    TFile *infSub = TFile::Open(TString::Format("../MassFit/%s/TagAndProbeHisto_SubdominantBkgs_newBins.root",year.Data()));

    TString regions[2] = {"hSRBTightAndProbe_", "hSRBTightAndSR_"};
    //get the histograms for the jetPt0 variable
    //[0] is denominator and [1] is numerator
    TH1F *hData[2], *hQCD[2], *hSub[2], *hTT[2];
    hData[0] = (TH1F*)infData->Get("hSRBTightAndProbe_jetPt0_expYield");
    hData[1] = (TH1F*)infData->Get("hSRBTightAndSR_jetPt0_expYield");

    hTT[0] = (TH1F*)infTT->Get("hSRBTightAndProbe_jetPt0_expYield");
    hTT[1] = (TH1F*)infTT->Get("hSRBTightAndSR_jetPt0_expYield");

    hSub[0] = (TH1F*)infSub->Get("hSRBTightAndProbe_jetPt0_expYield");
    hSub[1] = (TH1F*)infSub->Get("hSRBTightAndSR_jetPt0_expYield");

    hQCD[0] = (TH1F*)infQCD->Get("hSRBTightAndProbe_jetPt0_expYield");
    hQCD[1] = (TH1F*)infQCD->Get("hSRBTightAndSR_jetPt0_expYield");

    cout<<"hData Numerator: "<<hData[1]->Integral()<<endl;
    cout<<"hData denominator: "<<hData[0]->Integral()<<endl;

    cout<<"hSub Numerator: "<<hSub[1]->Integral()<<endl;
    cout<<"hSub denominator: "<<hSub[0]->Integral()<<endl;

    cout<<"hQCD Numerator: "<<hQCD[1]->Integral()<<endl;
    cout<<"hQCD denominator: "<<hQCD[0]->Integral()<<endl;

    cout<<"hTT Numerator: "<<hTT[1]->Integral()<<endl;
    cout<<"hTT denominator: "<<hTT[0]->Integral()<<endl;

    for(int i =0; i<2; i++)
    { 
    // scale QCD to its shape
    // get the NQCD 
    TFile *masFitResultsFile = TFile::Open(TString::Format("%s/Nominal/MassFitResults_%sSignalTemplates_.root",
                                                            year.Data(),  
                                                            regions[i].Data()));
    RooWorkspace *w = (RooWorkspace*)masFitResultsFile->Get("w");
    RooRealVar *value = (RooRealVar *)w->var("nFitQCD_2b");

    float val = value->getValV();
    float nqcd_error = value->getError();
    masFitResultsFile->Close();


    Double_t integral_error;
    float integral_qcd = hQCD[i]->IntegralAndError(1, hQCD[i]->GetNbinsX(), integral_error);
    for(int ibin=1; ibin<=hQCD[i]->GetNbinsX(); ibin++)
    {
        float error = hQCD[i]->GetBinError(ibin);
        float bin_value = hQCD[i]->GetBinContent(ibin);
        float new_error = TMath::Sqrt(TMath::Power((error*val)/integral_qcd,2) + 
                                TMath::Power((nqcd_error*bin_value)/integral_qcd, 2) + 
                                TMath::Power(bin_value*val*integral_error/TMath::Power(integral_qcd,2),2));
        hQCD[i]->SetBinError(ibin, new_error);
        hQCD[i]->SetBinContent(ibin, val*hQCD[i]->GetBinContent(ibin)/integral_qcd);
    }
    //hQCD[i]->Scale(val/integral_qcd);

    } 
    for(int i =0; i<2; i++)
    { 
        for(int ibin=1; ibin<=hData[i]->GetNbinsX(); ibin++)
        { 
            float new_error = TMath::Sqrt(TMath::Power(hQCD[i]->GetBinError(ibin) ,2)+ 
                                    TMath::Power(hSub[i]->GetBinError(ibin) ,2) +
                                    TMath::Power(hData[i]->GetBinError(ibin), 2));
            hData[i]->SetBinError(ibin, new_error);
            hData[i]->SetBinContent(ibin, hData[i]->GetBinContent(ibin) - hSub[i]->GetBinContent(ibin) - hQCD[i]->GetBinContent(ibin));
        }
    }

    //now measure the efficiency for data and ttbar files
    float eff_data = hData[1]->Integral() / hData[0]->Integral();
    float eff_tt   = hTT[1]->Integral() / hTT[0]->Integral();

    Double_t error_data[2], error_tt[2];
    Double_t intData[2], intTT[2];

    for(int i=0; i<2; i++)
    {
        intData[i] = hData[i]->IntegralAndError(1,hData[i]->GetNbinsX(),error_data[i]);
        intTT[i] = hTT[i]->IntegralAndError(1,hTT[i]->GetNbinsX(),error_tt[i]);
    }

    //calculate errors:
    float eff_data_error = TMath::Sqrt( TMath::Power(error_data[1]/intData[0],2) + TMath::Power(error_data[0] * intData[1]/ TMath::Power(intData[0],2),2));
    float eff_tt_error = TMath::Sqrt( TMath::Power(error_tt[1]/intTT[0],2) + TMath::Power(error_tt[0] * intTT[1]/ TMath::Power(intTT[0],2),2));

    // calculate top taggger SF inclusive
    // calculate the SF error here 
    float top_tagger_sf_error_incl = TMath::Sqrt( 
                                TMath::Power(eff_data_error/eff_tt,2) +
                                TMath::Power(eff_data*eff_tt_error/TMath::Power(eff_tt,2),2));
    float top_tagger_sf_incl = (eff_data/eff_tt);
    cout<<"THIS IS TOP TAGGER SF INCLUSIVE: "<<top_tagger_sf_incl<<endl;
    //now calculate it per pT region:
    //jetPt0 binning: {400,425,450,475,500,535,570,610,650,700,750,800,850,900,950,1025,1100,1200,1300,1400,1500}
    //[i][0] is the ptregion i, denominator and [i][1] is pt region i, numerator
    Double_t integralPtRegionsData[3][2],integralPtRegionsData_error[3][2];
    Double_t integralPtRegionsTT[3][2], integralPtRegionsTT_error[3][2];

    float eff_PtRegionsData[3];
    float eff_PtRegionsTT[3];
    //calculate errors in pt regions:
    float eff_data_errorPtRegions[3];
    float eff_tt_errorPtRegions[3];
    float eff_tt_errorPtRegions_systematic[3];
    // scale factor calculation per pt region:
    float top_tagger_sf_error[3];
    float top_tagger_sf[3];

    float start[3] = {1,5,6};
    float end[3] = {4,5,7};

    for(int ibin=1; ibin<=hData[0]->GetNbinsX(); ibin++)
    {
        cout<< "ibin "<< ibin <<": Content: "<< hData[0]->GetBinContent(ibin)<<endl;
    }
        
    for(int ipt=0; ipt<3;ipt++)
    {
        cout<< "integal"<< hData[0]->IntegralAndError(start[ipt],end[ipt],integralPtRegionsData_error[ipt][0])<<endl; 
        integralPtRegionsData[ipt][0] = hData[0]->IntegralAndError(start[ipt],end[ipt],integralPtRegionsData_error[ipt][0]);
        integralPtRegionsData[ipt][1] = hData[1]->IntegralAndError(start[ipt],end[ipt],integralPtRegionsData_error[ipt][1]);

        integralPtRegionsTT[ipt][0] = hTT[0]->IntegralAndError(start[ipt],end[ipt],integralPtRegionsTT_error[ipt][0]);
        integralPtRegionsTT[ipt][1] = hTT[1]->IntegralAndError(start[ipt],end[ipt],integralPtRegionsTT_error[ipt][1]);

        eff_PtRegionsData[ipt] = integralPtRegionsData[ipt][1]/integralPtRegionsData[ipt][0];
        eff_PtRegionsTT[ipt] = integralPtRegionsTT[ipt][1]/integralPtRegionsTT[ipt][0];
        //eff_PtRegionsTT[ipt] = ttbar_mc_tagnprobe_eff[year.Data()][start[ipt]];
        
        eff_data_errorPtRegions[ipt] = TMath::Sqrt( 
                        TMath::Power(integralPtRegionsData_error[ipt][1]/integralPtRegionsData[ipt][0],2) +  // 
                        TMath::Power(integralPtRegionsData_error[ipt][0] * integralPtRegionsData[ipt][1]/ TMath::Power(integralPtRegionsData[ipt][0],2),2));
        cout<<ipt<<"  "<<eff_PtRegionsData[ipt]<<endl;
        cout<<ipt<<"  "<<eff_PtRegionsTT[ipt]<<endl;
        // this is statistical error
        eff_tt_errorPtRegions[ipt] = TMath::Sqrt( TMath::Power(integralPtRegionsTT_error[ipt][1]/integralPtRegionsTT[ipt][0],2) + TMath::Power(integralPtRegionsTT_error[ipt][0] * integralPtRegionsTT[ipt][1]/ TMath::Power(integralPtRegionsTT[ipt][0],2),2));
        // this is systematic error
        eff_tt_errorPtRegions_systematic[ipt] = ttbar_mc_tagnprobe_eff_error[year.Data()][start[ipt]];

        // calculate the SF error here 
        top_tagger_sf_error[ipt] = TMath::Sqrt( 
                                    TMath::Power(eff_data_errorPtRegions[ipt]/eff_PtRegionsTT[ipt],2) +
                                    TMath::Power(eff_PtRegionsData[ipt]*eff_tt_errorPtRegions[ipt]/TMath::Power(eff_PtRegionsTT[ipt],2),2));
        top_tagger_sf[ipt] = (eff_PtRegionsData[ipt]/eff_PtRegionsTT[ipt]);  
        
    }



    // create a histogram 
    Float_t xbins[4] = {400, 600, 800, 1500};
    int NBINS = 3;
    TH1F *hTopSF = new TH1F("h_topTaggerSF","h_topTaggerSF", NBINS, xbins);
    
    for (int ibin=0; ibin<=hTopSF->GetNbinsX(); ibin++)
    {
        // use as a central value the inclusive SF
        hTopSF->SetBinContent(ibin+1, top_tagger_sf_incl);
        // error = central value * error_pct 
        // where error_pct = top_tagger_sf_error[ipt] / top_tagger_sf[ipt]
        // error = top_tagger_sf_incl * (top_tagger_sf_error[ipt]/top_tagger_sf[ipt])
        hTopSF->SetBinError(ibin+1, top_tagger_sf_incl * (top_tagger_sf_error[ibin]/top_tagger_sf[ibin]));
    }

    std::vector<float> y_up_values, y_down_values;
    std::vector<float> x_values;
    // Now fit the up/down values 
    for (int ibin=0; ibin<=hTopSF->GetNbinsX(); ibin = ibin+2)
    {
        float x = hTopSF->GetBinCenter(ibin+1);
        float yup = hTopSF->GetBinContent(ibin+1) + hTopSF->GetBinError(ibin+1);
        float ydown = hTopSF->GetBinContent(ibin+1) - hTopSF->GetBinError(ibin+1);
        cout<<"up: "<< x << " "<< yup<< endl;
        cout<<"down: "<< x << " "<< ydown<< endl;
        y_up_values.push_back(yup);
        y_down_values.push_back(ydown); 
        x_values.push_back(x);

    }
    // y = ax + b --> 1 for up and 1 for down
    float alpha_up = (y_up_values[1] - y_up_values[0]) / (x_values[1] - x_values[0]);
    float beta_up = y_up_values[0] - alpha_up * x_values[0];
    
    float alpha_down = (y_down_values[1] - y_down_values[0]) / (x_values[1] - x_values[0]);
    float beta_down = y_down_values[0] - alpha_down * x_values[0];
    
    // now do the tf1s
    TF1 *fup = new TF1("fup", "[0] *x + [1]", 400, 1500);
    TF1 *fdown = new TF1("fdown", "[0] *x + [1]", 400, 1500);

    fup->SetParameters(alpha_up, beta_up);
    fdown->SetParameters(alpha_down, beta_down);

    TCanvas *can = new TCanvas("effCanPt", "effCanPt", 800, 600);

    //testing 
    cout<<fup ->Eval(600)<<endl;
    cout<<fup ->Eval(1400)<<endl;
    cout<<fup ->Eval(5000)<<endl;

    cout<<fdown ->Eval(600)<<endl;
    cout<<fdown ->Eval(1400)<<endl;
    cout<<fdown ->Eval(5000)<<endl;


    TLegend *leg;
    //if (year.Contains("2016")) leg = new TLegend(0.5, 0.75, 0.7, 0.9);
    //else 
    leg = new TLegend(0.5, 0.15, 0.7, 0.3);

    int n = 3;
    float xData[] = {500,700,1150};
    std::vector<TString> x_labels = {"pT[400-500]", "pT[500-600]", "pT[600-Inf]"};
    float ex[]= {0,0,0};
    TGraphErrors *grData = new TGraphErrors(n,xData,top_tagger_sf,ex,top_tagger_sf_error);

    leg->AddEntry(grData,"pT dep. SF", "lep");
    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle(TString::Format("Top Tagger Scale Factors for %s; p_{T} (GeV); Scale Factor", year.Data()));

    grData->SetTitle("TGraphErrors Example");
    grData->SetMarkerColor(kBlue);
    grData->SetMarkerStyle(21);

    // SF central line 
    Double_t plotX_data[3] = {400, 700, 1500};
    Double_t plotY_data[3] = {top_tagger_sf_incl, top_tagger_sf_incl, top_tagger_sf_incl};
    Double_t plotX_data_error[3] = {0, 0, 0};
    Double_t plotY_data_error[3] = {top_tagger_sf_error_incl, top_tagger_sf_error_incl, top_tagger_sf_error_incl};
    TGraphErrors *data_central_err = new TGraphErrors(3,plotX_data,plotY_data, plotX_data_error, plotY_data_error);
    data_central_err->SetFillColor(kBlue);
    data_central_err->SetLineColor(kBlue);
    data_central_err->SetLineWidth(2.0);
    data_central_err->SetFillStyle(3354);

    leg->AddEntry(data_central_err,"Inclusice Top Tagger SF", "f");

    // now plot the tf1 
    TGraphAsymmErrors *fup_graph = new TGraphAsymmErrors(2);
    fup_graph->SetPoint(0, 400, top_tagger_sf_incl);
    fup_graph->SetPointEYhigh(0, fup->Eval(400) - top_tagger_sf_incl);
    fup_graph->SetPointEYlow(0, top_tagger_sf_incl - fdown->Eval(400));

    fup_graph->SetPoint(1, 1500, top_tagger_sf_incl);
    fup_graph->SetPointEYhigh(1, fup->Eval(1500) - top_tagger_sf_incl);
    fup_graph->SetPointEYlow(1, top_tagger_sf_incl - fdown->Eval(1500));
    fup_graph->SetLineColor(kRed);
    fup_graph->SetMarkerColor(kRed);
    fup_graph->SetFillColor(kRed);
    fup_graph->SetLineWidth(2.0);
    fup_graph->SetFillStyle(3345);

    leg->AddEntry(fup_graph,"pT dependent uncertainty", "f");

    mg->Add(grData, "AP");
    mg->Add(data_central_err, "A3L");
    mg->Add(fup_graph, "A3L");

    if (year.EqualTo("2016_preVFP"))
        mg->GetYaxis()->SetRangeUser(0.4, 1.4);
    else if (year.EqualTo("2016_postVFP"))
        mg->GetYaxis()->SetRangeUser(0.4, 1.5);  
    else 
        mg->GetYaxis()->SetRangeUser(0.5, 1.2);  
    mg->GetXaxis()->SetRangeUser(350, 1500);  
    mg->Draw("a");


    leg->Draw();
    can->Update();
    can->Print(TString::Format("%s/Nominal/plots/SF_perPtRegion_FitProcess.pdf",year.Data()),"pdf");


    // save the TF1s
    TFile *fit_file = new TFile(TString::Format("../VariationHandling/%s/TopTaggerUncertainty_FitValues.root",year.Data()),"RECREATE");
    fit_file->cd();
    fup->Write("fup");
    fdown->Write("fdown");

}
