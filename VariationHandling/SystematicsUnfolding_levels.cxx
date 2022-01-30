#include "TemplateConstants.h"
#include "../CMS_plots/tdrstyle.C"
#include "../CMS_plots/CMS_lumi.C"

//We split the uncertaintied in 5 groups and we create a histogram for each group
//1. Statistical which is basically the error bars of the unfolded histograms
//2. JES+JER+Pileup which is the variation (difference from the nominal restuls)
//   we get from these source
//3. Flavor tagging
//4. Parton Shower
//5. Hard scattering

//For every group we create a histogram that will containt the combined uncertainty
//that resulsts from the addition of the uncertainties from this gorup.
//The combination of the uncertainties comes by adding the in qudrature. Although
//we split the uncertainties in up and down depending on wheter they result in a
//upper or lower value for the XS in that bin we in the end present one combined
//which is given by 0.5 * 100 * (eUp + eDown)
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

void SystematicsUnfolding_levels(TString isParton = "Particle", bool isNorm = false)
{
  gStyle->SetOptStat(0);
  initFilesMapping();
  setTDRStyle();
  std::vector<TString> dirs = {"Nominal", "JES", "bTagVariation", "PSWeights", "PDFWeights", "ScaleWeights"};
  std::vector<TString> groups = {"Stat. Uncertainty", "JES+JER+Pileup", "Flavor Tagging", "Parton Shower", "Hard Scattering"};
  std::vector<int> groupColors = {kBlack, kRed, kGreen, kBlue, kOrange};

  //TString baseInputDir = "/afs/cern.ch/work/g/gbakas/public/HEP-NTUA/";
  TString baseInputDir = "/Users/georgebakas/Documents/HEP-NTUA_ul/VariationHandling/";
  TString outputDirectory = TString::Format("/SystematicBreakdown/results");
  //CheckAndCreateDirectory(outputDirectory);
  TString histo_name = "hUnfoldFinal_";
  //TString histo_name = "hSigInit_";
  
  const int NVAR = 10;

  TString fileName = TString::Format("%sUnfoldedCombined/Nominal/OutputFile%s.root", baseInputDir.Data(), isParton.Data());
  TFile *fNominal = TFile::Open(fileName);

  for (int i = 0; i<NVAR; i++)
  {
    std::cout << "Variable: " << vars[i]<< std::endl;

    std::vector<TH1F *> groupHistogramsUp;
    std::vector<TH1F *> groupHistogramsDown;
    std::vector<TH1F *> groupHistogramsSym;

    TH1F *hNominal = (TH1F *)fNominal->Get(TString::Format("%s%s", histo_name.Data(), vars[i].Data()));
    if(isNorm)
      hNominal->Scale(1./hNominal->Integral());

    //initialize group histograms
    for (int group = 0; group < groups.size(); group++)
    {
      TH1F *hSystematicsUp = (TH1F *)hNominal->Clone(TString::Format("%s %s Up", groups[group].Data(),
                                                                            vars[i].Data()));
      TH1F *hSystematicsDown = (TH1F *)hNominal->Clone(TString::Format("%s %s Down", groups[group].Data(),
                                                                              vars[i].Data()));
      TH1F *hSystematicsSym = (TH1F *)hNominal->Clone(TString::Format("%s %s Sym", groups[group].Data(),
                                                                              vars[i].Data()));
      hSystematicsUp->Reset();
      hSystematicsDown->Reset();
      hSystematicsSym->Reset();
      hSystematicsUp->SetLineColor(groupColors[group-1]);
      hSystematicsDown->SetLineColor(groupColors[group-1]);
      hSystematicsSym->SetLineColor(groupColors[group-1]);

      hSystematicsUp->SetLineWidth(3);
      hSystematicsDown->SetLineWidth(3);
      hSystematicsSym->SetLineWidth(3);

      groupHistogramsUp.push_back(hSystematicsUp);
      groupHistogramsDown.push_back(hSystematicsDown);
      groupHistogramsSym.push_back(hSystematicsSym);
    }

    //Start of nominal handling
    for (int bin = 0; bin < groupHistogramsUp[0]->GetNbinsX(); bin++)
    {
      double nominalValue = hNominal->GetBinContent(bin + 1);
      double nominalError = hNominal->GetBinError(bin + 1);
      if (nominalValue != 0)
      {
        groupHistogramsUp[0]->SetBinContent(bin + 1, (nominalError * nominalError) / (nominalValue * nominalValue));
        groupHistogramsDown[0]->SetBinContent(bin + 1, (nominalError * nominalError) / (nominalValue * nominalValue));
      }
      //std::cout << "nominalValue: " << nominalValue << " nominalError: " << nominalError << " " << (nominalError * nominalError) / (nominalValue * nominalValue) << std::endl;
    }
    //end of nominal handling

    //start variation handling

    for (int j = 0; j < dirs.size(); j++)
    {
      //std::cout << "Variation: " << dirs[j] << std::endl;
      TString variation = dirs[j];
      std::vector<TString> variationFiles = listFiles(TString::Format("UnfoldedCombined/%s/", dirs[j].Data()), isParton.Data());

      int group = -1;
      if (variation.Contains("Nominal"))
        group = 0;
      if (variation.Contains("JES") || variation.Contains("JER") || variation.Contains("Pileup"))
        group = 1;
      if (variation.Contains("bTagVariation"))
        group = 2;
      if (variation.Contains("PS"))
        group = 3;
      if (variation.Contains("PDF") || variation.Contains("Scale"))
        group = 4;
      if (variation.Contains("SystematicsFiles"))
        group = -1;

      cout<<variation<<" group: "<<group <<endl;

      if (group == -1)
        continue;

      for(int jfile=0; jfile<variationFiles.size(); jfile++)
      {
        /* fileName = TString::Format("UnfoldedCombined/%s/%s",
                                  variation.Data(),
                                   variationFiles[jfile].Data()); */

        fileName = TString::Format("UnfoldedCombined/%s/%s",
                                  variation.Data(),
                                  variationFiles[jfile].Data());
        
        if (!fileName.Contains("Def") && variation.Contains("PS"))
          continue;
        
        TFile *f = TFile::Open(fileName);
        
        cout<<fileName<<endl;
        /*if (fileName.EqualTo(TString::Format("UnfoldedCombined/PDFWeights/OutputFile%s_pdf_0.root", isParton.Data()))) continue;
        if (fileName.EqualTo(TString::Format("UnfoldedCombined/PDFWeights/OutputFile%s_pdf_1.root", isParton.Data()))) continue;
        if (fileName.EqualTo(TString::Format("UnfoldedCombined/PDFWeights/OutputFile%s_pdf_1.root", isParton.Data()))) continue;
        if (fileName.EqualTo(TString::Format("UnfoldedCombined/PDFWeights/OutputFile%s_pdf_56.root", isParton.Data()))) continue;
        if (fileName.EqualTo(TString::Format("UnfoldedCombined/PDFWeights/OutputFile%s_pdf_57.root", isParton.Data()))) continue;
        if (fileName.EqualTo(TString::Format("UnfoldedCombined/PDFWeights/OutputFile%s_pdf_35.root", isParton.Data()))) continue;
        if (fileName.EqualTo(TString::Format("UnfoldedCombined/PDFWeights/OutputFile%s_pdf_42.root", isParton.Data()))) continue;
        if (fileName.EqualTo(TString::Format("UnfoldedCombined/PDFWeights/OutputFile%s_pdf_76.root", isParton.Data()))) continue;
        if (fileName.EqualTo(TString::Format("UnfoldedCombined/PDFWeights/OutputFile%s_pdf_59.root", isParton.Data()))) continue;
        if (fileName.EqualTo(TString::Format("UnfoldedCombined/PDFWeights/OutputFile%s_pdf_58.root", isParton.Data()))) continue;
        if (fileName.EqualTo(TString::Format("UnfoldedCombined/PDFWeights/OutputFile%s_pdf_86.root", isParton.Data()))) continue;
        if (fileName.EqualTo(TString::Format("UnfoldedCombined/PDFWeights/OutputFile%s_pdf_93.root", isParton.Data()))) continue;*/
        if (fileName.EqualTo(TString::Format("UnfoldedCombined/PDFWeights/OutputFile%s_pdf_98.root", isParton.Data()))) continue;
        if (fileName.EqualTo(TString::Format("UnfoldedCombined/PDFWeights/OutputFile%s_pdf_99.root", isParton.Data()))) continue;
        if (fileName.EqualTo(TString::Format("UnfoldedCombined/PDFWeights/OutputFile%s_pdf_100.root", isParton.Data()))) continue;
        //if (fileName.EqualTo("UnfoldedCombined/ScaleWeights/OutputFileParton_scale_9.root")) continue;
        //if (fileName.EqualTo("UnfoldedCombined/ScaleWeights/OutputFileParton_scale_7.root")) continue;


        TH1F *hVariation = (TH1F *)f->Get(TString::Format("%s%s", histo_name.Data(), vars[i].Data()));
        if(isNorm)
          hVariation->Scale(1./hVariation->Integral());

        for (int bin = 0; bin < groupHistogramsUp[group]->GetNbinsX(); bin++)
        {
          double nominalValue = hNominal->GetBinContent(bin + 1);
          double nominalError = hNominal->GetBinError(bin + 1);
          double variationValue = hVariation->GetBinContent(bin + 1);
          double valuePull = (variationValue - nominalValue) / nominalValue;
          double variationErrorUp = groupHistogramsUp[group]->GetBinContent(bin + 1);
          double variationErrorDown = groupHistogramsDown[group]->GetBinContent(bin + 1);

          //std::cout << nominalValue << " " << nominalError << endl;
          //std::cout << " " << variationValue << " " << valuePull << endl;
          std::cout << fabs(variationValue - nominalValue) /nominalValue << endl;
          if (fabs(variationValue - nominalValue) /nominalValue > 0.7 ) 
            std::cout << fileName <<endl;
          
          //std::cout << " " << variationErrorUp << " " << variationErrorDown;
          if (nominalValue != 0)
          {
            if (valuePull >= 0)
            {
              //std::cout << " Variation up:";
              //std::cout << " " << groupHistogramsUp[group]->GetBinContent(bin + 1);
              //std::cout << " " << groupHistogramsDown[group]->GetBinContent(bin + 1);
              groupHistogramsUp[group]->SetBinContent(bin + 1, variationErrorUp + (valuePull * valuePull));
              //std::cout << " " << groupHistogramsUp[group]->GetBinContent(bin + 1) <<endl;
              //std::cout << " " << groupHistogramsDown[group]->GetBinContent(bin + 1);
            }
            else
            {
              //std::cout << " Variation down:";
              //std::cout << " " << groupHistogramsUp[group]->GetBinContent(bin + 1);
              //std::cout << " " << groupHistogramsDown[group]->GetBinContent(bin + 1);
              groupHistogramsDown[group]->SetBinContent(bin + 1, variationErrorDown + (valuePull * valuePull));
              //std::cout << " " << groupHistogramsUp[group]->GetBinContent(bin + 1) <<endl;
              //std::cout << " " << groupHistogramsDown[group]->GetBinContent(bin + 1);
            }
          }
          //std::cout << endl;
        }
        //std::cout << endl;
        f->Close();
      }//end of for loop (over specific variation files)
    } //end of for loop (over all variation types)


    for (int group = 0; group < groups.size(); group++)
    {
      std::cout << "Group: " << groups[group] << std::endl;
      for (int bin = 0; bin < groupHistogramsUp[group]->GetNbinsX(); bin++)
      {
        double errorUp = TMath::Sqrt(groupHistogramsUp[group]->GetBinContent(bin + 1));
        double errorDown = TMath::Sqrt(groupHistogramsDown[group]->GetBinContent(bin + 1));
        std::cout <<"bin "<<bin<<" " << errorUp << " " << errorDown;
        groupHistogramsSym[group]->SetBinContent(bin + 1, 0.5 * 100 * (errorUp + errorDown));
        std::cout << " " << groupHistogramsSym[group]->GetBinContent(bin + 1) << std::endl;
      }
    }
    TCanvas *c1 = new TCanvas(vars[i], vars[i], 800, 600);

    groupHistogramsSym[0]->SetFillColor(kGray);
    //groupHistogramsSym[0]->GetYaxis()->SetRangeUser(0, 120);
    groupHistogramsSym[0]->SetTitle("");
    groupHistogramsSym[0]->SetLineWidth(0);
    groupHistogramsSym[0]->GetXaxis()->SetTitle(vars[i]);
    groupHistogramsSym[0]->GetXaxis()->SetLabelSize(0.035);
    std::cout << groupHistogramsSym[0]->GetYaxis()->GetLabelSize() << std::endl;
    groupHistogramsSym[0]->Draw("hist");
    groupHistogramsSym[0]->GetYaxis()->SetRangeUser(0, 120);
    groupHistogramsSym[0]->GetYaxis()->SetTitle("Relative Uncertainty [%]");
    //groupHistogramsSym[0]->Draw();
    if (groupHistogramsSym[0]->GetMaximum() < 120)
    {
      groupHistogramsSym[0]->GetYaxis()->SetRangeUser(0, 120);
    }
    TLegend *leg = new TLegend(0.7, 0.75, 0.9, 0.9);
    leg->AddEntry(groupHistogramsSym[0], groups[0], "f");

    for (int group = 1; group < groups.size(); group++)
    {
      groupHistogramsSym[group]->Draw("hist same");
      leg->AddEntry(groupHistogramsSym[group], groups[group], "l");
    }

    leg->Draw();

    lumi_13TeV = TString::Format("%0.1f fb^{-1}", luminosity["luminosityAll"]/1000);
    int iPeriod = 13;
    int iPos = 0;
    writeExtraText=true;
    CMS_lumi(c1, iPeriod, iPos);
    c1->SaveAs(TString::Format("%s%s/Systematics%s%s_%s%s.pdf", 
                        baseInputDir.Data(), outputDirectory.Data(), 
                        histo_name.Data(), isParton.Data(), vars[i].Data(),
                        (isNorm ? "_normalized" : "")), "pdf");
  }
}
