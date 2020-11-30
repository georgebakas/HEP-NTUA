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
std::vector<TString> listFiles(const char *dirname="", const char *ext=".root")
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
       if (!file->IsDirectory() && fname.EndsWith(ext)) {
         //cout << fname.Data() << endl;
         list_of_files.push_back(fname.Data());
       }
     }
   }
   return list_of_files;
}


void Systematics(TString year)
{
  gStyle->SetOptStat(0);
  TString isParton = "Parton";
  initFilesMapping();
  TString baseDir = TString::Format("%s//results", year.Data());
  std::vector<TString> dirs = {"Nominal", "JES", "bTagVariation", "SystematicsFiles", "PSWeights", "PDFWeights"};
  std::vector<TString> groups = {"Stat. Uncertainty", "JES+JER+Pileup", "Flavor Tagging", "Parton Shower", "Hard Scattering"};
  std::vector<int> groupColors = {kBlack, kRed, kBlue, kGreen, kOrange};

  TString outputDirectory = TString::Format("%s/results/", year.Data());
  //CheckAndCreateDirectory(outputDirectory);
  const int NVAR = 7;
  TString fileName = TString::Format("../Unfolding/%s/%sMeasurements/Data/OutputFile.root",
                                     year.Data(),
                                     isParton.Data());
  TFile *fNominal = TFile::Open(fileName);
  for (int i = 0; i <1; i++)
  {
    std::cout << "Variable: " << vars[i]<< std::endl;

    std::vector<TH1F *> groupHistogramsUp;
    std::vector<TH1F *> groupHistogramsDown;
    std::vector<TH1F *> groupHistogramsSym;

    TH1F *hUnfoldedNominal = (TH1F *)fNominal->Get(TString::Format("hUnfold_%s", vars[i].Data()));

    //initialize group histograms
    for (int group = 0; group < groups.size(); group++)
    {
      TH1F *hSystematicsUp = (TH1F *)hUnfoldedNominal->Clone(TString::Format("%s %s Up", groups[group].Data(),
                                                                             vars[i].Data()));
      TH1F *hSystematicsDown = (TH1F *)hUnfoldedNominal->Clone(TString::Format("%s %s Down", groups[group].Data(),
                                                                               vars[i].Data()));
      TH1F *hSystematicsSym = (TH1F *)hUnfoldedNominal->Clone(TString::Format("%s %s Sym", groups[group].Data(),
                                                                              vars[i].Data()));
      hSystematicsUp->Reset();
      hSystematicsDown->Reset();
      hSystematicsSym->Reset();
      hSystematicsUp->SetLineColor(group);
      hSystematicsDown->SetLineColor(group);
      hSystematicsSym->SetLineColor(group);

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
      double nominalValue = hUnfoldedNominal->GetBinContent(bin + 1);
      double nominalError = hUnfoldedNominal->GetBinError(bin + 1);
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
      std::cout << "Variation: " << dirs[j] << std::endl;
      TString variation = dirs[j];
      std::vector<TString> variationFiles = listFiles(TString::Format("%s/Unfolding_%s/",year.Data(), dirs[j].Data()));
      int group = -1;
      if (variation.Contains("Nominal"))
        group = 0;
      if (variation.Contains("JESSrc") || variation.Contains("JER") || variation.Contains("Pileup"))
        group = 1;
      if (variation.Contains("bTagVariation"))
        group = 2;
      if (variation.Contains("SystematicsFiles"))
        group = 3;
      if (variation.Contains("PDF") || variation.Contains("PS") || variation.Contains("Scale"))
        group = 4;

      cout<<variation<<" group: "<<group <<endl;

      if (group == -1)
        continue;

      for(int jfile=0; jfile<variationFiles.size(); jfile++)
      {
        fileName = TString::Format("%s/Unfolding_%s/%s",
                                   year.Data(),
                                   variation.Data(),
                                   variationFiles[jfile].Data());

        TFile *f = TFile::Open(fileName);
        cout<<fileName<<endl;
        if(fileName.EqualTo("2018/Unfolding_PDFWeights/OutputFileParton_pdf_35.root")) continue;
        TH1F *hVariation = (TH1F *)f->Get(TString::Format("hUnfold_%s", vars[i].Data()));

        for (int bin = 0; bin < groupHistogramsUp[group]->GetNbinsX(); bin++)
        {
          double nominalValue = hUnfoldedNominal->GetBinContent(bin + 1);
          double nominalError = hUnfoldedNominal->GetBinError(bin + 1);
          double variationValue = hVariation->GetBinContent(bin + 1);
          double valuePull = (variationValue - nominalValue) / nominalValue;
          double variationErrorUp = groupHistogramsUp[group]->GetBinContent(bin + 1);
          double variationErrorDown = groupHistogramsDown[group]->GetBinContent(bin + 1);
          //std::cout << nominalValue << " " << nominalError;
          //std::cout << " " << variationValue << " " << valuePull;
          //std::cout << " " << variationErrorUp << " " << variationErrorDown;
          if (nominalValue != 0)
          {
            if (valuePull >= 0)
            {
              //std::cout << " Variation up:";
              //std::cout << " " << groupHistogramsUp[group]->GetBinContent(bin + 1);
              //std::cout << " " << groupHistogramsDown[group]->GetBinContent(bin + 1);
              groupHistogramsUp[group]->SetBinContent(bin + 1, variationErrorUp + (valuePull * valuePull));
              //std::cout << " " << groupHistogramsUp[group]->GetBinContent(bin + 1);
              //std::cout << " " << groupHistogramsDown[group]->GetBinContent(bin + 1);
            }
            else
            {
              //std::cout << " Variation down:";
              //std::cout << " " << groupHistogramsUp[group]->GetBinContent(bin + 1);
              //std::cout << " " << groupHistogramsDown[group]->GetBinContent(bin + 1);
              groupHistogramsDown[group]->SetBinContent(bin + 1, variationErrorDown + (valuePull * valuePull));
              //std::cout << " " << groupHistogramsUp[group]->GetBinContent(bin + 1);
              //std::cout << " " << groupHistogramsDown[group]->GetBinContent(bin + 1);
            }
          }
          //std::cout << endl;
        }
        //std::cout << endl;
        //f->Close();
      }//end of for loop (over specific variation files)
    } //end of for loop (over all variation types)


    for (int group = 0; group < groups.size(); group++)
    {
      //std::cout << "Group: " << groups[group] << std::endl;
      for (int bin = 0; bin < groupHistogramsUp[group]->GetNbinsX(); bin++)
      {
        double errorUp = TMath::Sqrt(groupHistogramsUp[group]->GetBinContent(bin + 1));
        double errorDown = TMath::Sqrt(groupHistogramsDown[group]->GetBinContent(bin + 1));
        //std::cout << " " << errorUp << " " << errorDown;
        groupHistogramsSym[group]->SetBinContent(bin + 1, 0.5 * 100 * (errorUp + errorDown));
        //std::cout << " " << groupHistogramsSym[group]->GetBinContent(bin + 1) << std::endl;
      }
    }
    TCanvas *c1 = new TCanvas(vars[i], vars[i], 800, 600);

    groupHistogramsSym[0]->SetFillColor(kGray);
    //groupHistogramsSym[0]->GetYaxis()->SetRangeUser(0, 120);
    groupHistogramsSym[0]->SetLineWidth(0);
    groupHistogramsSym[0]->GetXaxis()->SetTitle(vars[i]);
    groupHistogramsSym[0]->GetXaxis()->SetLabelSize(0.035);
    std::cout << groupHistogramsSym[0]->GetYaxis()->GetLabelSize() << std::endl;
    groupHistogramsSym[0]->Draw("hist");
    //groupHistogramsSym[0]->GetYaxis()->SetRangeUser(0, 120);
    groupHistogramsSym[0]->GetYaxis()->SetTitle("Relative Uncertainty [%]");
    //groupHistogramsSym[0]->Draw();
    if (groupHistogramsSym[0]->GetMaximum() < 120)
    {
      groupHistogramsSym[0]->GetYaxis()->SetRangeUser(0, 120);
    }
    TLegend *leg = new TLegend(0.2, 0.6, 0.5, 0.9);
    leg->AddEntry(groupHistogramsSym[0], groups[0], "f");

    for (int group = 1; group < groups.size(); group++)
    {
      groupHistogramsSym[group]->Draw("hist same");
      leg->AddEntry(groupHistogramsSym[group], groups[group], "l");
    }

    leg->Draw();

    lumi_13TeV = TString::Format("%0.1f fb^{-1}", luminosity["luminosity"+year]/1000);
    writeExtraText = true;
    CMS_lumi(c1, 4, 10);
    c1->SaveAs(TString::Format("%s/Systematics_%s.png", outputDirectory.Data(), vars[i].Data()));
  }
}
